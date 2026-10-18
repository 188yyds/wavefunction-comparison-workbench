import copy, csv, json, os, uuid, math, re
from decimal import Decimal, InvalidOperation
from pathlib import Path

ROOT = Path(os.path.abspath(__file__)).parent
VERSION = 1
COLORS = ['#2166ac', '#d55e00', '#009e73', '#7251a2', '#c74878', '#7b6500']
PARAMS = dict(t1=1.0, ratio=.4, dimer=1, nbar=.5, L=0, N=-1,
              dataset='psi', bulk_left=0, bulk_right=0, string_left=0,
              string_right=0, cut=0, max_order=6, entropy_coarse=1e-9,
              entropy_fine=1e-10, memory_mb=1024, blas_threads=1)
BASE_OBS = [
 ('纠缠量','half_chain_entropy','半链纠缠熵','S_half'),
 ('纠缠量','center_entropy','中心块纠缠熵（较耗内存）','S_block'),
 ('纠缠量','entanglement_gap','纠缠谱差','xi2 - xi1'),
 ('纠缠量','ES1','第一条纠缠能级','xi1'),('纠缠量','ES2','第二条纠缠能级','xi2'),
 ('非局域量','string_order','String order','O_string'),
 ('非局域量','string_order_Id','String Id / parity','O_Id'),
 ('非局域量','nonlocal_sum','Nonlocal sum（实部）','Re sum O_i'),
 ('能量与密度','energy','总能量','E'),('能量与密度','energy_per_site','每格点能量','E/L'),
 ('能量与密度','energy_variance_raw','能量方差（含相消舍入）','Var(H)'),
 ('能量与密度','total_number','总粒子数','N'),
 ('能量与密度','bulk_density','体区平均密度','mean n'),
 ('能量与密度','center_density','中心两格点平均密度','n_center'),
 ('能量与密度','edge_density','端点平均密度','n_edge'),
 ('能量与密度','bulk_local_variance','体区平均局域数涨落','mean Var(n)'),
 ('能量与密度','bulk_pair_occupation','体区平均成对占据','mean n(n-1)/2'),
 ('关联量','bulk_g1_nn','体区最近邻单粒子关联','mean Re <b_i† b_(i+1)>'),
 ('关联量','g1_bulk_endpoints','体区端点单粒子关联','Re <b_a† b_b>'),
 ('关联量','bulk_connected_nn','体区最近邻连接密度关联','mean C_nn'),
 ('关联量','connected_bulk_endpoints','体区端点连接密度关联','C_nn(a,b)'),
 ('数值诊断','max_highest_occupation_probability','最高占据概率（最大值）','max P(n=nmax)'),
 ('数值诊断','reflection_density_error','密度反射不对称程度','max |n_i-n_(L+1-i)|'),
 ('配对指标','fidelity','相对参考组的保真度','F'),
 ('配对指标','infidelity','不保真度','1-F'),
 ('配对指标','max_density_difference','最大逐格点密度差','max |delta n|')]
PAIR_KEYS = {'fidelity', 'infidelity', 'max_density_difference'}
class MappingConflict(ValueError): pass

def observable_specs(order=6):
    return BASE_OBS + [('普通累积量','kappa%d'%r,'%d 阶普通累积量'%r,'kappa%d'%r)
        for r in range(1,order+1)] + [('标准化累积量','standardized%d'%r,
        ('偏度' if r==3 else '超额峰度（不加3）' if r==4 else '%d 阶标准化累积量'%r),
        'kappa%d / kappa2^(%d/2)'%(r,r)) for r in range(3,order+1)]

def decimal_text(value):
    d=Decimal(str(value))
    if not d.is_finite(): raise ValueError('U 必须是有限数值')
    return format(d.normalize(),'f') if d else '0'

def grid(start,end,step):
    try: a,b,s=map(lambda v: Decimal(str(v)),(start,end,step))
    except InvalidOperation: raise ValueError('U 起点、终点、步长必须是数字')
    if not all(v.is_finite() for v in (a,b,s)) or s<=0 or b<a:
        raise ValueError('要求步长 > 0，终点 ≥ 起点，且所有值有限')
    count=(b-a)/s
    if count!=count.to_integral_value(): raise ValueError('终点与步长不整除，请调整终点或步长')
    if count>10000: raise ValueError('单组最多支持 10001 个 U 点')
    return [decimal_text(a+i*s) for i in range(int(count)+1)]

def new_group(index=0):
    return dict(id=uuid.uuid4().hex,name='波函数组 '+chr(65+index%26),
                color=COLORS[index%len(COLORS)],style=['-','--',':','-.'][index%4],notes='',grid=None,
                params={},files={},pending=[])

def new_project():
    groups=[new_group(0),new_group(1)]
    return dict(version=VERSION,title='新的波函数对比',grid=['1','2','0.5'],
                common=copy.deepcopy(PARAMS),groups=groups,reference=groups[0]['id'],
                observables=['half_chain_entropy','entanglement_gap','string_order','kappa2','fidelity'])

def group_grid(project,group): return grid(*(group.get('grid') or project['grid']))
def effective_params(project,group):
    p=copy.deepcopy(project['common']);p.update(group.get('params',{}));return p

def update_grid(project,values,group=None):
    new=grid(*values)
    affected=[group] if group else [g for g in project['groups'] if not g.get('grid')]
    for g in affected:
        g['pending'] += [p for u,p in g['files'].items() if u not in new]
        g['files']={u:p for u,p in g['files'].items() if u in new}
    if group: group['grid']=list(values)
    else: project['grid']=list(values)

def assign_files(project,group,paths,start=None):
    paths=[os.path.abspath(p) for p in paths]
    if not paths: return
    if any(not os.path.isfile(p) or Path(p).suffix.lower() not in ('.h5','.hdf5') for p in paths):
        raise ValueError('仅接受存在的 .h5 或 .hdf5 文件')
    us=group_grid(project,group)
    if start is not None: us=us[us.index(decimal_text(start)):]
    slots=[u for u in us if not group['files'].get(u)]
    if len(paths)>len(slots): raise ValueError('剩余 %d 个空位，拖入 %d 个文件；请调整网格或清空位置'%(len(slots),len(paths)))
    for u,p in zip(slots,paths): group['files'][u]=p

def replace_file(group,u,path):
    if not os.path.isfile(path) or Path(path).suffix.lower() not in ('.h5','.hdf5'): raise ValueError('仅接受存在的 .h5 或 .hdf5 文件')
    old=group['files'].get(u)
    if old and old!=path: group['pending'].append(old)
    group['files'][u]=os.path.abspath(path)

def validate_project(p):
    if p.get('version')!=VERSION: raise ValueError('方案版本不支持')
    grid(*p['grid'])
    if not p.get('groups'): raise ValueError('方案至少需要一组波函数')
    ids=[g['id'] for g in p['groups']]
    if any(not isinstance(i,str) or not re.fullmatch(r'[A-Za-z0-9_-]{1,64}',i) for i in ids):raise ValueError('组标识格式无效')
    if len(ids)!=len(set(ids)): raise ValueError('组标识重复')
    if p['reference'] not in ids: raise ValueError('参考组不存在')
    p['common']={**PARAMS,**p['common']}
    known={x[1] for x in observable_specs(int(p['common']['max_order']))}
    if not set(p['observables'])<=known: raise ValueError('方案含未知物理量')
    for g in p['groups']:
        if not isinstance(g['name'],str) or not g['name'].strip():raise ValueError('组名称不能为空')
        if g['style'] not in ('-','--',':','-.'):raise ValueError('线型无效')
        if not re.fullmatch(r'#[0-9a-fA-F]{6}',g['color']):raise ValueError('颜色须采用 #RRGGBB 格式')
        us=group_grid(p,g)
        if not set(g['files'])<=set(us): raise ValueError('方案中的文件对应关系超出 U 网格')
        g.setdefault('pending',[]);g.setdefault('params',{});g.setdefault('notes','')
    return p

def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
    os.replace(str(tmp),str(path))

def check_file(path,u,params):
    import h5py
    if not os.path.isfile(path): raise ValueError('文件不存在')
    with h5py.File(path,'r') as f:
        if params['dataset'] not in f: raise ValueError('没有数据集 '+params['dataset'])
        for field,expected in [('U',float(u)),('t1',params['t1']),('ratio',params['ratio']),
                               ('dimerization',params['dimer'])]:
            if field in f:
                v=float(f[field][()])
                if not math.isfinite(v) or abs(v-expected)>1e-8*max(1,abs(expected)):
                    raise MappingConflict('%s 冲突：手动=%s，文件内部=%s；请修正'%(field,expected,v))
        for field,key in [('L','L'),('particle_number','N')]:
            enabled=params[key]>0 if key=='L' else params[key]>=0
            if enabled and field in f and abs(float(f[field][()])-params[key])>1e-8:
                raise MappingConflict(field+' 与物理设置不一致')
    return '内部参数检查通过（文件名未参与匹配）'

def toml_value(v):
    if isinstance(v,bool): return 'true' if v else 'false'
    if isinstance(v,str): return json.dumps(v,ensure_ascii=False)
    if isinstance(v,(int,float)): return repr(v)
    if isinstance(v,list): return '['+', '.join(toml_value(x) for x in v)+']'
    raise TypeError(type(v))

def write_job(project,run_dir,force=False):
    run_dir=Path(run_dir);run_dir.mkdir(parents=True,exist_ok=True)
    points=[]
    for g in project['groups']:
        params=effective_params(project,g)
        for u in group_grid(project,g):
            path=g['files'].get(u)
            if path:
                try:check_file(path,u,params)
                except MappingConflict:raise
                except Exception:pass  # Worker records unreadable points while other points continue.
                points.append(dict(group_id=g['id'],U=u,path=os.path.abspath(path),**params))
    if not points: raise ValueError('请先导入波函数')
    if not project['observables']: raise ValueError('请至少选择一个物理量')
    points.sort(key=lambda p:(Decimal(p['U']),p['group_id']!=project['reference']))
    fields=dict(version=1,output_dir=os.path.abspath(str(run_dir)),
                cache_dir=str(ROOT/'data'/'cache'),reference=project['reference'],
                observables=project['observables'],force=force)
    lines=[k+' = '+toml_value(v) for k,v in fields.items()]
    for p in points:
        lines.append('\n[[points]]')
        lines += [k+' = '+toml_value(v) for k,v in p.items()]
    path=run_dir/'job.toml';path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    atomic_json(run_dir/'方案快照.json',project)
    return str(path),points

def read_results(path):
    if not Path(path).exists(): return []
    with open(str(path),encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))

def example_project():
    p=new_project();p['title']='L20 MTBL 与 DMRG 示例模板'
    for g,name in zip(p['groups'],['MTBL','DMRG']):
        g['name']=name;g['notes']='U=1、1.5、2；请按顺序导入自己的波函数'
    return p
