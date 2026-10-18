import os, sys, json, subprocess, ctypes
from pathlib import Path
root=Path(os.path.abspath(__file__)).parent
try:
    target=sys.executable
    config=root/'data'/'runtime.json'
    if config.exists():target=json.loads(config.read_text(encoding='utf-8'))['python']
    hidden=target.replace('python.exe','pythonw.exe')
    if os.path.isfile(hidden):target=hidden
    (root/'data').mkdir(exist_ok=True)
    log=open(str(root/'data'/'启动日志.txt'),'w',encoding='utf-8')
    environment=os.environ.copy();environment['PYTHONIOENCODING']='utf-8';environment['PYTHONDONTWRITEBYTECODE']='1'
    process=subprocess.Popen([target,str(root/'app.py')],cwd=str(root),stdout=log,stderr=log,
                             env=environment,creationflags=0x08000000)
    try:
        code=process.wait(timeout=3)
        if code:
            ctypes.windll.user32.MessageBoxW(None,'启动失败，请查看 data/启动日志.txt。删除 data/runtime.json 可恢复默认运行环境。','波函数对比工作台',16)
    except subprocess.TimeoutExpired:pass
except Exception as e:
    ctypes.windll.user32.MessageBoxW(None,str(e),'波函数对比工作台：启动错误',16)
