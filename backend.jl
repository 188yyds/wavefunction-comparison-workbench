using ITensors, HDF5, CSV, DataFrames, LinearAlgebra, Statistics, SHA, TOML, Base64
if !isdefined(ITensors,:MPS)
    @eval using ITensorMPS
end
include(joinpath(@__DIR__,"core","体区间累积量核心.jl"))
include(joinpath(@__DIR__,"core","scaling物理量定义.jl"))
const B=BulkNumberCumulants
const SD=ScalingDefinitions
const VERSION_HASH=bytes2hex(sha256(vcat(read(@__FILE__),read(joinpath(@__DIR__,"core","体区间累积量核心.jl")),read(joinpath(@__DIR__,"core","scaling物理量定义.jl")))))
const PAIR=Set(["fidelity","infidelity","max_density_difference"])
const NBAR=Ref(0.5)
function ITensors.op(::OpName"WCdn",::SiteType"Boson",s::Index)
    op("N",s)-NBAR[]*op("Id",s)
end
function ITensors.op(::OpName"WCphase",::SiteType"Boson",s::Index)
    exp(im*pi*(op("N",s)-NBAR[]*op("Id",s)))
end
function ITensors.op(::OpName"WChigh",::SiteType"Boson",d::Int)
    v=zeros(d,d);v[end,end]=1;v
end
function event(kind,g="",u="",key="",message="")
    println("WC\t",kind,"\t",g,"\t",u,"\t",key,"\t",base64encode(string(message)));flush(stdout)
end
function replace_file_atomic(tmp,path)
    if Sys.iswindows()
        source=transcode(UInt16,tmp*"\0");target=transcode(UInt16,path*"\0")
        for attempt in 1:100
            ok=GC.@preserve source target ccall((:MoveFileExW,"kernel32"),Cint,
                (Ptr{UInt16},Ptr{UInt16},UInt32),source,target,UInt32(1))
            ok!=0 && return
            sleep(.02)
        end
        error("结果文件持续被占用，无法更新：$(path)；请关闭占用它的程序")
    else
        Base.Filesystem.rename(tmp,path)
    end
end
function atomic_toml(path,data)
    mkpath(dirname(path));tmp=path*".$(getpid()).tmp"
    open(tmp,"w") do io;TOML.print(io,data);end
    replace_file_atomic(tmp,path)
end
function signature(p)
    statp=stat(p["path"])
    cfg=sort([(k,string(v)) for (k,v) in p if k!="group_id"];by=first)
    bytes2hex(sha256(string(VERSION_HASH,cfg,statp.size,statp.mtime)))
end
function cachepath(job,p)
    joinpath(job["cache_dir"],signature(p)*".toml")
end
function getcache(job,p)
    path=cachepath(job,p)
    if isfile(path)
        try;return TOML.parsefile(path);catch;end
    end
    Dict{String,Any}("values"=>Dict{String,Any}(),"info"=>Dict{String,Any}())
end
function geometry(p,L)
    automatic=B.subsystem_bounds(L;require_even_odd=false,require_bulk=true)
    a,b=(p["bulk_left"]==0 && p["bulk_right"]==0) ? automatic : (p["bulk_left"],p["bulk_right"])
    1<a<=b<L || error("体区间要求 1 < left ≤ right < L")
    st,en=(p["string_left"]==0 && p["string_right"]==0) ? SD.direct_bulk_string_pair(L) : (p["string_left"],p["string_right"])
    1<=st<en<=L || error("String 端点要求 1 ≤ left < right ≤ L")
    cut=p["cut"]==0 ? clamp(div(L,2),1,L-1) : p["cut"]
    1<=cut<L || error("切口要求 1 ≤ cut < L")
    a,b,st,en,cut
end
function load_validate(p)
    psi,meta=h5open(p["path"],"r") do f
        haskey(f,p["dataset"]) || error("没有指定的 MPS 数据集")
        m=Dict(k=>read(f,k) for k in keys(f) if k!=p["dataset"] && f[k] isa HDF5.Dataset)
        read(f,p["dataset"],MPS),m
    end
    for (k,v) in [("U",parse(Float64,p["U"])),("ratio",p["ratio"]),("t1",p["t1"]),("dimerization",p["dimer"])]
        if haskey(meta,k);isapprox(Float64(meta[k]),v;atol=1e-8,rtol=1e-8) || error("$k 与文件内部参数冲突，须修正对应关系");end
    end
    L=length(psi);L>=4 || error("此版本的体区间物理量要求 L ≥ 4")
    all(s->hastags(s,"Boson"),siteinds(psi)) || error("仅支持 Boson MPS")
    p["L"]>0 && p["L"]!=L && error("设置 L 与 MPS 不一致")
    raw=norm(psi);isfinite(raw) && raw>0 || error("波函数范数无效")
    normalize!(psi);n=Float64.(real.(expect(psi,"N")));N=sum(n)
    haskey(meta,"L") && Int(meta["L"])!=L && error("文件记录的 L 与 MPS 链长不一致")
    haskey(meta,"particle_number") && !isapprox(N,meta["particle_number"];atol=1e-7,rtol=1e-9) && error("文件记录的粒子数与 MPS 不一致")
    haskey(meta,"boson_dim") && !all(==(Int(meta["boson_dim"])),dim.(siteinds(psi))) && error("文件记录的局域维数与 MPS 不一致")
    p["N"]>=0 && !isapprox(N,p["N"];atol=1e-7,rtol=1e-9) && error("设置 N 与 MPS 不一致")
    a,b,st,en,cut=geometry(p,L)
    strong(x)=isodd(x)==(p["dimer"]==1) ? "strong" : "weak"
    info=Dict{String,Any}("L"=>L,"N"=>N,"dimensions"=>dim.(siteinds(psi)),"density"=>n,
        "raw_norm"=>raw,"maxlinkdim"=>maxlinkdim(psi),"has_qns"=>hasqns(psi),
        "bulk_left"=>a,"bulk_right"=>b,"string_left"=>st,"string_right"=>en,"cut"=>cut,
        "center_cut_type"=>strong(cut),"block_left_cut_type"=>strong(a-1),"block_right_cut_type"=>strong(b),
        "saved_accepted"=>string(get(meta,"accepted","not_recorded")),
        "saved_residual"=>Float64(get(meta,"energy_residual_norm",NaN)),
        "source_path"=>p["path"],"source_size"=>filesize(p["path"]),"source_mtime"=>mtime(p["path"]))
    psi,info
end
function hamiltonian(sites,p)
    os=OpSum();U=parse(Float64,p["U"])
    for i in 1:length(sites)-1
        t=p["t1"]*(isodd(i)==(p["dimer"]==1) ? 1.0 : p["ratio"])
        os+=-t,"Adag",i,"A",i+1;os+=-t,"Adag",i+1,"A",i
    end
    for i in eachindex(sites);os+=U/2,"N",i,"N",i;os+=-U/2,"N",i;end
    MPO(os,sites)
end
function center_svd(psi,cut)
    orthogonalize!(psi,cut)
    left=cut==1 ? (siteind(psi,cut),) : (linkind(psi,cut-1),siteind(psi,cut))
    _,S,_=svd(psi[cut],left);S
end
function projected_exterior(psi,a,b,cutoff)
    state=deepcopy(psi);losses=Float64[]
    for cut in (a-1,b)
        orthogonalize!(state,cut)
        left=cut==1 ? (siteind(state,cut),) : (linkind(state,cut-1),siteind(state,cut))
        U,S,V,spec=svd(state[cut],left;cutoff=cutoff)
        push!(losses,Float64(spec.truncerr));state[cut]=U;state[cut+1]=S*V*state[cut+1];normalize!(state)
    end
    p=clamp(1-prod(1 .- losses),0.0,1.0)
    h=p==0 ? 0.0 : -p*log(p)-(1-p)*log1p(-p)
    state,p,h+p*sum(log(dim(siteind(psi,i))) for i in a:b)
end
# Count actual QN output blocks before allocating a large interval contraction.
function contraction_bytes(A,C)
    ai=collect(inds(A));ci=collect(inds(C))
    shared=[i for i in eachindex(ai) if any(==(ai[i]),ci)]
    cshared=[findfirst(==(ai[i]),ci) for i in shared]
    ae=[i for i in eachindex(ai) if !(i in shared)];ce=[i for i in eachindex(ci) if !(i in cshared)]
    if !hasqns(A) || !hasqns(C)
        return 8.0*prod(Float64(dim(ai[i])) for i in ae)*prod(Float64(dim(ci[i])) for i in ce)
    end
    cb=Dict{Tuple,Vector{Tuple}}()
    for block in ITensors.eachnzblock(C)
        t=Tuple(block);k=Tuple(t[i] for i in cshared);push!(get!(cb,k,Tuple[]),t)
    end
    outputs=Set{Tuple}();elements=0.0
    for block in ITensors.eachnzblock(A)
        t=Tuple(block);k=Tuple(t[i] for i in shared)
        for ct in get(cb,k,Tuple[])
            ot=Tuple(vcat([t[i] for i in ae],[ct[i] for i in ce]))
            if !(ot in outputs)
                push!(outputs,ot)
                elements+=prod(Float64(ITensors.blockdim(ai[i],t[i])) for i in ae)*prod(Float64(ITensors.blockdim(ci[i],ct[i])) for i in ce)
            end
        end
    end
    elements*max(sizeof(eltype(A)),sizeof(eltype(C)))
end
function block_entropy(state,a,b,memory_mb)
    orthogonalize!(state,b)
    left=commonind(state[a],state[a-1]);right=commonind(state[b],state[b+1])
    bra=prime(dag(state),linkinds(state));limit=memory_mb*1024.0^2
    dl=dim(left);dr=dim(right);T=eltype(state[a])
    3.0*(dl*dr)^2*sizeof(T)>limit && error("块熵最终密度矩阵超过内存限额；未放宽精度")
    matrix=zeros(T,dl*dr,dl*dr);pieces=[];offset=0
    spaces=hasqns(left) ? space(left) : [nothing=>dl]
    for item in spaces
        q,d=item
        for start in 1:8:d
            count=min(8,d-start+1)
            small=hasqns(left) ? Index([q=>count],"EntropyLeft";dir=dir(left)) : Index(count,"EntropyLeft")
            W=ITensor(Float64,small,dag(left))
            for k in 1:count;W[small=>k,dag(left)=>offset+start+k-1]=1.0;end
            push!(pieces,(first=offset+start,count=count,index=small,ket=state[a]*W,bra=bra[a]*prime(dag(W))))
        end
        offset+=d
    end
    # Split exterior Schmidt sectors exactly; no additional truncation of the state.
    for (li,lpart) in enumerate(pieces),(ri,rpart) in enumerate(pieces)
        ri<li && continue
        rho=lpart.ket*rpart.bra
        for j in a+1:b
            next=nothing;s=siteind(state,j)
            for n in 1:dim(s)
                ket=state[j]*onehot(dag(s)=>n);norm(ket)==0 && continue
                bra_n=prime(dag(ket),linkinds(state))
                contraction_bytes(rho,ket)*3>limit && error("块熵分块收缩仍超过内存限额；未放宽精度")
                term=rho*ket
                contraction_bytes(term,bra_n)*3>limit && error("块熵分块收缩仍超过内存限额；未放宽精度")
                term=term*bra_n;next=isnothing(next) ? term : next+term
            end
            isnothing(next) && error("块熵收缩得到空环境")
            rho=next
        end
        data=Array(rho,lpart.index,right,prime(dag(rpart.index)),prime(dag(right)))
        rows=[lpart.first+k-1+(r-1)*dl for r in 1:dr for k in 1:lpart.count]
        cols=[rpart.first+k-1+(r-1)*dl for r in 1:dr for k in 1:rpart.count]
        block=reshape(data,length(rows),length(cols));matrix[rows,cols]=block
        ri!=li && (matrix[cols,rows]=adjoint(block))
    end
    spectrum=eigvals(Hermitian(matrix))
    probs=[abs(x) for x in spectrum if abs(x)>1e-12];isempty(probs) && error("块熵谱为空")
    probs/=sum(probs);-sum(p*log(p) for p in probs)
end
function entropy_checked(psi,p,info)
    a=info["bulk_left"];b=info["bulk_right"]
    coarse=p["entropy_coarse"];fine=p["entropy_fine"]
    0<fine<coarse<1 || error("要求 0 < 精细阈值 < 粗阈值 < 1")
    ss=Float64[];loss=Float64[];bounds=Float64[]
    for cutoff in (coarse,fine)
        state,pdrop,bound=projected_exterior(psi,a,b,cutoff)
        push!(ss,block_entropy(state,a,b,p["memory_mb"]));push!(loss,pdrop);push!(bounds,bound);GC.gc()
    end
    drift=abs(ss[2]-ss[1]);drift<=sum(bounds)+1e-9 || error("块熵两档检查未通过")
    info["entropy_discarded_weight"]=loss[2];info["entropy_error_bound"]=bounds[2]
    info["entropy_coarse_value"]=ss[1];info["entropy_precision_difference"]=drift
    ss[2]
end
function align_fidelity(psi,other)
    length(psi)==length(other) || error("L 不同")
    hasqns(psi)==hasqns(other) || error("QN 基底结构不同，暂不支持此保真度配对")
    x=deepcopy(psi);y=deepcopy(other)
    for i in 1:length(x)
        sx=siteind(x,i);sy=siteind(y,i);target=dim(sx)>=dim(sy) ? sx : sy
        vals=B.number_values(target)
        for state in (x,y)
            old=siteind(state,i);old==target && continue
            W=ITensor(Float64,target,dag(old))
            for (j,n) in enumerate(B.number_values(old))
                k=findfirst(==(n),vals);isnothing(k) && error("Fock 基底无法对齐")
                W[target=>k,dag(old)=>j]=1.0
            end
            state[i]=W*state[i]
        end
    end
    clamp(abs2(inner(x,y))/(real(inner(x,x))*real(inner(y,y))),0.0,1.0)
end
function pairpath(job,p,r)
    joinpath(job["cache_dir"],"pair_"*bytes2hex(sha256(signature(p)*signature(r)))*".toml")
end
function stored_pair(job,p,r)
    path=try;pairpath(job,p,r);catch;return nothing;end
    if !job["force"] && isfile(path)
        try
            values=TOML.parsefile(path)
            all(haskey(values,k) for k in job["observables"] if k in PAIR) && return values
        catch
        end
    end
    nothing
end

function main(jobfile,phase)
    job=TOML.parsefile(jobfile);points=job["points"];obs=String.(job["observables"])
    BLAS.set_num_threads(Int(first(points)["blas_threads"]))
    output=job["output_dir"];mkpath(output);mkpath(job["cache_dir"])
    resultpath=joinpath(output,"results.csv")
    rows=NamedTuple[]
    if phase=="center" && isfile(resultpath)
        for row in eachrow(CSV.read(resultpath,DataFrame;types=String))
            push!(rows,(group_id=String(row.group_id),U=String(row.U),observable=String(row.observable),
                value=String(row.value),status=String(row.status),message=ismissing(row.message) ? "" : String(row.message)))
        end
    end
    function output_row(p,k,v,status,msg="")
        push!(rows,(group_id=String(p["group_id"]),U=String(p["U"]),observable=k,value=string(v),status=status,message=string(msg)))
        tmp=resultpath*".$(getpid()).tmp"
        CSV.write(tmp,DataFrame(rows));replace_file_atomic(tmp,resultpath)
        event(status,p["group_id"],p["U"],k,msg)
    end
    refs=Dict(p["U"]=>p for p in points if p["group_id"]==job["reference"])
    infos=Dict{String,Any}();states=Dict{String,MPS}();activeu=""
    for p in points
        BLAS.set_num_threads(Int(p["blas_threads"]))
        u=p["U"];g=p["group_id"]
        u!=activeu && (empty!(states);empty!(infos);GC.gc();activeu=u)
        selected=phase=="center" ? filter(==("center_entropy"),obs) : filter(k->k!="center_entropy" && !(k in PAIR),obs)
        pairselected=phase=="center" ? String[] : filter(k->k in PAIR,obs)
        event("point",g,u,"",phase)
        try
            sig=signature(p);cachefile=joinpath(job["cache_dir"],sig*".toml")
            cache=getcache(job,p);vals=cache["values"];info=cache["info"];errors=Dict{String,String}()
            if job["force"]
                for k in selected;delete!(vals,k);end
            end
            needpair=(any(k->k in ("fidelity","infidelity"),pairselected) && haskey(refs,u) && any(q->q["U"]==u && q["group_id"]!=job["reference"] && isnothing(stored_pair(job,q,refs[u])),points))
            missingkeys=[k for k in selected if !haskey(vals,k)]
            psi=nothing
            if job["force"] || !isempty(missingkeys) || isempty(info) || needpair
                psi,info=load_validate(p);cache["info"]=info;states[g]=psi
            end
            infos[g]=info
            function family(keys,f)
                needed=[k for k in selected if k in keys && !haskey(vals,k)]
                isempty(needed) && return
                event("calculating",g,u,first(needed))
                try
                    got=f()
                    signature(p)==sig || error("计算期间源文件发生变化，请重新计算")
                    for (k,v) in got
                        if isfinite(v)
                            vals[k]=v
                        else
                            delete!(vals,k);errors[k]="$k 未得到有效有限值（检查谱秩、方差或数值稳定性）"
                        end
                    end
                    atomic_toml(cachefile,cache)
                    for k in needed
                        if haskey(vals,k);output_row(p,k,vals[k],"computed")
                        else;output_row(p,k,NaN,"failed",errors[k]);end
                    end
                catch err
                    for k in needed;errors[k]=sprint(showerror,err);output_row(p,k,NaN,"failed",errors[k]);end
                    GC.gc()
                end
            end
            n=info["density"];L=info["L"];a=info["bulk_left"];b=info["bulk_right"]
            family(["total_number","bulk_density","center_density","edge_density","reflection_density_error"],()->Dict(
                "total_number"=>sum(n),"bulk_density"=>mean(n[a:b]),
                "center_density"=>mean(n[clamp(div(L,2),1,L):clamp(div(L,2)+1,1,L)]),
                "edge_density"=>mean(n[[1,L]]),"reflection_density_error"=>maximum(abs.(n-reverse(n)))))
            family(["half_chain_entropy","entanglement_gap","ES1","ES2"],()->begin
                S=center_svd(psi,info["cut"]);gap,e1,e2=SD.entanglement_gap_from_singular_values(S)
                Dict("half_chain_entropy"=>SD.entropy_from_singular_values(S),"entanglement_gap"=>gap,"ES1"=>e1,"ES2"=>e2)
            end)
            ck=filter(k->startswith(k,"kappa") || startswith(k,"standardized"),selected)
            if any(k->!haskey(vals,k),ck)
                family(ck,()->begin
                    orders=[parse(Int,replace(k,r"\D"=>"")) for k in ck]
                    m=maximum(orders);cr=B.calculate_cumulants(psi,a,b;max_order=m,precision_mode=:checked)
                    got=Dict{String,Float64}()
                    info["cumulant_order_passed"]=collect(cr.order_passed)
                    info["cumulant_bits"]=cr.contraction_bits
                    info["cumulant_direction_difference"]=Float64.(cr.stability_difference)
                    for k in ck
                        r=parse(Int,replace(k,r"\D"=>""));required=startswith(k,"standardized") ? [2,r] : [r]
                        if !all(cr.order_passed[required])
                            got[k]=NaN;continue
                        end
                        v=cr.cumulants[r]
                        if startswith(k,"standardized")
                            if cr.cumulants[2]<=1e-10
                                got[k]=NaN;continue
                            end
                            v/=cr.cumulants[2]^(BigFloat(r)/2)
                        end
                        got[k]=Float64(v)
                    end
                    got
                end)
            end
            for k in ("string_order","string_order_Id")
                family([k],()->Dict(k=>-real(inner(psi',SD.string_mpo(info["string_left"],info["string_right"],siteinds(psi);endpoint_op=k=="string_order" ? "sigmaz_b" : "Id"),psi))))
            end
            family(["nonlocal_sum"],()->begin
                NBAR[]=p["nbar"];sites=siteinds(psi);total=0.0+0.0im
                for i in 1:L
                    os=OpSum();terms=("WCdn",i)
                    for j in i-1:-1:1;terms=(terms...,"WCphase",j);end
                    os+=terms;total+=inner(psi',MPO(os,sites),psi)
                end
                info["nonlocal_imag"]=imag(total);Dict("nonlocal_sum"=>real(total))
            end)
            family(["bulk_local_variance","bulk_pair_occupation"],()->begin
                n2=real.(expect(psi,"N * N"));Dict("bulk_local_variance"=>mean((n2-n.^2)[a:b]),"bulk_pair_occupation"=>mean((n2-n)[a:b])/2)
            end)
            family(["bulk_g1_nn","g1_bulk_endpoints"],()->begin
                c=real.(correlation_matrix(psi,"Adag","A"));Dict("bulk_g1_nn"=>mean([c[i,i+1] for i in a:b-1]),"g1_bulk_endpoints"=>c[a,b])
            end)
            family(["bulk_connected_nn","connected_bulk_endpoints"],()->begin
                c=real.(correlation_matrix(psi,"N","N"))-n*transpose(n)
                Dict("bulk_connected_nn"=>mean([c[i,i+1] for i in a:b-1]),"connected_bulk_endpoints"=>c[a,b])
            end)
            energykeys=["energy","energy_per_site","energy_variance_raw"]
            if any(k->k in missingkeys,energykeys)
                H=hamiltonian(siteinds(psi),p);E=real(inner(psi',H,psi))
                family(["energy","energy_per_site"],()->Dict("energy"=>E,"energy_per_site"=>E/L))
                family(["energy_variance_raw"],()->Dict("energy_variance_raw"=>real(inner(H,psi,H,psi))-E^2))
                info["energy_variance_roundoff_scale"]=eps(Float64)*max(E^2,1)
            end
            family(["max_highest_occupation_probability"],()->Dict("max_highest_occupation_probability"=>maximum(real.(expect(psi,"WChigh")))))
            family(["center_entropy"],()->Dict("center_entropy"=>entropy_checked(psi,p,info)))
            for k in selected
                if !(k in missingkeys) && haskey(vals,k);output_row(p,k,vals[k],"cached");end
            end
            signature(p)==sig || error("计算期间源文件变化，当前点结果无效")
            atomic_toml(cachefile,cache)
            atomic_toml(joinpath(output,"info_$(g)_$(u).toml"),info)
            if !isempty(pairselected)
                if !haskey(refs,u) || !haskey(infos,job["reference"])
                    for k in pairselected;output_row(p,k,NaN,"incompatible","参考组在此 U 没有有效状态");end
                else
                    r=refs[u];ri=infos[job["reference"]];pv=stored_pair(job,p,r);paircached=!isnothing(pv)
                    if isnothing(pv)
                        try
                            info["L"]==ri["L"] && isapprox(info["N"],ri["N"];atol=1e-7) || error("L 或粒子数扇区不同，配对指标不适用")
                            for (source,sourceinfo) in ((p,info),(r,ri))
                                sourceinfo["source_size"]==filesize(source["path"]) && sourceinfo["source_mtime"]==mtime(source["path"]) || error("配对期间源文件变化")
                            end
                            pairfile=pairpath(job,p,r)
                            pv=Dict{String,Any}("max_density_difference"=>maximum(abs.(n-ri["density"])))
                            if any(k->k in ("fidelity","infidelity"),pairselected)
                                f=g==job["reference"] ? 1.0 : align_fidelity(states[g],states[job["reference"]])
                                pv["fidelity"]=f;pv["infidelity"]=1-f
                            end
                            for (source,sourceinfo) in ((p,info),(r,ri))
                                sourceinfo["source_size"]==filesize(source["path"]) && sourceinfo["source_mtime"]==mtime(source["path"]) || error("配对期间源文件变化")
                            end
                            atomic_toml(pairfile,pv)
                        catch err
                            for k in pairselected;output_row(p,k,NaN,"incompatible",sprint(showerror,err));end
                        end
                    end
                    if !isnothing(pv)
                        for k in pairselected;output_row(p,k,pv[k],paircached ? "cached" : "computed");end
                    end
                end
            end
        catch err
            for k in vcat(selected,pairselected);output_row(p,k,NaN,"failed",sprint(showerror,err));end
        end
        event("point_done",g,u);GC.gc()
    end
    event("finished","","","",phase)
end
if abspath(PROGRAM_FILE)==@__FILE__
    main(ARGS[1],length(ARGS)>1 ? ARGS[2] : "normal")
end
