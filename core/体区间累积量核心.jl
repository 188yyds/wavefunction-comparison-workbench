# 共用函数；仅载入定义，不启动计算、画图或自检。
module BulkNumberCumulants

using ITensors
if !isdefined(ITensors, :MPS)
    @eval using ITensorMPS
end
using HDF5, JLD2, LinearAlgebra, Printf, Dates

const METHOD_VERSION = "bulk_number_taylor_v1"
const OBSERVABLE = "N_A = sum(N_j, j=left_site:right_site); ordinary cumulants, no length normalization"

# 扩展名沿用 .h5，但格式为 JLD2。直接指定读取器，避免 FileIO 按 .h5 错选 HDF5。
load_result(path) = JLD2.jldopen(path,"r") do f
    Dict(key => f[key] for key in keys(f))
end

function subsystem_bounds(L::Int; region=nothing, require_even_odd=true, require_bulk=true)
    if isnothing(region)
        q = div(L, 4)
        a, b = isodd(q) ? (q + 1, 3q) : (q + 2, 3q + 1)
    else
        length(region) == 2 || error("region 必须是 (起点,终点)")
        all(x -> x isa Integer, region) || error("区间端点必须是整数")
        a, b = Int.(region)
    end
    1 <= a <= b <= L || error("非法区间 $(a):$(b)，L=$L")
    require_bulk && !(1 < a <= b < L) && error("体区间须满足 1<a<=b<L")
    require_even_odd && !(iseven(a) && isodd(b)) && error("要求偶数起点、奇数终点，实际为 $a:$b")
    return a, b
end

function wavefunction_path(cfg, L, U, ratio)
    joinpath(cfg.wavefunction_dir, @sprintf(
        "dimer_%i_j_%.3f_U_%.3f_L_%i_ratio_%.3f_filling_%s_nsweep_%i_mps.h5",
        cfg.dimerization, cfg.t1, U, L, ratio, cfg.filling, cfg.nsweep))
end

function result_path(cfg, L, U, ratio, a, b)
    joinpath(cfg.high_order_data_dir, "cumulants", "order1-$(cfg.max_order)", @sprintf(
        "cumulants_ratio_%.3f_U_%.3f_dimer_%i_j_%.3f_L_%i_filling_%s_nsweep_%i_A_%i-%i.h5",
        ratio, U, cfg.dimerization, cfg.t1, L, cfg.filling, cfg.nsweep, a, b))
end

function number_values(s)
    # 读取真正的局域 N 算符，而非把数组下标或 QN 标签猜作粒子数。
    mat = Matrix(op("N", s), prime(s), dag(s))
    norm(mat - Diagonal(diag(mat))) <= 1e-13 || error("此方法要求局域 N 在站点基底中对角")
    maximum(abs, imag.(diag(mat))) <= 1e-13 || error("局域粒子数必须为实数")
    v = real.(diag(mat))
    all(x -> x >= 0 && abs(x-round(x)) <= 1e-13, v) || error("需要非负整数粒子数谱")
    return Int.(round.(v))
end

tensor_type(psi, ::Type{R}) where {R} = any(j -> eltype(psi[j]) <: Complex, 1:length(psi)) ? Complex{R} : R
matrix_as(::Type{T}, A) where {T} = eltype(A)==T ? A : T.(A)

function checked_real(z; tolerance=1e-10)
    isfinite(z) || throw(OverflowError("收缩产生非有限数；请提高精度/调整 moment_scale"))
    abs(imag(z)) <= tolerance * max(one(real(z)), abs(real(z))) || error("收缩存在不可忽略虚部: $z")
    return real(z)
end

function link_groups(i)
    if !hasqns(i)
        return (ids=[1], ranges=[1:dim(i)], dims=[dim(i)])
    end
    ids, ranges, dims = Int[], UnitRange{Int}[], Int[]
    labels = QN[]
    for k in 1:ITensors.nblocks(i)
        q = ITensors.qn(i,k)
        g = findfirst(==(q),labels)
        if isnothing(g)
            push!(labels,q)
            push!(dims,0)
            g = length(dims)
        end
        d = ITensors.blockdim(i,k)
        push!(ids,g)
        push!(ranges,dims[g]+1:dims[g]+d)
        dims[g] += d
    end
    return (ids=ids,ranges=ranges,dims=dims)
end

"""仅重排/复制已有非零块，不做 SVD 或截断；重复 QN 标签先合并以保留相干项。"""
function prepare_matrices(psi)
    network = []
    for j in 1:length(psi)
        A = psi[j]
        s = siteind(psi,j)
        left = j==1 ? (hasqns(A) ? Index([QN()=>1],"LeftBoundary") : Index(1,"LeftBoundary")) : commonind(A,psi[j-1])
        right = j==length(psi) ? (hasqns(A) ? Index([QN()=>1],"RightBoundary") : Index(1,"RightBoundary")) : commonind(A,psi[j+1])
        j==1 && (A *= onehot(left=>1))
        j==length(psi) && (A *= onehot(right=>1))
        A = permute(A,left,s,right) # 只对输入的 Float64/ComplexF64 张量重排
        lg,rg = link_groups(left),link_groups(right)
        T = eltype(A)
        matrices = Dict{NTuple{3,Int},Matrix{T}}()
        if hasqns(A)
            offset = cumsum(vcat(0,[ITensors.blockdim(s,k) for k in 1:ITensors.nblocks(s)]))
            for block in ITensors.eachnzblock(A)
                kl,ks,kr = Tuple(block)
                gl,gr = lg.ids[kl],rg.ids[kr]
                data = Array(ITensors.NDTensors.array(ITensors.NDTensors.blockview(ITensors.tensor(A),block)))
                for nlocal in axes(data,2)
                    n = offset[ks]+nlocal
                    mat = get!(matrices,(n,gl,gr)) do
                        zeros(T,lg.dims[gl],rg.dims[gr])
                    end
                    mat[lg.ranges[kl],rg.ranges[kr]] .= @view data[:,nlocal,:]
                end
            end
        else
            data = Array(A,left,s,right)
            for n in 1:dim(s)
                matrices[(n,1,1)] = Matrix(@view data[:,n,:])
            end
        end
        records = [(n=key[1],left=key[2],right=key[3],matrix=mat) for (key,mat) in sort!(collect(matrices);by=first)]
        # 固定 QN 的张量在给定 n 下必须一对一连接左右 QN 扇区。
        # 否则只存对角环境块会漏掉相干项，必须拒绝而非静默计算。
        for n in 1:dim(s)
            subset = filter(x -> x.n==n,records)
            length(unique(x.left for x in subset))==length(subset) || error("局域 QN 映射不是一对一")
            length(unique(x.right for x in subset))==length(subset) || error("局域 QN 映射不是一对一")
        end
        push!(network,(records=records,left_dims=lg.dims,right_dims=rg.dims,values=number_values(s)))
    end
    return network
end

function matrix_environment(E,site,::Type{T};reverse=false) where {T}
    dims = reverse ? site.left_dims : site.right_dims
    out = [zeros(T,d,d) for d in dims]
    for record in site.records
        A = matrix_as(T,record.matrix)
        if reverse
            out[record.left] .+= (A * E[record.right]) * adjoint(A)
        else
            out[record.right] .+= (adjoint(A) * E[record.left]) * A
        end
    end
    return out
end

"""
精确到舍入误差的截断幂级数收缩。没有 SVD、MPO 压缩、采样或数值差分。
始终从原始张量显式收缩区间外的左右环境；高精度路径也重算这些环境。
E[p+1] 记录 t^p 系数。局域权重 w[r+1]=(n-c_j)^r/(scale^r*r!)。
对每个 p,n 仅计算一次 A_n E_p A_n†，复用于全部 r>=p。
"""
function taylor_contract(psi, a, b, M, centers; realtype=Float64, scale=1.0, direction=:left,network=nothing)
    1 <= a <= b <= length(psi) || error("非法区间")
    M >= 1 || error("max_order 必须 >= 1")
    length(centers) == b-a+1 || error("中心数组长度错误")
    direction in (:left, :right) || error("direction 必须为 :left 或 :right")
    R = realtype
    T = tensor_type(psi, R)
    sc = R(scale)
    isfinite(sc) && sc > 0 || error("moment_scale 必须为有限正数")
    net = isnothing(network) ? prepare_matrices(psi) : network
    # 不使用正交规范假设；高精度路径从输入数据重新收缩全部环境。
    left, right = [ones(T,1,1)], [ones(T,1,1)]
    for j in 1:a-1
        left = matrix_environment(left,net[j],T)
    end
    for j in length(psi):-1:b+1
        right = matrix_environment(right,net[j],T;reverse=true)
    end
    start, finish = direction == :left ? (left, right) : (right, left)
    E = [k == 0 ? start : [zeros(T,size(mat)) for mat in start] for k in 0:M]
    sequence = direction == :left ? (a:b) : (b:-1:a)
    for j in sequence
        site = net[j]
        dims = direction==:left ? site.right_dims : site.left_dims
        nextE = [[zeros(T,d,d) for d in dims] for k in 0:M]
        for record in site.records
            A = matrix_as(T,record.matrix)
            w = ones(R, M+1)
            x = (R(site.values[record.n]) - R(centers[j-a+1])) / sc
            for k in 1:M
                w[k+1] = w[k] * x / R(k)
            end
            for p in 0:M
                if direction==:left
                    F = (adjoint(A) * E[p+1][record.left]) * A
                    destination = record.right
                else
                    F = (A * E[p+1][record.right]) * adjoint(A)
                    destination = record.left
                end
                for k in p:M
                    nextE[k+1][destination] .+= w[k-p+1] .* F
                end
            end
        end
        E = nextE
    end
    # tr(E*R)，不需要真的做一次矩阵乘法。这里不取共轭，因为 bra 已在环境中。
    z = [checked_real(sum((sum(e[g] .* transpose(finish[g])) for g in eachindex(e));init=zero(T))) for e in E]
    z[1] > 0 || error("波函数范数不为正")
    return (coefficients=z ./ z[1], norm2=z[1])
end

"""log Z 的代数递推；输出普通累积量。即使中心并非精确均值也成立。"""
function coefficients_to_cumulants(coefficients, centers, scale; bits=256)
    setprecision(BigFloat, bits) do
        a = BigFloat.(coefficients)
        a ./= a[1]
        M = length(a)-1
        b = zeros(BigFloat, M)
        condition = zeros(BigFloat, M)
        kappas = zeros(BigFloat, M)
        sc = BigFloat(scale)
        factorial_r = one(BigFloat)
        scale_r = one(BigFloat)
        for r in 1:M
            terms = [BigFloat(k) * b[k] * a[r-k+1] / BigFloat(r) for k in 1:r-1]
            b[r] = a[r+1] - sum(terms; init=zero(BigFloat))
            magnitude = abs(a[r+1]) + sum(abs, terms; init=zero(BigFloat))
            condition[r] = iszero(b[r]) ? (iszero(magnitude) ? one(BigFloat) : BigFloat(Inf)) : magnitude/abs(b[r])
            factorial_r *= r
            scale_r *= sc
            kappas[r] = factorial_r * scale_r * b[r]
        end
        kappas[1] += sum(BigFloat.(centers))
        all(isfinite, kappas) || throw(OverflowError("累积量非有限"))
        return (cumulants=kappas, log_coefficients=b, cancellation=condition)
    end
end

function evaluate(psi, a, b, M, centers; bits=53, direction=:left, scale=1.0,network=nothing)
    R = bits == 53 ? Float64 : BigFloat
    setprecision(BigFloat, max(256, bits)) do
        raw = if bits == 53
            taylor_contract(psi, a,b,M,centers; realtype=R,scale=scale,direction=direction,network=network)
        else
            setprecision(BigFloat, bits) do
                taylor_contract(psi, a,b,M,centers; realtype=R,scale=scale,direction=direction,network=network)
            end
        end
        post = coefficients_to_cumulants(raw.coefficients, centers, scale; bits=max(256, 2bits))
        return merge(raw, post, (contraction_bits=bits, direction=direction))
    end
end

function stability(first, second; atol, rtol)
    delta = abs.(first .- second)
    tolerance = BigFloat(atol) .+ BigFloat(rtol) .* max.(abs.(first), abs.(second))
    return (delta=delta, passed=delta .<= tolerance)
end

function calculate_cumulants(psi, a, b; max_order=6, precision_mode=:checked,
        precision_bits=128, confirmation_bits=256, atol=1e-10, rtol=1e-8,
        cancellation_limit=1e8, moment_scale=1.0)
    max_order isa Integer && max_order >= 1 || error("max_order 必须为正整数")
    1 <= a <= b <= length(psi) || error("非法区间")
    precision_mode in (:float64, :checked, :bigfloat) || error("未知 precision_mode")
    precision_bits isa Integer && confirmation_bits isa Integer &&
        64 <= precision_bits < confirmation_bits || error("精度必须为整数且满足 64<=precision_bits<confirmation_bits")
    atol isa Real && rtol isa Real && isfinite(atol) && isfinite(rtol) &&
        atol > 0 && rtol >= 0 || error("需要有限的 atol>0, rtol>=0")
    cancellation_limit isa Real && isfinite(cancellation_limit) && cancellation_limit >= 1 ||
        error("cancellation_limit 必须为 >=1 的有限数")
    moment_scale isa Real && isfinite(moment_scale) && moment_scale > 0 ||
        error("moment_scale 必须为有限正数")
    # expect 在副本上正交化，原始 MPS 张量不变；这里只用其结果选择平移中心。
    centers = real.(expect(psi, "N"; sites=a:b))
    all(isfinite, centers) || error("平均粒子数非有限")
    M = Int(max_order)
    network = prepare_matrices(psi)
    # 超过六阶的 checked 模式直接使用两档高精度，避免双精度下溢伪装成一致。
    escalated = precision_mode==:bigfloat || (precision_mode==:checked && M>6)
    base = nothing
    if !escalated
        try
            base = evaluate(psi,a,b,M,centers; bits=53,scale=moment_scale,network=network)
        catch err
            precision_mode==:checked && err isa OverflowError || rethrow()
            escalated = true
        end
    end
    selected = base
    difference = fill(BigFloat(NaN), M)
    passed = falses(M)
    status = fill("not_checked", M)
    low_bits = 53
    if precision_mode == :checked && !escalated
        try
            reverse_result = evaluate(psi,a,b,M,centers; bits=53,direction=:right,scale=moment_scale,network=network)
            check = stability(base.cumulants, reverse_result.cumulants; atol=atol,rtol=rtol)
            difference, passed = check.delta, check.passed
            # 相消指标仅作触发器，不宣称是严格误差界；接近零时也可能触发。
            escalated = !all(passed) || any(base.cancellation .> cancellation_limit) ||
                any(reverse_result.cancellation .> cancellation_limit)
            status = [ok ? "stable_float64_directions" : "unstable_float64" for ok in passed]
        catch err
            err isa OverflowError || rethrow()
            escalated = true
        end
    end
    if escalated
        println("    高精度复核：$(precision_bits) / $(confirmation_bits) bits（原始 MPS 全环境重算）")
        low = evaluate(psi,a,b,M,centers; bits=precision_bits,scale=moment_scale,network=network)
        selected = evaluate(psi,a,b,M,centers; bits=confirmation_bits,direction=:right,scale=moment_scale,network=network)
        check = stability(low.cumulants, selected.cumulants; atol=atol,rtol=rtol)
        difference, passed = check.delta, check.passed
        low_bits = precision_bits
        status = [ok ? "stable_bigfloat_precisions" : "unstable_bigfloat" for ok in passed]
    end
    # 平均数只作为独立低阶诊断，不强迫 κ1=c，不强迫奇数阶为零。
    mean_reference = sum(BigFloat.(centers))
    mean_difference = abs(selected.cumulants[1] - mean_reference)
    # 参考 expect 使用 Float64，不能拿它要求 1e-30 的精度；高精度精度由两档收缩比较。
    mean_tolerance = max(BigFloat(atol)+BigFloat(rtol)*abs(mean_reference),
        BigFloat(64eps(Float64))*max(one(BigFloat),sum(abs,BigFloat.(centers))))
    mean_ok = mean_difference <= mean_tolerance
    variance_ok = M < 2 || selected.cumulants[2] >= -BigFloat(atol)
    if !mean_ok || !variance_ok
        passed .= false
        status .= "failed_low_order_check"
    end
    return merge(selected, (
        centers=centers, scale=moment_scale, order_status=status,
        order_passed=passed, stability_difference=difference, comparison_bits=low_bits,
        mean_reference=mean_reference, mean_difference=mean_difference,mean_check_tolerance=mean_tolerance,
        variance_nonnegative=variance_ok, precision_mode=String(precision_mode),
        float64_cumulants=isnothing(base) ? fill(BigFloat(NaN),M) : base.cumulants,
        float64_difference=isnothing(base) ? fill(BigFloat(NaN),M) : abs.(selected.cumulants .- base.cumulants),
    ))
end

function old_comparison(cfg, L, U, ratio, a, b, result)
    # 保留旧画图脚本的 Float64 求和与公式，比较同一区域的普通累积量。
    cfg.old_compare_max_bytes isa Integer && cfg.old_compare_max_bytes > 0 ||
        error("old_compare_max_bytes 必须为正整数字节数")
    n = min(4, length(result.cumulants))
    m = fill(NaN, 4)
    reasons = fill("not_requested", 4)
    messages = fill("", 4)
    for r in 1:n
        reasons[r] = "missing"
        file = joinpath(cfg.high_order_data_dir, "order$r", @sprintf(
            "order%d_ratio_%.3f_U_%.3f_dimer_%i_j_%.3f_L_%i_filling_%s_nsweep_%i.h5",
            r, ratio,U,cfg.dimerization,cfg.t1,L,cfg.filling,cfg.nsweep))
        isfile(file) || continue
        if filesize(file) > cfg.old_compare_max_bytes
            reasons[r] = "skipped_file_size"
            continue
        end
        try
            # 仅用 HDF5 读取普通数值数组的形状/类型元数据，不加载数组。
            # 压缩文件可能很小，不能仅靠磁盘文件大小限制内存。
            bytes = h5open(file,"r") do f
                dataset = f["C_$(r)_all"]
                dataset isa HDF5.Dataset || error("旧关联对象不是数值数组")
                size(dataset) == ntuple(_ -> L,r) || error("旧关联文件维度错误")
                T = eltype(dataset)
                isbitstype(T) && T <: Number || error("旧关联数组不是定长数值类型")
                prod(big.(size(dataset))) * sizeof(T)
            end
            if bytes > cfg.old_compare_max_bytes
                reasons[r] = "skipped_array_size"
                continue
            end
            m[r] = JLD2.jldopen(file,"r") do f
                arr = f["C_$(r)_all"]
                idx = ntuple(_ -> a:b, r)
                checked_real(sum(@view arr[idx...]))
            end
            reasons[r] = "loaded"
        catch err
            (err isa InterruptException || err isa OutOfMemoryError) && rethrow()
            reasons[r] = "read_error"
            messages[r] = sprint(showerror,err)
            @warn "旧文件对照失败；继续保存新方法结果" file reason=messages[r]
        end
    end
    m1,m2,m3,m4 = m
    old = [m1, m2-m1^2, m3-3m2*m1+2m1^3,
           m4-4m3*m1-3m2^2+12m2*m1^2-6m1^4]
    return Dict("raw_moments"=>m, "cumulants"=>old,
        "difference"=>Float64.(abs.(result.cumulants[1:n] .- old[1:n])),
        "file_status"=>reasons, "file_messages"=>messages,
        "note"=>"旧关联文件不能证明来自完全相同的 MPS；差异不是新方法误差界。")
end

function signature(cfg, source, a, b)
    st = stat(source)
    return (method=METHOD_VERSION, source=abspath(source), bytes=st.size, mtime=st.mtime,
        a=a,b=b,max_order=cfg.max_order,scale=cfg.moment_scale,mode=cfg.precision_mode,
        bits=cfg.precision_bits,confirmation_bits=cfg.confirmation_bits,
        atol=cfg.atol,rtol=cfg.rtol,cancellation_limit=cfg.cancellation_limit)
end

function save_result(path, result, metadata; overwrite=false)
    isfile(path) && !overwrite && error("结果已存在，禁止覆盖: $path")
    mkpath(dirname(path))
    temp, io = mktemp(dirname(path))
    close(io)
    try
        # 同时保存便于画图的 Float64 和保留高精度的十进制文本。
        # 不把只剩双精度的数字冒充高精度保存。
        setprecision(BigFloat, max(256, 2result.contraction_bits)) do
            JLD2.jldopen(temp, "w") do f
                for (key,value) in metadata
                    f[string(key)] = value
                end
                f["observable"] = OBSERVABLE
                f["method_version"] = METHOD_VERSION
                f["orders"] = collect(1:length(result.cumulants))
                f["cumulants"] = Float64.(result.cumulants)
                f["cumulants_decimal"] = string.(result.cumulants)
                f["taylor_coefficients_decimal"] = string.(result.coefficients)
                f["log_coefficients_decimal"] = string.(result.log_coefficients)
                f["centers"] = result.centers
                f["moment_scale"] = result.scale
                f["contraction_bits"] = result.contraction_bits
                f["comparison_bits"] = result.comparison_bits
                f["postprocess_bits"] = max(256, 2result.contraction_bits)
                f["precision_mode"] = result.precision_mode
                f["order_passed"] = result.order_passed
                f["order_status"] = result.order_status
                f["stability_difference_decimal"] = string.(result.stability_difference)
                f["float64_difference_decimal"] = string.(result.float64_difference)
                f["cancellation_indicator"] = Float64.(result.cancellation)
                f["norm2_decimal"] = string(result.norm2)
                f["mean_reference"] = Float64(result.mean_reference)
                f["mean_difference"] = Float64(result.mean_difference)
                f["mean_check_tolerance"] = Float64(result.mean_check_tolerance)
                f["variance_nonnegative"] = result.variance_nonnegative
                f["wavefunction_convergence_checked"] = false
                f["created_at"] = string(now())
            end
        end
        mv(temp, path; force=overwrite)
    finally
        isfile(temp) && rm(temp)
    end
    return path
end

function run_point(cfg, spec, U, ratio)
    L = spec.system_size
    region = hasproperty(spec,:region) ? spec.region : nothing
    a,b = subsystem_bounds(L; region=region,require_even_odd=cfg.require_even_odd,require_bulk=cfg.require_bulk)
    source = wavefunction_path(cfg,L,U,ratio)
    if !isfile(source)
        cfg.skip_missing_files ? (@warn "缺少波函数，跳过" source; return nothing) : error("缺少波函数: $source")
    end
    path = result_path(cfg,L,U,ratio,a,b)
    sig = signature(cfg,source,a,b)
    if isfile(path) && !cfg.overwrite_results
        old = load_result(path)
        if get(old,"signature",nothing)==sig && all(get(old,"order_passed",[false]))
            println("已存在且来源/设置一致、校验通过，跳过: $path")
            return path
        end
        error("已有结果与当前来源/设置不一致或未通过校验；请改输出目录或显式设置 overwrite_results=true: $path")
    end
    println("L=$L U=$U ratio=$ratio A=$a:$b 长度=$(b-a+1) 阶数=1:$(cfg.max_order)")
    psi = h5open(source,"r") do f
        read(f,"psi",MPS)
    end
    length(psi)==L || error("文件内链长与参数不一致: $source")
    started = time()
    result = calculate_cumulants(psi,a,b; max_order=cfg.max_order,precision_mode=cfg.precision_mode,
        precision_bits=cfg.precision_bits,confirmation_bits=cfg.confirmation_bits,
        atol=cfg.atol,rtol=cfg.rtol,cancellation_limit=cfg.cancellation_limit,moment_scale=cfg.moment_scale)
    total_mean = sum(real.(expect(psi,"N")))
    expected_N = cfg.filling=="1010" ? L/2 : cfg.filling=="1111" ? L : nothing
    !isnothing(expected_N) && !isapprox(total_mean,expected_N;atol=1e-7,rtol=1e-9) &&
        error("MPS 总粒子数 $total_mean 与 filling=$(cfg.filling) 不一致")
    metadata = Dict{String,Any}(
        "signature"=>sig, "source_wavefunction"=>abspath(source),
        "system_size"=>L,"U"=>Float64(U),"ratio"=>Float64(ratio),
        "dimerization"=>cfg.dimerization,"t1"=>cfg.t1,"filling"=>cfg.filling,"nsweep"=>cfg.nsweep,
        "left_site"=>a,"right_site"=>b,"region_length"=>b-a+1,
        "cut_bonds"=>[a-1,b],"region_mode"=>(isnothing(region) ? "original_even_to_odd" : "manual"),
        "total_number_mean"=>total_mean,"has_qns"=>hasqns(psi),
        "maxlinkdim"=>maxlinkdim(psi),"local_dimensions"=>dim.(siteinds(psi)),
        "atol"=>cfg.atol,"rtol"=>cfg.rtol,"elapsed_seconds"=>time()-started,
    )
    cfg.compare_old_orders && (metadata["old_method_comparison"] = old_comparison(cfg,L,U,ratio,a,b,result))
    save_result(path,result,metadata;overwrite=cfg.overwrite_results)
    for r in 1:cfg.max_order
        println("  κ$r = $(Float64(result.cumulants[r]))  [$(result.order_status[r])]")
    end
    !all(result.order_passed) && @warn "部分阶数尚未通过数值稳定性检查；文件已保存状态，画图默认排除这些阶数的数据点" path
    println("已保存: $path")
    return path
end

function main(cfg)
    cfg.blas_threads isa Integer && cfg.blas_threads >= 1 || error("blas_threads 必须为正整数")
    BLAS.set_num_threads(cfg.blas_threads)
    # 外层顺序运行，避免多个大 MPS 同时驻留；BLAS 线程数由用户显式控制。
    for ratio in cfg.ratio_values, spec in cfg.scan_specs, U in spec.u_values
        run_point(cfg,spec,U,ratio)
        GC.gc(false)
    end
    return nothing
end

end # module
