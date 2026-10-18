module ScalingDefinitions
using ITensors, LinearAlgebra
if !isdefined(ITensors,:MPS)
    @eval using ITensorMPS
end
const string_offset_d=0
const string_selector=:bulk_pair
const string_bulk_fraction=0.25
const nonlocal_nbar=0.5
function ITensors.op(::OpName"sigmaz_b", ::SiteType"Boson", j::Index)
    return op("Id", j) - 2 * op("N", j)
end

function ITensors.op(::OpName"sz_b", ::SiteType"Boson", j::Index)
    return 0.5 * op("Id", j) - op("N", j)
end

function ITensors.op(::OpName"exp_im_pi", ::SiteType"Boson", j::Index)
    return exp(im * pi * (0.5 * op("Id", j) - op("N", j)))
end

function ITensors.op(::OpName"density_minus_nbar", ::SiteType"Boson", j::Index)
    return op("N", j) - nonlocal_nbar * op("Id", j)
end

function ITensors.op(::OpName"exp_i_pi_density_minus_nbar", ::SiteType"Boson", j::Index)
    return exp(im * pi * (op("N", j) - nonlocal_nbar * op("Id", j)))
end

function string_mpo(st_i::Int, end_i::Int, sites; endpoint_op::String="sigmaz_b")
    os = OpSum()
    ops = tuple()
    ops = (ops..., endpoint_op, end_i)
    for i = end_i - 1:-1:st_i + 1
        ops = (ops..., "exp_im_pi", i)
    end
    ops = (ops..., endpoint_op, st_i)
    os += ops
    return MPO(os, sites)
end

function nonlocal_term_mpo(sites::Vector{<:Index}, i::Int)
    os = OpSum()
    ops = tuple()
    ops = (ops..., "density_minus_nbar", i)
    for k = i - 1:-1:1
        ops = (ops..., "exp_i_pi_density_minus_nbar", k)
    end
    os += ops
    return MPO(os, sites)
end
function direct_bulk_string_pair(L::Int)#!string从偶数开始
    if isodd(div(L, 4))
        st = Int(div(L, 4) + 1) - 2*string_offset_d
        en = Int(div(L, 4)*3)
    else
        st = Int(div(L, 4))+2 - 2*string_offset_d
        en = Int(div(L, 4)*3+1)
    end

    return st, en
end

function fraction_string_pair(L::Int)
    st = clamp(round(Int, string_bulk_fraction * L), 1, div(L, 2))
    en = L - st + 1
    return st, en
end

function selected_string_pair(L::Int)
    if string_selector == :edge_pair
        return 1, L
    elseif string_selector == :bulk_pair
        return direct_bulk_string_pair(L)
    elseif string_selector == :fraction_pair
        return fraction_string_pair(L)
    elseif string_selector == :mean_bulk_pairs || string_selector == :max_abs
        return direct_bulk_string_pair(L)
    else
        error("Unknown string_selector = $string_selector")
    end
end
const entanglement_entropy_definition = "block_entropy_bulk_window"
const half_chain_entropy_definition = "single_center_bond_cut"
const entanglement_spectrum_definition = "single_bulk_bond_cut"

# CenterEntropy is the bulk block entropy on the middle [left:right] interval.
# HalfChainEntropy is the bipartition entropy across the center bond.
function entropy_from_singular_values(S; cutoff::Float64 = 1e-12)
    probs = Float64[]
    for n in 1:dim(S, 1)
        p = S[n, n]^2
        if p > cutoff
            push!(probs, p)
        end
    end
    if isempty(probs)
        return NaN
    end

    probs = probs ./ sum(probs)
    entropy = 0.0
    for p in probs
        entropy -= p * log(p)
    end
    return entropy
end

function central_bulk_interval(L::Int)
    q = div(L, 4)
    if isodd(q)
        left = q + 1
        right = 3 * q
    else
        left = q + 2
        right = 3 * q + 1
    end

    left = max(left, 2)
    right = min(right, L - 1)
    return left < right ? (left, right) : nothing
end

# This is the expensive middle-block density-matrix path. It is called only
# when observable_name == :center_entropy.
function center_entropy_from_mps(psi::MPS, L::Int)
    interval = central_bulk_interval(L)
    interval === nothing && return NaN
    left, right = interval

    orthogonalize!(psi, right)

    block_left_link = linkind(psi, left - 1)
    block_right_link = linkind(psi, right)
    if isnothing(block_left_link) || isnothing(block_right_link)
        return NaN
    end

    dagpsi = prime(dag(psi), linkinds(psi))
    rho = psi[left] * dagpsi[left]
    for j in (left + 1):right
        rho = rho * psi[j]
        rho = rho * dagpsi[j]
    end

    D, _ = eigen(rho, (block_left_link, block_right_link), (block_left_link', block_right_link'))
    probs = [abs(x) for x in diag(D) if abs(x) > 1e-12]
    isempty(probs) && return NaN
    probs = probs ./ sum(probs)

    center_entropy = 0.0
    for p in probs
        center_entropy -= p * log(p)
    end
    return center_entropy
end

# Both half-chain entropy and the entanglement gap start from the same center-
# bond SVD, but the caller requests only one observable at a time.
function center_bond_singular_values_from_mps(psi::MPS, L::Int)
    interval = central_bulk_interval(L)
    interval === nothing && return nothing
    left, right = interval
    spectrum_cut = clamp(div(L, 2), left, right - 1)
    orthogonalize!(psi, spectrum_cut)

    spectrum_left_link = linkind(psi, spectrum_cut - 1)
    isnothing(spectrum_left_link) && return nothing

    _, S, _ = svd(psi[spectrum_cut], (spectrum_left_link, siteind(psi, spectrum_cut)))
    return S
end

function entanglement_gap_from_singular_values(S; cutoff::Float64 = 1e-12)
    entanglement_spectrum = Float64[]
    for n in 1:dim(S, 1)
        lambda = S[n, n]
        if lambda > cutoff
            push!(entanglement_spectrum, -2.0 * log(lambda))
        end
    end

    sort!(entanglement_spectrum)
    es1 = length(entanglement_spectrum) >= 1 ? entanglement_spectrum[1] : NaN
    es2 = length(entanglement_spectrum) >= 2 ? entanglement_spectrum[2] : NaN
    entanglement_gap = isfinite(es1) && isfinite(es2) ? es2 - es1 : NaN

    return entanglement_gap, es1, es2
end
end
