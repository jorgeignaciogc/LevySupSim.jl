#= This code is a simple implementation of the paper
"Geometrically convergent simulation of the extrema of Lévy Processes"
by González Cázares, Mijatović, Uribe Bravo
=#

#= The main algorithm.

Input:
1. f: A sampler from a parametric distribution (parametrised by time t ∈ ℝ_+)
such that f(t) is distributed as X_t, where X is the Lévy process
(works better if vecotrised).
3. T: Time horizon T > 0.
3. n: Number of steps n ∈ ℕ.
4. N: (optional) Number of samples N ∈ ℕ.

Output:
A tuple (Y1,Y2,Y3) such that Y1, Y2, Y3 ∈ ℝ^N and such that the vectors
{(Y1[i], Y2[i], Y3[i])}_{i=1}^N are iid and, for every i,
the vector (Y1[i], Y2[i], Y3[i]) is approximately distributed as
(X_T, sup_{t<=T}X_t, inf{t>0: X_t ∨ X_{t-} = sup_{s<=T}X_s}).

Please note:
It is not recommended to take N > 10^6 as the RAM consumption might
kill julia's kernel. Please instead split N = N_1 + ... + N_k
and then run the function for all such N_i.
=#
function rand_levy_sup(f,T,n,N)
    Y1 = zeros(N)
    Y2 = zeros(N)
    Y3 = zeros(N)

    ℓ1 = zeros(N) .+ T
    ℓ = zeros(N)

    try
        for i = 1:(n-1)
            ℓ = ℓ1 .* rand(N)

            ξ = f(ℓ)
            Y1 .+= ξ
            Y2 .+= max.(ξ,0)
            Y3 .+= ℓ .* (ξ .> 0)

            ℓ1 = ℓ1 .- ℓ
        end
        ξ = f(ℓ1)
        Y1 .+= ξ
        Y2 .+= max.(ξ,0)
        Y3 .+= ℓ1 .* (ξ .> 0)
    catch
        for i = 1:(n-1)
            ℓ = ℓ1 .* rand(N)

            ξ = f.(ℓ)
            Y1 .+= ξ
            Y2 .+= max.(ξ,0)
            Y3 .+= ℓ .* (ξ .> 0)

            ℓ1 = ℓ1 .- ℓ
        end
        ξ = f.(ℓ1)
        Y1 .+= ξ
        Y2 .+= max.(ξ,0)
        Y3 .+= ℓ1 .* (ξ .> 0)
    end

    return (Y1,Y2,Y3)
end

#= PROCESSES
Examples of samplers f for the Levy processes:
Variance Gamma, Normal Inverse Gaussian and (Weakly) Stable
=#

using Distributions

#= Variance Gamma (Alg. 6.11 in Tankov (2004))

Input:
1. ℓ: a vector ℓ ∈ ℝ_+^n
2. Θ: a vector of parameters Θ = (σ,θ,κ,b) ∈ ℝ_+ × ℝ × ℝ_+ × ℝ

Output:
Vector X such that for each i (using the notation of Tankov (2004)),
X[i] is distributed as VarianceGamma(σ,θ,κ)_ℓ[i] + b*ℓ[i]
or equivalently,
𝔼[exp(i*u*X[i])] = exp(ℓ[i]*(b+1/κ) - ℓ[i]*√(1 - 2*i*u*θ*κ + κ*σ^2*u^2)/κ)

Remarks:
1. For exp(X_t-rt) to be a martingale it is required that
b = r + log(1 + θ*κ - σ^2*κ/2)-1)/κ
2. Calibration in Schoutens:
σ = 0.1213, θ = −0.1436, κ = 0.1686
=#
function VG(ℓ,Θ)
    (σ,θ,κ,b) = Θ
    Y = rand(Normal(0,σ),N)
    G = [rand(Gamma(ℓ[j]/κ,κ)) for j=1:N]
    return Y .* sqrt.(G) .+ θ .* G .+ b .* ℓ
end

#= Normal Inverse Gaussian (Algs 6.12 and 6.9 in Tankov (2004))

Input:
1. ℓ: a vector ℓ ∈ ℝ_+^n
2. Θ: a vector of parameters Θ = (σ,θ,κ,b) ∈ ℝ_+ × ℝ × ℝ_+ × ℝ

Output:
Vector X such that for each i (using the notation of Tankov (2004)),
X[i] is distributed as VarianceInverseGaussian(σ,θ,κ)_ℓ[i] + b*ℓ[i]
or equivalently,
𝔼[exp(i*u*X[i])] = exp(ℓ[i]*(b+1/κ) - ℓ[i]*√(1 - 2*i*u*θ*κ + κ*σ^2*u^2)/κ)

Remarks:
1. For exp(X_t-rt) to be a martingale (see Ex. 15.3 in Tankov (2004))
it is required that b = r + (√(1 - σ^2*κ - 2*θ*κ) - 1)/κ
2. Calibration in Schoutens:
σ = 0.1836, θ = −0.1313, κ = 1.2819
=#
function NIG(ℓ,Θ)
    (σ,θ,κ,b) = Θ
    G = rand(Normal(),N) .^ 2
    G = ℓ .+ (G .- sqrt.((4/κ).*ℓ.*G .+ G .^ 2)) .* (κ/2)
    U = rand(N) .< (ℓ ./ (G+ℓ))
    G = abs.(G .* U .+ (ℓ .^ 2 ./ G) .* (1 .- U))
    return rand(Normal(),N) .* σ .* sqrt.(G) .+ θ .* G .+ b .* ℓ
end

#= Weakly Stable (Chambers-Mellows-Stuck and Alg. 6.9 in Tankov (2004))

Input:
1. ℓ: a vector ℓ ∈ ℝ_+^n
2. Θ: a vector of parameters Θ = (α,β,γ,b) ∈ ℝ_+ × ℝ × ℝ_+ × ℝ

Output:
Vector X such that for each i (using Zolotarev's (C) form),
𝔼[exp(i*u*X[i])] = exp(ℓ[i]*(i*u*b* - γ*|u|^α*exp(-i*π*t*α*sgn(u)/2)),
where t = β*1_{α<=1}+ β*((α-2)/α)*1_{α>1}

Remarks:
1. Only has exponential moments if α=2 or if 1<α<2 and β=-1.
2. There is an additional weakly stable case for α=1 not covered here.
=#
function WS(ℓ,Θ)
    (α,β,γ,b) = Θ
    G = rand(Uniform(-pi/2,pi/2),N)
    G = rand(Exponential(),N) .^ (1-1/α) .* sin.(α .* (G + θ)) ./
        (cos.(G) .^ (1/α) .* cos.(G - α .* (G + θ)) .^ (1-1/α))
    return G .* (γ .* ℓ .^ (1/α)) .+ b .* ℓ
end
