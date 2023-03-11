# LevySupSim.jl
This code is an implementation of the algorithm presented in the paper

A Julia package for the approximate geometrically convergent simulation of
the supremum of a Lévy process over a finite interval. 
The output is the triplet of the state of the process, 
the approximate supremum, and the time at which the supremum is attained.
Mathematically, take a Lévy process $X$ and a positive time horizon $T$. 
Define 

$$\overline{X}_{T} = \sup_{t\leq T}X_{t}
\quad\text{and}\quad
\tau_{T}=\inf{\\{}t>0:X_{t}\vee X_{t-}=\overline{X}_{T}{\\}},$$

then the output takes $n$ 
steps and is an approximate simulation of 

$$(X_{T},\overline{X}_{T},\tau_{T})$$

such that:
1. the sample is exact in the first coordinate
2. converges in $L^1$
at rate 
$\eta^n$
for some $\eta\in[1/2,2/3]$
(depending on the Lévy characteristics) in the second coordinate
3. is off by at most 
$T2^{-n}$ (with probability 1) in the third coordinate.

## What is included in this code?
The main function is

    rand_levy_sup

which requires a sampler for the increments of the Levy process 
and returns the triplet in the form of a tuple of vectors (see the example below).
### Lévy increment samplers included
Fix some $t\in(0,\infty)^n$,
with $n\in\mathbb{N}$.

Variance Gamma

    VG(t,Θ)
where 

$$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_{+}\times\mathbb{R}\times\mathbb{R}_{+}\times\mathbb{R}.$$
    
Normal Inverse Gaussian

    NIG(t,Θ)
where 

$$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_{+}\times\mathbb{R}\times\mathbb{R}_{+}\times\mathbb{R}.$$

Weakly Stable

    WS(t,Θ)
where 

$$\Theta=(\alpha,\beta,\gamma,b)\in(0,2]\times[-1,1]\times\mathbb{R}_{+}\times\mathbb{R}.$$

See the code and the example below for details.

### Notes: StableSupremum
This distributions' implementation relies on a recent paper by the authors of the package. See the article for details at:  
Jorge I. González Cázares and Aleksandar Mijatović and Gerónimo Uribe Bravo, *Geometrically convergent simulation of the extrema of Lévy Processes*, Mathematics of Operations Research 47(2):1141-1168. doi:10.1287/moor.2021.1163 (2021) and arXiv:1810.11039v3 (2021).

## Example: Barrier Option + NIG
Consider the payoff 
$P(x,y,t)=e^{-rT} (x-K)^+ 1_{y < B}$

1. time horizon $T = 1$
2. initial stock price $S_0 = 100$
3. strike price $K = 100$
4. barrier level $B = 115$
5. risk-free interest rate $r = 0.05$

and an NIG process with parameters $\sigma = 0.1836, \theta = -0.1313, \kappa = 1.2819$.
Assume the stock price process is $S_t=S_0 e^{X_t}$.

    using Statistics, StatsBase
    # Parameters:
    T, r, S0, K, B = 1., 0.05, 100., 100., 115.
    σ, θ, κ = 0.1836, -0.1313, 1.2819
    b = r + (sqrt(1 - σ^2*κ - 2*θ*κ) - 1)/κ # to get a martingale stock price process
    Θ = (σ, θ, κ, b)

    # Payoff:
    g(X,T,r,S0,K,B) = exp(-r*T) .* max.(S0.*exp.(X[1]).-K,0) .* (exp.(X[2]).<B/S0)
    g(X) = g(X,T,r,S0,K,B)

    # Sampler:
    f(ℓ) = NIG(ℓ,Θ)

    # Number of steps and samples
    n, N = [20,50,100,200,400], 10^6
    
    for i = n
      t = convert(Float64,time_ns())
      X = rand_levy_sup(f,T,i,N)
      price = mean_and_std(g(X))
      t = convert(Float64,time_ns()) - t
      t = t*1e-9
      print("
        For n = $i,
        Estimate: 𝔼[g(X)]=$(price[1]),
        Correction term with confidence level 95%: $(1.96*price[2]/N),
        Time taken = $t seconds"
      )
    end
        

## Author and Contributor List
Jorge I. González Cázares

Aleksandar Mijatović  
Gerónimo Uribe Bravo
