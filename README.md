# LevySupSim.jl
This code is an implementation of the algorithm presented in the paper

A Julia package for exact simulation of the supremum of a stable process.
This is a simple algorithm that simulates approximately (geometrically quickly) 
the supremum of a Lévy process over a finite interval. 
The output is the triplet of the state of the process, 
the approximate supremum, and the time at which the supremum is attained.
Mathematically, take a Lévy process <a href="https://www.codecogs.com/eqnedit.php?latex=$X$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$X$" title="$X$" /></a> 
and a positive time horizon <a href="https://www.codecogs.com/eqnedit.php?latex=$T$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$T$" title="$T$" /></a>. 
Define 
<a href="https://www.codecogs.com/eqnedit.php?latex=$\overline{X}_T=\sup_{t\leq&space;T}X_t$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\overline{X}_T=\sup_{t\leq&space;T}X_t$" title="$\overline{X}_T=\sup_{t\leq T}X_t$" /></a>
and 
<a href="https://www.codecogs.com/eqnedit.php?latex=$\tau_T=\inf\{t>0:X_t\vee&space;X_{t-}=\overline{X}_T\}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\tau_T=\inf\{t>0:X_t\vee&space;X_{t-}=\overline{X}_T\}$" title="$\tau_T=\inf\{t>0:X_t\vee X_{t-}=\overline{X}_T\}$" /></a>
then the output takes 
<a href="https://www.codecogs.com/eqnedit.php?latex=$n$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$n$" title="$n$" /></a> 
steps and is an approximate simulation of 
<a href="https://www.codecogs.com/eqnedit.php?latex=$(X_T,\overline{X}_T,\tau_T)$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$(X_T,\overline{X}_T,\tau_T)$" title="$(X_T,\overline{X}_T,\tau_T)$" /></a>
such that:
1. is exact in the first coordinate
2. converges in <a href="https://www.codecogs.com/eqnedit.php?latex=$L^1$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$L^1$" title="$L^1$" /></a>
at rate 
<a href="https://www.codecogs.com/eqnedit.php?latex=$\eta^n$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\eta^n$" title="$\eta^n$" /></a>
for some <a href="https://www.codecogs.com/eqnedit.php?latex=$\eta\in[1/2,2/3]$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\eta\in[1/2,2/3]$" title="$\eta\in[1/2,2/3]$" /></a>
(depending on the Lévy characteristics) in the second coordinate
3. is off by at most 
<a href="https://www.codecogs.com/eqnedit.php?latex=$T2^{-n}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$T2^{-n}$" title="$T2^{-n}$" /></a>
(with probability 1) in the third coordinate.

## What is included in this code?
The main function is

    rand_levy_sup

which requires a sampler for the increments of the Levy process 
and returns the triplet in the form of a tuple of vectors (see the example below).
### Lévy increment samplers included
For some <a href="https://www.codecogs.com/eqnedit.php?latex=$t\in(0,\infty)^n$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$t\in(0,\infty)^n$" title="$t\in(0,\infty)^n$" /></a>,
with <a href="https://www.codecogs.com/eqnedit.php?latex=$n\in\mathbb{N}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$n\in\mathbb{N}$" title="$n\in\mathbb{N}$" /></a>.
Variance Gamma

    VG(ℓ,Θ)
where <a href="https://www.codecogs.com/eqnedit.php?latex=$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_&plus;\times\mathbb{R}\times\mathbb{R}_&plus;\times\mathbb{R}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_&plus;\times\mathbb{R}\times\mathbb{R}_&plus;\times\mathbb{R}$" title="$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_+\times\mathbb{R}\times\mathbb{R}_+\times\mathbb{R}$" /></a>.
    
Normal Inverse Gaussian

    NIG(ℓ,Θ)
where <a href="https://www.codecogs.com/eqnedit.php?latex=$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_&plus;\times\mathbb{R}\times\mathbb{R}_&plus;\times\mathbb{R}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_&plus;\times\mathbb{R}\times\mathbb{R}_&plus;\times\mathbb{R}$" title="$\Theta=(\sigma,\theta,\kappa,b)\in\mathbb{R}_+\times\mathbb{R}\times\mathbb{R}_+\times\mathbb{R}$" /></a>.

Weakly Stable

    WS(ℓ,Θ)
where <a href="https://www.codecogs.com/eqnedit.php?latex=$\Theta=(\alpha,\beta,\gamma,b)\in(0,2]\times[-1,1]\times\mathbb{R}_&plus;\times\mathbb{R}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\Theta=(\alpha,\beta,\gamma,b)\in(0,2]\times[-1,1]\times\mathbb{R}_&plus;\times\mathbb{R}$" title="$\Theta=(\alpha,\beta,\gamma,b)\in(0,2]\times[-1,1]\times\mathbb{R}_+\times\mathbb{R}$" /></a>.

See the code and the example below for details.

### Notes: StableSupremum
This distributions' implementation relies on a recent paper by the authors of the package. See the article for details at:  
Jorge I. González Cázares and Aleksandar Mijatović and Gerónimo Uribe Bravo, *Geometrically convergent simulation of the extrema of Lévy Processes* (2018).

## Example: Barrier Option + NIG
Consider the payoff 
<a href="https://www.codecogs.com/eqnedit.php?latex=$P(x,y,t)=e^{-rT}&space;(x-K)^&plus;&space;1_{y&space;<&space;B}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$P(x,y,t)=e^{-rT}&space;(x-K)^&plus;&space;1_{y&space;<&space;B}$" title="$P(x,y,t)=e^{-rT} (x-K)^+ 1_{y < B}$" /></a>

1. time horizon <a href="https://www.codecogs.com/eqnedit.php?latex=$T&space;=&space;1$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$T&space;=&space;1$" title="$T = 1$" /></a>
2. initial stock price <a href="https://www.codecogs.com/eqnedit.php?latex=$S_0&space;=&space;100$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$S_0&space;=&space;100$" title="$S_0 = 100$" /></a>
3. strike price <a href="https://www.codecogs.com/eqnedit.php?latex=$K&space;=&space;100$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$K&space;=&space;100$" title="$K = 100$" /></a>
4. barrier level <a href="https://www.codecogs.com/eqnedit.php?latex=$B&space;=&space;115$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$B&space;=&space;115$" title="$B = 115$" /></a>
5. risk-free interest rate <a href="https://www.codecogs.com/eqnedit.php?latex=r&space;=&space;0.05" target="_blank"><img src="https://latex.codecogs.com/gif.latex?r&space;=&space;0.05" title="r = 0.05" /></a>

and an NIG process with parameters <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma&space;=&space;0.1836,&space;\theta&space;=&space;-0.1313,&space;\kappa&space;=&space;1.2819" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma&space;=&space;0.1836,&space;\theta&space;=&space;-0.1313,&space;\kappa&space;=&space;1.2819" title="\sigma = 0.1836, \theta = -0.1313, \kappa = 1.2819" /></a>.
Assume the stock price process is <a href="https://www.codecogs.com/eqnedit.php?latex=S_t=S_0&space;e^{X_t}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S_t=S_0&space;e^{X_t}" title="S_t=S_0 e^{X_t}" /></a>.

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
