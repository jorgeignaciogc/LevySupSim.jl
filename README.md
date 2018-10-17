# LevySupSim
This code is an implementation of the algorithm presented in the paper

**Geometrically convergent simulation of the extrema of Lévy Processes**
by
- Jorge Ignacio González Cázares
- Aleksandar Mijatović
- Gerónimo Uribe Bravo

This is a simple algorithm that simulates approximately (geometrically quickly) 
the supremum of a Lévy process over a finite interval. 
The output is the triplet of the state of the process, 
the approximate supremum, and the time at which the supremum is attained.
Mathematically, take a Lévy process <a href="https://www.codecogs.com/eqnedit.php?latex=$X$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$X$" title="$X$" /></a> 
and a positive time horizon <a href="https://www.codecogs.com/eqnedit.php?latex=$T$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$T$" title="$T$" /></a>. 
Define 
<a href="https://www.codecogs.com/eqnedit.php?latex=$\overline{X}_T=\sup_{t\leq&space;T}X_t$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\overline{X}_T=\sup_{t\leq&space;T}X_t$" title="$\overline{X}_T=\sup_{t\leq T}X_t$" /></a>
and 
<a href="https://www.codecogs.com/eqnedit.php?latex=$\tau=\inf\{t>0:X_t\vee&space;X_{t-}=\overline{X}_T\}$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$\tau=\inf\{t>0:X_t\vee&space;X_{t-}=\overline{X}_T\}$" title="$\tau=\inf\{t>0:X_t\vee X_{t-}=\overline{X}_T\}$" /></a>
then the output takes 
<a href="https://www.codecogs.com/eqnedit.php?latex=$n$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$n$" title="$n$" /></a> 
steps and is an approximate simulation of 
<a href="https://www.codecogs.com/eqnedit.php?latex=$(X_T,\overline{X}_T,\tau)$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$(X_T,\overline{X}_T,\tau)$" title="$(X_T,\overline{X}_T,\tau)$" /></a>
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
