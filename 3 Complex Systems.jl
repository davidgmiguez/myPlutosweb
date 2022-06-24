### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d096a6be-65a1-428d-9bfb-da7fe89f4c19
using Plots 

# ╔═╡ f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
using DifferentialEquations

# ╔═╡ 6e728cf8-72e7-4567-90e9-a97f5f6a1f28
using ParameterizedFunctions

# ╔═╡ 481d6298-106a-4f6a-b822-7fda0a24bbe0
using PlutoUI

# ╔═╡ 6d9759c3-63d9-4166-8884-cfaa99ee33c4
using HypertextLiteral

# ╔═╡ 081daad5-b960-4c64-bdec-c0ecb0d6896b
md" # Simple models of population dynamics

## Unconstrained growth 

The most basic approach to population growth is to begin with the assumption that every individual produces two offspring in its lifetime, then dies, which would double the population size each generation. The same approach is for cells that divide and produce two identical cells, such as bacteria.

Stem cells also undergo a similar type of scheme: during the early stages of development, they undergo a phase of proliferation that expands the population, in a process that called proliferative divisions. Later on, they start to differentiate to produce the different types of progeny. A very important question in develpmental biology is to understand how these populations grow in time and how this growth may be affected by external stimuly. 

A first simplified population approach is to assume a constant relative growth rate. In a discrete model, this rate would represent growth over time intervals of some fixed length. Consider as an example the growth of a population of cells with an initial number of 'p(t=0)=p_0'. 

After a given time interval '∆t', a percentage of the population will reproduce, resulting in a number of new cells '∆p'. This number of newly produced cells has to be  proportional to the initial number of cells 'p(0)' (if a population of 20,000 cells produces 1200 new cells in 1 h, then a 4-fold bigger population of 80,000  will produce 4 times as many, i.e., 4800, new cells in 1 h). We can write this proportionality as a discrete system in the following way

```math
p_{t+1} = r \cdot p_t
```

where $p_t$ is the population at the start of hour $n$, $p_{t+1}$ is the number of cells at the end of the time interval, and $r$ is a fixed growth factor. Here, each term $p_t$ is simply multiplied by $r$ to produce the next term. Because the newly produced cells always add to the population, i.e. the number of cells in the system increases, it is straightforward to see that the number of cells in the system at succesive time imntervals will be

```math
\begin{align} 
p_0\\
p_1 &=r \cdot p_0 \\
p_2 &=r \cdot p_1 = r^2 \cdot p_0 \\
p_3 &=r \cdot p_2 = r^3 \cdot p_0 \\
···\\
\end{align}
```


So, in general 

```math
p_{t}=r^t \cdot p_0
```

The factor $r$ exceeds unity by the relative growth rate. For example, if the population increases by $6\%$ each hour then $r = 1.06$. In general, with a positive relative growth rate, the solution to is an exponential function with base $r > 1$
"

# ╔═╡ 4e2cead3-1c4f-48b4-8f5a-f78c4efaee2c
begin
	r_slide = @bind r html"<input type=range min=0.1 max=1.4 step=0.1>"
	
	md"""
	**Set the growth rate?**
	
	value of r: $(r_slide)
	
	"""
end

# ╔═╡ 16c39710-8f90-45c8-983a-25438019d90c
time=collect(1:1000) .* 0.01

# ╔═╡ bd8ee393-6193-4662-b199-edbe339ffc31
plot(time,r .^(time),ylims = (0,10))

# ╔═╡ 951f16b4-f2d4-4855-91b2-b2ab2cf7df43
md"we can do the same iteratively, calculating teh next value of p based on the previous value, using a for cycle"

# ╔═╡ 99d88df8-e6c8-4fe7-a9ab-61f908912e92
function unconstrained_growth_1(p, r; dt=0.01)
traj = []
for t in 1:10 # arbitrary, just leave enough for it to reach a steady state given the dt
	#p += dt * (p * (r - 1))	
	#n = n + (n * (r - 1))	
	p = p * r; #p = p * r; 
push!(traj,p)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 04bb79a4-5a72-4439-a3cd-b2b5ff985300
scatter([1:10],unconstrained_growth_1(1, r; dt=0.01),ylims = (0,10),xlabel=("Time"),ylabel=("Number of cells"),title=("Unconstrained growth, r=$r"))

# ╔═╡ 06cbef49-b654-41dd-ac10-13474ef316a7
md"we can do the same, to get more time points, using a dt
```math
\begin{align} 
p_{t+1}= r \cdot p_{t} = r \cdot p_{t} + p_{t} - p_{t} \\
p_{t+1}= p_{t} + (r \cdot p_{t} - p_{t})\\
p_{t+1}= p_{t} + p_{t} (r  - 1)
\end{align}
```
"

# ╔═╡ 9a68d323-e13d-42bb-970a-0ab5630155d5
function unconstrained_growth_2(p, r; dt=0.01)
traj = []
for t in 1:1000 # arbitrary, just leave enough for it to reach a steady state given the dt
	p += dt * (p * (r - 1))	
	#n = n + (n * (r - 1))	
	#p = p * r; #p = p * r; 
push!(traj,p)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 13aae8fd-f633-4073-a5ab-1e89f99e8568
plot(time,unconstrained_growth_2(1, r; dt=0.01),ylims = (0,10),xlabel=("Time"),ylabel=("Number of cells"),title=("Unconstrained growth, r=$r"))

# ╔═╡ 7b74322b-8481-4311-a77c-1a62cfb5b15c
md"
Although this population model might be fairly accurate in the short term, this type of exponential growth is not realistic in the long term. For instance, there will be potentially a limitation in environmental factors (size, nutrients...). A simple approximation but much more realistic, it is to assume that $r$ changes with the population size. 

```math
p_{n+1} = r(p_n) \cdot p_n
```

where the function $r(p)$ is set to decrease as the population p increases. From a mathematical standpoint, it is natural to begin an investigation of such models by considering the case of a linear function r(p).

Lets think of a population that initially grows with very little environmental constraint ($r>1$), a situation that might arise when a few members of a new species are introduced into an environment rich in nutrients and habitable area. Over time the population will increase until it approaches some maximum sustainable size, eventually reaching an equilibrium ($r=1$).

Introducing specific numbers, suppose that a population of p = 1 experiences a growth factor of r(1) = 1.1. Assume also that the equilibrium population size is 50, so that the corresponding growth factor is r(50) = 1.



"

# ╔═╡ a1ce2262-da8e-496a-81d8-10ff5246c17f
begin
	plot([0,100],[1.1,0.9],line = (:line, 4))
	title!("r for constrained growth")
	xlabel!("Number of cells")
	ylabel!("r")
end

# ╔═╡ 85defe33-1c87-4e81-a7d2-a363bee3e699
md"
The linear function r is thus determined to be

```math
Slope=\frac{1.0 - 1.1}{50-1}=\frac{- 0.1}{49}=-2 \cdot 10^{-3}\\
p(0)=1.1
```
Therefore, the function takes the form

```math
r(p) = 1.1 - 2 \cdot 10^{-3} p
```

and the difference equation is now
```math
p_{n+1} = p_{n} (1.1 -2 \cdot 10^{-3} p_{n})  = 
```

```math
p_{n+1} =  p_{n} -  p_{n} + p_{n} (1.1 -2 \cdot 10^{-3} p_{n})  =
```

```math
p_{n+1} =  p_{n} +  p_{n} ( p_{n} (1.1 -2 \cdot 10^{-3} p_{n}) -1) 
```

This is an instance of discrete logistic growth. 
"

# ╔═╡ b043d65b-a214-482a-96cb-e8c075814490
growth_factor(p) = 1.10 - 0.002 * p

# ╔═╡ 842e2c41-72a3-443c-82db-ab1ff3612a12
begin
	plot([0,100],[growth_factor(0),growth_factor(100)],line = (:line, 4))
	title!("r for constrained growth")
	xlabel!("Number of cells")
	ylabel!("r")
end

# ╔═╡ 79fdee9d-9e12-4d2c-ad61-bd2b4c6dbb1a
function constrained_growth_3(n, r,final_time)
traj = []
for t in 1:final_time # arbitrary, just leave enough for it to reach a steady state given the dt
n = n * (1.10 - 0.002 * n)
push!(traj,n)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 86804a6f-9178-40f3-ae82-cc6723412ae8
function unconstrained_growth_3(p, r; dt=0.01)
traj = []
for t in 1:100 # arbitrary, just leave enough for it to reach a steady state given the dt
	#p += dt * (p * ((1.10 - 0.002 * p) - 1))	
	p = p * (1.10 - 0.002 * p)
	#p = p * r; #p = p * r; 
push!(traj,p)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 8510e3df-ea93-4c25-ac3d-1069b067a62d
plot([1:100],unconstrained_growth_3(1, r; dt=0.01),ylims = (0,100),xlabel=("Time"),ylabel=("Number of cells"),title=("Unconstrained growth, r=$r"))

# ╔═╡ 753586ff-1721-48e3-befd-bce8f72ae285


# ╔═╡ 62542751-f6ae-4bf6-9a20-c362011f68d5
function Constrained_growth_4(n, rr; dt=0.01)
traj = []
for t in 1:1000 # arbitrary, just leave enough for it to reach a steady state given the dt
	n += dt * (n * (rr - 0.028 * n)) 
#n = n * (1.1 - 0.0004 * n) 
push!(traj,n)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 831a828f-166e-4467-a94e-3cc7981a3e16
begin
	rr_slide = @bind rr html"<input type=range min=0.1 max=1.4 step=0.1>"
	
	md"""
	**Set the growth rate?**
	
	value of r: $(rr_slide)
	
	"""
end

# ╔═╡ 1c127c36-c5fa-4a1f-beb1-8f3774abffaf
begin
	plot(Constrained_growth_4(1, rr; dt=0.015),ylims = (0,50))
			title!("Constrained growth")
			xlabel!("Time")
			ylabel!("Number of cells")
end

# ╔═╡ 3abe473c-67c7-4897-9345-1c39b5a2077c
function logistic_growth(n, rr )
traj = []
for t in 1:50 # arbitrary, just leave enough for it to reach a steady state given the dt
n += n * rr * (1. - n)
push!(traj,n)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 1dc48319-02d1-44d8-9873-2ebd28b358a8
function logistic_growth2(n, rr; dt=0.01)
traj = []
for t in 1:500 # arbitrary, just leave enough for it to reach a steady state given the dt
n += dt * n * rr * (1. - n) 
push!(traj,n)
end
return traj#[end-20:end] # this is sampling from the steady state
end

# ╔═╡ 8e2e2fec-f264-45df-95fa-556ac631a8f2
plot(logistic_growth(0.01, rr))

# ╔═╡ 7ba042d6-107d-4188-a877-c7f7d7b35c46
md"
## Unconstrained growth Continuos

Stem Cell Proliferation Growth In general stem cells undergo a phase of proliferation that expands the population, then they start to differentiate to produce the different types of progeny. Let’s calculate how a population of stem cells grow in time. Lets assume we have an initial number of cells `N` in a population 

After a given time we have a number of new stem cells `∆N` (number of newly produced cells) during a given time interval, `∆t`, is proportional to the initial  number of cells `N`. (If a population of 20,000 cells produces 1200 new cells in 1 h, a 4-fold bigger population of 80,000 cells of the same type of microorganism will produce 4 times as many, viz. 4800, new cells in 1 h): 

```math
\frac{\Delta N}{\Delta t} \sim N\tag{1}
```

Because the newly produced cells always add to the population, i.e. N increases steadily, we have to regard time intervals (and accordingly numbers of newly formed cells) that are as small as possible. Mathematically we thus have to deal with infinitesimal increases or differentials: 

```math
\frac{\mathrm{d} N}{\mathrm{d} t} \sim N\tag{2}
```

To get from proportionality  to an equation, a proportionality factor, μ, is introduced. This is the specific growth rate or often termed simply growth rate. Sorting of variables yields: 

```math
\frac{\mathrm{d} N}{N}  = \mu \mathrm{d} t\tag{3}
```

If we integrate this differential equation:

```math
\int \frac{\mathrm{d} N}{N}  = \int \mu \mathrm{d} t \tag{4}
```

 we obtain:
 
 ```math
\ln N = \mu t + C \tag{5}
```
 
 The undefined integration constant `C` can be fixed if we define the initial conditions: $N(t = 0) = N_0$

 ```math
\ln N_0 = \mu 0 + C =C \tag{6}
```
 
 so the equation becomes:
 
  ```math
\ln N = \mu t + \ln N_0 \tag{7}
```
  
  and after applying the properties of logatithms:
  
```math
\begin{align*}
\ln N - \ln N_0 &= \mu t  \tag{8}\\
\ln \frac{N}{N_0}&= \mu t  \tag{9}\\
\frac{N}{N_0} &= e^{\mu t}  \tag{10}\\
N &= N_0 e^{\mu t}  \tag{11}\\
\end{align*}
```
  
   If the time that has passed is exactly the average length of the cell cycle of the stem cell population, we can write that in $t = T$ the population has doubled, so $N_0$ has increased to $2 N_0$:

  
```math
\begin{align*}
2 N &= N_0 e^{\mu T}  \\
\frac{2 N_0}{N_0 }& =2= e^{\mu T}  \tag{13}\\
 \ln 2&= \mu T  \tag{14}\\
\end{align*}
```
   
   So, the proportionality constant is related to the cell cycle length as:
```math
\begin{align*}
 \mu = \frac{\ln 2}{T}  \tag{15}
\end{align*}
```


"

# ╔═╡ 85bc871e-6067-4a1d-bbb9-e3a9b0dfbd85
md"
# Logistic Growth 

We can introduce in the original equation a term that accounts for the limitation of resources 

```math
\begin{align*}
1 -\frac{N}{K} \tag{16}
\end{align*}
```

This factor is close to 1 (i.e., has no effect) when N is much smaller than K, and which is close to 0 when N is close to K. The resulting differential equation is called is the logistic growth model.


```math
\begin{align*}
\frac{\mathrm{d} N}{\mathrm{d} t}=\mu N(1 -\frac{N}{K}) \tag{17}
\end{align*}
```

separate variables

```math
\begin{align*}
\frac{\mathrm{d} N}{N(1 -\frac{N}{K})}=\mu \mathrm{d} t   \tag{18}
\end{align*}
```

Our next goal would be to integrate both sides of this equation, but the form of the right hand side doesn't look elementary and will require a partial fractions expansion. That is, we wish to write 


```math
\begin{align*}
\frac{1}{N(1 -\frac{N}{K})} = \frac{A}{N}+\frac{B}{1 -\frac{N}{K}} \tag{19}
\end{align*}
```

where $A$ and $B$ are unknown constants. If we multiply on the left and right hand sides by  $N \left( 1- \frac{N}{K} \right)$ (which is equivalent to putting the right hand side over a common denominator) we arrive at the equation 

```math
\begin{align*}
1 = A \ \left( 1 -\frac{N}{K}\right) + B \cdot N  = A + N (B - \frac{A}{K}) \tag{20}
\end{align*}
```

Since there is no term with $N$ on the left hand side, we see that 

```math
\begin{align*}
B - \frac{A}{K} = 0 \quad \mbox{ or } \quad B = \frac{A}{K} \tag{21}
\end{align*}
```

If we set $B = \frac{A}{K}$ then we are left with $A=1$, and thus the partial fraction decomposition is 

```math
\begin{align*}
\frac{1}{N(1 -\frac{N}{K})} = \frac{1}{N}+\frac{\frac{1}{K}}{1 -\frac{N}{K}} \tag{22}
\end{align*}
```

so the integral becomes

```math
\begin{align*}
\frac{\mathrm{d} N}{N}+\frac{\frac{\mathrm{d} N}{K}}{1 -\frac{N}{K}}=\mu \mathrm{d} t   \tag{23}
\end{align*}
```

the first part is simply:

```math
\begin{align*}
\int\frac{dN}{N} = \ln (N) \tag{24},
\end{align*}
```

For the second term, we must use a substitution $u=1-\frac{N}{K}$, which gives a differential  $du = \frac{-1}{K} \ dN$. Thus we may write the second term on the right hand side as:

```math
\begin{align*}
\int \frac{dN/K}{\left( 1-\frac{N}{K} \right)} = \int \frac{-du}{u} = -\ln (u) = -\ln (1-N/K) \tag{25}
\end{align*}
```

Putting all these terms together gives us:

```math
\begin{align*}
\mu t + c = \ln (N)-\ln (1-\frac{N}{K}) = \ln \left[\frac{N}{1-N/K} \right] \tag{26}
\end{align*}
```

Here we have used the property of logarithms to equate the difference of the logs with the log of the quotient. The additional term, $c$, on the left hand side is the free constant of integration, which will be determined by considering initial conditions to the differential equation. Exponentiating both sides of the equation gives 

```math
\begin{align*}
e^{\mu t + c} = \frac{N}{1-\frac{N}{K}} \tag{27},
\end{align*}
```

so

```math
\begin{align*}
 e^{\mu t} e^c = \frac{N}{1-\frac{N}{K}} \tag{28}
\end{align*}
```

```math
\begin{align*}
e^{\mu t} C = \frac{N}{1-\frac{N}{K}} \tag{29}
\end{align*}
```

to find $C$ we use teh inital condition that $N(t=0)=N_0$, and substituting gives 

```math
\begin{align*}
C = C e^0 = \frac{N_0}{1-\frac{N_0}{K}} = \frac{N_0}{1-\frac{N_0}{K}} \frac{K}{K} =\frac{K N_0}{K-N_0} \tag{30}
\end{align*}
```

Solving now for $P$, we first cross-multiply to arrive at 

```math
\begin{align*}
\left(1-\frac{N}{K} \right) C e^{\mu t} = N \tag{31}
\end{align*}
```

and putting all terms including $N$ on one side of the equation, 

```math
\begin{align*}
C e^{ \mu t} = N \left[1 + \frac{C e^{ \mu t}}{K} \right]  \tag{32}
\end{align*}
```

Solving now for $N$, 

 
```math
\begin{align*}
N = \frac{C e^{\mu t}}{1 + \frac{C e^{\mu t}}{K} } =
\frac{\frac{K \cdot N_0}{K-N_0} e^{\mu t}}{1 + \frac{\frac{K\cdot N_0}{K-N_0} e^{\mu t}}{K} }\tag{33}
\end{align*}
```

Simplifying this expression by multiplying numerator and denominator by  $(K-N_0) e^{-\mu t}$ gives 

```math
\begin{align*}
N = \frac{K N_0}{N_0 +(K-N_0) e^{- \mu t}} \tag{34}
\end{align*}
```

if we monitor hwo the value of $\mu$ changes 

```math
\begin{align*}
(K-N_0) e^{- \mu t} &=\frac{K N_0}{N} - N_0\\
(K-N_0) e^{- \mu t} &=\frac{N_0 (K - N)}{N}\\
 e^{- \mu t} &=\frac{N_0 (K - N)}{N (K-N_0)}\\
 - \mu t &= Log(\frac{N_0 (K - N)}{N (K-N_0)})\\
  \mu  &= \frac{1}{t}Log(\frac{N (K-N_0)}{N_0 (K - N)})\\
    \frac{Log 2}{T}  &= \frac{1}{t}Log(\frac{N (K-N_0)}{N_0 (K - N)})\\
    T  &=  t Log 2  Log(\frac{N_0 (K - N)}{N (K-N_0)})\\
\end{align*}
```
"

# ╔═╡ 61fe29b8-065c-430f-a28d-b463780f8b00
md"Fitness can fundamentally be achieved by two different strategies: long life (stability) or fast reproduction (multiplication, replication). These strategies are to some degree dependent: since no organism is immortal, a minimum amount of reproduction is needed to replace the organisms that have died; yet, in order to reproduce, the system must live long enough to reach the degree of development where it is able to reproduce. On the other hand, the two strategies cannot both be maximally pursued: the resources used for fast reproduction cannot be used for developing a system that will live long, and vice-versa. This means that all evolutionary systems are confronted with a development-reproduction trade-off: they must choose whether they invest more resources in the one or in the other.
How much a given system will invest in one strategy at the expense of the other one depends on the selective environment. In biology, this is called r-K selection: in an r-situation, organisms will invest in quick reproduction, in a K-situation they will rather invest in prolonged development and long life. Typical examples of r-species are mice, rabbits, weeds and bacteria, which have a lot of offspring, but a short life expectancy. Examples of organisms undergoing K-selection are tortoises, elephants, people, and sequoia trees: their offspring are few but long-lived. In summary, r-selection is selection for quantity, K-selection for quality of offspring."

# ╔═╡ a0bdad86-6cbd-4da1-b1ad-3c1f5e3c1ef3
md"r-organisms	K-organisms
short-lived	long-lived
small	large
weak	strong or well-protected
waste a lot of energy	energy efficient
less intelligent, experienced...	more intelligent, experienced...
have large litters	have small litters
reproduce at an early age	reproduce at a late age
fast maturation	slow maturation
little care for offspring	much care for offspring
strong sex drive	weak sex drive
small size at birth	large size at birth"

# ╔═╡ 7f30ad4d-9964-4140-8003-051653d5f1e4
md"
# Differentiation dynamics

what happens when we analyze the dynamcis of a population of cells that not only proliferates, but also differentiates into another type of cell? For instance, a population of stem cells, in a developing organ. We assume a developing organ as a population of cycling progenitors 'P' that cycle with an average cell cycle $T$. Some of these cells terminally differentiate, exit the cell cycle and acquire a given specialized phenotype 'D'. We start from a initial population of progenitors $P_0$ and differentiated $D_0$ cells. A common approac is top characterize the dynamics of the population focusing on the outcome of the cell division of the progenitors. In principle, each division of a 'P' cell can give two progenitors ('pp' division), two differentiated cells ('dd' division) and also an assymetric mode of divisin where a progenitor and a differentiated cell is generated ('pd' division). If we calculate the average amount of 'P' and 'D' generated after a single cell cycle (n=1) we can write the number of progenitor and differentiated cells as:

$$
\begin{eqnarray}
P_1&=&P_{0}(2pp+pd) \tag{35}\\
D_1&=&D_{0}+P_{0}(2dd+pd)\tag{36}
\end{eqnarray}
$$

where using the condition $pp+pd+dd=1$, 

$$
\begin{eqnarray}
P_1&=&P_{0} (1+pp-dd)\tag{37}\\
D_1&=&D_{0}+P_{0}(1+dd-pp)\tag{38}
\end{eqnarray}
$$

Therfore, for n=2,

$$
\begin{eqnarray}
P_2&=&P_{1} (1+pp-dd)\tag{39}\\
D_2&=&D_{1}+P_{1}(1+dd-pp)\tag{40}
\end{eqnarray}
$$

applying eqs. 37 and 38, we obtain

$$
\begin{eqnarray}
P_2&=&P_{0}(1+pp-dd)(1+pp-dd)=P_{0} (1+pp-dd)^2\tag{41}\\
D_2&=&D_{0}+P_{0}(1+dd-pp)+P_{0}(1+pp-dd)(1+dd-pp)\tag{42}
\end{eqnarray}
$$

and rearranging terms in eq. 42:
$$
\begin{eqnarray}
D_2&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd))\tag{43}\\
\end{eqnarray}
$$

Subsequently, for n=3
$$
\begin{eqnarray}
P_3&=&P_{2} (1+pp-dd)\tag{44}\\
D_3&=&D_{2}+P_{2}(1+dd-pp) \tag{45}
\end{eqnarray}
$$


applying eqs. 41 and 42, we obtain
$$
\begin{eqnarray}
P_3&=&P_{0} (1+pp-dd)^2 (1+pp-dd)= P_{0} (1+pp-dd)^3 \tag{46}\\
D_3&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd) + (1+pp-dd)^2) \tag{47}
\end{eqnarray}
$$

therefore, for $n$ steps, we obtain, 

$$
\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^n \tag{48}\\
D_n&=&D_{0}+P_{0}(1+dd-pp)(1+(1+pp-dd)+(1+pp-dd)^2+...+(1+pp-dd)^{n-1})\tag{49}
\end{eqnarray}
$$

where the second term in eq. 49 can be written as
$$
\begin{eqnarray}
1+(1+pp-dd)+(1+pp-dd)^2+...+(pp-dd)^{n-1}=\displaystyle\sum_{i=0}^{n-1} (1+pp-dd)^i\tag{50}
\end{eqnarray}
$$

which renaming $r=1+pp-dd$ is equivalent to 
$$
\begin{eqnarray}
\displaystyle\sum_{i=0}^{n-1} r^i=\frac{1-r^n}{1-r}\tag{51}
\end{eqnarray}
$$

therefore, eq. 50 can be written as
$$
\begin{eqnarray}
D_n&=&D_{0}+P_{0}(1+dd-pp)\frac{1-(1+pp-dd))^n}{1-(1+pp-dd)} \tag{51}
\end{eqnarray}
$$

which, after simplifying terms, can be rewritten as
$$
\begin{eqnarray}
D_n&=&D_{0}+P_{0}(1-(1+pp-dd)^n)\frac{1+dd-pp}{dd-pp}\tag{52}
\end{eqnarray}
$$

and taking into account eq. 48, we obtain the final equation for the number of progenitors 'P' and differentiated 'D' cells in a stem cell population that is growing and differentiating:
$$
\begin{eqnarray}
D_n&=&D_{0}+(P_{0}-P_{n})\frac{1+dd-pp}{dd-pp}=D_{0}+(P_{n}-P_{0})\frac{1+dd-pp}{pp-dd} \tag{53}
\end{eqnarray}
$$

Interestingly, both equations depend on the iteration step $n$ only via the number of progenitors at a given time in the system $P_n$. 

For simplicity, the system of equations has been derived for a situation of discrete $n=\Delta t/T$, (n=1,2,3...), i.e., with the time step equal to the average cell cycle $\Delta t=T$. If we instead consider the time step as half of the cell cycle ($\Delta t=T/2$), then $n=2 \Delta t/T$, (n=1,2,3...), and eq. 37 is now:
$$
\begin{eqnarray}
P_1&=&P_{0} (1+pp-dd)^{\frac{1}{2}}\tag{54}
\end{eqnarray}
$$

and following identical iteration steps we arrive at
$$
\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^{\frac{n}{2}}\tag{55}
\end{eqnarray}
$$

while eqs. for 'D' cells remain the same, since it does not depend explicitly on the iteration step. This way, for $\delta t=1$, and following the same process we obtain
$$

\begin{eqnarray}
P_n&=&P_{0} (1+pp-dd)^{\frac{n}{T}}\tag{56}
\end{eqnarray}
$$

which can be generalize for the continuum limit $n=t$,

$$
\begin{eqnarray}
P_{t}&=&P_0 (1+pp-dd)^{\frac{t}{T}}\tag{57}\\
D_{t}&=&D_{0}+(P_{t}-P_{0})\frac{1+dd-pp}{pp-dd}\tag{58}
\end{eqnarray}
$$

Finally, f we rewrite  $\Delta$$P=P_{t}-P_{0}$, we obtain the expression. 
$$
\begin{eqnarray}
P_{t}&=&P_0 (1+pp-dd)^{\frac{t}{T}}\tag{59}\\
D_{t}&=&D_{0}+\Delta P\frac{1+dd-pp}{pp-dd}\tag{60}
\end{eqnarray}
$$
"


# ╔═╡ 435e892c-e56d-4f20-be01-cbf462194882
pp=0.6
dd=0.2
P₀=100
D₀=50
T=24
t=collect(0:0.1:100)
P1=plot(t,t->P₀*(1+pp-dd)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd)^(t/T))-1)*((1-pp+dd)/(pp-dd)),label="D",seriestype=:line,ylims = (0,400))

pp=0.0001
dd=0.00001
P2=plot(t,t->P₀*(1+pp-dd)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd)^(t/T))-1)*((1-pp+dd)/(pp-dd)),label="D",seriestype=:line,ylims = (0,400))

pp=0.0001
dd=0.3
P3=plot(t,t->P₀*(1+pp-dd)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd)^(t/T))-1)*((1-pp+dd)/(pp-dd)),label="D",seriestype=:line,ylims = (0,400))

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ ae184b0d-f61d-4aea-853b-17cb597d1087
md"
to do:
    
    it will be cool to include a term of saturations as the logistic equation does. 
    
"

# ╔═╡ aecc1472-9167-4b09-a08c-9e701def7d54
md"The true power of these equatiosn is that they are analytical, in the sense that now we turn the equations down to see if we can predict the correct value of pp-dd and T"

# ╔═╡ 8e867da7-e8bc-4dca-8fe4-644f5277cb67
pp=0.6
dd=0.2
P₀=100
D₀=50
T_true=24
t=collect(0:0.1:100)

P=P₀.*(1+pp-dd).^(t./T_true)
D=D₀.+P₀.*(((1 .+pp-dd).^(t./T_true)).-1).*((1-pp+dd)/(pp-dd))

P1=plot(t,P,label="P",seriestype=:line,ylims = (0,4300))
plot!(t,D,label="D",seriestype=:line,ylims = (0,400))

P_=P[2:end]
D_=D[2:end]
t_=t[2:end]
P__=P[1:end-1]
D__=D[1:end-1];
t__=t[1:end-1];
gamma=1
pp_dd=(P_ .-P__) ./(P_ .-P__ .+D_ .-D__);
T=(t_ .-t__) .* log.(1 .+(gamma.*(abs.(pp_dd)))) ./ abs.(log.(P_ ./P__));
P2=plot(t_,pp_dd, ylims = (-1,1))
hline!([pp-dd])
P3=plot(t_,T,ylims = (23,25))
hline!([T_true])

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ d46dfe3b-f178-4d18-b95f-c961bc207466
pp=0.1
dd=0.101
P₀=100
D₀=50
T_true=20
t=collect(0:0.1:100)

P=P₀.*(1+pp-dd).^(t./T_true)
D=D₀.+P₀.*(((1 .+pp-dd).^(t./T_true)).-1).*((1-pp+dd)/(pp-dd))

P1=plot(t,P,label="P",seriestype=:line)
plot!(t,D,label="D",seriestype=:line)

P_=P[2:end]
D_=D[2:end]
t_=t[2:end]
P__=P[1:end-1]
D__=D[1:end-1];
t__=t[1:end-1];
gamma=1
pp_dd=(P_ .-P__) ./(P_ .-P__ .+D_ .-D__);
T=(t_ .-t__) .* log.(1 .+(gamma.*(abs.(pp_dd)))) ./ abs.(log.(P_ ./P__));
P2=plot(t_,pp_dd, ylims = (-1,1))
hline!([pp-dd])
P3=plot(t_,T,ylims = (0,25))
hline!([T_true])

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ ca6a3577-ac75-4993-b2e5-4975ef4aaf9f
md"The equations for the cell cycle start to fail whene we go to values lower than pp-dd=0"

# ╔═╡ 349b1edf-3956-4b47-9d2f-760df409d40e
pp=-0.1
dd=0.0
P₀=100
D₀=50
T_true=20
t=collect(0:0.1:100)

P=P₀.*(1+pp-dd).^(t./T_true)
D=D₀.+P₀.*(((1 .+pp-dd).^(t./T_true)).-1).*((1-pp+dd)/(pp-dd))

P1=plot(t,P,label="P",seriestype=:line)
plot!(t,D,label="D",seriestype=:line)

P_=P[2:end]
D_=D[2:end]
t_=t[2:end]
P__=P[1:end-1]
D__=D[1:end-1];
t__=t[1:end-1];
gamma=1
pp_dd=(P_ .-P__) ./(P_ .-P__ .+D_ .-D__);
T=(t_ .-t__) .* log.(1 .+(gamma.*(abs.(pp_dd)))) ./ abs.(log.(P_ ./P__));
P2=plot(t_,pp_dd, ylims = (-1,1),label="pp-dd predicted",line = (:blue, 1))
hline!([pp-dd],label="pp-dd true",line = (:red,:dash, 1))
P3=plot(t_,T,ylims = (0,25),label="T predicted",line = (:blue, 1))
hline!([T_true],label="T true",line = (:red,:dash, 1))

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ 3e8066bd-d8dc-47d9-98e5-0b8b0759e179
md" ### Including Apoptosis

what if we include apoptosis"

# ╔═╡ 2e9621c6-50c4-4d30-84f5-f26ea707a808
pp=0.6
dd=0.2
ø=0.1
P₀=100
D₀=50
T=24
t=collect(0:0.1:100)
P1=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))



pp=0.0001
dd=0.00001
ø=0.1
P2=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))


pp=0.0001
dd=0.3
ø=0.1
P3=plot(t,t->P₀*(1+pp-dd-ø)^(t/T),label="P",seriestype=:line,ylims = (0,4300))
plot!(t,t->D₀+P₀*(((1+pp-dd-ø)^(t/T))-1)*((1-pp+dd-ø)/(pp-dd-ø)),label="D",seriestype=:line,ylims = (0,400))
plot!(t,t->P₀*(((1+pp-dd-ø)^(t/T))-1)*(ø/(pp-dd-ø)),label="ø",seriestype=:line,ylims = (0,400))

plot(P1,P2,P3,layout=(1,3),legend=true,size = (800, 500))

# ╔═╡ 65983c74-85d8-11ec-0c78-c9585d7356d8
md" # 3. Computational Models

Complex systems are described by nonlinear equations, and therefore, very often, they cannot be solved analitically. Therefore, the most common aproach is the use of mathematical models. 

What type of models do we want? 

Models with all parameters and interactions, or simplified models?

_Models are more useful when they are wrong_ 

Discrete versus continuous systems

Logistic growth (Discrete)

Non-linear equations and steady-state in population growth: limited growth by available nutrients 

How to make a numerical model of a complex system. Because complex systems are basically composed of parts that are interacting, there is a clear analogy with  chemical reactions (composed of molecules that are reacting). Therefore, we can borrow many concepts from this field and applied them directly to the study of Complex Systems, such as: 

- Mass action Law
- Mass conservation law
- Chemical Equilibrium
- Stoichimestry 

## 3.1 Basics of Chemical Dynamics

Biological systems are composed of multiple entities (cell, molecuels, genes, organisms...) that interact at multiple time and spatial scales. We will study first how to deal mathematically with systems. Although not limited to chemistry, this is what the field of Chemical Dynamics has been perfecting for many years. Therefore, we will first borrow some ideas from Chemistry that will sound quite familar to some of the students of this course (but not to others). As a start, let's assume a system of interacting species (chemicals, cells, genes, proteins, organisms...) `A`, `B`, `C` and `D` that interact between them, so the amount of each one is allowed to change over time: "

# ╔═╡ 547db9e8-8de9-408d-80d4-48947b4aa1c2
md"Some possible schemes of interaction are, for instance:
```math
\begin{align*}
A + B  &\overset{k_1}{\longrightarrow} C  \tag{1} \\
2A + B  &\overset{k_2}{\longrightarrow} C  \tag{2} \\ 
C  &\overset{k_3}{\longrightarrow} A + B   \tag{3} \\ 
AB + C  &\overset{k_4}{\longrightarrow} AC + B   \tag{4} \\ 
AB + CD  &\overset{k_5}{\longrightarrow} AC + BD   \tag{5} \\ 
\end{align*}
```
"

# ╔═╡ 845a5477-29aa-4ebd-89a0-53ddeed9343d
md" which is the typycal notation of a chemical reaction, but it can be extrapolated to any other system where some entities interact. We know from chemistry that the speed of a chemical reaction (i.e., the change in the amount of a given species) is proportional to the amount of collisions (i.e., the amount of times that two given reactants collide into each other in the solution). It is easy to realize that this rate of collisions has to be proportional to the amount of reactants in the solution. For instance, let's illustrate this situation for the the first reaction:"

# ╔═╡ 5d00ba37-89a4-44a1-b8dc-c47394890b1f
begin
	dog_slide = @bind 🐶 html"<input type=range min=1 max=50 step=1>"
	cat_slide = @bind 🐱 html"<input type=range min=1 max=50 step=1>"
	
	md"""
	**How many molecules do you have?**
	
	Molecules of type A: $(dog_slide)
	
	Molecules of type B: $(cat_slide)
	"""
end

# ╔═╡ c4bb7c6a-1165-44d6-8c1e-ea6798297c72
begin
	p1 =scatter(fill(rand(🐶),1), fill(rand(🐶),1), fill(rand(🐶),1), m=(10, 0.8, :blues),framestyle=:box, title="Collisions",label="A")
	p1 =scatter!(fill(rand(🐱),1), fill(rand(🐱),1), fill(rand(🐱),1), m=(10, 0.8, :reds),framestyle=:box,label="B")

end

# ╔═╡ d4e5e672-4759-47ed-9a69-4d92efd85b05
md"It is clear that there will be more collisions in the situation on the right, so the reaction will take place faster. So, not only you will end up with more [C], you will also obtain [C] faster. Based on these basic facts, the speed of each reaction above can be matematically written as:

```math
\begin{align*}
speed~of~reaction~1~&\propto A \cdot B \tag{6} \\ 
speed~of~reaction~2~ &\propto A \cdot A \cdot B \tag{7} \\ 
speed~of~reaction~3 ~ &\propto C \tag{8} \\ 
speed ~of ~reaction~4 ~&\propto AB \cdot B    \tag{9} \\ 
speed ~of ~reaction~5 ~&\propto AB \cdot CD  \tag{10} \\ 
\end{align*}
```


and the proportionality constant is simply the value of the reaction rates $k_i$, which is a measure of how many of the collisions are effective, so:

```math
\begin{align*}
speed~of~reaction~1~&= k_1 \cdot A \cdot B \tag{11} \\ 
speed~of~reaction~2~ &= k_2 \cdot A \cdot A \cdot B \tag{12} \\ 
speed~of~reaction~3 ~ &= k_3 \cdot C \tag{13} \\ 
speed ~of ~reaction~4 ~&= k_4 \cdot AB \cdot B    \tag{14} \\ 
speed ~of ~reaction~5 ~&= k_5 \cdot AB \cdot CD  \tag{15} 
\end{align*}
```
"

# ╔═╡ 6b25dc03-352b-4351-b387-2c3d0dd80106
md" ### 1.1. Equilibrium 

In these system, the equilibrium (defined as the condition of no change in the amount of the species interacting) corresponds to the situation when the speed of the reaction is zero. In these type of __irreversible__ reactions, equilibrium occurs when the reaction is finished because one of the reactants (the limitant) has been fully consumed. On the contrary, in __reversible__ reactions, the equilibrium can be dynamic and does not mean explicitely that the reaction is stopped. For instance, let's asumme a very simple equilibrium between two species:

```math
\begin{align*}
2 NO_2   &\overset{k_1}{\underset{k_2}{\longleftrightarrow}} N_{2}O_4 \tag{16} \\
\end{align*}
```

This system is actually composed of two (not so) different reactions, 

```math
\begin{align*}
2 NO_2   &\overset{k_1}{\longrightarrow}  N_{2}O_4 \tag{17} \\
N_{2}O_4  &\overset{k_2}{\longrightarrow} 2 NO_2 \tag{18}
\end{align*}
```
and the speed of each reaction is:

```math
\begin{align*}
speed~of~reaction~16~&= k_1 \cdot NO_2 \cdot NO_2 = k_1 \cdot NO_{2}^{2} \tag{19} \\ 
speed~of~reaction~17~ &= k_2 \cdot N_{2}O_4 \tag{20}
\end{align*}
```

"

# ╔═╡ eb9eb1de-7b0f-4045-b438-ab7930243e5c
md"

For this reversible reaction, it is impossible to have pure $NO_2$ or  $N_{2}O_4$, (i.e., as soon as the amount of one species is approaching zero, the force towards the other direction in the reaction approaches infinite). In other words, the less of a molecule exists in solution, the faster will be generated. Therefore, equilibrium cannot be defined as the point where the reaction is finished (because it never finishes), or the point where one of the species has been fully consumed. In these conditions, the equilibrium is better defined as the situation when the speed of the two reactions is equal. 

```math
\begin{align*}
k_1 \cdot [NO_{2}]^{2}_{eq} = k_2 \cdot [N_{2}O_4]_{eq}  \tag{21} 
\end{align*}
```

In other words, equilibrium occurs when the concentration of the reactants do not change overtime. In our example, it means that you reach a value of $[N_{2}O_4]$ and $[NO_2]^2$ that is constant. As a consequence, in reversible reactions at equilibrium, the ratio $\frac{[N_{2}O_4]_{eq}}{[NO_2]^2_{_{eq}}}$ is a constant value that is proportional to the ratio between the reaction rates $k_1$ and $k_2$. 

```math
\frac{k_1}{k_2} = \frac{[N_{2}O_4]_{eq}}{[NO_{2}]^{2}_{eq}} \tag{23}\\
```

In consequence, that the ratio between the concentrations of reactants at equilibrium ($\frac{[N_{2}O_4]_{eq}}{[NO_2]^2_{eq}}$ in our example) does not depend on how much $[N_{2}O_4]$ or $[NO_2]^2$ you put intially in the system. This means that if does not matter if you start with zero concetration of one the reactants and tons of molecules of the other: at qeuilibrium the ratio between teh concentrations only depends on the ratio between the rates of the two reversible equations. 
"

# ╔═╡ 07735905-642b-4141-87fe-c9bc20ae04c2
md"Now let's practice with another example of chemical reaction: 

```math
Na_{2}CO_3 + CaCl_2 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} CaCO_3 + 2 \cdot NaCl \tag{22}
```

Following the rationale of the previous chemical reaction, we can write the same correspondence between 

```math
\frac{k_1}{k_2}
= \frac{[CaCO_3][NaCl]^2}{[Na_{2}CO_3][CaCl_2]} \tag{23}\\
```
"

# ╔═╡ 76691838-7086-4ca0-91d1-b5ef2c3e3b24
md"### 1.2 Order of reactions and Equilibrium Constant 

The order of a reaction refers to the power dependence of the rate on the concentration of each reactant. In brief, the order of the reaction indicates the correlation of its velocity with the amount of reactants.

- For a zero-order reaction, the rate does not depend on the concentration of any species.
- For a first-order reaction, the rate is dependent on the concentration of a single species.
- For a second-order reaction, the rate is dependent on the square of the concentration of a single reactant, or two reactants.

This way, the units of the rate constants will depend on the type of reaction taking place. 

A common characterization of a system of interacting species is the equilibrium constant $K_{eq}$, i.e, the ratio between the two reversible reaction rates. 

```math
\begin{align*}
K_{eq}=\frac{k_1}{k_2} \tag{22} 
\end{align*}
```

Here, we can ask ourselves what are the units of the equilibrium constant. From the equation above, since it is defined as the ratio between two kinetic constants of a reversible reaction, we can conclude that it has no dimensions. The correct aswer is that, since the units of the kinetic constants depend on the order of the reaction taking place, and that the two reactions that form a reversible reaction can have different orders, the units of the equilibrium constant will depend on each particular system.  

In this direction, we cannot compare the dynamics and the equilibrium state of two reactions simply based on the value of their kinetic constant. To illustrate this, we run below the numerical simulation of two different reactions that have the same value the kinetic constant. Since the two reversible reactions are formed by reactions of different order, the units of the kinetic constants are different, therefore, they cannot be compared.
"

# ╔═╡ 23c143d3-29a3-4824-877e-1f13d0818ab6
simpleODE1! = @ode_def abetterway begin
  da = -k1 * a + k2 * c 
  dc =  k1 * a - k2 * c 
    end k1 k2

# ╔═╡ 71777cd5-eedb-4dbf-a3a0-84ebe31ef602
simpleODE2! = @ode_def abetterway2 begin
  da = -k1 * a * b + k2 * c
  dc = k1 * a * b - k2 * c 
  db = -k1 * a * b + k2 * c 
    end k1 k2

# ╔═╡ aed09f8b-eb07-4a34-a0bf-dc9ac987e85f
begin
	b_slide = @bind b₀ html"<input type=range min=0.0 max=2 step=0.1>"

	md"""
	**Move the silder to change the initial concentration of b?**
	
	Concentracion inicial de b: $(b_slide)
	
	"""
end

# ╔═╡ 0c575ced-ef89-4bf8-8f8c-c784d7f2d33a
begin
		k1=2.3e-1;  # units 1/(Ms)
		k2=2.5e-1;  # units 1/(M M s)
		tspan = (0.0,10.0)
		p = (k1,k2)
		a₀=0.05; # units (M)
		#b₀=0.05; # units (M)
		c₀=0.00; # units (M)
		u₀=[a₀,c₀];
	
	prob = ODEProblem(simpleODE1!,u₀,tspan,p)
	sol1 = solve(prob);
	p3=plot(sol1,xlabel="Time [s]",ylabel="Concentration [M]",title="a <-> c",ylims = (0,0.05));
	
	u₀=[a₀,c₀,b₀];
	
	prob = ODEProblem(simpleODE2!,u₀,tspan,p)
	sol2 = solve(prob);
	p4=plot(sol2,xlabel="Time [s]",ylabel="Concentration [M]",title="a + b <-> c, b₀ = $b₀",ylims = (0,0.05));
	
	
	plot(p3,p4,layout=(1,2),legend=true)
end

# ╔═╡ 7b95dbaf-dc9d-45eb-9a0a-f2f75daf8765
md" 
You can see above two examples with the same value of the equilibrium constant, but completely different dynamics. The left one corresponds to a chemical reaction where forward and backward reaction have the same order, and therefore $K_{eq}$ is a nondimensional parameter. The right one corresponds to a chemical reaction where forward and backward reaction have different order, and the units of $K_{eq}$ is now [M]. You can move the slider to change the intial concentration of b to find a value where the two reactions have a similar dynamics. 





"

# ╔═╡ 2e15f653-b17b-4724-8fa9-a10606093c5d
md" 
## 2. General Formulation for a system of interacting entities
Let's now introduce a general notation for any interaction scheme between species. Let's define a system where species `A` and `B` react reversibly to give species `C` and `D`:

```math
\begin{align*}
aA + bB  &\overset{k_1}{\underset{k_2}{\longleftrightarrow}} cC + dD \tag{26} \\ 
\end{align*}
```

where `a`, `b`, `c`, `d` correspond to the stoichiometric coefficients for a balanced interaction. At any instant in time, we can define a ratio between the amounts of each species, such as:

```math
Q= \frac{[C]^c[D]^d}{[A]^a[B]^b} \tag{27}
```

where `Q` is defined in chemistry as the reaction quotient, and measures the relative amounts of the interacting species present during a reaction at a particular point in time. As the time evolves, the system moves towards its equilibrium, and the value of `Q` gradually approaches to the equilibrium constant $K_{eq}$. The general expression of this equilibrium constant is:

```math
K_{eq}= \frac{[C]^c[D]^d}{[A]^a[B]^b} \tag{28}
```

"


# ╔═╡ cd434441-2358-4d50-8744-98af4fc99176
md"### 3. Independence of ratio between equilibrium concentrations on initial conditions: 

We have seen how, in reversible reactions, the ratio between the concentrations of the reactants at equilibrium depends only on the ratio between the kinetic constants, and therefore it is __independent on the initial concentrations__ of the reactants. This apparent simple result has important implications, and therefore, it is important to think a bit about this feature. To illustrate that, we will show the solution of the system for two different initial conditions, and compare the dynamics and the value of ``Q(t)``."

# ╔═╡ 7a252b5a-13f1-489c-8dfd-ef325adeee56
simpleODE3! = @ode_def abetterway3 begin
  da = -k1 * a * b + k2 * c * d^2
  db = -k1 * a * b + k2 * c * d^2
  dc = k1 * a * b - k2 * c * d^2
  dd = 2 * k1 * a * b - 2 * k2 * c * d^2
    end k1 k2

# ╔═╡ 0a62b4ce-861b-44b2-bb56-c6f4cfce5fef
begin
	intial_a = @bind aa₀ html"<input type=range min=0.01 max=0.1 step=0.01>"
	
	intial_b = @bind bb₀ html"<input type=range min=0.01 max=0.1 step=0.01>"
	intial_c = @bind cc₀ html"<input type=range min=0.01 max=0.1 step=0.01>"
	intial_d = @bind dd₀ html"<input type=range min=0.01 max=0.1 step=0.01>"
	
	md"""
	**How many molecules do you have?**
	
	Initial concetration of a: $(intial_a)
	
	Initial concetration of b: $(intial_b)
	
	Initial concetration of c: $(intial_c)
	
	Initial concetration of d: $(intial_d)
	
	"""
end

# ╔═╡ e48547d5-f638-4885-8d2d-a8f70ad67dd5
begin

	prob3 = ODEProblem(simpleODE3!,[aa₀,bb₀,cc₀,dd₀],(0.0,50.0),(1.3e0,2.5e0))
	sol3 = solve(prob3);

	p6=plot(sol3,xlabel="Time [s]",ylabel="Concentration [M]",title="Dynamics",ylims = (0,0.2));

	Q1=(sol3[3,:].*sol3[4,:].^2)./(sol3[1,:].*sol3[2,:])
	p7=plot(sol3.t,Q1,title= "Quotient coefficient",xlabel="Time [s]",ylabel="Q [M]",ylims = (0,0.6))

	plot(p6,p7,layout=(1,2),legend=true)
end

# ╔═╡ ea66ecaf-96c5-4b54-92eb-55d0eda5a001
md"We see that the value of $Q$ does not change, and always aproaches the same final value at the end of the reaction, and this value is $K_{eq}$.

This apparently peculiar relationship between the amounts of reactants and products in an equilibrium (no matter how many reactants you start with) is based on the previous notion of dynamic equilibrium (two opposite reactions with the same speed). The values of `[A]`, `[B]`, `[C]` and `[D]` represent the amount of each species at equilibrium. If this amount is given as a concentration, lets say in moles/volume, `M`, the units of the equilibrium constant $K_{eq}$ are:



```math
[K_{eq}]=M^{c+d-a-b} \tag{29}
```
"

# ╔═╡ 442944b8-8f11-484e-988c-6b1bcbd13c8e
@htl("""

<div class='blue-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Computer Task 1: Working with functions"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

""")

# ╔═╡ 10d1f8cc-1a72-420b-b5cc-0286105850b3
md"A good way to write code to solve numerical models is to pack the scripts in separate functions that can be called with a single sentence. Modern programming languages are optimized to work with functions, resulting in much more efficient code. 

As a first task, write a simple computer program with two functions that take as input arguments, a vector of concentrations at equilibirum of two reactions reaction 16 and 22 and return their value $K_{eq}$ with the correct units"

# ╔═╡ 6f8ff973-6845-4517-8bca-1bec97fa4edf
@htl("""

<div class='red-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Solution Computer Task 1"

</script>

<style>
.red-background {
	padding: .5em;
	background: lightgreen;
	color: black;
}
</style>

""")

# ╔═╡ 313f484e-bea2-437d-ba91-d766d36dd248
function Calculate_Keq1(a)
	NO2=a[1]
	N2O4=a[2]
	NO2/N2O4^2
end

# ╔═╡ f8ebbcd7-61e6-46bd-91fd-546ba931f76c
function Calculate_Keq2(a)
	Na2CO3=a[1]
	CaCl2=a[2]
	CaCO3=a[3]
	NaCl=a[4]
	(CaCO3*NaCl^2)/(Na2CO3*CaCl2)
end

# ╔═╡ 7c053891-043a-4c0d-9afd-d03012d102c8
NO2=2; N2O4=3;

# ╔═╡ 5959c966-ff6e-4172-9fad-c05258050ff4
Calculate_Keq1([NO2,N2O4]) 

# ╔═╡ 29001c0e-9f2e-4fbb-ba8f-3b34de1d2a4b
Na2CO3=2;CaCl2=0.5;CaCO3=2;NaCl=1.2;

# ╔═╡ 10e57e8a-bc19-4abe-9b62-ba7b65c69c9d
Calculate_Keq2([Na2CO3,CaCl2,CaCO3,NaCl])

# ╔═╡ 6de3e149-e3e4-4c24-975b-e1e2259392f3
md"For reaction 16 we have:

For reaction 21 we have:

```math
[K_{eq}]=M^{1+0-2-0}=M^{-1}  \tag{32}\\
```

while for reaction 21 we have:

```math
[K_{eq}]=M^{1+2-1-1}=M  \tag{33}\\
```


"

# ╔═╡ 1525b27e-2613-4f27-9f13-171e24cd574e
md" ### 1.3 Effect of temperature 

We have seen that the speed of a reaction is proportional to the concentration of reactants, and to the kinetic rate constant. Another way to increase or decrease the speed of a given reaction is via changes in the _temperature_. At the molecular level, temperature is related to the random motions of the particles in matter. In other words, temperature in a solution is a measure of the average kinetic energy of the molecules involved, and therefore, changes in the temperature will result in changes in the number of collisions in a chemical reaction. 

The dependence between temperature and reaction rate constant is set by the The Arrhenius Equation:
```math
k = A e^{-E_a/RT} 
```
The parameter of $E_a$ is the activation energy, and $A$ is a parameter that relates to the frequency of collisions and the orientation of a favorable collision probability

The activation energy is the energy required for two molecules to interact, it is like a barrier of Energy. Lets plot below the dependence of the rate constant with the temperature. 

"

# ╔═╡ f55ef39d-192e-4e97-acf5-a5fad4c6bae8
begin
	E_slide = @bind E_a html"<input type=range min=1.0 max=20 step=1>"

	
	md"""
	**Move the silder to change the Activation Energy**
	
	Ea: $(E_slide)
	
	"""
end

# ╔═╡ 59fd13d8-e3bd-4756-a77f-a571876670a3
md"we can see that the kinetic rate constant increases with the temperature (more kinetic energy of the particles, will mean more efficient collisions). Also, as we increase $E_a$, the rate decreases, suggesting that the number of efficient collisions is reduced. Another way of looking at this is using the typical energy diagrams for chemical reactions. "

# ╔═╡ 2495b700-916f-4103-a822-1c085c357153
md"

Of course, in a reversible reaction, if the forward reaction is exothermic, the backward reaction is endothermic, and vice-versa. And the energy of activation is smaller for the exothermic than for the endothermic direction.  Lets see teh consequences of this difference.  To do that, lets play with the Arrienous equation, to get a straigh line

```math
\begin{align}  
\ln k &= \ln \left(Ae^{-E_a/RT} \right) \\
&= \ln A + \ln \left(e^{-E_a/RT}\right) \\
&= \left(\dfrac{-E_a}{R}\right) \left(\dfrac{1}{T}\right) + \ln A \\
&=\ln A - \dfrac{E_{a}}{RT}
\end{align}
```


Now let's plot the two kinetic rate constants for a reversible reaction. 
 "

# ╔═╡ c2cc318b-933a-43d4-be76-59a629d900c3
begin
	p5=plot(1 ./ T,log(A).-(E_a./2 ./ (0.082.*T)),ylims = (-2,0),label = "k1, Exothermic")
	p5=plot!(1 ./ T,log(A).-(E_a./(0.082.*T)),ylims = (-2,0),label = "k2, Endothermic")
	title!("Rate constant change with T in reversible reactions")
	ylabel!("log (k)")
	xlabel!("1/Temperature [1/K]")
end

# ╔═╡ d30be1ba-5658-4307-a69d-c6c6c97365b0
md"We can see that, both k1 and k2 increase with the temperature, but the dependence is stronger for the endothermic. Therefore, increasing the temperature  will result in an shift of the equilibirum towards the endothermic reaction (the reaction that consumes heat). On other words, more heat is available, the reaction that uses heat is favored). On the contrary a decrease in temperature, favors the reaction that increases the temperature.

As a conclusion: since the temperature affects differently both $ks$ of a reversible reaction, and their rate is the equilibrium constant, a change in the temperature will allways affect the equilibrium of a reversible reaction. 

For instance, lets look at the following reaction: 


```math
H_2+ I_2 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} 2IH \tag{24}
```

The equilibrium constant is $K_{eq} = \frac{[IH]^2}{[H_2][I_2]}$. If [$HI$]=0.75 M and [$H2$]= 0.20 M at equilibrium, then the concentration of [$I_2$] at equilibrium depends on the value of $K_{eq}$

```math
[I_2] = \frac{[0.75]^2}{[0.2][K_{eq}]} \tag{25}
```

"

# ╔═╡ 4409d520-1aae-4036-892c-27f2bfaec571
begin
	
	K_eq= LinRange(0.1,1,100)
	HI=0.75
	H2=0.2
	plot(K_eq,(HI^2)./(K_eq.*H2))
	title!("Equilibrium concentration of the reactant ")
	ylabel!("I_2 [M]")
	xlabel!("Equilibrium constant")
end

# ╔═╡ cb5bbeb3-1c07-4a20-ba2c-d4ea706609e0
md"As as general conclusion, __changes in the temperature of the system that result in an increase in $K_{eq}$ shift the equilibrium towards the products__, while a reduction in $K_{eq}$ shifts the concentrations at equilibrium favoring the reactants.  
"

# ╔═╡ 6aa8f4b5-c26e-476e-b00c-0b0b040bbd18
md" ### 1.4 Aplication to systems of interactions

Next, we can use these concepts from chemical systems to study any type of system. For instance, let's set a system of interactions where we can study the balance between single people and the formation of marriage couples


```math
2  a \overset{k_1}{\underset{k_2}{\longleftrightarrow}} b \tag{24}
```
In this very simple analogy, $k_1$ is the rate of marriages, $k_2$ is the rate of divorces. At any point, mariaages and divorces are taking place, but at a population level, the equilibrium can be reached. 

"

# ╔═╡ d385fd9c-b3a5-41a8-b03c-5c780f54d6d8
simpleODE4! = @ode_def abetterway4 begin
  da = - 2* k_1 * a^2 + 2 * k_2 * b 
  db =  k_1 * a^2- k_2 * b 
end k_1 k_2

# ╔═╡ 65ed9079-07bf-4722-bf4d-bd872d4de4ac
begin
	dimerization_slide = @bind 👍 html"<input type=range min=0.3e-1 max=5.3e-1 step=0.01>"
	release_slide = @bind 👎 html"<input type=range min=2.3e-1 max=8.3e-1 step=0.01>"
	
	md"""
	**Set the rates of marriage and divorce?**
	
	Rate of marriage: $(dimerization_slide)
	
	Rate of divorce: $(release_slide)
	"""
end

# ╔═╡ 86576d8a-8726-4125-b533-48b5d15ba021
begin
			p_ = (👍,👎)
			a_₀=0.5; # units (M)
			#b₀=0.05; # units (M)
			b_₀=0.00; # units (M)
			u₀2=[a_₀,b_₀];
	        prob4 = ODEProblem(simpleODE4!,u₀2,tspan,p_)
	        sol4 = solve(prob4);
	plot(sol4,xlabel="Time [s]",ylabel="Concentration [M]",title="2a <-> b");
end

# ╔═╡ e8b5f56e-7f44-46d6-9962-cec32eb237be
md"At equilibrium we can see the proportions by using the value of $K_{eq}$
```math
\begin{align*}
K_{eq}=\frac{k_1}{k_2} = \frac{[b_{eq}]}{[a_{eq}]^2}
\end{align*}
```
The value of $K_{eq}$= $(👍/👎)
"

# ╔═╡ 33730b5a-eab3-4d88-a96e-84bf7fa510ef
md"Now we use the Arrhenius Equation to see what sets the rates of marriage and divorce. :
```math
k_1 = A e^{-E_a/RT} 
```
In this analogy, we can see that, to perform a marriage, the activation energy can be identified as some sort of cost, perharps the cost of the wedding. While for the backwards reaction can be identified as the cost of a divorce. 

```math
k_2 = A e^{-E_a/RT} 
```
What is the temperature? It can be assumed as the wealth of the society. More wealth, more movement (easier to cross the energy barriers). 

Question, is this an exothermic or endothermic reaction ? Marriage consumes heat or releases heat. Depends on what is more expensive, to get married or to get divorce. Anoter way to see that is to explore what happens if the society gains wealth. If the amount of money avaliable increases, usually, the rate of divorce in a socienty increases.  

So, we are in the same situation as above, the rate of $k_2$ and $k_1$ both increase, but $k_2$ increases more rapidly. 
"


# ╔═╡ 0f5fd289-5f7b-4312-8610-337170f3e09c
plot(p5)

# ╔═╡ b7d6ee7d-8049-4e22-9ecf-307b939f040f
md"In conclusion, based on this, divorce is more expensive than marriage, becasue increasing wealth moves the equilibrium towards more single people. We can argue if this is true or not true, but it ios clear that probably teh model is too simple and we need to add more interactions (such as religious believes, divorce laws... )"

# ╔═╡ 4e711021-eddf-4cb2-a420-e6f17dad5542
md"Then the equilibrium constant can be calculated as:"

# ╔═╡ 54dd5337-c76f-4a23-8bec-a0cd57b433dd
exothermic_url = "https://cdn.kastatic.org/ka-perseus-images/dd9737ed130c0a9965efa9476715cd9084bc5a1d.svg"

# ╔═╡ e58e3af7-6b70-4e0b-964e-a670ac6923a2
md"""Another interesting graph is the progress in reaction energy for __exotermic reactions__ (reactions that release heat, as the reaction takes place). Despite being energetically favorable, we also have a Activation energy, $(Resource(exothermic_url))""" 

# ╔═╡ 5220090a-c8ef-4b9b-ae17-6692fc4b19e3
endothermic_url = "https://cdn.kastatic.org/ka-perseus-images/fad604021c159260b16946b55b2b3ae106c7f05f.svg"

# ╔═╡ 24e616ee-3312-411e-a214-a0ec39f7da9d
md"""  We can directly see in this type of energy diagrams that the $E_a$ represents the energy barrier of the transition state. In this case, the plot  represents an __endothermic reaction__ (energy of products is higher than energy of reactans, and a supply of energy in the form of temperature is required, which usually means that will be taken from the surroundings as the reaction takes place). 

In this type of reactions, if $E_a$ increases, more energy will be required to transit from reactants to products. These type of plots will be important later in the course, in the context of biochemical reactions (reactions that involve proteins and other biological molecules). $(Resource(endothermic_url))"""

# ╔═╡ b06bb803-6f5c-45ed-ba6a-faa88bc6e9b5
begin
	N₀=200
	t=collect(0:0.1:100)
	T=10
	μ=log(2)/T
	K=400
	plot(t,t->N₀*K/(N₀ +(K-N₀)*exp(-μ*t)),label="T=10",seriestype=:line,ylims = (200,400))
	T=5
	μ=log(2)/T
	plot!(t,t->N₀*K/(N₀ +(K-N₀)*exp(-μ*t)),label="T=5",seriestype=:line,ylims = (200,400))
	T=15
	μ=log(2)/T
	plot!(t,t->N₀*K/(N₀ +(K-N₀)*exp(-μ*t)),label="T=15",seriestype=:line,ylims = (200,400))
end

# ╔═╡ 1c24a2f1-2b15-48e6-b0d8-c3278103036d
begin
	N₀=200
	t=collect(0:0.1:10)
	μ1=log(2)/10
	μ2=log(2)/5
	μ3=log(2)/15
	plot(t,t->N₀*exp(μ1*t),label="T=10",seriestype=:line,ylims = (200,400))
	plot!(t,t->N₀*exp(μ2*t),label="T=5",seriestype=:line,ylims = (200,400))
	plot!(t,t->N₀*exp(μ3*t),label="T=15",seriestype=:line,ylims = (200,400))
end

# ╔═╡ c589178d-8b4b-488c-8d63-c7cc487848ec
begin
	A=1
	T= LinRange(100,300,100)
	plot(T,A.*exp.(-E_a./(0.082.*T)),ylims = (0,1))
	title!("Rate constant change with temperature")
	ylabel!("k")
	xlabel!("Temperature [K]")
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DifferentialEquations = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
ParameterizedFunctions = "65888b18-ceab-5e60-b2b9-181511a3b968"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
DifferentialEquations = "~7.1.0"
HypertextLiteral = "~0.9.3"
ParameterizedFunctions = "~5.13.1"
Plots = "~1.25.11"
PlutoUI = "~0.7.35"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "1ee88c4c76caa995a885dc2f22a5d548dfbbc0ba"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.2"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "56c347caf09ad8acb3e261fe75f8e09652b7b05b"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "0.7.10"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "ce68f8c2162062733f9b4c9e3700d5efc4a8ec47"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "0.16.11"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bijections]]
git-tree-sha1 = "705e7822597b432ebe152baa844b49f8026df090"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.3"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "5e98d6a6aa92e5758c4d58501b7bf23732699fa3"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.2"

[[deps.BoundaryValueDiffEq]]
deps = ["BandedMatrices", "DiffEqBase", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "NLsolve", "Reexport", "SparseArrays"]
git-tree-sha1 = "fe34902ac0c3a35d016617ab7032742865756d7d"
uuid = "764a87c0-6b3e-53db-9096-fe964310641d"
version = "2.7.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CPUSummary]]
deps = ["Hwloc", "IfElse", "Static"]
git-tree-sha1 = "849799453de85b55e78550fc7b0c8f442eb497ab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.8"

[[deps.CSTParser]]
deps = ["Tokenize"]
git-tree-sha1 = "6cc1759204bed5a4e2a5c2f00901fd5d90bc7a62"
uuid = "00ebfdb7-1f24-5e51-bd34-a7502290713f"
version = "3.3.1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "03dc838350fbd448fca0b99285ed4d60fc229b72"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "4cd7063c9bdebdbd55ede1af70f3c2f48fab4215"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.6"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DEDataArrays]]
deps = ["ArrayInterface", "DocStringExtensions", "LinearAlgebra", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "31186e61936fbbccb41d809ad4338c9f7addf7ae"
uuid = "754358af-613d-5f8d-9788-280bf1605d4c"
version = "0.2.0"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayDiffEq]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "LinearAlgebra", "Logging", "NonlinearSolve", "OrdinaryDiffEq", "Printf", "RecursiveArrayTools", "Reexport", "UnPack"]
git-tree-sha1 = "ceb3463f2913eec2f0af5f0d8e1386fb546fdd32"
uuid = "bcd4f6db-9728-5f36-b5f7-82caef46ccdb"
version = "5.34.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "ChainRulesCore", "DEDataArrays", "DataStructures", "Distributions", "DocStringExtensions", "FastBroadcast", "ForwardDiff", "FunctionWrappers", "IterativeSolvers", "LabelledArrays", "LinearAlgebra", "Logging", "MuladdMacro", "NonlinearSolve", "Parameters", "PreallocationTools", "Printf", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "StaticArrays", "Statistics", "SuiteSparse", "ZygoteRules"]
git-tree-sha1 = "433291c9e63dcfc1a0e42c6aeb6bb5d3e5ab1789"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.81.4"

[[deps.DiffEqCallbacks]]
deps = ["DataStructures", "DiffEqBase", "ForwardDiff", "LinearAlgebra", "NLsolve", "OrdinaryDiffEq", "Parameters", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "e57ecaf9f7875714c164ccca3c802711589127cf"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "2.20.1"

[[deps.DiffEqJump]]
deps = ["ArrayInterface", "Compat", "DataStructures", "DiffEqBase", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "StaticArrays", "TreeViews", "UnPack"]
git-tree-sha1 = "e30f058eb600407e3fd4ea082e2527e3a3671238"
uuid = "c894b116-72e5-5b58-be3c-e6d8d4ac2b12"
version = "8.2.1"

[[deps.DiffEqNoiseProcess]]
deps = ["DiffEqBase", "Distributions", "LinearAlgebra", "Optim", "PoissonRandom", "QuadGK", "Random", "Random123", "RandomNumbers", "RecipesBase", "RecursiveArrayTools", "Requires", "ResettableStacks", "SciMLBase", "StaticArrays", "Statistics"]
git-tree-sha1 = "d6839a44a268c69ef0ed927b22a6f43c8a4c2e73"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.9.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.DifferentialEquations]]
deps = ["BoundaryValueDiffEq", "DelayDiffEq", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffEqNoiseProcess", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEq", "Random", "RecursiveArrayTools", "Reexport", "SteadyStateDiffEq", "StochasticDiffEq", "Sundials"]
git-tree-sha1 = "3f3db9365fedd5fdbecebc3cce86dfdfe5c43c50"
uuid = "0c46a032-eb83-5123-abaf-570d42b7fbaa"
version = "7.1.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "5f5f0b750ac576bcf2ab1d7782959894b304923e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.9"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "7eb5d99577e478d23b1ba1faa9f8f6980d34d0a3"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.4"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "d7ab55febfd0907b285fbf8dc0c73c0825d9d6aa"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.3.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "3e1289d9a6a54791c1ee60da0850f4fd71188da6"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.11.0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FastBroadcast]]
deps = ["LinearAlgebra", "Polyester", "Static"]
git-tree-sha1 = "0f8ef5dcb040dbb9edd98b1763ac10882ee1ff03"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.1.12"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ec299fdc8f49ae450807b0cb1d161c6b76fd2b60"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.10.1"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "241552bc2209f0fa068b6415b1942cc0aa486bcc"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "3965a3216446a6b020f0d48f1ba94ef9ec01720d"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.6"

[[deps.Hwloc]]
deps = ["Hwloc_jll"]
git-tree-sha1 = "92d99146066c5c6888d5a3abc871e6a214388b91"
uuid = "0e44f5e4-bd66-52a0-8798-143a42290a1d"
version = "2.0.0"

[[deps.Hwloc_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d8bccde6fc8300703673ef9e1383b11403ac1313"
uuid = "e33a78d0-f292-5ffc-b300-72abe9b543c8"
version = "2.7.0+0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterativeSolvers]]
deps = ["LinearAlgebra", "Printf", "Random", "RecipesBase", "SparseArrays"]
git-tree-sha1 = "1169632f425f79429f245113b775a0e3d121457c"
uuid = "42fd0dbc-a981-5370-80f2-aaf504508153"
version = "0.9.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaFormatter]]
deps = ["CSTParser", "CommonMark", "DataStructures", "Pkg", "Tokenize"]
git-tree-sha1 = "fcfaddc61f766211b2c835d3eceaf999b6ea9555"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "0.22.4"

[[deps.KLU]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse_jll"]
git-tree-sha1 = "cae5e3dfd89b209e01bcd65b3a25e74462c67ee0"
uuid = "ef3ab10e-7fda-4108-b977-705223b18434"
version = "0.3.0"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "6333cc5b848295895f3b23eb763d020fc8e05867"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.7.12"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "0328ad9966ae29ccefb4e1b9bfd8c8867e4360df"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.3"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "LinearAlgebra", "MacroTools", "StaticArrays"]
git-tree-sha1 = "97e2adfcbe7ac07112ca79f03e34fc88cac6b9e7"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.7.2"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "6dd77ee76188b0365f7d882d674b95796076fa2c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "DocStringExtensions", "IterativeSolvers", "KLU", "Krylov", "KrylovKit", "LinearAlgebra", "RecursiveFactorization", "Reexport", "Requires", "SciMLBase", "Setfield", "SparseArrays", "SuiteSparse", "UnPack"]
git-tree-sha1 = "f27bb8e4eabdb93ed3703c55025b111e045ffe81"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "1.12.0"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "534aa24fae56f5f0956134d8789ab30d6fe2f615"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.102"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "0886d229caaa09e9f56bcf1991470bd49758a69f"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.3"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModelingToolkit]]
deps = ["AbstractTrees", "ArrayInterface", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqJump", "DiffRules", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "Graphs", "IfElse", "InteractiveUtils", "JuliaFormatter", "LabelledArrays", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "NaNMath", "NonlinearSolve", "RecursiveArrayTools", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SafeTestsets", "SciMLBase", "Serialization", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "Symbolics", "UnPack", "Unitful"]
git-tree-sha1 = "6d3dd18fbb1abf01894c5d064072285c6b863a98"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "8.5.1"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MuladdMacro]]
git-tree-sha1 = "c6190f9a7fc5d9d5915ab29f2134421b12d24a68"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.2"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "81b44a8cba10ff3cfb564da784bf92e5f834da0e"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.3"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.NonlinearSolve]]
deps = ["ArrayInterface", "FiniteDiff", "ForwardDiff", "IterativeSolvers", "LinearAlgebra", "RecursiveArrayTools", "RecursiveFactorization", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "UnPack"]
git-tree-sha1 = "b61c51cd5b9d8b197dfcbbf2077a0a4e1505278d"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "0.3.14"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "045d10789f5daff18deb454d5923c6996017c2f3"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.6.1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.OrdinaryDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "ExponentialUtilities", "FastClosures", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "Logging", "LoopVectorization", "MacroTools", "MuladdMacro", "NLsolve", "NonlinearSolve", "Polyester", "PreallocationTools", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "df82fa0f9f90f669cc3cf9e3f0400e431e0704ac"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.6.6"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7e2166042d1698b6072352c74cfd1fca2a968253"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.6"

[[deps.ParameterizedFunctions]]
deps = ["DataStructures", "DiffEqBase", "DocStringExtensions", "Latexify", "LinearAlgebra", "ModelingToolkit", "Reexport", "SciMLBase"]
git-tree-sha1 = "2f48f745e976dc5575bbc301e6c63b8fb5f12155"
uuid = "65888b18-ceab-5e60-b2b9-181511a3b968"
version = "5.13.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5c907bdee5966a9adb8a106807b7c387e51e4d6c"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.11"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.PoissonRandom]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "44d018211a56626288b5d3f8c6497d28c26dc850"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.0"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Requires", "Static", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "2232d3865bc9a098e664f69cbe340b960d48217f"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.6.6"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "dc11fa882240c43a875b48e21e6423704927d12f"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.4"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "LabelledArrays"]
git-tree-sha1 = "e4cb8d4a2edf9b3804c1fb2c2de57d634ff3f36e"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.2.3"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "ChainRulesCore", "DocStringExtensions", "FillArrays", "LinearAlgebra", "RecipesBase", "Requires", "StaticArrays", "Statistics", "ZygoteRules"]
git-tree-sha1 = "736699f42935a2b19b37a6c790e2355ca52a12ee"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.24.2"

[[deps.RecursiveFactorization]]
deps = ["LinearAlgebra", "LoopVectorization", "Polyester", "StrideArraysCore", "TriangularSolve"]
git-tree-sha1 = "7ad4c2ef15b7aecd767b3921c0d255d39b3603ea"
uuid = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
version = "0.2.9"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "256eeeec186fa7f26f2801732774ccf277f05db9"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.1.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "62c2da6eb66de8bb88081d20528647140d4daa0e"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "61a96d8b89083a53fb2b745f3b59a05359651bbe"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.30"

[[deps.SafeTestsets]]
deps = ["Test"]
git-tree-sha1 = "36ebc5622c82eb9324005cc75e7e2cc51181d181"
uuid = "1bc83da4-3b8d-516f-aca4-4fe02f6d838f"
version = "0.0.1"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "RecipesBase", "RecursiveArrayTools", "StaticArrays", "Statistics", "Tables", "TreeViews"]
git-tree-sha1 = "8ff1bf96965b3878ca5d235752ff1daf519e7a26"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.26.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SparseDiffTools]]
deps = ["Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays", "VertexSafeGraphs"]
git-tree-sha1 = "87efd1676d87706f4079e8e717a7a5f02b6ea1ad"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "1.20.2"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "6354dfaf95d398a1a70e0b28238321d5d17b2530"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.SteadyStateDiffEq]]
deps = ["DiffEqBase", "DiffEqCallbacks", "LinearAlgebra", "NLsolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "3e057e1f9f12d18cac32011aed9e61eef6c1c0ce"
uuid = "9672c7b4-1e72-59bd-8a11-6ac3964bc41f"
version = "1.6.6"

[[deps.StochasticDiffEq]]
deps = ["Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqJump", "DiffEqNoiseProcess", "DocStringExtensions", "FillArrays", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "Logging", "MuladdMacro", "NLsolve", "OrdinaryDiffEq", "Random", "RandomNumbers", "RecursiveArrayTools", "Reexport", "SparseArrays", "SparseDiffTools", "StaticArrays", "UnPack"]
git-tree-sha1 = "5f88440e7470baad99f559eed674a46d2b6b96f7"
uuid = "789caeaf-c7a9-5a7d-9973-96adeb23e2a0"
version = "6.44.0"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "ManualMemory", "Requires", "SIMDTypes", "Static", "ThreadingUtilities"]
git-tree-sha1 = "e0a02838565c4600ecd1d8874db8cfe263aaa6c7"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.2.12"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"

[[deps.Sundials]]
deps = ["CEnum", "DataStructures", "DiffEqBase", "Libdl", "LinearAlgebra", "Logging", "Reexport", "SparseArrays", "Sundials_jll"]
git-tree-sha1 = "76d881c22a2f3f879ad74b5a9018c609969149ab"
uuid = "c3572dad-4567-51f8-b174-8c6c989267f4"
version = "4.9.2"

[[deps.Sundials_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg", "SuiteSparse_jll"]
git-tree-sha1 = "04777432d74ec5bc91ca047c9e0e0fd7f81acdb6"
uuid = "fb77eaff-e24c-56d4-86b1-d163f2edb164"
version = "5.2.1+0"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "bfa211c9543f8c062143f2a48e5bcbb226fd790b"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.7"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "074e08aea1c745664da5c4b266f50b840e528b1c"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.3.0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "6dad289fe5fc1d8e907fa855135f85fb03c8fa7a"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.9"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "97e999be94a7147d0609d0b9fc9feca4bf24d76b"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.15"

[[deps.Tokenize]]
git-tree-sha1 = "0952c9cee34988092d73a5708780b3917166a0dd"
uuid = "0796e94c-ce3b-5d07-9a54-7f471281c624"
version = "0.5.21"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "1cda71cc967e3ef78aa2593319f6c7379376f752"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.72"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.TriangularSolve]]
deps = ["CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "LoopVectorization", "Polyester", "Static", "VectorizationBase"]
git-tree-sha1 = "5cbc1a4551fcf8afe8f80bb4f1f13e3271ee2656"
uuid = "d5829a12-d9aa-46ab-831f-fb7c9ab06edf"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "Hwloc", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "e9a35d501b24c127af57ca5228bcfb806eda7507"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.24"

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─081daad5-b960-4c64-bdec-c0ecb0d6896b
# ╟─4e2cead3-1c4f-48b4-8f5a-f78c4efaee2c
# ╠═16c39710-8f90-45c8-983a-25438019d90c
# ╠═bd8ee393-6193-4662-b199-edbe339ffc31
# ╟─951f16b4-f2d4-4855-91b2-b2ab2cf7df43
# ╠═99d88df8-e6c8-4fe7-a9ab-61f908912e92
# ╠═04bb79a4-5a72-4439-a3cd-b2b5ff985300
# ╠═06cbef49-b654-41dd-ac10-13474ef316a7
# ╠═9a68d323-e13d-42bb-970a-0ab5630155d5
# ╠═13aae8fd-f633-4073-a5ab-1e89f99e8568
# ╟─7b74322b-8481-4311-a77c-1a62cfb5b15c
# ╠═a1ce2262-da8e-496a-81d8-10ff5246c17f
# ╠═85defe33-1c87-4e81-a7d2-a363bee3e699
# ╠═b043d65b-a214-482a-96cb-e8c075814490
# ╠═842e2c41-72a3-443c-82db-ab1ff3612a12
# ╠═79fdee9d-9e12-4d2c-ad61-bd2b4c6dbb1a
# ╠═86804a6f-9178-40f3-ae82-cc6723412ae8
# ╠═8510e3df-ea93-4c25-ac3d-1069b067a62d
# ╠═753586ff-1721-48e3-befd-bce8f72ae285
# ╠═62542751-f6ae-4bf6-9a20-c362011f68d5
# ╟─831a828f-166e-4467-a94e-3cc7981a3e16
# ╠═1c127c36-c5fa-4a1f-beb1-8f3774abffaf
# ╠═3abe473c-67c7-4897-9345-1c39b5a2077c
# ╠═1dc48319-02d1-44d8-9873-2ebd28b358a8
# ╠═8e2e2fec-f264-45df-95fa-556ac631a8f2
# ╠═7ba042d6-107d-4188-a877-c7f7d7b35c46
# ╠═1c24a2f1-2b15-48e6-b0d8-c3278103036d
# ╟─85bc871e-6067-4a1d-bbb9-e3a9b0dfbd85
# ╠═b06bb803-6f5c-45ed-ba6a-faa88bc6e9b5
# ╠═61fe29b8-065c-430f-a28d-b463780f8b00
# ╠═a0bdad86-6cbd-4da1-b1ad-3c1f5e3c1ef3
# ╠═7f30ad4d-9964-4140-8003-051653d5f1e4
# ╠═435e892c-e56d-4f20-be01-cbf462194882
# ╠═ae184b0d-f61d-4aea-853b-17cb597d1087
# ╠═aecc1472-9167-4b09-a08c-9e701def7d54
# ╠═8e867da7-e8bc-4dca-8fe4-644f5277cb67
# ╠═d46dfe3b-f178-4d18-b95f-c961bc207466
# ╠═ca6a3577-ac75-4993-b2e5-4975ef4aaf9f
# ╠═349b1edf-3956-4b47-9d2f-760df409d40e
# ╠═3e8066bd-d8dc-47d9-98e5-0b8b0759e179
# ╠═2e9621c6-50c4-4d30-84f5-f26ea707a808
# ╟─65983c74-85d8-11ec-0c78-c9585d7356d8
# ╠═547db9e8-8de9-408d-80d4-48947b4aa1c2
# ╟─845a5477-29aa-4ebd-89a0-53ddeed9343d
# ╠═d096a6be-65a1-428d-9bfb-da7fe89f4c19
# ╠═c4bb7c6a-1165-44d6-8c1e-ea6798297c72
# ╟─5d00ba37-89a4-44a1-b8dc-c47394890b1f
# ╟─d4e5e672-4759-47ed-9a69-4d92efd85b05
# ╟─6b25dc03-352b-4351-b387-2c3d0dd80106
# ╟─eb9eb1de-7b0f-4045-b438-ab7930243e5c
# ╟─07735905-642b-4141-87fe-c9bc20ae04c2
# ╟─76691838-7086-4ca0-91d1-b5ef2c3e3b24
# ╠═23c143d3-29a3-4824-877e-1f13d0818ab6
# ╠═71777cd5-eedb-4dbf-a3a0-84ebe31ef602
# ╟─aed09f8b-eb07-4a34-a0bf-dc9ac987e85f
# ╠═0c575ced-ef89-4bf8-8f8c-c784d7f2d33a
# ╟─7b95dbaf-dc9d-45eb-9a0a-f2f75daf8765
# ╠═f0413f68-1eb4-4a3d-bdcf-61e8aa96c0e7
# ╠═6e728cf8-72e7-4567-90e9-a97f5f6a1f28
# ╟─2e15f653-b17b-4724-8fa9-a10606093c5d
# ╟─cd434441-2358-4d50-8744-98af4fc99176
# ╠═7a252b5a-13f1-489c-8dfd-ef325adeee56
# ╟─0a62b4ce-861b-44b2-bb56-c6f4cfce5fef
# ╟─e48547d5-f638-4885-8d2d-a8f70ad67dd5
# ╟─ea66ecaf-96c5-4b54-92eb-55d0eda5a001
# ╟─442944b8-8f11-484e-988c-6b1bcbd13c8e
# ╟─10d1f8cc-1a72-420b-b5cc-0286105850b3
# ╟─6f8ff973-6845-4517-8bca-1bec97fa4edf
# ╠═313f484e-bea2-437d-ba91-d766d36dd248
# ╠═f8ebbcd7-61e6-46bd-91fd-546ba931f76c
# ╠═7c053891-043a-4c0d-9afd-d03012d102c8
# ╠═5959c966-ff6e-4172-9fad-c05258050ff4
# ╠═29001c0e-9f2e-4fbb-ba8f-3b34de1d2a4b
# ╠═10e57e8a-bc19-4abe-9b62-ba7b65c69c9d
# ╟─6de3e149-e3e4-4c24-975b-e1e2259392f3
# ╟─1525b27e-2613-4f27-9f13-171e24cd574e
# ╟─f55ef39d-192e-4e97-acf5-a5fad4c6bae8
# ╠═481d6298-106a-4f6a-b822-7fda0a24bbe0
# ╠═c589178d-8b4b-488c-8d63-c7cc487848ec
# ╟─59fd13d8-e3bd-4756-a77f-a571876670a3
# ╠═24e616ee-3312-411e-a214-a0ec39f7da9d
# ╟─e58e3af7-6b70-4e0b-964e-a670ac6923a2
# ╟─2495b700-916f-4103-a822-1c085c357153
# ╟─c2cc318b-933a-43d4-be76-59a629d900c3
# ╟─d30be1ba-5658-4307-a69d-c6c6c97365b0
# ╠═4409d520-1aae-4036-892c-27f2bfaec571
# ╟─cb5bbeb3-1c07-4a20-ba2c-d4ea706609e0
# ╟─6aa8f4b5-c26e-476e-b00c-0b0b040bbd18
# ╠═d385fd9c-b3a5-41a8-b03c-5c780f54d6d8
# ╟─65ed9079-07bf-4722-bf4d-bd872d4de4ac
# ╠═86576d8a-8726-4125-b533-48b5d15ba021
# ╟─e8b5f56e-7f44-46d6-9962-cec32eb237be
# ╟─33730b5a-eab3-4d88-a96e-84bf7fa510ef
# ╠═0f5fd289-5f7b-4312-8610-337170f3e09c
# ╟─b7d6ee7d-8049-4e22-9ecf-307b939f040f
# ╟─4e711021-eddf-4cb2-a420-e6f17dad5542
# ╠═6d9759c3-63d9-4166-8884-cfaa99ee33c4
# ╠═54dd5337-c76f-4a23-8bec-a0cd57b433dd
# ╟─5220090a-c8ef-4b9b-ae17-6692fc4b19e3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
