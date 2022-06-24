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

# ╔═╡ 58a925ef-23db-4a78-9333-3aec9e53821d
using PlutoUI, HypertextLiteral

# ╔═╡ 0e64a8f0-e392-469a-ba9a-9e4805966b06
using Plots

# ╔═╡ 1c8134d8-8b40-11ec-20e5-03d3186dfbb3
md"# 7. Solving Differential Equations of Systems of Interactions Analitically 

Differential equations are central in the mathematical treatment of the dynamics of biological systems. Now that we know how to write the set of differential equations that govern any system of interacting species (using the differential form of Mass Action Law and the Mass Conservation Law), we will proceed now to integrate them and obtain an analytical solution and extract relevant information.

##  1. Irreversible First Order Reaction 

 We will start with the simplest example of one species  $X_1$ transforming into another $X_2$ with a kinetic constant $k_1$. The scheme of interactions is simply: 

```math
X_1 \overset{k_1}{\longrightarrow} X_2  \tag{1}
```

with initial concentrations $X_1(0)$ and $X_2(0)$ and stoichometric matrices

```math
A=\begin{bmatrix}
 1& 0\end{bmatrix} ;
B=\begin{bmatrix}
0 & 1  
\end{bmatrix} ; (B-A)^T= \begin{bmatrix}
 - 1\\ 1 \end{bmatrix} \tag{2}
```
 
in this particular case

```math
K=\begin{pmatrix}
 k\\ 
\end{pmatrix}
```

and 

```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^0
\end{pmatrix}
```

so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
 - 1\\ 1 \end{bmatrix} \begin{bmatrix}
k \end{bmatrix}\begin{pmatrix}
X_1^1\cdot X_2^0
\end{pmatrix} = \begin{bmatrix}
 - 1\\ 1 \end{bmatrix} \begin{pmatrix}
k \cdot X_1
\end{pmatrix} = \begin{bmatrix}
 - k \cdot X_1\\ k \cdot X_1 \end{bmatrix} 
\end{align}
```



The set of ODEs that define the system is

```math
\begin{align*}
   \frac{\mathrm{d} X_1}{\mathrm{d} t} &=  - k_1 X_1 \tag{3}\\
   \frac{\mathrm{d} X_2}{\mathrm{d} t} &=   k_1 X_1  \tag{4}\\
\end{align*}
```

It is called a first order reaction because the rate is proportional to the first power of the concentration of one reactant. Taking eq. 2, we reorganize variables and integrate.

```math
\begin{align*}
   \frac{\mathrm{d} X_1}{X_1} &=  - k_1 \mathrm{d}t \tag{5}\\
   \int  \frac{\mathrm{d} X_1}{X_1} &=  \int  - k_1 \mathrm{d}t \tag{6}\\
    \log(X_1(t)) &=  - k_1 t + C  \tag{7}\\
    X_1 (t)&=  e^{- k_1 t + C}  \tag{8}\\
    X_1 (t)&=  e^{- k_1 t}e^{C}  \tag{9}\\
\end{align*}
```

To calculate the integration constant, we use the intial condition at $t=0$.

```math
\begin{align*}
    X_1 (0)&=  e^{- k_1 0} e^{C}  \tag{10}\\
    X_1 (0)&=   e^{C}  \tag{11}\\
\end{align*}
```
so the final differential equation for $X_1$ takes the form:

```math
X_1 (t)= X_1 (0) e^{- k_1 t}  \tag{12}
```



Let's use now the conservation of mass to obtain the equation for the reactant $X_2$. Based on this law, if $C \cdot (B-A)^T =0$, then $\sum_{i=1}^{r} C_i \cdot X_i = cte $, so for the case of two variables and one single irreversible reaction 


```math
\begin{align*}
C_1 \cdot X_1 + C_2 \cdot X_2 &= cte \tag{13}\\
\end{align*}
```

so
```math
\begin{align*}
C \cdot (B-A)^T &=0 \tag{14}\\ 
\begin{pmatrix}C_1 & C_2\end{pmatrix} \begin{bmatrix}
 - 1\\ 1 \end{bmatrix}&=0 \tag{15}\\
 -C_1 + C_2 &= 0 \tag{16}\\
 C_1 &= C_2 \tag{17}
\end{align*}
```
Therefore,

```math
\begin{align*}
C_1 \cdot X_1(t) + C_1 \cdot X_2(t) &= cte \tag{18}\\
\end{align*}
```

for any value of $C_1$, so taking $C_1=1$, we have 

```math
X_1(t)+X_2(t)=cte \tag{19}
```

This constant value has to be true also at $t=0$,

```math
X_1(t)+X_2(t)= X_1(0)+X_2(0) = cte \tag{20}
```

using this conservation of mass, we can also obtain the analytical solution for $X_2$:

```math
\begin{align*}
X_2 (t) &= X_1 (0) + X_2 (0) - X_1 (t)\tag{21} \\
X_2 (t) &= X_1 (0) + X_2 (0) - X_1 (0) e^{- k_1 t}  \tag{22}\\
X_2 (t) &=  X_2 (0) + X_1 (0) (1-e^{- k_1 t})  \tag{23}
\end{align*}
```
"

# ╔═╡ 0a299414-2f9f-4182-ac44-1154230a67bc
@htl("""

<div class='blue-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Computer Exercise"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

""")

# ╔═╡ 5a28c78e-4d5e-402f-bddc-f5470f2c13af
md"Write a computer program that consists on a function that takes as input the time of the reaction, a vector of initial concentrations of the species $X_i$ and the reaction rate $k$ and plots the time evolution of the variables $X_i(t)$ for an irreversible first order reaction."

# ╔═╡ 9ca66203-c94b-4a51-80e7-a2c629cd0d75
gr()

# ╔═╡ e3ef157c-98d6-4579-9891-32bf00fb3b13
begin
	X1 =@bind 🐶 Slider(0.1:0.1:2.0)
	X2 =@bind 🐱 Slider(0.1:0.1:2.0)
	
	md"""
	**Use the sliders to set the initial conditions**
	
	X_{1}: $(X1)
	
	X_{2}: $(X2)
	"""
end

# ╔═╡ e0f47a77-3e55-40aa-b0a1-e1c41e02a041
md"""
X_{1}(0)= $(🐶), X_{2}(0)= $(🐱)
"""

# ╔═╡ 977a9c92-7777-41ef-b4aa-67b99fc640d5
begin
	t=collect(0:0.1:10)
	X₀ = [🐶,🐱]
	k_1 = 1.;
	plot(t,t->X₀[1]*exp(-k_1*t),label="\\ X_1",seriestype=:line)
	plot!(t,t->X₀[2]+X₀[1]*(1-exp(-k_1*t)),label="\\ X_2",seriestype=:line)
	title!("First order reaction")
	xaxis!("Time")
	yaxis!("Concentration")
end

# ╔═╡ 447fcd26-8908-4563-aa36-51cf69885bb1
md"## 2. Reversible First Order Reaction

Another simple example that can be solved analytically is the reversible first order reaction :


```math
X_1 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} X_2  \tag{24}
```
```math
A=\begin{bmatrix}
 1 & 0\\
 0 & 1 \end{bmatrix} ;
B=\begin{bmatrix}
0 & 1\\
 1 & 0  
\end{bmatrix} ; (B-A)^T= \begin{bmatrix}
-1 & 1\\
 1 & -1   \end{bmatrix} \tag{25}
```

 in this particular case

```math
K=\begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix} \tag{26}
```


```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^0\\
X_1^0\cdot X_2^1
\end{pmatrix} = \begin{pmatrix}
X_1\\
 X_2
\end{pmatrix}  \tag{27}
```

so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
-1 & 1\\
 1 & -1   \end{bmatrix} \begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix}\begin{pmatrix}
X_1\\
 X_2
\end{pmatrix}  
\end{align}
```

```math
\begin{align}
\begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}=  \begin{bmatrix}
-1 & 1\\
 1 & -1   \end{bmatrix}\begin{pmatrix}
k_1 \cdot X_1\\
  k_2 \cdot X_2
\end{pmatrix} =   \begin{bmatrix}
- k_1 \cdot X_1 +  k_2 \cdot X_2\\
 k_1 \cdot X_1 -k_2 \cdot X_2   \end{bmatrix}
\end{align} \tag{28}
```

The differential equations of the system can be written using the differential form of Mass Action Law as: 

```math
\begin{align}
\frac{dX_1}{dt} = -k_1 \cdot X_1 +k_2 \cdot X_2 \tag{29}\\
\frac{dX_2}{dt} = k_1 \cdot X_1 - k_2 \cdot X_2 \tag{30}\\
\end{align}
```

The mass conservation law is equivalent to the previous case (only $X_1$ and $X_2$ as reactants). We start by using mass conservation to de-couple the two equations. For the first one we have:

```math
\begin{align}
\frac{dX_1}{dt} = k_2\left( X_1(0)+X_2(0)-X_1\right)-k_1X_1 \tag{31}\\
\end{align}
```

which, when we rearrange terms becomes:
```math
\frac{dX_1}{dt}+\left(k_2+k_1\right)X_1=k_2\left( X_1(0)+X_2(0)\right) \tag{32}
```

To solve this diferential equation, we need to calculate the integration factor, $e^{ \int p(x)dx }$, that allows us to put the first term of the equation as one diferentiation. In this case $p(x)=k_2+k1 $, so $e^{ \int (k_2+k1)dt }=e^{(k_2+k1)t}$. This way, when we multiply both terms in eq. 18 by the integrating factor, we obtain:

```math
\frac{dX_1}{dt} e^{(k_2+k1)t} +\left(k_2+k_1\right) e^{(k_2+k1)t} X_1 = \left( X_1(0)+X_2(0)\right) k_2 e^{(k_2+k1)t} \tag{33}
```

The first term now is the derivative of $X_1 \cdot  e^{(k_2+k1)t}$ so we can reqwrite:

```math
\frac{d(X_1 e^{(k_2+k1)t})}{dt} = \left( X_1(0)+X_2(0)\right) k_2 e^{(k_2+k1)t} \tag{34}
```

we move the `dt` to the right side we integrate 

```math
\int d(X_1 e^{(k_2+k1)t}) = \left( X_1(0)+X_2(0)\right) k_2 \int e^{(k_2+k1)t} dt \tag{35}
```

solving

```math
X_1 e^{(k_2+k1)t} = \frac{\left( X_1(0)+X_2(0)\right) k_2}{k_2+k1}  e^{(k_2+k1)t} + C \tag{36}
```

rearranging terms becomes

```math
X_1 (t) = \frac{\left( X_1(0)+X_2(0)\right) k_2}{k_2+k1} + C \cdot e^{-(k_2+k1)t} \tag{37}
```

to calculate the integration constant we use initial value of $X_1$

```math
\begin{align}
X_1 (0) &= \frac{\left( X_1(0)+X_2(0)\right) k_2}{k_2+k1} + C \cdot e^{-(k_2+k1)0} \tag{38}\\
X_1 (0) &=\frac{\left( X_1(0)+X_2(0)\right) k_2}{k_2+k1} + C \tag{39}\\
C&= X_1 (0) - \frac{\left( X_1(0)+X_2(0)\right) k_2}{k_2+k1} \tag{40}\\
C&= \frac{X_1 (0)(k_2 + k 1) - X_1(0) k_2 - X_2(0) k_2 }{k_2+k1}\tag{41} \\
C&= \frac{X_1 (0) k 1 - X_2(0) k_2 }{k_2+k1}\tag{42} \\
\end{align}
```


Therefore the full solutions for $X_1$ is:
```math
X_1(t)=\frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2}+\frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{43}
```
"

# ╔═╡ aa4d3b18-62eb-4524-a736-9942e93ccbea

md"
In this case, the conservation of mass is exactly the same as the previous case of irreversible reaction, so:
```math
X_1(t)+X_2(t)=X_1(0)+X_2(0) \tag{44}
```


and, using the conservation of mass, the solution for $X_2$ is

```math
\begin{align}
X_2(t)&=X_1(0)+X_2(0)-\frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2}-\\
& - \frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{45}\\
X_2(t)&=\frac{(k_1+k_2)(X_1(0)+X_2(0))} {k_1+k_2}- \\ 
& - \frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2}-\frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{46}\\
X_2(t)&=\frac{k_1 (X_1(0)+X_2(0))} {k_1+k_2}+ \frac{k_2 (X_1(0)+X_2(0))}{k_1+k_2} - \\ 
& -\frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2}-\frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{47}
\end{align}
```

so the final set of equations is:
```math
\begin{align}
X_1(t)&=\frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2}+\frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{48}\\
X_2(t)&=\frac{k_1 (X_1(0)+X_2(0))} {k_1+k_2}-\frac{k_1 X_1(0)-k_2X_2(0) }{k_1+k_2}e^{-\left (k_1+k_2)\right )t} \tag{49}
\end{align}
```

where the equilibrium concentrations for $X_1$ and $X_2$ can be calculated using limits: 

```math
\begin{align}X_1(t->inf)\frac{k_2\left ( X_1(0)+X_2(0) \right )}{k_1+k_2} \tag{50}\\
X_2(t->inf)\frac{k_1\left ( X_1(0)+X_2(0) \right )}{k_1+k_2} \tag{51}
\end{align}
```"

# ╔═╡ 6224d853-94c1-4f11-b349-d003fa47ee2f
@htl("""

<div class='blue-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Computer Exercise"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

""")

# ╔═╡ ed9da8e1-723a-4e9b-a4c3-8a0bb925e01c
md"Write a computer program that consists on a function that takes as input the time of the reaction, a vector of initial concentrations of the species $X_i$ and the reaction rate $k$ and plots the time evolution of the variables $X_i(t)$ for a reversible first order reaction."

# ╔═╡ c81fc401-75de-42ec-9d3f-f04faacaadee
begin
	X3 =@bind 🔴 Slider(0.1:0.1:2.0)
	X4 =@bind 🟡 Slider(0.1:0.1:2.0)
	
	md"""
	**Use the sliders to set the initial conditions**
	
	X_{3}: $(X3)
	
	X_{4}: $(X4)
	"""
end

# ╔═╡ b673041f-9b22-4ede-bf5e-82c7fed0dbf7
md"""
X_{3}(0)= $(🔴), X_{4}(0)= $(🟡)
"""

# ╔═╡ 0d351d9e-32ab-4667-9563-a85b65e6cf0a
begin
	t₀=0
	t_final=10
	k_2 = 0.5
	X₀[1] = 🔴 # Initial concentrations
	X₀[2] = 🟡 # Initial concentrations
plot(t,t->((k_2*(X₀[1]+X₀[2]))/(k_1+k_2)) + ((k_1*X₀[1]-k_2*X₀[2])/(k_1+k_2))*exp(-t*(k_1+k_2)),label="\\ X_1",seriestype=:line)
plot!(t,t->X₀[1] + X₀[2] - ((k_2*(X₀[1]+X₀[2]))/(k_1+k_2)) - ((k_1*X₀[1]-k_2*X₀[2])/(k_1+k_2))*exp(-t*(k_1+k_2)),label="\\ X_2",seriestype=:line)
title!("Reversible first order")
xaxis!("Time")
yaxis!("Concentration")
end

# ╔═╡ bea913a2-ccf1-44e5-b710-8d0c8568605f
md"## 3. Second order reaction: irreversible dimerization

A special case of second order reaction is when conversion into product occur in pairs. This corresponds to a situation known as dimerization, very common in biology. In the next section, we will derive the solution for this scheme of interaction, where two identical molecules interact to form a complex. The general scheme of interaction is the following:

```math
2X_1 \overset{k_1}{\longrightarrow} X_2  \tag{52}
```

with initial concentrations $X_1(0)$ and $X_2(0)$ and stoichometric matrices

```math
A=\begin{bmatrix}
 2& 0\end{bmatrix} ;
B=\begin{bmatrix}
0 & 1  
\end{bmatrix} ; (B-A)^T= \begin{bmatrix}
 - 2\\ 1 \end{bmatrix} \tag{53}
```
 
in this particular case

```math
K=\begin{pmatrix}
 k_1\\ 
\end{pmatrix} \tag{54}
```

and 

```math
X^A=\begin{pmatrix}
X_1^2\cdot X_2^0
\end{pmatrix} \tag{55}
```

so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
 - 2\\ 1 \end{bmatrix} \begin{bmatrix}
k_1 \end{bmatrix}\begin{pmatrix}
X_1^2\cdot X_2^0
\end{pmatrix} = \begin{bmatrix}
 - 2\\ 1 \end{bmatrix} \begin{pmatrix}
k_1 \cdot X_1^2
\end{pmatrix} = \begin{bmatrix}
 - 2 k_1 X_1^2\\ k_1  X_1^2 \end{bmatrix}  \tag{56}
\end{align}
```

```math
\begin{align*}
   \frac{\mathrm{d} X_1}{\mathrm{d} t} &=  - 2 k_1 X_1^2 \tag{57}\\
   \frac{\mathrm{d} X_2}{\mathrm{d} t} &=    k_1 X_1^2  \tag{58}\\
\end{align*}
```

It is called a second order reaction because the rate is proportional to the second power of the concentration of reactant. Taking the first differential equation, we reorganize variables and integrate.
```math
\begin{align*}
   \frac{\mathrm{d} X_1}{X_1^2} &=  - 2k_1 \mathrm{d}t \tag{59}\\
   \int  \frac{\mathrm{d} X_1}{X_1^2} &=  \int  - 2k_1 \mathrm{d}t \tag{60}\\
     \frac{-1}{X_1} &=  - 2k_1 \cdot t + C \tag{61}\\
     X_1(t) &=\frac{1}{2k_1 \cdot  t - C}  \tag{62}\\
\end{align*}
```

To calculate the integration constant, we use the intial condition at $t=0$.

```math
\begin{align*}
     X_1(0) &=\frac{1}{2k_1  \cdot 0 - C} = - \frac{1}{C} \tag{63}
\end{align*}
```


so the final equation of is
```math
\begin{align}
X_1(t) &=\frac{1}{2k_1 \cdot t + \frac{1}{X_1(0)} } \tag{64}
\end{align}
```
so, multiplying numerator and denominator of the right side of the equation by $X_1(0)$ we obtain

```math
\begin{align}
X_1(t) &=\frac{X_1(0)}{(X_1(0))(2 k_1  \cdot t  + \frac{1}{X_1(0)}) }= \frac{X_1(0)}{2k_1 X_1(0) \cdot t + 1}\tag{65}
\end{align}
```


Let's use now the conservation of mass to obtain the equation for the reactant $X_2$. Remember, if $C \cdot (B-A)^T =0$, then $\sum_{i=1}^{r} C_i \cdot X_i = cte $, so

```math
\begin{align*}
C_1 \cdot X_1 + C_2 \cdot X_2 &= cte \tag{66}\\
\end{align*}
```

so
```math
\begin{align*}
C \cdot (B-A)^T &=0 \tag{67}\\ 
\begin{pmatrix}C_1 & C_2\end{pmatrix} \begin{bmatrix}
 - 2\\ 1 \end{bmatrix}&=0\tag{68}\\
 -2 C_1 + C_2 &= 0 \tag{69}\\
 2 C_1 &= C_2 \tag{70}
\end{align*}
```

Therefore,

```math
\begin{align*}
 C_1 \cdot X_1(t) + 2 C_1 \cdot X_2(t) &= cte \tag{71}\\
\end{align*}
```

for any value of $C_1$, so taking $C_1=1$, we have 

```math
X_1(t)+ 2 X_2(t)=cte \tag{72}
```

This constant value has to be true also at $t=0$,

```math
X_1(t)+2 X_2(t)= X_1(0)+ 2 X_2(0) = cte \tag{77}
```


using this restriction, we obtain for $X_2$:
```math
\begin{align}
X_2(t) &=  \frac{X_1(0)}{2}+ X_2(0) - \frac{X_1(0)}{4k_1 X_1(0) \cdot t + 1}\tag{78}\\
X_2(t) &= X_2(0) +  \frac{X_1(0)}{2} (1 - \frac{1}{2k_1 X_1(0) \cdot t + 1})\tag{79}\\
X_2(t) &= X_2(0) +  \frac{X_1(0)}{2} \frac{2k_1 X_1(0)t}{2k_1 X_1(0) \cdot t + 1}\tag{80}\\
X_2(t) &= X_2(0) +  \frac{k_1 X_1(0)^2 t}{2k_1 X_1(0) \cdot t + 1}\tag{81}\\
\end{align}
```



"

# ╔═╡ 04549d00-0d7f-472b-abb8-1b64d80da9c2
@htl("""

<div class='blue-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Computer Exercise"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

""")

# ╔═╡ d75faa1b-35f6-4650-b291-b1b7f798f0ad
md"Write a computer program that consists on a function that takes as input the time of the reaction, a vector of initial concentrations of the species $X_i$ and the reaction rate $k$ and plots the time evolution of the variables $X_i(t)$ for an irreversible dimerization reaction.

Compare the dynamics of consumption of $X_1$ between this second order reaction and the firts order irrevensible reaction. Why are they so different, even when we use the exact same values?"

# ╔═╡ b16367a6-cf40-47ff-9fb3-27129cbc4959
begin
	X5 =@bind 🟥 Slider(0.1:0.1:2.0)
	X6 =@bind 🟨 Slider(0.1:0.1:2.0)
	
	md"""
	**Use the sliders to set the initial conditions**
	
	X_{5}: $(X5)
	
	X_{6}: $(X6)
	"""
end

# ╔═╡ 3c1af85d-e895-4cf5-b9e3-0fece84ccd49
md"""
X_{1}(0)= $(🟥), X_{2}(0)= $(🟨)
"""

# ╔═╡ 5dd86c83-e59d-4f31-bc26-1a9ed0ce531f
begin
	X₀[1] = 🟥 # Initial concentrations
	X₀[2] = 🟨 # Initial concentrations
	
P3=plot(t,t->X₀[1]*exp(-k_1*t),label="monomer reaction \\  X_1",seriestype=:line)
P3=plot!(t,t->X₀[2]+X₀[1]*(1-exp(-k_1*t)),label="monomer reaction \\ X_2",seriestype=:line)
title!("Monomer reaction")
P4=plot(t,t->X₀[1]/((2*k_1*t*X₀[1])+1),label="dimer reaction \\ X_1 ")
P4=plot!(t,t->X₀[2]+(k_1*t*X₀[1]^2)/((2*k_1*t*X₀[1])+1),label="dimer reaction \\ X_2 ")
title!("Dimer reaction")


plot(P3,P4,layout=(1,2),legend=true)

xaxis!("Time")
yaxis!("Concentration")
end

# ╔═╡ 7bb7da1e-39d5-4ba8-bbc4-19d8b2edb932
md"We cannot compare because the units of the kinetic constants are diffrenet between first order [1/time] and second order [1/(time*concentration)] reactions. To compare the dynamcis we should adjust these two values. One way of comparing them is to divide de kinetic constant of the second order dimer reaction by the initial concentration of $X_1$ "

# ╔═╡ 08fbbe98-e7eb-4d02-ad29-267453509c8c
begin
	P1=plot(t,t->X₀[1]*exp(-k_1*t),label="monomer reaction \\  X_1",seriestype=:line)
	P1=plot!(t,t->X₀[2]+X₀[1]*(1-exp(-k_1*t)),label="monomer reaction \\ X_2",seriestype=:line)
	title!("Monomer reaction")
#	k_1 = 1.0 / X₀[1]
	P2=plot(t,t->X₀[1]/((2*k_1/ X₀[1]*t*X₀[1])+1),label="dimer reaction \\ X_1 ")
	P2=plot!(t,t->X₀[2]+(k_1/ X₀[1]*t*X₀[1]^2)/((2*k_1/ X₀[1]*t*X₀[1])+1),label="dimer reaction \\ X_2 ")
	title!("Dimer reaction")
	
	
	plot(P1,P2,layout=(1,2),legend=true)
	
	xaxis!("Time")
	yaxis!("Concentration")
end

# ╔═╡ 76865e9e-074b-4c3e-a811-d02492d69366
md"## Reversible Dimerization

Another simple reaction scheme is the reversible dimerization (a second order reaction forward and a first order reaction backwards):

```math
2 X_1 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} X_2  \tag{82} 
```

```math
A=\begin{bmatrix}
 2 & 0\\
 0 & 1 \end{bmatrix} ;
B=\begin{bmatrix}
0 & 1\\
 2 & 0  
\end{bmatrix} ; (B-A)^T= \begin{bmatrix}
-2 & 2\\
 1 & 1   \end{bmatrix}  \tag{83} 
```
 in this particular case

```math
K=\begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix} \tag{83} 
```

```math
X^A=\begin{pmatrix}
X_1^2 X_2^0\\
X_1^0 X_2^1
\end{pmatrix} = \begin{pmatrix}
X_1^2\\
 X_2
\end{pmatrix} \tag{84} 
```

so 
"

# ╔═╡ 13947e3c-f393-40af-8f4e-dc62bb64a3cb
md"

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
-2 & 2\\
 1 & -1   \end{bmatrix} \begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix}\begin{pmatrix}
X_1^2\\
 X_2
\end{pmatrix}  =  \begin{bmatrix}
-2 & 2\\
 1 & -1   \end{bmatrix}\begin{pmatrix}
k_1  X_1^2\\
  k_2 X_2
\end{pmatrix} 
\end{align}
```


"

# ╔═╡ 5df18e02-dc27-4851-9954-392bb6770d88
md"
```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix} =   \begin{bmatrix}
- 2 k_1  X_1^2 +  2 k_2  X_2\\
 k_1  X_1^2 - k_2  X_2   \end{bmatrix}\end{align} \tag{85} 
```


"

# ╔═╡ 3f0505c7-db5a-497c-9408-55b776d3f646
md" The differential equations of the system can be written using the differential form of Mass Action Law as: 

```math 
\begin{align}
\frac{dX_1}{dt} = - 2k_1 \cdot X_1^2 + 2k_2 \cdot X_2 \tag{86} \\
\frac{dX_2}{dt} = k_1 \cdot X_1^2 - k_2 \cdot X_2 \tag{87} \\
\end{align}
```

The mass conservation law is equivalent to the previous case (only $X_1$ and $X_2$ as reactants). $X_2(t) =  \frac{X_1(0)}{2}+ X_2(0) - \frac{X_1(t)}{2}$ We start by using mass conservation to de-couple the two equations. For the first one we have:

```math 
\begin{align}
\frac{dX_1(t)}{dt} =  k_2 \left(X_1(0)+ 2 X_2(0)- X_1(t)\right)- 2k_1X_1(t)^2 \tag{88} \\
\end{align}
```
We can solve the integral but we cannot express the value of $X_1$ directly as a function of $t$, so we have to solve it using other methods. We will learn how to do this in the next classes. 
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
HypertextLiteral = "~0.9.3"
Plots = "~1.26.0"
PlutoUI = "~0.7.36"
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

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

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

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

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
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

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

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

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

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "3f7cb7157ef860c637f3f4929c8ed5d9716933c6"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.7"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

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
git-tree-sha1 = "23d109aad5d225e945c813c6ebef79104beda955"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.26.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "2c87c85e397b7ffed5ffec054f532d4edd05d901"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.36"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

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

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

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

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

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
# ╟─1c8134d8-8b40-11ec-20e5-03d3186dfbb3
# ╠═58a925ef-23db-4a78-9333-3aec9e53821d
# ╟─0a299414-2f9f-4182-ac44-1154230a67bc
# ╟─5a28c78e-4d5e-402f-bddc-f5470f2c13af
# ╠═0e64a8f0-e392-469a-ba9a-9e4805966b06
# ╠═9ca66203-c94b-4a51-80e7-a2c629cd0d75
# ╠═e3ef157c-98d6-4579-9891-32bf00fb3b13
# ╠═e0f47a77-3e55-40aa-b0a1-e1c41e02a041
# ╠═977a9c92-7777-41ef-b4aa-67b99fc640d5
# ╟─447fcd26-8908-4563-aa36-51cf69885bb1
# ╟─aa4d3b18-62eb-4524-a736-9942e93ccbea
# ╟─6224d853-94c1-4f11-b349-d003fa47ee2f
# ╟─ed9da8e1-723a-4e9b-a4c3-8a0bb925e01c
# ╟─c81fc401-75de-42ec-9d3f-f04faacaadee
# ╟─b673041f-9b22-4ede-bf5e-82c7fed0dbf7
# ╠═0d351d9e-32ab-4667-9563-a85b65e6cf0a
# ╟─bea913a2-ccf1-44e5-b710-8d0c8568605f
# ╟─04549d00-0d7f-472b-abb8-1b64d80da9c2
# ╟─d75faa1b-35f6-4650-b291-b1b7f798f0ad
# ╟─b16367a6-cf40-47ff-9fb3-27129cbc4959
# ╟─3c1af85d-e895-4cf5-b9e3-0fece84ccd49
# ╠═5dd86c83-e59d-4f31-bc26-1a9ed0ce531f
# ╟─7bb7da1e-39d5-4ba8-bbc4-19d8b2edb932
# ╠═08fbbe98-e7eb-4d02-ad29-267453509c8c
# ╟─76865e9e-074b-4c3e-a811-d02492d69366
# ╟─13947e3c-f393-40af-8f4e-dc62bb64a3cb
# ╟─5df18e02-dc27-4851-9954-392bb6770d88
# ╟─3f0505c7-db5a-497c-9408-55b776d3f646
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
