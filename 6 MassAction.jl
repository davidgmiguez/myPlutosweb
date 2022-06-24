### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ fffc6cac-c044-4cd9-a3b7-084341878faf
using HypertextLiteral

# ╔═╡ dd3e39ac-8b2c-11ec-1774-7543d2573dc0
md"# 5. The Mass Action Law

## 5.1 Definition of the Law of Mass Action 

Now that we know how to obtain the stoichiometric matrices and the order of any given reaction, the next step is to write the differential equations that capture the dynamics, i.e., the changes in time of the amount of all the species involved.  To illustrate this, we will use a generic scheme of interaction in vector notation:

$$\begin{align*}
X_1 + X_2 \overset{k}{\longrightarrow} X_3  \tag{58}\\ 
\end{align*}$$
The stoichometrix matrix is: "

# ╔═╡ 7c7a696a-8716-4739-9a40-44a34c621485
md" The matrix A is [1,1,0]"

# ╔═╡ 30a79e11-2165-4908-94a3-d020cc6d2e16
md"The order of the reaction is $(sum([1,1,0]))"

# ╔═╡ e11754b1-6993-403b-a0ec-e0aaf86a5dd3
md"The Units of k are s^-1 x k^$(1-(sum([1,1,0])))"

# ╔═╡ ab0aa754-78e2-47b5-a398-a7e1ab8d6c58
md"We will write down the differential equations that account for the dynamics of each species. 

```math
\begin{align}
\frac{\mathrm{d} X_1}{\mathrm{d} t}  &=&  \tag{59}\\
\frac{\mathrm{d} X_2}{\mathrm{d} t} &=& \tag{60}\\
\frac{\mathrm{d} X_3}{\mathrm{d} t} &=&  \tag{61}
\end{align}
```

The __Mass Action Law__ describes how a chemical reaction takes place under ideal and constant conditions. In brief, it is based on the fact that the speed of any reaction is proportional to the concentrations of the reactants involved, which has been discussed previously in this course (the proportionalinty constant is the kinetic rate constant). From here, it follows that the rate of change of any reactant is also proportional to the concetrations of the reactants involved. This velocity depends on the affinity of the reaction and the concentrations of each of the reactants (in any scheme of interactions, the ones at the left of the arrow). The affinity of the reactions is defined by the rate constant `k` of the reaction. The `k` is the proportionality constant that correlates the speed fo the reaction with the concentrations of the species involved, which in this case are $X_1$ and $X_2$:

```math
\begin{align}
\frac{\mathrm{d} X_1}{\mathrm{d} t}  &= - k X_1 X_2 \tag{62}\\
\frac{\mathrm{d} X_2}{\mathrm{d} t} &=- k X_1 X_2 \tag{63}\\
\frac{\mathrm{d} X_3}{\mathrm{d} t} &=  k X_1 X_2\tag{64} 
\end{align}
```

Lets see if the units are correct:

```math
\begin{align}
\frac{[M]}{[s]}  &= - \frac{1}{[s][M]}[M] [M] \tag{65}\\
\frac{[M]}{[s]}  &= - \frac{1}{[s][M]}[M] [M] \tag{66}\\
\frac{[M]}{[s]}  &=  \frac{1}{[s][M]}[M] [M] \tag{67}\\
\end{align}
```


For the case where more than one molecule of reactant is involved, we will use the reaction of formation of water. 
"

# ╔═╡ d2e26cb1-f136-4b96-8ef3-d3ee30f61b3d
md" The matrix A is [2 1 0]"

# ╔═╡ 10f98ec7-fcad-4587-90f8-62f7b1b74f4d
md" The matrix B is [0 0 2]"

# ╔═╡ fee636f1-4944-46a4-bcbe-32b7273268c7
md"The order of the reaction is $(sum([2 1 0]))"

# ╔═╡ 56a68cf6-6e64-4875-b614-18987857bc06
md"The Units of k are s^-1 x k^$(1-(sum([2 1 0])))"

# ╔═╡ 05faae5d-fdb6-4f1d-a77f-0a1337e1666d
StMat=([0 0 2]-[2 1 0])'

# ╔═╡ 0d524884-afbf-4594-80ee-9fdfcdabbea0
md"we can decompose the equation as

$$\begin{align*}
X_1 +  X_1 + X_2 \overset{k}{\longrightarrow}  X_3 + X_3 \tag{68}\\ 
\end{align*}$$

We will write down the differential equations that account for the dynamics of each species. The speed of each reaction depends on the concentration of the species involved (what is before the arrow). In this case is two units of $X_1$ and one unit of $X_2$. In addition, since the reaction consumes 2 units of $X_1$ for each unit of $X_2$, the rate of consumption of $X_1$ is twice the rate of consumption $X_2$. Another more automaticic way of doing it is to use the matrix of stochoimetric coefficients $(B-A)'$. Therefore:

```math
\begin{align}
\frac{\mathrm{d} X_1}{\mathrm{d} t}  &= - 2k X_1^2 X_2 \tag{69}\\
\frac{\mathrm{d} X_2}{\mathrm{d} t} &=- k X_1^2 X_2 \tag{70}\\
\frac{\mathrm{d} X_3}{\mathrm{d} t} &=2  k X_1^2 X_2 \tag{71}
\end{align}
```

Lets check if units are correct:

```math
\begin{align}
\frac{[M]}{[s]}  &= - \frac{1}{[s][M]^2} [M]^2 [M] \tag{72}\\
\frac{[M]}{[s]}  &= - \frac{1}{[s][M]^2} [M]^2 [M] \tag{73}\\
\frac{[M]}{[s]}  &= \frac{1}{[s][M]^2} [M]^2 [M] \tag{74}
\end{align}
```


"


# ╔═╡ 61a1a25f-21b1-40a6-aafa-17fddf4a848c
md"## 5.2 General formulation of the Mass Action Law

When we have to write the differenetial equations of more complex systems of interactions, such as Eq 21, we can use a general method based omn the following equation: 
```math
\begin{align}
\frac{\mathrm{d} X}{\mathrm{d} t}&= (B-A)^T \cdot K \cdot X^A \tag{75}\\
\end{align}
```


where  `A` and `B` are the matrices with the stoichiometric coefficients, `X` is the state vector $X=[X_1,X_2,X_3]^T$, and `K` is a matrix in the form:

```math
K=\begin{pmatrix}
 k_1 & 0   & 0  & \dots & 0 \\ 
 0 &  k_2  & 0  & \dots & 0\\ 
 0 & 0  & k_3 &\dots  &0 \\ 
 \vdots &  \vdots   &  \vdots  &  \vdots  &  \vdots  \tag{76}\\ 
 0 & 0  & 0  &  \dots & k_r 
\end{pmatrix}
```

$X^A$ is a column vector, calculated as:


```math
X^A=\begin{pmatrix}
X_1^{A_{11}} \cdot X_2^{A_{12}} \cdot X_3^{A_{13}} \cdot  ... \cdot  X_s^{A_{1s}} \\ 
X_1^{A_{21}} \cdot X_2^{A_{22}} \cdot X_3^{A_{23}} \cdot  ... \cdot  X_s^{A_{2s}} \\ 
 \vdots  \,\,\,\,\,\,\,\,\,\,\,  \vdots  \,\,\,\,\,\,\,\,\,\,\,\,    \vdots \,\,\,\,\,\,\,\,\,   \vdots  \,\,\,\,\,\,\,\,   \vdots  \\ 
X_1^{A_{r1}} \cdot X_2^{A_{r2}} \cdot X_3^{A_{r3}} \cdot  ... \cdot  X_s^{A_{rs}} \\ 
\end{pmatrix} =\begin{pmatrix}
 \prod_{i=1}^{s} X_{i}^{A_{1i}}\\ 
\prod_{i=1}^{s} X_{i}^{A_{2i}}\\ 
   \vdots  \tag{77} \\ 
\prod_{i=1}^{s} X_{i}^{A_{ri}}\\ 
\end{pmatrix}
```

"

# ╔═╡ 2b58c174-c501-47ed-ade6-0ee9b88b76eb
md" ### 5.3 Example of differential form of Mass Action Law

As an example we will use the previous reaction for the formation of salt. 

```math
NaCO_3 + CaCl_2 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} CaCO_3 + 2 \cdot NaCl \tag{78}\\  
```



with two reactions $r=2$ and four species $s=4$, with the stoichiometric matrices `A` and `B` listed in Eq. 30"

# ╔═╡ 6c491ea3-10b9-46b3-a3ff-2c27d029efe0
begin
	A=[1 1 0 0; 0 0 1 2]
	B=[0 0 1 2; 1 1 0 0]
	stoichiometric_matrix= (B-A)'
	println("The stoichiometric matrix is $stoichiometric_matrix ")
end

# ╔═╡ eb79f46c-5b6d-4b36-a6e0-f92fe6938b46
md"in our case

```math
K=\begin{pmatrix}
 k_1 & 0   \tag{79}\\ 
 0 &  k_2  \\ 
\end{pmatrix}
```

and 


```math
X^A=\begin{pmatrix}
X_1\cdot X_2  \tag{80}\\ 
X_3\cdot X_4^2  
\end{pmatrix}
```

therefore

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \\ \frac{\mathrm{d} X_4}{\mathrm{d} t}
\end{bmatrix}&= \begin{bmatrix}
 -1   &  1 \\
    -1  &   1\\
     1  &  -1\\
     2  &  -2
\end{bmatrix}\begin{bmatrix}
k_1 & 0\\ 
 0& k_2
\end{bmatrix}\begin{bmatrix}
X_1 \cdot X_2\\ 
X_3 \cdot X_4^2
\end{bmatrix} \tag{81}
\end{align}
```

so calculating

```math
\begin{align}
\begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \\ \frac{\mathrm{d} X_4}{\mathrm{d} t}
\end{bmatrix}&= \begin{bmatrix}
 -1   &  1 \\
    -1  &   1\\
     1  &  -1\\
     2  &  -2
\end{bmatrix}\begin{bmatrix}
k_1 \cdot  X_1 \cdot X_2\\ 
k_2 \cdot X_3 \cdot X_4^2
\end{bmatrix} \tag{82}
\end{align}
```

In other words, the differential equations are obtained by multiplying the matrix $(B-A)'$ by a vector of speeds of the reactions involved. This way, the final set of equations is:

```math
\begin{align}        
            \frac{ \mathrm{d}X_1 }{\mathrm{d}t} &= - k_1 \cdot X_1 \cdot X_2 + k_2 \cdot X_3 \cdot X_4^2  \tag{83}\\ 
            \frac{ \mathrm{d}X_2 }{\mathrm{d}t} &= - k_1 \cdot X_1 \cdot X_2 + k_2 \cdot X_3 \cdot X_4^2   \tag{84}  \\
             \frac{ \mathrm{d} X_3 }{\mathrm{d}t} &=  k_1 \cdot X_1 \cdot X_2 - k_2 \cdot X_3 \cdot X_4^2  \tag{85}  \\
              \frac{ \mathrm{d} X_4 }{\mathrm{d}t} &= 2 k_1 \cdot X_1  \cdot X_2 - 2 k_2 \cdot X_3 \cdot X_4^2   \tag{86}  \\
\end{align}
```
"

# ╔═╡ 3f13dd7b-5430-41b2-a6c2-e185f924ba6b
@htl("""

<div class='blue-background'>
Hello!
</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Computer Task 4: Differential Equations '"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

""")

# ╔═╡ c5a8887c-6875-481d-97bc-cb7da256d2fb
md" 1. Find the order of the reactions, the units of the kinetic rate constants, and write the differential equations that govern the dynamics of the following system:

```math
X_1  \overset{k_1}{\underset{k_2}{\longleftrightarrow}} X_2 
```

2. Find the order of the reactions, the units of the kinetic rate constants, and write the differential equations that govern the dynamics of the following system:

```math 
X_1 + X_2 \overset{k_1}{\underset{k_2}{\longleftrightarrow}} 2 X_1 
```

3.  Find the order of the reactions, the units of the kinetic rate constants, and write the differential equations that govern the dynamics of the following system:

```math 
X_1\overset{k_1}{\longrightarrow}X_2\overset{k_2}{\longrightarrow}X_3
```
"

# ╔═╡ fba4e0b3-52a5-4d1a-b7da-6740f2c2a836
md"## Solution Task 1

Decomposing the reversible reaction into two irreversible reactions
```math
\begin{align*}
X_1  \overset{k_1}{\underset{}{\longrightarrow}} X_2 \\
X_2  \overset{k_2}{\underset{}{\longrightarrow}} X_1 
\end{align*}
```

let's put the stoichiometric values

```math
\begin{align*}
1 \cdot  X_1 + 0 \cdot X_2 \overset{k_1}{\underset{}{\longrightarrow}} 0 \cdot X_1+ 1 \cdot X_2 \\
0 \cdot X_1+ 1 \cdot X_2 \overset{k_2}{\underset{}{\longrightarrow}} 1 \cdot  X_1 + 0 \cdot X_2
\end{align*}
```



the stoichiometric matrices are:

"

# ╔═╡ 0731c83b-3b56-4b80-9654-1760dcbfca9b
AA=[1 0; 0 1];

# ╔═╡ 51ac8bbb-314d-4e82-ab3c-4d65be6f5e74
BB=[0 1; 1 0];

# ╔═╡ 7155b893-ec48-4213-b763-fe4d604efdaa
(BB-AA)'

# ╔═╡ e6a154b3-1f4f-4759-b0ab-62429b0c5200
function Calculate_Keq_stoichiometric_generic(A,B)

 for i=1:size(A)[1]
    println("The order of the reaction $i is ", sum(A[i,:]))
    println("The Units of k$i are s^-1 . k^$(1-(sum(A[i,:])))")
 end
    stoichiometric_matrix= (B-A)'
    println("The stoichiometric matrix is $stoichiometric_matrix ")
    
end

# ╔═╡ 7206af1d-1d79-4c18-82c2-7979d3baa2d5
Calculate_Keq_stoichiometric_generic(AA,BB)

# ╔═╡ b23e2b04-83a3-46f3-86e4-f49d8d42314d
md"The $K$ matrix and the $X^A$ vector required for the equation for the general form of the Mass action are:

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

	"

# ╔═╡ 0ca58de1-28bd-4df1-88ca-efa8d17250f6
md"## Solution Task 2

Decomposing the reversible reaction into two irreversible reactions
```math 
\begin{align}
X_1  + X_2 \overset{k_1}{\underset{}{\longrightarrow}} 2 X_1 \\
 2 X_1  \overset{k_2}{\underset{}{\longrightarrow}}  X_1  + X_2 
\end{align}
```


lest put the stoichiometric values

```math 
\begin{align}
1 \cdot  X_1 + 1 \cdot X_2 \overset{k_1}{\underset{}{\longrightarrow}} 2 \cdot X_1+ 0 \cdot X_2 \\
2 \cdot X_1+ 0 \cdot X_2 \overset{k_2}{\underset{}{\longrightarrow}} 1 \cdot  X_1 + 1 \cdot X_2 
\end{align}
```

the stoichiometric matrices are:"

# ╔═╡ 52a42411-35b9-472a-bd7e-99e7e6395015
AAA=[1 1; 2 0];

# ╔═╡ 637ae77a-a61c-4742-b737-aa4985ce6a32
BBB=[2 0; 1 1];

# ╔═╡ 45f11cc5-59bd-44a6-977b-78e1bea3d9c5
(BBB-AAA)'

# ╔═╡ a21e1203-57a8-4677-968b-c0a90cec9403
Calculate_Keq_stoichiometric_generic(AAA,BBB)

# ╔═╡ a796442b-85fa-45af-8c90-9ce9a8d4fcd2
md"The $K$ matrix and the $X^A$ vector required for the equation for the general form of the Mass action are:

```math
K=\begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix} \tag{26}
```


```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^1\\
X_1^2\cdot X_2^0
\end{pmatrix} = \begin{pmatrix}
X_1 \cdot X_2\\
 X_1^2
\end{pmatrix}  \tag{27}
```

so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
1 & -1\\
 -1 & 1   \end{bmatrix} \begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix}\begin{pmatrix}
X_1 \cdot X_2\\
 X_1^2
\end{pmatrix}  
\end{align}
```

```math
\begin{align}
\begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \end{bmatrix}=  \begin{bmatrix}
1 & -1\\
 -1 & 1   \end{bmatrix}\begin{pmatrix}
k_1 \cdot X_1 \cdot X_2\\
  k_2 \cdot X_1^2
\end{pmatrix} =   \begin{bmatrix}
 k_1 \cdot X_1 \cdot X_2 -  k_2 \cdot X_1^2\\
- k_1 \cdot X_1 \cdot X_2 + k_2 \cdot X_1^2   \end{bmatrix}
\end{align} \tag{28}
```



The differential equations of the system can be written using the differential form of Mass Action Law as: 

```math
\begin{align}
\frac{dX_1}{dt} =  k_1 \cdot X_1 \cdot X_2 -  k_2 \cdot X_1^2 \tag{29}\\
\frac{dX_2}{dt} = - k_1 \cdot X_1 \cdot X_2 + k_2 \cdot X_1^2 \tag{30}\\
\end{align}
```
"

# ╔═╡ 3ac754bf-e6de-417f-9093-35a73f036a1e
md"## Solution Task 3

Decomposing the reversible reaction into two irreversible reactions
```math
\begin{align}
X_1   \overset{k_1}{\underset{}{\longrightarrow}}  X_2 \\
X_2  \overset{k_2}{\underset{}{\longrightarrow}}  X_3  
\end{align}
```

lest put the stoichiometric values

```math
\begin{align}
1 \cdot  X_1 + 0 \cdot X_2 + 0 \cdot X_3 \overset{k_1}{\underset{}{\longrightarrow}} 0 \cdot X_1+ 1 \cdot X_2 + 0 \cdot X_3\\
 0 \cdot X_1+ 1 \cdot X_2 + 0 \cdot X_3 \overset{k_2}{\underset{}{\longrightarrow}}  0 \cdot X_1+ 0 \cdot X_2 + 1 \cdot X_3 
\end{align}
```





the stoichiometric matrices are:"

# ╔═╡ 6c027717-d272-4a3d-824c-54ea64f4fd92
AAAA=[1 0 0 ; 0 1 0];

# ╔═╡ d3950d79-5918-4ed7-bdad-dd04f22d3a30
BBBB=[0 1 0; 0 0 1];

# ╔═╡ 3895430e-9726-4c54-b0d2-9d7901f66f01
Calculate_Keq_stoichiometric_generic(AAAA,BBBB)

# ╔═╡ 3b18856d-27b8-4e1b-bdef-927f3e487e49
md"The $K$ matrix and the $X^A$ vector required for the equation for the general form of the Mass action are:

```math
K=\begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix} \tag{26}
```


```math
X^A=\begin{pmatrix}
X_1^1\cdot X_2^0 \cdot X_3^0\\
X_1^0\cdot X_2^1 \cdot X_3^0
\end{pmatrix} = \begin{pmatrix}
X_1 \\
 X_2
\end{pmatrix}  \tag{27}
```

so 

```math
\begin{align}
 \begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \end{bmatrix}= \begin{bmatrix}
-1 & 0\\
 1 & -1  \\
 0 & 1 \end{bmatrix} \begin{pmatrix}
 k_1 & 0\\
 0 & k_2
\end{pmatrix}\begin{pmatrix}
X_1 \\
 X_2
\end{pmatrix}  
\end{align}
```

```math
\begin{align}
\begin{bmatrix}
\frac{\mathrm{d} X_1}{\mathrm{d} t}\\ \frac{\mathrm{d} X_2}{\mathrm{d} t} \\ \frac{\mathrm{d} X_3}{\mathrm{d} t} \end{bmatrix}=  \begin{bmatrix}
-1 & 0\\
 1 & -1  \\
 0 & 1 \end{bmatrix} \begin{pmatrix}
k_1 X_1 \\
 k_2 X_2
\end{pmatrix}   =   \begin{bmatrix}
 - k_1 \cdot X_1 \\
 k_1 \cdot X_1 - k_2 \cdot X_2  \\
 k_2 \cdot X_2 \end{bmatrix}
\end{align} \tag{28}
```



The differential equations of the system can be written using the differential form of Mass Action Law as: 

```math
\begin{align}
\frac{dX_1}{dt} =  -k_1 \cdot X_1 \\
\frac{dX_2}{dt} =   k_1 \cdot X_1 \cdot X_2 - k_2 \cdot X_2 \tag{30}\\
\frac{dX_3}{dt} =   k_2 \cdot X_2 
\end{align}
```
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"

[compat]
HypertextLiteral = "~0.9.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"
"""

# ╔═╡ Cell order:
# ╟─dd3e39ac-8b2c-11ec-1774-7543d2573dc0
# ╠═7c7a696a-8716-4739-9a40-44a34c621485
# ╠═30a79e11-2165-4908-94a3-d020cc6d2e16
# ╠═e11754b1-6993-403b-a0ec-e0aaf86a5dd3
# ╟─ab0aa754-78e2-47b5-a398-a7e1ab8d6c58
# ╠═d2e26cb1-f136-4b96-8ef3-d3ee30f61b3d
# ╠═10f98ec7-fcad-4587-90f8-62f7b1b74f4d
# ╠═fee636f1-4944-46a4-bcbe-32b7273268c7
# ╠═56a68cf6-6e64-4875-b614-18987857bc06
# ╠═05faae5d-fdb6-4f1d-a77f-0a1337e1666d
# ╟─0d524884-afbf-4594-80ee-9fdfcdabbea0
# ╟─61a1a25f-21b1-40a6-aafa-17fddf4a848c
# ╟─2b58c174-c501-47ed-ade6-0ee9b88b76eb
# ╠═6c491ea3-10b9-46b3-a3ff-2c27d029efe0
# ╟─eb79f46c-5b6d-4b36-a6e0-f92fe6938b46
# ╠═fffc6cac-c044-4cd9-a3b7-084341878faf
# ╟─3f13dd7b-5430-41b2-a6c2-e185f924ba6b
# ╟─c5a8887c-6875-481d-97bc-cb7da256d2fb
# ╟─fba4e0b3-52a5-4d1a-b7da-6740f2c2a836
# ╠═0731c83b-3b56-4b80-9654-1760dcbfca9b
# ╠═51ac8bbb-314d-4e82-ab3c-4d65be6f5e74
# ╠═7155b893-ec48-4213-b763-fe4d604efdaa
# ╠═e6a154b3-1f4f-4759-b0ab-62429b0c5200
# ╠═7206af1d-1d79-4c18-82c2-7979d3baa2d5
# ╟─b23e2b04-83a3-46f3-86e4-f49d8d42314d
# ╟─0ca58de1-28bd-4df1-88ca-efa8d17250f6
# ╠═52a42411-35b9-472a-bd7e-99e7e6395015
# ╠═637ae77a-a61c-4742-b737-aa4985ce6a32
# ╠═45f11cc5-59bd-44a6-977b-78e1bea3d9c5
# ╠═a21e1203-57a8-4677-968b-c0a90cec9403
# ╟─a796442b-85fa-45af-8c90-9ce9a8d4fcd2
# ╟─3ac754bf-e6de-417f-9093-35a73f036a1e
# ╠═6c027717-d272-4a3d-824c-54ea64f4fd92
# ╠═d3950d79-5918-4ed7-bdad-dd04f22d3a30
# ╠═3895430e-9726-4c54-b0d2-9d7901f66f01
# ╟─3b18856d-27b8-4e1b-bdef-927f3e487e49
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
