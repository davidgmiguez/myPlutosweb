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

# ╔═╡ ce8ccbe5-2637-43df-906a-ce6a71e9b918
using DataFrames

# ╔═╡ 87e78056-fecc-4b15-a723-16f51fa19353
using Statistics

# ╔═╡ 7ffa98bf-633f-4307-8da3-c529febde4aa
using HypertextLiteral

# ╔═╡ 25fcbb57-f3ec-4cce-b22c-2503c9bb36e5
using PlutoUI

# ╔═╡ 0ebe3384-97aa-11ec-17b5-cbf23568df03
md" # 2 Mathematical characterization of systems

## 2.1 Ecological Systems

In fact, diversity is related to robustness and adaptation, because if some part o the system perform different tasks, the sayste itself may respond better to sudden environmental changes and laugh changes, as occurs in ecosystems


How do we measure this diversity? is a system composed of parts that are almost identical a diverse system? Is a bag of apples diverse? Depends of what are we looking for, Therefore the measure of diversity is going to depend on the the of system we are dealing with. 


Lets say we want to measure size diversity in a population of bacteria, We can just use standard deviation, Large std, very diverse, small std, no diversity. But lets say now that we focus on the genetic information of a system composed on foxes and  rabbits. ,, there is no way we can quantify the difference in the genes  based on standard deviation (standard deviation of what, the length of the genome, 

"

# ╔═╡ 943c55a4-00d2-4cb7-bba0-a6c46b30e0fa
md"Suppose a biologist wants to measure the diversity of species in a local forest. She collects the following data: "

# ╔═╡ 8ad237fa-8a09-4bae-80f4-ba1298c09581
begin
	dog_slide = @bind 🐶 html"<input type=range min=0 max=100 step=1>"
	mouse_slide = @bind 🐭 html"<input type=range min=0 max=100 step=1>"
	cat_slide = @bind 😺 html"<input type=range min=0 max=100 step=1>"
	fox_slide = @bind 🦊 html"<input type=range min=0 max=100 step=1>"
	bear_slide = @bind 🐻 html"<input type=range min=0 max=100 step=1>"
	
	md"""
	**How many molecules do you have?**
	
	Amount of dogs: $(dog_slide)
	
	Amount of dogs: $(mouse_slide)
	
	Amount of cats: $(cat_slide)
	
	Amount of foxes: $(fox_slide)
	
	Amount of bears: $(bear_slide)
	"""
end

# ╔═╡ bee1dd58-c8e6-4fa6-8eae-7cd4cca0a466
values_animals=DataFrame("Species" => ["Dog", "Mouse", "Cat",  "Fox",  "Bear"],"Frequency" => [🐶, 🐭, 😺, 🦊, 🐻])

# ╔═╡ 68713a25-a962-47e2-8830-c5e5c0511b43
md"Next, she can calculate the total number of organisms"

# ╔═╡ 279b1402-f4ba-4a3b-9d7d-82ac81b515e4
 total_number_of_organisms=sum(values_animals.Frequency)

# ╔═╡ 4190734c-4c6e-4662-9be9-2c19c4956a38
md"A first measure is the Richness $R$, i.e., how many species do we have "

# ╔═╡ 9981448a-edee-435f-9d2b-ec3271c4016a
R=size(values_animals.Frequency,1)

# ╔═╡ 0dbc209e-91dd-4569-afe9-b6407bdab995
md" For ecosystems, or systems composed of a number of individuals from species, we often use the Dominance Index:

```math
D=\frac{\sum n_i(n_i-1)}{N(N-1)}
```

where $n_i$ is the number of organisms in a species i, and N is the total number of organisms 

"

# ╔═╡ c3606f18-d301-44d8-8765-7746bb6493c0
N=sum(values_animals.Frequency)

# ╔═╡ 34a569b4-3f88-4b5f-928b-0edd40c4d6ef
md"the Dominace index is"

# ╔═╡ 3989ead3-2760-443f-80e8-21ac961e2ff5
D = (sum((values_animals.Frequency).* ((values_animals.Frequency) .- 1)) / (N*(N-1)))

# ╔═╡ 7085e737-82db-4879-9c2e-96c6f1f146f0
md"Another highly used idex is the Simpson’s Index of Diversity, which is simply 

```math
S=1-D=1-\frac{\sum n_i(n_i-1)}{N(N-1)}
```

"

# ╔═╡ 7bfac697-d85b-4a60-b90a-2a69c1155f5c
S = 1 - D

# ╔═╡ c71435ea-4bca-4bf7-ba90-576563385b33
md" Species __evenness__ refers to how close in numbers each species in an environment is. Mathematically it is defined as a diversity index, a measure of biodiversity which quantifies how equal the community is numerically. So if there are 40 foxes and 1000 dogs, the community is not very even.

Simpson evenness is closely related to Simpson diversity. It is actually just the Inverse Simpson Diversity with the richness added to the denominator, to mitigate the effect of species count on the range of the metric.
```math
E = \frac{1}{D * R}
```
"

# ╔═╡ 8f31e31d-1080-4b22-b508-b71bf30c2ec0
E= 1/ D / R

# ╔═╡ aa33422d-a8cb-46ed-8e11-d70cfd7caf36
@htl("""

<div class='blue-background'>
Exercise: Measure the Richness, the dominance, teh diversity and the evenness of the followin systems. Discuss and compare the values obtained between the three systems:

</div>

<script>
// more about selecting elements later!
currentScript.previousElementSibling.innerText = "Exercise: Measure the Richness, the dominance, teh diversity and the evenness of the followin systems. Discuss and compare the values obtained between the three systems:"

</script>

<style>
.blue-background {
	padding: .5em;
	background: lightblue;
	color: black;
}
</style>

- Which one is more diverse?

- Which one is more even?

- which one is more rich?

""")

# ╔═╡ 721ad120-8999-4f9f-9bd3-b9d966246436
values_animals_2=DataFrame("Species" => ["A", "B", "C",  "D",  "E", "F", "G", "H", "I", "J"],"Frequency" => [10,9, 11, 10, 8, 12, 10,11,10,9])

# ╔═╡ 8de6461d-e9ba-46df-a347-a7648b6f6e6a
values_animals_3=DataFrame("Species" => ["A", "B", "C",  "D",  "E", "F", "G", "H", "I", "J"],"Frequency" => [72,9, 11, 10, 8, 12, 10,11,10,9])

# ╔═╡ 8a2c17b6-3f5c-4f9f-9714-09df2944dd50
values_animals_4=DataFrame("Species" => ["A", "B", "C",  "D",  "E", "F", "G", "H", "I", "J"],"Frequency" => [33,33, 34, 0, 0, 0, 0, 0,0,0])

# ╔═╡ 02c540f7-cdf8-428d-88c0-64ae6cb905c5
fuji_url = "https://upload.wikimedia.org/wikipedia/commons/thumb/1/12/Mount_Fuji_and_Shinkansen_100_from_Fuji_River.jpg/1200px-Mount_Fuji_and_Shinkansen_100_from_Fuji_River.jpg"

# ╔═╡ 6be4813b-fb90-4fcb-bc93-314971a7b50d
ireland_url = "https://images.ireland.com/media/Images/Down/3e656c481c7a49a2b366d4f8696b17b2.jpg"

# ╔═╡ e25205ff-a4fa-4a54-8a6c-1f94c9dd9ac3
discrete_url = "https://reference.wolfram.com/language/ref/Files/DiscretePlot3D.en/O_3.png"

# ╔═╡ 52f9023e-fb67-43ca-8548-46dfeeb7fb86
continuum_url = "https://i.stack.imgur.com/jGmh5.png"

# ╔═╡ b20cdd19-08c4-4837-947f-97594b8822d2
md" ## 2.2 Graph Theory 
As we mentioned, the important feature of the system comes from the interactions between the parts. Now we will see some tools to study how the parts of a system are interacting and connected.

- Proteins can act as transcription factors for other genes, forming complex networks that ultimately set the identity of a cell.

- When some node of the network changes (for instance, due to mutations), the behavior of a cell can change (leading to genetic diseases, such as cancer).

- Therefore, to understand and treat genetic diseases, we need  to know the genes involved, how they interact with other genes, and how they ultimately affect cellular responses.


"

# ╔═╡ 73f4fd5b-4611-4d90-8f05-af7b279dd86a
md" ### - Type of Networks
The structure of a network determines many relevant aspects of its function.
Network architecture is established by topological and statistical analyses.
Networks can be studied based on its global properties, (size, connectivity, robustness…) 

Depending of the nature of the links between nodes ewe can classify networks in:

- Undirected: links between nodes have no direction (Facebook)
- Directed: links between nodes have direction (Twitter)
- Directed with self-links: Feedback loops (Gene regulatory networks)

[![EC1wSj.md.png](https://iili.io/EC1wSj.md.png)](https://freeimage.host/i/EC1wSj)
"

# ╔═╡ bfa2aecf-b7a8-4d68-aa1b-82d29032bd70
tree_url = "https://er.yuvayana.org/wp-content/uploads/sites/11/2016/02/Binary-Tree-network-model.jpg"

# ╔═╡ 40f52fcb-4be1-40c1-8720-ee48433473f3
md" ### - Size of Networks
The size of a network can refer to the number of nodes N or the number of links E.
The number of links is minimal in a tree network


$(Resource(tree_url)) 

"

# ╔═╡ af0276a5-322a-4a7f-9457-7e58524c5cbd
md"The minimal number of links E in a network of N nodes $E_{min}$ = N-$1$
The maximum number of links depends of the type of network
For undirected

```math
E_{max}=\binom{N}{2}=\frac{N!}{2!(N-2)!)}=\frac{N(N-1)(N-2)...1}{2\cdot 1 \cdot (N-2)(N-3)...1}=\frac{N(N-1)}{2}
```

For directed 

```math
\begin{align*}
E_{max}=2\binom{N}{2}=N+\frac{2 \cdot N!}{2!(N-2)!)}=\\
\\
E_{max}=\frac{2 \cdot N(N-1)(N-2)...1}{2\cdot 1 \cdot (N-2)(N-3)...1}=\\
\\
E_{max}=N(N-1)
\end{align*}
```



For directed with self-links

```math
\begin{align*}
E_{max}=N+2\binom{N}{2}=N+\frac{2 \cdot N!}{2!(N-2)!)}=\\
\\
E_{max}=N+\frac{2 \cdot N(N-1)(N-2)...1}{2\cdot 1 \cdot (N-2)(N-3)...1}=\\
\\
E_{max}=N+N(N-1)=N(1+N-1)=N^2
\end{align*}
```



"

# ╔═╡ d4c87e16-8dff-4818-9bfc-df89c372a498
md" ### - Density of Networks

The density D of a network of N nodes is defined as a ratio of the number of links E to the number of possible links $E_{max}$

```math
\begin{align*}
D=\frac{E-(N-1)}{E_{max}-(N-1)}
\end{align*}
```


"

# ╔═╡ 477ba0d8-1dd9-46c8-bff0-72e4bfa89238
md"Exercise, draw two networks (one directed and one undirected) and calculate the density, which one is mode dense? "

# ╔═╡ 5b7be6f9-18f1-46a0-8a90-1222cf3a4298
distance_to_1=[1,1,2,3,3,4,4]

# ╔═╡ f11db50b-6aa9-4b38-a1c5-7fc217c55e95
md"The average shortest distance to node 1 is"

# ╔═╡ 8a415590-f5d7-4110-b5d9-4b65c8de7ea4
mean(distance_to_1)

# ╔═╡ 1f025dba-f14c-437c-90ce-969bffa75352
md" Exercise: Calculate the average shortest path to node 2:

Solution:
2->1 = 1
2->0 = 2
2->3 = 3
2->7 = 4
2->4 = 4
2->5 = 5
2->6 = 5

0->2 = 2
0->1 = 1
0->3 = 1
0->7 = 2
0->4 = 2
0->4 = 2

"

# ╔═╡ 9179e620-b8ff-4559-b864-c9729a3a42ec
distance_to_2=[1,2,3,4,4,5,5]

# ╔═╡ 9c2a3452-8c7f-4f59-877f-a1c585107e0c
mean(distance_to_2)

# ╔═╡ 2fb50459-ddf7-479e-8a0f-2278c15345ae
md" Exercise: Calculate the average shortest path to node 0:

Solution:

0->2 = 2
0->1 = 1
0->3 = 1
0->7 = 2
0->4 = 2
0->6 = 3
0->5 = 3

"

# ╔═╡ a54318e0-582f-4857-8c7f-2dcd34348ce3
distance_to_0=[2,1,1,2,2,3,3]

# ╔═╡ 794ddca3-4a69-407a-9fe8-2d223775eadd
mean(distance_to_0)

# ╔═╡ 2688137a-c208-4515-a68b-447f03596746
md" Central nodes have a shortest distance to everyone else in teh network than nodes in the boundaries. It is likely that information travels through them."

# ╔═╡ 76ece02b-8170-487b-8e0b-8e32b579b131
md"## Adjacency Matrix: Matrix Representation of Networks

A = binary square matrix of $N$*$N$ elements. 
$A_{ij}$ = 1 if there exists an interaction between $i$ and $j$ 
$A_{ij}$ = 0 if not. 

- undirected graph $A_{ij}$ = $A_{ji}$, $A_{i} = 0$
- directed graph $A_{ij}$ $\isnot$ $A_{ji}$, $A_{i} = 0$,
- directed graph with self-links$A_{ij}$ $\isnot$ $A_{ji}$, $A_{i} $\isnot$ 0$,

Example 1: Draw and characterize the network that correspodns to the following graph,:

```math
\begin{equation*}
A_{i,j} = 
\begin{pmatrix}
0 & 1 & 0 & 1 \\
1 & 0 & 1 & 0 \\
0  & 1  & 0 & 1  \\
1 & 0 & 1 &0 
\end{pmatrix}
\end{equation*}
```

Example 2: Draw and characterize the the network that correspodns to the following graph:

```math
\begin{equation*}
A_{i,j} = 
\begin{pmatrix}
0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 \\
1  & 0  & 0 & 1  \\
0 & 1 & 0 & 0 
\end{pmatrix}
\end{equation*}
```

Example 3: Draw and characterize the the network that correspodns to the following graph:

```math
\begin{equation*}
A_{i,j} = 
\begin{pmatrix}
0 & 1 & 0 & 1 \\
1 & 0 & 1 & 1 \\
0  & 1  & 0 & 0  \\
1 & 1 & 0 & 1 
\end{pmatrix}
\end{equation*}
```

"

# ╔═╡ ddd62a2f-b19d-41b7-8b56-3c279f9c61ab
md" ### - Degree of Networks
Another more direct measure of how connect are the different nodes is the number of links they received. The degree $k$ of a node is the number of connections the node has to other nodes. If a network is directed, such as GRN, nodes have two different degrees:
- in-degree $k_{in}$: the number of incoming edges, and the 
- out-degree $k_{out}$:  the number of outgoing edges
- total-degree $k_{tot}$ = $k_{in}$ + $k_{out}$:

[![EC1wSj.md.png](https://iili.io/EC1wSj.md.png)](https://freeimage.host/i/EC1wSj)

_Undirected_:  
- degree of node 1 is 2
- degree of node 2 is 2  
- degree of node 3 is 3
- degree of node 4 is 1

In an undirected network the total number of links, E, can be expressed as the sum of the node degrees 

```math
\begin{align*}
E=\frac{1}{2}\sum_{i=1}^{N} k_i
\end{align*}
```

the 1/2 factor corrects for the fact that in the sum (2.1) each link is counted twice. For example, the link connecting the nodes 2 and 3 will be counted once in the degree of node 2 and once in the degree of node 3 

_Directed_:

- In-degree of node 1 is 0  
- In-degree of node 2 is 2
- Out-degree of node 1 is 2 
- Out degree of node 2 is 0

In a directed network the total number of links, E, can be expressed as the sum of all node in-degrees or out-degrees 

```math
\begin{align*}
E=\sum_{i=1}^{N} k_{in} = \sum_{i=1}^{N} k_{out}
\end{align*}
```

The __in__ degrees are:

```math
\begin{align*}
k^{in}_1 = 0\\
k^{in}_2 = 2\\
k^{in}_3 = 2\\
k^{in}_4 = 1
\end{align*}
```

The _out_ degrees are:
```math
\begin{align*}
k^{out}_1 = 2\\
k^{out}_2 = 0\\
k^{out}_3 = 2\\
k^{out}_4 = 1
\end{align*}
```

So, for this particular network:
```math
\begin{align*}
E=\sum_{i=1}^{N} k_{in} = \sum_{i=1}^{N} k_{out} =5
\end{align*}
```
"


# ╔═╡ b41588c9-c52c-4e85-ad1d-2775d143bd06
md"

Average degree = average number of links per node. Gives you a basic measurement of how the network is connected in average. 

In undirected networks:

```math
\begin{align*}
\left \langle k \right \rangle = \frac{k_1 + k_2+ k_3 + ... k_N}{N}=\frac{1}{N}\sum_{i=1}^{N}k_i=\frac{2E}{N}
\end{align*}
```

In directed networks:

```math
\begin{align*}
\left \langle k^{in}_i \right \rangle =\frac{1}{N}\sum_{i=1}^{N}k^{in}_i=\left \langle k^{out}_i \right \rangle =\frac{1}{N}\sum_{i=1}^{N}k^{out}_i=\frac{E}{N}
\end{align*}
```

using the notation of the adjacency matrix:
```math
\begin{align*}
\left \langle k \right \rangle =  \frac{ \sum_{i=j}^{N} \sum_{i=1}^{N} A_{ij}} {N}
\end{align*}
```

Exercise: find de average degree of the folliwng networks: 

[![EC1wSj.md.png](https://iili.io/EC1wSj.md.png)](https://freeimage.host/i/EC1wSj)


For the undirected: $\left \langle k \right \rangle = (2+2+3+1)/ 4 = 2$ 

For the directed: $\left \langle k_{in} \right \rangle = (2+2+3+1)/ 4 = 2$ 

"

# ╔═╡ 19f2c67b-fe8f-4913-868c-6fc6dfe32ffb


# ╔═╡ 40804fce-7a1e-420e-94cc-2b1dbdf6e357
md" ## 2.1 Fittness Landscapes

In many cases, the rather abstract and mathematically complex structure of a system of attractors and basins can be replaced by the more intuitive model of a fitness landscape. Under certain mathematical conditions, the deterministic dynamics of a system can be represented by a potential function F on the state space, which attaches a certain number to each state: F: S -> R: s -> F(s), such that trajectory of the system through the state space will always follow the path of steepest descent, i.e. move from a given state s to that neighboring state for which F is minimal. In mechanics, this function corresponds to the potential energy of the system. More generally, the function represents the degree to which a certain state is __preferable__ to another state: the lower the value of F, the __better__ or the more __fit__ the state. Thus, the potential can be seen as the negative of the fitness function. 

The fitness function transforms the state space into a fitness landscape, where every point in the space has a certain __height__ corresponding to its fitness value. This landscape has peaks and valleys. The attractors of the dynamics will now correspond to the local minima of the potential function (potential valleys), or, equivalently, to the maxima of the fitness function (fitness peaks). This can be understood by noting that the system will always move downward in the potential landscape. 

Landscapes are potential functions that define the behavior of the system. Plotting and analyzing landscapes is a great tool to approach a system (for instance , find maximums, minimums of energy, equilibrium points., local versus global maxima on minima...) 

z = z (x, y) 

the state of a system is its position in the _xy_ plane, 
the elevation _z_ is the value of a given property of the system (for instance, temperature…)

usually lanscapes are drawn in sucxh a way that moving up in _z_ (potential function) represents a change in state that requires energy, while moving down is spontaneous. 


## 2.2 Types of landscapes


$(Resource(fuji_url)) Mount Fuji landscapes: One peak
   
$(Resource(ireland_url)) Ripple landscapes: many peaks and valleys
   
$(Resource(discrete_url)) Discrete landscapes: a finite numer of states

$(Resource(continuum_url)) Continuum landscapes: infinite number of states

"

# ╔═╡ 64c9dbd0-edfd-4b1c-b14d-5e9875ea4f33
html"""<iframe width="700" height="400" src="https://www.youtube.com/embed/4pdiAneMMhU?start=83" title="static" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 8ad14cbb-2737-4142-b52b-2a22a3b8352a
md"Static landscapes: no not change z(x,y) do not depend on is not $z$"

# ╔═╡ 1822e04a-9986-42e0-9235-06ef0e486a74
html"""<iframe width="700" height="400" src="https://www.youtube.com/embed/4pdiAneMMhU?start=26" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 9db87f92-f361-465a-8841-fac3fc078343
md"Moving landscapes: z(x,y,z) feedback

A complex system is not the one that has several maxima, minima, are systems where the landscape moves, it depends on variables of the system. The landscape changes depending on the state of the system (Feedback). 

Imagine that the landscape is a mountain, that you want to climb, and that the peak moves to another location as you approach, this is a complex system. 

-  Feedback between you (position) and the mountain (value).
-  Adaptation, you need to change your direction if you want to go up. 

Examples of these landscapes that move:
- Living making decisions
- Evolution: evolutionary landscape
- Development: Fittness landscapes"

# ╔═╡ 42924ebd-0bc9-4bd2-9e22-daba975b19e2
html"""<iframe width="700" height="400" src="https://www.youtube.com/embed/iaq_Fpr4KZc?start=83" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 0e68a81f-d628-40e4-a56f-1dc386247946
md" ## Landscapes in Biology:

Exploration + Exploitation: 

Exploration = Mutations, Genetic  recombination (mix of genes of mother and father in sexual reproduction), (increased variability in a population of bacteria)

Exploitation: natural selection 
(Reduces de variability in a population of bacteria)

Imaging these with feedbacks, with moving landscapes

In moving landscapes, the system never stops exploring, if the system stops exploring, a shitf in the landscape, can convert a peak on a valley, and if you are not exploring you are done

So, complex systems have to always balance between exploring and exploiting

So, here it is the cool thing: complexity itself is an emerging property!
Generated by the need of a system to balance exploitation and exploration. 
"




# ╔═╡ 5b77331c-adfa-4df2-9fb8-c1d14785cc2b
html"""
<iframe width="700" height="400" src="https://www.youtube.com/embed/KUvuv74_E1U" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""


# ╔═╡ c87ffaca-cec4-4aef-8fa8-9d8608941284
md"In this second rendition of Wadington's 'epigenetic landscape', we introduce two variations which we believe reflect some important developmental processes. The first one is that the valleys and the hills, the actual landscape that cells will traverse, are not predermined, that the cells make them as they move along the developmental landscape. The second notion is that at every bifurcation of the landscape there is a binary decision (green or red)"

# ╔═╡ e23af05f-d125-4798-918f-8a20c0f55065
html"""<iframe width="700" height="400" src="https://www.youtube.com/embed/PRmWBBnWfdA" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>"""

# ╔═╡ 979d1f55-870f-424c-a488-71a369c1fd8b
degree_distribution_url = "https://www.researchgate.net/profile/Daniele-Condorelli/publication/242112667/figure/fig3/AS:669311584710678@1536587725639/An-important-graph-property-is-the-degree-distribution-function-P-k-that-describes_W640.jpg"

# ╔═╡ f9f7a01d-fbcf-4e5d-9cab-152eff279e5a
types_of_networks_url="https://www.researchgate.net/publication/343839881/figure/fig2/AS:928266965774337@1598327499590/Four-types-of-networks-in-the-scale-free-network-the-white-and-striped-nodes-represent_W640.jpg"

# ╔═╡ 25576ba0-c98c-41de-ae62-8c5ea6262a6b
md" Degree Distribution (probability of connecting), 

The degree distribution $p_k$ provides the probability that a randomly selected node in the network has degree $p_{k}$. 
Since $p_{k}$ is a probability, it must be normalized, i.e.

```math
\begin{align*}
\sum_{k=1}^{N} p_k =1
\end{align*}
```

P(k) is the probability distribution of the degrees over the whole network.  pk of a network is defined as the fraction of nodes in the network with degree $k_i$. Thus if there are N nodes in total in a network and M of them have degree k, we have $p_{k}$ = M/N.

The degree distribution is probably the most important feature of a network, since it determines many properties of how the network acts. So, the number of nodes in a network with degree k can be obtained by:


```math
\begin{align*}
N_k = N \cdot p_k
\end{align*}
```

And the average degree can be also calculated as:

```math
\begin{align*}
\left \langle k \right \rangle =\frac{1}{N}\sum_{i=1}^{N}k_i=\sum_{k=0}^{\infty }k \cdot p_k
\end{align*}
```
We can characterize the networks based on the degree distribution:
$(Resource(types_of_networks_url))

- __Regular networks and trees__ (usually man-made, have the lowest heterogeneity (e.g. the number of connections each node has is more or less the same) 

- __small-world networks__ higher clustering and almost the same average path than the random networks with the same number of nodes and edges. high modularity, high herogeneity

- __Scale-free networks__ a highly heterogeneous degree distribution (power-law degree distribution)

- __Random networks__, nodes paired with uniform probability (low heterogeneity , since most nodes have the same number of connections), the degree distribution will be a Gaussian bell-shaped curve.

$(Resource(degree_distribution_url))

The typical plot of the degree distribution of a random network as a maximum, which correspond to the characteristic scale of the network. For example, in this network, the characteristic scale is around 8, since there is a higher number of nodes with that number of connections. 

Biological gene regulatory networks have a degree distribution very different from random. There is no characteristic scale (no node with a more common number of connections, no maximum in the degree distribution plot). Many nodes poorly connected, very few nodes (hubs) highly connected.

Degree distributions are often represented in log-log plot. These types of networks are called power-law or scale free networks 

Some nonbiological networks also show a degree distribution that is close to scale free. Above we plot the degree distribution of some of them. 

(A) Actor collaboration graph with N = 212,250,〈k〉 = 28.78. 
(B) Small sample of the World Wide Web, N = 325,729, 〈k〉 = 5.46. 
(C) Power grid data, N = 4941, 〈k〉 = 2.67

What does it mean that so many real networks have a power law (scale free) degree distribution? 
As opposed to random networks, scale-free networks do not have a characteristic scale, meaning that there is no typical node in the network that represents the degree for the other nodes.

### - 2.3 Conclusion: 
The study of Networks is a very big, active and interesting field. The characteristics of networks of interactions are quite useful to understand key properties of large systems where their parts are interacting (big data…). Biological networks are a perfect systems to study using this formalism of network theory. It is well known that most of biologically occurring networks  have a degree distribution close to scale free, characterized by the presence of large hubs:

These power law distributions have the same functional form at all scales.
Why so many different systems have the same scale free degree distribution?

This is quite difficult question to try to explain in such a basic course. But the current agreement is that Scale Free networks often apear in  in systems where you have two properties combined:

 They grow in size (new nodes and links added). They started from small networks and became large by adding new nodes (as oppose to networks were all the links are established at the same time, i.e., a group in WhatsApp)
There is some sort of “preferential attachment”, meaning that there is some popular nodes that every one wants to get information from (rich gets richer)

Both social networks and GRN have these two properties. 


"

# ╔═╡ 6bd76627-d6f0-4adb-af46-203f21f95d27

annealing_url = "https://www.researchgate.net/profile/Omid_Ghasemalizadeh/publication/308786233/figure/fig1/AS:412751978614791@1475419148801/Simulated-Annealing-optimization-of-a-one-dimensional-objective-function.png"

# ╔═╡ 0aa7e21a-255c-4fcc-8ee6-8b3ef96a82b9
md" ## Dynamical systems in lanscape: 

Systems evolve dynamically in a landscape towards the point of minimal energy. They do so by exploring locally in small steps (using pertubrbation and noise).

- Small perturbations, (finds local minima)
- Large perturbations (global minimum, but risky, the system can move too fat away from optimal, and dies)

The process of moving in a landscape is often compared to the process of __Annealing__: a treatment that alters the physical and sometimes chemical properties of a material. It involves heating above recrystallization temperature,  followed by cooling. Heat increases the rate of diffusion, which results in redistributing and eradicating the dislocations in metals. Cooling is then perfomed gradually to reduce gradually the excursions and more local exploration. 

So intially the system can explore larger regions, and then smaller regions to find local minimum.

- $(Resource(annealing_url)) simulated annealing

 This process has inspired many algoritms and optimization methdos, called __Simulated anhealing__: a probabilistic method for finding the minimum of a complex functions using computers. 
"


# ╔═╡ b5174a24-428e-4a41-9b85-fccc7bcf8d32
path_url = "https://cdncontribute.geeksforgeeks.org/wp-content/uploads/exampleFigure-1.png"

# ╔═╡ 6551021e-e6d4-4d34-af80-efd8a7b3010d
md" ### - Distance in Networks

In a real network like the Internet, a short average path length facilitates the quick transfer of information and reduces costs 

- $(Resource(path_url))

Input: source node = $0$ and destination node is = $7$.

Output: Shortest path length is:$2$ -> [0 3 7]

        


Input: source vertex is = $2$ and destination vertex is = $6$ 

Output: Shortest path length is:$5$ -> [2 1 0 3 4 6] 

Exercise: Calculate the average shortest path to node 1:

Solution:
1->2 = 1 
1->0 = 1
1->3 = 2
1->7 = 3
1->4 = 3
1->6 = 4
1->5 = 4
"

# ╔═╡ 6c11b041-72c9-41c7-9fdb-f3257fea0fd9
waddington_url = "https://www.researchgate.net/profile/Vladimir-Kovac/publication/230718157/figure/fig1/AS:300417280954380@1448636469781/The-epigenetic-landscape-proposed-by-C-H-Waddington-1940-in-Organisers-genes.png"

# ╔═╡ 08ff363a-a0e9-4eb2-a118-5374869fc0da
md" ## Waddington epigenetic landscape: Stem cell differentiation

In 1957, Conrad Waddington described that mammalian development is unidirectional, which means that embryonic stem cells develop into a more mature differentiated state. This idea has been depicted as a ball rolling down from the top of Waddington’s mountain to the bottom of a valley (Figure 1a). It demonstrates visually the natural restriction of cell differentiation potential during normal development. 

$(Resource(waddington_url)) Waddington

In this short animation we take the original idea of the 'epigenetic landscape' by CH Waddington and put it into motion, as he probably conceived it in his mind. Waddington saw cells as pebbles rolling down a mountain with hills and valleys. The bottom of the valleys represents cell or differentiation states. Our rendition keeps close to the notions of Waddington except that it looks at the starting point not as a summit but rather as the crater of a volcano"

# ╔═╡ a44b48e7-5497-4d49-a4d7-64508ecb5c3a
adjacency_url="https://www.ebi.ac.uk/training/online/courses/network-analysis-of-protein-interaction-data-an-introduction/wp-content/uploads/sites/64/2020/08/new-fig-4.png"

# ╔═╡ 2db33033-cb6c-4e5b-b952-c9f06c98def0
md"Using the adjacency matrix, the degree is the sum of the row 

For undirected:
```math
\begin{align*}
k_1=\sum_{i=1}^{N} A_{ij}
\end{align*}
```

For directed 
```math
\begin{align*}
k_{1}^{out}=\sum_{i=1}^{N} A_{ij}\\
k_{1}^{in}=\sum_{j=1}^{N} A_{ij}
\end{align*}
```

For weighted 
```math
\begin{align*}
k_{1}^{out}=\sum_{i=1}^{N} A_{ij}\\
k_{1}^{in}=\sum_{j=1}^{N} A_{ij}
\end{align*}
```


- $(Resource(adjacency_url))
"

# ╔═╡ 2fdf8c3b-2c12-42fe-945e-6d91878259ca
https://math.libretexts.org/Bookshelves/Scientific_Computing_Simulations_and_Modeling/Book%3A_Introduction_to_the_Modeling_and_Analysis_of_Complex_Systems_(Sayama)/01%3A_Introduction_to_Modeling_and_Analysis/1.01%3A_Complex_Systems_in_a_Nutshell


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
DataFrames = "~1.3.2"
HypertextLiteral = "~0.9.3"
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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

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

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

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

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═0ebe3384-97aa-11ec-17b5-cbf23568df03
# ╟─943c55a4-00d2-4cb7-bba0-a6c46b30e0fa
# ╠═ce8ccbe5-2637-43df-906a-ce6a71e9b918
# ╠═8ad237fa-8a09-4bae-80f4-ba1298c09581
# ╟─bee1dd58-c8e6-4fa6-8eae-7cd4cca0a466
# ╟─68713a25-a962-47e2-8830-c5e5c0511b43
# ╟─279b1402-f4ba-4a3b-9d7d-82ac81b515e4
# ╟─4190734c-4c6e-4662-9be9-2c19c4956a38
# ╠═9981448a-edee-435f-9d2b-ec3271c4016a
# ╟─0dbc209e-91dd-4569-afe9-b6407bdab995
# ╠═c3606f18-d301-44d8-8765-7746bb6493c0
# ╟─34a569b4-3f88-4b5f-928b-0edd40c4d6ef
# ╠═3989ead3-2760-443f-80e8-21ac961e2ff5
# ╟─7085e737-82db-4879-9c2e-96c6f1f146f0
# ╠═7bfac697-d85b-4a60-b90a-2a69c1155f5c
# ╟─c71435ea-4bca-4bf7-ba90-576563385b33
# ╠═8f31e31d-1080-4b22-b508-b71bf30c2ec0
# ╟─aa33422d-a8cb-46ed-8e11-d70cfd7caf36
# ╠═721ad120-8999-4f9f-9bd3-b9d966246436
# ╠═8de6461d-e9ba-46df-a347-a7648b6f6e6a
# ╠═8a2c17b6-3f5c-4f9f-9714-09df2944dd50
# ╠═02c540f7-cdf8-428d-88c0-64ae6cb905c5
# ╠═6be4813b-fb90-4fcb-bc93-314971a7b50d
# ╠═e25205ff-a4fa-4a54-8a6c-1f94c9dd9ac3
# ╠═52f9023e-fb67-43ca-8548-46dfeeb7fb86
# ╠═b20cdd19-08c4-4837-947f-97594b8822d2
# ╠═73f4fd5b-4611-4d90-8f05-af7b279dd86a
# ╠═40f52fcb-4be1-40c1-8720-ee48433473f3
# ╟─bfa2aecf-b7a8-4d68-aa1b-82d29032bd70
# ╟─af0276a5-322a-4a7f-9457-7e58524c5cbd
# ╠═d4c87e16-8dff-4818-9bfc-df89c372a498
# ╟─477ba0d8-1dd9-46c8-bff0-72e4bfa89238
# ╠═6551021e-e6d4-4d34-af80-efd8a7b3010d
# ╠═5b7be6f9-18f1-46a0-8a90-1222cf3a4298
# ╠═87e78056-fecc-4b15-a723-16f51fa19353
# ╠═f11db50b-6aa9-4b38-a1c5-7fc217c55e95
# ╠═8a415590-f5d7-4110-b5d9-4b65c8de7ea4
# ╟─1f025dba-f14c-437c-90ce-969bffa75352
# ╠═9179e620-b8ff-4559-b864-c9729a3a42ec
# ╠═9c2a3452-8c7f-4f59-877f-a1c585107e0c
# ╟─2fb50459-ddf7-479e-8a0f-2278c15345ae
# ╠═a54318e0-582f-4857-8c7f-2dcd34348ce3
# ╠═794ddca3-4a69-407a-9fe8-2d223775eadd
# ╟─2688137a-c208-4515-a68b-447f03596746
# ╟─76ece02b-8170-487b-8e0b-8e32b579b131
# ╟─ddd62a2f-b19d-41b7-8b56-3c279f9c61ab
# ╟─2db33033-cb6c-4e5b-b952-c9f06c98def0
# ╠═b41588c9-c52c-4e85-ad1d-2775d143bd06
# ╠═25576ba0-c98c-41de-ae62-8c5ea6262a6b
# ╠═19f2c67b-fe8f-4913-868c-6fc6dfe32ffb
# ╠═40804fce-7a1e-420e-94cc-2b1dbdf6e357
# ╠═64c9dbd0-edfd-4b1c-b14d-5e9875ea4f33
# ╟─8ad14cbb-2737-4142-b52b-2a22a3b8352a
# ╟─1822e04a-9986-42e0-9235-06ef0e486a74
# ╟─9db87f92-f361-465a-8841-fac3fc078343
# ╟─0aa7e21a-255c-4fcc-8ee6-8b3ef96a82b9
# ╟─42924ebd-0bc9-4bd2-9e22-daba975b19e2
# ╟─0e68a81f-d628-40e4-a56f-1dc386247946
# ╟─08ff363a-a0e9-4eb2-a118-5374869fc0da
# ╟─5b77331c-adfa-4df2-9fb8-c1d14785cc2b
# ╟─c87ffaca-cec4-4aef-8fa8-9d8608941284
# ╟─e23af05f-d125-4798-918f-8a20c0f55065
# ╠═7ffa98bf-633f-4307-8da3-c529febde4aa
# ╠═25fcbb57-f3ec-4cce-b22c-2503c9bb36e5
# ╟─979d1f55-870f-424c-a488-71a369c1fd8b
# ╠═f9f7a01d-fbcf-4e5d-9cab-152eff279e5a
# ╟─6bd76627-d6f0-4adb-af46-203f21f95d27
# ╠═b5174a24-428e-4a41-9b85-fccc7bcf8d32
# ╠═6c11b041-72c9-41c7-9fdb-f3257fea0fd9
# ╠═a44b48e7-5497-4d49-a4d7-64508ecb5c3a
# ╠═2fdf8c3b-2c12-42fe-945e-6d91878259ca
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
