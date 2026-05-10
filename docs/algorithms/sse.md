# Stochastic Series Expansion

In this section, we provide a brief overview of the SSE algorithm. To this end, first consider the canonical ensemble partition function $Z$: 

$$
Z = \text{Tr} (e^{-\beta H}), \ \ \ \rho = e^{-\beta H}
$$

where $\beta = 1/T$ is the reciperocal temperature, $H$ is the Hamiltonian, and we $\rho$ is the imaginary time propagator. In statistical mechanics, we are interested in equilibrium properties which can be evaluated via the usual approach:

$$
\langle O \rangle = \frac{1}{Z} \text{Tr} (O e^{-\beta H}) 
\label{estiamtor_qmc}
$$

where $O$ is the operator representing the observable property of interest. QMC, in general, involves mapping the $d$ dimensional quantum problem to a $d+1$ dimensional classical problem and then constructing a Markov chain to statistically evaluate expection values. For SSE, this is done by Taylor expanding the propagator:

$$
Z = \sum_{\sigma} \sum_{n=0}^{\infty} \frac{(-\beta)^n}{n!} \langle \sigma|H|\sigma\rangle 
\label{taylor_rho_Z}
$$

Here, we have resolved the trace in $S_z$ basis and $\sigma = (\sigma_1, ..., \sigma_N)$ with $\sigma_i \in \{1,-1\} \ \forall \ i \in \{1,...,N\}$. Here, $N$ is the total number of sites in the lattice. Now, for the sake of brevity, we restrict ourselves to the long-range $XY$ model with power-law interactions:

$$
H = - J \sum_{i < j} \frac{1}{r_{ij}^\alpha} \left( S_i^x S_j^x + S_i^y S_j^y \right) = - J \sum_{b=1}^{N_b} \left(H_{1.b} + H_{2,b} \right) + JC
$$

Here, $b$ denotes a single bond in the lattice and $N_b$ represents the total number of bonds. The terms $H_{1,b}$ and $H_{2,b}$ are given by: 

$$
H_{1,b} = \frac{1}{2r_b^\alpha}, \ \ \ H_{2,b} = \frac{1}{2r_b^\alpha} \left(S_{i(b)}^+ S_{j(b)}^- + S_{i(b)}^- S_{j(b)}^+ \right)
$$

The constant offset $JC = \sum_b \frac{J}{2 r_b^\alpha}$ only shifts the spectrum so it can be ignored and added back to the energy at the end of the calculation. We can express each term in the Taylor expansion of the propagator above as a sum over all permutations (the order of the permutation group being determined by the order of the expansion). Then, we can impose a cutoff $n_{max} = M$ on the expansion order and combinatorially re-arrange terms to get:

$$
Z = \sum_{\{\sigma\}, S_M} \frac{(J \beta)^n (M-n)!}{M!} \prod_{k=1}^{M} \langle \sigma_k | H_{a_k, b_k} | \sigma_{k+1} \rangle
$$ 

In the above equation, $S_M = \{(a_1,b_1),...,(a_M,b_M)\}$ with $a_k \in \{0,1,2\}$ and $b_k \in \{1,...,N_b\}$. The operator string has been padded with identity operators wherever needed by introducing the notation $H_{0,b} = I$. Finally, we have also inserted $M-1$ resolutions of identity in the $S_z$ basis between each of the $M$ operators. The Taylor expansion and the insertion of identities effectively maps this to a classical theory in $(d+1)$ dimensions. Note that for each $k$, the many body indices $\sigma_k$ and $\sigma_{k+1}$ can either be the same or differ by a single flipped pair of spins in our construction. 

As with all QMC approaches, we treat the above sum as a set of classical configurations and use transition schemes (called moves) to propagate the system through the resulting state space. To ensure that the corresponding Markov chain converges to the desired distribution, these moves must satisfy detailed balance. In the SSE algorithm, we employ two types of moves: the diagonal transitions and the off-diagonal loop updates. 

## Diagonal Updates

The diagonal update allows us to vary the expansion order $n$ via insertions and removals of diagonal operators $H_{1,b}$ in the operator string. This can be done using the heatbath scheme, which allows us to pre-compute the probability tables and index them in the QMC update implementation efficiently so that the probabilities do not have to be evaluated on the fly. In particular for the long range XY model with power-law interactions, it can be shown that the transition probabilities are given by: 

$$ 
P(n \rightarrow n+1) = \text{min} \left(1, \frac{\beta C}{M-n} \right)
$$

$$
P(n+1 \rightarrow n) = \text{min} \left(1, \frac{M-n+1}{\beta C} \right)
$$

After an insertion update has been accepted for a particular time-slice $n$, the bond is determined by Gibbs sampling the discrete distribution defined by: $P(b) = 1/2r_b^\alpha$. 

## Off-Diagonal Updates 

The off-diagonal updates correspond to replacing non-identity diagonal operators with off-diagonal ones ($H_{2,b}$) probabilistically. This is done using the directed loop technique. In practice, this corresponds to mapping each matrix element in the operator list above, to a set of four vertex leg variables. Each leg denotes a classical spin with value $\sigma_l \in \{-1,1\}$. We then construct a loop in the $(d+1)$ configuration space by entering a vertex at one of its legs and exiting at another leg. This process flips both the entrance and exit legs, leaving the other two leg variables intact. 

We continue traversing the state space after exiting a vertex by moving our pointer to the next time-slice on the same site index with a non-identity operator. This new leg now becomes the entrance leg and same transition approach can be used again until the loop closes. The exit leg, given an entrance leg at each step is determined probabilistically. These probabilities, in general, depend on the model under question. Two common approaches for computing these probabilities are: the heatbath scheme and the directed loop scheme. Typically, the directed loop technique offers a more ergodic traversal of the state space by minimizing bounces (exit leg being the same as the entrance leg). In practice, these off-diagonal update schemes can be implemented efficiently using a linked list data structure, which is the approach used by the Heisenberg Ion package. 

For certain models such as the XY model, bounces can be entirely excluded, and it can be shown that the entire $(d+1)$ dimensional state space can be sub-divided into loops. Moreover, these loops can then be flipped independently with a probability of $1/2$ while respecting detailed balance. This approach is what we refer to as the deterministic algorithm in Heisenberg Ion package. In general, a set of transition probabilities is valid if it satisfies the detailed balance requirement defined by the following:

$$ 
P(s \rightarrow s') W(s) = P(s' \rightarrow s) W(s') 
$$

Here, $s$ and $s'$ are some arbitrary configurations of the entire state space with weights $W(s)$ and $W(s')$ respectively, and $P(s \rightarrow s')$ is the transition probability between those configurations. Both the heatbath and the directed loop schemes satisfy this requirement. 

## Estimators

SSE offers a variety of convenient estimators for common observable quantities of interest, such as the energy and spin stiffness. The equilibrium energy can be computed using the thermodynamic estimator: 

$$ 
E = -\frac{1}{Z} \frac{\partial Z}{\partial \beta}
$$

Using the Taylor expansion form for $Z$ given by Eq. $\eqref{taylor_rho_Z}$, we get the following estimator for the thermodynamic energy:

$$ 
E = \frac{1}{Z} \sum_{\sigma} \sum_{n=0}^{\infty} \left( \frac{-n}{\beta} \right) \frac{(-\beta)^n}{n!} \langle \sigma|H|\sigma\rangle = - \frac{\langle n \rangle}{\beta}
$$ 

The magnetization, defined by $M_z = \frac{1}{N}\sum_{i=1}^{N} \sigma_z^i$ can be computed directly using Eq. $\eqref{estiamtor_qmc}$: 

$$
\langle M_z \rangle = \frac{1}{N} \langle N_{0} - N_{1} \rangle
$$

where $N_0$ is the number of spins in the lattice with $\sigma = 1$ and $N_1$ is the number of spins in the lattice with $\sigma = -1$ at each simulation step for the first time-slice. Since the trace yields periodic boundaries in imaginary time, we could also average this estimator over all time-slices. 

The spin stiffnes is defined as the free energy response to a boundary phase twist $\theta$. The corresponding phase twist Hamiltonian is defined by: 

$$
H(\theta) = -\frac{J}{2} \sum_{i=0}^{N-1} \sum_{r=1}^{(N-1)/2} \frac{1}{r_b^\alpha} \left( e^{-i r \theta} S_i^{+} S_{i+r}^{-} + e^{i r \theta} S_{i}^{-} S_{i+r}^{+} \right)
$$

Note that we have assumed periodic boundaries here since a boundary phase twist would not necessarily be well-defined in the presence of open boundary conditions. The spin stiffness is then given by: 

$$
 \rho_s = \frac{1}{N}\frac{\partial^2 F}{\partial \theta^2} \Biggr|_{\theta = 0}
$$

In general, this stiffness can be related to the winding number $W$ as follows: 

$$
\rho_s = \frac{1}{\beta N^{2-d}} \langle W^2 \rangle
$$

The winding number has a very convenient estimator in SSE. It is given by accumulating the off-diagonal vertices across all time-slices, with a weight of $+r$ for each $S_i^- S_{i+r}^+$ operator, and a weight of $-r$ for each $S_i^{+} S_{i+r}^-$ operator, and $0$ otherwise. Therefore, the evaluation of the spin stiffness in SSE reduces to a simple operator counting estimator.  

In addition to these estimators, the Heisenberg Ion implementation of long range SSE can also be configured to produce shot data. This is produced by recording the entire list of spins at the first time-slice for each simulation step. This should, in principle, allow for the evaluation of all diagonal (in the $S_z$ basis) observables of interest. Finally, for the evaluation of ground state properties in QMC, we can evaluate our observables of interest for increasing values of $\beta$ until the results converge. The converged values of the observables will then correspond to negligible finite temperature systematic error. 