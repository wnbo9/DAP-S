# Simulations

## Simulation with normal distribution

To evaluate the performance of DAP-PIR, we conduct a simulation study with a normal distribution for the genotype data. Specifically, we consider a scenario with $n=500$ individuals and $p=10$ SNPs to relive the computational burden of exact computation of the posterior. Each SNP $g_{ij}$ is generated independently from a standard normal distribution. We randomly select one to five variants and specify them as causal variants with the setting:
$$
g_{ij} \sim N(0,1),
$$

$$
y_i = \sum_{j=1}^{10} \beta_j g_{ij} + \epsilon_i,
$$

$$
\epsilon_i \sim N(0, \tau^{-1}),
$$

$$
\beta_j \sim (1 - \gamma_j)\delta_0 + \gamma_j N\left(0, \frac{\sigma^2_0}{\tau}\right),
$$

where $\beta_j$ represents the effect size of the $j$-th SNP, $\epsilon_i$ is the residual term with varaince $\tau^{-1}$, and $\sigma_0^2/\tau$ is the prior effect size variance. The effect sizes $\beta_j$ are determined by a binary indicator $\gamma_j$: if $\gamma_j=0$, $\beta_j$ is set to zero and indicates a non-causal SNP; if $\gamma_j=1$, $\beta_j$ is drawn from a normal distribution with mean zero and variance $\sigma^2_0/\tau$, representing a causal SNP. The prior probability of each SNP being causal is set to $1/10$ by default, and we randomly assign $S$ causal variants. The indicator vector $\boldsymbol{\gamma} = \left(\gamma_1,\cdots,\gamma_{10}\right)$ is a 10-vector of binary variables such that only $S$ entries are 1 and the other entries are 0. In the simulation, we set $\tau=1$, $\sigma^2_0 = 0.6^2$, and choose $S \in \left\{1,2,3,4,5\right\}$. For each setting, we simulate 1000 times. Then we have simulated $1,000\times 5 = 5,000$ datasets for further investigation.




## Simulation on GTEx V8 data
Using the GTEx V8 data with $n=670$ individual, we select 1,000 genes and $p=1,000$ SNPs near the region of each selected gene. We simulate the phenotypes following SuSiE with various combinations of the number of effects, $S$, and proportion of variance explained (PVE) by genotypes, $\phi$. The simulation settings are as follows:
$$
y_i = \sum_{j=1}^{1000} \beta_j g_{ij}  + \epsilon_i,
$$

$$
\epsilon_i \sim N(0, \sigma^2),
$$

$$
\sigma^2 = \frac{1-\phi}{\phi}\text{Var}(\mathbf{G}\boldsymbol{\beta}),
$$

$$
\beta_j  \sim (1-\gamma_j)\delta_0 + \gamma_jN\left(0, 0.6^2\right).
$$
Specifically, for each gene, we randomly assign $S$ variants to be causal, while their effect sizes are independently drawn from $N(0, 0.6^2)$ and set the other effect sizes to zero. The variance of the error term $\sigma^2$ is given by $\frac{1-\phi}{\phi}\text{Var}(\mathbf{G}\boldsymbol{\beta})$ and the simulated phenotypes follow $N(\mathbf{G}\boldsymbol{\beta}, \sigma^2)$. We follow the SuSiE setting and generate data with pairwise combinations of $S \in \{1,2,3,4,5\}$ and $\phi \in \{0.05,0.1,0.2,0.4\}$. Then we have simulated $1,000\times 5 \times 4 = 20,000$ datasets for further investigation.