# evoSSS
Codes and intermediary files for reproducing figures from “Multi-Scale Modeling of Competitive Transmission and Variant Risk Prediction of Infectious Pathogens”, by Sisi Huang, Xiangyun Lu, Zhengyun Xiao, Yiqi Zhuang, Xin Wei, Zuyan Fan, Lingtong Huang, Shijie Zhao, Michael Lynch, Hangping Yao\*, Min Zheng\*, and Chao Jiang\*. (*co-corresponding authors)

## 1. Viral analysis of early COVID-19 patients

[1_Clinical_samples](./1_Clinical_samples/)

Intra-host diversity and mutation dynamics from a collection of 149 early clinical samples from 66 patients.

## 2. Infection and competition experiment

[2_Experiment](./2_Experiment/)

Infection and competition dynamics in multiple cell lines.

## 3. Epidemiological analysis

[3_Epidemiological_analysis](./3_Epidemiological_analysis/)

Epidemiological features of the lineages A and B using ~5.5 million high-quality sequences from the GISAID database.

## 4. evoSSS_model
[4_evoSSS_model](./4_evoSSS_model/)

The evoSSS model integrates the evoSIR model within the Seeding-Spreading framework.

### evoSIR model
The evoSIR model incorporates both the intra- and inter-host modeling, based on the classical Susceptible-Infected-Removed (SIR) model, to simulate the competitive dynamic of infectious diseases.

The model intra-host incorporates terms for the replication rate, competition during co-infection, and viral evolution. The intra-host model of multiple competing variants 1, 2, …, J follows the ordinary differential equations (ODEs):

$$\frac{dV_i}{dt} = r_i V_i \left(1 - \frac{V_i + \sum_{j=1}^{J} \alpha_{ij} V_j}{K} \right)$$

where $V_i$ and $V_j$ are the titers of variants $i$ and $j$, and $r_i$ is the replication rate of variant $i$. $\alpha_{ij}$ represents the inhibitory effect of variant $j$ on variant $i$, and $K$ is the total viral load capacity.

For the inter-host component, we extended the SIR model to include compartments for individuals infected with each variant. $S$ is the number of susceptible individuals correlated with the spreading capacity $Q$, which represents the susceptible contacts per infected individual; $I_j$ is the individuals infected by variant $j$, and $R$ is the number of removed individuals. $N$ is the total population ($N=S+I_1+I_2+R$). The inter-host model follows:

$$
\frac{dS}{dt} = -\sum_{j=1}^{J} \frac{\beta_j S I_j}{N},\ \frac{dI_j}{dt} = \sum_{j=1}^{J} \frac{\beta_j S I_j}{N} - \gamma I_j,\ \frac{dR}{dt} = \gamma \sum_{j=1}^{J} I_j
$$

where $\beta_j\propto V_j/\sum_{j=1}^J V_j$ is the transmission rate for variant $j$, $\gamma$ is the removal (recovery or death) rate for infected individuals. The expected number of new cases for variant $j$ is $\beta_j SI/N$ per day.

### The evoSIR-Seeding-Spreading (evoSSS) model
The evoSIR-Seeding-Spreading (evoSSS) model integrated a Seeding-Spreading framework with the evoSIR model. The transition from one cycle to the next is represented by the composition of the seeding process and the evoSIR model:

$$F=f_{\text{evoSIR}(Q,\beta)}\cdot f_\text{seed}$$

This composition function is applied to the matrix of infectious individuals, transforming the state from cycle $n$ to $n+1$:

$$F: I^n\rightarrow I^{n+1}$$

where $I^n$ is the matrix representing the time sequence of infectious individuals for variants $1,2,…,J$. The seeding process generates infectious seeds $D=[D_1,D_2,…,D_J]$ of cycle $n+1$ from infectious individuals $I=[I_1,I_2,…,I_J]$ of cycle $n$,

$$f_{\text{seed}}: I^n\rightarrow D^{n+1}$$

utilizes a mobility matrix $M=(m_{t_i,t_j})$ to simulate the generation of epidemic hotspots from the infectious seeds at different time points. The elements of the mobility matrix are defined as

$$
m(t_i, t_j) =
\begin{cases} 
0, & t_i < t_j \\
\in [0,1], & t_i \geq t_j
\end{cases}
$$

This process generates infectious seeds

$$D^{n+1}=M*I^n$$

These infectious seeds create epidemic hotspots at each time point. Each hotspot h undergoes the evoSIR modeling 

$$f_{\text{evoSIR}(Q,\beta)}: D^{n+1}\rightarrow I^{n+1,h}$$

where the spreading capacity $Q$ determines the susceptible individuals per seed, while transmission rates $\beta=[\beta_1,\beta_2,…,\beta_J]$ indicate the variant risk. The epidemic outcomes of these hotspots are then aggregated to represent the global pandemic state for the next epidemic cycle:

$$
I_j^{n+1}(t) = \sum_{h \in \text{Hotspots}} I_j^{n+1,h}(t), \quad j \in \{1,2,\dots,J\}.
$$

Specifically, the simulation in this paper used 200 days and an interspace of 30 days for each seeding-spreading cycle. 


## 5. evoSSS_application

[5_evoSSS_application](./5_evoSSS_application/)

To address potential future epidemic outbreaks, we aim to use the evoSSS model to simulate and predict the competition and transmission dynamics of infectious diseases under various scenarios. Spreading capacity, influenced by factors such as social interventions, seasonal outbreaks, and the nature of the virus, limits the number of susceptible individuals. Transmission rates, representing the competitive strength of the variants, serve as early warnings for potentially high-risk variants 



