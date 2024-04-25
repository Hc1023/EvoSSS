
## EvoSIR
"EvoSIR" model includes both within-host dynamics and across-host transmission.

### Within-Host Dynamics

We focus on the replication and competition of the viruses within the host. 

For two competing viral variants within a host, \( V_1 \) and \( V_2 \):

\[ \frac{dV_1}{dt} = r_1 V_1 \left(1 - \frac{V_1 + \alpha_{12} V_2}{K}\right)  \]
\[ \frac{dV_2}{dt} = r_2 V_2 \left(1 - \frac{V_2 + \alpha_{21} V_1}{K}\right)  \]

Where:
- \( V_1 \) and \( V_2 \) are the concentrations of variants 1 and 2 within a host.
- \( r_1 \) and \( r_2 \) are the replication rates of variants 1 and 2, respectively.
- \( K \) is the carrying capacity of the host for the viruses.
- \( \alpha_{12} \) represents the competitive effect of variant 2 on variant 1.
- \( \alpha_{21} \) represents the competitive effect of variant 1 on variant 2.

### Across-Host Dynamics

Extending the SIR model for two strains, \( I_1 \) and \( I_2 \), representing the number of individuals infected by variant 1 and 2, respectively:

\[ \frac{dS}{dt} = -\beta_1 S \frac{I_1}{N} - \beta_2 S \frac{I_2}{N} \]
\[ \frac{dI_1}{dt} = \beta_1 S \frac{I_1}{N} - \gamma_1 I_1 \]
\[ \frac{dI_2}{dt} = \beta_2 S \frac{I_2}{N} - \gamma_2 I_2 \]
\[ \frac{dR}{dt} = \gamma_1 I_1 + \gamma_2 I_2 \]

Where:
- \( S \) is the number of susceptible individuals.
- \( I_1 \) and \( I_2 \) are the number of individuals infected by variant 1 and 2.
- \( R \) is the number of removed individuals.
- \( \beta_1 \) and \( \beta_2 \) are the transmission rates for variants 1 and 2, influenced by within-host dynamics.
- \( \gamma_1 \) and \( \gamma_2 \) are the removing (recovery or death) rates for individuals infected with variants 1 and 2.
- \( N \) is the total population.

The transmission rates \( \beta_1 \) and \( \beta_2 \) could be functions of the within-host viral loads \( V_1 \) and \( V_2 \), reflecting how the prevalence of each variant within hosts affects its spread between hosts. The exact form of this dependency would be determined based on empirical data or theoretical considerations.

This coupled model, the EvoSIR, now provides a framework for understanding the evolution and spread of competing viral strains both within individual hosts and across the population.

### Bayesian parameter estimation
Use a Bayesian framework to update prior beliefs about within-host parameters based on the observed across-host data:

- Prior Distributions: Start with prior distributions for parameters based on laboratory and clinical data.
- Posterior Estimates: Update these priors with the likelihood of observing the actual epidemiological data under different parameter values, using Bayesian inference methods.
