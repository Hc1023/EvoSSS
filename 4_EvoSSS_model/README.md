
## EvoSIR
"EvoSIR" model includes both within-host dynamics and across-host transmission.

### Within-Host Dynamics

We focus on the replication and competition of the viruses within the host. 

For two competing viral variants within a host, \( V_1 \) and \( V_2 \):

\[ \frac{dV_1}{dt} = r_1 V_1 \left(1 - \frac{V_1 + \alpha_{12} V_2}{K}\right) - \mu V_1 \]
\[ \frac{dV_2}{dt} = r_2 V_2 \left(1 - \frac{V_2 + \alpha_{21} V_1}{K}\right)+  \mu V_1  \]

Where:
- \( V_1 \) and \( V_2 \) are the concentrations of variants 1 and 2 within a host.
- \( r_1 \) and \( r_2 \) are the replication rates of variants 1 and 2, respectively.
- \( K \) is the carrying capacity of the host for the viruses.
- \( \alpha_{12} \) represents the competitive effect of variant 2 on variant 1.
- \( \alpha_{21} \) represents the competitive effect of variant 1 on variant 2.

### Across-Host Dynamics

Extending the SIR model for two strains, \( I_1 \) and \( I_2 \), representing the number of individuals infected by variant 1 and 2, respectively:

\[ \frac{dS}{dt} = -\beta_1 S \frac{I_1}{N} - \beta_2 S \frac{I_2}{N} \]
\[ \frac{dI_1}{dt} = \beta_1 S \frac{I_1}{N} - \gamma I_1 \]
\[ \frac{dI_2}{dt} = \beta_2 S \frac{I_2}{N} - \gamma I_2 \]
\[ \frac{dR}{dt} = \gamma I_1 + \gamma I_2 \]

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


To effectively combine your within-host and across-host models into a cohesive framework, you need to establish connections that allow the dynamics of one model to influence the other. This will allow for a more realistic and integrated approach to modeling the spread and evolution of COVID-19. Here’s a step-by-step guide to help you combine the "EvoSIR" model components:

### Step 1: Define Interaction Points
First, determine how the outcomes from the within-host model (viral loads of different strains) affect the parameters in the across-host model (like transmission rates). This typically involves:
- **Transmission Dependency**: Making the transmission rates (\(\beta_1, \beta_2\)) dependent on the viral loads or proportions of each strain within a host at the time of transmission.
- **Recovery Rate Adjustment**: Adjusting recovery rates (\(\gamma_1, \gamma_2\)) based on the dynamics and interactions of the strains within the host.

### Step 2: Modify Across-Host Equations
Modify the across-host transmission equations to incorporate the outputs from the within-host model. For instance:
- **Dynamic Transmission Rates**: Define \(\beta_1\) and \(\beta_2\) as functions of the viral loads \(V_1\) and \(V_2\) from the within-host model. You might consider a simple proportional relationship or a more complex function based on empirical data:
  \[
  \beta_1 = f_1(V_1, V_2), \quad \beta_2 = f_2(V_1, V_2)
  \]
  Where \(f_1\) and \(f_2\) could incorporate factors like viral load peaks, competition effects, and host immune response.

### Step 3: Synchronize Time Scales and Initial Conditions
Ensure that the time scales and initial conditions in both models align:
- **Time Scales**: Both models need to operate over compatible time steps since the interactions and progression of the disease within and between hosts should correspond to the same temporal framework.
- **Initial Conditions**: Set the initial conditions for \(V_1\) and \(V_2\) based on typical initial viral loads observed in newly infected individuals. Similarly, ensure that the initial numbers of \(S, I_1, I_2,\) and \(R\) in the across-host model are realistic.

### Step 4: Numerical Integration and Coupling
Use numerical methods to solve the combined set of differential equations:
- **Coupling Mechanism**: At each time step, update the values of \(V_1\) and \(V_2\) using the within-host model, then feed these values into the across-host model to update \(\beta_1\) and \(\beta_2\), and consequently the states \(S, I_1, I_2,\) and \(R\).
- **Numerical Method**: Methods like Euler’s method or Runge-Kutta methods can be used for solving the coupled equations, adjusting the step size as needed for stability and accuracy.

### Step 5: Simulation and Validation
Run simulations to test the behavior of the combined model under various scenarios:
- **Scenario Testing**: Simulate different public health interventions, such as vaccination or social distancing, and observe how they affect the spread and evolution of the virus.
- **Validation**: Compare the model outputs with real-world data to validate its accuracy. Adjust parameters or model structure based on the fit to data and expert feedback.

By following these steps, you can develop a comprehensive model that not only captures the dynamics of COVID-19 within a host but also how these dynamics influence and are influenced by the transmission of the virus between hosts.