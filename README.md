# Covid_modelling
Repository of python code used to model the spread of the COVID-19 virus.

## Modelling Approach
At each point in time the population of a country is made up of three distinct groups: those who are currently infected ($I_t$); those who are susceptible ($S_t$) and those who have recovered ($$"R_t"$$). At each point in time only some fraction of those infected are tested and show a positive result. We denote the number of people infected at time $t$ by $I_t$. We distinguish within this group between those who have tested positive (denoted $I_{st}$) and those who were not tested or incorrectly tested negative ($I_{at}$)  such that $I_t=I_{st}+I_{at}$. (We use the subscripts s and a for these groups because those who were tested were disproportionately those with symptoms while those who were infected but not tested were likely to have had a higher proportion of the asymptomatic). The evolution of $S_t$, $I_t$ and $R_t$ in discrete time is given by the dynamic system:

<img src="https://render.githubusercontent.com/render/math?math=\Delta S_t = - \beta_t I_{t-1}\frac{S_{t-1}}{N}">

<img src="https://render.githubusercontent.com/render/math?math=\Delta R_t = \gamma I_{t-1}">

<img src="https://render.githubusercontent.com/render/math?math=\Delta I_t = \beta_t  I_{t-1} \frac{S_{t-1}}{N} - \gamma I_{t-1},">

$\Delta S$ is the change in the population of the susceptible; $N$ is the total population, $\beta_t$ is the transmission rate of the virus at a time $t$ (the mean number of people an infectious person will infect per unit time) and $\gamma$ is the rate of recovery. The initial infection rate over the infectious period, the reproduction number, is defined as $R_0 = \frac{\beta_t}{\gamma}$. Initially We shall assume that $\pi_a$ is a constant so that

<img src="https://render.githubusercontent.com/render/math?math=I_{st} = (1 - \pi_a) I_t">

and

<img src="https://render.githubusercontent.com/render/math?math=I_{at} = \pi_a I_t.">

The number of new cases at time $t$ ($y_t$) can be calculated as 

<img src="https://render.githubusercontent.com/render/math?math=y_{t} = \Delta I_t + \gamma I_{t-1}.">

New cases are the sum of the change in the number of outstanding cases plus the numbers recovered. The number of new recorded cases ($y_{st}$) is

<img src="https://render.githubusercontent.com/render/math?math=y_{st} = (1 - \pi_a) (\Delta I_t + \gamma I_{t-1}) = (1 - \pi_a) \bigg(\beta_t I_{t-1} \frac{S_{t-1}}{N}\bigg).">


