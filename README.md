# Covid_modelling
Repository of python code used to model the spread of the COVID-19 virus.

## Modelling Approach
At each point in time the population of a country is made up of three distinct groups: those who are currently infected (I<sub>t</sub>); those who are susceptible (S<sub>t</sub>) and those who have recovered (R<sub>t</sub>). At each point in time only some fraction of those infected are tested and show a positive result. We denote the number of people infected at time t by I<sub>t</sub>. We distinguish within this group between those who have tested positive (denoted I<sub>st</sub>) and those who were not tested or incorrectly tested negative (I<sub>at</sub>)  such that I<sub>t</sub>=I<sub>st</sub>+I<sub>at</sub>. (We use the subscripts s and a for these groups because those who were tested were disproportionately those with symptoms while those who were infected but not tested were likely to have had a higher proportion of the asymptomatic). The evolution of S<sub>t</sub>, I<sub>t</sub> and R<sub>t</sub> in discrete time is given by the dynamical system:

<img src="https://render.githubusercontent.com/render/math?math=\Delta S_t = - \beta_t I_{t-1}\frac{S_{t-1}}{N}">

<img src="https://render.githubusercontent.com/render/math?math=\Delta R_t = \gamma I_{t-1}">

<img src="https://render.githubusercontent.com/render/math?math=\Delta I_t = \beta_t  I_{t-1} \frac{S_{t-1}}{N} - \gamma I_{t-1},">

&Delta;S is the change in the population of the susceptible; N is the total population, "&beta"<sub>t</sub> is the transmission rate of the virus at a time t (the mean number of people an infectious person will infect per unit time) and &gamma; is the rate of recovery. The initial infection rate over the infectious period, the reproduction number, is defined as R<sub>0</sub> = &beta;<sub>t</sub>/&gamma;. Initially We shall assume that &pi;<sub>a</sub> is a constant so that

<img src="https://render.githubusercontent.com/render/math?math=I_{st} = (1 - \pi_a) I_t">

and

<img src="https://render.githubusercontent.com/render/math?math=I_{at} = \pi_a I_t.">

The number of new cases at time t (y<sub>at</sub>) can be calculated as 

<img src="https://render.githubusercontent.com/render/math?math=y_{t} = \Delta I_t + \gamma I_{t-1}.">

New cases are the sum of the change in the number of outstanding cases plus the numbers recovered. The number of new recorded cases (y<sub>st</sub>) is

<img src="https://render.githubusercontent.com/render/math?math=y_{st} = (1 - \pi_a) (\Delta I_t + \gamma I_{t-1}) = (1 - \pi_a) \bigg(\beta_t I_{t-1} \frac{S_{t-1}}{N}\bigg).">

The strategy that we pursue is to use the data on the numbers of new cases who test positive for the virus. We then seek the values of the parameters of the model - and in particular &pi;<sub>a</sub> - that give a predicted y<sub>st</sub> that matches the data. We use data on the numbers of those who test positive (in the UK and in other countries) as the variable we are trying to match. We use data on tests up the end of April 2020 by which time nearly 500,000 had been tested in the UK and around 175,000 had tested positive (according to data from the Office of National Statistics). It is over this period that new cases testing positive first rose dramatically and then began to fall sharply a few weeks after lockdown began.

To implement the estimation of the model we need to make assumptions about the transmission rate of the virus &beta;<sub>t</sub> and the recovery rate &gamma;. The transmission rate will not have been constant because of policy measures introduced to slow the spread of the infection and because of behavioural changes that were happening even in the runup to the lockdown. In the UK "lockdown", which began on March 23rd, was strict and social distancing was already happening just before this date; both will likely have brought &beta; down significantly. Similar policies were adopted at various times in March 2020 in other countries. We assume a constant value of &beta;<sub>t</sub> before the lockdown date (of &\beta;<sub>0</sub>), followed by a gradual reduction in the &beta;<sub>t</sub> value after this date to simulate the effect the measures have on transmission. The initial value of &beta;<sub>0</sub> is derived from the value of the initial reproduction rate R<sub>0</sub> and the recovery rate &gamma;, using the relation &beta;<sub>0</sub> = &gamma; R<sub>0</sub>.

We assume that after the lockdown date there is a lag until the value of &beta;<sub>t</sub> starts to change from &beta<sub>0</sub>. The lag is between the lockdown measures starting and the impact on the numbers testing positive for the virus. That lag reflects several distinct factors: it must include the lag in the impact on new infections, the lag before symptoms show, the lag before testing the symptomatic and finally the lag before results are known and recorded in the daily measure. We set the overall lag at 14 days, but also assess sensitivity of results to shorter lags in part because social distancing was already happening just before lockdown. After this lag, &beta; decays exponentially towards a value of &beta<sub>L</sub>. The time path for &beta;<sub>t</sub> can be expressed as


<img src="https://render.githubusercontent.com/render/math?math=\beta_t= \beta_0 \text{if}\ t \leq t^*">

<img src="https://render.githubusercontent.com/render/math?math=\beta_t= \beta_0 - (\beta_0 - \beta_L)(1 - e^{-(t-t^*)\lambda}) \text{if } t > t^*">
  
where t<sup>*</sup> is the lockdown time plus the 14 day lag period and &lambda; is the speed of adjustment in &beta; after lockdown measures begin to take effect. We assume that once the lockdown does begin to affect numbers testing positive it quite quickly reaches its full effectiveness, bringing the transmission rate down so that half of its long run impact on &beta; comes through in 3 days, implying that &lambda; = 0.231.

For given values of &gamma;, &beta;<sub>0</sub> and &lambda we search for the values of the two free parameters -  &beta;<sub>L</sub> and &pi;<sub>a</sub> - so as to maximise the fit of the model. We choose those two free parameters to minimise the root mean squared (RMS) deviation between the daily data on the numbers of new positive tests for the virus and the model prediction of that number y<sub>st</sub>. The parameters we fit are a measure of how effective the lockdown is in bringing down the infection rate (measured by how much lower &beta;<sub>L</sub> is relative to &beta;<sub>0</sub>) and the ratio of those infected but not tested to the total population of the infected (&pi;<sub>a</sub>).These key parameters are the ones which best match each country's data on test results. Separate estimation of these parameters for each country allow for cross-country differences in characteristics that might affect the spread of the virus. 

## Files

- Testing_data_UK.csv - Timeseries data on number of tests conducted.
- covid_data_JH.csv - Timeseries data on COVID-19 cases seperated by country.
- SRI_model_multiple_countries.py - SCript contating model and optimisation code designed to run optimisation for data on any country. Run using command "python SRI_model_multiple_countries.py --3 letter country code-- --numeric date of countries lockdown-- --population of country-- . The script produces RMS grids and best fit plots as well as significance testing for minima in RMS.
- SRI_model.ipynb - Jupyter notebook version of model and plotting scripts
- SRI_model_variable_pi.ipynb - version of the SRI that allows the asymptomatic fraction, &pi, vary with time.
  
  
 
