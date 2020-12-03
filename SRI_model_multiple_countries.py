
# coding: utf-8


import iris
import iris.coord_categorisation as coord_cat 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
import scipy
from scipy import stats
from collections import Counter
from scipy.fftpack import fft, fftfreq
import sys
import matplotlib.animation as animation
import csv
import pandas as pd
from datetime import datetime

#read in necessary constants
#country code
code = sys.argv[1]
#string of lockdown date
lockdown_date = sys.argv[2]
#population
N = float(sys.argv[3])


# Simplest way to read csv (which doesn't work for multi-level columns)
df = pd.read_csv('total-cases-covid-19.csv')

#isolate country that we want, work out when the 1st case was.
GB = df.loc[df['Code'] == code].values[:]
cumulative = GB[:,3]/N

#find index through the year of 1st date
date0 = GB[0,2]
datetime_object = datetime.strptime(date0, '%b %d, %Y')
first_index = int(datetime_object.timetuple().tm_yday)

#rate of cases and new daily cases
GB_new_cases = cumulative[1:] - cumulative[0:-1]

#convert the lockdown date into an index from the first case date
N_lockdown = np.where(GB[:,2] == lockdown_date)[0][0]


fig = plt.figure(figsize = (7,5), dpi = 100)
plt.plot(np.arange(first_index, first_index + len(GB_new_cases)), GB_new_cases*N)
plt.ylabel('new cases')
plt.xticks([1,32,61,92,122],['01/Jan','01/Feb', '01/Mar', '01/Apr', '01/May'])
plt.xlim(15,135)
plt.title('UK Daily New Cases', fontsize = 12)
plt.tight_layout()
#plt.show()
fig.savefig('./covid_optimisations/countries/'+code+'daily_new_cases.png', dpi = 200)



##functions
## function to perform 1 timestep of the SRI model given beta, gamma, N
## and values of populations from previous timestep, returns changes in each probability of population (S, R, I)
def Timestep(N, beta, gamma, S_m, R_m, I_m):
    deltaS = -1 * beta * I_m * S_m/N
    deltaR = gamma * I_m
    deltaI = beta * I_m * S_m/N - gamma * I_m
    
    return deltaS, deltaR, deltaI


#modified version of the run_sim function used to a step_function beta
def run_sim(N, beta0, gamma, N_days, s0, pia, false_symp, beta_step, N_lockdown, mode):
    #initial vals for populations
    R = np.zeros([1])
    symp0 = s0
    symp = np.array([symp0])
    I = np.array([symp[0]*N/(1-pia)])
    S = np.array([N - I[0]])
    symp_new_cases = np.array([symp0])
    #loop over number of timesteps, incrementing S,R,I and symp.
    k = np.log(2)/3
    R_number = np.array([R0*S[0]/N])
    
    for i in range(1, N_days):
        if i < N_lockdown + 14:
            beta = beta0
        else:
            if mode == 'step':
                beta = beta0 - beta_step
            if mode == 'decay':
                k = np.log(2)/3
                beta = beta0 - beta_step * (1 - np.exp(-(i - (N_lockdown + 14))*k))
        if i > len(GB_new_cases) + 10:
            beta = beta0 - beta_step/2.
        if split_betas == True:
            beta = beta*((1-pia) + 0.75*pia)
        deltaS, deltaR, deltaI = Timestep(N, beta, gamma, S[-1], R[-1], I[-1])
        S = np.append(S, S[-1] + deltaS)
        R = np.append(R, R[-1] + deltaR)
        I = np.append(I, I[-1] + deltaI)
        
        R_number = np.append(R_number, (S[i]*beta)/(N*gamma))

        symp = np.append(symp, (1-pia)*I[-1]/N + false_symp*(S[-1]/N + R[-1]/N))
        #calculate number of new symptomatic cases for this timestep
        symp_new_cases = np.append(symp_new_cases, symp[-1] - symp[-2] + (1 - pia)*deltaR/N)
        #print symp[-1], symp[-2], symp_new_cases[-1], deltaR
    return S,R,I,symp,symp_new_cases, R_number

### exploring the variability of the total symptomatic curve with pia
def pia_and_beta_vs_symp(beta0, pia_max, beta_step, pia_step, beta_drop_max, real_symp_new_cases, mode):
    pias = np.arange(0,pia_max,pia_step)
    beta_drops = np.arange(0,beta_drop_max + beta_step, beta_step)
    RMS = np.zeros([len(pias),len(beta_drops)])

    for i in range(len(pias)):
        for j in range(len(beta_drops)):
            print pias[i], beta_drops[j]
            S,R,I,symp,symp_new_cases, R_number = run_sim(N, beta0, gamma, N_days, s0, pias[i], false_symp, beta_drops[j], N_lockdown, mode)
            RMS[i,j] = np.sqrt(np.mean((real_symp_new_cases - symp_new_cases)**2))
    return pias, beta_drops, RMS

#parameters converted to day units of time
half_life = float(sys.argv[5])
gamma = 1 - 0.5**(1/half_life)

N_days = len(GB_new_cases)
s0 = cumulative[0]
false_symp = 0
pia_max = 1
split_betas = False
R0 =  float(sys.argv[4])
beta0 = gamma*R0

#run optimisation
pias, beta_drops, RMS = pia_and_beta_vs_symp(beta0, 1, 0.005, 0.001, beta0, GB_new_cases, 'decay')

#plotting
fig = plt.figure()
cs = plt.contourf( beta0 - beta_drops, pias, RMS, levels = np.linspace(np.min(RMS),np.max(RMS),20), cmap = 'Blues')
min_pia = pias[np.where(RMS == np.min(RMS))[0]]
min_beta = beta0 - beta_drops[np.where(RMS == np.min(RMS))[1]]
bottom5_pia = pias[np.where(RMS < 1.1*np.min(RMS))[0]]
bottom5_beta = beta0 - beta_drops[np.where(RMS < 1.1*np.min(RMS))[1]]
plt.scatter(min_beta, min_pia, color = 'red', s = 4, label = 'min RMS')
plt.scatter(bottom5_beta, bottom5_pia, color = 'g', alpha = 0.6, s = 1.5, label = 'within 10% of min RMS')
plt.ylim(0,1.02)
plt.xlim(0,beta0)
plt.xlabel('Post Lockdown $\\beta_L$')
plt.title('$R_0$ = ' + str(R0))
plt.ylabel('$\pi_a$', fontsize = 15)
plt.colorbar(cs, label = 'RMS')
plt.legend(loc = 3, fontsize = 7)
plt.tight_layout()
#plt.show()
fig.savefig('./covid_optimisations/countries/'+code+'gamma_='+ str(round(gamma, 4)) +'_RMS_R0_=_' + str(R0) +'latest_data.png', dpi = 200)

fig = plt.figure(figsize = (7,4), dpi = 100)
S,R,I,symp,symp_new_cases,R_number = run_sim(N, beta0, gamma, N_days, s0, min_pia, false_symp, beta0 - min_beta, N_lockdown, 'decay')
plt.plot(np.arange(first_index, first_index + len(GB_new_cases)), GB_new_cases)
plt.plot(np.arange(first_index, first_index + len(symp_new_cases)), symp_new_cases)
plt.ylabel('$y_{st}$', fontsize = 15)
plt.xticks([1,32,61,92,122],['01/Jan','01/Feb', '01/Mar', '01/Apr', '01/May'])
plt.xlim(15,135)
plt.title('Best fit: $\pi_a$ = ' + str(round(min_pia[0], 4)) + ', post lockdown $\\beta_L$ = '+ str(round(min_beta[0], 4)))
plt.tight_layout()
#plt.show()
fig.savefig('./covid_optimisations/countries/'+code+'gamma_=_'+ str(round(gamma, 4)) +'_best_fit_R0_=_' + str(R0) +'latest_data.png', dpi = 200)



#now use the optimal values for beta step down and pia to run the model forward in time past present day
n_projection = 180
S,R,I,symp,symp_new_cases,R_number = run_sim(N, beta0, gamma, n_projection, s0, min_pia, false_symp, beta0 - min_beta, N_lockdown, 'decay')
alpha = 0.8
fig = plt.figure(figsize = (7,5), dpi = 100)
plt.plot(np.arange(first_index, n_projection + first_index), symp_new_cases, color = 'g', label = '$y_{st}$')
plt.plot(np.arange(first_index, first_index + len(GB_new_cases)), GB_new_cases, color = 'blue', label = 'actual data')
plt.ylabel('$y_{st}$')
plt.xticks([1,32,61,92,122,152,183],['01/Jan','01/Feb', '01/Mar', '01/Apr', '01/May', '01/Jun', '01/Jul'])
plt.legend(loc = 1)
plt.title('projections with best fit $\pi_a$ and $\\beta_L$ values', fontsize = 12)
plt.tight_layout()
#plt.show()
fig.savefig('./covid_optimisations/countries/'+code+'gamma_=_'+ str(round(gamma, 4)) +'_forward_projection_R0_=_' + str(R0) +'_with_beta_increase_may10th_latest_data.png', dpi = 200)


#now use the optimal values for beta step down and pia to run the model forward in time past present day
n_projection = 180
S,R,I,symp,symp_new_cases,R_number = run_sim(N, beta0, gamma, n_projection, s0, min_pia, false_symp, beta0 - min_beta, N_lockdown, 'decay')
alpha = 0.8
fig = plt.figure(figsize = (7,5), dpi = 100)
plt.plot(np.arange(first_index, n_projection + first_index), R_number, color = 'g', label = '$R_t$')
plt.ylabel('$R_t$')
plt.xticks([1,32,61,92,122,152,183],['01/Jan','01/Feb', '01/Mar', '01/Apr', '01/May', '01/Jun', '01/Jul'])
plt.legend(loc = 1)
plt.title('$R_t$ projections with best fit $\pi_a$ and $\\beta_L$ values', fontsize = 12)
plt.tight_layout()
#plt.show()
fig.savefig('./covid_optimisations/countries/'+code+'gamma_=_'+ str(round(gamma, 4)) +'_forward_projection_of_R_t_R0_=_' + str(R0) +'_with_beta_increase_may10th_latest_data.png', dpi = 200)


#plot just real datan_projection = 180
#### function to perform a bootstrap stats test on two sets of squared residuals data
#### this is done by pooling both residual datasets, randomly selecting
#### 2 random subsets of the same size as original data
def bootstrap_sig_test(beta_drop_min, pia_min, beta_drop_other, pia_other):
    
    #run simulation with min and other parameter combinations
    symp_new_cases_min = run_sim(N, beta0, gamma, N_days, s0, pia_other, false_symp, beta_drop_min, N_lockdown, 'decay')[4]
    symp_new_cases_other = run_sim(N, beta0, gamma, N_days, s0, pia_other, false_symp, beta_drop_other, N_lockdown, 'decay')[4]
    
    #compute square of residuals for both simulations
    sq_residual_best = (symp_new_cases_min - GB_new_cases)**2
    sq_residual_other = (symp_new_cases_other - GB_new_cases)**2
    
    ### real difference in SSW rates
    real_dif =  np.mean(sq_residual_other) - np.mean(sq_residual_best)

    ### pool residual squared datasets
    residual_pooled = np.append(sq_residual_best, sq_residual_other)

    ## loop over bootstrapping, partition into 2 subsets (same size as original SSW arrays)
    ## store difference in SSW rates in dif array
    dif = np.empty(0)
    for j in range(10000):

        ind = np.random.choice(range(residual_pooled.shape[0]), size=N_days, replace=False)
        rest = np.array([i for i in range(0,residual_pooled.shape[0]) if i not in ind])
        
        dummy1 = residual_pooled[ind]
        dummy2 = residual_pooled[rest]
        
        dif = np.append(dif, np.mean(dummy1) - np.mean(dummy2))

    percentiles = np.percentile(dif, np.arange(0,100.001,0.001))
    
    idx = (np.abs(percentiles - real_dif)).argmin()
    
    return dif, real_dif, idx



#pick other point on RMS grid to test against min
pia_other = 0.9
drop = beta_drops[np.where(RMS[int(pia_other*1000),:] == np.min(RMS[int(pia_other*1000),:]))]
print RMS.shape
print RMS[int(pia_other*1000),:]/np.min(RMS)

beta_drop_other = beta0 - drop
dif, real_dif, percentile_index = bootstrap_sig_test(beta0 - min_beta, min_pia, beta_drop_other, pia_other)

fig = plt.figure()
plt.hist(dif, bins=30)  # `density=False` would make counts
plt.ylabel('Frequency')
plt.plot([real_dif,real_dif], [0,1000], label = 'actual difference')
plt.xlabel('Sqaured Residuals Difference')
plt.title(' ')
plt.legend(loc = 2)
fig.savefig('./covid_optimisations/figures_for_paper/bootstrap_RMS_difs' + code + '.png', dpi = 150)
plt.show()

print np.arange(0,100.001,0.001)[percentile_index]