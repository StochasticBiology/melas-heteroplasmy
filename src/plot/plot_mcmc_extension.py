from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

thin_orig = 1
burnin = 0#int(0.5e7/thin_orig)

report = 1000

data_dir = Path('../simulate')

lag_range = 1000
print('Loading posterior_samples...')
posterior_samples = np.loadtxt(data_dir/'posterior_samples.dat', delimiter = ',')
print('Loading lps...')
log_posteriors = np.loadtxt(data_dir/'lps.dat')
print('Loading acceptance rates')
acc_rates = np.loadtxt(data_dir/'acc_rates.dat', delimiter = ',')

###########################################
# Perform analysis on the samples
###########################################

def acf(theta, lags):
	N = len(theta) # number of observations 
	acf = np.zeros(len(lags)) # container for the acf, for each value of lag
	for i, lag in enumerate(lags): # for a given value of lag
		lag = int(lag)
		theta_t = theta[0:(N-lag)] - np.mean(theta)
		theta_lag = theta[lag:N] - np.mean(theta)
		acf[i] = np.sum(theta_t * theta_lag)/(len(theta_t) * np.var(theta, ddof = 1))
	return acf

print('Plotting...')

mean_names = ['h_int', 'k_mrna', 'f_m', 'k_m', 'deg_b_p', 'ko', 'kg', 'zeta', 'c1', 'm2', 'k_gr', 'k_d', 'k_r']
sigma_names = ['sigma_gly', 'sigma_etc', 'sigma_petc', 'sigma_vol', 'sigma_g', 'sigma_rmax']
all_params = mean_names + sigma_names	

plt.close('all')
size = 10
#plt.rc('text', usetex=True)
plt.rc('font',**{'size':size,'family':'sans-serif','sans-serif':['Helvetica']})
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

plt.locator_params(nbins=5)


lags = np.linspace(0,lag_range, 50+1)

for i in range(0,6):
	print(i)
	# pick parameters to plot in sets of 3
	params_to_plot = all_params[3*i:3*i+3] # param names 
	indices_to_plot = range(3*i,3*i+3) # param indices in all_params

	fig, axs = plt.subplots(nrows=3, ncols=3)
	axs = np.reshape(axs,(1,9))[0]

	for j, index in enumerate(indices_to_plot): # for each index to plot

		param_name = all_params[index].replace('_','')

		ax1 = axs[3*j]
		ax2 = axs[3*j+1]
		ax3 = axs[3*j+2]

		ax1.plot(range(0,len(posterior_samples[burnin:])), posterior_samples[burnin:,index], '-k') # plot value against iterations
		ax1.set_title(param_name)
		ax1.set_xlabel('Iteration')
		ax1.set_ylabel('Value')
		
		ax2.plot(lags, acf(posterior_samples[burnin:,index], lags), '-k') # plot acf
		ax2.set_title(param_name)
		ax2.set_xlabel('Lag')
		ax2.set_ylabel('ACF')
		
		ax3.hist(posterior_samples[burnin:,index], bins = 25, color = 'black', density = True) # plot value against iterations
		ax3.set_title(param_name)
		ax3.set_xlabel('Value')
		ax3.set_ylabel('Density')
	fig.tight_layout()	
	plt.savefig('out{}.png'.format(i))

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(range(0,len(log_posteriors[burnin:])), log_posteriors[burnin:], '-k')
ax.set_xlabel('Iteration')
ax.set_ylabel('Log posterior')
fig.tight_layout()	
plt.savefig('log_pos.png')

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(acc_rates[int(burnin/report):,0], acc_rates[int(burnin/report):,1], '-k')
ax.set_xlabel('Iteration')
ax.set_ylabel('Acceptance rate')
fig.tight_layout()	
plt.savefig('acc_rate.png')
