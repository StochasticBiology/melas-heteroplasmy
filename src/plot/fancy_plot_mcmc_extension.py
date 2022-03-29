import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from itertools import cycle
import matplotlib.lines as mlines

plt.close('all')
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

def acf(theta, lags):
	N = len(theta) # number of observations 
	acf = np.zeros(len(lags)) # container for the acf, for each value of lag
	for i, lag in enumerate(lags): # for a given value of lag
		lag = int(lag)
		theta_t = theta[0:(N-lag)] - np.mean(theta)
		theta_lag = theta[lag:N] - np.mean(theta)
		acf[i] = np.sum(theta_t * theta_lag)/(len(theta_t) * np.var(theta, ddof = 1))
	return acf

demo = int(sys.argv[1])

thin_orig = 1
if demo == 1:
  burnin = 0 
else:
  burnin = int(0.5e7/thin_orig)

report = 1000

pos_dir = '../simulate/'

lag_range = 1000
print('Loading posterior_samples...')
posterior_samples = np.loadtxt(pos_dir+'posterior_samples.dat', delimiter = ',')
print('Loading lps...')
log_posteriors = np.loadtxt(pos_dir+'lps.dat')
#%%
##################################
# Histograms
##################################

print('Plotting Histograms...')

# h_int, k_mrna, f_m, k_m, deg_b_p, ko, kg, zeta, c1, m2, k_gr, k_r, sigma_m_etc, sigma_p_plus, sigma_m_gly, sigma_v, sigma_g, sigma_r_max

mean_names = ['h^*', '\ln(k_{mRNA})', 'f_m', '\ln(k_m)', '\ln(\delta_p)', '\ln(k_o)', '\ln(k_g)', '\ln(\zeta)', 'c_1', 'm_2', 'k_{gr}', '\ln(k_p)']
sigma_names = ['\sigma_{gly}', '\sigma_{ETC}', '\sigma_{P}', '\sigma_V', '\sigma_G', '\sigma_{R}']
all_params = mean_names + sigma_names	

size = 18
plt.rc('font',**{'size':size,'family':'sans-serif','sans-serif':['Helvetica']})

plt.locator_params(nbins=5)


lags = np.linspace(0,lag_range, 50+1)


fig, axs = plt.subplots(nrows=5, ncols=4, figsize=(20,15))
axs = np.reshape(axs,(1,20))[0]

for i, index in enumerate(all_params): # for each index to plot
	ax = axs[i]
	if index == 'm_2':
		ax.hist(posterior_samples[:,i], bins = 100,color = 'red', density = True, alpha = 0.3) # plot value against iterations
		ax.set_xlim([0,20])
	else:
		ax.hist(posterior_samples[:,i], bins = 25,color = 'red', density = True, alpha = 0.3) # plot value against iterations
	

	ax.set_xlabel(r'${}$'.format(all_params[i]), fontsize = 22)
	ax.set_ylabel('Frequency')

##############################
# Plot priors for sigma
##############################
reg_scales = np.loadtxt('../../data/regularization/regularizer_scales.dat')
omega = 2 # reg_scale in met_functions.h

p2 = mlines.Line2D([], [], color='red')
p1 = plt.Rectangle((0, 0), 0.1, 0.1, fc="red", alpha = 0.3, ec = 'red')
leg_size = 9

for i in range(6):
	ax = axs[i+12]
	x = np.linspace(0,max(posterior_samples[:,i+12]))
	lambda_i = reg_scales[i]  * omega
	ax.plot(x, lambda_i*np.exp(-lambda_i*x), '-r')
	ax.legend([p1,p2], ['Marginal Posterior', 'Prior'], prop={'size':leg_size})
##############################


axs[-1].axis('off')
axs[-2].axis('off')

fig.tight_layout()	
plt.savefig('marginal_posteriors.png', dpi = 600)
plt.savefig('marginal_posteriors.pdf')

size = 20
plt.rc('font',**{'size':size,'family':'sans-serif','sans-serif':['Helvetica']})
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(posterior_samples[:,0], bins = 25,color = 'red',density = True, alpha = 0.3) # plot value against iterations
ax.set_xlabel(r'Critical heteroplasmy, ${}$'.format(all_params[0]))
ax.set_ylabel('Density')
fig.tight_layout()	
plt.savefig('marginal_posterior_hint_main.png', dpi = 600)
plt.savefig('marginal_posterior_hint_main.pdf')
#%%

print('Plotting lps...')

# Wrap the array into a 2D array of chunks, truncating the last chunk if 
# chunksize isn't an even divisor of the total size.
# (This part won't use _any_ additional memory)
leg_size = 12

fig, ax = plt.subplots()
chunksize = 100
x = np.array(range(0,len(log_posteriors)))
numchunks = log_posteriors.size // chunksize 
ychunks = log_posteriors[:chunksize*numchunks].reshape((-1, chunksize))
xchunks = x[:chunksize*numchunks].reshape((-1, chunksize))

# Calculate the max, min, and means of chunksize-element chunks...
max_env = ychunks.max(axis=1)
min_env = ychunks.min(axis=1)

ystd = ychunks.std(axis=1, ddof = 1)


ycenters = ychunks.mean(axis=1)
xcenters = xchunks.mean(axis=1)

# Now plot the bounds and the mean...
ax.fill_between(xcenters, min_env, max_env, color='gray', 
                edgecolor='none', alpha=0.5, label = 'range')

ax.fill_between(xcenters, ycenters + ystd,  ycenters - ystd, color='red', 
                edgecolor='none', alpha=0.5 )
ax.plot(xcenters, ycenters, '-k', label = 'mean')

ax.set_xlabel('Iterations')
ax.set_ylabel('Log Posterior')
p1 = plt.Rectangle((0, 0), 1, 1, fc="gray", alpha = 0.5)
p2 = plt.Rectangle((0, 0), 1, 1, fc="red", alpha = 0.5)
p3 = mlines.Line2D([], [], color='black')

plt.legend([p1, p2, p3], ['range', r'$\pm$ std', 'mean'], loc = 'lower right', prop={'size':leg_size})

fig.savefig('log_pos_main.png', dpi=600)
fig.savefig('log_pos_main.pdf')


###################################
# Plot ACFs
###################################
print('Plotting ACFs...')
fig = plt.figure()
ax = fig.add_subplot(111)

lines = ["-","--",":"]
linecycler = cycle(lines)

leg_size = 12
for i, index in enumerate(all_params): # for each index to plot
	ax.plot(lags, acf(posterior_samples[:,i], lags), next(linecycler), label = r'${}$'.format(all_params[i])) # plot acf

ax.set_xlabel('Lag')
ax.set_ylabel('ACF')
ax.legend(loc = 'upper right',prop={'size':leg_size}, ncol = 3)
ax.set_ylim([-0.1,1])

fig.tight_layout()
plt.savefig('acf_main.png', dpi = 600)
plt.savefig('acf_main.pdf')
