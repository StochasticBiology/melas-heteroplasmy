import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

data_dir = '../../data/from-raw/'

pos_dir = '../simulate/'

trunc = 1

demo = int(sys.argv[1])

if demo == 1:
  # TODO: These have been tweaked for "demo" mode. Please see comment for original values
  thin_orig = 1 # was originally 10
  burnin = 0  # 0.5e7/thin_orig # posterior burnin
  thin_again = 1 # was originally 10  # a second round of thinning for the summary stats
else:
  thin_orig = 10
  burnin = 0.5e7/thin_orig
  thin_again = 10

'''
f = open(pos_dir+'map.dat')
l= f.readlines()
l=l[0]
l=l.split(',')
l[0]=l[0].replace('map: ','')
l[-1]=l[-1].replace('\n','')
map_params = [eval(x) for x in l]
'''

f = open(pos_dir + 'map.dat')
l = f.readlines()
l = l[0]
l = l.split(',')
l[0] = l[0].replace('map: ', '')
l[-1] = l[-1].replace('\n', '')

# TODO: This is the original code, but it's broken. Will hard-code the actual MAP.
#       Need to fix met_hast.c to generate MAP parameters. Could be because of the "demo" mode. Perhaps
#       running with the original parameters fixes this.
# map_params = [eval(x) for x in l]
map_params = [0.517104,0.500965,0.004852,4.036712,-0.741138,-0.722366,-0.964056,0.128937,0.920223,3.192715,0.673631,0.070328,0.200925,0.203106,0.291619,0.063939,0.002678,0.171003]

print('Loading posterior_samples...')
posterior_samples = np.loadtxt(pos_dir + 'posterior_samples.dat', delimiter=',')

print('Loading lps...')
lps = np.loadtxt(pos_dir + 'lps.dat')

################################
# Remove burnin and thin
################################
# posterior_samples = posterior_samples[burnin:]
# posterior_samples = posterior_samples[::thin_again] # true thinning is thus thin_again * thin_orin


logged_params = ['k_mrna',
                 'k_m',
                 'deg_b_p',
                 'ko',
                 'kg',
                 'zeta',
                 'k_r']  # parameters to walk in log space

mean_names = ['h_int', 'k_mrna', 'f_m', 'k_m', 'deg_b_p', 'ko', 'kg', 'zeta', 'c1', 'm2', 'k_gr', 'k_r']
mean_names_dict = {}
for i in range(len(mean_names)):
    mean_names_dict[mean_names[i]] = i


def eval_model_mean(_parameters, **kwargs):
    # {'h_int':0, 'k_mrna':1, 'f_m':2, 'k_m':3, 'deg_b_p':4, 'ko':5, 'kg':6, 'zeta':7, 
    # 'c1':8, 'm2':9, 'k_gr':10, 'k_r':11}
    logged_params = kwargs['logged_params']

    _p = [x for x in _parameters]  # create a local variable, so we don't actually alter _parameters
    #######################################################################
    # Exponentiate some of the variables, so they wonder in log space
    #######################################################################
    for l_param in logged_params:
        log_index = mean_names_dict[l_param]
        _p[log_index] = np.exp(_p[log_index])

    n_plus = (1.0 - h)

    # n_plus_den = np.interp(h, mtDNA_data['Heteroplasmy']/100.0, n_good_den)

    ##################################
    # Shape ETC mRNA degradation
    ##################################
    # Use logistic function
    h_mid = _p[0] - np.log((1.0 - _p[2]) / _p[2]) / _p[3]  # Calculate the corresponding value of the midpoint
    deg_a_m = 1.0 / (1.0 + np.exp(_p[3] * (h - h_mid)))

    ##################################
    # glycolysis mRNA model
    ##################################
    index_h_int = np.argmin(abs(h - _p[0]))
    m_gly = _p[8] * np.ones(len(h))
    c2 = _p[8] - _p[9] * _p[0]
    m_gly[index_h_int:] = _p[9] * h[index_h_int:] + c2

    ##################################
    # Evaluate mRNA, tRNA, protein
    ##################################
    m_etc = _p[7] * n_plus / (deg_a_m + 1.0 / _p[1])
    p_plus = m_etc * n_plus / _p[4]

    ##########################
    # Volume model
    ##########################
    v = _p[5] * p_plus + _p[6] * m_gly

    ##########################
    # Growth model
    ##########################
    g = _p[10] / v

    ##########################
    # Max respiratory rate
    ##########################
    r_max = np.zeros(len(h))
    for i in range(0, len(h)):
        r_max[i] = _p[11] * p_plus[i]

    # return m_etc, p_plus, m_gly, v, g, r_max
    return {'m_etc': m_etc, 'p_plus': p_plus, 'm_gly': m_gly, 'v': v, 'g': g, 'r_max': r_max,
            'deg_a_m': deg_a_m * _p[1]}


################################################
# Import data
################################################

print
'Importing other data...'

# gly_data = pd.read_csv('/home/juvid/Dropbox/Work/Mit_and_Metabolism/Wallace_MELAS/Data/EDA/Transcripts/Rel_GLY.csv')


gly_data = pd.read_csv(data_dir + 'gly_tfm.dat', sep=',')
vol = pd.read_csv(data_dir + 'vol_tfm.dat', sep=',')
etc_data = pd.read_csv(data_dir + 'etc_tfm.dat', sep=',')
petc_data = pd.read_csv(data_dir + 'petc_tfm.dat', sep=',')
growth_data = pd.read_csv(data_dir + 'growth_rate_tfm_fix.dat', sep=',')
resp_data = pd.read_csv(data_dir + 'resp_tfm.dat', sep=',')

data_errors = {'m_gly': 'gly_x_v_err', 'm_etc': 'etc_x_v_err', 'p_plus': 'petc_x_v_err', 'v': 'Err_v'}

values = {'m_gly': 'gly_x_v', 'm_etc': 'etc_x_v', 'p_plus': 'petc_x_v', 'v': 'Cell_Volume', 'g': 'growth',
          'r_max': 'rmax'}  # the column name of the data, in the corresponding dataframe

data = {'m_gly': gly_data[:-1], 'm_etc': etc_data[:-1], 'p_plus': petc_data[:-2], 'v': vol[:-1], 'g': growth_data[:-2],
        'r_max': resp_data[:-2]}
h = np.linspace(0, 0.9, 181)

mtDNA_data = pd.read_csv(data_dir + 'mtDNA_nDNA_f1i.txt', sep=',', skiprows=2,
                         names=['Heteroplasmy', 'mtDNA', 'mtDNA_err'])
mtDNA_data['Heteroplasmy'] = mtDNA_data['Heteroplasmy'] / 100.0
mtDNA_data['n_plus'] = mtDNA_data['mtDNA'] * (1.0 - mtDNA_data['Heteroplasmy'])
mtDNA_data['n_plus_err'] = mtDNA_data['mtDNA_err'] * (1.0 - mtDNA_data['Heteroplasmy'])
merge = pd.merge(mtDNA_data, vol, on='Heteroplasmy', how='inner')
n_good_den = mtDNA_data['n_plus'] / merge['Cell_Volume']
n_good_den_err = np.sqrt(
    mtDNA_data['n_plus_err'] ** 2 / merge['Cell_Volume'] ** 2 + mtDNA_data['n_plus'] ** 2 / merge['Cell_Volume'] ** 4 *
    merge['Err_v'] ** 2)

N_data = np.interp(h, mtDNA_data['Heteroplasmy'], mtDNA_data['mtDNA'])

'''
###################
# Add Dumas data
###################
v_int = np.interp(h, vol['Heteroplasmy']/100.0, vol['Cell_Volume'])
h_0_7_index = np.argmin(abs(h-0.7))
v_0_7 = v_int[h_0_7_index]
dumas = pd.read_csv(data_dir + 'Dumas_CI_CIII_CIV.csv', sep = ',')
dumas['r_CI'] = dumas['CI']/dumas['a-Tubulin']
dumas['r_CIII'] = dumas['CIII']/dumas['a-Tubulin']
dumas['r_CIV'] = dumas['CIV']/dumas['a-Tubulin']

# Note: this is different to what I do for the ETC transcripts, which are added, then normalised. Here I normalise, then add, then normalise again.
# I only do this to be consistent with the data I have access to, in the figure.
# Note: Not able to put an error on this, as we don't know an error on the cell volume at h=0.7
dumas['sum_prot'] = (dumas['r_CI']/dumas['r_CI'][0] + dumas['r_CIII']/dumas['r_CIII'][0] + dumas['r_CIV']/dumas['r_CIV'][0])/3.0
dumas_x_v = v_0_7 * dumas['sum_prot'][1]
'''

############################################################
# Find posterior percentiles
############################################################

print('Finding posterior realisations...')

m_etc_realisations = np.zeros((len(posterior_samples), len(h)))
p_plus_realisations = np.zeros((len(posterior_samples), len(h)))
m_gly_realisations = np.zeros((len(posterior_samples), len(h)))
v_realisations = np.zeros((len(posterior_samples), len(h)))

g_realisations = np.zeros((len(posterior_samples), len(h)))
r_max_realisations = np.zeros((len(posterior_samples), len(h)))

deg_realisations = np.zeros((len(posterior_samples), len(h)))

n_plus_den_realisations = np.zeros((len(posterior_samples), len(h)))

for i, theta_sample in enumerate(posterior_samples):  # for every sample of posterior
    if i % 10000 == 0:
        print(i)
    model = eval_model_mean(theta_sample, logged_params=logged_params, N_data=N_data)  # evaluate the model of the mean

    m_etc = model['m_etc']
    p_plus = model['p_plus']
    m_gly = model['m_gly']
    v = model['v']
    g = model['g']
    r_max = model['r_max']

    deg_a_m = model['deg_a_m']

    m_etc_realisations[i, :] = m_etc
    p_plus_realisations[i, :] = p_plus
    m_gly_realisations[i, :] = m_gly
    v_realisations[i, :] = v
    g_realisations[i, :] = g
    r_max_realisations[i, :] = r_max

    deg_realisations[i, :] = deg_a_m

    n_plus_den_realisations[i, :] = N_data * (1 - h) / v

print('Finding percentiles')

m_etc_percentiles = np.zeros((2, len(h)))
p_plus_percentiles = np.zeros((2, len(h)))
m_gly_percentiles = np.zeros((2, len(h)))
v_percentiles = np.zeros((2, len(h)))
g_percentiles = np.zeros((2, len(h)))
r_max_percentiles = np.zeros((2, len(h)))
deg_percentiles = np.zeros((2, len(h)))

n_plus_den_percentiles = np.zeros((2, len(h)))

m_etc_std = np.zeros(len(h))
p_plus_std = np.zeros(len(h))
m_gly_std = np.zeros(len(h))
v_std = np.zeros(len(h))
g_std = np.zeros(len(h))
r_max_std = np.zeros(len(h))

deg_std = np.zeros(len(h))
n_plus_den_std = np.zeros(len(h))

m_etc_mean = np.zeros(len(h))
p_plus_mean = np.zeros(len(h))
m_gly_mean = np.zeros(len(h))
v_mean = np.zeros(len(h))
g_mean = np.zeros(len(h))
r_max_mean = np.zeros(len(h))
n_plus_den_mean = np.zeros(len(h))

deg_mean = np.zeros(len(h))
deg_median = np.zeros(len(h))

q_low = 25  # lower percentile to plot
q_high = 75  # upper percentile to plot

for i, h_val in enumerate(h):

    if i % 10 == 0: print(i)

    m_etc_percentiles[0, i] = np.percentile(m_etc_realisations[:, i], q_low)
    m_etc_percentiles[1, i] = np.percentile(m_etc_realisations[:, i], q_high)

    p_plus_percentiles[0, i] = np.percentile(p_plus_realisations[:, i], q_low)
    p_plus_percentiles[1, i] = np.percentile(p_plus_realisations[:, i], q_high)

    m_gly_percentiles[0, i] = np.percentile(m_gly_realisations[:, i], q_low)
    m_gly_percentiles[1, i] = np.percentile(m_gly_realisations[:, i], q_high)

    v_percentiles[0, i] = np.percentile(v_realisations[:, i], q_low)
    v_percentiles[1, i] = np.percentile(v_realisations[:, i], q_high)

    g_percentiles[0, i] = np.percentile(g_realisations[:, i], q_low)
    g_percentiles[1, i] = np.percentile(g_realisations[:, i], q_high)

    r_max_percentiles[0, i] = np.percentile(r_max_realisations[:, i], q_low)
    r_max_percentiles[1, i] = np.percentile(r_max_realisations[:, i], q_high)

    deg_percentiles[0, i] = np.percentile(deg_realisations[:, i], q_low)
    deg_percentiles[1, i] = np.percentile(deg_realisations[:, i], q_high)

    n_plus_den_percentiles[0, i] = np.percentile(n_plus_den_realisations[:, i], q_low)
    n_plus_den_percentiles[1, i] = np.percentile(n_plus_den_realisations[:, i], q_high)

    m_etc_std[i] = np.std(m_etc_realisations[:, i], ddof=1)
    p_plus_std[i] = np.std(p_plus_realisations[:, i], ddof=1)
    m_gly_std[i] = np.std(m_gly_realisations[:, i], ddof=1)
    v_std[i] = np.std(v_realisations[:, i], ddof=1)
    g_std[i] = np.std(g_realisations[:, i], ddof=1)
    r_max_std[i] = np.std(r_max_realisations[:, i], ddof=1)

    deg_std = np.std(deg_realisations[:, i], ddof=1)
    n_plus_den_std[i] = np.std(n_plus_den_realisations[:, i], ddof=1)

    m_etc_mean[i] = np.mean(m_etc_realisations[:, i])
    p_plus_mean[i] = np.mean(p_plus_realisations[:, i])
    m_gly_mean[i] = np.mean(m_gly_realisations[:, i])
    v_mean[i] = np.mean(v_realisations[:, i])
    g_mean[i] = np.mean(g_realisations[:, i])
    r_max_mean[i] = np.mean(r_max_realisations[:, i])

    deg_mean[i] = np.mean(deg_realisations[:, i])
    n_plus_den_mean[i] = np.mean(n_plus_den_realisations[:, i])

    deg_median[i] = np.median(deg_realisations[:, i])

#####################
# Get the map
#####################
map_model = eval_model_mean(map_params, logged_params=logged_params, N_data=N_data)
m_etc_map = map_model['m_etc']
p_plus_map = map_model['p_plus']
m_gly_map = map_model['m_gly']
v_map = map_model['v']
g_map = map_model['g']
r_max_map = map_model['r_max']

n_plus_den_map = N_data * (1 - h) / v_map

#####################
# Plot
#####################
leg_size = 5

plt.close('all')
size = 12
plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

fig = plt.figure()

ax = fig.add_subplot(321)
plt.errorbar(etc_data['Heteroplasmy'], etc_data['etc_x_v'], yerr=etc_data['etc_x_v_err'], fmt='kx', markersize=8,
             label='Data')
plt.plot(h, m_etc_map, '-r', label='MAP')
plt.plot(h, m_etc_mean, '-k', label='Posterior Mean')
plt.plot(h, m_etc_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, m_etc_percentiles[1, :], '--r')
plt.fill_between(h, m_etc_percentiles[1, :], m_etc_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC mRNA, $M_{ETC}$')
ax.legend(loc='upper right', prop={'size': leg_size})
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(322)
plt.errorbar(petc_data['Heteroplasmy'][:-1], petc_data[values['p_plus']][:-1],
             yerr=petc_data[data_errors['p_plus']][:-1], fmt='kx', markersize=8, label='Data')
plt.plot(h, p_plus_map, '-r', label='MAP')
plt.plot(h, p_plus_mean, '-k', label='Posterior Mean')
plt.plot(h, p_plus_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, p_plus_percentiles[1, :], '--r')
# plt.plot(0.7, dumas_x_v, 'k*', label = 'Desquiret-Dumas')
plt.fill_between(h, p_plus_percentiles[1, :], p_plus_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC Protein, $P^+$')
ax.legend(loc='upper right', prop={'size': leg_size})
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(323)
plt.errorbar(gly_data['Heteroplasmy'], gly_data[values['m_gly']], yerr=gly_data[data_errors['m_gly']], fmt='kx',
             markersize=8, label='Data')
plt.plot(h, m_gly_map, '-r', label='MAP')
plt.plot(h, m_gly_mean, '-k', label='Posterior Mean')
plt.plot(h, m_gly_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, m_gly_percentiles[1, :], '--r')
plt.fill_between(h, m_gly_percentiles[1, :], m_gly_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Glycolysis mRNA, $M_{gly}$')
ax.legend(loc='upper left', prop={'size': leg_size})
ax.set_xlim([-0.10, 1.1])

ax = fig.add_subplot(324)
plt.errorbar(vol['Heteroplasmy'], vol['Cell_Volume'], yerr=vol['Err_v'], fmt='kx', markersize=8, label='Data')
plt.plot(h, v_map, '-r', label='MAP')
plt.plot(h, v_mean, '-k', label='Posterior Mean')
plt.plot(h, v_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, v_percentiles[1, :], '--r')
plt.fill_between(h, v_percentiles[1, :], v_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Cell volume, $V$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='lower right', prop={'size': leg_size})

ax = fig.add_subplot(325)
plt.errorbar(growth_data['Heteroplasmy'], growth_data['growth'], yerr=growth_data['growth_err'], fmt='kx', markersize=8,
             label='Data')
plt.plot(h, g_map, '-r', label='MAP')
plt.plot(h[:-1], g_mean[:-1], '-k', label='Posterior Mean')
plt.plot(h, g_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, g_percentiles[1, :], '--r')
plt.fill_between(h, g_percentiles[1, :], g_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Growth, $G$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='upper right', prop={'size': leg_size})

ax = fig.add_subplot(326)
plt.plot(resp_data['Heteroplasmy'][:-1], resp_data[values['r_max']][:-1], 'kx', markersize=8, label='Data')
plt.plot(h, r_max_map, '-r', label='MAP')
plt.plot(h, r_max_mean, '-k', label='Posterior Mean')
plt.plot(h, r_max_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, r_max_percentiles[1, :], '--r')
plt.fill_between(h, r_max_percentiles[1, :], r_max_percentiles[0, :], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Max Resp, $R_{max}$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='upper right', prop={'size': leg_size})

fig.tight_layout()

plt.savefig('posterior_model_main.png')
plt.savefig('posterior_model_main.pdf')

################################################################################
# Just the data
################################################################################


fig = plt.figure()

ax = fig.add_subplot(321)
plt.errorbar(etc_data['Heteroplasmy'], etc_data['etc_x_v'], yerr=etc_data['etc_x_v_err'], fmt='kx', markersize=8,
             label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC mRNA, $M_{ETC}$')
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(322)
plt.errorbar(petc_data['Heteroplasmy'][:-1], petc_data[values['p_plus']][:-1],
             yerr=petc_data[data_errors['p_plus']][:-1], fmt='kx', markersize=8, label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC Protein, $P^+$')
ax.set_xlim([0, 1.1])

ax = fig.add_subplot(323)
plt.errorbar(gly_data['Heteroplasmy'], gly_data[values['m_gly']], yerr=gly_data[data_errors['m_gly']], fmt='kx',
             markersize=8, label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Glycolysis mRNA, $M_{gly}$')
ax.set_xlim([-0.10, 1.1])

ax = fig.add_subplot(324)
plt.errorbar(vol['Heteroplasmy'], vol['Cell_Volume'], yerr=vol['Err_v'], fmt='kx', markersize=8, label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Cell volume, $V$')
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(325)
plt.errorbar(growth_data['Heteroplasmy'], growth_data['growth'], yerr=growth_data['growth_err'], fmt='kx', markersize=8,
             label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Growth, $G$')
ax.set_xlim([-0.1, 1.1])
ax.set_ylim([0.6, 1.4])

ax = fig.add_subplot(326)
plt.plot(resp_data['Heteroplasmy'][:-1], resp_data[values['r_max']][:-1], 'kx', markersize=8, label='Data')
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Max Resp, $R_{max}$')
ax.set_xlim([-0.1, 1.1])

fig.tight_layout()

plt.savefig('picard_data.png')
plt.savefig('picard_data.pdf')

############################
# Mean and std
############################

plt.close('all')
size = 12
plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

fig = plt.figure()

ax = fig.add_subplot(321)
plt.errorbar(etc_data['Heteroplasmy'], etc_data['etc_x_v'], yerr=etc_data['etc_x_v_err'], fmt='kx', markersize=8,
             label='Data')
plt.plot(h, m_etc_map, '-k', label='MAP')
plt.plot(h, m_etc_mean, '-r', label='posterior mean')
plt.plot(h, m_etc_mean + m_etc_std, '--r', label=r'$\pm$sd')
plt.plot(h, m_etc_mean - m_etc_std, '--r')
plt.fill_between(h, m_etc_mean + m_etc_std, m_etc_mean - m_etc_std, color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC mRNA, $M_{ETC}$')
ax.legend(loc='upper right', prop={'size': leg_size})
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(322)
plt.errorbar(petc_data['Heteroplasmy'][:-1], petc_data[values['p_plus']][:-1],
             yerr=petc_data[data_errors['p_plus']][:-1], fmt='kx', markersize=8, label='Data')
plt.plot(h, p_plus_map, '-k', label='MAP')
plt.plot(h, p_plus_mean, '-r', label='posterior mean')
plt.plot(h, p_plus_mean + p_plus_std, '--r', label=r'$\pm$sd')
plt.plot(h, p_plus_mean - p_plus_std, '--r')
# plt.plot(0.7, dumas_x_v, 'k*', label = 'Desquiret-Dumas')
plt.fill_between(h, p_plus_mean + p_plus_std, p_plus_mean - p_plus_std, color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'ETC Protein, $P^+$')
ax.legend(loc='upper right', prop={'size': leg_size})
ax.set_xlim([-0.1, 1.1])

ax = fig.add_subplot(323)
plt.errorbar(gly_data['Heteroplasmy'], gly_data[values['m_gly']], yerr=gly_data[data_errors['m_gly']], fmt='kx',
             markersize=8, label='Data')
plt.plot(h, m_gly_map, '-k', label='MAP')
plt.plot(h, m_gly_mean, '-r', label='posterior mean')
plt.plot(h, m_gly_mean + m_gly_std, '--r', label=r'$\pm$sd')
plt.plot(h, m_gly_mean - m_gly_std, '--r')
plt.fill_between(h, m_gly_mean + m_gly_std, m_gly_mean - m_gly_std, color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Glycolysis mRNA, $M_{gly}$')
ax.legend(loc='upper left', prop={'size': leg_size})
ax.set_xlim([-0.10, 1.1])

ax = fig.add_subplot(324)
plt.errorbar(vol['Heteroplasmy'], vol['Cell_Volume'], yerr=vol['Err_v'], fmt='kx', markersize=8, label='Data')
plt.plot(h, v_map, '-k', label='MAP')
plt.plot(h, v_mean, '-r', label='posterior mean')
plt.plot(h, v_mean + v_std, '--r', label=r'$\pm$sd')
plt.plot(h, v_mean - v_std, '--r')
plt.fill_between(h, v_mean + v_std, v_mean - v_std, color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Cell volume, $V$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='lower right', prop={'size': leg_size})

ax = fig.add_subplot(325)
plt.errorbar(growth_data['Heteroplasmy'], growth_data['growth'], yerr=growth_data['growth_err'], fmt='kx', markersize=8,
             label='Data')
plt.plot(h, g_map, '-k', label='MAP')
plt.plot(h, g_mean, '-r', label='posterior mean')
plt.plot(h[:-1], g_mean[:-1] + g_std[:-1], '--r', label=r'$\pm$sd')
plt.plot(h[:-1], g_mean[:-1] - g_std[:-1], '--r')
plt.fill_between(h[:-1], g_mean[:-1] + g_std[:-1], g_mean[:-1] - g_std[:-1], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Growth, $G$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='upper right', prop={'size': leg_size})

ax = fig.add_subplot(326)
plt.plot(resp_data['Heteroplasmy'][:-1], resp_data[values['r_max']][:-1], 'kx', markersize=8, label='Data')
plt.plot(h, r_max_map, '-k', label='MAP')
plt.plot(h, r_max_mean, '-r', label='posterior mean')
plt.plot(h[:-1], r_max_mean[:-1] + r_max_std[:-1], '--r', label=r'$\pm$sd')
plt.plot(h[:-1], r_max_mean[:-1] - r_max_std[:-1], '--r')
plt.fill_between(h[:-1], r_max_mean[:-1] + r_max_std[:-1], r_max_mean[:-1] - r_max_std[:-1], color='red', alpha=0.3)
ax.set_xlabel(r'Heteroplasmy, $h$')
ax.set_ylabel(r'Max Resp, $R_{max}$')
ax.set_xlim([-0.1, 1.1])
ax.legend(loc='upper right', prop={'size': leg_size})

fig.tight_layout()

plt.savefig('mean_std_profile.png')
plt.savefig('mean_std_profile.pdf')

##########################
# degradation
###########################

plt.close('all')
size = 20
leg_size = 12
plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

fig = plt.figure()

plt.plot(h, deg_median, '-k', label='posterior median')
plt.plot(h, deg_percentiles[0, :], '--r', label='Posterior 25-75\% CI')
plt.plot(h, deg_percentiles[1, :], '--r')
plt.fill_between(h, deg_percentiles[1, :], deg_percentiles[0, :], color='red', alpha=0.3)
plt.xlabel(r'Heteroplasmy, $h$')
plt.ylabel(r'mRNA degradation, $\delta_{m}^a$')
plt.xlim([-0.1, 1.1])
plt.legend(loc='upper right', prop={'size': leg_size})

fig.tight_layout()

plt.savefig('deg_profile.png')
plt.savefig('deg_profile.pdf')

##########################
# n_plus_den
###########################


plt.close('all')
size = 20
leg_size = 12
plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
#plt.rc('text', usetex=True)
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

fig = plt.figure()
# plt.errorbar(mtDNA['Heteroplasmy'], n_good_den, yerr=n_good_den_err, fmt='kx', markersize=12, label='Data')  # TODO: Undefined param, needs detective work JA20220325
plt.plot(h, n_plus_den_mean, '-k', label='posterior mean')
plt.plot(h, n_plus_den_map, '-r', label='MAP')
plt.plot(h, n_plus_den_mean + n_plus_den_std, '--r', label=r'$\pm$std')
plt.plot(h, n_plus_den_mean - n_plus_den_std, '--r')
plt.fill_between(h, n_plus_den_mean + n_plus_den_std, n_plus_den_mean - n_plus_den_std, color='red', alpha=0.3)
plt.xlabel(r'Heteroplasmy, $h$')
plt.ylabel(r'Posterior Wt-mtDNA Density, $N^+/V$')
plt.xlim([-0.1, 1.1])
plt.legend(loc='upper right', prop={'size': leg_size})

fig.tight_layout()

plt.savefig('mtdna_den.png')
plt.savefig('mtdna_den.pdf')
