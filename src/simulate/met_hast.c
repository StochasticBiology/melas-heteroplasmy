/*
Metropolis Hastings with a constant covariance matrix. Requires the Cholesky factor (A), AA^T = Sigma. Fixes bug where prior scales with data.

gcc -o3 met_hast.c -lm -o met_hast.ce
./met_hast.ce

// Parameter ordering convention:
// h_int, k_mrna, f_m, k_m, deg_b_p, ko, kg, zeta, c1, m2, k_gr, k_d, k_r, sigma_m_etc, sigma_p_plus, sigma_m_gly, sigma_v, sigma_g, sigma_r_max
	
Sub-model ordering convention:
m_etc
p_plus
m_gly
v
g
r_max

Author: Juvid Aryaman
Date Created: 24/02/16

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// For 64-bit integer
#include <inttypes.h>
// For limits on value ranges
#include <limits.h>

// Enter debugging mode
//#define VERBOSE

// What machine is this? 0 = cluster, 1 = work, 2 = home
#define MACHINE 1

// Number of iterations for Met Hastings
#define NUM_STEPS 1e10

// Period to thin
#define THINNING 1e5

// Period to print to file. Should be <= NUM_STEPS/THINNING. 
// THINNING * BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define BUFFER_LEN 1e1

// Period to report the acceptance rate (ignores THINNING). Should be < NUM_STEPS
// Typically want to report less often than store params
#define REPORT 1e6

// Period to print acceptance rate to file
// REPORT * REPORT_BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define REPORT_BUFFER_LEN 100

/*
// Number of iterations for Met Hastings
#define NUM_STEPS 1e6

// Period to thin
#define THINNING 1e1

// Period to print to file. Should be <= NUM_STEPS/THINNING. 
// THINNING * BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define BUFFER_LEN 1e4

// Period to report the acceptance rate (ignores THINNING). Should be < NUM_STEPS
// Typically want to report less often than store params
#define REPORT 1e4

// Period to print acceptance rate to file
// REPORT * REPORT_BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define REPORT_BUFFER_LEN 100
*/

/*
// Number of iterations for Met Hastings
#define NUM_STEPS 1e7

// Period to thin
#define THINNING 1e2

// Period to print to file. Should be <= NUM_STEPS/THINNING. 
// THINNING * BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define BUFFER_LEN 1e5

// Period to report the acceptance rate (ignores THINNING). Should be < NUM_STEPS
// Typically want to report less often than store params
#define REPORT 1e4

// Period to print acceptance rate to file
// REPORT * REPORT_BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define REPORT_BUFFER_LEN 100
*/

/*
// Number of iterations for Met Hastings
#define NUM_STEPS 1e4

// How often to thin
#define THINNING 10

// How frequently to print to file. Should be <= NUM_STEPS/THINNING
#define BUFFER_LEN 1e1

// How frequently to report the acceptance rate, given thinning
#define REPORT 1e2

// Period to print acceptance rate to file
// REPORT * REPORT_BUFFER_LEN = Period to print in iterations. About 10^7 is sensible
#define REPORT_BUFFER_LEN 100
*/



// Define a random number
#define RND drand48()

// Number of locations to evaluate h
#define NUM_H 6

// Number of models (m_etc, m_gly, etc...)
#define NUM_MOD 6

// A small number for floating point precision
#define EPSILON 1e-8

//Scale factor for proposal
#define SCALE_FAC 5.76

// Define global variables
double reg_scale;

#include "met_functions.h"

int main( int argc, char *argv[]){
	// Check command line arguments
	int my_seed;


	if(argc != 3){
		printf("Requires (R, seed)\n");
		return 0;
	}
	else{
		reg_scale = atof(argv[1]);		
		printf("Regulariser scale = %f\n", reg_scale);
		my_seed = atoi(argv[2]);		
		printf("Seed = %d\n", my_seed);
		
	}

	int i;
	srand48(my_seed); //randomise the seed of drand48() using command line argument

	// Check that various deinitions do not overflow
	uint64_t num_prints;
	num_prints = NUM_STEPS / THINNING;

	if(num_prints % (uint64_t) BUFFER_LEN != 0){
		printf("Error: num prints indivisible by buffer length\n");
		return 0;
	}

	if(NUM_STEPS > UINT64_MAX){
		printf("NUM_STEPS too big!\n");
		return 0;
	}
	
	if(BUFFER_LEN > LONG_MAX){
		printf("BUFFER_LEN too big!\n");
		return 0;
	}

	if(REPORT > LONG_MAX){
		printf("REPORT too big!\n");
		return 0;
	}

	if( (uint64_t)(NUM_STEPS/REPORT) % (uint64_t)REPORT_BUFFER_LEN != 0){
		printf("Report buffer incorrect!\n");
		return 0;
	}

	if( (uint64_t)(NUM_STEPS/THINNING) > 5e6){
		printf("Too much data!\n");
		return 0;
	}

	printf("Definitions ok...\n");
	printf("Report every %f iterations\n", REPORT);

	///////////////////////////////////////////////////////////////////
	// READ DATA
	///////////////////////////////////////////////////////////////////

	char data_dir[500];
	// Directory to data
	if(MACHINE == 0) { // cluster
		snprintf(data_dir, sizeof(data_dir), "/scratchcomp16/ja1109/Wallace_Data/Wallace_Data_Final/Processed_times_vol/");
	}
	else if (MACHINE == 1) { // work
		snprintf(data_dir, sizeof(data_dir),"../../data/reprocessed");
	}
	else if (MACHINE == 2) { // home
		snprintf(data_dir, sizeof(data_dir), "/media/hdd/home/juvid/Dropbox/Work/Mit_and_Metabolism/Wallace_MELAS/Data/Wallace_Data_Final/Processed_times_vol/");
	}
	else {
		printf("Machine not known!\n"); return 0;
	}
	
	
	

	char chol_dir[] = "../../data/cholesky";
	FILE *ch_p;
	char chol_pth[500];
	snprintf(chol_pth, sizeof(chol_pth), "%schol_flat.dat", chol_dir); 
	printf("Cholesky Factor is in: %s\n", chol_pth);

	ch_p = fopen(chol_pth, "r");

	printf("Data is in: %s\n", data_dir);

	// Pointers for each data file
	FILE *g_p, *m_etc_p, *m_gly_p, *p_plus_p, *r_max_p, *v_p;
	
	// Arrays for the whole path to file. Must be >= length of string
	char g_dir[500];
	char m_etc_dir[500];
	char m_gly_dir[500];
	char p_plus_dir[500];
	char r_max_dir[500];
	char v_dir[500];

	// Concatenate the data directory and the file name
	snprintf(g_dir, sizeof(g_dir), "%sg_fix.dat", data_dir);
	snprintf(m_etc_dir, sizeof(m_etc_dir), "%sm_etc.dat", data_dir);
	snprintf(m_gly_dir, sizeof(m_gly_dir), "%sm_gly.dat", data_dir);
	snprintf(p_plus_dir, sizeof(p_plus_dir), "%sp_plus.dat", data_dir);
	snprintf(r_max_dir, sizeof(r_max_dir), "%sr_max.dat", data_dir);
	snprintf(v_dir, sizeof(v_dir), "%sv.dat", data_dir);

	// Open pointers to each file
	g_p = fopen(g_dir, "r");
	m_etc_p = fopen(m_etc_dir, "r");
	m_gly_p = fopen(m_gly_dir, "r");
	p_plus_p = fopen(p_plus_dir, "r");
	r_max_p = fopen(r_max_dir, "r");
	v_p = fopen(v_dir, "r");

	// Number of lines in each file = (amount of data * 2)
	int g_len = 12;
	int m_etc_len = 14;
	int m_gly_len = 14;
	int p_plus_len = 10;
	int r_max_len = 12;
	int v_len = 14;

	// Truncate datasets by not reading in h = 1 data
	g_len  = g_len - 2;
	m_etc_len  = m_etc_len - 2;
	m_gly_len  = m_gly_len - 2;
	p_plus_len  = p_plus_len - 2;
	r_max_len  = r_max_len - 2;
	v_len  = v_len - 2;

	// Containers for data
	double g_data[g_len];
	double m_etc_data[m_etc_len];
	double m_gly_data[m_gly_len];
	double p_plus_data[p_plus_len];
	double r_max_data[r_max_len];
	double v_data[v_len];

	
	// Read data into arrays line by line	
	for (i=0; i < g_len; i++){
		fscanf(g_p, "%lf\n", &g_data[i]);
	}
	for (i=0; i < m_etc_len; i++){
		fscanf(m_etc_p, "%lf\n", &m_etc_data[i]);
	}
	for (i=0; i < m_gly_len; i++){
		fscanf(m_gly_p, "%lf\n", &m_gly_data[i]);
	}
	for (i=0; i < p_plus_len; i++){
		fscanf(p_plus_p, "%lf\n", &p_plus_data[i]);
	}
	for (i=0; i < r_max_len; i++){
		fscanf(r_max_p, "%lf\n", &r_max_data[i]);
	}
	for (i=0; i < v_len; i++){
		fscanf(v_p, "%lf\n", &v_data[i]);
	}

	#ifdef VERBOSE
		// Print data
		printf("g\n");	
		for (i=0; i < g_len; i++){
			printf("%lf\n", g_data[i]);
		}
		printf("\n");
		printf("m_etc\n");		
		for (i=0; i < m_etc_len; i++){
			printf("%lf\n", m_etc_data[i]);
		}
		printf("\n");
		printf("m_gly\n");		
		for (i=0; i < m_gly_len; i++){
			printf("%lf\n", m_gly_data[i]);
		}
		printf("\n");
		printf("p_plus\n");		
		for (i=0; i < p_plus_len; i++){
			printf("%lf\n", p_plus_data[i]);
		}
		printf("\n");
		printf("r_max\n");		
		for (i=0; i < r_max_len; i++){
			printf("%lf\n", r_max_data[i]);
		}
		printf("\n");
		printf("v\n");		
		for (i=0; i < v_len; i++){
			printf("%lf\n", v_data[i]);
		}
		printf("\n");
	#endif
	
	// Close the pointers
	fclose(m_etc_p);
	fclose(p_plus_p);
	fclose(m_gly_p);
	fclose(v_p);
	fclose(g_p);
	fclose(r_max_p);
	
	printf("Data read ok...\n");

	//////////////////////////////////////////////////////////////////////////
	// INITIALISE VARIABLES
	//////////////////////////////////////////////////////////////////////////

	double h[] = {0.0, 0.2, 0.3, 0.5, 0.6, 0.9}; // all possible discrete values of h
	
	// Parameter ordering convention:
	// h_int, k_mrna, f_m, k_m, deg_b_p, ko, kg, zeta, c1, m2, k_gr, k_r, sigma_m_etc, sigma_p_plus, sigma_m_gly, sigma_v, sigma_g, sigma_r_max
	
	// Parameter initialisation
	double params[] = {0.37319628,  0.00434872,  0.34545602,  4.06544628, -0.49812794,
       -0.89170928, -0.91006802,  0.60002103,  0.94414452,  1.98640131,
        0.65711522,  0.34257764,  0.35816268,  0.22491422,  0.10311401,
        0.06982591,  0.01275206,  0.26357679
        }; // initial values of parameters
    int num_params = sizeof(params) / sizeof(params[0]);

	/*
	logged_params = ['k_mrna',
    'k_m',
    'deg_b_p',
    'ko',
    'kg',
    'zeta',    
    'k_gr',
    'k_d',
    'k_r']
	*/

    int logged_params[] = {1, 3, 4, 5, 6, 7, 11};
    int num_log_params = sizeof(logged_params) / sizeof(logged_params[0]);


	// Priors in normal space
	double mean_min_priors[] = {
	0.0, // h_int
	1e-10, // k_mrna
	0.0, // f_m
	1e-10, // k_m
	1e-10,  // deg_b_p
	1e-10, // ko
	1e-10, // kg
	1e-10, // zeta
	-10.0, // c1
	0.0, // m2
	1e-10, // k_gr
	1e-10 // k_r
	};

	double mean_max_priors[] = {
	1.0, // h_int
	100.0, // k_mrna
	1.0, // f_m
	100.0, // k_m
	100.0, // deg_b_p
	100.0, // ko
	100.0, // kg
	100.0, // zeta
	10.0, // c1
	100.0, // m2
	100.0, // k_gr
	100.0 // k_r
	};

	// Log transform boundaries of the appropriate parameters, to make uniforms in
	// log space
	for(i = 0; i < num_log_params; i++){
		mean_min_priors[logged_params[i]] = log(mean_min_priors[logged_params[i]]);
		mean_max_priors[logged_params[i]] = log(mean_max_priors[logged_params[i]]);
	}
  
	// Create a model structure, to evaluate the model at discrete values of h, corresponding to data  
	Model model_current;

	// Split (h,y) values of data
	double m_etc_vals[m_etc_len/2], m_etc_h[m_etc_len/2];
	for(i = 0; i < m_etc_len/2; i++){m_etc_vals[i] = m_etc_data[2*i+1];m_etc_h[i] = m_etc_data[2*i];}

	double p_plus_vals[p_plus_len/2], p_plus_h[p_plus_len/2];
	for(i = 0; i < p_plus_len/2; i++){p_plus_vals[i] = p_plus_data[2*i+1];p_plus_h[i] = p_plus_data[2*i];}

	double m_gly_vals[m_gly_len/2], m_gly_h[m_gly_len/2];
	for(i = 0; i < m_gly_len/2; i++){m_gly_vals[i] = m_gly_data[2*i+1];m_gly_h[i] = m_gly_data[2*i];}

	double v_vals[v_len/2], v_h[v_len/2];
	for(i = 0; i < v_len/2; i++){v_vals[i] = v_data[2*i+1];v_h[i] = v_data[2*i];}

	double g_vals[g_len/2], g_h[g_len/2];
	for(i = 0; i < g_len/2; i++){g_vals[i] = g_data[2*i+1];g_h[i] = g_data[2*i];}

	double r_max_vals[r_max_len/2], r_max_h[r_max_len/2];
	for(i = 0; i < r_max_len/2; i++){r_max_vals[i] = r_max_data[2*i+1];r_max_h[i] = r_max_data[2*i];}

	// All of the y-values for data
	double *all_vals[] = {m_etc_vals, p_plus_vals, m_gly_vals, v_vals, g_vals, r_max_vals};
	// All of the h-values for data
	double *all_hs[] = {m_etc_h, p_plus_h, m_gly_h, v_h, g_h, r_max_h};
	// Number of data points, per sub-model
	double all_vals_lens[] = {m_etc_len/2, p_plus_len/2,	m_gly_len/2, v_len/2, g_len/2, r_max_len/2};
	
	printf("Model initialised ok...\n");

	//////////////////////////////////////////////////////////////////
	// Determine regulariser scales by calculating value range
	//////////////////////////////////////////////////////////////////
	double all_reg_scales[NUM_MOD];
	
	for (i = 0; i < NUM_MOD; i++)
	{
		all_reg_scales[i] = 1.0/(Max_array(all_vals[i], all_vals_lens[i]) - Min_array(all_vals[i], all_vals_lens[i]))
;	}
	

	///////////////////////////////////////////
	// Read in Cholesky factor as 1D array
	//////////////////////////////////////////
	
	double chol[num_params*num_params];
	for (i = 0; i < num_params*num_params; i++)
	{
		fscanf(ch_p, "%lf\n", &chol[i]);
	}
	close(ch_p);


	#ifdef VERBOSE
		for (i = 0; i < num_params*num_params; i++)
		{
			printf("%lf\n", chol[i]);
		}
	#endif

	// Scale Chol
	double scale = SCALE_FAC / num_params;
	
	for(i = 0; i < num_params * num_params; i++){		
		chol[i] = chol[i] * scale;
	}

	//////////////////////////////////////////////////////////////////////////////
	// Initialise Met-Hast variables
	//////////////////////////////////////////////////////////////////////////////
	
	#ifdef VERBOSE
		double sample[num_params];
		Rand_multivar_n(sample, chol, num_params);
		Print_array(sample, num_params);
	#endif

	// Evaluate model on a discrete h-lattice, given a parameter vector
	Eval_model(h, params, &model_current, logged_params, num_params, num_log_params);
	
	// Evaluate log posterior for the model
	double lp_current, lp_proposal, ll_current, ll_proposal, lp_map, ll_map;
	
	double map[num_params]; // MAP

	lp_current = log_posterior(params, num_params, &model_current, all_vals, all_vals_lens, all_reg_scales, all_hs, mean_min_priors, mean_max_priors, &ll_current);
	lp_map = lp_current;
	ll_map = ll_current;
	for(i=0;i<num_params;i++){map[i] = params[i];}
	
	printf("Initial lp: %f\n", lp_current);

	uint64_t iter = 0; // define an unsigned 64-bit integer, as iterator
	

	double dtheta[num_params]; // Innovation
	double proposal[num_params]; // New proposal for params

	Model model_proposal; // Resultant model proposal


	// malloc for big containers, put into heap
	double** samples; // a 2D array
	samples = (double**) malloc(BUFFER_LEN * sizeof(double*)); // Create an array of pointers
	for (i = 0; i < (long)BUFFER_LEN; i++){
		samples[i] = (double*) malloc(num_params * sizeof(double)); // Give a pointer to each pointer
	}

	
	double *log_posteriors = malloc(BUFFER_LEN * sizeof(double));
	double *log_likelihoods = malloc(BUFFER_LEN * sizeof(double));

	unsigned long long_i;

	// Make a REPORT_BUFFER_LEN x 2 array, for (iter, acc_count)
	uint64_t** acc_counts;
	acc_counts = (uint64_t**) malloc(REPORT_BUFFER_LEN * sizeof(uint64_t*)); // Create array of pointers
	for(long_i = 0; long_i < REPORT_BUFFER_LEN; long_i++){
		acc_counts[long_i] = (uint64_t*) malloc(2 * sizeof(uint64_t)); // Give a pointer to each pointer
	}
	

	long acc_counter = 0; // counter for acceptance rate
	long iter_print = 0; // counter for when to save data for printing
	long iter_acc = 0; // counter for when to save acceptance rate for printing
	

	FILE *samples_fp;
	FILE *lps_fp;
	FILE *acc_fp;
	FILE *ll_fp;
	FILE *map_fp;
	samples_fp = fopen("posterior_samples.dat","w"); fclose(samples_fp);
	lps_fp = fopen("lps.dat","w"); fclose(lps_fp);
	ll_fp = fopen("ll.dat","w"); fclose(ll_fp);
	acc_fp = fopen("acc_rates.dat","w"); fclose(acc_fp);
	map_fp = fopen("map.dat","w"); fclose(map_fp);

	int j;	

	printf("Met hast initialised ok...\n");

	//////////////////////////////////////////////////////////////////////////////
	// Perform Met-Hast with correlated multidimensional Gaussian proposal
	//////////////////////////////////////////////////////////////////////////////
	
	for(iter = 0; iter < NUM_STEPS; iter++){

		// Save the parameters, if on a thinning iteration
		if((iter % (long)THINNING) == 0){			
			for(i = 0; i < num_params; i++){ samples[iter_print][i] = params[i]; } // save params
			log_posteriors[iter_print] = lp_current; // save posteriors
			log_likelihoods[iter_print] = ll_current;
			iter_print++;						
		}
		// When we reach the buffer length, write to file
		if(iter_print == (long)BUFFER_LEN){
			printf("Writing params to file...\n");		
			samples_fp = fopen("posterior_samples.dat","a");
			lps_fp = fopen("lps.dat","a"); 
			ll_fp = fopen("ll.dat","a");
			for(long_i = 0; long_i < BUFFER_LEN ; long_i++){
				for(j = 0; j < num_params; j++){
					if(j < num_params - 1){fprintf(samples_fp, "%f,", samples[long_i][j]);}
					else{fprintf(samples_fp, "%f\n", samples[long_i][j]);}
				}
				fprintf(lps_fp, "%f\n", log_posteriors[long_i]);
				fprintf(ll_fp, "%f\n", log_likelihoods[long_i]);
			}
			fclose(samples_fp);
			fclose(lps_fp);
			fclose(ll_fp);

			// Reset the buffer to start overwriting
			iter_print = 0;

			// Print MAP
			map_fp = fopen("map.dat","w");
			fprintf(map_fp, "map: ");
			for(i=0; i<num_params; i++){
				if(i == num_params - 1) fprintf(samples_fp, "%f\n", map[i]);			
				else fprintf(samples_fp, "%f,", map[i]);		
			}
			fprintf(map_fp, "lp_map: %f\n", lp_map);
			fprintf(map_fp, "ll_map: %f\n", ll_map);
			fclose(map_fp);
		}

		// Report 			
		if (iter % (long)REPORT == 0){			
			printf("%" PRIu64 ", %f, %f\n",iter, (acc_counter/((float)REPORT)), lp_current );
			acc_counts[iter_acc][0] = iter;
			acc_counts[iter_acc][1] = acc_counter;
			acc_counter = 0;
			iter_acc++;
		}
		// Send report to file once buffer is full
		if(iter_acc == (long)REPORT_BUFFER_LEN){
			printf("Writing report to file...\n");
			acc_fp = fopen("acc_rates.dat","a");
			for(long_i=0; long_i < REPORT_BUFFER_LEN; long_i++){
				fprintf(acc_fp, "%" PRIu64 ",%ld\n", acc_counts[long_i][0], acc_counts[long_i][1]);
			}
			fclose(acc_fp);

			// Reset the buffer to start overwriting
			iter_acc = 0;
		}

		//////////////////////////////////////////////////////////////////////////////
		// Do the Metropolis-Hastings step
		//////////////////////////////////////////////////////////////////////////////

		Rand_multivar_n(dtheta, chol, num_params); // Draw from proposal		

		for(i = 0; i < num_params; i++){proposal[i] = params[i] + dtheta[i];} // Perturb params
		
		// Evaluate the model for the proposal
		Eval_model(h, proposal, &model_proposal, logged_params, num_params, num_log_params);

		// Evaluate the log posterior for the proposal
		lp_proposal = log_posterior(proposal, num_params, &model_proposal, all_vals, all_vals_lens, all_reg_scales, all_hs, mean_min_priors, mean_max_priors, &ll_proposal);

		if (lp_proposal - lp_current > log(RND)){
			// Update the log posterior, likelihood and parameter values
			lp_current = lp_proposal;
			ll_current = ll_proposal;
			for(i=0; i < num_params; i++){params[i] = proposal[i];}
			acc_counter++; // increment acceptance counter
		} // Otherwise, don't update	

		// Greedy algorithm for finding the MAP
		if(lp_current > lp_map){
			lp_map = lp_current;
			ll_map = ll_current;
			for(i=0;i<num_params;i++){map[i] = params[i];}
		}

				
		
	}

	

	printf("Exiting...\n");
	return 0;
	
}




