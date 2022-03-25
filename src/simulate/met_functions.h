// Create a new variable type, called Model
typedef struct _model {
	double m_etc[NUM_H];
	double p_plus[NUM_H];
	double m_gly[NUM_H];
	double v[NUM_H];
	double g[NUM_H];
	double r_max[NUM_H];

	double all_mods[NUM_H][NUM_H] ;
} Model;



// A test function to show struct parsing syntax
void Print_value(Model *mod){
	printf("value: %f\n", (*mod).m_etc[0]); // print an element of the struct
	(*mod).m_etc[0] = (*mod).m_etc[0] + 1.0; // augment an element of the struct

}
// A test function to show array parsing syntax
void Print_value2(double *a){
	printf("value: %f\n", a[0]);
	a[0] = a[0] + 1.0;
}

void Print_array(double *a, int num_elements){
	int i;
	for (i = 0; i < num_elements; i++)
	{
		printf("%f\n", a[i]);
	}
	printf("\n");
}

void Print_model(Model *m){
	int i;
	printf("m_etc\n");
	for(i = 0; i < NUM_H; i++){
		printf("%f\n", (*m).m_etc[i]); 
	}
	printf("\n");
	printf("p_plus\n");
	for(i=0; i< NUM_H; i++){
		printf("%f\n", (*m).p_plus[i]);
	}
	printf("\n");
	printf("m_gly\n");
	for(i=0; i< NUM_H; i++){
		printf("%f\n", (*m).m_gly[i]);
	}
	printf("\n");
	printf("v\n");
	for(i=0; i< NUM_H; i++){
		printf("%f\n", (*m).v[i]);
	}
	printf("\n");
	printf("g\n");
	for(i=0; i< NUM_H; i++){
		printf("%f\n", (*m).g[i]);
	}
	printf("\n");
	printf("r_max\n");
	for(i=0; i< NUM_H; i++){
		printf("%f\n", (*m).r_max[i]);
	}
	printf("\n");
}

void Update_model_iterator(Model *m){ // Given up-to-date sub-models, update the 2D array inside the Model structure
	int i;

	for(i=0; i < NUM_H; i++){
	(*m).all_mods[0][i] = (*m).m_etc[i];
	}
	for(i=0; i < NUM_H; i++){	
	(*m).all_mods[1][i] = (*m).p_plus[i];
	}
	for(i=0; i < NUM_H; i++){	
	(*m).all_mods[2][i] = (*m).m_gly[i];
	}
	for(i=0; i < NUM_H; i++){	
	(*m).all_mods[3][i] = (*m).v[i];
	}
	for(i=0; i < NUM_H; i++){	
	(*m).all_mods[4][i] = (*m).g[i];
	}
	for(i=0; i < NUM_H; i++){	
	(*m).all_mods[5][i] = (*m).r_max[i];
	}
}

void Eval_model(double *h, double *p, Model *mod, int *logged_params, int num_params, int num_log_params){
	// Parameter ordering convention:
	// {'h_int':0, 'k_mrna':1, 'f_m':2, 'k_m':3, 'deg_b_p':4, 'ko':5, 'kg':6, 'zeta':7, 
    // 'c1':8, 'm2':9, 'k_gr':10, 'k_r':11}
	// Note: It is important NOT to modify *p in this function
	double pt[num_params];	// Create temporary container for params
	int i;
	// copy over params to local variable
	for (i = 0; i < num_params; i++){pt[i] = p[i];}  
	// exponentiate appropriate params
	for (i = 0; i < num_log_params; i++){pt[logged_params[i]] = exp(pt[logged_params[i]]);} 

	/////////////////////////////
	// Evaluate model
	/////////////////////////////
	int hin; // h iterator
	double h_mid, deg_a_m, c2, n_plus;

	for (hin = 0; hin < NUM_H; hin++){
		n_plus = (1. - h[hin]);
		//////////////////////////////////////
		// Shape ETC mRNA degradation
		//////////////////////////////////////
 		h_mid = pt[0] - log((1.0 - pt[2])/pt[2])/pt[3]; // Calculate the corresponding value of the midpoint
 		deg_a_m = 1.0 / (1.0 + exp(pt[3]*(h[hin] - h_mid)));
 		
 		//////////////////////////////////////
 		// m_gly 
 		//////////////////////////////////////
 		c2 = pt[8] - pt[9]*pt[0];
 		if (h[hin] < pt[0]){ // if h < h_int
 			(*mod).m_gly[hin] = pt[8];	
 		}
 		else{
 			(*mod).m_gly[hin] = pt[9] * h[hin] + c2;	
 		}
 		
 		//////////////////////////////////////
 		// m_etc
 		//////////////////////////////////////
 		(*mod).m_etc[hin] = pt[7] * n_plus / (deg_a_m + 1.0 / pt[1]);
 		
 		//////////////////////////////////////
 		// p_plus
 		//////////////////////////////////////
 		(*mod).p_plus[hin] = (*mod).m_etc[hin] * n_plus / pt[4]; 

 		//////////////////////////////////////
 		// v
 		//////////////////////////////////////
 		(*mod).v[hin] = pt[5] * (*mod).p_plus[hin] + pt[6] * (*mod).m_gly[hin];

 		//////////////////////////////////////
 		// g
 		//////////////////////////////////////
 		(*mod).g[hin] = pt[10] / (*mod).v[hin];

 		//////////////////////////////////////
 		// r_max
 		//////////////////////////////////////
 		(*mod).r_max[hin] = pt[11] * (*mod).p_plus[hin];

 		// Update the iterator
 		Update_model_iterator(mod);

	}
}

double Max_array(double *a, int num_elements){
	double max = a[0];
	int i;	
	for (i = 1; i < num_elements; i++)
	{
		if(a[i] > max){max = a[i];}
	}
	return max;
}

double Min_array(double *a, int num_elements){
	double min = a[0];
	int i;	
	for (i = 1; i < num_elements; i++)
	{
		if(a[i] < min){min = a[i];}
	}
	return min;
}

// Generate a single N(0,1) variable using Box-Muller transformation
double Rand_n(){ 
	double U1, U2;
	static double Z0, Z1; // static retains value between calls
	static unsigned int call = 0; // Variable to flip-flop upon each call
	//printf("call: %d\n", call);
	if(call == 1){ // spit out stored normal
		call = !call; // switch the state of call, so next time we calculate
		return Z1;
	}
	U1 = RND;
	U2 = RND;
	Z0 = sqrt( -2 * log(U1)) * cos(2*M_PI*U2);
	Z1 = sqrt( -2 * log(U1)) * sin(2*M_PI*U2);
	//printf("Generation:\n");
	//printf("%f, %f\n",Z0, Z1 );
	call = !call;
	return Z0;
}

// Draw from a multivariate Gaussian, N(0, Sigma). Requires the Cholsky factor AA^T = Sigma
// doi:10.1007/978-0-387-98144-4
void Rand_multivar_n(double *ans, double *chol, int num_rows){
	int i, j;
	double z[num_rows];
	double temp;
	for (i = 0; i < num_rows; i++)
	{
		z[i] = Rand_n(); // assign independent N(0,1)
		
	}
	// Perform matrix multiplication Az
	for(i = 0; i < num_rows; i++){
		temp = 0; // reset the value of the ith value
		for(j=0; j < num_rows; j++){
			temp += chol[num_rows*i + j] * z[j];
		}
		ans[i] = temp;
	}

}

int Get_h_index(double h_val){
	double h[] = {0.0, 0.2, 0.3, 0.5, 0.6, 0.9};
	int i;
	for (i = 0; i < NUM_H; i++)
	{
		if(fabs(h_val - h[i]) < EPSILON){
			return i;
		}		
	}
	return -1; // Value of h was not near the allowed values
}

double log_posterior(double *params, int num_params, Model *m, double *all_vals[], 
					 double *all_vals_lens, double *all_reg_scales, double *all_hs[], 
					 double *mean_min_priors, double *mean_max_priors, double *ll_current){
	int i, j, h_index;
	double y, mu, total_lp, res, ll, total_ll;

	// Initialise the total log posterior to zero
	total_lp = 0.0;
	total_ll = 0.0; // and the total log likelihood

	int num_mean_params;
	num_mean_params = num_params - NUM_MOD;

	// Get the sigmas
	double sigmas[NUM_MOD];
	for (i = 0; i < NUM_MOD; i++){sigmas[i] = params[num_mean_params + i];}
	
	// Check priors
	for(i=0; i<num_mean_params; i++) if(params[i] < mean_min_priors[i] || params[i] > mean_max_priors[i]) return -INFINITY;
	for(i=0; i<NUM_MOD; i++) if(sigmas[i] < 0) return -INFINITY;
	
	// For debugging

	for(i=0; i < NUM_MOD; i++){ // for each model
		for(j=0; j<all_vals_lens[i]; j++){ // for all the data of the submodel

			// Get the index of the h, of the data
			h_index = Get_h_index(all_hs[i][j]);
			if (h_index == -1){printf("h error!\n");exit(99);}
			
			// Get the theory value, corresponding to the h_index
			mu = (*m).all_mods[i][h_index];
			// Get data value, corresponding to h_index
			y =  all_vals[i][j];

			res = (y - mu); // residual
			
			ll = -log(sigmas[i]) - 0.5*res*res/(sigmas[i] * sigmas[i]); // log likelihood


			total_ll = total_ll + ll; // sum up log likelihood 
			total_lp = total_lp + ll; // sum up log posterior
			#ifdef VERBOSE
				printf("h: %f, th: %f, dat: %f \n",all_hs[i][j], mu, y);
			#endif 
		} 

		// Add prior. Unform priors cancel.
		total_lp = total_lp - reg_scale * all_reg_scales[i] * sigmas[i]; // log posterior
        

	}
	*ll_current = total_ll; // modify ll_current
	return total_lp;

}