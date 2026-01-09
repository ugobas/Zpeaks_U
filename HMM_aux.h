// Internal parameters assigned in HMM_Clusters
//int RMAX;          // Maximum number of initial conditions
//int IT_MAX;        // Maximum number of iterations
//float EPS;         // Threshold for convergence
//float lambda0;     // Correction to singular matrix 
//double lik0;

float max(float x1, float x2);
float Rel_diff(float x, float y);
void Copy_matrix_f(float **m1, float **m2, int n1, int n2);
float Rel_diff_matrix_f(float **m1, float **m2, int n1, int n2, int sym);

float Get_lcorr(float **data, int N, int Nvar);
float AIC(float lik, int Npara, int Nsam);
float BIC(float lik, int Npara, int Nsam);
float Corr_coeff(float *slope, float *offset, float *xx, float *yy, int n);
float **Transpose_matrix(float **m, int n1, int n2);
float Determinant(float **A, int N);
double Compute_logGauss(float *dx, float *x,
			float *mu, float logdet,
			float **sighalf, float **sigvec,
			float *sigval, float *logval,
			int singular, int Nvar);
float Cosine(float *v1, float *v2, int N);
void Copy_mu(float **Mu, float **mu1, int ncl, int Nvar);
float Compare_mu(float ***Mu, int nround,  int ncl, int Nvar);

float EM_Clusters(int *clustout,
		  // Variables for storing optimal parameters:
		  float **mu3, float ***sig3, float *tau3,		  
		  float **x_VarSam,
		  int N, int ncl, int Nvar, int first);
double Assign_clusters(int *cluster, float *Pstate, int *numcl, // Output 
		       float **log_Pxk, int N, int Nvar, int ncl);
