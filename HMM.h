// Parameters internally set in expmax.c and also used in expmax_order.c :
extern float lambda0;   // Correction to singular matrix
extern int IT_MAX;      // Maximum number of iterations
extern float EPS;       // Threshold for convergence
extern int RMAX;        // Maximum number of initial conditions
extern int COMPACT;
extern int ORDER;

struct Para{
  float *tau;    // Frequency of hidden states
  float **trans; // Transition probabilities of hidden states
  float **mu;    // Mean of Gaussian distributions or max of exp distr
  float ***sig;  // Covariance matrix of Gaussian distributions
  float **scale_pos; // P(x)=(a/sp)*exp(-(x-mu)/sp) x > mu
  float **scale_neg; // P(x)=((1-a)/sn)*exp(-(mu-x)/sn) x<mu
  float *log_scale; // = -log(sp+sn)
};


// Main routines
float HMM_Clusters(int *cluster, float *Pstate,
		   int ncl, struct Para *Par, float **x, int N, int Nvar,
		   int first, char *model);
// Common auxiliary routines
int Set_parameters(struct Para *Par,
		   float **w1_ik, float ***w2_ikk1,
		   float **x_SamVar, int N, int Nvar, int ncl, int *clus);
int Initialize_clusters_2(int *clus, float **x_SamVar, int N, int Round,
			  float Thr);
int Initialize_clusters(int *clus,
			float **x_SamVar, int N, int Nvar, int ncl,
			int Round, int *cluster);

int Print_parameters(int *numcl, struct Para *Par,
		     int ncl, int Nvar, int end);
int Print_parameters_f(int *numcl, struct Para *Par,
		       int ncl, int Nvar, FILE *file_out);
float Check_convergence(struct Para *Par1, struct Para *Par2,
			int ncl, int Nvar);
void Copy_parameters(struct Para *Par1, struct Para *Par2,
		     int ncl, int Nvar);
int Gaussian_parameters(struct Para *Par,
			float **wik, float ***wikk1,
			float **x_SamVar,
			int N, int Nvar, int ncl, int iter);
int Gaussian_parameters_1(struct Para *Par,
			  float **w_ik, float **x_SamVar,
			  int N, int Nvar, int ncl);
void Empty_para(struct Para *Par, int Nvar, int ncl);
void Allocate_Para(struct Para *Par, int ncl, int Nvar, char *model);

double Assign_clusters(int *cluster, float *Pstate, int *numcl, // Output
		       float **Pxk, int N, int Nvar, int ncl);
int Process_parameters(float *logdet, float ***sighalf,  // Output
		       float ***sigvec, float  **sigval, 
		       float  **logval, int *singular,
		       int *print_sing, float ***sig, // Input
		       int Nvar, int ncl, int iter);
double Compute_logGauss(float *dx, float *x,
			float *mu, float logdet,
			float **sighalf, float **sigvec,
			float *sigval, float *logval,
			int singular, int Nvar);
float Compute_log_exp(float *x, float *mu, float *scale_pos, float *scale_neg,
		      float log_scale, int Nvar);
int Compactify_clusters(int *cluster, int *numcl, double *lik, // Output 
			float **w, float **mu, float **x_SamVar,
			int N, int Nvar, int ncl);

int Sort(int *rank, float *score, int N);
float Quad_form_choldc(float *dx, float **sighalf, int Nvar);
float Determinant(float **A, int N);
double d_Determinant(double **A, int N);
float Scalar_prod(float *u, float *v, int n);
void Empty_variance_parameters(float *logdet, float ***sighalf,
			       float ***sigvec,
			       float  **sigval, float  **logval,
			       int *singular, int *print_sing,
			       int Nvar, int ncl);


float **Transpose_matrix(float **m, int n1, int n2);
float Dist2(float *x1, float *y1, int n);
