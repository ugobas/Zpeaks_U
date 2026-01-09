int Cluster_peaks(struct peak *peaks, int N_peak, char *name, int makeclust);
float K_means(int k_clus, int d, int *cluster, int *nc, double **Mean,
	      float **x, int N);
float EM(int *cluster, int *nc, double **Mean, //output
	 int nclus, int d, int SIGMA, float **x, int N); //input
