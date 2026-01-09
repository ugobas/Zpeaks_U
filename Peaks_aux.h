#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PEAKSIZE_DEF 125  // Peaks located at xo +- PEAKSIZE
#define SIZEMAX 10000
#define NCHAR 400

//int PEAKSIZE;
//int ZSCORE;

struct peak{
  int chr;
  long ini, end;
  int size;
  int c;
  float y;  // Peak score
  float yout; // Peak score outside the peak
  struct peak *next;
  int ifrag;
  // int joined;
  // int removed;
  int match;
  int mult;
  int xo, x1, x2;
  //double xsum, norm;
  //float d_l, d_r;
};

void Count_matches_ZP(int *n_match1, int *n_match2, float DTOL,
		      struct peak *peak1, int N_peak1,
		      struct peak *peak2, int N_peak2);
float Dist_peaks(struct peak *peak1, int N_peak1,
		 struct peak *peak2, int N_peak2);
int Join(struct peak *peak2, struct peak *peak3);
struct peak *Extract_peaks(struct peak *peak_old, int N_peak, int Ncr,
			   long *len);

void Read_name(int *nc, char **file, char *PATH, char *string, int nmax);
int Read_file(long **x, float **y, char *file);
int Read_file_3f(long **x, float **y, float **z, char *file);
struct peak *Read_peaks(int *npeak, char *file, int shift);
int Read_wig(long *nn, long **x, float **y, short **chr,
	     char ***chr_name, char *file);
int Read_wig_2(long **x, float **y, short **chr,
	       char *file, int *nch, int Nchr);
int Read_wig_1(long *n, long *x, float *y, long nn, 
	       int *nc, long *lc, int *num_chr, char *file, int NCMAX);
int Read_chroms_wig(long *nch,char **chr_name,int *step,char *file,int Ncmax);
int Read_file_3(long **kk, float **x1, float **x2, char *file);
void Set_peak(struct peak *peak, int next);

void Output_name(char *nameout, char *namein, char *ext);
FILE *Open_file_w(char *name);

void Shift_coord(long **x_new, long *x_old, long n);
void Sizediff_stat(struct peak *peak, double *Prof_peak, double **Prof_score,
		   double *Discriminant_score, int N_prof, char **name_prof,
		   char *nameout, int PRINTALL, char *string, FILE *file_all);
int Make_cluster_1(struct peak *peak, int N_peak, float DCLUST);
void Make_cluster_4(struct peak *peak, int N_peak, float DCLUST);
int Make_cluster_check(struct peak *peak, int N_peak,
		       float DCLUST, long *xx, float *zz, long nn,
		       int check, float Thr);
void Profile_score(double *Prof_peak, double **Prof_score,
		   double *Discriminant_score,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_prof, int N_prof,
		   struct peak *peak, int N_peak);
void Plot_profile(struct peak *peak, int p,
		  long **xprof, float **yprof, long *nprof, int nchr,
		  char *name_prof, char *nameout);
int Remove_fragments(struct peak *peak, int N_peak, float Size_min);

void Print_Peaks(struct peak *peak, char *nameout);
void Print_Peaks_nomatch(struct peak *peak_old, int N_peak_old, char *nameout);
void Print_Peaks_new(struct peak *peak, char *nameout, char *nameold);
int Select_peaks_match(struct peak *peak2, struct peak *peak, int N, int match);

#define PMAX 1000000 // Maximum number of peaks
