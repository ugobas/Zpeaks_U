#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PEAKSIZE_DEF 150  // Origins located at xo +- PEAKSIZE
#define SIZEMAX 30000
extern int PEAKSIZE, MINSIZE;
#define BIDIR 100

#define NCHAR 400

/*struct fragment{
  int chr;
  long ini, end;
  float value;
  short strand;
  }; */

struct fragshort{
  int chr;
  long ini, end;
  float value;
  short strand;
  short multi;
};

struct frag{
  int chr;
  long ini, end;
  float value;
  short multi;
  short strand;
  //
  int match;
  int size;
  long x0, x1, x2;
  int lab;
  struct frag *next;
  //float w;
  //double xsum, norm;
  //float d_l, d_r;
};

// Input
void Read_name(int *nc, char **file, char *PATH, char *string, int nmax);
int Read_file(long **x, float **y, char *file);
int Read_file_3f(long **x, float **y, float **z, char *file);
int Read_file_3(long **kk, float **x1, float **x2, char *file);
int Read_bed_file(long **x, float **y, long *nn, char *file, int Nchr);
struct frag *Read_frags(long *Nfrag, int *icol,
			int *nchr, char **chr_name, int *chr_num, char *file);
struct fragshort *Read_peaks(long *Nfrag, int *icol,
			     int *Nchr, char **chr_name, int *chr_num,
			     char *file);
struct fragshort *Read_peaks_nochr(int *nori, char *file);
int Read_multi_chr(long **x, float **y, long *nn, int **strand,
		   int *Nchr, int nchr_in, char **chr_name, int *chr_num,
		   char *file);
char *Extension(char *string);
char *Remove_dir(char *s);

// Output
void Output_name(char *nameout, char *namein, char *ext);
FILE *Open_file_w(char *name);
void Print_peak(FILE *file_out, struct frag ori);
void Print_ave(double *ave, double *dev,
	       double sum1, double sum2, double norm, char *name,
	       FILE *file_out);
double Mean_size(struct frag *peaks, int N_peak);

void Plot_profile(struct frag *ori, int p,
		  long **xprof, float **yprof, long *nprof,
		  char *name_prof, char *nameout);
int Print_metaplot(double *Prof, char *name_prof, int p,
		   int NBIN, float BINSIZE, char *nameout);
void Print_peak(FILE *file_out, struct frag ori);

// Fragment operations
struct fragshort *Fragshort(struct frag *frags, int N);
int Overlap(struct frag *ori1, struct frag *ori2, int DTOL);
int Join(struct frag *ori2, struct frag *ori3, int *joined);
int Make_cluster_1(struct frag *ori, int N_ori, int *joined, float DCLUST);
void Make_cluster_4(struct frag *ori, int N_ori, int *joined, float DCLUST);
struct frag *Extract_peaks(struct frag *ori_old, int N_ori, int Ncr);
void Count_matches(int *n_match1, int *n_match2, int DTOL,
		   struct frag *ori1, int N_ori1,
		   struct frag *ori2, int N_ori2);
int Sum_match_d(int *n_matched, float *d_ave, int *same_strand,
		struct frag *ori1, int N_ori1,
		struct frag *ori2, int DTOL);
float Dist_peaks(struct frag *ori1, int N_ori1,
		 struct frag *ori2, int N_ori2);
int Remove_fragments(struct frag *ori, int N_ori, int *removed,
		     float Size_min);
void Set_ori(struct frag *orii, int i);
int Orient_peaks(struct fragshort *peaks, int N,
		 long **xp, int **strand, long *nn, int Nchr);
struct frag *Genome_regions(long *N_regions, long *Chrsize,
			      int PEAKSIZE, int Nchr);

// Coordinate operations
void Shift_coord(long **x_new, long *x_old, long n);

// Statistics of fragments
void Sizediff_stat(struct frag *ori, double *Prof_ori, double **Prof_score,
		   double *Discriminant_score, int N_prof, char **name_prof,
		   char *nameout, int PRINTALL);
/*void Compute_profile(double *Prof, double *Prof2, int Missing,
		     int NBIN, float BINSIZE,
		     struct fragshort *frags, long N_peaks, 
		     long **xprof, float **yprof, long *npro,
		     char *name);
//Moved into Metaplots.c
*/

// Statistics of vectors
void Profile_score(double *Prof_ori, double **Prof_score,
		   double *Prof_genome, double *Prof2_genome,
		   double *Discriminant_score,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_p, int N_prof,
		   struct fragshort *ori, int N_ori);
int Get_peaks(double *Discriminant_score,
	      double **Prof_score, double *Prof_ori,
	      long **x_scr, float **y_scr, long *nn, int Ncr,
	      struct frag *ori_old, int N_ori_old,
	      long ***xprof, float ***yprof, long **nprof,
	      char **name_prof, int N_prof,
	      int DCLUST, float Size_Min, int WINHALF,
	      int WINSIZE, int Stat_Comp, int DTOL, float MC,
	      char *file_old, float Thr, int PRINT);
long Genome_statistics(double *Ave, double *Dev, float **y, long *nn, int Nchr);
void Profile_peaks(double **Prof_peaks,
		   double *Ave_all, double *Dev_all, long *N_bin,
		   double *Ave_peaks, double *Sem_peaks, 
		   double *Ave_nopeak, double *Dev_nopeak,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_p, int N_prof,
		   struct fragshort *peak, int N_PEAKS);

// General
float Corr_coeff(float *slope, float *offset, float *xx, float *yy, int n);
float Corr_coeff_old(float *xx, float *yy, int n);


#define ORIMAX 1000000 // Maximum number of origins
