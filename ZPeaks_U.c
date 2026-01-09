#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Peaks_aux.h"
#include "HMM_aux.h"
#include "HMM.h"
#include "allocate.h"
#include "cluster_score.h"
#include "optimization.h"
#include "cluster_peaks.h"

int SYMMETRIC=0; // Do not require symmetric maxima
int EXCLUDE_PEAK=1; // Do not average across peaks y>y_mid
int LMIN_excl=3; // Minimum number of segments for excluding a peak 3 5
int LMIN_peak=3; // Minimum number of segments for accepting a peak 3 5
int LOCAL=1; // Normalize Z scor by local or global s.dev?
int MAKECLUST=0;
float THR=2.0;
float THR_STEP=0.5;
float THR_MIN=1.5;
float EXP_MIN=1.5;
int WP_MIN=2;
float Y_MAX=0.0;  //y_thr=Y_MAX*y_max+(1-Y_MAX)*y_min
int CENTER=0;
int ZTNUM=0; // threshold on Z score depends on number of fragments
float PENALTY=0.05; // if ZTNUM thr*=(1+PENALTY*log(num[c]/num[cmin]))
float DAMP_MIN=10, DAMP_STEP=20, DAMP_MAX=100;
float EPS=0.2; // Damp=100
//int RANGE_NUC=300;
float R_MIN=10, RANGE=0;
float R_MAX=200;
int R_OUT=10000;
int N_RANGE=20;
int PRINT_CLASS=0;
int SEPARATE_NORM=1;
int SEP_SMOOTH=1;
int SCORE=1;
char SCORE_R='N'; // Use non-monotonic score to optimize RANGE
float OUTL=3;
int OPEN_FILE=1;

//int ZRANGE=200;
//int WGT_INI=8, WGT_END=20, WGT_STEP=4;
//int KMEAN=0; 
// Outdated:
int DCLUST=0, DTOL=0;
int SIZE_MIN=0;
int SEP_KMEAN=0;

int RANDOM=1;
int OMIT_MITO=1; // Omit mitochondrion and chloroplast?
int NORM_PEAK_SCORE=1; // Normalize the score of each peak dividing by length?

// Classes of maxima
int NC=6, NC1; //NC-1;
//float ymax_c[9]={-1,-0.5,0.0,0.5,1.0,100}; // Classes of peak score
//float ymax_c[6]={-10,-1,-0.25,0.25,1.0}; // Classes of peak score
//float ymax_c[7]={-10,-0.75,0.0,0.75,10}; // Classes of peak score
//int lab_cl[6]={-3,-2,-1,1,2,3}; // Label of classes
float ymax_c[7]={-100,-0.7,-0.1,0.5,1.5,100}; // Classes of peak score
//float ymax_c[7]={-10,-0.6,-0.2,0.2,0.6,10}; // Classes of peak score
int lab_cl[7]={-3,-2,-1,0,1,2,3}; // Label of classes
float *thr_WP; int N_WP=4, N_WP1;
// We require lower Z score if WP is larger

void Get_threshold(float Thr, float *zt, float *y_scr, long nn, int NC);
int Get_maxima(int *lpeak, float *yout_peak, double *S, float *y1,
	       float *y_scr, float *y_exper, long nn,
	       float ythr, int wp, float *Z_thr, float Sdy,
	       float Sqr_w, int get_y);
extern int Read_multi_chr(long **x, float **y, long *nn, int **strand,
			  int *Nchr, int Nchr_in, char **chr_name, int *chr_num,
			  char *file);
extern float Corr_coeff(float *slope, float *offset,
			float *x, float *y, int n);
extern struct peak *Read_peaks_ZP(int *npeak, char *file, int shift);
extern void Print_peak_score(struct peak *peak, double **Peak_score,
			     int N_prof, char **name_prof, char *nameout);
extern void Profile_score_ZP(double **Peak_score, FILE *file_out, 
			     long ***xprof, float ***yprof, long **nprof,
			     int *nchr_prof,
			     double **ave_prof,double **dev_prof,int N_prof,  
			     struct peak *peak, int N_PEAK);

extern struct peak *Extract_peaks_ZP(struct peak *peak_old, int N_peak,
				     int Nchr, long *len_chr);
extern void Count_matches_ZP(int *n_match1, int *n_match2, float DTOL,
			     struct peak *peak1, int N_peak1,
			     struct peak *peak2, int N_peak2);
static void Ave_prof(double **ave_prof, double **dev_prof, int N_prof, int Nchr,
		     float ***yprof, long **nprof);

void AT_norm_control(char *file_c, long *x_coord, float *y_contr,
		     char **chromosome, char **chr_name, long *Lchr,
		     long *nnch, int Nchr);
long Read_sequence(char **seq, char *file_chr);
void Get_chr_name(char **chr_lab, char *string);
static char Maiuscule(char s);

extern int ZSCORE, PEAKSIZE, IT_MAX, RMAX;
extern float EPS, lambda0; 
extern double lik0;

#define CODE "ZPeaks_U"
#define NCHRMAX 80   // Maximum number of chromosomes
#define NPROF 100 // Maximum number of profiles to read
#define VERBOSE 1

char NAMEOUT[80]="Peaks_U";

// Control variables
int STAT_COMP=1;    // Make statistics of comparison?
int PRSCORE=0;      // Print scores in file?

// Global parameters
char *file_e[NCHRMAX], *file_c[NCHRMAX];
int SIZE_INI=200; // For removing small fragments
int SIZE_END=300;
int SIZE_STEP=50;
float Thr;
int WIN=50;

// HMM parameters
int ncl=2;
int *numclus, *numdom;
char namelog[200];
FILE *filelog;

// data
int Nchr, CONTR=0;
long nnch[NCHRMAX], nn, *len_chr ;
long inichr[NCHRMAX], endchr[NCHRMAX];
long *x_coord;
float *y_contr=NULL, *y_exper=NULL, *y_diff=NULL;
char *chr_name[NCHRMAX];
int chr_num[NCHRMAX];
int *Omit_chr=NULL;

// Comparison
struct peak *peak_ref=NULL;
char *file_ref=NULL;
int N_peak_ref=0;

// Profiles
int N_prof;
int nchr_prof[NPROF];
char  **file_prof[NPROF], *name_prof[NPROF];

int Get_peaks(int *N, double *S, int *lpeak, float *y_scr, int *WIN,
	      float Damp, float Thr, float l_EPS, int S_MIN);
//int K_means_chr(int *N, double *S, double **Mean, int *cluster, float *y_scr,
//		int *WIN, float Damp, float l_EPS, int NC);
int K_means_old(int *Nc, double *Mean, double *Sd, float *ythr_clus,
		int *cluster, float *y_scr, int nk, int NC);
struct peak *Peaks2Frags(int *N_peak, int *lpeak, float *yout_peak,
			 float *y_scr, long nn, int c);
struct peak *CopyPeaks(int *N_peak, struct peak *peaks_all,
		       int N_peak_all, int CLUS);
struct peak *Center_peaks(int *N_peak, struct peak *peaks,
			  float *y_scr, long *x_coord, long nn);

float Smooth(float *yb, float *y, long N, float *Sqr_w,
	     float DAMP, int *WIN, float l_EPS);
void Smooth_score(float *y_scr, long nn, float Damp, int *WIN, float l_EPS);
int Get_ave_sd_outliers(double *Y1, double *Y2, float *y, int nk, float OUTL);

float Peak_likelihood(int *lpeak, float *y, int nn, int ncl);
int Get_frags(int *lpeak, double *S, float Thr, float *y_scr, int SIZEMIN);
float Get_likelihood(float Thr, int *lpeak, float *y, long nn, int SIZE_MIN);
static void Count_peaks(int *numclus, int *numdom, int *lpeak, long *endchr,
			int N,int ncl);
static int Set_control(float *yc, long nn, float Mean);
static int Get_input(int *CONTR, char **file_c, char **file_e,
		     // control and experiment
		     char **file_ref, // reference peaks
		     char **file_genome, char **chr_lab, int *NORM_AT,// genome
		     char ***file_prof, int *N_prof, // profiles
		     char **name_prof, int *nchr_prof,
		     // parameters: 
		     int *NC, float *THR, float *THR_STEP, float *THR_MIN,
		     float *EXP_MIN, float *Y_MAX,
		     int *LOCAL, int *ZTNUM, 
		     int *SEPARATE_NORM, int *SEP_SMOOTH,
		     int *WP_MIN, int *SCORE, char *SCORE_R, float *OUTL,
		     float *DAMP_MIN, float *DAMP_STEP, float *DAMP_MAX,
		     float *EPS, int *LMIN_excl, int *LMIN_peak,
		     float *R_MIN, float *R_MAX, int *N_RANGE,
		     int *PRINT_CLASS,int *PRSCORE, char *NAMEOUT,
		     //int *KMEAN, int *SEP_KMEAN, float *THR_MIN,
		     //int *SIZE_MIN, int *DCLUST, int *DTOL, 
		     //int *ZRANGE, int *WGT_INI, int *WGT_END, int *WGT_STEP,
		     int argc, char **argv);
static void help(char *name);

static void Copy_prof(long *nprof, long **xprof, float **yprof,
		      long *x, float *y, long *nnch, int Nchr);
void Comp_Z_score(float *Z_scr, long nn, float *ye_w, float *yc_w, float sd);
static void Rescale_counts(float *yy, long nn, float scale);
static float Normalize_counts(float *y, long n, float AVE);
static float Mean_counts(float *sd, float *yy, long nn);
static void Print_score(long *x, float *y, long *nnch, long nn,
			char **chr_name, char *NAMEOUT);
static int Read_contr_exp(char **chr_name, long *nn, long *nnch,
			  long **x_coord, float **y_contr, float **y_exper,
			  float *r,
			  char **file_e, int nexp, char **file_c, int ncontr);
static int Count_chromosomes(char **chr_name, long *lchr,  int *step,
			     char **file, int nf);
static void Read_chromosomes(long *x, float *y, long *lch, long nn,
			     int *num_chr, char **file, int nf);
extern char *Extension(char *string);
void Examine_peaks(struct peak *peaks, int N_peak,
		   struct peak *peak_ref, int N_peak_ref, int Nchr,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_prof, char **name_prof,
		   double **ave_prof, double **dev_prof, int N_prof,
		   char *NAMEOUT,
		   //int WINDOW, int DCLUST, int SIZE_MIN,
		   int WINDOW, float Thr, char *Param,
		   char *what, char *file_ref,
		   float *y1, int clus);
int Close_peak(struct peak **peak, long x_coord, float y, float yout,
	       int *numch, double *sizech);
/*void Compute_scores(int *lpeak, short *chr, float **x_VarSam,
		    int Nsam, int Nvar, int ncl, float lik,
		    struct Para Par, FILE *filelog);*/

//float *Mean_control, *Mean_exp, *sd_control, *sd_exp;
int BOX1; // Size of the wig files

int main(int argc, char **argv)
{
  ZSCORE=1; // Print Z score of properties (1) or raw properties (0)?
  //strcpy(MODEL,"E"); // E=exponential G=Gaussian
  int p, k, j; long i;

  // Read input
  for(p=0; p<NPROF; p++){
    file_prof[p]=malloc(NCHRMAX*sizeof(char *));
    for(k=0; k<NCHRMAX; k++)file_prof[p][k]=NULL;
    nchr_prof[p]=0;
  }
  char *file_genome[NCHRMAX]; // Chromosome sequences
  char *chr_lab[NCHRMAX];  // Character labelling chromosome
  int NORM_AT=0; // Normalize control?
  for(k=0; k<NCHRMAX; k++){
    file_genome[k]=NULL; chr_lab[k]=NULL;
  }

  int ncontr,nexp=
    Get_input(&ncontr, file_c, file_e, &file_ref,
	      file_genome, chr_lab, &NORM_AT,
	      file_prof, &N_prof, name_prof, nchr_prof,
	      &NC, &THR, &THR_STEP, &THR_MIN, &EXP_MIN, &Y_MAX, &LOCAL, &ZTNUM,
	      &SEPARATE_NORM, &SEP_SMOOTH,
	      &WP_MIN, &SCORE, &SCORE_R, &OUTL,
	      &DAMP_MIN, &DAMP_STEP, &DAMP_MAX, &EPS, &LMIN_excl, &LMIN_peak,
	      &R_MIN, &R_MAX, &N_RANGE, &PRINT_CLASS, &PRSCORE, NAMEOUT,
	      // &KMEAN, &SEP_KMEAN, &THR_MIN, &SIZE_MIN, &DCLUST, &DTOL
	      // &ZRANGE, &WGT_INI, &WGT_END, &WGT_STEP,
	      argc, argv);

  if(THR<THR_MIN){
    printf("WARNING, THR= %.3f is lower than minimum allowed value, "
	   "setting it to %.3f\n");
    THR=THR_MIN;
  }

  // Read experiment and control
  float corr_exp_contr=0;
  Nchr=Read_contr_exp(chr_name, &nn, nnch, &x_coord, &y_contr, &y_exper, 
		      &corr_exp_contr, file_e, nexp, file_c, ncontr);

  if(OMIT_MITO){
    printf("Eliminating mitochondria and chloroplast, if any\n");
    Omit_chr=malloc(Nchr*sizeof(int));
    int no=0, imax=0; long nk=0;
    for(i=0; i<Nchr; i++){
      if((strncmp(chr_name[i], "mito", 4)==0)||
	 (strncmp(chr_name[i], "chloro", 6)==0)){
	Omit_chr[i]=1; no++;
      }else{
	Omit_chr[i]=0;
	imax=i+1; nk+=nnch[i];
      }
    }
    printf("%d chromosomes eliminated, %d remaining\n", no, Nchr-no);
    if((imax+no)==Nchr){Nchr=imax; nn=nk;}
    else{
      printf("WARNING, the eliminated chromosomes are not the last ones\n");
    }
  }

  // Read chromosomes
  char *chromosome[Nchr]; long Lchr[Nchr], Lmax=0;
  for(k=0; k<Nchr; k++){chromosome[k]=NULL; Lchr[k]=0;}
  if(NORM_AT && file_genome[0] && y_contr){
    printf("Normalizing control with T score read in %s\n",file_genome[0]);
    { // Check if control name contains AT_norm
      char *s=file_c[0];
      while(*s!='\0'){
	if(strncmp(s, "ATnorm", 6)==0 || strncmp(s, "AT_norm", 7)==0){
	  printf("ERROR, you try to normalize already normalized control %s\n",
		 file_c[0]); exit(8);
	}
	s++;
      }
    }

    for(k=0; k<Nchr; k++){
      printf("Reading chromosome %s in %s\n", chr_name[k], file_genome[k]);
      Lchr[k]=Read_sequence(&(chromosome[k]), file_genome[k]);
      if(Lchr[k]>Lmax)Lmax=Lchr[k];
      printf("L= %ld\n", Lchr[k]);
    }
    AT_norm_control(file_c[0], x_coord, y_contr, chromosome, chr_name,
		    Lchr, nnch, Nchr);
  }

  // Compare genome and control
  len_chr=malloc(Nchr*sizeof(long));
  long ini=0;
  for(k=0; k<Nchr; k++){
    inichr[k]=ini; endchr[k]=ini+nnch[k]; ini=endchr[k];
    len_chr[k]=x_coord[endchr[k]-1];
    if(Lchr[k]){
      if(abs(Lchr[k]-len_chr[k])>200){
	printf("WARNING, different chromosome length in control (%ld) and "
	       "in genome (%ld)\n", len_chr[k], Lchr[k]);
      }
      if(Lchr[k]>len_chr[k]){len_chr[k]=Lchr[k];}
    }
    printf("chr %s (%d): ini= %ld length= %ld bp nn=%ld\n",
	   chr_name[k], k+1, x_coord[inichr[k]], len_chr[k], nnch[k]);
  }

  BOX1=x_coord[1]-x_coord[0];
  if(y_contr)CONTR=1; 
  printf("First bin: %d %d\n", x_coord[0], x_coord[1]);
  printf("Box size= %d nucleotides %d chromosomes\n", BOX1, Nchr);
  //RANGE_NUC/(float)BOX1; 
  if(R_MIN<5)R_MIN=5;
  if(RANGE<R_MIN)RANGE=R_MIN;
  R_OUT/=BOX1;
  if(R_MIN>R_MAX){
    int tmp=R_MIN; R_MIN=R_MAX; R_MAX=tmp;
  }
  printf("RANGE= %.0f bins, %d bp BOX1= %d\n", RANGE, RANGE*BOX1, BOX1);
  printf("%d chromosomes\n", Nchr);

  // Read profiles, if any
  N_prof++; // Profile N_prof-1 is Peakscore and it is set later
  long  **nprof=malloc(N_prof*sizeof(long *));
  long  ***xprof=malloc(N_prof*sizeof(long **));
  float ***yprof=malloc(N_prof*sizeof(float **));
  for(i=0; i<Nchr; i++)chr_num[i]=i+1;
  int ip=0, skipped=0;
  for(p=0; p<(N_prof-1); p++){
    printf("Reading profile %s %d files\n", name_prof[p], nchr_prof[p]);
    nprof[ip]=malloc(Nchr*sizeof(long));
    xprof[ip]=malloc(Nchr*sizeof(long *));
    yprof[ip]=malloc(Nchr*sizeof(float *));
    if(nchr_prof[p]==Nchr){
      for(k=0; k<Nchr; k++){
	nprof[ip][k]=Read_file(&xprof[ip][k], &yprof[ip][k], file_prof[ip][k]);
	printf("Reading %d lines %s\n", nprof[ip][k], file_prof[p][k]);
      }
    }else if(nchr_prof[p]==1){
      int nchr=Nchr; 
      Read_multi_chr(xprof[ip], yprof[ip], nprof[ip], NULL, &nchr, Nchr,
		     chr_name, chr_num, file_prof[p][0]);
    }else{
      printf("WARNING, profile %s in %d files, neither %d nor 1 skipping it\n",
	     name_prof[p], nchr_prof[p], Nchr);
      free(nprof[ip]); free(xprof[ip]); free(yprof[ip]); ip--;
      skipped++;
    }
    ip++;
  }
  N_prof-=skipped;
  printf("%d profiles stored, %d skipped\n", N_prof, skipped);

  // Read reference peaks, if any
  char namefound[200];
  char **lab_p=NULL;
  if(file_ref!=NULL){
    printf("Reading %s\n", file_ref);
    peak_ref=Read_peaks_ZP(&N_peak_ref, file_ref, -1); // Chr=chr-1
    printf("%d reference peaks\n", N_peak_ref);
    sprintf(namefound, "%s_found.bed", file_ref);
  }

  // Normalize to one read per nucleotide
  float AVE=1.0;
  // Normalize both control and experiment per each chromosome
  if(CONTR==0){
    y_contr=malloc(nn*sizeof(float));
    Set_control(y_contr, nn, AVE);
    printf("Using uniform control\n");
  }else{
    printf("Using control\n");
  }
  printf("%d chromosomes\n", Nchr);

  double YC1=0, YC2=0;
  {float *yc=y_contr; for(i=0; i<nn; i++){YC1+=*yc; YC2+=(*yc)*(*yc); yc++;}}
  YC1/=nn; YC2=YC2-nn*YC1*YC1; YC2=sqrt(YC2/(nn-1));
  printf("Mean control: %.3g s.dev.: %.3g\n", YC1, YC2);

  if(SEPARATE_NORM){
    long n0=0;
    for(i=0; i<Nchr; i++){
      long ni=nnch[i];
      Normalize_counts(y_exper+n0, ni, AVE);
      if(ncontr)Normalize_counts(y_contr+n0, ni, AVE);
      n0+=ni;
    }
  }else{
    Normalize_counts(y_exper, nn, AVE);
    if(ncontr)Normalize_counts(y_contr, nn, AVE);  
  }

  YC1=0; YC2=0;
  {float *yc=y_contr; for(i=0; i<nn; i++){YC1+=*yc; YC2+=(*yc)*(*yc); yc++;}}
  YC1/=nn; YC2=YC2-nn*YC1*YC1; YC2=sqrt(YC2/(nn-1));
  printf("Mean control: %.3g s.dev.: %.3g\n", YC1, YC2);
  if(isnan(YC1) || isnan(YC2))exit(8);
  printf("Exclude bin from ave if y>y_min(peak)+ %.3f*(y_max-y_min)(peak)\n",
	 Y_MAX);

  y_diff=malloc(nn*sizeof(float));
  for(i=0; i<nn; i++)y_diff[i]=(y_exper[i]-y_contr[i])/YC2;

  // Set Peakscore
  p=N_prof-1;
  nchr_prof[p]=Nchr;
  name_prof[p]=malloc(100*sizeof(char));
  strcpy(name_prof[p], "Peakscore");
  nprof[p]=malloc(Nchr*sizeof(long));
  xprof[p]=malloc(Nchr*sizeof(long *));
  yprof[p]=malloc(Nchr*sizeof(float *));
  for(k=0; k<Nchr; k++){
    nprof[p][k]=nnch[k];
    xprof[p][k]=malloc(nnch[k]*sizeof(long));
    yprof[p][k]=malloc(nnch[k]*sizeof(float));
  }

  NC1=NC-1; //NC1=2 is the label of peaks
  sprintf(namelog,"%s.log", NAMEOUT);
  filelog=fopen(namelog, "w");
  fprintf(filelog, "Corr(exp,contr)= %.3f\n", corr_exp_contr);

  int *lpeak=malloc(nn*sizeof(int));
  float *yout_peak=malloc(nn*sizeof(float));
  float *y_scr=malloc(nn*sizeof(float));
  numdom=malloc(NC*sizeof(int));
  numclus=malloc(NC*sizeof(int));
  float *y1=malloc((NC+1)*sizeof(float));

  float ythr=EXP_MIN*(AVE/BOX1); // Minimum value of y_exper for peak

  thr_WP=malloc(N_WP*sizeof(float)); N_WP1=N_WP-1; Thr=THR;
  for(int i=0; i<N_WP; i++){
    if(Thr < THR_MIN){Thr=THR_MIN;}
    thr_WP[i]=Thr; Thr-=THR_STEP;
    printf("WP= %d thr= %.2g\n", WP_MIN+i, thr_WP[i]);
    fprintf(filelog, "#WP= %d thr= %.2g\n", WP_MIN+i, thr_WP[i]);
  } 
  Thr=THR;
  float zt[NC]; for(int c=0; c<NC; c++)zt[c]=1;

  // Find optimal damping factor and EPS
  int N_EPS=20;
  float l_EPS_MIN=0.7, l_EPS_MAX=10;
  float l_EPS_INI=-log(EPS);
  if(l_EPS_INI<l_EPS_MIN){l_EPS_INI=l_EPS_MIN;}
  if(l_EPS_INI>l_EPS_MAX){l_EPS_INI=l_EPS_MAX;}
  float lE_STEP=log(1.5), l_EPS=l_EPS_INI-lE_STEP;
  float eps[4], Sc_e[4];
  // Damp
  int N_Damp=20;
  float Damp, Damp_ini=DAMP_MAX*0.5;
  float D[4], Sc_d[4]; //int wgt_o[4], wp_o[4], range_o[4];
  // Range
  float R_STEP=(R_MAX-R_MIN)/4; //R_MAX*0.1; if(R_STEP<5)R_STEP=5;
  float RANGE_INI=R_MAX/2-R_STEP; //2*R_MIN; //0.5*R_MAX; 
  if(N_RANGE==1)RANGE_INI=R_MAX;
  if(RANGE_INI<R_MIN)RANGE_INI=R_MIN;
  float rr[4], Sc_r[4], N_r[4], S_r[4];

  printf("Damp= %.3g-%.3g %.3g\n\n", DAMP_MIN, DAMP_MAX, DAMP_STEP);

  float S_opt=-100000, Score_opt=0;
  int N_opt=0, W_opt=0;
  float range_opt=0, lE_opt=0, Damp_opt=0;

  int NPARA_MAX=1000, Npara=0;
  float *Para[NPARA_MAX], thr_para=0.001;
  for(int i=0; i<NPARA_MAX; i++)Para[i]=malloc(3*sizeof(float));
  int failure_r=0, kr, ir, n_failure=3;
  RANGE=RANGE_INI-R_STEP; //R_MAX-3*R_STEP;
  for(kr=0; kr<N_RANGE; kr++){
    if(kr<2){
      ir=kr;
      RANGE+=R_STEP;
    }else if(kr==2){
      ir=2;
      if(SCORE_R=='X'){
	if(S_r[1]<S_r[0]){ // Decreasing function, use S to optimize range
	  SCORE_R='S';
	  printf("## The score is a decreasing function of range\n" 
		 "## Changing criterion to optimize range from N to score\n");
	  Sc_r[0]=S_r[0];
	  Sc_r[1]=S_r[1];
	}else if(N_r[1]<N_r[0]){ // Decreasing function, use N to optimize range
	  SCORE_R='N';
	  printf("## N is a decreasing function of range\n" 
		 "## Using criterion N to optimize range\n");
	  Sc_r[0]=N_r[0];
	  Sc_r[1]=N_r[1];
	}else{
	  printf("## Both score and N are increasing function of range\n");
	  printf("## Fix range by maximizing Score/");
	  float tmp0=S_r[0]/sqrt(N_r[0]), tmp1=S_r[1]/sqrt(N_r[1]);
	  if(tmp1<tmp0){
	    SCORE_R='D';
	    printf("sqrt(N)\n");
	    Sc_r[0]=tmp0; Sc_r[1]=tmp1;
	  }else{
	    SCORE_R='s';
	    printf("N\n");
	    Sc_r[0]=S_r[0]/N_r[0]; Sc_r[1]=S_r[1]/N_r[1];
	  } 
	}
      }
      printf("## S[range=%.0f]: %.4g\n", rr[0], Sc_r[0]);
      printf("## S[range=%.0f]: %.4g\n", rr[1], Sc_r[1]);
      if(Sc_r[1]>Sc_r[0]){
	RANGE=rr[1]+R_STEP;
	S_opt=Sc_r[1]; N_opt=N_r[1];
      }else{ // Decreasing function 
	//RANGE=R_MIN;
	RANGE=rr[0]-R_STEP; 
	S_opt=Sc_r[0]; N_opt=N_r[0];
      }

    }else{
      ir=3;
      float R_OLD=RANGE;
      RANGE=Find_max_quad(rr, Sc_r, R_MIN, R_MAX);
      if(isnan(RANGE)){printf("Range is nan, breaking\n"); break;}
      if(fabs(RANGE-R_OLD)<R_OLD*0.002){
	printf("RANGE almost repeated: %.0f %.0f %.0f < =.f, breaking\n",
	       RANGE, R_OLD, fabs(RANGE-R_OLD), R_OLD*0.001);
	break;
      }
      printf("New value of range: %.0f in %.0f %.0f\n", RANGE, R_MIN, R_MAX);
    }
    if(RANGE<R_MIN)RANGE=R_MIN;
    if(RANGE>R_MAX)RANGE=R_MAX;
    if(ir==3 && (RANGE==rr[0] || RANGE==rr[1] || RANGE==rr[2])){
      printf("Range= %g repeated, breaking\n", RANGE);
      printf("Previous values= %g %g %g\n", rr[0], rr[1], rr[2]);
      break;
    }
    //printf("l_EPS= %.3g Damp=%.3g range=%d\n", l_EPS, Damp, RANGE);
    
    float Sc_opt_r=-10000, Sc_opt_e=-10000;
    int N_opt_r=0, W_opt_r=0;
    float lE_opt_r=0, Damp_opt_r=0;

    // Loop over eps
    int failure_e=0, ke, ie=0;
    if(kr==0){l_EPS=l_EPS_INI;}
    else{
      l_EPS= lE_opt-lE_STEP; if(l_EPS>0.5*l_EPS_MAX)l_EPS=0.5*l_EPS_MAX;
    }
    for(ke=0; ke<N_EPS; ke++){
      if(ke<2){
	ie=ke; l_EPS+=lE_STEP;
      }else if(ke==2){
	ie=2;
	if(Sc_e[0]>Sc_e[1]){l_EPS=eps[0]-lE_STEP;}
	else{l_EPS=eps[1]+lE_STEP;}
      }else{
	ie=3;
	l_EPS=Find_max_quad(eps, Sc_e, l_EPS_MIN, l_EPS_MAX);
	if(isnan(l_EPS)){printf("EPS is nan, breaking\n"); break;}
      }
      if(l_EPS>l_EPS_MAX)l_EPS=l_EPS_MAX;
      if(l_EPS<l_EPS_MIN)l_EPS=l_EPS_MIN;
      if(ie==3 && (l_EPS==eps[0] || l_EPS==eps[1] || l_EPS==eps[2])){
	printf("l_EPS= %.3g repeated, breaking\n", l_EPS);
	break;
      }

      int failure_d=0, kd, id;
      float Sc_opt_d=-10000;
      if(kr==0){Damp=DAMP_MIN;}
      else{
	Damp=Damp_opt-DAMP_STEP; if(Damp>Damp_ini)Damp=Damp_ini;
	//Damp=DAMP_MIN;
      }
      for(kd=0; kd<N_Damp; kd++){
	if(kd<2){
	  id=kd; Damp+=DAMP_STEP;
	}else if(kd==2){
	  id=2;
	  if(Sc_d[0]>Sc_d[1]){Damp=D[0]-DAMP_STEP;}
	  else{Damp=D[1]+DAMP_STEP;}
	}else{
	  id=3;
	  Damp=Find_max_quad(D, Sc_d, 1, DAMP_MAX);
	  if(isnan(Damp)){printf("Damp is nan, breaking\n"); break;}
	}
	if(Damp<1)Damp=1;
	if(Damp>DAMP_MAX)Damp=DAMP_MAX;
	if(id==3 && (Damp==D[0] || Damp==D[1] || Damp==D[2])){
	  printf("Damp=%.3g repeated, breaking\n", Damp); 
	  break;
	}

	float Sqr_w, Sdy=Smooth(y_scr, y_diff, nn, &Sqr_w, Damp, &WIN, l_EPS);
	if(ZTNUM)Get_threshold(Thr, zt, y_scr, nn, NC);

	double S=0; int N=0;
	for(int i=0; i<Npara; i++){
	  if(fabs(Para[i][0]-RANGE)<fabs(Para[i][0]*thr_para) &&
	     fabs(Para[i][1]-WIN)<fabs(Para[i][1]*thr_para) &&
	     fabs(Para[i][2]-Damp) <fabs(Para[i][2]*thr_para)){
	    printf("Parameters are repeated, leaving\n");
	    goto next_d;
	  }
	}

	N=Get_maxima(lpeak,yout_peak, &S,y1,y_scr,y_exper,
		     nn,ythr,WP_MIN,zt,Sdy,Sqr_w,0);
	printf("%d peaks RANGE= %.0f WIN=%d Damp=%.3g score= %.0f\n",
	       N, RANGE, WIN, Damp, S);
	Para[Npara][0]=RANGE;
	Para[Npara][1]=WIN;
	Para[Npara][2]=Damp;
	if(Npara<(NPARA_MAX-1)){Npara++;}

	Sc_d[id]=S;
	D[id]=Damp;
	if((Sc_d[id] > Sc_opt_d) || (id==0 && ir==0)){
	  if(Sc_d[id] > Sc_opt_d*1.001){failure_d=0;}
	  else{failure_d++;}
	  Sc_opt_d=Sc_d[id];
	  if(Sc_opt_d > Sc_opt_r){
	    Sc_opt_r=Sc_opt_d;
	    N_opt_r=N;
	    Damp_opt_r=Damp;
	    lE_opt_r=l_EPS;
	    W_opt_r=WIN;
	  }
	}else if(id==3){
	  failure_d++;
	}
	if(failure_d==n_failure)break;
	if(kd>=3){
	  if(Rearrange_points(D, Sc_d, D[3], Sc_d[3])<0)break;
	}
      next_d: continue;
      } // End loop of Damping
      if(kd==N_Damp && failure_d==0)
	printf("WARNING, optimal damping not found in %d steps\n",kd);
      
      eps[ie]=l_EPS;
      Sc_e[ie]=Sc_opt_d; 
      if((Sc_e[ie]>Sc_opt_e)||(ie==0)){
	if(Sc_e[ie]>Sc_opt_e*1.001){failure_e=0;}
	else{failure_e++;}
	Sc_opt_e=Sc_e[ie];
      }else if(ie==3){
	failure_e++;
      }
      if(failure_e==n_failure)break;
      if(ke>=3){
	if(Rearrange_points(eps, Sc_e, eps[3], Sc_e[3])<0)break;
      }

    } // End loop of l_EPS
    if(ke==N_EPS && failure_e==0){
      printf("WARNING, optimal l_EPS not found in %d steps\n",ke);
    }
    rr[ir]=RANGE;
    if(SCORE_R=='N'){
      Sc_r[ir]=N_opt_r;  // Determine RANGE by optimizing number of peaks
    }else if(SCORE_R=='S' || SCORE_R=='X'){ 
      Sc_r[ir]=Sc_opt_r;
    }else if(SCORE_R=='s'){
      Sc_r[ir]=Sc_opt_r/N_opt_r;
    }else if(SCORE_R=='D'){
      Sc_r[ir]=Sc_opt_r/sqrt(N_opt_r);
    }
    N_r[ir]=N_opt_r;
    S_r[ir]=Sc_opt_r;

    printf("### Best: Score= %.0f (%d) Np= %d "
	   "W=%d Damp=%.0f l_EPS=%.3g range= %.0f "
	   "wp= %d Z_thr=%.2g local= %d S/sqrt(Np)= %.4g Score[range]=%.4g\n\n",
	   Sc_opt_r, SCORE, N_opt_r, (2*W_opt_r+1)*BOX1, Damp_opt_r, lE_opt_r,
	   RANGE, WP_MIN,THR,LOCAL,Sc_opt_r/sqrt(N_opt_r),Sc_r[ir]);


    if(kr==0 || Sc_r[ir]>S_opt){
      if(Sc_r[ir]>S_opt*1.001){failure_r=0;}
      else{failure_r++;}
      N_opt=N_opt_r; 
      Score_opt=Sc_opt_r;
      W_opt=W_opt_r;
      lE_opt=lE_opt_r;
      Damp_opt=Damp_opt_r;
      range_opt=RANGE;
      S_opt=Sc_r[ir];
    }else if(ir==3){
      failure_r++;
    }
    if(failure_r==n_failure)break;
    if(ir>=3){
      if(Rearrange_points(rr, Sc_r, rr[3], Sc_r[3])<0)break;
    }
  } // end loop of range

  if(kr==N_RANGE && N_RANGE>4 && failure_r==0){
    printf("WARNING, optimal range not found in %d steps\n",kr);
  }

  Damp=Damp_opt;
  l_EPS=lE_opt;
  RANGE=range_opt;
  printf("\nBest: Score= %.0f (%d) Np= %d "
	 "W=%d Damp=%.0f l_EPS=%.3g range= %.0f "
	 "wp= %d ythr=%.2g*(AVE/BOX1) Z_thr=%.2g local= %d ztnum= %d "
	 "S/sqrt(N)=%.4g\n",
	 Score_opt, SCORE, N_opt, (2*W_opt+1)*BOX1, Damp_opt, lE_opt,range_opt,
	 WP_MIN,EXP_MIN,THR,LOCAL,ZTNUM,Score_opt/sqrt(N_opt));

  // Record score in wig file and in profile N_prof-1
  float Sqr_w, Sdy=Smooth(y_scr, y_diff, nn, &Sqr_w, Damp, &WIN, l_EPS);
  if(PRSCORE)Print_score(x_coord, y_scr, nnch, nn, chr_name, NAMEOUT);
  // Copy profile score
  p=N_prof-1;
  Copy_prof(nprof[p], xprof[p], yprof[p], x_coord, y_scr, nnch, Nchr);

  double **ave_prof=Allocate_mat2_d(N_prof, Nchr);
  double **dev_prof=Allocate_mat2_d(N_prof, Nchr);
  Ave_prof(ave_prof, dev_prof, N_prof, Nchr, yprof, nprof);

  if(ZTNUM)Get_threshold(Thr, zt, y_scr, nn, NC);
  double S=0;
  int N=Get_maxima(lpeak,yout_peak,&S,y1,y_scr,y_exper,
		   nn,ythr,WP_MIN,zt,Sdy,Sqr_w,1);
  printf("%d peaks score= %.0f win= %d\n", N, S, WIN);

  // Transform into fragments
  int N_peak[NC+1]; struct peak *peaks[NC+1];
  peaks[0]=NULL; N_peak[NC]=0;
  Count_peaks(numclus, N_peak, lpeak, endchr, nn, NC);
  for(int i=1; i<NC; i++){
    printf("lpeak %d: %d peaks < %.2g\n", i, N_peak[i], ymax_c[i]);
    if(i)N_peak[NC]+=N_peak[i];
  }
  printf("All peaks: %d\n\n", N_peak[NC]);

  char Param[20000], tmp[200];
  sprintf(Param, "# Pars: SCORE=%d LOCAL=%d Y_MAX=%.3g EXP_MIN=%.3g WP_MIN=%d "
	  "THR=%.3g THR_MIN=%.3g LMIN_peak=%d LMIN_excl=%d SYMMETRIC=%d\n"
	  "# (the threshold decreases of %.3g for each unit of WP "
	  "until the minimum value min= %.3g)\n",
	  SCORE, LOCAL, Y_MAX, EXP_MIN, WP_MIN, THR,
	  THR_MIN,LMIN_peak,LMIN_excl, SYMMETRIC,THR_STEP,THR_MIN);
  sprintf(tmp, "# SEPARATE_NORM=%d SEP_SMOOTH=%d ZTNUM=%d\n",
	  SEPARATE_NORM, SEP_SMOOTH, ZTNUM);
  strcat(Param, tmp);
  strcat(Param, "# The score for selecting parameters is the sum of ");
  if(SCORE==0){sprintf(tmp, "bins in Peaks\n");}
  else if(SCORE==2){sprintf(tmp, "Z-score of bins / sqrt(num. Peaks)\n");}
  else{sprintf(tmp, "Z-score of bins in Peaks\n");}
  strcat(Param, tmp);
  sprintf(tmp, "# Optimized Damp factor: "
	  "%.0f (min=%.0f max= %.0f step=%.0f %d iter)\n",
	  Damp, DAMP_MIN, DAMP_MAX, DAMP_STEP, N_Damp);
  strcat(Param, tmp);
  sprintf(tmp, "# Optimized -log(EPS): "
	  "%.2f (min=%.2f max=%.2g step=%.2g %d iter)\n",
	  l_EPS, l_EPS_MIN, l_EPS_MAX, lE_STEP, N_EPS);
  strcat(Param, tmp);

  sprintf(tmp, "# Resulting window for smoothing: "
	  "%d bins %d nucleotides\n",(2*W_opt+1),(2*W_opt+1)*BOX1);
  strcat(Param, tmp);
  
  strcat(Param, "# The score for selecting range is ");
  if(SCORE_R=='S'){sprintf(tmp, "the same score as for other parameters\n");}
  else if(SCORE_R=='N'){sprintf(tmp, "the number of peaks\n");}
  else if(SCORE_R=='s'){sprintf(tmp, "the score/number of peaks\n");}
  else{sprintf(tmp, "the score / sqrt(number of peaks)\n");}
  strcat(Param, tmp);
 
  sprintf(tmp, "# Optimized Range for local average: "
    "%.0f (min=%.0f max=%.0f step=%.0f %d iter)\n",
    RANGE, R_MIN, R_MAX, R_STEP, N_RANGE);
  strcat(Param, tmp);

  N=N_peak[NC];
  int Nm=0; for(int i=0; i<nn; i++){if(lpeak[i])Nm++;}
  sprintf(tmp, "# Final score: %.0f %d Peaks S=%.3g "
	  "S/sqrt(N)=%.4g normalized=%.4g\n#\n",
	  S_opt, N, Score_opt, Score_opt/sqrt(N),
	  Score_opt*sqrt(N)/(float)Nm);
  strcat(Param, tmp);
  fprintf(filelog, "%s", Param);

  char string[1000];
  sprintf(string, "#c Npeaks size ");
  for(p=0; p<N_prof; p++)sprintf(string, "%s %s", string, name_prof[p]);
  fprintf(filelog, "%s\n", string);

  if(N_peak_ref){
    sprintf(string,
	    "%s d_12 d_21 over_1 over_2 prod o1/ran o2/ran prod/ran", string);
  }

  int WINDOW=(2*WIN+1)*BOX1;
  peaks[NC]=Peaks2Frags(N_peak+NC, lpeak, yout_peak, y_scr, nn, NC);
  char out_name[200];
  sprintf(out_name, "%s_T%.2g", NAMEOUT, Thr);
  if(MAKECLUST)Cluster_peaks(peaks[NC], N_peak[NC], out_name, MAKECLUST);

  for(int c=1; c<=NC; c++){
    fprintf(filelog, "%d", c);
    if(c<NC){peaks[c]=CopyPeaks(N_peak+c, peaks[NC], N_peak[NC], c);}
    if(peaks[c]==NULL){continue;}
    struct peak *peak=peaks[c];
    for(i=0; i<N_peak[c]; i++){Set_peak(peak, 0); peak++;}
    
    // Peaks are printed here!
    struct peak *cpeaks;
    if(CENTER){
      // Center peaks in such a way that the center coincides with the maximum
      // of the not-smoothed score and the peak is included in the previous one
      cpeaks=Center_peaks(N_peak+c, peaks[c], y_scr, x_coord, nn);
    }else{
      cpeaks=peaks[c];
    }
    Examine_peaks(cpeaks, N_peak[c], peak_ref, N_peak_ref, Nchr, //centered
		  xprof, yprof, nprof, nchr_prof, name_prof,
		  ave_prof, dev_prof, N_prof,
		  NAMEOUT, WINDOW, Thr, Param, "", file_ref, y1, c);
    if(0 && CENTER){
      printf("\nExamine not centered peaks\n"); //not centered
      Examine_peaks(peaks[c], N_peak[c], peak_ref, N_peak_ref, Nchr,
		    xprof, yprof, nprof, nchr_prof, name_prof,
		    ave_prof, dev_prof, N_prof,
		    NAMEOUT, WINDOW, Thr, Param, "NotCentered_",
		    file_ref, y1, c);
    }

    if(CENTER)free(cpeaks);
  }
    
  printf("\nBest: Score= %.0f (%d) Np= %d "
	 "W=%d Damp=%.0f l_EPS=%.3g range= %.0f "
	 "wp= %d ythr=%.2g*(AVE/BOX1) Z_thr=%.2g local= %d ztnum= %d "
	 "S/sqrt(N)=%.4g\n",
	 Score_opt, SCORE, N_opt, (2*W_opt+1)*BOX1, Damp_opt, lE_opt,range_opt,
	 WP_MIN,EXP_MIN,THR,LOCAL,ZTNUM,Score_opt/sqrt(N_opt));

  return(0);
}

long Read_sequence(char **seq, char *file_chr)
{
  long l=0;
  FILE *file_in=fopen(file_chr, "r");
  if(file_in==NULL){
    printf("ERROR, file %s not found\n", file_chr); exit(8);
  }
  char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>')continue;
    //l+=Count_char(string);
    {    // Count_char
      int n=0; char *s1=string;
      while(*s1!='\n'){s1++; n++;}
      l+=n;
    }
  }
  fclose(file_in);
  printf("%ld residues read in %s\n", l, file_chr);

  (*seq)=malloc((l+1)*sizeof(char)); l=0;
  file_in=fopen(file_chr, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>')continue;
    { //l+=Copy_char(string, *seq+l);
      int n=0; char *s1=string, *s2=*seq+l;
      while(*s1!='\n'){*s2=Maiuscule(*s1); s1++; s2++; n++;}
      l+=n;
    }
  }
  (*seq)[l]='\0';
  fclose(file_in);
  return(l);
}

void AT_norm_control(char *file_c, long *x_coord, float *y_contr,
		     char **chromosome, char **chr_name, long *Lchr,
		     long *nnch, int Nchr)
{
  long n=0; double AVE=0;
  int step=x_coord[1]-x_coord[0];
  for(int k=0; k<Nchr; k++){ // for each chromosome
    long nk=nnch[k], n1=nk-1, ini=0, L_thr=Lchr[k]+2000;
    char *nuc=chromosome[k];
    if(x_coord[n+nk-1]>L_thr){
      printf("ERROR, too many bp in control in chrom %d: %ld > %ld\n",
	     k, x_coord[n+nk-1], Lchr[k]); exit(8);
    }
    for(long i=0; i<nk; i++){ // for each bin
      long end; if(i<n1){end=x_coord[n];}else{end=ini+step;}
      int A=0, T=0; 
      for(long j=ini; j<end; j++){ // for each nucleotide;
	if(*nuc=='A'){A++;}else if(*nuc=='T'){T++;}
	nuc++;
      }
      if(A>T){y_contr[n]*=A;}else{y_contr[n]*=T;}
      AVE+=y_contr[n];
      ini=end; n++;
    }
  }
  AVE/=n;
  for(long i=0; i<n; i++){y_contr[i]/=AVE;}

  // Write oututput
  char out[200];
  // Look for last dot and chr
  /*if(0){
    char *s=file_c, *s2=NULL, tmp;
    while(*s!='\0'){
      if(strncmp(s, "chr", 3)==0 || strncmp(s, "Chr", 3)==0){s2=s;}
      s++;
    }
    while(*s!='\0'){s++;}
    s--; while(*s!='.'){s--;}
    if(s2){if(*(s2-1)=='_' || *(s2-1)=='.'){s2--;}  tmp=*s2; *s2='\0';}
    if(*s=='.'){
      *s='\0'; sprintf(out, "%s.AT_norm.wig", file_c); *s='.';
    }else{
      sprintf(out, "%s.AT_norm.wig", file_c);
    }
    if(s2){*s2=tmp;}
    }else{ */
  sprintf(out, "Control.AT_norm.wig", file_c); 
  //}

  printf("Writing normalized control in %s\n", out);  
  FILE *file_out=fopen(out, "w");
  int k=-1; long nk=0;
  for(long i=0; i<n; i++){
    if(k<0 || i==nk){
      k++; if(k<Nchr){nk+=nnch[k];}else{nk*=5;}
      fprintf(file_out, "fixedStep chrom=%s start=1 step=%d span=%d\n",
	      chr_name[k], step, step);
    }  
    fprintf(file_out, "%.4f\n", y_contr[i]);
  }
  fclose(file_out);

}

void Examine_peaks(struct peak *peaks, int N_peak,
		   struct peak *peak_ref, int N_peak_ref, int Nchr,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_prof, char **name_prof,
		   double **ave_prof, double **dev_prof, int N_prof,
		   char *NAMEOUT, int WINDOW, float Thr,
		   char *Param, char *what, char *file_ref,
		   float *y1, int clus)
{
  int PRINT=1; // Print peaks in bed format
  // Statistics of peaks
  int cl=lab_cl[clus]; // Index of cluster

  char name_all[200], out_name[200];
  sprintf(out_name, "%s%s_T%.2g_W%d", NAMEOUT, what, Thr, WINDOW);
  sprintf(name_all, "Properties_ave_%s.dat", out_name);
  FILE *file_all;

  char header[10000], header2[10000]="", tmp[100];
  sprintf(header, "#### For each variable: x=ave x+1=s.e.m. x+2=t\n");
  strcat(header, "#1=class 2=num 3=size 4=sd"); int a=5;
  strcpy(header2, "#1=class 2=num 3=size");
  for(int i=0; i<N_prof; i++){
      sprintf(tmp," %d=%s", a, name_prof[i]);
      strcat(header, tmp);
      sprintf(tmp," %d=%s", i+4, name_prof[i]);
      strcat(header2, tmp);
      a+=3;
  }
  if(OPEN_FILE){ 
    file_all=fopen(name_all, "w");
    fprintf(file_all, "%s", Param);
    fprintf(file_all, "#### %s:\n", NAMEOUT);
    fprintf(file_all, "%s\n", header);
    OPEN_FILE=0;
  }else{
    file_all=fopen(name_all, "a");
  }
  fprintf(file_all, "#class %d y= %.3g\n", cl, y1[clus]);
 
  double *Peak_score[N_prof]; int p;
  if(clus==NC){
    for(p=0; p<N_prof; p++)Peak_score[p]=malloc(N_peak*sizeof(double));
  }else{
    for(p=0; p<N_prof; p++)Peak_score[p]=NULL;
  }

  fprintf(file_all, "%d\t", cl);
  printf("%s\n%d\t", header2, cl);
  Profile_score_ZP(Peak_score, file_all,
		   xprof, yprof, nprof, nchr_prof, ave_prof, dev_prof, N_prof,
		   peaks, N_peak);

   char name_out[200];
  // Compare with references, if any
  if(N_peak_ref){
    
    /*float dist=Dist_peaks(peaks,N_peak,peak_ref,N_peak_ref);
      printf(" %.0f", dist);
      dist=Dist_peaks(peak_ref,N_peak_ref,peaks,N_peak);
      printf(" %.0f", dist); */
    
    // Statistics of peaks present / not present in previous set
    int N_peak_max=N_peak; if(N_peak_ref>N_peak_max)N_peak_max=N_peak_ref;
    struct peak *peak2=malloc(N_peak_max*sizeof(struct peak));
    int n_match=0, n_match_ref=0;

    // Overlap with random peaks
    peak2=Extract_peaks_ZP(peak_ref, N_peak_ref, Nchr, len_chr);
    Count_matches_ZP(&n_match,&n_match_ref,DTOL,peaks,N_peak,peak2,N_peak_ref);
    float ov1r=(float)n_match/N_peak; 
    float ov2r=(float)n_match_ref/N_peak_ref;
    if(ov1r==0){ov1r=0.5/(float)N_peak;}
    
    // Overlap with reference
    Count_matches_ZP(&n_match,&n_match_ref,DTOL,peaks,N_peak,
		     peak_ref,N_peak_ref);
    float ov1=(float)n_match/N_peak;
    float ov2=(float)n_match_ref/N_peak_ref;
    //int Np2=Select_peaks_match(peak2, peaks, N_peak, 1); // matched

    printf("Overlap: %d peaks %.3f %.3f",  n_match, ov1, ov2);
    printf(" /ran: %.1f %.1f\n", ov1/ov1r, ov2/ov2r);
    fprintf(file_all, "#overlap: %d %.3f %.1f\n", n_match, ov1, ov1/ov1r);
    
    if(clus<NC){goto end_reference;}
    
    fprintf(file_all, "#### ReferencePeaks %s:\n", file_ref);
    fprintf(file_all, "%d\t", cl);
    printf("\n#### ReferencePeaks:\n%s\n%d\t", header2, cl);
    Profile_score_ZP(NULL, file_all,
		     xprof, yprof, nprof, nchr_prof, ave_prof, dev_prof, N_prof,
		     peak_ref,N_peak_ref);
    
    if(0){
      int Np2=Select_peaks_match(peak2, peaks, N_peak, 0); // not matched
      sprintf(name_out, "Properties_Common%s.dat", NAMEOUT);
      fprintf(file_all, "#### CommonPeaks:\n");	
      fprintf(file_all, "%d\t", cl);
      printf("#### CommonPeaks:\n%s\n%d\t", header2, cl);
      Profile_score_ZP(NULL, file_all,
		       xprof, yprof, nprof,nchr_prof,ave_prof,dev_prof,N_prof,
		       peak2, Np2);
    }
    
    if(0){
      int Np2=Select_peaks_match(peak2, peaks, N_peak, 0); // not matched
      fprintf(file_all, "#### NewPeaks:\n");
      fprintf(file_all, "%d\t", cl);
      printf("#### NewPeaks:\n%s\n%d\t", header2, cl);
      Profile_score_ZP(NULL, file_all,
		       xprof,yprof,nprof,nchr_prof,ave_prof,dev_prof, N_prof,
		       peak2, Np2);
    }
    
    if(0){
      int Np2=Select_peaks_match(peak2, peak_ref, N_peak_ref, 0); // not matched
      fprintf(file_all, "#### UnconfirmedPeaks:\n");
      fprintf(file_all, "%d\t", cl);
      printf("#### UnconfirmedPeaks:\n%s\n%d\t", header2, cl);
      Profile_score_ZP(NULL, file_all,
		       xprof,yprof,nprof,nchr_prof,ave_prof,dev_prof, N_prof,
		       peak2, Np2);
    }

    if(PRINT && 0){
      sprintf(name_out, "%s_notfound.bed", file_ref);
      Print_Peaks_nomatch(peak_ref, N_peak_ref, name_out);
      
      char nameold[200];
      sprintf(nameold, "%s_OLD.bed", out_name);
      sprintf(name_out, "%s_NEW.bed", out_name);
      Print_Peaks_new(peaks, name_out, nameold);
    }
    
  end_reference:
    free(peak2);
  } // end reference

  if(0 && RANDOM){
    // Random peaks
    struct peak *peakran=Extract_peaks_ZP(peaks, N_peak, Nchr, len_chr);
    fprintf(file_all, "#### Random:\n-1\t");
    printf("\n#### Random:\n%s\n%d\t", header2, cl);
    Profile_score_ZP(Peak_score, file_all,
		     xprof, yprof, nprof, nchr_prof, ave_prof, dev_prof, N_prof,
		     peakran,N_peak);
    free(peakran);
  }

  printf("\n");
  fclose(file_all);

  // Print peaks
  if(PRINT && (PRINT_CLASS || clus==NC)){
    if(clus<NC){
      sprintf(name_out, "%s_c%d.bed", out_name, cl);
    }else{
      sprintf(name_out, "%s_All.bed", out_name);
    }
    Print_Peaks(peaks, name_out);
  }
 
  // Print properties of peaks
  if(clus==NC){
    sprintf(name_out, "Properties_Peaks_%s.dat", out_name);
    Print_peak_score(peaks, Peak_score, N_prof, name_prof, name_out);
  }

 
  /************* Metaplots ******************/
  if(0){
    sprintf(name_out, "Metaplots_%s.dat", NAMEOUT); //, SIZE_MIN
    printf("\n");
    for(p=0; p<N_prof; p++){
      Plot_profile(peaks, p, xprof[p], yprof[p], nprof[p], nchr_prof[p],
		   name_prof[p], name_out);
      printf("Metaplot of %s\n", name_prof[p]);
    }
    if(N_peak_ref){
      for(p=0; p<N_prof; p++){
	Plot_profile(peak_ref, p, xprof[p], yprof[p], nprof[p], nchr_prof[p],
		     name_prof[p], "Metaplots_reference.dat");
      }
    }
  }
  if(clus==NC){printf("Writing %s\n", name_all);}

}

int Get_input(int *nc, char **file_c, char **file_e, // control and exper
	      char **file_ref, // reference peaks
	      char **file_genome, char **chr_lab, int *NORM_AT, // genome
	      char ***file_prof, int *N_prof, // profiles
	      char **name_prof, int *nfile_prof,
	      // parameters: 
	      int *NC, float *THR, float *thr_step, float *thr_min,
	      float *EXP_MIN, float *Y_MAX,
	      int *LOCAL, int *ZTNUM, 
	      int *SEPARATE_NORM, int *SEP_SMOOTH,
	      int *WP_MIN, int *SCORE, char *SCORE_R, float *OUTL,
	      float *DAMP_MIN, float *DAMP_STEP, float *DAMP_MAX,
	      float *EPS, int *LMIN_excl, int *LMIN_peak,
	      float *R_MIN, float *R_MAX, int *N_RANGE,
	      int *PRINT_CLASS, int *PRSCORE, char *NAMEOUT,
	      //int *KMEAN, int *SEP_KMEAN, float *THR_MIN,
	      //int *SIZE_MIN, int *DCLUST, int *DTOL, 
	      //int *ZRANGE, int *WGT_INI, int *WGT_END, int *WGT_STEP,
	      int argc, char **argv)      // input
{
  int i; char file[200];
  if(argc<2)help(argv[0]);
  for(i=0; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){help(argv[0]);}
    else if(strncmp(argv[i], "-file", 5)==0){strcpy(file, argv[i+1]);}
    //if(i>1)printf("WARNING, unrecognized option %s\n", argv[i]);
  }

  // Default

  // Reading parameter file
  if(file[0]!='\0'){printf("Reading %s\n", file);}
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }

  char string[1000], dumm[80]; *nc=0;
  int ne=0, np=0, n1=0, npf=0, ng=0;
  int re=0, rc=0, rp=0, r1=0, rg=0;
  char PATH[NCHAR]; PATH[0]='\0';
  double norm, NSAM; int  wd; float x;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'){
      continue; // Comment line
    }else if(strncmp(string, "END", 3)==0){
      rp=0; r1=0; re=0; rc=0; rg=0; PATH[0]='\0';
      if(n1){nfile_prof[npf]=n1; npf++; n1=0;}
    }else if(strncmp(string, "DIR", 3)==0){
      char *s=string+4; int k=0;
      while(*s!='\n' && *s!=' '){PATH[k]=*s; k++; s++;}
      *s--; if(*s!='/'){PATH[k]='/'; k++;} PATH[k]='\0';
    }else if(strncmp(string, "EXPER", 5)==0){
      re=1; ne=0;
    }else if(re){
      Read_name(&ne, file_e, PATH, string, NCHRMAX);
    }else if(strncmp(string, "GENOME", 6)==0){
      rg=1; ng=0;
    }else if(rg){
      Get_chr_name(&(chr_lab[ng]), string);
      Read_name(&ng, file_genome, PATH, string, NCHRMAX);
    }else if(strncmp(string, "CONTROL:", 8)==0){
      rc=1; *nc=0; 
    }else if(rc){
      Read_name(nc, file_c, PATH, string, NCHRMAX);
    }else if(strncmp(string, "PREDICTION", 8)==0){
      rp=1;
    }else if(rp){
      if(np){
	printf("ERROR, only one reference file allowed\n"); exit(8);
      }
      Read_name(&np, file_ref, PATH, string, NCHRMAX);
    }else if(strncmp(string, "PROF", 4)==0){
      r1=1;
      name_prof[npf]=malloc(NCHAR*sizeof(char));
      sscanf(string+5, "%s", name_prof[npf]);
    }else if(r1){
      Read_name(&n1, file_prof[npf], PATH, string, NCHRMAX);
    }else if(strncmp(string, "NAME", 4)==0){
      sscanf(string+5, "%s", NAMEOUT);
    }else if(strncmp(string, "THR_STEP=", 9)==0){
      sscanf(string+9, "%f", thr_step);
      printf("THR_STEP= %.2f\n", *thr_step);
    }else if(strncmp(string, "THR_MIN=", 8)==0){
      sscanf(string+8, "%f", thr_min);
      printf("THR_MIN= %.2f\n", *thr_min);
    }else if(strncmp(string, "THR=", 4)==0){
      sscanf(string+4, "%f", THR);
      printf("THR= %.2f\n", *THR);
    }else if(strncmp(string, "SCORE_RANGE", 11)==0){
      char tmp; sscanf(string+12, "%c", &tmp);
      if(tmp=='X' || tmp=='S' || tmp=='s'|| tmp=='N'|| tmp=='D'){
	*SCORE_R=tmp; printf("Score_Range= %c\n", tmp);
      }else{
	printf("WARNING, SCORE_RANGE %c not allowed, using default %c\n",
	       tmp,*SCORE_R);
      }
    }else if(strncmp(string, "SCORE", 5)==0){
      int tmp; sscanf(string+6, "%d", &tmp);
      if(tmp==0 || tmp==1 || tmp==2){*SCORE=tmp;}
      else{
	printf("WARNING, SCORE %d not allowed, using default %d\n",tmp,*SCORE);
      }

    }else if(strncmp(string, "NORM_AT", 7)==0){
      int yes; sscanf(string+8, "%d", &yes);
      if(yes){*NORM_AT=1;}
    }else if(strncmp(string, "OUTL", 4)==0){
      sscanf(string+5, "%f", OUTL);
    }else if(strncmp(string, "SEPARATE_NORM", 13)==0){
      sscanf(string+14, "%d", SEPARATE_NORM);
    }else if(strncmp(string, "SEP_SMOOTH", 10)==0){
      sscanf(string+11, "%d", SEP_SMOOTH);
    }else if(strncmp(string, "NC", 2)==0){
      int dumm; sscanf(string+3, "%d", &dumm);
      if((dumm<2)||(dumm>11)){
	printf("ERROR: Not allowed value NC=%d %s", dumm, string);
	printf("Using default %d\n", *NC);
      }else{
	*NC=dumm;
      }
    }else if(strncmp(string, "EXP_MIN", 7)==0){
      sscanf(string+8, "%f", EXP_MIN);
    }else if(strncmp(string, "Y_MAX", 5)==0){
      float dumm; sscanf(string+6, "%f", &dumm);
      if(dumm<0 || dumm >1){printf("WARNING, Y_MAX= %.3g not allowed\n",dumm);}
      else{*Y_MAX=dumm;}
    }else if(strncmp(string, "ZTNUM", 5)==0){
      sscanf(string+6, "%d", ZTNUM);
    }else if(strncmp(string, "LOCAL", 5)==0){
      sscanf(string+6, "%d", LOCAL);
    }else if(strncmp(string, "WP_MIN", 6)==0){
      sscanf(string+7, "%d", WP_MIN);
    }else if(strncmp(string, "DAMP_MIN", 8)==0){
      sscanf(string+9, "%f", DAMP_MIN);
    }else if(strncmp(string, "DAMP_MAX", 8)==0){
      sscanf(string+9, "%f", DAMP_MAX);
    }else if(strncmp(string, "DAMP_STEP", 9)==0){
      sscanf(string+10, "%f", DAMP_STEP);
    }else if(strncmp(string, "EPS", 3)==0){
      sscanf(string+4, "%f", EPS);
    }else if(strncmp(string, "LMIN_e", 6)==0){
      sscanf(string+7, "%d", LMIN_excl);
    }else if(strncmp(string, "LMIN_p", 6)==0){
      sscanf(string+7, "%d", LMIN_peak);
    }else if(strncmp(string, "RANGE_MAX", 9)==0){
      sscanf(string+10, "%f\n", R_MAX);
      printf("R_MAX= %f\n", *R_MAX);
    }else if(strncmp(string, "RANGE_MIN", 9)==0){
      sscanf(string+10, "%f", R_MIN);
      printf("R_MIN= %f\n", *R_MIN);
    }else if(strncmp(string, "N_RANGE", 7)==0){
      sscanf(string+8, "%d", N_RANGE);
      printf("N_RANGE= %d\n", *N_RANGE);
    }else if(strncmp(string, "PRINT_CLASS", 11)==0){
      sscanf(string+12, "%d", PRINT_CLASS);
    }else if(strncmp(string, "PRINT_SCORE", 11)==0){
      sscanf(string+12, "%d", PRSCORE);

      //}else if(strncmp(string, "WGT_INI", 7)==0){
      //sscanf(string+8, "%d", WGT_INI);
      //}else if(strncmp(string, "WGT_END", 7)==0){
      //sscanf(string+8, "%d", WGT_END);
      //}else if(strncmp(string, "WGT_STEP", 8)==0){
      //sscanf(string+9, "%d", WGT_STEP);
      //}else if(strncmp(string, "ZRANGE", 6)==0){
      //sscanf(string+7, "%d", ZRANGE);

      //}else if(strncmp(string, "WP_STEP", 7)==0){
      //sscanf(string+8, "%d", WP_STEP);
      //}else if(strncmp(string, "KMEAN", 5)==0){
      //sscanf(string+6, "%d", KMEAN);
      //}else if(strncmp(string, "SEP_KMEAN", 9)==0){
      //sscanf(string+10, "%d", SEP_KMEAN);
      //}else if(strncmp(string, "THR_MIN", 7)==0){
      //sscanf(string+8, "%f", THR_MIN);
      //printf("THR_MIN= %.2f\n", *THR_MIN);

      //}else if(strncmp(string, "SIZE_MIN", 8)==0){
      //sscanf(string+9, "%d", SIZE_MIN);
      //printf("Minimum size required for calling a peak= %d\n", *SIZE_MIN);
      //}else if(strncmp(string, "DCLUST=", 7)==0){
      //sscanf(string+7, "%d", DCLUST);
      //}else if(strncmp(string, "DTOL=", 5)==0){
      //sscanf(string+5, "%d", DTOL);
      //printf("DTOL= %d\n", *DTOL);
      /*}else if(strncmp(string, "SIZE_END", 8)==0){
      sscanf(string+9, "%d", SIZE_END);
      printf("Minimum size required (largest)= %d\n", *SIZE_END);
    }else if(strncmp(string, "SIZE_STEP", 9)==0){
      sscanf(string+10, "%d", SIZE_STEP);
      printf("Minimum size required (step)= %d\n", *SIZE_STEP);*/
      /*}else if(strncmp(string, "MODEL", 5)==0){
      sscanf(string+6, "%s", dumm);
      if((strcmp(dumm, "E")!=0)&&(strcmp(dumm, "G")!=0)){
	printf("WARNING, model %s not implemented\n", dumm);
	printf("ALLowed models are E (exponential) and G (Gaussian)\n");
	printf("Default is %s\n", MODEL);
      }else{
	printf("Model changed from %s to %s\n", MODEL, dumm);
	strcpy(MODEL, dumm);
	}*/
    }else if(string[0]!='\n'){
      printf("WARNING, unrecognized line: %s\n", string);
    }
  }
  if(ng){printf("%d genome files to read\n", ng);}


  for(i=0; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){
      help(argv[0]);
    }else if(strncmp(argv[i], "-file", 5)==0){
      i++; strcpy(file, argv[i]);
    }else if(strncmp(argv[i], "-at_norm", 8)==0){
      *NORM_AT=1;
    }else{printf("WARNING, unrecognized option %s\n", argv[i]);}
  }


  *N_prof=npf;
  if((ne==0)){
    printf("ERROR, no input files specified\n"); help(argv[0]);
  }
  /*if((ne>1)||(*nc>1)){
    printf("ERROR, only one wig file allowed for experiment and control,");
    printf(" found %d and %d\n", ne, nc); exit(8);
    }*/
  printf("%d files for experiment and %d for control\n", ne, *nc);
  if(*nc>1 && *nc!=ne){
    printf("ERROR, if control is present we expect the same number of files ");
    printf("as the experiment\n"); exit(8);
  }
  if(*nc==0){
    printf("WARNING, no control files specified\n");
    printf("Using mean experiment box as control\n");
  }

  printf("Executing %s\n", CODE);
  return(ne);  
}

void Get_chr_name(char **chr_lab, char *string){
  *chr_lab=malloc(20*sizeof(char));
  char *s1=string, *s2=*chr_lab;
  while(*s1!='\0'){
    if(strncmp(s1, "chr", 3)==0 || strncmp(s1, "Chr", 3)==0){
      s1+=3; break;
    }
    s1++;
  }
  while(*s1!='.' && *s1!='\0'){*s2=*s1; s1++; s2++;}
  *s2='\0';
}

char Maiuscule(char s){
  int i=(int)s;
  if(i>96){return((char)(i-32));}else{return(s);}
}

void help(char *name){
  printf("\nPROGRAM %s\n", name);
  printf("Author: Ugo Bastolla,\n"
	 "Centro de Biologia Molecular Severo Ochoa\n"
	 "(CSIC-UAM), Madrid Spain\n"
	 "<ubastolla@cbm.csic.es>\n\n"
	 //
	 "Mandatory argument: -file <parameter file>\n"
	 "FORMAT of parameter file:\n"
	 "################################\n"
	 "EXPER:\n"
	 "DIR=<path of exper files>\n"
	 "<exper file> (one for each chromosome, MANDATORY)\n"
	 "END\n"
	 "CONTROL:\n"
	 "DIR=<path of control files>\n"
	 "<control file> (one for each chromosome, optional)\n"
	 "END\n"
	 "GENOME:\n"
	 "DIR=<path of genome files>\n"
	 "<genome file> (one for each chromosome, needed for AT norm)\n"
	 "END\n"
	 "PREDICTION:\n"
	 "DIR=<path of reference file>\n"
	 "<reference file> (only one, optional)\n"
	 "END\n"
	 "PROF <profile name>\n"
	 "DIR=<path of prof files>\n"
	 "<prof file> (only one, or one for each chromosome, optional)\n"
	 "END\n"
	 "(Same for all experimental profiles)\n");
  printf("############ Parameters:\n");
  printf("NORM_AT=0 or 1 (normalize control with T content)\n");
  printf("WP_MIN= <Min.length of bell-shaped region> def: %d\n", WP_MIN);
  printf("THR=<Threshold Z score for minimum WP> (default %.2f)\n", THR);
  printf("THR_STEP=<Decrease of threshold for unit of WP> (default %.2f)"
	 "until min. value %.2f\n", THR_STEP, THR_MIN);
  printf("LOCAL= <Normalize z score by local s.d.?> def: %d\n", LOCAL);
  printf("Y_MAX=%.3f ! "
	 "compute <y> excluding y>y_thr=y_min+Y_MAX*(y_max-y_min)\n", Y_MAX);
  printf("EXP_MIN=%.3f ! Exclude peak if exp/<exp> < %.3f\n",  EXP_MIN);
  printf("ZTNUM= <Threshold for Z score depends on number of frags> def: %d\n",
	 ZTNUM);

  printf("SEPARATE_NORM= <> ! Normalize chrosmosomes separately? 1=yes "
	 "def: %d\n", SEPARATE_NORM);
  printf("SEP_SMOOTH= <>   ! Smooth chrosmosomes separately? 1=yes "
	 "def: %d\n", SEP_SMOOTH);
  printf("SCORE= <>   ! Score for determining optimal window: number =0, sum score=1 /srt(N)=2 "
	 "def: %d\n", SCORE);
  printf("OUTL= <>    ! Compute score s.dev excluding outliers > OUTL "
	 "(OUT<0: disable) def: %d\n", OUTL);
  printf("EXP_MIN= <Min. reads per bp at peaks> def: %.2g\n", EXP_MIN);
  printf("DAMP_MIN= <Min. damping factor> def: %d\n", DAMP_MIN);
  printf("DAMP_STEP=<Increment of damping factor> def: %d\n", DAMP_STEP);
  printf("EPS= <Threshold for averaging bins> def: %.2g\n", EPS);
  printf("RANGE= <Min. number of nuc for averaging score> def: %.0f\n",RANGE);
  printf("RANGE_MAX= <Max. number of bins for averaging score> def: %d\n",
	 R_MAX);
  printf("N_RANGE= <Number of iterations for optimizing RANGE> def: %d\n",
	 N_RANGE);
  printf("PRINT_CLASS= <Print bed file for each class?> def: %d\n",PRINT_CLASS);
  printf("PRINT_SCORE= <Print Peaks score as wig file?> def: %d\n", PRSCORE);

  //printf("WP_STEP= <increment of tested value of wp> def: %d\n", WP_STEP);
  //printf("WGT_INI= <Min.tested value of wgt> def: %d\n", WGT_INI);
  //printf("WGT_END= <Max.tested value of wgt> def: %d\n", WGT_END);
  //printf("WGT_STEP= <increment of tested value of wgt> def: %d\n", WGT_STEP);
  //printf("ZRANGE= <half width for average local Z score> def: %d\n", ZRANGE);

  //printf("SIZE_MIN=<Minimum size for calling a peak>\n");
  //printf("DCLUST=<Distance threshold for joining fragments>\n");
  //printf("DTOL=<Tolerance for comparison>\n");
  //
  /*printf("SIZE_INI=<Minimum of minimum size for calling a peak>\n");
  printf("SIZE_END=<Maximum of minimum size for calling a peak>\n");
  printf("SIZE_INI=<Step of minimum size for calling a peak>\n");
  printf("MODEL=E ! Score distribution used in the HMM\n");
  printf("# Allowed: E (exponential) and G (Gaussian)\n");*/

  printf("Command line options:\n"
         "   -at_norm   ! Normalize control with AT content\n"
	 "\n"
	 );
  exit(8);
}

void Copy_prof(long *nprof, long **xprof, float **yprof,
	       long *x_scr, float *y_scr, long *nnch, int Nchr)
{
  int k; int i;
  long *xs=x_scr; float *ys=y_scr;
  for(k=0; k<Nchr; k++){
    long *xp=xprof[k]; float *yp=yprof[k];
    for(i=0; i<nnch[k]; i++){
      *xp=*xs; xp++; xs++;
      *yp=*ys; yp++; ys++;
    }
  }
}

void Comp_Z_score(float *Z_scr, long nn, float *ye_box, float *yc_box, float sd)
{
  double Z1=0, Z2=0; long i;
  float *ycb=yc_box, *yeb=ye_box, *Z=Z_scr, yc;
  for(i=0; i<nn; i++){
    float ZZ=(*yeb-*ycb);
    //float yc=0.5*(*ycb+1.0), ZZ=(*yeb-yc);
    //float ZZ=(*yeb-*ycb)/(*ycb);
    // WARNING, the last normalization presents a bias towards AT rich
    *Z=ZZ; Z1+=ZZ; Z2+=ZZ*ZZ;
    Z++; yeb++; ycb++;
  }
  Z1/=nn; Z2=sqrt(Z2/nn-Z1*Z1);
  if(Z2>0){
    for(i=0; i<nn; i++)Z_scr[i]=(Z_scr[i]-Z1)/Z2;
  }
}

int Set_control(float *y, long nn, float Mean){
  long i; float *yi=y;
  for(i=0; i<nn; i++){*yi=Mean; yi++;}
  return(0);
}

float Normalize_counts(float *y, long n, float AVE){
  float sd, Mean=Mean_counts(&sd, y, n);
  Rescale_counts(y, n, AVE/Mean);
  printf("Ave= %.3g sd= %.3g n=%ld\n", Mean, sd, n);
  return(sd*AVE/Mean);
}

float Mean_counts(float *sd, float *yy, long nn){
  double Y1=0, Y2=0; float *y=yy; long i;
  for(i=0; i<nn; i++){
    Y1+=*y; Y2+=(*y)*(*y); y++;
  }
  Y1/=nn; Y2=sqrt((Y2-nn*Y1*Y1)/(nn-1));
  *sd=Y2;
  return(Y1);
}

void Rescale_counts(float *yy, long nn, float scale){
  long i; float *y=yy;
  for(i=0; i<nn; i++){*yy *= scale; yy++;}
}

void Print_score(long *x, float *y, long *nnch, long nn,
		 char **chr_name, char *NAMEOUT)
{
  char nameout[200];
  sprintf(nameout, "%s_score.wig", NAMEOUT);
  printf("Writing %s\n", nameout);
  FILE *file_out=fopen(nameout, "w");
  fprintf(file_out, "track type=wiggle_0\n");
  int step=x[1]-x[0];
  int k=0; long i, n0=0;
  for(i=0; i<nn; i++){
    if(i==n0){
      fprintf(file_out, "fixedStep chrom=%s start=1 step=%d span=%d\n",
	      chr_name[k], step, step);
      n0+=nnch[k]; k++;
    }
    fprintf(file_out, "%.4f\n", y[i]);
  }
  fclose(file_out);
  printf("Written\n");
}


void Count_peaks(int *numclus, int *numdom, int *lpeak, long *endchr,
		 int N, int ncl)
{
  int i, k, k0=-1, ch=-1; long n0=0;
  for(k=0; k<ncl; k++){numclus[k]=0; numdom[k]=0;}
  for(i=0; i<N; i++){
    k=lpeak[i];
    if((k<0)||(k>=ncl)){
      printf("WARNING!!!! wrong lpeak identifier %d (%d lpeaks expected)\n",
	     k, ncl); continue;
    }
    if(ch<0 || i>=endchr[ch]){ch++; k0=-1;}
    if(k!=k0){numdom[k]++; k0=k;}
    numclus[k]++;
  }
}

int Close_peak(struct peak **peak, long x_coord, float y, float yout,
	       int *numch, double *sizech)
{
  numch[(*peak)->chr]++;
  sizech[(*peak)->chr]+=(x_coord-(*peak)->ini);
  (*peak)->end=x_coord-1;
  (*peak)->y=y;
  (*peak)->yout=yout;
  (*peak)->next=(*peak)+1;
  (*peak)++;
  return(1);
}

/*
void Compute_scores(int *lpeak, short *chr, float **x_VarSam,
		    int Nsam, int Nvar, int ncl, float lik,
		    struct Para Par, FILE *filelog)
{
  int numclus[ncl], numdom[ncl], k;
  Count_peaks(numclus, numdom, lpeak, chr, Nsam, ncl);
  printf("Peak sizes (peaks/no peaks): ");
  for(k=0; k<ncl; k++)printf(" %d", numdom[k]); printf(" domains ");
  for(k=0; k<ncl; k++)printf(" %d", numclus[k]); printf(" elements\n");

  // Compute likelihood and other scores
 float cscore=Peak_score(lpeak, ncl, x_VarSam, Nsam, Nvar);
 // Effective number of variables
 float lcorr=Get_lcorr(x_VarSam, Nsam, 0);
 float Nsam_eff=Nsam;
 if((lcorr > 0)&&(lcorr < Nsam))Nsam_eff/=lcorr;
 float norm_lik1=Nsam_eff/Nsam, norm_lik2=Nsam_eff;
 lik*=norm_lik1;
 // Compute number of parameters
 int N_para = Nvar + 1; // mu, tau (per lpeak)
 if(Par.sig){N_para+= Nvar*(Nvar+1)/2;}
 else if(Par.scale_pos){
   N_para+=Nvar;
   if(Rel_diff(Par.scale_neg[0][0],Par.scale_pos[0][0])>0.02)N_para++;
   if(Rel_diff(Par.scale_neg[1][0],Par.scale_pos[1][0])>0.02)N_para++;
 }
 if(Par.trans)N_para+=(ncl-1);
 N_para=ncl*N_para-1;

 float aic=AIC(lik, N_para, Nsam_eff);
 float bic=BIC(lik, N_para, Nsam_eff);
 lik/=norm_lik2; aic/=norm_lik2; bic/=norm_lik2;
 // Separation score
 double chi=-1;
 if(Par.sig){
   chi=Chi2(Par.mu[0], Par.mu[1], Par.sig[0], Par.sig[1], Nvar);
 }else if(Par.scale_pos){
   float sig0=max(Par.scale_pos[0][0], Par.scale_neg[0][0]);
   float sig1=max(Par.scale_pos[1][0], Par.scale_neg[1][0]);
   chi=Chi2(Par.mu[0], Par.mu[1], &sig0, &sig1, 1);
 }

 // Print
 char txt[400];
 sprintf(txt,
	 "# %d lpeaks: lik= %.4f AIC=%.4f BIC= %.4f (/NPC) score=%.3f\n",
	 ncl, lik, aic, bic, cscore);
 sprintf(txt,"%s# Discriminative power chi2= %.3f\n", txt, chi);
 printf("%s\n", txt);
 fprintf(filelog, "%s", txt);
 Print_parameters_f(numclus, &Par, ncl, Nvar, filelog);
}
*/

//
/*int K_means_chr(int *N, double *S, double **Mean, int *cluster, float *y_scr,
		//int *WIN, float Damp, float l_EPS
		int NC)
{
  // Smooth over window, compute Z score and assign clusters
  //Smooth_score(y_scr, nn, Damp, WIN, l_EPS);

  if(SEP_KMEAN){
    long n0=0; int chr; *N=0; *S=0; 
    for(chr=0; chr<Nchr; chr++){
      long nk=nnch[chr]; double S_chr=0;
      int N_chr=K_means_old(&S_chr, Mean[chr], cluster+n0, y_scr+n0, nk, NC);
      n0+=nk; (*N)+=N_chr; (*S)+=S_chr;
    }
  }else{ // All chromosomes together
    *N=K_means_old(S, Mean[0], cluster, y_scr, nn);
  }
  return(0);
}
*/


int K_means_old(int *Nc, double *Mean, double *Sd, float *ythr_clus,
		int *cluster, float *y_scr, int nk, int NC)
{
  int itmax=100, iter, NC1=NC-1, k, i;

  // Initialize
  float ymax=-100, ymin=100;
  for(i=0; i<nk; i++){
    if(y_scr[i]<ymin){ymin=y_scr[i];}else if(y_scr[i]>ymax){ymax=y_scr[i];}
  }
  float step=(ymax-ymin)/NC, yy=step;
  for(k=0; k<NC; k++){ythr_clus[0]=yy; yy+=step;}

  int Nc_old[NC]; double Mean_old[NC];
  for(k=0; k<NC;  k++){Mean_old[k]=0; Nc_old[k]=0;}
  for(iter=0; iter<100; iter++){
    for(k=0; k<NC;  k++){Mean[k]=0; Nc[k]=0;}
    for(i=0; i<nk; i++){
      for(k=0; k<NC1; k++)if(y_scr[i]<ythr_clus[k])break;
      Nc[k]++; Mean[k]+=y_scr[i];
    }
    for(k=0; k<NC; k++)if(Nc[k])Mean[k]/=Nc[k];
    for(k=0; k<NC; k++){
      if(Nc[k]!=Nc_old[k] || Mean[k]!=Mean_old[k])break;
    }
    if(k==NC)break; // convergence
    for(k=0; k<NC; k++){Nc_old[k]=Nc[k]; Mean_old[k]=Mean[k];}
    for(k=0; k<NC1; k++){ythr_clus[k]=(Mean[k]+Mean[k+1])/2;}
  }
  if(iter==itmax)
    printf("WARNING, K_means did not converge after %d steps\n",iter);

  for(k=0; k<NC;  k++){Nc[k]=0; Mean[k]=0; Sd[k]=0;}
  for(i=0; i<nk; i++){
    for(k=0; k<NC1; k++)if(y_scr[i]<ythr_clus[k])break;
    cluster[i]=k; Nc[k]++; Mean[k]+=y_scr[i]; Sd[k]+=y_scr[i]*y_scr[i]; 
  }
  for(k=0; k<NC; k++){
    if(Nc[k]>1){
      Mean[k]/=Nc[k]; Sd[k]=(Sd[k]-Nc[k]*Mean[k]*Mean[k])/(Nc[k]-1);
    }
    Sd[k]=sqrt(Sd[k]);
    printf("Cluster %d %d fragments Mean: %.3g S.d.: %.2g Max: %.2g\n",
	   k, Nc[k], Mean[k], Sd[k], ythr_clus[k]);
  }
  return(0);
}

int Get_peaks(int *N, double *S, int *lpeak, float *y_scr, int *WIN,
	      float Damp, float Thr, float l_EPS, int SIZE_MIN)
{
  // Smooth over window, compute Z score and assign clusters
  Smooth_score(y_scr, nn, Damp, WIN, l_EPS);

  // Remove small domains
  *N=Get_frags(lpeak, S, Thr, y_scr, SIZE_MIN);
  printf("Thr=%.2f W=%d D=%.0f N=%d (S>%d)\n",
	 Thr, *WIN, Damp, *N, SIZE_MIN);
  return(0);
}

void Smooth_score(float *y_scr, long nn, float Damp, int *WIN, float l_EPS)
{
  double Z1=0, Z2=0;
  long n0=0; int i, k, nst=0;
  float Sqr_w;

  for(k=0; k<Nchr; k++){
    long nk=nnch[k];
    if(k==0)printf("W= %d Damp=%.1f\n", *WIN, Damp);
    if((Omit_chr)&&(Omit_chr[k])){n0+=nk; continue;} 
    Smooth(y_scr+n0, y_diff+n0, nk, &Sqr_w, Damp, WIN, l_EPS);
    double Y1=0, Y2=0; float *yy=y_scr+n0;
    for(i=0; i<nk; i++){Y1+=*yy; Y2+=(*yy)*(*yy); yy++;}
    Z1+=Y1; Z2+=Y2; nst+=nk;
    if(SEP_SMOOTH){
      int ns=nk; yy=y_scr+n0;
      if(OUTL>0)ns=Get_ave_sd_outliers(&Y1, &Y2, y_scr+n0, nk, OUTL);
      Y1/=ns; Y2=sqrt(Y2/ns-Y1*Y1);
      yy=y_scr+n0; for(i=0; i<nk; i++)yy[i]=(yy[i]-Y1)/Y2;
      printf("chr %d Av. diff. smoothed: %.3f s.d. %.3g", k+1,Y1,Y2);
      if(OUTL>0)printf(" (removing outliers > %.2f)", OUTL);
      printf("\n");
    }
    n0+=nk;
  }
  if(SEP_SMOOTH==0){
    if(OUTL>0)nst=Get_ave_sd_outliers(&Z1, &Z2, y_scr, nn, OUTL);
    Z1/=nst; Z2=sqrt(Z2/nst-Z1*Z1);
    for(i=0; i<nn; i++)y_scr[i]=(y_scr[i]-Z1)/Z2;
    printf("Ave diff smoothed exper-control: %.2g s.d.: %.2g", Z1, Z2);
    if(OUTL>0)printf(" (removing outliers > %.2f)", OUTL);
    printf("\n");
  }
}

int Get_ave_sd_outliers(double *Y1, double *Y2, float *y, int nk, float OUTL)
{
  if(OUTL<=0)return(nk);
  float ave=*Y1/nk, sd=sqrt(*Y2/nk-ave*ave), thr=sd*OUTL;
  int ns=0, i; float *yy=y;
  *Y1=0; *Y2=0;
  for(i=0; i<nk; i++){
    if(fabs(*yy-ave)<=thr){(*Y1)+=*yy; (*Y2)+=(*yy)*(*yy); ns++;} yy++;
  }
  return(ns);
}

float Smooth(float *y_weight,
	     float *y, long N, float *Sqr_w, float DAMP, int *WIN, float l_EPS)
{
  /* Weighted average of the number of reads in neighboring boxes
   */
  double Y1=0, Y2=0;

  float DampFact=exp(-BOX1/DAMP);
  *WIN=l_EPS*DAMP/BOX1;

  float *Damp=malloc((*WIN+1)*sizeof(float));
  Damp[0]=1.0;
  (*Sqr_w)=1.0; int k;
  for(k=1; k<= *WIN; k++){
    Damp[k]=Damp[k-1]*DampFact;
    (*Sqr_w)+=2*Damp[k];
  }
  (*Sqr_w)=sqrt(*Sqr_w);

  float *y_w=y_weight, *y_center=y;
  for(long i=0; i<N; i++){
    float ww=1, ysum=(*y_center), *w, *yy;
    w=Damp+1; yy=y_center+1; long j=i+1;
    for(k=1; k<=*WIN; k++){
      if(j>=N)break;
      ysum+=(*yy)*(*w); ww+=(*w); w++; j++; yy++;
    }
    w=Damp+1; yy=y_center-1; j=i-1; 
    for(k=1; k<=*WIN; k++){
      if(j<0)break;
      ysum+=(*yy)*(*w); ww+=(*w); w++; j--; yy--;
    }
    *y_w=(ysum/ww); Y1+=*y_w; Y2+=(*y_w)*(*y_w);
    y_w++; y_center++;
  }
  free(Damp);
  Y1/=N; Y2=(Y2-N*Y1*Y1); Y2=sqrt(Y2/(N-1));
  return(Y2);
}

float Get_likelihood(float Thr, int *lpeak, float *y, long nn, int SIZE_MIN)
{
  double S; int Np=Get_frags(lpeak, &S, Thr, y, SIZE_MIN);
  float lik=Peak_likelihood(lpeak, y, nn, 3)/nn;
  printf("Thr=%.2f likelihood/N: %.4f\n", Thr, lik);
  return(lik);
}

int Get_frags(int *lpeak, double *Score, float Thr, float *y_scr, int SIZE_MIN)
{
  int Np=0, i; *Score=0;
  int S=SIZE_MIN/BOX1, m; long n0=0, j;
  for(i=0; i<Nchr; i++){
    long ni=nnch[i];
    if((Omit_chr)&&(Omit_chr[i])){
      for(j=n0; j<(n0+ni); j++)lpeak[j]=0;
      n0+=ni; continue;
    }
    for(j=n0; j<(n0+ni); j++){
      if(y_scr[j]>Thr){Np++; *Score+=y_scr[j]; lpeak[j]=NC1;}
      else{lpeak[j]=0;}
    }
    int k0=-1, k, j0=n0;
    for(j=n0; j<(n0+ni); j++){
      k=lpeak[j];
      if(k!=k0){
	if(k==0){ // Peak ends, check whether size >= S
	  if((j-j0)<S)
	    for(m=j0; m<j; m++){lpeak[m]=0; *Score-=y_scr[m]; Np--;}
	}else{ // Peak starts
	  j0=j;
	}
	k0=k;
      }
    }
    if((k==NC1)&&((j-j0)<S))
      for(m=j0; m<j; m++){lpeak[m]=0; *Score-=y_scr[m]; Np--;}
    n0+=ni;
  }
  return(Np);
}

float Peak_likelihood(int *lpeak, float *y, int nn, int ncl)
{
  printf("Computing likelihood\n");
  double *y1=malloc(ncl*sizeof(double));
  double *y2=malloc(ncl*sizeof(double));
  double *ll=malloc(ncl*sizeof(double));
  int *nc=malloc(ncl*sizeof(int)), i, k;
  for(k=0; k<ncl; k++){y1[k]=0; y2[k]=0; nc[k]=0;}
  for(i=0; i<nn; i++){
    k=lpeak[i];
    if((k<0)||(k>=ncl))printf("ERROR, cluster %d out of bound\n", k);
    y1[k]+=y[i]; y2[k]+=y[i]*y[i]; nc[k]++;
  }
  for(k=0; k<ncl; k++){
    y1[k]/=nc[k]; y2[k]=y2[k]/nc[k]-y1[k]*y1[k]; ll[k]=0;
  }
  for(i=0; i<nn; i++){
    k=lpeak[i]; float z=y[i]-y1[k];
    ll[k]+=z*z;
  }
  double lik=0;
  for(k=0; k<ncl; k++){
    lik+=ll[k]/y2[k]+nc[k]*log(y2[k]);
  }
  lik *= (-0.5);
  lik -= 0.5*nn*log(6.283);
  free(y1); free(y2); free(ll); free(nc);
  return(lik);
}

struct peak *Center_peaks(int *N_peak, struct peak *peaks,
			  float *y_scr, long *x_coord, long nn)
{
  struct peak *cpeaks=malloc(*N_peak*sizeof(struct peak));
  struct peak *peak=peaks, *cpeak=cpeaks, *c1=NULL; int n=*N_peak;
  int i, B=BOX1/2; long last_j=endchr[0];
  int chr=0; int no_max=0, joined=0;
  for(i=0; i<*N_peak; i++){
    // Find maximum
    int jmax=peak->ifrag;
    long x0; float y_max;
    while(jmax>=last_j){chr++; last_j=endchr[chr];}
    if(jmax<0 || jmax>=last_j){
      printf("ERROR, peak %d has ifrag= %d chr=%d %d jmax=%d last_j=%d\n",
	     i,peak->ifrag,peak->chr,chr+1,jmax, last_j);
      no_max++;
      x0=(peak->ini+peak->end)/2; y_max=0;
    }else{
      // Find maximum
      int j=jmax+1;
      while(1){
	if(y_scr[j]>y_scr[jmax]){jmax=j;}
	j++;
	if((j>=last_j)||(x_coord[j]>=peak->end))break;
      }
      x0=x_coord[jmax]+B; // Half the way
      y_max=y_scr[jmax];
    }
    int d=x0-peak->ini, d2=peak->end-x0; if(d2 < d)d=d2;
    if(d<B)d=B;
    /*
    // Find maximum distance at which score is > DAMP*max in both dir.
    float y_thr=y_max*DAMP;
    int j1; for(j1=jmax-1; j1>=0; j1--)if(y_scr[j1]<y_thr)break; j1++;
    int j2; for(j2=jmax+1; j2<nn; j2++)if(y_scr[j2]<y_thr)break; j2--;
    int jd=jmax-j1, d2=j2-jmax; if(d2 < jd)jd=d2;*/

    // Set peak
    //cpeak->ini=x_coord[jmax-jd];
    //cpeak->end=x_coord[jmax+jd]+BOX1-1;
    long ini=x0-d; if(ini<1)ini=1;
    long end=x0+d; if(end>len_chr[peak->chr])end=len_chr[peak->chr];
    if(c1 && c1->chr==peak->chr && ini< c1->end){
      // The peaks overlap, join them
      cpeak=c1; n--; joined++;
      if(ini < cpeak->ini)cpeak->ini=ini;
      if(end > cpeak->end)cpeak->end=end;
      if(y_max>cpeak->y)cpeak->y=y_max;
    }else{
      cpeak->chr=peak->chr;
      cpeak->ini=ini;
      cpeak->end=end;
      cpeak->y=y_max;
    }
    cpeak->size=cpeak->end-cpeak->ini+1;
    cpeak->xo=(cpeak->ini+cpeak->end)/2;
    c1=cpeak;
    cpeak->next=cpeak+1;
    if(0){
      printf("peak: %d %ld %ld %d %.1f  ",
	     peak->chr, peak->ini, peak->end, peak->size, peak->y); 
      printf("cpeak: %d %ld %ld %d %.1f",
	     cpeak->chr, cpeak->ini, cpeak->end, cpeak->size, cpeak->y);
      printf("   %d %d\n", peak->size-cpeak->size, peak->xo-cpeak->xo);
    }
  next_p: 
    peak=peak->next;
    cpeak=cpeak->next;
  }
  (cpeak-1)->next=NULL;
  if(n!=*N_peak){
    printf("WARNING number of peaks reduced from %d to %d after centering\n",
	   *N_peak, n);
    printf("%d peaks without maximum and %d peaks joined\n", no_max, joined);
    *N_peak=n;
  }
  return(cpeaks);
}

struct peak *Peaks2Frags(int *N_peak, int *lpeak, float *yout_peak,
			 float *y_scr, long nn, int CLUS)
{
  int n=0, N_frag=numclus[CLUS];
  struct peak *peaks=malloc(*N_peak*sizeof(struct peak)), *peak=peaks;
  int i, m=0, k, k1=-1, nch=0; double y=0, yout=0;
  int NCH=1000, numch[NCH]; double sizech[NCH];
  for(i=0; i<NCH; i++){numch[i]=0; sizech[i]=0;}
  for(i=0; i<nn; i++){
    if(i>=endchr[nch]){ // New chromosome starts
      //if(k1==CLUS){
      if(k1>0){
	n+=Close_peak(&peak, x_coord[i-1]+BOX1, y/m, yout/m, numch, sizech);
	k1=-1;
      }
      nch++;
      if(nch>=Nchr){
	printf("ERROR too many chromosomes i=%d n=%ld ch=%d\n", i, nn, nch);
	exit(8);
      }
    } // end new chromosome
    k=lpeak[i];
    //if(k==CLUS){
    if(k1!=k){ // k=1, k1=0 : start peak
      //if(k1==CLUS){ // k=0, k1=1: close peak
      if(k1>0){ // k=0, k1=1: close peak
	n+=Close_peak(&peak, x_coord[i], y/m, yout/m, numch, sizech);
      }
      if(k>0){
	peak->ifrag=i;
	peak->ini=x_coord[i];
	peak->c=k;
	peak->chr=nch;
	y=0; yout=0; m=0;
      }
      k1=k;
    }
    if(k>0){
      y+=y_scr[i]; yout+=yout_peak[i]; m++;
    }

  }
  //if(k1==CLUS)
  if(k1>0)
    n+=Close_peak(&peak,x_coord[i-1]+BOX1,y/m,yout/m,numch,sizech);
  (peak-1)->next=NULL;
  printf("%d peaks in %d chromosomes in Peak2Frag %d\n", n, nch+1, CLUS);
  *N_peak=n;
  for(i=0; i<=nch; i++){
    if(numch[i])sizech[i]/=(float)numch[i];
  }
  return(peaks);
}

struct peak *CopyPeaks(int *N_peak, struct peak *peaks_all,
		       int N_peak_all, int CLUS)
{
  int nc=0, i;
  for(i=0; i<N_peak_all; i++){if(peaks_all[i].c==CLUS){nc++;}}
  printf("%d peaks found in cluster %d\n", nc, CLUS);
  if(nc==0){
    printf("WARNING, no peaks found in cluster %d\n", CLUS);
    return(NULL);
  }
  *N_peak=nc;

  struct peak *peaks=malloc((*N_peak)*sizeof(struct peak));
  struct peak *prev=NULL, *p=peaks;
  for(i=0; i<N_peak_all; i++){
    if(peaks_all[i].c==CLUS){
      *p=peaks_all[i];
      if(prev){prev->next=p;}
      prev=p; p++;
    }
  }
  prev->next=NULL;
  return(peaks);
}


int Read_contr_exp(char **chr_name, long *nn, long *nnch,
		   long **x_coord, float **y_contr, float **y_exper, float *r,
		   char **file_e, int nexp, char **file_c, int ncontr)
{
  int k;
  long lch1[NCHRMAX], lch2[NCHRMAX];
  char *chr_name1[NCHRMAX], *chr_name2[NCHRMAX];
  int step1=0, step2=0;

  printf("Counting chromosomes for experiment and control\n");
  for(k=0; k<NCHRMAX; k++){lch1[k]=0; lch2[k]=0;}
  int Nch1=Count_chromosomes(chr_name1, lch1, &step1, file_e, nexp);
  int Nch2=Count_chromosomes(chr_name2, lch2, &step2, file_c, ncontr);
  int Nch=Nch1; if(Nch2>Nch1)Nch=Nch2;

  int num_chr1[Nch], num_chr2[Nch];
  for(k=0; k<Nch; k++){
    num_chr1[k]=-1; num_chr2[k]=-1;
    chr_name[k]=malloc(20*sizeof(char));
  }

  if(Nch2==0){ // Control is not present
    for(k=0; k<Nch; k++){
      num_chr1[k]=k;
      strcpy(chr_name[k], chr_name1[k]);
      nnch[k]=lch1[k];
    }
  }else{ // Control is present, use it for numbering

    if(step2!=step1){
      printf("ERROR, different step\n");
      printf("Experiment: %s %d chromosomes step=%d\n",file_e[0],Nch1,step1);
      printf("Control:    %s %d chromosomes step=%d\n",file_c[0],Nch2,step2);
      exit(8);
    }
    if(Nch2 != Nch1){
      printf("WARNING, different number of chromosomes\n");
      printf("Experiment: %s %d chromosomes step=%d\n",file_e[0],Nch1,step1);
      printf("Control:    %s %d chromosomes step=%d\n",file_c[0],Nch2,step2);
    }
    if(Nch2<Nch1){
      printf("ERROR, fewer chromosomes in control than experiment\n");
      exit(8);
    }
    int j;
    for(k=0; k<Nch2; k++){
      strcpy(chr_name[k], chr_name2[k]);
      num_chr2[k]=k; 
      for(j=0; j<Nch1; j++)if(strcmp(chr_name2[k],chr_name1[j])==0)break;
      if(j==Nch1){
	printf("WARNING, chromosome %s not present in exp\n", chr_name2[k]);
	nnch[k]=lch2[k];
      }else{
	if(lch1[j]>lch2[k]){nnch[k]=lch1[j];}else{nnch[k]=lch2[k];}
	if(lch1[j]!=lch2[k]){
	  printf("WARNING, different number of lines in exp (%ld) "
		 "and contr (%ld) chrom %s\n", lch1[j], lch2[k], chr_name[k]);
	  if(abs(lch1[j]-lch2[k])>200){
	    printf("Exiting\n"); exit(8);
	  }
	}
	num_chr1[j]=k; 
      }
    }
  }
  printf("%d chromosomes found:\n", Nch);
  for(k=0; k<Nch; k++){printf("%d %d\n", num_chr1[k], num_chr2[k]);}
  if(Nch<1){printf("Too few chromosomes, exiting\n"); exit(8);}
  (*nn)=0; for(k=0; k<Nch; k++)(*nn)+=nnch[k];

  /*printf("%d chromosomes in experiment (%d files, step=%d):\n",
	 Nch1,nexp,step1);
  for(k=0; k<Nch1; k++){
    printf("%s %d %d %d\n",
	   chr_name1[k],nnch[k],lch1[k],step1*lch1[k]);
  }
  printf("%d chromosomes in control (%d files, step=%d):\n",
	 Nch2,ncontr,step2);
  for(k=0; k<Nch2; k++){
    printf("%s %d %d %d\n",
	   chr_name2[k],nnch[k],lch2[k],step2*lch2[k]);
  }
  printf("%ld lines total\n", *nn);*/

  // initialize coordinates
  *x_coord=malloc(*nn*sizeof(long));
  *y_exper=malloc(*nn*sizeof(long));
  *y_contr=NULL; long i;
  for(i=0; i<*nn; i++){(*x_coord)[i]=-1; (*y_exper)[i]=0;}
  printf("Reading chromosomes for experiment and control\n");
  Read_chromosomes(*x_coord, *y_exper, lch1, *nn, num_chr1, file_e, nexp);
  if(Nch2){
    *y_contr=malloc(*nn*sizeof(long));
    for(i=0; i<*nn; i++)(*y_contr)[i]=0;
    Read_chromosomes(*x_coord, *y_contr, lch2, *nn, num_chr2, file_c, ncontr);
    float slope, offset;
    *r=Corr_coeff(&slope,&offset,(*y_contr),(*y_exper),(int)(*nn));
    printf("Correlation between experiment and control: %.3f\n", *r);
  }

  printf("%d chromosomes in experiment (%d files, step=%d):\n",
	 Nch1,nexp,step1);
  for(k=0; k<Nch1; k++){
    printf("%s %d %d %d\n",
	   chr_name1[k],nnch[k],lch1[k],step1*lch1[k]);
  }
  printf("%d chromosomes in control (%d files, step=%d):\n",
	 Nch2,ncontr,step2);
  for(k=0; k<Nch2; k++){
    printf("%s %d %d %d\n",
	   chr_name2[k],nnch[k],lch2[k],step2*lch2[k]);
  }  
  printf("%d chromosomes, %ld lines total\n", Nch, *nn);

  return(Nch);
}


int Count_chromosomes(char **chr_name, long *lchr, int *step,
		      char **file, int nf)
{
  int Nchr=0, k;
  for(k=0; k<nf; k++){
    printf("Reading %s\n", file[k]);
    char *file_ext=Extension(file[k]);
    if(strncmp(file_ext, "wig",3)==0){
      int n=Read_chroms_wig(lchr+Nchr, chr_name+Nchr, step, file[k], NCHRMAX);
      printf("file %s %d chromosomes read\n",file[k], n);
      Nchr+=n;
    }else if(strncmp(file_ext, "gr",2)==0){
      FILE *file_in=fopen(file[k], "r"); int n=0;
      char string[1000];
      while(fgets(string, sizeof(string), file_in)!=NULL){
	if(string[0]!='#')n++;
	if(n==1)sscanf(string, "%d", step);
      }
      chr_name[Nchr]=malloc(80*sizeof(char));
      sprintf(chr_name[Nchr], "chr%d", k+1);
      lchr[Nchr]=n; Nchr++; 
    }else{
      printf("ERROR, unknown file extension %s\n", file[k]);
      exit(8);
    }
  }
  return(Nchr);
}

void Read_chromosomes(long *x, float *y, long *lch, long nn,
		      int *num_chr, char **file, int nf)
{
  int Nchr=0, k, nc=0; long n=0, nk;
  for(k=0; k<nf; k++){
    printf("Reading %s\n", file[k]);
    char *file_ext=Extension(file[k]);
    if(strncmp(file_ext, "wig",3)==0){
      Read_wig_1(&n, x+n, y+n, nn, &nc, lch+nc, num_chr, file[k], NCHRMAX);
    }else if(strncmp(file_ext,"gr",2)==0){
      if(num_chr[k]<0)continue;
      FILE *file_in=fopen(file[k], "r");
      char string[1000];
      while(fgets(string, sizeof(string), file_in)!=NULL){
	if(string[0]=='#')continue;
	sscanf(string, "%ld %f", x, y); x++; y++; n++;
      }
      fclose(file_in);
    }
  }
  if(n != nn){
    printf("ERROR, different number of lines: %ld %ld\n", n, nn);
    exit(8);
  }
}

int Get_maxima(int *lpeak, float *yout_peak, double *S, float *y_ave,
	       float *y_scr, float *y_exper, long nn,
	       float ythr, int wp_min, float *zt, float Sdy, float Sqr_w,
	       int get_y)
{
  float s=Sdy, Sd2=Sdy*Sdy;
  int N=0, Nm=0, wp2_min=2*wp_min+2; long n1=0; *S=0;
  for(int i=0; i<nn; i++){lpeak[i]=0; yout_peak[i]=0;}
  int num[NC+1]; for(int i=0; i<=NC; i++){y_ave[i]=0; num[i]=0;}

  for(int k=0; k<Nchr; k++){
    long nk=nnch[k], n2=n1+nk;
    if((Omit_chr)&&(Omit_chr[k])){n1=n2; continue;}
    if(n2>nn){printf("ERROR in Get_maxima\n"); exit(8);}
    for(int i=n1; i<n2; i++){
      if(lpeak[i])continue;
      if(y_exper[i]<=ythr)continue;
      // Verify that the peak is bell-shaped for i-wp<=j<=i+wp
      int wp=0, i1=i, i2=i;
      if(SYMMETRIC){
	while((i1 >=n1) && i2<n2){
	  if((y_scr[i1]<y_scr[--i1])||(y_scr[i2]<y_scr[++i2]))break;
	  wp++;
	}
        if(wp<wp_min)continue;
      }else{
	int w1=0, w2=0;
	while(i1 >=n1){
	  if(y_scr[i1]<y_scr[--i1]){break;} w1++;
	}
	if(w1==0)continue;
	while(i2 <n2){
	  if(y_scr[i2]<y_scr[++i2]){break;} w2++;
	}
	if(w2==0)continue;
        if((w1+w2)<wp2_min)continue;
	wp=(w1+w2)/2;
      }

      // Determine mid value of y at peak
      i1++; i2--;
      float y_max=y_scr[i], ym1=y_scr[i-wp_min]; //[i1]
      float ym2=y_scr[i+wp_min], y_min; // i2
      if(ym1<y_max){
	y_min=ym1; if(ym2<y_max && ym2>y_min){y_min=ym2;}
      }else{
	y_min=ym2;
      }
      float y_thr=y_min+Y_MAX*(y_max-y_min); // 0.5<=W<=1
      // Compute local Z score, impose it is > Z_thr 
      double y1=0, y2=0; float *yy;
      int ns=0, j=i1, ini=0; yy=y_scr+j;
      while(j>=n1){
	//if(m>=RANGE)break;
	if(ns>=RANGE)break;
	if(EXCLUDE_PEAK &&  *yy >y_thr){
	  int jj=j-1; float *yp=yy-1;
	  while(*yp>y_thr && jj>=n1){jj--; yp--;}
	  int l=j-jj;
	  if(l>=LMIN_excl || ini){yy=yp; j=jj;}
	  if(j<n1){goto end_left;}
	}
	ns++; y1+=(*yy); if(LOCAL)y2+=(*yy)*(*yy);
	if(ini){ini=0;}
	//next1:
	j--; yy--;
      }
    end_left:
      int n=ns; ns=0; j=i2; yy=y_scr+j;
      while(j<n2){
	//if(m>=RANGE)break;
	if(ns>=RANGE)break;
	if(EXCLUDE_PEAK &&  *yy >y_thr){
	  int jj=j+1, ini=0; float *yp=yy+1;
	  while(*yp>y_thr && jj<n2){jj++; yp++;}
	  int l=jj-j;
	  if(l>=LMIN_excl || ini){yy=yp; j=jj;}
	  if(j>=n2){goto end_right;}
	}
	ns++; y1+=(*yy); if(LOCAL)y2+=(*yy)*(*yy);
	//next2:
	j++; yy++;
 	if(ini){ini=0;}
      }
    end_right:
      n+=ns;
      //if(n<=m || n<2)continue; // m: size of peak n: range for computing <y>
      if(n<5)continue;

      //int j2=j-1; // where Z=thr
      y1/=n;
      if(LOCAL){
	y2=(y2-n*y1*y1)/(n-1);
	if(y2>Sd2){s=sqrt(y2);}else{s=Sdy;} // Choose the maximum
	//if(y2>0){s=pow(Sd2*y2,0.25);}else{s=Sdy;} // geometric mean
	//if(y2>Sd2){s=pow(Sd2*y2,0.25);}else{s=Sdy;} // Choose the maximum
      }
      // determine class of ymax and y1 in [1,NC1]
      j=wp-wp_min; if(j>N_WP1){j=N_WP1;}else if(j<0){j=0;}
      float t=thr_WP[j];

      double y_out=0;
      if(get_y){
	int j=i1; yy=y_scr+j; int ns1=0;
	while(j>=n1 && ns1<R_OUT){y_out+=*yy; ns1++; j--; yy--;}
	int n_out=ns1;
	j=i2; yy=y_scr+j; ns1=0;
	while(j<n2 && ns1<R_OUT){y_out+=*yy; ns1++; j++; yy++;}
	n_out+=ns1;
	y_out/=n_out;
      }else{
	y_out=y1;
      }

      int c; for(c=1; c<NC1; c++){if(y_out<ymax_c[c])break;} // y_max

      if(ZTNUM){t*=zt[c];}
      float thr=y1+t*s; // y<thr' => (y-y1)/s<thr
      if(y_scr[i]<thr){continue;}

      int m=1, j1=i, j2=i;
      if(SYMMETRIC){
	while(y_scr[j1]>=thr && y_scr[j2]>=thr){
	  m+=2; j1--; j2++;
	  if(j1<i1 || j2>i2)break; // i1 i2
	}
      }else{
	while(y_scr[j1]>=thr){
	  m++; j1--; if(j1<i1)break; // i1 
	}
	while(y_scr[j2]>=thr){
	  m++; j2++; if(j2>i2)break; //i2
	}
      }
      if(m<LMIN_peak)continue; // Exclude short fragments

      y_ave[c]+=y_out; num[c]++;
      float sc=0; m=0;
      if(lpeak[i]==0){
	if(SCORE){sc+=y_scr[i];} m++; //Z
	N++; // Otherwise, it was already counted
      }
      lpeak[i]=c;
      yout_peak[i]=y_out;
      j1=i; j2=i;
      if(SYMMETRIC){
	while(y_scr[j1]>=thr && y_scr[j2]>=thr){
	  if(lpeak[j1]==0){m++; if(SCORE){sc+=y_scr[j1];}}
	  lpeak[j1]=c; yout_peak[j1]=y_out;
	  if(lpeak[j2]==0){m++; if(SCORE){sc+=y_scr[j2];}}
	  lpeak[j2]=c; yout_peak[j2]=y_out;
	  j1--; j2++; if(j1<i1 || j2>i2)break; //if(j1<n1 || j2==n2)break; 
	}
      }else{
	while(y_scr[j1]>=thr){
	  if(lpeak[j1]==0){m++; if(SCORE){sc+=y_scr[j1];}}
	  lpeak[j1]=c; yout_peak[j1]=y_out;
	  j1--; if(j1<i1)break; //if(j1<n1)break; 
	}
	while(y_scr[j2]>=thr){
	  if(lpeak[j2]==0){m++; if(SCORE){sc+=y_scr[j2];}}
	  lpeak[j2]=c; yout_peak[j2]=y_out;
	  j2++; if(j2>i2)break; //if(j2==n2)break; 
	}
      }
      if(SCORE){(*S)+=(sc-m*y1)/s;} // (sc/m-y1)/s
      Nm+=m;
    }
    n1+=nk; // next chromosome
  }
  if(SCORE==0){*S=Nm;} // Number of bins
  else{
    // *S/=Sqr_w; // Normalize by bin width
    if(SCORE==2){(*S)*=(N/Nm);}
  }
  // Compute mean y
  for(int i=0; i<NC; i++){
    y_ave[NC]+=y_ave[i];
    num[NC]+=num[i];
    if(num[i])y_ave[i]/=num[i];
  }
  if(num[NC])y_ave[NC]/=num[NC];
  return(N);
}


void Get_threshold(float Thr, float *zt, float *y_scr, long nn, int NC)
{
  long num[NC], i, c, NC1=NC-1; for(c=0; c<NC; c++)num[c]=0;
  for(i=0; i<nn; i++){
    for(c=1; c<NC1; c++){if(y_scr[i]<ymax_c[c])break;}
    num[c]++;
  }
  int c0=1; for(c=2; c<NC; c++){if(num[c]>0 && num[c]<num[c0]){c0=c;}}

  for(c=1; c<NC; c++){
    if(num[c]==0){zt[c]=1;}
    else{zt[c]=1+PENALTY*log(num[c]/(float)num[c0]);}
    //printf("cluster %d has size %d Z= %.3g\n", c, num[c], zt[c]);
  }
}

void Ave_prof(double **ave_prof, double **dev_prof, int N_prof, int Nchr,
	      float ***yprof, long **nprof)
{
 for(int p=0; p<N_prof; p++){
   for(int chr=0; chr<Nchr; chr++){
     double y1=0, y2=0;
     long nj=nprof[p][chr];
     float *yp=yprof[p][chr];
     for(long j=0; j<nj; j++){y1+=(*yp); y2+=(*yp)*(*yp); yp++;}
     y1/=nj;
     y2=(y2/nj-y1*y1);
     ave_prof[p][chr]=y1;
     dev_prof[p][chr]=sqrt(y2);
   }
 }
}
