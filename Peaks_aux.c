#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "Peaks_aux.h"
#include "random3.h"

#define NCHAR 400
#define VERBOSE 0

static unsigned long randomgenerator(void);
static int INIRAN;
static float RANFACTOR;

extern char *Extension(char *string);
extern int Get_step_wig(char *string);
extern int Get_chr_wig(int *span, char *chr_name, char *string);
extern int Count_columns(char *string);
int ZSCORE, PEAKSIZE;

static int Overlap(struct peak peak1, struct peak peak2, float DTOL);
static void Complete_chromosome(long **xx,float **yy, short **ch,
				int xc, int step, int n, int nch, int ichr);
float Get_ave_sd(double *dev, double sum1, double sum2, double norm);

float Sum_score(long *k2, long ini, long end,
		float *z, long *x, long k, long m);
static void Print_peak(FILE *file_out, struct peak peak);


void Profile_score_ZP(double **Peak_score, FILE *file_out, 
		      long ***xprof, float ***yprof, long **nprof,
		      int *nchr_prof,
		      double **ave_prof, double **dev_prof, int N_prof,	      
		      struct peak *peak, int N_PEAK)
{

  int Nchr=100, num_chr[Nchr]; for(int i=0; i<Nchr; i++)num_chr[i]=0;
  int NP=0; struct peak *peak1=peak; Nchr=0;
  double size_1=0, size_2=0;
  while(peak1!=NULL){
    int chr=peak1->chr;
    float s=peak1->end-peak1->ini+1;
    num_chr[chr]++; if(chr>Nchr)Nchr=chr;
    size_1+=s; size_2+=s*s; NP++;
    peak1=peak1->next;
  }
  Nchr++;
  size_1=Get_ave_sd(&size_2, size_1, size_2, NP);
  char buffer[1000], buffer2[1000];
  sprintf(buffer, "%d\t%.0f\t%.0f", NP, size_1, size_2);
  sprintf(buffer2, "%d\t%.0f", NP, size_1);

  // Print distribution across chromosomes
  char str_char[Nchr*20], tmp[20];
  sprintf(str_char, "# Npeaks: %d (tot), ", NP);
  for(int i=0; i<Nchr; i++){
    sprintf(tmp, " %d (ch%d)", num_chr[i],i+1); strcat(str_char, tmp);
  }

  double ave_peak[N_prof], dev_peak[N_prof], Discriminant_score[N_prof];
  float sn=sqrt(NP);
  int warn=1, nwarn=0;
  for(int p=0; p<N_prof; p++){
    //printf("Profile= %d over %d\n", p, N_prof);
    // average score across whole genome
    double y_tot_1=0, y_tot_2=0; long n_tot=0;
    for(int chr=0; chr<Nchr; chr++){
      long nj=nprof[p][chr];
      n_tot+=nj;
      y_tot_1+=nj*ave_prof[p][chr];
      y_tot_2+=nj*(ave_prof[p][chr]*ave_prof[p][chr]+
		   dev_prof[p][chr]*dev_prof[p][chr]);
    }

    double *Ps=NULL; if(Peak_score){Ps=Peak_score[p];}
    double sum_peak_1=0, sum_peak_2=0;
    double y_peak_1=0, y_peak_2=0; long l_peak=0;
    int chr=-1, nchr=nchr_prof[p], j=0, nj=0, N_peak=0;
    long *xp; float *yp;
    struct peak *peak1=peak;
    while(peak1!=NULL){
      if(peak1->chr !=chr){
	if(peak1->chr >= nchr)break;
	chr=peak1->chr; j=0;
	nj=nprof[p][chr];
	xp=xprof[p][chr];
	yp=yprof[p][chr];
      }
      while((j<nj)&&(*xp<peak1->ini)){  //ini
	j++; xp++; yp++;
      }
      // average across peak
      double y_peak=0; int norm_peak=0;
      while((j<nj)&&(*xp<=peak1->end)){ //end
	y_peak+=(*yp);
	sum_peak_2+=(*yp)*(*yp);
	norm_peak++;
	j++; xp++; yp++;
      }
      sum_peak_1+=y_peak;
      l_peak+=norm_peak;
      
      if(norm_peak){y_peak/=norm_peak;}
      y_peak_1+=y_peak;
      y_peak_2+=y_peak*y_peak;
      if(Ps && N_peak<N_PEAK){
	if(ZSCORE){Ps[N_peak]=(y_peak-ave_prof[p][chr])/dev_prof[p][chr];}
	else{Ps[N_peak]=y_peak;}
      }

      peak1=peak1->next;
      N_peak++;
      if((N_peak==N_PEAK)&&(peak1)){
	nwarn++;
	if(warn){
	  printf("WARNING, too many peaks in profile score\n");
	  printf("p[%d]: %d %ld %ld\n", N_peak, peak1->chr,
		 peak1->ini, peak1->end);
	  warn=0;
	}
      }
    }
    
    // Mean score of no-peaks
    double y_nopeak_1=y_tot_1-sum_peak_1;
    double y_nopeak_2=y_tot_2-sum_peak_2;
    long n_nopeak=n_tot-l_peak;
    y_nopeak_1=Get_ave_sd(&y_nopeak_2, y_nopeak_1, y_nopeak_2, n_nopeak);

    // Mean score of peaks
    if(0){y_peak_1=Get_ave_sd(&y_peak_2, y_peak_1, y_peak_2, N_peak);}
    y_peak_1=Get_ave_sd(&y_peak_2, sum_peak_1, sum_peak_2, l_peak);

    // Mean score of all genome
    y_tot_1=Get_ave_sd(&y_tot_2, y_tot_1, y_tot_2, n_tot);

    ave_peak[p]=(y_peak_1-y_tot_1)/y_tot_2; //*sqrt((float)N_peak)
    dev_peak[p]=y_tot_2/sn;
    // Normalize by standard dev. of peaks. It has drawbacks.
    Discriminant_score[p]=
      sn*(y_peak_1-y_nopeak_1)/y_tot_2;
      //sqrt(y_peak_2*y_peak_2+y_nopeak_2*y_nopeak_2);

    char tmp[90];
    sprintf(tmp, "\t%.3g\t%.2g\t%.2g",
	    ave_peak[p], dev_peak[p], Discriminant_score[p]);
    strcat(buffer, tmp);
    sprintf(tmp, "\t%.3g",ave_peak[p]); //Discriminant_score[p]
    strcat(buffer2, tmp);

  } // end p
  fprintf(file_out, "%s\n", buffer);
  printf("%s\n", buffer2);
  //fprintf(file_out,"%s\n", str_char);
  printf("%s\n", str_char);
  if(nwarn)printf("WARNING, %.0f extra peaks in profile score\n",
		  nwarn/N_prof);

}

void Print_peak_score(struct peak *peak, double **Peak_score,
		      int N_prof, char **name_prof, char *nameout)
{
  //int MM = 3+2*N_prof;
  int M0=1; //2; 1=size 2=dist
  int MM = M0+N_prof, M=MM, i, j;  

  // Print to file
  FILE *file_out=Open_file_w(nameout);
  fprintf(file_out, "### 1=size ");
  if(M0==2){fprintf(file_out, " 2=min_dist/kb");}
  for(i=0; i<N_prof; i++)fprintf(file_out," %d=%s",i+M0+1, name_prof[i]);
  fprintf(file_out, " chr x_peak\n");
  fprintf(file_out, "### "); for(i=0; i<M; i++)fprintf(file_out, " 1");
  fprintf(file_out, " 0 0\n");
  fprintf(file_out, "### "); for(i=0; i<M; i++)fprintf(file_out, " 1");
  fprintf(file_out, " 0 0\n");

  int ipeak=0;
  struct peak *peak1=NULL, *peak2=peak, *peak3=peak->next;
  while(peak2!=NULL){
    peak3=peak2->next;
    if(M0==2){
      float d, d1, d2;
      if((peak1!=NULL)&&(peak1->chr==peak2->chr)){d1=peak2->ini-peak1->end;}
      else{d1=1000000000;}
      if((peak3!=NULL)&&(peak3->chr==peak2->chr)){d2=peak3->ini-peak2->end;}
      else{d2=1000000000;}
      if(d1<d2){d=d1;}else{d=d2;}
      i=1;
      fprintf(file_out, "%.2g ",d/1000);
    }
    i=0;
    float s=peak2->size, w=s;
    fprintf(file_out, "%.0f ", s);
    if(peak2->size>SIZEMAX)fprintf(file_out, "#");
    for(j=0; j<N_prof; j++){
      s=Peak_score[j][ipeak]; i++;
      fprintf(file_out, "%.4g ", s);
    }
    fprintf(file_out, " %d %.0f\n",
	    peak2->chr, 0.5*(peak2->end+peak2->ini));
    peak1=peak2; peak2=peak2->next; ipeak++;
    if(peak2==NULL)break;
  }
  fclose(file_out);
}

struct peak *Read_peaks_ZP(int *npeak, char *file, int shift)
{
  if(file==NULL)return(0);
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  // Count lines and allocate
  int n=0; char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue; n++;
  }
  fclose(file_in);
  if(n==0)return(NULL); *npeak=n;

  // Read
  struct peak *peak=malloc(n*sizeof(struct peak)), *peaki=peak; 
  int i=0, s=0, mult=0;
  file_in=fopen(file, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue; char *ptr=string;
    if((string[0]=='c')||(string[0]=='C'))ptr+=3;
    s=sscanf(ptr, "%d%ld%ld%d",
	     &(peaki->chr), &(peaki->ini), &(peaki->end), &mult);
    if(mult>0){peaki->mult=mult;}else{peaki->mult=1;}
    peaki->chr += shift;
    //Set_peak(peaki, 1);
    peaki->size=peaki->end-peaki->ini+1;
    peaki->xo=(peaki->end+peaki->ini)/2;
    peaki->x1=peaki->xo-PEAKSIZE;
    peaki->x2=peaki->xo+PEAKSIZE;
    peaki->match=0;
    peaki->next=peaki+1;
    i++;  peaki++;
  }
  peak[n-1].next=NULL;
  fclose(file_in);
  return(peak); 
}

struct peak *Extract_peaks_ZP(struct peak *peak_old, int N_peak, int Nchr,
			      long *len_chr)
{
  unsigned long iran;
  if(INIRAN==0){
    INIRAN=1;
    iran=randomgenerator();
    InitRandom( (RANDOMTYPE)iran);
    RANFACTOR=pow(12.0,1/3.0);
  }
 
  struct peak peak[N_peak];
  struct peak chr_peak[Nchr]; int p_chr[Nchr], ichr=0, i;
  for(i=0; i<Nchr; i++){chr_peak[i].next=NULL; p_chr[i]=0;}

  for(i=0; i<N_peak; i++){
    if(peak_old[i].chr!=ichr)ichr=peak_old[i].chr;
    p_chr[ichr]++;
    struct peak *peaki=peak+i;
    peaki->chr=ichr;
    peaki->size=peak_old[i].size;
    peaki->ini=RandomFloating()*len_chr[ichr];
    peaki->end=peaki->ini+peaki->size-1;
    // Sort
    struct peak *peakj1=chr_peak+ichr, *peakj2=peakj1->next;
    while(peakj2!=NULL){
      if(peaki->ini<peakj2->ini){
	peakj1->next=peaki; peaki->next=peakj2; break;
      }
      peakj1=peakj2; peakj2=peakj2->next;
    }
    if(peakj2==NULL){peakj1->next=peaki; peaki->next=peakj2;}
  }
  // Copy
  struct peak *peaknew=malloc(N_peak*sizeof(struct peak));
  struct peak *p1=peaknew;
  int n=0;
  for(ichr=0; ichr<Nchr; ichr++){
    struct peak *peaki=chr_peak[ichr].next;
    while(peaki!=NULL){
      *p1=*peaki; p1->next=p1+1;
      peaki=peaki->next; p1++; n++;
    }
    //printf("Chromosome %d %d peaks extracted\n", ichr+1, p_chr[ichr]);
  }
  //printf("%d peaks extracted over %d\n", n, N_peak);
  (p1-1)->next=NULL;
  return(peaknew);
}

unsigned long randomgenerator(void){
     
     unsigned long tm;
     time_t seconds;
     
     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   // waiting time equal 1 second 
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     
     return((unsigned long) (rand()));

}

void Count_matches_ZP(int *n_match1, int *n_match2, float DTOL,
		      struct peak *peak1, int N_peak1,
		      struct peak *peak2, int N_peak2)
{
  int i1, i2=0, match=0;
  struct peak *peak1_p=peak1, *peak2_p=peak2;
  for(i1=0; i1<N_peak1; i1++)peak1[i1].match=0;
  for(i2=0; i2<N_peak2; i2++)peak2[i2].match=0;
  while(1){
    while((peak1_p->chr<peak2_p->chr)||
	  ((peak1_p->chr==peak2_p->chr)&&(peak1_p->end+DTOL<peak2_p->ini))){
      peak1_p=peak1_p->next; if(peak1_p==NULL)goto count;
    }
    while((peak2_p->chr<peak1_p->chr)||
	  ((peak2_p->chr==peak1_p->chr)&&(peak2_p->end+DTOL<peak1_p->ini))){
      peak2_p=peak2_p->next; if(peak2_p==NULL)goto count;
    }
    if(Overlap(*peak1_p, *peak2_p, DTOL)){
      peak1_p->match=1; peak2_p->match=1;
      if(peak1_p->end<peak2_p->end){
	peak1_p=peak1_p->next; if(peak1_p==NULL)goto count;
      }else{
	peak2_p=peak2_p->next; if(peak2_p==NULL)goto count;
      }
    }
  }
 count:
  (*n_match1)=0; struct peak *p=peak1;
  while(p!=NULL){if(p->match)(*n_match1)++; p=p->next;}
  (*n_match2)=0; p=peak2;
  while(p!=NULL){if(p->match)(*n_match2)++; p=p->next;}
}

int Overlap(struct peak peak1, struct peak peak2, float DTOL)
{
  if(peak1.chr!=peak2.chr)return(0);
  if((peak2.ini>peak1.end+DTOL)||(peak1.ini>peak2.end+DTOL))return(0);
  //iif((peak2.x1>peak1.x2+DTOL)||(peak1.x1>peak2.x2+DTOL))return(0);
  return(1);
}


int Read_wig(long *nn, long **x, float **y, short **chr,
	     char ***chr_name, char *file)
{
  long n=0; char string[1000];
  FILE *file_in=fopen(file, "r");
  int WIG=0, step=0, ncol=1, i;
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  if(strncmp(Extension(file), "wig", 3)!=0){
    printf("ERROR, file %s must be a wig file\n", file); exit(8);
  }
  // Count lines and allocate
  int nchr=0;
  fgets(string, sizeof(string), file_in); // Discard first line
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "track", 5)==0)continue;
    if(strncmp(string, "fixed", 5)==0){
      step=Get_step_wig(string); nchr++; WIG=1; continue;
    }
    n++; if(n==2)ncol=Count_columns(string);
  }
  fclose(file_in);
  (*nn)=n;
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));
  (*chr)=malloc(n*sizeof(short));
  (*chr_name)=malloc(nchr*sizeof(char *));
  for(i=0; i<nchr; i++)(*chr_name)[i]=malloc(20*sizeof(char *));

  // Report
  printf("Reading file %s ", file);
  if(WIG)printf("(wig format), step= %d ", step);
  printf("%d columns, %d lines %d chromosomes\n", ncol, n, nchr);
  if((ncol==1)&&(step==0)){
    printf("ERROR, only one column but step not reported\n"); exit(8);
  }

  // Read
  n=0;
  long *xx=*x, xxcount=1, xxn, xxo=0;
  float *yy=*y; short *ch=*chr;
  char word[100];
  file_in=fopen(file, "r");
  int ichr=-1;
  if(WIG){
    //fgets(string, sizeof(string), file_in);
    fgets(string, sizeof(string), file_in);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "fixed", 5)==0){
      ichr++; int span=0;
      Get_chr_wig(&span, (*chr_name)[ichr], string);
      xxcount=1; continue;
    }
    *ch=ichr;
    if(ncol==1){
      sscanf(string, "%f", yy);
      *xx=xxcount; xxcount+=step; 
    }else{
      sscanf(string, "%ld%s", &xxn, word);
      if(word[0]=='n'){*yy=0;}
      else{sscanf(word, "%f", yy);}
      *xx=(xxn+xxo+1)/2; xxo=xxn;
    }
    n++; xx++; yy++; ch++;
  }
  fclose(file_in);
  printf("x[0]= %ld x_last=%ld\n", (*x)[0], *(xx-1));
  return(nchr);
}

int Read_wig_1(long *n, long *x, float *y, long nn, 
	       int *nc, long *lc, int *num_chr, char *file, int NCMAX)
{
  long *xx=x, xxcount=0, xxi=0, xxn;
  char string[1000], word[80];
  FILE *file_in=fopen(file, "r");
  int step=0, ncol=-1, i, read=0;
  float *yy=y; *n=0;
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  if(strncmp(Extension(file), "wig", 3)!=0){
    printf("ERROR, file %s must be a wig file\n", file); exit(8);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "track", 5)==0)continue;
    if(strncmp(string, "fixed", 5)==0){
      step=Get_step_wig(string);
      //if(*nc){lc[num_chr[(*nc)-1]]=*n;}
      //if((*nc)>1){lc[num_chr[(*nc)-1]]-=lc[num_chr[(*nc)-2]];}
      if(*nc){lc[(*nc)-1]=*n;}
      if((*nc)>1){lc[(*nc)-1]-=lc[(*nc)-2];}

      (*nc)++;
      if(*nc> NCMAX){
	printf("ERROR, > %d chromosomes, exiting\n", NCMAX); exit(8);
      }
      xxcount=step; xxi=0; read=1;
      continue;
    }
    if(read==0)continue;
    if(ncol<0)ncol=Count_columns(string);
    if(ncol==1){
      sscanf(string, "%f", yy);
      *xx=xxcount; xxcount+=step; 
    }else{
      sscanf(string, "%ld%s", &xxn, word);
      if(word[0]=='n'){*yy=0;}
      else{sscanf(word, "%f", yy);}
      *xx=(xxn+xxi+1)/2; xxi=xxn;
    }
    (*n)++; xx++; yy++; 
    if(*n > nn){
      printf("ERROR, too many lines > %ld\n", nn);
      for(int i=0; i<*nc; i++){printf("chr%d: %ld\n", i+1, lc[i]);}
      exit(8);
    }
  }
  fclose(file_in);
  //if(*nc){lc[num_chr[(*nc)-1]]=*n;}
  //if((*nc)>1){lc[num_chr[(*nc)-1]]-=lc[num_chr[(*nc)-2]];}
  if(*nc){lc[(*nc)-1]=*n;}
  if((*nc)>1){lc[(*nc)-1]-=lc[(*nc)-2];}
  printf("x[0]= %ld x_last=%ld\n", x[0], *(xx-1));
  return(step);
}


int Read_wig_2(long **x, float **y, short **chr,
	       char *file, int *nch, int Nchr)
{
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  if(strncmp(Extension(file), "wig", 3)!=0){
    printf("ERROR, file %s must be a wig file\n", file); exit(8);
  }

  long n=0; int i;
  for(i=0; i<Nchr; i++)n+=nch[i];
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));
  (*chr)=malloc(n*sizeof(short));

  char string[1000];
  int WIG=0, step=0, ncol=1;
  // Read
  n=0;
  long *xx=*x, xxcount=1, xxn, xxo=0;
  float *yy=*y; short *ch=*chr;
  char word[100];
  file_in=fopen(file, "r");
  int ichr=-1;
  fgets(string, sizeof(string), file_in); // Eliminate first line
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(strncmp(string, "track", 5)==0)continue;
    if(strncmp(string, "fixed", 5)==0){
      step=Get_step_wig(string);
      if(ichr>=0)
	Complete_chromosome(&xx,&yy,&ch,xxcount,step,n,nch[ichr],ichr);
      ichr++; n=0; xxcount=1; WIG=1;
      continue;
    }
    *ch=ichr;
    if(ncol==1){
      sscanf(string, "%f", yy);
      *xx=xxcount; xxcount+=step; 
    }else{
      sscanf(string, "%ld%s", &xxn, word);
      if(word[0]=='n'){*yy=0;}
      else{sscanf(word, "%f", yy);}
      *xx=(xxn+xxo+1)/2; xxo=xxn;
    }
    n++; xx++; yy++; ch++;
  }
  Complete_chromosome(&xx,&yy,&ch,xxcount,step,n,nch[ichr],ichr);
  fclose(file_in);
  printf("x[0]= %ld x_last=%ld\n", (*x)[0], *(xx-1));
  return(ichr+1);
}

void Complete_chromosome(long **xx,float **yy, short **ch,
			 int xc, int step, int n, int nch, int ichr)
{
  int i, xxcount=xc;
  for(i=n; i<nch; i++){
    **ch=ichr; **yy=0; **xx=xxcount; xxcount+=step;
    (*ch)++; (*xx)++; (*yy)++;
  }
}

int Read_chroms_wig(long *nch, char **chr_name, int *step,
		    char *file, int Ncmax)
{
  char string[1000];
  FILE *file_in=fopen(file, "r");
  int WIG=0, ncol=1, i;
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  if(strncmp(Extension(file), "wig", 3)!=0){
    printf("ERROR, file %s must be a wig file\n", file); exit(8);
  }
  // Count lines per chromosome
  int nchr=-1; long n=0;
  //(*chr_name)=malloc(Ncmax*sizeof(char *));
  //fgets(string, sizeof(string), file_in); // Discard first line
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(n==0 && strncmp(string, "track", 5)==0)continue;
    if(strncmp(string, "fixed", 5)==0){
      if(nchr >= Ncmax){
	printf("ERROR in %s, too many chromosomes > %d\n",
	       file, Ncmax); exit(8);
      }
      if(nchr >= 0)nch[nchr]=n; n=0; nchr++;
      *step=Get_step_wig(string); int span=0; 
      chr_name[nchr]=malloc(20*sizeof(char *));
      Get_chr_wig(&span, chr_name[nchr], string);
      WIG=1; continue;
    }
    n++; if(n==2)ncol=Count_columns(string);
  }
  nch[nchr]=n; nchr++;
  fclose(file_in);

  // Report
  printf("Reading file %s ", file);
  if(WIG)printf("(wig format), step= %d ", *step);
  printf("%d columns, %d chromosomes\n", ncol, nchr);
  if((ncol==1)&&(step==0)){
    printf("ERROR, only one column but step not reported\n"); exit(8);
  }
  return(nchr);
}

void Set_peak(struct peak *peak, int next){
  peak->size=peak->end-peak->ini+1;
  peak->xo=(peak->end+peak->ini)/2;
  peak->x1=peak->xo-PEAKSIZE;
  peak->x2=peak->xo+PEAKSIZE;
  peak->match=0;
  if((next)&&(peak->next==NULL))peak->next=peak+1;
}


float Get_ave_sd(double *dev, double sum1, double sum2, double norm)
{
  float ave=sum1;
  if(norm){ave/=norm;}
  if(norm==1 || norm==0){
    *dev=fabs(ave);
  }else{
    *dev=sqrt((sum2-norm*ave*ave)/(norm-1));
  }
  return(ave);
}

int Make_cluster_check(struct peak *peak, int N_peak,
		       float DCLUST, long *xx, float *zz, long nn,
		       int check, float Thr)
{
  long *x=xx, k=0, k2=0;
  float *z=zz, d1, zj; int i, n=N_peak;
  //for(i=0; i<N_peak; i++)peak[i].joined=0;
  struct peak *peak1=peak, *peak2=peak1->next;
  while(peak2!=NULL){
    d1=peak2->ini-peak1->end;
    if(d1<DCLUST){
      // join 2 and 1
      if((k2>k)&&(k2<peak1->ini))k=k2; 
      zj=Sum_score(&k2, peak1->ini, peak2->end, z, x, k, nn);
      if((check==0)||(zj>Thr)){
	Join(peak1, peak2);
	peak1->y=zj; n--;
	goto next;
      }
    }
    peak1=peak2; // Not joined
  next:
    peak2=peak2->next;
  }
  return(n);
}

float Sum_score(long *k2, long ini, long end,
		float *z, long *x, long k, long m)
{
  float zj=0; int n=0; long i;
  for(i=k; i<m; i++){
    if(x[i]>end)break;
    if(x[i]>ini)zj+=z[i]; n++;
  }
  *k2=i;
  return(zj);
}

// Printing

void Print_peak(FILE *file_out, struct peak peak){
  //if(peak.size>SIZEMAX)fprintf(file_out, "#");
  fprintf(file_out, "%d\t%ld\t%ld\t%.3g",
	  peak.chr+1, peak.ini, peak.end, peak.y);
  if(peak.yout){fprintf(file_out, "\t%.3g", peak.yout);}
  fprintf(file_out, "\n");
}

void Print_Peaks(struct peak *peak, char *nameout)
{
  FILE *file_out=Open_file_w(nameout);
  //fprintf(file_out, "%s", parameters);
  struct peak *peak1=peak; int Np=0;
  while(peak1!=NULL){
    Print_peak(file_out, *peak1); Np++; peak1=peak1->next;
  }
  fclose(file_out);
  printf("Printing %d peaks in %s\n", Np, nameout);
}

void Print_Peaks_nomatch(struct peak *peak_old, int N_peak_old, char *nameout)
{ 
  int i, Np=0;
  FILE *file_out=Open_file_w(nameout);
  for(i=0; i<N_peak_old; i++){
    if(peak_old[i].match==0){
      Print_peak(file_out, peak_old[i]); Np++;
    }
  }
  fclose(file_out);
  printf("Printing %d previous not confirmed peaks in %s\n",Np, nameout);
}

void Print_Peaks_new(struct peak *peak, char *nameout, char *nameold)
{
  FILE *file_out=Open_file_w(nameold);
  FILE *file_new=Open_file_w(nameout);
  struct peak *peak1=peak; int Np=0, Np_new=0;
  while(peak1!=NULL){
    if(peak1->match==0){
      Print_peak(file_new, *peak1); Np_new++;
    }else{
      Print_peak(file_out, *peak1); Np++;
    }
    peak1=peak1->next;
  }
  fclose(file_out); fclose(file_new); 
  printf("Printing %d previous confirmed peaks in %s\n",
	 Np, nameold);
  printf("Printing %d new peaks in %s\n",
	 Np_new, nameout);
}

int Select_peaks_match(struct peak *peak2, struct peak *peak, int N, int match)
{
  int N_match=0;
  struct peak *peak1=peak;
  struct peak *p2=peak2;
  while(peak1!=NULL){
    if(peak1->match==match){
      *p2=*peak1; p2->next=p2+1;
      p2++; N_match++;
    }
    peak1=peak1->next;
  }
  (p2-1)->next=NULL;
  return(N_match);
}


/*

void Intersect_peaks(struct peak *peak,
		     struct peak *peak1,
		     struct peak *peak2);

int Read_file(long **x, float **y, char *file){
  long n=0; char string[1000];
  FILE *file_in=fopen(file, "r");
  int WIG=0, step=0, ncol=1;
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  if(strncmp(Extension(file), "wig", 3)==0){
    while(step==0){
      fgets(string, sizeof(string), file_in);
      if(strncmp(string, "fixed", 5)==0){
	step=Get_step_wig(string); WIG=1;
      }
    }
    if(step==0){
      printf("ERROR, no step found in wig file %s\n", file); exit(8);
    }
  }
  // Count lines and allocate
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue; n++;
    if(n==2)ncol=Count_columns(string);
  }
  fclose(file_in);
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));

  // Report
  printf("Reading file %s ", file);
  if(WIG)printf("(wig format), step= %d ", step);
  printf("%d columns, %d lines\n", ncol, n);
  if((ncol==1)&&(step==0)){
    printf("ERROR, only one column but step not reported\n"); exit(8);
  }

  // Read
  n=0;
  long *xx=*x, xxcount=(1-step)/2, xxn, xxo=0;
  float *yy=*y;
  char word[100];
  file_in=fopen(file, "r");
  if(WIG){
    //fgets(string, sizeof(string), file_in);
    fgets(string, sizeof(string), file_in);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    if(ncol==1){
      sscanf(string, "%f", yy);
      xxcount+=step; *xx=xxcount;
    }else{
      sscanf(string, "%ld%s", &xxn, word);
      if(word[0]=='n'){*yy=0;}
      else{sscanf(word, "%f", yy);}
      *xx=(xxn+xxo+1)/2; xxo=xxn;
    }
    n++; xx++; yy++;
  }
  fclose(file_in);
  printf("x[0]= %ld x_last=%ld\n", (*x)[0], *(xx-1));
  return(n);
}

void Read_name(int *n, char **file, char *PATH, char *string, int nmax)
{
  char name[NCHAR];
  if(*n>nmax){
    printf("ERROR, too many input files, maximum %d\n", nmax);
    exit(8);
  }
  sscanf(string, "%s", name);
  file[*n]=malloc(NCHAR*sizeof(char));
  sprintf(file[*n], "%s%s", PATH, name);
  //printf("File to read: %s\n", file[*n]);
  (*n)++;
}

void Output_name(char *nameout, char *namein, char *ext){
  char *s=namein, *si=namein;
  while((*s!='\0')&&(*s!='\n')){
    if(*s=='/'){si=s+1;}
    else if(*s==' '){*s='-';}
    s++;
  }
  sprintf(nameout, "%s%s",si, ext);
}

void Shift_coord(long **x_new, long *x_old, long n)
{
  int i;
  (*x_new)=malloc(n*sizeof(long));
  long x0=0, *x1=*x_new, *x2=x_old;
  for(i=0; i<n; i++){
    *x1=0.5*(x0+*x2); x0=*x2; x1++; x2++;
  }
}

FILE *Open_file_w(char *name)
{
  FILE *file_out=fopen(name, "w");
  printf("\nWriting %s\n", name);
  return(file_out);
}


int Unify_peaks(struct peak *peak,
		struct peak *peak1, int N_peak1,
		struct peak *peak2, int N_peak2,
		float DTOL)
{
  int Np=0, i, match=0;
  struct peak *peak1_p=peak1, *peak2_p=peak2, *peak_p=peak-1;
  for(i=0; i<(N_peak1+N_peak2); i++){
    peak[i].ini=-1; peak[i].end=-1; peak[i].match=0;
  }
  for(i=0; i<N_peak1; i++)peak1[i].match=0;
  for(i=0; i<N_peak2; i++)peak2[i].match=0;
  printf("Mult 1: ");
  for(i=0; i<60; i++)printf("%d", peak1[i].mult); printf("\n");
  printf("Mult 2: ");
  for(i=0; i<60; i++)printf("%d", peak2[i].mult); printf("\n");

  while(1){
    while((peak1_p->chr<peak2_p->chr)||
	  ((peak1_p->chr==peak2_p->chr)&&(peak1_p->end+DTOL<peak2_p->ini))){
      if(peak1_p->match==0){peak_p++; *peak_p=*peak1_p; Np++;}
      peak1_p=peak1_p->next; if(peak1_p==NULL)goto count_2;
    }
    while((peak2_p->chr<peak1_p->chr)||
	  ((peak2_p->chr==peak1_p->chr)&&(peak2_p->end+DTOL<peak1_p->ini))){
      if(peak2_p->match==0){peak_p++; *peak_p=*peak2_p; Np++;}
      peak2_p=peak2_p->next; if(peak2_p==NULL)goto count_2;
    }
    if(Overlap(*peak1_p, *peak2_p, DTOL)){
      if((peak1_p->match==0)&&(peak2_p->match==0)){peak_p++; Np++;}
      Intersect_peaks(peak_p, peak1_p, peak2_p);
      peak_p->match=1; peak1_p->match=1; peak2_p->match=1;
      if(peak1_p->end<peak2_p->end){
	peak1_p=peak1_p->next; if(peak1_p==NULL)goto count_2;
      }else{
	peak2_p=peak2_p->next; if(peak2_p==NULL)goto count_2;
      }
    }
  }
 count_2:
  while(peak1_p!=NULL){
    if(peak1_p->match==0){peak_p++; *peak_p=*peak1_p; Np++;}
    peak1_p=peak1_p->next;
  }
  while(peak2_p!=NULL){
    if(peak2_p->match==0){peak_p++; *peak_p=*peak2_p; Np++;}
    peak2_p=peak2_p->next;
  }
  for(i=0; i<Np; i++)peak[i].next=peak+i+1;
  peak[Np-1].next=NULL;
  return(Np);
}

void Intersect_peaks(struct peak *peak,
		     struct peak *peak1,
		     struct peak *peak2)
{
  peak->chr=peak1->chr;
  if(peak1->ini<peak2->ini){
    if((peak->ini<0)||(peak->ini>peak1->ini))peak->ini=peak1->ini;
  }else{
    if((peak->ini<0)||(peak->ini>peak2->ini))peak->ini=peak2->ini;
  }
  if(peak1->end>peak2->end){
    if((peak->end<0)||(peak->end<peak1->end))peak->end=peak1->end;
  }else{
    if((peak->end<0)||(peak->end<peak2->end))peak->end=peak2->end;
  }
  if(peak->match==0){peak->mult=peak1->mult+peak2->mult;}
}

float Dist_peaks(struct peak *peak1, int N_peak1,
		 struct peak *peak2, int N_peak2)
{
  double d_sum=0; float d1, d2;
  struct peak *peak=peak1, *peak_l=peak2, *peak_r;
  while(peak!=NULL){
    while((peak_l!=NULL)&&(peak_l->chr<peak->chr))peak_l=peak_l->next;
    if(peak_l!=NULL)peak_r=peak_l->next;
    while((peak_l!=NULL)&&(peak_l->chr==peak->chr)&&
	  (peak_r!=NULL)&&(peak_r->chr==peak->chr)&&(peak_r->xo < peak->xo)){
      peak_l=peak_r; peak_r=peak_l->next;
    }
    if((peak_l!=NULL)&&(peak_l->chr==peak->chr))
      {d1=fabs(peak_l->xo-peak->xo);}else{d1=10000000;}
    if((peak_r!=NULL)&&(peak_r->chr==peak->chr))
      {d2=fabs(peak_r->xo-peak->xo);}else{d2=10000000;}
    if(d1<d2){d_sum+=d1;}else{d_sum+=d2;}
    peak=peak->next;
  }
  return(d_sum/N_peak1);
}

int Join(struct peak *peak1, struct peak *peak2) 
{
  //peak1->joined=1;
  peak1->end=peak2->end;
  peak1->next=peak2->next;
  //if(peak1->norm){
    //peak1->norm+=peak2->norm;
    //peak1->xsum+=peak2->xsum;
    //peak1->xo=peak1->xsum/peak1->norm;
    //if((peak1->xo>=peak1->end)||(peak1->xo<=peak1->ini)){
    //printf("ERROR in join, wrong mid point %.0f (%.0f-%.0f)\n",
    //peak1->xo, peak1->ini, peak1->end); exit(8);
    //}
    //}else{
  peak1->xo=(peak1->end+peak1->ini)/2.;
  //}
  peak1->x1=peak1->xo-PEAKSIZE;
  peak1->x2=peak1->xo+PEAKSIZE;
  peak1->size=peak1->end-peak1->ini+1;
  return(1);
}


void Make_cluster_4(struct peak *peak, int N_peak, float DCLUST)
{
  int i;
  //for(i=0; i<N_peak; i++)peak[i].joined=0;
  int njoin=1, inichr=1;
  float d1, d2, d3;
  while(njoin){
    njoin=0;
    struct peak *peak1=NULL, *peak2=NULL, *peak3=NULL, *peak4=peak;
    while(peak4!=NULL){
      if(peak1==NULL){
	inichr=1; peak1=peak;
      }else if(peak4->chr!=peak3->chr){
	inichr=1; peak1=peak4;
	// Join at end of chromosome
	if((d2<DCLUST)&&(d2<d1)){
	  njoin+=Join(peak2, peak3);
	}
      }
      if(inichr){
	inichr=0;
	peak2=peak1->next; d1=peak2->xo-peak1->xo;
	peak3=peak2->next; d2=peak3->xo-peak2->xo;
	if((d1<DCLUST)&&(d1<d2)){
	  njoin+=Join(peak2, peak3);
	  peak2=peak1->next; d1=peak2->xo-peak1->xo;
	  peak3=peak2->next; d2=peak3->xo-peak2->xo;
	}
	peak4=peak3->next; d3=peak4->xo-peak3->xo;
      }
      if((d2<DCLUST)&&(d2<d1)&&(d2<d3)){
	// join 2 and 3
	njoin+=Join(peak2, peak3);
	peak3=peak2->next; d2=peak3->xo-peak2->xo;
	peak4=peak3->next; if(peak4==NULL)break;
	d3=peak4->xo-peak3->xo;
	
      }
      peak1=peak2; peak2=peak3; peak3=peak4; peak4=peak3->next;
      if(peak4==NULL)break;
      d1=d2; d2=d3; d3=peak4->xo-peak3->xo;
    }
  }
}

int Make_cluster_1(struct peak *peak, int N_peak, float DCLUST)
{
  int i, n=2;   float d1;
  //for(i=0; i<N_peak; i++)peak[i].joined=0;
  struct peak *peak1=peak, *peak2=peak1->next;
  while(peak2!=NULL){
    if(peak2->chr!=peak1->chr){
      peak1=peak2; peak2=peak2->next;
      if(peak2==NULL)break;
    }
    d1=peak2->ini-peak1->end;
    if(d1<DCLUST){
      // join 2 and 1
      Join(peak1, peak2);
    }else{
      peak1=peak2; n++;
    }
    peak2=peak2->next;
  }
  return(n);
}


int Remove_fragments(struct peak *peak, int N_peak, float Size_min)
{
  int n=N_peak, i;
  struct peak *peak1=peak, *peak2;
  while((peak1!=NULL)&&(peak1->size<Size_min)){
    //peak1->removed=1;
    peak1=peak1->next; n--;
  }
  if(peak1==NULL)return(0);
  if(peak1!=peak){*peak=*peak1; peak1=peak;}
  peak2=peak1->next;
  while(peak2!=NULL){
    if(peak2->size<Size_min){
      //peak2->removed=1; 
      peak1->next=peak2->next; n--;
    }else{
      //peak2->removed=0;
      peak1=peak2;
    }
    peak2=peak2->next;
  }
  printf("Removing fragments with size < %.0f ", Size_min);
  printf("obtaining %d from %d peaks\n", n, N_peak);
  return(n);
}

int Get_chr_wig(char *chr_name, char *string){
  char *s=string; int i, l=sizeof(string), step;
  while(*s!='\n'){
    if(strncmp(s, "chrom=", 6)==0){
      sscanf(s+6, "%s", chr_name); return(0);
    }
    s++;
  }
  printf("ERROR, chromosome name not found in wig file, %s", string);
  exit(8);
}

int Get_step_wig(char *string){
  char *s=string; int i, l=sizeof(string), step;
  while(*s!='\n'){
    if((strncmp(s, "step=", 5)==0)||
       (strncmp(s, "Step=", 5)==0)||
       (strncmp(s, "span=", 5)==0)){
      sscanf(s+5, "%d", &step);
      return(step);
    }
    s++;
  }
  printf("ERROR, step not found in wig file, %s", string);
  exit(8);
}

int Count_columns(char *string){
  int m=0; char *ptr=string;
  while(*ptr!='\n'){
    if((*ptr!=' ')&&(*ptr!='\t')){
      m++;
      while((*ptr!=' ')&&(*ptr!='\t')&&(*ptr!='\n'))ptr++;
    }else{
      ptr++;
    }
  }
  return(m);
}
void Plot_profile(struct peak *peak, int p,
		  long **xprof, float **yprof,
		  long *nprof, int nchr, char *name_prof,
		  char *nameout)
{

  int NBIN=80;
  float BINSIZE=50;
  int LENGTH=(NBIN*BINSIZE)/2;
  int chr=-1, i, j, nj, N_peak=0; 
  double Prof[NBIN]; int norm[NBIN];
  for(i=0; i<NBIN; i++){Prof[i]=0; norm[i]=0;}
  long *xp; float *yp, w;
  double Prof_ave=0, Prof_dev=0, normtot=0;

  struct peak *peak1=peak;
  while(peak1!=NULL){
    if(peak1->chr !=chr){
      if(peak1->chr >= nchr)break;
      chr=peak1->chr; j=0;
      nj=nprof[chr]; xp=xprof[chr]; yp=yprof[chr];
    }
    long x0=0.5*(peak1->end+peak1->ini);
    long xmin=x0-LENGTH, xmax=x0+LENGTH;
    while((j<nj)&&(*xp<xmin)){  //ini
      Prof_ave+=*yp; Prof_dev+=(*yp)*(*yp); normtot++;
      j++; xp++; yp++;
    }
    while((j<nj)&&(*xp<=xmax)){ //end
      int bin=(*xp-xmin)/BINSIZE;
      if(bin < NBIN){
	Prof[bin]+=(*yp); norm[bin]++;
      }
      Prof_ave+=*yp; Prof_dev+=(*yp)*(*yp); normtot++;
      j++; xp++; yp++;
    }
    peak1=peak1->next;
  }

  FILE *file_out;
  if(p==0){
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
  }else{
    file_out=fopen(nameout, "a");
  }
  fprintf(file_out, "# %s\n", name_prof);
  Prof_ave/=normtot;
  Prof_dev=sqrt(Prof_dev/normtot-Prof_ave*Prof_ave);
  for(i=0; i<NBIN; i++){
    fprintf(file_out, "%.0f %.3g\n",
	    i*BINSIZE-LENGTH,(Prof[i]/norm[i]-Prof_ave)/Prof_dev);
  }
  fprintf(file_out, "&\n");
  fclose(file_out);
}

char *Extension(char *string){
  char *s=string, *s0=s; int i=0, l=sizeof(string);
  while(*s!='\0'){if(*s=='.')s0=s+1; s++;}
  return(s0);
}


*/
