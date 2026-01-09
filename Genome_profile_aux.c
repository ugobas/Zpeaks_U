#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include "Genome_profile_aux.h"
#include "random3.h"

#define ZSCORE 1 // Print Z score of properties (1) or raw properties (0)?
#define NCHAR 400
#define VERBOSE 0
#define DELOG 0  // Exponentiate logarithm
#define W1 0.0 //0.25


int NCMAX=200;

unsigned long randomgenerator(void);
int INIRAN;
float RANFACTOR;
extern int PEAKSIZE;
static int Read_file_gr(long **x, float **y, char *file, int *LOG);
static int Read_file_wig(long **x, float **y, char *file, int *LOG);
struct fragshort *Read_frag_wig(long *Nfrag, int *span,
				int *Nchr, char **chr_name, int *chr_num,
				char *file);
static int Get_strand(char *);

void sort(int n, float *x);
void mysort(float *x, int n);

int Count_columns(char *string);

int Sum_matches(struct frag *ori1, int N_ori1, struct frag *ori2, int DTOL);
struct frag *Smaller_ini(struct frag **f_r, struct frag *f1, struct frag *f2);
struct frag *Smaller_end(struct frag **f_r, struct frag *f1, struct frag *f2);
int Overlap(struct frag *ori1, struct frag *ori2, int DTOL);
int Overlap_fx(struct fragshort *frag, long *xx, long k);

char *Extension(char *string);
int Get_step_wig(char *string);
int Get_chr_wig(int *span, char *scr, char *string);
int Get_chr_number(char *scr, char **chr_name, int *chr_num,
		   int *nchr, char *file);
int Chromosome_fragments(long *Lmax, long *Nf, struct fragshort *frags,
			 long Nfrag, int Nchr);
void Intersect_peaks(struct frag *peak,
		     struct frag *peak1,
		     struct frag *peak2);

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

int Read_file(long **x, float **y, char *file)
{
  // Check if exists
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  fclose(file_in);

  int LOG=0; long i, n;
  char *ext=Extension(file);
  if(strncmp(ext, "wig", 3)==0){
    n=Read_file_wig(x, y, file, &LOG);
  }else if(strncmp(ext, "gr", 2)==0 ||
	   strncmp(ext, "dat", 3)==0){
    n=Read_file_gr(x, y, file, &LOG);
  }else{
    printf("Reading file %s\n", file);
    printf("ERROR, unknown file type .%s\n", ext);
    exit(8);
  }

  // Exponentiate logarithm
  if(DELOG && LOG){
    float *yy=*y; for(i=0; i<n; i++){*yy=exp(*yy); yy++;}
  }
  return(n);
}

int Read_file_3f(long **x, float **y, float **z, char *file){
  long n=0; char string[1000];
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  // Count lines and allocate
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue; n++;
  }
  fclose(file_in);
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));
  (*z)=malloc(n*sizeof(float));

  // Read
  n=0;
  file_in=fopen(file, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    sscanf(string, "%ld%f%f", &((*x)[n]), &((*y)[n]), &((*z)[n])); n++; 
  }
  fclose(file_in);
  return(n);
}

int Read_bed_file(long **x, float **y, long *nn, char *file, int Nchr)
{
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }

  if(strncmp(Extension(file), "bed", 3)!=0){
    printf("WARNING, file %s has not bed extention\n", file);
  }
  fclose(file_in);

  int Nfrag, i; 
  struct fragshort *frags=Read_peaks_nochr(&Nfrag, file);

  // Shortest fragment 
  int Lmin=frags[0].end-frags[0].ini+1, k, l, chr, imin;
  long *Lmax=malloc(Nchr*sizeof(long));
  for(k=0; k<Nchr; k++)Lmax[k]=0;
  struct fragshort *frag=frags;
  for(i=0; i<Nfrag; i++){
    l=frag->end-frag->ini+1;
    if(l<Lmin){Lmin=l; imin=i;}
    chr=frag->chr;
    if((chr<0)||(chr>=Nchr)){
      printf("ERROR file %s, number of chromosome %d not in [1,%d]\n",
	     file, frag->chr, Nchr); exit(8);
    }
    if(frag->end > Lmax[chr])Lmax[chr]=frag->end;
    frag++;
  }
  printf("Chromosome sizes: ");
  for(k=0; k<Nchr; k++)printf(" %ld", Lmax[k]);
  printf("\nfile=%s\n", file);
  printf("imin: %d (over %d) %d %ld %ld\n", imin, Nfrag,
	 frags[imin].chr, frags[imin].ini, frags[imin].end);

  // Representative points
  int step=100;
  int start=step/2, c;
  frag=frags; i=0;
  for(c=0; c<Nchr; c++){
    if(Lmax[c]==0){
      printf("ERROR, no data found for chromosome %d in file %s\n",
	     c+1, file); exit(8);
    }
    long n=Lmax[c]/step+1; nn[c]=n;
    long *xx=malloc(n*sizeof(long)); x[c]=xx;
    float *yy=malloc(n*sizeof(float)); y[c]=yy;
    long k, xk=start; int n1=Nfrag-1;
    for(k=0; k<n; k++){
      xx[k]=xk; chr=c+1;
      while((chr < frag->chr)&&(i<n1)){i++; frag++;}
      if(chr==frag->chr){
	if(xk < frag->ini){
	  yy[k]=0;
	}else{
	  if(xk <frag->end){
	    yy[k]=1;
	  }else{
	    while((xk >=frag->end)&&(i<n1)){i++; frag++;}
	  }
	}
      }else{
	yy[k]=0;
      }
      xk+=step;
    }
  }

  free(frags);
  return(0);
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

int Read_multi_chr(long **x, float **y, long *nn, int **strand,
		   int *Nchr, int Nchr_in, char **chr_name, int *chr_num,
		   char *file)
{
  int BIN_SIZE=100;
  FILE *file_in=fopen(file, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  fclose(file_in);


  int icol=3, start=0, step=0;
  long Nfrag=0; struct fragshort *frags=NULL;
  if(strncmp(Extension(file), "wig", 3)==0){
    frags=Read_frag_wig(&Nfrag, &step, Nchr, chr_name, chr_num, file);
  }else{
    if(strncmp(Extension(file), "bed", 3)!=0){
      printf("WARNING, file %s has no known extension, treating it as .bed\n",
	     file);
    }
    //frags=Read_peaks(&Nfrag, &icol, Nchr, chr_name, chr_num, file);
    struct frag *ff=Read_frags(&Nfrag, &icol, Nchr, chr_name, chr_num, file);
    frags=malloc(Nfrag*sizeof(struct fragshort));
    struct fragshort *fs=frags; struct frag *f=ff; 
    for(int i=0; i<Nfrag; i++){
      fs->value=f->value; fs->strand=f->strand;
      fs->ini=f->ini; fs->end=f->end; fs->chr=f->chr;
      fs++; f++;
    }
    free(ff);
  }
  if(*Nchr > Nchr_in){
    printf("ERROR, %d chromosomes expected from configuration file but "
	   "%d found\n", Nchr_in, *Nchr); exit(8);
  }

  // Delog
  if(DELOG){
    int LOG=0; long i;
    for(i=0; i<Nfrag; i++)if(frags[i].value<0){LOG=1; break;}
    if(LOG){
      for(i=0; i<Nfrag; i++)frags[i].value=exp(frags[i].value);
    }
  }

  // Set step
  if(step==0){ // No wig file is present
    step=BIN_SIZE; printf("Setting step to %d\n", step);
  }

  // Number of fragments per chromosome
  long Lmax[Nchr_in], Nf[Nchr_in];
  for(int c=0; c<Nchr_in; c++){Lmax[c]=0; Nf[c]=0;}
  int Lmin=Chromosome_fragments(Lmax, Nf, frags, Nfrag, *Nchr);

  long Ltot=0; int c;
  printf("Chromosome Mb windows: (bin size=%d)\n", step);
  for(c=0; c<*Nchr; c++){
    if(Lmax[c]){nn[c]=(Lmax[c]/step)+1;}else{nn[c]=0;}
    printf("%s %.3f %ld\n", chr_name[c], Lmax[c]*0.000001, nn[c]);
    Ltot+=nn[c];
  }
  printf("Total windows: %ld\n", Ltot);

  if(strand){
    if(icol<5){
      printf("WARNING, strand information not present in file %s\n",file);
      for(c=0; c<*Nchr; c++)strand[c]=NULL;
    }else{
      printf("Obtaining strand information from file %s\n", file);
    }
  }

  // Make bins
  for(c=0; c<*Nchr; c++){
    if(Lmax[c]==0){
      printf("WARNING, no data found for chromosome %s in file %s\n",
	     chr_name[c], file);
      nn[c]=0; x[c]=NULL; y[c]=NULL; if(strand)strand[c]=NULL;
    }else{
      long n=nn[c], k=1;
      long *xx=malloc(n*sizeof(long));
      x[c]=xx; *xx=step; // xx[k] end point of fragment k
      while((*xx)<Lmax[c]){
	*(xx+1)=*xx+step; xx++; k++;
      }
      nn[c]=k;
      y[c]=malloc(nn[c]*sizeof(float));
      if((strand)&&(icol>=5)){
	strand[c]=malloc(nn[c]*sizeof(int));
	int *st=strand[c]; for(k=0; k<nn[c]; k++){*st=0; st++;}
      }
    }
  }

  float *yy=NULL;
  long nm=0, np=0, nd=0, nover=0;
  long *xx=NULL, n=0, k=0, n1=0;
  int *stra=NULL, printerr=0;
  struct fragshort *frag=frags; 
  c=-1;
  for(long i=0; i<Nfrag; i++){
    if(frag->chr!=c){
      c=frag->chr;
      xx=x[c]; yy=y[c]; n=nn[c];
      if(strand)stra=strand[c];
      for(k=0; k<n; k++){yy[k]=0;}
      n1=n-1; k=0;
    }
    int m=0;
    while((*(xx)  <= frag->ini)&&(k<n1)){xx++; k++;}
    while((*(xx-1)>=frag->end)&&(k>0)){xx--; k--;}
    while((*(xx)>frag->ini)&&((k==0)||(*(xx-1)<=frag->end))){
      int w=Overlap_fx(frag, xx, k);
      if(w){   
	yy[k]+=w*frag->value; m=1;
	if(stra){
	  if(frag->strand<0){
	    nm++;
	    if(stra[k]>0){stra[k]=BIDIR; nd++;}
	    else{stra[k]=-1;}
	  }else{
	    np++;
	    if(stra[k]<0){stra[k]=BIDIR; nd++;}
	    else{stra[k]=1;}
	  }
	}
      }
      xx++; k++;
    }
    if(m)nover++;
    if((printerr==0)&&(m==0)){
      printerr=1;
      printf("WARNING, fragment %d %ld %ld not assigned ",
	     frag->chr, frag->ini, frag->end);
      printf("chr= %d x-1=%ld x=%ld x+1=%ld k=%ld\n",
	     c+1, *(xx-1), *(xx), *(xx+1), k);
    }
    frag++;
  }
  printf("%ld fragments over %ld assigned\n", nover, Nfrag);
  //int mid=step/2;
  // Normalize to average value per nucleotide
  for(c=0; c<*Nchr; c++){
    //xx=x[c]; for(k=0; k<nn[c]; k++){*xx-=mid; xx++;}
    yy=y[c]; for(k=0; k<nn[c]; k++){*yy/=step; yy++;}
  }

  if(strand && strand[0]){
    printf("%ld coord. in negative strand %ld in positive %ld bi-directional\n",
	   nm,np,nd);
    // Assign strand as closest fragment
    long n0=0, n1=0, k; int l, s1, s2;
    for(c=0; c<*Nchr; c++){
      long n=nn[c];
      int *stra=strand[c];
      for(k=0; k<n; k++){
	if(stra[k]!=0)continue;
	n0++;
	for(l=0; l<10; l++){
	  int kl1=k+l, kl2=k-l;
	  if(kl2<0 || kl1 >=n)continue;
	  s1=fabs(stra[kl1]);
	  s2=fabs(stra[kl2]);
	  if(s1==1){
	    if((s2==0)||(s2==BIDIR)){stra[k]=stra[kl1]/2; n1++;}
	    break;
	  }else if(s2==1){
	    if((s1==0)||(s1==BIDIR)){stra[k]=stra[kl2]/2; n1++;}
	    break;
	  }
	} 
      }
    }
    printf("strand assigned to %ld coordinates out of %ld\n", n1, n0);
  }
  free(frags);
  /*Ltot=0;
  for(c=0; c<*Nchr; c++){
    printf("%s %.3f %ld\n", chr_name[c], Lmax[c]*0.000001, nn[c]);
    Ltot+=nn[c];
  }
  printf("Total windows: %ld\n", Ltot);*/


  return(0);
}

int Overlap_fx(struct fragshort *frag, long *xx, long k)
{
  long xmax, xmin;
  if(*xx<=frag->end){xmax=*xx;}else{xmax=frag->end+1;}
  if((k==0)||(*(xx-1)<frag->ini)){xmin=frag->ini;}else {xmin=*(xx-1);}
  return(xmax-xmin);
}

int Read_file_bed(long **x, float **y, char *file, int *LOG,
		  int *Nchr, char **chr_name, int *chr_num)
{
  int icol; long Nfrag;
  struct fragshort *frags=
    Read_peaks(&Nfrag, &icol, Nchr, chr_name, chr_num, file);

  // Shortest fragment 
  int Lmin=frags[0].end-frags[0].ini, i, l;
  long Lmax=frags[0].end;
  for(i=1; i<Nfrag; i++){
    l=frags[i].end-frags[i].ini;
    if(l<Lmin)Lmin=l;
    if(frags[i].end > Lmax)Lmax=frags[i].end;
  }

  // Representative points
  int start=Lmin/4, step=Lmin/2;
  long k, xk=start, n=Lmax/step+1; i=0;
  struct fragshort *frag=frags;
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));
  for(k=0; k<n; k++){
    (*x)[k]=xk;
    if(xk < frag->ini){
      (*y)[k]=0;
    }else{
      if(xk <frags->end){
	 (*y)[k]=1;
      }else{
	while(xk >=frags->end){
	  i++; frag++;
	  if(i>=Nfrag){frag--; break;}
	}
      }
    }
    xk+=step;
  }


  free(frags);
  return(n);
}

int Read_file_wig(long **x, float **y, char *file, int *LOG)
{
  long n=0; char string[1000];
  FILE *file_in=fopen(file, "r");
  fgets(string, sizeof(string), file_in);
  fgets(string, sizeof(string), file_in);
  int step=Get_step_wig(string);

  // Count lines and allocate
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    n++;
  }
  fclose(file_in);
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));

  // Read
  n=0; *LOG=0;
  long *xx=*x, xxcount=0; *xx=step;
  float *yy=*y;
  file_in=fopen(file, "r"); // if(WIG)
  fgets(string, sizeof(string), file_in);
  fgets(string, sizeof(string), file_in);
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    // if(WIG){
    sscanf(string, "%f", yy);
    xxcount+=step; *xx=xxcount;
    if((DELOG)&&(*LOG==0)&&(*yy<0))*LOG=1;
    n++; xx++; yy++;
  }
  fclose(file_in);
  return(n);
}

int Read_file_gr(long **x, float **y, char *file, int *LOG)
{
  FILE *file_in=fopen(file, "r");
  long n=0; char string[1000];
  // Count lines and allocate
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue; n++;
  }
  fclose(file_in);
  printf("%ld lines in file %s\n", n, file);
  (*x)=malloc(n*sizeof(long));
  (*y)=malloc(n*sizeof(float));

  // Read
  n=0; *LOG=0;
  long *xx=*x;
  float *yy=*y;
  char word[100];
  file_in=fopen(file, "r"); 
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    sscanf(string, "%ld%s", xx, word);
    if(word[0]=='n'){*yy=0;}
    else{sscanf(word, "%f", yy);}
    if((DELOG)&&(*LOG==0)&&(*yy<0))*LOG=1;
    n++; xx++; yy++;
  }
  fclose(file_in);
  return(n);
}


int Read_file_3(long **kk, float **x1, float **x2, char *file)
{
  if(file==NULL)return(0);
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  // Count lines and allocate
  long n=0; char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    n++;
  }
  fclose(file_in);
  (*x1)=malloc(n*sizeof(float));
  (*x2)=malloc(n*sizeof(float));
  (*kk)=malloc(n*sizeof(long));

  // Read
  n=0;
  file_in=fopen(file, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    sscanf(string, "%ld%f%f", &((*kk)[n]), &((*x1)[n]), &((*x2)[n])); n++; 
  }
  fclose(file_in);
  return(n); 
}

int Read_window(char *file){
  int win=1; char string[150];
  FILE *file_in=fopen(file, "r");
  fgets(string, sizeof(string), file_in);
  if(strncmp(string, "#WINDOW=",8)==0){
    sscanf(string+8, "%d", &win);
  }else{
    printf("WARNING, window size not found in %s\n", file);
    printf("Setting default value %d\n", win);
  }
  fclose(file_in);
  return(win);
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

float Corr_coeff(float *slope, float *offset, float *xx, float *yy, int n)
{
  float *x=xx, *y=yy; int i;
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  for(i=0; i<n; i++){
    x1+=(*x); x2+=(*x)*(*x);
    y1+=(*y); y2+=(*y)*(*y);
    xy+=(*x)*(*y); x++; y++;
  }
  x1/=n; y1/=n;
  x2=(x2-n*x1*x1);
  y2=(y2-n*y1*y1);
  if(slope)*slope = xy/x2;
  if(offset)*offset=(y1-(*slope)*x1);
  xy=(xy-n*x1*y1)/sqrt(x2*y2);
  return(xy);
}


float Corr_coeff_old(float *xx, float *yy, int n)
{
  float *x=xx, *y=yy; int i;
  double x1=0, x2=0, y1=0, y2=0, xy=0;
  for(i=0; i<n; i++){
    x1+=(*x); x2+=(*x)*(*x);
    y1+=(*y); y2+=(*y)*(*y);
    xy+=(*x)*(*y); x++; y++;
  }
  x1/=n; y1/=n;
  xy=(xy-n*x1*y1)/sqrt((x2-n*x1*x1)*(y2-n*y1*y1));
  return(xy);
}

FILE *Open_file_w(char *name)
{
  FILE *file_out=fopen(name, "w");
  printf("\nWriting %s\n", name);
  return(file_out);
}


struct frag *Read_frags(long *Nfrag, int *icol,
			int *nchr, char **chr_name, int *chr_num,
			char *file)
{
  if(file==NULL)return(0);
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }

  // Count lines and allocate
  long nline=0; char string[1000]; *icol=-1;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    if(nline==0)*icol=Count_columns(string);
    nline++;
  }
  fclose(file_in);
  printf("%ld lines and %d columns found in %s\n", nline, *icol, file);
  if(nline==0)return(NULL);
  if((*icol<3)||(*icol>6)){
    printf("WARNING, this program needs a bed file ");
    printf("with 3 to 6 columns\n");
    if(*icol<3){printf("Exiting\n"); return(NULL);}
  }

  // Allocate
  struct frag *frags=NULL; int box=1; long n=nline;
  while(frags==NULL){
    frags=malloc(n*sizeof(struct frag));
    if(frags==NULL){
      printf("WARNING!!!!!! No memory available to read file %s, %ld lines\n",
	     file, n); n/=2; box*=2;
    }
  }
  printf("Allocating %ld lines for file %s, box=%d\n", n, file, box);
  struct frag *frag=frags;

  int nchr_old=*nchr, i;

  // Read
  int multi=0, c=0, ib, ichr=-1; 
  long k, n1, ns=0;
  char scr[80], chr_name_act[80]="       ";
  char strand[8], name[40], nameold[40], value[40]=".";
  file_in=fopen(file, "r"); nline=0; int nexcl=0;
  double Sum1=0, Sum2=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    if(*icol==3){
      sscanf(string,"%s%ld%ld",scr, &(frag->ini), &(frag->end));
    }else if(*icol==4){
      sscanf(string,"%s%ld%ld%s",scr,&(frag->ini),&(frag->end),value);
    }else if(*icol==5){
      sscanf(string, "%s%ld%ld%s%s", scr,
	     &(frag->ini), &(frag->end), value, strand);
    }else if(*icol>=6){
      sscanf(string, "%s%ld%ld%s%s%s", scr,
	     &(frag->ini), &(frag->end), name, value, strand);
    }
    frag->strand=Get_strand(strand);
    if(frag->strand<0)ns++;

    if(*icol>3){
      if(nline==0){
	strcpy(nameold, name); n1=nline; multi=1;
      }else if(strcmp(name, nameold)!=0){
	for(i=n1; i<nline; i++)frags[i].multi=multi;
	strcpy(nameold, name); n1=nline; multi=1; 
      }else{
	multi++;
      }
    }
    if(isdigit(name[0])){
      sscanf(name, "%f", &(frag->value));
      if(isnan(frag->value))(frag->value)=1;
    }else{
      frag->value=1;
    }
    if(isdigit(value[0])){
      sscanf(value, "%f", &(frag->value));
      if(isnan(frag->value))(frag->value)=1;
    }else{
      frag->value=1;
    }
    if(*icol!=3){Sum1+=frag->value; Sum2+=frag->value*frag->value;}

    if(strcmp(scr, chr_name_act)!=0 || ichr<0){
      ichr=Get_chr_number(scr, chr_name, chr_num, nchr, file);
      strcpy(chr_name_act, scr);
    }

    frag->chr=ichr;
    nline++;  frag++;
    //for(ib=1; ib<box; ib++)
    //if(fgets(string,sizeof(string),file_in)==NULL)goto close;
  }
 close:
  fclose(file_in);
  *Nfrag=nline;
  printf("%d chromosomes with %ld fragments found in %s:\n",
	 *nchr, nline, file);

  if(nchr_old && *nchr!=nchr_old){
    printf("WARNING, different number of chr: %d input, %d in %s\n",
	   nchr_old, *nchr, file);
  }


  //printf("%d columns\n", *icol);
  if(*icol>=5)printf("%ld fragments in negative strand %ld in positive\n",
		     ns, nline-ns);
  if(nexcl)printf("%d fragments omitted\n", nexcl);
  if(*icol!=3){
    Sum1/=*Nfrag;
    Sum2-=(*Nfrag)*Sum1*Sum1; Sum2/=(*Nfrag-1);
    printf("Mean: %.3g S.dev: %.3g\n", Sum1, sqrt(Sum2));
    if((fabs(Sum1)<0.00001)&&(Sum2<0.0001)){
      // value is not specified, setting it to 1
      printf("WARNING, all values are 0, setting them to 1\n");
      frag=frags;
      for(i=0; i<*Nfrag; i++){frag->value=1; frag++;}
    }
  }
  frag=frags;
  for(i=0; i<*Nfrag; i++){frag->next=frag+1; frag=frag->next;}
  frags[(*Nfrag)-1].next=NULL;
  /*frag=frags;
  for(i=0; i<*Nfrag; i++){
    printf("%d %d %d\n", frag->chr, frag->ini, frag->end); frag++;
    }*/
  return(frags); 
}


struct fragshort *Read_peaks_nochr(int *nori, char *file)
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
    if(string[0]=='#'||string[0]=='\n')continue;
    n++;
  }
  fclose(file_in);
  if(n==0)return(NULL); *nori=n;

  // Read
  struct fragshort *ori=malloc(n*sizeof(struct fragshort)), *orii=ori; 
  int s=0, mult=0, chr=0; char schr[5], schr_old[5]="     ";
  long i=0;
  file_in=fopen(file, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    char *ptr=string;
    if((string[0]=='c')||(string[0]=='C'))ptr+=3;
    s=sscanf(ptr, "%s%ld%ld%d", schr, &(orii->ini), &(orii->end), &mult);
    if(strcmp(schr,schr_old)!=0){strcpy(schr_old, schr); chr++;}
    (orii->chr)=chr;
    if(mult>0){orii->multi=mult;}else{orii->multi=1;}
    //Set_ori(orii, i);
    i++;  orii++;
  }
  //ori[n-1].next=NULL;
  fclose(file_in);
  printf("%d chromosomes found Peaksize= %d\n", chr, PEAKSIZE);
  return(ori); 
}


struct fragshort *Read_peaks(long *Nfrag, int *icol,
			     int *nchr, char **chr_name, int *chr_num,
			     char *file)
{
  if(file==NULL)return(0);
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }
  // Count lines and allocate
  long nline=0; char string[1000];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    if(nline==0)*icol=Count_columns(string);
    nline++;
  }
  fclose(file_in);
  printf("%ld lines and %d columns found in %s\n", nline, *icol, file);
  if(nline==0)return(NULL);
  if((*icol!=3)&&(*icol!=4)&&(*icol!=6)){
    printf("WARNING, this program only recognizes bed files ");
    printf("with 3, 4 or 6 columns\n");
    if(*icol<3){printf("Exiting\n"); return(NULL);}
  }

  // Allocate
  struct fragshort *frags=NULL; int box=1; long n=nline;
  while(frags==NULL){
    frags=malloc(n*sizeof(struct fragshort));
    if(frags==NULL){
      printf("WARNING!!!!!! No memory available to read file %s, %ld lines\n",
	     file, n); n/=2; box*=2;
    }
  }
  printf("Allocating %ld lines for file %s, box=%d\n", n, file, box);
  struct fragshort *frag=frags;

  int nchr_old=*nchr;

  // Read
  int multi=0, ib, ichr, mchr=0;
  long i, n1, ns=0;
  char crn[80], scr_old[80]="       ", scr_act[80]="";
  char strand[8], name[40], nameold[40], value[40]=".";
  file_in=fopen(file, "r"); nline=0;
  double Sum1=0, Sum2=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    frag->strand=0;
    if(*icol==3){
      sscanf(string,"%s%ld%ld",crn, &(frag->ini), &(frag->end));
    }else if(*icol==4){
      sscanf(string,"%s%ld%ld%f",crn,&(frag->ini),&(frag->end),value);
    }else if(*icol==5){
      sscanf(string,"%s%ld%ld%f%s",crn,&(frag->ini),&(frag->end),value,strand);
      frag->strand=Get_strand(strand);
    }else if(*icol==6){
      sscanf(string, "%s%ld%ld%s%s%s", crn,
	     &(frag->ini), &(frag->end), name, value, strand);
      frag->strand=Get_strand(strand);
      if(nline==0){
	strcpy(nameold, name); n1=nline; multi=1;
      }else if(strcmp(name, nameold)!=0){
	for(i=n1; i<nline; i++)frags[i].multi=multi;
	strcpy(nameold, name); n1=nline; multi=1; 
      }else{
	multi++;
      }
    }
    if(isdigit(value[0])){
      sscanf(value, "%f", &(frag->value));
      if(isnan(frag->value))frag->value=1;
    }else{
      frag->value=1;
    }

    char *scr=crn; if((*scr=='c')||(*scr=='C'))scr+=3;
    int ichr=mchr;
    if(strcmp(scr,scr_act)==0)goto found;
    for(ichr=0; ichr<*nchr; ichr++)
      if(strcmp(scr, chr_name[ichr])==0)goto found;

    // New chromosome
    if(strcmp(scr, scr_act)!=0){
      printf("New chromosome: %d %s\n", *nchr, scr);
      strcpy(scr_act, scr);
      mchr=*nchr; ichr=mchr; (*nchr)++;
      if(*nchr>=NCMAX){
	printf("ERROR, too many chromosomes (more than %d)\n", NCMAX);
	printf("If correct, modify NMCAX in Genome_profile_aux.c\n");
	exit(8);
      }
      if(chr_name[mchr]==NULL){chr_name[mchr]=malloc(80*sizeof(char));}
      strcpy(chr_name[mchr], scr);
    }

  found:
    if(ichr<0){
      printf("ERROR, chromosome number %d (%s) not valid\n",ichr, scr);
      exit(8);
    }

    frag->chr=ichr;
    if(*icol!=3){Sum1+=frag->value; Sum2+=frag->value*frag->value;}
    if(frag->strand<0)ns++;  
    nline++;  frag++;
    for(ib=1; ib<box; ib++)
      if(fgets(string,sizeof(string),file_in)==NULL)goto close;
  }
 close:
  fclose(file_in);
  *Nfrag=nline;
  printf("%d chromosomes with %ld fragments found in %s", *nchr, nline, file);
  if(*nchr != nchr_old){
    printf("WARNING, different number of chromosomes in input (%d) "
	   "and in %s (%d)\n", nchr_old, file, *nchr);
  }

  if(*icol==6)printf(" %ld of them in negative strand", ns);
  printf("\n");
  if(*icol!=3){
    Sum1/=*Nfrag;
    Sum2-=(*Nfrag)*Sum1*Sum1; Sum2/=(*Nfrag-1);
    printf("Mean: %.3g S.dev: %.3g\n", Sum1, sqrt(Sum2));
    if((fabs(Sum1)<0.00001)&&(Sum2<0.0001)){
      // value is not specified, setting it to 1
      printf("WARNING, all values are 0, setting them to 1\n");
      frag=frags;
      for(i=0; i<*Nfrag; i++){frag->value=1; frag++;}
    }
  }
  return(frags); 
}

struct fragshort *Read_frag_wig(long *Nfrag, int *span,
				int *nchr, char **chr_name, int *chr_num,
				char *file)
{
  if(file==NULL)return(0);
  FILE *file_in=fopen(file, "r");
  // Check if exists
  if(file_in==NULL){
    printf("ERROR, file %s does not exist\n", file); exit(8);
  }

  // Count lines and allocate
  long nline=0; char string[1000]; int icol;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#'||string[0]=='\n')continue;
    if(strncmp(string,"track",5)==0 || strncmp(string,"fixed",5)==0)continue;
    if(nline==0)icol=Count_columns(string);
    nline++;
  }
  fclose(file_in);
  printf("%ld lines %d columns found in %s\n", nline, icol, file);
  if(nline==0)return(NULL);

  // Allocate
  struct fragshort *frags=NULL; int box=1; long n=nline;
  while(frags==NULL){
    frags=malloc(n*sizeof(struct fragshort));
    if(frags==NULL){
      printf("WARNING!!!!!! No memory available to read file %s, %ld lines\n",
	     file, n); n/=2; box*=2;
    }
  }
  printf("Allocating %ld lines for file %s, box=%d\n", n, file, box);
  struct fragshort *frag=frags;

  // Read
  int variable=0, span1, c, shift=-1;
  long x=1;
  char scr[7];
  file_in=fopen(file, "r"); nline=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "track", 5)==0)continue;
    if(string[0]=='#'||string[0]=='\n')continue;
    if((strncmp(string, "variable", 8)==0)||
       (strncmp(string, "fixed", 5)==0)){
      if(strncmp(string, "fixed", 5)==0){variable=0;}
      else{variable=1;}
      Get_chr_wig(span, scr, string); span1=*span-1;
      c=Get_chr_number(scr, chr_name, chr_num, nchr, file);
      if(shift<0){if(c){shift=1;}else{shift=0;}}
      c-=shift;
      x=-1;
      continue;
    }else if(variable){
      sscanf(string, "%ld%f", &(frag->ini), &(frag->value));
      frag->end=frag->ini+span1; frag->strand=0; 
    }else{
      sscanf(string, "%f", &(frag->value));
      frag->ini=x;  frag->end=x+span1; frag->strand=0;
      x=frag->end+1;
    }
    (frag->chr)=c;
    //Set_ori(orii, nline);
    nline++;  frag++;
  }
  fclose(file_in);
  *Nfrag=nline;
  printf("%d chromosomes with %ld fragments found in %s\n", *nchr, nline, file);
  return(frags);
}


void Set_ori(struct frag *orii, int i){
  orii->size=orii->end-orii->ini;
  orii->x0=(orii->end+orii->ini)/2;
  orii->x1=orii->x0-PEAKSIZE;
  orii->x2=orii->x0+PEAKSIZE;
  orii->match=0;
  orii->lab=i;
  orii->next=orii+1;
  //orii->w=1;
}

void Count_matches(int *n_match1, int *n_match2, int DTOL,
		   struct frag *ori1, int N_ori1,
		   struct frag *ori2, int N_ori2)
{
  *n_match1=Sum_matches(ori1, N_ori1, ori2, DTOL);
  *n_match2=Sum_matches(ori2, N_ori2, ori1, DTOL);
}

int Sum_matches(struct frag *ori1, int N_ori1, struct frag *ori2, int DTOL)
{
  int n=0, m=0;
  struct frag *o=ori1, *p_l=ori2, *p_r=p_l->next;
  while(o){
    o->match=0;
    while(p_r && ((p_r->chr<o->chr)||((p_r->chr==o->chr)&&(p_r->end<o->ini)))){
      p_l=p_r; p_r=p_r->next;
    }
    if(Overlap(p_l, o, DTOL)||Overlap(p_r, o, DTOL)){o->match=1; m++;}
    o=o->next; n++;
  }
  printf("%d match over %d\n", m, n);
  return(m);
}


int Sum_match_d(int *n_matches, float *d_ave, int *same_strand,
		struct frag *ori1, int N_ori1,
		struct frag *ori2, int DTOL)
{
  // For each fragment in O1 count the number of frag in O2 that match them

  int n=0, m=0; // All frags in O1, frags in O1 with at least 1 match
  *n_matches=0; // Total matches
  *d_ave=0;     // Mean of min dist of frags in O1 to O2
  if(same_strand)*same_strand=0;
  struct frag *o=ori1, *p_l=ori2, *p_r=p_l->next;
  while(o){
    o->match=0;
    while(p_r && ((p_r->chr<o->chr)||((p_r->chr==o->chr)&&(p_r->end<o->ini)))){
      p_l=p_r; p_r=p_l->next;
    }
    float d_min;
    if(p_l && p_l->chr==o->chr){d_min=fabs(o->x0-p_l->x0);}
    else{d_min=1000000000;}
    while(p_r && ((p_r->chr==o->chr)&&(p_r->ini<o->end))){
      float d=fabs(o->x0-p_r->x0); if(d<d_min)d_min=d;
      if(Overlap(o, p_r, DTOL)){
	if(same_strand && (o->strand==p_r->strand))(*same_strand)++;
	o->match=1;
	(*n_matches)++;
      }
      if(p_r->end<=o->end){p_r=p_r->next;}
      else{break;}
    }
    if(o->match){m++;}
    (*d_ave)+=d_min; n++;
    p_l=p_r;
    o=o->next;
  }
  (*d_ave)/=n;
  //printf(" %d matches over %d\n", m, n);
  return(m);
}


int Overlap(struct frag *f1, struct frag *f2, int DTOL)
{
  if((f1==NULL)||(f2==NULL)||(f1->chr!=f2->chr))return(0);
  if((f1->end+DTOL<f2->ini)||(f2->end+DTOL<f1->ini))return(0);
  return(1);
}

int Unify_peaks(struct frag *ori_new,
		struct frag *ori1, int N_ori1,
		struct frag *ori2, int N_ori2,
		int DTOL)
{
  struct frag *ori=ori_new;
  struct frag *ori1_p=ori1, *ori2_p=ori2;
  int n=0, m=0;
  while(1){
    if((ori1_p==NULL)&&(ori2_p==NULL))break;
    // New ori is equal to smaller ori
    struct frag *ori_r, *ori_l=Smaller_ini(&ori_r, ori1_p, ori2_p);
    ori->chr=ori_l->chr; ori->ini=ori_l->ini;
    ori->end=ori_l->end; ori->next=ori+1;
    ori->match=0;
    // set label and update copied ori
    if(ori_l==ori1_p){ori->lab=1; ori1_p=ori1_p->next;}
    else{ori->lab=2; ori2_p=ori2_p->next;}

    while(Overlap(ori_l,ori_r,DTOL)){
      // overlap
      ori_l->match=1; ori_r->match=1; 
      ori_l=Smaller_end(&ori_r, ori_l, ori_r);
      ori_l=ori_l->next;
      ori->end=ori_r->end;
      ori->match++;
      ori->lab=3; m++;
    }
    if(ori->match){
      while(ori1_p &&(ori1_p->chr==ori->chr)&&(ori1_p->end<=ori->end))
	ori1_p=ori1_p->next;
      while(ori2_p &&(ori2_p->chr==ori->chr)&&(ori2_p->end<=ori->end))
	ori2_p=ori2_p->next;
    }
    ori++; n++;
  }
  (ori-1)->next=NULL;
  printf("%d overlapping fragments %d union\n", m, n);
  return(n);
}

struct frag *Smaller_ini(struct frag **f_r, struct frag *f1, struct frag *f2)
{
  if(f1==NULL){*f_r=f1; return(f2);}
  if(f2==NULL){*f_r=f2; return(f1);}
  if((f1->chr<f2->chr)||((f1->chr==f2->chr)&&(f1->ini<=f2->ini))){
    *f_r=f2; return(f1);
  }else{
    *f_r=f1; return(f2);
  }
}

struct frag *Smaller_end(struct frag **f_r, struct frag *f1, struct frag *f2)
{
  if(f1==NULL){*f_r=f1; return(f2);}
  if(f2==NULL){*f_r=f2; return(f1);}
  if((f1->chr<f2->chr)||((f1->chr==f2->chr)&&(f1->end<=f2->end))){
    *f_r=f2; return(f1);
  }else{
    *f_r=f1; return(f2);
  }
}
  
void Intersect_peaks(struct frag *peak,
		     struct frag *peak1,
		     struct frag *peak2)
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
  if(peak->match==0){peak->multi=peak1->multi+peak2->multi;}
}

float Dist_peaks(struct frag *ori1, int N_ori1,
		 struct frag *ori2, int N_ori2)
{
  double d_sum=0; float d1, d2; int n=0;
  struct frag *ori=ori1, *ori_l=ori2, *ori_r;
  while(ori!=NULL){
    while((ori_l!=NULL)&&(ori_l->chr<ori->chr))ori_l=ori_l->next;
    if(ori_l!=NULL)ori_r=ori_l->next;
    while((ori_l!=NULL)&&(ori_l->chr==ori->chr)&&
	  (ori_r!=NULL)&&(ori_r->chr==ori->chr)&&(ori_r->x0 < ori->x0)){
      ori_l=ori_r; ori_r=ori_l->next;
    }
    if((ori_l!=NULL)&&(ori_l->chr==ori->chr))
      {d1=fabs(ori_l->x0-ori->x0);}else{d1=10000000;}
    if((ori_r!=NULL)&&(ori_r->chr==ori->chr))
      {d2=fabs(ori_r->x0-ori->x0);}else{d2=10000000;}
    if(d1<d2){d_sum+=d1;}else{d_sum+=d2;} n++;
    ori=ori->next;
  }
  return(d_sum/n);
}


struct fragshort *Fragshort(struct frag *frags, int N){
  struct fragshort *fragnew=malloc(N*sizeof(struct fragshort));
  struct frag *f=frags; struct fragshort *x=fragnew; 
  int i;
  for(i=0; i<N; i++){
    x->chr=f->chr;
    x->ini=f->ini;
    x->end=f->end;
    x++; f++;
  }
  return(fragnew);
}

struct frag *Extract_peaks(struct frag *frag_old, int Nf, int Nchr)
{
  struct frag *frag=malloc(Nf*sizeof(struct frag));

  unsigned long iran;
  if(INIRAN==0){
    INIRAN=1;
    iran=randomgenerator();
    InitRandom( (RANDOMTYPE)iran);
    RANFACTOR=pow(12.0,1/3.0);
  }

  int i, chr; 
  int num_chr[Nchr]; long *size=malloc(Nchr*sizeof(long));
  for(i=0; i<Nchr; i++){size[i]=0; num_chr[i]=0;}
  struct frag *f_old=frag_old;
  for(i=0; i<Nf; i++){
    chr=f_old->chr;
    if((chr<0)||(chr>=Nchr))printf("frag %d chr %d\n", i, f_old->chr);
    num_chr[chr]++;
    if(f_old->end > size[chr])size[chr]=f_old->end;
    f_old++;
  }


  //printf("Extracting %d frags\n", Nf);
  struct frag *f=frag; int nn=0;
  for(chr=0; chr<Nchr; chr++){
    int n=num_chr[chr], k;
    if(n==0)continue;
    long L=size[chr];
    float ini_chr[n];
    for(i=0; i<n; i++)ini_chr[i]=RandomFloating()*L;
    //mysort(ini_chr, n);
    struct frag *f_old=frag_old;
    for(k=0; k<n; k++){
      int imin=-1;
      for(i=0; i<n; i++){
	if((ini_chr[i]>=0)&&((imin<0)||(ini_chr[i]<ini_chr[imin])))imin=i;
      }
      if(imin<0){
	printf("ERROR in random, imin not found: chr= %d i= %d over %d\n",
	       chr, k, n); exit(8);
      }
      f->ini=ini_chr[imin]; ini_chr[imin]=-1;
      //f->ini=ini_chr[k];
      while(f_old->chr!=chr)f_old++;
      if(f_old->chr!=chr){
	printf("ERROR in random, f_old not found: chr= %d i= %d over %d\n",
	       chr, k, n); exit(8);
      }
      f->chr=chr;
      f->size=f_old->size;
      f->end=f->ini+f->size;
      f->x0=(f->end+f->ini)/2;
      f->x1=f->x0-PEAKSIZE;
      f->x2=f->x0+PEAKSIZE;
      f->next=f+1;
      f++; nn++;
    }
  }
  (f-1)->next=NULL;
  if(nn!=Nf){
    printf("ERROR, %d instead of %d extracted fragments\n", nn, Nf);
    for(chr=0; chr<Nchr; chr++)
      printf("chr%d %d\n", chr+1, num_chr[chr]);
    for(i=0; i<Nf; i++){
      chr=(frag_old+i)->chr;
      if((chr<0)||(chr>=Nchr))
	printf("frag %d chr %d\n", i, (frag_old+i)->chr);
    }
    //exit(8);
  }


  return(frag);
}


void mysort(float *x, int n)
{
  struct rank{
    float x;
    struct rank *next;
  };
  struct rank rank_ptr[n+1], *rr=rank_ptr+1;
  rank_ptr->next=rr;
  rr->x=x[0]; rr->next=NULL; rr++;
  int i;
  for(i=1; i<n; i++){
    rr->x=x[i];
    struct rank *rp=rank_ptr, *rn=rp->next;
    while(rn != NULL){
      if(rr->x < rn->x){
	rp->next=rr;
	rr->next=rn;
	break;
      }
      rp=rn; rn=rn->next;
    }
    if(rn==NULL){rp->next=rr; rr->next=NULL;}
    rr++;
  }
  rr=rank_ptr->next;
  for(i=0; i<n; i++){
    x[i]=rr->x; rr=rr->next;
    
  }
  /*printf("Ranked values: \n");
    for(i=0; i<n; i++)printf("%.0f ", x[i]);
    printf("\n");*/

}


void sort(int n, float *ra)
{
  int l,j,ir,i;
  float rra;
  
  l=(n >> 1)+1;
  ir=n-1;
  for (;;) {
    if (l){
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[0];
      if (--ir == 0) {
	ra[0]=rra;
	return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
  }
}


unsigned long randomgenerator(void){
     
     unsigned long tm;
     time_t seconds;
     
     time(&seconds);
     srand((unsigned)(seconds % 65536));
     do   /* waiting time equal 1 second */
       tm= clock();
     while (tm/CLOCKS_PER_SEC < (unsigned)(1));
     
     return((unsigned long) (rand()));

}

int Join(struct frag *ori1, struct frag *ori2, int *joined)
{
  joined[ori1->lab]=1;
  joined[ori2->lab]=1;
  ori1->end=ori2->end;
  ori1->next=ori2->next;
  /*if(ori1->norm){
    ori1->norm+=ori2->norm;
    ori1->xsum+=ori2->xsum;
    ori1->x0=ori1->xsum/ori1->norm;
    if((ori1->x0>=ori1->end)||(ori1->x0<=ori1->ini)){
      printf("ERROR in join, wrong mid point %.0f (%.0f-%.0f)\n",
	     ori1->x0, ori1->ini, ori1->end); exit(8);
    }
    }else{*/
    ori1->x0=(ori1->end+ori1->ini)/2.;
    //}
  ori1->x1=ori1->x0-PEAKSIZE;
  ori1->x2=ori1->x0+PEAKSIZE;
  ori1->size=ori1->end-ori1->ini;
  return(1);
}

void Sizediff_stat(struct frag *ori, double *Prof_ori, double **Prof_score,
		   double *Discriminant_score,
		   int N_prof, char **name_prof, char *nameout,
		   int PRINTALL)
{
  //int MM = 3+2*N_prof;
  int M0=2; int MM = M0+N_prof;
  char *namex[MM]; int p, i, j, M;
  struct frag *ori1=NULL, *ori2=ori, *ori3=ori->next;
  struct frag *ori_next;
  double x_1[MM], x_2[MM], ave, dev, norm=0;
  float s, d, d1, d2;
  if(Prof_score){M=MM;}else{M=MM-N_prof;} //-2*N_prof;
  for(i=0; i<MM; i++){
    x_1[i]=0; x_2[i]=0; namex[i]=malloc(40*sizeof(char));
  }
  strcpy(namex[0], "size/10^3");
  strcpy(namex[1], "min_dist/10^4");
  for(i=0; i<N_prof; i++)strcpy(namex[M0+i],  name_prof[i]);

  // Print to file
  FILE *file_out;
  if(PRINTALL){file_out=Open_file_w(nameout);}
  else{file_out=fopen(nameout, "w");}
  fprintf(file_out, "### ");
  for(i=0; i<M; i++)fprintf(file_out," %s",namex[i]);
  fprintf(file_out, " chr x_peak\n");
  fprintf(file_out, "### ");
  for(i=0; i<2; i++)fprintf(file_out, " 0");
  for(i=2; i<M; i++)fprintf(file_out, " 1");
  fprintf(file_out, " 0 0\n");
  fprintf(file_out, "### "); for(i=0; i<M; i++)fprintf(file_out, " 1");
  fprintf(file_out, " 0 0\n");

  long n=0;
  while(ori2!=NULL){
    norm++;
    if((ori1!=NULL)&&(ori1->chr==ori2->chr)){d1=ori2->ini-ori1->end;}
    else{d1=1000000000;}
    ori3=ori2->next;
    if((ori3!=NULL)&&(ori3->chr==ori2->chr)){d2=ori3->ini-ori2->end;}
    else{d2=1000000000;}
    if(d1<d2){ori_next=ori1; d=d1;}else{ori_next=ori3; d=d2;}
    /*if(ori2->size>ori_next->size){omax=ori2; omin=ori_next;}
    else{omin=ori2; omax=ori_next;}
    xx[2]=omin->size/omax->size;*/

    if(ori2->size>SIZEMAX)fprintf(file_out, "#");
    i=0; s=ori2->size;
    fprintf(file_out, "%.4g ", s/1000); x_1[i]+=s; x_2[i]+=s*s;
    i=1; s=d;
    fprintf(file_out, "%.4g ",s/10000); x_1[i]+=s; x_2[i]+=s*s;
    for(j=0; j<N_prof; j++){
      s=Prof_score[j][n]; i++;
      fprintf(file_out, "%.4g ", s); x_1[i]+=s; x_2[i]+=s*s;
    }
    fprintf(file_out, " %d %ld\n", ori2->chr, ori2->x0);
    ori1=ori2; ori2=ori1->next; n++;
  }

  // Print to file
  fprintf(file_out, "# Profile sqrt(N_ori)(<S>_ori-<S>_noori)/s_noori");
  fprintf(file_out," (<S>_ori-<S>_noori)/s_ori\n");
  for(p=0; p<N_prof; p++){
    fprintf(file_out, "# %15s: %6.3f  %.3g\n",
	    name_prof[p], Prof_ori[p],  x_1[p+M0]/norm); 
  }   //Discriminant_score[p]
  fclose(file_out);

  // Print to screen
  if((PRINTALL)&&(VERBOSE)){
    printf("num ");
    for(i=0; i<M; i++)printf(" %s", namex[i]);
    printf("\n");
    printf("%5.0f  ", norm);
    for(i=0; i<M; i++){
      Print_ave(&ave, &dev, x_1[i], x_2[i], norm, namex[i], file_out);
      printf("  %.3g %.2g", ave, dev);
    }
    printf("\n");
    for(i=0; i<M; i++)printf(" %s", namex[i]); printf("\n");
  }
  /*printf("PEAKS  ");
    for(i=0; i<M; i++)printf(" %s", namex[i]); printf("\n");*/
  printf("%5.0f  ", norm);
  for(i=0; i<M; i++)printf(" %.3g", x_1[i]/norm); // i<M0
  //for(i=0; i<N_prof; i++)printf("%.3f ", Discriminant_score[i]);
  //printf("\n");
}

void Print_ave(double *ave, double *dev,
	       double sum1, double sum2, double norm, char *name,
	       FILE *file_out)
{
  *ave=sum1/norm;
  *dev=sum2/norm-(*ave)*(*ave); *dev=sqrt(*dev/norm);
  char string[1000];
  sprintf(string, "# %10s: ave=%.3g %.3g\n", name, *ave, *dev);
  fprintf(file_out, "%s", string);
  //printf("%s", string);
}

void Make_cluster_4(struct frag *ori, int N_ori, int *joined, float DCLUST)
{
  int i;
  for(i=0; i<N_ori; i++)joined[i]=0;
  int njoin=1, inicr=1;
  float d1, d2, d3;
  while(njoin){
    njoin=0;
    struct frag *ori1=NULL, *ori2=NULL, *ori3=NULL, *ori4=ori;
    while(ori4!=NULL){
      if(ori1==NULL){
	inicr=1; ori1=ori;
      }else if(ori4->chr!=ori3->chr){
	inicr=1; ori1=ori4;
	// Join at end of chromosome
	if((d2<DCLUST)&&(d2<d1)){
	  njoin+=Join(ori2, ori3, joined);
	}
      }
      if(inicr){
	inicr=0;
	ori2=ori1->next; d1=ori2->x0-ori1->x0;
	ori3=ori2->next; d2=ori3->x0-ori2->x0;
	if((d1<DCLUST)&&(d1<d2)){
	  njoin+=Join(ori2, ori3, joined);
	  ori2=ori1->next; d1=ori2->x0-ori1->x0;
	  ori3=ori2->next; d2=ori3->x0-ori2->x0;
	}
	ori4=ori3->next; d3=ori4->x0-ori3->x0;
      }
      if((d2<DCLUST)&&(d2<d1)&&(d2<d3)){
	// join 2 and 3
	njoin+=Join(ori2, ori3, joined);
	ori3=ori2->next; d2=ori3->x0-ori2->x0;
	ori4=ori3->next; if(ori4==NULL)break;
	d3=ori4->x0-ori3->x0;
	
      }
      ori1=ori2; ori2=ori3; ori3=ori4; ori4=ori3->next;
      if(ori4==NULL)break;
      d1=d2; d2=d3; d3=ori4->x0-ori3->x0;
    }
  }
}

int Make_cluster_1(struct frag *ori, int N_ori, int *joined, float DCLUST)
{
  int i, n=2;   float d1;
  for(i=0; i<N_ori; i++)joined[i]=0;
  struct frag *ori1=ori, *ori2=ori1->next;
  while(ori2!=NULL){
    if(ori2->chr!=ori1->chr){
      ori1=ori2; ori2=ori2->next;
      if(ori2==NULL)break;
    }
    d1=ori2->ini-ori1->end;
    if(d1<DCLUST){
      // join 2 and 1
      Join(ori1, ori2, joined);
    }else{
      ori1=ori2; n++;
    }
    ori2=ori2->next;
  }
  return(n);
}


void Profile_peaks(double **Prof_peaks,
		   double *Ave_all, double *Dev_all, long *N_bin,
		   double *Ave_peaks, double *Sem_peaks, 
		   double *Ave_nopeak, double *Sem_nopeak,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_p, int N_prof,
		   struct fragshort *peaks, int N_peaks)
{
  int p, i;
  float Ps[N_peaks]; int norm[N_peaks];
  for(p=0; p<N_prof; p++){
    if(Prof_peaks[p]){
      for(i=0; i<N_peaks; i++){Ps[i]=0; norm[i]=0;}
    }
    long Norm_peak=0; int Num_nopeak=0;
    double Prof_1_peak=0, Prof_2_peak=0;

    struct fragshort *peak=peaks;
    int chr=-1;
    long *xp, nj, j=0, k=0;
    float *yp;
    while(k<N_peaks){
      if(peak->chr!=chr){
	chr=peak->chr; j=0;
	nj=nprof[p][chr];
	xp=xprof[p][chr];
	yp=yprof[p][chr];
	Num_nopeak++; 
      }
      while(j<nj){
	if(j && *(xp-1)>peak->end){
	  break;
	}else if(*xp>peak->ini){
	  if(Prof_peaks[p]){Ps[k]+=*yp; norm[k]++;}
	  Prof_1_peak +=(*yp); //wpy;
	  Prof_2_peak+=(*yp)*(*yp); //wpy
	  Norm_peak++; //=wp;
	}
	xp++; yp++; j++;
      } // end bins of peak k
      k++; peak++;
    } // end peaks

    // Mean score of no-peaks
    Num_nopeak+=N_peaks; 
    double Prof_1_nopeak=N_bin[p]*Ave_all[p]-Prof_1_peak;
    double Prof_2_nopeak=N_bin[p]*(Dev_all[p]+Ave_all[p]*Ave_all[p])
      -Prof_2_peak;
    long Norm_nopeak=N_bin[p]-Norm_peak;
    Ave_nopeak[p]=Prof_1_nopeak/Norm_nopeak;
    Sem_nopeak[p]=Prof_2_nopeak/Norm_nopeak-Ave_nopeak[p]*Ave_nopeak[p];
    if(Sem_nopeak[p] > 0.00001){
      Sem_nopeak[p]=sqrt(Sem_nopeak[p]/Num_nopeak);
    }else{
      Sem_nopeak[p]=0;
    }

    // Mean score of peaks
    Ave_peaks[p]=Prof_1_peak/Norm_peak;
    Sem_peaks[p]=sqrt((Prof_2_peak/Norm_peak-Ave_peaks[p]*Ave_peaks[p])/
		       N_peaks);
    if(Prof_peaks[p]){
      for(i=0; i<N_peaks; i++){
	if(norm[i])Ps[i]/=norm[i];
	if(ZSCORE){Ps[i]=(Ps[i]-Ave_all[p])/Dev_all[p];}
	Prof_peaks[p][i]=Ps[i];
      }
    }
    //printf("%d Ave_peak= %.5f Ave_all= %.5f N= %ld %ld\n",
    //	   p, Ave_peaks[p], Ave_all[p], Norm_peak, N_bin[p]);

  }
  //exit(8);
}


void Profile_score(double *Prof_peak, double **Prof_score,
		   double *Ave_prof, double *Dev_prof,
		   double *Discriminant_score,
		   long ***xprof, float ***yprof, long **nprof,
		   int *nchr_p, int N_prof,
		   struct fragshort *peak, int N_PEAKS)
{
  int p, i;
  int *norm_peak=malloc(N_PEAKS*sizeof(int));
  float *Ps=malloc(N_PEAKS*sizeof(float));
  for(p=0; p<N_prof; p++){
    long Norm_peak=0, Norm_nopeak=0;
    double Prof_1_peak=0, Prof_2_peak=0;
    double Prof_1_nopeak=0, Prof_2_nopeak=0;
    int chr, chr_old=1;
    long j=0, k=0, nj=nprof[p][0], ww=nj, npk;
    long *xp=xprof[p][0]; float *yp=yprof[p][0];
    if(Prof_score){
      for(i=0; i<N_PEAKS; i++)norm_peak[i]=0;
      for(i=0; i<N_PEAKS; i++)Ps[i]=0;
    }
    struct fragshort *peak1=peak;
    while(k<N_PEAKS){
      while(j<nj){
	if(*xp<peak1->ini){
	  Prof_1_nopeak +=(*yp);
	  Prof_2_nopeak+=(*yp)*(*yp);
	  Norm_nopeak ++; j++; xp++; yp++;
	}else if(*xp<peak1->end){
	  if(Prof_score){
	    Ps[k]+=(*yp); norm_peak[k]++;
	  }
	  Prof_1_peak+=(*yp);
	  Prof_2_peak+=(*yp)*(*yp);
	  Norm_peak++; j++; xp++; yp++;
	}else{
	  break;
	}
      }
      //peak1=peak1->next; k++;
      peak1++; k++;
      if((peak1==NULL)||(peak1->chr != chr)){
	while(j<nj){
	  Prof_1_nopeak +=(*yp);
	  Prof_2_nopeak+=(*yp)*(*yp);
	  Norm_nopeak ++; j++; xp++; yp++;
	}
	if((k==N_PEAKS)||(peak1==NULL))break;
	chr=peak1->chr;
	if(chr<0 || chr>=nchr_p[p])break;
	nj=nprof[p][chr]; if(nj<=0)break;
	xp=xprof[p][chr]; yp=yprof[p][chr];
	if(xp==NULL || yp==NULL)break;
	ww+=nj; j=0;

      }
    }
    npk=k;

    // Mean score of all genome
    double ave_all=(Prof_1_peak+Prof_1_nopeak)/(Norm_peak+Norm_nopeak);
    double dev_all=(Prof_2_peak+Prof_2_nopeak)/(Norm_peak+Norm_nopeak);
    dev_all=sqrt(dev_all-ave_all*ave_all);

    // Mean score of no-peaks
    if(Norm_nopeak){
      Prof_1_nopeak/=Norm_nopeak;
      Prof_2_nopeak=sqrt(Prof_2_nopeak/Norm_nopeak-Prof_1_nopeak*Prof_1_nopeak);
    }
    // Mean score of peaks
    Prof_1_peak/=Norm_peak;
    Prof_2_peak=sqrt(Prof_2_peak/Norm_peak-Prof_1_peak*Prof_1_peak);
 
   // Score of individual peaks
    double ave_peak=0, dev_peak=0;
    if(Prof_score){
      double *Pp=Prof_score[p], w, w1=W1, w2=w1*w1;
      float z=(float)npk/(float)Norm_peak;
      for(k=0; k<npk; k++){
	if(W1){
	  double norm=norm_peak[k], sum=Ps[k];
	  j=k-2; if(j>=0) {w=w2; sum+=w*Ps[j]; norm+=w*norm_peak[j];}
	  j=k-1; if(j>=0) {w=w1; sum+=w*Ps[j]; norm+=w*norm_peak[j];}
	  j=k+1; if(j<npk){w=w1; sum+=w*Ps[j]; norm+=w*norm_peak[j];}
	  j=k+2; if(j<npk){w=w2; sum+=w*Ps[j]; norm+=w*norm_peak[j];}
	  if(norm){sum/=norm;}else{sum=ave_all;}
	  Pp[k]=sum;
	}else{
	  Pp[k]=Ps[k]; if(norm_peak[k])Pp[k]/=norm_peak[k];
	}
	ave_peak+=Pp[k]; dev_peak+=Pp[k]*Pp[k];
      }
      ave_peak/=npk;
      dev_peak=sqrt((dev_peak-npk*ave_peak*ave_peak)/(npk-1));
      if(Ave_prof){
	Ave_prof[p]=ave_peak;
	Dev_prof[p]=dev_peak;
      }
      if(ZSCORE){
	for(k=0; k<npk; k++)Pp[k]=(Pp[k]-ave_all)/dev_all;
      }
    }
    Prof_peak[p]=(ave_peak-ave_all)*sqrt((float)npk)/dev_peak;
    // Normalize by standard dev. of peaks. It has drawbacks.
    Discriminant_score[p]=(ave_peak-Prof_1_nopeak)/
      sqrt(dev_peak*Prof_2_nopeak);
    //sqrt((float)N_ori)
    // Discriminant_score[p]=(Prof_ori_1-Prof_1)/Prof_2;
    if(VERBOSE)
      printf("Profile %2d pk: %.3f %7ld nopk: %.3f %7ld "
	     "all: %.3f %7ld np: %ld\n",
	     p, ave_peak, Norm_peak, Prof_1_nopeak, Norm_nopeak, ave_all,
	     Norm_peak+Norm_nopeak, npk);
  }
  free(norm_peak); free(Ps);
}


int Remove_fragments(struct frag *ori, int N_ori, int *removed,
		     float Size_min)
{
  printf("Removing fragments with size < %.0f\n", Size_min);
  int n=N_ori;
  struct frag *ori1=ori, *ori2=ori1->next;
  while(ori2!=NULL){
    removed[ori2->lab]=0;
    if(ori2->size<Size_min){ori1->next=ori2->next; removed[ori2->lab]=1; n--;}
    else{ori1=ori2;}
    ori2=ori2->next;
  }
  printf("Obtaining %d from %d origins\n", n, N_ori);
  return(n);
}

char *Extension(char *string){
  char *s=string, *s0=s;
  while(*s!='\0'){if(*s=='.')s0=s+1; s++;}
  return(s0);
}

int Get_step_wig(char *string){
  char *s=string; int step;
  while(*s!='\n'){
    if(strncmp(s, "step=", 5)==0){
      sscanf(s+5, "%d", &step);
      return(step);
    }
    s++;
  }
  printf("ERROR, step not found in wig file, %s", string); exit(8);
}

void Plot_profile(struct frag *ori, int p,
		  long **xprof, float **yprof,
		  long *nprof, char *name_prof,
		  char *nameout)
{
  int cr, i, j, nj; long *xp; float *yp;
  struct frag *ori1=ori;
  int NBIN=80;
  float BINSIZE=50;
  int LENGTH=(NBIN*BINSIZE)/2;
  double *Prof=malloc(NBIN*sizeof(double));
  int *norm=malloc(NBIN*sizeof(int));
  for(i=0; i<NBIN; i++){Prof[i]=0; norm[i]=0;}
  double Prof_ave=0, Prof_dev=0, normtot=0;

  while(ori1!=NULL){
    if(ori1->chr !=cr){
      cr=ori1->chr; j=0;
      nj=nprof[cr]; xp=xprof[cr]; yp=yprof[cr];
    }
    long xmin=ori1->x0-LENGTH, xmax=ori1->x0+LENGTH;
    while((j<nj)&&(*xp<xmin)){  //ini
      Prof_ave+=*yp; Prof_dev+=(*yp)*(*yp); normtot++;
      j++; xp++; yp++;
    }
    while((j<nj)&&(*xp<=xmax)){ //end
      int bin=(*xp-ori1->x0+LENGTH)/BINSIZE;
      if(bin < NBIN){
	Prof[bin]+=(*yp); norm[bin]++;
      }
      Prof_ave+=*yp; Prof_dev+=(*yp)*(*yp); normtot++;
      j++; xp++; yp++;
    }
    ori1=ori1->next;
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
  fclose(file_out); free(Prof); free(norm);

}

int Get_chr_number(char *scr, char **chr_name, int *chr_num, int *nchr,
		   char *file)
{
  int c; char *s=scr;
  if(strncmp(scr, "chr", 3)==0 || strncmp(scr, "Chr", 3)==0 ){
    s=scr+3;
  }
  int ic=0; char *cn=chr_name[0];
  if(cn && (strncmp(cn, "chr", 3)==0 || strncmp(cn, "Chr", 3)==0)){ic=3;}
  for(c=0; c<*nchr; c++){
    if(strcmp(s, chr_name[c]+ic)==0)return(c);
  }

  for(c=0; c<*nchr; c++){
    int chr=chr_num[c];
    if(chr >=0 && chr<*nchr && strcmp(s, chr_name[chr]+ic)==0)return(chr);
  }

  // Not found in list, store it
  c=*nchr;
  /*int num=1; char *cc=s;
  while(*cc!='\0'){if(isdigit(*cc)==0){num=0; break;} cc++;}
  if(num){sscanf(s, "%d", &c); c--;}*/
  printf("Found new chromosome %d: %s\n", *nchr, scr);
  if(c<0){
    printf("ERROR, wrong chromosome number %d (%d)\n", s, c+1); exit(8);
  }
  if(*nchr>=NCMAX){
    printf("ERROR, too many chromosomes (more than %d)\n", NCMAX);
    printf("If correct, modify NMCAX in Genome_profile_aux.c\n");
    exit(8);
  }
  if(chr_num[*nchr]>=0){
    printf("ERROR, chromosome number %d already assigned: %d\n",
	   *nchr, chr_num[*nchr]); exit(8);
  }
  if(chr_name[*nchr]==NULL){chr_name[*nchr]=malloc(80*sizeof(char));}
  strcpy(chr_name[*nchr], s);
  chr_num[*nchr]=*nchr;
  (*nchr)++;
  return(c);
}

int Get_chr_wig(int *span, char *scr, char *string){
  char *s=string; int c=0, sp=0; 
  while(*s!='\n'){
    if(strncmp(s, "span=", 5)==0){
      sscanf(s+5, "%d", span); sp=1;
    }else if(strncmp(s, "chrom=", 6)==0){
      sscanf(s+6, "%s", scr); c=1;
    }
    s++;
  }
  if(c==0){
    printf("ERROR, chromosome name not found, line=%s\n",string);
    exit(8);
  }
  if(sp==0){
    printf("ERROR, span not found, line=%s\n",string);
    exit(8);
  }
  return(0);
}

void Print_peak(FILE *file_out, struct frag ori){
  /*if(PRINT==0){
    fprintf(file_out, "%d\t%.0f\t%.0f\n", ori.chr, ori.ini, ori.end);
    }else{*/
  if(ori.size>SIZEMAX)fprintf(file_out, "#");
  fprintf(file_out, "%d\t%d\t%d\n",
	  ori.chr, ori.ini, ori.end);
}

int Print_metaplot(double *Prof, char *name_prof, int p,
		   int NBIN, float BINSIZE, char *nameout)
{
  int LENGTH=((NBIN-1)*BINSIZE)/2, i;
  FILE *file_out;
  if(p==0){
    file_out=fopen(nameout, "w");
    printf("Writing %s\n", nameout);
  }else{
    file_out=fopen(nameout, "a");
  }
  fprintf(file_out, "# %s\n", name_prof);
  for(i=0; i<NBIN; i++){
    fprintf(file_out, "%.0f %.3g\n", i*BINSIZE-LENGTH, Prof[i]);
  }
  fprintf(file_out, "&\n");
  fclose(file_out);
  return(0);
}

struct frag *Genome_regions(long *N_regions, long *Chrsize,
			    int PEAKSIZE, int Nchr)
{
  struct frag *regions=NULL;
  long n=*N_regions; int box=1;
  while(regions==NULL){
    regions=malloc(n*sizeof(struct frag));
    if(regions==NULL){
      printf("WARNING!!! No memory available for genomic regions, %ld lines\n",
	     n);
      n/=2; box*=2;
    }
  }
  printf("Allocating %ld lines for genomic regions, box=%d\n", n, box);
  int size=PEAKSIZE*box;

  struct frag *region=regions, *next=regions+1;
  long i, x1=0, x0=size, x2=2*size;
  int chr=0, Frag_size=2*size;
  for(i=0; i<n; i++){
    if(x0 > Chrsize[chr]){chr++; x1=0, x0=PEAKSIZE, x2=2*PEAKSIZE;}
    if(chr>=Nchr)break;
    region->ini=x1; region->end=x2; region->size=Frag_size;
    region->x0=x0; region->x1=x1; region->x2=x2;
    region->lab=i; region->next=next;
    region->chr=chr;
    x1=x2; x0=x1+size; x2=x0+size;
    region=next; next++;
  }
  regions[i-1].next=NULL;
  *N_regions=i;
  return(regions);
}

char *Remove_dir(char *s){
  char *s1=s, *s2=s;
  while((*s1!='\0')&&(*s1!=' ')&&(*s1!='\n')){
    if(*s1=='/')s2=s1+1; s1++;
  }
  return(s2);
}

int Remove_file_extension(char *filename){
  char *s=filename;
  while(*s!='\0')s++;
  while(s!=filename){s--; if(*s=='.'){*s='\0'; break;}}
  return(0);
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

int Chromosome_fragments(long *Lmax, long *Nf, struct fragshort *frags,
			 long Nfrag, int Nchr)
{
  int chr, c; 
  long i, l, imin=0;
  long Lmin=frags[0].end-frags[0].ini+1;
  struct fragshort *frag=frags;
  int shift=0; if(frag->chr){shift=1;}
  for(i=0; i<Nfrag; i++){
    l=frag->end-frag->ini+1;
    if(l<Lmin){Lmin=l; imin=i;}
    chr=frag->chr; if(shift){chr--;} // CHANGE
    if((chr<0)||(chr>Nchr)){
      printf("ERROR frag %ld of %ld ", i, Nfrag);
      printf("(%d %ld %ld %f)", frag->chr,frag->ini, frag->end, frag->value);
      printf(", number of chromosome %d not in [1,%d]\n", frag->chr, Nchr);
      frag--;
      printf("previous: %d %ld %ld %f\n",
	     frag->chr,frag->ini,frag->end,frag->value);
      if(i<Nfrag-1){
	frag=frags+i+1;
	printf("next: %d %ld %ld %f\n",
	       frag->chr,frag->ini,frag->end,frag->value);
      }
      exit(8);
    } // end error
    Nf[chr]++;
    if(frag->end > Lmax[chr])Lmax[chr]=frag->end;
    frag++;
  }
  printf("Lmin= %ld\n", Lmin);
  /*if(Lmin<5)printf(" imin: line %ld (over %ld) chr %d %ld %ld", imin,
		   Nfrag,frags[imin].chr,frags[imin].ini,frags[imin].end);
		   printf("\n");*/
  return(Lmin);
}

int Orient_peaks(struct fragshort *peak, int N,
		 long **xpr, int **strand, long *nn, int Nchr)
{
  struct fragshort *peak1=peak;
  long j=0, k=0, nj=nn[0], *xp=xpr[0];
  int *stra=strand[0], chr=-1, ns=0, nm=0;
  int bidir=0;
  while(k<N){
    while(j<nj){
      if(*xp<peak1->ini){
	j++; xp++; if(stra)stra++;
      }else if(*xp<peak1->end){
	if(stra){
	  if(*stra==BIDIR){
	    peak1->strand=BIDIR;
	    bidir++;
	  }else if(*stra>0){
	    if(peak1->strand==0)peak1->strand=1;
	  }else if(*stra<0){
	    if(peak1->strand==0)peak1->strand=-1;
	    nm++;
	  }
	  stra++; ns++;
	}
	j++; xp++; 
      }else{
	break;
      }
    }
    peak1++; k++;
    if(k==N)break;
    if(peak1->chr != chr){
      chr=peak1->chr;
      j=0; nj=nn[chr];
      xp=xpr[chr]; stra=strand[chr];
    }
  }
  printf("strand assigned to %d fragments of %ld, %d negative %d positive\n",
	 ns, k, nm, ns-nm-bidir);
  return(bidir);
}

long Genome_statistics(double *Ave, double *Dev, float **y, long *nn, int Nchr)
{
  int Ntot=0; // Nvalue=0;
  (*Ave)=0; (*Dev)=0;
  for(int c=0; c<Nchr; c++){
    float *yy=y[c];
    for(int i=0; i<nn[c]; i++){ //nn[c]
      if(isnan(*yy))exit(8);
      (*Ave)+=(*yy); (*Dev)+=(*yy)*(*yy); //Nvalue++;
      yy++; //Ntot++;
    }
    Ntot+=nn[c];
  }
  (*Ave)=(*Ave)/Ntot;
  (*Dev)=sqrt(((*Dev)-Ntot*(*Ave)*(*Ave))/(Ntot-1));
  return(Ntot);
}

double Mean_size(struct frag *peaks, int N_peak){
  double size_tot=0; int i;
  struct frag *f=peaks;
  for(i=0; i<N_peak; i++){
    f->size=f->end-f->ini+1; size_tot+=f->size; f++;
  }
  return(size_tot/N_peak);
}

void Rank_frags_new(struct frag *frags, int n, int Nchr)
{
  int i, k, chr, m=0;
  // First fragment
  struct frag *f=frags, *fmin=f;
  for(i=0; i<n; i++){
    if((f->chr<=fmin->chr)&&(f->ini<fmin->ini))fmin=f;
    f->next=NULL; f++;
  }
  struct frag tmp=frags[0]; frags[0]=*fmin; *fmin=tmp;
  frags[0].next=frags+1;
  struct frag *frag_prev=frags;
  for(chr=1; chr<=Nchr; chr++){
    for(k=0; k<n; k++){
      struct frag *f=frags+1, *fmin=NULL;
      for(i=0; i<n; i++){
	if((f->chr==chr)&&(f->next==NULL)&&((fmin==NULL)||(f->ini<fmin->ini)))
	  fmin=f;
	f++;
      }
      if(fmin==NULL)break;
      frag_prev->next=fmin; frag_prev=fmin; fmin->next=fmin+1; m++;
    }
  }
  if(m!=n){
    printf("ERROR, %d frags but only %d ranked\n", n, m);
    exit(8);
  }
  frag_prev->next=NULL;
}

void Rank_frags(struct frag *frags, int n, int Nchr)
{
  int i, k, chr, m=0;
  struct frag tmp[n], *ftmp=tmp;
  for(chr=1; chr<=Nchr; chr++){
    for(k=0; k<n; k++){
      struct frag *f=frags, *fmin=NULL;
      for(i=0; i<n; i++){
	if((f->chr==chr)&&(f->ini >=0)&&((fmin==NULL)||(f->ini<fmin->ini)))
	  fmin=f;
	f++;
      }
      if(fmin==NULL)break;
      *ftmp=*fmin; fmin->ini=-1; ftmp++; m++;
    }
  }
  if(m!=n){
    printf("ERROR, %d frags but only %d ranked\n", n, m);
    printf("1= %d %ld %ld\n", tmp[0].chr, tmp[0].ini, tmp[0].end);
    printf("%d= %d %ld %ld\n", m,tmp[m-1].chr,tmp[m-1].ini,tmp[m-1].end);
    exit(8);
  }
  ftmp=tmp; struct frag *f=frags;
  for(k=0; k<n; k++){*f=*ftmp; f->next=f+1; f++; ftmp++;}
  (f-1)->next=NULL;
}

int Get_strand(char *strand)
{
  if(strand[0]=='-'){return(-1);}
  else if(strand[0]=='+'){return(1);}
  else{return(0);}
}
