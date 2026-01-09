#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Peaks_aux.h"
#include "cluster_peaks.h"
#include "allocate.h"


static int Rank_clusters(int *cluster_out, int *nc, double **Mean_out,
			 double ***S2inv, double *logdet, float *logtau,
			 int *cluster, double **Mean, double ***Sig2,
			 int N, int nclus, int d, int SIGMA);

int Cluster_peaks(struct peak *peaks, int N_peak, char *name, int makeclust){
  // Output: cluster number peak->c
  // Input: peak->y, peak->yout

  char nameout[200]; strcpy(nameout, name); strcat(nameout, ".clus");
  FILE *file_out=fopen(nameout, "w");
  printf("Writing clustering results in file %s\n", nameout);

  int i, clus_opt[N_peak];;
  float **x=malloc(N_peak*sizeof(float *));
  for(i=0; i<N_peak; i++){
    x[i]=malloc(2*sizeof(float));
    x[i][0]=peaks[i].yout;
    x[i][1]=peaks[i].y;
  }
  int k_min=1, k_max=8;
  int cluster[N_peak];
  int nc[k_max];
  float S_min=100000, score;
  double *Mean[k_max];
  for(int k=0; k<k_max; k++){Mean[k]=malloc(2*sizeof(double));}
  for(int d=1; d<=2; d++){
    for(int SIGMA=1; SIGMA>=-1; SIGMA--){
      for(int k=k_min; k<=k_max; k++){
	printf("nclus= %d d= %d SIGMA= %d\n", k, d, SIGMA);
	if(SIGMA>=0){
	  score=EM(cluster, nc, Mean, k, d, SIGMA, x, N_peak);
	}else{
	  score=K_means(k, d, cluster, nc, Mean, x, N_peak); score=-score;
	}
	int opt=0;
	//if(silhouette>S_max){S_max=silhouette;} //opt=1;
	if(score<S_min){S_min=score; opt=1;}
	if(k==5 && d==1 && SIGMA==1){opt=1;}else{opt=0;}
	if(opt){for(i=0; i<N_peak; i++){clus_opt[i]=cluster[i];}}
	fprintf(file_out, "d= %d SIGMA=%d k= %d ", d, SIGMA, k);
	if(SIGMA>=0){fprintf(file_out, "bic");}
	else{fprintf(file_out, "sil"); score=-score;}
	fprintf(file_out, "= %.4g clus:", score);
	for(int kk=0; kk<k; kk++){
	  fprintf(file_out, " %d (%.2g", nc[kk], Mean[kk][0]);
	  for(int j=1; j<d; j++)fprintf(file_out, ",%.2g", Mean[kk][j]);
	  fprintf(file_out, ");");
	}
	fprintf(file_out, "\n");
      }
    }
  }

  for(i=0; i<N_peak; i++){free(x[i]);}
  if(makeclust){
    for(i=0; i<N_peak; i++){peaks[i].c=clus_opt[i]+1;}
   }
  fclose(file_out);
}


float EM(int *cluster_out, int *nc, double **Mean_out, //output
	 int nclus, int d, int SIGMA, float **x, int N) //input
{
  // Uses the fact that dimension d can be only 1 or 2
  int k, i, j, j2;
  float bic=0, log2pi=log(6.283185);

  // Compute mean
  float x0max[d], x0min[d], x1[d], x2[d];
  for(j=0; j<d; j++){
    x1[j]=0; x2[j]=0;
    for(i=0; i<N; i++){
      float *xi=x[i];
      x1[j]+=xi[j]; x2[j]+=xi[j]*xi[j];
      if(j==0){
	if(i==0){x0min[j]=xi[j]; x0max[j]=xi[j];}
	else if(xi[j]<x0min[j]){x0min[j]=xi[j];}
	else if(xi[j]>x0max[j]){x0max[j]=xi[j];}
      }
    }
    x1[j]/=N;
    x2[j]=sqrt((x2[j]-N*x1[j]*x1[j])/(N-1));
  }
    
  // Initialize clusters
  float tau[nclus], tau_new[nclus], logtau[nclus];
  int nc_new[nclus], cluster[N], cluster_new[N];
  double **Sig2[nclus], **S2inv[nclus], *S2kj, **Sig2_new[nclus],
    *Mean_new[nclus], *Mean[nclus], logdet[nclus];

  for(k=0; k<nclus;  k++){
    nc[k]=0;
    nc_new[k]=0;
    Mean[k]=malloc(d*sizeof(double));
    Mean_new[k]=malloc(d*sizeof(double));
    if(SIGMA==0){
      Sig2[k]=NULL; S2inv[k]=NULL; Sig2_new[k]=NULL;
    }else{
      Sig2[k]=Allocate_mat2_d(d, d);
      S2inv[k]=Allocate_mat2_d(d, d);
      Sig2_new[k]=Allocate_mat2_d(d, d);
    }
  }
  float peq=1/(float)nclus, lpeq=log(peq);
  for(k=0; k<nclus; k++){
    tau_new[k]=peq; logtau[k]=lpeq;
  }
  for(j=0; j<d; j++){
    float step=(x0max[j]-x0min[j])/nclus, yy=x0min[j]+step/2;
    for(k=0; k<nclus; k++){
      Mean_new[k][j]=yy; yy+=step; if(SIGMA)Sig2_new[k][j][j]=1;
    }
  }
  
  float eps=0.000001; int iter, itmax=100;
  double negloglik_prev=-N*100000;
  for(iter=0; iter<itmax; iter++){

    // Initialize mean, variance and size of clusters
    for(k=0; k<nclus;  k++){
      tau[k]=tau_new[k];
      tau_new[k]=0;
      if(tau[k]>0){logtau[k]=log(tau[k]);}
      nc[k]=nc_new[k];
      nc_new[k]=0;
      float det=0;
      for(j=0; j<d; j++){
	Mean[k][j]=Mean_new[k][j];
	Mean_new[k][j]=0;
      }
      if(SIGMA){
	double **S2k=Sig2[k], **S2invk=S2inv[k], **S2k_new=Sig2_new[k];
	for(j=0; j<d; j++){
	  for(j2=0; j2<=j; j2++){
	    S2k[j][j2]=S2k_new[j][j2];
	    S2k_new[j][j2]=0;
	  }
	}
	if(d==1){
	  det=S2k[0][0];
	  S2invk[0][0]=1/S2k[0][0];
	}else if(d==2){
	  det=S2k[0][0]*S2k[1][1]-S2k[1][0]*S2k[1][0];
	  S2invk[0][0]=S2k[1][1]/det;
	  S2invk[1][0]=-2*S2k[1][0]/det; // double, to sum only on j2<=j
	  S2invk[1][1]=S2k[0][0]/det;
	}
	// General case: compute determinant and inverse
	if(det>0){logdet[k]=log(det);}
      }
    }

    // Compute mean and variance
    double negloglik_all=0;
    for(i=0; i<N; i++){
      float pi[nclus], lp[nclus];
      float *xi=x[i];
      cluster[i]=cluster_new[i];
      cluster_new[i]=-1;
      int kk=-1; float lpmin=100000, Z=0;
      for(k=0; k<nclus; k++){
	float diff[d], negloglik=0;
	if(SIGMA){negloglik+=logdet[k];}
	for(j=0; j<d; j++){
	  diff[j]=Mean[k][j]-x[i][j];
	  if(SIGMA){
	    double *Sinvkj=S2inv[k][j];
	    for(j2=0; j2<=j; j2++){
	      negloglik+=diff[j]*Sinvkj[j2]*diff[j2];
	    }
	  }else{
	    negloglik+=diff[j]*diff[j];
	  }
	}
	negloglik=0.5*negloglik-logtau[k];
	if(k==0 || negloglik<lpmin){lpmin=negloglik; kk=k;}
	lp[k]=negloglik;
	pi[k]=exp(-negloglik); Z+=pi[k];
      }
      cluster_new[i]=kk;
      nc_new[kk]++;
      if(Z<=0 || isnan(Z)){
	printf("WARNING, i=%d Z=%.2g\n", i, Z); Z=0;
	for(k=0; k<nclus; k++){lp[k]-=lpmin; pi[k]=exp(-lp[k]); Z+=pi[k];}
	if(Z==0 || isnan(Z)){printf("ERROR, i=%d Zk=%.2g\n", i, Z); goto end;}
      }

      //Z*=N; // sum_i p_i/Z= 1/N
      for(k=0; k<nclus; k++){
	pi[k]/=Z;
	tau[k]+=pi[k];
	negloglik_all+=pi[k]*lp[k];
	for(j=0; j<d; j++){
	  float px=pi[k]*xi[j];
	  Mean_new[k][j]+=px;
	  if(SIGMA){
	    S2kj=Sig2_new[k][j];
	    for(j2=0; j2<=j; j2++)S2kj[j2]+=px*xi[j2];
	  }
	}
      } //end k
    } // end i
    // Normalization
    for(k=0; k<nclus; k++){
      for(j=0; j<d; j++){
	Mean_new[k][j]/=tau[k];
	if(SIGMA){
	  double *S2kj=Sig2_new[k][j];
	  for(j2=0; j2<=j; j2++){
	    S2kj[j2]=S2kj[j2]/tau[k]-Mean_new[k][j]*Mean_new[k][j2];
	  }
	}
      } // end j
      tau[k]/=N;
    } // end k
    // Test convergence
    for(k=0; k<nclus; k++){
      if(nc[k]!=nc_new[k])break;
      for(j=0; j<d; j++){
	if(fabs(Mean[k][j]-Mean_new[k][j])>eps)break;
      }
    }
    if(k==nclus){
      printf("Convergence in %d iterations\n", iter); break;
    }

    negloglik_all=negloglik_all+0.5*d*log2pi;
    printf("%.4g\n", negloglik_all);
    if(isnan(negloglik_all)){printf("nan, exiting\n"); break;}
    if(iter && negloglik_all>negloglik_prev){
      printf("lik down, exiting\n");break;
    }
    negloglik_prev=negloglik_all;
  }
  if(iter==itmax)
    printf("WARNING, EM did not converge after %d steps\n",iter);


  // Order clusters based on x[0]
  Rank_clusters(cluster_out, nc, Mean_out, S2inv, logdet, logtau,
		cluster, Mean, Sig2, N, nclus, d, SIGMA);

  double negloglik=0;
  for(i=0; i<N; i++){
    int k=cluster_out[i];
    float diff[d];
    for(j=0; j<d; j++){
      diff[j]=Mean_out[k][j]-x[i][j];
      if(SIGMA){
	double *Sinvkj=S2inv[k][j];
	for(j2=0; j2<=j; j2++){
	  negloglik+=diff[j]*Sinvkj[j2]*diff[j2];
	}
      }else{
	negloglik+=diff[j]*diff[j];
      }
    }
    if(SIGMA){negloglik+=logdet[k];}
  }
  int nclus_act=0; for(k=0; k<nclus;  k++){if(nc[k])nclus_act++;}
  int param=nclus_act*(d+1);
  if(SIGMA){
    param+=nclus_act*d*(d+1)/2;
    negloglik+=N*d*log2pi;
  }
  float logdetmin=0;
  for(k=0; k<nclus;  k++){
    if(nc[k]==0)continue;
    if(logdetmin==0 || logdet[k]<logdetmin)logdetmin=logdet[k];
  }
  negloglik-=N*logdetmin;
  bic=negloglik+param*log((float)N);

 end:
  for(k=0; k<nclus;  k++){
    free(Mean_new[k]);
    free(Mean[k]);
    if(SIGMA){
      Empty_matrix_d(Sig2[k],d);
      Empty_matrix_d(S2inv[k],d);
      Empty_matrix_d(Sig2_new[k],d);
    }
  }

  return(bic/N);
}
float K_means(int nclus, int d, int *cluster, int *nc, double **Mean,
	      float **x, int N)
	    //int *nc, double *Mean, double *Sd, float *ythr_clus, int *cluster,
	    //float *y_scr, int nk, int NC)
{
  // Use the fact that dimension d can be only 1 or 2
  int k, i, j;

  // Compute mean
  float x1[d], x2[d], x0max=-100, x0min=100;
  for(j=0; j<d; j++){
    x1[j]=0; x2[j]=0;
    for(i=0; i<N; i++){
      float *xi=x[i];
      x1[j]+=xi[j]; x2[j]+=xi[j]*xi[j];
      if(j==0){
	if(xi[0]<x0min){x0min=xi[0];}else if(xi[0]>x0max){x0max=xi[0];}
      }
    }
    x1[j]/=N;
    x2[j]=sqrt((x2[j]-N*x1[j]*x1[j])/(N-1));
  }
    
  // Initialize clusters
  int error=0;
  float x_thr[nclus], step=(x0max-x0min)/nclus, yy=x0min;
  for(k=0; k<nclus; k++){x_thr[k]=yy; yy+=step;}
  for(i=0; i<N; i++){
    for(k=nclus-1; k>=0; k--){
      if(x[i][0]>=x_thr[k])break;
    }
    if(k<0){
      printf("WARNING, %.3f not classified min: %.3f\n",x[i][0], x0min);
      error++; k=0;
    }
    cluster[i]=k;
  }
  if(error){
    printf("%d errors in initialization, leaving\n", error); exit(8);
  }

  int nc_old[nclus];
  double *Mean_old[nclus];
  for(k=0; k<nclus;  k++){
    nc_old[k]=0;
    Mean_old[k]=malloc(d*sizeof(double));
    for(j=0; j<d; j++){Mean_old[k][j]=0;}
  }

  float eps=0.00001; int iter, itmax=100;
  for(iter=0; iter<itmax; iter++){
    // Compute mean
    for(k=0; k<nclus;  k++){
      nc[k]=0; for(j=0; j<d; j++){Mean[k][j]=0;}
    }
    for(i=0; i<N; i++){
      k=cluster[i];
      if(k<0 || k>=nclus){
	printf("ERROR, wrong cluster index %d max: %d", k, nclus-1);
	exit(8);
      }
      nc[k]++; float *xi=x[i];
      for(j=0; j<d; j++){Mean[k][j]+=xi[j];}
    }
    for(k=0; k<nclus; k++){
      if(nc[k]){
	for(j=0; j<d; j++){Mean[k][j]/=nc[k];}
      }
    }

    // Test convergence
    for(k=0; k<nclus; k++){
      if(nc[k]!=nc_old[k])break;
      for(j=0; j<d; j++){if(fabs(Mean[k][j]-Mean_old[k][j])>eps)break;}
    }
    if(k==nclus){
      printf("Convergence in %d iterations\n", iter); break;
    }
    for(k=0; k<nclus; k++){
      nc_old[k]=nc[k]; nc[k]=0;
      for(j=0; j<d; j++)Mean_old[k][j]=Mean[k][j];
    }

    // Set clusters
    for(i=0; i<N; i++){
      int kk=-1; float d2min=1000;
      for(k=0; k<nclus; k++){
	float d2=0;
	for(j=0; j<d; j++){
	  float dd=Mean[k][j]-x[i][j]; d2+=dd*dd;
	}
	if(k==0 || d2<d2min){d2min=d2; kk=k;}
      }
      if(kk<0){
	printf("ERROR, wrong cluster %d\n", kk); exit(8);
      }
      cluster[i]=kk;
      nc[kk]++;
    }
  }
  if(iter==itmax)
    printf("WARNING, K_means did not converge after %d steps\n",iter);

  int cluster_old[N]; for(i=0; i<N; i++){cluster_old[i]=cluster[i];}
  Rank_clusters(cluster, nc, Mean, NULL, NULL, NULL,
		cluster_old, Mean_old, NULL, N, nclus, d, 0);

  for(k=0; k<nclus;  k++){free(Mean_old[k]);}

  float silhouette=0;
  for(i=0; i<N; i++){
    int kk=cluster[i]; if(nc[kk]<=1){continue;} // silhouette=0
    float dist[nclus]; for(k=0; k<nclus; k++){dist[k]=0;}
    float *xi=x[i];
    for(int i2=0; i2<N; i2++){
      if(i2==i)continue;
      int k2=cluster[i2];
      float d2=0;
      for(j=0; j<d; j++){
	float dd=x[i2][j]-xi[j]; d2+=dd*dd;
      }
      dist[k2]+=sqrt(d2);
    }
    float dmin=1000; int j=-1;
    for(k=0; k<nclus; k++){
      if(k==kk || nc[k]<=0 ){continue;}
      dist[k]/=nc[k]; if(j<0 || dist[k]<dist[j]){j=k;}
    }
    if(nc[kk]>1)dist[kk]/=(nc[kk]-1);
    if(j>=0){
      float sil=(dist[j]-dist[kk]);
      if(dist[j]>dist[kk]){sil/=dist[j];}else{sil/=dist[kk];}
      silhouette+=sil;
    }
  }
  return(silhouette/N);
}

int Rank_clusters(int *cluster_out, int *nc, double **Mean_out,
		  double ***S2inv, double *logdet, float *logtau,
		  int *cluster, double **Mean, double ***Sig2,
		  int N, int nclus, int d, int SIGMA)
{
  int i, j, k;
  
  int copied[nclus];
  for(k=0; k<nclus;  k++){copied[k]=0; nc[k]=0;}
  double logtau_old[nclus];
  if(logtau){for(k=0; k<nclus;  k++)logtau_old[k]=logtau[k];}
  for(i=0; i<N; i++){cluster_out[i]=-1;}
  for(k=0; k<nclus;  k++){
    int kk=-1;
    for(j=0; j<nclus; j++){
      if(copied[j])continue;
      if(kk<0 || Mean[j][0]<Mean[kk][0]){kk=j;}
    }
    if(kk<0){
      printf("ERROR ranking clusters, rank %d not found\n", k); exit(8);
    }
    copied[kk]=1;
    for(i=0; i<N; i++){
      if(cluster[i]==kk){cluster_out[i]=k; nc[k]++;}
    }
    for(j=0; j<d; j++){
      Mean_out[k][j]=Mean[kk][j];
    }
    if(logtau){logtau[k]=logtau_old[kk];}
    if(SIGMA){
      double **S2k=Sig2[kk], **S2invk=S2inv[k], det=0;
      if(d==1){
	det=S2k[0][0];
	S2invk[0][0]=1/S2k[0][0];
      }else if(d==2){
	det=S2k[0][0]*S2k[1][1]-S2k[1][0]*S2k[1][0];
	S2invk[0][0]=S2k[1][1]/det;
	S2invk[1][0]=-2*S2k[1][0]/det; // double, to sum only on j2<=j
	S2invk[1][1]=S2k[0][0]/det;
      }
      // General case: compute determinant and inverse
      /*for(j=0; j<d; j++){
	for(j2=0; j2<=j; j2++){
	  S2k_new[j][j2]=S2k[j][j2];
	}
	}*/
      if(det>0){logdet[k]=log(det);}
    }
  }
  return(0);
}
