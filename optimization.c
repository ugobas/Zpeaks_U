#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "optimization.h"

float Find_max_quad(float *x, float *y,
		     float MIN, float MAX)
{
  // Rank x 
  int i0=-1, i1=-1, i2=-1;
  if(x[0]<=x[1] && x[0]<=x[2]){
    i0=0; if(x[1]<=x[2]){i1=1; i2=2;}else{i1=2; i2=1;}
  }else if(x[1]<=x[0] && x[1]<=x[2]){
    i0=1; if(x[0]<=x[2]){i1=0; i2=2;}else{i1=2; i2=0;}
  }else if(x[2]<=x[0] && x[2]<=x[1]){
    i0=2; if(x[0]<=x[1]){i1=0; i2=1;}else{i1=1; i2=0;}
  }
  if(i0<0 || i1<0 || i2<0){
    printf("ERROR in Find_max_quad2, x= %.2g %.2g %.2g i0=%d i1=%d i2=%d\n",
	   x[0], x[1], x[2], i0, i1, i2); exit(8);
  }

  float xnew=0;
  // Minimum
  if((y[i1]<y[i0])&&(y[i1]<y[i2])){
    if(y[i0]>y[i2]){
      xnew=2*x[i0]-x[i1]; goto done;
      //return((x[i0]+MIN)/2);
    }else{
      xnew=2*x[i2]-x[i1]; goto done;
      //return((x[i2]+MAX)/2);
    }
  }
  if(x[0]==x[1] || x[0]==x[2] || x[1]==x[2]){
    int imin=0;
    if(y[1]<y[imin]){imin=1;}
    if(y[2]<y[imin]){imin=2;}
    xnew=(x[0]+x[1]+x[2]-x[imin])/2;
    goto done;
  }

  float yprime_1=((y[i0]-y[i1])/(x[i0]-x[i1]));
  float yprime_2=((y[i1]-y[i2])/(x[i1]-x[i2]));
  float a=(yprime_1-yprime_2)/(x[i0]-x[i2]);
  if(a>=0){
    // Second derivative is positive, no maximum!
    if(y[i0]>y[i2]){
      xnew=2*x[i0]-x[i1]; goto done;
      //return((x[i0]+MIN)/2);
    }else{
      xnew=2*x[i2]-x[i1]; goto done;
      //return((x[i2]+MAX)/2);
    }
  }
  float b=yprime_1-a*(x[i0]+x[i1]);
  xnew=-b/(2*a);
 done:
  if(xnew<MIN){xnew=MIN;}
  else if(xnew>MAX){xnew=MAX;}
  return(xnew);
}

int Rearrange_points(float *x, float *y,
		     float x0, float y0)
{
  int i0=0; // Find minimum value of y
  if(y[1]<y[i0]){i0=1;}
  if(y[2]<y[i0]){i0=2;}
  if(y0<=y[i0]){
    // Check if x0 is in the interval. Substitute with closest extremum 
    int imin=0, imax=0;
    if(x[1]<x[imin]){imin=1;}
    else if(x[1]>x[imax]){imax=1;}
    if(x[2]<x[imin]){imin=2;}
    else if(x[2]>x[imax]){imax=2;}
    if(x0<=x[imin] || x0>=x[imax]){
      printf("y=%g < %g not inserted because x=%.3g not in %.3g-%.3g\n",
	     y0, y[i0], x0, x[imin], x[imax]);
      return(-1);
    }
    if((x0-x[imin])<(x[imax]-x0)){i0=imin;}
    else{i0=imax;}
  }
  if(fabs(x0-x[i0])<0.002*fabs(x0)){
    printf("y=%g not inserted because x=%g ~ %g y= %g\n",
	   y0, x0, x[i0], y[i0]);
    return(-1);
  }
  x[i0]=x0; y[i0]=y0;
  return(0);

  /*if(x0<x[0]){
    x[2]=x[1]; y[2]=y[1];
    x[1]=x[0]; y[1]=y[0];
    x[0]=x0;   y[0]=y0;
  }else if(x0<x[1]){
    x[2]=x[1]; y[2]=y[1];
    x[1]=x0;   y[1]=y0;
  }else if(x0<x[2]){
    x[0]=x[1]; y[0]=y[1];
    x[1]=x0;   y[1]=y0;
  }else{
    x[0]=x[1]; y[0]=y[1];
    x[1]=x[2]; y[1]=y[2];
    x[2]=x0;  y[2]=y0;
    }*/
}


float Find_max_quad_old(float x1, float x2, float x3,
			float y1, float y2, float y3,
			float MIN, float MAX)
{
  if((y2<y1)&&(y2<y3)){
    if(y1>y3){
      return((x1+MIN)/2);
    }else{
      return((x3+MAX)/2);
    }
  }
  float yprime_1=((y1-y2)/(x1-x2));
  float yprime_2=((y2-y3)/(x2-x3));
  float a=(yprime_1-yprime_2)/(x1-x3);
  if(a>0){
    // Second derivative is positive, no maximum!
    if(y1>y3){
      return((x1+MIN)/2);
    }else{
      return((x3+MAX)/2);
    }
  }
  float b=yprime_1-a*(x1+x2);
  float x0=-b/(2*a);
  if(x0<MIN){x0=MIN;}
  else if(x0>MAX){x0=MAX;}
  return(x0);
}
  
void Rearrange_points_old(float *x1, float *x2, float *x3,
			  float *y1, float *y2, float *y3,
			  float x0, float y0)
{
  if(x0<*x1){
    *x3=*x2; *y3=*y2;
    *x2=*x1; *y2=*y1;
    *x1=x0;  *y1=y0;
  }else if(x0<*x2){
    *x3=*x2; *y3=*y2;
    *x2=x0;  *y2=y0;
  }else if(x0<*x3){
    *x1=*x2; *y1=*y2;
    *x2=x0; *y2=y0;
  }else{
    *x1=*x2; *y1=*y2;
    *x2=*x3; *y2=*y3;
    *x3=x0;  *y3=y0;
  }
}



