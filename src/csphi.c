#include <math.h>
#include <stdlib.h>
#include "csphi.h"


#define B(n,v) (((n)-(v))*((n)+(v)+1.0))

double ynnm(int n, int m)
{
	int pm,k;
	double a,ynnm;

	a=1.0/sqrt(4.0*PI);
	pm = abs(m);

	if(n < pm){
		ynnm=0;
	}
	else if (n==0){
		ynnm=a;
	}
	else{
		ynnm=a;
		for(k=1; k<=n; k++){
			ynnm=sqrt((2.0*k+1)/(8.0*k))*ynnm;
		}
		if(n!=pm){
			for(k=n-1; k>=pm; k--){
			ynnm=sqrt(1.0*(n+k+1.0)/(n-k))*ynnm;
			}
		}
	}
	return ynnm;
}

void ynunm(int en,int em,SFLOAT* y,int len)
{
	int k;

	for(k=0;k<len;k++)
		*(y + k)=0;

	if(abs(em) <= en){
		*(y + en) = ynnm(en,em);
		k=en-2;
		if(k >= 0){
			*(y + k) = (B(en,k+1.0)+B(en,k+2.0) - 4.0*em*em)*(*(y + k + 2))/B(en,k);	
			for(k=k-2;k>=0;k-=2){
				(*(y + k))=((B(en,k+1.0)+B(en,k+2.0)-4.0*em*em)*(*(y + k + 2))-B(en,k+3.0)*(*(y + k + 4)))/B(en,k);
			}
		}
	}
}
