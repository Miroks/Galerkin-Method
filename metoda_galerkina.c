#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "ode.h"

float pt(float x){
	//return (x*x*x+1)*cos(x);
	return (x*x+1.0f)*cos(x);
	//return 2.0f;
}

float qt(float x){
	return sin(x)-3.0f;
	//return 1.0f;
}

float ft(float x){
	//return 2+2*x*(x*x*x+1)*cos(x)+(sin(x)-5)*(x*x+1);
	return 6.0f*x+(x*x+1.0f)*3.0f*x*x*cos(x)+(sin(x)-3.0f)*(x*x*x+1.0f);
	//return x*x+4*x+3.0f;
}

void okresl_komponenty_ode(ODE *s){
	int i,j,m;
	//s->a=-1.0f;
	s->a=-2.0f;
	//s->b=3.0f;
	s->b=1.0f;
	//s->ua=2.0f;
	s->ua=-7.0f;
	//s->ub=10.0f;
	s->ub=2.0f;
	s->N=200;
	s->h=(s->b-s->a)/(s->N+1);
	s->x=vector(0,s->N+1);
	s->u=vector(0,s->N+1);
	s->S=vector(1,s->N);
	s->R=matrix(1,s->N,1,s->N);
	s->C1=matrix(0,s->N+1,0,s->N+1);
	s->C2=matrix(0,s->N+1,0,s->N+1);
	s->C3=f3tensor(0,s->N+1,0,s->N+1,0,s->N+1);
	s->C4=f3tensor(0,s->N+1,0,s->N+1,0,s->N+1);
	for(i=0;i<=s->N+1;i++){
		s->x[i]=s->a+i*s->h;
	}
	s->p=pt;
	s->q=qt;
	s->f=ft;
	for(j=0;j<=s->N+1;j++){
		for(m=0;m<=s->N+1;m++){
			if (fabsf(j-m) > 1){
				s->C1[j][m]=0.0f;
				s->C2[j][m]=0.0f;
			}
			else if(m==j){
				s->C1[j][m]=(2.0f*s->h)/3.0f;
				s->C2[j][m]=2.0f/(s->h);
			}
			else if(m!=j){
				s->C1[j][m]=(s->h)/6.0f;
				s->C2[j][m]=-1.0f/(s->h);
			}
			for(i=0;i<=s->N+1;i++){
				if ((fabsf(j-i) > 1) || (fabsf(j-m) > 1) || (fabsf(i-m) > 1)){
					s->C3[j][m][i]=0.0f;
					s->C4[j][m][i]=0.0f;
				}
				if((j==i)&&(i==m)){
					s->C3[j][m][i]=(s->h)/2.0f;
					s->C4[j][m][i]=0.0f;
				}
				if((j==m)&&(m==i-1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)1.0f/3.0f;
				}
				if((j==m)&&(m==i+1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)-1.0f/3.0f;
				}
				if((j==i)&&(i==m-1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)-1.0f/6.0f;
				}
				if((j==i)&&(i==m+1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)1.0f/6.0f;
				}
				if((j==m+1)&&(j==i+1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)-1.0f/6.0f;
				}
				if((j==m-1)&&(j==i-1)){
					s->C3[j][m][i]=(s->h)/12.0f;
					s->C4[j][m][i]=(float)1.0f/6.0f;
				}
			}
		}
	}
	
	
	s->S[1]= (*s->f)(s->x[0])*s->C1[1][1] + (*s->f)(s->x[1])*s->C1[1][1] + (*s->f)(s->x[2])*s->C1[1][2] + s->ua*(s->C2[1][0] - (*s->p)(s->x[0])*s->C4[1][0][0] - (*s->p)(s->x[1])*s->C4[1][1][0] - (*s->q)(s->x[0])*s->C3[1][0][0] - (*s->q)(s->x[1])*s->C3[1][1][0]);
	for(j=2;j<=s->N-1;j++){
		s->S[j]=(*s->f)(s->x[j-1])*s->C1[j][j-1]+(*s->f)(s->x[j])*s->C1[j][j]+(*s->f)(s->x[j+1])*s->C1[j][j+1];
	}
	s->S[s->N]=(*s->f)(s->x[s->N-1])*s->C1[s->N][s->N-1]+(*s->f)(s->x[s->N])*s->C1[s->N][s->N]+(*s->f)(s->x[s->N+1])*s->C1[s->N][s->N+1]+s->ub*(s->C2[s->N][s->N+1]-(*s->p)(s->N)*s->C4[s->N][s->N][s->N+1] - (*s->p)(s->x[s->N+1])*s->C4[s->N][s->N+1][s->N+1] - (*s->q)(s->x[s->N])*s->C3[s->N][s->N][s->N+1] - (*s->q)(s->x[s->N+1])*s->C3[s->N][s->N+1][s->N+1]);
	
	for(i=1;i<=s->N;i++){
		for(j=1;j<=s->N;j++){
			s->R[i][j]=0.0f;
		}
	}
	s->R[1][1]=-s->C2[1][1]+(*s->p)(s->x[0]) * s->C4[1][0][1]+(*s->q)(s->x[0])*s->C3[1][0][1]+(*s->p)(s->x[1])*s->C4[1][1][1]+(*s->q)(s->x[1])*s->C3[1][1][1]+(*s->p)(s->x[2])*s->C4[1][2][1]+(*s->q)(s->x[2])*s->C3[1][2][1];
	s->R[1][2]=-s->C2[1][2]+(*s->p)(s->x[1])*s->C4[1][1][2]+(*s->q)(s->x[1])*s->C3[1][1][2]+(*s->p)(s->x[2])*s->C4[1][2][2]+(*s->q)(s->x[2])*s->C3[1][2][2];
	for(j=2;j<=s->N-1;j++){
		s->R[j][j-1]=-s->C2[j][j-1]+(*s->p)(s->x[j-1])*s->C4[j][j-1][j-1]+(*s->q)(s->x[j-1])*s->C3[j][j-1][j-1]+(*s->p)(s->x[j])*s->C4[j][j][j-1]+(*s->q)(s->x[j])*s->C3[j][j][j-1];
		s->R[j][j]=-s->C2[j][j]+(*s->p)(s->x[j-1])*s->C4[j][j-1][j]+(*s->q)(s->x[j-1])*s->C3[j][j-1][j]+(*s->p)(s->x[j])*s->C4[j][j][j]+(*s->q)(s->x[j])*s->C3[j][j][j]+(*s->p)(s->x[j+1])*s->C4[j][j+1][j]+(*s->q)(s->x[j+1])*s->C3[j][j+1][j];
		s->R[j][j+1]=-s->C2[j][j+1]+(*s->p)(s->x[j])*s->C4[j][j][j+1]+(*s->q)(s->x[j])*s->C3[j][j][j+1]+(*s->p)(s->x[j+1])*s->C4[j][j+1][j+1]+(*s->q)(s->x[j+1])*s->C3[j][j+1][j+1];
	}
	s->R[s->N][s->N-1]=-s->C2[s->N][s->N-1]+(*s->p)(s->x[s->N-1])*s->C4[s->N][s->N-1][s->N-1]+(*s->q)(s->x[s->N-1])*s->C3[s->N][s->N-1][s->N-1]+(*s->p)(s->x[s->N])*s->C4[s->N][s->N][s->N-1]+(*s->q)(s->x[s->N])*s->C3[s->N][s->N][s->N-1];
	s->R[s->N][s->N]=-s->C2[s->N][s->N]+(*s->p)(s->x[s->N-1])*s->C4[s->N][s->N-1][s->N]+(*s->q)(s->x[s->N-1])*s->C3[s->N][s->N-1][s->N]+(*s->p)(s->x[s->N])*s->C4[s->N][s->N][s->N]+(*s->q)(s->x[s->N])*s->C3[s->N][s->N][s->N]+(*s->p)(s->x[s->N+1])*s->C4[s->N][s->N+1][s->N]+(*s->q)(s->x[s->N+1])*s->C3[s->N][s->N+1][s->N];
	
}


int main(void){
	ODE *s = (ODE*)malloc(2*sizeof(ODE));
	okresl_komponenty_ode(s);
	int i,j;
	int n=s->N-1;
	int *indx;
	float *d;
	indx=ivector(1,s->N);
	d=vector(1,s->N);
	ludcmp(s->R,n,indx,d);
	lubksb(s->R,n,indx,s->S);
	s->u[0]=s->ua;
	s->u[s->N+1]=s->ub;
	for(i=1;i<=s->N;i++){
		s->u[i]=s->S[i-1];
	}
	for(i=0;i<=s->N+1;i++){
		printf("x[%i]=%.3f -> %f\n",i,s->x[i],s->u[i]);
	}
	rysunek(s->x,s->u,s->N);
	free_vector(s->x,0,s->N+1);
	free_vector(s->u,0,s->N+1);
	free_vector(s->S,0,s->N-1);
	free_matrix(s->R,0,s->N-1,0,s->N-1);
	free_matrix(s->C1,0,s->N+1,0,s->N+1);
	free_matrix(s->C2,0,s->N+1,0,s->N+1);
	free_f3tensor(s->C3,0,s->N+1,0,s->N+1,0,s->N+1);
	free_f3tensor(s->C4,0,s->N+1,0,s->N+1,0,s->N+1);
	free_ivector(indx,1,s->N);
	free_vector(d,1,s->N);
}





