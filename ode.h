#include <stdio.h>



typedef struct ode{
	float a,b,ua,ub;
	float (*p)(float), (*q)(float), (*f)(float);
	int N;
	float *x, *u, *S, ***C3, ***C4, **C1, **C2, **R;
	float h;
}ODE;
