#include <stdio.h>
#include <math.h>
#include <dislin.h>

void rysunek(float *x, float *y, int n) {	
	int i, ic;
	float rys_x[n], rys_y1[n], rys_y2[n];
	
	for(i=0; i<=n; i++) {
		rys_x[i] = (float) x[i];
		rys_y1[i] = (float) y[i];
		rys_y2[i] = (float) x[i]*x[i]*x[i] + 1;
	}
	
	metafl ("xwin");
	scrmod ("revers");
	disini ();
	pagera ();
	complx ();
	axspos (450, 1800);
	axslen (2200, 1200);
	
	name   ("x", "x");
	name   ("y", "y");
	
	labdig (-3, "x");
	ticks  (9, "x");
	ticks  (10, "y");
	
	titlin ("Wykres rownania:", 1);
	
	ic =   intrgb (0.95,0.95,0.95);
	axsbgd (ic);
	
	graf   (-3.0, 2.0, -3.0, 1, -14.0, 13.0, -14.0, 2);
	setrgb (0.7, 0.7, 0.7);
	grid   (1, 1);
	
	color  ("fore");
	height (50);
	title  ();
	
	color  ("red");
	curve  (rys_x, rys_y1, n+1);
	color  ("green");
	curve  (rys_x, rys_y2, n+1);
	disfin ();
}
