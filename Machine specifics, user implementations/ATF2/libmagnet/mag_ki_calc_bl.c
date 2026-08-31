//F!==============================================================================
//F!
//F	integer*4 function mag_ki_calc
//F	1	( mode, n, cur, fld, energy, eff_len, kvalue, field, current )
//F!
//F!	make linear interpolation from the current-field data.
//F!	
//F!	(supplied)
//F!	    mode;	integer*4   ; b_to_i = 1 or i_to_b = 2.
//F!	    n;		integer*4   ; array size of cur and fld.
//F!	    cur;	real*4	    ; array of the current data.
//F!	    fld;	real*4	    ; array of the field data.
//F!	    energy;	real*4	    ; energy in gev.
//F!	    eff_len;	real*4	    ; effective length in meter.
//F!
//F!	(supplied and returned)
//F!	    kvalue;	real*4	    ; k value.
//F!	    current;	real*4	    ; current (a).
//F!
//F!	(possible return code)
//F!	    mag_ki_calc = 0
//F!	    mag_ki_calc = 1
//F!
//F!	author: n.terunuma (kek)
//F!	
//F!	history:    8-dec-1997 first version.
//F!
//F!==============================================================================
#include <stdio.h>
#include <math.h>

int mag_ki_calc( int mode, int n, float cur[], float fld[], float energy, float eff_len, 
	float *kvalue, float *field, float *current )
{
	int i;
	const int mode_b_to_i = 1;
	const int mode_i_to_b =2;
	double xin, yout, x[n], y[n];

	if( eff_len <= 0.0 ) return(0);
	if( energy <= 0.0 ) return(0);

	if( mode == mode_b_to_i ){
	    *field = energy * fabsf(*kvalue) / eff_len / (double)0.3;
	    xin = *field;

		if( mag_fld_convert(mode,ndata,cur,fld,ncoil[2],field,current) )


 )	// Array of current

	    for( i=0; i<n; i++ ){
		x[i] = fld[i];
		y[i] = cur[i];
	    }

	}else if( mode  == mode_i_to_b ) {
	    xin = *current;
	    for( i=0; i<n; i++ ){
		x[i] = cur[i];
		y[i] = fld[i];
	    }

	}else {
	    return(0);
	}

	if( xin == x[0] ){
	    yout = y[0];
	    goto G8000;
	}else if( xin < x[0] ){
	    yout = y[0];
	    goto G9000;
	}

	for( i=0; i<(n-1); i++ ){
	    if( (xin > x[i]) && (xin <= x[i+1] ) ) goto G1000;
	}

	yout = y[n-1];
	
	goto G9000;

G1000:
	yout = y[i] + (y[i+1]-y[i]) / (x[i+1]-x[i]) * (xin-x[i]);

G8000:

G9000:
	if( mode == mode_b_to_i ){
	     *current = yout;
	}else if( mode == mode_i_to_b ){
	     *field  = yout;
	     *kvalue = (*field) / energy * eff_len * (double)0.3;
	}

	return(1);
}
