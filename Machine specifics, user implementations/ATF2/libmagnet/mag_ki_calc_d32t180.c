//F!==============================================================================
//F!
//F	integer*4 function mag_ki_calc_d32t180
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
//F!	    mag_ki_calc_d32t180 = 0
//F!	    mag_ki_calc_d32t180 = 1
//F!
//F!	history:    16-oct-2006 first version.
//F!
//F!==============================================================================
//F	implicit none
#include <stdlib.h>
#include "magnet.h"
#include <math.h>
#include <stdio.h>

int mag_ki_calc_d32t180( int mode, int n, float cur[], float fld[], 
	float energy, float eff_len, float *kvalue, float *field, float *current )
{
    // for compile warning
    if( eff_len < 0.0 ) printf( "eff_len<0\n");
//F	integer*4 n, i, mode
	int i;
//F	integer*4 mode_efflen / 0 /, mode_b_to_i / 1 /, mode_i_to_b / 2 /
	const int mode_b_to_i = 1;
	const int mode_i_to_b = 2;
//F	real*4    cur(n), fld(n)
//F	real*4    energy, kvalue, field, current, eff_len
//F	real*8    xin, yout, x(n), y(n)
	float xin, yout, x[n], y[n];

	int status=0;

//F	mag_ki_calc_d32t180 = 0

//F! 2006-oct-13	  eff_len not used.
//Fc	if( eff_len .le. 0.0 ) return
//F!!!

//F	if( energy  .le. 0.0 ) return
	if( energy <= 0.0 ) return(status);

//F	if( mode .eq. mode_b_to_i ) then
	if( mode == mode_b_to_i ) {

//F! 2006-oct-13	  eff_len not used.
//Fc	    field = energy * abs(kvalue) / eff_len / 0.3d0
//F	    field = energy * abs(kvalue) / 0.3d0
        *field = energy * fabs(*kvalue) / (double)0.3;
      //  float kvv = *kvalue;
      //  *field = energy * fabsf(kvv) / (double)0.3;
//F!!!

//F	    xin   = field
	    xin   = *field;
//F	    do i=1, n
	    for( i=0; i<n; i++ ){
//F		x(i) = fld(i)
		x[i] = fld[i];
//F		y(i) = cur(i)
		y[i] = cur[i];
//F	    end do
	    }
//F	else if( mode .eq. mode_i_to_b ) then
	}else if( mode == mode_i_to_b ){
//F	    xin = current
	    xin = *current;
//F	    do i=1, n
	    for( i=0; i<n; i++ ){
//F		x(i) = cur(i)
		x[i] = cur[i];
//F		y(i) = fld(i)
		y[i] = fld[i];
//F	    end do
	    }
//F	else
	}else {
//F	    return
	    return( status );
//F	end if
	}

//F	if ( xin .eq. x(1) ) then
	if( xin == x[1-1] ){
//F	    yout = y(1)
	    yout = y[1-1];
//F	    goto 8000
	    goto G8000;
//F	else if( xin .lt. x(1) ) then
	}else if( xin < x[1-1] ){
//F	    yout = y(1)
	    yout = y[1-1];
//F	    goto 9000
	    goto G9000;
//F	end if
	}

//F	do i=1, n-1
	for( i=0; i<n-1; i++ ){
//F	    if( xin.gt.x(i) .and. xin.le.x(i+1) ) goto 1000
	    if( xin > x[i] && xin <= x[i+1] ) goto G1000;
//F	end do
	}

//F	yout = y(n)
	yout = y[n-1];
//F	goto 9000
	goto G9000;

G1000:
//F1000	yout = y(i) + (y(i+1)-y(i)) / (x(i+1)-x(i)) * (xin-x(i))
	yout = y[i] + (y[i+1]-y[i]) / (x[i+1]-x[i]) * (xin-x[i]);

G8000:
//F8000	mag_ki_calc_d32t180 = 1
	status = 1;

G9000:
//F9000	if ( mode .eq. mode_b_to_i ) then
	if( mode == mode_b_to_i ){
//F	     current = yout
	     *current = yout;
//F	else if( mode .eq. mode_i_to_b ) then
	}else if( mode == mode_i_to_b ){
//F	     field  = yout
	     *field = yout;

//F! 2006-oct-13	  eff_len not used.
//Fccc	     kvalue = field / energy * eff_len * 0.3d0
//F	     kvalue = field / energy * 0.3d0
	     *kvalue = *field / energy * (double)0.3;
//F!!!
//F	end if
	}

//F	end
	return(status);
}

