//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )
//F!	----------------------------------------------------------------------
//F!
//F!	this routine called from mag_ki_z_***, mag_ki_q_***, ...
//F!
//F!	(parameter)
//F!	    mode:	i4, = 1	mode_k_to_i
//F!			    = 2 mode_i_to_k
//F!			    = 3 mode_b_to_i
//F!			    = 4 mode_i_to_b
//F!	    energy:	r4, energy in gev.
//F!	    slope:	r4, field/current, b(t) = slope * i(a).
//F!	    efflen:	r4, effective length in meter.
//F!	    kvalue:	r4, k value.
//F!	    field:	r4, field data in tesla... .
//F!	    current:	r4, magnet current in ampare.
//F!
//F!	(history)
//F!	    13-mar-1998	n.terunuma, created.
//F!
//F!==============================================================================
#include <stdio.h>
#include "magnet.h"

int mag_ki_calc_liner( int mode, float energy, float slope, 
	float efflen, float *kvalue, float *current, float *field )
{

//F	implicit    none
//F	integer*4   mode
//F	real*4      energy, slope, efflen, kvalue, field, current

//F	integer*4   mode_k_to_i / 1 /, mode_i_to_k / 2 /
//F	integer*4   mode_b_to_i / 3 /, mode_i_to_b / 4 /, mode_efflen / 5 /
	#define MODE_K_TO_I	1
	#define MODE_I_TO_K	2
	#define MODE_B_TO_I	3
	#define MODE_I_TO_B	4
	#define MODE_EFFLEN	5
	
//F	mag_ki_calc_liner = 0

//F	if( mode .eq. mode_k_to_i ) then
	if( mode == MODE_K_TO_I ){

//F		if( slope  .le. 0.0 ) return
		if( slope <= 0.0 ) return(0);
//F		if( efflen .le. 0.0 ) return
		if( efflen <= 0.0 ) return(0);
//F		field = energy * kvalue / efflen / 0.3d0
		*field = energy * (*kvalue) / efflen / (double)0.3;
//F		current = field / slope
		*current = (*field) / slope;

//F	else if( mode .eq. mode_i_to_k ) then
	}else if( mode == MODE_I_TO_K ){

//F		if( energy .le. 0.0 ) return
		if( energy <= 0.0 ) return(0);
//F		field = current * slope
		*field =(*current) * slope;
//F		kvalue = field / energy * efflen * 0.3d0
		*kvalue = (*field) / energy * efflen * (double)0.3;

//F	else if( mode .eq. mode_b_to_i ) then
	}else if( mode == MODE_B_TO_I ){

//F		if( slope  .le. 0.0 ) return
		if( slope  <= 0.0 ) return(0);
//F		current = field / slope
		*current = (*field) / slope;

//F	else if( mode .eq. mode_i_to_b ) then
	}else if( mode == MODE_I_TO_B ){

//F		field = current * slope
		*field = (*current) * slope;

//F	else if( mode .eq. mode_efflen ) then
	}else if( mode == MODE_EFFLEN ){
//F	else
	}else{
//F		return
		return(0);
//F	end if
	}

//F	mag_ki_calc_liner = 1
	return(1);
//F	end
}
