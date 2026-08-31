//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_z_nkk_cv
//F	1           ( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	magnetic field data of vertical correctors made by nkk co. ltd.
//F!	id number is cv**.
//F!
//F!	(parameter)
//F!	    mode:	i4, =1 for b to i, =2 for i to b.
//F!	    energy:	r4, energy in gev.
//F!	    efflen:	r4, effective length in meter.
//F!	    kvalue:	r4, k value.
//F!	    field:	r4, field data in tesla.
//F!	    current:	r4, magnet current in ampare.
//F!
//F!	(field data)
//F!	    measured by nkk co. ltd.
//F!	    only three magnets measured in this type.
//F!	    average of data shows that the slope is 112.01 and offset is 1.7g.
//F!	    we ignore this very small offset.
//F!
//F!	    b(g) = 112.0(g)/1(a) * i(a)
//F!	    b(t) = 112e-4(t)/1(a) * i(a)
//F!
//F!	(effective length)
//F!	    use data of nkk but it has no tail part data. cv05
//F!	    we assume this part by extending other part. the accuracy of the 
//F!	    effective length is less than ??%.
//F!
//F!	(history)
//F!	    13-mar-1998	n.terunuma, created. not finished, (efflen)
//F!
//F!==============================================================================
#include <stdio.h>
#include "magnet.h"

int mag_ki_z_nkk_cv( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field )
{

//F	implicit    none
//F	integer*4   mode, mag_ki_calc_liner
//F	real*4      energy, slope, efflen, kvalue, field, current
	float slope;

//F	slope  = 112.0e-4
	slope = 112.0e-4;
//F	efflen = 0.1248
	*efflen = 0.1248;

//F	mag_ki_z_nkk_cv = mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )

	return mag_ki_calc_liner( mode, energy, slope, *efflen, kvalue, current, field );
//F	end
}
