//F!===============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_zv100r
//F	1		    ( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	magnet field data for horizontal correctors made by ???.
//F!	id number is ???.
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
//F!	    measured by ???.
//F!
//F!	(effective length)
//F!
//F!	(history)
//F!	    05-dec-2007	created.
//F!
//F!	----------------------------------------------------------------------
//F!
//F!	        bl      /    b    =  efflen
//F!	    0.005924151 / 0.0427  = 0.138738
//F!	                          = 0.13874
//F!
//F!	        b       /  efflen =   slope
//F!	    0.0011815   / 0.13874 = 0.0085159
//F!	                          = 8.5e-3
//F!
//F!===============================================================================
#include <stdio.h>
#include "magnet.h"

int mag_ki_zv100r( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field )
{
//F	implicit    none
//F	integer*4   mode, mag_ki_calc_liner
//F	real*4      energy, slope, efflen, kvalue, field, current


//F	slope  = 8.5e-3
	float slope  = 8.5e-3;
//F	efflen = 0.13874
	*efflen = 0.13874;

//F	mag_ki_zv100r = mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )

	return mag_ki_calc_liner( mode, energy, slope, *efflen, kvalue, current, field );

//F	end
}
