//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_z_tecno_zh82983
//F	1		    ( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	magnet field data for horizontal correctors made by tecno elect.
//F!	id number is 58283.
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
//F!	    measured by tecno electric co. ltd.
//F!	    only one magnet measured in this type.
//F!
//F!	    b(g) = 1028(g)/10(a) * i(a)
//F!	    b(t) = 1028e-4(t)/10(a) * i(a)
//F!
//F!	(effective length)
//F!	    use data of tecno. but it has no tail part data.
//F!	    we assume this part by extending other part. the accuracy of the 
//F!           effective length is around 1%. it was estimated by changing the end
//F!           of tail data with a reasonable range.
//F!
//F!	(history)
//F!	    03-dec-2007	y.tsukada, created.
//F!
//F!==============================================================================
#include <stdio.h>
#include "magnet.h"

int mag_ki_z_tecno_zh82983( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field )
{

//F	implicit    none
//F	integer*4   mode, mag_ki_calc_liner
//F	real*4      energy, slope, efflen, kvalue, field, current

//F	slope  = 5.9e-3
	float slope  = 5.9e-3;
//F	efflen = 0.1679
	*efflen = 0.1679;

//F	mag_ki_z_tecno_zh82983 = mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )
	return mag_ki_calc_liner( mode, energy, slope, *efflen, kvalue, current, field );

//F	end
}
