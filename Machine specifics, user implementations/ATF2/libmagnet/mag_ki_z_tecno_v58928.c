//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_z_tecno_v58928
//F	1		    ( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	magnet field data for vertical correctors at a wiggler section.
//F!	this magnet was made by tecno elect. id number is 58928.
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
//F!	    b(g) = 510(g)/5(a) * i(a)
//F!	    b(t) = 510e-4(t)/5(a) * i(a)
//F!
//F!	(effective length)
//F!	    use data of tecno. but it has no tail part data.
//F!	    we measured field on the vacuum chamber and normalize to that of
//F!	    tecno.
//F!	    we assume this part by extending other part. the accuracy of the 
//F!	    effective length is less than ??%.
//F!
//F!	(history)
//F!	    13-mar-1998	n.terunuma, created. not finished, (efflen)
//F!
//F!==============================================================================
//F	implicit    none
#include <stdio.h>
#include "magnet.h"

int mag_ki_z_tecno_v58928( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field )
{
	
//F	integer*4   mode, mag_ki_calc_liner
//F	real*4      energy, slope, efflen, kvalue, field, current

//F	slope  = 102.0e-4
	float slope  = 102.0e-4;
//F	efflen = 0.172
	*efflen = 0.172;

//F	mag_ki_z_tecno_v58928 = mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )
	return mag_ki_calc_liner( mode, energy, slope, *efflen, kvalue, current, field );

//F	end
}
