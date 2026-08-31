//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_z_nkk_ch
//F	1	    ( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	magnetic field data of horizontal correctors made by nkk co. ltd.
//F!	id number is ch01 to ch40.
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
//F!	    excitation curve of all magnets were measured in this type.
//F!	    average of data shows that the slope is 108.04 and offset is 0.6g.
//F!	    we ignore this very small offset.
//F!
//F!	    b(g) = 108.0(g)/1(a) * i(a)
//F!	    b(t) = 108e-4(t)/1(a) * i(a)
//F!
//F!	(effective length)
//F!	    s-direction field distribution was measured only one magnet ch03
//F!	    by nkk co. ltd. it has no tail part.
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

int mag_ki_z_nkk_ch( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field )
{

//F	integer*4   mode, mag_ki_calc_liner
//F	real*4      energy, slope, efflen, kvalue, field, current

//F	slope  = 108.0e-4
	float slope  = 108.0e-4;
//F	efflen = 0.119210
	*efflen = 0.119210;

//F	mag_ki_z_nkk_ch = mag_ki_calc_liner
//F	1	    ( mode, energy, slope, efflen, kvalue, current, field )
	return mag_ki_calc_liner( mode, energy, slope, *efflen, kvalue, current, field );

//F	end
}

