//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_idx_skew
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	(parameter)
//F!	    coil:	i4, main coil = 1, trim coil = 2.
//F!	    mode:	i4, =1 for b to i, =2 for i to b.
//F!	    energy:	r4, energy in gev.
//F!	    efflen:	r4, effective length in meter.
//F!	    kvalue:	r4, k value.
//F!	    field:	r4, field data in t/m.
//F!	    current:	r4, magnet current in ampare.
//F!
//F!	(field data)
//F!	    measured by idx co. ltd.
//F!	    only one magnet measured in this type.
//F!	    measurement was done at x=10mm. 
//F!
//F!	(effective length)
//F!	    effective length had not been measwred,
//F!	    then we used the effective length of 0.07867,
//F!	    which was as well as that of hitachi type2
//F!	    quadrapole magnet.
//F!
//F!	parameter
//F!	    coils:  based on the drawings #bq1069b-001, 99-07-01
//F!
//F!	history
//F!	    01-nov-1999 j.ozawa, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_idx_skew( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 11;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 30
	const int ncoil_main = 30;
//F        parameter ncoil_trim = 0
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867
	const float eff_length = 0.07867;

//F	data cur(01) /  -20.00000 /, fld(01) /   -5.66200 /
//F	data cur(02) /  -16.00000 /, fld(02) /   -4.52200 /
//F	data cur(03) /  -12.00000 /, fld(03) /   -3.38500 /
//F	data cur(04) /   -8.00000 /, fld(04) /   -2.24700 /
//F	data cur(05) /   -4.00000 /, fld(05) /   -1.11500 /
//F	data cur(06) /    0.00000 /, fld(06) /    0.00000 /
//F	data cur(07) /    4.00000 /, fld(07) /    1.11500 /
//F	data cur(08) /    8.00000 /, fld(08) /    2.24700 /
//F	data cur(09) /   12.00000 /, fld(09) /    3.38500 /
//F	data cur(10) /   16.00000 /, fld(10) /    4.52200 /
//F	data cur(11) /   20.00000 /, fld(11) /    5.66200 /
//F	
	cur[ 1-1] =  -20.00000; fld[ 1-1] =   -5.66200; 
	cur[ 2-1] =  -16.00000; fld[ 2-1] =   -4.52200; 
	cur[ 3-1] =  -12.00000; fld[ 3-1] =   -3.38500; 
	cur[ 4-1] =   -8.00000; fld[ 4-1] =   -2.24700; 
	cur[ 5-1] =   -4.00000; fld[ 5-1] =   -1.11500; 
	cur[ 6-1] =    0.00000; fld[ 6-1] =    0.00000; 
	cur[ 7-1] =    4.00000; fld[ 7-1] =    1.11500; 
	cur[ 8-1] =    8.00000; fld[ 8-1] =    2.24700; 
	cur[ 9-1] =   12.00000; fld[ 9-1] =    3.38500; 
	cur[10-1] =   16.00000; fld[10-1] =    4.52200; 
	cur[11-1] =   20.00000; fld[11-1] =    5.66200; 

//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_idx_skew = mag_ki_control
//F	1		    (
//F	1		     mode,	    ! access mode
//F	1		     mx,	    ! number of data for excitation.
//F	1		     cur,	    ! current data array
//F	1		     fld,	    ! field data array
//F	1		     energy,	    ! energy in gev
//F	1		     efflen,	    ! effective length (m)
//F	1		     ncoil_main,    ! number of turn for main coil.
//F	1		     ncoil_trim,    ! number of turn for trim coil.
//F	1		     kvalue,	    ! array of k values.
//F	1		     field,	    ! array of field
//F	1		     current	    ! array of current
//F	1		    )
	return mag_ki_control(
			     mode,	    //! access mode
			     mx,	    //! number of data for excitation.
			     cur,	    //! current data array
			     fld,	    //! field data array
			     energy,	    //! energy in gev
			     *efflen,	    //! effective length (m)
			     ncoil_main,    //! number of turn for main coil.
			     ncoil_trim,    //! number of turn for trim coil.
			     kvalue,	    //! array of k values.
			     field,	    //! array of field
			     current	    //! array of current
		);

//F	end
}

