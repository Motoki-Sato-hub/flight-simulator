//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_hitachi_3
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
//F!	    measured by hitachi co. ltd.
//F!	    only one magnet measured in this type.
//F!	    measurement was done at x=10mm. 
//F!
//F!	(effective length)
//F!	    integrated field strength = 844234 gauss*mm.
//F!	    field strength at z=0 is 4253.61 gauss*mm, then eff=0.19847 meter.
//F!
//F!	parameter
//F!	    coils:  based on the drawings #310q431-359, 96-01-11.
//F!
//F!	history
//F!	    13-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_hitachi_3( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 11;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 29
	const int ncoil_main = 29;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.19847
	const float eff_length = 0.19847;
//F	
//F	data cur(01) /    0.00000 /, fld(01) /    0.00000 /
//F	data cur(02) /   20.40000 /, fld(02) /    5.62800 /
//F	data cur(03) /   40.00000 /, fld(03) /   11.10900 /
//F	data cur(04) /   60.00000 /, fld(04) /   16.65600 /
//F	data cur(05) /   80.20000 /, fld(05) /   22.24300 /
//F	data cur(06) /  100.80000 /, fld(06) /   27.86600 /
//F	data cur(07) /  120.20000 /, fld(07) /   33.12200 /
//F	data cur(08) /  140.00000 /, fld(08) /   38.48700 /
//F	data cur(09) /  160.20000 /, fld(09) /   43.82800 /
//F	data cur(10) /  180.20000 /, fld(10) /   48.57500 /
//F	data cur(11) /  200.20000 /, fld(11) /   51.97100 /
//F	
	cur[ 1-1] =    0.00000; fld[ 1-1] =    0.00000; 
	cur[ 2-1] =   20.40000; fld[ 2-1] =    5.62800; 
	cur[ 3-1] =   40.00000; fld[ 3-1] =   11.10900; 
	cur[ 4-1] =   60.00000; fld[ 4-1] =   16.65600; 
	cur[ 5-1] =   80.20000; fld[ 5-1] =   22.24300; 
	cur[ 6-1] =  100.80000; fld[ 6-1] =   27.86600; 
	cur[ 7-1] =  120.20000; fld[ 7-1] =   33.12200; 
	cur[ 8-1] =  140.00000; fld[ 8-1] =   38.48700; 
	cur[ 9-1] =  160.20000; fld[ 9-1] =   43.82800; 
	cur[10-1] =  180.20000; fld[10-1] =   48.57500; 
	cur[11-1] =  200.20000; fld[11-1] =   51.97100; 

//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_hitachi_3 = mag_ki_control
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
	//printf("mode' = %d\n", mode );
	//printf("energy' = %f\n", energy );
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

