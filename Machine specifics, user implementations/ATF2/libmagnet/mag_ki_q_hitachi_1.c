//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_hitachi_1
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	(parameter)
//F!	    mode:	i4, =1 for b to i, =2 for i to b.
//F!	    energy:	r4, energy in gev.
//F!	    kvalue:	r4, k value.
//F!	    current:	r4, magnet current in ampare.
//F!	    efflen:	r4, effective length in meter.
//F!	    field:	r4, field data in t/m.
//F!
//F!	(field data)
//F!	    measured by hitachi co. ltd.
//F!	    only one magnet measured in this type.
//F!	    measurement was done at x=10mm. 
//F!
//F!	(effective length)
//F!	    integrated field strength = 145652 gauss.mm at x=8mm offset.
//F!	    field strength at z=0 is 1849.19 gauss, then the effective length
//F!	    becomes 0.078765 meter.
//F!
//F!	documents
//F!	    drawings: #310q431-357, 96-01-11.
//F!
//F!	history
//F!	    13-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_hitachi_1( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 11;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 17
	const int ncoil_main = 17;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.078765
	const float eff_length = 0.078765;

//F	data cur( 1) /   0.0 /, fld( 1) /  0.000 /
//F	data cur( 2) /  14.2 /, fld( 2) /  2.350 /
//F	data cur( 3) /  28.6 /, fld( 3) /  4.712 /
//F	data cur( 4) /  42.2 /, fld( 4) /  6.981 /
//F	data cur( 5) /  56.0 /, fld( 5) /  9.256 /
//F	data cur( 6) /  70.2 /, fld( 6) / 11.604 /
//F	data cur( 7) /  84.2 /, fld( 7) / 13.892 /
//F	data cur( 8) /  98.2 /, fld( 8) / 16.167 /
//F	data cur( 9) / 112.2 /, fld( 9) / 18.410 /
//F	data cur(10) / 126.0 /, fld(10) / 20.629 /
//F	data cur(11) / 140.2 /, fld(11) / 22.898 /

	cur[ 1-1] =   0.0; fld[ 1-1] =  0.000; 
	cur[ 2-1] =  14.2; fld[ 2-1] =  2.350; 
	cur[ 3-1] =  28.6; fld[ 3-1] =  4.712; 
	cur[ 4-1] =  42.2; fld[ 4-1] =  6.981; 
	cur[ 5-1] =  56.0; fld[ 5-1] =  9.256; 
	cur[ 6-1] =  70.2; fld[ 6-1] = 11.604; 
	cur[ 7-1] =  84.2; fld[ 7-1] = 13.892; 
	cur[ 8-1] =  98.2; fld[ 8-1] = 16.167; 
	cur[ 9-1] = 112.2; fld[ 9-1] = 18.410; 
	cur[10-1] = 126.0; fld[10-1] = 20.629; 
	cur[11-1] = 140.2; fld[11-1] = 22.898; 

//F	
//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_hitachi_1 = mag_ki_control
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

