//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_hitachi_2
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
//F!	    measured by M.Masuzawa(KEK)
//F!
//F!	(effective length)
//F!
//F!	parameter
//F!
//F!	history
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_qea29( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 16;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 49;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.07867
	const float eff_length = 0.19849;

//F	data cur(01) /    0.00000 /, fld(01) /    0.00000 /
//F	data cur(02) /   10.20000 /, fld(02) /    3.83200 /
//F	data cur(03) /   20.20000 /, fld(03) /    7.64600 /
//F	data cur(04) /   30.20000 /, fld(04) /   11.49800 /
//F	data cur(05) /   40.00000 /, fld(05) /   15.20200 /
//F	data cur(06) /   50.00000 /, fld(06) /   18.92700 /
//F	data cur(07) /   60.20000 /, fld(07) /   22.65700 /
//F	data cur(08) /   70.20000 /, fld(08) /   26.28500 /
//F	data cur(09) /   80.20000 /, fld(09) /   29.89600 /
//F	data cur(10) /   90.20000 /, fld(10) /   33.43800 /
//F	data cur(11) /  100.20000 /, fld(11) /   36.75800 /
//F	

	cur[1-1]= 0.000;	
	cur[2-1]= 3.333;
	cur[3-1]= 6.667;
	cur[4-1]= 10.000;
	cur[5-1]= 13.333;
	cur[6-1]= 16.667;
	cur[7-1]= 20.000;
	cur[8-1]= 23.333;
	cur[9-1]= 26.667;
	cur[10-1]= 30.000;
	cur[11-1]= 33.333;
	cur[12-1]= 36.667;
	cur[13-1]= 40.000;
	cur[14-1]= 43.333;
	cur[15-1]= 46.667;
	cur[16-1]= 50.000;

	fld[ 1-1]= 0.000; 
	fld[ 2-1]= 1.604; 
	fld[ 3-1]= 3.136; 
	fld[ 4-1]= 4.682; 
	fld[ 5-1]= 6.237; 
	fld[ 6-1]= 7.790; 
	fld[ 7-1]= 9.341; 
	fld[ 8-1]= 10.891; 
	fld[ 9-1]= 12.479; 
	fld[10-1]= 14.023; 
	fld[11-1]= 15.560; 
	fld[12-1]= 17.098; 
	fld[13-1]= 18.634; 
	fld[14-1]= 20.165; 
	fld[15-1]= 21.694; 
	fld[16-1]= 23.230;


//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_hitachi_2 = mag_ki_control
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

