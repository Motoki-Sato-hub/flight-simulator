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

int mag_ki_q_qea02( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
//	const int mx = 16;
//2011/12/06 16 -> 31
//T.Yamauchi
	const int mx = 31;
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

//T.Yamauchi
/*
	cur[1-1]= 0.0;
	cur[2-1]= 10.0;
	cur[3-1]= 20.0;
	cur[4-1]= 30.0;
	cur[5-1]= 40.0;
	cur[6-1]= 50.0;
	cur[7-1]= 60.0;
	cur[8-1]= 70.0;
	cur[9-1]= 80.0;
	cur[10-1]= 90.0;
	cur[11-1]= 100.0;
	cur[12-1]= 110.0;
	cur[13-1]= 120.0;
	cur[14-1]= 130.0;
	cur[15-1]= 140.0;
	cur[16-1]= 150.0;

	fld[ 1-1]= 0.000; 
	fld[ 2-1]= 4.733;
	fld[ 3-1]= 9.372;
	fld[ 4-1]= 14.080;
	fld[ 5-1]= 18.754;
	fld[ 6-1]= 23.454;
	fld[ 7-1]= 28.081;
	fld[ 8-1]= 32.677;
	fld[ 9-1]= 37.260;
	fld[10-1]= 41.710;
	fld[11-1]= 45.856;
	fld[12-1]= 49.188;
	fld[13-1]= 51.605;
	fld[14-1]= 53.509;
	fld[15-1]= 55.080;
	fld[16-1]= 56.453;
*/
	cur[ 1-1]= -150.0;
	cur[ 2-1]= -140.0;
	cur[ 3-1]= -130.0;
	cur[ 4-1]= -120.0;
	cur[ 5-1]= -110.0;
	cur[ 6-1]= -100.0;
	cur[ 7-1]= -90.0;
	cur[ 8-1]= -80.0;
	cur[ 9-1]= -70.0;
	cur[10-1]= -60.0;
	cur[11-1]= -50.0;
	cur[12-1]= -40.0;
	cur[13-1]= -30.0;
	cur[14-1]= -20.0;
	cur[15-1]= -10.0;
	cur[16-1]= 0.0;
	cur[17-1]= 10.0;
	cur[18-1]= 20.0;
	cur[19-1]= 30.0;
	cur[20-1]= 40.0;
	cur[21-1]= 50.0;
	cur[22-1]= 60.0;
	cur[23-1]= 70.0;
	cur[24-1]= 80.0;
	cur[25-1]= 90.0;
	cur[26-1]= 100.0;
	cur[27-1]= 110.0;
	cur[28-1]= 120.0;
	cur[29-1]= 130.0;
	cur[30-1]= 140.0;
	cur[31-1]= 150.0;

	fld[ 1-1]= -56.453;
	fld[ 2-1]= -55.080;
	fld[ 3-1]= -53.509;
	fld[ 4-1]= -51.605;
	fld[ 5-1]= -49.188;
	fld[ 6-1]= -45.856;
	fld[ 7-1]= -41.710;
	fld[ 8-1]= -37.260;
	fld[ 9-1]= -32.677;
	fld[10-1]= -28.081;
	fld[11-1]= -23.454;
	fld[12-1]= -18.754;
	fld[13-1]= -14.080;
	fld[14-1]= -9.372;
	fld[15-1]= -4.733;
	fld[16-1]= 0.000;
	fld[17-1]= 4.733;
	fld[18-1]= 9.372;
	fld[19-1]= 14.080;
	fld[20-1]= 18.754;
	fld[21-1]= 23.454;
	fld[22-1]= 28.081;
	fld[23-1]= 32.677;
	fld[24-1]= 37.260;
	fld[25-1]= 41.710;
	fld[26-1]= 45.856;
	fld[27-1]= 49.188;
	fld[28-1]= 51.605;
	fld[29-1]= 53.509;
	fld[30-1]= 55.080;
	fld[31-1]= 56.453;

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

