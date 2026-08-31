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
#include <time.h>
#include "magnet.h"

int mag_ki_q_qd0( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
////	const int mx = 16;
	const int mx = 18;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 11;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867
	const float eff_length = 0.475;
	int new_flg = 0;
	int t_ncoil_main, t_ncoil_trim;

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

/*	cur[1-1]= 0.0;
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
*/
	cur[1-1]= 0.0;
	cur[2-1]= 50.011;
	cur[3-1]= 60.028;
	cur[4-1]= 70.086;
	cur[5-1]= 80.072;
	cur[6-1]= 90.169;
	cur[7-1]= 100.03;
	cur[8-1]= 110.09;
	cur[9-1]= 120.07;
	cur[10-1]= 130.14;
	cur[11-1]= 140.15;
	cur[12-1]= 145.22;
	cur[13-1]= 150.22;
	cur[14-1]= 170.19;
	cur[15-1]= 190.26;
	cur[16-1]= 200.09;
	cur[17-1]= 210.14;
	cur[18-1]= 220.18;

/*	fld[ 1-1]= 0.000;
	fld[ 2-1]= 0.93;
	fld[ 3-1]= 1.86;
	fld[ 4-1]= 2.79;
	fld[ 5-1]= 3.72;
	fld[ 6-1]= 4.65;
	fld[ 7-1]= 5.57;
	fld[ 8-1]= 6.5;
	fld[ 9-1]= 7.43;
	fld[10-1]= 8.36;
	fld[11-1]= 9.29;
	fld[12-1]= 10.22;
	fld[13-1]= 11.15;
	fld[14-1]= 12.08;
	fld[15-1]= 13.01;
	fld[16-1]= 13.94;
*/
	fld[ 1-1]= 0.000;
	fld[ 2-1]= 4.7354;
	fld[ 3-1]= 5.6817;
	fld[ 4-1]= 6.6332;
	fld[ 5-1]= 7.5776;
	fld[ 6-1]= 8.5317;
	fld[ 7-1]= 9.463;
	fld[ 8-1]= 10.411;
	fld[ 9-1]= 11.35;
	fld[10-1]= 12.296;
	fld[11-1]= 13.236;
	fld[12-1]= 13.709;
	fld[13-1]= 14.179;
	fld[14-1]= 16.079;
	fld[15-1]= 17.955;
	fld[16-1]= 18.869;
	fld[17-1]= 19.802;
	fld[18-1]= 20.731;

	//F        efflen = eff_length
        *efflen = eff_length;
	t_ncoil_main = ncoil_main;
	t_ncoil_trim = ncoil_trim;
// 2023.10.31 -
	time_t t1 = 1698731902;// 2023.10.31 14:58
	// check 2023.10.31 before or after.
	time_t t = time(NULL);

	if( t1 < t ) new_flg = 1;
	
	const int n_mx = 5;
	const float n_eff_length = 0.48005;

	int t_mx;

	if( new_flg == 1 ){
	    cur[0] = 0;
	    cur[1] = 50.0;
	    cur[2] = 100.0;
	    cur[3] = 150.0;
	    cur[4] = 200.0;

	    fld[0] = 0.0000;
	    fld[1] = 3.1247;
	    fld[2] = 6.2493;
	    fld[3] = 9.3740;
	    fld[4] = 12.4987;

	    *efflen = n_eff_length;
	    t_mx    = n_mx;
	    t_ncoil_main = 18;
	    t_ncoil_trim = 20;
	}else{

	    t_mx = mx;
	}

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
			     t_mx,	    //! number of data for excitation.
			     cur,	    //! current data array
			     fld,	    //! field data array
			     energy,	    //! energy in gev
			     *efflen,	    //! effective length (m)
			     t_ncoil_main,    //! number of turn for main coil.
			     t_ncoil_trim,    //! number of turn for trim coil.
			     kvalue,	    //! array of k values.
			     field,	    //! array of field
			     current	    //! array of current
		);
//F	end
}

