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

int mag_ki_q_qf1( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
////	const int mx = 11;
// 2012-11-26　奥木さんからの指示で変更
//	const int mx = 18;
//	const int mx = 14;
	const int mx = 21;

//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 11;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867

// 2012-11-26　奥木さんからの指示で変更
//	const float eff_length = 0.475;
// 2012-11-29　奥木さんからの指示で変更
//	const float eff_length = 0.43;
//	const float eff_length = 0.4441;
// 2026-03-11 奥木さんからの指示で変更
	const float eff_length = 0.48005;

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
*/
// 2012-11-26　奥木さんからの指示で変更
/*
	cur[1-1]= 0.0;
	cur[2-1]= 50.017;
	cur[3-1]= 60.025;
	cur[4-1]= 70.096;
	cur[5-1]= 80.072;
	cur[6-1]= 90.162;
	cur[7-1]= 100.03;
	cur[8-1]= 110.09;
	cur[9-1]= 120.06;
	cur[10-1]= 130.13;
	cur[11-1]= 140.15;
	cur[12-1]= 145.21;
	cur[13-1]= 150.21;
	cur[14-1]= 170.18;
	cur[15-1]= 190.26;
	cur[16-1]= 200.09;
	cur[17-1]= 210.13;
	cur[18-1]= 220.18;

	cur[ 1-1]= 0.0;
	cur[ 2-1]= 19.96615;
	cur[ 3-1]= 39.95245;
	cur[ 4-1]= 59.92775;
	cur[ 5-1]= 79.99465;
	cur[ 6-1]= 89.99305;
	cur[ 7-1]= 99.98775;
	cur[ 8-1]= 109.9852;
	cur[ 9-1]= 119.95675;
	cur[10-1]= 124.9367;
	cur[11-1]= 130.04185;
	cur[12-1]= 134.97095;
	cur[13-1]= 140.0183;
	cur[14-1]= 150.0187;
*/
/*	fld[ 1-1]= 0.0;
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
*/
// 2012-11-26　奥木さんからの指示で変更
/*
	fld[ 1-1]= 0.0;
	fld[ 2-1]= 4.73;
	fld[ 3-1]= 5.6761;
	fld[ 4-1]= 6.627;
	fld[ 5-1]= 7.5699;
	fld[ 6-1]= 8.5221;
	fld[ 7-1]= 9.4514;
	fld[ 8-1]= 10.398;
	fld[ 9-1]= 11.335;
	fld[10-1]= 12.278;
	fld[11-1]= 13.217;
	fld[12-1]= 13.69;
	fld[13-1]= 14.159;
	fld[14-1]= 16.065;
	fld[15-1]= 17.939;
	fld[16-1]= 18.852;
	fld[17-1]= 19.784;
	fld[18-1]= 20.711;
*/
// 2012-11-29　奥木さんからの指示で変更
/*
 	fld[ 1-1]= 0;
	fld[ 2-1]= 1.210960465;
	fld[ 3-1]= 2.406481395;
	fld[ 4-1]= 3.603847674;
	fld[ 5-1]= 4.807048837;
	fld[ 6-1]= 5.406362791;
	fld[ 7-1]= 6.005002326;
	fld[ 8-1]= 6.60365814;
	fld[ 9-1]= 7.199184884;
	fld[10-1]= 7.496866279;
	fld[11-1]= 7.801623256;
	fld[12-1]= 8.09494186;
	fld[13-1]= 8.394638372;
	fld[14-1]= 8.985097674;

        fld[ 1-1]= 0;
        fld[ 2-1]= 1.210960465 * 43 / 44.41;
        fld[ 3-1]= 2.406481395 * 43 / 44.41;
        fld[ 4-1]= 3.603847674 * 43 / 44.41;
        fld[ 5-1]= 4.807048837 * 43 / 44.41;
        fld[ 6-1]= 5.406362791 * 43 / 44.41;
        fld[ 7-1]= 6.005002326 * 43 / 44.41;
        fld[ 8-1]= 6.60365814  * 43 / 44.41;
        fld[ 9-1]= 7.199184884 * 43 / 44.41;
        fld[10-1]= 7.496866279 * 43 / 44.41;
        fld[11-1]= 7.801623256 * 43 / 44.41;
        fld[12-1]= 8.09494186  * 43 / 44.41;
        fld[13-1]= 8.394638372 * 43 / 44.41;
        fld[14-1]= 8.985097674 * 43 / 44.41;
*/
	cur[0] = 0.00000;
	cur[1] = 10.00000;
	cur[2] = 20.00000;
	cur[3] = 30.00000;
	cur[4] = 40.00000;
	cur[5] = 50.00000;
	cur[6] = 60.00000;
	cur[7] = 70.00000;
	cur[8] = 80.00000;
	cur[9] = 90.00000;
	cur[10] = 100.00000;
	cur[11] = 110.00000;
	cur[12] = 120.00000;
	cur[13] = 130.00000;
	cur[14] = 140.00000;
	cur[15] = 150.00000;
	cur[16] = 160.00000;
	cur[17] = 170.00000;
	cur[18] = 180.00000;
	cur[19] = 190.00000;
	cur[20] = 200.00000;

	fld[0] = 0.00000;
	fld[1] = 0.68278;
	fld[2] = 1.36555;
	fld[3] = 2.04833;
	fld[4] = 2.73111;
	fld[5] = 3.41388;
	fld[6] = 4.09666;
	fld[7] = 4.77943;
	fld[8] = 5.46221;
	fld[9] = 6.14499;
	fld[10] = 6.82776;
	fld[11] = 7.51054;
	fld[12] = 8.19332;
	fld[13] = 8.87609;
	fld[14] = 9.55887;
	fld[15] = 10.24164;
	fld[16] = 10.92442;
	fld[17] = 11.60720;
	fld[18] = 12.28997;
	fld[19] = 12.97275;
	fld[20] = 13.65553;


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

