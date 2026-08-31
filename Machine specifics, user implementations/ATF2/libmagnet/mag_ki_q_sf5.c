#include <stdio.h>
#include "magnet.h"

int mag_ki_q_sf5( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	parameter   mx = 11
	//const int mx = 17;
	//const int mx = 11;
	const int mx=21;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 1;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867

	// 2015-08-31 Y.Tsukada(KIS)
	// SF5マグネットの入れ替えのため、データ変更
	//const float eff_length = 0.10078;
	const float eff_length = 0.07018;
/*
	cur[1-1]= 0.0;
	cur[2-1]= 0.396;
	cur[3-1]= 0.829;
	cur[4-1]= 1.399;
	cur[5-1]= 1.832;
	cur[6-1]= 2.302;
	cur[7-1]= 4.97;
	cur[8-1]= 7.442;
	cur[9-1]= 9.909;
	cur[10-1]= 15.065;
	cur[11-1]= 19.962;
	cur[12-1]= 24.961;
	cur[13-1]= 30.055;
	cur[14-1]= 35.066;
	cur[15-1]= 39.912;
	cur[16-1]= 45.063;
	cur[17-1]= 50.004;

		
	fld[ 1-1]= 0.0;
	fld[ 2-1]= 19.89;
	fld[ 3-1]= 29.094;
	fld[ 4-1]= 41.236;
	fld[ 5-1]= 50.425;
	fld[ 6-1]= 60.517;
	fld[ 7-1]= 118.87;
	fld[ 8-1]= 173.08;
	fld[ 9-1]= 227.71;
	fld[10-1]= 343.04;
	fld[11-1]= 453.28;
	fld[12-1]= 565.99;
	fld[13-1]= 680.01;
	fld[14-1]= 791.82;
	fld[15-1]= 899.54;
	fld[16-1]= 1013.7;
	fld[17-1]= 1122.6;
*/

	cur[1-1]  = -10.0;
	cur[2-1]  = -9.0;
	cur[3-1]  = -8.0;
	cur[4-1]  = -7.0;
	cur[5-1]  = -6.0;
	cur[6-1]  = -5.0;
	cur[7-1]  = -4.0;
	cur[8-1]  = -3.0;
	cur[9-1] = -2.0;
	cur[10-1] = -1.0;
	cur[11-1]  = 0.0;
	cur[12-1]  = 1.0;
	cur[13-1]  = 2.0;
	cur[14-1]  = 3.0;
	cur[15-1]  = 4.0;
	cur[16-1]  = 5.0;
	cur[17-1]  = 6.0;
	cur[18-1]  = 7.0;
	cur[19-1]  = 8.0;
	cur[20-1] = 9.0;
	cur[21-1] = 10.0;


	fld[1-1]  = -319.3;
	fld[2-1]  = -286.9;
	fld[3-1]  = -254.6;
	fld[4-1]  = -223.1;
	fld[5-1]  = -191.0;
	fld[6-1]  = -159.1;
	fld[7-1]  = -127.0;
	fld[8-1]  = -95.27;
	fld[9-1] = -64.29;
	fld[10-1] = -33.47;

	//fld[1-1]  = 3.38;
	fld[11-1]  = 0.0;
	fld[12-1]  = 33.47;
	fld[13-1]  = 64.29;
	fld[14-1]  = 95.27;
	fld[15-1]  = 127.0;
	fld[16-1]  = 159.1;
	fld[17-1]  = 191.0;
	fld[18-1]  = 223.1;
	fld[19-1]  = 254.6;
	fld[20-1] = 286.9;
	fld[21-1] = 319.3;

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

