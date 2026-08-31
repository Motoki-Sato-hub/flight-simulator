#include <stdio.h>
#include "magnet.h"

int mag_ki_q_sd4( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	parameter   mx = 11
	const int mx = 17;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 1;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867
	const float eff_length = 0.10078;

	cur[1-1]= 0.0;
	cur[2-1]= 0.397;
	cur[3-1]= 0.83;
	cur[4-1]= 1.399;
	cur[5-1]= 1.832;
	cur[6-1]= 2.302;
	cur[7-1]= 4.971;
	cur[8-1]= 7.443;
	cur[9-1]= 9.911;
	cur[10-1]= 15.068;
	cur[11-1]= 19.962;
	cur[12-1]= 24.962;
	cur[13-1]= 30.055;
	cur[14-1]= 35.066;
	cur[15-1]= 39.91;
	cur[16-1]= 45.061;
	cur[17-1]= 50.002;
	
	fld[ 1-1]= 0.0;
	fld[ 2-1]= 20.158;
	fld[ 3-1]= 29.394;
	fld[ 4-1]= 41.579;
	fld[ 5-1]= 50.795;
	fld[ 6-1]= 60.921;
	fld[ 7-1]= 119.34;
	fld[ 8-1]= 173.74;
	fld[ 9-1]= 228.64;
	fld[10-1]= 344.14;
	fld[11-1]= 454.79;
	fld[12-1]= 568.03;
	fld[13-1]= 682.74;
	fld[14-1]= 794.88;
	fld[15-1]= 902.77;
	fld[16-1]= 1017.2;
	fld[17-1]= 1126.3;

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

