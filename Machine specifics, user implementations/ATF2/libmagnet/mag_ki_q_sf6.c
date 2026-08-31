#include <stdio.h>
#include "magnet.h"

int mag_ki_q_sf6( int mode, float energy, 
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
	cur[2-1]= 0.398;
	cur[3-1]= 0.829;
	cur[4-1]= 1.399;
	cur[5-1]= 1.832;
	cur[6-1]= 2.301;
	cur[7-1]= 4.972;
	cur[8-1]= 7.443;
	cur[9-1]= 9.912;
	cur[10-1]= 15.069;
	cur[11-1]= 19.962;
	cur[12-1]= 24.962;
	cur[13-1]= 30.055;
	cur[14-1]= 35.066;
	cur[15-1]= 39.908;
	cur[16-1]= 45.059;
	cur[17-1]= 49.999;
	
	fld[ 1-1]= 0.0;
	fld[ 2-1]= 22.481;
	fld[ 3-1]= 31.634;
	fld[ 4-1]= 43.738;
	fld[ 5-1]= 52.886;
	fld[ 6-1]= 62.946;
	fld[ 7-1]= 121.29;
	fld[ 8-1]= 175.41;
	fld[ 9-1]= 230.08;
	fld[10-1]= 345.5;
	fld[11-1]= 456.04;
	fld[12-1]= 569.96;
	fld[13-1]= 685.74;
	fld[14-1]= 799.28;
	fld[15-1]= 908.37;
	fld[16-1]= 1023.8;
	fld[17-1]= 1133.9;

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

