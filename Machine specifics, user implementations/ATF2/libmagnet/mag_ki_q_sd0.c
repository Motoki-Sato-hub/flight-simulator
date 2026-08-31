#include <stdio.h>
#include "magnet.h"

int mag_ki_q_sd0( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	parameter   mx = 11
	const int mx = 5;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 39
	const int ncoil_main = 1;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 0;
//F        parameter eff_length = 0.07867
	const float eff_length = 0.1;

	cur[1-1]= 0.0;
	cur[2-1]= 2.9621;
	cur[3-1]= 5.9571;
	cur[4-1]= 8.9525;
	cur[5-1]= 11.964;
	
	fld[ 1-1]= 0.0;
	fld[ 2-1]= 85.983;
	fld[ 3-1]= 174.79;
	fld[ 4-1]= 263.39;
	fld[ 5-1]= 351.37;

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

