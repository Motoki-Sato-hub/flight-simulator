//F!==============================================================================
//F!
//F!
//F	integer*4 function mag_ki_control
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
//F!
//F!
//F!	supplied parameter
//F!		mode,	    ! access mode, =1 for k_to_i, =2 for i_to_k
//F!		mx,	    ! number of data for excitation.
//F!		cur,	    ! current data array
//F!		fld,	    ! field data array
//F!		energy,	    ! energy in gev
//F!		efflen,	    ! effective length (m)
//F!		ncoil_main, ! number of turn for main coil.
//F!		ncoil_trim, ! number of turn for trim coil.
//F!
//F!	supplied or returned parameter
//F!		kvalue,	    ! array of k values.
//F!		field,	    ! array of field
//F!		current	    ! array of current
//F!
//F!	history
//F!	    13-apr-1998 n.terunuma, created to use main and trim coils.
//F!
//F!==============================================================================
#include <stdio.h>
#include "magnet.h"

int mag_ki_control(
	int     mode,	        // Access mode
	int     mx,	        // Number of data for excitation.
	float    cur[],	        // Current data array
	float    fld[],	        // Field data array
	float    energy,	// Energy in gev
	float    eff_len,	        // Effective length (m)
	float    ncoil_main,   // Number of turn for main coil.
	float    ncoil_trim,    // Number of turn for trim coil.
	float    kvalue[2],     // Array of k values.
	float    field[2],        // Array of field
	float    current[2]    // Array of current
){
	const int mode_k_to_i = 1;
	const int mode_i_to_k = 2;

	int ncoil[2];
	int status;
	float flddata[2], dc[2];
	float cur_main, cur_trim;

	ncoil[0] = ncoil_main;
	ncoil[1] = ncoil_trim;

	if( mode == mode_k_to_i ){

		flddata[0] = energy * fabsf(kvalue[0]) / eff_len / 0.3;
//		flddata[1] = energy * fabsf(*kvalue[1]) / eff_len / 0.3;
		flddata[1] = energy *      (kvalue[1]) / eff_len / 0.3;

		status = mag_fld_convert( mode, mx, cur, fld, ncoil, flddata, dc );
		if( status ){
			current[0] = dc[0];
			current[1] = dc[1];
			field[0]   = flddata[0];
			field[1]   = flddata[1];
		}

	}else if ( mode == mode_i_to_k ){

		dc[0] = cur_main = current[0];
		dc[1] = cur_trim = current[1];

		status = mag_fld_convert( mode, mx, cur, fld, ncoil, flddata, dc );
		if( status ){
			field[0]  = flddata[0];
			field[1]  = flddata[1];
			kvalue[0] = flddata[0] / energy * eff_len * 0.3;
			kvalue[1] = flddata[1] / energy * eff_len * 0.3;
		}
	}

	return( status );
}
