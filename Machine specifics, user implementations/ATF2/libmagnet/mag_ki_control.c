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
	float    efflen,	        // Effective length (m)
	float    ncoil_main,   // Number of turn for main coil.
	float    ncoil_trim,    // Number of turn for trim coil.
	float    kvalue[2],     // Array of k values.
	float    field[2],        // Array of field
	float    current[2]    // Array of current
){
//F	implicit none
//F	integer*4   mode, mx, ncoil_main, ncoil_trim, mag_ki_calc
//F        integer*4   mode_k_to_i / 1 /, mode_i_to_k / 2 /
	const int mode_k_to_i = 1;
	const int mode_i_to_k = 2;

//F        real*4      cur(100), fld(100), energy, efflen
//F        real*4      kvalue(2), field(2), current(2)

//F        real*4      cur_main,  kvl_main,  fld_main
	float cur_main,  kvl_main,  fld_main;
//F        real*4      cur_trim,  kvl_trim,  fld_trim
	float cur_trim,  kvl_trim,  fld_trim;
//F        real*4      cur_total, kvl_total, fld_total
	float cur_total, kvl_total, fld_total;

	int status;
//F!--------------------------------------------------------------------------
//F    kvl_main = kvalue(1)
	kvl_main = kvalue[0];
//F    kvl_trim = kvalue(2)
	kvl_trim = kvalue[1];
//F    cur_main = current(1)
	cur_main = current[0];
//F    cur_trim = current(2)
	cur_trim = current[1];
//F! main coil control

//F	mag_ki_control = mag_ki_calc( mode, mx, cur, fld, energy, efflen,
//F	1			 kvl_main, fld_main, cur_main )
//F				 if( .not. mag_ki_control ) goto 9000
	status = mag_ki_calc( mode, mx, cur, fld, energy, efflen,
			&kvl_main, &fld_main, &cur_main );
	if( status%2 != 1 ) return(status);

//F! trim coil control

//F	if( mode .eq. mode_k_to_i ) then
	if( mode == mode_k_to_i ){

//F	    kvl_total = kvl_main + kvl_trim
	    kvl_total = kvl_main + kvl_trim;

//F	    mag_ki_control = mag_ki_calc(
//F	1			    mode, mx, cur, fld, energy, efflen,
//F	1			    kvl_total, fld_total, cur_total )
//F				    if( .not. mag_ki_control ) goto 9000
	    status = mag_ki_calc( mode, mx, cur, fld, energy, efflen,
				&kvl_total, &fld_total, &cur_total );
	    if( status%2 != 1 ) return(status);

//F	    fld_trim = fld_total - fld_main
	    fld_trim = fld_total - fld_main;
//F	    cur_trim = cur_total - cur_main
	    cur_trim = cur_total - cur_main;

//Fc 99/12/02 j.ozawa
//Fc	    cur_trim = real(ncoil_main)/real(ncoil_trim)*cur_trim

//F	    if( ncoil_trim .ne. 0 ) then
	    if( ncoil_trim != 0 ){
//F		cur_trim = real(ncoil_main)/real(ncoil_trim)*cur_trim
		cur_trim = (float)ncoil_main / (float)ncoil_trim * cur_trim;
//F	    else
	    }else{
//F		cur_trim = 0.0
		cur_trim = 0.0;
//F	    endif
	    }
//Fcc

//F	    current(1) = cur_main
	    current[0] = cur_main;
//F	    current(2) = cur_trim
	    current[1] = cur_trim;
//F	    field  (1) = fld_main
	    field[0] = fld_main;
//F	    field  (2) = fld_trim
	    field[1] = fld_trim;


//F	elseif( mode .eq. mode_i_to_k ) then
	}else if ( mode == mode_i_to_k ){

//F	    cur_total = cur_main + real(ncoil_trim)/real(ncoil_main)*cur_trim
	    cur_total = cur_main + (float)ncoil_trim / (float)ncoil_main * cur_trim;

//F	    mag_ki_control = mag_ki_calc( 
//F	1			    mode, mx, cur, fld, energy, efflen,
//F	1			    kvl_total, fld_total, cur_total )
//F				    if( .not. mag_ki_control ) goto 9000
	    status = mag_ki_calc( mode, mx, cur, fld, energy, efflen,
				&kvl_total, &fld_total, &cur_total );
	    //printf("mag_ki_calc status = %d\n", status );
	    if( status%2 != 1 ) return(status);
		
//F	    kvl_trim = kvl_total - kvl_main
	    kvl_trim = kvl_total - kvl_main;
//F	    fld_trim = fld_total - fld_main
	    fld_trim = fld_total - fld_main;

//F	    kvalue(1) = kvl_main
	    kvalue[0] = kvl_main;
//F	    kvalue(2) = kvl_trim
	    kvalue[1] = kvl_trim;
//F	    field (1) = fld_main
	    field [0] = fld_main;
//F	    field (2) = fld_trim
	    field [1] = fld_trim;

//F	end if
	}

//F9000	continue
	return( status );
//F	end
}
