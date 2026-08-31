//F!==============================================================================
//F	integer*4 function mag_ki_q_ecube_skew( mode, kvalue, current )
//F!
//F!	supplied
//F!	    mode    : 1 for k_to_i mode or
//F!		      2 for i_to_k mode.
//F!
//F!	supplied & returned
//F!	    kvalue  : k value
//F!	    current : current
//F!
//F!	2001.05.09 created by f.takagi
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_ecube_skew( int mode, float *kvalue, float *current )
{

//F	parameter mode_k_to_i=1
	const int mode_k_to_i=1;
//F	parameter mode_i_to_k=2
	const int mode_i_to_k=2;

//F	integer*4 mode
//F	real*4    kvalue, current

//F	if( mode .eq. mode_k_to_i ) then	    ! k-->i transfer mode
	if( mode == mode_k_to_i ){
//F	    current = 4.0*kvalue/0.003
	    *current = 4.0*(*kvalue)/0.003;
//F	else if( mode .eq. mode_i_to_k ) then	    ! i-->k transfer mode
	}else if( mode == mode_i_to_k ){
//F	    kvalue  = 0.003*current/4.0
	    *kvalue  = 0.003*(*current)/4.0;
//F	else
	}else{
//F	    return
	    return(0);
//F	end if
	}

//F	mag_ki_q_ecube_skew = 1
	return(1);
//F	end
}

