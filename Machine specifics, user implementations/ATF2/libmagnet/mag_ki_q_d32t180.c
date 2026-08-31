//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_d32t180
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!
//F!
//F!	field data for quadrupole magnet of qea-d32t180 18 cm.
//F!	type number is 13.
//F!
//F!	(parameter)
//F!           mode:       i4, = 1 mode_k_to_i
//F!                           = 2 mode_i_to_k
//F!                           = 3 mode_b_to_i
//F!                           = 4 mode_i_to_b
//F!	    coil:	i4, main coil = 1, trim coil = 2.
//F!	    energy:	r4, energy in gev.
//F!	    kvalue:	r4, k value.
//F!	    current:	r4, magnet current in ampare.
//F!	    field:	r4, field data in t/m.
//F!	    efflen:	r4, effective length in meter.
//F!
//F!	(magnet specifications)
//F!	    manufacture	    : 
//F!	    type number	    :  13
//F!	    pole length	    : 180 mm
//F!	    bore diameter   :     mm
//F!
//F!	history
//F!	    16-oct-2006	created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_d32t180( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode
//F	parameter   mx = 14
	const int mx = 14;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F	integer*4   mag_ki_control_d32t180
//F	integer*4   mag_ki_calc_d32t180

//F!	parameter ncoil_main = 26
//F	parameter ncoil_main = 49	    ! 2007-03-07    change 
	const int ncoil_main = 49;
//F	parameter ncoil_trim = 20
	const int ncoil_trim = 20;

//F! 2006-oct-13 not used.
//F	parameter eff_length = 0
	const float eff_length = 0;
//F!!!

//F	data cur(01) /   0.0 /, fld(01) /  0.0 /
//F	data cur(02) /   5.0 /, fld(02) /  0.469959 /
//F	data cur(03) /   8.0 /, fld(03) /  0.748402 /
//F	data cur(04) /  10.0 /, fld(04) /  0.932288 /
//F	data cur(05) /  12.0 /, fld(05) /  1.119242 /
//F	data cur(06) /  14.0 /, fld(06) /  1.304396 /
//F	data cur(07) /  16.0 /, fld(07) /  1.490985 /
//F	data cur(08) /  20.0 /, fld(08) /  1.864267 /
//F	data cur(09) /  25.0 /, fld(09) /  2.330425 /
//F	data cur(10) /  30.0 /, fld(10) /  2.795264 /
//F	data cur(11) /  35.0 /, fld(11) /  3.253318 /
//F	data cur(12) /  40.0 /, fld(12) /  3.711283 /
//F	data cur(13) /  45.0 /, fld(13) /  4.168370 /
//F	data cur(14) /  50.0 /, fld(14) /  4.624610 /

	cur[ 1-1] =   0.0; fld[ 1-1] =  0.0; 
	cur[ 2-1] =   5.0; fld[ 2-1] =  0.469959; 
	cur[ 3-1] =   8.0; fld[ 3-1] =  0.748402; 
	cur[ 4-1] =  10.0; fld[ 4-1] =  0.932288; 
	cur[ 5-1] =  12.0; fld[ 5-1] =  1.119242; 
	cur[ 6-1] =  14.0; fld[ 6-1] =  1.304396; 
	cur[ 7-1] =  16.0; fld[ 7-1] =  1.490985; 
	cur[ 8-1] =  20.0; fld[ 8-1] =  1.864267; 
	cur[ 9-1] =  25.0; fld[ 9-1] =  2.330425; 
	cur[10-1] =  30.0; fld[10-1] =  2.795264; 
	cur[11-1] =  35.0; fld[11-1] =  3.253318; 
	cur[12-1] =  40.0; fld[12-1] =  3.711283; 
	cur[13-1] =  45.0; fld[13-1] =  4.168370; 
	cur[14-1] =  50.0; fld[14-1] =  4.624610; 

//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_d32t180 = mag_ki_control_d32t180
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
	return mag_ki_control_d32t180(
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

