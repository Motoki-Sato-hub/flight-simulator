//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_tokin_3393
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	field data for quadrupole magnet of tokin 6cm. type number is 3393.
//F!
//F!	(parameter)
//F!           mode:       i4, = 1 mode_k_to_i
//F!                           = 2 mode_i_to_k
//F!                           = 3 mode_b_to_i
//F!                           = 4 mode_i_to_b
//F!	    energy:	r4, energy in gev.
//F!	    kvalue:	r4, k value.
//F!	    current:	r4, magnet current in ampare.
//F!	    field:	r4, field data in t/m.
//F!	    efflen:	r4, effective length in meter.
//F!
//F!	(magnet specifications)
//F!	    manufacture	    : tokin co. ltd,.
//F!	    type number	    : 3393
//F!	    pole length	    : 60 mm
//F!	    bore diameter   : 32 mm
//F!
//F!	    main coil data
//F!		measured by tokin co. ltd,. by a hole probe.
//F!		all magnets were measured in this type. we use no.1 magnet data.
//F!		measurement was done at a offset of x=20mm. 
//F!		no numerical data exist, so we read data from b-i graph.
//F!		current value of the maximum point was 139a and field was 4.6kg.
//F!		linerity of the data was very good.
//F!		so we set,    4.6kg @20mm --> 23.0 t/m (at 139a).
//F!			      b(t/m) = 0.16547(t/m/a) * i(a)
//F!	    trim coil data
//F!		no numerical data exist, so we read data from b-i graph.
//F!		taken at main coil 139.0a. linearity of data was good. 
//F!		data shows the field difference was 200g/5a.
//F!		it becomes 200g @20mm --> 1 t/m (at 5a),
//F!				b(t/m) = 0.2(t/m/a) * i(a)
//F!	    effective length
//F!		no numerical data exist, so we read data from b-z graph.
//F!		integrate the area of graph and divide it by peak value.
//F!		it shows the effective length is 0.079148 meter.
//F!
//F!	documents
//F!
//F!	history
//F!	    14-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_tokin_3390a( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 8
	const int mx = 17;
//F        real*4      cur(mx), fld(mx), energy, efflen
//F        real*4      kvalue(2), field(2), current(2)
	float cur[mx], fld[mx];

//F        parameter ncoil_main = 17
	const int ncoil_main = 67;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.07890677
	//const float eff_length = 0.07890677;
	const float eff_length = 0.19742664;

//F	data cur( 1) /   0.0 /, fld( 1) /  0.0 /
//F	data cur( 2) /  20.0 /, fld( 2) /  3.4 /
//F	data cur( 3) /  40.0 /, fld( 3) /  6.8 /
//F	data cur( 4) /  60.0 /, fld( 4) / 10.0 /
//F	data cur( 5) /  80.0 /, fld( 5) / 13.2 /
//F	data cur( 6) / 100.0 /, fld( 6) / 16.6 /
//F	data cur( 7) / 120.0 /, fld( 7) / 19.8 /
//F	data cur( 8) / 139.0 /, fld( 8) / 23.0 /

	cur[ 1-1] =   0.0; fld[ 1-1] =  0.0; 
	cur[ 2-1] =  14.44; fld[ 2-1] =  2.00*5.0; 
	cur[ 3-1] =  20.00; fld[ 3-1] =  2.70*5.0; 
	cur[ 4-1] =  29.62; fld[ 4-1] = 4.00*5.0; 
	cur[ 5-1] =  40.00; fld[ 5-1] = 5.32*5.0; 
	cur[ 6-1] = 45.18; fld[ 6-1] = 6.00*5.0; 
	cur[ 7-1] = 60.00; fld[ 7-1] = 7.93*5.0; 
	cur[ 8-1] = 65.92; fld[ 8-1] = 8.68*5.0; 
	cur[ 9-1] = 71.11; fld[ 9-1] = 9.28*5.0; 
	cur[ 10-1] = 75.92; fld[ 10-1] = 9.77*5.0; 
	cur[ 11-1] = 80.00; fld[ 11-1] = 10.10*5.0; 
	cur[ 12-1] = 84.07; fld[ 12-1] = 10.40*5.0; 
	cur[ 13-1] = 87.77; fld[ 13-1] = 10.60*5.0; 
	cur[ 14-1] = 91.11; fld[ 14-1] = 10.80*5.0; 
	cur[ 15-1] = 94.44; fld[ 15-1] = 11.00*5.0; 
	cur[ 16-1] = 98.15; fld[ 16-1] = 11.20*5.0; 
	cur[ 17-1] = 100.00; fld[ 17-1] = 11.30*5.0; 

//F	
//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_tokin_3393 = mag_ki_control
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

