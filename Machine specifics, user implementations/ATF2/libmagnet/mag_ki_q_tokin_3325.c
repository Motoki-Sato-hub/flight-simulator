//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_tokin_3325
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!
//F!	field data for quadrupole magnet of tokin 18 cm. type number is 3325.
//F!
//F!	(parameter)
//F!           mode:       i4, = 1 mode_k_to_i
//F!                           = 2 mode_i_to_k
//F!                           = 3 mode_b_to_i
//F!                           = 4 mode_i_to_b
//F!	    energy:	r4, energy in gev.
//F!	    kvalue:	r4, k value.
//F!	    current:	r4, magnet current in ampare.
//F!	    efflen:	r4, effective length in meter.
//F!	    field:	r4, field data in t/m.
//F!
//F!	(magnet specifications)
//F!	    manufacture	    : tokin co. ltd,.
//F!	    type number	    : 3325
//F!	    pole length	    : 180 mm
//F!	    bore diameter   :  32 mm
//F!
//F!	    main coil data
//F!		measured by tokin co. ltd,. by a hole probe.
//F!		all magnets were measured in this type. we use no.1 magnet data.
//F!		measurement was done at a offset of x=20mm. 
//F!		no numerical data exist, so we read data from b-i graph.
//F!		most part of data have a linearity but higher current part
//F!		has saturation. we use multi-point data to control the magnets.
//F!			      b(t/m) = liner-interpolation (t/m/a) * i(a).
//F!	    trim coil data
//F!		no numerical data exist, so we read data from b-i graph.
//F!		taken at main coil 510.0a. linearity of data was good. 
//F!		data shows the field difference was 200g/5a.
//F!		it becomes 200g @20mm --> 1 t/m (at 5a),
//F!				b(t/m) = 0.2(t/m/a) * i(a)
//F!	    effective length
//F!		no numerical data exist, so we read data from b-z graph.
//F!		integrate the area of graph and divide it by peak value.
//F!		it shows the effective length is 0.21563 meter.
//F!		but, this value seems to be wrong. we set the effective length
//F!		as same as that of an average for hitachi 180mm thickness 
//F!		magnets, 0.19886 meter.
//F!
//F!	parameter
//F!
//F!	history
//F!	    13-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_tokin_3325( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 11;
//F        real*4      cur(mx), fld(mx), energy, efflen
//F        real*4      kvalue(2), field(2), current(2)
	float cur[mx], fld[mx];

//F        parameter ncoil_main = 11
	const int ncoil_main = 11;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.19886 ! temporary 
	const float eff_length = 0.19886;
	

//F	data cur(01) /    0.00 /, fld(01) /  0.00 /
//F	data cur(02) /  100.00 /, fld(02) / 10.80 /
//F	data cur(03) /  200.00 /, fld(03) / 21.45 /
//F	data cur(04) /  300.00 /, fld(04) / 32.05 /
//F	data cur(05) /  400.00 /, fld(05) / 42.65 /
//F	data cur(06) /  417.24 /, fld(06) / 44.70 /
//F	data cur(07) /  434.48 /, fld(07) / 46.45 /
//F	data cur(08) /  451.72 /, fld(08) / 48.10 /
//F	data cur(09) /  468.97 /, fld(09) / 49.70 /
//F	data cur(10) /  485.21 /, fld(10) / 51.05 /
//F	data cur(10) /  503.45 /, fld(10) / 52.35 /
//F	data cur(11) /  512.07 /, fld(11) / 52.80 /

	cur[ 1-1] =    0.00; fld[ 1-1] =  0.00; 
	cur[ 2-1] =  100.00; fld[ 2-1] = 10.80; 
	cur[ 3-1] =  200.00; fld[ 3-1] = 21.45; 
	cur[ 4-1] =  300.00; fld[ 4-1] = 32.05; 
	cur[ 5-1] =  400.00; fld[ 5-1] = 42.65; 
	cur[ 6-1] =  417.24; fld[ 6-1] = 44.70; 
	cur[ 7-1] =  434.48; fld[ 7-1] = 46.45; 
	cur[ 8-1] =  451.72; fld[ 8-1] = 48.10; 
	cur[ 9-1] =  468.97; fld[ 9-1] = 49.70; 
	cur[10-1] =  485.21; fld[10-1] = 51.05; 
	cur[10-1] =  503.45; fld[10-1] = 52.35; 
	cur[11-1] =  512.07; fld[11-1] = 52.80; 


//F	
//F        efflen = eff_length 
	*efflen = eff_length; 

//F	mag_ki_q_tokin_3325 = mag_ki_control
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

