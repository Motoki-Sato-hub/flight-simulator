//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_tokin_3581
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!
//F!
//F!	field data for quadrupole magnet of tokin 6 cm. type number is 3581.
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
//F!	    manufacture	    : tokin co. ltd,.
//F!	    type number	    : 3581
//F!	    pole length	    :  60 mm
//F!	    bore diameter   :  42 mm
//F!
//F!	    main coil data
//F!		measured by tokin co. ltd,. by a hole probe.
//F!		all magnets were measured in this type. we use no.1 magnet data.
//F!		measurement was done at a offset of x=10mm. 
//F!		no numerical data exist, so we read data from b-i graph.
//F!		most part of data have a linearity but higher current part
//F!		has saturation. we use multi-point data to control the magnets.
//F!			      b(t/m) = liner-interpolation (t/m/a) * i(a).
//F!	    trim coil data
//F!		no numerical data exist, so we read data from b-i graph.
//F!		taken at main coil 245.0a. linearity of data was not so good. 
//F!		this may be caused by a saturation of the magnet.
//F!		data shows the field difference was 20g/5a.
//F!		it becomes 20g @10mm --> 0.2 t/m (at 5a),
//F!				b(t/m) = 0.04(t/m/a) * i(a)
//F!	    effective length
//F!		no numerical data exist, so we read data from b-z graph.
//F!		integrate the area of graph and divide it by peak value.
//F!		it shows the effective length is ??? meter.

//F!	parameter
//F!	    coils:  based on the drawings #3581-s100-00, 98-10-03.
//F!
//F!	history
//F!	    13-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_tokin_3581( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 11
	const int mx = 11;

//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F        parameter ncoil_main = 26
	const int ncoil_main = 26;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//F        parameter eff_length = 0.084339
	const float eff_length = 0.084339;

//F	data cur(01) /   0.0 /, fld(01) /  0.0 /
//F	data cur(02) /  50.0 /, fld(02) /  7.5 /
//F	data cur(03) /  80.0 /, fld(03) / 11.8 /
//F	data cur(04) / 100.0 /, fld(04) / 14.8 /
//F	data cur(05) / 130.0 /, fld(05) / 19.1 /
//F	data cur(06) / 150.0 /, fld(06) / 21.8 /
//F	data cur(07) / 170.0 /, fld(07) / 24.4 /
//F	data cur(08) / 190.0 /, fld(08) / 26.4 /
//F	data cur(09) / 210.0 /, fld(09) / 28.0 /
//F	data cur(10) / 230.0 /, fld(10) / 29.3 /
//F	data cur(11) / 245.0 /, fld(11) / 30.1 /

	cur[ 1-1] =   0.0; fld[ 1-1] =  0.0; 
	cur[ 2-1] =  50.0; fld[ 2-1] =  7.5; 
	cur[ 3-1] =  80.0; fld[ 3-1] = 11.8; 
	cur[ 4-1] = 100.0; fld[ 4-1] = 14.8; 
	cur[ 5-1] = 130.0; fld[ 5-1] = 19.1; 
	cur[ 6-1] = 150.0; fld[ 6-1] = 21.8; 
	cur[ 7-1] = 170.0; fld[ 7-1] = 24.4; 
	cur[ 8-1] = 190.0; fld[ 8-1] = 26.4; 
	cur[ 9-1] = 210.0; fld[ 9-1] = 28.0; 
	cur[10-1] = 230.0; fld[10-1] = 29.3; 
	cur[11-1] = 245.0; fld[11-1] = 30.1; 

//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_tokin_3581 = mag_ki_control
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

