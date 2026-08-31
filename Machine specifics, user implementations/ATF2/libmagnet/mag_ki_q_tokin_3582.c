//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_q_tokin_3582
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	field data for quadrupole magnet of tokin 18 cm. type number is 3582.
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
//F!	    type number	    : 3582
//F!	    pole length	    : 180 mm
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
//F!		taken at main coil 245.0a. linearity of data was good. 
//F!		data shows the field difference was 50g/5a.
//F!		it becomes 50g @10mm --> 0.5 t/m (at 5a),
//F!				b(t/m) = 0.1(t/m/a) * i(a)
//F!	    effective length
//F!		no numerical data exist, so we read data from b-z graph.
//F!		integrate the area of graph and divide it by peak value.
//F!		it shows the effective length is 0.202628 meter.
//F!
//F!
//F!	documents
//F!	    drawings: #310q431-357, 96-01-11.
//F!
//F!	history
//F!	    13-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_q_tokin_3582( int mode, float energy, 
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
//F        parameter eff_length = 0.202628
	const float eff_length = 0.202628;
//F	
//F	data cur(01) /   0.0 /, fld(01) /  0.0 /
//F	data cur(02) /  50.0 /, fld(02) /  7.6 /
//F	data cur(03) / 100.0 /, fld(03) / 14.9 /
//F	data cur(04) / 150.0 /, fld(04) / 22.2 /
//F	data cur(05) / 170.0 /, fld(05) / 25.0 /
//F	data cur(06) / 190.0 /, fld(06) / 28.0 /
//F	data cur(07) / 210.0 /, fld(07) / 30.8 /
//F	data cur(08) / 220.0 /, fld(08) / 32.1 /
//F	data cur(09) / 230.0 /, fld(09) / 33.3 /
//F	data cur(10) / 240.0 /, fld(10) / 34.2 /
//F	data cur(11) / 245.0 /, fld(11) / 34.7 /

	cur[ 1-1] =   0.0; fld[ 1-1] =  0.0; 
	cur[ 2-1] =  50.0; fld[ 2-1] =  7.6; 
	cur[ 3-1] = 100.0; fld[ 3-1] = 14.9; 
	cur[ 4-1] = 150.0; fld[ 4-1] = 22.2; 
	cur[ 5-1] = 170.0; fld[ 5-1] = 25.0; 
	cur[ 6-1] = 190.0; fld[ 6-1] = 28.0; 
	cur[ 7-1] = 210.0; fld[ 7-1] = 30.8; 
	cur[ 8-1] = 220.0; fld[ 8-1] = 32.1; 
	cur[ 9-1] = 230.0; fld[ 9-1] = 33.3; 
	cur[10-1] = 240.0; fld[10-1] = 34.2; 
	cur[11-1] = 245.0; fld[11-1] = 34.7; 

//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_q_tokin_3582 = mag_ki_control
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

