//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_b_shi_c
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	(parameter)
//F!	    mode:	i4, =1 for b to i, =2 for i to b.
//F!	    energy:	r4, energy in gev.
//F!	    kvalue:	r4, k value.
//F!	    current:	r4, magnet current in ampare.
//F!	    efflen:	r4, effective length in meter.
//F!	    field:	r4, field data in t/m.
//F!
//F!	specifications
//F!	    pole length along the beam orbit, 1350 mm.
//F!	    main coils: 80 turns.
//F!	    trim coils: 20 turns.
//F!
//F!	field data
//F!	    there are no data prepared by sumitomo campany so we measure them
//F!	    by a hole probe. (terunuma,imai).
//F!
//F!	history
//F!	    15-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int  mag_ki_b_shi_c( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control

//F	parameter   mx = 17
	const int mx=17;
//F        real*4      cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F        real*4      kvalue(2), field(2), current(2)

//F	parameter ncoil_main = 80
	const int ncoil_main = 80;
//F	parameter ncoil_trim = 20
	const int ncoil_trim = 20;

//Fc 99/11/12 j.ozawa & okugi
//Fc	parameter eff_length = 1.37296		! meter
//F	parameter eff_length = 1.34296		! meter
	const float eff_length = 1.34296;

//F	data cur( 1) /   0.0 /, fld( 1) / 0.0000 /
//F	data cur( 2) /  50.0 /, fld( 2) / 0.1560 /
//F	data cur( 3) / 100.0 /, fld( 3) / 0.3095 /
//F	data cur( 4) / 150.0 /, fld( 4) / 0.4632 /
//F	data cur( 5) / 200.0 /, fld( 5) / 0.6145 /
//F	data cur( 6) / 250.0 /, fld( 6) / 0.7631 /
//F	data cur( 7) / 300.0 /, fld( 7) / 0.9092 /
//F	data cur( 8) / 320.0 /, fld( 8) / 0.9635 /
//F	data cur( 9) / 340.0 /, fld( 9) / 1.0110 /
//F	data cur(10) / 360.0 /, fld(10) / 1.0492 /
//F	data cur(11) / 380.0 /, fld(11) / 1.0806 /
//F	data cur(12) / 400.0 /, fld(12) / 1.1079 /
//F	data cur(13) / 420.0 /, fld(13) / 1.1313 /
//F	data cur(14) / 440.0 /, fld(14) / 1.1530 /
//F	data cur(15) / 460.0 /, fld(15) / 1.1727 /
//F	data cur(16) / 480.0 /, fld(16) / 1.1911 /
//F	data cur(17) / 500.0 /, fld(17) / 1.2000 /  ! add by interpolation

	cur[ 1-1] =   0.0; fld[ 1-1] = 0.0000; //
	cur[ 2-1] =  50.0; fld[ 2-1] = 0.1560; //
	cur[ 3-1] = 100.0; fld[ 3-1] = 0.3095; //
	cur[ 4-1] = 150.0; fld[ 4-1] = 0.4632; //
	cur[ 5-1] = 200.0; fld[ 5-1] = 0.6145; //
	cur[ 6-1] = 250.0; fld[ 6-1] = 0.7631; //
	cur[ 7-1] = 300.0; fld[ 7-1] = 0.9092; //
	cur[ 8-1] = 320.0; fld[ 8-1] = 0.9635; //
	cur[ 9-1] = 340.0; fld[ 9-1] = 1.0110; //
	cur[10-1] = 360.0; fld[10-1] = 1.0492; //
	cur[11-1] = 380.0; fld[11-1] = 1.0806; //
	cur[12-1] = 400.0; fld[12-1] = 1.1079; //
	cur[13-1] = 420.0; fld[13-1] = 1.1313; //
	cur[14-1] = 440.0; fld[14-1] = 1.1530; //
	cur[15-1] = 460.0; fld[15-1] = 1.1727; //
	cur[16-1] = 480.0; fld[16-1] = 1.1911; //
	cur[17-1] = 500.0; fld[17-1] = 1.2000; //  ! add by interpolation

//F	
//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_b_shi_c = mag_ki_control
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

