//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_b_melco
//F	1		( mode, energy, kvalue, current, efflen, field )
//F!	----------------------------------------------------------------------
//F!
//F!	(parameter)
//F!	    coil:	i4, main coil = 1, trim coil = 2.
//F!	    mode:	i4, =1 for b to i, =2 for i to b.
//F!	    energy:	r4, energy in gev.
//F!	    efflen:	r4, effective length in meter.
//F!	    kvalue:	r4, k value.
//F!	    field:	r4, field data in t/m.
//F!	    current:	r4, magnet current in ampare.
//F!
//F!	(field data)
//F!	    measured by mitsubishi electoric co. ltd. (melco).
//F!	    only one magnet measured in this type.
//F!	    measurement was done at r=5730mm. 
//F!
//F!	(effective length)
//F!	    measurement was done at r=5730mm. i = 1000.2(a). 
//F!	    result from melco shows that the effective length is 1022.70mm.
//F!
//F!	(specifications)
//F!	    main coils: 1000a, 24turns.
//F!	    trim coils: 5a, 40turns.
//F!
//F!	(history)
//F!	    09-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include "magnet.h"

int mag_ki_b_melco( int mode, float energy,
		float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_calc
	int status;

//F	parameter   mx = 17
	const int mx=17;

//F	integer*4   mode_k_to_i / 1 /, mode_i_to_k / 2 /
	const int mode_k_to_i = 1;
	const int mode_i_to_k = 2;

//F	real*4	    cur(mx), fld(mx), energy, efflen
	float cur[mx], fld[mx];
//F	real*4	    kvalue(2), field(2), current(2)

//F	real*4	    cur_main, kvl_main, fld_main
	float cur_main, kvl_main, fld_main;

//F	real*4	    cur_trim, kvl_trim, fld_trim
	float cur_trim, kvl_trim, fld_trim;

//F	real*4	    cur_total, kvl_total, fld_total
	float cur_total, kvl_total, fld_total;

//F	data cur(01) /    0.000 /, fld(01) / 0.000000 /
//F	data cur(02) /    0.240 /, fld(02) / 0.000571 /
//F	data cur(03) /   99.807 /, fld(03) / 0.092213 /
//F	data cur(04) /  201.411 /, fld(04) / 0.186149 /
//F	data cur(05) /  299.124 /, fld(05) / 0.276556 /
//F	data cur(06) /  413.427 /, fld(06) / 0.382121 /
//F	data cur(07) /  507.789 /, fld(07) / 0.468801 /
//F	data cur(08) /  599.295 /, fld(08) / 0.552041 /
//F	data cur(09) /  650.544 /, fld(09) / 0.598576 /
//F	data cur(10) /  702.777 /, fld(10) / 0.646138 /
//F	data cur(11) /  749.325 /, fld(11) / 0.688338 /
//F	data cur(12) /  799.914 /, fld(12) / 0.734115 /
//F	data cur(13) /  855.612 /, fld(13) / 0.784023 /
//F	data cur(14) /  901.839 /, fld(14) / 0.825531 /
//F	data cur(15) /  956.214 /, fld(15) / 0.874221 /
//F	data cur(16) /  999.819 /, fld(16) / 0.912321 /
//F	data cur(17) / 1052.862 /, fld(17) / 0.958653 /

	cur[ 1-1] =    0.000; fld[ 1-1] = 0.000000; //
	cur[ 2-1] =    0.240; fld[ 2-1] = 0.000571; //
	cur[ 3-1] =   99.807; fld[ 3-1] = 0.092213; //
	cur[ 4-1] =  201.411; fld[ 4-1] = 0.186149; //
	cur[ 5-1] =  299.124; fld[ 5-1] = 0.276556; //
	cur[ 6-1] =  413.427; fld[ 6-1] = 0.382121; //
	cur[ 7-1] =  507.789; fld[ 7-1] = 0.468801; //
	cur[ 8-1] =  599.295; fld[ 8-1] = 0.552041; //
	cur[ 9-1] =  650.544; fld[ 9-1] = 0.598576; //
	cur[10-1] =  702.777; fld[10-1] = 0.646138; //
	cur[11-1] =  749.325; fld[11-1] = 0.688338; //
	cur[12-1] =  799.914; fld[12-1] = 0.734115; //
	cur[13-1] =  855.612; fld[13-1] = 0.784023; //
	cur[14-1] =  901.839; fld[14-1] = 0.825531; //
	cur[15-1] =  956.214; fld[15-1] = 0.874221; //
	cur[16-1] =  999.819; fld[16-1] = 0.912321; //
	cur[17-1] = 1052.862; fld[17-1] = 0.958653; //

//F	parameter ncoil_main = 24
	const int ncoil_main = 24;

//F	parameter ncoil_trim = 40
	const int ncoil_trim = 40;

//F	parameter eff_length = 1.02270		! meter
	const float eff_length = 1.02270;

//F	efflen = eff_length
	*efflen = eff_length;

//F	kvl_main = kvalue(1)
	kvl_main = kvalue[1-1];

//F	kvl_trim = kvalue(2)
	kvl_trim = kvalue[2-1];

//F	cur_main = current(1)
	cur_main = current[1-1];

//F	cur_trim = current(2)
	cur_trim = current[2-1];

//F! main coil control

//F1000	mag_ki_b_melco = mag_ki_calc(
//F	1		    mode, mx, cur, fld, energy, efflen, 
//F	1		    kvl_main, fld_main, cur_main )
	status = mag_ki_calc(
			mode, mx, cur, fld, energy, *efflen,
			&kvl_main, &fld_main, &cur_main );

//F! trim coil control

//F2000	if( mode .eq. mode_k_to_i ) then
	if( mode == mode_k_to_i ){

//F	    kvl_total = kvl_main + kvl_trim
	    kvl_total = kvl_main + kvl_trim;

//F	    mag_ki_b_melco = mag_ki_calc( 
//F	1			mode, mx, cur, fld, energy, efflen, 
//F	1			kvl_total, fld_total, cur_total )
	    status = mag_ki_calc( mode, mx, cur, fld, energy, *efflen,
				&kvl_total, &fld_total, &cur_total );

//F	    fld_trim = fld_total - fld_main
	    fld_trim = fld_total - fld_main;

//F	    cur_trim = cur_total - cur_main
	    cur_trim = cur_total - cur_main;

//F	    cur_trim = real(ncoil_main)/real(ncoil_trim)*cur_trim
	    cur_trim = (float)ncoil_main / (float)ncoil_trim * cur_trim;

//F	    current(1) = cur_main
	    current[1-1] = cur_main;

//F	    current(2) = cur_trim
	    current[2-1] = cur_trim;

//F	    field  (1) = fld_main
	    field[1-1] = fld_main;

//F	    field  (2) = fld_trim
	    field[2-1] = fld_trim;

//F	elseif( mode .eq. mode_i_to_k ) then
	} else if( mode == mode_i_to_k ){

//F	    cur_total = cur_main + real(ncoil_trim)/real(ncoil_main)*cur_trim
	    cur_total = cur_main + (float)ncoil_trim / (float)ncoil_main * cur_trim;

//F	    mag_ki_b_melco = mag_ki_calc( 
//F	1			mode, mx, cur, fld, energy, efflen, 
//F	1			kvl_total, fld_total, cur_total )
	    status = mag_ki_calc( mode, mx, cur, fld, energy, *efflen,
				&kvl_total, &fld_total, &cur_total );

//F	    kvl_trim = kvl_total - kvl_main
	    kvl_trim = kvl_total - kvl_main;

//F	    fld_trim = fld_total - fld_main
	    fld_trim = fld_total - fld_main;

//F	    kvalue(1) = kvl_main
	    kvalue[1-1] = kvl_main;

//F	    kvalue(2) = kvl_trim
	    kvalue[2-1] = kvl_trim;

//F	    field (1) = fld_main
	    field[1-1] = fld_main;

//F	    field (2) = fld_trim
	    field[2-1] = fld_trim;

//F	end if
	}

	return( status );
//F	end
}
