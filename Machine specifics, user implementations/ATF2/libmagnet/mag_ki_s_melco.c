//F!==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_s_melco
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
//F!
//F!	(effective length)
//F!
//F!	(specifications)
//F!
//F!	(history)
//F!	    09-apr-1998	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include <math.h>
#include "magnet.h"

int mag_ki_s_melco( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{
//F	integer*4   mode, mag_ki_control, i

//F	parameter   mx = 17
	const int mx=17;
//F        real*4      cur(mx), fld(mx), energy, efflen
	static float cur[17], fld[17];
//F        real*4      kvalue(2), field(2), current(2)
//F	logical	    firsttime / .true. /
	static char firsttime = 1;
//F! 
//F! s_mitsubishi.dat
//F!
//F        parameter ncoil_main = 17
	const int ncoil_main = 17;
//F        parameter ncoil_trim = 20
	const int ncoil_trim = 20;
//Fc        parameter eff_length = 0.06678	  ! by fit
//F        parameter eff_length = 0.07077	  ! changed 21-apr-98
	const float eff_length = 0.07077;

//F	data cur(01) /    0.00000 /, fld(01) /    0.00000 /
//F	data cur(02) /    0.07500 /, fld(02) /    0.00050 /
//F	data cur(03) /   25.01000 /, fld(03) /    0.00986 /
//F	data cur(04) /   49.98700 /, fld(04) /    0.01959 /
//F	data cur(05) /   75.00100 /, fld(05) /    0.02927 /
//F	data cur(06) /  124.99800 /, fld(06) /    0.04816 /
//F	data cur(07) /  149.99100 /, fld(07) /    0.05737 /
//F	data cur(08) /  175.00000 /, fld(08) /    0.06623 /
//F	data cur(09) /  200.00700 /, fld(09) /    0.07398 /
//F	data cur(10) /  224.98399 /, fld(10) /    0.07947 /
//F	data cur(11) /  249.98300 /, fld(11) /    0.08339 /
//F	data cur(12) /  275.00299 /, fld(12) /    0.08645 /
//F	data cur(13) /  299.98700 /, fld(13) /    0.08891 /
//F	data cur(14) /  324.99500 /, fld(14) /    0.09093 /
//F	data cur(15) /  350.00000 /, fld(15) /    0.09261 /
//F	data cur(16) /  374.98999 /, fld(16) /    0.09404 /
//F	data cur(17) /  399.99600 /, fld(17) /    0.09527 /

//F	if( firsttime ) then
//F	    firsttime = .false.
//F	    do i=1, mx
//F		fld(i)=2.*fld(i)/(6e-3)**2
//F	    end do
//F	end if
	if( firsttime ){
		cur[ 1-1] =    0.00000; fld[ 1-1] =    0.00000; 
		cur[ 2-1] =    0.07500; fld[ 2-1] =    0.00050; 
		cur[ 3-1] =   25.01000; fld[ 3-1] =    0.00986; 
		cur[ 4-1] =   49.98700; fld[ 4-1] =    0.01959; 
		cur[ 5-1] =   75.00100; fld[ 5-1] =    0.02927; 
		cur[ 6-1] =  124.99800; fld[ 6-1] =    0.04816; 
		cur[ 7-1] =  149.99100; fld[ 7-1] =    0.05737; 
		cur[ 8-1] =  175.00000; fld[ 8-1] =    0.06623; 
		cur[ 9-1] =  200.00700; fld[ 9-1] =    0.07398; 
		cur[10-1] =  224.98399; fld[10-1] =    0.07947; 
		cur[11-1] =  249.98300; fld[11-1] =    0.08339; 
		cur[12-1] =  275.00299; fld[12-1] =    0.08645; 
		cur[13-1] =  299.98700; fld[13-1] =    0.08891; 
		cur[14-1] =  324.99500; fld[14-1] =    0.09093; 
		cur[15-1] =  350.00000; fld[15-1] =    0.09261; 
		cur[16-1] =  374.98999; fld[16-1] =    0.09404; 
		cur[17-1] =  399.99600; fld[17-1] =    0.09527; 

		int i;
		for( i=0; i<mx; i++ ){
			fld[i] = 2.*fld[i]/pow((6e-3),2);
		}
	}

//F	
//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_s_melco = mag_ki_control
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

