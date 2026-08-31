//F! ==============================================================================
//F!
//F!	----------------------------------------------------------------------
//F	integer*4 function mag_ki_s_slac_213
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
//F!	    measured by slac and reported by mark woodley.
//F!	    for type "2.13s3.00" (sx1x and sx2x)
//F!
//F!	(effective length)
//F!
//F!	(specifications)
//F!
//F!	(history)
//F!	    07-mar-2005	n.terunuma, created.
//F!
//F!==============================================================================
//F	implicit none
#include <stdio.h>
#include <math.h>
#include "magnet.h"

int mag_ki_s_slac_213( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] )
{

//F	integer*4   mode, mag_ki_control, i

//F	parameter mx = 21
	const int mx = 21;

//F        real*4      cur(mx), fld(mx), energy, efflen
	static float cur[21], fld[21];
//F        real*4      kvalue(2), field(2), current(2)
//F	logical	    firsttime / .true. /
	static char firsttime =1;
//F! 
//F!"2.13s3.00" (sd1x and sd2x)

//F! parameter mx = 21
//F	parameter ncoil_main = 87
	const int ncoil_main = 87;
//F	parameter ncoil_trim = 0
	const int ncoil_trim = 0;
//F	parameter eff_length = 0.1
	const float eff_length = 0.1;
//F	parameter r0 = 0.0254*(2.13/2.0)
	const float r0 = 0.0254*(2.13/2.0);

//F	data cur(01) /  0.0 /, fld(01) / 0.000000 /
//F	data cur(02) /  0.5 /, fld(02) / 0.004694 /
//F	data cur(03) /  1.0 /, fld(03) / 0.010154 /
//F	data cur(04) /  1.5 /, fld(04) / 0.015612 /
//F	data cur(05) /  2.0 /, fld(05) / 0.021068 /
//F	data cur(06) /  2.5 /, fld(06) / 0.026522 /
//F	data cur(07) /  3.0 /, fld(07) / 0.031973 /
//F	data cur(08) /  3.5 /, fld(08) / 0.037422 /
//F	data cur(09) /  4.0 /, fld(09) / 0.042868 /
//F	data cur(10) /  4.5 /, fld(10) / 0.048312 /
//F	data cur(11) /  5.0 /, fld(11) / 0.053753 /
//F	data cur(12) /  5.5 /, fld(12) / 0.059190 /
//F	data cur(13) /  6.0 /, fld(13) / 0.064624 /
//F	data cur(14) /  6.5 /, fld(14) / 0.070054 /
//F	data cur(15) /  7.0 /, fld(15) / 0.075480 /
//F	data cur(16) /  7.5 /, fld(16) / 0.080900 /
//F	data cur(17) /  8.0 /, fld(17) / 0.086314 /
//F	data cur(18) /  8.5 /, fld(18) / 0.091722 /
//F	data cur(19) /  9.0 /, fld(19) / 0.097121 /
//F	data cur(20) /  9.5 /, fld(20) / 0.102511 /
//F	data cur(21) / 10.0 /, fld(21) / 0.107890 /

//F	if( firsttime ) then
//F	    firsttime = .false.
//F	    do i=1, mx
//Fc		fld(i)=2.*fld(i)/(6e-3)**2
//F		fld(i)=2.*fld(i)/r0**2
//F	    end do
//F	end if

	if( firsttime ){
		cur[ 1-1] =  0.0; fld[ 1-1] = 0.000000; 
		cur[ 2-1] =  0.5; fld[ 2-1] = 0.004694; 
		cur[ 3-1] =  1.0; fld[ 3-1] = 0.010154; 
		cur[ 4-1] =  1.5; fld[ 4-1] = 0.015612; 
		cur[ 5-1] =  2.0; fld[ 5-1] = 0.021068; 
		cur[ 6-1] =  2.5; fld[ 6-1] = 0.026522; 
		cur[ 7-1] =  3.0; fld[ 7-1] = 0.031973; 
		cur[ 8-1] =  3.5; fld[ 8-1] = 0.037422; 
		cur[ 9-1] =  4.0; fld[ 9-1] = 0.042868; 
		cur[10-1] =  4.5; fld[10-1] = 0.048312; 
		cur[11-1] =  5.0; fld[11-1] = 0.053753; 
		cur[12-1] =  5.5; fld[12-1] = 0.059190; 
		cur[13-1] =  6.0; fld[13-1] = 0.064624; 
		cur[14-1] =  6.5; fld[14-1] = 0.070054; 
		cur[15-1] =  7.0; fld[15-1] = 0.075480; 
		cur[16-1] =  7.5; fld[16-1] = 0.080900; 
		cur[17-1] =  8.0; fld[17-1] = 0.086314; 
		cur[18-1] =  8.5; fld[18-1] = 0.091722; 
		cur[19-1] =  9.0; fld[19-1] = 0.097121; 
		cur[20-1] =  9.5; fld[20-1] = 0.102511; 
		cur[21-1] = 10.0; fld[21-1] = 0.107890; 

		int i;
		for( i=0; i<mx; i++ ){
			fld[i] = 2.*fld[i]/pow(r0,2);
		}

		firsttime = 0;
	}

//F	
//F        efflen = eff_length 
	*efflen = eff_length;

//F	mag_ki_s_slac_213 = mag_ki_control
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

