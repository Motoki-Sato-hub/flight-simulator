//F!==============================================================================
//F!
//F	integer*4 function mag_ki_calc
//F	1	( mode, n, cur, fld, energy, eff_len, kvalue, field, current )
//F!
//F!	make linear interpolation from the current-field data.
//F!	
//F!	(supplied)
//F!	    mode;	integer*4   ; b_to_i = 1 or i_to_b = 2.
//F!	    n;		integer*4   ; array size of cur and fld.
//F!	    cur;	real*4	    ; array of the current data.
//F!	    fld;	real*4	    ; array of the field data.
//F!	    energy;	real*4	    ; energy in gev.
//F!	    eff_len;	real*4	    ; effective length in meter.
//F!
//F!	(supplied and returned)
//F!	    kvalue;	real*4	    ; k value.
//F!	    current;	real*4	    ; current (a).
//F!
//F!	(possible return code)
//F!	    mag_ki_calc = 0
//F!	    mag_ki_calc = 1
//F!
//F!	author: n.terunuma (kek)
//F!	
//F!	history:    8-dec-1997 first version.
//F!	            9-Feb-2010 Checked for C. by N.Terunuma
//F!
//F!==============================================================================
#include <stdio.h>
#include <math.h>

int mag_ki_calc( int mode, int n, float cur[], float fld[], float energy, float eff_len, 
		 float *kvalue, float *field, float *current )
{
  int i;
  
  const int mode_b_to_i = 1;
  const int mode_i_to_b = 2;
  
  double xin;
  
  
  // Reject negative effective length
  
  if( eff_len <= 0.0 ) return(0);
  
  // Reject negative beam energy
  
  if( energy <= 0.0 ) return(0);
  
  if( mode == mode_b_to_i ){
    
    //------------------------------------------------------------------------
    
    *field = energy * fabsf(*kvalue) / (double)eff_len / (double)0.3;
    xin = *field;
    
    if( xin <= fld[0] ){
      *current = cur[0];
      return(1);
    }else if( xin >= fld[n-1] ){
      *current = cur[n-1];
      return(1);
    }
    
    for( i=0; i<(n-1); i++ ){
      if( (xin>fld[i]) && (xin<=fld[i+1]) ){
	*current = cur[i] + (cur[i+1]-cur[i]) / (fld[i+1]-fld[i]) * (xin-fld[i]);
	return(1);
      }
    }
    
    *current = 0.0;
    return(0);
    
    //------------------------------------------------------------------------
    
  }else if( mode  == mode_i_to_b ){
    
    //------------------------------------------------------------------------
    
    xin = *current;
    
    if( xin <= cur[0] ){
      *field = fld[0];
      goto G1000;
    }else if( xin >= cur[n-1] ){
      *field = fld[n-1];
      goto G1000;
    }
    
    for( i=0; i<(n-1); i++ ){
      if( (xin>cur[i]) && (xin<=cur[i+1]) ){
	*field = fld[i] + (fld[i+1]-fld[i]) / (cur[i+1]-cur[i]) * (xin-cur[i]);
	goto G1000;
      }
    }
    
    *field = 0.0;
    *kvalue = 0.0;
    return(0);
    
  G1000:
    *kvalue = (*field) / energy * eff_len * (double)0.3;
    return(1);
    
    //------------------------------------------------------------------------
    
  }else{
    
    return(0);
    
  }
}
