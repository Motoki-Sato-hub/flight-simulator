//F!============================================================================
//F!
//F	integer*4 function mag_ki_main
//F	1   ( mode, magname0, energy, kvalue, current, efflen, field )
//F!
//F!
//F!	search magnet data routines and get data
//F!
//F!	supplied
//F!	    mode    integer*4,	    access mode, =1 k-to-i, =2 i-to-k
//F!	    magname character*(*),  magnet name for a main coil.
//F!	    energy  real*4,	    beam energy in [gev].
//F!
//F!	returned
//F!	    efflen  real*4,	    effective length in [meter].
//F!
//F!	supplied or returned
//F!	    kvalue(2)	real*4,	    array for k values.
//F!	    current(2)	real*4,	    array for the setting current.
//F!	    field(2)	real*4,	    array for the field strength.
//F!
//F!	    for above array, 1 for main coil, 2 for trim coil. 
//F!
//F!	99/11/05 j.ozawa    add.    mag_ki_q_idx_skew
//F!       01/05/09 f.takagi   add.    mag_ki_q_ecube_skew
//F!============================================================================
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include "atf_common.h"
#include "magnet.h"

int mag_ki_main( int mode, const char *magname0, 
		 float energy, float kvalue[2],  float current[2], float *efflen, float field[2] )
{
  char magname[80];
  int status;
  
//  str_upcase( magname, sizeof(magname), magname0 );
    for( int i=0; i<(int)strlen(magname0)+1; i++ ) magname[i] = toupper(magname0[i]);
  printf("mode = %d\n", mode );
  printf("magname = %s\n", magname );
  printf("current[0] = %f\n", current[0] );
  printf("current[1] = %f\n", current[1] );
  printf("kvalue[0] = %f\n", kvalue[0] );
  printf("kvalue[1] = %f\n", kvalue[1] );
  printf("energy     = %f\n", energy );

  if( magname[0] == 'Z' ) goto G1000;
  if( magname[0] == 'B' ) goto G2000;
  if( magname[0]=='S' )
    {
      if( magname[3]=='R' ) goto SEXT_DR;
      goto SEXT_FF;
    }
  if( magname[0]=='Q' )
    {
      if( strncmp(magname,"QK",2)==0 ) goto SKEW_X;
      if( strncmp(magname,"QS",2)==0 ) goto SKEW_X;
      
      if( magname[3]=='R' ) goto QUAD_DR1;
      if( magname[4]=='R' ) goto QUAD_DR2;
      if( magname[3]=='X' ) goto QUAD_EXT1;
      if( magname[4]=='X' ) goto QUAD_EXT2;
      if( magname[3]=='F' ) goto QUAD_FF;
      if( magname[4]=='F' ) goto QUAD_FF;
      if( magname[5]=='F' ) goto QUAD_FF;
      if( magname[4]=='T' ) goto QUAD_BT;
    }
  
  goto UNDEFINED;
  
  //F!--------------------------------------------------------------------------
  
  //F! 2007-apr-04 add, zh101r
  //F! 2007-may-02 add, zh102r
  //F! 2007-dec-03 del, zh101r
  
 G1000:
  //F1000	if( magname(1:6).eq.'zh102r') then
  if( strncmp(magname,"ZH102R",6) == 0 ){
    
    //F		    mag_ki_main = mag_ki_z_nkk_cv
    //F	1           ( mode, energy, kvalue, current, efflen, field )
    status = mag_ki_z_nkk_cv( mode, energy, kvalue, current, efflen, field );
    //if( status !=1 ) atferr( "mag_ki_z_nkk_cv" );
    //F! 2007-dec-03 add, zh100r
    
    //F	else if( magname(1:6).eq.'zh100r') then
  }else if( strncmp(magname,"ZH100R",6) == 0 ){
    
    //F		    mag_ki_main = mag_ki_z_tecno_zh82982
    //F	1           ( mode, energy, kvalue, current, efflen, field )
    status = mag_ki_z_tecno_zh82982( mode, energy, kvalue, current, efflen, field );
    //if( status !=1 ) atferr( "mag_ki_z_tecno_zh82982" );
    //F! 2007-dec-03 add, zh101r
    
    //F	else if( magname(1:6).eq.'zh101r') then
  }else if( strncmp(magname,"ZH101R",6) == 0 ){
    
    //F		    mag_ki_main = mag_ki_z_tecno_zh82983
    //F	1           ( mode, energy, kvalue, current, efflen, field )
    status = mag_ki_z_tecno_zh82983( mode, energy, kvalue, current, efflen, field );
    //if( status !=1 ) atferr( "mag_ki_z_tecno_zh82983" );
    //F! 2007-dec-06 add, zv100r
    
    //F	else if( magname(1:6).eq.'zv100r') then
  }else if( strncmp(magname,"ZV100R",6) == 0 ){
    
    //F		    mag_ki_main = mag_ki_zv100r
    //F	1           ( mode, energy, kvalue, current, efflen, field )
    status = mag_ki_zv100r ( mode, energy, kvalue, current, efflen, field );
    //if( status !=1 ) atferr( "mag_ki_zv100r" );
    
    //F	else if( magname(2:2) .eq. 'v' ) then
  }else if( magname[1] == 'V' ){
    //F		 
    //F	    if( magname(3:5).eq.'10r' .or. 
    //F	1	magname(3:5).eq.'11r' .or.
    //F	1	magname(3:5).eq.'12r' .or.
    //F	1	magname(3:5).eq.'13r' .or.
    //F	1	magname(3:5).eq.'14r' .or.
    //F	1	magname(3:5).eq.'37r' .or.
    //F	1	magname(3:5).eq.'38r' .or.
    //F	1	magname(3:5).eq.'39r' .or.
    //F	1	magname(3:5).eq.'40r' .or.
    //F	1	magname(3:5).eq.'41r' ) then
    if( strncmp( &magname[2], "10R",3)==0 ||
	strncmp( &magname[2], "11R",3)==0 ||
	strncmp( &magname[2], "12R",3)==0 ||
	strncmp( &magname[2], "13R",3)==0 ||
	strncmp( &magname[2], "14R",3)==0 ||
	strncmp( &magname[2], "37R",3)==0 ||
	strncmp( &magname[2], "38R",3)==0 ||
	strncmp( &magname[2], "39R",3)==0 ||
	strncmp( &magname[2], "40R",3)==0 ||
	strncmp( &magname[2], "41R",3)==0 ){
      
      //F		    mag_ki_main = mag_ki_z_tecno_v58928
      //F	1           ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_tecno_v58928 ( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_tecno_v58928" );
      
      //F	    else if( 
      //F	1	magname(3:5).eq.'15r' .or. 
      //F	1	magname(3:5).eq.'34r' .or.
      //F	1	magname(3:5).eq.'35r' .or.
      //F	1	magname(3:5).eq.'36r' .or.
      //F	1	magname(3:5).eq.'42r' ) then
    }else if( 
	     strncmp( &magname[2], "15R",3)==0 ||
	     strncmp( &magname[2], "34R",3)==0 ||
	     strncmp( &magname[2], "35R",3)==0 ||
	     strncmp( &magname[2], "36R",3)==0 ||
	     strncmp( &magname[2], "42R",3)==0 ){
      
      //F		    mag_ki_main = mag_ki_z_tecno_v58284
      //F	1           ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_tecno_v58284( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_tecno_v58284" );
      //F	    else
    }else{
      //F		    mag_ki_main = mag_ki_z_nkk_cv
      //F	1           ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_nkk_cv( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_nkk_cv" );
      
      //F	    end if
    }
    //F	else if( magname(2:2) .eq. 'h' ) then
  }else if( magname[1] == 'H' ){
    
    //F	    if( magname(3:5).eq.'10r' .or. 
    //F	1	magname(3:5).eq.'11r' .or.
    //F	1	magname(3:5).eq.'12r' .or.
    //F	1	magname(3:5).eq.'13r' .or.
    //F	1	magname(3:5).eq.'14r' .or.
    
    //F	1	magname(3:5).eq.'34r' .or.
    //F	1	magname(3:5).eq.'35r' .or.
    //F	1	magname(3:5).eq.'36r' .or.
    //F	1	magname(3:5).eq.'37r' .or.
    //F	1	magname(3:5).eq.'38r' .or.
    //F	1	magname(3:5).eq.'39r' ) then
    if(    strncmp( &magname[2], "10R",3)==0  ||
	   strncmp( &magname[2], "11R",3)==0  ||
	   strncmp( &magname[2], "12R",3)==0  ||
	   strncmp( &magname[2], "13R",3)==0  ||
	   strncmp( &magname[2], "14R",3)==0  ||
	   strncmp( &magname[2], "34R",3)==0  ||
	   strncmp( &magname[2], "35R",3)==0  ||
	   strncmp( &magname[2], "36R",3)==0  ||
	   strncmp( &magname[2], "37R",3)==0  ||
	   strncmp( &magname[2], "38R",3)==0  ||
	   strncmp( &magname[2], "39R",3)==0  ){
      //F	
      //F		    mag_ki_main = mag_ki_z_tecno_h58283
      //F	1           ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_tecno_h58283( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_tecno_h58283" );

      //F! 2005-nov-18 add,
      
      //F	    else if( magname(3:6).eq.'100r' ) then
    }else if( strncmp( &magname[2], "100R",4)==0 ){
      
      //F		    mag_ki_main = mag_ki_z_tecno_h59789
      //F	1	    ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_tecno_h59789( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_tecno_h59789" );
      
      //F	    else
    }else{
      //F		    mag_ki_main = mag_ki_z_nkk_ch
      //F	1           ( mode, energy, kvalue, current, efflen, field )
      status = mag_ki_z_nkk_ch( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_z_nkk_ch" );
      
      //F	    end if
    }
    
    //F	end if
  }
  //F	return
  return( status );
  
  //F!--------------------------------------------------------------------------
  //F! main bending dipoles, 'b***'
  //F!--------------------------------------------------------------------------
  
 G2000:
  //F2000	if( magname(2:4) .eq. 'h1r' ) then
  if( strncmp( &magname[1], "H1R", 3 ) == 0 ){
    
    //F		mag_ki_main = mag_ki_b_melco
    //F	1		    ( mode, energy, kvalue, current, efflen, field )
    status = mag_ki_b_melco( mode, energy, kvalue, current, efflen, field );
    //if( status !=1 ) atferr( "mag_ki_b_melco" );
    
    //F	else
  }else{
    
    //F	    if( magname(2:2).eq.'v'   .or.
    //F	1	magname(2:4).eq.'h1x' ) then
    if( magname[1]=='V' ||  strncmp(&magname[1],"H1X",3)==0 ){
      //F	
      //Fc 99/12/02 j.ozawa
      //Fc		mag_ki_main = mag_ki_b_shi_h
      //Fc	1		    ( mode, energy, kvalue, current, efflen, field )
      
      //F		if(magname(1:2).eq.'bv') then
      if( strncmp(magname,"BV",2) == 0 ){
	//F		    mag_ki_main = mag_ki_b_shi_r
	//F	1		    ( mode, energy, kvalue, current, efflen, field )
	status = mag_ki_b_shi_r( mode, energy, kvalue, current, efflen, field );
	//if( status !=1 ) atferr( "mag_ki_b_shi_r" );
	//F		else
      }else{
	//F		    mag_ki_main = mag_ki_b_shi_h
	//F	1		    ( mode, energy, kvalue, current, efflen, field )
	status = mag_ki_b_shi_h( mode, energy, kvalue, current, efflen, field );
	//if( status !=1 ) atferr( "mag_ki_b_shi_h" );
	//F		endif
      }
      //Fcc
      //F	    else
    }else{
      
      //F		mag_ki_main = mag_ki_b_shi_c
      //F	1		    ( mode, energy, kvalue, current, efflen, field )
      status =mag_ki_b_shi_c( mode, energy, kvalue, current, efflen, field );
      //if( status !=1 ) atferr( "mag_ki_b_shi_c" );
      //F	    end if
    }
    
    //F	end if
  }
  
  //F	return
  return( status );
 
  //==================================================================================================
  
 SEXT_DR:
  
  status = mag_ki_s_melco( mode, energy, kvalue, current, efflen, field );   
  return( status );
  
 
 SEXT_FF:

  //2010/02/09 T.Yamauchi(KIS)
  //	add SF6FF, SF5FF, SD4FF, SF1FF, SD0FF      
  
  if( strncmp(magname,"SF6FF",5) == 0 ){
    
    status = mag_ki_q_sf6( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"SF5FF",5) == 0 ){
    
    status = mag_ki_q_sf5( mode, energy, kvalue, current, efflen, field );
    printf("kvalue[0] = %f\n", kvalue[0] );
  }else if( strncmp(magname,"SD4FF",5) == 0 ){
    
    status = mag_ki_q_sd4( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"SF1FF",5) == 0 ){
    
    status = mag_ki_q_sf1( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"SD0FF",5) == 0 ){
    
    status = mag_ki_q_sd0( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
  
  if( strncmp(magname,"SD",2) == 0 ) goto DEFOCUS;
  return( status );

 QUAD_BT:

  if( strncmp(magname,"QD10T",5)==0 || strncmp(magname,"QF10T",5 )==0 ||
      strncmp(magname,"QF30T_S",7)==0 || strncmp(magname,"QF31T_S",7)==0 ||
      strncmp(magname,"QD32T_S",7)==0 || strncmp(magname,"QF40T",5)==0 ||
      strncmp(magname,"QD40T",5)==0 || strncmp(magname,"QF41T",5)==0 ||
      strncmp(magname,"QD50T",5)==0 || strncmp(magname,"QF50T",5)==0 ||
      strncmp(magname,"QD51T",5)==0 || strncmp(magname,"QF11T",5)==0 ){
    status = mag_ki_q_tokin_3390a( mode, energy, kvalue, current, efflen, field );
  }else if( strncmp(magname,"QD11T",5)==0 || strncmp(magname,"QD30T_S",7)==0 ||
      strncmp(magname,"QD31T_S",7)==0 || strncmp(magname,"QD41T",5)==0 ||
      strncmp(magname,"QF42T",5)==0 || strncmp(magname,"QD42T",5)==0 ||
      strncmp(magname,"QF43T_S",7)==0 || strncmp(magname,"QF51T",5)==0 ||
      strncmp(magname,"QD52T",5)==0 || strncmp(magname,"QF53T",5)==0 ||
      strncmp(magname,"QD12T",5)==0 ){
    status = mag_ki_q_tokin_3390b( mode, energy, kvalue, current, efflen, field );
  }else if( strncmp(magname,"QF20T",5)==0 || strncmp(magname,"QF21T_S",7)==0 ||
      strncmp(magname,"QF52T",5)==0 ){
    status = mag_ki_q_tokin_3390c( mode, energy, kvalue, current, efflen, field );
  }else{
    goto UNDEFINED;
  }

  return( status );
  
  //==================================================================================================
  // Quadrupoles
  //==================================================================================================
  
 QUAD_DR1:
  if( strncmp(magname,"QF1R",4)==0 ){
      
    status = mag_ki_q_hitachi_1( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF2R",4) == 0 ){
    
    status = mag_ki_q_tokin_3325( mode, energy, kvalue, current, efflen, field );

  }else if( strncmp(magname,"QM1R",4)==0 ){
    //printf("energy''=%f\n", energy );
    status = mag_ki_q_hitachi_3( mode, energy, kvalue, current, efflen, field );
    //printf("status = %d\n", status );
  }else if( strncmp(magname,"QM2R",4)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;
    
  }else if( strncmp(magname,"QM3R",4)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;

  }else if( strncmp(magname,"QM4R",4)==0 ){
    
    status = mag_ki_q_hitachi_3( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM5R",4)==0 ){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM6R",4)==0 ){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;
    
  }else if( strncmp(magname,"QM7R",4)==0 ||
	    strncmp(magname,"QM8R",4)==0 ){
    
    status = mag_ki_q_tokin_3393( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM9R",4)==0){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;

  }else{
    goto UNDEFINED;
  }
    
  return( status );
  

 QUAD_DR2: 
  
  if( strncmp(magname,"QM7AR",5) == 0 ){
    
    status = mag_ki_q_tokin_3581( mode, energy, kvalue, current, efflen, field );  

  }else if( strncmp(magname,"QM9R",4)==0 ||
	    strncmp(magname,"QM10R",5)==0){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM11R",5)==0 ){
    
    status = mag_ki_q_tokin_3582( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;
    
  }else if( strncmp(magname,"QM12R",5)==0 ){
    
    status = mag_ki_q_qea13( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM13R",5)==0 ){
    
    status = mag_ki_q_qea19( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;
  
  }else if( strncmp(magname,"QM14R",5)==0 ){
    
    status = mag_ki_q_qea21( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM15R",5)==0 ){
    
    status = mag_ki_q_tokin_3582( mode, energy, kvalue, current, efflen, field );
    goto DEFOCUS;

  }else if( strncmp(magname,"QM16R",5)==0 ||
	    strncmp(magname,"QM17R",5)==0 ||
	    strncmp(magname,"QM18R",5)==0 ||
	    strncmp(magname,"QM19R",5)==0 ){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    if( strncmp(magname,"QM17R",5)==0 ||
	strncmp(magname,"QM18R",5)==0 ) goto DEFOCUS;
    
  }else if( strncmp(magname,"QM20R",5)==0 ){
    
    status = mag_ki_q_hitachi_3( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM21R",5)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
     goto DEFOCUS;
  
  }else if( strncmp(magname,"QM22R",5)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
     goto DEFOCUS;
 
  }else if( strncmp(magname,"QM23R",5)==0 ){
    
    status = mag_ki_q_hitachi_3( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
    
  return( status );
  

 SKEW_X: 
  
  if( strncmp(magname,"QK1X",4)==0 ||
      strncmp(magname,"QK2X",4)==0 ||
      strncmp(magname,"QK3X",4)==0 ||
      strncmp(magname,"QK4X",4)==0 ||
      strncmp(magname,"QS1X",4)==0 ||
      strncmp(magname,"QS2X",4)==0 ){
    
    status = mag_ki_q_idx_skew( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
    
  return( status );
  

 QUAD_EXT1: 
  
  if( strncmp(magname,"QF1X",4)==0 ||
      strncmp(magname,"QD2X",4)==0 ||
      strncmp(magname,"QF3X",4)==0 ||
      strncmp(magname,"QF4X",4)==0 ||
      strncmp(magname,"QD5X",4)==0 ||
      strncmp(magname,"QF6X",4)==0 ||
      strncmp(magname,"QD8X",4)==0 ||
      strncmp(magname,"QF9X",4)==0 ){
    
    status = mag_ki_q_hitachi_5( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF7X",4)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
  
  if( strncmp(magname,"QD",2)==0 ) goto DEFOCUS;
  return( status );
  

 QUAD_EXT2: 
  
  if( strncmp(magname,"QD10X",5) == 0 ){
    
    status = mag_ki_q_qea12( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF11X",5) == 0 ){
    
    status = mag_ki_q_qea11( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD12X",5) == 0 ){
    
    status = mag_ki_q_qea16( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF13X",5)==0 ||
	    strncmp(magname,"QD14X",5)==0 ||
	    strncmp(magname,"QF15X",5)==0 ){
    
    status = mag_ki_q_hitachi_4 ( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD16X",5) == 0 ){
    
    status = mag_ki_q_qea03( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF17X",5) == 0 ){
    
    status = mag_ki_q_qea06( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD18X",5) == 0 ){
    
    status = mag_ki_q_qea08( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF19X",5) == 0 ){
    
    status = mag_ki_q_qea10( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD20X",5)==0 ||
	    strncmp(magname,"QF21X",5)==0 ){
    
    status = mag_ki_q_hitachi_2( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
  
  if( strncmp(magname,"QD",2)==0 ) goto DEFOCUS;
  return( status );
  
  
 QUAD_FF: 
  
  if( strncmp(magname,"QM16FF",6) == 0 ){
    
    status = mag_ki_q_qea09( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM15FF",6) == 0 ){
    
    status = mag_ki_q_qea02( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM14FF",6) == 0 ){
    
    status = mag_ki_q_qea01( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM13FF",6) == 0 ){
    
    status = mag_ki_q_qea05( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM12FF",6) == 0 ){
    
    status = mag_ki_q_qea15( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QM11FF",6) == 0 ){
    
    status = mag_ki_q_qea04( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD10BFF",7) == 0 ){
    
    status = mag_ki_q_qea25( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD10AFF",7) == 0 ){
    
    status = mag_ki_q_qea22( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF9BFF",6) == 0 ){
    
    status = mag_ki_q_qea23( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF9AFF",6) == 0 ){
    
    status = mag_ki_q_qea24( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD8FF",5) == 0 ){
    
    status = mag_ki_q_qea20( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF7FF",5) == 0 ){
    
    status = mag_ki_q_qea27( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD6FF",5) == 0 ){
    
    status = mag_ki_q_qea26( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF5BFF",6) == 0 ){
    
    status = mag_ki_q_qea29( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF5AFF",6) == 0 ){
    
    status = mag_ki_q_qea28( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD4BFF",6) == 0 ){
    
    status = mag_ki_q_qea33( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD4AFF",6) == 0 ){
    
    status = mag_ki_q_qea30( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF3FF",5) == 0 ){
    
    status = mag_ki_q_qea31( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD2BFF",6) == 0 ){
    
    status = mag_ki_q_qea34( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD2AFF",6) == 0 ){
    
    status = mag_ki_q_qea32( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QF1FF",5) == 0 ){
    
    status = mag_ki_q_qf1( mode, energy, kvalue, current, efflen, field );
    
  }else if( strncmp(magname,"QD0FF",5) == 0 ){
    
    status = mag_ki_q_qd0( mode, energy, kvalue, current, efflen, field );
    
  }else{
    goto UNDEFINED;
  }
   
  if( strncmp(magname,"QD",2)==0 ) goto DEFOCUS;
  return( status );
  
  //================================================================================
  
 DEFOCUS:
  
  // Mode = I to K  
  if( mode == 2 ){
    kvalue[0] = -kvalue[0];
    kvalue[1] = -kvalue[1];
  }
  
  return( status );
  
 UNDEFINED:
  
  printf( "Magnet %s is not defined in mag_ki_main.\n", magname );
  return( 990 );
  
}
