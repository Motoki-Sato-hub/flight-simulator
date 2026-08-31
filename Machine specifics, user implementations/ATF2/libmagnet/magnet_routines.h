#ifndef MAGNET_ROUTINES_H
#define MAGNET_ROUTINES_H
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
int mag_ki_calc( int mode, int n, float cur[], float fld[], float energy, float eff_len, 
        float *kvalue, float *field, float *current );

int mag_ki_main( int mode, const char *magname0, 
	float energy, float kvalue[2],  float current[2], float *efflen, float field[2] );
int mag_ki_z_nkk_cv( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_zh82982( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_zh82983( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_zv100r( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_v58928( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_v58284( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_h58283( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_tecno_h59789( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_z_nkk_ch( int mode, float energy, 
	float *kvalue, float *current, float *efflen, float *field );
int mag_ki_calc_liner( int mode, float energy, float slope, float efflen, 
			float *kvalue, float *current, float *field );
int mag_ki_b_shi_r( int mode, int energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_b_shi_h( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int  mag_ki_b_shi_c( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_s_slac_138( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_s_slac_213( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_s_melco( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3325( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3393( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3582( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3581( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_hitachi_1( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_hitachi_2( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_hitachi_3( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_hitachi_4( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_hitachi_5( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_idx_skew( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3390a( int mode, float energy,
        float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3390b( int mode, float energy,
        float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_tokin_3390c( int mode, float energy,
        float kvalue[2], float current[2], float *efflen, float field[2] );

int mag_ki_control(
	int     mode,	        // Access mode
	int     mx,	        // Number of data for excitation.
	float    cur[],	        // Current data array
	float    fld[],	        // Field data array
	float    energy,	// Energy in gev
	float    efflen,	        // Effective length (m)
	float    ncoil_main,   // Number of turn for main coil.
	float    ncoil_trim,    // Number of turn for trim coil.
	float    kvalue[2],     // Array of k values.
	float    field[2],        // Array of field
	float    current[2]    // Array of current
);
int mag_ki_q_d32t180( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_control_d32t180(
	     int mode,	    //! access mode
	     int mx,	    //! number of data for excitation.
	     float cur[],	    //! current data array
	     float fld[],	    //! field data array
	     float energy,	    //! energy in gev
	     float efflen,	    //! effective length (m)
	     int ncoil_main,    //! number of turn for main coil.
	     int ncoil_trim,    //! number of turn for trim coil.
	     float kvalue[2],	    //! array of k values.
	     float field[2],	    //! array of field
	     float current[2]	    //! array of current
);
int mag_ki_calc_d32t180( int mode, int n, float cur[], float fld[], 
	float energy, float eff_len, float *kvalue, float *field, float *current );
int mag_ki_q_ecube_skew( int mode, float *kvalue, float *current );

int mag_ki_b_melco( int mode, float energy,
		float kvalue[2], float current[2], float *efflen, float field[2] );

int mag_ki_q_qea12( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea11( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea16( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea03( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea06( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea08( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea10( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea09( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea02( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea01( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea05( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea15( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea04( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea25( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea22( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea23( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea24( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea20( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea27( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea26( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea29( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea28( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea33( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea30( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea31( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea34( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea32( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qf1( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qd0( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea13( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea19( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_qea21( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );

//2010/02/09 T.Yamauchi(KIS)
//	add SF6FF, SF5FF, SD4FF, SF1FF, SD0FF
int mag_ki_q_sf6( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_sf5( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_sd4( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_sf1( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
int mag_ki_q_sd0( int mode, float energy, 
	float kvalue[2], float current[2], float *efflen, float field[2] );
#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif

