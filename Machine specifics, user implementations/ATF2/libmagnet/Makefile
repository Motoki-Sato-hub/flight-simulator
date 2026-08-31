OPT=-Wall -g

OBJS += mag_ki_main.o
OBJS += mag_ki_z_nkk_cv.o
OBJS += mag_ki_z_tecno_zh82982.o
OBJS += mag_ki_z_tecno_zh82983.o
OBJS += mag_ki_zv100r.o
OBJS += mag_ki_z_tecno_v58928.o
OBJS += mag_ki_z_tecno_v58284.o
OBJS += mag_ki_z_tecno_h58283.o
OBJS += mag_ki_z_tecno_h59789.o
OBJS += mag_ki_z_nkk_ch.o
OBJS += mag_ki_b_melco.o
OBJS += mag_ki_b_shi_r.o
OBJS += mag_ki_b_shi_h.o
OBJS += mag_ki_b_shi_c.o
OBJS += mag_ki_s_slac_138.o
OBJS += mag_ki_s_slac_213.o
OBJS += mag_ki_s_melco.o
OBJS += mag_ki_q_tokin_3325.o
OBJS += mag_ki_q_tokin_3393.o
OBJS += mag_ki_q_d32t180.o
OBJS += mag_ki_control_d32t180.o
OBJS += mag_ki_calc_d32t180.o
OBJS += mag_ki_q_tokin_3582.o
OBJS += mag_ki_q_hitachi_1.o
OBJS += mag_ki_q_hitachi_2.o
OBJS += mag_ki_q_hitachi_3.o
OBJS += mag_ki_q_hitachi_4.o
OBJS += mag_ki_q_hitachi_5.o
OBJS += mag_ki_q_idx_skew.o
OBJS += mag_ki_q_ecube_skew.o
OBJS += mag_ki_calc_liner.o
OBJS += mag_ki_control.o
OBJS += mag_ki_calc.o
OBJS += mag_ki_q_tokin_3581.o
OBJS += mag_ki_q_tokin_3390a.o
OBJS += mag_ki_q_tokin_3390b.o
OBJS += mag_ki_q_tokin_3390c.o
OBJS += mag_ki_q_qea13.o
OBJS += mag_ki_q_qea01.o
OBJS += mag_ki_q_qea02.o
OBJS += mag_ki_q_qea03.o
OBJS += mag_ki_q_qea04.o
OBJS += mag_ki_q_qea05.o
OBJS += mag_ki_q_qea06.o
OBJS += mag_ki_q_qea07.o
OBJS += mag_ki_q_qea08.o
OBJS += mag_ki_q_qea09.o
OBJS += mag_ki_q_qea10.o
OBJS += mag_ki_q_qea11.o
OBJS += mag_ki_q_qea12.o
OBJS += mag_ki_q_qea15.o
OBJS += mag_ki_q_qea16.o
OBJS += mag_ki_q_qf1.o
OBJS += mag_ki_q_qd0.o
OBJS += mag_ki_q_qea14.o
OBJS += mag_ki_q_qea17.o
OBJS += mag_ki_q_qea18.o
OBJS += mag_ki_q_qea19.o
OBJS += mag_ki_q_qea20.o
OBJS += mag_ki_q_qea21.o
OBJS += mag_ki_q_qea22.o
OBJS += mag_ki_q_qea23.o
OBJS += mag_ki_q_qea24.o
OBJS += mag_ki_q_qea25.o
OBJS += mag_ki_q_qea26.o
OBJS += mag_ki_q_qea27.o
OBJS += mag_ki_q_qea28.o
OBJS += mag_ki_q_qea29.o
OBJS += mag_ki_q_qea30.o
OBJS += mag_ki_q_qea31.o
OBJS += mag_ki_q_qea32.o
OBJS += mag_ki_q_qea33.o
OBJS += mag_ki_q_qea34.o
OBJS += mag_ki_q_sf6.o
OBJS += mag_ki_q_sf5.o
OBJS += mag_ki_q_sd4.o
OBJS += mag_ki_q_sf1.o
OBJS += mag_ki_q_sd0.o
OBJS += qkilnbt.o

all: libmagnet.a libmagnet.so

.c.o:
	gcc $(OPT) -c $< $(CFLAG)

.f.o:
	gfortran $(OPT) -c $< $(CFLAG)

libmagnet.a: $(OBJS)

	rm -f libmagnet.a
	ar q libmagnet.a $(OBJS)
	ranlib libmagnet.a

libmagnet.so: $(OBJS)
	gfortran -shared -fPIC -o libmagnet.so *.o

clean:
	rm -f *.o *.a *.so

