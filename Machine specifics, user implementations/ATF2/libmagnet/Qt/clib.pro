#---------------------------------------------------------------
#
# Renewed on 2023 Aug.
#
#---------------------------------------------------------------

QT += core gui widgets xml network
CONFIG += c++17

TARGET = magnet
TEMPLATE = lib
CONFIG += staticlib

#LIBS += -lgfortran

INCLUDEPATH += $(ATF_OP_BASE)/common/atflib/

SOURCES += \
    ../mag_ki_b_melco.c \
    ../mag_ki_b_shi_c.c \
    ../mag_ki_b_shi_h.c \
    ../mag_ki_b_shi_r.c \
    ../mag_ki_calc.c \
    ../mag_ki_calc_d32t180.c \
    ../mag_ki_calc_liner.c \
    ../mag_ki_control.c \
    ../mag_ki_control_d32t180.c \
    ../mag_ki_main.c \
    ../mag_ki_q_d32t180.c \
    ../mag_ki_q_ecube_skew.c \
    ../mag_ki_q_hitachi_1.c \
    ../mag_ki_q_hitachi_2.c \
    ../mag_ki_q_hitachi_3.c \
    ../mag_ki_q_hitachi_4.c \
    ../mag_ki_q_hitachi_5.c \
    ../mag_ki_q_idx_skew.c \
    ../mag_ki_q_qd0.c \
    ../mag_ki_q_qd0_old.c \
    ../mag_ki_q_qea01.c \
    ../mag_ki_q_qea02.c \
    ../mag_ki_q_qea03.c \
    ../mag_ki_q_qea04.c \
    ../mag_ki_q_qea05.c \
    ../mag_ki_q_qea06.c \
    ../mag_ki_q_qea07.c \
    ../mag_ki_q_qea08.c \
    ../mag_ki_q_qea09.c \
    ../mag_ki_q_qea10.c \
    ../mag_ki_q_qea11.c \
    ../mag_ki_q_qea12.c \
    ../mag_ki_q_qea13.c \
    ../mag_ki_q_qea14.c \
    ../mag_ki_q_qea15.c \
    ../mag_ki_q_qea16.c \
    ../mag_ki_q_qea17.c \
    ../mag_ki_q_qea18.c \
    ../mag_ki_q_qea19.c \
    ../mag_ki_q_qea20.c \
    ../mag_ki_q_qea21.c \
    ../mag_ki_q_qea22.c \
    ../mag_ki_q_qea23.c \
    ../mag_ki_q_qea24.c \
    ../mag_ki_q_qea25.c \
    ../mag_ki_q_qea26.c \
    ../mag_ki_q_qea27.c \
    ../mag_ki_q_qea28.c \
    ../mag_ki_q_qea29.c \
    ../mag_ki_q_qea30.c \
    ../mag_ki_q_qea31.c \
    ../mag_ki_q_qea32.c \
    ../mag_ki_q_qea33.c \
    ../mag_ki_q_qea34.c \
    ../mag_ki_q_qf1.c \
    ../mag_ki_q_sd0.c \
    ../mag_ki_q_sd4.c \
    ../mag_ki_q_sf1.c \
    ../mag_ki_q_sf5.c \
    ../mag_ki_q_sf6.c \
    ../mag_ki_q_tokin_3325.c \
    ../mag_ki_q_tokin_3390a.c \
    ../mag_ki_q_tokin_3390b.c \
    ../mag_ki_q_tokin_3390c.c \
    ../mag_ki_q_tokin_3393.c \
    ../mag_ki_q_tokin_3581.c \
    ../mag_ki_q_tokin_3582.c \
    ../mag_ki_s_melco.c \
    ../mag_ki_s_slac_138.c \
    ../mag_ki_s_slac_213.c \
    ../mag_ki_z_nkk_ch.c \
    ../mag_ki_z_nkk_cv.c \
    ../mag_ki_z_tecno_h58283.c \
    ../mag_ki_z_tecno_h59789.c \
    ../mag_ki_z_tecno_v58284.c \
    ../mag_ki_z_tecno_v58928.c \
    ../mag_ki_z_tecno_zh82982.c \
    ../mag_ki_z_tecno_zh82983.c \
    ../mag_ki_zv100r.c  #\
    #../qkilnbt.f

HEADERS += \
    ../magnet.h \
    ../magnet_routines.h

