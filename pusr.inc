c      *******************************************************
c      ** 
c      ** pointers into iuser & user
c      ** and problem dimensions
c      **
c      *******************************************************

c     set constant values

      integer TYPE_M, TYPE_MU, TYPE_U, TYPE_NU, TYPE_C, TYPE_PWL
      parameter (TYPE_M=1, TYPE_MU=2, TYPE_C=3, TYPE_U=4,
     &           TYPE_NU=5, TYPE_PWL=6)
      integer CTYPE_M, CTYPE_C, CTYPE_SUM, CTYPE_NURAT, CTYPE_OBJ,
     &        CTYPE_U, CTYPE_NU, CTYPE_SPECIAL, CTYPE_NONLINEAR,
     &        CTYPE_PWL
      parameter (CTYPE_OBJ=99, CTYPE_M=1, CTYPE_C=2, CTYPE_SUM=3,
     &           CTYPE_U=4, CTYPE_NU=5, CTYPE_NURAT=6, CTYPE_SPECIAL=7,
     &           CTYPE_NONLINEAR=8, CTYPE_PWL=9)

c     ------------- problem dimensions ------------------------

      integer N_RMINC
      integer N_NURATP
      integer N_NURATM
      integer N_RMAVA
      integer N_MU
      integer N_U
      integer N_NU
      integer N_PWL
      integer N_VR
      integer N_CS
      integer N_MAXA
      integer N_HESS
      integer NX_NTY
      integer NX_DNTY
      
c     ---------------- pointer into user ---------------------

c     begin-1 of c_i (prices for raw materials) (set to 0)

c
c     begin-1 of r_il (raw_mat-nutrient table) 
      integer PU_R

c     begin-1 of tonnages for demand
      integer PU_D0

c     begin-1 of constraint data in user array
      integer PU_CS

c       ------------ pointers into iuser ----------------------

      integer PI_CS
c     ... at iuser(pi_cs_pi+i) is start of int information about cs-i
      integer PI_CS_PI
c     ... at iuser(pi_cs_pu+i) is start of real information about cs-i
      integer PI_CS_PU
      integer PI_CS_DA
      integer PU_CS_DA

c     ---------------- pointers for bin-network --------------------

      integer BN_N_SRB   
      integer BN_N_SRBB        
      integer BN_P_SRB   
      integer BN_L_SRB   
      integer BN_N_TGB   
      integer BN_N_TGBB  
      integer BN_P_TGB   
      integer BN_L_TGB   
      integer BN_FOB  

c     ---------------- static variables for reg_cs/reg_vr ------------
c     ... flag to indicate that reg_cs has not been called yet
      integer REG_FST_CL
c     ... next constraint to be registered 
      integer REG_NX_CS
c     ... next variable to be registered
      integer REG_NX_VR

      common /pusrc/ N_RMINC, N_NURATP, N_NURATM, 
     &     N_RMAVA, N_MU, N_U, N_NU, N_VR, N_CS, N_MAXA, N_HESS,
     &     NX_NTY, NX_DNTY, 
     &     PU_R, PU_D0, PU_CS, 
     &     PI_CS, PI_CS_PI, PI_CS_PU, PI_CS_DA, PU_CS_DA,
     &     BN_N_SRB, BN_N_SRBB, BN_P_SRB, BN_L_SRB, BN_N_TGB, 
     &     BN_N_TGBB, BN_P_TGB, BN_L_TGB, BN_FOB,
     &     REG_FST_CL, REG_NX_CS, REG_NX_VR
