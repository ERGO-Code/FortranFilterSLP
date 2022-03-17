*     This include file and the associated SLPCOMMON is for global
*     variables that are used within SLP. It is not to be used for
*     an interface with the outside world
*
*
*---------------
*  Parameters  *
*---------------
*
        INTEGER MAX_RMS
c     &          MAX_SPS,
c     &          MAX_RELATION,
c     &          MAX_GIG,
c     &          MAX_REL,
c     &          MAX_INFEASIBLES,
c     &          MAX_NUTS
*        
        PARAMETER (MAX_RMS=500)         ! Defined in NTGMIX.INC
c        PARAMETER (MAX_SPS=200)         ! Defined in NTGMIX.INC
c        PARAMETER (MAX_RELATION=100)    ! Defined in NTGMIX.INC
c        PARAMETER (MAX_GIG=MAX_SPS*MAX_RELATION)
c        PARAMETER (MAX_REL=MAX_SPS*MAX_RELATION*2)
c        PARAMETER (MAX_INFEASIBLES=15)  ! Defined in NTGMIX.INC
c        PARAMETER (MAX_NUTS=1000)       ! Defined in NTGMIX.INC

*
        CHARACTER*1 COMMA
        PARAMETER (COMMA=',')
*
        INTEGER MIN_DSPACE
        PARAMETER (MIN_DSPACE=1000)
*
        INTEGER BAD_EMSERROR
        PARAMETER (BAD_EMSERROR=300)

*
*-------------------------
*  SLP file name common  *
*-------------------------
*
        COMMON / SLPCOMMON /  dirname,
     &                        DIM_FILE,
     &                        LMTS_FILE,
     &                        INGRID_FILE,
     &                        BRIAN_FILE,
     &                        RAWMAT_FILE,
     &                        STRAIGHTS_FILE,
     &                        TONS_FILE,
     &                        PRM_FILE,
     &                        PRMS_FILE,
     &                        MIDELMT_FILE,
     &                        RMINC_FILE,
     &                        NURATP_FILE,
     &                        NURATM_FILE,
     &                        RMLMTS_FILE,
     &                        SPECS_FILE,
*
     &                        RNDSD_FILE,
     &                        EMSCT_FILE,
     &                        SLPCT_FILE,
     &                        SUMMARY_FILE,
     &                        SPECIAL_FILE,
     &                        CSTRPT_FILE,
     &	                      NONLINEAR_FILE
*
        CHARACTER*80 dirname,
     &               DIM_FILE,
     &               LMTS_FILE,
     &               INGRID_FILE,
     &               BRIAN_FILE,
     &               RAWMAT_FILE,
     &               STRAIGHTS_FILE,
     &               TONS_FILE,
     &               PRM_FILE,
     &               PRMS_FILE,
     &               MIDELMT_FILE,
     &               RMINC_FILE,
     &               NURATP_FILE,
     &               NURATM_FILE,
     &               SPECIAL_FILE,
     &	             NONLINEAR_FILE,
     &               RMLMTS_FILE,
     &               SPECS_FILE,
*
     &               RNDSD_FILE,
     &               EMSCT_FILE,
     &               SLPCT_FILE,
     &               SUMMARY_FILE,
     &               CSTRPT_FILE
*
*---------------------------
*  SLP problem dimensions  *
*---------------------------
*
        COMMON / SLPCOMMON / NUMBER_NUTRIENTS,
     &                        NUMBER_RAWMATERIALS,
     &                        NUMBER_SPECIFICATIONS,
     &                        NUMBER_PREMIXES
*
        INTEGER NUMBER_NUTRIENTS,
     &          NUMBER_RAWMATERIALS,
     &          NUMBER_SPECIFICATIONS,
     &          NUMBER_PREMIXES
*
*--------------------------
*  SLP control variables  *
*--------------------------
*
        COMMON / SLPCOMMON / EMSOL_ERRORMESSAGE,
     &	                      nout,
     &                        NUMBER_STARTS,
     &                        PROB_GLOBAL
        INTEGER nout,
     &          NUMBER_STARTS,
     &          PROB_GLOBAL
        CHARACTER*80 EMSOL_ERRORMESSAGE
