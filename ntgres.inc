c     This common is to pass global variables between different
c     bits of Integra in order to gather all the results for
c     the NTGRESULT array.

      integer   MAX_SPS,
     &          MAX_RELATION,
     &          MAX_GIG,
     &          MAX_REL,
     &          MAX_INFEASIBLES,
     &          MAX_NUTS

      PARAMETER (MAX_SPS=200)   ! Defined in NTGMIX.INC
      PARAMETER (MAX_RELATION=100) ! Defined in NTGMIX.INC
      PARAMETER (MAX_GIG=MAX_SPS*MAX_RELATION)
      PARAMETER (MAX_REL=MAX_SPS*MAX_RELATION*2)
      PARAMETER (MAX_INFEASIBLES=15) ! Defined in NTGMIX.INC
      PARAMETER (MAX_NUTS=1000) ! Defined in NTGMIX.INC

c     these are here to pass information beween different parts of 
c     the code
      integer NUM_GIG_RESULTS,
     &     NUM_RATIO_RESULTS
      double precision FLOATING_BLANK
      PARAMETER (FLOATING_BLANK=1.0D-30)
      
c     These are here to ensure persistence across several calls to
c     getsensitivity
      
      integer GLO_INDEX_FLAG,
     &        REL_INDEX_FLAG,
     &        RATIO_NUTR_INDEX_FLAG

      integer NUM_RELATIONSHIPS
      integer GLOINDEX(MAX_GIG,2),
     &        RINDEX(MAX_REL),
     &        RATIO_NUTR_INDEX(MAX_NUTS,2)


      common /ntgrescommon/ NUM_GIG_RESULTS,
     &                      NUM_RATIO_RESULTS,
     &                      GLO_INDEX_FLAG,
     &                      REL_INDEX_FLAG,
     &                      RATIO_NUTR_INDEX_FLAG,
     &                      GLOINDEX,
     &                      RINDEX,
     &                      RATIO_NUTR_INDEX,
     &                      NUM_RELATIONSHIPS

