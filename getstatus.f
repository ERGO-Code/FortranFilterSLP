        SUBROUTINE GETGIGSTATUS(ifail,
     &                          SENSITIVITY,
     &                          NTGRESULT,
     &                          GIGRESULT,
     &                          NRESULTS,
     &                          LEVEL,
     &                          PROD,
     &                          INGR,
     &                          NUMVIOLATED,
     &                          VIOLTYPE,
     &                          VIOLPROD,
     &                          VIOLNTRM,
     &                          VIOLVALUE)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : GETGIGSTATUS                         |
*       | Description: Set the STATUS for GIGs              |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*-------------------
*  Include common  *
*-------------------
*
        INCLUDE 'SLPCOMM.INC'
        INCLUDE 'pusr.inc'
*
*------------------
*  Indirect data  *
*------------------
*
        DIMENSION        NTGRESULT(*),
     &                   VIOLTYPE (*),
     &                   VIOLPROD (*),
     &                   VIOLNTRM (*),
     &                   VIOLVALUE(*)
*
        DOUBLE PRECISION SENSITIVITY,NTGRESULT,LEVEL,VIOLVALUE
*
*---------------
*  Local data  *
*---------------
*
        INTEGER GIGINDEX
        DOUBLE PRECISION DLOWER,DUPPER,DIFFLEVEL
*
*---------------
*  Initialise  *
*---------------
*
        TEMPSTATUS=0
*
*--------------------------------------------
*  Get Min and Max constraints for the GIG  *
*--------------------------------------------
*
        NDX=GIGINDEX(GIGRESULT,2,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
        DLOWER=NTGRESULT(NDX)
*                  
        NDX=GIGINDEX(GIGRESULT,3,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
        DUPPER=NTGRESULT(NDX)
*
*----------------------------------
*  Correct the level for the GIG  *
*----------------------------------
*
        DIFFLEVEL=0D0
*
*  Check if the GIG is violated.
*
        DO I=1,NUMVIOLATED
          IF(VIOLTYPE(I).EQ.TYPE_NU.AND.
     &       VIOLPROD(I).EQ.INGR.AND.
     &       VIOLNTRM(I).EQ.PROD) THEN
            DIFFLEVEL=VIOLVALUE(I)
          END IF
        END DO
*
        LEVEL=LEVEL+DIFFLEVEL
*
*-----------------------------
*  Set STATUS when feasible  *
*-----------------------------
*
        IF(ifail.EQ.0) THEN
          IF(SENSITIVITY.EQ.0D0) THEN
*  Within limits
            TEMPSTATUS=1
          ELSE IF(SENSITIVITY.GT.0D0) THEN
*  On minimum limit
            TEMPSTATUS=2
          ELSE IF(SENSITIVITY.LT.0D0) THEN
*  On maximum limit
            TEMPSTATUS=3
          END IF
*
*  In some cases the above method is not correct.....
*
          CALL COMPAREFLOATINGNUMBERS(DLOWER,DUPPER,ISIGN)
*
*  If min/max constraints are different, correct where level equals min/max
*
          IF(ISIGN.NE.0D0) THEN
            CALL COMPAREFLOATINGNUMBERS(LEVEL,DLOWER,ISIGN)
            IF(ISIGN.EQ.0D0) THEN
              TEMPSTATUS=2
            END IF
            CALL COMPAREFLOATINGNUMBERS(LEVEL,DUPPER,ISIGN)
            IF(ISIGN.EQ.0D0) THEN
              TEMPSTATUS=3
            END IF
          END IF
        ELSE
*
*-----------------------------
*  Set STATUS if infeasible  *
*-----------------------------
*
          IF(SENSITIVITY.EQ.0D0) THEN
*  Within limits
            TEMPSTATUS=1
          ELSE
*
*  Put the level into RESULT
*
            NDX=GIGINDEX(GIGRESULT,1,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                   NRESULTS)
            NTGRESULT(NDX)=LEVEL
*            
            IF(SENSITIVITY.GT.0D0) THEN
*  On minimum limit
              TEMPSTATUS=2
*
            ELSE IF(SENSITIVITY.LT.0D0) THEN
*  On maximum limit
              TEMPSTATUS=3
            END IF
*
*  In all cases
*            
            CALL COMPAREFLOATINGNUMBERS(LEVEL,DLOWER,ISIGN)
*  Under minimum limit
            IF(ISIGN.LT.0) THEN
              TEMPSTATUS=4
            END IF
            CALL COMPAREFLOATINGNUMBERS(LEVEL,DUPPER,ISIGN)
*  Above maximum limit
            IF(ISIGN.GT.0) THEN
              TEMPSTATUS=5
            END IF
          END IF
        END IF
*
*  Put the status into RESULT
*
        IF(TEMPSTATUS.EQ.0) THEN
          write(*,*)'FATAL ERROR !!. STATUS not found for a GIG!'
          STOP
        END IF
*        
        NDX=GIGINDEX(GIGRESULT,8,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
        NTGRESULT(NDX)=TEMPSTATUS
*
*  Put the level into RESULT
*
        NDX=GIGINDEX(GIGRESULT,1,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
        NTGRESULT(NDX)=LEVEL
        
*        
        RETURN
*        
        END





        SUBROUTINE GETRATIOSTATUS(ifail,
     &                            SENSITIVITY,
     &                            NTGRESULT,
     &                            RATIORESULT,
     &                            NRESULTS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Function   : GETRATIOSTATUS                       |
*       | Description: Set the STATUS for ratios            |
*       | Revision   : November 2003                        |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*-------------------
*  Include common  *
*-------------------
*
        INCLUDE 'SLPCOMM.INC'
        INCLUDE 'pusr.inc'
*
*------------------
*  Indirect data  *
*------------------
*
        DIMENSION        NTGRESULT(*)
*
        DOUBLE PRECISION SENSITIVITY,NTGRESULT
*
*---------------
*  Local data  *
*---------------
*
        DOUBLE PRECISION LEVEL,DLOWER,DUPPER,LEVEL1,LEVEL2,DIFFLEVEL
*
*-----------------------------
*  Set STATUS when feasible  *
*-----------------------------
*
        IF(ifail.EQ.0) THEN
          IF(SENSITIVITY.EQ.0D0) THEN
*  Within limits
            TEMPSTATUS=1
          ELSE IF(SENSITIVITY.GT.0D0) THEN
*  On minimum limit
            TEMPSTATUS=2
          ELSE IF(SENSITIVITY.LT.0D0) THEN
*  On maximum limit
            TEMPSTATUS=3
          END IF
          NDX=RATIOINDEX(RATIORESULT,8,nu,rm,mi,de,NRESULTS)
          IF(TEMPSTATUS.GT.NTGRESULT(NDX)) THEN
            NTGRESULT(NDX)=TEMPSTATUS
          END IF
*  Set sensitivity
          NDX=RATIOINDEX(RATIORESULT,4,nu,rm,mi,de,NRESULTS)
          NTGRESULT(NDX)=SENSITIVITY
*
*-----------------------------
*  Set STATUS if infeasible  *
*-----------------------------
*
        ELSE
*  Set odd status and let Integra-Mix sort the status out
          NDX=RATIOINDEX(RATIORESULT,8,nu,rm,mi,de,NRESULTS)
          NTGRESULT(NDX)=9999D0
*  Set sensitivity
          NDX=RATIOINDEX(RATIORESULT,4,nu,rm,mi,de,NRESULTS)
          NTGRESULT(NDX)=SENSITIVITY
        END IF
        RETURN
*        
        END
