      SUBROUTINE COMPAREFLOATINGNUMBERS(FIRSTVALUE,SECONDVALUE,
     &                                    STATUS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '03. All rights reserved |
*       +---------------------------------------------------+
*       | Subroutine :  COMPAREFLOATINGNUMBERS              |
*       | Description:  Compare for floating equality       |
*       | Revision   :  November 2003                       |
*       +---------------------------------------------------+
*
*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT INTEGER (A-Z)
*
*------------------
*  Indirect data  *
*------------------
*
        DOUBLE PRECISION FIRSTVALUE,SECONDVALUE
*
*---------------
*  Local data  *
*---------------
*
        DOUBLE PRECISION DIFF,RATIO,ETOL
*       DATA ETOL /0.000001D0/
        DATA ETOL /0.0000000001D0/
*
*----------------------------------
*             STATUS              *
*----------------------------------
*                                 *
*  -1 = Firstvalue < Secondvalue  *
*   0 = Firstvalue = Secondvalue  *
*   1 = Firstvalue > Secondvalue  *
*------------------------------------
*
*-----------------------------
*  Test if values are equal  *
*-----------------------------
*
        IF(FIRSTVALUE.EQ.SECONDVALUE) THEN
          STATUS=0
          RETURN
        END IF
*--------------------------------
*  Test if first value is zero  *
*--------------------------------
*
        IF(FIRSTVALUE.EQ.0D0) THEN
          IF(DABS(SECONDVALUE).LT.ETOL) THEN
            STATUS=0
          ELSE IF(SECONDVALUE.LT.0D0) THEN
            STATUS=1
          ELSE
            STATUS=-1
          END IF
          RETURN
        END IF
*
*---------------------------------
*  Test if second value is zero  *
*---------------------------------
*
        IF(SECONDVALUE.EQ.0D0) THEN
          IF(DABS(FIRSTVALUE).LT.ETOL) THEN
            STATUS=0
          ELSE IF(FIRSTVALUE.LT.0D0) THEN
            STATUS=-1
          ELSE
            STATUS=1
          END IF
          RETURN
        END IF
*
*---------------------------
*  Perform magnitude test  *
*---------------------------
*
        DIFF=DABS(FIRSTVALUE-SECONDVALUE)
        IF(DIFF.LE.DABS(FIRSTVALUE)) THEN
          RATIO=DIFF/DABS(SECONDVALUE)
          IF(RATIO.LT.ETOL) THEN
            STATUS=0
            RETURN
          END IF
        END IF
*
*---------------------------
*  Magnitude test failure  *
*---------------------------
*
        IF(FIRSTVALUE.LT.SECONDVALUE) THEN
          STATUS=-1
        ELSE
          STATUS=1
        END IF
        RETURN
*        
        END

