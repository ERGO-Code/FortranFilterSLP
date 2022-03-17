        SUBROUTINE GETMIXERINCLUSION(MIXER,n,sn,rm,x,vr_reg,
     &                               INCLUSION)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int '08. All rights reserved |
*       +---------------------------------------------------+
*       | Subroutine :  GETMIXERINCLUSION                   |
*       | Description:  Extract the mixer inclusion         |
*       | Revision   :  June 2008                           |
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
        INTEGER          vr_reg(2+2*sn+3*n)
        DOUBLE PRECISION INCLUSION
        DOUBLE PRECISION x(n)
*
*---------------------
*  Local parameters  *
*---------------------
*
        PARAMETER (TYPE_U=4)
*
*----------------------------
*  Look for mixe inclusion  *
*----------------------------
*
        INCLUSION=0D0
        DO i=1,n
          IF (vr_reg(2+2*sn+i).eq.TYPE_U) THEN
            j = vr_reg(2+2*sn+n+i)
            IF (j-rm.EQ.MIXER) THEN
*
*  Total premix consumption
* 
              INCLUSION=x(i)
              RETURN
            END IF
          END IF
        END DO
        RETURN
        END
