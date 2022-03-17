      subroutine ranging(dspace, ispace, nspace, ncols, nrows, status,
     .     dcup, drup, dclo, drlo, dobj, result, ifail)

      implicit none

        INCLUDE 'SLPCOMM.INC'

c     ------------ declaration of passed parameters -----------------

      integer nspace, ncols, nrows, ifail
      
      integer ispace(2*nspace), status(ncols+nrows)
      double precision dspace(nspace), dcup(ncols), drup(nrows),
     &     result(7,nrows+ncols), dclo(ncols), drlo(nrows),
     &     dobj(ncols)
      
c     ------------- declaration of internal variables ---------------

      double precision dzero, volume, dsmall
      parameter (dzero = 1d-10)
      parameter (volume = 1.d0)
      parameter (dsmall = 0.d0)

      integer i, j, k, l
      integer istat, jstat, kstat, lstat, rtcod
      double precision dlevel, dlowercost, dupperlevel, duppercost,
     .     dlowerlevel

      include 'EMSOLI.INC'
      include 'EMSOLN.INC'
      include 'EMSOLR.INC'

c     --------------------- subroutine body -------------------------
*
*-----------------------------
*  Get pointers to EMS data  *
*-----------------------------
*
      call ems_iget(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
*
*--------------------
*  Perform ranging  *
*--------------------
*
      call ems_rgda(rtcod, dspace)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rgda: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
c            CALL GET_EMSERROR('RGDA',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
*
*-----------------------------
*  Get pointers to EMS data  *
*-----------------------------
*
      call ems_nget(rtcod, dspace, emsoln, emsolnln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_nget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('NGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
*
*----------------------------------------
*  Start extracting column status bits  *
*----------------------------------------
*
      DO 160 J=1,NCOLS
*
*------------------------
*  Extract status bits  *
*------------------------
*     
         CALL TEST_BIT(ispace(emsoln(10)+j-1),1,ISTAT)
         CALL TEST_BIT(ispace(emsoln(10)+j-1),2,JSTAT)
         CALL TEST_BIT(ispace(emsoln(10)+j-1),3,KSTAT)
         CALL TEST_BIT(ispace(emsoln(10)+j-1),4,LSTAT)

c         write (*,'(I3,4I2)') j, LSTAT, KSTAT, JSTAT, ISTAT
*
*--------------------
*  Feasible column  *
*--------------------
*
c     Status = 1 : Variable is basic (in between bounds)
c              2 : Variable at lower bound
c              3 : Variable at upper bound
c              4 : Variable below lower bound
c              5 : Variable above upper bound
c
         IF (LSTAT.EQ.0) THEN
            IF (ISTAT.NE.0) THEN
               STATUS(J)=1
            ELSE IF (JSTAT.NE.0) THEN
               STATUS(J)=3
            ELSE IF (KSTAT.NE.0) THEN
               STATUS(J)=2
            ELSE
               IF (dspace(emsoln(9)+J-1).LT.0D0) THEN
                  STATUS(J)=3
               ELSE
                  STATUS(J)=2
               END IF
            END IF
*
*----------------------
*  Infeasible column  *
*----------------------
*
         ELSE
            IF(JSTAT.NE.0) THEN
               STATUS(J)=5
            ELSE IF(KSTAT.NE.0) THEN
               STATUS(J)=4
            ELSE
               IF(dspace(emsoln(7)+J-1).GT.
     &              dcup(J)) THEN
                  STATUS(J)=5
               ELSE
                  STATUS(J)=4
               END IF
            END IF
         END IF
*
*-----------------------------------------
*  Finish extracting column status bits  *
*-----------------------------------------
*
 160  CONTINUE
*
*-------------------------------------
*  Start extracting row status bits  *
*-------------------------------------
*
      DO 170 I=1,NROWS
*
*------------------------
*  Extract status bits  *
*------------------------
*
        K=NCOLS+I
        CALL TEST_BIT(ispace(emsoln(5)+i-1),1,ISTAT)
        CALL TEST_BIT(ispace(emsoln(5)+i-1),2,JSTAT)
        CALL TEST_BIT(ispace(emsoln(5)+i-1),3,KSTAT)
        CALL TEST_BIT(ispace(emsoln(5)+i-1),4,LSTAT)
*
*-----------------
*  Feasible row  *
*-----------------
*
        IF(LSTAT.EQ.0) THEN
           IF(ISTAT.NE.0) THEN
              STATUS(K)=1
           ELSE IF(JSTAT.NE.0) THEN
              STATUS(K)=3
           ELSE IF(KSTAT.NE.0) THEN
              STATUS(K)=2
           ELSE
              IF(dspace(emsoln(4)+I-1).GE.0D0) THEN
                 STATUS(K)=3
              ELSE
                 STATUS(K)=2
              END IF
           END IF
*
*-------------------
*  Infeasible row  *
*-------------------
*
        ELSE
           IF(JSTAT.NE.0) THEN
              STATUS(K)=5
           ELSE IF(KSTAT.NE.0) THEN
              STATUS(K)=4
           ELSE
              IF(dspace(emsoln(2)+I-1).GT.
     &             drup(i)) THEN
                 STATUS(K)=5
              ELSE
                 STATUS(K)=4
              END IF
           END IF
        END IF
*     
*--------------------------------------
*  Finish extracting row status bits  *
*--------------------------------------
*     
 170  CONTINUE
*
*------------------------------------
*  Start extracting column results  *
*------------------------------------
*
      DO 200 J=1,NCOLS
*     
*----------------------------
*  Extract column activity  *
*----------------------------
*
         DLEVEL=dspace(emsoln(7)+J-1)
*
         IF(STATUS(J).LE.3) THEN
            DLEVEL=MIN(MAX(DLEVEL, dclo(J)), dcup(J))
            IF(ABS(DLEVEL).LE.DZERO/VOLUME) THEN
               IF(STATUS(J).EQ.2) THEN
                  DLEVEL=0D0
               ELSE
                  IF(DLEVEL.EQ.0D0) THEN
                     DLEVEL=DSMALL/VOLUME
                  END IF
               END IF
            END IF
         END IF
*     
         RESULT(1,J)=DLEVEL
         RESULT(2,J)=dclo(J)
         RESULT(3,J)=dcup(J)
*
*------------------------------
*  Test if have done ranging  *
*------------------------------
*
c        IF(IFALL.NE.1) THEN
c                RESULT(4,J)=0D0
c                RESULT(5,J)=0D0
c                RESULT(6,J)=0D0
c                RESULT(7,J)=0D0
c                GOTO 200
c        END IF
*
*---------------------
*  Get ranging data  *
*---------------------
*
c     EMSOL vectors:
c        32 : lowest cost to maintain current basis
c        31 : highest cost to maintain current basis
c        39 : highest upper bound to maintain current basis
c        40 : lowest lower bound to maintain current basis
         
c     DLOWERCOST : maximal change of cost down
c     DUPPERCOST : maximal change of cost up
c     DLOWERLEVEL: maximal change of variable value up
c     DUPPERLEVEL: maximal change of variable value down

         DLOWERCOST = dobj(J)-dspace(emsoln(32)+J-1)
         DUPPERLEVEL= dspace(emsoln(39)+J-1)-DLEVEL
         DUPPERCOST = dspace(emsoln(31)+J-1)-dobj(J)
         DLOWERLEVEL= DLEVEL-dspace(emsoln(40)+J-1)
*
*-------------------
*  Basic variable  *  (in between bounds)
*-------------------
*
         IF(STATUS(J).EQ.1) THEN
            RESULT(4,J)=DLOWERCOST
            RESULT(5,J)=DUPPERLEVEL
            RESULT(6,J)=DUPPERCOST
            RESULT(7,J)=DLOWERLEVEL
*
*----------------
*  Lower bound  *
*----------------
*
         ELSE IF(STATUS(J).EQ.2) THEN
            RESULT(4,J)=DLOWERCOST
            RESULT(5,J)=DUPPERLEVEL
            RESULT(6,J)=DLOWERLEVEL
            RESULT(7,J)=0D0
*
*----------------
*  Upper bound  *
*----------------
*
         ELSE
            RESULT(4,J)=DUPPERCOST
            RESULT(5,J)=DLOWERLEVEL
            RESULT(6,J)=DUPPERLEVEL
            RESULT(7,J)=0D0
         END IF
*
*-------------------------------
*  Check if results too small  *
*-------------------------------
*
        DO 180 L=4,7
        IF (ABS(RESULT(L,J)).LE.DZERO/VOLUME) THEN
           RESULT(L,J)=0D0
        END IF
 180  CONTINUE
*
*-------------------------------------
*  Finish extracting column results  *
*-------------------------------------
*
 200  CONTINUE
*
*---------------------------------
*  Start extracting row results  *
*---------------------------------
*
      DO 220 I=1,NROWS
*
*-------------------------
*  Extract row activity  *
*-------------------------
*
         K=NCOLS+I
         DLEVEL=dspace(emsoln(2)+I-1)
         IF(STATUS(K).LE.3) THEN
            DLEVEL=MIN(MAX(DLEVEL,drlo(i)),drup(i))
         END IF
         RESULT(1,K)=DLEVEL
         RESULT(2,K)=drlo(i)
         RESULT(3,K)=drup(i)
*     
*------------------------------
*  Test if have done ranging  *
*------------------------------
*
c        IF(IFALL.NE.1) THEN
c                RESULT(4,K)=0D0
c                RESULT(5,K)=0D0
c                RESULT(6,K)=0D0
c                RESULT(7,K)=0D0
c                GOTO 220
c        END IF
*
*---------------------
*  Get ranging data  *
*---------------------
*
c     EMSOL ARRAYS:
c       32 : continuation of lowests costs to maintain current basis
c       31 : continuation of highest costs to maintain current basis
c       47 : highest row bound to maintain current basis
c       48 : lowest row bound to maintain current basis

c     DLOWERCOST : continuation of lowest cost to maintain basis
c     DUPPERCOST : continuation of highest cost to maintain basis
c     DLOWERLEVEL: maximum change of bounds down
c     DUPPERLEVEL: maximum change of bounds up


         DLOWERCOST = dspace(emsoln(32)+K-1)
         DUPPERLEVEL= dspace(emsoln(47)+I-1)-DLEVEL
         DUPPERCOST = dspace(emsoln(31)+K-1)
         DLOWERLEVEL= DLEVEL-dspace(emsoln(48)+I-1)
*
*-------------------
*  Basic variable  *
*-------------------
*
         IF(STATUS(K).EQ.1) THEN
            RESULT(4,K)=-DLOWERCOST
            RESULT(5,K)=DUPPERLEVEL
            RESULT(6,K)=DUPPERCOST
            RESULT(7,K)=DLOWERLEVEL
*     
*----------------
*  Lower bound  *
*----------------
*
         ELSE IF(STATUS(K).EQ.2) THEN
            RESULT(4,K)=-DLOWERCOST
            RESULT(5,K)=DUPPERLEVEL
            RESULT(6,K)=DLOWERLEVEL
            RESULT(7,K)=0D0
*
*----------------
*     Upper bound  *
*----------------
*     
         ELSE
            RESULT(4,K)=DUPPERCOST
            RESULT(5,K)=DLOWERLEVEL
            RESULT(6,K)=DUPPERLEVEL
            RESULT(7,K)=0D0
         END IF
         
*
*-------------------------------
*  Check if results too small  *
*-------------------------------
*
         DO 210 L=4,7
            IF(ABS(RESULT(L,K)).LE.DZERO) THEN
               RESULT(L,K)=0D0
            END IF
 210     CONTINUE
*
*----------------------------------
*  Finish extracting row results  *
*----------------------------------
*
 220  CONTINUE

c      do i=1,ncols
c         print '(I3,3G15.8,I4)', i, (result(j, i), j=1,3), status(i)
c      end do
      
      ifail = 0

      return
      
      
         
 999  continue

      ifail = 1

      return

      end
