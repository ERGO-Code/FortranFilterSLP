      subroutine do_range(dspace, ispace, nspace, n, m, status, dcup,
     &     drup, dclo, drlo, dobj, result, rifail, oldbnd, s, rho,
     &     x, c, ngrad, mrow, mcol, dels, sv_bs_r, sv_bs_c, iprint,
     &     oldrobjvalue, sv_t_fi)

      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'EMSOLI.INC'
      include 'EMSOLN.INC'
      include 'EMSOLR.INC'

c     ------------ declaration of passed parameters -----------------

      integer nspace, n, m, rifail, ngrad, iprint
      double precision rho, oldrobjvalue
      double precision range_tol

c     dChange
c     for new shape phase one 
c
c     change number of columns from n to n+2*m in
c     sv_bs_c(n+2*m), dcup(n+2*m), dclo
c
c     change total row/col info from  2*ngrad+2*m to ngrad+2*m in
c     mrow, mcol,  dels
c
c     (don't change result or status)
c     endChange

      integer ispace(2*nspace), status(n+m), mrow(ngrad+2*m),
     &     mcol(ngrad+2*m), sv_bs_r(m), sv_bs_c(n+2*m)
      double precision dspace(nspace), dcup(n+2*m), drup(m),
     &     result(7,n+m), dclo(n+2*m), drlo(m),
     &     dobj(n+2*m), x(n), c(m), oldbnd(2*n+2*m), s(n+m),
     &     dels(ngrad+2*m)
      logical sv_t_fi
c     -------------- declaration of local variables ------------------

      integer i, rtcod

      integer    bt30_32
      parameter (bt30_32 = 2**31+2**30+2**29)
c     -------------- declaration of external functions ---------------
      character*32 bin32
c     ==================== subroutine body ===========================

c     prepare for ranging:
c     ... setup last problem with original bounds

      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

c      write(nout,*) RTOLPINF, RTOLDINF
c      RTOLPINF = 100*RTOLPINF
c      RTOLDINF = 100*RTOLDINF

      call ems_rset(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

c     --- this assumes that last problem was Phase II
c
c     this is the Phase II problem (including trust region bounds)
c     but removing bounds that were tightened due to the presolve
c     oldbnd() are the bounds before presolve
c
c     AGR 23/04/13
c     This should now be redundant as untightened bounds are now
c     already substituted at the end of the final SLP solve
c
c     There are, however, still occasional ranging failures,
c     this is strange, since exactly the same problem as
c     at the end of slpmain is solved. This basically seems to
c     come down to a LP warmstart failure, even though the warmstart
c     has been successful in the last few SLP iterations

      do i=1,n
         dclo(i) = max(-rho*s(i), oldbnd(i)-x(i))
         dcup(i) = min( rho*s(i), oldbnd(n+m+i)-x(i))
c         if (max(-rho*s(i), oldbnd(i)-x(i)).ne.dclo(i)) then
c            print *, i, dclo(i), oldbnd(i)-x(i), -rho*s(i)
c         end if
c         if (min(rho*s(i), oldbnd(n+m+i)-x(i)).ne.dcup(i)) then
c            print *, i, dcup(i), oldbnd(n+m+i)-x(i), rho*s(i)
c         end if
      end do
      do i=1,m
         drlo(i) = oldbnd(n+i)-c(i)
         drup(i) = oldbnd(n+m+n+i)-c(i)
c         if (oldbnd(n+i)-c(i).ne.drlo(i)) then
c            print *, i, drlo(i), oldbnd(n+i)-c(i)
c         end if
c         if (oldbnd(n+m+n+i)-c(i).ne.drup(i)) then
c            print *, i, drup(i), oldbnd(n+m+n+i)-c(i)
c         end if
      end do

      call ems_lmdl(rtcod, dspace, 1, m, n, ngrad, dobj, drlo, drup,
     &     dclo, dcup, mrow, mcol, dels)
      if (rtcod.ne.0) then
         if (iprint.ge.3) WRITE(nout,*)'error in ems_lmdl: ',rtcod
         rifail = 1
         goto 99
      end if

c     ----- write model to file --------

      if (sv_t_fi) then
         open(12,file='lprange.mps')
         call ems_bcdo(rtcod, dspace, 12, 1, 2)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_bcdo: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('BCDO',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            stop
         end if
         close(12)
         print *, 'finished writing to file'

         open(12,file='asrange.dat')
         do i=1,n
            write(12,'(A)') bin32(sv_bs_c(i))
         end do
         do i=1,m
            write(12,'(A)') bin32(sv_bs_r(i))
         end do
          close(12)
         print *, 'finished writing to file'
      end if
c     -------------------- set warmstart information ----------------

      call ems_nget(rtcod, dspace, emsoln, emsolnln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_nget: ',rtcod
         rifail = 1
         goto 99
      end if

      do i=1,m
         ispace(nrowstat+i-1) =
     &        iand(ispace(nrowstat+i-1),not(bt30_32)) + sv_bs_r(i)
         if (.not.btest(sv_bs_r(i),31)) then
            if (.not.btest(sv_bs_r(i),30))
     &           dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
            if (.not.btest(sv_bs_r(i),29))
     &           dspace(nrowacts+i-1) = dspace(nrowupper+i-1)
         else
            dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
         end if
      end do
      do i=1,n
c     WRITE(nout,'(B33.32)') ispace(ncolstat+i-1)
         ispace(ncolstat+i-1) =
     &        iand(ispace(ncolstat+i-1),not(bt30_32)) + sv_bs_c(i)
c     WRITE(nout,'(2B33.32)') ispace(ncolstat+i-1), sv_bs_c(i)
         if (.not.btest(sv_bs_c(i),31)) then
c     c           ... variable is active
            if (.not.btest(sv_bs_c(i),30))
     &           dspace(ncolsol+i-1) = dclo(i)
            if (.not.btest(sv_bs_c(i),29))
     &           dspace(ncolsol+i-1) = dcup(i)
         else
            dspace(ncolsol+i-1) = dclo(i)
         end if
      end do

c     ------------- re-solve the modified last LP problem ------------

      call ems_sslv(rtcod, dspace, 1, 0)
      IF(RTCOD.GE.BAD_EMSERROR) THEN
         if (iprint.ge.3) WRITE(nout,*)'error in ems_sslv: ',rtcod
         rifail = 1
         goto 99
      end if

      call ems_iget(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iget: ',rtcod
         rifail = 1
         goto 99
      end if
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         rifail = 1
         goto 99
      end if

c     ----- write model to file --------

      if (sv_t_fi) then
         do i=1,n
            sv_bs_c(i) = iand(ispace(ncolstat+i-1), bt30_32)
         end do
         do i=1,m
            sv_bs_r(i) = iand(ispace(nrowstat+i-1), bt30_32)
         end do
         open(12,file='asrangesol.dat')
         do i=1,n
            write(12,'(A)') bin32(sv_bs_c(i))
         end do
         do i=1,m
            write(12,'(A)') bin32(sv_bs_r(i))
         end do
          close(12)
         print *, 'finished writing to file'
      end if

c      WRITE(nout,*) 'LP exit code (iters): ',iprobstat, Iiternum
c      write(nout,*) RTOLPINF, RTOLDINF
c      do i=1,n
c         write(*,'(I3,3G15.8,B5,1X,4L1)') i, 
c     &        dspace(ncollower+i-1), dspace(ncolsol+i-1),
c     &        dspace(ncolupper+i-1),
c     &        ishft(ispace(ncolstat+i-1),-28),
c     &        btest(ispace(ncolstat+i-1),28),
c     &        btest(ispace(ncolstat+i-1),29),
c     &        btest(ispace(ncolstat+i-1),30),
c     &        btest(ispace(ncolstat+i-1),31)
c      end do
c      print *, emsoln(10), ncolstat

c     max expected change in optimal value: the vertex shouldn't change
c     but the basis may. This seems the correct estimate of the possible
c     change (in primal/dual vectors rather than objective)
      range_tol = 2.*(RTOLPINF+RTOLDINF)

c      WRITE(nout,*)range_tol
      if (abs(robjvalue-oldrobjvalue).gt.range_tol) then
         WRITE(nout,*)'WARNING: Change of value in last LP'
         WRITE(nout,*)'=> Ranging data might be wrong'
         WRITE(nout,*)robjvalue, oldrobjvalue
         WRITE(nout,*)'Difference, tol', abs(robjvalue-oldrobjvalue), 
     &        range_tol
         rifail = 1
         goto 99
      end if
c      print *, robjvalue


c     ... call the ranging routine
      call ranging(dspace, ispace, nspace, n, m, status, dcup,
     &     drup, dclo, drlo, dobj, result, rifail)


 99    continue

      end
