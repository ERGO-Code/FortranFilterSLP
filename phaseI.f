c     phaseI solves the minimize total infeasibilities problem
c
c                min sum_i s^+_i+s^-_i,   
c                s.t. bl_i <= c(x) + s_i^+ - s_i^- <= bu_i   
c                            c(x) + s_i >= bl_i
c
c     again by SLP. That is at each iteration it solves the linearizat'n
c
c          min sum_i s^+_i + s^-_i,  
c          s.t.  
c             bl_i- c_(x) <= c'(x)*dx +s^+_i - s^-_i <= bu_i - c_i(x)
c
c
c     Phase I terminates if the objective value of a subproblem is 0
c     (infact <csr_eps) in which case we have a current point x and a 
c     step dx that satisfies
c
c               bl_i <= c_i(x) + c'(x)*dx <= bu_i
c
c     i.e. the step dx is feasible for the Phase II subproblem
c
c     Compared to the old way of doing this, this avoids having two
c     copies of A in the constraint matrix. On the dowside it introduced
c     additional variables.
c
c                      old            new
c     nvar             n+m           n+2*m   
c     ncon             2*m             m
c     nonzeros         2*ngrad+2*m   ngrad+2*m
c
c     dChange
c     also get sent Bf,Bhc,Bx,Bsv_bs_c,Bsv_bs_r
c     from slpmain.f
      subroutine phaseI(n, m, maxa, ngrad, maxlws, iprint,
     &     ifail, maxf, rho, f, x, xnew, c, cnew,
     &     drlo, drup, dclo, dcup, a, la, mrow, mcol, dels,
     &     dobj, d, ff, fc, lws, user, iuser,
     &     istat, rstat, blo, bup, lam, dspace,
     &     ispace, pr_to_fi, nspace, iter, s,
     &     sv_bs_c, sv_bs_r, msg, max_iter, hc, sv_t_fi,
     &     se_mx_lp_it, normd, Bf, Bhc, Bx, Bsv_bs_c, Bsv_bs_r)
c     end dChange
      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'msg.inc'

c     ----------- start of phase I --------------------------------
c     return codes
c
c     ifail = 0  - no problems
c     ifail = 1  - convergence to infeasible point
c     ifail = 2  - max number of slp iterations
c     ifail = 3  - computational error
c     ifail = 4  - max number of lp iterations
c     ifail = 5  - max number of filter entries
c     ifail = 6  - not enough workspace

c     ------------ declaration of passed parameters -----------------

      integer n, m, maxa, ngrad, maxlws, iprint, ifail,
     &        maxf, nspace, iter, max_iter

      double precision rho, f, hc, normd

      logical pr_to_fi, fst_lp, sv_t_fi, se_mx_lp_it

      character*20 msg



c     dChange
c     for new shape phase one 
c
c     change number of columns from n+m to n+2*m in
c     sv_bs_c, dclo, dcup, dobj
c
c     change number of rows from 2*m to m in
c     sv_bs_r, drlo, drup
c
c     change total row/col info from  2*ngrad+2*m to ngrad+2*m in
c     mrow, mcol,  dels
c
c     (don't change result or status)
c     endChange

      integer la(0:maxa+m+2), mrow(ngrad+2*m), mcol(ngrad + 2*m),
     &        lws(maxlws), istat(12), iuser(*),
     &        sv_bs_c(n+2*m), sv_bs_r(m),
     &        ispace(2*nspace)

      double precision x(n), xnew(n), c(m), cnew(m), drlo(m),
     &             drup(m), dclo(n+2*m), dcup(n+2*m), a(maxa),
     &             user(*), dels(ngrad+2*m),
     &             rstat(3), blo(n+m), bup(n+m), lam(n+m), d(n),
     &             ff(maxf), fc(maxf), dobj(n+2*m), dspace(nspace),
     &             s(n+m)

      integer        seed, bigiter
      common /summc/ seed, bigiter


c     ------------- declaration of internal variables ---------------

      integer i, j, flag, sizef, rtcod, pjp, k_el, mxa, nr_rsv, nr_lp

      double precision q, hcnew, fnew, ratiof, ratioc, ratio,
     &                 maxmu

      logical accept

      character ch

      integer trim77
      integer ems_time
      character*4 int2char
      character*32 bin32
      double precision ems_drand

      integer tt0, tt1
      common /ttc/ tt0, tt1

      integer    bt30_32
      parameter (bt30_32 = 2**31+2**30+2**29)

      double precision     infty, eps
      common /NLP_eps_inf/ infty, eps

      double precision stp_eps, csr_eps
      common /tol/     stp_eps, csr_eps

c     dChange
c     extra variables for various things
c     Bf, Bhc, Bx, Bsv_bs_c, Bsv_bs_r save the best point yet
c     bt30, bt31 are for warm starting basis lower and upper bounds
      double precision Bf, Bhc, Bx(n)
      integer  Bsv_bs_c(n), Bsv_bs_r(m)
      integer    bt30, bt31
      parameter (bt30 = 2**29, bt31=2**30) 
c     end dChange

      include 'EMSOLI.INC'
      include 'EMSOLN.INC'
      include 'EMSOLR.INC'

c     ... flush filter and set for phase I
      sizef = 1
      ff(sizef) = infty
      fc(sizef) = infty
      hc = infty
      fst_lp = .true.
      nr_rsv = 0

c     ... set new EMSOL tolerance
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
      RTOLPINF = csr_eps
      call ems_rset(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

 105  continue

      call gradient(n, m, mxa, x, a, la, maxa, user, iuser, flag)

c     dChange 
c     alter formation of phase one problem for new shape
      do i=1,n
         dclo(i) = max(-s(i)*rho, blo(i)-x(i))
         dcup(i) = min( s(i)*rho, bup(i)-x(i))
      end do
      do i=1,m
         dclo(i+n) = 0.
         dcup(i+n) = infty
         dclo(i+m+n) = 0
         dcup(i+m+n) = infty
      end do

      do i=1,m
         drlo(i) = blo(n+i)-c(i)
         drup(i) = bup(n+i)-c(i)
      end do

      do i=1,n
         dobj(i) = 0.
      end do
      do i=1,m
         dobj(n+i) = 1.
         dobj(n+m+i) = 1.
      end do

      k_el = 1
      pjp = la(0)
      do j=1,m
         do i=la(pjp+j), la(pjp+j+1)-1
            dels(k_el) = a(i)
            mcol(k_el) = la(i)
            mrow(k_el) = j
            k_el = k_el + 1
         end do
         dels(k_el) = -1.
         mcol(k_el) = n+j
         mrow(k_el) = j
         k_el = k_el + 1
         dels(k_el) = 1.
         mcol(k_el) = n+m+j
         mrow(k_el) = j
         k_el = k_el + 1
      end do
     
c      WRITE(nout,'(5(2I5,G15.5))')
c     &        (mrow(i), mcol(i), dels(i), i=1,ngrad+2*m)

      call ems_lmdl(rtcod, dspace, 1, m, n+2*m, ngrad+2*m, dobj, drlo,
     &              drup, dclo, dcup, mrow, mcol, dels)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_lmdl: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('LMDL',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 999
      end if

c     end dChange




c     dChange
c     minor rearrangement, this moved up out of set warmstart bracket
c     so gets done whether or not fst_lp is on
c     (as I always warmstart, either from phase one to phase one
c     or loading problem from phase two to phase one)
      call ems_nget(rtcod, dspace, emsoln, emsolnln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_nget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('NGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
c     end dChange

c     ---------------    set warmstart information -----------------

c     The 3 final bits of ispace(nrowstat/ncolstat) give the activity 
c     status of a row/column. For a warmstart this status needs to be set
c     to the correct code AND the corresponding vector of activities
c     (ispace(nrowact/ncolsol)) need to be set to the correct value
c     The activity status from the last LP is stored in sv_bs_r/c
c     The meaning of bits 30-32 is as follows
c         32:  1=basic, 0=nonbasic (for basic bits 30/31 are meaningless)
c         31:  1=may go down (i.e at upper bound and not fixed)
c              0=may not go down (i.e. at lower bound)
c         30:  1=may go up (i.e. at lower bound and not fixed)
c              0=may not go up (i.e. at upper bound)

c     dChange 
c     new look warmstart with new dimensions

c     This is a standard warmstart (Phase I -> Phase I)
      if (.not.fst_lp) then

         do i=1,m
c           clear bits 30-32 and replace by sv_bs_r
            ispace(nrowstat+i-1) =
     &           iand(ispace(nrowstat+i-1),not(bt30_32)) + sv_bs_r(i)
c           set row activities (to lower or upper bounds) if nonbasic
            if (.not.btest(sv_bs_r(i),31)) then
               if (.not.btest(sv_bs_r(i),30))
     &              dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
               if (.not.btest(sv_bs_r(i),29))
     &              dspace(nrowacts+i-1) = dspace(nrowupper+i-1)
            end if
         end do
c     do nothing much for columns
         do i=1,n+2*m
c           clear bits 30-32 and replace by sv_bs_c
            ispace(ncolstat+i-1) =
     &          iand(ispace(ncolstat+i-1),not(bt30_32)) + sv_bs_c(i)
c           keep activities where they are (why not for rows?)
         end do
      else


c     --- warmstart for the phase II -> phase I case ---
         if (nr_rsv==0) then
c     this is new warm start from phase two into phase one
c     we call eg the first block of columns c1
c     set 'c1', first block of columns 
c      (x variables from their final value in phase II)
            do i = 1,n        
               ispace(ncolstat+i-1) = 
     &              iand(ispace(ncolstat+i-1),not(bt30_32)) + sv_bs_c(i)
               if (.not.btest(sv_bs_c(i),31)) then 
                  if (.not.btest(sv_bs_c(i),30)) 
     &                 dspace(ncolsol+i-1) = dclo(i)
                  if (.not.btest(sv_bs_c(i),29))
     &                 dspace(ncolsol+i-1) = dcup(i) 
               end if
            end do
c     warmstart status for rows and the two slack columns are set at the same
c     time. Decision based on status of row in phase II problem.
c     set 'r1' rows, and 'c2' 'c3' identity columns
            do i = 1,m
c     nonbasic (active) and may not go up (i.e at upper bound)
               if (.not.btest(sv_bs_r(i),31) .and. 
     &              .not.btest(sv_bs_r(i),29)) then 
c               row nb at ub (or fixed): set r1 at ub, c2 at 0, c3 at 0
                  ispace(nrowstat+i-1) = 
     &              iand(ispace(nrowstat+i-1),not(bt30_32))+bt31
                  dspace(nrowacts+i-1)= dspace(nrowupper+i-1)
c                  set slack columns (both to lower bnd active at 0)
                  ispace(ncolstat+n+i-1) = 
     &              iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+i-1)=0
                  ispace(ncolstat+n+m+i-1) = 
     &              iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+m+i-1)=0
               end if
c     nonbasic (active) and may not go down (i.e at lower bound)
               if (.not.btest(sv_bs_r(i),31) .and. 
     &              .not.btest(sv_bs_r(i),30)) then 
c               row nb at lb (or fixed): set r1 at lb, c2 at 0, c3 at 0
                  ispace(nrowstat+i-1) = 
     &                 iand(ispace(nrowstat+i-1),not(bt30_32))+bt30
                  dspace(nrowacts+i-1)= dspace(nrowlower+i-1)
c                  set slack columns (both to lower bnd active at 0)
                  ispace(ncolstat+n+i-1) = 
     &                 iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+i-1)=0
                 ispace(ncolstat+n+m+i-1) = 
     &                 iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30
                 dspace(ncolsol+n+m+i-1)=0
               end if
c     nonbasic and may go neither up or down (fixed)
               if (.not.btest(sv_bs_r(i),31) .and. 
     &              .not.btest(sv_bs_r(i),30) .and.
     &              .not.btest(sv_bs_r(i),29)) then 
c                row fixed: set r1 at lb, c2 at 0, c3 at 0
                  ispace(nrowstat+i-1) = 
     &                 iand(ispace(nrowstat+i-1),not(bt30_32))
                  dspace(nrowacts+i-1)= dspace(nrowlower+i-1)
c                 set slack columns (both to lower bnd active at 0)
                  ispace(ncolstat+n+i-1) = 
     &                 iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+i-1)=0
                  ispace(ncolstat+n+m+i-1) = 
     &                 iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+m+i-1)=0
               end if
c     basic and above upper bound
               if (btest(sv_bs_r(i),31) .and. 
     &              dspace(nrowacts+i-1)>drup(i)+csr_eps) then 
c                 row basic above ub: set r1 at ub, c2 basic, c3 at 0
                  ispace(nrowstat+i-1) =
     &                 iand(ispace(nrowstat+i-1),not(bt30_32))+bt31
                  dspace(nrowacts+i-1) = dspace(nrowupper+i-1)
c                 set slack columns (s+ basic, s- active at 0)
                  ispace(ncolstat+n+i-1) = 
     &                 iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30_32
                  ispace(ncolstat+n+m+i-1) = 
     &                 iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+m+i-1)=0
               end if
c     basic and below lower bound
               if (btest(sv_bs_r(i),31) .and. 
     &              dspace(nrowacts+i-1)<drlo(i)-csr_eps) then 
c                 row basic below ub: set r1 at lb, c2 0, c3 basic
                  ispace(nrowstat+m+i-1) =
     &                 iand(ispace(nrowstat+m+i-1),not(bt30_32))+bt30   
                  dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
c                 set slack columns (s+ active at 0, s- basic)
                  ispace(ncolstat+n+i-1) = 
     &                 iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+i-1)=0
                  ispace(ncolstat+n+m+i-1) = 
     &               iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30_32
               end if
c     basic within bounds
               if (btest(sv_bs_r(i),31) .and. 
     &              dspace(nrowacts+i-1)<=drup(i)+csr_eps .and. 
     &              (dspace(nrowacts+i-1)>=drlo(i)-csr_eps) ) then 
c                   row basic in bounds: set r1 basic, c2 0, c3 0
                  ispace(nrowstat+i-1) = 
     &                 iand(ispace(nrowstat+i-1),not(bt30_32))+bt30_32
c                 set slack columns (both to lower bnd active at 0)
                  ispace(ncolstat+n+i-1) = 
     &                 iand(ispace(ncolstat+n+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+i-1)=0
                  ispace(ncolstat+n+m+i-1) = 
     &                 iand(ispace(ncolstat+n+m+i-1),not(bt30_32))+bt30
                  dspace(ncolsol+n+m+i-1)=0
               end if
            end do              
         end if                 
      end if                    
c     end dChange

c     -----------  end set warmstart information --------------------

c     ----- write model to file --------
      if (sv_t_fi) then
         write(nout,*) 
     &        'writing model to file: lp'//int2char(iter+1)//'.mps'
         open(12,file='lp'//int2char(iter+1)//'.mps')
c         open(12,file='lastlp.mps')
         call ems_bcdo(rtcod, dspace, 12, 1, 2)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_bcdo: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('BCDO',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 999
         end if
         close(12)
         open(12,file='lastas.dat')
         write(12,*) fst_lp
         do i=1,n+2*m
            write(12,'(A)') bin32(sv_bs_c(i))
         end do
         do i=1,m
            write(12,'(A)') bin32(sv_bs_r(i))
         end do
         close(12)
      end if

      if (fst_lp) fst_lp = .false.

c      call ems_mset(rtcod, dspace, 1, 0, 256, 0, 0, 9999, 1)
c      if (rtcod.ne.0) then
c         WRITE(nout,*)'error in ems_mset: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('MSET',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 999
c      end if

      call ems_scal(rtcod, dspace)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_scal: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('SCAL',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

c      call ems_mset(rtcod, dspace, 1, 0, -1, 0, 0, 9999, 1)
c      if (rtcod.ne.0) then
c         WRITE(nout,*)'error in ems_mset: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('MSET',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 999
c      end if

      call ems_sslv(rtcod, dspace, 1, 0)
      if (rtcod.ne.0) then
         if (iprint.ge.3) WRITE(nout,*)'error in ems_sslv: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('SSLV',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 999
      end if
      tt1 = ems_time()
c      WRITE(nout,*)tt1-tt0
      tt0=tt1


      call ems_nget(rtcod, dspace, emsoln, emsolnln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_nget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('NGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

c     --------------- get warmstart information ----------------------
c      do i=0,n+m-1
c         if (i.le.2*m-1) then
c            WRITE(nout,'(I,G20.10, F10.5, 3G20.10,F10.5,3G20.10)') i+1,
c     &    dspace(nobjective+i), dspace(ncolaux+i), dspace(ncollower+i),
c     &    dspace(ncolsol+i), dspace(ncolupper+i),
c     &    dspace(nrowaux+i), dspace(nrowlower+i), dspace(nrowacts+i),
c     &    dspace(nrowupper+i)
c         else
c            WRITE(nout,'(I,G20.10,F10.5,3G20.10)') i,
c     &    dspace(nobjective+i), dspace(ncolaux+i), dspace(ncollower+i),
c     &    dspace(ncolsol+i), dspace(ncolupper+i)
c         end if
c      end do

c     dChange
c     alter dimensions of basis save for new shape phase one
      do i=1,n+2*m
         sv_bs_c(i) = iand(ispace(ncolstat+i-1), bt30_32)
c         WRITE(nout,'(2b33.32)') sv_bs_c(i), ispace(ncolstat+i-1)
      end do
      do i=1,m
         sv_bs_r(i) = iand(ispace(nrowstat+i-1), bt30_32)
      end do
c     end dChange

c      do i=1,n+m
c         if (i.gt.2*m) then
c         WRITE(nout,'(I5,4G20.10)') i, dspace(nobjective+i-1),
c     &        dspace(ncollower+i-1), dspace(ncolupper+i-1),
c     &        dspace(ncolaux+i-1)
c         else
c         WRITE(nout,'(I5,7G20.10)') i, dspace(nobjective+i-1),
c     &        dspace(ncollower+i-1), dspace(ncolupper+i-1),
c     &        dspace(ncolaux+i-1),
c     &        dspace(nrowlower+i-1), dspace(nrowupper+i-1),
c     &        dspace(nrowaux+i-1)
c         end if
c      end do

c     -------------- END get warmstart information ------------------

      call ems_iget(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

      if (iprint.ge.2) WRITE(nout,*) 'LP exit code: ',iprobstat
      nr_lp = Iiternum
      
      if (iprobstat.eq.3) then
         WRITE(nout,*)' Maximal LP iterations reached', Iiternum
         WRITE(nout,*)
     &        ' Number of primal/dual infeasibilities and values:'
         WRITE(nout,*)'  primal:',inumpinf, rsumpinf
         WRITE(nout,*)'  dual  :',inumdinf, rsumdinf

c     dChange
cAGR  what does ems_itru.f do?
c     if we had ems_itru.f (which we dont here)
c     it is possible we have the message that max LP iterations are 
c     reached as we have reached a close to optimal point and ems_itru
c     has caused an exit with iprobstat=3
c     therefore we are more willing to carry on here than before
         if (2==0) then
            if (((rsumpinf<=50*rtolpinf).and.(rsumdinf<=50*rtoldinf)) 
     &           .or. (inumpinf.eq.0)) then
               goto 110
            end if
         else
            if (inumpinf.eq.0) then
c     ... boldly go...
               goto 110
            end if
         end if
c     end dChange

         msg = ' max lp iter'
         ifail = 4
         goto 118
      else if (iprobstat.eq.6) then
         WRITE(nout,*)' Not enough workspace for EMSOL'
         msg = ' out of memory'
         ifail = 6
         goto 118
c      else if (iprobstat.eq.2) then
c         write(*,*)'Here is the stop, that makes Integra freeze.....!'
c         stop
      elseif (iprobstat.ne.0) then
         WRITE(nout,*) ' Unexpected ifail from LP: ',iprobstat
         WRITE(nout,*) 
     &        ' Number of primal/dual infeasibilities and values:'
         WRITE(nout,*)'  primal:',inumpinf, rsumpinf
         WRITE(nout,*)'  dual  :',inumdinf, rsumdinf
         if (inumpinf.eq.0) then
            if ((rsumdinf.le.10.*rtoldinf).or.(nr_rsv.eq.2)) then
c           ... boldly go...
               goto 110
            end if
         end if
         if (nr_rsv.eq.0) then
            WRITE(nout,*) ' 1. resolve - coldstart'
            nr_rsv = 1
            goto 105
         end if
         if (nr_rsv.eq.1) then
            WRITE(nout,*)' 2. resolve - disturbing point'
            nr_rsv = 2
            do i=1,n
               x(i) = x(i) + 1.d-10*ems_drand(0)
            end do
            goto 105
         end if
         WRITE(nout,*) ' Maximal number of unsuccesful resolves!'
         WRITE(nout,*) ' STOP!'
         goto 999
      end if

c     ... LP solved succesfully (sort of)
 110  continue

      nr_rsv = 0

      if (se_mx_lp_it) then
         call ems_iget(rtcod, dspace, emsoli, emsoliln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_iget: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 999
         end if
         IMAXITER = 2*IITERNUM
         se_mx_lp_it = .false.
         call ems_iset(rtcod, dspace, emsoli, emsoliln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_iset: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('ISET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 999
         end if
      end if


      do i=1,n
         d(i) = dspace(ncolsol+i-1)
         lam(i) = dspace(ncolrcosts+i-1)
      end do
      do i=1,m
c         print '(I3,3F20.12)', i, dspace(nrowduals+i-1),
c     .        dspace(nrowduals+m+i-1), dspace(ncolsol+n+i-1)
         lam(n+i) = dspace(nrowduals+i-1)+ dspace(nrowduals+m+i-1)
c         if (dspace(ncolsol+n+i-1).le.1d-10) lam(n+i) = 0.d0
      end do


      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

      q = robjvalue

      iter = iter + 1

      normd = 0.
      do i=1,n
         xnew(i) = x(i) + d(i)
         normd = max(normd,abs(d(i)/s(i)))
      end do

      maxmu = 0.
      do i=1,n
         if (lam(i).le.0) then
            if (abs(dcup(i)-s(i)*rho).le.1.e-12)
     &                   maxmu = max(maxmu, abs(lam(i)))
         else
            if (abs(dclo(i)+s(i)*rho).le.1.e-12)
     &                   maxmu = max(maxmu, abs(lam(i)))
         end if
      end do

      if (iprint.ge.3)
     &    WRITE(nout,*) 'maximal multiplier on trust region: ',maxmu

      if (pr_to_fi) then
c         open(26,file='/home/agr/project/slp/slp'//int2char(iter)//'.m')
c         write(26, *) 'x = ['
c         do i=1,n
c            write(26,*) x(i)
c         end do
c         write(26,*) '];'
c         write(26, *) 'd = ['
c         do i=1,n
c            write(26,*) d(i)
c         end do
c         write(26,*) '];'
c         write(26,*) ' rho = ', rho
c         close(26)

         open(26,file='lastiter.dat')
         write(26,*) rho
         do i=1,n
            write(26,*) x(i)
         end do
         close(26)
      end if

      call objfun(xnew, n, fnew, user, iuser, flag)
      call confun(xnew, n, m, cnew, a, la, user, iuser, flag)

      hcnew = 0.
      do i=1,m
         hcnew = hcnew +
     &      max(cnew(i)-bup(n+i),blo(n+i)-cnew(i),0.d0)/s(n+i)
      end do

c     dChange
c     see if best point yet
      call checkBestPoint(n,m, csr_eps, fnew, hcnew, xnew, 
     &     ispace(ncolstat:ncolstat+n-1),ispace(nrowstat:nrowstat+m-1)
     &     ,Bf,Bhc,Bx,sv_bs_c,Bsv_bs_r)
c     end dChange

c     ----------------- check for convergence ---------------

      if (q.le.csr_eps) then
         if (iprint.ge.2) WRITE(nout,*) 'convergence'
         goto 115
      end if
      if ((rho.le.stp_eps).or.(q.eq.hc)) then
         ifail = 1
         goto 118
      end if
      if (iter.gt.max_iter) then
         if (iprint.ge.2) WRITE(nout,*)' maximal number of iterations'
         msg = msg(1:trim77(msg))//' max iter'
         ifail = 2
         goto 118
      end if
      if (iprint.ge.2) WRITE(nout,*) 'model hc:', q

      accept = .true.
c      write(nout,*) 'curr', hcnew, q
      do i=1,sizef
c         write(nout,*) i, fc(i), ff(i)
         if ((hcnew.ge.fc(i)).and.(q.ge.ff(i))) accept = .false.
      end do

      ratioc = (hcnew - hc)/(q - hc)
      ratiof = 1.
      ratio = max(ratiof, ratioc)

      if (iprint.ge.2) then
         WRITE(nout,*) 'reductions in model & reality + ratio:'
         WRITE(nout,'(A,3G20.10)') 'hc: ',q - hc, hcnew - hc, ratioc
      end if

      if (accept) then
         if (iprint.ge.2) WRITE(nout,*) 'step acceptable'
         if (iprint.ge.1) then
            WRITE(nout,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &           iter, normd, rho, hcnew,nr_lp, '+12'
            
            CALL LOGOUT(1,dspace,bigiter,
     &           iter,normd,rho,0D0,hcnew,nr_lp,+12)
         end if

         IF (msg_wrt_sumdat) THEN

            write(22,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &           iter, normd, rho, hcnew,nr_lp, '+12'

         END IF

         do i=1,sizef
            if ((hcnew.le.fc(i)).and.(q.le.ff(i))) then
               do j=i,sizef-1
                  fc(j) = fc(j+1)
                  ff(j) = ff(j+1)
               end do
               sizef = sizef - 1
            end if
         end do
         if (sizef.ge.maxf) then
            WRITE(nout,*)' Maximal filter size reached'
            msg = ' max filter'
            ifail = 5
            goto 118
         end if
         sizef = sizef + 1
         fc(sizef) = hcnew
         ff(sizef) = q
         hc = hcnew
         do i=1,n
            x(i) = xnew(i)
         end do
         do i=1,m
            c(i) = cnew(i)
         end do
         if ((ratio.gt.0.75).and.(abs(normd-rho).le.1.e-10)) then
            rho = 2.*rho
            if (iprint.ge.2) WRITE(nout,*) 'New trust region: ',rho
         end if
      else
         if (iprint.ge.2) WRITE(nout,*) 'step not acceptable'
         if (iprint.ge.1) then
            WRITE(nout,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &           iter, normd, rho, hcnew,nr_lp, '-12'

            CALL LOGOUT(1,dspace,bigiter,
     &           iter,normd,rho,0D0,hcnew,nr_lp,-12)
         end if
         
         IF (msg_wrt_sumdat) THEN

            write(22,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &           iter, normd, rho, hcnew,nr_lp, '-12'

         END IF

         rho = min(rho, normd)/4.
         if (iprint.ge.2)
     &      WRITE(nout,*) 'resolve with smaller trust region: ',rho
         if (rho.le.stp_eps) then
            ifail = 1
            goto 118
         end if
      end if

      goto 105

 115  continue

      if (iprint.ge.2) WRITE(nout,*) ' Feasible LP found => PHASE II'

      call objfun(x, n, f, user, iuser, flag)
      call confun(x, n, m, c, a, la, user, iuser, flag)
      hc = 0.
      do i=1,m
         hc = hc +
     &       max(c(i)-bup(n+i),blo(n+i)-c(i),0.d0)/s(n+i)
      end do
      if (iprint.ge.1) then
         WRITE(nout,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &        iter, normd, rho, hc,nr_lp, '+20'
         
         CALL LOGOUT(1,dspace,bigiter,
     &        iter,normd,rho,0D0,hc,nr_lp,+20)
      end if

      IF (msg_wrt_sumdat) THEN

         write(22,'(I3,2G16.8,14X,G14.8,I5,A5)')
     &        iter, normd, rho, hc,nr_lp, '+20'
         
      END IF


c     ... reset EMSOL tolerance
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if
      RTOLPINF = csr_eps*10.
      call ems_rset(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 999
      end if

      do i=1,ngrad+m
         if (mcol(i).le.n) then
            drup(mrow(i)) = drup(mrow(i)) - dels(i)*d(mcol(i))
         end if
      end do
c      WRITE(nout,*) csr_eps, q
      do i=1,m
c         WRITE(nout,*) i, drup(i),
c     . dspace(nrowacts+i-1)-dspace(nrowupper+i-1)+dspace(ncolsol+i+n-1),
c     .dspace(nrowacts+i+m-1)-dspace(nrowupper+i-1)-dspace(ncolsol+i+n-1)
         drup(1) = max(abs(drup(i)),drup(1))
      end do
c      WRITE(nout,*)'max c/s violation: ',drup(1)


      ifail = 0

      return

 118  continue
      WRITE(nout,*) 'Stop after ',iter,' Iterations.'
      WRITE(nout,*) 'No feasible point could be found.'

c      open(12,file='diet.sol')
c      write(12,*) 'values of variables'
c      do i=1,n
c         write(12,'(I4,G25.10)') i,x(i)
c      end do
c      write(12,*)
c      write(12,*) 'residuals of constraints'
c      do i=1,m
c         write(12,'(I4,G25.10)') i, c(i)
c      end do
c      close(12)
      hc = hcnew
      f = 0.d0
      do i=1,m
         f = f +
     &       max(c(i)-bup(n+i),blo(n+i)-c(i),0.d0)
      end do

      return

 999  continue
      ifail = 3
      hc = hcnew
      return

      end


      integer function trim77(string)

      character*(*) string
      integer l, i

      l = len(string)
      i = l

 5    continue
      if (string(i:i).ne.achar(32).and.string(i:i).ne.achar(9))
     &     goto 10
      i = i-1
      if (i.eq.0) goto 10
      goto 5

 10   continue

      trim77 = i

      end


      character*32 function bin32(n)

      integer i, n_rmd

      do i=1,32
         if (btest(n,i)) then
            bin32(i:i) = '1'
         else
            bin32(i:i) = '0'
         end if
      end do

      end
