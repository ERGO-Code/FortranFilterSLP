      subroutine slp_main(n, m, maxa, ngrad, nspace, maxlws, iprint,
     &              ifail, rifail, maxf, rho, f, hc, x, xnew, c, cnew,
     &              drlo, drup, dclo, dcup, a, la,
     &              mrow, mcol, dels, dobj,
     &              d, ff, fc, lws, dspace,
     &              user, iuser, max_iter, istat, rstat, blo, bup, lam,
     &              ispace, sv_bs_r, sv_bs_c, s, sv_t_fi,
     &              se_mx_lp_it, wr_msg, status, result, oldbnd, 
     &              do_rng)

      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'msg.inc'
c     dChange 
c     add this to have access to pointers into user/iuser (for maxCost)
        include 'pusr.inc'
c     end dChange

c     ------------ declaration of passed parameters -----------------

      integer n, m, maxa, ngrad, nspace, maxlws, iprint, ifail,
     &        max_iter, maxf, rifail
      double precision rho, f, hc

c     dChange
c     new shape for phase one problem needs to be declared here
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
      
      integer la(0:maxa+m+2), mrow(ngrad+2*m), mcol(ngrad + 2*m),
     &        lws(maxlws), istat(12), iuser(*),
     &        ispace(2*nspace), sv_bs_c(n+2*m), sv_bs_r(m),
     &        status(n+m)
      double precision x(n), xnew(n), c(m), cnew(m), drlo(m),
     &             drup(m), dclo(n+2*m), dcup(n+2*m), a(maxa),
     &             dspace(nspace), user(*), dels(ngrad+2*m),
     &             rstat(3), blo(n+m), bup(n+m), lam(n+m), d(n),
     &             ff(maxf), fc(maxf), dobj(n+m), s(n+m),
     &             result(7, n+m), oldbnd(2*n+2*m)
c     end dChange

      logical sv_t_fi, se_mx_lp_it, wr_msg, do_rng
      integer f_lp_it
c     ------------- declaration of internal variables ---------------

      integer i, j, k, flag, mxa, iter, sizef, rtcod, pjp, INFifail,
     &        maxprt, k_el, nr_lp, nr_rsv

      integer n_two, n_xb

      double precision q, hcnew, fnew, ratiof, ratioc, ratio, normd,
     &                 penalty, mdl_f, mdl_hc, maxmu, kt_rsdu, mx_inf,
     &                 oldrobjvalue

      logical accept, f_prob, pr_to_fi, set_stat, fi_xst, f_lp_af_phI,
     &        scnd_try, rlx_bnds, f_rng_it

      character ch
      character*20 msg, ch_fnew
      character*4 int2char
      integer trim77
      integer ems_time
      character*32 bin32
      double precision ems_drand

      integer tt0, tt1, tt00
      common /ttc/ tt0, tt1

      integer    bt30_32
      parameter (bt30_32 = 2**31+2**30+2**29)


      double precision stp_eps, csr_eps
      common /tol/     stp_eps, csr_eps

      integer drv_ct(3)

      double precision     infty, eps
      common /NLP_eps_inf/ infty, eps

      logical wr_msg_ag
      common /wr_msg_ag_cn/ wr_msg_ag

      integer        seed, bigiter
      common /summc/ seed, bigiter

c     dChange
c     extra variables needed for various things:

c     maximal magnitude of objective coefficients and pointers into 
c     iuser/user needed to access this information on the problem
      double precision maxCost
      integer n_lt, n_nt, pcsi, p_ltd, p_ntd

c     kt_tol is estimate of sensible threshold on kt_rsdu
c     optimalityConditionMet is logical to say if converged
      double precision kt_tol
      logical optimalityConditionMet
      double precision oldKTmeasure

c     Storage for best point yet (with properties and optimal basis)
      double precision Bf, Bhc, Bx(n)
      integer  Bsv_bs_c(n), Bsv_bs_r(m)
c     end dChange
      double precision sv_obj


      include 'EMSOLI.INC'
      include 'EMSOLN.INC'
      include 'EMSOLR.INC'

c     --------------------- subroutine body -------------------------

      INTEGER NINF
      DOUBLE PRECISION SUMPRIMALINF,SUMDUALINF,SUMRELATIVEINF

c     dChange
c     size of initial trust region set to typically much larger than 5
         if (1==1) then
            rho = 0
            do i=1,n
               if ( ((bup(i)-blo(i))/s(i)) .ge. rho) then
                  rho = (bup(i)-blo(i))/s(i)
               end if
            end do
c     if too big just use 5
            if (rho>1e8) then 
               if (iprint>=1) print*, 'initial rho set too big at', rho,
     &              'reset to ', 5
               rho=5
            else
               rho = ceiling(rho) 
            end if
         end if
c     end dChange

         NINF=-1
         rifail = 1

         DO I=1,MIN(MIN_DSPACE,NSPACE)          ! Inserted JCH 7/3/2000
            DSPACE(I)=0D0                       ! Inserted JCH 7/3/2000
         END DO                                 ! Inserted JCH 7/3/2000

      iter = 0
      tt0 = ems_time()
      tt00 = tt0

      do i=1,n+m
         if (blo(i).gt.bup(i)) then
c            write(nout,*) i, blo(i), bup(i)
            fnew = 0.d0
            WRITE(nout,*) 'Inconsistent bounds'
            msg = ' incons bnds'
            normd = 0.d0
            hcnew = infty
            ifail = 7
            goto 15
         end if
      end do
      do i=1,n
c         if (x(i).gt.bup(i)+1d-7.or.x(i).lt.blo(i)-1d-7) then
c            print *, i, blo(i), x(i), bup(i)
c         end if
         x(i) = min(x(i), bup(i))
         x(i) = max(x(i), blo(i))
      end do

c      equivalence(ispace(1), dspace(1))

c      WRITE(nout,*)'EMSOL workspace :',nspace

      wr_msg_ag = .false.
      f_lp_af_phI = .false.
      scnd_try = .false.
      rlx_bnds = .false.
      f_rng_it = .false.

      msg = '.'
      maxprt = -1
      if (wr_msg) maxprt = 256
      ch = 'n'                  ! Corrected JCH 7/1/2000
      pr_to_fi = .false.
      if (ch.eq.'y') pr_to_fi = .true.
      
      if (iprint.gt.0) then
         WRITE(nout,*)
         WRITE(nout,*)' start solving problem by SLP method (v3_0)'
         WRITE(nout,*)'  n     = ',n
         WRITE(nout,*)'  m     = ',m
c     WRITE(nout,*)'  ngrad = ',ngrad
      end if

      flag = 0
      f_prob = .true.
      set_stat = .false.
      nr_rsv = 0
      ifail = 0

      fnew = infty

      call ems_mset(rtcod, dspace, 1, 0, maxprt, 0, 0, 9999, 1)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_mset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('MSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

CJCH      if (.not.f_prob) then

         call ems_init(rtcod, dspace)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_init: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('INIT',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if

CJCH      else
CJCH         f_prob = .false.
CJCH      end if

      call ems_iget(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if
      call ems_rd_ct_vr(rtcod, emsoliln, emsolrln, emsoli, emsolr)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rd_ct_vr: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('rd_ct_vr',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if
 
c     dChange
c     potential change to  EMSOL before first LP solve
c     so only calls ems_itru.f at start and when believed to be optimal
      if (2==0) then
         Iiterufreq = 999999       
      end if
c     end dChange

      call ems_iset(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('ISET',RTCOD,EMSOL_ERRORMESSAGE)
          END IF
         goto 99
      end if
      call ems_rset(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

c      call ems_rd_ct_vr(rtcod, dspace, drv_ct)
c      if (rtcod.ne.0) then
c         WRITE(nout,*)'error in ems_rd_ct_vr: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('rd_ct_vr',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 99
c      end if

      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         if (iprint.ge.3) WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if
c      WRITE(nout,*)'emsol_tol = ',csr_eps
      RTOLPINF = csr_eps *10.
      call ems_rset(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rset: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if


      call ems_dsca(rtcod, dspace, nspace, 1)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_dsca: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('DSCA',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 99
      end if

      call objfun(x, n, f, user, iuser, flag)
      call confun(x, n, m, c, a, la, user, iuser, flag)

      sizef = 1
      ff(sizef) = f
      hc = 0.
      do i=1,m
         hc = hc +
     &         max(c(i)-bup(n+i),blo(n+i)-c(i),0.d0)/s(n+i)
c         if (c(i).gt.bup(n+i)+1d-7.or.c(i).lt.blo(n+i)-1d-7) then
c            print '(A,I5,3g17.11)', 'C',i, blo(n+i), c(i), bup(n+i)
c         end if

      end do
c      stop

      fc(sizef) = hc

c     dChange
c     initialisation of best point vectors
      Bf = f
      Bhc = hc
      Bx = x
c     end dChange


c     dChange
c     get maxCost (i.e. the maximal objective coefficient) for use later
c     in KT test and in setting rpweight
c     Calculated from problem data => for linear objectives this bound 
c     is valid through all iterations. For bilinear objectives the 
c     calulated bound assumes that all variables are below 1.

c     .. get pointer into iuser to start of integer information about 
c        objective
      pcsi = iuser(pi_cs_pi+n_cs+1)
c     ... this slicing up explained in register.f:
c     number of linear and bilinear terms
      n_lt  = iuser(pcsi)
      n_nt  = iuser(pcsi+1)
c     pointer into user for real info on linear/bilinear terms
      p_ltd = iuser(pcsi+4) 
      p_ntd = iuser(pcsi+5)
      maxCost = 0            
      do i=1,n_lt
c     ... user(p_ltd+i-1) are the objective coefficients 
c         (for linear terms)
         if (abs(user(p_ltd+i-1))>maxCost) then
            maxCost=abs(user(p_ltd+i-1))
         end if
      end do
      do i=1,n_nt
c     ... user(p_ltd+i-1) are the objective coefficients 
c         (for bilinear terms) (are there any?)
         if (abs(user(p_ntd+i-1))>maxCost) then
            maxCost=abs(user(p_ntd+i-1))
         end if
      end do               
c     end dChange

     

      if (iprint.eq.1) then
         WRITE(nout,*)
         WRITE(nout,*) 'iter    |d|          rho           f   '//
     &     '        |c(x)|       lp_iter   st'
         WRITE(nout,*) 
     &      '------------------------------------------------'//
     &      '--------------------------'
      end if

      IF (msg_wrt_sumdat) THEN

         write(22,*) 'iter    |d|          rho           f   '//
     &        '        |c(x)|       lp_iter   st'
         write(22,*) '-----------------------------------------------'//
     &        '---------------------------'
         
      END IF


c     dChange
c     print out 0th iteration giving initial objective and infeas
      if (iprint.eq.1) then
         WRITE(nout,'(I3,2G16.8,2G14.8,I5,A5)')
     &        iter, 0.0, 0.0, f, hc, 0, ' 00'
      end if
      CALL LOGOUT(2,dspace,bigiter,
     &     iter,0,0,fnew,hcnew,0,+01)
c     end dChange

c      do i=1,n
c         WRITE(nout,*) i, x(i)
c      end do

      call gradient(n, m, mxa, x, a, la, maxa, user, iuser, flag)





 5    continue

c     ----------------- setup phase II LP ------------------------

      do i=1,n
         dclo(i) = max(-rho*s(i), blo(i)-x(i))
         dcup(i) = min( rho*s(i), bup(i)-x(i))
c         write(nout,*) i, blo(i), x(i), bup(i)
      end do
      do i=1,m
         drlo(i) = blo(n+i)-c(i)
         drup(i) = bup(n+i)-c(i)
c         print *, i, c(i)
      end do

      do i=1,n
         dobj(i) = a(i)
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
      end do

c      call pmdl(rtcod, dspace, 1, m, n, ngrad, dobj, drlo, drup,
c     &              dclo, dcup, mrow, mcol, dels)
      call ems_lmdl(rtcod, dspace, 1, m, n, ngrad, dobj, drlo, drup,
     &              dclo, dcup, mrow, mcol, dels)
      if (rtcod.ne.0) then
         if (iprint.ge.3) WRITE(nout,*)'error in ems_lmdl: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('LMDL',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         msg = msg(1:trim77(msg))//' bad model detected: data error?'
         goto 99
      end if

c     -------------------- set warmstart information ----------------
      if (set_stat) then
         call ems_nget(rtcod, dspace, emsoln, emsolnln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_nget: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('NGET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if

         do i=1,m
            ispace(nrowstat+i-1) =
     &          iand(ispace(nrowstat+i-1),not(bt30_32)) + sv_bs_r(i)
            if (.not.btest(sv_bs_r(i),31)) then
               if (.not.btest(sv_bs_r(i),30))
     &              dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
               if (.not.btest(sv_bs_r(i),29))
     &              dspace(nrowacts+i-1) = dspace(nrowupper+i-1)
            else
               dspace(nrowacts+i-1) = dspace(nrowlower+i-1)
            end if
         end do
         do i=1,n
c            WRITE(nout,'(B33.32)') ispace(ncolstat+i-1)
            ispace(ncolstat+i-1) =
     &          iand(ispace(ncolstat+i-1),not(bt30_32)) + sv_bs_c(i)
c            WRITE(nout,'(2B33.32)') ispace(ncolstat+i-1), sv_bs_c(i)
            if (.not.btest(sv_bs_c(i),31)) then
cc           ... variable is active
               if (.not.btest(sv_bs_c(i),30))
     &              dspace(ncolsol+i-1) = dclo(i)
               if (.not.btest(sv_bs_c(i),29))
     &              dspace(ncolsol+i-1) = dcup(i)
            else
               if (f_lp_af_phI) then
                  dspace(ncolsol+i-1) = x(i)+d(i)
               else
                  dspace(ncolsol+i-1) = dclo(i)
               end if
            end if
         end do
      end if

c     ----- write model to file --------

      if (sv_t_fi) then
c         open(12,file='lastlp.mps')
         open(12,file='lp'//int2char(iter+1)//'.mps')
         write(nout, *) 'write model to file: lp'//
     &        int2char(iter+1)//'.mps'
         call ems_bcdo(rtcod, dspace, 12, 1, 2)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_bcdo: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('BCDO',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if
         close(12)

c         open(12,file='lastas.dat')
         write(nout,*) 
     &        'write warmstart info to file: as'//int2char(iter+1)//
     &        '.dat   ',set_stat
         open(12,file='as'//int2char(iter+1)//'.dat')
         write(12,*) set_stat
         if (set_stat) then
            do i=1,n
               write(12,'(A)') bin32(sv_bs_c(i))
            end do
            do i=1,m
               write(12,'(A)') bin32(sv_bs_r(i))
            end do
         end if
         close(12)
      end if

      if (.not.set_stat) set_stat = .true.

c     ----------------- solve phase II LP -----------------------

c      call ems_mset(rtcod, dspace, 1, 0, 256, 0, 0, 9999, 1)
c      if (rtcod.ne.0) then
c         WRITE(nout,*)'error in ems_mset: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('MSET',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 99
c      end if

      call ems_scal(rtcod, dspace)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_scal: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('SCAL',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

c      call ems_mset(rtcod, dspace, 1, 0, -1, 0, 0, 9999, 1)
c      if (rtcod.ne.0) then
c         WRITE(nout,*)'error in ems_mset: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('MSET',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 99
c      end if

      call ems_sslv(rtcod, dspace, 1, 0)
      if (rtcod.ne.0) then
         if (iprint.ge.3) WRITE(nout,*)'error in ems_sslv: ',rtcod
C         IF(RTCOD.GE.BAD_EMSERROR) THEN
C            CALL GET_EMSERROR('SSLV',RTCOD,EMSOL_ERRORMESSAGE)
C         END IF
c         goto 99
      end if
      tt1 = ems_time()
c      WRITE(nout,*)tt1-tt0
      tt0=tt1

                CALL EMS_CA_G_N_SU_MX_PR_IFS(NINF,
     &                                       SUMPRIMALINF,
     &                                       SUMDUALINF,
     &                                       SUMRELATIVEINF,
     &                                       DSPACE,ISPACE)

      call ems_nget(rtcod, dspace, emsoln, emsolnln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_nget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('NGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

c     ---------------- get warmstart information --------------------
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


      do i=1,n
         sv_bs_c(i) = iand(ispace(ncolstat+i-1), bt30_32)
c         WRITE(nout,'(2b33.32)') sv_bs_c(i), ispace(ncolstat+i-1)
      end do
      do i=1,m
         sv_bs_r(i) = iand(ispace(nrowstat+i-1), bt30_32)
      end do

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

c      open(11,file='lpas'//int2char(iter)//'.dat')
c      do i=1,n
c         write(11,'(b33.32)') sv_bs_c(i)
c      end do
c      do i=1,m
c         write(11,'(b33.32)') sv_bs_r(i)
c      end do
c      close(11)


c     -------------- END get warmstart information ------------------


      call ems_iget(rtcod, dspace, emsoli, emsoliln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_iget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if
      call ems_rget(rtcod, dspace, emsolr, emsolrln)
      if (rtcod.ne.0) then
         WRITE(nout,*)'error in ems_rget: ',rtcod
         IF(RTCOD.GE.BAD_EMSERROR) THEN
            CALL GET_EMSERROR('RGET',RTCOD,EMSOL_ERRORMESSAGE)
         END IF
         goto 99
      end if

      if (iprint.ge.2) WRITE(nout,*) 'LP exit code: ',iprobstat
      nr_lp = Iiternum

      if (iprobstat.eq.1) then
         if (iprint.ge.2)
     &        WRITE(nout,*)' Infeasible LP detected => PHASE I'
         mx_inf = 0.d0
         do i=1,n
            mx_inf = max(mx_inf, max(
     &   (dspace(ncolsol+i-1)-dcup(i))*dspace(ncolscales+i-1),
     &   (dclo(i)-dspace(ncolsol+i-1))*dspace(ncolscales+i-1)))
c            WRITE(nout,*)i, dspace(ncolscales+i-1)
         end do
         do i=1,m
            mx_inf = max(mx_inf, max(
     &   (dspace(nrowacts+i-1)-drup(i))*dspace(nrowscales+i-1),
     &   (drlo(i)-dspace(nrowacts+i-1))*dspace(nrowscales+i-1)))
c            WRITE(nout,*) i, dspace(nrowscales+i-1)
         end do
c         WRITE(nout,*)'max c/s violation: ',mx_inf

         if (f_lp_af_phI) then
            if (scnd_try) then
               msg = msg(1:trim77(msg))//' failed feas'
               ch_fnew = ' '
               hcnew = hc
               ifail = 3
               RTOLPINF = csr_eps*10.
               goto 17
            end if
            WRITE(nout,*)'First Lp after phase I not feasible'
            WRITE(nout,*)'Relax tolerance to ',mx_inf*1.1

            RTOLPINF = mx_inf*1.1
            call ems_rset(rtcod, dspace, emsolr, emsolrln)
            if (rtcod.ne.0) then
               WRITE(nout,*)'error in ems_rset: ',rtcod
               IF(RTCOD.GE.BAD_EMSERROR) THEN
                  CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
               END IF
               goto 99
            end if
            do i=1,n
               sv_bs_c(i) = iand(ispace(ncolstat+i-1), bt30_32)
               d(i) = dspace(ncolsol+i-1)
            end do
            do i=1,m
               sv_bs_r(i) = iand(ispace(nrowstat+i-1), bt30_32)
            end do

            scnd_try = .true.
            goto 5
         end if

c     dChange
c     do five things after p2 solve point before entering p1
c     
c     get info on the point
c     see if best point yet
c     increment iteration counter
c     print out iteration with 01 code
c     check if SLP iteration limit reached

         if (1==1) then

c     get info on the point
            normd = 0.d0
            do i = 1,n
               d(i) = dspace(ncolsol+i-1)
               lam(i) = dspace(ncolrcosts+i-1) 
               xnew(i) = min(max(blo(i), x(i) + d(i)), bup(i))
               normd = max(normd,abs(d(i)/s(i))) 
            end do
            do i = 1,m
               lam(n+i) = dspace(nrowduals+i-1) 
            end do
            call objfun(xnew, n, fnew, user, iuser, flag)
            call confun(xnew, n, m, cnew, a, la, user, iuser, flag)
            hcnew = 0.
            do i = 1,m
               hcnew = hcnew + 
     .              max(cnew(i)-bup(n+i),blo(n+i)-cnew(i),0.d0)/s(n+i)
            end do  
            
c     see if best point yet
            call checkBestPoint(n,m, csr_eps, fnew, hcnew, xnew, 
     &  ispace(ncolstat:ncolstat+n-1),ispace(nrowstat:nrowstat+m-1)
     &           ,Bf,Bhc,Bx,sv_bs_c,Bsv_bs_r)
            
c     increment iteration counter
            iter = iter+1
            
c     print out iteration with 01 code
            if (iprint.ge.1) then
               WRITE(nout,'(I3,2G16.8,2G14.8,I5,A5)')
     &              iter, normd, rho, fnew, hcnew,nr_lp, ' 01'
            end if
            CALL LOGOUT(2,dspace,bigiter,
     &           iter,normd,rho,fnew,hcnew,nr_lp,+01)

c     check if SLP iteration limit reached
            if (iter>=max_iter) then
               if (iprint>=2) 
     &              write(nout,*) ' maximal number of iterations'
               msg = msg(1:trim77(msg))//' max iter'
               ifail = 2
               goto 15
            end if

         end if
c     end dChange

c     dChange
c     also send Bf, Bhc, Bx, Bsv_bs_c, Bsv_bs_r to phaseI
c     for keeping track of best point
         call phaseI(n, m, maxa, ngrad, maxlws, iprint,
     &        ifail, maxf, rho, f, x, xnew, c, cnew,
     &        drlo, drup, dclo, dcup, a, la, mrow, mcol, dels,
     &        dobj, d, ff, fc, lws, user, iuser,
     &        istat, rstat, blo, bup, lam,dspace, ispace,
     &        pr_to_fi, nspace, iter, s,
     &        sv_bs_c, sv_bs_r, msg, max_iter, hc, sv_t_fi,
     &        se_mx_lp_it, normd, Bf, Bhc, Bx, Bsv_bs_c, Bsv_bs_r)
c     end dChange

         if (ifail.ne.0) then
c     ... phaseI returns in f the best obtained constraint violation
c         and in hc the scaled constraint violation
            msg = msg(1:trim77(msg))//' infeasible'
            ch_fnew = ''
c     ... remember best found c/s violation (unscaled) for reporting
            hc = f
            hcnew = hc
chc   ... Do ranging even if infeasible
            call ranging(dspace, ispace, nspace, n, m, status, dcup,
     &           drup, dclo, drlo, dobj, result, rifail)
            goto 17
         end if

c     dChange
c     add multiplier to LP phase one for SLP phase two
c     size of multiplier is 0.01 scaled by maxCost
      if (1==1) then
         rpweight = 0.01/maxCost 
         call ems_rset(rtcod, dspace, emsolr, emsolrln)  
         if (iprint>=2) print*, 'setting rpweight', rpweight
      end if
c     end dChange

         f_lp_af_phI = .true.
         scnd_try = .false.
         n_two = 0
c     ... set status vector for phase II hotstart:

c     dChange
c     with new phase one have a different warm start procedure
c     for reentering phase two from phase one
c     basically dont need to do anything:
c       If a constraint is active (by necessity with zero slacks) in 
c       Phase I it is also active in the corresponding phase II problem.
         do i=1,m
            sv_bs_r(i) = ibset(sv_bs_r(i),31)
         end do
c     end dChange

         set_stat = .true.
         sizef = 1
         ff(sizef) = f
         fc(sizef) = hc
         goto 5
      elseif (iprobstat.eq.3) then
         WRITE(nout,*)' Maximal LP iterations reached'
         WRITE(nout,*)
     &        ' Number of primal/dual infeasibilities and values:'
         WRITE(nout,*)'  primal:',inumpinf, rsumpinf
         WRITE(nout,*)'  dual  :',inumdinf, rsumdinf

c     dChange
c     if we have ems_itru.f (but we dont here)
c     it is possible we have the message that max LP iterations are 
c     reached as we have reached a close to optimal point and ems_itru
c     has caused an exit with iprobstat=3
c     if that is the case then we continue problem just as before.
c     In that case we are likely more likely to reach here
c     end dChange

         if ((inumpinf.eq.0).or.(rsumpinf.le.50.*rtolpinf)) then
c        ... boldly go...
            goto 10
         end if
         msg = ' max lp iter'
         ifail = 4
         goto 15
      elseif (iprobstat.eq.6) then
         WRITE(nout,*)' Not enough workspace'
         msg = ' out of memory'
         ifail = 6
         goto 15
      elseif (iprobstat.ne.0) then
         WRITE(nout,*) ' Unexpected ifail from LP: ',iprobstat
         WRITE(nout,*) 
     &        ' Number of primal/dual infeasibilities and values:'
         WRITE(nout,*)'  primal:',inumpinf, rsumpinf
         WRITE(nout,*)'  dual  :',inumdinf, rsumdinf
         if ((inumpinf.eq.0).or.(rsumpinf.le.50.*rtolpinf)) then
            if ((rsumdinf.le.50.*rtoldinf).or.(nr_rsv.eq.2)) then
c           ... boldly go...
               goto 10
            end if
         end if
         if (nr_rsv.eq.0) then
            WRITE(nout,*) '1. resolve - coldstart'
            nr_rsv = 1
            set_stat = .false.
            goto 5
         end if
         if (nr_rsv.eq.1) then
            WRITE(nout,*)' 2. resolve - disturbing point'
            nr_rsv = 2
            set_stat = .false.
            do i=1,n
               x(i) = x(i) + 1.d-10*ems_drand(0)
            end do
            goto 5
         end if
         WRITE(nout,*) ' Maximal number of unsuccesful resolves!'
         WRITE(nout,*) ' STOP!'
         goto 99
      end if

c     ... LP solved succesfully (sort of)

 10   continue

      f_lp_af_phI = .false.
      if (scnd_try) then
c        ... reset tolerance to original value
         WRITE(nout,*)'reset tolerance to ',10.*csr_eps

         RTOLPINF = csr_eps*10.
         call ems_rset(rtcod, dspace, emsolr, emsolrln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_rset: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('RSET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if
         scnd_try = .false.
      end if

      if (se_mx_lp_it) then
         call ems_iget(rtcod, dspace, emsoli, emsoliln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_iget: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('IGET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if
         IMAXITER = 2*IITERNUM
         se_mx_lp_it = .false.
         call ems_iset(rtcod, dspace, emsoli, emsoliln)
         if (rtcod.ne.0) then
            WRITE(nout,*)'error in ems_iset: ',rtcod
            IF(RTCOD.GE.BAD_EMSERROR) THEN
               CALL GET_EMSERROR('ISET',RTCOD,EMSOL_ERRORMESSAGE)
            END IF
            goto 99
         end if
      end if

      nr_rsv = 0

      normd = 0.d0
      do i=1,n
         d(i) = dspace(ncolsol+i-1)
         if (abs(d(i))-s(i)*rho.gt.1e-5) then
            WRITE(nout,*) i, dclo(i), d(i), dcup(i)
         end if
         lam(i) = dspace(ncolrcosts+i-1)
         xnew(i) = min(max(blo(i), x(i) + d(i)), bup(i))
         normd = max(normd,abs(d(i)/s(i)))
      end do
      do i=1,m
         lam(n+i) = dspace(nrowduals+i-1)
      end do

      q = robjvalue
      oldrobjvalue = robjvalue
      iter = iter + 1

c      do i=1,n
c         write(6,*) i, dspace(ncollower+i-1), d(i),
c     &        dspace(ncolupper+i-1)
c      end do
c      write(6,*) robjvalue

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
     &   WRITE(nout,*) 'maximal multiplier on trust region: ',maxmu

c     ---------------- calculate KT residual ----------------------

c      pjp = la(0)
c      d = a(1:n)
c      do j=1,m
c         do i=la(pjp+j),la(pjp+j+1)-1
c            d(la(i)) = d(la(i)) - a(i)*lam(n+j)
c         end do
c      end do
c      kt_rsdu = 0.d0
c      do i=1,n
c         d(i) = d(i) - lam(i)
c         kt_rsdu = max(kt_rsdu,abs(d(i)))
c      end do

c      WRITE(nout,*)'LP KT res = ',kt_rsdu

      if (pr_to_fi) then
         open(26,file='lastiter.dat')
         write(26,*) rho
         do i=1,n
            write(26,*) xnew(i)
         end do
         close(26)
      end if

      call objfun(xnew, n, fnew, user, iuser, flag)
      call confun(xnew, n, m, cnew, a, la, user, iuser, flag)
      call gradient(n, m, mxa, xnew, a, la, maxa, user, iuser, flag)

      hcnew = 0.
      do i=1,m
         hcnew = hcnew +
     &       max(cnew(i)-bup(n+i),blo(n+i)-cnew(i),0.d0)/s(n+i)
      end do

c     dChange
c     see if best point yet
      call checkBestPoint(n,m, csr_eps, fnew, hcnew, xnew, 
     &     ispace(ncolstat:ncolstat+n-1),ispace(nrowstat:nrowstat+m-1)
     &     ,Bf,Bhc,Bx,sv_bs_c,Bsv_bs_r)
c     end dChange

c     ---------------- calculate KT residual ----------------------

c     The KKT-residual, for the problem 
c              min f(x) s.t. cl<= c(x) <= cu, bl <= x <= bu
c     is
c              Df(x) + A'y - lam    
c     where A = Dc(x), and y, lam are multipliers (in the last LP) on
c     constraint and variable bounds respectively
c     FIXME: is Danny right about the sign of +A'y???
c
c     In order to get the KKT-residual of the underlying NLP, 
c     multipliers lam should only be used if they refer to proper (not 
c     TR) bounds.

c     dChange
c     two changes to calculating kt resisual
c     1) correct error: change -a(i)*lam(n+j) to +a(i)*lam(n+j)
c     2) only subtract column dual if that column is at a real bound
c     and not if at a trust region bound
c     makes convergence test harder to satisfy

      pjp = la(0)
      do i=1,n
         d(i) = a(i)
      end do
      do j=1,m
         do i=la(pjp+j),la(pjp+j+1)-1
            d(la(i)) = d(la(i)) + a(i)*lam(n+j)
       !      write(nout,'(a,i5,a,i5,a,g12.5,a,g12.5)') 'var',
       !&           la(i), ' con', j, ' cgrad', a(i), ' row dual', lam(n+j)
         end do
      end do
      kt_rsdu = 0.d0

cAGR  the old KT residual (including TR bounds) is still computed for 
c     comp

      oldKTmeasure = 0 
      do i=1,n
         oldKTmeasure = max(oldKTmeasure,abs( d(i) - lam(i)))
      end do

      do i=1,n            
         if ( (abs(xnew(i)-blo(i))>=1d-12) .and.
     &        (abs(xnew(i)-bup(i))>=1d-12) .and. 
     &        ((abs(dspace(ncolsol+i-1)-dclo(i))<=1d-12) .or.
     &        (abs(dspace(ncolsol+i-1)-dcup(i))<=1d-12)) ) then
c     at trust bound not actual bound so do nothing
      !write(nout,'(a,i5,a,g12.5,a,g12.5,a,g12.5,a,g12.5,a,g12.5,a,g12.5)') 
      !&              'var', i, ' obj', a(i), ' vgrad', lam(i), 
      !&              ' val', xnew(i), 
      !&              ' old kt', d(i)-lam(i), ' new kt*', d(i)
         else
c     not at trust regino bound so here remove col dual as usual
      ! write(nout,'(a,i5,a,g12.5,a,g12.5,a,g12.5,a,g12.5,a,g12.5,a,g12.5)') 
      !&              'var', i, ' obj', a(i), ' vgrad', lam(i), 
      !&              ' val', xnew(i), 
      !&              ' old kt', d(i)-lam(i), ' new kt', d(i)-lam(i)
            d(i) = d(i) - lam(i) 
         end if
         kt_rsdu = max(kt_rsdu,abs(d(i)))
      end do
c     end dChange

      !if (kt_rsdu .ne. oldKTmeasure) 
      !&     WRITE(nout,*)'SLP KT res = ',kt_rsdu, ' old kt', oldKTmeasure



c     ---------------- check for convergence ----------------------

c     dChange
c     get kt_tol for use with new optimality condition
c     make it one times  biggest cost times csr_eps
      kt_tol =  maxCost*csr_eps
c     end dChange


c     dChange 
c     see if optimal by whether or not optimalityConditionMet
c     actual test same except use kt_tol instead of csr_eps
c     this typically laxer
      optimalityConditionMet=.false.
      if ((((normd.le.stp_eps).or.(kt_rsdu.le.kt_tol)).and.
     &     (hcnew.le.csr_eps*10.))) optimalityConditionMet=.true. 
c     end dChange

c      print *, optimalityConditionMet, do_rng, rlx_bnds
      if (optimalityConditionMet.and.do_rng.and.(.not.rlx_bnds)) then
         optimalityConditionMet = .false.
         rlx_bnds = .true.
         f_rng_it = .true.
         sv_obj = fnew
         print *,'Ranging: relax bounds and carry on'
         do i=1,n
            blo(i) = oldbnd(i)
            bup(i) = oldbnd(n+m+i)
         end do
         do i=1,m
            blo(n+i) = oldbnd(n+i)
            bup(n+i) = oldbnd(n+n+m+i)
         end do
            
      end if
      if (optimalityConditionMet)   then
         if (rlx_bnds) then
            if (abs(sv_obj-fnew)/fnew.gt.1d-8) then
               write(nout, *) 'WARNNING!'
               write(nout, *) 'Change of objective after bounds relaxed'
               write(nout, '(2G16.8)') sv_obj, fnew
               write(nout, *) 'This should not happen!'
               stop
            end if
         end if
         if (iprint.ge.1) then
            WRITE(nout,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '+22'

            CALL LOGOUT(2,dspace,bigiter,
     &                        iter,normd,rho,fnew,hcnew,nr_lp,+22)
         end if
         
         IF (msg_wrt_sumdat) THEN

            write(22,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '+22'

            CALL LOGOUT(2,dspace,bigiter,
     &                        iter,normd,rho,fnew,hcnew,nr_lp,+22)
            
         END IF

         if (iprint.ge.2) WRITE(nout,*) ' convergence!'

         if (do_rng) then
            call do_range(dspace, ispace, nspace, n, m, status, dcup,
     &           drup, dclo, drlo, dobj, result, rifail, oldbnd, s,
     &           rho, x, c, ngrad, mrow, mcol, dels, 
     &           sv_bs_r, sv_bs_c, iprint, oldrobjvalue, sv_t_fi)

c           call ranging(dspace, ispace, nspace, n, m, status, dcup,
c     &          drup, dclo, drlo, dobj, result, rifail)
            if (rifail.eq.1) then
               if (iprint.ge.1) WRITE(nout,*)' ranging failed'
               msg = msg(1:trim77(msg))//' range fail'
c              ifail = 3
            end if
         end if


c     ... set x = xnew, c = cnew
         f = fnew
         hc = hcnew
         do i=1,n
            x(i) = xnew(i)
         end do
         do i=1,m
            c(i) = cnew(i)
         end do


         goto 15
      end if
      if (iter.ge.max_iter) then
         if (iprint.ge.2) WRITE(nout,*)' maximal number of iterations'
         msg = msg(1:trim77(msg))//' max iter'
         ifail = 2
         goto 15
      end if

c     dChange
c     change to filter
c     consider a point dominated if both a) and b) hold
c     a) some other point has lower feas, or is feas 
c     b) that same other point has lower obj
      accept = .true.
      if (1==1) then
         do i = 1,sizef
cAGR        I DONT BELIEVE IN THE CHANGE OF THE FILTER!
cAGR        The first change is sensible (once a point is feasible within the
cAGR        tolerance, only objective function value should be compared)
cAGR        but WHY reject if not at least 0.1 better in objective? This
cAGR        sounds arbitrary 
            if ( ( (hcnew>=fc(i)) .or.(10*csr_eps>=fc(i)) )
     &           .and.(fnew+1e-01>=ff(i)) )  accept = .false.
         end do
      else 
         do i=1,sizef
            if ((hcnew.ge.fc(i)).and.(fnew.ge.ff(i))) accept = .false.
         end do 
      end if

c     end dChange
      if (f_rng_it) then
c        ... if we carry on with relaxed bounds after convergence
c            then always accept first step
         accept = .true.
         f_rng_it = .false.
      end if

c     dChange 
c     a few changes to ratio test
c     corrected sign error in ratiof 
c     added divided by 0 condition for ratioc
c     dont use ratioc if new point feasible

      if (q.lt.0) then
         ratiof = (fnew - f)/q
      else
c        to prevent division by zero
         ratiof = 0.0
      end if

      if (hc==0) then 
c        again to prevent division by zero
         ratioc = 0.0
      else
         ratioc = -(hcnew - hc)/hc
      end if
c     Once new point feasible witin tolerance, set ratioc = 1.0
      if (hcnew<=10*csr_eps) ratioc = 1.0
c     end dChange

      ratio = min(ratiof, ratioc)
         
      if (iprint.ge.2) then
         WRITE(nout,*) 'reductions in model & reality + ratio:'
         WRITE(nout,'(A,3G20.10)') ' f: ',q, fnew - f, ratiof
         WRITE(nout,'(A,3G20.10)') 'hc: ',-hc, hcnew - hc, ratioc

         WRITE(nout,'(A,G20.10)') ' => ratio = ',ratio
      end if

      if (accept) then
         if (iprint.ge.2) WRITE(nout,*) 'step acceptable'
         if (iprint.ge.1) then
            WRITE(nout,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '+22'
         end if

         IF (msg_wrt_sumdat) THEN

            write(22,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '+22'

         END IF

         do i=1,sizef
            if ((hcnew.le.fc(i)).and.(fnew.le.ff(i))) then
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
            goto 15
         end if
         sizef = sizef + 1
         fc(sizef) = hcnew
         ff(sizef) = fnew
         f = fnew
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
            WRITE(nout,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '-22'
            
            CALL LOGOUT(2,dspace,bigiter,
     &                        iter,normd,rho,fnew,hcnew,nr_lp,-22)
         end if

         IF (msg_wrt_sumdat) THEN

            write(22,'(I3,2G16.8,2G14.8,I5,A5)')
     &           iter, normd, rho, fnew, hcnew,nr_lp, '-22'

         END IF

         rho = min(rho, normd)/4.
         if (iprint.ge.2)
     &     WRITE(nout,*) 'resolve with smaller trust region: ',rho

      end if

c     dChange
c     see if point unbounded, but shouldn't ever be
c     unbounded if a bound is +- infty and a var is at that bound
      if (1==1) then
         do i = 1,n
            if ( (blo(i)==-infty .and. (abs(x(i)-blo(i))<=1)) .or.
     &           (bup(i)== infty .and. (abs(x(i)-bup(i))<=1)) ) then
               msg = ' unbounded '
               ifail = 8
               goto 15
            end if
         end do
      end if
c     end dChange

      goto 5

 99   continue
c     computational error
      ifail = 3
      msg = msg(1:trim77(msg))//' slp error'
 15   continue
 
      write(ch_fnew,'(G16.8)') fnew
 17   continue

      IF (msg_wrt_sumdat) THEN

         inquire(file=SUMMARY_FILE, exist=fi_xst)
         if (fi_xst) then
c     open(23, file=SUMMARY_FILE, access='APPEND')
            open(23, file=SUMMARY_FILE, position='APPEND')
         else
            open(23, file=SUMMARY_FILE)
            write(23,'(a)') ' Values at the end of each slp run'
            write(23,'(3a)') ' nr  iter    step           trust',
     &           '         objective       csr_error',
     &           '        seed'
         end if

c     dChange
c     at this stage of reporting and recording solution
c     could potentially change so recording 
c     Bf and Bhc instead of fnew and hcnew 
c     end dChange

         write(23,'(2I3,2G16.8,a16,G16.8,i11,i10,a)')
     &        bigiter, iter, normd, rho, ch_fnew, hcnew,
     &        seed, tt1-tt00, msg(1:trim77(msg))
         close(23)
         
      END IF
      tt1 = ems_time()

c      IF(IFAIL.EQ.0.OR.IFAIL.EQ.1) THEN
c         NUMBER_SOLUTIONS=NUMBER_SOLUTIONS+1
c         WRITE(20,'(3(I3,A),4(G16.8,A))') NINF,',',
c     &                                           bigiter,',',
c     &                                           iter,',',
c     &                                           normd,',',
c     &                                           rho,',',
c     &                                           fnew,',',
c     &                                           hcnew,','
c      END IF

      end subroutine slp_main


c     **************************************************************
c     *
c     * checkBestPoint
c     *
c     **************************************************************

        subroutine checkBestPoint(n,m, csr_eps, fTest, hcTest, xTest, 
     &       sv_bs_cTest, sv_bs_rTest,
     &       Bf,Bhc,Bx,Bsv_bs_c,Bsv_bs_r)

c     n, m     problem dimensions
c     csr_eps  feasibility tolerance (hc<10*csr_eps is feasible)
c        the new point that may be the best one yet:
c     fTest, hcTest, xTest, sv_bs_c/rTest
c        the best point so far
c     Bf, Bhc, Bx, Bsv_bs_c/r

        implicit none
c     general problem variables
        integer n,m
        double precision csr_eps
c     best storage
        double precision Bf
        double precision Bhc
        double precision Bx(n) 
        integer  Bsv_bs_c(n), Bsv_bs_r(m)
c     input vars
        double precision fTest, hcTest, xTest(n)
        integer  sv_bs_cTest(n), sv_bs_rTest(m)
c     local vars
        logical bestPoint
        integer i
        integer    bt30_32
        parameter (bt30_32 = 2**31+2**30+2**29)
      
c     best point is most feas if no feas point found else
c     best obj among feasible ponits

      bestPoint=.false.
c     if we already have a feasible point, then make sure new point is 
c     feasible and compare objective values
      if ( Bhc<=10*csr_eps) then 
         if ((hcTest<=10*csr_eps).and. fTest<Bf) then
            bestPoint=.true.
         end if
c     otherwise just compare feasibilities
      else
         if (hcTest<Bhc ) then
            bestPoint=.true.    
         end if
      end if
      if (bestPoint) then
         !print*, 'old best', Bf, Bhc
         !print*, 'new best', fTest, hcTest
c     save new best point
         Bf = fTest
         Bhc = hcTest
         Bx = xTest
         do i = 1,n
            Bsv_bs_c(i) = iand(sv_bs_cTest(i), bt30_32)
         end do
         do i = 1,m
            Bsv_bs_r(i) = iand(sv_bs_rTest(i), bt30_32)
         end do
      end if
         

      end subroutine checkBestPoint
c     end dChange

c     retruns codes
c
c     ifail = 0  - no problems
c     ifail = 1  - convergence to infeasible point
c     ifail = 2  - max number of slp iterations
c     ifail = 3  - computational error
c     ifail = 4  - max number of lp iterations
c     ifail = 5  - max number of filter entries
c     ifail = 6  - not enough workspace
c     ifail = 7  - inconsistent bounds
c     dChange
c     ifail = 8  - unbounded
c     end dChange
