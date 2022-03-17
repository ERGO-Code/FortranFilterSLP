c     SLP driver: interface routines to the integra2 SLP solver:
c     -----------------------------------------------------------
c     
c     There are two routines: slp_driver and get_dimension
c     get_dimension should be called first:
c     
c     get_dimension(filename)
c     
c     where 'filename' is the path to a sonoco.def file which 
c     holds information about the *.WRK files defining the problem
c     
c     The information read in by get_dimensions is stored in 
c     the module interface_mod (which needs to be USE'd). This 
c     provides access to the problem size parameters:
c     
c     nu, rm, mi, de: basic problem dimensions
c     ntgdim        : size of the ntgresults array
c
c     The calling method then needs to allocate ntgresults and the 
c     real, integer workspace in ws, lws and call
c
c     slp_driver(ws, sz_ws, lws, sz_lws, ntgresults, nresults, ifail)
c
c     which will solve the SLP problem and return information on
c     the best solution found in ntgresults.
c
c     The behaviour of slp_driver can be changed by setting some
c     control variables which are defined in COMMONs:
c
c     in SLPCOMM.INC:
c      - PROB_GLOBAL: The required probability of finding the 
c                     global solution (default: 95)
c      - NUMBER_STARTS: maximal number of restarts (default: 50)
c                    (less may be done depending on the setting
c                     of PROB_GLOBAL, but NUMBER_STARTS is an upper
c                     bound)
c     in msg.inc:
c      - Various control parameters (msg_xxxx) to switch debugging
c        ouput on (all set to .false. by default)
c     
c     Information on the solution status is returned in ifail
c      - ifail = 0: ok
c              = 1: Problem infeasible (in all runs)
c              = 20: Insufficient Memory passed
c              = 21: Insufficient Memory (dynamic allocation)
c              = 30: Some hardwired limit exceeded
c              = 50: Wrong usage (some data file is inconsistent)
c              = 99: Solver internal error
c     
c     More information on the solution status/error messages can be
c     obtained by including 'return.inc' and inspecting the variables
c      - global_ifail (which should equal ifail above), 
c      - global_err_msg,
c      - global_int_arg/global_double_arg
c
c


c     -----------------------------------------------------------------
c     SLP_DRIVER
c     -----------------------------------------------------------------
      subroutine slp_driver(ws,sz_ws,lws,sz_lws, ntgresult, nresults,
     &     ifail)
c     
c     This is the main interface routine into SLP. 
c     
c     PARAMETERS(IN)
cD      ws(sz_ws)      double workspace 
cI      lws(sz_lws)    integer workspace
cI      nresults
c     PARAMATERS(OUT)
cD      ntgresults(*)
c          Array to return results to Integra-GUI. Carries information
c          about Integras variables. 
c          The array is sliced up into blocks of nresult(=8) values,
c          their meaning is as follows
c            1 - value
c            2 - lower bound
c            3 - upper bound
c            4 - sensitivity
c            5 - ranging 1
c            6 - ranging 2
c            7 - ranging 3
c            8 - status (1=free, 2=on min,3=on max,4= <min, 5= >max)
c          The Integra variables are ordered as follows
c            [(MI+DE)*(NU+RM+MI)]
c                         nutrient/rawmat/mixer content of mixer/demand
c            [RM+MI]         total usage
c            [NURATM/P]      nutrient ratios
c            [RMINC]         raw mat inclusions (GIGs)
c            [SPECIAL]       special constraints
c     OTHER
c     

      use interface_mod
      use types_module
      implicit none

        INCLUDE 'SLPCOMM.INC'
        include 'pusr.inc'
        include 'msg.inc'
        include 'return.inc'
cHM
        include 'ntgres.inc'

      integer, parameter :: MAX_CLASSES = 50

      integer sz_ws, sz_lws, nresults, ifail

CJCH *** parameter (sz_ws = 1000000, sz_lws = 100000)

      integer lws(sz_lws)
      double precision ws(sz_ws), ntgresult(*)


cHM

      integer n, m, max_f, m1, sn

      integer iprint, rifail, maxiter, error, iflag


      double precision iz_rho, f, hc, rho, mx_iz_x

      integer istat(13)

      double precision rstat(3)

      character*4 int2char
      integer trim77
      integer ems_time
      type(special_da):: SPD
      type(integra_da):: DAT

      double precision, allocatable :: oldbnd(:)

      double precision ems_drand, dummy
      logical rd_rnd_sd, fi_xst, sv_t_fi, se_mx_lp_it, wr_msg, do_rng
      integer sz_user, sz_iuser, max_cs_terms
      integer sv_t_fi_int, se_mx_lp_it_int, wr_msg_int

      integer p_blo, p_bup, p_x, p_c, p_s, p_lam, p_user, p_wws,
     &     p_bns, n_bns, n_lws, sz_wws, sz_llws, p_llws, p_nx
      integer p_iuser, p_status, p_result, p_bdscl, p_vrreg

      double precision     infty, eps
      common /NLP_eps_inf/ infty, eps

      double precision stp_eps, csr_eps
      common /tol/     stp_eps, csr_eps

      character ch
      logical file_exists
      integer i, j, mx_big_it, slp_ws, nspace, mxwk, mxiwk, err
      double precision totmu

      integer        seed, bigiter
      common /summc/ seed, bigiter

c     --- variables for the stopping logic

      double precision BESTCOST, DTEMP, val
      integer BESTCOUNT, BESTSEED, BESTIFAIL 

      integer k, p
      logical fd, resolve
      integer :: kobj = 0                !number different objs seen
      double precision vobj(MAX_CLASSES) !values of different objectives
      integer iobj(MAX_CLASSES)          !ifail for different sols
      integer nobj(MAX_CLASSES)          !number of times sol i is seen
      integer :: perclevel(4) = (/ 90, 95, 98, 99/)
      integer :: nbar(9,4) = reshape(
     &  (/   4,   7,  11,  15,  19,  24,  28,  33,  38,
     &       5,   9,  13,  18,  23,  28,  33,  39,  50,
     &       6,  11,  17,  22,  28,  34,  40,  50,  50,
     &       7,  13,  19,  25,  32,  38,  46,  50,  50/), (/9, 4/))


c     ========================  procedure body  ======================
c     - read in problem dimensions & allocate space -

c     ... set p to the correct level for the required probability
c         of having found the global optimum
        p = 1;
        do i=2,4
           if (prob_global.ge.perclevel(i)) p = i
        end do


c     dChange
c     some problems eg sl27Jun solve an awful lot faster if we put
c     mx_iz_x=1 to force a start point near zero
c     so now this option read in by SLPCT_FILE 
                                !mx_iz_x = 10000.d0
c     end dChange
     
      global_ifail = -1
       

c     flag to make sure that the GLOINDEX array is initialised once
c     (at the start of the 50 runs?)
      GLO_INDEX_FLAG = 0
c     flag to make sure that the RINDEX array is initialised once
      REL_INDEX_FLAG = 0
c     flag to make sure that the RINDEX array is initialised once
      RATIO_NUTR_INDEX_FLAG = 0

c     ---------------- read control variable file -------------------

      inquire(file=SLPCT_FILE,exist=file_exists)
      if (.not.file_exists) then
         maxiter = 500

CJCH         iz_rho = 5000.
         iz_rho = 5.                 ! Was this in previous version

         eps   = 1d-6   !  set tolerances and infinity here
         infty = 1d31   !    depending on problem dimension
         stp_eps = 2.d-4     ! originally stp_eps = 1.e-5
         csr_eps = 1.d-6     ! originally csr_eps = 1.e-7
         iprint = 1     ! Input iprint for LP--subproblems
         nout = 6           ! output channel (>= 6) (6 = screen)
         nspace = 5000000  ! EMSOL workspace
         max_f  = 100
         se_mx_lp_it = .false.
         sv_t_fi = .false.
         wr_msg = .false.
         mx_iz_x = 10000
      else
         open(31,file=SLPCT_FILE)
         read(31,*) maxiter
         read(31,*) iz_rho
         read(31,*) eps
         read(31,*) infty
         read(31,*) stp_eps
         read(31,*) csr_eps
         read(31,*) iprint
         read(31,*) nout
         read(31,*) nspace
         read(31,*) max_f
         read(31,*) se_mx_lp_it_int
         read(31,*) sv_t_fi_int
         read(31,*) wr_msg_int
         read(31,*) mx_iz_x
         close(31)
         if (sv_t_fi_int.eq.1) then
            sv_t_fi = .true.
         else
            sv_t_fi = .false.
         end if
         if (se_mx_lp_it_int.eq.1) then
            se_mx_lp_it = .true.
         else
            se_mx_lp_it = .false.
         end if
         if (wr_msg_int.eq.1) then
            wr_msg = .true.
         else
            wr_msg = .false.
         end if
      end if
      if (msg_log) then
         write(nout,*) ' Maxiter       = ', maxiter
         write(nout,*) ' iz_rho        = ', iz_rho
         write(nout,*) ' eps           = ', eps
         write(nout,*) ' infty         = ', infty
         write(nout,*) ' stp_eps       = ', stp_eps
         write(nout,*) ' csr_eps       = ', csr_eps
         write(nout,*) ' iprint        = ', iprint
         write(nout,*) ' nout          = ', nout
         write(nout,*) ' nspace        = ', nspace
         write(nout,*) '  WARNING request for nspace can not be obeyed'
         write(nout,*) ' max_f         = ', max_f
         write(nout,*) ' se_mx_lp_it   = ', se_mx_lp_it
         write(nout,*) ' sv_t_fi       = ', sv_t_fi
         write(nout,*) ' wr_msg        = ', wr_msg
         write(nout,*) ' mx_iz_x       = ', mx_iz_x
      end if

c     -------------------- get problem dimensions --------------------


      p_bns = 1
      p_llws = sz_lws/2
      
      n_bns = p_llws - p_bns
      n_lws = sz_lws - p_llws

      call rd_prob_dim(sn, lws(p_bns), n_bns, lws(p_llws), 
     &     n_lws, ws, SPD, DAT, iflag)
      if (iflag.gt.0) goto 99

      if (msg_log) then
         WRITE(nout,*) ' number of nutrients    : ',nu
         WRITE(nout,*) ' number of raw materials: ',rm
         WRITE(nout,*) ' number of demands      : ',de
         WRITE(nout,*) ' number of premixes     : ',mi
      end if

      p_nx = p_bns + 4*(mi+de) + 2*n_mu + 3*(rm+mi)  

c      WRITE(nout,*) ' number of rmava c/s    : ',n_u
c      WRITE(nout,*) ' number of rminc c/s    : ',n_nu
c      WRITE(nout,*) ' number of nurat c/s    : ',n_nuratp+n_nuratm

c     ------------- set memory-map for rd_prob_da ----------------

      n = n_vr
      m = n_cs
      m1 = m+1

      p_bdscl  = 1
      p_blo    = 1
      p_bup    = p_blo    + n+m1
      p_s      = p_bup    + n+m1
      p_user   = p_bdscl  + 3*(n+m1)
      sz_wws   = 2+3*nu + 2*(mi+max(de,rm))
     &     + max(2*nu, 2*(rm+mi), 2*n_rminc)+1
      p_wws    = sz_ws - sz_wws
      sz_user  = p_wws - p_user

      p_vrreg  = p_nx
      p_iuser  = p_vrreg  + 4*n + 2*sn + 2
      sz_llws  = mi+max(de,rm) +
     &     max(2*n_rminc+rm, n_nuratp, n_nuratm)+1
      p_llws   = sz_lws - sz_llws
      sz_iuser = p_llws - p_iuser

      if (msg_log) then
         write(nout,*) 
     &        'real/int workspace for rd_prob_da: ',sz_user, sz_iuser
      end if

      max_cs_terms = max(max(max(40, rm), de), nu)

      allocate(oldbnd(2*m+2*n), stat=err)
      if (err.ne.0) then
         global_ifail = 20
         write(global_err_msg, '(A,I5,A)') 
     &        'Failed allocation of ',2*m+2*n,' bytes.'
         global_int_arg(1) = 2*m+2*n
         goto 99
      end if

 10   continue
      call rd_prob_da(n, m, sn, ws(p_user), lws(p_iuser), sz_user,
     &     sz_iuser, lws(p_vrreg), ws(p_bdscl), oldbnd, ws(p_wws), 
     &     sz_wws, lws(p_llws), sz_llws, lws(p_bns), n_bns, 
     &     max_cs_terms, SPD, DAT, iflag, NTGRESULT,NRESULTS)
    

      if (iflag.eq.3) then
c     ... max_cs_terms too low
         max_cs_terms = 10*max_cs_terms
         goto 10
      end if
      call free_special(SPD)
      if (iflag.eq.1) then
         if (msg_err) then
            write(nout,*)
     &           'FATAL ERROR: not enough real or int memory in reg_cs'
         end if
c     ... global_ifail/global_err_msg is set in reg_cs
         goto 99
      end if
      if (ifail.eq.2) then
         global_ifail = 1
cHM         ...Set i rd_prob_da.f for ifail = 2
         if (msg_err) then
            write(nout,'(a)') 'Infeasibility detected in presolve'
         end if
         goto 99
      end if
      if (ifail.ne.0) then
         if (msg_err) then
            write(nout,*) 'FATAL ERROR: read problem data failed'
         end if
         global_err_msg = 'ERROR: read problem data failed.'
         global_ifail = 99
         goto 99
      end if



      p_x      = p_user + nx_dnty + 1
      p_c      = p_x    + n
      p_lam    = p_c    + m
      p_result = p_lam  + n+m
      p_wws    = p_result + 7*(n+m)

      p_status = p_iuser + nx_nty + 1
      p_llws    = p_status + n+m



      mxwk = sz_ws - p_wws + 1
      mxiwk = sz_lws - p_llws + 1
      if (msg_log) then
         write(nout,*)
         write(nout,*) ' double workspace for SLP: ',mxwk
         write(nout,*) 'integer workspace for SLP: ',mxiwk
      end if
      if (mxwk.le.0.or.mxiwk.le.0) then
         if (msg_err) then
            write(nout,*) 'FATAL ERROR: not enough workspace'
            WRITE(nout,*)mxwk, mxiwk
         end if
         if (mxwk.le.0) then
            write(global_err_msg, '(A,I8,A,I8)') 
     &           'ERROR: not enough real workspace to call SLP, need ',
     &           p_wws, ' ,got ',sz_ws
         end if
         if (mxiwk.le.0) then
cHM
            write(global_err_msg, '(A)') 
     &           'ERROR: not enough int workspace to call SLP, need  ',
     &           'xxxxxxxxxx got xxxxxxxxxx.'
            global_int_arg(1) = p_llws
            global_int_arg(2) = sz_lws
         end if
         global_ifail = 21
c            write(global_err_msg, '(A,I,A,I)') 
c     &           'ERROR: not enough int workspace to call SLP, need ',
c     &           p_llws, ' ,got ',sz_lws
c         end if
c         global_ifail = 3
         goto 99
      end if

        CALL LOGSET(m,n,MAX(mxwk/131072,1))

      inquire(file=RNDSD_FILE,EXIST= fi_xst)
      if (fi_xst) then

         open(25,file=RNDSD_FILE)
         read(25,*) seed
         close(25)
      else
         seed = ems_time()
c         seed = mod(seed,30081)
      end if


c     Fixed seed is better than random for debugging, but not for production
c     code. If a fixed seed is desired this should be done through the 
c     rnd_sd.dat file
c      seed = 36966
c      seed = 30856
c      seed = 30862
c      seed = 1359538442
c       seed = 1366723166

      OPEN(27,FILE='COST_PER_STEP.TXT')
c      write(27,*)' Cost Repeat %          : ',COSTREPETE
c      write(27,*)' Minimum Steps          : ',MINSTEPS
c      write(27,*)' Minimum Cheapest Costs : ',MINCHEAPESTCOSTS
      write(27,*)' '
      write(27,*)' Step    Seed            Cost'//
     &     '           Cheapest Cost    Count  Inf.   %'
      write(27,*)' ----- --------- -------------------'//
     &     ' ------------------- ----- ----- -----'
      CLOSE(27)      
cHM end ----------------------------------- 


c     ===============================================================
c       Start multiple random start iterations
c     ===============================================================

      bigiter = 0
      BESTCOST      = 1.d+20
      BESTCOUNT = 0
      resolve = .false.
      do_rng = .false.

 5    continue
c
chm Append to log file
      OPEN(27,FILE='COST_PER_STEP.TXT',POSITION='APPEND')

      bigiter = bigiter + 1

      if (bigiter.gt.1.and.(.not.resolve)) then
         seed = seed +1
      end if

c     seed = ems_time()
c     seed = mod(seed,30081)
      if (msg_log) then
         write(nout,*) 'Initial seed = ', seed
      end if
      dummy = ems_drand(seed)

      do i=1,n+m
         ws(p_lam+i-1) = 0.d0
      end do

c     ... get random number
      do i=1,n
         ws(p_wws+i-1) = ems_drand(0)
      end do

c     ... set random starting point
      do i=1,n
         ws(p_x+i-1) = max(ws(p_blo+i-1),-mx_iz_x) 
     &        + ws(p_wws+i-1)*(min(ws(p_bup+i-1),mx_iz_x)
     &        -max(ws(p_blo+i-1),-mx_iz_x))
      end do

      rho = iz_rho

c      open(21,file='start.dat')
c      do i=1,n
c         read(21,*) ws(p_x+i-1)
c      end do
c      close(21)

      IF (msg_wrt_sumdat) THEN
c     open(unit=22,file=dirname(1:trim77(dirname))//'/sum'//
c     &                 int2char(bigiter)//'.dat')
         open(22,file='sum'//
     &        int2char(bigiter)//'.dat')
         
         write(22,*) seed
      END IF

c     ... call the main SLP routine
      call filterSLP(n, m, n_maxa, max_f, mxwk, mxiwk,
     &            iprint, iflag, rifail, rho, ws(p_x), ws(p_c), f, hc,
     &            ws(p_blo), ws(p_bup), ws(p_s), ws(p_wws), lws(p_llws),
     &            ws(p_lam), ws(p_user), lws(p_iuser), maxiter,istat,
     &            rstat, sv_t_fi, se_mx_lp_it, wr_msg, lws(p_status),
     &            ws(p_result), oldbnd, do_rng)

      if(iflag.eq.0.or.iflag.eq.1) then
         if (msg_solprt) 
     &      call sol_print(n, sn, ws(p_x), ws(p_blo), ws(p_bup), f, 
     &        lws(p_vrreg), bigiter, lws(p_llws), mxiwk,
     &        ws(p_wws), mxwk)
         if (msg_lamprt)
     &     call lam_print(n, sn, m, ws(p_x), lws(p_vrreg), ws(p_user), 
     &        lws(p_iuser), ws(p_blo), ws(p_bup), ws(p_lam), f, ws(p_c),
     &        mxiwk, lws(p_llws), lws(p_status))
         if (msg_ranprt.and.(rifail.eq.0))
     &        call ran_print(n, m, sn, ws(p_result), lws(p_vrreg), 
     &              lws(p_iuser))

c     write default error message: if everything is fine, the 
c     message from the final (best) run is reported

         if (iflag.eq.0) then
            write(global_err_msg,'(A)') 'Solution found OK'
            global_ifail = 0
         else
            write(global_err_msg,'(A)') 'Problem infeasible'
            global_ifail = 1
         end if
      else
c     ... all other errors should have been picked up before, this
c         should not happen
         if (msg_err) then
            WRITE(nout,'(A,I2,A)')
     &           'Could not solve SLP run ',bigiter," succesfully"
         end if
         write(global_err_msg, '(A,I3,A)')
     &        'Could not solve SLP run ',bigiter,' succesfully'
         global_ifail = 99
         
      end if

      IF (msg_wrt_sumdat) THEN
         close(22)
      END IF

c      open(21,file='bestdbl.dat')
c      do i=1,n
c         write(21,'(d25.16)') ws(p_x+i-1)
c      end do
c      close(21)
c      stop


c     ----------- Danny/AGR logic to decide when to stop --------------
      
c     ifail, f, hc are returned from filterSLP

c     The idea is to keep track of how many different solutions (k)  
c     we have seen in n runs and stop when the chance of having found
c     the best solution is >= p%

c     nbar(k, p) is the maximal number of iterations necessary if
c       k different classes have been seen so far.

c     need to keep track of different solutions seen and how often
c     they are obtained


cHM   ----------- HMs logic to find the best solution yet -------------

c     bestcost   = best objective found so far
c     bestseed   = seed leading to best objective 
c     bestcount  = number of times the best objective has been found
c     dtemp          = percentage that best sol is found
c
c


chm        write(*,*) seed
chm        pause

c     If this solve is feasible, check if this is better than the 
c     best cost found so far. If so remember value and seed.
c     counts the number of times the cheapest value so far has been
c     found.
      
      if (.not.resolve) then
         val = f
         if (iflag.eq.1) val = hc
         if (msg_multstart)
     &        print *, 'Last solution:', iflag, val
         
c     ... check if solution has been found before
         fd = .false.
         if (msg_multstart) then
            print *, 'Currently know ',kobj,' different solutions:'
            do k=1,kobj
               print *, k, iobj(k), vobj(k), nobj(k)
            end do
         end if
         do k=1,kobj
            if ((iflag.eq.iobj(k)).and.(dabs(f-vobj(k))).le.1.d-3) then
               if (msg_multstart)
     &          print *, ' found solution ',k,'again:', iobj(k), vobj(k)
               fd = .true.
               nobj(k) = nobj(k) + 1
               exit
            end if
         end do
         if (.not.fd.and.(kobj.lt.MAX_CLASSES)) then
            kobj = kobj +1
            iobj(kobj) = iflag
            vobj(kobj) = val
            nobj(kobj) = 1
         end if

         mx_big_it = NUMBER_STARTS
         if (kobj.le.9) then
            mx_big_it = nbar(kobj, p)
         end if
      
         if (msg_multstart) then
            print *, 'found ',kobj,' classes in ',bigiter,' trials'
            print *, 'should stop after ',mx_big_it,'iterations' 
            if (bigiter.ge.mx_big_it) then
               print *, 'Should stop now!'
            end if
         end if
      
c     if first iteration initialise 
         IF(BESTCOUNT.EQ.0) THEN
            BESTCOST=val
            BESTIFAIL = iflag
            BESTCOUNT=1
            BESTSEED=seed
         ELSE
c     found a better solution
            if ((iflag.lt.BESTIFAIL).or.
     &           (iflag.eq.BESTIFAIL.and.(f.lt.bestcost-1.d-6))) then
               BESTCOST=val
               BESTCOUNT=1
               BESTSEED=seed
               BESTIFAIL = iflag
            END IF
*     Found a repeated cheap cost     
            if (iflag.eq.BESTIFAIL.and.dabs(val-BESTCOST).lt.1.d-6) then
               BESTCOUNT=BESTCOUNT+1
            END IF
         END IF
      

c     BESTCOST     : best value found so far
c     BESTCOUNT: number of times it has been found
c     BESTSEED  : seed leading to the best value
         
*     
         DTEMP    =BESTCOUNT*100.0D0/bigiter
*     
         WRITE(27,'(I6,I10,2F20.6,3I6)') bigiter,seed,f,BESTCOST,
     &        BESTCOUNT,iflag, INT(DTEMP)
*     
      end if
      CLOSE(27)

c     If not on resolve iteration 
      if ((.not.resolve).and.(bigiter.lt.mx_big_it)) then
         goto 5
      end if
      
c     if done all the iterations, resolve the best one so far
      if (.not.resolve) then
         seed=BESTSEED
         resolve = .true.
         do_rng = .true.
         goto 5
      END IF

c     ===============================================================
c       End multiple random start iterations
c     ===============================================================

c      global_ifail = 0

 99   continue

      ifail = global_ifail
      if (msg_err) then
         write(nout,'(A)') global_err_msg
      end if

c
chc  Close the cost per step output file
c
      if((iflag.eq.0.or.iflag.eq.1).and.bigiter.ge.mx_big_it) then
c
chc  Get sensitivity data back into the NTGRESULT table     
c
         call GETSENSITIVITY(n, sn, m, ws(p_x),lws(p_vrreg),
     &        lws(p_iuser), ws(p_blo), ws(p_bup), ws(p_lam), ws(p_c),
     &        mxiwk, lws(p_llws), ws(p_result), lws(p_status),ifail,
     &        NTGRESULT,NRESULTS)
       end if
*

chm   ... deallocate memory for 'oldbnd'
      deallocate(oldbnd)

      end



c     =============================================================
c     get_dimensions
c     =============================================================
      
      subroutine get_dimensions(filename, ifail)
c      subroutine get_dimensions(filename, nu, rm, mi, de,nnurat,nrminc,
c     &     nspec, nnl)

c     This subroutines reads the data files defining the problem
c     and returns the problem dimensions. Can be used in order to
c     allocate work-arrays of appropriate sizes before the solver
c     proper is called
c     
c     PARAMETERS
cCI    filename         file giving names of the datafiles that define
c                         the problem
cIO    nu, rm, mi, de   
cIO    nnurat
cIO    nrminc
cIO    nspec
cIO    nnl

      use interface_mod
      use types_module

      implicit none

      include 'SLPCOMM.INC'
      include 'msg.inc'
      include 'return.inc'
      
      character(*) filename
      integer nnurat, nrminc, nspec, nnl, ifail

      integer iflag
      logical file_exists

      type(special_da):: SPD
      
      ifail = 0
c     locate the sonoco.def file and read the filenames
      call read_filenames(filename, ifail)
      if (ifail.ne.0) goto 99

      inquire(file=DIM_FILE, exist=file_exists)
      if (.not.file_exists) then
         WRITE(nout,*) 'Problem dimension file ',DIM_FILE,' not found.'
         write(global_err_msg, '(A,A,A)') 
     &        'Problem dimension file ',DIM_FILE,' not found.'
         global_ifail = 50
         ifail = 1
         goto 99
      end if
      if (msg_log) then
         WRITE(nout,*) ' reading problem dimensions:',DIM_FILE
      end if
      open(31,file=DIM_FILE)
      read(31,*) nu, rm, mi, de
      close(31)

      ntgdim = (MI+DE)*(NU+RM+MI) + (RM+MI) + 500
      
c      NUMBER_NUTRIENTS = nu
c      NUMBER_RAWMATERIALS = rm
c      NUMBER_SPECIFICATIONS = de
c      NUMBER_PREMIXES = mi

      call rd_prob_dim2(filename, SPD, iflag)

 99   continue

      end subroutine

c     =============================================================
c     read_filenames
c     =============================================================
c     
c     Given the input filename (which may include a path) this 
c     routine reads the pointed to file for the names of the
c     other problem definition files. 
c     If the initial filename contains a path, this path is
c     added to all the filenames.
c
c     INPUT: 
c         filename
c
c     EFFECTS: 
c       - The filenames in SLPCOMM.INC are set

      subroutine read_filenames(filename, ifail)

      implicit none
      
      character(*) filename
      integer ifail

      CHARACTER*80     LOGOUTFILE,
     &                 RESULTSFILE,
     &                 STRING
      integer iunit, pos
      logical file_exists
      character ch

      data iunit /1/

      include 'SLPCOMM.INC'
      include 'return.inc'

*-----------------------------------
*  Read the SONOCO parameter file  *
*-----------------------------------
*
      ifail = 0
      
      inquire(file=filename, exist=file_exists)
      if (.not.file_exists) then
         WRITE(nout,*) ' Problem definition file ',filename,' not found'
         write(global_err_msg, '(A,A,A)')
     &        ' Problem definition file ',filename,' not found'
         global_ifail = 50
         ifail = 1
         return
      end if

c     ... scan if filename contains a '/', if so remember dirname

      pos = index(filename, "/", .true.)
c      fd_slash = .false.
c      do i=trim77(filename),1, -1
c         if (filename(i:i)=='/') then
c            fd_slash = .true.
c            break
c         end if
c      end do
c      if (fd_slash) then
      if (pos>0) then
         DIRNAME = filename(1:pos)
      else
         DIRNAME = "./"
      end if

      OPEN(IUNIT,FILE=filename)
*
      READ(IUNIT,'(A)') STRING  ! Skip the straights cost penalty
*
      READ(IUNIT,'(A)') STRING ! Problem dimemsions file
      DIM_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Specification NT constraints
      LMTS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Premix RM constraints
      INGRID_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Premix NT constraints
      BRIAN_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Raw material analyses
      RAWMAT_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Final product added RM constraints
      STRAIGHTS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Specification tonnes
      TONS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Raw Material prices
      PRM_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Raw material penalty prices
      PRMS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Premix inclusion limits
      MIDELMT_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Global RM inclusion limits
      RMINC_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Product nutrient ratios
      NURATP_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)') STRING ! Premix nutrient ratios
      NURATM_FILE = trim(dirname)//trim(STRING)
*     
      READ(IUNIT,'(A)') STRING ! Iteration log file
      LOGOUTFILE = trim(dirname)//trim(STRING)
*     
      READ(IUNIT,'(A)') STRING ! Results file (machine readable)
      RESULTSFILE = trim(dirname)//trim(STRING)
*
      READ(IUNIT,'(A)') STRING  ! Skip the solution list file
      READ(IUNIT,'(A)') STRING  ! Skip the solution expost file
*
c     ... set default names: these lines are optional in sonoco.def
      RMLMTS_FILE   ='RMLMTS.WRK'
      SPECS_FILE    ='SPECS.WRK'
      SPECIAL_FILE  ='SPECIAL.WRK'
*
      READ(IUNIT,'(A)',END=100) STRING ! Overall RM limits file
      RMLMTS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)',END=100) STRING ! Formula nutrient specs (unused)
      SPECS_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)',END=100) STRING ! Summary file (in debug mode)
      SUMMARY_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)',END=100) STRING ! Relationships file
      SPECIAL_FILE = trim(dirname)//trim(STRING)
      READ(IUNIT,'(A)',END=100) STRING ! Nonlinear constraints def'tion
      NONLINEAR_FILE = trim(dirname)//trim(STRING)
*     
 100  CLOSE(IUNIT)

c     ... setup constant filenames
      SUMMARY_FILE  ='SUMMARY.WRK'
      SLPCT_FILE = 'slp_ct.dat'
      RNDSD_FILE = trim(dirname)//'/rnd_sd.dat'
      inquire(file=RNDSD_FILE,EXIST= file_exists)
      if (file_exists) then
         print*, ' read random seed from file ? '
         read *, ch
c         if (ch.ne.'y') RNDSD_FILE = ' '
         if (ch.ne.'y') 
     &        RNDSD_FILE = trim(dirname)//'/dummy'
      end if

c      open(31,file=dirname(1:trim77(dirname))//'/newdim.dat')
c      read(31,*) nu, rm, mi, de
c      close(31)

      end

c     -----------------------------------------------------------------
c     BLOCKDATA
c     -----------------------------------------------------------------
c     This blockdata block is used to provide default settings
c     for SLP control variables. 
c     The defaults can be changed by including the corresonding
c     COMMON from the calling program and changing the value


      blockdata setup
      
      include 'msg.inc'
      include 'SLPCOMM.INC'
      data msg_reg, msg_net, msg_rep_cs, msg_tght_cs,msg_warn_edge, 
     &     msg_err, msg_log, msg_solprt, msg_lamprt, msg_ranprt,
     &     msg_nonlin, msg_multstart, msg_wrt_sumdat
     &     /13 * .false./
      data NUMBER_STARTS, PROB_GLOBAL /50,95/

      end blockdata
