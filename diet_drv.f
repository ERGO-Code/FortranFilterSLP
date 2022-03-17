c     This is the driver routine. Everything done here is 
c     the same that Henrik would do in his GUI
c     - call to get_dimensions
c     - allocate arrays
c     - call slp_driver
c

      program diet_drv

      use interface_mod

      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'msg.inc'
      include 'return.inc'

c     ... size of EMSOL workarrays: need to be passed to slpdriver
      integer sz_ws, sz_lws
      parameter (sz_ws = 15000000, sz_lws = 5100000)
c      parameter (sz_ws = 15000000, sz_lws = 5000000)

      integer nresults
      parameter (nresults = 8)

      integer lws(sz_lws)
      double precision ws(sz_ws)
      double precision, allocatable ::  ntgresult(:)

      integer mx_big_it

      integer err, ifail

      msg_reg = .false.        !log from c/s, var registering
      msg_net = .false.        !print bin network
      msg_rep_cs = .false.     !report all constraints
      msg_tght_cs = .false.    !report on constraint tightening
      msg_warn_edge = .false. !warn for missing edge (and other warnings)
      msg_err       = .true.   !all error messages 
      msg_log       = .false.  !log from rd_prob_da: reading file
c      msg_solprt  = .false.  !print solution
      msg_solprt  = .true.  !print solution
      msg_lamprt  = .false.  !print sensitivity
c      msg_lamprt  = .true.  !print sensitivity
      msg_ranprt  = .false.  !print ranging
      msg_nonlin  = .false.  !debug from nonlinear reading
      msg_multstart = .false. !log from multistart logic
      msg_wrt_sumdat = .true. !write iteration log in sumXXX.dat files

      print *,' Number of random starts : '
      read *, mx_big_it
      NUMBER_STARTS = mx_big_it

cf90      write(6,'(A)',advance='no') ' Name of problem directory : '
      print*, ' Name of problem directory : '
      read '(a20)', dirname
c      dirname='jan99'


c     ... call to get_dimensions. This reads sonoco.def to get the
c     filenames of the problem definition files and calls rdprobdi
c     to setup problem dimensions. 
c     All the data read in is stored in interface_mod

      call get_dimensions(trim(dirname)//'/sonoco.def', ifail)


      allocate(ntgresult(nresults*ntgdim), stat=err)
      if (err.ne.0) stop

     
      call slp_driver(ws, sz_ws, lws, sz_lws, ntgresult, nresults,
     &     ifail)

      
      end


      subroutine get_emserror(name, rtcod, errmessage)
      implicit none

      character*8 name
      character*80 errmessage
      integer rtcod

      end

      subroutine logout(two, dspace, bigiter, iter, normd, rho, fnew,
     .     hcnew, nr_lp, st)
      implicit none
      integer two, bigiter, iter, nr_lp, st
      double precision dspace(*), normd, rho, fnew, hcnew

      end

      subroutine logset(m, n, d)
      implicit none
      integer m, n, d
      end


      subroutine test_bit(no, bit, res)
      implicit none

      integer no, bit, res

      res = 0
      if (btest(no,32-bit)) res = 1

      end


