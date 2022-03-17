c     ==============================================================
c
c     subroutines for variable and constraint registration
c
c     ==============================================================
c
c     Every variable that is used by the model must be registered 
c     first. This can be done anytime before it is first referred
c     to. A variable is identified by three integer values:
c        1)   type:    MU, M, C, U, NU
c        2/3) ix1/ix2: index/indices, vars with only one index 
c                      have ix2=0
c     All the variable registering information is stored in the
c     array vr_reg(2+2*sn+4*n)  n = max number of vars
c                              sn = sqrt(n)
c     
c     Constraints can be registered by using the 3-int identification
c     of the variables as described above.
c     ixl, ixn1, ixn2 are (3,...) arrays that hold the 3-int information
c     for all the varioables in question.


c     structure of c/s registration entries
c     
c     in iuser at pi_cs_pi/pi_cs_pu are pointers to the start of 
c     integer/real information for each constraint. Each constraint has
c     the following information:
c     IUSER: 0              number of linear terms (n_lt)
c            1              number of non-linear (bilinear) terms (n_nt)
c            2              start of integer info for linear terms
c            3              start of integer info for bilinear terms
c            4              start of real info for linear terms
c            5              start of real info for bilinear terms
c            6              type of constraint
c            7              ix 1
c            8              ix 2
c            9              info on linear terms (1 entry per term)
c            9+n_lt         info on nonlinear terms (2 entries per term)
c            9+n_lt+2*n_nt  start of next constraint
c     USER:  0              start of info for linear terms (1 entry per term)
c            n_lt           start of info for bilin terms (1 entry per term)
c            n_lt + n_nt    start of info for next c/s

      subroutine reg_cs(iuser, user, sz_iuser, sz_user, type, ix1, ix2, 
     &     n, m, sn, vr_reg, bdscl, n_lt, n_nt, ixl, ixn1, ixn2, col, 
     &     con, lb, ub, scl, ifail)

c     RETURN VALUES:
c       ifail = 0: okay
c             = 1: out of memory
c
      implicit none

      integer sz_iuser, sz_user, ifail
      integer type, ix1, ix2, n_lt, n_nt, n, sn, m
      integer ixl(3,n_lt), ixn1(3,n_nt), ixn2(3,n_nt), vr_reg(*)
      double precision col(n_lt), con(n_nt), lb, ub, scl
      integer iuser(sz_iuser)
      double precision user(sz_user), bdscl(3*(n+m+1))
      integer i, i1, i2

      include 'pusr.inc'
      include 'msg.inc'
      include 'return.inc'
      include 'SLPCOMM.INC'


      ifail = 0

      if (REG_FST_CL.eq.0) then
         REG_FST_CL = 1
         pi_cs_pi = 0
         pi_cs_pu = m+1
         pi_cs_da = 2*(m+1)
         pu_cs_da = 0
         nx_nty = pi_cs_da+1
         nx_dnty = pu_cs+1
      end if

      REG_NX_CS = REG_NX_CS + 1

      if (msg_reg) then
         write(nout, '(A,I5,3I4,A,2I4)') ' Register c/s no ',
     &        REG_NX_CS, type, ix1, ix2,' : ',n_lt, n_nt
         write(nout, *) ' Starting at position ',nx_nty, nx_dnty
      end if


c     check for enough space to register constraints
      if ((REG_NX_CS.gt.m).and.(type.ne.CTYPE_OBJ)) then
         print '(A,I5,A,I5)', 
     &        'Attempt to register c/s ',REG_NX_CS,' of ',m
         stop
      end if


      if (nx_nty + 9+n_lt+2*n_nt.gt.sz_iuser) then
         if (msg_err) then
            write(nout, *) 'ERROR: out of integer memory in reg_cs.'
         end if
         global_err_msg = 'ERROR: out of integer memory in reg_cs.'
         global_ifail = 21
         ifail = 1
         goto 99
      end if
      if (nx_dnty + n_lt + n_nt .gt. sz_user) then
         if (msg_err) then
            write(nout,*) 'ERROR: out of real memory in reg_cs.'
         end if
         global_err_msg = 'ERROR: out of real memory in reg_cs.'
         global_ifail = 21
         ifail = 1
         goto 99
      end if

      iuser(pi_cs_pi+REG_NX_CS) = nx_nty
      iuser(pi_cs_pu+REG_NX_CS) = nx_dnty
      iuser(nx_nty+0) = n_lt
      iuser(nx_nty+1) = n_nt
c     ... set start of integer info for linear/bilinear terms
      iuser(nx_nty+2) = nx_nty+9
      iuser(nx_nty+3) = nx_nty+9+n_lt
c     ... set start of real info for linear/bilinear terms
      iuser(nx_nty+4) = nx_dnty
      iuser(nx_nty+5) = nx_dnty+n_lt

      iuser(nx_nty+6) = type
      iuser(nx_nty+7) = ix1
      iuser(nx_nty+8) = ix2
c     ... set bounds
      bdscl(n+REG_NX_CS)           = lb
      bdscl(n+m+1+n+REG_NX_CS)     = ub
      bdscl(2*(n+m+1)+n+REG_NX_CS) = scl
      
c     ... set linear terms
      nx_nty = nx_nty + 9
      do i=1,n_lt
         call fd_vr(n, sn, vr_reg, ixl(1,i), ixl(2,i), ixl(3,i), i1)
         if (i1.lt.0) then
            write(nout,*) 
     &           'Could not find var:',ixl(1,i), ixl(2,i), ixl(3,i)
            stop
         end if
         iuser(nx_nty) = i1
         user(nx_dnty) = col(i)
         nx_nty = nx_nty +1
         nx_dnty = nx_dnty +1
      end do
c     ... set bilinear terms
      do i=1,n_nt
         call fd_vr(n, sn, vr_reg, ixn1(1,i), ixn1(2,i), ixn1(3,i), i1)
         call fd_vr(n, sn, vr_reg, ixn2(1,i), ixn2(2,i), ixn2(3,i), i2)
         if (i1.lt.0) then
            write(nout, *)
     &           'Could not find var:',ixn1(1,i),ixn1(2,i),ixn1(3,i)
            stop
         end if
         if (i2.lt.0) then
            write(nout,*) 
     &           'Could not find var:',ixn2(1,i),ixn2(2,i),ixn2(3,i)
            stop
         end if
         iuser(nx_nty)   = i1
         iuser(nx_nty+1) = i2
         user(nx_dnty)   = con(i)
         nx_nty = nx_nty +2
         nx_dnty = nx_dnty +1
      end do

c      print *,'end of reg_cs'
 99   continue

      end


c     ----------------------------------------------------------------
c     reg_vr
c     ----------------------------------------------------------------

      subroutine reg_vr(n, m, sn, vr_reg, bdscl, type, ix1, ix2, lb, 
     &     ub, scl)

c     structure of vr_reg array:
c               1   n_vars
c               1   empty
c               sn  pointer to first element with hashcode
c               sn  pointer to most recent element with hashcode
c               n   type
c               n   ix1
c               n   ix2
c               n   hashcode

      implicit none

      integer n, m, sn, vr_reg(*), type, ix1, ix2
      integer hashcode, i
      double precision bdscl(3*(n+m+1)), lb, ub, scl

      include 'msg.inc'
      include 'pusr.inc'
      include 'SLPCOMM.INC'

      if (REG_NX_VR.eq.0) then
c     ... setup arrays
         do i=1,sn
            vr_reg(2+i) = 0
            vr_reg(2+sn+i) = 0
         end do
         do i=1,n
            vr_reg(2+2*sn+3*n+i) = 0
         end do
      end if

      REG_NX_VR = REG_NX_VR + 1

c      vr_reg(1) = vr_reg(1) + 1
      hashcode = mod(mod(type+ix1+ix2,sn)+sn,sn) + 1
      if (msg_reg) 
     &     write(nout, '(a,6I6,2g10.3)') 'reg var:',
     &     REG_NX_VR, n, type, ix1, ix2, hashcode, lb, ub

      if (REG_NX_VR.gt.n) then
         print '(A,I5,A,I5)', 
     &        'Attempt to register var ',REG_NX_VR,' of ',n
         stop
      end if

      vr_reg(2+2*sn+0*n+REG_NX_VR) = type 
      vr_reg(2+2*sn+1*n+REG_NX_VR) = ix1
      vr_reg(2+2*sn+2*n+REG_NX_VR) = ix2
      
      if (vr_reg(2+hashcode).eq.0) then
c     ... first entry with this hashcode
         vr_reg(2+hashcode) = REG_NX_VR
         vr_reg(2+sn+hashcode) = REG_NX_VR
      else
         vr_reg(2+2*sn+3*n+vr_reg(2+sn+hashcode)) = REG_NX_VR
         vr_reg(2+sn+hashcode) = REG_NX_VR
      end if
      
      bdscl(REG_NX_VR)           = lb 
      bdscl(n+m+1+REG_NX_VR)     = ub 
      bdscl(2*(n+m+1)+REG_NX_VR) = scl
c      call prt_vr(n, m, sn, vr_reg) 

      end

c     ----------------------------------------------------------------
c     fd_vr
c     ----------------------------------------------------------------

      subroutine fd_vr(n, sn, vr_reg, type, ix1, ix2, ix)
      
      integer n, sn, vr_reg(*), type, ix1, ix2, ix
      integer hashcode, p

      hashcode = mod(mod(type+ix1+ix2,sn)+sn,sn)+1
      p = vr_reg(2+hashcode)

 100  continue
      if (p.eq.0) then
         ix = -1
         return
      end if

      if (vr_reg(2+2*sn+p).eq.type.and.vr_reg(2+2*sn+n+p).eq.ix1.and.
     &     vr_reg(2+2*sn+2*n+p).eq.ix2) then
         ix = p
         return
      end if
      p = vr_reg(2+2*sn+3*n+p)
      goto 100

      end


c     ****************************************************************
c
c     subroutine tighten_bounds_for_constraints:
c        finds all constraints of class type (M, C, NU, U)
c        and uses the information to derive tighter bounds for the
c        variable defined by this constraint
c
c     ****************************************************************
c
c     this routine assumes that all weights (linear, nonlinear) have to sum 
c     up to 1. i.e. constraints is of the form
c      
c           -x + sum_i mu_i col_i + sum_j nu_j con_j y_j = 0
c       and                      sum_i mu_i + sum_j nu_j = 1
c
c     col_i, con_j    coefficients (constants)
c     mu_i, nu_j   linear, nonlinear weights, can have bounds outside [0,1]
c
c     Assumptions are justified for M, C constraints, but NOT U, NU-cs
c     
c     The routine scans for bounds on mu_i, nu_j that are outside [0,1] 
c     and sets   smneg = sum (-lb(mu_i,nu_j))+
c                smpos = sum (ub(mu_i, nu_j)-1)+
c     further scans for lv, uv (=lowest, highest) coefficient 
c     con*lb(y), con*ub(y) in case of nonlinear terms and sets new 
c     bounds
c
c        lb = lv - min(smneg,smpos)*(uv-lv)
c        ub = uv + min(smneg,smpos)*(uv-lv)

c     dChange
c     boundMax passed here from rd_prod_da.f
      subroutine tgt_bnd_cs(type, iuser, user, n, sn, vr_reg, bdscl, 
     &     ifail, boundMax)
c     end dChange

c     PARAMETERS:
c      boundMax(1:nu): maximal level of nutrient i in any raw material
c      boundMax(nu+1): sum of all product demands (total production)
c      boundMwx(nu+2): maximal cost for any raw material
c      bdscl:   array of lower/upper bounds and scale factors for all
c               variables and constraints:
c               1         - n+m    :  lb (first vars then cons)
c               n+m+1     - 2*(n+m):  lu (first vars then cons)
c               2*(n+m)+1 - 3*(n+m):  scale (first vars then cons)
c
c     RETURN VALUES
c       ifail = 0: okay
c             = 1: tightened bounds are now infeasible
c     
      use interface_mod

      implicit none

      include 'pusr.inc'
      include 'msg.inc'
      include 'SLPCOMM.INC'
      include 'return.inc'

      character*2 LFCR
      
      integer type, n, sn, ifail
      integer iuser(*), vr_reg(*)
      double precision user(*), bdscl(*)

      integer m, i, j, p0i, p0d, pi1, pi2, pd1, pd2, nl, nn, vr,
     &     ix, ix1, ix2
      double precision lb, ub, smpos, smneg, uv, lv, co, tol, slb

c     dChange
c     boundMax
c     boundMax is maximum amounts of nutrient,mass,cost
      double precision boundMax(0:nu+2) 
c     end dChange

chm      LFCR = char(13)//char(11)
      LFCR = char(13)//char(10)
      ifail = 0
      m = n_cs
      
c     ... search all constraints
      do i=1,m
         p0i = iuser(pi_cs_pi+i)
         p0d = iuser(pi_cs_pu+i)
         
         if (iuser(p0i+6).eq.type) then
            if (msg_tght_cs)
     &           write(nout, '(A,I4,3I3,A,2I3)') 
     &           'try to tighten using c/s:',i, iuser(p0i+6),
     &           iuser(p0i+7), iuser(p0i+8),' : ',
     &           iuser(p0i), iuser(p0i+1)
c        ... found constraint - identify variable to tighten
            pi1 = iuser(p0i+2)
            pd1 = iuser(p0i+4)
            pi2 = iuser(p0i+3)
            pd2 = iuser(p0i+5)
            nl  = iuser(p0i+0)
            nn  = iuser(p0i+1)
            vr = iuser(pi1)
            
            if (msg_tght_cs)
     &           write(nout,'(A,I4,A,3I3)')
     &           'try to tighten bounds on variable:',
     &           vr,':',vr_reg(2+2*sn+vr), vr_reg(2+2*sn+n+vr),
     &           vr_reg(2+2*sn+2*n+vr)
            
            if (abs(user(pd1)+1.d0).gt.1d-6) then
               write(nout,*) 'cannot tighten bound: '//
     &              'weight on first term \= -1:',user(pd1)
               stop
            end if

c     The defined variable is a (almost) convex combination of fixed
c     terms (linear terms) and variable terms (bilinear terms).
c
c     We need:   - sum of negative allowable weights             smneg
c                - sum of allowable weights above 1.0            smpos
c                - lowest value to contribute to weighted sum    lv  
c                - highst value to contribute to weighted sum    uv  
c
c     The bounds are then calculated as
c             tol = min(smneg,smpos)
c             lb = lv - tol*(uv-lv)
c             ub = uv + tol*(uv-lv)
c     also calculate slb (simple lb) which is the sum over lower bounds
c     of all terms involved. For some constraints this may give tighter
c     bounds. An analog sub does not make sense, since a value of 1
c     for *all* lambda will give a large overestimate
            
            smneg = 0.d0
            smpos = 0.d0
            lv = 1.d10
            uv = -1.d10
            slb = 0.d0

c           ... go through linear terms
            do j=2,nl
               ix = iuser(pi1+j-1)
               co = user(pd1+j-1)
               if (co.gt.uv) uv = co
               if (co.lt.lv) lv = co
               if (bdscl(ix).lt.0.d0) smneg = smneg - bdscl(ix)
               if (bdscl(n+m+1+ix).gt.1.d0) smpos = smpos 
     &              + bdscl(n+m+1+ix) - 1.d0
               if (co.gt.0.d0) then
                  slb = slb + co*bdscl(ix)
               else
                  slb = slb + co*bdscl(n+m+1+ix)
               end if
            end do
            
c           ... go through bilinear terms
            do j=1,nn
               ix1 = iuser(pi2+2*j-2)
               ix2 = iuser(pi2+2*j-1)
               co = user(pd2+j-1)
               if (co.lt.0.d0) then
                  write(nout, *) 
     &                 'negative bilinear weight: cannot get bound'
                  stop
               end if

               if (co*bdscl(ix2).lt.lv) lv = co*bdscl(ix2)
               if (co*bdscl(n+m+1+ix2).gt.uv) uv = co*bdscl(n+m+1+ix2)
               if (bdscl(ix1).lt.0.d0) smneg = smneg - bdscl(ix1)
               if (bdscl(n+m+1+ix1).gt.1.d0) smpos = smpos 
     &              + bdscl(n+m+1+ix1) - 1.d0
c     ... add slb term
               if (bdscl(ix1).ge.0.d0.and.bdscl(ix2).ge.0d0) then
                  slb = slb + bdscl(ix1)*bdscl(ix2)
               else
                  if (bdscl(ix1)*bdscl(n+m+1+ix2).lt.
     &                 bdscl(ix2)*bdscl(n+m+1+ix1)) then
                     slb = slb + bdscl(ix1)*bdscl(n+m+1+ix2)
                  else
                     slb = slb + bdscl(ix2)*bdscl(n+m+1+ix1)
                  end if
               end if
               
            end do
            
            if (msg_tght_cs) then
               write(nout, '(A,2F8.3)')
     &              'min/max allowable violation of convex bounds:',
     &              smneg, smpos
               write(nout, *) 'Possible cvx range of values:',lv, uv
            end if

            tol = min(smneg,smpos)
            lb = lv - tol*(uv-lv)
            ub = uv + tol*(uv-lv)
            if (msg_tght_cs) then
               write(nout, *) 'Possible range of values    :',lb, ub
               write(nout, *) 'Simple lower bound          :',slb
               write(nout, *) 'Old bounds                  :',
     &              bdscl(vr), bdscl(n+m+1+vr)  
               write(nout, *)
            end if
            
c     ... use slb if that is tighter
            if (slb>lb) lb = slb

c     ... check that bounds are feasible

            if (max(lb, bdscl(vr)) > min(ub, bdscl(n+m+1+vr))) then
c     .. if variable to be tightened is M-Type
               if (vr_reg(2+2*sn+vr).eq.1) then
chm                  write(global_err_msg, FMT=1001) LFCR, 
chm     &                 vr_reg(2+2*sn+2*n+vr),
chm     &                    vr_reg(2+2*sn+n+vr), LFCR, lb, ub, LFCR, 
chm     &                    bdscl(vr), bdscl(n+m+1+vr)
                  write(global_err_msg, '(A)')
     &   'ERROR: Infeasible problem specification.'//LFCR
     & //'Nutrient         : xxxxxxxxxx'//LFCR
     & //'Product          : xxxxxxxxxx'//LFCR
     & //'Original Min/Max : xxxxxxxxxxxxxxxx  xxxxxxxxxxxxxxxx'//LFCR
     & //'Possible Minimum : xxxxxxxxxxxxxxxx  xxxxxxxxxxxxxxxx'
                  global_int_arg(1)    = vr_reg(2+2*sn+2*n+vr)
                  global_int_arg(2)    = vr_reg(2+2*sn+n+vr)
                  global_double_arg(1) = lb
                  global_double_arg(2) = ub
                  global_double_arg(3) = bdscl(vr)
                  global_double_arg(4) = bdscl(n+m+1+vr)
               else
                  write(global_err_msg, '(A)')
     &                 'ERROR: infeasible problem specification.'
               end if
               if (msg_err) then
                  write(nout,'(A)') global_err_msg
               end if
c                  write(nout, *) 'ERROR: infeasible problem specification'
cc     .. if variable to be tightened is M-Type
c                  if (vr_reg(2+2*sn+vr).eq.1) then
c                     print '(A,I3,A,I3)','   Nutrient level NU-', 
c     &                    vr_reg(2+2*sn+2*n+vr),
c     &                    ', MI/DE-', vr_reg(2+2*sn+n+vr)
c                     print '(A,2G15.8)',
c     &                    '   Bounds from brian.dat are:         ',
c     &                    bdscl(vr), bdscl(n+m+1+vr)
c                     print '(A,2G15.8)',
c     &                    '   Derived bounds from inclusion are: ',
c     &                    lb, ub
c                     write(6,FMT=1001) vr_reg(2+2*sn+2*n+vr),
c     &                    vr_reg(2+2*sn+n+vr), lb, ub, bdscl(vr), 
c     &                    bdscl(n+m+1+vr)
c                  end if
c               end if
               ifail = 1
chm               global_ifail = 1
               global_ifail = 5
               return
            end if

c     ... tighted bound

            if (lb.gt.bdscl(vr))       bdscl(vr)       = lb
            if (ub.lt.bdscl(n+m+1+vr)) bdscl(n+m+1+vr) = ub


c     dChange
c     perform extra bound tightening usnig boundMax
c     For the case of CTYPE_M/CTYPE_C (nutrient level and cost of mixes) 
c     constraints we can derive tighter bounds using values in the boundMax
c     array:
c     - Level of nutrient i in any mix has to be in [-maxL, 2*maxL] where
c       maxL is max level of the nutrient in any raw-material
c     - Cost of mix j has to be in [-maxC, 2*maxC] where
c       maxC is the max cost of any raw-material
c     AGR: are these ever active. Seems that these bounds should be picked up
c          by analysing the constraints as above.
c          Also would be good to have better bounds on possible weights than
c          [-1, 2].
            if (type==CTYPE_M) then
               if (1==1) then
c                 tighten upper bound on nutrient level in mix:
c                 has to be <= 2*(max level of nu in raw_materials) 
c                 ... bdscl(n+m+1+vr) is ub(vr)
                  if (bdscl(n+m+1+vr).gt.
     &                 2*boundMax(vr_reg(2+2*sn+2*n+vr))) then
                     if (msg_tght_cs) then
                        write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                 'type vr', type, vr,' old ub', bdscl(n+m+1+vr), 
     &                 ' new ub', 2*boundMax(vr_reg(2+2*sn+2*n+vr))
                     end if
                     bdscl(n+m+1+vr) = 2*boundMax(vr_reg(2+2*sn+2*n+vr))
                  end if
               end if
               if (1==1) then
c                 tighten lower bound on nutrient level in mix:
c                 has to be => -(max level of nu in raw_materials) 
c                 ... bdscl(vr) is lb(vr)
                 if (bdscl(vr).lt.-boundMax(vr_reg(2+2*sn+2*n+vr))) then
                     if (msg_tght_cs) then
                        write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                       'type vr', type, vr,
     &                       ' old lb', bdscl(vr),' new lb',
     &                       -boundMax(vr_reg(2+2*sn+2*n+vr))
                     end if
                     bdscl(vr) =  -boundMax(vr_reg(2+2*sn+2*n+vr))
                  end if
               end if
            end if
c     get costs between -max and 2*max
            if (type==CTYPE_C) then
               if (1==1) then
c                 tighten upper bound on cost of mix:
c                 has to be <= 2*(maxCost of raw_materials) 
c                 ... bdscl(n+m+1+vr) is ub(vr)
                  if (bdscl(n+m+1+vr).gt.2*boundMax(nu+2)) then 
                     if (msg_tght_cs) then    
                        write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                 'type vr', type, vr,' old ub', bdscl(n+m+1+vr), 
     &                 ' new ub', 2*boundMax(nu+2)
                     end if
                     bdscl(n+m+1+vr) = 2*boundMax(nu+2)
                  end if
               end if
               if (1==1) then
c                 tighten lower bound on cost of mix:
c                 has to be >= -(maxCost of raw_materials) 
c                 ... bdscl(vr) is lb(vr)
                  if (bdscl(vr).lt.-boundMax(nu+2) ) then
                     if (msg_tght_cs) then
                        write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                   'type vr', type, vr,
     &                   ' old lb', bdscl(vr),' new lb', -boundMax(nu+2)
                     end if
                     bdscl(vr) =  -boundMax(nu+2)
                  end if
               end if
            end if
c     end dChange



         end if
      end do

 1001 FORMAT("ERROR: infeasible problem specification",A2,
     &     "Nutrient ",I3, 
     &     "  in MI/DE ",I3,A2,
     &     "Derived bounds: Min ",G15.8," Max ",G15.8,A2,
     &     "Set bounds:     Min ",G15.8," Max ",G15.8)

      end


c     ****************************************************************
c
c     subroutine tighten_bounds_for_constraints:
c        finds all constraints of class type (M, C, NU, U)
c        and uses the information to derive tighter bounds for the
c        variable defined by this constraint
c
c     ****************************************************************
c
c     this routine just assumes that the constraint is of the form
c
c       -x + sum_i mu_i col_i + sum_j nu_j con_j y_j = 0
c       (weights don't have to sum to 1)
c
c     and is applicable to NU and U constraints. But gives weaker bounds
c     than the alternative tgt_bnd_cs
c
c     simply sets   lb = sum_i lb(mu_i) col_i + sum_j lb(nu_j*y_j) con_j
c                   ub = sum_i ub(mu_i) col_i + sum_j ub(nu_j*y_j) con_j
c

c     dChange
c     boundMax passed here from rd_prob_da.f
      subroutine tgt_bnd_cs_smp(type, iuser, user, n, sn, vr_reg, bdscl,
     &     boundMax)
c     end dChange

c     PARAMETERS:
c      boundMax(nu+1): sum of all product demands (total production)

      use interface_mod

      implicit none

      include 'pusr.inc'
      include 'SLPCOMM.INC'
      include 'msg.inc'

      integer type, n, sn
      integer iuser(*), vr_reg(*)
      double precision user(*), bdscl(*)

      integer m, i, j, p0i, p0d, pi1, pi2, pd1, pd2, nl, nn, vr,
     &     ix, ix1, ix2
      double precision lb, ub, smpos, smneg, uv, lv, co, tol, 
     &     l1, l2, u1, u2
      
c     dChange
c     declare boundMax
c     boundMax is maximum amounts of nutrient,mass,cost
      double precision boundMax(0:nu+2) 
c     end dChange

c     ... search all constraints
      m = n_cs
      do i=1,m
         p0i = iuser(pi_cs_pi+i)
         p0d = iuser(pi_cs_pu+i)
         
         if (iuser(p0i+6).eq.type) then
            if (msg_tght_cs)
     &           write(nout, '(A,I5,3I4,A,2I4)') 
     &           'try to tighten using c/s:',i, iuser(p0i+6),
     &           iuser(p0i+7), iuser(p0i+8),' : ',
     &           iuser(p0i), iuser(p0i+1)
c        ... found constraint - identify variable to tighten
            pi1 = iuser(p0i+2)
            pd1 = iuser(p0i+4)
            pi2 = iuser(p0i+3)
            pd2 = iuser(p0i+5)
            nl  = iuser(p0i+0)
            nn  = iuser(p0i+1)
            vr = iuser(pi1)
            
            if (msg_tght_cs)
     &       write(nout, '(A,I5,A,3I4)') 
     &           'try to tighten bounds on variable:',
     &           vr,':',vr_reg(2+2*sn+vr), vr_reg(2+2*sn+n+vr),
     &           vr_reg(2+2*sn+2*n+vr)
            
            if (abs(user(pd1)+1.d0).gt.1d-6) then
               write(nout, *) 'cannot tighten bound: '//
     &              'weight on first term \= -1:',user(pd1)
               stop
            end if

c     do simple bound finding: for each term calculate lowest and highest
c     possible value and sum them up

            lb = 0.d0
            ub = 0.d0

c           ... go through linear terms
c               co*mu: co*lb(mu), cu*ub(mu) are bounds
            do j=2,nl
               ix = iuser(pi1+j-1)
               co = user(pd1+j-1)
               lb = lb + co*bdscl(ix)
               ub = ub + co*bdscl(n+m+1+ix)
            end do

c           ... go through bilinear terms: co*mu*y: 
c                co*min{lb(mu)lb(y),ub(mu)lb(y),lb(mu)ub(y),ub(mu)ub(y)}
            do j=1,nn
               ix1 = iuser(pi2+2*j-2)
               ix2 = iuser(pi2+2*j-1)
               co = user(pd2+j-1)
               l1 = bdscl(ix1)
               l2 = bdscl(ix2)
               u1 = bdscl(n+m+1+ix1)
               u2 = bdscl(n+m+1+ix2)

c               print *, j, co, l1, l2, u1, u2
               lv = min(min(l1*l2, l1*u2), min(u1*l2, u1*u2))
               uv = max(max(l1*l2, l1*u2), max(u1*l2, u1*u2))

               lb = lb + lv
               ub = ub + uv
            end do

            if (msg_tght_cs) then
               write(nout, *) 'Possible range of values    :',lb, ub
               write(nout, *) 'Old bounds                  :',
     &              bdscl(vr), bdscl(n+m+1+vr)  
               write(nout, *)
            end if

c     ... tighted bound

            if (lb.gt.bdscl(vr))       bdscl(vr)       = lb
            if (ub.lt.bdscl(n+m+1+vr)) bdscl(n+m+1+vr) = ub


c     dChange
c     perform extra tightening with boundMax for type four cons
c     put masses between -max and 2*max
c     For constraints of type CTYPE_U, tighten bounds on U variable
c     (total flow through rm/bin): has to be within [-max, 2*max],
c     where max is the total production (sum of all demands)
c     For constraints of type CTYPE_NU, tighten bounds on NU variable
c     (total mass of raw material i in bin/product) to be within [-max, 2*max],
c     where max is as above
c     AGR: This seems true, but unlikely ever to be activated. Lower bound
c     on all these variables should be 0 anyhow. For upper bound it would
c     be good to get a better factor than 2*.
            if (1==1) then
               if (bdscl(n+m+1+vr).gt.2*boundMax(nu+1))then 
                  if (msg_tght_cs) then
                     write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                 'type vr', type, vr,' old ub', bdscl(n+m+1+vr), 
     &                    ' new ub', 2*boundMax(nu+1) 
                  end if
                  bdscl(n+m+1+vr) = 2*boundMax(nu+1)   
               end if
            end if
            if (1==1) then
               if (bdscl(vr).lt.-boundMax(nu+1) ) then
                  if (msg_tght_cs) then
                     write(nout, '(a,2i6,a,g12.4,a,g12.4)') 
     &                    'type vr', type, vr,' old lb', bdscl(vr), 
     &                    ' new lb', -boundMax(nu+1) 
                  end if
                  bdscl(vr) = -boundMax(nu+1)
               end if
            end if
c     end dChange


         end if
      end do
      

      end



      subroutine prt_cs(n, sn, iuser, user, vr_reg)

      implicit none

      include 'pusr.inc'
      include 'SLPCOMM.INC'

      integer n, sn
      integer iuser(*), vr_reg(*)
      double precision user(*)

      integer m, i, j, ix1, ix2, p0i, p0d, pi1, pi2, pd1, pd2, nl, nn

      write(nout, *)
     &     '------------------------------------------------------'
      write(nout, *) 'Total no of c/s:', n_cs
      m = n_cs
      do i=1,m
         p0i = iuser(pi_cs_pi+i)
         p0d = iuser(pi_cs_pu+i)
         pi1 = iuser(p0i+2)
         pd1 = iuser(p0i+4)
         pi2 = iuser(p0i+3)
         pd2 = iuser(p0i+5)
         nl  = iuser(p0i+0)
         nn  = iuser(p0i+1)

c         write(nout, *)  i, p0i, p0d, pi1, pi2, pd1, pd2
  
         write(nout, *) '    #  TP   I  J: LIN NL-terms' 
         write(nout, '(A,I4,A,3I3,A,2I3,A,2G15.8)') 'cs',i, ' -',
     &        iuser(p0i+6), iuser(p0i+7), iuser(p0i+8),':',
     &        nl, nn
         do j=1,nl
            ix1 = iuser(pi1+j-1)
            write(nout, '(F8.3,I4,A,3I4)') user(pd1+j-1), ix1,':', 
     &           vr_reg(2+2*sn+ix1), 
     &           vr_reg(2+2*sn+n+ix1), vr_reg(2+2*sn+2*n+ix1)
         end do
         do j=1,nn
            ix1 = iuser(pi2+2*j-2)
            ix2 = iuser(pi2+2*j-1)
            write(nout, '(F8.3,2I4,A,3I4,A,3I4)')
     &           user(pd2+j-1),ix1,ix2,':', 
     &   vr_reg(2+2*sn+ix1),vr_reg(2+2*sn+n+ix1),vr_reg(2+2*sn+2*n+ix1),
     &           '|',
     &   vr_reg(2+2*sn+ix2),vr_reg(2+2*sn+n+ix2),vr_reg(2+2*sn+2*n+ix2)
         end do

         write(nout, *)
      end do

      end
      
c     ----------------------------------------------------------------
c     prt_vr
c     ----------------------------------------------------------------
      
      subroutine prt_vr(n, m, sn, vr_reg) 
      
      implicit none

      integer n, m, sn, vr_reg(*)
      
      integer nvr, i, type, ix1, ix2, p, hash

      include 'pusr.inc'
      include 'SLPCOMM.INC'
      nvr = REG_NX_VR
      
      write(nout, *)  'HASHCODEs:  HASH   first   last'

      do i=1,sn
         write(nout, '(3I5)') i, vr_reg(2+i), vr_reg(2+sn+i)
      end do

      write(nout, *) 'VAR: # TYPE  IX1  IX2  hash next-same-hash'
      do i=1,nvr
         type = vr_reg(2+2*sn+i)
         ix1 = vr_reg(2+2*sn+n+i)
         ix2 = vr_reg(2+2*sn+2*n+i)
         p   = vr_reg(2+2*sn+3*n+i)
         hash = mod(mod(type+ix1+ix2,sn)+sn,sn) + 1
         write(nout, '(6I5)') 
     &        i, type, ix1, ix2, hash, p
      end do

      end
