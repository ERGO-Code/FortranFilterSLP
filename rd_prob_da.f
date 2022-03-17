      subroutine rd_prob_da(n, m, sn, user, iuser, sz_user, sz_iuser, 
     &     vr_reg, bdscl, oldbnd, ws, nws, lws, nlws, bns, n_bns, 
     &     max_cs_terms, SPD, DAT, ifail, NTGRESULT,NRESULTS)

c
c     rd_prob_da needs  ws[ 2+3*nu + max{2*nu, 2*(rm+mi), 2*nrst} ]
c                      lws[ mi+max(de, rm)+ 
c                          max{2*nrst+rm, nurat_p, nurat_m} ]
c
c     PARAMETERS:
c      bdscl:   array of lower/upper bounds and scale factors for all
c               variables and constraints:
c               1         - n+m    :  lb (first vars then cons)
c               n+m+1     - 2*(n+m):  lu (first vars then cons)
c               2*(n+m)+1 - 3*(n+m):  scale (first vars then cons)
c
c     RETURN VALUES:
c      ifail = 0: read ok 
c            = 1: out of memory in reg_cs
c            = 2: infeasible in presolve (tgt_bnd_cs)
c            = 3: exceeded hardwired limit on number of c/s terms
c            = 99: fatal: bail out with error message and
c                  global_ifail set here
c
c
c
c
c     ------------------------------------------------------------------
c     NEGATIVE WEIGHTS:
c
c     Raw materials into bins and raw mats into products (straights)
c     are allowed to have negative weights
c
c     All negative weights below -limnwt (-2.d0) are rejected
cHM?  All negative weights below -limnwt (-9.d0) are rejected
c     Bins cannot have negative costs
c
c     ------------------------------------------------------------------
c
c     The model uses the following variables:
c     TYPE_M (nu,mi+de)     nutrient content of bins/demands    
c     TYPE_MU(rm+mi, mi+de) weight of rm/bin in bin/demand
c     TYPE_U (rm+mi)        total flow through rm/bin 
c                            (sparse, only when needed) 
c     TYPE_NU(rm, mi+de)    total weight of rm in mi/de
c                             (sparse: only when needed, i.e mi/de  
c                              is in source path of de and mi is not
c                              first order bin (no bin sources))
c     TYPE_C(mi+de)         price of bin/demand
c
c     and the following constraints:
c     CTYPE_M(nu,mi+de):
c              m_ij = sum_{k \in sources} mu_jk*m_ik
c     CTYPE_C(mi+de):
c              c_j  = sum_{k \in sources} mu_jk*c_k
c     CTYPE_SUM(mi+de):
c              sum_{k \in sources} mu_jk = 1  \forall j 
c     CTYPE_NURAT(l, k):
c              nutrient ratio constraint #l on mi/de k
c              l is couted separately for NURATM/NURATP
c     CTYPE_U(rm+mi): 
c              define U: sparse only generated when needed
c     CTYPE_NU
c              define NU: sparse only generated when needed
c
c     objective is (CTYPE_OBJ):
c              sum_{k \in demands} tons_k*C_k

cFIXME: check sizes of user, iuser array
      
      use interface_mod
      use types_module
      use nonlin_types
      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'pusr.inc'
      include 'msg.inc'
chm
      include 'return.inc'
      include 'ntgres.inc'

      type(special_da):: SPD
      type(integra_da):: DAT
c      type(nl_info):: NLD
      integer n, m, sn, sz_user, sz_iuser, n_bns, nws, nlws, ifail,
     &     max_cs_terms, cntlin, cntnl
      double precision scl
      integer iuser(sz_iuser), lws(nlws), vr_reg(4*n+2*sn+2), bns(n_bns)
      double precision bdscl(3*n+3*m+3), user(sz_user), ws(nws)
      double precision oldbnd(2*n+2*m)
      integer NRESULTS
      double precision NTGRESULT(*)
      integer trim77
c      
chc   New parameters
c
      integer NDX,RESULTINDEX
      integer RATIOCOUNT

      integer i, j, k, l, nznu, i0, i1, i2, i3, nrst, pcsi, ix, os
      integer srb, srbb, srbr, p1, fifo1, fifo2, tgb, tgbb, tgbd

      integer p_scl, p_maxnu, p_da, p_minnu, p_rmi_lb, p_rmi_ub,
     &     p_rmi_fg, p_lda
      double precision tol, com_tol, totde, avde
      double precision dummy, maxv, minv, limnwt, lb, ub
      double precision dp1, dp2, dp1l, dp1u, dp2l, dp2u

      logical file_exists
      logical edge_exists

      logical fi_xst, adv_it, lp_sol, found
      character ch
      character(180) line
      integer err

      double precision     infty, eps
      common /NLP_eps_inf/ infty, eps

c     to build sparse constraint structure array
      integer lt_ix1(3,max_cs_terms), nt_ix1(3,max_cs_terms), 
     &     nt_ix2(3,max_cs_terms)
      double precision lt_co(max_cs_terms), nt_co(max_cs_terms)
      integer fifo(2*MAX_RMS)

c     automatic arrays to carry temporary information
c     ... bounds on rm availability (set by rmlmts, used later)
      double precision rmava_lb(rm+mi), rmava_ub(rm+mi)

c     dChange 
c     boundMax is maximum amounts of nutrient,mass,cost
      double precision boundMax(0:nu+2) 
c     end dChange

      ifail = 0

c     ... zero static variables for reg_cs/reg_vr
      REG_FST_CL = 0
      REG_NX_CS = 0
      REG_NX_VR = 0

c     ...pointer into user
      pu_r      = rm
      pu_d0     = pu_r  + nu*rm
      pu_cs     = pu_d0 + de

c     ... pointer into iuser
      pi_cs     = 0
      pi_cs_pi  = pi_cs    + 0
      pi_cs_pu  = pi_cs_pi + m + 1
      pi_cs_da  = pi_cs_pu + m + 1

cc     ... map of sources for each bin/demand
c      bn_n_srb   = 0
c      bn_n_srbb  = bn_n_srb   + mi+de
c      bn_p_srb   = bn_n_srbb  + mi+de
c      bn_l_srb   = bn_p_srb   + mi+de
cc     ... map of targets for each rm/bin
c      bn_n_tgb   = bn_l_srb   + n_mu
c      bn_n_tgbb  = bn_n_tgb   + rm+mi
c      bn_p_tgb   = bn_n_tgbb  + rm+mi
c      bn_l_tgb   = bn_p_tgb   + rm+mi
cc     ... list of first order bins
c      bn_fob     = bn_l_tgb   + n_mu
c                              +mi+de

c     ... pointers into workspace ws
c        ...avarage value of nutrient & average price
      p_scl       = 0
c        ... maximal value of nutrient
      p_maxnu     = p_scl    + nu+2
c        ... min of values of nutrients
      p_minnu     = p_maxnu  + nu
c        ... lb/ub for rmava/rminc constraints
      p_rmi_lb    = p_minnu  + nu
      p_rmi_ub    = p_rmi_lb  + mi + max(de,rm)
c        ... data vector to read in values from file
      p_da        = p_rmi_ub + mi + max(de,rm) 
c                         + max(2*nu, rm, 2*mi)
c     pointers into workspace lws
      p_rmi_fg    = 0

      p_lda       = p_rmi_fg + mi + max(de,rm)


c      n_vr = 0
c      n_cs = 0
      n_maxa = 0


      limnwt = 2.d0             ! Stick to this 22/7-05
chm   limnwt = 9.d0      ! Causes high instability and float point error
chm   limnwt = 4.d0      ! Does work on the SDF 3603 but not preferred.


c     limnwt sets the maximum acceptable negative weight
c     used to calculate bounds on lam, mu, m, etc

      lp_sol = .false.

c     ================ subroutine body =====================
c     --------------- read raw material prices: PRM.DAT ---------

c     user[1:rm] holds raw material prices.
c     ws[p_scl+nu+1] holds average price.

      ws(p_scl+nu+1) = 0.
      if (msg_log) then
         WRITE(nout,*)' reading raw material prices: ',
     &        PRM_FILE(1:trim77(PRM_FILE))
      end if
      open(10,file=PRM_FILE)
      do i=1,rm
c         DAT%is_piecewise_price(i) = .false.
         read(10,'(a)') line
c         print *, i, line
         read(line,*, iostat=err) DAT%price1(i), DAT%price2(i)
c         print *, i, err, DAT%price1(i), DAT%price2(i)
c        .. this here works: if there is only one record: err = -1
c           if there are two records then err = 0
         if (err.eq.0.and.DAT%price2(i).lt.500000) then
c         if (err.eq.0) then
            if (DAT%is_piecewise_price(i) .neqv. .true.) then
               print *,'previously identified as needing PWL'
               stop
            end if
         else
            if (DAT%is_piecewise_price(i) .neqv. .false.) then
               print *,'previously identified as not needing PWL'
               stop
            end if

         end if

         user(i) = DAT%price1(i)
         ws(p_scl+nu+1) = ws(p_scl+nu+1) + user(i)
         ws(p_scl+nu+1) = ws(p_scl+nu+1) + user(i)
      end do
      close(10)
      ws(p_scl+nu+1) = ws(p_scl+nu+1)/rm
      if (abs(ws(p_scl+nu+1)).le.1d-8) ws(p_scl+nu+1) = 1.d0


c     - - - - - register cost variables for PWL raw mat - - - - - - 
      do i=1,rm
         if (DAT%is_piecewise_price(i)) then
            lb = 0.d0
            ub = 1.d10
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_C, -i, 0, lb, ub,
     &           ws(p_scl+nu+1))
         end if
      end do

c     ------------- read raw material specs: RAW_MAT.DAT ----------

c     user[pu_r+1:pu_r+rm*nu] holds raw material definition
c     ws[p_da+1:p_da+rm] used as temporary array 

      if (msg_log) then
         WRITE(nout,*) 
     &        ' reading nutrient composition of raw materials: ',
     &        RAWMAT_FILE(1:trim77(RAWMAT_FILE))
      end if
      open(10,file=RAWMAT_FILE)
      do l=1,nu
         read(10,*) (ws(p_da+i), i=1,rm)
c      WRITE(nout,*) (ws(p_da+i), i=1,rm)

         do i=1,rm
            user(pu_r+(i-1)*nu+l) = ws(p_da+i)
         end do
      end do

      
      close(10)

c     ... calculate maxnu & minnu

      do l=1,nu
         ws(p_maxnu+l) = 0.
         ws(p_minnu+l) = 1.d10
         do i=1,rm
            ws(p_maxnu+l) = max(ws(p_maxnu+l),user(pu_r+(i-1)*nu+l))
            ws(p_minnu+l) = min(ws(p_minnu+l),user(pu_r+(i-1)*nu+l))
         end do
      end do

c     -------------- read tonnages of demand: TONS.DAT ---------------

      open(10,file=TONS_FILE)
      if (msg_log) then
         WRITE(nout,*) ' reading tonnages for demand: ',
     &        TONS_FILE(1:trim77(TONS_FILE))
      end if
      totde = 0.d0
      do j=1,de
         read(10,*) user(pu_d0+j)
         totde = totde + user(pu_d0+j)
c         WRITE(nout,*) user(pu_d0+j)
      end do
      close(10)
      avde = totde/de

c     ----------------------- calculate scales --------------
c
c     scl(i), 1<=i<=nu is average level of nu i in rm's
c     scl(nu+1)        is average price for rm's
c     scl(nu+2)        is average demand for prod's

      ws(p_scl+nu+2) = avde
      do l=1,nu
         nznu = 0
         ws(p_scl+l) = 0.
         do i=1,rm
            if (user(pu_r+(i-1)*nu+l).gt.1.e-10) then
               nznu = nznu + 1
               ws(p_scl+l) = ws(p_scl+l) + user(pu_r+(i-1)*nu+l)
            end if
         end do
         if (nznu.eq.0.or.(abs(ws(p_scl+l)).le.1d-8)) then
            ws(p_scl+l) = 1.d0
         else
            ws(p_scl+l) = ws(p_scl+l)/nznu
         end if
      end do
c      WRITE(nout,*) ' Print scales:'
c      do i=1,nu+1
c         WRITE(nout,'(I3,F15.5)') i,ws(p_scl+i)
c      end do

c     -------- read bounds on nutrients in bins: BRIAN.DAT ---------
c
c     Every premix is described by two lines of BRIAN.DAT:
c       first gives lower bounds for all nutrients
c       second gives upper bounds for all nutrients

      inquire(file=BRIAN_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*)
     &           ' reading bounds on nut. content of premixes: ',
     &           BRIAN_FILE(1:trim77(BRIAN_FILE))
         end if
         open(10,file=BRIAN_FILE)
         do j=1,mi
            read(10,*) (ws(p_da+i), i=1,nu)

chc         remember the original bounds to return to GUI
            do i=1,nu
              NDX=RESULTINDEX(j,i,2,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+i)
            end do
c
            read(10,*) (ws(p_da+nu+i), i=1,nu)
chc
            do i=1,nu
              NDX=RESULTINDEX(j,i,3,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+nu+i)
            end do
c
            do i=1,nu
               lb = ws(p_da+i)
               ub = ws(p_da+nu+i)
               
               call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_M, j, i, 
     &              lb, ub, ws(p_scl+i))
            end do
            lb = -1.d10
            ub = 1.d10
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_C, j, 0, lb, ub,
     &             ws(p_scl+nu+1))
         end do
         close(10)
      else
c     ... register M, C variables
         do j=1,mi
            do i=1,nu
               lb = -1.d10
               ub = 1.d10
               call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_M, j, i, 
     &              lb, ub, ws(p_scl+i))
            end do
            lb = -1.d10
            ub = 1.d10
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_C, j, 0, lb, ub,
     &             ws(p_scl+nu+1))
         end do
      end if



c     ----- read nutrient specs of demands: SPECS.DAT/LMTS.DAT -------

      inquire(file=LMTS_FILE,exist=file_exists)
      if (file_exists) then
         com_tol = 0.0
         if (msg_log) then
            WRITE(nout,*)' reading bounds for demands: ',
     &           LMTS_FILE(1:trim77(LMTS_FILE))
            WRITE(nout,*)'  tolerance on bounds is:',com_tol
         end if
         open(10,file=LMTS_FILE)

         do j=1,de
            read(10,*) (ws(p_da+i), i=1,2*nu)

chc         remember the original bounds to return to GUI
            k=0
            do i=1,nu
              k=k+1
c     Max constraint               
              NDX=RESULTINDEX(j+mi,i,3,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+k)
              k=k+1
c  Min constraint               
              NDX=RESULTINDEX(j+mi,i,2,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+k)
            end do

            do l=1,nu
               tol = max(abs(ws(p_da+2*l)),abs(ws(p_da+2*l-1)))*com_tol
               call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_M, mi+j, l,
     &              ws(p_da+2*l)-tol,  ws(p_da+2*l-1)+tol,
     &              ws(p_scl+l))
c               bup(pv_dk+(j-1)*nu+l) = ws(p_da+2*l-1)*(1.+com_tol)
c               blo(pv_dk+(j-1)*nu+l) = ws(p_da+2*l)*(1.-com_tol)
            end do
            lb = -1.d10
            ub = 1.d10
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_C, mi+j, 0, 
     &           lb, ub, ws(p_scl+nu+1))

         end do
         close(10)
      else
         if (msg_log) then
            WRITE(nout,*)' reading specifications for demands: ',
     &           SPECS_FILE(1:trim77(SPECS_FILE))
            WRITE(nout,*) '  tolerance for specs: '
         end if
         open(10,file=SPECS_FILE)
         read(10,*) com_tol
         do j=1,de
            read(10,*) (ws(p_da+i), i=1,nu)
            do i=1,nu
               tol = abs(ws(p_da+i))*com_tol
               call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_M, mi+j, l,
     &             ws(p_da+i)-tol,ws(p_da+i)+tol, ws(p_scl+i))
c               blo(pv_dk+(j-1)*nu+i) = ws(p_da+i)*(1.-com_tol)
c               bup(pv_dk+(j-1)*nu+i) = ws(p_da+i)*(1.+com_tol)
            end do
            lb = -1.d10
            ub = 1.d10
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_C, mi+j, 0, 
     &           lb, ub, ws(p_scl+nu+1))
         end do
         close(10)
      end if


c     ------------ read bounds on bin-inclusions: INGRID.DAT ---------

c     INGRID.DAT now contains bounds on all weights:
c                bin/raw_mat maekup of all bins/products

c                The old INGRID.DAT/STRAIGHTS.DAT/MIDELMT.DAT
c                are now redundant.
c
c     INGRID.DAT now also contains information on the sparsity of the 
c                inclusion graph: Edges with LB=UB=0 simply do not
c                exsist. 
c
c     A variable is registered (with appropriate bounds) for every edge
c     that does exist.
c     
c     A description of the network has already been generated in
c     rd_prob_dim: it set up the following arrays:
c     (bn_???? are pointers to the first elements to the array in bns)
c
c     bn_n_srb   [mi+de]       #sources for bin/demand
c     bn_n_srbb  [mi+de]          of which other bins
c     bn_p_srb   [mi+de]       pointer into list of src for each bn/de
c     bn_l_srb   [n_mu]        list of sources
c     bn_n_tgb   [rm+mi]       #targets for rm/bin
c     bn_n_tgbb  [rm+mi]          of which other bins
c     bn_p_tgb   [rm+mi]       pointer into list of trg for each rm/bn
c     bn_l_tgb   [n_mu]        list of targets
c     bn_fob     [mi+de]       first order bins (no bin sources)

      inquire(file=INGRID_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*) ' reading bounds on bin inclusions: ',
     &           INGRID_FILE(1:trim77(INGRID_FILE))
         end if
         open(10,file=INGRID_FILE)
         do j=1,mi+de
c           ... read lower bounds
            read(10,*) (ws(p_da+i), i=1,rm+mi)
chc
            do i=1,rm+mi
              NDX=RESULTINDEX(j,i+nu,2,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+i)
            end do

c           ... read upper bounds
            read(10,*) (ws(p_da+rm+mi+i), i=1,rm+mi)
chc
            do i=1,rm+mi
              NDX=RESULTINDEX(j,i+nu,3,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=ws(p_da+rm+mi+i)
            end do

            do i=1,rm+mi
               if (max(abs(ws(p_da+rm+mi+i)),abs(ws(p_da+i))).gt.1e-8) 
     &              then
                  lb = min(max(ws(p_da+i),-limnwt),1.d0+limnwt)
                  ub = min(max(ws(p_da+rm+mi+i),-limnwt),1.d0+limnwt)

                  call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_MU, i, j, 
     &                 lb, ub, 1.d0) 
               end if
            end do
         end do
         close(10)
      else
c     ... no INGRID.DAT exists: take default
         print *, INGRID_FILE
         print *,INGRID_FILE(1:trim77(INGRID_FILE))//
     &        ' does not exists: No predefined action'
         stop
      end if

c     ---------------- generate M, C, SUM constraints ------------------
c
c     Now that all weight variables and the structure of the network
c     is in place, we can generate the M, C, SUM constraints for 
c     every BIN/PRODUCT
c

c     ... print network if wanted
      if (msg_net) call prt_bn(bns, n_bns)


      do i=1,mi+de
c      ... for all bins + products
         srb  = bns(bn_n_srb+i)
         srbb = bns(bn_n_srbb+i)
         srbr = srb-srbb
         if (srb+1.gt.max_cs_terms) then
            ifail = 3
            return
         end if
            
         p1 = bns(bn_p_srb+i)

c        ... set nutrient content constraint for each nutrient
         do j=1,nu
c           ... set linear terms
            lt_ix1(1,1) = TYPE_M
            lt_ix1(2,1) = i
            lt_ix1(3,1) = j
            lt_co(1)  = -1.d0
            do k=1,srbr
               i1 = bns(p1+k)
               lt_ix1(1,1+k) = TYPE_MU
               lt_ix1(2,1+k) = i1
               lt_ix1(3,1+k) = i
               lt_co(1+k) = user(pu_r+(i1-1)*nu+j)
            end do
c           ... set nonlinear terms
            do k=1,srbb
               i1 = bns(p1+srbr+k)-rm
               nt_ix1(1,k) = TYPE_MU
               nt_ix1(2,k) = i1+rm 
               nt_ix1(3,k) = i
               nt_ix2(1,k) = TYPE_M
               nt_ix2(2,k) = i1 
               nt_ix2(3,k) = j
               nt_co(k)  = 1.d0
            end do
            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_M, i, j, 
     &           n, m, sn, vr_reg, 
     &           bdscl, srbr+1, srbb, lt_ix1, nt_ix1, nt_ix2, lt_co,
     $           nt_co,0.d0, 0.d0, ws(p_scl+j), ifail)
            if (ifail.ne.0) return
         end do
c        >>price constraint
c        ... set linear terms
         lt_ix1(1,1) = TYPE_C
         lt_ix1(2,1) = i
         lt_ix1(3,1) = 0
         lt_co(1)  = -1.d0
         cntlin = 1
         cntnl = 0
c     ... loop over all raw-material sources of this bin
         do k=1,srbr
            i1 = bns(p1+k)
c           ... check if raw-mat has a piecewise linear price
            if (DAT%is_piecewise_price(i1)) then
c              ... set of a nonlinear term
               cntnl = cntnl + 1
               nt_ix1(1,cntnl) = TYPE_MU
               nt_ix1(2,cntnl) = i1
               nt_ix1(3,cntnl) = i
               nt_ix2(1,cntnl) = TYPE_C
               nt_ix2(2,cntnl) = -i1
               nt_ix2(3,cntnl) = 0
               nt_co(cntnl)  = 1.d0
            else
               cntlin = cntlin + 1
               lt_ix1(1,cntlin) = TYPE_MU
               lt_ix1(2,cntlin) = i1
               lt_ix1(3,cntlin) = i
               lt_co(cntlin) = user(i1)
            end if
         end do
c        ... set nonlinear terms
         do k=1,srbb
            i1 = bns(p1+srbr+k)-rm
            cntnl = cntnl + 1
            nt_ix1(1,cntnl) = TYPE_MU
            nt_ix1(2,cntnl) = rm+i1
            nt_ix1(3,cntnl) = i
            nt_ix2(1,cntnl) = TYPE_C
            nt_ix2(2,cntnl) = i1
            nt_ix2(3,cntnl) = 0
            nt_co(cntnl)  = 1.d0
         end do
         call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_C, i, 0, n, 
     &        m, sn, vr_reg, 
     &        bdscl, cntlin, cntnl, lt_ix1, nt_ix1, nt_ix2, lt_co, 
     &        nt_co ,0.d0, 0.d0, ws(p_scl+nu+1), ifail)
         if (ifail.ne.0) return
c        >>sum = 1 constraint     
         do k=1,srb
            i1 = bns(p1+k)
            lt_ix1(1,k) = TYPE_MU
            lt_ix1(2,k) = i1
            lt_ix1(3,k) = i
            lt_co(k)  = 1.d0
         end do
         call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_SUM, i, 0, n,
     &        m, sn, vr_reg, 
     &        bdscl, srb, 0, lt_ix1, nt_ix1, nt_ix2, lt_co, nt_co, 
     &        1.d0, 1.d0, 1.d0, ifail)
         if (ifail.ne.0) return
      end do
      
c     ==============================================================
c     
c     This has been the basic model. The rest involves variables 
c     that are setup in a sparse manner:
c      
c      - nu : total inclusion of rm in a certain bin/prod
c      - u  : total flow through a certain rm/bin
c
c      - NU can be specified for [rm][bin/prod]
c           if any NU is specified (due to a rminc constrint = bound, 
c           or a special constraint) then all NU for this [rm] that
c           might lead to the current node in the tree must also
c           be specified
c      - u  can be specified for [rm/bin]
c           if any U is specified (due to a rmlmts constraint = bound,
c           or a special constraint) then all U above the current node
c           in the tree must also be specified.
c

c     ------------------- read rmlmts.dat -------------------------

c     FORMAT:
c     [no of raw mats/bins with availability constraints]
c     [ix of rm/bin] [lb] [ub]     (first rm are raw-mat then bins)
c     [ix of rm/bin] [lb] [ub]
c     ...
      
c      vars: rm_mi_lb(rm+mi), rm_mi_ub(rm+mi), rm_mi_fg(rm+mi)
c            fifo(rm+mi), fifo1, fifo2
c
c     set up lists: rm_mi_fg(rm+mi):   1 if rmava c/s needed
c                   rm_mi_lb(rm+mi):   lb if needed
c                   rm_mi_ub(rm+mi):   ub if needed


      do i=1,rm+mi
         lws(p_rmi_fg+i) = 0
         rmava_lb(i) = -0.05*totde
         rmava_ub(i) = 1.01d0*totde
      end do
      fifo1 = 1
      fifo2 = 1

      inquire(file=RMLMTS_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*) 
     &           ' reading bounds on raw material availability: ',
     &           RMLMTS_FILE(1:trim77(RMLMTS_FILE))
         end if
         open(10,file=RMLMTS_FILE)
         read (10,*) nrst

         do i=1,nrst
            read (10, '(a)') line
            read (line, *, iostat=err) i1, dp1l, dp1u, dp2l, dp2u
c            read (10,*) i1, dp1, dp2
            if (i1.le.rm+mi) then
               rmava_lb(i1) = dp1l
               rmava_ub(i1) = dp1u

c               print *, i1, dp1l, dp1u
               DAT%rmava1lb(i1) = dp1l
               DAT%rmava1ub(i1) = dp1u
               if (err.eq.0) then
                  if (i1.le.rm) then
                     if (.not.DAT%is_piecewise_price(i1)) then
                        print '(A,I4,A)', 'Raw material ',i1,
     &                       ': two sets of bounds but not two prices' 
c                        stop
                     end if
                     DAT%rmava2lb(i1) = dp2l
                     DAT%rmava2ub(i1) = dp2u
c                    print *, i1, dp2l, dp2u
                     rmava_lb(i1) = dp1l+dp2l
                     rmava_ub(i1) = dp1u+dp2u
                  else
                     if (msg_warn_edge) then
                        print '(A,I4,A)', 
     &               'WARNING: two set of lb/ub given for rm/mi ',i1, 
     &               ' which is a bin'
                     end if
                     
                  end if
               end if   
               lws(p_rmi_fg+i1) = 1
               fifo(fifo2) = i1
               fifo2 = fifo2 + 1
chc         ... remember original bounds to return to GUI
               NDX=RESULTINDEX(de+mi+1,i,2,nu+rm+mi,NRESULTS)
               NTGRESULT(NDX)=dp1
               NDX=RESULTINDEX(de+mi+1,i,3,nu+rm+mi,NRESULTS)
               NTGRESULT(NDX)=dp2
            else
               print *, 'No raw material ',i1
               stop
            end if
c
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c     ... scan special file for any U variables mentioned
         
         do k=1,SPD%n
            do i=SPD%p_cs(k),SPD%p_cs(k)+SPD%n_terms(k)-1
               if (SPD%vtype(i).eq.TYPE_U) then
                  i1 = SPD%spec1(i)
                  if (lws(i1).eq.0) then
                     lws(i1) = 1
                     fifo(fifo2) = i1
                     fifo2 = fifo2 + 1
                  end if
                  if (SPD%n_terms(k).eq.1) then
c                     print *,'BND>', SPD%lb(k), SPD%coeff(i), SPD%ub(k)
                     rmava_lb(i1) = max(rmava_lb(i1),
     &                    SPD%lb(k)/SPD%coeff(i))
                     rmava_ub(i1) = min(rmava_ub(i1),
     &                    SPD%ub(k)/SPD%coeff(i))
                  end if
               end if
            end do
         end do

c     ...   flag all bins/rm that need implied rmava constraints 
 10      continue
         i1 = fifo(fifo1)
         i2 = bns(bn_p_tgb+i1)
         fifo1 = fifo1 + 1
         do i=1,bns(bn_n_tgb+i1)
            i3 = bns(i2+i)+rm
            if (i3.le.rm+mi) then
               if (lws(p_rmi_fg+i3).eq.0) then
                  lws(p_rmi_fg+i3) = 1
                  fifo(fifo2) = i3
                  fifo2 = fifo2 + 1
               end if
            end if
         end do
         if (fifo1.lt.fifo2) goto 10
c     .... register all vars
         
         do i=1,rm+mi
            if (lws(p_rmi_fg+i).eq.1) then
               call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_U, i, 0, 
     &              rmava_lb(i), rmava_ub(i), ws(p_scl+nu+2))
            end if
         end do
         
c     ... register constraints
         do i=1,rm+mi
            if (lws(p_rmi_fg+i).eq.1) then
               
               tgb = bns(bn_n_tgb+i)
               tgbb = bns(bn_n_tgbb+i)
               tgbd = tgb-tgbb
               if (tgb+1.gt.max_cs_terms) then
                  ifail = 3
                  return
               end if
               
               p1 = bns(bn_p_tgb+i)
c               ...set bilinear terms
               do k=1,tgbb
                  i1 = bns(p1+k)
                  nt_ix1(1,k) = TYPE_MU
                  nt_ix1(2,k) = i 
                  nt_ix1(3,k) = i1
                  nt_ix2(1,k) = TYPE_U
                  nt_ix2(2,k) = i1+rm 
                  nt_ix2(3,k) = 0
                  nt_co(k)  = 1.d0
               end do
c              ... set linear terms
               lt_ix1(1,1) = TYPE_U
               lt_ix1(2,1) = i
               lt_ix1(3,1) = 0
               lt_co(1)  = -1.d0
               do k=1,tgbd
                  i1 = bns(p1+tgbb+k)
                  lt_ix1(1,1+k) = TYPE_MU
                  lt_ix1(2,1+k) = i
                  lt_ix1(3,1+k) = i1
c                  lt_co(1+k) = ws(i1-mi)
                  lt_co(1+k) = user(pu_d0+i1-mi)
               end do

               call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_U, i, 
     &              0, n, m, sn, vr_reg, bdscl,tgbd+1, tgbb, lt_ix1, 
     &              nt_ix1, nt_ix2, lt_co, nt_co,0.d0, 0.d0, 
     &              ws(p_scl+nu+2), ifail)
               if (ifail.ne.0) return
            end if
         end do

      end if

c     ------------------- read rminc.dat -------------------------

c     FORMAT:
c
c     [no of c/s of type RMINC]
c     [ix of bin/product] [ix of raw_mat] [lower bound] [upper bound]
c     [ix of bin/product] [ix of raw_mat] [lower bound] [upper bound]
c     [ix of bin/product] [ix of raw_mat] [lower bound] [upper bound]
c     [ix of bin/product] [ix of raw_mat] [lower bound] [upper bound]
c     ...

cFIXME: saved local info in ws(p_da), lws(p_lda)     

c
c     work out which bin/prod-rm combinations need a NU (total ingred)
c     variable
c

      nrst = 0
      do i=1,rm
         lws(p_lda+i) = 0
      end do
      inquire(file=RMINC_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*)
     &           ' reading bounds on raw_mat inclusion in products: ',
     &           RMINC_FILE(1:trim77(RMINC_FILE))
         end if
c        ... read file into temporary arrays:
c            lws(p_lda+       i)  bin/de affected by cs-i
c            lws(p_lda+  nrst+i)  rm affected by cs-i
c            lws(p_lda+2*nrst+i)  no of c/s affecting rm-i
c             ws(p_da+       i)  lower bound of cs-i
c             ws(p_da+  nrst+i)  upper bnd of cs-i

         open(10,file=RMINC_FILE)
         read(10,*) nrst

chc      ... remember to pass to getsensitivity via common
         NUM_GIG_RESULTS=nrst
c         
         do i=1,rm
            lws(p_lda+2*nrst+i) = 0
         end do
         do i=1,nrst
            read(10,*) lws(p_lda+i), lws(p_lda+nrst+i), 
     &           ws(p_da+i), ws(p_da+nrst+i)
c        ... count no of constraints for given raw-mat
            lws(p_lda+2*nrst+lws(p_lda+nrst+i)) = 
     &           lws(p_lda+2*nrst+lws(p_lda+nrst+i))+1
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c        ... for all raw_materials that have constraints
         do i=1,rm
c          ... get list of all bin/de affected by particular rm
            do j=1,mi+de
               lws(p_rmi_fg+j) = 0        !flag if NU variable needed
c     FIXME: these bounds should be tightened by working through
c            from the rm level how much of the rm could end up in
c            the bin in question. I'm not doing this now, there should
c            be no problem from fairly lax bounds             (26/10/06)
               ws(p_rmi_lb+j) = -1.d0          ! and bounds for it
               ws(p_rmi_ub+j) = 2.d0
            end do
c            for each rm make list of bins/de that have constraints
c            lws(p_rmi_fg+j): =1 if bin/de needs NU var for this rm
c             ws(p_rmi_lb+j): bounds on NU variable
c             ws(p_rmi_ub+j): 
            fifo1 = 1
            fifo2 = 1
c           fifo(fifo1:fifo2) is FIFO stack of bin/de with rminc constraints
            do j=1,nrst
               if (lws(p_lda+nrst+j).eq.i) then
c              ... mark all bin/de affected directly by constraints
                  k = lws(p_lda+j)
                  lws(p_rmi_fg+k) = 1
                  ws(p_rmi_lb+k) = ws(p_da+j)
                  ws(p_rmi_ub+k) = ws(p_da+nrst+j)
c                 ... and put the on the stack
                  fifo(fifo2) = k
                  fifo2 = fifo2 + 1
               end if
            end do
c     ... scan special file for any NU variables for this rm
            do k=1,SPD%n
               do j=SPD%p_cs(k),SPD%p_cs(k)+SPD%n_terms(k)-1
                  if ((SPD%vtype(j).eq.TYPE_NU).and.
     &                 SPD%spec1(j).eq.i) then
                     i1 = SPD%spec2(j)
                     if (lws(p_rmi_fg+i1).eq.0) then
                        lws(p_rmi_fg+i1) = 1
                        fifo(fifo2) = i1
                        fifo2 = fifo2 + 1
                     end if
                     if (SPD%n_terms(k).eq.1) then
c                     print *,'BND>', SPD%lb(k), SPD%coeff(j), SPD%ub(k)
                        ws(p_rmi_lb+i1) = max(ws(p_rmi_lb+i1),
     &                       SPD%lb(k)/SPD%coeff(j))
                        ws(p_rmi_ub+i1) = min(ws(p_rmi_ub+i1),
     &                       SPD%ub(k)/SPD%coeff(j))
                     end if
                  end if
               end do
            end do
c           ... and also include all bin-sources
            do
               if (fifo1.ge.fifo2) exit
c                get bin/product of the stack and find all sources of it
               i1 = fifo(fifo1)
c                i1 is a bin/product taken off the stack
               i2 = bns(bn_p_srb+i1)
c                i2 is a pointer (into bns) to a list of sources of i1  
               srb = bns(bn_n_srb+i1)
               srbb =  bns(bn_n_srbb+i1)
               srbr = srb-srbb
c                srb: total #sources, srbb: of which other bins
               fifo1 = fifo1 + 1
c                loop through all bin sources of i1
               do j=srbr+1,srb
                  i3 = bns(i2+j)-rm
                  if (lws(p_rmi_fg+i3).eq.0) then
                     lws(p_rmi_fg+i3) = 1
                     fifo(fifo2) = i3
                     fifo2 = fifo2 + 1
                  end if
               end do
            end do

c     at this point all bin/products with lws(p_rmi_fg+j), j=1,..,mi+de
c     need a NU variable (total raw mat) for ram material i

c           ... register all vars
            do j=1,mi+de
               if (lws(p_rmi_fg+j).eq.1) then
c                 ... if bin-j has bin sources, need NU variable
c                     (otherwise just use MU)
                  if (bns(bn_n_srbb+j).ne.0) then
                     call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_NU, 
     &                    i, j, ws(p_rmi_lb+j), 
     &                    ws(p_rmi_ub+j), 1.d0)
                  else
c        if there actually is a bound on NU in a first order bin, 
c        then this bound could be set throgh ingrid.dat
c        Nevertheless, a corresponding MU exists, which is tightened
c        if necessary
                     if (ws(p_rmi_lb+j).gt.0.d0.or.
     &                    ws(p_rmi_ub+j).lt.1.d0) then
c                  ... if first-order-bin then there is a corresponding
c                    MU variable => see if bounds need to be tightened
                        call fd_vr(n, sn, vr_reg, TYPE_MU, i, j, ix)
chm
                        if (ix.ge.0) then
c                         .. mu variable exists -> this should have been
c                            used
                           if (msg_err) then
                            write(nout,*)'WARNING: trying to set a'
     &   //' bound through rminc that should be set in ingrid'
                            write(nout,*) i, j, ws(p_rmi_lb+j), 
     &                                     ws(p_rmi_ub+j)
                          end if
                          write(global_err_msg, '(A)')
     &     ' WARNING: A constrained GIG is not allowed at this level. '
     &   //' Product: xxxxxxxxxx, Component: xxxxxxxxxx'
                          global_int_arg(1)    = i
                          global_int_arg(2)    = j
                          global_double_arg(1) = ws(p_rmi_lb+j)
                          global_double_arg(2) = ws(p_rmi_ub+j)
                          global_ifail = 50
                          ifail = 99
                          return
                        end if
                        if (ix.lt.0) then
c                         ... mu variable does not exits: this is
c                           an error in the logic that should be flagged
                          if (msg_err) then
                            write(nout,*)'Variable MU ',i,j,
     &    ' that should exists because of'                            
     &  //' NU involving first order bin, doesn''t exists'
                          end if
                          write(global_err_msg, '(A)')
     &    'Variable MU xxxxxxxxxx  xxxxxxxxxx that should exists'
     &  //' because of NU involving first order bin, doesn''t exists.'
                          global_int_arg(1) = i
                          global_int_arg(2) = j
                          global_ifail = 99
                          ifail = 99
                          return
                        else
                           if (ws(p_rmi_lb+j).gt.bdscl(ix))
     &                          bdscl(ix) = ws(p_rmi_lb+j)
                           if (ws(p_rmi_ub+j).lt.bdscl(n+m+1+ix))
     &                          bdscl(n+m+1+ix) = ws(p_rmi_ub+j)
                        end if
                     end if
                  end if
               end if
            end do
c              ... register constraints

            do j=1,mi+de
               if (lws(p_rmi_fg+j).eq.1
     &              .and.bns(bn_n_srbb+j).ne.0) then
                  srb = bns(bn_n_srb+j)
                  srbb =  bns(bn_n_srbb+j)
                  srbr = srb-srbb
c                    p1 = pointer to start of source for current mi/de
                  p1 = bns(bn_p_srb+j)
c                 ... set linear/bilinear terms
                     
                  lt_ix1(1,1) = TYPE_NU
                  lt_ix1(2,1) = i
                  lt_ix1(3,1) = j
                  lt_co(1) = -1.d0
                  i1 = 1
                  i2 = 0
                  do k=1,srb
c                 i3 counts through all sources
                     i3 = bns(p1+k)
c                    i3 is a source, it can be 
c                      - == i (the current raw_material)
c                      - > rm (a bin source)
                     if (i3.eq.i) then
                        i1 = i1 + 1
                        if (i1.gt.max_cs_terms) then
                           ifail = 3
                           return
                        end if
                        
                        lt_ix1(1,i1) = TYPE_MU
                        lt_ix1(2,i1) = i
                        lt_ix1(3,i1) = j
                        lt_co(i1) = 1.d0
                     else if (i3.gt.rm) then
                        if (bns(bn_n_srbb+i3-rm).eq.0) then
                           if (edge_exists(i, i3-rm, bns)) then
                              i2 = i2 + 1
                              if (i2.gt.max_cs_terms) then
                                 ifail = 3
                                 return
                              end if
                              nt_ix1(1,i2) = TYPE_MU
                              nt_ix1(2,i2) = i3
                              nt_ix1(3,i2) = j
                              nt_ix2(1,i2) = TYPE_MU
                              nt_ix2(2,i2) = i
                              nt_ix2(3,i2) = i3-rm
                              nt_co(i2) = 1.d0
                           else
                              if (msg_warn_edge) then
c     FIXME: what exactly is happening here and why can warning be
c     ignored?
                                 print *,'WARNING: MU ',i,i3-rm
                                 print *,'Left out in NU-c/s because'
     &                                //' no such edge exists'
                              end if
                           end if
                        else 
                           i2 = i2 + 1
                           if (i2.gt.max_cs_terms) then
                              ifail = 3
                              return
                           end if
                           nt_ix1(1,i2) = TYPE_MU
                           nt_ix1(2,i2) = i3
                           nt_ix1(3,i2) = j
                           nt_ix2(1,i2) = TYPE_NU
                           nt_ix2(2,i2) = i
                           nt_ix2(3,i2) = i3-rm
                           nt_co(i2) = 1.d0
                        end if
                     end if
                  end do
                  
                  call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_NU, 
     &                 i, j, n, m, sn,
     $                 vr_reg, bdscl, i1, i2, lt_ix1, nt_ix1, nt_ix2,
     $                 lt_co,nt_co,0.d0, 0.d0, 1.d0, ifail)
                  if (ifail.ne.0) return
               end if
            end do
            
         end do
      end if

c     --------------- register special constraints -------------------
c     see common_types.f for an explanation of the format

      do i=1,SPD%n
         if (SPD%n_terms(i).gt.max_cs_terms) then
            ifail = 3
            return
         end if
         if (SPD%n_terms(i).gt.1) then
            os = SPD%p_cs(i)-1
            do j=1,SPD%n_terms(i)
               lt_ix1(1,j) = SPD%vtype(os+j)
               lt_ix1(2,j) = SPD%spec1(os+j)
               lt_ix1(3,j) = SPD%spec2(os+j)
               lt_co(j) = SPD%coeff(os+j)
            end do
            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_SPECIAL, 
     &           i, 0, n, m, sn, vr_reg, bdscl, SPD%n_terms(i), 0, 
     &           lt_ix1, nt_ix1, nt_ix2, lt_co, nt_co, SPD%lb(i), 
     &           SPD%ub(i), 1.d0, ifail)
            if (ifail.ne.0) return
         end if
      end do
         

c     ------------------- read nurat.dat -------------------------

c     FORMAT:
c     [no of c/s of type NURAT]
c     [ix of prod] [no of coeff] [lb] [ub] [ix of nut 1] [coeff 1] [xi of n2] ...
c     [ix of prod] [no of coeff] [lb] [ub] [ix of nut 1] [coeff 1] [xi of n2] ...
c     [ix of prod] [no of coeff] [lb] [ub] [ix of nut 1] [coeff 1] [xi of n2] ...
c     ...

chm
      NUM_RATIO_RESULTS=0
      inquire(file=NURATP_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*) 
     &           ' reading bounds on nutrient ratios in products: ',
     &           NURATP_FILE(1:trim77(NURATP_FILE))
         end if
         open(10,file=NURATP_FILE)
chm Added RATIOCOUNT
c        ... NUM_RATIO_RESULTS counts how many NURAT type entries should
c        be reported back in NTGRESULT. The problem is that a constraint
c        on nutrient ratios may be represented by 1 or 2 linear 
c        constraints on nutrients of type NURAT. 
c     FIXME: this is pretty inelegant.
c     
         read(10,*) n_nuratp, RATIOCOUNT
         NUM_RATIO_RESULTS=NUM_RATIO_RESULTS+RATIOCOUNT
c
         do l=1,n_nuratp
            read(10,*) k, lws(p_lda+l)
         end do
         close(10)
         open(10,file=NURATP_FILE)
         read(10,*) nrst
         do l=1,n_nuratp
            if (lws(p_lda+l).gt.max_cs_terms) then
               ifail = 3
               return
            end if
            
            read(10,*) k, nrst, ws(p_da+1), ws(p_da+2),
     &           (lt_ix1(3,i),lt_co(i), i=1,lws(p_lda+l))

            scl = 0.d0
c     ... register constraint
            do i=1,lws(p_lda+l)
               lt_ix1(1,i) = TYPE_M 
               lt_ix1(2,i) = k+mi
c               lt_ix1(3,i) = 
c               lt_co(i)    = 
               scl = scl + abs(lt_co(i))*ws(p_scl+lt_ix1(3,i))
            end do
            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_NURAT, l, 
     &           k+mi, n, m, sn, 
     &           vr_reg, bdscl, lws(p_lda+l), 0, lt_ix1, nt_ix1, nt_ix2,
     $           lt_co, nt_co,ws(p_da+1), ws(p_da+2), scl, ifail)
            if (ifail.ne.0) return

         end do
         close(10)
      end if

      inquire(file=NURATM_FILE,exist=file_exists)
      if (file_exists) then
         if (msg_log) then
            WRITE(nout,*) 
     &           ' reading bounds on nutrient ratios in premixes: ',
     &           NURATM_FILE(1:trim77(NURATM_FILE))
         end if
         open(10,file=NURATM_FILE)
chm Added RATIOCOUNT
c        ... see above (NURATP) for RATIOCOUNT
         read(10,*) n_nuratm,RATIOCOUNT
         NUM_RATIO_RESULTS=NUM_RATIO_RESULTS+RATIOCOUNT
c
         do l=1,n_nuratm
            read(10,*) k, lws(p_lda+l)
         end do
         close(10)
         open(10,file=NURATM_FILE)
         read(10,*) nrst
         do l=1,n_nuratm
            if (lws(p_lda+l).gt.max_cs_terms) then
               ifail = 3
               return
            end if
            read(10,*) k, nrst, ws(p_da+1), ws(p_da+2),
     &           (lt_ix1(3,i),lt_co(i), i=1,lws(p_lda+l))

            scl = 0.d0
c     ... register constraint
            do i=1,lws(p_lda+l)
               lt_ix1(1,i) = TYPE_M 
               lt_ix1(2,i) = k
c              lt_ix1(3,i) = 
c              lt_co(i)    = 
               scl = scl + abs(lt_co(i))*ws(p_scl+lt_ix1(3,i))
            end do
            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_NURAT, l, 
     &           k, n, m, sn, 
     &           vr_reg, bdscl, lws(p_lda+l), 0, lt_ix1, nt_ix1, nt_ix2,
     $           lt_co, nt_co,ws(p_da+1), ws(p_da+2), scl, ifail)
            if (ifail.ne.0) return
         end do
         close(10)
      end if

c     ----------------- scan nonlinear cs file ------------------------

      call read_nonlinear(NLD, ifail, n, sn, vr_reg)
c     make sure that it is read in correctly
      if (ifail.ne.0) then 
c     ... bail out, global_ifail set
         ifail = 99 
         return
      end if

      if (msg_nonlin) call print_nl_data(n, sn, vr_reg)
c      print *, 'Read ',NLD%n_cs,' nonlinear constraints:'
c      do i=1,NLD%n_cs
c         print *, print_nl_cs(NLD%root(i)%p, ch, n, sn, vr_reg)
c      end do
c      stop
 
c     register dummy constraintsL with 0 linear and nonlinear terms
c     but so that the objective function is registered in the correct
c     location
      do i=1,NLD%n_cs
         call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_NONLINEAR, 
     &        i, 0, n, m, sn, vr_reg, bdscl, 0, 0, 
     &        lt_ix1, nt_ix1, nt_ix2, lt_co, nt_co, 
     &        NLD%lb(i), NLD%ub(i), 1.d0, ifail)
      end do

c     ------ set up the remainder of the piecewise linear costs ----
c     
c     Modelling of piecewise linear costs:
c     - raw material costs are given by two components, each with their
c       own availability limits. 
c     - This is modelled by introducing a TYPE_U variable that keeps track of
c       how much of each raw material is used (such a variable is generated
c       for every rawmat with rmlmts constraints anyhow)
c     - Introduce 5 additional variables:
c        + lam1, lam2: (TYPE_PWL, i, 1)/(TYPE_PWL, i, 2):  weights
c        + u_i1, u_i2: (TYPE_U, i, 1)/(TYPE_U, i, 2):      flows
c        + p_i (TYPE_C, -i):                               price
c     - and introcuce 4 additional constraints
c        + lam1*p1 + lam2*p2 = p1
c        + lam1 + lam2 = 1
c        + u_i1 = lam1*u_i
c        + u_i2 = lam2*u_i
c     - Overall problem dimensions (rm, mi) stay the same as before!

c
      do i=1,rm
         if (DAT%is_piecewise_price(i)) then
c     ... need constraints of the form
c            p_i = lam1*p1 + lam2*p2
c            lam1 + lam2 = 1
c         and
c            f1 = lam1*f_i, f2 = lam2*f_i
c     ...  f_i is the total flow though raw_mat i; TYPE_U: i, 0

            lb = 0.d0
            ub = 1.d0
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_PWL, i, 1, 
     &              lb, ub, 1.d0)
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_PWL, i, 2, 
     &              lb, ub, 1.d0)

            lb = DAT%rmava1lb(i)
            ub = DAT%rmava1ub(i)
c           ... p_scl+nu+2 is avde (average demand)
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_U, i, 1, 
     &              lb, ub, ws(p_scl+nu+2))

            lb = DAT%rmava2lb(i)
            ub = DAT%rmava2ub(i)
            call reg_vr(n, m, sn, vr_reg, bdscl, TYPE_U, i, 2, 
     &              lb, ub, ws(p_scl+nu+2))
            
c     ... and now register the 4  constraints
c     ...  - p_i + lam1*p1 + lam2*p2 = 0
            lt_ix1(1, 1) = TYPE_C
            lt_ix1(2, 1) = -i
            lt_ix1(3, 1) = 0
            lt_co(1) = -1.d0
            lt_ix1(1, 2) = TYPE_PWL
            lt_ix1(2, 2) = i
            lt_ix1(3, 2) = 1
            lt_co(2) = DAT%price1(i)
            lt_ix1(1, 3) = TYPE_PWL
            lt_ix1(2, 3) = i
            lt_ix1(3, 3) = 2
            lt_co(3) = DAT%price2(i)

            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_PWL, i, 0,
     &           n, m, sn, vr_reg, bdscl, 3, 0, lt_ix1, nt_ix1, nt_ix2, 
     &           lt_co, nt_co, 0.d0, 0.d0, ws(p_scl+nu+1), ifail)


c     ...  - lam1 + lam2 = 1.0
            lt_ix1(1, 1) = TYPE_PWL
            lt_ix1(2, 1) = i
            lt_ix1(3, 1) = 1
            lt_co(1) = 1.d0
            lt_ix1(1, 2) = TYPE_PWL
            lt_ix1(2, 2) = i
            lt_ix1(3, 2) = 2
            lt_co(2) = 1.d0

            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_SUM,-i, 0,
     &           n, m, sn, vr_reg, bdscl, 2, 0, lt_ix1, nt_ix1, nt_ix2, 
     &           lt_co, nt_co, 1.d0, 1.d0, ws(p_scl+nu+1), ifail)

c     ...  - -f1 + lam1*fl = 0
            lt_ix1(1, 1) = TYPE_U
            lt_ix1(2, 1) = i
            lt_ix1(3, 1) = 1
            lt_co(1) = -1.d0
            nt_ix1(1, 1) = TYPE_PWL
            nt_ix1(2, 1) = i
            nt_ix1(3, 1) = 1
            nt_ix2(1, 1) = TYPE_U
            nt_ix2(2, 1) = i
            nt_ix2(3, 1) = 0
            nt_co(1) = 1.d0

            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_U, i, 1,
     &           n, m, sn, vr_reg, bdscl, 1, 1, lt_ix1, nt_ix1, nt_ix2, 
     &           lt_co, nt_co, 0.d0, 0.d0, ws(p_scl+nu+2), ifail)

c     ...  - -f2 + lam2*f2 = 0
            lt_ix1(1, 1) = TYPE_U
            lt_ix1(2, 1) = i
            lt_ix1(3, 1) = 2
            lt_co(1) = -1.d0
            nt_ix1(1, 1) = TYPE_PWL
            nt_ix1(2, 1) = i
            nt_ix1(3, 1) = 2
            nt_ix2(1, 1) = TYPE_U
            nt_ix2(2, 1) = i
            nt_ix2(3, 1) = 0
            nt_co(1) = 1.d0

            call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_U, i, 2,
     &           n, m, sn, vr_reg, bdscl, 1, 1, lt_ix1, nt_ix1, nt_ix2, 
     &           lt_co, nt_co, 0.d0, 0.d0, ws(p_scl+nu+2), ifail)

         end if
      end do

c     ----------------- register objective ------------------------

      do i=1,de
         lt_ix1(1,i) = TYPE_C
         lt_ix1(2,i) = mi+i
         lt_ix1(3,i) = 0
         lt_co(i)    = user(pu_d0+i)
      end do
      
      call reg_cs(iuser, user, sz_iuser, sz_user, CTYPE_OBJ, 0, 0, n, m,
     &     sn, vr_reg, bdscl,
     &     de, 0, lt_ix1, nt_ix1, nt_ix2, lt_co, nt_co, 
     &     0.d0, 0.d0, 0.d0, ifail)
      if (ifail.ne.0) return
c     ------------------- tighten variable bounds -----------------

c      do i=1,n
c         print *, i, bdscl(i), bdscl(n+m+i)
c      end do

c     save old untightend bounds
      do i=1,n+m
         oldbnd(i) = bdscl(i)
         oldbnd(n+m+i) = bdscl(n+m+1+i)
      end do

      if (msg_rep_cs) call prt_cs(n, sn, iuser, user, vr_reg)

c     dChange
c     create matrix boundMax here
c     position 1:nu is max nutrient in any raw material
c     position nu+1 is sum of product masses , nu+2 max cost in any raw
      do l = 1,nu
         boundMax(l)=0
         do i = 1,rm
            if (user(pu_r+(i-1)*nu+l)>boundMax(l)) then
               boundMax(l) = user(pu_r+(i-1)*nu+l)
               j = i
            end if
         end do
      end do
      boundMax(nu+1)=0
      do k = 1,de
         boundMax(nu+1) = boundMax(nu+1) + user(pu_d0+k)
      end do
      boundMax(nu+2)=0
      do i = 1,rm
         if (user(i)>boundMax(nu+2)) then
            boundMax(nu+2) = user(i)
            j = i
         end if
      end do  
c     end dChange

     
c     dChange
c     send boundMax to these four subroutines
      call tgt_bnd_cs(CTYPE_M, iuser, user, n, sn, vr_reg, bdscl, ifail,
     &     boundMax)
      if (ifail.ne.0) then
         ifail = 2
         return
      end if
      call tgt_bnd_cs(CTYPE_C, iuser, user, n, sn, vr_reg, bdscl, ifail,
     &     boundMax)
      if (ifail.ne.0) then
         ifail = 2
         return
      end if
      call tgt_bnd_cs_smp(CTYPE_NU, iuser, user, n, sn, vr_reg, bdscl,
     &     boundMax)
      call tgt_bnd_cs_smp(CTYPE_U, iuser, user, n, sn, vr_reg, bdscl,
     &     boundMax)
c     end dChange

c     ------------------ count constraint nonzeros ----------------

      n_maxa = n
      n_hess = 0
      do i=1,m-NLD%n_cs
         pcsi = iuser(pi_cs_pi+i)
c         print *,i, pcsi
         n_maxa = n_maxa + iuser(pcsi) + 2*iuser(pcsi+1)
         n_hess = n_hess + iuser(pcsi+1)
      end do
      do i=1,NLD%n_cs
         n_maxa = n_maxa + NLD%n_nz_cs(i)
      end do

      pcsi = iuser(pi_cs_pi+m+1)
      n_hess = n_hess + iuser(pcsi+1)
      


      end

