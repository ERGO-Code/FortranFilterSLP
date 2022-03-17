c     =============================================================
c     rd_prob_dim
c     =============================================================
      subroutine rd_prob_dim(sn, bns, n_bns, lws, n_lws, ws, 
     &     SPD, DAT, ifail)
c     
c     subroutine to analyse data files to obtain
c     
c      - n_mu, n_u, n_nu, n_nuratp, n_nuratm
c      - problem dimensions n_vr, n_cs
c      - bin-network bns

c     iflag = 0: OK
c     iflag = 1: Error, global_ifail set

c     lws needs size  rm+mi, 2*n_rminc + rm + mi + de
c     ws  needs size  2*(rm+mi) 

      use interface_mod
      use types_module
      use nonlin_types
      implicit none

        INCLUDE 'SLPCOMM.INC'
        include 'pusr.inc'
        include 'msg.inc'
        include 'return.inc'

      integer sn, n_bns, n_lws, ifail
      integer bns(n_bns), lws(n_lws)
      double precision ws(*)

cFIXME: size of fifo assumes mi,rm < MAX_RMS
      integer fifo(2*MAX_RMS)

      integer i, j, k, i1, i2, i3, fifo1, fifo2
      integer srb, srbb, srbr, nrst, p_nx, p_rmi_fg
c     ... variables for reading SPECIAL file
      type(special_da):: SPD
      type(integra_da):: DAT
c      type(nl_info):: NLD

      integer n_special, n_nl_cs, n_term, n_term_tt
      integer sp_sp1, sp_sp2
      character sp_type
      character(180) line
      integer err
      double precision sp_coeff

      integer trim77

      logical file_exists


c     ... search for dim.dat file

c      nu = NUMBER_NUTRIENTS
c      rm = NUMBER_RAWMATERIALS
c      de = NUMBER_SPECIFICATIONS
c      mi = NUMBER_PREMIXES

      ifail = 0

      if (rm.gt.MAX_RMS) then
         if (msg_err) then
            print *,'too many raw materials'
         end if
chm
         write(global_err_msg, '(A)')
     &     'Number of raw materials (xxxxxx) '//
     &     'bigger than max allowed (xxxxxx)'
chm         write(global_err_msg, '(A,I,A,I,A)') 
chm     &        'Number of raw materials (',rm,
chm     &        ') bigger than max allowed (',MAX_RMS,')'
         write(global_err_msg(25:30), '(I6)') rm
         write(global_err_msg(58:63), '(I6)') MAX_RMS
         global_ifail = 30
         ifail = 1
         goto 99
      end if

c     ------------------------------------------------------------------
c      Analyse Bin-network
c     ------------------------------------------------------------------

      inquire(file=INGRID_FILE,exist=file_exists)
      if (file_exists) then
c     --------- count total number of allowed inclusions: n_mu ---------

         open(10,file=INGRID_FILE)
         n_mu = 0
         do j=1,mi+de
            read (10,*) (ws(i), i=1,rm+mi) 
            read (10,*) (ws(rm+mi+i), i=1,rm+mi) 
            do i=1,rm+mi
               if (max(abs(ws(i)),abs(ws(rm+mi+i))).gt.1d-8) then
                  n_mu=n_mu+1
               end if
            end do
         end do
         close(10)
         
c     Pointers into bns
c     ... map of sources for each bin/demand
         bn_n_srb   = 0
         bn_n_srbb  = bn_n_srb   + mi+de
         bn_p_srb   = bn_n_srbb  + mi+de
         bn_l_srb   = bn_p_srb   + mi+de
c     ... map of targets for each rm/bin
         bn_n_tgb   = bn_l_srb   + n_mu
         bn_n_tgbb  = bn_n_tgb   + rm+mi
         bn_p_tgb   = bn_n_tgbb  + rm+mi
         bn_l_tgb   = bn_p_tgb   + rm+mi
c     ... list of first order bins
         bn_fob     = bn_l_tgb   + n_mu
c                                +mi+de
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
c                               NOT USED: use bn_n_srbb instead
         if (bn_fob+mi+de.gt.n_bns) then
            if (msg_err) then
               print *,'Out of workspace for bin-structure array'
               print *, bn_fob+mi+de, n_bns
            end if
            write(global_err_msg, '(A)')
     &           'Out of workspace in rd_prob_dim: Need xxxxxxxxxx'//
     &           ' ,got xxxxxxxxxx'
chm            write(global_err_msg, '(A,I,A,I)')
chm     &           'Out of workspace in rd_prob_dim: Need ',bn_fob+mi+de,
chm     &           ' ,got ',n_bns
            write(global_err_msg(38:47), '(I10)') bn_fob+mi+de
            write(global_err_msg(54:63), '(I10)') n_bns
            ifail = 1
            global_ifail = 20
            goto 99
         end if

         if (msg_log) then
            WRITE(nout,*) ' reading bounds on bin inclusions: ',
     &           INGRID_FILE(1:trim77(INGRID_FILE)),
     &           ' (analyse bin network)'
         end if
         open(10,file=INGRID_FILE)
         p_nx = 0
         do i=1,rm+mi
            bns(bn_n_tgb+i) = 0
            bns(bn_n_tgbb+i) = 0
         end do
         do i=1,mi+de
            bns(bn_n_srb+i) = 0
            bns(bn_n_srbb+i) = 0
         end do

         do j=1,mi+de
c           ... read lower bounds
            read(10,*) (ws(i), i=1,rm+mi)
c           ... read upper bounds
            read(10,*) (ws(rm+mi+i), i=1,rm+mi)
            bns(bn_p_srb+j) = bn_l_srb + p_nx 

            do i=1,rm+mi

               if (max(abs(ws(rm+mi+i)),abs(ws(i))).gt.1d-8) 
     &              then
c                  print *,'found link:',i, j
                  bns(bn_n_srb+j) = bns(bn_n_srb+j) + 1
                  bns(bn_n_tgb+i) = bns(bn_n_tgb+i) + 1
                  if (i.gt.rm) bns(bn_n_srbb+j) = bns(bn_n_srbb+j)+1
                  if (j.le.mi) bns(bn_n_tgbb+i) = bns(bn_n_tgbb+i)+1

                  p_nx = p_nx+1
                  bns(bn_l_srb+p_nx) = i
               end if
            end do
         end do
         close(10)
      else
c     ... no INGRID.DAT exists: take default
         print *,INGRID_FILE(1:trim77(INGRID_FILE))//
     &        '  does not exists: No predefined action'
         stop
      end if

c     ... convert source table to target table:
c     ... set pointers to list
      bns(bn_p_tgb+1) = bn_l_tgb
      lws(1) = bn_l_tgb+1
      do i=1,rm+mi-1
         bns(bn_p_tgb+i+1) = bns(bn_p_tgb+i) + bns(bn_n_tgb+i)
         lws(i+1) = bns(bn_p_tgb+i+1)+1
      end do
c      do i=1,rm+mi
c         print *, i, lws(i)
c      end do
c     ... fill list
      do i=1,mi+de
         do j=bns(bn_p_srb+i)+1,bns(bn_p_srb+i)+bns(bn_n_srb+i)
            bns(lws(bns(j))) = i
            lws(bns(j)) = lws(bns(j)) + 1
         end do
      end do
c      do i=1,rm+mi
c         print *, i, lws(i)
c      end do

c       call prt_bn(bns, n_bns)
c     ------------------------------------------------------------------
c      END: Analyse Bin-network
c     ------------------------------------------------------------------

      call read_special(SPD, bns, n_bns, ifail)
c     ... any error just bail out, global_ifail/err_msg has been set
      if (ifail.ne.0) goto 99

      call get_n_nl_cs(n_nl_cs)

c     ------------------------------------------------------------------
c      count number of piecewise linear raw material costs: N_PWL
c     ------------------------------------------------------------------

      allocate(DAT%is_piecewise_price(rm), 
     &     DAT%price1(rm), DAT%price2(rm), DAT%rmava1lb(rm+mi), 
     &     DAT%rmava1ub(rm+mi), DAT%rmava2lb(rm), DAT%rmava2ub(rm)) 

      N_PWL = 0
      open(10,file=PRM_FILE)
      do i=1,rm
         DAT%is_piecewise_price(i) = .false.
         read(10,'(a)') line
         read(line,*, iostat=err) DAT%price1(i), DAT%price2(i)
c        .. if there is only one record: err = -1
c           if there are two records then err = 0

         if (err.eq.0.and.DAT%price2(i).lt.500000) then
            DAT%is_piecewise_price(i) = .true.
            N_PWL = N_PWL + 1
         end if
      end do
      close(10)

c     ------------------------------------------------------------------
c      count number of implied rmava constraints: N_U 
c     ------------------------------------------------------------------
c     Every rm with an availability constraints implies an
c     availability constraint along its targets tree
c
c     fifo implements fifo stack: 
c        fifo1 = start of stack
c        fifo2 = next free element of stack
c     elements are added to end and taken from start

c     lws(i) = flag if rm/mi needs a rmava constraint
      do i=1,rm+mi
         lws(i) = 0
      end do
      fifo1 = 1
      fifo2 = 1
      inquire(file=RMLMTS_FILE,exist=file_exists)
      N_U = 0
      if (file_exists) then
         open(10,file=RMLMTS_FILE)
         read (10,*) nrst
c     read all rm with rmava constraints, mark in lws(*) and put on stck
         do i=1,nrst
            read (10,*) i1
            lws(i1) = 1
            fifo(fifo2) = i1
            fifo2 = fifo2 + 1
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c     ... scan special file for any U variables mentioned
         
         do i=1,SPD%n_terms_tt
            if (SPD%vtype(i).eq.TYPE_U) then
               if (lws(SPD%spec1(i)).eq.0) then
                  lws(SPD%spec1(i)) = 1
                  fifo(fifo2) = SPD%spec1(i)
                  fifo2 = fifo2 + 1
               end if
            end if
         end do
         
c         do i=fifo1,fifo2-1
c            print *, i, fifo(i)
c         end do

c     ... flag all bins/rm that need implied rmava constraints 
 10      continue
         i1 = fifo(fifo1)
         i2 = bns(bn_p_tgb+i1)
         fifo1 = fifo1 + 1
         do i=1,bns(bn_n_tgb+i1)
            i3 = bns(i2+i)+rm
            if (i3.le.rm+mi) then
               if (lws(i3).eq.0) then
                  lws(i3) = 1
                  fifo(fifo2) = i3
                  fifo2 = fifo2 + 1
               end if
            end if
         end do
         if (fifo1.lt.fifo2) goto 10

         do i=1,rm+mi
            if (lws(i).eq.1) then
               N_U = N_U + 1
c               print *, i
            end if
         end do
      end if

c     ------------------------------------------------------------------
c           count number of implied rm-inc constraints: N_NU 
c     ------------------------------------------------------------------
c     each rminc c/s in a demand implies a rminc c/s along its 
c     source tree. (NOT NEEDED for bins that are directly fed from
c     the rm in question - first-order-bins (fob))
c     ADDENDUM: the concept of first-order-bin does depend on the rm
c               in question: NOT NEEDED for bins that are only fed 
c               directly from the rm in question. This is now checked 
c               directly
c     
c     this analysis needs to be done for every rm

      nrst = 0
      p_rmi_fg = rm
      do i=1,rm
         lws(i) = 0
      end do
      N_NU = 0
      inquire(file=RMINC_FILE,exist=file_exists)
      if (file_exists) then
c        ... read file into temporary arrays
         open(10,file=RMINC_FILE)
         read(10,*) nrst

         if (2*nrst+rm+mi+de.gt.n_lws) then
            if (msg_err) then
               print *,'Not enough integer workspace in rd_prob_di' 
            end if
            write(global_err_msg, '(A)')
     &           'Out of workspace in rd_prob_dim: Need xxxxxxxxxx'//
     &           ' ,got xxxxxxxxxx'
chm            write(global_err_msg, '(A,I,A,I)')
chm     &           'Out of workspace in rd_prob_dim: Need ',
chm     &           2*nrst+rm+mi+de,' ,got ',n_lws
            write(global_err_msg(38:47), '(I10)') 2*nrst+rm+mi+de
            write(global_err_msg(54:63), '(I10)') n_lws
            ifail = 1
            global_ifail = 20
            goto 99
         end if

c     map in lws:
c     lws(1)        [nrst]     list of rminc constraints: de
c     lws(1+nrst)   [nrst]                "               rm
c     lws(1+2*nrst) [rm]       #of rminc affecting rm-i
c     lws(p_rmi_fg) [mi+de]    flag if mi/de need rminc c/s for this rm

         n_rminc = nrst

c     - - - - - - -  read list of rminc & count no per rm  - - - - - - 

         do i=1,rm
            lws(2*nrst+i) = 0
         end do
         p_rmi_fg = 2*nrst+rm
         do i=1,nrst
            read(10,*) lws(i), lws(nrst+i)
c        ... count no of constraints for given ram-mat
            lws(2*nrst+lws(nrst+i)) = lws(2*nrst+lws(nrst+i))+1
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c     - - - - - - -  Analyse source tree for all affected rm  - - - -
c        ... for all raw_materials that have constraints
         do i=1,rm
c        ... get list of all bin/de affected by particular rm
c            flag them and put them on the stack
            do j=1,mi+de
               lws(p_rmi_fg+j) = 0
            end do
            fifo1 = 1
            fifo2 = 1
            do j=1,nrst
               if (lws(nrst+j).eq.i) then
                  k = lws(j)
                  lws(p_rmi_fg+k) = 1
                  fifo(fifo2) = k
                  fifo2 = fifo2 + 1
               end if
            end do
c     ... scan special file for any NU variables for this rm
            do j=1,SPD%n_terms_tt
               if ((SPD%vtype(j).eq.TYPE_NU).and.
     &              SPD%spec1(j).eq.i) then
                  i1 = SPD%spec2(j)
                  if (lws(p_rmi_fg+i1).eq.0) then
                     lws(p_rmi_fg+i1) = 1
                     fifo(fifo2) = i1
                     fifo2 = fifo2 + 1
                  end if
               end if
            end do
            
c           ... and also include all bin-sources
            do
               if (fifo1.ge.fifo2) exit
               i1 = fifo(fifo1)
               i2 = bns(bn_p_srb+i1)
               srb = bns(bn_n_srb+i1)
               srbb = bns(bn_n_srbb+i1)
               srbr = srb-srbb
               fifo1 = fifo1 + 1
               do j=srbr+1,srb
                  i3 = bns(i2+j)-rm
                  if (lws(p_rmi_fg+i3).eq.0) then
                     lws(p_rmi_fg+i3) = 1
                     fifo(fifo2) = i3
                     fifo2 = fifo2 + 1
                  end if
               end do
            end do
c     ...            count
            do j=1,mi+de
c              ... if this bin/prod is marked then setup NU variable
               if (lws(p_rmi_fg+j).eq.1) then
                  if (bns(bn_n_srbb+j).ne.0) then
c                     print *,'Nu in ',i,j
                     N_NU = N_NU +1
                  end if
               end if
            end do
         end do
      end if

c     ------------------------------------------------------------------
c     ... search for nutrient ratio file (products): 'nuratp.dat'
c     ------------------------------------------------------------------
      n_nuratp = 0
      inquire(file=NURATP_FILE,exist=file_exists)
      if (file_exists) then
         open(31,file=NURATP_FILE)
         read(31,*) n_nuratp
         close(31)
      end if

c     ------------------------------------------------------------------
c     ... search for nutrient ratio file (premixes): 'nuratm.dat'
c     ------------------------------------------------------------------
      n_nuratm = 0
      inquire(file=NURATM_FILE,exist=file_exists)
      if (file_exists) then
         open(31,file=NURATM_FILE)
         read(31,*) n_nuratm
         close(31)
      end if

c     ------------------------------------------------------------------
c      SET: n_vr, n_cs, sn 
c     ------------------------------------------------------------------

c                      M  C    U     MU      NU
      n_vr = (mi+de)*(nu+1) + n_u + n_mu + n_nu
      n_vr = n_vr + 5*n_pwl
      sn = int(sqrt(real(n_vr)))+1
c                      M C S    U     NU    NURAT
      n_cs = (mi+de)*(nu+1+1)+ n_u + n_nu + n_nuratp+n_nuratm
     &     + SPD%n_cs + n_nl_cs
      n_cs = n_cs + 4*n_pwl

c      print *, 'N_U = ',n_u
c      print *, 'N_NU = ',n_nu
c      print *, 'N_NURATP = ',n_nuratp
c      print *, 'N_NURATM = ',n_nuratm
c      print *, 'N_SPECIAL = ',SPD%n_cs
c      print *, 'N_PWL =',n_pwl
c      print *, 'n_vr =',n_vr
c      stop

 99   continue

      end

c     =============================================================
c     rd_prob_dim2
c     =============================================================
      subroutine rd_prob_dim2(filename, ifail)
c   
c     Currently assumes that nu, rm, mi, de are already set in SLPCOMMON
c     and can just be read out.
c     dim.dat is currently read in diet_drv (and copied into SLPCOMMON)
c     
c     subroutine to analyse data files to obtain
c     
c      - n_mu, n_u, n_nu, n_nuratp, n_nuratm
c      - problem dimensions n_vr, n_cs
c      - bin-network bns

c     ifail returns
c     0 ok
c     1 inconsistent data files
c     3 out of memory
c     4 exceeding hard wired limit

c     lws needs size  rm+mi, 2*n_rminc + rm + mi + de
c     ws  needs size  2*(rm+mi) 

      use interface_mod
      use types_module
      use nonlin_types
      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'pusr.inc'
      include 'msg.inc'
      include 'return.inc'

      character(*) filename
      integer ifail
      type(special_da):: SPD
      type(integra_da):: DAT
  
      integer sn, n_bns, n_lws, iflag
      integer, allocatable:: bns(:)

      integer, allocatable :: lws(:)
      double precision, allocatable :: ws(:)

cFIXME: size of fifo assumes mi,rm < MAX_RMS
      integer fifo(2*MAX_RMS)

      integer i, j, k, i1, i2, i3, fifo1, fifo2
      integer srb, srbb, srbr, nrst, p_nx, p_rmi_fg
c     ... variables for reading SPECIAL file
c      type(nl_info):: NLD

      integer n_special, n_nl_cs, n_term, n_term_tt
      integer sp_sp1, sp_sp2
      character sp_type
      character(180) line
      double precision sp_coeff
      integer err

      integer trim77

      logical file_exists

      
c     ... search for dim.dat file

c      nu = NUMBER_NUTRIENTS
c      rm = NUMBER_RAWMATERIALS
c      de = NUMBER_SPECIFICATIONS
c      mi = NUMBER_PREMIXES

      ifail = 0

      if (rm.gt.MAX_RMS) then
         if (msg_err) then
            print *,'too many raw materials'
         end if
chm
         write(global_err_msg, '(A)')
     &     'Number of raw materials (xxxxxx) '//
     &     'bigger than max allowed (xxxxxx)'
chm         write(global_err_msg, '(A,I,A,I,A)') 
chm     &        'Number of raw materials (',rm,
chm     &        ') bigger than max allowed (',MAX_RMS,')'

         global_ifail = 30
         write(global_err_msg(25:30), '(I6)') rm
         write(global_err_msg(58:63), '(I6)') MAX_RMS
         ifail = 1
         goto 99
      end if

c     ------------------------------------------------------------------
c      Analyse Bin-network
c     ------------------------------------------------------------------

      inquire(file=INGRID_FILE,exist=file_exists)
      if (file_exists) then
c     --------- count total number of allowed inclusions: n_mu ---------
         allocate(ws(2*(rm+mi)), stat=err)
         if (err.ne.0) then
            ifail = 1
            write(global_err_msg, '(A)') 
     &           'Allocation failed in rdprobdi.f'
            global_int_arg(1) = 2*(rm+mi)
            global_ifail = 21
            goto 99
         end if
         
         open(10,file=INGRID_FILE)
         n_mu = 0
         do j=1,mi+de
            read (10,*) (ws(i), i=1,rm+mi) 
            read (10,*) (ws(rm+mi+i), i=1,rm+mi) 
            do i=1,rm+mi
               if (max(abs(ws(i)),abs(ws(rm+mi+i))).gt.1d-8) then
                  n_mu=n_mu+1
               end if
            end do
         end do
         close(10)
         
         deallocate(ws)

c     Pointers into bns
c     ... map of sources for each bin/demand
         bn_n_srb   = 0
         bn_n_srbb  = bn_n_srb   + mi+de
         bn_p_srb   = bn_n_srbb  + mi+de
         bn_l_srb   = bn_p_srb   + mi+de
c     ... map of targets for each rm/bin
         bn_n_tgb   = bn_l_srb   + n_mu
         bn_n_tgbb  = bn_n_tgb   + rm+mi
         bn_p_tgb   = bn_n_tgbb  + rm+mi
         bn_l_tgb   = bn_p_tgb   + rm+mi
c     ... list of first order bins
         bn_fob     = bn_l_tgb   + n_mu
c                                +mi+de
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
c                               NOT USED: use bn_n_srbb instead

         n_bns = bn_fob+mi+de
         allocate(bns(bn_fob+mi+de), stat=err)
         if (err.ne.0) then
            write(global_err_msg, '(A,I8)')
     &           'Failure to allocate memory: Need', bn_fob+mi+de
            global_ifail = 21
            global_int_arg(1) = bn_fob+mi+de
            ifail = 1
            goto 99
         end if

         if (msg_log) then
            WRITE(nout,*) ' reading bounds on bin inclusions: ',
     &           INGRID_FILE(1:trim77(INGRID_FILE)),
     &           ' (analyse bin network)'
         end if
         open(10,file=INGRID_FILE)
         p_nx = 0
         do i=1,rm+mi
            bns(bn_n_tgb+i) = 0
            bns(bn_n_tgbb+i) = 0
         end do
         do i=1,mi+de
            bns(bn_n_srb+i) = 0
            bns(bn_n_srbb+i) = 0
         end do

         allocate(ws(2*(rm+mi)), stat=err)
         if (err.ne.0) then
            write(global_err_msg, '(A,I8)')
     &           'Failure to allocate in rdprobdi: Need', 2*(rm+mi)
            global_ifail = 21
            global_int_arg(1) = 2*(rm+mi)
            ifail = 1
            goto 99
         end if

         do j=1,mi+de
c           ... read lower bounds
            read(10,*) (ws(i), i=1,rm+mi)
c           ... read upper bounds
            read(10,*) (ws(rm+mi+i), i=1,rm+mi)
            bns(bn_p_srb+j) = bn_l_srb + p_nx 

            do i=1,rm+mi

               if (max(abs(ws(rm+mi+i)),abs(ws(i))).gt.1d-8) 
     &              then
c                  print *,'found link:',i, j
                  bns(bn_n_srb+j) = bns(bn_n_srb+j) + 1
                  bns(bn_n_tgb+i) = bns(bn_n_tgb+i) + 1
                  if (i.gt.rm) bns(bn_n_srbb+j) = bns(bn_n_srbb+j)+1
                  if (j.le.mi) bns(bn_n_tgbb+i) = bns(bn_n_tgbb+i)+1

                  p_nx = p_nx+1
                  bns(bn_l_srb+p_nx) = i
               end if
            end do
         end do
         close(10)

         deallocate(ws)

      else
c     ... no INGRID.DAT exists: take default
         print *,INGRID_FILE(1:trim77(INGRID_FILE))//
     &        '  does not exists: No predefined action'
         stop
      end if


      allocate(lws(rm+mi), stat=err)
      if (err.ne.0) then
         write(global_err_msg, '(A,I8)')
     &        'Failure to allocate in rdprobdi: Need', rm+mi
         global_ifail = 21
         global_int_arg(1) = rm+mi
         ifail = 1
         goto 99
      end if

c     ... convert source table to target table:
c     ... set pointers to list
      bns(bn_p_tgb+1) = bn_l_tgb
      lws(1) = bn_l_tgb+1
      do i=1,rm+mi-1
         bns(bn_p_tgb+i+1) = bns(bn_p_tgb+i) + bns(bn_n_tgb+i)
         lws(i+1) = bns(bn_p_tgb+i+1)+1
      end do
c      do i=1,rm+mi
c         print *, i, lws(i)
c      end do
c     ... fill list
      do i=1,mi+de
         do j=bns(bn_p_srb+i)+1,bns(bn_p_srb+i)+bns(bn_n_srb+i)
            bns(lws(bns(j))) = i
            lws(bns(j)) = lws(bns(j)) + 1
         end do
      end do
c      do i=1,rm+mi
c         print *, i, lws(i)
c      end do

c       call prt_bn(bns, n_bns)
c     ------------------------------------------------------------------
c      END: Analyse Bin-network
c     ------------------------------------------------------------------

      call read_special(SPD, bns, n_bns, iflag)
      if (iflag.ne.0) goto 99

      call get_n_nl_cs(n_nl_cs)

c     ------------------------------------------------------------------
c      count number of piecewise linear raw material costs: N_PWL
c     ------------------------------------------------------------------

      allocate(DAT%is_piecewise_price(rm), 
     &     DAT%price1(rm), DAT%price2(rm), DAT%rmava1lb(rm+mi), 
     &     DAT%rmava1ub(rm+mi), DAT%rmava2lb(rm), DAT%rmava2ub(rm)) 

      N_PWL = 0
      open(10,file=PRM_FILE)
      do i=1,rm
         DAT%is_piecewise_price(i) = .false.
         read(10,'(a)') line
         read(line,*, iostat=err) DAT%price1(i), DAT%price2(i)
c        .. if there is only one record: err = -1
c           if there are two records then err = 0

         if (err.eq.0) then
            DAT%is_piecewise_price(i) = .true.
            N_PWL = N_PWL + 1
         end if
      end do
      close(10)

c     ------------------------------------------------------------------
c      count number of implied rmava constraints: N_U 
c     ------------------------------------------------------------------
c     Every rm with an availability constraints implies an
c     availability constraint along its targets tree
c
c     fifo implements fifo stack: 
c        fifo1 = start of stack
c        fifo2 = next free element of stack
c     elements are added to end and taken from start

c     lws(i) = flag if rm/mi needs a rmava constraint
      do i=1,rm+mi
         lws(i) = 0
      end do
      fifo1 = 1
      fifo2 = 1
      inquire(file=RMLMTS_FILE,exist=file_exists)
      N_U = 0
      if (file_exists) then
         open(10,file=RMLMTS_FILE)
         read (10,*) nrst
c     read all rm with rmava constraints, mark in lws(*) and put on stck
         do i=1,nrst
            read (10,*) i1
            lws(i1) = 1
            fifo(fifo2) = i1
            fifo2 = fifo2 + 1
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c     ... scan special file for any U variables mentioned
         
         do i=1,SPD%n_terms_tt
            if (SPD%vtype(i).eq.TYPE_U) then
               if (lws(SPD%spec1(i)).eq.0) then
                  lws(SPD%spec1(i)) = 1
                  fifo(fifo2) = SPD%spec1(i)
                  fifo2 = fifo2 + 1
               end if
            end if
         end do
         
c         do i=fifo1,fifo2-1
c            print *, i, fifo(i)
c         end do

c     ... flag all bins/rm that need implied rmava constraints 
 10      continue
         i1 = fifo(fifo1)
         i2 = bns(bn_p_tgb+i1)
         fifo1 = fifo1 + 1
         do i=1,bns(bn_n_tgb+i1)
            i3 = bns(i2+i)+rm
            if (i3.le.rm+mi) then
               if (lws(i3).eq.0) then
                  lws(i3) = 1
                  fifo(fifo2) = i3
                  fifo2 = fifo2 + 1
               end if
            end if
         end do
         if (fifo1.lt.fifo2) goto 10

         do i=1,rm+mi
            if (lws(i).eq.1) then
               N_U = N_U + 1
c               print *, i
            end if
         end do
      end if

      deallocate(lws)

c     ------------------------------------------------------------------
c           count number of implied rm-inc constraints: N_NU 
c     ------------------------------------------------------------------
c     each rminc c/s in a demand implies a rminc c/s along its 
c     source tree. (NOT NEEDED for bins that are directly fed from
c     the rm in question - first-order-bins (fob))
c     ADDENDUM: the concept of first-order-bin does depend on the rm
c               in question: NOT NEEDED for bins that are only fed 
c               directly from the rm in question. This is now checked 
c               directly
c     
c     this analysis needs to be done for every rm

      nrst = 0
      p_rmi_fg = rm
      N_NU = 0
      inquire(file=RMINC_FILE,exist=file_exists)
      if (file_exists) then
c        ... read file into temporary arrays
         open(10,file=RMINC_FILE)
         read(10,*) nrst

         allocate(lws(2*nrst+rm+mi+de), stat=err)
         if (err.ne.0) then
            write(global_err_msg, '(A,I8)')
     &           'Failure to allocate in rdprobdi: Need', 
     &           2*nrst+rm+mi+de
            global_ifail = 21
            global_int_arg(1) = 2*nrst+rm+mi+de
            ifail = 1
            goto 99
         end if
         do i=1,rm
            lws(i) = 0
         end do
         
c         if (2*nrst+rm+mi+de.gt.n_lws) then
c            if (msg_err) then
c               print *,'Not enough integer workspace in rd_prob_di' 
c            end if
c            write(global_err_msg, '(A)')
c     &           'Out of workspace in rd_prob_dim: Need xxxxxxxxxx'//
c     &           ' ,got xxxxxxxxxx'
cchm            write(global_err_msg, '(A,I,A,I)')
cchm     &           'Out of workspace in rd_prob_dim: Need ',
cchm     &           2*nrst+rm+mi+de,' ,got ',n_lws
c            write(global_err_msg(38:47), '(I10)') 2*nrst+rm+mi+de
c            write(global_err_msg(54:63), '(I10)') n_lws
c            iflag = 1
c            ifail = 3
cchm            global_ifail = 3
c            goto 99
c         end if

c     map in lws:
c     lws(1)        [nrst]     list of rminc constraints: de
c     lws(1+nrst)   [nrst]                "               rm
c     lws(1+2*nrst) [rm]       #of rminc affecting rm-i
c     lws(p_rmi_fg) [mi+de]    flag if mi/de need rminc c/s for this rm

         n_rminc = nrst

c     - - - - - - -  read list of rminc & count no per rm  - - - - - - 

         do i=1,rm
            lws(2*nrst+i) = 0
         end do
         p_rmi_fg = 2*nrst+rm
         do i=1,nrst
            read(10,*) lws(i), lws(nrst+i)
c        ... count no of constraints for given ram-mat
            lws(2*nrst+lws(nrst+i)) = lws(2*nrst+lws(nrst+i))+1
         end do
         close(10)
      end if

      if (file_exists.or.SPD%n.ne.0) then
c     - - - - - - -  Analyse source tree for all affected rm  - - - -
c        ... for all raw_materials that have constraints
         do i=1,rm
c        ... get list of all bin/de affected by particular rm
c            flag them and put them on the stack
            do j=1,mi+de
               lws(p_rmi_fg+j) = 0
            end do
            fifo1 = 1
            fifo2 = 1
            do j=1,nrst
               if (lws(nrst+j).eq.i) then
                  k = lws(j)
                  lws(p_rmi_fg+k) = 1
                  fifo(fifo2) = k
                  fifo2 = fifo2 + 1
               end if
            end do
c     ... scan special file for any NU variables for this rm
            do j=1,SPD%n_terms_tt
               if ((SPD%vtype(j).eq.TYPE_NU).and.
     &              SPD%spec1(j).eq.i) then
                  i1 = SPD%spec2(j)
                  if (lws(p_rmi_fg+i1).eq.0) then
                     lws(p_rmi_fg+i1) = 1
                     fifo(fifo2) = i1
                     fifo2 = fifo2 + 1
                  end if
               end if
            end do
            
c           ... and also include all bin-sources
            do
               if (fifo1.ge.fifo2) exit
               i1 = fifo(fifo1)
               i2 = bns(bn_p_srb+i1)
               srb = bns(bn_n_srb+i1)
               srbb = bns(bn_n_srbb+i1)
               srbr = srb-srbb
               fifo1 = fifo1 + 1
               do j=srbr+1,srb
                  i3 = bns(i2+j)-rm
                  if (lws(p_rmi_fg+i3).eq.0) then
                     lws(p_rmi_fg+i3) = 1
                     fifo(fifo2) = i3
                     fifo2 = fifo2 + 1
                  end if
               end do
            end do
c     ...            count
            do j=1,mi+de
c              ... if this bin/prod is marked then setup NU variable
               if (lws(p_rmi_fg+j).eq.1) then
                  if (bns(bn_n_srbb+j).ne.0) then
c                     print *,'Nu in ',i,j
                     N_NU = N_NU +1
                  end if
               end if
            end do
         end do
      end if

      deallocate(lws, bns)

c     ------------------------------------------------------------------
c     ... search for nutrient ratio file (products): 'nuratp.dat'
c     ------------------------------------------------------------------
      n_nuratp = 0
      inquire(file=NURATP_FILE,exist=file_exists)
      if (file_exists) then
         open(31,file=NURATP_FILE)
         read(31,*) n_nuratp
         close(31)
      end if

c     ------------------------------------------------------------------
c     ... search for nutrient ratio file (premixes): 'nuratm.dat'
c     ------------------------------------------------------------------
      n_nuratm = 0
      inquire(file=NURATM_FILE,exist=file_exists)
      if (file_exists) then
         open(31,file=NURATM_FILE)
         read(31,*) n_nuratm
         close(31)
      end if

c     ------------------------------------------------------------------
c      SET: n_vr, n_cs, sn 
c     ------------------------------------------------------------------

c                      M  C    U     MU      NU
      n_vr = (mi+de)*(nu+1) + n_u + n_mu + n_nu
      n_vr = n_vr + 5*n_pwl
      sn = int(sqrt(real(n_vr)))+1



c                      M C S    U     NU    NURAT
      n_cs = (mi+de)*(nu+1+1)+ n_u + n_nu + n_nuratp+n_nuratm
     &     + SPD%n_cs + n_nl_cs
      n_cs = n_cs + 4*n_pwl
c      print *, 'N_U = ',n_u
c      print *, 'N_NU = ',n_nu
c      print *, 'N_NURATP = ',n_nuratp
c      print *, 'N_NURATM = ',n_nuratm
c      print *, 'N_SPECIAL = ',SPD%n_cs


 99   continue

      end



c     =============================================================
c     tgtnm
c     =============================================================

      character*5 function tgtnm(i)
      
      use interface_mod

      implicit none
      integer i
      character*4 int2char
      character*5 name

c      include 'pusr.inc'

      if (i.le.mi) then
         name(2:5) = int2char(i)
         name(1:3) = ' mi'
      else
         name(2:5) = int2char(i-mi)
         name(1:3) = ' de'
      end if
      
      tgtnm = name

      end

c     =============================================================
c     srcmn
c     =============================================================
      character*5 function srcnm(i)
      
      use interface_mod
      
      implicit none
      integer i
      character*4 int2char
      character*5 name
c      include 'pusr.inc'

      if (i.le.rm) then
         name(2:5) = int2char(i)
         name(1:3) = ' rm'
      else
         name(2:5) = int2char(i-rm)
         name(1:3) = ' mi'
      end if
      
      srcnm = name

      end

c     =============================================================
c     edge_exists
c     =============================================================
      logical function edge_exists(src, tgt, bns)
      implicit none

      integer src, tgt
      integer bns(*)

      integer i, p

      include 'pusr.inc'
      

      edge_exists = .false.
      p = bns(bn_p_tgb+src)
      do i=1,bns(bn_n_tgb+src)
         if (bns(p+i).eq.tgt) edge_exists = .true.
      end do

      end


c     =============================================================
c     prt_bn
c     =============================================================
      subroutine prt_bn(bns, n_bns)

      use interface_mod

      implicit none

      include 'pusr.inc'
      
      integer n_bns
      integer bns(n_bns)

      integer i, j, srb, ll, p
      character*5 srcnm, tgtnm
      
      
c      ll = bn_n_tgb - bn_l_srb
c      print *,'Total links/space in list: ',n_mu, ll
c      do i=0,(ll-mod(ll,10))/10
c         print '(I3,A,10I3)', bn_l_srb+i*10, ':',
c     &        (bns(bn_l_srb+i*10+j), j=1,10)
c      end do

      print *,'Bin networks by sources: bins:',mi,' Demands: ',de
      print *
      print *,'  bin/de     #sources   #bin sources'
      do i=1,mi+de
         print *, tgtnm(i), bns(bn_n_srb+i), bns(bn_n_srbb+i),
     &        bns(bn_p_srb+i)
      end do
      print *
      print *,'Sources for each bin/de'
      do i=1,mi+de
         srb = bns(bn_n_srb+i)
         p = bns(bn_p_srb+i)
         print *,tgtnm(i), ' : ',(srcnm(bns(p+j)), j=1,srb)
      end do
      print *
c      ll = bn_fob - bn_l_tgb
c      print *,'Total links/space in list: ',n_mu, ll
c      do i=0,(ll-mod(ll,10))/10
c         print '(I3,A,10I3)', bn_l_tgb+i*10, ':',
c     &        (bns(bn_l_tgb+i*10+j), j=1,10)
c      end do
c      print *
      print *,'Bin networks by targets: rm:',rm,' bins: ',mi
      print *
      print *,'  rm/bin     #targets   #bin targets'
      do i=1,rm+mi
         print *, srcnm(i), bns(bn_n_tgb+i), bns(bn_n_tgbb+i),
     &        bns(bn_p_tgb+i)
      end do
      print *
      print *,'Targets for each rm/bin'
      do i=1,rm+mi
         srb = bns(bn_n_tgb+i)
         p = bns(bn_p_tgb+i)
         print *,srcnm(i), ' : ',(tgtnm(bns(p+j)), j=1,srb)
      end do
 
      end


c     =============================================================
c     setup_datafiles
c     =============================================================
      subroutine setup_datafiles(filename, ifail)
      
      implicit none 

      INCLUDE 'SLPCOMM.INC'
      include 'msg.inc'
      include 'return.inc'

      character(len=*) filename
      integer ifail

      logical file_exists
      integer pos
      CHARACTER*80     LOGOUTFILE, RESULTSFILE, STRING

      inquire(file=filename, exist=file_exists)
      if (.not.file_exists) then
         ifail = 1
         global_ifail = 50
         write(global_err_msg, '(3A)') 
     &        'Problem definition file ',filename,' not found.'
         return
      end if
      
      open(10, file=filename)

      READ(10,'(A)') STRING  ! Skip the straights cost penalty
*     
      READ(10,'(A)') DIM_FILE ! Problem dimemsions file
      READ(10,'(A)') LMTS_FILE ! Specification NT constraints
      READ(10,'(A)') INGRID_FILE ! Premix RM constraints
      READ(10,'(A)') BRIAN_FILE ! Premix NT constraints
      READ(10,'(A)') RAWMAT_FILE ! Raw material analyses
      READ(10,'(A)') STRAIGHTS_FILE ! Final product added RM constraints
      READ(10,'(A)') TONS_FILE ! Specification tonnes
      READ(10,'(A)') PRM_FILE ! Raw Material prices
      READ(10,'(A)') PRMS_FILE ! Raw material penalty prices
      READ(10,'(A)') MIDELMT_FILE ! Premix inclusion limits
      READ(10,'(A)') RMINC_FILE ! Global RM inclusion limits
      READ(10,'(A)') NURATP_FILE ! Product nutrient ratios
      READ(10,'(A)') NURATM_FILE ! Premix nutrient ratios
*     
      READ(10,'(A)') LOGOUTFILE ! Iteration log file
*     
      READ(10,'(A)') RESULTSFILE ! Results file (machine readable)
*     
      READ(10,'(A)') STRING  ! Skip the solution list file
      READ(10,'(A)') STRING  ! Skip the solution expost file
*     
      RMLMTS_FILE   ='RMLMTS.WRK'
      SPECS_FILE    ='SPECS.WRK'
      SUMMARY_FILE  ='SUMMARY.WRK'
      SPECIAL_FILE  ='SPECIAL.WRK'
*     
*     
c     the following names are read in. Should no more lines exist
c     reading will stop (which is not satisfactory!)

      READ(10,'(A)',END=100) RMLMTS_FILE ! Overall RM limits file
      READ(10,'(A)',END=100) SPECS_FILE ! Formula nutrient specs
      READ(10,'(A)',END=100) SUMMARY_FILE ! Summary file (in debug mode)
      READ(10,'(A)',END=100) SPECIAL_FILE ! Relationships file
*     
 100  CLOSE(10)
*
*---------------------------
*  Define other SLP files  *
*---------------------------
*
      RNDSD_FILE ='RND_SD.DAT'  ! Starting random number seed
      EMSCT_FILE ='CT_VR.DAT'   ! EMSOL control file
      SLPCT_FILE ='SLP_CT.DAT'  ! SLP control file
      CSTRPT_FILE='CST_REPT.DAT' ! The cost repete % file
*
*

c     scan for a directory deliminator in the filename
      pos = scan(filename, '/\', back=.true.)
      if (pos.gt.0) then
         dirname = filename(1:pos-1)
      else
         dirname =''
      end if


      end subroutine
      
