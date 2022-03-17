      module types_module
      
      implicit none
      save

      type special_da
c     description of all constraints
      integer :: n               !total lineas in special.dat
      integer :: n_cs            !  of which constraints
      integer :: n_bnd           !  of which bounds
      integer :: n_terms_tt      !total terms in all bounds & c/s
      integer, pointer :: p_cs(:)
      integer, pointer :: n_terms(:)
      double precision, pointer :: lb(:), ub(:)

c     description of all terms
      integer, pointer :: vtype(:)
      integer, pointer :: spec1(:), spec2(:)
      double precision, pointer :: coeff(:)

      end type special_da

c     integra data is a data type that should store all data read in
c     from the data files. This should eventually replace the user/iuser
c     ws/lws arrays where this information is currently stored

      type integra_da

c     ... indicate if rm has just one price or more components
      logical, pointer :: is_piecewise_price(:)
c     ... prices for raw materials (price2 only where present)
      double precision, pointer :: price1(:), price2(:)
c     ... lower/upper availability for rm at price1
      double precision, pointer :: rmava1lb(:), rmava1ub(:)
c     ... lower/upper availability for rm at price2
      double precision, pointer :: rmava2lb(:), rmava2ub(:)


      end type integra_da

      contains
      
c     -----------------------------------------------------------------
c     read_special: read file special.dat
c     -----------------------------------------------------------------
c
c     SPECIAL.DAT provides an interface to pose any general linear
c     constraint involving variables that are present in the model
c
c     its format is as follows
c
c     [no of special constraints] 
c     [1: no of terms in cs] [1:lower bnd] [1:upper bnd] 
c          followed by a list of all terms in the c/s: 
c     [TYPE] [coeff] [spec-1] [spec-2]
c        : : ; : 
c     [TYPE] [coeff] [spec-1] [spec-2] 
c           followed by specification of special constraint 2, 3, ...
c
c     The TYPE, spec-1 spec-2 refer to the variables as follows
c
c     ID  Type                       spec-1           spec-2
c     R   direct raw-mat inclusion   rm/bin involved  bin/prod involved 
c     T   total raw-mat inclusion    rm               bin/prod
c     N   nutrient content           bin/prod         nu involved 
c     U   total usage of rm/bin      rm/bin           0

      subroutine read_special(SPD, bns, n_bns, ifail)
      
      implicit none
      INCLUDE 'SLPCOMM.INC'
      include 'pusr.inc'
      include 'msg.inc'
      include 'return.inc'

      type(special_da) :: SPD
      integer ifail, n_bns, bns(n_bns)


      integer i, j, k, err, d1, c_cs
      logical fi_xst, fd_src
      character(len=1) :: ch
      double precision d2
      integer n_ign_terms
 
      ifail = 0
      SPD%n_terms_tt = 0
      SPD%n = 0
      SPD%n_cs = 0
      SPD%n_bnd = 0

      inquire(file=SPECIAL_FILE, exist=fi_xst)
      if (.not.fi_xst) then

c     dChange
c     allocate all these with size 0 just to stop error
c     when deallocating them
         allocate(SPD%p_cs(0), SPD%n_terms(0),
     &     SPD%lb(0), SPD%ub(0), STAT = err)
         allocate(SPD%vtype(0), SPD%spec1(0),
     &        SPD%spec2(0), SPD%coeff(0), 
     &        STAT = err)
         if (err.ne.0) goto 99
c     end dChange
         return

      end if

      if (msg_log) then
         print *,'reading general constraints: special.dat'
      end if

      open(10, file=SPECIAL_FILE)

      read(10,*) SPD%n
      allocate(SPD%p_cs(SPD%n+1), SPD%n_terms(SPD%n),
     &     SPD%lb(SPD%n), SPD%ub(SPD%n), STAT = err)
      if (err.ne.0) goto 99

      do i=1,SPD%n
         read(10,*) SPD%n_terms(i) 
         SPD%n_terms_tt = SPD%n_terms_tt + SPD%n_terms(i)
         if (SPD%n_terms(i).eq.1) then
            SPD%n_bnd = SPD%n_bnd + 1
         else
            SPD%n_cs = SPD%n_cs + 1
         end if
         do j=1,SPD%n_terms(i)
            read(10,*)
         end do
      end do
      close(10)
      
      allocate(SPD%vtype(SPD%n_terms_tt), SPD%spec1(SPD%n_terms_tt),
     &     SPD%spec2(SPD%n_terms_tt), SPD%coeff(SPD%n_terms_tt), 
     &     STAT = err)
      if (err.ne.0) goto 99

      
      c_cs = 1
      open(10, file=SPECIAL_FILE)
      read(10, *) d1
      do i=1,SPD%n
         read(10,*) d1, SPD%lb(i), SPD%ub(i)
         SPD%p_cs(i) = c_cs
         n_ign_terms = 0
         do j=1,SPD%n_terms(i)
            read(10,*) ch, SPD%coeff(c_cs), 
     &           SPD%spec1(c_cs), SPD%spec2(c_cs)
            if (SPD%n_terms(i).eq.1.and.SPD%coeff(c_cs).le.1e-10) then
               print *,'ERROR: when setting bounds in SPECIAL.DAT'//
     &              ' coeff has to be > 0'
               stop
c               goto 99
            end if
            select case(ch)
            case ('R','r')
c              .. direct raw-mat/bin inclusion in bin/prod
               SPD%vtype(c_cs) = TYPE_MU
            case ('T','t')
c              ... total raw-mat inclusion in product
               SPD%vtype(c_cs) = TYPE_NU
c               ... for first order bins T is the same as R

c              ... check if no bin sources => Total is equal R
c               print *, n_bns, bn_n_srbb+SPD%spec2(c_cs)
               if (bns(bn_n_srbb+SPD%spec2(c_cs)).eq.0) then
                  if (msg_warn_edge) then
                     print '(A,I3,A,I3,A)',
     &                    'defined NU-type variable for rm ',
     &                    SPD%spec1(c_cs),' in bin/prod ',
     &                    SPD%spec2(c_cs),' which is first order'
                     print '(A)', 'changed to MU-type variable' 
                  end if
                  SPD%vtype(c_cs) = TYPE_MU
c                 ... also check that this edge exists in the graph
c                  print *,'Debug:'
c                  print *,'Bin ',SPD%spec2(c_cs),' has ',
c     &                 bns(bn_n_srb+SPD%spec2(c_cs)),'sources'
                  fd_src = .false.
                  do k=1,bns(bn_n_srb+SPD%spec2(c_cs))
c                     print *, bns(bns(bn_p_srb+SPD%spec2(c_cs))+k)
                     if (bns(bns(bn_p_srb+SPD%spec2(c_cs))+k)
     &                    .eq.SPD%spec1(c_cs)) fd_src = .true.
                  end do
                  if (.not.fd_src) then
                     if (msg_warn_edge) then
                        print '(A,I3,A,I3)',
     &     'Defined Total ingredient constraint on RM',SPD%spec1(c_cs),
     &     'that does not feed into bin',SPD%spec2(c_cs)
                        print *, 'Will ignore this term'
                     end if
                     n_ign_terms = n_ign_terms +1
                     cycle
                  end if
               end if
            case ('N','n')
c              ... nutrient content of unit of bin/prod
               SPD%vtype(c_cs) = TYPE_M
            case ('U','u')
c              ... total useage of raw-mat/bin
               SPD%vtype(c_cs) = TYPE_U
               SPD%spec2(c_cs) = 0
            case default
               write(nout,*) 'ERROR: Unknown constraint type in special'
               stop
c               goto 99
            end select
c            print *,'RDSPECIAL: read record:'
c            print *, c_cs, SPD%vtype(c_cs), SPD%coeff(c_cs), 
c     &           SPD%spec1(c_cs), SPD%spec2(c_cs)
            c_cs = c_cs + 1
         end do
         SPD%n_terms(i) = SPD%n_terms(i) - n_ign_terms
      end do
      close(10)

      return
 99   continue
      
      if (msg_err) then
         write(nout,*) 'ERROR allocating memory in read_special'
      end if
      write(global_err_msg, '(A)') 
     &     'ERROR allocating memory in read_special.'
      ifail = 1
      global_ifail = 21
      return

      end subroutine
      

      subroutine free_special(SPD)

      implicit none
      type(special_da) :: SPD

      deallocate(SPD%p_cs, SPD%n_terms, SPD%lb, SPD%ub, SPD%vtype, 
     &     SPD%spec1, SPD%spec2, SPD%coeff)
      
      end subroutine 

      end module
