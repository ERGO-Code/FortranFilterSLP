      module nonlin_types
      
      implicit none
      save  
c     This makes all values persistent between accesses to the module



c     ----------------------------------------------------------------
c     NL_ATOM is a node in the tree describing a nonlinear constraint
c     an NL_ATOM has a type (*,+,/,exp,VAL,VAR, etc) and can have one
c     or two children. A VAL is a double value. A VAR is a variable. 
c     
      type nl_atom
c     
      character :: atype
c     number of children and pointers to them
      integer :: n_chd
      type(nl_atom), pointer :: chd1, chd2
c     the value for a VALUE type atom
      double precision :: value
c     index of the variable in case of a VAR type atom
      integer :: ix

      end type nl_atom

c     ----------------------------------------------------------------
c     PATOM is a pointer to nl_atom
c     This is needed since FORTRAN90/95 does not otherwise allow
c     arrays of pointers. This workaround suggested in Ellis et al, p576
      type patom
      type(nl_atom), pointer :: p
      end type patom


c     The following two data type set up a linked list of integers. 
c     It is used to keep track of indices of all variables used in
c     a particular constraint. The linked list is used so that 
c     variables can be added dynamically (since every constraint can
c     have a different number of variables and the number of variables
c     is not known before the parsing). Using this we can also generate
c     an ordered list of indices for each constraint
c     ----------------------------------------------------------------
c     INT_LL_EL is an element in a linked list of integers
      type int_ll_el
      integer ix
      type(int_ll_el), pointer :: next
      type(int_ll_el), pointer :: last
      end type int_ll_el

c     ... same trick as above: 
c           so that we can have an array of list heads
      type pint_ll_el
      type(int_ll_el), pointer :: p
      end type pint_ll_el
      
c     ----------------------------------------------------------------
c     NL_INFO stores the information about nonlinear constraints

      type nl_info
      
      integer :: n_cs
      type(patom), allocatable :: root(:)
c     ... the number of variables used in the c/s definition
      integer, allocatable :: n_nz_cs(:)
      type(pint_ll_el), allocatable :: ix_nz_cs(:)
      double precision, pointer :: lb(:), ub(:)

      end type nl_info

c     ------------- Global variables -------------
      type(nl_info) :: NLD
      integer :: cu_cs
      
      contains

      
c     -----------------------------------------------------------------
c     get_n_nl_cs: get number of nonlinear constraints
c     -----------------------------------------------------------------
      subroutine get_n_nl_cs(n_cs)

      implicit none
      INCLUDE 'SLPCOMM.INC'

      integer n_cs
      logical fi_xst

      n_cs = 0
      inquire(file=NONLINEAR_FILE, exist=fi_xst)
      if (.not.fi_xst) then
         return
      end if

      open(10, file=NONLINEAR_FILE)

      read(10,*) n_cs
      close(10)

      end subroutine

c     -----------------------------------------------------------------
c     read_nonlinear: read file nonlinear.dat
c     -----------------------------------------------------------------
c
c     NONLINEAR.DAT provides an interface to pose any general nonlinear
c     constraint involving variables that are present in the model
c     
c
c     its format is as follows
c
c     [no of nonlinear constraints] 
c     [Definition of constraint-1]
c     [lb-1] [ub-1]
c     [Definition of constraint-2]
c     [lb-2] [ub-2]
c     
c     The definition of constraints is given as a normal formula involving
c     +, -, *, /, ^, exp().
c     It can refer to variables in the model as N(i,j), where the meaning
c     of the name and indices is as for the SPECIAL constrains, i.e.
c
c     ID  Type                       spec-1           spec-2
c     R   direct raw-mat inclusion   rm/bin involved  bin/prod involved 
c     T   total raw-mat inclusion    rm               bin/prod
c     N   nutrient content           bin/prod         nu involved 
c     U   total usage of rm/bin      rm/bin           0

      subroutine read_nonlinear(NLD, ifail, n, sn, vr_reg)
      
      implicit none
      INCLUDE 'SLPCOMM.INC'
      include 'pusr.inc'
      include 'msg.inc'
      include 'return.inc'    
      
      type(nl_info) :: NLD
      integer ifail
      integer n, sn, vr_reg(*)

c     n, sn, vr_reg: needed for fd_vr: to find the index of a variable
c                    defined by TYPE and two indices


      logical fi_xst
      integer i, j, err
      character(len=200) line
      character(len=20) str1, str2
      
      ifail = 0
      NLD%n_cs = 0
      inquire(file=NONLINEAR_FILE, exist=fi_xst)
      if (.not.fi_xst) then
         return
      end if

      if (msg_log) then
         print *,'reading general nonlinear constraints: nonlinear.dat'
      end if

      open(10, file=NONLINEAR_FILE)

      read(10,*) NLD%n_cs
c      print *, NLD%n_cs
      allocate(NLD%root(NLD%n_cs),NLD%lb(NLD%n_cs), NLD%ub(NLD%n_cs), 
     &    NLD%n_nz_cs(NLD%n_cs), NLD%ix_nz_cs(NLD%n_cs), STAT = err)
      if (err.ne.0) goto 99
      
      do i=1,NLD%n_cs
         NLD%n_nz_cs(i) = 0
      end do

      do i=1,NLD%n_cs
         read(10,'(A)') line
c         print *, line

c        Read in lower and upper bounds as string and then convert them
c        to double. If the bound is specified as '.', then set to 
c        +/- 1d+30
         read(10,*) str1, str2
c         print *, str1, str2
         if (str1.eq.'.') then
            NLD%lb(i) = -1d31
         else
            read(str1, *) NLD%lb(i)
         end if
         if (str2.eq.'.') then
            NLD%ub(i) = 1d31
         else
            read(str2, *) NLD%ub(i)
         end if
c         print *, NLD%lb(i), NLD%ub(i)

c        parse the line defining the constraint
         cu_cs = i
         ifail = 0
         call parse_line(line, NLD%root(i)%p, n, sn, vr_reg, ifail)
      end do
      close(10)

      return
 99   continue
      
      if (msg_err) then
         write(6,*) 'ERROR allocating memory in read_nonlinear'
      end if
      write(global_err_msg, '(A)') 
     &     'ERROR allocating memory in read_nonlin'
      ifail = 1
      global_ifail = 21
      return

      end subroutine

c     ------------------------------------------------------------------
c     parse_line
c     ------------------------------------------------------------------
c     Recursive routine to parse a line declaring a nonlinear constraint
c     and to set up the corresponding expression tree
c
c     PARAMETERS:
c     IN: 
c      character(*) line: the (remainder) of the string to be parsed
c      integer n 
c      integer sn 
c      integer(*) vr_reg: the store of registered variables. Passed 
c                    on to fd_vr in order to obtain the index of
c                    the variables that is refered to (and to store
c                    the index on the nl_atom
c     IN/OUT
c      nl_atom root: A (as yet unassigned) pointer to a nl_atom.
c                    The routine will create the node (via
c                    declare_atom). On return this will be the
c                    root node to the created expression tree.
c     out
c      integer ifail: failure code
      recursive subroutine parse_line(line, root, n, sn, vr_reg, ifail)

      include 'msg.inc'
      include 'return.inc'

      type(nl_atom), pointer :: root
      character(len=*) line
      integer n, sn, vr_reg(*), ifail
      double precision val

      integer pos, ln, ix, ps, pe, iflag
c     parse line will try to find operators in the line set up the 
c     appropriate atom and call itself recursively on the remainder
      
c     remove matching brackets and set ps,pe to proper start of 
c     expression
      call remove_brackets(line, ps, pe, ifail)
      if (ifail.ne.0) return
c      print *,'return remove_brackets:>'//line//'<'
      if (ps.gt.1.or.pe.lt.len(line)) then
         call parse_line(line(ps:pe), root, n, sn, vr_reg, ifail)
         return
      end if

c     first look for '+' or '-' type clauses
      pos = scan_outside_brackets(line, '+', ifail)
      if (ifail.ne.0) return
      if (pos.ge.0) then
c         print *,'Found +:',line(1:pos-1),line(pos+1:len(line))
         call declare_atom(root, '+', 2, ifail)
         if (ifail.ne.0) return
         call parse_line(line(1:pos-1), root%chd1, n, sn, vr_reg, ifail)
         if (ifail.ne.0) return
         call parse_line(line(pos+1:len(line)), root%chd2, n, sn, 
     &        vr_reg, ifail)
         return
      end if
      pos = scan_outside_brackets(line, '-', ifail)
      if (ifail.ne.0) return
c     ... we do not allow a match for a '-' at start -> this is a sign 
      if (pos.ge.2) then
         call declare_atom(root, '-', 2, ifail)
         if (ifail.ne.0) return
         call parse_line(line(1:pos-1), root%chd1, n, sn, vr_reg, ifail)
         if (ifail.ne.0) return
         call parse_line(line(pos+1:len(line)), root%chd2, n, sn,
     &        vr_reg, ifail)
         return
      end if
      
c     next scan for '*' or '/' clauses
      pos = scan_outside_brackets(line, '*', ifail)
      if (ifail.ne.0) return
      if (pos.ge.0) then
c         print *,'Found *:',line(1:pos-1),line(pos+1:len(line))
         call declare_atom(root, '*', 2, ifail)
         if (ifail.ne.0) return
         call parse_line(line(1:pos-1), root%chd1, n, sn, vr_reg, ifail)
         if (ifail.ne.0) return
c         print *, 'parsed before *'
         call parse_line(line(pos+1:len(line)), root%chd2, n, sn,
     &        vr_reg, ifail)
c         print *, 'parsed after *'
         return
      end if
      pos = scan_outside_brackets(line, '/', ifail)
      if (ifail.ne.0) return
      if (pos.ge.0) then
         call declare_atom(root, '/', 2, ifail)
         if (ifail.ne.0) return
         call parse_line(line(1:pos-1), root%chd1, n, sn, vr_reg, ifail)
         if (ifail.ne.0) return
         call parse_line(line(pos+1:len(line)), root%chd2, n, sn,
     &        vr_reg, ifail)
         return
      end if
c     next scan for '^'
      pos = scan_outside_brackets(line, '^', ifail)
      if (ifail.ne.0) return
      if (pos.ge.0) then
         call declare_atom(root, '^', 2, ifail)
         if (ifail.ne.0) return
         call parse_line(line(1:pos-1), root%chd1, n, sn, vr_reg, ifail)
         if (ifail.ne.0) return
         call parse_line(line(pos+1:len(line)), root%chd2, n, sn,
     &        vr_reg, ifail)
         return
      end if

c     Ok, then it could be a literal number, a variable or a function
c     function and variable only differ by the number of arguments 
c     and the name
      
c      print *,'No operator found:>'//line//'<'
      ln = len(line)
      if (line(ln:ln).eq.')') then
c         print *,'scan for (',line
         pos = index(line,'(')
c         print *, 'test for var:',line
         ix = test_var(line(1:pos-1),line(pos+1:ln-1), n, sn, vr_reg,
     &        ifail)
         if (ifail.ne.0) return
         if (ix>0) then 
c            print *,ix
            call declare_atom(root, 'V',0, ifail)
            if (ifail.ne.0) return
            root%ix = ix
c           ... put this variable on the list of variables for cs
            call put_vr_lst(ix, iflag)
c           this is a variable
         else 
c           this must be a function
            if ((line(1:pos-1).eq.'exp')
     &           .or.(line(1:pos-1).eq.'EXP')) then
               call declare_atom(root, 'e', 1, ifail)
               if (ifail.ne.0) return
               call parse_line(line(pos+1:ln-1), root%chd1, n, sn,
     &              vr_reg, ifail)
               return
            end if
            if (msg_err) then
               print *, 'Parse error: substring ends on ) but not'//
     &              ' a variable or function'
               print *, '>'//line//'<'
            end if
            global_ifail = 50
            write(global_err_msg,'(A)') 'Parse error:'//
     &           ' substring ends on ) but not a variable or function'
            ifail = 1
            return
         end if
      else
c        does not end in ')', I assume this must be a float literal
         read (line,*) val
         call declare_atom(root, '1', 0, ifail)
         if (ifail.ne.0) return
         root%value = val
         return
      end if


      end subroutine


c     ------------------------------------------------------------------
c     declare atom
c     ------------------------------------------------------------------
      subroutine declare_atom(atom, type, nchd, ifail)

      include 'msg.inc'
      include 'return.inc'

      type(nl_atom), pointer :: atom
      character*1 type
      integer nchd, ifail

      integer err

      ifail = 0

      allocate(atom, stat=err)
      if (err.ne.0) then
         if (msg_err)
     &        print *, 'Error allocating memory for nonlinear atom'
         write(global_err_msg,'(A)') 
     &        'Error allocating memory for nonlinear atom'
         global_ifail = 21
         ifail = 1
         return
      end if
      atom%n_chd = nchd
      atom%atype = type
      end subroutine 

c     ------------------------------------------------------------------
c     remove_brackets
c     ------------------------------------------------------------------
      subroutine remove_brackets(line, p1, p2, ifail)
      
      include 'msg.inc'
      include 'return.inc'

      character(len=*) line

      integer p1,p2, pe, ifail

c      print *,'remove_brackets called with:>'//line//'<'
c     remove any trailing and leading blanks
c      line = trim(adjustl(line))

c     p1,p2 are pointer to beginning and end of current section
      p1 = 1
      p2 = len(line)
c     skip any leading blanks
      do while(line(p1:p1).eq.' ')
         p1 = p1+1
      end do
c     skip any trailing blanks
      do while(line(p2:p2).eq.' ')
         p2 = p2-1
      end do

c     need to be careful when removing brackets: only matching brackets 
c     can be removed. Beware of "(1+M(1,2)) - M(2,3)"!
      do 
         if (line(p1:p1).eq.'(') then
c           ... find where the closing bracket is
            pe = find_match_bracket(line, p1, p2)
            if (pe.le.0) then
               if (msg_err)
     &           print *, 'Nonlin parse: Could not find closing bracket'
               ifail = 1
               global_ifail = 50
               write(global_err_msg,'(A)') 
     &              'Nonlin parse: Could not find closing bracket'
               return
            end if
            if (pe.eq.p2) then
               p1 = p1+1
               p2 = p2-1
            else
c              ... no brackets to remove if match is within string
               exit;
            end if
         else
c           ... no brackets to remove if first is not bracket
            exit;
         end if
      end do
c     skip any leading blanks
      do while(line(p1:p1).eq.' ')
         p1 = p1+1
      end do
c     skip any trailing blanks
      do while(line(p2:p2).eq.' ')
         p2 = p2-1
      end do

      end subroutine

c     ------------------------------------------------------------------
c     find_match_bracket
c     ------------------------------------------------------------------
c     scan the expression 'line' for the maching bracket for the
c     opening bracket at p1. If there is no matching bracket return -1
c     do not scan past p2

      integer function find_match_bracket(line, p1, p2)

      implicit none
      character(len=*) line
      integer p1, p2

      integer pe, lvl
      pe = p1+1
      lvl = 1
      find_match_bracket = -1
      do 
         if (line(pe:pe).eq.'(') then
            lvl = lvl +1
         end if
         if (line(pe:pe).eq.')') then
            if (lvl.eq.1) then
c              ... this is the matching bracket
               find_match_bracket = pe
               exit;
            else
               lvl = lvl-1
            end if
         end if
         if (pe.eq.p2) exit;
         pe = pe+1
      end do
      

      end function
c     ------------------------------------------------------------------
c     scan_outside_brackets
c     ------------------------------------------------------------------
      integer function scan_outside_brackets(line, ch, ifail)
      
      implicit none

      include 'msg.inc'
      include 'return.inc'


      character(len=*) line
      character ch
      integer ifail

c     position of symbol found
      integer pos
c     number of levels down in bracketing expressions
      integer brk_lev

      integer i

      pos = -1
      brk_lev = 0

      do i=1,len(line)
c        adjust bracketing level
         if (line(i:i).eq.'(') then
            brk_lev = brk_lev + 1
         end if
         if (line(i:i).eq.')') then
            brk_lev = brk_lev - 1
         end if
         if (brk_lev.lt.0) then
            if (msg_err) then
               print *, 'Nonlin parse:'//
     &              ' Term has more closing than opening brackets'
               print *, line
            end if
            ifail = 1
            global_ifail = 50
            write(global_err_msg,'(A)') 
     &       'Nonlin parse: Term has more closing than opening brackets'
            return
         end if
         if ((line(i:i).eq.ch).and.(brk_lev.eq.0)) then
            pos = i
            exit
         end if
      end do

      scan_outside_brackets = pos 
      
      end function

c     ------------------------------------------------------------------
c     test_var
c     ------------------------------------------------------------------
c
c     This function tests if 'name' is an identifer for a variable
c     if so it look in indices which should be a comma separated list
c     of integers and will identify the exact variable
c     
c     test_var returns the index of the variable or -1
c
      integer function test_var(name, indices, n, sn, vr_reg, ifail)
      character(len=*) name, indices

      character(len=10) str1
      integer pcom, i, j, ix, n, sn, ifail
      integer vr_reg(*)

      include 'pusr.inc'
      include 'msg.inc'
      include 'return.inc'

      test_var = -1
c      print *, 'name =',name
c      print *, 'indices =',indices
      if (name.eq.'N') then
         pcom = index(indices, ',')
         read (indices(1:pcom-1),*) i
         read (indices(pcom+1:len(indices)),*) j
c         print *, i, j
         call fd_vr(n, sn, vr_reg, TYPE_M, i, j, ix)
         if (ix<0) then
            if (msg_err) then
               print '(A,A1,I3,I3)','nonlinear term refers to variable',
     &              'N', i,j
               print *, 'fd_vr cannot find it'
            end if
            ifail = 1
            global_ifail = 99
            write(global_err_msg,'(A)') 
     &           'Nonlin parse: nonlinear variable ref not found'
            return
         end if
         test_var = ix
         return
      end if
      if (name.eq.'R') then
         pcom = index(indices, ',')
         read (indices(1:pcom-1),*) i
         read (indices(pcom+1:len(indices)),*) j
c         print *, i, j
         call fd_vr(n, sn, vr_reg, TYPE_MU, i, j, ix)
         if (ix<0) then
            if (msg_err) then
               print '(A,A1,I3,I3)','nonlinear term refers to variable',
     &              'R', i,j
               print *, 'fd_vr cannot find it'
            end if
            ifail = 1
            global_ifail = 99
            write(global_err_msg,'(A)') 
     &           'Nonlin parse: nonlinear variable ref not found'
            return
         end if
         test_var = ix
         return
      end if
      if (name.eq.'T') then
         pcom = index(indices, ',')
         read (indices(1:pcom-1),*) i
         read (indices(pcom+1:len(indices)),*) j
c         print *, i, j
         call fd_vr(n, sn, vr_reg, TYPE_NU, i, j, ix)
         if (ix<0) then
            if (msg_err) then
               print '(A,A1,I3,I3)','nonlinear term refers to variable',
     &              'T', i,j
               print *, 'fd_vr cannot find it'
            end if
            ifail = 1
            global_ifail = 99
            write(global_err_msg,'(A)') 
     &           'Nonlin parse: nonlinear variable ref not found'
            return
         end if
         test_var = ix
         return
      end if
      if (name.eq.'U') then
         read (indices,*) i
         call fd_vr(n, sn, vr_reg, TYPE_U, i, i, ix)
         if (ix<0) then
            if (msg_err) then
               print '(A,A1,I3)','nonlinear term refers to variable',
     &              'U', i
               print *, 'fd_vr cannot find it'
            end if
            ifail = 1
            global_ifail = 99
            write(global_err_msg,'(A)') 
     &           'Nonlin parse: nonlinear variable ref not found'
            return
         end if
         test_var = ix
         return
      end if

      return

      end function

c     ------------------------------------------------------------------
c     put_vr_lst
c     ------------------------------------------------------------------
c     This subroutine adds index 'ix' to the linked list of variable 
c     indices refered by a constraint. The current constraint is
c     stored in global variables cu_cs
c
      subroutine put_vr_lst(ix, ifail)

      implicit none
      include 'msg.inc'
      include 'return.inc'
      
      integer ix, ifail

      integer err
      type(int_ll_el), pointer :: el, pel

      ifail = 0
c     if this is the first variable, initialise list
      if (NLD%n_nz_cs(cu_cs).eq.0) then
         allocate(el, stat=err)
         if (err.ne.0) goto 99
         nullify(el%next)
         nullify(el%last)
         el%ix = ix
         NLD%ix_nz_cs(cu_cs)%p => el
         NLD%n_nz_cs(cu_cs) = 1
      else
c        traverse list until we find an element with larger index or
c        the end of the list
         pel => NLD%ix_nz_cs(cu_cs)%p
         do while (pel%ix<ix .and.associated(pel%next))
            pel => pel%next
         end do
c        ... found variable again -> do nothing
         if (pel%ix.eq.ix) return

c        ...create new element for linked list
         allocate(el, stat=err)
         if (err.ne.0) goto 99
         el%ix = ix
         NLD%n_nz_cs(cu_cs) = NLD%n_nz_cs(cu_cs)+1
c         print *,'n_nz_cs is now',NLD%n_nz_cs(cu_cs)
            
c        ... first variable already larger -> add new head
         if (.not.associated(pel%last)) then
c        ... could also be that there is only one element
            if (ix<pel%ix) then
               el%next => pel
               nullify(el%last)
               pel%last => el
               NLD%ix_nz_cs(cu_cs)%p => el
            else
c             ... add behind
               nullify(el%next)
               el%last => pel
               pel%next => el
            end if
            return
         end if
c        ... if final node and ix bigger
         if ((.not.associated(pel%next)).and.ix>pel%ix) then
            pel%next => el
            el%last => pel
            nullify(el%next)
            return
         end if
c        ... finally this is a normal node: add before pel
         el%next => pel
         el%last => pel%last
         pel%last%next => el
         pel%last => el
         return
      end if

      return
 99   continue
      
      if (msg_err) then
         write(6,*) 'ERROR allocating memory in put_vr_lst'
      end if
      write(global_err_msg, '(A)') 
     &     'ERROR allocating memory in read_nonlin'
      ifail = 1
      return

c     ------ here follows the error handling in nonlin.f ---------
 199  continue

      

      end subroutine

c     ------------------------------------------------------------------
c     print_nl_cs
c     ------------------------------------------------------------------
c
c     assemble the expression given by the atom root and return it
c     in the function return value
c     FIXED: the bug below is almost fixed. At least for 
c     interaction of '*'/'/' with '+'/'-'
c     BUG: this does not indicate the order of calculations
c     i.e. in bracketed expressions 
c           c*(a+b), 
c     the brackets are discarded during the parse and not inserted 
c     during the print. The expression will thus be printed as
c           c*b+a
c     nevertheless, the calculations are correct, since b+a is
c     evaluated before c*(...)
c     FEATURE request: rather than printing 'var(   ix)' should 
c     attempt to print the proper variable name (not sure how
c     easy this is).
c     

      recursive character(len=200) 
     &     function print_nl_cs(root, op, n, sn, vr_reg) result(res)
      
      type(nl_atom), pointer :: root
      character op
c     these are needed so that we can look up variable type and indices
      integer n, sn, vr_reg(*)

      integer vr_tp
      character op1, op2
      character(len=200) str1, str2
      
      if (root%atype.eq.'*') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
         if (op1.eq.'+'.or.op1.eq.'-') then
            str1 = '('//trim(str1)//')'
         end if
         str2 = print_nl_cs(root%chd2, op2, n, sn, vr_reg)
         if (op2.eq.'+'.or.op2.eq.'-') then
            str2 = '('//trim(str2)//')'
         end if
c         print_nl_cs=trim(str1)//'*'//trim(str2)
         res=trim(str1)//'*'//trim(str2)
         op = '*'
         return
      end if
      if (root%atype.eq.'/') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
         if (op1.eq.'+'.or.op1.eq.'-') then
            str1 = '('//trim(str1)//')'
         end if
         str2 = print_nl_cs(root%chd2, op2, n, sn, vr_reg)
         if (op2.eq.'+'.or.op2.eq.'-') then
            str2 = '('//trim(str2)//')'
         end if
c         print_nl_cs=trim(str1)//'/'//trim(str2)
         res=trim(str1)//'/'//trim(str2)
         op = '/'
         return
      end if
      if (root%atype.eq.'+') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
         str2 = print_nl_cs(root%chd2, op2, n, sn, vr_reg)
c          print_nl_cs=trim(str1)//'+'//trim(str2)
          res=trim(str1)//'+'//trim(str2)
          op = '+'
          return
      end if
      if (root%atype.eq.'-') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
         str2 = print_nl_cs(root%chd2, op2, n, sn, vr_reg)
c          print_nl_cs=trim(str1)//'-'//trim(str2)
          res=trim(str1)//'-'//trim(str2)
          op = '-'
          return
      end if
      if (root%atype.eq.'^') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
         str2 = print_nl_cs(root%chd2, op2, n, sn, vr_reg)
c         print_nl_cs=trim(str1)//'^'//trim(str2)
         res=trim(str1)//'^'//trim(str2)
         op = '^'
         return
      end if
      if (root%atype.eq.'e') then
         str1 = print_nl_cs(root%chd1, op1, n, sn, vr_reg)
c          print_nl_cs='exp('//trim(str1)//')'
          res='exp('//trim(str1)//')'
          op = 'f'
          return
      end if

c     and now treat terminal symbols
      if (root%atype.eq.'1') then
         write(str1,*) root%value
c         print_nl_cs=trim(str1)
         res=trim(str1)
         op = '1'
         return
      end if
      if (root%atype.eq.'V') then
c        if I know the index of a variable I can get the following info
c        vr_reg(2+2*sn    +ix): type
c        vr_reg(2+2*sn  +n+ix): ix1
c        vr_reg(2+2*sn+2*n+ix): ix2
c        also (from pusr.inc): M=1,MU=2,C=3,U=4,NU=5
         

         write(str1,'("  (",I3,",",I3,")")') 
     &        vr_reg(2+2*sn+n+root%ix),vr_reg(2+2*sn+2*n+root%ix)
         vr_tp = vr_reg(2+2*sn+root%ix)
         if (vr_tp.eq.1) str1(1:2) = ' M'
         if (vr_tp.eq.2) str1(1:2) = 'MU'
         if (vr_tp.eq.3) str1(1:2) = ' C'
         if (vr_tp.eq.4) str1(1:2) = ' U'
         if (vr_tp.eq.5) str1(1:2) = 'NU '
c         print_nl_cs=trim(str1)
         res=trim(str1)
         op = 'V'
         return
      end if
   
      print *, 'trying to print unknown atom type: ',root%atype
      stop

      end function


      subroutine print_nl_data(n, sn, vr_reg)
      
      implicit none
      integer n, sn, vr_reg(*)

      integer i
      character ch


      type(int_ll_el), pointer :: el
      
      print *, '------------- nonlinear data -----------------'
      print *, ' There are ',NLD%n_cs,' constraints:'
      do i=1,NLD%n_cs
         print *, 'Constraint: ',i
         print *, 'Variables refered: ',NLD%n_nz_cs(i)
         el => NLD%ix_nz_cs(i)%p
         do
            print *, el%ix
            if (.not.associated(el%next)) exit
            el=>el%next
         end do
         print *, print_nl_cs(NLD%root(i)%p, ch, n, sn, vr_reg)
      end do

      end subroutine

c     ------------------------------------------------------------------
c     eval_nl_cs
c     ------------------------------------------------------------------
c     Evaluate a nonlinear constraint function around a given point
c
      recursive double precision function eval_nl_cs(root, x) 
     &     result(res)
      
      implicit none
      type(nl_atom), pointer:: root
      double precision x(*)
      
      select case(root%atype)
      case ('+')
c         eval_nl_cs = 
         res = 
     &        eval_nl_cs(root%chd1, x) + eval_nl_cs(root%chd2, x)
      case ('*')
c         eval_nl_cs = 
         res = 
     &        eval_nl_cs(root%chd1, x) * eval_nl_cs(root%chd2, x)
      case ('/')
c         eval_nl_cs = 
         res = 
     &        eval_nl_cs(root%chd1, x) / eval_nl_cs(root%chd2, x)
      case ('-')
c         eval_nl_cs = 
         res = 
     &        eval_nl_cs(root%chd1, x) - eval_nl_cs(root%chd2, x)
      case ('^')
c         eval_nl_cs = 
         res = 
     &        eval_nl_cs(root%chd1, x) ** eval_nl_cs(root%chd2, x)
      case ('e')
         res = exp(eval_nl_cs(root%chd1, x))
c         eval_nl_cs = exp(eval_nl_cs(root%chd1, x))
      case ('1')
         res = root%value
c         eval_nl_cs = root%value
      case ('V')
         res = x(root%ix)
c         eval_nl_cs = x(root%ix)
      case default
         print *, 'EVAL_NL: encountered unknown node:',root%atype
         stop
      end select

      end function


c     ------------------------------------------------------------------
c     eval_nl_cs
c     ------------------------------------------------------------------
      recursive subroutine get_gradient(root, n, x, c, g)
      
      type(nl_atom), pointer :: root
      integer n
      double precision x(n), g(n), c

      double precision c1, c2, g1(n), g2(n)

c     calculate value and gradient for up to two children
      if (root%n_chd.ge.1) then
         call get_gradient(root%chd1, n, x, c1, g1)
      end if
      if (root%n_chd.ge.2) then
         call get_gradient(root%chd2, n, x, c2, g2)
      end if

c     and combine them depending on the type of operation
      select case(root%atype)
      case ('+')
         g = g1+g2
         c = c1+c2
      case ('-')
         g = g1-g2
         c = c1-c2
      case ('*')
c         print *, c1, c2
c         print *,'g1= ',g1
c         print *,'g2= ',g2
         g = c1*g2+c2*g1
         c = c1*c2
c         print *, 'c = ',c
c         print *, 'g = ',g
      case ('/')
         g = (c2*g1-c1*g2)/(c2*c2)
         c = c1/c2
      case ('^')
c        g1=0 => d/dx (f1(x)^n) = n*(d/dx f1(x))*f1(x)^(n-1) 
c        d/dx (f1(x)^f2(x)) = f1^f2 *(d/dx f2 *log (f1) + f2*d/dx f1/f1) 
c         g = c2*c1**(c2-1d0)*g1
         g = (c1**c2)*(g2*log(c1)+c2*g1/c1) 
         c = c1**c2
      case ('e')
c        d/dx exp(f1(x)) = g1(x) exp(f1(x))
         c = exp(c1)
         g = c*g1
      case('1')
         g = 0.d0
         c = root%value
      case('V')
         c = x(root%ix)
         g = 0.d0
         g(root%ix) = 1.d0
      end select
      end subroutine

      end module
