      subroutine filterSLP(n,m,maxa,max_f,mxwk, mxiwk, iprint,
     .            ifail, rifail, rho, x, c, f, blo, bup, s, 
     .            ws, lws, lam, user, iuser, max_iter,
     .            istat, rstat, sv_t_fi, se_mx_lp_it, wr_msg,
     .            status, result, oldbnd)


      implicit none

c     ... declaration of passed parameters -- scalars
      integer n, m, maxa, max_f, mxwk, mxiwk, iprint,
     .        rifail, ifail, max_iter
      double precision    rho, f
      double precision a(11000)
      integer la(0:11000)
      integer status(n+m)

c     ... declaration of passed parameters -- arrays
      integer lws(mxiwk), iuser(*), istat(11)
      double precision  blo(n+m), bup(n+m), x(n), c(m), 
     .                  lam(n+m), s(n+m), ws(mxwk), user(*), rstat(3),
     .                  result(7, n+m), oldbnd(2*n+2*m)
      logical sv_t_fi, se_mx_lp_it, wr_msg


c     ... declaration of internal variables
      double precision rand(1000),f1,f2
      double precision c1(1000), c2(1000), eps, dummy, ems_drand
      integer i, j, k, pjp, mark, nr, seed, flag, mxa
      integer ems_time

      if (maxa+2+m.gt.11000.or.n.gt.1000.or.m.gt.1000) then
         print *, n, m, maxa
         print *,'Some dimension too big'
         stop
      end if
      flag = 0

      print *, ' eps :'
      read *, eps
      
      seed = ems_time()
      dummy = ems_drand(seed)

c      call random_seed(nr)
c      do i=1,nr
c         call system_clock(seed(i))
c         seed(i) = mod(seed(i),30081)
c      end do  
c      call random_seed(nr,seed(1:nr))

      do i=1,n
         rand(i) = ems_drand(0)
      end do

      do i=1,n
         x(i) = max(blo(i),-1000.d0) + 
     .             rand(i)*(min(bup(i),1000.d0)-max(blo(i),-1000.d0))
      end do

      print *,'x-values'
      do i=1,n
         print '(I5,3F10.4)', i,blo(i),x(i),bup(i)
      end do
      
      call gradient(n, m, mxa, x, a, la, maxa, user, iuser, flag)
      pjp = la(0)
      
      print *, la(0)
      print *, (la(i), i=1,mxa)
      print *, (la(pjp+i), i=0,m+1)
      print *,' Objective function'
      call objfun(x, n, f1, user, iuser, flag)
      do i=1,n
         x(i)=x(i)+eps
         call objfun(x, n, f2, user, iuser, flag)
         x(i)=x(i)-eps
         print *,i,' -- ', a(i), (f2 - f1)/eps
      end do 
      
      do i=1,m 
         print *,'---------------------------------------------------'
         print *,'Equation ',i
         call confun(x, n, m, c1, a, la, user, iuser, flag)
         do j=la(pjp+i),la(pjp+i+1)-1     
            k = la(j)
            x(k) = x(k) + eps
            call confun(x, n, m, c2, a, la, user, iuser, flag)
            x(k) = x(k) - eps
            if (abs(a(j)/(c2(i) - c1(i))*eps-1.).ge. 0.05) then
               if (a(j).ne.0.or.c2(i).ne.c1(i)) 
     .                             print *,'+++++++++++++++'
            end if
            print *,'    Var ',k,' -- ', a(j), (c2(i) - c1(i))/eps
         end do
      end do
         
      print *,'----------------------------------------------------'
      print *,' other derivatives '
      do i=1,m
         do j=1,n
            mark = 0
            do k=la(pjp+i),la(pjp+i+1)-1
               if (la(k).eq.j) mark = 1
            end do
            if (mark.eq.0) then
               x(j) = x(j) + eps
               call confun(x, n, m, c2, a, la, user, iuser, flag)
               x(j) = x(j) - eps
               if (c2(i).ne.c1(i)) then
                  print '(A,I3,A,I3,A,F25.15)',
     .                 'Error in Eq.',i,' Var.',j,':',c2(i)-c1(i) 
               end if
            end if
         end do
      end do

      stop

      end 

      
