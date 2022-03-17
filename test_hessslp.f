c     *****************************************************************
c     **
c     **   Check routine for the Hessian of a function suite
c     **
c     **   This routine now uses the 'sqp_9_june' convention of
c     **   calling filterSQP and the subroutines
c     **
c     *****************************************************************

      subroutine filterSLP(n,m,maxa,max_f, mxwk, mxiwk, iprint,
     .            ifail, rho, x, c, f, blo, bup, s, 
     .            ws, lws, lam, user, iuser, max_iter,
     .            istat, rstat, sv_t_fi, se_mx_lp_it, wr_msg)


      implicit none

c     ... declaration of passed parameters -- scalars
      integer n, m, maxa, max_f, mxwk, mxiwk, iprint, 
     .        ifail, max_iter
      double precision    rho, f
      double precision a(11000)
      integer la(0:11000)

c     ... declaration of passed parameters -- arrays
      integer lws(mxiwk), iuser(*), istat(11)
      double precision blo(n+m), bup(n+m), x(n), c(m), 
     .                 lam(n+m), s(n+m), ws(mxwk), user(*), rstat(3)
      logical sv_t_fi, se_mx_lp_it, wr_msg


c     ... declaration of internal variables

      double precision rand(2000),a1(11000),a2(11000)
      double precision ei(2000), v(2000), eps, c1, c2, dummy, ems_drand

      integer i, j, k, pjp, mark, nr, seed, flag, ii, mxa, l_hess,
     .        li_hess
      integer ems_time

      flag = 0

      print *,' eps: '
      read *, eps

      seed = ems_time()
      dummy = ems_drand(seed)

      
      do i=1,n
         rand(i) = ems_drand(0)
      end do
      
      do i=1,n
         x(i) = max(blo(i),-1000.d0) + 
     .        rand(i)*(min(bup(i),1000.d0)-max(blo(i),-1000.d0))
         ei(i) = 0.d0
      end do
      
      do i=1,n+m
         lam(i) = 0.d0
      end do
      
      
      call gradient(n, m, mxa, x, a1, la, maxa, user, iuser, flag)  
      pjp=la(0)
      
      print *,'Objective Function'
      l_hess = mxwk
      li_hess = mxiwk
      call hessian(x, n, m, 2, lam, ws, lws, user, iuser, l_hess,
     .             li_hess, flag)

      do j=1,n
         x(j) = x(j) + eps
         call gradient(n, m, mxa, x, a2, la, maxa, user, iuser, flag)
         x(j) = x(j) - eps
         ei(j) = 1.
         call gdotx(n, ei, ws, lws, v) 
         ei(j) = 0.
         do k = 1,n
            c1 = v(k)
            c2 = (a2(k)-a1(k))/eps
c            print *,j,k,c1, c2
            if ((((abs(c1).gt.1.e-9).and.(abs(c2).gt.1.e-9))
     .           .and.(abs((c1-c2)/c1).gt..01))
     .         .or.((abs(c1).gt.1.e-9).and.(abs(c2).le.1.e-9))
     .         .or.((abs(c1).le.1.e-9).and.(abs(c2).gt.1.e-9))) then
               print *,'+++++++++++++++++++++++++++++'
               read *
            end if
         end do
      end do
      
      do i=300,m
         print *,'=================================================='
         print *,'Equation ',i
         lam(i+n) = - 1.
         l_hess = mxwk
         li_hess = mxiwk
         call hessian(x, n, m, 1, lam, ws, lws, user, iuser, l_hess,
     .                li_hess, flag)
         do j=la(pjp+i),la(pjp+i+1)-1    
            x(la(j)) = x(la(j)) + eps
            call gradient(n, m, mxa, x, a2, la, maxa, user, iuser,flag)
            x(la(j)) = x(la(j)) - eps
            ei(la(j)) = 1.
            call gdotx(n, ei, ws, lws, v) 
            ei(la(j)) = 0.
            do k=la(pjp+i),la(pjp+i+1)-1
               c1 = v(la(k))
               c2 = (a2(k)-a1(k))/eps
               print *,la(j),la(k),c1, c2
               if ((((abs(c1).gt.1.e-9).and.(abs(c2).gt.1.e-9))
     .           .and.(abs((c1-c2)/c1).gt..01))
     .         .or.((abs(c1).gt.1.e-9).and.(abs(c2).le.1.e-9))
     .         .or.((abs(c1).le.1.e-9).and.(abs(c2).gt.1.e-9))) then
                  print *,'+++++++++++++++++++++++++++'
                  read *
               end if
            end do
         end do
         lam(i+n) = 0.
      end do
      
      end 
      subroutine gdotx(n, x, ws, lws, v)
 
      implicit none
      
      integer             phl, phr, phc, phe, pvs, pcs
      common /hess_scale/ phl, phr, phc, phe, pvs, pcs
      
      integer n, lws(0:*)
      double precision x(n), v(n), ws(-phe+1:*)
      
      integer i, hl, row, col
      
      do i=1,n
         v(i) = 0.
      end do
      
      hl = lws(0)
      do i=1,hl
         row = lws(i)
         col = lws(hl+i)
         v(row) = v(row) + ws(i) * x(col)
         if (row.ne.col) then
            v(col) = v(col) + ws(i) * x(row)
         end if
      end do
      
      end 
