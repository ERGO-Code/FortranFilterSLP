      subroutine gradient(n, m, mxa, x, a, la, maxa, user, iuser, flag)

      use nonlin_types
      implicit none

        INCLUDE 'SLPCOMM.INC'

      integer n, m, mxa, maxa, iuser(*), la(0:*)
      integer flag
      double precision x(n), a(*), user(*)

      integer i, j, k, l, pjp, ll, p_a

      integer n_lt, n_nt, p_lti, p_nti, p_ltd, p_ntd, pcsi, pcsu,
     &     ix1, ix2, nx_pos

c     ... temporary gradient (for the nonlinear constraint)
      double precision c, tmpg(n)
      type(int_ll_el), pointer :: el

      include 'pusr.inc'

c     --------------- gradient of objective -----------------------

      mxa = maxa
      pjp = maxa+1
      la(pjp) = 1
      la(pjp+1) = n+1
      la(0) = pjp
      do i=1,n
         la(i) = i
         a(i) = 0.d0
      end do

      pcsi = iuser(pi_cs_pi+m+1)
      pcsu = iuser(pi_cs_pu+m+1)

      n_lt  = iuser(pcsi)
      n_nt  = iuser(pcsi+1)
      p_lti = iuser(pcsi+2) 
      p_nti = iuser(pcsi+3) 
      p_ltd = iuser(pcsi+4) 
      p_ntd = iuser(pcsi+5) 

      do i=1,n_lt
         ix1 = iuser(p_lti+i-1)
         a(ix1) = a(ix1) + user(p_ltd+i-1)
      end do

      do i=1,n_nt
         ix1 = iuser(p_nti+2*i-2)
         ix2 = iuser(p_nti+2*i-1)
         a(ix1) = a(ix1) + x(ix2)*user(p_ntd+i-1)
         a(ix2) = a(ix2) + x(ix1)*user(p_ntd+i-1)
      end do

c     -------------------- all other gradients -----------------------

      nx_pos = n+1
      
      do i=1,m-NLD%n_cs
         
         pcsi = iuser(pi_cs_pi+i)
         pcsu = iuser(pi_cs_pu+i)

         n_lt  = iuser(pcsi)
         n_nt  = iuser(pcsi+1)
         p_lti = iuser(pcsi+2) 
         p_nti = iuser(pcsi+3) 
         p_ltd = iuser(pcsi+4) 
         p_ntd = iuser(pcsi+5) 
         
         do j=1,n_lt
            ix1 = iuser(p_lti+j-1)
            a(nx_pos) = user(p_ltd+j-1)
            la(nx_pos) = ix1
            nx_pos = nx_pos + 1
         end do
         
         do j=1,n_nt
            ix1 = iuser(p_nti+2*j-2)
            ix2 = iuser(p_nti+2*j-1)
            a(nx_pos) = user(p_ntd+j-1)*x(ix2)
            la(nx_pos) = ix1
            nx_pos = nx_pos + 1
            a(nx_pos) = user(p_ntd+j-1)*x(ix1)
            la(nx_pos) = ix2
            nx_pos = nx_pos + 1
         end do
         la(pjp+i+1) = nx_pos 

      end do

      do i=1,NLD%n_cs
c        we should fill in the next locations after nx_pos in the 
c        la/a arrays. I guess the parser needs to count the number
c        of variables and their position.
         call get_gradient(NLD%root(i)%p, n, x, c, tmpg)
         
         el => NLD%ix_nz_cs(i)%p
         la(nx_pos) = el%ix
         a(nx_pos) = tmpg(el%ix)
         nx_pos = nx_pos + 1
         do j=2,NLD%n_nz_cs(i)
            el => el%next
            la(nx_pos) = el%ix
            a(nx_pos) = tmpg(el%ix)
            nx_pos = nx_pos + 1
         end do
         la(pjp+m-NLD%n_cs+i+1) = nx_pos
      end do


      end
