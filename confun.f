      subroutine confun(x, n, m, res, a, la, user, iuser, flag)

      use nonlin_types
      implicit none

        INCLUDE 'SLPCOMM.INC'

      integer n, m, iuser(*), la(0:*)
      double precision x(n), res(m), a(*), user(*)
      integer flag

      integer i, j
      integer n_lt, n_nt, p_lti, p_nti, p_ltd, p_ntd, pcsi, pcsu,
     &     ix1, ix2

      include 'pusr.inc'


c     for all constraints apart of the nonlinear ones
      do i=1,m-NLD%n_cs
         res(i) = 0.d0
         
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
            res(i) = res(i) + user(p_ltd+j-1)*x(ix1)
         end do
         
         do j=1,n_nt
            ix1 = iuser(p_nti+2*j-2)
            ix2 = iuser(p_nti+2*j-1)
            res(i) = res(i) + user(p_ntd+j-1)*x(ix1)*x(ix2)
         end do
         
      end do

      do i=1,NLD%n_cs
         res(m-NLD%n_cs+i) = eval_nl_cs(NLD%root(i)%p, x)
      end do

c      do i=1,m
c         print *, i, res(i)
c      end do

      end
