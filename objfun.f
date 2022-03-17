      subroutine objfun(x, n, f, user, iuser, flag)

      implicit none

        INCLUDE 'SLPCOMM.INC'

      integer n, iuser(*)
      integer flag
      double precision f, x(n), user(*)

      integer i, n_lt, n_nt, p_lti, p_nti, p_ltd, p_ntd, pcsi, pcsu,
     &     m, ix1, ix2

      include 'pusr.inc'

      f = 0.d0
      
      m = n_cs+1
      pcsi = iuser(pi_cs_pi+m)
      pcsu = iuser(pi_cs_pu+m)

      if (iuser(pcsi+6).ne.CTYPE_OBJ) then
         print *,'Fatal Error: objective not where expected!'
         stop
      end if

      n_lt  = iuser(pcsi)
      n_nt  = iuser(pcsi+1)
      p_lti = iuser(pcsi+2) 
      p_nti = iuser(pcsi+3) 
      p_ltd = iuser(pcsi+4) 
      p_ntd = iuser(pcsi+5) 

      do i=1,n_lt
         ix1 = iuser(p_lti+i-1)
         f = f + user(p_ltd+i-1)*x(ix1)
      end do

      do i=1,n_nt
         ix1 = iuser(p_nti+2*i-2)
         ix2 = iuser(p_nti+2*i-1)
         f = f + user(p_ntd+i-1)*x(ix1)*x(ix2)
      end do

      end
