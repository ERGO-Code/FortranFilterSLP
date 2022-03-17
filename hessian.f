      subroutine hessian(x,n,m, phase, lam, ws, lws, user, iuser,
     &                   l_hess, li_hess, flag)

      implicit none

        INCLUDE 'SLPCOMM.INC'

      integer    phe
      parameter (phe=0)

      integer n,m,phase,lws(0:*),iuser(*)
      integer flag
      double precision x(n),lam(-n+1:m),ws(-phe+1:*),user(*)

      integer l, i, j, k, col

      integer l_hess, li_hess
      integer n_lt, n_nt, p_lti, p_nti, p_ltd, p_ntd, pcsi, pcsu,
     &     ix1, ix2, nx_pos

      include 'pusr.inc'

      lws(0) = n_hess
      col = lws(0)

      if ((lws(0).gt.l_hess).or.(2*lws(0)+1.gt.li_hess)) then
         WRITE(nout,*)' Not enough storage space for hessian'
         WRITE(nout,*) 'STOP!'
         stop
      end if
      l_hess = lws(0)
      li_hess = 2*lws(0)+1

      nx_pos = 1
      
      pcsi = iuser(pi_cs_pi+m+1)
      pcsu = iuser(pi_cs_pu+m+1)

      n_lt  = iuser(pcsi)
      n_nt  = iuser(pcsi+1)
      p_lti = iuser(pcsi+2) 
      p_nti = iuser(pcsi+3) 
      p_ltd = iuser(pcsi+4) 
      p_ntd = iuser(pcsi+5) 

      do i=1,n_nt
         ix1 = iuser(p_nti+2*i-2)
         ix2 = iuser(p_nti+2*i-1)
         lws(nx_pos) = ix1
         lws(col+nx_pos) = ix2
         if (phase.eq.2) then
            ws(nx_pos) = user(p_ntd+i-1)
         else
            ws(nx_pos) = 0.d0
         end if
         nx_pos = nx_pos + 1
      end do

c     -------------------- all other gradients -----------------------

      do i=1,m
         
         pcsi = iuser(pi_cs_pi+i)
         pcsu = iuser(pi_cs_pu+i)

         n_nt  = iuser(pcsi+1)
         p_nti = iuser(pcsi+3) 
         p_ntd = iuser(pcsi+5) 
         
         do j=1,n_nt
            ix1 = iuser(p_nti+2*j-2)
            ix2 = iuser(p_nti+2*j-1)
            
            lws(nx_pos) = ix1
            lws(col+nx_pos) = ix2
            ws(nx_pos) = -user(p_ntd+j-1)*lam(i)
            nx_pos = nx_pos + 1
         end do

      end do

c     set hessian to 0

c      do j=1,n_maxa
c         ws(j) = 0.
c      end do

      end
