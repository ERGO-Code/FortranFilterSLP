      subroutine ran_print(n, m, sn, result, vr_reg, iuser)
      
      implicit none

      include 'pusr.inc'
      include 'SLPCOMM.INC'


      integer n, m, sn
      double precision result(7, n+m)
      integer vr_reg(2+2*sn+3*n), iuser(*)

      integer i, tp, j, k, l, pi

c     ==================== subroutine body =======================

      write(nout,*)
      write(nout,*)
      write(nout,*)  ' ============= RANGING ================='
      write(nout,*)
      write(nout,*)

      write(nout,*) ' Information for variables:'
      write(nout,*) 'I TP J K   RANGING INFORMATION '
      do i=1,n
         tp = vr_reg(2+2*sn+i)
         j = vr_reg(2+2*sn+n+i)
         k = vr_reg(2+2*sn+2*n+i)
         write(nout,'(I3,3I2,7G10.2)') i, tp, j, k, (result(l,i), l=1,7)
      end do

      write(nout,*) ' Information for constraints'
      write(nout,*) 'I TP J K   RANGING INFORMATION '
      do i=1,m
         pi = iuser(pi_cs_pi+i)
         tp = iuser(pi+6)
         j = iuser(pi+7)
         k = iuser(pi+8)
         write(nout,'(I3,3I2,7G10.2)') i, tp, j, k, 
     &        (result(l,n+i), l=1,7)
      end do

      end
