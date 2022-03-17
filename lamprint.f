      subroutine lam_print(n, sn, m, x, vr_reg, user, iuser, blo, bup, 
     &     lam, f, c, nlws, lws, status)

      use interface_mod
      use nonlin_types
      implicit none

      include 'pusr.inc'
      include 'SLPCOMM.INC'

      integer n, sn, m, nlws
      integer iuser(*), vr_reg(2+2*sn+3*n), lws(nlws), status(n+m)

      double precision f
      double precision x(n), user(*), blo(n+m), bup(n+m), 
     &     lam(n+m), c(m)

c     .... local variables

      integer i, ix, j, k, l, pi
      double precision viol
      
c     ==================== subroutine body =======================

      if (max((mi+de)*(rm+mi), (nu+1)*(mi+de)).gt.nlws) then
         write (nout,*) 'Not enough integer workspace in lamprint'
         stop
      end if

      write(nout,*)
      write(nout,*)
      write(nout,*)  ' ============= SENSITIVITY ================='
      write(nout,*)
      write(nout,*)  ' f = ',f
      write(nout,*)
      write(nout,*) 'Composition of mixers (in terms of raw materials):'
      write(nout,*)

c     ... get mu's
      do i=1,(mi+de)*(rm+mi)
         lws(i) = 0
      end do
      
      do i=1,n
         if (vr_reg(2+2*sn+i).eq.TYPE_MU) then
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
            lws((j-1)*(mi+de)+k) = i
         end if
      end do

      do j=1,mi
         write(nout,*) 'Mixer ',j
         do i=1,rm
            ix = lws((i-1)*(mi+de)+j)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'rm-',i, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
         do i=1,mi
            ix = lws((rm+i-1)*(mi+de)+j)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'mi-',i, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
      end do
      write(nout,*)
      write(nout,*)  'Composition of products (in terms of rm/mi)'
      write(nout,*)
      do j=1,de
         write(nout,*) 'Product ',j
         do i=1,rm
            ix = lws((i-1)*(mi+de)+mi+j)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'rm-',i, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
         do i=1,mi
            ix = lws((i+rm-1)*(mi+de)+mi+j)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'mi-',i, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
      end do
      
c     ... scan for TYPE_M variables

      do i=1,nu*(mi+de)
         lws(i) = 0
      end do
      do i=1,n
         if (vr_reg(2+2*sn+i).eq.TYPE_M) then
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
            lws((j-1)*(nu+1)+k) = i
         end if
         if (vr_reg(2+2*sn+i).eq.TYPE_C) then
            j = vr_reg(2+2*sn+n+i)
            if (j.ge.1)
     &           lws((j-1)*(nu+1)+nu+1) = i
         end if
      end do

      write(nout,*)
      write(nout,*)  'Composition of mixers (in terms of nutrients):'
      write(nout,*)
      do j=1,mi
         write(nout,*) 'Mixer ',j
         do i=1,nu
            ix = lws((j-1)*(nu+1)+i)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'nu-',i, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
      end do
      write(nout,*)
      write(nout,*)  
     &     'Specification of products (in terms of nutrients):'
      write(nout,*)
      do k=1,de
         write(nout,*) 'Product ',k
         do l=1,nu
            ix = lws((k+mi-1)*(nu+1)+l)
            if (ix.ne.0) write(nout,'(A,I3,4F15.8,I4)')
     &           'nu-',l, blo(ix), x(ix), bup(ix), lam(ix), status(ix)
         end do
      end do
      
      write(nout,*)
      write(nout,*)  'Raw material percentage in Products'
      write(nout,*)
c     ... scan for TYPE_NU

      do i=1,n
         if (vr_reg(2+2*sn+i).eq.TYPE_NU) then
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
            if (j.le.rm) then
               if (k.le.mi) then
                  write(nout,'(A,I3,A,I3,4F15.8,I4)')
     &          'rm-',j, ': mi-',k, blo(i),x(i),bup(i),lam(i),status(ix)
               else
                  write(nout,'(A,I3,A,I3,4F15.8,I4)')
     &       'rm-',j, ': de-',k-mi, blo(i),x(i),bup(i),lam(i),status(ix)
               end if
            else
               if (k.le.mi) then
                  write(nout,'(A,I3,A,I3,4F15.8,I4)')
     &       'mi-',j-rm, ': mi-',k, blo(i),x(i),bup(i),lam(i),status(ix)
               else
                  write(nout,'(A,I3,A,I3,4F15.8,I4)')
     &    'mi-',j-rm, ': de-',k-mi, blo(i),x(i),bup(i),lam(i),status(ix)
               end if
            end if
         end if
      end do

      write(nout,*)
      write(nout,*)  'Constraints on raw material availability/'//
     &     'bin capacity'
      write(nout,*)
      do i=1,n
         if (vr_reg(2+2*sn+i).eq.TYPE_U) then
            j = vr_reg(2+2*sn+n+i)
            if (j.le.rm) then
               write(nout,'(A,I3,4F15.8)')
     &              'rm-',j, blo(i),x(i),bup(i),lam(i)
            else
               write(nout,'(A,I3,4F15.8)')
     &              'mi-',j-rm, blo(i),x(i),bup(i),lam(i)
            end if
         end if
      end do

      write(nout,*)
      write(nout,*)  'General linear constraints'
      write(nout,*)
      
      do i=1,m
         pi = iuser(pi_cs_pi+i)
         if (iuser(pi+6).eq.CTYPE_SPECIAL) then
            j = iuser(pi+7)
               write(nout,'(A,I3,4F15.8)')
     &             'no:',j, blo(n+i),c(i),bup(n+i),lam(n+i)
         end if
      end do

      write(nout,*)
      write(nout,*)  'Linear constraints on mixer/product composition'
      write(nout,*)
      
      do i=1,m
         pi = iuser(pi_cs_pi+i)
         if (iuser(pi+6).eq.CTYPE_NURAT) then
            j = iuser(pi+7)
            k = iuser(pi+8)
            if (k.lt.mi) then
               write(nout,'(A,I3,A,I3,4F15.8)')
     &             'mi-',k,' no:',j, blo(n+i),c(i),bup(n+i),lam(n+i)
            else
               write(nout,'(A,I3,A,I3,4F15.8)')
     &             'de-',k-mi,' no:',j, blo(n+i),c(i),bup(n+i),lam(n+i)
            end if
         end if
      end do

      write(nout,*)
      write(nout,*)  'General nonlinear constraints'
      write(nout,*)

      do i=m-NLD%n_cs+1,m
         write(nout,'(A,I3,4F15.8)')
     &        'no:',i+NLD%n_cs-m, blo(n+i),c(i),bup(n+i),lam(n+i)
      end do

      
      write(nout,*)
      write(nout,*)  'Violated constraints are:'
      write(nout,*)
      write(nout,*) ' #    TP  I  J       blo        c(i)         bup'

      do i=1,m
         pi = iuser(pi_cs_pi+i)
         viol = max(c(i)-bup(n+i), blo(n+i)-c(i))
         if (viol.gt.1.d-6) then
            write(nout,'(I4,A,3I3,3F12.5)') i,': ',
     &           iuser(pi+6), iuser(pi+7), iuser(pi+8),
     &           blo(i+n), c(i), bup(i+n)
         end if
      end do


      end

