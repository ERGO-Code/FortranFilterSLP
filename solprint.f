      subroutine sol_print(n, sn, x, blo, bup, f, vr_reg, k_it, 
     &     lws, nlws, ws, nws)
c     
c     sol_print: prints the solution in sol####.dat and sol####.m
c                in subdirectory <dirn>
c
c     f           optimal value
c     n           #vars
c     x(n)        value of vars
c     blo(n)      bounds on vars
c     bup(n)
c     vr_reg      var registering information
c     lws(nlws)   workspace
c     dirn        directory to write file to
c     k_it        #of solve (this is used for the #### in filename)
c
      use interface_mod

      implicit none

      INCLUDE 'SLPCOMM.INC'
      include 'pusr.inc'

      integer n, k_it, nlws, sn, nws
      double precision f

      double precision x(n), blo(n), bup(n), ws(nws)
      integer lws(nlws), vr_reg(2+2*sn+3*n)

      integer i, j, k, l, demin, demax

      character*4 int2char
      integer trim77

      integer        seed, bigiter
      common /summc/ seed, bigiter

        DIMENSION ITEMP(2)
        DOUBLE PRECISION DTEMP
        INTEGER ITEMP
        EQUIVALENCE (ITEMP,DTEMP)

c     ==================== subroutine body =======================

      if (max((mi+de)*(rm+mi), (nu+1)*(mi+de)).gt.nlws) then
         write(nout, *) 'Not enough integer workspace in solprint'
         stop
      end if
      if (mi+de.gt.nws) then
         write(nout, *) 'Not enough double workspace in solprint'
         stop
      end if


      open(12, file='sol'//int2char(k_it)//'_.dat')
      do i=1,n
         write(12,*) i, x(i)
      end do
      close(12)


      open(12,file='sol'//int2char(k_it)//'.dat')

      open(13,file='sol'//int2char(k_it)//'.m')

      write(12,*)  ' ============= SOLUTION ================='
      write(12,*)
      write(12,*)  ' f = ',f
      write(13,*) 'f = ',f
      write(12,*)
      write(12,*)  'Composition of mixers (in terms of raw mat & bins):'
      write(12,*)
      
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

      do i=1,rm
         do j=1,mi
            if (lws((i-1)*(mi+de)+j).gt.0) then
               ws(j) = x(lws((i-1)*(mi+de)+j))
            else
               ws(j) = 0.d0
            end if
         end do
         write(12,'(A,I3,400F8.3)') 'Rm-',i,(ws(j), j=1,mi) 
      end do
      write(12,*) '---------------------------------------------------'
      do i=1,mi
         do j=1,mi
            if (lws((i+rm-1)*(mi+de)+j).gt.0) then
               ws(j) = x(lws((i+rm-1)*(mi+de)+j))
            else
               ws(j) = 0.d0
            end if
         end do
         write(12,'(A,I3,400F8.3)') 'Mi-',i,(ws(j), j=1,mi)
      end do
      write(12,*)
      write(12,*)  'Composition of demand products (in terms of mi+rm):'
      write(12,*)
      k = 1
 10   continue
      demin = (k-1)*9+1
      demax = min(de,demin+8)
      do i=1,mi
         do j=demin,demax
            if (lws((i+rm-1)*(mi+de)+mi+j).gt.0) then
               ws(j) = x(lws((i+rm-1)*(mi+de)+mi+j))
            else
               ws(j) = 0.d0
            end if
         end do
         write(12,'(A,i3,10F8.3)')
     .             'Mi-',i,(ws(j),j=demin,demax)
      end do
      write(12,*) '---------------------------------------------------'
      do i=1,rm
         do j=demin,demax
            if (lws((i-1)*(mi+de)+mi+j).gt.0) then
               ws(j) = x(lws((i-1)*(mi+de)+mi+j))
            else
               ws(j) = 0.d0
            end if
         end do
         write(12,'(A,i3,10F8.3)')
     .             'RM-',i,(ws(j),j=demin,demax)
      end do
      write(12,*)
      k = k+1
      if (demax.lt.de) goto 10

      write(13,*) 'mu = ['
      do i=1,rm+mi
         do j=1,mi+de
            if (lws((i-1)*(mi+de)+j).gt.0) then
               ws(j) = x(lws((i-1)*(mi+de)+j))
            else
               ws(j) = 0.d0
            end if
         end do
         write(13,'(400F20.10)') (ws(j), j=1,mi+de)
      end do
      write(13,*) '];'
      write(13,*)

c     ... scan for TYPE_M variables

c     ... create a look up table:
c         lws((j-1)*(nu+1)+k     is idx of M(j, k)
c         lws((j-1)*(nu+1)+nu+1  is idx of C(j)
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
c            print *, nu, j, (j-1)*(nu+1)+nu+1
            if (j.ge.1) 
     &           lws((j-1)*(nu+1)+nu+1) = i
         end if
      end do

      write(12,*)
      write(12,*)  'Composition of mixers in terms of nutrients'
      do i=1,nu
         write(12,'(A,I3,400F8.3)') 'Nu-',i,(x(lws((j-1)*(nu+1)+i)), 
     &        j=1,mi)
      end do
      write(12,*)
      write(12,*)  'Cost of mixers'
      write(12,'(A,10F8.3)') '      ',(x(lws((j-1)*(nu+1)+nu+1)), 
     &     j=1,mi)
      write(12,*)
      write(12,*)  'Composition of demand products by nutrients'
      do k=1,de
         write(12,'(A,I3,A,F10.5)') 'demand product ',k,
     .                       ' Cost =',x(lws((k+mi-1)*(nu+1)+nu+1))
         do l=1,nu
            write(12,'(A,I3,3F10.5)') 'Nu-',l,
     .       blo(lws((k+mi-1)*(nu+1)+l)),x(lws((k+mi-1)*(nu+1)+l)),
     .       bup(lws((k+mi-1)*(nu+1)+l))
         end do
         write(12,*)
      end do


      write(13,*) 'm = ['
      do i=1,mi+de
         write(13,'(400F20.10)') (x(lws((i-1)*(nu+1)+j)), j=1,nu)
      end do
      write(13,*) '];'
      write(13,*)
      write(13,*) 'mlb = ['
      do i=1,mi+de
         write(13,'(400F20.10)') (blo(lws((i-1)*(nu+1)+j)), j=1,nu)
      end do
      write(13,*) '];'
      write(13,*)
      write(13,*) 'mub = ['
      do i=1,mi+de
         write(13,'(400F20.10)') (bup(lws((i-1)*(nu+1)+j)), j=1,nu)
      end do
      write(13,*) '];'
      write(13,*)
      write(13,*) 'c = ['
      write(13,'(400F20.10)') (x(lws((j-1)*(nu+1)+nu+1)), j=1,mi+de)
      write(13,*) '];'
      write(13,*)
      

      WRITE(nout,*) ' '
      WRITE(nout,*) ' write solution to file: sol'//int2char(k_it)//'.m'
      WRITE(nout,*) ' '

CJCH      open(12,file=dirname(1:trim77(dirname))//'/sol'//
CJCH     .                                  int2char(k_it)//'.m')

c      do j=1,mi
c         ws(j) = 0.
c         do k=1,de
c            ws(j) = ws(j) + x(pv_lamm+(k-1)*mi+j)*user(rm+nu*rm+k)
c         end do
c      end do
c      write(12,*) 'tonmi= ['
c      do i=1,mi
c         write(12,'(F20.10)') ws(i)
c      end do
c      write(12,*) '];'
c      write(12,*)
c      write(12,*) 'totnumi= ['
c      do i=1,nu
c         write(12,'(20F20.10)') (x(pv_m+(j-1)*nu+i)*ws(j), j=1,mi)
c      end do
c      write(12,*) '];'
c      close(12)

C----------------------- CODE INSERTED BY JCH ------------------------
c
c100     DO 140 J=1,DE
c        WRITE(20,'(I5,A,A)') J,COMMA,'Product'
c           DO 120 I=1,RM
c           DTEMP=x(pv_lamr+(j-1)*rm+i)
c           WRITE(20,'(2Z8,1X,G16.8,A,I5,A,A)')
c     &                                  ITEMP,
c     &                                  x(pv_lamr+(j-1)*rm+i),COMMA,
c     &                                  I,COMMA,'Rawmat'
c120        CONTINUE
c           DO 130 I=1,MI
c           DTEMP=x(pv_lamm+(j-1)*mi+i)
c           WRITE(20,'(2Z8,1X,G16.8,A,I5,A,A)')
c     &                                  ITEMP,
c     &                                  x(pv_lamm+(j-1)*mi+i),COMMA,
c     &                                  I,COMMA,'Premix'
c130        CONTINUE
c140     CONTINUE
c
c        DO 160 J=1,MI
c        WRITE(20,'(I5,A,A)') J,COMMA,'Premix'
c           DO 150 I=1,RM
c           DTEMP=x(pv_mu+(j-1)*rm+i)
c           WRITE(20,'(2Z8,1X,G16.8,A,I5,A,A)')
c     &                                  ITEMP,
c     &                                  x(pv_mu+(j-1)*rm+i),COMMA,
c     &                                  I,COMMA,'Rawmat'
c150        CONTINUE
c160     CONTINUE
c
cC----------------------------------------------------------------------

      end


