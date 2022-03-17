      subroutine pmdl(rtcod, dspace, ix, m, n, ngrad, dobj, drlo, drup,
     &              dclo, dcup, mrow, mcol, dels)
      
      implicit none

      integer rtcod, m, n, ix, ngrad
      double precision dspace(*)
      double precision drlo(m), drup(m), dclo(n), dcup(n)
      integer mrow(ngrad), mcol(ngrad)
      double precision dobj(n), dels(ngrad)

      integer i, j

      do i=1,n
         print *, i, dclo(i), dcup(i), dobj(i)
      end do
      do i=1,m
         print *, i, drlo(i), drup(i)
      end do
      do i=1,ngrad
         print *,i, mrow(i), mcol(i), dels(i)
      end do

      end
      
