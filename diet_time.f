      integer function ems_time()
c
c     This is a system dependent function time, that returns
c     time elapsed in seconds since some run-independent event
c
c     This subroutine is for SUN systems
      implicit none
      
      integer time

      ems_time = time()

      end
