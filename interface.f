      module interface_mod

c     This provides the interface (and routines to access it)
c     to communicate information on the problem size and the
c     splicing of the primal/dual vector to the caller
c     
c     It basically stores all the information that is read in (and set
c     up) by the call to rd_prob_dim. 


      implicit none

c     ... save so that all information is stored between accesses to
c     the module
      save 

c      type interface_type

      integer nu, rm, mi, de
      integer ntgdim

c      integer nnurat, nrminc, nspec, nnl

c      end type interface_type


      end module
