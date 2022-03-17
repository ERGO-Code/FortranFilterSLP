      subroutine ems_itru(dspace,ispace,reason,userrtcd)
      implicit none
c     read in only
      double precision dspace(*)
      integer ispace(*)
      integer reason
c     read in maybe change
      integer userrtcd
c     extra for calling iget and rget
      integer rtcod
c     extra for judging feasibility
      double precision stp_eps, csr_eps
      common /tol/     stp_eps, csr_eps
     
      include 'EMSOLI.INC'
      include 'EMSOLR.INC'

c     this called every time emsol resets
c     checks to see if we are approx feasible and if so aborts
c     by telling emsol max LP iterations have been reached
c     saves iterations on hard problems where otherwise emsol keeps getting
c     close to feas point, refactoring finding its not feas and carrying on

        if (1==1) then

           call ems_rget(rtcod, dspace, emsolr, emsolrln)
           call ems_iget(rtcod, dspace, emsoli, emsoliln)
           
           if ( (rsumpinf<=50*rtolpinf) .and. (rsumdinf<=50*rtoldinf) 
     &          .and. (inumpinf>0 .or. inumdinf>0)
     &          .and. (Iiternum>1000) ) then
              userrtcd = 3      !make emsol exit as if max its
              print '(a,i5,a,i5,a,g15.7,a,g15.7,a,2g15.7)',  
     &             'inumpinf', inumpinf, ' inumdinf',  inumdinf,
     &             ' rsumpinf', rsumpinf, ' rsumdinf ' ,rsumdinf 
     &             , ' pd tols', rtolpinf, rtoldinf
           end if   
        end if

      end subroutine ems_itru
