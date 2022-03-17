      subroutine filterSLP (n,m,maxa,maxf,mxwk,mxiwk,iprint,
     &                     ifail,rifail,rho,x,c,f,hc, blo,bup,s,ws,
     &                     lws,lam,user,iuser,max_iter,istat,
     &                     rstat, sv_t_fi, se_mx_lp_it, wr_msg,
     &                     status, result, oldbnd, do_rng)

      implicit none

      INCLUDE 'SLPCOMM.INC'

c     ----------- declaration of passed parameters -------------------

      integer n, m, kmax, maxa, maxf, mxwk, mxiwk, iprint,
     &        ifail, rifail, max_iter

      double precision rho, f, hc

      integer lws(mxiwk), iuser(*), istat(12), status(n+m)


c     dChange
c     for new shape phase one 
c
c     change number of columns from n+m to n+2*m in
c     p_sv_bs_c, p_dclo, p_dcup, p_dobj
c
c     change number of rows from 2*m to m in
c     p_sv_bs_r, p_drlo, p_drup
c
c     change total row/col info from  2*ngrad+2*m to ngrad+2*m in
c     p_mrow, p_mcol,  p_dels
c
c     endChange

      double precision x(n), c(m), blo(n+m), bup(n+m), s(n+m),
     &                 ws(mxwk), lam(n+m), user(*), rstat(3),
     &                 result(7,n+m), oldbnd(2*n+2*m)


      logical sv_t_fi, se_mx_lp_it, wr_msg, do_rng

c     ------------ declaration of internal parameters ----------------

      integer p_drlo, p_drup, p_dclo, p_dcup, p_mrow,
     &        p_ws, p_lws, p_mcol, p_dels, p_dobj,
     &        p_d, p_xnew, p_cnew, p_ff, p_fc, p_sv_bs_r, p_sv_bs_c,
     &        p_a, p_la

      integer maxws, maxlws
      integer ngrad, i

      ngrad = maxa - n

c     workspace: SLP needs
c        10*m + 6*n + 3*ngrad + 2*maxf  real workspace, rest is passed
c                                       as dspace/ispace to EMSOL
c        8*m + 2*n + 5*ngrad + 3        integer workspace, rest is
c                                       passed as lws to slpmain
c                                       (and never used)


c     storage map for integer data

      p_la   = 1
      p_mrow = p_la + maxa+m+3
      p_mcol = p_mrow + ngrad+2*m
      p_sv_bs_r = p_mcol + ngrad + 2*m
      p_sv_bs_c = p_sv_bs_r + m
      p_lws     = p_sv_bs_c + n+2*m

c     storage map for double precision data

      p_a    = 1
      p_drlo = p_a    + maxa
      p_drup = p_drlo + m
      p_dclo = p_drup + m
      p_dcup = p_dclo + n+2*m
      p_dels = p_dcup + n+2*m
      p_dobj = p_dels + ngrad + 2*m
      p_xnew = p_dobj + n+2*m
      p_cnew = p_xnew + n
      p_d    = p_cnew + m
      p_ff   = p_d    + n
      p_fc   = p_ff   + maxf
      p_ws   = p_fc   + maxf

      maxws = mxwk - p_ws +1
      maxlws = mxiwk - p_lws +1

      if (maxws.le.0.or.maxlws.le.0) then
         WRITE(nout,*) 'FATAL ERROR: not enough workspace in FilterSLP'
         WRITE(nout,*) 'maxws = ',maxws,' maxlws = ',maxlws
         ifail = 6
         goto 99
      end if


c      do i=1,m
c         print *, i, blo(n+i), bup(n+i)
c      end do

      
      call slp_main(n, m, maxa, ngrad, maxws, maxlws, iprint,
     &              ifail, rifail, maxf, rho, f, hc, x, ws(p_xnew), c,
     &              ws(p_cnew), ws(p_drlo), ws(p_drup), ws(p_dclo),
     &              ws(p_dcup), ws(p_a), lws(p_la),
     &              lws(p_mrow), lws(p_mcol), ws(p_dels), ws(p_dobj),
     &              ws(p_d), ws(p_ff), ws(p_fc),
     &              lws(p_lws), ws(p_ws), user, iuser, max_iter,
     &              istat, rstat, blo, bup, lam, ws(p_ws),
     &              lws(p_sv_bs_r), lws(p_sv_bs_c), s,
     &              sv_t_fi, se_mx_lp_it, wr_msg, status, result,
     &              oldbnd, do_rng)

 99   continue

      end
