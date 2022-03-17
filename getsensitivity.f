        SUBROUTINE GETSENSITIVITY(n, sn, m, x, vr_reg, iuser, 
     &                            blo, bup, lam, c, nlws, lws,
     &                            result, status, ifail,
     &                            NTGRESULT, NRESULTS)
*
*       +---------------------------------------------------+
*       | Copyright (c) FORMAT Int.'03. All rights reserved |
*       +---------------------------------------------------+
*       | Subroutine :  GETSENSITIVITY                      |
*       | Description:  Get the sensitivity data            |
*       | Revision   :  September 2003                      |
*       +---------------------------------------------------+
*
*------------------------
*  Use Statements       *
*------------------------
*

      use interface_mod

*------------------------
*  Implicit definition  *
*------------------------
*
        IMPLICIT NONE
*
*-------------------
*  Include common  *
*-------------------
*
        INCLUDE 'SLPCOMM.INC'
        INCLUDE 'pusr.inc'
        INCLUDE 'ntgres.inc'
*
*------------------
*  Indirect data  *
*------------------
*
        INTEGER          n, sn, m, nlws, ifail, NRESULTS
*        
        INTEGER          iuser   (*), 
     &                   vr_reg  (2+2*sn+3*n),
     &                   lws     (nlws),
     &                   status  (n+m)
*
        DOUBLE PRECISION x        (n), 
     &                   blo      (n+m),
     &                   bup      (n+m), 
     &                   lam      (n+m),
     &                   c        (m),
     &                   result   (7,n+m),
     &                   NTGRESULT(*)
*
*---------------
*  Local data  *
*---------------
*
        DIMENSION        VIOLTYPE(m),
     &                   VIOLPROD(m),
     &                   VIOLNTRM(m),
     &                   VIOLVALUE(m),
     &                   REMSTART(m),
     &                   REMPRODUCT(m),
     &                   REMNUTRIENT(m),
     &                   REMVIOLATED(m)
*
        INTEGER          i, ix, j, k, l, pi, tp, p, NDX,
     &                   NUMVIOLATED,VIOLTYPE,VIOLPROD,VIOLNTRM,
     &                   CHANGETYPE,SLPLINE,RATIORESULT,NUTNDX1,
     &                   NUTNDX2,TEMPSTATUS,PROD,INGR,GIGRESULT,
     &                   NUMREM,REMSTART,REMPRODUCT,REMNUTRIENT,
     &                   REMVIOLATED,REMCOUNT
*     
        INTEGER          RESULTINDEX,RATIOINDEX,GIGINDEX,RELINDEX,
     &                   INFINDEX
*        
        DOUBLE PRECISION viol,VIOLVALUE,LEVEL,DLOWER,DUPPER,
     &                   LEVEL1,LEVEL2,DIFFLEVEL,SENSITIVITY,
     &                   DMIXINCL

*
*---------------
*  Initialise  *
*---------------
*
*
*---------------------------------
*  Not enough integer workspace  *
*---------------------------------
*
        IF(max((mi+de)*(rm+mi), (nu+1)*(mi+de)).gt.nlws) then
          PRINT *,'Not enough integer workspace in lamprint'
          STOP
        END IF
*
*-------------
*  Get mu's  *
*-------------
*
        DO i=1,(mi+de)*(rm+mi)
          lws(i) = 0
        END DO
*      
        DO i=1,n
          IF (vr_reg(2+2*sn+i).eq.TYPE_MU) THEN
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
            lws((j-1)*(mi+de)+k) = i
          END IF
        END DO
*
*----------------------------
*  Fill the GIG index once  *
*----------------------------
*
        IF(GLO_INDEX_FLAG.EQ.0) THEN
          DO i=1,NUM_GIG_RESULTS
            NDX =GIGINDEX(i,6,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
            GLOINDEX(i,1)=NTGRESULT(NDX)
            NDX =GIGINDEX(i,7,nu,rm,mi,de,NUM_RATIO_RESULTS,NRESULTS)
            GLOINDEX(i,2)=NTGRESULT(NDX)
          END DO
          GLO_INDEX_FLAG=1
        END IF
*
*----------------------------
*  Fill the REL index once  *
*----------------------------
*
        IF(REL_INDEX_FLAG.EQ.0) THEN
          l=1
          DO i=1,m
            pi = iuser(pi_cs_pi+i)
            IF (iuser(pi+6).eq.CTYPE_SPECIAL) THEN
              j = iuser(pi+7)
*
*  Get number of constraints
*
              NDX=RELINDEX(l,1,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                     NUM_GIG_RESULTS,NRESULTS)
              p=NTGRESULT(NDX)
*            
*  Keep the same NDX if have 2 constraints per relationship
*
              RINDEX(j)=l
*            
              IF(p.EQ.1) THEN
                l=l+1
              ELSE
                NTGRESULT(NDX)=1D0
              END IF
            END IF
          END DO
          REL_INDEX_FLAG   =1
          NUM_RELATIONSHIPS=l-1
        END IF
*
*-----------------------------------
*  Fill the NUTR ratio index once  *
*-----------------------------------
*
        IF(RATIO_NUTR_INDEX_FLAG.EQ.0) THEN
          DO i=1,NUM_RATIO_RESULTS
            NDX=RATIOINDEX(i,4,nu,rm,mi,de,NRESULTS)
            RATIO_NUTR_INDEX(i,1)=NTGRESULT(NDX)
            NDX=RATIOINDEX(i,5,nu,rm,mi,de,NRESULTS)
            RATIO_NUTR_INDEX(i,2)=NTGRESULT(NDX)
          END DO
          RATIO_NUTR_INDEX_FLAG=1
        END IF
*
*-------------------------
*  Violated constraints  *
*-------------------------
*
*  #    TP  I  J       blo        c(i)         bup
        NUMVIOLATED=0
        DO i=1,m
          pi   = iuser(pi_cs_pi+i)
          viol = max(c(i)-bup(n+i), blo(n+i)-c(i))
          IF (viol.gt.1.d-6) THEN
            NUMVIOLATED=NUMVIOLATED+1
            VIOLTYPE(NUMVIOLATED) =iuser(pi+6)
            VIOLPROD(NUMVIOLATED) =iuser(pi+7)
            VIOLNTRM(NUMVIOLATED) =iuser(pi+8)
            VIOLVALUE(NUMVIOLATED)=c(i)
          END IF
        END DO
*
*------------------------------------------------------
*  Composition of mixers (in terms of raw materials)  *
*------------------------------------------------------
*
        DO j=1,mi
          DO i=1,rm
            ix = lws((i-1)*(mi+de)+j)
            IF (ix.ne.0) THEN
*
*  Inclusion
* 
              NDX=RESULTINDEX(j,i+nu,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
*              
              NDX=RESULTINDEX(j,i+nu,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'rm-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
          DO i=1,mi
            ix = lws((rm+i-1)*(mi+de)+j)
            IF (ix.ne.0) THEN
*
*  Inclusion
* 
              NDX=RESULTINDEX(j,i+nu+rm,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
* 
              NDX=RESULTINDEX(j,i+nu+rm,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'mi-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
        END DO
*
*------------------------------------------------
*  Composition of products (in terms of rm/mi)  *
*------------------------------------------------
*
        DO j=1,de
          DO i=1,rm
            ix = lws((i-1)*(mi+de)+mi+j)
            IF (ix.ne.0) THEN
*
*  Inclusion
* 
              NDX=RESULTINDEX(j+mi,i+nu,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
* 
              NDX=RESULTINDEX(j+mi,i+nu,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'rm-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
          DO i=1,mi
            ix = lws((i+rm-1)*(mi+de)+mi+j)
            IF (ix.ne.0) THEN
*
*  Inclusion
* 
              NDX=RESULTINDEX(j+mi,i+nu+rm,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
* 
              NDX=RESULTINDEX(j+mi,i+nu+mi,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'mi-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
        END DO
*      
*-------------------------------------
*     ... scan for TYPE_M variables  *
*-------------------------------------
*
        DO i=1,nu*(mi+de)
          lws(i) = 0
        END DO
*        
        DO i=1,n
          IF (vr_reg(2+2*sn+i).eq.TYPE_M) THEN
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
            lws((j-1)*(nu+1)+k) = i
          END IF
          IF (vr_reg(2+2*sn+i).eq.TYPE_C) THEN
            j = vr_reg(2+2*sn+n+i)
            if (j.gt.0) then
c       ... make sure this is an TYPE_C: i, 0 entry and not T_C: -1, 0
               lws((j-1)*(nu+1)+nu+1) = i
            end if
         END IF
        END DO
*
*--------------------------------------------------
*  Composition of mixers (in terms of nutrients)  *
*--------------------------------------------------
*
        DO j=1,mi
*  Mixer ',j
          DO i=1,nu
            ix = lws((j-1)*(nu+1)+i)
            IF (ix.ne.0) THEN
*
*  Level
* 
              NDX=RESULTINDEX(j,i,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
* 
              NDX=RESULTINDEX(j,i,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'nu-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
        END DO
*
*------------------------------------------------------
*  Specification of products (in terms of nutrients)  *
*------------------------------------------------------
*
        DO j=1,de
*  Product ',j
          DO i=1,nu
            ix = lws((j+mi-1)*(nu+1)+i)
            IF (ix.ne.0) THEN
*
*  Level
* 
              NDX=RESULTINDEX(j+mi,i,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(ix)
*
*  Sensitivity
* 
              NDX=RESULTINDEX(j+mi,i,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(ix)
*              write(*,'(A,I3,4F15.8)')
*     &           'nu-',i, blo(ix), x(ix), bup(ix), lam(ix)
            END IF
          END DO
        END DO
*
*----------------------------------------
*  Raw material percentage in Products  *
*     ... scan for TYPE_NU              *
*----------------------------------------
*
*  j: Ingredient
*  k: Produkt (mi/de)
*
        DO i=1,n
          IF (vr_reg(2+2*sn+i).eq.TYPE_NU) THEN
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
*
            DO l=1,NUM_GIG_RESULTS
              PROD=GLOINDEX(l,1)
              INGR=GLOINDEX(l,2)
*              
              GIGRESULT=l
*     
              IF(PROD.EQ.k.AND.INGR.EQ.J) THEN
*
*  Set STATUS and level
*
                CALL GETGIGSTATUS(ifail,
     &                            lam(i),
     &                            NTGRESULT,
     &                            GIGRESULT,
     &                            NRESULTS,
     &                            x(i),
     &                            PROD,
     &                            INGR,
     &                            NUMVIOLATED,
     &                            VIOLTYPE,
     &                            VIOLPROD,
     &                            VIOLNTRM,
     &                            VIOLVALUE)
*                write(*,'(A,I3,A,I3,4F15.8)')
*     &              'rm-',j, ': mi-',k, blo(i),x(i),bup(i),lam(i)
              END IF
            END DO
          END IF
        END DO
*
*----------------------------------------------------------
*  Constraints on raw material availability/bin capacity  *
*----------------------------------------------------------
*
        DO i=1,n
          IF (vr_reg(2+2*sn+i).eq.TYPE_U) THEN
            j = vr_reg(2+2*sn+n+i)
            IF (j.le.rm) THEN
*
*  Total raw material consumption
* 
              NDX=RESULTINDEX(mi+de+1,j,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(i)
*
*  Sensibility
* 
              NDX=RESULTINDEX(mi+de+1,j,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(i)
*              write(*,'(A,I3,4F15.8)')
*     &              'rm-',j, blo(i),x(i),bup(i),lam(i)
            ELSE
*
*  Total premix consumption
* 
              NDX=RESULTINDEX(mi+de+1,j,1,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=x(i)
              CALL CORRECTSENSITIVITY(j-rm,mi,rm,nu,NRESULTS,
     &                                NTGRESULT,x(i))
*
*  Sensibility
* 
              NDX=RESULTINDEX(mi+de+1,j,4,nu+rm+mi,NRESULTS)
              NTGRESULT(NDX)=lam(i)
*              write(*,'(A,I3,4F15.8)')
*     &              'mi-',j-rm, blo(i),x(i),bup(i),lam(i)
            END IF
          END IF
        END DO
*
*-----------------
*  Reset STATUS  *
*-----------------
*
        DO l=1,NUM_RATIO_RESULTS
          NDX=RATIOINDEX(l,8,nu,rm,mi,de,NRESULTS)
          NTGRESULT(NDX)=FLOATING_BLANK
        END DO
*
*--------------------------
*  Reset REL sensitivity  *
*--------------------------
*
        DO i=1,m
          pi = iuser(pi_cs_pi+i)
          IF (iuser(pi+6).eq.CTYPE_SPECIAL) THEN
            j = iuser(pi+7)
*
            l=RINDEX(j)
*
*  Reset the sensitivity to only get the last set
*
            NDX=RELINDEX(l,4,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                   NUM_GIG_RESULTS,NRESULTS)
            NTGRESULT(NDX)=FLOATING_BLANK
          END IF
        END DO
*
*-----------------------------------------------
*  General linear constraints (Relationships)  *
*-----------------------------------------------
*
        DO i=1,m
          pi = iuser(pi_cs_pi+i)
          IF (iuser(pi+6).eq.CTYPE_SPECIAL) THEN
            j = iuser(pi+7)
*
*  Get the sensitivity
*
            l=RINDEX(j)
            NDX=RELINDEX(l,4,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                   NUM_GIG_RESULTS,NRESULTS)
            IF(NTGRESULT(NDX).EQ.FLOATING_BLANK) THEN
*  First time round
              NTGRESULT(NDX)=lam(n+i)
            ELSE IF(NTGRESULT(NDX).EQ.0D0) THEN
*  Second time round
              NTGRESULT(NDX)=lam(n+i)
            END IF
          END IF
        END DO
*
*----------------------------------------------------
*  Linear constraints on mixer/product composition  *
*----------------------------------------------------
*
        CHANGETYPE=0
        DO i=1,m
          pi = iuser(pi_cs_pi+i)
          IF (iuser(pi+6).eq.CTYPE_NURAT) THEN
            j = iuser(pi+7)
            k = iuser(pi+8)
*            if (k.lt.mi) then
*               write(*,'(A,I3,A,I3,4F15.8)')
*     &             'mi-',k,' no:',j, blo(n+i),c(i),bup(n+i),lam(n+i)
*            else
*               write(*,'(A,I3,A,I3,4F15.8)')
*     &             'de-',k-mi,' no:',j, blo(n+i),c(i),bup(n+i),lam(n+i)
*            end if
*
*  k: de/mi index
*  j: sequence number for ratio constraint
*
            IF(CHANGETYPE.EQ.0) THEN
*
*  Products
*
              DO l=1,NUM_RATIO_RESULTS
                NDX        =RATIOINDEX(l,6,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j) THEN
                  GOTO 100
                END IF
*                  
                NDX        =RATIOINDEX(l,7,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j) THEN
                  GOTO 100
                END IF
              END DO
              write(*,*)'FATAL ERROR !!. '//
     &                  'Product RESULT index not found for ratio!'
              STOP
*
*  Set STATUS for a product
*
100           CALL GETRATIOSTATUS(ifail,
     &                            lam(n+i),
     &                            NTGRESULT,
     &                            RATIORESULT,
     &                            NRESULTS)
              IF(j.EQ.n_nuratp) THEN
                CHANGETYPE=1
              END IF
            ELSE
*
*  Premixes
*
              DO l=1,NUM_RATIO_RESULTS
                NDX        =RATIOINDEX(l,6,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j+n_nuratp) THEN
                  GOTO 110
                END IF
*
                NDX        =RATIOINDEX(l,7,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j+n_nuratp) THEN
                  GOTO 110
                END IF
              END DO
              write(*,*)'FATAL ERROR !!. '//
     &                  'Premix RESULT index not found for ratio!'
              STOP
*
*  Set STATUS for a premix
*
110           CALL GETRATIOSTATUS(ifail,
     &                            lam(n+i),
     &                            NTGRESULT,
     &                            RATIORESULT,
     &                            NRESULTS)
            END IF
          END IF
        END DO
*
*----------------------
*  RANGING VARIABLES  *
*----------------------
*
*  Information for variables:
*  I TP J K   RANGING INFORMATION 
*  TP : Type
*  j  : Spec.     if type = 1 (TYPE_M).  Mixes first then products
*  j  : Component if type = 2 (TYPE_MU). Ingredients first then products
*  j  : Component if type = 4 (TYPE_U).  Ingredients first then products
*  k  : Nutrient  if type = 1 (TYPE_M).  Nutrient number
*  k  : Spec.     if type = 2 (TYPE_MU). Mixes first then products
*  k  : NA        if type = 4 (TYPE_U).  Not in use
*
        NUMREM=0
        DO i=1,n
          tp = vr_reg(2+2*sn+i)
          j = vr_reg(2+2*sn+n+i)
          k = vr_reg(2+2*sn+2*n+i)
*
*---------------------------------------------
*  Nutrient. Get Ranging 1, 2, 3 and STATUS  *
*---------------------------------------------
*
          IF(tp.EQ.TYPE_M) THEN
*
*  Status
*
            NDX=RESULTINDEX(j,k,8,nu+rm+mi,NRESULTS)
            NTGRESULT(NDX)=status(i)
            NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
            LEVEL=NTGRESULT(NDX)
*  Get minimum constraint
            NDX=RESULTINDEX(j,k,2,nu+rm+mi,NRESULTS)
            DLOWER=NTGRESULT(NDX)
*  Get maximum constraint
            NDX=RESULTINDEX(j,k,3,nu+rm+mi,NRESULTS)
            DUPPER=NTGRESULT(NDX)
*      write(*,*)'GetSensitivity: ifail,Level,Min,Max,Status : ',ifail,
*     & LEVEL,DLOWER,DUPPER,status(i)
            IF(ifail.EQ.0) THEN
*
*  Get mixer inclusion
*
              CALL GETMIXERINCLUSION(j,n,sn,rm,x,vr_reg,DMIXINCL)
*
*  Sensitivity when unconstrained
*
*              IF(status(i).EQ.1) THEN
                NDX=RESULTINDEX(j,k,4,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=result(4,i)
*
*  Correct the sensitivity
*
                IF(DMIXINCL.NE.0D0) THEN
                  NTGRESULT(NDX)=NTGRESULT(NDX)/DMIXINCL
                END IF
*              END IF
*
              DO l=5,7
                NDX=RESULTINDEX(j,k,l,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=result(l,i)
*
*  Correct the ranging
*
                IF(l.EQ.6.AND.DMIXINCL.NE.0D0) THEN
                  NTGRESULT(NDX)=NTGRESULT(NDX)/DMIXINCL
                END IF
              END DO
            END IF
*
*  Correct the status if not feasible
*
            IF(ifail.NE.0) THEN
              DO l=1,NUMVIOLATED
                IF(VIOLTYPE(l).EQ.TYPE_M.AND.
     &             VIOLPROD(l).EQ.j.AND.
     &             VIOLNTRM(l).EQ.k) THEN
*  Get level
                  NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
                  LEVEL=NTGRESULT(NDX)
*  Get minimum constraint
                  NDX=RESULTINDEX(j,k,2,nu+rm+mi,NRESULTS)
                  DLOWER=NTGRESULT(NDX)
*  Get maximum constraint
                  NDX=RESULTINDEX(j,k,3,nu+rm+mi,NRESULTS)
                  DUPPER=NTGRESULT(NDX)
*                  
                  IF(LEVEL+VIOLVALUE(l).LT.DLOWER) THEN
                    NDX=RESULTINDEX(j,k,8,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=4
*  Set level = level+violvalue
                    NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=LEVEL+VIOLVALUE(l)
                  ELSE IF(LEVEL+VIOLVALUE(l).GT.DUPPER) THEN
                    NDX=RESULTINDEX(j,k,8,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=5
*  Set level = level+violvalue
                    NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=LEVEL+VIOLVALUE(l)
                  ELSE
                    NUMREM=NUMREM+1
                    REMSTART(NUMREM)   =j
                    REMPRODUCT(NUMREM) =j
                    REMNUTRIENT(NUMREM)=k
                    REMVIOLATED(NUMREM)=l
                  END IF
                END IF
              END DO
            END IF
          END IF
*
*-----------------------------------------------
*  Ingredient. Get Ranging 1, 2, 3 and STATUS  *
*-----------------------------------------------
*
          IF(tp.EQ.TYPE_MU) THEN
*
*  Status
* 
            NDX=RESULTINDEX(k,j+nu,8,nu+rm+mi,NRESULTS)
            NTGRESULT(NDX)=status(i)

            NDX=RESULTINDEX(k,j+nu,1,nu+rm+mi,NRESULTS)
            LEVEL=NTGRESULT(NDX)

            NDX=RESULTINDEX(k,j+nu,2,nu+rm+mi,NRESULTS)
            DLOWER=NTGRESULT(NDX)

            NDX=RESULTINDEX(k,j+nu,3,nu+rm+mi,NRESULTS)
            DUPPER=NTGRESULT(NDX)
*
            IF(ifail.EQ.0) THEN
*
*  Get mixer inclusion
*
              CALL GETMIXERINCLUSION(k,n,sn,rm,x,vr_reg,DMIXINCL)
*
*  Sensitivity when unconstrained
*
              IF(status(i).EQ.1) THEN
                NDX=RESULTINDEX(k,j+nu,4,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=result(4,i)
*
*  Correct the sensitivity
*
                IF(DMIXINCL.NE.0D0) THEN
                  NTGRESULT(NDX)=NTGRESULT(NDX)/DMIXINCL
                END IF
              END IF
*
              DO l=5,7
                NDX=RESULTINDEX(k,j+nu,l,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=result(l,i)
*
*  Correct the ranging
*
                IF(l.EQ.6.AND.DMIXINCL.NE.0D0) THEN
                  NTGRESULT(NDX)=NTGRESULT(NDX)/DMIXINCL
                END IF
              END DO
            END IF
          END IF
*
*-----------------------------------------------------------------
*  Total ingredient consumption. Get Ranging 1, 2, 3 and STATUS  *
*-----------------------------------------------------------------
*
          IF(tp.EQ.TYPE_U) THEN
            IF(ifail.EQ.0) THEN
              DO l=5,7
                NDX=RESULTINDEX(mi+de+1,j,l,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=result(l,i)
              END DO
            END IF
            NDX=RESULTINDEX(mi+de+1,j,8,nu+rm+mi,NRESULTS)
            NTGRESULT(NDX)=status(i)
*
*  Correct the status if not feasible
*
            IF(ifail.NE.0) THEN
              DO l=1,NUMVIOLATED
                IF(VIOLTYPE(l).EQ.TYPE_U.AND.
     &             VIOLPROD(l).EQ.j) THEN
*  Get minimum constraint
                  NDX=RESULTINDEX(mi+de+1,j,2,nu+rm+mi,NRESULTS)
                  DLOWER=NTGRESULT(NDX)
*  Get maximum constraint
                  NDX=RESULTINDEX(mi+de+1,j,3,nu+rm+mi,NRESULTS)
                  DUPPER=NTGRESULT(NDX)
*  Get level
                  NDX=RESULTINDEX(mi+de+1,j,1,nu+rm+mi,NRESULTS)
                  LEVEL=NTGRESULT(NDX)
*  Add violated value to level
                  LEVEL=LEVEL+VIOLVALUE(l)
*  Store level
                  NTGRESULT(NDX)=LEVEL
*                  
                  IF(LEVEL.LT.DLOWER) THEN
                    NDX=RESULTINDEX(mi+de+1,j,8,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=4
                  ELSE IF(LEVEL.GT.DUPPER) THEN
                    NDX=RESULTINDEX(mi+de+1,j,8,nu+rm+mi,NRESULTS)
                    NTGRESULT(NDX)=5
                  END IF
                END IF
              END DO
            END IF
          END IF
*
*------------------------------
*  GIGs. Get Ranging 1, 2, 3  *
*------------------------------
*
          IF(tp.EQ.TYPE_NU) THEN
            IF(ifail.EQ.0) THEN
*              
              DO l=1,NUM_GIG_RESULTS
                PROD=GLOINDEX(l,1)
                INGR=GLOINDEX(l,2)
*              
                IF(PROD.EQ.k.AND.INGR.EQ.J) THEN
*
*  GIG found
*
                  DO p=5,7
                    NDX=GIGINDEX(l,p,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                           NRESULTS)
                    NTGRESULT(NDX)=result(p,i)
                  END DO
                END IF
              END DO
            END IF
          END IF
*          
*          write(nout,'(I3,3I2,1X,7(F14.6,1X))')i,tp,j,k,
*     &         (result(l,i),l=1,7)
        END DO
*
*------------------------
*  RANGING CONSTRAINTS  *
*------------------------
*
*  Information for constraints'
*  I TP J K   RANGING INFORMATION '
        CHANGETYPE=0
        DO i=1,m
          pi = iuser(pi_cs_pi+i)
          tp = iuser(pi+6)
          j  = iuser(pi+7)
          k  = iuser(pi+8)
*
*  RELs. Get Ranging 1, 2, 3
*
          IF(tp.EQ.CTYPE_SPECIAL.AND.ifail.EQ.0) THEN
*
*  Get the sensitivity
*
            l=RINDEX(j)
            NDX=RELINDEX(l,4,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                   NUM_GIG_RESULTS,NRESULTS)
            LEVEL=NTGRESULT(NDX)
*
*  Get ranging info for constrained term
*
            IF(LEVEL.NE.0D0) THEN
              DO p=5,7
                NDX=RELINDEX(l,p,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                       NUM_GIG_RESULTS,NRESULTS)
                NTGRESULT(NDX)=result(p,n+i)
              END DO
            END IF
          END IF
*
*  Nutrient ratios. Get Ranging 1, 2, 3 and STATUS
*
          IF(tp.EQ.CTYPE_NURAT.AND.ifail.EQ.0) THEN
*
*  k: de/mi index
*  j: sequence number for ratio constraint
*
            IF(CHANGETYPE.EQ.0) THEN
*
*  Products
*
              DO l=1,NUM_RATIO_RESULTS
                NDX        =RATIOINDEX(l,6,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j) GOTO 130
*                  
                NDX        =RATIOINDEX(l,7,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j) GOTO 130
              END DO
              write(*,*)'FATAL ERROR !!. '//
     &                  'Product RESULT index not found for ratio!'
              STOP
*
*  Set nutrient ratio sensitivity for a product
*
130           NDX=RATIOINDEX(l,8,nu,rm,mi,de,NRESULTS)
              IF(NTGRESULT(NDX).EQ.1D0) THEN
                NDX=RATIOINDEX(l,4,nu,rm,mi,de,NRESULTS)
                NTGRESULT(NDX)=result(4,n+i)
              END IF
*
              IF(j.EQ.n_nuratp) THEN
                CHANGETYPE=1
              END IF
            ELSE
*
*  Premixes
*
              DO l=1,NUM_RATIO_RESULTS
                NDX        =RATIOINDEX(l,6,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j+n_nuratp) GOTO 135
*
                NDX        =RATIOINDEX(l,7,nu,rm,mi,de,NRESULTS)
                SLPLINE    =NTGRESULT(NDX)
                RATIORESULT=l
                IF(SLPLINE.EQ.j+n_nuratp) GOTO 135
              END DO
              write(*,*)'FATAL ERROR !!. '//
     &                  'Premix RESULT index not found for ratio!'
              STOP
*
*  Set nutrient ratio sensitivity for a premix
*
135           NDX=RATIOINDEX(l,8,nu,rm,mi,de,NRESULTS)
              IF(NTGRESULT(NDX).EQ.1D0) THEN
                NDX=RATIOINDEX(l,4,nu,rm,mi,de,NRESULTS)
                NTGRESULT(NDX)=result(4,n+i)
              END IF
*         write(*,'(I3,3I2,7(F14.6,1X))') i, tp, j, k, 
*     &        (result(l,n+i), l=1,7)
            END IF
          END IF
        END DO
*
*--------------------------------------------------------------------
*  Find the correct product for the infeasible nutrient constraint  *
*--------------------------------------------------------------------
*
        DO 140 P=1,NUMREM
          REMCOUNT=0
*
*  Try the next premix/product
*
120       REMSTART(P)=REMSTART(P)+1
*          
          IF(REMSTART(P).GT.mi+de) REMSTART(P)=1
*
*  Stop when comming round
*            
          IF(REMSTART(P).EQ.REMPRODUCT(P)) THEN
            NDX=RESULTINDEX(REMPRODUCT(P),REMNUTRIENT(P),8,nu+rm+mi,
     &                      NRESULTS)
            NTGRESULT(NDX)=1
            GOTO 140
          END IF
*
*  Do not correct forever
*
          REMCOUNT=REMCOUNT+1
          IF(REMCOUNT.GT.mi+de) THEN
            NDX=RESULTINDEX(REMPRODUCT(P),REMNUTRIENT(P),8,nu+rm+mi,
     &                      NRESULTS)
            NTGRESULT(NDX)=1
            GOTO 140
          END IF
*          
          DO i=1,n
            tp = vr_reg(2+2*sn+i)
            j = vr_reg(2+2*sn+n+i)
            k = vr_reg(2+2*sn+2*n+i)
*
*  Nutrient.
*
            IF(tp.EQ.TYPE_M.AND.
     &         j.EQ.REMSTART(P).AND.
     &         k.EQ.REMNUTRIENT(P)) THEN
*  Get level
              NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
              LEVEL=NTGRESULT(NDX)
*  Get minimum constraint
              NDX=RESULTINDEX(j,k,2,nu+rm+mi,NRESULTS)
              DLOWER=NTGRESULT(NDX)
*  Get maximum constraint
              NDX=RESULTINDEX(j,k,3,nu+rm+mi,NRESULTS)
              DUPPER=NTGRESULT(NDX)
*  Get sensitivity
              NDX=RESULTINDEX(j,k,4,nu+rm+mi,NRESULTS)
              SENSITIVITY=NTGRESULT(NDX)
*              
              IF(LEVEL+VIOLVALUE(REMVIOLATED(P)).LT.DLOWER.AND.
     &           SENSITIVITY.GT.0D0) THEN
                NDX=RESULTINDEX(j,k,8,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=4
*  Set level = level+violvalue
                NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=LEVEL+VIOLVALUE(REMVIOLATED(P))
              ELSE IF(LEVEL+VIOLVALUE(REMVIOLATED(P)).GT.DUPPER.AND.
     &                SENSITIVITY.LT.0D0) THEN
                NDX=RESULTINDEX(j,k,8,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=5
*  Set level = level+violvalue
                NDX=RESULTINDEX(j,k,1,nu+rm+mi,NRESULTS)
                NTGRESULT(NDX)=LEVEL+VIOLVALUE(REMVIOLATED(P))
              ELSE
                GOTO 120
              END IF
            END IF
          END DO
140     CONTINUE
*
*------------------------------------
*  Store the infeasibles in RESULT  *
*------------------------------------
*
        DO I=1,NUMVIOLATED
          IF(I.LE.MAX_INFEASIBLES) THEN
            IF(VIOLTYPE(I).EQ.TYPE_M.OR.
     &         VIOLTYPE(I).EQ.TYPE_MU.OR.
     &         VIOLTYPE(I).EQ.TYPE_U.OR.
     &         VIOLTYPE(I).EQ.TYPE_NU) THEN
*  TP in index n,1              
              NDX=INFINDEX(I,1,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                     NUM_GIG_RESULTS,NUM_RELATIONSHIPS,NRESULTS)
              NTGRESULT(NDX)=VIOLTYPE(I)
*  mi/de in index n,2
              NDX=INFINDEX(I,2,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                     NUM_GIG_RESULTS,NUM_RELATIONSHIPS,NRESULTS)
              NTGRESULT(NDX)=VIOLPROD(I)
*  rm/nut in index n,3
              NDX=INFINDEX(I,3,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                     NUM_GIG_RESULTS,NUM_RELATIONSHIPS,NRESULTS)
              NTGRESULT(NDX)=VIOLNTRM(I)
*  Diff in index n,4
              NDX=INFINDEX(I,4,nu,rm,mi,de,NUM_RATIO_RESULTS,
     &                     NUM_GIG_RESULTS,NUM_RELATIONSHIPS,NRESULTS)
              NTGRESULT(NDX)=VIOLVALUE(I)
            END IF
          END IF
        END DO
*
        END

