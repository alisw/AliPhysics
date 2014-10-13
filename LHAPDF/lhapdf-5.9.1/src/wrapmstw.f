! -*- F90 -*-

!-- LHAPDF version of MSTW interpolation code and alphaS routines.
!-- 22/01/2009 by Graeme Watt <watt(at)hep.ucl.ac.uk>

!-- Modified to allow possibility of different heavy-quark masses.
!-- 25/06/2010 by Graeme Watt <Graeme.Watt(at)cern.ch>

!-- Fix "NaN" bug for q <= m_c when m_c^2 < 1.25 GeV^2.
!-- 25/01/2011 by Graeme Watt <Graeme.Watt(at)cern.ch>

subroutine MSTWevolve(x,Q,xpdf,xphoton)
  implicit none
  !
  include 'parmsetup.inc'
  character*16 name(nmxset)
  integer iset,mem,nset,nmem(nmxset),ndef(nmxset),mmem,imem
  COMMON/NAME/name,nmem,ndef,mmem
  double precision xpdf(-6:6),xvalence(6),xphoton, &
       &     x,Q,alfas,Qalfa,Eorder,Q2fit,MSTWALPHAS
  logical warn,fatal
  parameter(warn=.false.,fatal=.true.)
  !   Set warn=.true. to turn on warnings when extrapolating.
  !   Set fatal=.false. to return zero instead of terminating when
  !    invalid input values of x and q are used.
  integer ih,f,nhess,nx,nq,np,nqc0,nqb0,n,m,ip,io, &
       &     nExtraFlavours
  double precision xmin,xmax,qsqmin,qsqmax,mc2,mb2,eps, &
       &     qsq,xlog,qsqlog,res,res1,anom,ExtrapolateMSTWPDF, &
       &     InterpolateMSTWPDF
  parameter(nx=64,nq=48,np=12)
  parameter(xmin=1d-6,xmax=1d0,qsqmin=1d0,qsqmax=1d9,eps=1d-6)
  parameter(nhess=2*23)
  character dummyChar,dummyWord*50
  double precision gridx(nmxgridx),gridq(nmxgridq)
  integer ngridx,ngridq,jx,jq
  double precision ff(np,nx,nq)
  double precision qqorig(nq),qq(nq),xx(nx),cc(np,0:nhess,nx,nq,4,4,nmxset)
  double precision xxl(nx),qql(nq,nmxset,0:nhess)
  !   Store values of distance and tolerance along each eigenvector,
  !   heavy-quark masses and alphaS parameters in COMMON block.
  !   Allow the possibility of different alphaS parameters and
  !   heavy-quark masses for each member of the set.
  double precision distance(nmxset,0:nhess),tolerance(nmxset,0:nhess), &
       &     mCharm(nmxset,0:nhess),mBottom(nmxset,0:nhess), &
       &     alphaSQ0(nmxset,0:nhess),alphaSMZ(nmxset,0:nhess)
  integer alphaSorder(nmxset,0:nhess),alphaSnfmax(nmxset,0:nhess)
  common/mstwCommon/distance,tolerance, &
       &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
  double precision mCharmSave,mBottomSave,mTopSave
  double precision cmass(nmxset),bmass(nmxset),tmass(nmxset)
  common/masses_LHA/cmass,bmass,tmass
  data xx/1d-6,2d-6,4d-6,6d-6,8d-6, &
       &     1d-5,2d-5,4d-5,6d-5,8d-5, &
       &     1d-4,2d-4,4d-4,6d-4,8d-4, &
       &     1d-3,2d-3,4d-3,6d-3,8d-3, &
       &     1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
       &     .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
       &     .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
       &     .5d0,.525d0,.55d0,.575d0,.6d0,.625d0,.65d0,.675d0, &
       &     .7d0,.725d0,.75d0,.775d0,.8d0,.825d0,.85d0,.875d0, &
       &     .9d0,.925d0,.95d0,.975d0,1d0/
  data qqorig/1.d0, &
       &     1.25d0,1.5d0,0.d0,0.d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0, &
       &     1d1,1.2d1,0.d0,0.d0,2.6d1,4d1,6.4d1,1d2, &
       &     1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
       &     1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
       &     1.8d6,3.2d6,5.6d6,1d7,1.8d7,3.2d7,5.6d7,1d8, &
       &     1.8d8,3.2d8,5.6d8,1d9/
  save

  call getnset(iset)
  call getnmem(iset,imem)
  ih = imem

  qsq=q*q
  mc2=mCharm(iset,ih)**2
  mb2=mBottom(iset,ih)**2
  !   If mc2 < qsq < mc2+eps, then qsq = mc2+eps.
  if (qsq.gt.mc2.and.qsq.lt.mc2+eps) then
     qsq = mc2+eps
  end if
  !   If mb2 < qsq < mb2+eps, then qsq = mb2+eps.
  if (qsq.gt.mb2.and.qsq.lt.mb2+eps) then
     qsq = mb2+eps
  end if

  xlog=log10(x)
  qsqlog=log10(qsq)

  do f = 0, 13 ! loop over flavours

     res = 0.d0

     if (f.eq.0) then          ! gluon
        ip = 1
     else if (f.ge.1.and.f.le.5) then ! quarks
        ip = f+1
     else if (f.le.-1.and.f.ge.-5) then ! antiquarks
        ip = -f+1
     else if (f.ge.7.and.f.le.11) then ! valence quarks
        ip = f
     else if (f.eq.13) then    ! photon
        ip = 12
     else if (abs(f).ne.6.and.f.ne.12) then
        if (warn.or.fatal) print *,"Error in MSTWevolve: f = ",f
        if (fatal) stop
     end if

     if (x.le.0.d0.or.x.gt.xmax.or.q.le.0.d0) then

        if (warn.or.fatal) print *,"Error in MSTWevolve: x,qsq = ", &
             &        x,qsq
        if (fatal) stop

     else if (abs(f).eq.6.or.f.eq.12) then ! set top quarks to zero

        res = 0.d0

     else if (qsq.lt.qsqmin) then ! extrapolate to low Q^2

        if (warn) then
           print *, "Warning in MSTWevolve, extrapolating: f = ",f, &
                &           ", x = ",x,", q = ",q
        end if

        if (x.lt.xmin) then    ! extrapolate to low x

           res = ExtrapolateMSTWPDF(ip,np,ih,nhess,xlog, &
                & log10(qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           res1 = ExtrapolateMSTWPDF(ip,np,ih,nhess,xlog, &
                & log10(1.01D0*qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
              res = res - ExtrapolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                   & log10(qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
              res1 = res1 - ExtrapolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                   & log10(1.01D0*qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           end if

        else                   ! do usual interpolation

           res = InterpolateMSTWPDF(ip,np,ih,nhess,xlog, &
                & log10(qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           res1 = InterpolateMSTWPDF(ip,np,ih,nhess,xlog, &
                & log10(1.01D0*qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
              res = res - InterpolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                   & log10(qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
              res1 = res1 - InterpolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                   & log10(1.01D0*qsqmin),nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
           end if

        end if

        !   Calculate the anomalous dimension, dlog(xf)/dlog(qsq),
        !   evaluated at qsqmin.  Then extrapolate the PDFs to low
        !   qsq < qsqmin by interpolating the anomalous dimenion between
        !   the value at qsqmin and a value of 1 for qsq << qsqmin.
        !   If value of PDF at qsqmin is very small, just set
        !   anomalous dimension to 1 to prevent rounding errors.
        !   Impose minimum anomalous dimension of -2.5.
        if (abs(res).ge.1.D-5) then
           anom = max( -2.5D0, (res1-res)/res/0.01D0 )
        else
           anom = 1.D0
        end if
        res = res*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin)

     else if (x.lt.xmin.or.qsq.gt.qsqmax) then ! extrapolate

        if (warn) then
           print *, "Warning in MSTWevolve, extrapolating: f = ",f, &
                &           ", x = ",x,", q = ",q
        end if

        res = ExtrapolateMSTWPDF(ip,np,ih,nhess,xlog, &
             & qsqlog,nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))

        if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
           res = res - ExtrapolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                & qsqlog,nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
        end if

     else                      ! do usual interpolation

        res = InterpolateMSTWPDF(ip,np,ih,nhess,xlog, &
             & qsqlog,nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))

        if (f.le.-1.and.f.ge.-5) then ! antiquark = quark - valence
           res = res - InterpolateMSTWPDF(ip+5,np,ih,nhess,xlog, &
                & qsqlog,nx,nq,xxl,qql(1,iset,ih),cc(1,0,1,1,1,1,iset))
        end if

     end if

     if (f.ge.7.and.f.le.12) then
        xvalence(f-6) = res
     else if (f.eq.13) then
        xphoton = res
     else
        xpdf(f) = res
     end if

  end do

  do f = 1, 6
     xpdf(-f) = xpdf(f) - xvalence(f) ! antiquarks
  end do

  return 
  !                                                                       
  entry MSTWgetgrid(nset,ngridx,ngridq,gridx,gridq)
      do jx=1,nx
          gridx(jx)=xx(jx)
      enddo
      do jq=1,nq
          gridq(jq)=qq(jq)
      enddo
      ngridx=nx
      ngridq=nq        
  return

  entry MSTWread(nset)

  read(1,*) nmem(nset),ndef(nset)

  do ih = 0, nmem(nset)
     
     !   Read header containing heavy-quark masses and alphaS values.
     read(1,*)
     read(1,*)
     read(1,*) dummyChar,dummyWord,dummyWord,dummyChar, &
          &        distance(nset,ih),tolerance(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyChar,mCharm(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyChar,mBottom(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyChar,alphaSQ0(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyChar,alphaSMZ(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyWord,dummyChar, &
          &        alphaSorder(nset,ih),alphaSnfmax(nset,ih)
     read(1,*) dummyChar,dummyWord,dummyChar,nExtraFlavours
     read(1,*)
     read(1,*)

     !   Heavy-quark masses and alphaS values
     !   are stored in a COMMON block.
     !  WRITE(6,*) "mCharm = ",mCharm(nset,ih), &
     !       &     ", mBottom = ",mBottom(nset,ih)
     !  WRITE(6,*) "alphaS(Q0) = ",alphaSQ0(nset,ih),", alphaS(MZ) = ", &
     !       &     alphaSMZ(nset,ih),", alphaSorder = ",alphaSorder(nset,ih), &
     !       &     ", alphaSnfmax = ",alphaSnfmax(nset,ih)

     mc2=mCharm(nset,ih)**2
     mb2=mBottom(nset,ih)**2

     !   Check that the heavy-quark masses are sensible.
     !   Redistribute grid points if not in usual range.
     do m=1,nq
        qq(m) = qqorig(m)
     end do
     if (mc2.le.qq(1).or.mc2+eps.ge.qq(8)) then
        print *,"Error in MSTWevolve: invalid mCharm = ",mCharm(nset,ih)
        stop
     else if (mc2.lt.qq(2)) then
        nqc0=2
        qq(4)=qq(2)
        qq(5)=qq(3)
     else if (mc2.lt.qq(3)) then
        nqc0=3
        qq(5)=qq(3)
     else if (mc2.lt.qq(6)) then
        nqc0=4
     else if (mc2.lt.qq(7)) then
        nqc0=5
        qq(4)=qq(6)
     else
        nqc0=6
        qq(4)=qq(6)
        qq(5)=qq(7)
     end if
     if (mb2.le.qq(12).or.mb2+eps.ge.qq(17)) then
        print *,"Error in MSTWevolve: invalid mBottom = ",mBottom(nset,ih)
        stop
     else if (mb2.lt.qq(13)) then
        nqb0=13
        qq(15)=qq(13)
     else if (mb2.lt.qq(16)) then
        nqb0=14
     else
        nqb0=15
        qq(14)=qq(16)
     end if
     qq(nqc0)=mc2
     qq(nqc0+1)=mc2+eps
     qq(nqb0)=mb2
     qq(nqb0+1)=mb2+eps

     !   The nExtraFlavours variable is provided to aid compatibility
     !   with future grids where, for example, a photon distribution
     !   might be provided (cf. the MRST2004QED PDFs).
     if (nExtraFlavours.lt.0.or.nExtraFlavours.gt.1) then
        print *,"Error in MSTWevolve: invalid nExtraFlavours = ", &
             &           nExtraFlavours
        stop
     end if

     !   Now read in the grids from the grid file.
     do n=1,nx-1
        do m=1,nq
           if (nExtraFlavours.gt.0) then
              if (alphaSorder(nset,ih).eq.2) then ! NNLO
                 read(1,'(12(1pe12.4))',iostat=io) &
                      & (ff(ip,n,m),ip=1,12)
              else          ! LO or NLO
                 ff(10,n,m) = 0.d0 ! = chm-cbar
                 ff(11,n,m) = 0.d0 ! = bot-bbar
                 read(1,'(10(1pe12.4))',iostat=io) &
                      & (ff(ip,n,m),ip=1,9),ff(12,n,m)
              end if
           else             ! nExtraFlavours = 0
              if (alphaSorder(nset,ih).eq.2) then ! NNLO
                 ff(12,n,m) = 0.d0 ! = photon
                 read(1,'(11(1pe12.4))',iostat=io) &
                      & (ff(ip,n,m),ip=1,11)
              else          ! LO or NLO
                 ff(10,n,m) = 0.d0 ! = chm-cbar
                 ff(11,n,m) = 0.d0 ! = bot-bbar
                 ff(12,n,m) = 0.d0 ! = photon
                 read(1,'(9(1pe12.4))',iostat=io) &
                      & (ff(ip,n,m),ip=1,9)
              end if
           end if
           if (io.ne.0) then
              print *,"Error in MSTWevolve reading file"
              stop
           end if
        enddo
     enddo

     !   PDFs are identically zero at x = 1.
     do m=1,nq
        do ip=1,np
           ff(ip,nx,m)=0d0
        enddo
     enddo

     do n=1,nx
        xxl(n)=log10(xx(n))
     enddo
     do m=1,nq
        qql(m,nset,ih)=log10(qq(m))
     enddo

     !   Initialise all parton flavours.
     do ip=1,np
        call InitialiseMSTWPDF(ip,np,ih,nhess,nx,nq,nqc0,nqb0, &
             & xxl,qql(1,nset,ih),ff,cc(1,0,1,1,1,1,nset))
     enddo

  end do

  return 
  !                                                                       
  entry MSTWalfa(alfas,Qalfa)
  call getnset(iset)
  call getnmem(iset,imem)

  ! Check value of Qalfa is below appropriate thresholds.
  ! If not, redefine thresholds and reinitialise alphaS.
  IF (alphaSnfmax(iset,imem).EQ.5) THEN ! maximum of five flavours
     IF (Qalfa.GT.mTopSave) THEN
        mTopSave = 2.D0*Qalfa
        CALL MSTWINITALPHAS(alphaSorder(iset,imem), &
             &     1.D0,1.D0,alphaSQ0(iset,imem), &
             &     mCharmSave,mBottomSave,mTopSave)
     END IF
  ELSE IF (alphaSnfmax(iset,imem).EQ.4) THEN ! maximum of four flavours
     IF (Qalfa.GT.mBottomSave) THEN
        mBottomSave = 2.D0*Qalfa
        mTopSave = mBottomSave/0.9D0
        CALL MSTWINITALPHAS(alphaSorder(iset,imem), &
             &     1.D0,1.D0,alphaSQ0(iset,imem), &
             &     mCharmSave,mBottomSave,mTopSave)
     END IF
  ELSE IF (alphaSnfmax(iset,imem).EQ.3) THEN ! maximum of three flavours
     IF (Qalfa.GT.mCharmSave) THEN
        mCharmSave = 2.D0*Qalfa
        mBottomSave = mCharmSave*0.9D0/0.8D0
        mTopSave = mBottomSave/0.9D0
        CALL MSTWINITALPHAS(alphaSorder(iset,imem), &
             &     1.D0,1.D0,alphaSQ0(iset,imem), &
             &     mCharmSave,mBottomSave,mTopSave)
     END IF
  END IF
  alfas = MSTWALPHAS(Qalfa)
  return
  !
  entry MSTWinit(nset,Eorder,Q2fit)
  return
  !                                                                       
  entry MSTWpdf(mem)
  call getnset(iset)
  call setnmem(iset,mem)

  ! Call the initialisation routine with alpha_S(Q_0).
  IF (alphaSnfmax(iset,mem).EQ.5) THEN ! maximum of five flavours
     mCharmSave = mCharm(iset,mem)
     mBottomSave = mBottom(iset,mem)
     mTopSave = 1.D10
  ELSE IF (alphaSnfmax(iset,mem).EQ.4) THEN ! maximum of four flavours
     mCharmSave = mCharm(iset,mem)
     mBottomSave = 0.9D10
     mTopSave = 1.D10
  ELSE IF (alphaSnfmax(iset,mem).EQ.3) THEN ! maximum of three flavours
     mCharmSave = 0.8D10
     mBottomSave = 0.9D10
     mTopSave = 1.D10
  END IF

  ! Update masses stored in common/masses_LHA/.
  cmass(iset) = mCharm(iset,mem)
  bmass(iset) = mBottom(iset,mem)
  tmass(iset) = mTopSave
  
  CALL MSTWINITALPHAS(alphaSorder(iset,mem),1.D0,1.D0,alphaSQ0(iset,mem), &
       &     mCharmSave,mBottomSave,mTopSave)
  
  !  !   Check calculated value of alpha_S(M_Z) matches stored value.
  !  WRITE(6,'(" alphaS(MZ) = ",F7.5," = ",F7.5)') &
  !       &     MSTWALPHAS(91.1876D0),alphaSMZ(iset,mem)

  return 
  !                                                                       
END subroutine MSTWevolve



!----------------------------------------------------------------------

subroutine InitialiseMSTWPDF(ip,np,ih,nhess,nx,my,myc0,myb0, &
     &     xx,yy,ff,cc)
  implicit none
  integer nhess,ih,nx,my,myc0,myb0,j,k,l,m,n,ip,np
  double precision xx(nx),yy(my),ff(np,nx,my), &
       &     ff1(nx,my),ff2(nx,my),ff12(nx,my),ff21(nx,my), &
       &     yy0(4),yy1(4),yy2(4),yy12(4),z(16), &
       &     cl(16),cc(np,0:nhess,nx,my,4,4),iwt(16,16), &
       &     polderiv1,polderiv2,polderiv3,d1,d2,d1d2,xxd

      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
     &     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, &
     &     -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0, &
     &     2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0, &
     &     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, &
     &     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, &
     &     0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1, &
     &     0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1, &
     &     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0, &
     &     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0, &
     &     9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2, &
     &     -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2, &
     &     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0, &
     &     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0, &
     &     -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1, &
     &     4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/

      do m=1,my
         ff1(1,m)=polderiv1(xx(1),xx(2),xx(3), &
     &        ff(ip,1,m),ff(ip,2,m),ff(ip,3,m))
         ff1(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx), &
     &        ff(ip,nx-2,m),ff(ip,nx-1,m),ff(ip,nx,m))
         do n=2,nx-1
            ff1(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1), &
     &           ff(ip,n-1,m),ff(ip,n,m),ff(ip,n+1,m))
         enddo
      enddo

!--   Calculate the derivatives at qsq=mc2,mc2+eps,mb2,mb2+eps
!--   in a similar way as at the endpoints qsqmin and qsqmax.
      do n=1,nx
         do m=1,my
            if (myc0.eq.2.and.m.eq.1) then
               ff2(n,m)=(ff(ip,n,m+1)-ff(ip,n,m))/(yy(m+1)-yy(m))
            else if (myc0.eq.2.and.m.eq.2) then
               ff2(n,m)=(ff(ip,n,m)-ff(ip,n,m-1))/(yy(m)-yy(m-1))
            else if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff2(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2), &
     &              ff(ip,n,m),ff(ip,n,m+1),ff(ip,n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff2(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m), &
     &              ff(ip,n,m-2),ff(ip,n,m-1),ff(ip,n,m))
            else
               ff2(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1), &
     &              ff(ip,n,m-1),ff(ip,n,m),ff(ip,n,m+1))
            end if
         end do
      end do

!--   Calculate the cross derivatives (d/dx)(d/dy).
      do m=1,my
         ff12(1,m)=polderiv1(xx(1),xx(2),xx(3), &
     &        ff2(1,m),ff2(2,m),ff2(3,m))
         ff12(nx,m)=polderiv3(xx(nx-2),xx(nx-1),xx(nx), &
     &        ff2(nx-2,m),ff2(nx-1,m),ff2(nx,m))
         do n=2,nx-1
            ff12(n,m)=polderiv2(xx(n-1),xx(n),xx(n+1), &
     &           ff2(n-1,m),ff2(n,m),ff2(n+1,m))
         enddo
      enddo

!--   Calculate the cross derivatives (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            if (myc0.eq.2.and.m.eq.1) then
               ff21(n,m)=(ff1(n,m+1)-ff1(n,m))/(yy(m+1)-yy(m))
            else if (myc0.eq.2.and.m.eq.2) then
               ff21(n,m)=(ff1(n,m)-ff1(n,m-1))/(yy(m)-yy(m-1))
            else if (m.eq.1.or.m.eq.myc0+1.or.m.eq.myb0+1) then
               ff21(n,m)=polderiv1(yy(m),yy(m+1),yy(m+2), &
     &              ff1(n,m),ff1(n,m+1),ff1(n,m+2))
            else if (m.eq.my.or.m.eq.myc0.or.m.eq.myb0) then
               ff21(n,m)=polderiv3(yy(m-2),yy(m-1),yy(m), &
     &              ff1(n,m-2),ff1(n,m-1),ff1(n,m))
            else
               ff21(n,m)=polderiv2(yy(m-1),yy(m),yy(m+1), &
     &              ff1(n,m-1),ff1(n,m),ff1(n,m+1))
            end if
         end do
      end do

!--   Take the average of (d/dx)(d/dy) and (d/dy)(d/dx).
      do n=1,nx
         do m = 1, my
            ff12(n,m)=0.5*(ff12(n,m)+ff21(n,m))
         end do
      end do

      do n=1,nx-1
         do m=1,my-1
            d1=xx(n+1)-xx(n)
            d2=yy(m+1)-yy(m)
            d1d2=d1*d2
            
            yy0(1)=ff(ip,n,m)
            yy0(2)=ff(ip,n+1,m)
            yy0(3)=ff(ip,n+1,m+1)
            yy0(4)=ff(ip,n,m+1)
            
            yy1(1)=ff1(n,m)
            yy1(2)=ff1(n+1,m)
            yy1(3)=ff1(n+1,m+1)
            yy1(4)=ff1(n,m+1)
            
            yy2(1)=ff2(n,m)
            yy2(2)=ff2(n+1,m)
            yy2(3)=ff2(n+1,m+1)
            yy2(4)=ff2(n,m+1)
            
            yy12(1)=ff12(n,m)
            yy12(2)=ff12(n+1,m)
            yy12(3)=ff12(n+1,m+1)
            yy12(4)=ff12(n,m+1)
            
            do k=1,4
               z(k)=yy0(k)
               z(k+4)=yy1(k)*d1
               z(k+8)=yy2(k)*d2
               z(k+12)=yy12(k)*d1d2
            enddo
            
            do l=1,16
               xxd=0.d0
               do k=1,16
                  xxd=xxd+iwt(k,l)*z(k)
               enddo
               cl(l)=xxd
            enddo
            l=0
            do k=1,4
               do j=1,4
                  l=l+1
                  cc(ip,ih,n,m,k,j)=cl(l)
               enddo
            enddo
         enddo
      enddo
      return
      end

!----------------------------------------------------------------------

      double precision function InterpolateMSTWPDF(ip,np,ih,nhess,x,y, &
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,MSTWlocx,l,m,n,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4), &
     &     x,y,z,t,u

      n=MSTWlocx(xx,nx,x)
      m=MSTWlocx(yy,my,y)
      
      t=(x-xx(n))/(xx(n+1)-xx(n))
      u=(y-yy(m))/(yy(m+1)-yy(m))
      
      z=0.d0
      do l=4,1,-1
         z=t*z+((cc(ip,ih,n,m,l,4)*u+cc(ip,ih,n,m,l,3))*u &
     &        +cc(ip,ih,n,m,l,2))*u+cc(ip,ih,n,m,l,1)
      enddo

      InterpolateMSTWPDF = z

      return
      end

!----------------------------------------------------------------------

      double precision function ExtrapolateMSTWPDF(ip,np,ih,nhess,x,y, &
     &     nx,my,xx,yy,cc)
      implicit none
      integer ih,nx,my,nhess,MSTWlocx,n,m,ip,np
      double precision xx(nx),yy(my),cc(np,0:nhess,nx,my,4,4), &
     &     x,y,z,f0,f1,z0,z1,InterpolateMSTWPDF

      n=MSTWlocx(xx,nx,x)           ! 0: below xmin, nx: above xmax
      m=MSTWlocx(yy,my,y)           ! 0: below qsqmin, my: above qsqmax

!   If extrapolation in small x only:
      if (n.eq.0.and.m.gt.0.and.m.lt.my) then
         f0 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(1),y,nx,my,xx,yy,cc)
         f1 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(2),y,nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
         end if
!   If extrapolation into large q only:
      else if (n.gt.0.and.m.eq.my) then
         f0 = InterpolateMSTWPDF(ip,np,ih,nhess,x,yy(my),nx,my,xx,yy,cc)
         f1 = InterpolateMSTWPDF(ip,np,ih,nhess,x,yy(my-1),nx,my,xx,yy,cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
     &           (y-yy(my)))
         else
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
!   If extrapolation into large q AND small x:
      else if (n.eq.0.and.m.eq.my) then
         f0 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(1),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(1),yy(my-1),nx,my,xx,yy, &
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
     &           (y-yy(my)))
         else
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         f0 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(2),yy(my),nx,my,xx,yy,cc)
         f1 = InterpolateMSTWPDF(ip,np,ih,nhess,xx(2),yy(my-1),nx,my,xx,yy, &
     &        cc)
         if (f0.gt.1.d-3.and.f1.gt.1.d-3) then
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
     &           (y-yy(my)))
         else
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
         end if
         if (z0.gt.1.d-3.and.z1.gt.1.d-3) then
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
         else
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
         end if
      else
         print *,"Error in ExtrapolateMSTWPDF"
         stop
      end if

      ExtrapolateMSTWPDF = z      

      return
      end

!----------------------------------------------------------------------

      integer function MSTWlocx(xx,nx,x)
!   returns an integer j such that x lies inbetween xx(j) and xx(j+1).
!   nx is the length of the array with xx(nx) the highest element.
      implicit none
      integer nx,jl,ju,jm
      double precision x,xx(nx)
      if(x.eq.xx(1)) then
         MSTWlocx=1
         return
      endif
      if(x.eq.xx(nx)) then
         MSTWlocx=nx-1  
         return
      endif
      ju=nx+1
      jl=0
    1 if((ju-jl).le.1) goto 2
      jm=(ju+jl)/2
      if(x.ge.xx(jm)) then
         jl=jm
      else
         ju=jm
      endif
      goto 1
    2 MSTWlocx=jl
      return
      end

!----------------------------------------------------------------------

      double precision function polderiv1(x1,x2,x3,y1,y2,y3)
!--   returns the estimate of the derivative at x1 obtained by a
!--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv1=(x3*x3*(y1-y2)+2.d0*x1*(x3*(-y1+y2)+x2*(y1-y3)) &
     &     +x2*x2*(-y1+y3)+x1*x1*(-y2+y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv2(x1,x2,x3,y1,y2,y3)
!--   returns the estimate of the derivative at x2 obtained by a
!--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv2=(x3*x3*(y1-y2)-2.d0*x2*(x3*(y1-y2)+x1*(y2-y3)) &
     &     +x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

      double precision function polderiv3(x1,x2,x3,y1,y2,y3)
!--   returns the estimate of the derivative at x3 obtained by a
!--   polynomial interpolation using the three points (x_i,y_i).
      implicit none
      double precision x1,x2,x3,y1,y2,y3
      polderiv3=(x3*x3*(-y1+y2)+2.d0*x2*x3*(y1-y3)+x1*x1*(y2-y3) &
     &     +x2*x2*(-y1+y3)+2.d0*x1*x3*(-y2+y3))/ &
     &     ((x1-x2)*(x1-x3)*(x2-x3))
      return
      end

!----------------------------------------------------------------------

!----------------------------------------------------------------------
!--   Stand-alone code for alpha_s cannibalised (with permission)
!--   from Andreas Vogt's QCD-PEGASUS package (hep-ph/0408244).
!--   The running coupling alpha_s is obtained at N^mLO (m = 0,1,2,3)
!--   by solving the renormalisation group equation in the MSbar scheme
!--   by a fourth-order Runge-Kutta integration.  Transitions from
!--   n_f to n_f+1 flavours are made when the factorisation scale
!--   mu_f equals the pole masses m_h (h = c,b,t).  At exactly
!--   the thresholds m_{c,b,t}, the number of flavours n_f = {3,4,5}.
!--   The top quark mass should be set to be very large to evolve with
!--   a maximum of five flavours.  The factorisation scale mu_f may be
!--   a constant multiple of the renormalisation scale mu_r.  The input
!--   factorisation scale mu_(f,0) should be less than or equal to
!--   the charm quark mass.  However, if it is greater than the
!--   charm quark mass, the value of alpha_s at mu_(f,0) = 1 GeV will
!--   be found using a root-finding algorithm.
!--
!--   Example of usage.
!--   First call the initialisation routine (only needed once):
!--
!--    IORD = 2                  ! perturbative order (N^mLO,m=0,1,2,3)
!--    FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2
!--    MUR = 1.D0                ! input mu_r in GeV
!--    ASMUR = 0.5D0             ! input value of alpha_s at mu_r
!--    MC = 1.4D0                ! charm quark mass
!--    MB = 4.75D0               ! bottom quark mass
!--    MT = 1.D10                ! top quark mass
!--    CALL MSTWINITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
!--
!--   Then get alpha_s at a renormalisation scale mu_r with:
!--
!--    MUR = 100.D0              ! renormalisation scale in GeV
!--    ALFAS = MSTWALPHAS(MUR)
!--
!----------------------------------------------------------------------
!--   Comments to Graeme Watt <watt(at)hep.ucl.ac.uk>
!----------------------------------------------------------------------

      SUBROUTINE MSTWINITALPHAS(IORD, FR2, MUR, ASMUR, MC, MB, MT)
!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   MUR = input renormalisation scale (in GeV) for alpha_s.
!--   ASMUR = input value of alpha_s at the renormalisation scale MUR.
!--   MC,MB,MT = heavy-quark masses in GeV.
      IMPLICIT NONE
      INTEGER IORD,IORDc,MAXF,MODE
      DOUBLE PRECISION FR2,MUR,ASMUR,MC,MB,MT,EPS,A,B,MSTWDZEROX, &
     &     R0c,FR2c,MURc,ASMURc,MCc,MBc,MTc,FINDALPHASR0,R0,ASI
      COMMON / MSTWDZEROXcommon / FR2c,MURc,ASMURc,MCc,MBc,MTc,R0c,IORDc
      PARAMETER(EPS=1.D-10,MAXF=10000,MODE=1)
      EXTERNAL FINDALPHASR0

      IF (MUR*sqrt(FR2).LE.MC) THEN ! Check that MUF <= MC.
         R0 = MUR
         ASI = ASMUR
      ELSE                      ! Solve for alpha_s at R0 = 1 GeV.
!--   Copy variables to common block.
         R0c = 1.D0/sqrt(FR2)
         IORDc = IORD
         FR2c = FR2
         MURc = MUR
         ASMURc = ASMUR
         MCc = MC
         MBc = MB
         MTc = MT
!--   Now get alpha_s(R0) corresponding to alpha_s(MUR).
         A = 0.02D0              ! lower bound for alpha_s(R0)
         B = 2.00D0              ! upper bound for alpha_s(R0)
         R0 = R0c
         ASI = MSTWDZEROX(A,B,EPS,MAXF,FINDALPHASR0,MODE)
      END IF

      CALL MSTWINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)

      RETURN
      END

!----------------------------------------------------------------------

!--   Find the zero of this function using MSTWDZEROX.
      DOUBLE PRECISION FUNCTION FINDALPHASR0(ASI)
      IMPLICIT NONE
      INTEGER IORD
      DOUBLE PRECISION FR2, R0, ASI, MC, MB, MT, MUR, ASMUR, MSTWALPHAS
      COMMON / MSTWDZEROXcommon / FR2, MUR, ASMUR, MC, MB, MT, R0, IORD

      CALL MSTWINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)
      FINDALPHASR0 = MSTWALPHAS(MUR) - ASMUR ! solve equal to zero

      RETURN
      END
      
!----------------------------------------------------------------------

      SUBROUTINE MSTWINITALPHASR0(IORD, FR2, R0, ASI, MC, MB, MT)
!--   IORD = 0 (LO), 1 (NLO), 2 (NNLO), 3 (NNNLO).
!--   FR2 = ratio of mu_f^2 to mu_r^2 (must be a fixed value).
!--   R0 = input renormalisation scale (in GeV) for alphas_s.
!--   ASI = input value of alpha_s at the renormalisation scale R0.
!--   MC,MB,MT = heavy-quark masses in GeV.
!--   Must have R0*sqrt(FR2) <= MC to call this subroutine.
      IMPLICIT NONE
      INTEGER IORD,NAORD,NASTPS,IVFNS,NFF
      DOUBLE PRECISION FR2,R0,ASI,MC,MB,MT,LOGFR,R20, &
     &     PI,ZETA,CF,CA,TR,AS0,M20,MC2,MB2,MT2
      PARAMETER(PI = 3.14159265358979D0)

      COMMON / RZETA  / ZETA(6)
      COMMON / COLOUR / CF, CA, TR
      COMMON / ASINP  / AS0, M20
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / VARFLV / IVFNS
      COMMON / NFFIX  / NFF
      COMMON / FRRAT  / LOGFR

!
! ..QCD colour factors
!
      CA = 3.D0
      CF = 4./3.D0
      TR = 0.5D0
!
! ..The lowest integer values of the Zeta function
!
      ZETA(1) = 0.57721566490153D0
      ZETA(2) = 1.644934066848226D0
      ZETA(3) = 1.202056903159594D0
      ZETA(4) = 1.082323233711138D0
      ZETA(5) = 1.036927755143370D0
      ZETA(6) = 1.017343061984449D0

      IVFNS = 1                 ! variable flavour-number scheme (VFNS)
!      IVFNS = 0                 ! fixed flavour-number scheme (FFNS)
      NFF = 4                   ! number of flavours for FFNS
      NAORD = IORD              ! perturbative order of alpha_s
      NASTPS = 20               ! num. steps in Runge-Kutta integration
      R20 = R0**2               ! input renormalisation scale
      MC2 = MC**2               ! mu_f^2 for charm threshold
      MB2 = MB**2               ! mu_f^2 for bottom threshold
      MT2 = MT**2               ! mu_f^2 for top threshold
      LOGFR = LOG(FR2)          ! log of ratio of mu_f^2 to mu_r^2
      M20 = R20 * FR2           ! input factorisation scale

!
! ..Stop some nonsense
!
      IF ( (IVFNS .EQ. 0) .AND. (NFF .LT. 3) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 0) .AND. (NFF .GT. 5) ) THEN
         WRITE (6,*) 'Wrong flavour number for FFNS evolution. STOP'
         STOP
      END IF
!     
      IF ( NAORD .GT. 3 ) THEN
         WRITE (6,*) 'Specified order in a_s too high. STOP' 
         STOP
      END IF
!
      IF ( (IVFNS .NE. 0) .AND. (FR2 .GT. 4.001D0) ) THEN
         WRITE (6,*) 'Too low mu_r for VFNS evolution. STOP'
         STOP
      END IF
!
      IF ( (IVFNS .EQ. 1) .AND. (M20 .GT. MC2) ) THEN
         WRITE (6,*) 'Too high mu_0 for VFNS evolution. STOP'
         STOP
      END IF
!     
      IF ( (ASI .GT. 2.D0) .OR. (ASI .LT. 2.D-2) ) THEN
         WRITE (6,*) 'alpha_s out of range. STOP'
         STOP
      END IF
!     
      IF ( (IVFNS .EQ. 1) .AND. (MC2 .GT. MB2) ) THEN
         WRITE (6,*) 'Wrong charm-bottom mass hierarchy. STOP'
         STOP
      END IF
      IF ( (IVFNS .EQ. 1) .AND. (MB2 .GT. MT2) ) THEN
         WRITE (6,*) 'Wrong bottom-top mass hierarchy. STOP'
         STOP
      END IF
!

!--   Store the beta function coefficients in a COMMON block.
      CALL BETAFCT

!--   Store a_s = alpha_s(mu_r^2)/(4 pi) at the input scale R0.
      AS0 = ASI / (4.D0* PI)

!--   Store alpha_s at the heavy flavour thresholds in a COMMON block.
       IF (IVFNS .NE. 0) THEN
          CALL EVNFTHR (MC2, MB2, MT2)
       END IF

      RETURN
      END

!----------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION MSTWALPHAS(MUR)
      IMPLICIT NONE
      INTEGER NFF,IVFNS,NF
      DOUBLE PRECISION PI,LOGFR,AS0,M20,ASC,M2C,ASB,M2B,AST,M2T,M2,MUR, &
     &     R2,ASI,ASF,R20,R2T,R2B,R2C,AS
      PARAMETER ( PI = 3.14159265358979D0 )
!
! ..Input common blocks 
! 
       COMMON / NFFIX  / NFF
       COMMON / VARFLV / IVFNS 
       COMMON / FRRAT  / LOGFR
       COMMON / ASINP  / AS0, M20
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

       R2 = MUR**2
       M2 = R2 * EXP(+LOGFR)
       IF (IVFNS .EQ. 0) THEN
!
!   Fixed number of flavours
!
          NF  = NFF
          R20 = M20 * R2/M2
          ASI = AS0
          ASF = AS (R2, R20, AS0, NF)
!
       ELSE
!
! ..Variable number of flavours
!
          IF (M2 .GT. M2T) THEN
             NF = 6
             R2T = M2T * R2/M2
             ASI = AST
             ASF = AS (R2, R2T, AST, NF)
!
          ELSE IF (M2 .GT. M2B) THEN
             NF = 5
             R2B = M2B * R2/M2
             ASI = ASB
             ASF = AS (R2, R2B, ASB, NF)
!     
          ELSE IF (M2 .GT. M2C) THEN
             NF = 4
             R2C = M2C * R2/M2
             ASI = ASC
             ASF = AS (R2, R2C, ASC, NF)
!     
          ELSE
             NF = 3
             R20 = M20 * R2/M2
             ASI = AS0
             ASF = AS (R2, R20, AS0, NF)
!       
          END IF
!
       END IF
!
! ..Final value of alpha_s
!
       MSTWALPHAS = 4.D0*PI*ASF

       RETURN
       END
!
! =================================================================av==


! =====================================================================
!
! ..The threshold matching of the QCD coupling in the MS(bar) scheme,  
!    a_s = alpha_s(mu_r^2)/(4 pi),  for  NF -> NF + 1  active flavours 
!    up to order a_s^4 (NNNLO).
!
! ..The value  ASNF  of a_s for NF flavours at the matching scale, the 
!    logarithm  LOGRH = ln (mu_r^2/m_H^2) -- where m_H is the pole mass
!    of the heavy quark -- and  NF  are passed as arguments to the 
!    function  ASNF1.  The order of the expansion  NAORD  (defined as 
!    the 'n' in N^nLO) is provided by the common-block  ASPAR.
!
! ..The matching coefficients are inverted from Chetyrkin, Kniehl and
!    Steinhauser, Phys. Rev. Lett. 79 (1997) 2184. The QCD colour
!    factors have been hard-wired in these results. The lowest integer 
!    values of the Zeta function are given by the common-block  RZETA.
!
! =====================================================================
!
!
      DOUBLE PRECISION FUNCTION ASNF1 (ASNF, LOGRH, NF)
!
      IMPLICIT NONE
      INTEGER NF, NAORD, NASTPS, PRVCLL, K1, K2
      DOUBLE PRECISION ASNF,LOGRH,ZETA,CMC,CMCI30,CMCF30,CMCF31, &
     &     CMCI31,ASP,LRHP

      DIMENSION CMC(3,0:3)
!
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
      COMMON / ASPAR  / NAORD, NASTPS
      COMMON / RZETA  / ZETA(6)
!
! ..Variables to be saved for the next call
!
      SAVE CMC, CMCI30, CMCF30, CMCF31, CMCI31, PRVCLL
!
! ---------------------------------------------------------------------
!
! ..The coupling-constant matching coefficients (CMC's) up to NNNLO 
!   (calculated and saved in the first call of this routine)
!
       IF (PRVCLL .NE. 1) THEN
!
         CMC(1,0) =  0.D0
         CMC(1,1) =  2./3.D0
!
         CMC(2,0) = 14./3.D0
         CMC(2,1) = 38./3.D0
         CMC(2,2) =  4./9.D0  
!
         CMCI30 = + 80507./432.D0 * ZETA(3) + 58933./1944.D0 &
     &            + 128./3.D0 * ZETA(2) * (1.+ DLOG(2.D0)/3.D0)
         CMCF30 = - 64./9.D0 * (ZETA(2) + 2479./3456.D0)
         CMCI31 =   8941./27.D0
         CMCF31 = - 409./27.D0
         CMC(3,2) = 511./9.D0
         CMC(3,3) = 8./27.D0
!
         PRVCLL = 1
!
       END IF
!
! ---------------------------------------------------------------------
!
! ..The N_f dependent CMC's, and the alpha_s matching at order NAORD 
!
       CMC(3,0) = CMCI30 + NF * CMCF30
       CMC(3,1) = CMCI31 + NF * CMCF31
!
       ASNF1 = ASNF
       IF (NAORD .EQ. 0) GOTO 1
       ASP   = ASNF
!
       DO 11 K1 = 1, NAORD 
         ASP = ASP * ASNF
         LRHP = 1.D0
!
       DO 12 K2 = 0, K1
         ASNF1 = ASNF1 + ASP * CMC(K1,K2) * LRHP
         LRHP = LRHP * LOGRH
!
  12   CONTINUE
  11   CONTINUE
!
! ---------------------------------------------------------------------
!
   1   RETURN
       END

!
! =================================================================av==
!
! ..The subroutine  EVNFTHR  for the evolution of  a_s = alpha_s/(4 pi)
!    from a three-flavour initial scale to the four- to six-flavour
!    thresholds (identified with the squares of the corresponding quark
!    masses).  The results are written to the common-block  ASFTHR.
!
! ..The input scale  M20 = mu_(f,0)^2  and the corresponding value 
!    AS0  of a_s  are provided by  ASINP.  The fixed scale logarithm
!    LOGFR = ln (mu_f^2/mu_r^2) is specified in  FRRAT.  The alpha_s
!    matching is done by the function ASNF1.
!
! =====================================================================
!
!
       SUBROUTINE EVNFTHR (MC2, MB2, MT2)
!
       IMPLICIT NONE
       DOUBLE PRECISION MC2, MB2, MT2, M20, M2C, M2B, M2T, R20, R2C, &
     &                  R2B, R2T, AS, ASNF1, AS0, ASC, ASB, AST, &
     &                  ASC3, ASB4, AST5, LOGFR, SC, SB, ST
!
! ---------------------------------------------------------------------
! 
! ..Input common blocks
!  
       COMMON / ASINP  / AS0, M20
       COMMON / FRRAT  / LOGFR
!
! ..Output common blocks
!
       COMMON / ASFTHR / ASC, M2C, ASB, M2B, AST, M2T

! ---------------------------------------------------------------------
!
! ..Coupling constants at and evolution distances to/between thresholds
! 
       R20 = M20 * EXP(-LOGFR)
!
! ..Charm
!
       M2C  = MC2
       R2C  = M2C * R20/M20
       ASC3 = AS (R2C, R20, AS0, 3)
       SC   = LOG (AS0 / ASC3)
       ASC  = ASNF1 (ASC3, -LOGFR, 3)
!
! ..Bottom 
!
       M2B  = MB2
       R2B  = M2B * R20/M20
       ASB4 = AS (R2B, R2C, ASC, 4)
       SB   = LOG (ASC / ASB4)
       ASB  = ASNF1 (ASB4, -LOGFR, 4)
!
! ..Top
!
       M2T  = MT2
       R2T  = M2T * R20/M20
       AST5 = AS (R2T, R2B, ASB, 5)
       ST   = LOG (ASB / AST5)
       AST  = ASNF1 (AST5, -LOGFR, 5)

       RETURN
       END

!
! =================================================================av==
!
! ..The running coupling of QCD,  
!
!         AS  =  a_s  =  alpha_s(mu_r^2)/(4 pi),
!
!    obtained by integrating the evolution equation for a fixed number
!    of massless flavours  NF.  Except at leading order (LO),  AS  is 
!    obtained using a fourth-order Runge-Kutta integration.
!
! ..The initial and final scales  R20  and  R2,  the value  AS0  at
!    R20, and  NF  are passed as function arguments.  The coefficients 
!    of the beta function up to  a_s^5 (N^3LO)  are provided by the 
!    common-block  BETACOM.  The order of the expansion  NAORD (defined
!    as the 'n' in N^nLO) and the number of steps  NASTPS  for the 
!    integration beyond LO are given by the common-block  ASPAR.
!
! =====================================================================
!
!
      DOUBLE PRECISION FUNCTION AS (R2, R20, AS0, NF)
!
      IMPLICIT NONE
      INTEGER NFMIN, NFMAX, NF, NAORD, NASTPS, K1
      DOUBLE PRECISION R2, R20, AS0, SXTH, BETA0, BETA1, BETA2, BETA3, &
     &     FBETA1,FBETA2,FBETA3,A,LRRAT,DLR,XK0,XK1,XK2,XK3
      PARAMETER (NFMIN = 3, NFMAX = 6)
      PARAMETER ( SXTH = 0.166666666666666D0 )
!
! ---------------------------------------------------------------------
!
! ..Input common-blocks 
!
       COMMON / ASPAR  / NAORD, NASTPS
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX), &
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)
!
! ..The beta functions FBETAn at N^nLO for n = 1, 2, and 3
!
       FBETA1(A) = - A**2 * ( BETA0(NF) + A *   BETA1(NF) )
       FBETA2(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF) &
     &                        + A * BETA2(NF) ) )
       FBETA3(A) = - A**2 * ( BETA0(NF) + A * ( BETA1(NF) &
     &                        + A * (BETA2(NF) + A * BETA3(NF)) ) )
!
! ---------------------------------------------------------------------
!
! ..Initial value, evolution distance and step size
!
       AS = AS0
       LRRAT = LOG (R2/R20)
       DLR = LRRAT / NASTPS
!
! ..Solution of the evolution equation depending on  NAORD
!   (fourth-order Runge-Kutta beyond the leading order)
!
       IF (NAORD .EQ. 0) THEN
!
         AS = AS0 / (1.+ BETA0(NF) * AS0 * LRRAT)
!
       ELSE IF (NAORD .EQ. 1) THEN
!
       DO 2 K1 = 1, NASTPS
         XK0 = DLR * FBETA1 (AS)
         XK1 = DLR * FBETA1 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA1 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA1 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  2    CONTINUE
!
       ELSE IF (NAORD .EQ. 2) THEN
!
       DO 3 K1 = 1, NASTPS
         XK0 = DLR * FBETA2 (AS)
         XK1 = DLR * FBETA2 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA2 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA2 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  3    CONTINUE
!  
       ELSE IF (NAORD .EQ. 3) THEN
!
       DO 4 K1 = 1, NASTPS
         XK0 = DLR * FBETA3 (AS)
         XK1 = DLR * FBETA3 (AS + 0.5 * XK0)
         XK2 = DLR * FBETA3 (AS + 0.5 * XK1)
         XK3 = DLR * FBETA3 (AS + XK2)
         AS = AS + SXTH * (XK0 + 2.* XK1 + 2.* XK2 + XK3)
  4    CONTINUE
       END IF
!
! ---------------------------------------------------------------------
!
       RETURN
       END

!
! =================================================================av==
!
! ..The subroutine BETAFCT for the coefficients  BETA0...BETA3  of the 
!    beta function of QCD up to order alpha_s^5 (N^3LO), normalized by 
!
!        d a_s / d ln mu_r^2  =  - BETA0 a_s^2 - BETA1 a_s^3 - ... 
!
!    with  a_s = alpha_s/(4*pi). 
!
! ..The MSbar coefficients are written to the common-block BETACOM for 
!   NF = 3...6  (parameters NFMIN, NFMAX) quark flavours.
!
! ..The factors CF, CA and TF  are taken from the common-block  COLOUR.
!    Beyond NLO the QCD colour factors are hard-wired in this routine,
!    and the numerical coefficients are truncated to six digits.
!
! =====================================================================
!
!
       SUBROUTINE BETAFCT
!
       IMPLICIT DOUBLE PRECISION (A - Z)
       INTEGER NFMIN, NFMAX, NF
       PARAMETER (NFMIN = 3, NFMAX = 6)
!
! ---------------------------------------------------------------------
!
! ..Input common-block
!
       COMMON / COLOUR / CF, CA, TR
!
! ..Output common-block
!
       COMMON / BETACOM   / BETA0 (NFMIN:NFMAX), BETA1 (NFMIN:NFMAX), &
     &                   BETA2 (NFMIN:NFMAX), BETA3 (NFMIN:NFMAX)

!
! ---------------------------------------------------------------------
!
! ..The full LO and NLO coefficients 
!
       B00 =  11./3.D0 * CA
       B01 =  -4./3.D0 * TR
       B10 =  34./3.D0 * CA**2
       B11 = -20./3.D0 * CA*TR - 4.* CF*TR
!
! ..Flavour-number loop and output to the array
!
       DO 1 NF = NFMIN, NFMAX
!
       BETA0(NF) = B00 + B01 * NF
       BETA1(NF) = B10 + B11 * NF
!
       BETA2(NF) = 1428.50 - 279.611 * NF + 6.01852 * NF**2
       BETA3(NF) = 29243.0 - 6946.30 * NF + 405.089 * NF**2 &
     &             + 1.49931 * NF**3
!
! ---------------------------------------------------------------------
!
  1    CONTINUE
!
       RETURN
       END
!
! =================================================================av==


!--   G.W. MSTWDZEROX taken from CERNLIB to find the zero of a function.
      DOUBLE PRECISION FUNCTION MSTWDZEROX(A0,B0,EPS,MAXF,F,MODE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!     Based on
!
!        J.C.P. Bus and T.J. Dekker, Two Efficient Algorithms with
!        Guaranteed Convergence for Finding a Zero of a Function,
!        ACM Trans. Math. Software 1 (1975) 330-345.
!
!        (MODE = 1: Algorithm M;    MODE = 2: Algorithm R)
      CHARACTER*80 ERRTXT
      LOGICAL LMT
      DIMENSION IM1(2),IM2(2),LMT(2)
      PARAMETER (Z1 = 1, HALF = Z1/2)
      DATA IM1 /2,3/, IM2 /-1,3/
      MSTWDZEROX = 0.D0             ! G.W. to prevent compiler warning
      IF(MODE .NE. 1 .AND. MODE .NE. 2) THEN
       C=0
       WRITE(ERRTXT,101) MODE
       WRITE(6,*) ERRTXT
       GOTO 99
      ENDIF
      FA=F(B0)
      FB=F(A0)
      IF(FA*FB .GT. 0) THEN
       C=0
       WRITE(ERRTXT,102) A0,B0
       WRITE(6,*) ERRTXT
       GOTO 99
      ENDIF
      ATL=ABS(EPS)
      B=A0
      A=B0
      LMT(2)=.TRUE.
      MF=2
    1 C=A
      FC=FA
    2 IE=0
    3 IF(ABS(FC) .LT. ABS(FB)) THEN
       IF(C .NE. A) THEN
        D=A
        FD=FA
       END IF
       A=B
       B=C
       C=A
       FA=FB
       FB=FC
       FC=FA
      END IF
      TOL=ATL*(1+ABS(C))
      H=HALF*(C+B)
      HB=H-B
      IF(ABS(HB) .GT. TOL) THEN
       IF(IE .GT. IM1(MODE)) THEN
        W=HB
       ELSE
        TOL=TOL*SIGN(Z1,HB)
        P=(B-A)*FB
        LMT(1)=IE .LE. 1
        IF(LMT(MODE)) THEN
         Q=FA-FB
         LMT(2)=.FALSE.
        ELSE
         FDB=(FD-FB)/(D-B)
         FDA=(FD-FA)/(D-A)
         P=FDA*P
         Q=FDB*FA-FDA*FB
        END IF
        IF(P .LT. 0) THEN
         P=-P
         Q=-Q
        END IF
        IF(IE .EQ. IM2(MODE)) P=P+P
        IF(P .EQ. 0 .OR. P .LE. Q*TOL) THEN
         W=TOL
        ELSEIF(P .LT. HB*Q) THEN
         W=P/Q
        ELSE
         W=HB
        END IF
       END IF
       D=A
       A=B
       FD=FA
       FA=FB
       B=B+W
       MF=MF+1
       IF(MF .GT. MAXF) THEN
        WRITE(6,*) "Error in MSTWDZEROX: TOO MANY FUNCTION CALLS"
        GOTO 99
       ENDIF
       FB=F(B)
       IF(FB .EQ. 0 .OR. SIGN(Z1,FC) .EQ. SIGN(Z1,FB)) GOTO 1
       IF(W .EQ. HB) GOTO 2
       IE=IE+1
       GOTO 3
      END IF
      MSTWDZEROX=C
   99 CONTINUE
      RETURN
  101 FORMAT('Error in MSTWDZEROX: MODE = ',I3,' ILLEGAL')
  102 FORMAT('Error in MSTWDZEROX: F(A) AND F(B) HAVE THE SAME SIGN, A = ', &
     &     1P,D15.8,', B = ',D15.8)
      END

! ---------------------------------------------------------------------
