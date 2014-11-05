! -*- F90 -*-


      subroutine MRSTqedevolve(x,Q,pdf,photon)
      implicit real*8(a-h,o-z)
      include 'parmsetup.inc'
      character*16 name(nmxset)
      character*80 line
      character*512 setpath
      integer nmem(nmxset),ndef(nmxset),mmem
!      integer member(nmxset)
      integer nset,iset
      common/NAME/name,nmem,ndef,mmem
      parameter(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26,          &
     &nhess=0)
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
      real*8 pdf(-6:6),photon,phot
      real*8 f1(nx,nq)                                                  &
     &,f2(nx,nq)                                                        &
     &,f3(nx,nq)                                                        &
     &,f4(nx,nq)                                                        &
     &,f5(nx,nq)                                                        &
     &,f6(nx,nq)                                                        &
     &,f7(nx,nq)                                                        &
     &,f8(nx,nq)                                                        &
     &,f9(nx,nq)                                                        &
     &,fc(nx,nqc),fb(nx,nqb)
      real*8 qq(nq),xx(nx),                                             &
     &cc1(0:nhess,nx,nq,4,4,nmxset),cc2(0:nhess,nx,nq,4,4,nmxset),      &
     &cc3(0:nhess,nx,nq,4,4,nmxset),cc4(0:nhess,nx,nq,4,4,nmxset),      &
     &cc6(0:nhess,nx,nq,4,4,nmxset),cc8(0:nhess,nx,nq,4,4,nmxset),      &
     &cc9(0:nhess,nx,nq,4,4,nmxset),                                    &
     &ccc(0:nhess,nx,nqc,4,4,nmxset),ccb(0:nhess,nx,nqb,4,4,nmxset)
      real*8 xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
      data xx/1d-5,2d-5,4d-5,6d-5,8d-5,                                 &
     &              1d-4,2d-4,4d-4,6d-4,8d-4,                           &
     &              1d-3,2d-3,4d-3,6d-3,8d-3,                           &
     &              1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2,               &
     &           .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0,     &
     &           .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0,     &
     &           .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0,        &
     &           .8d0,.9d0,1d0/
      data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1,       &
     &        1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2,                          &
     &        1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4,          &
     &        1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6,              &
     &        1.8d6,3.2d6,5.6d6,1d7/
      data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
      save

      xsave=x
      qsq = q*q
      q2save=qsq

      xlog=dlog(x)
      qsqlog=dlog(qsq)

      call getnset(iset)
!      imem=member(iset)
      call getnmem(iset,imem)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc1(0,1,1,1,1,iset),upv)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc2(0,1,1,1,1,iset),dnv)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc3(0,1,1,1,1,iset),glu)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc4(0,1,1,1,1,iset),usea)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc6(0,1,1,1,1,iset),str)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc8(0,1,1,1,1,iset),dsea)
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,cc9(0,1,1,1,1,iset),photon)

      chm=0.d0
      if(qsq.gt.emc2) then
      call jeppe2(0,xlog,qsqlog,nx,nqc,xxl,qqlc,ccc(0,1,1,1,1,iset),chm)
      endif

      bot=0.d0
      if(qsq.gt.emb2) then
      call jeppe2(0,xlog,qsqlog,nx,nqb,xxl,qqlb,ccb(0,1,1,1,1,iset),bot)
      endif

      pdf(0)  = glu
      pdf(1)  = dnv+dsea
      pdf(-1) = dsea
      pdf(2)  = upv+usea
      pdf(-2) = usea
      pdf(3)  = str
      pdf(-3) = str
      pdf(4)  = chm
      pdf(-4) = chm
      pdf(5)  = bot
      pdf(-5) = bot
      pdf(6) = 0.0d0
      pdf(-6) = 0.0d0

      x=xsave
      qsq=q2save
      return
 
      entry MRSTqedgetgrid(nset,ngridx,ngridq,gridx,gridq)
      do jx=1,nx
          gridx(jx)=xx(jx)
      enddo
      do jq=1,nq
          gridq(jq)=qq(jq)
      enddo
      ngridx=nx
      ngridq=nq        
      return
!
      entry MRSTqedread(nset)
! - dummy read in to get to End: (stream 1 is still open)
      read(1,*)nmem(nset),ndef(nset)

      do i=0,nmem(nset)
        do n=1,nx-1
          do m=1,nq
            read(1,'(a)'),line
          enddo
        enddo
      enddo

      do n=1,nx
      xxl(n)=dlog(xx(n))
      enddo
      do m=1,nq
      qql(m)=dlog(qq(m))
      enddo

      return
!
      entry MRSTqedpdf(mem)
        call getnset(iset)
        call setnmem(iset,mem)

! have to reopen stream 1
        call getsetpath(setpath)
        open(1,file=setpath(1:len_trim(setpath)),action='READ')
        line = ''
        do while (line(2:11).ne.'Evolution:')
           read(1,'(a)'),line
        enddo
        read(1,'(a)'),line
        read(1,'(a)'),line

        read(1,*)nmem(iset),ndef(iset)

! - dummy read up to the member requested
      do i=0,mem-1
        do n=1,nx-1
          do m=1,nq
            read(1,'(a)')line
          enddo
        enddo
      enddo


!- read in the data of the requested member
        do 20 n=1,nx-1
        do 20 m=1,nq
        read(1,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),                      &
     &                  f5(n,m),f7(n,m),f6(n,m),f8(n,m),f9(n,m)
! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea 9-photon
   20 continue
!      write(*,*)'PDF set ',nm,' first element ',f1(1,1)

      do 40 m=1,nq
      f1(nx,m)=0.d0
      f2(nx,m)=0.d0
      f3(nx,m)=0.d0
      f4(nx,m)=0.d0
      f5(nx,m)=0.d0
      f6(nx,m)=0.d0
      f7(nx,m)=0.d0
      f8(nx,m)=0.d0
      f9(nx,m)=0.d0
   40 continue

      call jeppe1(0,nx,nq,xxl,qql,f1,cc1(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f2,cc2(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f3,cc3(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f4,cc4(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f6,cc6(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f8,cc8(0,1,1,1,1,iset))
      call jeppe1(0,nx,nq,xxl,qql,f9,cc9(0,1,1,1,1,iset))

      emc2=2.045
      emb2=18.5

      do 44 m=1,nqc
      qqlc(m)=qql(m+nqc0)
      do 44 n=1,nx
      fc(n,m)=f5(n,m+nqc0)
   44 continue
      qqlc(1)=dlog(emc2)
      call jeppe1(0,nx,nqc,xxl,qqlc,fc,ccc(0,1,1,1,1,iset))

      do 45 m=1,nqb
      qqlb(m)=qql(m+nqb0)
      do 45 n=1,nx
      fb(n,m)=f7(n,m+nqb0)
   45 continue
      qqlb(1)=dlog(emb2)
      call jeppe1(0,nx,nqb,xxl,qqlb,fb,ccb(0,1,1,1,1,iset))

      close(1)
   50 format(9f10.5)


      return
!
      entry MRSTqedalfa(nflav,alfas,Qalfa)
        call getnset(iset)
        call alphamrs(nflav,alfas,Qalfa)
      return
!
      entry MRSTqedinit(Eorder,Q2fit)
      return
!
      END
