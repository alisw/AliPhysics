      subroutine ABM11evolve(xb,Q,fset)
      implicit none
      include 'parmsetup.inc'
      integer npdf,npar,kschem,i,k,kk,n,m,kx,nxbb,nxb1,nxb2
      integer nxb,nq,np,nvar,pdfmem,imem,nvar2
      parameter(nxb=99,nq=20,np=9,nvar=0,nvar2=25)
      character*16 name(nmxset)
      character*80 line 
      character*512 setpath 
      integer nmem(nmxset),ndef(nmxset),mem
      common/NAME/name,nmem,ndef,mem
      double precision gridx(nmxgridx),gridq(nmxgridq)
      integer ngridx,ngridq,jx,jq
!      integer iset,iimem
!      common/SET/iset,iimem
      integer nset,iset
      real*8 f(0:nvar,nxb,nq+1,0:np),xx(nxb)
      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 bsp(0:nvar,nxb,nq+1,0:np),csp(0:nvar,nxb,nq+1,0:np) &
     &      ,dsp(0:nvar,nxb,nq+1,0:np)
      real*8 x1,xd,del,dels,delx,delx1,xlog1
      real*8 pdfs(0:np),fset(-6:6),alfas
      real*8 x,Q,xb,q2,xmin,xmax,xmax1,qsq,qsqmin,qsqmax,b,ss
      real*8 aa,f0,fp,fm
      real*8 pdfs0(0:np),pdfsp(0:np),pdfsm(0:np),q2plus,q2minus,q2zero
      real*8 q2mat,delq2,f2p,f2pp,aq2,bq2,cq2
      integer ios
      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/
      data q2mat /1d0/

      save f,npdf,npar,pdfmem,dels,delx,x1,delx1,xlog1,nxbb,xx &
     &,fsp,bsp,csp,dsp,nxb1,nxb2,xmax1
!      save 
      
      q2=Q*Q
      if(q2.gt.qsqmax) print 99,q2
      if(xb.lt.xmin.or.xb.gt.xmax)    print 98,xb

  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ',g12.3)
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ',g12.3)

      x=max(xb,xmin)
      x=min(x,xmax)
      qsq=min(q2,qsqmax)

!      if (x.gt.x1) then
!        xd=(1d0-x1)**2-(1d0-x)**2
!        n=int(xd/delx1)+nxbb
!      else
!        xd=dlog(x)-xlog1
!        n=nxbb+int(xd/DELX)-1
!      end if

      do n=1,nxb-1
        if (x.lt.xx(n+1)) goto 300
        end do
 300  aa=x-xx(n)

!      print *,x,xx(n),xx(n+1),aa
      !k=pdfmem
      k=0

      if (qsq.ge.q2mat) then 
        ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
        m=int(ss/dels)+1
        b=ss/dels-dble(m)+1.d0
 
        do i=1,npdf
          f0=f(k,n,m,i) + aa*bsp(k,n,m,i) + aa**2*csp(k,n,m,i) & 
     &                  + aa**3*dsp(k,n,m,i)
          fp=f(k,n,m+1,i) + aa*bsp(k,n,m+1,i) + aa**2*csp(k,n,m+1,i) &
     &                    + aa**3*dsp(k,n,m+1,i)
          if (m.ge.2) then 
            fm=f(k,n,m-1,i) + aa*bsp(k,n,m-1,i) + aa**2*csp(k,n,m-1,i) &
     &                      +aa**3*dsp(k,n,m-1,i)
            pdfs(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
          else 
            pdfs(i)= f0*(1d0-b) + fp*b
           end if
        end do
      else 
        delq2=q2mat*0.1

        do i=1,npdf

          q2plus=q2mat+delq2
          ss=dlog(dlog(q2plus/0.04d0))-dlog(dlog(qsqmin/0.04d0))
          m=int(ss/dels)+1
          b=ss/dels-dble(m)+1.d0

          f0=f(k,n,m,i) + aa*bsp(k,n,m,i) + aa**2*csp(k,n,m,i)  &
     &                  + aa**3*dsp(k,n,m,i)
          fp=f(k,n,m+1,i) + aa*bsp(k,n,m+1,i) + aa**2*csp(k,n,m+1,i) &
     &                    + aa**3*dsp(k,n,m+1,i)
          if (m.ge.2) then 
            fm=f(k,n,m-1,i) + aa*bsp(k,n,m-1,i) + aa**2*csp(k,n,m-1,i) &
     &                      +aa**3*dsp(k,n,m-1,i)
            pdfsp(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
          else 
            pdfsp(i)= f0*(1d0-b) + fp*b
           end if

          q2minus=q2mat-delq2
          ss=dlog(dlog(q2minus/0.04d0))-dlog(dlog(qsqmin/0.04d0))
          m=int(ss/dels)+1
          b=ss/dels-dble(m)+1.d0

          f0=f(k,n,m,i) + aa*bsp(k,n,m,i) + aa**2*csp(k,n,m,i) & 
     &                  + aa**3*dsp(k,n,m,i)
          fp=f(k,n,m+1,i) + aa*bsp(k,n,m+1,i) + aa**2*csp(k,n,m+1,i) &
     &                    + aa**3*dsp(k,n,m+1,i)
          if (m.ge.2) then 
            fm=f(k,n,m-1,i) + aa*bsp(k,n,m-1,i) + aa**2*csp(k,n,m-1,i) &
     &                      +aa**3*dsp(k,n,m-1,i)
            pdfsm(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
          else 
            pdfsm(i)= f0*(1d0-b) + fp*b
          end if

          q2zero=q2mat
          ss=dlog(dlog(q2zero/0.04d0))-dlog(dlog(qsqmin/0.04d0))
          m=int(ss/dels)+1
          b=ss/dels-dble(m)+1.d0

          f0=f(k,n,m,i) + aa*bsp(k,n,m,i) + aa**2*csp(k,n,m,i) & 
     &                  + aa**3*dsp(k,n,m,i)
          fp=f(k,n,m+1,i) + aa*bsp(k,n,m+1,i) + aa**2*csp(k,n,m+1,i) &
     &                    + aa**3*dsp(k,n,m+1,i)
          if (m.ge.2) then 
            fm=f(k,n,m-1,i) + aa*bsp(k,n,m-1,i) + aa**2*csp(k,n,m-1,i) &
     &                      +aa**3*dsp(k,n,m-1,i)
            pdfs0(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
          else 
            pdfs0(i)= f0*(1d0-b) + fp*b
           end if

          f2p=(pdfsp(i)-pdfsm(i))/2./delq2
          f2pp=(pdfsp(i)+pdfsm(i)-2*pdfs0(i))/delq2**2
          cq2=(f2pp*q2mat**2-2*f2p*q2mat+2*pdfs0(i))/2/q2mat**3
          bq2=f2pp/2.-3*cq2*q2mat
          aq2=f2p-2*bq2*q2mat-3*cq2*q2mat**2
          pdfs(i)=aq2*q2 + bq2*q2**2 + cq2*q2**3
!          print *,i,xb,q2,aq2,bq2,cq2
        end do
      end if

!      fset(-6)=pdfs(9)
!      fset(-5)=pdfs(8)
!      fset(-4)=pdfs(7)
      fset(-3)=pdfs(5)
!--reversed mrs 7/7/04 due
      fset(-1)=pdfs(6)
      fset(-2)=pdfs(4)
      fset(0)=pdfs(3)
!--reversed mrs 7/7/04 due
      fset(2)=pdfs(1)+pdfs(4)
      fset(1)=pdfs(2)+pdfs(6)
      fset(3)=pdfs(5)
!      fset(4)=pdfs(7)
!      fset(5)=pdfs(8)
!      fset(6)=pdfs(9)

      do kk=7,npdf
        fset(kk-3)=pdfs(kk)
        fset(-kk+3)=pdfs(kk)
      end do
      do kk=npdf+1,9
        fset(kk-3)=0.
        fset(-kk+3)=0.
      end do

      return
!
      entry ABM11alfa(alfas,Q)
      q2=Q*Q
      if(q2.gt.qsqmax) print 99,q2
      qsq=max(q2,q2mat)
      qsq=min(qsq,qsqmax)

      ss=log(log(qsq/0.04))-log(log(qsqmin/0.04))
      m=int(ss/dels)+1
      b=ss/dels-m+1

      k=pdfmem

      fp=f(k,1,m+1,0)
      f0=f(k,1,m,0)
      if (m.ge.2) then 
        fm=f(k,1,m-1,0)
        alfas=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
      else 
        alfas= f0*(1d0-b) + fp*b
      end if

!      print *,q2,q2mat,m,alfas


      return
!
      entry ABM11getgrid(nset,ngridx,ngridq,gridx,gridq)
     
      ngridx=nxb
      do jx=1,nxb
          gridx(jx)=xx(jx)
      enddo

      ngridq=nq
      do jq=1,ngridq
          gridq(jq)= 0.04*exp(exp( &
     & log(log(qsqmin/0.04))+(float(jq-1)/19)*( log(log(qsqmax/0.04))-log(log(qsqmin/0.04)) ) & 
     & ))
      enddo 
       
      return
!
      entry ABM11read(nset)
! following fix because members are 0-nvar
!      nmem = nvar + 1

      nmem(nset) = nvar
      ndef(nset) = 0

! - dummy read in to get to End: (stream 1 is still open)               
      read(1,fmt="(i2,i2)") npdf,npar
!      print *,'number of members',npdf,npar
      if (npar.eq.0) npar=nvar2

      do k=0,npar
         do n=1,nxb-1
            do m=1,nq
               read(1,*) (f(k,n,m,i),i=0,npdf)
!	       print 100,k,n,m,i,(f(k,n,m,i),i=0,npdf)
            enddo
         enddo
      enddo
  100 format (4i4,13f11.5)

      return
!
      entry ABM11init

      dels=(dlog(dlog(qsqmax/0.04d0))- &
     &      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)


      return
!
      entry ABM11pdf(imem)
      pdfmem=imem
      if ((pdfmem.lt.0).or.(pdfmem.gt.npar)) then
         write(*,*) 'ABM11 PDF set:'
         write(*,*) 'PDF member out of range:'
         write(*,*) 'member = ',pdfmem,'    member range = (0,',npar,')'
         stop
      endif

      call getnset(iset)
      nmem(iset) = npar 
      ndef(iset) = 0 
      
! have to reopen stream 1                                               
      call getsetpath(setpath) 
      open(1,file=setpath(1:len_trim(setpath)),action='READ') 
      line = '' 
      do while (line(2:11).ne.'Evolution:') 
         read(1,'(a)') line 
      enddo 

      read(1,'(a)'),line 
      read(1,'(a)'),line 

      read(1,fmt="(i2,i2)") npdf,npar
      if (npar.eq.0) npar=nvar2

      do k=0,imem-1
         do n=1,nxb-1
            do m=1,nq
               read(1,*) (f(0,n,m,i),i=0,npdf)
            enddo
         enddo
      enddo
      
      do n=1,nxb-1 
         do m=1,nq 
            read(1,*) (f(0,n,m,i),i=0,npdf) 
         enddo 
      enddo 

!...X GRID

      nxb1=30
      nxb2=nxb-nxb1
      x1=0.3d0
      xmax1=0.99d0
      delx=(log(x1)-log(xmin))/(nxb1-1)
      DELX1=(log(1-xmax1)-log(1-x1))/(nxb2-1)

      do kx=1,nxb1
        xx(kx)=exp(log(xmin)+delx*(kx-1))
      end do

      do kx=nxb-1,nxb1,-1
        xx(kx)=1-exp(log(1-xmax1)+delx1*(kx+1-nxb))
      end do

!      nxbb=nxb/2
!      x1=0.3d0
!      xlog1=dlog(x1)
!      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
!      DELX1=(1.d0-x1)**2/dble(nxbb+1)
!      do kx=1,nxbb
!        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
!      end do
!      do kx=nxbb+1,nxb-1
!        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
!      end do

      xx(nxb)=1d0

      do i=0,npdf
        do m=1,nq
          if (i.ne.0) then 
            f(0,nxb,m,i)=0d0
          else 
            f(0,nxb,m,i)=f(0,nxb-1,m,i)
          end if
          do n=1,nxb-1
          end do
          do n=1,nxb
            fsp(n)=f(0,n,m,i)
          end do
          call ABM11spline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            bsp(0,n,m,i)=bs(n)
            csp(0,n,m,i)=cs(n)
            dsp(0,n,m,i)=ds(n)
          end do
        end do
      end do
      close(1)
      return
      end

! ---------------------------------------------------------------------
      SUBROUTINE ABM11SPLINE(N,X,Y,B,C,D)
! ---------------------------------------------------------------------
! CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
! INTERPOLATION SUBROUTINES ARE TAKEN FROM
! G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
! COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
!
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION X(N), Y(N), B(N), C(N), D(N)
!
      NM1=N-1
      IF(N.LT.2) RETURN
      IF(N.LT.3) GO TO 250
      D(1)=X(2)-X(1)
      C(2)=(Y(2)-Y(1))/D(1)
      DO 210 K=2,NM1
         D(K)=X(K+1)-X(K)
         B(K)=2.0D0*(D(K-1)+D(K))
         C(K+1)=(Y(K+1)-Y(K))/D(K)
         C(K)=C(K+1)-C(K)
  210 CONTINUE
      B(1)=-D(1)
      B(N)=-D(N-1)
      C(1)=0.0D0
      C(N)=0.0D0
      IF(N.EQ.3) GO TO 215
      C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
      C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
      C(1)=C(1)*D(1)**2.0D0/(X(4)-X(1))
      C(N)=-C(N)*D(N-1)**2.0D0/(X(N)-X(N-3))
 215  CONTINUE
      DO 220 K=2,N
         T=D(K-1)/B(K-1)
         B(K)=B(K)-T*D(K-1)
         C(K)=C(K)-T*C(K-1)
 220  CONTINUE
      C(N)=C(N)/B(N)
      DO 230 IB=1,NM1
         K=N-IB
         C(K)=(C(K)-D(K)*C(K+1))/B(K)
 230  CONTINUE
      B(N)=(Y(N)-Y(NM1))/D(NM1) &
     &     +D(NM1)*(C(NM1)+2.0D0*C(N))
      DO 240 K=1,NM1
         B(K)=(Y(K+1)-Y(K))/D(K) &
     &        -D(K)*(C(K+1)+2.0D0*C(K))
         D(K)=(C(K+1)-C(K))/D(K)
         C(K)=3.0D0*C(K)
 240  CONTINUE
      C(N)=3.0D0*C(N)
      D(N)=D(N-1)
      RETURN
 250  CONTINUE
      B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0D0
      D(1)=0.0D0
      B(2)=B(1)
      C(2)=0.0D0
      D(2)=0.0D0
      RETURN
      END
