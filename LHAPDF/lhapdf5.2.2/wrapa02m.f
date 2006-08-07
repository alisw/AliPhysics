      subroutine A02Mevolve(xb,Q,fset)
      implicit none
      include 'parmsetup.inc'
      integer npdf,npar,kschem,i,k,n,m,kx,nxbb
      integer nxb,nq,np,nvar,pdfmem,imem
      parameter(nxb=99,nq=20,np=9,nvar=17)
      integer nexp(0:np)
      data nexp / 0, 3, 4, 5, 5, 5, 5, 5, 5, 5 /
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem
      common/NAME/name,nmem,ndef,mem
c      integer iset,iimem
c      common/SET/iset,iimem
      integer nset
      real*8 f(0:nvar,nxb,nq+1,0:np),xx(nxb)
      real*8 fsp(nxb),bs(nxb),cs(nxb),ds(nxb)
      real*8 bsp(0:nvar,nxb,nq+1,0:np),csp(0:nvar,nxb,nq+1,0:np)
     +      ,dsp(0:nvar,nxb,nq+1,0:np)
      real*8 x1,xd,del,dels,delx,delx1,xlog1
      real*8 pdfs(0:np),fset(-6:6),alfas
      real*8 x,Q,xb,q2,xmin,xmax,qsq,qsqmin,qsqmax,b,ss
      real*8 aa,f0,fp,fm
      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/

      save f,npdf,npar,pdfmem,dels,delx,x1,delx1,xlog1,nxbb,xx
     +,fsp,bsp,csp,dsp
c      save 
      
      q2=Q*Q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      if(xb.lt.xmin.or.xb.gt.xmax)       print 98,x
  99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ')
  98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ')

      x=max(xb,xmin)
      x=min(xb,xmax)
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)

      if (x.gt.x1) then
        xd=(1d0-x1)**2-(1d0-x)**2
        n=int(xd/delx1)+nxbb
      else
        xd=dlog(x)-xlog1
        n=nxbb+int(xd/DELX)-1
      end if
      aa=x-xx(n)

      ss=dlog(dlog(qsq/0.04d0))-dlog(dlog(qsqmin/0.04d0))
      m=int(ss/dels)+1
      b=ss/dels-dble(m)+1.d0
 
      k=pdfmem

      do i=0,npdf
        f0=f(k,n,m,i) + aa*bsp(k,n,m,i) + aa**2*csp(k,n,m,i) 
     +              + aa**3*dsp(k,n,m,i)
        fp=f(k,n,m+1,i) + aa*bsp(k,n,m+1,i) + aa**2*csp(k,n,m+1,i)
     +                + aa**3*dsp(k,n,m+1,i)
        if (m.ge.2) then 
          fm=f(k,n,m-1,i) + aa*bsp(k,n,m-1,i) + aa**2*csp(k,n,m-1,i)
     +                   +aa**3*dsp(k,n,m-1,i)
          pdfs(i)=fm*b*(b-1d0)/2d0 + f0*(1d0-b**2) + fp*b*(b+1d0)/2d0
        else 
          pdfs(i)= f0*(1d0-b) + fp*b
         end if
        pdfs(i) = pdfs(i)*(1d0-x)**dble(nexp(i))
      end do

      fset(-6)=pdfs(9)
      fset(-5)=pdfs(8)
      fset(-4)=pdfs(7)
      fset(-3)=pdfs(5)
c--reversed mrs 7/7/04 due
      fset(-1)=pdfs(6)
      fset(-2)=pdfs(4)
      fset(0)=pdfs(3)
c--reversed mrs 7/7/04 due
      fset(2)=pdfs(1)+pdfs(4)
      fset(1)=pdfs(2)+pdfs(6)
      fset(3)=pdfs(5)
      fset(4)=pdfs(7)
      fset(5)=pdfs(8)
      fset(6)=pdfs(9)
      return
*
      entry A02Malfa(alfas,Q)
      q2=Q*Q
      if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
      qsq=max(q2,qsqmin)
      qsq=min(q2,qsqmax)

      ss=log(log(qsq/0.04))-log(log(qsqmin/0.04))
      m=int(ss/dels)+1
      b=ss/dels-m+1

      k=pdfmem
      alfas=(1d0-b)*f(k,1,m,0)+b*f(k,1,m+1,0)
      return
*
      entry A02Mread(nset)
c following fix because members are 0-nvar
c      nmem = nvar + 1

      nmem(nset) = nvar
      ndef(nset) = 0


      read(1,*) npdf
c      print *,'number of members',npdf
      npar=nvar
      do k=0,npar
         do n=1,nxb-1
            do m=1,nq
               read(1,100) (f(k,n,m,i),i=0,npdf)
c	       print 100,(f(k,n,m,i),i=0,npdf)
            enddo
         enddo
      enddo
  100 format (13f11.5)

      nxbb=nxb/2
      x1=0.3d0
      xlog1=dlog(x1)
      delx=(dlog(x1)-dlog(xmin))/dble(nxbb-1)
      DELX1=(1.d0-x1)**2/dble(nxbb+1)

*...X GRID
      do kx=1,nxbb
        xx(kx)=dexp(dlog(xmin)+delx*dble(kx-1))
      end do
      do kx=nxbb+1,nxb-1
        xx(kx)=1.d0-dsqrt(dabs((1.d0-x1)**2-delx1*dble(kx-nxbb)))
      end do
      xx(nxb)=1d0

      do k=0,npar
      do i=0,npdf
        do m=1,nq
          if (i.ne.0) then 
            f(k,nxb,m,i)=0d0
          else 
            f(k,nxb,m,i)=f(k,nxb-1,m,i)
          end if
          do n=1,nxb-1
            f(k,n,m,i)=f(k,n,m,i)/(1d0-xx(n))**nexp(i)
          end do
          do n=1,nxb
            fsp(n)=f(k,n,m,i)
          end do
          call a02mspline (nxb,xx,fsp,bs,cs,ds)
          do n=1,nxb
            bsp(k,n,m,i)=bs(n)
            csp(k,n,m,i)=cs(n)
            dsp(k,n,m,i)=ds(n)
          end do
        end do
      end do
      end do

      return
*
      entry A02Minit

      dels=(dlog(dlog(qsqmax/0.04d0))-
     +      dlog(dlog(qsqmin/0.04d0)))/dble(nq-1)


      return
*
      entry A02Mpdf(imem)
      pdfmem=imem
      if ((pdfmem.lt.0).or.(pdfmem.gt.npar)) then
         write(*,*) 'A02M PDF set:'
         write(*,*) 'PDF member out of range:'
         write(*,*) 'member = ',pdfmem,'    member range = (0,',npar,')'
         stop
      endif
      return
      end

* ---------------------------------------------------------------------
      SUBROUTINE A02MSPLINE(N,X,Y,B,C,D)
* ---------------------------------------------------------------------
* CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
* INTERPOLATION SUBROUTINES ARE TAKEN FROM
* G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
* COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*
      IMPLICIT REAL*8(A-H,O-Z)

      DIMENSION X(N), Y(N), B(N), C(N), D(N)
*
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
      B(N)=(Y(N)-Y(NM1))/D(NM1)
     1     +D(NM1)*(C(NM1)+2.0D0*C(N))
      DO 240 K=1,NM1
         B(K)=(Y(K+1)-Y(K))/D(K)
     1        -D(K)*(C(K+1)+2.0D0*C(K))
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
