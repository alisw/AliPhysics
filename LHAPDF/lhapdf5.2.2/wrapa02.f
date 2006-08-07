      subroutine A02evolve(xb,Q,fset)
      implicit none
      include 'parmsetup.inc'
      integer npdf,npar,kschem,i,k,n,m,nxbb
      integer nxb,nq,np,nvar,pdfmem,imem
      parameter(nxb=99,nq=20,np=9,nvar=15)
      character*16 name(nmxset)
      integer nmem(nmxset),ndef(nmxset),mem
      common/NAME/name,nmem,ndef,mem
      integer nset
      real*4 f(0:nvar,nxb,nq+1,0:np)
      real*8 x1,xd,del,dels,delx,delx1,xlog1
      real*8 pdfs(np),fset(-6:6),alfas
      real*8 x,Q,xb,q2,xmin,xmax,qsq,qsqmin,qsqmax,a,b,ss
      data xmin,xmax,qsqmin,qsqmax/1d-7,1d0,0.8d0,2d8/

      save f,npdf,npar,pdfmem,dels,delx,x1,delx1,xlog1,nxbb

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
        xd=(1-x1)**2-(1-x)**2
        n=int(xd/delx1)+nxbb
        a=xd/delx1-n+nxbb
      else
        xd=log(x)-xlog1
        n=nxbb+int(xd/DELX)-1
        a=xd/delx-n+nxbb
      end if

      ss=log(log(qsq/0.04))-log(log(qsqmin/0.04))
      m=int(ss/dels)+1
      b=ss/dels-m+1

      k=pdfmem

      do i=1,npdf
        if (n.gt.1.and.m.gt.1.and.n.ne.49) then 
        pdfs(i)= f(k,n,m,i)*(1+a*b-a**2-b**2) 
     +         + f(k,n+1,m+1,i)*a*b
     +         + f(k,n+1,m,i)*a*(a-2*b+1)/2. 
     +         + f(k,n,m+1,i)*b*(b-2*a+1)/2.
     +         +  f(k,n-1,m,i)*a*(a-1)/2. 
     +         + f(k,n,m-1,i)*b*(b-1)/2.
        else
        pdfs(i)= (1-a)*(1-b)*f(k,n,m,i) 
     .         + (1-a)*b*f(k,n,m+1,i)
     .         + a*(1-b)*f(k,n+1,m,i) 
     .         + a*b*f(k,n+1,m+1,i)
        end if
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
      entry A02alfa(alfas,Q)
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
      entry A02read(nset)
c following fix because members are 0-nvar
c      nmem = nvar + 1
      nmem(nset) = nvar
      ndef(nset) = 0
      read(1,*) npdf
      npar=nvar
      do k=0,npar
         do n=1,nxb-1
            do m=1,nq
               read(1,100) (f(k,n,m,i),i=0,npdf)
            enddo
         enddo
      enddo
  100 format (13f11.5)
      return
*
      entry A02init

      dels=(log(log(qsqmax/0.04))-log(log(qsqmin/0.04)))/(nq-1)

      nxbb=nxb/2
      x1=0.3
      xlog1=log(x1)
      delx=(log(x1)-log(xmin))/(nxbb-1)
      DELX1=(1-x1)**2/(nxbb+1)

      do i=1,npdf
         do m=1,nq
            do k=1,npar
               f(k,nxb,m,i)=0d0
            end do
         end do
         do m=1,nq
            do k=1,npar
               f(k,nxb,m,0)=f(k,nxb-1,m,0)
            end do
         end do
      end do
      return
*
      entry A02pdf(imem)
      pdfmem=imem
      if ((pdfmem.lt.0).or.(pdfmem.gt.npar)) then
         write(*,*) 'A02 PDF set:'
         write(*,*) 'PDF member out of range:'
         write(*,*) 'member = ',pdfmem,'    member range = (0,',npar,')'
         stop
      endif
      return
      end
