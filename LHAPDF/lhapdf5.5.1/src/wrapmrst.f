! -*- F90 -*-


      subroutine MRSTevolve(x,Q,pdf) 
      implicit real*8(a-h,o-z) 
      include 'parmsetup.inc' 
      character*16 name(nmxset) 
      integer nmem(nmxset),ndef(nmxset),mmem 
!      integer member(nmxset)                                           
      integer nset,iset 
      common/NAME/name,nmem,ndef,mmem 
      parameter(nx=49,nq=37,np=8,nqc0=2,nqb0=11,nqc=35,nqb=26,          &
     &nhess=30)                                                         
      real*8 pdf(-6:6) 
      real*8 f1(nx,nq)                                                  &
     &,f2(nx,nq)                                                        &
     &,f3(nx,nq)                                                        &
     &,f4(nx,nq)                                                        &
     &,f5(nx,nq)                                                        &
     &,f6(nx,nq)                                                        &
     &,f7(nx,nq)                                                        &
     &,f8(nx,nq)                                                        &
     &,fc(nx,nqc),fb(nx,nqb)                                            
      real*8 qq(nq),xx(nx),                                             &
     &cc1(0:nhess,nx,nq,4,4,nmxset),cc2(0:nhess,nx,nq,4,4,nmxset),      &
     &cc3(0:nhess,nx,nq,4,4,nmxset),cc4(0:nhess,nx,nq,4,4,nmxset),      &
     &cc6(0:nhess,nx,nq,4,4,nmxset),cc8(0:nhess,nx,nq,4,4,nmxset),      &
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
                                                                        
                                ! G.W. 24/04/2008                       
      if (qsq.lt.qsqmin) then 
         qsqlog=dlog(1.01D0*qsqmin) 
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &   cc1(0,1,1,1,1,iset),upv1)                                      
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &         cc2(0,1,1,1,1,iset),dnv1)                                
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &         cc3(0,1,1,1,1,iset),glu1)                                
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &   cc4(0,1,1,1,1,iset),usea1)                                     
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &   cc6(0,1,1,1,1,iset),str1)                                      
         call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                    &
     &   cc8(0,1,1,1,1,iset),dsea1)                                     
         qsqlog=dlog(qsqmin) 
      end if 
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc1(0,1,1,1,1,iset),upv)                                          
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc2(0,1,1,1,1,iset),dnv)                                          
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc3(0,1,1,1,1,iset),glu)                                          
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc4(0,1,1,1,1,iset),usea)                                         
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc6(0,1,1,1,1,iset),str)                                          
      call jeppe2(imem,xlog,qsqlog,nx,nq,xxl,qql,                       &
     &cc8(0,1,1,1,1,iset),dsea)                                         
                                ! G.W. 24/04/2008                       
      if (qsq.lt.qsqmin) then 
!--   Calculate the anomalous dimension, dlog(xf)/dlog(qsq),            
!--   evaluated at qsqmin.  Then extrapolate the PDFs to low            
!--   qsq < qsqmin by interpolating the anomalous dimenion between      
!--   the value at qsqmin and a value of 1 for qsq << qsqmin.           
!--   If value of PDF at qsqmin is very small, just set                 
!--   anomalous dimension to 1 to prevent rounding errors.              
         anom = (upv1-upv)/upv/0.01D0 
         if (abs(upv).lt.1.D-4) anom = 1.D0 
         upv = upv*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
         anom = (dnv1-dnv)/dnv/0.01D0 
         if (abs(dnv).lt.1.D-4) anom = 1.D0 
         dnv = dnv*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
         anom = (glu1-glu)/glu/0.01D0 
         if (abs(glu).lt.1.D-4) anom = 1.D0 
         glu = glu*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
         anom = (usea1-usea)/usea/0.01D0 
         if (abs(usea).lt.1.D-4) anom = 1.D0 
         usea = usea*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
         anom = (str1-str)/str/0.01D0 
         if (abs(str).lt.1.D-4) anom = 1.D0 
         str = str*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
         anom = (dsea1-dsea)/dsea/0.01D0 
         if (abs(dsea).lt.1.D-4) anom = 1.D0 
         dsea = dsea*(qsq/qsqmin)**(anom*qsq/qsqmin+1.D0-qsq/qsqmin) 
      end if 
                                                                        
                                                                        
                                                                        
      chm=0.d0 
      if(qsq.gt.emc2) then 
      call jeppe2(imem,xlog,qsqlog,nx,nqc,xxl,qqlc,                     &
     &ccc(0,1,1,1,1,iset),chm)                                          
      endif 
                                                                        
      bot=0.d0 
      if(qsq.gt.emb2) then 
      call jeppe2(imem,xlog,qsqlog,nx,nqb,xxl,qqlb,                     &
     &ccb(0,1,1,1,1,iset),bot)                                          
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
      pdf(6)  = 0.0d0 
      pdf(-6) = 0.0d0 
                                                                        
      x=xsave 
      qsq=q2save 
      return 
!                                                                       
      entry MRSTread(nset) 
      read(1,*)nmem(nset),ndef(nset) 
!      print *,nmem(nset),ndef(nset)                                    
!      do nm = 0,nmem-1                                                 
      do nm = 0,nmem(nset) 
        do 20 n=1,nx-1 
        do 20 m=1,nq 
        read(1,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m),                      &
     &                  f5(n,m),f7(n,m),f6(n,m),f8(n,m)                 
! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea         
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
   40 continue 
                                                                        
      do n=1,nx 
      xxl(n)=dlog(xx(n)) 
      enddo 
      do m=1,nq 
      qql(m)=dlog(qq(m)) 
      enddo 
                                                                        
      call jeppe1(nm,nx,nq,xxl,qql,f1,cc1(0,1,1,1,1,nset)) 
      call jeppe1(nm,nx,nq,xxl,qql,f2,cc2(0,1,1,1,1,nset)) 
      call jeppe1(nm,nx,nq,xxl,qql,f3,cc3(0,1,1,1,1,nset)) 
      call jeppe1(nm,nx,nq,xxl,qql,f4,cc4(0,1,1,1,1,nset)) 
      call jeppe1(nm,nx,nq,xxl,qql,f6,cc6(0,1,1,1,1,nset)) 
      call jeppe1(nm,nx,nq,xxl,qql,f8,cc8(0,1,1,1,1,nset)) 
                                                                        
      emc2=2.045 
      emb2=18.5 
                                                                        
      do 44 m=1,nqc 
      qqlc(m)=qql(m+nqc0) 
      do 44 n=1,nx 
      fc(n,m)=f5(n,m+nqc0) 
   44 continue 
      qqlc(1)=dlog(emc2) 
      call jeppe1(nm,nx,nqc,xxl,qqlc,fc,ccc(0,1,1,1,1,nset)) 
                                                                        
      do 45 m=1,nqb 
      qqlb(m)=qql(m+nqb0) 
      do 45 n=1,nx 
      fb(n,m)=f7(n,m+nqb0) 
   45 continue 
      qqlb(1)=dlog(emb2) 
      call jeppe1(nm,nx,nqb,xxl,qqlb,fb,ccb(0,1,1,1,1,nset)) 
                                                                        
                                                                        
      enddo 
   50 format(8f10.5) 
      return 
!                                                                       
                                                                        
      entry MRSTalfa(nflav,alfas,Qalfa) 
        call getnset(iset) 
!        mem = member(iset)                                             
!        call setnmem(member(iset))                                     
        call alphamrs(nflav,alfas,Qalfa) 
      return 
!                                                                       
      entry MRSTinit(Eorder,Q2fit) 
      return 
!                                                                       
      entry MRSTpdf(mem) 
!      if(mem.eq.0) mem=ndef                                            
!      imem = mem                                                       
        call getnset(iset) 
!        member(iset)=mem                                               
        call setnmem(iset,mem) 
                                                                        
      return 
!                                                                       
      END                                           
                                                                        
      subroutine jeppe1(imem,nx,my,xx,yy,ff,cc) 
      implicit real*8(a-h,o-z) 
      parameter(nnx=49,mmy=37,nhess=30) 
      dimension xx(nx),yy(my),ff(nnx,mmy),                              &
     &ff1(nnx,mmy),ff2(nnx,mmy),                                        &
     &ff12(nnx,mmy),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16),       &
     &cl(16),cc(0:nhess,nx,my,4,4),iwt(16,16)                           
                                                                        
      data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,                         &
     &                  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,                &
     &                  -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0,             &
     &                  2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0,               &
     &                  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,                &
     &                  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,                &
     &                  0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1,             &
     &                  0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1,               &
     &                  -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,             &
     &                  0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,             &
     &                  9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2,          &
     &                  -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2,      &
     &                  2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,               &
     &                  0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,               &
     &                  -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1,      &
     &                  4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/          
                                                                        
                                                                        
      do 42 m=1,my 
      dx=xx(2)-xx(1) 
      ff1(1,m)=(ff(2,m)-ff(1,m))/dx 
      dx=xx(nx)-xx(nx-1) 
      ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx 
      do 41 n=2,nx-1 
      ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m),        &
     &ff(n+1,m))                                                        
   41 continue 
   42 continue 
                                                                        
      do 44 n=1,nx 
      dy=yy(2)-yy(1) 
      ff2(n,1)=(ff(n,2)-ff(n,1))/dy 
      dy=yy(my)-yy(my-1) 
      ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy 
      do 43 m=2,my-1 
      ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m),        &
     &ff(n,m+1))                                                        
   43 continue 
   44 continue 
                                                                        
      do 46 m=1,my 
      dx=xx(2)-xx(1) 
      ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx 
      dx=xx(nx)-xx(nx-1) 
      ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx 
      do 45 n=2,nx-1 
      ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m),     &
     &ff2(n+1,m))                                                       
   45 continue 
   46 continue 
                                                                        
      do 53 n=1,nx-1 
      do 52 m=1,my-1 
      d1=xx(n+1)-xx(n) 
      d2=yy(m+1)-yy(m) 
      d1d2=d1*d2 
                                                                        
      yy0(1)=ff(n,m) 
      yy0(2)=ff(n+1,m) 
      yy0(3)=ff(n+1,m+1) 
      yy0(4)=ff(n,m+1) 
                                                                        
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
                                                                        
      do 47 k=1,4 
      z(k)=yy0(k) 
      z(k+4)=yy1(k)*d1 
      z(k+8)=yy2(k)*d2 
      z(k+12)=yy12(k)*d1d2 
   47 continue 
                                                                        
      do 49 l=1,16 
      xxd=0. 
      do 48 k=1,16 
      xxd=xxd+iwt(k,l)*z(k) 
   48 continue 
      cl(l)=xxd 
   49 continue 
      l=0 
      do 51 k=1,4 
      do 50 j=1,4 
      l=l+1 
      cc(imem,n,m,k,j)=cl(l) 
   50 continue 
   51 continue 
   52 continue 
   53 continue 
      return 
      END                                           
                                                                        
      subroutine jeppe2(i,x,y,nx,my,xx,yy,cc,z) 
!--   G.W. 02/07/2007 Allow extrapolation to small x and large q.       
      implicit real*8(a-h,o-z) 
      parameter(nhess=30) 
      dimension xx(nx),yy(my),cc(0:nhess,nx,my,4,4) 
                                                                        
      n=locx(xx,nx,x) 
      m=locx(yy,my,y) 
                                                                        
      if (n.gt.0.and.n.lt.nx.and.m.gt.0.and.m.lt.my) then 
!--   Do usual interpolation.                                           
         t=(x-xx(n))/(xx(n+1)-xx(n)) 
         u=(y-yy(m))/(yy(m+1)-yy(m)) 
         z=0.d0 
         do l=4,1,-1 
            z=t*z+((cc(i,n,m,l,4)*u+cc(i,n,m,l,3))*u                    &
     &           +cc(i,n,m,l,2))*u+cc(i,n,m,l,1)                        
         enddo 
                                                                        
      else if (n.eq.0.and.m.gt.0.and.m.lt.my) then 
!--   Extrapolate to small x.                                           
         call jeppe3(i,xx(1),y,nx,my,xx,yy,cc,f0) 
         call jeppe3(i,xx(2),y,nx,my,xx,yy,cc,f1) 
         if (f0.gt.0.d0.and.f1.gt.0.d0) then 
            z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1))) 
         else 
            z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1)) 
         end if 
                                                                        
      else if (n.gt.0.and.m.eq.my) then 
!--   Extrapolate to large q.                                           
         call jeppe3(i,x,yy(my),nx,my,xx,yy,cc,f0) 
         call jeppe3(i,x,yy(my-1),nx,my,xx,yy,cc,f1) 
         if (f0.gt.0.d0.and.f1.gt.0.d0) then 
            z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*        &
     &           (y-yy(my)))                                            
         else 
            z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my)) 
         end if 
                                                                        
      else if (n.eq.0.and.m.eq.my) then 
!--   Extrapolate to small x AND large q.                               
         call jeppe3(i,xx(1),yy(my),nx,my,xx,yy,cc,f0) 
         call jeppe3(i,xx(1),yy(my-1),nx,my,xx,yy,cc,f1) 
         if (f0.gt.0.d0.and.f1.gt.0.d0) then 
            z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*       &
     &           (y-yy(my)))                                            
         else 
            z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my)) 
         end if 
         call jeppe3(i,xx(2),yy(my),nx,my,xx,yy,cc,f0) 
         call jeppe3(i,xx(2),yy(my-1),nx,my,xx,yy,cc,f1) 
         if (f0.gt.0.d0.and.f1.gt.0.d0) then 
            z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))*       &
     &           (y-yy(my)))                                            
         else 
            z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my)) 
         end if 
         if (z0.gt.0.d0.and.z1.gt.0.d0) then 
            z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1))) 
         else 
            z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1)) 
         end if 
                                                                        
      else 
!--   Set parton distribution to zero otherwise.                        
         z = 0.d0 
                                                                        
      end if 
                                                                        
      return 
      END                                           
                                                                        
!--   G.W. 02/07/2007 Copy of the original jeppe2,                      
!--   only used for extrapolation.                                      
      subroutine jeppe3(i,x,y,nx,my,xx,yy,cc,z) 
      implicit real*8(a-h,o-z) 
      parameter(nhess=30) 
      dimension xx(nx),yy(my),cc(0:nhess,nx,my,4,4) 
      n=locx(xx,nx,x) 
      m=locx(yy,my,y) 
      t=(x-xx(n))/(xx(n+1)-xx(n)) 
      u=(y-yy(m))/(yy(m+1)-yy(m)) 
      z=0.d0 
      do l=4,1,-1 
         z=t*z+((cc(i,n,m,l,4)*u+cc(i,n,m,l,3))*u                       &
     &        +cc(i,n,m,l,2))*u+cc(i,n,m,l,1)                           
      enddo 
      return 
      END                                           
                                                                        
      integer function locx(xx,nx,x) 
      implicit real*8(a-h,o-z) 
      dimension xx(nx) 
!$$$      if(x.le.xx(1)) then                                           
                          ! G.W. 02/07/2007                             
      if(x.eq.xx(1)) then 
      locx=1 
      return 
      endif 
!$$$      if(x.ge.xx(nx)) then                                          
                           ! G.W. 02/07/2007                            
      if(x.eq.xx(nx)) then 
      locx=nx-1 
      return 
      endif 
      ju=nx+1 
      jl=0 
    1 if((ju-jl).le.1) go to 2 
      jm=(ju+jl)/2 
      if(x.ge.xx(jm)) then 
      jl=jm 
      else 
      ju=jm 
      endif 
      go to 1 
    2 locx=jl 
      return 
      END                                           
                                                                        
      real*8 function  polderiv(x1,x2,x3,y1,y2,y3) 
      implicit real*8(a-h,o-z) 
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*                    &
     &(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))   
      return 
      END                                           
