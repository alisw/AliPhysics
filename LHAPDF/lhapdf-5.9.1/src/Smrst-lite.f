! -*- F90 -*-
                                                                        
      subroutine jeppe1(imem,nx,my,xx,yy,ff,cc) 
      implicit real*8(a-h,o-z) 
      parameter(nnx=49,mmy=37,nhess=0) 
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
      parameter(nhess=0) 
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
      parameter(nhess=0) 
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
    1 if((ju-jl).le.1) goto 2 
      jm=(ju+jl)/2 
      if(x.ge.xx(jm)) then 
      jl=jm 
      else 
      ju=jm 
      endif 
      goto 1 
    2 locx=jl 
      return 
      END                                           
                                                                        
      real*8 function  polderiv(x1,x2,x3,y1,y2,y3) 
      implicit real*8(a-h,o-z) 
      polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1*                    &
     &(y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))   
      return 
      END                                           
