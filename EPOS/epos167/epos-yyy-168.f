c----------------------------------------------------------------------
      subroutine hnbaaa156(ip,iret)
c----------------------------------------------------------------------
c  microcanonical decay of cluster ip via loop over hnbmet
c----------------------------------------------------------------------
      include 'epos.inc'
      common/cxyzt/xptl(mxptl),yptl(mxptl),zptl(mxptl),tptl(mxptl)
     *,optl(mxptl),uptl(mxptl),sptl(mxptl),rptl(mxptl,3)
      parameter(maxp=500)
      common/confg/np,amass(maxp),ident(maxp),pcm(5,maxp),wtxlog,wtlog
      common/citer/iter,itermx
      double precision tpro,zpro,ttar,ztar,ttaus,detap,detat
      common /cttaus/  tpro,zpro,ttar,ztar,ttaus,detap,detat 
      integer jc(nflav,2)
      double precision p(5),c(5)
      parameter(maxit=50000)
      common/count/nacc,nrej,naccit(maxit),nptot,npit(maxit)
      dimension be(4),pe(5),pa(5)
      common/yradx/yrad(maxp),phirad(maxp)
      common/xxxspecsy/ndrop(-4:4,-4:4,-4:4)
      common/cdelzet/delzet,delsgr /cvocell/vocell
      call utpri('hnbaaa',ish,ishini,4)

      if(ish.ge.3)then
      write(ifch,140)sngl(ttaus)
  140 format(/' ----------------------------------'/
     *'    droplet decay  at tau =',f6.2/
     *' ----------------------------------')
      write(ifch,*)'droplet:'
      call alist('&',ip,ip)
      endif
      
      iret=0
      do j=1,5
      c(j)=pptl(j,ip)
      enddo
                                    
      call idquac(ip,nqi,nsi,nai,jc)
      keu=jc(1,1)-jc(1,2)
      ked=jc(2,1)-jc(2,2)
      kes=jc(3,1)-jc(3,2)
      kec=jc(4,1)-jc(4,2)
      keb=jc(5,1)-jc(5,2)
      ket=jc(6,1)-jc(6,2)
      !print*,'droplet uds=',keu,ked,kes,'   E=',pptl(5,ip)

      volu=4./3.*pi*radptl(ip)**3

      if(volu.le.0)call utstop('hnbaaa: volume = 0&')
      if(volu.lt.0.01)then
        call utmsg('hnbaaa')
        write(ifch,*)'*****  very small volume:',volu
        write(ifch,*)
     *  'id:',idptl(ip),' r:',radptl(ip),' m:',pptl(5,ip)
        call utmsgf
      endif
   
    !~~~~~redefine energy in case of radial flow       
      amin=utamnu(keu,ked,kes,kec,keb,ket,4)   !utamnu(...,4) and not utamnu(...,5) 
      aumin=amuseg                             !otherwise droplet from remnant decay 
      ipo=ip                                   !could be too light after flow
      if(ityptl(ip).eq.60)ipo=iorptl(ip)
      tecmor=pptl(5,ipo)
      if(iappl.eq.4.or.iorsdf.ne.3
     &.or.ityptl(ip).eq.40.or.ityptl(ip).eq.50)then !not for droplets from remnants 
        yrmax=0                            
      else        	
        yrmax=max(0.0,yradmx+yradmi*alog10(engy/200.))
        	if(maproj.eq.1.and.matarg.eq.1)then
          yrmax=max(0.0,yradpp+yradpi*alog10(engy/1800.))
        	  !aumin=amin
        	  !if(yrmax.gt.0.2)print*,'===',tecmor,aamin,yrmax
        	endif  
      endif        	
      fradflo=1.
      if(yrmax.ne.0.)
     &fradflo=1./((sinh(yrmax)*yrmax-cosh(yrmax)+1.)/(yrmax**2/2.)) 
      tecm=pptl(5,ip)
      tecmxx=tecm
      if(yrmax.gt.0..and.tecmor.gt.aumin
     &  .and.tecm*fradflo.gt.amin) then
        ! redefine energy to account for collective flow
        ! \int_0^yrmax f(y) d2y = E_new (effective mass)  
        ! \int_0^yrmax cosh(y) f(y) d2y = E_old
        tecm=tecm*fradflo
        if(tecm.lt.amin)stop'aaahnb: small mass. should not happen.   '
      else 
        yrmax=0.  
      endif
 
    !~~~~~redefine energy in case of long coll flow  
      if(iappl.eq.4.or.iorsdf.ne.3
     &.or.ityptl(ip).eq.40.or.ityptl(ip).eq.50)then !not for droplets from remnants 
        yco=0
      else
       if(ylongmx.lt.0.)then
        yco=delzet * 1.75
       else
        yco=ylongmx
       endif
      endif 
      tecmx=tecm
      if(yco.gt.0..and.tecmor.gt.aumin) then
        tecm=tecm/sinh(yco)*yco     
      else 
        yco=0.  
      endif
      !print*,'========= cluster energy: ',pptl(5,ip),tecmx,tecm

    !~~~~~~~~~redefine volume~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ccc      vocri=tecm/epscri(ioclude)   !??????????????????????????????
ccc      volu=max(vocri,vocell)       !????????????????????????????

    !~~~~~~~~~decay~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      call hnbini(iret)
      !if(iret.ne.0)write(ifch,*)'***** unsucessfull hnbini *****'
      if(iret.ne.0)goto1000
      if(ioinct.ge.1)goto1
        
      do iter=1,itermx
        naccit(iter)=0
        call hnbmet
      enddo
     
1     continue

      if(ioceau.eq.1.and.iappl.eq.1)call xhnbte(ip)

    !~~~~~~~~~~long coll flow -> particles~~~~~~~~~~~~~~~~
      if(yco.gt.0.) then
        errlim=0.0001
        tecm=tecmx
 611    energ=0.
        do i=1,np
          yrad(i)=(2*rangen()-1)*yco
          be(3)=sinh(yrad(i))
          be(4)=cosh(yrad(i))
          energ=energ+be(4)*pcm(4,i)-be(3)*pcm(3,i)
        enddo
        if(abs(energ-tecm).gt.0.1) goto 611
        	!print*,'===== energy after flow boosts',energ,'   soll: ',tecm
        do j=1,4
        	  pe(j)=0.
        	enddo
        do i=1,np 
          be(1)= 0
          be(2)= 0
          be(3)= sinh(yrad(i))
          be(4)= cosh(yrad(i))
          call utlob3(1,be(1),be(2),be(3),be(4),1e0
     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          do j=1,4
          pe(j)=pe(j)+pcm(j,i)
          enddo
        enddo
        	pe(5)=sqrt(pe(4)**2-pe(3)**2-pe(2)**2-pe(1)**2)
        	!write(6,'(a,5e11.3)')'flow boosts',pe
        do j=1,4
        	  pa(j)=0.
        	enddo
        do i=1,np 
          call utlob3(1,pe(1),pe(2),pe(3),pe(4),pe(5)
     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          do j=1,4
        	    pa(j)=pa(j)+pcm(j,i)
        	  enddo
        enddo
        	pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
        	!write(6,'(a,5e11.3)')' cms boost ',pa
        esoll=tecm
        scal=1.
        do ipass=1,200
          sum=0.
          do  j=1,np
            do k=1,3
              pcm(k,j)=scal*pcm(k,j)
            enddo
            pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2+pcm(3,j)**2
     *           +amass(j)**2)
            sum=sum+pcm(4,j)
          enddo          
          scal=esoll/sum
          !write(6,*)'ipass,scal,e,esoll:'
          !    $         ,ipass,scal,sum,esoll 
          if(abs(scal-1.).le.errlim) goto301
        enddo
 301    continue      
        do j=1,4
        	  pa(j)=0.
        	enddo
        do i=1,np 
          do j=1,4
        	    pa(j)=pa(j)+pcm(j,i)
        	  enddo
        enddo
        pa(5)=sqrt(pa(4)**2-pa(3)**2-pa(2)**2-pa(1)**2)
        !write(6,'(a,5e11.3)')' rescaling ',pa
      endif

    !~~~~~~~~~~radial flow -> particles~~~~~~~~~~~~~~~~~~
      if(yrmax.gt.0.) then
        fecc=0
        	aa=1
        	bb=0
        	cc=0
        	dd=1
        if(ityptl(ip).eq.60)then        	
        	  ipo=iorptl(ip)
        	  xx=uptl(ipo)   ! <x**2>
        	  yy=optl(ipo)   ! <y**2>
        	  xy=desptl(ipo) ! <x*y>
        	  dta=0.5*abs(xx-yy)
        	  ev1=(xx+yy)/2+sqrt(dta**2+xy**2)
        	  ev2=(xx+yy)/2-sqrt(dta**2+xy**2)
        	  if(xy.lt.0..and.dta.ne.0.)then
        	    theta=0.5*atan(-xy/dta)
        	  elseif(xy.gt.0..and.dta.ne.0.)then
        	    theta=-0.5*atan(xy/dta)
        	  else
        	    theta=0
        	  endif
        	  !eccx=(yy-xx)/(yy+xx)
          yy=ev1
        	  xx=ev2
             	  ecc=(yy-xx)/(yy+xx)
        	  !print*,eccx,ecc,theta
        	  fecc=facecc*ecc
        	  fecc=min(0.3,fecc)
        	  phiclu=phievt+theta
          aa=cos(phiclu)
          bb=sin(phiclu)
        	  cc=-sin(phiclu)
        	  dd=cos(phiclu)
        endif
        	errlim=0.0001
        tecm=tecmxx
 610    energ=0.
        do i=1,np
          yrad(i)=sqrt(rangen())
          phirad(i)=2.*pi*rangen()
          bex=-dsinh(dble(yrad(i)*yrmax))*cos(phirad(i))*(1+fecc) 
          bey=-dsinh(dble(yrad(i)*yrmax))*sin(phirad(i))*(1-fecc)
          be(1)=aa*bex+cc*bey
          be(2)=bb*bex+dd*bey
          be(3)=-0d0    
          be(4)=sqrt(1+be(1)**2+be(2)**2)
          bp=0d0
          do k=1,3
            bp=bp+pcm(k,i)*be(k) 
          enddo
          en=be(4)*pcm(4,i)+bp
          energ=energ+en
        enddo
        if(abs(energ-tecm).gt.0.1) goto 610
        energ=0.
        do i=1,np 
          bex=dsinh(dble(yrad(i)*yrmax))*cos(phirad(i))*(1+fecc)
          bey=dsinh(dble(yrad(i)*yrmax))*sin(phirad(i))*(1-fecc)
          be(1)=aa*bex+cc*bey
          be(2)=bb*bex+dd*bey
          be(3)=0d0           
          be(4)=sqrt(1+be(1)**2+be(2)**2)
          call utlob3(1,be(1),be(2),be(3),be(4),1e0
     *         , pcm(1,i), pcm(2,i), pcm(3,i), pcm(4,i))
          energ=energ+pcm(4,i)
        enddo
        esoll=tecm
        scal=1.
        do ipass=1,200
          sum=0.
          do  j=1,np
            do k=1,3
              pcm(k,j)=scal*pcm(k,j)
            enddo
            pcm(4,j)=sqrt(pcm(1,j)**2+pcm(2,j)**2+pcm(3,j)**2
     *           +amass(j)**2)
            sum=sum+pcm(4,j)
          enddo          
          scal=esoll/sum
          !write(6,*)'ipass,scal,e,esoll:'
          ! $         ,ipass,scal,sum,esoll 
          if(abs(scal-1.).le.errlim) goto300
        enddo
 300    continue      
      else
        do n=1,np
          yrad(n)=0.
          phirad(n)=0.
        enddo
      endif
    !~~~~~~~~~~~~~~~

      nptlb=nptl
      do n=1,np
        nptl=nptl+1
        if(nptl.gt.mxptl)call utstop('hnbptl: mxptl too small&')
        idptl(nptl)=ident(n)
        do j=1,4
          p(j)=pcm(j,n)
        enddo
        p(5)=amass(n)
        call utlob2(-1,c(1),c(2),c(3),c(4),c(5),p(1),p(2),p(3),p(4),10)
        do j=1,5
          pptl(j,nptl)=p(j)
        enddo
        if(tecmor.gt.aumin)then
          ityptl(nptl)=60
        else
          ityptl(nptl)=19
        endif
        if(ityptl(ip).eq.60)then
          if(ityptl(nptl).eq.60)then
            ipo=iorptl(ip)
            xx=uptl(ipo)        ! <x**2>
            yy=optl(ipo)        ! <y**2>
            rini=sqrt(5./3.*(xx+yy)) !<r**2>=3/5*R**2 for sphere of radius R
            r=1.15*rini*yrad(n) !yrad=y/ymax
            tau=2.25/sqrt(yrad(n)**2+0.04)-0.75
            z=xorptl(3,ipo)
            t=xorptl(4,ipo)
           !zeta=0.5*log((t+z)/(t-z))-0.5*delzet+2*0.5*delzet*rangen()
            zeta=0.5*log((p(4)+p(3))/(p(4)-p(3)))
            z=tau*sinh(zeta)
            t=tau*cosh(zeta)
            xorptl(1,nptl)=xorptl(1,ipo)+r*cos(phirad(n))
            xorptl(2,nptl)=xorptl(2,ipo)+r*sin(phirad(n))
            xorptl(3,nptl)=z
            xorptl(4,nptl)=t
          else
            xorptl(1,nptl)=xorptl(1,ip)
            xorptl(2,nptl)=xorptl(2,ip)
            xorptl(3,nptl)=xorptl(3,ip)
            xorptl(4,nptl)=xorptl(4,ip)
          endif 
        endif
      enddo 

      if(ish.ge.3)then
        write(ifch,*)'decay products:'
        call alist('&',nptlb+1,nptl)
        if(ish.ge.5)then
          write(ifch,*)'momentum sum:'
          do kk=1,5
            pptl(kk,nptl+1)=0
            do ii=nptlb+1,nptl
              pptl(kk,nptl+1)=pptl(kk,nptl+1)+pptl(kk,ii)
            enddo
            pptl(kk,nptl+2)=c(kk)
          enddo
          call alist('&',nptl+1,nptl+2)
        endif
      endif

 1000 continue
      call utprix('hnbaaa',ish,ishini,4)
      return
      end
