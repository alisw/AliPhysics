*******************************************************************
      SUBROUTINE eloss_qsimp(sr,sd)
      REAL*8           xr,length,noc,qq,nnorm,rri,fracc
      COMMON /input/   xr,length,noc,qq,nnorm,rri,fracc
      INTEGER JMAX
      REAL*8 a1,b1,diff1,d1,sr,sd,EPS,fu1r,fu1d
      EXTERNAL eloss_func1, eloss_func2
      PARAMETER (EPS=1.e-9, JMAX=12)
      REAL*8 osr,ostr
      ostr=-1.e10
      osr=-1.e10
      a1=0.0
      b1=(0.99999-xr)*fracc
      diff1 = b1-a1
      d1 = 0.5*diff1
*
      res=dgauss(eloss_func2,a1,b1,1.d-6)
      call eloss_func1(a1,fu1r,fu1d)
      sd=fu1d
      sr=res
      END

      function eloss_func2(yy)
      implicit double precision (a-h,o-z)
      call eloss_func1(yy,fu1r,fu1d)
      eloss_func2=fu1r
      return
      end


**************************************************************
*
      SUBROUTINE eloss_func1(yy,funr,fund)
*
      REAL*8           funr,yy,fund
      REAL*8           xr,length,noc,qq,nnorm,rri,fracc
      COMMON /input/   xr,length,noc,qq,nnorm,rri,fracc
      EXTERNAL         eloss_lookup
      REAL*8           cont, disc, wwt, tepsi
*
      tepsi = yy
      wwt = tepsi
      if(wwt.ge.1.3) then
         call eloss_lookup(rri,1.d0,cont,disc)
         funr=0.0
         fund=disc
      else
         call eloss_lookup(rri,wwt,cont,disc)
         funr = cont*eloss_fragm(xr/(1.0-tepsi/fracc),qq)
     .         /(1.0-tepsi/fracc)
         fund = disc
      endif
      END
*******************************************************************
      SUBROUTINE eloss_lookup(rrrr,xxxx,continuous,discrete)
*
      REAL*8           xx(400), da(30), ca(30,260), rrr(30)
      COMMON /data/    xx, da, ca, rrr
      REAL*8           rrrr,xxxx, continuous, discrete
      REAL*8           rrin, xxin
      INTEGER          nrlow, nrhigh, nxlow, nxhigh
      REAL*8           rrhigh, rrlow, rfraclow, rfrachigh
      REAL*8           xfraclow, xfrachigh
      REAL*8           clow, chigh
*
      rrin = rrrr
      xxin = xxxx
*
*    determine the tabulated values xx(nxlow), xx(nxhigh)
*    rrlow, rrhigh such that
*    xx(nxlow) < xxin <  xx(nxhigh)
*    rrlow < rrin < rrhigh
*
      nxlow = int(xxin/0.005) + 1
      nxhigh = nxlow + 1
      xfraclow = (xx(nxhigh)-xxin)/0.005
      xfrachigh = (xxin - xx(nxlow))/0.005
*
      do 666, nr=1,30
         if (rrin.lt.rrr(nr)) then
            rrhigh = rrr(nr)
         else
            rrhigh = rrr(nr-1)
            rrlow = rrr(nr)
            nrlow = nr
            nrhigh = nr-1
            goto 665
         endif
 666     enddo
 665     continue
*
      rfraclow = (rrhigh-rrin)/(rrhigh-rrlow)
      rfrachigh = (rrin-rrlow)/(rrhigh-rrlow)
*
      clow = xfraclow*ca(nrlow,nxlow)+xfrachigh*ca(nrlow,nxhigh)
      chigh = xfraclow*ca(nrhigh,nxlow)+xfrachigh*ca(nrhigh,nxhigh)
      continuous = rfraclow*clow + rfrachigh*chigh
      discrete = rfraclow*da(nrlow) + rfrachigh*da(nrhigh)
*
      END

***************************************************************

c	BKK FF

      FUNCTION eloss_fragmbkk(xxx,qqq)
      REAL*8   alphav, betav, gammav, nv
      REAL*8   alphas, betas, gammas, ns
      REAL*8   sbar, xx, qq, xxx, qqq, lambda, fragv, frags
*
      xx = xxx
      qq = qqq
      lambda = 1.0
      sbar=log(log(qq*qq/(lambda*lambda))/log(4.0/(lambda*lambda)))
*
      alphav = -1.0 - 0.0272*sbar
      betav = 1.2 + 0.67*sbar
      gammav = -0.393*sbar
      nv = 0.551 - 0.053*sbar - 0.032*sbar*sbar
*
*      alphav = -1.0 - 0.059*sbar
*      betav = 1.2 + 0.6*sbar
*      gammav = -0.163*sbar
*      nv = 0.338 - 0.064*sbar - 0.0105*sbar*sbar
*
      alphas = -1.0 + 0.447*sbar - 0.266*sbar*sbar
      betas = 4.7 - 2.88*sbar + 2.05*sbar*sbar
      gammas = -9.01*sbar + 4.36*sbar*sbar
      ns = 1.23 + 2.85*sbar - 1.6*sbar*sbar
*
*      alphas = -1.0 + 0.757*sbar - 0.537*sbar*sbar
*      betas = 5.26 - 5.22*sbar + 3.62*sbar*sbar
*      gammas = -13.6*sbar + 8.17*sbar*sbar
*      ns = 1.19 + 4.20*sbar - 2.86*sbar*sbar
*
      fragv = nv*(xx**alphav)*((1.0-xx)**betav)*((1.0+xx)**gammav)
      frags = ns*(xx**alphas)*((1.0-xx)**betas)*((1.0+xx)**gammas)
c      fragmbkk = (fragv+frags)
      eloss_fragmbkk = fragv
      END

***************************************************************

      FUNCTION eloss_fragm(xxx,qqq)
      REAL*8   alpha, beta, gamma, n
      REAL*8   sbar, xx, qq, xxx, qqq, lambda
*
      xx = xxx
      qq = qqq
      lambda = 0.088
      sbar=log(log(qq*qq/(lambda*lambda))/log(2.0/(lambda*lambda)))
*

c	u or d -> pi+ + pi-

      n=0.54610-0.22946*sbar-0.22594*sbar**2+0.21119*sbar**3
      alpha=-1.46616-0.45404*sbar-0.12684*sbar**2+0.27646*sbar**3
      beta=1.01864+0.95367*sbar-1.09835*sbar**2+0.74657*sbar**3
      gamma=-0.01877*sbar+0.02949*sbar**2

c	g -> pi+ + pi-

      n=6.04510-6.61523*sbar-1.64978*sbar**2+2.68223*sbar**3
      alpha=-.71378+0.14705*sbar-1.08423*sbar**2-.43182*sbar**3
      beta=2.92133+1.48429*sbar+1.32887*sbar**2-1.78696*sbar**3
      gamma=0.23086*sbar-0.29182*sbar**2

      eloss_fragm = n*xx**alpha*(1.-xx)**beta*(1.+gamma/xx)/2.

      END

