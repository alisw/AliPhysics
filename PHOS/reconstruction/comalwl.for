	parameter (nalw=1,ialw=40000,ialws=9,mxck=1000)
	common /wlihal/nt,mm(ialw),d(3,ialw),x(3,ialw),
     ,  ns(ialw),nks(ialw),ks(ialws,ialw)
c        common /iter/niter,st(ialw),sm(ialw),sqm(ialw),wm(ialw),
c     ,  cfi(ialw),wcf(ialw),ndl(ialw),fndlsv(ialw),alfa(ialw)
	common /wall/wl(ialw)
	common /walt/wlt(ialw),swlt(ialw),swlt0(ialw),swltmx(ialw)
	common /wlcl/cf(ialw)
	common /walgl/nglw(nalw),ipwgl(nalw),ncls(nalw),ipwlcl(nalw)

c	data fndlsv /ialw*1000./,alfa /ialw*0.000/
c ped
c	parameter (nchp=32)
c	common /pdihis/stat,npih(nchp,ialw),nout(ialw),max(ialw),
c     ,  nmax(ialw),nsum(ialw)
        COMMON /PEDIHP/STh,IHTHR(ialw),NST(ialw),NSTGT(ialw)
        common /gamthr/gmthr(ialw)
c        DATA IHTHR /ialw*0/
	common /wlclus/nck(nalw),etck(nalw),
     ,  emxck(mxck,nalw),eck(mxck,nalw),
     ,  xmick(3,mxck,nalw),xmack(3,mxck,nalw),xavck(3,mxck,nalw)
