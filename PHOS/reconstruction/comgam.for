	parameter (NGp=1000,nsps=10,nvertmax=1000)
        COMMON /RCGAMMA/KG,MW(ngp),ID(ngp),JD(ngp),E(ngp),E4(ngp),
     ,  XW(ngp),YW(ngp),ES(nsps,ngp),ET(nsps,ngp),ISsd(ngp),
     ,  IGDEV(ngp),ZGDEV(ngp),sigexy(3,ngp),Emimx(2,nsps,ngp),
     ,  kgfix,igfix(ngp),cgfix(3,ngp),sgfix(3,ngp),hiw(ngp),
     ,  wsw(nsps,ngp),h1w(ngp),h0w(ngp),raxay(5,ngp),
     ,  sigmaes0(nsps,ngp),dispeces(nsps,ngp),
     ,  igamvert(ngp)
c      COMMON/MASS/ EFMAS(ngp,ngp),ang(ngp,ngp),MES(ngp,ngp),
c     ,himes(ngp,ngp,3)
      COMMON/P4/P4(4,ngp),PTOT(4),EFMTOT,sigp4(3,ngp)
      COMMON/P4f/P4f(4,ngp),PTOTf(4),EFMfit
      COMMON/comvertx/ntotvertx,gamvertex(3,nvertmax)
