	parameter (Npg=1000,nssp=10)
        COMMON /GAMMAx/KGx,MWx(npg),IDx(npg),JDx(npg),Ex(npg),E4x(npg),
     ,  XWx(npg),YWx(npg),ESx(nssp,npg),ETx(nssp,npg),ISsdx(npg),
     ,  IGDEVx(npg),ZGDEVx(npg)
	COMMON /GAMx/NGx,EGx(npg),CGx(3,npg),PGx(4,npg),hgx(3,npg),    
     ,ihgw(npg),label(npg),egxr0(npg),pgxr0(3,npg),hgxr0(3,npg)
	COMMON /trax/Ntrx,ichrg(npg),Ctrx(3,npg),Ptrx(4,npg),
     ,htrx(3,npg),ihtrw(npg)
