	parameter (ncelmax=10)
	parameter (narrmax=10,nholmax=10)
	parameter (ncompmax=2,ncomparmax=5)
	parameter (mmmax=40000)

	common /comcells/ncells,idcells(ncelmax),cellsize(3,ncelmax),
     ,  idcelmat(ncelmax)

	common /comarray/narray,ijarray(2,narrmax),idarray(narrmax),
     ,  idcellar(narrmax),isizarray(narrmax),
     ,  sizarray(3,narrmax),
     ,  nholes(narrmax),poshole(3,nholmax,narrmax),
     ,  sizhole(3,nholmax,narrmax)

	common /comcompose/ncompose,poscomp(3,ncompmax),
     ,  ncomparr(ncompmax),
     ,  poscompar(3,ncomparmax,ncompmax),idcomparr(ncomparmax,ncompmax),
     ,  mmaxcompar(ncomparmax,ncompmax),mmincompar(ncomparmax,ncompmax),
     ,  mmincomp(ncompmax),mmaxcomp(ncompmax)

	common /comspace/ipmmm,mmcells(mmmax),cellpos(3,mmmax),
     ,  neibortyp(mmmax)
	
	
	
