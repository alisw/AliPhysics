      subroutine USERevolve(x,Q,f)
      implicit none
      integer nset
      real*8 f(-6:6)
      real*8 x,Q
      real*8 alfas
      real*8 Eorder,Q2fit	
      return
*
      entry USERread(nset)
      return
*
      entry USERalfa(alfas,Q)
      return
*
      entry USERinit(nset,Eorder,Q2fit)
      return
*
      entry USERpdf(nset)
      return
*
      end
*
