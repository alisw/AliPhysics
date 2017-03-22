c -*-Fortran-*-

      integer maxpup
      parameter(maxpup=100)
      integer idbmup,pdfgup,pdfsup,idwtup,nprup,lprup
      double precision ebmup,xsecup,xerrup,xmaxup
      common /heprup/ idbmup(2),ebmup(2),pdfgup(2),pdfsup(2),
     &                idwtup,nprup,xsecup(maxpup),xerrup(maxpup),
     &                xmaxup(maxpup),lprup(maxpup)
      integer maxnup
      parameter (maxnup=500)
      integer nup,idprup,idup,istup,mothup,icolup
      double precision xwgtup,scalup,aqedup,aqcdup,pup,vtimup,spinup
      common/hepeup/nup,idprup,xwgtup,scalup,aqedup,aqcdup,
     &              idup(maxnup),istup(maxnup),mothup(2,maxnup),
     &              icolup(2,maxnup),pup(5,maxnup),vtimup(maxnup),
     &              spinup(maxnup)
      save /hepeup/
