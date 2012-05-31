//
// This macro is useful to check the correct corrispondence
// between the TOF volume GEANT numbering
// (i.e. sector [0;17], module [0;4], strip [0;14/18], padZ [0;1], padX [0;47])
// and the equipment identification
// (i.e. DDL [0;71], TRM [3/4;12], TDC [0;14], chain [0;1]; channel [0;7]).
//
// Author: A. De Caro (decaro@sa.infn.it - annalisa.de.caro@cern.ch)
//

#if !defined( __CINT__) || defined(__MAKECINT__)
#include <Riostream.h>

#include "AliTOFRawStream.h"
#endif

void CheckEQID2GEANT2EQID(Int_t drmCheck=0, Int_t trmCheck=3);
void CheckGEANT2EQID2GEANT(Int_t selectedSector=0);

void CheckEQID2GEANT2EQID(Int_t drmCheck, Int_t trmCheck)
{

  Int_t iDDL    = -1;
  Int_t iSector = -1;
  Int_t iPlate  = -1;
  Int_t iStrip  = -1;
  Int_t iPadZ   = -1;
  Int_t iPadX   = -1;

  Int_t vol[5] = {-1, -1, -1, -1, -1};

  Int_t sDDL = -1;
  Int_t sTRM = -1;
  Int_t sTDC = -1;
  Int_t sChain = -1;
  Int_t sChannel = -1;

  Int_t ddlInf = drmCheck;
  Int_t ddlSup = drmCheck+1;
  Int_t trmInf = trmCheck;
  Int_t trmSup = trmCheck+1;

  if (drmCheck==-1) ddlInf = 0, ddlSup = 72;
  if (trmCheck==-1) trmInf = 3, trmSup = 12;
  for (Int_t nDDL=ddlInf; nDDL<ddlSup; nDDL++) {
    iSector = AliTOFRawStream::GetSectorNumber(nDDL);

    iDDL = AliTOFRawStream::GetDDLnumberPerSector(nDDL);

    for (Int_t nTRM=trmInf; nTRM<=trmSup; nTRM++) {
      for (Int_t nTDC=0 ; nTDC<15; nTDC++) {
      if (
	  (nDDL%2==1 && nTRM==3 && nTDC>2)
	  ||
	  (nDDL%2==0 && nTRM==3)
	  )
	continue;

	iPlate = AliTOFRawStream::Equip2VolNplate(iDDL, nTRM, nTDC);
	iStrip = AliTOFRawStream::Equip2VolNstrip(iDDL, nTRM, nTDC);

	for (Int_t nChain=0; nChain<2; nChain++) {
	  for (Int_t nChannel=0 ; nChannel<8; nChannel++) {

	    //iPad = AliTOFRawStream::Equip2VolNpad(iDDL, nChain, nTDC, nChannel);
	    iPadX = AliTOFRawStream::Equip2VolNpadX(iDDL, nChain, nTDC, nChannel);
	    iPadZ = AliTOFRawStream::Equip2VolNpadZ(iDDL, nChain, nTDC, nChannel);

	    vol[0] = iSector, vol[1] = iPlate, vol[2] = iStrip, vol[3] = iPadZ, vol[4] = iPadX;

	    sDDL = AliTOFRawStream::Geant2DDL(vol);
	    sTRM = AliTOFRawStream::Geant2TRM(vol);
	    sTDC = AliTOFRawStream::Geant2TDC(vol);
	    sChain = AliTOFRawStream::Geant2Chain(vol);
	    sChannel = AliTOFRawStream::Geant2Channel(vol);

	    if (
		sDDL!=nDDL ||
		sTRM!=nTRM ||
		sTDC!=nTDC ||
		sChain!=nChain ||
		sChannel!=nChannel
		)
	      {
		printf(" %2i, %2i, %2i, %1i, %1i ---> ",
		       nDDL,
		       nTRM,
		       nTDC,
		       nChain,
		       nChannel
		       );
		printf(" %2i %1i %2i %1i %2i ---> ",
		       iSector, iPlate, iStrip, iPadZ, iPadX);
		printf(" %2i, %2i, %2i, %1i, %1i\n",
		       sDDL,
		       sTRM,
		       sTDC,
		       sChain,
		       sChannel
		       );
	      }

	  }
	}
      }
    }

  }


}



void CheckGEANT2EQID2GEANT(Int_t selectedSector)
{

  Int_t iSector = -1;
  Int_t iPlate  = -1;
  Int_t iStrip  = -1;
  Int_t iPadZ   = -1;
  Int_t iPadX   = -1;

  Int_t vol[5] = {-1, -1, -1, -1, -1};

  Int_t nDDL    = -1;
  Int_t iDDL = -1;
  Int_t iTRM = -1;
  Int_t iTDC = -1;
  Int_t iChain = -1;
  Int_t iChannel = -1;

  Int_t infSector= selectedSector;
  Int_t supSector= selectedSector+1;

  if (selectedSector==-1) infSector = 0, supSector = 18;

  for (Int_t nSector=infSector; nSector<supSector; nSector++) {
    for (Int_t nPlate=0; nPlate<5; nPlate++) {
      if ((nSector==13 || nSector==14 || nSector==15) && nPlate==2) continue;
      for (Int_t nStrip=0; nStrip<AliTOFGeometry::NStrip(nPlate); nStrip++) {
	for (Int_t nPadZ=0; nPadZ<2; nPadZ++) {
	  for (Int_t nPadX=0; nPadX<48; nPadX++) {

	    vol[0] = nSector, vol[1] = nPlate, vol[2] = nStrip, vol[3] = nPadZ, vol[4] = nPadX;

	    nDDL = AliTOFRawStream::Geant2DDL(vol);
	    iDDL = AliTOFRawStream::GetDDLnumberPerSector(nDDL);

	    iTRM = AliTOFRawStream::Geant2TRM(vol);
	    iTDC = AliTOFRawStream::Geant2TDC(vol);
	    iChain = AliTOFRawStream::Geant2Chain(vol);
	    iChannel = AliTOFRawStream::Geant2Channel(vol);

	    iSector = AliTOFRawStream::GetSectorNumber(nDDL);
	    iPlate = AliTOFRawStream::Equip2VolNplate(iDDL, iTRM, iTDC);
	    iStrip = AliTOFRawStream::Equip2VolNstrip(iDDL, iTRM, iTDC);
	    iPadX = AliTOFRawStream::Equip2VolNpadX(iDDL, iChain, iTDC, iChannel);
	    iPadZ = AliTOFRawStream::Equip2VolNpadZ(iDDL, iChain, iTDC, iChannel);

	    if (
		nSector!= iSector ||
		nPlate!= iPlate  ||
		nStrip!= iStrip ||
		nPadZ!= iPadZ  ||
		nPadX!= iPadX
		)
	      {
		printf(" %2i %1i %2i %1i %2i ---> ",
		       iSector, iPlate, iStrip, iPadZ, iPadX);
		printf(" %2i, %2i, %2i, %1i, %1i ---> ",
		       iDDL,
		       iTRM,
		       iTDC,
		       iChain,
		       iChannel
		       );
		printf(" %2i %1i %2i %1i %2i\n    ",
		       nSector, nPlate, nStrip, nPadZ, nPadX);
	      }

	  }
	}
      }
    }
  }


}
