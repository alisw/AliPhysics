// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 * 
 **************************************************************************/
#include "AliEveTOFStrip.h"

#include <TEveManager.h>

#include <AliTOFdigit.h>
#include <AliTOFGeometry.h>

#include <TStyle.h>

using namespace std;

Bool_t       AliEveTOFStrip::fgStaticInitDone = kFALSE;
TEveFrameBox*    AliEveTOFStrip::fgTOFstripFrameBox = 0;
TEveRGBAPalette* AliEveTOFStrip::fgTOFstripPalette  = 0;

//_______________________________________________________
ClassImp(AliEveTOFStrip)

/* ************************************************************************ */

AliEveTOFStrip::AliEveTOFStrip(const Text_t* n, const Text_t* t) :
  TEveQuadSet(n, t),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(0),
  fSector(-1), fPlate(-1), fStrip(-1),
  fDx(0), fDz(0)
{

  fGeoManager = (TGeoManager*)gEve->GetGeometry("$REVESYS/alice-data/alice_fullgeo.root");
  if (!fGeoManager) printf("ERROR: no TGeo\n");

}
/* ************************************************************************ */

AliEveTOFStrip::AliEveTOFStrip(TGeoManager *localGeoManager,
		   Int_t nSector, Int_t nPlate, Int_t nStrip)
  :
  TEveQuadSet(Form("Strip%i",nStrip)),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(0),
  fSector(nSector), fPlate(nPlate), fStrip(nStrip),
  fDx(0), fDz(0),
  fGeoManager(localGeoManager)
{
  
  //if (!fGeoManager) printf("ERROR: no TGeo\n");

  InitModule();

}
/* ************************************************************************ */

AliEveTOFStrip::AliEveTOFStrip(TGeoManager *localGeoManager,
		   Int_t nSector, Int_t nPlate, Int_t nStrip,
		   TClonesArray *tofArray)
  :
  TEveQuadSet(Form("Strip%i",nStrip)),
  fTOFgeometry(new AliTOFGeometry()),
  fTOFarray(tofArray),
  fSector(nSector), fPlate(nPlate), fStrip(nStrip),
  fDx(0), fDz(0),
  fGeoManager(localGeoManager)
{

  InitModule();

}
/* ************************************************************************ */

AliEveTOFStrip::~AliEveTOFStrip()
{

  fGeoManager = 0x0;
  delete fGeoManager;

  fTOFarray = 0x0;
  delete fTOFarray;

}

/* ************************************************************************ */
/*
void AliEveTOFStrip::SetDigitsInfo(AliEveTOFDigitsInfo* info)
{
  if(fInfo) fInfo->DecRefCount();
  fInfo = info;
  if(fInfo) fInfo->IncRefCount();

}
*/
/* ************************************************************************ */
void AliEveTOFStrip::InitStatics()
{
  if (fgStaticInitDone) return;

  Float_t dx = 2.5*48;
  Float_t dz = 3.5*2;
  fgTOFstripFrameBox = new TEveFrameBox();

  fgTOFstripFrameBox->SetAAQuadXZ(-dx*0.5, 0, -dz*0.5, dx, dz);
  fgTOFstripFrameBox->SetFrameColor((Color_t) 32);//31);

  //fgTOFstripPalette  = new TEveRGBAPalette(0, 2048); // TOT
  fgTOFstripPalette  = new TEveRGBAPalette(0, 8192); // TDC

  fgStaticInitDone = kTRUE;
}

/* ************************************************************************ */
void AliEveTOFStrip::InitModule()
{

  fDx = fTOFgeometry->XPad()*fTOFgeometry->NpadX();
  fDz = fTOFgeometry->ZPad()*fTOFgeometry->NpadZ();

  if (!fgStaticInitDone) InitStatics();

  SetFrame(fgTOFstripFrameBox);
  SetPalette(fgTOFstripPalette);
  //fFrame   = fgTOFstripFrameBox;
  //fPalette = fgTOFstripPalette;

  LoadQuads();  
  ComputeBBox();
  SetTrans();

}

/* ************************************************************************ */
void AliEveTOFStrip::LoadQuads()
{

  //Int_t n_col = gStyle->GetNumberOfColors();

  Int_t iPadX = -1;
  Int_t iPadZ = -1;
  Int_t tdc = -1;
  Int_t tot = -1;
  Float_t x = -1;
  Float_t z = -1;

  Reset(kQT_RectangleXZFixedY, kFALSE, 32);

  AliTOFdigit *tofDigit;

  //printf(" fTOFarray->GetEntries() = %4i \n",fTOFarray->GetEntries());

  for (Int_t ii=0; ii<fTOFarray->GetEntries(); ii++) {

    tofDigit = (AliTOFdigit*)fTOFarray->UncheckedAt(ii);

    iPadX = tofDigit->GetPadx();
    iPadZ = tofDigit->GetPadz();

    tot = tofDigit->GetToT();
    tdc = tofDigit->GetTdc();

    //if (fSector==4 && fPlate==2 && fStrip==0) printf(" %2i  %1i\n", iPadX, iPadZ);
    //if (iPadX==23 || iPadX==24) printf(" %2i  %1i %2i \n", fSector, fPlate, fStrip);

    fTOFgeometry->DetToStripRF(iPadX, iPadZ, x, z);

    AddQuad(x, z, 2.5, 3.5);
    //AddQuad(-2.5*0.5, -3.5*0.5, 2.5, 3.5);

    // In principle could have color based on number of neigbours. We
    // can insert the time-of-flight value for each pad
    //QuadValue((Int_t)tot);
    QuadValue((Int_t)tdc);

    //if (fSector==4 && fPlate==2  && fStrip==0) printf("  %1i   %2i    %f  %f \n", iPadZ, iPadX, x, z);

  }

  RefitPlex();
 
}

/* ************************************************************ */
void AliEveTOFStrip::SetTrans()
{

  fHMTrans.UnitTrans();

  //Int_t det[5] = {fSector, fPlate, fStrip, -1, -1};
  Char_t path[100];
  //fTOFgeometry->GetVolumePath(det,path);
  fTOFgeometry->GetVolumePath(fSector, fPlate, fStrip, path);

  fGeoManager->cd(path);
  TGeoHMatrix global = *fGeoManager->GetCurrentMatrix();
  Double_t *rotMat = global.GetRotationMatrix();
  
  /*
  // ok till 19 April 2007
  fHMTrans.SetBaseVec(1, rotMat[0], rotMat[1], rotMat[2]);
  fHMTrans.SetBaseVec(2, rotMat[3], rotMat[4], rotMat[5]);
  fHMTrans.SetBaseVec(3, rotMat[6], rotMat[7], rotMat[8]);
  */

  fHMTrans.SetBaseVec(1, rotMat[0], rotMat[3], rotMat[6]);
  fHMTrans.SetBaseVec(2, rotMat[1], rotMat[4], rotMat[7]);
  fHMTrans.SetBaseVec(3, rotMat[2], rotMat[5], rotMat[8]);

  Double_t *tr = global.GetTranslation();
  fHMTrans.SetBaseVec(4, tr);

}
