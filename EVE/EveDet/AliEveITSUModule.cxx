
/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include <AliEveITSUModule.h>
#include <TGeoMatrix.h>
#include <TEveTrans.h>
#include <TClonesArray.h>
#include <TStyle.h>
#include "AliLog.h"

#include <TEveManager.h>
#include <TGeoManager.h>
#include <AliEveEventManager.h>
#include <AliGeomManager.h>

#include <AliITSUGeomTGeo.h>
#include <AliITSMFTSegmentationPix.h>
#include <AliITSMFTDigitPix.h>
//______________________________________________________________________________
//
// Visualization of an ITS Upgrade module.

ClassImp(AliEveITSUModule)

Bool_t AliEveITSUModule::fgStaticInitDone = 0;

TEveFrameBox*    AliEveITSUModule::fgITSUFrameBox     = 0;
TEveFrameBox*    AliEveITSUModule::fgITSUFrameBoxDead = 0;
TEveRGBAPalette* AliEveITSUModule::fgITSUPalette  = 0;

AliITSUGeomTGeo* fGM                 = 0;
const AliITSMFTSegmentationPix* fSegm      = 0;

/******************************************************************************/

AliEveITSUModule::AliEveITSUModule(const Text_t* n, const Text_t* t) :
  TEveQuadSet(n, t),
  fID(0),
  fkLayer(0),
  fkLadder(0),
  fkDetector(0),
  fDpx(0), fDpz(0),
  fAtLeastOneDigit(kFALSE)
{
  // Constructor.
  
}

AliEveITSUModule::AliEveITSUModule(AliITSUGeomTGeo *gm, Int_t id, Int_t layer, Int_t ladder, Int_t detector) :
  TEveQuadSet(Form("ITSU module %d; (lay,lad,det)=(%d,%d,%d)", id,layer,ladder,detector),Form("%d",id)),
  fID(id),
  fkLayer(layer),
  fkLadder(ladder),
  fkDetector(detector),
  fDpx(0), fDpz(0),
  fAtLeastOneDigit(kFALSE)
{

  // 
  // constructor
  //
  fGM = gm; // ITSU Geometry Manager
  fgStaticInitDone = kFALSE; 
  fSegm = fGM->GetSegmentation(layer);
  fDpx = fSegm->Dpx(0);  // pixel pitch in x
  fDpz = fSegm->Dpz(0);  // pixel pitch in z
  SetID(id);
  //
}

AliEveITSUModule::~AliEveITSUModule()
{
  // Destructor.

 
}

/******************************************************************************/

void AliEveITSUModule::InitStatics()
{
  // Initialize static variables.
  //
  // Warning all sensor sizes are cm
  // In Eve half-lengths/widths are used, hence 1/2.

  if (fgStaticInitDone) return;
  fgStaticInitDone = kTRUE;

  Float_t dx =  fSegm->Dx(); // dimension in x in cm
  Float_t dz =  fSegm->Dz(); // dimension in y in cm
  Float_t dy =  0;// ? eventuelly a few 100 micron, right?

  {
    fgITSUFrameBox = new TEveFrameBox();
    fgITSUFrameBox->SetAAQuadXZ(-dx/2, dy, -dz/2, dx, dz);
    fgITSUFrameBox->SetFrameColor(kBlue-4);
    fgITSUFrameBox->SetFrameFill(kTRUE);
    fgITSUFrameBox->IncRefCount();

    fgITSUPalette  = new TEveRGBAPalette(0,1);
    fgITSUPalette->IncRefCount();

    fgITSUFrameBoxDead = new TEveFrameBox();
    fgITSUFrameBoxDead->SetAAQuadXZ(-dx/2, dy, -dz/2, dx, dz);
    fgITSUFrameBoxDead->SetFrameColor(kRed);
    fgITSUFrameBoxDead->SetFrameFill(kTRUE);
    fgITSUFrameBoxDead->IncRefCount();
  }


}


/******************************************************************************/

void AliEveITSUModule::SetID(Int_t gid, Bool_t trans)
{
  // Set detector id.

  static const TEveException kEH("AliEveITSUModule::SetID ");

  fID = gid;
  if (!fgStaticInitDone)
  {
    InitStatics();
  }

  SetFrame(fgITSUFrameBox);
  SetPalette(fgITSUPalette);

  RefitPlex();
  ComputeBBox();
  InitMainTrans();
  if (trans)
    SetTrans();

}

/******************************************************************************/

void AliEveITSUModule::SetDigitInQuad(AliITSMFTDigitPix *pDig)
{
  // Sets a digit from source in a visualization structure - called quads.


  if (!fAtLeastOneDigit) {
    Reset(kQT_RectangleXZFixedY, kFALSE, 32);
    fAtLeastOneDigit = kTRUE;
  }

  Float_t x,z;
  fSegm->DetToLocal(pDig->GetCoord2(),pDig->GetCoord1(),x,z);

  AddQuad(x-fDpx/2, z-fDpz/2, fDpx, fDpz);
  QuadId(pDig);

  Int_t intSignal = pDig->GetSignalPix();
  QuadValue(intSignal);
  if (fgITSUPalette->GetMaxVal()<intSignal) {
    fgITSUPalette->SetMax(intSignal);
    fgITSUPalette->MinMaxValChanged();
  }


}

/******************************************************************************/

void AliEveITSUModule::SetTrans()
{
  // Set transformation matrix 
   
  const TGeoHMatrix *mat = fGM->GetMatrixSens(fID);
  fMainTrans->SetFrom(*mat);

}


/******************************************************************************/

void AliEveITSUModule::Print(Option_t* ) const
{
  // Print object summary information.

  printf("AliEveITSUModule: ModuleId: %d, layer %d, ladder %d, detector %d\n",
         fID, fkLayer, fkLadder, fkDetector);

}

/******************************************************************************/

void AliEveITSUModule::DigitSelected(Int_t idx)
{
  // Override secondary select (alt-click) from TEveQuadSet.
 
  //  for (Int_t i=0;i<7;i++) {
  //    idx=i;
  DigitBase_t *qb  = GetDigit(idx);
  TObject     *obj = GetId(idx);
  AliITSMFTDigitPix *pDig  = dynamic_cast<AliITSMFTDigitPix*>(obj);
  printf("AliEveITSUModule::QuadSelected "); 
  printf("  idx=%d, value=%d, obj=0x%lx, digit=0x%lx\n",
	 idx, qb->fValue, (ULong_t)obj, (ULong_t)pDig);
  if (pDig) { 
    Float_t x,z;
    fSegm->DetToLocal(pDig->GetCoord2(),pDig->GetCoord1(),x,z);
    printf(" Digit info: mod|lay/lad/det=%d|%d/%d/%d; row/col=%3d/%4d; \n",
	   fID,fkLayer,fkLadder,fkDetector,
	   pDig->GetCoord2(),pDig->GetCoord1());
    printf(" local (x,z)=(%.4lf,%.4lf)cm; signal:%5d e-;  generated by tracks ",
	   x,z,pDig->GetSignalPix()); 
    for (int itr=0;itr<pDig->GetNTracks();itr++)  
      if (pDig->GetTrack(itr)>=0) printf(" %5d",pDig->GetTrack(itr)); printf("\n");
  }
  
} 
