//*************************************************************************
// EMCAL event display
// Visualization of an EMCAL super module.
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008
//*************************************************************************

#include <Riostream.h>
#include <vector>

#include <TEveTrans.h>
#include <TEveElement.h>
#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEvePointSet.h>
#include <TClonesArray.h>
#include <TVectorT.h>
#include <TStyle.h>
#include <TBuffer3DTypes.h>
#include <TBuffer3D.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>
#include <TEveRGBAPalette.h>

#include "AliEveEMCALData.h"
#include "AliEveEMCALSModule.h"
#include "AliEveEMCALSModuleData.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"


ClassImp(AliEveEMCALSModule)

TEveFrameBox*    AliEveEMCALSModule::fFrameBigBox = 0;
TEveFrameBox*    AliEveEMCALSModule::fFrameSmallBox = 0;
TEveRGBAPalette* AliEveEMCALSModule::fFrameDigPalette = 0;
TEveRGBAPalette* AliEveEMCALSModule::fFrameCluPalette = 0;

AliEveEMCALSModule::AliEveEMCALSModule(Int_t smid, const Text_t* n, const Text_t* t) :
  TEveElement(fFrameColor),
  TNamed(n,t),
  fEMCALData(0),
  fEMCALSModuleData(0),
  fFrameColor((Color_t)10),
  fRTS(1),
  fSModuleID(smid),
  fQuadSet(new TEveQuadSet(n,t)),
  fQuadSet2(new TEveQuadSet(n,t)),
  fPointSet(new TEvePointSet(n)),
  fThreshold(0),
  fMaxVal(4096),
  fClusterSize(5),
  fHitSize(5),
  fColorArray(0),
  fDebug(0)
{
  // Constructor.
  Char_t name[256];
  if (smid < 10) {
    sprintf(name,"Full Super Module %02d",smid);
  } else {
    sprintf(name,"Half Super Module %02d",smid);
  }
  SetName(name);

  for(Int_t i=0; i<3; i++) fSMBigBBox[i] = 0.;
  for(Int_t i=0; i<3; i++) fSMSmallBBox[i] = 0.;
  for(Int_t i=0; i<3; i++) fSMBBoxCenter[i] = 0.;

  // Hits
  fPointSet->IncDenyDestroy();
  AddElement(fPointSet);
  // Digits
  fQuadSet->IncDenyDestroy();
  AddElement(fQuadSet);
  // Clusters
  fQuadSet2->IncDenyDestroy();
  AddElement(fQuadSet2);

}

AliEveEMCALSModule::AliEveEMCALSModule(const AliEveEMCALSModule &esm) :
  TEveElement(fFrameColor),
  TNamed(),
  fEMCALData(esm.fEMCALData),
  fEMCALSModuleData(esm.fEMCALSModuleData),
  fFrameColor(esm.fFrameColor),
  fRTS(esm.fRTS),
  fSModuleID(esm.fSModuleID),
  fQuadSet(esm.fQuadSet),
  fQuadSet2(esm.fQuadSet2),
  fPointSet(esm.fPointSet),
  fThreshold(esm.fThreshold),
  fMaxVal(esm.fMaxVal),
  fClusterSize(esm.fClusterSize),
  fHitSize(esm.fHitSize),
  fColorArray(esm.fColorArray),
  fDebug(esm.fDebug)
{
  // Constructor.
  Char_t name[256];
  if (fSModuleID < 10) {
    sprintf(name,"Full Super Module %02d",fSModuleID);
  } else {
    sprintf(name,"Half Super Module %02d",fSModuleID);
  }
  SetName(name);

  for(Int_t i=0; i<3; i++) fSMBigBBox[i] = 0.;
  for(Int_t i=0; i<3; i++) fSMSmallBBox[i] = 0.;
  for(Int_t i=0; i<3; i++) fSMBBoxCenter[i] = 0.;

  // Hits
  fPointSet->IncDenyDestroy();
  AddElement(fPointSet);
  // Digits
  fQuadSet->IncDenyDestroy();
  AddElement(fQuadSet);
  // Clusters
  fQuadSet2->IncDenyDestroy();
  AddElement(fQuadSet2);

}

AliEveEMCALSModule::~AliEveEMCALSModule()
{
  //
  // Destructor.
  //

  fPointSet->DecDenyDestroy();
  fQuadSet->DecDenyDestroy();

  if(fEMCALData) fEMCALData->DecRefCount();
  if(fFrameBigBox)   fFrameBigBox->Delete();
  if(fFrameSmallBox) fFrameSmallBox->Delete();
  if(fFrameDigPalette) fFrameDigPalette->Delete();
  if(fFrameCluPalette) fFrameCluPalette->Delete();

}

//______________________________________________________________________________
void AliEveEMCALSModule::DropData()
{
  //
  // release the sm data
  //

//   fNDigits   = 0;
//   fNClusters = 0;
//   fNHits     = 0;

  return;

}

//______________________________________________________________________________
void AliEveEMCALSModule::ComputeBBox()
{
  //
  // bounding box
  //

  fEMCALSModuleData->GetSModuleBigBox(fSMBigBBox[0],fSMBigBBox[1], fSMBigBBox[2]);
  fEMCALSModuleData->GetSModuleSmallBox(fSMSmallBBox[0],fSMSmallBBox[1], fSMSmallBBox[2]);

//           if (fgStaticInitDone) return;
//           fgStaticInitDone = kTRUE;

  fFrameBigBox = new TEveFrameBox();
  fFrameBigBox->SetAABoxCenterHalfSize(0, 0, 0, fSMBigBBox[0], fSMBigBBox[1], fSMBigBBox[2]);
  fFrameBigBox->SetFrameColor((Color_t)10);
  fFrameDigPalette = new TEveRGBAPalette(0,512);
  fFrameDigPalette->SetLimits(0, 1024);
  fFrameDigPalette->IncRefCount();

  fFrameSmallBox = new TEveFrameBox();
  fFrameSmallBox->SetAABoxCenterHalfSize(0, 0, 0, fSMSmallBBox[0], fSMSmallBBox[1], fSMSmallBBox[2]);
  fFrameSmallBox->SetFrameColor((Color_t)10);
  fFrameCluPalette  = new TEveRGBAPalette(0,512);
  fFrameCluPalette->SetLimits(0, 1024);
  fFrameCluPalette->IncRefCount();

  BBoxInit();

  fBBox[0] = - 2*fSMBigBBox[0];
  fBBox[1] = + 2*fSMBigBBox[0];
  fBBox[2] = - 2*fSMBigBBox[1];
  fBBox[3] = + 2*fSMBigBBox[1];
  fBBox[4] = - 2*fSMBigBBox[2];
  fBBox[5] = + 2*fSMBigBBox[2];

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetThreshold(Short_t t)
{
  //
  // digits amplitude threshold
  //

  fThreshold = TMath::Min(t, (Short_t)(fMaxVal - 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetMaxVal(Int_t mv)
{
  //
  // digits amplitude maximum value
  //

  fMaxVal = TMath::Max(mv, (Int_t)(fThreshold + 1));
  ClearColorArray();
  IncRTS();

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetClusterSize(Int_t size)
{
  //
  // cluster point size
  //

  fClusterSize = TMath::Max(1, size);
  IncRTS();

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetHitSize(Int_t size)
{
  //
  // hit point size
  //

  fHitSize = TMath::Max(1, size);
  IncRTS();

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetupColor(Int_t val, UChar_t* pixel) const
{
  //
  // RGBA color for amplitude "val"
  //

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  TEveUtil::TEveUtil::ColorFromIdx(gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin)), pixel);

}

//______________________________________________________________________________
Int_t AliEveEMCALSModule::ColorIndex(Int_t val) const
{
  //
  // index color
  //

  if(val < fThreshold) val = fThreshold;
  if(val > fMaxVal)    val = fMaxVal;

  Float_t div  = TMath::Max(1, fMaxVal - fThreshold);
  Int_t   nCol = gStyle->GetNumberOfColors();
  Int_t   cBin = (Int_t) TMath::Nint(nCol*(val - fThreshold)/div);

  return gStyle->GetColorPalette(TMath::Min(nCol - 1, cBin));

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetupColorArray() const
{
  //
  // build array of colors
  //

  if(fColorArray)
    return;

  fColorArray = new UChar_t [4 * (fMaxVal - fThreshold + 1)];
  UChar_t* p = fColorArray;
  for(Int_t v=fThreshold; v<=fMaxVal; ++v, p+=4)
    SetupColor(v, p);

}

//______________________________________________________________________________
void AliEveEMCALSModule::ClearColorArray()
{
  //
  // delete array of colors
  //

  if(fColorArray) {
    delete [] fColorArray;
    fColorArray = 0;
  }
}

//______________________________________________________________________________
void AliEveEMCALSModule::SetDataSource(AliEveEMCALData* data)
{
  // Set source of data.

  if (data == fEMCALData) return;
  if(fEMCALData) fEMCALData->DecRefCount();
  fEMCALData = data;
  if(fEMCALData) fEMCALData->IncRefCount();

  fEMCALSModuleData = GetSModuleData();

  IncRTS();
}

//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALSModule::GetSModuleData() const
{
  // Return source of data.

  return fEMCALData ? fEMCALData->GetSModuleData(fSModuleID) : 0;
}

//______________________________________________________________________________
void AliEveEMCALSModule::UpdateQuads()
{
  // Update hit/digit/cluster representation.

  vector< vector<Float_t> > bufferDigit;
  vector< vector<Float_t> > bufferCluster;
  vector< vector<Float_t> > bufferHit;
  Float_t x0, y0, z, w, h, clsq;
  Int_t charge, cathode, nDigits, nClusters, nHits, oldSize, ic1, ic2;
  Double_t clsX, clsY, clsZ;
  Float_t hitX, hitY, hitZ;
  Int_t smId = fEMCALSModuleData->GetSmId();

  //--------------------------
  // Hits from runloader
  //--------------------------
  fPointSet->Reset();

  TEvePointSet* points = fEMCALData->GetPointSetData();
  char form[1000];
  if(points){
    sprintf(form,"N=%d", points->Size());
    points->SetTitle(form);
    points->SetMarkerSize(.5);
    points->SetMarkerColor((Color_t)2);
    fPointSet->AddElement(points);
  }
  else {printf("There is no hits in Runloader \n"); }
  
  if (fEMCALSModuleData != 0) {
    
    // digits ------------------------
    fQuadSet->SetOwnIds(kTRUE);
    fQuadSet->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix());
    fQuadSet->SetPalette(fFrameDigPalette);
    if(smId<fEMCALSModuleData->GetNsmf()) 
      fQuadSet->SetFrame(fFrameBigBox);
    else fQuadSet->SetFrame(fFrameSmallBox);

    bufferDigit = fEMCALSModuleData->GetDigitBuffer();
    if(!bufferDigit.empty())
      {
	nDigits = fEMCALSModuleData->GetNDigits();
	if(fDebug>1) cout << "nDigits: " << nDigits << endl;
	// loop over digits
	for (Int_t id = 0; id < nDigits; id++) {
	  if(fDebug>1) {
	    cout << "bufferDigit[" << id << "][0]: " << bufferDigit[id][0] << endl;
	    cout << "bufferDigit[" << id << "][1]: " << bufferDigit[id][1] << endl;
	    cout << "bufferDigit[" << id << "][2]: " << bufferDigit[id][2] << endl;
	    cout << "bufferDigit[" << id << "][3]: " << bufferDigit[id][3] << endl;
	    cout << "bufferDigit[" << id << "][4]: " << bufferDigit[id][4] << endl;
	  }
	  Int_t iid = bufferDigit[id][0];
	  Int_t isupMod = bufferDigit[id][1];
	  Float_t iamp = bufferDigit[id][2];
	  Float_t ix = bufferDigit[id][3];
	  Float_t iy = bufferDigit[id][4];
	  Float_t iz = bufferDigit[id][5];
	  
	  fQuadSet->AddQuad(iy, iz);
	  fQuadSet->QuadValue(iamp);
      //      fQuadSet->QuadId(iid);

	} // end digits loop
      }
    else { printf("There is no digits in SM %d \n", smId); }

    // hits --------------------------
    bufferHit = fEMCALSModuleData->GetDigitBuffer();
    if(!bufferHit.empty())
      {
	nHits = fEMCALSModuleData->GetNHits();
	if(fDebug>1) cout << "nHits: " << nHits << endl;
	oldSize = fPointSet->GrowFor(nHits);
	// Loop over hits
	for (Int_t ih = 0; ih < nHits; ih++) {
	  hitX = bufferHit[ih][2];
	  hitY = bufferHit[ih][3];
	  hitZ = bufferHit[ih][4];
	  fPointSet->SetPoint(ih,hitX,hitY,hitZ);
	}
      }
    else {printf("There is no hits in SM %d \n", smId); }

    // clusters ------------------------
    fQuadSet2->SetOwnIds(kTRUE);
    fQuadSet2->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet2->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet2->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet2->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix());
    fQuadSet2->SetPalette(fFrameCluPalette);
    if(smId<fEMCALSModuleData->GetNsmf()) 
      fQuadSet2->SetFrame(fFrameBigBox);
    else fQuadSet2->SetFrame(fFrameSmallBox);

    bufferCluster = fEMCALSModuleData->GetClusterBuffer();
    if(!bufferCluster.empty())
      {
	nClusters = fEMCALSModuleData->GetNClusters();
	if(fDebug>1) cout << "nClusters: " << nClusters << endl;
	// loop over clusters
	for (Int_t id = 0; id < nClusters; id++) {
	  if(fDebug>1) {
	    cout << "bufferCluster[" << id << "][0]: " << bufferCluster[id][0] << endl;
	    cout << "bufferCluster[" << id << "][1]: " << bufferCluster[id][1] << endl;
	    cout << "bufferCluster[" << id << "][2]: " << bufferCluster[id][2] << endl;
	    cout << "bufferCluster[" << id << "][3]: " << bufferCluster[id][3] << endl;
	    cout << "bufferCluster[" << id << "][4]: " << bufferCluster[id][4] << endl;
	  }
	  Int_t isupMod = bufferCluster[id][0];
	  Float_t iamp = bufferCluster[id][1];
	  Float_t ix = bufferCluster[id][2];
	  Float_t iy = bufferCluster[id][3];
	  Float_t iz = bufferCluster[id][4];
	  
	  fQuadSet2->AddQuad(iy, iz);
	  fQuadSet2->QuadValue(iamp);
	  //      fQuadSet2->QuadId(iid);

	} // end clusters loop
      }
    else { printf("There is no clusters in SM %d \n", smId); }

  } // end if (fEMCALSModuleData != 0)

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetSModuleID(Int_t id)
{
  // Set id of chamber to display.

  if (id <  0) id = 0;
  if (id > 12) id = 12;

  fSModuleID = id;
  IncRTS();
}
