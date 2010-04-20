// EMCAL event display
// Visualization of an EMCAL super module.
//
//  Author: Magali Estienne (magali.estienne@cern.ch)
//  June 30 2008

#include "AliEveEMCALSModule.h"

//#include <vector>

#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>


class Riostream;
//class vector;
class TEveTrans;
class TEveElement;
class TClonesArray;
class TStyle;
class TBuffer3DTypes;
class TBuffer3D;
class TVirtualPad;
class TVirtualViewer3D;
class AliEveEMCALData;
class AliEMCALHit;
class AliEMCALDigit;
class AliEveEMCALSModuleData;

ClassImp(AliEveEMCALSModule)

Bool_t           AliEveEMCALSModule::fgStaticInit = kFALSE;
Float_t          AliEveEMCALSModule::fgSMBigBBox[3];
Float_t          AliEveEMCALSModule::fgSMSmallBBox[3];
TEveFrameBox*    AliEveEMCALSModule::fgFrameBigBox = 0;
TEveFrameBox*    AliEveEMCALSModule::fgFrameSmallBox = 0;
TEveRGBAPalette* AliEveEMCALSModule::fgFrameDigPalette = 0;
TEveRGBAPalette* AliEveEMCALSModule::fgFrameCluPalette = 0;

AliEveEMCALSModule::AliEveEMCALSModule(Int_t smid, const Text_t* n, const Text_t* t) :
  TEveElement(fFrameColor),
  TNamed(n,t),
  TAtt3D(),
  fEMCALData(0),
  fEMCALSModuleData(0),
  fFrameColor((Color_t)10),
  fSModuleID(smid),
  fQuadSet(new TEveQuadSet(n,t)),
  fQuadSet2(new TEveQuadSet(n,t)),
  fPointSet(new TEvePointSet(n)),
  fClusterSize(5),
  fHitSize(5),
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
  TAtt3D(),
  fEMCALData(esm.fEMCALData),
  fEMCALSModuleData(esm.fEMCALSModuleData),
  fFrameColor(esm.fFrameColor),
  fSModuleID(esm.fSModuleID),
  fQuadSet(esm.fQuadSet),
  fQuadSet2(esm.fQuadSet2),
  fPointSet(esm.fPointSet),
  fClusterSize(esm.fClusterSize),
  fHitSize(esm.fHitSize),
  fDebug(esm.fDebug)
{
  // Copy constructor.
  Char_t name[256];
  if (fSModuleID < 10) {
    sprintf(name,"Full Super Module %02d",fSModuleID);
  } else {
    sprintf(name,"Half Super Module %02d",fSModuleID);
  }
  SetName(name);

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
  fQuadSet2->DecDenyDestroy();

  if(fEMCALData) fEMCALData->DecRefCount();
}

//______________________________________________________________________________
void AliEveEMCALSModule::DropData() const
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
void AliEveEMCALSModule::InitStatics(AliEveEMCALSModuleData* md)
{
  //
  // Bounding box, Framebox and Palette
  //

  if (fgStaticInit) return;
  fgStaticInit = kTRUE;

  md->GetSModuleBigBox(fgSMBigBBox[0],fgSMBigBBox[1], fgSMBigBBox[2]);
  md->GetSModuleSmallBox(fgSMSmallBBox[0],fgSMSmallBBox[1], fgSMSmallBBox[2]);

  fgFrameBigBox = new TEveFrameBox();
  fgFrameBigBox->SetAABoxCenterHalfSize(0, 0, 0, fgSMBigBBox[0], fgSMBigBBox[1], fgSMBigBBox[2]);
  fgFrameBigBox->SetFrameColor((Color_t)10);
  fgFrameBigBox->IncRefCount();
  fgFrameDigPalette = new TEveRGBAPalette(0,512);
  fgFrameDigPalette->SetLimits(0, 1024);
  fgFrameDigPalette->IncRefCount();

  fgFrameSmallBox = new TEveFrameBox();
  fgFrameSmallBox->SetAABoxCenterHalfSize(0, 0, 0, fgSMSmallBBox[0], fgSMSmallBBox[1], fgSMSmallBBox[2]);
  fgFrameSmallBox->SetFrameColor((Color_t)10);
  fgFrameSmallBox->IncRefCount();
  fgFrameCluPalette  = new TEveRGBAPalette(0,512);
  fgFrameCluPalette->SetLimits(0, 1024);
  fgFrameCluPalette->IncRefCount();
}

//______________________________________________________________________________
void AliEveEMCALSModule::SetClusterSize(Int_t size)
{
  //
  // Cluster point size
  //

  fClusterSize = TMath::Max(1, size);
}

//______________________________________________________________________________
void AliEveEMCALSModule::SetHitSize(Int_t size)
{
  //
  // hit point size
  //

  fHitSize = TMath::Max(1, size);
}

//______________________________________________________________________________
void AliEveEMCALSModule::SetDataSource(AliEveEMCALData* const data)
{
  //
  // Set source of data.
  //

  if (data == fEMCALData) return;
  if(fEMCALData) fEMCALData->DecRefCount();
  fEMCALData = data;
  if(fEMCALData) fEMCALData->IncRefCount();

  // Get pointer on SM data
  fEMCALSModuleData = GetSModuleData();
}

//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALSModule::GetSModuleData() const
{
  //
  // Return source of data.
  //

  return fEMCALData ? fEMCALData->GetSModuleData(fSModuleID) : 0;
}

//______________________________________________________________________________
void AliEveEMCALSModule::UpdateQuads()
{
  //
  // Update hit/digit/cluster representation.
  //

  vector< vector<Double_t> > bufferDigit;
  vector< vector<Double_t> > bufferCluster;
  vector< vector<Float_t> > bufferHit;
  Int_t nDigits, nClusters, nHits, oldSize;
  Float_t hitX, hitY, hitZ;
  Int_t smId = fEMCALSModuleData->GetSmId();

  //--------------------------
  // Hits from runloader
  //--------------------------
  fPointSet->Reset();

  /*
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
  
  */

  if (fEMCALSModuleData != 0) {

    if (!fgStaticInit)
      InitStatics(fEMCALSModuleData);

    // digits ------------------------

    // Define TEveQuadSet for digits
    fQuadSet->SetOwnIds(kTRUE);
    fQuadSet->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix());
    fQuadSet->SetPalette(fgFrameDigPalette);
    if(smId<fEMCALSModuleData->GetNsmf()) 
      fQuadSet->SetFrame(fgFrameBigBox);
    else fQuadSet->SetFrame(fgFrameSmallBox);

    // Get the digit information from the buffer
    bufferDigit = fEMCALSModuleData->GetDigitBuffer();
    if(!bufferDigit.empty())
      {
	nDigits = fEMCALSModuleData->GetNDigits();
	if(fDebug>1) cout << "nDigits: " << nDigits << endl;
	// loop over digits
	for (Int_t id = 0; id < nDigits; id++) {
	  //	  Int_t iid = (Int_t)bufferDigit[id][0];
	  //	  Int_t isupMod = (Int_t)bufferDigit[id][1];
	  Double_t iamp = bufferDigit[id][2];
	  Int_t amp = (Int_t)(iamp+0.5);
	  //	  Double_t ix = bufferDigit[id][3];
	  Double_t iy = bufferDigit[id][4];
	  Double_t iz = bufferDigit[id][5];
	  
	  // Add digit information to the TEveQuadSet
	  fQuadSet->AddQuad(iy, iz);
	  fQuadSet->QuadValue(amp);
	} // end digits loop
      }
    else { if (fDebug) printf("There is no digits in SM %d \n", smId); }

    // hits --------------------------
    bufferHit = fEMCALSModuleData->GetHitBuffer();
    if(!bufferHit.empty())
      {
	char form[1000];
	nHits = fEMCALSModuleData->GetNHits();
	if(fDebug>1) cout << "nHits: " << nHits << endl;
	oldSize = fPointSet->GrowFor(nHits);
	// Loop over hits
	for (Int_t ih = 0; ih < nHits; ih++) {
	  hitX = bufferHit[ih][3];
	  hitY = bufferHit[ih][4];
	  hitZ = bufferHit[ih][5];
	  fPointSet->SetPoint(ih,hitX,hitY,hitZ);
	  sprintf(form,"N=%d", fPointSet->Size());
	  fPointSet->SetTitle(form);
	  fPointSet->SetMarkerSize(.5);
	  fPointSet->SetMarkerColor((Color_t)2);
	}
      }
    else { if (fDebug) printf("There is no hits in SM %d \n", smId); }

    // clusters ------------------------

    // Define TEveQuadSet for clusters
    fQuadSet2->SetOwnIds(kTRUE);
    fQuadSet2->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet2->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet2->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet2->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix());
    fQuadSet2->SetPalette(fgFrameCluPalette);
    if(smId<fEMCALSModuleData->GetNsmf()) 
      fQuadSet2->SetFrame(fgFrameBigBox);
    else fQuadSet2->SetFrame(fgFrameSmallBox);

    // Get the cluster information from the buffer
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
	  //	  Int_t isupMod = (Int_t)bufferCluster[id][0];
	  Double_t iamp = bufferCluster[id][1];
	  Int_t amp = (Int_t)(iamp+0.5);
	  //	  Double_t ix = bufferCluster[id][2];
	  Double_t iy = bufferCluster[id][3];
	  Double_t iz = bufferCluster[id][4];
	  
	  // Add cluster information to the TEveQuadSet
	  fQuadSet2->AddQuad(iy, iz);
	  fQuadSet2->QuadValue(amp);
	  //      fQuadSet2->QuadId(iid);

	} // end clusters loop
      }
    else { if (fDebug) printf("There is no clusters in SM %d \n", smId); }

  } // end if (fEMCALSModuleData != 0)

}

//______________________________________________________________________________
void AliEveEMCALSModule::SetSModuleID(Int_t id)
{
  //
  // Set id of the SM to display.
  //

  if (id <  0) id = 0;
  if (id > 12) id = 12;

  fSModuleID = id;
}
