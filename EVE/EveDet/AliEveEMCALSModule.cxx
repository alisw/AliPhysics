/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <iostream>
#include <vector>

#include <TString.h>
#include <TEveFrameBox.h>
#include <TEveQuadSet.h>
#include <TEvePointSet.h>
#include <TEveRGBAPalette.h>

#include "AliEveEMCALSModule.h"

class TEveTrans;
class TEveElement;
class TClonesArray;
class TStyle;
class TBuffer3DTypes;
class TBuffer3D;
class TVirtualPad;
class TVirtualViewer3D;

class AliEveEMCALData;
class AliEveEMCALSModuleData;

/// \cond CLASSIMP
ClassImp(AliEveEMCALSModule) ;
/// \endcond

Bool_t           AliEveEMCALSModule::fgStaticInit = kFALSE;

Float_t          AliEveEMCALSModule::fgSMBigBBox   [3];
Float_t          AliEveEMCALSModule::fgSMSmallBBox [3];
Float_t          AliEveEMCALSModule::fgSMDCalBBox  [3];
Float_t          AliEveEMCALSModule::fgSMSmallDBBox[3];

TEveFrameBox*    AliEveEMCALSModule::fgFrameBigBox    = 0;
TEveFrameBox*    AliEveEMCALSModule::fgFrameSmallBox  = 0;
TEveFrameBox*    AliEveEMCALSModule::fgFrameDCalBox   = 0;
TEveFrameBox*    AliEveEMCALSModule::fgFrameSmallDBox = 0;

TEveRGBAPalette* AliEveEMCALSModule::fgFrameDigPalette = 0;
TEveRGBAPalette* AliEveEMCALSModule::fgFrameCluPalette = 0;

///
/// Constructor
///
//______________________________________________________________________________
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
  Char_t name[256];
  
  if (smid < 10) 
    snprintf(name,256,"Full Super Module %02d",smid);
  else 
    snprintf(name,256,"Half Super Module %02d",smid);
  
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

///
/// Copy constructor
///
//______________________________________________________________________________
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
  Char_t name[256];
  if (fSModuleID < 10) 
    snprintf(name,256,"Full Super Module %02d",fSModuleID);
  else 
    snprintf(name,256,"Half Super Module %02d",fSModuleID);
  
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

///
/// Destructor.
///
AliEveEMCALSModule::~AliEveEMCALSModule()
{
  fPointSet->DecDenyDestroy();
  fQuadSet ->DecDenyDestroy();
  fQuadSet2->DecDenyDestroy();

  if(fEMCALData) fEMCALData->DecRefCount();
}

///
/// Bounding box, Framebox and Palette
///
//______________________________________________________________________________
void AliEveEMCALSModule::InitStatics(AliEveEMCALSModuleData* md)
{
  if (fgStaticInit) return;
  
  fgStaticInit = kTRUE;

  md->GetSModuleBigBox   (fgSMBigBBox   [0], fgSMBigBBox   [1], fgSMBigBBox   [2]);
  md->GetSModuleSmallBox (fgSMSmallBBox [0], fgSMSmallBBox [1], fgSMSmallBBox [2]);  
  md->GetSModuleDCalBox  (fgSMDCalBBox  [0], fgSMDCalBBox  [1], fgSMDCalBBox  [2]);
  md->GetSModuleSmallDBox(fgSMSmallDBBox[0], fgSMSmallDBBox[1], fgSMSmallDBBox[2]);

  fgFrameBigBox = new TEveFrameBox();
  fgFrameBigBox->SetAABoxCenterHalfSize   (0, 0, 0, fgSMBigBBox   [0], fgSMBigBBox   [1], fgSMBigBBox   [2]);
  fgFrameBigBox->SetFrameColor((Color_t)10);
  fgFrameBigBox->IncRefCount();

  fgFrameSmallBox = new TEveFrameBox();
  fgFrameSmallBox->SetAABoxCenterHalfSize (0, 0, 0, fgSMSmallBBox [0], fgSMSmallBBox [1], fgSMSmallBBox [2]);
  fgFrameSmallBox->SetFrameColor((Color_t)10);
  fgFrameSmallBox->IncRefCount();

  fgFrameDCalBox = new TEveFrameBox();
  fgFrameDCalBox->SetAABoxCenterHalfSize  (0, 0, 0, fgSMDCalBBox  [0], fgSMDCalBBox  [1], fgSMDCalBBox  [2]);
  fgFrameDCalBox->SetFrameColor((Color_t)10);
  fgFrameDCalBox->IncRefCount();
  
  fgFrameSmallDBox = new TEveFrameBox();
  fgFrameSmallDBox->SetAABoxCenterHalfSize(0, 0, 0, fgSMSmallDBBox[0], fgSMSmallDBBox[1], fgSMSmallDBBox[2]);
  fgFrameSmallDBox->SetFrameColor((Color_t)10);
  fgFrameSmallDBox->IncRefCount();  
  
  fgFrameDigPalette = new TEveRGBAPalette(0,512);
  fgFrameDigPalette->SetLimits(0, 1024);
  fgFrameDigPalette->IncRefCount();
  
  fgFrameCluPalette  = new TEveRGBAPalette(0,512);
  fgFrameCluPalette->SetLimits(0, 1024);
  fgFrameCluPalette->IncRefCount();
}

///
/// Cluster point size
///
//______________________________________________________________________________
void AliEveEMCALSModule::SetClusterSize(Int_t size)
{
  fClusterSize = TMath::Max(1, size);
}

///
/// Hit point size
///
//______________________________________________________________________________
void AliEveEMCALSModule::SetHitSize(Int_t size)
{
  fHitSize = TMath::Max(1, size);
}

///
/// Set source of data.
///
//______________________________________________________________________________
void AliEveEMCALSModule::SetDataSource(AliEveEMCALData * data)
{
  if ( data == fEMCALData ) return;
  
  if ( fEMCALData ) fEMCALData->DecRefCount();
  
  fEMCALData = data;
  
  if ( fEMCALData ) fEMCALData->IncRefCount();

  // Get pointer on SM data
  fEMCALSModuleData = GetSModuleData();
}

///
/// Return source of data.
///
//______________________________________________________________________________
AliEveEMCALSModuleData* AliEveEMCALSModule::GetSModuleData() const
{
  return fEMCALData ? fEMCALData->GetSModuleData(fSModuleID) : 0;
}

///
/// Update hit/digit/cluster representation.
///
//______________________________________________________________________________
void AliEveEMCALSModule::UpdateQuads(Bool_t iHits, Bool_t iDigits, Bool_t iClusters)
{  
  Int_t smId = fEMCALSModuleData->GetSmId();
  
  if (!fgStaticInit)
    InitStatics(fEMCALSModuleData);
  
  //--------------------------------
  // hits --------------------------
  //--------------------------------
  
  if(iHits)
  {    
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
    
    Float_t hitX = 0, hitY = 0, hitZ = 0;

    std::vector< std::vector<Float_t> >  bufferHit;
    
    bufferHit = fEMCALSModuleData->GetHitBuffer();
    
    if(!bufferHit.empty())
    {
      char form[1000];
      
      Int_t nHits = fEMCALSModuleData->GetNHits();
      
      AliDebug(1,Form("nHits: %d", nHits));
      
      Int_t oldSize = fPointSet->GrowFor(nHits);
      
      // Loop over hits
      for (Int_t ih = 0; ih < nHits; ih++) 
      {
        hitX = bufferHit[ih][3];
        hitY = bufferHit[ih][4];
        hitZ = bufferHit[ih][5];
        
        fPointSet->SetPoint(ih,hitX,hitY,hitZ);
        
        snprintf(form,1000,"N=%d", fPointSet->Size());
        fPointSet->SetTitle(form);
        fPointSet->SetMarkerSize(.5);
        fPointSet->SetMarkerColor((Color_t)2);
      }
    }
    else AliDebug(1,Form("There are no hits in SM %d", smId));
  }  
  
  //--------------------------------
  // digits ------------------------
  //--------------------------------
  
  if(iDigits)
  {
    std::vector< std::vector<Double_t> > bufferDigit;
    
    // Define TEveQuadSet for digits
    fQuadSet->SetOwnIds(kTRUE);
    fQuadSet->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix(smId));
    fQuadSet->SetPalette(fgFrameDigPalette);
    
    if     (smId < 10) 
      fQuadSet->SetFrame(fgFrameBigBox   );
    else if(smId < 12)  
      fQuadSet->SetFrame(fgFrameSmallBox );
    else if(smId < 18)  
      fQuadSet->SetFrame(fgFrameDCalBox  );
    else if(smId < 20)  
      fQuadSet->SetFrame(fgFrameSmallDBox);
    
    // Get the digit information from the buffer
    bufferDigit = fEMCALSModuleData->GetDigitBuffer();
    
    if(!bufferDigit.empty())
    {
      Int_t nDigits = fEMCALSModuleData->GetNDigits();
      
      AliDebug(1,Form("nDigits: %d", nDigits) );
      
      // loop over digits
      for (Int_t id = 0; id < nDigits; id++) 
      {
        //	  Int_t iid = (Int_t)bufferDigit[id][0];
        //	  Int_t isupMod = (Int_t)bufferDigit[id][1];
        Double_t iamp = bufferDigit[id][2];
        //Int_t amp = (Int_t)(iamp+0.5); // Why? Let it be float.
        
        //	  Double_t ix = bufferDigit[id][3];
        Double_t iy = bufferDigit[id][4];
        Double_t iz = bufferDigit[id][5];
        
        // Add digit information to the TEveQuadSet
        fQuadSet->AddQuad(iy, iz);
        fQuadSet->QuadValue(iamp);
      } // end digits loop
    }
    else AliDebug(1,Form("There are no digits in SM %d", smId)); 
  }
  
  //----------------------------------
  // clusters ------------------------
  //----------------------------------
  
  if(iClusters)
  {
    std::vector< std::vector<Double_t> > bufferCluster;
    
    // Define TEveQuadSet for clusters
    fQuadSet2->SetOwnIds(kTRUE);
    fQuadSet2->Reset(TEveQuadSet::kQT_RectangleYZFixedDimX, kFALSE, 32);
    fQuadSet2->SetDefWidth (fEMCALSModuleData->GetPhiTileSize());
    fQuadSet2->SetDefHeight(fEMCALSModuleData->GetEtaTileSize());
    fQuadSet2->RefMainTrans().SetFrom(*fEMCALSModuleData->GetSModuleMatrix(smId));
    fQuadSet2->SetPalette(fgFrameCluPalette);
    
    if     (smId < 10) 
      fQuadSet2->SetFrame(fgFrameBigBox   );
    else if(smId < 12)  
      fQuadSet2->SetFrame(fgFrameSmallBox );
    else if(smId < 18)  
      fQuadSet2->SetFrame(fgFrameDCalBox  );
    else if(smId < 20)  
      fQuadSet2->SetFrame(fgFrameSmallDBox);
    
    // Get the cluster information from the buffer
    bufferCluster = fEMCALSModuleData->GetClusterBuffer();
    if(!bufferCluster.empty())
    {
      Int_t nClusters = fEMCALSModuleData->GetNClusters();
      
      AliDebug(1, Form("nClusters: %d", nClusters) );
      
      // loop over clusters
      for (Int_t id = 0; id < nClusters; id++) 
      {
        AliDebug(1,Form("bufferCluster[%d][0]: %f",id, bufferCluster[id][0]));
        AliDebug(1,Form("bufferCluster[%d][1]: %f",id, bufferCluster[id][1]));
        AliDebug(1,Form("bufferCluster[%d][2]: %f",id, bufferCluster[id][2]));
        AliDebug(1,Form("bufferCluster[%d][3]: %f",id, bufferCluster[id][3]));
        AliDebug(1,Form("bufferCluster[%d][4]: %f",id, bufferCluster[id][4]));
        
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
    else AliDebug(1,Form("There are no clusters in SM %d", smId));
  }
}

///
/// Set id of the SM to display.
///
//______________________________________________________________________________
void AliEveEMCALSModule::SetSModuleID(Int_t id)
{
  if (id <  0) id = 0;
  if (id > 20) id = 20;

  fSModuleID = id;
}
