/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  base class for the ITS elements in the HLT EVE display

#include "AliHLTEveITS.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliEveHOMERManager.h"
#include "TEveManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEvePointSet.h"
#include "ITS/AliHLTITSClusterDataFormat.h"
#include "AliITSRecPoint.h"
#include "TCanvas.h"

ClassImp(AliHLTEveITS)

AliHLTEveITS::AliHLTEveITS(TString name) : 
AliHLTEveBase(),
  fName(name), 
  fPointSet(NULL)
{
  // Constructor.

  SetDetector(fName);
}

AliHLTEveITS::~AliHLTEveITS()
{
  //Destructor
  if(fPointSet)
    delete fPointSet;
  fPointSet = NULL;
}

void AliHLTEveITS::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
 
  if ( block->GetDataType().CompareTo("ROOTHIST") == 0 ) {
    if(!fCanvas) {  
      fCanvas = CreateCanvas(Form("%s QA",fName.Data()), Form("%s QA", fName.Data()));
      fCanvas->Divide(3, 3);
      SetMaxHistograms(9);
    }
    AddHistogramsToCanvas( block , fCanvas, fHistoCount);
  } 

  else if ( block->GetDataType().CompareTo("CLUSTERS") == 0 ) {
    if(!fPointSet) {
      fPointSet = CreatePointSet(fName);
      fEventManager->GetEveManager()->AddElement(fPointSet);
    }
    ProcessClusters(block, fPointSet);
  }
}


void AliHLTEveITS::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  if(fPointSet) fPointSet->ElementChanged();
}

void AliHLTEveITS::ResetElements() {
  //See header file for documentation
  fHistoCount = 0;
  if(fPointSet) fPointSet->Reset();
}

TEvePointSet * AliHLTEveITS::CreatePointSet(TString name) {
  //See header file for documentation
  TEvePointSet * ps = new TEvePointSet(name.Data());
  SetUpPointSet(ps);
  return ps;
} 

void AliHLTEveITS::SetUpPointSet(TEvePointSet * ps ) {
  //See header file for documentation
  ps->SetMainColor(kBlack);
  ps->SetMarkerStyle((Style_t)kFullDotMedium);
}

void AliHLTEveITS::ProcessClusters(AliHLTHOMERBlockDesc * block, TEvePointSet * cont ) {
  //See header file for documentation
  AliHLTITSClusterData *cd = reinterpret_cast<AliHLTITSClusterData*> (block->GetData());
  UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);
  
  if ( cd->fSpacePointCnt != 0 ) {
    for (UInt_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTITSSpacePointData)) {
      AliHLTITSSpacePointData *sp = reinterpret_cast<AliHLTITSSpacePointData*> (data);
  
      Int_t lab[4]   = {0,0,0,0};
      Float_t hit[6] = {0,0,0,0,0,0};
      Int_t info[3]  = {0,0,0};
 				 
      lab[0]  = sp->fTracks[0];
      lab[1]  = sp->fTracks[1];
      lab[2]  = sp->fTracks[2];
      lab[3]  = sp->fIndex;
      hit[0]  = sp->fY;
      hit[1]  = sp->fZ;
      hit[2]  = sp->fSigmaY2;
      hit[3]  = sp->fSigmaZ2;
      hit[4]  = sp->fQ;
      hit[5]  = sp->fSigmaYZ;
      info[0] = sp->fNy;
      info[1] = sp->fNz;
      info[2] = sp->fLayer;
      
      Float_t xyz[3];
      AliITSRecPoint recpoint(lab,hit,info);
      recpoint.GetGlobalXYZ(xyz);

      cont->SetNextPoint(xyz[0], xyz[1], xyz[2]);
    }
  }
}
