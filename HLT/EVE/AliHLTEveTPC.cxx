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
/// @brief  TPC processor for the HLT EVE display

#include "AliHLTEveTPC.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliHLTEveBase.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TColor.h"
#include "TMath.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "TH1F.h"

ClassImp(AliHLTEveTPC)

AliHLTEveTPC::AliHLTEveTPC() : 
  AliHLTEveBase("TPC Clusters"), 
  fEveClusters(NULL),
  fEveColClusters(NULL),
  fNColorBins(15), 
  fHistCharge(NULL), 
  fHistQMax(NULL), 
  fHistQMaxOverCharge(NULL)
{
  // Constructor.
}

AliHLTEveTPC::~AliHLTEveTPC()
{
  //Destructor
  if(fEveColClusters)
    delete fEveColClusters;
  fEveColClusters = NULL;

  if(fEveClusters)
    delete fEveClusters;
  fEveClusters = NULL;

  if(fHistQMaxOverCharge)
    delete fHistQMaxOverCharge;
  fHistQMaxOverCharge = NULL;

  if(fHistQMax)
    delete fHistQMax;
  fHistQMax = NULL;

  if(fHistCharge)
    delete fHistCharge;
  fHistCharge = NULL;

}


void AliHLTEveTPC::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {

    if(!fEveClusters){	  
      fEveClusters = CreatePointSet();
      //fEventManager->GetEveManager()->AddElement(fEveClusters);
    } 
    
    if(!fEveColClusters){
      fEveColClusters = CreatePointSetArray();
      AddElement(fEveColClusters);
    } 
    
    ProcessClusters(block, fEveClusters, fEveColClusters);
   
  }
  
//   else if ( ! block->GetDataType().CompareTo("HWCL_ALT") ) {
//     if(!gTPCTestClusters){	  
      
//       gTPCTestClusters = new TEvePointSet("TPC Clusters Test");
//       //ggTPCTestClusters->ApplyVizTag("TPC Clusters");
//       gTPCTestClusters->SetMainColor(kBlue);
//       gTPCTestClusters->SetMarkerStyle((Style_t)kFullDotSmall);
//       gEve->AddElement(gTPCTestClusters);
//     }
    
//     processTPCClusters(block, gTPCTestClusters);
//     gTPCTestClusters->ElementChanged();
//   }
  
  
}

TEvePointSet * AliHLTEveTPC::CreatePointSet() {
  //See header file for documentation

  TEvePointSet * ps = new TEvePointSet("TPC Clusters");
  ps->SetMainColor(kRed);
  ps->SetMarkerStyle((Style_t)kFullDotSmall);
 
  return ps;

}

TEvePointSetArray * AliHLTEveTPC::CreatePointSetArray(){
  //See header file for documentation

  TEvePointSetArray * cc = new TEvePointSetArray("TPC Clusters Colorized");
  cc->SetMainColor(kRed);
  cc->SetMarkerStyle(4); // antialiased circle
  cc->SetMarkerSize(0.4);
  cc->InitBins("Cluster Charge", fNColorBins, 0., fNColorBins*20.);
  
  const Int_t nCol = TColor::GetNumberOfColors();
  
  for (Int_t ii = 0; ii < fNColorBins + 1; ++ii) {
    cc->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (fNColorBins+2)));
  }

  return cc;
     
}


void AliHLTEveTPC::UpdateElements() {
  //See header file for documentation

  if(fEveClusters) fEveClusters->ResetBBox();

  if (fHistQMax || fHistQMaxOverCharge || fHistCharge )
    // DrawHistograms();
 
  if(fCanvas) fCanvas->Update();
 
  // if(fEveColClusters){
  //  for (Int_t ib = 0; ib <= fNColorBins + 1; ++ib) {
  //    fEveColClusters->GetBin(ib)->ResetBBox();
  //  }
  // }

}

void AliHLTEveTPC::ResetElements(){
  //See header file for documentation

  if(fEveClusters) fEveClusters->Reset();
  if(fEveColClusters){
    for (Int_t ib = 0; ib <= fNColorBins + 1; ++ib) {
      fEveColClusters->GetBin(ib)->Reset();
    }
  }

}

Int_t AliHLTEveTPC::ProcessClusters( AliHLTHOMERBlockDesc * block, TEvePointSet * cont, TEvePointSetArray * contCol ){
  //See header file for documentation


  // if(!fHistCharge) fHistCharge = new TH1F("ClusterCharge","ClusterCharge",100,0,500);
  // if(!fHistQMax) fHistQMax = new TH1F("QMax","QMax",50,0,250);
  // if(!fHistQMaxOverCharge) fHistQMaxOverCharge = new TH1F("QMaxOverCharge","QMaxOverCharge",50,0,1);


  Int_t   slice = block->GetSubDetector();
  Float_t phi   = ( slice + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
  
  AliHLTTPCClusterData *cd = reinterpret_cast<AliHLTTPCClusterData*> (block->GetData());
  UChar_t *data            = reinterpret_cast<UChar_t*> (cd->fSpacePoints);

  if ( cd->fSpacePointCnt != 0 ) {
    for (UInt_t iter = 0; iter < cd->fSpacePointCnt; ++iter, data += sizeof(AliHLTTPCSpacePointData)) {
      AliHLTTPCSpacePointData *sp = reinterpret_cast<AliHLTTPCSpacePointData*> (data);
      cont->SetNextPoint(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ);
      if (contCol)
	contCol->Fill(cos*sp->fX - sin*sp->fY, sin*sp->fX + cos*sp->fY, sp->fZ, sp->fCharge);


      // fHistCharge->Fill(sp->fCharge);
      // fHistQMax->Fill(sp->fQMax);
      // fHistQMaxOverCharge->Fill(((Float_t)sp->fQMax)/((Float_t)sp->fCharge));
    }
  }


  cont->ElementChanged();

  if(contCol) contCol->ElementChanged();
    
  return 0;  


}

void AliHLTEveTPC::DrawHistograms() {
  //See header file for documentation
  if (!fCanvas) {
    fCanvas = CreateCanvas("TPC Cl QA", "TPC Cluster QA");
    fCanvas->Divide(2, 2);
  }
  
  Int_t icd = 1;
  fCanvas->cd(icd++);
  fHistCharge->Draw();
  fCanvas->cd(icd++);
  fHistQMax->Draw();
  fCanvas->cd(icd++);
  fHistQMaxOverCharge->Draw();
  fCanvas->cd();
  

}
