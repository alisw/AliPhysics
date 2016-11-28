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

#include "TClonesArray.h"
#include "AliHLTEveTRD.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "AliHLTEveBase.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"
#include "TEvePointSet.h"
#include "TColor.h"
#include "TMath.h"
#include "TH1F.h"
#include "AliTRDcluster.h"

ClassImp(AliHLTEveTRD)

AliHLTEveTRD::AliHLTEveTRD() : 
  AliHLTEveBase("TRD"), 
  fEveClusters(NULL),
  fEveColClusters(NULL),
  fNColorBins(15),
  fClusterArray(NULL)
{
  // Constructor.
  fClusterArray = new TClonesArray("AliTRDcluster");
}

AliHLTEveTRD::~AliHLTEveTRD()
{
  //Destructor, not implemented
  if(fEveColClusters)
    delete fEveColClusters;
  fEveColClusters = NULL;
  
  if(fEveClusters)
    delete fEveClusters;
  fEveClusters = NULL;

  fClusterArray->Delete();
  delete fClusterArray;
  fClusterArray = NULL;
}


void AliHLTEveTRD::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation

  if ( ! block->GetDataType().CompareTo("CLUSTERS") ) {
    
    if(!fEveColClusters){
      fEveColClusters = CreatePointSetArray();
      AddElement(fEveColClusters);
    } 
    
    ProcessClusters(block, fEveColClusters);
   
  } else if ( ! block->GetDataType().CompareTo("ROOTHIST") ) {
  
    if(!fCanvas) {
      fCanvas = CreateCanvas("TRD QA", "TRD QA");
      fCanvas->Divide(3, 3);
    }

    AddHistogramsToCanvas(block, fCanvas, fHistoCount);
		   
  }
  
}


void AliHLTEveTRD::AddHistogramsToCanvas(AliHLTHOMERBlockDesc* block, TCanvas * canvas, Int_t &cdCount ) {
  //See header file for documentation
  if ( ! block->GetClassName().CompareTo("TH1F")) {
    TH1F* histo = reinterpret_cast<TH1F*>(block->GetTObject());
    ++cdCount;
  
    TVirtualPad* pad = canvas->cd(cdCount);
    histo->Draw();
    pad->SetGridy();
    pad->SetGridx();

    if ( ! strcmp(histo->GetName(), "nscls") ) {
      histo->GetXaxis()->SetRangeUser(0.,15.);
    }

    if ( ! strcmp(histo->GetName(),"sclsdist") ||
	 ! strcmp(histo->GetName(),"evSize") )
      pad->SetLogy();
  }

}



TEvePointSet * AliHLTEveTRD::CreatePointSet() {
  //See header file for documentation
  TEvePointSet * ps = new TEvePointSet("TRD Clusters");
  ps->SetMainColor(kBlue);
  ps->SetMarkerStyle((Style_t)kFullDotSmall);
 
  return ps;

}

TEvePointSetArray * AliHLTEveTRD::CreatePointSetArray(){

  Int_t i = 0;

  cout << i++ << endl;
  //See header file for documentation
  TEvePointSetArray * cc = new TEvePointSetArray("TRD Clusters Colorized");
  cout << i++ << endl;
  cc->SetMainColor(kRed);
  cout << i++ << endl;
  cc->SetMarkerStyle(4); // antialiased circle
  cout << i++ << endl;
  cc->SetMarkerSize(0.4);
  cout << i++ << endl;
  cc->InitBins("Cluster Charge", fNColorBins, 0., fNColorBins*100.);
  cout << i++ << endl;
  
  const Int_t nCol = TColor::GetNumberOfColors();
  cout << i++ << endl;
  for (Int_t ii = 0; ii < fNColorBins + 1; ++ii) {
    cc->GetBin(ii)->SetMainColor(TColor::GetColorPalette(ii * nCol / (fNColorBins+2)));
  }
  cout << i++ << endl;

  return cc;
     
}


void AliHLTEveTRD::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  // if(fEveClusters) fEveClusters->ResetBBox();

  if(fEveColClusters) {
    for (Int_t ib = 0; ib <= fNColorBins + 1; ++ib) {
      fEveColClusters->GetBin(ib)->ResetBBox();
    }
  }
}

void AliHLTEveTRD::ResetElements(){
  //See header file for documentation
  // if(fEveClusters) fEveClusters->Reset();
 
  if(fEveColClusters){
    for (Int_t ib = 0; ib <= fNColorBins + 1; ++ib) {
      fEveColClusters->GetBin(ib)->Reset();
    }
  }

  fHistoCount = 0;
}

Int_t AliHLTEveTRD::ProcessClusters( AliHLTHOMERBlockDesc * block, TEvePointSetArray * contCol ){
  //See header file for documentation

  Int_t iResult = 0;

  Int_t sm = block->GetSubDetector();
  if ( sm == 6 ) sm = 7;
  
  Float_t phi   = ( sm + 0.5 ) * TMath::Pi() / 9.0;  
  Float_t cos   = TMath::Cos( phi );
  Float_t sin   = TMath::Sin( phi );
  
  Byte_t* ptrData = reinterpret_cast<Byte_t*>(block->GetData());
  UInt_t ptrSize = block->GetSize();
  Int_t unused;

  /*AliHLTTRDUtils::ReadClusters(fClusterArray, ptrData, ptrSize, &unused);

  for(int num=fClusterArray->GetEntriesFast(); num--;){
    AliTRDcluster* trdCluster = (AliTRDcluster*)fClusterArray->At(num);
    contCol->Fill(cos*trdCluster->GetX() - sin*trdCluster->GetY(), 
		  sin*trdCluster->GetX() + cos*trdCluster->GetY(), 
		  trdCluster->GetZ(),
		  trdCluster->GetQ() );  
  }*/
  
  return iResult;
}


