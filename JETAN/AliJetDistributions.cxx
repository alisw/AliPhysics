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
 
//---------------------------------------------------------------------
// JetDistributions class 
// Get different basic distributions
// Authors: mercedes.lopez.noriega@cern.ch
//---------------------------------------------------------------------
 
#include "AliJetDistributions.h"
ClassImp(AliJetDistributions)
 
////////////////////////////////////////////////////////////////////////
// includes
#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
  
#include "AliJetProductionDataPDC2004.h"
#include "AliJet.h"
#include "AliJetReaderHeader.h"
#include "AliLeading.h"
  
////////////////////////////////////////////////////////////////////////
// constructor/destructor
  
AliJetDistributions::AliJetDistributions():
  fReaderHeader(0x0),
  fDirectory(0x0),
  fFile(0x0),
  fEventMin(0),
  fEventMax(-1),
  fRunMin(0),
  fRunMax(11),
  fPercentage(-1.0),
  fPartPtCut(0.0),
  fPythia(kFALSE),
  fDoPart(kTRUE),
  fDoGenJ(kTRUE),
  fDoRecJ(kTRUE),
  fPart(0),
  fGenJ(0),
  fRecJ(0),
  fRetaH(0),
  fRphiH(0),
  fRptH(0),
  fRetaphiH(0),
  fMultH(0)
{
  // Default constructor
  fFile = "jets.root";   
  SetReaderHeader();
  
}

////////////////////////////////////////////////////////////////////////
// define histogrames 

void AliJetDistributions::DefineHistograms()
{
  // Define histograms to be used
  fRetaH = new TH1F("fRetaH","Reconstructed eta",140,-1.2,1.2);
  SetProperties(fRetaH,"#eta","entries");
  fRphiH = new TH1F("fRphiH","Reconstructed phi",18,0,2.0*TMath::Pi());
  SetProperties(fRphiH,"#phi","entries");
  fRptH = new TH1F("fRptH","Reconstructed pt",150,0,150);
  SetProperties(fRptH,"p_{T} (GeV/c)","entries");

  fRetaphiH = new TH2F("fRetaphiH","Reconstructed eta vs. phi",140,-1.2,1.2,18,0.,2.0*TMath::Pi());
  SetProperties(fRetaphiH,"#eta","#phi");

  fMultH = new TH1F("fMultH","Reconstructed Multiplicity",1000,0,10000);
  SetProperties(fMultH,"Multiplicity","entries");
}

void AliJetDistributions::SetProperties(TH1* h,const char* x, const char* y) const
{
  // Properties of histograms (x title, y title and error propagation)
  h->SetXTitle(x);
  h->SetYTitle(y);
  h->Sumw2();
}
 
////////////////////////////////////////////////////////////////////////
// fill histogrames 

void AliJetDistributions::FillHistograms()
{
  // Run data 
  AliJetProductionDataPDC2004* runData = new AliJetProductionDataPDC2004();
  
  // Loop over runs
  TFile* jFile = 0x0;
  for (Int_t iRun = fRunMin; iRun <= fRunMax; iRun++) {
    char fn[20];
    sprintf(fn,"%s/%s.root",fDirectory,(runData->GetRunTitle(iRun)).Data());
    jFile = new TFile(fn);
    printf("  Analyzing run: %d %s\n", iRun,fn);	
    
    // Get reader header and events to be looped over
    AliJetReaderHeader *jReaderH = (AliJetReaderHeader*)(jFile->Get(fReaderHeader));
    if (fEventMin == -1) fEventMin =  jReaderH->GetFirstEvent();
    if (fEventMax == -1) {
      fEventMax =  jReaderH->GetLastEvent();
    } else {
      fEventMax = TMath::Min(fEventMax, jReaderH->GetLastEvent());
    }

    //AliUA1JetHeader *jH = (AliUA1JetHeader *) (jFile->Get("AliUA1JetHeader"));
    //jH->PrintParameters();

    
    // Loop over events
    for (Int_t i = fEventMin; i < fEventMax; i++) {
      //printf("  Analyzing run: %d  Event %d / %d \n",iRun, i, fEventMax);

      // Get next tree with AliJet
      char nameT[100];
      sprintf(nameT, "TreeJ%d",i);
      TTree *jetT =(TTree *)(jFile->Get(nameT));
      if (fDoRecJ) jetT->SetBranchAddress("FoundJet",    &fRecJ);
      if (fDoGenJ) jetT->SetBranchAddress("GenJet",      &fGenJ);
      if (fDoPart) jetT->SetBranchAddress("LeadingPart", &fPart);
      
      jetT->GetEntry(0);

      FillDistributions(fRecJ);

      delete jetT;
    } // end loop over events in one file
    if (jFile) jFile->Close();
    delete jFile;
  } // end loop over files
}

void AliJetDistributions::FillDistributions(AliJet *j)
{
  // Fill histrograms
  TArrayI inJet = j->GetInJet();
  TArrayF etain = j->GetEtaIn();
  TArrayF ptin = j->GetPtIn();
  TArrayF phiin = j->GetPhiIn();

  fMultH->Fill(inJet.GetSize(),1);
  for(Int_t i=0;i<inJet.GetSize();i++) {
    fRetaH->Fill(etain[i],1);
    fRphiH->Fill(phiin[i],1);
    fRptH->Fill(ptin[i],1);
    fRetaphiH->Fill(etain[i],phiin[i],1);
  }
}

////////////////////////////////////////////////////////////////////////
// Plot histogrames 

void AliJetDistributions::PlotHistograms()
{
  // improved!
  fMultH->Draw();
  fRetaH->Draw();
  fRphiH->Draw();
  fRetaphiH->Draw();
  fRptH->Draw();

}


////////////////////////////////////////////////////////////////////////
// Save histogrames 

void AliJetDistributions::SaveHistograms()
{
  // Save histogramas
  TFile *fOut = new TFile(fFile,"recreate");
  fOut->cd();
  fMultH->Write();
  fRetaphiH->Write();
  fRetaH->Write();
  fRphiH->Write();
  fRptH->Write();
  fOut->Close();
}

// main Analysis function

void AliJetDistributions::Analyze()
{
  // Do the analysis
  DefineHistograms();
  FillHistograms();
  PlotHistograms();
  SaveHistograms();
}


