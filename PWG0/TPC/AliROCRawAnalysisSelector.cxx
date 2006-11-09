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

/* $Id$ */


#include "AliROCRawAnalysisSelector.h"

#include <AliLog.h>

#include <AliRawEvent.h>
#include <AliRawReaderRoot.h>
#include <AliRawEventHeaderBase.h>
#include <AliTPCRawStream.h>
#include <AliTPCParamSR.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TTimeStamp.h>


ClassImp(AliROCRawAnalysisSelector)

AliROCRawAnalysisSelector::AliROCRawAnalysisSelector() :
  TSelector(),
  fRawEvent(0),
  fTree(0)
{
  
  //
  // Constructor. Initialization of pointers
  //
  
  AliDebug(AliLog::kInfo, "Constructor....");
}

AliROCRawAnalysisSelector::~AliROCRawAnalysisSelector()
{
  //
  // Destructor
  //
}

void AliROCRawAnalysisSelector::SlaveBegin(TTree* tree)
{
  //
  
  AliDebug(AliLog::kInfo, "SlaveBegin....");

  TSelector::SlaveBegin(tree);

  if (tree != 0)
    Init(tree);
} 

void AliROCRawAnalysisSelector::Init(TTree* tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses of the tree
  // will be set. It is normaly not necessary to make changes to the
  // generated code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running with PROOF.

  fTree = tree;

  TSelector::Init(tree);

  // Set branch address
  if (tree) {
    AliDebug(AliLog::kInfo, "INFO: Tree found");

    tree->SetBranchAddress("rawevent", &fRawEvent);
  }

}

Bool_t AliROCRawAnalysisSelector::Process(Long64_t entry)
{
  //
  // Implement your analysis here. Do not forget to call the parent class Process by
  // if (TSelector::Process(entry) == kFALSE)
  //   return kFALSE;
  //

  fTree->GetTree()->GetEntry(entry);
  
  AliRawReaderRoot* rawReader        = new AliRawReaderRoot(fRawEvent);
  AliRawEventHeaderBase* eventHeader = (AliRawEventHeaderBase*)rawReader->GetEventHeader();


  if (eventHeader) {
    eventHeader->Print();
    
    UInt_t timeStamp = eventHeader->Get("Timestamp");
    UInt_t eventType = eventHeader->Get("Type");
    
    printf("Time stamp: %s, event type %d\n", TTimeStamp(timeStamp).AsString(), eventType);
  }           
  
  AliTPCRawStream* tpcRawStream = new AliTPCRawStream(rawReader);
     
  AliTPCParamSR* fParam = new AliTPCParamSR;
     
  const Int_t kNIS = fParam->GetNInnerSector();
  const Int_t kNOS = fParam->GetNOuterSector();
  const Int_t kNS = kNIS + kNOS;
  
  Float_t fSign;
    
  // kNS
  for(Int_t fSector = 0; fSector < kNS; fSector++) {
    printf("*** Looking at sector %d ***\n", fSector);
            
    Int_t nRows = 0;
    Int_t nDDLs = 0, indexDDL = 0;
    
    if (fSector < kNIS) {
      nRows = fParam->GetNRowLow();
      fSign = (fSector < kNIS/2) ? 1 : -1;
      nDDLs = 2;
      indexDDL = fSector * 2;
    }
    else {
      nRows = fParam->GetNRowUp();
      fSign = ((fSector-kNIS) < kNOS/2) ? 1 : -1;
      nDDLs = 4;
      indexDDL = (fSector-kNIS) * 4 + kNIS * 2;
    }
    
    // Loas the raw data for corresponding DDLs
    rawReader->Reset();
    tpcRawStream->SetOldRCUFormat(kTRUE);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);
    Int_t digCounter=0;
    // Begin loop over altro data
    
    printf("Selected DDLs %d ... %d\n", indexDDL,indexDDL+nDDLs-1);
    
//     Int_t count = 0;
    
//     // hist for this sector
//     TString title;
//     title.Form("sector_%d", fSector);
//     gROOT->cd();
//     TH3F* hist = dynamic_cast<TH3F*> (gROOT->FindObject(title)); 
//     if (!hist)
//       hist = new TH3F(title, Form("%s;row;pad;time", title.Data()), 90, -0.5, 89.5, 120, -0.5, 119.5, 100, 0, 1200);
    
//     title.Form("sector_%d_signal", fSector);
//     TH1F* signal = dynamic_cast<TH1F*> (gROOT->FindObject(title)); 
//     if (!signal)
//       signal = new TH1F(title, title, 200, 0, 2000);
//     //TProfile3D* hist = new TProfile3D(title, title, 90, -0.5, 89.5, 200, -0.5, 199.5, 100, 0, 1500); 
    
//     while (tpcRawStream->Next())
//       {
// 	if (tpcRawStream->GetSector() != fSector)
// 	  {
// 	    printf("Sector index mismatch ! Expected (%d), but got (%d) !\n",fSector,tpcRawStream->GetSector());
// 	    return;
// 	  }
	
	
// 	if ((count++ % 100000) == 0)
// 	  printf("Found %d. digit in sector %d: row %d, pad %d, time %d, signal %d\n", count, tpcRawStream->GetSector(), tpcRawStream->GetRow(), tpcRawStream->GetPad(), tpcRawStream->GetTime(), tpcRawStream->GetSignal());
	
// 	if (tpcRawStream->GetSignal() > 200)
// 	  hist->Fill(tpcRawStream->GetRow(), tpcRawStream->GetPad(), tpcRawStream->GetTime(), tpcRawStream->GetSignal());
// 	signal->Fill(tpcRawStream->GetSignal());
// 	//if (++count == 2)
// 	//    break;
//       }
    
//     if (count > 0 && event == Nevents-1)
//       {
// 	TCanvas* canvas = new TCanvas(hist->GetName(), hist->GetName(), 900, 450);
// 	canvas->Divide(2, 1);
// 	canvas->cd(1);
// 	hist->Draw();
// 	canvas->cd(2);
// 	signal->Draw();
//       }    
//     //else
//     //    delete hist;
    
  }
   
  return kTRUE;
}

void AliROCRawAnalysisSelector::SlaveTerminate()
{
  //
  AliDebug(AliLog::kInfo, "SlaveTerminate....");
  
  //for (Int_t i=0; i<kTPCSectors; i++)
  // if (fClusterHistograms[i])
  //    fOutput->Add(fClusterHistograms[i]);
} 

void AliROCRawAnalysisSelector::Terminate()
{
  AliDebug(AliLog::kInfo, "Terminate....");

  // TODO read from output list for PROOF
    
  TFile* file = TFile::Open("rocRaw.root", "RECREATE");
  
  //  for (Int_t i=0; i<kTPCSectors; i++)
  //  if (fClusterHistograms[i])
  //    fClusterHistograms[i]->SaveHistograms();

  file->Close();
} 
