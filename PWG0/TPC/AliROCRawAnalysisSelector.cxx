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
#include <TROOT.h>

#include <TPC/AliTPCRawHistograms.h>


ClassImp(AliROCRawAnalysisSelector)

AliROCRawAnalysisSelector::AliROCRawAnalysisSelector() :
  TSelector(),
  fRawEvent(0),
  fTree(0),
  fParam(0)  
{
  //
  // Constructor. Initialization of pointers
  //

  for (Int_t i=0; i<kTPCSectors; i++)
    fHistograms[i] = 0;

  fParam = new AliTPCParamSR;
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

  // Set branch address
  if (tree) 
  {
    AliDebug(AliLog::kInfo, "INFO: Tree found");

    tree->SetBranchAddress("rawevent", &fRawEvent);
  }
}

Bool_t AliROCRawAnalysisSelector::Process(Long64_t entry)
{
  //
  //
  //

  AliDebug(AliLog::kInfo, Form("Processing event %lld", entry));
  
  fTree->GetTree()->GetEntry(entry);

  if (!fRawEvent)
  {
    AliDebug(AliLog::kError, "fRawEvent empty");
    return kFALSE;
  }
  
  AliRawReaderRoot* rawReader = new AliRawReaderRoot(fRawEvent);

  const AliRawEventHeaderBase* eventHeader = dynamic_cast<const AliRawEventHeaderBase*> (rawReader->GetEventHeader());
  if (eventHeader) 
  {
    eventHeader->Print();
    
    UInt_t timeStamp = eventHeader->Get("Timestamp");
    UInt_t eventType = eventHeader->Get("Type");
    
    AliDebug(AliLog::kInfo, Form("Time stamp: %s, event type %d", TTimeStamp(timeStamp).AsString(), eventType));
  }           
  
  AliTPCRawStream* tpcRawStream = new AliTPCRawStream(rawReader);
     
  const Int_t kNIS = fParam->GetNInnerSector();
  const Int_t kNOS = fParam->GetNOuterSector();
  const Int_t kNS = kNIS + kNOS;
  
  for (Int_t sector = 0; sector < kNS; sector++) 
  {
    AliDebug(AliLog::kInfo, Form("*** Looking at sector %d ***", sector));
            
    Int_t nRows = 0;
    Int_t nDDLs = 0, indexDDL = 0;
    
    if (sector < kNIS) 
    {
      nRows = fParam->GetNRowLow();
      nDDLs = 2;
      indexDDL = sector * 2;
    }
    else 
    {
      nRows = fParam->GetNRowUp();
      nDDLs = 4;
      indexDDL = (sector-kNIS) * 4 + kNIS * 2;
    }
    
    // Loas the raw data for corresponding DDLs
    rawReader->Reset();
    tpcRawStream->SetOldRCUFormat(kTRUE);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);
    
    AliDebug(AliLog::kDebug, Form("Selected DDLs %d ... %d", indexDDL, indexDDL+nDDLs-1));
    
    Int_t count = 0;
    
    while (tpcRawStream->Next())
    {
   	  if (tpcRawStream->GetSector() != sector)
      {
 	    AliDebug(AliLog::kError, Form("Sector index mismatch ! Expected (%d), but got (%d) !",sector,tpcRawStream->GetSector()));
 	    return kFALSE;
      }
      
      if ((count++ % 100000) == 0)
 	    AliDebug(AliLog::kDebug, Form("Found %d. digit in sector %d: row %d, pad %d, time %d, signal %d", count, 
            tpcRawStream->GetSector(), tpcRawStream->GetRow(), tpcRawStream->GetPad(), tpcRawStream->GetTime(), tpcRawStream->GetSignal()));

      if (!fHistograms[sector])
      {
        // not sure if this is still needed, should prevent creation of the histogram in the opened root file
        gROOT->cd();
        fHistograms[sector] = new AliTPCRawHistograms(sector);
      }
      
      fHistograms[sector]->FillDigit(tpcRawStream);
    }
  }
  
  delete fRawEvent;
  fRawEvent = 0;
   
  return kTRUE;
}

void AliROCRawAnalysisSelector::SlaveTerminate()
{
  //
  for (Int_t i=0; i<kTPCSectors; i++)
   if (fHistograms[i])
     fOutput->Add(fHistograms[i]);
} 

void AliROCRawAnalysisSelector::Terminate()
{
  AliDebug(AliLog::kInfo, "Terminate....");

  // TODO read from output list for PROOF
    
  TFile* file = TFile::Open("rocRaw.root", "RECREATE");
  
  for (Int_t i=0; i<kTPCSectors; i++)
    if (fHistograms[i])
      fHistograms[i]->SaveHistograms();

  file->Close();

  for (Int_t i=0; i<kTPCSectors; i++)
    if (fHistograms[i])
      fHistograms[i]->DrawHistograms();
} 
