
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------
//         AliSpectraAODHistoManager class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisDataContainer.h"
#include "AliSpectraAODHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODHistoManager)


using namespace AliSpectraNameSpace;

AliSpectraAODHistoManager::AliSpectraAODHistoManager(const char *name): TNamed(name, "AOD Spectra Histo Manager"), fOutputList(0)
{
  // ctor
   fOutputList = new TList;
   fOutputList->SetOwner();
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   for (Int_t ihist  = 0; ihist < kNHist ; ihist++)
   {
      if (ihist <= kNPtHist) BookPtHistogram(kHistName[ihist]);  // PT histos
      if (ihist > kNPtHist) BookPIDHistogram(kHistName[ihist]);  // PID histos
   }
   TH1::AddDirectory(oldStatus);

}
//_______________________________________________________
TH1F* AliSpectraAODHistoManager::BookPtHistogram(const char * name)
{
   // Return a pt histogram with predefined binning, set the ID and add it to the output list
   AliInfo(Form("Booking pt histogram %s", name));

   TH1F * hist = new TH1F(name, Form("P_{T} distribution (%s)", name), 40, 0, 1.2);
   hist->GetXaxis()->SetTitle("P_{T} (GeV / c)");
   hist->GetYaxis()->SetTitle("No. of events");
   hist->SetMarkerStyle(kFullCircle);
   hist->Sumw2();
   fOutputList->Add(hist);

   return hist;
}
//_____________________________________________________________________________

TH2F* AliSpectraAODHistoManager::BookPIDHistogram(const char * name)
{
   // Return a pt histogram with predefined binning, set the ID and add it to the output list
   AliInfo(Form("Booking pt histogram %s", name));

   TH2F * hist = new TH2F(name, Form("Particle Identification (%s)", name), 100, 0, 1.2, 1000, 0, 1000);
   hist->GetXaxis()->SetTitle("P (GeV / c)");
   hist->GetYaxis()->SetTitle("TPC");
//  hist->Sumw2();
   fOutputList->Add(hist);

   return hist;
}


Long64_t AliSpectraAODHistoManager::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODHistoManager objects with this.
  // Returns the number of merged objects (including this).

  //  AliInfo("Merging");

  if (!list)
    return 0;

  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;

  // collections of all histograms
  TList collections;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraAODHistoManager* entry = dynamic_cast<AliSpectraAODHistoManager*> (obj);
    if (entry == 0) 
      continue;

    TList * hlist = entry->GetOutputList();      
    collections.Add(hlist);
    count++;
  }
  
  fOutputList->Merge(&collections);
  
  delete iter;

  return count+1;
}

