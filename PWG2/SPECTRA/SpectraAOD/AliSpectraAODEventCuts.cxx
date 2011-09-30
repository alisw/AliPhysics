
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
//         AliSpectraAODEventCuts class
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TLegend.h"
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
#include "AliSpectraAODEventCuts.h"
#include "AliSpectraAODHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODEventCuts)

AliSpectraAODEventCuts::AliSpectraAODEventCuts(const char *name) : TNamed(name, "AOD Event Cuts"), fAOD(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fHistoCuts(0)
{
   // Constructor
   fHistoCuts = new TH1I("fEventCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
   fCentralityCutMin = 0.0;      // default value of centrality cut minimum, 0 ~ no cut
   fCentralityCutMax = 10000.0;  // default value of centrality cut maximum,  ~ no cut

}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::IsSelected(AliAODEvent * aod)
{
// Returns true if Event Cuts are selected and applied
   fAOD = aod;
   fIsSelected = (CheckVtxRange() && CheckCentralityCut());
   if(fIsSelected)  fHistoCuts->Fill(kAcceptedEvents);
   return fIsSelected;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckVtxRange()
{
  // reject events outside of range
   AliAODVertex * vertex = fAOD->GetPrimaryVertex();
   if (!vertex)
   {
      fHistoCuts->Fill(kVtxNoEvent);
      return kFALSE;
   }
   if (TMath::Abs(vertex->GetZ()) < 10)
   {
      return kTRUE;
   }
   fHistoCuts->Fill(kVtxRange);
   return kFALSE;
}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::CheckCentralityCut()
{
   // Check centrality cut
   if ( (fAOD->GetCentrality()->GetCentralityPercentile("V0M") <= fCentralityCutMax)  &&  (fAOD->GetCentrality()->GetCentralityPercentile("V0M") >= fCentralityCutMin) )  return kTRUE;   
   fHistoCuts->Fill(kVtxCentral);
   return kFALSE;
}

//______________________________________________________
void AliSpectraAODEventCuts::PrintCuts()
{
    // print info about event cuts
    cout << "Event Stats" << endl;
    cout << " > Number of accepted events: " << fHistoCuts->GetBinContent(kAcceptedEvents + 1) << endl;
    cout << " > Vertex out of range: " << fHistoCuts->GetBinContent(kVtxRange + 1) << endl;
    cout << " > Events cut by centrality: " << fHistoCuts->GetBinContent(kVtxCentral + 1) << endl;
    cout << " > Events without vertex: " << fHistoCuts->GetBinContent(kVtxNoEvent + 1) << endl;
    }
//______________________________________________________

Long64_t AliSpectraAODEventCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODEventCuts objects with this.
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
    AliSpectraAODEventCuts* entry = dynamic_cast<AliSpectraAODEventCuts*> (obj);
    if (entry == 0) 
      continue;

    TH1I * histo = entry->GetHistoCuts();      
    collections.Add(histo);
    count++;
  }
  
  fHistoCuts->Merge(&collections);
  
  delete iter;

  return count+1;
}

