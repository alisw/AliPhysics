
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
//         AliSpectraAODTrackCuts class
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
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODTrackCuts)


AliSpectraAODTrackCuts::AliSpectraAODTrackCuts(const char *name) : TNamed(name, "AOD Track Cuts"), fIsSelected(0), fTrackBits(0), fEtaCut(0), fDCACut(0), fPCut(0), fPtCut(0), fHistoCuts(0), fTrack(0)

{
   // Constructor
   fHistoCuts = new TH1I("fTrkCuts", "Track Cuts", kNTrkCuts, -0.5, kNTrkCuts - 0.5);
   fEtaCut = 100000.0; // default value of eta cut ~ no cut
   fDCACut = 100000.0; // default value of dca cut ~ no cut
   fPCut = 100000.0; // default value of p cut ~ no cut
   fPtCut = 100000.0; // default value of pt cut ~ no cut 
   
}

//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::IsSelected(AliAODTrack * track)
{
// Returns true if Track Cuts are selected and applied
   if (!track)
   {
      printf("ERROR: Could not receive track");
      return kFALSE;
   }
   fTrack = track;
   fIsSelected = (CheckTrackType() && CheckEtaCut() && CheckDCACut() && CheckPCut() && CheckPtCut());
   return fIsSelected ;
}
//_________________________________________________________

Bool_t AliSpectraAODTrackCuts::CheckTrackType()
{
   // Check track cuts
   if (fTrack->TestFilterBit(fTrackBits)) return kTRUE;
   fHistoCuts->Fill(kTrkBit);
   return kFALSE;
}
//________________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckEtaCut()
{
   // Check eta cut
   if (fTrack->Eta() < fEtaCut && fTrack->Eta() > - fEtaCut) return kTRUE;
   fHistoCuts->Fill(kTrkEta);
   return kFALSE;
}
//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckDCACut()
{
   // Check DCA cut
   if (fTrack->DCA() < fDCACut) return kTRUE;
   fHistoCuts->Fill(kTrkDCA);
   return kFALSE;
}
//________________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckPCut()
{
   // Check P cut
   if (fTrack->P() < fPCut) return kTRUE;
   fHistoCuts->Fill(kTrkP);
   return kFALSE;
}
//_______________________________________________________
Bool_t AliSpectraAODTrackCuts::CheckPtCut()
{
    // check Pt cut
//    if ((fTrack->Pt() < fPtCut) && (fTrack->Pt() > 0.3 )) return kTRUE;
   if (fTrack->Pt() < fPtCut) return kTRUE;
    fHistoCuts->Fill(kTrkPt);
    return kFALSE;
}
//_______________________________________________________
void AliSpectraAODTrackCuts::PrintCuts() const
{
  // Print cuts
    cout << "Track Cuts" << endl;
    cout << " > TrackBit\t" << fTrackBits << endl;
    cout << " > Eta cut\t" << fEtaCut << endl;
    cout << " > DCA cut\t" << fDCACut << endl;
    cout << " > P cut\t" << fPCut << endl;
    cout << " > Pt cut \t" << fPtCut << endl;
}
//_______________________________________________________
void AliSpectraAODTrackCuts::SetTrackType(UInt_t bit)
{
   // Set the type of track to be used. The argument should be the bit number. The mask is produced automatically.
   fTrackBits = (0x1 << (bit - 1));
}
//_______________________________________________________

Long64_t AliSpectraAODTrackCuts::Merge(TCollection* list)
{
  // Merge a list of AliSpectraAODTrackCuts objects with this.
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
    AliSpectraAODTrackCuts* entry = dynamic_cast<AliSpectraAODTrackCuts*> (obj);
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

