
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
#include "AliSpectraAODTrackCuts.h"
#include "AliSpectraAODHistoManager.h"
#include <iostream>

using namespace std;

ClassImp(AliSpectraAODEventCuts)

AliSpectraAODEventCuts::AliSpectraAODEventCuts(const char *name) : TNamed(name, "AOD Event Cuts"), fAOD(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0)
{
  // Constructor
  fHistoCuts = new TH1I("fEventCuts", "Event Cuts", kNVtxCuts, -0.5, kNVtxCuts - 0.5);
  fHistoVtxBefSel = new TH1F("fHistoVtxBefSel", "Vtx distr before event selection",500,-15,15);
  fHistoVtxAftSel = new TH1F("fHistoVtxAftSel", "Vtx distr after event selection",500,-15,15);
  fHistoEtaBefSel = new TH1F("fHistoEtaBefSel", "Eta distr before event selection",500,-2,2);
  fHistoEtaAftSel = new TH1F("fHistoEtaAftSel", "Eta distr after event selection",500,-2,2);
  fHistoNChAftSel = new TH1F("fHistoNChAftSel", "NCh distr after event selection",3000,-0.5,2999.5);
  fCentralityCutMin = 0.0;      // default value of centrality cut minimum, 0 ~ no cut
  fCentralityCutMax = 10000.0;  // default value of centrality cut maximum,  ~ no cut

}

//______________________________________________________
Bool_t AliSpectraAODEventCuts::IsSelected(AliAODEvent * aod,AliSpectraAODTrackCuts     *trackcuts)
{
  // Returns true if Event Cuts are selected and applied
  fAOD = aod;
  fTrackCuts = trackcuts;
  fHistoCuts->Fill(kProcessedEvents);
  //loop on tracks, before event selection, filling QA histos
  AliAODVertex * vertex = fAOD->GetPrimaryVertex();//FIXME vertex is recreated
  if(vertex)fHistoVtxBefSel->Fill(vertex->GetZ());
  fIsSelected = (CheckVtxRange() && CheckCentralityCut());
  if(fIsSelected){
    fHistoCuts->Fill(kAcceptedEvents);
    if(vertex)fHistoVtxAftSel->Fill(vertex->GetZ());
  }
  Int_t Nch=0;
  for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
    AliAODTrack* track = fAOD->GetTrack(iTracks);
    if (!fTrackCuts->IsSelected(track)) continue;
    fHistoEtaBefSel->Fill(track->Eta());
    if(fIsSelected){
      fHistoEtaAftSel->Fill(track->Eta());
      Nch++;
    }
  }
  if(fIsSelected)fHistoNChAftSel->Fill(Nch);
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
  cout << " > Number of processed events: " << fHistoCuts->GetBinContent(kProcessedEvents + 1) << endl;
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
  TList collections;//FIXME we should use only 1 collection
  TList collections_histoVtxBefSel;
  TList collections_histoVtxAftSel;
  TList collections_histoEtaBefSel;
  TList collections_histoEtaAftSel;
  TList collections_histoNChAftSel;

  Int_t count = 0;

  while ((obj = iter->Next())) {
    AliSpectraAODEventCuts* entry = dynamic_cast<AliSpectraAODEventCuts*> (obj);
    if (entry == 0) 
      continue;

    TH1I * histo = entry->GetHistoCuts();      
    collections.Add(histo);
    TH1F * histo_histoVtxBefSel = entry->GetHistoVtxBefSel();      
    collections_histoVtxBefSel.Add(histo_histoVtxBefSel);
    TH1F * histo_histoVtxAftSel = entry->GetHistoVtxAftSel();      
    collections_histoVtxAftSel.Add(histo_histoVtxAftSel);
    TH1F * histo_histoEtaBefSel = entry->GetHistoEtaBefSel();      
    collections_histoEtaBefSel.Add(histo_histoEtaBefSel);
    TH1F * histo_histoEtaAftSel = entry->GetHistoEtaAftSel();      
    collections_histoEtaAftSel.Add(histo_histoEtaAftSel);
    TH1F * histo_histoNChAftSel = entry->GetHistoNChAftSel();      
    collections_histoNChAftSel.Add(histo_histoNChAftSel);
    count++;
  }
  
  fHistoCuts->Merge(&collections);
  fHistoVtxBefSel->Merge(&collections_histoVtxBefSel);
  fHistoVtxAftSel->Merge(&collections_histoVtxAftSel);
  fHistoEtaBefSel->Merge(&collections_histoEtaBefSel);
  fHistoEtaAftSel->Merge(&collections_histoEtaAftSel);
  fHistoNChAftSel->Merge(&collections_histoNChAftSel);
  
  delete iter;

  return count+1;
}

/// FIXME: Q vector
// //Selection on QVector, before ANY other selection on the event
// //Spectra MUST be normalized wrt events AFTER the selection on Qvector
// // Can we include this in fEventCuts
// Double_t Qx2EtaPos = 0, Qy2EtaPos = 0;
// Double_t Qx2EtaNeg = 0, Qy2EtaNeg = 0;
// Int_t multPos = 0;
// Int_t multNeg = 0;
// for(Int_t iT = 0; iT < fAOD->GetNumberOfTracks(); iT++) {
//   AliAODTrack* aodTrack = fAOD->GetTrack(iT);
//   if (!fTrackCuts->IsSelected(aodTrack)) continue;
//   if (aodTrack->Eta() >= 0){
//     multPos++;
//     Qx2EtaPos += TMath::Cos(2*aodTrack->Phi()); 
//     Qy2EtaPos += TMath::Sin(2*aodTrack->Phi());
//   } else {
//     multNeg++;
//     Qx2EtaNeg += TMath::Cos(2*aodTrack->Phi()); 
//     Qy2EtaNeg += TMath::Sin(2*aodTrack->Phi());
//   }
// } 
// Double_t qPos=-999;
// if(multPos!=0)qPos= TMath::Sqrt((Qx2EtaPos*Qx2EtaPos + Qy2EtaPos*Qy2EtaPos)/multPos);
// Double_t qNeg=-999;
// if(multNeg!=0)qNeg= TMath::Sqrt((Qx2EtaNeg*Qx2EtaNeg + Qy2EtaNeg*Qy2EtaNeg)/multNeg);
  
// if((qPos>fTrackCuts->GetQvecMin() && qPos<fTrackCuts->GetQvecMax()) || (qNeg>fTrackCuts->GetQvecMin() && qNeg<fTrackCuts->GetQvecMax())){

//fill q distributions vs centrality, after all event selection
// fHistMan->GetqVecHistogram(kHistqVecPos)->Fill(qPos,fAOD->GetCentrality()->GetCentralityPercentile("V0M"));  // qVector distribution
// fHistMan->GetqVecHistogram(kHistqVecNeg)->Fill(qNeg,fAOD->GetCentrality()->GetCentralityPercentile("V0M"));  // qVector distribution
