
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

AliSpectraAODEventCuts::AliSpectraAODEventCuts(const char *name) : TNamed(name, "AOD Event Cuts"), fAOD(0),fIsMC(0), fIsSelected(0), fCentralityCutMin(0), fCentralityCutMax(0), fHistoCuts(0),fHistoVtxBefSel(0),fHistoVtxAftSel(0),fHistoEtaBefSel(0),fHistoEtaAftSel(0),fHistoNChAftSel(0)
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
  Double_t cent=0;
  if(fIsMC)cent=fAOD->GetCentrality()->GetCentralityPercentile("V0M");
  else cent=ApplyCentralityPatchAOD049();
  
  if ( (cent <= fCentralityCutMax)  &&  (cent >= fCentralityCutMin) )  return kTRUE;   
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

Double_t AliSpectraAODEventCuts::ApplyCentralityPatchAOD049()
{
   //
   //Apply centrality patch for AOD049 outliers
   //
  // if (fCentralityType!="V0M") {
  //   AliWarning("Requested patch forAOD049 for wrong value (not centrality from V0).");
  //   return -999.0;
  // }
  AliCentrality *centrality = fAOD->GetCentrality();
   if (!centrality) {
     AliWarning("Cannot get centrality from AOD event.");
     return -999.0;
   }
   
   Float_t cent = (Float_t)(centrality->GetCentralityPercentile("V0M"));
   /*
     Bool_t isSelRun = kFALSE;
     Int_t selRun[5] = {138364, 138826, 138828, 138836, 138871};
     if(cent<0){
     Int_t quality = centrality->GetQuality();
     if(quality<=1){
     cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
     } else {
     Int_t runnum=aodEvent->GetRunNumber();
     for(Int_t ir=0;ir<5;ir++){
     if(runnum==selRun[ir]){
     isSelRun=kTRUE;
     break;
     }
     }
     if((quality==8||quality==9)&&isSelRun) cent=(Float_t)centrality->GetCentralityPercentileUnchecked("V0M");
     }
     }
   */
   if(cent>=0.0) {
     Float_t v0 = 0.0;
     AliAODEvent *aodEvent = (AliAODEvent *)fAOD;
     AliAODVZERO *aodV0 = (AliAODVZERO *) aodEvent->GetVZEROData();
     v0+=aodV0->GetMTotV0A();
     v0+=aodV0->GetMTotV0C();
     if ( (cent==0) && (v0<19500) ) {
       AliDebug(3, Form("Filtering issue in centrality -> cent = %5.2f",cent));
       return -999.0;
     }
     Float_t tkl = (Float_t)(aodEvent->GetTracklets()->GetNumberOfTracklets());
     Float_t val = 1.30552 +  0.147931 * v0;
     
     Float_t tklSigma[101] = {176.644, 156.401, 153.789, 153.015, 142.476, 137.951, 136.127, 129.852, 127.436, 124.86,
			      120.788, 115.611, 113.172, 110.496, 109.127, 104.421, 102.479, 99.9766, 97.5152, 94.0654,
			      92.4602, 89.3364, 87.1342, 83.3497, 82.6216, 81.1084, 78.0793, 76.1234, 72.9434, 72.1334,
			      68.0056, 68.2755, 66.0376, 62.9666, 62.4274, 59.65, 58.3776, 56.6361, 54.5184, 53.4224,
			      51.932, 50.8922, 48.2848, 47.912, 46.5717, 43.4114, 43.2083, 41.3065, 40.1863, 38.5255,
			      37.2851, 37.5396, 34.4949, 33.8366, 31.8043, 31.7412, 30.8392, 30.0274, 28.8793, 27.6398,
			      26.6488, 25.0183, 25.1489, 24.4185, 22.9107, 21.2002, 21.6977, 20.1242, 20.4963, 19.0235,
			      19.298, 17.4103, 16.868, 15.2939, 15.2939, 16.0295, 14.186, 14.186, 15.2173, 12.9504, 12.9504,
			      12.9504, 15.264, 12.3674, 12.3674, 12.3674, 12.3674, 12.3674, 18.3811, 13.7544, 13.7544,
			      13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544, 13.7544
     };
     
     if ( TMath::Abs(tkl-val) > 6.*tklSigma[(Int_t)cent] )  {
       AliDebug(3, Form("Outlier event in centrality -> cent = %5.2f",cent));
       return -999.0;
     }
   } else {
     //force it to be -999. whatever the negative value was
     cent = -999.;
   }
   return cent;
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
