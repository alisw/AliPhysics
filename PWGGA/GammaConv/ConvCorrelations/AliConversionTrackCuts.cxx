/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Author: Svein Lindal <slindal@fys.uio.no>                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliConversionTrackCuts.cxx
/// @author Svein Lindal
/// @brief  Base class for analysation of conversion particle - track correlations


#include "AliConversionTrackCuts.h"
//#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include <TFormula.h>
#include <iostream>


using namespace std;
ClassImp(AliConversionTrackCuts)

//________________________________________________________________________
AliConversionTrackCuts::AliConversionTrackCuts() : 
AliAnalysisCuts(),
  fFlagsOn(0x0),
  fFlagsOff(0x0),
  fRejectKinkDaughters(kTRUE),
  fDCARfixed(kTRUE),
  fDCARptFormula(""),
  fDCARmax(1E20),
  fDCAZfixed(kTRUE),
  fDCAZptFormula(""),
  fDCAZmax(1E20),
  fSPDminNClusters(0),
  fITSminNClusters(0),
  fITSmaxChi2(1E20),
  fTPCminNClusters(0),
  fTPCmaxChi2(1E20),
  fAODTestFilterBit(-1)
{
  //Constructor
}
//________________________________________________________________________
AliConversionTrackCuts::AliConversionTrackCuts(TString name, TString title = "title") : 
  AliAnalysisCuts(name, title),
  fFlagsOn(0x0),
  fFlagsOff(0x0),
  fRejectKinkDaughters(kTRUE),
  fDCARfixed(kTRUE),
  fDCARptFormula(""),
  fDCARmax(1E20),
  fDCAZfixed(kTRUE),
  fDCAZptFormula(""),
  fDCAZmax(1E20),
  fSPDminNClusters(0),
  fITSminNClusters(0),
  fITSmaxChi2(1E20),
  fTPCminNClusters(0),
  fTPCmaxChi2(1E20),
  fAODTestFilterBit(-1)

{
  //Constructor
}


//________________________________________________________________________________
// AliConversionTrackCuts::~AliConversionTrackCuts() {
//   ///destructor
// }

Bool_t AliConversionTrackCuts::AcceptTrack(AliAODTrack * track, AliAODEvent * aodEvent) {
  // Check an AOD track.
// This is done doing directly all checks, since there is not
// an equivalend checker for AOD tracks
//

   // try to retrieve the reference AOD event
   // AliAODEvent *aodEvent = 0x0;
   // if (fEvent) aodEvent = fEvent->GetRefAOD();
   // if (!aodEvent) {
   //    AliError("AOD reference event is not initialized!");
   //    return kFALSE;
   // }


   // step #0: check SPD and ITS clusters
   Int_t nSPD = 0;
   nSPD  = TESTBIT(track->GetITSClusterMap(), 0);
   nSPD += TESTBIT(track->GetITSClusterMap(), 1);
   if (nSPD < fSPDminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Not enough SPD clusters in this track. Rejected");
      return kFALSE;
   }

   // step #1: check number of clusters in TPC
   if (track->GetTPCNcls() < fTPCminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Too few TPC clusters. Rejected");
      return kFALSE;
   }

   if (track->GetITSNcls() < fITSminNClusters) {
      AliDebug(AliLog::kDebug + 2, "Too few ITS clusters. Rejected");
      return kFALSE;
   }
 
   // step #2: check chi square
   if (track->Chi2perNDF() > fTPCmaxChi2) {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
   }

   if (track->Chi2perNDF() > fITSmaxChi2) {
      AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
      return kFALSE;
   }


   // step #3: reject kink daughters
   AliAODVertex *vertex = track->GetProdVertex();
   if (vertex && fRejectKinkDaughters) {
      if (vertex->GetType() == AliAODVertex::kKink) {
         AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
         return kFALSE;
      }
   }


   // step #4: DCA cut (transverse)
   Double_t b[2], cov[3];
   vertex = aodEvent->GetPrimaryVertex();
   if (!vertex) {
      AliDebug(AliLog::kDebug + 2, "NULL vertex");
      return kFALSE;
   }

   if (!track->PropagateToDCA(vertex, aodEvent->GetMagneticField(), kVeryBig, b, cov)) {
      AliDebug(AliLog::kDebug + 2, "Failed propagation to vertex");
      return kFALSE;
   }

   // if the DCA cut is not fixed, compute current value
   if (!fDCARfixed) {
      static TString str(fDCARptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaXY(Form("%s_dcaXY", GetName()), str.Data());
      fDCARmax = dcaXY.Eval(track->Pt());
   }

   // check the cut
   if (TMath::Abs(b[0]) > fDCARmax) {
      AliDebug(AliLog::kDebug + 2, "Too large transverse DCA");
      return kFALSE;
   }

 
   // step #5: DCA cut (longitudinal)
   // the DCA has already been computed above
   // if the DCA cut is not fixed, compute current value
   if (!fDCAZfixed) {
      static TString str(fDCAZptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaZ(Form("%s_dcaXY", GetName()), str.Data());
      fDCAZmax = dcaZ.Eval(track->Pt());
   }

   // check the cut
  if (TMath::Abs(b[1]) > fDCAZmax) {
      AliDebug(AliLog::kDebug + 2, "Too large longitudinal DCA");
      return kFALSE;
   }

 
   // step #6: check eta/pt range
   if (track->Eta() < fEta[0] || track->Eta() > fEta[1]) {
      AliDebug(AliLog::kDebug + 2, "Outside ETA acceptance");
      return kFALSE;
   }

   // if (track->Pt() < fPt[0] || track->Pt() > fPt[1]) {
   //    AliDebug(AliLog::kDebug + 2, "Outside PT acceptance");
   //    return kFALSE;
   // }

   // if we are here, all cuts were passed and no exit point was got
   return kTRUE;
}

//_________________________________________________________________________________________________
void AliConversionTrackCuts::Print(const Option_t *) const
{
//
// Print information on this cut
//

   AliInfo(Form("Cut name                : %s", GetName()));
   AliInfo(Form("Required flags (off, on): %lx %lx", fFlagsOn, fFlagsOff));
   AliInfo(Form("Ranges in eta, pt       : %.2f - %.2f, %.2f - %.2f", fEta[0], fEta[1], fPt[0], fPt[1]));
   AliInfo(Form("Kink daughters are      : %s", (fRejectKinkDaughters ? "rejected" : "accepted")));
   AliInfo(Form("TPC requirements        : min. cluster = %d, max chi2 = %f", fTPCminNClusters, fTPCmaxChi2));
   AliInfo(Form("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f", fITSminNClusters, fSPDminNClusters, fITSmaxChi2));

   if (fDCARfixed) {
	 AliInfo(Form("DCA r cut               : fixed to %f cm", fDCARmax));
   } else {
	 AliInfo(Form("DCA r cut formula       : %s", fDCARptFormula.Data()));
   }
   
   if (fDCAZfixed) {
	 AliInfo(Form("DCA z cut               : fixed to %f cm", fDCAZmax));
   } else {
	 AliInfo(Form("DCA z cut formula       : %s", fDCAZptFormula.Data()));
   }
   
   
}







 
//________________________________________________________________________________
// void AliConversionTrackCuts::CreateHistograms() {
//   CreateBaseHistograms();
// }


// ///________________________________________________________________________________
// void AliConversionTrackCuts::SetUpDefaultBins() {
//   //Set up default bins
//   Double_t ptbins[19] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.5, 15, 20, 25, 30, 50, 100};
//   fAxisdEta.Set(160, -1.6, 1.6);
//   fAxisdEta.SetNameTitle("dEta", "delta eta");

//   fAxisdPhi.Set(64, -TMath::PiOver2(), 3*TMath::PiOver2());
//   fAxisdPhi.SetNameTitle("dPhi", "delta Phi");

//   fAxistPt.Set(18, ptbins);
//   fAxistPt.SetNameTitle("tPt", "trigger Pt");

//   fAxiscPt.Set(18, ptbins);
//   fAxiscPt.SetNameTitle("cPt", "track Pt");

//   fAxisIso.Set(3, -0.5, 2.5);
//   fAxisIso.SetNameTitle("iso", "isolation");

//   fAxesList.AddAt(&fAxisdEta, 0);
//   fAxesList.AddAt(&fAxisdPhi, 1);
//   fAxesList.AddAt(&fAxistPt, 2);
//   fAxesList.AddAt(&fAxiscPt, 3);
//   fAxesList.AddAt(&fAxisIso, 4);

//   fAxisMEEta.Set(160, -0.8, 0.8);
//   fAxisMEEta.SetNameTitle("eta", "eta");
  
//   fAxisMEPhi.Set(64, 0, TMath::TwoPi());
//   fAxisMEPhi.SetNameTitle("phi", "phi");

//   fTrackAxisList.AddAt(&fAxisMEEta, 0);
//   fTrackAxisList.AddAt(&fAxisMEPhi, 1);
//   fTrackAxisList.AddAt(&fAxistPt, 2);
//   fTrackAxisList.AddAt(&fAxiscPt, 3);
//   fTrackAxisList.AddAt(&fAxisIso, 4);

//   fTrigAxisList.AddAt(&fAxisMEEta, 0);
//   fTrigAxisList.AddAt(&fAxisMEPhi, 1);
//   fTrigAxisList.AddAt(&fAxistPt, 2);
//   fTrigAxisList.AddAt(&fAxisIso, 3);

//   for(int iIso = 0; iIso < 2; iIso++) {
//     fHNTriggers[iIso] = NULL;
//   }
// }


// //________________________________________________________________________
// void AliConversionTrackCuts::CreateBaseHistograms() {
//   //Create histograms add, to outputlis

//   cout << "Creating histograms for "<< GetName() << endl;

//   fHistograms = new TList();
//   fHistograms->SetOwner(kTRUE);
//   fHistograms->SetName(fName);

//   for(int iIso = 0; iIso < 2; iIso++) {

//     fHNTriggers[iIso] = new TH1F(Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
// 								 Form("%s_%s_fNTriggers", fName.Data(), (iIso==0)?"nonIso":"isolated"), 
// 								 fAxistPt.GetNbins(), fAxistPt.GetXbins()->GetArray());
//     fHNTriggers[iIso]->Sumw2();
//     fHistograms->Add(fHNTriggers[iIso]);
  
//   }

//   fCorrSparse = CreateSparse(GetName(), GetTitle(), &fAxesList);
//   fHistograms->Add(fCorrSparse);

//   fTrackSparse = CreateSparse(Form("%s_%s", GetName(), "METrack"), Form("%s %s", GetTitle(), "ME Tracks"), &fTrackAxisList);
//   fHistograms->Add(fTrackSparse);

//   fTrigSparse = CreateSparse(Form("%s_%s", GetName(), "METrig"), Form("%s %s", GetTitle(), "ME Triggers"), &fTrigAxisList);
//   fHistograms->Add(fTrigSparse);

// }

// ///________________________________________________________________________
// THnSparseF * AliConversionTrackCuts::CreateSparse(TString nameString, TString titleString, TList * axesList) {
//   //Create sparse
//   const Int_t dim = axesList->GetSize();
  
//   cout << nameString << " " << titleString << " " <<   "    dimesion: " << dim << endl;

//   TAxis * axes[dim];
//   Int_t   bins[dim];
//   Double_t min[dim];
//   Double_t max[dim];

//   for(Int_t i = 0; i<dim; i++) {
// 	TAxis * axis = dynamic_cast<TAxis*>(axesList->At(i));
// 	if(axis) axes[i] = axis;
// 	else {
// 	  cout << "AliAnalysisTaskdPhi::CreateSparse: Error error, all the axes are not present in axis list" << endl;
// 	  return NULL;
// 	}
//   }

//   for(Int_t i = 0; i<dim; i++) {
// 	cout << axes[i]->GetTitle() << endl;
// 	bins[i] = axes[i]->GetNbins(); 
// 	min[i] = axes[i]->GetBinLowEdge(1);
// 	max[i] = axes[i]->GetBinUpEdge(axes[i]->GetNbins());
//   }

//   THnSparseF * sparse = new THnSparseF(Form("%s", nameString.Data()), 
// 									   Form("%s", titleString.Data()), 
// 									   dim, bins, min, max);
  
//   for(Int_t i = 0; i<dim; i++) {
// 	sparse->GetAxis(i)->SetNameTitle(axes[i]->GetName(), axes[i]->GetTitle() );
// 	if(axes[i]->GetXbins()->GetSize() > 0) {
// 	  sparse->SetBinEdges(i, axes[i]->GetXbins()->GetArray() );
// 	}
//   }
//   return sparse;
// }


// ///____________________________________________________________________________
// // void AliConversionTrackCuts::FillTriggerCounters(Float_t tPt, Bool_t isolated){ 
// //   //Fill histogram with trigger counters

// //   fHNTriggers[0]->Fill(tPt);
  
// //   if(isolated) {
// //     fHNTriggers[isolated]->Fill(tPt);
    
// //   }
// // }

// // ///_____________________________________________________________________________
// // void AliConversionTrackCuts::FillHistograms(Float_t tPt, Float_t cPt, Float_t dPhi, Float_t dEta, Bool_t isolated) {
// //   //Fill histograms

// //   if(dEta) { ;}
// //   //fHdPhi[0]->Fill(tPt, cPt, dPhi);
// //   if(isolated) {
// //     //fHdPhi[isolated]->Fill(tPt, cPt, dPhi);
// //   }
// // }

// //_______________________________________________________________________________

// void AliConversionTrackCuts::PrintStatistics()  { 
//   //Print some statistics between each file
//   for(Int_t i = 1; i <= fHNTriggers[0]->GetNbinsX(); i++) {
//     Int_t nTrig = (Int_t) fHNTriggers[0]->GetBinContent(i+1);
//     cout << "triggers: " << nTrig << endl;

//   }
// }


// //_______________________________________________________________________________
// void AliConversionTrackCuts::FillTriggerCounters(const AliAODConversionParticle * particle, Bool_t leading) {
//   fHNTriggers[leading]->Fill(particle->Pt());
// }


// //________________________________________________________________
// void AliConversionTrackCuts::CorrelateWithTracks(AliAODConversionParticle * particle, TObjArray * tracks, Int_t const tIDs[4], Int_t isolated = 0) {
//   //Correlate particle with tracks

//   FillTriggerCounters(particle, isolated);

//   Int_t nDim = fAxesList.GetSize();
//   Double_t dphivalues[nDim];
//   Double_t trackValues[nDim];
//   Double_t trigValues[nDim - 1];

//   for(int ij = 0; ij < tracks->GetEntriesFast(); ij++) {
// 	AliVTrack * track = static_cast<AliVTrack*>(tracks->UncheckedAt(ij));
// 	Int_t tid = track->GetID();

// 	if((tid > 0) && (tid == tIDs[0] || tid == tIDs[1] || tid == tIDs[2] || tid == tIDs[3]) ) {
// 	  continue;
// 	}
	
// 	dphivalues[0] = particle->Eta() - track->Eta();
// 	dphivalues[1] = GetDPhi(particle->Phi() - track->Phi());
// 	dphivalues[2] = particle->Pt();
// 	dphivalues[3] = track->Pt();
// 	dphivalues[4] = isolated;

// 	trackValues[0] = track->Eta();
// 	trackValues[1] = track->Phi();
// 	trackValues[2] = particle->Pt();
// 	trackValues[3] = track->Pt();
// 	trackValues[4] = isolated;

// 	trigValues[0] = particle->Eta();
// 	trigValues[1] = particle->Phi();
// 	trigValues[2] = particle->Pt();
// 	trigValues[4] = isolated;

// 	if(nDim > 4) {
// 	  dphivalues[5] = particle->M();
// 	  trackValues[5] = particle->M();
// 	  trigValues[4] = particle->M();
// 	}
	
// 	fCorrSparse->Fill(dphivalues);
// 	fTrackSparse->Fill(trackValues);
// 	fTrigSparse->Fill(trigValues);
//   }
// }

