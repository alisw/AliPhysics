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
#include "TH2F.h"


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
  fAODTestFilterBit(-1),
  fhPhi(NULL),
  fHistograms(NULL)
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
  fAODTestFilterBit(-1),
  fhPhi(NULL),
  fHistograms(NULL)
{
  //Constructor
}


//________________________________________________________________________________
 AliConversionTrackCuts::~AliConversionTrackCuts() {
   ///destructor
   // if(fHistograms)
   // 	 delete fHistograms;
   // fHistograms = NULL;

}

TList * AliConversionTrackCuts::CreateHistograms() {

  if(!fHistograms) fHistograms = new TList();

  fHistograms->SetOwner(kTRUE);
  fHistograms->SetName("trackCuts");

  fhPhi = new TH2F("phi", "phi", 20, -0.5, 19.5, 128, 0, TMath::TwoPi());
  fHistograms->Add(fhPhi);

  fhPhiPt = new TH2F("phipt", "phipt", 80, 0, 100, 128, 0, TMath::TwoPi());
  fHistograms->Add(fhPhiPt);


  return fHistograms;
}


void AliConversionTrackCuts::FillHistograms(Int_t cutIndex, AliVTrack * track) {
  
  fhPhi->Fill(cutIndex, track->Phi());
  fhPhiPt->Fill(track->Pt(), track->Phi());
}

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

  Int_t cutIndex = 0;
  
  FillHistograms(cutIndex, track);
  cutIndex++;
  // step #0: check SPD and ITS clusters
  Int_t nSPD = 0;
  nSPD  = TESTBIT(track->GetITSClusterMap(), 0);
  nSPD += TESTBIT(track->GetITSClusterMap(), 1);
  if (nSPD < fSPDminNClusters) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Not enough SPD clusters in this track. Rejected");
	return kFALSE;
  }
  cutIndex++;
  
  // step #1: check number of clusters in TPC
  if (track->GetTPCNcls() < fTPCminNClusters) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Too few TPC clusters. Rejected");
	return kFALSE;
  }
  cutIndex++;
  
  if (track->GetITSNcls() < fITSminNClusters) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Too few ITS clusters. Rejected");
	return kFALSE;
  }
  cutIndex++;
  
  // step #2: check chi square
  if (track->Chi2perNDF() > fTPCmaxChi2) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
	return kFALSE;
  }
  cutIndex++;

  if (track->Chi2perNDF() > fITSmaxChi2) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
	return kFALSE;
  }
  cutIndex++;


   // step #3: reject kink daughters
   AliAODVertex *vertex = track->GetProdVertex();
   if (vertex && fRejectKinkDaughters) {
      if (vertex->GetType() == AliAODVertex::kKink) {
  FillHistograms(cutIndex, track);
         AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
         return kFALSE;
      }
   }
  cutIndex++;


   // step #4: DCA cut (transverse)
   Double_t b[2], cov[3];
   vertex = aodEvent->GetPrimaryVertex();
   if (!vertex) {
  FillHistograms(cutIndex, track);
      AliDebug(AliLog::kDebug + 2, "NULL vertex");
      return kFALSE;
   }
  cutIndex++;

   if (!track->PropagateToDCA(vertex, aodEvent->GetMagneticField(), kVeryBig, b, cov)) {
      AliDebug(AliLog::kDebug + 2, "Failed propagation to vertex");
  FillHistograms(cutIndex, track);
      return kFALSE;
   }
  cutIndex++;

   // if the DCA cut is not fixed, compute current value
   if (!fDCARfixed) {
	 FillHistograms(cutIndex, track);
      static TString str(fDCARptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaXY(Form("%s_dcaXY", GetName()), str.Data());
      fDCARmax = dcaXY.Eval(track->Pt());
   }
  cutIndex++;

   // check the cut
   if (TMath::Abs(b[0]) > fDCARmax) {
       FillHistograms(cutIndex, track);
	   AliDebug(AliLog::kDebug + 2, "Too large transverse DCA");
      return kFALSE;
   }
  cutIndex++;

 
   // step #5: DCA cut (longitudinal)
   // the DCA has already been computed above
   // if the DCA cut is not fixed, compute current value
   if (!fDCAZfixed) {
        FillHistograms(cutIndex, track);
		static TString str(fDCAZptFormula);
      str.ReplaceAll("pt", "x");
      static const TFormula dcaZ(Form("%s_dcaXY", GetName()), str.Data());
      fDCAZmax = dcaZ.Eval(track->Pt());
   }
  cutIndex++;

   // check the cut
  if (TMath::Abs(b[1]) > fDCAZmax) {
       FillHistograms(cutIndex, track);
	   AliDebug(AliLog::kDebug + 2, "Too large longitudinal DCA");
      return kFALSE;
   }

  cutIndex++;
 
   // step #6: check eta/pt range
   if (track->Eta() < fEta[0] || track->Eta() > fEta[1]) {
      FillHistograms(cutIndex, track);
	  AliDebug(AliLog::kDebug + 2, "Outside ETA acceptance");
      return kFALSE;
   }
  cutIndex++;

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






