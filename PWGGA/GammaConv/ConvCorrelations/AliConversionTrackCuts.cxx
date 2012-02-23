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
  fDCAXYmax(1E20),
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
  
  if (track->GetTPCNcls() < fTPCminNClusters) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Too few TPC clusters. Rejected");
	return kFALSE;
  }
  cutIndex++;
  
  if (track->Chi2perNDF() > fTPCmaxChi2) {
  FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Bad chi2. Rejected");
	return kFALSE;
  }
  cutIndex++;

  AliAODVertex *vertex = track->GetProdVertex();
  if (vertex && fRejectKinkDaughters) {
	if (vertex->GetType() == AliAODVertex::kKink) {
	  FillHistograms(cutIndex, track);
	  AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
	  return kFALSE;
	}
  }
  cutIndex++;

  if(track->ZAtDCA() > fDCAZmax) {
	FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
	return kFALSE;
  }
  cutIndex++;

  Float_t xatdca = track->XAtDCA();
  Float_t yatdca = track->YAtDCA();
  
  if(xatdca*xatdca * yatdca*yatdca > fDCAXYmax) {
	FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
	return kFALSE;
  }
  cutIndex++;



  ULong_t status = track->GetStatus();
  if ((status&AliESDtrack::kTPCrefit) == 0) {
	FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
	return kFALSE;
  }
  cutIndex++;

  if ((status&AliESDtrack::kITSrefit) == 0) {
	FillHistograms(cutIndex, track);
	AliDebug(AliLog::kDebug + 2, "Kink daughter. Rejected");
	return kFALSE;
  }
  cutIndex++;



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






