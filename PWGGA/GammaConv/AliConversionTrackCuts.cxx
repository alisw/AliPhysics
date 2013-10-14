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
#include "AliESDtrackCuts.h"

using namespace std;
ClassImp(AliConversionTrackCuts)


const char* AliConversionTrackCuts::fgkCutNames[AliConversionTrackCuts::kNCuts] = {
  "nClusTPC", 
  "FoundFindable", 
  "Chi2PerNDF", 
  "Kink", 
  "DCA_Z", 
  "DCA_XY", 
  "TPCRefit"
  "kAccTracks"
};



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
  fDCAXYmax(1E20),
  fSPDminNClusters(0),
  fITSminNClusters(0),
  fITSmaxChi2(1E20),
  fTPCminNClusters(0),
  fTPCClusOverFindable(0.0),
  fTPCmaxChi2(1E20),
  fAODTestFilterBit(-1),
  fRequireTPCRefit(kFALSE),
  fESDCuts(NULL),
  fhPhi(NULL),
  fhPt(NULL),
  fhPhiPt(NULL),
  fhdcaxyPt(NULL),
  fhdcazPt(NULL),
  fhdca(NULL),
  fhnclpt(NULL),
  fhnclsfpt(NULL),
  fHistograms(NULL) {
  //Constructor
  //SetUpAxes();
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
  fTPCClusOverFindable(0.0),
  fTPCmaxChi2(1E20),
  fAODTestFilterBit(-1),
  fRequireTPCRefit(kFALSE),
  fESDCuts(NULL),
  fhPhi(NULL),  
  fhPt(NULL),
  fhPhiPt(NULL),
  fhdcaxyPt(NULL),
  fhdcazPt(NULL),
  fhdca(NULL),
  fhnclpt(NULL),
  fhnclsfpt(NULL),
  fHistograms(NULL)
{
  //Constructor
//  SetUpAxes();
}


//________________________________________________________________________________
 AliConversionTrackCuts::~AliConversionTrackCuts() {
   ///destructor
   // if(fHistograms)
   // 	 delete fHistograms;
   // fHistograms = NULL;

   if(fESDCuts)
     delete fESDCuts;
   fESDCuts = NULL;
}

TList * AliConversionTrackCuts::CreateHistograms() {
  //Create the histograms

  if(!fHistograms) fHistograms = new TList();

  fHistograms->SetOwner(kTRUE);
  fHistograms->SetName("trackCuts");

  fhPhi = new TH2F("phi", "phi", kNCuts+2, kPreCut -0.5, kNCuts + 0.5, 
				   128, 0, TMath::TwoPi());
  TAxis * xax = fhPhi->GetXaxis();
  for(Int_t i = 0; i < kNCuts; i++){
	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  }
  fHistograms->Add(fhPhi);
  


  fhPt = new TH2F("pt", "pt", kNCuts+2, kPreCut -0.5, kNCuts + 0.5, 
				  100, 0., 100.);
  xax = fhPt->GetXaxis();
  for(Int_t i = 0; i < kNCuts; i++){
	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  }
  fHistograms->Add(fhPt);

  //  fhPhiPt = new TH2F("phipt", "phipt", 100, 0, 100, 64, 0, TMath::TwoPi());
  //fHistograms->Add(fhPhiPt);

  fhdcaxyPt = new TH2F("dcaxypt", "dcaxypt", 20, 0, 20, 50, 0, 5);
  fHistograms->Add(fhdcaxyPt);

  fhdcazPt = new TH2F("dcazpt", "dcazpt", 20, 0, 20, 50, 0, 5);
  fHistograms->Add(fhdcazPt);

  fhdca = new TH2F("dca", "dca", 60, -3, 3, 60, -3, 3);
  fHistograms->Add(fhdca);

  fhnclpt = new TH2F("nclstpcvspt", "nclstpcvspt", 20, 0, 20, 100, 0, 200);
  fHistograms->Add(fhnclpt);

  fhnclsfpt = new TH2F("nclsfpt", "nclsfpt", 20, 0, 20, 100, 0, 1.2);
  fHistograms->Add(fhnclsfpt);
  
  return fHistograms;
}


void AliConversionTrackCuts::FillHistograms(Int_t cutIndex, AliVTrack * track, Bool_t passed = kFALSE) {
  //Fill histograms
  (void)passed;
  fhPhi->Fill(cutIndex, track->Phi());
  fhPt->Fill(cutIndex, track->Pt());
  //if(passed) fhPhiPt->Fill(track->Pt(), track->Phi());

}

void AliConversionTrackCuts::FillDCAHist(Float_t dcaz, Float_t dcaxy, AliVTrack * track) {
  fhdcaxyPt->Fill(track->Pt(), dcaxy);
  fhdcazPt->Fill(track->Pt(), dcaz);
  fhdca->Fill(dcaz, dcaxy);
}


Bool_t AliConversionTrackCuts::AcceptTrack(AliESDtrack * track) {
  //Check esd track
  if(!fESDCuts)
    fESDCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  
  FillHistograms(kPreCut, track);

  if( !fESDCuts->IsSelected(track)) {
    return kFALSE;
    
  }

  fhnclpt->Fill(track->Pt(), track->GetTPCNcls());
  if(track->GetTPCNclsF() > 0) fhnclsfpt->Fill(track->Pt(), ((Double_t) track->GetTPCNcls())/track->GetTPCNclsF());
  FillHistograms(kPreCut + 1, track);

  ///Get impact parameters
  Double_t extCov[15];
  track->GetExternalCovariance(extCov);
  Float_t b[2];
  Float_t bCov[3];
  track->GetImpactParameters(b,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  
  Float_t dcaToVertexXY = b[0];
  Float_t dcaToVertexZ = b[1];
  FillDCAHist(dcaToVertexZ, dcaToVertexXY, track);
  return kTRUE;
}

Bool_t AliConversionTrackCuts::AcceptTrack(AliAODTrack * track) {
  //Check aod track
  FillHistograms(kPreCut, track);
  if(!track->TestFilterBit(128)) {
    return kFALSE;
  }

  if (track->GetTPCNcls() < fTPCminNClusters) return kFALSE;
  FillHistograms(kCutNcls, track);

  if (track->Chi2perNDF() > fTPCmaxChi2) return kFALSE;
  FillHistograms(kCutNDF, track);

  AliAODVertex *vertex = track->GetProdVertex();
  if (vertex && fRejectKinkDaughters) {
    if (vertex->GetType() == AliAODVertex::kKink) {
      return kFALSE;
    }
  }
  FillHistograms(kCutKinc, track);
  
  if(TMath::Abs(track->ZAtDCA()) > fDCAZmax) {
    return kFALSE;
  }
  FillHistograms(kCutDCAZ, track);

  Float_t xatdca = track->XAtDCA();
  Float_t yatdca = track->YAtDCA();
  Float_t xy = xatdca*xatdca + yatdca*yatdca;
  if(xy > fDCAXYmax) {
    return kFALSE;
  }

  FillHistograms(kCutDCAXY, track);


  fhnclpt->Fill(track->Pt(), track->GetTPCNcls());
  if(track->GetTPCNclsF() > 0) fhnclsfpt->Fill(track->Pt(), ((Double_t) track->GetTPCNcls())/track->GetTPCNclsF());
  FillDCAHist(track->ZAtDCA(), TMath::Sqrt(track->XAtDCA()*track->XAtDCA() + track->YAtDCA()*track->YAtDCA()), track);


  return kTRUE;
}

//_________________________________________________________________________________________________
void AliConversionTrackCuts::Print(const Option_t *) const
{
//
// Print information on this cut
//

  printf("Cut name                : %s \n", GetName());
  printf("Kink daughters are      : %s \n", (fRejectKinkDaughters ? "rejected" : "accepted"));
  printf("TPC requirements        : clusters/findable %f, min. cluster = %d, max chi2 = %f, %s require refit\n", fTPCClusOverFindable, fTPCminNClusters, fTPCmaxChi2, (fRequireTPCRefit) ? "" : "Don't");
  printf("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f \n", fITSminNClusters, fSPDminNClusters, fITSmaxChi2);
  printf("DCA z cut               : fixed to %f cm \n", fDCAZmax);
  printf("DCA xy cut              : fixed to %f cm \n", fDCAXYmax);
}
 
//_________________________________________________________________________________________________

Bool_t AliConversionTrackCuts::IsSelected(TObject * object ) {
  AliAODTrack * aodtrack = dynamic_cast<AliAODTrack*>(object);
  if (aodtrack) {
    return AcceptTrack(aodtrack);
  } else {
    AliESDtrack * track = dynamic_cast<AliESDtrack*>(object);
    if (track)
      return AcceptTrack(track);
  }

return kFALSE;
} 




