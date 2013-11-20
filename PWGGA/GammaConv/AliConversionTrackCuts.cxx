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
  fEsdTrackCuts(NULL),
  fEsdTrackCutsExtra1(NULL),
  fEsdTrackCutsExtra2(NULL),
  fEvent(NULL),
  fFilterBit(2048),
  fDCAZmax(-1),
  fDCAXYmax(-1),
  fOwnedTracks(),
  fInitialized(kFALSE),
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
  fOwnedTracks.SetOwner(kTRUE);
}
//________________________________________________________________________
AliConversionTrackCuts::AliConversionTrackCuts(TString name, TString title = "title") : 
  AliAnalysisCuts(name, title),
  fEsdTrackCuts(NULL),
  fEsdTrackCutsExtra1(NULL),
  fEsdTrackCutsExtra2(NULL),
  fEvent(NULL),
  fFilterBit(2048),
  fDCAZmax(-1),
  fDCAXYmax(-1),
  fOwnedTracks(),
  fInitialized(kFALSE),
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
  fOwnedTracks.SetOwner(kTRUE);
}


//________________________________________________________________________________
 AliConversionTrackCuts::~AliConversionTrackCuts() {
   ///destructor
   // if(fHistograms)
   // 	 delete fHistograms;
   // fHistograms = NULL;

   if(fEsdTrackCuts)
     delete fEsdTrackCuts;
   fEsdTrackCuts = NULL;
   
   if(fEsdTrackCutsExtra1)
     delete fEsdTrackCutsExtra1;
   fEsdTrackCutsExtra1 = NULL;
   
   if(fEsdTrackCutsExtra2)
     delete fEsdTrackCutsExtra2;
   fEsdTrackCutsExtra2 = NULL;

   fOwnedTracks.Delete();
}

//______________________________________________________________________________
void AliConversionTrackCuts::DefineESDCuts() {
  // Reproduces the cuts of the corresponding bit in the ESD->AOD filtering
  // (see $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C)
  ///Copied from alianalyseleadingue
  const Int_t filterbit = fFilterBit;

  if (filterbit == 128) {
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(70);
    
  }  else if (filterbit == 256) {
    // syst study
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(80);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(3);
    fEsdTrackCuts->SetMaxDCAToVertexZ(2.7);
    fEsdTrackCuts->SetMaxDCAToVertexXY(1.9);
    
  }  else if (filterbit == 512) {
    // syst study
    fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    fEsdTrackCuts->SetMinNClustersTPC(60);
    fEsdTrackCuts->SetMaxChi2PerClusterTPC(5);
    fEsdTrackCuts->SetMaxDCAToVertexZ(3.7);
    fEsdTrackCuts->SetMaxDCAToVertexXY(2.9);
    
  } else if (filterbit == 1024) {
    if(!fEsdTrackCuts) {
      fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fEsdTrackCuts->SetMinNClustersTPC(-1);
      fEsdTrackCuts->SetMinNCrossedRowsTPC(70);
      fEsdTrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    }
  } else if (filterbit == 2048)  {
    // mimic hybrid tracks 
    // correspond to esdTrackCutsHTG, but WITHOUT spd constraint. this is checked with the next object
    if(!fEsdTrackCuts) {
      fEsdTrackCuts = new AliESDtrackCuts();
      TFormula *f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
      fEsdTrackCuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep, 100);
      fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
      fEsdTrackCuts->SetRequireTPCStandAlone(kTRUE);
      fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
      fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
      fEsdTrackCuts->SetMaxFractionSharedTPCClusters(0.4);
      
      fEsdTrackCuts->SetMaxDCAToVertexXY(2.4);
      fEsdTrackCuts->SetMaxDCAToVertexZ(3.2);
      fEsdTrackCuts->SetDCAToVertex2D(kTRUE);
	
      fEsdTrackCuts->SetMaxChi2PerClusterITS(36);
      fEsdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
	
      fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);
	
      fEsdTrackCuts->SetEtaRange(-0.9, 0.9);
      fEsdTrackCuts->SetPtRange(0.1, 1000000.0);
	
      fEsdTrackCuts->SetRequireITSRefit(kFALSE); //not here, n
    }
    // Add SPD requirement 
    fEsdTrackCutsExtra1 = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
    fEsdTrackCutsExtra1->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    fEsdTrackCutsExtra1->SetRequireITSRefit(kTRUE);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra1 corresponds to esdTrackCutsHTG
    
    fEsdTrackCutsExtra2 = new AliESDtrackCuts("No_SPD", "Reject tracks with cluster in SPD");
    fEsdTrackCutsExtra2->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    // A track passing fEsdTrackCuts and fEsdTrackCutsExtra2 corresponds to esdTrackCutsHTGC and needs to be constrained
    
    
  }
}

//______________________________________________________________________________
Bool_t AliConversionTrackCuts::AcceptTrack(AliESDtrack * track) {
  //Check esd track

  if(!fInitialized) {
    DefineESDCuts();
    if(fDCAXYmax > 0) {
      if(fEsdTrackCuts) fEsdTrackCuts->SetMaxDCAToVertexXY(fDCAXYmax);
    }
    if(fDCAZmax > 0) {
      if(fEsdTrackCuts) fEsdTrackCuts->SetMaxDCAToVertexZ(fDCAZmax);
    }
  
    fInitialized = kTRUE;
  }

  FillHistograms(kPreCut, track);


  
  if( !fEsdTrackCuts->IsSelected(track)) return kFALSE;


  ///If only one track cuts then it has passed the cuts
  if( !(fEsdTrackCutsExtra1 && fEsdTrackCutsExtra2)) {
    FillHistograms(1, track);
    return kTRUE;
  }
  ///If passing extra
  if (fEsdTrackCutsExtra1 && fEsdTrackCutsExtra1->IsSelected(track)) {
    FillHistograms(2, track);
    FillHistograms(1, track);
    return kTRUE;
  } 

  ///If passing extra2
  if (fEsdTrackCutsExtra2 && fEsdTrackCutsExtra2->IsSelected(track)) {
    const AliExternalTrackParam * param = track->GetConstrainedParam();
    if(param) {
      AliESDtrack* esdTrack = new AliESDtrack(track);
      esdTrack->CopyFromVTrack(param);
      track = esdTrack;
      fOwnedTracks.Add(track);

      FillHistograms(3, track);
      FillHistograms(1, track);
      return kTRUE;
    } else {
      return kFALSE;
    }
  } else {
    return kFALSE;
  }

  cout << "error error, should not be herer!"<<endl;
  return kFALSE;

  // FillHistograms(kPreCut + 1, track);
  // return kTRUE;

  // fhnclpt->Fill(track->Pt(), track->GetTPCNcls());
  // if(track->GetTPCNclsF() > 0) fhnclsfpt->Fill(track->Pt(), ((Double_t) track->GetTPCNcls())/track->GetTPCNclsF());
  // FillHistograms(kPreCut + 1, track);

  // ///Get impact parameters
  // Double_t extCov[15];
  // track->GetExternalCovariance(extCov);
  // Float_t b[2];
  // Float_t bCov[3];
  // track->GetImpactParameters(b,bCov);
  // if (bCov[0]<=0 || bCov[2]<=0) {
  //   AliDebug(1, "Estimated b resolution lower or equal zero!");
  //   bCov[0]=0; bCov[2]=0;
  // }
  
  // Float_t dcaToVertexXY = b[0];
  // Float_t dcaToVertexZ = b[1];
  // FillDCAHist(dcaToVertexZ, dcaToVertexXY, track);
  // return kTRUE;
}

Bool_t AliConversionTrackCuts::AcceptTrack(AliAODTrack * track) {
  //Check aod track
  FillHistograms(kPreCut, track);
  if(!track->TestFilterBit(fFilterBit)) {
    return kFALSE;
  }

  ///Do dca xy cut!

  return kTRUE;


  // if (track->GetTPCNcls() < fTPCminNClusters) return kFALSE;
  // FillHistograms(kCutNcls, track);

  // if (track->Chi2perNDF() > fTPCmaxChi2) return kFALSE;
  // FillHistograms(kCutNDF, track);

  // AliAODVertex *vertex = track->GetProdVertex();
  // if (vertex && fRejectKinkDaughters) {
  //   if (vertex->GetType() == AliAODVertex::kKink) {
  //     return kFALSE;
  //   }
  // }
  // FillHistograms(kCutKinc, track);
  
  // if(TMath::Abs(track->ZAtDCA()) > fDCAZmax) {
  //   return kFALSE;
  // }
  // FillHistograms(kCutDCAZ, track);

  // Float_t xatdca = track->XAtDCA();
  // Float_t yatdca = track->YAtDCA();
  // Float_t xy = xatdca*xatdca + yatdca*yatdca;
  // if(xy > fDCAXYmax) {
  //   return kFALSE;
  // }

  // FillHistograms(kCutDCAXY, track);


  // fhnclpt->Fill(track->Pt(), track->GetTPCNcls());
  // if(track->GetTPCNclsF() > 0) fhnclsfpt->Fill(track->Pt(), ((Double_t) track->GetTPCNcls())/track->GetTPCNclsF());
  // FillDCAHist(track->ZAtDCA(), TMath::Sqrt(track->XAtDCA()*track->XAtDCA() + track->YAtDCA()*track->YAtDCA()), track);


}






TList * AliConversionTrackCuts::CreateHistograms() {
  //Create the histograms

  if(!fHistograms) fHistograms = new TList();

  fHistograms->SetOwner(kTRUE);
  fHistograms->SetName("trackCuts");

  fhPhi = new TH2F("phi", "phi", 5, -0.5, 4.5, 32, 0, TMath::TwoPi());
  // TAxis * xax = fhPhi->GetXaxis();
  // for(Int_t i = 0; i < kNCuts; i++){
  // 	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  // }
  fHistograms->Add(fhPhi);
  


  // fhPt = new TH2F("pt", "pt", kNCuts+2, kPreCut -0.5, kNCuts + 0.5, 
  // 				  20, 0., 20.);
  // xax = fhPt->GetXaxis();
  // for(Int_t i = 0; i < kNCuts; i++){
  // 	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  // }
  // fHistograms->Add(fhPt);

  //  fhPhiPt = new TH2F("phipt", "phipt", 100, 0, 100, 64, 0, TMath::TwoPi());
  //fHistograms->Add(fhPhiPt);

  // fhdcaxyPt = new TH2F("dcaxypt", "dcaxypt", 20, 0, 20, 50, 0, 5);
  // fHistograms->Add(fhdcaxyPt);

  // fhdcazPt = new TH2F("dcazpt", "dcazpt", 20, 0, 20, 50, 0, 5);
  // fHistograms->Add(fhdcazPt);

  // fhdca = new TH2F("dca", "dca", 60, -3, 3, 60, -3, 3);
  // fHistograms->Add(fhdca);

  // fhnclpt = new TH2F("nclstpcvspt", "nclstpcvspt", 20, 0, 20, 50, 0, 100);
  // fHistograms->Add(fhnclpt);

  // fhnclsfpt = new TH2F("nclsfpt", "nclsfpt", 20, 0, 20, 60, 0, 1.2);
  // fHistograms->Add(fhnclsfpt);
  
  return fHistograms;
}

void AliConversionTrackCuts::FillHistograms(Int_t cutIndex, AliVTrack * track) {

  //Fill histograms
  if(fhPhi) fhPhi->Fill(cutIndex, track->Phi());
  //  if(fhPt) fhPt->Fill(cutIndex, track->Pt());
  //if(passed) fhPhiPt->Fill(track->Pt(), track->Phi());

}

void AliConversionTrackCuts::FillDCAHist(Float_t dcaz, Float_t dcaxy, AliVTrack * track) {
  if(fhdcaxyPt) fhdcaxyPt->Fill(track->Pt(), dcaxy);
  if(fhdcazPt) fhdcazPt->Fill(track->Pt(), dcaz);
  if(fhdca) fhdca->Fill(dcaz, dcaxy);
}









//_________________________________________________________________________________________________
void AliConversionTrackCuts::Print(const Option_t *) const
{
//
// Print information on this cut
//

  // printf("Cut name                : %s \n", GetName());
  // printf("Kink daughters are      : %s \n", (fRejectKinkDaughters ? "rejected" : "accepted"));
  // printf("TPC requirements        : clusters/findable %f, min. cluster = %d, max chi2 = %f, %s require refit\n", fTPCClusOverFindable, fTPCminNClusters, fTPCmaxChi2, (fRequireTPCRefit) ? "" : "Don't");
  // printf("ITS requirements        : min. cluster = %d (all), %d (SPD), max chi2 = %f \n", fITSminNClusters, fSPDminNClusters, fITSmaxChi2);
  // printf("DCA z cut               : fixed to %f cm \n", fDCAZmax);
  // printf("DCA xy cut              : fixed to %f cm \n", fDCAXYmax)
    ;
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




