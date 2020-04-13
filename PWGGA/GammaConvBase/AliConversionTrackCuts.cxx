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
/// @brief  Base class for analysation of conversion particle - track correlations - track cuts


#include "AliConversionTrackCuts.h"
//#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include <TFormula.h>
#include <iostream>
#include "TH2F.h"
#include "AliESDtrackCuts.h"
#include "THn.h"

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
  fDCAZmax(3.2*3.2),
  fDCAXYmax(2.4*2.4),
  fInitialized(kFALSE),
  fhPhi(NULL),
  //  fhPt(NULL),
  //fhPhiPt(NULL),
  fhdcaxyPt(NULL),
  fhdcazPt(NULL),
  fhdca(NULL),
  fhnclpt(NULL),
  fhnclsfpt(NULL),
  fhEtaPhi(NULL),
  fhTrackEff(NULL),
  fkCreateTrackEff(kFALSE),
  fHistograms(NULL) 
{
  //Constructor
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
  fInitialized(kFALSE),
  fhPhi(NULL),  
  //fhPt(NULL),
  //fhPhiPt(NULL),
  fhdcaxyPt(NULL),
  fhdcazPt(NULL),
  fhdca(NULL),
  fhnclpt(NULL),
  fhnclsfpt(NULL),
  fhEtaPhi(NULL),
  fhTrackEff(NULL),
  fkCreateTrackEff(kFALSE),
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

   if(fEsdTrackCuts)
     delete fEsdTrackCuts;
   fEsdTrackCuts = NULL;
   
   if(fEsdTrackCutsExtra1)
     delete fEsdTrackCutsExtra1;
   fEsdTrackCutsExtra1 = NULL;
   
   if(fEsdTrackCutsExtra2)
     delete fEsdTrackCutsExtra2;
   fEsdTrackCutsExtra2 = NULL;

}

//______________________________________________________________________________
void AliConversionTrackCuts::DefineESDCuts() {
  // Reproduces the cuts of the corresponding bit in the ESD->AOD filtering
  // (see $ALICE_ROOT/ANALYSIS/macros/AddTaskESDFilter.C)
  ///Copied from alianalyseleadingue
  const Int_t filterbit = fFilterBit;

  if (filterbit == 128) {
    if(!fEsdTrackCuts) {
      fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fEsdTrackCuts->SetMinNClustersTPC(70);
    }
  }  else if (filterbit == 256) {
    if(!fEsdTrackCuts) {
      // syst study
      fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fEsdTrackCuts->SetMinNClustersTPC(80);
      fEsdTrackCuts->SetMaxChi2PerClusterTPC(3);
      fEsdTrackCuts->SetMaxDCAToVertexZ(2.7);
      fEsdTrackCuts->SetMaxDCAToVertexXY(1.9);
    }
  }  else if (filterbit == 512) {
    if(!fEsdTrackCuts) {
      // syst study
      fEsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
      fEsdTrackCuts->SetMinNClustersTPC(60);
      fEsdTrackCuts->SetMaxChi2PerClusterTPC(5);
      fEsdTrackCuts->SetMaxDCAToVertexZ(3.7);
      fEsdTrackCuts->SetMaxDCAToVertexXY(2.9);
    }
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
Bool_t AliConversionTrackCuts::AcceptTrack(const AliESDtrack * track) {
  //Check esd track
  FillHistograms(kPreCut, track);

  if( fFilterBit == 256) {

    ///Standalone tpc tracks constrained
    const AliExternalTrackParam * param = track->GetConstrainedParam();
    if(param) {
      AliESDtrack esdTrack(track);
      esdTrack.CopyFromVTrack(param);

      if( !fEsdTrackCuts->IsSelected(&esdTrack)) return kFALSE;

      FillHistograms(1, track);

      Double_t dca[2];
      GetDCA(&esdTrack, dca);
      
      FillDCAHist(dca[1], dca[0], &esdTrack);
      if(fhEtaPhi) fhEtaPhi->Fill(esdTrack.Eta(), esdTrack.Phi());
      return kTRUE;
    } else {
      return kFALSE;
    }

    return kFALSE;
  }


  if(!fInitialized) {
    DefineESDCuts();
    // if(fDCAXYmax > 0) {
    //   if(fEsdTrackCuts) fEsdTrackCuts->SetMaxDCAToVertexXY(fDCAXYmax);
    // }
    // if(fDCAZmax > 0) {
    //   if(fEsdTrackCuts) fEsdTrackCuts->SetMaxDCAToVertexZ(fDCAZmax);
    // }
  
    fInitialized = kTRUE;
  }


  Double_t dca[2];
  GetDCA(track, dca);

  
  ///If only one track cuts then it has passed the cuts
  if( !(fEsdTrackCutsExtra1 && fEsdTrackCutsExtra2)) {
    FillHistograms(1, track);
    FillDCAHist(dca[1], dca[0], track);
    if(fhEtaPhi) fhEtaPhi->Fill(track->Eta(), track->Phi());
    return kTRUE;
  }

  ///If passing extra
  AliESDtrack internaltrack(track); // local copy used in order to keep constantnes of the input object
  if (fEsdTrackCutsExtra1 && fEsdTrackCutsExtra1->IsSelected(&internaltrack)) {
    FillHistograms(1, &internaltrack);
    FillHistograms(2, &internaltrack);

    FillDCAHist(dca[1], dca[0], &internaltrack);
    if(fhEtaPhi) fhEtaPhi->Fill(internaltrack.Eta(), internaltrack.Phi());
    
    return kTRUE;
  } 

  ///If passing extra2
  if (fEsdTrackCutsExtra2 && fEsdTrackCutsExtra2->IsSelected(&internaltrack)) {
    const AliExternalTrackParam * param = internaltrack.GetConstrainedParam();
    if(param) {
      AliESDtrack esdTrack(internaltrack);
      esdTrack.CopyFromVTrack(param);

      FillHistograms(3, &esdTrack);
      FillHistograms(1, &esdTrack);

      FillDCAHist(dca[1], dca[0], &esdTrack);
      if(fhEtaPhi) fhEtaPhi->Fill(esdTrack.Eta(), esdTrack.Phi());

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
  // return kTRUE;
}

Bool_t AliConversionTrackCuts::AcceptTrack(const AliAODTrack * track) {
  //Check aod track
  
  FillHistograms(kPreCut, track);
  
  if (fFilterBit == 768) {
    if(!track->IsHybridGlobalConstrainedGlobal()) return kFALSE;
      
    if (!(track->GetStatus() & AliVTrack::kITSrefit)) {
      return kFALSE;
    }
      
    //The cluster sharing cut can be done with:
    Double_t frac = Double_t(track->GetTPCnclsS()) / Double_t(track->GetTPCncls());
    if (frac > 0.4) return kFALSE;
      
    ///Do dca xy cut!
    FillHistograms(1, track);
      
    ///DCA
    Double_t dca[2] = { -999, -999};
    //Bool_t dcaok = 
    GetDCA(track, dca);
    FillDCAHist(dca[1], dca[0], track);
      
      
    if(track->IsGlobalConstrained()) {
      FillHistograms(3, track);
    } else {
      FillHistograms(2, track);
    }
      
    if(fhEtaPhi) fhEtaPhi->Fill(track->Eta(), track->Phi());
      
    return kTRUE;
    
    ////////////////////////////////
    //// Standalone
    ////////////////////////////////
  } else  if(fFilterBit == 256) {
    if(!track->IsTPCConstrained()) return kFALSE;



    ///DCA
    Double_t dca[2] = { -999, -999};
    GetDCA(track, dca);

    if( (dca[0]*dca[0]/fDCAXYmax + dca[1]*dca[1]/fDCAZmax) > 1 ) {
      FillHistograms(3, track);
      return kFALSE;
    }

    if(track->GetTPCncls() < 70) {
      FillHistograms(4, track);
      return kFALSE;
    }

    AliAODVertex * vtx = track->GetProdVertex();
    if (vtx->GetType() == AliAODVertex::kKink ) {
      FillHistograms(5, track);
      return kFALSE;
    }

    if(track->Chi2perNDF() > 36) {
      FillHistograms(6, track);
      return kFALSE;
    }
    if(track->Chi2perNDF() > 26) {
      FillHistograms(7, track);
      return kFALSE;
    }
    if(track->Chi2perNDF() > 16) {
      FillHistograms(8, track);
      return kFALSE;
    }
    if(track->Chi2perNDF() > 4) {
      FillHistograms(9, track);
      return kFALSE;
    }



    FillDCAHist(dca[1], dca[0], track);

    FillHistograms(2, track);
    if(fhEtaPhi) fhEtaPhi->Fill(track->Eta(), track->Phi());
    return kTRUE;

  }
  return kFALSE;
}


///______________________________________________________________________________
Bool_t AliConversionTrackCuts::GetDCA(const AliESDtrack *track, Double_t dcaxyz[2]) {
  ///Get track dca esd trck
  Float_t dca[2];
  Float_t bCov[3];
  track->GetImpactParameters(dca,bCov);
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
    return kFALSE;
  }

  dcaxyz[0] = dca[0];
  dcaxyz[1] = dca[1];
  
  return kTRUE;
}

///_____________________________________________________________________________
Bool_t AliConversionTrackCuts::GetDCA(const AliAODTrack *track, Double_t dca[2]) {
  ///Get track dca aod trck
  if(track->TestBit(AliAODTrack::kIsDCA)){
    dca[0]=track->DCA();
    dca[1]=track->ZAtDCA();
    return kTRUE;
  }
  
  Bool_t ok=kFALSE;
  if(fEvent) {
    Double_t covdca[3];
    //AliAODTrack copy(*track);
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);
    
    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      dca[0]=-999.;
      dca[1]=-999.;
    //printf("This method can be used only for propagation inside the beam pipe \n");
    return kFALSE;
    }


    AliAODVertex *vtx =(AliAODVertex*)(fEvent->GetPrimaryVertex());
    Double_t fBzkG = fEvent->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,dca,covdca);
    //ok = copy.PropagateToDCA(vtx,fBzkG,kVeryBig,dca,covdca);
  }
  if(!ok){
    dca[0]=-999.;
    dca[1]=-999.;
  }
  return ok;
}



TList * AliConversionTrackCuts::CreateHistograms() {
  //Create the histograms

  if(!fHistograms) fHistograms = new TList();

  fHistograms->SetOwner(kTRUE);
  fHistograms->SetName("trackCuts");

  fhPhi = new TH2F(Form("phi_%s", GetName()), Form("phi_%s", GetTitle()), 5, -0.5, 4.5, 32, 0, TMath::TwoPi());
  // TAxis * xax = fhPhi->GetXaxis();
  // for(Int_t i = 0; i < kNCuts; i++){
  // 	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  // }
  fHistograms->Add(fhPhi);
  

  fhEtaPhi = new TH2F(Form("etaphi_%s",GetName()), Form("etaphi_%s", GetTitle()), 36, -0.9, 0.9, 32, 0, TMath::TwoPi());
  fHistograms->Add(fhEtaPhi);

  // fhPt = new TH2F("pt", "pt", kNCuts+2, kPreCut -0.5, kNCuts + 0.5, 
  // 				  20, 0., 20.);
  // xax = fhPt->GetXaxis();
  // for(Int_t i = 0; i < kNCuts; i++){
  // 	xax->SetBinLabel(xax->FindFixBin(i), fgkCutNames[i]);
  // }
  // fHistograms->Add(fhPt);

  //  fhPhiPt = new TH2F("phipt", "phipt", 100, 0, 100, 64, 0, TMath::TwoPi());
  //fHistograms->Add(fhPhiPt);

  fhdcaxyPt = new TH2F(Form("dcaxypt_%s", GetName()),  Form("dcaxypt_%s", GetTitle()), 20, 0, 20, 50, -2.5, 2.5);
  fHistograms->Add(fhdcaxyPt);

  fhdcazPt = new TH2F(Form("dcazpt_%s", GetName()),  Form("dcazpt_%s", GetTitle()), 20, 0, 20, 70, -3.5, 3.5);
  fHistograms->Add(fhdcazPt);

  fhdca = new TH2F(Form("dca_%s", GetName()),  Form("dca_%s", GetTitle()), 70, -3.5, 3.5, 50, -2.5, 2.5);
  fhdca->SetXTitle("dca z");
  fhdca->SetYTitle("dca xy");

  
  fHistograms->Add(fhdca);

  // fhnclpt = new TH2F("nclstpcvspt", "nclstpcvspt", 20, 0, 20, 50, 0, 100);
  // fHistograms->Add(fhnclpt);

  // fhnclsfpt = new TH2F("nclsfpt", "nclsfpt", 20, 0, 20, 60, 0, 1.2);
  // fHistograms->Add(fhnclsfpt);

  

  if (fkCreateTrackEff) {
    const Double_t ptbins[23] = {0.5,  0.6, 0.7,  0.8, 0.9,  1.0, 1.2, 1.4,  1.6, 1.8, 
			  2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 
			  4.5, 4.75, 5.0};
    
    const Int_t bins[4] = { 12,  22, 36, 32};
    const Double_t min[4] = { -12, 0.5, -.9,  0};
    const Double_t max[4] = {  12,   5,  .9,  2*TMath::Pi() };
    
    fhTrackEff = new THnF(Form("hTrackEff_%s", GetName()), "hTrackEff", 4, bins, min, max);
    fhTrackEff->SetBinEdges(1, ptbins);
    fHistograms->Add(fhTrackEff);
  }

  
  return fHistograms;
}

void AliConversionTrackCuts::FillHistograms(Int_t cutIndex, const AliVTrack * track) {

  //Fill histograms
  if(fhPhi) fhPhi->Fill(cutIndex, track->Phi());
  //  if(fhPt) fhPt->Fill(cutIndex, track->Pt());
  //if(passed) fhPhiPt->Fill(track->Pt(), track->Phi());

}

void AliConversionTrackCuts::FillDCAHist(Float_t dcaz, Float_t dcaxy, const AliVTrack * track) {
  if(fhdcaxyPt) fhdcaxyPt->Fill(track->Pt(), dcaxy);
  if(fhdcazPt) fhdcazPt->Fill(track->Pt(), dcaz);
  if(fhdca) fhdca->Fill(dcaz, dcaxy);
  
  if(fhTrackEff) {
    Double_t val[4];
    val[0] = fEvent->GetPrimaryVertex()->GetZ();
    val[1] =  track->Pt();
    val[2] =  track->Eta();
    val[3] =  track->Phi();
    
    fhTrackEff->Fill(val);
  }
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




