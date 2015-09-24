/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, proviyaded that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purapose. It is         *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// Analysis task for the systematic study of the uncertainties related to   //
// the tracking and ITS-TPC matching efficiency for different particle      //
// species.                                                                 //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


#include "Riostream.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3F.h"
#include "THn.h"
#include "TList.h"
#include "TMath.h"
#include "TParticlePDG.h"
//
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliPID.h"
#include "AliESDtrackCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDUtils.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"
//
#include "AliAnalysisTrackingUncertainties.h"


ClassImp(AliAnalysisTrackingUncertainties)

// histogram constants
const double AliAnalysisTrackingUncertainties::kMaxChi2 = 200;

//________________________________________________________________________
AliAnalysisTrackingUncertainties::AliAnalysisTrackingUncertainties() 
  : AliAnalysisTaskSE("TaskTestPA"), 
    fESD(0),
    fESDpid(0),
    fUtils(0),
    fMCtrue(0),
    fListHist(0),
    fESDtrackCuts(0),
    fMatchTr(),
    fMatchChi(),
    fExcludeMomFromChi2ITSTPC(0) 
{
  // default Constructor
  /* fast compilation test

     gSystem->Load("libANALYSIS");
     gSystem->Load("libANALYSISalice");
     .L AliAnalysisTrackingUncertainties.cxx++
   */
}


//________________________________________________________________________
AliAnalysisTrackingUncertainties::AliAnalysisTrackingUncertainties(const char *name) 
  : AliAnalysisTaskSE(name),
    fESD(0),
    fESDpid(0),
    fUtils(0),
    fMCtrue(0),
    fListHist(0),
    fESDtrackCuts(0),
    fMatchTr(),
    fMatchChi(),
    fExcludeMomFromChi2ITSTPC(0) 
{
  //
  // standard constructur which should be used
  //
  fMCtrue = kTRUE; 
  //
  // create track cuts
  //
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
  fESDtrackCuts->SetEtaRange(-1., 1.);
  //
  // analysis utils if needed
  //
  fUtils = new AliAnalysisUtils();
  //

  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());

}


//________________________________________________________________________
void AliAnalysisTrackingUncertainties::UserCreateOutputObjects() 
{
  //
  // Create histograms
  // Called once
  //
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);
  //
  // (1.) basic QA and statistics histograms
  //
  TH2F * histVertexSelection = new TH2F("histVertexSelection", "vertex selection; vertex z (cm); accepted/rejected", 100, -50., 50., 2, -0.5, 1.5);
  fListHist->Add(histVertexSelection);
  //
  // (2.) track cut variation histograms
  //
  InitializeTrackCutHistograms();
  //
  // (3.) ITS -> TPC matching histograms
  //
  //  Int_t    binsMatch[kNumberOfAxes] = { 10,   50,    20,            18,  6};
  //  Double_t minMatch[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  //  Double_t maxMatch[kNumberOfAxes]  = {200,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  //  TString axisNameMatch[kNumberOfAxes]  = {"matchChi2","pT","eta","phi","pid"};
  //  TString axisTitleMatch[kNumberOfAxes] = {"matchChi2","pT","eta","phi","pid"};
  //
  //  THnF * hBestMatch = new THnF("hBestMatch","ITS -> TPC matching ",kNumberOfAxes, binsMatch, minMatch, maxMatch);
  //  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
  //    hBestMatch->GetAxis(iaxis)->SetName(axisNameMatch[iaxis]);
  // hBestMatch->GetAxis(iaxis)->SetTitle(axisTitleMatch[iaxis]);
  //  }
  //  BinLogAxis(hBestMatch, 1);
  //  fListHist->Add(hBestMatch);
  //
  //
  //
  const int nbPt=40;
  const double ptMax=5;
  //
  TH2F * hNMatch    = new TH2F("hNMatch","N Matches",nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5);
  TH2F * hBestMatch = new TH2F("hBestMatch","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2); // OB
  TH2F * hBestMatch_cuts = new TH2F("hBestMatch_cuts","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatch  = new TH2F("hAllMatch","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatchGlo  = new TH2F("hAllMatchGlo","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hPtCorr_ITSTPC = new TH2F("hPtCorr_ITSTPC","PtCorr",nbPt,0,ptMax,nbPt,0,ptMax);
  TH2F * hdPtRel_ITSTPC = new TH2F("hdPtRel_ITSTPC","dPt/pt",nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax);
  TH2F * hdInvPtRel_ITSTPC = new TH2F("hdInvPtRel_ITSTPC","pt*dPt^{-1}",nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax);
  //
  fListHist->Add(hNMatch);
  fListHist->Add(hBestMatch);
  fListHist->Add(hBestMatch_cuts);
  fListHist->Add(hAllMatch);
  fListHist->Add(hAllMatchGlo);
  fListHist->Add(hPtCorr_ITSTPC);
  fListHist->Add(hdPtRel_ITSTPC);
  fListHist->Add(hdInvPtRel_ITSTPC);
  //
  //
  TH2F * hNMatchBg    = new TH2F("hNMatchBg","N Matches",nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5);
  TH2F * hBestMatchBg = new TH2F("hBestMatchBg","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hBestMatchBg_cuts = new TH2F("hBestMatchBg_cuts","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2); // OB
  TH2F * hAllMatchBg  = new TH2F("hAllMatchBg","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatchGloBg  = new TH2F("hAllMatchGloBg","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hdPtRelBg_ITSTPC = new TH2F("hdPtRelBg_ITSTPC","dPt/pt",nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax);
  TH2F * hdInvPtRelBg_ITSTPC = new TH2F("hdInvPtRelBg_ITSTPC","pt*dPt^{-1}",nbPt,0,ptMax,2*nbPt+1,-0.4*ptMax,0.4*ptMax);

  //cout<<" here 0 : hdPtRelBg_ITSTPC "<<hdPtRelBg_ITSTPC<<" hdInvPtRelBg_ITSTPC "<<hdInvPtRelBg_ITSTPC<<endl;

  fListHist->Add(hNMatchBg);
  fListHist->Add(hBestMatchBg);
  fListHist->Add(hBestMatchBg_cuts);
  fListHist->Add(hAllMatchBg);
  fListHist->Add(hAllMatchGloBg);
   fListHist->Add(hdPtRelBg_ITSTPC);
  fListHist->Add(hdInvPtRelBg_ITSTPC);
  //add default track cuts in the output list
  fListHist->Add(fESDtrackCuts);
  //
  // post data
  //
  PostData(1, fListHist);


}



//________________________________________________________________________
void AliAnalysisTrackingUncertainties::UserExec(Option_t *) 
{
  //
  // main event loop
  //
  fESD = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!fESD) {
    PostData(1, fListHist);
    return;
  }
  //
  if (!fESDpid) fESDpid = ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();
  //
  // Check Monte Carlo information and other access first:
  //
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    fMCtrue = kFALSE;
  }
  //
  // extract generated particles information
  //
  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent) {
    if (fMCtrue) return;
  }
  if (fMCtrue) {
    //
    stack = mcEvent->Stack();
    if (!stack) return;
    //
    for(Int_t i = 0; i < stack->GetNtrack(); i++) {
      /* at the moment nothing is needed here
      TParticle * trackMC = stack->Particle(i);
      Double_t rap = trackMC->Eta();
      Double_t y = trackMC->Y();
      Double_t pT  = trackMC->Pt();
      */
      //
    }

  }
  //
  if (!fESDtrackCuts) {
    PostData(1, fListHist);
    return;
  }
  //
  // monitor vertex position and selection
  //
  TH2F * histVertexSelection = (TH2F *) fListHist->FindObject("histVertexSelection");
  //
  Float_t vertexZ = 0.;
  if (IsVertexAccepted(fESD, vertexZ)) {
    histVertexSelection->Fill(vertexZ, 0);
  } else {
    histVertexSelection->Fill(vertexZ, 1);
    return;
  }
  //
  // fill track cut variation histograms
  //
  ProcessTrackCutVariation();
  //
  // fill ITS->TPC matching histograms
  //
  ProcessItsTpcMatching();
  
  // Post output data  
  PostData(1, fListHist);
  
}      


//________________________________________________________________________
void  AliAnalysisTrackingUncertainties::ProcessItsTpcMatching(){
  //
  // check how many its-sa tracks get matched to TPC
  //
  int ntr = fESD->GetNumberOfTracks();
  //
  // initialize histograms
  //
  TH2F * hNMatch         = (TH2F*) fListHist->FindObject("hNMatch");
  TH2F * hBestMatch      = (TH2F*) fListHist->FindObject("hBestMatch");
  TH2F * hBestMatch_cuts = (TH2F*) fListHist->FindObject("hBestMatch_cuts");
  TH2F * hAllMatch       = (TH2F*) fListHist->FindObject("hAllMatch");
  TH2F * hAllMatchGlo    = (TH2F*) fListHist->FindObject("hAllMatchGlo");  
  TH2F * hPtCorr_ITSTPC  = (TH2F*) fListHist->FindObject("hPtCorr_ITSTPC");
  TH2F * hdPtRel_ITSTPC  = (TH2F*) fListHist->FindObject("hdPtRel_ITSTPC");
  TH2F * hdInvPtRel_ITSTPC = (TH2F*) fListHist->FindObject("hdInvPtRel_ITSTPC");

  //
  TH2F * hNMatchBg          = (TH2F*) fListHist->FindObject("hNMatchBg");
  TH2F * hBestMatchBg       = (TH2F*) fListHist->FindObject("hBestMatchBg");
  TH2F * hBestMatchBg_cuts  = (TH2F*) fListHist->FindObject("hBestMatchBg_cuts");
  TH2F * hAllMatchBg        = (TH2F*) fListHist->FindObject("hAllMatchBg");
  TH2F * hAllMatchGloBg     = (TH2F*) fListHist->FindObject("hAllMatchGloBg");    
  TH2F * hdPtRelBg_ITSTPC    = (TH2F*) fListHist->FindObject("hdPtRelBg_ITSTPC");
  TH2F * hdInvPtRelBg_ITSTPC = (TH2F*) fListHist->FindObject("hdInvPtRelBg_ITSTPC");

  //cout<<" here 1: hdPtRelBg_ITSTPC "<<hdPtRelBg_ITSTPC<<" hdInvPtRelBg_ITSTPC "<<hdInvPtRelBg_ITSTPC<<endl;

  //
  for (int it=0;it<ntr;it++) {
    AliESDtrack* trSA = fESD->GetTrack(it);
    if (!trSA->IsOn(AliESDtrack::kITSpureSA) || !trSA->IsOn(AliESDtrack::kITSrefit)) continue;
    double pt = trSA->Pt();

    // OB - fiducial eta and pt cuts
    Double_t etaSA = trSA->Eta();
    // std::cout<<" etaSA "<<etaSA<<std::endl; // eta range up to +/- 1.4

    Double_t ptSA  = trSA->Pt();

    if(TMath::Abs(etaSA)>0.8) continue;

    //
    Int_t nmatch = 0;
    for (int i=kMaxMatch;i--;) {fMatchChi[i]=0; fMatchTr[i]=0;}
    for (int it1=0;it1<ntr;it1++){

      //std::cout<<" here 0, it1 "<<it1<<" it "<<it<<std::endl;

      if (it1==it) continue;
      
      //std::cout<<" here 2, it1 "<<it1<<" it "<<it<<std::endl;

      AliESDtrack* trESD = fESD->GetTrack(it1);
      if (!trESD->IsOn(AliESDtrack::kTPCrefit)) continue;
      
      //std::cout<<" call match: it1 "<<it1<<" it "<<it<<" nmatch "<<nmatch<<std::endl;
      Match(trSA,trESD, nmatch, fExcludeMomFromChi2ITSTPC);
      //std::cout<<" left match: it1 "<<it1<<" it "<<it<<" nmatch "<<nmatch<<std::endl;
    }
    //
    
    // std::cout<<" if "<<it<<" filling nmatch "<<nmatch<<" best chi2 "<<fMatchChi[0]<<std::endl;
    
    hNMatch->Fill(pt,nmatch);
    if (nmatch>0){
      hBestMatch->Fill(pt,fMatchChi[0]);
      hPtCorr_ITSTPC->Fill(pt,fMatchTr[0]->Pt()); 
      hdPtRel_ITSTPC->Fill(pt,(pt-fMatchTr[0]->Pt())/pt); 
      hdInvPtRel_ITSTPC->Fill(pt,pt*( 1/pt - (1/fMatchTr[0]->Pt()) )); 
    }
    
    if (nmatch>0 && fESDtrackCuts){
      
      if(fESDtrackCuts->AcceptTrack(fMatchTr[0])){
	hBestMatch_cuts->Fill(pt,fMatchChi[0]);
      }
    }
    
    //
    for (int imt=nmatch;imt--;) {
      hAllMatch->Fill(pt,fMatchChi[imt]);
      if (fMatchTr[imt]->IsOn(AliESDtrack::kITSrefit)) hAllMatchGlo->Fill(pt,fMatchChi[imt]);
    }
    //
    nmatch = 0;
    for (int i=kMaxMatch;i--;) {fMatchChi[i]=0; fMatchTr[i]=0;}
    for (int it1=0;it1<ntr;it1++) {
      if (it1==it) continue;
      AliESDtrack* trESD = fESD->GetTrack(it1);
      if (!trESD->IsOn(AliESDtrack::kTPCrefit)) continue;
      Match(trSA,trESD, nmatch, fExcludeMomFromChi2ITSTPC, TMath::Pi());
    }
    //
    hNMatchBg->Fill(pt,nmatch);
    if (nmatch>0){
      hBestMatchBg->Fill(pt,fMatchChi[0]);
      hdPtRelBg_ITSTPC->Fill(pt,(pt-fMatchTr[0]->Pt())/pt); 
      hdInvPtRelBg_ITSTPC->Fill(pt,pt*( 1/pt - (1/fMatchTr[0]->Pt()) )); 
    }

    if (nmatch>0 && fESDtrackCuts){
      if(fESDtrackCuts->AcceptTrack(fMatchTr[0])){
	hBestMatchBg_cuts->Fill(pt,fMatchChi[0]);
      }
    }

    for (int imt=nmatch;imt--;) {
      hAllMatchBg->Fill(pt,fMatchChi[imt]);
      if (fMatchTr[imt]->IsOn(AliESDtrack::kITSrefit)) hAllMatchGloBg->Fill(pt,fMatchChi[imt]);
    }
    //
  }


}


void AliAnalysisTrackingUncertainties::Match(const AliESDtrack* tr0, const AliESDtrack* tr1, Int_t& nmatch, 
					     Bool_t excludeMom, Double_t rotate) {
  //
  // check if two tracks are matching, possible rotation for combinatoric backgr.
  // 
  Float_t bField = fESD->GetMagneticField();
  //
  const AliExternalTrackParam* trtpc0 = tr1->GetInnerParam();
  if (!trtpc0) return;
  AliExternalTrackParam trtpc(*trtpc0);
  //
  if (TMath::Abs(rotate)>1e-5) {
    const double *par = trtpc.GetParameter();
    const double *cov = trtpc.GetCovariance();
    double alp = trtpc.GetAlpha() + rotate;
    trtpc.Set(trtpc.GetX(),alp,par,cov);
  }
  //
  if (!trtpc.Rotate(tr0->GetAlpha())) return;
  if (!trtpc.PropagateTo(tr0->GetX(),bField)) return;
  double chi2 = tr0->GetPredictedChi2(&trtpc);

  //std::cout<<" in Match, nmatch "<<nmatch<<" par[4] before "<<trtpc.GetParameter()[4]<<" chi2 "<<chi2<<endl;

  // OB chi2 excluding pt 
  if(excludeMom){
    ((double*)trtpc.GetParameter())[4] = tr0->GetParameter()[4]; // set ITS mom equal TPC mom
    chi2 = tr0->GetPredictedChi2(&trtpc);

    //std::cout<<" in Match, nmatch "<<nmatch<<" par[4] after "<<trtpc.GetParameter()[4]<<" tr0 mom "<<tr0->GetParameter()[4]
    //	     <<" chi2 "<<chi2<<std::endl;
  }


  if (chi2>kMaxChi2) return;

  // std::cout<<" found good match, tr1 "<<tr1<<" chi2 "<<chi2<<std::endl;
  // std::cout<<" before: fMatchChi[0]  "<<fMatchChi[0]<<" [1] "<<fMatchChi[1]
  // 	   <<" [2]  "<<fMatchChi[2]<<" [3] "<<fMatchChi[3]
  // 	   <<" [4]  "<<fMatchChi[4]<<std::endl; 

  // std::cout<<" before: fMatchTr[0]  "<<fMatchTr[0]<<" [1] "<<fMatchTr[1]
  // 	   <<" [2]  "<<fMatchTr[2]<<" [3] "<<fMatchTr[3]
  // 	   <<" [4]  "<<fMatchTr[4]<<std::endl; 

  //
  int ins;
  for (ins=0;ins<nmatch;ins++) if (chi2<fMatchChi[ins]) break;
  if (ins>=kMaxMatch) return;
  
  for (int imv=nmatch;imv>ins;imv--) {
    if (imv>=kMaxMatch) continue;
    fMatchTr[imv]  = fMatchTr[imv-1];
    fMatchChi[imv] = fMatchChi[imv-1];
  }
  fMatchTr[ins] = tr1;
  fMatchChi[ins] = chi2;
  nmatch++;
  if (nmatch>=kMaxMatch) nmatch = kMaxMatch;
  //
}


//________________________________________________________________________
void AliAnalysisTrackingUncertainties::ProcessTrackCutVariation() {
  //
  // fill track cut variation histograms - undo cuts step-by-step and fill histograms
  //
  //
  // initialize histograms
  //
  THnF * histNcl         = (THnF *) fListHist->FindObject("histNcl");
  THnF * histChi2Tpc     = (THnF *) fListHist->FindObject("histChi2Tpc");
  THnF * histDcaZ        = (THnF *) fListHist->FindObject("histDcaZ");
  THnF * histSpd         = (THnF *) fListHist->FindObject("histSpd");
  THnF * histNcr         = (THnF *) fListHist->FindObject("histNcr");
  THnF * histCRoverFC    = (THnF *) fListHist->FindObject("histCRoverFC");
  THnF * histChi2Its     = (THnF *) fListHist->FindObject("histChi2Its");
  THnF * histTpcLength   = (THnF *) fListHist->FindObject("histTpcLength");
  THnF * histTpcItsMatch = (THnF *) fListHist->FindObject("histTpcItsMatch");
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  //
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    //
    AliESDtrack *track =fESD->GetTrack(i);
    //
    // relevant variables
    //
    //Double_t pid        = Double_t(GetPid(track));
    //
    Int_t nclsTPC       = track->GetTPCncls();
    Float_t pT          = track->Pt();
    Float_t eta         = track->Eta();
    Float_t phi         = track->Phi();
    Float_t chi2TPC     = track->GetTPCchi2();
    Float_t ncrTPC      = track->GetTPCCrossedRows();
    Int_t nclsTPCF      = track->GetTPCNclsF(); 
    Float_t nCRoverFC   = track->GetTPCCrossedRows();
    Double_t chi2ITS    = track->GetITSchi2();
    Int_t nclsITS       = track->GetITSclusters(0);
    Float_t tpcLength   = 0.;

    if (track->GetInnerParam() && track->GetESDEvent()) {
      tpcLength = track->GetLengthInActiveZone(1, 1.8, 220, track->GetESDEvent()->GetMagneticField());
    }

    if (nclsTPC != 0) {
      chi2TPC /= nclsTPC; 
    } else {
      chi2TPC = 999.;
    }

    if (nclsTPCF !=0) {
      nCRoverFC /= nclsTPCF;
    } else {
      nCRoverFC = 999.;
    }

    if (nclsITS != 0){
      chi2ITS /= nclsITS;
    }else {
      chi2ITS = 999.;
    }
    //
    track->GetImpactParameters(dca, cov);
    //
    // (1.) fill number of clusters histogram
    //
    Int_t minNclsTPC = fESDtrackCuts->GetMinNClusterTPC();
    fESDtrackCuts->SetMinNClustersTPC(0);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistNcl[kNumberOfAxes] = {static_cast<Double_t>(nclsTPC), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histNcl->Fill(vecHistNcl);
      }
    }
    fESDtrackCuts->SetMinNClustersTPC(minNclsTPC);
    //
    // (2.) fill chi2 TPC histogram
    //
    Float_t maxChi2 = fESDtrackCuts->GetMaxChi2PerClusterTPC();
    fESDtrackCuts->SetMaxChi2PerClusterTPC(999.);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistChi2Tpc[kNumberOfAxes] = {chi2TPC, pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histChi2Tpc->Fill(vecHistChi2Tpc);
      }
    }
    fESDtrackCuts->SetMaxChi2PerClusterTPC(maxChi2);
    //
    // (3.) fill dca_z histogram
    //
    Float_t maxDcaZ = fESDtrackCuts->GetMaxDCAToVertexZ();
    fESDtrackCuts->SetMaxDCAToVertexZ(999.);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistDcaZ[kNumberOfAxes] = {TMath::Abs(dca[1]), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histDcaZ->Fill(vecHistDcaZ);
      }
    }
    fESDtrackCuts->SetMaxDCAToVertexZ(maxDcaZ);
    //
    // (4.) fill hit in SPD histogram
    //
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    if (fESDtrackCuts->AcceptTrack(track)) {
      Int_t hasPoint = 0;
      if (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) hasPoint = 1;
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistSpd[kNumberOfAxes] = {static_cast<Double_t>(hasPoint), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histSpd->Fill(vecHistSpd);
      }
    }
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    //
    // (5.) fill number of crossed rows histogram
    //
    //Int_t minNcrTPC = fESDtrackCuts->GetMinNCrossedRowsTPC(); //wait for getter in ESDtrackCuts
    Int_t minNcrTPC = 0;  //for now use standard value from 2010 !!
    fESDtrackCuts->SetMinNCrossedRowsTPC(0);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistNcr[kNumberOfAxes] = {static_cast<Double_t>(ncrTPC), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histNcr->Fill(vecHistNcr);
      }
    }
    fESDtrackCuts->SetMinNCrossedRowsTPC(minNcrTPC);
    //
    // (6.) fill crossed rows over findable clusters histogram
    //
    //Int_t minCRoverFC = fESDtrackCuts->GetMinRatioCrossedRowsOverFindableClustersTPC(); //wait for getter in ESDtrackCuts
    Int_t minCRoverFC = 0.; //for now use standard value from 2010 !!
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.);
    if (fESDtrackCuts->AcceptTrack(track)) {
       for(Int_t iPid = 0; iPid < 6; iPid++) {
	 Double_t vecHistCRoverFC[kNumberOfAxes] = {static_cast<Double_t>(nCRoverFC), pT, eta, phi, static_cast<Double_t>(iPid)};
	 if (IsConsistentWithPid(iPid, track)) histCRoverFC->Fill(vecHistCRoverFC);
       }
    }
    fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(minCRoverFC);
    //
    // (7.) fill chi2 ITS histogram
    //
    Float_t maxChi2ITS = fESDtrackCuts->GetMaxChi2PerClusterITS();
    fESDtrackCuts->SetMaxChi2PerClusterITS(999.);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistChi2ITS[kNumberOfAxes] = {chi2ITS, pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histChi2Its->Fill(vecHistChi2ITS);
      }
    }
    fESDtrackCuts->SetMaxChi2PerClusterITS(maxChi2ITS);
    //
    // (8.) fill active length in TPC histogram
    //
    Int_t minTpcLength = fESDtrackCuts->GetMinLengthActiveVolumeTPC();
    fESDtrackCuts->SetMinLengthActiveVolumeTPC(0);
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistTpcLength[kNumberOfAxes] = {tpcLength, pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histTpcLength->Fill(vecHistTpcLength);
      }
    }
    fESDtrackCuts->SetMinLengthActiveVolumeTPC(minTpcLength);
    //
    // (9.) fill TPC->ITS matching efficiency histogram
    //
    Bool_t isMatched = kFALSE;
    // remove all ITS requirements
    //
    // Leonardo and Emilia: 
    //  -> if MC is available: fill it only for true primaries, 
    //        --to be done for every cut?
    //  -> Postprocessing: plot histogram with 1 divided by histogram with 0 as a function of pT/eta/phi
    //  -> Do we want to remove the DCA cut?
    Bool_t refit=fESDtrackCuts->GetRequireITSRefit();
    Float_t chi2tpc= fESDtrackCuts->GetMaxChi2TPCConstrainedGlobal();
    Float_t chi2its= fESDtrackCuts->GetMaxChi2PerClusterITS();
    //TString str = fESDtrackCuts->GetMaxDCAToVertexXYPtDep();
    
    fESDtrackCuts->SetRequireITSRefit(kFALSE);
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(99999.);
    fESDtrackCuts->SetMaxChi2PerClusterITS(999999.);
	//TString str = fESDtrackCuts->GetMaxDCAToVertexXYPtDep();
    //fESDtrackCuts->SetMaxDCAToVertexXYPtDep();
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistTpcItsMatch[kNumberOfAxes] = {static_cast<Double_t>(isMatched), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histTpcItsMatch->Fill(vecHistTpcItsMatch); // fill with 1 here
      }
    }
    //apply back the cuts
    fESDtrackCuts->SetRequireITSRefit(refit);
    fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(chi2tpc);
    fESDtrackCuts->SetMaxChi2PerClusterITS(chi2its);
    //fESDtrackCuts->SetMaxDCAToVertexXYPtDep(str.Data());
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    //set is matched
    isMatched=kTRUE;
    if (fESDtrackCuts->AcceptTrack(track)) {
      for(Int_t iPid = 0; iPid < 6; iPid++) {
	Double_t vecHistTpcItsMatch[kNumberOfAxes] = {static_cast<Double_t>(isMatched), pT, eta, phi, static_cast<Double_t>(iPid)};
	if (IsConsistentWithPid(iPid, track)) histTpcItsMatch->Fill(vecHistTpcItsMatch); // fill with 0 here
      }
    }
    
  } // end of track loop


}



//________________________________________________________________________
void AliAnalysisTrackingUncertainties::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query


}


//________________________________________________________________________
void AliAnalysisTrackingUncertainties::BinLogAxis(const THn *h, Int_t axisNumber) {
  //
  // Method for the correct logarithmic binning of histograms
  //
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
  
}


//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsVertexAccepted(AliESDEvent * esd, Float_t &vertexZ) {
  //
  // function to check if a proper vertex is reconstructed and write z-position in vertexZ
  //
  vertexZ = -999.;
  Bool_t vertexOkay = kFALSE;
  const AliESDVertex *vertex = esd->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    //
    vertex = esd->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) {
      vertexOkay = kFALSE; }
    else {
      vertexOkay = kTRUE;
    }
    //
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vertex->IsFromVertexerZ() && (zRes>0.25)) vertexOkay = kFALSE;
  }
  else {
    vertexOkay = kTRUE;
  }

  vertexZ = vertex->GetZ();  
  return vertexOkay;

}

//________________________________________________________________________
AliAnalysisTrackingUncertainties::ESpecies_t AliAnalysisTrackingUncertainties::GetPid(const AliESDtrack * const tr, Bool_t useTPCTOF) const {
    //
    // Determine particle species for a given track
    // Two approaches can be used: As default the selection is done using TPC-only, in addition
    // the TOF usage is optional. In case of TPC-TOF, a valid TOF signal has to be provided for 
    // the given track. The identification is delegated to helper function for each species. 
    // Tracks which are selected as more than one species (ambiguous decision) are rejected.
    //
    // @Return: Particles species (kUndef in case no identification is possible)
    //
    if(!fESDpid) return kUndef;
    if(useTPCTOF && !(tr->GetStatus() & AliVTrack::kTOFpid)) return kUndef;

    Bool_t isElectron(kFALSE), isPion(kFALSE), isKaon(kFALSE), isProton(kFALSE);
    Int_t nspec(0);
    if((isElectron = IsElectron(tr, useTPCTOF))) nspec++;
    if((isPion = IsPion(tr, useTPCTOF))) nspec++;
    if((isKaon = IsKaon(tr, useTPCTOF))) nspec++;
    if((isProton = IsProton(tr,useTPCTOF))) nspec++;
    if(nspec != 1) return kUndef;   // No decision or ambiguous decision;
    if(isElectron) return kSpecElectron;
    if(isPion) return kSpecPion;
    if(isProton) return kSpecProton;
    if(isKaon) return kSpecKaon;
    return kUndef;
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsElectron(const AliESDtrack * const tr, Bool_t useTPCTOF) const {
    //
    // Selection of electron candidates using the upper half of the TPC sigma band, starting at 
    // the mean ignoring its shift, and going up to 3 sigma above the mean. In case TOF information 
    // is available, tracks which are incompatible with electrons within 3 sigma are rejected. If 
    // no TOF information is used, the momentum regions where the kaon and the proton line cross 
    // the electron line are cut out using a 3 sigma cut around the kaon or proton line.
    //

    Float_t nsigmaElectronTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kElectron);
    if(nsigmaElectronTPC < 0 || nsigmaElectronTPC  > 3) return kFALSE;

    if(useTPCTOF){
        Float_t nsigmaElectronTOF = fESDpid->NumberOfSigmasTOF(tr, AliPID::kElectron);
        if(TMath::Abs(nsigmaElectronTOF) > 3) return kFALSE;
        else return kTRUE;
    } else {
        Float_t nsigmaKaonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kKaon),
                nsigmaProtonTPC =fESDpid->NumberOfSigmasTPC(tr, AliPID::kProton);
        if(TMath::Abs(nsigmaKaonTPC < 3) || TMath::Abs(nsigmaProtonTPC < 3)) return kFALSE;
        else return kTRUE;
    }
}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsConsistentWithPid(Int_t type, const AliESDtrack * const tr) {
  //
  // just check if the PID is consistent with a given hypothesis in order to 
  // investigate effects which are only dependent on the energy loss.
  //
  if (type == kSpecPion)     return IsPion(tr);
  if (type == kSpecKaon)     return IsKaon(tr);
  if (type == kSpecProton)   return IsProton(tr);
  if (type == kSpecElectron) return IsElectron(tr);
  if (type == kAll)          return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsPion(const AliESDtrack * const tr, Bool_t /*useTPCPTOF*/) const{
  //
  // Selectron of pion candidates
  // @TODO: To be implemented
  //
  Float_t nsigmaPionTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kPion);
  if (TMath::Abs(nsigmaPionTPC) < 3) return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsKaon(const AliESDtrack * const tr, Bool_t /*useTPCPTOF*/) const {
  //
  // Selection of kaon candidates
  // @TODO: To be implemented
  //
  Float_t nsigmaKaonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kKaon);
  if (TMath::Abs(nsigmaKaonTPC) < 3) return kTRUE;
  return kFALSE;

}

//________________________________________________________________________
Bool_t AliAnalysisTrackingUncertainties::IsProton(const AliESDtrack * const tr, Bool_t /*useTPCPTOF*/) const{
  // 
  // Selection of proton candidates
  // @TODO: To be implemented
  //
  Float_t nsigmaProtonTPC = fESDpid->NumberOfSigmasTPC(tr, AliPID::kProton);
  if (TMath::Abs(nsigmaProtonTPC) < 3) return kTRUE;
  return kFALSE;
  
}

//________________________________________________________________________
void AliAnalysisTrackingUncertainties::InitializeTrackCutHistograms() {
  //
  // create histograms for the track cut studies
  //
  //
  // (1.) number of clusters       
  //                               0-ncl, 1-pt, 2-eta,         3-phi, 4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all)
  Int_t    binsNcl[kNumberOfAxes] = { 40,   50,    20,            18,  6};
  Double_t minNcl[kNumberOfAxes]  = { 20,  0.1,    -1,             0, -0.5};
  Double_t maxNcl[kNumberOfAxes]  = {160,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameNcl[kNumberOfAxes]  = {"ncl","pT","eta","phi","pid"};
  TString axisTitleNcl[kNumberOfAxes] = {"ncl","pT","eta","phi","pid"};
  //
  THnF * histNcl = new THnF("histNcl","number of clusters histogram",kNumberOfAxes, binsNcl, minNcl, maxNcl);
  BinLogAxis(histNcl, 1);
  fListHist->Add(histNcl);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histNcl->GetAxis(iaxis)->SetName(axisNameNcl[iaxis]);
    histNcl->GetAxis(iaxis)->SetTitle(axisTitleNcl[iaxis]);
  }
  //
  // (2.) chi2/cls-TPC            
  //                                  0-chi2, 1-pt, 2-eta,         3-phi, 4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all)
  Int_t    binsChi2Tpc[kNumberOfAxes] = { 40,   50,    20,            18,  6};
  Double_t minChi2Tpc[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  Double_t maxChi2Tpc[kNumberOfAxes]  = {  8,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameChi2Tpc[kNumberOfAxes]  = {"chi2tpc","pT","eta","phi","pid"};
  TString axisTitleChi2Tpc[kNumberOfAxes] = {"chi2tpc","pT","eta","phi","pid"};
  //
  THnF * histChi2Tpc = new THnF("histChi2Tpc","chi2 per cls. in TPC",kNumberOfAxes, binsChi2Tpc, minChi2Tpc, maxChi2Tpc);
  BinLogAxis(histChi2Tpc, 1);
  fListHist->Add(histChi2Tpc);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histChi2Tpc->GetAxis(iaxis)->SetName(axisNameChi2Tpc[iaxis]);
    histChi2Tpc->GetAxis(iaxis)->SetTitle(axisTitleChi2Tpc[iaxis]);
  }
  //
  // (3.) dca_z
  //                               0-dcaZ, 1-pt, 2-eta,         3-phi, 4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all)
  Int_t    binsDcaZ[kNumberOfAxes] = { 20,   50,    20,            18,  6};
  Double_t minDcaZ[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  Double_t maxDcaZ[kNumberOfAxes]  = {  4,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameDcaZ[kNumberOfAxes]  = {"dcaZ","pT","eta","phi","pid"};
  TString axisTitleDcaZ[kNumberOfAxes] = {"dcaZ","pT","eta","phi","pid"};
  //
  THnF * histDcaZ = new THnF("histDcaZ","dca_z to prim. vtx.",kNumberOfAxes, binsDcaZ, minDcaZ, maxDcaZ);
  BinLogAxis(histDcaZ, 1);
  fListHist->Add(histDcaZ);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histDcaZ->GetAxis(iaxis)->SetName(axisNameDcaZ[iaxis]);
    histDcaZ->GetAxis(iaxis)->SetTitle(axisTitleDcaZ[iaxis]);
  }
  //
  // (4.) hit in SPD layer
  //                              0-spdHit, 1-pt, 2-eta,         3-phi, 4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all)
  Int_t    binsSpd[kNumberOfAxes] = {    2,   50,    20,            18,  6};
  Double_t minSpd[kNumberOfAxes]  = { -0.5,  0.1,    -1,             0, -0.5};
  Double_t maxSpd[kNumberOfAxes]  = {  1.5,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameSpd[kNumberOfAxes]  = {"spdHit","pT","eta","phi","pid"};
  TString axisTitleSpd[kNumberOfAxes] = {"spdHit","pT","eta","phi","pid"};
  //
  THnF * histSpd = new THnF("histSpd","hit in SPD layer or not",kNumberOfAxes, binsSpd, minSpd, maxSpd);
  BinLogAxis(histSpd, 1);
  fListHist->Add(histSpd);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histSpd->GetAxis(iaxis)->SetName(axisNameSpd[iaxis]);
    histSpd->GetAxis(iaxis)->SetTitle(axisTitleSpd[iaxis]);
  }
  //
  // (5.) number of crossed rows       
  //                               0-ncr, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsNcr[kNumberOfAxes] = { 40,   50,    20,            18,  6};
  Double_t minNcr[kNumberOfAxes]  = { 20,  0.1,    -1,             0, -0.5};
  Double_t maxNcr[kNumberOfAxes]  = {160,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameNcr[kNumberOfAxes]  = {"Ncr","pT","eta","phi","pid"};
  TString axisTitleNcr[kNumberOfAxes] = {"Ncr","pT","eta","phi","pid"};
  //
  THnF * histNcr = new THnF("histNcr","number of crossed rows TPC histogram",kNumberOfAxes, binsNcr, minNcr, maxNcr);
  BinLogAxis(histNcr, 1);
  fListHist->Add(histNcr);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histNcr->GetAxis(iaxis)->SetName(axisNameNcr[iaxis]);
    histNcr->GetAxis(iaxis)->SetTitle(axisTitleNcr[iaxis]);
  }
  //
  // (6.) ratio crossed rows over findable clusters       
  //                                0-CRoverFC, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsCRoverFC[kNumberOfAxes] = {  26,   50,    20,            18,  6};
  Double_t minCRoverFC[kNumberOfAxes]  = { 0.4,  0.1,    -1,             0, -0.5};
  Double_t maxCRoverFC[kNumberOfAxes]  = { 1.8,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameCRoverFC[kNumberOfAxes]  = {"CRoverFC","pT","eta","phi","pid"};
  TString axisTitleCRoverFC[kNumberOfAxes] = {"CRoverFC","pT","eta","phi","pid"};
  //
  THnF * histCRoverFC = new THnF("histCRoverFC","number of crossed rows over findable clusters histogram",kNumberOfAxes, binsCRoverFC, minCRoverFC, maxCRoverFC);
  BinLogAxis(histCRoverFC, 1);
  fListHist->Add(histCRoverFC);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histCRoverFC->GetAxis(iaxis)->SetName(axisNameCRoverFC[iaxis]);
    histCRoverFC->GetAxis(iaxis)->SetTitle(axisTitleCRoverFC[iaxis]);
  }
  //
  // (7.) max chi2 / ITS cluster       
  //                               0-Chi2Its, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsChi2Its[kNumberOfAxes] = { 25,   50,    20,            18,  6};
  Double_t minChi2Its[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  Double_t maxChi2Its[kNumberOfAxes]  = { 50,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameChi2Its[kNumberOfAxes]  = {"Chi2Its","pT","eta","phi","pid"};
  TString axisTitleChi2Its[kNumberOfAxes] = {"Chi2Its","pT","eta","phi","pid"};
  //
  THnF * histChi2Its = new THnF("histChi2Its","number of crossed rows TPC histogram",kNumberOfAxes, binsChi2Its, minChi2Its, maxChi2Its);
  BinLogAxis(histChi2Its, 1);
  fListHist->Add(histChi2Its);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histChi2Its->GetAxis(iaxis)->SetName(axisNameChi2Its[iaxis]);
    histChi2Its->GetAxis(iaxis)->SetTitle(axisTitleChi2Its[iaxis]);
  }
  //
  // (8.) tpc active volume length       
  //                                0-TpcLength, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsTpcLength[kNumberOfAxes] = {  40,   50,    20,            18,  6};
  Double_t minTpcLength[kNumberOfAxes]  = {   0,  0.1,    -1,             0, -0.5};
  Double_t maxTpcLength[kNumberOfAxes]  = { 170,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameTpcLength[kNumberOfAxes]  = {"TpcLength","pT","eta","phi","pid"};
  TString axisTitleTpcLength[kNumberOfAxes] = {"TpcLength","pT","eta","phi","pid"};
  //
  THnF * histTpcLength = new THnF("histTpcLength","number of crossed rows TPC histogram",kNumberOfAxes, binsTpcLength, minTpcLength, maxTpcLength);
  BinLogAxis(histTpcLength, 1);
  fListHist->Add(histTpcLength);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histTpcLength->GetAxis(iaxis)->SetName(axisNameTpcLength[iaxis]);
    histTpcLength->GetAxis(iaxis)->SetTitle(axisTitleTpcLength[iaxis]);
  }
  //
  // (9.) match TPC->ITS
  //                                  0-is matched, 1-pt, 2-eta,         3-phi, 4-pid(0-3 -> electron-proton, 4 -> undef, 5 -> all)
  Int_t    binsTpcItsMatch[kNumberOfAxes] = {    2,   50,    20,            18,  6};
  Double_t minTpcItsMatch[kNumberOfAxes]  = { -0.5,  0.1,    -1,             0, -0.5};
  Double_t maxTpcItsMatch[kNumberOfAxes]  = {  1.5,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameTpcItsMatch[kNumberOfAxes]  = {"isMatched","pT","eta","phi","pid"};
  TString axisTitleTpcItsMatch[kNumberOfAxes] = {"isMatched","pT","eta","phi","pid"};
  //
  THnF * histTpcItsMatch = new THnF("histTpcItsMatch","TPC -> ITS matching",kNumberOfAxes, binsTpcItsMatch, minTpcItsMatch, maxTpcItsMatch);
  BinLogAxis(histTpcItsMatch, 1);
  fListHist->Add(histTpcItsMatch);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histTpcItsMatch->GetAxis(iaxis)->SetName(axisNameTpcItsMatch[iaxis]);
    histTpcItsMatch->GetAxis(iaxis)->SetTitle(axisTitleTpcItsMatch[iaxis]);
  }




}

