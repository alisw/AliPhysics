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
const Int_t kNumberOfAxes = 5;

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
    fMatchChi()
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
    fMatchChi()
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
  // basic QA and statistics histograms
  //
  TH2F * histVertexSelection = new TH2F("histVertexSelection", "vertex selection; vertex z (cm); accepted/rejected", 100, -50., 50., 2, -0.5, 1.5);
  fListHist->Add(histVertexSelection);
  //
  // track cut variation histograms
  InitializeTrackCutHistograms();
  //
  //
  // matching histograms
  //
  const int nbPt=40;
  const double ptMax=5;
  //
  TH2F * hNMatch    = new TH2F("hNMatch","N Matches",nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5);
  TH2F * hBestMatch = new TH2F("hBestMatch","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatch  = new TH2F("hAllMatch","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatchGlo  = new TH2F("hAllMatchGlo","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  fListHist->Add(hNMatch);
  fListHist->Add(hBestMatch);
  fListHist->Add(hAllMatch);
  fListHist->Add(hAllMatchGlo);
  //
  TH2F * hNMatchBg    = new TH2F("hNMatchBg","N Matches",nbPt,0,ptMax,kMaxMatch+1,-0.5,kMaxMatch+0.5);
  TH2F * hBestMatchBg = new TH2F("hBestMatchBg","Best Match Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatchBg  = new TH2F("hAllMatchBg","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  TH2F * hAllMatchGloBg  = new TH2F("hAllMatchGloBg","All Matches Chi2",nbPt,0,ptMax,2*int(TMath::Max(1.1,kMaxChi2)),0,kMaxChi2);
  fListHist->Add(hNMatchBg);
  fListHist->Add(hBestMatchBg);
  fListHist->Add(hAllMatchBg);
  fListHist->Add(hAllMatchGloBg);
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
  TH2F * hNMatch      = (TH2F*) fListHist->FindObject("hNMatch");
  TH2F * hBestMatch   = (TH2F*) fListHist->FindObject("hBestMatch");
  TH2F * hAllMatch    = (TH2F*) fListHist->FindObject("hAllMatch");
  TH2F * hAllMatchGlo = (TH2F*) fListHist->FindObject("hAllMatchGlo");  
  //
  TH2F * hNMatchBg      = (TH2F*) fListHist->FindObject("hNMatchBg");
  TH2F * hBestMatchBg   = (TH2F*) fListHist->FindObject("hBestMatchBg");
  TH2F * hAllMatchBg    = (TH2F*) fListHist->FindObject("hAllMatchBg");
  TH2F * hAllMatchGloBg = (TH2F*) fListHist->FindObject("hAllMatchGloBg");    
  //
  for (int it=0;it<ntr;it++) {
    AliESDtrack* trSA = fESD->GetTrack(it);
    if (!trSA->IsOn(AliESDtrack::kITSpureSA) || !trSA->IsOn(AliESDtrack::kITSrefit)) continue;
    double pt = trSA->Pt();
    //
    Int_t nmatch = 0;
    for (int i=kMaxMatch;i--;) {fMatchChi[i]=0; fMatchTr[i]=0;} // reset array
    for (int it1=0;it1<ntr;it1++) { 
      if (it1==it) continue;
      AliESDtrack* trESD = fESD->GetTrack(it1);
      if (!trESD->IsOn(AliESDtrack::kTPCrefit)) continue;
      Match(trSA,trESD, nmatch);
    }
    //
    hNMatch->Fill(pt,nmatch);
    if (nmatch>0) hBestMatch->Fill(pt,fMatchChi[0]);
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
      Match(trSA,trESD, nmatch, TMath::Pi());
    }
    //
    hNMatchBg->Fill(pt,nmatch);
    if (nmatch>0) hBestMatchBg->Fill(pt,fMatchChi[0]);
    for (int imt=nmatch;imt--;) {
      hAllMatchBg->Fill(pt,fMatchChi[imt]);
      if (fMatchTr[imt]->IsOn(AliESDtrack::kITSrefit)) hAllMatchGloBg->Fill(pt,fMatchChi[imt]);
    }
    //
  }


}


void AliAnalysisTrackingUncertainties::Match(const AliESDtrack* tr0, const AliESDtrack* tr1, Int_t &nmatch, Double_t rotate) {
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
  if (chi2>kMaxChi2) return;
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
  THnF * histNcl     = (THnF *) fListHist->FindObject("histNcl");
  THnF * histChi2Tpc = (THnF *) fListHist->FindObject("histChi2Tpc");
  THnF * histDcaZ    = (THnF *) fListHist->FindObject("histDcaZ");
  THnF * histSpd     = (THnF *) fListHist->FindObject("histSpd");
  //
  Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
  //
  for (Int_t i=0;i<fESD->GetNumberOfTracks();++i) {
    //
    AliESDtrack *track =fESD->GetTrack(i);
    //
    // relevant variables
    //
    Int_t nclsTPC = track->GetTPCncls();
    Float_t pT  = track->Pt();
    Float_t eta = track->Eta();
    Float_t phi = track->Phi();
    Float_t chi2TPC = track->GetTPCchi2();
    if (nclsTPC != 0) {
      chi2TPC /= nclsTPC; 
    } else {
      chi2TPC = 999.;
    }
    //
    track->GetImpactParameters(dca, cov);
    //
    // (1.) fill number of clusters histogram
    //
    Int_t minNclsTPC = fESDtrackCuts->GetMinNClusterTPC();
    fESDtrackCuts->SetMinNClustersTPC(0);
    if (fESDtrackCuts->AcceptTrack(track)) {
      Double_t vecHistNcl[kNumberOfAxes] = {nclsTPC, pT, eta, phi, 0.};
      histNcl->Fill(vecHistNcl);
    }
    fESDtrackCuts->SetMinNClustersTPC(minNclsTPC);
    //
    // (2.) fill chi2 TPC histogram
    //
    Float_t maxChi2 = fESDtrackCuts->GetMaxChi2PerClusterTPC();
    fESDtrackCuts->SetMaxChi2PerClusterTPC(999.);
    if (fESDtrackCuts->AcceptTrack(track)) {
      Double_t vecHistChi2Tpc[kNumberOfAxes] = {chi2TPC, pT, eta, phi, 0.};
      histChi2Tpc->Fill(vecHistChi2Tpc);
    }
    fESDtrackCuts->SetMaxChi2PerClusterTPC(maxChi2);
    //
    // (3.) fill dca_z histogram
    //
    Float_t maxDcaZ = fESDtrackCuts->GetMaxDCAToVertexZ();
    fESDtrackCuts->SetMaxDCAToVertexZ(999.);
    if (fESDtrackCuts->AcceptTrack(track)) {
      Double_t vecHistDcaZ[kNumberOfAxes] = {TMath::Abs(dca[1]), pT, eta, phi, 0.};
      histDcaZ->Fill(vecHistDcaZ);
    }
    fESDtrackCuts->SetMaxDCAToVertexZ(maxDcaZ);
    //
    // (4.) fill hit in SPD histogram
    //
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kOff);
    if (fESDtrackCuts->AcceptTrack(track)) {
      Int_t hasPoint = 0;
      if (track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) hasPoint = 1;
      Double_t vecHistSpd[kNumberOfAxes] = {hasPoint, pT, eta, phi, 0.};
      histSpd->Fill(vecHistSpd);
    }
    fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);


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
    TString vtxTyp = vertex->GetTitle();
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) vertexOkay = kFALSE;
  }
  else {
    vertexOkay = kTRUE;
  }

  vertexZ = vertex->GetZ();  
  return vertexOkay;

}


//________________________________________________________________________
void AliAnalysisTrackingUncertainties::InitializeTrackCutHistograms() {
  //
  // create histograms for the track cut studies
  //
  //
  // (1.) number of clusters       
  //                               0-ncl, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
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
  //                                  0-chi2, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsChi2Tpc[kNumberOfAxes] = { 40,   50,    20,            18,  6};
  Double_t minChi2Tpc[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  Double_t maxChi2Tpc[kNumberOfAxes]  = {  8,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameChi2Tpc[kNumberOfAxes]  = {"chi2tpc","pT","eta","phi","pid"};
  TString axisTitleChi2Tpc[kNumberOfAxes] = {"chi2tpc","pT","eta","phi","pid"};
  //
  THnF * histChi2Tpc = new THnF("histChi2Tpc","number of clusters histogram",kNumberOfAxes, binsChi2Tpc, minChi2Tpc, maxChi2Tpc);
  BinLogAxis(histChi2Tpc, 1);
  fListHist->Add(histChi2Tpc);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histChi2Tpc->GetAxis(iaxis)->SetName(axisNameChi2Tpc[iaxis]);
    histChi2Tpc->GetAxis(iaxis)->SetTitle(axisTitleChi2Tpc[iaxis]);
  }
  //
  // (3.) dca_z
  //                               0-dcaZ, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsDcaZ[kNumberOfAxes] = { 20,   50,    20,            18,  6};
  Double_t minDcaZ[kNumberOfAxes]  = {  0,  0.1,    -1,             0, -0.5};
  Double_t maxDcaZ[kNumberOfAxes]  = {  4,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameDcaZ[kNumberOfAxes]  = {"dcaZ","pT","eta","phi","pid"};
  TString axisTitleDcaZ[kNumberOfAxes] = {"dcaZ","pT","eta","phi","pid"};
  //
  THnF * histDcaZ = new THnF("histDcaZ","number of clusters histogram",kNumberOfAxes, binsDcaZ, minDcaZ, maxDcaZ);
  BinLogAxis(histDcaZ, 1);
  fListHist->Add(histDcaZ);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histDcaZ->GetAxis(iaxis)->SetName(axisNameDcaZ[iaxis]);
    histDcaZ->GetAxis(iaxis)->SetTitle(axisTitleDcaZ[iaxis]);
  }
  //
  // (4.) hit in SPD layer
  //                              0-spdHit, 1-pt, 2-eta,         3-phi, 4-pid(0,unid,etc.)
  Int_t    binsSpd[kNumberOfAxes] = {    2,   50,    20,            18,  6};
  Double_t minSpd[kNumberOfAxes]  = { -0.5,  0.1,    -1,             0, -0.5};
  Double_t maxSpd[kNumberOfAxes]  = {  1.5,   20,    +1, 2*TMath::Pi(),  5.5};
  //
  TString axisNameSpd[kNumberOfAxes]  = {"spdHit","pT","eta","phi","pid"};
  TString axisTitleSpd[kNumberOfAxes] = {"spdHit","pT","eta","phi","pid"};
  //
  THnF * histSpd = new THnF("histSpd","number of clusters histogram",kNumberOfAxes, binsSpd, minSpd, maxSpd);
  BinLogAxis(histSpd, 1);
  fListHist->Add(histSpd);
  //
  for (Int_t iaxis=0; iaxis<kNumberOfAxes;iaxis++){
    histSpd->GetAxis(iaxis)->SetName(axisNameSpd[iaxis]);
    histSpd->GetAxis(iaxis)->SetTitle(axisTitleSpd[iaxis]);
  }




}

