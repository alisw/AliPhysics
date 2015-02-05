// $Id$
//
// Jet model task to merge to existing branches
// only implemented for track branches
//
// Author: M. Verweij

#include "AliJetFastSimulation.h"

#include <TClonesArray.h>
#include <TFolder.h>
#include <TLorentzVector.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TProfile.h>
#include <TGrid.h>
#include <TFile.h>
#include <TF1.h>
#include "AliAnalysisManager.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALRecPoint.h"
#include "AliGenerator.h"
#include "AliHeader.h"
#include "AliLog.h"
#include "AliPicoTrack.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliStack.h"
#include "AliVCluster.h"
#include "AliVEvent.h"

ClassImp(AliJetFastSimulation)

//________________________________________________________________________
AliJetFastSimulation::AliJetFastSimulation() : 
AliAnalysisTaskEmcal("AliJetFastSimulation",kTRUE),
  fTracksOutName(""),
  fTracksOut(0x0),
  fNTrackClasses(2),
  fRandom(0),
  fEfficiencyFixed(1.1),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrPtResolutionSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fDiceEfficiencyMinPt(-1.),
  fUncertEfficiency(1),
  fUseTrPtResolutionFromOADB(kFALSE),
  fUseTrEfficiencyFromOADB(kFALSE),
  fPathTrPtResolution(""),
  fPathTrEfficiency(""),
  fHistPtDet(0),
  fh2PtGenPtSmeared(0),
  fp1Efficiency(0),
  fp1PtResolution(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetFastSimulation::AliJetFastSimulation(const char *name) : 
  AliAnalysisTaskEmcal(name,kTRUE),
  fTracksOutName(""),
  fTracksOut(0x0),
  fNTrackClasses(2),
  fRandom(0),
  fEfficiencyFixed(1.1),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrPtResolutionSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fDiceEfficiencyMinPt(-1.),
  fUncertEfficiency(1),
  fUseTrPtResolutionFromOADB(kFALSE),
  fUseTrEfficiencyFromOADB(kFALSE),
  fPathTrPtResolution(""),
  fPathTrEfficiency(""),
  fHistPtDet(0),
  fh2PtGenPtSmeared(0),
  fp1Efficiency(0),
  fp1PtResolution(0)
{
  // Standard constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetFastSimulation::~AliJetFastSimulation()
{
  // Destructor

  delete fRandom;
}

//________________________________________________________________________
void AliJetFastSimulation::ExecOnce() 
{
  // Exec only once.

  AliAnalysisTaskEmcal::ExecOnce();

  if(!fRandom) fRandom = new TRandom3(0);

  if (!fTracksOutName.IsNull()) {
    fTracksOut = new TClonesArray("AliPicoTrack");
    fTracksOut->SetName(fTracksOutName);
    if (InputEvent()->FindListObject(fTracksOutName)) {
      AliFatal(Form("%s: Collection %s is already present in the event!", GetName(), fTracksOutName.Data()));
      return;
    }
    else {
      InputEvent()->AddObject(fTracksOut);
    }
  }
}
//________________________________________________________________________
void AliJetFastSimulation::LocalInit() {
  //initialize track response
  if(fUseTrPtResolutionFromOADB) LoadTrPtResolutionRootFileFromOADB();
  if(fUseTrEfficiencyFromOADB)   LoadTrEfficiencyRootFileFromOADB();
  FitMomentumResolution();
}

//________________________________________________________________________
void AliJetFastSimulation::UserCreateOutputObjects() 
{
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    } else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 1.0;
    }
  }

  fHistPtDet = new TH1F("fHistpt","fHistPtDet;#it{p}_{T};N",nBinPt,binLimitsPt);
  fOutput->Add(fHistPtDet);

  fh2PtGenPtSmeared = new TH2F("fh2PtGenPtSmeared","fh2PtGenPtSmeared",nBinPt,binLimitsPt,nBinPt,binLimitsPt);
  fOutput->Add(fh2PtGenPtSmeared);

  fp1Efficiency = new TProfile("fp1Efficiency","fp1Efficiency",nBinPt,binLimitsPt);
  fOutput->Add(fp1Efficiency);

  fp1PtResolution = new TProfile("fp1PtResolution","fp1PtResolution",nBinPt,binLimitsPt);
  fOutput->Add(fp1PtResolution);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
Bool_t AliJetFastSimulation::Run() 
{
  //Check if information is provided detector level effects
  if(!fMomResH1 && fNTrackClasses>0 ) 
    fUseTrPtResolutionSmearing = 0;
  if(!fMomResH2 && fNTrackClasses>1 ) 
    fUseTrPtResolutionSmearing = 0;
  if(!fMomResH3 && fNTrackClasses>2 ) 
    fUseTrPtResolutionSmearing = 0;

  if(fEfficiencyFixed < 1.)
    fUseDiceEfficiency = 1; // 1 is the default; 2 can be set by user but not implemented
  else {
    if(!fhEffH1 && fNTrackClasses>0 ) 
      fUseDiceEfficiency = 0;
    if(!fhEffH2 && fNTrackClasses>1 ) 
      fUseDiceEfficiency = 0;
    if(!fhEffH3 && fNTrackClasses>2 ) 
      fUseDiceEfficiency = 0;
  }

  fTracksOut->Delete();
  SimulateTracks();
  return kTRUE;
}

//________________________________________________________________________
void AliJetFastSimulation::SimulateTracks()
{
  //Apply toy detector simulation to tracks
  Int_t it = 0;
  const Int_t nTracks = fTracks->GetEntriesFast();
   for (Int_t i = 0; i < nTracks; ++i) {
    AliPicoTrack *picotrack = static_cast<AliPicoTrack*>(fTracks->At(i));
    if (!picotrack)
      continue;

    Bool_t accept = kTRUE;
    Double_t eff[3] = {0};
    Double_t rnd = fRandom->Uniform(1.);  
    if(fUseDiceEfficiency) accept = DiceEfficiency(picotrack,eff,rnd);
    if(!accept) continue;

    AliPicoTrack *track = NULL;
    if(fUseTrPtResolutionSmearing) {
      track = SmearPt(picotrack,eff,rnd);
      (*fTracksOut)[it] = track;
    } else
      track = new ((*fTracksOut)[it]) AliPicoTrack(*picotrack);

    track->SetBit(TObject::kBitMask,1);
    fHistPtDet->Fill(track->Pt());
    it++;
   }
}

//________________________________________________________________________
Bool_t AliJetFastSimulation::DiceEfficiency(AliPicoTrack *vp, Double_t eff[3], Double_t rnd)
{
  // Dice to decide if particle is kept or not - toy  model for efficiency
  //
  Double_t sumEff = 0.;
  Double_t pT = 0.;
 
  if(fEfficiencyFixed<1.)
    sumEff = fEfficiencyFixed;
  else {
    pT = vp->Pt();
    Double_t pTtmp = pT;
    if(pT>10.) pTtmp = 10.;
    if(fhEffH1) eff[0] = fhEffH1->GetBinContent(fhEffH1->FindBin(pTtmp));
    if(fhEffH2) eff[1] = fhEffH2->GetBinContent(fhEffH2->FindBin(pTtmp));
    if(fhEffH3) eff[2] = fhEffH3->GetBinContent(fhEffH3->FindBin(pTtmp));

    sumEff = eff[0]+eff[1]+eff[2];
  }
  if(fUncertEfficiency!=1) sumEff=sumEff+fUncertEfficiency;
  fp1Efficiency->Fill(vp->Pt(),sumEff);
  if(rnd>sumEff && pT > fDiceEfficiencyMinPt) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
AliPicoTrack* AliJetFastSimulation::SmearPt(AliPicoTrack *vp, Double_t eff[3], Double_t rnd)
{
  //Smear momentum of generated particle
  Double_t smear = 1.;
  Double_t  pT = vp->Pt();
  //Select hybrid track category

  //Sort efficiencies from large to small
  Int_t cat[3] = {0};
  TMath::Sort(3,eff,cat);
  if(rnd<=eff[cat[2]])
    smear = GetMomentumSmearing(cat[2],pT);
  else if(rnd<=(eff[cat[2]]+eff[cat[1]]))
    smear = GetMomentumSmearing(cat[1],pT);
  else
    smear = GetMomentumSmearing(cat[0],pT);
  
  fp1PtResolution->Fill(vp->Pt(),smear);

  Double_t sigma = vp->Pt()*smear;
  Double_t pTrec = fRandom->Gaus(vp->Pt(),sigma);
  fh2PtGenPtSmeared->Fill(vp->Pt(),pTrec);

  AliPicoTrack *picotrack = new AliPicoTrack(pTrec,
					     vp->Eta(),
					     vp->Phi(),
					     vp->Charge(),
					     vp->GetLabel(),
					     AliPicoTrack::GetTrackType(vp),
					     vp->GetTrackEtaOnEMCal(),
					     vp->GetTrackPhiOnEMCal(),
					     vp->GetTrackPtOnEMCal(),
					     vp->IsEMCAL(),
					     0.13957); //assume pion mass
    return picotrack;
}

//________________________________________________________________________
Double_t AliJetFastSimulation::GetMomentumSmearing(Int_t cat, Double_t pt) {

  //
  // Get smearing on generated momentum
  //

  TProfile *fMomRes = 0x0;
  if(cat==1 && fMomResH1) fMomRes = (TProfile*)fMomResH1->Clone("fMomRes");
  if(cat==2 && fMomResH2) fMomRes = (TProfile*)fMomResH2->Clone("fMomRes");
  if(cat==3 && fMomResH3) fMomRes = (TProfile*)fMomResH3->Clone("fMomRes");

  if(!fMomRes)
    return 0.;

  Double_t smear = 0.;
  if(pt>20.) {
    if(cat==1 && fMomResH1Fit) smear = fMomResH1Fit->Eval(pt);
    if(cat==2 && fMomResH2Fit) smear = fMomResH2Fit->Eval(pt);
    if(cat==3 && fMomResH3Fit) smear = fMomResH3Fit->Eval(pt);
  }
  else {
    Int_t bin = fMomRes->FindBin(pt);
    smear = fRandom->Gaus(fMomRes->GetBinContent(bin),fMomRes->GetBinError(bin));
  }
  
  if(fMomRes) delete fMomRes;

  return smear;
}

//________________________________________________________________________
void AliJetFastSimulation::LoadTrPtResolutionRootFileFromOADB() {

  if (!gGrid && fPathTrPtResolution.Contains("alien://")) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

  TFile *f = TFile::Open(fPathTrPtResolution.Data());
  if(!f)return;
  TProfile *fProfPtPtSigma1PtGlobSt     = NULL;
  TProfile *fProfPtPtSigma1PtGlobCnoSPD = NULL;
  TProfile *fProfPtPtSigma1PtGlobCnoITS = NULL;
  if(fNTrackClasses>0) fProfPtPtSigma1PtGlobSt     = (TProfile*)f->Get("fProfPtPtSigma1PtGlobSt");
  if(fNTrackClasses>1) fProfPtPtSigma1PtGlobCnoSPD = (TProfile*)f->Get("fProfPtPtSigma1PtGlobCnoSPD");
  if(fNTrackClasses>2) fProfPtPtSigma1PtGlobCnoITS = (TProfile*)f->Get("fProfPtPtSigma1PtGlobCnoITS");

  SetSmearResolution(kTRUE);
  SetMomentumResolutionHybrid(fProfPtPtSigma1PtGlobSt,fProfPtPtSigma1PtGlobCnoSPD,fProfPtPtSigma1PtGlobCnoITS);
}

//________________________________________________________________________
void AliJetFastSimulation::LoadTrEfficiencyRootFileFromOADB() {

   if (!gGrid && fPathTrPtResolution.Contains("alien://")) {
     AliInfo("Trying to connect to AliEn ...");
     TGrid::Connect("alien://");
   }

  TFile *f = TFile::Open(fPathTrEfficiency.Data());
  if(!f)return;
  TH1D *hEffPosGlobSt     = NULL;
  TH1D *hEffPosGlobCnoSPD = NULL;
  TH1D *hEffPosGlobCnoITS = NULL;
  if(fNTrackClasses>0) hEffPosGlobSt = (TH1D*)f->Get("hEffPosGlobSt");
  if(fNTrackClasses>1) hEffPosGlobCnoSPD = (TH1D*)f->Get("hEffPosGlobCnoSPD");
  if(fNTrackClasses>2) hEffPosGlobCnoITS = (TH1D*)f->Get("hEffPosGlobCnoITS");

  SetDiceEfficiency(kTRUE);
  SetEfficiencyHybrid(hEffPosGlobSt,hEffPosGlobCnoSPD,hEffPosGlobCnoITS);
}

//________________________________________________________________________
void AliJetFastSimulation::SetMomentumResolutionHybrid(TProfile *p1, TProfile *p2, TProfile *p3) {
  //
  // set mom res profiles
  //
  if(fMomResH1) delete fMomResH1;
  if(fMomResH2) delete fMomResH2;
  if(fMomResH3) delete fMomResH3;
  
  if(p1) fMomResH1 = new TProfile(*p1);//(TProfile*)p1->Clone("fMomResH1");
  if(p2) fMomResH2 = new TProfile(*p2);//(TProfile*)p2->Clone("fMomResH2");
  if(p3) fMomResH3 = new TProfile(*p3);//(TProfile*)p3->Clone("fMomResH3");
}

//________________________________________________________________________
void AliJetFastSimulation:: SetEfficiencyHybrid(TH1 *h1, TH1 *h2, TH1 *h3) {
  //
  // set tracking efficiency histos
  //
  if(h1) fhEffH1 = (TH1*)h1->Clone("fhEffH1");
  if(h2) fhEffH2 = (TH1*)h2->Clone("fhEffH2");
  if(h3) fhEffH3 = (TH1*)h3->Clone("fhEffH3");
}

//________________________________________________________________________
void AliJetFastSimulation::FitMomentumResolution() {
  //
  // Fit linear function on momentum resolution at high pT
  //

  if(!fMomResH1Fit && fMomResH1) {
    fMomResH1Fit = new TF1("fMomResH1Fit","[0]+[1]*x",0.,200.);
    fMomResH1->Fit(fMomResH1Fit,"LL V0","",5.,30.);
    fMomResH1Fit ->SetRange(5.,100.);
  }

  if(!fMomResH2Fit && fMomResH2) {
    fMomResH2Fit = new TF1("fMomResH2Fit","[0]+[1]*x",0.,200.);
    fMomResH2->Fit(fMomResH2Fit,"LL V0","",5.,30.);
    fMomResH2Fit ->SetRange(5.,100.);
  }

  if(!fMomResH3Fit && fMomResH3) {
    fMomResH3Fit = new TF1("fMomResH3Fit","[0]+[1]*x",0.,200.);
    fMomResH3->Fit(fMomResH3Fit,"LL V0","",5.,30.);
    fMomResH3Fit ->SetRange(5.,100.);
  }

}


