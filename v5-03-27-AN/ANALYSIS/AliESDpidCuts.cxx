/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
#include <TCanvas.h>
#include <TClass.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TMath.h>
#include <TIterator.h>
#include <TString.h>
#include <TList.h>

#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDpid.h"

#include "AliESDpidCuts.h"

ClassImp(AliESDpidCuts)

const Int_t AliESDpidCuts::kNcuts = 3;

//_____________________________________________________________________
AliESDpidCuts::AliESDpidCuts(const Char_t *name, const Char_t *title):
    AliAnalysisCuts(name, title)
  , fESDpid(NULL)
  , fTPCsigmaCutRequired(0)
  , fTOFsigmaCutRequired(0)
  , fCutTPCclusterRatio(0.)
  , fMinMomentumTOF(0.5)
  , fHcutStatistics(NULL)
  , fHcutCorrelation(NULL)
{
  //
  // Default constructor
  //
  
  fESDpid = new AliESDpid;
  memset(fCutTPCnSigma, 0, sizeof(Float_t)* AliPID::kSPECIES * 2);
  memset(fCutTOFnSigma, 0, sizeof(Float_t)* AliPID::kSPECIES * 2);

  memset(fHclusterRatio, 0, sizeof(TH1F *) * 2);
  memset(fHnSigmaTPC, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);
  memset(fHnSigmaTOF, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);
}

//_____________________________________________________________________
AliESDpidCuts::AliESDpidCuts(const AliESDpidCuts &ref):
    AliAnalysisCuts(ref)
  , fESDpid(NULL)
  , fTPCsigmaCutRequired(ref.fTPCsigmaCutRequired)
  , fTOFsigmaCutRequired(ref.fTOFsigmaCutRequired)
  , fCutTPCclusterRatio(ref.fCutTPCclusterRatio)
  , fMinMomentumTOF(ref.fMinMomentumTOF)
  , fHcutStatistics(NULL)
  , fHcutCorrelation(NULL)
{
  //
  // Copy constructor
  //
  fESDpid = new AliESDpid(*ref.fESDpid);
  memcpy(fCutTPCnSigma, ref.fCutTPCnSigma, sizeof(Float_t) * AliPID::kSPECIES * 2);
  memcpy(fCutTOFnSigma, ref.fCutTOFnSigma, sizeof(Float_t) * AliPID::kSPECIES * 2);
  
  if(ref.fHcutStatistics) fHcutStatistics = dynamic_cast<TH1I *>(ref.fHcutStatistics->Clone());
  if(ref.fHcutCorrelation) fHcutCorrelation = dynamic_cast<TH2I *>(ref.fHcutCorrelation->Clone());
  for(Int_t imode = 0; imode < 2; imode++){
    if(ref.fHclusterRatio[imode]) fHclusterRatio[imode] = dynamic_cast<TH1F *>(ref.fHclusterRatio[imode]->Clone());
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      if(fHnSigmaTPC[ispec][imode]) fHnSigmaTPC[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTPC[ispec][imode]->Clone());
      if(fHnSigmaTOF[ispec][imode]) fHnSigmaTOF[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTPC[ispec][imode]->Clone());
    }
  }
}

//_____________________________________________________________________
AliESDpidCuts &AliESDpidCuts::operator=(const AliESDpidCuts &ref){
  //
  // Assignment operator
  //
  if(this != &ref)
    ref.Copy(*this);
  return *this;
}

//_____________________________________________________________________
AliESDpidCuts::~AliESDpidCuts(){
  //
  // Destructor
  //
  delete fESDpid;

  delete fHcutStatistics;
  delete fHcutCorrelation;
  for(Int_t imode = 0; imode < 2; imode++){
    delete fHclusterRatio[imode];
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      delete fHnSigmaTPC[ispec][imode];
      delete fHnSigmaTOF[ispec][imode];
    }
  }
}

//_____________________________________________________________________
Bool_t AliESDpidCuts::IsSelected(TObject *obj){
  //
  // Select Track
  // 
  AliESDtrack * trk = dynamic_cast<AliESDtrack*>(obj);
  if(!trk){
    AliError("Provided object is not AliESDtrack!");
    return kFALSE;
  }
  const AliESDEvent* evt = trk->GetESDEvent();
  if(!evt){
    AliError("No AliESDEvent!");
    return kFALSE;
  }
  return AcceptTrack(trk, evt);
}

//_____________________________________________________________________
void AliESDpidCuts::Copy(TObject &c) const {
  //
  // Copy function
  //
  AliESDpidCuts &target = dynamic_cast<AliESDpidCuts &>(c);

  target.fESDpid = new AliESDpid(*fESDpid);

  target.fCutTPCclusterRatio = fCutTPCclusterRatio;
  target.fMinMomentumTOF = fMinMomentumTOF;

  target.fTPCsigmaCutRequired = fTPCsigmaCutRequired;
  target.fTOFsigmaCutRequired = fTOFsigmaCutRequired;
  
  if(fHcutStatistics) target.fHcutStatistics = dynamic_cast<TH1I *>(fHcutStatistics->Clone());
  if(fHcutCorrelation) target.fHcutCorrelation = dynamic_cast<TH2I *>(fHcutCorrelation->Clone());
  for(Int_t imode = 0; imode < 2; imode++){
    if(fHclusterRatio[imode]) target.fHclusterRatio[imode] = dynamic_cast<TH1F *>(fHclusterRatio[imode]->Clone());
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      if(fHnSigmaTPC[ispec][imode]) target.fHnSigmaTPC[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTPC[ispec][imode]->Clone());
      if(fHnSigmaTOF[ispec][imode]) target.fHnSigmaTOF[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTOF[ispec][imode]->Clone());
    }
  }
 
  memcpy(target.fCutTPCnSigma, fCutTPCnSigma, sizeof(Float_t) * AliPID::kSPECIES * 2);
  memcpy(target.fCutTOFnSigma, fCutTOFnSigma, sizeof(Float_t) * AliPID::kSPECIES * 2);
 
  AliESDpidCuts::Copy(c);
}

//_____________________________________________________________________
Long64_t AliESDpidCuts::Merge(TCollection *coll){
  //
  // Merge Cut objects
  //
  if(!coll) return 0;
  if(coll->IsEmpty())   return 1;
  if(!HasHistograms())  return 0;
  
  TIterator *iter = coll->MakeIterator();
  TObject *o = NULL; 
  AliESDpidCuts *ref = NULL;
  Int_t counter = 0;
  while((o = iter->Next())){
    ref = dynamic_cast<AliESDpidCuts *>(o);
    if(!ref) continue;
    if(!ref->HasHistograms()) continue;

    fHcutStatistics->Add(ref->fHcutStatistics);
    fHcutCorrelation->Add(ref->fHcutCorrelation);
    for(Int_t imode = 0; imode < 2; imode++){
      fHclusterRatio[imode]->Add(ref->fHclusterRatio[imode]);
      for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
        fHnSigmaTPC[ispec][imode]->Add(ref->fHnSigmaTPC[ispec][imode]);
        fHnSigmaTOF[ispec][imode]->Add(ref->fHnSigmaTOF[ispec][imode]);
      }
    }
    ++counter;
  }
  return ++counter;
}

//_____________________________________________________________________
void AliESDpidCuts::DefineHistograms(Color_t color){
  //
  // Swich on QA and create the histograms
  //
  SetBit(kHasHistograms, kTRUE);
  fHcutStatistics = new TH1I("fHcutStatistics", "Cut Statistics", kNcuts, 0, kNcuts);
  fHcutStatistics->SetLineColor(color);
  fHcutCorrelation = new TH2I("fHcutCorrelation", "Cut Correlation", kNcuts, 0, kNcuts, kNcuts, 0, kNcuts);
  TString cutname[kNcuts] = {"TPCclusterRatio", "TPC sigma", "TOF sigma"};
  for(Int_t icut = 0; icut < kNcuts; icut++){
    fHcutStatistics->GetXaxis()->SetBinLabel(fHcutStatistics->GetXaxis()->GetFirst() + icut, cutname[icut].Data());
    fHcutCorrelation->GetXaxis()->SetBinLabel(fHcutCorrelation->GetXaxis()->GetFirst() + icut, cutname[icut].Data());
    fHcutCorrelation->GetYaxis()->SetBinLabel(fHcutCorrelation->GetYaxis()->GetFirst() + icut, cutname[icut].Data());
  }
  Char_t hname[256], htitle[256];
  for(Int_t imode = 0; imode < 2; imode++){
    snprintf(hname, 256, "fHclusterRatio%s", imode ? "After" : "Before");
    snprintf(htitle, 256, "TPC cluster Ratio %s cuts;Ratio;Entries", imode ? "after" : "before");
    fHclusterRatio[imode] = new TH1F(hname, htitle, 20, 0., 1.);
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      snprintf(hname, 256, "fHnSigma%sTPC%s", AliPID::ParticleName(ispec), imode ? "after" : "before");
      snprintf(htitle, 256, "TPC sigma for %s %s cuts;sigma;Entries", AliPID::ParticleName(ispec), imode ? "after" : "before");
      fHnSigmaTPC[ispec][imode] = new TH1F(hname, htitle, 200, -10., 10.);
      snprintf(hname, 256, "fHnSigma%sTOF%s", AliPID::ParticleName(ispec), imode ? "after" : "before");
      snprintf(htitle, 256, "TOF sigma for %s %s cuts;sigma;Entries", AliPID::ParticleName(ispec), imode ? "after" : "before");
      fHnSigmaTOF[ispec][imode] = new TH1F(hname, htitle, 200, -10., 10.);
    }
  }
}

//_____________________________________________________________________
Bool_t AliESDpidCuts::AcceptTrack(const AliESDtrack *track, const AliESDEvent *event){
  //
  // Check whether the tracks survived the cuts
  //
  enum{
    kCutClusterRatioTPC,
    kCutNsigmaTPC,
    kCutNsigmaTOF
  };
  Long64_t cutRequired=0, cutFullfiled = 0;
  if(fTOFsigmaCutRequired && event == 0)  {
    AliError("No event pointer. Need event pointer for T0 for TOF cut");
    return (0);
  }
  Double_t clusterRatio = track->GetTPCNclsF() ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t>(track->GetTPCNclsF()) : 1.;
  if(fCutTPCclusterRatio > 0.){
    SETBIT(cutRequired, kCutClusterRatioTPC);
    if(clusterRatio >= fCutTPCclusterRatio) 
      SETBIT(cutFullfiled, kCutClusterRatioTPC);
  }
  // check TPC nSigma cut
  Float_t nsigmaTPC[AliPID::kSPECIES], nsigma;   // need all sigmas for QA plotting
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    nsigmaTPC[ispec] = nsigma = fESDpid->NumberOfSigmasTPC(track,(AliPID::EParticleType)ispec);
    if(!(fTPCsigmaCutRequired & 1 << ispec)) continue;
    SETBIT(cutRequired, kCutNsigmaTPC); // We found at least one species where the n-Sigma Cut is required
    if(nsigma >= fCutTPCnSigma[2*ispec] && nsigma <= fCutTPCnSigma[2*ispec+1]) SETBIT(cutFullfiled, kCutNsigmaTPC);    // Fullfiled for at least one species
  }
  // check TOF nSigma cut
  Float_t nsigmaTOF[AliPID::kSPECIES];    // see above
  Bool_t hasTOFpid = track->GetStatus() & AliESDtrack::kTOFpid; // only apply TOF n-sigma cut when PID Status Bit is set
  Double_t times[AliPID::kSPECIES];
  track->GetIntegratedTimes(times);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    
    if(hasTOFpid && event) nsigmaTOF[ispec] = nsigma = fESDpid->NumberOfSigmasTOF(track,(AliPID::EParticleType)ispec, event->GetT0());
    if(!(fTOFsigmaCutRequired && 1 << ispec)) continue;
    SETBIT(cutRequired, kCutNsigmaTOF);
    if(track->GetOuterParam()->P() >= fMinMomentumTOF){
      if(hasTOFpid && nsigma <= fCutTOFnSigma[2*ispec] && nsigma >= fCutTOFnSigma[2*ispec+1]) SETBIT(cutFullfiled, kCutNsigmaTOF);
    }
  }

  // Fill Histograms before cuts
  if(HasHistograms()){
    fHclusterRatio[0]->Fill(clusterRatio);
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      fHnSigmaTPC[ispec][0]->Fill(nsigmaTPC[ispec]);
      if(hasTOFpid) fHnSigmaTOF[ispec][0]->Fill(nsigmaTOF[ispec]);
    }
  }
  if(cutRequired != cutFullfiled){
    // Fill cut statistics
    if(HasHistograms()){
      for(Int_t icut = 0; icut < kNcuts; icut++){
	if(TESTBIT(cutRequired, icut) && !TESTBIT(cutFullfiled, icut)){
	  // cut not fullfiled
	  fHcutStatistics->Fill(icut);
	  for(Int_t jcut = 0; jcut <= icut; jcut++)
	    if(TESTBIT(cutRequired, jcut) && !TESTBIT(cutFullfiled, jcut)) fHcutCorrelation->Fill(jcut, icut);
	}
      }
    }
    return kFALSE;    // At least one cut is not fullfiled
  }

  // Fill Histograms after cuts
  if(HasHistograms()){
    fHclusterRatio[1]->Fill(clusterRatio);
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      fHnSigmaTPC[ispec][1]->Fill(nsigmaTPC[ispec]);
      if(hasTOFpid) fHnSigmaTOF[ispec][1]->Fill(nsigmaTOF[ispec]);
    }
  }

  return kTRUE;
}

//_____________________________________________________________________
void AliESDpidCuts::SaveHistograms(const Char_t * location){
  //
  // Save the histograms to a file
  //
  if(!HasHistograms()){
    AliError("Histograms not on - Exiting");
    return;
  }
  if(!location) location = GetName();
  gDirectory->mkdir(location);
  gDirectory->cd(location);
  fHcutStatistics->Write();
  fHcutCorrelation->Write();

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");

  gDirectory->cd("before_cuts");
  fHclusterRatio[0]->Write();
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    fHnSigmaTPC[ispec][0]->Write();
    fHnSigmaTOF[ispec][0]->Write();
  }

  gDirectory->cd("../after_cuts");
  fHclusterRatio[1]->Write();
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    fHnSigmaTPC[ispec][1]->Write();
    fHnSigmaTOF[ispec][1]->Write();
  }

  gDirectory->cd("..");
}

//_____________________________________________________________________
void AliESDpidCuts::DrawHistograms(){
  //
  // Draw the Histograms
  //
  TCanvas *stat = new TCanvas("cutStat", "Cut Statistics", 640, 480);
  stat->cd();
  fHcutStatistics->SetStats(kFALSE);
  fHcutStatistics->Draw();
  stat->SaveAs(Form("%s_%s.gif", GetName(), stat->GetName()));

  TCanvas *correl = new TCanvas("cutCorrelation", "Cut Correlation", 640, 480);
  correl->cd();
  fHcutCorrelation->SetStats(kFALSE);
  fHcutCorrelation->Draw("colz");
  correl->SaveAs(Form("%s_%s.gif", GetName(), correl->GetName()));

  TCanvas *cRatio = new TCanvas("ClusterRatioTPC", "TPC cluster Ratio", 640, 480);
  cRatio->cd();
  fHclusterRatio[0]->SetLineColor(kRed);
  fHclusterRatio[0]->SetStats(kFALSE);
  fHclusterRatio[0]->Draw();
  fHclusterRatio[1]->SetLineColor(kBlue);
  fHclusterRatio[1]->SetStats(kFALSE);
  fHclusterRatio[1]->Draw("same");
  cRatio->SaveAs(Form("%s_%s.gif",  GetName(), cRatio->GetName()));

  TCanvas *cNsigTPC = new TCanvas("NsigmaTPC", "TPC n-sigma", 640, 480);
  cNsigTPC->Divide(3,2);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    cNsigTPC->cd(ispec + 1);
    fHnSigmaTPC[ispec][0]->SetLineColor(kRed);
    fHnSigmaTPC[ispec][0]->SetStats(kFALSE);
    fHnSigmaTPC[ispec][0]->Draw();
    fHnSigmaTPC[ispec][1]->SetLineColor(kBlue);
    fHnSigmaTPC[ispec][1]->SetStats(kFALSE);
    fHnSigmaTPC[ispec][1]->Draw("same");
  }
  cNsigTPC->SaveAs(Form("%s_%s.gif", GetName(), cNsigTPC->GetName()));

  TCanvas *cNsigTOF = new TCanvas("NsigmaTOF", "TOF n-sigma", 640, 480);
  cNsigTOF->Divide(3,2);
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    cNsigTOF->cd(ispec + 1);
    fHnSigmaTOF[ispec][0]->SetLineColor(kRed);
    fHnSigmaTOF[ispec][0]->SetStats(kFALSE);
    fHnSigmaTOF[ispec][0]->Draw();
    fHnSigmaTOF[ispec][1]->SetLineColor(kBlue);
    fHnSigmaTOF[ispec][1]->SetStats(kFALSE);
    fHnSigmaTOF[ispec][1]->Draw("same");
  }
  cNsigTOF->SaveAs(Form("%s_%s.gif", GetName(), cNsigTOF->GetName()));
}

