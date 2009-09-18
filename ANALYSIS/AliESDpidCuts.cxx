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

#include "AliESDtrack.h"
#include "AliLog.h"
#include "AliTPCpidESD.h"
#include "AliTOFpidESD.h"

#include "AliESDpidCuts.h"

ClassImp(AliESDpidCuts)

const Int_t AliESDpidCuts::kNcuts = 3;

//_____________________________________________________________________
AliESDpidCuts::AliESDpidCuts(const Char_t *name, const Char_t *title):
    AliAnalysisCuts(name, title)
  , fTPCpid(NULL)
  , fTOFpid(NULL)
  , fCutTPCclusterRatio(0.)
  , fMinMomentumTOF(0.5)
  , fHcutStatistics(NULL)
  , fHcutCorrelation(NULL)
{
  //
  // Default constructor
  //
  
  fTPCpid = new AliTPCpidESD;
  fTOFpid = new AliTOFpidESD;
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    fCutTPCnSigma[ispec] = -1;
    fCutTOFnSigma[ispec] = -1;
  }
  //memset(fCutTPCnSigma, -1, sizeof(Float_t) * AliPID::kSPECIES);
  //memset(fCutTOFnSigma, -1, sizeof(Float_t) * AliPID::kSPECIES);

  memset(fHclusterRatio, 0, sizeof(TH1F *) * 2);
  memset(fHnSigmaTPC, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);
  memset(fHnSigmaTOF, 0, sizeof(TH1F *) * AliPID::kSPECIES * 2);
}

//_____________________________________________________________________
AliESDpidCuts::AliESDpidCuts(const AliESDpidCuts &ref):
    AliAnalysisCuts(ref)
  , fTPCpid(NULL)
  , fTOFpid(NULL)
  , fCutTPCclusterRatio(ref.fCutTPCclusterRatio)
  , fMinMomentumTOF(ref.fMinMomentumTOF)
  , fHcutStatistics(NULL)
  , fHcutCorrelation(NULL)
{
  //
  // Copy constructor
  //
  fTPCpid = new AliTPCpidESD(*ref.fTPCpid);
  fTOFpid = new AliTOFpidESD(*ref.fTOFpid);
  memcpy(fCutTPCnSigma, ref.fCutTPCnSigma, sizeof(Float_t) * AliPID::kSPECIES);
  memcpy(fCutTOFnSigma, ref.fCutTOFnSigma, sizeof(Float_t) * AliPID::kSPECIES);
  
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
  delete fTPCpid;
  delete fTOFpid;

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
Bool_t AliESDpidCuts::IsSelected(TObject *o){
  //
  // Select Track
  // 
  if(TString(o->IsA()->GetName()).CompareTo("AliESDtrack")){
    Char_t errormessage[256];
    sprintf(errormessage, "Provided object not an AliESDtrack: Type %s", o->IsA()->GetName());
    AliError(errormessage);
    return kFALSE;
  }
  return AcceptTrack(const_cast<const AliESDtrack *>(dynamic_cast<AliESDtrack *>(o)));
}

//_____________________________________________________________________
void AliESDpidCuts::Copy(TObject &c) const {
  //
  // Copy function
  //
  AliESDpidCuts &target = dynamic_cast<AliESDpidCuts &>(c);

  target.fTPCpid = new AliTPCpidESD;
  target.fTOFpid = new AliTOFpidESD;

  target.fCutTPCclusterRatio = fCutTPCclusterRatio;
  target.fMinMomentumTOF = fMinMomentumTOF;
  
  if(fHcutStatistics) target.fHcutStatistics = dynamic_cast<TH1I *>(fHcutStatistics->Clone());
  if(fHcutCorrelation) target.fHcutCorrelation = dynamic_cast<TH2I *>(fHcutCorrelation->Clone());
  for(Int_t imode = 0; imode < 2; imode++){
    if(fHclusterRatio[imode]) target.fHclusterRatio[imode] = dynamic_cast<TH1F *>(fHclusterRatio[imode]->Clone());
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      if(fHnSigmaTPC[ispec][imode]) target.fHnSigmaTPC[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTPC[ispec][imode]->Clone());
      if(fHnSigmaTOF[ispec][imode]) target.fHnSigmaTOF[ispec][imode] = dynamic_cast<TH1F *>(fHnSigmaTOF[ispec][imode]->Clone());
    }
  }
 
  memcpy(target.fCutTPCnSigma, fCutTPCnSigma, sizeof(Float_t) * AliPID::kSPECIES);
  memcpy(target.fCutTOFnSigma, fCutTOFnSigma, sizeof(Float_t) * AliPID::kSPECIES);
 
  AliESDpidCuts::Copy(c);
}

//_____________________________________________________________________
Long64_t AliESDpidCuts::Merge(TCollection *coll){
  //
  // Merge Cut objects
  //
  if(coll) return 0;
  if(coll->IsEmpty()) return 1;
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
    sprintf(hname, "fHclusterRatio%s", imode ? "After" : "Before");
    sprintf(htitle, "TPC cluster Ratio %s cuts;Ratio;Entries", imode ? "after" : "before");
    fHclusterRatio[imode] = new TH1F(hname, htitle, 20, 0., 1.);
    for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      sprintf(hname, "fHnSigma%sTPC%s", AliPID::ParticleName(ispec), imode ? "after" : "before");
      sprintf(htitle, "TPC sigma for %s %s cuts;sigma;Entries", AliPID::ParticleName(ispec), imode ? "after" : "before");
      fHnSigmaTPC[ispec][imode] = new TH1F(hname, htitle, 200, -10., 10.);
      sprintf(hname, "fHnSigma%sTOF%s", AliPID::ParticleName(ispec), imode ? "after" : "before");
      sprintf(htitle, "TOF sigma for %s %s cuts;sigma;Entries", AliPID::ParticleName(ispec), imode ? "after" : "before");
      fHnSigmaTOF[ispec][imode] = new TH1F(hname, htitle, 200, -10., 10.);
    }
  }
}

//_____________________________________________________________________
Bool_t AliESDpidCuts::AcceptTrack(const AliESDtrack *track){
  //
  // Check whether the tracks survived the cuts
  //
  enum{
    kCutClusterRatioTPC,
    kCutNsigmaTPC,
    kCutNsigmaTOF
  };
  Long64_t cutRequired=0, cutFullfiled = 0;
  Double_t clusterRatio = track->GetTPCNclsF() ? static_cast<Float_t>(track->GetTPCNcls())/static_cast<Float_t>(track->GetTPCNclsF()) : 1.;
  if(fCutTPCclusterRatio > 0.){
    SETBIT(cutRequired, kCutClusterRatioTPC);
    if(clusterRatio >= fCutTPCclusterRatio) 
      SETBIT(cutFullfiled, kCutClusterRatioTPC);
  }
  // check TPC nSigma cut
  Float_t nsigmaTPC[AliPID::kSPECIES];   // need all sigmas for QA plotting
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    nsigmaTPC[ispec] = fTPCpid->GetNumberOfSigmas(track, static_cast<AliPID::EParticleType>(ispec));
    if(fCutTPCnSigma[ispec] < 0) continue;
    SETBIT(cutRequired, kCutNsigmaTPC); // We found at least one species where the n-Sigma Cut is required
    if(TMath::Abs(nsigmaTPC[ispec]) <= fCutTPCnSigma[ispec]) SETBIT(cutFullfiled, kCutNsigmaTPC);    // Fullfiled for at least one species
  }
  // check TOF nSigma cut
  Float_t nsigmaTOF[AliPID::kSPECIES];    // see above
  Bool_t hasTOFpid = track->GetStatus() & AliESDtrack::kTOFpid; // only apply TOF n-sigma cut when PID Status Bit is set
  for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
    if(hasTOFpid) nsigmaTOF[ispec] = fTOFpid->GetNumberOfSigmas(track, static_cast<AliPID::EParticleType>(ispec));
    if(fCutTOFnSigma[ispec] < 0) continue;
    SETBIT(cutRequired, kCutNsigmaTOF);
    if(track->GetOuterParam()->P() >= fMinMomentumTOF){
      if(hasTOFpid && (TMath::Abs(nsigmaTOF[ispec]) <= fCutTOFnSigma[ispec])) SETBIT(cutFullfiled, kCutNsigmaTOF);
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

