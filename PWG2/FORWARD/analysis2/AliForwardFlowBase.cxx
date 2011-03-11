//
// Class used to handle the input from AODs and put it into histograms
// the Forward Flow tasks can run on
//
#include <iostream>
#include "AliForwardFlowBase.h"
#include "AliAODCentralMult.h"
#include "AliAODMCParticle.h"
#include "AliLog.h"
ClassImp(AliForwardFlowBase)

//_____________________________________________________________________
AliForwardFlowBase::AliForwardFlowBase() : fList(0) {} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardFlowBase::AliForwardFlowBase(TList* OutputList) : fList(0)
{ 
  // 
  // Constructor
  //
  // Parameters:
  //  OutputList: list of histograms for flow analysis
  //
  fList = OutputList;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::AODCheck(AliAODForwardMult* aodfm)
{
  // 
  // Function to check that and AOD event meets the cuts
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  if (!aodfm->IsTriggerBits(AliAODForwardMult::kInel)) return kFALSE;
  if (!aodfm->HasIpZ()) return kFALSE;
  if (!aodfm->InRange(-10,10)) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopAODFMD(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODFowardMult object and fill histograms provided
  // by flow task.
  //
  AliAODForwardMult* AODFMult = static_cast<AliAODForwardMult*>(AODevent->FindListObject("Forward"));
  if (!AODFMult) return kFALSE;
  if (!AODCheck(AODFMult)) return kFALSE;

  // Memory leak starts here!
  TH2D* hdNdetadphi = static_cast<TH2D*>(AODFMult->GetHistogram().Clone("d2ndetadphi"));
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSE");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSE");

  Double_t Mult = 0;
  Double_t Eta = 0;
  Double_t Phi = 0;
  Double_t Weight = 0;
  for(Int_t eta = 1; eta<=hdNdetadphi->GetNbinsX(); eta++) {
    Eta = hdNdetadphi->GetXaxis()->GetBinCenter(eta);
    for(Int_t phi = 1; phi<=hdNdetadphi->GetNbinsY()+1; phi++) {
      Phi    = hdNdetadphi->GetYaxis()->GetBinCenter(phi);
      Weight = hdNdetadphi->GetBinContent(eta, phi);
      dNdetadphi->Fill(Eta, Phi, Weight);
      dNdphi->Fill(Phi, Weight);
      Mult += Weight;
    }
  }
  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE;
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopAODSPD(AliAODEvent* AODevent) 
{
  // 
  // Loop over AliAODCentralMult object and fill histograms
  // provided by flow task
  //
  AliAODCentralMult* AODCMult = static_cast<AliAODCentralMult*>(AODevent->FindListObject("CentralClusters"));
  if (!AODCMult) return kFALSE;

  // Memory leak starts here! 
  TH2D* hdNdetadphi = static_cast<TH2D*>(AODCMult->GetHistogram().Clone("central"));
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSE");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSE");

  Double_t Mult = 0;
  Double_t Eta = 0;
  Double_t Phi = 0;
  Double_t Weight = 0;
  for(Int_t eta = 1; eta<=hdNdetadphi->GetNbinsX(); eta++) {
    Eta = hdNdetadphi->GetXaxis()->GetBinCenter(eta);
    for(Int_t phi = 1; phi<=hdNdetadphi->GetNbinsY()+1; phi++) {
      Phi = hdNdetadphi->GetYaxis()->GetBinCenter(phi);
      Weight = hdNdetadphi->GetBinContent(eta, phi);
      dNdetadphi->Fill(Eta, Phi, Weight);
      dNdphi->Fill(Phi, Weight);
      Mult += Weight;
    }
  }
  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE;
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopAODFMDandSPD(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODForwardMult and AliAODCentralMult onject and fill
  // histograms provided by flow task
  //
  AliAODForwardMult* AODFMult = static_cast<AliAODForwardMult*>(AODevent->FindListObject("Forward"));
  AliAODCentralMult* AODCMult = static_cast<AliAODCentralMult*>(AODevent->FindListObject("CentralClusters"));
  if (!AODCMult || !AODFMult) return kFALSE;
  if (!AODCheck(AODFMult)) return kFALSE;

  TH2D* hdNdetadphiF = static_cast<TH2D*>(AODFMult->GetHistogram().Clone("d2ndetadphi"));
  TH2D* hdNdetadphiC = static_cast<TH2D*>(AODCMult->GetHistogram().Clone("central"));
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSE");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSE");

  Double_t Mult = 0;
  Double_t Eta = 0;
  Double_t Phi = 0;
  Double_t Weight = 0;
  Double_t Forward = 0;
  Double_t Central = 0;
  for(Int_t eta = 1; eta<=hdNdetadphiF->GetNbinsX(); eta++) {
    Eta = hdNdetadphiF->GetXaxis()->GetBinCenter(eta);
    for(Int_t phi = 1; phi<=hdNdetadphiF->GetNbinsY(); phi++) {
      Phi = hdNdetadphiF->GetYaxis()->GetBinCenter(phi);
      Forward = hdNdetadphiF->GetBinContent(eta, phi);
      Central = hdNdetadphiC->GetBinContent(eta, phi);
      if (Forward && Central) Weight = 0.5*Forward + 0.5*Central;
      else Weight = Forward + Central;
      dNdetadphi->Fill(Eta, Phi, Weight);
      dNdphi->Fill(Phi, Weight);
      Mult += Weight;
    }
  }
  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE;
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopAODmc(AliAODEvent* AODevent) 
{
  // 
  // Loop over AliAODParticle branch and fill flow histograms
  //

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)AODevent->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSEMC");
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Float_t Mult = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (TMath::Abs(particle->Eta()) < 1.7) continue;
    if (particle->Phi() < TMath::Pi() / 9.) continue;
    if (particle->Eta() > -3.4 && particle->Eta() < 5) {
      dNdphi->Fill(particle->Phi());
      dNdetadphi->Fill(particle->Eta(), particle->Phi());
      Mult++;
    }
  }

  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE; 

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopAODtrrefHits(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODForwardMult object, get MC track ref information
  // and fill flow histograms
  //

  AliAODForwardMult* AODFMult = static_cast<AliAODForwardMult*>(AODevent->FindListObject("ForwardMC"));
  if (!AODFMult) return kFALSE;
 // if (!AODCheck(AODFMult)) return kFALSE;

  TH2D* hdNdetadphi = static_cast<TH2D*>(AODFMult->GetHistogram().Clone("d2ndetadphiMC"));
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSETrRef");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSETrRef");

  Float_t Mult = 0;
  for(Int_t eta = 1; eta<=hdNdetadphi->GetNbinsX(); eta++) {
    Float_t Eta = hdNdetadphi->GetXaxis()->GetBinCenter(eta);
    for(Int_t phi = 1; phi<=hdNdetadphi->GetNbinsY()+1; phi++) {
      Float_t Phi = hdNdetadphi->GetYaxis()->GetBinCenter(phi);
      Float_t Weight = hdNdetadphi->GetBinContent(eta, phi);
      dNdetadphi->Fill(Eta, Phi, Weight);
      dNdphi->Fill(Phi, Weight);
      Mult += Weight;
    }
  }
  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopMCaddptFlow(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODParticle branch and add pt dependant flow, 
  // fill flow histograms
  //

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)AODevent->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSEMC");
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Double_t Mult = 0;
  Double_t Weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (TMath::Abs(particle->Eta()) < 1.7) continue;
    if (particle->Eta() > -3.4 && particle->Eta() < 5) {
      Weight = 0.5*particle->Pt();
      Weight *= 0.5*TMath::Cos(2*particle->Phi());
      Weight += 1;
      dNdphi->Fill(particle->Phi(), Weight);
      dNdetadphi->Fill(particle->Eta(), particle->Phi(), Weight);
      Mult++;
    }
  }

  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE; 

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopMCaddpdgFlow(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODParticle branch and add pid dependant flow, 
  // fill flow histograms
  //

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)AODevent->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSEMC");
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Double_t Mult = 0;
  Double_t Weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (TMath::Abs(particle->Eta()) < 1.7) continue;
    if (particle->Eta() > -3.4 && particle->Eta() < 5) {
      Weight = 0.07;
      if (TMath::Abs(particle->PdgCode()) ==  211)
        Weight = 0.1;
      if (TMath::Abs(particle->PdgCode()) ==  2212)
        Weight = 0.05;
      Weight *= 2.*TMath::Cos(2*particle->Phi());
      Weight += 1;
      dNdphi->Fill(particle->Phi(), Weight);
      dNdetadphi->Fill(particle->Eta(), particle->Phi(), Weight);
      Mult++;
    }
  }

  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE; 

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowBase::LoopMCaddetaFlow(AliAODEvent* AODevent) 
{
  //
  // Loop over AliAODParticle branch and add eta dependant flow, 
  // fill flow histograms
  //

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)AODevent->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSEMC");
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Double_t Mult = 0;
  Double_t Weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (TMath::Abs(particle->Eta()) < 1.7) continue;
    if (particle->Eta() > -3.4 && particle->Eta() < 5) {
      Weight = 0.03 + 0.07 * (1 - TMath::Abs(particle->Eta()) / 6);
      Weight *= 2.*TMath::Cos(2*particle->Phi());
      Weight += 1;
      dNdphi->Fill(particle->Phi(), Weight);
      dNdetadphi->Fill(particle->Eta(), particle->Phi(), Weight);
      Mult++;
    }
  }

  dNdphi->SetBinContent(0, Mult);
  if (Mult < 100) return kFALSE; 

  return kTRUE;
}
//_____________________________________________________________________
//
//
// EOF
