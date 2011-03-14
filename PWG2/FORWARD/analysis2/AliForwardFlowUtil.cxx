//
// Class used to handle the input from AODs and put it into histograms
// the Forward Flow tasks can run on.
// It can also add flow to AliAODMCParticles. 
//
#include <iostream>
#include "AliForwardFlowUtil.h"
#include "AliAODCentralMult.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliLog.h"
#include "AliAODForwardMult.h"
#include "AliAODEvent.h"

ClassImp(AliForwardFlowUtil)

//_____________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil() : 
  fList(0),
  fZvertex(0) {} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil(TList* outputList) :
  fList(0),
  fZvertex(0)
{ 
  // 
  // Constructor
  //
  // Parameters:
  // TList: list containing histograms for flow analysis
  //
  fList = outputList;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::AODCheck(const AliAODForwardMult* aodfm) const
{
  // 
  // Function to check that and AOD event meets the cuts
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  if (!aodfm->IsTriggerBits(AliAODForwardMult::kInel)) return kFALSE;
  if (!aodfm->HasIpZ()) return kFALSE;
  if (!aodfm->InRange(-fZvertex,fZvertex)) return kFALSE;

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODFMD(const AliAODEvent* aodevent) const
{
  //
  // Loop over AliAODFowardMult object and fill histograms provided
  // by flow task.
  //
  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
  if (!aodfmult) return kFALSE;
  if (!AODCheck(aodfmult)) return kFALSE;

  TH2D hdNdetadphi = aodfmult->GetHistogram();
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSE");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSE");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t mult = 0;
  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  for (Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    for (Int_t phiBin = 1; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      dNdetadphi->Fill(eta, phi, weight);
      if (eta < 4.) dNdphi->Fill(phi, weight);
      mult += weight;
    }
  }
  dNdphi->SetBinContent(0, mult);
//  if (aodfmult->HasCentrality()) dNdphi->SetBinContent(dNdphi->GetNbinsX()+1, aodfmult->GetCentrality());
//  else dNdphi->SetBinContent(dNdphi->GetNbinsX()+1, -1);

  dNdphi->SetBinContent(dNdphi->GetNbinsX()+1, GetCentFromMC(aodevent));
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODSPD(const AliAODEvent* aodevent) const
{
  // 
  // Loop over AliAODCentralMult object and fill histograms
  // provided by flow task
  //
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(aodevent->FindListObject("CentralClusters"));
  if (!aodcmult) return kFALSE;

  TH2D hdNdetadphi = aodcmult->GetHistogram();
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSESPD");
  dNdetadphi->Reset();

  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  for (Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    for (Int_t phiBin = 1; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      dNdetadphi->Fill(eta, phi, weight);
    }
  }
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODtrrefHits(const AliAODEvent* aodevent) const 
{
  //
  // Loop over AliAODForwardMult object, get MC track ref information
  // and fill flow histograms
  //

  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("ForwardMC"));
  if (!aodfmult) return kFALSE;
 // if (!AODCheck(aodfmult)) return kFALSE;

  TH2D hdNdetadphi = aodfmult->GetHistogram();
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSETrRef");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSETrRef");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t mult = 0;
  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  for(Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    for(Int_t phiBin = 1; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      dNdetadphi->Fill(eta, phi, weight);
      dNdphi->Fill(phi, weight);
      mult += weight;
    }
  }
  dNdphi->SetBinContent(0, mult);

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODmc(const AliAODEvent* aodevent, 
                                     TString addFlow = "", 
                                     Int_t type = 0, 
                                     Int_t order = 2) const 
{
  // 
  // Loop over AliAODParticle branch and fill flow histograms
  //

  //retreive MC particles from event
  TClonesArray* mcArray = (TClonesArray*)aodevent->FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  Double_t rp = 0;
  if (addFlow.Length() > 1) {
    AliAODMCHeader* header = (AliAODMCHeader*)aodevent->FindListObject(AliAODMCHeader::StdBranchName());
    rp = header->GetReactionPlaneAngle();
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH1D* dNdphi = (TH1D*)fList->FindObject("hdNdphiSEMC");
  dNdphi->Reset();
  dNdetadphi->Reset();
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Double_t mult = 0;
  Double_t weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    weight = 0;
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    if (particle->Eta() > -3.4 && particle->Eta() < 5) {
      // Add flow if it is in the argument
      if (addFlow.Length() > 1) {
        if (addFlow.Contains("pt"))
          weight += AddptFlow(particle->Pt(), type);
        if (addFlow.Contains("pid"))
          weight += AddpidFlow(particle->PdgCode(), type);
        if (addFlow.Contains("eta"))
          weight += AddetaFlow(particle->Eta(), type);
        if (addFlow.Contains("flat"))
          weight = 0.1*type;
        weight *= 2.*TMath::Cos((Double_t)order*(particle->Phi()-rp)); 
      }
      weight += 1;
      
      dNdphi->Fill(particle->Phi(), weight);
      dNdetadphi->Fill(particle->Eta(), particle->Phi(), weight);
      mult += weight;
    }
  }

  dNdphi->SetBinContent(0, mult);

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::AddptFlow(Double_t Pt = 0, Int_t type = 0) const 
{
  //
  // Add pt dependent flow factor
  //
  Double_t weight = 0;

  if (type == 1) weight = 0.125*Pt;
      
  if (type == 2) {
    weight = 0.2*(Pt/2.0);
    if (Pt > 2.0) weight = 0.2;
  }

  if (type == 3) {
      weight = 0.05*(Pt/2.0);
      if (Pt < 2.0) weight = 0.05;
  } 

  if (type == 4) {
      weight = 0.2*(Pt/1.0);
      if (Pt > 1.0) weight = 0.2;
  }

  if (type == 5) { 
      weight = 0.05*(Pt/1.0);
      if (Pt < 1.0) weight = 0.05;
  }                                                      

  if (type == 6) {
      weight = 0.2*(Pt/3.0);
      if (Pt > 3.0) weight = 0.2;
  }

  return weight;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::AddpidFlow(Int_t ID = 0, Int_t type = 0) const 
{
  //
  // Add pid dependent flow factor 
  //
  Double_t weight = 0;

  if (type == 1) {
    weight = 0.07;
    if (TMath::Abs(ID) ==  211) // pion flow
      weight = 0.1;
    if (TMath::Abs(ID) ==  2212) // proton flow
      weight = 0.05;
  }

  if (type == 2) {
    weight = 0.06;
    if (TMath::Abs(ID) ==  211) // pion flow
      weight = 0.1;
    if (TMath::Abs(ID) ==  2212) // proton flow
      weight = 0.08;
  }

  if (type == 3) {
    weight = 0.05;
    if (TMath::Abs(ID) ==  211) // pion flow
      weight = 0.1;
    if (TMath::Abs(ID) ==  2212) // proton flow
      weight = 0.07;
  }

  if (type == 4) {
    weight = 0.07;
    if (TMath::Abs(ID) ==  211) // pion flow
      weight = 0.1;
    if (TMath::Abs(ID) ==  2212) // proton flow
      weight = 0.085;
  }

  return weight;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::AddetaFlow(Double_t eta = 0, Int_t type = 0) const 
{
  //
  // Add eta dependent flow factor 
  //
  Double_t weight = 0;
  
  if (type == 1) weight = 0.03 + 0.07 * (1 - TMath::Abs(eta) / 6.);

  if (type == 2) weight = 0.07 * (1 - TMath::Abs(eta) / 6.);
  
  if (type == 3) weight = 0.07 * (1 - TMath::Abs(eta) / 8.);
  
  if (type == 4) weight = 0.07 * (1 - TMath::Abs(eta) / 10.);
  
  if (type == 5) weight = 0.07 * (1 - TMath::Abs(eta) / 12.); 

  if (type == 6) weight = 0.07 * (1 - TMath::Abs(eta) / 14.); 

  if (type == 7) weight = 0.07 * (1 - TMath::Abs(eta) / 16.); 

  return weight;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::GetCentFromMC(const AliAODEvent* aodevent) const
{
  //
  // Get centrality from MC impact parameter.
  // Values taken from: https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentStudies
  //
  Double_t cent = -1.;
  Double_t b = -1.;
  AliAODMCHeader* header = (AliAODMCHeader*)aodevent->FindListObject(AliAODMCHeader::StdBranchName());
  if (!header) return cent;
  b = header->GetImpactParameter();
  
  if ( 0.00 <= b && b <  3.50)
    cent = 2.5;
  if ( 3.50 <= b && b <  4.95)
    cent = 7.5;
  if ( 4.95 <= b && b <  6.98)
    cent = 15.;
  if ( 6.98 <= b && b <  8.55)
    cent = 25.;
  if ( 8.55 <= b && b <  9.88)
    cent = 35.;
  if ( 9.88 <= b && b < 11.04)
    cent = 45.;
  if (11.04 <= b && b < 12.09)
    cent = 55.;
  if (12.09 <= b && b < 13.06)
    cent = 65.;
  if (13.06 <= b && b < 13.97)
    cent = 75.;
  if (13.97 <= b && b < 19.61)
    cent = 90.; 

  return cent;
}
//_____________________________________________________________________
//
//
// EOF

