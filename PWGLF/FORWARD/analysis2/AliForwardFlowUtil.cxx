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
#include "TProfile2D.h"
#include "TGraph.h"
#include "TF1.h"

ClassImp(AliForwardFlowUtil)

//_____________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil() : 
  fList(0),
  fCent(-1),
  fVertex(1111),
  fAliceCent4th(),
  fAlicePt2nd4050(),
  fAlicePt4th3040(),
  fAlicePt4th4050(),
  fImpactParToCent()
   {} 
  //
  // Default Constructor
  //
//_____________________________________________________________________
AliForwardFlowUtil::AliForwardFlowUtil(TList* outputList) :
  fList(0),
  fCent(-1),
  fVertex(1111),
  fAliceCent4th(),
  fAlicePt2nd4050(),
  fAlicePt4th3040(),
  fAlicePt4th4050(),
  fImpactParToCent()
{ 
  // 
  // Constructor
  //
  // Parameters:
  // TList: list containing histograms for flow analysis
  //
  fList = outputList;

  Double_t xCumulant4thTPCrefMultTPConlyAll[] = {2.5,7.5,15,25,35,45,55,65};
  Double_t yCumulant4thTPCrefMultTPConlyAll[] = {0.017855,0.032440,0.055818,0.073137,0.083898,0.086690,0.082040,0.077777};
  Int_t nPointsCumulant4thTPCrefMultTPConlyAll = sizeof(xCumulant4thTPCrefMultTPConlyAll)/sizeof(Double_t);
  fAliceCent4th = new TGraph(nPointsCumulant4thTPCrefMultTPConlyAll,xCumulant4thTPCrefMultTPConlyAll,
                                                        yCumulant4thTPCrefMultTPConlyAll);
 
  Double_t xCumulant2nd4050ALICE[] = {0.000000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
  Double_t yCumulant2nd4050ALICE[] = {0.000000,0.043400,0.059911,0.073516,0.089756,0.105486,0.117391,0.128199,0.138013,
  0.158271,0.177726,0.196383,0.208277,0.216648,0.242954,0.249961,0.240131,0.269006,0.207796};
  Int_t nPointsCumulant2nd4050ALICE = sizeof(xCumulant2nd4050ALICE)/sizeof(Double_t);                                      
  fAlicePt2nd4050 = new TGraph(nPointsCumulant2nd4050ALICE,xCumulant2nd4050ALICE,yCumulant2nd4050ALICE);

  Double_t xCumulant4th3040ALICE[] = {0.00000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000,
  5.500000,7.000000,9.000000};
  Double_t yCumulant4th3040ALICE[] = {0.000000,0.037071,0.048566,0.061083,0.070910,0.078831,0.091396,0.102026,0.109691,
  0.124449,0.139819,0.155561,0.165701,0.173678,0.191149,0.202015,0.204540,0.212560,0.195885,
  0.000000,0.000000,0.000000};
  Int_t nPointsCumulant4th3040ALICE = sizeof(xCumulant4th3040ALICE)/sizeof(Double_t);                                      
  fAlicePt4th3040 = new TGraph(nPointsCumulant4th3040ALICE,xCumulant4th3040ALICE,yCumulant4th3040ALICE);

  Double_t xCumulant4th4050ALICE[] = {0.000000,0.250000,0.350000,0.450000,0.550000,0.650000,0.750000,0.850000,0.950000,
  1.100000,1.300000,1.500000,1.700000,1.900000,2.250000,2.750000,3.250000,3.750000,4.500000};
  Double_t yCumulant4th4050ALICE[] = {0.000000,0.038646,0.049824,0.066662,0.075856,0.081583,0.099778,0.104674,0.118545,
  0.131874,0.152959,0.155348,0.169751,0.179052,0.178532,0.198851,0.185737,0.239901,0.186098};
  Int_t nPointsCumulant4th4050ALICE = sizeof(xCumulant4th4050ALICE)/sizeof(Double_t);   
  fAlicePt4th4050 = new TGraph(nPointsCumulant4th4050ALICE, xCumulant4th4050ALICE, yCumulant4th4050ALICE);

  Double_t impactParam[] = {0.,1.75,4.225,5.965,7.765,9.215,10.46,11.565,12.575,13.515,16.679};
  Double_t centrality[] = {0.,2.5,7.5,15,25,35,45,55,65,75,90};

  Int_t nPoints = sizeof(impactParam)/sizeof(Double_t);
  fImpactParToCent = new TGraph(nPoints, impactParam, centrality);

}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::AODCheck(const AliAODForwardMult* aodfm) 
{
  // 
  // Function to check that and AOD event meets the cuts
  //
  // Parameters: 
  //  AliAODForwardMult: forward mult object with trigger and vertex info
  //
  fCent = -1;
  fVertex = 1111;

  if (!aodfm->IsTriggerBits(AliAODForwardMult::kOffline)) return kFALSE;
  fCent = (Double_t)aodfm->GetCentrality();
  if (0. >= fCent || fCent >= 100.) return kFALSE;
  TH1D* vertex = (TH1D*)fList->FindObject("VertexAll");
  fVertex = aodfm->GetIpZ();
  vertex->Fill(fVertex);
  if (TMath::Abs(fVertex) >= 5.) return kFALSE;
  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODFMD(const AliAODEvent* aodevent)
{
  //
  // Loop over AliAODFowardMult object and fill histograms provided
  // by flow task.
  //

  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("Forward"));
  if (!aodfmult) return kFALSE;
  if (!AODCheck(aodfmult)) return kFALSE;

  TH2D hdNdetadphi = aodfmult->GetHistogram();
  TH2D* vertex = (TH2D*)fList->FindObject("CoverageVsVertex");
  TH2D* dNdphi = (TH2D*)fList->FindObject("hdNdphiSEFMD");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEFMD");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t mult = 0;
  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  for (Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    if (TMath::Abs(eta) < 1.75) continue;
    for (Int_t phiBin = 0; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      if (phiBin == 0) {
        vertex->Fill(eta, fVertex, weight);
        continue;
      }
      dNdetadphi->Fill(eta, phi, weight);
      dNdphi->Fill(eta, phi, weight);
      dNdphi->Fill(-1.*eta, phi, weight);
      mult += weight;
    }
  }
//  fCent = GetCentFromMC(aodevent);
//  fCent = 0.5;
 
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
  TH2D* vertex = (TH2D*)fList->FindObject("CoverageVsVertex");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSESPD");
  TH2D* dNdphi = (TH2D*)fList->FindObject("hdNdphiSESPD");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  Double_t mult = 0;
  for (Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    if (TMath::Abs(eta) > 1.75) continue;
    for (Int_t phiBin = 0; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      if (phiBin == 0) {
        vertex->Fill(eta, fVertex, weight);
        continue;
      }
      dNdetadphi->Fill(eta, phi, weight);
      if (TMath::Abs(eta) < 0.5) continue;
      dNdphi->Fill(eta, phi, weight);
//      dNdphi->Fill(-1.*eta, phi, weight);
      mult += weight;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODFMDtrrefHits(const AliAODEvent* aodevent) const 
{
  //
  // Loop over AliAODForwardMult object, get MC track ref information
  // and fill flow histograms
  //

  AliAODForwardMult* aodfmult = static_cast<AliAODForwardMult*>(aodevent->FindListObject("ForwardMC"));
  if (!aodfmult) return kFALSE;

  TH2D hdNdetadphi = aodfmult->GetHistogram();
  TH2D* dNdphi = (TH2D*)fList->FindObject("hdNdphiSEFMDTR");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEFMDTR");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t mult = 0;
  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  for(Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    if (TMath::Abs(eta) < 1.75) continue;
    for(Int_t phiBin = 1; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      dNdetadphi->Fill(eta, phi, weight);
      dNdphi->Fill(eta, phi, weight);
      dNdphi->Fill(-1.*eta, phi, weight);
      mult += weight;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Bool_t AliForwardFlowUtil::LoopAODSPDtrrefHits(const AliAODEvent* aodevent) const
{
  // 
  // Loop over AliAODCentralMult object and fill histograms
  // provided by flow task
  //
  AliAODCentralMult* aodcmult = static_cast<AliAODCentralMult*>(aodevent->FindListObject("CentralClustersMC"));
  if (!aodcmult) return kFALSE;

  TH2D hdNdetadphi = aodcmult->GetHistogram();
  TH2D* dNdphi = (TH2D*)fList->FindObject("hdNdphiSESPDTR");
  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSESPDTR");
  dNdphi->Reset();
  dNdetadphi->Reset();

  Double_t eta = 0;
  Double_t phi = 0;
  Double_t weight = 0;
  Double_t mult = 0;
  for (Int_t etaBin = 1; etaBin<=hdNdetadphi.GetNbinsX(); etaBin++) {
    eta = hdNdetadphi.GetXaxis()->GetBinCenter(etaBin);
    if (TMath::Abs(eta) > 1.75) continue;
    for (Int_t phiBin = 1; phiBin<=hdNdetadphi.GetNbinsY(); phiBin++) {
      phi = hdNdetadphi.GetYaxis()->GetBinCenter(phiBin);
      weight = hdNdetadphi.GetBinContent(etaBin, phiBin);
      dNdetadphi->Fill(eta, phi, weight);
      dNdphi->Fill(eta, phi, weight);
      dNdphi->Fill(-1.*eta, phi, weight);
      mult += weight;
    }
  }

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
//    AliWarning("No MC array found in AOD. Try making it again.");
    return kFALSE;
  }

  Double_t rp = 0;
  AliAODMCHeader* header = dynamic_cast<AliAODMCHeader*>(aodevent->FindListObject(AliAODMCHeader::StdBranchName()));
  if (!header) {
    AliWarning("No header file found.");
  }
  else {
    rp = header->GetReactionPlaneAngle();
  }

  TH2D* dNdetadphi = (TH2D*)fList->FindObject("hdNdetadphiSEMC");
  TH2D* dNdphi = (TH2D*)fList->FindObject("hdNdphiSEMC");
  TProfile2D* mcTruth = (TProfile2D*)fList->FindObject("pMCTruth");
  dNdphi->Reset();
  dNdetadphi->Reset();
  
  Int_t ntracks = mcArray->GetEntriesFast();

  // Track loop: chek how many particles will be accepted
  Double_t mult = 0;
  Double_t weight = 0;
  for (Int_t it = 0; it < ntracks; it++) {
    weight = 1;
    AliAODMCParticle* particle = (AliAODMCParticle*)mcArray->At(it);
    if (!particle) {
      AliError(Form("Could not receive track %d", it));
      continue;
    }
    if (!particle->IsPhysicalPrimary()) continue;
    if (particle->Charge() == 0) continue;
    Double_t pT = particle->Pt();
//    if (pT <= 0.2 || pT >= 5.) continue;
    Double_t eta = particle->Eta();
    Double_t phi = particle->Phi();
    if (TMath::Abs(eta) < 6.) {
      // Add flow if it is in the argument
      if (addFlow.Length() > 1) {
        if (addFlow.Contains("pt"))
          weight *= AddptFlow(pT, type);
//          weight *= AddptFlow(pT, 1);
        if (addFlow.Contains("pid"))
          weight *= AddpidFlow(particle->PdgCode(), 1);
        if (addFlow.Contains("eta"))
          weight *= AddetaFlow(eta, type);
//          weight *= AddetaFlow(eta, 2);
        if (addFlow.Contains("cent")) 
          weight *= fAliceCent4th->Eval(fCent)/fAliceCent4th->Eval(45);
        
        weight *= 20*2.*TMath::Cos((Double_t)order*(phi-rp)); 
        weight += 1;
      }
      dNdphi->Fill(eta, phi, weight);
      dNdphi->Fill(-1.*eta, phi, weight);
      dNdetadphi->Fill(eta, phi, weight);
      mcTruth->Fill(eta, fCent, weight*TMath::Cos(2.*(phi-rp)));
      mult += weight;
    }
  }

  return kTRUE;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::AddptFlow(Double_t pt = 0, Int_t type = 0) const 
{
  //
  // Add pt dependent flow factor
  //
  Double_t weight = 0;

  switch(type) 
  {
    case 1: weight = fAlicePt2nd4050->Eval(pt)*0.5+fAlicePt4th4050->Eval(pt)*0.5;
            break;
    case 2: weight = fAlicePt2nd4050->Eval(pt);
            break;
    case 3: weight = fAlicePt4th3040->Eval(pt);
            break;
    case 4: weight = fAlicePt4th4050->Eval(pt);
            break;
  } 
  
  return weight;
}
//_____________________________________________________________________
Double_t AliForwardFlowUtil::AddpidFlow(Int_t id = 0, Int_t type = 0) const 
{
  //
  // Add pid dependent flow factor 
  //
  Double_t weight = 0;

  switch(type)
  {
    case 1: if (TMath::Abs(id) ==  211) // pion flow
              weight = 1.;
            else if (TMath::Abs(id) ==  2212) // proton flow
              weight = 1.3;
            else 
              weight = 0.7;
            break;
    case 2: weight = 1.207;
            break;
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

  TF1 gaus = TF1("gaus", "gaus", -6, 6);

  switch(type)
  {
     case 1: gaus.SetParameters(0.1, 0., 9);
             break;
     case 2: gaus.SetParameters(0.1, 0., 3);
             break;
     case 3: gaus.SetParameters(0.1, 0., 15);
             break;
  }

  weight = gaus.Eval(eta);

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

  cent = fImpactParToCent->Eval(b);

  return cent;
}
//_____________________________________________________________________
//
//
// EOF

