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

/* $Id: AliAnalysisTaskPi0V2.cxx 55404 2012-03-29 10:10:19Z fca $ */

/* AliAnalysisTaskPi0V2.cxx
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Instructions for adding histograms can be found below, starting with NEW HISTO
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#include "AliAnalysisTaskPi0V2.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliVCluster.h"

#include "AliEventplane.h"
#include "AliEMCALGeometry.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TString.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPi0V2)

//________________________________________________________________________
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutput(0),
    fESD(0),
    fTracksName("PicoTrack"), fV1ClusName("CaloCluster"), fV2ClusName("CaloCluster"),
    fTrigClass("CVLN_|CSEMI_|CCENT|CVHN"),
    fTracks(0), fV1Clus(0), fV2Clus(0),
    fRunNumber(-999.),
    fEvtSelect(1),
    fVtxCut(15.),
    fNcellCut(2), fECut(1), fEtaCut(0.65), fM02Cut(0.5), fPi0AsyCut(0), isV1Clus(1),
    fCentrality(99.),
    fEPTPC(-999.),
    fEPTPCreso(0.), 
    fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.), fEPV0Ar(-999.), fEPV0Cr(-999.), fEPV0r(-999.),
    fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.), fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
    hEvtCount(0), hAllcentV0(0), hAllcentV0r(0), hAllcentV0A(0), hAllcentV0C(0), hAllcentTPC(0),
    h2DcosV0r(0), h2DsinV0r(0), h2DcosV0A(0), h2DsinV0A(0), h2DcosV0C(0), h2DsinV0C(0), h2DcosTPC(0), h2DsinTPC(0), 
    hEPTPC(0), hresoTPC(0),
    hEPV0(0), hEPV0A(0), hEPV0C(0), hEPV0Ar(0), hEPV0Cr(0), hEPV0r(0), hEPV0AR4(0), hEPV0AR7(0), hEPV0CR0(0), hEPV0CR3(0),
    hdifV0Ar_V0Cr(0), hdifV0A_V0CR0(0), hdifV0A_V0CR3(0), hdifV0ACR0_V0CR3(0), hdifV0C_V0AR4(0), hdifV0C_V0AR7(0), hdifV0AR4_V0AR7(0),
    hdifV0A_V0C(0), hdifV0A_TPC(0), hdifV0C_TPC(0), hdifV0C_V0A(0), 
    hM02vsPtA(0), hM02vsPtB(0), hClusDxDZA(0), hClusDxDZB(0),
    hdifEMC_EPV0(0), hdifEMC_EPV0A(0), hdifEMC_EPV0C(0), hdifful_EPV0(0), hdifful_EPV0A(0), hdifful_EPV0C(0), 
    hdifout_EPV0(0), hdifout_EPV0A(0), hdifout_EPV0C(0), hdifEMC_EPTPC(0), hdifful_EPTPC(0), hdifout_EPTPC(0),
    hdifClus_EPV0(0), hdifClus_EPV0A(0), hdifClus_EPV0C(0), hdifClus_EPTPC(0),
    fHEPV0r(0), fHEPV0A(0), fHEPV0C(0), fHEPTPC(0)

{
    // Dummy constructor ALWAYS needed for I/O.
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2() // All data members should be initialised here
   :AliAnalysisTaskSE("default_name"),
    fOutput(0),
    fESD(0),
    fTracksName("PicoTracks"),
    fTrigClass("CVLN_|CSEMI_|CCENT|CVHN"),
    fTracks(0),
    fRunNumber(-999.),
    fEvtSelect(1),
    fVtxCut(15.),
    fNcellCut(2), fECut(1), fEtaCut(0.65), fM02Cut(0.5), fPi0AsyCut(0), isV1Clus(1),
    fCentrality(99.),
    fEPTPC(-999.),
    fEPTPCreso(0.),
    fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.), fEPV0Ar(-999.), fEPV0Cr(-999.), fEPV0r(-999.),
    fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.), fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
    hEvtCount(0), hAllcentV0(0), hAllcentV0r(0), hAllcentV0A(0), hAllcentV0C(0), hAllcentTPC(0),
    h2DcosV0r(0), h2DsinV0r(0), h2DcosV0A(0), h2DsinV0A(0), h2DcosV0C(0), h2DsinV0C(0), h2DcosTPC(0), h2DsinTPC(0), 
    hEPTPC(0), hresoTPC(0),
    hEPV0(0), hEPV0A(0), hEPV0C(0), hEPV0Ar(0), hEPV0Cr(0), hEPV0r(0), hEPV0AR4(0), hEPV0AR7(0), hEPV0CR0(0), hEPV0CR3(0),
    hdifV0Ar_V0Cr(0), hdifV0A_V0CR0(0), hdifV0A_V0CR3(0), hdifV0ACR0_V0CR3(0), hdifV0C_V0AR4(0), hdifV0C_V0AR7(0), hdifV0AR4_V0AR7(0),
    hdifV0A_V0C(0), hdifV0A_TPC(0), hdifV0C_TPC(0), hdifV0C_V0A(0),  
    hM02vsPtA(0), hM02vsPtB(0), hClusDxDZA(0), hClusDxDZB(0),
    hdifEMC_EPV0(0), hdifEMC_EPV0A(0), hdifEMC_EPV0C(0), hdifful_EPV0(0), hdifful_EPV0A(0), hdifful_EPV0C(0), 
    hdifout_EPV0(0), hdifout_EPV0A(0), hdifout_EPV0C(0), hdifEMC_EPTPC(0), hdifful_EPTPC(0), hdifout_EPTPC(0),
    hdifClus_EPV0(0), hdifClus_EPV0A(0), hdifClus_EPV0C(0), hdifClus_EPTPC(0),
    fHEPV0r(0), fHEPV0A(0), fHEPV0C(0), fHEPTPC(0)
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskPi0V2::~AliAnalysisTaskPi0V2()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
     delete fOutput;
}
//_____________________________________________________________________
Double_t AliAnalysisTaskPi0V2::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = 0;
  if (fESD)
    cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  Double_t maxe = 0;
  const Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}
//_____________________________________________________________________
Double_t AliAnalysisTaskPi0V2::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax) const
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = 0;
  if (fESD)
    cells = fESD->GetEMCALCells();
  if (!cells)
    return 0;

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0.;

  geom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
        (aphidiff==0 && aetadiff==1) ) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}
//_____________________________________________________________________
Bool_t AliAnalysisTaskPi0V2::IsWithinFiducialVolume(Short_t id) const
{
  // Check if cell is within given fiducial volume.

  Double_t fNFiducial = 1;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;

  Bool_t okrow = kFALSE;
  Bool_t okcol = kFALSE;

  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if (!geom)
    return kFALSE;

  Int_t cellAbsId = id;
  geom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
  geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);

  // Check rows/phi
  if (iSupMod < 10) {
    if (iphi >= fNFiducial && iphi < 24-fNFiducial)
      okrow = kTRUE;
  } else {
    if (iphi >= fNFiducial && iphi < 12-fNFiducial)
      okrow = kTRUE;
  }
  // Check columns/eta
  Bool_t noEMCALBorderAtEta0 = kTRUE;
  if (!noEMCALBorderAtEta0) {
    if (ieta > fNFiducial && ieta < 48-fNFiducial)
      okcol = kTRUE;
  } else {
    if (iSupMod%2==0) {
      if (ieta >= fNFiducial)
        okcol = kTRUE;
    } else {
      if (ieta < 48-fNFiducial)
        okcol = kTRUE;
    }
  }
  if (okrow && okcol)
     return kTRUE;

  return kFALSE;
}
//______________________________________________________________________
Bool_t AliAnalysisTaskPi0V2::IsGoodCluster(const AliVCluster *c) const
{

  if(!c)
    return kFALSE;

  if(c->GetNCells() < fNcellCut)
   return kFALSE;

  if(c->E() < fECut)
   return kFALSE;

  Short_t id = -1;
  Double_t maxE = GetMaxCellEnergy(c, id); 
     if((1. - double(GetCrossEnergy(c,id))/maxE) > 0.97)
    return kFALSE;


  Float_t pos1[3];
  c->GetPosition(pos1);
  TVector3 clsPos(pos1);
  Double_t eta = clsPos.Eta();

  if(TMath::Abs(eta) > fEtaCut)
    return kFALSE;  

  if (!IsWithinFiducialVolume(id))
    return kFALSE;

  if(c->GetM02() >fM02Cut)
    return kFALSE;


  return kTRUE;

}
//________________________________________________________________________________________________
Bool_t AliAnalysisTaskPi0V2::IsGoodClusterV1(const AliVCluster *c) const
{

  if(!c)
    return kFALSE;

  if(c->GetNCells() < fNcellCut)
   return kFALSE;

  if(c->E() < fECut)
   return kFALSE;

  Short_t id = -1;
  Double_t maxE = GetMaxCellEnergy(c, id);
     if((1. - double(GetCrossEnergy(c,id))/maxE) > 0.97)
    return kFALSE;


  Float_t pos1[3];
  c->GetPosition(pos1);
  TVector3 clsPos(pos1);
  Double_t eta = clsPos.Eta();

  if(TMath::Abs(eta) > fEtaCut)
    return kFALSE;

  if (!IsWithinFiducialVolume(id))
    return kFALSE;

  if(c->GetM02() <0.5)
    return kFALSE;

  Double_t dr = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx() + c->GetTrackDz()*c->GetTrackDz());
  if(dr>0.025)
    return kFALSE;

  return kTRUE;

}
//_____________________________________________________________________
Bool_t AliAnalysisTaskPi0V2::IsGoodPion(const TLorentzVector &p1, const TLorentzVector &p2) const
{
  // Is good pion?

  if(fPi0AsyCut){
    Double_t asym = TMath::Abs(p1.E()-p2.E())/(p1.E()+p2.E());
    if (asym>0.7)
      return kFALSE;
  }
  TLorentzVector pion;
  pion = p1 + p2;
  Double_t eta = pion.Eta();
  if(TMath::Abs(eta) > fEtaCut)
    return kFALSE;

  return kTRUE;
}
//_______________________________________________________________________
void AliAnalysisTaskPi0V2::FillPion(const TLorentzVector& p1, const TLorentzVector& p2, Double_t EPV0r, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC)
{
  // Fill histogram.

  if (!IsGoodPion(p1,p2))
    return;
  TLorentzVector pion;
  pion = p1 + p2;

  Double_t mass = pion.M();
  Double_t pt   = pion.Pt();
  Double_t phi  = pion.Phi();

  Double_t dphiV0   = phi-EPV0r;
  Double_t dphiV0A  = phi-EPV0A;
  Double_t dphiV0C  = phi-EPV0C;
  Double_t dphiTPC  = phi-EPTPC;

  Double_t cos2phiV0  = TMath::Cos(2.*(dphiV0));
  Double_t cos2phiV0A = TMath::Cos(2.*(dphiV0A));
  Double_t cos2phiV0C = TMath::Cos(2.*(dphiV0C));
  Double_t cos2phiTPC = TMath::Cos(2.*(dphiTPC));

  dphiV0  = TVector2::Phi_0_2pi(dphiV0);  if(dphiV0  >TMath::Pi())  dphiV0 -= TMath::Pi();
  dphiV0A = TVector2::Phi_0_2pi(dphiV0A); if(dphiV0A >TMath::Pi())  dphiV0A -= TMath::Pi();
  dphiV0C = TVector2::Phi_0_2pi(dphiV0C); if(dphiV0C >TMath::Pi())  dphiV0C -= TMath::Pi();
  dphiTPC = TVector2::Phi_0_2pi(dphiTPC); if(dphiTPC >TMath::Pi())  dphiTPC -= TMath::Pi();

  Double_t xV0[5]; // Match ndims in fH  V0 EP
  xV0[0]       = mass;
  xV0[1]       = pt;
  xV0[2]       = fCentrality;
  xV0[3]       = dphiV0;
  xV0[4]       = cos2phiV0;
  fHEPV0r->Fill(xV0);

  Double_t xV0A[5]; // Match ndims in fH V0A EP
  xV0A[0]       = mass;
  xV0A[1]       = pt;
  xV0A[2]       = fCentrality;
  xV0A[3]       = dphiV0A;
  xV0A[4]       = cos2phiV0A;
  fHEPV0A->Fill(xV0A);

  Double_t xV0C[5]; // Match ndims in fH V0C EP
  xV0C[0]       = mass;
  xV0C[1]       = pt;
  xV0C[2]       = fCentrality;
  xV0C[3]       = dphiV0C;
  xV0C[4]       = cos2phiV0C;
  fHEPV0C->Fill(xV0C);

  Double_t xTPC[5]; // Match ndims in fH TPC EP
  xTPC[0]       = mass;
  xTPC[1]       = pt;
  xTPC[2]       = fCentrality;
  xTPC[3]       = dphiTPC;
  xTPC[4]       = cos2phiTPC;
  fHEPTPC->Fill(xTPC);


}
//________________________________________________________________________________________________________________________________
void AliAnalysisTaskPi0V2::FillCluster(const TLorentzVector& p1, Double_t EPV0r, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC)
{
  //cluster(photon) v2 method
  Double_t Pt   = p1.Pt();
  Double_t Phi  = p1.Phi();

  Double_t difClusV0 = TVector2::Phi_0_2pi(Phi-EPV0r);   if(difClusV0 >TMath::Pi()) difClusV0  -= TMath::Pi();
  Double_t difClusV0A = TVector2::Phi_0_2pi(Phi-EPV0A);  if(difClusV0A >TMath::Pi()) difClusV0A -= TMath::Pi();
  Double_t difClusV0C = TVector2::Phi_0_2pi(Phi-EPV0C);  if(difClusV0C >TMath::Pi()) difClusV0C -= TMath::Pi();
  Double_t difClusTPC = TVector2::Phi_0_2pi(Phi-EPTPC);  if(difClusTPC >TMath::Pi()) difClusTPC -= TMath::Pi();

  hdifClus_EPV0->Fill(fCentrality,  difClusV0, Pt);
  hdifClus_EPV0A->Fill(fCentrality, difClusV0A, Pt);
  hdifClus_EPV0C->Fill(fCentrality, difClusV0C, Pt);
  hdifClus_EPTPC->Fill(fCentrality, difClusTPC, Pt);


}
//_________________________________________________________________________________________________
void AliAnalysisTaskPi0V2::GetMom(TLorentzVector& p, const AliVCluster *c, Double_t *vertex)
{
  // Calculate momentum.
  Float_t posMom[3];
  c->GetPosition(posMom);
  TVector3 clsPos2(posMom);

  Double_t e   = c->E();
  Double_t r   = clsPos2.Perp();
  Double_t eta = clsPos2.Eta();
  Double_t phi = clsPos2.Phi();

  TVector3 pos;
  pos.SetPtEtaPhi(r,eta,phi);

  if (vertex) { //calculate direction relative to vertex
    pos -= vertex;
  }

  Double_t rad = pos.Mag();
  p.SetPxPyPzE(e*pos.x()/rad, e*pos.y()/rad, e*pos.z()/rad, e);

}
//________________________________________________________________________
void AliAnalysisTaskPi0V2::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
        
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!

    hEvtCount = new TH1F("hEvtCount", " Event Plane", 10, 0.5, 10.5);
    hEvtCount->GetXaxis()->SetBinLabel(1,"All");
    hEvtCount->GetXaxis()->SetBinLabel(2,"Evt Cut");
    hEvtCount->GetXaxis()->SetBinLabel(3,"Trg Class");
    hEvtCount->GetXaxis()->SetBinLabel(4,"Vtx");
    hEvtCount->GetXaxis()->SetBinLabel(5,"Cent");
    hEvtCount->GetXaxis()->SetBinLabel(5,"EPtask");
    hEvtCount->GetXaxis()->SetBinLabel(7,"EPvalue");
    hEvtCount->GetXaxis()->SetBinLabel(8,"Pass");
    fOutput->Add(hEvtCount);
    
    hEPTPC   = new TH2F("hEPTPC",   "EPTPC     vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hresoTPC = new TH2F("hresoTPC", "TPc reso  vs cent", 100, 0., 100., 100, 0., 1.);
    hEPV0    = new TH2F("hEPV0",    "EPV0      vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0A   = new TH2F("hEPV0A",   "EPV0A     vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0C   = new TH2F("hEPV0C",   "EPV0C     vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0Ar  = new TH2F("hEPV0Ar",  "EPV0Ar    vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0Cr  = new TH2F("hEPV0Cr",  "EPV0Cr    vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0r   = new TH2F("hEPV0r",   "EPV0r     vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0AR4 = new TH2F("hEPV0AR4", "EPV0AR4   vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0AR7 = new TH2F("hEPV0AR7", "EPV0AR7   vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0CR0 = new TH2F("hEPV0CR0", "EPV0CR0   vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    hEPV0CR3 = new TH2F("hEPV0CR3", "EPV0CR3   vs cent", 100, 0., 100., 100, 0., TMath::Pi());
    fOutput->Add(hEPTPC);
    fOutput->Add(hresoTPC);
    fOutput->Add(hEPV0);
    fOutput->Add(hEPV0A);
    fOutput->Add(hEPV0C);
    fOutput->Add(hEPV0Ar);
    fOutput->Add(hEPV0Cr);
    fOutput->Add(hEPV0r);
    fOutput->Add(hEPV0AR4);
    fOutput->Add(hEPV0AR7);
    fOutput->Add(hEPV0CR0);
    fOutput->Add(hEPV0CR3);

    hdifV0Ar_V0Cr    = new TH2F("hdifV0Ar_V0Cr",    "EP Ar-Cr ", 100, 0., 100., 100, -1., 1.);    
    hdifV0A_V0CR0    = new TH2F("hdifV0A_V0CR0",    "EP A-R0 ",  100, 0., 100., 100, -1., 1.);    
    hdifV0A_V0CR3    = new TH2F("hdifV0A_V0CR3",    "EP A-R3 ",  100, 0., 100., 100, -1., 1.);    
    hdifV0ACR0_V0CR3 = new TH2F("hdifV0ACR0_V0CR3", "EP R0-R3 ", 100, 0., 100., 100, -1., 1.);    
    hdifV0C_V0AR4    = new TH2F("hdifV0C_V0AR4",    "EP C-R4 ",  100, 0., 100., 100, -1., 1.);    
    hdifV0C_V0AR7    = new TH2F("hdifV0C_V0AR7",    "EP C-R7 ",  100, 0., 100., 100, -1., 1.);    
    hdifV0AR4_V0AR7  = new TH2F("hdifV0AR4_V0AR7",  "EP R4-R7 ", 100, 0., 100., 100, -1., 1.);    
    fOutput->Add(hdifV0Ar_V0Cr);
    fOutput->Add(hdifV0A_V0CR0);
    fOutput->Add(hdifV0A_V0CR3);
    fOutput->Add(hdifV0ACR0_V0CR3);
    fOutput->Add(hdifV0C_V0AR4);
    fOutput->Add(hdifV0C_V0AR7);
    fOutput->Add(hdifV0AR4_V0AR7);

    hdifV0A_V0C = new TH2F("hdifV0A_V0C", "EP A-C  ", 100, 0., 100., 100, -1., 1.);
    hdifV0A_TPC = new TH2F("hdifV0A_TPC", "EP A-TPC", 100, 0., 100., 100, -1., 1.);
    hdifV0C_TPC = new TH2F("hdifV0C_TPC", "EP C-TPC", 100, 0., 100., 100, -1., 1.);
    hdifV0C_V0A = new TH2F("hdifV0C_V0A", "EP C-A  ", 100, 0., 100., 100, -1., 1.);
    fOutput->Add(hdifV0A_V0C);
    fOutput->Add(hdifV0A_TPC);
    fOutput->Add(hdifV0C_TPC);
    fOutput->Add(hdifV0C_V0A);



    hdifEMC_EPV0  = new TH3F("hdifEMC_EPV0",  "dif phi in EMC with EP",  100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifEMC_EPV0A = new TH3F("hdifEMC_EPV0A", "dif phi in EMC with EP",  100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifEMC_EPV0C = new TH3F("hdifEMC_EPV0C", "dif phi in EMC with EP",  100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hdifEMC_EPV0);
    fOutput->Add(hdifEMC_EPV0A);
    fOutput->Add(hdifEMC_EPV0C);

    hdifful_EPV0 = new TH3F("hdifful_EPV0",    "dif phi in full with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifful_EPV0A = new TH3F("hdifful_EPV0A",  "dif phi in full with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifful_EPV0C = new TH3F("hdifful_EPV0C",  "dif phi in full with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hdifful_EPV0);
    fOutput->Add(hdifful_EPV0A);
    fOutput->Add(hdifful_EPV0C);

    hdifout_EPV0  = new TH3F("hdifout_EPV0",  "dif phi NOT in EMC with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifout_EPV0A = new TH3F("hdifout_EPV0A", "dif phi NOT in EMC with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifout_EPV0C = new TH3F("hdifout_EPV0C", "dif phi NOT in EMC with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hdifout_EPV0);
    fOutput->Add(hdifout_EPV0A);
    fOutput->Add(hdifout_EPV0C);

    hdifEMC_EPTPC = new TH3F("hdifEMC_EPTPC", "dif phi in EMC with EP",  100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifful_EPTPC = new TH3F("hdifful_EPTPC",  "dif phi in full with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifout_EPTPC = new TH3F("hdifout_EPTPC", "dif phi NOT in EMC with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hdifEMC_EPTPC);
    fOutput->Add(hdifful_EPTPC);
    fOutput->Add(hdifout_EPTPC);

    hdifClus_EPV0 = new TH3F("hdifClus_EPV0", "dif phi in EMC Clus with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifClus_EPV0A = new TH3F("hdifClus_EPV0A", "dif phi in EMC Clus with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifClus_EPV0C = new TH3F("hdifClus_EPV0C", "dif phi in EMC Clus with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    hdifClus_EPTPC = new TH3F("hdifClus_EPTPC", "dif phi in EMC Clus with EP", 100, 0., 100., 100, 0., TMath::Pi(), 15, 0., 15.);
    fOutput->Add(hdifClus_EPV0);
    fOutput->Add(hdifClus_EPV0A);
    fOutput->Add(hdifClus_EPV0C);
    fOutput->Add(hdifClus_EPTPC);

    hAllcentV0  = new TH1F("hAllcentV0",  "All cent EP V0",  100, 0., TMath::Pi());
    hAllcentV0r = new TH1F("hAllcentV0r", "All cent EP V0r", 100, 0., TMath::Pi());
    hAllcentV0A = new TH1F("hAllcentV0A", "All cent EP V0A", 100, 0., TMath::Pi());
    hAllcentV0C = new TH1F("hAllcentV0C", "All cent EP V0C", 100, 0., TMath::Pi());
    hAllcentTPC = new TH1F("hAllcentTPC", "All cent EP TPC", 100, 0., TMath::Pi());
    fOutput->Add(hAllcentV0);
    fOutput->Add(hAllcentV0r);
    fOutput->Add(hAllcentV0A);
    fOutput->Add(hAllcentV0C);
    fOutput->Add(hAllcentTPC);

    h2DcosV0r = new TH2F("h2DcosV0r", "cos(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DsinV0r = new TH2F("h2DsinV0r", "sin(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DcosV0A = new TH2F("h2DcosV0A", "cos(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DsinV0A = new TH2F("h2DsinV0A", "sin(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DcosV0C = new TH2F("h2DcosV0C", "cos(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DsinV0C = new TH2F("h2DsinV0C", "sin(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DcosTPC = new TH2F("h2DcosTPC", "cos(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    h2DsinTPC = new TH2F("h2DsinTPC", "sin(Phi) V0r vs Run NUmber", 200, 0, 200, 100, -1, 1);
    fOutput->Add(h2DcosV0r);
    fOutput->Add(h2DsinV0r);
    fOutput->Add(h2DcosV0A);
    fOutput->Add(h2DsinV0A);
    fOutput->Add(h2DcosV0C);
    fOutput->Add(h2DsinV0C);
    fOutput->Add(h2DcosTPC);
    fOutput->Add(h2DsinTPC);

    hM02vsPtA = new TH2F("hM02vsPtA", "M02 vs Et before cut", 5000, 0, 50, 400, 0, 4.);
    hM02vsPtB = new TH2F("hM02vsPtB", "M02 vs Et before cut", 5000, 0, 50, 400, 0, 4.);
    fOutput->Add(hM02vsPtA);
    fOutput->Add(hM02vsPtB);

    hClusDxDZA = new TH2F("hClusDxDZA", "clus Dx vs Dz", 1000, -1., 1., 1000, -1., 1);
    hClusDxDZB = new TH2F("hClusDxDZB", "clus Dx vs Dz", 1000, -1., 1., 1000, -1., 1);
    fOutput->Add(hClusDxDZA);
    fOutput->Add(hClusDxDZB);

    const Int_t ndims = 5;
    Int_t nMgg=500, nPt=40, nCent=20, nDeltaPhi=315,  ncos2phi=500;
    Int_t bins[ndims] = {nMgg, nPt, nCent, nDeltaPhi, ncos2phi};
    Double_t xmin[ndims] = { 0,   0.,  0,   0.,     -1.};
    Double_t xmax[ndims] = { 0.5, 20., 100, 3.15,   1.};
    fHEPV0r  = new THnSparseF("fHEPV0r",  "Flow histogram EPV0",  ndims, bins, xmin, xmax);
    fHEPV0A = new THnSparseF("fHEPV0A",   "Flow histogram EPV0A", ndims, bins, xmin, xmax);
    fHEPV0C = new THnSparseF("fHEPV0C",   "Flow histogram EPV0C", ndims, bins, xmin, xmax);
    fHEPTPC = new THnSparseF("fHEPTPC",   "Flow histogram EPTPC", ndims, bins, xmin, xmax);
    fOutput->Add(fHEPV0r);
    fOutput->Add(fHEPV0A);
    fOutput->Add(fHEPV0C);
    fOutput->Add(fHEPTPC);

    

    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPi0V2::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
        
    // Create pointer to reconstructed event
   AliVEvent *event = InputEvent();
   if (!event) { Printf("ERROR: Could not retrieve event"); return; }
   // create pointer to event
   fESD = dynamic_cast<AliESDEvent*>(event);
   if (!fESD) {
     AliError("Cannot get the ESD event");
     return;
   }
   hEvtCount->Fill(1);
    
  Int_t AbsRunNumber = fESD->GetRunNumber();
  fRunNumber = ConvertToInternalRunNumber(AbsRunNumber);

  Bool_t isSelected =0;      
  if(fEvtSelect == 1){  //MB+SemiCentral
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kAnyINT | AliVEvent::kSemiCentral));
  } else if (fEvtSelect == 2){  //MB+Central+SemiCentral
    isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kAnyINT | AliVEvent::kSemiCentral | AliVEvent::kCentral));
  } else if(fEvtSelect == 3){  //MB
 isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kAnyINT ));
  }

  if(!isSelected )
        return; 

  hEvtCount->Fill(2);
  if(!fTrigClass.IsNull()){
    TString fired;
    fired = fESD->GetFiredTriggerClasses();
    if (!fired.Contains("-B-"))
      return;
    TObjArray *arr = fTrigClass.Tokenize("|");
    if (!arr)
      return;
    Bool_t match = 0;
    for (Int_t i=0;i<arr->GetEntriesFast();++i) {
      TObject *obj = arr->At(i);
      if (!obj)
	continue;
      if (fired.Contains(obj->GetName())) {
	match = 1;
	break;
      }
    }
    delete arr;
    if (
	!match && //select by Trigger classes in KCentral and KSemiCentral
        !(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & (AliVEvent::kAnyINT ))              // always accept MB
	) 
      return; //Not match skip this event
  }

    hEvtCount->Fill(3);
    const AliESDVertex* fvertex = fESD->GetPrimaryVertex();
    if(TMath::Abs(fvertex->GetZ())>fVtxCut)
      return;
    Double_t vertex[3] = {fvertex->GetX(), fvertex->GetY(), fvertex->GetZ()};

    hEvtCount->Fill(4);

    if(fESD->GetCentrality()) {
      fCentrality = 
	fESD->GetCentrality()->GetCentralityPercentile("CL1"); //spd vertex
    } else{
	   return;
    }

    hEvtCount->Fill(5);
    AliEventplane *ep = fESD->GetEventplane();
    if (ep) {
      if (ep->GetQVector())
	fEPTPC    = ep->GetQVector()->Phi()/2. ;
      else
	fEPTPC = -999.;
      if (ep->GetQsub1()&&ep->GetQsub2())
	fEPTPCreso  = TMath::Cos(2.*(ep->GetQsub1()->Phi()/2.-ep->GetQsub2()->Phi()/2.));
      else
	fEPTPCreso = -1;

      fEPV0    = ep->GetEventplane("V0",  fESD);
      fEPV0A   = ep->GetEventplane("V0A", fESD);
      fEPV0C   = ep->GetEventplane("V0C", fESD);
      Double_t qx=0, qy=0;
      Double_t qxr=0, qyr=0;
      fEPV0Ar  = ep->CalculateVZEROEventPlane(fESD, 4, 5, 2, qxr, qyr);
      fEPV0Cr  = ep->CalculateVZEROEventPlane(fESD, 2, 3, 2, qx,  qy);
      qxr += qx;
      qyr += qy;
      fEPV0r   = TMath::ATan2(qyr,qxr)/2.;
      fEPV0AR4 = ep->CalculateVZEROEventPlane(fESD, 4, 2, qx, qy);
      fEPV0AR5 = ep->CalculateVZEROEventPlane(fESD, 5, 2, qx, qy);
      fEPV0AR6 = ep->CalculateVZEROEventPlane(fESD, 6, 2, qx, qy);
      fEPV0AR7 = ep->CalculateVZEROEventPlane(fESD, 7, 2, qx, qy);
      fEPV0CR0 = ep->CalculateVZEROEventPlane(fESD, 0, 2, qx, qy);
      fEPV0CR1 = ep->CalculateVZEROEventPlane(fESD, 1, 2, qx, qy);
      fEPV0CR2 = ep->CalculateVZEROEventPlane(fESD, 2, 2, qx, qy);
      fEPV0CR3 = ep->CalculateVZEROEventPlane(fESD, 3, 2, qx, qy);
    }
    FillEPQA(); //Fill the EP QA

    hEvtCount->Fill(6);

    if( fEPV0A<-2. || fEPV0C<-2. || fEPTPC<-2. || fEPV0r<-2.) 
      return;

    hEvtCount->Fill(7);

    fEPV0   = TVector2::Phi_0_2pi(fEPV0);    if(fEPV0>TMath::Pi())   fEPV0  = fEPV0  - TMath::Pi();
    fEPV0r  = TVector2::Phi_0_2pi(fEPV0r);   if(fEPV0r>TMath::Pi())  fEPV0r = fEPV0r - TMath::Pi();
    fEPV0A  = TVector2::Phi_0_2pi(fEPV0A);   if(fEPV0A>TMath::Pi())  fEPV0A = fEPV0A - TMath::Pi();
    fEPV0C  = TVector2::Phi_0_2pi(fEPV0C);   if(fEPV0C>TMath::Pi())  fEPV0C = fEPV0C - TMath::Pi();
    fEPV0Ar = TVector2::Phi_0_2pi(fEPV0Ar);  if(fEPV0Ar>TMath::Pi()) fEPV0Ar = fEPV0Ar - TMath::Pi();
    fEPV0Cr = TVector2::Phi_0_2pi(fEPV0Cr);  if(fEPV0Cr>TMath::Pi()) fEPV0Cr = fEPV0Cr - TMath::Pi();
    fEPV0AR4   = TVector2::Phi_0_2pi(fEPV0AR4);    if(fEPV0AR4>TMath::Pi())   fEPV0AR4  = fEPV0AR4 - TMath::Pi();
    fEPV0AR7   = TVector2::Phi_0_2pi(fEPV0AR7);    if(fEPV0AR7>TMath::Pi())   fEPV0AR7  = fEPV0AR7 - TMath::Pi();
    fEPV0CR0   = TVector2::Phi_0_2pi(fEPV0CR0);    if(fEPV0CR0>TMath::Pi())   fEPV0CR0  = fEPV0CR0 - TMath::Pi();
    fEPV0CR3   = TVector2::Phi_0_2pi(fEPV0CR3);    if(fEPV0CR3>TMath::Pi())   fEPV0CR3  = fEPV0CR3 - TMath::Pi();

   if(fEPTPC != -999.)
   hEPTPC->Fill(fCentrality,  fEPTPC); 
   if(fEPTPCreso!=-1) hresoTPC->Fill(fCentrality, fEPTPCreso);
   hEPV0->Fill(fCentrality,   fEPV0);
   hEPV0A->Fill(fCentrality,  fEPV0A);
   hEPV0C->Fill(fCentrality,  fEPV0C);
   hEPV0Ar->Fill(fCentrality, fEPV0Ar);
   hEPV0Cr->Fill(fCentrality, fEPV0Cr);
   hEPV0r->Fill(fCentrality,  fEPV0r);
   hEPV0AR4->Fill(fCentrality, fEPV0AR4);
   hEPV0AR7->Fill(fCentrality, fEPV0AR7);
   hEPV0CR0->Fill(fCentrality, fEPV0CR0);
   hEPV0CR3->Fill(fCentrality, fEPV0CR3);

   hAllcentV0->Fill(fEPV0);
   hAllcentV0r->Fill(fEPV0r);
   hAllcentV0A->Fill(fEPV0A);
   hAllcentV0C->Fill(fEPV0C);  
   hAllcentTPC->Fill(fEPTPC);

   hdifV0Ar_V0Cr->Fill(fCentrality, TMath::Cos(2.*(fEPV0Ar - fEPV0Cr)));
   hdifV0A_V0CR0->Fill(fCentrality, TMath::Cos(2.*(fEPV0A - fEPV0CR0)));
   hdifV0A_V0CR3->Fill(fCentrality, TMath::Cos(2.*(fEPV0A - fEPV0CR3)));
   hdifV0ACR0_V0CR3->Fill(fCentrality, TMath::Cos(2*(fEPV0CR0 - fEPV0CR3)));
   hdifV0C_V0AR4->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0AR4)));
   hdifV0C_V0AR7->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0AR7)));
   hdifV0AR4_V0AR7->Fill(fCentrality, TMath::Cos(2*(fEPV0AR4 - fEPV0AR7)));
        
   hdifV0A_V0C->Fill(fCentrality, TMath::Cos(2*(fEPV0A - fEPV0C)));
   hdifV0A_TPC->Fill(fCentrality, TMath::Cos(2*(fEPV0A - fEPTPC)));
   hdifV0C_TPC->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPTPC)));
   hdifV0C_V0A->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0A)));
    // Cluster loop for reconstructed event

   if (!fV2ClusName.IsNull() && !fV2Clus) {
     fV2Clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV2ClusName));
     if (!fV2Clus) {
       AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fV2ClusName.Data()));
       return;
     }
  }


//   Int_t nCluster =  fESD->GetNumberOfCaloClusters(); 
   Int_t nCluster =  fV2Clus->GetEntries(); 
   for(Int_t i=0; i<nCluster; ++i){
//     AliESDCaloCluster *c1 = fESD->GetCaloCluster(i);
     AliVCluster *c1 = static_cast<AliVCluster*>(fV2Clus->At(i));      
     hClusDxDZA->Fill(c1->GetTrackDz(), c1->GetTrackDx());
     Float_t clsPosEt[3] = {0,0,0};
     c1->GetPosition(clsPosEt);
     TVector3 clsVec(clsPosEt);
     Double_t Et = c1->E()*TMath::Sin(clsVec.Theta());
     hM02vsPtA->Fill(Et, c1->GetM02());
     if(!c1->IsEMCAL()) continue;
     if(!IsGoodCluster(c1)) continue;
     hM02vsPtB->Fill(Et, c1->GetM02());	
     hClusDxDZB->Fill(c1->GetTrackDz(), c1->GetTrackDx());
     TLorentzVector p1;
     GetMom(p1, c1, vertex);
     for(Int_t j=i+1; j<nCluster; ++j){
//       AliESDCaloCluster *c2 = fESD->GetCaloCluster(j);
       AliVCluster *c2 = static_cast<AliVCluster*>(fV2Clus->At(i));      
       if(!c2->IsEMCAL()) continue;
       if(!IsGoodCluster(c2)) continue;
       TLorentzVector p2;
       GetMom(p2, c2, vertex);
       FillPion(p1, p2, fEPV0r, fEPV0A, fEPV0C, fEPTPC);
     }
   }

  if(isV1Clus){
    if (!fV2ClusName.IsNull() && !fV1Clus) {
      fV1Clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV1ClusName));
      if (!fV1Clus) {
	AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fV1ClusName.Data()));
	return;
      }
    }
    Int_t nClusterV1 = fV1Clus->GetEntries();
    for(Int_t i=0; i<nClusterV1; ++i){
//      AliESDCaloCluster *c3 = fESD->GetCaloCluster(i);      
      AliVCluster *c3 = static_cast<AliVCluster*>(fV1Clus->At(i));      
      if(!c3->IsEMCAL()) continue;
      if(!IsGoodClusterV1(c3)) continue;
      TLorentzVector p3;
      GetMom(p3, c3, vertex);
      FillCluster(p3, fEPV0r, fEPV0A, fEPV0C, fEPTPC);
    }
  }


   if (!fTracksName.IsNull() && !fTracks) {
     fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
     if (!fTracks) {
       AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fTracksName.Data())); 
       return;
     }
  }

   Int_t ntracks = fTracks->GetEntries();
   for(Int_t i=0; i<ntracks; ++i){
     AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(i));
     Double_t tPhi = track->Phi();
     Double_t tPt  = track->Pt();

     Double_t difTrackV0  = TVector2::Phi_0_2pi(tPhi-fEPV0);   if(difTrackV0  >TMath::Pi()) difTrackV0  -= TMath::Pi();
     Double_t difTrackV0A = TVector2::Phi_0_2pi(tPhi-fEPV0A);  if(difTrackV0A >TMath::Pi()) difTrackV0A -= TMath::Pi();
     Double_t difTrackV0C = TVector2::Phi_0_2pi(tPhi-fEPV0C);  if(difTrackV0C >TMath::Pi()) difTrackV0C -= TMath::Pi();
     Double_t difTrackTPC = TVector2::Phi_0_2pi(tPhi-fEPTPC);  if(difTrackTPC >TMath::Pi()) difTrackTPC -= TMath::Pi();
     if(track->IsEMCAL()){	
       hdifEMC_EPV0->Fill(fCentrality, difTrackV0, tPt);
       hdifEMC_EPV0A->Fill(fCentrality, difTrackV0A, tPt);
       hdifEMC_EPV0C->Fill(fCentrality, difTrackV0C, tPt);
       hdifEMC_EPTPC->Fill(fCentrality, difTrackTPC, tPt);
     }else{
       hdifout_EPV0->Fill(fCentrality, difTrackV0, tPt);
       hdifout_EPV0A->Fill(fCentrality, difTrackV0A, tPt);
       hdifout_EPV0C->Fill(fCentrality, difTrackV0C, tPt);
       hdifout_EPTPC->Fill(fCentrality, difTrackTPC, tPt);
     }
     hdifful_EPV0->Fill(fCentrality,    difTrackV0, tPt);
     hdifful_EPV0A->Fill(fCentrality,   difTrackV0A, tPt);
     hdifful_EPV0C->Fill(fCentrality,   difTrackV0C, tPt);
     hdifful_EPTPC->Fill(fCentrality,   difTrackTPC, tPt);
    }
    hEvtCount->Fill(8);

    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
}
//____________________________________________________________________
Int_t AliAnalysisTaskPi0V2::ConvertToInternalRunNumber(Int_t n)
{
    switch(n){
    case  170593 : return 179 ;
    case  170572 : return 178 ;
    case  170556 : return 177 ;
    case  170552 : return 176 ;
    case  170546 : return 175 ;
    case  170390 : return 174 ;
    case  170389 : return 173 ;
    case  170388 : return 172 ;
    case  170387 : return 171 ;
    case  170315 : return 170 ;
    case  170313 : return 169 ;
    case  170312 : return 168 ;
    case  170311 : return 167 ;
    case  170309 : return 166 ;
    case  170308 : return 165 ;
    case  170306 : return 164 ;
    case  170270 : return 163 ;
    case  170269 : return 162 ;
    case  170268 : return 161 ;
    case  170267 : return 160 ;
    case  170264 : return 159 ;
    case  170230 : return 158 ;
    case  170228 : return 157 ;
    case  170208 : return 156 ;
    case  170207 : return 155 ;
    case  170205 : return 154 ;
    case  170204 : return 153 ;
    case  170203 : return 152 ;
    case  170195 : return 151 ;
    case  170193 : return 150 ;
    case  170163 : return 149 ;
    case  170162 : return 148 ;
    case  170159 : return 147 ;
    case  170155 : return 146 ;
    case  170152 : return 145 ;
    case  170091 : return 144 ;
    case  170089 : return 143 ;
    case  170088 : return 142 ;
    case  170085 : return 141 ;
    case  170084 : return 140 ;
    case  170083 : return 139 ;
    case  170081 : return 138 ;
    case  170040 : return 137 ;
    case  170038 : return 136 ;
    case  170036 : return 135 ;
    case  170027 : return 134 ;
    case  169981 : return 133 ;
    case  169975 : return 132 ;
    case  169969 : return 131 ;
    case  169965 : return 130 ;
    case  169961 : return 129 ;
    case  169956 : return 128 ;
    case  169926 : return 127 ;
    case  169924 : return 126 ;
    case  169923 : return 125 ;
    case  169922 : return 124 ;
    case  169919 : return 123 ;
    case  169918 : return 122 ;
    case  169914 : return 121 ;
    case  169859 : return 120 ;
    case  169858 : return 119 ;
    case  169855 : return 118 ;
    case  169846 : return 117 ;
    case  169838 : return 116 ;
    case  169837 : return 115 ;
    case  169835 : return 114 ;
    case  169683 : return 113 ;
    case  169628 : return 112 ;
    case  169591 : return 111 ;
    case  169590 : return 110 ;
    case  169588 : return 109 ;
    case  169587 : return 108 ;
    case  169586 : return 107 ;
    case  169584 : return 106 ;
    case  169557 : return 105 ;
    case  169555 : return 104 ;
    case  169554 : return 103 ;
    case  169553 : return 102 ;
    case  169550 : return 101 ;
    case  169515 : return 100 ;
    case  169512 : return 99 ;
    case  169506 : return 98 ;
    case  169504 : return 97 ;
    case  169498 : return 96 ;
    case  169475 : return 95 ;
    case  169420 : return 94 ;
    case  169419 : return 93 ;
    case  169418 : return 92 ;
    case  169417 : return 91 ;
    case  169415 : return 90 ;
    case  169411 : return 89 ;
    case  169238 : return 88 ;
    case  169236 : return 87 ;
    case  169167 : return 86 ;
    case  169160 : return 85 ;
    case  169156 : return 84 ;
    case  169148 : return 83 ;
    case  169145 : return 82 ;
    case  169144 : return 81 ;
    case  169143 : return 80 ;
    case  169138 : return 79 ;
    case  169099 : return 78 ;
    case  169094 : return 77 ;
    case  169091 : return 76 ;
    case  169045 : return 75 ;
    case  169044 : return 74 ;
    case  169040 : return 73 ;
    case  169035 : return 72 ;
    case  168992 : return 71 ;
    case  168988 : return 70 ;
    case  168984 : return 69 ;
    case  168826 : return 68 ;
    case  168777 : return 67 ;
    case  168514 : return 66 ;
    case  168512 : return 65 ;
    case  168511 : return 64 ;
    case  168467 : return 63 ;
    case  168464 : return 62 ;
    case  168461 : return 61 ;
    case  168460 : return 60 ;
    case  168458 : return 59 ;
    case  168362 : return 58 ;
    case  168361 : return 57 ;
    case  168356 : return 56 ;
    case  168342 : return 55 ;
    case  168341 : return 54 ;
    case  168325 : return 53 ;
    case  168322 : return 52 ;
    case  168318 : return 51 ;
    case  168311 : return 50 ;
    case  168310 : return 49 ;
    case  168213 : return 48 ;
    case  168212 : return 47 ;
    case  168208 : return 46 ;
    case  168207 : return 45 ;
    case  168206 : return 44 ;
    case  168205 : return 43 ;
    case  168204 : return 42 ;
    case  168203 : return 41 ;
    case  168181 : return 40 ;
    case  168177 : return 39 ;
    case  168175 : return 38 ;
    case  168173 : return 37 ;
    case  168172 : return 36 ;
    case  168171 : return 35 ;
    case  168115 : return 34 ;
    case  168108 : return 33 ;
    case  168107 : return 32 ;
    case  168105 : return 31 ;
    case  168104 : return 30 ;
    case  168103 : return 29 ;
    case  168076 : return 28 ;
    case  168069 : return 27 ;
    case  168068 : return 26 ;
    case  168066 : return 25 ;
    case  167988 : return 24 ;
    case  167987 : return 23 ;
    case  167986 : return 22 ;
    case  167985 : return 21 ;
    case  167921 : return 20 ;
    case  167920 : return 19 ;
    case  167915 : return 18 ;
    case  167909 : return 17 ;
    case  167903 : return 16 ;
    case  167902 : return 15 ;
    case  167818 : return 14 ;
    case  167814 : return 13 ;
    case  167813 : return 12 ;
    case  167808 : return 11 ;
    case  167807 : return 10 ;
    case  167806 : return 9 ;
    case  167713 : return 8 ;
    case  167712 : return 7 ;
    case  167711 : return 6 ;
    case  167706 : return 5 ;
    case  167693 : return 4 ;
    case  166532 : return 3 ;
    case  166530 : return 2 ;
    case  166529 : return 1 ;

    default : return 199;
    }
}
//_______________________________________________________________________
void AliAnalysisTaskPi0V2::FillEPQA()
{
  
  h2DcosV0r->Fill(fRunNumber, TMath::Cos(fEPV0r));
  h2DsinV0r->Fill(fRunNumber, TMath::Sin(fEPV0r));
  h2DcosV0A->Fill(fRunNumber, TMath::Cos(fEPV0A));
  h2DsinV0A->Fill(fRunNumber, TMath::Sin(fEPV0A));
  h2DcosV0C->Fill(fRunNumber, TMath::Cos(fEPV0C));
  h2DsinV0C->Fill(fRunNumber, TMath::Sin(fEPV0C));
  h2DcosTPC->Fill(fRunNumber, TMath::Cos(fEPTPC));
  h2DsinTPC->Fill(fRunNumber, TMath::Sin(fEPTPC));


}
//________________________________________________________________________
void AliAnalysisTaskPi0V2::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
//    fOutput = dynamic_cast<TList*> (GetOutputData(1));
 //   if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
        
    // Get the physics selection histograms with the selection statistics
    //AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    //AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    //TH2F *histStat = (TH2F*)inputH->GetStatistics();
   
   
    // NEW HISTO should be retrieved from the TList container in the above way,
    // so it is available to draw on a canvas such as below

}
