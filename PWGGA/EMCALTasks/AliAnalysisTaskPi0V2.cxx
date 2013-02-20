/* $Id: AliAnalysisTaskPi0V2.cxx 55404 2012-03-29 10:10:19Z fca $ */

#include "AliAnalysisTaskPi0V2.h"

#include <Riostream.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include "AliAODEvent.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPID.h"
#include "AliCaloTrackReader.h"
#include "AliCalorimeterUtils.h"
#include "AliEMCALGeometry.h"
#include "AliEPFlattener.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliMCEvent.h"
#include "AliOADBContainer.h"
#include "AliStack.h"
#include "AliVCluster.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskPi0V2)

//________________________________________________________________________
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2(const char *name) :
  AliAnalysisTaskSE(name),
  fOutput(0),
  fESD(0),fAOD(0),
  fTracksName("PicoTrack"), fV1ClusName("CaloCluster"), fV2ClusName("CaloCluster"),
  fTrigClass("CVLN_|CSEMI_|CCENT|CVHN"),
  fTracks(0), fV1Clus(0), fV2Clus(0),
  fRunNumber(-999),fInterRunNumber(-999),
  fVtxCut(15.),
  fNcellCut(2.), fECut(1.), fEtaCut(0.65), fM02Cut(0.5),fDrCut(0.025), fPi0AsyCut(0), isV1Clus(1), isPhosCali(0),
  fCentrality(99.),
  fEPTPC(-999.),
  fEPTPCreso(0.), 
  fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.), fEPV0Ar(-999.), fEPV0Cr(-999.), fEPV0r(-999.),
  fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.), fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
  hEvtCount(0), 
  h2DcosV0A(0), h2DsinV0A(0), h2DcosV0C(0), h2DsinV0C(0), h2DcosTPC(0), h2DsinTPC(0), 
  hEPTPC(0), hresoTPC(0),
  hEPV0(0), hEPV0A(0), hEPV0C(0), hEPV0Ar(0), hEPV0Cr(0), hEPV0r(0), hEPV0AR4(0), hEPV0AR7(0), hEPV0CR0(0), hEPV0CR3(0),
  hEPTPCCor(0), hEPV0ACor(0), hEPV0CCor(0),
  hdifV0Ar_V0Cr(0), hdifV0A_V0CR0(0), hdifV0A_V0CR3(0), hdifV0ACR0_V0CR3(0), hdifV0C_V0AR4(0), hdifV0C_V0AR7(0), hdifV0AR4_V0AR7(0),
  hdifV0A_V0C(0), hdifV0A_TPC(0), hdifV0C_TPC(0), hdifV0C_V0A(0), 
  hM02vsPtA(0), hM02vsPtB(0), hClusDxDZA(0), hClusDxDZB(0),
  hdifEMC_EPV0(0), hdifEMC_EPV0A(0), hdifEMC_EPV0C(0), hdifful_EPV0(0), hdifful_EPV0A(0), hdifful_EPV0C(0), 
  hdifout_EPV0(0), hdifout_EPV0A(0), hdifout_EPV0C(0), 
  fEPcalibFileName("$ALICE_ROOT/OADB/PHOS/PHOSflat.root"), fTPCFlat(0x0), fV0AFlat(0x0),  fV0CFlat(0x0),
  fClusterPbV0(0), fClusterPbV0A(0), fClusterPbV0C(0), fClusterPbTPC(0),    
  fHEPV0r(0x0), fHEPV0A(0x0), fHEPV0C(0x0), fHEPTPC(0x0)
{
  // Dummy constructor ALWAYS needed for I/O.
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskPi0V2::AliAnalysisTaskPi0V2() :
  AliAnalysisTaskSE("default_name"),
  fOutput(0),
  fESD(0),fAOD(0),
  fTracksName("PicoTrack"), fV1ClusName("CaloCluster"), fV2ClusName("CaloCluster"),
  fTrigClass("CVLN_|CSEMI_|CCENT|CVHN"),
  fTracks(0), fV1Clus(0), fV2Clus(0),
  fRunNumber(-999),fInterRunNumber(-999),
  fVtxCut(15.),
  fNcellCut(2.), fECut(1.), fEtaCut(0.65), fM02Cut(0.5), fDrCut(0.025), fPi0AsyCut(0), isV1Clus(1),isPhosCali(0),
  fCentrality(99.),
  fEPTPC(-999.),
  fEPTPCreso(0.),
  fEPV0(-999.), fEPV0A(-999.), fEPV0C(-999.), fEPV0Ar(-999.), fEPV0Cr(-999.), fEPV0r(-999.),
  fEPV0AR4(-999.), fEPV0AR5(-999.), fEPV0AR6(-999.), fEPV0AR7(-999.), fEPV0CR0(-999.), fEPV0CR1(-999.), fEPV0CR2(-999.), fEPV0CR3(-999.),
  hEvtCount(0), 
  h2DcosV0A(0), h2DsinV0A(0), h2DcosV0C(0), h2DsinV0C(0), h2DcosTPC(0), h2DsinTPC(0),
  hEPTPC(0), hresoTPC(0),
  hEPV0(0), hEPV0A(0), hEPV0C(0), hEPV0Ar(0), hEPV0Cr(0), hEPV0r(0), hEPV0AR4(0), hEPV0AR7(0), hEPV0CR0(0), hEPV0CR3(0),
  hEPTPCCor(0), hEPV0ACor(0), hEPV0CCor(0),
  hdifV0Ar_V0Cr(0), hdifV0A_V0CR0(0), hdifV0A_V0CR3(0), hdifV0ACR0_V0CR3(0), hdifV0C_V0AR4(0), hdifV0C_V0AR7(0), hdifV0AR4_V0AR7(0),
  hdifV0A_V0C(0), hdifV0A_TPC(0), hdifV0C_TPC(0), hdifV0C_V0A(0),
  hM02vsPtA(0), hM02vsPtB(0), hClusDxDZA(0), hClusDxDZB(0),
  hdifEMC_EPV0(0), hdifEMC_EPV0A(0), hdifEMC_EPV0C(0), hdifful_EPV0(0), hdifful_EPV0A(0), hdifful_EPV0C(0),
  hdifout_EPV0(0), hdifout_EPV0A(0), hdifout_EPV0C(0), 
  fEPcalibFileName("$ALICE_ROOT/OADB/PHOS/PHOSflat.root"), fTPCFlat(0x0), fV0AFlat(0x0),  fV0CFlat(0x0),
  fClusterPbV0(0), fClusterPbV0A(0), fClusterPbV0C(0), fClusterPbTPC(0),    
  fHEPV0r(0x0), fHEPV0A(0x0), fHEPV0C(0x0), fHEPTPC(0x0)
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
  if (fTPCFlat) 
    delete fTPCFlat;  
  fTPCFlat=0x0;
  if (fV0AFlat) 
    delete fV0AFlat;  
  fV0AFlat=0x0;
  if (fV0CFlat) 
    delete fV0CFlat;  
  fV0CFlat=0x0;
  delete fOutput;
}

//_____________________________________________________________________
Double_t AliAnalysisTaskPi0V2::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = 0;
  if (fESD) {
    cells = fESD->GetEMCALCells();
  } else {
    cells = fAOD->GetEMCALCells();
  }
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

  AliVCaloCells *cells;
  if (fESD) {
    cells = fESD->GetEMCALCells();
  } else {
    cells = fAOD->GetEMCALCells();
  }
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
  if (!c)
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

  if (TMath::Abs(eta) > fEtaCut)
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
  if (!c)
    return kFALSE;

  if (c->GetNCells() < fNcellCut)
   return kFALSE;

  if (c->E() < fECut)
   return kFALSE;

  Short_t id = -1;
  Double_t maxE = GetMaxCellEnergy(c, id);
  if((1. - double(GetCrossEnergy(c,id))/maxE) > 0.97)
    return kFALSE;


  Float_t pos1[3];
  c->GetPosition(pos1);
  TVector3 clsPos(pos1);
  Double_t eta = clsPos.Eta();

  if (TMath::Abs(eta) > fEtaCut)
    return kFALSE;

  if (!IsWithinFiducialVolume(id))
    return kFALSE;

  if (c->GetM02() <fM02Cut)
    return kFALSE;

  Double_t dr = TMath::Sqrt(c->GetTrackDx()*c->GetTrackDx() + c->GetTrackDz()*c->GetTrackDz());
  if(dr<fDrCut)
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
void AliAnalysisTaskPi0V2::FillPion(const TLorentzVector& p1, const TLorentzVector& p2, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC)
{
  // Fill histogram.

  if (!IsGoodPion(p1,p2))
    return;
  TLorentzVector pion;
  pion = p1 + p2;

  Double_t mass = pion.M();
  Double_t pt   = pion.Pt();
  Double_t phi  = pion.Phi();

  Double_t dphiV0A  = phi-EPV0A;
  Double_t dphiV0C  = phi-EPV0C;
  Double_t dphiTPC  = phi-EPTPC;

  dphiV0A = TVector2::Phi_0_2pi(dphiV0A); if(dphiV0A >TMath::Pi())  dphiV0A -= TMath::Pi();
  dphiV0C = TVector2::Phi_0_2pi(dphiV0C); if(dphiV0C >TMath::Pi())  dphiV0C -= TMath::Pi();
  dphiTPC = TVector2::Phi_0_2pi(dphiTPC); if(dphiTPC >TMath::Pi())  dphiTPC -= TMath::Pi();

  Double_t xV0A[4]; // Match ndims in fH V0A EP
  xV0A[0]       = mass;
  xV0A[1]       = pt;
  xV0A[2]       = fCentrality;
  xV0A[3]       = dphiV0A;
  fHEPV0A->Fill(xV0A);

  Double_t xV0C[4]; // Match ndims in fH V0C EP
  xV0C[0]       = mass;
  xV0C[1]       = pt;
  xV0C[2]       = fCentrality;
  xV0C[3]       = dphiV0C;
  fHEPV0C->Fill(xV0C);

  Double_t xTPC[4]; // Match ndims in fH TPC EP
  xTPC[0]       = mass;
  xTPC[1]       = pt;
  xTPC[2]       = fCentrality;
  xTPC[3]       = dphiTPC;
  fHEPTPC->Fill(xTPC);
}

//________________________________________________________________________________________________________________________________
void AliAnalysisTaskPi0V2::FillCluster(const TLorentzVector& p1, Double_t EPV0A, Double_t EPV0C, Double_t EPTPC, AliVCluster *c)
{
  // Cluster(photon) v2 method

  Double_t Et   = p1.Et();
  Double_t Phi  = p1.Phi();
  Double_t M02  = c->GetM02();
  Double_t DxClus = c->GetTrackDx();
  Double_t DzClus = c->GetTrackDz();
  Double_t dr = TMath::Sqrt(DxClus*DxClus + DzClus*DzClus);

  Double_t difClusV0A = TVector2::Phi_0_2pi(Phi-EPV0A);  if(difClusV0A >TMath::Pi()) difClusV0A -= TMath::Pi();
  Double_t difClusV0C = TVector2::Phi_0_2pi(Phi-EPV0C);  if(difClusV0C >TMath::Pi()) difClusV0C -= TMath::Pi();
  Double_t difClusTPC = TVector2::Phi_0_2pi(Phi-EPTPC);  if(difClusTPC >TMath::Pi()) difClusTPC -= TMath::Pi();

  Double_t DataV0A[5];
  DataV0A[0] = Et;
  DataV0A[1] = M02;
  DataV0A[2] = fCentrality;
  DataV0A[3] = difClusV0A;
  DataV0A[4] = dr;
  fClusterPbV0A->Fill(DataV0A);

  Double_t DataV0C[5];
  DataV0C[0] = Et;
  DataV0C[1] = M02;
  DataV0C[2] = fCentrality;
  DataV0C[3] = difClusV0C;
  DataV0C[4] = dr;
  fClusterPbV0C->Fill(DataV0C);

  Double_t DataTPC[5];
  DataTPC[0] = Et;
  DataTPC[1] = M02;
  DataTPC[2] = fCentrality;
  DataTPC[3] = difClusTPC;
  DataTPC[4] = dr;
  fClusterPbTPC->Fill(DataTPC);
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

  hEvtCount = new TH1F("hEvtCount", " Event Plane", 9, 0.5, 9.5);
  hEvtCount->GetXaxis()->SetBinLabel(1,"All");
  hEvtCount->GetXaxis()->SetBinLabel(2,"Evt");
  hEvtCount->GetXaxis()->SetBinLabel(3,"Trg Class");
  hEvtCount->GetXaxis()->SetBinLabel(4,"Vtx");
  hEvtCount->GetXaxis()->SetBinLabel(5,"Cent");
  hEvtCount->GetXaxis()->SetBinLabel(6,"EPtask");
  hEvtCount->GetXaxis()->SetBinLabel(7,"ClusterTask");
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

  hEPTPCCor  = new TH2F("hEPTPCCor",   "EPTPC  vs cent after PHOS Correct", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0ACor  = new TH2F("hEPV0ACor",   "EPV0A  vs cent after PHOS Correct", 100, 0., 100., 100, 0., TMath::Pi());
  hEPV0CCor  = new TH2F("hEPV0CCor",   "EPV0C  vs cent after PHOS Correct", 100, 0., 100., 100, 0., TMath::Pi());
  fOutput->Add(hEPTPCCor);
  fOutput->Add(hEPV0ACor);
  fOutput->Add(hEPV0CCor);

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

  if (isV1Clus) {
    //  Et   M02  spdcent DeltaPhi    Dr  
    Int_t    bins[5] = {  500, 350,  100,     100,      100}; // binning
    Double_t min[5]  = {  0.0, 0.0,    0,     0.0,      0 }; // min x
    Double_t max[5]  = { 50.0, 3.5,  100,  TMath::Pi(), 0.1}; // max x

    fClusterPbV0A = new THnSparseF("fClusterPbV0A","",5,bins,min,max);
    fClusterPbV0A->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); 
    fClusterPbV0A->GetAxis(1)->SetTitle("M02"); 
    fClusterPbV0A->GetAxis(2)->SetTitle("V0M Centrality");
    fClusterPbV0A->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); 
    fClusterPbV0A->GetAxis(4)->SetTitle("Dr"); 
    fOutput->Add(fClusterPbV0A);

    fClusterPbV0C = new THnSparseF("fClusterPbV0C","",5,bins,min,max);
    fClusterPbV0C->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); 
    fClusterPbV0C->GetAxis(1)->SetTitle("M02"); 
    fClusterPbV0C->GetAxis(2)->SetTitle("V0M Centrality");
    fClusterPbV0C->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); 
    fClusterPbV0C->GetAxis(4)->SetTitle("Dr");
    fOutput->Add(fClusterPbV0C);

    fClusterPbTPC = new THnSparseF("fClusterPbTPC","",5,bins,min,max);
    fClusterPbTPC->GetAxis(0)->SetTitle("Transverse Energy [GeV]"); 
    fClusterPbTPC->GetAxis(1)->SetTitle("M02"); 
    fClusterPbTPC->GetAxis(2)->SetTitle("V0M Centrality");
    fClusterPbTPC->GetAxis(3)->SetTitle("Delta(#phi) [rad]"); 
    fClusterPbTPC->GetAxis(4)->SetTitle("Dr");
    fOutput->Add(fClusterPbTPC);
  }
  
  h2DcosV0A = new TProfile("h2DcosV0A", "cos(Phi) V0r vs Run NUmber", 200, 0., 200.);
  h2DsinV0A = new TProfile("h2DsinV0A", "sin(Phi) V0r vs Run NUmber", 200, 0., 200.);
  h2DcosV0C = new TProfile("h2DcosV0C", "cos(Phi) V0r vs Run NUmber", 200, 0., 200.);
  h2DsinV0C = new TProfile("h2DsinV0C", "sin(Phi) V0r vs Run NUmber", 200, 0., 200.);
  h2DcosTPC = new TProfile("h2DcosTPC", "cos(Phi) V0r vs Run NUmber", 200, 0., 200.);
  h2DsinTPC = new TProfile("h2DsinTPC", "sin(Phi) V0r vs Run NUmber", 200, 0., 200.);
  fOutput->Add(h2DcosV0A);
  fOutput->Add(h2DsinV0A);
  fOutput->Add(h2DcosV0C);
  fOutput->Add(h2DsinV0C);
  fOutput->Add(h2DcosTPC);
  fOutput->Add(h2DsinTPC);

  if (isV1Clus) {
    hM02vsPtA = new TH2F("hM02vsPtA", "M02 vs Et before cut", 5000, 0, 50, 400, 0, 4.);
    hM02vsPtB = new TH2F("hM02vsPtB", "M02 vs Et before cut", 5000, 0, 50, 400, 0, 4.);
    fOutput->Add(hM02vsPtA);
    fOutput->Add(hM02vsPtB);
  }
  hClusDxDZA = new TH2F("hClusDxDZA", "clus Dx vs Dz", 1000, -1., 1., 1000, -1., 1);  
  hClusDxDZB = new TH2F("hClusDxDZB", "clus Dx vs Dz", 1000, -1., 1., 1000, -1., 1);
  fOutput->Add(hClusDxDZA);
  fOutput->Add(hClusDxDZB);
    
  if (!isV1Clus) {
    const Int_t ndims = 4;
    Int_t nMgg=500, nPt=40, nCent=20, nDeltaPhi=315;
    Int_t binsv1[ndims] = {nMgg, nPt, nCent, nDeltaPhi};
    Double_t xmin[ndims] = { 0,   0.,  0,   0.};
    Double_t xmax[ndims] = { 0.5, 20., 100, 3.15};
    fHEPV0A = new THnSparseF("fHEPV0A",   "Flow histogram EPV0A", ndims, binsv1, xmin, xmax);
    fHEPV0C = new THnSparseF("fHEPV0C",   "Flow histogram EPV0C", ndims, binsv1, xmin, xmax);
    fHEPTPC = new THnSparseF("fHEPTPC",   "Flow histogram EPTPC", ndims, binsv1, xmin, xmax);
    fHEPV0A->GetAxis(0)->SetTitle("m_{#gamma#gamma} "); 
    fHEPV0A->GetAxis(1)->SetTitle("p_{T}[GeV]"); 
    fHEPV0A->GetAxis(2)->SetTitle("centrality");
    fHEPV0A->GetAxis(3)->SetTitle("#delta #phi");
    fHEPV0C->GetAxis(0)->SetTitle("m_{#gamma#gamma} "); 
    fHEPV0C->GetAxis(1)->SetTitle("p_{T}[GeV]"); 
    fHEPV0C->GetAxis(2)->SetTitle("centrality");
    fHEPV0C->GetAxis(3)->SetTitle("#delta #phi");
    fHEPTPC->GetAxis(0)->SetTitle("m_{#gamma#gamma} "); 
    fHEPTPC->GetAxis(1)->SetTitle("p_{T}[GeV]"); 
    fHEPTPC->GetAxis(2)->SetTitle("centrality");
    fHEPTPC->GetAxis(3)->SetTitle("#delta #phi");
    fOutput->Add(fHEPV0A);
    fOutput->Add(fHEPV0C);
    fOutput->Add(fHEPTPC);
  }
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskPi0V2::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  hEvtCount->Fill(1);
  // Create pointer to reconstructed event

  AliVEvent *event = InputEvent();
  if (!event) { 
    AliError("Could not retrieve event"); 
    return; 
  }

  // create pointer to event
  TString type = AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()->GetDataType();
  if (!type){
    AliError("Cannot get the event");
    return;
  }

  if (type=="ESD"){
    fESD = dynamic_cast<AliESDEvent*>(event);
    if (!fESD) {
      AliError("Cannot get the ESD event");
      return;
    }
  } else if (type=="AOD"){
    fAOD = dynamic_cast<AliAODEvent*>(event);
    if (!fAOD) {
      AliError("Cannot get the AOD event");
      return;
    }
  } else {
    AliError("Cannot happen");
    return;
  }

  if (fESD) 
    fRunNumber = fESD->GetRunNumber();
  else 
    fRunNumber = fAOD->GetRunNumber();

  fInterRunNumber = ConvertToInternalRunNumber(fRunNumber);

  hEvtCount->Fill(2);
  if (!fTrigClass.IsNull()) {
    TString fired;
    if (fESD) {
      fired = fESD->GetFiredTriggerClasses();
    } else {
      fired = fAOD->GetFiredTriggerClasses();
    }
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
    if (!match)
      return;
  }
  hEvtCount->Fill(3);

  if (isPhosCali) {
    if(fESD){
     if( fRunNumber != fESD->GetRunNumber()) 
      SetFlatteningData();
    } else{
    if( fRunNumber != fAOD->GetRunNumber())
      SetFlatteningData();
    }
  }

  const AliVVertex* fvertex;
  if (fESD){
    fvertex = fESD->GetPrimaryVertex();
  } else {
    fvertex = fAOD->GetPrimaryVertex();
  }

  if (TMath::Abs(fvertex->GetZ())>fVtxCut)
    return;
  Double_t vertex[3] = {fvertex->GetX(), fvertex->GetY(), fvertex->GetZ()};

  hEvtCount->Fill(4);

  fCentrality = event->GetCentrality()->GetCentralityPercentile("CL1"); //spd vertex

  hEvtCount->Fill(5);

  AliEventplane *ep = event->GetEventplane();
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

  fEPV0   = TVector2::Phi_0_2pi(fEPV0);    
  if (fEPV0>TMath::Pi())   
    fEPV0  = fEPV0 - TMath::Pi();
  fEPV0r  = TVector2::Phi_0_2pi(fEPV0r);   
  if (fEPV0r>TMath::Pi())  
    fEPV0r = fEPV0r - TMath::Pi();
  fEPV0A  = TVector2::Phi_0_2pi(fEPV0A);   
  if (fEPV0A>TMath::Pi())  
    fEPV0A = fEPV0A - TMath::Pi();
  fEPV0C  = TVector2::Phi_0_2pi(fEPV0C);   
  if (fEPV0C>TMath::Pi())  
    fEPV0C = fEPV0C - TMath::Pi();
  fEPV0Ar = TVector2::Phi_0_2pi(fEPV0Ar);  
  if (fEPV0Ar>TMath::Pi()) 
    fEPV0Ar = fEPV0Ar - TMath::Pi();
  fEPV0Cr = TVector2::Phi_0_2pi(fEPV0Cr);  
  if (fEPV0Cr>TMath::Pi()) 
    fEPV0Cr = fEPV0Cr - TMath::Pi();
  fEPV0AR4   = TVector2::Phi_0_2pi(fEPV0AR4);    
  if (fEPV0AR4>TMath::Pi())   
    fEPV0AR4  = fEPV0AR4 - TMath::Pi();
  fEPV0AR7   = TVector2::Phi_0_2pi(fEPV0AR7);    
  if (fEPV0AR7>TMath::Pi())   
    fEPV0AR7  = fEPV0AR7 - TMath::Pi();
  fEPV0CR0   = TVector2::Phi_0_2pi(fEPV0CR0);    
  if (fEPV0CR0>TMath::Pi())   
    fEPV0CR0  = fEPV0CR0 - TMath::Pi();
  fEPV0CR3   = TVector2::Phi_0_2pi(fEPV0CR3);    
  if (fEPV0CR3>TMath::Pi())   
    fEPV0CR3  = fEPV0CR3 - TMath::Pi();
  if (fEPTPC != -999. &&  fEPTPC != -1)
    hEPTPC->Fill(fCentrality,  fEPTPC); 
  if (fEPTPCreso!=-1) 
    hresoTPC->Fill(fCentrality, fEPTPCreso);
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

  if (isPhosCali) {
    // PHOS Flattening
    fEPV0A = ApplyFlatteningV0A(fEPV0A, fCentrality); //V0A after Phos flatten
    fEPV0C = ApplyFlatteningV0C(fEPV0C, fCentrality); //V0C after Phos flatten
    fEPTPC = ApplyFlattening(fEPTPC, fCentrality);    //TPC after Phos flatten
  }

  if (!isPhosCali) { 
    if(fESD){
     if( fRunNumber != fESD->GetRunNumber())
      SetFlatteningData();
    } else{
    if( fRunNumber != fAOD->GetRunNumber())
      SetFlatteningData();
    }
    hEPTPCCor->Fill(fCentrality, ApplyFlattening(fEPTPC, fCentrality));
    hEPV0ACor->Fill(fCentrality, ApplyFlatteningV0A(fEPV0A, fCentrality));
    hEPV0CCor->Fill(fCentrality, ApplyFlatteningV0C(fEPV0C, fCentrality));
  } else {
    hEPTPCCor->Fill(fCentrality, fEPTPC);
    hEPV0ACor->Fill(fCentrality, fEPV0A);
    hEPV0CCor->Fill(fCentrality, fEPV0C);
  } 

  hdifV0Ar_V0Cr->Fill(fCentrality, TMath::Cos(2.*(fEPV0Ar - fEPV0Cr)));
  hdifV0A_V0CR0->Fill(fCentrality, TMath::Cos(2.*(fEPV0A - fEPV0CR0)));
  hdifV0A_V0CR3->Fill(fCentrality, TMath::Cos(2.*(fEPV0A - fEPV0CR3)));
  hdifV0ACR0_V0CR3->Fill(fCentrality, TMath::Cos(2*(fEPV0CR0 - fEPV0CR3)));
  hdifV0C_V0AR4->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0AR4)));
  hdifV0C_V0AR7->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0AR7)));
  hdifV0AR4_V0AR7->Fill(fCentrality, TMath::Cos(2*(fEPV0AR4 - fEPV0AR7)));
        
  hdifV0A_V0C->Fill(fCentrality, TMath::Cos(2*(fEPV0A - fEPV0C)));
  if (fEPTPC!=-1 && fEPTPC!=-999.){
    hdifV0A_TPC->Fill(fCentrality, TMath::Cos(2*(fEPV0A - fEPTPC)));
    hdifV0C_TPC->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPTPC)));
  }
  hdifV0C_V0A->Fill(fCentrality, TMath::Cos(2*(fEPV0C - fEPV0A)));

  // Cluster loop for reconstructed event

  //================ for v2 clusterize analysis==============================================
  if (!isV1Clus) {
    if (!fV2ClusName.IsNull() && !fV2Clus) {
      fV2Clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV2ClusName));
      if (!fV2Clus) {
	AliError(Form("%s: Could not retrieve v2 cluster name %s!", GetName(), fV2ClusName.Data()));
	return;
      }
    }
    Int_t nCluster =  fV2Clus->GetEntries(); 
    for (Int_t i=0; i<nCluster; ++i) {
      AliVCluster *c1 = static_cast<AliVCluster*>(fV2Clus->At(i));      
      if (!c1) 
	continue;
      hClusDxDZA->Fill(c1->GetTrackDz(), c1->GetTrackDx());
      if (!c1->IsEMCAL()) 
	continue;
      if (!IsGoodCluster(c1)) 
	continue;
      hClusDxDZB->Fill(c1->GetTrackDz(), c1->GetTrackDx());
      TLorentzVector p1;
      GetMom(p1, c1, vertex);
      for (Int_t j=i+1; j<nCluster; ++j) {
	AliVCluster *c2 = static_cast<AliVCluster*>(fV2Clus->At(j));      
	if (!c2) 
	  continue;
	if (!c2->IsEMCAL()) 
	  continue;
	if (!IsGoodCluster(c2)) 
	  continue;
	TLorentzVector p2;
	GetMom(p2, c2, vertex);
	FillPion(p1, p2, fEPV0A, fEPV0C, fEPTPC);
      }
    }
  }

  //================ for v1 clusterize analysis==============================================
  if (isV1Clus) {
    if (!fV2ClusName.IsNull() && !fV1Clus) {
      fV1Clus = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fV1ClusName));
      if (!fV1Clus) {
	AliError(Form("%s: Could not retrieve v1 cluster name %s!", GetName(), fV1ClusName.Data()));
	return;
      }
    }
    Int_t nClusterV1 = fV1Clus->GetEntries();
    for (Int_t i=0; i<nClusterV1; ++i) {
      AliVCluster *c3 = dynamic_cast<AliVCluster*>(fV1Clus->At(i));      
      if (!c3) 
	continue;
      if (!c3->IsEMCAL()) 
	continue;
      Double_t M02c3 = c3->GetM02();
      Double_t Dxc3  = c3->GetTrackDx();
      Double_t Dzc3  = c3->GetTrackDz(); 

      hClusDxDZA->Fill(Dzc3, Dxc3);
      Float_t clsPosEt[3] = {0,0,0};
      c3->GetPosition(clsPosEt);
      TVector3 clsVec(clsPosEt);
      Double_t Et = c3->E()*TMath::Sin(clsVec.Theta());
      hM02vsPtA->Fill(Et, M02c3);
      if (!IsGoodClusterV1(c3)) 
	continue;
      hM02vsPtB->Fill(Et, M02c3);
      hClusDxDZB->Fill(Dzc3, Dxc3);
      TLorentzVector p3;
      GetMom(p3, c3, vertex);
      FillCluster(p3, fEPV0A, fEPV0C, fEPTPC, c3);
    }
  }

  hEvtCount->Fill(7);

  if (!fTracksName.IsNull() && !fTracks) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Could not retrieve tracks %s!", GetName(), fTracksName.Data())); 
      return;
    }
  }

  Int_t ntracks = fTracks->GetEntries();
  for (Int_t i=0; i<ntracks; ++i){
    AliVTrack *track = static_cast<AliVTrack*>(fTracks->At(i));
    if (!track) 
      continue;
    Double_t tPhi = track->Phi();
    Double_t tPt  = track->Pt();
    Double_t Eta  = track->Eta();
    
    Double_t difTrackV0  = TVector2::Phi_0_2pi(tPhi-fEPV0);   
    if (difTrackV0  >TMath::Pi()) 
      difTrackV0  -= TMath::Pi();
    Double_t difTrackV0A = TVector2::Phi_0_2pi(tPhi-fEPV0A);  
    if (difTrackV0A >TMath::Pi()) 
      difTrackV0A -= TMath::Pi();
    Double_t difTrackV0C = TVector2::Phi_0_2pi(tPhi-fEPV0C);  
    if (difTrackV0C >TMath::Pi()) 
      difTrackV0C -= TMath::Pi();
    Double_t difTrackTPC = TVector2::Phi_0_2pi(tPhi-fEPTPC);  
    if (difTrackTPC >TMath::Pi()) 
      difTrackTPC -= TMath::Pi();
    if (tPhi*TMath::RadToDeg()>80. && tPhi*TMath::RadToDeg()<180. && Eta <0.7 && Eta >(-0.7)){	
      hdifEMC_EPV0->Fill(fCentrality, difTrackV0, tPt);
      hdifEMC_EPV0A->Fill(fCentrality, difTrackV0A, tPt);
      hdifEMC_EPV0C->Fill(fCentrality, difTrackV0C, tPt);
    } else {
      hdifout_EPV0->Fill(fCentrality, difTrackV0, tPt);
      hdifout_EPV0A->Fill(fCentrality, difTrackV0A, tPt);
      hdifout_EPV0C->Fill(fCentrality, difTrackV0C, tPt);
    }
    hdifful_EPV0->Fill(fCentrality,    difTrackV0, tPt);
    hdifful_EPV0A->Fill(fCentrality,   difTrackV0A, tPt);
    hdifful_EPV0C->Fill(fCentrality,   difTrackV0C, tPt);
  } 
  hEvtCount->Fill(8);

  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
}

//____________________________________________________________________
Int_t AliAnalysisTaskPi0V2::ConvertToInternalRunNumber(Int_t n)
{
  switch(n) {
    case  170593 : return 179;
    case  170572 : return 178;
    case  170556 : return 177;
    case  170552 : return 176;
    case  170546 : return 175;
    case  170390 : return 174;
    case  170389 : return 173;
    case  170388 : return 172;
    case  170387 : return 171;
    case  170315 : return 170;
    case  170313 : return 169;
    case  170312 : return 168;
    case  170311 : return 167;
    case  170309 : return 166;
    case  170308 : return 165;
    case  170306 : return 164;
    case  170270 : return 163;
    case  170269 : return 162;
    case  170268 : return 161;
    case  170267 : return 160;
    case  170264 : return 159;
    case  170230 : return 158;
    case  170228 : return 157;
    case  170208 : return 156;
    case  170207 : return 155;
    case  170205 : return 154;
    case  170204 : return 153;
    case  170203 : return 152;
    case  170195 : return 151;
    case  170193 : return 150;
    case  170163 : return 149;
    case  170162 : return 148;
    case  170159 : return 147;
    case  170155 : return 146;
    case  170152 : return 145;
    case  170091 : return 144;
    case  170089 : return 143;
    case  170088 : return 142;
    case  170085 : return 141;
    case  170084 : return 140;
    case  170083 : return 139;
    case  170081 : return 138;
    case  170040 : return 137;
    case  170038 : return 136;
    case  170036 : return 135;
    case  170027 : return 134;
    case  169981 : return 133;
    case  169975 : return 132;
    case  169969 : return 131;
    case  169965 : return 130;
    case  169961 : return 129;
    case  169956 : return 128;
    case  169926 : return 127;
    case  169924 : return 126;
    case  169923 : return 125;
    case  169922 : return 124;
    case  169919 : return 123;
    case  169918 : return 122;
    case  169914 : return 121;
    case  169859 : return 120;
    case  169858 : return 119;
    case  169855 : return 118;
    case  169846 : return 117;
    case  169838 : return 116;
    case  169837 : return 115;
    case  169835 : return 114;
    case  169683 : return 113;
    case  169628 : return 112;
    case  169591 : return 111;
    case  169590 : return 110;
    case  169588 : return 109;
    case  169587 : return 108;
    case  169586 : return 107;
    case  169584 : return 106;
    case  169557 : return 105;
    case  169555 : return 104;
    case  169554 : return 103;
    case  169553 : return 102;
    case  169550 : return 101;
    case  169515 : return 100;
    case  169512 : return 99;
    case  169506 : return 98;
    case  169504 : return 97;
    case  169498 : return 96;
    case  169475 : return 95;
    case  169420 : return 94;
    case  169419 : return 93;
    case  169418 : return 92;
    case  169417 : return 91;
    case  169415 : return 90;
    case  169411 : return 89;
    case  169238 : return 88;
    case  169236 : return 87;
    case  169167 : return 86;
    case  169160 : return 85;
    case  169156 : return 84;
    case  169148 : return 83;
    case  169145 : return 82;
    case  169144 : return 81;
    case  169143 : return 80;
    case  169138 : return 79;
    case  169099 : return 78;
    case  169094 : return 77;
    case  169091 : return 76;
    case  169045 : return 75;
    case  169044 : return 74;
    case  169040 : return 73;
    case  169035 : return 72;
    case  168992 : return 71;
    case  168988 : return 70;
    case  168984 : return 69;
    case  168826 : return 68;
    case  168777 : return 67;
    case  168514 : return 66;
    case  168512 : return 65;
    case  168511 : return 64;
    case  168467 : return 63;
    case  168464 : return 62;
    case  168461 : return 61;
    case  168460 : return 60;
    case  168458 : return 59;
    case  168362 : return 58;
    case  168361 : return 57;
    case  168356 : return 56;
    case  168342 : return 55;
    case  168341 : return 54;
    case  168325 : return 53;
    case  168322 : return 52;
    case  168318 : return 51;
    case  168311 : return 50;
    case  168310 : return 49;
    case  168213 : return 48;
    case  168212 : return 47;
    case  168208 : return 46;
    case  168207 : return 45;
    case  168206 : return 44;
    case  168205 : return 43;
    case  168204 : return 42;
    case  168203 : return 41;
    case  168181 : return 40;
    case  168177 : return 39;
    case  168175 : return 38;
    case  168173 : return 37;
    case  168172 : return 36;
    case  168171 : return 35;
    case  168115 : return 34;
    case  168108 : return 33;
    case  168107 : return 32;
    case  168105 : return 31;
    case  168104 : return 30;
    case  168103 : return 29;
    case  168076 : return 28;
    case  168069 : return 27;
    case  168068 : return 26;
    case  168066 : return 25;
    case  167988 : return 24;
    case  167987 : return 23;
    case  167986 : return 22;
    case  167985 : return 21;
    case  167921 : return 20;
    case  167920 : return 19;
    case  167915 : return 18;
    case  167909 : return 17;
    case  167903 : return 16;
    case  167902 : return 15;
    case  167818 : return 14;
    case  167814 : return 13;
    case  167813 : return 12;
    case  167808 : return 11;
    case  167807 : return 10;
    case  167806 : return 9;
    case  167713 : return 8;
    case  167712 : return 7;
    case  167711 : return 6;
    case  167706 : return 5;
    case  167693 : return 4;
    case  166532 : return 3;
    case  166530 : return 2;
    case  166529 : return 1;

    default : return 199;
  }
}

//_______________________________________________________________________
void AliAnalysisTaskPi0V2::FillEPQA()
{
  h2DcosV0A->Fill(fInterRunNumber, TMath::Cos(fEPV0A));
  h2DsinV0A->Fill(fInterRunNumber, TMath::Sin(fEPV0A));
  h2DcosV0C->Fill(fInterRunNumber, TMath::Cos(fEPV0C));
  h2DsinV0C->Fill(fInterRunNumber, TMath::Sin(fEPV0C));
  h2DcosTPC->Fill(fInterRunNumber, TMath::Cos(fEPTPC));
  h2DsinTPC->Fill(fInterRunNumber, TMath::Sin(fEPTPC));
}

//_________________________________________________________________________________
void AliAnalysisTaskPi0V2::SetFlatteningData()
{
  //Read objects with flattening parameters
  AliOADBContainer flatContainer("phosFlat");
  flatContainer.InitFromFile(fEPcalibFileName.Data(),"phosFlat");
  TObjArray *maps = (TObjArray*)flatContainer.GetObject(fRunNumber,"phosFlat");
  if (!maps) {
    AliError(Form("Can not read Flattening for run %d. \n From file >%s<\n",fRunNumber,fEPcalibFileName.Data())) ;    
  } else {
    AliInfo(Form("Setting PHOS flattening with name %s \n",maps->GetName())) ;
    AliEPFlattener * h = (AliEPFlattener*)maps->At(0) ;  
    if(fTPCFlat) delete fTPCFlat ;
    fTPCFlat = new AliEPFlattener();
    fTPCFlat = h ;
    h = (AliEPFlattener*)maps->At(1);  
    if(fV0AFlat) delete fV0AFlat ;
    fV0AFlat = new AliEPFlattener();
    fV0AFlat = h ;
    h = (AliEPFlattener*)maps->At(2);  
    if(fV0CFlat) delete fV0CFlat ;
    fV0CFlat = new AliEPFlattener();
    fV0CFlat = h;
  }    
}
 
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0V2::ApplyFlattening(Double_t phi, Double_t c)
{
  if(fTPCFlat)
    return fTPCFlat->MakeFlat(phi,c);
  return phi;
}

//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0V2::ApplyFlatteningV0A(Double_t phi, Double_t c)
{
  if(fV0AFlat)
    return fV0AFlat->MakeFlat(phi,c);
  return phi;
}

//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0V2::ApplyFlatteningV0C(Double_t phi, Double_t c){
 
  if(fV0CFlat)
    return fV0CFlat->MakeFlat(phi,c);
  return phi;
}

//________________________________________________________________________
void AliAnalysisTaskPi0V2::Terminate(Option_t *) 
{
}
