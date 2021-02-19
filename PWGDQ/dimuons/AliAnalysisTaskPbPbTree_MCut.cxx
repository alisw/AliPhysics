/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: AliAnalysisTaskPbPbTree_MCut.cxx $ */

//-----------------------------------------------------------------------------
// Analysis task to compute muon/dimuon kinematic distributions
// The output is a list of histograms.
// The macro class can run on AOD or in the train with the ESD filter.
// R. Arnaldi, Luca Micheletti
//
//-----------------------------------------------------------------------------

//#ifndef AliAnalysisTaskPbPbTree_MCut_CXX
//#define AliAnalysisTaskPbPbTree_MCut_CXX

// ROOT includes
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TList.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include "TGrid.h"
#include "TRandom3.h"

#include "AliInputEventHandler.h"
#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"
#include "AliMultSelection.h"
#include "AliOADBContainer.h"
#include "AliAODv0.h"

#include "AliAnalysisTaskPbPbTree_MCut.h"
#include "AliAODZDC.h"
#include "AliTriggerAnalysis.h"
#include "AliVMultiplicity.h"
#include "AliAODTracklets.h"

// STL includes
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

TH1D*        fMultV0;                        // profile from V0 multiplicity
TH1D*        fQx2mV0A[18];                    // <Qxn> V0A
TH1D*        fQy2mV0A[18];                    // <Qyn> V0A
TH1D*        fQx2sV0A[18];                    // sigma Qxn V0A
TH1D*        fQy2sV0A[18];                    // sigma Qyn V0A

TH1D*        fQx2mV0C[18];                    // <Qxn> V0C
TH1D*        fQy2mV0C[18];                    // <Qyn> V0C
TH1D*        fQx2sV0C[18];                    // sigma Qxn V0C
TH1D*        fQy2sV0C[18];                    // sigma Qyn V0C

TH1D*        fQx2mTrk[18];                    // <Qxn> tracklets
TH1D*        fQy2mTrk[18];                    // <Qyn> tracklets
TH1D*        fQx2sTrk[18];                    // sigma Qxn tracklets
TH1D*        fQy2sTrk[18];                    // sigma Qyn tracklets


TH1D*        fQx3mV0A[18];                    // <Qxn> V0A
TH1D*        fQy3mV0A[18];                    // <Qyn> V0A
TH1D*        fQx3sV0A[18];                    // sigma Qxn V0A
TH1D*        fQy3sV0A[18];                    // sigma Qyn V0A

TH1D*        fQx3mV0C[18];                    // <Qxn> V0C
TH1D*        fQy3mV0C[18];                    // <Qyn> V0C
TH1D*        fQx3sV0C[18];                    // sigma Qxn V0C
TH1D*        fQy3sV0C[18];                    // sigma Qyn V0C

TH1D*        fQx3mTrk[18];                    // <Qxn> tracklets
TH1D*        fQy3mTrk[18];                    // <Qyn> tracklets
TH1D*        fQx3sTrk[18];                    // sigma Qxn tracklets
TH1D*        fQy3sTrk[18];                    // sigma Qyn tracklets

Double_t CostHE_PbPb(AliAODTrack*, AliAODTrack*);
Double_t CostCS_PbPb(AliAODTrack*, AliAODTrack*);
Double_t PhiHE_PbPb(AliAODTrack*, AliAODTrack*);
Double_t PhiCS_PbPb(AliAODTrack*, AliAODTrack*);
Double_t CostEPnB_PbPb(AliAODTrack*, AliAODTrack*, Double_t);

ClassImp(AliAnalysisTaskPbPbTree_MCut)
//__________________________________________________________________________
AliAnalysisTaskPbPbTree_MCut::AliAnalysisTaskPbPbTree_MCut() :
  AliAnalysisTaskSE(),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fMassCut(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fPercentV0M(0x0),
  fPercentCL0(0x0),
  fPercentCL1(0x0),
  fPercentV0A(0x0),
  fPercentV0C(0x0),
  fPercentZNA(0x0),
  fPercentZNC(0x0),
  fIsPhysSelected(0x0),
  fPsi2Trkl(0x0),
  fPsi3Trkl(0x0),
  fPsi2RP(0x0),
  fAODEvent(0x0),
//  fTrigClass(0x0),
  finpmask(0),
  fNTracks(0x0)
{
  //
  //Default ctor
  //
  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "StandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;
  }
  for(Int_t i=0; i<400;i++){
    fDimuPt[i]=999.;
    fDimuPx[i]=999.;
    fDimuPy[i]=999.;
    fDimuPz[i]=999.;
    fDimuY[i]=999.;
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    fDimuCostHE[i] = 999;
    fDimuPhiHE[i] = 999;
    fDimuCostCS[i] = 999;
    fDimuPhiCS[i] = 999;
    fDimuCostEPnB[i] = 999;
    fDimuCostRPnB[i] = 999;
    fDimuPhi[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  }
}

//__________________________________________________________________________
AliAnalysisTaskPbPbTree_MCut::AliAnalysisTaskPbPbTree_MCut(const char *name) :
  AliAnalysisTaskSE(name),
  fOutputTree(0x0),
  fNevt(0x0),
  fBeamEnergy(0.),
  fMassCut(0.),
  fkAnalysisType(0x0),
  fPeriod(0x0),
  fCountTotEv(0x0),
  fCountTrigger(0x0),
  fCountCINT7(0x0),
  fCountCMUL7(0x0),
  fCountCMLL7(0x0),
  fCountCMSL7(0x0),
  fCountCMSH7(0x0),
  fNMuons(0x0),
  fNTracklets(0x0),
  fNContributors(0x0),
  fNDimu(0x0),
  fPercentV0M(0x0),
  fPercentCL0(0x0),
  fPercentCL1(0x0),
  fPercentV0A(0x0),
  fPercentV0C(0x0),
  fPercentZNA(0x0),
  fPercentZNC(0x0),
  fIsPhysSelected(0x0),
  fPsi2Trkl(0x0),
  fPsi3Trkl(0x0),
  fPsi2RP(0x0),
  fAODEvent(0x0),
//  fTrigClass(0x0),
  finpmask(0),
  fNTracks(0x0)
{
  //
  // Constructor. Initialization of Inputs and Outputs
  //
  Info("AliAnalysisTaskPbPbTree_MCut","Calling Constructor");

  fMuonTrackCuts = new AliMuonTrackCuts("StandardMuonTrackCuts", "TestStandardMuonTrackCuts");
  fMuonTrackCuts->SetFilterMask(AliMuonTrackCuts::kMuPdca);
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);

  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999;
  }
  for(Int_t i=0; i<400;i++){
    fDimuPt[i]=999.;
    fDimuPx[i]=999.;
    fDimuPy[i]=999.;
    fDimuPz[i]=999.;
    fDimuY[i]=999.;
    fDimuMass[i]=999.;
    fDimuCharge[i]=999;
    fDimuMatch[i]=999;
    fDimuCostHE[i] = 999;
    fDimuPhiHE[i] = 999;
    fDimuCostCS[i] = 999;
    fDimuPhiCS[i] = 999;
    fDimuCostEPnB[i] = 999;
    fDimuCostRPnB[i] = 999;
    fDimuPhi[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  }

  DefineOutput(1,TTree::Class());
  DefineOutput(2,TH1D::Class());
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_MCut& AliAnalysisTaskPbPbTree_MCut::operator=(const AliAnalysisTaskPbPbTree_MCut& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c) ;
  }
  return *this;
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_MCut::AliAnalysisTaskPbPbTree_MCut(const AliAnalysisTaskPbPbTree_MCut& c) :
  AliAnalysisTaskSE(c),
  fOutputTree(c.fOutputTree),
  fNevt(c.fNevt),
  fBeamEnergy(c.fBeamEnergy),
  fMassCut(c.fMassCut),
  fkAnalysisType(c.fkAnalysisType),
  fPeriod(c.fPeriod),
  fCountTotEv(c.fCountTotEv),
  fCountTrigger(c.fCountTrigger),
  fCountCINT7(c.fCountCINT7),
  fCountCMUL7(c.fCountCMUL7),
  fCountCMLL7(c.fCountCMLL7),
  fCountCMSL7(c.fCountCMSL7),
  fCountCMSH7(c.fCountCMSH7),
  fNMuons(c.fNMuons),
  fNTracklets(c.fNTracklets),
  fNContributors(c.fNContributors),
  fNDimu(c.fNDimu),
  fPercentV0M(c.fPercentV0M),
  fPercentCL0(c.fPercentCL0),
  fPercentCL1(c.fPercentCL1),
  fPercentV0A(c.fPercentV0A),
  fPercentV0C(c.fPercentV0C),
  fPercentZNA(c.fPercentZNA),
  fPercentZNC(c.fPercentZNC),
  fIsPhysSelected(c.fIsPhysSelected),
  fPsi2Trkl(c.fPsi2Trkl),
  fPsi3Trkl(c.fPsi3Trkl),
  fPsi2RP(c.fPsi2RP),
  fAODEvent(c.fAODEvent),
//  fTrigClass(c.fTrigClass),
  finpmask(c.finpmask),
  fMuonTrackCuts(c.fMuonTrackCuts),
  fNTracks(c.fNTracks)

 {
  //
  // Copy Constructor
  //
}

//___________________________________________________________________________
AliAnalysisTaskPbPbTree_MCut::~AliAnalysisTaskPbPbTree_MCut() {
  //
  //destructor
  //
  Info("~AliAnalysisTaskPbPbTree_MCut","Calling Destructor");
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis) delete fOutputTree;
}

//___________________________________________________________________________
void AliAnalysisTaskPbPbTree_MCut::NotifyRun()
{
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->SetRun((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()));
}

//___________________________________________________________________________
void AliAnalysisTaskPbPbTree_MCut::UserCreateOutputObjects(){

  if (fOutputTree) return;

  OpenFile(1,"RECREATE");
  fOutputTree = new TTree("PbPbTree","Data Tree");

  fOutputTree->Branch("FiredTriggerClasses",fTrigClass,"FiredTriggerClasses/C");
  fOutputTree->Branch("inpmask",&finpmask,"inpmask/i");

  fOutputTree->Branch("NMuons",&fNMuons,"NMuons/I");
  fOutputTree->Branch("Vertex",fVertex,"Vertex[3]/D");
  fOutputTree->Branch("PercentV0M",&fPercentV0M,"PercentV0M/F");
  fOutputTree->Branch("PercentCL0",&fPercentCL0,"PercentCL0/F");
  fOutputTree->Branch("PercentCL1",&fPercentCL1,"PercentCL1/F");
  fOutputTree->Branch("PercentV0A",&fPercentV0A,"PercentV0A/F");
  fOutputTree->Branch("PercentV0C",&fPercentV0C,"PercentV0C/F");
  fOutputTree->Branch("PercentZNA",&fPercentZNA,"PercentZNA/F");
  fOutputTree->Branch("PercentZNC",&fPercentZNC,"PercentZNC/F");
  fOutputTree->Branch("NTracks",&fNTracks,"NTracks/I");

  fOutputTree->Branch("Pt",fPt,"Pt[NMuons]/D");
  fOutputTree->Branch("E",fE,"E[NMuons]/D");
  fOutputTree->Branch("Px",fPx,"Px[NMuons]/D");
  fOutputTree->Branch("Py",fPy,"Py[NMuons]/D");
  fOutputTree->Branch("Pz",fPz,"Pz[NMuons]/D");
  fOutputTree->Branch("Y",fY,"Y[NMuons]/D");
  fOutputTree->Branch("Eta",fEta,"Eta[NMuons]/D");
  fOutputTree->Branch("Phi",fPhi,"Phi[NMuons]/D");
  fOutputTree->Branch("MatchTrig",fMatchTrig,"MatchTrig[NMuons]/I");
  fOutputTree->Branch("TrackChi2",fTrackChi2,"TrackChi2[NMuons]/D");
  fOutputTree->Branch("MatchTrigChi2",fMatchTrigChi2,"MatchTrigChi2[NMuons]/D");
  fOutputTree->Branch("Charge",fCharge,"Charge[NMuons]/I");
  fOutputTree->Branch("RAtAbsEnd",fRAtAbsEnd,"RAtAbsEnd[NMuons]/D");
  fOutputTree->Branch("pDCA",fpDCA,"pDCA[NMuons]/I");

  fOutputTree->Branch("NDimu",&fNDimu,"NDimu/I");
  fOutputTree->Branch("DimuMu",fDimuMu,"DimuMu[NDimu][2]/I");
  fOutputTree->Branch("DimuPt",fDimuPt,"DimuPt[NDimu]/D");
  fOutputTree->Branch("DimuPx",fDimuPx,"DimuPx[NDimu]/D");
  fOutputTree->Branch("DimuPy",fDimuPy,"DimuPy[NDimu]/D");
  fOutputTree->Branch("DimuPz",fDimuPz,"DimuPz[NDimu]/D");
  fOutputTree->Branch("DimuY",fDimuY,"DimuY[NDimu]/D");
  fOutputTree->Branch("DimuMass",fDimuMass,"DimuMass[NDimu]/D");
  fOutputTree->Branch("DimuCharge",fDimuCharge,"DimuCharge[NDimu]/I");
  fOutputTree->Branch("DimuMatch",fDimuMatch,"DimuMatch[NDimu]/I");
  fOutputTree->Branch("DimuCostHE",fDimuCostHE,"DimuCostHE[NDimu]/D");
  fOutputTree->Branch("DimuPhiHE",fDimuPhiHE,"DimuPhiHE[NDimu]/D");
  fOutputTree->Branch("DimuCostCS",fDimuCostCS,"DimuCostCS[NDimu]/D");
  fOutputTree->Branch("DimuPhiCS",fDimuPhiCS,"DimuPhiCS[NDimu]/D");
  fOutputTree->Branch("DimuCostEPnB",fDimuCostEPnB,"DimuCostEPnB[NDimu]/D");
  fOutputTree->Branch("DimuCostRPnB",fDimuCostRPnB,"DimuCostRPnB[NDimu]/D");
  fOutputTree->Branch("IsPhysSelected",&fIsPhysSelected,"IsPhysSelected/O");

  fOutputTree -> Branch("Psi2Trkl",&fPsi2Trkl,"Psi2Trkl/D");
  fOutputTree -> Branch("Psi3Trkl",&fPsi3Trkl,"Psi3Trkl/D");
  fOutputTree -> Branch("Psi2RP",&fPsi2RP,"Psi2RP/D");
  fOutputTree -> Branch("DimuPhiEP",fDimuPhi,"DimuPhiEP[NDimu]/D");

  fOutputTree->ls();

  PostData(1,fOutputTree);

  //trigger summary
  fhNEv = new TH1D("fhNEv","hNEv",8,0.,8.);
  TString namelabel1[8]={"TotEv","CINT7","CMUL7","CMLL7","CMSL7","CMSH7","CINT7_CENT","CINT7ZAC_CENT"};
  for(int k=0;k<8;k++) fhNEv->GetXaxis()->SetBinLabel(k+1,namelabel1[k]);

  PostData(2,fhNEv);

}

//_________________________________________________
void AliAnalysisTaskPbPbTree_MCut::UserExec(Option_t *)
{

  fNTracks=0;
  fNMuons=0;
  fNTracklets=-1;
  fNContributors=-1;
  fNDimu=0;
  fPercentV0M=-1.;
  fPercentCL0=-1.;
  fPercentCL1=-1.;
  fPercentV0A=-1.;
  fPercentV0C=-1.;
  fPercentZNA=-1.;
  fPercentZNC=-1.;
  fPsi2Trkl=-999.;
  fPsi3Trkl=-999.;
  fPsi2RP=-999.;
  fVertex[0]=999.; fVertex[1]=999.; fVertex[2]=999.;
  for(Int_t i=0; i<1500;i++){
    fPt[i]=999.;
    fE[i]=999.;
    fPx[i]=999;
    fPy[i]=999;
    fPz[i]=999;
    fY[i]=999.;
    fEta[i]=999.;
    fPhi[i]=999.;
    fMatchTrig[i]=999.;
    fTrackChi2[i]=999.;
    fMatchTrigChi2[i]=999.;
    fDCA[i]=999.;
    fCharge[i]=999;
    fRAtAbsEnd[i]=999;
    fpDCA[i] = 999.;
  }
  for(Int_t i=0; i<400;i++){
    fDimuPt[i]=999.;
    fDimuPx[i]=999.;
    fDimuPy[i]=999.;
    fDimuPz[i]=999.;
    fDimuY[i]=999.;
    fDimuMass[i]=999.;
    fDimuCharge[i]=999.;
    fDimuMatch[i]=0;
    fDimuCostHE[i] = 999.;
    fDimuPhiHE[i] = 999.;
    fDimuCostCS[i] = 999.;
    fDimuPhiCS[i] = 999.;
    fDimuCostEPnB[i] = 999;
    fDimuCostRPnB[i] = 999;
    fDimuPhi[i] = 999;
    for(Int_t k=0;k<2;k++) fDimuMu[i][k]=999;
  }

//
// Execute analysis for current event
//
  fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if ( ! fAODEvent ) {
    AliError ("AOD event not found. Nothing done!");
    return;
  }

   AliAODHeader *aodheader=dynamic_cast<AliAODHeader*>(fAODEvent->GetHeader());
    TString firedtrigger = aodheader->GetFiredTriggerClasses();
    sprintf(fTrigClass,"%s",firedtrigger.Data());

   finpmask = aodheader->GetL0TriggerInputs();

  //   to apply physics selection
    UInt_t fSelectMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    fIsPhysSelected = fSelectMask & (AliVEvent::kMuonUnlikeLowPt7 | AliVEvent::kMuonLikeLowPt7 | AliVEvent::kMuonSingleLowPt7 |   AliVEvent::kMuonSingleHighPt7  | AliVEvent::kINT7inMUON  | AliVEvent::kINT7);


   Bool_t TriggerSelected=kFALSE;
   Bool_t TriggerSelected_CINT7=kFALSE;
   Bool_t TriggerSelected_CMUL7=kFALSE;
   Bool_t TriggerSelected_CMLL7=kFALSE;
   Bool_t TriggerSelected_CMSL7=kFALSE;
   Bool_t TriggerSelected_CMSH7=kFALSE;
   Bool_t TriggerSelected_CINT7_CENT=kFALSE;
   Bool_t TriggerSelected_CINT7ZAC_CENT=kFALSE;

  if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
  else TriggerSelected = kFALSE;
  if(firedtrigger.Contains("CINT7-B-NOPF-MUFAST")) TriggerSelected_CINT7 = kTRUE;  // o devo guardare CENT??? o quello con ZN?
										  // guardare quello con ZN anche per il calcolo di fnorm cint7_zac
									       // check se fermarsi a 80%
  else TriggerSelected_CINT7 = kFALSE;
  if(firedtrigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected_CMUL7 = kTRUE;
  else TriggerSelected_CMUL7 = kFALSE;
  if(firedtrigger.Contains("CMLL7-B-NOPF-MUFAST")) TriggerSelected_CMLL7 = kTRUE;
  else TriggerSelected_CMLL7 = kFALSE;
  if(firedtrigger.Contains("CMSL7-B-NOPF-MUFAST")) TriggerSelected_CMSL7 = kTRUE;
  else TriggerSelected_CMSL7 = kFALSE;
  if(firedtrigger.Contains("CMSH7-B-NOPF-MUFAST")) TriggerSelected_CMSH7 = kTRUE;
  else TriggerSelected_CMSH7 = kFALSE;
  if(firedtrigger.Contains("CINT7-B-NOPF-CENT")) TriggerSelected_CINT7_CENT = kTRUE;
  if(firedtrigger.Contains("CINT7ZAC-B-NOPF-CENT")) TriggerSelected_CINT7ZAC_CENT = kTRUE;

  Double_t DeltaCh=0.;
  fhNEv->Fill(0.+DeltaCh);
  if (TriggerSelected_CINT7) fhNEv->Fill(1.+DeltaCh);
  if (TriggerSelected_CMUL7) fhNEv->Fill(2.+DeltaCh);
  if (TriggerSelected_CMLL7) fhNEv->Fill(3.+DeltaCh);
  if (TriggerSelected_CMSL7) fhNEv->Fill(4.+DeltaCh);
  if (TriggerSelected_CMSH7) fhNEv->Fill(5.+DeltaCh);
  if (TriggerSelected_CINT7_CENT) fhNEv->Fill(16.+DeltaCh);
  if (TriggerSelected_CINT7ZAC_CENT) fhNEv->Fill(17.+DeltaCh);

  // Evaluation of Event-Plane variables
  Int_t run = fAODEvent->GetRunNumber();
  if(run != fRun){
    // Load the calibrations run dependent
    OpenInfoCalbration(run);
    fRun = run;
  }

  Float_t zvtx = GetVertex(fAODEvent);
  fVtxCut = 14.;

  if(zvtx< -990){

  } else {

    if (TMath::Abs(zvtx) < fVtxCut){
      //Centrality
      Float_t v0Centr    = -100.;
      Float_t cl1Centr   = -100.;
      Float_t cl0Centr   = -100.;

      AliMultSelection* MultSelection = 0x0;
      MultSelection = (AliMultSelection*) fAODEvent->FindListObject("MultSelection");
      if( !MultSelection) {
          AliWarning("AliMultSelection object not found!");
          return;
      } else {
          v0Centr  = MultSelection -> GetMultiplicityPercentile("V0M");
          cl1Centr = MultSelection -> GetMultiplicityPercentile("CL1");
          cl0Centr = MultSelection -> GetMultiplicityPercentile("CL0");
      }

      if (v0Centr >= 90. || cl1Centr >= 90 || cl0Centr >= 90)
          return;

        Int_t iCentV0 = Int_t(v0Centr);
        if (iCentV0 >= 90)
            return;


        Int_t iCentSPD = Int_t(cl1Centr);
        if (iCentSPD >= 90)
            return;


        Short_t zvt = GetVertexZ(zvtx);
        if (zvt < 0)
            return;

        AliAODTracklets* aodTrkl = (AliAODTracklets*)fAODEvent->GetTracklets();
        Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();

        //Tracklets
        Double_t Qx2trkcut = 0, Qy2trkcut = 0;
        Double_t Qx3trkcut = 0, Qy3trkcut = 0;
        Int_t multtrkcut = 0;


        for (Int_t it = 0; it < nITSTrkls; it++){
            Double_t corPhi = CalcCorPhi(aodTrkl->GetPhi(it), aodTrkl->GetDeltaPhi(it));

            Qx2trkcut += TMath::Cos(2.*corPhi);
            Qy2trkcut += TMath::Sin(2.*corPhi);

            Qx3trkcut += TMath::Cos(3.*corPhi);
            Qy3trkcut += TMath::Sin(3.*corPhi);

            multtrkcut++;
        }

        Double_t Qytr2Cor = (Qy2trkcut - fQy2mTrk[zvt]->GetBinContent(iCentV0+1))/fQy2sTrk[zvt]->GetBinContent(iCentV0+1);
        Double_t Qxtr2Cor = (Qx2trkcut - fQx2mTrk[zvt]->GetBinContent(iCentV0+1))/fQx2sTrk[zvt]->GetBinContent(iCentV0+1);

        Double_t Qytr3Cor = (Qy3trkcut - fQy3mTrk[zvt]->GetBinContent(iCentV0+1))/fQy3sTrk[zvt]->GetBinContent(iCentV0+1);
        Double_t Qxtr3Cor = (Qx3trkcut - fQx3mTrk[zvt]->GetBinContent(iCentV0+1))/fQx3sTrk[zvt]->GetBinContent(iCentV0+1);

        Double_t psi2Trkl = TMath::ATan2(Qytr2Cor, Qxtr2Cor)/2.;
        Double_t psi3Trkl = TMath::ATan2(Qytr3Cor, Qxtr3Cor)/3.;

        fPsi2Trkl = psi2Trkl;
        fPsi3Trkl = psi3Trkl;
        fPsi2RP = (-TMath::Pi()/2.) + TMath::Pi()*(gRandom -> Rndm());
    }
  }
  // centrality determination
  // https://twiki.cern.ch/twiki/bin/viewauth/ALICE/PACentStudiesRun2
  Float_t PercV0M = 300;
  Float_t PercCL0 = 300;
  Float_t PercCL1 = 300;
  Float_t PercV0A = 300;
  Float_t PercV0C = 300;
  Float_t PercZNA = 300;
  Float_t PercZNC = 300;

  AliMultSelection *MultSelection = 0x0;
  MultSelection = (AliMultSelection * ) fAODEvent -> FindListObject("MultSelection");
  if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
  }else{
        PercV0M = MultSelection -> GetMultiplicityPercentile("V0M"); //second argument in kFALSE
        PercCL0 = MultSelection -> GetMultiplicityPercentile("CL0");
        PercCL1 = MultSelection -> GetMultiplicityPercentile("CL1");
        PercV0A = MultSelection -> GetMultiplicityPercentile("V0A");
        PercV0C = MultSelection -> GetMultiplicityPercentile("V0C");
        PercZNA = MultSelection -> GetMultiplicityPercentile("ZNA");
        PercZNC = MultSelection -> GetMultiplicityPercentile("ZNC");
  }
  fPercentV0M = PercV0M;
  fPercentCL0 = PercCL0;
  fPercentCL1 = PercCL1;
  fPercentV0A = PercV0A;
  fPercentV0C = PercV0C;
  fPercentZNA = PercZNA;
  fPercentZNC = PercZNC;


  AliAODVertex *PrimVertex =  fAODEvent->GetPrimaryVertex();
  fVertex[0]=PrimVertex->GetX();
  fVertex[1]=PrimVertex->GetY();
  fVertex[2]=PrimVertex->GetZ();

  //------------------------------------------------------
  // build dimuon object starting from single muons
  // Keep only dimuons with
  // -- mass > mcut
  // -- 2.5<y<4
  // -- matchtrigger = 2
  //------------------------------------------------------

   Int_t numdimu = 0;
   Int_t nummu = 0;
   Int_t ntracks = fAODEvent->GetNumberOfTracks();
   if(ntracks!=0) {
    fNTracks = ntracks;

   Int_t LabelOld1[1500];
   Int_t LabelOld2[1500];
   Bool_t GoodMuon[1500]={kFALSE};

   //if(ntracks>1500) return;    //skip events with a huge number of tracks --> for pp collisions

   for (Int_t i=0;i<ntracks;i++){

     AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
     if (!mu0->IsMuonTrack()) continue;

     for(Int_t j=i+1;j<ntracks;j++){

       AliAODTrack *mu1=(AliAODTrack*) fAODEvent->GetTrack(j);
       if (!mu1->IsMuonTrack()) continue;

       AliAODDimuon *dimu=new AliAODDimuon(mu0,mu1);

       if(dimu->Mass()>fMassCut && (mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) && dimu->Y()>-4 && dimu->Y()<-2.5){

	 fDimuMass[numdimu] = dimu->Mass();
	 fDimuPt[numdimu] = dimu->Pt();
	 fDimuPx[numdimu] = dimu->Px();
	 fDimuPy[numdimu] = dimu->Py();
	 fDimuPz[numdimu] = dimu->Pz();
	 fDimuY[numdimu] = dimu->Y();
	 fDimuCharge[numdimu]= dimu->Charge();
	 if(mu0->GetMatchTrigger()>1 && mu1->GetMatchTrigger()>1) fDimuMatch[numdimu]=2;
	 else fDimuMatch[numdimu]=-1;
   fDimuPhi[numdimu] = dimu->Phi();

	 fDimuCostHE[numdimu] = CostHE_PbPb(mu0,mu1);
	 fDimuPhiHE[numdimu] = PhiHE_PbPb(mu0,mu1);
	 fDimuCostCS[numdimu] = CostCS_PbPb(mu0,mu1);
	 fDimuPhiCS[numdimu] = PhiCS_PbPb(mu0,mu1);
   fDimuCostEPnB[numdimu] = CostEPnB_PbPb(mu0,mu1,fPsi2Trkl);
   fDimuCostRPnB[numdimu] = CostEPnB_PbPb(mu0,mu1,fPsi2RP);

	 LabelOld1[numdimu]=i;
	 LabelOld2[numdimu]=j;
	 GoodMuon[i]=kTRUE;
	 GoodMuon[j]=kTRUE;
	 numdimu++;
	} // close loop on dimuon cuts
	delete dimu;
      } // close loop on second dimuon tracks
  } // close loop on first dimuon track

 // loop on single muons to keep only muons belonging to a dimuon surviving cuts
  for(Int_t i=0;i<ntracks;i++){
      AliAODTrack *mu0=(AliAODTrack*)fAODEvent->GetTrack(i);
      if (!mu0->IsMuonTrack()) continue;

     if(GoodMuon[i]) {
	fCharge[nummu] = mu0->Charge();
	fPt[nummu] = mu0->Pt();
	fPx[nummu] = mu0->Px();
	fPy[nummu] = mu0->Py();
	fPz[nummu] = mu0->Pz();
	fPt[nummu] = mu0->Pt();
	fY[nummu]  = mu0->Y();
	fEta[nummu]= mu0->Eta();
	fPhi[nummu]= mu0->Phi();
	fE[nummu] = mu0->E();
	fMatchTrig[nummu]   = mu0->GetMatchTrigger();
	fMatchTrigChi2[nummu]= mu0->GetChi2MatchTrigger();
	fRAtAbsEnd[nummu]=mu0->GetRAtAbsorberEnd();
     if(fMuonTrackCuts -> IsSelected(mu0)) fpDCA[nummu] = 1;
     nummu++;
     }
   }

 // reassign labels to muons belonging to a dimuon surviving cuts  [labels need to be reassigned since only selected tracks are saved and not all of them]
  for(Int_t i = 0;i<numdimu;i++){
    Int_t LabelNew1 = 0;
    Int_t LabelNew2 = 0;
    for(int j=0;j<LabelOld1[i];j++){
      if(GoodMuon[j]) LabelNew1++;
    }
    for(int j=0;j<LabelOld2[i];j++){
      if(GoodMuon[j]) LabelNew2++;
    }
    fDimuMu[i][0]=LabelNew1;
    fDimuMu[i][1]=LabelNew2;
   }

  fNMuons =nummu;
  fNDimu=numdimu;
  } // end loop on ntracks !=0
  //
  // keep only events where there is at least one dimuon
  if(fNDimu>0){
    fOutputTree->Fill();
    PostData(1,fOutputTree);
  }

  PostData(2,fhNEv);

}
//_____________________________________________________________________________
Float_t AliAnalysisTaskPbPbTree_MCut::GetVertex(AliAODEvent* aod) const
{

  Float_t vtxz = -999;

  const AliAODVertex* trkVtx = aod->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0)
    return vtxz;

    //add extra cuts

  vtxz = trkVtx->GetZ();

  return vtxz;
}
//__________________________________________________________
Short_t AliAnalysisTaskPbPbTree_MCut::GetVertexZ(Float_t vtxZ) const
{
    Short_t zvtx = -10;

    if (vtxZ >= -14. && vtxZ < -12.)
    zvtx = 0;
    if (vtxZ >= -12. && vtxZ < -10.)
    zvtx = 1;
    if (vtxZ >= -10. && vtxZ < -8.)
    zvtx = 2;
    if (vtxZ >= -8. && vtxZ < -6.)
    zvtx = 3;
    if (vtxZ >= -6. && vtxZ < -4.)
    zvtx = 4;
    if (vtxZ >= -4. && vtxZ < -3.)
    zvtx = 5;
    if (vtxZ >= -3. && vtxZ < -2.)
    zvtx = 6;
    if (vtxZ >= -2. && vtxZ < -1.)
    zvtx = 7;
    if (vtxZ >= -1. && vtxZ < 0)
    zvtx = 8;
    if (vtxZ >= 0 && vtxZ < 1.)
    zvtx = 9;
    if (vtxZ >= 1. && vtxZ < 2.)
    zvtx = 10;
    if (vtxZ >= 2. && vtxZ < 3.)
    zvtx = 11;
    if (vtxZ >= 3. && vtxZ < 4.)
    zvtx = 12;
    if (vtxZ >= 4. && vtxZ < 6.)
    zvtx = 13;
    if (vtxZ >= 6. && vtxZ < 8.)
    zvtx = 14;
    if (vtxZ >= 8. && vtxZ < 10.)
    zvtx = 15;
    if (vtxZ >= 10. && vtxZ < 12.)
    zvtx = 16;
    if (vtxZ >= 12. && vtxZ <= 14.)
    zvtx = 17;

    return zvtx;

}
//_____________________________________________________________________________
void AliAnalysisTaskPbPbTree_MCut::OpenInfoCalbration(Int_t run)
{
    if (!gGrid) {
        TGrid::Connect("alien://");
        printf("CONNECTING TO GRID ...\n");
    }

    // Searching the period name reading the corresponding run list
    std::string line;
    std::string strRun = std::to_string(run);
    std::string tmpPerName; // string for the period name

    printf("Opening calibV0TrklNoEtaCutRun215oVtx14MergedRuns.root\n");
    TFile* tmpFile215o = 0;
    tmpFile215o = TFile::Open("alien:///alice/cern.ch/user/l/lmichele/Event_Plane_calibration_files/calibV0TrklNoEtaCutRun215oVtx14MergedRuns.root");

    AliOADBContainer* cont215o = (AliOADBContainer*) tmpFile215o->Get("hMultV0BefCorPfpx");
    if((cont215o->GetObject(run))){tmpPerName = "215oVtx14MergedRuns";}
    else{printf("run %i does not belong to LHC15o\n", run);}
    tmpFile215o->Close();

    printf("Opening calibV0TrklNoEtaCutRun218qVtx14MRP2New.root\n");
    TFile* tmpFile218q = 0;
    tmpFile218q = TFile::Open("alien:///alice/cern.ch/user/l/lmichele/Event_Plane_calibration_files/calibV0TrklNoEtaCutRun218qVtx14MRP2New.root");
    AliOADBContainer* cont218q = (AliOADBContainer*) tmpFile218q->Get("hMultV0BefCorPfpx");
    if((cont218q->GetObject(run))){tmpPerName = "218qVtx14MRP2New";}
    else{printf("run %i does not belong to LHC18q\n", run);}
    tmpFile218q->Close();

    printf("Opening calibV0TrklNoEtaCutRun218rVtx14MRP2New.root\n");
    TFile* tmpFile218r = 0;
    tmpFile218r = TFile::Open("alien:///alice/cern.ch/user/l/lmichele/Event_Plane_calibration_files/calibV0TrklNoEtaCutRun218rVtx14MRP2New.root");
    AliOADBContainer* cont218r = (AliOADBContainer*) tmpFile218r->Get("hMultV0BefCorPfpx");
    if((cont218r->GetObject(run))){tmpPerName = "218rVtx14MRP2New";}
    else{printf("run %i does not belong to LHC18r\n", run);}
    tmpFile218r->Close();

    printf("RUN %i \n",run);
    printf("PERIOD NAME = %s\n",tmpPerName.c_str());

    TFile* foadb = 0;
    printf("Opening alien:///alice/cern.ch/user/l/lmichele/Event_Plane_calibration_files/calibV0TrklNoEtaCutRun%s.root",tmpPerName.c_str());
    foadb = TFile::Open(Form("alien:///alice/cern.ch/user/l/lmichele/Event_Plane_calibration_files/calibV0TrklNoEtaCutRun%s.root",tmpPerName.c_str()));
    printf("Reading calibration files...\n");

    AliOADBContainer* cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorPfpx");
    if(!cont){
        printf("OADB object hMultV0BefCorPfpx is not available in the file\n");
        return;
    }
    if(!(cont->GetObject(run))){
        printf("OADB object hMultV0BefCorPfpx is not available for run %i\n", run);
        return;
    }
    fMultV0 = ((TH1D*) cont->GetObject(run));

    for (Int_t k = 0; k < 18; k++){

        AliOADBContainer* contQx2am = (AliOADBContainer*) foadb->Get(Form("fqxa2m_%d", k));
        if(!contQx2am){
            printf("OADB object fqxa2m is not available in the file\n");
            return;
        }
        if(!(contQx2am->GetObject(run))){
            printf("OADB object fqxa2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0A[k]= ((TH1D*) contQx2am->GetObject(run));


        AliOADBContainer* contQy2am = (AliOADBContainer*) foadb->Get(Form("fqya2m_%d", k));
        if(!contQy2am){
            printf("OADB object fqya2m is not available in the file\n");
            return;
        }
        if(!(contQy2am->GetObject(run))){
            printf("OADB object fqya2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0A[k]= ((TH1D*) contQy2am->GetObject(run));



        AliOADBContainer* contQx2as = (AliOADBContainer*) foadb->Get(Form("fqxa2s_%d", k));
        if(!contQx2as){
            printf("OADB object fqxa2s is not available in the file\n");
            return;
        }
        if(!(contQx2as->GetObject(run))){
            printf("OADB object fqxa2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0A[k]= ((TH1D*) contQx2as->GetObject(run));


        AliOADBContainer* contQy2as = (AliOADBContainer*) foadb->Get(Form("fqya2s_%d", k));
        if(!contQy2as){
            printf("OADB object fqya2s is not available in the file\n");
            return;
        }
        if(!(contQy2as->GetObject(run))){
            printf("OADB object fqya2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0A[k]= ((TH1D*) contQy2as->GetObject(run));




        AliOADBContainer* contQx2cm = (AliOADBContainer*) foadb->Get(Form("fqxc2m_%d", k));
        if(!contQx2cm){
            printf("OADB object fqxc2m is not available in the file\n");
            return;
        }
        if(!(contQx2cm->GetObject(run))){
            printf("OADB object fqxc2m is not available for run %i\n", run);
            return;
        }
        fQx2mV0C[k]= ((TH1D*) contQx2cm->GetObject(run));


        AliOADBContainer* contQy2cm = (AliOADBContainer*) foadb->Get(Form("fqyc2m_%d", k));
        if(!contQy2cm){
            printf("OADB object fqyc2m is not available in the file\n");
            return;
        }
        if(!(contQy2cm->GetObject(run))){
            printf("OADB object fqyc2m is not available for run %i\n", run);
            return;
        }
        fQy2mV0C[k]= ((TH1D*) contQy2cm->GetObject(run));


        AliOADBContainer* contQx2cs = (AliOADBContainer*) foadb->Get(Form("fqxc2s_%d", k));
        if(!contQx2cs){
            printf("OADB object fqxc2s is not available in the file\n");
            return;
        }
        if(!(contQx2cs->GetObject(run))){
            printf("OADB object fqxc2s is not available for run %i\n", run);
            return;
        }
        fQx2sV0C[k]= ((TH1D*) contQx2cs->GetObject(run));


        AliOADBContainer* contQy2cs = (AliOADBContainer*) foadb->Get(Form("fqyc2s_%d", k));
        if(!contQy2cs){
            printf("OADB object fqyc2s is not available in the file\n");
            return;
        }
        if(!(contQy2cs->GetObject(run))){
            printf("OADB object fqyc2s is not available for run %i\n", run);
            return;
        }
        fQy2sV0C[k]= ((TH1D*) contQy2cs->GetObject(run));




        AliOADBContainer* contQx2trm = (AliOADBContainer*) foadb->Get(Form("fqxtr2m_%d", k));
        if(!contQx2trm){
            printf("OADB object fqxtr2m is not available in the file\n");
            return;
        }
        if(!(contQx2trm->GetObject(run))){
            printf("OADB object fqxtr2m is not available for run %i\n", run);
            return;
        }
        fQx2mTrk[k]= ((TH1D*) contQx2trm->GetObject(run));


        AliOADBContainer* contQy2trm = (AliOADBContainer*) foadb->Get(Form("fqytr2m_%d", k));
        if(!contQy2trm){
            printf("OADB object fqytr2m is not available in the file\n");
            return;
        }
        if(!(contQy2trm->GetObject(run))){
            printf("OADB object fqytr2m is not available for run %i\n", run);
            return;
        }
        fQy2mTrk[k]= ((TH1D*) contQy2trm->GetObject(run));


        AliOADBContainer* contQx2trs = (AliOADBContainer*) foadb->Get(Form("fqxtr2s_%d", k));
        if(!contQx2trs){
            printf("OADB object fqxtr2s is not available in the file\n");
            return;
        }
        if(!(contQx2trs->GetObject(run))){
            printf("OADB object fqxtr2s is not available for run %i\n", run);
            return;
        }
        fQx2sTrk[k]= ((TH1D*) contQx2trs->GetObject(run));


        AliOADBContainer* contQy2trs = (AliOADBContainer*) foadb->Get(Form("fqytr2s_%d", k));
        if(!contQy2trs){
            printf("OADB object fqytr2s is not available in the file\n");
            return;
        }
        if(!(contQy2trs->GetObject(run))){
            printf("OADB object fqytr2s is not available for run %i\n", run);
            return;
        }
        fQy2sTrk[k]= ((TH1D*) contQy2trs->GetObject(run));




        AliOADBContainer* contQx3am = (AliOADBContainer*) foadb->Get(Form("fqxa3m_%d", k));
        if(!contQx3am){
            printf("OADB object fqxa3m is not available in the file\n");
            return;
        }
        if(!(contQx3am->GetObject(run))){
            printf("OADB object fqxa3m is not available for run %i\n", run);
            return;
        }
        fQx3mV0A[k]= ((TH1D*) contQx3am->GetObject(run));


        AliOADBContainer* contQy3am = (AliOADBContainer*) foadb->Get(Form("fqya3m_%d", k));
        if(!contQy3am){
            printf("OADB object fqya3m is not available in the file\n");
            return;
        }
        if(!(contQy3am->GetObject(run))){
            printf("OADB object fqya3m is not available for run %i\n", run);
            return;
        }
        fQy3mV0A[k]= ((TH1D*) contQy3am->GetObject(run));


        AliOADBContainer* contQx3as = (AliOADBContainer*) foadb->Get(Form("fqxa3s_%d", k));
        if(!contQx3as){
            printf("OADB object fqxa3s is not available in the file\n");
            return;
        }
        if(!(contQx3as->GetObject(run))){
            printf("OADB object fqxa3s is not available for run %i\n", run);
            return;
        }
        fQx3sV0A[k]= ((TH1D*) contQx3as->GetObject(run));


        AliOADBContainer* contQy3as = (AliOADBContainer*) foadb->Get(Form("fqya3s_%d", k));
        if(!contQy3as){
            printf("OADB object fqya3s is not available in the file\n");
            return;
        }
        if(!(contQy3as->GetObject(run))){
            printf("OADB object fqya3s is not available for run %i\n", run);
            return;
        }
        fQy3sV0A[k]= ((TH1D*) contQy3as->GetObject(run));




        AliOADBContainer* contQx3cm = (AliOADBContainer*) foadb->Get(Form("fqxc3m_%d", k));
        if(!contQx3cm){
            printf("OADB object fqxc3m is not available in the file\n");
            return;
        }
        if(!(contQx3cm->GetObject(run))){
            printf("OADB object fqxc3m is not available for run %i\n", run);
            return;
        }
        fQx3mV0C[k]= ((TH1D*) contQx3cm->GetObject(run));


        AliOADBContainer* contQy3cm = (AliOADBContainer*) foadb->Get(Form("fqyc3m_%d", k));
        if(!contQy3cm){
            printf("OADB object fqyc3m is not available in the file\n");
            return;
        }
        if(!(contQy3cm->GetObject(run))){
            printf("OADB object fqyc3m is not available for run %i\n", run);
            return;
        }
        fQy3mV0C[k]= ((TH1D*) contQy3cm->GetObject(run));


        AliOADBContainer* contQx3cs = (AliOADBContainer*) foadb->Get(Form("fqxc3s_%d", k));
        if(!contQx3cs){
            printf("OADB object fqxc3s is not available in the file\n");
            return;
        }
        if(!(contQx3cs->GetObject(run))){
            printf("OADB object fqxc3s is not available for run %i\n", run);
            return;
        }
        fQx3sV0C[k]= ((TH1D*) contQx3cs->GetObject(run));


        AliOADBContainer* contQy3cs = (AliOADBContainer*) foadb->Get(Form("fqyc3s_%d", k));
        if(!contQy3cs){
            printf("OADB object fqyc3s is not available in the file\n");
            return;
        }
        if(!(contQy3cs->GetObject(run))){
            printf("OADB object fqyc3s is not available for run %i\n", run);
            return;
        }
        fQy3sV0C[k]= ((TH1D*) contQy3cs->GetObject(run));




        AliOADBContainer* contQx3trm = (AliOADBContainer*) foadb->Get(Form("fqxtr3m_%d", k));
        if(!contQx3trm){
            printf("OADB object fqxtr3m is not available in the file\n");
            return;
        }
        if(!(contQx3trm->GetObject(run))){
            printf("OADB object fqxtr3m is not available for run %i\n", run);
            return;
        }
        fQx3mTrk[k]= ((TH1D*) contQx3trm->GetObject(run));


        AliOADBContainer* contQy3trm = (AliOADBContainer*) foadb->Get(Form("fqytr3m_%d", k));
        if(!contQy3trm){
            printf("OADB object fqytr3m is not available in the file\n");
            return;
        }
        if(!(contQy3trm->GetObject(run))){
            printf("OADB object fqytr3m is not available for run %i\n", run);
            return;
        }
        fQy3mTrk[k]= ((TH1D*) contQy3trm->GetObject(run));


        AliOADBContainer* contQx3trs = (AliOADBContainer*) foadb->Get(Form("fqxtr3s_%d", k));
        if(!contQx3trs){
            printf("OADB object fqxtr3s is not available in the file\n");
            return;
        }
        if(!(contQx3trs->GetObject(run))){
            printf("OADB object fqxtr3s is not available for run %i\n", run);
            return;
        }
        fQx3sTrk[k]= ((TH1D*) contQx3trs->GetObject(run));


        AliOADBContainer* contQy3trs = (AliOADBContainer*) foadb->Get(Form("fqytr3s_%d", k));
        if(!contQy3trs){
            printf("OADB object fqytr3s is not available in the file\n");
            return;
        }
        if(!(contQy3trs->GetObject(run))){
            printf("OADB object fqytr3s is not available for run %i\n", run);
            return;
        }
        fQy3sTrk[k]= ((TH1D*) contQy3trs->GetObject(run));

    }
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskPbPbTree_MCut::CalcCorPhi(Double_t phi, Double_t dPhi) const
{

    phi += 39./34.*dPhi;
    if (phi < 0.) phi += 2.*TMath::Pi();
    if (phi > 2.*TMath::Pi()) phi -= 2.*TMath::Pi();
    return phi;

}
//______________________________________________________________________________
Double_t CostHE_PbPb(AliAODTrack* Mu0, AliAODTrack* Mu1){
  Double_t EBeam = 2510*208;
  Double_t mp = 195.323567174;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjLab(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargLab(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the LAB frame
  //
  TLorentzVector pMu1Lab(pla10,pla11,pla12,e1);
  TLorentzVector pMu2Lab(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the LAB frame
  //
  TLorentzVector pDimuLab = pMu1Lab + pMu2Lab;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuLab.E())*pDimuLab.Vect();
  TLorentzVector pMu1Dimu = pMu1Lab;
  TLorentzVector pMu2Dimu = pMu2Lab;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);

  TVector3 zaxis;
  zaxis=(pDimuLab.Vect()).Unit();
  //
  // --- Calculation of the polarization angle (Helicity)
  // (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  } else {
    cost = zaxis.Dot((pMu2Dimu.Vect()).Unit());
  }
  return cost;
}
//______________________________________________________________________________
//______________________________________________________________________________
Double_t PhiHE_PbPb(AliAODTrack* Mu0, AliAODTrack* Mu1){
  // Calculation the Helicity aimuthal angle (adapted from code by R. Arnaldi)
  Double_t EBeam = 2510*208;
  if(EBeam <= 0){
    printf("Can not compute phiHE with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 195.323567174;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we consider no crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the muon parameters in the dimuon rest frame
  //
  TVector3 zaxis = (pDimuCM.Vect()).Unit();
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);

  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);

  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxis),(pMu2Dimu.Vect()).Dot(xaxis));

   return phi;
}
//______________________________________________________________________________
Double_t CostCS_PbPb(AliAODTrack* Mu0, AliAODTrack* Mu1){
  Double_t EBeam = 2510*208;
  Double_t mp = 195.323567174;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0 -> Py();
  Double_t pla12 = Mu0 -> Pz();
  Double_t e1 = Mu0 -> E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  Double_t mu2Charge = Mu1 -> Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  Double_t cost;
  if(mu1Charge > 0) {
    cost = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
    // Theta CS is not properly defined for Like-Sign muons
    if(mu2Charge > 0 && cost<0) cost = -cost;
  } else {
    // Theta CS is not properly defined for Like-Sign muons
    cost = zaxisCS.Dot((pMu2Dimu.Vect()).Unit());
    if(mu2Charge < 0 && cost<0) cost = -cost;
  }
  return cost;
}
//______________________________________________________________________________
Double_t PhiCS_PbPb(AliAODTrack* Mu0, AliAODTrack* Mu1){
  // Cosinus of the Collins-Soper polar decay angle
  Double_t EBeam = 2510*208;
  if(EBeam <= 0){
    printf("Can not compute phiCS with EBeam=%f\n",EBeam);
    return -999999999;
  }
  Double_t mp = 195.323567174;
  Double_t pbeam = TMath::Sqrt(EBeam*EBeam - mp*mp);
  Double_t pla10 = Mu0 -> Px();
  Double_t pla11 = Mu0->Py();
  Double_t pla12 = Mu0->Pz();
  Double_t e1 = Mu0->E();
  Double_t mu1Charge = Mu0 -> Charge();
  Double_t pla20 = Mu1 -> Px();
  Double_t pla21 = Mu1 -> Py();
  Double_t pla22 = Mu1 -> Pz();
  Double_t e2 = Mu1 -> E();
  //Double_t mu2Charge=Mu1->Charge();
  if(pla10==0 && pla11==0 && pla12==0 && e1==0 && mu1Charge==0 && pla20==0 && pla21==0 && pla22==0 && e2==0.){return -666.;}

  // Fill the Lorentz vector for projectile and target
  // For the moment we do not consider the crossing angle
  // Projectile runs towards the MUON arm
  TLorentzVector pProjCM(0.,0.,-pbeam,EBeam); // projectile
  TLorentzVector pTargCM(0.,0., pbeam,EBeam); // target
  //
  // --- Get the muons parameters in the CM frame
  //
  TLorentzVector pMu1CM(pla10,pla11,pla12,e1);
  TLorentzVector pMu2CM(pla20,pla21,pla22,e2);
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;
  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  TLorentzVector pMu1Dimu = pMu1CM;
  TLorentzVector pMu2Dimu = pMu2CM;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the CS angle
  //
  TVector3 zaxisCS = (((pProjDimu.Vect()).Unit()) - ((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
   TVector3 yaxisCS = (((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
   TVector3 xaxisCS = (yaxisCS.Cross(zaxisCS)).Unit();

   Double_t phi;
   if(mu1Charge>0) phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
   else phi = TMath::ATan2((pMu2Dimu.Vect()).Dot(yaxisCS),((pMu2Dimu.Vect()).Dot(xaxisCS)));

   return phi;
}
//________________________________________________________________________
Double_t CostEPnB_PbPb(AliAODTrack* Mu0, AliAODTrack* Mu1, Double_t Psi){
  //printf("CosTheta in the Event-Plane reference frame --> No Boost of the Event-Plane vector in the J/psi rest frame \n");
  Double_t PxMu0      = Mu0 -> Px();
  Double_t PyMu0      = Mu0 -> Py();
  Double_t PzMu0      = Mu0 -> Pz();
  Double_t EMu0       = Mu0 -> E();
  Double_t ChargeMu0  = Mu0 -> Charge();
  Double_t PxMu1      = Mu1 -> Px();
  Double_t PyMu1      = Mu1 -> Py();
  Double_t PzMu1      = Mu1 -> Pz();
  Double_t EMu1       = Mu1 -> E();
  Double_t ChargeMu1  = Mu1 -> Charge();

  if(PxMu0 == 0 && PyMu0 == 0 && PzMu0 == 0 && EMu0 == 0 && ChargeMu0 == 0 && PxMu1 == 0 && PyMu1 == 0 && PzMu1 == 0 && EMu1 == 0.){return -666.;}

  //
  // --- Get the muons parameters in the laboratory frame
  //
  TLorentzVector pMu1CM(PxMu0,PyMu0,PzMu0,EMu0);
  TLorentzVector pMu2CM(PxMu1,PyMu1,PzMu1,EMu1);

  //
  // --- Calculation of the Q-vector and its orthogonal
  //
  TVector3 Qtr2vect(TMath::Cos(Psi),TMath::Sin(Psi),0.);
  //printf("Event Plane vector : \n");
  //Qtr2vect.Print();
  TVector3 Qtr2vectOrtLab = Qtr2vect.Orthogonal();
  //printf("Event Plane vector orthogonal : \n");
  //Qtr2vectOrtLab.Print();
  //
  // --- Obtain the dimuon parameters in the CM frame
  //
  TLorentzVector pDimuCM = pMu1CM + pMu2CM;

  //
  // --- Translate the dimuon parameters in the dimuon rest frame
  //
  TVector3 beta = (-1./pDimuCM.E())*pDimuCM.Vect();
  //printf("Beta : \n");
  //beta.Print();
  TLorentzVector pMu1Dimu = pMu1CM;
  //printf("4-vector mu1 : \n");
  //pMu1Dimu.Print();
  TLorentzVector pMu2Dimu = pMu2CM;
  //printf("4-vector mu2 : \n");
  //pMu2Dimu.Print();
  pMu1Dimu.Boost(beta);
  //printf("4-vector mu1 boosted : \n");
  //pMu1Dimu.Print();
  pMu2Dimu.Boost(beta);
  //printf("4-vector mu2 boosted : \n");
  //pMu2Dimu.Print();

  //
  // --- Determine the z axis for the EP angle
  //
  TVector3 zaxisEP = Qtr2vectOrtLab.Unit();
  //printf("zaxis Event Plane unitary \n");
  //zaxisEP.Print();

  //
  // --- Determine the EP angle (angle between mu+ and the z axis defined above)
  //
  Double_t CosTheta_EP;
  if(ChargeMu0 > 0) {
    CosTheta_EP = zaxisEP.Dot((pMu1Dimu.Vect()).Unit());
  } else {
    CosTheta_EP = zaxisEP.Dot((pMu2Dimu.Vect()).Unit());
  }
  //printf("CosTheta Event Plane = %f\n",CosTheta_EP);
  return CosTheta_EP;
}
//________________________________________________________________________
void AliAnalysisTaskPbPbTree_MCut::Terminate(Option_t *)
{

 }
