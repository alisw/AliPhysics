/*************************************************************************
* Copyright(c) 1998-2008,ALICE Experiment at CERN,All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use,copy,modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee,provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez  (cperez@cern.ch)
//          Pawel Debski  (pdebski@cern.ch)
//////////////////////////////////////////////////////

#include "TChain.h"
#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TFile.h"

#include "TRandom3.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliVVertex.h"
#include "AliVVZERO.h"
#include "AliStack.h"
#include "AliMCEvent.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"
#include "AliAODv0.h"
#include "AliAODTracklets.h"
#include "AliAODHeader.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include "TMath.h"
#include "TObjArray.h"
#include "AliFlowCandidateTrack.h"

#include "AliFlowTrackCuts.h"
#include "AliFlowEventCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowBayesianPID.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowVector.h"

#include "AliAnalysisTaskFlowStrange.h"

ClassImp(AliAnalysisTaskFlowStrange)

//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange() :
  AliAnalysisTaskSE(),
  fPIDResponse(NULL),
  fFB1(NULL),
  fFB1024(NULL),
  fTPCevent(NULL),
  fVZEevent(NULL),
  fCandidates(NULL),
  fList(NULL),
  fRunNumber(-1),
  fDebug(0),
  fQAlevel(0),
  fReadESD(kFALSE),
  fReadMC(kFALSE),
  fPostMatched(0),
  fAvoidExec(kFALSE),
  fSkipSelection(kFALSE),
  fSkipFlow(kFALSE),
  fSkipDHcorr(kTRUE),
  fUseFP(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fExtraEventRejection(kFALSE),
  fCentMethod("V0MTRK"),
  fCentPerMin(0),
  fCentPerMax(100),
  fThisCent(-1.0),
  fVertexZcut(10.0),
  fExcludeTPCEdges(kFALSE),
  fSpecie(0),
  fOnline(kFALSE),
  fHomemade(kFALSE),
  fWhichPsi(1),
  fVZEsave(kFALSE),
  fVZEload(NULL),
  fVZEResponse(NULL),
  fVZEmb(kFALSE),
  fVZEByDisk(kTRUE),
  fVZECa(0),
  fVZECb(3),
  fVZEAa(0),
  fVZEAb(3),
  fVZEQA(NULL),
  fPsi2(0.0),
  fMCEP(0.0),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fMinMassX(-1.0),
  fMaxMassX(-1.0),
  fRFPFilterBit(1),
  fRFPminPt(0.2),
  fRFPmaxPt(5.0),
  fRFPminEta(-0.8),
  fRFPmaxEta(+0.8),
  fRFPTPCsignal(10.0),
  fRFPmaxIPxy(2.4),
  fRFPmaxIPz(3.2),
  fRFPTPCncls(70),
  fDecayMass(0.0),
  fDecayPhi(0.0),
  fDecayEta(0.0),
  fDecayPt(0.0),
  fDecayDCAdaughters(0.0),
  fDecayCosinePointingAngleXY(0.0),
  fDecayRadXY(0.0),
  fDecayDecayLength(0.0),
  fDecayQt(0.0),
  fDecayAlpha(0.0),
  fDecayRapidity(0.0),
  fDecayProductIPXY(0.0),
  fDecayIPneg(0.0),
  fDecayIPpos(0.0),
  fDecayXneg(0.0),
  fDecayXpos(0.0),
  fDecayIDneg(-1),
  fDecayIDpos(-1),
  fDecayID(-1),
  fDecayMatchOrigin(0.0),
  fDecayMatchPhi(0.0),
  fDecayMatchEta(0.0),
  fDecayMatchPt(0.0),
  fDecayMatchRadXY(0.0),
  fDecayMinEta(0.0),
  fDecayMaxEta(0.0),
  fDecayMinPt(0.0),
  fDecayMaxDCAdaughters(0.0),
  fDecayMinCosinePointingAngleXY(0.0),
  fDecayMinQt(0.0),
  fDecayAPCutPie(kTRUE),
  fDecayMinRadXY(0.0),
  fDecayMaxDecayLength(0.0),
  fDecayMaxProductIPXY(0.0),
  fDecayMaxRapidity(0.0),
  fDaughterPhi(0.0),
  fDaughterEta(0.0),
  fDaughterPt(0.0),
  fDaughterNClsTPC(0),
  fDaughterCharge(0),
  fDaughterNFClsTPC(0),
  fDaughterNSClsTPC(0),
  fDaughterChi2PerNClsTPC(0.0),
  fDaughterXRows(0.0),
  fDaughterImpactParameterXY(0.0),
  fDaughterImpactParameterZ(0.0),
  fDaughterStatus(0),
  fDaughterITScm(0),
  fDaughterNSigmaPID(0.0),
  fDaughterKinkIndex(0),
  fDaughterMatchPhi(0.0),
  fDaughterMatchEta(0.0),
  fDaughterMatchPt(0.0),
  fDaughterMatchImpactParameterXY(0.0),
  fDaughterMatchImpactParameterZ(0.0),
  fDaughterUnTag(kTRUE),
  fDaughterMinEta(0.0),
  fDaughterMaxEta(0.0),
  fDaughterMinPt(0.0),
  fDaughterMinNClsTPC(0),
  fDaughterMinXRows(0),
  fDaughterMaxChi2PerNClsTPC(0.0),
  fDaughterMinXRowsOverNClsFTPC(0.0),
  fDaughterMinImpactParameterXY(0.0),
  fDaughterMaxNSigmaPID(0.0) {
  //ctor
  for(Int_t i=0; i!=6; ++i) fDaughterITSConfig[i]=-1;
}
//=======================================================================
AliAnalysisTaskFlowStrange::AliAnalysisTaskFlowStrange(const char *name) :
  AliAnalysisTaskSE(name),
  fPIDResponse(NULL),
  fFB1(NULL),
  fFB1024(NULL),
  fTPCevent(NULL),
  fVZEevent(NULL),
  fCandidates(NULL),
  fList(NULL),
  fRunNumber(-1),
  fDebug(0),
  fQAlevel(0),
  fReadESD(kFALSE),
  fReadMC(kFALSE),
  fPostMatched(0),
  fAvoidExec(kFALSE),
  fSkipSelection(kFALSE),
  fSkipFlow(kFALSE),
  fSkipDHcorr(kTRUE),
  fUseFP(kFALSE),
  fRunOnpA(kFALSE),
  fRunOnpp(kFALSE),
  fExtraEventRejection(kFALSE),
  fCentMethod("V0MTRK"),
  fCentPerMin(0),
  fCentPerMax(100),
  fThisCent(-1.0),
  fVertexZcut(10.0),
  fExcludeTPCEdges(kFALSE),
  fSpecie(0),
  fOnline(kFALSE),
  fHomemade(kFALSE),
  fWhichPsi(1),
  fVZEsave(kFALSE),
  fVZEload(NULL),
  fVZEResponse(NULL),
  fVZEmb(kFALSE),
  fVZEByDisk(kTRUE),
  fVZECa(0),
  fVZECb(3),
  fVZEAa(0),
  fVZEAb(3),
  fVZEQA(NULL),
  fPsi2(0.0),
  fMCEP(0.0),
  fMassBins(0),
  fMinMass(0.0),
  fMaxMass(0.0),
  fMinMassX(-1.0),
  fMaxMassX(-1.0),
  fRFPFilterBit(1),
  fRFPminPt(0.2),
  fRFPmaxPt(5.0),
  fRFPminEta(-0.8),
  fRFPmaxEta(+0.8),
  fRFPTPCsignal(10.0),
  fRFPmaxIPxy(2.4),
  fRFPmaxIPz(3.2),
  fRFPTPCncls(70),
  fDecayMass(0.0),
  fDecayPhi(0.0),
  fDecayEta(0.0),
  fDecayPt(0.0),
  fDecayDCAdaughters(0.0),
  fDecayCosinePointingAngleXY(0.0),
  fDecayRadXY(0.0),
  fDecayDecayLength(0.0),
  fDecayQt(0.0),
  fDecayAlpha(0.0),
  fDecayRapidity(0.0),
  fDecayProductIPXY(0.0),
  fDecayIPneg(0.0),
  fDecayIPpos(0.0),
  fDecayXneg(0.0),
  fDecayXpos(0.0),
  fDecayIDneg(-1),
  fDecayIDpos(-1),
  fDecayID(-1),
  fDecayMatchOrigin(0.0),
  fDecayMatchPhi(0.0),
  fDecayMatchEta(0.0),
  fDecayMatchPt(0.0),
  fDecayMatchRadXY(0.0),
  fDecayMinEta(0.0),
  fDecayMaxEta(0.0),
  fDecayMinPt(0.0),
  fDecayMaxDCAdaughters(0.0),
  fDecayMinCosinePointingAngleXY(0.0),
  fDecayMinQt(0.0),
  fDecayAPCutPie(kTRUE),
  fDecayMinRadXY(0.0),
  fDecayMaxDecayLength(0.0),
  fDecayMaxProductIPXY(0.0),
  fDecayMaxRapidity(0.0),
  fDaughterPhi(0.0),
  fDaughterEta(0.0),
  fDaughterPt(0.0),
  fDaughterNClsTPC(0),
  fDaughterCharge(0),
  fDaughterNFClsTPC(0),
  fDaughterNSClsTPC(0),
  fDaughterChi2PerNClsTPC(0.0),
  fDaughterXRows(0.0),
  fDaughterImpactParameterXY(0.0),
  fDaughterImpactParameterZ(0.0),
  fDaughterStatus(0),
  fDaughterITScm(0),
  fDaughterNSigmaPID(0.0),
  fDaughterKinkIndex(0),
  fDaughterMatchPhi(0.0),
  fDaughterMatchEta(0.0),
  fDaughterMatchPt(0.0),
  fDaughterMatchImpactParameterXY(0.0),
  fDaughterMatchImpactParameterZ(0.0),
  fDaughterUnTag(kTRUE),
  fDaughterMinEta(0.0),
  fDaughterMaxEta(0.0),
  fDaughterMinPt(0.0),
  fDaughterMinNClsTPC(0),
  fDaughterMinXRows(0),
  fDaughterMaxChi2PerNClsTPC(0.0),
  fDaughterMinXRowsOverNClsFTPC(0.0),
  fDaughterMinImpactParameterXY(0.0),
  fDaughterMaxNSigmaPID(0.0) {
  //ctor
  for(Int_t i=0; i!=6; ++i) fDaughterITSConfig[i]=-1;
  DefineInput( 0,TChain::Class());
  DefineOutput(1,TList::Class());
  DefineOutput(2,AliFlowEventSimple::Class()); // TPC object
  DefineOutput(3,AliFlowEventSimple::Class()); // VZE object
}
//=======================================================================
AliAnalysisTaskFlowStrange::~AliAnalysisTaskFlowStrange() {
  //dtor
  if (fCandidates) delete fCandidates;
  if (fTPCevent)   delete fTPCevent;
  if (fVZEevent)   delete fVZEevent;
  if (fList)       delete fList;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserCreateOutputObjects() {
  //UserCreateOutputObjects
  fList=new TList();
  fList->SetOwner();
  AddQAEvents();
  AddQACandidates();

  if(fReadESD) MakeFilterBits();

  AliFlowCommonConstants *cc = AliFlowCommonConstants::GetMaster();
  cc->SetNbinsMult(3000); cc->SetMultMin(0);   cc->SetMultMax(30000);
  cc->SetNbinsPt(200); cc->SetPtMin(0.0);   cc->SetPtMax(20.0);
  cc->SetNbinsPhi(100);  cc->SetPhiMin(0.0);  cc->SetPhiMax(TMath::TwoPi());
  cc->SetNbinsEta(100);  cc->SetEtaMin(-5.0); cc->SetEtaMax(+5.0);
  cc->SetNbinsQ(100);    cc->SetQMin(0.0);    cc->SetQMax(3.0);
  cc->SetNbinsMass(fMassBins);
  cc->SetMassMin(fMinMass);
  cc->SetMassMax(fMaxMass);

  if(fMinMassX<0) {
    if(fSpecie==0) {
      fMinMassX = 0.494;
      fMaxMassX = 0.502;
      //      double lowEdge[14]={0.398, 0.420, 0.444, 0.468, 0.486,
      //			  0.490, 0.494, 0.498, 0.502, 0.506, 
      //			  0.524, 0.548, 0.572, 0.598};
    } else {
      fMinMassX = 1.114;
      fMaxMassX = 1.118;
      //      double lowEdge[13]={1.084, 1.094, 1.104, 1.110, 1.114,
      //			  1.116, 1.118, 1.122, 1.128, 1.138,
      //			  1.148, 1.158, 1.168};} else {
    }
  }

  //loading pid response
  AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  fTPCevent = new AliFlowEvent(100);
  fVZEevent = new AliFlowEvent(100);

  //array of candidates
  fCandidates = new TObjArray(100);
  fCandidates->SetOwner();

  PostData(1,fList);
  if(fUseFP) { // for connection to the flow package
    PostData(2,fTPCevent);
    PostData(3,fVZEevent);
  }

  gRandom->SetSeed();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQAEvents() {
  // function to add event qa
  TH1D *tH1D;
  TProfile *tProfile;
  TList *tQAEvents=new TList();
  tQAEvents->SetName("Event");
  tQAEvents->SetOwner();
  tH1D = new TH1D("Events","Number of Events",6,0,6); tQAEvents->Add(tH1D);
  tH1D->GetXaxis()->SetBinLabel(1,"exec");
  tH1D->GetXaxis()->SetBinLabel(2,"userexec");
  tH1D->GetXaxis()->SetBinLabel(3,"reached");
  tH1D->GetXaxis()->SetBinLabel(4,"selected");
  tH1D->GetXaxis()->SetBinLabel(5,"rejectedByLowQw");
  tH1D->GetXaxis()->SetBinLabel(6,"rejectedByErrorLoadVZEcal");
  tProfile = new TProfile("Configuration","Configuration",20,0,20); tQAEvents->Add(tProfile);
  tProfile->Fill( 0.5,fCentPerMin,1); tProfile->GetXaxis()->SetBinLabel( 1,"fCentPerMin");
  tProfile->Fill( 1.5,fCentPerMax,1); tProfile->GetXaxis()->SetBinLabel( 2,"fCentPerMax");
  tProfile->Fill( 2.5,fDaughterMinEta,1);               tProfile->GetXaxis()->SetBinLabel( 3,"fDaughterMinEta");
  tProfile->Fill( 3.5,fDaughterMaxEta,1);               tProfile->GetXaxis()->SetBinLabel( 4,"fDaughterMaxEta");
  tProfile->Fill( 4.5,fDaughterMinPt,1);                tProfile->GetXaxis()->SetBinLabel( 5,"fDaughterMinPt");
  tProfile->Fill( 5.5,fDaughterMinNClsTPC,1);           tProfile->GetXaxis()->SetBinLabel( 6,"fDaughterMinNClsTPC");
  tProfile->Fill( 6.5,fDaughterMaxChi2PerNClsTPC,1);    tProfile->GetXaxis()->SetBinLabel( 7,"fDaughterMaxChi2PerNClsTPC");
  tProfile->Fill( 7.5,fDaughterMinXRowsOverNClsFTPC,1); tProfile->GetXaxis()->SetBinLabel( 8,"fDaughterMinXRowsOverNClsFTPC");
  tProfile->Fill( 8.5,fDaughterMinImpactParameterXY,1); tProfile->GetXaxis()->SetBinLabel( 9,"fDaughterMinImpactParameterXY");
  tProfile->Fill( 9.5,fDaughterMaxNSigmaPID,1);         tProfile->GetXaxis()->SetBinLabel(10,"fDaughterMaxNSigmaPID");
  tProfile->Fill(10.5,fDecayMaxDCAdaughters,1);          tProfile->GetXaxis()->SetBinLabel(11,"fDecayMaxDCAdaughters");
  tProfile->Fill(11.5,fDecayMinCosinePointingAngleXY,1); tProfile->GetXaxis()->SetBinLabel(12,"fDecayMinCosinePointingAngleXY");
  tProfile->Fill(12.5,fDecayMinQt,1);                    tProfile->GetXaxis()->SetBinLabel(13,"fDecayMinQt");
  tProfile->Fill(13.5,fDecayMinRadXY,1);                 tProfile->GetXaxis()->SetBinLabel(14,"fDecayMinRadXY");
  tProfile->Fill(14.5,fDecayMaxDecayLength,1);           tProfile->GetXaxis()->SetBinLabel(15,"fDecayMaxDecayLength");
  tProfile->Fill(15.5,fDecayMaxProductIPXY,1);           tProfile->GetXaxis()->SetBinLabel(16,"fDecayMaxProductIPXY");
  tProfile->Fill(16.5,fDecayMaxRapidity,1);              tProfile->GetXaxis()->SetBinLabel(17,"fDecayMaxRapidity");
  tProfile->Fill(17.5,fDecayMinEta,1);                   tProfile->GetXaxis()->SetBinLabel(18,"fDecayMinEta");
  tProfile->Fill(18.5,fDecayMaxEta,1);                   tProfile->GetXaxis()->SetBinLabel(19,"fDecayMaxEta");
  tProfile->Fill(19.5,fDecayMinPt,1);                    tProfile->GetXaxis()->SetBinLabel(20,"fDecayMinPt");

  tH1D = new TH1D("POI","POIs;multiplicity",800,0,800);         tQAEvents->Add(tH1D);
  tH1D = new TH1D("UNTAG","UNTAG;Untagged Daughters",800,0,800);tQAEvents->Add(tH1D);
  tH1D = new TH1D("RealTime","RealTime;LogT sec",2000,-10,+10); tQAEvents->Add(tH1D);
  fList->Add(tQAEvents);
  AddEventSpy("EventsReached");
  AddEventSpy("EventsSelected");
  AddMakeQSpy();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddEventSpy(TString name) {
  TH1D *tH1D;
  TH2D *tH2D;
  TList *tList=new TList();
  tList->SetName(name.Data());
  tList->SetOwner();
  tH2D = new TH2D("VTXZ","VTXZ;Global||SPD;SPD",60,-25,+25,60,-25,+25); tList->Add( tH2D );
  tH2D = new TH2D("CCCC","CCCC;V0M;TRK",60,-10,110,60,-10,110);         tList->Add( tH2D );
  tH2D = new TH2D("REFM","REFM;TPC;GLOBAL",100,0,3000,100,0,3000);      tList->Add( tH2D );
  if(fReadMC) {
    tH1D = new TH1D("MCCC","MCCC;Xsection",100,-10,110); tList->Add( tH1D );
    tH1D = new TH1D("MCEP","MCEP;MCEP",100,-TMath::TwoPi(),TMath::TwoPi()); tList->Add( tH1D );
  }
  fList->Add(tList);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddMakeQSpy() {
  if(fSkipFlow) return;
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPF1;
  TList *tList=new TList();
  tList->SetName("MakeQSpy");
  tList->SetOwner();
  tH1D = new TH1D("RFPTPC","TPC Refrence Multiplicity;multiplicity",3000,0,3000);     tList->Add( tH1D );
  tH1D = new TH1D("RFPVZE","VZERO Reference Multiplicity;multiplicity",3000,0,30000); tList->Add( tH1D );
  tH1D = new TH1D("QmTPC","TPC Normalized Q vector;|Q|/M",3000,0,1);   tList->Add( tH1D );
  tH1D = new TH1D("QmVZE","VZERO Normalized Q vector;|Q|/M",3000,0,1); tList->Add( tH1D );
  tH2D = new TH2D("TPCAllPhiEta","TPCall;Phi;Eta",180,0,TMath::TwoPi(),80,-0.9,+0.9); tList->Add( tH2D );
  tH2D = new TH2D("VZEAllPhiEta","VZEall;Phi;Eta",20,0,TMath::TwoPi(),40,-4.0,+6.0);  tList->Add( tH2D );
  tH1D = new TH1D("TPCPSI","TPCPSI;PSI",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("TPCPSIA","TPCPSIA;PSIA",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("TPCPSIB","TPCPSIB;PSIB",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSI","VZEPSI;PSI",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSIA","VZEPSIA;PSIA",72,0,TMath::Pi()); tList->Add( tH1D );
  tH1D = new TH1D("VZEPSIB","VZEPSIB;PSIB",72,0,TMath::Pi()); tList->Add( tH1D );
  tPF1 = new TProfile("TPCQ","TPCQ",6,0.5,6.5,"s");
  tPF1->GetXaxis()->SetBinLabel(1,"Qay"); tPF1->GetXaxis()->SetBinLabel(2,"Qax");
  tPF1->GetXaxis()->SetBinLabel(3,"Qby"); tPF1->GetXaxis()->SetBinLabel(4,"Qbx");
  tPF1->GetXaxis()->SetBinLabel(5,"Qy");  tPF1->GetXaxis()->SetBinLabel(6,"Qx");
  tList->Add( tPF1 );
  tPF1 = new TProfile("VZEQ","VZEQ",6,0.5,6.5,"s");
  tPF1->GetXaxis()->SetBinLabel(1,"Qay"); tPF1->GetXaxis()->SetBinLabel(2,"Qax");
  tPF1->GetXaxis()->SetBinLabel(3,"Qby"); tPF1->GetXaxis()->SetBinLabel(4,"Qbx");
  tPF1->GetXaxis()->SetBinLabel(5,"Qy");  tPF1->GetXaxis()->SetBinLabel(6,"Qx");
  tList->Add( tPF1 );
  
  fList->Add(tList);
  if(!fSkipFlow) {
    tList=new TList(); tList->SetName("TPCRFPall"); tList->SetOwner(); AddTPCRFPSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("TPCRFPsel"); tList->SetOwner(); AddTPCRFPSpy(tList); fList->Add(tList);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddQACandidates() {
  // function to add histogramming for candidates
  if(fSkipSelection) return;

  TList *tList;
  TH1D *tH1D;
  TH2D *tH2D;
  TH3D *tH3D;

  //reconstruction
  if(fReadESD) {
    tList=new TList(); tList->SetName("TrkAll"); tList->SetOwner(); AddTracksSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkSel"); tList->SetOwner(); AddTracksSpy(tList); fList->Add(tList);
    tH2D = new TH2D("NPAIR", "NPAIR;NPOS;NNEG",1000,0,5000,1000,0,5000); tList->Add(tH2D);
    tH2D = new TH2D("PtIPXY","PtIPXY;Pt;IPxy", 100,0,10,200,-10,+10); tList->Add(tH2D);
  }
  //candidates
  tList=new TList(); tList->SetName("RecAll"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  tH2D = new TH2D("V0SADC","V0S AFTER DAUGHTER CUTS;V0ALL;V0IMW",100,0,1000,100,0,1000); tList->Add(tH2D);
  tList=new TList(); tList->SetName("RecSel"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  //daughters
  tList=new TList(); tList->SetName("TrkDau"); tList->SetOwner(); AddTracksSpy(tList); fList->Add(tList);
  if(!fSkipDHcorr) {
    //corr
    tList=new TList(); tList->SetName("DHCORR"); tList->SetOwner(); 
    tH3D = new TH3D("DPHI","DPHI;dPT;dPHI;dETA", 20, -1, +1, 120, -TMath::TwoPi(), TMath::TwoPi(), 16, -1.6, +1.6 ); tList->Add(tH3D);
    fList->Add(tList);
  }
  if(fQAlevel>1) {
    // IN-OUT
    tList=new TList(); tList->SetName("RecAllIP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("RecAllOP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("RecSelIP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("RecSelOP"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
  }
  //match
  if(fReadMC) {
    tList=new TList(); tList->SetName("STATMC"); tList->SetOwner(); fList->Add(tList);
    tH1D = new TH1D("Events", "Events",5,0.5,5.5); tList->Add(tH1D);
    tH1D->GetXaxis()->SetBinLabel(1,"Selected events");
    tH1D->GetXaxis()->SetBinLabel(2,"Stack found");
    tH1D->GetXaxis()->SetBinLabel(3,"Daughters in stack");
    tH1D->GetXaxis()->SetBinLabel(4,"Correspond to decay");
    tH1D->GetXaxis()->SetBinLabel(5,"Decay has mother");

    tList=new TList(); tList->SetName("RecMth"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkMth"); tList->SetOwner(); AddTracksSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("RecUNMth"); tList->SetOwner(); AddCandidatesSpy(tList,false); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkUNMth"); tList->SetOwner(); AddTracksSpy(tList,false); fList->Add(tList);

    tList=new TList(); tList->SetName("NegNegMth"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("NegPosMth"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkNegMth"); tList->SetOwner(); AddTracksSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("TrkPosMth"); tList->SetOwner(); AddTracksSpy(tList,true); fList->Add(tList);
    tList=new TList(); tList->SetName("MthFDW"); tList->SetOwner(); AddCandidatesSpy(tList,true); fList->Add(tList);

  }
  //stack
  if(fReadMC) {
    tList=new TList(); tList->SetName("GenTru"); tList->SetOwner(); AddCandidatesSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTK0sGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTLdaGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTPhiGenAcc"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTXiGenAcc");  tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTK0s"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
    tList=new TList(); tList->SetName("MCTLda"); tList->SetOwner(); AddMCParticleSpy(tList); fList->Add(tList);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Exec(Option_t* option) {
  // bypassing ::exec (needed because of AMPT)
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(0);
  if(fAvoidExec) {
    AliAnalysisTaskFlowStrange::UserExec(option);
  } else {
    AliAnalysisTaskSE::Exec(option);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::UserExec(Option_t *option) {
  // bridge
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(1);
  AliAnalysisTaskFlowStrange::MyUserExec(option);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyNotifyRun() {
  if(fQAlevel>5 && !fReadESD) AddVZEQA();
  if(fVZEsave) AddVZEROResponse();
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::CalibrateEvent() {
  if(fVZEsave) SaveVZEROResponse();
  if(fQAlevel>5 && !fReadESD) SaveVZEROQA(); // 2BIMPROVED
  Bool_t okay=kTRUE;
  if(fVZEload) {
    LoadVZEROResponse();
    if(!fVZEResponse) okay = kFALSE;
  }
  return okay;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptAAEvent(AliESDEvent *tESD) {
  Double_t acceptEvent=kTRUE;
  Double_t tTPCVtxZ = tESD->GetPrimaryVertexTPC()->GetZ();
  if(tESD->GetPrimaryVertexTPC()->GetNContributors()<=0) return kFALSE;
  Double_t tSPDVtxZ = tESD->GetPrimaryVertexSPD()->GetZ();
  if(tESD->GetPrimaryVertexSPD()->GetNContributors()<=0) return kFALSE;
  // EventCuts
  AliCentrality *cent = tESD->GetCentrality();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  TString mycent = fCentMethod;
  if(fCentMethod.Contains("V0MTRK")) {
    acceptEvent = TMath::Abs(cc1-cc2)>5.0?kFALSE:acceptEvent; // a la Alex
    mycent = "V0M";
  }
  fThisCent = cent->GetCentralityPercentile( mycent );
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ-tSPDVtxZ)>0.5?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  // EndOfCuts
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  }
  return acceptEvent;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptAAEvent(AliAODEvent *tAOD) {
  Double_t acceptEvent=kTRUE;
  //=>Pile-up rejection (hardcoded)
  Double_t tVtxZ = tAOD->GetPrimaryVertex()->GetZ();
  if(tAOD->GetPrimaryVertex()->GetNContributors()<=0) return kFALSE;
  Double_t tSPDVtxZ = tAOD->GetPrimaryVertexSPD()->GetZ();
  if(tAOD->GetPrimaryVertexSPD()->GetNContributors()<=0) return kFALSE;
  Int_t tpc = RefMultTPC();
  Int_t glo = RefMultGlobal();
  if(fExtraEventRejection) {
    TString name = tAOD->GetPrimaryVertex()->GetTitle();
    if( !name.Contains("VertexerTracks") ) return kFALSE;
    if( ( Float_t(tpc) < -40.3+1.22*glo ) || ( Float_t(tpc)>(32.1+1.59*glo) ) ) return kFALSE;
  }
  // EventCuts
  AliCentrality *cent = tAOD->GetHeader()->GetCentralityP();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  TString mycent = fCentMethod;
  if(fCentMethod.Contains("V0MTRK")) {
    acceptEvent = TMath::Abs(cc1-cc2)>5.0?kFALSE:acceptEvent;
    mycent = "V0M";
  }
  fThisCent = cent->GetCentralityPercentile( mycent );

  Double_t xsec=0;
  if(fReadMC) {
    AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(tAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if (!mcHeader) {
      return kFALSE;
    }
    xsec = mcHeader->GetCrossSection();
    //fMCEP = mcHeader->GetReactionPlaneAngle();
    fMCEP = mcHeader->GetReactionPlaneAngle() + (gRandom->Rndm()>0.5)*TMath::Pi();
  }

  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ-tSPDVtxZ)>0.5?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("REFM"))->Fill( tpc, glo );
  if(fReadMC) {
    ((TH1D*)((TList*)fList->FindObject("EventsReached"))->FindObject("MCCC"))->Fill( xsec );
    ((TH1D*)((TList*)fList->FindObject("EventsReached"))->FindObject("MCEP"))->Fill( fMCEP );
  }
  // EndOfCuts
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("REFM"))->Fill( tpc, glo );
    if(fReadMC) {
      ((TH1D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("MCCC"))->Fill( xsec );
      ((TH1D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("MCEP"))->Fill( fMCEP );
    }
  }
  return acceptEvent;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptPPEvent(AliAODEvent *tAOD) {
  Double_t acceptEvent=kTRUE;
  Double_t tVtxZ = tAOD->GetPrimaryVertex()->GetZ();
  if(tAOD->GetPrimaryVertex()->GetNContributors()<=0) return kFALSE;
  Double_t tSPDVtxZ = tAOD->GetPrimaryVertexSPD()->GetZ();
  if(tAOD->GetPrimaryVertexSPD()->GetNContributors()<=0) return kFALSE;
  // EventCuts
  AliCentrality *cent = tAOD->GetHeader()->GetCentralityP();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  fThisCent = GetReferenceMultiplicity();
  //for pp i use fCentPerXXX to select on multiplicity
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ-tSPDVtxZ)>0.5?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  // EndOfCuts
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  }
  return acceptEvent;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::GetReferenceMultiplicity() { //toberefined
  AliAODEvent *tAOD = (AliAODEvent *) InputEvent();
  if(!tAOD) return -1;
  AliAODTrack *track;
  Int_t rawN = tAOD->GetNumberOfTracks();
  Int_t ref=0;
  for(Int_t id=0; id!=rawN; ++id) {
    track = tAOD->GetTrack(id);
    if(!track->TestFilterBit(fRFPFilterBit)) continue;
    ++ref;
  }
  return ref;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptPAEvent(AliAODEvent *tAOD) {
  //if(aod->GetHeader()->GetEventNumberESDFile() == 0) return; //rejecting first chunk NOT NEEDED ANYMORE
  Int_t bc2 = tAOD->GetHeader()->GetIRInt2ClosestInteractionMap();
  if(bc2!=0) return kFALSE;
  Int_t bc1 = tAOD->GetHeader()->GetIRInt1ClosestInteractionMap();
  if(bc1!=0) return kFALSE;
  Short_t isPileup = tAOD->IsPileupFromSPD(5);
  if(isPileup!=0) return kFALSE;
  if(tAOD->GetHeader()->GetRefMultiplicityComb08()<0) return kFALSE;

  const AliAODVertex* spdVtx = tAOD->GetPrimaryVertexSPD();
  if(!spdVtx) return kFALSE;
  if(spdVtx->GetNContributors()<=0) return kFALSE;

  const AliAODVertex* tpcVtx=NULL;
  Int_t nVertices = tAOD->GetNumberOfVertices();
  for(Int_t iVertices = 0; iVertices < nVertices; iVertices++){
    const AliAODVertex* vertex = tAOD->GetVertex(iVertices);
    if (vertex->GetType() != AliAODVertex::kMainTPC) continue;
    tpcVtx = vertex;
  }
  if(!tpcVtx) return kFALSE;
  if(tpcVtx->GetNContributors()<=0) return kFALSE;
  Double_t tTPCVtxZ = tpcVtx->GetZ();
  Double_t tSPDVtxZ = spdVtx->GetZ();
  if (TMath::Abs(tSPDVtxZ - tTPCVtxZ)>2.0) return kFALSE;
  if(plpMV(tAOD)) return kFALSE;

  Double_t acceptEvent=kTRUE;
  // EventCuts
  AliCentrality *cent = tAOD->GetHeader()->GetCentralityP();
  Double_t cc1, cc2;
  cc1 = cent->GetCentralityPercentile("V0M");
  cc2 = cent->GetCentralityPercentile("TRK");
  if(fCentMethod.Contains("V0MTRK")) fCentMethod = "V0M";
  fThisCent = cent->GetCentralityPercentile( fCentMethod );
  acceptEvent = (fThisCent<fCentPerMin||fThisCent>fCentPerMax)?kFALSE:acceptEvent;
  acceptEvent = TMath::Abs(tTPCVtxZ)>fVertexZcut?kFALSE:acceptEvent;
  // EndOfCuts
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
  ((TH2D*)((TList*)fList->FindObject("EventsReached"))->FindObject("CCCC"))->Fill( cc1, cc2 );
  if(acceptEvent) {
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("VTXZ"))->Fill( tTPCVtxZ, tSPDVtxZ );
    ((TH2D*)((TList*)fList->FindObject("EventsSelected"))->FindObject("CCCC"))->Fill( cc1, cc2 );    
  }
  return acceptEvent;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MyUserExec(Option_t *) {
  // MAIN ROUTINE
  TStopwatch tTime;
  tTime.Start();
  fCandidates->SetLast(-1);
  AliESDEvent *tESD=dynamic_cast<AliESDEvent*>(InputEvent());
  AliAODEvent *tAOD=dynamic_cast<AliAODEvent*>(InputEvent());
  Int_t thisRun = fRunNumber;
  //=>check event
  Bool_t acceptEvent=kFALSE;
  if(fReadESD) {
    if(!tESD) {ResetContainers(); Publish(); return;}
    acceptEvent = fRunOnpp?kFALSE:fRunOnpA?kFALSE:AcceptAAEvent(tESD);
    thisRun = tESD->GetRunNumber();
  } else {
    if(!tAOD) {ResetContainers(); Publish(); return;}
    acceptEvent = fRunOnpp?AcceptPPEvent(tAOD):fRunOnpA?AcceptPAEvent(tAOD):AcceptAAEvent(tAOD);
    thisRun = tAOD->GetRunNumber();
  }
  if(thisRun!=fRunNumber) {
    fRunNumber = thisRun;
    MyNotifyRun();
  }
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(2);
  //=>does the event clear?
  if(!acceptEvent) {ResetContainers(); Publish(); return;}
  // healthy event incomming
  if( !CalibrateEvent() ) { // saves/retrieves/qas VZEROCAL
    ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(5);
    ResetContainers(); Publish(); return; // issue retrieving callibration
  }
  if(!fSkipFlow) {
    MakeQVectors();
    if(fPsi2<-0.1) {
      ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(4);
      ResetContainers(); Publish(); return;
    }
  }
  //=>great, lets do our stuff!
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("Events"))->Fill(3);
  //=>load candidates
  if(!fSkipSelection) {
    if(fReadESD) {
      ReadFromESD(tESD);
    } else {
      if(fSpecie<10) ReadFromAODv0(tAOD);
      else ChargeParticles(tAOD);
    }
    if(fUseFP) {
      if(!fSkipDHcorr) MakeDHcorr();
      AddCandidates();
    }
    //=>flow
    //=>done
  }
  tTime.Stop();
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("RealTime"))->Fill( TMath::Log( tTime.RealTime() ) );
  Publish();
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Publish() {
  PostData(1,fList);
  if(fUseFP) {
    PostData(2,fTPCevent);
    PostData(3,fVZEevent);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromESD(AliESDEvent *tESD) {
  AliStack *stack=NULL;
  if(fReadMC) {
    AliMCEvent *mcevent=NULL;
    mcevent = MCEvent();
    if(mcevent) stack = mcevent->Stack();
  }

  Int_t num = tESD->GetNumberOfTracks();
  AliESDtrack *myTrack;
  Int_t plist[3000], nlist[3000], np=0, nn=0;
  Double_t pd0[3000], nd0[3000];
  for (Int_t i=0; i!=num; ++i) {
    myTrack = (AliESDtrack*) tESD->GetTrack(i);
    if(!myTrack) continue;
    LoadTrack(myTrack);
    FillTrackSpy("TrkAll");
    if(!AcceptDaughter()) continue;
    FillTrackSpy("TrkSel");
    ((TH2D*)((TList*)fList->FindObject("TrkSel"))->FindObject("PtIPXY" ))->Fill( myTrack->Pt(), fDaughterImpactParameterXY );
    if( myTrack->Charge()>0 ) {
      pd0[np] = fDaughterImpactParameterXY;
      plist[np++] = i;
    } else {
      nd0[nn] = fDaughterImpactParameterXY;
      nlist[nn++] = i;
    }
  }
  ((TH1D*)((TList*)fList->FindObject("TrkSel"))->FindObject("NPAIR" ))->Fill( np,nn );
  const AliESDVertex *vtx = tESD->GetPrimaryVertex();
  AliESDtrack *pT, *nT;
  for(int p=0; p!=np; ++p) {
    pT = (AliESDtrack*) tESD->GetTrack( plist[p] );
    for(int n=0; n!=nn; ++n) {
      nT = (AliESDtrack*) tESD->GetTrack( nlist[n] );
      fDecayProductIPXY = pd0[p]*nd0[n];
      AliExternalTrackParam pETP(*pT), nETP(*nT);
      Double_t xa, xb;
      pETP.GetDCA(&nETP,tESD->GetMagneticField(),xa,xb);
      fDecayDCAdaughters = pETP.PropagateToDCA(&nETP,tESD->GetMagneticField());
      AliESDv0 vertex( nETP,nlist[n], pETP,plist[p] );
      fDecayCosinePointingAngleXY = CosThetaPointXY( &vertex, vtx );
      fDecayRadXY = DecayLengthXY( &vertex, vtx );
      fDecayPt = vertex.Pt();
      fDecayPhi = vertex.Phi();
      fDecayEta = vertex.Eta();
      Double_t pmx, pmy, pmz, nmx, nmy, nmz;
      vertex.GetNPxPyPz(nmx,nmy,nmz);
      vertex.GetPPxPyPz(pmx,pmy,pmz);
      TVector3 mom1(pmx,pmy,pmz), mom2(nmx,nmy,nmz), mom(vertex.Px(),vertex.Py(),vertex.Pz());
      Double_t qlpos = mom1.Dot(mom)/mom.Mag();
      Double_t qlneg = mom2.Dot(mom)/mom.Mag();
      fDecayQt = mom1.Perp(mom);
      fDecayAlpha = (qlpos-qlneg)/(qlpos+qlneg);
      Double_t mpi = 0.13957018;
      if(fSpecie==0) {
        Double_t eppi = TMath::Sqrt( mpi*mpi + pmx*pmx + pmy*pmy + pmz*pmz );
        Double_t enpi = TMath::Sqrt( mpi*mpi + nmx*nmx + nmy*nmy + nmz*nmz );
        fDecayMass = TMath::Sqrt( mpi*mpi + mpi*mpi + 2*(eppi*enpi - pmx*nmx - pmy*nmy - pmz*nmz ) );
        fDecayRapidity = vertex.RapK0Short();
      } else {
        Double_t mpr = 0.938272013;
        Double_t epi, epr;
        if(fDecayAlpha>0) {
          epr = TMath::Sqrt( mpr*mpr + pmx*pmx + pmy*pmy + pmz*pmz );
          epi = TMath::Sqrt( mpi*mpi + nmx*nmx + nmy*nmy + nmz*nmz );
        } else {
          epi = TMath::Sqrt( mpi*mpi + pmx*pmx + pmy*pmy + pmz*pmz );
          epr = TMath::Sqrt( mpr*mpr + nmx*nmx + nmy*nmy + nmz*nmz );
        }
        fDecayMass = TMath::Sqrt( mpi*mpi + mpr*mpr + 2*(epi*epr - pmx*nmx - pmy*nmy - pmz*nmz ) );
        fDecayRapidity = vertex.RapLambda();
      }
      Double_t energy = TMath::Sqrt( fDecayMass*fDecayMass + vertex.Px()*vertex.Px() + vertex.Py()*vertex.Py() + vertex.Pz()*vertex.Pz() );
      Double_t gamma = energy/fDecayMass;
      fDecayDecayLength = DecayLength( &vertex, vtx )/gamma;
      Double_t dPHI = fDecayPhi;
      Double_t dDPHI = dPHI - fPsi2;
      if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
      if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
      if(fQAlevel>1) {
        if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("RecAllOP");
        else FillCandidateSpy("RecAllIP");
      }
      FillCandidateSpy("RecAll");
      ((TH2D*)((TList*)fList->FindObject("RecAll"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
      ((TH2D*)((TList*)fList->FindObject("RecAll"))->FindObject("XPOSXNEG"))->Fill( xa, xb );
      if(!AcceptCandidate()) continue;
      if(fDecayMass<fMinMass) continue;
      if(fDecayMass>fMaxMass) continue;
      // PID missing
      if(fQAlevel>1) {
        if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("RecSelOP");
        else FillCandidateSpy("RecSelIP");
      }
      FillCandidateSpy("RecSel");
      ((TH2D*)((TList*)fList->FindObject("RecSel"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
      ((TH2D*)((TList*)fList->FindObject("RecSel"))->FindObject("XPOSXNEG"))->Fill( xa, xb );

      fDecayIDneg = nT->GetID();
      fDecayIDpos = pT->GetID();
      MakeTrack();
      LoadTrack(pT); FillTrackSpy("TrkDau");
      LoadTrack(nT); FillTrackSpy("TrkDau");

      //===== BEGIN OF MCMATCH
      if(stack) {
        bool matched = false;
        Int_t labelpos = pT->GetLabel();
        Int_t labelneg = nT->GetLabel();
        Double_t rOri=-1;
        if( labelpos>0 && labelneg>0 ) {
          TParticle *mcpos = stack->Particle( labelpos );
          TParticle *mcneg = stack->Particle( labelneg );
          Int_t pdgRecPos = mcpos->GetPdgCode();
          Int_t pdgRecNeg = mcneg->GetPdgCode();
          if( pdgRecPos==211&&pdgRecNeg==-211 ) if(mcpos->GetMother(0)>0) {
            if( mcpos->GetMother(0)==mcneg->GetMother(0) ) {
              TParticle *mcmot = stack->Particle( mcpos->GetMother(0) );
              rOri = TMath::Sqrt( mcmot->Vx()*mcmot->Vx() + mcmot->Vy()*mcmot->Vy() );
              if( TMath::Abs(mcmot->GetPdgCode())==310) {
                if(mcmot->GetNDaughters()==2) matched=true;
              }
            }
          }
        }
        if(matched) {
          FillCandidateSpy("RecMth");
          ((TH2D*)((TList*)fList->FindObject("RecMth"))->FindObject("D0PD0N"))->Fill( pd0[p],nd0[n] );
          ((TH2D*)((TList*)fList->FindObject("RecMth"))->FindObject("XPOSXNEG"))->Fill( xa, xb );
          ((TH1D*)((TList*)fList->FindObject("RecMth"))->FindObject("MCOrigin"))->Fill( rOri );
          LoadTrack(pT); FillTrackSpy("TrkMth");
          LoadTrack(nT); FillTrackSpy("TrkMth");
        }
      }
      //===== END OF MCMATCH
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadStack(TClonesArray* mcArray) {
  if(!mcArray) return;
  AliAODMCParticle *myMCTrack, *iMCDau, *jMCDau;
  for(int i=0; i!=mcArray->GetEntriesFast(); ++i) {
    myMCTrack = dynamic_cast<AliAODMCParticle*>(mcArray->At( i ));
    if(!myMCTrack) continue;
    int tPDG=310;
    if(fSpecie>0) tPDG = 3122;
    if( TMath::Abs(myMCTrack->PdgCode())==tPDG )
      if( myMCTrack->GetNDaughters() == 2 ) {
        Int_t iDau = myMCTrack->GetDaughter(0);
        Int_t jDau = myMCTrack->GetDaughter(1);
        AliAODMCParticle *posDau=NULL;
        AliAODMCParticle *negDau=NULL;
        if(iDau>0&&jDau>0) {
          iMCDau = dynamic_cast<AliAODMCParticle*>(mcArray->At( iDau ));
          jMCDau = dynamic_cast<AliAODMCParticle*>(mcArray->At( jDau ));
          if(iMCDau) {
            if(iMCDau->Charge()>0) posDau=iMCDau;
            else negDau=iMCDau;
          }
          if(jMCDau) {
            if(jMCDau->Charge()>0) posDau=jMCDau;
            else negDau=jMCDau;
          }
        } //got two daughters
        if(posDau&&negDau) {
          Double_t dx = myMCTrack->Xv() - posDau->Xv();
          Double_t dy = myMCTrack->Yv() - posDau->Yv();
          Double_t dz = myMCTrack->Zv() - posDau->Zv();
          fDecayRadXY = TMath::Sqrt( dx*dx + dy*dy );
          TVector3 momPos(posDau->Px(),posDau->Py(),posDau->Pz());
          TVector3 momNeg(negDau->Px(),negDau->Py(),negDau->Pz());
          TVector3 momTot(myMCTrack->Px(),myMCTrack->Py(),myMCTrack->Pz());
          Double_t qlpos = momPos.Dot(momTot)/momTot.Mag();
          Double_t qlneg = momNeg.Dot(momTot)/momTot.Mag();
          fDecayQt = momPos.Perp(momTot);
          fDecayAlpha = 1.-2./(1.+qlpos/qlneg);
          fDecayMass = myMCTrack->GetCalcMass();
          Double_t energy = myMCTrack->E();
          Double_t gamma = energy/fDecayMass;
          fDecayDecayLength = TMath::Sqrt(dx*dx+dy*dy+dz*dz)/gamma;
          fDecayPt = myMCTrack->Pt();
          fDecayPhi = myMCTrack->Phi();
          fDecayEta = myMCTrack->Eta();
          fDecayRapidity = myMCTrack->Y();
          fDecayDCAdaughters = 0;
          fDecayCosinePointingAngleXY = 1;
          fDecayProductIPXY = -1;
          if(AcceptCandidate()) FillCandidateSpy("GenTru");
        }
      } // k0/lda with two daughters
    //==== BEGIN TRACK CUTS
    if(myMCTrack->Eta()<-0.8) continue;
    if(myMCTrack->Eta()>+0.8) continue;
    //==== END TRACK CUTS
    switch( TMath::Abs(myMCTrack->PdgCode()) ) {
    case (310): //k0s
      FillMCParticleSpy( "MCTK0s", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTK0sGenAcc", myMCTrack );
      break;
    case (3122): //lda
      FillMCParticleSpy( "MCTLda", myMCTrack );
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTLdaGenAcc", myMCTrack );
      break;
    case (333): //phi
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTPhiGenAcc", myMCTrack );
      break;
    case (3312): //xi
      if( myMCTrack->IsPrimary() )
        FillMCParticleSpy( "MCTXiGenAcc", myMCTrack );
      break;
    }
    
  }
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::CosThetaPointXY(AliESDv0 *me, const AliVVertex *vtx) {
  TVector3 mom( me->Px(), me->Py(), 0 );
  TVector3 fli( me->Xv()-vtx->GetX(), me->Yv()-vtx->GetY(), 0 );
  Double_t ctp = mom.Dot(fli) / mom.Mag() / fli.Mag();
  return ctp;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::CosThetaPointXY(AliAODv0 *me, const AliVVertex *vtx) {
  TVector3 mom( me->Px(), me->Py(), 0 );
  TVector3 fli( me->Xv()-vtx->GetX(), me->Yv()-vtx->GetY(), 0 );
  Double_t ctp = mom.Dot(fli) / mom.Mag() / fli.Mag();
  return ctp;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLengthXY(AliESDv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLengthXY(AliAODv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLength(AliESDv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dz = me->Zv()-vtx->GetZ();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy + dz*dz );
  return dxy;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::DecayLength(AliAODv0 *me, const AliVVertex *vtx) {
  Double_t dx = me->Xv()-vtx->GetX();
  Double_t dy = me->Yv()-vtx->GetY();
  Double_t dz = me->Zv()-vtx->GetZ();
  Double_t dxy = TMath::Sqrt( dx*dx + dy*dy + dz*dz );
  return dxy;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ReadFromAODv0(AliAODEvent *tAOD) {
  TClonesArray* mcArray=NULL;
  if(fReadMC) {
    mcArray = dynamic_cast<TClonesArray*>(tAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    ReadStack(mcArray);
  }

  Int_t nV0s = tAOD->GetNumberOfV0s();
  AliAODv0 *myV0;
  Int_t v0all=0, v0imw=0;
  for (Int_t i=0; i!=nV0s; ++i) {
    myV0 = (AliAODv0*) tAOD->GetV0(i);
    if(!myV0) continue;
    if(!fOnline) if(myV0->GetOnFlyStatus() ) continue;
    if(fOnline) if(!myV0->GetOnFlyStatus() ) continue;
    AliAODTrack *iT, *jT;
    AliAODVertex *vtx = tAOD->GetPrimaryVertex();
    Double_t pos[3],cov[6];
    vtx->GetXYZ(pos);
    vtx->GetCovarianceMatrix(cov);
    const AliESDVertex vESD(pos,cov,100.,100);
    // TESTING CHARGE
    int iPos, iNeg;
    iT=(AliAODTrack*) myV0->GetDaughter(0);
    if(iT->Charge()>0) {
      iPos = 0; iNeg = 1;
    } else {
      iPos = 1; iNeg = 0;
    }
    // END OF TEST

    iT=(AliAODTrack*) myV0->GetDaughter(iPos); // positive
    AliESDtrack ieT( iT );
    ieT.SetTPCClusterMap( iT->GetTPCClusterMap() );
    ieT.SetTPCSharedMap( iT->GetTPCSharedMap() );
    ieT.SetTPCPointsF( iT->GetTPCNclsF() );
    ieT.PropagateToDCA(&vESD, tAOD->GetMagneticField(), 100);
    LoadTrack(&ieT,iT->Chi2perNDF());
    Float_t ip[2];
    ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    fDecayIPpos = fDaughterImpactParameterXY; //ieT.GetD(pos[0], pos[1], tAOD->GetMagneticField());
    if(!AcceptDaughter()) continue;

    jT=(AliAODTrack*) myV0->GetDaughter(iNeg); // negative
    AliESDtrack jeT( jT );
    jeT.SetTPCClusterMap( jT->GetTPCClusterMap() );
    jeT.SetTPCSharedMap( jT->GetTPCSharedMap() );
    jeT.SetTPCPointsF( jT->GetTPCNclsF() );
    jeT.PropagateToDCA(&vESD, tAOD->GetMagneticField(), 100);
    LoadTrack(&jeT,jT->Chi2perNDF());
    jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    fDecayIPneg = fDaughterImpactParameterXY; //jeT.GetD(pos[0], pos[1], tAOD->GetMagneticField());
    if(!AcceptDaughter()) continue;

    if( fExcludeTPCEdges ) {
      if( IsAtTPCEdge(iT->Phi(),iT->Pt(),+1,tAOD->GetMagneticField()) ) continue;
      if( IsAtTPCEdge(jT->Phi(),jT->Pt(),-1,tAOD->GetMagneticField()) ) continue;
    }
    ieT.GetDCA(&jeT,tAOD->GetMagneticField(),fDecayXpos,fDecayXneg);
    /*
    // cutting out population close to TPC edges :: strange excess saw in 2010
    if( fExcludeTPCEdges ) {
    Double_t phimod = myV0->Phi();
    int sectors[6] = {5,6,9,10,11,12};
    for(int ii=0; ii!=6; ++ii)
    if( (phimod<(sectors[ii]+1)*TMath::Pi()/9.0) && (phimod>sectors[ii]*TMath::Pi()/9.0) )
    return 0;
    }
    */
    if(fSpecie==0)
      fDecayRapidity = myV0->RapK0Short();
    else
      fDecayRapidity = myV0->RapLambda();
    fDecayDCAdaughters = myV0->DcaV0Daughters();
    fDecayCosinePointingAngleXY = CosThetaPointXY( myV0, vtx );
    fDecayRadXY = DecayLengthXY( myV0, vtx );
    fDecayPt = myV0->Pt();
    fDecayPhi = myV0->Phi();
    fDecayEta = myV0->Eta();
    fDecayProductIPXY = fDecayIPpos*fDecayIPneg;
    fDecayQt = myV0->PtArmV0();
    fDecayAlpha = myV0->AlphaV0(); // AlphaV0 -> AODRecoDecat::Alpha -> return 1.-2./(1.+QlProng(0)/QlProng(1));
    if(myV0->ChargeProng(iPos)<0) fDecayAlpha = -fDecayAlpha; // protects for a change in convention
    fDecayPt = myV0->Pt();
    fDecayEta = myV0->Eta();
    if( fSpecie==0 ) {
      fDecayMass = myV0->MassK0Short();
    } else {
      if(fDecayAlpha>0) fDecayMass = myV0->MassLambda();
      else fDecayMass = myV0->MassAntiLambda();
    }
    v0all++;
    if(fDecayMass<fMinMass) continue;
    if(fDecayMass>fMaxMass) continue;
    v0imw++;
    Double_t energy = TMath::Sqrt( fDecayMass*fDecayMass + myV0->Px()*myV0->Px() + myV0->Py()*myV0->Py() + myV0->Pz()*myV0->Pz() );
    Double_t gamma = energy/fDecayMass;
    fDecayDecayLength = DecayLength( myV0, vtx )/gamma;
    Double_t dPHI = fDecayPhi;
    Double_t dDPHI = dPHI - fPsi2;
    if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
    if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
    if(fQAlevel>1) {
      if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("RecAllOP");
      else FillCandidateSpy("RecAllIP");
    }
    FillCandidateSpy("RecAll");
    if(!AcceptCandidate()) continue;
    //PID for lambda::proton
    if( fSpecie>0 )
      if(fDecayPt<1.2) {
        if(fDecayAlpha>0) {
          if( !PassesPIDCuts(&ieT) ) continue; //positive track
        } else {
          if( !PassesPIDCuts(&jeT) ) continue; //negative track
        }
      }
    if(fQAlevel>1) {
      if( (dDPHI>TMath::PiOver4()) && (dDPHI<3*TMath::PiOver4()) ) FillCandidateSpy("RecSelOP");
      else FillCandidateSpy("RecSelIP");
    }
    FillCandidateSpy("RecSel");
    // ============================
    // Posting for FlowAnalysis
    if(!fPostMatched) {
      fDecayIDneg = iT->GetID();
      fDecayIDpos = jT->GetID();
      MakeTrack();
    }
    // ============================

    LoadTrack(&ieT,iT->Chi2perNDF());
    ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    FillTrackSpy("TrkDau");
    LoadTrack(&jeT,jT->Chi2perNDF()); 
    jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    FillTrackSpy("TrkDau");
    //===== BEGIN OF MCMATCH
    if(fReadMC) ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 1 ); // Selected event
    if(mcArray) {
      ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 2 ); // Stack found
      bool matched = false;
      bool feeddown = false;
      Int_t labelpos = iT->GetLabel();
      Int_t labelneg = jT->GetLabel();
      AliAODMCParticle *mcpos = (AliAODMCParticle*) mcArray->At( TMath::Abs(labelpos) );
      AliAODMCParticle *mcneg = (AliAODMCParticle*) mcArray->At( TMath::Abs(labelneg) );
      if( mcpos && mcneg ) {
        ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 3 ); // Daughters in stack
        Int_t pdgRecPos = mcpos->GetPdgCode();
        Int_t pdgRecNeg = mcneg->GetPdgCode();
        int pospdg=211, negpdg=211;
        int mompdg=310, fdwpdg=333;
        if(fSpecie>0) {
          mompdg=3122;
          fdwpdg=3312;
          if(fDecayAlpha>0) {
            pospdg=2212; negpdg=211;
          } else {
            negpdg=2212; pospdg=211;
          }
        }
        if( TMath::Abs(pdgRecPos)==pospdg&&TMath::Abs(pdgRecNeg)==negpdg )
          if(mcpos->GetMother()>-1)
            if( mcpos->GetMother()==mcneg->GetMother() ) {
              AliAODMCParticle *mcmot = (AliAODMCParticle*) mcArray->At( mcpos->GetMother() );
              fDecayMatchOrigin = TMath::Sqrt( mcmot->Xv()*mcmot->Xv() + mcmot->Yv()*mcmot->Yv() );
              fDecayMatchPt = mcmot->Pt();
              fDecayMatchEta = mcmot->Eta();
              fDecayMatchPhi = mcmot->Phi();
              if( TMath::Abs(mcmot->GetPdgCode())==mompdg) {
                if(mcmot->GetNDaughters()==2) {
		  ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 4 ); // Correspond to decay
                  matched=true;
                  Double_t dx = mcmot->Xv() - mcpos->Xv();
                  Double_t dy = mcmot->Yv() - mcpos->Yv();
                  fDecayMatchRadXY = TMath::Sqrt( dx*dx + dy*dy );
                }
                if(mcmot->GetMother()>-1) {
		  ((TH1D*)((TList*)fList->FindObject("STATMC"))->FindObject("Events"))->Fill( 5 ); // Decay has mother
                  AliAODMCParticle *mcfdw = (AliAODMCParticle*) mcArray->At( mcmot->GetMother() );
                  if( TMath::Abs(mcfdw->GetPdgCode())==fdwpdg)
                    feeddown=true;
                } // k0/lda have mother
              } // mother matches k0/lda
            } // both have same mother
      }
      if(matched) {
        FillCandidateSpy("RecMth",true);
	if(fPostMatched>0) {
	  fDecayIDneg = iT->GetID();
	  fDecayIDpos = jT->GetID();
	  MakeTrack();
	}
	if(labelpos<0&&labelneg<0) {
	  FillCandidateSpy("NegNegMth",true);
	} else if(labelpos*labelneg<0) {
	  FillCandidateSpy("NegPosMth",true);
	}
        LoadTrack(&ieT,iT->Chi2perNDF());
	fDaughterMatchPhi=mcpos->Phi();
	fDaughterMatchEta=mcpos->Eta();
	fDaughterMatchPt=mcpos->Pt();
        ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("TrkMth",true);
	if(labelpos<0||labelneg<0) FillTrackSpy("TrkNegMth",true);
	else FillTrackSpy("TrkPosMth",true);
        LoadTrack(&jeT,jT->Chi2perNDF());
	fDaughterMatchPhi=mcneg->Phi();
	fDaughterMatchEta=mcneg->Eta();
	fDaughterMatchPt=mcneg->Pt();
        jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("TrkMth",true);
	if(labelpos<0||labelneg<0) FillTrackSpy("TrkNegMth",true);
	else FillTrackSpy("TrkPosMth",true);
      } else {
        FillCandidateSpy("RecUNMth",false);
	if(fPostMatched<0) {
	  fDecayIDneg = iT->GetID();
	  fDecayIDpos = jT->GetID();
	  MakeTrack();
	}
        LoadTrack(&ieT,iT->Chi2perNDF());
        ieT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("TrkUNMth",false);
        LoadTrack(&jeT,jT->Chi2perNDF());
        jeT.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
        fDaughterImpactParameterXY = ip[0];
        fDaughterImpactParameterZ = ip[1];
        FillTrackSpy("TrkUNMth",false);
      }
      if(feeddown) {
        FillCandidateSpy("MthFDW",true);
      }
    }
    //===== END OF MCMATCH
  }
  ((TH2D*)((TList*)fList->FindObject("RecAll"))->FindObject("V0SADC"))->Fill( v0all,v0imw );
  return;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesPIDCuts(AliESDtrack *myTrack, AliPID::EParticleType pid) {
  Bool_t pass=kTRUE;
  if(fPIDResponse) {
    fDaughterNSigmaPID = fPIDResponse->NumberOfSigmasTPC(myTrack,pid);
    if( TMath::Abs(fDaughterNSigmaPID) > fDaughterMaxNSigmaPID )
      pass = kFALSE;
  }
  return pass;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ChargeParticles(AliAODEvent *tAOD) {
  //benchmark purposes (for the ultimate measurement for LHC10h see Alex, Francesco)
  if(!tAOD) return;
  // (for the moment) only compatible with SPVZE (no untagging)
  for(int i=0; i!=tAOD->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = tAOD->GetTrack( i );
    if(!t) continue;
    if( !t->TestFilterBit(1024) ) continue;
    fDecayMass=0.0; // using mass as nsigmas control plot
    if(fPIDResponse) { // PID
      switch(fSpecie) { // TPC PID only
      case(kPION):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kPion);
        break;
      case(kKAON):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kKaon);
        break;
      case(kPROTON):
        fDecayMass = fPIDResponse->NumberOfSigmasTPC(t,AliPID::kProton);
        break;
      }
    }
    Bool_t pass = kTRUE;
    if( TMath::Abs(fDecayMass) > 3.0 ) pass=kFALSE;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) pass=kFALSE;
    if( t->Pt()<0.2 || t->Pt()>20.0 ) pass=kFALSE;
    if( t->GetTPCsignal()<10.0 ) pass=kFALSE;
    fDecayPt=t->Pt();
    fDecayPhi=t->Phi();
    fDecayEta=t->Eta();
    fDecayID=t->GetID();

    FillCandidateSpy("RecAll");
    if(!pass) continue;

    AliESDtrack et( t );
    et.SetTPCClusterMap( t->GetTPCClusterMap() );
    et.SetTPCSharedMap( t->GetTPCSharedMap() );
    et.SetTPCPointsF( t->GetTPCNclsF() );
    Float_t ip[2];
    LoadTrack(&et,t->Chi2perNDF()); 
    AliAODVertex *vtx = tAOD->GetPrimaryVertex();
    Double_t pos[3];
    vtx->GetXYZ(pos);
    et.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    fDaughterImpactParameterXY = ip[0];
    fDaughterImpactParameterZ = ip[1];
    FillTrackSpy("TrkDau");
    FillCandidateSpy("RecSel");

    MakeTrack();
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::Terminate(Option_t *) {
  //terminate
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeTrack() {
  // create track for flow tasks
  if(fCandidates->GetLast()+5>fCandidates->GetSize()) {
    fCandidates->Expand( fCandidates->GetSize()+20 );
  }
  Bool_t overwrite = kTRUE;
  AliFlowCandidateTrack *oTrack = (static_cast<AliFlowCandidateTrack*> (fCandidates->At( fCandidates->GetLast()+1 )));
  if( !oTrack ) { // creates new
    oTrack = new AliFlowCandidateTrack();
    overwrite = kFALSE;
  } else { // overwrites
    oTrack->ClearMe();
  }
  oTrack->SetMass(fDecayMass);
  oTrack->SetPt(fDecayPt);
  oTrack->SetPhi(fDecayPhi);
  oTrack->SetEta(fDecayEta);
  if(fSpecie<10) {
    oTrack->AddDaughter(fDecayIDpos);
    oTrack->AddDaughter(fDecayIDneg);
  } else {
    oTrack->SetID( fDecayID );
  }
  oTrack->SetForPOISelection(kTRUE);
  oTrack->SetForRPSelection(kFALSE);
  if(overwrite) {
    fCandidates->SetLast( fCandidates->GetLast()+1 );
  } else {
    fCandidates->AddLast(oTrack);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidates() {
  // adds candidates to flow events (untaging if necessary)
  if(fDebug) printf("FlowEventTPC %d tracks | %d RFP | %d POI\n",fTPCevent->NumberOfTracks(),fTPCevent->GetNumberOfRPs(),fTPCevent->GetNumberOfPOIs());
  if(fDebug) printf("FlowEventVZE %d tracks | %d RFP | %d POI\n",fVZEevent->NumberOfTracks(),fVZEevent->GetNumberOfRPs(),fVZEevent->GetNumberOfPOIs());
  if(fDebug) printf("I received %d candidates\n",fCandidates->GetEntriesFast());
  Int_t untagged=0;
  Int_t poi=0;
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = static_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    cand->SetForPOISelection(kTRUE);
    cand->SetForRPSelection(kFALSE);
    poi++;
    if(fDebug) printf(" >Checking at candidate %d with %d daughters: mass %f\n",
                      iCand,cand->GetNDaughters(),cand->Mass());
    if(fSpecie<10) { // DECAYS
      // untagging ===>
      if(fDaughterUnTag) {
	for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) {
	  if(fDebug) printf("  >Daughter %d with fID %d", iDau, cand->GetIDDaughter(iDau));
	  for(int iRPs=0; iRPs!=fTPCevent->NumberOfTracks(); ++iRPs ) {
	    AliFlowTrack *iRP = static_cast<AliFlowTrack*>(fTPCevent->GetTrack( iRPs ));
	    if(!iRP) continue;
	    if(!iRP->InRPSelection()) continue;
	    if(cand->GetIDDaughter(iDau) == iRP->GetID()) {
	      if(fDebug) printf(" was in RP set");
	      ++untagged;
	      iRP->SetForRPSelection(kFALSE);
	      fTPCevent->SetNumberOfRPs( fTPCevent->GetNumberOfRPs() -1 );
	    }
	  }
	  if(fDebug) printf("\n");
	}
      }
      // <=== untagging 
      fTPCevent->InsertTrack( ((AliFlowTrack*) cand) );
    } else {  // CHARGED
      // adding only new tracks and tagging accordingly ===>
      Bool_t found=kFALSE;
      for(int iRPs=0; iRPs!=fTPCevent->NumberOfTracks(); ++iRPs ) {
        AliFlowTrack *iRP = static_cast<AliFlowTrack*>(fTPCevent->GetTrack( iRPs ));
        if(!iRP) continue;
        if(!iRP->InRPSelection()) continue;
        if(cand->GetID() == iRP->GetID()) {
          if(fDebug) printf("  >charged track (%d) was also found in RP set (adding poi tag)\n",cand->GetID());
          iRP->SetForPOISelection(kTRUE);
          found = kTRUE;
        }
      }
      if(!found) // not found adding track
        fTPCevent->InsertTrack( ((AliFlowTrack*) cand) );
    }
    fVZEevent->InsertTrack( ((AliFlowTrack*) cand) );
  } //END OF LOOP
  fTPCevent->SetNumberOfPOIs( poi );
  fVZEevent->SetNumberOfPOIs( poi );
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("POI"))->Fill( poi );
  ((TH1D*)((TList*)fList->FindObject("Event"))->FindObject("UNTAG"))->Fill( untagged );
  if(fDebug) printf("FlowEventTPC %d tracks | %d RFP | %d POI\n",fTPCevent->NumberOfTracks(),fTPCevent->GetNumberOfRPs(),fTPCevent->GetNumberOfPOIs());
  if(fDebug) printf("FlowEventVZE %d tracks | %d RFP | %d POI\n",fVZEevent->NumberOfTracks(),fVZEevent->GetNumberOfRPs(),fVZEevent->GetNumberOfPOIs());
}
//=======================================================================
void AliAnalysisTaskFlowStrange::PushBackFlowTrack(AliFlowEvent *flowevent, Double_t pt, Double_t phi, Double_t eta, Double_t w, Int_t id) {
  AliFlowTrack rfp;
  rfp.SetPt(pt);
  rfp.SetPhi(phi);
  rfp.SetEta(eta);
  rfp.SetWeight(w);
  rfp.SetForRPSelection(kTRUE);
  rfp.SetForPOISelection(kFALSE);
  rfp.SetID(id);
  flowevent->InsertTrack( &rfp );
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::IsAtTPCEdge(Double_t phi,Double_t pt,Int_t charge,Double_t b) {
  // Origin: Alex Dobrin
  // Implemented by Carlos Perez
  TF1 cutLo("cutLo", "-0.01/x+pi/18.0-0.015", 0, 100);
  TF1 cutHi("cutHi", "0.55/x/x+pi/18.0+0.03", 0, 100);
  Double_t phimod = phi;
  if(b<0) phimod = TMath::TwoPi()-phimod;  //for negatve polarity field
  if(charge<0) phimod = TMath::TwoPi()-phimod; //for negatve charge
  phimod += TMath::Pi()/18.0;
  phimod = fmod(phimod, TMath::Pi()/9.0);
  if( phimod<cutHi.Eval(pt) && phimod>cutLo.Eval(pt) )
    return kTRUE;

  return kFALSE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQVectors() {
  //computes event plane and updates fPsi2
  //if there is a problem fPsi->-1
  fPsi2=-1;
  Double_t vzeqax, vzeqay, vzeqaw, vzeqbx, vzeqby, vzeqbw;
  Double_t tpcqax, tpcqay, tpcqaw, tpcqbx, tpcqby, tpcqbw;
  //=>loading event
  MakeQVZE(InputEvent(),vzeqax,vzeqay,vzeqaw,vzeqbx,vzeqby,vzeqbw);
  MakeQTPC(InputEvent(),tpcqax,tpcqay,tpcqaw,tpcqbx,tpcqby,tpcqbw);
  if(fReadMC) {
    fVZEevent->SetMCReactionPlaneAngle( fMCEP );    
    fTPCevent->SetMCReactionPlaneAngle( fMCEP );    
  }
  //=>computing psi
  //VZERO
  Double_t vqx, vqy;
  vqx=vqy=0;
  Double_t psivzea,psivzeb,psivze,vzew;
  psivzea = ( TMath::Pi()+TMath::ATan2(-vzeqay,-vzeqax) )/2.0;
  psivzeb = ( TMath::Pi()+TMath::ATan2(-vzeqby,-vzeqbx) )/2.0;
  vqx = vzeqax + vzeqbx;
  vqy = vzeqay + vzeqby;
  vzew = vzeqaw + vzeqbw;
  psivze = ( TMath::Pi()+TMath::ATan2(-vqy,-vqx) )/2.0;
  //TPC
  Double_t tqx, tqy;
  tqx=tqy=0;
  Double_t psitpca,psitpcb,psitpc,tpcw;
  psitpca = ( TMath::Pi()+TMath::ATan2(-tpcqay,-tpcqax) )/2.0;
  psitpcb = ( TMath::Pi()+TMath::ATan2(-tpcqby,-tpcqbx) )/2.0;
  tqx = tpcqax + tpcqbx;
  tqy = tpcqay + tpcqby;
  tpcw = tpcqaw + tpcqbw;
  psitpc = ( TMath::Pi()+TMath::ATan2(-tqy,-tqx) )/2.0;
  //=>does the event clear?
  switch(fWhichPsi) {
  case(1): //VZERO
    if( (vzeqaw<2)||(vzeqbw<2) ) return;
    fPsi2 = psivze;
    break;
  case(2): //TPC
    if( (tpcqaw<2)||(tpcqbw<2) ) return;
    fPsi2 = psitpc;
    break;
  }
  //=>great! recording
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSI"))->Fill( psivze );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSIA"))->Fill( psivzea );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEPSIB"))->Fill( psivzeb );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("RFPVZE"))->Fill( vzew );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmVZE"))->Fill( TMath::Sqrt(vqx*vqx+vqy*vqy)/vzew );
  //------
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSI"))->Fill( psitpc );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSIA"))->Fill( psitpca );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCPSIB"))->Fill( psitpcb );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("RFPTPC"))->Fill( tpcw );
  ((TH1D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("QmTPC"))->Fill( TMath::Sqrt(tqx*tqx+tqy*tqy)/tpcw );
  //------
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 1., tpcqay/tpcqaw, tpcqaw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 2., tpcqax/tpcqaw, tpcqaw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 3., tpcqby/tpcqbw, tpcqbw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 4., tpcqbx/tpcqbw, tpcqbw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 5., tqy/tpcw, tpcw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCQ"))->Fill( 6., tqx/tpcw, tpcw );
  //------
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 1., vzeqay/vzeqaw, vzeqaw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 2., vzeqax/vzeqaw, vzeqaw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 3., vzeqby/vzeqbw, vzeqbw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 4., vzeqbx/vzeqbw, vzeqbw );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 5., vqy/vzew, vzew );
  ((TProfile*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEQ"))->Fill( 6., vqx/vzew, vzew );

  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQVZE(AliVEvent *tevent,
                                          Double_t &qxa,Double_t &qya,Double_t &qwa,
                                          Double_t &qxb,Double_t &qyb,Double_t &qwb) {
  //=>cleaning
  if(fUseFP) fVZEevent->ClearFast();
  //=>external weights
  Double_t extW[64];
  for(int i=0;i!=64;++i) extW[i]=1;
  if((!fVZEsave)&&(fVZEResponse)) {
    Double_t minC = fCentPerMin, maxC = fCentPerMax;
    if(fVZEmb) {
      minC = 0;
      maxC = 80;
    }
    Int_t ybinmin = fVZEResponse->GetYaxis()->FindBin(minC+1e-6);
    Int_t ybinmax = fVZEResponse->GetYaxis()->FindBin(maxC-1e-6);
    for(int i=0;i!=64;++i) extW[i] = fVZEResponse->Integral(i+1,i+1,ybinmin,ybinmax)/(maxC-minC);
    //ring-wise normalization
    Double_t ring[8];
    for(int j=0; j!=8; ++j) {
      ring[j]=0;
      for(int i=0;i!=8;++i) ring[j] += extW[j*8+i]/8;
    }
    //disk-wise normalization
    Double_t disk[2];
    int xbinmin, xbinmax;
    xbinmin = 1+8*fVZECa;
    xbinmax = 8+8*fVZECb;
    disk[0] = fVZEResponse->Integral(xbinmin,xbinmax,ybinmin,ybinmax)/(maxC-minC)/(xbinmax-xbinmin+1);
    xbinmin = 33+8*fVZEAa;
    xbinmax = 40+8*fVZEAb;
    disk[1] = fVZEResponse->Integral(xbinmin,xbinmax,ybinmin,ybinmax)/(maxC-minC)/(xbinmax-xbinmin+1);
    //for(int i=0;i!=64;++i) printf("CELL %d -> W = %f ||",i,extW[i]);

    if(fVZEByDisk) {
      for(int i=0;i!=64;++i) extW[i] = disk[i/32]/extW[i];
    } else {
      for(int i=0;i!=64;++i) extW[i] = ring[i/8]/extW[i];
    }
    //for(int i=0;i!=64;++i) printf(" W = %f \n",extW[i]);
  }
  //=>computing
  qxa=qya=qwa=qxb=qyb=qwb=0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //v0c -> qa
  for(int id=fVZECa*8;id!=8+fVZECb*8;++id) {
    eta = -3.45+0.5*(id/8);
    phi = TMath::PiOver4()*(0.5+id%8);
    w = tevent->GetVZEROEqMultiplicity(id)*extW[id];
    qxa += w*TMath::Cos(2*phi);
    qya += w*TMath::Sin(2*phi);
    qwa += w;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fVZEevent,0,phi,eta,w,0);
  }
  //v0a -> qb
  for(int id=32+fVZEAa*8;id!=40+fVZEAb*8;++id) {
    eta = +4.8-0.6*((id/8)-4);
    phi = TMath::PiOver4()*(0.5+id%8);
    w = tevent->GetVZEROEqMultiplicity(id)*extW[id];
    qxb += w*TMath::Cos(2*phi);
    qyb += w*TMath::Sin(2*phi);
    qwb += w;
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("VZEAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fVZEevent,0,phi,eta,w,0);
  }
  if(fUseFP) fVZEevent->SetNumberOfRPs(rfp);
  if(fDebug>0&&fUseFP) printf("Inserted tracks in FlowEventVZE %d ==> %.1f\n",fVZEevent->NumberOfTracks(),qwa+qwb);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddTPCRFPSpy(TList *me) {
  TH1D *tH1D;
  tH1D = new TH1D("PT",   "PT",           50,0,5);     me->Add(tH1D);
  tH1D = new TH1D("PHI",  "PHI", 90,0,TMath::TwoPi()); me->Add(tH1D);
  tH1D = new TH1D("ETA",  "ETA",          40,-1,+1);   me->Add(tH1D);
  tH1D = new TH1D("TPCS", "TPC Signal",   100,0,500);  me->Add(tH1D);
  tH1D = new TH1D("IPXY", "IPXY",         100,-2,+2);  me->Add(tH1D);
  tH1D = new TH1D("IPZ",  "IPZ",          100,-2,+2);  me->Add(tH1D);
  // TPC
  tH1D = new TH1D("TPCNCLS", "NCLS", 170,-0.5,+169.5);   me->Add(tH1D);
  tH1D = new TH1D("TPCSHCL", "NSCLS / NCLS", 100,0,1);   me->Add(tH1D);
  tH1D = new TH1D("TPCFICL", "NCLS1I / NCLS",100,0,1);   me->Add(tH1D);
  tH1D = new TH1D("TPCXRNF", "XROW / NFCLS", 100,0,1.5); me->Add(tH1D);
  tH1D = new TH1D("TPCRCHI", "CHI2 / NCLS",  50,0,5);    me->Add(tH1D);
  // ITS
  tH1D = new TH1D("ITSNCLS", "NCLS",   7,-0.5,+6.5); me->Add(tH1D);
  tH1D = new TH1D("ITSRCHI", "CHI2 / NCLS", 50,0,5); me->Add(tH1D);

}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesRFPTPCCuts(AliESDtrack *track, Double_t aodchi2cls, Float_t aodipxy, Float_t aodipz) {
  if(track->GetKinkIndex(0)>0) return kFALSE;
  if( (track->GetStatus()&AliESDtrack::kTPCrefit)==0 ) return kFALSE;
  Double_t pt = track->Pt();
  Double_t phi = track->Phi();
  Double_t eta = track->Eta();
  Double_t tpcs = track->GetTPCsignal();
  Float_t ipxy, ipz;
  track->GetImpactParameters(ipxy,ipz);
  Int_t cls = track->GetTPCclusters(0);
  Double_t xrows, findcls, chi2;
  findcls = track->GetTPCNclsF();
  xrows = track->GetTPCCrossedRows();
  chi2 = track->GetTPCchi2();
  Double_t rchi2 = chi2/cls;
  if(!fReadESD) {
    rchi2 = aodchi2cls;
    ipxy = aodipxy;
    ipz = aodipz;
  }
  Double_t xrnfcls = xrows/findcls;
  Double_t scls, cls1i, itschi2;
  Int_t itscls;
  cls1i = track->GetTPCNclsIter1();
  scls = track->GetTPCnclsS();
  itscls = track->GetITSclusters(0);
  itschi2 = track->GetITSchi2();
  Double_t shcl = scls/cls;
  Double_t ficl = cls1i/cls;
  Double_t itsrchi2 = itscls/itschi2;
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("PT"))->Fill( pt );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("PHI"))->Fill( phi );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ETA"))->Fill( eta );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCS"))->Fill( tpcs );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("IPXY"))->Fill( ipxy );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("IPZ"))->Fill( ipz );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCNCLS"))->Fill( cls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCSHCL"))->Fill( shcl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCFICL"))->Fill( ficl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCXRNF"))->Fill( xrnfcls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("TPCRCHI"))->Fill( rchi2 );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ITSNCLS"))->Fill( itscls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPall"))->FindObject("ITSRCHI"))->Fill( itsrchi2 );
  if(pt<fRFPminPt) return kFALSE; //0.2
  if(pt>fRFPmaxPt) return kFALSE; //5.0
  if(eta<fRFPminEta) return kFALSE; //-0.8
  if(eta>fRFPmaxEta) return kFALSE; //+0.8
  if(tpcs<fRFPTPCsignal) return kFALSE; //10.0
  if( TMath::Sqrt(ipxy*ipxy/fRFPmaxIPxy/fRFPmaxIPxy+ipz*ipz/fRFPmaxIPz/fRFPmaxIPz)>1 ) return kFALSE; //2.4 3.2
  if(cls<fRFPTPCncls) return kFALSE; //70
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("PT"))->Fill( pt );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("PHI"))->Fill( phi );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ETA"))->Fill( eta );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCS"))->Fill( tpcs );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("IPXY"))->Fill( ipxy );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("IPZ"))->Fill( ipz );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCNCLS"))->Fill( cls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCSHCL"))->Fill( shcl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCFICL"))->Fill( ficl );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCXRNF"))->Fill( xrnfcls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("TPCRCHI"))->Fill( rchi2 );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ITSNCLS"))->Fill( itscls );
  ((TH1D*)((TList*)fList->FindObject("TPCRFPsel"))->FindObject("ITSRCHI"))->Fill( itsrchi2 );
  return kTRUE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliVEvent *tevent,
                                          Double_t &qxa,Double_t &qya,Double_t &qwa,
                                          Double_t &qxb,Double_t &qyb,Double_t &qwb) {
  AliESDEvent *tESD = (AliESDEvent*) (tevent);
  AliAODEvent *tAOD = (AliAODEvent*) (tevent);
  if(fReadESD) {
    if(!tESD) return;
    MakeQTPC(tESD,qxa,qya,qwa,qxb,qyb,qwb);
  } else {
    if(!tAOD) return;
    MakeQTPC(tAOD,qxa,qya,qwa,qxb,qyb,qwb);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliAODEvent *tAOD,
                                          Double_t &qxa,Double_t &qya,Double_t &qwa,
                                          Double_t &qxb,Double_t &qyb,Double_t &qwb) {
  //=>cleaning
  if(fUseFP) fTPCevent->ClearFast();
  qxa=qya=qwa=qxb=qyb=qwb=0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //=>aod stuff
  AliAODVertex *vtx = tAOD->GetPrimaryVertex();
  Double_t pos[3],cov[6];
  vtx->GetXYZ(pos);
  vtx->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);
  AliAODTrack *track;
  //=>looping
  Int_t rawN = tAOD->GetNumberOfTracks();
  for(Int_t id=0; id!=rawN; ++id) {
    track = tAOD->GetTrack(id);
    //=>cuts
    if(!track->TestFilterBit(fRFPFilterBit)) continue;
    if( fExcludeTPCEdges )
      if( IsAtTPCEdge( track->Phi(), track->Pt(), track->Charge(), tAOD->GetMagneticField() ) )	continue;
    AliESDtrack etrack( track );
    etrack.SetTPCClusterMap( track->GetTPCClusterMap() );
    etrack.SetTPCSharedMap( track->GetTPCSharedMap() );
    etrack.SetTPCPointsF( track->GetTPCNclsF() );
    Float_t ip[2];
    etrack.GetDZ(pos[0], pos[1], pos[2], tAOD->GetMagneticField(), ip);
    if(!PassesRFPTPCCuts(&etrack,track->Chi2perNDF(),ip[0],ip[1])) continue;
    //=>collecting info
    phi = track->Phi();
    eta = track->Eta();
    w = 1;
    //=>subevents
    if(eta<0.5*(fRFPminEta+fRFPmaxEta)) {
      qxa += w*TMath::Cos(2*phi);
      qya += w*TMath::Sin(2*phi);
      qwa += w;
    } else {
      qxb += w*TMath::Cos(2*phi);
      qyb += w*TMath::Sin(2*phi);
      qwb += w;
    }
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fTPCevent,track->Pt(),phi,eta,w,track->GetID());
  }
  if(fUseFP) fTPCevent->SetNumberOfRPs(rfp);
  if(fDebug) printf("Inserted tracks in FlowEventTPC %d ==> %.1f\n",fTPCevent->NumberOfTracks(),qwa+qwb);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeQTPC(AliESDEvent *tESD,
                                          Double_t &qxa,Double_t &qya,Double_t &qwa,
                                          Double_t &qxb,Double_t &qyb,Double_t &qwb) {
  //=>cleaning
  if(fUseFP) fTPCevent->ClearFast();
  qxa=qya=qwa=qxb=qyb=qwb=0;
  Int_t rfp=0;
  Double_t eta, phi, w;
  //=>looping
  AliESDtrack *track;
  Int_t rawN = tESD->GetNumberOfTracks();
  for(Int_t id=0; id!=rawN; ++id) {
    track = tESD->GetTrack(id);
    //=>cuts
    if( fExcludeTPCEdges )
      if( IsAtTPCEdge( track->Phi(), track->Pt(), track->Charge(), tESD->GetMagneticField() ) )	continue;
    if(!PassesFilterBit(track)) continue;
    if(!PassesRFPTPCCuts(track)) continue;
    //=>collecting info
    phi = track->Phi();
    eta = track->Eta();
    w = 1;
    //=>subevents
    if(eta<0) {
      qxa += w*TMath::Cos(2*phi);
      qya += w*TMath::Sin(2*phi);
      qwa += w;
    } else {
      qxb += w*TMath::Cos(2*phi);
      qyb += w*TMath::Sin(2*phi);
      qwb += w;
    }
    ((TH2D*)((TList*)fList->FindObject("MakeQSpy"))->FindObject("TPCAllPhiEta"))->Fill( phi, eta, w );
    rfp++;
    if(fUseFP) PushBackFlowTrack(fTPCevent,track->Pt(),phi,eta,w,track->GetID());
  }
  if(fUseFP) fTPCevent->SetNumberOfRPs(rfp);
  if(fDebug) printf("Inserted tracks in FlowEventTPC %d ==> %.1f\n",fTPCevent->NumberOfTracks(),qwa+qwb);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddMCParticleSpy(TList *me) {
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPro;
  tH1D = new TH1D("Pt",   "Pt",   100,0.0,20);  me->Add(tH1D);
  tH1D = new TH1D("Phi",  "Phi",  100,0,TMath::TwoPi()); me->Add(tH1D);
  tH1D = new TH1D("Eta",  "Eta",  100,-1,+1);   me->Add(tH1D);
  tH1D = new TH1D("Rad2", "Rad2", 1000,0,+100); me->Add(tH1D);
  tH2D = new TH2D("Dphi", "phi-MCEP;pt;dphi",100,0,20, 72,0,TMath::Pi()); me->Add(tH2D);
  tPro = new TProfile("Cos2dphi","Cos2dphi",100,0,20); me->Add(tPro);
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillMCParticleSpy(TString listName, AliAODMCParticle *p) {
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Pt" ))->Fill( p->Pt() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Eta" ))->Fill( p->Eta() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Phi" ))->Fill( p->Phi() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Rad2" ))->Fill( TMath::Sqrt( p->Xv()*p->Xv() +
												 p->Yv()*p->Yv() ) );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Dphi" ))->Fill( p->Pt(), GetMCDPHI(p->Phi()) );
  ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("Cos2dphi" ))->Fill( p->Pt(), TMath::Cos( 2*GetMCDPHI(p->Phi()) ), 1 );
  return;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::GetMCDPHI(Double_t phi) {
  Double_t dDPHI = phi - fMCEP;
  if( dDPHI < 0 ) dDPHI += TMath::TwoPi();
  if( dDPHI > TMath::Pi() ) dDPHI = TMath::TwoPi()-dDPHI;
  return dDPHI;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillMCParticleSpy(TString listName, TParticle *p) {
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Pt" ))->Fill( p->Pt() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Eta" ))->Fill( p->Eta() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Phi" ))->Fill( p->Phi() );
  ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("Rad2" ))->Fill( TMath::Sqrt( p->Vx()*p->Vx() +
												 p->Vy()*p->Vy() ) );
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddCandidatesSpy(TList *me,Bool_t res) {
  TH1D *tH1D;
  TH2D *tH2D;
  TProfile *tPro;
  TProfile2D *tPro2;
  tH2D = new TH2D("PhiEta",  "PhiEta;Phi;Eta", 100,0,TMath::TwoPi(),100,-1,+1); me->Add(tH2D);
  tH2D = new TH2D("PtRAP",   "PtRAP;Pt;Y",          100,0,20,100,-2.0,+2.0);  me->Add(tH2D);
  tH2D = new TH2D("PtDCA",   "PtDCA;Pt;DCA",        100,0,20,100,0,10);       me->Add(tH2D);
  tH2D = new TH2D("PtCTP",   "PtCTP;Pt;CTP",        100,0,20,100,-1,+1);      me->Add(tH2D);
  //tH2D = new TH2D("PtCTP",   "PtCTP;Pt;CTP",        100,0,10,100,0.90,+1);    me->Add(tH2D);
  tH2D = new TH2D("PtD0D0",  "PtD0D0;Pt;D0D0",      100,0,20,100,-5,+5);      me->Add(tH2D);
  tH2D = new TH2D("PtRad2",  "PtRad2;Pt;RadXY",     100,0,20,100,0,+50);      me->Add(tH2D);
  tH2D = new TH2D("PtDL",    "PtDL;Pt;DL",          100,0,20,100,0,+50);      me->Add(tH2D);
  tH2D = new TH2D("PtMASS",  "PtMASS;Pt;MASS", 100,0,20,fMassBins,fMinMass,fMaxMass); me->Add(tH2D);
  tH2D = new TH2D("APPOS",   "APPOS;alphaPOS;QtPOS",100,-2,+2,100,0,0.3);     me->Add(tH2D);
  tH2D = new TH2D("D0PD0N",  "D0PD0N;D0P;D0N",      200,-10,+10,200,-10,+10); me->Add(tH2D);
  tH2D = new TH2D("XPOSXNEG","XPOSXNEG;XPOS;XNEG",  200,-50,+50,200,-50,+50); me->Add(tH2D);

  if(fReadMC) {
    if(res) {
      tH1D = new TH1D("MCOrigin", "MCOrigin;Rad2",1000,0,50); me->Add(tH1D);
      tH2D = new TH2D("PHIRes","PHIRes;PHI;MC-DAT", 72,   0, TMath::TwoPi(),100,-0.12,+0.12); me->Add(tH2D);
      tH2D = new TH2D("ETARes","ETARes;ETA;MC-DAT", 16,-0.8,           +0.8,100,-0.2,+0.2);   me->Add(tH2D);
      tH2D = new TH2D("PTRes", "PTRes;Pt;MC-DAT",  100,   0,             20,100,-0.4,+0.4);   me->Add(tH2D);
      tH2D = new TH2D("RXYRes","RXYRes;RXY;MC-DAT",100,   0,             50,100,-4.0,+4.0);   me->Add(tH2D);
    }
    tH2D = new TH2D("PTDPHIMC","PtDPHIMC;Pt;PHI-MCEP",100,0,20,72,0,TMath::Pi()); me->Add(tH2D);
    tPro = new TProfile("Cos2dphiMC",  "Cos2dphiMC",100,0,20); me->Add(tPro);
    tPro = new TProfile("Cos2dphiMCPK","Cos2dphiMC PK",100,0,20); me->Add(tPro);
    tPro2=new TProfile2D("C2DPHIMCMASS","C2DPHIMCMASS",100,0,20,fMassBins,fMinMass,fMaxMass); me->Add(tPro2);
  }
  return;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillCandidateSpy(TString listName, Bool_t fillRes) {
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PhiEta"))->Fill( fDecayPhi, fDecayEta );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtRAP" ))->Fill( fDecayPt, fDecayRapidity );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtDCA" ))->Fill( fDecayPt, fDecayDCAdaughters );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtCTP" ))->Fill( fDecayPt, fDecayCosinePointingAngleXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtD0D0"))->Fill( fDecayPt, fDecayProductIPXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtRad2"))->Fill( fDecayPt, fDecayRadXY );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtDL"  ))->Fill( fDecayPt, fDecayDecayLength );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PtMASS"))->Fill( fDecayPt, fDecayMass );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("APPOS" ))->Fill( fDecayAlpha, fDecayQt );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("D0PD0N"))->Fill( fDecayIPpos, fDecayIPneg );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("XPOSXNEG"))->Fill( fDecayXpos, fDecayXneg );
  if(fReadMC) {
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTDPHIMC" ))->Fill( fDecayPt, GetMCDPHI( fDecayPhi ) );
    ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("Cos2dphiMC" ))->Fill( fDecayPt, TMath::Cos( 2*GetMCDPHI(fDecayPhi) ), 1 );
    if( fDecayMass>fMinMassX && fDecayMass<fMaxMassX ) {
      ((TProfile*)((TList*)fList->FindObject(listName.Data()))->FindObject("Cos2dphiMCPK" ))->Fill( fDecayPt, TMath::Cos( 2*GetMCDPHI(fDecayPhi) ), 1 );
    }
    ((TProfile2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("C2DPHIMCMASS" ))->Fill(fDecayPt,fDecayMass,TMath::Cos(2*GetMCDPHI(fDecayPhi)), 1 );
    if(fillRes) {
      ((TH1D*)((TList*)fList->FindObject(listName.Data()))->FindObject("MCOrigin"))->Fill( fDecayMatchOrigin );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PHIRes"))->Fill( fDecayPhi, fDecayMatchPhi-fDecayPhi );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("ETARes"))->Fill( fDecayEta, fDecayMatchEta-fDecayEta );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTRes"))->Fill( fDecayPt, fDecayMatchPt-fDecayPt );
      ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("RXYRes"))->Fill( fDecayRadXY, fDecayMatchRadXY-fDecayRadXY );
    }
  }
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptCandidate() {
  if(fDecayEta<fDecayMinEta) return kFALSE;
  if(fDecayEta>fDecayMaxEta) return kFALSE;
  if(fDecayPt<fDecayMinPt) return kFALSE;
  if(fDecayProductIPXY>fDecayMaxProductIPXY) return kFALSE;
  if(fDecayDCAdaughters>fDecayMaxDCAdaughters) return kFALSE;
  if(fDecayCosinePointingAngleXY<fDecayMinCosinePointingAngleXY) return kFALSE;
  if(fDecayRadXY<fDecayMinRadXY) return kFALSE;
  if(TMath::Abs(fDecayRapidity)>fDecayMaxRapidity) return kFALSE;
  if(fSpecie==0) {
    if(fDecayAPCutPie) {
      if(fDecayQt/TMath::Abs(fDecayAlpha)<fDecayMinQt) return kFALSE;
    } else {
      if(fDecayQt<fDecayMinQt) return kFALSE;
    }
    if(fDecayDecayLength>fDecayMaxDecayLength*2.6842) return kFALSE;
  } else {
    if(fDecayDecayLength>fDecayMaxDecayLength*7.89) return kFALSE;
  }
  if(fSpecie==1) if(fDecayAlpha>0) return kFALSE;
  if(fSpecie==2) if(fDecayAlpha<0) return kFALSE;
  return kTRUE;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddTracksSpy(TList *me,Bool_t res) {
  TH2D *tH2D;
  tH2D = new TH2D("PHIETA",       "PHIETA;PHI;ETA",               100,0,TMath::TwoPi(),100,-2,2); me->Add(tH2D);
  tH2D = new TH2D("IPXYIPZ",      "IPXYIPZ;IPXY;IPZ",             1000,-20,+20,1000,-20,+20); me->Add(tH2D);
  tH2D = new TH2D("PTTPCNCLS",    "PTTPCNCLS;PT;NCLS",            100,0,20,170,0,170);  me->Add(tH2D);
  tH2D = new TH2D("PTITSLAY",     "PTITSLAY;PT;ITSLAYER",         100,0,20,6,-0.5,+5.5);me->Add(tH2D);
  tH2D = new TH2D("PTITSTPCrefit","PTITSTPCrefit;PT",             100,0,20,2,-0.5,+1.5);me->Add(tH2D);
  tH2D->GetYaxis()->SetBinLabel(1,"ITS refit");
  tH2D->GetYaxis()->SetBinLabel(2,"TPC refit");

  tH2D = new TH2D("POSTPCNCLCHI2","POSTPCNCLCHI2;NCLS;CHI2/NCLS", 170,0,170,100,0,8);   me->Add(tH2D);
  tH2D = new TH2D("POSTPCNFCLNXR","POSTPCNFCLNXR;NFCLS;NXR",      170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("POSTPCNCLNFCL","POSTPCNCLNFCL;NCLS;NFCLS",     170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("POSTPCNCLNSCL","POSTPCNCLNSCL;NCLS;NSCLS",     170,0,170,170,0,170); me->Add(tH2D);

  tH2D = new TH2D("NEGTPCNCLCHI2","NEGTPCNCLCHI2;NCLS;CHI2/NCLS", 170,0,170,100,0,8);   me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNFCLNXR","NEGTPCNFCLNXR;NFCLS;NXR",      170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNCLNFCL","NEGTPCNCLNFCL;NCLS;NFCLS",     170,0,170,170,0,170); me->Add(tH2D);
  tH2D = new TH2D("NEGTPCNCLNSCL","NEGTPCNCLNSCL;NCLS;NSCLS",     170,0,170,170,0,170); me->Add(tH2D);

  if(res) {
    tH2D = new TH2D("PHIRes", "PHIRes;PHI;MC-DAT",   72,   0, TMath::TwoPi(),100,-0.12,+0.12); me->Add(tH2D);
    tH2D = new TH2D("ETARes", "ETARes;ETA;MC-DAT",   16,-0.8,           +0.8,100,-0.2,+0.2);   me->Add(tH2D);
    tH2D = new TH2D("PTRes",  "PTRes;Pt;MC-DAT",    100,   0,             20,100,-0.4,+0.4);   me->Add(tH2D);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::FillTrackSpy(TString listName,Bool_t res) {
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PHIETA" ))->Fill( fDaughterPhi, fDaughterEta );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "IPXYIPZ" ))->Fill( fDaughterImpactParameterXY, fDaughterImpactParameterZ );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTTPCNCLS" ))->Fill( fDaughterPt, fDaughterNClsTPC );
  if( TESTBIT(fDaughterITScm,0) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 0 );
  if( TESTBIT(fDaughterITScm,1) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 1 );
  if( TESTBIT(fDaughterITScm,2) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 2 );
  if( TESTBIT(fDaughterITScm,3) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 3 );
  if( TESTBIT(fDaughterITScm,4) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 4 );
  if( TESTBIT(fDaughterITScm,5) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSLAY" ))->Fill( fDaughterPt, 5 );
  if( (fDaughterStatus&AliESDtrack::kITSrefit) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSTPCrefit" ))->Fill( fDaughterPt, 0 );
  if( (fDaughterStatus&AliESDtrack::kTPCrefit) ) ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( "PTITSTPCrefit" ))->Fill( fDaughterPt, 1 );

  TString ch="NEG";
  if(fDaughterCharge>0) ch="POS";
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLCHI2",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterChi2PerNClsTPC );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNFCLNXR",ch.Data()) ))->Fill( fDaughterNFClsTPC, fDaughterXRows );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLNFCL",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterNFClsTPC );
  ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject( Form("%sTPCNCLNSCL",ch.Data()) ))->Fill( fDaughterNClsTPC, fDaughterNSClsTPC );

  if(res) {
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PHIRes"))->Fill( fDaughterPhi, fDaughterMatchPhi-fDaughterPhi );
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("ETARes"))->Fill( fDaughterEta, fDaughterMatchEta-fDaughterEta );
    ((TH2D*)((TList*)fList->FindObject(listName.Data()))->FindObject("PTRes"))->Fill( fDaughterPt, fDaughterMatchPt-fDaughterPt );
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::LoadTrack(AliESDtrack *myTrack, Double_t aodChi2NDF) {
  fDaughterCharge = myTrack->Charge();
  fDaughterXRows = myTrack->GetTPCCrossedRows();
  fDaughterNFClsTPC = myTrack->GetTPCNclsF();
  fDaughterNSClsTPC = myTrack->GetTPCnclsS();
  fDaughterNClsTPC = myTrack->GetTPCclusters(0);
  if(fReadESD) {
    if(fDaughterNClsTPC>0) fDaughterChi2PerNClsTPC = myTrack->GetTPCchi2()/fDaughterNClsTPC;
  } else {
    fDaughterChi2PerNClsTPC = aodChi2NDF;
  }
  myTrack->GetImpactParameters(fDaughterImpactParameterXY,fDaughterImpactParameterZ);
  fDaughterStatus = myTrack->GetStatus();
  fDaughterITScm = myTrack->GetITSClusterMap();
  fDaughterPhi = myTrack->Phi();
  fDaughterEta = myTrack->Eta();
  fDaughterPt = myTrack->Pt();
  fDaughterKinkIndex = myTrack->GetKinkIndex(0);
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::AcceptDaughter() {
  if(fDaughterKinkIndex>0) return kFALSE;
  if( (fDaughterStatus&AliESDtrack::kTPCrefit)==0 ) return kFALSE;
  if(fDaughterNFClsTPC<1) return kFALSE;
  if(fDaughterPt<fDaughterMinPt) return kFALSE;
  if(fDaughterEta<fDaughterMinEta) return kFALSE;
  if(fDaughterEta>fDaughterMaxEta) return kFALSE;
  if(fDaughterNClsTPC<fDaughterMinNClsTPC) return kFALSE;
  if(fDaughterXRows<fDaughterMinXRows) return kFALSE;
  if(fDaughterChi2PerNClsTPC>fDaughterMaxChi2PerNClsTPC) return kFALSE;
  if(TMath::Abs(fDaughterImpactParameterXY)<fDaughterMinImpactParameterXY) return kFALSE;
  if(fDaughterXRows<fDaughterMinXRowsOverNClsFTPC*fDaughterNFClsTPC) return kFALSE;
  for(Int_t lay=0; lay!=6; ++lay)
    if(fDaughterITSConfig[lay]>-0.5) {
      if(fDaughterITSConfig[lay]) {
	if(!TESTBIT(fDaughterITScm,lay)) return kFALSE;
      } else {
	if(TESTBIT(fDaughterITScm,lay)) return kFALSE;
      }
    }
  return kTRUE;
}
//=======================================================================
Double_t AliAnalysisTaskFlowStrange::GetWDist(const AliVVertex* v0, const AliVVertex* v1) {
  // calculate sqrt of weighted distance to other vertex
  if (!v0 || !v1) {
    printf("One of vertices is not valid\n");
    return 0;
  }
  static TMatrixDSym vVb(3);
  double dist = -1;
  double dx = v0->GetX()-v1->GetX();
  double dy = v0->GetY()-v1->GetY();
  double dz = v0->GetZ()-v1->GetZ();
  double cov0[6],cov1[6];
  v0->GetCovarianceMatrix(cov0);
  v1->GetCovarianceMatrix(cov1);
  vVb(0,0) = cov0[0]+cov1[0];
  vVb(1,1) = cov0[2]+cov1[2];
  vVb(2,2) = cov0[5]+cov1[5];
  vVb(1,0) = vVb(0,1) = cov0[1]+cov1[1];
  vVb(0,2) = vVb(1,2) = vVb(2,0) = vVb(2,1) = 0.;
  vVb.InvertFast();
  if (!vVb.IsValid()) {printf("Singular Matrix\n"); return dist;}
  dist = vVb(0,0)*dx*dx + vVb(1,1)*dy*dy + vVb(2,2)*dz*dz
    +    2*vVb(0,1)*dx*dy + 2*vVb(0,2)*dx*dz + 2*vVb(1,2)*dy*dz;
  return dist>0 ? TMath::Sqrt(dist) : -1;
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::plpMV(const AliVEvent *event) {
  // check for multi-vertexer pile-up
  const AliAODEvent *aod = (const AliAODEvent*)event;
  const AliESDEvent *esd = (const AliESDEvent*)event;
  //
  const int    kMinPlpContrib = 5;
  const double kMaxPlpChi2 = 5.0;
  const double kMinWDist = 15;
  //
  if (!aod && !esd) {
    printf("Event is neither of AOD nor ESD\n");
    exit(1);
  }
  //
  const AliVVertex* vtPrm = 0;
  const AliVVertex* vtPlp = 0;
  int nPlp = 0;
  //
  if (aod) {
    if ( !(nPlp=aod->GetNumberOfPileupVerticesTracks()) ) return kFALSE;
    vtPrm = aod->GetPrimaryVertex();
    if (vtPrm == aod->GetPrimaryVertexSPD()) return kTRUE; // there are pile-up vertices but no primary
  }
  else {
    if ( !(nPlp=esd->GetNumberOfPileupVerticesTracks())) return kFALSE;
    vtPrm = esd->GetPrimaryVertexTracks();
    if (((AliESDVertex*)vtPrm)->GetStatus()!=1) return kTRUE; // there are pile-up vertices but no primary
  }
  //int bcPrim = vtPrm->GetBC();
  //
  for (int ipl=0;ipl<nPlp;ipl++) {
    vtPlp = aod ? (const AliVVertex*)aod->GetPileupVertexTracks(ipl) : (const AliVVertex*)esd->GetPileupVertexTracks(ipl);
    //
    if (vtPlp->GetNContributors() < kMinPlpContrib) continue;
    if (vtPlp->GetChi2perNDF() > kMaxPlpChi2) continue;
    //  int bcPlp = vtPlp->GetBC();
    //  if (bcPlp!=AliVTrack::kTOFBCNA && TMath::Abs(bcPlp-bcPrim)>2) return kTRUE; // pile-up from other BC
    //
    double wDst = GetWDist(vtPrm,vtPlp);
    if (wDst<kMinWDist) continue;
    //
    return kTRUE; // pile-up: well separated vertices
  }
  //
  return kFALSE;
  //
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeFilterBits() {
  //FilterBit 1
  fFB1 = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  //FilterBit1024
  fFB1024 = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,1);
  fFB1024->SetMinNCrossedRowsTPC(120);
  fFB1024->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
  fFB1024->SetMaxChi2PerClusterITS(36);
  fFB1024->SetMaxFractionSharedTPCClusters(0.4);
  fFB1024->SetMaxChi2TPCConstrainedGlobal(36);
  fFB1024->SetEtaRange(-0.9,0.9);
  fFB1024->SetPtRange(0.15, 1e10);
}
//=======================================================================
Bool_t AliAnalysisTaskFlowStrange::PassesFilterBit(AliESDtrack *track) {
  Bool_t ret=kFALSE;
  switch(fRFPFilterBit) {
    case(1024):
      ret = fFB1024->AcceptTrack(track);
      break;
    default:
      ret = fFB1->AcceptTrack(track);
  }
  return ret;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::LoadVZEROResponse() {
  if(fVZEResponse) {
    TString run = fVZEResponse->GetTitle();
    if( run.Atoi() == fRunNumber ) return;
    fVZEResponse = NULL;
  }
  //==>loading
  fVZEResponse = dynamic_cast<TH2D*> (fVZEload->FindObject( Form("%d",fRunNumber) ));
  if(fVZEResponse) {
    printf("New VZE calibration: run %d || %s -> Entries %.0f\n",fRunNumber, fVZEResponse->GetTitle(),fVZEResponse->GetEntries());
  } else {
    printf("VZE calibration: request but not found!!!\n");
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddVZEQA() {
  fVZEQA = new TList();
  fVZEQA->SetName( Form("VZEQA%d",fRunNumber) );
  fVZEQA->SetOwner();
  if(fQAlevel>0) {
    TProfile2D *prof = new TProfile2D("LINP","LINP;VZEcell;VZEmult;SPDtrkl", 64,0,64,500,0,700,0,10000); fVZEQA->Add( prof );
    prof = new TProfile2D("MULP","MULP;VZEcell;CENTR;VZEmult", 64,0,64,100,0,100,0,10000); fVZEQA->Add( prof );
    TH3D *tH3D = new TH3D("EQU","EQU;VZEeqmult;VZEmult",100,0,700,100,0,700,64,0,64); fVZEQA->Add( tH3D );
  }
  fList->Add(fVZEQA);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SaveVZEROQA() {
  AliAODEvent *event = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!event) return;
  AliVVZERO *vzero = event->GetVZEROData();
  AliAODTracklets *tracklets = event->GetTracklets();
  if(!vzero) return;
  if(!tracklets) return;
  Double_t mult, eqmult;
  Int_t trkl;
  for(int id=0; id!=64; ++id) {
    trkl = tracklets->GetNumberOfTracklets();
    mult = vzero->GetMultiplicity(id);
    eqmult = event->GetVZEROEqMultiplicity(id);
    ((TProfile2D*) fVZEQA->FindObject( "LINP" ))->Fill(id,mult,trkl,1);
    ((TProfile2D*) fVZEQA->FindObject( "MULP" ))->Fill(id,fThisCent,mult,1);
    ((TH3D*) fVZEQA->FindObject("EQU"))->Fill(eqmult,mult,id);
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::AddVZEROResponse() {
  fVZEResponse = NULL;
  AliVEvent *event = InputEvent();
  if(!event) return;
  Int_t thisrun = event->GetRunNumber();
  fVZEResponse = new TH2D( Form("%d",thisrun), Form("%d;cell;CC",thisrun), 64,0,64, 100, 0, 100);
  fList->Add(fVZEResponse);
}
//=======================================================================
void AliAnalysisTaskFlowStrange::SaveVZEROResponse() {
  if(!fVZEResponse) return;
  AliVEvent *event = InputEvent();
  if(!event) return;
  Double_t w;
  for(int id=0; id!=64; ++id) {
    w = event->GetVZEROEqMultiplicity(id);
    fVZEResponse->Fill(id,fThisCent,w);
  }
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::RefMultTPC() {
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!ev) return -1;
  Int_t found = 0;
  for(int i=0; i!=ev->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = ev->GetTrack( i );
    if(!t) continue;
    if( !t->TestFilterBit(1) ) continue;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) continue;
    if( t->Pt()<0.2 || t->Pt()>5.0 ) continue;
    if( t->GetTPCNcls()<70 ) continue;
    if( t->GetTPCsignal()<10.0 ) continue;
    if( t->Chi2perNDF()<0.2 ) continue;    
    ++found;
  }
  return found;
}
//=======================================================================
Int_t AliAnalysisTaskFlowStrange::RefMultGlobal() {
  AliAODEvent *ev = dynamic_cast<AliAODEvent*> (InputEvent());
  if(!ev) return -1;
  Int_t found = 0;
  for(int i=0; i!=ev->GetNumberOfTracks(); ++i) {
    AliAODTrack *t = ev->GetTrack( i );
    if(!t) continue;
    if( !t->TestFilterBit(16) ) continue;
    if( t->Eta()<-0.8 || t->Eta()>+0.8 ) continue;
    if( t->Pt()<0.2 || t->Pt()>5.0 ) continue;
    if( t->GetTPCNcls()<70 ) continue;
    if( t->GetTPCsignal()<10.0 ) continue;
    if( t->Chi2perNDF()<0.1 ) continue;    
    Double_t b[3], bcov[3];
    if( !t->PropagateToDCA(ev->GetPrimaryVertex(),ev->GetMagneticField(),100,b,bcov) ) continue;
    if( b[0]>+0.3 || b[0]<-0.3 || b[1]>+0.3 || b[1]<-0.3) continue;
    ++found;
  }
  return found;
}
//=======================================================================
void AliAnalysisTaskFlowStrange::MakeDHcorr() {
  // Adds DH corr
  for(int iCand=0; iCand!=fCandidates->GetEntriesFast(); ++iCand ) {
    AliFlowCandidateTrack *cand = static_cast<AliFlowCandidateTrack*>(fCandidates->At(iCand));
    if(!cand) continue;
    for(int iRPs=0; iRPs!=fTPCevent->NumberOfTracks(); ++iRPs ) {
      AliFlowTrack *iRP = static_cast<AliFlowTrack*>(fTPCevent->GetTrack( iRPs ));
      if(!iRP) continue;
      if(!iRP->InRPSelection()) continue;
      if(cand->GetID() == iRP->GetID()) continue; //avoid autocorr
      for(int iDau=0; iDau!=cand->GetNDaughters(); ++iDau) //if it is a decay
        if(cand->GetIDDaughter(iDau) == iRP->GetID()) continue;
      //corr
      Double_t dDPHI = iRP->Phi() - cand->Phi();
      Double_t dDETA = iRP->Eta() - cand->Eta();
      Double_t dDPT  = iRP->Pt()  - cand->Pt();
      ((TH3D*)((TList*)fList->FindObject("DHCORR"))->FindObject("DPHI"))->Fill( dDPT, dDPHI, dDETA );
      //end of corr
    }
  }
}
//=======================================================================
void AliAnalysisTaskFlowStrange::ResetContainers() {
  fTPCevent->ClearFast();
  fVZEevent->ClearFast();
}
