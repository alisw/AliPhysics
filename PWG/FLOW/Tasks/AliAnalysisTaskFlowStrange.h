/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez (cperez@cern.ch)
//          Pawel Debski (pdebski@cern.ch)
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICExperiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowStrange_H
#define AliAnalysisTaskFlowStrange_H

#include "AliAnalysisTaskSE.h"

class TList;
class TH2D;
class TObjArray;
class TClonesArray;
class AliAODMCParticle;
class AliESDtrackCuts;
class AliFlowEventCuts;
class AliPIDResponse;
class AliESDEvent;
class AliAODEvent;
class AliAODv0;
class AliESDv0;
class AliVVertex;
class AliFlowBayesianPID;
class AliAODVertex;

class AliAnalysisTaskFlowStrange : public AliAnalysisTaskSE {
 public:
  enum Especie {kKZE=0,kLDA=1,kLDABAR=2,kLDAALL=3,kCHARGED=90,kPION=91,kKAON=92,kPROTON=93};
  AliAnalysisTaskFlowStrange();
  AliAnalysisTaskFlowStrange(const Char_t *name);
  virtual ~AliAnalysisTaskFlowStrange();
  virtual void UserCreateOutputObjects();
  virtual void Exec(Option_t*);
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();

  void MyUserExec(Option_t *);
  void SetDebug(Int_t val=1) {fDebug = val;}
  void SetQAlevel(Int_t qa) {fQAlevel = qa;}

  void SetMass(Int_t n, Double_t m, Double_t M) {fMassBins=n;fMinMass=m;fMaxMass=M;}

  void SetpA() {fRunOnpA = kTRUE;  fRunOnpp = kFALSE;}
  void Setpp() {fRunOnpA = kFALSE; fRunOnpp = kTRUE; }
  void SetK0L0(Int_t specie) {fSpecie=specie;}
  void SetOnline(Bool_t val) {fOnline=val;}
  void SetHomemade(Bool_t val) {fHomemade=val;}
  void SetExcludeTPCEdges(Bool_t value) {fExcludeTPCEdges=value;}
  void SetCentralityRange(TString val, Int_t m, Int_t M) {fCentMethod=val; fCentPerMin=m; fCentPerMax=M;}
  void SetReadESD(Bool_t val) {fReadESD=val;}
  void SetReadMC(Bool_t val) {fReadMC=val;}
  void SetAvoidExec(Bool_t val) {fAvoidExec=val;}
  void SetSkipSelection(Bool_t val) {fSkipSelection=val;}
  void SetSkipFlow(Bool_t val) {fSkipFlow=val;}
  void SetUseFlowPackage(Bool_t val) {fUseFP=val;}
  void SetExtraEventRejection(Bool_t val) {fExtraEventRejection=val;}

  void SetWhichPsi(Int_t val) {fWhichPsi=val;}
  void SetStoreVZEResponse(Bool_t val) {fVZEsave=val;}
  void LoadVZEResponse(TList *val) {fVZEload=val;}
  
  void SetRFPFilterBit(Int_t val) {fRFPFilterBit=val;}
  void SetRFPMinPt(Double_t val) {fRFPminPt=val;}
  void SetRFPMaxPt(Double_t val) {fRFPmaxPt=val;}
  void SetRFPMinEta(Double_t val) {fRFPminEta=val;}
  void SetRFPMaxEta(Double_t val) {fRFPmaxEta=val;}
  void SetRFPTPCSignal(Double_t val) {fRFPTPCsignal=val;}
  void SetRFPMaxIPxy(Double_t val) {fRFPmaxIPxy=val;}
  void SetRFPMaxIPz(Double_t val) {fRFPmaxIPz=val;}
  void SetRFPMinTPCCls(Int_t val) {fRFPTPCncls=val;}

  void SetDauMinNClsTPC(Int_t val) {fDaughterMinNClsTPC=val;}
  void SetDauMinXRows(Int_t val) {fDaughterMinXRows=val;}
  void SetDauMaxChi2PerNClsTPC(Double_t val) {fDaughterMaxChi2PerNClsTPC=val;}
  void SetDauMinXRowsOverNClsFTPC(Double_t val) {fDaughterMinXRowsOverNClsFTPC=val;}
  void SetDauMinEta(Double_t val) {fDaughterMinEta=val;}
  void SetDauMaxEta(Double_t val) {fDaughterMaxEta=val;}
  void SetDauMinPt(Double_t val) {fDaughterMinPt=val;}
  void SetDauMinImpactParameterXY(Double_t val) {fDaughterMinImpactParameterXY=val;}
  void SetDauMaxNSigmaPID(Double_t val) {fDaughterMaxNSigmaPID=val;}

  void SetMaxRapidity(Double_t val) {fDecayMaxRapidity=val;}
  void SetMinEta(Double_t val) {fDecayMinEta=val;}
  void SetMaxEta(Double_t val) {fDecayMaxEta=val;}
  void SetMinPt(Double_t val) {fDecayMinPt=val;}
  void SetMaxDCAdaughters(Double_t val) {fDecayMaxDCAdaughters=val;}
  void SetMinCosinePointingAngleXY(Double_t val) {fDecayMinCosinePointingAngleXY=val;}
  void SetMinQt(Double_t val, Bool_t val2=kTRUE) {fDecayMinQt=val; fDecayAPCutPie=val2;}
  void SetMinRadXY(Double_t val) {fDecayMinRadXY=val;}
  void SetMaxDecayLength(Double_t val) {fDecayMaxDecayLength=val;}
  void SetMaxProductIPXY(Double_t val) {fDecayMaxProductIPXY=val;}

 private:
  AliAnalysisTaskFlowStrange(const AliAnalysisTaskFlowStrange& analysisTask);
  AliAnalysisTaskFlowStrange& operator=(const AliAnalysisTaskFlowStrange& analysisTask);
  void AddQAEvents();
  void AddQACandidates();

  void AddEventSpy();
  Bool_t AcceptAAEvent(AliESDEvent *tESD);
  Bool_t AcceptAAEvent(AliAODEvent *tAOD);
  Bool_t AcceptPPEvent(AliAODEvent *tAOD);
  Bool_t AcceptPAEvent(AliAODEvent *tAOD);
  Int_t GetReferenceMultiplicity();

  void ReadStack(TClonesArray* mcArray);
  void ReadFromESD(AliESDEvent *tESD);
  void ReadFromAODv0(AliAODEvent *tAOD);

  void ChargeParticles(AliAODEvent *tAOD);

  void ComputePsi2(AliVEvent *event);
  void AddMakeQSpy();
  void MakeQVZE(AliVEvent *event,Double_t &qxa,Double_t &qya,Double_t &qwa,Double_t &qxb,Double_t &qyb,Double_t &qwb);
  void MakeQTPC(AliESDEvent *event,Double_t &qxa,Double_t &qya,Double_t &qwa,Double_t &qxb,Double_t &qyb,Double_t &qwb);
  void MakeQTPC(AliAODEvent *event,Double_t &qxa,Double_t &qya,Double_t &qwa,Double_t &qxb,Double_t &qyb,Double_t &qwb);
  void AddTPCRFPSpy(TList *val);
  Bool_t PassesRFPTPCCuts(AliESDtrack *myTrack, Double_t aodChi2NDF=0, Float_t aodipxy=0, Float_t aodipz=0);
  void MakeQVectors();

  void AddCandidates();
  void ReadEventPlanesFromAOD(AliAODEvent *tAOD);

  Double_t CosThetaPointXY(AliESDv0 *me, const AliVVertex *vtx);
  Double_t CosThetaPointXY(AliAODv0 *me, const AliVVertex *vtx);
  Double_t DecayLengthXY(AliESDv0 *me, const AliVVertex *vtx);
  Double_t DecayLengthXY(AliAODv0 *me, const AliVVertex *vtx);
  Double_t DecayLength(AliESDv0 *me, const AliVVertex *vtx);
  Double_t DecayLength(AliAODv0 *me, const AliVVertex *vtx);

  void AddMCParticleSpy(TList *val);
  void FillMCParticleSpy(TString listName, AliAODMCParticle *par);
  void FillMCParticleSpy(TString listName, TParticle *par);

  void AddCandidatesSpy(TList *val);
  void FillCandidateSpy(TString listName);

  void AddTracksSpy(TList *val);
  void FillTrackSpy(TString listName);

  void MakeFilterBits();
  Bool_t PassesFilterBit(AliESDtrack *me);

  void LoadTrack(AliESDtrack *myTrack, Double_t aodChi2NDF=0);
  Bool_t AcceptDaughter();
  Bool_t AcceptCandidate();
  Bool_t PassesPIDCuts(AliESDtrack *myTrack, AliPID::EParticleType pid=AliPID::kProton);

  Bool_t IsAtTPCEdge(Double_t phi,Double_t pt,Int_t charge,Double_t b);

  void MakeTrack();
  void PushBackFlowTrack(AliFlowEvent *event, Double_t pt, Double_t phi, Double_t eta, Double_t we, Int_t id);

  Double_t GetWDist(const AliVVertex* v0, const AliVVertex* v1);
  Bool_t plpMV(const AliVEvent *event);

  void LoadVZEROResponse();
  void AddVZEROResponse();
  void SaveVZEROResponse();
  void AddVZEQA();
  void SaveVZEROQA();

  Int_t RefMultTPC();
  Int_t RefMultGlobal();

  AliPIDResponse *fPIDResponse; //! PID response object
  AliESDtrackCuts *fFB1;        // filterbit cut equivalent
  AliESDtrackCuts *fFB1024;     // filterbit cut equivalent
  AliFlowEvent   *fTPCevent;    // flow event (needed here due to ev selection)
  AliFlowEvent   *fVZEevent;    // flow event (needed here due to ev selection)
  TObjArray      *fCandidates;  // array of selected candidates
  TList          *fList;        // stores the final list of output histograms

  Int_t fDebug;   // debug level
  Int_t fQAlevel; // QA plots

  Bool_t fReadESD;       // move back to ESD
  Bool_t fReadMC;        // read MC files
  Bool_t fAvoidExec;     // avoids Exec
  Bool_t fSkipSelection; // skip decay finder
  Bool_t fSkipFlow;      // skip flow-wise code
  Bool_t fUseFP;         // flow package?
  Bool_t fRunOnpA;       // make task compatible with pA event selection
  Bool_t fRunOnpp;       // make task compatible with pp event selection
  Bool_t fExtraEventRejection; // to reject pile up
  TString  fCentMethod; // CC
  Int_t    fCentPerMin; // CC
  Int_t    fCentPerMax; // CC
  Double_t fThisCent;   // CC

  Bool_t fExcludeTPCEdges; // exclude TPC edges from single track selection

  Int_t  fSpecie;   // K0=>0 L0=>1
  Bool_t fOnline;   // change into online v0 finder
  Bool_t fHomemade; // homemade v0 finder

  Int_t fWhichPsi;  // detector for Psi2
  Bool_t  fVZEsave; // make vze response
  TList  *fVZEload; // adress to calibration file
  TH2D   *fVZEResponse; // vze response vs centrality class
  TList  *fVZEQA;   // adress to qalist
  Double_t fPsi2;   // best estimation of Psi2

  Int_t    fMassBins; // opens
  Double_t fMinMass;  // mass
  Double_t fMaxMass;  // window

  Int_t fRFPFilterBit;    // RFP TPC
  Double_t fRFPminPt;     // RFP TPC
  Double_t fRFPmaxPt;     // RFP TPC
  Double_t fRFPminEta;    // RFP TPC
  Double_t fRFPmaxEta;    // RFP TPC
  Double_t fRFPTPCsignal; // RFP TPC
  Double_t fRFPmaxIPxy;   // RFP TPC
  Double_t fRFPmaxIPz;    // RFP TPC
  Int_t fRFPTPCncls;      // RFP TPC

  Double_t fDecayMass;                  // DECAY
  Double_t fDecayPhi;                   // DECAY
  Double_t fDecayEta;                   // DECAY
  Double_t fDecayPt;                    // DECAY
  Double_t fDecayDCAdaughters;          // DECAY
  Double_t fDecayCosinePointingAngleXY; // DECAY
  Double_t fDecayRadXY;                 // DECAY
  Double_t fDecayDecayLength;           // DECAY
  Double_t fDecayQt;                    // DECAY
  Double_t fDecayAlpha;                 // DECAY
  Double_t fDecayRapidity;              // DECAY
  Double_t fDecayProductIPXY;           // DECAY
  Int_t    fDecayIDneg;                 // DECAY
  Int_t    fDecayIDpos;                 // DECAY

  Double_t fDecayMinEta;                   // DECAY CUTS
  Double_t fDecayMaxEta;                   // DECAY CUTS
  Double_t fDecayMinPt;                    // DECAY CUTS
  Double_t fDecayMaxDCAdaughters;          // DECAY CUTS
  Double_t fDecayMinCosinePointingAngleXY; // DECAY CUTS
  Double_t fDecayMinQt;                    // DECAY CUTS
  Bool_t   fDecayAPCutPie;                 // DECAY CUTS
  Double_t fDecayMinRadXY;                 // DECAY CUTS
  Double_t fDecayMaxDecayLength;           // DECAY CUTS
  Double_t fDecayMaxProductIPXY;           // DECAY CUTS
  Double_t fDecayMaxRapidity;              // DECAY CUTS

  Double_t fDaughterPhi;               // DAUGHTER
  Double_t fDaughterEta;               // DAUGHTER
  Double_t fDaughterPt;                // DAUGHTER
  Int_t    fDaughterNClsTPC;           // DAUGHTER
  Int_t    fDaughterCharge;            // DAUGHTER
  Int_t    fDaughterNFClsTPC;          // DAUGHTER
  Int_t    fDaughterNSClsTPC;          // DAUGHTER
  Double_t fDaughterChi2PerNClsTPC;    // DAUGHTER
  Double_t fDaughterXRows;             // DAUGHTER
  Float_t  fDaughterImpactParameterXY; // DAUGHTER
  Float_t  fDaughterImpactParameterZ;  // DAUGHTER
  UInt_t   fDaughterStatus;            // DAUGHTER
  Double_t fDaughterNSigmaPID;         // DAUGHTER
  Int_t    fDaughterKinkIndex;         // DAUGHTER

  Double_t fDaughterMinEta;               // DAUGHTER CUTS
  Double_t fDaughterMaxEta;               // DAUGHTER CUTS
  Double_t fDaughterMinPt;                // DAUGHTER CUTS
  Int_t    fDaughterMinNClsTPC;           // DAUGHTER CUTS
  Int_t    fDaughterMinXRows;             // DAUGHTER CUTS
  Double_t fDaughterMaxChi2PerNClsTPC;    // DAUGHTER CUTS
  Double_t fDaughterMinXRowsOverNClsFTPC; // DAUGHTER CUTS
  Double_t fDaughterMinImpactParameterXY; // DAUGHTER CUTS
  Double_t fDaughterMaxNSigmaPID;         // DAUGHTER CUTS

  ClassDef(AliAnalysisTaskFlowStrange, 5);
};
#endif
