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

#ifndef AliAnalysisTaskFlowStrangee_H
#define AliAnalysisTaskFlowStrangee_H

#include "AliAnalysisTaskSE.h"

class TList;
class TH2D;
class TObjArray;
class TClonesArray;
class AliAODMCParticle;
class AliESDtrackCuts;
class AliFlowEventCuts;
//class AliFlowTrackCuts;
class AliPIDResponse;
class AliESDEvent;
class AliAODEvent;
class AliAODv0;
class AliESDv0;
class AliVVertex;
class AliFlowBayesianPID;
class AliAODVertex;

class AliAnalysisTaskFlowStrangee : public AliAnalysisTaskSE {
 public:
  enum Especie {kKZE=0,kLDA=1,kLDABAR=2,kLDAALL=3,kCHARGED=90,kPION=91,kKAON=92,kPROTON=93};
  enum Econfig {kSpecie=1,kHarmonic,kReadMC,kSkipSelection};
  AliAnalysisTaskFlowStrangee();
  AliAnalysisTaskFlowStrangee(const Char_t *name);
  virtual ~AliAnalysisTaskFlowStrangee();
  virtual void UserCreateOutputObjects();
  virtual void Exec(Option_t*);
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  virtual void MyUserExec(Option_t *);
  virtual void MyUserCreateOutputObjects();
  virtual void MyPrintConfig();
  virtual void PrintConfig();

  void SetHarmonic(Int_t val) {fHarmonic= val;}

  void SetOutputList(TList *lst) {fList=lst;}
  TList* GetOutputList() {return fList;}
  TList* RunTerminateAgain(TList *lst);

  void SetDebug(Int_t val=1) {fDebug = val;}
  void SetQAlevel(Int_t qa) {fQAlevel = qa;}

  void SetpA() {fRunOnpA = kTRUE;  fRunOnpp = kFALSE;}
  void Setpp() {fRunOnpA = kFALSE; fRunOnpp = kTRUE; }
  void SetReadESD(Bool_t val) {fReadESD=val;}
  void SetReadMC(Bool_t val) {fReadMC=val;}

  void SetAvoidExec(Bool_t val) {fAvoidExec=val;}
  void SetVertexZcut(Double_t val) {fVertexZcut=val;}
  void SetSkipCentralitySelection(Bool_t val) {fSkipCentralitySelection=val;}
  void SetCentralityRange(TString val, Int_t m, Int_t M) {fCentMethod=val; fCentPerMin=m; fCentPerMax=M;}
  void SetExtraEventRejection(Bool_t val) {fExtraEventRejection=val;}
  void SetSkipTerminate(Bool_t val) {fSkipTerminate=val;}

  void SetAddPiToMCReactionPlane(Bool_t val) {fAddPiToMCReactionPlane=val;}
  void SetUseFlowPackage(Bool_t val) {fUseFP=val;}
  void SetWhichPsi(Int_t val) {fWhichPsi=val;}
  void SetStoreVZEResponse(Bool_t val) {fVZEsave=val;}
  void LoadVZEResponse(TList *val, Bool_t val2=kFALSE, Bool_t val3=kTRUE) {fVZEload=val;fVZEmb=val2;fVZEByDisk=val3;}
  void SetRFPFilterBit(Int_t val) {fRFPFilterBit=val;}
  void SetRFPMinPt(Double_t val) {fRFPminPt=val;}
  void SetRFPMaxPt(Double_t val) {fRFPmaxPt=val;}
  void SetRFPAMinEta(Double_t val) {fRFPAminEta=val;}
  void SetRFPAMaxEta(Double_t val) {fRFPAmaxEta=val;}
  void SetRFPCMinEta(Double_t val) {fRFPCminEta=val;}
  void SetRFPCMaxEta(Double_t val) {fRFPCmaxEta=val;}
  void SetRFPTPCSignal(Double_t val) {fRFPTPCsignal=val;}
  void SetRFPMaxIPxy(Double_t val) {fRFPmaxIPxy=val;}
  void SetRFPMaxIPz(Double_t val) {fRFPmaxIPz=val;}
  void SetRFPMinTPCCls(Int_t val) {fRFPTPCncls=val;}
  void SetRFPVZERingRange(Int_t val1, Int_t val2, Int_t val3, Int_t val4)
    {fVZECa=val1;fVZECb=val2;fVZEAa=val3;fVZEAb=val4;}
  //void SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }

  void SetSkipSelection(Bool_t val) {fSkipSelection=val;}
  void SetSkipVn(Bool_t val) {fSkipVn=val;}
  void SetPostMatched(Int_t val) {fPostMatched=val;}
  void SetK0L0(Int_t specie) {fSpecie=specie;}
  void SetMass(Int_t n, Double_t m, Double_t M) {fMassBins=n;fMinMass=m;fMaxMass=M;}
  void SetPtEdges(Int_t n, Double_t *p);
  void SetOnline(Bool_t val) {fOnline=val;}
  void SetHomemade(Bool_t val) {fHomemade=val;}
  void SetExcludeTPCEdges(Bool_t value) {fExcludeTPCEdges=value;}
  void SetMaxRapidity(Double_t val) {fDecayMaxRapidity=val;}
  void SetMinEta(Double_t val) {fDecayMinEta=val;}
  void SetMaxEta(Double_t val) {fDecayMaxEta=val;}
  void SetMinPt(Double_t val) {fDecayMinPt=val;}
  void SetMaxDCAdaughters(Double_t val) {fDecayMaxDCAdaughters=val;}
  void SetMinCosinePointingAngleXY(Double_t val) {fDecayMinCosinePointingAngleXY=val;}
  void SetMinQt(Double_t val, Bool_t val2=kTRUE) {fDecayMinQt=val; fDecayAPCutPie=val2;}
  void SetStopPIDAtPt(Double_t val) {fDecayStopPIDAtPt=val;}
  void SetMinRadXY(Double_t val) {fDecayMinRadXY=val;}
  void SetMaxDecayLength(Double_t val) {fDecayMaxDecayLength=val;}
  void SetMaxProductIPXY(Double_t val) {fDecayMaxProductIPXY=val;}

  void SetDauMinNClsTPC(Int_t val) {fDaughterMinNClsTPC=val;}
  void SetDauMinNClsITS(Int_t val) {fDaughterMinNClsITS=val;}
  void SetDauMinXRows(Int_t val) {fDaughterMinXRows=val;}
  void SetDauMaxChi2PerNClsTPC(Double_t val) {fDaughterMaxChi2PerNClsTPC=val;}
  void SetDauMinXRowsOverNClsFTPC(Double_t val) {fDaughterMinXRowsOverNClsFTPC=val;}
  void SetDauITSLayer(Int_t layer, Int_t config) {fDaughterITSConfig[layer]=config;}
  void SetDauMinEta(Double_t val) {fDaughterMinEta=val;}
  void SetDauMaxEta(Double_t val) {fDaughterMaxEta=val;}
  void SetDauMinPt(Double_t val) {fDaughterMinPt=val;}
  void SetDauMinImpactParameterXY(Double_t val) {fDaughterMinImpactParameterXY=val;}
  void SetDauMaxNSigmaPID(Double_t val) {fDaughterMaxNSigmaPID=val;}
  void SetDauUnTagProcedure(Bool_t val) {fDaughterUnTag=val;}
  void SetDauSPDRequireAny(Bool_t val) {fDaughterSPDRequireAny=val;}
  void SetDauITSrefit(Bool_t val) {fDaughterITSrefit=val;}

  //newITScuts
  void SetMaxSharedITSCluster(Int_t maxITSclusterShared){fmaxSharedITSCluster = maxITSclusterShared;}
  void SetMaxChi2perITSCluster(Double_t maxITSChi2){fMaxchi2perITSCluster = maxITSChi2;}
  void OpenToyModel();
  void MakeToyEvent(Int_t seed=0, Int_t m_decay = 30, Double_t v_decay = 0.05,
		    Double_t mass_decay_mu = 0.497648, Double_t mass_decay_sg = 0.01,
		    Int_t m_bgr = 30, Double_t v_bgr = 0.08,
		    Int_t mtpc_a = 300, Double_t v_tpca = 0.10, Int_t mtpc_c = 300, Double_t v_tpcc = 0.10,
		    Int_t mvze_a = 300, Double_t v_vzea = 0.10, Int_t mvze_c = 300, Double_t v_vzec = 0.10 );
  void CloseToyModel();
  TList* RebinDecayVn(Int_t nbins, Int_t *bins);

 private:
  AliAnalysisTaskFlowStrangee(const AliAnalysisTaskFlowStrangee& analysisTask);
  AliAnalysisTaskFlowStrangee& operator=(const AliAnalysisTaskFlowStrangee& analysisTask);
  void AddQAEvents();
  void AddQACandidates();

  void MyNotifyRun();
  Bool_t CalibrateEvent();
  void Publish();
  
  void AddEventSpy(TString name);
  void FillEventSpy(TString name);

  Bool_t MinimumRequirementsAA(AliAODEvent *tAOD);
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
  void FillMakeQSpy();
  void ComputeChi2VZERO();
  void MakeQVZE(AliVEvent *event);
  void MakeQTPC(AliVEvent *event);
  void MakeQTPC(AliESDEvent *event);
  void MakeQTPC(AliAODEvent *event);
  void AddTPCRFPSpy(TList *val);
  Bool_t PassesRFPTPCCuts(AliESDtrack *myTrack, Double_t aodChi2NDF=0, Float_t aodipxy=0, Float_t aodipz=0);
  void MakeQVectors();
  void ResetContainers();

  void AddCandidates();
  TList* RebinDecayVn(TList *tList,Int_t nbins, Int_t *bins);

  Double_t GetMCDPHI(Double_t phi);

  Double_t CosThetaPointXY(AliESDv0 *me, const AliVVertex *vtx);
  Double_t CosThetaPointXY(AliAODv0 *me, const AliVVertex *vtx);
  Double_t DecayLengthXY(AliESDv0 *me, const AliVVertex *vtx);
  Double_t DecayLengthXY(AliAODv0 *me, const AliVVertex *vtx);
  Double_t DecayLength(AliESDv0 *me, const AliVVertex *vtx);
  Double_t DecayLength(AliAODv0 *me, const AliVVertex *vtx);

  void AddMCParticleSpy(TList *val);
  void FillMCParticleSpy(TString listName, AliAODMCParticle *par);
  void FillMCParticleSpy(TString listName, TParticle *par);

  void AddCandidatesSpy(TList *val, Bool_t fillRes=kFALSE);
  void FillCandidateSpy(TString listName, Bool_t fillRes=kFALSE);

  void AddTrackSpy(TList *val, Bool_t fillRes=kFALSE);
  void FillTrackSpy(TString listName, Bool_t fillRes=kFALSE);

  void AddDecayVn(TList *val);
  void FillDecayVn(TString listName,Double_t ms,Double_t pt,Double_t phi,Double_t eta,Int_t fid1,Int_t fid2);
  void QCStoreDecayVn(TString name);
  void ComputeDecayVn(TString listName);

  void AddTrackVn(TList *val);
  void FillTrackVn(TString listName,Double_t pt,Double_t phi,Double_t eta,Int_t fid);
  void QCStoreTrackVn(TString name);
  void ComputeTrackVn(TString listName);
  Bool_t InQTPC(Int_t id);

  void MakeFilterBits();
  Bool_t PassesFilterBit(AliESDtrack *me);

  void LoadTrack(AliESDtrack *myTrack, Double_t aodChi2NDF=0);
  Bool_t AcceptDaughter(Bool_t strongITS=kTRUE,  Bool_t newITScuts=kFALSE);
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
  void FillVZEQA();
  void FillVZEQA(AliAODEvent *tAOD);

  Int_t RefMult(AliAODEvent *tAOD, Int_t fb);
  Int_t RefMultTPC();
  Int_t RefMultGlobal();

  AliPIDResponse *fPIDResponse; //! PID response object
  AliESDtrackCuts *fFB1;        // filterbit cut equivalent
  AliESDtrackCuts *fFB1024;     // filterbit cut equivalent
  AliFlowEvent   *fTPCevent;    // flow event (needed here due to ev selection)
  AliFlowEvent   *fVZEevent;    // flow event (needed here due to ev selection)
  //  AliFlowTrackCuts* fCutsRP;

  TObjArray      *fCandidates;  // array of selected candidates
  TList          *fList;        // stores the final list of output histograms

  Int_t fRunNumber; // current run number

  Int_t fDebug;   // debug level
  Int_t fQAlevel; // QA plots

  Bool_t fReadESD;       // move back to ESD
  Bool_t fReadMC;        // read MC files
  Bool_t fAddPiToMCReactionPlane; // add pi randomly (MCTUNED)
  Int_t fPostMatched;    // post only (un)matched particles
  Bool_t fAvoidExec;     // avoids Exec
  Bool_t fSkipSelection; // skip decay finder
  Bool_t fSkipVn;        // skip flow computation
  Bool_t fUseFP;         // flow package?
  Bool_t fRunOnpA;       // make task compatible with pA event selection
  Bool_t fRunOnpp;       // make task compatible with pp event selection
  Bool_t fExtraEventRejection; // to reject pile up
  Bool_t fSkipCentralitySelection; // to skip centrality
  TString  fCentMethod; // CC
  Int_t    fCentPerMin; // CC
  Int_t    fCentPerMax; // CC
  Double_t fThisCent;   // CC
  Double_t fV0M; // V0M CC
  Double_t fTRK; // TRK CC
  Double_t fPriVtxZ; // vtxZ
  Double_t fSPDVtxZ; // vtxZ
  Int_t fSPDtracklets; // spd tracklets
  Float_t fVZETotM; // vzero total multiplicity
  Int_t fRefMultTPC; // tpc only multiplicity
  Int_t fRefMultHyb; // hybrid multiplicity

  Double_t fVertexZcut; // cut on main vertex Z

  Bool_t fExcludeTPCEdges; // exclude TPC edges from single track selection

  Int_t  fSpecie;   // K0=>0 L0=>1
  Bool_t fOnline;   // change into online v0 finder
  Bool_t fHomemade; // homemade v0 finder

  Int_t fWhichPsi;  // detector for Psi2

  Bool_t  fVZEsave; // make vze response
  TList  *fVZEload; // adress to calibration file
  TH2D   *fVZEResponse; // vze response vs centrality class
  Double_t fVZEextW[64]; // vze weights
  Bool_t  fVZEmb;   // integrate response (linearity)
  Bool_t  fVZEByDisk; // normalized by disk
  Int_t   fVZECa;   // start of V0C (ring number 0-3)
  Int_t   fVZECb;   // end of V0C (ring number 0-3)
  Int_t   fVZEAa;   // start of V0A (ring number 0-3)
  Int_t   fVZEAb;   // end of V0A (ring number 0-3)
  TList  *fVZEQA;   // address to qalist

  Int_t fHarmonic;  // flow angle order
  Double_t fPsi2;   // best estimation of Psi2
  Double_t fMCEP;   // stores MC EP (when available)
  // VZE QVector
  Double_t fQVZEACos;
  Double_t fQVZEASin;
  Double_t fQVZECCos;
  Double_t fQVZECSin;
  Double_t fQVZEA;
  Double_t fQVZEC;
  Bool_t fVZEWarning;
  // TPC QVector
  Double_t fQTPCACos;
  Double_t fQTPCASin;
  Double_t fQTPCCCos;
  Double_t fQTPCCSin;
  Double_t fQTPC2hCos;
  Double_t fQTPC2hSin;
  Double_t fQTPCA;
  Double_t fQTPCC;
  Int_t fQTPCA_nTracks;
  Int_t fQTPCC_nTracks;
  Int_t fQTPCA_fID[2000];
  Int_t fQTPCC_fID[2000];
  Bool_t fSkipTerminate;

  Int_t    fMassBins; // opens
  Double_t fMinMass;  // mass
  Double_t fMaxMass;  // window
  Int_t fPtBins;        // to shrink
  Double_t fPtBinEdge[100]; // output

  Int_t fRFPFilterBit;    // RFP TPC
  Double_t fRFPminPt;     // RFP TPC
  Double_t fRFPmaxPt;     // RFP TPC
  Double_t fRFPAminEta;   // RFP TPC
  Double_t fRFPAmaxEta;   // RFP TPC
  Double_t fRFPCminEta;   // RFP TPC
  Double_t fRFPCmaxEta;   // RFP TPC
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
  Double_t fDecayDecayLengthLab;        // DECAY
  Double_t fDecayQt;                    // DECAY
  Double_t fDecayAlpha;                 // DECAY
  Double_t fDecayRapidity;              // DECAY
  Double_t fDecayProductIPXY;           // DECAY
  Double_t fDecayIPneg;                 // DECAY
  Double_t fDecayIPpos;                 // DECAY
  Double_t fDecayXneg;                  // DECAY
  Double_t fDecayXpos;                  // DECAY
  Int_t    fDecayIDneg;                 // DECAY
  Int_t    fDecayIDpos;                 // DECAY
  Int_t    fDecayID;                    // DECAY

  Double_t fDecayMatchOrigin;  // MC DECAY
  Double_t fDecayMatchPhi;     // MC DECAY
  Double_t fDecayMatchEta;     // MC DECAY
  Double_t fDecayMatchPt;      // MC DECAY
  Double_t fDecayMatchRadXY;   // MC DECAY

  Double_t fDecayMinEta;                   // DECAY CUTS
  Double_t fDecayMaxEta;                   // DECAY CUTS
  Double_t fDecayMinPt;                    // DECAY CUTS
  Double_t fDecayMaxDCAdaughters;          // DECAY CUTS
  Double_t fDecayMinCosinePointingAngleXY; // DECAY CUTS
  Double_t fDecayMinQt;                    // DECAY CUTS
  Bool_t   fDecayAPCutPie;                 // DECAY CUTS
  Double_t fDecayStopPIDAtPt;              // DECAY CUTS
  Double_t fDecayMinRadXY;                 // DECAY CUTS
  Double_t fDecayMaxDecayLength;           // DECAY CUTS
  Double_t fDecayMaxProductIPXY;           // DECAY CUTS
  Double_t fDecayMaxRapidity;              // DECAY CUTS

  Double_t fDaughterPhi;               // DAUGHTER
  Double_t fDaughterEta;               // DAUGHTER
  Double_t fDaughterPt;                // DAUGHTER
  Int_t    fDaughterNClsTPC;           // DAUGHTER
  Int_t    fDaughterNClsITS;           // DAUGHTER
  Int_t    fDaughterITSConfig[6];      // DAUGHTER
  Int_t    fDaughterCharge;            // DAUGHTER
  Int_t    fDaughterNFClsTPC;          // DAUGHTER
  Int_t    fDaughterNSClsTPC;          // DAUGHTER
  Double_t fDaughterChi2PerNClsTPC;    // DAUGHTER
  Double_t fDaughterXRows;             // DAUGHTER
  Float_t  fDaughterImpactParameterXY; // DAUGHTER
  Float_t  fDaughterImpactParameterZ;  // DAUGHTER
  UInt_t   fDaughterStatus;            // DAUGHTER
  UChar_t  fDaughterITScm;             // DAUGHTER
  Double_t fDaughterNSigmaPID;         // DAUGHTER
  Int_t    fDaughterKinkIndex;         // DAUGHTER
  Double_t fDaughterAtSecPhi;          // DAUGHTER
  Double_t fDaughterAtSecEta;          // DAUGHTER
  Double_t fDaughterAtSecPt;           // DAUGHTER

 //newITScuts
  Int_t    fsharedITSCluster;          // DAUGHTER 
  Double_t fchi2perClusterITS;         // DAUGHTER 
  Int_t    fcounterForSharedCluster;   // DAUGHTER

  Double_t fDaughterMatchPhi;               // MC DAUGHTER
  Double_t fDaughterMatchEta;               // MC DAUGHTER
  Double_t fDaughterMatchPt;                // MC DAUGHTER
  Float_t  fDaughterMatchImpactParameterXY; // MC DAUGHTER
  Float_t  fDaughterMatchImpactParameterZ;  // MC DAUGHTER

  Bool_t   fDaughterUnTag;             // UNTAG PROCEDURE

  Double_t fDaughterMinEta;               // DAUGHTER CUTS
  Double_t fDaughterMaxEta;               // DAUGHTER CUTS
  Double_t fDaughterMinPt;                // DAUGHTER CUTS
  Int_t    fDaughterMinNClsTPC;           // DAUGHTER CUTS
  Int_t    fDaughterMinNClsITS;           // DAUGHTER CUTS
  Int_t    fDaughterMinXRows;             // DAUGHTER CUTS
  Double_t fDaughterMaxChi2PerNClsTPC;    // DAUGHTER CUTS
  Double_t fDaughterMinXRowsOverNClsFTPC; // DAUGHTER CUTS
  Double_t fDaughterMinImpactParameterXY; // DAUGHTER CUTS
  Double_t fDaughterMaxNSigmaPID;         // DAUGHTER CUTS
  Bool_t   fDaughterSPDRequireAny;        // DAUGHTER CUTS
  Bool_t   fDaughterITSrefit;             // DAUGHTER CUTS
 
  //newITScuts
  Double_t fMaxchi2perITSCluster;        // DAUGHTER CUTS
  Int_t fmaxSharedITSCluster;            // DAUGHTER CUTS

  ClassDef(AliAnalysisTaskFlowStrangee, 6);
};
#endif
