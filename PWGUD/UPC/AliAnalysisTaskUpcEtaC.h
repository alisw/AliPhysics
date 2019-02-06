/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKUPCETAC_H
#define ALIANALYSISTASKUPCETAC_H

class TClonesArray;
class TTree;
class TH1;
class TH2;
class TList;
class AliPIDResponse;
class AliAODEvent;
class AliESDEvent;
class AliTOFTriggerMask;

#define ntrg 17
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskUpcEtaC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcEtaC();
  AliAnalysisTaskUpcEtaC(const char *name);
  virtual ~AliAnalysisTaskUpcEtaC();

  virtual void Init();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void RunAODtrig();
  virtual void RunAODhist();
  virtual void RunAODtree();
  virtual void RunAODMC(AliAODEvent *aod);
  virtual void RunAODsystematics(AliAODEvent *aod);
  virtual void RunESDtrig();
  virtual void RunESDhist();
  virtual void RunESDtree();
  virtual void RunESDMC(AliESDEvent *esd);
  virtual void Terminate(Option_t *);
  void SetRunTree(Bool_t runTree){fRunTree = runTree;}
  void SetRunHist(Bool_t runHist){fRunHist = runHist;}
  void SetRunSyst(Bool_t runSyst){fRunSystematics = runSyst;}
  void SetTracking(Int_t tracking){fTracking = tracking;}
  void SetIsMC(Bool_t MC){isMC = MC;}
  void InitSystematics();
  Double_t GetMedian(Double_t *daArray);
  Bool_t CheckMeritCutWinner(int cutChoice, double oldPars[3], double newPars[3]);
  void BoostCut(TLorentzVector d1, TLorentzVector d2, TLorentzVector parent, Double_t *boostInfo);


 private:
  Int_t fType; // 0 - ESD, 1 - AOD
  Int_t fTracking; //0 - Global, 1 - ITSsa
  Bool_t isMC;
  Bool_t fRunTree; 
  Bool_t fRunHist;
  Bool_t fRunSystematics;
  
  int fMeritCutChoice; //set to 4 by default

  AliPIDResponse *fPIDResponse;
  
  //event tree
  TTree *fEtaCK0sChannelTree;
  TTree *fEtaCTree;
  //tree variables
  Int_t fRunNum;
  UInt_t fPerNum, fOrbNum;
  //trigger
  Bool_t fTrigger[ntrg];
  Bool_t fTriggerInputsMC[ntrg];
  UInt_t fL0inputs, fL1inputs;
  AliTOFTriggerMask *fTOFmask;
  Bool_t fIsPhysicsSelected;

  //PID for EtC->Pi+Pi-K+K- Channel  
  Double_t fPIDTPCMuon[7];
  Double_t fPIDTPCElectron[7];
  Double_t fPIDTPCPion[7];
  Double_t fPIDTPCKaon[7];
  Double_t fPIDTPCProton[7];
  
  Double_t fPIDTOFMuon[7];
  Double_t fPIDTOFElectron[7];
  Double_t fPIDTOFPion[7];
  Double_t fPIDTOFKaon[7];
  Double_t fPIDTOFProton[7];

  //PID for EtC->K0s K+- Pi-+ Channel
  Double_t fPIDTPCMuonPos[7];
  Double_t fPIDTPCElectronPos[7];
  Double_t fPIDTPCPionPos[7];
  Double_t fPIDTPCKaonPos[7];
  Double_t fPIDTPCProtonPos[7];
  
  Double_t fPIDTOFMuonPos[7];
  Double_t fPIDTOFElectronPos[7];
  Double_t fPIDTOFPionPos[7];
  Double_t fPIDTOFKaonPos[7];
  Double_t fPIDTOFProtonPos[7];
 
  Double_t fPIDTPCMuonNeg[7];
  Double_t fPIDTPCElectronNeg[7];
  Double_t fPIDTPCPionNeg[7];
  Double_t fPIDTPCKaonNeg[7];
  Double_t fPIDTPCProtonNeg[7];
  
  Double_t fPIDTOFMuonNeg[7];
  Double_t fPIDTOFElectronNeg[7];
  Double_t fPIDTOFPionNeg[7];
  Double_t fPIDTOFKaonNeg[7];
  Double_t fPIDTOFProtonNeg[7];
  
  Int_t fVtxContrib;
  Double_t fVtxPos[3];
  Double_t fMCVtxPos[3];
  Double_t fVtxErr[3];
  Double_t fVtxChi2,fVtxNDF;
  Double_t fKfVtxPos[3];
  Int_t fSpdVtxContrib;
  Double_t fSpdVtxPos[3];
  
  Bool_t fIsVtxContributor[7];
  
  UShort_t fBCrossNum, fNtracklets, fNLooseTracks;
  //vzero, zdc
  Double_t fZNAenergy, fZNCenergy;
  Double_t fZPAenergy, fZPCenergy;
  Double_t fZDCAtime, fZDCCtime;
  Int_t fV0Adecision, fV0Cdecision;
  Int_t fADAdecision, fADCdecision;
  //input data
  TObjString *fDataFilnam;
  Short_t fRecoPass;
  Long64_t fEvtNum;
  //tracks
  TClonesArray *fJPsiAODTracks;
  TClonesArray *fJPsiESDTracks; 
  TClonesArray *fEtaCAODTracks;
  TClonesArray *fEtaCESDTracks;
    //mc
  TClonesArray *fGenPart;
  
  TList *fListTrig;
  TH1D *fHistCcup4TriggersPerRun;
  TH1D *fHistCcup7TriggersPerRun;
  TH1D *fHistCcup2TriggersPerRun;
  TH1D *fHistCint1TriggersPerRun;
  TH1D *fHistCint6TriggersPerRun;
  TH1D *fHistC0tvxAndCint1TriggersPerRun;
  TH1D *fHistZedTriggersPerRun;
  TH1D *fHistCvlnTriggersPerRun;
  TH1D *fHistMBTriggersPerRun;
  TH1D *fHistCentralTriggersPerRun;
  TH1D *fHistSemiCentralTriggersPerRun;
  
  TH1D *fHistCTest58TriggersPerRun;
  TH1D *fHistCTest59TriggersPerRun;
  TH1D *fHistCTest60TriggersPerRun;
  TH1D *fHistCTest61TriggersPerRun;
  
  TH1D *fHistCcup8TriggersPerRun;
  TH1D *fHistCcup9TriggersPerRun;
  TH1D *fHistCcup10TriggersPerRun;
  TH1D *fHistCcup11TriggersPerRun;
  TH1D *fHistCcup12TriggersPerRun;
  TH1D *fHistCtrueTriggersPerRun;

  //My histos and stuff
  TList *fListHist;
  TList *fListHistKstar;
  TList *fListHist2Rho4Pion;
  TList *fListHistK0s3PiPi4K;
  TList *fListHistZDC;

  //New ZDC histos
  TH1D *fHistZDCAenergy;
  TH1D *fHistZDCCenergy;
  TH1D *fHistZDCAtime;
  TH1D *fHistZDCCtime;
  TH1D *fHistZDCImpactParameter;
  TH1D *fHistZDCAImpactParameter;
  TH1D *fHistZDCCImpactParameter;

  TH1D *fHistNeventsEtaC; //Count potential EtaC events at each step
  TH2D *fMPiKvsMPiK; //Dalitz Plot, Mass first PiK combo vs Mass second PiK combo
  //2 Kstar case
  TH1D *f2KstarPtPiPlus;
  TH1D *f2KstarPtPiMinus;
  TH1D *f2KstarPtKPlus;
  TH1D *f2KstarPtKMinus;
  TH2D *f2KstarTPCsignalPion; //dEdx correlation Pions
  TH2D *f2KstarTPCsignalKaon; //dEdx correlation Kaons
  TH2D *f2KstarDedxVsPtPion;
  TH2D *f2KstarDedxVsPtKaon;
  TH2D *f2KstarTPCsignalVsQPtPion; //dEdx Pions vs Q*Pt
  TH2D *f2KstarTPCsignalVsQPtKaon; //dEdx Kaons vs Q*Pt
  TH2D *f2KstarPtVsMinvFirstKstar; 
  TH2D *f2KstarPtVsMinvSecondKstar;
  //    TH1D *f2KstarMinvFirstKstar;
  //    TH1D *f2KstarMinvSecondKstar;
  TH2D *f2KstarPtVsMinvEtaC;
  TH2D *f2KstarEtaVsMinvEtaC;
  TH2D *f2KstarEtaVsMinvEtaC400MeVPtMax;
  TH2D *f2KstarEtaVsMinvEtaC100MeVPtMax;
  TH2D *f2KstarSumPzVsMinvEtaC;
  TH1D *f2KstarScalarSumP;
  TH1D *f2KstarVectorSumPt;
  //    TH1D *f2KstarMinvEtaC;
  //1 Kstar case
  TH1D *f1KstarPtPiPlus;
  TH1D *f1KstarPtPiMinus;
  TH1D *f1KstarPtKPlus;
  TH1D *f1KstarPtKMinus;
  TH2D *f1KstarTPCsignalPion; //dEdx correlation Pions
  TH2D *f1KstarTPCsignalKaon; //dEdx correlation Kaons
  TH2D *f1KstarDedxVsPtPion;
  TH2D *f1KstarDedxVsPtKaon;
  TH2D *f1KstarTPCsignalVsQPtPion; //dEdx Pions vs Q*Pt
  TH2D *f1KstarTPCsignalVsQPtKaon; //dEdx Kaons vs Q*Pt
  TH2D *f1KstarPtVsMinvKstar;
  TH2D *f1KstarPtVsMinvOtherPiKcombo;
  //    TH1D *f1KstarMinvKstar;
  //    TH1D *f1KstarMinvOtherPiKcombo;
  TH2D *f1KstarPtVsMinvEtaC;
  TH2D *f1KstarEtaVsMinvEtaC;
  TH2D *f1KstarEtaVsMinvEtaC400MeVPtMax;
  TH2D *f1KstarEtaVsMinvEtaC100MeVPtMax;
  TH2D *f1KstarSumPzVsMinvEtaC;
  TH1D *f1KstarScalarSumP;
  TH1D *f1KstarVectorSumPt;
  //    TH1D *f1KstarMinvEtaC;
  //0 Kstar case
  TH1D *f0KstarPtPiPlus;
  TH1D *f0KstarPtPiMinus;
  TH1D *f0KstarPtKPlus;
  TH1D *f0KstarPtKMinus;
  TH2D *f0KstarTPCsignalPion; //dEdx correlation Pions
  TH2D *f0KstarTPCsignalKaon; //dEdx correlation Kaons
  TH2D *f0KstarDedxVsPtPion;
  TH2D *f0KstarDedxVsPtKaon;
  TH2D *f0KstarTPCsignalVsQPtPion; //dEdx Pions vs Q*Pt
  TH2D *f0KstarTPCsignalVsQPtKaon; //dEdx Kaons vs Q*Pt
  TH2D *f0KstarPtVsMinvFirstPiKcombo;
  TH2D *f0KstarPtVsMinvSecondPiKcombo;
  //    TH1D *f0KstarMinvFirstPiKcombo;
  //    TH1D *f0KstarMinvSecondPiKcombo;
  TH2D *f0KstarPtVsMinvEtaC;
  TH2D *f0KstarEtaVsMinvEtaC;
  TH2D *f0KstarEtaVsMinvEtaC400MeVPtMax;
  TH2D *f0KstarEtaVsMinvEtaC100MeVPtMax;
  TH2D *f0KstarSumPzVsMinvEtaC;
  TH1D *f0KstarScalarSumP;
  TH1D *f0KstarVectorSumPt;
  //    TH1D *f0KstarMinvEtaC;

  //K0s Channel
  TH1D *fHistK0sCandidatesPerEvent; //Track number of K0s candidates per event
  TH1D *fK0sPosDaughterPt;
  TH1D *fK0sNegDaughterPt;
  TH2D *fK0sPosVsNegDaughterPt;
  TH1D *fK0sPionPt;
  TH1D *fK0sKaonPt;
  TH2D *fK0sPtVsMinvK0s;
  //  TH1D *fK0sMinv;
  TH2D *fKPiPtVsMinvK0sChannel;
  //  TH1D *fKPiMinvK0sChannel;
  TH2D *fM2K0sVsM2KPiK0sChannel; //Dalitz Plot, Mass K0s vs Mass of PiK combo
  TH2D *fM2K0sPiVsM2KPiK0sChannel; //Dalitz Plot, Mass K0sPi vs Mass of PiK combo
  TH2D *fM2K0sKVsM2KPiK0sChannel; //Dalitz Plot, Mass K0sK vs Mass of PiK combo
  TH2D *fK0sPtVsMinvEtaC;
  //  TH1D *fEtaCMinvK0sChannel;
  TH1D *fK0sDecayLength;
  TH2D *fK0sEtaVsMinvEtaC;
  TH2D *fK0sEtaVsMinvEtaC400MeVPtMax;
  TH2D *fK0sEtaVsMinvEtaC100MeVPtMax;
  TH2D *fK0sSumPzVsMinvEtaC;
  TH1D *fK0sScalarSumP;
  TH1D *fK0sVectorSumPt;

  TH2D *fHistEtaCMassVsPt;
  TH1D *fHistEtaCMassCoherent;

  TH1D *fHistNeventsEtaCK0sChannel; //Count potential EtaC events at each step

  //Diagnostic histos to understand what kind of events are being analyzed.
  TH1D *fHistNpion;
  TH1D *fHistNK0sPion;
  TH1D *fHistNkaon;
  TH1D *fHistPiMinusK;

  //New diagnostic histos to investigate alternative PID approaches.
  TH2D *fNSigmaPionTPCvsNSigmaPionTOFLowPt;
  TH2D *fNSigmaPionTPCvsNSigmaPionTOFMidPt;
  TH2D *fNSigmaPionTPCvsNSigmaPionTOFHighPt;
  TH2D *fNSigmaKaonTPCvsNSigmaKaonTOFLowPt;
  TH2D *fNSigmaKaonTPCvsNSigmaKaonTOFMidPt;
  TH2D *fNSigmaKaonTPCvsNSigmaKaonTOFHighPt;
  TH2D *fTPCdEdxVsTOFbetaAll;
  TH2D *fTPCdEdxVsTOFbetaPionsWithPID;
  TH2D *fTPCdEdxVsTOFbetaKaonsWithPID;
  TH2D *fTOFTimeVsTPCdEdxAll;
  TH2D *fTOFTimeVsTPCdEdxPionsWithPID;
  TH2D *fTOFTimeVsTPCdEdxKaonsWithPID;
  TH2D *fTOFbetaVsPtAll;
  TH2D *fTOFbetaVsPtPionsWithPID;
  TH2D *fTOFbetaVsPtKaonsWithPID;
  TH1D *fNTracksWithTOFPIDPerEvent;
  TH1D *fNTracksMissingDueToTOFPerEvent;

  TH1D *fV0DaughterDca;
  TH1D *fK0sDcaToPrimVertex;
  TH1D *fK0sDaughterDcaToPrimVertex;
  TH1D *fK0sMassDistribution;
  TH1D *fV0DecayLength;
  TH1D *fV0Eta;
  TH1D *fCosPointingAngle;

  //RhoRho Channel histos.
  TH1D *fHistNeventsEtaCRhoChannel;
  TH2D *f2RhoPtVsMinvRho;
  TH2D *f4PionPtVsMinvRho;
  TH2D *f2RhoPtVsMinvEtaC;
  TH2D *f4PionPtVsMinvEtaC;

  TH2D *f2RhoPtVsMinvOtherRho;
  TH2D *f2RhoPtVsMinvNonRhoPairs;
  TH2D *f4PiVs2PiMinv;
  TH2D *f4PiVs2PiMinvSquared;
  TH2D *fM2PiPiVsM2PiPi;
  TH2D *f2RhoEtaVsMinvEtaC;
  TH2D *f4PionEtaVsMinvEtaC;
  TH2D *f2RhoEtaVsMinvEtaC400MeVPtMax;
  TH2D *f4PionEtaVsMinvEtaC400MeVPtMax;
  TH2D *f2RhoEtaVsMinvEtaC100MeVPtMax;
  TH2D *f4PionEtaVsMinvEtaC100MeVPtMax;
  TH2D *f2RhoSumPzVsMinvEtaC;
  TH2D *f4PionSumPzVsMinvEtaC;
  TH1D *f2RhoScalarSumP;
  TH1D *f4PionScalarSumP;
  TH1D *f2RhoVectorSumPt;
  TH1D *f4PionVectorSumPt;

  //3PiPi Channel histos
  TH1D *fHistNeventsEtaC3PiPiChannel;
  TH2D *f3PiPiPtVsMinvEtaC;
  TH2D *f3PiPiEtaVsMinvEtaC;
  TH2D *f3PiPiEtaVsMinvEtaC400MeVPtMax;
  TH2D *f3PiPiEtaVsMinvEtaC100MeVPtMax;
  TH2D *f3PiPiSumPzVsMinvEtaC;
  TH1D *f3PiPiScalarSumP;
  TH1D *f3PiPiVectorSumPt;
    
  //Helicity cut histos
  TH1D *fKstarParentPx;
  TH1D *fKstarParentPy;
  TH1D *fKstarParentPz;
  TH1D *fKstarDaughterParentAngle;
  TH1D *fKstarDaughterParentCosAngle;
  TH1D *fKstarDaughterDaughterAngle;
  TH1D *fKstarDaughterDaughterCosAngle;
  TH1D *fKstarDaughterPtotal;
  TH1D *fKstarDaughterPtotalNorm;

  //Helicity cut histos - Check histos
  TH1D *fKstarParentPxCheck;
  TH1D *fKstarParentPyCheck;
  TH1D *fKstarParentPzCheck;
  TH1D *fKstarDaughterParentAngleCheck;
  TH1D *fKstarDaughterParentCosAngleCheck;
  TH1D *fKstarDaughterDaughterAngleCheck;
  TH1D *fKstarDaughterDaughterCosAngleCheck;
  TH1D *fKstarDaughterPtotalCheck;
  TH1D *fKstarDaughterPtotalNormCheck;

  //2Rho0 channel Helicity cut histos
  TH1D *f2RhoParentPx;
  TH1D *f2RhoParentPy;
  TH1D *f2RhoParentPz;
  TH1D *f2RhoDaughterParentAngle;
  TH1D *f2RhoDaughterParentCosAngle;
  TH1D *f2RhoDaughterDaughterAngle;
  TH1D *f2RhoDaughterDaughterCosAngle;
  TH1D *f2RhoDaughterPtotal;

  //2Rho0 channel Helicity cut histos - Check histos
  TH1D *f2RhoParentPxCheck;
  TH1D *f2RhoParentPyCheck;
  TH1D *f2RhoParentPzCheck;
  TH1D *f2RhoDaughterParentAngleCheck;
  TH1D *f2RhoDaughterParentCosAngleCheck;
  TH1D *f2RhoDaughterDaughterAngleCheck;
  TH1D *f2RhoDaughterDaughterCosAngleCheck;
  TH1D *f2RhoDaughterPtotalCheck;

  //4 Kaon channel
  TH1D *fHistNeventsEtaC4KaonChannel;
  TH2D *f4KaonPtVsMinvEtaC;
  TH2D *f4KaonPtVsMinvKK;
  TH2D *f4KVs2KMinv;
  TH2D *f4KVs2KMinvSquared;
  TH2D *fM2KKVsM2KK;
  TH2D *f4KaonEtaVsMinvEtaC;
  TH2D *f4KaonEtaVsMinvEtaC400MeVPtMax;
  TH2D *f4KaonEtaVsMinvEtaC100MeVPtMax;
  TH2D *f4KaonSumPzVsMinvEtaC;
  TH1D *f4KaonScalarSumP;
  TH1D *f4KaonVectorSumPt;

  TH1D *fHistZDCCuts;
  
  TList *fListSystematics;
  TList *fListJPsiLoose;
  TList *fListJPsiTight;
  TList *fListEtaCLoose;
  TList *fListEtaCTight;
  
  AliAnalysisTaskUpcEtaC(const AliAnalysisTaskUpcEtaC&); //not implemented
  AliAnalysisTaskUpcEtaC& operator =(const AliAnalysisTaskUpcEtaC&); //not implemented
  
  ClassDef(AliAnalysisTaskUpcEtaC, 5); 
};

#endif





