/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */
//Full code, 18r.pbpb_Legotrain_01
#ifndef AliAnalysisTaskUpcFourPionsGNP_H 
#define AliAnalysisTaskUpcFourPionsGNP_H

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

class AliAnalysisTaskUpcFourPionsGNP : public AliAnalysisTaskSE {
public:
	AliAnalysisTaskUpcFourPionsGNP();
	AliAnalysisTaskUpcFourPionsGNP(const char *name);
	virtual ~AliAnalysisTaskUpcFourPionsGNP();

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
	void SetRunTree(Bool_t runTree) { fRunTree = runTree; }
	void SetRunHist(Bool_t runHist) { fRunHist = runHist; }
	void SetRunSyst(Bool_t runSyst) { fRunSystematics = runSyst; }
	void SetTracking(Int_t tracking) { fTracking = tracking; }
	void SetIsMC(Bool_t MC) { isMC = MC; }
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
	Double_t fVtxChi2, fVtxNDF;
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
	TH1D *fHistCcup29TriggersPerRun;
	TH1D *fHistCcup30TriggersPerRun;
	TH1D *fHistCcup31TriggersPerRun;

	TH1D *fHistZedTriggersPerRun;
	TH1D *fHistMBTriggersPerRun;
	TH1D *fHistCentralTriggersPerRun;
	TH1D *fHistSemiCentralTriggersPerRun;
	TH1D *fHistCtrueTriggersPerRun;

	//My histos and stuff
	TList *fListHist;
	TList *fListHistKstar;
	TList *fListHist2Rho4Pion;
	TList *fListHistK0s3PiPi4K;
	TList *fListHistZDC;

	TList *fListHistPID;
	TList *fListHistUnusedPID;
	TList *fListHistPIDcuts;
	TList *fList4TrackPID;
	TList *fList6TrackPID;

	TList *fList2KStar;
	TList *fList1KStar;
	TList *fList0KStar;
	TList *fListHelicityCuts;
	TList *fList2KStarDiagnostic;
	TList *fList1KStarDiagnostic;
	TList *fList0KStarDiagnostic;
	TList *fList2KStarEtaC;
	TList *fList1KStarEtaC;
	TList *fList0KStarEtaC;

	TList *fList2RhoCandidates;
	TList *fList2Rho;
	TList *fList1Rho;
	TList *fList0Rho;
	TList *fList2RhoDiagnostic;
	TList *fList1RhoDiagnostic;
	TList *fList0RhoDiagnostic;
	TList *fList2RhoEtaC;
	TList *fList1RhoEtaC;
	TList *fList0RhoEtaC;

	TList *fList3PiPi;
	TList *fList4K;
	TList *fListK0Short;
	TList *fList3PiPiDiagnostic;
	TList *fList4KDiagnostic;
	TList *fListK0ShortDiagnostic;
	TList *fList3PiPiEtaC;
	TList *fList4KEtaC;
	TList *fListK0ShortEtaC;
	TList *fListK0ShortPID;
	TList *fList2K4Pi;
	TList *fList2K4PiEtaC;
	TList *fList2K4PiDiagnostic;

	//New ZDC histos
	TH1D *fHistZDCAenergy;
	TH1D *fHistZDCCenergy;
	TH1D *fHistZDCAtime;
	TH1D *fHistZDCCtime;
	TH1D *fHistZDCImpactParameter;
	TH1D *fHistZDCAImpactParameter;
	TH1D *fHistZDCCImpactParameter;

	//Diagnostic histos to understand what kind of events are being analyzed.
	TH1D *fEtaCCandidatesPerChannel;
	TH1D *fEtaCLowPtCandidatesPerChannel;
	TH2D *fAllPtVsMinvEtaC;
	TH1D *fAllMinvEtaCLowPt;
	TH2D *fChannelVsMinvEtaC;
	TH1D *fHistNTracks;
	TH1D *fHistFourTracksNpion;
	TH1D *fHistSixTracksNpion;
	TH1D *fHistNK0sPion;
	TH1D *fHistFourTracksNkaon;
	TH1D *fHistSixTracksNkaon;
	TH1D *fHistPIDCheck;
	TH1D *fFourTrackMissing;
	TH1D *fSixTrackMissing;
	TH1D *fHistFourTracksNboth;
	TH1D *fHistFourTracksNneither;
	TH1D *fHistSixTracksNboth;
	TH1D *fHistSixTracksNneither;
	TH1D *fListNKaonChange;
	TH1D *fListNPionChange;
	TH1D *fHistPostFourTracksNkaon;
	TH1D *fHistPostFourTracksNpion;
	TH2D *f4KnKaonVsnewKaon;
	TH1D *f4PiEventCandidates;
	TH2D *f4PinPionVsnewPion;

	//New diagnostic histos to investigate alternative PID approaches.
	TH2D *fPionTPCvsPionTOFLowP;
	TH2D *fPionTPCvsPionTOFMidP;
	TH2D *fPionTPCvsPionTOFHighP;
	TH2D *fKaonTPCvsKaonTOFLowP;
	TH2D *fKaonTPCvsKaonTOFMidP;
	TH2D *fKaonTPCvsKaonTOFHighP;
	TH2D *fPionTPCvsKaonTPCLowP;
	TH2D *fPionTPCvsKaonTPCMidP;
	TH2D *fPionTPCvsKaonTPCHighP;
	TH2D *fPionTOFvsKaonTOFMidP;
	TH2D *fPionTOFvsKaonTOFPionMidP;
	TH2D *fPionTOFvsKaonTOFKaonMidP;
	TH2D *fTPCdEdxVsTOFbetaAll;
	TH2D *fTPCdEdxVsTOFbetaPionsWithPID;
	TH2D *fTPCdEdxVsTOFbetaKaonsWithPID;
	TH1D *fTOFIntegratedLength;
	TH2D *fTOFbetaVsPAll;
	TH2D *fTOFbetaVsPPion;
	TH2D *fTOFbetaVsPKaon;
	TH2D *fTOFbetaVsPBoth;
	TH2D *fDedxVsPAll;
	TH2D *fDedxVsPPion;
	TH2D *fDedxVsPKaon;
	TH2D *fDedxVsPBoth;

	TH1D *fHistNeventsEtaC; //Count potential EtaC events at each step

							//2 Kstar case
	TH1D *fHistNeventsEtaCKstarChannel;
	TH1D *fKstarEventCandidates;
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
	TH2D *f2KstarPtVsMinvEtaC;
	TH2D *f2KstarEtaVsMinvEtaC;
	TH2D *f2KstarEtaVsMinvEtaC400MeVPtMax;
	TH2D *f2KstarEtaVsMinvEtaC100MeVPtMax;
	TH2D *f2KstarSumPzVsMinvEtaC;
	TH1D *f2KstarScalarSumP;
	TH1D *f2KstarVectorSumPt;

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
	TH2D *f1KstarPtVsMinvEtaC;
	TH2D *f1KstarEtaVsMinvEtaC;
	TH2D *f1KstarEtaVsMinvEtaC400MeVPtMax;
	TH2D *f1KstarEtaVsMinvEtaC100MeVPtMax;
	TH2D *f1KstarSumPzVsMinvEtaC;
	TH1D *f1KstarScalarSumP;
	TH1D *f1KstarVectorSumPt;

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
	TH2D *f0KstarPtVsMinvEtaC;
	TH2D *f0KstarEtaVsMinvEtaC;
	TH2D *f0KstarEtaVsMinvEtaC400MeVPtMax;
	TH2D *f0KstarEtaVsMinvEtaC100MeVPtMax;
	TH2D *f0KstarSumPzVsMinvEtaC;
	TH1D *f0KstarScalarSumP;
	TH1D *f0KstarVectorSumPt;
	TH2D *fMPiKvsMPiK; //Dalitz Plot, Mass first PiK combo vs Mass second PiK combo
	TH2D *fKstarMPiKvsMPiK;

	TH1D *fHistNeventsEtaC2K4PiChannel;
	TH1D *f2K4PiEventCandidates;
	TH2D *f2K4PiPtVsMinvEtaC;
	TH2D *f2K4PiEtaVsMinvEtaC;
	TH2D *f2K4PiEtaVsMinvEtaC400MeVPtMax;
	TH2D *f2K4PiEtaVsMinvEtaC100MeVPtMax;
	TH2D *f2K4PiSumPzVsMinvEtaC;
	TH1D *f2K4PiScalarSumP;
	TH1D *f2K4PiVectorSumPt;

					   //K0s Channel
	TH1D *fHistNeventsEtaCK0sChannel; //Count potential EtaC events at each step
	TH1D *fK0sEventCandidates;
	TH1D *fHistK0sCandidatesPerEvent; //Track number of K0s candidates per event
	TH1D *fK0sPosDaughterPt;
	TH1D *fK0sNegDaughterPt;
	TH2D *fK0sPosVsNegDaughterPt;
	TH1D *fK0sPionPt;
	TH1D *fK0sKaonPt;
	TH2D *fK0sPtVsMinvK0s;
	TH2D *fK0sPtVsMinvKPi;
	TH2D *fK0sM2K0sVsM2KPi; //Dalitz Plot, Mass K0s vs Mass of PiK combo
	TH2D *fK0sM2K0sPiVsM2KPi; //Dalitz Plot, Mass K0sPi vs Mass of PiK combo
	TH2D *fK0sM2K0sKVsM2KPi; //Dalitz Plot, Mass K0sK vs Mass of PiK combo
	TH2D *fK0sPtVsMinvEtaC;
	TH1D *fK0sDecayLength;
	TH2D *fK0sEtaVsMinvEtaC;
	TH2D *fK0sEtaVsMinvEtaC400MeVPtMax;
	TH2D *fK0sEtaVsMinvEtaC100MeVPtMax;
	TH2D *fK0sSumPzVsMinvEtaC;
	TH1D *fK0sScalarSumP;
	TH1D *fK0sVectorSumPt;

	TH1D *fK0DaughterDca;
	TH1D *fK0sDcaToPrimVertex;
	TH1D *fK0sDaughterDcaToPrimVertex;
	TH1D *fK0sMassDistribution;
	TH1D *fV0DecayLength;
	TH1D *fV0Eta;
	TH1D *fV0CosPointingAngle;
	TH1D *fV0sMassK0s;

	TH1D *fHistNProngFound;
	TH1D *fK0GoodTracks;
	TH2D *fK0sTOFbetaVsPAll;
	TH2D *fK0sTOFbetaVsPPion;
	TH2D *fK0sDedxVsPAll;
	TH2D *fK0sDedxVsPPion;

	//RhoRho Channel histos.

	TH1D *fHistNeventsEtaCRhoChannel;
  TH1D *fHistNeventsFourPicharge;
	TH2D *f2Rho2PairPtVsMinvFirstPair;
	TH2D *f2Rho1PairPtVsMinvRhoPair;
	TH2D *f2Rho1PairPtVsMinvNonRhoPair;
	TH2D *f4PionPtVsMinvRho;
	TH2D *f2Rho2PairPtVsMinvEtaC;
	TH2D *f2Rho1PairPtVsMinvEtaC;
	TH2D *f4PionPtVsMinvEtaC;

	TH2D *f2Rho2PairPtVsMinvSecondPair;
	TH2D *f2PairRhoM2VsPiPiM2;
	TH2D *f1PairRhoM2VsPiPiM2;
	TH2D *f2RhoM2VsPiPiM2;
	TH2D *f2Rho2PairEtaVsMinvEtaC;
	TH2D *f2Rho1PairEtaVsMinvEtaC;
	TH2D *f4PionEtaVsMinvEtaC;
	TH2D *f2Rho2PairEtaVsMinvEtaC400MeVPtMax;
	TH2D *f2Rho1PairEtaVsMinvEtaC400MeVPtMax;
	TH2D *f4PionEtaVsMinvEtaC400MeVPtMax;
	TH2D *f2Rho2PairEtaVsMinvEtaC100MeVPtMax;
	TH2D *f2Rho1PairEtaVsMinvEtaC100MeVPtMax;
	TH2D *f4PionEtaVsMinvEtaC100MeVPtMax;
	TH2D *f2Rho2PairSumPzVsMinvEtaC;
	TH2D *f2Rho1PairSumPzVsMinvEtaC;
	TH2D *f4PionSumPzVsMinvEtaC;
	TH1D *f2Rho2PairScalarSumP;
	TH1D *f2Rho1PairScalarSumP;
	TH1D *f4PionScalarSumP;
	TH1D *f2Rho2PairVectorSumPt;
	TH1D *f2Rho1PairVectorSumPt;
	TH1D *f4PionVectorSumPt;

	//3PiPi Channel histos
	TH1D *fHistNeventsEtaC3PiPiChannel;
	TH1D *f6PiEventCandidates;
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
	TH1D *f4KEventCandidates;

	TH1D *fHistZDCCuts;

	TList *fListSystematics;
	TList *fListJPsiLoose;
	TList *fListJPsiTight;
	TList *fListEtaCLoose;
	TList *fListEtaCTight;

	AliAnalysisTaskUpcFourPionsGNP(const AliAnalysisTaskUpcFourPionsGNP&); //not implemented
	AliAnalysisTaskUpcFourPionsGNP& operator =(const AliAnalysisTaskUpcFourPionsGNP&); //not implemented

	ClassDef(AliAnalysisTaskUpcFourPionsGNP, 5);
};

#endif
