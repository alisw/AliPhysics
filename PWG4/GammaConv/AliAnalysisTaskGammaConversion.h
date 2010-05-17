#ifndef ALIANALYSISTASKGAMMACONVERSION_H
#define ALIANALYSISTASKGAMMACONVERSION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class used to do analysis on conversion pairs
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskSE.h"
#include <vector>
#include "AliV0Reader.h"
//#include "AliCFManager.h"  // for CF
//#include "AliCFContainer.h"   // for CF

class AliGammaConversionAODObject;
class TNtuple;
class AliGammaConversionHistograms;
class AliESDv0;
class AliKFParticle;
class AliESDInputHandler;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class TList;
class AliStack;
class AliESDtrackCuts;
class AliCFManager; // for CF
class AliCFContainer; // for CF

class AliAnalysisTaskGammaConversion : public AliAnalysisTaskSE
{
	
 public:
  AliAnalysisTaskGammaConversion();
  AliAnalysisTaskGammaConversion(const char* name);
  virtual ~AliAnalysisTaskGammaConversion() ;// virtual destructor
		
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  //virtual void ConnectInputData(Option_t * option);
		
  void ProcessMCData();
  void ProcessV0sNoCut();
  void ProcessV0s();
  void ProcessGammasForNeutralMesonAnalysis();
  void ProcessGammasForOmegaMesonAnalysis();
  void ProcessConvPHOSGammasForNeutralMesonAnalysis();
		
  // for CF
  void SetCFManager(AliCFManager * const io) {fCFManager = io;};
  AliCFManager *GetCFManager() const {return fCFManager;}
		
		
  // AOD
  TString GetAODBranchName() const {return  fAODBranchName;}
  void SetAODBranchName(TString name)  {fAODBranchName = name ;}	
  void FillAODWithConversionGammas();
  // end AOD
		
  static Bool_t IsGoodImpPar(AliESDtrack *const track);
	
  // for GammaJetAnalysis
  void ProcessGammasForGammaJetAnalysis();
  void CreateListOfChargedParticles();
  Double_t GetMinimumDistanceToCharge(Int_t indexHighestPtGamma);
  void CalculateJetCone(Int_t gammaIndex);
  Int_t GetIndexHighestPtGamma();
  void SetESDtrackCuts();
  // end of Gamma Jet
		
  void SetMinPtForGammaJet(Double_t minPtForGammaJet){fMinPtForGammaJet=minPtForGammaJet;}
  void SetMinIsoConeSize(Double_t minIsoConeSize){fMinIsoConeSize=minIsoConeSize;}
  void SetMinPtIsoCone(Double_t minPtIsoCone){fMinPtIsoCone=minPtIsoCone;}
  void SetMinPtGamChargedCorr(Double_t minPtGamChargedCorr){fMinPtGamChargedCorr=minPtGamChargedCorr;}
  void SetMinPtJetCone(Double_t minPtJetCone){fMinPtJetCone=minPtJetCone;}

  void SetLowPtMapping(Double_t lowPtMapping){fLowPtMapping=lowPtMapping;}
  void SetHighPtMapping(Double_t highPtMapping){fHighPtMapping=highPtMapping;}

		
  void SetHistograms(AliGammaConversionHistograms *const histograms){fHistograms=histograms;}
  void SetTriggerCINT1B(Bool_t flag){fTriggerCINT1B=flag;}
  void SetDoMCTruth(Bool_t flag){fDoMCTruth=flag;}
  void SetDoNeutralMeson(Bool_t flag){fDoNeutralMeson=flag;}
  void SetDoOmegaMeson(Bool_t flag){fDoOmegaMeson=flag;}
  void SetDoNeutralMesonV0MCCheck(Bool_t flag){fDoNeutralMesonV0MCCheck=flag;}
  void SetDoJet(Bool_t flag){fDoJet=flag;}
  void SetDoChic(Bool_t flag){fDoChic=flag;}
		
  void SetElectronMass(Double_t electronMass){fElectronMass = electronMass;}
  void SetGammaMass(Double_t gammaMass){fGammaMass = gammaMass;}
  void SetGammaWidth(Double_t gammaWidth){fGammaWidth = gammaWidth;}
  void SetPi0Mass(Double_t pi0Mass){fPi0Mass = pi0Mass;}
  void SetPi0Width(Double_t pi0Width){fPi0Width = pi0Width;}
  void SetEtaMass(Double_t etaMass){fEtaMass = etaMass;}
  void SetEtaWidth(Double_t etaWidth){fEtaWidth = etaWidth;}
  void SetMinOpeningAngleGhostCut(Double_t ghostCut){fMinOpeningAngleGhostCut = ghostCut;}
  void SetV0Reader(AliV0Reader* const reader){fV0Reader=reader; fV0Reader->SetESDtrackCuts(fEsdTrackCuts);}
  void SetCalculateBackground(Bool_t bg){fCalculateBackground=bg;}
  void CalculateBackground();
  void SetWriteNtuple(Bool_t writeNtuple){fWriteNtuple = writeNtuple;}
  void FillNtuple();
  Double_t GetMCOpeningAngle(TParticle* const daughter0, TParticle* const daughter1) const;
  void CheckV0Efficiency();

		
  //////////////////Chi_c Analysis////////////////////////////
  void GetPID(AliESDtrack *track, Stat_t &pid, Stat_t &weight);	
  double GetSigmaToVertex(AliESDtrack* t);
  void ElectronBackground(TString hBg, TClonesArray e);
  void FillAngle(TString histoName,TClonesArray const tlVeNeg, TClonesArray const tlVePos);
  void FillElectronInvMass(TString histoName, TClonesArray const negativeElectron, TClonesArray const positiveElectron);
  void FillGammaElectronInvMass(TString histoMass,TString histoDiff, TClonesArray const fKFGammas, TClonesArray const tlVeNeg,TClonesArray const tlVePos);
  void CleanWithAngleCuts(TClonesArray const negativeElectrons, TClonesArray const positiveElectrons, TClonesArray const gammas);
  TClonesArray GetTLorentzVector(TClonesArray* esdTrack);	
  void ProcessGammaElectronsForChicAnalysis();
  ///////////////////////////////////////////////////////////////

  void SetDoCF(Bool_t flag){fDoCF = flag;}
		
 private:
  AliAnalysisTaskGammaConversion(const AliAnalysisTaskGammaConversion&); // Not implemented
  AliAnalysisTaskGammaConversion& operator=(const AliAnalysisTaskGammaConversion&); // Not implemented
		
  // for CF
  enum{
    kStepGenerated = 0,
    kStepReconstructable = 1, 
    kStepGetOnFly = 2,
    kStepLikeSign = 3,
    kStepTPCRefit = 4,
    kStepKinks = 5,
    kStepdEdx_electronselection = 6,
    kStepdEdx_pionrejection = 7,
    kStepNContributors = 8,
    kStepTPCPID = 9,
    kStepR = 10,
    kStepLine = 11,
    kStepZ = 12,
    kStepMinClsTPC = 13,
    kStepSinglePt = 14,
    kStepNDF = 15,
    kStepChi2 = 16,
    kStepEta = 17,
    kStepPt = 18,
    kStepTrueGamma = 19
  };
  
  AliV0Reader* fV0Reader; // The V0 reader object 
		
  AliStack * fStack; // pointer to the MC particle stack
  AliMCEventHandler *fMCTruth; // for CF   pointer to MCTruth
  AliMCEvent *fGCMCEvent;  // for CF    pointer to the MC Event
  AliESDEvent* fESDEvent; //pointer to the ESDEvent
  TList * fOutputContainer; // Histogram container
  AliCFManager *fCFManager;  // for CF
  //  AliCFContainer *container;   // for CF

		
		
  AliGammaConversionHistograms *fHistograms; // Pointer to the histogram handling class
  Bool_t fTriggerCINT1B; //Flag to select trigger CINT1B
  Bool_t fDoMCTruth; // Flag to switch on/off MC truth 
  Bool_t fDoNeutralMeson; // flag
  Bool_t fDoOmegaMeson; // flag
  Bool_t fDoJet; // flag
  Bool_t fDoChic; // flag
		
  TClonesArray * fKFReconstructedGammasTClone; //! transient
  TClonesArray * fKFReconstructedPi0sTClone; //! transient
  TClonesArray * fCurrentEventPosElectronTClone; //! transient
  TClonesArray * fCurrentEventNegElectronTClone; //! transient
  TClonesArray * fKFReconstructedGammasCutTClone; //! transient
  TClonesArray * fPreviousEventTLVNegElectronTClone; //! transient
  TClonesArray * fPreviousEventTLVPosElectronTClone; //! transient
		
  //  vector<AliKFParticle> fKFReconstructedGammas; // vector containing all reconstructed gammas
  vector<Int_t> fElectronv1; // vector containing index of electron 1
  vector<Int_t> fElectronv2; // vector containing index of electron 2
		
  vector<Int_t> fGammav1; // vector containing index of gamma 1
  vector<Int_t> fGammav2; // vector containing index of gamma 2
		


  ///////Chi_c Analysis///////////////////////////
  //  vector<AliESDtrack*> fCurrentEventPosElectron;       // comment here
  //  vector<AliESDtrack*> fCurrentEventNegElectron;       // comment here
  //  vector<AliKFParticle> fKFReconstructedGammasCut;     // comment here
  //  vector<TLorentzVector> fPreviousEventTLVNegElectron; // comment here
  //  vector<TLorentzVector> fPreviousEventTLVPosElectron; // comment here
  //////////////////////////////////////////////////	
		
  //mass defines
  Double_t fElectronMass; //electron mass
  Double_t fGammaMass;    //gamma mass
  Double_t fPi0Mass;      //pi0mass
  Double_t fEtaMass;      //eta mass
		
  // width defines
  Double_t fGammaWidth; //gamma width cut
  Double_t fPi0Width;   // pi0 width cut
  Double_t fEtaWidth;   // eta width cut
		
  Double_t fMinOpeningAngleGhostCut; // minimum angle cut
		
  AliESDtrackCuts* fEsdTrackCuts;           // Object containing the parameters of the esd track cuts
		
  Bool_t fCalculateBackground; //flag to set backgrount calculation on/off
  Bool_t fWriteNtuple;         // flag to set if writing to ntuple on/off
  TNtuple *fGammaNtuple;       // Ntuple for gamma values
  TNtuple *fNeutralMesonNtuple;// NTuple for mesons
		
  Int_t fTotalNumberOfAddedNtupleEntries; // number of added ntuple entries
		
  TClonesArray* fChargedParticles;  //! transient
  vector<Int_t> fChargedParticlesId;  //! transient
		
  Double_t fGammaPtHighest;  //! transient
  Double_t fMinPtForGammaJet;  //! transient
  Double_t fMinIsoConeSize; //! transient
  Double_t fMinPtIsoCone; //! transient
  Double_t fMinPtGamChargedCorr; //! transient
  Double_t fMinPtJetCone; //! transient
  Int_t    fLeadingChargedIndex; //! transient
  Double_t fLowPtMapping; //! transient
  Double_t fHighPtMapping; //! transient
  Bool_t fDoCF; //! transient
		
  TClonesArray* fAODBranch ;        //! selected particles branch
  TString fAODBranchName; // New AOD branch name

  Bool_t fDoNeutralMesonV0MCCheck; //flag

  vector<Int_t>fKFReconstructedGammasV0Index; // index of the reconstructed v0s

  ClassDef(AliAnalysisTaskGammaConversion, 8); // Analysis task for gamma conversions
};

#endif //ALIANALYSISTASKGAMMA_H
