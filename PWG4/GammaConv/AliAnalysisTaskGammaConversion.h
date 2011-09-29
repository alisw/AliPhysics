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
#include "AliGammaConversionBGHandler.h"
#include "TRandom3.h"
#include "TF1.h"
#include "AliMultiplicity.h"
//#include "AliCFManager.h"  // for CF
//#include "AliCFContainer.h"   // for CF
#include "AliAODConversionMother.h"


class AliAODPWG4Particle;
class AliAODConversionPhoton;
class AliAODConversionMother;
class AliKFConversionPhoton;
class AliKFConversionMother;
class TNtuple;
class AliGammaConversionHistograms;
class AliESDv0;
class AliV0;
class AliKFParticle;
class AliKFVertex;
class AliESDInputHandler;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class TList;
class AliStack;
class AliESDtrackCuts;
class AliTriggerAnalysis;
class AliCFManager; // for CF
class AliCFContainer; // for CF
class TRandom3;
class TF1;

class AliAnalysisTaskGammaConversion : public AliAnalysisTaskSE
{
	
 public:
  typedef enum { kProcSD, kProcDD, kProcND, kProcUnknown, kNProcs } ProcType_t; 
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
  void CheckMesonProcessTypeEventQuality(Int_t evtQ);
  Int_t GetProcessType(const AliMCEvent * mcEvt) ;
  void ProcessMCData();
  void ProcessV0sNoCut();
  void ProcessV0s();
  void ProcessGammasForNeutralMesonAnalysis();
  void ProcessGammasForOmegaMesonAnalysis();
  //  void ProcessConvPHOSGammasForNeutralMesonAnalysis();
//  void RecalculateV0ForGamma();
  // for CF
  void SetCFManager(AliCFManager * const io) {fCFManager = io;};
  AliCFManager *GetCFManager() const {return fCFManager;}
		
		
  // AOD
  TString GetAODBranchName() const {return  fAODBranchName;}
  void SetAODBranchName(TString name)  {fAODBranchName = name ;}	
  void SetForceAOD(Bool_t forceAOD ) { fKFForceAOD = forceAOD; }
  void FillAODWithConversionGammas();
  void AddGammaToAOD(AliKFConversionPhoton * kfParticle);
 // void AddPionToAOD(AliKFConversionMother * kfParticle);
 // void AddOmegaToAOD(AliKFParticle * kfParticle, Int_t daughter1, Int_t daughter2);
  void TagDaughter(Int_t gammaIndex);

  // end AOD
		
  static Bool_t IsGoodImpPar(const AliESDtrack *const track);
	
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
  void SetRecalculateV0ForGamma(Bool_t flag){fRecalculateV0ForGamma=flag;}		

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
  Double_t GetMCOpeningAngle(const TParticle* const daughter0,const TParticle* const daughter1) const;
  void CheckV0Efficiency();
  void SetDeltaAODFileName(TString fn) { fKFDeltaAODFileName = fn; };
  void SetCreateAOD(Bool_t doAod) { fKFCreateAOD = doAod; };
  TString GetDeltaAODFileName() const { return fKFDeltaAODFileName; };
  //////////////////Chi_c Analysis////////////////////////////
  void GetPID(const AliESDtrack *track, Stat_t &pid, Stat_t &weight);	
  double GetSigmaToVertex(const AliESDtrack* t);
  void ElectronBackground(TString hBg, TClonesArray e);
  void FillAngle(TString histoName,TClonesArray const tlVeNeg, TClonesArray const tlVePos);
  void FillElectronInvMass(TString histoName, TClonesArray const negativeElectron, TClonesArray const positiveElectron);
  void FillGammaElectronInvMass(TString histoMass,TString histoDiff, TClonesArray const fKFGammas, TClonesArray const tlVeNeg,TClonesArray const tlVePos);
  void CleanWithAngleCuts(TClonesArray const negativeElectrons, TClonesArray const positiveElectrons, TClonesArray const gammas);
  TClonesArray GetTLorentzVector(TClonesArray* esdTrack);	
  void ProcessGammaElectronsForChicAnalysis();
  ///////////////////////////////////////////////////////////////

  void SetDoCF(Bool_t flag){fDoCF = flag;}

  void SetUseChargedTracksMultiplicityForBG(Bool_t flag){fUseTrackMultiplicityForBG = flag;}		

  void SetMoveParticleAccordingToVertex(Bool_t flag){fMoveParticleAccordingToVertex = flag;}

  void MoveParticleAccordingToVertex(AliKFParticle * particle,const AliGammaConversionBGHandler::GammaConversionVertex *vertex);

  void SetApplyChi2Cut(Bool_t flag){fApplyChi2Cut = flag;}

  void SetDoRotation(Bool_t flag){fDoRotation = flag;}

  void SetPMDegreesBG(Int_t deg){fNDegreesPMBackground=deg;}

  void SetNumberOfRotationsBG(Int_t nRot){fNRandomEventsForBG=nRot;}

  void RotateKFParticle(AliKFParticle * kfParticle,Double_t angle);

  void SetCheckBGProbability(Bool_t flag){fCheckBGProbability = flag;}

  void SetRemovePileUp(Bool_t removePileUp) { fRemovePileUp = removePileUp; }

  void SetSelectV0AND(Bool_t selectV0AND) { fSelectV0AND = selectV0AND; }
  void SetUseMultiplicity(Int_t useMultiplicity) {fUseMultiplicity=useMultiplicity;}
  void SetUseMultiplicityBin(Int_t useMultiplicityBin) {fUseMultiplicityBin=useMultiplicityBin;}
  void SetUseHBTMultiplicity(Int_t useHBTMultiplicity) {fUseHBTMultiplicity=useHBTMultiplicity;}
  void SetUseHBTMultiplicityBin(Int_t useHBTMultiplicityBin) {fUseHBTMultiplicityBin=useHBTMultiplicityBin;}

  Int_t CalculateMultiplicityBin();
  void SetUseCentrality(Int_t useCentrality) {fUseCentrality=useCentrality;}
  void SetUseCentralityBin(Int_t useCentralityBin) {fUseCentralityBin=useCentralityBin;}




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
    kStepdEdxElectronSelection = 6,
    kStepdEdxPpionRejection = 7,
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
  Bool_t fRecalculateV0ForGamma;//flag
		
  TClonesArray * fKFReconstructedGammasTClone; //! transient
  TClonesArray * fKFReconstructedPi0sTClone; //! transient
  TClonesArray * fKFRecalculatedGammasTClone; //! transient
  TClonesArray * fCurrentEventPosElectronTClone; //! transient
  TClonesArray * fCurrentEventNegElectronTClone; //! transient
  TClonesArray * fKFReconstructedGammasCutTClone; //! transient
  TClonesArray * fPreviousEventTLVNegElectronTClone; //! transient
  TClonesArray * fPreviousEventTLVPosElectronTClone; //! transient
		
  //  vector<AliKFParticle> fKFReconstructedGammas; // vector containing all reconstructed gammas
		
  //  AliESDpid * fESDpid; // esd pid


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
		
  TClonesArray * fAODGamma; //TClonesArray for gammas to put in AOD
  //TClonesArray * fAODPi0; //TTClonesArray for Pi0s to put in AOD
  //TClonesArray * fAODOmega; //TTClonesArray for omegas to put in AOD
  TString fAODBranchName; // New AOD branch name
  Bool_t fKFCreateAOD; //Create the AOD tclones? (regardless if storing or not)
  
  Bool_t fKFForceAOD;  //Set the Analysis Manager FillAOD variable to true every event
  TString fKFDeltaAODFileName; //! File name for delta AOD (if any)
  Bool_t fDoNeutralMesonV0MCCheck; //flag
  Bool_t fUseTrackMultiplicityForBG; //flag
  Bool_t fMoveParticleAccordingToVertex; //flag
  Bool_t fApplyChi2Cut; //flag
  Int_t fNRandomEventsForBG; //number of random events to use in rotation method
  Int_t fNDegreesPMBackground; // number of degree window to rotate particles for rotation method
  Bool_t fDoRotation; //flag
  Bool_t fCheckBGProbability; //flag
//  vector<Int_t>fKFReconstructedGammasV0Index; // index of the reconstructed v0s
  Bool_t fRemovePileUp;                 // Remove Pile Up
  Bool_t fSelectV0AND;                 // Select V0AND
  AliTriggerAnalysis *fTriggerAnalysis; //! Trigger Analysis for Normalisation
  Int_t fMultiplicity;
  Int_t fUseMultiplicity;
  Int_t fUseMultiplicityBin;
  Int_t fUseHBTMultiplicity;
  Int_t fUseHBTMultiplicityBin;
  Int_t fUseCentrality;
  Int_t fUseCentralityBin;
  TRandom3 fRandom;
  ClassDef(AliAnalysisTaskGammaConversion, 20); // Analysis task for gamma conversions
};

#endif //ALIANALYSISTASKGAMMA_H
