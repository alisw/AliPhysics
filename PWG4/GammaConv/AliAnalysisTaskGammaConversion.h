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

class TNtuple;
class AliGammaConversionHistograms;
class AliESDv0;
class AliKFParticle;
class AliESDInputHandler;
class AliESDEvent;
class AliAODEvent;
class TList;
class AliStack;

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
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void ConnectInputData(Option_t *);
	
  void ProcessMCData();
  void ProcessV0sNoCut();
  void ProcessV0s();
  void ProcessGammasForNeutralMesonAnalysis();
  void SetHistograms(AliGammaConversionHistograms *const histograms){fHistograms=histograms;}
  void SetDoMCTruth(Bool_t flag){fDoMCTruth=flag;}
  void SetElectronMass(Double_t electronMass){fElectronMass = electronMass;}
  void SetGammaMass(Double_t gammaMass){fGammaMass = gammaMass;}
  void SetGammaWidth(Double_t gammaWidth){fGammaWidth = gammaWidth;}
  void SetPi0Mass(Double_t pi0Mass){fPi0Mass = pi0Mass;}
  void SetPi0Width(Double_t pi0Width){fPi0Width = pi0Width;}
  void SetEtaMass(Double_t etaMass){fEtaMass = etaMass;}
  void SetEtaWidth(Double_t etaWidth){fEtaWidth = etaWidth;}
  void SetMinOpeningAngleGhostCut(Double_t ghostCut){fMinOpeningAngleGhostCut = ghostCut;}
  void SetV0Reader(AliV0Reader* const reader){fV0Reader=reader;}
  void SetCalculateBackground(Bool_t bg){fCalculateBackground=bg;}
  void CalculateBackground();
  void SetWriteNtuple(Bool_t writeNtuple){fWriteNtuple = writeNtuple;}
  void FillNtuple();
  Double_t GetMCOpeningAngle(TParticle* const daughter0, TParticle* const daughter1) const;
  void CheckV0Efficiency();


  //////////////////Chi_c Analysis////////////////////////////
  void GetPID(AliESDtrack *track, Stat_t &pid, Stat_t &weight);	
  double GetSigmaToVertex(AliESDtrack* t);
  void ElectronBackground(TString hBg, vector <TLorentzVector> e);
  void FillAngle(TString histoName,vector <TLorentzVector> tlVeNeg, vector <TLorentzVector> tlVePos);
  void FillElectronInvMass(TString histoName, vector <TLorentzVector> negativeElectron, 
	vector <TLorentzVector> positiveElectron);
  void FillGammaElectronInvMass(TString histoMass,TString histoDiff,vector <AliKFParticle> fKFGammas,
        vector <TLorentzVector> tlVeNeg,vector<TLorentzVector> tlVePos);
  void CleanWithAngleCuts(vector <AliESDtrack*> negativeElectrons,
	vector <AliESDtrack*> positiveElectrons, vector <AliKFParticle> gammas);
  vector <TLorentzVector> GetTLorentzVector(vector <AliESDtrack*> esdTrack);	
  void ProcessGammaElectronsForChicAnalysis();
  ///////////////////////////////////////////////////////////////



 private:
  AliAnalysisTaskGammaConversion(const AliAnalysisTaskGammaConversion&); // Not implemented
  AliAnalysisTaskGammaConversion& operator=(const AliAnalysisTaskGammaConversion&); // Not implemented

  AliV0Reader* fV0Reader; // The V0 reader object 

  AliStack * fStack; // pointer to the MC particle stack
  AliESDEvent* fESDEvent; //pointer to the ESDEvent
  TList * fOutputContainer ; // Histogram container

  AliGammaConversionHistograms *fHistograms; // Pointer to the histogram handling class

  Bool_t fDoMCTruth; // Flag to switch on/off MC truth 
    
  vector<TParticle*> fMCAllGammas; // vector containing all MC gammas
  vector<TParticle*> fMCPi0s; //vector containing all MC Pi0s
  vector<TParticle*> fMCEtas; //vector containing all MC Etas
  vector<TParticle*> fMCGammaChic; //vector containing all MC Chi_c's

  vector<AliKFParticle> fKFReconstructedGammas; // vector containing all reconstructed gammas
  vector<Bool_t> fIsTrueReconstructedGammas;    // vector containing information if this was a true gamma or not (follows the index of fKFReconstructedGammas)
  vector<Int_t> fElectronv1; // vector containing index of electron 1
  vector<Int_t> fElectronv2; // vector containing index of electron 2

  ///////Chi_c Analysis///////////////////////////
  vector<AliESDtrack*> fCurrentEventPosElectron;       // comment here
  vector<AliESDtrack*> fPreviousEventPosElectron;      //comment here
  vector<AliESDtrack*> fCurrentEventNegElectron;       //comment here
  vector<AliESDtrack*> fPreviousEventNegElectron;      //comment here
  vector<AliKFParticle> fKFReconstructedGammasCut;     //comment here
  vector<TLorentzVector> fPreviousEventTLVNegElectron; //comment here
  vector<TLorentzVector> fPreviousEventTLVPosElectron; //comment here
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

  Bool_t fCalculateBackground; //flag to set backgrount calculation on/off
  Bool_t fWriteNtuple;         // flag to set if writing to ntuple on/off
  TNtuple *fGammaNtuple;       // Ntuple for gamma values
  TNtuple *fNeutralMesonNtuple;// NTuple for mesons

  Int_t fTotalNumberOfAddedNtupleEntries; // number of added ntuple entries

  ClassDef(AliAnalysisTaskGammaConversion, 3); // Analysis task for gamma conversions
};
 
#endif //ALIANALYSISTASKGAMMA_H
