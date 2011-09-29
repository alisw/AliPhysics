#ifndef ALIANALYSISTASKGAMMACONVDALITZ_H
#define ALIANALYSISTASKGAMMACONVDALITZ_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for pi0->e+e-gamma (Dalitz decay)

#include "AliAnalysisTaskSE.h"

class AliESDInputHandler;
class AliMCEventHandler;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDpidCuts;
class AliV0Reader;
class AliGammaConversionHistograms;

class AliAnalysisTaskGammaConvDalitz: public AliAnalysisTaskSE
{
  public:

	AliAnalysisTaskGammaConvDalitz();
	AliAnalysisTaskGammaConvDalitz( const char* name );
	virtual ~AliAnalysisTaskGammaConvDalitz();

	virtual void UserExec(Option_t *option);
	virtual void UserCreateOutputObjects();
	virtual void ConnectInputData(Option_t *option);
	virtual void Terminate(Option_t *option);
	
	enum TrackSelectionCriteria { kITSsaTrack=0, kGlobalTrack=1, kITSsaGlobalTrack=2 };
	
	void SetRunStandalone( Bool_t flag=kFALSE ) { fStandalone = flag; }
	void SetComputeBackground( Bool_t flag=kTRUE ) { fComputeBkg = flag; }
	void SetUseBayesPID( Bool_t flag=kTRUE ) { fUseBayesPID = flag; }
	void SetUseESDtrackIndexCut( Bool_t flag=kTRUE) { fUseTrackIndexCut = flag; }
	void SetUsePsiPairCut(Bool_t flag=kTRUE) { fUsePsiPairCut = flag; }
	void SetUseMassCut(Bool_t flag=kTRUE) { fUseMassCut = flag; }
	void SetUseGammaCut(Bool_t flag=kTRUE) { fUseGammaCut = flag; }
	void SetUseAliKF(Bool_t flag=kFALSE) { fUseAliKF = flag; }
	void SetTrackSelectionCriteria(AliAnalysisTaskGammaConvDalitz::TrackSelectionCriteria sel=kGlobalTrack) { fTrkSelectionCriteria = sel; }
	
	void SetPsiPairCut(Double_t psi=0.45, Double_t phiMin=0., Double_t phiMax=0.12, Bool_t readMagFieldSgn=kTRUE){fPsiPairCut = psi; fDeltaPhiCutMin = phiMin; fDeltaPhiCutMax = phiMax; fReadMagFieldSign = readMagFieldSgn;}
	void SetMassCut(Double_t min, Double_t max) {fMassCutMin = min; fMassCutMax = max; }

	void SetNSigmasElecTPC( Double_t min, Double_t max) { fNSigmaBelowElecTPCbethe = min; fNSigmaAboveElecTPCbethe = max; }
	void SetNSigmasPionTPC( Double_t max )   {  fNSigmaAbovePionTPCbethe = max; }
	void SetNSigmasKaonTPC( Double_t max )   {  fNSigmaAboveKaonTPCbethe	= max; }
	void SetNSigmasProtonTPC( Double_t max ) {  fNSigmaAboveProtonTPCbethe  = max; }
	
	void SetV0Reader( AliV0Reader* reader ) { fV0Reader = reader; }
	void SetDoMC(Bool_t flag) { fDoMC = flag; }
	void SetBGHandler( AliGammaConversionBGHandler* BG ) {  fBGEventHandler = BG; }
	
	void AdoptHistograms( AliGammaConversionHistograms* histograms ) { fHistograms = histograms; }
	void AdoptITSsaTrackCuts( AliESDtrackCuts* esdCuts = 0 );
	void AdoptESDtrackCuts( AliESDtrackCuts* esdCuts = 0 );
	void AdoptESDpidCuts( AliESDpidCuts* esdPIDCuts = 0 );
	
  private:

	void ProcessMCData();
	void CreateListOfDalitzPairCandidates();
	void ProcessGammaElectronsForDalitzAnalysis();
	
	void ESDtrackIndexCut(vector<Int_t>& pos, vector<Int_t>& neg, const TClonesArray* gamma);
	void PsiPairCut(vector<Int_t>& pos, vector<Int_t>& neg);
	void MassCut(vector<Int_t>& pos, vector<Int_t>& neg);
	void CleanArray(vector<Int_t>& x, const vector<Bool_t>& tag);
	
	TClonesArray* IndexToAliKFParticle(const vector<Int_t>& v, Int_t PDG);
	TClonesArray* FindElectronFromPi0Dalitz(const vector<Int_t>& candidates, const Int_t PDG);
	TClonesArray* FindGammaFromPi0Dalitz(const TClonesArray* candidates, const vector<Int_t>& pos, const vector<Int_t>& neg);
	TClonesArray* FindGamma(const TClonesArray* candidates, const vector<Int_t>& pos, const vector<Int_t>& neg);
	TClonesArray* FindDalitzPair(const TClonesArray* pos, const TClonesArray* neg);
	TClonesArray* FindDalitzPair(const vector<Int_t>& pos, const vector<Int_t>& neg,Int_t motherOpc);
        TClonesArray* FindJpsi(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx,Int_t motherOpc);
	TClonesArray* FindParticleDalitz(const TClonesArray* pos, const TClonesArray* neg, const TClonesArray* gamma,Int_t opc);
	TClonesArray* FindParticleDalitz(const vector<Int_t>& pos, const vector<Int_t>& neg, const TClonesArray* gamma, const vector<Int_t>& posGam, const vector<Int_t>& negGam,Int_t motherOpc);
        TClonesArray* FindParticleChic(const vector<Int_t>& posIdx, const vector<Int_t>& negIdx, const TClonesArray* gamma, const vector<Int_t>& posGam, const vector<Int_t>& negGam,Int_t motherOpc);

	void SetGammaPoolMaxSize(UInt_t maxSize=10) { fPoolMaxSize = maxSize; }
	void UpdateGammaPool(const TClonesArray* gamma);
	void UpdateElectronPool(TClonesArray* elec);
	TClonesArray* GammasFromBGHandler() const;
	TClonesArray* ElectronFromBGHandler() const;

	Bool_t IsPi0DalitzDaughter( Int_t label ) const;
	Bool_t IsDalitzPair( Int_t labelPos, Int_t labelNeg, Int_t motherOpc ) const;
	Bool_t IsFromGammaConversion( Int_t labelPos, Int_t labelNeg  ) const;
	Bool_t IsFromGammaConversion( Double_t psiPair, Double_t deltaPhi ) const;
	Bool_t HaveSameMother( Int_t label1, Int_t label2 ) const;
	
	Double_t GetPsiPair( const AliKFParticle* pos, const AliKFParticle* neg ) const;
	Double_t GetPsiPair( const TLorentzVector* pos, const TLorentzVector* neg ) const;
	Double_t GetPsiPair( const AliESDtrack* trackPos, const AliESDtrack* trackNeg ) const;
	
	Int_t GetMonteCarloPid(const AliESDtrack* t) const;
	Int_t GetBayesPid(const AliESDtrack* t, Int_t trackType ) const;
	Int_t GetNSigmaPid(const AliESDtrack* t, Int_t trackType ) const;
	
	void GetGammaCandidates(TClonesArray*& gamma, vector<Int_t>& posIndex, vector<Int_t>& negIndex);
	void AngleEposEnegGammaCut( const vector<Int_t>& xPosIndex, const vector<Int_t>& yNegIndex, const TClonesArray* zGamma, TClonesArray*& gamma, vector<Int_t>& posIndex, vector<Int_t>& negIndex);
	void FillPsiPair(const TClonesArray* pos, const TClonesArray* neg, const TString& hName);
	void FillAngle(const TClonesArray* x, const TClonesArray* y, const TString& hName);
	Double_t Rapidity(const TParticle* p) const;
	void FillPidTable(const TParticle* p, Int_t pid);

	
  // protected:
  private:

	AliStack*   fStack;                 //! MC particle stack
	AliMCEvent* fGCMCEvent;               //! for CF pointer to the MC Event

	AliESDEvent* fESDEvent;             //! ESD event

	vector<Int_t> fEposCandidateIndex;  //! track indexes of e+ candidates
	vector<Int_t> fEnegCandidateIndex;  //! track indexes of e- candidates
	vector<Int_t> fGammaCandidatePosIndex; //! track indexes for gamma candidates positive track
	vector<Int_t> fGammaCandidateNegIndex; //! track indexes for gamma candidates negative track

	TClonesArray* fGammaCandidates;     //! AliKFParticle gamma candidates
	TClonesArray* fGammaPool;           //! AliKFParticle gamma pool of previous events
	Int_t fPoolMaxSize; // size of the gamma pool
	Int_t fGamPoolPos; // Posisiton of last added gamma in the pool
	
	AliGammaConversionBGHandler* fBGEventHandler; // Background event handler

	TList* fOutputContainer;           // Histogram container
	AliMCEventHandler* fMCTruth;       // for CF pointer to MCTruth
	AliV0Reader* fV0Reader;            // The V0 reader object
	AliESDpid* fESDpid;                // for dEdx cut based on nSigma to a particle line
	AliESDtrackCuts* fESDtrackCuts;    // ESD global track cuts
	AliESDtrackCuts* fITSsaTrackCuts;  // ITS standalone ESD track cuts
	AliESDpidCuts* fESDpidCuts;        // ESD PID cuts
	AliGammaConversionHistograms* fHistograms;  // histogram container
	
	Bool_t fStandalone;        // Run the task as standalone for the V0reader
	Bool_t fDoMC;              // process montecarlo simulation
	Bool_t fComputeBkg;        // Compute combinatorial background
	Bool_t fUseBayesPID;       // use bayesian pid
	Bool_t fUseTrackIndexCut;  // use esd track index cut
	Bool_t fUsePsiPairCut;     // use psi pair cut
	Bool_t fUseMassCut;        // use mass cut
	Bool_t fUseGammaCut;       // use e+e- plane angle gamma cut
	Bool_t fReadMagFieldSign;  // Read the magnetic field sign from the ESD for Psi pair cut
	Bool_t fUseAliKF;          // Use AliKFParticle to reconstruct the pi0 instead of TLorentzVector class
	
	Int_t fMagFieldSign;                 // Magnetic field sign
	const Double_t fkElectronMass;        // Electron mass
	Double_t fPsiPairCut;                // Psi pair cut value
	Double_t fDeltaPhiCutMin;              // Delta_Phi minimum cut value
	Double_t fDeltaPhiCutMax;              // Delta_Phi maximum cut value
	Double_t fMassCutMin;                // Minimum value of e+e- mass (GeV/c2)
	Double_t fMassCutMax;                // Maximum value of e+e- mass (GeV/c2)
	Double_t fNSigmaBelowElecTPCbethe;   // Number of sigmas below the electron BB in the TPC
	Double_t fNSigmaAboveElecTPCbethe;   // Number of sigmas above the electron BB in the TPC
	Double_t fNSigmaAbovePionTPCbethe;   // Number of sigmas above the Pion BB in the TPC
	Double_t fNSigmaAboveKaonTPCbethe;   // Number of sigmas above the Kaon BB in the TPC
	Double_t fNSigmaAboveProtonTPCbethe; // Number of sigmas above the Proton BB in the TPC

	TrackSelectionCriteria fTrkSelectionCriteria;      // Selected criteria for track cuts

 private:
 
	AliAnalysisTaskGammaConvDalitz( const AliAnalysisTaskGammaConvDalitz& ); // Not implemented
	AliAnalysisTaskGammaConvDalitz& operator=( const AliAnalysisTaskGammaConvDalitz& ); // Not implemented
	
	ClassDef( AliAnalysisTaskGammaConvDalitz, 2 );
};

#endif // ALIANALYSISTASKGAMMACONVDALITZ_H
