#ifndef ALIANALYSISTASKB2AOD_H
#define ALIANALYSISTASKB2AOD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for B2 (AOD)
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <AliAnalysisTaskSE.h>

class AliAODEvent;
class AliAODTrack;
class AliAODMCParticle;
class TString;
class TList;
class TProfile;
class AliLnAODtrackCuts;
class AliLnHistoMap;
class AliLnID;

class AliAnalysisTaskB2AOD: public AliAnalysisTaskSE
{
  public:
	AliAnalysisTaskB2AOD();
	AliAnalysisTaskB2AOD(const char* name);
	
	virtual ~AliAnalysisTaskB2AOD();
	
	virtual void UserCreateOutputObjects();
	virtual void UserExec(Option_t* option);
	virtual void Terminate(Option_t *);
	
	void SetV0ANDtrigger(Bool_t flag=1) { fV0AND = flag; }
	void SetNoFastOnlyTrigger(Bool_t flag=1) { fNoFastOnly = flag; }
	void SetNtrkMultTrigger(Bool_t flag=1) { fNtrkMultTrigger = flag; }
	void SetKNOmultInterval(Double_t min, Double_t max) { fMinKNOmult = min; fMaxKNOmult = max; }
	void SetCentralityInterval(Double_t min, Double_t max) { fMinCentrality = min; fMaxCentrality = max; }
	
	void SetMeanNtrk(Double_t mean) { fMeanNtrk = mean; }
	
	void SetVertexXInterval(Double_t min, Double_t max) { fMinVx = min; fMaxVx = max; }
	void SetVertexYInterval(Double_t min, Double_t max) { fMinVy = min; fMaxVy = max; }
	void SetVertexZInterval(Double_t min, Double_t max) { fMinVz = min; fMaxVz = max; }
	
	void SetDCAxyInterval(Double_t min, Double_t max) { fMinDCAxy = min; fMaxDCAxy = max; }
	void SetDCAzInterval(Double_t min, Double_t max) { fMinDCAz = min; fMaxDCAz = max; }
	void SetMaxNSigmaToVertex(Double_t max ) { fMaxNSigma = max; }
	
	void SetEtaInterval(Double_t min, Double_t max) { fMinEta = min; fMaxEta = max; }
	void SetRapidityInterval(Double_t min, Double_t max) { fMinY = min; fMaxY = max; }
	
	void SetSimulation(Bool_t flag = kTRUE) { fSimulation = flag; }
	void SetHeavyIons(Bool_t flag = kTRUE) { fHeavyIons = flag; }
	
	void SetMomentumCorrection(Bool_t flag = kTRUE) { fMomentumCorrection = flag; }
	void SetMomentumCorrectionProfile(TProfile* pfx) { fMoCpfx = pfx; }
	
	void SetHistogramMap(AliLnHistoMap* map) { fHistoMap = map; }
	
	void SetTrackCuts(AliLnAODtrackCuts* trackCuts) { fTrackCuts = trackCuts; }
	
	void SetPID(AliLnID* lnID) { fLnID = lnID; }
	void SetMaxNSigmaITS(Double_t max) { fMaxNSigmaITS = max; }
	void SetMaxNSigmaTPC(Double_t max) { fMaxNSigmaTPC = max; }
	void SetMaxNSigmaTOF(Double_t max) { fMaxNSigmaTOF = max; }
	
	void SetM2Interval(Double_t min, Double_t max) { fMinM2 = min; fMaxM2 = max; };
	
	void SetParticleSpecies(const TString& species);
	
    private:

	AliAnalysisTaskB2AOD(const AliAnalysisTaskB2AOD& other);
	AliAnalysisTaskB2AOD& operator=(const AliAnalysisTaskB2AOD& other);
	
	AliLnHistoMap* CreateHistograms();
	
	Int_t GetParticles();
	Int_t GetTracks();
	
	Bool_t IsV0AND(UInt_t triggerBits) const;
	Bool_t IsFastOnly(UInt_t triggerBits) const;
	Bool_t IsMB(UInt_t triggerBits) const;
	
	AliAODMCParticle* GetParticle(const AliAODTrack* trk) const;
	
	Int_t GetChargedMultiplicity(Double_t maxEta) const;
	
	Bool_t IsFakeTrack(const AliAODTrack* trk) const;
	
	Double_t GetPhi(Double_t p[3]) const;
	Double_t GetRapidity(Double_t p, Double_t pz, Double_t m) const;
	Double_t GetBeta(const AliAODTrack* trk) const;
	Double_t GetMassSquared(Double_t p, Double_t beta) const;
	Int_t    GetITSnClusters(const AliAODTrack* trk) const;
	Int_t    GetITSnPointsPID(const AliAODTrack* trk) const;
	
	Int_t GetPidCode(const TString& species) const;
	
	Double_t GetMomentumCorrection(Double_t ptrec) const;
	
	Double_t GetM2Difference(Double_t beta, Double_t p, Double_t m) const;
	Double_t GetExpectedTime(const AliAODTrack* trk, Double_t m) const;
	
  private:
 
	TString fSpecies; // particle species for the analysis
	Int_t fPartCode; // particle species code
	
	Bool_t fHeavyIons; // analysis of heavy ions data
	Bool_t fSimulation; // analysis of MC simulation
	
	Bool_t fMultTriggerFired; //! track multiplicity trigger flag
	Bool_t fCentTriggerFired; //! centrality trigger flag
	Bool_t fTriggerFired; //! trigger flag
	Bool_t fGoodVertex; //! good vertex flag
	Bool_t fPileUpEvent; //! pile-up flag
	
	Bool_t fV0AND; // V0AND trigger flag
	Bool_t fNoFastOnly; // No kFastOnly trigger flag
	Bool_t fNtrkMultTrigger; // enable combined multiplicity trigger
	Double_t fMinKNOmult; // minimum KNO track multiplicity scaling
	Double_t fMaxKNOmult; // maximum KNO track multiplicity scaling
	Double_t fMinCentrality; // minimum centrality for HI
	Double_t fMaxCentrality; // maximum centrality for HI
	
	Double_t fNch; //! current charged multipicity
	Double_t fNtrk; //! current track multipicity
	Double_t fMeanNtrk; // average track multiplicity
	Double_t fKNOmult; //! KNO track multiplicity scaling
	
	Double_t fMinVx; // vertex low X value
	Double_t fMaxVx; // vertex high X value
	Double_t fMinVy; // vertex low Y value
	Double_t fMaxVy; // vertex high Y value
	Double_t fMinVz; // vertex low Z value
	Double_t fMaxVz; // vertex high Z value
	
	Double_t fMinDCAxy; // minimum DCAxy
	Double_t fMaxDCAxy; // maximum DCAxy
	Double_t fMinDCAz; // minimum DCAz
	Double_t fMaxDCAz; // maximum DCAz
	Double_t fMaxNSigma; // maximum number of sigmas to primary vertex
	
	Double_t fMinEta; // minimum pseudorapidity
	Double_t fMaxEta; // maximum pseudorapidity
	Double_t fMinY; // minimum rapidity
	Double_t fMaxY; // maximum rapidity
	
	AliAODEvent* fAODevent; //! AOD event
	
	TList* fOutputContainer; // output container
	AliLnHistoMap*  fHistoMap; // histogram map (defined somewhere else)
	AliLnAODtrackCuts* fTrackCuts; // track cuts (defined somewhere else)
	AliLnID* fLnID; // PID for light nuclei (defined somewhere else)
	
	Double_t fMaxNSigmaITS; // maximum number of sigmas to dEdx in the ITS
	Double_t fMaxNSigmaTPC; // maximum number of sigmas to dEdx in the TPC
	Double_t fMaxNSigmaTOF; // maximum number of sigmas to dEdx in the TOF
	
	Double_t fMinM2; // minimum m2 for TPC+TOF pid
	Double_t fMaxM2; // maximum m2 for TPC+TOF pid
	
	Bool_t fMomentumCorrection; // enable momentum correction
	TProfile* fMoCpfx; // momentum correction from simulation
	
	ClassDef(AliAnalysisTaskB2AOD, 1)
};

#endif // ALIANALYSISTASKB2AOD_H
