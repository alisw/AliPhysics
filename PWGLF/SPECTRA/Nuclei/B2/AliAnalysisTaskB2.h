#ifndef ALIANALYSISTASKB2_H
#define ALIANALYSISTASKB2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis task for B2
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <AliAnalysisTask.h>
#include <AliPID.h>

class AliESDtrack;
class AliMCEvent;
class AliESDEvent;
class TString;
class AliESDtrackCuts;
class AliLnHistoMap;
class AliLnID;
class TParticle;
class TList;

class AliAnalysisTaskB2: public AliAnalysisTask
{
  public:
	AliAnalysisTaskB2();
	AliAnalysisTaskB2(const char* name);
	
	virtual ~AliAnalysisTaskB2();
	
	virtual void ConnectInputData(Option_t *);
	virtual void CreateOutputObjects();
	virtual void Exec(Option_t* option);
	virtual void Terminate(Option_t *);
	
	void SetV0ANDtrigger(Bool_t flag=1) { fV0AND = flag; }
	void SetNoFastOnlyTrigger(Bool_t flag=1) { fNoFastOnly = flag; }
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
	
	void SetHistogramMap(AliLnHistoMap* map) { fHistoMap = map; }
	
	void SetESDtrackCuts(AliESDtrackCuts* esdTrackCuts) {fESDtrackCuts = esdTrackCuts; }
	
	void SetPID(AliLnID* lnID) { fLnID = lnID; }
	void SetMaxNSigmaITS(Double_t max) { fMaxNSigmaITS = max; }
	void SetMaxNSigmaTPC(Double_t max) { fMaxNSigmaTPC = max; }
	void SetMaxNSigmaTOF(Double_t max) { fMaxNSigmaTOF = max; }
	
	void SetTOFmatch(Bool_t flag = kTRUE) { fTOFmatch = flag; }
	Bool_t AcceptTOFtrack(const AliESDtrack* trk) const;
	
	void SetM2Interval(Double_t min, Double_t max) { fMinM2 = min; fMaxM2 = max; };
	
	void SetParticleSpecies(const TString& species);
	
    private:

	AliAnalysisTaskB2(const AliAnalysisTaskB2& other);
	AliAnalysisTaskB2& operator=(const AliAnalysisTaskB2& other);
	
	AliLnHistoMap* CreateHistograms();
	
	Int_t GetParticles();
	Int_t GetTracks();
	
	Bool_t IsV0AND() const;
	Bool_t IsFastOnly(UInt_t triggerBits) const;
	Bool_t IsMB(UInt_t triggerBits) const;
	
	Bool_t TOFmatch(const AliESDtrack* trk) const;
	
	TParticle* GetParticle(const AliESDtrack* trk) const;
	
	Int_t GetChargedMultiplicity(Double_t maxEta) const;
	
	Bool_t IsFakeTrack(const AliESDtrack* trk) const;
	Bool_t IsPhysicalPrimary(const TParticle* prt) const;
	Bool_t IsFromMaterial(const TParticle* prt) const;
	Bool_t IsFromWeakDecay(const TParticle* prt) const;
	
	Double_t GetSign(TParticle* prt) const;
	
	Double_t GetPhi(const AliESDtrack* trk) const;
	Double_t GetTheta(const AliESDtrack* trk) const;
	Double_t GetRapidity(const AliESDtrack* trk, Int_t pid) const;
	Double_t GetITSmomentum(const AliESDtrack* trk) const;
	Double_t GetTOFmomentum(const AliESDtrack* trk) const;
	Double_t GetBeta(const AliESDtrack* trk) const;
	Double_t GetMassSquare(const AliESDtrack* trk) const;
	Double_t GetTimeOfFlight(const AliESDtrack* trk) const;
	Double_t GetITSchi2PerCluster(const AliESDtrack* trk) const;
	Int_t    GetITSnClusters(const AliESDtrack* trk) const;
	Int_t    GetITSnPointsPID(const AliESDtrack* trk) const;
	
	Int_t GetPidCode(const TString& species) const;
	
  private:
 
	TString fSpecies; // particle species for the analysis
	Int_t fPartCode; // particle species code
	
	Bool_t fHeavyIons; // analysis of heavy ions data
	Bool_t fSimulation; // analysis of MC simulation
	
	Bool_t fMultTrigger; //! track multiplicity trigger flag
	Bool_t fCentTrigger; //! centrality trigger flag
	Bool_t fTriggerFired; //! trigger flag
	Bool_t fGoodVertex; //! good vertex flag
	Bool_t fPileUpEvent; //! pile-up flag
	
	Bool_t fV0AND; // V0AND trigger flag
	Bool_t fNoFastOnly; // No kFastOnly trigger flag
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
	
	AliMCEvent* fMCevent;   //! monte carlo event
	AliESDEvent* fESDevent; //! ESD event
	
	TList* fOutputContainer; // output container
	AliLnHistoMap*  fHistoMap; // histogram map (defined somewhere else)
	AliESDtrackCuts* fESDtrackCuts; // track cuts (defined somewhere else)
	AliLnID* fLnID; // PID for light nuclei (defined somewhere else)
	
	Double_t fMaxNSigmaITS; // maximum number of sigmas to dEdx in the ITS
	Double_t fMaxNSigmaTPC; // maximum number of sigmas to dEdx in the TPC
	Double_t fMaxNSigmaTOF; // maximum number of sigmas to dEdx in the TOF
	
	class AliTriggerAnalysis* fTrigAna; //! to access trigger information
	class AliESDpid* fESDpid; //! ESD pid
	Bool_t fIsPidOwner; // whether we own fESDpid
	Int_t fTimeZeroType;  // time zero type
	
	Bool_t fTOFmatch; // TOF match flag
	
	Double_t fMinM2; // minimum m2 for TPC+TOF pid
	Double_t fMaxM2; // maximum m2 for TPC+TOF pid
	
	ClassDef(AliAnalysisTaskB2, 1)
};

#endif // ALIANALYSISTASKB2_H
