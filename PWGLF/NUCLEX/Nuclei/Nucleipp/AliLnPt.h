#ifndef ALILNPT_H
#define ALILNPT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// reconstructed pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;
class TH1D;
class TH2D;
class RooWorkspace;
class TF1;

class AliLnPt: public TObject
{
  public:
	
	AliLnPt(const TString& particle, Double_t trigEff, const TString& inputFilename, const TString& outputFilename, const TString& otag, const TString& corrFilename, const TString& corrtag);
	virtual ~AliLnPt();
	
	TH1D* PID(const TH1D* hPt, const TH2D* hPidPt, Double_t ptmin, Double_t ptmax, const TString& name);
	TH1D* Secondaries(const TH1D* hPt, const TH1D* hFracMatPt, const TH1D* hFracFdwnPt, const TString& name) const;
	TH1D* Secondaries(const TH1D* hPt, TF1* fncMatPt, TF1* fncFdwnPt, const TString& name) const;

	TH1D* Efficiency(const TH1D* hPt, const TH1D* hEffVtxPt, const TH1D* hEffAccTrkPt, const TString& name) const;
	
	Double_t GetVertexCorrection(const TH1D* hData, const TH1D* hMC) const;
	
	Int_t Exec();
	
	void SetParticle(const TString& particle) { fParticle= particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	void SetCorrectionTag(const TString& tag) { fCorrTag = tag; }
	
	void SetPtInterval(Double_t min, Double_t max) { fPtMin = min; fPtMax = max; };
	
	void SetPidProcedure(Int_t proc) { fPidProc=proc; }
	void SetPidPt(Double_t ptpid) { fPtPid = ptpid; };
	void SetBkgInterval(Double_t min, Double_t max) { fBkgMin = min; fBkgMax = max; };
	void SetPidInterval(Double_t min, Double_t max) { fIntMin = min; fIntMax = max; };
	
	void SetPid(Bool_t flag=1) { fPid = flag; }
	void SetSecondaries(Bool_t flag=1) { fSecondaries = flag; }
	void SetEfficiency(Bool_t flag=1) { fEfficiency = flag; }
	
	void SetOnlyGeneration(Bool_t flag=1) { fIsOnlyGen = flag; }
	
	void SetMakeStats(Bool_t flag=1) { fMakeStats = flag; }
	
	void SetMCtoINEL(Bool_t flag=1) { fMCtoINEL=flag;}
	void SetFitFractionCorr(Bool_t flag=1) { fFitFrac=flag; }
	void SetFeedDownCorr(Bool_t flag=1) { fFdwnCorr=flag; }
	void SetSameFeedDownCorr(Bool_t flag=1) { fSameFdwn=flag; }
	
	void SetPidEfficiency(Double_t eff) { fPidEff=eff; }
	
	void SetDebugLevel(Int_t level) { fDebugLevel = level; }
	
	enum { kMassSquared=0, kMassDifference, kTimeOfFlight };
	
  private:
 
	AliLnPt(const AliLnPt& other);
	AliLnPt& operator=(const AliLnPt& other);
	
	void GetPtFromPid(TH1D* hPt, Double_t ptmin, Double_t ptmax, const TH2D* hM2Pt, Double_t m2min, Double_t m2max) const;
	Double_t GetM2Width(Double_t pt, Double_t mass) const;
	RooWorkspace* GetM2Model(Double_t pt, Double_t m, const TString& name, Double_t max) const;
	RooWorkspace* GetToFModel(Double_t pt, const TString& name) const;
	Double_t GetExpectedDT(Double_t pt, Double_t m) const;
	
  private:
	
	TString fParticle; // particle name
	
	Double_t fTrigEff; // trigger efficiency
	
	TString fInputFilename; // input filename
	
	TString fOutputFilename; // output for the result
	TString fOutputTag; // tag for the ouput
	
	TString fCorrFilename; // file with the corrections
	TString fCorrTag; // tag for the correction file
	
	Double_t fPtMin; // minimum pt value in GeV/c
	Double_t fPtMax; // maximum pt value in GeV/c
	
	Bool_t fPid; // enable additional PID
	Bool_t fSecondaries; // remove secondaries
	Bool_t fEfficiency; // correct by efficiency
	
	Bool_t fIsOnlyGen; // if the simulation is only generation
	
	Double_t fPtPid; // minimum pt value for pid correction
	
	Double_t fBkgMin; // minimum value to remove background
	Double_t fBkgMax; // maximum value to remove background
	
	Double_t fIntMin; // minimum value for integration
	Double_t fIntMax; // maximum value for integration
	
	Bool_t fMakeStats; // make event stats
	
	Bool_t fMCtoINEL; // MC to extrapolate to inel or for triggering events
	Double_t fVtxFactor; // vertex correction factor between simulation and data
	
	Bool_t fFitFrac; // fit for fraction of secondaries
	Bool_t fFdwnCorr; // enable feed-down correction
	Bool_t fSameFdwn; // same feed-down correction for positives and negatives
	
	Int_t fPidProc;  // selected pid procedure on the pt distribution
	
	Double_t fPidEff;  // pid efficiency for all pt bins
	
	Int_t fDebugLevel; // 0 no verbose, > 1 verbose
	
	ClassDef(AliLnPt,1)
};

#endif // ALILNPT_H
