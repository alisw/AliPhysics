#ifndef ALILNPT_H
#define ALILNPT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// reconstructed pt
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

#ifndef HAVE_ROOUNFOLD
//#define HAVE_ROOUNFOLD
#endif

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
	
	TH1D* PID(const TH1D* hPt, const TH2D* hM2Pt, Int_t lowbin, Int_t hibin, const TString& name);
#ifdef HAVE_ROOUNFOLD
	TH1D* Unfolding(const TH1D* hPt, const TH1D* hMeasuredPt, const TH1D* hTruePt, const TH2D* ResponseMtx, const TString& name, Int_t iterations=4) const;
#endif
	TH1D* Secondaries(const TH1D* hPt, const TH1D* hFracMatPt, const TH1D* hFracFdwnPt, const TString& name) const;
	TH1D* Secondaries(const TH1D* hPt, TF1* fncMatPt, TF1* fncFdwnPt, const TString& name) const;
	TH1D* Efficiency(const TH1D* hPt, const TH1D* hEffVtxPt, const TH1D* hEffAccTrkPt, const TString& name) const;
	
	Double_t GetVertexCorrection(const TH1D* hData, const TH1D* hMC) const;
	
	Int_t Exec();
	
	void SetRapidityInterval(Double_t ymin, Double_t ymax) { fYMin = ymin; fYMax = ymax; }
	
	void SetParticle(const TString& particle) { fParticle= particle; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	void SetCorrectionTag(const TString& tag) { fCorrTag = tag; }
	
	void SetPtBinInterval(Int_t lowbin, Int_t hibin) { fLowPtBin = lowbin; fHiPtBin = hibin; };
	void SetM2BinInterval(Int_t lowbin, Int_t hibin) { fLowM2Bin = lowbin; fHiM2Bin = hibin; };
	
	void SetM2BkgInterval(Double_t min, Double_t max) { fMinM2Bkg = min; fMaxM2Bkg = max; };
	void SetM2TPCInterval(Double_t min, Double_t max) { fMinM2tpc = min; fMaxM2tpc = max; };
	
	void SetPidM2(Bool_t flag=1) { fPidM2 = flag; }
	void SetUnfolding(Bool_t flag=1, Int_t niter=4) { fUnfolding = flag; fNIter=niter; }
	void SetSecondaries(Bool_t flag=1) { fSecondaries = flag; }
	void SetEfficiency(Bool_t flag=1) { fEfficiency = flag; }
	
	void SetOnlyGeneration(Bool_t flag=1) { fIsOnlyGen = flag; }
	
	void SetMakeStats(Bool_t flag=1) { fMakeStats = flag; }
	
	void SetMCtoINEL(Bool_t flag=1) { fMCtoINEL=flag;}
	void SetFitFractionCorr(Bool_t flag=1) { fFitFrac=flag; }
	void SetFeedDownCorr(Bool_t flag=1) { fFdwnCorr=flag; }
	void SetSameFeedDownCorr(Bool_t flag=1) { fSameFdwn=flag; }
	
	
  private:
 
	AliLnPt(const AliLnPt& other);
	AliLnPt& operator=(const AliLnPt& other);
	
	void GetPtFromM2(TH1D* hPt, Int_t lowbin, Int_t hibin, const TH2D* hM2Pt, Double_t m2min, Double_t m2max) const;
	Double_t GetM2Width(Double_t pt, Double_t mass) const;
	RooWorkspace* GetM2Model(Double_t pt, Double_t m, const TString& name, Double_t max) const;
	
  private:
	
	TString fParticle; // particle name
	
	Double_t fTrigEff; // trigger efficiency
	
	TString fInputFilename; // input filename
	
	TString fOutputFilename; // output for the result
	TString fOutputTag; // tag for the ouput
	
	TString fCorrFilename; // file with the corrections
	TString fCorrTag; // tag for the correction file
	
	Int_t fLowPtBin; // low bin for the pt in GeV/c/A
	Int_t fHiPtBin; // high bin for the pt in GeV/c/A
	
	Bool_t fPidM2; // correct contamination of m2 pid
	Bool_t fUnfolding; // unfold the pt correction
	Int_t  fNIter; // number of iterations for Bayesian unfolding
	Bool_t fSecondaries; // remove secondaries
	Bool_t fEfficiency; // correct by efficiency
	
	Bool_t fIsOnlyGen; // if the simulation is only generation
	
	Double_t fYMin; // rapidity interval min y
	Double_t fYMax; // rapidity interval max y
	
	Int_t fLowM2Bin; // low pt bin for pid correction
	Int_t fHiM2Bin; // high pt bin for pid correction
	
	Double_t fMinM2Bkg; // minimum m2 to remove background
	Double_t fMaxM2Bkg; // maximum m2 to remove background
	
	Double_t fMinM2tpc; // minimum m2 for integration
	Double_t fMaxM2tpc; // maximum m2 for integration
	
	Bool_t fMakeStats; // make event stats
	
	Bool_t fMCtoINEL; // MC to extrapolate to inel or for triggering events
	Double_t fVtxFactor; // vertex correction factor between simulation and data
	
	Bool_t fFitFrac; // fit for fraction of secondaries
	Bool_t fFdwnCorr; // enable feed-down correction
	Bool_t fSameFdwn; // same feed-down correction for positives and negatives
	
	ClassDef(AliLnPt,1)
};

#endif // ALILNPT_H
