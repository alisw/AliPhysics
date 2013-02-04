#ifndef ALILNDRIVER_H
#define ALILNDRIVER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// driver for computing the pt and spectra
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#include <TObject.h>
#include <TString.h>

class TString;
class TFileMerger;

class AliLnDriver: public TObject
{
  public:
	
	AliLnDriver();
	virtual ~AliLnDriver();
	
	void MakePtCorr() const;
	void MakePt() const;
	void MakeRatio() const;
	void MakeSpectra() const;
	
	void PrintFilenames() const;
	
	Int_t Run() const;
	
	TString GetSpecies() const { return fSpecies; }
	TString GetOutputTag() const { return fOutputTag; }
	TString GetOutputCorrTag() const { return fOutputCorTag; }
	
	TString GetPtCorrDebugFilename() const { return fOutputPtCorrDebug; }
	TString GetPtDebugFilename() const { return fOutputPtDebug; }
	
	void SetInputFilenames(const TString& data, const TString& simu, const TString& simuFix, const TString& ptcorr);
	void SetOutputFilenames(const TString& pt, const TString& ratio, const TString& spectra);
	void SetDebugFilenames(const TString& corrdebug, const TString& ptdebug) { fOutputPtCorrDebug = corrdebug; fOutputPtDebug = ptdebug; }
	
	void SetSpecies(const TString& species) { fSpecies = species; }
	
	void SetOutputTag(const TString& tag) { fOutputTag = tag; }
	void SetOutputCorTag(const TString& tag) { fOutputCorTag = tag; }
	
	void SetTriggerEfficiency(Double_t eff[3]);
	void SetInelXSection(Double_t xsec[3]);
	void SetExtrapolateToINEL(Bool_t flag=1) { fINEL = flag; }
	void SetOnlyGeneration(Bool_t flag=1) { fIsOnlyGen = flag; }
	
	void SetMakeCorrections(Bool_t flag=1) { fMakeCorr = flag; }
	void SetMakePt(Bool_t flag=1) { fMakePt = flag; }
	void SetMakeRatio(Bool_t flag=1) { fMakeRatio = flag; }
	void SetMakeSpectra(Bool_t flag=1) { fMakeSpectra = flag; }
	void SetMakeStats(Bool_t flag=1) { fMakeStats = flag; }
	
	void SetRapidityInterval(Double_t ymin, Double_t ymax) { fYMin = ymin; fYMax = ymax; }
	
	void SetPtBinInterval(Int_t lowbin, Int_t hibin) { fLowPtBin = lowbin; fHiPtBin = hibin; }
	void SetM2BinInterval(Int_t lowbin, Int_t hibin) { fLowM2Bin = lowbin; fHighM2Bin = hibin; }
	void SetM2BkgInterval(Double_t min, Double_t max) { fMinM2Bkg = min; fMaxM2Bkg = max; }
	void SetM2TPCInterval(Double_t min, Double_t max) { fMinM2tpc = min; fMaxM2tpc = max; }
	void SetPidM2(Bool_t flag=1) { fPidM2 = flag; }
	void SetUnfolding(Bool_t flag=1, Int_t niter=4) { fUnfolding = flag; fNIter=niter; }
	void SetSecondaries(Bool_t flag=1) { fSecondaries = flag; }
	void SetSecProd(Int_t prod) { fSecProd = prod; }
	void SetAntiNucleusAsTemplate(Bool_t flag=1) { fANucTemplate = flag; }
	void SetMatDCAxyModel(Int_t model=1) { fMatDCAxyMod = model; }
	void SetNBin(Int_t nbin) { fNbin = nbin; }
	void SetDCAxyInterval(Double_t dcamin, Double_t dcamax) { fMinDCAxy = dcamin; fMaxDCAxy = dcamax; }
	void SetEfficiency(Bool_t flag=1, Bool_t g3Fluka=0) { fEfficiency = flag; fG3Fluka=g3Fluka; }
	void SetScalingFactors(Double_t mat, Double_t fd) { fScMat=mat; fScFd=fd; }
	
	void SetMCtoINEL(Bool_t flag=1) { fMCtoINEL = flag; }
	void SetFitFractionCorr(Bool_t flag=1) { fFitFrac=flag; }
	void SetFeedDownCorr(Bool_t flag=1) { fFdwnCorr=flag; }
	void SetSameFeedDownCorr(Bool_t flag=1) { fSameFdwn = flag; }
	
	void SetAddFakeTracks(Bool_t flag=1) { fAddFakeTracks = flag; }
	
	void SetSysErr( Double_t pos, Double_t neg) { fSysPos = pos; fSysNeg = neg; }
	
  private:
 
	AliLnDriver(const AliLnDriver& other);
	AliLnDriver& operator=(const AliLnDriver& other);
	
  private:
 
	TString fSpecies;       // particle species
	TString fOutputTag;     // tag for output file
	TString fOutputCorTag;  // tag for correction file
	
	Double_t fTrigEff[3];   // trigger efficiency, stat. and syst. errors
	Double_t fXsec[3];      // total inelastic cross section, stat. and syst. errors
	Bool_t  fIsOnlyGen;     // if it is only generation
	Bool_t  fINEL;          // extrapolate to inelastic events
	
	Bool_t  fMakeCorr;      // make corrections
	Bool_t  fMakePt;        // make pt
	Bool_t  fMakeRatio;     // make antiparticle/particle ratio
	Bool_t  fMakeSpectra;   // make spectra
	Bool_t  fMakeStats;     // make event stats
	
	Int_t    fLowPtBin;     // low pt bin
	Int_t    fHiPtBin;      // high pt bin
	Bool_t   fPidM2;        // enable m2 pid correction
	Int_t    fLowM2Bin;     // low m2 bin for pid contamination
	Int_t    fHighM2Bin;    // high m2 bin for pid contamination
	Bool_t   fUnfolding;    // unfolding correction
	Int_t    fNIter;        // number of iterations for Bayesian unfolding
	Bool_t   fSecondaries;  // correction of secondaries
	Int_t    fSecProd;      // procedure for estimating fractions
	Int_t    fMatDCAxyMod;  // DCAxy model for correction of secondaries
	Bool_t   fANucTemplate; // enable antinucleus as template for primaries
	Int_t    fNbin;         // rebin of DCAxy distribution
	Double_t fYMin;         // min rapidity
	Double_t fYMax;         // max rapidity
	Double_t fMinDCAxy;     // min DCAxy
	Double_t fMaxDCAxy;     // max DCAxy
	Double_t fMinM2Bkg;     // min M2 for removing background
	Double_t fMaxM2Bkg;     // max M2 for removing background
	Double_t fMinM2tpc;     // min M2 for integration
	Double_t fMaxM2tpc;     // max M2 for integration
	Bool_t   fEfficiency;   // efficiency correction
	Bool_t   fG3Fluka;      // enable G3/Fluka correction for TPC
	Double_t fScMat;        // scaling factor for material fraction
	Double_t fScFd;         // scaling factor for feed-down fraction
	
	Double_t fSysPos;       // variation for positives
	Double_t fSysNeg;       // variation for negatives
	
	TString fInputData;     // input data filename
	TString fInputSimu;     // input simulation filename
	TString fInputSimuFix;  // input fixed simulation filename
	
	TString fOutputPtCorr;  // output filename for pt corrections
	TString fOutputPt;      // output filename for pt
	TString fOutputRatio;   // output filename for antiparticle/particle ratio
	TString fOutputSpectra; // output filename for differential yields
	
	TString fOutputPtCorrDebug; // output filename for debugging pt corrections
	TString fOutputPtDebug;     // output filename for debugging pt
	
	Bool_t fFitFrac;        // fit for fraction of secondaries
	Bool_t fFdwnCorr;       // enable feed-down correction
	Bool_t fSameFdwn;       // same feed-down correction for positives and negatives
	Bool_t fMCtoINEL;       // MC to extrapolate to inel or for triggering events
	
	Bool_t fAddFakeTracks;  // include fake tracks in the efficiency and templates
	
	ClassDef(AliLnDriver,3)
};

#endif // ALILNDRIVER_H
