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
	void SetExtrapolateToINEL(Bool_t flag=1) { fINEL = flag; }
	void SetOnlyGeneration(Bool_t flag=1) { fIsOnlyGen = flag; }
	
	void SetMakeCorrections(Bool_t flag=1) { fMakeCorr = flag; }
	void SetMakePt(Bool_t flag=1) { fMakePt = flag; }
	void SetMakeRatio(Bool_t flag=1) { fMakeRatio = flag; }
	void SetMakeSpectra(Bool_t flag=1) { fMakeSpectra = flag; }
	void SetMakeStats(Bool_t flag=1) { fMakeStats = flag; }
	
	void SetRapidityInterval(Double_t ymin, Double_t ymax) { fYMin = ymin; fYMax = ymax; }
	
	void SetPtInterval(Double_t min, Double_t max) { fPtMin = min; fPtMax = max; }
	void SetPidPt(Double_t ptpid) { fPidPt = ptpid; }
	void SetBkgInterval(Double_t min, Double_t max) { fBkgMin = min; fBkgMax = max; }
	void SetPidInterval(Double_t min, Double_t max) { fIntMin = min; fIntMax = max; }
	void SetPid(Bool_t flag=1) { fPid = flag; }
	void SetSecondaries(Bool_t flag=1) { fSecondaries = flag; }
	void SetSecProcedure(Int_t proc) { fSecProc = proc; }
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
	
	void SetPidProcedure(Int_t proc) { fPidProc=proc; }
	void SetPidEfficiency(Double_t eff) { fPidEff=eff; }
	void SetAddFakeTracks(Bool_t flag=1) { fAddFakeTracks = flag; }
	
	void SetDebugLevel(Int_t level) { fDebugLevel = level; }
	
  private:
 
	AliLnDriver(const AliLnDriver& other);
	AliLnDriver& operator=(const AliLnDriver& other);
	
  private:
 
	TString fSpecies;       // particle species
	TString fOutputTag;     // tag for output file
	TString fOutputCorTag;  // tag for correction file
	
	Double_t fTrigEff[3];   // trigger efficiency, stat. and syst. errors
	Bool_t  fIsOnlyGen;     // if it is only generation
	Bool_t  fINEL;          // extrapolate to inelastic events
	
	Bool_t  fMakeCorr;      // make corrections
	Bool_t  fMakePt;        // make pt
	Bool_t  fMakeRatio;     // make antiparticle/particle ratio
	Bool_t  fMakeSpectra;   // make spectra
	Bool_t  fMakeStats;     // make event stats
	
	Double_t fPtMin;        // minimum pt value
	Double_t fPtMax;        // maximum pt value
	Bool_t   fPid;          // enable pid correction
	Double_t fPidPt;        // minimum pt value for pid correction
	Bool_t   fSecondaries;  // enable correction of secondaries
	Int_t    fSecProc;      // procedure to estimate fractions
	Int_t    fMatDCAxyMod;  // DCAxy model for correction of secondaries
	Bool_t   fANucTemplate; // enable antinucleus as template for primaries
	Int_t    fNbin;         // rebin of DCAxy distribution
	Double_t fYMin;         // min rapidity
	Double_t fYMax;         // max rapidity
	Double_t fMinDCAxy;     // min DCAxy
	Double_t fMaxDCAxy;     // max DCAxy
	Double_t fBkgMin;       // lower limit for removing background
	Double_t fBkgMax;       // upper limit for removing background
	Double_t fIntMin;       // lower limit for integration
	Double_t fIntMax;       // upper limit for integration
	Bool_t   fEfficiency;   // enable efficiency correction
	Bool_t   fG3Fluka;      // enable G3/Fluka correction for TPC
	Double_t fScMat;        // scaling factor for material fraction
	Double_t fScFd;         // scaling factor for feed-down fraction
	
	TString fInputData;     // input data filename
	TString fInputSimu;     // input simulation filename
	TString fInputSimuFix;  // input fixed simulation filename
	
	TString fOutputPtCorr;  // output filename for pt corrections
	TString fOutputPt;      // output filename for pt
	TString fOutputRatio;   // output filename for antiparticle/particle ratio
	TString fOutputSpectra; // output filename for differential yields
	
	TString fOutputPtCorrDebug; // output filename for debugging pt corrections
	TString fOutputPtDebug;     // output filename for debugging pt
	
	Bool_t fFitFrac;        // enable fit to fractions
	Bool_t fFdwnCorr;       // enable feed-down correction
	Bool_t fSameFdwn;       // same feed-down correction for positives and negatives
	Bool_t fMCtoINEL;       // MC to extrapolate to inel or for triggering events
	
	Bool_t fAddFakeTracks;  // include fake tracks in the efficiency and templates
	
	Int_t fPidProc;  // pid procedure for the pt distribution
	
	Double_t fPidEff;  // pid efficiency for all pt bins
	
	Int_t fDebugLevel; // 0 no verbose, > 1 verbose
	
	ClassDef(AliLnDriver,3)
};

#endif // ALILNDRIVER_H
