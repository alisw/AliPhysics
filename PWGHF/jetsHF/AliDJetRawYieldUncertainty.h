/**
 * \file AliDJetRawYieldUncertainty.h
 * \brief Declaration of class AliDJetRawYieldUncertainty
 *
 * In this header file the class AliDJetRawYieldUncertainty is declared.
 * Class to extract jet Pt spectrum yield uncertainty via multi-trial approach.
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Feb 28, 2016
 */

#ifndef ALIDJETRAWYIELDUNCERTAINTY_H
#define ALIDJETRAWYIELDUNCERTAINTY_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <vector>

class TH1D;
class TH2D;
class TH1F;
class TCanvas;
class AliHFMultiTrials;

class AliDJetVReader;

/**
 * \class AliDJetRawYieldUncertainty
 * \brief Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach.
 *
 * Implementation of a class to extract jet Pt spectrum yield uncertainty via multi-trial approach.
 */
class AliDJetRawYieldUncertainty : public TObject {

public:

  enum EDMesonSpecies_t {
    kUnknownMeson,
    kD0toKpi,
    kDStarD0pi
  };

  enum EYieldMethod_t {
    kUnknownMethod,
    kEffScale,
    kSideband
  };

  /**
   * \struct SBResults
   * \brief Implementation of a struct used to hold results of the SB uncertainty evalualtion
   *
   * Implementation of a struct used to hold results of the SB uncertainty evalualtion
   */
  struct SBResults {
    Bool_t     fSuccess                  ; //!<!Status of the last run
    TH1D      *fJetYieldCentral          ; //!<!Central values of the yield of jet spectrum + syst yield uncertainty
    TH1D      *fJetYieldUnc              ; //!<!Yield uncertainty vs jet pT bin
    TH1F     **fJetSpectrSBVars          ; //!<!Array of jet spectrum histograms, one per variation (sideband approach)
    TH1F      *fJetSpectrSBDef           ; //!<!Array of jet spectrum histograms, default trial (sideband approach)
    TH1F     **fJetBinYieldDistribution  ; //!<!Array of histograms with yield distributions from the trials for each pT(jet)
  };

  AliDJetRawYieldUncertainty();
  AliDJetRawYieldUncertainty(const AliDJetRawYieldUncertainty &source);
  virtual ~AliDJetRawYieldUncertainty();

  Bool_t SetDmesonSpecie(EDMesonSpecies_t k);
  Bool_t SetYieldMethod(EYieldMethod_t meth);
  void SetAllowRepetitionOfTrialExtraction(Bool_t allow)           { fAllowRepetitions       = allow ; }
  void SetPtBinEdgesForMassPlot(Double_t ptmin, Double_t ptmax)    { fpTmin                  = ptmin ; fpTmax                  = ptmax ; }
  void SetZedges(Double_t zmin, Double_t zmax)                     { fzmin                   = zmin  ; fzmax                   = zmax  ; }
  void SetSigmaForSignalRegion(Double_t nsigma)                    { fnSigmaSignReg          = nsigma; }
  void SetSigmaSideBandLeft(Double_t nsig1, Double_t nsig2)        { fnSigmaSideBandLeft1    = nsig1 ; fnSigmaSideBandLeft2    = nsig2 ; }
  void SetSigmaSideBandRight(Double_t nsig1, Double_t nsig2)       { fnSigmaSideBandRight1   = nsig1 ; fnSigmaSideBandRight2   = nsig2 ; }
  void SetMaxNTrialsForSidebandMethod(Int_t nmax)                  { fnMaxTrials             = nmax  ; }
  void SetSaveInvMassFitCanvases(Bool_t s)                         { fSaveInvMassFitCanvases = s     ; }
  void SetChi2Cut(Double_t cut)                                    { fChi2Cut                = cut   ; }
  void SetDebugLevel(Int_t d)                                      { fDebug                  = d     ; }
  void SetDJetReader(AliDJetVReader* r)                            { fDJetReader             = r     ; }

  // Reflections
  void SetFitReflections(Bool_t refl)                              { fFitRefl                = refl  ; }
  void SetReflFilename(TString fname)                              { fReflFilenameInput      = fname ; }
  void SetMCSigFilename(TString fname)                             { fSigMCFilenameInput     = fname ; }
  void SetReflHistoname(TString hname)                             { fReflHistoName          = hname ; }
  void SetMCSigHistoname(TString hname)                            { fSigMCHistoName         = hname ; }
  void SetUseCoherentChoice(Bool_t c)                              { fCoherentChoice         = c     ; }
  void SetUseBkgInBinEdges(Bool_t b)                               { fUseBkgInBinEdges       = b     ; }

  void SetValueOfReflOverSignal(Double_t ratio, Double_t minrange=1.7, Double_t maxrange=2.1) { fFixRiflOverS = ratio; fReflRangeL = minrange; fReflRangeR = maxrange; }

  void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetJetzBins(Int_t nbins=0, Double_t* zedges=0x0);
  void SetDmesonEfficiency(Double_t* effvalues=0x0);

  void SetSigmaToFixDPtBins(Double_t* sigmafix);
  void SetSigmaToFixJetPtBins(Double_t* sigmafix);

  void SetMeanSigmaVariations(Bool_t* cases);
  void SetBkgVariations(Bool_t* cases);
  void SetRebinSteps(Int_t nsteps, Int_t* cases);
  void SetMinMassSteps(Int_t nsteps, Double_t* cases);
  void SetMaxMassSteps(Int_t nsteps, Double_t* cases);
  void SetSigmaBinCounting(Int_t nsteps, Double_t* cases=0x0);
  void SetMaskOfVariations(Int_t ncases, Bool_t* cases);

  static void FitReflDistr(Int_t nPtBins, TString inputfile, TString fitType = "DoubleGaus");

  AliHFMultiTrials* RunMultiTrial();
  Bool_t CombineMultiTrialOutcomes();

  Bool_t ExtractInputMassPlot();

  Bool_t EvaluateUncertainty();
  Bool_t EvaluateUncertaintyEffScale();
  SBResults EvaluateUncertaintySideband(TString obs, Int_t nJetBins, Double_t* jetBinEdges);
  Bool_t EvaluateUncertaintySideband();

  Bool_t Success() const { return fSuccess; }

  Bool_t GenerateJetSpectrum(TH2* hInvMassJetObs, Double_t mean, Double_t sigma, Double_t bkg, Int_t iDbin, TH1* hjetobs, TH1* hjetobs_s, TH1* hjetobs_s1, TH1* hjetobs_s2);

  virtual void ClearObjects();

protected:

  static const Int_t fgkNSigmaVar             = 6; //!<!Number of mean/sigma variations
  static const Int_t fgkNBkgVar               = 8; //!<!Number of bkg variations

  AliDJetVReader    *fDJetReader                 ; ///< Pointer to the object that reads the output of the analysis task and produces the histograms
  Bool_t             fSaveInvMassFitCanvases     ; ///< Switch to save the invariant mass fit canvases
  EDMesonSpecies_t   fDmesonSpecie               ; ///< D meson specie
  TString            fDmesonLabel                ; ///< D meson label
  EYieldMethod_t     fYieldApproach              ; ///< Method to extract jet pT spectrum
  TString            fMethodLabel                ; ///< Method label
  Double_t           fpTmin                      ; ///< pT lower edge of mass plot to evaluate variations of yields
  Double_t           fpTmax                      ; ///< pT upper edge of mass plot to evaluate variations of yields
  Double_t           fzmin                       ; ///< z minimum value to extract jet pT spectrum
  Double_t           fzmax                       ; ///< z maximum value to extract jet pT spectrum
  Int_t              fnDbins                     ; ///< Number of D-meson pT bins (for eff scaling)
  Double_t          *fDbinpTedges                ; ///< D-meson pt bin edges values
  Int_t              fnJetPtbins                 ; ///< Number of jet pT bins to be used for spectrum
  Double_t          *fJetPtBinEdges              ; ///< Jet pT bin edges to be used for spectrum
  Int_t              fnJetzbins                  ; ///< Number of jet z bins to be used for spectrum
  Double_t          *fJetzBinEdges               ; ///< Jet z bin edges to be used for spectrum
  Double_t          *fDEffValues                 ; ///< D-meson efficiency values

  Double_t           fnSigmaSignReg              ; ///< Number of sigma for signal region
  Double_t           fnSigmaSideBandLeft1        ; ///< Number of sigma for sideband left region (upper limit)
  Double_t           fnSigmaSideBandLeft2        ; ///< Number of sigma for sideband left region (lower limit)
  Double_t           fnSigmaSideBandRight1       ; ///< Number of sigma for sideband right region (upper limit)
  Double_t           fnSigmaSideBandRight2       ; ///< Number of sigma for sideband right region (lower limit)
  Double_t          *fSigmaToFixDPtBins          ; ///< Value of fixed sigma for MultiTrial
  Double_t          *fSigmaToFixJetPtBins        ; ///< Value of fixed sigma for MultiTrial
  Bool_t             fMeanSigmaVar[fgkNSigmaVar] ; ///< Array of bools for mean/sigma variations
  Bool_t             fBkgVar[fgkNBkgVar]         ; ///< Array of bools for bkg variations
  Int_t              fnRebinSteps                ; ///< Number of steps for rebin
  Int_t             *fRebinSteps                 ; ///< Values of rebin steps
  Int_t              fnMinMassSteps              ; ///< Number of steps for low mass edges
  Double_t          *fMinMassSteps               ; ///< Values of low mass edges
  Int_t              fnMaxMassSteps              ; ///< Number of steps for up mass edges
  Double_t          *fMaxMassSteps               ; ///< Values of up mass edges
  Int_t              fnSigmaBC                   ; ///< Number of steps for BC (different sigma)
  Double_t          *fSigmaBC                    ; ///< Values of sigmas for BC
  Int_t              fnMask                      ; ///< Number of elements of array of variations to be kept (only for sigma/mean and bkg config)
  Bool_t            *fMask                       ; ///< Array of variations to be kept (only for sigma/mean and bkg config)
  Double_t           fChi2Cut                    ; ///< Maximum value of allowed chi2
  Int_t              fnMaxTrials                 ; ///< Max number of random trials for each pT(D) bin to build pT(jet) spectrum variations (sideband approach)
  Bool_t             fAllowRepetitions           ; ///< Allow repetitions in the extraction of trials in a give pT(D) bin, for sideband approach

  Bool_t             fFitRefl                    ; ///< Include reflection template in the mass fit
  TString            fReflFilenameInput          ; ///< Name of input file for reflection template
  TString            fSigMCFilenameInput         ; ///< Name of input file for MC signal
  TString            fReflHistoName              ; ///< Name of reflection template histogram
  TString            fSigMCHistoName             ; ///< Name of signal MC histogram
  Double_t           fFixRiflOverS               ; ///< Refl/signMC value for reflections
  Double_t           fReflRangeL                 ; ///< Lower range of reflection template (for refl/signMC)
  Double_t           fReflRangeR                 ; ///< Upper range of reflection template (for refl/signMC)

  Bool_t             fCoherentChoice             ; ///< Use coherent choice for the SB method

  Bool_t             fUseBkgInBinEdges           ; ///< Whether the background should be computed between bin edges rather than exact limits (must be false if the binning of the histogram used for the projections in the SB method is different compared to the histogram used for the invariant mass fit)

  Int_t              fDebug                      ; ///< Debug level

  TH1D              *fMassPlot                   ; //!<!Mass spectra to be fitted
  TH2D              *fMassVsJetPtPlot            ; //!<!Mass vs jet pt (SB method)
  TH2D              *fMassVsJetzPlot             ; //!<!Mass vs jet z (SB method)

  TH1D              *fJetPtYieldCentral          ; //!<!Central values of the yield of jet spectrum + syst yield uncertainty
  TH1D              *fJetPtYieldUnc              ; //!<!Yield uncertainty vs jet pT bin
  TH1F             **fJetPtSpectrSBVars          ; //!<!Array of jet spectrum histograms, one per variation (sideband approach)
  TH1F              *fJetPtSpectrSBDef           ; //!<!Array of jet spectrum histograms, default trial (sideband approach)
  TH1F             **fJetPtBinYieldDistribution  ; //!<!Array of histograms with yield distributions from the trials for each pT(jet)

  TH1D              *fJetzYieldCentral           ; //!<!Central values of the yield of jet spectrum + syst yield uncertainty
  TH1D              *fJetzYieldUnc               ; //!<!Yield uncertainty vs jet pT bin
  TH1F             **fJetzSpectrSBVars           ; //!<!Array of jet spectrum histograms, one per variation (sideband approach)
  TH1F              *fJetzSpectrSBDef            ; //!<!Array of jet spectrum histograms, default trial (sideband approach)
  TH1F             **fJetzBinYieldDistribution   ; //!<!Array of histograms with yield distributions from the trials for each pT(jet)

  Bool_t             fSuccess                    ; //!<!Status of the last run

  std::vector<TCanvas*>
                     fCanvases                   ; //!<!Canvases created by this class

private:
  ClassDef(AliDJetRawYieldUncertainty,1);
};

#endif
