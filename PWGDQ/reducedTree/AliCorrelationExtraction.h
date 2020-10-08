/**************************************************************
 *                                                            *
 * class used for the extraction of J/psi-hadron correlations *
 *                                                            *
 * authors: Lucas Altenkamper (lucas.altenkamper@cern.ch)     *
 *          Antoine Lardeux   (antoine.lardeux@cern.ch)       *
 *                                                            *
 * 12/08/2018                                                 *
 *                                                            *
 **************************************************************/

/*

 Brief usage guide: ---------------------------------------------------------------------------------------------

 This class is supposed to be used as a toolbox for the J/psi - hadron correlation analysis. It provides
 all necessary functionality to extract the correlation signal based on e+e- - hadron correlation
 distributions (provided as THnF or THnSparseF) obtained with the reduced tree framework. The class uses
 the AliResonanceFits class to extract the J/psi signal from the e+e- invariant mass distributions. The
 AliResonanceFits object is supposed to be configured in a macro and provided to the AliCorrelationExtraction
 class. A brief overview of the main functionality of the class is provided in the following. For a more
 detailed description of the different functions check AliCorrelationExtraction.cxx.

 1) Input data --------------------------------------------------------------------------------------------------

  The main input distributions are the same-event opposite sign (fSEOS) and mixed-event opposite sign (fMEOS)
  THnF (or THnSparseF) distributions obtained from the reducedTree framework, namely the class
  AliReducedAnalysisJpsi2eeCorrelations. Hereby, opposite sign refers to the e+e- pairs used in the
  e+e- - hadron correlation histograms. The input distributions are set using SetSEOSHistogram() and
  SetMEOSHistogram(), respectively.

  Additionally, depending on the correlation signal extraction method that is used, like sign distributions
  are required. Hereby, correlation distributions (fSEMM, fSEPP, fMEMM and fMEPP) and electron/positron pair
  distributions (fSEPPPair, fSEMMPair, fMEPPPair and fMEMMPair) are set with the functions SetSELSHistogram(),
  SetMELSHistogram(), SetSELSPairHistogram() and SetMELSPairHistogram(), respectively.

  All input THnF (or THnSparseF) can also be set at the same time using SetHistograms() for the correlation
  and SetPairHistograms() for the pair distributions.

  In order to correct the J/psi - hadron correlation function for the associated hadron efficiency, an
  1-dimensional efficiency distribution (fHadronEff) can be set via SetHadronEfficiencyHistogram(). Note that
  the variable/dimension (see 2) in which the efficiency is provided has to be set.

 2) Dimensions and variables ------------------------------------------------------------------------------------

  There are 3 "types" of variables (THnF/THnSparseF dimensions) that can be (2.1 and 2.3) or need to be (2.2)
  defined by the user in order to do the correlation signal extraction with this class:

  2.1) General variables, which correspond to the input THnF/THnSparseF dimensions (i.e. axes), can be set
    using AddVariable() or AddVariables(). Hereby, variable ranges can be set using SetVarRange() in order to
    cut on the axes of the input distributions.
    Not all dimensions have to be provided and those which are not provided by the user or are not given a
    specific variable range will be integrated over.

  2.2) Correlation variables, i.e. the delta phi, delta eta and e+e- pair mass variables, are required for the
    correlation signal extraction and can be set using SetDeltaPhiVariable(), SetDeltaEtaVariable() and
    SetMassVariable(), respectively. In addition to the specific setting using the 3 functions quoted above,
    the corresponding dimension of the input THnF/THnSparseF has to be defined using AddVariable() or
    AddVariables() as described in 2.1. Variable ranges, i.e. cuts, can also be defined for the 3 correlation
    variables using SetVarRange().
    In addition, if the like sign distributions (see 1) are used for the correlation signal extraction, the
    index (i.e. the axis) of the mass variable in the like sign electron pair input distributions has to be
    defined using SetMassVariablePairIndex().

  2.3) Mixing variables, i.e. event variables which are used to define the event mixing pool, can be provided
    using AddMixingVariable() or AddMixingVariables(). These have to be dimensions of the input THnF/THnSparseF
    and, if provided, will be used during the calculation of the inclusive e+e- - hadron correlation extraction.
    Hereby, the inclusive correlation function will be extracted in bins of the mixing variables and summed
    afterwards.
    Typical examples for these variables are the z_vtx. position or the event multiplicity. These variables are
    not required for the correlation signal extraction and a test on 2016 pp 13 TeV data showed no significant
    difference between using the mixing variables (z_vtx, mult, and both) or not.

 3) User options ------------------------------------------------------------------------------------------------

  There are different options for the correlation signal extraction which can be set by the user.

  3.1) The background method, all of which can be found in the enumerated list BackgroundMethods, that is used
    to determine the non-J/psi background in the e+e- - hadron correlation function can be set via
    SetBackgroundMethod(). By default, no background extraction method is set and the class will not run the
    processing.

  3.2) If a fitting method (kBkgFitting) is requested for the background determination, a fit function must be
    provided via SetBackgroundFitFunction(). There is no default fit function defined. A mass exclusion range
    that should not be considered for the background fit is required as well as there is no default. It must be
    set via SetMassExclusionRange().

  3.3) The e+e- pair mass signal window must be defined using SetSignalMassWindow(). Again, there is not default
    setting and the processing will be stopped if none is provided.

  3.4) The e+e- pair background mass window(s) can be set with SetBackgroundMassWindows() and define the mass
    windows used for the correlation background extraction with the various background extraction methods. Note
    that some methods (e.g. kBkgSideband) require only one background window while others (e.g. kBkgFitting or
    kBkgSuperposition) require more. Check Initialize() or the different background extraction functions in order
    to find out how many windows you need.

  3.5) Additionally, a verbose mode can be used by calling SetVerbose(). Per default, verbose mode is switched off.

 4) Processing --------------------------------------------------------------------------------------------------

  The processing, i.e. the correlation signal extraction, is done by calling Process() in the user defined macro.
  This function returns a flag (fProcessDone) that can be tested to see if the processing was successful. The
  processing works as follows:

  4.1) Initialize() is called first and will stop the processing in case of erroneous inputs or conflicting user
    options.

  4.2) The J/psi signal is extracted using the user provided AliResonanceFits object.

  4.3) The inclusive e+e- - hadron correlation function is calculated in the user provided mass signal window using
    the (private) function CalculateInclusiveCorrelation().

  4.4) Depending on the user options for the background extraction method, a specific function to determined the
    background correlation function (CalculateBackgroundCorrelation...) is called. For some background extraction
    methods (kBkgSuperposition and kBkgFitting), the signal correlation function (i.e. J/psi - hadron) is already
    calculated here and 4.5 is skipped.

  4.5) The signal correlation function is calculated using the superposition principle from the inclusive
    e+e- - hadron and the background correlation function using CalculateSignalCorrelation().

  4.6) There are two ways an efficiency correction can be applied, depending on if the J/psi should be efficiency
    corrected or not: 1) If no J/psi efficicency correction is required, a hadron efficiency map (1D) can be set
    using SetHadronEfficiencyHistogram(). The signal correlation will then be corrected for the hadron
    efficiency after the signal extraction has been performed. 2) In case a J/psi efficiency correction is required
    in addition to the hadron efficiency correction, the input correlation histograms already have to be corrected
    for the efficiency using a 1/(eff_jpsi x eff_hadron) weight when filling the input histograms. If this is the case,
    the J/psi efficiency must be provided via SetJpsiEfficiency() in order to properly correct the trigger normalization.

 5) Output ------------------------------------------------------------------------------------------------------

  After the processing is done, the final and some intermediate distributions can be obtained using the
  provided getter functions. The most important output objects hereby are the inclusive (all e+e-), background
  (non-J/psi e+e-) and signal (J/psi) correlation functions which can be obtained using the GetInclusiveCF...(),
  GetBackgroundCF...() and GetSignalCF...() functions, respectively. For all other possible output objects,
  please consult the class in order to understand in detail what they represent.

  Additionally, the trigger (i.e. e+e- pair) values caluclated using the AliResonanceFits object can be obtained
  in the signal and background windows using GetTriggerValuesSignal() and GetTriggerValuesBackground(),
  respectively. The trigger values are provided as arrays with the elements according to the enumerated list
  TriggerValues, including e.g. the signal (J/psi) and background counts.

*/

#ifndef ALICORRELATIONEXTRACTION_H
#define ALICORRELATIONEXTRACTION_H

#include <TObject.h>
#include <THn.h>
#include <THnSparse.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include "AliResonanceFits.h"

class TMinuit;
class TH1;

//_______________________________________________________________________________
class AliCorrelationExtraction : public TObject {
  
  public:

    AliCorrelationExtraction();
    virtual ~AliCorrelationExtraction();
  
    enum BackgroundMethods {
      kBkgNone=-1,
      kBkgFitting,
      kBkgSideband,
      kBkgLikeSign,
      kBkgInterpolation,
      kBkgSuperposition,
      kNBackgroundMethods=5
    };
    enum TriggerValues {
      kSig=0,
      kSigErr,
      kBkg,
      kBkgErr,
      kSplusB,
      kSplusBErr,
      kNTrigSEPP,         // no. trigger SE-PP (N_++)
      kNTrigSEPPErr,
      kNTrigSEMM,         // no. trigger SE-MM (N_--)
      kNTrigSEMMErr,
      kNTrigMEOS,         // no. trigger ME-OS (M_+-)
      kNTrigMEOSErr,
      kNTrigMEPP,         // no. trigger ME-PP (M_++)
      kNTrigMEPPErr,
      kNTrigMEMM,         // no. trigger ME-MM (M_--)
      kNTrigMEMMErr,
      kR,                 // M_+-/2*sqrt(M_++*M_--)
      kRErr,
      kBkgComb,           // 2*R*sqrt(N_++*N_--)
      kBkgCombErr,
      kSigFrac,           // S/S+B
      kSigFracErr,
      kBkgFrac,           // B/S+B
      kBkgFracErr,
      kBkgCombFrac,       // B_comb/S+B
      kBkgCombFracErr,
      kNTriggerValues=26
    };

    // maximum array dimensions
    static const Int_t kNMaxVariables             = 10;   // variables (i.e. THnF dimensions)
    static const Int_t kNMaxMixingVariables       = 10;   // mixing variables (i.e. THnF dimensions)
    static const Int_t kNMaxBackgroundMassRanges  = 10;   // background mass windows
    static const Int_t kNMaxMixingVarBins         = 100;  // mixing variable bins in THnF
    static const Int_t kNMaxDeltaPhiBins          = 100;  // delta phi bins in THnF
    static const Int_t kNMaxDeltaEtaBins          = 100;  // delta eta bins in THnF
  
    // setters
    void SetHistograms(THnF* seos, THnF* sepp=0x0, THnF* semm=0x0, THnF* meos=0x0, THnF* mepp=0x0, THnF* memm=0x0);
    void SetHistograms(THnSparseF* seos, THnSparseF* sepp=0x0, THnSparseF* semm=0x0, THnSparseF* meos=0x0, THnSparseF* mepp=0x0, THnSparseF* memm=0x0);
    void SetPairHistograms(THnF* sepp=0x0, THnF* semm=0x0, THnF* meos=0x0, THnF* mepp=0x0, THnF* memm=0x0);
    void SetSEOSHistogram(THnF* h) {fSEOS = h; fProcessDone = kFALSE;}
    void SetSEOSHistogram(THnSparseF* h) {fSEOSSparse = h; fProcessDone = kFALSE;}
    void SetMEOSHistogram(THnF* h) {fMEOS = h; fProcessDone = kFALSE;}
    void SetMEOSHistogram(THnSparseF* h) {fMEOSSparse = h; fProcessDone = kFALSE;}
    void SetSELSHistogram(THnF* hpp, THnF* hmm) {fSEPP = hpp; fSEMM = hmm; fProcessDone = kFALSE;}
    void SetSELSHistogram(THnSparseF* hpp, THnSparseF* hmm) {fSEPPSparse = hpp; fSEMMSparse = hmm; fProcessDone = kFALSE;}
    void SetMELSHistogram(THnF* hpp, THnF* hmm) {fMEPP = hpp; fMEMM = hmm; fProcessDone = kFALSE;}
    void SetMELSHistogram(THnSparseF* hpp, THnSparseF* hmm) {fMEPPSparse = hpp; fMEMMSparse = hmm; fProcessDone = kFALSE;}
    void SetMEOSPairHistogram(THnF* h) {fMEOSPair = h; fProcessDone = kFALSE;}
    void SetSELSPairHistogram(THnF* hpp, THnF* hmm) {fSEPPPair = hpp; fSEMMPair = hmm; fProcessDone = kFALSE;}
    void SetMELSPairHistogram(THnF* hpp, THnF* hmm) {fMEPPPair = hpp; fMEMMPair = hmm; fProcessDone = kFALSE;}
    void SetJpsiEfficiency(Double_t eff, Double_t effErr) {fJpsiEff = eff; fJpsiEffErr = effErr; fUseJpsiEfficiency = kTRUE; fProcessDone = kFALSE;}
    void SetHadronEfficiencyHistogram(TH1D* h, Int_t var) {fHadronEff = h; fHadronEfficiencyVariable = var; fProcessDone = kFALSE;}
    void SetAliResonanceFitsObject(AliResonanceFits* resonanceFits) {fResonanceFits = (AliResonanceFits*)resonanceFits->Clone("ResonanceFits"); fProcessDone = kFALSE;}
    void SetBackgroundMethod(Int_t method, Bool_t integrateDeltaEta=kFALSE);
    void SetBackgroundFitFunction(TF1* fitFunc) {fBkgFitFunction = (TF1*)fitFunc->Clone("BkgFitFunction"); fProcessDone = kFALSE;}
    void SetFitPrecision(Double_t prec=1.e-06) {fFitPrecision = prec;}
    void SetSignalMassWindow(Double_t min, Double_t max) {fMassSignalRange[0] = min; fMassSignalRange[1] = max; fProcessDone = kFALSE;}
    void SetBackgroundMassWindows(Int_t n, Double_t* min, Double_t* max);
    void SetMassExclusionRange(Double_t min, Double_t max) {fMassExclusionRange[0] = min, fMassExclusionRange[1] = max;}
    void SetVerbose() {fVerboseFlag = kTRUE;}
  
    // add variables and set ranges on the THnF
    void AddVariables(Int_t nVars, Int_t* vars, Int_t* indices);
    void AddVariable(Int_t var, Int_t index);
    void SetVarRange(Int_t var, Double_t* lims);
    void SetVarRange(Int_t var, Double_t min, Double_t max);
    void AddMixingVariables(Int_t nVars, Int_t* vars, Int_t* indices);
    void AddMixingVariable(Int_t var, Int_t index);

    // indicate mass, delta phi and delta eta variables
    void SetMassVariable(Int_t var) {fMassVariable = var; fProcessDone = kFALSE;}
    void SetDeltaPhiVariable(Int_t var) {fDeltaPhiVariable = var; fProcessDone = kFALSE;}
    void SetDeltaEtaVariable(Int_t var) {fDeltaEtaVariable = var; fProcessDone = kFALSE;}
    void SetMassVariablePairIndex(Int_t index) {fMassVariableIndexPair = index; fProcessDone = kFALSE;}

    // getters
    TH1D*             GetInclusiveCF1D() const {return (fProcessDone ? fInclusiveCF1D : 0x0);}
    TH1D*             GetInclusiveCF1D(Int_t massWindow) const {return (fProcessDone ? fInclusiveCF1DBackgroundMassWindow[massWindow] : 0x0);}
    TH1D*             GetInclusiveCF1D(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMass[phiBin][etaBin] : 0x0);}
    TH1D*             GetInclusiveCF1DBackgroundRange(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMassBackgroundRange[phiBin][etaBin] : 0x0);}
    TH2D*             GetInclusiveCF2D() const {return (fProcessDone ? fInclusiveCF2D : 0x0);}
    TH2D*             GetInclusiveCF2D(Int_t massWindow) const {return (fProcessDone ? fInclusiveCF2DBackgroundMassWindow[massWindow] : 0x0);}
    TH3D*             GetInclusiveCF3D() const {return (fProcessDone ? fInclusiveCF3D : 0x0);}
    TH1D*             GetBackgroundCF1D() const {return (fProcessDone ? fBackgroundCF1D : 0x0);}
    TH2D*             GetBackgroundCF2D() const {return (fProcessDone ? fBackgroundCF2D : 0x0);}
    TH1D*             GetCombinatorialBackgroundCF1D(Int_t massWindow) const {return (fProcessDone ? fCombinatorialBackgroundCF1D[massWindow] : 0x0);}
    TH2D*             GetCombinatorialBackgroundCF2D(Int_t massWindow) const {return (fProcessDone ? fCombinatorialBackgroundCF2D[massWindow] : 0x0);}
    TH1D*             GetSignalCF1D() const {return (fProcessDone ? fSignalCF1D : 0x0);}
    TH1D*             GetSignalCF1DEfficiencyCorrected() const {return (fProcessDone ? fSignalCF1DEffCorr : 0x0);}
    TH2D*             GetSignalCF2D() const {return (fProcessDone ? fSignalCF2D : 0x0);}
    TH2D*             GetSignalCF2DEfficiencyCorrected() const {return (fProcessDone ? fSignalCF2DEffCorr : 0x0);}
    AliResonanceFits* GetAliResonanceFitsObject() const {return (fProcessDone ? fResonanceFits : 0x0);}
    Int_t             GetBackgroundMethod() const {return (fProcessDone ? fOptionBkgMethod : 0x0);}
    TF1*              GetBackgroundFitFunction() const {return (fProcessDone ? fBkgFitFunction : 0x0);}
    TF1*              GetBackgroundFitFunction(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fBackgroundCF1DInvMassFit[phiBin][etaBin] : 0x0);}
    TF1*              GetGlobalFitFunction(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMassFit[phiBin][etaBin] : 0x0);}
    TGraphErrors*     GetGlobalFitFunctionCI(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fInclusiveCF1DInvMassFitCI[phiBin][etaBin] : 0x0);}
    TH1*              GetSoverBMC(Int_t phiBin, Int_t etaBin) const {return (fProcessDone ? fSoverBMC[phiBin][etaBin] : 0x0);}
    TH1D*             GetPairInvMassSEPP() const {return (fProcessDone ? fSEPPPairInvMass : 0x0);}
    TH1D*             GetPairInvMassSEMM() const {return (fProcessDone ? fSEMMPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEOS() const {return (fProcessDone ? fMEOSPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEPP() const {return (fProcessDone ? fMEPPPairInvMass : 0x0);}
    TH1D*             GetPairInvMassMEMM() const {return (fProcessDone ? fMEMMPairInvMass : 0x0);}
    TH2D*             GetSEOS() const {return (fProcessDone ? fSEOSNorm : 0x0);}
    TH2D*             GetSEOS(Int_t massWindow) const {return (fProcessDone ? fSEOSNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetSEPP(Int_t massWindow) const {return (fProcessDone ? fSEPPNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetSEMM(Int_t massWindow) const {return (fProcessDone ? fSEMMNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEOS() const {return (fProcessDone ? fMEOSNorm : 0x0);}
    TH2D*             GetMEOS(Int_t massWindow) const {return (fProcessDone ? fMEOSNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEPP(Int_t massWindow) const {return (fProcessDone ? fMEPPNormBackgroundMassWindow[massWindow] : 0x0);}
    TH2D*             GetMEMM(Int_t massWindow) const {return (fProcessDone ? fMEMMNormBackgroundMassWindow[massWindow] : 0x0);}
    const Double_t*   GetTriggerValuesSignal() const {return (fProcessDone ? fTrigValSig : 0x0);}
    const Double_t*   GetTriggerValuesBackground(Int_t massWindow) const {return (fProcessDone ? fTrigValBkg[massWindow] : 0x0);}

    // processing
    Bool_t  Process();
    void    PrintUserOptions();  // print summary of user options
  
  private:
  
    // input histograms
    THnF*       fSEOS;
    THnSparseF* fSEOSSparse;
    THnF*       fSEPP;
    THnSparseF* fSEPPSparse;
    THnF*       fSEMM;
    THnSparseF* fSEMMSparse;
    THnF*       fMEOS;
    THnSparseF* fMEOSSparse;
    THnF*       fMEPP;
    THnSparseF* fMEPPSparse;
    THnF*       fMEMM;
    THnSparseF* fMEMMSparse;
    THnF*       fSEPPPair;
    THnF*       fSEMMPair;
    THnF*       fMEOSPair;
    THnF*       fMEPPPair;
    THnF*       fMEMMPair;
  
    TH1D*       fHadronEff;

    static TH1* fSoverBMC[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];

    // output histograms
    TH2D* fSEOSNorm;
    TH2D* fSEOSNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fSEPPNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fSEMMNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH1D* fSEPPPairInvMass;
    TH1D* fSEMMPairInvMass;
    TH1D* fMEOSPairInvMass;
    TH1D* fMEPPPairInvMass;
    TH1D* fMEMMPairInvMass;
    TH2D* fMEOSNorm;
    TH2D* fMEOSNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fMEPPNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fMEMMNormBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH1D* fInclusiveCF1D;
    TH1D* fInclusiveCF1DInvMass[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH1D* fInclusiveCF1DInvMassBackgroundRange[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH1D* fInclusiveCF1DBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH2D* fInclusiveCF2D;
    TH2D* fInclusiveCF2DBackgroundMassWindow[kNMaxBackgroundMassRanges];
    TH3D* fInclusiveCF3D;
    TH1D* fBackgroundCF1D;
    static TF1*   fBackgroundCF1DInvMassFit[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TF1*          fInclusiveCF1DInvMassFit[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TGraphErrors* fInclusiveCF1DInvMassFitCI[kNMaxDeltaPhiBins][kNMaxDeltaEtaBins];
    TH2D* fBackgroundCF2D;
    TH1D* fCombinatorialBackgroundCF1D[kNMaxBackgroundMassRanges];
    TH2D* fCombinatorialBackgroundCF2D[kNMaxBackgroundMassRanges];
    TH1D* fSignalCF1D;
    TH1D* fSignalCF1DEffCorr;
    TH2D* fSignalCF2D;
    TH2D* fSignalCF2DEffCorr;

    // variables
    Int_t     fNVariables;                                    // number of variables to be handled
    Int_t     fVariables[kNMaxVariables];                     // list of variables
    Double_t  fVarLimits[kNMaxVariables][2];                  // variable limits (i.e. "cuts" or selection ranges)
    Int_t     fVarIndices[kNMaxVariables];                    // indices of variables in THnF
    Int_t     fNMixingVariables;                              // number of variables used in event mixing
    Int_t     fMixingVariables[kNMaxMixingVariables];         // list of mixing variables
    Int_t     fMixingVarIndices[kNMaxMixingVariables];        // indices of mixing variables in THnF
    Int_t     fNMixingVarBins[kNMaxMixingVariables];          // number of bins per mixing variable
    Double_t  fMixingVarBinLimits[kNMaxMixingVariables][kNMaxMixingVarBins][2]; // mixing variable bin limits
    Int_t     fMassVariable;
    Int_t     fMassVariableIndex;
    Int_t     fMassVariableIndexPair;                         // mass variable index in LS pair histogram: fSEPPPair and fSEMMPair
    Int_t     fDeltaPhiVariable;
    Int_t     fDeltaPhiVariableIndex;
    Int_t     fDeltaEtaVariable;
    Int_t     fDeltaEtaVariableIndex;
    Int_t     fHadronEfficiencyVariable;                      // variable the efficiency (fHadronEff) is a function of
    Int_t     fHadronEfficiencyVariableIndex;                 // index of fHadronEfficiencyVariable in fSEOS, etc.

    // J/psi efficiency
    Double_t  fJpsiEff;
    Double_t  fJpsiEffErr;

    // user options
    Bool_t            fVerboseFlag;
    Bool_t            fUseMixingVars;
    Bool_t            fIntegrateDeltaEta[kNBackgroundMethods];
    Bool_t            fUseJpsiEfficiency;
    AliResonanceFits* fResonanceFits;
    Int_t             fOptionBkgMethod;
    TF1*              fBkgFitFunction;
    Double_t          fFitPrecision;
  
    // mass ranges
    Int_t     fNBackgroundMassRanges;                               // number of background mass windows
    Double_t  fBackgroundMassRanges[kNMaxBackgroundMassRanges][2];  // background mass windows
    Double_t  fMassSignalRange[2];                                  // signal mass window
    Double_t  fMassExclusionRange[2];                               // mass exclusion range for fit method
  
    // values
    Double_t fTrigValSig[kNTriggerValues];                            // trigger values signal mass range
    Double_t fTrigValBkg[kNMaxBackgroundMassRanges][kNTriggerValues]; // trigger values background mass ranges
  
    // process flag
    Bool_t fProcessDone;
  
    // member functions
    void    ApplyUserRanges(THnBase* h);
    Bool_t  Initialize();
    Bool_t  NormalizeToNearSidePeak(TH2D* h);
    Bool_t  IsBackgroundRange(Double_t min, Double_t max, Int_t& index);
    TH1D*   ProjectToDeltaPhi(TH2D* hIn, TString name);
    static Double_t GlobalFitFunction(Double_t* x, Double_t* par);
    Bool_t  CalculateInclusiveCorrelationInMixingBins(Int_t currentVar, Int_t& nCalls,
                                                      THnBase* seos, THnBase* meos, TH2D* (&inclCF));
    Bool_t  CalculateInclusiveCorrelation(Double_t minMass, Double_t maxMass, Bool_t isSignalRange,
                                          TH2D* (&seos), TH2D* (&meos),
                                          TH2D* (&incl2D), TH1D* (&incl1D));
    Bool_t  CalculateBackgroundCorrelationFitting();
    Bool_t  CalculateBackgroundCorrelationSideband();
    Bool_t  CalculateBackgroundCorrelationLikeSign();
    Bool_t  CalculateBackgroundCorrelationInterpolation();
    Bool_t  CalculateBackgroundCorrelationSuperposition();
    Bool_t  CalculateSignalCorrelation();
    Bool_t  HadronEfficiencyCorrection();
  
  ClassDef(AliCorrelationExtraction, 7);
};

#endif
