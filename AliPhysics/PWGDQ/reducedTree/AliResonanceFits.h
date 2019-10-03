// Class used for extracting resonance yields from invariant mass distributions
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 
/*
   Brief usage guide
   
   1) Input data
       THnF histograms tipically obtained from reducedTree analyses, but generic THnF can be used as well.
       The histograms can be set via the SetHistograms() methods.
       An initialized AliResonanceFits can be reused multiple times, and these histograms can be reset.
   2) Specifying dimensions / variables
       The AliResonanceFits needs to know what kind of variable corresponds to a given axis in the THnF, but not all axes
       have to be described. See below how the variables are used.
       
       Event - wise variables
       AliResonanceFits is designed to do matching / fitting in each possible event category. 
       The event categories are defined by all possible combinations of (E1, E2, ... Ei, ... En), where Ei is an event - wise variable
       which spans a range [xmin_i -> xmax_i] and is binned in n_i non-equal width bins.
       So, the signal (S), background (B) and signal+background (S+B) distributions are obtained by summing
       the contributions for each event category:
       e.g. B = sum_1 sum_2 ... sum_i ... sum_n  b(e1,e2,...ei,...en) 
       Here b(e1,e2,...,en) is the scaled background mass distribution in a given event category. 
       The sums sum_i is the sum over the bins of category i. The scaling is done in each event category using various methods (described below).
       There is no minimum or maximum of how many event categories should the user define. All the THnF axes that are not defined
       via the AddVariable() are simply integrated over.
       For example, if no event wise variable is defined (even if they are available in the THnF itself) this means that
       the user want to perform the matching / fitting just using the fully integrated S+B and B distributions.
       
       Mass variable:
       It is compulsory to define at least the invariant mass variable. When adding a variable via AddVariable(var, index), the "var"
       argument is an identifier of the variable. If the reducedTree framework is used, then it is convenient to use the variables from 
       AliReducedVarManager. The index argument is the index of the corresponding dimension in the THnF. 
       At the initialization time, there is a check whether the mass variable is defined, by looking in the list of all defined var's
       and comparing them to fMassVariable. If a match is not found, then the the whole procedure fails. By default, fMassVariable is AliReducedVarManager::kMass, so if one is using reducedTree conventions then this is simple. 
       Otherwise, a custom mass variable can be specified via: SetMassVariable().
       
       Pt variable:
       It is not compulsory to have the pt variable in THnF or add it as variable to AliResonanceFits.
       However, if it is added, this has a special treatment to be distinguished from the event variables on which summations are performed.
       The presence of the pt variable is however strongly advisable. If it is defined, the pt variable can be used for
       a) kinematic cuts on the signal
       b) kinematic cuts on the pt range used for matching / fitting (this range can be different wrt to the pt range used to count the signal)
       c) 2D matching: performing (m,pt) matching instead of 1D (matching using only mass distributions)
       The presence of the pt variable among the defined variables is checked at the initialization and by default is considered to be
       AliReducedVarManager::kPt. A custom pt variable can be redefined by SetPtVariable().
       
       Custom ranges for all defined variables can be applied using SetVarRange()
       
   3) User options
       For all the options there are defaults which would give you most of the time something reasonable if all the prerequisites are met.
       The user can specify:
       a) the background type: mixed-event, like-sign, fit (TODO: fit is not currently implemented)
       b) the matching option for scaling the bkg: scale to SE-OS in sidebands around the peak, scale to like-sign
       c) whether to use 1D or 2D matching (NOTE: 2D is a really experimental feature, use with care)
       d) scaling option: using the counts, weighted average of S/B distribution, minuit fit
       e) weighted average power: default is 2 (standard statistics) but this can be changed in certain situations
       f) Minuit fit option if its the case: chi2 minimization or likelihood maximization
       g) Use zero counts bins as significant during chi2 estimation (count is zero with error of 1 count)
       h) Like-sign method when needed: default is using geometric mean; can be changed to use arithmetic mean
       i) Use R-factor correction for the like-sign bkg
       j) after the summation of event categories is finished, an additional matching procedure can be optionally run.
       k) fitting ranges (independent of what is being set via SetVarRange())
         mass fit range
         pt fit range
         mass exclusion range (tipically we exclude the peak region if we match bkg to the SE-OS).
         
   4) Output 
      After setting up all the required inputs, the user should run the Process() function
      If running of Process() is succesful, the S, B, S+B and S/B histograms become available.
      At this point, the function ComputeOutputValues() can be called to calculate integrated signal, S/B, signif, etc.
      The user must provide as arguments, the mass range where the signal should be counted / integrated, and optionally the pt range.
      The return values of ComputeOutputValues() is an array of doubles which contains the information specified in the FitValues enumeration.
      NOTE: If any of the options, input data or ranges of AliResonanceFits are changed, the output values become unavailble 
                 and the Process() function must be called again. So, the output information must be always in sync with the active options.
                 
   5) Plotting
      There is currently no plotting function included. This can be done for the moment by the user in an external function
      using the output of the matching process.
      TODO: Add a plotting solution
       
 */

#ifndef ALIRESONANCEFITS_H
#define ALIRESONANCEFITS_H

#include <TObject.h>
#include <THn.h>
#include <TF1.h>
#include <TFitResult.h>

class TMinuit;
class TH1;

//_____________________________________________________________________
class AliResonanceFits : public TObject {
 
 public:
  
  enum Constants {
     kBkgMixedEvent,    // opposite-sign mixed event bkg
     kBkgLikeSign,         // same-event like-sign bkg
     kBkgMixedEventAndResidualFit,   // fit of residual bkg with a user function for bkg and MC signal shape for signal
     kBkgFitFunction,     // fit of the SE-OS with a user function for bkg and MC signal shape for signal
     kMatchSEOS,          // match to same-event opposite-sign outside signal region (side bands)
     kMatchSELS,           // match to same-event like-sign
     kScaleEntries,         // scale using the bin counts
     kScaleWeightedAverage,  // scale using weighted average of ratios in individual bins
     kScaleFit,                // scale using fitting
     kLSArithmeticMean,    // build LS bkg using arithmetic mean
     kLSGeometricMean,    // build LS bkg using geometric mean
     kMinuitMethodChi2,                // chi2 minimization method
     kMinuitMethodLikelihood,      // log-likelihood maximization
     kDebug,
     kNMaxVariables=10
  };
  
  enum FitValues {
    kSig,                          // signal counts
    kSigErr,                     // signal counts error
    kBkg,                         // bkg counts
    kBkgErr,                    // bkg counts error
    kBkgResidual,           // residual bkg obtained after fitting the combinatorial bkg subtracted minv distribution
    kBkgResidualErr,
    kSplusB,                    // S+B counts
    kSplusBerr,               // S+B counts error
    kSoverB,                   // S/B
    kSoverBerr,               // S/B error
    kSignif,                      // significance
    kSignifErr,                 // significance error
    kChisqSideBands,     // chi2/ndf in the side bands around peak
    kChisqMCPeak,          // chi2/ndf wrt MC signal shape in the peak region
    kChisqMCTotal,     // chi2/ndf wrt MC signal shape
    kFitProbability,        // probability from TFitResult::Prob()
    kMCYieldFraction,     // yield fraction in counting window
    kBkgScale,                     // scale factor for the bkg
    kBkgScaleErr,              // error on the scale factor
    kNFitValues
  };
   
  AliResonanceFits();
  virtual ~AliResonanceFits();            

  
  // setters
  //void Reset();
  // User input
  void SetHistograms(THnF* seos, THnF* meos = 0x0, 
		     THnF* selsLeg1=0x0, THnF* selsLeg2=0x0, THnF* melsLeg1=0x0, THnF* melsLeg2=0x0);
  void SetSEOSHistogram(THnF* hist) {fSEOS = hist; fMatchingIsDone = kFALSE;};
  void SetSELSHistograms(THnF* hLeg1, THnF* hLeg2) {fSELSleg1 = hLeg1; fSELSleg2 = hLeg2; fMatchingIsDone = kFALSE;};
  void SetMEOSHistogram(THnF* hist) {fMEOS = hist; fMatchingIsDone = kFALSE;};
  void SetMELSHistograms(THnF* hLeg1, THnF* hLeg2) {fMELSleg1 = hLeg1; fMELSleg2 = hLeg2; fMatchingIsDone = kFALSE;}
  void SetSEOSMCHistogram(THnF* hist) {fSEOS_MCtruth = hist;}
  void SetSignalMCshape(TH1* shape) {fSignalMCshape = shape;}
  
  // add variables and set ranges on the THnF
  void AddVariables(Int_t nVars, Int_t* vars, Int_t* indices);
  void AddVariable(Int_t var, Int_t index);
  void SetVarRange(Int_t var, Double_t* lims);
  void SetVarRange(Int_t var, Double_t min, Double_t max);
  
  // indicate the special mass and pt variables
  void SetMassVariable(Int_t var) {fMassVariable = var; fMatchingIsDone = kFALSE;}
  void SetPtVariable(Int_t var) {fPtVariable = var; fMatchingIsDone = kFALSE;}
  
  // set various options (see also defaults)
  void SetBkgMethod(Int_t method) {fOptionBkgMethod = method; fMatchingIsDone = kFALSE;}
  void SetMEMatchingMethod(Int_t option) {fgOptionMEMatching = option; fMatchingIsDone = kFALSE;}
  void SetUseRfactorCorrection(Bool_t use=kTRUE) {fOptionUseRfactorCorrection = use; fMatchingIsDone = kFALSE;}
  void SetUse2DMatching(Bool_t flag=kTRUE) {fgOptionUse2DMatching = flag; fMatchingIsDone = kFALSE;}
  void SetScalingOption(Int_t option) {fOptionScale = option; fMatchingIsDone = kFALSE;}
  void SetLSmethod(Int_t option) {fOptionLSmethod = option; fMatchingIsDone = kFALSE;}
  void SetWeightedAveragePower(Double_t power) {fWeightedAveragePower = power; fMatchingIsDone = kFALSE;}
  void SetMinuitFitOption(Float_t option) {fOptionMinuit = option; fMatchingIsDone = kFALSE;}
  void SetUseSignificantZero(Bool_t option) {fgOptionUseSignificantZero = option; fMatchingIsDone = kFALSE;}
  void SetScaleSummedBkg(Bool_t option) {fOptionScaleSummedBkg = option; fMatchingIsDone = kFALSE;}
  void SetDebugMode(Bool_t option) {fOptionDebug=option; fMatchingIsDone = kFALSE;}
  
  // set various ranges
  void SetMassFitRange(Double_t min, Double_t max) {fgMassFitRange[0] = min+1.0e-6; fgMassFitRange[1] = max-1.0e-6; fUserEnabledMassFitRange = kTRUE; fMatchingIsDone = kFALSE;}
  void SetPtFitRange(Double_t min, Double_t max) {fgPtFitRange[0] = min+1.0e-6; fgPtFitRange[1] = max-1.0e-6; fUserEnabledPtFitRange = kTRUE; fMatchingIsDone = kFALSE;}
  void AddMassExclusionRange(Double_t min, Double_t max) {
     if(fgNMassExclusionRanges==10) return;        // maximum 10 mass exclusion ranges
     fgMassExclusionRanges[fgNMassExclusionRanges][0] = min + 1.0e-6; fgMassExclusionRanges[fgNMassExclusionRanges][1] = max -1.0e-6;
     fgNMassExclusionRanges++;
     fMatchingIsDone = kFALSE;
  }
  void SetMassExclusionRange(Double_t min, Double_t max) {
     fgMassExclusionRanges[0][0] = min + 1.0e-6; fgMassExclusionRanges[0][1] = max -1.0e-6;
     fgNMassExclusionRanges = 1;
     fMatchingIsDone = kFALSE;
  }
  void SetResidualFitFunction(TF1* fitFunc) {fResidualFitFunc = (TF1*)fitFunc->Clone("ResidualFitFunction");}
  void SetBkgFitFunction(TF1* fitFunc, Bool_t forceNew=kFALSE) {
     if(!forceNew && fBkgFitFunction) return;
     if(fBkgFitFunction) delete fBkgFitFunction;
     fBkgFitFunction = (TF1*)fitFunc->Clone("BkgFitFunction");
  }
  
  Bool_t Process();
  Double_t* ComputeOutputValues(Double_t minMass, Double_t maxMass, Double_t minPt=-1., Double_t maxPt=-1.);
  void Print();        // print a summary of all user options
  void PrintFitValues(); 
  
  /*void SetEffHistogram(TH2D* eff)  {fEffVsPtCent = eff;}
  void SetWeightHistogram(TH1D* weights)  {fWeightVsCent = weights;}
  void SetEventsHistogram(TH1F* events) {fEventVsCent = events;}
  
  */
  // Getters
  TH1* GetSplusB() const {return (fMatchingIsDone ? fSplusB : 0x0);}
  TH1* GetBkg() const {return (fMatchingIsDone ? fBkg : 0x0);}
  TH1* GetSignal() const {return (fMatchingIsDone ? fSig : 0x0);}
  TH1* GetSoverB(Bool_t fromMCshape=kFALSE) const {return (fMatchingIsDone ? (fromMCshape ? fSoverBfromMCshape : fSoverB) : 0x0);}
  TH1* GetSplusResidualBkg() const {return (fMatchingIsDone ? fSplusResidualBkg : 0x0);}
  TH1* GetBkgCombinatorial() const {return (fMatchingIsDone ? fBkgCombinatorial : 0x0);}
  TH1* GetResidualBkg() const {return (fMatchingIsDone ? fBkgResidual : 0x0);}
  TH1* GetSignalMC() const {return (fMatchingIsDone ? fSignalMCshape : 0x0);}
  
  Int_t GetBkgMethod() const {return fOptionBkgMethod;}
  Int_t GetScalingOption() const {return fOptionScale;}
  Int_t GetMEMatchingMethod() const {return fgOptionMEMatching;}
  Int_t GetMinuitFitOption() const {return fOptionMinuit;}
  Double_t* GetMassFitRange() const {return fgMassFitRange;}
  Int_t     GetNMassExclusionRanges() const {return fgNMassExclusionRanges;}
  Double_t* GetMassExclusionRange(Int_t i=0) const {return (i<fgNMassExclusionRanges ? fgMassExclusionRanges[i] : 0x0);}
  const Double_t* GetFitValues() const {return fFitValues;}
  TF1* GetResidualFitFunction() const {return fResidualFitFunc;}
  TF1* GetBkgFitFunction() const {return fBkgFitFunction;}
  TF1* GetGlobalFitFunction() const {return fGlobalFitFunction;}
  Bool_t GetDebugMode() const {return fOptionDebug;}
   
 private:
   // User input data ------------------------------------------------------------------------------------------------
   THnF*    fSEOS;
   THnF*    fSELSleg1;
   THnF*    fSELSleg2;
   THnF*    fMEOS;
   THnF*    fMELSleg1;
   THnF*    fMELSleg2;
   THnF*    fSEOS_MCtruth;                          // signal histogram
  
   Int_t  fNVariables;                                                                           // number of variables to be handled
   Int_t  fVariables[kNMaxVariables];                                                  // list of variables
   Double_t fVarLimits[kNMaxVariables][2];                                       // variable limits used for integrating and fitting the mass (or mass-pt) distribution
                                                                                                          // NOTE: these limits are the most inclusive, such that both signal counting, 
                                                                                                          //          plotting and fit ranges are included
   Int_t fVarBinLimits[kNMaxVariables][2];
   Int_t  fVarIndices[kNMaxVariables];                                                // indices of variables in the THnF
  
   Int_t fMassVariable;             // the mass variable among the fNVariables  (defaults to AliReducedVarManager::kMass)
   Int_t fPtVariable;                  // the transverse momentum variable among the fNVariables (defaults to -1 -> not set)
  
   // Temporary variables
   Int_t fNLoopingVariables;     
   Int_t fCurrentVariable;
   Int_t fIter[kNMaxVariables];
   static TH1*  fgTempSignal;             // pointer to temporary signal histogram used during fitting
   static TH1*  fgTempBkg;             // pointer to temporary bkg histogram used during fitting
  
   // User options --------------------------------------------------------------------------------------------------
   static Bool_t     fgOptionUse2DMatching;        // FALSE (default): match inv.mass projections;  TRUE: match (m,pt) projections
             Int_t        fOptionBkgMethod;               // either one of these: kBkgMixedEvent (default), kBkgLikeSign, kBkgFunction
   static Int_t        fgOptionMEMatching;               // either one of these: kMatchSEOS (default), kMatchSELS
             Bool_t     fOptionUseRfactorCorrection;     // if true apply R-factor correction; default is FALSE
             Int_t        fOptionScale;                         // either one of these: kScaleEntries (default), kScaleWeightedAverage, kScaleFit
             Int_t        fOptionLSmethod;                 // either one of : kLSGeometricMean (default), kLSArithmeticMean (used for low stat situations)
             Double_t fWeightedAveragePower;    // (default is 2.0) power of the inverse statistical error used as weights for the weighted average
             Int_t        fOptionMinuit;                // either kMinuitMethodChi2 (default) or kMinuitMethodLikelihood
   static Bool_t      fgOptionUseSignificantZero;    // if true, assume zero entries as significant and error of 1 during the chi2 calculation
             Bool_t      fOptionScaleSummedBkg;       // if true, run the matching procedure on the summed S+B and bkg (default is false)
             Bool_t      fOptionDebug;                       // if true, construct all possible distributions
             
   // Matching / fit ranges
   // NOTE: Mass and pt ranges used for matching / fitting can in principle be different (a sub-interval only) wrt ranges in fVarLimits
   //            If the dedicated setter function are not called by user, these ranges will be made same as in fVarLimits at Initialize() time
   static Double_t fgMassFitRange[2];             // mass range used in the bkg to signal matching or in the fit procedure
   Bool_t     fUserEnabledMassFitRange;   // default false, enabled when SetMassFitRange() is called
   static Double_t fgPtFitRange[2];                  // pt range used in the bkg to signal matching or in the fit procedure
   Bool_t     fUserEnabledPtFitRange;        // default false, enabled when SetPtFitRange() is called
   //static Double_t fgMassExclusionRange[2];         // mass exclusion range, used in matching / fitting
   static Double_t fgMassExclusionRanges[10][2];     // mass exclusion range, used in matching / fitting
   static Int_t fgNMassExclusionRanges;                // number of mass exclusion ranges
  
   // Utility data members
   TH1* fSplusB;               //  total signal + bkg projection
   TH1* fBkg;                    //  background projection
   TH1* fSig;                    // signal projection
   
   TH1* fBkgLikeSign;
   TH1* fBkgLikeSignLeg1;
   TH1* fBkgLikeSignLeg2;
   TH1* fBkgMixedEvent;
   
   static TF1* fBkgFitFunction;
   static TF1* fGlobalFitFunction;
   TFitResultPtr fFitResult;                // fit result of the residual fit
   
   /////////////////////////////////
   TH1* fSplusResidualBkg;    //  combinatorial bkg subtracted minv distribution (signal + residual bkg)   
   TH1* fSplusBblind;       //  bkg minv distribution; signal blind (area around signal excluded)   
   TH1* fBkgCombinatorial;    //  combinatorial bkg (used when the residual bkg fit option is switched on)
   TH1* fBkgResidual;             // residual bkg obtained after fitting the combinatorial bkg subtracted distribution
   ////////////////////////////////
   
   TH1* fSoverB;              // S/B projection
   TH1* fSoverBfromMCshape;   // S/B projection using the scaled MC shape for the signal
   static TH1* fSignalMCshape;    // MC truth signal shape
   Double_t fFitValues[kNFitValues];       // array used to store information on the signal fit
   Bool_t fMatchingIsDone;                  // set to true if the matching procedure was succesfully run; false if the object is in any other state
   
   TMinuit* fMinuitFitter;                    // used if fit option is required
   ///////////////////////////////////////////////////
   TF1*      fResidualFitFunc;            // fit function used to fit the combinatorial bkg subtracted minv distribution
   
   ////////////////////////////////////////////////////
   
   // Private utility functions
   Bool_t Initialize();            // returns true if all prerequisites for signal extraction are met
   void Slice();
   void ApplyUserRanges(THnF* h);
   void AddSlice();
   TH1* BuildLSbkg(TH1* selsLeg1, TH1* selsLeg2, TH1* meos=0x0, TH1* melsLeg1=0x0, TH1* melsLeg2=0x0);
   void SqrtTH1(TH1* h, Bool_t is2D=kFALSE);
   void  ComputeEntryScale(TH1* signal, TH1* bkg);
   void  ComputeWeightedScale(TH1* sig, TH1* bkg);
   void  FitScale(TH1* sig, TH1* bkg, Bool_t fixScale=kFALSE);
   void  ComputeScale(TH1* scaleHist, TH1* bkgHist);
   static void Fcn(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t);
   static Double_t Chi2(TH1* sig, TH1* bkg, Double_t scale, Double_t scaleError=0.0);
   void FitInvMass();
   void FitResidualBkg();
   static Double_t GlobalFitFunction(Double_t *x, Double_t* par);

   ClassDef(AliResonanceFits, 6);
};

#endif
