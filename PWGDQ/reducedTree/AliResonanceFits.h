// Class used for extracting resonance yields from invariant mass distributions
// Author: Ionut-Cristian Arsene (iarsene@cern.ch)
// 
#ifndef ALIRESONANCEFITS_H
#define ALIRESONANCEFITS_H

#include <TObject.h>
#include <THn.h>

class TH1;
class TH1D;
class TH1F;
class TH2D;
class TVirtualPad;

//_____________________________________________________________________
class AliResonanceFits : public TObject {
 
 public:
  enum Variables {
    kMass=0,
    kPt,
    kRapidity,
    kCentrality,
    kVertexZ,
    kEP2,
    kNVariables
  };
  enum FitValues {
    kSig,
    kSigErr,
    kBkg,
    kBkgErr,
    kSplusB,
    kSplusBerr,
    kSoverB,
    kSoverBerr,
    kSignif,
    kSignifErr,
    kChisq,
    kChisqMC,
    kBkgScale,
    kBkgScaleErr,
    kNFitValues
  };
   
  AliResonanceFits();
  virtual ~AliResonanceFits();            

  // setters
  void Reset();
  // User input
  void SetHistograms(THnF* seos, THnF* meos, 
		     THnF* selsLeg1=0x0, THnF* selsLeg2=0x0, THnF* melsLeg1=0x0, THnF* melsLeg2=0x0);
  void SetVars(Int_t* vars);
  void SetEffHistogram(TH2D* eff)  {fEffVsPtCent = eff;}
  void SetWeightHistogram(TH1D* weights)  {fWeightVsCent = weights;}
  void SetEventsHistogram(TH1F* events) {fEventVsCent = events;}
  void SetSignalMCshape(TH1D* shape) {fSignalMCshape = shape;}
  // User options
  void SetFitRange(Float_t min, Float_t max)   {fFitRange[0]=min; fFitRange[1]=max;}
  void SetExclusionRange(Float_t min, Float_t max)  {fExclusionRange[0]=min; fExclusionRange[1]=max;}
  void SetMassRange(Float_t min, Float_t max)  {fMassRange[0]=min; fMassRange[1]=max;}
  void SetSignalRange(Float_t min, Float_t max)   {fSignalRange[0]=min; fSignalRange[1]=max;}
  void SetCentralityRange(Float_t min, Float_t max) {fCentralityRange[0]=min; fCentralityRange[1]=max; fCentralitySelection = kTRUE;}
  void SetVertexZRange(Float_t min, Float_t max) {fVertexZRange[0]=min; fVertexZRange[1]=max; fVertexSelection = kTRUE;}
  void SetEP2Range(Float_t min, Float_t max) {fEP2Range[0]=min; fEP2Range[1]=max; fEPSelection = kTRUE;}
  void SetPtRange(Float_t min, Float_t max) {fPtRange[0]=min; fPtRange[1]=max; fPtSelection = kTRUE;}
  
  void SetBkgMethod(Int_t method) {fBkgMethod = method;}
  void SetUse2DMatching(Bool_t flag=kTRUE) {fUse2DMatching = flag;}
  void SetMatchingOption(Int_t option) {fMatchingOption = option;}
  void SetWeightedAveragePower(Float_t power) {fWeightedAveragePower = power;}
  void SetMinuitFitOption(Float_t option) {fMinuitFitOption = option;}
  void SetFixBkgScale(Float_t scale) {fFixScale=kTRUE; fFitValues[kBkgScale]=scale;}
  void SetMEMatchingBkg(Int_t option) {fMEMatchOption = option;}
  void SetLSmethod(Int_t method) {fLSmethod = method;}
  void SetPlotingOption(Int_t option) {fPlottingOption=option;}
  
  // Getters
  TH1* GetSignal() const {if(!fIsProcessed) return 0x0; return fHistBkgSubtracted;}
  TH1* GetSplusB() const {if(!fIsProcessed) return 0x0; return fHistSEOS;}
  TH1* GetMEbkg()  const {if(!fIsProcessed) return 0x0; return fHistMEbkg;}
  TH1* GetLSbkg()  const {if(!fIsProcessed) return 0x0; return fHistSELSbkg;}
  TH1* GetBkg()    const {if(!fIsProcessed) return 0x0; if(fBkgMethod==1) return fHistMEbkg; else return fHistSELSbkg;}
  const Double_t* GetFitValues() const {if(!fIsProcessed) return 0; return fFitValues;}
  
  // Main working functions
  void Process();
  void Process(Int_t bkgMethod, Int_t matchOption, Float_t waPower, Float_t minuitOption, Float_t fixBkgScale, 
	       Int_t meMatchBkg, Int_t lsMethod);
  void ExtractSignal(TH1* signal, TH1* bkg, Bool_t fixScale=kFALSE);
  void DrawSignalExtraction(Bool_t save=kFALSE, const Char_t* name="", const Char_t* outputDir="",
                            Bool_t makeNicePlot=kFALSE, 
                            TVirtualPad* externalPad=0x0, Bool_t noYlabels=kFALSE, Bool_t noLegends=kFALSE);
  
  // Utility stuff  
  
 private:
  // User input data ------------------------------------------------------------------------------------------------
  THnF*    fSEOS;
  THnF*    fSELSleg1;
  THnF*    fSELSleg2;
  THnF*    fMEOS;
  THnF*    fMELSleg1;
  THnF*    fMELSleg2;
  Bool_t   fUsedVars[kNVariables];
  Int_t    fVarIndex[kNVariables];
  TH2D*    fEffVsPtCent;             // Efficiency vs centrality and pt to be applied to the invariant mass spectrum 
                                     // TODO: Add more dimensions (e.g. pt)
  TH1D*    fWeightVsCent;            // Weights vs centrality to be applied to the invariant mass spectrum
                                     // Needed for example when the centrality distribution in data is not flat
  TH1F*    fEventVsCent;             // Distribution of events over centrality                                   
  TH1D*    fSignalMCshape;           // Signal shape from MC
  
  // User options --------------------------------------------------------------------------------------------------
  static Float_t  fFitRange[2];             // mass range used for the sig extraction procedure
  static Float_t  fExclusionRange[2];       // mass range excluded from matching (e.g. range around resonance peak)
  Float_t  fMassRange[2];            // mass range used for plotting
  Float_t  fSignalRange[2];          // mass range used to count the signal
  Float_t  fCentralityRange[2];
  Bool_t   fCentralitySelection;
  Float_t  fVertexZRange[2];
  Bool_t   fVertexSelection;
  Float_t  fEP2Range[2];
  Bool_t   fEPSelection;
  static Float_t fPtRange[2];
  Bool_t   fPtSelection;
  
  Int_t    fPlottingOption;          // (default is 0) 0 - mass projection; 1 - pt projection; 2 - (mass,pt) projection
  Int_t    fBkgMethod;               // (default is 1) 1 - Mixed event; 2 - Like sign
  static Bool_t fUse2DMatching;      // (default is false - match inv.mass projections)
                                     // if toggled use (m,pt) projections for matching
  Int_t    fMatchingOption;          // (default is 1) 1 - weighted average; 2 - fit; 3 - signal entries
  Float_t  fWeightedAveragePower;    // (default is 2.0) power of the inverse statistical error used as weights for the weighted average
  Float_t  fMinuitFitOption;         // (default is 1.0) fit option for TMinuit: 1.0 - chisquare minimization; 0.5 - log-likelihood
  Bool_t   fFixScale;                // (default is kFALSE) Fix the scale in the Minuit fit (the point is just to get the chisquare)
  Int_t    fMEMatchOption;           // (default is 2)
                                     // 1 - match to the SEOS outside signal region; 
                                     // 2 - match to the SELS R-factor corrected bkg
  Int_t    fLSmethod;                // (default is 1) 1 - arithmetic mean; 2 - geometric mean
    
  // Utility data members ---------------------------------------------------------------------------------------------
  Bool_t       fIsProcessed;
  static Int_t fNdf;
  Double_t     fFitValues[kNFitValues];
  TH1*         fHistBkgSubtracted;      // TH2D/TH1D depending on whether fUse2DMatching is toggled or not 
  TH1*         fHistSEOS;               // TH2D/TH1D depending on whether fUse2DMatching is toggled or not
  TH1*         fHistSELSbkg;            // TH2D/TH1D depending on whether fUse2DMatching is toggled or not
  TH1*         fHistMEbkg;              // TH2D/TH1D depending on whether fUse2DMatching is toggled or not
  static TH1*  fTempHistSignal;
  static TH1*  fTempHistBkgnd;
  
  // private utility functions
  static Double_t Chi2(TH1* signal, TH1* background, Double_t scale, Double_t scaleError=0.0);
  static void Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
  void     ComputeWeightedScale(TH1* signal, TH1* bkg);
  void     FitScale(TH1* signal, TH1* bkg);
  void     ComputeEntryScale(TH1* signal, TH1* bkg);
  void     MakeSEOS();
  void     MakeLSbkg();
  void     MakeMEbkg();
  TH1*     GetLScorrected(TH1* selsLeg1, TH1* selsLeg2, TH1* meos, TH1* melsLeg1, TH1* melsLeg2);
  void     SetDefaultRanges();
  void     DrawMassProjection(TH1* seos, TH1* bkg, TH1* signal, Bool_t save=kFALSE, const Char_t* name="", const Char_t* outputDir="",
                              Bool_t makeNicePlot=kFALSE, TVirtualPad* externalPad=0x0, Bool_t noYlabels=kFALSE, Bool_t noLegends=kFALSE);
  
  ClassDef(AliResonanceFits, 2);
};

#endif
