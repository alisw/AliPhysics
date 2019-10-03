// Class for advanced fitting methods (template fits, (weighted) LL fits, regularised fits, simultaneous fits)
//
// Based on original code of: 
// Xianguo Lu
// lu@physi.uni-heidelberg.de 
//
// Modified by:
// Benjamin Hess
// Benjamin-Andreas.Hess@Uni-Tuebingen.de

#ifndef ALITPCPIDMATHFIT_H
#define ALITPCPIDMATHFIT_H

class TString;
class TMinuit;
class TH1;

class AliTPCPIDmathFit: public TObject
{
public:
  typedef Double_t (*FitFunc_t)(const Double_t *, const Double_t *);
  
  static AliTPCPIDmathFit* Instance(const Int_t numXbinsRegularisation = -1, const Int_t numSimultaneousFits = -1,
                                    const Int_t maxDataPoints = 200000);
  
  Int_t AddRefHisto(TH1* refHisto);
  void ClearRefHistos();
  
  Bool_t GetApplyPatching() const { return fApplyPatching; };
  Int_t GetDebugLevel() const { return fDebugLevel; }
  Double_t GetEpsilon() const { return fEpsilon;  }
  Int_t GetNrefHistos() const { return fNrefHistos;  }
  TH1* GetRefHisto(Int_t index) const;
  Int_t GetIndexParametersToRegularise(Int_t index) const;
  Int_t GetMaxCalls() const { return fMaxCalls; }
  const TString GetMinimisationStrategy() const { return fMinimisationString; }; 
  Int_t GetNumParametersPerXbin() const { return fNumParametersPerXbin; };
  Int_t GetNumParametersToRegularise() const { return fNumParametersToRegularise; };
  Int_t GetNumSimultaneousFits() const { return fNumSimultaneousFits; };
  Int_t GetNumXbinsRegularisation() const { return fNumXbinsRegularisation; };
  const Double_t* GetParAdditional() const { return fParAdditional; };
  Int_t GetXbinIndex() const { return fXbinIndex; };
  const Double_t* GetXstatisticalWeight() const { return fXstatisticalWeight; };
  const Double_t* GetXstatisticalWeightError() const { return fXstatisticalWeightError; };
  const Double_t* GetXstatisticalWeight2() const { return fXstatisticalWeight2; };
  const Double_t* GetXstatisticalWeightError2() const { return fXstatisticalWeightError2; };
  const Double_t* GetXvaluesForRegularisation() const { return fXvaluesForRegularisation; };
  Int_t GetRegularisation() const { return fRegularisation; };
  Double_t GetRegularisationFactor() const { return fRegularisationFactor; };
  Double_t GetScaleFactorError() const { return fScaleFactorError; };
  Bool_t GetUseLogLikelihood() const { return fUseLogLikelihood; }
  Bool_t GetUseRegularisation() const { return fRegularisation > 0 && fNumParametersToRegularise > 0 && fNumXbinsRegularisation > 1; };
  Bool_t GetUseWeightsForLogLikelihood() const { return fUseWeightsForLoglikelihood; }
  
  void InputData(const TH1 *hh, const Int_t indexXbinRegularisation, const Int_t indexSimultaneousFit, const Double_t leftBoundary,
                 const Double_t rightBoundary, const Double_t threshold, const Bool_t kXerr);
  Int_t MinuitFit(FitFunc_t *fn, FitFunc_t *fnError, const Int_t nPar, Double_t *par, Double_t *err, Double_t *covMatrix, 
                  Double_t &ChiSquareOrLogLikelihood, Int_t &NDF, const Double_t *stepSize = 0x0, const Double_t *lowLim = 0x0, 
                  const Double_t *hiLim = 0x0);
  
  void SetApplyPatching(const Int_t applyPatching) { fApplyPatching = applyPatching; };
  void SetDebugLevel(const Int_t level = 1) { fDebugLevel = level; }
  void SetEpsilon(const Double_t eps);
  void SetMaxCalls(const Int_t maxCalls);
  void SetMinimisationStrategy(TString strategy) { fMinimisationString = strategy; };
  Bool_t SetParametersToRegularise(const Int_t numParams, const Int_t numParamsPerXbin, const Int_t* indexParametersToRegularise,
                                   const Int_t* lastNotFixedIndexOfParameters, const Double_t* xValuesForRegularisation,
                                   const Double_t* xStatisticalWeight, const Double_t* xStatisticalWeightError,
                                   const Double_t* xStatisticalWeight2, const Double_t* xStatisticalWeightError2,
                                   const Double_t* parAdditional);
  void SetRegularisation(const Int_t reg, const Double_t regFactor) { fRegularisation = reg; fRegularisationFactor = regFactor; };
  void SetScaleFactorError(const Double_t scaleFactorError) { fScaleFactorError = scaleFactorError; };
  void SetUseLogLikelihood(const Bool_t useLogLikelihood = kTRUE) { fUseLogLikelihood = useLogLikelihood; }
  void SetUseWeightsForLogLikelihood(const Bool_t useWeightsForLogLikelihood = kTRUE) 
    { fUseWeightsForLoglikelihood = useWeightsForLogLikelihood; }
  
 private:
  AliTPCPIDmathFit(const Int_t numXbinsRegularisation = 1, const Int_t numSimultaneousFits = 1, const Int_t maxDataPoints = 200000);
  virtual ~AliTPCPIDmathFit();
  
  void AddPenaltyTermForRegularisation(const Bool_t useLogLikelihood);
  inline Double_t EvalLog(Double_t x) const;
  Int_t Fit(const Double_t *inPar = 0x0, Double_t *covMatrix = 0x0, const Double_t *stepSize = 0x0, const Double_t *lowLim = 0x0,
            const Double_t *hiLim = 0x0);
  static void FitFCN(Int_t &nPar, Double_t *input, Double_t &chiSquareOrLogLikelihood, Double_t *par, Int_t flag);
  Double_t GetChiSquare();
  Double_t GetNegLogLikelihood();
  void GetCovarianceMatrix(Double_t* covMatrix) const;
  inline static Double_t PatchParameter(Double_t par, Double_t statWeightTot1, Double_t statWeightPar2, Double_t statWeightTotSummed);
  Bool_t PolynomialInterpolation(Double_t xa[], Double_t ya[], Int_t n, Double_t x, Double_t* y, Double_t* dy);

  //______________________________________________
  //______________________________________________
  
  static AliTPCPIDmathFit* fgInstance;   //  Instance of this class (singleton implementation)
  
  const Double_t fkDoubleEpsilonLimit; // Limit for comparing double with zero
  
  //------ Initialized in MinuitFit(), func, nPar, par and err passed in from outside ------
  //------ all reset after Fit() ------
  FitFunc_t *fFunc;      // Function used for fitting
  FitFunc_t *fErrorFunc; // Function used for error estimation while fitting
  Int_t fNpar;           // Number of fit parameters
  Double_t *fPar;        // Fit parameters
  Double_t *fErr;        // Error of fit parameters
  
  //chi2, fNegLogLikelihood, ndf only used internally, since if it stays, it may be confused by different MinuitFit.
  //Reset after Fit() and only passed out after a Fit()
  Double_t fChi2;             // Chi^2 of fit
  Double_t fNegLogLikelihood; // Negative log-likelihood of fit
  Int_t fNDF;                 // NDF of fit

  Int_t fNstep;               // Current fitting step
  
  //------ initialized via InputData() ------
  Int_t **fNdata;             // Array containing number of data points
  
  Double_t ***fX;             // Array containing x values of data points
  Double_t ***fY;             // Array containing y values of data points
  
  Double_t ***fXerr;          // Array containing x errors of data points
  Double_t ***fYerr;          // Array containing y errors of data points
  
  //------ Set via setter functions -----
  Int_t fMaxCalls;            // Maximum number of calls for fitting
  Double_t fEpsilon;          // Epsilon value for fit convergence
  Bool_t fUseLogLikelihood;                          // Flag indicating whether to use LL fit or chi^2 fit
  Bool_t fUseWeightsForLoglikelihood;                // Flag indicating whether to use weighted LL fit or normal one (global)
  Bool_t fUseWeightsForLoglikelihoodForCurrentFit;   // Flag indicating whether to use weighted LL fit or normal one for the current fit
  
  Double_t fScaleFactorError; // Scale errors with this factor
  
  Int_t fRegularisation;      // Do regularisation with +/- this number of bins 
  Int_t fXbinIndex;           // Index of current x bin
  
  Int_t fNumParametersToRegularise;       // Number of fit parameter that will be regularised
  Int_t fNumParametersPerXbin;            // Number of fit parameters for each x bin
  Int_t *fIndexParametersToRegularise;    // Indicices of fit parameters that are to be regularised
  Int_t *fLastNotFixedIndexOfParameters;  // Indices of last free parameter in the corresponding x bin
  
  Double_t *fXvaluesForRegularisation;    // Array with x values that will be used for regularisation
  
  Double_t *fXstatisticalWeight;          // Statistical weights of each x bin
  Double_t *fXstatisticalWeightError;     // Errors of statistical weights of each x bin
  
  Bool_t    fApplyPatching;               // Turn on/off patching of parameters with following additional weights 
  Double_t *fParAdditional;               // Additional contribution to fit parameters
  Double_t *fXstatisticalWeight2;         // Statistical weights of each x bin - additional weight from other source
  Double_t *fXstatisticalWeightError2;    // Errors of statistical weights of each x bin - additional weight error from other source
  
  Double_t fRegularisationFactor;         // Relative weighting factor for regularisation (1 = weighted with error)
  
  TString fMinimisationString;            // Minimisation method

  //------ Initialized via constructor -----
  const Int_t fkMaxNdata;                 // Maximum number of data points
  Int_t fNumSimultaneousFits;             // Number of simultaneous fits
  Int_t fNumXbinsRegularisation;          // Number of x bins for the regularisation
  
  //------ initialized and destroyed in Fit()
  TMinuit *fMinuit;                       // Pointer to minuit object
  
  Double_t *fPolIntX;                     // Used for internal interpolation during regularisation
  Double_t *fPolIntY;                     // Used for internal interpolation during regularisation
  Double_t *fPolIntC;                     // Used for internal interpolation during regularisation
  Double_t *fPolIntD;                     // Used for internal interpolation during regularisation
  
  //------
  Int_t fDebugLevel;                      // Debug level
  
  
  // Reference histos used for the fitting
  TH1** fRefHistos;                       // Array with histos (templates) used for fitting
  Int_t fNrefHistos;                      // Number of histos (templates) for fitting
  
  AliTPCPIDmathFit(const AliTPCPIDmathFit&);            // Not implemented
  AliTPCPIDmathFit& operator=(const AliTPCPIDmathFit&); // Not implemented
  
  ClassDef(AliTPCPIDmathFit, 1);
};

#endif
