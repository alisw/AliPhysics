#include "TH1.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TMinuit.h"
#include "TStopwatch.h"
#include "TString.h"

#include "AliLog.h"
#include <iostream>
#include <limits>

#include "AliTPCPIDmathFit.h"

// Class for advanced fitting methods (template fits, (weighted) LL fits, regularised fits, simultaneous fits)
//
// Based on original code of: 
// Xianguo Lu
// lu@physi.uni-heidelberg.de 
//
// Modified by:
// Benjamin Hess
// Benjamin-Andreas.Hess@Uni-Tuebingen.de


ClassImp(AliTPCPIDmathFit);

AliTPCPIDmathFit* AliTPCPIDmathFit::fgInstance = NULL;

//***Public functions***
//______________________________________________________
Int_t AliTPCPIDmathFit::AddRefHisto(TH1* refHisto)
{
  // Add refHisto to the list of reference histos. The index of the added histo is returned.

  TH1** newRefHistos = new TH1*[fNrefHistos + 1];
  for (Int_t i = 0; i < fNrefHistos; i++)   {
    newRefHistos[i] = (TH1*)fRefHistos[i]->Clone();
    newRefHistos[i]->SetName(Form("%s_%d", newRefHistos[i]->GetName(), i));
    newRefHistos[i]->SetDirectory(0);
    delete fRefHistos[i];
  }
  
  newRefHistos[fNrefHistos] = (TH1*)refHisto->Clone();
  
  delete [] fRefHistos;
  fRefHistos = newRefHistos;
  
  return fNrefHistos++;
}


//______________________________________________________
void AliTPCPIDmathFit::ClearRefHistos()
{
  // Clear list of reference histos
  
  for (Int_t i = 0; i < fNrefHistos; i++)
    delete fRefHistos[i];
  
  delete [] fRefHistos;
  fRefHistos = 0x0;
  fNrefHistos = 0;
}


//______________________________________________________
TH1* AliTPCPIDmathFit::GetRefHisto(Int_t index) const
{
  if (index < 0 || index >= fNrefHistos || !fRefHistos) 
    return 0x0;
  
  return fRefHistos[index];
}


//______________________________________________________
Int_t AliTPCPIDmathFit::GetIndexParametersToRegularise(Int_t index) const
{
  if (index < 0 || index >= fNumParametersToRegularise || !fIndexParametersToRegularise)
    return -1;
  
  return fIndexParametersToRegularise[index];
}


//______________________________________________________
AliTPCPIDmathFit* AliTPCPIDmathFit::Instance(const Int_t numXbinsRegularisation, const Int_t numSimultaneousFits,
                                             const Int_t maxDataPoints)
{
  // Get instance of AliTPCPIDmathFit. If numSimultaneousFits is > 0 and there is already an existing instance,
  // the instance will be deleted and a new one with the new number of simultaneous fits will be created.
  // The same applies applies to numXbinsRegularisation.
  
  if(!fgInstance) {
    fgInstance = new AliTPCPIDmathFit(numXbinsRegularisation, numSimultaneousFits, maxDataPoints);
  }

  if (numSimultaneousFits > 0 && numSimultaneousFits != fgInstance->GetNumSimultaneousFits()) {
    printf("Different number of simultaneous fits. Creating new instance with %d instead of %d simultaneous fits...\n",
           numSimultaneousFits, fgInstance->GetNumSimultaneousFits());
    delete fgInstance;
    fgInstance = new AliTPCPIDmathFit(numXbinsRegularisation, numSimultaneousFits, maxDataPoints);
  }
  
  if (numXbinsRegularisation > 0 && numXbinsRegularisation != fgInstance->GetNumXbinsRegularisation()) {
    printf("Different number of x bins for regularisation. Creating new instance with %d instead of %d such bins...\n",
           numXbinsRegularisation, fgInstance->GetNumXbinsRegularisation());
    delete fgInstance;
    fgInstance = new AliTPCPIDmathFit(numXbinsRegularisation, numSimultaneousFits, maxDataPoints);
  }
   
  return fgInstance;
}


//______________________________________________________
void AliTPCPIDmathFit::InputData(const TH1 *hh, const Int_t indexXbinRegularisation, const Int_t indexSimultaneousFit,
                                 const Double_t leftBoundary, const Double_t rightBoundary, const Double_t threshold,
                                 const Bool_t kXerr)
{
  // Take index of simultaneousFit and xBinRegularisation and fill the corresponding data points from the histo
  
  if (indexSimultaneousFit < 0 || indexSimultaneousFit >= fNumSimultaneousFits) {
    printf("Error: Index %d of simultaneous fit does not exist. Must be 0 - %d!\n", indexSimultaneousFit, fNumSimultaneousFits - 1);
    return;
  }
  
  if (indexXbinRegularisation < 0 || indexXbinRegularisation >= fNumXbinsRegularisation) {
    printf("Error: Index %d of x bin regularisation does not exist. Must be 0 - %d!\n", indexXbinRegularisation,
           fNumXbinsRegularisation - 1);
    return;
  }
  
  fNdata[indexXbinRegularisation][indexSimultaneousFit] = 0;
  
  for (Int_t i = 0; i < fkMaxNdata; i++) {
    fX[indexXbinRegularisation][indexSimultaneousFit][i] = -999;
    fY[indexXbinRegularisation][indexSimultaneousFit][i] = -999;
    fXerr[indexXbinRegularisation][indexSimultaneousFit][i] = -999;
    fYerr[indexXbinRegularisation][indexSimultaneousFit][i] = -999;
  }

  //------
  
  //const Double_t hmax =  hh->GetBinContent(hh->GetMaximumBin());
  const TAxis * ax = hh->GetXaxis();
  
  // Only fill bins inside the desired range and only bins with content > threshold
  for (Int_t id = ax->GetFirst(); id <= ax->GetLast(); id++) {
    if (hh->GetBinCenter(id) < leftBoundary)  continue;
    if (hh->GetBinCenter(id) > rightBoundary)  break;
     
    const Double_t yy = hh->GetBinContent(id);
    if (yy <= threshold)//"=" important to exclude 0-bins
        continue;

    fX[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = hh->GetBinCenter(id);
    fY[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = yy;

    fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = hh->GetBinError(id);
    
    // Small yErr only matters for chi^2 fit. For likelihood fit, just set error to 1, if no error is set.
    if (fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] 
        < fkDoubleEpsilonLimit) {
      if (fUseLogLikelihood) {
        fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = 1.0;
      }
      else if (fY[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] 
               < fkDoubleEpsilonLimit) {
        Int_t firstNonZeroBin = hh->FindFirstBinAbove(0.0);
        Int_t lastNonZeroBin = hh->FindLastBinAbove(0.0);
        
        // Bins with small content usually sit at the edges of the histo
        Double_t errorOfBinsWithSmallestContent = 0.0;
        if (firstNonZeroBin > 0 && lastNonZeroBin > 0) {
          errorOfBinsWithSmallestContent = (hh->GetBinError(firstNonZeroBin) + hh->GetBinError(lastNonZeroBin)) / 2.0;
          
          
          if (fDebugLevel > 0)
            printf("MathFit::InputData strange yerr %e for y %e -> Setting arbitrary error %e!\n", 
                  fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]],
                  fY[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]],
                  errorOfBinsWithSmallestContent);
            
          fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] =
            errorOfBinsWithSmallestContent;
        }
        else {
          // There is no information, exclude from fit
          fX[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = -999;
          fY[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = -999;
          fXerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = -999;
          fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = -999;
          continue;
        }

        //// For empty bins assume bin content = 1, i.e. error = 1
        //fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = 1.0;
        /*
        const Double_t ww = 10.;
        fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = hmax / ww;
        //printf("MathFit::InputData yerr asigned arbitrary as %f = hamx %f / %f for 0-bin!!\n", 
                 fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]], hmax, ww);
        */
      }
      else {
        printf("AliTPCPIDmathFit::InputData fYerr = %e for fY = %e and fNdata = %d\n",
               fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]],
               fY[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]], fNdata[indexXbinRegularisation][indexSimultaneousFit]);
        exit(1);
      }
    }
    /*
    if (abs(fYerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]]) 
        < std::numeric_limits<double>::min()){
      printf("AliTPCPIDmathFit::InputData fYerr = 0!! %d\n", id);
      exit(1);
    }*/
    
    if (kXerr) 
      fXerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = 
        hh->GetBinWidth(id) / 2.0;
    else       
      fXerr[indexXbinRegularisation][indexSimultaneousFit][fNdata[indexXbinRegularisation][indexSimultaneousFit]] = 0.0;

    fNdata[indexXbinRegularisation][indexSimultaneousFit]++;
    if(fNdata[indexXbinRegularisation][indexSimultaneousFit] >= fkMaxNdata) {
      printf("AliTPCPIDmathFit::InputData fNdata >= fkMaxNdata!! %d\n", fkMaxNdata);
      exit(1);
    }
  }  
}


//______________________________________________________
Int_t AliTPCPIDmathFit::MinuitFit(FitFunc_t *fn, FitFunc_t *fnError, const Int_t nPar, Double_t *par, Double_t *err,
                                  Double_t *covMatrix, Double_t &ChiSquareOrLogLikelihood, Int_t &NDF, 
                                  const Double_t *stepSize, const Double_t *lowLim, const Double_t *hiLim)
{
  // Setup minuit object with in put parameters and internally stored things (like data points) and do the fitting.
  
  Bool_t hasData = kFALSE;
  for (Int_t indexXbinRegularisation = 0; indexXbinRegularisation < fNumXbinsRegularisation; indexXbinRegularisation++) {
    for (Int_t indexSimultaneousFit = 0; indexSimultaneousFit < fNumSimultaneousFits; indexSimultaneousFit++) {
      if (fNdata[indexXbinRegularisation][indexSimultaneousFit] > 0) {
        hasData = kTRUE;
        break;
      }
    }
  }
  
  if (!hasData) {
    printf("AliTPCPIDmathFit::MinuitFit no data yet, call InputData() first!!\n");
    return -1;
  }
  
  if (fRegularisation > 0 && fNumSimultaneousFits <= 1) {
    printf("Regularisation only implemented for simultaneous fit!\n");
    return -1;
  }
  
  for (Int_t indexSimultaneousFit = 0; indexSimultaneousFit < fNumSimultaneousFits; indexSimultaneousFit++) {
    fFunc[indexSimultaneousFit] = fn ? fn[indexSimultaneousFit] : 0x0;
    fErrorFunc[indexSimultaneousFit] = fnError ? fnError[indexSimultaneousFit] : 0x0;
  }
    
  fNpar = nPar;
  fPar = par;
  fErr = err;
  
  const Double_t *inpar = fPar;
  Int_t errFlag = 0;
    
  // Start the stopwatch. If reset is kTRUE reset the stopwatch before starting it (including the stopwatch counter). 
  // Use kFALSE to continue timing after a Stop() without resetting the stopwatch.
  TStopwatch stwa;
  stwa.Start(kTRUE);
  errFlag = Fit(inpar, covMatrix, stepSize, lowLim, hiLim);
  stwa.Stop();
  stwa.Print("m");//"m": milli sec; "u": micro sec
  
  //Store Chi^2/Loglikelihood and NDF
  if (fUseLogLikelihood)
    ChiSquareOrLogLikelihood = fNegLogLikelihood;
  else
    ChiSquareOrLogLikelihood = fChi2;
  NDF = fNDF;

  // Reset internal variables after fit
  for (Int_t indexSimultaneousFit = 0; indexSimultaneousFit < fNumSimultaneousFits; indexSimultaneousFit++) {
    fFunc[indexSimultaneousFit] = 0x0;
    fErrorFunc[indexSimultaneousFit] = 0x0;
  }
  
  fNpar = 0;
  fPar = 0x0;
  fErr = 0x0;

  fNstep = 0;

  fChi2 = 0;
  fNegLogLikelihood = 0;
  fNDF = 0;
  
  return errFlag;
}


//______________________________________________________
void AliTPCPIDmathFit::SetEpsilon(const Double_t eps)  
{ 
  // Set epsilon defining convergence of fit
  
  if (eps > 0)
    fEpsilon = eps;
  else 
    AliError("epsilon must be > 0");
}


//______________________________________________________
void AliTPCPIDmathFit::SetMaxCalls(const Int_t maxCalls)  
{ 
  // Set maximum number of calls for fitting
  
  if (maxCalls > 0)
    fMaxCalls = maxCalls;
  else 
    AliError("maxCalls must be > 0");
}


//______________________________________________________
Bool_t AliTPCPIDmathFit::SetParametersToRegularise(const Int_t numParams, const Int_t numParamsPerXbin, 
                                                   const Int_t* indexParametersToRegularise, const Int_t* lastNotFixedIndexOfParameters,
                                                   const Double_t* xValuesForRegularisation, const Double_t* xStatisticalWeight,
                                                   const Double_t* xStatisticalWeightError, const Double_t* xStatisticalWeight2,
                                                   const Double_t* xStatisticalWeightError2, const Double_t* parAdditional)
{
  // Sets the number of parameters to the new value and saves the indices of the parameters to regularise in an internal array.
  // Also, the number of parameters per xBin is set, the x values of the individual x bins and the last index of the x bin
  // in which the corresponding parameter is not fixed plus the statistical weights and errors for each xBin.
  // kFALSE is returned and the internal values are reset, if invalid values are provided.
  // If patching is to be used (and it must be set before), the corresponding statistical weights and parAddional can be forwarded here.
  // Note that parAdditional is the statistical weight per parameter. The size of parAdditional must be equal to the number of parameters 
  // used in the fit later!
  
  fNumParametersPerXbin = 0;
  fNumParametersToRegularise = 0;
  
  delete [] fIndexParametersToRegularise;
  fIndexParametersToRegularise = 0x0;
  
  delete [] fLastNotFixedIndexOfParameters;
  fLastNotFixedIndexOfParameters = 0x0;
  
  delete [] fXvaluesForRegularisation;
  fXvaluesForRegularisation = 0x0;
  
  delete [] fXstatisticalWeight;
  fXstatisticalWeight = 0x0;
  
  delete [] fXstatisticalWeightError;
  fXstatisticalWeightError = 0x0;
  
  delete [] fXstatisticalWeight2;
  fXstatisticalWeight2 = 0x0;
  
  delete [] fXstatisticalWeightError2;
  fXstatisticalWeightError2 = 0x0;
  
  delete [] fParAdditional;
  fParAdditional = 0x0;
  
  if (!indexParametersToRegularise || !xValuesForRegularisation || numParamsPerXbin <= 0 || numParams <= 0
      || !xStatisticalWeight || !xStatisticalWeightError || !lastNotFixedIndexOfParameters) {
    printf("Error: Cannot set parameters to regularise. Invalid input!\n");
    
    return kFALSE;
  }
  
  if (fApplyPatching && (!xStatisticalWeight2 || !xStatisticalWeightError2 || !parAdditional)) {
    printf("Error: Patching is switched on, but input is missing!\n");
    
    return kFALSE;
  }
  
  if (!fApplyPatching && (xStatisticalWeight2 || xStatisticalWeightError2 || parAdditional)) {
    printf("Error: Patching is switched off, but input is there! Please enable patching before!\n");
    
    return kFALSE;
  }
  
  fIndexParametersToRegularise = new Int_t[numParams];
  for (Int_t i = 0; i < numParams; i++)
    fIndexParametersToRegularise[i] = indexParametersToRegularise[i];
  
  fLastNotFixedIndexOfParameters = new Int_t[numParamsPerXbin];
  for (Int_t i = 0; i < numParamsPerXbin; i++)
    fLastNotFixedIndexOfParameters[i] = lastNotFixedIndexOfParameters[i];
  
  fXvaluesForRegularisation = new Double_t[fNumXbinsRegularisation];
  fXstatisticalWeight = new Double_t[fNumXbinsRegularisation];
  fXstatisticalWeightError = new Double_t[fNumXbinsRegularisation];
  
  const Int_t numAdditionalParams = numParamsPerXbin * fNumXbinsRegularisation;
  
  if (fApplyPatching) {
    fXstatisticalWeight2 = new Double_t[fNumXbinsRegularisation];
    fXstatisticalWeightError2 = new Double_t[fNumXbinsRegularisation];
    fParAdditional = new Double_t[numAdditionalParams];
  }
  
  Bool_t success = kTRUE;
  
  Bool_t xValuesDecreasing = kFALSE;
  if (fNumXbinsRegularisation >= 2) {
    if (xValuesForRegularisation[1] < xValuesForRegularisation[0]) {
      xValuesDecreasing = kTRUE;
      
      printf("Found decreasing ordering of x values for regularisation!\n");
    }
  }
  
  for (Int_t i = 0; i < fNumXbinsRegularisation; i++) {
    fXvaluesForRegularisation[i] = xValuesForRegularisation[i];
    
    if (i > 0) {
      if ((!xValuesDecreasing && fXvaluesForRegularisation[i] <= fXvaluesForRegularisation[i - 1]) ||
          (xValuesDecreasing  && fXvaluesForRegularisation[i] >= fXvaluesForRegularisation[i - 1])) {
        printf("Error: x values for regularisation must be strictly increasing or decreasing (assumed %s here)!\n",
              xValuesDecreasing ? "decreasing" : "increasing");
        success = kFALSE;
        break;
      }
    }
    
    if (xStatisticalWeight[i] < 0 || xStatisticalWeightError[i] < 0) {
      printf("Error: Invalid statistical weight and/or error for bin %d: value %e, error %e\n",
             i, xStatisticalWeight[i], xStatisticalWeightError[i]);
      success = kFALSE;
      break;
    }
    
    fXstatisticalWeight[i] = xStatisticalWeight[i];
    fXstatisticalWeightError[i] = xStatisticalWeightError[i];
    
    if (fApplyPatching) {
      if (xStatisticalWeight2[i] < 0 || xStatisticalWeightError2[i] < 0) {
        printf("Error: Invalid statistical weight 2 and/or error for bin %d: value %e, error %e\n",
              i, xStatisticalWeight2[i], xStatisticalWeightError2[i]);
        success = kFALSE;
        break;
      }
      
      fXstatisticalWeight2[i] = xStatisticalWeight2[i];
      fXstatisticalWeightError2[i] = xStatisticalWeightError2[i];
    }
  }

  if (fApplyPatching && success) {
    for (Int_t i = 0; i < numAdditionalParams ; i++) {
      if (parAdditional[i] < 0) {
        printf("Error: Invalid parAdditional error for bin %d: %e\n", i, parAdditional[i]);
        success = kFALSE;
        break;
      }
      
      fParAdditional[i] = parAdditional[i];
    }    
  }
  
  if (!success) {
    fNumParametersPerXbin = 0;
    fNumParametersToRegularise = 0;
    
    delete [] fIndexParametersToRegularise;
    fIndexParametersToRegularise = 0x0;
    
    delete [] fLastNotFixedIndexOfParameters;
    fLastNotFixedIndexOfParameters = 0x0;
  
    delete [] fXvaluesForRegularisation;
    fXvaluesForRegularisation = 0x0;
    
    delete [] fXstatisticalWeight;
    fXstatisticalWeight = 0x0;
    
    delete [] fXstatisticalWeightError;
    fXstatisticalWeightError = 0x0;
    
    delete [] fXstatisticalWeight2;
    fXstatisticalWeight2 = 0x0;
    
    delete [] fXstatisticalWeightError2;
    fXstatisticalWeightError2 = 0x0;
    
    delete [] fParAdditional;
    fParAdditional = 0x0;
  
    return kFALSE;
  }
  
  fNumParametersToRegularise = numParams;
  fNumParametersPerXbin = numParamsPerXbin;
  
  return kTRUE;
}



//***Private functions***
//______________________________________________________
AliTPCPIDmathFit::AliTPCPIDmathFit(const Int_t numXbinsRegularisation, const Int_t numSimultaneousFits, const Int_t maxDataPoints): 
fkDoubleEpsilonLimit(2. * std::numeric_limits<double>::min()),
fFunc(0x0), 
fErrorFunc(0x0),
fNpar(0), 
fPar(0x0), 
fErr(0x0), 
fChi2(0), 
fNegLogLikelihood(0),
fNDF(0), 
fNstep(0), 
fNdata(0x0), 
fX(0x0),
fY(0x0),
fXerr(0x0),
fYerr(0x0),
fMaxCalls(2000),
fEpsilon(1.e-4),
fUseLogLikelihood(kFALSE),
fUseWeightsForLoglikelihood(kFALSE),
fUseWeightsForLoglikelihoodForCurrentFit(kFALSE),
fScaleFactorError(1.0),
fRegularisation(0),
fXbinIndex(0),
fNumParametersToRegularise(0),
fNumParametersPerXbin(0),
fIndexParametersToRegularise(0x0),
fLastNotFixedIndexOfParameters(0x0),
fXvaluesForRegularisation(0x0),
fXstatisticalWeight(0x0),
fXstatisticalWeightError(0x0),
fApplyPatching(kFALSE),
fParAdditional(0x0),
fXstatisticalWeight2(0x0),
fXstatisticalWeightError2(0x0),
fRegularisationFactor(1),
fMinimisationString("MINIMIZE"),
fkMaxNdata(maxDataPoints), 
fNumSimultaneousFits(numSimultaneousFits),
fNumXbinsRegularisation(numXbinsRegularisation),
fMinuit(0x0), 
fPolIntX(0x0),
fPolIntY(0x0),
fPolIntC(0x0),
fPolIntD(0x0),
fDebugLevel(0),
fRefHistos(0x0),
fNrefHistos(0)
{
  // Private ctor
  
  if (fNumSimultaneousFits < 1) {
    printf("Error: numSimultaneousFits must be >= 1! Is now set to 1!\n\n");
    fNumSimultaneousFits = 1;
  }
  
  fX = new Double_t**[fNumXbinsRegularisation];
  fY = new Double_t**[fNumXbinsRegularisation];
  fXerr = new Double_t**[fNumXbinsRegularisation];
  fYerr = new Double_t**[fNumXbinsRegularisation];
  
  fNdata = new Int_t*[fNumXbinsRegularisation];
  
  fFunc = new FitFunc_t[fNumSimultaneousFits];
  fErrorFunc = new FitFunc_t[fNumSimultaneousFits];

  for (Int_t i = 0; i < fNumXbinsRegularisation; i++) {
    fX[i] = new Double_t*[fNumSimultaneousFits];
    fY[i] = new Double_t*[fNumSimultaneousFits];
    fXerr[i] = new Double_t*[fNumSimultaneousFits];
    fYerr[i] = new Double_t*[fNumSimultaneousFits];
    
    fNdata[i] = new Int_t[fNumSimultaneousFits];

    for (Int_t j = 0; j < fNumSimultaneousFits; j++) {
      fX[i][j] = new Double_t[fkMaxNdata];
      fY[i][j] = new Double_t[fkMaxNdata];
      fXerr[i][j] = new Double_t[fkMaxNdata];
      fYerr[i][j] = new Double_t[fkMaxNdata];

      for (Int_t k = 0; k < fkMaxNdata; k++) {
        fX[i][j][k] = -999;
        fY[i][j][k] = -999;
        fXerr[i][j][k] = -999;
        fYerr[i][j][k] = -999;
      }
      
      fNdata[i][j] = 0;
    }
  }
  
  for (Int_t i = 0; i < fNumSimultaneousFits; i++) {
    fFunc[i] = 0x0;
    fErrorFunc[i] = 0x0;
  }
}


//______________________________________________________
AliTPCPIDmathFit::~AliTPCPIDmathFit()
{
  // dtor
  
  ClearRefHistos();
  
  for (Int_t i = 0; i < fNumXbinsRegularisation; i++) {
    for (Int_t j = 0; j < fNumSimultaneousFits; j++) {
      delete [] fX[i][j];
      delete [] fY[i][j];
      delete [] fXerr[i][j];
      delete [] fYerr[i][j];
      
      fX[i][j]  = 0x0;
      fY[i][j]  = 0x0;
      fXerr[i][j]  = 0x0;
      fYerr[i][j]  = 0x0;
    }
    
    delete [] fX[i];
    delete [] fY[i];
    delete [] fXerr[i];
    delete [] fYerr[i];
    
    fX[i]  = 0x0;
    fY[i]  = 0x0;
    fXerr[i]  = 0x0;
    fYerr[i]  = 0x0;
  }
        
        
  delete [] fX;
  fX = 0x0;
  
  delete [] fY;
  fY = 0x0;
  
  delete [] fXerr;
  fXerr = 0x0;
  
  delete [] fYerr;
  fYerr = 0x0;
  
  
  for (Int_t i = 0; i < fNumXbinsRegularisation; i++) {
    delete [] fNdata[i];

    fNdata[i]  = 0x0;
  }

  delete [] fNdata;
  fNdata = 0x0;
  
  
  
  delete [] fFunc;
  fFunc = 0x0;
  
  delete [] fErrorFunc;
  fErrorFunc = 0x0;
  
  delete [] fIndexParametersToRegularise;
  fIndexParametersToRegularise = 0x0;
  
  delete [] fLastNotFixedIndexOfParameters;
  fLastNotFixedIndexOfParameters = 0x0;
  
  delete [] fXvaluesForRegularisation;
  fXvaluesForRegularisation = 0;
  
  delete fMinuit;
  fMinuit = 0x0;
  
  delete [] fPar;
  fPar = 0x0;
  
  delete [] fErr;
  fErr = 0x0;
  
  delete [] fXstatisticalWeight;
  fXstatisticalWeight = 0x0;
  
  delete [] fXstatisticalWeightError;
  fXstatisticalWeightError = 0x0;
  
  delete [] fParAdditional;
  fParAdditional = 0x0;
  
  delete [] fXstatisticalWeight2;
  fXstatisticalWeight2 = 0x0;
  
  delete [] fXstatisticalWeightError2;
  fXstatisticalWeightError2 = 0x0;
  
  delete [] fPolIntX;
  fPolIntX = 0x0;
  
  delete [] fPolIntY;
  fPolIntY = 0x0;
  
  delete [] fPolIntC;
  fPolIntC = 0x0;
  
  delete [] fPolIntD;
  fPolIntD = 0x0;
  
  fNumXbinsRegularisation = 0;
  fNumSimultaneousFits = 0;
  
  fgInstance = 0x0;
}


//______________________________________________________
inline Double_t AliTPCPIDmathFit::EvalLog(Double_t x) const
{ 
  // Safely evaluate logarithms -> Adapted root util.h 
  
  if (x < fkDoubleEpsilonLimit)
    return x / fkDoubleEpsilonLimit + TMath::Log(fkDoubleEpsilonLimit) - 1.; 
  else      
    return TMath::Log(x);
}
  

//______________________________________________________
Int_t AliTPCPIDmathFit::Fit(const Double_t *inPar, Double_t *covMatrix, const Double_t *stepSize, const Double_t *lowLim,
                            const Double_t *hiLim)
{
  // Fit the currently loaded data, intialise the parameters with inPar (stepSize, lowLim and hiLim).
  // Save covariance matrix to covMatrix.
  
  Int_t fNdataTotal = 0;
  
  for (Int_t iNumXbinRegularisation = 0; iNumXbinRegularisation < fNumXbinsRegularisation; iNumXbinRegularisation++) {
    for (Int_t iNumSimultaneousFits = 0; iNumSimultaneousFits < fNumSimultaneousFits; iNumSimultaneousFits++) {
      fNdataTotal += fNdata[iNumXbinRegularisation][iNumSimultaneousFits];
      
      for (Int_t idata = 0; idata < fNdata[iNumXbinRegularisation][iNumSimultaneousFits]; idata++) {
        if (fDebugLevel > 1)
          printf("%3d %2d %7d %13.6e %13.6e %10.2e %10.2e\n", iNumXbinRegularisation, iNumSimultaneousFits, 
                 idata, fX[iNumXbinRegularisation][iNumSimultaneousFits][idata],
                 fY[iNumXbinRegularisation][iNumSimultaneousFits][idata], fXerr[iNumXbinRegularisation][iNumSimultaneousFits][idata], 
                 fYerr[iNumXbinRegularisation][iNumSimultaneousFits][idata]);
      }
    }
  }
  
  AliTPCPIDmathFit* ptr = AliTPCPIDmathFit::Instance();
  printf("\n================= Starting MathFit useLogLikelihood %d ndata %d =================\n\n", ptr->GetUseLogLikelihood(),
         fNdataTotal);
  
  fMinuit = new TMinuit(fNpar);
  fMinuit->SetFCN(FitFCN);

  // Standard error is for chi^2 -> Need factor 1/2 for log likelihood
  fMinuit->SetErrorDef(ptr->GetUseLogLikelihood() ? 0.5 : 1);

  // Define the start parameters
  for (Int_t i = 0; i < fNpar; i++)  {
    const Double_t start = inPar ? inPar[i] : 0.5;
    // Here is only the INITIAL step size
    const Double_t ss = stepSize ? stepSize[i] : 0.01;
    const Double_t low = lowLim ? lowLim[i] : 0;
    const Double_t hi = hiLim ? hiLim[i] : 0;

    fMinuit->DefineParameter(i, Form("p%d", i), start, ss, low, hi);
  }
  
  fMinuit->SetMaxIterations(10000);
  
  Double_t arg[2] = {fMaxCalls, fEpsilon};
  Int_t errFlag = 0;
  Int_t hesseFlag = -999;
  
  
  
  // If regularisation is to be used, initialise the function for interpolation.
  // Use a pol(2*n-1) for the interpolation, if +-n bins (set by fRegularisation) are taken into account,
  // i.e. 2 * fRegularisation points are used.
  if (GetUseRegularisation()) {
    if (!fXvaluesForRegularisation) {
      printf("Cannot calculate regularisation: x value array not initialised!\n");
      exit(-1);
    }
    
    // NOTE: The interpolation routine is based on self-defined vectors, which have indices from 1 to n.
    // In order not to mess something up, just add one dummy entry at index zero
    fPolIntX = new Double_t[fRegularisation * 2 + 1];
    fPolIntY = new Double_t[fRegularisation * 2 + 1];
    fPolIntC = new Double_t[fRegularisation * 2 + 1];
    fPolIntD = new Double_t[fRegularisation * 2 + 1];
    
    for (Int_t i = 0; i < fRegularisation * 2 + 1; i++) {
      fPolIntX[i] = 0.;
      fPolIntY[i] = 0.;
      fPolIntC[i] = 0.;
      fPolIntD[i] = 0.;
    }
    
    printf("Regularisation enabled! Interpolating with pol%d between +-%d bins and regularisation factor %f.\n\n",
           2 * fRegularisation - 1, fRegularisation, fRegularisationFactor);
  }
  else 
     printf("Regularisation disabled!\n\n");
  
  // Using "MINIMIZE" means to try "MIGRAD". Only if MIGRAD fails,
  // SIMPLEX and then MIGRAD are applied. This yields to better convergence,
  // but in case of MIGRAD failing the first time, errors are not reliable and
  // also the minimum might not be correct. But this is most likely still better
  // than a failed fit.
  
  if (fUseWeightsForLoglikelihood && fUseLogLikelihood) {
    
    Double_t covMatrixTemp[fNpar][fNpar];
    Double_t hesseMatrixTemp[fNpar][fNpar];
    Double_t matrixTemp[fNpar][fNpar];
    Double_t covMatrixFinal[fNpar][fNpar];
    for (Int_t i = 0; i < fNpar; i++) {
      for (Int_t j = 0; j < fNpar; j++) {
        covMatrixTemp[i][j] = 0;
        hesseMatrixTemp[i][j] = 0;
        matrixTemp[i][j] = 0;
        covMatrixFinal[i][j] = 0;
      }
    }
    
    
    fUseWeightsForLoglikelihoodForCurrentFit = kFALSE;
    
    // Firstly, minimise w/o weighting and obtain covariance matrix
    fMinuit->mnexcm(fMinimisationString.Data(), arg, 2, errFlag);
    GetCovarianceMatrix(&covMatrixTemp[0][0]);
    
    // Save results (to be used if further steps fail)
    for(Int_t ipar = 0; ipar < fNpar; ipar++)
      fMinuit->GetParameter(ipar, fPar[ipar], fErr[ipar]);
    
    if (covMatrix) {
      GetCovarianceMatrix(covMatrix);
      
      
      // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
      for (Int_t i = 0; i < fNpar; ++i) {
        for (Int_t j = 0; j < fNpar; ++j) {
            covMatrix[i * fNpar + j] *= fScaleFactorError * fScaleFactorError; 
        }
      }
    }
    
    // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
    for (Int_t i = 0; i < fNpar; ++i)
      fErr[i] *= fScaleFactorError;
    
    
    // Secondly, enable weithing and obtain HESSE matrix
    fUseWeightsForLoglikelihoodForCurrentFit = kTRUE;
    
    
    // TODO THIS CALL IS NOT NECESSARILY NEEDED fMinuit->mnexcm(fMinimisationString.Data(), arg, 2, errFlag2);
    
    //improve error bars a lot!
    //only 1 argument: fMaxCalls
    fMinuit->mnexcm("HESSE", arg, 1, hesseFlag);
    
    if (hesseFlag == 0) {
      // If HESSE call was successful, save results (to be used if further steps fail)
      for(Int_t ipar = 0; ipar < fNpar; ipar++)
        fMinuit->GetParameter(ipar, fPar[ipar], fErr[ipar]);
      
      if (covMatrix) {
        GetCovarianceMatrix(covMatrix);
      
        // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
        for (Int_t i = 0; i < fNpar; ++i) {
          for (Int_t j = 0; j < fNpar; ++j) {
              covMatrix[i * fNpar + j] *= fScaleFactorError * fScaleFactorError; 
          }
        }
      }
      
      // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
      for (Int_t i = 0; i < fNpar; ++i)
        fErr[i] *= fScaleFactorError;
    }
    
    
    // Hesse matrix is the inverse of the error matrix after the hesse step
    
    Int_t nfree = fMinuit->GetNumFreePars();
    TMatrixDSym orig(nfree);   
    fMinuit->mnemat(orig.GetMatrixArray(), nfree); 
    
    TMatrixDSym inv(orig);
    inv.Invert();
    
    // If the inversion failed, the result is equal to the original
    Bool_t success = kFALSE;
    for (Int_t i = 0; i < nfree; ++i) { 
      for (Int_t j = 0; j < nfree; ++j) {
        if (TMath::Abs(inv(i,j) - orig(i,j)) > 1e-9) {
          success = kTRUE;
          break;
        }
      }
      
      if (success)
        break;
    }
    
    
    // Only re-weight errors, if the matrix inversion and the fits have been successful.
    if ((errFlag + hesseFlag) == 0 && success) {
      Int_t l = 0; 
      for (Int_t i = 0; i < fNpar; ++i) {       
        if (fMinuit->fNiofex[i] > 0) {  // not fixed ?
          Int_t m = 0; 
          for (Int_t j = 0; j <= i; ++j) { 
            if (fMinuit->fNiofex[j] > 0) {  //not fixed
                hesseMatrixTemp[i][j] =  inv(l, m);
                hesseMatrixTemp[j][i] = hesseMatrixTemp[i][j]; 
                m++;
            }
          }
          l++;
        }
      }  
      
      // We want cov * hes * cov
      // -> Do hes * cov first
      for (Int_t i = 0; i < fNpar; ++i) {
        for (Int_t j = 0; j < fNpar; ++j) {
            for (Int_t k = 0; k < fNpar; ++k) 
              matrixTemp[i][j] += hesseMatrixTemp[i][k] * covMatrixTemp[k][j];
        }
      }

      // Now do cov * tmp
      for (Int_t i = 0; i < fNpar; ++i) {
        for (Int_t j = 0; j < fNpar; ++j) {
          for (Int_t k = 0; k < fNpar; ++k) {
            covMatrixFinal[i][j] += covMatrixTemp[i][k] * matrixTemp[k][j];
          }
          
          // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
          covMatrixFinal[i][j] *= fScaleFactorError * fScaleFactorError; 
          
          if (covMatrix)
            covMatrix[i * fNpar + j] = covMatrixFinal[i][j];
        }
      }
      
      // Obtain the parameters (and old errors)
      for(Int_t ipar = 0; ipar < fNpar; ipar++)
        fMinuit->GetParameter(ipar, fPar[ipar], fErr[ipar]);
      
      // Set the new parameter errors
      for(Int_t ipar = 0; ipar < fNpar; ipar++)
        fErr[ipar] = TMath::Sqrt(covMatrixFinal[ipar][ipar]);
    }
    else {
      // An error occurred -> Leave the errors as they are (still better than doing re-weighting with problematic matrices)
      printf("Error: Could not re-weight errors (MIGRAD %d, HESSE %d, matrix inversion %d) -> Errors and covariance matrix not reliable!\n", 
             errFlag, hesseFlag, success);
    }
  }
  else {
    fMinuit->mnexcm(fMinimisationString.Data(), arg, 2, errFlag);

    //improve error bars a lot!
    //only 1 argument: fMaxCalls
    fMinuit->mnexcm("HESSE", arg, 1, hesseFlag);
    
    for(Int_t ipar = 0; ipar < fNpar; ipar++)
      fMinuit->GetParameter(ipar, fPar[ipar], fErr[ipar]);
    
    if (covMatrix) {
      GetCovarianceMatrix(covMatrix);
    
      // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
      for (Int_t i = 0; i < fNpar; ++i) {
        for (Int_t j = 0; j < fNpar; ++j) {
            covMatrix[i * fNpar + j] *= fScaleFactorError * fScaleFactorError; 
        }
      }
    }
    
    // If the parameter errors (NOT(!) the errors for the data points in fY!!!) are to be scaled, we need to do this here!
    for (Int_t i = 0; i < fNpar; ++i)
      fErr[i] *= fScaleFactorError;
  }
  
  fNDF = fNdataTotal - fNpar;
  Double_t fChi2orLogLikelihood = 0;
  
  if (!fUseLogLikelihood) {
    GetChiSquare();
    fChi2orLogLikelihood = fChi2;
  }
  else {
    GetNegLogLikelihood();
    fChi2orLogLikelihood = fNegLogLikelihood;
  }
  
  
  if ((errFlag + hesseFlag) != 0) printf("*****Warning: Abnormal termination: MIGRAD flag %d, HESSE flag %d\n", errFlag, hesseFlag);

  
  printf("======> ChiSquare/NDF: %f / %d%s\nerrFlag %d (MIGRAD %d, HESSE %d)\n", fChi2orLogLikelihood, fNDF, 
         fNDF > 0 ? Form(" = %f", fChi2orLogLikelihood / fNDF) : "", errFlag + hesseFlag, errFlag, hesseFlag);
  errFlag += hesseFlag;

  delete fMinuit;
  fMinuit = 0x0;
  
  delete [] fPolIntX;
  fPolIntX = 0x0;
  
  delete [] fPolIntY;
  fPolIntY = 0x0;
  
  delete [] fPolIntC;
  fPolIntC = 0x0;
  
  delete [] fPolIntD;
  fPolIntD = 0x0;
  
  fXbinIndex = 0;
  
  
  return errFlag;
}


//______________________________________________________
void AliTPCPIDmathFit::FitFCN(Int_t &nPar, Double_t *, Double_t &chiSquareOrLogLikelihood, Double_t *par, Int_t)
{
  // Perform fitting step
  
  AliTPCPIDmathFit* ptr = AliTPCPIDmathFit::Instance();
  ptr->fNstep++;
  
  for (Int_t ii = 0; ii < ptr->fNpar; ii++)
    ptr->fPar[ii] = par[ii];

  if (ptr->GetUseLogLikelihood()) {
    chiSquareOrLogLikelihood = ptr->GetNegLogLikelihood();
  }
  else {
    chiSquareOrLogLikelihood = ptr->GetChiSquare();
  }

  if (ptr->GetDebugLevel() >= 3) {
    printf("%d %f %d --- ", ptr->fNstep, chiSquareOrLogLikelihood, nPar);
    
    for (Int_t i = 0; i < nPar; i++)
      printf("%f ", par[i]);
    
    printf("-----------------\n");
    
    if (ptr->fNstep == 2 && ptr->GetDebugLevel() >= 4)
       exit(1);
  }
}


//______________________________________________________
Double_t AliTPCPIDmathFit::GetChiSquare()
{
  // Calculates chi^2 and sets this value for the internal variable fChi2. In addition, the calculated value is returned.
  
  fChi2 = 0;
  
  for (Int_t xBin = 0; xBin < fNumXbinsRegularisation; xBin++) {
    // This parameter tells the following fFunc call which x bin (and which parameters) are to be used for the fit
    // of the current x bin
      
    fXbinIndex = xBin; 
    
    for (Int_t iSimultaneousFit = 0; iSimultaneousFit < fNumSimultaneousFits; iSimultaneousFit++) {
      for(Int_t idata = 0; idata < fNdata[xBin][iSimultaneousFit]; idata++) {
        const Double_t evaly = fFunc[iSimultaneousFit](&(fX[xBin][iSimultaneousFit][idata]), fPar);
        const Double_t dy = fY[xBin][iSimultaneousFit][idata] - evaly;
        
        Double_t dfdx = 0;
        //calculated df/dx if necessary
        if (fXerr[xBin][iSimultaneousFit][idata] != 0)  {
          const Double_t derix0 = fX[xBin][iSimultaneousFit][idata] - fXerr[xBin][iSimultaneousFit][idata];
          const Double_t derix1 = fX[xBin][iSimultaneousFit][idata] + fXerr[xBin][iSimultaneousFit][idata];
          const Double_t deriy0 = fFunc[iSimultaneousFit](&derix0, fPar);
          const Double_t deriy1 = fFunc[iSimultaneousFit](&derix1, fPar);
          if (abs(derix1 - derix0) < fkDoubleEpsilonLimit) {
            AliError("derix1-derix0 is 0 -> Devision by zero!");
            dfdx = 0;
          }
          else  
            dfdx = (deriy1 - deriy0) / (derix1 - derix0);
        }
        Double_t den = TMath::Power(fYerr[xBin][iSimultaneousFit][idata], 2) +
                       TMath::Power(dfdx * fXerr[xBin][iSimultaneousFit][idata], 2);
        
        if (fErrorFunc[iSimultaneousFit])  {
          den += fErrorFunc[iSimultaneousFit](&(fX[xBin][iSimultaneousFit][idata]), fPar);
        }
        
        Double_t part = 0;
        
        if (abs(den) < fkDoubleEpsilonLimit) {
          AliError("den is 0 -> Devision by zero!");
          part = 0;
        }
        else 
          part = dy * dy / den;
        
        if (fDebugLevel >= 4)  {
          for(Int_t ip = 0; ip < fNpar; ip++) {
            printf("par%d %e\n", ip, fPar[ip]);
          }
          printf("%d %d: %e %e %e -- %e %e -- %e %e %e -- %e\n", iSimultaneousFit, idata, fX[xBin][iSimultaneousFit][idata], 
                 fY[xBin][iSimultaneousFit][idata], fYerr[xBin][iSimultaneousFit][idata], evaly, dy, dfdx, den, part, fChi2);
        }
        
        fChi2 += part;
      }
    }
  }
  
  AddPenaltyTermForRegularisation(kFALSE);
  
  return fChi2;
}


//______________________________________________________
Double_t AliTPCPIDmathFit::GetNegLogLikelihood()
{
  // Implementation based on ROOT's FitUtil::EvaluatePoissonLogL in $ROOTSYS/math/mathcore/src/FitUtil.cxx.
  // Especially the weighted branch is based on the root version.
  
  fNegLogLikelihood = 0;
  
  if (fUseWeightsForLoglikelihoodForCurrentFit) {
    for (Int_t xBin = 0; xBin < fNumXbinsRegularisation; xBin++) {
      // This parameter tells the following fFunc call which x bin (and which parameters) are to be used for the fit
      // of the current x bin
      
      fXbinIndex = xBin; 
      
      for (Int_t iSimultaneousFit = 0; iSimultaneousFit < fNumSimultaneousFits; iSimultaneousFit++) {
        for(Int_t idata = 0; idata < fNdata[xBin][iSimultaneousFit]; idata++) {
          
          // Need to apply weight correction. Effective weight is error^2 / y
          // and expected events in bins is evalY/weight.
          // One can apply the correction only when y is not zero otherwise weight is undefined
          // (in case of weighted likelihood one doesn't care about the constant term due to 
          // the saturated model).
          
          if (TMath::Abs(fY[xBin][iSimultaneousFit][idata]) > fkDoubleEpsilonLimit) {
            Double_t evalY = fFunc[iSimultaneousFit](&(fX[xBin][iSimultaneousFit][idata]), fPar);
            evalY = TMath::Max(evalY, 0.0);  // Avoid negative or too small values 
            
            // Bin effective weight
            const Double_t weight = (fYerr[xBin][iSimultaneousFit][idata] * fYerr[xBin][iSimultaneousFit][idata]) /
                                    fY[xBin][iSimultaneousFit][idata];  
            fNegLogLikelihood += weight * evalY - weight * fY[xBin][iSimultaneousFit][idata] * EvalLog(evalY);
          }
        }
      }
    }
  }
  else {
    // Binned Poissonian likelihood chi-square as in http://www.scribd.com/doc/65048561/1983-S-Baker-R-D-Cousins-NIM-221-437-442
    // -> Download http://www.sciencedirect.com/science/article/pii/0167508784900164#
    for (Int_t xBin = 0; xBin < fNumXbinsRegularisation; xBin++) {
      // This parameter tells the following fFunc call which x bin (and which parameters) are to be used for the fit
      // of the current x bin
      fXbinIndex = xBin; 
      
      for (Int_t iSimultaneousFit = 0; iSimultaneousFit < fNumSimultaneousFits; iSimultaneousFit++) {
        for(Int_t idata = 0; idata < fNdata[xBin][iSimultaneousFit]; idata++) {
          Double_t evalY = fFunc[iSimultaneousFit](&(fX[xBin][iSimultaneousFit][idata]), fPar);
          // Since the fit MINIMIZES, we need the negative value of the log likelihood here,
          // i.e. minus sign already taken into account in the following!
          evalY = TMath::Max(evalY, 0.0);  // Avoid negative or too small values 

          Double_t part = evalY - fY[xBin][iSimultaneousFit][idata];
          if (fY[xBin][iSimultaneousFit][idata] > 0.0) {
            part += fY[xBin][iSimultaneousFit][idata] * (EvalLog(fY[xBin][iSimultaneousFit][idata]) - EvalLog(evalY));  
          }
          
          if(fDebugLevel >= 4) 
            printf("%d %d: %e -- %e %e %e\n", iSimultaneousFit, idata, fX[xBin][iSimultaneousFit][idata], evalY, part,
                   fNegLogLikelihood);
          
          // NO factor 2 as in paper. This is because the error definition is set to 0.5, which exactly absorbs the parameter
          // 2 here (error means that changing the parameters by the error, the negLogLikelihood changes by the error definition)
          // -> So changing the negLogLikelihood as defined here by 1 sigma produces 0.5 of the change of the real negLoglikelihood
          // with the factor 2 in front.
          fNegLogLikelihood += part; 
        }
      }
    }
  }
  
  AddPenaltyTermForRegularisation(kTRUE);
  
  return fNegLogLikelihood;
}


//______________________________________________________
void AliTPCPIDmathFit::AddPenaltyTermForRegularisation(const Bool_t useLogLikelihood)
{
  // In case of regularisation, add the corresponding penalty term to chi^2 or negLogLikelihood.
  
  if (GetUseRegularisation()) {
    // In the following, the roles of x and y are different from the actual fit!
    // In particular, y is now the value of a certain parameter!
    
    // Loop over all x bins to obtain the penalty term
    for (Int_t iPar = 0; iPar < fNumParametersToRegularise; iPar++) {
      const Int_t currParIndex = fIndexParametersToRegularise[iPar];
      
      // currParIndex < 0 means that this parameter in this xBin should NOT contribute to the penalty term
      if (currParIndex < 0) {
        if (fDebugLevel > 1)
          printf("Skipping for reg: %d - parInBin %d, xBin %d\n", -currParIndex, -currParIndex % fNumParametersPerXbin, (Int_t)(-currParIndex / fNumParametersPerXbin));
        continue;
      }
      
      // the xBin of the current parameter can be obtained from the parameter index
      // and the number of parameters per xBin
      const Int_t xBin = (Int_t)(currParIndex / fNumParametersPerXbin);
      
      const Double_t xStatisticalWeightSummed = fApplyPatching ? (fXstatisticalWeight[xBin] + fXstatisticalWeight2[xBin])
                                                               :  fXstatisticalWeight[xBin];
      const Double_t xStatisticalWeightErrorSquaredSummed = fApplyPatching
                                                                ? fXstatisticalWeightError[xBin]*fXstatisticalWeightError[xBin] +  
                                                                  fXstatisticalWeightError2[xBin]*fXstatisticalWeightError2[xBin]
                                                                :  fXstatisticalWeightError[xBin]*fXstatisticalWeightError[xBin];
      if (xStatisticalWeightSummed < fkDoubleEpsilonLimit)
          continue; // If there is no statistics in this bin, one can safely skip it (the term cannot be calculated anyway)
      
      // Determine the number of points used for the regularisation (= order of interpolation polynomial)
      const Int_t nRegBins = 2 * fRegularisation; // Regularise with +-fRegularisation bins
      
      // No regularisation at the edges, i.e. less than fRegularisation bins left or right of the current bin
      // => This is IMPORTANT. Consider for example fRegularisation = 1, i.e. +-1 bins. If the left bin is missing,
      // the "inter"polation polynomial would be a pol0, i.e. a constant. But for low pT there is some finite slope
      // at the edges. So, extrapolating to the edges would bias the results.
      // The edge points are constraint by their neighbours in the sense that they provide an interpolation value.
      // A bin is also an "edge bin", if the right neighbour bin is fixed
      if (xBin - fRegularisation < 0) 
        continue; // Points left of the current bin are missing
      if (xBin + fRegularisation > fNumXbinsRegularisation - 1)
        continue; // Points right of the current bin are missing
      if (xBin + fRegularisation > fLastNotFixedIndexOfParameters[currParIndex % fNumParametersPerXbin])
        continue; // Points right of the current bin cannot be used, since they are fixed
      
      /* OLD
      Int_t nRegBins = 2 * fRegularisation; // Regularise with +-fRegularisation bins
      if (xBin - fRegularisation < 0) 
        nRegBins -= fRegularisation - xBin; // Points left of the current bin are missing
      if (xBin + fRegularisation > fNumXbinsRegularisation - 1)
        nRegBins -= (xBin + fRegularisation) - (fNumXbinsRegularisation - 1); // Points right of the current bin are missing
      
      // To perform the interpolation, the interpolation function must be well determined.
      // This is only true for nRegBins >= DOF(pol(2 * fRegularisation - 1)) = 2 * fRegularisation.
      // Therefore, bins at the edges will not be regularised.
        
      if (nRegBins < 2 * fRegularisation)
        continue;
      
      
      
      Int_t xBinTemp = TMath::Max(0, xBin - fRegularisation);
      */
      
      // Bin effective weight required for weighted log likelihood
      const Double_t weight = fUseWeightsForLoglikelihoodForCurrentFit 
                                ? (xStatisticalWeightErrorSquaredSummed / xStatisticalWeightSummed)
                                : 1.0;
      
      // Current values of this parameter
      const Double_t currYvalue = fApplyPatching ? PatchParameter(fPar[currParIndex], fXstatisticalWeight[xBin], fParAdditional[currParIndex],
                                                                  xStatisticalWeightSummed)
                                                 : fPar[currParIndex];
      
      // Set the data points for the interpolation
      
      // Reset first; i = 0 is not used!
      for (Int_t i = 1; i < fRegularisation * 2 + 1; i++) {
        fPolIntX[i] = 0.;
        fPolIntY[i] = 0.;
      }

      Int_t xBinTemp = xBin - fRegularisation; // Above if statement already guarantees xBinTemp >= 0
      for (Int_t i = 0; i < nRegBins; xBinTemp++) {
        if (xBinTemp == xBin) // Do not add considered bin to interpolation
          continue;
        
        const Int_t currParIndexTemp = currParIndex % fNumParametersPerXbin + xBinTemp * fNumParametersPerXbin;
        // NOTE: +1 to ignore dummy index 0!
        fPolIntX[i + 1] = fXvaluesForRegularisation[xBinTemp];
        fPolIntY[i + 1] = fApplyPatching ? PatchParameter(fPar[currParIndexTemp], fXstatisticalWeight[xBinTemp], 
                                                          fParAdditional[currParIndexTemp],
                                                          fXstatisticalWeight[xBinTemp] + fXstatisticalWeight2[xBinTemp])
                                         : fPar[currParIndexTemp];
        
        i++;
      }
      
      // Get the interpolated y value
      Double_t yInterpolated = 0., yInterPolatedErr = 0;
      if (!PolynomialInterpolation(fPolIntX, fPolIntY, nRegBins, fXvaluesForRegularisation[xBin], &yInterpolated, &yInterPolatedErr)) 
      {
        printf("Polynomial interpolation failed! Regularisation term for current x bin (%d) ignored\n", xBin);
        continue;
      }
      
      // Calculate the penalty term...
      const Double_t yDiff = currYvalue - yInterpolated;
      
      // For the error, just take the statistical error (assuming the statistical weight for this bin to have no error, since
      // this is only a number of measured counts, which is a fact). A weighting correction will be applied for the
      // penalty term.
      
      const Double_t yError2 = (0.5 * (currYvalue + yInterpolated)) / xStatisticalWeightSummed;
      
      Double_t regPenaltyTerm = 0.;
      if (yError2 > fkDoubleEpsilonLimit)
        regPenaltyTerm = fRegularisationFactor * weight * (yDiff * yDiff) / yError2;
      //else {
        // Should only happen, if currYvalue ~ yInterpolates ~ 0. But then, yDiff is ~ 0 also. So, take regPenaltyTerm = 0 in this case.
      //}
        
      // Each fit of the simultaneous fit group is a fit of its own and has its own regularisation. Thus,
      // the penalty term needs to be added for every individual fit in order to get the same strength for
      // the reg. as for non-simultaneous fit. This has been checked and yields consistent results with only
      // one single fit.
      regPenaltyTerm *= fNumSimultaneousFits;
      
      // ...and add it to the negLogLikelihood or the chi^2 depending on the fit method.
      // The factor 0.5 for the log likelihood is to take the proper definition of negLogLikelihood
      // (which is 0.5), i.e. changing by yError should change the penalty term by 0.5.
      if (useLogLikelihood)
        fNegLogLikelihood += 0.5 * regPenaltyTerm;
      else
        fChi2 += regPenaltyTerm;
    }
  }
}


//______________________________________________________
void AliTPCPIDmathFit::GetCovarianceMatrix(Double_t* covMatrix) const
{
  // Retrieve covariance matrix.
  
  if (!fMinuit || !covMatrix)
    return;

  // Take special care in case of fixed parameters, i.e. fill in zeros, compare
  // http://root.cern.ch/root/html/src/TMinuitMinimizer.cxx.html#Ilx.Z

  const Int_t nfree = fMinuit->GetNumFreePars();

  TMatrixDSym matfree(nfree); 
  fMinuit->mnemat(matfree.GetMatrixArray(), nfree);

  // r: row, c: column
  Int_t rfree = 0;
  for (Int_t rr = 0; rr < fNpar; rr++) {
    if (fMinuit->fNiofex[rr] <= 0) // fMinuit->fNiofex[rr] is zero in case of fix parameter
      continue;

    Int_t cfree = 0;
    for (Int_t cc = 0; cc < fNpar; cc++) {
      if (fMinuit->fNiofex[cc] <= 0) // fMinuit->fNiofex[cc] is zero in case of fix parameter
        continue;
      
      covMatrix[rr * fNpar + cc] = matfree[rfree][cfree];
      cfree++;
    }
    rfree++;
  }
}


//______________________________________________________
inline Double_t AliTPCPIDmathFit::PatchParameter(Double_t par, Double_t statWeightTot1, Double_t statWeightPar2, Double_t statWeightTotSummed)
{
  // Patch parameter - NOTE: statWeightTotSummed is assumed to be statWeightTot1 + statWeightTot2!!!
  // E.g. TOF patching: par = fraction_species, statWeightTot1 = TPC_yield_tot, statWeightTot2 = TOF_yield_tot,
  // statWeightPar2 = TOF_yield_species
  if (statWeightTotSummed > 0.)
    return (par * statWeightTot1 + statWeightPar2) / statWeightTotSummed;
  
  return 0.;
}


//______________________________________________________
Bool_t AliTPCPIDmathFit::PolynomialInterpolation(Double_t xa[], Double_t ya[], Int_t n, Double_t x, Double_t* y, Double_t* dy)
{
  // Code from "Numerical Recipes in C" with slight modifications.
  // In the following, n must be <= 2*fRegularisation!
  // Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y and
  // an error estimate dy. If P(x) is the polynomial of degree N − 1 such that P(xai) = yai, i =
  // 1, . . . , n, then the returned value y = P(x).
  
  
  Int_t i, m, ns = 1;
  Double_t den, dif, dift, ho, hp, w;
  
  dif = fabs(x - xa[1]);
  
  for (Int_t j = 1; j <= n; j++) {
    fPolIntC[j] = 0.;
    fPolIntD[j] = 0.;
  }
  
  for (i = 1; i <= n; i++) {
    // Find the index ns of the closest table entry
    if ( (dift = fabs(x - xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    
    // Initialise tableau of c's and d's
    fPolIntC[i] = ya[i];
    fPolIntD[i] = ya[i];
  }
  
  *y = ya[ns--]; // This is the initial approximation to y.
  for (m = 1; m < n; m++) { // For each column of the tableau,
    for (i = 1; i <= n - m; i++) { // loop over the current c’s and d’s and update them
      ho = xa[i] - x;
      hp = xa[i + m] - x;
      w = fPolIntC[i + 1] - fPolIntD[i];
      
      if ( fabs((den = ho - hp)) < fkDoubleEpsilonLimit) {
        // This error can occur only if two input xa’s are (to within roundoff) identical.
        printf("Error in PolynomialInterpolation!\n");
        
        return kFALSE;
      }
      
      den = w / den;
      
      // Update c's and d's
      fPolIntD[i] = hp * den;
      fPolIntC[i] = ho * den;
    }
    
    // After each column in the tableau is completed, we decide which correction, c or d,
    // we want to add to our accumulating value of y, i.e., which path to take through the
    // tableau—forking up or down. We do this in such a way as to take the most “straight
    // line” route through the tableau to its apex, updating ns accordingly to keep track of
    // where we are. This route keeps the partial approximations centered (insofar as possible)
    // on the target x. The last dy added is thus the error indication.
    *y += (*dy = (2 * ns < (n-m) ? fPolIntC[ns + 1] : fPolIntD[ns--]));

  }
  
  return kTRUE;
}
