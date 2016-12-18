#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
//#include "Math/RootFinderAlgorithms.h"


//________________________________________________________
Double_t getWeightedMean(Double_t* x, Double_t* par)
{
  // Estimate the weighted mean with the chi^2 method from:
  // J. Phys A: Math. Gen. 25 (1992) 1967-1979
  
  const Double_t sysError = x[0];
  const Int_t nPoints = (Int_t)par[1];
  
  const Double_t* xMean = &par[2];
  const Double_t* xSigma = &par[2 + nPoints];
  
  const Double_t sysError2 = sysError * sysError;
  
  // Calculate the weighted mean x^
  Double_t errTot2Inverse = 0;
  Double_t weightedMean = 0;
  
  Double_t totErr2OnePoint = 0;
  for (Int_t i = 0; i < nPoints; i++) {
    totErr2OnePoint = xSigma[i] * xSigma[i] + sysError2;
    errTot2Inverse += totErr2OnePoint > 0 ? 1. / totErr2OnePoint : 1e12;
    
    weightedMean += totErr2OnePoint > 0 ? xMean[i] / totErr2OnePoint : xMean[i] * 1e12;
  }
  
  if (errTot2Inverse > 0)
    weightedMean /= errTot2Inverse;
  else
    weightedMean *= 1e12;
  
  return weightedMean;
}


//________________________________________________________
Double_t chiSquareEstimateOfSysError(Double_t* x, Double_t* par)
{
  // Estimate the systematic error with the chi^2 method from:
  // J. Phys A: Math. Gen. 25 (1992) 1967-1979
  
  const Double_t sysError = x[0];
  const Double_t chi2Expected = par[0];
  const Int_t nPoints = (Int_t)par[1];
  const Bool_t ignoreSigmaErrors = (Bool_t)par[2 * nPoints + 2];
  
  const Double_t* xMean = &par[2];
  const Double_t* xSigma = &par[2 + nPoints];
  
  const Double_t sysError2 = sysError * sysError;
  
  // Calculate the weighted mean x^
  Double_t weightedMean = getWeightedMean(x, par);
  
  Double_t totErr2OnePoint = 0;
  Double_t chi2 = 0;
  for (Int_t i = 0; i < nPoints; i++) {
    totErr2OnePoint = xSigma[i] * xSigma[i] + sysError2;
    
    chi2 += totErr2OnePoint > 0 ? (xMean[i] - weightedMean) * (xMean[i] - weightedMean) / totErr2OnePoint
                                : (xMean[i] - weightedMean) * (xMean[i] - weightedMean) * 1e12;
    if (!(totErr2OnePoint > 0) && !ignoreSigmaErrors && xMean[i] > 0) {
      printf("Warning: xSgima2[%d] = %f, but xMean[%d] = %f\n", i, xSigma[i], i, xMean[i]);
    }
    
    //printf("%d: xMean %f, xSigma %f, weightedMean %f, errTot2Inverse %f, sysError %f, chi2 - (N - 1) %e\n", 
    //       i, xMean[i], xSigma[i], weightedMean, errTot2Inverse, sysError, chi2 - (nPoints - 1));
  }
  
  return chi2 - chi2Expected;
}


//________________________________________________________
Double_t findSystematicError(Int_t nPoints, Double_t* xMean, Double_t* xSigma, Bool_t ignoreSigmaErrors,
                             Double_t minError = 0.0, Double_t maxError = 1.0, Double_t minMean = 0., Double_t maxMean = 1.,
                             Double_t* xMeanResult = 0x0) 
{
  Bool_t success = kTRUE;
  Double_t sysError[3] = {999., 999., 999.};
  
  // Root of function is estimate for systematic error
  // Means lie all within [minMean, maxMean]
  TF1* f = new TF1("findSystematicError", chiSquareEstimateOfSysError, minMean, maxMean, 2 * nPoints + 2 + 1);
  f->SetParameter(2 * nPoints + 2, ignoreSigmaErrors);
  f->SetParameter(1, nPoints);
  
  for (Int_t i = 0; i < nPoints; i++) {
    f->SetParameter(i + 2, xMean[i]);
    f->SetParameter(i + (2 + nPoints), ignoreSigmaErrors ? 0. : xSigma[i]);
  }
    
  for (Int_t iter = 0; iter < 3; iter++) {
    // Mean expected chi^2
    Double_t chi2Percentile = 0.5;
    
    // Mean expected chi^2 not exactly known - allow for +- 1 sigma change
    if (iter == 1)
      chi2Percentile = 0.16;
    else if (iter == 2)
      chi2Percentile = 0.84;
    
    
    f->SetParameter(0, nPoints - 1);
    //OLD: more accurate, but different from common definition of sys. error 
    //f->SetParameter(0, TMath::ChisquareQuantile(chi2Percentile, nPoints - 1)); // NDF = nPoints - 1
    
    
    // Wrap function and use Brent's method to find the root
    ROOT::Math::WrappedTF1 wf(*f);
    ROOT::Math::BrentRootFinder brf;

    // In case of particle fractions, error must be within [0, 1], in case of yields the
    // error should be: 0 <= error <= yield
    brf.SetFunction(wf, minError, maxError);
    brf.SetNpx(1000);
    brf.SetDefaultNSearch(30);
    Bool_t currSuccess = brf.Solve(1000);
    
    if (iter == 0 && xMeanResult) {
      // Means lie all within [minMean, maxMean]
      TF1* fMean = new TF1("findMean", getWeightedMean, minMean, maxMean, 2 * nPoints + 2 + 1);
      fMean->SetParameters(f->GetParameters());
      *xMeanResult = fMean->Eval(brf.Root());
      delete fMean;
    }
    
    success = success && currSuccess;
    
    if (currSuccess)
      sysError[iter] = TMath::Max(0., brf.Root());
    
    break;//TODO iter > 0 not used
  }

  delete f;
  
  // TODO OLD As the final estimate of the systematic error, take the error for the expected chi^2 (iter 0)
  // and add the sigma of the 3 estimates (i.e. chi^2 varied by +-1sigma around mean)
  //TODO Double_t sigmaSysErrorEstimates = TMath::RMS(3, sysError);
  
  Double_t sysErrorEstimateFinal = sysError[0];//TMath::Sqrt(sysError[0] * sysError[0] + sigmaSysErrorEstimates * sigmaSysErrorEstimates); 
  
  if (!success) {
    printf("Error: Failed to find root!\n");
    sysErrorEstimateFinal = 999.;
  }
  
  /*
  if (sysErrorEstimateFinal > 900) {
    printf("%e (nPoints %d)\nminMean %e, maxMean %e, minError %e, maxError %e\n", sysErrorEstimateFinal, nPoints, minMean, maxMean, minError, maxError);
    
    for (Int_t i = 0; i < nPoints; i++)
      printf("mean %d: %e\n", i, xMean[i]);
    printf("\n\n");
  }*/
  

  return TMath::Max(0., sysErrorEstimateFinal);
}


//________________________________________________________
Bool_t extractSystematicError(const Int_t nHistos, TH2D** hInput, TH2D* hResults, const Bool_t setMean,
                              const Bool_t addErrorsQuadratically, const Bool_t ignoreSigmaErrors) 
{
  if (!hInput || !hResults)
    return kFALSE;
  
  Double_t meansForFit[nHistos];
  Double_t sigmasForFit[nHistos];
  
  const Int_t nBinsX = hResults->GetNbinsX();
  const Int_t nBinsY = hResults->GetNbinsY();
  
  for (Int_t binY = 1; binY <= nBinsY; binY++) {
    for (Int_t binX = 1; binX <= nBinsX; binX++) {
      Double_t maxMean = -1.;
      Double_t minMean = -1.;
      Bool_t minMaxSet = kFALSE;
      
      for (Int_t j = 0; j < nHistos; j++) {
        meansForFit[j] = hInput[j]->GetBinContent(binX, binY);
        sigmasForFit[j] = hInput[j]->GetBinError(binX, binY);
        
        if (!minMaxSet || meansForFit[j] > maxMean)
          maxMean = meansForFit[j];
        
        if (!minMaxSet || meansForFit[j] < minMean)
          minMean = meansForFit[j];
        
        minMaxSet = kTRUE;
      }
      
      /* Calculate sigma directly
      Double_t weightedMean = 0.;
      for (Int_t i = 0; i < nHistos; i++) 
        weightedMean += meansForFit[i];

      weightedMean /= nHistos;

      Double_t sysError = 0;
      for (Int_t i = 0; i < nHistos; i++) {
        sysError += (meansForFit[i] - weightedMean) * (meansForFit[i] - weightedMean);
      }

      sysError /= TMath::ChisquareQuantile(0.5, nHistos - 1);//nHistos - 1;
      
      sysError = TMath::Sqrt(sysError);
      */
      
      Double_t weightedMean = 0;
      // The root finder needs some given range. For yields, the error should lie within [0, maxMean - minMean]
      Double_t sysError = findSystematicError(nHistos, meansForFit, sigmasForFit, ignoreSigmaErrors, 0, maxMean - minMean, minMean, maxMean,
                                              &weightedMean);
      
      if (setMean)
        hResults->SetBinContent(binX, binY, weightedMean);
      
      if (addErrorsQuadratically) {
        Double_t currBinError = hResults->GetBinError(binX, binY);
        hResults->SetBinError(binX, binY, TMath::Sqrt(currBinError * currBinError + sysError * sysError));
      }
      else
        hResults->SetBinError(binX, binY, sysError);
      
      /*
      TH1D* h = new TH1D(Form("h_%s", hResults->GetName()), "", TMath::Max(TMath::Ceil(maxMean), 1000.), 0, maxMean);
      for (Int_t j = 0; j < nHistos; j++)
        h->Fill(meansForFit[j]);
        
      new TCanvas();
      h->Draw();
      printf("%s: RMS %e, sysError %e, RMS/mean %e, sysError/mean %e\n", h->GetName(), h->GetRMS(), sysError, h->GetRMS()/weightedMean,
              sysError/weightedMean);
      */
    }
  }
  
  return kTRUE;
}
