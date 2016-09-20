#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TSystem.h"

#include "AliPID.h"

#include <iostream>

#include "THnSparseDefinitions.h"

//#include "/hera/alice/bhess/trunk/bhess_PID/AliAnalysisTaskPID.h"
#include "FFs/bhess_PID/AliAnalysisTaskPID.h"
#include "SystematicErrorUtils.h"

//enum axisDataProj { kPidPtProj = 0, kPidJetPtProj = 1, kPidZProj = 2, kPidXiProj = 3, kNDimsProj = 4 };
enum axisDataProj { kPidPtProj = 0, kPidZProj = 1, kPidXiProj = 2, kNDimsProj = 3 };

const TString yieldAxisTitlePt = "1/N_{Jets} dN/dp_{T} (GeV/c)^{-1}";
const TString yieldAxisTitleZ  = "1/N_{Jets} dN/dz";
const TString yieldAxisTitleXi = "1/N_{Jets} dN/d#xi";


//___________________________________________________________________
void setupHist(TH1* h, TString histName, TString histTitle, TString xAxisTitle, TString yAxisTitle, Int_t color, Bool_t isMC)
{
  if (!h)
    return;
  
  if (histName != "")
    h->SetName(histName.Data());
  h->SetTitle(histTitle.Data());
  
  if (xAxisTitle != "")
    h->GetXaxis()->SetTitle(xAxisTitle.Data());
  if (yAxisTitle != "")
    h->GetYaxis()->SetTitle(yAxisTitle.Data());
  
  h->SetMarkerStyle(isMC ? 24 : 20);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  
  h->SetStats(kFALSE);
}


//____________________________________________________________________________________________________________________
void normaliseYieldHist2D(TH2* hData, TH2* hNumJets, const Int_t lowerCentralityBinLimit, const Int_t upperCentralityBinLimit)
{
  // Normalise to 1/numJets and bin width. NOTE: jetPt binning of hData and hNumJets assumed to be the same!
  
  if (!hData)
    return;
  
  for (Int_t binJetPt = 0; binJetPt <= hData->GetNbinsY() + 1; binJetPt++) {
    // No normalisation to numJets, if no histo is provided
    const Double_t numJets = hNumJets ? hNumJets->Integral(lowerCentralityBinLimit, upperCentralityBinLimit, binJetPt, binJetPt) : 1.;
    Bool_t noJets = numJets < 1e-13;
    
    for (Int_t binObs = 0; binObs <= hData->GetNbinsX() + 1; binObs++) {
      if (noJets) {
        if (hData->GetBinContent(binObs, binJetPt) > 0.) {
          printf("Error: No jets for jetPt ~ %f, but found content %f at y-coord %f!\n", hData->GetYaxis()->GetBinCenter(binJetPt),
                 hData->GetBinContent(binObs, binJetPt),  hData->GetXaxis()->GetBinCenter(binObs));
        }
        continue;
      }
      const Double_t dObservable = hData->GetXaxis()->GetBinWidth(binObs);
      const Double_t normFactor = 1. / (numJets * dObservable);
      
      hData->SetBinContent(binObs, binJetPt, hData->GetBinContent(binObs, binJetPt) * normFactor);
      hData->SetBinError(binObs, binJetPt, hData->GetBinError(binObs, binJetPt) * normFactor);
    }
  }
}

//____________________________________________________________________________________________________________________
void normaliseYieldHist(TH1* h, Double_t numJets)
{
  // Normalise to 1/numJets and bin width
  
  if (!h)
    return;
  
  if (numJets <= 0) // Do not normalise
    numJets = 1.;
  
  for (Int_t bin = 0; bin <= h->GetNbinsX() + 1; bin++) {
    const Double_t dObservable = h->GetBinWidth(bin);
    const Double_t normFactor = 1. / (numJets * dObservable);
    h->SetBinContent(bin, h->GetBinContent(bin) * normFactor);
    h->SetBinError(bin, h->GetBinError(bin) * normFactor);
  }
}


//___________________________________________________________________
void normaliseWeightingFactor(TH2* hWeightingFactor, TH1* hTotalWeight)
{
  // Normalise the weighting factors in hWeightingFactor such that every row sums up to 1.
  // In this case, hTotalWeight is assumed to contain just the projection on the y-axis of
  // hWeightingFactor (especially this implies the same bins in y).
  
  if (!hWeightingFactor || !hTotalWeight) {
    printf("Error normaliseWeightingFactor: Missing histos!\n");
    return;
  }
  
  for (Int_t binY = 1; binY <= hWeightingFactor->GetNbinsY(); binY++) {
    const Double_t totalWeight = hTotalWeight->GetBinContent(binY);
    if (totalWeight > 0.) {
      const Double_t normFactor = 1. / totalWeight;
      
      for (Int_t binPt = 1; binPt <= hWeightingFactor->GetNbinsX(); binPt++) {
        hWeightingFactor->SetBinContent(binPt, binY, hWeightingFactor->GetBinContent(binPt, binY) * normFactor);
        hWeightingFactor->SetBinError(binPt, binY, hWeightingFactor->GetBinError(binPt, binY) * normFactor);
      }
    }
  }
}


//___________________________________________________________________
void normaliseErrorWeightingFactor(TH2* hErrorWeightingFactor, TH1* hTotalErrorWeight)
{
  // Normalise the weighting factors in hErrorWeightingFactor such that every column sums up to 1.
  // In this case, hTotalErrorWeight is assumed to contain just the projection on the x-axis of
  // hErrorWeightingFactor (especially this implies the same bins in pT).
  
  if (!hErrorWeightingFactor || !hTotalErrorWeight) {
    printf("Error normaliseErrorWeightingFactor: Missing histos!\n");
    return;
  }
  
  for (Int_t binPt = 1; binPt <= hErrorWeightingFactor->GetNbinsX(); binPt++) {
    const Double_t totalWeight = hTotalErrorWeight->GetBinContent(binPt);
    if (totalWeight > 0.) {
      const Double_t normFactor = 1. / totalWeight;
      
      for (Int_t binY = 1; binY <= hErrorWeightingFactor->GetNbinsY(); binY++) {
        hErrorWeightingFactor->SetBinContent(binPt, binY, hErrorWeightingFactor->GetBinContent(binPt, binY) * normFactor);
        hErrorWeightingFactor->SetBinError(binPt, binY, hErrorWeightingFactor->GetBinError(binPt, binY) * normFactor);
      }
    }
  }
}

/*TODO Try with iterative approach (seems to fail or convergence is questionable!)
//___________________________________________________________________
void calculateWeightedMean(TH1D** hFractionIDFFvar, TH1D** hFractionIDFFtrackPt, TH2* hWeightingVarvsPt, TH2* hErrorWeightingVarvsPt,
                           Bool_t zerothIter, Double_t alpha)
{
  // Loop over all var (=z, xi) bins and calculate for the mean fraction: Fraction(var) = Sum_Pt (fraction(Pt) * weight(Pt, var)).
  // For the error, calculate:
  // FractionError(var)^2 = Sum_Pt ((fractionError(Pt) * sqrt(1/errorWeight(Pt, var)))^2 * weight(Pt, var)^2)
  //                      = Sum_Pt ((fractionError(Pt) * weight(Pt, var))^2 / errorWeight(Pt, var));
  // here, the fraction error in pT is weighted with a weight along the pT column (note that for FRACTIONS the error scales with
  // 1/sqrt(N) for not too small N).
  // NOTE: The denominator (sum (weights) in case of the mean and (sum (weights))^2 in case of the error) is one
  // by construction (normalised weights).
  //
  // NOTE: The same code can also be used for the to-pi-ratios (replace "fraction" by "ratio" in the following),
  // but NOT for yields (because then the scaling of sigma is sqrt(N) and not 1/sqrt(N) and also one does not want
  // to average the yields of different bins, but to add them).
  // In case of the yields, one rather takes the fractions (with errors defined as above) and scales (i.e. no error on the
  // multiplicative factor) with the yield in the var bin.
  
  if (!hFractionIDFFvar || !hFractionIDFFtrackPt || !hWeightingVarvsPt || !hErrorWeightingVarvsPt) {
    printf("Error calculateWeightedMean: Missing histos!\n");
    return;
  }
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hFractionIDFFvar[species])
      continue; // Will happen for to-pi-ratio in case of species = AliPID::kPion
    
    
    //TODO adapt errors?!
    //TODO adapt comment?!
    // Calculate estimate of probability density matrix
    const Int_t nBinsPt = hFractionIDFFtrackPt[species]->GetNbinsX();
    const Int_t nBinsVar = hFractionIDFFvar[species]->GetNbinsX();
    
    
    // Add 1 in each dimension to take into account that axes start at index 1
    Double_t probDensity[nBinsPt + 1][nBinsVar + 1];
    
    if (zerothIter) {
      // For zeroth iteration, guess of frac(var) extremely bad - do ordinary weighted mean first without any freaky weighting procedure
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          probDensity[binPt][binVar] = hFractionIDFFtrackPt[species]->GetBinContent(binPt);
        }
      }
    }
    else {
      // weightVarInv factors
      Double_t weightVarInv[nBinsPt + 1];
      for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
        weightVarInv[binVar] = 0;
        Double_t sumVarWeightMatrix = 0;
        for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
          weightVarInv[binVar] += hWeightingVarvsPt->GetBinContent(binPt, binVar) * hFractionIDFFtrackPt[species]->GetBinContent(binPt); 
          sumVarWeightMatrix += hWeightingVarvsPt->GetBinContent(binPt, binVar);
        }
        if (weightVarInv[binVar] <= 0) {
          // This case should only happen for e.g. a xi bin beyond the range. But this means that sumVarWeightMatrix is already zero.
          // If not, something is wrong
          if (sumVarWeightMatrix > 0) {
            printf("Error (species %s): weightVarInv[binVar] = %f <= 0 and sumVarWeightMatrix = %f > 0!\n", 
                  AliPID::ParticleShortName(species), weightVarInv[binVar], sumVarWeightMatrix);
            return;
          }
        }
      }
      
      // weightPtInv factors
      Double_t weightPtInv[nBinsPt + 1];
      weightPtInv[0] = 0;
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        weightPtInv[binPt] = 0;
        Double_t sumPtWeightMatrix = 0;
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          weightPtInv[binPt] += hErrorWeightingVarvsPt->GetBinContent(binPt, binVar) * hFractionIDFFvar[species]->GetBinContent(binVar); 
          sumPtWeightMatrix += hErrorWeightingVarvsPt->GetBinContent(binPt, binVar);
        }
        if (weightPtInv[binPt] <= 0) {
          // This case should only happen for e.g. a xi bin beyond the range. But this means that sumPtWeightMatrix is already zero.
          // If not, something is wrong
          if (sumPtWeightMatrix > 0) {
            printf("Error (species %s): weightPtInv[binPt] = %f <= 0 and sumPtWeightMatrix = %f > 0!\n", 
                  AliPID::ParticleShortName(species), weightPtInv[binPt], sumPtWeightMatrix);
            return;
          }
        }
      }
      
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          if (hWeightingVarvsPt->GetBinContent(binPt, binVar) > 0) {
            // weightPtInv[binPt] == 0 doesn't matter, if sumPtWeightMatrix is zero because matrix * weight is then still zero
            Double_t weightTotal = 0.;
            weightTotal = TMath::Sqrt(1./weightPtInv[binPt]*1./weightPtInv[binPt] + 1./weightVarInv[binVar]*1./weightVarInv[binVar]);
            //if ((weightPtInv[binPt] + weightVarInv[binVar]) > 0)  {
            //  weightTotal = 1. / (weightPtInv[binPt] + weightVarInv[binVar]);
            //}
            //else {
            //  printf("Error: Total weight is <= 0 for binPt/binVar = %d/%d\n", binPt, binVar);
            //}
            probDensity[binPt][binVar] = hFractionIDFFvar[species]->GetBinContent(binVar)
                                        * hFractionIDFFtrackPt[species]->GetBinContent(binPt)
                                        * weightTotal;
          }
          else
            probDensity[binPt][binVar] = 0.;
        }
      }
      
      // Normalise probDensity such that the overall (all pT, var) species fraction equals the one measured in all pT bins
      
      Double_t overallFractionSpecies = 0;
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          overallFractionSpecies += hErrorWeightingVarvsPt->GetBinContent(binPt, binVar)
                                    * hFractionIDFFtrackPt[species]->GetBinContent(binPt);
        }
      }
      
      Double_t overallFractionSpeciesRec = 0;
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          overallFractionSpeciesRec += probDensity[binPt][binVar] * hErrorWeightingVarvsPt->GetBinContent(binPt, binVar);
        }
      }

      if (overallFractionSpeciesRec <= 0) {
        printf("Yield for species %s is zero!\n", AliPID::ParticleShortName(species));
        return;
      }

      const Double_t normFactor = overallFractionSpecies / overallFractionSpeciesRec; 

      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
          probDensity[binPt][binVar] = probDensity[binPt][binVar] * normFactor;
        }
      }
    }
  
    for (Int_t binVar = 1; binVar <= nBinsVar; binVar++) {
      Double_t fraction = 0.;
      Double_t fractionError2 = 0.;
      
      // Binning of pT was adjusted such that all histos have the same
      for (Int_t binPt = 1; binPt <= nBinsPt; binPt++) {
        const Double_t weight = hWeightingVarvsPt->GetBinContent(binPt, binVar);
        fraction += probDensity[binPt][binVar] * weight;
        
        const Double_t errorWeight = hErrorWeightingVarvsPt->GetBinContent(binPt, binVar);
        if (errorWeight > 0) {
          fractionError2 += TMath::Power(hFractionIDFFtrackPt[species]->GetBinError(binPt) * weight, 2)
                            / errorWeight;
        }
        else {
          if (weight > 0) {
            printf("ERROR: errorWeight is %f. This should not happen if a finite yield (for any species) was found in this bin, i.e. there is a finite weight: %f!\n",
                   errorWeight, weight);
          }
        }
      }
      
      const Double_t fractionError = TMath::Sqrt(fractionError2);
      
      
      //// Smoothly change fraction w.r.t. to RELATIVE value: old + alpha * (new/old - 1) * old
      //const Double_t oldFraction = hFractionIDFFvar[species]->GetBinContent(binVar);
      //hFractionIDFFvar[species]->SetBinContent(binVar, oldFraction + alpha * (fraction - oldFraction));
      //hFractionIDFFvar[species]->SetBinError(binVar, fractionError);
      
      // Smoothly change fraction, i.e. only alpha * difference (new - old) to old
      const Double_t oldFraction = hFractionIDFFvar[species]->GetBinContent(binVar);
      hFractionIDFFvar[species]->SetBinContent(binVar, (1. - alpha) * oldFraction + alpha * fraction);
      hFractionIDFFvar[species]->SetBinError(binVar, fractionError);
    }
  }
}


//___________________________________________________________________
void calculateWeightedMeanIteratively(TH1D** hFractionIDFFvar, TH1D** hFractionIDFFtrackPt, TH2* hWeightingVarvsPt,
                                      TH2* hErrorWeightingVarvsPt, TFile* saveFile)
{
  // Calculate weighted mean iteratively: Firstly, flat priors vs. var are used to calculate the weighted mean.
  // In the next iterations, the probability matrix using result(var) * result(pT) * some weighting is calculated
  // and used to find the result(var) of the next iteration, which is then again the input of the following iteration.
  
  if (!hFractionIDFFvar || !hFractionIDFFtrackPt || !hWeightingVarvsPt || !hErrorWeightingVarvsPt) {
    printf("Error calculateWeightedMeanIteratively: Missing histos!\n");
    return;
  }
  
  // Initialise flat priors - actual value does not matter as long as it is constant vs. var
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hFractionIDFFvar[species])
      continue; // Will happen for to-pi-ratio in case of species = AliPID::kPion
    
    Double_t initialGuess = 1. / AliPID::kSPECIES;
    for (Int_t binVar = 1; binVar <= hFractionIDFFvar[species]->GetNbinsX(); binVar++) {
      hFractionIDFFvar[species]->SetBinContent(binVar, initialGuess);
      hFractionIDFFvar[species]->SetBinError(binVar, 1); // no information - large error
    }
  }
  
  const Int_t numIter = 3;
  const Double_t alpha = 0.5;
  for (Int_t iter = 0; iter <= numIter; iter++) {
    // Always take as input the results (vs. var) from the last iteration/initial guess
    
    // Set alpha to 1 for iter 0 because initial guess is much worse than that after iter 0
    calculateWeightedMean(hFractionIDFFvar, hFractionIDFFtrackPt, hWeightingVarvsPt, hErrorWeightingVarvsPt,
                          (iter == 0), (iter == 0) ? 1. : alpha);
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      // Save result of some iterations to file
      if (((Int_t)(iter * alpha * 100)) % 10 != 0)
        continue;
      
      if (hFractionIDFFvar[species]) {
        if (saveFile) {
          const TString subDir = "iterations";
          if (!saveFile->FindObject(subDir.Data()))
            saveFile->mkdir(subDir.Data());
          saveFile->cd(subDir.Data());
          
          TH1D* hFractionIDFFvarIter = new TH1D(*hFractionIDFFvar[species]);
          hFractionIDFFvarIter->SetName(Form("%s_iter%d", hFractionIDFFvar[species]->GetName(), iter));
          hFractionIDFFvarIter->Write();
          
          saveFile->cd();
        }
      }
    }
  }
  //TODO store var histos for each iteration (outside) => Check, whether procedure converges for all jetPt bins,
  // if so, develop convergence criterium
  //TODO Fix number of iterations at first => Adapt comment
}
*/

//___________________________________________________________________
void calculateWeightedMean(TH1D** hFractionIDFFvar, TH1D** hFractionIDFFtrackPt, TH2* hWeightingVarvsPt, TH2* hErrorWeightingVarvsPt)
{
  // Loop over all var (=z, xi) bins and calculate for the mean fraction: Fraction(var) = Sum_Pt (fraction(Pt) * weight(Pt, var)).
  // For the error, calculate:
  // FractionError(var)^2 = Sum_Pt ((fractionError(Pt) * sqrt(1/errorWeight(Pt, var)))^2 * weight(Pt, var)^2)
  //                      = Sum_Pt ((fractionError(Pt) * weight(Pt, var))^2 / errorWeight(Pt, var));
  // here, the fraction error in pT is weighted with a weight along the pT column (note that for FRACTIONS the error scales with
  // 1/sqrt(N) for not too small N).
  // NOTE: The denominator (sum (weights) in case of the mean and (sum (weights))^2 in case of the error) is one
  // by construction (normalised weights).
  //
  // NOTE: The same code can also be used for the to-pi-ratios (replace "fraction" by "ratio" in the following),
  // but NOT for yields (because then the scaling of sigma is sqrt(N) and not 1/sqrt(N) and also one does not want
  // to average the yields of different bins, but to add them).
  // In case of the yields, one rather takes the fractions (with errors defined as above) and scales (i.e. no error on the
  // multiplicative factor) with the yield in the var bin.
  
  if (!hFractionIDFFvar || !hFractionIDFFtrackPt || !hWeightingVarvsPt || !hErrorWeightingVarvsPt) {
    printf("Error calculateWeightedMean: Missing histos!\n");
    return;
  }
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    if (!hFractionIDFFvar[species])
      continue; // Will happen for to-pi-ratio in case of species = AliPID::kPion
    
    
    for (Int_t binVar = 1; binVar <= hFractionIDFFvar[species]->GetNbinsX(); binVar++) {
      Double_t fraction = 0.;
      Double_t fractionError2 = 0.;
      
      // Binning of pT was adjusted such that all histos have the same
      for (Int_t binPt = 1; binPt <= hFractionIDFFtrackPt[species]->GetNbinsX(); binPt++) {
        const Double_t weight = hWeightingVarvsPt->GetBinContent(binPt, binVar);
        fraction += hFractionIDFFtrackPt[species]->GetBinContent(binPt) * weight;
        
        const Double_t errorWeight = hErrorWeightingVarvsPt->GetBinContent(binPt, binVar);
        if (errorWeight > 0) {
          fractionError2 += TMath::Power(hFractionIDFFtrackPt[species]->GetBinError(binPt) * weight, 2)
                            / errorWeight;
        }
        else {
          if (weight > 0) {
            printf("ERROR: errorWeight is %f. This should not happen if a finite yield (for any species) was found in this bin, i.e. there is a finite weight: %f!\n",
                   errorWeight, weight);
          }
        }
      }
      
      const Double_t fractionError = TMath::Sqrt(fractionError2);
      hFractionIDFFvar[species]->SetBinContent(binVar, fraction);
      hFractionIDFFvar[species]->SetBinError(binVar, fractionError);
    }
  }
}


//___________________________________________________________________
void translateFractionToYield(TH1D** hIDFFvar, TH1D** hFractionIDFFvar, TH1D* hTotalWeight)
{
  // Take the fractions for var and just scale (i.e. no error for scale parameter) the
  // content of each bin with the total yield in this var bin.
  // The error of the fractions is already constructed such that it includes the error of the yield
  // (via log-likelihood fit vs. pT and then proper propagation from pT to var).
  //
  // NOTE: The binning must be the same for all histograms!
  
  if (!hIDFFvar || !hFractionIDFFvar || !hTotalWeight) {
    printf("Error translateFractionToYield: Missing histos!\n");
    return;
  }
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    for (Int_t binVar = 1; binVar <= hIDFFvar[species]->GetNbinsX(); binVar++) {
      const Double_t totalWeight = hTotalWeight->GetBinContent(binVar);
      const Double_t fraction = hFractionIDFFvar[species]->GetBinContent(binVar);
      const Double_t fractionError = hFractionIDFFvar[species]->GetBinError(binVar);
      
      hIDFFvar[species]->SetBinContent(binVar, fraction * totalWeight);
      hIDFFvar[species]->SetBinError(binVar, fractionError * totalWeight);
    }
  }
}


  /*
//___________________________________________________________________
void GenerateParticleYields(THnSparse* hDataProj, AliAnalysisTaskPID* pidTask, const Double_t centrality,
                            TH2D** hIDFFtrackPt, TH2D** hIDFFz, TH2D** hIDFFxi,
                            const Bool_t setMean, const Bool_t addErrorsQuadratically, const Bool_t smearByError,
                            const Bool_t uniformSystematicError, const Bool_t takeIntoAccountSysError, Int_t nGenerations)
{
  
  if (!hDataProj || !pidTask || !hIDFFtrackPt || !hIDFFz || !hIDFFxi) {
    printf("Cannot generate particle fractions - missing input!\n");
    return;
  }
  
  if (takeIntoAccountSysError && smearByError) {
    printf("It doesn't make sense to correlate statistical and systematic errors!\n");
    return;
  }
  
  const Int_t nDimProj = hDataProj->GetNdimensions();
  const Long64_t nBinsTHnSparseProj = hDataProj->GetNbins();
  Double_t binContent = 0., binError2 = 0.;
  Int_t binCoord[nDimProj];
  Double_t binCentre[nDimProj];
  Double_t prob[AliPID::kSPECIES];
  
  Bool_t success = kTRUE;
  
  // NOTE: In the following, the bin error is divided according to the fraction. The idea is to divide bin content (and error)
  // into independent sub-samples according to the fractions. But then, adding up the errors of the sub-samples quadratically
  // should recover the original error (think of 100 tracks, error sqrt(100), divided into 10 samples a 10 tracks => error of
  // samples should be sqrt(10) then).
  // This means, that one needs to add error^2 * fraction and NOT error^2 * fraction^2!!
  // However, the errors are ignored anyway at the moment....
  
  if (!smearByError && !takeIntoAccountSysError) {
    nGenerations = 1;
    // Just calculate fractions/spectra w/o any smearing
    const Int_t smearSpeciesByError = -1;
    const Int_t takeIntoAccountSysErrorOfSpecies = -1;
  
    for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
      binContent = hDataProj->GetBinContent(bin, binCoord);
      binError2  = hDataProj->GetBinError2(bin);
      
      // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
      if (binContent < 2. * std::numeric_limits<double>::min())
        continue;
      
      for (Int_t dim = 0; dim < nDimProj; dim++) 
        binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
      
      // Since there is no smearing, the fraction will stay the same and there is no need to save the result for one iteration in a
      // histogram (cf. case smearByError >= 0).
      success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob,
                                              smearSpeciesByError, takeIntoAccountSysErrorOfSpecies, uniformSystematicError);
      
      if (success) {
        for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
          // Track pT
          hIDFFtrackPt[species]->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                              hIDFFtrackPt[species]->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[species]);
          Double_t tempError = hIDFFtrackPt[species]->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
          hIDFFtrackPt[species]->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
          
          // z
          hIDFFz[species]->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                        hIDFFz[species]->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                        + binContent * prob[species]);
          tempError = hIDFFz[species]->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
          hIDFFz[species]->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                      TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
          
          // xi
          hIDFFxi[species]->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                          hIDFFxi[species]->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                          + binContent * prob[species]);
          tempError = hIDFFxi[species]->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
          hIDFFxi[species]->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                        TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
        }
      }
      else 
        printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f\n", binCentre[kPidPtProj], binCentre[kPidJetPtProj],
              centrality);
    }
  }
  else {
    // Calculate fractions/spectra w/ smearing "nGenerations" times for every species and obtain the statistical or systematic error 
    // from the spread of the results.
    // Note: For a given species, only the spread of the variation of exactly that species is considered, i.e. one is not
    // looking at the change of the fraction if another species is varied (this would bias the result towards smaller errors
    // by construction since the fraction of the considered species changes weighted with its statistical/systematic error)
    
    TH2D* hIDFFtrackPtVaried[AliPID::kSPECIES][nGenerations];
    TH2D* hIDFFzVaried[AliPID::kSPECIES][nGenerations];
    TH2D* hIDFFxiVaried[AliPID::kSPECIES][nGenerations];
    
    // Create this histo and set all entries to -1 (= not calculated yet) in every loop
    TH2D* hPartFractionsSmeared = hDataProj->Projection(kPidJetPtProj, kPidPtProj);
    hPartFractionsSmeared->SetName("hPartFractionsSmeared");
    
    // Generate results with varying fractions
    for (Int_t smearSpeciesByError = 0; smearSpeciesByError < AliPID::kSPECIES; smearSpeciesByError++) {
      for (Int_t i = 0; i < nGenerations; i++) {
        hIDFFtrackPtVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFtrackPt[smearSpeciesByError]);
        hIDFFtrackPtVaried[smearSpeciesByError][i]->SetName(Form("hIDFFtrackPtVaried_%d_%d", smearSpeciesByError, i));
        hIDFFtrackPtVaried[smearSpeciesByError][i]->Reset();
        
        hIDFFzVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFz[smearSpeciesByError]);
        hIDFFzVaried[smearSpeciesByError][i]->SetName(Form("hIDFFzVaried_%d_%d", smearSpeciesByError, i));
        hIDFFzVaried[smearSpeciesByError][i]->Reset();
        
        hIDFFxiVaried[smearSpeciesByError][i] = new TH2D(*(TH2D*)hIDFFxi[smearSpeciesByError]);
        hIDFFxiVaried[smearSpeciesByError][i]->SetName(Form("hIDFFxiVaried_%d_%d", smearSpeciesByError, i));
        hIDFFxiVaried[smearSpeciesByError][i]->Reset();
        
        // In each iteration (moving over all bins), one want only ONE fixed particle fraction per bin in the spectra map.
        // Otherwise, one would assing different fractions to e.g. the same trackPt bin, if only z and xi is different
        // (but jetPt bin is still the same). This would then cause the fluctuations to cancel partially.
        // Thus: Calculate the smeared fraction for each bin in the spectra map only once and store it in a histo. 
        // If it is requested again (i.e. histoEntry >= 0), use the value from the histo.
        
        // Set all entries to -1 (= not calculated yet)
        for (Int_t iX = 0; iX <= hPartFractionsSmeared->GetNbinsX(); iX++) {
          for (Int_t iY = 0; iY <= hPartFractionsSmeared->GetNbinsY(); iY++) {
            hPartFractionsSmeared->SetBinContent(iX, iY, -1);
          }
        }
        
        for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
          binContent = hDataProj->GetBinContent(bin, binCoord);
          binError2  = hDataProj->GetBinError2(bin);
          
          // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
          if (binContent < 2. * std::numeric_limits<double>::min())
            continue;
          
          for (Int_t iS = 0; iS < AliPID::kSPECIES; iS++)
            prob[iS] = 0.;
          
          for (Int_t dim = 0; dim < nDimProj; dim++) 
            binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
          
          if (hPartFractionsSmeared->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj]) >= 0) {
            success = kTRUE;
            prob[smearSpeciesByError] = hPartFractionsSmeared->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
          }
          else {
            const Int_t smearSpeciesByStatisticalError = smearByError ? smearSpeciesByError : -1;
            const Int_t takeIntoAccountSysErrorOfSpecies = takeIntoAccountSysError ? smearSpeciesByError : -1;
            
            success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob, 
                                                    smearSpeciesByStatisticalError, takeIntoAccountSysErrorOfSpecies, 
                                                    uniformSystematicError);
            
            if (success)
              hPartFractionsSmeared->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], prob[smearSpeciesByError]);
          }
          
          // NOTE: Since only the bin content of the current species is stored in hPartFractionsSmeared,
          // only prob[smearSpeciesByError] can be used in the following!!!!
          
          if (success) {
            // To make things readable...
            TH2D* hIDFFtrackPtVariedCurr = hIDFFtrackPtVaried[smearSpeciesByError][i];
            TH2D* hIDFFzVariedCurr       = hIDFFzVaried[smearSpeciesByError][i];
            TH2D* hIDFFxiVariedCurr      = hIDFFxiVaried[smearSpeciesByError][i];
            
            // Track pT
            hIDFFtrackPtVariedCurr->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                  hIDFFtrackPtVariedCurr->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                                  + binContent * prob[smearSpeciesByError]);
            Double_t tempError = hIDFFtrackPtVariedCurr->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
            hIDFFtrackPtVariedCurr->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
            
            // z
            hIDFFzVariedCurr->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                            hIDFFzVariedCurr->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                            + binContent * prob[smearSpeciesByError]);
            tempError = hIDFFzVariedCurr->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
            hIDFFzVariedCurr->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                          TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
            
            // xi
            hIDFFxiVariedCurr->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                              hIDFFxiVariedCurr->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[smearSpeciesByError]);
            tempError = hIDFFxiVariedCurr->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
            hIDFFxiVariedCurr->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[smearSpeciesByError]));
          }
          else 
            printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f, smearSpeciesByError %d\n", 
                   binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, smearSpeciesByError);
        }
      }
    }
    
    delete hPartFractionsSmeared;
    
    const Int_t nHistos = nGenerations;
    
    */
    /*OLD error for each species by variation of ALL species (will bias towards smaller errors!)
    // TODO Still does not store all the values in histogram to avoid different fractions for same map bin in same iteration.
    // => If this is build in, one needs a histogram for EACH species!!!

    // Calculate fractions/spectra w/ smearing "nGenerations" times for every species and obtain the statistical error from
    // the spread of the results
    
    TH2D* hIDFFtrackPtVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    TH2D* hIDFFzVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    TH2D* hIDFFxiVaried[AliPID::kSPECIES][AliPID::kSPECIES * nGenerations];
    
    // Generate results with varying fractions
    for (Int_t smearSpeciesByError = 0; smearSpeciesByError < AliPID::kSPECIES; smearSpeciesByError++) {
      for (Int_t i = 0; i < nGenerations; i++) {
        for (Int_t consideredSpecies = 0; consideredSpecies < AliPID::kSPECIES; consideredSpecies++) {
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFtrackPt[consideredSpecies]);
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFtrackPtVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFtrackPtVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFz[consideredSpecies]);
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFzVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFzVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i] = 
            new TH2D(*(TH2D*)hIDFFxi[consideredSpecies]);
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->Reset();
          hIDFFxiVaried[consideredSpecies][smearSpeciesByError * nGenerations + i]->SetName(Form("hIDFFxiVaried_%d_%d",
                                                                                          consideredSpecies,
                                                                                          smearSpeciesByError * nGenerations + i));
        }
        
        for (Long64_t bin = 0; bin < nBinsTHnSparseProj; bin++) {
          binContent = hDataProj->GetBinContent(bin, binCoord);
          binError2  = hDataProj->GetBinError2(bin);
          
          // If the bin is empty, do not compute the particle fractions -> This bin contributes nothing anyway
          if (binContent < 2. * std::numeric_limits<double>::min())
            continue;
          
          for (Int_t dim = 0; dim < nDimProj; dim++) 
            binCentre[dim] = hDataProj->GetAxis(dim)->GetBinCenter(binCoord[dim]);
          
          success = pidTask->GetParticleFractions(binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, prob, 
                                                  smearSpeciesByError);
          
          if (success) {
            for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
              // To make things readable...
              TH2D* hIDFFtrackPtVariedCurr = hIDFFtrackPtVaried[species][smearSpeciesByError * nGenerations + i];
              TH2D* hIDFFzVariedCurr       = hIDFFzVaried[species][smearSpeciesByError * nGenerations + i];
              TH2D* hIDFFxiVariedCurr      = hIDFFxiVaried[species][smearSpeciesByError * nGenerations + i];
              
              // Track pT
              hIDFFtrackPtVariedCurr->SetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                    hIDFFtrackPtVariedCurr->GetBinContent(binCoord[kPidPtProj], binCoord[kPidJetPtProj])
                                                    + binContent * prob[species]);
              Double_t tempError = hIDFFtrackPtVariedCurr->GetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj]);
              hIDFFtrackPtVariedCurr->SetBinError(binCoord[kPidPtProj], binCoord[kPidJetPtProj], 
                                                  TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
              
              // z
              hIDFFzVariedCurr->SetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj], 
                                              hIDFFzVariedCurr->GetBinContent(binCoord[kPidZProj], binCoord[kPidJetPtProj])
                                              + binContent * prob[species]);
              tempError = hIDFFzVariedCurr->GetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj]);
              hIDFFzVariedCurr->SetBinError(binCoord[kPidZProj], binCoord[kPidJetPtProj],
                                            TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
              
              // xi
              hIDFFxiVariedCurr->SetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj], 
                                               hIDFFxiVariedCurr->GetBinContent(binCoord[kPidXiProj], binCoord[kPidJetPtProj])
                                               + binContent * prob[species]);
              tempError = hIDFFxiVariedCurr->GetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj]);
              hIDFFxiVariedCurr->SetBinError(binCoord[kPidXiProj], binCoord[kPidJetPtProj],
                                             TMath::Sqrt(tempError * tempError + binError2 * prob[species]));
            }
          }
          else 
            printf("Problem obtaining fractions for: trackPt %f, jetPt %f, cent %f, smearSpeciesByError %d\n", 
                   binCentre[kPidPtProj], binCentre[kPidJetPtProj], centrality, smearSpeciesByError);
        }
      }
    }
    
    const Int_t nHistos = AliPID::kSPECIES * nGenerations;
    */
    /*
    // Compare results to obtain error
    const Bool_t ignoreSigmaErrors = kTRUE;
    
    for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
      if (!extractSystematicError(nHistos, hIDFFtrackPtVaried[species], hIDFFtrackPt[species], setMean, addErrorsQuadratically,
                                  ignoreSigmaErrors))
        printf("Failed to determine systematic error for trackPt, species %d\n", species);
      
      if (!extractSystematicError(nHistos, hIDFFzVaried[species], hIDFFz[species], setMean, addErrorsQuadratically, ignoreSigmaErrors))
        printf("Failed to determine systematic error for z, species %d\n", species);
      
      if (!extractSystematicError(nHistos, hIDFFxiVaried[species], hIDFFxi[species], setMean, addErrorsQuadratically, 
                                  ignoreSigmaErrors))
        printf("Failed to determine systematic error for xi, species %d\n", species);
      
      for (Int_t i = 0; i < nHistos; i++) {
        delete hIDFFtrackPtVaried[species][i];
        hIDFFtrackPtVaried[species][i] = 0x0;
        
        delete hIDFFzVaried[species][i];
        hIDFFzVaried[species][i] = 0x0;
        
        delete hIDFFxiVaried[species][i];
        hIDFFxiVaried[species][i] = 0x0;
      }
    }
  }
}
*/
  
  
//___________________________________________________________________
Int_t extractFFs(TString pathNameData, TString listName /*= ""*/, TString pathNameFractionsAndYields,
                 Int_t chargeMode /*kNegCharge = -1, kAllCharged = 0, kPosCharge = 1*/,
                 Double_t lowerCentrality = -2, Double_t upperCentrality = -2, Double_t lowerJetPt = -1, Double_t upperJetPt = -1,
                 Int_t rebinZ = 1, Int_t rebinXi = 1, Bool_t onlyUseRelevantMCIDforMatrix = kFALSE)
                 // NOTE: rebinPt makes no sense, since one takes the pT histos from the file (histos already fixed w.r.t. binning)
{
  TObjArray* histList = 0x0;
  
  if (listName == "") {
    listName = pathNameData;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  // Load pT fractions, yields and to-pi-ratios (data and MC)
  TFile* fFractionsAndYields = TFile::Open(pathNameFractionsAndYields.Data());
  if (!fFractionsAndYields)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameFractionsAndYields.Data() << "\"!" << std::endl;
    return -1;
  }
  
  
  TH1D* hFractionIDFFtrackPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hFractionIDFFtrackPtMC[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hIDFFtrackPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hIDFFtrackPtMC[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hRatioToPiIDFFtrackPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiIDFFtrackPtMC[AliPID::kSPECIES] = { 0x0, };

  Int_t numMCHistsFound = 0;
  
  // In case of MC also retrieve the MC (raw) yields and fractions
  // NOTE: The MC histograms are also there for data, but then use the "most probable PID".
  // In case of data, just ignore the MC histos.
  Bool_t hasMC = kFALSE;
  TH1D* hTestMC = (TH1D*)fFractionsAndYields->Get("hFractionComparisonPions");
  if (hTestMC) {
    TString yTitle = hTestMC->GetYaxis()->GetTitle();
    if (yTitle.Contains("Most Probable PID", TString::kIgnoreCase) == kFALSE)
      hasMC = kTRUE;
  }
  
  //printf("TODO: MC manually set to be present!\n");
  //hasMC = kTRUE; //TODO
  
  if (hasMC) 
    printf("Fraction file seems to contain MC histos -> Trying do extract them...\n\n");
  else
    printf("Fraction file does not seem to contain MC histos -> Ignoring them...\n\n");
    
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    TString speciesName = AliPID::ParticleName(species);
    TString firstLetter = speciesName(0);
    firstLetter.ToUpper();
    speciesName.Replace(0, 1, firstLetter.Data());
    TString histName = Form("hYield%ss", speciesName.Data());
    hIDFFtrackPt[species] = (TH1D*)fFractionsAndYields->Get(histName.Data());
    if (!hIDFFtrackPt[species]) {
      printf("Failed to load hist \"%s\"\n", histName.Data());
      return -1;
    }
    hIDFFtrackPt[species]->SetFillStyle(0);
    hIDFFtrackPt[species]->SetName(Form("hIDFFtrackPt_%s", AliPID::ParticleShortName(species)));
    
    TString histNameFraction = histName;
    histNameFraction.ReplaceAll("Yield", "Fraction");
    hFractionIDFFtrackPt[species] = (TH1D*)fFractionsAndYields->Get(histNameFraction.Data());
    if (!hFractionIDFFtrackPt[species]) {
      printf("Failed to load hist \"%s\"\n", histNameFraction.Data());
      return -1;
    }
    hFractionIDFFtrackPt[species]->SetFillStyle(0);
    hFractionIDFFtrackPt[species]->SetName(Form("hFractionIDFFtrackPt_%s", AliPID::ParticleShortName(species)));
    
    
    TString histNameRatioToPi = Form("hRatioToPi%ss", speciesName.Data());
    if (species != AliPID::kPion) {
      hRatioToPiIDFFtrackPt[species] = (TH1D*)fFractionsAndYields->Get(histNameRatioToPi.Data());
      if (!hRatioToPiIDFFtrackPt[species]) {
        printf("Failed to load hist \"%s\"\n", histNameRatioToPi.Data());
        //return -1; // For old data, these histos are not there, just go on without them
      }
      else {
        hRatioToPiIDFFtrackPt[species]->SetFillStyle(0);
        hRatioToPiIDFFtrackPt[species]->SetName(Form("hRatioToPiIDFFtrackPt_%s", AliPID::ParticleShortName(species)));
      }
    }
    
    if (hasMC) {
      TString histNameMC = Form("%sMC", histName.Data());
      hIDFFtrackPtMC[species] = (TH1D*)fFractionsAndYields->Get(histNameMC.Data());
      
      TString histNameFractionMC = Form("%sMC", histNameFraction.Data());
      hFractionIDFFtrackPtMC[species] = (TH1D*)fFractionsAndYields->Get(histNameFractionMC.Data());
      
      if (species != AliPID::kPion) {
        TString histNameRatioToPiMC = Form("%sMC", histNameRatioToPi.Data());
        hRatioToPiIDFFtrackPtMC[species] = (TH1D*)fFractionsAndYields->Get(histNameRatioToPiMC.Data());
      }
      
      if (hFractionIDFFtrackPtMC[species] && hIDFFtrackPtMC[species] &&
          (hRatioToPiIDFFtrackPtMC[species] || !hRatioToPiIDFFtrackPt[species] || species == AliPID::kPion)) {
        numMCHistsFound++;
        
        hIDFFtrackPtMC[species]->SetFillStyle(0);
        hIDFFtrackPtMC[species]->SetName(Form("hIDFFtrackPtMC_%s", AliPID::ParticleShortName(species)));
        
        hFractionIDFFtrackPtMC[species]->SetFillStyle(0);
        hFractionIDFFtrackPtMC[species]->SetName(Form("hFractionIDFFtrackPtMC_%s", AliPID::ParticleShortName(species)));
      
        if (species != AliPID::kPion && hRatioToPiIDFFtrackPtMC[species]) {
          hRatioToPiIDFFtrackPtMC[species]->SetFillStyle(0);
          hRatioToPiIDFFtrackPtMC[species]->SetName(Form("hRatioToPiIDFFtrackPtMC_%s", AliPID::ParticleShortName(species)));
        }
      }
    }
  }
  
  if (numMCHistsFound > 0 && numMCHistsFound != AliPID::kSPECIES) {
    printf("Error: Unable to retrieve all MC histos! Got %d.\n", numMCHistsFound);
    return -1;
  }
  
  
  
  // Extract the data histogram
  TFile* f = TFile::Open(pathNameData.Data());
  if (!f)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\"!" << std::endl;
    return -1;
  }
  
  histList = (TObjArray*)(f->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  THnSparse* hPIDdata = dynamic_cast<THnSparse*>(histList->FindObject("hPIDdataAll"));
  if (!hPIDdata) {
    std::cout << std::endl;
    std::cout << "Failed to load data histo!" << std::endl;
    return -1;
  }
  
  // Set proper errors, if not yet calculated
  if (!hPIDdata->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << hPIDdata->GetName() << "..." << std::endl;
    hPIDdata->Sumw2();
    Long64_t nBinsTHnSparse = hPIDdata->GetNbins();
    Double_t binContent = 0;
    
    for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
      binContent = hPIDdata->GetBinContent(bin);
      hPIDdata->SetBinError(bin, TMath::Sqrt(binContent));
    }
  }
  
  // If desired, restrict centrality axis
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2; // Integral(lowerCentBinLimit, uppCentBinLimit) will not be restricted if these values are kept
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hPIDdata->GetAxis(kPidCentrality)->FindBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= hPIDdata->GetAxis(kPidCentrality)->GetNbins()) {
      actualLowerCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hPIDdata->GetAxis(kPidCentrality)->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
  }
  else {
    std::cout << "All" << std::endl;
  }
    
  if (restrictCentralityAxis) {
    hPIDdata->GetAxis(kPidCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  }
  
  // If desired, restrict jetPt axis
  Int_t lowerJetPtBinLimit = -1;
  Int_t upperJetPtBinLimit = -1;
  Bool_t restrictJetPtAxis = kFALSE;
  Double_t actualLowerJetPt = -1.;
  Double_t actualUpperJetPt = -1.;
  
  if (lowerJetPt >= 0 && upperJetPt >= 0) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(lowerJetPt + 0.001);
    upperJetPtBinLimit = hPIDdata->GetAxis(kPidJetPt)->FindBin(upperJetPt - 0.001);
    
    // Check if the values look reasonable
    if (lowerJetPtBinLimit <= upperJetPtBinLimit && lowerJetPtBinLimit >= 1 && upperJetPtBinLimit <= hPIDdata->GetAxis(kPidJetPt)->GetNbins()) {
      actualLowerJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinLowEdge(lowerJetPtBinLimit);
      actualUpperJetPt = hPIDdata->GetAxis(kPidJetPt)->GetBinUpEdge(upperJetPtBinLimit);

      restrictJetPtAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested jet pT range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  std::cout << "jet pT: ";
  if (restrictJetPtAxis) {
    std::cout << actualLowerJetPt << " - " << actualUpperJetPt << std::endl;
  }
  else {
    std::cout << "All" << std::endl;
  }
  
  if (restrictJetPtAxis) {
    hPIDdata->GetAxis(kPidJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
  }
  
  // If desired, restrict charge axis
  std::cout << "Charge selection: ";
  if (chargeMode == kAllCharged)
    std::cout << "All charged particles" << std::endl;
  else if (chargeMode == kNegCharge)
    std::cout << "Negative particles only" << std::endl;
  else if (chargeMode == kPosCharge)
    std::cout << "Positive particles only" << std::endl;
  else {
    std::cout << "Unknown -> ERROR" << std::endl;
    return -1;
  }
  
  const Bool_t restrictCharge = (chargeMode != kAllCharged);
  
  const Int_t indexChargeAxisData = GetAxisByTitle(hPIDdata, "Charge (e_{0})");
  if (indexChargeAxisData < 0 && restrictCharge) {
    std::cout << "Error: Charge axis not found for data histogram!" << std::endl;
    return -1;
  }
  Int_t lowerChargeBinLimitData = -1;
  Int_t upperChargeBinLimitData = -2;
  Double_t actualLowerChargeData = -999;
  Double_t actualUpperChargeData = -999;
  
  if (restrictCharge) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    if (chargeMode == kNegCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(-1. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. - 0.001);
    }
    else if (chargeMode == kPosCharge) {
      lowerChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(0. + 0.001);
      upperChargeBinLimitData = hPIDdata->GetAxis(indexChargeAxisData)->FindBin(1. - 0.001);
    }
    
    // Check if the values look reasonable
    if (lowerChargeBinLimitData <= upperChargeBinLimitData && lowerChargeBinLimitData >= 1
        && upperChargeBinLimitData <= hPIDdata->GetAxis(indexChargeAxisData)->GetNbins()) {
      actualLowerChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinLowEdge(lowerChargeBinLimitData);
      actualUpperChargeData = hPIDdata->GetAxis(indexChargeAxisData)->GetBinUpEdge(upperChargeBinLimitData);
      
      std::cout << "Charge range data: " << actualLowerChargeData << " - " << actualUpperChargeData << std::endl;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested charge range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
    
    hPIDdata->GetAxis(indexChargeAxisData)->SetRange(lowerChargeBinLimitData, upperChargeBinLimitData);
  }
  
  // Just take one arbitrary selectSpecies to avoid multiple counting
  hPIDdata->GetAxis(kPidSelectSpecies)->SetRange(1, 1);
  
  // If desired, throw away under- and overflow bins from MC ID
  if (onlyUseRelevantMCIDforMatrix) {
    hPIDdata->GetAxis(kPidMCpid)->SetRange(1, hPIDdata->GetAxis(kPidMCpid)->GetNbins());
  }
  
  // Get projection on dimensions relevant for FFs
  const Int_t nDimProj = kNDimsProj;
  Int_t dimProj[nDimProj];
  dimProj[kPidPtProj] = kPidPt;
  dimProj[kPidZProj] = kPidZ;
  dimProj[kPidXiProj] =kPidXi;
  
  THnSparse* hDataProj = hPIDdata->Projection(nDimProj, dimProj, "e");
  
  // If desired, rebin axes
  if (rebinZ != 1 || rebinXi != 1) {
    Int_t group[hDataProj->GetNdimensions()];
    
    for (Int_t i = 0; i < hDataProj->GetNdimensions(); i++) {
       if (i == kPidZProj)
        group[i] = rebinZ;
      else if (i == kPidXiProj)
        group[i] = rebinXi;
      else
        group[i] = 1;
    }
    
    THnSparse* hTemp = hDataProj->Rebin(group);
    hDataProj->SetName("temp");
    delete hDataProj;

    hDataProj = hTemp;
  }
  
  // Take pT binning from pion fraction (binning for all species the same) and create a new THnSparse with this new binning
  TH1D hDummyPt(*hFractionIDFFtrackPt[AliPID::kPion]);
  hDummyPt.SetName("hDummyPt");
  TAxis* axisPt = hDummyPt.GetXaxis();
  
  Int_t binsProj[nDimProj];
  Double_t xminProj[nDimProj];
  Double_t xmaxProj[nDimProj];
  
  for (Int_t iDim = 0; iDim < nDimProj; iDim++) {
    if (iDim == kPidPtProj) {
      binsProj[iDim] = axisPt->GetNbins();
      xminProj[iDim] = axisPt->GetBinLowEdge(1);
      xmaxProj[iDim] = axisPt->GetBinUpEdge(axisPt->GetNbins());
    }
    else {
      binsProj[iDim] = hDataProj->GetAxis(iDim)->GetNbins();
      xminProj[iDim] = hDataProj->GetAxis(iDim)->GetBinLowEdge(1);
      xmaxProj[iDim] = hDataProj->GetAxis(iDim)->GetBinUpEdge(hDataProj->GetAxis(iDim)->GetNbins());
    }
  }
  
  THnSparse* hDataProjRebinned = new THnSparseD("hDataProjRebinned","", nDimProj, binsProj, xminProj, xmaxProj);
  hDataProjRebinned->SetName("hDataProjRebinned");
  hDataProjRebinned->Sumw2();
  
  for (Int_t iDim = 0; iDim < nDimProj; iDim++) {
    if (iDim == kPidPtProj) {
      if (axisPt->GetXbins()->fN != 0)
        hDataProjRebinned->SetBinEdges(iDim, axisPt->GetXbins()->fArray);
    }
    else {
      if (hDataProj->GetAxis(iDim)->GetXbins()->fN != 0)
        hDataProjRebinned->SetBinEdges(iDim, hDataProj->GetAxis(iDim)->GetXbins()->fArray);
    }
    
    hDataProjRebinned->GetAxis(iDim)->SetTitle(hDataProj->GetAxis(iDim)->GetTitle());
  }
  
  // Now just fill the THnSparse with the original data. RebinnedAdd already takes into account the different binning properly
  // (was also tested explicitely!).
  hDataProjRebinned->RebinnedAdd(hDataProj, 1.);
  
  // Delete the old THnSparse and set the pointer to the new one
  delete hDataProj;
  hDataProj = hDataProjRebinned;
  
  TH2D* hNjetsGen = 0x0;
  TH2D* hNjetsRec = 0x0;
  
  TH1D* hMCgenPrimYieldPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimYieldZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimYieldXi[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hMCgenPrimYieldTotalPt = 0x0;
  TH1D* hMCgenPrimYieldTotalZ = 0x0;
  TH1D* hMCgenPrimYieldTotalXi = 0x0;
  
  TH1D* hMCgenPrimFractionPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimFractionZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimFractionXi[AliPID::kSPECIES] = { 0x0, };
  
  hNjetsGen = (TH2D*)histList->FindObject("fh2FFJetPtGen");
  hNjetsRec = (TH2D*)histList->FindObject("fh2FFJetPtRec");
  
  if (!hNjetsRec) {
    printf("Failed to load number of jets (rec) histo!\n");
    
    // For backward compatibility (TODO REMOVE IN FUTURE): Load info from fixed AnalysisResults file (might be wrong, if other
    // period is considered; also: No multiplicity information)
    TFile* fBackward = TFile::Open("finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/AnalysisResults.root");
    
    TString dirDataInFile = "";
    TDirectory* dirData = fBackward ? (TDirectory*)fBackward->Get(fBackward->GetListOfKeys()->At(0)->GetName()) : 0x0;
  
    TList* list = dirData ? (TList*)dirData->Get(dirData->GetListOfKeys()->At(0)->GetName()) : 0x0;

    TH1D* hFFJetPtRec = list ? (TH1D*)list->FindObject("fh1FFJetPtRecCutsInc") : 0x0;
    
    if (hFFJetPtRec) {
      printf("***WARNING: For backward compatibility, using file \"finalCuts/pp/7TeV/LHC10e.pass2/corrected/finalisedSplines/finalMapsAndTail/Jets/noCutOn_ncl_or_liav/AnalysisResults.root\" to get number of jets. BUT: Might be wrong period and has no mult info!***\n");
      printf("ALSO: Using Njets for inclusive jets!!!!\n");
      
      hNjetsRec = new TH2D("fh2FFJetPtRec", "", 1, -1, 1,  hPIDdata->GetAxis(kPidJetPt)->GetNbins(),
                          hPIDdata->GetAxis(kPidJetPt)->GetXbins()->GetArray());
      
      for (Int_t iJet = 1; iJet <= hNjetsRec->GetNbinsY(); iJet++) {
        Int_t lowerBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinLowEdge(iJet) + 1e-3);
        Int_t upperBin = hFFJetPtRec->FindFixBin(hNjetsRec->GetYaxis()->GetBinUpEdge(iJet) - 1e-3);
        hNjetsRec->SetBinContent(1, iJet, hFFJetPtRec->Integral(lowerBin, upperBin));
      }
    }
    
    if (!hNjetsRec)
      return -1;
  }
  
  THnSparse* hMCgeneratedYieldsPrimaries = (THnSparse*)histList->FindObject("fhMCgeneratedYieldsPrimaries");
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hMCgenPrimYieldPt[i] = 0x0;
    hMCgenPrimYieldZ[i] = 0x0;
    hMCgenPrimYieldXi[i] = 0x0;
  }
  
  if (hMCgeneratedYieldsPrimaries && hMCgeneratedYieldsPrimaries->GetEntries() > 0) {
    // Set proper errors, if not yet calculated
    if (!hMCgeneratedYieldsPrimaries->GetCalculateErrors()) {
      std::cout << "Re-calculating errors of " << hMCgeneratedYieldsPrimaries->GetName() << "..." << std::endl;
      
      hMCgeneratedYieldsPrimaries->Sumw2();
      
      Long64_t nBinsTHnSparseGenYield = hMCgeneratedYieldsPrimaries->GetNbins();
      Double_t binContentGenYield = 0;
      for (Long64_t bin = 0; bin < nBinsTHnSparseGenYield; bin++) {
        binContentGenYield = hMCgeneratedYieldsPrimaries->GetBinContent(bin);
        hMCgeneratedYieldsPrimaries->SetBinError(bin, TMath::Sqrt(binContentGenYield));
      }
    }
    
    // If desired, rebin axes
    if (rebinZ != 1 || rebinXi != 1) {
      Int_t group[hMCgeneratedYieldsPrimaries->GetNdimensions()];
      
      for (Int_t k = 0; k < hMCgeneratedYieldsPrimaries->GetNdimensions(); k++) {
        if (k == kPidGenYieldZ)
          group[k] = rebinZ;
        else if (k == kPidGenYieldXi)
          group[k] = rebinXi;
        else
          group[k] = 1;
      }
      
      THnSparse* hTemp = hMCgeneratedYieldsPrimaries->Rebin(group);
      hMCgeneratedYieldsPrimaries->SetName("temp");
      delete hMCgeneratedYieldsPrimaries;
      
      hMCgeneratedYieldsPrimaries = hTemp;
    }


    if (restrictCentralityAxis)
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldCentrality)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    
    if (restrictJetPtAxis) 
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldJetPt)->SetRange(lowerJetPtBinLimit, upperJetPtBinLimit);
    
    if (restrictCharge) {
      const Int_t indexChargeAxisGenYield = GetAxisByTitle(hMCgeneratedYieldsPrimaries, "Charge (e_{0})");
      if (indexChargeAxisGenYield < 0) {
        std::cout << "Error: Charge axis not found for gen yield histogram!" << std::endl;
        return -1;
      }
  
      Int_t lowerChargeBinLimitGenYield = -1;
      Int_t upperChargeBinLimitGenYield = -2;
      Double_t actualLowerChargeGenYield = -999;
      Double_t actualUpperChargeGenYield = -999;
  
      // Add subtract a very small number to avoid problems with values right on the border between to bins
      if (chargeMode == kNegCharge) {
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(-1. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. - 0.001);
      }
      else if (chargeMode == kPosCharge) {
        lowerChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(0. + 0.001);
        upperChargeBinLimitGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->FindBin(1. - 0.001);
      }
      
      // Check if the values look reasonable
      if (lowerChargeBinLimitGenYield <= upperChargeBinLimitGenYield && lowerChargeBinLimitGenYield >= 1
          && upperChargeBinLimitGenYield <= hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetNbins()) {
        actualLowerChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinLowEdge(lowerChargeBinLimitGenYield);
        actualUpperChargeGenYield = hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->GetBinUpEdge(upperChargeBinLimitGenYield);
        
        if (TMath::Abs(actualLowerChargeGenYield - actualLowerChargeData) > 1e-4 ||
            TMath::Abs(actualUpperChargeGenYield - actualUpperChargeData) > 1e-4) {
          std::cout << std::endl;
          std::cout << "Error: Charge range gen yield: " << actualLowerChargeGenYield << " - " << actualUpperChargeGenYield
                    << std::endl << "differs from that of data: " << actualLowerChargeData << " - " << actualUpperChargeData
                    << std::endl;
          return -1;
        }
      }
      else {
        std::cout << std::endl;
        std::cout << "Requested charge range (gen yield) out of limits or upper and lower limit are switched!" << std::endl;
        return -1;
      }
      
      hMCgeneratedYieldsPrimaries->GetAxis(indexChargeAxisGenYield)->SetRange(lowerChargeBinLimitGenYield, 
                                                                              upperChargeBinLimitGenYield);
    }
  
    for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(MCid + 1, MCid + 1);
      
      hMCgenPrimYieldPt[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldPt, "e");
      setupHist(hMCgenPrimYieldPt[MCid],
                Form("hMCgenYieldsPrimPt_%s", AliPID::ParticleShortName(MCid)),
                Form("%s", AliPID::ParticleName(MCid)),
                     Form("%s", hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldPt)->GetTitle()),
                     yieldAxisTitlePt.Data(), getLineColorAliPID(MCid), kTRUE);
      hMCgenPrimYieldPt[MCid]->SetMarkerStyle(28);
      
      hMCgenPrimYieldZ[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldZ, "e");
      setupHist(hMCgenPrimYieldZ[MCid],
                Form("hMCgenYieldsPrimZ_%s", AliPID::ParticleShortName(MCid)),
                Form("%s", AliPID::ParticleName(MCid)),
                     Form("%s", hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldZ)->GetTitle()),
                     yieldAxisTitleZ.Data(), getLineColorAliPID(MCid), kTRUE);
      hMCgenPrimYieldZ[MCid]->SetMarkerStyle(28);
      
      hMCgenPrimYieldXi[MCid] = hMCgeneratedYieldsPrimaries->Projection(kPidGenYieldXi, "e");
      setupHist(hMCgenPrimYieldXi[MCid],
                Form("hMCgenYieldsPrimXi_%s", AliPID::ParticleShortName(MCid)),
                Form("%s", AliPID::ParticleName(MCid)),
                Form("%s", hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldXi)->GetTitle()),
                yieldAxisTitleXi.Data(), getLineColorAliPID(MCid), kTRUE);
      hMCgenPrimYieldXi[MCid]->SetMarkerStyle(28);
      
      hMCgeneratedYieldsPrimaries->GetAxis(kPidGenYieldMCpid)->SetRange(0, -1);
    }
    
    // Total generated yields as sum of identified yields (i.e. no other species)
    hMCgenPrimYieldTotalPt = new TH1D(*hMCgenPrimYieldPt[0]);
    hMCgenPrimYieldTotalPt->SetName("hMCgenYieldsPrimPt_total");
    hMCgenPrimYieldTotalPt->SetTitle("Total");
    
    hMCgenPrimYieldTotalZ = new TH1D(*hMCgenPrimYieldZ[0]);
    hMCgenPrimYieldTotalZ->SetName("hMCgenYieldsPrimZ_total");
    hMCgenPrimYieldTotalZ->SetTitle("Total");
    
    hMCgenPrimYieldTotalXi = new TH1D(*hMCgenPrimYieldXi[0]);
    hMCgenPrimYieldTotalXi->SetName("hMCgenYieldsPrimXi_total");
    hMCgenPrimYieldTotalXi->SetTitle("Total");
    
    // Errors are correct since MC yields truly independent
    for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
      hMCgenPrimYieldTotalPt->Add(hMCgenPrimYieldPt[MCid]);
      hMCgenPrimYieldTotalZ->Add(hMCgenPrimYieldZ[MCid]);
      hMCgenPrimYieldTotalXi->Add(hMCgenPrimYieldXi[MCid]);
    }
    
    // Calculate fractions: Binomial error, since numerator is subset of denominator
    for (Int_t MCid = 0; MCid < AliPID::kSPECIES; MCid++) {
      hMCgenPrimFractionPt[MCid] = new TH1D(*hMCgenPrimYieldPt[MCid]);
      hMCgenPrimFractionPt[MCid]->SetName(Form("hMCgenFractionPrimPt_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimFractionPt[MCid]->GetYaxis()->SetTitle("Particle Fraction");
      hMCgenPrimFractionPt[MCid]->Divide(hMCgenPrimYieldPt[MCid], hMCgenPrimYieldTotalPt, 1., 1., "B");
      
      hMCgenPrimFractionZ[MCid] = new TH1D(*hMCgenPrimYieldZ[MCid]);
      hMCgenPrimFractionZ[MCid]->SetName(Form("hMCgenFractionPrimZ_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimFractionZ[MCid]->GetYaxis()->SetTitle("Particle Fraction");
      hMCgenPrimFractionZ[MCid]->Divide(hMCgenPrimYieldZ[MCid], hMCgenPrimYieldTotalZ, 1., 1., "B");
      
      hMCgenPrimFractionXi[MCid] = new TH1D(*hMCgenPrimYieldXi[MCid]);
      hMCgenPrimFractionXi[MCid]->SetName(Form("hMCgenFractionPrimXi_%s", AliPID::ParticleShortName(MCid)));
      hMCgenPrimFractionXi[MCid]->GetYaxis()->SetTitle("Particle Fraction");
      hMCgenPrimFractionXi[MCid]->Divide(hMCgenPrimYieldXi[MCid], hMCgenPrimYieldTotalXi, 1., 1., "B");
    }
  }
  
  // Create empty histos vs. z and xi
  TH1D* hFractionIDFFz[AliPID::kSPECIES] = { 0x0, };
  TH1D* hFractionIDFFxi[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hIDFFz[AliPID::kSPECIES] = { 0x0, };
  TH1D* hIDFFxi[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hRatioToPiIDFFz[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiIDFFxi[AliPID::kSPECIES] = { 0x0, };
  
  
  TH1D* hFractionIDFFzMC[AliPID::kSPECIES] = { 0x0, };
  TH1D* hFractionIDFFxiMC[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hIDFFzMC[AliPID::kSPECIES] = { 0x0, };
  TH1D* hIDFFxiMC[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hRatioToPiIDFFzMC[AliPID::kSPECIES] = { 0x0, };
  TH1D* hRatioToPiIDFFxiMC[AliPID::kSPECIES] = { 0x0, };
  
  TH1D* hMCgenPrimRatioToPiPt[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimRatioToPiZ[AliPID::kSPECIES] = { 0x0, };
  TH1D* hMCgenPrimRatioToPiXi[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hIDFFz[i] = hDataProj->Projection(kPidZProj);
    hIDFFz[i]->Reset();
    if (hIDFFz[i]->GetSumw2N() <= 0)
      hIDFFz[i]->Sumw2();
    setupHist(hIDFFz[i], Form("hIDFFz_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
              hDataProj->GetAxis(kPidZProj)->GetTitle(), yieldAxisTitleZ.Data(), getLineColorAliPID(i), kFALSE);
    hFractionIDFFz[i] = new TH1D(*hIDFFz[i]);
    hFractionIDFFz[i]->SetName(Form("hFractionIDFFz_%s", AliPID::ParticleShortName(i)));
    hFractionIDFFz[i]->GetYaxis()->SetTitle("Particle Fraction");
    
    hIDFFxi[i] = hDataProj->Projection(kPidXiProj);
    hIDFFxi[i]->Reset();
    if (hIDFFxi[i]->GetSumw2N() <= 0)
      hIDFFxi[i]->Sumw2();    
    setupHist(hIDFFxi[i], Form("hIDFFxi_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
              hDataProj->GetAxis(kPidXiProj)->GetTitle(), yieldAxisTitleXi.Data(), getLineColorAliPID(i), kFALSE);
    hFractionIDFFxi[i] = new TH1D(*hIDFFxi[i]);
    hFractionIDFFxi[i]->SetName(Form("hFractionIDFFxi_%s", AliPID::ParticleShortName(i)));
    hFractionIDFFxi[i]->GetYaxis()->SetTitle("Particle Fraction");
    
    // If MC histos for pT are available, also create those vs. z,xi
    if (hFractionIDFFtrackPtMC[i]) {
      hIDFFzMC[i] = new TH1D(*hIDFFz[i]);
      setupHist(hIDFFzMC[i], Form("hIDFFzMC_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
                "", "", getLineColorAliPID(i), kTRUE);
      hFractionIDFFzMC[i] = new TH1D(*hIDFFzMC[i]);
      hFractionIDFFzMC[i]->SetName(Form("hFractionIDFFzMC_%s", AliPID::ParticleShortName(i)));
      hFractionIDFFzMC[i]->GetYaxis()->SetTitle("Particle Fraction");
      
      hIDFFxiMC[i] = new TH1D(*hIDFFxi[i]);
      setupHist(hIDFFxiMC[i], Form("hIDFFxiMC_%s", AliPID::ParticleShortName(i)), Form("%s", AliPID::ParticleShortName(i)),
               "", "", getLineColorAliPID(i), kTRUE);
      hFractionIDFFxiMC[i] = new TH1D(*hIDFFxiMC[i]);
      hFractionIDFFxiMC[i]->SetName(Form("hFractionIDFFxiMC_%s", AliPID::ParticleShortName(i)));
      hFractionIDFFxiMC[i]->GetYaxis()->SetTitle("Particle Fraction");
    }
    
    
    
    // If to-pi-ratios are available for pT, also create those for z, xi.
    // Additionally, calculate the generated true to-pi-ratios.
    if (hRatioToPiIDFFtrackPt[i]) {
      hRatioToPiIDFFz[i] = new TH1D(*hIDFFz[i]);
      hRatioToPiIDFFz[i]->SetName(Form("hRatioToPiIDFFz_%s", AliPID::ParticleShortName(i)));
      hRatioToPiIDFFz[i]->SetTitle(hRatioToPiIDFFtrackPt[i]->GetTitle());
      TString tempTitle = hRatioToPiIDFFtrackPt[i]->GetYaxis()->GetTitle();
      tempTitle.ReplaceAll("P_{T}", "z");
      tempTitle.ReplaceAll("p_{T}", "z");
      hRatioToPiIDFFz[i]->GetYaxis()->SetTitle(tempTitle.Data());
      
      hRatioToPiIDFFxi[i] = new TH1D(*hIDFFxi[i]);
      hRatioToPiIDFFxi[i]->SetName(Form("hRatioToPiIDFFxi_%s", AliPID::ParticleShortName(i)));
      hRatioToPiIDFFxi[i]->SetTitle(hRatioToPiIDFFtrackPt[i]->GetTitle());
      tempTitle = hRatioToPiIDFFtrackPt[i]->GetYaxis()->GetTitle();
      tempTitle.ReplaceAll("P_{T}", "#xi");
      tempTitle.ReplaceAll("p_{T}", "#xi");
      hRatioToPiIDFFxi[i]->GetYaxis()->SetTitle(tempTitle.Data());
      
      // Also for MC, if available
      if (hRatioToPiIDFFtrackPtMC[i]) {
        hRatioToPiIDFFzMC[i] = new TH1D(*hRatioToPiIDFFz[i]);
        setupHist(hRatioToPiIDFFzMC[i], Form("hRatioToPiIDFFzMC_%s", AliPID::ParticleShortName(i)), hRatioToPiIDFFz[i]->GetTitle(),
                  "", "", getLineColorAliPID(i), kTRUE);
        
        hRatioToPiIDFFxiMC[i] = new TH1D(*hRatioToPiIDFFxi[i]);
        setupHist(hRatioToPiIDFFxiMC[i], Form("hRatioToPiIDFFxiMC_%s", AliPID::ParticleShortName(i)), hRatioToPiIDFFxi[i]->GetTitle(),
                  "", "", getLineColorAliPID(i), kTRUE);
      }
      
      // Generated true to-pi-ratios
      if (hMCgenPrimYieldPt[i]) {
        hMCgenPrimRatioToPiPt[i] = new TH1D(*hMCgenPrimYieldPt[i]);
        setupHist(hMCgenPrimRatioToPiPt[i],
                  Form("hMCgenPrimRatioToPiPt_%s", AliPID::ParticleShortName(i)), hRatioToPiIDFFtrackPt[i]->GetTitle(),  "",
                  hRatioToPiIDFFtrackPt[i]->GetYaxis()->GetTitle(), getLineColorAliPID(i), kTRUE);
        hMCgenPrimRatioToPiPt[i]->SetMarkerStyle(hMCgenPrimYieldPt[i]->GetMarkerStyle());
        
        // True yield -> Statistically independent, just divide
        hMCgenPrimRatioToPiPt[i]->Divide(hMCgenPrimYieldPt[i], hMCgenPrimYieldPt[AliPID::kPion]);
      }
      
      if (hMCgenPrimYieldZ[i]) {
        hMCgenPrimRatioToPiZ[i] = new TH1D(*hMCgenPrimYieldZ[i]);
        setupHist(hMCgenPrimRatioToPiZ[i],
                  Form("hMCgenPrimRatioToPiZ_%s", AliPID::ParticleShortName(i)), hRatioToPiIDFFz[i]->GetTitle(),  "",
                  hRatioToPiIDFFz[i]->GetYaxis()->GetTitle(), getLineColorAliPID(i), kTRUE);
        hMCgenPrimRatioToPiZ[i]->SetMarkerStyle(hMCgenPrimYieldZ[i]->GetMarkerStyle());
        
        // True yield -> Statistically independent, just divide
        hMCgenPrimRatioToPiZ[i]->Divide(hMCgenPrimYieldZ[i], hMCgenPrimYieldZ[AliPID::kPion]);
      }
      
      if (hMCgenPrimYieldXi[i]) {
        hMCgenPrimRatioToPiXi[i] = new TH1D(*hMCgenPrimYieldXi[i]);
        setupHist(hMCgenPrimRatioToPiXi[i],
                  Form("hMCgenPrimRatioToPiXi_%s", AliPID::ParticleShortName(i)), hRatioToPiIDFFxi[i]->GetTitle(),  "",
                  hRatioToPiIDFFxi[i]->GetYaxis()->GetTitle(), getLineColorAliPID(i), kTRUE);
        hMCgenPrimRatioToPiXi[i]->SetMarkerStyle(hMCgenPrimYieldXi[i]->GetMarkerStyle());
        
        // True yield -> Statistically independent, just divide
        hMCgenPrimRatioToPiXi[i]->Divide(hMCgenPrimYieldXi[i], hMCgenPrimYieldXi[AliPID::kPion]);
      }
    }
  }
  
  // Obtain maps for change of variable (pT -> z,xi) which in the end contain the weighting factors
  TH2D* hWeightingXivsPt = hDataProj->Projection(kPidXiProj, kPidPtProj, "e");
  hWeightingXivsPt->SetName("hWeightingXivsPt");
  hWeightingXivsPt->SetTitle("Weighting factors for p_{T} -> #xi");
  
  TH2D* hErrorWeightingXivsPt = new TH2D(*hWeightingXivsPt);
  hErrorWeightingXivsPt->SetName("hErrorWeightingXivsPt");
  hErrorWeightingXivsPt->SetTitle("Error weighting factors for p_{T} -> #xi");
  
  TH2D* hWeightingZvsPt = hDataProj->Projection(kPidZProj, kPidPtProj, "e");
  hWeightingZvsPt->SetName("hWeightingZvsPt");
  hWeightingZvsPt->SetTitle("Weighting factors for p_{T} -> z");
  
  TH2D* hErrorWeightingZvsPt = new TH2D(*hWeightingZvsPt);
  hErrorWeightingZvsPt->SetName("hErrorWeightingZvsPt");
  hErrorWeightingZvsPt->SetTitle("Error weighting factors for p_{T} -> z");
  
  // Get the total weights in each row (i.e. fixed z or xi).
  // NOTE: Setting the projection from 1 to nBinsX already EXCLUDES underflow and overflow bins!
  // Anyway this should make no difference since one anyway demands pT > 0.15 GeV/c (or this is even included in the cuts)
  // and pT > 50 GeV/c has more or less no statistics (only relevant for jetPt > 50 GeV/c, which means almost zero jets).
  TH1D* hTotalWeightXi = hWeightingXivsPt->ProjectionY("hTotalWeightXi", 1, hWeightingXivsPt->GetNbinsX(), "e");
  hTotalWeightXi->SetTitle("Total weight for p_{T} -> #xi");
  
  TH1D* hTotalWeightZ  = hWeightingZvsPt->ProjectionY("hTotalWeightZ", 1, hWeightingZvsPt->GetNbinsX(), "e"); 
  hTotalWeightZ->SetTitle("Total weight for p_{T} -> z");
  
  // Normalise the weighting factors such that every row sums up to 1
  normaliseWeightingFactor(hWeightingXivsPt, hTotalWeightXi);
  normaliseWeightingFactor(hWeightingZvsPt, hTotalWeightZ);
  
  
  // Prepare save of results to file
  TString chargeString = "";
  if (chargeMode == kPosCharge)
    chargeString = "_posCharge";
  else if (chargeMode == kNegCharge)
    chargeString = "_negCharge";
  
  TString saveFileName = pathNameFractionsAndYields;
  saveFileName.Replace(0, pathNameFractionsAndYields.Last('/') + 1, "");
  
  TString savePath = pathNameFractionsAndYields;
  savePath.ReplaceAll(Form("/%s", saveFileName.Data()), "");
  
  saveFileName.Prepend("output_extractedFFs_");
  saveFileName.ReplaceAll(".root", Form("__centrality_%s%s%s.root",
                                        restrictCentralityAxis ? Form("%.0f_%.0f.root", actualLowerCentrality, actualUpperCentrality)
                                                               : "all",
                                        restrictJetPtAxis ? Form("_jetPt%.1f_%.1f", actualLowerJetPt, actualUpperJetPt) : "",
                                        chargeString.Data()));
  
  TString saveFilePathName = Form("%s/%s", savePath.Data(), saveFileName.Data());
  TFile* saveFile = TFile::Open(saveFilePathName.Data(), "RECREATE");
  
  if (!saveFile) {
    printf("Failed to save results to file \"%s\"!\n", saveFilePathName.Data());
    return -1;
  }
  
  
  
  // Get the total weights in each column (i.e. fixed pT).
  // NOTE: Just take the projection of z on pT. Is the same as for xi.
  TH1D* hTotalErrorWeightPt = hErrorWeightingZvsPt->ProjectionX("hTotalErrorWeightPt", 1, hErrorWeightingZvsPt->GetNbinsY(), "e");
  hTotalErrorWeightPt->SetTitle("Total error weight for p_{T} -> z,#xi");
  
  // Normalise the weighting factors such that every column sums up to 1
  normaliseErrorWeightingFactor(hErrorWeightingXivsPt, hTotalErrorWeightPt);
  normaliseErrorWeightingFactor(hErrorWeightingZvsPt, hTotalErrorWeightPt);
  
  
  // Calculate the weighted means for the fractions
  calculateWeightedMean(&hFractionIDFFz[0],  &hFractionIDFFtrackPt[0], hWeightingZvsPt,  hErrorWeightingZvsPt);
  calculateWeightedMean(&hFractionIDFFxi[0], &hFractionIDFFtrackPt[0], hWeightingXivsPt, hErrorWeightingXivsPt);
  
  // Obtain the yields
  translateFractionToYield(&hIDFFz[0], &hFractionIDFFz[0], hTotalWeightZ);
  translateFractionToYield(&hIDFFxi[0], &hFractionIDFFxi[0], hTotalWeightXi);
  
  // Calculate the weighted means for the to-pi-ratios
  calculateWeightedMean(&hRatioToPiIDFFz[0],  &hRatioToPiIDFFtrackPt[0], hWeightingZvsPt,  hErrorWeightingZvsPt);
  calculateWeightedMean(&hRatioToPiIDFFxi[0], &hRatioToPiIDFFtrackPt[0], hWeightingXivsPt, hErrorWeightingXivsPt);
  
  
  // Same for MC, if available. The weights are the same, only the fraction and fraction errors change (MC ID used)
  if (hFractionIDFFtrackPtMC[AliPID::kPion]) {
    // Calculate the weighted means for the fractions
    calculateWeightedMean(&hFractionIDFFzMC[0],  &hFractionIDFFtrackPtMC[0], hWeightingZvsPt,  hErrorWeightingZvsPt);
    calculateWeightedMean(&hFractionIDFFxiMC[0], &hFractionIDFFtrackPtMC[0], hWeightingXivsPt, hErrorWeightingXivsPt);
    
    // Obtain the yields
    translateFractionToYield(&hIDFFzMC[0], &hFractionIDFFzMC[0], hTotalWeightZ);
    translateFractionToYield(&hIDFFxiMC[0], &hFractionIDFFxiMC[0], hTotalWeightXi);
    
    // Calculate the weighted means for the to-pi-ratios
    calculateWeightedMean(&hRatioToPiIDFFzMC[0],  &hRatioToPiIDFFtrackPtMC[0], hWeightingZvsPt,  hErrorWeightingZvsPt);
    calculateWeightedMean(&hRatioToPiIDFFxiMC[0], &hRatioToPiIDFFtrackPtMC[0], hWeightingXivsPt, hErrorWeightingXivsPt);
  }
  
  
  
  
  /*
  // First iteration: Just take the default fractions to obtain the "mean" of the yields
  Bool_t setMean = kTRUE;
  Bool_t addErrorsQuadratically = kFALSE;
  Bool_t smearByError = kFALSE;
  Bool_t takeIntoAccountSysError = kFALSE;
  // setMean and all the follwoing parameters are anyway irrelevant, if smearByError = kFALSE and takeIntoAccountSysError = kFALSE
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPt, hIDFFz, hIDFFxi, setMean, addErrorsQuadratically,
                         smearByError, uniformSystematicError, takeIntoAccountSysError, 1);
  

  // Next iteration: Vary the PID map within statistical errors several times and calculate the resulting statistical error
  // of the spectra from this. But since the mean of the fraction is some kind of "best estimate" of the truth, leave the means
  // as they have been set during the last step.
  
  
  // NOTE: Do NOT add the statistical errors from the last step, since the yields/fractions (rel. error is the same) already
  // contain the stat. uncertainty of the yield of this species in the corresponding bin! I.e. one just takes some kind of weighted
  // mean of the fractions (mean and error) (equivalent of taking the yields with corresponding error).
  setMean = kFALSE;
  addErrorsQuadratically = kFALSE;
  smearByError = kTRUE;
  takeIntoAccountSysError = kFALSE;
  const Int_t nGenerations = 5000; //TODO set to 5000 in the end, after all testing was successful
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPt, hIDFFz, hIDFFxi, setMean, addErrorsQuadratically,
                         smearByError, uniformSystematicError, takeIntoAccountSysError, nGenerations);
  
  
  // Next iteration: Vary the PID map within systematic errors several times and calculate the resulting systematic error
  // of the spectra from this. But since the mean of the fraction is some kind of "best estimate" of the truth, leave the means
  // as they have been set during the last step.
  setMean = kFALSE;
  addErrorsQuadratically = kFALSE;
  smearByError = kFALSE;
  takeIntoAccountSysError = kTRUE;
  // Clone histograms with final statistical errors -> Only set systematic errors, but leave mean as it is
  TH2D* hIDFFtrackPtSysError[AliPID::kSPECIES] = { 0x0, };
  TH2D* hIDFFzSysError[AliPID::kSPECIES] = { 0x0, };
  TH2D* hIDFFxiSysError[AliPID::kSPECIES] = { 0x0, };
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    hIDFFtrackPtSysError[i] = new TH2D(*hIDFFtrackPt[i]);
    hIDFFtrackPtSysError[i]->SetName(Form("%s_sysError", hIDFFtrackPt[i]->GetName()));
    hIDFFtrackPtSysError[i]->SetFillStyle(0);
    
    hIDFFzSysError[i] = new TH2D(*hIDFFz[i]);
    hIDFFzSysError[i]->SetName(Form("%s_sysError", hIDFFz[i]->GetName()));
    hIDFFzSysError[i]->SetFillStyle(0);
    
    hIDFFxiSysError[i] = new TH2D(*hIDFFxi[i]);
    hIDFFxiSysError[i]->SetName(Form("%s_sysError", hIDFFxi[i]->GetName()));
    hIDFFxiSysError[i]->SetFillStyle(0);
  }
  
  GenerateParticleYields(hDataProj, pidTask, centrality, hIDFFtrackPtSysError, hIDFFzSysError, hIDFFxiSysError, setMean,
                         addErrorsQuadratically, smearByError, uniformSystematicError, takeIntoAccountSysError, nGenerations);
  
  delete pidTask;
 */
  
  // Normalise properly
  const Double_t numJetsRec = hNjetsRec ? hNjetsRec->Integral(lowerCentralityBinLimit, upperCentralityBinLimit,
                                                              lowerJetPtBinLimit, upperJetPtBinLimit) : 0.;
  const Double_t numJetsGen = hNjetsGen ? hNjetsGen->Integral(lowerCentralityBinLimit, upperCentralityBinLimit,
                                                              lowerJetPtBinLimit, upperJetPtBinLimit) : 0.;
  Bool_t noJets = numJetsRec < 1e-13;
  
  if (noJets) {
    printf("Error: No jets in desired range!\n");
    return -1;
  }
  
  for (Int_t species = 0; species < AliPID::kSPECIES; species++) {
    // NOTE: Pt is already properly normalised (except for generated true yield!)
    
    //normaliseYieldHist(hIDFFtrackPt[species], numJetsRec);
    normaliseYieldHist(hIDFFz[species], numJetsRec);
    normaliseYieldHist(hIDFFxi[species], numJetsRec);
    
    // This is RECONSTRUCTED, but identified via MC label -> Normalise to number of RECONSTRUCTED jets!
    //normaliseYieldHist(hIDFFtrackPtMC[species], numJetsRec);
    normaliseYieldHist(hIDFFzMC[species], numJetsRec);
    normaliseYieldHist(hIDFFxiMC[species], numJetsRec);
    
    // GENERATED yield must be normalised to number of GENERATED jets!
    normaliseYieldHist(hMCgenPrimYieldPt[species], numJetsGen);
    normaliseYieldHist(hMCgenPrimYieldZ[species], numJetsGen);
    normaliseYieldHist(hMCgenPrimYieldXi[species], numJetsGen);
  }
  
  normaliseYieldHist(hMCgenPrimYieldTotalPt, numJetsGen);
  normaliseYieldHist(hMCgenPrimYieldTotalZ, numJetsGen);
  normaliseYieldHist(hMCgenPrimYieldTotalXi, numJetsGen);
  
  
  saveFile->cd();
  
  if (hWeightingXivsPt)
    hWeightingXivsPt->Write();
  
  if (hTotalWeightXi)
    hTotalWeightXi->Write();
  
  if (hErrorWeightingXivsPt)
    hErrorWeightingXivsPt->Write();
  

  if (hWeightingZvsPt)
    hWeightingZvsPt->Write();
  
  if (hTotalWeightZ)
    hTotalWeightZ->Write();
  
  if (hErrorWeightingZvsPt)
    hErrorWeightingZvsPt->Write();
  
  
  if (hTotalErrorWeightPt)
    hTotalErrorWeightPt->Write();
  
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    if (hIDFFtrackPt[i])
      hIDFFtrackPt[i]->Write();
    
    if (hIDFFz[i])
      hIDFFz[i]->Write();
    
    if (hIDFFxi[i])
      hIDFFxi[i]->Write();
    
    
    if (hFractionIDFFtrackPt[i])
      hFractionIDFFtrackPt[i]->Write();
    
    if (hFractionIDFFz[i])
      hFractionIDFFz[i]->Write();
    
    if (hFractionIDFFxi[i])
      hFractionIDFFxi[i]->Write();
    
    
    if (hRatioToPiIDFFtrackPt[i])
      hRatioToPiIDFFtrackPt[i]->Write();
    
    if (hRatioToPiIDFFz[i])
      hRatioToPiIDFFz[i]->Write();
    
    if (hRatioToPiIDFFxi[i])
      hRatioToPiIDFFxi[i]->Write();
    
    
    if (hIDFFtrackPtMC[i])
      hIDFFtrackPtMC[i]->Write();
    
    if (hIDFFzMC[i])
      hIDFFzMC[i]->Write();
    
    if (hIDFFxiMC[i])
      hIDFFxiMC[i]->Write();
    
    
    if (hFractionIDFFtrackPtMC[i])
      hFractionIDFFtrackPtMC[i]->Write();
    
    if (hFractionIDFFzMC[i])
      hFractionIDFFzMC[i]->Write();
    
    if (hFractionIDFFxiMC[i])
      hFractionIDFFxiMC[i]->Write();
    
    
    if (hRatioToPiIDFFtrackPtMC[i])
      hRatioToPiIDFFtrackPtMC[i]->Write();
    
    if (hRatioToPiIDFFzMC[i])
      hRatioToPiIDFFzMC[i]->Write();
    
    if (hRatioToPiIDFFxiMC[i])
      hRatioToPiIDFFxiMC[i]->Write();
    
    
    if (hMCgenPrimYieldPt[i])
      hMCgenPrimYieldPt[i]->Write();
    
    if (hMCgenPrimYieldZ[i])
      hMCgenPrimYieldZ[i]->Write();
    
    if (hMCgenPrimYieldXi[i])
      hMCgenPrimYieldXi[i]->Write();
    
    
    if (hMCgenPrimFractionPt[i])
      hMCgenPrimFractionPt[i]->Write();
    
    if (hMCgenPrimFractionZ[i])
      hMCgenPrimFractionZ[i]->Write();
    
    if (hMCgenPrimFractionXi[i])
      hMCgenPrimFractionXi[i]->Write();
    
    
    if (hMCgenPrimRatioToPiPt[i])
      hMCgenPrimRatioToPiPt[i]->Write();
    
    if (hMCgenPrimRatioToPiZ[i])
      hMCgenPrimRatioToPiZ[i]->Write();
    
    if (hMCgenPrimRatioToPiXi[i])
      hMCgenPrimRatioToPiXi[i]->Write();
  }
  
  
  if (hMCgenPrimYieldTotalPt)
      hMCgenPrimYieldTotalPt->Write();
    
  if (hMCgenPrimYieldTotalZ)
    hMCgenPrimYieldTotalZ->Write();
  
  if (hMCgenPrimYieldTotalXi)
    hMCgenPrimYieldTotalXi->Write();
  
  
  if (hNjetsGen)
    hNjetsGen->Write();
  
  if (hNjetsRec)
    hNjetsRec->Write();
  

  TNamed* settings = new TNamed(
      Form("Settings: Data file \"%s\", File with fractions and yields \"%s\", lowerCentrality %.0f, upperCentrality %.0f, chargeMode %d, rebinZ %d, rebinXi %d, onlyUseRelevantMCIDforMatrix %d\n",
           pathNameData.Data(), pathNameFractionsAndYields.Data(), lowerCentrality, upperCentrality, chargeMode, rebinZ, rebinXi, 
           onlyUseRelevantMCIDforMatrix), "");
  settings->Write();
  
  saveFile->Close();
  
  
  
  /*
  Int_t t1 = hIDFFtrackPt[AliPID::kPion]->GetYaxis()->FindBin(5.1);
  Int_t t2 = hIDFFtrackPt[AliPID::kPion]->GetYaxis()->FindBin(9.9);
  
  printf("t1 %d, t2 %d\n", t1, t2);
  new TCanvas();
  TH1D* hTemp = hIDFFtrackPt[AliPID::kPion]->ProjectionX("_pfy1", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFtrackPtSysError[AliPID::kPion]->ProjectionX("sys_pfy1", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldPt[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldPt[AliPID::kPion]->ProjectionX("MC_pfy1", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFz[AliPID::kPion]->ProjectionX("_pfy2", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFzSysError[AliPID::kPion]->ProjectionX("sys_pfy2", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldZ[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldZ[AliPID::kPion]->ProjectionX("MC_pfy2", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFxi[AliPID::kPion]->ProjectionX("_pfy3", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFxiSysError[AliPID::kPion]->ProjectionX("sys_pfy3", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldXi[AliPID::kPion]) {
    hTemp = hMCgenPrimYieldXi[AliPID::kPion]->ProjectionX("MC_pfy3", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  
  
  new TCanvas();
  hTemp = hIDFFtrackPt[AliPID::kProton]->ProjectionX("_pfy4", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFtrackPtSysError[AliPID::kProton]->ProjectionX("sys_pfy4", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldPt[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldPt[AliPID::kProton]->ProjectionX("MC_pfy4", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFz[AliPID::kProton]->ProjectionX("_pfy5", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFzSysError[AliPID::kProton]->ProjectionX("sys_pfy5", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldZ[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldZ[AliPID::kProton]->ProjectionX("MC_pfy5", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }
  
  
  new TCanvas();
  hTemp = hIDFFxi[AliPID::kProton]->ProjectionX("_pfy", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("");
  
  hTemp = hIDFFxiSysError[AliPID::kProton]->ProjectionX("sys_pfy", t1, t2, "e");
  hTemp->SetFillStyle(0);
  hTemp->Draw("E2same");
  
  if (hMCgenPrimYieldXi[AliPID::kProton]) {
    hTemp = hMCgenPrimYieldXi[AliPID::kProton]->ProjectionX("MC_pfy", t1, t2, "e");
    hTemp->SetFillStyle(0);
    hTemp->SetMarkerStyle(23);
    hTemp->Draw("same");
  }*/
  
  f->Close();
  fFractionsAndYields->Close();
  
  return 0;
}