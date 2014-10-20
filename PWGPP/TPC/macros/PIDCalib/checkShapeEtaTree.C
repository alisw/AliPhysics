#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH3D.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVirtualFitter.h"
#include "TLinearFitter.h"
#include "TList.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMath.h"
#include "TDatime.h"
#include "TSystem.h"
#include "TSpline.h"

#include "AliPID.h"

#include <iostream>
#include <iomanip>

#include "THnSparseDefinitions.h"
#include "ProgressBar.h"

enum mapType {
  kNoNormalisation = 1,
  kSmallEtaNormalisation = 2,
  kLargeEtaNormalisation = 3
};
enum type {
  kTypeDelta = 1,
  kTypeDeltaPrime = 2,
  kTypeSigmaDeltaPrime = 3
};

const Double_t binContentThreshold = 1e-4;
const Double_t massProton = AliPID::ParticleMass(AliPID::kProton);
const Double_t massPion = AliPID::ParticleMass(AliPID::kPion);
const TString suffixMapType[kLargeEtaNormalisation + 1] = {"", "NoNormalisation", "SmallEtaNormalisation", "LargeEtaNormalisation"};

TF1 fGauss("fGauss", "[0]*TMath::Gaus(x, [1], [2], 1)", 0.6, 1.6);


//TODO NOW getBin... functions needed? Only, if correction (momentum) necessary
//_________________________________________________________________________________________
Int_t getBinX(TH2D* h, Double_t tanTheta)
{
  Int_t binX = h->GetXaxis()->FindBin(tanTheta);
  
  if (binX <= 0)
    binX = 1;
  if (binX > h->GetXaxis()->GetNbins())
    binX = h->GetXaxis()->GetNbins();

    return binX;
}


//_________________________________________________________________________________________
Int_t getBinY(TH2D* h, Double_t dEdxInv)
{
  Int_t binY = h->GetYaxis()->FindBin(dEdxInv);
  
  if (binY <= 0)
    binY = 1;
  if (binY > h->GetYaxis()->GetNbins())
    binY = h->GetYaxis()->GetNbins();

  return binY;
}


//__________________________________________________________________________________________
Bool_t FindFitRange(TH1* h, Double_t& lowThreshold, Double_t& highThreshold, Double_t fractionForRange = 0.1)
{
    lowThreshold = 0.6;
    highThreshold = 1.6;
    
    if (!h)
      return kFALSE;
    
    // Average around maximum bin -> Might catch some outliers
    Int_t maxBin = h->GetMaximumBin();
    Double_t maxVal = h->GetBinContent(maxBin);
    UChar_t usedBins = 1;
    if (maxBin > 1) {
      maxVal += h->GetBinContent(maxBin - 1);
      usedBins++;
    }
    if (maxBin < h->GetNbinsX()) {
      maxVal += h->GetBinContent(maxBin + 1);
      usedBins++;
    }
    maxVal /= usedBins;
    
    Double_t thresholdFraction = fractionForRange * maxVal; 
    Int_t lowThrBin = TMath::Max(1, h->FindFirstBinAbove(thresholdFraction));
    Int_t highThrBin = TMath::Min(h->GetNbinsX(), h->FindLastBinAbove(thresholdFraction));
    
    lowThreshold = h->GetBinCenter(lowThrBin);
    highThreshold = h->GetBinCenter(highThrBin);
    
    return kTRUE;
}


//__________________________________________________________________________________________
void SetHistAxis(TH2D* h, Int_t type)
{
  h->GetXaxis()->SetTitle("tan(#Theta)");
  h->GetYaxis()->SetTitle("1/(dE/dx) (arb. unit)");
  if (type == kTypeDelta)
    h->GetZaxis()->SetTitle("#Delta (arb. unit)");
  else if (type == kTypeDeltaPrime) 
    h->GetZaxis()->SetTitle("#Delta' (arb. unit)");
  else if (type == kTypeSigmaDeltaPrime) 
    h->GetZaxis()->SetTitle("Par1(#sigma_{rel}(#Delta'))");
  else
    printf("SetHistAxis: Unknown type %d!\n", type);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(2.2);
  h->GetXaxis()->SetLabelSize(0.03);
  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(2.6);
  h->GetYaxis()->SetLabelSize(0.03);
  h->GetZaxis()->SetTitleSize(0.04);
  h->GetZaxis()->SetTitleOffset(1.3);
}

//__________________________________________________________________________________________
Double_t getMedianOfNonZeros(Double_t* input, const Int_t dim)
{
  Double_t values[dim];
  for (Int_t i = 0; i < dim; i++)
    values[i] = 0.0;
    
  Int_t numNonZero = 0;
  
  for (Int_t i = 0; i < dim; i++) {
    if (input[i] > 0) {
      values[numNonZero] = input[i];
      numNonZero++;
    }
  }

  return ((numNonZero > 0) ? TMath::Median(numNonZero, values) : 0);
}


//__________________________________________________________________________________________
void normaliseHisto(TH2D* h, Double_t lowerLimit, Double_t upperLimit, Int_t type)
{
  if (lowerLimit >= upperLimit) {
    printf("!!!Error normaliseHisto: lowerLimit >= upperLimit!\n");
    return;
  }
  
  if (lowerLimit < 0 || upperLimit < 0) {
    printf("!!!Error normaliseHisto: lowerLimit, upperLimit >= 0 required!\n");
    return;
  }
  
  Int_t binLow = h->GetXaxis()->FindBin(lowerLimit);
  if (binLow < 1)
    binLow = 1;
  
  Int_t binHigh = h->GetXaxis()->FindBin(upperLimit);
  if (binHigh > h->GetXaxis()->GetNbins())
    binHigh = h->GetXaxis()->GetNbins();
  
  Int_t binLow2 = h->GetXaxis()->FindBin(-upperLimit);
  if (binLow2 < 1)
    binLow2 = 1;
  
  Int_t binHigh2 = h->GetXaxis()->FindBin(-lowerLimit);
  if (binHigh2 > h->GetXaxis()->GetNbins())
    binHigh2 = h->GetXaxis()->GetNbins();
  
  // If 0 is included in range, it might happen that both ranges overlap -> Remove one of the double counted binsPforMap
    if (binHigh2 >= binLow) {
      binHigh2 = binLow - 1;
      if (binHigh2 < 1)   {
        printf("!!!Error: binHigh2 = binLow - 1 is < 1\n");
        return;
      }
    }
    
    if (binLow2 > binHigh2)
      binLow2 = binHigh2;
    
    // Normalise with respect to some given tan(theta)
      // To be robust against fluctuations: Take median of a few tan(theta) bins around the desired value for scaling
      const Int_t nThetaBins = (binHigh - binLow + 1) + (binHigh2 - binLow2 + 1);
      
      if (nThetaBins <= 0)  {
        printf("Cannot renormalise histo due to bad limits for reference bins: %f -> %f\n", lowerLimit, upperLimit);
        return;
      }
      
      Double_t* values;
      values = new Double_t[nThetaBins];
      
      for (Int_t i = 1; i <= h->GetYaxis()->GetNbins(); i++) {
        // Reset values
        Int_t binIndex = 0;
        for (Int_t thetaBin = 1; thetaBin <= nThetaBins; thetaBin++)
          values[thetaBin - 1] = 0;
        
        for (Int_t thetaBin = binLow; thetaBin <= binHigh; thetaBin++)  {
          Double_t temp = h->GetBinContent(thetaBin, i);
          values[binIndex] = (temp > 0 || type == kTypeDelta) ? temp : 0;
          binIndex++;
        }
        for (Int_t thetaBin = binLow2; thetaBin <= binHigh2; thetaBin++)  {
          Double_t temp = h->GetBinContent(thetaBin, i);
          values[binIndex] = (temp > 0 || type == kTypeDelta) ? temp : 0;
          binIndex++;
        }
        
        Double_t temp = 0;
        Double_t scale = 0;
        
        if (type == kTypeDeltaPrime)  {
          temp = getMedianOfNonZeros(values, nThetaBins);
          if (temp <= 0)
            continue;
          scale = 1. / temp;
        }
        else if (type == kTypeDelta)  {
          scale = TMath::Median(nThetaBins, values);
        }
        else  {
          printf("Error: Type %d not supported for normalisation!\n", type);
          return;
        }
        
        
        for (Int_t thetaBin = 1; thetaBin <= h->GetXaxis()->GetNbins(); thetaBin++)  {
          if (type == kTypeDeltaPrime)  {
            h->SetBinContent(thetaBin, i, h->GetBinContent(thetaBin, i) * scale);
            h->SetBinError(thetaBin, i, h->GetBinError(thetaBin, i) * scale);
          }
          else if (type == kTypeDelta)  {
            h->SetBinContent(thetaBin, i, h->GetBinContent(thetaBin, i) - scale);
            h->SetBinError(thetaBin, i, h->GetBinError(thetaBin, i) - scale);
          }
        }
      }  
      delete values;
}


//__________________________________________________________________________________________
void eliminateNonPositivePointsInHisto(TH2D* h)
{
  const Int_t nNeighbours = 24;
  Double_t* values;
  values = new Double_t[nNeighbours];

  // Search for bins with content <= 0 (bins with fit errors). Then take all surrounding points in +-2 rows and columns
  // and take their median (only those bins without fit errors!) as the new bin content of the considered bin.
  
  for (Int_t binX = 1; binX <= h->GetXaxis()->GetNbins(); binX++) {
    Int_t firstBinLeft = TMath::Max(1, binX - 2);
    Int_t lastBinRight = TMath::Min(h->GetXaxis()->GetNbins(), binX + 2);
    
    for (Int_t binY = 1; binY <= h->GetYaxis()->GetNbins(); binY++) {
      if (h->GetBinContent(binX, binY) <= 0) {
        Int_t firstBinBelow = TMath::Max(1, binY - 2);
        Int_t lastBinAbove = TMath::Min(h->GetYaxis()->GetNbins(), binY + 2);
        
        // Reset values
        Int_t binIndex = 0;
        for (Int_t i = 0; i < nNeighbours; i++)
          values[i] = 0;
        
        for (Int_t binX2 = firstBinLeft; binX2 <= lastBinRight; binX2++) {
          for (Int_t binY2 = firstBinBelow; binY2 <= lastBinAbove; binY2++) {
            if (binX2 == binX && binY2 == binY)
              continue; // skip bin that is currently under consideration
            
            // Only take values from (hopefully) valid fits
            if (h->GetBinContent(binX2, binY2) > 0) {
              values[binIndex] = h->GetBinContent(binX2, binY2);
              binIndex++;
            }
          }
        }
        
        Double_t temp = getMedianOfNonZeros(values, nNeighbours);
        if (temp <= 0) {
          // Only print error message, if there is at least one positive value
          if (binIndex > 0)
            printf("Error: Could not eliminate values <= 0 for bin at (%f, %f)!\n",
                  h->GetXaxis()->GetBinCenter(binX), h->GetYaxis()->GetBinCenter(binY));
          temp = -1;
        }
        
        printf("Eliminated non-positive value at bin (%f, %f): %f -> %f!\n",
               h->GetXaxis()->GetBinCenter(binX), h->GetYaxis()->GetBinCenter(binY),
               h->GetBinContent(binX, binY), temp);
        h->SetBinContent(binX, binY, temp);
        h->SetBinError(binX, binY, 1000);
      }
    }
  }
  
  delete values;
}


//__________________________________________________________________________________________
void eliminateOutliers(TH2D* h, Double_t threshold)
{
  const Int_t nNeighbours = 8;
  Double_t* values;
  values = new Double_t[nNeighbours];

  // Search for outliers (most likely bad fits). Then take all surrounding points in +-1 rows and columns
  // and take their median (only those bins without (known) fit errors, i.e. bin content > 0).
  // If the current bin content deviates by more than "threshold" from the median, take the median as the new bin content
  // of the considered bin.
  
  for (Int_t binY = 1; binY <= h->GetYaxis()->GetNbins(); binY++) {
    Int_t firstBinBelow = TMath::Max(1, binY - 1);
    Int_t lastBinAbove = TMath::Min(h->GetYaxis()->GetNbins(), binY + 1);
  
    for (Int_t binX = 1; binX <= h->GetXaxis()->GetNbins(); binX++) {
      Int_t firstBinLeft = TMath::Max(1, binX - 1);
      Int_t lastBinRight = TMath::Min(h->GetXaxis()->GetNbins(), binX + 1);
    
      Double_t binContent = h->GetBinContent(binX, binY);
    
    
      
      // Reset values
      Int_t binIndex = 0;
      for (Int_t i = 0; i < nNeighbours; i++)
        values[i] = 0;
      
      for (Int_t binX2 = firstBinLeft; binX2 <= lastBinRight; binX2++) {
        for (Int_t binY2 = firstBinBelow; binY2 <= lastBinAbove; binY2++) {
          if (binX2 == binX && binY2 == binY)
            continue; // skip bin that is currently under consideration
          
          // Only take values from (hopefully) valid fits
          if (h->GetBinContent(binX2, binY2) > 0) {
            values[binIndex] = h->GetBinContent(binX2, binY2);
            binIndex++;
          }
        }
      }
      
      Double_t temp = getMedianOfNonZeros(values, nNeighbours);
      if (temp <= 0) {
        // Only print error message, if there is at least one positive value
        if (binIndex > 0)
          printf("Error: Could not eliminate outlier for bin at (%f, %f)!\n",
                 h->GetXaxis()->GetBinCenter(binX), h->GetYaxis()->GetBinCenter(binY));
        temp = -1;
      }
      else {
        if (TMath::Abs(binContent - temp) > threshold) {
          printf("Eliminated outlier at bin (%f, %f): %f -> %f!\n",
               h->GetXaxis()->GetBinCenter(binX), h->GetYaxis()->GetBinCenter(binY),
               binContent, temp);
          h->SetBinContent(binX, binY, temp);
          h->SetBinError(binX, binY, 1000);
        }
      }
    }
  }
  
  delete values;
}


//__________________________________________________________________________________________
void addPointToHyperplane(TH2D* h, TLinearFitter* linExtrapolation, Int_t binX, Int_t binY)
{
  if (h->GetBinContent(binX, binY) <= binContentThreshold)
    return; // Reject bins without content (within some numerical precision) or with strange content
    
  Double_t coord[2] = {0, 0};
  coord[0] = h->GetXaxis()->GetBinCenter(binX);
  coord[1] = h->GetYaxis()->GetBinCenter(binY);
  Double_t binError = h->GetBinError(binX, binY);
  if (binError <= 0) {
    printf("Warning: Trying to add bin in addPointToHyperplane with error (%e) not set (%f, %f), but bin content %f. Maybe fit problem -> Setting large error\n",
           binError, coord[0], coord[1], h->GetBinContent(binX, binY));
    binError = 1000;
  }
  linExtrapolation->AddPoint(coord, h->GetBinContent(binX, binY, binError));
}

//__________________________________________________________________________________________
TH2D* refineHistoViaLinearExtrapolation(TH2D* h, Double_t refineFactorX, Double_t refineFactorY, Int_t mapType, 
                                        TFile* fSave, Bool_t sigmaMap = kFALSE, Bool_t skipBinsAtBoundaries = kFALSE)
{
  if (!h)
    return 0x0;
  
  if (h->GetEntries() <= 0)
    return 0x0;
  
  Int_t nBinsX = h->GetXaxis()->GetNbins();
  Int_t nBinsY = h->GetYaxis()->GetNbins();

  Int_t firstYbin = 1;
  
  if (skipBinsAtBoundaries) {
    // Remove the two first and the last bin on y-axis
    nBinsY -= 3;
    firstYbin = 3; 
  }
  
  // Extrapolate map to larger 1/dEdx using the last 3 points (taking into account skipBinsAtBoundaries). Then: Refine
  
  Int_t nBinsYextrapolated = nBinsY + TMath::Ceil((1. / 20. - h->GetYaxis()->GetBinUpEdge(nBinsY)) / h->GetYaxis()->GetBinWidth(nBinsY));
  // Make sure that the "old" bins stay the same and only some additional bins are added (i.e. bin width should NOT change)
  Double_t yUpBoundExtrapolated =  h->GetYaxis()->GetBinUpEdge(nBinsY) +  h->GetYaxis()->GetBinWidth(nBinsY) * (nBinsYextrapolated - nBinsY);
  
  TH2D* hExtrapolated = new TH2D(Form("%s_%s_extrapolated", h->GetName(), suffixMapType[mapType].Data()), Form("%s (refined)",h->GetTitle()),
                                 nBinsX, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nBinsX),
                                 nBinsYextrapolated, h->GetYaxis()->GetBinLowEdge(1), yUpBoundExtrapolated);
  
  // Copy the old bins and take into account skipBinsAtBoundaries
  // Get rid of most non-positive bin contents
  eliminateNonPositivePointsInHisto(h);
  
  for (Int_t binX = 1; binX <= nBinsX; binX++)  {
    for (Int_t binY = 1; binY <= nBinsY; binY++)    {
      if (h->GetBinContent(binX, binY) > 0 && h->GetBinError(binX, binY) > 0)   {
        hExtrapolated->SetBinContent(binX, binY, h->GetBinContent(binX, binY));
        hExtrapolated->SetBinError(binX, binY, h->GetBinError(binX, binY));
      }
    }
    // Default value for extrapolated points is 1
    for (Int_t binY = nBinsY + 1; binY <= nBinsYextrapolated; binY++)    {
      hExtrapolated->SetBinContent(binX, binY, 1); 
      hExtrapolated->SetBinError(binX, binY, 1000); 
    }
  }
  
  
  // Do the extrapolation
  
  TLinearFitter* linExtrapolation = new TLinearFitter(2, "hyp2", "");
  
  for (Int_t binX = 1; binX <= nBinsX; binX++)  {
    linExtrapolation->ClearPoints();
    
    // Last 3 points for fitting for eta map,
    // but only last 2 points for sigma map, since the resolution
    // usually levels off around MIP and 3 points used for fitting
    // would cause an underestimation of the resolution.
    
    for (Int_t binY = (sigmaMap ? (nBinsY - 1) : (nBinsY - 2)); binY <= nBinsY; binY++)   {
      Double_t centerX = hExtrapolated->GetXaxis()->GetBinCenter(binX);
      Double_t centerY = hExtrapolated->GetYaxis()->GetBinCenter(binY);
      
      Int_t oldBinX = h->GetXaxis()->FindBin(centerX);
      if (oldBinX < 1)  
        oldBinX = 1;
      if (oldBinX > nBinsX)
        oldBinX = nBinsX;
      
      Int_t oldBinY = h->GetYaxis()->FindBin(centerY);
      if (oldBinY < firstYbin)  
        oldBinY = firstYbin;
      if (oldBinY > nBinsY)
        oldBinY = nBinsY;
      
      // Neighbours left column
        if (oldBinX >= 2) {
          if (oldBinY >= 2) {
            addPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY - 1);
          }
          
          addPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY);
          
          if (oldBinY < nBinsY) {
            addPointToHyperplane(h, linExtrapolation, oldBinX - 1, oldBinY + 1);
          }
        }
        
        // Neighbours (and point itself) same column
        if (oldBinY >= 2) {
          addPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY - 1);
        }
        
        addPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY);
        
        if (oldBinY < nBinsY) {
          addPointToHyperplane(h, linExtrapolation, oldBinX, oldBinY + 1);
        }
        
        // Neighbours right column
        if (oldBinX < nBinsX) {
          if (oldBinY >= 2) {
            addPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY - 1);
          }
          
          addPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY);
          
          if (oldBinY < nBinsY) {
            addPointToHyperplane(h, linExtrapolation, oldBinX + 1, oldBinY + 1);
          }
        }
    } 
    // Fit 2D-hyperplane
    if (linExtrapolation->GetNpoints() <= 0)
      continue;
    
    if (linExtrapolation->Eval() != 0)// EvalRobust -> Takes much, much, [...], much more time (~hours instead of seconds)
        continue;
        
    // ....extrapolation to the rest
    for (Int_t binY = nBinsY + 1; binY <= nBinsYextrapolated; binY++)    {
      Double_t centerX = hExtrapolated->GetXaxis()->GetBinCenter(binX);
      Double_t centerY = hExtrapolated->GetYaxis()->GetBinCenter(binY);
      
      // Fill the bin of the refined histogram with the extrapolated value
      Double_t extrapolatedValue = linExtrapolation->GetParameter(0) + linExtrapolation->GetParameter(1) * centerX
      + linExtrapolation->GetParameter(2) * centerY;
      
      // Do not allow too small or even negative values
      extrapolatedValue = TMath::Max(binContentThreshold * 2, extrapolatedValue);
      hExtrapolated->SetBinContent(binX, binY, extrapolatedValue);
      hExtrapolated->SetBinError(binX, binY, 50); // High, but constant, value => should almost not influence not-extrapolated data      
    }
  }
  
  eliminateNonPositivePointsInHisto(hExtrapolated);
  
  // Normalise map on demand
  
  // Use kNoNormalisation for final QA
  if (mapType == kSmallEtaNormalisation) {
    normaliseHisto(hExtrapolated, 0., 0.11, kTypeDeltaPrime);
  }
  else if (mapType == kLargeEtaNormalisation) {
    normaliseHisto(hExtrapolated, 0.81, 0.99, kTypeDeltaPrime);
  }
  
  // Interpolate to finer map
  Int_t nBinsXrefined = nBinsX * refineFactorX;
  Int_t nBinsYrefined = nBinsYextrapolated * refineFactorY; 
  
  TH2D* hRefined = new TH2D(Form("%s_%s_refined", h->GetName(), suffixMapType[mapType].Data()),  Form("%s (refined)",h->GetTitle()),
                            nBinsXrefined, hExtrapolated->GetXaxis()->GetBinLowEdge(1), hExtrapolated->GetXaxis()->GetBinUpEdge(nBinsX),
                            nBinsYrefined, hExtrapolated->GetYaxis()->GetBinLowEdge(firstYbin), 
                            hExtrapolated->GetYaxis()->GetBinUpEdge(nBinsYextrapolated));
  
  for (Int_t binX = 1; binX <= nBinsXrefined; binX++)  {
    Double_t centerX = hRefined->GetXaxis()->GetBinCenter(binX);
    
    for (Int_t binY = 1; binY <= nBinsYrefined; binY++)  {

      hRefined->SetBinContent(binX, binY, 1); // Default value is 1
      
      Double_t centerY = hRefined->GetYaxis()->GetBinCenter(binY);
      
      /*OLD
      linExtrapolation->ClearPoints();
      
      // For interpolation: Just take the corresponding bin from the old histo.
      // For extrapolation: take the last available bin from the old histo.
      // If the boundaries are to be skipped, also skip the corresponding bins
      Int_t oldBinX = hExtrapolated->GetXaxis()->FindBin(centerX);
      if (oldBinX < 1)  
        oldBinX = 1;
      if (oldBinX > nBinsX)
        oldBinX = nBinsX;
      
      Int_t oldBinY = hExtrapolated->GetYaxis()->FindBin(centerY);
      if (oldBinY < firstYbin)  
        oldBinY = firstYbin;
      if (oldBinY > nBinsYextrapolated)
        oldBinY = nBinsYextrapolated;
      
      // Neighbours left column
      if (oldBinX >= 2) {
        if (oldBinY >= 2) {
          addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX - 1, oldBinY - 1);
        }
        
        addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX - 1, oldBinY);
        
        if (oldBinY < nBinsYextrapolated) {
          addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX - 1, oldBinY + 1);
        }
      }
      
      // Neighbours (and point itself) same column
      if (oldBinY >= 2) {
        addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX, oldBinY - 1);
      }
        
      addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX, oldBinY);
        
      if (oldBinY < nBinsYextrapolated) {
        addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX, oldBinY + 1);
      }
      
      // Neighbours right column
      if (oldBinX < nBinsX) {
        if (oldBinY >= 2) {
          addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX + 1, oldBinY - 1);
        }
        
        addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX + 1, oldBinY);
        
        if (oldBinY < nBinsYextrapolated) {
          addPointToHyperplane(hExtrapolated, linExtrapolation, oldBinX + 1, oldBinY + 1);
        }
      }
      
      
      // Fit 2D-hyperplane
      if (linExtrapolation->GetNpoints() <= 0)
        continue;
        
      if (linExtrapolation->Eval() != 0)// EvalRobust -> Takes much, much, [...], much more time (~hours instead of seconds)
        continue;
      
      // Fill the bin of the refined histogram with the extrapolated value
      Double_t interpolatedValue = linExtrapolation->GetParameter(0) + linExtrapolation->GetParameter(1) * centerX
                                 + linExtrapolation->GetParameter(2) * centerY;
      */
      Double_t interpolatedValue = hExtrapolated->Interpolate(centerX, centerY);
      hRefined->SetBinContent(binX, binY, interpolatedValue);      
    }
  } 
  
  
  // Problem: Interpolation does not work before/beyond center of first/last bin (as the name suggests).
  // Therefore, for each row in dEdx: Take last bin from old map and interpolate values from center and edge.
  // Assume line through these points and extropolate to last bin of refined map
  const Double_t firstOldXbinUpEdge = hExtrapolated->GetXaxis()->GetBinUpEdge(1);
  const Double_t firstOldXbinCenter = hExtrapolated->GetXaxis()->GetBinCenter(1);
  
  const Double_t oldXbinHalfWidth = firstOldXbinUpEdge - firstOldXbinCenter;
  
  const Double_t lastOldXbinLowEdge = hExtrapolated->GetXaxis()->GetBinLowEdge(hExtrapolated->GetNbinsX());
  const Double_t lastOldXbinCenter = hExtrapolated->GetXaxis()->GetBinCenter(hExtrapolated->GetNbinsX());
  
  for (Int_t binY = 1; binY <= nBinsYrefined; binY++)  {
    Double_t centerY = hRefined->GetYaxis()->GetBinCenter(binY);
    
    const Double_t interpolatedCenterFirstXbin = hExtrapolated->Interpolate(firstOldXbinCenter, centerY);
    const Double_t interpolatedUpEdgeFirstXbin = hExtrapolated->Interpolate(firstOldXbinUpEdge, centerY);
    
    const Double_t extrapolationSlopeFirstXbin = (interpolatedUpEdgeFirstXbin - interpolatedCenterFirstXbin) / oldXbinHalfWidth;
    const Double_t extrapolationOffsetFirstXbin = interpolatedCenterFirstXbin;
    
    
    const Double_t interpolatedCenterLastXbin = hExtrapolated->Interpolate(lastOldXbinCenter, centerY);
    const Double_t interpolatedLowEdgeLastXbin = hExtrapolated->Interpolate(lastOldXbinLowEdge, centerY);
    
    const Double_t extrapolationSlopeLastXbin = (interpolatedCenterLastXbin - interpolatedLowEdgeLastXbin) / oldXbinHalfWidth;
    const Double_t extrapolationOffsetLastXbin = interpolatedCenterLastXbin;

    for (Int_t binX = 1; binX <= nBinsXrefined; binX++)  {
      Double_t centerX = hRefined->GetXaxis()->GetBinCenter(binX);
     
      if (centerX < firstOldXbinCenter) {
        Double_t extrapolatedValue = extrapolationOffsetFirstXbin + (centerX - firstOldXbinCenter) * extrapolationSlopeFirstXbin;
        hRefined->SetBinContent(binX, binY, extrapolatedValue);      
      }
      else if (centerX <= lastOldXbinCenter) {
        continue;
      }
      else {
        Double_t extrapolatedValue = extrapolationOffsetLastXbin + (centerX - lastOldXbinCenter) * extrapolationSlopeLastXbin;
        hRefined->SetBinContent(binX, binY, extrapolatedValue);     
      }
    }
  } 
  
  delete linExtrapolation;
  
  
  if (fSave)    {
    fSave->cd();
    hExtrapolated->Write();
  }
  
  delete hExtrapolated;
  
  return hRefined;
}


//__________________________________________________________________________________________
TFitResult*  doubleGaussFit(TH1D* h, Double_t currentMeanMomentum, TSpline3* splPr, TSpline3* splPi, TString fitOption = "QNS") 
{
  Double_t lowThreshold = 0.6;
  Double_t highThreshold = 1.6;
  FindFitRange(h, lowThreshold, highThreshold);
  TFitResultPtr result = h->Fit("gaus", "QNS", "", lowThreshold, highThreshold);
  
  Double_t contaminationPeakMean = splPi->Eval(currentMeanMomentum / massPion) / splPr->Eval(currentMeanMomentum / massProton);

  if (contaminationPeakMean < h->GetXaxis()->GetBinLowEdge(1) ||
      contaminationPeakMean > h->GetXaxis()->GetBinUpEdge(h->GetNbinsX())) {
    return ((Int_t)result == 0) ? (TFitResult*)result.Get()->Clone() : 0x0;
  }
 
  // Estimate parameters for doubleGauss fit
  Double_t estimatedMean = 0;
  Double_t estimatedSigma = 0;
  Double_t estimatedYield = 0;
  
  Double_t chi2oneGauss = 0;
  Int_t NDFoneGauss = 0;
  Double_t reducedChi2oneGauss = 0;
  
  if ((Int_t) result == 0) {
    estimatedMean = result->GetParams()[1];
    estimatedSigma = result->GetParams()[2];
    estimatedYield = result->GetParams()[0];
    
    chi2oneGauss = result->Chi2();
    NDFoneGauss = result->Ndf();
    reducedChi2oneGauss = (NDFoneGauss > 0) ? chi2oneGauss / NDFoneGauss : -1;
  }
  else {
    estimatedMean = h->GetMean();
    estimatedSigma = h->GetRMS();
    estimatedYield = (estimatedSigma > 0) ? (h->Integral("width") / estimatedSigma) : h->GetEntries();
  }
  
  TF1* doubleGaus = new TF1("doubleGaus", "[0]*TMath::Gaus(x,[1],[2],0)+[3]*TMath::Gaus(x,[4],[2],0)", 0.6, 1.6);

  estimatedMean = TMath::Max(0.6, estimatedMean);
  estimatedYield = TMath::Max(1., estimatedYield);
  
  Double_t newPars[5] = { estimatedYield, estimatedMean, estimatedSigma, estimatedYield / 10., contaminationPeakMean };
  doubleGaus->SetParameters(newPars);
  doubleGaus->SetParLimits(0, 0.01 * estimatedYield, 100. * estimatedYield);//TODO 0., 100. * estimatedYield);
  doubleGaus->SetParLimits(1, estimatedMean - 0.05, estimatedMean + 0.05);//TODO 0.6, 1.6);
  doubleGaus->SetParLimits(2, 0, 100. * estimatedSigma);
  doubleGaus->SetParLimits(3, 0, 0.8 * estimatedYield);
  doubleGaus->SetParLimits(4, contaminationPeakMean - 0.15, contaminationPeakMean + 0.15);//TODO doubleGaus->SetParLimits(4, 0.6, 1.6);
  TFitResultPtr result2 = h->Fit(doubleGaus, fitOption.Data());

  printf("\n%f -> %f\n", currentMeanMomentum, contaminationPeakMean);//TODO NOW NOW NOW
  printf("%f, %f, %f, %f, %f\n%f, %f, %f\n", result2->GetParams()[0], result2->GetParams()[1], result2->GetParams()[2], result2->GetParams()[3], result2->GetParams()[4], result->GetParams()[0], result->GetParams()[1], result->GetParams()[2]); printedSomething = kTRUE;//TODO NOW NOW NOW
  if ((Int_t)result2 == 0) {
    Double_t chi2doubleGauss = 0;
    Int_t NDFdoubleGauss = 0;
    Double_t reducedChi2doubleGauss = 0;
    
    chi2doubleGauss = result2->Chi2();
    NDFdoubleGauss = result2->Ndf();
    reducedChi2doubleGauss = (NDFdoubleGauss > 0) ? chi2doubleGauss / NDFdoubleGauss : -1;
    
    printf("%f / %d = %f\n %f / %d = %f\n\n", chi2oneGauss, NDFoneGauss, reducedChi2oneGauss, chi2doubleGauss, NDFdoubleGauss, reducedChi2doubleGauss); printedSomething = kTRUE;//TODO NOW NOW NOW
    // Only take double gauss result, if (valid) reduced chi2 is better than that of the one gauss fit
    if (reducedChi2oneGauss < 0 || (reducedChi2doubleGauss > 0 && reducedChi2doubleGauss <= reducedChi2oneGauss))
      return (TFitResult*)result2.Get()->Clone();
  }
  
  // If fit failed, return results of standard fit instead
  return ((Int_t)result == 0) ? (TFitResult*)result.Get()->Clone() : 0x0;
}


//__________________________________________________________________________________________
TFitResult* extractClusterDependence(TH3D** hPreMapClusterResolved, Int_t numClusterBins, Int_t clusterLowBound, Int_t clusterUpBound,
                                     Int_t binX, Int_t binY, Int_t binYhigh, Double_t c0, TFile* fSave, TSpline3* splPr, TSpline3* splPi)
{
  TH1D* hClusMean = new TH1D(Form("hClusMean_tanTheta_%f_pTPC_%f", hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                  hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                             Form("#Delta' vs. ncl (tan(#Theta) %f, p_{TPC} %f);ncl;#Delta'", hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                  hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                             numClusterBins, clusterLowBound, clusterUpBound);
  
  TH1D* hClusSigmaRelGauss = new TH1D(Form("hClusSigmaRelGauss_tanTheta_%f_pTPC_%f", hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                           hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                                      Form("#sigma_{rel, Gauss} vs. ncl (tan(#Theta) %f, p_{TPC} %f);ncl;#sigma_{rel, Gauss}",
                                           hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                           hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                                      numClusterBins, clusterLowBound, clusterUpBound);

  TH1D* hClusSigmaRel = new TH1D(Form("hClusSigmaRel_tanTheta_%f_pTPC_%f", hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                           hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                                      Form("#sigma_{rel} vs. ncl (tan(#Theta) %f, p_{TPC} %f);ncl;#sigma_{rel}",
                                           hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                           hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY)),
                                      numClusterBins, clusterLowBound, clusterUpBound);
  
  /*//TODO ADDED 03.07.2013
  TH1D* hTest = hPreMapClusterResolved[0]->ProjectionZ("hTest", binX, binX, binY, binYhigh);
  Long64_t nEntries = hTest->GetEntries();
  Double_t clusterMean = (clusterLowBound + (clusterUpBound - clusterLowBound) / numClusterBins * 0) * nEntries;
  
  for (Int_t clusterBin = 1; clusterBin < numClusterBins; clusterBin++) {
    TH1D* hTempProjectionZ = hPreMapClusterResolved[clusterBin]->ProjectionZ("hTempProjectionZ", binX, binX, binY, binYhigh);

    if (hTempProjectionZ->GetEntries() < 10) {
      delete hTempProjectionZ;
      continue;
    }
    
    hTest->Add(hTempProjectionZ);
    nEntries += hTempProjectionZ->GetEntries();
    clusterMean += (clusterLowBound + (clusterUpBound - clusterLowBound) / numClusterBins * clusterBin) * hTempProjectionZ->GetEntries();
    delete hTempProjectionZ;
  }
  if (nEntries > 0)
    clusterMean /= nEntries;
  
  Double_t lowThreshold = 0.6;
  Double_t highThreshold = 1.6;
  FindFitRange(hTest, lowThreshold, highThreshold);
  TFitResultPtr fitResPtr3 = hTest->Fit("gaus", "QNS", "", lowThreshold, highThreshold);
  TFitResult* fitRes3 = ((Int_t)fitResPtr3 == 0) ? (TFitResult*)fitResPtr3.Get()->Clone() : 0x0;
  
  if (!fitRes3 || fitRes3->GetParams()[1] <= 0) {
    printf("ERROR: Fit test failed\n");
    return 0x0;
  }
  
  Double_t mean = fitRes3->GetParams()[1];
  Double_t sigma = fitRes3->GetParams()[2];
  
  Double_t relSigma = sigma / mean;
  
  TGraphErrors * grClusSigmaRel2 = new TGraphErrors(1);
  grClusSigmaRel2->SetPoint(0, clusterMean, relSigma);
  
  TF1* fStats2 = new TF1("fStats2", "TMath::Sqrt([0] * [0] + [1] * [1] / x)", clusterLowBound, clusterUpBound);
  fStats2->SetLineColor(kBlue);
  fStats2->FixParameter(0, c0);
  //fStats2->SetParameter(1, 0.1);
  //fStats2->SetParLimits(1, 0, 2);
  fStats2->SetNpx(300);
  
  TFitResultPtr fitRes14 = grClusSigmaRel2->Fit(fStats2, "QS", "", clusterLowBound, clusterUpBound);
  
  delete fStats2;
  delete grClusSigmaRel2;
  
  delete hTest;

  if (((Int_t)fitRes14) != 0) {
    printf("\nError: Cluster fit: Dependence fit failed (%d) for (%f, %f - %f): skipped\n", (Int_t)fitRes14,
           hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY),
           hPreMapClusterResolved[0]->GetYaxis()->GetBinUpEdge(binYhigh));
    printedSomething = kTRUE;
    return 0x0;
  }
  
  return fitRes14.Get();
  //TODO END*/
  
  Double_t minLowThreshold = 0.6;
  Double_t maxHighThreshold = 1.6;
  fGauss.SetRange(minLowThreshold, maxHighThreshold);
  
  for (Int_t clusterBin = 0; clusterBin < numClusterBins; clusterBin++) {
    TH1D* hTempProjectionZ = hPreMapClusterResolved[clusterBin]->ProjectionZ("hTempProjectionZ", binX, binX, binY, binYhigh);

    if (hTempProjectionZ->GetEntries() < 10) {
      /*
      printf("\nWarning: Cluster fit: No entries for (%f, %f - %f, ncl %f): skipped\n",
             hPreMapClusterResolved[clusterBin]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinLowEdge(binY),
             hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinUpEdge(binYhigh), hClusMean->GetXaxis()->GetBinCenter(clusterBin + 1));
      printedSomething = kTRUE;

      */
      delete hTempProjectionZ;
      continue;
    }

    TFitResult* fitRes = 0x0;
    
    Double_t mean = -999.;
    Double_t errorMean = 999.;

    Double_t sigma = 999.;
    Double_t errorSigma = 999.;
      
    if (splPr && splPi) { // do DoubleGaussFit
      fitRes = doubleGaussFit(hTempProjectionZ, 
                              0.5 * (hPreMapClusterResolved[0]->GetYaxis()->GetBinCenter(binY) + 
                                     hPreMapClusterResolved[0]->GetYaxis()->GetBinCenter(binYhigh)),
                              splPr, splPi, "QNS");
      
      if (fitRes) {
        mean = fitRes->GetParams()[1];
        errorMean = fitRes->GetErrors()[1];

        sigma = fitRes->GetParams()[2];
        errorSigma = fitRes->GetErrors()[2];
      }
    }
    else {
      Double_t lowThreshold = minLowThreshold;
      Double_t highThreshold = maxHighThreshold;
      
      // First, fit with the desired restriction via FindFitRange to fix the mean,
      // but to fit the width, fit the full range (otherwise one typically gets biased
      // to larger sigmas). However, fitting the full range would also result in a different
      // mean. Thus, to be consistent, these two steps are needed.
      // NOTE: If the restricted fit fails, just take the mean of the fit over the full range (if it doesn't fail),
      // since this is still better than no information
      
      FindFitRange(hTempProjectionZ, lowThreshold, highThreshold);
      TFitResultPtr fitResPtrFirstStep = hTempProjectionZ->Fit("gaus", "QNS", "", lowThreshold, highThreshold);
      
      if ((Int_t)fitResPtrFirstStep == 0) {
        fGauss.SetParameter(0, fitResPtrFirstStep->GetParams()[0]);
        fGauss.SetParError(0, fitResPtrFirstStep->GetErrors()[0]);
        fGauss.SetParameter(2, fitResPtrFirstStep->GetParams()[2]);
        fGauss.SetParError(2, fitResPtrFirstStep->GetErrors()[2]);
        
        fGauss.FixParameter(1, fitResPtrFirstStep->GetParams()[1]);
        TFitResultPtr fitResPtrSecondStep = hTempProjectionZ->Fit(&fGauss, "QNS", "", minLowThreshold, maxHighThreshold);
        
        fitRes = ((Int_t)fitResPtrSecondStep == 0) ? (TFitResult*)fitResPtrSecondStep.Get()->Clone() : 0x0;
        
        if (fitRes) {
          // Mean and error from first step (error is zero in second step because mean is fixed there)
          mean = fitResPtrFirstStep->GetParams()[1];
          errorMean = fitResPtrFirstStep->GetErrors()[1];

          sigma = fitRes->GetParams()[2];
          errorSigma = fitRes->GetErrors()[2];
        }
      }
      else {
        TFitResultPtr fitResPtrSecondStep = hTempProjectionZ->Fit("gaus", "QNS", "", minLowThreshold, maxHighThreshold);
        
        fitRes = ((Int_t)fitResPtrSecondStep == 0) ? (TFitResult*)fitResPtrSecondStep.Get()->Clone() : 0x0;
        
        if (fitRes) {
          mean = fitRes->GetParams()[1];
          errorMean = fitRes->GetErrors()[1];

          sigma = fitRes->GetParams()[2];
          errorSigma = fitRes->GetErrors()[2];
        }
      }
    }
            
    if (!fitRes) {
      printf("\nWarning: Cluster fit: Gauss fit failed for (%f, %f - %f, ncl %f): skipped\n",
             hPreMapClusterResolved[clusterBin]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinLowEdge(binY),
             hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinUpEdge(binYhigh), hClusMean->GetXaxis()->GetBinCenter(clusterBin + 1));
      printedSomething = kTRUE;

      delete hTempProjectionZ;
      continue;
    }
    else {
      hClusMean->SetBinContent(clusterBin + 1, mean);
      hClusMean->SetBinError(clusterBin + 1, errorMean);
      
      hClusSigmaRelGauss->SetBinContent(clusterBin + 1, sigma);
      hClusSigmaRelGauss->SetBinError(clusterBin + 1, errorSigma);

      if (mean > 1e-4)   {
        hClusSigmaRel->SetBinContent(clusterBin + 1, sigma / mean);
        hClusSigmaRel->SetBinError(clusterBin + 1,
                                   TMath::Sqrt(TMath::Power(errorSigma / mean, 2) + TMath::Power(errorMean * sigma / (mean * mean), 2)));
      }
      else {
        printf("\nWarning: Cluster fit: Gauss fit with strange mean value (%e) for (%f, %f - %f, ncl %f): skipped\n", mean,
             hPreMapClusterResolved[clusterBin]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinLowEdge(binY),
             hPreMapClusterResolved[clusterBin]->GetYaxis()->GetBinUpEdge(binYhigh), hClusMean->GetXaxis()->GetBinCenter(clusterBin + 1));
        printedSomething = kTRUE;
      }

      delete hTempProjectionZ;
    }
  }

  TF1* fStats = new TF1("fStats", "TMath::Sqrt([0] * [0] + [1] * [1] / x)", clusterLowBound, clusterUpBound);
  fStats->SetLineColor(kBlue);
  fStats->FixParameter(0, c0);
  //fStats->SetParameter(1, 0.1);
  //fStats->SetParLimits(1, 0, 2);
  fStats->SetNpx(300);
  
  TGraphErrors * grClusSigmaRel = new TGraphErrors(hClusSigmaRel);
  for(Int_t ip = 0; ip < grClusSigmaRel->GetN(); ip ++) {
    Bool_t removePoint = grClusSigmaRel->GetY()[ip] <= 0 || grClusSigmaRel->GetEY()[ip] / TMath::Max(1e-8, grClusSigmaRel->GetY()[ip]) > 0.20; //TODO NOW was > 0.10
                         //|| grClusSigmaRel->GetX()[ip] > 140; // TODO (was 160) For some reason it is better to skip the last clusters 
    if (removePoint) {
      grClusSigmaRel->RemovePoint(ip);
      ip--;
      continue;
    }
  }
  grClusSigmaRel->SetMarkerStyle(20);
  grClusSigmaRel->SetMarkerColor(kMagenta);
  
  if (grClusSigmaRel->GetN() <= 0) {
    printf("\nError: Cluster fit: No good data points for dependence fit for (%f, %f - %f): skipped\n", 
           hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY),
           hPreMapClusterResolved[0]->GetYaxis()->GetBinUpEdge(binYhigh));
    printedSomething = kTRUE;
    return 0x0;
  }
  
  TFitResultPtr fitRes = grClusSigmaRel->Fit(fStats, "QS", "", clusterLowBound, clusterUpBound);
  //printf("\nchi^2/NDF = %f / %d = %f\n\n", fitRes->Chi2(), fitRes->Ndf(), (fitRes->Ndf() > 0) ? fitRes->Chi2() / fitRes->Ndf() : -1);

  TCanvas* canvCluster = new TCanvas(Form("canvCluster_tanTheta_%f_pTPC_%f_%f", 
                                          hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                          hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY),
                                          hPreMapClusterResolved[0]->GetYaxis()->GetBinUpEdge(binYhigh)),
                                     Form("canvCluster_tanTheta_%f_pTPC_%f_%f", 
                                          hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX),
                                          hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY),
                                          hPreMapClusterResolved[0]->GetYaxis()->GetBinUpEdge(binYhigh)),
                                     100,10,1380,800);
  canvCluster->Divide(3,1);

  canvCluster->cd(1);
  hClusMean->DrawCopy("e");
  
  canvCluster->cd(2);
  hClusSigmaRelGauss->DrawCopy("e");

  canvCluster->cd(3);
  hClusSigmaRel->DrawCopy("e");
  grClusSigmaRel->Draw("same p");
  
  if (fSave)  {
    fSave->cd("clusterQA");
    canvCluster->Write();
    //hClusMean->Write();
    //hClusSigmaRelGauss->Write();
    //hClusSigmaRel->Write();
  }

  delete fStats;
  delete grClusSigmaRel;
  
  delete hClusMean;
  delete hClusSigmaRel;
  delete hClusSigmaRelGauss;

  delete canvCluster;

  if (((Int_t)fitRes) != 0) {
    printf("\nError: Cluster fit: Dependence fit failed (%d) for (%f, %f - %f): skipped\n", (Int_t)fitRes,
           hPreMapClusterResolved[0]->GetXaxis()->GetBinCenter(binX), hPreMapClusterResolved[0]->GetYaxis()->GetBinLowEdge(binY),
           hPreMapClusterResolved[0]->GetYaxis()->GetBinUpEdge(binYhigh));
    printedSomething = kTRUE;
    return 0x0;
  }
  
  return fitRes.Get();
}

//__________________________________________________________________________________________
// NOTE: Use smoothing only for creation of sigma map. For normal eta map this pulls down the edges too much <-> trade-off between
// (normally very small) jumps in the map and good description at the edges.
// For the sigma map, there are not that sharp edges, but the fluctuations are by far stronger. Here it makes sense to use the smoothing.
TH2* checkShapeEtaTree(TTree* tree = 0x0, Double_t resolutionFactorX = 1/*2 or 2.1*/, Double_t resolutionFactorY = 1/*2 or 2.1*/, 
                       Double_t pThresholdTOFcut = -1/*0.6*/, Int_t nMergeBinsAroundThreshold = 3, Double_t pThresholdV0cut = 1.2,
                       Double_t pThresholdV0plusTOFcut = 2.0, Double_t outlierThreshold = 0.04,
                       Bool_t doSmoothing = kFALSE, Double_t c0 = 0.0232,  Int_t nclCut = 60,
                       TString path = ".", TString suffixForFileName = "", Bool_t useDoubleGaussFit = kFALSE, 
                       TString pathNameSplinesFile = "TPCPIDResponse.root", TString prSplinesName = "TSPLINE3_DATA_PROTON_LHC10H_PASS2_PBPB_MEAN", 
                       Int_t returnMapType = kNoNormalisation, TString treeName = "fTree") 
{ 
  TString fileName = Form("bhess_PIDetaTree%s.root", suffixForFileName.Data());
  
  if (resolutionFactorX <= 0 || resolutionFactorY <= 0)  {
	  printf("Factor for lower resolution <= 0 is not allowed!\n");
	  return 0x0;
  }
  
  TFile* f = 0x0;
  
  Bool_t saveResults = kFALSE;

  if (!tree)  {
    f = TFile::Open(Form("%s/%s", path.Data(), fileName.Data()));
    if (!f)  {
      std::cout << "Failed to open file \"" << Form("%s/%s", path.Data(), fileName.Data()) << "\"!" << std::endl;
      
      
      return 0x0;
    }
    
    // Extract the data Tree
    tree = dynamic_cast<TTree*>(f->Get(treeName.Data()));
    if (!tree) {
      std::cout << "Failed to load data tree!" << std::endl;
      return 0x0;
    }
    
    saveResults = kTRUE;
  }
  
  
  
  //TODO NOW start
  /*
  //TString pathNameMomentumMap = "pionTree_11b2/outputCheckShapeEtaTreePions_2012_10_16__11_38.root";
  //TString pathNameMomentumMap = "pionTree_10d1/outputCheckShapeEtaTreePions_2012_10_17__07_38.root";
  TString pathNameMomentumMap = "finalCuts/uncorrected/0.6GeVcut/pp/7TeV/LHC10d.pass2/correctedV0splines_newV0tree/outputCheckShapeEtaTreePions_2012_10_23__16_02_correctedWith_10_38.root";
  TFile* fPiMap = TFile::Open(pathNameMomentumMap.Data());
  if (!fPiMap)  {
    std::cout << "Failed to open pion momentum map file \"" << pathNameMomentumMap.Data() << "\"!" << std::endl;
    return 0x0;
  }

  TH2D* hPiMap = 0x0;
  
  if (fPiMap) {
    hPiMap = dynamic_cast<TH2D*>(fPiMap->Get("hRefinedSmallEtaNormalisation")); //TODO NOW NOW NOWForm(hRefined%s", suffixMapType[returnMapType].Data())));
    if (!hPiMap) {
      std::cout << "Failed to load pion momentum map!" << std::endl;
      return 0x0;
    }
  }*/
  //TODO NOW end
  
  
  
  
  TSpline3* splPr = 0x0;
  TSpline3* splPi = 0x0;
  
  if (useDoubleGaussFit) {
    std::cout << "Using double gauss fit!" << std::endl << std::endl;
  
    std::cout << "Loading splines from file \"" << pathNameSplinesFile.Data() << "\"" << std::endl;
    
    TFile* fSpl = TFile::Open(pathNameSplinesFile.Data());
    if (!fSpl) {
      std::cout << "Failed to open spline file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
      return 0x0;
    }
    
    TObjArray* TPCPIDResponse = (TObjArray*)fSpl->Get("TPCPIDResponse");
    if (!TPCPIDResponse) {
      std::cout << "Failed to load object array from spline file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
      return 0x0;
    }
    
    splPr = (TSpline3*)TPCPIDResponse->FindObject(prSplinesName.Data());
    TString piSplinesName = prSplinesName.ReplaceAll("PROTON", "PION");
    splPi = (TSpline3*)TPCPIDResponse->FindObject(piSplinesName.Data());
    
    if (!splPr || !splPi) {
      std::cout << "Failed to load splines from file \"" << pathNameSplinesFile.Data() << "\"!" << std::endl;
      return 0x0;
    }
  }
  else
    std::cout << "NOT using double gauss fit!" << std::endl << std::endl;

  
  // Output file
  
  TFile* fSave = 0x0;
  TString saveFileName = "";
  if (saveResults)  {
    TDatime daTime;
    
    TString saveSuffix = "";
    if (suffixForFileName.Length() > 0)
      saveSuffix = Form("_%s", suffixForFileName.Data());
    
    saveFileName = Form("outputCheckShapeEtaTree_%04d_%02d_%02d__%02d_%02d%s.root", daTime.GetYear(), daTime.GetMonth(),
                        daTime.GetDay(), daTime.GetHour(), daTime.GetMinute(), saveSuffix.Data());
    
    fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
    if (!fSave) {
      std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
      return 0x0;
    }
    
    fSave->mkdir("clusterQA");
  }
  
  printf("\nProcessing file %s\n", f->GetName());
  
  // Only activate the branches of interest to save processing time
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  //tree->SetBranchStatus("pT", 1);
  //tree->SetBranchStatus("phiPrime", 1);
  tree->SetBranchStatus("dEdx", 1);
  tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  //if (isPbPb)
  //  tree->SetBranchStatus("fMultiplicity", 1);
  //tree->SetBranchStatus("sinAlpha", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
  tree->SetBranchStatus("pidType", 1);
  
  TString multiplicitySelectionString = "";
// if (isPbPb)// Usually 20 bins from 0 to 4000 => bin width ~ 200
//    multiplicitySelectionString = "&& fMultiplicity > 3000 && fMultiplicity <= 3200";
  
  // Bins along tanTheta
  const Int_t numTanThetaBins = 30; //WARNING: Also change in AliPIDResponse, if changed here!
  Int_t binsX = numTanThetaBins; 
  
  if (resolutionFactorX != 1) 
    binsX /= resolutionFactorX;
  
  // Do not mix tan(theta) > 0 with tan(theta) < 0!
  // Since the eta range is symmetric around zero (i.e. |lowerTanTheta| = |upperTanTheta|),
  // this means that the number of bins for tan(theta) < 0 must equal that of tan(theta) > 0.
  // Thus, in order to have a bin edge at tan(theta) = 0, binsX must be even.
  if (binsX % 2 != 0) {
    binsX++;
    // Set resolutionFactorX correspondingly to the really used value
    resolutionFactorX = ((Double_t)numTanThetaBins) / ((Double_t)binsX);
    printf("\nNOTE: In order not to mix up tan(theta) < 0 with > 0, an even number of tanTheta-bins is required => ResolutionFactorX set to %f!\n\n",
           resolutionFactorX);
  }
  
  // ---------------------------------------------------------------------------------
  // Use momentum bins and take expected dEdx to translate from momentum to dEdx in map

  const Int_t clusterUpBound = 160; 
  //TODO was 70 and 10; changed on 24.06.13
  const Int_t clusterLowBound = 60; // Other possibility: 80. Fits look better, but obviously the results for the PID-task are worse
  const Int_t clusterBinWidth = 10; // Other possibility: 20. Fits look better, but obviously the results for the PID-task are worse
  // NOTE: Too large cluster bins are bad because you fit a superposition of Gauss curves with different widths. This potentially
  // biases your fit to larger values. The effect for the resulting dependence 1/sqrt(ncl) is typically few permille for deltaNcl = 20
  // compared to deltaNcl = 10
  const Int_t numClusterBins = (clusterUpBound - clusterLowBound) / clusterBinWidth;
  
  const Double_t pBoundLow = 0.15;
  const Double_t pBoundUp = 5.0;
  const Double_t lowerTanTheta = -1.0;
  const Double_t upperTanTheta = 1.0;
  
  
  
  
  const Int_t nPbinsForMap = 110;
  Double_t binsPforMap[nPbinsForMap + 1];

  const Double_t factorBinning = TMath::Power(pBoundUp/pBoundLow, 1./nPbinsForMap);
  
  // Log binning
  binsPforMap[0] = pBoundLow;
  for (Int_t i = 0 + 1; i <= nPbinsForMap; i++) {
    binsPforMap[i] = factorBinning * binsPforMap[i - 1];
  }
  binsPforMap[nPbinsForMap] = pBoundUp;

  /*TODO changed on 23.06.13
  // 68 bins from 0.15 to 1.0 (width 0.0125)
  // 20 bins from 1.0 to 1.5 (width 0.025)
  // 10 bins from 1.5 to 2.0 (width 0.05)
  //  5 bins from 2.0 to 2.5 (width 0.1)
  //  4 bins from 2.5 to 3.5 (width 0.25)
  //  3 bins from 3.5 to 5.0 (width 0.5)
  const Int_t nPbinsForMap = 68 + 20 + 10 + 5 + 4 + 3;
  Double_t binsPforMap[nPbinsForMap + 1];
  for (Int_t i = 0; i < 68; i++)  {
    binsPforMap[i] = pBoundLow + i * 0.0125;
  }
  for (Int_t i = 68, j = 0; i < 68 + 20; i++, j++)  {
    binsPforMap[i] = 1.0 + j * 0.025;
  }
  for (Int_t i = 68 + 20, j = 0; i < 68 + 20 + 10; i++, j++)  {
    binsPforMap[i] = 1.5 + j * 0.05;
  }
  for (Int_t i = 68 + 20 + 10, j = 0; i < 68 + 20 + 10 + 5; i++, j++)  {
    binsPforMap[i] = 2.0 + j * 0.1;
  }
  for (Int_t i = 68 + 20 + 10 + 5, j = 0; i < 68 + 20 + 10 + 5 + 4; i++, j++)  {
    binsPforMap[i] = 2.5 + j * 0.25;
  }
  for (Int_t i = 68 + 20 + 10 + 5 + 4, j = 0; i < 68 + 20 + 10 + 5 + 4 + 3; i++, j++)  {
    binsPforMap[i] = 3.5 + j * 0.5;
  }
  binsPforMap[nPbinsForMap] = pBoundUp;
  */

  progressbar(0.);
  
  const Int_t nDim = 3;
  
  Int_t thnBins[nDim] = {  binsX, nPbinsForMap, 200  };
  Double_t xmin[nDim] = {  lowerTanTheta, pBoundLow, 0.6  };
  Double_t xmax[nDim] = {  upperTanTheta, pBoundUp, 1.6  };
  THnSparseD* hDummy = new THnSparseD("hDummy", "", nDim, thnBins, xmin, xmax);
  hDummy->SetBinEdges(1, binsPforMap);
  TH3D* hPreMap = hDummy->Projection(0, 1, 2);
  hPreMap->SetName("hPreMap");
  hPreMap->SetTitle("#Delta' map for p_{TPC_Inner} vs. tan(#theta)");   
  hPreMap->GetXaxis()->SetTitle("tan(#theta)");
  hPreMap->GetYaxis()->SetTitle("p_{TPC_inner} (GeV/c)");
  hPreMap->GetZaxis()->SetTitle("#Delta'");
  
  delete hDummy;
  
  // Use the expected dEdx to translate from pTPC to 1/dEdx -> Create a 2D histo with the correct momentum binning for this translation  
  const Int_t nDim2 = 2;
  
  Int_t thnBins2[nDim2] = {  nPbinsForMap, 200000  };
  Double_t xmin2[nDim2] = {  pBoundLow, 0.  };
  Double_t xmax2[nDim2] = {  pBoundUp, 2000.  };
  
  hDummy = new THnSparseD("hDummy", "", nDim2, thnBins2, xmin2, xmax2);
  hDummy->SetBinEdges(0, binsPforMap);
  TH2D* hMomTranslation = hDummy->Projection(1, 0);
  hMomTranslation->SetName("hMomTranslation");
  hMomTranslation->SetTitle("Momentum to dEdxExpected translation map)");   
  hMomTranslation->GetXaxis()->SetTitle("p_{TPC_inner} (GeV/c)");
  hMomTranslation->GetYaxis()->SetTitle("expected(dE/dx) (arb. unit)");
  
  delete hDummy;  
  
  // Maps for different number of clusters - used to extract dependence of ncl and create corresponding map for paramatrisation  
  TH3D* hPreMapClusterResolved[numClusterBins];
  
  hDummy = new THnSparseD("hDummy", "", nDim, thnBins, xmin, xmax);
  hDummy->SetBinEdges(1, binsPforMap);
  
  for (Int_t numClusters = clusterLowBound, clusterBin = 0; numClusters + clusterBinWidth <= clusterUpBound;
       numClusters += clusterBinWidth, clusterBin++)   {
    TH3D* hTemp = 0x0;
    hTemp = hDummy->Projection(0, 1, 2);
    hTemp->SetName("hTemp");
    hTemp->SetTitle(Form("#Delta' map for p_{TPC_Inner} vs. tan(#theta) for ncl around %d", numClusters));
    hTemp->GetXaxis()->SetTitle("tan(#theta)");
    hTemp->GetYaxis()->SetTitle("p_{TPC_inner} (GeV/c)");
    hTemp->GetZaxis()->SetTitle("#Delta'");
    
    hTemp->SetName(Form("hPreMapClusterResolved_%d", clusterBin));
    hPreMapClusterResolved[clusterBin] = hTemp;
  }
       
  delete hDummy;
  
  
  
  Long64_t nTreeEntries = tree->GetEntriesFast();
  
  Double_t dEdx = 0.; // Measured dE/dx
  Double_t dEdxExpected = 0.; // Expected dE/dx according to parametrisation
  Double_t tanTheta = 0.; // Tangens of (local) theta at TPC inner wall
  Double_t pTPC = 0.; // Momentum at TPC inner wall
  UShort_t tpcSignalN = 0; // Number of clusters used for dEdx
  //Int_t fMultiplicity = 0;
  UChar_t pidType = 0;
  //Double_t phiPrime = 0;
  
  
  tree->SetBranchAddress("dEdx", &dEdx);
  tree->SetBranchAddress("dEdxExpected", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  tree->SetBranchAddress("pTPC", &pTPC);
  tree->SetBranchAddress("pidType", &pidType);
  //tree->SetBranchAddress("phiPrime", &phiPrime);
  //if (isPbPb)
  //  tree->SetBranchAddress("fMultiplicity", &fMultiplicity);
  
  
  /*TODO NOW
  TF1* fShapeSmallP = new TF1("fShapeSmallP", "pol5", -0.4, 0.4);
  fShapeSmallP->SetParameters(1.01712, -0.0202725, -0.260692, 0.261623, 0.671854, -1.14014);
  */
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);
    
    if (dEdx <= 0 || dEdxExpected <= 0)
      continue;
    
    if (nclCut >= 0 && tpcSignalN <= nclCut)
      continue;

    
    //if (pidType == kV0idNoTOF || pidType == kV0idPlusTOFaccepted || pidType == kV0idPlusTOFrejected) continue; //For testing to get the old result
    
    // Use different PID techniques depending on the momentum - except for MC
    if (pidType != kMCid) {
      if ((pTPC < 0.6 && pidType != kTPCid)               // No V0s/TOF below 0.6 due to bad statistics
        //TODO NOW Maybe don't use this line since the tails get fewer V0's than the main peak, so the tails are overall reduced || (pTPC >= 0.6 && pTPC < pThresholdV0cut && pidType != kTPCandTOFid) || // Do not mix V0's and TOF in this range
          || (pTPC >= pThresholdV0cut && pTPC < pThresholdV0plusTOFcut &&
              pidType != kV0idNoTOF && pidType != kV0idPlusTOFaccepted) // Only V0's above pThresholdV0cut due to contamination
          || (pTPC >= pThresholdV0plusTOFcut && pidType != kV0idPlusTOFaccepted)) // Additionally require TOF
        continue;
    }
    
    //Double_t pT = pTPC*TMath::Sin(-TMath::ATan(tanTheta)+TMath::Pi()/2.0);
    //if ((phiPrime > 0.072/pT+TMath::Pi()/18.0-0.035 && phiPrime < 0.07/pT/pT+0.1/pT+TMath::Pi()/18.0+0.035)) 
    //  continue;
    
    /*TODO 
    if (TMath::Abs(tanTheta) <= 0.4) {
      Double_t p0 = fShapeSmallP->Eval(tanTheta) - 1.0; // Strength of the correction
      Double_t p1 = -9.0; // How fast the correction is turned off
      Double_t p2 = -0.209; // Turn off correction around 0.2 GeV/c
      Double_t p3 = 1.0; // Delta' for large p should be 1

      Double_t corrFactor = TMath::Erf((pTPC + p2) * p1) * p0 + p3 + p0; // Add p0 to have 1 for p3 = 1 and large pTPC
      dEdxExpected *= corrFactor;
    }*/
    
    // TODO Old unsuccessful try
    //Double_t thetaGlobalTPC = -TMath::ATan(tanTheta) + TMath::Pi() / 2.;
    //Double_t pTtpc = pTPC * TMath::Sin(thetaGlobalTPC);
    //Double_t pTtpcInv = (pTtpc > 0) ? 1. / pTtpc : 0;
    //Double_t p0 = 1.0;
    //Double_t p1 = 1./ 0.5;//TODO 2.0;
    //Double_t p2 = -0.05;//TODO 0.1
    //Double_t pTcorrFactor = p0 + (pTtpcInv > p1) * p2 * (pTtpcInv - p1);
    //
    //dEdxExpected *= pTcorrFactor;
    
    
    // Correct for pure momentum dependence (azimuthal dependence comes into play)
    /*TODO NOW
    Double_t correctionFactor = hPiMap->GetBinContent(getBinX(hPiMap, tanTheta), getBinY(hPiMap, pTPC));
    correctionFactor -= 1.;
    correctionFactor *= 0.5*(TMath::Erf((0.5 - pTPC) / 0.05) + 1.); // Only correct up to 0.5 GeV/c and switch off within about 2*0.05 GeV/c
    //correctionFactor *=  0.5* splPr->Eval(pTPC / massProton) / splPi->Eval(pTPC / massPion);//TODO NOW
    correctionFactor += 1.;
    dEdx /= correctionFactor;*/
    
    
    // TODO Maybe needed(?) - maybe not for the maps, since this should be valid for all particles!: "acceptedByPhiPrimeCut == 1"
    //if (!(TMath::Abs(tpcSignalN - 150) <= 5 && (phiPrime < 0.072/pT+pi/18.0-0.035 || phiPrime > 0.07/pT/pT+0.1/pT+pi/18.0+0.035)))
    //  continue
    
    hPreMap->Fill(tanTheta, pTPC, dEdx / dEdxExpected);
    hMomTranslation->Fill(pTPC, dEdxExpected);
    
    for (Int_t numClusters = clusterLowBound, clusterBin = 0; numClusters + clusterBinWidth <= clusterUpBound;
         numClusters += clusterBinWidth, clusterBin++)   {
      if (tpcSignalN > numClusters && tpcSignalN <= numClusters + clusterBinWidth) {
        hPreMapClusterResolved[clusterBin]->Fill(tanTheta, pTPC, dEdx / dEdxExpected);
        break;
      }
    }
        
    if (i % 1000 == 0)
      progressbar(100. * (((Double_t)i) / nTreeEntries));
  }
  
  progressbar(100.);
  
  
  // Go to MIP region (around p = 3 GeV/c, tan(theta) small) and take the mean value as the end point of the map
  // Since the contamination with kaons already starts, it is safer to go only to up to p = 2.7 GeV/c. Or: Use V0s
  printf("\nDetermining upper and lower 1/dEdx bound of map...\n\n");
  ///*TODO NOW NOW
  tree->Project("hAdjustMIP(200000,0,100)", "dEdxExpected", 
                Form("pTPC > 1.5 && dEdx > 0 && dEdxExpected > 0 %s", multiplicitySelectionString.Data()));

  TH1D* hAdjustMIP = (TH1D*)gDirectory->Get("hAdjustMIP");
  
  Int_t binMIP = hAdjustMIP->FindFirstBinAbove(0.0);
  
  Double_t tempMean = (binMIP > 0) ? hAdjustMIP->GetBinCenter(binMIP) : -1.;
  if (tempMean <= 0)  {
    printf("ERROR: Failed to determine expected dEdx in MIP region!\n");
    return 0x0;
  }
  Double_t upperMapBound = 1. / tempMean;
  
  TProfile* hMomTranslationProfile = hMomTranslation->ProfileX();
	hMomTranslationProfile->GetXaxis()->SetRangeUser(0.7, 4.5);
  Int_t mipBin = hMomTranslationProfile->GetMinimumBin();
  Double_t mipMomentum = hMomTranslationProfile->GetXaxis()->GetBinUpEdge(mipBin);
  
  printf("Upper 1/dEdx bound of map set to: 1./%f = %f (pTPC ~ %.2f GeV/c)\n", tempMean, upperMapBound, mipMomentum);
  
  delete hAdjustMIP;
  delete hMomTranslationProfile;
  
  // Take pre-map and loop through tanTheta-rows in momentum bins. Take the first row with sufficient statistics as
  // the lower bound and translate to the corresponding 1./dEdx
  TH2D* hPreMap_p_tanTheta = (TH2D*)hPreMap->Project3D("yx");
  Int_t lowMomentumBinWithSufficientStatistics = 1;
  for (Int_t binY = 1; binY <= hPreMap_p_tanTheta->GetNbinsY(); binY++) {
    Bool_t sufficientStatistics = kTRUE;
    for (Int_t binX = 1; binX <= hPreMap_p_tanTheta->GetNbinsX(); binX++) {
      // TODO Try to estimate the statistics with resolutionFactorY
      if (hPreMap_p_tanTheta->GetBinContent(binX, binY) < 100 / resolutionFactorY) {
        sufficientStatistics = kFALSE;
        break;
      }
    }
    if (sufficientStatistics) {
      lowMomentumBinWithSufficientStatistics = binY;
      break;
    }    
  }
  TH1D* hDeDxExpectedLowPbound = hMomTranslation->ProjectionY("_pyLow", lowMomentumBinWithSufficientStatistics, lowMomentumBinWithSufficientStatistics);
    
  Double_t meanDeDxExpectedLowPbound = hDeDxExpectedLowPbound->GetMean();  
  delete hDeDxExpectedLowPbound;
      
  if (meanDeDxExpectedLowPbound <= 0)  {
    printf("\nCould not determine mean dEdx_expected for lower p bound (pTPC = %f)!\n", 
           hPreMap_p_tanTheta->GetYaxis()->GetBinCenter(lowMomentumBinWithSufficientStatistics));
    return 0x0;
  }     
  Double_t lowerMapBound = 1. / meanDeDxExpectedLowPbound;
  printf("Lower 1/dEdx bound of map set to: 1./%f = %f (pTPC ~ %.2f GeV/c)\n", meanDeDxExpectedLowPbound, lowerMapBound, 
         hPreMap_p_tanTheta->GetYaxis()->GetBinCenter(lowMomentumBinWithSufficientStatistics));
  //*/
  
  /*
  printf("TODO: Fixed upper map bound to 0.02!!\n");//TODO NOW NOW
  Double_t upperMapBound = 0.02;//TODO NOW NOW
  printf("TODO: Fixed lower map bound to 0.002!!\n");//TODO NOW NOW
  Double_t lowerMapBound = 0.002;//TODO NOW NOW
  */
  //Double_t lowerMapBound = 0.0016;//TODO NOW NOW Should be used for PbPb due to bad splines
  //printf("\nTODO NOW: lowerMapBound fixed to 0.0016 - only needed for bad splines in PbPb!!!\n\n");//TODO NOW NOW
  
  
  // Binning was found to yield good results, if 40 bins are chosen for the range 0.0016 to 0.02. For the new variable range,
  // scale the number of bins correspondingly
  
  Int_t binsY = TMath::Nint((upperMapBound - lowerMapBound) / (0.02 - 0.0016) * 40);//WARNING: Also change in AliPIDResponse, if changed here!
  
  if (resolutionFactorY != 1)
    binsY /= resolutionFactorY;
  
  printf("Used number of y bins: %d\n", binsY);

  printf("\n\nExtracting histograms from tree. This may take some time...\n");
  
  TH2D* hRes3DprimeFit = new TH2D("hRes3DprimeFit", 
                                  "#Delta' map for 1./(dE/dx) vs. tan(#theta) via Gauss fit;tan(#theta);1./(dE/dx) (arb. unit);#Delta'",
                                  binsX, lowerTanTheta, upperTanTheta, binsY, lowerMapBound, upperMapBound);
  TH2D* hSigmaPar1 = (TH2D*)hRes3DprimeFit->Clone("hSigmaPar1");
  hSigmaPar1->SetTitle("Parameter 1 of #sigma_{rel}(#Delta') map for 1./(dE/dx) vs. tan(#theta) via Gauss fit");
  hSigmaPar1->GetZaxis()->SetTitle("Par1(#sigma_{rel}(#Delta'))");
  
  TH2D* hRes3DprimeProfile = (TH2D*)hRes3DprimeFit->Clone("hRes3DprimeProfile");
  hRes3DprimeProfile->SetTitle("#Delta' map for 1./(dE/dx) vs. tan(#theta) via mean");

  

  printf("\n\nStarted map creation....\n");
  
  
  const Int_t numTotalBins = hPreMap->GetXaxis()->GetNbins() * hPreMap->GetYaxis()->GetNbins();
  Bool_t specialHandlingApplied = kFALSE;
  
  // If special handling for the threshold of the TOF cut is requested: Determine the p-Bin which contains the threshold
  Int_t binWithThreshold = -1;
  if (pThresholdTOFcut > 0) {
    for (Int_t binY = 1; binY <= hPreMap->GetYaxis()->GetNbins(); binY++)  {
      if (hPreMap->GetYaxis()->GetBinLowEdge(binY) <= pThresholdTOFcut &&
          hPreMap->GetYaxis()->GetBinUpEdge(binY) > pThresholdTOFcut)  {
        binWithThreshold = binY;
        break;
      }
    }
  }
  
  for (Int_t binY = 1; binY <= hPreMap->GetYaxis()->GetNbins(); binY++) {
     // TODO added on 24.06.13
    // Do not go beyond MIP momentum
    if (hPreMap->GetYaxis()->GetBinLowEdge(binY) >= mipMomentum) {
      printf("\nReached momentum of MIPs. Stopping...\n");
      printedSomething = kTRUE;
      break;
    }
    
    // Use the expected dEdx to translate from pTPC to 1/dEdx   
    TH1D* hDeDxExpected = hMomTranslation->ProjectionY("_py", binY, binY);
    
    Double_t meanDeDxExpected = hDeDxExpected->GetMean();  
    delete hDeDxExpected;
      
    if (meanDeDxExpected <= 0)  {
      printf("\nCould not determine mean dEdx_expected for (p = %f) - skipping!\n", hPreMap->GetYaxis()->GetBinCenter(binY));
      printedSomething = kTRUE;
      continue;
    }     
    
    Int_t dEdxBin = hRes3DprimeFit->GetYaxis()->FindBin(1. / meanDeDxExpected);
    if (dEdxBin < 1 || dEdxBin > hRes3DprimeFit->GetYaxis()->GetNbins())  {
      printf("\nSkipped bin (out of map range): p = %f GeV/c, 1./<dEdxExpected> = %f\n", hPreMap->GetYaxis()->GetBinCenter(binY), 
             1. / meanDeDxExpected);
      printedSomething = kTRUE;
      continue;
    }
    
    // If there are subsequent bins with same 1/dEdxExpected-bin, merge them to increase statistics.
    Int_t nMergeBins = 1;
    for (; ; nMergeBins++) {
      if (binY + nMergeBins > hPreMap->GetYaxis()->GetNbins())  {
        nMergeBins--;
        break;
      }
      
      // TODO added on 24.06.13
      // Do not go beyond MIP momentum
      if (hPreMap->GetYaxis()->GetBinLowEdge(binY + nMergeBins) >= mipMomentum) {
        nMergeBins--;
        break;
      }
        
      hDeDxExpected = hMomTranslation->ProjectionY("_py", binY + nMergeBins, binY + nMergeBins);
      Double_t meanDeDxExpectedTemp = hDeDxExpected->GetMean();
      delete hDeDxExpected;
      
      if (meanDeDxExpectedTemp <= 0 || hRes3DprimeFit->GetYaxis()->FindBin(1. / meanDeDxExpectedTemp) != dEdxBin) {
        nMergeBins--;
        break;
      }
    }
    
    Int_t binYhigh = binY + nMergeBins;
    
    //if (binYhigh != binY) {
    //  printf("****Merged: %d for %f - %f\n", nMergeBins, hPreMap->GetYaxis()->GetBinLowEdge(binY),
    //         hPreMap->GetYaxis()->GetBinUpEdge(binYhigh));
    //}
    
    // Average over some bins to get rid of threshold effects around the cut
    //Double_t referenceY = hPreMap->GetYaxis()->GetBinLowEdge(binYhigh);
    //Double_t binWidthY = hPreMap->GetYaxis()->GetBinWidth(binYhigh);
    
    Int_t nMapRows = 1; // If no special handling: Only one row of the map filled

    if (pThresholdTOFcut > 0) {
      // TOF cut causes step down in statistics behind the treshold towards larger p. Therefore: Merge the requested number of bins
      // above threshold - 1 (or only a part of them, if due to the same 1/dEdx-Bin some bins already got merged above)
      if (binYhigh >= binWithThreshold - 1 && binYhigh <= binWithThreshold - 1 + nMergeBinsAroundThreshold - 1) {
        // Check, if the next momentum bin after the merging, i.e. binYhigh + nMergeBinsAroundThreshold + 1, is apart in the row
        // of the map. If so, fill all rows before (usually at most 1 row) with the same value.
        // Take the bin AFTER the last merged bin because this bin should have better statistics
        hDeDxExpected = hMomTranslation->ProjectionY("_py", binWithThreshold + nMergeBinsAroundThreshold + 1, 
                                                            binWithThreshold + nMergeBinsAroundThreshold + 1);
        Double_t meanDeDxExpectedTemp = hDeDxExpected->GetMean();
        delete hDeDxExpected;
      
        if (meanDeDxExpectedTemp <= 0)  {
          printf("\nError in special handling: Failed to determine mean dEdxExpected of next row after \"the gap\"!\n");
          printedSomething = kTRUE;
        }
        else  {
          // nMapRows should be at least 1. If the row minus 1 of the bin after the map coincides with the "usual" row, then the second
          // argument in the following call should yield 1. For sanity, take Max(1,x) instead of x, although this should make no
          // difference, if everything is fine.
          nMapRows = TMath::Max(1, (hRes3DprimeFit->GetYaxis()->FindBin(1. / meanDeDxExpectedTemp) - 1) - dEdxBin + 1);
        }
        
        Int_t nAdditionalMergeBins = binWithThreshold + nMergeBinsAroundThreshold - binYhigh;
        nMergeBins += nAdditionalMergeBins;
        binYhigh += nAdditionalMergeBins;
        
        // Special handling is only really applied, if bins get merged due to it
        if (nAdditionalMergeBins > 0) {
          if (specialHandlingApplied) {
            printf("\nWarning in special handling: It is applied twice, which should not happen by design!\n");
            printedSomething = kTRUE;
          }
          
          printf("\nAdditional merging due to special handling of threshold for TOF cut: %d p bins and %d 1/dEdx rows\n", 
                 nAdditionalMergeBins, nMapRows);
          printedSomething = kTRUE;
          
          specialHandlingApplied = kTRUE;
        }
      }
      /*OLD
      if (referenceY + nMergeBinsAroundThreshold / 4. * binWidthY >= pThresholdTOFcut && referenceY < pThresholdTOFcut) {
        // Check, if the next momentum bin after the merging, i.e. binYhigh + nMergeBinsAroundThreshold + 1, is appart in the row
        // of the map. If so, fill all rows before (usually at most 1 row) with the same value.
        // Take the bin AFTER the last merged bin because this bin should have better statistics
        hDeDxExpected = hMomTranslation->ProjectionY("_py", binYhigh + nMergeBinsAroundThreshold + 1, 
                                                            binYhigh + nMergeBinsAroundThreshold + 1);
        Double_t meanDeDxExpectedTemp = hDeDxExpected->GetMean();
        delete hDeDxExpected;
      
        if (meanDeDxExpectedTemp <= 0)  {
          printf("\nError in special handling: Failed to determine mean dEdxExpected of next row after \"the gap\"!\n");
          printedSomething = kTRUE;
        }
        else  {
          // nMapRows should be at least 1. If the row minus 1 of the bin after the map coincides with the "usual" row, then the second
          // argument in the following call should yield 1. For sanity, take Max(1,x) instead of x, although this should make no
          // difference, if everything is fine.
          nMapRows = TMath::Max(1, (hRes3DprimeFit->GetYaxis()->FindBin(1. / meanDeDxExpectedTemp) - 1) - dEdxBin + 1);
        }
        
        nMergeBins += nMergeBinsAroundThreshold;
        binYhigh += nMergeBinsAroundThreshold;
        specialHandlingApplied = kTRUE;
      }*/
    }
    
    for (Int_t binX = 1; binX <= hPreMap->GetXaxis()->GetNbins(); binX++) {  
      TH1D* hTempProjectionZ = 0x0;
      
      hTempProjectionZ = hPreMap->ProjectionZ("hTempProjectionZ", binX, binX, binY, binYhigh);
      
      if (hTempProjectionZ->GetEntries() < 10) {
        printf("\nWarning: Too few entries for (%f, %f - %f (->%f)): skipped\n",
               hPreMap->GetXaxis()->GetBinCenter(binX), hPreMap->GetYaxis()->GetBinLowEdge(binY),
               hPreMap->GetYaxis()->GetBinUpEdge(binYhigh), hRes3DprimeFit->GetYaxis()->GetBinCenter(dEdxBin)); 
        printedSomething = kTRUE;

        delete hTempProjectionZ;
        continue;
      }  
      
      Double_t meanWithoutFit = hTempProjectionZ->GetMean();

      // Do not overwrite with data from lower statistics (larger momentum) -> Most of the cases should be caught by the merging!
      for (Int_t i = 0; i < nMapRows; i++)  {
        if (hRes3DprimeProfile->GetBinContent(binX, dEdxBin + i) <= 0) { 
          hRes3DprimeProfile->SetBinContent(binX, dEdxBin + i, meanWithoutFit);
          hRes3DprimeProfile->SetBinError(binX, dEdxBin + i, (meanWithoutFit <= 0) ? 1000 : 1);
        }
        else  {
          printf("\nWarning: Data for mean (%f) already set (%f, %f (->%f) to %f - will NOT be overwritten\n",
                  meanWithoutFit, hPreMap->GetXaxis()->GetBinCenter(binX), 
                  hPreMap->GetYaxis()->GetBinCenter(binY), hRes3DprimeFit->GetYaxis()->GetBinCenter(dEdxBin + i),
                  hRes3DprimeProfile->GetBinContent(binX, dEdxBin + i));
          printedSomething = kTRUE;
        }
      }
      
      
      TFitResult* fitRes = 0x0;
      if (useDoubleGaussFit) {
        fitRes = doubleGaussFit(hTempProjectionZ, 0.5 *(hPreMap->GetYaxis()->GetBinCenter(binY) + hPreMap->GetYaxis()->GetBinCenter(binYhigh)),
                                splPr, splPi, "QNS") ;
      }
      else {
        Double_t lowThreshold = 0.6;
        Double_t highThreshold = 1.6;
        FindFitRange(hTempProjectionZ, lowThreshold, highThreshold);
        TFitResultPtr fitResPtr = hTempProjectionZ->Fit("gaus", "QNS", "", lowThreshold, highThreshold);
        fitRes = ((Int_t)fitResPtr == 0) ? (TFitResult*)fitResPtr.Get()->Clone() : 0x0;
      }
        
      Double_t mean = 0;
      Double_t meanError = 0;
                
      // If the fit failed, use the mean of the histogram as an approximation instead
      if (!fitRes) {
        mean = meanWithoutFit;
        meanError = 1000;
        
        printf("\nWarning: Fit failed for (%f, %f (->%f)), entries: %.0f -> Using mean instead (%f +- %f)\n",
               hPreMap->GetXaxis()->GetBinCenter(binX), hPreMap->GetYaxis()->GetBinCenter(binY), 
               hRes3DprimeFit->GetYaxis()->GetBinCenter(dEdxBin), hTempProjectionZ->GetEntries(), mean, meanError); 
        printedSomething = kTRUE;
      }
      else {
        mean = fitRes->GetParams()[1];
        meanError = fitRes->GetErrors()[1];
      }
      
      // Do not overwrite with data from lower statistics (larger momentum) -> Most of the cases should be caught by the merging!
      for (Int_t i = 0; i < nMapRows; i++)  {
        if (hRes3DprimeFit->GetBinContent(binX, dEdxBin + i) <= 0)  {
          hRes3DprimeFit->SetBinContent(binX, dEdxBin + i, mean);
          hRes3DprimeFit->SetBinError(binX, dEdxBin + i, meanError);
        }
        else  {
          printf("\nWarning: Data for fit (%f) already set (%f, %f (->%f) to %f - will NOT be overwritten\n",
                  mean, hPreMap->GetXaxis()->GetBinCenter(binX), 
                  hPreMap->GetYaxis()->GetBinCenter(binY), hRes3DprimeFit->GetYaxis()->GetBinCenter(dEdxBin + i),
                  hRes3DprimeFit->GetBinContent(binX, dEdxBin + i));
          printedSomething = kTRUE;
        }
      }
      
      delete fitRes;
      delete hTempProjectionZ;

      //progressbar(100. * (((Double_t)(((binY - 1.) * hPreMap->GetXaxis()->GetNbins()) + binX)) / numTotalBins)); continue; //TODO NOW NOW
      // Extract the dependence on the number of clusters and fill the map with the resolution parameter
      TFitResult* fNclStats = extractClusterDependence(&hPreMapClusterResolved[0], numClusterBins, clusterLowBound, clusterUpBound, 
                                                       binX, binY, binYhigh, c0, fSave, splPr, splPi);

      
      Double_t c1 = -1; 
      Double_t c1Error = 1000;

      if (fNclStats) {
        c1 = fNclStats->GetParams()[1];
        c1Error = fNclStats->GetErrors()[1];
      }

      // Do not overwrite with data from lower statistics (larger momentum) -> Most of the cases should be caught by the merging!
      for (Int_t i = 0; i < nMapRows; i++)  {
        if (hSigmaPar1->GetBinContent(binX, dEdxBin + i) <= 0)  {
          hSigmaPar1->SetBinContent(binX, dEdxBin + i, c1);
          hSigmaPar1->SetBinError(binX, dEdxBin + i, c1Error);
        }
        else  {
          printf("\nWarning: Data (sigma par 1) for fit (%f) already set (%f, %f (->%f) to %f - will NOT be overwritten\n",
                 c1, hPreMap->GetXaxis()->GetBinCenter(binX), hPreMap->GetYaxis()->GetBinCenter(binY),
                 hSigmaPar1->GetYaxis()->GetBinCenter(dEdxBin + i),
                 hSigmaPar1->GetBinContent(binX, dEdxBin + i));
          printedSomething = kTRUE;
        }
      }
      
      //delete fNclStats;
      
      
      progressbar(100. * (((Double_t)(((binY - 1.) * hPreMap->GetXaxis()->GetNbins()) + binX)) / numTotalBins));
    }
    
    binY += nMergeBins;
  }
  
  progressbar(100);
  printf("\n");
  
  /*OLD  
   // ---------------------------------------------------------------------------------
   // Use dEdx bins directly to create the map
   
   //upperMapBound usually = 0.02
   tree->Project(Form("hRes3Dprime(%d,-1,1, %d,0.004,%f, 200,0.6,1.6)", binsX, binsY, upperMapBound), 
   "dEdx/dEdxExpected:1./dEdx:tanTheta", "dEdx > 0 && dEdxExpected > 0 && pTPC < 3");
   //tree->Draw(Form("dEdx/dEdxExpected:1./dEdx:tanTheta>>hRes3Dprime(%d,-1,1, %d,0.004,%f, 200,0.6,1.6)", binsX, binsY, upperMapBound), 
   //           "dEdx > 0 && dEdxExpected > 0 && pTPC < 3");
   //tree->Draw("dEdx/dEdxExpected:1./dEdx:tanTheta>>hRes3Dprime(40,-1,1, 40,0.002,0.0225, 200,0.6,1.6)");// "acceptedByPhiPrimeCut == 0");
   //tree->Draw("dEdx/dEdxExpected:dEdx:tanTheta>>hRes3Dprime(30,-1,1, 100,45.,600., 200,0.6,1.6)");
   TH3F* hRes3Dprime = (TH3F*)gDirectory->Get("hRes3Dprime");
   TProfile2D* hTemp = hRes3Dprime->Project3DProfile("yx");
   hTemp->SetName("hTemp");
   TH2D* hRes3DprimeProfile = hTemp->ProjectionXY();
   hRes3DprimeProfile->SetName("hRes3DprimeProfile");
   normaliseHisto(hRes3DprimeProfile, -0.1, 0.1, kTypeDeltaPrime);
   hRes3DprimeProfile->GetZaxis()->SetRangeUser(0.95, 1.15);
   SetHistAxis(hRes3DprimeProfile, kTypeDeltaPrime);
   hRes3DprimeProfile->DrawCopy("surf1");
   canv3Dprime->cd(2);
   //OLD
   //hRes3Dprime->FitSlicesZ();
   //TH2D *hRes3DprimeFit = (TH2D*)gDirectory->Get("hRes3Dprime_1");
   
   
   TH2D *hRes3DprimeFit = (TH2D*)hRes3DprimeProfile->Clone("hRes3DprimeFit");
   hRes3DprimeFit->Reset();
   
   const Int_t numBinsSpecialHandling = 2;
   for (Int_t binY = 1; binY <= hRes3Dprime->GetYaxis()->GetNbins(); binY++) {
     Bool_t specialHandling = kFALSE;
     Double_t referenceY = hRes3Dprime->GetYaxis()->GetBinLowEdge(binY);
     Double_t binWidthY = hRes3Dprime->GetYaxis()->GetBinWidth(binY);
     
     // Average over some bins to get rid of threshold effects around the cut
     if (TOFcutDeDx > 0) {
       if (referenceY + numBinsSpecialHandling / 2. * binWidthY >= 1. / TOFcutDeDx && referenceY < 1. / TOFcutDeDx)
         specialHandling = kTRUE;
     }
     
     for (Int_t binX = 1; binX <= hRes3Dprime->GetXaxis()->GetNbins(); binX++) {  
       TH1D* hTempProjectionZ = 0x0;
       if (specialHandling)
         hTempProjectionZ = hRes3Dprime->ProjectionZ("hTempProjectionZ", binX, binX, binY, binY + numBinsSpecialHandling, "e");
       else
         hTempProjectionZ = hRes3Dprime->ProjectionZ("hTempProjectionZ", binX, binX, binY, binY, "e");
       
       Int_t binZ = 0;
       Int_t binZmax = 0;
       
       // Start from the very right side and go to the left -> Remove some of the first bins with non-zero content from the fitting.
       // This is important, since these bin might be affected from binning effects and spoil the gaussian shape.
       for (binZ = hTempProjectionZ->GetXaxis()->GetNbins(); binZ >= 1; binZ--) {
         if (hTempProjectionZ->GetBinContent(binZ) > 0)  {
           binZmax = binZ;
           
           break;
         }
       }
       
       binZmax = TMath::Max(1, binZmax - 2);
       
       TFitResultPtr fitRes = hTempProjectionZ->Fit("gaus", "QS", "", hTempProjectionZ->GetXaxis()->GetBinLowEdge(1),
       hTempProjectionZ->GetXaxis()->GetBinUpEdge(binZmax));
       Double_t mean = 0;
       
       // If the fit failed, use the mean of the histogram as an approximation instead
       if (((Int_t)fitRes) != 0) {
         printf("Warning: Fit failed for (%f, %f), entries: %.0f -> Using mean instead (%f)\n",
         hRes3DprimeProfile->GetXaxis()->GetBinCenter(binX), hRes3DprimeProfile->GetYaxis()->GetBinCenter(binY),
         fitRes->GetParams()[1], hTempProjectionZ->GetEntries(), hTempProjectionZ->GetMean()); 
         mean = hTempProjectionZ->GetMean();
       }
       else
         mean = fitRes->GetParams()[1];
       hRes3DprimeFit->SetBinContent(binX, binY, mean);
       
       if (specialHandling)  {
         for (Int_t i = 1; i <= numBinsSpecialHandling && i + binY <= hRes3DprimeFit->GetYaxis()->GetNbins(); i++)  
           hRes3DprimeFit->SetBinContent(binX, binY + i, mean);
       }
       
       
       delete hTempProjectionZ;
   }
   
   if (specialHandling)  {
     binY += numBinsSpecialHandling;
   }
   }
   
   normaliseHisto(hRes3DprimeFit, -0.1, 0.1, kTypeDeltaPrime);
   hRes3DprimeFit->GetZaxis()->SetRangeUser(0.95, 1.15);
   SetHistAxis(hRes3DprimeFit, kTypeDeltaPrime);
   hRes3DprimeFit->DrawCopy("surf1");
   
   delete hTemp;
   delete hRes3Dprime;
   */
  
  
  
  
  TH2D* hPureResults = (TH2D*)hRes3DprimeFit->Clone("hPureResults");
  TH2D* hPureResultsSigmaPar1 = (TH2D*)hSigmaPar1->Clone("hPureResultsSigmaPar1");
  
  printf("\n\nEliminating outliers....\n\n");
  eliminateOutliers(hRes3DprimeFit, outlierThreshold);
  
  printf("\n\nSame for the sigma map....\n\n");
  eliminateOutliers(hSigmaPar1, outlierThreshold); 
  
  TH2D* hRes3DprimeDiffSmooth = 0x0;
  
  if (doSmoothing) {
    // --------------------------------------------------
    // Determine difference between smoothed and pure map (gauss fit)    
    hRes3DprimeDiffSmooth = (TH2D*)hRes3DprimeFit->Clone("hRes3DprimeDiffSmooth");
    hRes3DprimeDiffSmooth->SetTitle("tan(#theta) and 1/(dE/dx) dependence of #Delta': Smooth(Gauss) - Gauss");
    
    // Smooth the histos
    
    // If the smoothing is to be applied, the bins with fit errors (content <= 0) must be eliminated in order not to bias the smoothing.
    // In case of no smoothing, the refinement will take care of such bins....
    eliminateNonPositivePointsInHisto(hRes3DprimeProfile);
    eliminateNonPositivePointsInHisto(hRes3DprimeFit);
    eliminateNonPositivePointsInHisto(hSigmaPar1);
    
    
    hRes3DprimeProfile->Smooth();
    hRes3DprimeFit->Smooth();
        
    hSigmaPar1->Smooth();
    
    hRes3DprimeDiffSmooth->Add(hRes3DprimeFit, -1);
    //hRes3DprimeDiffSmooth->GetZaxis()->SetRangeUser(-0.05, 0.05);
    SetHistAxis(hRes3DprimeDiffSmooth, kTypeDeltaPrime);
  }
  
  /* Normalisation is now done after the extrapolation
   *   // Use kNoNormalisation for final QA
   *   if (mapType == kSmallEtaNormalisation) {
   *     normaliseHisto(hRes3DprimeProfile, 0., 0.11, kTypeDeltaPrime);
   *     normaliseHisto(hRes3DprimeFit, 0., 0.11, kTypeDeltaPrime);
   }
   else if (mapType == kLargeEtaNormalisation) {
    normaliseHisto(hRes3DprimeProfile, 0.81, 0.99, kTypeDeltaPrime);
  normaliseHisto(hRes3DprimeFit, 0.81, 0.99, kTypeDeltaPrime);
   }
   */
  
  
  //hRes3DprimeProfile->GetZaxis()->SetRangeUser(0.95, 1.15);
  SetHistAxis(hRes3DprimeProfile, kTypeDeltaPrime);
  
  //hRes3DprimeFit->GetZaxis()->SetRangeUser(0.95, 1.15);
  SetHistAxis(hRes3DprimeFit, kTypeDeltaPrime);
  
  SetHistAxis(hSigmaPar1, kTypeSigmaDeltaPrime);
  
  
  if (saveResults)  {
    fSave->cd();
    hPreMap->Write();
    hMomTranslation->Write();
  }
  
  delete hPreMap;
  delete hMomTranslation;
  
  
  // --------------------------------------------------
  // Determine difference between gauss fit and pure mean
  TH2D* hRes3DprimeDiff = (TH2D*)hRes3DprimeFit->Clone("hRes3DprimeDiff");
  hRes3DprimeDiff->SetTitle("tan(#theta) and 1/(dE/dx) dependence of #Delta': Gauss - Profile");
  hRes3DprimeDiff->Add(hRes3DprimeProfile, -1);
  hRes3DprimeDiff->GetZaxis()->SetRangeUser(-0.05, 0.05);
  SetHistAxis(hRes3DprimeDiff, kTypeDeltaPrime);
  
  // --------------------------------------------------
  // Refine the map
  printf("\n\nMap created. Started to refine....\n\n");
  
  Double_t factorX = 6;
  Double_t factorY = 6;
  
  if (resolutionFactorX != 1)
    factorX *= resolutionFactorX;
  if (resolutionFactorY != 1)
    factorY *= resolutionFactorY;
  
  
  
  Bool_t skipBinsAtBoundaries = kFALSE; // Diese Bins vielleicht doch nicht rausnehmen. Der Anstieg in Sigma bei extrem kleinen Impulsen ist insofern richtig, da hier die Splines abweichen und somit eine große Streuung auftritt. Und bei großen Impulsen sollten im Moment die V0s Abhilfe von fehlerhaften Abweichungen schaffen. Beim Mean gilt dasselbe.
  TH2D* hRefinedNoNorm = refineHistoViaLinearExtrapolation(hRes3DprimeFit, factorX, factorY, kNoNormalisation, fSave, kFALSE,
                                                           skipBinsAtBoundaries);
  if (hRefinedNoNorm) {
    hRefinedNoNorm->SetName(Form("hRefined%s", suffixMapType[kNoNormalisation].Data()));
    SetHistAxis(hRefinedNoNorm, kTypeDeltaPrime);
  }
  
  TH2D* hRefinedSmallEtaNorm = refineHistoViaLinearExtrapolation(hRes3DprimeFit, factorX, factorY, kSmallEtaNormalisation, fSave, kFALSE,
                                                                 skipBinsAtBoundaries);
  if (hRefinedSmallEtaNorm)  {
    hRefinedSmallEtaNorm->SetName(Form("hRefined%s", suffixMapType[kSmallEtaNormalisation].Data()));
    SetHistAxis(hRefinedSmallEtaNorm, kTypeDeltaPrime);
  }
  
  TH2D* hRefinedLargeEtaNorm = refineHistoViaLinearExtrapolation(hRes3DprimeFit, factorX, factorY, kLargeEtaNormalisation, fSave, kFALSE,
                                                                 skipBinsAtBoundaries);
  if (hRefinedLargeEtaNorm) {
    hRefinedLargeEtaNorm->SetName(Form("hRefined%s", suffixMapType[kLargeEtaNormalisation].Data()));
    SetHistAxis(hRefinedLargeEtaNorm, kTypeDeltaPrime);
  }
  
  
  
  printf("\nSame for the sigma map.....\n\n");
  
  // Normalisation makes no sense for sigma map (will anyway be created for corrected data only)
  TH2D* hRefinedSigmaPar1 = refineHistoViaLinearExtrapolation(hSigmaPar1, factorX, factorY, kNoNormalisation, fSave, kTRUE,
                                                              skipBinsAtBoundaries);
  if (hRefinedSigmaPar1) {
    hRefinedSigmaPar1->SetName("hThetaMapSigmaPar1");
    SetHistAxis(hRefinedSigmaPar1, kTypeSigmaDeltaPrime);
  }
  
  
  printf("\nDone!\n\n");
  
  if (resolutionFactorX != 1 || resolutionFactorY != 1) {
    printf("\n***WARNING***: Resolution factor != 1 used for map creation! Resolution scaled by factor (%f, %f)\n\n", 
           resolutionFactorX, resolutionFactorY);
  }
  
  if (pThresholdTOFcut > 0) {
    printf("\n***INFO***: Requested special handling for map around p = %f GeV/c to take care of potential TOF cut around this value!", 
           pThresholdTOFcut);
    if (specialHandlingApplied) 
      printf(" -> Special handling was used!\n\n");
    else
      printf(" -> Special handling was not used (seemed to be not necessary)!\n\n");
  }
  else  {
    printf("\n***WARNING***: No special care taken of the TOF cut. Fine for MC, but not for data. Please check, if it is correct!\n\n");
  }
  
  printf("WARNING: Remember to switch on/off ncl cut w.r.t. geo cut!\n");
  
  if (saveResults)  {
    TNamed* c0Info = new TNamed("c0", Form("%f", c0));
    
    fSave->cd();
    hPureResults->Write();
    hRes3DprimeProfile->Write();
    hRes3DprimeFit->Write();
    hRes3DprimeDiff->Write();
    if (doSmoothing)
      hRes3DprimeDiffSmooth->Write();
    
    hPureResultsSigmaPar1->Write();
    hSigmaPar1->Write();
    
    if (hRefinedNoNorm)
      hRefinedNoNorm->Write();
    if (hRefinedSmallEtaNorm)
      hRefinedSmallEtaNorm->Write();
    if (hRefinedLargeEtaNorm)
      hRefinedLargeEtaNorm->Write();
    
    if (hRefinedSigmaPar1)
      hRefinedSigmaPar1->Write();
    
    c0Info->Write();
    
    TNamed* settings = new TNamed(
      Form("Settings for map: resolutionFactor (%f, %f),  pThresholdTOFcut %f, nMergeBinsAroundThreshold %d, pThresholdV0cut %f, pThresholdV0plusTOFcut %f, outlierThreshold %f, doSmoothing %d, c0 %f, nclCut %d, useDoubleGaussFit %d\n",
           resolutionFactorX, resolutionFactorY, pThresholdTOFcut, nMergeBinsAroundThreshold, pThresholdV0cut, pThresholdV0plusTOFcut,
           outlierThreshold, doSmoothing, c0, nclCut, useDoubleGaussFit), "");
    settings->Write();
    
    f->Close();
    fSave->Close();
    
    delete c0Info;
    
    gSystem->Exec(Form("echo \"%s\" >> %s/settingsForMapCreation.txt",
                       TString(settings->GetName()).ReplaceAll("Settings for map",
                                                               Form("Settings for map %s", saveFileName.Data())).Data(), path.Data()));
    return 0x0;
  }
  
  if (returnMapType == kSmallEtaNormalisation)
    return hRefinedSmallEtaNorm;
  else if (returnMapType == kLargeEtaNormalisation)
    return hRefinedLargeEtaNorm;
  else if (returnMapType == kNoNormalisation)
    return hRefinedNoNorm;
  
  printf("Map type %d is not supported. Returning 0x0...\n", returnMapType);
  
  return 0x0;
}
