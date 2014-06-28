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
const TString suffixMapType[kLargeEtaNormalisation + 1] = {"", "NoNormalisation", "SmallEtaNormalisation", "LargeEtaNormalisation"};

//__________________________________________________________________________________________
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
              
            values[binIndex] = h->GetBinContent(binX2, binY2);
            binIndex++;
          }
        }
        
        Double_t temp = getMedianOfNonZeros(values, nNeighbours);
        if (temp <= 0) {
          printf("Error: Could no eliminate values <= 0 for bin at (%f, %f)!\n",
                 h->GetXaxis()->GetBinCenter(binX), h->GetYaxis()->GetBinCenter(binY));
          temp = -1;
        }
        
        h->SetBinContent(binX, binY, temp);
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
    binError = 1000; // Should not happen because bins without content are rejected
    printf("ERROR: This should never happen: Trying to add bin in addPointToHyperplane with error not set....\n");
  }
  linExtrapolation->AddPoint(coord, h->GetBinContent(binX, binY, binError));
}

//__________________________________________________________________________________________
TH2D* refineHistoViaLinearExtrapolation(TH2D* h, Double_t refineFactorX, Double_t refineFactorY, Int_t mapType, 
                                        TFile* fSave, Bool_t skipBinsAtBoundaries = kFALSE)
{
  if (!h)
    return 0x0;
  
  Int_t nBinsX = h->GetXaxis()->GetNbins();
  Int_t nBinsY = h->GetYaxis()->GetNbins();

  Int_t firstYbin = 1;
  
  if (skipBinsAtBoundaries) {
    // Remove the two first and the last bin on y-axis
    nBinsY -= 3;
    firstYbin = 3; 
  }
    
  // Normalise map on demand
  
  // Use kNoNormalisation for final QA
  if (mapType == kSmallEtaNormalisation) {
    normaliseHisto(h, 0., 0.11, kTypeDeltaPrime);
  }
  else if (mapType == kLargeEtaNormalisation) {
    normaliseHisto(h, 0.81, 0.99, kTypeDeltaPrime);
  }
  
  // Interpolate to finer map
  TLinearFitter* linExtrapolation = new TLinearFitter(2, "hyp2", "");

  Int_t nBinsXrefined = nBinsX * refineFactorX;
  Int_t nBinsYrefined = nBinsY * refineFactorY; 
  
  TH2D* hRefined = new TH2D(Form("%s_%s_refined", h->GetName(), suffixMapType[mapType].Data()),  Form("%s (refined)",h->GetTitle()),
                            nBinsXrefined, h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(nBinsX),
                            nBinsYrefined, h->GetYaxis()->GetBinLowEdge(firstYbin), h->GetYaxis()->GetBinUpEdge(nBinsY));
  
  for (Int_t binX = 1; binX <= nBinsXrefined; binX++)  {
    for (Int_t binY = 1; binY <= nBinsYrefined; binY++)  {
      
      linExtrapolation->ClearPoints();
      
      hRefined->SetBinContent(binX, binY, 1); // Default value is 1
      
      Double_t centerX = hRefined->GetXaxis()->GetBinCenter(binX);
      Double_t centerY = hRefined->GetYaxis()->GetBinCenter(binY);
      
      // For interpolation: Just take the corresponding bin from the old histo.
      // For extrapolation: take the last available bin from the old histo.
      // If the boundaries are to be skipped, also skip the corresponding bins
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
      
      
      // Fit 2D-hyperplane
      if (linExtrapolation->GetNpoints() <= 0)
        continue;
        
      if (linExtrapolation->Eval() != 0)// EvalRobust -> Takes much, much, [...], much more time (~hours instead of seconds)
        continue;
      
      // Fill the bin of the refined histogram with the extrapolated value
      Double_t extrapolatedValue = linExtrapolation->GetParameter(0) + linExtrapolation->GetParameter(1) * centerX
                                 + linExtrapolation->GetParameter(2) * centerY;
                                 
      hRefined->SetBinContent(binX, binY, extrapolatedValue);      
    }
  } 
  
  delete linExtrapolation;
  
  return hRefined;
}


//__________________________________________________________________________________________
TFitResult*  doubleGaussFit(TH1D* h, Double_t currentMeanMomentum, TSpline3* splPr, TSpline3* splPi, TString fitOption = "QNS") 
{
  TFitResultPtr result = h->Fit("gaus", fitOption.Data());
  
  Double_t contaminationPeakMean = splPi->Eval(currentMeanMomentum / AliPID::ParticleMass(AliPID::kPion)) 
                                 / splPr->Eval(currentMeanMomentum / AliPID::ParticleMass(AliPID::kProton));

  if (contaminationPeakMean < h->GetXaxis()->GetBinLowEdge(1) ||
      contaminationPeakMean > h->GetXaxis()->GetBinUpEdge(h->GetNbinsX())) {
    return ((Int_t)result == 0) ? result.Get() : 0x0;
  }
  
  // Estimate parameters for doubleGauss fit
  Double_t estimatedMean = 0;
  Double_t estimatedSigma = 0;
  Double_t estimatedYield = 0;
  
  if ((Int_t) result == 0) {
    estimatedMean = result->GetParams()[1];
    estimatedSigma = result->GetParams()[2];
    estimatedYield = result->GetParams()[0];
  }
  else {
    estimatedMean = h->GetMean();
    estimatedSigma = h->GetRMS();
    estimatedYield = (estimatedSigma > 0) ? (h->Integral("width") / estimatedSigma) : h->GetEntries();
  }
  
  TF1* doubleGaus = new TF1("doubleGaus", "[0]*TMath::Gaus(x,[1],[2],0)+[3]*TMath::Gaus(x,[4],[2],0)", 0.6, 1.6);

  Double_t newPars[5] = { estimatedYield, estimatedMean, estimatedSigma, estimatedYield / 10., contaminationPeakMean };
  doubleGaus->SetParameters(newPars);
  doubleGaus->SetParLimits(0, 0., 100. * estimatedYield);
  doubleGaus->SetParLimits(1, 0.6, 1.6);
  doubleGaus->SetParLimits(2, 0, 100. * estimatedSigma);
  doubleGaus->SetParLimits(3, 0, 0.8 * estimatedYield);
  doubleGaus->SetParLimits(4, contaminationPeakMean - 0.1, contaminationPeakMean + 0.1);//TODO doubleGaus->SetParLimits(4, 0.6, 1.6);
  TFitResultPtr result2 = h->Fit(doubleGaus, fitOption.Data());

  // If fit failed, return results of standard fit instead
  if ((Int_t)result2 == 0) {
    return result2.Get();
  }
  
  return ((Int_t)result == 0) ? result.Get() : 0x0;
}


//__________________________________________________________________________________________
// NOTE: Use smoothing only for creation of sigma map. For normal eta map this pulls down the edges too much <-> trade-off between
// (normally very small) jumps in the map and good description at the edges.
// For the sigma map, there are not that sharp edges, but the fluctuations are by far stronger. Here it makes sense to use the smoothing.
TH2* checkShapeEtaTreePions(TTree* tree = 0x0, Double_t lowResolutionFactorX = 1/*2 or 2.1*/, Double_t lowResolutionFactorY = 1/*2 or 2.1*/, 
                       Bool_t doSmoothing = kFALSE, TString path = ".", TString suffixForFileName = "", Bool_t useDoubleGaussFit = kFALSE, 
                       TString pathNameThetaMap = "finalCuts/uncorrected/0.6GeVcut/pp/7TeV/LHC10d.pass2/outputCheckShapeEtaTree_2012_07_18__14_04.root", 
                       TString mapSuffix = "NoNormalisation",
                       TString pathNameSplinesFile = "DefaultSplines.root", TString prSplinesName = "TSPLINE3_DATA_PION_LHC10D_PASS2_PP_MEAN", 
                       Int_t returnMapType = kNoNormalisation, TString treeName = "fTreePions") 
{ 
  TString fileName = Form("bhess_PIDetaTreePions%s.root", suffixForFileName.Data());
  
  if (lowResolutionFactorX <= 0 || lowResolutionFactorY <= 0)  {
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
  
  TSpline3* splPr = 0x0;
  TSpline3* splPi = 0x0;
  
  // Extract the correction maps
  TFile* fPrMap = TFile::Open(pathNameThetaMap.Data());
  if (!fPrMap)  {
    std::cout << "Failed to open proton map file \"" << pathNameThetaMap.Data() << "\"!" << std::endl;
    return 0x0;
  }

  TH2D* hPrMap = 0x0;
  
  if (fPrMap) {
    hPrMap = dynamic_cast<TH2D*>(fPrMap->Get(Form("hRefined%s", mapSuffix.Data())));
    if (!hPrMap) {
      std::cout << "Failed to load proton theta map!" << std::endl;
      return 0x0;
    }
  }
  
  
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
    hPiMap = dynamic_cast<TH2D*>(fPiMap->Get("hRefinedSmallEtaNormalisation")); //TODO NOW NOW NOWForm("hRefined%s", mapSuffix.Data())));
    if (!hPiMap) {
      std::cout << "Failed to load pion momentum map!" << std::endl;
      return 0x0;
    }
  }
  
  // Output file
  
  TFile* fSave = 0x0;
  TString saveFileName = "";
  if (saveResults)  {
    TDatime daTime;
    
    TString saveSuffix = "";
    if (suffixForFileName.Length() > 0)
      saveSuffix = Form("_%s", suffixForFileName.Data());
    
    saveFileName = Form("outputCheckShapeEtaTreePions_%04d_%02d_%02d__%02d_%02d%s.root", daTime.GetYear(), daTime.GetMonth(),
                        daTime.GetDay(), daTime.GetHour(), daTime.GetMinute(), saveSuffix.Data());
    
    fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
    if (!fSave) {
      std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
      return 0x0;
    }
  }
  
  printf("\nProcessing file %s\n", f->GetName());
  
  
  // Only activate the branches of interest to save processing time
  tree->SetBranchStatus("*", 0); // Disable all branches
  tree->SetBranchStatus("pTPC", 1);
  //tree->SetBranchStatus("pT", 1);
  tree->SetBranchStatus("dEdx", 1);
  tree->SetBranchStatus("dEdxExpected", 1);
  tree->SetBranchStatus("tanTheta", 1);
  //if (isPbPb)
  //  tree->SetBranchStatus("fMultiplicity", 1);
  tree->SetBranchStatus("tpcSignalN", 1);
  
  TString multiplicitySelectionString = "";
// if (isPbPb)// Usually 20 bins from 0 to 4000 => bin width ~ 200
//    multiplicitySelectionString = "&& fMultiplicity > 3000 && fMultiplicity <= 3200";
  
  
  Int_t binsX = 30 / lowResolutionFactorX;
  Int_t binsY = 40 / lowResolutionFactorY;
  const Double_t pBoundLow = 0.15;// 1. / 0.6;
  const Double_t pBoundUp = 0.6;// 1. / 0.15;
  const Double_t lowerTanTheta = -1.0;
  const Double_t upperTanTheta = 1.0;
  
  progressbar(0.);
  
  TH2D* hMap = new TH2D("hMap", "#Delta' map for pTPC vs. tan(#theta) via Gauss fit;tan(#theta);p_{TPC} (GeV/c);#Delta'",
                        binsX, lowerTanTheta, upperTanTheta, binsY, pBoundLow, pBoundUp);
  
  const Int_t nDim = 3;
  
  Int_t thnBins[nDim] = {  binsX, binsY, 200  };
  Double_t xmin[nDim] = {  lowerTanTheta, pBoundLow, 0.6  };
  Double_t xmax[nDim] = {  upperTanTheta, pBoundUp, 1.6  };
  THnSparseD* hDummy = new THnSparseD("hDummy", "", nDim, thnBins, xmin, xmax);
  TH3D* hPreMap = hDummy->Projection(0, 1, 2);
  hPreMap->SetName("hPreMap");
  hPreMap->SetTitle("#Delta' map for p_{TPC_Inner} vs. tan(#theta)");   
  hPreMap->GetXaxis()->SetTitle("tan(#theta)");
  hPreMap->GetYaxis()->SetTitle("p_{TPC_inner} (GeV/c)");
  hPreMap->GetZaxis()->SetTitle("#Delta'");
  
  delete hDummy;
  
  Long64_t nTreeEntries = tree->GetEntriesFast();
  
  Double_t dEdx = 0.; // Measured dE/dx
  Double_t dEdxExpected = 0.; // Expected dE/dx according to parametrisation
  Double_t tanTheta = 0.; // Tangens of (local) theta at TPC inner wall
  Double_t pTPC = 0.; // Momentum at TPC inner wall
  Double_t pT = 0;
  UShort_t tpcSignalN = 0; // Number of clusters used for dEdx
  //Int_t fMultiplicity = 0;
  
  
  tree->SetBranchAddress("dEdx", &dEdx);
  tree->SetBranchAddress("dEdxExpected", &dEdxExpected);
  tree->SetBranchAddress("tanTheta", &tanTheta);
  tree->SetBranchAddress("tpcSignalN", &tpcSignalN);
  tree->SetBranchAddress("pTPC", &pTPC);
  //tree->SetBranchAddress("pT", &pT);
  //if (isPbPb)
  //  tree->SetBranchAddress("fMultiplicity", &fMultiplicity);
  
  Double_t correctionFactor = 1.0;
  Double_t dEdxOrig = 0;
  
  for (Long64_t i = 0; i < nTreeEntries; i++) {
    tree->GetEntry(i);
    
    if (dEdx <= 0 || dEdxExpected <= 0 || tpcSignalN <= 60) 
      continue;
    
    dEdxOrig = dEdx;
    
    // Correct for pure momentum dependence (azimuthal dependence comes into play)
    correctionFactor = hPiMap->GetBinContent(getBinX(hPiMap, tanTheta), getBinY(hPiMap, pTPC));
    correctionFactor -= 1.;
    correctionFactor *= 0.5*(TMath::Erf((0.5 - pTPC) / 0.05) + 1.); // Only correct up to 0.5 GeV/c and switch off within about 2*0.05 GeV/c
    correctionFactor += 1.;
    //TODO NOW dEdx /= correctionFactor;
    
  
    // Correct eta dependence
    correctionFactor = hPrMap->GetBinContent(getBinX(hPrMap, tanTheta), getBinY(hPrMap, 1./dEdxOrig));
    dEdx /= correctionFactor; 
  
    hPreMap->Fill(tanTheta, pTPC, dEdx / dEdxExpected);
    
    if (i % 1000 == 0)
    progressbar(100. * (((Double_t)i) / nTreeEntries));
  }
  
  progressbar(100.);

  

  printf("\n\nStarted map creation....\n");
  
  
  const Int_t numTotalBins = hPreMap->GetXaxis()->GetNbins() * hPreMap->GetYaxis()->GetNbins();
 
  for (Int_t binY = 1; binY <= hPreMap->GetYaxis()->GetNbins(); binY++) {
    for (Int_t binX = 1; binX <= hPreMap->GetXaxis()->GetNbins(); binX++) {  
      TH1D* hTempProjectionZ = 0x0;
      
      hTempProjectionZ = hPreMap->ProjectionZ("hTempProjectionZ", binX, binX, binY, binY);
      
      if (hTempProjectionZ->GetEntries() < 10) {
        printf("\nWarning: Too few entries for (%f, %f - %f): skipped\n",
               hPreMap->GetXaxis()->GetBinCenter(binX), hPreMap->GetYaxis()->GetBinLowEdge(binY),
               hPreMap->GetYaxis()->GetBinUpEdge(binY)); 
        printedSomething = kTRUE;

        delete hTempProjectionZ;
        continue;
      }  
      
      Double_t meanWithoutFit = hTempProjectionZ->GetMean();

      // Do not overwrite with data from lower statistics (larger momentum) -> Most of the cases should be caught by the merging!
      TFitResult* fitRes = 0x0;
      TFitResultPtr fitResPtr = hTempProjectionZ->Fit("gaus", "QSN", "");//TODO removed option L
      fitRes = ((Int_t)fitResPtr == 0) ? fitResPtr.Get() : 0x0;
      
      Double_t mean = 0;
      Double_t meanError = 0;
                
      // If the fit failed, use the mean of the histogram as an approximation instead
      if (!fitRes) {
        mean = meanWithoutFit;
        meanError = 1000;
        
        printf("\nWarning: Fit failed for (%f, %f), entries: %.0f -> Using mean instead (%f +- %f)\n",
               hPreMap->GetXaxis()->GetBinCenter(binX), hPreMap->GetYaxis()->GetBinCenter(binY), 
               hTempProjectionZ->GetEntries(), mean, meanError); 
        printedSomething = kTRUE;
      }
      else {
        mean = fitRes->GetParams()[1];
        meanError = fitRes->GetErrors()[1];
      }

      hMap->SetBinContent(binX, binY, mean);
      hMap->SetBinError(binX, binY, meanError);
      
      delete hTempProjectionZ;

      progressbar(100. * (((Double_t)(((binY - 1.) * hPreMap->GetXaxis()->GetNbins()) + binX)) / numTotalBins)); 
    }
  }
  
  progressbar(100);
  
  
  
  TH2D* hPureResults = (TH2D*)hMap->Clone("hPureResults");  
  TH2D* hRes3DprimeDiffSmooth = 0x0;
  
  if (doSmoothing) {
    // --------------------------------------------------
    // Determine difference between smoothed and pure map (gauss fit)    
    hRes3DprimeDiffSmooth = (TH2D*)hMap->Clone("hRes3DprimeDiffSmooth");
    hRes3DprimeDiffSmooth->SetTitle("tan(#theta) and pTPC dependence of #Delta': Smooth(Gauss) - Gauss");
    
    // Smooth the histos
    
    // If the smoothing is to be applied, the bins with fit errors (content <= 0) must be eliminated in order not to bias the smoothing.
    // In case of no smoothing, the refinement will take care of such bins....
    eliminateNonPositivePointsInHisto(hMap);
    
    
    hMap->Smooth();
        
    
    hRes3DprimeDiffSmooth->Add(hMap, -1);
    //hRes3DprimeDiffSmooth->GetZaxis()->SetRangeUser(-0.05, 0.05);
    SetHistAxis(hRes3DprimeDiffSmooth, kTypeDeltaPrime);
  }
  
  //hMap->GetZaxis()->SetRangeUser(0.95, 1.15);
  SetHistAxis(hMap, kTypeDeltaPrime);
    
  
  if (saveResults)  {
    fSave->cd();
    hPreMap->Write();
  }
  
  delete hPreMap;
  
  
  // --------------------------------------------------
  // Refine the map
  printf("\n\nMap created. Started to refine....\n\n");
  
  Double_t factorX = 6;
  Double_t factorY = 6;
  
  if (lowResolutionFactorX != 1)
    factorX *= lowResolutionFactorX;
  if (lowResolutionFactorY != 1)
    factorY *= lowResolutionFactorY;
  
  
  
  Bool_t skipBinsAtBoundaries = kFALSE; // TODO NOW Diese Bins vielleicht doch nicht rausnehmen. Der Anstieg in Sigma bei extrem kleinen Impulsen ist insofern richtig, da hier die Splines abweichen und somit eine große Streuung auftritt. Und bei großen Impulsen sollten im Moment die V0s Abhilfe von fehlerhaftenAbweichungen schaffen. Beim Mean gilt dasselbe.
  TH2D* hRefinedNoNorm = refineHistoViaLinearExtrapolation(hMap, factorX, factorY, kNoNormalisation, fSave, skipBinsAtBoundaries);
  hRefinedNoNorm->SetName(Form("hRefined%s", suffixMapType[kNoNormalisation].Data()));
  SetHistAxis(hRefinedNoNorm, kTypeDeltaPrime);
  
  TH2D* hRefinedSmallEtaNorm = refineHistoViaLinearExtrapolation(hMap, factorX, factorY, kSmallEtaNormalisation, fSave, skipBinsAtBoundaries);
  hRefinedSmallEtaNorm->SetName(Form("hRefined%s", suffixMapType[kSmallEtaNormalisation].Data()));
  SetHistAxis(hRefinedSmallEtaNorm, kTypeDeltaPrime);
  
  TH2D* hRefinedLargeEtaNorm = refineHistoViaLinearExtrapolation(hMap, factorX, factorY, kLargeEtaNormalisation, fSave, skipBinsAtBoundaries);
  hRefinedLargeEtaNorm->SetName(Form("hRefined%s", suffixMapType[kLargeEtaNormalisation].Data()));
  SetHistAxis(hRefinedLargeEtaNorm, kTypeDeltaPrime);
  
  printf("\nDone!\n\n");
  
  if (lowResolutionFactorX != 1 || lowResolutionFactorY != 1) {
    printf("\n***WARNING***: Low resolution used for map creation! Resolution scaled down by factor (%f, %f)\n\n", 
           lowResolutionFactorX, lowResolutionFactorY);
  }
  
  
  
  TH2D* hPureResultsSmallEtaNorm = (TH2D*)hPureResults->Clone("hPureResultsSmallEtaNorm");  
  normaliseHisto(hPureResultsSmallEtaNorm, 0., 0.11, kTypeDeltaPrime);
  
  if (saveResults)  {    
    fSave->cd();
    hPureResults->Write();
    hPureResultsSmallEtaNorm->Write();
    hMap->Write();
    if (doSmoothing)
      hRes3DprimeDiffSmooth->Write();
    
    hRefinedNoNorm->Write();
    hRefinedSmallEtaNorm->Write();
    hRefinedLargeEtaNorm->Write();
    
    TNamed* settings = new TNamed(
      Form("Settings for map: lowResolutionFactor (%f, %f), doSmoothing %d", lowResolutionFactorX, lowResolutionFactorY, doSmoothing), "");
    settings->Write();
    
    f->Close();
    fSave->Close();
        
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
