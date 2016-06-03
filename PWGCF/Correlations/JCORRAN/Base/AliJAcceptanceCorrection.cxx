/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// A class for calculating acceptance correction

#include "AliJAcceptanceCorrection.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>
#include <TF1.h>
#include "AliJAcceptanceFunctions.h"

/*
 * Default constructor.
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection() :
    fCard(),
    fDEtaNearAcceptance(),
    fDEtaDPhiNearAcceptance(),
    fDEtaDPhi3DNearAcceptance(),
    fMinCountsPerBinInclusive(1000),
    fDEtaNearLoaded(false),
    fDEtaDPhiNearLoaded(false),
    fDEtaDPhi3DNearLoaded(false),
    fLeadingParticleCorrelation(true)
{
  // default constructor
  Generate3DAcceptanceCorrection();
}

/*
 * Constructor that should be used for correct functioning of the class.
 * Requires definition of the JCard.
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection(AliJCard *inputCard) :
    fCard(inputCard),
    fDEtaNearAcceptance(),
    fDEtaDPhiNearAcceptance(),
    fDEtaDPhi3DNearAcceptance(),
    fMinCountsPerBinInclusive(1000),
    fDEtaNearLoaded(false),
    fDEtaDPhiNearLoaded(false),
    fDEtaDPhi3DNearLoaded(false),
    fLeadingParticleCorrelation(true)
{
  // Constructor with JCard
  Generate3DAcceptanceCorrection();
}

/*
 * Copy constructor
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection(const AliJAcceptanceCorrection& a):
    fCard(a.fCard),
    fDEtaNearAcceptance(a.fDEtaNearAcceptance),
    fDEtaDPhiNearAcceptance(a.fDEtaDPhiNearAcceptance),
    fDEtaDPhi3DNearAcceptance(a.fDEtaDPhi3DNearAcceptance),
    fDEtaDPhi3DNearAcceptanceCalculation(a.fDEtaDPhi3DNearAcceptanceCalculation),
    fMinCountsPerBinInclusive(a.fMinCountsPerBinInclusive),
    fDEtaNearLoaded(a.fDEtaNearLoaded),
    fDEtaDPhiNearLoaded(a.fDEtaDPhiNearLoaded),
    fDEtaDPhi3DNearLoaded(a.fDEtaDPhi3DNearLoaded),
    fLeadingParticleCorrelation(a.fLeadingParticleCorrelation)
{
  //copy constructor
}

/*
 * Destructor
 */
AliJAcceptanceCorrection::~AliJAcceptanceCorrection(){
  // destructor
  
  delete fDEtaDPhi3DNearAcceptanceCalculation;
  
}

/*
 * Equal sign operator
 */
AliJAcceptanceCorrection&  AliJAcceptanceCorrection::operator=(const AliJAcceptanceCorrection& a){
  //operator = 
  if(this != &a){
    fCard = a.fCard;
    fDEtaNearAcceptance = a.fDEtaNearAcceptance;
    fDEtaDPhiNearAcceptance = a.fDEtaDPhiNearAcceptance;
    fDEtaDPhi3DNearAcceptance = a.fDEtaDPhi3DNearAcceptance;
    fDEtaDPhi3DNearAcceptanceCalculation = a.fDEtaDPhi3DNearAcceptanceCalculation;
    fMinCountsPerBinInclusive = a.fMinCountsPerBinInclusive;
    fDEtaNearLoaded = a.fDEtaNearLoaded;
    fDEtaDPhiNearLoaded = a.fDEtaDPhiNearLoaded;
    fDEtaDPhi3DNearLoaded = a.fDEtaDPhi3DNearLoaded;
    fLeadingParticleCorrelation = a.fLeadingParticleCorrelation;
  }
  return *this;
}

/*
 * Method for generating 3D near side acceptance correction
 */
void AliJAcceptanceCorrection::Generate3DAcceptanceCorrection(){
  // Generate 3D near side acceptance correction
  
  // Read the eta range from card
  double etaRange = fCard->Get("EtaRange");
  
  // Define the function for calculating the acceptance correction
  AliJAcceptanceFunctions *functionWrapper = new AliJAcceptanceFunctions();
  TF1 *acceptanceCorrectionFunction = new TF1("acceptanceCorrectionFunction",functionWrapper,&AliJAcceptanceFunctions::AcceptanceCorrection3DNearSide,-2*etaRange,2*etaRange,2,"AliJAcceptanceFunctions","AcceptanceCorrection3DNearSide");
  acceptanceCorrectionFunction->SetParameter(1,etaRange);
  
  // Create the two dimensional histogram to store the values
  // This needs to have a lot of bins to have a decent accuracy
  fDEtaDPhi3DNearAcceptanceCalculation = new TH2D("acceptanceCalculation","acceptanceCalculation",800*etaRange,-2*etaRange,2*etaRange,1280,-TMath::Pi(),TMath::Pi());
  
  // Fill the histogram with values from the function
  // The reason why the function is not used directly is the time it takes to make numerical
  // solving needed to get a number out of the function.
  double binCenterX;
  double binCenterY;
  double binContent;
  int nBinsX = fDEtaDPhi3DNearAcceptanceCalculation->GetNbinsX();
  int nBinsY = fDEtaDPhi3DNearAcceptanceCalculation->GetNbinsY();
  for(int binX = 1; binX <= nBinsX; binX++){
    for(int binY = 1; binY <= nBinsY; binY++){
      binCenterX = fDEtaDPhi3DNearAcceptanceCalculation->GetXaxis()->GetBinCenter(binX);
      binCenterY = fDEtaDPhi3DNearAcceptanceCalculation->GetYaxis()->GetBinCenter(binY);
      acceptanceCorrectionFunction->SetParameter(0,binCenterY);
      binContent = acceptanceCorrectionFunction->Eval(binCenterX);
      fDEtaDPhi3DNearAcceptanceCalculation->SetBinContent(binX,binY,binContent);
    }
  }
  
  delete functionWrapper;
  delete acceptanceCorrectionFunction;
  
}

/*
 *  Method for reading acceptance correction histograms from the file.
 *  This method tries to read histograms for all possible corrections.
 *  If some of the histograms are missing, histogram manager gives an
 *  empty histogram at their place.
 *
 *  const char* fileName = file name for the file, in which the histogram are located
 */
void AliJAcceptanceCorrection::ReadMixedEventHistograms(const char *fileName){
  // read the mixed event histograms
  
  TPMERegexp sep("::");
  int ncol = sep.Split( fileName );
  TString filename = sep[0];
  
  if (TString(fileName).BeginsWith("alien:"))  TGrid::Connect("alien:");
  TFile *inclusiveFile = TFile::Open(filename);
  TDirectory * dir =  (TDirectory*) inclusiveFile;
  if( ncol > 1 ) dir = (TDirectory*)( inclusiveFile->Get(sep[1]));
  if( !dir ) {
    std::cout << "Wrong file name or directory name given for AliJAcceptanceCorrection::ReadMixedEventHistograms!" << std::endl;
    std::cout << "Inclusive histograms are not loaded!" << std::endl;
    return;
  }

  AliJHistManager *histogramReader = new AliJHistManager("hst",sep[1]);
  histogramReader->LoadConfig();
  
  // Check if the new acceptance histogram exists in the inclusive file.
  // This check is done for backwards compatibility. If the histogram is
  // found, flag this so we know it in fursther algorithms.
  if(histogramReader->HistogramExists("hDEtaNearMixAcceptance")){
    fDEtaNearAcceptance = histogramReader->GetTH1D("hDEtaNearMixAcceptance");
    NormalizeAcceptanceTraditional(fDEtaNearAcceptance, kAssocType);
    fDEtaNearLoaded = true;
  } else {
    std::cout << "Could not find histogram: hDEtaNearMixAcceptance" << std::endl;
    std::cout << "Inclusive histograms for traditional 1D acceptance correction are not loaded! " << std::endl;
  }
  
  // Load the mixed event deltaPhi deltaEta histogram in the traditional near side
  if(histogramReader->HistogramExists("hDphiDetaPta")){
    fDEtaDPhiNearAcceptance = histogramReader->GetTH2D("hDphiDetaPta");
    NormalizeAcceptanceTraditionalInclusive(fDEtaDPhiNearAcceptance, kAssocType);
    fDEtaDPhiNearLoaded = true;
  } else {
    std::cout << "Could not find histogram: hDphiDetaPta" << std::endl;
    std::cout << "Inclusive histograms for traditional 2D acceptance correction are not loaded! " << std::endl;
  }
  
  // Load the mixed event deltaPhi deltaEta histogram in the 3D near side
  if(histogramReader->HistogramExists("hDphiDetaXlong")){
    fDEtaDPhi3DNearAcceptance = histogramReader->GetTH2D("hDphiDetaXlong");
    NormalizeAcceptance3DNearSideInclusive(fDEtaDPhi3DNearAcceptance, kXeType);
    fDEtaDPhi3DNearLoaded = true;
  } else {
    std::cout << "Could not find histogram: hDphiDetaXlong" << std::endl;
    std::cout << "Inclusive histograms for 3D near side acceptance correction are not loaded! " << std::endl;
  }
  
  delete histogramReader;
  
}

/*
 * Traditional method for normalizing the acceptance histogram
 *
 * In this method the histogram is first rebinned and then normalized to interval [0,1]
 *
 *  AliJTH1D &acceptanceHisto = Custom histogram array containing the histogram to be normalized
 *  corrType assocType = Associated particle binning type, mainly pTa or xLong
 */
void AliJAcceptanceCorrection::NormalizeAcceptanceTraditional(AliJTH1D &acceptanceHisto, corrType assocType){
  // Method for normalizing and rebinning the inclusive acceptance histograms
  
  // Find the correct binning
  int numCent  = fCard->GetNoOfBins(kCentrType);
  int numPtt   = fCard->GetNoOfBins(kTriggType);
  int numAssoc = fCard->GetNoOfBins(assocType);
  
  // Variables for the loop
  int rebin;
  double counts;
  double maxValue;
  
  // Get the number of bins in the histograms
  const int nBins = acceptanceHisto[0][0][0]->GetNbinsX();
  
  // Loop over the input histograms and find the correct normalization
  for (int iCent = 0; iCent < numCent; iCent++) {
    for (int iPtt = 0; iPtt < numPtt; iPtt++){
      for (int iAssoc = 0; iAssoc < numAssoc; iAssoc++){
        
        // Do rebinning requiring at least fMinCountsPerBinInclusive counts in every bin in the histogram.  Maximum rebin is 16
        counts  = acceptanceHisto[iCent][iPtt][iAssoc]->Integral();
        rebin = GetRebin(counts,nBins,1);
        acceptanceHisto[iCent][iPtt][iAssoc]->Rebin(rebin);
        maxValue = acceptanceHisto[iCent][iPtt][iAssoc]->GetMaximum();
        if(maxValue > 0) acceptanceHisto[iCent][iPtt][iAssoc]->Scale(1.0/maxValue);
        
      }
    }
  }
}

/*
 * Do a proper normalization for the acceptance histograms
 * pTa or xLong bins are summed over to remove possible biases and to get better statistics
 * The  histograms arethen scaled to a given peak value
 *
 *  AliJTH2D &acceptanceHisto = Custom histogram array containing the histogram to be normalized
 *  corrType assocType = Associated particle binning type, mainly pTa or xLong
 */
void AliJAcceptanceCorrection::NormalizeAcceptanceInclusive(AliJTH2D &acceptanceHisto, corrType assocType, double peakValue){
  // Method for normalizing and rebinning the inclusive acceptance histograms from 3D near side
  
  // Find the correct binning
  const int numCent    = fCard->GetNoOfBins(kCentrType);
  const int numPtt     = fCard->GetNoOfBins(kTriggType);
  const int numAssoc   = fCard->GetNoOfBins(assocType);
  const int numZvertex = fCard->GetNoOfBins(kZVertType);
  
  // Loop over the input histograms and find the correct normalization
  for (int iCent = 0; iCent < numCent; iCent++) {
    for (int iPtt = 0; iPtt < numPtt; iPtt++){
      for (int iZVertex = 0; iZVertex < numZvertex; iZVertex++){
        
        // xLong bins bias the deltaEta deltaPhi distribution, so we need to integrate over them
        // Also pTa bins are integrated over to get a good statistics
        for (int iAssoc = 1; iAssoc < numAssoc; iAssoc++){
          // Check the leading particle condition for the histograms
          if(assocType == kAssocType && fCard->Get(kTriggType,iPtt) < fCard->Get(kAssocType,iAssoc) && fLeadingParticleCorrelation) continue;
          acceptanceHisto[1][iCent][iZVertex][iPtt][0]->Add(acceptanceHisto[1][iCent][iZVertex][iPtt][iAssoc]);
        }
        
        // Sum over z-vertex bins and put the histograms on first array index 0
        if(iZVertex == 0) acceptanceHisto[0][iCent][0][iPtt][1]->Reset();
        acceptanceHisto[0][iCent][0][iPtt][1]->Add(acceptanceHisto[1][iCent][iZVertex][iPtt][0]);
        
        // Rebin and normalize the acceptance histograms to the interval [0,peakValue]
        RebinAndNormalize(acceptanceHisto[1][iCent][iZVertex][iPtt][0],peakValue);
        
        // Note: The correction for away side effect for 3D near side done in GetAcceptanceCorrection3DNearSidelInclusive method
      } // z-vertex bins
      
      // Rebin and normalize the acceptance histograms integrated over z-vertices to the interval [0,peakValue]
      RebinAndNormalize(acceptanceHisto[0][iCent][0][iPtt][1],peakValue);
      
    } // Trigger pT bins
  } // Centrality bins
}


/*
 * Do a proper normalization for the acceptance histograms for traditional near side
 * These histograms we can just normalize to interval [0,1]
 *
 *  AliJTH2D &acceptanceHisto = Custom histogram array containing the histogram to be normalized
 *  corrType assocType = Associated particle binning type, mainly pTa or xLong
 */
void AliJAcceptanceCorrection::NormalizeAcceptanceTraditionalInclusive(AliJTH2D &acceptanceHisto, corrType assocType){
  // Method for normalizing and rebinning the inclusive acceptance histograms from 3D near side
  
  NormalizeAcceptanceInclusive(acceptanceHisto,assocType,1);
}

/*
 * Do a proper normalization for the acceptance histograms for 3D near side
 * The idea is to first scale the histogram to the length of constant deltaEta lines
 * Then the 3D away side effect is corrected by a fraction of pairs that are left
 * outside of the acceptance.
 *
 *  AliJTH2D &acceptanceHisto = Custom histogram array containing the histogram to be normalized
 *  corrType assocType = Associated particle binning type, mainly pTa or xLong
 */
void AliJAcceptanceCorrection::NormalizeAcceptance3DNearSideInclusive(AliJTH2D &acceptanceHisto, corrType assocType){
  // Method for normalizing and rebinning the inclusive acceptance histograms from 3D near side
  
  // Find the acceptance eta range
  const double etaRange = fCard->Get("EtaRange");

  // Do the normalization
  const double peakValue = 2*sqrt(2)*etaRange;
  NormalizeAcceptanceInclusive(acceptanceHisto,assocType,peakValue);
  
}

/*
 * Calculate the acceptance correction assuming triangular shape in deltaEta
 *
 *  double deltaEta = deltaEta for the particle pair
 *
 *  return = Acceptance correction based on input deltaEta
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTriangle(double deltaEta){
  // calculate acceptance correction on pseudorapidity triangle
  
  double absDEta = fabs(deltaEta);
  double etaRange = fCard->Get("EtaRange");

  // Triangle, f(0) = 1; f(etaRange) = 0 = f(-etaRange)
  double denominator = 1 - absDEta/(2*etaRange);

  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
}

/*
 * Return a calculated acceptance correction for 3D near side
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection3DNearSideCalculation(double deltaEta, double deltaPhi){
  // return the acceptance correction from the pre-calculated surface

  double denominator = fDEtaDPhi3DNearAcceptanceCalculation->GetBinContent(fDEtaDPhi3DNearAcceptanceCalculation->FindBin(deltaEta,deltaPhi));
  
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
}

/*
 * Calculate the acceptance correction from input histogram given in ReadInclusiveHistos
 *
 *  double deltaEta = deltaEta for the particle pair
 *  int centralityBin = bin index for centrality
 *  int triggerBin = bin index for trigger particle transverse momentum
 *  int assocBin = bin index for associated particle binning
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTriangleInclusive(double deltaEta, int centralityBin, int triggerBin, int assocBin){
  // Inclusive acceptance correction
  
  // If the inclusive histograms are not found from the file, return correction from triangle
  if(!fDEtaNearLoaded) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // If the given bin is negative, return correction from triangle
  if(assocBin < 0) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // The acceptance histogram comes always from pTa bins since 1D correction for 3D near side is meaningless
  TH1D *acceptanceHistogram = fDEtaNearAcceptance[centralityBin][triggerBin][assocBin];
  
  // Get the number of bins in the histograms
  const int nBins = acceptanceHistogram->GetNbinsX();
  
  // If there is less than defined number of entries per bin, just use triangle instead of inclusive
  if(acceptanceHistogram->GetEntries() < nBins*fMinCountsPerBinInclusive) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // Use the value in the bin corresponding to deltaEta as acceptance correction
  int bin =  acceptanceHistogram->FindBin(deltaEta);
  double denominator  =  acceptanceHistogram->GetBinContent(bin);
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
  
}

/*
 * Calculate the acceptance correction from 2D input histogram given in ReadInclusiveHistos
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int zVertexBin = z-vertex bin
 *  int triggerBin = trigger pT bin
 *  int firstBin = 0 for z-vertex summed correction, 1 for z-vertex binned correction
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTraditionalInclusiveBin(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin, int firstBin){
  // Inclusive acceptance correction from two dimensional histogram
  
  // If the inclusive histograms are not found from the file, return correction from triangle
  if(!fDEtaDPhiNearLoaded) return GetAcceptanceCorrectionTriangle(deltaEta);

  // The acceptance histogram comes always from 2D distribution
  // If z-vertex bin not specified in the argument list, return sum over z-vertex bins
  TH2D *acceptanceHistogram = fDEtaDPhiNearAcceptance[firstBin][centralityBin][zVertexBin][triggerBin][1-firstBin];
  
  // If there are less than fMinCountsPerBinInclusive entries per bin use calculation instead of histogram
  const int nBinsEta = acceptanceHistogram->GetNbinsX();
  const int nBinsPhi = acceptanceHistogram->GetNbinsY();
  const int nBins = nBinsEta*nBinsPhi;
  
  if(acceptanceHistogram->GetEntries() < nBins*fMinCountsPerBinInclusive) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // Use the value in the bin corresponding to deltaEta as acceptance correction
  double denominator  =  acceptanceHistogram->GetBinContent(acceptanceHistogram->FindBin(deltaEta,deltaPhi));
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
  
}

/*
 * Return a calculated acceptance correction for 3D near side. The value from the histogram
 * needs to be corrected for the away side effect before returning the correction
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int triggerBin = trigger pT bin
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTraditionalInclusive(double deltaEta, double deltaPhi, int centralityBin, int triggerBin){
  return GetAcceptanceCorrectionTraditionalInclusiveBin(deltaEta, deltaPhi, centralityBin, 0, triggerBin, 0);
}

/*
 * Return a calculated acceptance correction for 3D near side. The value from the histogram
 * needs to be corrected for the away side effect before returning the correction
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int zVertexBin = z-vertex bin
 *  int triggerBin = trigger pT bin
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTraditionalInclusive(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin){
  return GetAcceptanceCorrectionTraditionalInclusiveBin(deltaEta, deltaPhi, centralityBin, zVertexBin, triggerBin, 1);
}

/*
 * Return a calculated acceptance correction for 3D near side. The value from the histogram
 * needs to be corrected for the away side effect before returning the correction
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int zVertexBin = z-vertex bin
 *  int triggerBin = trigger pT bin
 *  int firstBin = 0 for z-vertex summed correction, 1 for z-vertex binned correction
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection3DNearSideInclusiveBin(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin, int firstBin){
  // Return the acceptance correction found from the inclusive near side deltaEta deltaPhi distributions
  
  // If the inclusive histograms are not found from the file, return correction from calculation
  if(!fDEtaDPhi3DNearLoaded) return GetAcceptanceCorrection3DNearSideCalculation(deltaEta,deltaPhi);
  
  // If there are less than fMinCountsPerBinInclusive entries per bin use calculation instead of histogram
  // This can happen only if after rebin of 16 the entries per bin are still low
  // If z-vertex bin not specified in the argument list, return sum over z-vertex bins
  TH2D *acceptanceHistogram = fDEtaDPhi3DNearAcceptance[firstBin][centralityBin][zVertexBin][triggerBin][1-firstBin];
  const int nBinsEta = acceptanceHistogram->GetNbinsX();
  const int nBinsPhi = acceptanceHistogram->GetNbinsY();
  const int nBins = nBinsEta*nBinsPhi;
  
  if(acceptanceHistogram->GetEntries() < nBins*fMinCountsPerBinInclusive) return GetAcceptanceCorrection3DNearSideCalculation(deltaEta,deltaPhi);
  
  // First find the length of the deltaEta line in the near side
  double nearSideLength = acceptanceHistogram->GetBinContent(acceptanceHistogram->FindBin(deltaEta,deltaPhi));
  
  // Then calculate the length of the deltaEta line outside of the acceptance
  double outsideAcceptance = sqrt(pow(deltaEta,2)+pow(deltaEta,2));
  
  // Calculate the correction to histogram based on these results
  if(nearSideLength + outsideAcceptance < 1e-6) return 0;
  double denominator = nearSideLength / (nearSideLength + outsideAcceptance);
  
  // Return the correction
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
}

/*
 * Return a calculated acceptance correction for 3D near side. The value from the histogram
 * needs to be corrected for the away side effect before returning the correction
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int triggerBin = trigger pT bin
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection3DNearSideInclusive(double deltaEta, double deltaPhi, int centralityBin, int triggerBin){
  return GetAcceptanceCorrection3DNearSideInclusiveBin(deltaEta, deltaPhi, centralityBin, 0, triggerBin, 0);
}

/*
 * Return a calculated acceptance correction for 3D near side. The value from the histogram
 * needs to be corrected for the away side effect before returning the correction
 *
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = centrality bin for the particle pair
 *  int zVertexBin = z-vertex bin
 *  int triggerBin = trigger pT bin
 *
 *  return = Acceptance correction based on input deltaEta and deltaPhi
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection3DNearSideInclusive(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin){
  return GetAcceptanceCorrection3DNearSideInclusiveBin(deltaEta, deltaPhi, centralityBin, zVertexBin, triggerBin, 1);
}

/*
 * Get the acceptance correction for traditional near side definition either from triangle
 * or from 2D deltaEta deltaPhi mixed event distribution based on the sampling method flag
 *
 *  int samplingMethod = 0 for triangle, 1 for mixed event
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = bin index for centrality
 *  int triggerBin = bin index for trigger particle transverse momentum
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionTraditional(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin){
  
  if(samplingMethod == 0){
    return GetAcceptanceCorrectionTriangle(deltaEta);
  } else {
    return GetAcceptanceCorrectionTraditionalInclusive(deltaEta,deltaPhi,centralityBin,triggerBin);
  }
  
}

/*
 * Get the acceptance correction for 3D near side definition either from triangle
 * or from 2D deltaEta deltaPhi mixed event distribution based on the sampling method flag
 *
 *  int samplingMethod = 0 for triangle, 1 for mixed event
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = bin index for centrality
 *  int triggerBin = bin index for trigger particle transverse momentum
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection3DNearSide(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin){
  
  if(samplingMethod == 0){
    return GetAcceptanceCorrection3DNearSideCalculation(deltaEta,deltaPhi);
  } else {
    return GetAcceptanceCorrection3DNearSideInclusive(deltaEta,deltaPhi,centralityBin,triggerBin);
  }
  
}

/*
 * Get the acceptance correction for traditional or 3D near side either from triangle
 * or from 2D deltaEta deltaPhi mixed event distribution based on the sampling method
 * near side definition flags
 *
 *  int nearSideDefinition = 0 or 1 for 3D near side, anything else for traditional
 *  int samplingMethod = 0 for triangle, 1 for mixed event
 *  double deltaEta = deltaEta for the particle pair
 *  double deltaPhi = deltaPhi for the particle pair
 *  int centralityBin = bin index for centrality
 *  int triggerBin = bin index for trigger particle transverse momentum
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrection(int nearSideDefinition, int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin){
  
  if(nearSideDefinition == 0 || nearSideDefinition == 1){
    return GetAcceptanceCorrection3DNearSide(samplingMethod,deltaEta,deltaPhi,centralityBin,triggerBin);
  } else {
    return GetAcceptanceCorrectionTraditional(samplingMethod,deltaEta,deltaPhi,centralityBin,triggerBin);
  }
  
}

/*
 * Get the rebinning factor for histograms. The idea is that each bin in the histogram
 * has to have at least fMinCountsPerBinInclusive counts. Maximum value for rebin is 16.
 *
 *  double counts = number of entries in the histogram
 *  int nBins = number of bins in the histogram
 *  int dimension = dimension of the histogram
 */
int AliJAcceptanceCorrection::GetRebin(double counts, int nBins, int dimension){
  
  // Do rebinning requiring at least fMinCountsPerBinInclusive counts in every bin in the histogram.  Maximum rebin is 16
  int rebin = 1;
  if(counts<nBins*fMinCountsPerBinInclusive) rebin=2;
  if(counts<(nBins/pow(2.0,dimension))*fMinCountsPerBinInclusive) rebin=4;
  if(counts<(nBins/pow(4.0,dimension))*fMinCountsPerBinInclusive) rebin=8;
  if(counts<(nBins/pow(8.0,dimension))*fMinCountsPerBinInclusive) rebin=10;
  if(counts<(nBins/pow(10.0,dimension))*fMinCountsPerBinInclusive) rebin=16;
  return rebin;
  
}

/*
 * Method for rebinning and normalizing two dimensional acceptance histograms
 * 
 *  TH2 *histogram = histogram to be rebinned and normalized
 *  double peakValue = maximum value of histogram after normalization
 */
void AliJAcceptanceCorrection::RebinAndNormalize(TH2 *histogram, double peakValue){
  
  // First, find out the number of bins in the histogram
  const int nBinsEta = histogram->GetNbinsX();
  const int nBinsPhi = histogram->GetNbinsY();
  const int nBins = nBinsEta*nBinsPhi;
  
  // Find out the number of entries in the histogram. Before normalization, this can be obtained from the integral.
  const double counts  = histogram->Integral();
  
  // Set the correct number of entries to histogram. This way the information is not lost when the histogram is normalized.
  histogram->SetEntries(counts);
  
  // Do rebinning requiring at least fMinCountsPerBinInclusive counts in every bin in the histogram.  Maximum rebin is 16
  int rebin = GetRebin(counts,nBins,2);
  histogram->Rebin2D(rebin,rebin);
  
  // Normalize the rebinned histogram. GetMaximum is probably good enough for finding maximum value.
  const double maxValue = histogram->GetMaximum();
  if(maxValue > 0) histogram->Scale(peakValue/maxValue);
}