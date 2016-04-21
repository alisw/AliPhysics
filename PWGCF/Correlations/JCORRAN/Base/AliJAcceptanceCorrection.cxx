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
//
// Note (19.4.2016): 3D near side correction and 2D correction under development
//                   Only traditional 1D correction functional at the moment

#include "AliJAcceptanceCorrection.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGrid.h>


//ClassImp(AliJAcceptanceCorrection)  // I do not know if this is needed

/*
 * Default constructor.
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection() :
    fCard(),
    fDEtaNearAcceptance(),
    fDEta3DNearAcceptance(),
    fDEtaDPhiNearAcceptance(),
    fDEtaDPhi3DNearAcceptance(),
    fInclusiveHistogramsNotFound(false)
{
  // default constructor
}

/*
 * Constructor that should be used for correct functioning of the class.
 * Requires definition of the JCard.
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection(AliJCard *inputCard) :
    fCard(inputCard),
    fDEtaNearAcceptance(),
    fDEta3DNearAcceptance(),
    fDEtaDPhiNearAcceptance(),
    fDEtaDPhi3DNearAcceptance(),
    fInclusiveHistogramsNotFound(false)
{
  // Constructor with JCard
}

/*
 * Copy constructor
 */
AliJAcceptanceCorrection::AliJAcceptanceCorrection(const AliJAcceptanceCorrection& a):
    fCard(a.fCard),
    fDEtaNearAcceptance(a.fDEtaNearAcceptance),
    fDEta3DNearAcceptance(a.fDEta3DNearAcceptance),
    fDEtaDPhiNearAcceptance(a.fDEtaDPhiNearAcceptance),
    fDEtaDPhi3DNearAcceptance(a.fDEtaDPhi3DNearAcceptance),
    fInclusiveHistogramsNotFound(a.fInclusiveHistogramsNotFound)
{
  //copy constructor
}


/*
 * Equal sign operator
 */
AliJAcceptanceCorrection&  AliJAcceptanceCorrection::operator=(const AliJAcceptanceCorrection& a){
  //operator = 
  if(this != &a){
    fCard = a.fCard;
    fDEtaNearAcceptance = a.fDEtaNearAcceptance;
    fDEta3DNearAcceptance = a.fDEta3DNearAcceptance;
    fDEtaDPhiNearAcceptance = a.fDEtaDPhiNearAcceptance;
    fDEtaDPhi3DNearAcceptance = a.fDEtaDPhi3DNearAcceptance;
    fInclusiveHistogramsNotFound = a.fInclusiveHistogramsNotFound;
  }
  return *this;
}

/*
 *  Method for reading acceptance correction histograms from the file.
 *  This method tries to read histograms for all possible corrections.
 *  If some of the histograms are missing, histogram manager gives an
 *  empty histogram at their place.
 *
 *  const char* inclusFileName = file name for the file, in which the histogram are located
 */
void AliJAcceptanceCorrection::ReadInclusiveHistos(const char *inclusFileName){
  // read inclusive histos
  
  TPMERegexp sep("::");
  int ncol = sep.Split( inclusFileName );
  TString filename = sep[0];
  
  if (TString(inclusFileName).BeginsWith("alien:"))  TGrid::Connect("alien:");
  TFile *inclusFile = TFile::Open(filename);
  TDirectory * dir =  (TDirectory*) inclusFile;
  if( ncol > 1 ) dir = (TDirectory*)( inclusFile->Get(sep[1]));
  if( !dir ) {
    std::cout << " ReadInclusiveHistos wrong file name or dirname !!!!" << std::endl;
    return;
  }

  AliJHistManager *histogramReader = new AliJHistManager("hst",sep[1]);
  
  // Check if the new acceptance histogram exists in the inclusive file.
  // This check is done for backwards compatibility. If the histogram does not exist,
  // mark that with a boolean flag so that we can use the information in the future.
  if(histogramReader->HistogramExists("hDEtaNearMixAcceptance")){
    fDEtaNearAcceptance = histogramReader->GetTH1D("hDEtaNearMixAcceptance");
    NormalizeAcceptanceTraditional(fDEtaNearAcceptance, kAssocType);
  } else {
    fInclusiveHistogramsNotFound = true;
  }
  
  if(histogramReader->HistogramExists("hDEta3DNearMixAcceptance")){
    fDEta3DNearAcceptance = histogramReader->GetTH1D("hDEta3DNearMixAcceptance");
    NormalizeAcceptanceTraditional(fDEta3DNearAcceptance, kXeType);
  }
  
  if(histogramReader->HistogramExists("hDPhiDEtaPta")){
    fDEtaDPhiNearAcceptance = histogramReader->GetTH2D("hDPhiDEtaPta");
    //NormalizeAcceptanceHistos(fhDEta3DNearMixFromFile, kXeType);
  }
  
  if(histogramReader->HistogramExists("hDPhiDEtaXlong")){
    fDEtaDPhi3DNearAcceptance = histogramReader->GetTH2D("hDPhiDEtaXlong");
    //NormalizeAcceptanceHistos(fhDEta3DNearMixFromFile, kXeType);
  }
  
}

/*
 * Traditional method for normalizing the acceptance histogram
 *
 * In this method the histogram is firts rebinned and then normalized to interval [0,1]
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
  
  // Loop over the input histograms and find the correct normalization
  for (int iCent = 0; iCent < numCent; iCent++) {
    for (int iPtt = 0; iPtt < numPtt; iPtt++){
      for (int iAssoc = 0; iAssoc < numAssoc; iAssoc++){
        
        // Rebin and normalize to the interval [0,1]
        counts  = acceptanceHisto[iCent][iPtt][iAssoc]->Integral();
        rebin = 4;
        if(counts<5000) rebin=8;
        if(counts<3000) rebin=10;
        if(counts<1000) rebin=16;
        acceptanceHisto[iCent][iPtt][iAssoc]->Rebin(rebin);
        maxValue = acceptanceHisto[iCent][iPtt][iAssoc]->GetMaximum();
        if(maxValue > 0) acceptanceHisto[iCent][iPtt][iAssoc]->Scale(1.0/maxValue);
        
      }
    }
  }
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
 * Calculate the acceptance correction from input histogram given in ReadInclusiveHistos
 *
 *  double deltaEta = deltaEta for the particle pair
 *  int centralityBin = bin index for centrality
 *  int triggerBin = bin index for trigger particle transverse momentum
 *  int assocBin = bin index for associated particle binning
 *  int assocType = index for associated particle type (0 = kLong, 1 = xLong, 2 = pTa)
 */
double AliJAcceptanceCorrection::GetAcceptanceCorrectionInclusive(double deltaEta, int centralityBin, int triggerBin, int assocBin, int assocType){
  // Inclusive acceptance correction
  
  // If the inclusive histograms are not found from the file, return correction from triangle
  if(fInclusiveHistogramsNotFound) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // If the given bin is negative, return correction from triangle
  if(assocBin < 0) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // Define the acceptance histogram
  TH1D *acceptanceHistogram;
  
  // Choose different acceptance histogram for xlong bins and other bins
  if(assocType == 1){
    acceptanceHistogram = fDEta3DNearAcceptance[centralityBin][triggerBin][assocBin];
  } else {
    acceptanceHistogram = fDEtaNearAcceptance[centralityBin][triggerBin][assocBin];
  }
  
  // If there is less than 1000 entries, just use triangle instead of inclusive
  if(acceptanceHistogram->GetEntries()<1000) return GetAcceptanceCorrectionTriangle(deltaEta);
  
  // Use the value in the bin corresponding to deltaEta as acceptance correction
  int bin =  acceptanceHistogram->FindBin(deltaEta);
  double denominator  =  acceptanceHistogram->GetBinContent(bin);
  if(denominator > 1e-6)
    return 1.0/denominator;
  else
    return 0;
  
}