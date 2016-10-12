#include "THnSparseDefinitions.h"

void runPIDiterative(TString fileName, Double_t deta, Double_t pLow, Double_t pHigh, Bool_t isMCdataSet, Int_t fitMethod, Int_t 
                     muonFractionHandlingParameter, Bool_t useIdentifiedGeneratedSpectra, Bool_t plotIdentifiedSpectra, Int_t mode,
                     Int_t chargeMode, Double_t lowerCentrality, Double_t upperCentrality, 
                     Double_t lowerJetPt, Double_t upperJetPt,
                     Int_t rebin, Int_t rebinDeltaPrime,
                     TString listName, Bool_t useLogLikelihood,
                     Bool_t useWeightsForLogLikelihood, Int_t regularisation, Double_t regularisationFactor,
                     Bool_t applyTOFpatching, Int_t binTypePt = kPtBinTypeJets, Double_t yieldThresholdForFitting = 100)
{
  // Load class for fit
  gROOT->LoadMacro("histFitting/AliTPCPIDmathFit.cxx+g");
  gROOT->LoadMacro("PID.C+");
  
  TString filePathNameResultsWithoutReg = "";
  
  // First estimate w/o regularisation
  Int_t returnValue = PID(fileName, deta, pLow, pHigh, isMCdataSet, fitMethod, muonFractionHandlingParameter,
                          useIdentifiedGeneratedSpectra, plotIdentifiedSpectra, mode, chargeMode, 
                          lowerCentrality, upperCentrality, lowerJetPt, upperJetPt, rebin, rebinDeltaPrime, listName, useLogLikelihood, 
                          useWeightsForLogLikelihood, 0, regularisationFactor, applyTOFpatching, "", &filePathNameResultsWithoutReg, binTypePt,
                          yieldThresholdForFitting);
  
  // In case of success, use the result to run w/ regularisation.
  // In case of error, there is little chance that the fitting with regularisation will succeed, but just try without initial values
  if (returnValue == 0) {
    Int_t returnValue2 = PID(fileName, deta, pLow, pHigh, isMCdataSet, fitMethod, muonFractionHandlingParameter, useIdentifiedGeneratedSpectra,
                             plotIdentifiedSpectra, mode, chargeMode, lowerCentrality, upperCentrality, lowerJetPt, upperJetPt, rebin, 
                             rebinDeltaPrime, listName, useLogLikelihood, useWeightsForLogLikelihood, regularisation, regularisationFactor, 
                             applyTOFpatching, filePathNameResultsWithoutReg, 0x0, binTypePt, yieldThresholdForFitting);
    
    // Sometimes the fit crashes if initial values are used. In that case, try again without initial values (it sometimes works!)
    PID(fileName, deta, pLow, pHigh, isMCdataSet, fitMethod, muonFractionHandlingParameter, useIdentifiedGeneratedSpectra,
        plotIdentifiedSpectra, mode, chargeMode, lowerCentrality, upperCentrality, lowerJetPt, upperJetPt, rebin, rebinDeltaPrime,
        listName, useLogLikelihood, useWeightsForLogLikelihood, regularisation, regularisationFactor, applyTOFpatching, "", 0x0, binTypePt,
        yieldThresholdForFitting);
  }
  else
    PID(fileName, deta, pLow, pHigh, isMCdataSet, fitMethod, muonFractionHandlingParameter, useIdentifiedGeneratedSpectra,
        plotIdentifiedSpectra, mode, chargeMode, lowerCentrality, upperCentrality, lowerJetPt, upperJetPt, rebin, rebinDeltaPrime,
        listName, useLogLikelihood, useWeightsForLogLikelihood, regularisation, regularisationFactor, applyTOFpatching, "", 0x0, binTypePt,
        yieldThresholdForFitting);
}