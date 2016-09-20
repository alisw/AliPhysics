#include "THnSparseDefinitions.h"

void runPID(TString fileName, Double_t deta, Double_t pLow, Double_t pHigh, Bool_t isMCdataSet, Int_t fitMethod,
            Int_t muonFractionHandlingParameter,
            Bool_t useIdentifiedGeneratedSpectra, Bool_t plotIdentifiedSpectra, Int_t mode, Int_t chargeMode,
            Double_t lowerCentrality, Double_t upperCentrality, 
            Double_t lowerJetPt, Double_t upperJetPt,
            Int_t rebin, Int_t rebinDeltaPrime,
            TString listName, Bool_t useLogLikelihood, Bool_t useWeightsForLogLikelihood,
            Int_t regularisation,
            Double_t regularisationFactor,
            Bool_t applyTOFpatching,
            TString filePathNameFileWithInititalFractions,
            Int_t binTypePt = kPtBinTypeJets)
{
  // Load class for fit
  gROOT->LoadMacro("histFitting/AliTPCPIDmathFit.cxx+g");
  gROOT->LoadMacro("PID.C+");
  PID(fileName, deta, pLow, pHigh, isMCdataSet, fitMethod, muonFractionHandlingParameter, useIdentifiedGeneratedSpectra,
      plotIdentifiedSpectra, mode, chargeMode, lowerCentrality, upperCentrality, lowerJetPt, upperJetPt, rebin, rebinDeltaPrime,
      listName, useLogLikelihood, useWeightsForLogLikelihood, regularisation, regularisationFactor, applyTOFpatching,
      filePathNameFileWithInititalFractions, 0x0, binTypePt);
}