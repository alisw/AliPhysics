#include "TChain.h"
#include "TFile.h"
#include "TF1.h"
#include "TAxis.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"

#include "AliMCParticle.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliInputEventHandler.h"

#include "AliVVertex.h"
#include "AliAnalysisFilter.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

#include "AliAnalysisTaskPID.h"

/*
This task collects PID output from different detectors.
Only tracks fulfilling some standard quality cuts are taken into account.
At the moment, only data from TPC and TOF is collected. But in future,
data from e.g. HMPID is also foreseen.

Contact: bhess@cern.ch
*/

ClassImp(AliAnalysisTaskPID)

//________________________________________________________________________
AliAnalysisTaskPID::AliAnalysisTaskPID()
  : AliAnalysisTaskPIDV0base()
  , fPIDcombined(new AliPIDCombined())
  , fInputFromOtherTask(kFALSE)
  , fStoreCentralityPercentile(kFALSE)
  , fStoreAdditionalJetInformation(kFALSE)
  , fTakeIntoAccountMuons(kFALSE)
  , fUseITS(kFALSE)
  , fUseTOF(kFALSE)
  , fUsePriors(kFALSE)
  , fTPCDefaultPriors(kFALSE)
  , fUseMCidForGeneration(kTRUE)
  , fUseConvolutedGaus(kFALSE) 
  , fkConvolutedGausNPar(3)
  , fAccuracyNonGaussianTail(1e-8)
  , fkDeltaPrimeLowLimit(0.02)
  , fkDeltaPrimeUpLimit(40.0)
  , fConvolutedGausDeltaPrime(0x0)
  , fEtaAbsCutLow(0.0)
  , fEtaAbsCutUp(0.9)
  , fDoAnySystematicStudiesOnTheExpectedSignal(kFALSE)
  , fSystematicScalingSplines(1.0)
  , fSystematicScalingEtaCorrectionMomentumThr(0.35)
  , fSystematicScalingEtaCorrectionLowMomenta(1.0)
  , fSystematicScalingEtaCorrectionHighMomenta(1.0)
  , fSystematicScalingEtaSigmaPara(1.0)
  , fSystematicScalingMultCorrection(1.0)
  , fCentralityEstimator("V0M")
  , fhPIDdataAll(0x0)
  , fhGenEl(0x0)
  , fhGenKa(0x0)
  , fhGenPi(0x0)
  , fhGenMu(0x0)
  , fhGenPr(0x0)
  , fGenRespElDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  /*
  , fGenRespElDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPr(new Double_t[fgkMaxNumGenEntries])
  */
  , fhEventsProcessed(0x0)
  , fhSkippedTracksForSignalGeneration(0x0)
  , fhMCgeneratedYieldsPrimaries(0x0)
  , fh2FFJetPtRec(0x0)
  , fh2FFJetPtGen(0x0)
  , fh1Xsec(0x0)
  , fh1Trials(0x0)
  , fContainerEff(0x0)
  , fOutputContainer(0x0)
{
  // default Constructor
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskPID", AliLog::kInfo); 
  
  fConvolutedGausDeltaPrime = new TF1("convolutedGausDeltaPrime", this, &AliAnalysisTaskPID::ConvolutedGaus,
                                      fkDeltaPrimeLowLimit, fkDeltaPrimeUpLimit,
                                      fkConvolutedGausNPar, "AliAnalysisTaskPID", "ConvolutedGaus");
  
  // Initialisation of translation parameters is time consuming.
  // Therefore, default values will only be initialised if they are really needed.
  // Otherwise, it is left to the user to set the parameter properly.
  fConvolutedGaussTransitionPars[0] = -999;
  fConvolutedGaussTransitionPars[1] = -999;
  fConvolutedGaussTransitionPars[2] = -999;
  
  // Fraction histos
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    fFractionHists[i] = 0x0;
    fFractionSysErrorHists[i] = 0x0;
  }
}

//________________________________________________________________________
AliAnalysisTaskPID::AliAnalysisTaskPID(const char *name)
  : AliAnalysisTaskPIDV0base(name)
  , fPIDcombined(new AliPIDCombined())
  , fInputFromOtherTask(kFALSE)
  , fStoreCentralityPercentile(kFALSE)
  , fStoreAdditionalJetInformation(kFALSE)
  , fTakeIntoAccountMuons(kFALSE)
  , fUseITS(kFALSE)
  , fUseTOF(kFALSE)
  , fUsePriors(kFALSE)
  , fTPCDefaultPriors(kFALSE)
  , fUseMCidForGeneration(kTRUE)
  , fUseConvolutedGaus(kFALSE) 
  , fkConvolutedGausNPar(3)
  , fAccuracyNonGaussianTail(1e-8)
  , fkDeltaPrimeLowLimit(0.02)
  , fkDeltaPrimeUpLimit(40.0)
  , fConvolutedGausDeltaPrime(0x0)
  , fEtaAbsCutLow(0.0)
  , fEtaAbsCutUp(0.9)
  , fDoAnySystematicStudiesOnTheExpectedSignal(kFALSE)
  , fSystematicScalingSplines(1.0)
  , fSystematicScalingEtaCorrectionMomentumThr(0.35)
  , fSystematicScalingEtaCorrectionLowMomenta(1.0)
  , fSystematicScalingEtaCorrectionHighMomenta(1.0)
  , fSystematicScalingEtaSigmaPara(1.0)
  , fSystematicScalingMultCorrection(1.0)
  , fCentralityEstimator("V0M")
  , fhPIDdataAll(0x0)
  , fhGenEl(0x0)
  , fhGenKa(0x0)
  , fhGenPi(0x0)
  , fhGenMu(0x0)
  , fhGenPr(0x0)
  , fGenRespElDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimeEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimeKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimePi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPrimePr(new Double_t[fgkMaxNumGenEntries])
  /*
  , fGenRespElDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespElDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespKaDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPiDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespMuDeltaPr(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaEl(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaKa(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPi(new Double_t[fgkMaxNumGenEntries])
  , fGenRespPrDeltaPr(new Double_t[fgkMaxNumGenEntries])
  */
  , fhEventsProcessed(0x0)
  , fhSkippedTracksForSignalGeneration(0x0)
  , fhMCgeneratedYieldsPrimaries(0x0)
  , fh2FFJetPtRec(0x0)
  , fh2FFJetPtGen(0x0)
  , fh1Xsec(0x0)
  , fh1Trials(0x0)
  , fContainerEff(0x0)
  , fOutputContainer(0x0)
{
  // Constructor
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskPID", AliLog::kInfo);
  
  fConvolutedGausDeltaPrime = new TF1("convolutedGausDeltaPrime", this, &AliAnalysisTaskPID::ConvolutedGaus,
                                      fkDeltaPrimeLowLimit, fkDeltaPrimeUpLimit,
                                      fkConvolutedGausNPar, "AliAnalysisTaskPID", "ConvolutedGaus");
  
  // Initialisation of translation parameters is time consuming.
  // Therefore, default values will only be initialised if they are really needed.
  // Otherwise, it is left to the user to set the parameter properly.
  fConvolutedGaussTransitionPars[0] = -999;
  fConvolutedGaussTransitionPars[1] = -999;
  fConvolutedGaussTransitionPars[2] = -999;
  
  // Fraction histos
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    fFractionHists[i] = 0x0;
    fFractionSysErrorHists[i] = 0x0;
  }
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TObjArray::Class());
  
  DefineOutput(2, AliCFContainer::Class());
}


//________________________________________________________________________
AliAnalysisTaskPID::~AliAnalysisTaskPID()
{
  // dtor
  
  CleanupParticleFractionHistos();
  
  delete fOutputContainer;
  fOutputContainer = 0;

  delete fConvolutedGausDeltaPrime;
  fConvolutedGausDeltaPrime = 0;
  
  delete [] fGenRespElDeltaPrimeEl;
  delete [] fGenRespElDeltaPrimeKa;
  delete [] fGenRespElDeltaPrimePi;
  delete [] fGenRespElDeltaPrimePr;
  
  fGenRespElDeltaPrimeEl = 0x0;
  fGenRespElDeltaPrimeKa = 0x0;
  fGenRespElDeltaPrimePi = 0x0;
  fGenRespElDeltaPrimePr = 0x0;
  
  delete [] fGenRespKaDeltaPrimeEl;
  delete [] fGenRespKaDeltaPrimeKa;
  delete [] fGenRespKaDeltaPrimePi;
  delete [] fGenRespKaDeltaPrimePr;
  
  fGenRespKaDeltaPrimeEl = 0x0;
  fGenRespKaDeltaPrimeKa = 0x0;
  fGenRespKaDeltaPrimePi = 0x0;
  fGenRespKaDeltaPrimePr = 0x0;
  
  delete [] fGenRespPiDeltaPrimeEl;
  delete [] fGenRespPiDeltaPrimeKa;
  delete [] fGenRespPiDeltaPrimePi;
  delete [] fGenRespPiDeltaPrimePr;
  
  fGenRespPiDeltaPrimeEl = 0x0;
  fGenRespPiDeltaPrimeKa = 0x0;
  fGenRespPiDeltaPrimePi = 0x0;
  fGenRespPiDeltaPrimePr = 0x0;
  
  delete [] fGenRespMuDeltaPrimeEl;
  delete [] fGenRespMuDeltaPrimeKa;
  delete [] fGenRespMuDeltaPrimePi;
  delete [] fGenRespMuDeltaPrimePr;
  
  fGenRespMuDeltaPrimeEl = 0x0;
  fGenRespMuDeltaPrimeKa = 0x0;
  fGenRespMuDeltaPrimePi = 0x0;
  fGenRespMuDeltaPrimePr = 0x0;
  
  delete [] fGenRespPrDeltaPrimeEl;
  delete [] fGenRespPrDeltaPrimeKa;
  delete [] fGenRespPrDeltaPrimePi;
  delete [] fGenRespPrDeltaPrimePr;
  
  fGenRespPrDeltaPrimeEl = 0x0;
  fGenRespPrDeltaPrimeKa = 0x0;
  fGenRespPrDeltaPrimePi = 0x0;
  fGenRespPrDeltaPrimePr = 0x0;
  
  /*OLD with deltaSpecies 
  delete [] fGenRespElDeltaEl;
  delete [] fGenRespElDeltaKa;
  delete [] fGenRespElDeltaPi;
  delete [] fGenRespElDeltaPr;
  
  fGenRespElDeltaEl = 0x0;
  fGenRespElDeltaKa = 0x0;
  fGenRespElDeltaPi = 0x0;
  fGenRespElDeltaPr = 0x0;
  
  delete [] fGenRespKaDeltaEl;
  delete [] fGenRespKaDeltaKa;
  delete [] fGenRespKaDeltaPi;
  delete [] fGenRespKaDeltaPr;
  
  fGenRespKaDeltaEl = 0x0;
  fGenRespKaDeltaKa = 0x0;
  fGenRespKaDeltaPi = 0x0;
  fGenRespKaDeltaPr = 0x0;
  
  delete [] fGenRespPiDeltaEl;
  delete [] fGenRespPiDeltaKa;
  delete [] fGenRespPiDeltaPi;
  delete [] fGenRespPiDeltaPr;
  
  fGenRespPiDeltaEl = 0x0;
  fGenRespPiDeltaKa = 0x0;
  fGenRespPiDeltaPi = 0x0;
  fGenRespPiDeltaPr = 0x0;
  
  delete [] fGenRespMuDeltaEl;
  delete [] fGenRespMuDeltaKa;
  delete [] fGenRespMuDeltaPi;
  delete [] fGenRespMuDeltaPr;
  
  fGenRespMuDeltaEl = 0x0;
  fGenRespMuDeltaKa = 0x0;
  fGenRespMuDeltaPi = 0x0;
  fGenRespMuDeltaPr = 0x0;
  
  delete [] fGenRespPrDeltaEl;
  delete [] fGenRespPrDeltaKa;
  delete [] fGenRespPrDeltaPi;
  delete [] fGenRespPrDeltaPr;
  
  fGenRespPrDeltaEl = 0x0;
  fGenRespPrDeltaKa = 0x0;
  fGenRespPrDeltaPi = 0x0;
  fGenRespPrDeltaPr = 0x0;
  */
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpPIDcombined()
{
  // Initialise the PIDcombined object
  
  if (!fPIDcombined) {
    AliFatal("No PIDcombined object!\n");
    return;
  }
  
  fPIDcombined->SetSelectedSpecies(AliPID::kSPECIESC);
  fPIDcombined->SetEnablePriors(fUsePriors);
 
  if (fTPCDefaultPriors)
    fPIDcombined->SetDefaultTPCPriors();
  
  //TODO use individual priors...
  
  // Change detector mask (TPC,TOF,ITS)
  Int_t detectorMask = AliPIDResponse::kDetTPC;
  
  // Reject mismatch mask - mismatch only relevant for TOF at the moment - other detectors do not use it
  Int_t rejectMismatchMask = AliPIDResponse::kDetTPC;
  
  
  if (fUseITS) {
    detectorMask = detectorMask | AliPIDResponse::kDetITS;
    rejectMismatchMask = rejectMismatchMask | AliPIDResponse::kDetITS;
  }
  if (fUseTOF) {
    detectorMask = detectorMask | AliPIDResponse::kDetTOF;
    rejectMismatchMask = rejectMismatchMask | AliPIDResponse::kDetTOF;
  }
  
  fPIDcombined->SetDetectorMask(detectorMask);
  fPIDcombined->SetRejectMismatchMask(rejectMismatchMask);
}


//________________________________________________________________________
void AliAnalysisTaskPID::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  SetUpPIDcombined();

  // Input handler
  AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  
  if (!inputHandler)
    AliFatal("Input handler needed");
  else {
    // PID response object
    fPIDResponse = inputHandler->GetPIDResponse();
    if (!fPIDResponse)
      AliFatal("PIDResponse object was not created");
  }
  
  OpenFile(1);
    
  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetName(GetName()) ;
  fOutputContainer->SetOwner(kTRUE);
  
  const Int_t nPtBins = 68;
  Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
           0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
           1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
           2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
           4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
           26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
  
  const Int_t nCentBins = 12;
  //-1 for pp; 90-100 has huge electro-magnetic impurities
  Double_t binsCent[nCentBins+1] = {-1, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 }; 

  const Int_t nJetPtBins = 11;
  Double_t binsJetPt[nJetPtBins+1] = {0, 2, 5, 10, 15, 20, 30, 40, 60, 80, 120, 200};
  
  const Int_t nChargeBins = 2;
  const Double_t binsCharge[nChargeBins+1] = { -1.0 - 1e-4, 0.0, 1.0 + 1e-4 };
  
  const Int_t nBinsJets = kDataNumAxes;
  const Int_t nBinsNoJets = nBinsJets - fgkNumJetAxes;
  
  const Int_t nBins = fStoreAdditionalJetInformation ? nBinsJets : nBinsNoJets;
 
  // deltaPrimeSpecies binning
  const Int_t deltaPrimeNBins = 600;
  Double_t deltaPrimeBins[deltaPrimeNBins + 1];

  const Double_t fromLow = fkDeltaPrimeLowLimit;
  const Double_t toHigh = fkDeltaPrimeUpLimit;
  const Double_t factor = TMath::Power(toHigh/fromLow, 1./deltaPrimeNBins);

  // Log binning for whole deltaPrime range
  deltaPrimeBins[0] = fromLow;
  for (Int_t i = 0 + 1; i <= deltaPrimeNBins; i++) {
    deltaPrimeBins[i] = factor * deltaPrimeBins[i - 1];
  }
  
  const Int_t nMCPIDbins = 5;
  const Double_t mcPIDmin = 0.;
  const Double_t mcPIDmax = 5.;
  
  const Int_t nSelSpeciesBins = 4;
  const Double_t selSpeciesMin = 0.;
  const Double_t selSpeciesMax = 4.;
  
  const Int_t nZBins = 20;
  const Double_t zMin = 0.;
  const Double_t zMax = 1.;
  
  const Int_t nXiBins = 70;
  const Double_t xiMin = 0.;
  const Double_t xiMax = 7.;
  
  // MC PID, SelectSpecies, pT, deltaPrimeSpecies, centrality percentile, jet pT, z = track_pT/jet_pT, xi = log(1/z)
  Int_t binsNoJets[nBinsNoJets] =    { nMCPIDbins, nSelSpeciesBins,         nPtBins,                 deltaPrimeNBins,  
                                       nCentBins, nChargeBins};
  Int_t binsJets[nBinsJets]     =    { nMCPIDbins, nSelSpeciesBins,         nPtBins,                 deltaPrimeNBins,           
                                       nCentBins,                      nJetPtBins, nZBins, nXiBins, nChargeBins };
  Int_t *bins = fStoreAdditionalJetInformation ? &binsJets[0] : &binsNoJets[0];
  
  Double_t xminNoJets[nBinsNoJets] = {   mcPIDmin,   selSpeciesMin,       binsPt[0],               deltaPrimeBins[0],                       
                                       binsCent[0], binsCharge[0]};
  Double_t xminJets[nBinsJets] =     {   mcPIDmin,   selSpeciesMin,       binsPt[0],               deltaPrimeBins[0],                       
                                       binsCent[0],                  binsJetPt[0],   zMin,   xiMin, binsCharge[0] };
  Double_t *xmin = fStoreAdditionalJetInformation? &xminJets[0] : &xminNoJets[0];

  Double_t xmaxNoJets[nBinsNoJets] = {   mcPIDmax,   selSpeciesMax, binsPt[nPtBins], deltaPrimeBins[deltaPrimeNBins], 
                                       binsCent[nCentBins], binsCharge[nChargeBins] };
  Double_t xmaxJets[nBinsJets] =     {   mcPIDmax,   selSpeciesMax, binsPt[nPtBins], deltaPrimeBins[deltaPrimeNBins], 
                                       binsCent[nCentBins], binsJetPt[nJetPtBins],   zMax,   xiMax, binsCharge[nChargeBins] };
  Double_t *xmax = fStoreAdditionalJetInformation? &xmaxJets[0] : &xmaxNoJets[0];
  
  /*OLD with TOF, p_TPC_Inner and p_vertex and deltaSpecies
  const Int_t nBins = 8;
  //TODO In case of memory trouble: Remove deltaTOFspecies and p(Vertex) (can be added later, if needed)?
  //TODO IF everything is working fine: p(TPC_inner) and p(Vertex) can be removed, since everything in the analysis is only a 
  // function of pT -> Reduces memory consumption significantly!!!
  
  // MC PID, SelectSpecies, P(TPC_inner), pT, p(Vertex), deltaSpecies, deltaPrimeSpecies, deltaTOFspecies
  const Int_t deltaPrimeNBins = 201;
  
  const Int_t deltaNBins = 801;
  const Double_t deltaLowLimit = -200.;
  const Double_t deltaUpLimit = 200.;
  
  Int_t bins[nBins] = 
    {  5, 4, nPtBins, nPtBins, nPtBins, deltaNBins, deltaPrimeNBins, 501 }; 
  Double_t xmin[nBins] = 
    {  0., 0., 0., 0., 0., deltaLowLimit, fkDeltaPrimeLowLimit, -5000.};
  Double_t xmax[nBins] = 
    {  5., 4., 50.0, 50.0, 50.0, deltaUpLimit, fkDeltaPrimeUpLimit, 5000.};
  */
  
  fConvolutedGausDeltaPrime->SetNpx(deltaPrimeNBins);

  fhPIDdataAll = new THnSparseD("hPIDdataAll","", nBins, bins, xmin, xmax);
  SetUpHist(fhPIDdataAll, binsPt, deltaPrimeBins, binsCent, binsJetPt);
  fOutputContainer->Add(fhPIDdataAll);
  
  // Generated histograms (so far, bins are the same as for primary THnSparse)
  const Int_t nGenBins = fStoreAdditionalJetInformation ? nBinsJets : nBinsNoJets;
  // MC PID, SelectSpecies, Pt, deltaPrimeSpecies, jet pT, z = track_pT/jet_pT, xi = log(1/z)
  
  Int_t *genBins = fStoreAdditionalJetInformation ? &binsJets[0] : &binsNoJets[0];
  Double_t *genXmin = fStoreAdditionalJetInformation? &xminJets[0] : &xminNoJets[0];
  Double_t *genXmax = fStoreAdditionalJetInformation? &xmaxJets[0] : &xmaxNoJets[0];

  fhGenEl = new THnSparseD("hGenEl", "", nGenBins, genBins, genXmin, genXmax);
  SetUpGenHist(fhGenEl, binsPt, deltaPrimeBins, binsCent, binsJetPt);
  fOutputContainer->Add(fhGenEl);
  
  fhGenKa = new THnSparseD("hGenKa", "", nGenBins, genBins, genXmin, genXmax);
  SetUpGenHist(fhGenKa, binsPt, deltaPrimeBins, binsCent, binsJetPt);
  fOutputContainer->Add(fhGenKa);
  
  fhGenPi = new THnSparseD("hGenPi", "", nGenBins, genBins, genXmin, genXmax);
  SetUpGenHist(fhGenPi, binsPt, deltaPrimeBins, binsCent, binsJetPt);
  fOutputContainer->Add(fhGenPi);
  
  if (fTakeIntoAccountMuons) {
    fhGenMu = new THnSparseD("hGenMu", "", nGenBins, genBins, genXmin, genXmax);
    SetUpGenHist(fhGenMu, binsPt, deltaPrimeBins, binsCent, binsJetPt);
    fOutputContainer->Add(fhGenMu);
  }
  
  fhGenPr = new THnSparseD("hGenPr", "", nGenBins, genBins, genXmin, genXmax);
  SetUpGenHist(fhGenPr, binsPt, deltaPrimeBins, binsCent, binsJetPt);
  fOutputContainer->Add(fhGenPr);
  
  
  fhEventsProcessed = new TH1D("fhEventsProcessed","Number of processed events;Centrality percentile", nCentBins, binsCent);
  fhEventsProcessed->Sumw2();
  fOutputContainer->Add(fhEventsProcessed);
  
  fhSkippedTracksForSignalGeneration = new TH2D("fhSkippedTracksForSignalGeneration",
                                                "Number of tracks skipped for the signal generation;P_{T}^{gen} (GeV/c);TPC signal N", 
                                                nPtBins, binsPt, 161, -0.5, 160.5);
  fhSkippedTracksForSignalGeneration->Sumw2();
  fOutputContainer->Add(fhSkippedTracksForSignalGeneration);
  
  
  // Generated yields within acceptance
  const Int_t nBinsGenYields = fStoreAdditionalJetInformation ? kGenYieldNumAxes : kGenYieldNumAxes - 3;
  Int_t genYieldsBins[kGenYieldNumAxes]    = { nMCPIDbins,         nPtBins,           nCentBins,            nJetPtBins, nZBins, nXiBins,
                                               nChargeBins };
  genYieldsBins[GetIndexOfChargeAxisGenYield()] = nChargeBins;
  Double_t genYieldsXmin[kGenYieldNumAxes] = {   mcPIDmin,       binsPt[0],         binsCent[0],          binsJetPt[0],   zMin,   xiMin,
                                               binsCharge[0] };
  genYieldsXmin[GetIndexOfChargeAxisGenYield()] = binsCharge[0];
  Double_t genYieldsXmax[kGenYieldNumAxes] = {   mcPIDmax, binsPt[nPtBins], binsCent[nCentBins], binsJetPt[nJetPtBins],   zMax,   xiMax, 
                                               binsCharge[nChargeBins] };
  genYieldsXmax[GetIndexOfChargeAxisGenYield()] = binsCharge[nChargeBins];
  
  fhMCgeneratedYieldsPrimaries = new THnSparseD("fhMCgeneratedYieldsPrimaries", 
                                                "Generated yields w/o reco and cuts inside acceptance (physical primaries)", 
                                                nBinsGenYields, genYieldsBins, genYieldsXmin, genYieldsXmax);
  SetUpGenYieldHist(fhMCgeneratedYieldsPrimaries, binsPt, binsCent, binsJetPt);
  fOutputContainer->Add(fhMCgeneratedYieldsPrimaries);
  
  // Container with several process steps (generated and reconstructed level with some variations)
  OpenFile(2);
 
  // Array for the number of bins in each dimension
  // Dimensions: MC-ID, trackPt, trackEta, trackCharge, cenrality percentile, jetPt, z, xi TODO phi???
  const Int_t nEffDims = fStoreAdditionalJetInformation ? kEffNumAxes : kEffNumAxes - 3; // Number of dimensions for the efficiency
  
  const Int_t nMCIDbins = AliPID::kSPECIES;
  Double_t binsMCID[nMCIDbins];
  
  for(Int_t i = 0; i <= nMCIDbins; i++) {
    binsMCID[i]= i; 
  }
  
  const Int_t nEtaBins = 18;
  const Double_t binsEta[nEtaBins+1] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                                         0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
  
  const Int_t nEffBins[kEffNumAxes] = { nMCIDbins, nPtBins, nEtaBins, nChargeBins, nCentBins, nJetPtBins, nZBins, nXiBins };
  
  fContainerEff = new AliCFContainer("containerEff", "Reconstruction Efficiency x Acceptance x Resolution and Secondary Correction",
                                     kNumSteps, nEffDims, nEffBins);

  // Setting the bin limits
  fContainerEff->SetBinLimits(kEffMCID, binsMCID);
  fContainerEff->SetBinLimits(kEffTrackPt, binsPt);
  fContainerEff->SetBinLimits(kEffTrackEta, binsEta);
  fContainerEff->SetBinLimits(kEffTrackCharge, binsCharge);
  fContainerEff->SetBinLimits(kEffCentrality, binsCent);
  if (fStoreAdditionalJetInformation) {
    fContainerEff->SetBinLimits(kEffJetPt, binsJetPt);
    fContainerEff->SetBinLimits(kEffZ, zMin, zMax);
    fContainerEff->SetBinLimits(kEffXi, xiMin, xiMax);
  }
  
  fContainerEff->SetVarTitle(kEffMCID,"MC ID");
  fContainerEff->SetVarTitle(kEffTrackPt,"P_{T} (GeV/c)");
  fContainerEff->SetVarTitle(kEffTrackEta,"#eta");
  fContainerEff->SetVarTitle(kEffTrackCharge,"Charge (e_{0})");
  fContainerEff->SetVarTitle(kEffCentrality, "Centrality Percentile");
  if (fStoreAdditionalJetInformation) {
    fContainerEff->SetVarTitle(kEffJetPt, "P_{T}^{jet} (GeV/c)");
    fContainerEff->SetVarTitle(kEffZ, "z = P_{T}^{track} / P_{T}^{jet}");
    fContainerEff->SetVarTitle(kEffXi, "#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  }
  
  // Define clean MC sample
  fContainerEff->SetStepTitle(kStepGenWithGenCuts, "Particle level, cuts on particle level");
  // For Acceptance x Efficiency correction of primaries
  fContainerEff->SetStepTitle(kStepRecWithGenCuts, "Detector level (rec) with cuts on particle level"); 
  // For (pT) resolution correction
  fContainerEff->SetStepTitle(kStepRecWithGenCutsMeasuredObs, 
                              "Detector level (rec) with cuts on particle level with measured observables");
  // For secondary correction
  fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObs, "Detector level, all cuts on detector level with measured observables");
  fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObsPrimaries, 
                              "Detector level, all cuts on detector level with measured observables, only MC primaries");
  fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObsStrangenessScaled, 
                              "Detector level (strangeness scaled), all cuts on detector level with measured observables");
  
  // Generated jets
  fh2FFJetPtRec = new TH2D("fh2FFJetPtRec", "Number of reconstructed jets;Centrality Percentile;P_{T}^{jet} (GeV/c)",
                           nCentBins, binsCent, nJetPtBins, binsJetPt);
  fh2FFJetPtRec->Sumw2();
  fOutputContainer->Add(fh2FFJetPtRec);
  fh2FFJetPtGen = new TH2D("fh2FFJetPtGen", "Number of generated jets;Centrality Percentile;P_{T}^{jet} (GeV/c)",
                           nCentBins, binsCent, nJetPtBins, binsJetPt);
  fh2FFJetPtGen->Sumw2();
  fOutputContainer->Add(fh2FFJetPtGen);
  
  // Pythia information
  fh1Xsec = new TProfile("fh1Xsec", "xsec from pyxsec.root", 1, 0, 1);
  fh1Xsec->Sumw2();
  fh1Xsec->GetXaxis()->SetBinLabel(1, "<#sigma>");
  fh1Trials = new TH1D("fh1Trials", "trials from pyxsec.root", 1, 0, 1);
  fh1Trials->Sumw2();
  fh1Trials->GetXaxis()->SetBinLabel(1, "#sum{ntrials}");
  
  fOutputContainer->Add(fh1Xsec);
  fOutputContainer->Add(fh1Trials);
  
  PostOutputData();
}


//________________________________________________________________________
void AliAnalysisTaskPID::UserExec(Option_t *)
{
  // Main loop
  // Called for each event
  
  // No processing of event, if input is fed in directly from another task
  if (fInputFromOtherTask)
    return;

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: fEvent not available");
    return;
  }
  
  fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  
  if (!fPIDResponse || !fPIDcombined)
    return;
  
  if (!GetVertexIsOk(fEvent))
    return;
  
  fESD = dynamic_cast<AliESDEvent*>(fEvent);
  const AliVVertex* primaryVertex = fESD ? fESD->GetPrimaryVertexTracks() : fEvent->GetPrimaryVertex(); 
  if (!primaryVertex)
    return;
  
  if(primaryVertex->GetNContributors() <= 0) 
    return;
  
  Double_t magField = fEvent->GetMagneticField();
  
  //OLD with DeltaSpecies const Bool_t usePureGausForDelta = kTRUE;
  

  Double_t centralityPercentile = -1;
  if (fStoreCentralityPercentile)
    centralityPercentile = fEvent->GetCentrality()->GetCentralityPercentile(fCentralityEstimator.Data());
  
  if (fMC) {
    for (Int_t iPart = 0; iPart < fMC->GetNumberOfTracks(); iPart++) { 
      AliMCParticle *mcPart  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(iPart));
      
      if (!mcPart)
          continue;
      
      // Define clean MC sample with corresponding particle level track cuts:
      // - MC-track must be in desired eta range
      // - MC-track must be physical primary
      // - Species must be one of those in question (everything else goes to the overflow bin of mcID)
      
      // Geometrie should be the same as on the reconstructed level -> By definition analysis within this eta interval
      if (TMath::Abs(mcPart->Eta()) < fEtaAbsCutLow || TMath::Abs(mcPart->Eta()) > fEtaAbsCutUp)  continue;
      
      Int_t mcID = PDGtoMCID(mcPart->PdgCode());
      
      // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
      Double_t chargeMC = mcPart->Charge() / 3.;
      
      if (TMath::Abs(chargeMC) < 0.01)
        continue; // Reject neutral particles (only relevant, if mcID is not used)
      
      Double_t valuesGenYield[kGenYieldNumAxes] = { mcID, mcPart->Pt(), centralityPercentile, -1, -1, -1, -1 };
      valuesGenYield[GetIndexOfChargeAxisGenYield()] = chargeMC;
      
      if (!fMC->IsPhysicalPrimary(iPart)) 
          continue;
      
      fhMCgeneratedYieldsPrimaries->Fill(valuesGenYield);
      
      
      Double_t valueEff[kEffNumAxes] = { mcID, mcPart->Pt(), mcPart->Eta(), chargeMC, centralityPercentile,
                                         -1, -1, -1 };
      
      fContainerEff->Fill(valueEff, kStepGenWithGenCuts);    
    }
  }
  
  // Track loop to fill a Train spectrum
  for (Int_t iTracks = 0; iTracks < fEvent->GetNumberOfTracks(); iTracks++) {
    AliVTrack* track = dynamic_cast<AliVTrack*>(fEvent->GetTrack(iTracks));
    if (!track) {
      Printf("ERROR: Could not retrieve track %d", iTracks);
      continue;
    }
    
    
    // Apply detector level track cuts
    //TODO NOW if (track->GetTPCsignalN() < 60)
    //TODO NOW  continue;//TODO was removed for a while
    
    
    if(fTrackFilter && !fTrackFilter->IsSelected(track))
      continue;
    
    /*
    if (fESD) {
      AliESDtrack* trackESD = fESD->GetTrack(iTracks);
      if (trackESD) {
        if (trackESD->GetLengthInActiveZone(1, 1.8, 220, magField) < 120)
           continue; //TODO NOW TEST
      }
    }*/
    
    if(fUsePhiCut) {
      if (!PhiPrimeCut(track, magField))
        continue; // reject track
    }
    
    Double_t dEdxTPC = fPIDResponse->IsTunedOnData() ? fPIDResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
    
    if (dEdxTPC <= 0)
      continue;
    
    Int_t pdg =  0; // = 0 indicates data for the moment
    AliMCParticle* mcTrack = 0x0;
    Int_t mcID = AliPID::kUnknown;
    Int_t label = 0;
    
    if (fMC) {
      label = track->GetLabel();
    
      //if (label < 0)
      //  continue;

      mcTrack = dynamic_cast<AliMCParticle*>(fMC->GetTrack(TMath::Abs(label)));
      if (!mcTrack) {
        Printf("ERROR: Could not retrieve mcTrack with label %d for track %d", label, iTracks);
        continue;
      }
      
      pdg = mcTrack->PdgCode();
      mcID = PDGtoMCID(mcTrack->PdgCode());
      
      // For efficiency: Reconstructed track has survived all cuts on the detector level (excluding acceptance)
      // and has an associated MC track which is a physical primary and was generated inside the acceptance
      if (fMC->IsPhysicalPrimary(TMath::Abs(label)) &&
          (TMath::Abs(mcTrack->Eta()) >= fEtaAbsCutLow && TMath::Abs(mcTrack->Eta()) <= fEtaAbsCutUp)) {
        
        // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        Double_t value[kEffNumAxes] = { mcID, mcTrack->Pt(), mcTrack->Eta(), mcTrack->Charge() / 3., centralityPercentile,
                                        -1, -1, -1 };
        fContainerEff->Fill(value, kStepRecWithGenCuts);    
          
        Double_t valueMeas[kEffNumAxes] = { mcID, track->Pt(), track->Eta(), track->Charge(), centralityPercentile,
                                            -1, -1, -1 };
        fContainerEff->Fill(valueMeas, kStepRecWithGenCutsMeasuredObs);    
      }
    }
   
    // Only process tracks inside the desired eta window    
    if (TMath::Abs(track->Eta()) < fEtaAbsCutLow || TMath::Abs(track->Eta()) > fEtaAbsCutUp)  continue;
   
    ProcessTrack(track, pdg, centralityPercentile, -1); // No jet information in this case -> Set jet pT to -1
    
    if (mcTrack) {
      Double_t valueRecAllCuts[kEffNumAxes] = { mcID, track->Pt(), track->Eta(), track->Charge(), centralityPercentile,
                                                -1, -1, -1 };
      fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObs);
      
      Double_t weight = IsSecondaryWithStrangeMotherMC(fMC, TMath::Abs(label)) ? GetMCStrangenessFactorCMS(fMC, mcTrack) : 1.0;
      fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsStrangenessScaled, weight);
      
      if (fMC->IsPhysicalPrimary(TMath::Abs(label)))
        fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsPrimaries);
    }
  } //track loop 
  
  IncrementEventsProcessed(centralityPercentile);

  PostOutputData();
}      

//________________________________________________________________________
void AliAnalysisTaskPID::Terminate(const Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::CheckDoAnyStematicStudiesOnTheExpectedSignal()
{
  // Check whether at least one scale factor indicates the ussage of systematic studies
  // and set internal flag accordingly.
  
  fDoAnySystematicStudiesOnTheExpectedSignal = kFALSE;
  
  if (TMath::Abs(fSystematicScalingSplines - 1.0) > fgkEpsilon) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if ((TMath::Abs(fSystematicScalingEtaCorrectionLowMomenta - 1.0) > fgkEpsilon) ||
      (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - 1.0) > fgkEpsilon)) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if (TMath::Abs(fSystematicScalingEtaSigmaPara - 1.0) > fgkEpsilon) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if (TMath::Abs(fSystematicScalingMultCorrection - 1.0) > fgkEpsilon) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
}


//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::PDGtoMCID(Int_t pdg)
{
  // Returns the corresponding AliPID index to the given pdg code.
  // Returns AliPID::kUnkown if pdg belongs to a not considered species.
  
  Int_t absPDGcode = TMath::Abs(pdg);
  if (absPDGcode == 211) {//Pion
    return AliPID::kPion;
  }
  else if (absPDGcode == 321) {//Kaon
    return AliPID::kKaon;
  }
  else if (absPDGcode == 2212) {//Proton
    return AliPID::kProton;
  }
  else if (absPDGcode == 11) {//Electron      
    return AliPID::kElectron;
  }
  else if (absPDGcode == 13) {//Muon
    return AliPID::kMuon;
  }
  
  return AliPID::kUnknown;
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::GetJetTrackObservables(const Double_t trackPt, const Double_t jetPt, Double_t& z, Double_t& xi)
{
  // Uses trackPt and jetPt to obtain z and xi.
  
  z = (jetPt > 0 && trackPt >= 0) ? (trackPt / jetPt) : -1;
  xi = (z > 0) ? TMath::Log(1. / z) : -1;
  
  if(trackPt > (1. - 1e-06) * jetPt && trackPt < (1. + 1e-06) * jetPt) { // case z=1 : move entry to last histo bin <1
    z  = 1. - 1e-06;
    xi = 1e-06;
  }
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::CleanupParticleFractionHistos()
{
  // Delete histos with particle fractions
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    delete fFractionHists[i];
    fFractionHists[i] = 0x0;
    
    delete fFractionSysErrorHists[i];
    fFractionSysErrorHists[i] = 0x0;
  }
}


//_____________________________________________________________________________
Double_t AliAnalysisTaskPID::ConvolutedGaus(const Double_t* xx, const Double_t* par) const
{
  // Convolutes gauss with an exponential tail which describes dEdx-response better than pure gaussian

  const Double_t mean = par[0];
  const Double_t sigma = par[1];
  const Double_t lambda = par[2];
  
  return lambda/sigma*TMath::Exp(-lambda/sigma*(xx[0]-mean)+lambda*lambda*0.5)*0.5*TMath::Erfc((-xx[0]+mean+sigma*lambda)/sigma*fgkOneOverSqrt2);
}


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::FastGaus(const Double_t x, const Double_t mean, const Double_t sigma) const
{
  // Calculate an unnormalised gaussian function with mean and sigma.

  if (sigma < fgkEpsilon)
    return 1.e30;
  
  const Double_t arg = (x - mean) / sigma;
  return exp(-0.5 * arg * arg);
}


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::FastNormalisedGaus(const Double_t x, const Double_t mean, const Double_t sigma) const
{
  // Calculate a normalised (divided by sqrt(2*Pi)*sigma) gaussian function with mean and sigma.

  if (sigma < fgkEpsilon)
    return 1.e30;
  
  const Double_t arg = (x - mean) / sigma;
  const Double_t res = exp(-0.5 * arg * arg);
  return res / (2.50662827463100024 * sigma); //sqrt(2*Pi)=2.50662827463100024
}


//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::FindBinWithinRange(TAxis* axis, Double_t value) const
{
  // Find the corresponding bin of the axis. Values outside the range (also under and overflow) will be set to the first/last
  // available bin
  
  Int_t bin = axis->FindFixBin(value);
  
  if (bin <= 0)
    bin = 1;
  if (bin > axis->GetNbins())
    bin = axis->GetNbins();

  return bin;
}

  
//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::FindFirstBinAboveIn3dSubset(const TH3* hist, const Double_t threshold, const Int_t yBin,
                                                      const Int_t zBin) const
{
  // Kind of projects a TH3 to 1 bin combination in y and z
  // and looks for the first x bin above a threshold for this projection.
  // If no such bin is found, -1 is returned.
  
  if (!hist)
    return -1;
    
  Int_t nBinsX = hist->GetNbinsX();
  for (Int_t xBin = 1; xBin <= nBinsX; xBin++) {
    if (hist->GetBinContent(xBin, yBin, zBin) > threshold)
      return xBin;
  }
  
  return -1;
}


//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::FindLastBinAboveIn3dSubset(const TH3* hist, const Double_t threshold, const Int_t yBin,
                                                     const Int_t zBin) const
{
  // Kind of projects a TH3 to 1 bin combination in y and z 
  // and looks for the last x bin above a threshold for this projection.
  // If no such bin is found, -1 is returned.
  
  if (!hist)
    return -1;
    
  Int_t nBinsX = hist->GetNbinsX();
  for (Int_t xBin = nBinsX; xBin >= 1; xBin--) {
    if (hist->GetBinContent(xBin, yBin, zBin) > threshold)
      return xBin;
  }
  
  return -1;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::GetParticleFraction(const Double_t trackPt, const Double_t jetPt, const Double_t centralityPercentile,
                                               const AliPID::EParticleType species,
                                               Double_t& fraction, Double_t& fractionErrorStat, Double_t& fractionErrorSys) const
{
  // Computes the particle fraction for the corresponding species for the given trackPt, jetPt and centrality.
  // Use jetPt = -1 for inclusive spectra and centralityPercentile = -1 for pp.
  // On success (return value kTRUE), fraction contains the particle fraction, fractionErrorStat(Sys) the sigma of its
  // statistical (systematic) error
  
  fraction = -999.;
  fractionErrorStat = 999.;
  fractionErrorSys = 999.;
  
  if (species > AliPID::kProton || species < AliPID::kElectron) {
    AliError(Form("Only fractions for species index %d to %d availabe, but not for the requested one: %d", 0, AliPID::kProton, species));
    return kFALSE;
  }
  
  if (!fFractionHists[species]) {
    AliError(Form("Histo with particle fractions for species %d not loaded!", species));
    
    return kFALSE;
  }
  
  Int_t jetPtBin = FindBinWithinRange(fFractionHists[species]->GetYaxis(), jetPt);
  Int_t centBin  = FindBinWithinRange(fFractionHists[species]->GetZaxis(), centralityPercentile);
  
  // The following interpolation takes the bin content of the first/last available bin,
  // if requested point lies beyond bin center of first/last bin.
  // The interpolation is only done for the x-axis (track pT), i.e. jetPtBin and centBin are fix,
  // because the analysis will anyhow be run in bins of jetPt and centrality and
  // it is not desired to correlate different jetPt bins via interpolation.
  
  // The same procedure is used for the error of the fraction
  TAxis* xAxis = fFractionHists[species]->GetXaxis();
  
  // No interpolation to values beyond the centers of the first/last bins (we don't know exactly where the spectra start or stop,
  // thus, search for the first and last bin above 0.0 to constrain the range
  Int_t firstBin = TMath::Max(1, FindFirstBinAboveIn3dSubset(fFractionHists[species], 0.0, jetPtBin, centBin));
  Int_t lastBin  = TMath::Min(fFractionHists[species]->GetNbinsX(), 
                              FindLastBinAboveIn3dSubset(fFractionHists[species], 0.0, jetPtBin, centBin));
  
  if (trackPt <= xAxis->GetBinCenter(firstBin)) {
    fraction = fFractionHists[species]->GetBinContent(firstBin, jetPtBin, centBin);
    fractionErrorStat = fFractionHists[species]->GetBinError(firstBin, jetPtBin, centBin);
    fractionErrorSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(firstBin, jetPtBin, centBin) : 0.;
  }
  else if (trackPt >= xAxis->GetBinCenter(lastBin)) {
    fraction = fFractionHists[species]->GetBinContent(lastBin, jetPtBin, centBin);
    fractionErrorStat = fFractionHists[species]->GetBinError(lastBin, jetPtBin, centBin);
    fractionErrorSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(lastBin, jetPtBin, centBin) : 0.;
  }
  else {
    Double_t x0 = 0., x1 = 0., y0 = 0., y1 = 0.;
    Double_t y0errStat = 0., y1errStat = 0., y0errSys = 0., y1errSys = 0.;
    Int_t trackPtBin = xAxis->FindBin(trackPt);
    
    // Linear interpolation between nearest neighbours in trackPt
    if (trackPt <= xAxis->GetBinCenter(trackPtBin)) {
        y0 = fFractionHists[species]->GetBinContent(trackPtBin - 1, jetPtBin, centBin);
        y0errStat = fFractionHists[species]->GetBinError(trackPtBin - 1, jetPtBin, centBin);
        y0errSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(trackPtBin - 1, jetPtBin, centBin)
                                                   : 0.;
        x0 = xAxis->GetBinCenter(trackPtBin - 1);
        y1 = fFractionHists[species]->GetBinContent(trackPtBin, jetPtBin, centBin);
        y1errStat = fFractionHists[species]->GetBinError(trackPtBin, jetPtBin, centBin);
        y1errSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(trackPtBin, jetPtBin, centBin)
                                                   : 0.;
        x1 = xAxis->GetBinCenter(trackPtBin);
    }
    else {
        y0 = fFractionHists[species]->GetBinContent(trackPtBin, jetPtBin, centBin);
        y0errStat = fFractionHists[species]->GetBinError(trackPtBin, jetPtBin, centBin);
        y0errSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(trackPtBin, jetPtBin, centBin)
                                                   : 0.;
        x0 = xAxis->GetBinCenter(trackPtBin);
        y1 = fFractionHists[species]->GetBinContent(trackPtBin + 1, jetPtBin, centBin);
        y1errStat = fFractionHists[species]->GetBinError(trackPtBin + 1, jetPtBin, centBin);
        y1errSys = fFractionSysErrorHists[species] ? fFractionSysErrorHists[species]->GetBinError(trackPtBin + 1, jetPtBin, centBin)
                                                   : 0.;
        x1 = xAxis->GetBinCenter(trackPtBin + 1);
    }
    
    // Per construction: x0 < trackPt < x1
    fraction = y0 + (trackPt - x0) * ((y1 - y0) / (x1 - x0));
    fractionErrorStat = y0errStat + (trackPt - x0) * ((y1errStat - y0errStat) / (x1 - x0));
    fractionErrorSys = fFractionSysErrorHists[species] ? (y0errSys + (trackPt - x0) * ((y1errSys - y0errSys) / (x1 - x0))) : 0.;
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::GetParticleFractions(const Double_t trackPt, const Double_t jetPt, const Double_t centralityPercentile,
                                                Double_t* prob, const Int_t smearSpeciesByError,
                                                const Int_t takeIntoAccountSpeciesSysError, const Bool_t uniformSystematicError) const
{
  // Fills the particle fractions for the given trackPt, jetPt and centrality into "prob".
  // Use jetPt = -1 for inclusive spectra and centralityPercentile = -1 for pp.
  // If smearSpeciesByError is >= 0 && < AliPID::kSPECIES, the returned fractions will be a random number distributed
  // with a gauss with mean being the corresponding particle fraction and sigma it's error for the considered species
  // "smearSpeciesByError".
  // Note that in this case the fractions for all species will NOT sum up to 1!
  // Thus, all other species fractions will be re-scaled weighted with their corresponding statistical error.
  // A similar procedure is used for "takeIntoAccountSpeciesSysError":  The systematic error of the corresponding species
  // is used to generate a random number with uniform distribution in [mean - sysError, mean + sysError] for the new mean
  // (in cace of uniformSystematicError = kTRUE, otherwise it will be a gaus(mean, sysError)),
  // then the other species will be re-scaled according to their systematic errors.
  // First, the systematic error uncertainty procedure will be performed (that is including re-scaling), then the statistical
  // uncertainty procedure.
  // On success, kTRUE is returned.
  
  if (!prob || smearSpeciesByError >= AliPID::kSPECIES || takeIntoAccountSpeciesSysError >= AliPID::kSPECIES)
    return kFALSE;
  
  Double_t probTemp[AliPID::kSPECIES];
  Double_t probErrorStat[AliPID::kSPECIES];
  Double_t probErrorSys[AliPID::kSPECIES];
  
  Bool_t success = kTRUE;
  success = success && GetParticleFraction(trackPt, jetPt, centralityPercentile, AliPID::kElectron,
                                           probTemp[AliPID::kElectron], probErrorStat[AliPID::kElectron], 
                                           probErrorSys[AliPID::kElectron]);
  success = success && GetParticleFraction(trackPt, jetPt, centralityPercentile, AliPID::kMuon,
                                           probTemp[AliPID::kMuon], probErrorStat[AliPID::kMuon], probErrorSys[AliPID::kMuon]);
  success = success && GetParticleFraction(trackPt, jetPt, centralityPercentile, AliPID::kPion,
                                           probTemp[AliPID::kPion], probErrorStat[AliPID::kPion], probErrorSys[AliPID::kPion]);
  success = success && GetParticleFraction(trackPt, jetPt, centralityPercentile, AliPID::kKaon,
                                           probTemp[AliPID::kKaon], probErrorStat[AliPID::kKaon], probErrorSys[AliPID::kKaon]);
  success = success && GetParticleFraction(trackPt, jetPt, centralityPercentile, AliPID::kProton,
                                           probTemp[AliPID::kProton], probErrorStat[AliPID::kProton], probErrorSys[AliPID::kProton]);
  
  if (!success)
    return kFALSE;
  
  // If desired, take into account the systematic error of the corresponding species and re-generate probTemp accordingly
  if (takeIntoAccountSpeciesSysError >= 0) {
    // Generate random fraction of the considered species "smearSpeciesByError" according to mean and sigma
    Double_t generatedFraction = uniformSystematicError 
                                   ? fRandom->Rndm() * 2. * probErrorSys[takeIntoAccountSpeciesSysError]
                                     - probErrorSys[takeIntoAccountSpeciesSysError]
                                     + probTemp[takeIntoAccountSpeciesSysError]
                                   : fRandom->Gaus(probTemp[takeIntoAccountSpeciesSysError],
                                                   probErrorSys[takeIntoAccountSpeciesSysError]);
    
    // Catch cases with invalid fraction (can happen for large errors), i.e. fraction < 0 or > 1
    if (generatedFraction < 0.)
      generatedFraction = 0.;
    else if (generatedFraction > 1.)
      generatedFraction = 1.;
    
    // Calculate difference from original fraction (original fractions sum up to 1!)
    Double_t deltaFraction = generatedFraction - probTemp[takeIntoAccountSpeciesSysError];
    
    // Fractions must (including errors) lie inside [0,1] -> Adapt weights accordingly by setting the errors
    if (deltaFraction > 0) {
      // Some part will be SUBTRACTED from the other fractions
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (probTemp[i] - probErrorSys[i] < 0) 
          probErrorSys[i] = probTemp[i];
      }
    }
    else {
      // Some part will be ADDED to the other fractions
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (probTemp[i] + probErrorSys[i] > 1) 
          probErrorSys[i] = 1. - probTemp[i];
      }
    }
    
    // Compute summed weight of all fractions except for the considered one
    Double_t summedWeight = 0.;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != takeIntoAccountSpeciesSysError)
        summedWeight += probErrorSys[i]; 
    }
    
    // Compute the weight for the other species
    /*
    if (summedWeight <= 1e-13) {
      // If this happens for some reason (it should not!), just assume flat weight
      printf("Error: summedWeight (sys error) ~ 0 for trackPt %f, jetPt %f, centralityPercentile %f. Setting flat weight!\n",
             trackPt, jetPt, centralityPercentile);
    }*/
    
    Double_t weight[AliPID::kSPECIES];
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != takeIntoAccountSpeciesSysError) {
        if (summedWeight > 1e-13)
          weight[i] = probErrorSys[i] / summedWeight;
        else
          weight[i] = probErrorSys[i] / (AliPID::kSPECIES - 1);
      }
    }
    
    // For the final generated fractions, set the generated value for the considered species
    // and the generated value minus delta times statistical weight
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != takeIntoAccountSpeciesSysError) 
        probTemp[i] = probTemp[i] - weight[i] * deltaFraction;
      else
        probTemp[i] = generatedFraction;
    }
  }
  
  // Using the values of probTemp (either the original ones or those after taking into account the systematic error),
  // calculate the final fractions - if the statistical error is to be taken into account, smear the corresponding
  // fraction. If not, just write probTemp to the final result array.
  if (smearSpeciesByError >= 0) {
    // Generate random fraction of the considered species "smearSpeciesByError" according to mean and sigma
    Double_t generatedFraction = fRandom->Gaus(probTemp[smearSpeciesByError], probErrorStat[smearSpeciesByError]);
    
    // Catch cases with invalid fraction (can happen for large errors), i.e. fraction < 0 or > 1
    if (generatedFraction < 0.)
      generatedFraction = 0.;
    else if (generatedFraction > 1.)
      generatedFraction = 1.;
    
    // Calculate difference from original fraction (original fractions sum up to 1!)
    Double_t deltaFraction = generatedFraction - probTemp[smearSpeciesByError];
    
    // Fractions must (including errors) lie inside [0,1] -> Adapt weights accordingly by setting the errors
    if (deltaFraction > 0) {
      // Some part will be SUBTRACTED from the other fractions
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (probTemp[i] - probErrorStat[i] < 0) 
          probErrorStat[i] = probTemp[i];
      }
    }
    else {
      // Some part will be ADDED to the other fractions
      for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
        if (probTemp[i] + probErrorStat[i] > 1) 
          probErrorStat[i] = 1. - probTemp[i];
      }
    }
    
    // Compute summed weight of all fractions except for the considered one
    Double_t summedWeight = 0.;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != smearSpeciesByError)
        summedWeight += probErrorStat[i]; 
    }
    
    // Compute the weight for the other species
    /*
    if (summedWeight <= 1e-13) {
      // If this happens for some reason (it should not!), just assume flat weight
      printf("Error: summedWeight (stat error) ~ 0 for trackPt %f, jetPt %f, centralityPercentile %f. Setting flat weight!\n",
             trackPt, jetPt, centralityPercentile);
    }*/
    
    Double_t weight[AliPID::kSPECIES];
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != smearSpeciesByError) {
        if (summedWeight > 1e-13)
          weight[i] = probErrorStat[i] / summedWeight;
        else
          weight[i] = probErrorStat[i] / (AliPID::kSPECIES - 1);
      }
    }
    
    // For the final generated fractions, set the generated value for the considered species
    // and the generated value minus delta times statistical weight
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (i != smearSpeciesByError) 
        prob[i] = probTemp[i] - weight[i] * deltaFraction;
      else
        prob[i] = generatedFraction;
    }
  }
  else {
    // Just take the generated values
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      prob[i] = probTemp[i];
  }
  
  
  // Should already be normalised, but make sure that it really is:
  Double_t probSum = 0.;
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    probSum += prob[i];
  }
  
  if (probSum <= 0)
    return kFALSE;
  
  if (TMath::Abs(probSum - 1.0) > 1e-4) {
    printf("Warning: Re-normalising sum of fractions: Sum is %e\n", probSum);
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      prob[i] /= probSum;
    }
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
const TH3D* AliAnalysisTaskPID::GetParticleFractionHisto(const Int_t species, const Bool_t sysError) const
{
  if (species < AliPID::kElectron || species > AliPID::kProton)
    return 0x0;
  
  return sysError ? fFractionSysErrorHists[species] : fFractionHists[species];
}


//_____________________________________________________________________________
Double_t AliAnalysisTaskPID::GetMCStrangenessFactorCMS(Int_t motherPDG, Double_t motherGenPt)
{
  // Strangeness ratio MC/data as function of mother pt from CMS data in |eta|<2.0
  // -> Based on function in PWGJE/AliAnalysisTaskFragmentationFunction, which uses
  // the following data from CMS pp @ 7 TeV inclusive (JHEP 05 (2011) 064)
  
  Double_t fac = 1;
  
  const Int_t absMotherPDG = TMath::Abs(motherPDG);
  
  if (absMotherPDG == 310 || absMotherPDG == 321) { // K0s / K+ / K-
    if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.768049;
    else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.732933;
    else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.650298;
    else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.571332;
    else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.518734;
    else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.492543;
    else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.482704;
    else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.488056;
    else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.488861;
    else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.492862;
    else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.504332;
    else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.501858;
    else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.512970;
    else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.524131;
    else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.539130;
    else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.554101;
    else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.560348;
    else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.568869;
    else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.583310;
    else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.604818;
    else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.632630;
    else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.710070;
    else if(6.00 <= motherGenPt && motherGenPt < 8.00) fac = 0.736365;
    else if(8.00 <= motherGenPt && motherGenPt < 10.00) fac = 0.835865;
  }
  
  if (absMotherPDG == 3122) { // Lambda
    if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.645162;
    else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.627431;
    else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.457136;
    else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.384369;
    else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.330597;
    else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.309571;
    else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.293620;
    else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.283709;
    else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.282047;
    else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.277261;
    else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.275772;
    else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.280726;
    else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.288540;
    else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.288315;
    else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.296619;
    else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.302993;
    else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.338121;
    else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.349800;
    else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.356802;
    else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.391202;
    else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.422573;
    else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.573815;
    else if(6.00 <= motherGenPt && motherGenPt < 8.00) fac = 0.786984;
    else if(8.00 <= motherGenPt && motherGenPt < 10.00) fac = 1.020021;
  } 
  
  if (absMotherPDG == 3312 || absMotherPDG == 3322) { // xi 
    if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.666620;
    else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.575908;
    else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.433198;
    else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.340901;
    else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.290896;
    else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.236074;
    else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.218681;
    else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.207763;
    else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.222848;
    else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.208806;
    else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.197275;
    else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.183645;
    else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.188788;
    else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.188282;
    else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.207442;
    else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.240388;
    else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.241916;
    else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.208276;
    else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.234550;
    else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.251689;
    else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.310204;
    else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.343492;  
  }
  
  const Double_t weight = 1. / fac;
  
  return weight;
}


//_____________________________________________________________________________
Double_t AliAnalysisTaskPID::GetMCStrangenessFactorCMS(AliMCEvent* mcEvent, AliMCParticle* daughter)
{
  // Strangeness ratio MC/data as function of mother pt from CMS data in |eta|<2.0
  // -> Based on function in PWGJE/AliAnalysisTaskFragmentationFunction

  if (!mcEvent)
    return 1.;

  AliMCParticle* currentMother   = daughter;
  AliMCParticle* currentDaughter = daughter;


  // find first primary mother K0s, Lambda or Xi   
  while(1) {
    Int_t daughterPDG = currentDaughter->PdgCode();  

    Int_t motherLabel = currentDaughter->GetMother();
    if(motherLabel >= mcEvent->GetNumberOfTracks()){ // protection
      currentMother = currentDaughter; 
      break; 
    }

    currentMother = (AliMCParticle*)mcEvent->GetTrack(motherLabel);

    if (!currentMother) { 
      currentMother = currentDaughter; 
      break; 
    }

    Int_t motherPDG = currentMother->PdgCode();  
 
    // phys. primary found ?    
    if (mcEvent->IsPhysicalPrimary(motherLabel))
      break; 

    if (TMath::Abs(daughterPDG) == 321) { 
      // K+/K- e.g. from phi (ref data not feeddown corrected)
      currentMother = currentDaughter;
      break; 
    }   
    if (TMath::Abs(motherPDG) == 310) { 
      // K0s e.g. from phi (ref data not feeddown corrected)
      break; 
    }   
    if (TMath::Abs(motherPDG) == 3212 && TMath::Abs(daughterPDG) == 3122) { 
      // Mother Sigma0, daughter Lambda (this case not included in feeddown corr.)
      currentMother = currentDaughter;
      break; 
    }

    currentDaughter = currentMother;
  }


  Int_t motherPDG = currentMother->PdgCode();  
  Double_t motherGenPt = currentMother->Pt(); 

  return GetMCStrangenessFactorCMS(motherPDG, motherGenPt);
}


// _________________________________________________________________________________
Bool_t AliAnalysisTaskPID::IsSecondaryWithStrangeMotherMC(AliMCEvent* mcEvent, Int_t partLabel)
{
  // Check whether particle is a secondary with strange mother, i.e. returns kTRUE if a strange mother is found
  // and the particle is NOT a physical primary. In all other cases kFALSE is returned
  
  if (!mcEvent || partLabel < 0)
    return kFALSE;
  
  AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(partLabel);
  
  if (!part)
    return kFALSE;
  
  if (mcEvent->IsPhysicalPrimary(partLabel))
    return kFALSE;
  
  Int_t iMother = part->GetMother();
  if (iMother < 0)
    return kFALSE;
  
  
  AliMCParticle* partM = (AliMCParticle*)mcEvent->GetTrack(iMother);
  if (!partM) 
   return kFALSE;
  
  Int_t codeM = TMath::Abs(partM->PdgCode());
  Int_t mfl = Int_t(codeM / TMath::Power(10, Int_t(TMath::Log10(codeM))));
  if (mfl == 3 && codeM != 3) // codeM = 3 is for s quark
    return kTRUE;
  
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::SetParticleFractionHisto(const TH3D* hist, const Int_t species, const Bool_t sysError)
{
  // Store a clone of hist (containing the particle fractions of the corresponding species with statistical error (sysError = kFALSE)
  // or systematic error (sysError = kTRUE), respectively), internally 
  
  if (species < AliPID::kElectron || species > AliPID::kProton) {
    AliError(Form("Only fractions for species index %d to %d can be set, but not for the requested one: %d", 0,
                  AliPID::kProton, species));
    return kFALSE;
  }
  
  if (sysError) {
    delete fFractionSysErrorHists[species];
    
    fFractionSysErrorHists[species] = new TH3D(*hist);
  }
  else {
    delete fFractionHists[species];
  
    fFractionHists[species] = new TH3D(*hist);
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::SetParticleFractionHistosFromFile(const TString filePathName, const Bool_t sysError)
{
  // Loads particle fractions for all species from the desired file and returns kTRUE on success.
  // The maps are assumed to be of Type TH3D, to sit in the main directory and to have names 
  // Form("hFraction_%e", AliPID::ParticleName(i)) for sysError = kFALSE and
  // Form("hFractionSysError_%e", AliPID::ParticleName(i)) for sysError = kTRUE.
  
  TFile* f = TFile::Open(filePathName.Data());
  if (!f)  {
    std::cout << "Failed to open file with particle fractions \"" << filePathName.Data() << "\"!" << std::endl;
    return kFALSE;
  }

  TH3D* hist = 0x0;
  for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
    TString histName = Form("hFraction%s_%s", sysError ? "SysError" : "", AliPID::ParticleName(i));
    hist = dynamic_cast<TH3D*>(f->Get(histName.Data()));
    if (!hist) {
      std::cout << "Failed to load particle fractions for " << histName.Data() << "!";
      std::cout << std::endl << "Cleaning up particle fraction histos!" << std::endl;
      CleanupParticleFractionHistos();
      return kFALSE;
    }
    
    if (!SetParticleFractionHisto(hist, i, sysError)) {
      std::cout << "Failed to load particle fractions for " << histName.Data() << "!";
      std::cout << std::endl << "Cleaning up particle fraction histos!" << std::endl;
      CleanupParticleFractionHistos();
      return kFALSE;
    }
  }
  
  delete hist;
  
  return kTRUE;
  
}


//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::GetRandomParticleTypeAccordingToParticleFractions(const Double_t trackPt, const Double_t jetPt,
                                                                            const Double_t centralityPercentile, 
                                                                            const Bool_t smearByError,
                                                                            const Bool_t takeIntoAccountSysError) const
{
  // Uses the stored histograms with the particle fractions to generate a random particle type according to these fractions.
  // In case of problems (e.g. histo missing), AliPID::kUnknown is returned.
  // If smearByError is kTRUE, the used fractions will be random numbers distributed with a gauss with mean
  // being the corresponding particle fraction and sigma it's error.
  // Note that in this case only the fraction of a random species is varied in this way. The other fractions
  // will be re-normalised according their statistical errors.
  // The same holds for the systematic error of species "takeIntoAccountSpeciesSysError", but the random number will be
  // uniformly distributed within [mean - sys, mean + sys] and the re-normalisation will be weighted with the systematic errors.
  // Note that the fractions will be calculated first with only the systematic error taken into account (if desired), including
  // re-normalisation. Then, the resulting fractions will be used to calculate the final fractions - either with statistical error
  // or without. The species, for which the error will be used for smearing, is the same for sys and stat error.
  
  Double_t prob[AliPID::kSPECIES];
  Int_t randomSpecies = (smearByError || takeIntoAccountSysError) ? (Int_t)(fRandom->Rndm() * AliPID::kSPECIES) : -1;
  Bool_t success = GetParticleFractions(trackPt, jetPt, centralityPercentile, prob, randomSpecies, randomSpecies);
  
  if (!success)
    return AliPID::kUnknown;
  
  Double_t rnd = fRandom->Rndm(); // Produce uniformly distributed floating point in ]0, 1]
  
  if (rnd <= prob[AliPID::kPion])
    return AliPID::kPion;
  else if (rnd <= prob[AliPID::kPion] + prob[AliPID::kKaon])
    return AliPID::kKaon;
  else if (rnd <= prob[AliPID::kPion] + prob[AliPID::kKaon] + prob[AliPID::kProton])
    return AliPID::kProton;
  else if (rnd <= prob[AliPID::kPion] + prob[AliPID::kKaon] + prob[AliPID::kProton] + prob[AliPID::kElectron])
    return AliPID::kElectron;
  
  return AliPID::kMuon;  //else it must be a muon (only species left)
}


//_____________________________________________________________________________
AliAnalysisTaskPID::ErrorCode AliAnalysisTaskPID::GenerateDetectorResponse(const AliAnalysisTaskPID::ErrorCode errCode, 
                                                                           const Double_t mean, const Double_t sigma,
                                                                           Double_t* responses, const Int_t nResponses, 
                                                                           const Bool_t usePureGaus)
{
  // Generate detector response. If a previous generation was not successful or there is something wrong with this signal generation,
  // the function will return kFALSE  
  if (!responses)
    return kError;
  
  // Reset response array
  for (Int_t i = 0; i < nResponses; i++)
    responses[i] = -999;
  
  if (errCode == kError)
    return kError;
  
  ErrorCode ownErrCode = kNoErrors;
  
  if (fUseConvolutedGaus && !usePureGaus) {
    // In case of convoluted gauss, calculate the probability density only once to save a lot of time!
    
    TH1* hProbDensity = 0x0;
    ownErrCode = SetParamsForConvolutedGaus(mean, sigma);
    if (ownErrCode == kError)
      return kError;
    
    hProbDensity = fConvolutedGausDeltaPrime->GetHistogram();
    
    for (Int_t i = 0; i < nResponses; i++) {
      responses[i] = hProbDensity->GetRandom();
      //responses[i] fConvolutedGausDeltaPrime->GetRandom(); // MUCH slower than using the binned version via the histogram
    }
  }
  else {
    for (Int_t i = 0; i < nResponses; i++) {
      responses[i] = fRandom->Gaus(mean, sigma);
    }
  }
  
  // If forwarded error code was a warning (error case has been handled before), return a warning
  if (errCode == kWarning)
    return kWarning;
  
  return ownErrCode; // Forward success/warning
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::PrintSettings(Bool_t printSystematicsSettings) const
{
  // Print current settings.
  
  printf("\n\nSettings for task %s:\n", GetName());
  printf("Is pPb/Pbp: %d -> %s\n", GetIsPbpOrpPb(), GetIsPbpOrpPb() ? "Adapting vertex cuts" : "Using standard vertex cuts");
  printf("Track cuts: %s\n", fTrackFilter ? fTrackFilter->GetTitle() : "-");
  printf("Eta cut: %.2f <= |eta| <= %.2f\n", GetEtaAbsCutLow(), GetEtaAbsCutUp());
  printf("Phi' cut: %d\n", GetUsePhiCut());
  
  printf("\n");
  
  printf("Centrality estimator: %s\n", GetCentralityEstimator().Data());
  
  printf("\n");
  
  printf("Use MC-ID for signal generation: %d\n", GetUseMCidForGeneration());
  printf("Use ITS: %d\n", GetUseITS());
  printf("Use TOF: %d\n", GetUseTOF());
  printf("Use priors: %d\n", GetUsePriors());
  printf("Use TPC default priors: %d\n", GetUseTPCDefaultPriors());
  printf("Use convoluted Gauss: %d\n", GetUseConvolutedGaus());
  printf("Accuracy of non-Gaussian tail: %e\n", GetAccuracyNonGaussianTail());
  printf("Take into account muons: %d\n", GetTakeIntoAccountMuons());
  printf("\nParams for transition from gauss to asymmetric shape:\n");
  printf("[0]: %e\n", GetConvolutedGaussTransitionPar(0));
  printf("[1]: %e\n", GetConvolutedGaussTransitionPar(1));
  printf("[2]: %e\n", GetConvolutedGaussTransitionPar(2));
  
  printf("\n");

  printf("Input from other task: %d\n", GetInputFromOtherTask());
  printf("Store additional jet information: %d\n", GetStoreAdditionalJetInformation());
  printf("Store centrality percentile: %d", GetStoreCentralityPercentile());
  
  if (printSystematicsSettings)
    PrintSystematicsSettings();
  else
    printf("\n\n\n");
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::PrintSystematicsSettings() const
{
  // Print current settings for systematic studies.
  
  printf("\n\nSettings for systematics for task %s:\n", GetName());
  printf("Splines:\t%f\n", GetSystematicScalingSplines());
  printf("EtaCorrMomThr:\t%f\n", GetSystematicScalingEtaCorrectionMomentumThr());
  printf("EtaCorrLowP:\t%f\n", GetSystematicScalingEtaCorrectionLowMomenta());
  printf("EtaCorrHighP:\t%f\n", GetSystematicScalingEtaCorrectionHighMomenta());
  printf("SigmaPara:\t%f\n", GetSystematicScalingEtaSigmaPara());
  printf("MultCorr:\t%f\n", GetSystematicScalingMultCorrection());
  
  printf("\n\n");
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::ProcessTrack(const AliVTrack* track, Int_t particlePDGcode, Double_t centralityPercentile,
                                        Double_t jetPt) 
{
  // Process the track (generate expected response, fill histos, etc.).
  // particlePDGcode == 0 means data. Otherwise, the corresponding MC ID will be assumed.
  
  //Printf("Debug: Task %s is starting to process track: dEdx %f, pTPC %f, eta %f, ncl %d\n", GetName(), track->GetTPCsignal(), track->GetTPCmomentum(),
  //       track->Eta(), track->GetTPCsignalN());
  
  const Bool_t isMC = (particlePDGcode == 0) ? kFALSE : kTRUE;
  
  Int_t binMC = -1;
  
  if (isMC) {
    if (TMath::Abs(particlePDGcode) == 211) {//Pion
      binMC = 3;
    }
    else if (TMath::Abs(particlePDGcode) == 321) {//Kaon
      binMC = 1;
    }
    else if (TMath::Abs(particlePDGcode) == 2212) {//Proton
      binMC = 4;
    }
    else if (TMath::Abs(particlePDGcode) == 11) {//Electron      
      binMC = 0;
    }
    else if (TMath::Abs(particlePDGcode) == 13) {//Muon
      binMC = 2;
    }
    else // In MC-ID case, set to underflow bin such that the response from this track is only used for unidentified signal generation
         // or signal generation with PID response and the track is still there (as in data) - e.g. relevant w.r.t. deuterons.
         // This is important to be as much as possible consistent with data. And the tracks can still be removed by disabling the
         // underflow bin for the projections
      binMC = -1; 
  }
  
  // Momenta
  //Double_t p = track->GetP();
  //Double_t pTPC = track->GetTPCmomentum();
  Double_t pT = track->Pt();
  
  Double_t z = -1, xi = -1;
  GetJetTrackObservables(pT, jetPt, z, xi);
  
  
  Double_t trackCharge = track->Charge();
  
  // TPC signal
  Double_t dEdxTPC = fPIDResponse->IsTunedOnData() ? fPIDResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
  
  if (dEdxTPC <= 0) {
    Printf("Skipping track with strange dEdx value: dEdx %f, pTPC %f, eta %f, ncl %d\n", track->GetTPCsignal(), track->GetTPCmomentum(),
           track->Eta(), track->GetTPCsignalN());
    return kFALSE;
  }
  
  
  
  
  Double_t dEdxEl, dEdxKa, dEdxPi, dEdxMu, dEdxPr;
  Double_t sigmaEl, sigmaKa, sigmaPi, sigmaMu, sigmaPr;
  
  if (fDoAnySystematicStudiesOnTheExpectedSignal) {
    // Get the uncorrected signal first and the corresponding correction factors.
    // Then modify the correction factors and properly recalculate the corrected dEdx
    
    // Get pure spline values for dEdx_expected, without any correction
    dEdxEl = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, kFALSE, kFALSE);
    dEdxKa = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, kFALSE, kFALSE);
    dEdxPi = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, kFALSE, kFALSE);
    dEdxMu = !fTakeIntoAccountMuons ? -1 :
                            fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, kFALSE, kFALSE);
    dEdxPr = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, kFALSE, kFALSE);
    
    // Scale splines, if desired
    if (TMath::Abs(fSystematicScalingSplines - 1.0) > fgkEpsilon) {
      dEdxEl *= fSystematicScalingSplines;
      dEdxKa *= fSystematicScalingSplines;
      dEdxPi *= fSystematicScalingSplines;
      dEdxMu *= fTakeIntoAccountMuons ? fSystematicScalingSplines : 1.;
      dEdxPr *= fSystematicScalingSplines;
    }
    
    // Get the eta correction factors for the (modified) expected dEdx
    Double_t etaCorrEl = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxEl) : 1.;
    Double_t etaCorrKa = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxKa) : 1.;
    Double_t etaCorrPi = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxPi) : 1.;
    Double_t etaCorrMu = fTakeIntoAccountMuons && !fPIDResponse->UseTPCEtaCorrection() ? 
                            fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxMu) : 1.;
    Double_t etaCorrPr = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxPr) : 1.;

    // Scale eta correction factors, if desired (and eta correction maps are to be used, otherwise it is not possible!)
    if (fPIDResponse->UseTPCEtaCorrection() &&
        (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - 1.0) > fgkEpsilon ||
         TMath::Abs(fSystematicScalingEtaCorrectionLowMomenta - 1.0) > fgkEpsilon)) {
      // Since we do not want to scale the splines with this, but only the eta variation, only scale the deviation of the correction factor!
      // E.g. if we would have a flat eta depence fixed at 1.0, we would shift the whole thing equal to shifting the splines by the same factor!
      
      
      // Due to additional azimuthal effects, there is an additional eta dependence for low momenta which is not corrected successfully so far.
      // One can assign a different (higher) systematic scale factor for this low-p region and a threshold which separates low- and high-p.
      // An ERF will be used to get (as a function of P_TPC) from one correction factor to the other within roughly 0.2 GeV/c
      Double_t usedSystematicScalingEtaCorrection = fSystematicScalingEtaCorrectionHighMomenta;
      
      if (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - fSystematicScalingEtaCorrectionLowMomenta) > fgkEpsilon) {
        const Double_t pTPC = track->GetTPCmomentum();
        const Double_t fractionHighMomentumScaleFactor = 0.5 * (1. + TMath::Erf((pTPC - fSystematicScalingEtaCorrectionMomentumThr) / 0.1));
        usedSystematicScalingEtaCorrection = fSystematicScalingEtaCorrectionLowMomenta * (1 - fractionHighMomentumScaleFactor)
                                             + fSystematicScalingEtaCorrectionHighMomenta * fractionHighMomentumScaleFactor;
      }
      
      etaCorrEl = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrEl - 1.0);
      etaCorrKa = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrKa - 1.0);
      etaCorrPi = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrPi - 1.0);
      etaCorrMu = fTakeIntoAccountMuons ? (1.0 + usedSystematicScalingEtaCorrection * (etaCorrMu - 1.0)) : 1.0;
      etaCorrPr = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrPr - 1.0);
    }
    
    // Get the multiplicity correction factors for the (modified) expected dEdx
    const Int_t currEvtMultiplicity = fPIDResponse->GetTPCResponse().GetCurrentEventMultiplicity();
    
    Double_t multiplicityCorrEl = fPIDResponse->UseTPCMultiplicityCorrection() ? fPIDResponse->GetTPCResponse().GetMultiplicityCorrectionFast(track,
                                    dEdxEl * etaCorrEl, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrKa = fPIDResponse->UseTPCMultiplicityCorrection() ? fPIDResponse->GetTPCResponse().GetMultiplicityCorrectionFast(track,
                                    dEdxKa * etaCorrKa, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrPi = fPIDResponse->UseTPCMultiplicityCorrection() ? fPIDResponse->GetTPCResponse().GetMultiplicityCorrectionFast(track,
                                    dEdxPi * etaCorrPi, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrMu = fTakeIntoAccountMuons && fPIDResponse->UseTPCMultiplicityCorrection() ? 
                                    fPIDResponse->GetTPCResponse().GetMultiplicityCorrectionFast(track, dEdxMu * etaCorrMu, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrPr = fPIDResponse->UseTPCMultiplicityCorrection() ? fPIDResponse->GetTPCResponse().GetMultiplicityCorrectionFast(track,
                                    dEdxPr * etaCorrPr, currEvtMultiplicity) : 1.;
    
    Double_t multiplicityCorrSigmaEl = fPIDResponse->UseTPCMultiplicityCorrection() ? 
                                        fPIDResponse->GetTPCResponse().GetMultiplicitySigmaCorrectionFast(dEdxEl * etaCorrEl, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrSigmaKa = fPIDResponse->UseTPCMultiplicityCorrection() ? 
                                        fPIDResponse->GetTPCResponse().GetMultiplicitySigmaCorrectionFast(dEdxKa * etaCorrKa, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrSigmaPi = fPIDResponse->UseTPCMultiplicityCorrection() ?
                                        fPIDResponse->GetTPCResponse().GetMultiplicitySigmaCorrectionFast(dEdxPi * etaCorrPi, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrSigmaMu = fTakeIntoAccountMuons && fPIDResponse->UseTPCMultiplicityCorrection() ? 
                                        fPIDResponse->GetTPCResponse().GetMultiplicitySigmaCorrectionFast(dEdxMu * etaCorrMu, currEvtMultiplicity) : 1.;
    Double_t multiplicityCorrSigmaPr = fPIDResponse->UseTPCMultiplicityCorrection() ? 
                                        fPIDResponse->GetTPCResponse().GetMultiplicitySigmaCorrectionFast(dEdxPr * etaCorrPr, currEvtMultiplicity) : 1.;
    
    // Scale multiplicity correction factors, if desired (and multiplicity correction functions are to be used, otherwise it is not possible!)
    if (fPIDResponse->UseTPCMultiplicityCorrection() && TMath::Abs(fSystematicScalingMultCorrection - 1.0) > fgkEpsilon) {
      // Since we do not want to scale the splines with this, but only the multiplicity variation, only scale the deviation of the correction factor!
      // E.g. if we would have a flat mult depence fix at 1.0, we would shift the whole thing equal to shifting the splines by the same factor!
      
      multiplicityCorrEl = 1.0 + fSystematicScalingMultCorrection * (multiplicityCorrEl - 1.0);
      multiplicityCorrKa = 1.0 + fSystematicScalingMultCorrection * (multiplicityCorrKa - 1.0);
      multiplicityCorrPi = 1.0 + fSystematicScalingMultCorrection * (multiplicityCorrPi - 1.0);
      multiplicityCorrMu = fTakeIntoAccountMuons ? (1.0 + fSystematicScalingMultCorrection * (multiplicityCorrMu - 1.0)) : 1.0;
      multiplicityCorrPr = 1.0 + fSystematicScalingMultCorrection * (multiplicityCorrPr - 1.0);
    }
    
    // eta correction must be enabled in order to use the new sigma parametrisation maps. Since this is the absolute sigma
    // for a track calculated with the unscaled paramaters, we have to devide out dEdxExpectedEtaCorrected and then need
    // to scale with the multiplicitySigmaCorrFactor * fSystematicScalingEtaSigmaPara. In the end, one has to scale with the
    // (modified) dEdx to get the absolute sigma
    // This means there is no extra parameter for the multiplicitySigmaCorrFactor, but only for the sigma map itself.
    // This is valid, since it appears only as a product. One has to assume a larger systematic shift in case of additional
    // multiplicity dependence....
    Double_t sigmaRelEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          / fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          * fSystematicScalingEtaSigmaPara * multiplicityCorrSigmaEl;
    
    Double_t sigmaRelKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          / fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          * fSystematicScalingEtaSigmaPara * multiplicityCorrSigmaKa;
    
    Double_t sigmaRelPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          / fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          * fSystematicScalingEtaSigmaPara * multiplicityCorrSigmaPi;
    
    Double_t sigmaRelMu = fTakeIntoAccountMuons ? 
                            fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                            / fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                            * fSystematicScalingEtaSigmaPara * multiplicityCorrSigmaMu
                            : 999.;
    
    Double_t sigmaRelPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          / fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                          * fSystematicScalingEtaSigmaPara * multiplicityCorrSigmaPr;
    
    // Now scale the (possibly modified) spline values with the (possibly modified) correction factors
    dEdxEl *= etaCorrEl * multiplicityCorrEl;
    dEdxKa *= etaCorrKa * multiplicityCorrKa;
    dEdxPi *= etaCorrPi * multiplicityCorrPi;
    dEdxMu *= etaCorrMu * multiplicityCorrMu;
    dEdxPr *= etaCorrPr * multiplicityCorrPr;
    
    // Finally, get the absolute sigma
    sigmaEl = sigmaRelEl * dEdxEl;
    sigmaKa = sigmaRelKa * dEdxKa;
    sigmaPi = sigmaRelPi * dEdxPi;
    sigmaMu = sigmaRelMu * dEdxMu;
    sigmaPr = sigmaRelPr * dEdxPr;
    
  }
  else {
    // No systematic studies on expected signal - just take it as it comve from the TPCPIDResponse
    dEdxEl = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    dEdxKa = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    dEdxPi = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    dEdxMu = !fTakeIntoAccountMuons ? -1 :
                            fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, 
                                                                             fPIDResponse->UseTPCEtaCorrection(),
                                                                             fPIDResponse->UseTPCMultiplicityCorrection());
    dEdxPr = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    
    sigmaEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    sigmaKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    sigmaPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    sigmaMu = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection());
    sigmaPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault, 
                                                              fPIDResponse->UseTPCEtaCorrection(),
                                                              fPIDResponse->UseTPCMultiplicityCorrection()); 
  }
  /*OLD with deltaSpecies
  Double_t deltaElectron = dEdxTPC - dEdxEl;
  Double_t deltaKaon = dEdxTPC - dEdxKa;
  Double_t deltaPion = dEdxTPC - dEdxPi;
  Double_t deltaProton = dEdxTPC - dEdxPr;
  */
  
  Double_t deltaPrimeElectron = (dEdxEl > 0) ? dEdxTPC / dEdxEl : -1;
  if (dEdxEl <= 0)  {
    Printf("Error: Expected TPC signal <= 0 for electron hypothesis");
    return kFALSE;
  }
    
  Double_t deltaPrimeKaon = (dEdxKa > 0) ? dEdxTPC / dEdxKa : -1;
  if (dEdxKa <= 0)  {
    Printf("Error: Expected TPC signal <= 0 for kaon hypothesis");
    return kFALSE;
  }
  
  Double_t deltaPrimePion = (dEdxPi > 0) ? dEdxTPC / dEdxPi : -1;
  if (dEdxPi <= 0)  {
    Printf("Error: Expected TPC signal <= 0 for pion hypothesis");
    return kFALSE;
  }
  
  Double_t deltaPrimeProton = (dEdxPr > 0) ? dEdxTPC / dEdxPr : -1;
  if (dEdxPr <= 0)  {
    Printf("Error: Expected TPC signal <= 0 for proton hypothesis");
    return kFALSE;
  }
      
  /*TODO for TOF 
  // TOF signal
  Double_t times[AliPID::kSPECIES];
  track->GetIntegratedTimes(times);
  AliTOFPIDResponse tofPIDResponse = fPIDResponse->GetTOFResponse();
  Float_t electronDeltaTOF = GetDeltaTOF(track, &tofPIDResponse, times, AliPID::kElectron);
  Float_t pionDeltaTOF = GetDeltaTOF(track, &tofPIDResponse, times, AliPID::kPion);
  Float_t kaonDeltaTOF = GetDeltaTOF(track, &tofPIDResponse, times, AliPID::kKaon);
  Float_t protonDeltaTOF = GetDeltaTOF(track, &tofPIDResponse, times, AliPID::kProton);
  */
  
  
  // Use probabilities to weigh the response generation for the different species.
  // Also determine the most probable particle type.
  Double_t prob[AliPID::kSPECIESC];
  for (Int_t i = 0; i < AliPID::kSPECIESC; i++)
    prob[i] = 0;

  fPIDcombined->ComputeProbabilities(track, fPIDResponse, prob);
  
  // Bug: One can set the number of species for PIDcombined, but PIDcombined will call PIDresponse, which writes without testing
  // the probs for kSPECIESC (including light nuclei) into the array.
  // In this case, when only kSPECIES are considered, the probabilities have to be rescaled!
  for (Int_t i = AliPID::kSPECIES; i < AliPID::kSPECIESC; i++)
    prob[i] = 0;
  
  // If muons are not to be taken into account, just set their probability to zero and normalise the remaining probabilities
  if (!fTakeIntoAccountMuons)
    prob[AliPID::kMuon] = 0;
  
  Double_t probSum = 0.;
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    probSum += prob[i];
  
  if (probSum > 0) {
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      prob[i] /= probSum;
  }
  
  if (!isMC) {
    // If there is no MC information, take the most probable species for the ID
    Float_t max = 0.;
    Int_t maxIndex = -1;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      if (prob[i] > max) {
        max = prob[i];
        maxIndex = i;
      }
    }
          
    // Translate from AliPID numbering to numbering of this class
    if (max > 0) {
      if (maxIndex == AliPID::kElectron)
        binMC = 0;
      else if (maxIndex == AliPID::kKaon)
        binMC = 1;
      else if (maxIndex == AliPID::kMuon)
        binMC = 2;
      else if (maxIndex == AliPID::kPion)
        binMC = 3;
      else if (maxIndex == AliPID::kProton)
        binMC = 4;
      else
        binMC = -1;
    }
    else {
      // Only take track into account for expectation values, if valid pid response is available.. Otherwise: Set to underflow bin.
      binMC = -1;
    }
  }
  
  /*
  //For testing: Swap p<->pT to analyse pure p-dependence => Needs to be removed later
  Double_t temp = pT;
  pT = pTPC;
  pTPC = temp;
  */
  
  
  Double_t entry[fStoreAdditionalJetInformation ? kDataNumAxes : kDataNumAxes - fgkNumJetAxes];
  entry[kDataMCID] = binMC;
  entry[kDataSelectSpecies] = 0;
  entry[kDataPt] = pT;
  entry[kDataDeltaPrimeSpecies] = deltaPrimeElectron;
  entry[kDataCentrality] = centralityPercentile;
  
  if (fStoreAdditionalJetInformation) {
    entry[kDataJetPt] = jetPt;
    entry[kDataZ] = z;
    entry[kDataXi] = xi;
  }
  
  entry[GetIndexOfChargeAxisData()] = trackCharge;
  
  /*OLD with TOF, p_TPC_Inner and p_vertex and deltaSpecies
  // MC PID, SelectSpecies, P(TPC_inner), pT, p(Vertex), deltaSpecies, deltaPrimeSpecies, deltaTOFspecies
  Double_t entry[8] = { binMC, 0, pTPC, pT, p, deltaElectron, deltaPrimeElectron, electronDeltaTOF };
  */
  
  fhPIDdataAll->Fill(entry);
  
  entry[kDataSelectSpecies] = 1;
  //OLD with deltaSpecies entry[5] = deltaKaon;
  entry[kDataDeltaPrimeSpecies] = deltaPrimeKaon;
  //TODO for TOF entry[7] = kaonDeltaTOF;
  fhPIDdataAll->Fill(entry);
    
  entry[kDataSelectSpecies] = 2;
  //OLD with deltaSpecies entry[5] = deltaPion;
  entry[kDataDeltaPrimeSpecies] = deltaPrimePion;
  //TODO for TOF entry[7] = pionDeltaTOF;
  fhPIDdataAll->Fill(entry);
    
  entry[kDataSelectSpecies] = 3;
  //OLD with deltaSpecies entry[5] = deltaProton;
  entry[kDataDeltaPrimeSpecies] = deltaPrimeProton;
  //TODO for TOF entry[7] = protonDeltaTOF;
  fhPIDdataAll->Fill(entry);
  
  
  // Construct the expected shape for the transition p -> pT
  
  Double_t genEntry[fStoreAdditionalJetInformation ? kGenNumAxes : kGenNumAxes - fgkNumJetAxes];
  genEntry[kGenMCID] = binMC;
  genEntry[kGenSelectSpecies] = 0;
  genEntry[kGenPt] = pT;
  genEntry[kGenDeltaPrimeSpecies] = -999;
  genEntry[kGenCentrality] = centralityPercentile;
  
  if (fStoreAdditionalJetInformation) {
    genEntry[kGenJetPt] = jetPt;
    genEntry[kGenZ] = z;
    genEntry[kGenXi] = xi;
  }
  
  genEntry[GetIndexOfChargeAxisGen()] = trackCharge;
  
  //OLD with deltaSpecies Double_t genEntry[5] = { binMC, 0, pT, -999, -999 }; // MC PID, SelectSpecies, pT, deltaSpecies, deltaPrimeSpecies
  
  // Generate numGenEntries "responses" with fluctuations around the expected values.
  // The higher the (transverse) momentum, the more "responses" will be generated in order not to run out of statistics too fast.
  Int_t numGenEntries = 10;
  if (pT >= 5) 
    numGenEntries *= 5;
  else if (pT >= 2)
    numGenEntries *= 2;
  
  // Jets have even worse statistics, therefore, scale numGenEntries further
  if (jetPt >= 40)
    numGenEntries *= 20;
  else if (jetPt >= 20)
    numGenEntries *= 10;
  else if (jetPt >= 10)
    numGenEntries *= 2;
  
  
  // Do not generate more entries than available in memory!
  if (numGenEntries > fgkMaxNumGenEntries)// fgkMaxNumGenEntries = 1000
    numGenEntries = fgkMaxNumGenEntries;
      
  ErrorCode errCode = kNoErrors;
  
  // Electrons
  errCode = GenerateDetectorResponse(errCode, 1.,              sigmaEl / dEdxEl, fGenRespElDeltaPrimeEl, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxEl / dEdxKa, sigmaEl / dEdxKa, fGenRespElDeltaPrimeKa, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxEl / dEdxPi, sigmaEl / dEdxPi, fGenRespElDeltaPrimePi, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxEl / dEdxPr, sigmaEl / dEdxPr, fGenRespElDeltaPrimePr, numGenEntries);

  // Kaons
  errCode = GenerateDetectorResponse(errCode, dEdxKa / dEdxEl, sigmaKa / dEdxEl, fGenRespKaDeltaPrimeEl, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, 1.,              sigmaKa / dEdxKa, fGenRespKaDeltaPrimeKa, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxKa / dEdxPi, sigmaKa / dEdxPi, fGenRespKaDeltaPrimePi, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxKa / dEdxPr, sigmaKa / dEdxPr, fGenRespKaDeltaPrimePr, numGenEntries);

  // Pions
  errCode = GenerateDetectorResponse(errCode, dEdxPi / dEdxEl, sigmaPi / dEdxEl, fGenRespPiDeltaPrimeEl, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxPi / dEdxKa, sigmaPi / dEdxKa, fGenRespPiDeltaPrimeKa, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, 1.,              sigmaPi / dEdxPi, fGenRespPiDeltaPrimePi, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxPi / dEdxPr, sigmaPi / dEdxPr, fGenRespPiDeltaPrimePr, numGenEntries);

  // Muons, if desired
  if (fTakeIntoAccountMuons) {
    errCode = GenerateDetectorResponse(errCode, dEdxMu / dEdxEl, sigmaMu / dEdxEl, fGenRespMuDeltaPrimeEl, numGenEntries);
    errCode = GenerateDetectorResponse(errCode, dEdxMu / dEdxKa, sigmaMu / dEdxKa, fGenRespMuDeltaPrimeKa, numGenEntries);
    errCode = GenerateDetectorResponse(errCode, dEdxMu / dEdxPi, sigmaMu / dEdxPi, fGenRespMuDeltaPrimePi, numGenEntries);
    errCode = GenerateDetectorResponse(errCode, dEdxMu / dEdxPr, sigmaMu / dEdxPr, fGenRespMuDeltaPrimePr, numGenEntries);
  }
  
  // Protons
  errCode = GenerateDetectorResponse(errCode, dEdxPr / dEdxEl, sigmaPr / dEdxEl, fGenRespPrDeltaPrimeEl, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxPr / dEdxKa, sigmaPr / dEdxKa, fGenRespPrDeltaPrimeKa, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, dEdxPr / dEdxPi, sigmaPr / dEdxPi, fGenRespPrDeltaPrimePi, numGenEntries);
  errCode = GenerateDetectorResponse(errCode, 1.,              sigmaPr / dEdxPr, fGenRespPrDeltaPrimePr, numGenEntries);
  
  
  /*OLD with deltaSpecies 
  // Delta
  
  // Electrons
  errCode = GenerateDetectorResponse(errCode, 0.,              sigmaEl, fGenRespElDeltaEl, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxEl - dEdxKa, sigmaEl, fGenRespElDeltaKa, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxEl - dEdxPi, sigmaEl, fGenRespElDeltaPi, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxEl - dEdxPr, sigmaEl, fGenRespElDeltaPr, numGenEntries, usePureGausForDelta);
  
  // Kaons
  errCode = GenerateDetectorResponse(errCode, dEdxKa - dEdxEl, sigmaKa, fGenRespKaDeltaEl, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, 0.,              sigmaKa, fGenRespKaDeltaKa, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxKa - dEdxPi, sigmaKa, fGenRespKaDeltaPi, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxKa - dEdxPr, sigmaKa, fGenRespKaDeltaPr, numGenEntries, usePureGausForDelta);
  
  // Pions
  errCode = GenerateDetectorResponse(errCode, dEdxPi - dEdxEl, sigmaPi, fGenRespPiDeltaEl, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxPi - dEdxKa, sigmaPi, fGenRespPiDeltaKa, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, 0.,              sigmaPi, fGenRespPiDeltaPi, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxPi - dEdxPr, sigmaPi, fGenRespPiDeltaPr, numGenEntries, usePureGausForDelta);
  
  // Muons
  errCode = GenerateDetectorResponse(errCode, dEdxMu - dEdxEl, sigmaMu, fGenRespMuDeltaEl, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxMu - dEdxKa, sigmaMu, fGenRespMuDeltaKa, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxMu - dEdxPi, sigmaMu, fGenRespMuDeltaPi, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxMu - dEdxPr, sigmaMu, fGenRespMuDeltaPr, numGenEntries, usePureGausForDelta);
  
  // Protons
  errCode = GenerateDetectorResponse(errCode, dEdxPr - dEdxEl, sigmaPr, fGenRespPrDeltaEl, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxPr - dEdxKa, sigmaPr, fGenRespPrDeltaKa, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, dEdxPr - dEdxPi, sigmaPr, fGenRespPrDeltaPi, numGenEntries, usePureGausForDelta);
  errCode = GenerateDetectorResponse(errCode, 0.,              sigmaPr, fGenRespPrDeltaPr, numGenEntries, usePureGausForDelta);
  */
  
  if (errCode != kNoErrors) {
    if (errCode == kWarning) {
      //Printf("Warning: Questionable detector response for track, most likely due to very low number of PID clusters! Debug output (dEdx_expected, sigma_expected):");
    }
    else 
      Printf("Error: Failed to generate detector response for track - skipped! Debug output (dEdx_expected, sigma_expected):");
    
    /*
    Printf("Pr: %e, %e", dEdxPr, sigmaPr);
    Printf("Pi: %e, %e", dEdxPi, sigmaPi);
    Printf("El: %e, %e", dEdxEl, sigmaEl);
    Printf("Mu: %e, %e", dEdxMu, sigmaMu);
    Printf("Ka: %e, %e", dEdxKa, sigmaKa);
    Printf("track: dEdx %f, pTPC %f, eta %f, ncl %d\n", track->GetTPCsignal(), track->GetTPCmomentum(), track->Eta(), 
           track->GetTPCsignalN());
    */
    
    if (errCode != kWarning) {
      fhSkippedTracksForSignalGeneration->Fill(track->GetTPCmomentum(), track->GetTPCsignalN());
      return kFALSE;// Skip generated response in case of error
    }
  }
  
  for (Int_t n = 0; n < numGenEntries; n++)  {
    if (!isMC || !fUseMCidForGeneration) {
      // If no MC info is available or shall not be used, use weighting with priors to generate entries for the different species
      Double_t rnd = fRandom->Rndm(); // Produce uniformly distributed floating point in ]0, 1]
      
      // Consider generated response as originating from...
      if (rnd <= prob[AliPID::kElectron])
        genEntry[kGenMCID] = 0; // ... an electron
      else if (rnd <= prob[AliPID::kElectron] + prob[AliPID::kKaon])
        genEntry[kGenMCID] = 1;  // ... a kaon
      else if (rnd <=  prob[AliPID::kElectron] + prob[AliPID::kKaon] + prob[AliPID::kMuon])
        genEntry[kGenMCID] = 2;  // ... a muon -> NOTE: prob[AliPID::kMuon] = 0 in case of fTakeIntoAccountMuons = kFALSE
      else if (rnd <= prob[AliPID::kElectron] + prob[AliPID::kKaon] + prob[AliPID::kMuon] + prob[AliPID::kPion])
        genEntry[kGenMCID] = 3;  // ... a pion
      else
        genEntry[kGenMCID] = 4;  // ... a proton
    }
    
    // Electrons
    genEntry[kGenSelectSpecies] = 0;
    //OLD with deltaSpecies genEntry[3] = fGenRespElDeltaEl[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimeEl[n];
    fhGenEl->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 1;
    //OLD with deltaSpecies genEntry[3] = fGenRespElDeltaKa[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimeKa[n];
    fhGenEl->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 2;
    //OLD with deltaSpecies genEntry[3] = fGenRespElDeltaPi[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimePi[n];
    fhGenEl->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 3;
    //OLD with deltaSpecies genEntry[3] = fGenRespElDeltaPr[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimePr[n];
    fhGenEl->Fill(genEntry);
    
    // Kaons
    genEntry[kGenSelectSpecies] = 0;
    //OLD with deltaSpecies genEntry[3] = fGenRespKaDeltaEl[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimeEl[n];
    fhGenKa->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 1;
    //OLD with deltaSpecies genEntry[3] = fGenRespKaDeltaKa[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimeKa[n];
    fhGenKa->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 2;
    //OLD with deltaSpecies genEntry[3] = fGenRespKaDeltaPi[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimePi[n];
    fhGenKa->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 3;
    //OLD with deltaSpecies genEntry[3] = fGenRespKaDeltaPr[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimePr[n];
    fhGenKa->Fill(genEntry);
    
    // Pions
    genEntry[kGenSelectSpecies] = 0;
    //OLD with deltaSpecies genEntry[3] = fGenRespPiDeltaEl[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimeEl[n];
    fhGenPi->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 1;
    //OLD with deltaSpecies genEntry[3] = fGenRespPiDeltaKa[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimeKa[n];
    fhGenPi->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 2;
    //OLD with deltaSpecies genEntry[3] = fGenRespPiDeltaPi[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimePi[n];
    fhGenPi->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 3;
    //OLD with deltaSpecies genEntry[3] = fGenRespPiDeltaPr[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimePr[n];
    fhGenPi->Fill(genEntry);
    
    if (fTakeIntoAccountMuons) {
      // Muons, if desired
      genEntry[kGenSelectSpecies] = 0;
      //OLD with deltaSpecies genEntry[3] = fGenRespMuDeltaEl[n];
      genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimeEl[n];
      fhGenMu->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 1;
      //OLD with deltaSpecies genEntry[3] = fGenRespMuDeltaKa[n];
      genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimeKa[n];
      fhGenMu->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 2;
      //OLD with deltaSpecies genEntry[3] = fGenRespMuDeltaPi[n];
      genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimePi[n];
      fhGenMu->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 3;
      //OLD with deltaSpecies genEntry[3] = fGenRespMuDeltaPr[n];
      genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimePr[n];
      fhGenMu->Fill(genEntry);
    }
    
    // Protons
    genEntry[kGenSelectSpecies] = 0;
    //OLD with deltaSpecies genEntry[3] = fGenRespPrDeltaEl[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimeEl[n];
    fhGenPr->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 1;
    //OLD with deltaSpecies genEntry[3] = fGenRespPrDeltaKa[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimeKa[n];
    fhGenPr->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 2;
    //OLD with deltaSpecies genEntry[3] = fGenRespPrDeltaPi[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimePi[n];
    fhGenPr->Fill(genEntry);
    
    genEntry[kGenSelectSpecies] = 3;
    //OLD with deltaSpecies genEntry[3] = fGenRespPrDeltaPr[n];
    genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimePr[n];
    fhGenPr->Fill(genEntry);
  }
  
  return kTRUE;
}


//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::SetConvolutedGaussLambdaParameter(Double_t lambda)
{
  // Set the lambda parameter of the convolution to the desired value. Automatically
  // calculates the parameters for the transition (restricted) gauss -> convoluted gauss.
  fConvolutedGaussTransitionPars[2] = lambda;
  
  // Save old parameters and settings of function to restore them later:
  Double_t* oldFuncParams = new Double_t[fkConvolutedGausNPar];
  for (Int_t i = 0; i < fkConvolutedGausNPar; i++)
    oldFuncParams[i] = fConvolutedGausDeltaPrime->GetParameter(i);
  Int_t oldFuncNpx = fConvolutedGausDeltaPrime->GetNpx();
  Double_t oldFuncRangeLow = 0, oldFuncRangeUp = 100;
  fConvolutedGausDeltaPrime->GetRange(oldFuncRangeLow, oldFuncRangeUp);
    
  // Choose some sufficiently large range
  const Double_t rangeStart = 0.5;
  const Double_t rangeEnd = 2.0;
  
  // To get the parameters for the transition, just choose arbitrarily input parameters for mu and sigma
  // (it makes sense to choose typical values). The ratio sigma_gauss / sigma_convolution is independent
  // of mu and as well the difference mu_gauss - mu_convolution.
  Double_t muInput = 1.0;
  Double_t sigmaInput = fgkSigmaReferenceForTransitionPars;
  
  
  // Step 1: Generate distribution with input parameters
  const Int_t nPar = 3;
  Double_t inputPar[nPar] = { muInput, sigmaInput, lambda };

  TH1D* hInput = new TH1D("hInput", "Input distribution", 2000, rangeStart, rangeEnd);

  fConvolutedGausDeltaPrime->SetParameters(inputPar);
  fConvolutedGausDeltaPrime->SetRange(rangeStart, rangeEnd);
  fConvolutedGausDeltaPrime->SetNpx(2000);

  /*OLD
  // The resolution and mean of the AliPIDResponse are determined in finite intervals
  // of ncl (also second order effects due to finite dEdx and tanTheta).
  // Take this into account for the transition parameters via assuming an approximately flat
  // distritubtion in ncl in this interval.
  // NOTE: The ncl interval should be the same as the one used for the sigma map creation!
  const Int_t nclStart = 151;
  const Int_t nclEnd = 160;
  const Int_t nclSteps = (nclEnd - nclStart) + 1;
  for (Int_t ncl = nclStart; ncl <= nclEnd; ncl++) {
    // Resolution scales with 1/sqrt(ncl)
    fConvolutedGausDeltaPrime->SetParameter(1, inputPar[1] * sqrt(nclEnd) / sqrt(ncl));
    TH1* hProbDensity = fConvolutedGausDeltaPrime->GetHistogram();
    
    for (Int_t i = 0; i < 50000000 / nclSteps; i++)  
      hInput->Fill(hProbDensity->GetRandom());
  }
  */
  
  TH1* hProbDensity = fConvolutedGausDeltaPrime->GetHistogram();
  
  for (Int_t i = 0; i < 50000000; i++)
    hInput->Fill(hProbDensity->GetRandom());
  
  // Step 2: Fit generated distribution with restricted gaussian
  Int_t maxBin = hInput->GetMaximumBin();
  Double_t max = hInput->GetBinContent(maxBin);
  
  UChar_t usedBins = 1;
  if (maxBin > 1) {
    max += hInput->GetBinContent(maxBin - 1);
    usedBins++;
  }
  if (maxBin < hInput->GetNbinsX()) {
    max += hInput->GetBinContent(maxBin + 1);
    usedBins++;
  }
  max /= usedBins;
  
  // NOTE: The range (<-> fraction of maximum) should be the same
  // as for the spline and eta maps creation
  const Double_t lowThreshold = hInput->GetXaxis()->GetBinLowEdge(hInput->FindFirstBinAbove(0.1 * max));
  const Double_t highThreshold = hInput->GetXaxis()->GetBinUpEdge(hInput->FindLastBinAbove(0.1 * max));
  
  TFitResultPtr fitResGaussFirstStep = hInput->Fit("gaus", "QNRS", "", lowThreshold, highThreshold);
  
  TFitResultPtr fitResGauss;
  
  if ((Int_t)fitResGaussFirstStep == 0) {
    TF1 fGauss("fGauss", "[0]*TMath::Gaus(x, [1], [2], 1)", rangeStart, rangeEnd);
    fGauss.SetParameter(0, fitResGaussFirstStep->GetParams()[0]);
    fGauss.SetParError(0, fitResGaussFirstStep->GetErrors()[0]);
    fGauss.SetParameter(2, fitResGaussFirstStep->GetParams()[2]);
    fGauss.SetParError(2, fitResGaussFirstStep->GetErrors()[2]);

    fGauss.FixParameter(1, fitResGaussFirstStep->GetParams()[1]);
    fitResGauss = hInput->Fit(&fGauss, "QNS", "s", rangeStart, rangeEnd);
  }
  else {
    fitResGauss = hInput->Fit("gaus", "QNRS", "same", rangeStart, rangeEnd);
  }
  //OLD TFitResultPtr fitResGauss = hInput->Fit("gaus", "QNRS", "", hInput->GetXaxis()->GetBinLowEdge(hInput->FindFirstBinAbove(0.1 * max)),
  //                                        hInput->GetXaxis()->GetBinUpEdge(hInput->FindLastBinAbove(0.1 * max)));
  
  
  // Step 3: Use parameters from gaussian fit to obtain parameters for the transition "restricted gauss" -> "convoluted gauss"
  
  // 3.1 The ratio sigmaInput / sigma_gaussFit ONLY depends on lambda (which is fixed per period) -> Calculate this first
  // for an arbitrary (here: typical) sigma. The ratio is then ~the same for ALL sigma for given lambda!
  if ((Int_t)fitResGauss != 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Gauss Fit failed!\n");
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete oldFuncParams;
    
    return kFALSE; 
  }
  
  if (fitResGauss->GetParams()[2] <= 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Sigma of gauss fit <= 0!\n");
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete oldFuncParams;
    
    return kFALSE;
  }
  
  // sigma correction factor
  fConvolutedGaussTransitionPars[1] = sigmaInput / fitResGauss->GetParams()[2];
  
  // 3.2 Now that sigma und lambda are determined, one can calculate mu by shifting the maximum to the desired position,
  // i.e. the maximum of the input distribution should coincide with that of the re-calculated distribution,
  // which is achieved by getting the same mu for the same sigma.
  // NOTE: For fixed lambda, the shift is proportional to sigma and independent of mu!
  // So, one can calculate the shift for an arbitrary fixed (here: typical)
  // sigma and then simply use this shift for any other sigma by scaling it correspondingly!!!
  
  // Mu shift correction:
  // Shift in mu (difference between mean of gauss and mean of convolution) is proportional to sigma!
  // Thus, choose a reference sigma (typical -> 0.05), such that for arbitrary sigma one can simple scale 
  // this shift correction with sigma / referenceSigma.
  fConvolutedGaussTransitionPars[0] = (fitResGauss->GetParams()[1] - muInput);
  
  
  /*Changed on 03.07.2013
  
  // Maximum of fConvolutedGausDeltaPrime should agree with maximum of input
  Double_t par[nPar] = { fitResGauss->GetParams()[1], // just as a guess of the maximum position
                         sigmaInput,
                         lambda };
                         
  fConvolutedGausDeltaPrime->SetParameters(par);
  
  Double_t maxXInput = fConvolutedGausDeltaPrime->GetMaximumX(TMath::Max(0.001, muInput - 3. * sigmaInput),
                                                              muInput + 10. * sigmaInput,
                                                              0.0001);
                                                              
  // Maximum shifts proportional to sigma and is linear in mu (~mean of gauss)
  // Maximum should be typically located within [gaussMean, gaussMean + 3 gaussSigma]. 
  // -> Larger search range for safety reasons (also: sigma and/or mean might be not 100% accurate).
  Double_t maxXconvoluted = fConvolutedGausDeltaPrime->GetMaximumX(TMath::Max(0.001, 
                                                                              fitResGauss->GetParams()[1] - 3. * fitResGauss->GetParams()[2]),
                                                                   fitResGauss->GetParams()[1] + 10. * fitResGauss->GetParams()[2],
                                                                   0.0001);
  if (maxXconvoluted <= 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Maximum of fConvolutedGausDeltaPrime <= 0!\n");
    
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete oldFuncParams;
    
    return kFALSE;
  }
  
  // maxX perfectly shifts as par[0] (scaled by sigma) -> Can shift maxX to input value.
  // Mu shift correction:
  fConvolutedGaussTransitionPars[0] = maxXconvoluted - maxXInput;
  */
  
  
  
  fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
  fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
  fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
  
  delete hInput;
  delete oldFuncParams;

  return kTRUE;
}


//_____________________________________________________________________________
AliAnalysisTaskPID::ErrorCode AliAnalysisTaskPID::SetParamsForConvolutedGaus(const Double_t gausMean, const Double_t gausSigma) 
{
  // Set parameters for convoluted gauss using parameters for a pure gaussian.
  // If SetConvolutedGaussLambdaParameter has not been called before to initialise the translation parameters,
  // some default parameters will be used and an error will show up.
  
  if (fConvolutedGaussTransitionPars[1] < -998) {
    AliError("Transition parameters not initialised! Default parameters will be used. Please call SetConvolutedGaussLambdaParameter(...) before any calculations!");
    SetConvolutedGaussLambdaParameter(2.0);
    AliError(Form("Parameters set to:\n[0]: %f\n[1]: %f\n[2]: %f\n", fConvolutedGaussTransitionPars[0],           
                  fConvolutedGaussTransitionPars[1], fConvolutedGaussTransitionPars[2]));
  }
  
  Double_t par[fkConvolutedGausNPar];
  par[2] = fConvolutedGaussTransitionPars[2];
  par[1] = fConvolutedGaussTransitionPars[1] * gausSigma;
  // maxX perfectly shifts as par[0] (scaled by sigma) -> Can shift maxX so that it sits at the right place.
  par[0] = gausMean - fConvolutedGaussTransitionPars[0] * par[1] / fgkSigmaReferenceForTransitionPars;
  
  ErrorCode errCode = kNoErrors;
  fConvolutedGausDeltaPrime->SetParameters(par);
  
  fConvolutedGausDeltaPrime->SetNpx(20); // Small value speeds up following algorithm (valid, since extrema far apart)

  // Accuracy of 10^-5 is enough to get 0.1% precise peak for MIPS w.r.t. to dEdx = 2000 of protons
  // (should boost up the algorithm, because 10^-10 is the default value!)
  Double_t maxX= fConvolutedGausDeltaPrime->GetMaximumX(TMath::Max(0.001, gausMean - 2. * gausSigma), 
                                                        gausMean + 6. * gausSigma, 1.0E-5);
  
  const Double_t maximum = fConvolutedGausDeltaPrime->Eval(maxX);
  const Double_t maximumFraction = maximum * fAccuracyNonGaussianTail;
  
  // Estimate lower boundary for subsequent search:
  Double_t lowBoundSearchBoundLow = TMath::Max(1e-4, maxX - 5. * gausSigma);
  Double_t lowBoundSearchBoundUp = maxX;
  
  Bool_t lowerBoundaryFixedAtZero = kFALSE;
  
  while (fConvolutedGausDeltaPrime->Eval(lowBoundSearchBoundLow) >= maximumFraction) {
    if (lowBoundSearchBoundLow <= 0) {
      // This should only happen to low dEdx particles with very few clusters and therefore large sigma, such that the gauss goes below zero deltaPrime
      if (maximum <= 0) { // Something is weired
        printf("Error generating signal: maximum is <= 0!\n");
        return kError;
      }
      else {
        const Double_t valueAtZero = fConvolutedGausDeltaPrime->Eval(0);
        if (valueAtZero / maximum > 0.05) {
          // Too large fraction below zero deltaPrime. Signal generation cannot be reliable in this case
          printf("Error generating signal: Too large fraction below zero deltaPrime: convGauss(0) / convGauss(max) =  %e / %e = %e!\n",
                 valueAtZero, maximum, valueAtZero / maximum);
          return kError;
        }
      }
      
      /*
      printf("Warning: LowBoundSearchBoundLow gets smaller zero -> Set left boundary to zero! Debug output: maximumFraction * fAccuracyNonGaussianTail = %e * %e = %e maxX %f, par[0] %f, par[1] %f, par[2] %f, gausMean %f, gausSigma %f\n",
             fConvolutedGausDeltaPrime->Eval(maxX), fAccuracyNonGaussianTail, maximumFraction, maxX, par[0], par[1], par[2], gausMean, gausSigma);
      */
      
      lowerBoundaryFixedAtZero = kTRUE;
      
      if (errCode != kError)
        errCode = kWarning;
      
      break;
    }
    
    lowBoundSearchBoundUp -= gausSigma;
    lowBoundSearchBoundLow -= gausSigma;
    
    if (lowBoundSearchBoundLow < 0) {
      lowBoundSearchBoundLow = 0;
      lowBoundSearchBoundUp += gausSigma;
    }
  }
  
  // Determine lower boundary inside estimated range. For small values of the maximum: Need more precision, since finer binning!
  Double_t rangeStart = lowerBoundaryFixedAtZero ? 0 :
                        fConvolutedGausDeltaPrime->GetX(maximumFraction, lowBoundSearchBoundLow, lowBoundSearchBoundUp, (maxX < 0.4) ? 1e-5 : 0.001);
  
  // .. and the same for the upper boundary
  Double_t rangeEnd = 0;
  // If distribution starts beyond upper boundary, everything ends up in the overflow bin. So, just reduce range and Npx to minimum
  if (rangeStart > fkDeltaPrimeUpLimit) {
    rangeEnd = rangeStart + 0.00001;
    fConvolutedGausDeltaPrime->SetRange(rangeStart,rangeEnd);
    fConvolutedGausDeltaPrime->SetNpx(4);
  }
  else {
    // Estimate upper boundary for subsequent search:
    Double_t upBoundSearchBoundUp = maxX + 5 * gausSigma;
    Double_t upBoundSearchBoundLow = maxX;
    while (fConvolutedGausDeltaPrime->Eval(upBoundSearchBoundUp) >= maximumFraction) {
      upBoundSearchBoundUp += gausSigma;
      upBoundSearchBoundLow += gausSigma;
    }
  
    //  For small values of the maximum: Need more precision, since finer binning!
    rangeEnd = fConvolutedGausDeltaPrime->GetX(maximumFraction, upBoundSearchBoundLow, upBoundSearchBoundUp, (maxX < 0.4) ? 1e-5 : 0.001);
    
    fConvolutedGausDeltaPrime->SetRange(rangeStart,rangeEnd);
    fConvolutedGausDeltaPrime->SetNpx(fhPIDdataAll->GetAxis(kDataDeltaPrimeSpecies)->FindBin(rangeEnd)
                                      - fhPIDdataAll->GetAxis(kDataDeltaPrimeSpecies)->FindBin(rangeStart) + 1);
    //fConvolutedGausDeltaPrime->SetNpx((rangeEnd - rangeStart) / fDeltaPrimeBinWidth + 1);
  }
  
  return errCode;
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpGenHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kGenPt, binsPt);
  hist->SetBinEdges(kGenDeltaPrimeSpecies, binsDeltaPrime);
  hist->SetBinEdges(kGenCentrality, binsCent);
  
  if (fStoreAdditionalJetInformation)
    hist->SetBinEdges(kGenJetPt, binsJetPt);
                          
  // Set axes titles
  hist->GetAxis(kGenMCID)->SetTitle("MC PID");
  hist->GetAxis(kGenMCID)->SetBinLabel(1, "e");
  hist->GetAxis(kGenMCID)->SetBinLabel(2, "K");
  hist->GetAxis(kGenMCID)->SetBinLabel(3, "#mu");
  hist->GetAxis(kGenMCID)->SetBinLabel(4, "#pi");
  hist->GetAxis(kGenMCID)->SetBinLabel(5, "p");
  
  hist->GetAxis(kGenSelectSpecies)->SetTitle("Select Species");
  hist->GetAxis(kGenSelectSpecies)->SetBinLabel(1, "e");
  hist->GetAxis(kGenSelectSpecies)->SetBinLabel(2, "K");
  hist->GetAxis(kGenSelectSpecies)->SetBinLabel(3, "#pi");
  hist->GetAxis(kGenSelectSpecies)->SetBinLabel(4, "p");
  
  hist->GetAxis(kGenPt)->SetTitle("P_{T} (GeV/c)");
  
  hist->GetAxis(kGenDeltaPrimeSpecies)->SetTitle("TPC #Delta'_{species} (arb. unit)");
  
  hist->GetAxis(kGenCentrality)->SetTitle(Form("Centrality Percentile (%s)", fCentralityEstimator.Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kGenJetPt)->SetTitle("P_{T}^{jet} (GeV/c)");
    
    hist->GetAxis(kGenZ)->SetTitle("z = P_{T}^{track} / P_{T}^{jet}");
    
    hist->GetAxis(kGenXi)->SetTitle("#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisGen())->SetTitle("Charge (e_{0})");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpGenYieldHist(THnSparse* hist, Double_t* binsPt, Double_t* binsCent, Double_t* binsJetPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kGenYieldPt, binsPt);
  hist->SetBinEdges(kGenYieldCentrality, binsCent);
  if (fStoreAdditionalJetInformation)
    hist->SetBinEdges(kGenYieldJetPt, binsJetPt);
  
  for (Int_t i = 0; i < 5; i++)
    hist->GetAxis(kGenYieldMCID)->SetBinLabel(i + 1, AliPID::ParticleLatexName(i));
  
  // Set axes titles
  hist->GetAxis(kGenYieldMCID)->SetTitle("MC PID");
  hist->GetAxis(kGenYieldPt)->SetTitle("P_{T}^{gen} (GeV/c)");
  hist->GetAxis(kGenYieldCentrality)->SetTitle(Form("Centrality Percentile (%s)", fCentralityEstimator.Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kGenYieldJetPt)->SetTitle("P_{T}^{jet, gen} (GeV/c)");
    
    hist->GetAxis(kGenYieldZ)->SetTitle("z = P_{T}^{track} / P_{T}^{jet}");
    
    hist->GetAxis(kGenYieldXi)->SetTitle("#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisGenYield())->SetTitle("Charge (e_{0})");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kDataPt, binsPt);
  hist->SetBinEdges(kDataDeltaPrimeSpecies, binsDeltaPrime);
  hist->SetBinEdges(kDataCentrality, binsCent);
  
  if (fStoreAdditionalJetInformation)
    hist->SetBinEdges(kDataJetPt, binsJetPt);
  
  // Set axes titles
  hist->GetAxis(kDataMCID)->SetTitle("MC PID");
  hist->GetAxis(kDataMCID)->SetBinLabel(1, "e");
  hist->GetAxis(kDataMCID)->SetBinLabel(2, "K");
  hist->GetAxis(kDataMCID)->SetBinLabel(3, "#mu");
  hist->GetAxis(kDataMCID)->SetBinLabel(4, "#pi");
  hist->GetAxis(kDataMCID)->SetBinLabel(5, "p");
  
  hist->GetAxis(kDataSelectSpecies)->SetTitle("Select Species");
  hist->GetAxis(kDataSelectSpecies)->SetBinLabel(1, "e");
  hist->GetAxis(kDataSelectSpecies)->SetBinLabel(2, "K");
  hist->GetAxis(kDataSelectSpecies)->SetBinLabel(3, "#pi");
  hist->GetAxis(kDataSelectSpecies)->SetBinLabel(4, "p");
  
  hist->GetAxis(kDataPt)->SetTitle("P_{T} (GeV/c)");
    
  hist->GetAxis(kDataDeltaPrimeSpecies)->SetTitle("TPC #Delta'_{species} (arb. unit)");
  
  hist->GetAxis(kDataCentrality)->SetTitle(Form("Centrality Percentile (%s)", fCentralityEstimator.Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kDataJetPt)->SetTitle("P_{T}^{jet} (GeV/c)");
  
    hist->GetAxis(kDataZ)->SetTitle("z = P_{T}^{track} / P_{T}^{jet}");
  
    hist->GetAxis(kDataXi)->SetTitle("#xi = ln(P_{T}^{jet} / P_{T}^{track})");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisData())->SetTitle("Charge (e_{0})");
  
  /*OLD with TOF, p_TPC_Inner and p_vertex
  // MC PID, SelectSpecies, P(TPC_inner), pT, p(Vertex), deltaSpecies, deltaPrimeSpecies, deltaTOFspecies
  hist->SetBinEdges(2, binsPt);
  hist->SetBinEdges(3, binsPt);
  hist->SetBinEdges(4, binsPt);
                          
  // Set axes titles
  hist->GetAxis(0)->SetTitle("MC PID");
  hist->GetAxis(0)->SetBinLabel(1, "e");
  hist->GetAxis(0)->SetBinLabel(2, "K");
  hist->GetAxis(0)->SetBinLabel(3, "#mu");
  hist->GetAxis(0)->SetBinLabel(4, "#pi");
  hist->GetAxis(0)->SetBinLabel(5, "p");
  
  hist->GetAxis(1)->SetTitle("Select Species");
  hist->GetAxis(1)->SetBinLabel(1, "e");
  hist->GetAxis(1)->SetBinLabel(2, "K");
  hist->GetAxis(1)->SetBinLabel(3, "#pi");
  hist->GetAxis(1)->SetBinLabel(4, "p");
  
  hist->GetAxis(2)->SetTitle("P_{TPC_inner} (GeV/c)");
  hist->GetAxis(3)->SetTitle("P_{T} (GeV/c)");
  hist->GetAxis(4)->SetTitle("P_{vertex} (GeV/c)");
  
  hist->GetAxis(5)->SetTitle("TPC #Delta_{species} (arb. unit)");
  
  hist->GetAxis(6)->SetTitle("TPC #Delta'_{species} (arb. unit)");
  
  hist->GetAxis(7)->SetTitle("#Delta TOF_{species} (ps)");
  */
}
