
#include <iostream>

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

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisFilter.h"
#include "AliInputEventHandler.h"

#include "AliVVertex.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"

#include "AliPPVsMultUtils.h"

#include "AliAnalysisTaskPID.h"

/*
This task collects PID output from different detectors.
Only tracks fulfilling some standard quality cuts are taken into account.
At the moment, only data from TPC and TOF is collected. But in future,
data from e.g. HMPID is also foreseen.

Contact: bhess@cern.ch
*/

ClassImp(AliAnalysisTaskPID)

const Int_t AliAnalysisTaskPID::fgkNumJetAxes = 5; // Number of additional axes for jets
const Double_t AliAnalysisTaskPID::fgkEpsilon = 1e-8; // Double_t threshold above zero
const Int_t AliAnalysisTaskPID::fgkMaxNumGenEntries = 500; // Maximum number of generated detector responses per track and delta(Prime) and associated species

const Double_t AliAnalysisTaskPID::fgkOneOverSqrt2 = 0.707106781186547462; // = 1. / TMath::Sqrt2();

const Double_t AliAnalysisTaskPID::fgkSigmaReferenceForTransitionPars = 0.05; // Reference sigma chosen to calculate transition

AliAnalysisTaskPID::EventGenerator AliAnalysisTaskPID::fgEventGenerator = AliAnalysisTaskPID::kPythia6Perugia0;

//________________________________________________________________________
AliAnalysisTaskPID::AliAnalysisTaskPID()
  : AliAnalysisTaskPIDV0base()
  , fRun(-1)
  , fPIDcombined(new AliPIDCombined())
  , fPPVsMultUtils(new AliPPVsMultUtils())
  , fInputFromOtherTask(kFALSE)
  , fDoPID(kTRUE)
  , fDoEfficiency(kTRUE)
  , fDoPtResolution(kFALSE)
  , fDoDeDxCheck(kFALSE)
  , fDoBinZeroStudy(kFALSE)
  , fStoreCentralityPercentile(kFALSE)
  , fStoreAdditionalJetInformation(kFALSE)
  , fTakeIntoAccountMuons(kFALSE)
  , fUseITS(kFALSE)
  , fUseTOF(kFALSE)
  , fStoreTOFInfo(kTRUE)
  , fUsePriors(kFALSE)
  , fTPCDefaultPriors(kFALSE)
  , fStoreCharge(kTRUE)
  , fUseMCidForGeneration(kTRUE)
  , fUseConvolutedGaus(kFALSE) 
  , fkConvolutedGausNPar(3)
  , fAccuracyNonGaussianTail(1e-8)
  , fkDeltaPrimeLowLimit(0.02)
  , fkDeltaPrimeUpLimit(40.0)
  , fConvolutedGausDeltaPrime(0x0)
  , fTOFmode(1)
  , fEtaAbsCutLow(0.0)
  , fEtaAbsCutUp(0.9)
  , fDoAnySystematicStudiesOnTheExpectedSignal(kFALSE)
  , fSystematicScalingSplinesThreshold(50.)
  , fSystematicScalingSplinesBelowThreshold(1.0)
  , fSystematicScalingSplinesAboveThreshold(1.0)
  , fSystematicScalingEtaCorrectionMomentumThr(0.35)
  , fSystematicScalingEtaCorrectionLowMomenta(1.0)
  , fSystematicScalingEtaCorrectionHighMomenta(1.0)
  , fSystematicScalingEtaSigmaParaThreshold(250.)
  , fSystematicScalingEtaSigmaParaBelowThreshold(1.0)
  , fSystematicScalingEtaSigmaParaAboveThreshold(1.0)
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
  , fDeltaPrimeAxis(0x0)
  , fhMaxEtaVariation(0x0)
  , fhEventsProcessed(0x0)
  , fhEventsTriggerSel(0x0)
  , fhEventsTriggerSelVtxCut(0x0) 
  , fhEventsProcessedNoPileUpRejection(0x0)
  , fChargedGenPrimariesTriggerSel(0x0)
  , fChargedGenPrimariesTriggerSelVtxCut(0x0)
  , fChargedGenPrimariesTriggerSelVtxCutZ(0x0)
  , fChargedGenPrimariesTriggerSelVtxCutZPileUpRej(0x0)
  , fhMCgeneratedYieldsPrimaries(0x0)
  , fh2FFJetPtRec(0x0)
  , fh2FFJetPtGen(0x0)
  , fh1Xsec(0x0)
  , fh1Trials(0x0)
  , fh1EvtsPtHardCut(0x0)
  , fContainerEff(0x0)
  , fQASharedCls(0x0)
  , fDeDxCheck(0x0)
  , fOutputContainer(0x0)
  , fQAContainer(0x0)
  , fIsUEPID(kFALSE)
  , fh2UEDensity(0x0)
  , fh1JetArea(0x0)
{
  // default Constructor
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskPID", AliLog::kInfo); 
  
  fConvolutedGausDeltaPrime = new TF1("convolutedGausDeltaPrime", this, &AliAnalysisTaskPID::ConvolutedGaus,
                                      fkDeltaPrimeLowLimit, fkDeltaPrimeUpLimit,
                                      fkConvolutedGausNPar, "AliAnalysisTaskPID", "ConvolutedGaus");
  
  // Set some arbitrary parameteres, such that the function call will not crash
  // (although it should not be called with these parameters...)
  fConvolutedGausDeltaPrime->SetParameter(0, 0);
  fConvolutedGausDeltaPrime->SetParameter(1, 1);
  fConvolutedGausDeltaPrime->SetParameter(2, 2);
  
  
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
    
    fPtResolution[i] = 0x0;
  }
}

//________________________________________________________________________
AliAnalysisTaskPID::AliAnalysisTaskPID(const char *name)
  : AliAnalysisTaskPIDV0base(name)
  , fRun(-1)
  , fPIDcombined(new AliPIDCombined())
  , fPPVsMultUtils(new AliPPVsMultUtils())
  , fInputFromOtherTask(kFALSE)
  , fDoPID(kTRUE)
  , fDoEfficiency(kTRUE)
  , fDoPtResolution(kFALSE)
  , fDoDeDxCheck(kFALSE)
  , fDoBinZeroStudy(kFALSE)
  , fStoreCentralityPercentile(kFALSE)
  , fStoreAdditionalJetInformation(kFALSE)
  , fTakeIntoAccountMuons(kFALSE)
  , fUseITS(kFALSE)
  , fUseTOF(kFALSE)
  , fStoreTOFInfo(kTRUE)
  , fUsePriors(kFALSE)
  , fTPCDefaultPriors(kFALSE)
  , fStoreCharge(kTRUE)
  , fUseMCidForGeneration(kTRUE)
  , fUseConvolutedGaus(kFALSE) 
  , fkConvolutedGausNPar(3)
  , fAccuracyNonGaussianTail(1e-8)
  , fkDeltaPrimeLowLimit(0.02)
  , fkDeltaPrimeUpLimit(40.0)
  , fConvolutedGausDeltaPrime(0x0)
  , fTOFmode(1)
  , fEtaAbsCutLow(0.0)
  , fEtaAbsCutUp(0.9)
  , fDoAnySystematicStudiesOnTheExpectedSignal(kFALSE)
  , fSystematicScalingSplinesThreshold(50.)
  , fSystematicScalingSplinesBelowThreshold(1.0)
  , fSystematicScalingSplinesAboveThreshold(1.0)
  , fSystematicScalingEtaCorrectionMomentumThr(0.35)
  , fSystematicScalingEtaCorrectionLowMomenta(1.0)
  , fSystematicScalingEtaCorrectionHighMomenta(1.0)
  , fSystematicScalingEtaSigmaParaThreshold(250.)
  , fSystematicScalingEtaSigmaParaBelowThreshold(1.0)
  , fSystematicScalingEtaSigmaParaAboveThreshold(1.0)
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
  , fDeltaPrimeAxis(0x0)
  , fhMaxEtaVariation(0x0)
  , fhEventsProcessed(0x0)
  , fhEventsTriggerSel(0x0)
  , fhEventsTriggerSelVtxCut(0x0) 
  , fhEventsProcessedNoPileUpRejection(0x0)
  , fChargedGenPrimariesTriggerSel(0x0)
  , fChargedGenPrimariesTriggerSelVtxCut(0x0)
  , fChargedGenPrimariesTriggerSelVtxCutZ(0x0)
  , fChargedGenPrimariesTriggerSelVtxCutZPileUpRej(0x0)
  , fhMCgeneratedYieldsPrimaries(0x0)
  , fh2FFJetPtRec(0x0)
  , fh2FFJetPtGen(0x0)
  , fh1Xsec(0x0)
  , fh1Trials(0x0)
  , fh1EvtsPtHardCut(0x0)
  , fContainerEff(0x0)
  , fQASharedCls(0x0)
  , fDeDxCheck(0x0)
  , fOutputContainer(0x0)
  , fQAContainer(0x0)
  , fIsUEPID(kFALSE)
  , fh2UEDensity(0x0)
  , fh1JetArea(0x0)
{
  // Constructor
  
  AliLog::SetClassDebugLevel("AliAnalysisTaskPID", AliLog::kInfo);
  
  fConvolutedGausDeltaPrime = new TF1("convolutedGausDeltaPrime", this, &AliAnalysisTaskPID::ConvolutedGaus,
                                      fkDeltaPrimeLowLimit, fkDeltaPrimeUpLimit,
                                      fkConvolutedGausNPar, "AliAnalysisTaskPID", "ConvolutedGaus");
  
  // Set some arbitrary parameteres, such that the function call will not crash
  // (although it should not be called with these parameters...)
  fConvolutedGausDeltaPrime->SetParameter(0, 0);
  fConvolutedGausDeltaPrime->SetParameter(1, 1);
  fConvolutedGausDeltaPrime->SetParameter(2, 2);
  
  
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
    
    fPtResolution[i] = 0x0;
  }
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());

  DefineOutput(1, TObjArray::Class());
  
  DefineOutput(2, AliCFContainer::Class());
  
  DefineOutput(3, TObjArray::Class());
}


//________________________________________________________________________
AliAnalysisTaskPID::~AliAnalysisTaskPID()
{
  // dtor
  
  CleanupParticleFractionHistos();
  
  delete fPIDcombined;
  fPIDcombined = 0x0;
  
  delete fPPVsMultUtils;
  fPPVsMultUtils = 0x0;
  
  delete fOutputContainer;
  fOutputContainer = 0x0;
  
  delete fQAContainer;
  fQAContainer = 0x0;

  delete fConvolutedGausDeltaPrime;
  fConvolutedGausDeltaPrime = 0x0;
  
  delete fDeltaPrimeAxis;
  fDeltaPrimeAxis = 0x0;
  
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
  
  delete fhMaxEtaVariation;
  fhMaxEtaVariation = 0x0;
  
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
  
  if (!fDoPID && !fDoDeDxCheck)
    return;
  
  if(fDebug > 1)
    printf("File: %s, Line: %d: SetUpPIDcombined\n", (char*)__FILE__, __LINE__);
  
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
  
  if(fDebug > 1)
    printf("File: %s, Line: %d: SetUpPIDcombined done\n", (char*)__FILE__, __LINE__);
}


//________________________________________________________________________
Bool_t AliAnalysisTaskPID::CalculateMaxEtaVariationMapFromPIDResponse()
{
  // Calculate the maximum deviation from unity of the eta correction factors for each row in 1/dEdx(splines)
  // from the eta correction map of the TPCPIDResponse. The result is stored in fhMaxEtaVariation.
  
  if (!fPIDResponse) {
    AliError("No PID response!");
    return kFALSE;
  }
  
  delete fhMaxEtaVariation;
  
  const TH2D* hEta = fPIDResponse->GetTPCResponse().GetEtaCorrMap();
  if (!hEta) {
    AliError("No eta correction map!");
    return kFALSE;
  }
  
  // Take binning from hEta in Y for fhMaxEtaVariation
  fhMaxEtaVariation = hEta->ProjectionY("hMaxEtaVariation");
  fhMaxEtaVariation->SetDirectory(0);
  fhMaxEtaVariation->Reset();
  
  // For each bin in 1/dEdx, loop of all tanTheta bins and find the maximum deviation from unity.
  // Store the result in fhMaxEtaVariation
  
  for (Int_t binY = 1; binY <= fhMaxEtaVariation->GetNbinsX(); binY++) {
    Double_t maxAbs = -1;
    for (Int_t binX = 1; binX <= hEta->GetNbinsX(); binX++) {
      Double_t curr = TMath::Abs(hEta->GetBinContent(binX, binY) - 1.);
      if (curr > maxAbs)
        maxAbs = curr;
    }
    
    if (maxAbs < 1e-12) {
      AliError(Form("Maximum deviation from unity is zero for 1/dEdx = %f (bin %d)", hEta->GetYaxis()->GetBinCenter(binY), binY));
      delete fhMaxEtaVariation;
      return kFALSE;
    }
    
    fhMaxEtaVariation->SetBinContent(binY, maxAbs);
  }
  
  printf("AliAnalysisTaskPID: Calculated max eta variation from map \"%s\".\n", hEta->GetTitle());
  
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskPID::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  if(fDebug > 1)
    printf("File: %s, Line: %d: UserCreateOutputObjects\n", (char*)__FILE__, __LINE__);
  
  // Setup basic things, like PIDResponse
  AliAnalysisTaskPIDV0base::UserCreateOutputObjects();
  
  if (!fPIDResponse)
    AliFatal("PIDResponse object was not created");
  
  SetUpPIDcombined();
  
  OpenFile(1);
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: UserCreateOutputObjects -> OpenFile(1) successful\n", (char*)__FILE__, __LINE__);
  
  fOutputContainer = new TObjArray(1);
  fOutputContainer->SetName(GetName());
  fOutputContainer->SetOwner(kTRUE);
  
  if (fDebug > 2)
    std::cout << "OutputContainer successfully created" << std::endl;
  
  /* Old binning (coarser in the intermediate region and a real subset of the new binning)
  const Int_t nPtBins = 68;
  Double_t binsPt[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4,  0.45,
           0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
           1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
           2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
           4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
           11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
           26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };*/

      const Int_t nPtBins = 73;
    Double_t binsPt[nPtBins + 1] = {0.,  0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
                                  0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
                                   1.,  1.1, 1.2,  1.3, 1.4,  1.5, 1.6,  1.7, 1.8,  1.9,
                                   2.,  2.1, 2.2,  2.3, 2.4,  2.5, 2.6,  2.7, 2.8,  2.9, 
                                   3.,  3.2, 3.4,  3.6, 3.8,   4., 4.5,   5., 5.5,   6.,
                                  6.5,   7.,  8.,   9., 10.,  11., 12.,  13., 14.,  15.,
                                  16.,  18., 20.,  22., 24.,  26., 28.,  30., 32.,  34., 
                                  36.,  40., 45.,  50. };

 
//   const Int_t nPtBins = 59;
//   Double_t binsPt[nPtBins + 1] = {0.01, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0};
  
  const Bool_t useITSTPCtrackletsCentEstimatorWithNewBinning = fCentralityEstimator.CompareTo("ITSTPCtracklets", TString::kIgnoreCase) == 0
                                                               && fStoreCentralityPercentile;
  
  const Int_t nCentBinsGeneral = 12;
  const Int_t nCentBinsNewITSTPCtracklets = 16;
  
  const Int_t nCentBins = useITSTPCtrackletsCentEstimatorWithNewBinning ? nCentBinsNewITSTPCtracklets : nCentBinsGeneral;

  Double_t binsCent[nCentBins+1];
  for (Int_t i = 0; i < nCentBins + 1; i++)
    binsCent[i] = -1;
  
  //-1 for pp (unless explicitely requested); 90-100 has huge electro-magnetic impurities
  Double_t binsCentV0[nCentBinsGeneral+1] = {-1, 0,  5, 10, 20, 30, 40, 50, 60, 70, 80,  90, 100 };
  
  // These centrality estimators deal with integers! This implies that the ranges are always [lowlim, uplim - 1]
  Double_t binsCentITSTPCTracklets[nCentBinsNewITSTPCtracklets+1] = { 0, 1, 4, 7, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 9999 };
  Double_t binsCentITSTPCTrackletsOldPreliminary[nCentBinsGeneral+1] = { 0, 7, 13, 20, 29, 40, 50, 60, 72, 83, 95, 105, 115 };
  
  // Special centrality binning for pp
  Double_t binsCentpp[nCentBinsGeneral+1] =   { 0, 0.01, 0.1, 1, 5, 10, 15, 20, 30, 40, 50, 70, 100};
  
  if (fCentralityEstimator.CompareTo("ITSTPCtrackletsOldPreliminaryBinning", TString::kIgnoreCase) == 0 && fStoreCentralityPercentile) {
    // Special binning for this centrality estimator; but keep number of bins!
    for (Int_t i = 0; i < nCentBinsGeneral+1; i++)
      binsCent[i] = binsCentITSTPCTrackletsOldPreliminary[i];
  }
  else if (fCentralityEstimator.CompareTo("ITSTPCtracklets", TString::kIgnoreCase) == 0 && fStoreCentralityPercentile) {
    // Special binning for this centrality estimator and different number of bins!
    for (Int_t i = 0; i < nCentBinsNewITSTPCtracklets+1; i++)
      binsCent[i] = binsCentITSTPCTracklets[i];
  }
  else if (fCentralityEstimator.Contains("ppMult", TString::kIgnoreCase) && fStoreCentralityPercentile) {
    // Special binning for this pp centrality estimator; but keep number of bins!
    for (Int_t i = 0; i < nCentBinsGeneral+1; i++)
      binsCent[i] = binsCentpp[i];
  }
  else {
    // Take default binning for VZERO
    for (Int_t i = 0; i < nCentBinsGeneral+1; i++)
      binsCent[i] = binsCentV0[i];
  }

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
  
  fDeltaPrimeAxis = new TAxis(deltaPrimeNBins, deltaPrimeBins);
  
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
  
  const Int_t nTOFpidInfoBins = kNumTOFpidInfoBins;
  const Double_t tofPIDinfoMin = kNoTOFinfo;
  const Double_t tofPIDinfoMax = kNoTOFinfo + kNumTOFpidInfoBins;
  
  const Int_t nDistanceBins = 30;
  const Double_t distanceBinsMin = 0.;
  const Double_t distanceBinsMax = 0.6;
  
  // jT binning - to have binning down to zero and log binning at the same time,
  // use first bin from zero extending to real start of log binning
  const Int_t nJtBins = 30 + 1;
  Double_t binsJt[nJtBins + 1];

  const Double_t fromLowJt = 0.05;
  const Double_t toHighJt = 10.;
  const Double_t factorJt = TMath::Power(toHighJt/fromLowJt, 1./(nJtBins-1));

  // Log binning for whole jT range
  binsJt[0] = 0.;
  binsJt[1] = fromLowJt;
  for (Int_t i = 0 + 1 + 1; i <= nJtBins; i++) {
    binsJt[i] = factorJt * binsJt[i - 1];
  }
  
  
  // MC PID, SelectSpecies, pT, deltaPrimeSpecies, centrality percentile, jet pT, z = track_pT/jet_pT, xi = log(1/z)
  Int_t binsNoJets[nBinsNoJets] =    { nMCPIDbins,
                                       nSelSpeciesBins,
                                       nPtBins,
                                       deltaPrimeNBins,
                                       nCentBins,
                                       nChargeBins,
                                       nTOFpidInfoBins };
  
  Int_t binsJets[nBinsJets]     =    { nMCPIDbins,
                                       nSelSpeciesBins,
                                       nPtBins,
                                       deltaPrimeNBins,
                                       nCentBins,
                                       nJetPtBins,
                                       nZBins,
                                       nXiBins,
                                       nChargeBins,
                                       nTOFpidInfoBins,
                                       nDistanceBins,
                                       nJtBins };
  
  Int_t *bins = fStoreAdditionalJetInformation ? &binsJets[0] : &binsNoJets[0];
  
  Double_t xminNoJets[nBinsNoJets] = { mcPIDmin,  
                                       selSpeciesMin,
                                       binsPt[0],
                                       deltaPrimeBins[0],
                                       binsCent[0],
                                       binsCharge[0],
                                       tofPIDinfoMin };
  
  Double_t xminJets[nBinsJets] =     { mcPIDmin,
                                       selSpeciesMin,
                                       binsPt[0],
                                       deltaPrimeBins[0],
                                       binsCent[0],
                                       binsJetPt[0],
                                       zMin,
                                       xiMin,
                                       binsCharge[0],
                                       tofPIDinfoMin,
                                       distanceBinsMin,
                                       binsJt[0] };
  
  Double_t *xmin = fStoreAdditionalJetInformation? &xminJets[0] : &xminNoJets[0];

  Double_t xmaxNoJets[nBinsNoJets] = { mcPIDmax,
                                       selSpeciesMax,
                                       binsPt[nPtBins],
                                       deltaPrimeBins[deltaPrimeNBins], 
                                       binsCent[nCentBins],
                                       binsCharge[nChargeBins],
                                       tofPIDinfoMax };
  
  Double_t xmaxJets[nBinsJets] =     { mcPIDmax,
                                       selSpeciesMax,
                                       binsPt[nPtBins],
                                       deltaPrimeBins[deltaPrimeNBins], 
                                       binsCent[nCentBins],
                                       binsJetPt[nJetPtBins],
                                       zMax,
                                       xiMax,
                                       binsCharge[nChargeBins],
                                       tofPIDinfoMax,
                                       distanceBinsMax,
                                       binsJt[nJtBins] };
  
  Double_t *xmax = fStoreAdditionalJetInformation? &xmaxJets[0] : &xmaxNoJets[0];
  
  fConvolutedGausDeltaPrime->SetNpx(deltaPrimeNBins);

  if (fDoPID) {
    fhPIDdataAll = new THnSparseD("hPIDdataAll","", nBins, bins, xmin, xmax);
    SetUpHist(fhPIDdataAll, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhPIDdataAll);
  }
  
  // Generated histograms (so far, bins are the same as for primary THnSparse)
  const Int_t nGenBins = fStoreAdditionalJetInformation ? nBinsJets : nBinsNoJets;
  // MC PID, SelectSpecies, Pt, deltaPrimeSpecies, jet pT, z = track_pT/jet_pT, xi = log(1/z)
  
  Int_t *genBins = fStoreAdditionalJetInformation ? &binsJets[0] : &binsNoJets[0];
  Double_t *genXmin = fStoreAdditionalJetInformation? &xminJets[0] : &xminNoJets[0];
  Double_t *genXmax = fStoreAdditionalJetInformation? &xmaxJets[0] : &xmaxNoJets[0];

  if (fDoPID) {
    fhGenEl = new THnSparseD("hGenEl", "", nGenBins, genBins, genXmin, genXmax);
    SetUpGenHist(fhGenEl, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhGenEl);
    
    fhGenKa = new THnSparseD("hGenKa", "", nGenBins, genBins, genXmin, genXmax);
    SetUpGenHist(fhGenKa, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhGenKa);
    
    fhGenPi = new THnSparseD("hGenPi", "", nGenBins, genBins, genXmin, genXmax);
    SetUpGenHist(fhGenPi, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhGenPi);
    
    if (fTakeIntoAccountMuons) {
      fhGenMu = new THnSparseD("hGenMu", "", nGenBins, genBins, genXmin, genXmax);
      SetUpGenHist(fhGenMu, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
      fOutputContainer->Add(fhGenMu);
    }
    
    fhGenPr = new THnSparseD("hGenPr", "", nGenBins, genBins, genXmin, genXmax);
    SetUpGenHist(fhGenPr, binsPt, deltaPrimeBins, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhGenPr);
  }
  
    if (fDebug > 2)
    std::cout << "Adding Event Histograms" << std::endl;
  
  fhEventsProcessed = new TH1D("fhEventsProcessed",
                               "Number of events passing trigger selection, vtx and zvtx cuts and pile-up rejection;Centrality Percentile", 
                               nCentBins, binsCent);
  fhEventsProcessed->Sumw2();
  fOutputContainer->Add(fhEventsProcessed);
  
  fhEventsTriggerSelVtxCut = new TH1D("fhEventsTriggerSelVtxCut",
                                      "Number of events passing trigger selection and vtx cut;Centrality Percentile", 
                                      nCentBins, binsCent);
  fhEventsTriggerSelVtxCut->Sumw2();
  fOutputContainer->Add(fhEventsTriggerSelVtxCut);
  
  fhEventsTriggerSel = new TH1D("fhEventsTriggerSel",
                                "Number of events passing trigger selection;Centrality Percentile", 
                                nCentBins, binsCent);
  fOutputContainer->Add(fhEventsTriggerSel);
  fhEventsTriggerSel->Sumw2();
  
  
  fhEventsProcessedNoPileUpRejection = new TH1D("fhEventsProcessedNoPileUpRejection",
                                                "Number of events passing trigger selection, vtx and zvtx cuts;Centrality Percentile", 
                                                nCentBins, binsCent);
  fOutputContainer->Add(fhEventsProcessedNoPileUpRejection);
  fhEventsProcessedNoPileUpRejection->Sumw2();
  
  if (fDebug > 2)
    std::cout << "Event Histograms added" << std::endl;  
  
  // Generated yields within acceptance
  const Int_t nBinsGenYields = fStoreAdditionalJetInformation ? kGenYieldNumAxes : kGenYieldNumAxes - 5;
  Int_t genYieldsBins[kGenYieldNumAxes]    = { nMCPIDbins,         nPtBins,           nCentBins,            nJetPtBins, nZBins, nXiBins,
                                               nChargeBins,   nDistanceBins,   nJtBins };
  genYieldsBins[GetIndexOfChargeAxisGenYield()] = nChargeBins;
  Double_t genYieldsXmin[kGenYieldNumAxes] = {   mcPIDmin,       binsPt[0],         binsCent[0],          binsJetPt[0],   zMin,   xiMin,
                                               binsCharge[0], distanceBinsMin, binsJt[0] };
  genYieldsXmin[GetIndexOfChargeAxisGenYield()] = binsCharge[0];
  Double_t genYieldsXmax[kGenYieldNumAxes] = {   mcPIDmax, binsPt[nPtBins], binsCent[nCentBins], binsJetPt[nJetPtBins],   zMax,   xiMax, 
                                               binsCharge[nChargeBins], distanceBinsMax, binsJt[nJtBins] };
  genYieldsXmax[GetIndexOfChargeAxisGenYield()] = binsCharge[nChargeBins];
  
  if (fDoPID) {
    fhMCgeneratedYieldsPrimaries = new THnSparseD("fhMCgeneratedYieldsPrimaries", 
                                                  "Generated yields w/o reco and cuts inside acceptance (physical primaries)", 
                                                  nBinsGenYields, genYieldsBins, genYieldsXmin, genYieldsXmax);
    SetUpGenYieldHist(fhMCgeneratedYieldsPrimaries, binsPt, binsCent, binsJetPt, binsJt);
    fOutputContainer->Add(fhMCgeneratedYieldsPrimaries);
  }
  
  if (fIsUEPID) {
    fh2UEDensity = new TH2D("fh2UEDensity", "p_{T} density of the Underlying event;Centrality Percentile;UE p_{T}/Event", nCentBins, binsCent, 10, 0.0, 4.0);
    fOutputContainer->Add(fh2UEDensity);
    fh1JetArea = new TH1D("fh1JetArea", "Jet Area (real)/#Pi*#Rho^{2} of the Jets used in the UE;Centrality Percentile",nCentBins,binsCent);
    fOutputContainer->Add(fh1JetArea);
  }
  
  // Container with several process steps (generated and reconstructed level with some variations)
  if (fDoEfficiency) {
    
    if (fDebug > 2)
    std::cout << "Try OpenFile(2)" << std::endl;  
    OpenFile(2);
    
    if(fDebug > 2)
    printf("File: %s, Line: %d: UserCreateOutputObjects -> OpenFile(2) successful\n", (char*)__FILE__, __LINE__);
  
    // Array for the number of bins in each dimension
    // Dimensions: MC-ID, trackPt, trackEta, trackCharge, cenrality percentile, jetPt, z, xi, distance, jT
    const Int_t nEffDims = fStoreAdditionalJetInformation ? kEffNumAxes : kEffNumAxes - 5; // Number of dimensions for the efficiency
    
    const Int_t nMCIDbins = AliPID::kSPECIES;
    Double_t binsMCID[nMCIDbins + 1];
    
    for(Int_t i = 0; i <= nMCIDbins; i++) {
      binsMCID[i]= i; 
    }
    
    const Int_t nEtaBins = 18;
    const Double_t binsEta[nEtaBins+1] = {-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
                                          0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    
    const Int_t nEffBins[kEffNumAxes] = { nMCIDbins, nPtBins, nEtaBins, nChargeBins, nCentBins, nJetPtBins, nZBins, nXiBins,
                                          nDistanceBins, nJtBins };
    
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
      fContainerEff->SetBinLimits(kEffDistance, distanceBinsMin, distanceBinsMax);
      fContainerEff->SetBinLimits(kEffJt, binsJt);
    }
    
    fContainerEff->SetVarTitle(kEffMCID,"MC ID");
    fContainerEff->SetVarTitle(kEffTrackPt,"p_{T} (GeV/c)");
    fContainerEff->SetVarTitle(kEffTrackEta,"#eta");
    fContainerEff->SetVarTitle(kEffTrackCharge,"Charge (e_{0})");
    fContainerEff->SetVarTitle(kEffCentrality, "Centrality Percentile");
    if (fStoreAdditionalJetInformation) {
      fContainerEff->SetVarTitle(kEffJetPt, "p_{T}^{jet} (GeV/c)");
      fContainerEff->SetVarTitle(kEffZ, "z = p_{T}^{track} / p_{T}^{jet}");
      fContainerEff->SetVarTitle(kEffXi, "#xi = ln(p_{T}^{jet} / p_{T}^{track})");
      fContainerEff->SetVarTitle(kEffDistance, "R");
      fContainerEff->SetVarTitle(kEffJt, "j_{T} (GeV/c)");
    }
    
    // Define clean MC sample
    fContainerEff->SetStepTitle(kStepGenWithGenCuts, "Particle level, cuts on particle level");
    // For Acceptance x Efficiency correction of primaries
    fContainerEff->SetStepTitle(kStepRecWithGenCuts, "Detector level (rec) with cuts on particle level"); 
    // For (pT) resolution correction
    fContainerEff->SetStepTitle(kStepRecWithGenCutsMeasuredObs, 
                                "Detector level (rec) with cuts on particle level with measured observables");
    // For secondary correction
    fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObs, 
                                "Detector level, all cuts on detector level with measured observables");
    fContainerEff->SetStepTitle(kStepRecWithRecCutsPrimaries, 
                                "Detector level, all cuts on detector level, only MC primaries");
    fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObsPrimaries, 
                                "Detector level, all cuts on detector level with measured observables, only MC primaries");
    fContainerEff->SetStepTitle(kStepRecWithRecCutsMeasuredObsStrangenessScaled, 
                                "Detector level (strangeness scaled), all cuts on detector level with measured observables");
  }
  
  if (fDoPID || fDoEfficiency) {
    // Generated jets
    fh2FFJetPtRec = new TH2D("fh2FFJetPtRec", "Number of reconstructed jets;Centrality Percentile;p_{T}^{jet} (GeV/c)",
                             nCentBins, binsCent, nJetPtBins, binsJetPt);
    fh2FFJetPtRec->Sumw2();
    fOutputContainer->Add(fh2FFJetPtRec);
    fh2FFJetPtGen = new TH2D("fh2FFJetPtGen", "Number of generated jets;Centrality Percentile;p_{T}^{jet} (GeV/c)",
                             nCentBins, binsCent, nJetPtBins, binsJetPt);
    fh2FFJetPtGen->Sumw2();
    fOutputContainer->Add(fh2FFJetPtGen);
  }
  
  // Pythia information
  fh1Xsec = new TProfile("fh1Xsec", "xsec from pyxsec.root", 1, 0, 1);
  fh1Xsec->Sumw2();
  fh1Xsec->GetXaxis()->SetBinLabel(1, "<#sigma>");
  fh1Trials = new TH1D("fh1Trials", "trials from pyxsec.root", 1, 0, 1);
  fh1Trials->Sumw2();
  fh1Trials->GetXaxis()->SetBinLabel(1, "#sum{ntrials}");
  
  fh1EvtsPtHardCut = new TH1F("fh1EvtsPtHardCut", "#events before and after MC #it{p}_{T,hard} cut;;Events",2,0,2);
  fh1EvtsPtHardCut->Sumw2();
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(1, "All");
  fh1EvtsPtHardCut->GetXaxis()->SetBinLabel(2, "#it{p}_{T,hard}");
  
  fOutputContainer->Add(fh1Xsec);
  fOutputContainer->Add(fh1Trials);
  fOutputContainer->Add(fh1EvtsPtHardCut);
  
  if (fDoDeDxCheck || fDoPtResolution) {
    OpenFile(3);
    fQAContainer = new TObjArray(1);
    fQAContainer->SetName(Form("%s_QA", GetName()));
    fQAContainer->SetOwner(kTRUE);
    
    if(fDebug > 2)
      printf("File: %s, Line: %d: UserCreateOutputObjects -> OpenFile(3) successful\n", (char*)__FILE__, __LINE__);
  }
  
  if (fDoPtResolution) {
    const Int_t nPtBinsRes = 100;
    Double_t pTbinsRes[nPtBinsRes + 1];

    const Double_t fromLowPtRes = 0.15;
    const Double_t toHighPtRes = 50.;
    const Double_t factorPtRes = TMath::Power(toHighPtRes/fromLowPtRes, 1./nPtBinsRes);
    // Log binning for whole pT range
    pTbinsRes[0] = fromLowPtRes;
    for (Int_t i = 0 + 1; i <= nPtBinsRes; i++) {
      pTbinsRes[i] = factorPtRes * pTbinsRes[i - 1];
    }
    
    const Int_t nBinsPtResolution = kPtResNumAxes;
    Int_t ptResolutionBins[kPtResNumAxes]    = { nJetPtBins,                       nPtBinsRes,            nPtBinsRes,             
                                                 nChargeBins,             nCentBins };
    Double_t ptResolutionXmin[kPtResNumAxes] = { binsJetPt[0],                   pTbinsRes[0],          pTbinsRes[0],           
                                                 binsCharge[0],           binsCent[0] };
    Double_t ptResolutionXmax[kPtResNumAxes] = { binsJetPt[nJetPtBins], pTbinsRes[nPtBinsRes], pTbinsRes[nPtBinsRes], 
                                                 binsCharge[nChargeBins], binsCent[nCentBins] };

    for (Int_t i = 0; i < AliPID::kSPECIES; i++) {
      fPtResolution[i] = new THnSparseD(Form("fPtResolution_%s", AliPID::ParticleShortName(i)), 
                                        Form("Pt resolution for primaries, %s", AliPID::ParticleLatexName(i)),
                                        nBinsPtResolution, ptResolutionBins, ptResolutionXmin, ptResolutionXmax);
      SetUpPtResHist(fPtResolution[i], pTbinsRes, binsJetPt, binsCent);
      fQAContainer->Add(fPtResolution[i]);
    }
    
    
    // Besides the pT resolution, also perform check on shared clusters
    const Int_t nBinsQASharedCls = kQASharedClsNumAxes;
    Int_t qaSharedClsBins[kQASharedClsNumAxes]    = { nJetPtBins, nPtBinsRes, 160, 160 };
    Double_t qaSharedClsXmin[kQASharedClsNumAxes] = { binsJetPt[0], pTbinsRes[0], 0, -1 };
    Double_t qaSharedClsXmax[kQASharedClsNumAxes] = { binsJetPt[nJetPtBins], pTbinsRes[nPtBinsRes], 160, 159 };
    
    fQASharedCls = new THnSparseD("fQASharedCls", "QA shared clusters", nBinsQASharedCls, qaSharedClsBins, qaSharedClsXmin, qaSharedClsXmax);
    
    SetUpSharedClsHist(fQASharedCls, pTbinsRes, binsJetPt);
    fQAContainer->Add(fQASharedCls);
  }
  
  
  
  if (fDoDeDxCheck) {
    const Int_t nEtaAbsBins = 9;
    const Double_t binsEtaAbs[nEtaAbsBins+1] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };
    
    const Double_t dEdxMin = 20;
    const Double_t dEdxMax = 110;
    const Int_t nDeDxBins = (Int_t) ((dEdxMax - dEdxMin) / 0.02);
    const Int_t nBinsDeDxCheck= kDeDxCheckNumAxes;
    Int_t dEdxCheckBins[kDeDxCheckNumAxes]    = { nSelSpeciesBins, nPtBins, nJetPtBins, nEtaAbsBins, nDeDxBins };
    Double_t dEdxCheckXmin[kDeDxCheckNumAxes] = { selSpeciesMin, binsPt[0], binsJetPt[0], binsEtaAbs[0], dEdxMin };
    Double_t dEdxCheckXmax[kDeDxCheckNumAxes] = { selSpeciesMax, binsPt[nPtBins], binsJetPt[nJetPtBins], binsEtaAbs[nEtaAbsBins], dEdxMax };
    
    fDeDxCheck = new THnSparseD("fDeDxCheck", "dEdx check", nBinsDeDxCheck, dEdxCheckBins, dEdxCheckXmin, dEdxCheckXmax);
    SetUpDeDxCheckHist(fDeDxCheck, binsPt, binsJetPt, binsEtaAbs);
    fQAContainer->Add(fDeDxCheck);
  }
  
  if (fDoBinZeroStudy) {
    const Double_t etaLow = -0.9;
    const Double_t etaUp = 0.9;
    const Int_t nEtaBins = 18;
    
    const Int_t nBinsBinZeroStudy = kBinZeroStudyNumAxes;
    Int_t binZeroStudyBins[nBinsBinZeroStudy]    = { nCentBins,                   nPtBins, nEtaBins };
    Double_t binZeroStudyXmin[nBinsBinZeroStudy] = { binsCent[0],               binsPt[0], etaLow   };
    Double_t binZeroStudyXmax[nBinsBinZeroStudy] = { binsCent[nCentBins], binsPt[nPtBins], etaUp    };
    
    fChargedGenPrimariesTriggerSel = new THnSparseD("fChargedGenPrimariesTriggerSel", "Trigger sel.", nBinsBinZeroStudy, binZeroStudyBins,
                                                    binZeroStudyXmin, binZeroStudyXmax);
    SetUpBinZeroStudyHist(fChargedGenPrimariesTriggerSel, binsCent, binsPt);
    fOutputContainer->Add(fChargedGenPrimariesTriggerSel);
    
    fChargedGenPrimariesTriggerSelVtxCut = new THnSparseD("fChargedGenPrimariesTriggerSelVtxCut", "Vertex cut", nBinsBinZeroStudy,
                                                          binZeroStudyBins, binZeroStudyXmin, binZeroStudyXmax);
    SetUpBinZeroStudyHist(fChargedGenPrimariesTriggerSelVtxCut, binsCent, binsPt);
    fOutputContainer->Add(fChargedGenPrimariesTriggerSelVtxCut);
    
    fChargedGenPrimariesTriggerSelVtxCutZ = new THnSparseD("fChargedGenPrimariesTriggerSelVtxCutZ", "Vertex #it{z} cut", nBinsBinZeroStudy,
                                                          binZeroStudyBins, binZeroStudyXmin, binZeroStudyXmax);
    SetUpBinZeroStudyHist(fChargedGenPrimariesTriggerSelVtxCutZ, binsCent, binsPt);
    fOutputContainer->Add(fChargedGenPrimariesTriggerSelVtxCutZ);
    
    fChargedGenPrimariesTriggerSelVtxCutZPileUpRej = new THnSparseD("fChargedGenPrimariesTriggerSelVtxCutZPileUpRej", "Vertex #it{z} cut", 
                                                                    nBinsBinZeroStudy, binZeroStudyBins, binZeroStudyXmin, binZeroStudyXmax);
    SetUpBinZeroStudyHist(fChargedGenPrimariesTriggerSelVtxCutZPileUpRej, binsCent, binsPt);
    fOutputContainer->Add(fChargedGenPrimariesTriggerSelVtxCutZPileUpRej);
  }
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: UserCreateOutputObjects -> Posting output data\n", (char*)__FILE__, __LINE__);
  
  PostData(1, fOutputContainer);
  if (fDoEfficiency)
    PostData(2, fContainerEff);
  if (fDoDeDxCheck || fDoPtResolution) 
    PostData(3, fQAContainer);
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: UserCreateOutputObjects -> Done\n", (char*)__FILE__, __LINE__);
}


//________________________________________________________________________
void AliAnalysisTaskPID::UserExec(Option_t *)
{
  // Main loop
  // Called for each event

  if(fDebug > 1)
    printf("File: %s, Line: %d: UserExec\n", (char*)__FILE__, __LINE__);
  
  // No processing of event, if input is fed in directly from another task
  if (fInputFromOtherTask)
    return;
  
  if(fDebug > 1)
    printf("File: %s, Line: %d: UserExec -> Processing started\n", (char*)__FILE__, __LINE__);

  fEvent = dynamic_cast<AliVEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: fEvent not available");
    return;
  }
  
  ConfigureTaskForCurrentEvent(fEvent);
  
  fMC = dynamic_cast<AliMCEvent*>(MCEvent());
  
  if (!fPIDResponse || !fPIDcombined)
    return;
  
  Double_t centralityPercentile = -1;
  Double_t centralityPercentileNoEventSelection = -1; 
  
  if (fStoreCentralityPercentile) {
    if (fCentralityEstimator.Contains("ITSTPCtracklets", TString::kIgnoreCase)) {
      // Special pp centrality estimator
      centralityPercentile = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fEvent, kTRUE);
      centralityPercentileNoEventSelection = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fEvent, kFALSE);
      //centralityPercentile = AliESDtrackCuts::GetReferenceMultiplicity(esdEvent, AliESDtrackCuts::kTrackletsITSTPC, fEtaAbsCutUp);// NOTE: Needs esd event!
    }
    else if (fCentralityEstimator.Contains("ppMult", TString::kIgnoreCase)) {
      // Another special pp centrality estimator
      centralityPercentile = fAnaUtils->GetMultiplicityPercentile(fEvent, GetPPCentralityEstimator().Data());
      centralityPercentileNoEventSelection = fPPVsMultUtils->GetMultiplicityPercentile(fEvent, GetPPCentralityEstimator().Data(), kFALSE);
    }
    else {
      // Ordinary centrality estimator
      centralityPercentile = fEvent->GetCentrality()->GetCentralityPercentile(fCentralityEstimator.Data());              //.Data() converts to const char
      centralityPercentileNoEventSelection = centralityPercentile; // Event selection not really implemented for this....
    }
  }
  
  const Bool_t nonNegativeCentralityPercentileNoEventSelection = centralityPercentileNoEventSelection >= 0;
  const Bool_t nonNegativeCentralityPercentile = centralityPercentile >= 0;
  
  // MB
  Bool_t passedVertexSelectionMult, passedVertexZSelectionMult;
  Bool_t isPileUpMult;
  Bool_t passedDAQCheck, passedTrackClustCut;
  
  //Checks for all run modes. Checks are different->Done in Functions
  // Check if vertex is ok, but don't apply cut on z position
  const Bool_t passedVertexSelectionMB = GetVertexIsOk(fEvent, kFALSE);
  // Now check again, but also require z position to be in desired range
  const Bool_t passedVertexZSelectionMB = GetVertexIsOk(fEvent, kTRUE);
  // Check pile-up
  const Bool_t isPileUpMB = GetIsPileUp(fEvent);
  
  passedVertexSelectionMult = passedVertexZSelectionMult = kFALSE;
  isPileUpMult = kTRUE;
  passedDAQCheck = passedTrackClustCut = kFALSE;
  
  if (GetRunMode() == kJetPIDMode) {
    // Mult (check only needed for non-negative centrality percentile, otherwise, event is not used anyway
    // Check if is INEL > 0 (slight abuse of notation with "vertex selection"....)
    passedVertexSelectionMult = nonNegativeCentralityPercentileNoEventSelection ? AliPPVsMultUtils::IsINELgtZERO(fEvent) : kFALSE;
    // Check z position of vertex
    passedVertexZSelectionMult = nonNegativeCentralityPercentileNoEventSelection ? AliPPVsMultUtils::IsAcceptedVertexPosition(fEvent) : kFALSE;
    // Check pile-up  (and also consistency between SPD and track vertex, which is again a z cut, but is a check for pile-up!)
    isPileUpMult = nonNegativeCentralityPercentileNoEventSelection ? (!AliPPVsMultUtils::IsNotPileupSPDInMultBins(fEvent)  ||
                                                                                 !AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices(fEvent)) : kTRUE;
  }
  if (GetRunMode() == kLightFlavorMode) {
    // Check Incomplete DAQ rejection
    passedDAQCheck = !fEvent->IsIncompleteDAQ();
    //Check Tracklet vs. Cluster Cut
    passedTrackClustCut = !fAnaUtils->IsSPDClusterVsTrackletBG(fEvent);
  }
  
  if (fDoBinZeroStudy && fMC) {
    for (Int_t iPart = 0; iPart < fMC->GetNumberOfTracks(); iPart++) { 
      AliMCParticle *mcPart  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(iPart));
      
      if (!mcPart)
          continue;
      
      if (!fMC->IsPhysicalPrimary(iPart)) 
          continue;
      
      const Double_t etaGen = mcPart->Eta();
      const Double_t ptGen = mcPart->Pt();
      
      Double_t values[kBinZeroStudyNumAxes] = { 0. };
      values[kBinZeroStudyGenPt] = ptGen;
      values[kBinZeroStudyGenEta] = etaGen;
      
      // For multiplicity selection:
      if (nonNegativeCentralityPercentileNoEventSelection) {
        values[kBinZeroStudyCentrality] = centralityPercentileNoEventSelection;
        fChargedGenPrimariesTriggerSel->Fill(values);
        if (passedVertexSelectionMult) {
            fChargedGenPrimariesTriggerSelVtxCut->Fill(values);
          if (passedVertexZSelectionMult) {
              fChargedGenPrimariesTriggerSelVtxCutZ->Fill(values);
            if (!isPileUpMult && nonNegativeCentralityPercentile) {
              // If nonNegativeCentralityPercentile is kFALSE, but nonNegativeCentralityPercentileNoEventSelection was true,
              // then this should only be due to pile-up
              values[kBinZeroStudyCentrality] = centralityPercentile;
              fChargedGenPrimariesTriggerSelVtxCutZPileUpRej->Fill(values);
            }
          }
        }
      }
      
      // For MB selection:
      values[kBinZeroStudyCentrality] = -13;
      fChargedGenPrimariesTriggerSel->Fill(values);
      if (passedVertexSelectionMB) {
          fChargedGenPrimariesTriggerSelVtxCut->Fill(values);
        if (passedVertexZSelectionMB) {
            fChargedGenPrimariesTriggerSelVtxCutZ->Fill(values);
          if (!isPileUpMB) {
            fChargedGenPrimariesTriggerSelVtxCutZPileUpRej->Fill(values);
          }
        }
      }
    }
  }
  
  //Increment event counters (trigger selection, vertex cuts and pile-up rejection) for MB and mult:
  //MB
  //Flags that indicate whether the event passed all selections for MB and/or mult
  Bool_t isMBSelected = kFALSE;
  Bool_t isMultSelected = kFALSE;
  
  
  if(GetRunMode() == AliAnalysisTaskPID::kJetPIDMode) {
    IncrementEventCounter(centralityPercentile, kTriggerSel);
    if (passedVertexSelectionMB) {
      IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCut); 
      if (passedVertexZSelectionMB) {
        IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
        if (!isPileUpMB) {
  /*      ATTENTION: Is this the right place for the pile-up rejection? Important to have still the proper bin-0 correction,
          which is done solely with sel and selVtx, since the zvtx selection does ~not change the spectra. The question is whether the pile-up
          rejection changes the spectra. If not, then it is perfectly fine to put it here and keep the usual histo for the normalisation to 
          number of events. But if it does change the spectra, this must somehow be corrected for.
          NOTE: multiplicity >= 0 usually implies a properly reconstructed vertex. Hence, the bin-0 correction cannot be done in multiplicity 
          bins. Furthermore, there seems to be no MC simulation with pile-up rejection, so the bin-0 correction cannot be extracted with it. 
          Pile-up rejection has only a minor impact, so maybe there is no need to dig further.*/
          IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCutAndZvtxCut);
          isMBSelected = kTRUE;
        }
      }
    }
  
//         Mult (again, only centrality percentile >= 0 considered)
    if (nonNegativeCentralityPercentileNoEventSelection) {
      IncrementEventCounter(centralityPercentileNoEventSelection, kTriggerSel);
      if (passedVertexSelectionMult) {
        IncrementEventCounter(centralityPercentileNoEventSelection, kTriggerSelAndVtxCut);
        if (passedVertexZSelectionMult) {
          IncrementEventCounter(centralityPercentileNoEventSelection, kTriggerSelAndVtxCutAndZvtxCutNoPileUpRejection);
          if (!isPileUpMult && nonNegativeCentralityPercentile) {
            /*NOTE: Same comment as for MB
            If nonNegativeCentralityPercentile is kFALSE, but nonNegativeCentralityPercentileNoEventSelection was true,
            then this should only be due to pile-up*/
            IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCutAndZvtxCut);
            isMultSelected = kTRUE;
          }
        }
      }
    }
  }
  
  if (GetRunMode() == AliAnalysisTaskPID::kLightFlavorMode) {
    if (passedDAQCheck && passedTrackClustCut && !isPileUpMB) {
      IncrementEventCounter(centralityPercentile, kTriggerSel);
      if (passedVertexSelectionMB) {
        IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCut); 
        if (passedVertexZSelectionMB) {
          IncrementEventCounter(centralityPercentile, kTriggerSelAndVtxCutAndZvtxCut);
	  isMBSelected = kTRUE;
        }
      }
    }
  }
  
  
//   Done, if neither MB, nor mult requirements fulfilled.
  if (!isMBSelected && !isMultSelected)
    return;
  
  
  
  Double_t magField = fEvent->GetMagneticField();
  
  if (fMC) {
    if (fDoPID || fDoEfficiency) {
      for (Int_t iPart = 0; iPart < fMC->GetNumberOfTracks(); iPart++) { 
        AliMCParticle *mcPart  = dynamic_cast<AliMCParticle*>(fMC->GetTrack(iPart));
        
        if (!mcPart)
            continue;
        
        // Define clean MC sample with corresponding particle level track cuts:
        // - MC-track must be in desired eta range
        // - MC-track must be physical primary
        // - Species must be one of those in question (everything else goes to the overflow bin of mcID)
        
        // Geometrie should be the same as on the reconstructed level -> By definition analysis within this eta interval
        if (!IsInAcceptedEtaRange(TMath::Abs(mcPart->Eta()))) continue;
        
        Int_t mcID = PDGtoMCID(mcPart->PdgCode());
        
        // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        Double_t chargeMC = mcPart->Charge() / 3.;
        
        if (TMath::Abs(chargeMC) < 0.01)
          continue; // Reject neutral particles (only relevant, if mcID is not used)
        
        if (!fMC->IsPhysicalPrimary(iPart)) 
            continue;
        
        if (fDoPID) {
          Double_t valuesGenYield[kGenYieldNumAxes] = {  static_cast<Double_t>(mcID), mcPart->Pt(), centralityPercentile, -1, -1, -1, -1 };
          valuesGenYield[GetIndexOfChargeAxisGenYield()] = chargeMC;
          
          if (isMultSelected)
            fhMCgeneratedYieldsPrimaries->Fill(valuesGenYield);
          
          if (isMBSelected) {
            valuesGenYield[kGenYieldCentrality] = -13;
            fhMCgeneratedYieldsPrimaries->Fill(valuesGenYield);
          }
        }
        
        
        if (fDoEfficiency) {
          Double_t valueEff[kEffNumAxes] = {  static_cast<Double_t>(mcID), mcPart->Pt(), mcPart->Eta(), chargeMC, centralityPercentile,
                                            -1, -1, -1, -1, -1 };
          
          if (isMultSelected)            fContainerEff->Fill(valueEff, kStepGenWithGenCuts);
          
          if (isMBSelected) {
            valueEff[kEffCentrality] = -13;
            fContainerEff->Fill(valueEff, kStepGenWithGenCuts);
          }
        }
      }
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
    const Bool_t tuneOnDataTPC = fPIDResponse->IsTunedOnData() &&
                                 ((fPIDResponse->GetTunedOnDataMask() & AliPIDResponse::kDetTPC) == AliPIDResponse::kDetTPC);
    Double_t dEdxTPC = tuneOnDataTPC ? fPIDResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
    if (dEdxTPC <= 0)
      continue;
    
    if(fTrackFilter && !fTrackFilter->IsSelected(track))
      continue;
    
    if (GetUseTPCCutMIGeo()) {
      if (!TPCCutMIGeo(track, fEvent))
        continue;
    }
    else if (GetUseTPCnclCut()) {
      if (!TPCnclCut(track))
        continue;
    }
    
    if(fUsePhiCut) {
      if (!PhiPrimeCut(track, magField))
        continue; // reject track
    }
    
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
      
      if (fDoEfficiency) {
        // For efficiency: Reconstructed track has survived all cuts on the detector level (excluding acceptance)
        // and has an associated MC track which is a physical primary and was generated inside the acceptance
        if (fMC->IsPhysicalPrimary(TMath::Abs(label)) &&
            IsInAcceptedEtaRange(TMath::Abs(mcTrack->Eta()))) {
          
          // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
          Double_t value[kEffNumAxes] = {  static_cast<Double_t>(mcID), mcTrack->Pt(), mcTrack->Eta(), mcTrack->Charge() / 3.,
                                           centralityPercentile, -1, -1, -1, -1, -1 };
          if (isMultSelected)
            fContainerEff->Fill(value, kStepRecWithGenCuts);
          
          if (isMBSelected) {
            value[kEffCentrality] = -13;
            fContainerEff->Fill(value, kStepRecWithGenCuts);
          }
            
          Double_t valueMeas[kEffNumAxes] = {  static_cast<Double_t>(mcID), track->Pt(), track->Eta(),  static_cast<Double_t>(track->Charge()), 
                                               centralityPercentile, -1, -1, -1, -1, -1 };
          if (isMultSelected)
            fContainerEff->Fill(valueMeas, kStepRecWithGenCutsMeasuredObs);
          
          if (isMBSelected) {
            valueMeas[kEffCentrality] = -13;
            fContainerEff->Fill(valueMeas, kStepRecWithGenCutsMeasuredObs);
          }
        }
      }
    }
   
    // Only process tracks inside the desired eta window    
    if (!IsInAcceptedEtaRange(TMath::Abs(track->Eta()))) continue;
   
    if (fDoPID || fDoDeDxCheck || fDoPtResolution) 
      ProcessTrack(track, pdg, centralityPercentile, -1, isMBSelected, isMultSelected); // No jet information in this case -> Set jet pT to -1
    
    if (fDoPtResolution) {
      if (mcTrack && fMC->IsPhysicalPrimary(TMath::Abs(label))) {
        // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        Double_t valuePtRes[kPtResNumAxes] = { -1, mcTrack->Pt(), track->Pt(), mcTrack->Charge() / 3., centralityPercentile };
        
        if (isMultSelected)
          FillPtResolution(mcID, valuePtRes);
        
        if (isMBSelected) {
          valuePtRes[kPtResCentrality] = -13;
          FillPtResolution(mcID, valuePtRes);
        }
      }
    }
    
    if (fDoEfficiency) {
      if (mcTrack) {
        Double_t valueRecAllCuts[kEffNumAxes] = {  static_cast<Double_t>(mcID), track->Pt(), track->Eta(), static_cast<Double_t>(track->Charge()), 
                                                   centralityPercentile, -1, -1, -1, -1, -1 };
        Double_t weight = IsSecondaryWithStrangeMotherMC(fMC, TMath::Abs(label)) ? GetMCStrangenessFactorCMS(fMC, mcTrack) : 1.0;
        
        if (isMultSelected) {
          fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObs);
          fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsStrangenessScaled, weight);
        }
        
        if (isMBSelected) {
          valueRecAllCuts[kEffCentrality] = -13;
          fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObs);
          fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsStrangenessScaled, weight);
        }
        
        
        // AliMCParticle->Charge() calls TParticlePDG->Charge(), which returns the charge in units of e0 / 3
        Double_t valueGenAllCuts[kEffNumAxes] = {  static_cast<Double_t>(mcID), mcTrack->Pt(), mcTrack->Eta(), mcTrack->Charge() / 3., 
                                                   centralityPercentile, -1, -1, -1, -1, -1 };
        if (fMC->IsPhysicalPrimary(TMath::Abs(label))) {
          if (isMultSelected) {
            valueRecAllCuts[kEffCentrality] = centralityPercentile;
            fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsPrimaries);
            fContainerEff->Fill(valueGenAllCuts, kStepRecWithRecCutsPrimaries);
          }
          if (isMBSelected) {
            valueRecAllCuts[kEffCentrality] = -13;
            valueGenAllCuts[kEffCentrality] = -13;
            fContainerEff->Fill(valueRecAllCuts, kStepRecWithRecCutsMeasuredObsPrimaries);
            fContainerEff->Fill(valueGenAllCuts, kStepRecWithRecCutsPrimaries);
          }
        }
      }
    }
  } //track loop 
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: UserExec -> Processing done\n", (char*)__FILE__, __LINE__);
  
  PostOutputData();
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: UserExec -> Done\n", (char*)__FILE__, __LINE__);
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
  
  
  if ((TMath::Abs(fSystematicScalingSplinesBelowThreshold - 1.0) > fgkEpsilon) ||
      (TMath::Abs(fSystematicScalingSplinesAboveThreshold - 1.0) > fgkEpsilon)) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if ((TMath::Abs(fSystematicScalingEtaCorrectionLowMomenta - 1.0) > fgkEpsilon) ||
      (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - 1.0) > fgkEpsilon)) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if ((TMath::Abs(fSystematicScalingEtaSigmaParaBelowThreshold - 1.0) > fgkEpsilon) ||
      (TMath::Abs(fSystematicScalingEtaSigmaParaAboveThreshold - 1.0) > fgkEpsilon)) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
  
  if (TMath::Abs(fSystematicScalingMultCorrection - 1.0) > fgkEpsilon) {
    fDoAnySystematicStudiesOnTheExpectedSignal = kTRUE;
    return;
  }
}


//_____________________________________________________________________________
void AliAnalysisTaskPID::ConfigureTaskForCurrentEvent(AliVEvent* event)
{
  // Configure the task for the current event. In particular, this is needed if the run number changes
  
  if (!event) {
    AliError("Could not set up task: no event!");
    return;
  }
  
  Int_t run = event->GetRunNumber();
  
  if (run != fRun){
    // If systematics on eta is investigated, need to calculate the maxEtaVariationMap
    if ((TMath::Abs(fSystematicScalingEtaCorrectionLowMomenta - 1.0) > fgkEpsilon) ||
        (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - 1.0) > fgkEpsilon)) {
      if (!CalculateMaxEtaVariationMapFromPIDResponse())
        AliFatal("Systematics on eta correction requested, but failed to calculate max eta varation map!");
    }
  }
  
  fRun = run;
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
void AliAnalysisTaskPID::GetJetTrackObservables(Double_t trackPt, Double_t jetPt, Double_t& z, Double_t& xi, Bool_t storeXi)
{
  // Uses trackPt and jetPt to obtain z and xi.
  
  z = (jetPt > 0 && trackPt >= 0) ? (trackPt / jetPt) : -1;
  if (storeXi) {
    xi = (z > 0) ? TMath::Log(1. / z) : -1;
  }
  else {
    xi = -1;
  }
  
  if(trackPt > (1. - 1e-06) * jetPt && trackPt < (1. + 1e-06) * jetPt) { // case z=1 : move entry to last histo bin <1
    z  = 1. - 1e-06;
    xi = storeXi ? 1e-06 : -1;
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
  
  if(fDebug > 5)
    printf("File: %s, Line: %d: ConvolutedGaus: mean %e, sigma %e, lambda %e\n", (char*)__FILE__, __LINE__, mean, sigma, lambda);
  
  return lambda/sigma*TMath::Exp(-lambda/sigma*(xx[0]-mean)+lambda*lambda*0.5)*0.5*TMath::Erfc((-xx[0]+mean+sigma*lambda)/sigma*fgkOneOverSqrt2);
}


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::FastGaus(Double_t x, Double_t mean, Double_t sigma) const
{
  // Calculate an unnormalised gaussian function with mean and sigma.

  if (sigma < fgkEpsilon)
    return 1.e30;
  
  const Double_t arg = (x - mean) / sigma;
  return exp(-0.5 * arg * arg);
}


//_____________________________________________________________________________
inline Double_t AliAnalysisTaskPID::FastNormalisedGaus(Double_t x, Double_t mean, Double_t sigma) const
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
Int_t AliAnalysisTaskPID::FindFirstBinAboveIn3dSubset(const TH3* hist, Double_t threshold, Int_t yBin,
                                                      Int_t zBin) const
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
Int_t AliAnalysisTaskPID::FindLastBinAboveIn3dSubset(const TH3* hist, Double_t threshold, Int_t yBin,
                                                     Int_t zBin) const
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
Bool_t AliAnalysisTaskPID::GetParticleFraction(Double_t trackPt, Double_t jetPt, Double_t centralityPercentile,
                                               AliPID::EParticleType species,
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
    Int_t trackPtBin = xAxis->FindFixBin(trackPt);
    
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
Bool_t AliAnalysisTaskPID::GetParticleFractions(Double_t trackPt, Double_t jetPt, Double_t centralityPercentile,
                                                Double_t* prob, Int_t smearSpeciesByError,
                                                Int_t takeIntoAccountSpeciesSysError, Bool_t uniformSystematicError) const
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
const TH3D* AliAnalysisTaskPID::GetParticleFractionHisto(Int_t species, Bool_t sysError) const
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
  
  if (GetEventGenerator() == kPythia6Perugia2011) {
    // Values from mcplots.cern.ch for Pythia 6.425, tune 350 (= Perugia 2011) 
  
    if (absMotherPDG == 310 || absMotherPDG == 321) { // K0s / K+ / K-
      if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.873424;
      else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.854657;
      else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.800455;
      else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.738324;
      else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.687298;
      else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.650806;
      else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.629848;
      else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.619261;
      else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.610045;
      else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.601626;
      else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.605392;
      else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.596221;
      else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.607561;
      else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.604021;
      else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.600392;
      else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.603259;
      else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.619247;
      else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.614940;
      else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.632294;
      else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.633038;
      else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.646769;
      else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.700733;
      else if(6.00 <= motherGenPt && motherGenPt < 8.00) fac = 0.664591;
      else if(8.00 <= motherGenPt && motherGenPt < 10.00) fac = 0.683853;
    }

    if (absMotherPDG == 3122) { // Lambda
    //if (absMotherPDG == 3122 || absMotherPDG == 3112 || absMotherPDG == 3222) { // Lambda / Sigma- / Sigma+
      if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.871675;
      else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.892235;
      else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.705598;
      else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.630633;
      else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.552697;
      else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.505789;
      else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.461067;
      else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.433770;
      else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.422565;
      else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.398517;
      else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.393404;
      else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.394656;
      else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.390861;
      else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.380383;
      else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.396162;
      else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.388568;
      else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.429110;
      else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.427236;
      else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.437851;
      else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.470140;
      else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.509113;
      else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.616101;
      else if(6.00 <= motherGenPt && motherGenPt < 8.00) fac = 0.832494;
      else if(8.00 <= motherGenPt && motherGenPt < 10.00) fac = 0.997015;
    }

    if (absMotherPDG == 3312 || absMotherPDG == 3322) { // xi
      if (0.00 <= motherGenPt && motherGenPt < 0.20) fac = 0.946902;
      else if(0.20 <= motherGenPt && motherGenPt < 0.40) fac = 0.885799;
      else if(0.40 <= motherGenPt && motherGenPt < 0.60) fac = 0.712161;
      else if(0.60 <= motherGenPt && motherGenPt < 0.80) fac = 0.595333;
      else if(0.80 <= motherGenPt && motherGenPt < 1.00) fac = 0.531432;
      else if(1.00 <= motherGenPt && motherGenPt < 1.20) fac = 0.424845;
      else if(1.20 <= motherGenPt && motherGenPt < 1.40) fac = 0.378739;
      else if(1.40 <= motherGenPt && motherGenPt < 1.60) fac = 0.347243;
      else if(1.60 <= motherGenPt && motherGenPt < 1.80) fac = 0.366527;
      else if(1.80 <= motherGenPt && motherGenPt < 2.00) fac = 0.332981;
      else if(2.00 <= motherGenPt && motherGenPt < 2.20) fac = 0.314722;
      else if(2.20 <= motherGenPt && motherGenPt < 2.40) fac = 0.287927;
      else if(2.40 <= motherGenPt && motherGenPt < 2.60) fac = 0.285158;
      else if(2.60 <= motherGenPt && motherGenPt < 2.80) fac = 0.296892;
      else if(2.80 <= motherGenPt && motherGenPt < 3.00) fac = 0.291956;
      else if(3.00 <= motherGenPt && motherGenPt < 3.20) fac = 0.347254;
      else if(3.20 <= motherGenPt && motherGenPt < 3.40) fac = 0.348836;
      else if(3.40 <= motherGenPt && motherGenPt < 3.60) fac = 0.287240;
      else if(3.60 <= motherGenPt && motherGenPt < 3.80) fac = 0.325536;
      else if(3.80 <= motherGenPt && motherGenPt < 4.00) fac = 0.364368;
      else if(4.00 <= motherGenPt && motherGenPt < 5.00) fac = 0.405848;
      else if(5.00 <= motherGenPt && motherGenPt < 6.00) fac = 0.439721;
    }  
  }
  else if (GetEventGenerator() == kPythia6Perugia0) {
    // Values from mcplots.cern.ch for Pythia 6 (Perugia 0) 
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
    //if (absMotherPDG == 3122 || absMotherPDG == 3112 || absMotherPDG == 3222) { // Lambda / Sigma- / Sigma+
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
AliAnalysisTaskPID::TOFpidInfo AliAnalysisTaskPID::GetTOFType(const AliVTrack* track, Int_t tofMode) const
{
  // Get the (locally defined) particle type judged by TOF
  if (!fStoreTOFInfo) {
    return kNoTOFinfo;
  }
  
  if (!fPIDResponse) {
    Printf("ERROR: fPIDResponse not available -> Cannot determine TOF type!");
    return kNoTOFinfo;
  }
  
  /*TODO still needs some further thinking how to implement it....
  // TOF PID status (kTOFout, kTIME) automatically checked by return value of ComputeTPCProbability;
  // also, probability array will be set there (no need to initialise here)
  Double_t p[AliPID::kSPECIES];
  const AliPIDResponse::EDetPidStatus tofStatus = fPIDResponse->ComputeTPCProbability(track, AliPID::kSPECIES, p);
  if (tofStatus != AliPIDResponse::kDetPidOk)
    return kNoTOFinfo;
  
  // Do not consider muons
  p[AliPID::kMuon] = 0.;
  
  // Probabilities are not normalised yet
  Double_t sum = 0.;
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    sum += p[i];
  
  if (sum <= 0.)
    return kNoTOFinfo;
  
  for (Int_t i = 0; i < AliPID::kSPECIES; i++)
    p[i] /= sum;
  
  Double_t probThreshold = -999.;
  
  // If there is only one distribution, the threshold corresponds to...
  if (tofMode == 0) {
    probThreshold = ;
  }
  else if (tofMode == 1) { // default
    probThreshold = 0.9973; // a 3-sigma inclusion cut
  }
  else if (tofMode == 2) {
    inclusion = 3.;
    exclusion = 3.5;
  }
  else {
    Printf("ERROR: Bad TOF mode: %d!", tofMode);
    return kNoTOFinfo;
  }
  
  */
  
  ///* OLD: cut with n-sigmas (ATTENTION: Does not take into account properly the TOF tail!)
  // Check kTOFout, kTIME, mismatch
  const AliPIDResponse::EDetPidStatus tofStatus = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
  if (tofStatus != AliPIDResponse::kDetPidOk)
    return kNoTOFinfo;
  
  Double_t nsigma[kNumTOFspecies + 1] = { -999., -999., -999., -999. };
  nsigma[kTOFpion]   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
  nsigma[kTOFkaon]   = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon);
  nsigma[kTOFproton] = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
  
  Double_t inclusion = -999;
  Double_t exclusion = -999;
  
  if (tofMode == 0) {
    inclusion = 3.;
    exclusion = 2.5;
  }
  else if (tofMode == 1) { // default
    inclusion = 3.;
    exclusion = 3.;
  }
  else if (tofMode == 2) {
    inclusion = 3.;
    exclusion = 3.5;
  }
  else {
    Printf("ERROR: Bad TOF mode: %d!", tofMode);
    return kNoTOFinfo;
  }
  
  // Do not cut on nSigma electron because this would also remove pions in a large pT range.
  // The electron contamination of the pion sample can be treated by TPC dEdx fits afterwards.
  if (TMath::Abs(nsigma[kTOFpion]) < inclusion && TMath::Abs(nsigma[kTOFkaon]) > exclusion && TMath::Abs(nsigma[kTOFproton]) > exclusion)
    return kTOFpion;
  if (TMath::Abs(nsigma[kTOFpion]) > exclusion && TMath::Abs(nsigma[kTOFkaon]) < inclusion && TMath::Abs(nsigma[kTOFproton]) > exclusion)
    return kTOFkaon;
  if (TMath::Abs(nsigma[kTOFpion]) > exclusion && TMath::Abs(nsigma[kTOFkaon]) > exclusion && TMath::Abs(nsigma[kTOFproton]) < inclusion)
    return kTOFproton;
  //*/
  
  // There are no TOF electrons selected because the purity is rather bad, even for small momenta
  // (also a small mismatch probability significantly affects electrons because their fraction
  // is small). There is no need for TOF electrons anyway, since the dEdx distribution of kaons and 
  // protons in a given pT bin (also at the dEdx crossings) is very different from that of electrons.
  // This is due to the steeply falling dEdx of p and K with momentum, whereas the electron dEdx stays
  // rather constant.
  // As a result, the TPC fit yields a more accurate electron fraction than the TOF selection can do.
  
  return kNoTOFpid;
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
Bool_t AliAnalysisTaskPID::SetParticleFractionHisto(const TH3D* hist, Int_t species, Bool_t sysError)
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
Bool_t AliAnalysisTaskPID::SetParticleFractionHistosFromFile(const TString filePathName, Bool_t sysError)
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
Double_t AliAnalysisTaskPID::GetMaxEtaVariation(Double_t dEdxSplines)
{
  // Returns the maximum eta variation (i.e. deviation of eta correction factor from unity) for the
  // given (spline) dEdx
  
  if (dEdxSplines < 1. || !fhMaxEtaVariation) {
    Printf("Error GetMaxEtaVariation: No map or invalid dEdxSplines (%f)!", dEdxSplines);
    return 999.;
  } 
  
  Int_t bin = fhMaxEtaVariation->GetXaxis()->FindFixBin(1. / dEdxSplines);
  
  if (bin == 0) 
    bin = 1;
  if (bin > fhMaxEtaVariation->GetXaxis()->GetNbins())
    bin = fhMaxEtaVariation->GetXaxis()->GetNbins();
  
  return fhMaxEtaVariation->GetBinContent(bin);
}


//_____________________________________________________________________________
Int_t AliAnalysisTaskPID::GetRandomParticleTypeAccordingToParticleFractions(Double_t trackPt, Double_t jetPt,
                                                                            Double_t centralityPercentile, 
                                                                            Bool_t smearByError,
                                                                            Bool_t takeIntoAccountSysError) const
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
AliAnalysisTaskPID::ErrorCode AliAnalysisTaskPID::GenerateDetectorResponse(AliAnalysisTaskPID::ErrorCode errCode, 
                                                                           Double_t mean, Double_t sigma,
                                                                           Double_t* responses, Int_t nResponses, 
                                                                           Bool_t usePureGaus)
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
  printf("TPCCutMIGeo: %d\n", GetUseTPCCutMIGeo());
  if (GetUseTPCCutMIGeo()) {
    printf("GetCutGeo: %f\n", GetCutGeo());
    printf("GetCutNcr: %f\n", GetCutNcr());
    printf("GetCutNcl: %f\n", GetCutNcl());
  }
  printf("TPCnclCut: %d\n", GetUseTPCnclCut());
  if (GetUseTPCnclCut()) {
    printf("GetCutPureNcl: %d\n", GetCutPureNcl());
  }
  
  printf("\n");
  
  printf("Centrality estimator: %s\n", GetCentralityEstimator().Data());
  
  printf("\n");
  
  printf("Pile up rejection type: %d\n", (Int_t)GetPileUpRejectionType());
  
  printf("\n");
  
  printf("Use MC-ID for signal generation: %d\n", GetUseMCidForGeneration());
  printf("Use ITS: %d\n", GetUseITS());
  printf("Use TOF: %d\n", GetUseTOF());
  printf("Store TOF: %d\n", GetStoreTOFInfo());
  printf("Use priors: %d\n", GetUsePriors());
  printf("Use TPC default priors: %d\n", GetUseTPCDefaultPriors());
  printf("Store Charge: %d\n", GetStoreCharge());
  printf("Use convoluted Gauss: %d\n", GetUseConvolutedGaus());
  printf("Accuracy of non-Gaussian tail: %e\n", GetAccuracyNonGaussianTail());
  printf("Take into account muons: %d\n", GetTakeIntoAccountMuons());
  printf("TOF mode: %d\n", GetTOFmode());
  printf("\nParams for transition from gauss to asymmetric shape:\n");
  printf("[0]: %e\n", GetConvolutedGaussTransitionPar(0));
  printf("[1]: %e\n", GetConvolutedGaussTransitionPar(1));
  printf("[2]: %e\n", GetConvolutedGaussTransitionPar(2));
  
  printf("\n");
  
  std::cout << "Run Mode: ";
  if (GetRunMode() == AliAnalysisTaskPID::kJetPIDMode) {
    std::cout << "Jet PID Mode";
  } 
  else if (GetRunMode() == AliAnalysisTaskPID::kLightFlavorMode) {
    std::cout << "Light Flavor Mode";
  }
  else {
    std::cout << "No recognized Run mode";
  }
  std::cout << std::endl;
  
  std::cout << "Store Charge: " << std::endl;  
  if (GetStoreCharge())
    std::cout << "true";
  else
    std::cout << "false";
  
  std::cout << std::endl;
  printf("Do PID: %d\n", fDoPID);
  printf("Do Efficiency: %d\n", fDoEfficiency);
  printf("Do PtResolution: %d\n", fDoPtResolution);
  printf("Do dEdxCheck: %d\n", fDoDeDxCheck);
  printf("Do binZeroStudy: %d\n", fDoBinZeroStudy);
  if (fDoEfficiency) {
    if (GetEventGenerator() == kPythia6Perugia2011) {
      std::cout << "Use Strangeness weighting factors for Pythia 6 (Perugia 2011)" << std::endl;
    }
    else if (GetEventGenerator() == kPythia6Perugia0) {
      std::cout << "Use Strangeness weighting factors for Pythia 6 (Perugia 0)" << std::endl;
    }
    else {
      std::cout << "Using kUnknown Strangeness weighting factors" << std::endl;
    }
  }
  
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
  printf("SplinesThr:\t%f\n", GetSystematicScalingSplinesThreshold());
  printf("SplinesBelowThr:\t%f\n", GetSystematicScalingSplinesBelowThreshold());
  printf("SplinesAboveThr:\t%f\n", GetSystematicScalingSplinesAboveThreshold());
  printf("EtaCorrMomThr:\t%f\n", GetSystematicScalingEtaCorrectionMomentumThr());
  printf("EtaCorrLowP:\t%f\n", GetSystematicScalingEtaCorrectionLowMomenta());
  printf("EtaCorrHighP:\t%f\n", GetSystematicScalingEtaCorrectionHighMomenta());
  printf("SigmaParaThr:\t%f\n", GetSystematicScalingEtaSigmaParaThreshold());
  printf("SigmaParaBelowThr:\t%f\n", GetSystematicScalingEtaSigmaParaBelowThreshold());
  printf("SigmaParaAboveThr:\t%f\n", GetSystematicScalingEtaSigmaParaAboveThreshold());
  printf("MultCorr:\t%f\n", GetSystematicScalingMultCorrection());
  printf("TOF mode: %d\n", GetTOFmode());
  
  printf("\n\n");
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPID::ProcessTrack(const AliVTrack* track, Int_t particlePDGcode, Double_t centralityPercentile,
                                        Double_t jetPt, Bool_t isMBSelected, Bool_t isMultSelected, Bool_t storeXi,
                                        Double_t radialDistanceToJet, Double_t jT) 
{
  // Process the track (generate expected response, fill histos, etc.).
  // particlePDGcode == 0 means data. Otherwise, the corresponding MC ID will be assumed.
  
  //Printf("Debug: Task %s is starting to process track: dEdx %f, pTPC %f, eta %f, ncl %d\n", GetName(), track->GetTPCsignal(), track->GetTPCmomentum(),
  //       track->Eta(), track->GetTPCsignalN());
  
  if(fDebug > 1)
    printf("File: %s, Line: %d: ProcessTrack\n", (char*)__FILE__, __LINE__);
  
  if (!isMBSelected && !isMultSelected)
    return kTRUE; // Obviously, this event was not selected at all
  
  if (!fDoPID && !fDoDeDxCheck && !fDoPtResolution)
    return kFALSE;
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: ProcessTrack -> Processing started\n", (char*)__FILE__, __LINE__);
  
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
  const Double_t pTPC = track->GetTPCmomentum();
  Double_t pT = track->Pt();
  
  Double_t z = -1, xi = -1;
  GetJetTrackObservables(pT, jetPt, z, xi, storeXi);
  
  Double_t trackCharge = fStoreCharge ? track->Charge() : -2;
  
  // TPC signal
  const Bool_t tuneOnDataTPC = fPIDResponse->IsTunedOnData() &&
                               ((fPIDResponse->GetTunedOnDataMask() & AliPIDResponse::kDetTPC) == AliPIDResponse::kDetTPC);
  Double_t dEdxTPC = tuneOnDataTPC ? fPIDResponse->GetTPCsignalTunedOnData(track) : track->GetTPCsignal();
  
  if (dEdxTPC <= 0) {
    if (fDebug > 1)
      Printf("Skipping track with strange dEdx value: dEdx %f, pTPC %f, eta %f, ncl %d\n", track->GetTPCsignal(), pTPC,
             track->Eta(), track->GetTPCsignalN());
    return kFALSE;
  }
  
  // Completely remove tracks from light nuclei defined by el and pr <dEdx> as function of pT.
  // A very loose cut to be sure not to throw away electrons and/or protons
  Double_t nSigmaPr = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
  Double_t nSigmaEl = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
  
  if ((pTPC >= 0.3 && (nSigmaPr > 10 && nSigmaEl > 10)) ||
      (pTPC <  0.3 && (nSigmaPr > 15 && nSigmaEl > 15))) {
    if (fDebug > 1)
      Printf("Skipping track which seems to be a light nucleus: dEdx %f, pTPC %f, pT %f, eta %f, ncl %d, nSigmaPr %f, nSigmaEl %f\n",
             track->GetTPCsignal(), pTPC, pT, track->Eta(), track->GetTPCsignalN(), nSigmaPr, nSigmaEl);
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
    if ((TMath::Abs(fSystematicScalingSplinesBelowThreshold - 1.0) > fgkEpsilon) ||
        (TMath::Abs(fSystematicScalingSplinesAboveThreshold - 1.0) > fgkEpsilon)) {
       
      // Tune turn-on of correction for pions (not so relevant for the other species, since at very large momenta)
      Double_t bg = 0;
      Double_t scaleFactor = 1.;
      Double_t usedSystematicScalingSplines = 1.;
        
      bg = pTPC / AliPID::ParticleMass(AliPID::kElectron);
      scaleFactor = 0.5 * (1. + TMath::Erf((bg - fSystematicScalingSplinesThreshold) / 4.));
      usedSystematicScalingSplines = fSystematicScalingSplinesBelowThreshold * (1 - scaleFactor)
                                            + fSystematicScalingSplinesAboveThreshold * scaleFactor;
      dEdxEl *= usedSystematicScalingSplines;
      
      bg = pTPC / AliPID::ParticleMass(AliPID::kKaon);
      scaleFactor = 0.5 * (1. + TMath::Erf((bg - fSystematicScalingSplinesThreshold) / 4.));
      usedSystematicScalingSplines = fSystematicScalingSplinesBelowThreshold * (1 - scaleFactor)
                                            + fSystematicScalingSplinesAboveThreshold * scaleFactor;
      dEdxKa *= usedSystematicScalingSplines;
      
      bg = pTPC / AliPID::ParticleMass(AliPID::kPion);
      scaleFactor = 0.5 * (1. + TMath::Erf((bg - fSystematicScalingSplinesThreshold) / 4.));
      usedSystematicScalingSplines = fSystematicScalingSplinesBelowThreshold * (1 - scaleFactor)
                                            + fSystematicScalingSplinesAboveThreshold * scaleFactor;
      dEdxPi *= usedSystematicScalingSplines;
      
      if (fTakeIntoAccountMuons) {
        bg = pTPC / AliPID::ParticleMass(AliPID::kMuon);
        scaleFactor = 0.5 * (1. + TMath::Erf((bg - fSystematicScalingSplinesThreshold) / 4.));
        usedSystematicScalingSplines = fSystematicScalingSplinesBelowThreshold * (1 - scaleFactor)
                                              + fSystematicScalingSplinesAboveThreshold * scaleFactor;
        dEdxMu *= usedSystematicScalingSplines;
      }
      
      
      bg = pTPC / AliPID::ParticleMass(AliPID::kProton);
      scaleFactor = 0.5 * (1. + TMath::Erf((bg - fSystematicScalingSplinesThreshold) / 4.));
      usedSystematicScalingSplines = fSystematicScalingSplinesBelowThreshold * (1 - scaleFactor)
                                            + fSystematicScalingSplinesAboveThreshold * scaleFactor;
      dEdxPr *= usedSystematicScalingSplines;
    }
    
    // Get the eta correction factors for the (modified) expected dEdx
    Double_t etaCorrEl = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxEl) : 1.;
    Double_t etaCorrKa = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxKa) : 1.;
    Double_t etaCorrPi = fPIDResponse->UseTPCEtaCorrection() ? fPIDResponse->GetTPCResponse().GetEtaCorrectionFast(track, dEdxPi) : 1.;
    Double_t etaCorrMu = fTakeIntoAccountMuons && fPIDResponse->UseTPCEtaCorrection() ? 
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
      // An ERF will be used to get (as a function of p_TPC) from one correction factor to the other within roughly 0.2 GeV/c
      Double_t usedSystematicScalingEtaCorrection = fSystematicScalingEtaCorrectionHighMomenta;
      
      if (TMath::Abs(fSystematicScalingEtaCorrectionHighMomenta - fSystematicScalingEtaCorrectionLowMomenta) > fgkEpsilon) {
        const Double_t fractionHighMomentumScaleFactor = 0.5 * (1. + TMath::Erf((pTPC - fSystematicScalingEtaCorrectionMomentumThr) / 0.1));
        usedSystematicScalingEtaCorrection = fSystematicScalingEtaCorrectionLowMomenta * (1 - fractionHighMomentumScaleFactor)
                                             + fSystematicScalingEtaCorrectionHighMomenta * fractionHighMomentumScaleFactor;
      }
      
      Double_t maxEtaVariationEl = GetMaxEtaVariation(dEdxEl);
      etaCorrEl = etaCorrEl * (1.0 + (usedSystematicScalingEtaCorrection - 1.) * (etaCorrEl - 1.0) / maxEtaVariationEl);
      
      Double_t maxEtaVariationKa = GetMaxEtaVariation(dEdxKa);
      etaCorrKa = etaCorrKa * (1.0 + (usedSystematicScalingEtaCorrection - 1.) * (etaCorrKa - 1.0) / maxEtaVariationKa);
      
      Double_t maxEtaVariationPi = GetMaxEtaVariation(dEdxPi);
      etaCorrPi = etaCorrPi * (1.0 + (usedSystematicScalingEtaCorrection - 1.) * (etaCorrPi - 1.0) / maxEtaVariationPi);
      
      if (fTakeIntoAccountMuons) {
        Double_t maxEtaVariationMu = GetMaxEtaVariation(dEdxMu);
        etaCorrMu = etaCorrMu * (1.0 + (usedSystematicScalingEtaCorrection - 1.) * (etaCorrMu - 1.0) / maxEtaVariationMu);
      }
      else
        etaCorrMu = 1.0;
      
      Double_t maxEtaVariationPr = GetMaxEtaVariation(dEdxPr);
      etaCorrPr = etaCorrPr * (1.0 + (usedSystematicScalingEtaCorrection - 1.) * (etaCorrPr - 1.0) / maxEtaVariationPr);
      
      
      /*OLD
      etaCorrEl = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrEl - 1.0);
      etaCorrKa = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrKa - 1.0);
      etaCorrPi = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrPi - 1.0);
      etaCorrMu = fTakeIntoAccountMuons ? (1.0 + usedSystematicScalingEtaCorrection * (etaCorrMu - 1.0)) : 1.0;
      etaCorrPr = 1.0 + usedSystematicScalingEtaCorrection * (etaCorrPr - 1.0);
      */
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
    
    Bool_t doSigmaSystematics = (TMath::Abs(fSystematicScalingEtaSigmaParaBelowThreshold - 1.0) > fgkEpsilon) ||
                                (TMath::Abs(fSystematicScalingEtaSigmaParaAboveThreshold - 1.0) > fgkEpsilon);
    
    
    Double_t dEdxElExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,
                                                                               kTRUE, kFALSE);
    Double_t systematicScalingEtaSigmaParaEl = 1.;
    if (doSigmaSystematics) {
      Double_t scaleFactor = 0.5 * (1. + TMath::Erf((dEdxElExpected - fSystematicScalingEtaSigmaParaThreshold) / 25.));
      systematicScalingEtaSigmaParaEl = fSystematicScalingEtaSigmaParaBelowThreshold * (1 - scaleFactor)
                                        + fSystematicScalingEtaSigmaParaAboveThreshold * scaleFactor;
    }
    Double_t sigmaRelEl = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kElectron, AliTPCPIDResponse::kdEdxDefault,
                                                                          kTRUE, kFALSE)
                          / dEdxElExpected * systematicScalingEtaSigmaParaEl * multiplicityCorrSigmaEl;
    
    
    Double_t dEdxKaExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault, 
                                                                               kTRUE, kFALSE);
    Double_t systematicScalingEtaSigmaParaKa = 1.;
    if (doSigmaSystematics) {
      Double_t scaleFactor = 0.5 * (1. + TMath::Erf((dEdxKaExpected - fSystematicScalingEtaSigmaParaThreshold) / 25.));
      systematicScalingEtaSigmaParaKa = fSystematicScalingEtaSigmaParaBelowThreshold * (1 - scaleFactor)
                                        + fSystematicScalingEtaSigmaParaAboveThreshold * scaleFactor;
    }
    Double_t sigmaRelKa = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kKaon, AliTPCPIDResponse::kdEdxDefault,
                                                                          kTRUE, kFALSE)
                          / dEdxKaExpected * systematicScalingEtaSigmaParaKa * multiplicityCorrSigmaKa;
    
    
    Double_t dEdxPiExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault,
                                                                               kTRUE, kFALSE);
    Double_t systematicScalingEtaSigmaParaPi = 1.;
    if (doSigmaSystematics) {
      Double_t scaleFactor = 0.5 * (1. + TMath::Erf((dEdxPiExpected - fSystematicScalingEtaSigmaParaThreshold) / 25.));
      systematicScalingEtaSigmaParaPi = fSystematicScalingEtaSigmaParaBelowThreshold * (1 - scaleFactor)
                                        + fSystematicScalingEtaSigmaParaAboveThreshold * scaleFactor;
    }
    Double_t sigmaRelPi = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kPion, AliTPCPIDResponse::kdEdxDefault,
                                                                          kTRUE, kFALSE)
                          / dEdxPiExpected * systematicScalingEtaSigmaParaPi * multiplicityCorrSigmaPi;
    
    
    Double_t sigmaRelMu = 999.;
    if (fTakeIntoAccountMuons) {
      Double_t dEdxMuExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, 
                                                                                 kTRUE, kFALSE);
      Double_t systematicScalingEtaSigmaParaMu = 1.;
      if (doSigmaSystematics) {
        Double_t scaleFactor = 0.5 * (1. + TMath::Erf((dEdxMuExpected - fSystematicScalingEtaSigmaParaThreshold) / 25.));
        systematicScalingEtaSigmaParaMu = fSystematicScalingEtaSigmaParaBelowThreshold * (1 - scaleFactor)
                                          + fSystematicScalingEtaSigmaParaAboveThreshold * scaleFactor;
      }
      sigmaRelMu = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kMuon, AliTPCPIDResponse::kdEdxDefault, kTRUE, kFALSE)
                   / dEdxMuExpected * systematicScalingEtaSigmaParaMu * multiplicityCorrSigmaMu;
    }
    
    
    Double_t dEdxPrExpected = fPIDResponse->GetTPCResponse().GetExpectedSignal(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,
                                                                               kTRUE, kFALSE);
    Double_t systematicScalingEtaSigmaParaPr = 1.;
    if (doSigmaSystematics) {
      Double_t scaleFactor = 0.5 * (1. + TMath::Erf((dEdxPrExpected - fSystematicScalingEtaSigmaParaThreshold) / 25.));
      systematicScalingEtaSigmaParaPr = fSystematicScalingEtaSigmaParaBelowThreshold * (1 - scaleFactor)
                                        + fSystematicScalingEtaSigmaParaAboveThreshold * scaleFactor;
    }
    Double_t sigmaRelPr = fPIDResponse->GetTPCResponse().GetExpectedSigma(track, AliPID::kProton, AliTPCPIDResponse::kdEdxDefault,
                                                                          kTRUE, kFALSE)
                          / dEdxPrExpected * systematicScalingEtaSigmaParaPr * multiplicityCorrSigmaPr;
    
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
    // No systematic studies on expected signal - just take it as it comes from the TPCPIDResponse
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
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: ProcessTrack -> Compute probabilities\n", (char*)__FILE__, __LINE__);
  
  
  // Use probabilities to weigh the response generation for the different species.
  // Also determine the most probable particle type.
  Double_t prob[AliPID::kSPECIESC];
  for (Int_t i = 0; i < AliPID::kSPECIESC; i++)
    prob[i] = 0;
  
//   //For Testing
//   
  
  fPIDcombined->ComputeProbabilities(track, fPIDResponse, prob);
  
  
  //Testing
  
//   if (track->GetTPCsignal()<=50 && track->GetTPCmomentum()>=2.5) {
//     std::cout << "Centrality in PID-Task: " << fPIDResponse->GetCurrentCentrality() << std::endl;
//     std::cout << "Track: dEdx " << track->GetTPCsignal() << " pTPC " << track->GetTPCmomentum() << std::endl;
//     Double_t priors[AliPID::kSPECIESC];
//     memset(priors,0,AliPID::kSPECIESC*sizeof(Double_t));  
//     fPIDcombined->GetPriors(track,priors,fPIDResponse->GetCurrentCentrality(),kTRUE);
//     for (Int_t i=0;i<AliPID::kSPECIESC;i++) {
//       std::cout << "Probabilities: " << prob[i] << std::endl;
//     }
//     for (Int_t i=0;i<AliPID::kSPECIESC;i++) {
//       std::cout << "Priors: " << priors[i] << std::endl;
//     }  
//     std::cout << std::endl;
//   }
//   
  //End of Test
  
  // If muons are not to be taken into account, just set their probability to zero and normalise the remaining probabilities later
  if (!fTakeIntoAccountMuons)
    prob[AliPID::kMuon] = 0;
  
  // Priors might cause problems in case of mismatches, i.e. mismatch particles will be identified as pions.
  // Can be easily caught for low p_TPC by just setting pion (and muon) probs to zero, if dEdx is larger than
  // expected dEdx of electrons and kaons
  if (pTPC < 1. && dEdxTPC > dEdxEl && dEdxTPC > dEdxKa) {
    prob[AliPID::kMuon] = 0;
    prob[AliPID::kPion] = 0;
  }
  
  
  // Bug: One can set the number of species for PIDcombined, but PIDcombined will call PIDresponse, which writes without testing
  // the probs for kSPECIESC (including light nuclei) into the array.
  // In this case, when only kSPECIES are considered, the probabilities have to be rescaled!
  
  // Exceptions:
  // 1. If the most probable PID is a light nucleus and above expected dEdx + 5 sigma of el to be sure not to mix up things in the
  //    high-pT region (also contribution there completely negligable!)
  // 2. If dEdx is larger than the expected dEdx + 5 sigma of el and pr, it is most likely a light nucleus
  //    (if dEdx is more than 5 sigma away from expectation, flat TPC probs are set... and the light nuclei splines are not
  //     too accurate)
  // In these cases:
  // The highest expected dEdx is either for pr or for el. Assign 100% probability to the species with highest expected dEdx then.
  // In all other cases: Throw away light nuclei probs and rescale other probs to 100%.
  
  // Find most probable species for the ID 
  Int_t maxIndex = TMath::LocMax(AliPID::kSPECIESC, prob);

  if ((prob[maxIndex] > 0 && maxIndex >= AliPID::kSPECIES && dEdxTPC > dEdxEl + 5. * sigmaEl) ||
      (dEdxTPC > dEdxEl + 5. * sigmaEl && dEdxTPC > dEdxPr + 5. * sigmaPr)) {
    for (Int_t i = 0; i < AliPID::kSPECIESC; i++)
      prob[i] = 0;
    
    if (dEdxEl > dEdxPr) {
      prob[AliPID::kElectron] = 1.;
      maxIndex = AliPID::kElectron;
    }
    else {
      prob[AliPID::kProton] = 1.;
      maxIndex = AliPID::kProton;
    }
  }
  else {
    for (Int_t i = AliPID::kSPECIES; i < AliPID::kSPECIESC; i++)
      prob[i] = 0;
    
    Double_t probSum = 0.;
    for (Int_t i = 0; i < AliPID::kSPECIES; i++)
      probSum += prob[i];
    
    if (probSum > 0) {
      for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        prob[i] /= probSum;
    }
    
    // If all probabilities are equal (should only happen if no priors are used), set most probable PID to pion
    // because LocMax returns just the first maximal value (i.e. AliPID::kElectron)
    if (TMath::Abs(prob[maxIndex] - 1. / AliPID::kSPECIES) < 1e-5)
      maxIndex = AliPID::kPion;
    
    // In all other cases, the maximum remains untouched from the scaling (and is < AliPID::kSPECIES by construction)
  }
  
  if (fDoDeDxCheck) {
    // Generate the expected responses in DeltaPrime and translate to real dEdx. Only take responses for exactly that species
    // (i.e. with pre-PID)
    
    Int_t numGenEntries = fgkMaxNumGenEntries; // fgkMaxNumGenEntries = 500
    
    ErrorCode errCode = kNoErrors;
  
    if(fDebug > 2)
      printf("File: %s, Line: %d: ProcessTrack -> Generate Responses for dEdx check\n", (char*)__FILE__, __LINE__);
    
    // Electrons
    errCode = GenerateDetectorResponse(errCode, 1., sigmaEl / dEdxEl, fGenRespElDeltaPrimeEl, numGenEntries);

    // Kaons
    errCode = GenerateDetectorResponse(errCode, 1., sigmaKa / dEdxKa, fGenRespKaDeltaPrimeKa, numGenEntries);

    // Pions
    errCode = GenerateDetectorResponse(errCode, 1., sigmaPi / dEdxPi, fGenRespPiDeltaPrimePi, numGenEntries);

    // Protons
    errCode = GenerateDetectorResponse(errCode, 1., sigmaPr / dEdxPr, fGenRespPrDeltaPrimePr, numGenEntries);
    
    if (errCode == kNoErrors) {
      Double_t genEntry[kDeDxCheckNumAxes];
      genEntry[kDeDxCheckJetPt] = jetPt;
      genEntry[kDeDxCheckEtaAbs] = TMath::Abs(track->Eta());
      genEntry[kDeDxCheckP] = pTPC;
      
        
      for (Int_t n = 0; n < numGenEntries; n++)  {
        // If no MC info is available or shall not be used, use weighting with priors to generate entries for the different species
        Double_t rnd = fRandom->Rndm(); // Produce uniformly distributed floating point in ]0, 1]
        
        // Consider generated response as originating from...
        if (rnd <= prob[AliPID::kElectron]) {
          genEntry[kDeDxCheckPID] = 0; // ... an electron
          genEntry[kDeDxCheckDeDx] = fGenRespElDeltaPrimeEl[n] * dEdxEl;
        }
        else if (rnd <= prob[AliPID::kElectron] + prob[AliPID::kKaon]) {
          genEntry[kDeDxCheckPID] = 1;  // ... a kaon
          genEntry[kDeDxCheckDeDx] = fGenRespKaDeltaPrimeKa[n] * dEdxKa;
        }
        else if (rnd <= prob[AliPID::kElectron] + prob[AliPID::kKaon] + prob[AliPID::kMuon] + prob[AliPID::kPion]) {
          genEntry[kDeDxCheckPID] = 2;  // ... a pion (or a muon)
          genEntry[kDeDxCheckDeDx] = fGenRespPiDeltaPrimePi[n] * dEdxPi;
        }
        else {
          genEntry[kDeDxCheckPID] = 3;  // ... a proton
          genEntry[kDeDxCheckDeDx] = fGenRespPrDeltaPrimePr[n] * dEdxPr;
        }
     
        fDeDxCheck->Fill(genEntry);
      }
    }
    
    if(fDebug > 2)
      printf("File: %s, Line: %d: ProcessTrack -> Generate Responses for dEdx check done\n", (char*)__FILE__, __LINE__);
  }
  
  if (fDoPtResolution) {
    // Check shared clusters, which is done together with the pT resolution
    Double_t qaEntry[kQASharedClsNumAxes];
    qaEntry[kQASharedClsJetPt] = jetPt;
    qaEntry[kQASharedClsPt] = pT;
    qaEntry[kDeDxCheckP] = pTPC;
    qaEntry[kQASharedClsNumSharedCls] = track->GetTPCSharedMapPtr()->CountBits();
    
    Int_t iRowInd = -1;
    // iRowInd == -1 for "all rows w/o multiple counting"
    qaEntry[kQASharedClsPadRow] = iRowInd;
    fQASharedCls->Fill(qaEntry);

    // Fill hist for every pad row with shared cluster
    for (iRowInd = 0; iRowInd < 159; iRowInd++) {
      if (track->GetTPCSharedMapPtr()->TestBitNumber(iRowInd)) {
        qaEntry[kQASharedClsPadRow] = iRowInd;
        fQASharedCls->Fill(qaEntry);
      }
    }
  }
  
  if (!fDoPID)
    return kTRUE;
  
  if (!isMC) {
    // Clean up the most probable PID for low momenta (not an issue for the probabilities, since fractions small,
    // but to get clean (and nice looking) templates (if most probable PID is used instead of generated responses), it should be done).
    // Problem: If more than 5 sigma away (and since priors also included), most probable PID for dEdx around protons could be pions.
    // Idea: Only clean selection for K and p possible. So, just take for the most probable PID the probs from TPC only calculated
    // by hand without restriction to 5 sigma and without priors. This will not affect the template generation, but afterwards one
    // can replace the K and p templates with the most probable PID distributions, so all other species templates remain the same.
    // I.e. it doesn't hurt if the most probable PID is distributed incorrectly between el, pi, mu.
    Bool_t maxIndexSetManually = kFALSE;
    
    if (pTPC < 1.) {
      Double_t probManualTPC[AliPID::kSPECIES];
       for (Int_t i = 0; i < AliPID::kSPECIES; i++)
        probManualTPC[i] = 0;
      
      probManualTPC[AliPID::kElectron] = TMath::Exp(-0.5*(dEdxTPC-dEdxEl)*(dEdxTPC-dEdxEl)/(sigmaEl*sigmaEl))/sigmaEl;
      // Muons are not important here, just ignore and save processing time
      probManualTPC[AliPID::kMuon] = 0.;//TMath::Exp(-0.5*(dEdxTPC-dEdxMu)*(dEdxTPC-dEdxMu)/(sigmaMu*sigmaMu))/sigmaMu;
      probManualTPC[AliPID::kPion] = TMath::Exp(-0.5*(dEdxTPC-dEdxPi)*(dEdxTPC-dEdxPi)/(sigmaPi*sigmaPi))/sigmaPi;
      probManualTPC[AliPID::kKaon] = TMath::Exp(-0.5*(dEdxTPC-dEdxKa)*(dEdxTPC-dEdxKa)/(sigmaKa*sigmaKa))/sigmaKa;
      probManualTPC[AliPID::kProton] = TMath::Exp(-0.5*(dEdxTPC-dEdxPr)*(dEdxTPC-dEdxPr)/(sigmaPr*sigmaPr))/sigmaPr;
      
      const Int_t maxIndexManualTPC = TMath::LocMax(AliPID::kSPECIES, probManualTPC);
      // Should always be > 0, but in case the dEdx is far off such that the machine accuracy sets every prob to zero,
      // better take the "old" result
      if (probManualTPC[maxIndexManualTPC] > 0.)
        maxIndex = maxIndexManualTPC;
      
      maxIndexSetManually = kTRUE;
    }
    
    
    // Translate from AliPID numbering to numbering of this class
    if (prob[maxIndex] > 0 || maxIndexSetManually) {
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
  
  TOFpidInfo tofPIDinfo = GetTOFType(track, fTOFmode);
  
  Double_t entry[fStoreAdditionalJetInformation ? kDataNumAxes : kDataNumAxes - fgkNumJetAxes];
  entry[kDataMCID] = binMC;
  entry[kDataPt] = pT;
  
  if (fStoreAdditionalJetInformation) {
    entry[kDataJetPt] = jetPt;
    entry[kDataZ] = z;
    entry[kDataXi] = xi;
    entry[kDataDistance] = radialDistanceToJet;
    entry[kDataJt] = jT;
  }
  
  entry[GetIndexOfChargeAxisData()] = trackCharge;
  entry[GetIndexOfTOFpidInfoAxisData()] = tofPIDinfo;
  
  // Now consider both cases MB and mult and set centrality to negative value for MB
  for (Int_t k = 0; k < 2; k++) {
    if (k == 0) {
      if (!isMultSelected)
        continue;
      
      entry[kDataCentrality] = centralityPercentile;
    }
    else {
      if (!isMBSelected)
        continue;
      
      entry[kDataCentrality] = -13;
    }
    
    entry[kDataSelectSpecies] = 0;
    entry[kDataDeltaPrimeSpecies] = deltaPrimeElectron;
    fhPIDdataAll->Fill(entry);
    
    entry[kDataSelectSpecies] = 1;
    entry[kDataDeltaPrimeSpecies] = deltaPrimeKaon;
    fhPIDdataAll->Fill(entry);
    
    entry[kDataSelectSpecies] = 2;
    entry[kDataDeltaPrimeSpecies] = deltaPrimePion;
    fhPIDdataAll->Fill(entry);
    
    entry[kDataSelectSpecies] = 3;
    entry[kDataDeltaPrimeSpecies] = deltaPrimeProton;
    fhPIDdataAll->Fill(entry);
  }
  
  
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
    genEntry[kGenDistance] = radialDistanceToJet;
    genEntry[kGenJt] = jT;
  }
  
  genEntry[GetIndexOfChargeAxisGen()] = trackCharge;
  genEntry[GetIndexOfTOFpidInfoAxisGen()] = tofPIDinfo;
  
  // Generate numGenEntries "responses" with fluctuations around the expected values.
  // fgkMaxNumGenEntries = 500 turned out to give reasonable templates even for highest track pT in 15-20 GeV/c jet pT bin.
  Int_t numGenEntries = fgkMaxNumGenEntries; // fgkMaxNumGenEntries = 500
  
  /*OLD: Different number of responses depending on pT and jet pT for fgkMaxNumGenEntries = 1000
   * => Problem: If threshold to higher number of responses inside a bin (or after rebinning), then the template
   * is biased to the higher pT.
  // Generate numGenEntries "responses" with fluctuations around the expected values.
  // The higher the (transverse) momentum, the more "responses" will be generated in order not to run out of statistics too fast.
  Int_t numGenEntries = 10;
 
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
  */
  
  
  ErrorCode errCode = kNoErrors;
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: ProcessTrack -> Generate Responses\n", (char*)__FILE__, __LINE__);
  
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
  
  if (errCode != kNoErrors) {
    if (errCode == kWarning && fDebug > 1) {
      Printf("Warning: Questionable detector response for track, most likely due to very low number of PID clusters! Debug output (dEdx_expected, sigma_expected):");
    }
    else 
      Printf("Error: Failed to generate detector response for track - skipped! Debug output (dEdx_expected, sigma_expected):");
    
    if (fDebug > 1) {
      Printf("Pr: %e, %e", dEdxPr, sigmaPr);
      Printf("Pi: %e, %e", dEdxPi, sigmaPi);
      Printf("El: %e, %e", dEdxEl, sigmaEl);
      Printf("Mu: %e, %e", dEdxMu, sigmaMu);
      Printf("Ka: %e, %e", dEdxKa, sigmaKa);
      Printf("track: dEdx %f, pTPC %f, eta %f, ncl %d\n", track->GetTPCsignal(), track->GetTPCmomentum(), track->Eta(), 
            track->GetTPCsignalN());
    }
    
    if (errCode != kWarning) {
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
    
    // Now consider both cases MB and mult and set centrality to negative value for MB
    for (Int_t k = 0; k < 2; k++) {
      if (k == 0) {
        if (!isMultSelected)
          continue;
        
        genEntry[kGenCentrality] = centralityPercentile;
      }
      else {
        if (!isMBSelected)
          continue;
        
        genEntry[kGenCentrality] = -13;
      }
      
      
      // Electrons
      genEntry[kGenSelectSpecies] = 0;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimeEl[n];
      fhGenEl->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 1;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimeKa[n];
      fhGenEl->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 2;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimePi[n];
      fhGenEl->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 3;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespElDeltaPrimePr[n];
      fhGenEl->Fill(genEntry);
      
      // Kaons
      genEntry[kGenSelectSpecies] = 0;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimeEl[n];
      fhGenKa->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 1;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimeKa[n];
      fhGenKa->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 2;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimePi[n];
      fhGenKa->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 3;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespKaDeltaPrimePr[n];
      fhGenKa->Fill(genEntry);
      
      // Pions
      genEntry[kGenSelectSpecies] = 0;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimeEl[n];
      fhGenPi->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 1;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimeKa[n];
      fhGenPi->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 2;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimePi[n];
      fhGenPi->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 3;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPiDeltaPrimePr[n];
      fhGenPi->Fill(genEntry);
      
      if (fTakeIntoAccountMuons) {
        // Muons, if desired
        genEntry[kGenSelectSpecies] = 0;
        genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimeEl[n];
        fhGenMu->Fill(genEntry);
        
        genEntry[kGenSelectSpecies] = 1;
        genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimeKa[n];
        fhGenMu->Fill(genEntry);
        
        genEntry[kGenSelectSpecies] = 2;
        genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimePi[n];
        fhGenMu->Fill(genEntry);
        
        genEntry[kGenSelectSpecies] = 3;
        genEntry[kGenDeltaPrimeSpecies] = fGenRespMuDeltaPrimePr[n];
        fhGenMu->Fill(genEntry);
      }
      
      // Protons
      genEntry[kGenSelectSpecies] = 0;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimeEl[n];
      fhGenPr->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 1;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimeKa[n];
      fhGenPr->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 2;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimePi[n];
      fhGenPr->Fill(genEntry);
      
      genEntry[kGenSelectSpecies] = 3;
      genEntry[kGenDeltaPrimeSpecies] = fGenRespPrDeltaPrimePr[n];
      fhGenPr->Fill(genEntry);
    }
  }
  
  if(fDebug > 2)
    printf("File: %s, Line: %d: ProcessTrack -> Done\n", (char*)__FILE__, __LINE__);
  
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
    delete [] oldFuncParams;
    
    return kFALSE; 
  }
  
  if (fitResGauss->GetParams()[2] <= 0) {
    AliError("Not able to calculate parameters for the transition \"restricted gauss\" -> \"convoluted gauss\": Sigma of gauss fit <= 0!\n");
    fConvolutedGausDeltaPrime->SetParameters(oldFuncParams);
    fConvolutedGausDeltaPrime->SetNpx(oldFuncNpx);
    fConvolutedGausDeltaPrime->SetRange(oldFuncRangeLow, oldFuncRangeUp);
    
    delete hInput;
    delete [] oldFuncParams;
    
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
    delete [] oldFuncParams;
    
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
  delete [] oldFuncParams;

  return kTRUE;
}


//_____________________________________________________________________________
AliAnalysisTaskPID::ErrorCode AliAnalysisTaskPID::SetParamsForConvolutedGaus(Double_t gausMean, Double_t gausSigma) 
{
  // Set parameters for convoluted gauss using parameters for a pure gaussian.
  // If SetConvolutedGaussLambdaParameter has not been called before to initialise the translation parameters,
  // some default parameters will be used and an error will show up.
  
  if(fDebug > 1)
    printf("File: %s, Line: %d: SetParamsForConvolutedGaus: mean %e, sigma %e\n", (char*)__FILE__, __LINE__, gausMean, gausSigma);
  
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
  
  if(fDebug > 3)
    printf("File: %s, Line: %d: SetParamsForConvolutedGaus -> Parameters set to: %e, %e, %e (transition pars: %e, %e, %e, %e)\n",
           (char*)__FILE__, __LINE__, par[0], par[1], par[2], fConvolutedGaussTransitionPars[0], fConvolutedGaussTransitionPars[1], 
           fConvolutedGaussTransitionPars[2], fgkSigmaReferenceForTransitionPars);
  
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
    fConvolutedGausDeltaPrime->SetNpx(fDeltaPrimeAxis->FindFixBin(rangeEnd) - fDeltaPrimeAxis->FindFixBin(rangeStart) + 1);
    
    fConvolutedGausDeltaPrime->SetNpx(fDeltaPrimeAxis->FindFixBin(rangeEnd) - fDeltaPrimeAxis->FindFixBin(rangeStart) + 1);
    //fConvolutedGausDeltaPrime->SetNpx(fhPIDdataAll->GetAxis(kDataDeltaPrimeSpecies)->FindFixBin(rangeEnd)
    //                                  - fhPIDdataAll->GetAxis(kDataDeltaPrimeSpecies)->FindFixBin(rangeStart) + 1);
    //fConvolutedGausDeltaPrime->SetNpx((rangeEnd - rangeStart) / fDeltaPrimeBinWidth + 1);
  }
  
  if(fDebug > 3)
    printf("File: %s, Line: %d: SetParamsForConvolutedGaus -> range %f - %f, error code %d\n", (char*)__FILE__, __LINE__,
           rangeStart, rangeEnd, errCode);
  
  return errCode;
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpGenHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt,
                                      Double_t* binsJt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kGenPt, binsPt);
  hist->SetBinEdges(kGenDeltaPrimeSpecies, binsDeltaPrime);
  hist->SetBinEdges(kGenCentrality, binsCent);
  
  if (fStoreAdditionalJetInformation) {
    hist->SetBinEdges(kGenJetPt, binsJetPt);
    hist->SetBinEdges(kGenJt, binsJt);
  }
  
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
  
  hist->GetAxis(kGenPt)->SetTitle("p_{T} (GeV/c)");
  
  hist->GetAxis(kGenDeltaPrimeSpecies)->SetTitle("TPC #Delta'_{species} (arb. unit)");
  
  hist->GetAxis(kGenCentrality)->SetTitle(Form("Centrality Percentile (%s)", GetPPCentralityEstimator().Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kGenJetPt)->SetTitle("p_{T}^{jet} (GeV/c)");
    
    hist->GetAxis(kGenZ)->SetTitle("z = p_{T}^{track} / p_{T}^{jet}");
    
    hist->GetAxis(kGenXi)->SetTitle("#xi = ln(p_{T}^{jet} / p_{T}^{track})");
    
    hist->GetAxis(kGenDistance)->SetTitle("R");
  
    hist->GetAxis(kGenJt)->SetTitle("j_{T} (GeV/c)");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisGen())->SetTitle("Charge (e_{0})");
  
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetTitle("TOF PID Info");
  // Offset is (TMath::Abs(kNoTOFinfo) + 1), such that bin 1 (= first label) corresponds to kNoTOFinfo (< 0)
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetBinLabel(kNoTOFinfo + (TMath::Abs(kNoTOFinfo) + 1), "No TOF Info");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetBinLabel(kNoTOFpid + (TMath::Abs(kNoTOFinfo) + 1), "No TOF PID");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetBinLabel(kTOFpion + (TMath::Abs(kNoTOFinfo) + 1), "#pi");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetBinLabel(kTOFkaon + (TMath::Abs(kNoTOFinfo) + 1), "K");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisGen())->SetBinLabel(kTOFproton + (TMath::Abs(kNoTOFinfo) + 1), "p");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpGenYieldHist(THnSparse* hist, Double_t* binsPt, Double_t* binsCent, Double_t* binsJetPt, Double_t* binsJt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kGenYieldPt, binsPt);
  hist->SetBinEdges(kGenYieldCentrality, binsCent);
  
  if (fStoreAdditionalJetInformation) {
    hist->SetBinEdges(kGenYieldJetPt, binsJetPt);
    hist->SetBinEdges(kGenYieldJt, binsJt);
  }
  
  for (Int_t i = 0; i < 5; i++)
    hist->GetAxis(kGenYieldMCID)->SetBinLabel(i + 1, AliPID::ParticleLatexName(i));
  
  // Set axes titles
  hist->GetAxis(kGenYieldMCID)->SetTitle("MC PID");
  hist->GetAxis(kGenYieldPt)->SetTitle("p_{T}^{gen} (GeV/c)");
  hist->GetAxis(kGenYieldCentrality)->SetTitle(Form("Centrality Percentile (%s)", GetPPCentralityEstimator().Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kGenYieldJetPt)->SetTitle("p_{T}^{jet, gen} (GeV/c)");
    
    hist->GetAxis(kGenYieldZ)->SetTitle("z = p_{T}^{track} / p_{T}^{jet}");
    
    hist->GetAxis(kGenYieldXi)->SetTitle("#xi = ln(p_{T}^{jet} / p_{T}^{track})");
    
    hist->GetAxis(kGenYieldDistance)->SetTitle("R");
  
    hist->GetAxis(kGenYieldJt)->SetTitle("j_{T} (GeV/c)");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisGenYield())->SetTitle("Charge (e_{0})");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpHist(THnSparse* hist, Double_t* binsPt, Double_t* binsDeltaPrime, Double_t* binsCent, Double_t* binsJetPt,
                                   Double_t* binsJt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kDataPt, binsPt);
  hist->SetBinEdges(kDataDeltaPrimeSpecies, binsDeltaPrime);
  hist->SetBinEdges(kDataCentrality, binsCent);
  
  if (fStoreAdditionalJetInformation) {
    hist->SetBinEdges(kDataJetPt, binsJetPt);
    hist->SetBinEdges(kDataJt, binsJt);
  }
  
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
  
  hist->GetAxis(kDataPt)->SetTitle("p_{T} (GeV/c)");
    
  hist->GetAxis(kDataDeltaPrimeSpecies)->SetTitle("TPC #Delta'_{species} (arb. unit)");
  
  hist->GetAxis(kDataCentrality)->SetTitle(Form("Centrality Percentile (%s)", GetPPCentralityEstimator().Data()));
  
  if (fStoreAdditionalJetInformation) {
    hist->GetAxis(kDataJetPt)->SetTitle("p_{T}^{jet} (GeV/c)");
  
    hist->GetAxis(kDataZ)->SetTitle("z = p_{T}^{track} / p_{T}^{jet}");
  
    hist->GetAxis(kDataXi)->SetTitle("#xi = ln(p_{T}^{jet} / p_{T}^{track})");
    
    hist->GetAxis(kDataDistance)->SetTitle("R");
  
    hist->GetAxis(kDataJt)->SetTitle("j_{T} (GeV/c)");
  }
  
  hist->GetAxis(GetIndexOfChargeAxisData())->SetTitle("Charge (e_{0})");
  
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetTitle("TOF PID Info");
  // Offset is (TMath::Abs(kNoTOFinfo) + 1), such that bin 1 (= first label) corresponds to kNoTOFinfo (< 0)
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetBinLabel(kNoTOFinfo + (TMath::Abs(kNoTOFinfo) + 1), "No TOF Info");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetBinLabel(kNoTOFpid + (TMath::Abs(kNoTOFinfo) + 1), "No TOF PID");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetBinLabel(kTOFpion + (TMath::Abs(kNoTOFinfo) + 1), "#pi");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetBinLabel(kTOFkaon + (TMath::Abs(kNoTOFinfo) + 1), "K");
  hist->GetAxis(GetIndexOfTOFpidInfoAxisData())->SetBinLabel(kTOFproton + (TMath::Abs(kNoTOFinfo) + 1), "p");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpPtResHist(THnSparse* hist, Double_t* binsPt, Double_t* binsJetPt, Double_t* binsCent) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kPtResJetPt, binsJetPt);
  hist->SetBinEdges(kPtResGenPt, binsPt);
  hist->SetBinEdges(kPtResRecPt, binsPt);
  hist->SetBinEdges(kPtResCentrality, binsCent);
  
  // Set axes titles
  hist->GetAxis(kPtResJetPt)->SetTitle("p_{T}^{jet, rec} (GeV/c)");
  hist->GetAxis(kPtResGenPt)->SetTitle("p_{T}^{gen} (GeV/c)");
  hist->GetAxis(kPtResRecPt)->SetTitle("p_{T}^{rec} (GeV/c)");  
  
  hist->GetAxis(kPtResCharge)->SetTitle("Charge (e_{0})");
  hist->GetAxis(kPtResCentrality)->SetTitle(Form("Centrality Percentile (%s)", GetPPCentralityEstimator().Data()));
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpSharedClsHist(THnSparse* hist, Double_t* binsPt, Double_t* binsJetPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  
  hist->SetBinEdges(kQASharedClsJetPt, binsJetPt);
  hist->SetBinEdges(kQASharedClsPt, binsPt);
  
  // Set axes titles
  hist->GetAxis(kQASharedClsJetPt)->SetTitle("#it{p}_{T}^{jet} (GeV/#it{c})");
  hist->GetAxis(kQASharedClsPt)->SetTitle("#it{p}_{T} (GeV/#it{c})");
  hist->GetAxis(kQASharedClsNumSharedCls)->SetTitle("#it{N}_{shared}^{cls}");
  hist->GetAxis(kQASharedClsPadRow)->SetTitle("Pad row");
  
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpDeDxCheckHist(THnSparse* hist, const Double_t* binsPt, const Double_t* binsJetPt, const Double_t* binsEtaAbs) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  hist->SetBinEdges(kDeDxCheckP, binsPt);
  hist->SetBinEdges(kDeDxCheckJetPt, binsJetPt);
  hist->SetBinEdges(kDeDxCheckEtaAbs, binsEtaAbs);
  
  
  // Set axes titles
  hist->GetAxis(kDeDxCheckPID)->SetTitle("PID");
  hist->GetAxis(kDeDxCheckPID)->SetBinLabel(1, "e");
  hist->GetAxis(kDeDxCheckPID)->SetBinLabel(2, "K");
  hist->GetAxis(kDeDxCheckPID)->SetBinLabel(3, "#pi");
  hist->GetAxis(kDeDxCheckPID)->SetBinLabel(4, "p");
  
  hist->GetAxis(kDeDxCheckJetPt)->SetTitle("p_{T}^{jet} (GeV/c)");
  hist->GetAxis(kDeDxCheckEtaAbs)->SetTitle("|#eta|");  
  hist->GetAxis(kDeDxCheckP)->SetTitle("p_{TPC} (GeV/c)"); 
  hist->GetAxis(kDeDxCheckDeDx)->SetTitle("TPC dE/dx (arb. unit)");
}


//________________________________________________________________________
void AliAnalysisTaskPID::SetUpBinZeroStudyHist(THnSparse* hist, const Double_t* binsCent, const Double_t* binsPt) const
{
  // Sets bin limits for axes which are not standard binned and the axes titles.
  hist->SetBinEdges(kBinZeroStudyCentrality, binsCent);
  hist->SetBinEdges(kBinZeroStudyGenPt, binsPt);
  
  // Set axes titles
  hist->GetAxis(kBinZeroStudyCentrality)->SetTitle(Form("Centrality Percentile (%s)", GetPPCentralityEstimator().Data()));
  hist->GetAxis(kBinZeroStudyGenEta)->SetTitle("#it{#eta}^{gen}");  
  hist->GetAxis(kBinZeroStudyGenPt)->SetTitle("#it{p}_{T}^{gen} (GeV/#it{c})"); 
}


//________________________________________________________________________
void AliAnalysisTaskPID::FillUEDensity(Double_t cent, Double_t UEpt) {
  if (fIsUEPID) {
    fh2UEDensity->Fill(cent, UEpt);
  }
}

void AliAnalysisTaskPID::FillJetArea(Double_t cent, Double_t area) {
  if (fIsUEPID) {
    fh1JetArea->Fill(cent,area);
  }
}

void AliAnalysisTaskPID::NormalizeJetArea(Double_t jetParameter) {
  if (!fIsUEPID)
    return;
  
  if (!fh2FFJetPtRec || !fh1JetArea) {
    std::cout << "Cannot normalize Jet Area" << std::endl;
    return;
  }
  
  fh1JetArea->Scale(1.0/(fh2FFJetPtRec->GetEntries() * jetParameter * jetParameter * TMath::Pi()));
}
