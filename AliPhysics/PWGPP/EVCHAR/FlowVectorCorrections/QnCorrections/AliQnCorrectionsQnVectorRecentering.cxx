/**************************************************************************************************
 *                                                                                                *
 * Package:       FlowVectorCorrections                                                           *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch                              *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com                             *
 *                Víctor González, UCM, victor.gonzalez@cern.ch                                   *
 *                Contributors are mentioned in the code where appropriate.                       *
 * Development:   2012-2016                                                                       *
 *                                                                                                *
 * This file is part of FlowVectorCorrections, a software package that corrects Q-vector          *
 * measurements for effects of nonuniform detector acceptance. The corrections in this package    *
 * are based on publication:                                                                      *
 *                                                                                                *
 *  [1] "Effects of non-uniform acceptance in anisotropic flow measurements"                      *
 *  Ilya Selyuzhenkov and Sergei Voloshin                                                         *
 *  Phys. Rev. C 77, 034904 (2008)                                                                *
 *                                                                                                *
 * The procedure proposed in [1] is extended with the following steps:                            *
 * (*) alignment correction between subevents                                                     *
 * (*) possibility to extract the twist and rescaling corrections                                 *
 *      for the case of three detector subevents                                                  *
 *      (currently limited to the case of two “hit-only” and one “tracking” detectors)            *
 * (*) (optional) channel equalization                                                            *
 * (*) flow vector width equalization                                                             *
 *                                                                                                *
 * FlowVectorCorrections is distributed under the terms of the GNU General Public License (GPL)   *
 * (https://en.wikipedia.org/wiki/GNU_General_Public_License)                                     *
 * either version 3 of the License, or (at your option) any later version.                        *
 *                                                                                                *
 **************************************************************************************************/

/// \file AliQnCorrectionsQnVectorRecentering.cxx
/// \brief Implementation of procedures for Qn vector recentering.
#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileComponents.h"
#include "AliQnCorrectionsHistogramSparse.h"
#include "AliQnCorrectionsDetector.h"
#include "AliLog.h"
#include "AliQnCorrectionsQnVectorRecentering.h"

const Int_t AliQnCorrectionsQnVectorRecentering::fDefaultMinNoOfEntries = 2;
const char *AliQnCorrectionsQnVectorRecentering::szCorrectionName = "Recentering and width equalization";
const char *AliQnCorrectionsQnVectorRecentering::szKey = "CCCC";
const char *AliQnCorrectionsQnVectorRecentering::szSupportHistogramName = "Qn";
const char *AliQnCorrectionsQnVectorRecentering::szCorrectedQnVectorName = "rec";
const char *AliQnCorrectionsQnVectorRecentering::szQANotValidatedHistogramName = "Rec NvE";
const char *AliQnCorrectionsQnVectorRecentering::szQAQnAverageHistogramName = "Rec Qn avg ";

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsQnVectorRecentering);
/// \endcond

/// Default constructor
/// Passes to the base class the identity data for the recentering and width equalization correction step
AliQnCorrectionsQnVectorRecentering::AliQnCorrectionsQnVectorRecentering() :
    AliQnCorrectionsCorrectionOnQvector(szCorrectionName, szKey) {
  fInputHistograms = NULL;
  fCalibrationHistograms = NULL;
  fQANotValidatedBin = NULL;
  fQAQnAverageHistogram = NULL;
  fApplyWidthEqualization = kFALSE;
  fMinNoOfEntriesToValidate = fDefaultMinNoOfEntries;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsQnVectorRecentering::~AliQnCorrectionsQnVectorRecentering() {
  if (fInputHistograms != NULL)
    delete fInputHistograms;
  if (fCalibrationHistograms != NULL)
    delete fCalibrationHistograms;
  if (fQANotValidatedBin != NULL)
    delete fQANotValidatedBin;
  if (fQAQnAverageHistogram != NULL)
    delete fQAQnAverageHistogram;
}

/// Asks for support data structures creation
///
/// Creates the recentered Qn vector
void AliQnCorrectionsQnVectorRecentering::CreateSupportDataStructures() {

  Int_t nNoOfHarmonics = fDetectorConfiguration->GetNoOfHarmonics();
  Int_t *harmonicsMap = new Int_t[nNoOfHarmonics];
  fDetectorConfiguration->GetHarmonicMap(harmonicsMap);
  fCorrectedQnVector = new AliQnCorrectionsQnVector(szCorrectedQnVectorName, nNoOfHarmonics, harmonicsMap);
  fInputQnVector = fDetectorConfiguration->GetPreviousCorrectedQnVector(this);
  delete [] harmonicsMap;
}

/// Asks for support histograms creation
///
/// Allocates the histogram objects and creates the calibration histograms.
/// The histograms are constructed with standard deviation error calculation
/// for the proper behavior of optional gain equalization step.
///
/// Process concurrency requires Calibration Histograms creation for all c
/// concurrent processes but not for Input Histograms so, we delete previously
/// allocated ones.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorRecentering::CreateSupportHistograms(TList *list) {

  TString histoNameAndTitle = TString::Format("%s %s ",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());

  if (fInputHistograms != NULL) delete fInputHistograms;
  fInputHistograms = new AliQnCorrectionsProfileComponents((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      fDetectorConfiguration->GetEventClassVariablesSet(), "s");
  fInputHistograms->SetNoOfEntriesThreshold(fMinNoOfEntriesToValidate);
  fCalibrationHistograms = new AliQnCorrectionsProfileComponents((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      fDetectorConfiguration->GetEventClassVariablesSet(), "s");

  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = fDetectorConfiguration->GetNoOfHarmonics();
  Int_t *harmonicsMap = new Int_t[nNoOfHarmonics];
  fDetectorConfiguration->GetHarmonicMap(harmonicsMap);
  fCalibrationHistograms->CreateComponentsProfileHistograms(list,nNoOfHarmonics, harmonicsMap);
  delete [] harmonicsMap;
  return kTRUE;
}

/// Attaches the needed input information to the correction step
/// \param list list where the inputs should be found
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorRecentering::AttachInput(TList *list) {

  if (fInputHistograms->AttachHistograms(list)) {
    AliInfo(TString::Format("Recentering on %s going to be applied", fDetectorConfiguration->GetName()).Data());
    fState = QCORRSTEP_applyCollect;
    return kTRUE;
  }
  return kFALSE;
}

/// Asks for QA histograms creation
///
/// Allocates the histogram objects and creates the QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorRecentering::CreateQAHistograms(TList *list) {

  fQAQnAverageHistogram = new AliQnCorrectionsProfileComponents(
      TString::Format("%s %s", szQAQnAverageHistogramName, fDetectorConfiguration->GetName()).Data(),
      TString::Format("%s %s", szQAQnAverageHistogramName, fDetectorConfiguration->GetName()).Data(),
      fDetectorConfiguration->GetEventClassVariablesSet());

  /* get information about the configured harmonics to pass it for histogram creation */
  Int_t nNoOfHarmonics = fDetectorConfiguration->GetNoOfHarmonics();
  Int_t *harmonicsMap = new Int_t[nNoOfHarmonics];
  fDetectorConfiguration->GetHarmonicMap(harmonicsMap);
  fQAQnAverageHistogram->CreateComponentsProfileHistograms(list,nNoOfHarmonics, harmonicsMap);
  delete [] harmonicsMap;
  return kTRUE;
}

/// Asks for non validated entries QA histograms creation
///
/// Allocates the histogram objects and creates the non validated entries QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorRecentering::CreateNveQAHistograms(TList *list) {

  fQANotValidatedBin = new AliQnCorrectionsHistogramSparse(
      TString::Format("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()).Data(),
      TString::Format("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()).Data(),
      fDetectorConfiguration->GetEventClassVariablesSet());
  fQANotValidatedBin->CreateHistogram(list);
  return kTRUE;
}

/// Processes the correction step
///
/// Pure virtual function
/// \return kTRUE if the correction step was applied
Bool_t AliQnCorrectionsQnVectorRecentering::ProcessCorrections(const Float_t *variableContainer) {
  Int_t harmonic;
  switch (fState) {
  case QCORRSTEP_calibration:
    /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
    /* we have not perform any correction yet */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the correction if the current Qn vector is good enough */
    AliInfo(TString::Format("Recentering process in detector %s: applying correction.", fDetectorConfiguration->GetName()).Data());
    if (fDetectorConfiguration->GetCurrentQnVector()->IsGoodQuality()) {
      /* we get the properties of the current Qn vector but its name */
      fCorrectedQnVector->Set(fDetectorConfiguration->GetCurrentQnVector(),kFALSE);
      harmonic = fDetectorConfiguration->GetCurrentQnVector()->GetFirstHarmonic();

      /* let's check the correction histograms */
      Long64_t bin = fInputHistograms->GetBin(variableContainer);
      if (fInputHistograms->BinContentValidated(bin)) {
        /* correction information validated */
        while (harmonic != -1) {
          Float_t widthX = 1.0;
          Float_t widthY = 1.0;
          if (fApplyWidthEqualization) {
            widthX = fInputHistograms->GetXBinError(harmonic, bin);
            widthY = fInputHistograms->GetYBinError(harmonic, bin);
          }
          fCorrectedQnVector->SetQx(harmonic, (fDetectorConfiguration->GetCurrentQnVector()->Qx(harmonic)
              - fInputHistograms->GetXBinContent(harmonic, bin))
              / widthX);
          fCorrectedQnVector->SetQy(harmonic, (fDetectorConfiguration->GetCurrentQnVector()->Qy(harmonic)
              - fInputHistograms->GetYBinContent(harmonic, bin))
              / widthY);
          harmonic = fDetectorConfiguration->GetCurrentQnVector()->GetNextHarmonic(harmonic);
        }
      } /* correction information not validated, we leave the Q vector untouched */
      else {
        if (fQANotValidatedBin != NULL) fQANotValidatedBin->Fill(variableContainer, 1.0);
      }
    }
    else {
      /* not done! input vector with bad quality */
      fCorrectedQnVector->SetGood(kFALSE);
    }
    /* and update the current Qn vector */
    fDetectorConfiguration->UpdateCurrentQnVector(fCorrectedQnVector);
    break;
  default:
    /* we are in passive state waiting for proper conditions, no corrections applied */
    return kFALSE;
  }
  /* if we reached here is because we applied the correction */
  return kTRUE;
}

/// Processes the correction step data collection
///
/// Pure virtual function
/// \return kTRUE if the correction step was applied
Bool_t AliQnCorrectionsQnVectorRecentering::ProcessDataCollection(const Float_t *variableContainer) {
  Int_t harmonic;
  switch (fState) {
  case QCORRSTEP_calibration:
    AliInfo(TString::Format("Recentering process in detector %s: collecting data.", fDetectorConfiguration->GetName()).Data());
    /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
    if (fInputQnVector->IsGoodQuality()) {
      harmonic = fInputQnVector->GetFirstHarmonic();
      while (harmonic != -1) {
        fCalibrationHistograms->FillX(harmonic,variableContainer,fInputQnVector->Qx(harmonic));
        fCalibrationHistograms->FillY(harmonic,variableContainer,fInputQnVector->Qy(harmonic));
        harmonic = fInputQnVector->GetNextHarmonic(harmonic);
      }
    }
    /* we have not perform any correction yet */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    AliInfo(TString::Format("Recentering process in detector %s: collecting data.", fDetectorConfiguration->GetName()).Data());
    /* collect the data needed to further produce correction parameters if the current Qn vector is good enough */
    if (fInputQnVector->IsGoodQuality()) {
      harmonic = fInputQnVector->GetFirstHarmonic();
      while (harmonic != -1) {
        fCalibrationHistograms->FillX(harmonic,variableContainer,fInputQnVector->Qx(harmonic));
        fCalibrationHistograms->FillY(harmonic,variableContainer,fInputQnVector->Qy(harmonic));
        harmonic = fInputQnVector->GetNextHarmonic(harmonic);
      }
    }
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the correction if the current Qn vector is good enough */
    /* provide QA info if required */
    if (fQAQnAverageHistogram != NULL) {
      harmonic = fCorrectedQnVector->GetFirstHarmonic();
      while (harmonic != -1) {
        fQAQnAverageHistogram->FillX(harmonic, variableContainer, fCorrectedQnVector->Qx(harmonic));
        fQAQnAverageHistogram->FillY(harmonic, variableContainer, fCorrectedQnVector->Qy(harmonic));
        harmonic = fCorrectedQnVector->GetNextHarmonic(harmonic);
      }
    }
    break;
  default:
    /* we are in passive state waiting for proper conditions, no corrections applied */
    return kFALSE;
  }
  /* if we reached here is because we applied the correction */
  return kTRUE;
}

/// Clean the correction to accept a new event
void AliQnCorrectionsQnVectorRecentering::ClearCorrectionStep() {

  fCorrectedQnVector->Reset();
}

/// Reports if the correction step is being applied
/// Returns TRUE if in the proper state for applying the correction step
/// \return TRUE if the correction step is being applied
Bool_t AliQnCorrectionsQnVectorRecentering::IsBeingApplied() const {
  switch (fState) {
  case QCORRSTEP_calibration:
    /* we are collecting */
    /* but not applying */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* we are collecting */
  case QCORRSTEP_apply:
    /* and applying */
    return kTRUE;
    break;
  default:
    break;
  }
  return kFALSE;
}

/// Report on correction usage
/// Correction step should incorporate its name in calibration
/// list if it is producing information calibration in the ongoing
/// step and in the apply list if it is applying correction in
/// the ongoing step.
/// \param calibrationList list containing the correction steps producing calibration information
/// \param applyList list containing the correction steps applying corrections
/// \return kTRUE if the correction step is being applied
Bool_t AliQnCorrectionsQnVectorRecentering::ReportUsage(TList *calibrationList, TList *applyList) {
  switch (fState) {
  case QCORRSTEP_calibration:
    /* we are collecting */
    calibrationList->Add(new TObjString(szCorrectionName));
    /* but not applying */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* we are collecting */
    calibrationList->Add(new TObjString(szCorrectionName));
  case QCORRSTEP_apply:
    /* and applying */
    applyList->Add(new TObjString(szCorrectionName));
    break;
  default:
    return kFALSE;
    break;
  }
  return kTRUE;
}

