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

/// \file AliQnCorrectionsQnVectorAlignment.cxx
/// \brief Implementation of procedures for Qn vector alignment correction.
#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileCorrelationComponents.h"
#include "AliQnCorrectionsProfileComponents.h"
#include "AliQnCorrectionsHistogramSparse.h"
#include "AliQnCorrectionsDetector.h"
#include "AliQnCorrectionsManager.h"
#include "AliLog.h"
#include "AliQnCorrectionsQnVectorAlignment.h"

const Int_t AliQnCorrectionsQnVectorAlignment::fDefaultMinNoOfEntries = 2;
const char *AliQnCorrectionsQnVectorAlignment::szCorrectionName = "Alignment";
const char *AliQnCorrectionsQnVectorAlignment::szKey = "EEEE";
const char *AliQnCorrectionsQnVectorAlignment::szSupportHistogramName = "QnQn";
const char *AliQnCorrectionsQnVectorAlignment::szCorrectedQnVectorName = "align";
const char *AliQnCorrectionsQnVectorAlignment::szQANotValidatedHistogramName = "Align NvE";
const char *AliQnCorrectionsQnVectorAlignment::szQAQnAverageHistogramName = "Align Qn avg ";


/// \cond CLASSIMP
ClassImp(AliQnCorrectionsQnVectorAlignment);
/// \endcond

/// Default constructor
/// Passes to the base class the identity data for the recentering and width equalization correction step
AliQnCorrectionsQnVectorAlignment::AliQnCorrectionsQnVectorAlignment() :
    AliQnCorrectionsCorrectionOnQvector(szCorrectionName, szKey),
    fDetectorConfigurationForAlignmentName() {
  fInputHistograms = NULL;
  fCalibrationHistograms = NULL;
  fQANotValidatedBin = NULL;
  fQAQnAverageHistogram = NULL;
  fHarmonicForAlignment = -1;
  fDetectorConfigurationForAlignment = NULL;
  fMinNoOfEntriesToValidate = fDefaultMinNoOfEntries;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsQnVectorAlignment::~AliQnCorrectionsQnVectorAlignment() {
  if (fInputHistograms != NULL)
    delete fInputHistograms;
  if (fCalibrationHistograms != NULL)
    delete fCalibrationHistograms;
  if (fQANotValidatedBin != NULL)
    delete fQANotValidatedBin;
  if (fQAQnAverageHistogram != NULL)
    delete fQAQnAverageHistogram;
}

/// Set the detector configuration used as reference for alignment
/// The detector configuration name is stored for further use.
/// \param name the name of the reference detector configuration
void AliQnCorrectionsQnVectorAlignment::SetReferenceConfigurationForAlignment(const char *name) {
  AliInfo(TString::Format("Reference name: %s, attached to detector configuration: %s",
      name,
      ((fDetectorConfiguration != NULL) ? "yes" : "no")).Data());

  fDetectorConfigurationForAlignmentName = name;

  /* we could be in different situations of framework attachment */
  /* so, we do nothing for the time being */
}

/// Informs when the detector configuration has been attached to the framework manager
/// Basically this allows interaction between the different framework sections at configuration time
void AliQnCorrectionsQnVectorAlignment::AttachedToFrameworkManager() {
  AliInfo(TString::Format("Attached! reference for alignment: %s", fDetectorConfigurationForAlignmentName.Data()).Data());

}

/// Asks for support data structures creation
///
/// Locates the reference detector configuration for alignment if its name has been previously stored
/// Creates the recentered Qn vector
void AliQnCorrectionsQnVectorAlignment::CreateSupportDataStructures() {

  /* now, definitely, we should have the reference detector configurations */
  if (fDetectorConfigurationForAlignmentName.Length() != 0) {
    if (fDetectorConfiguration->GetCorrectionsManager()->FindDetectorConfiguration(fDetectorConfigurationForAlignmentName.Data()) != NULL) {
      fDetectorConfigurationForAlignment = fDetectorConfiguration->GetCorrectionsManager()->FindDetectorConfiguration(fDetectorConfigurationForAlignmentName.Data());
    }
    else {
      AliFatal(Form("Wrong reference detector configuration %s for %s alignment correction step",
          fDetectorConfigurationForAlignmentName.Data(),
          fDetectorConfiguration->GetName()));
    }
  }
  else {
    AliFatal(Form("Missing reference detector configuration for %s alignment correction step",
        fDetectorConfiguration->GetName()));
  }

  Int_t nNoOfHarmonics = fDetectorConfiguration->GetNoOfHarmonics();
  Int_t *harmonicsMap = new Int_t[nNoOfHarmonics];
  /* make sure the alignment harmonic processing is active */
  fDetectorConfiguration->ActivateHarmonic(fHarmonicForAlignment);
  /* in both configurations */
  fDetectorConfigurationForAlignment->ActivateHarmonic(fHarmonicForAlignment);
  /* and now create the corrected Qn vector */
  fDetectorConfiguration->GetHarmonicMap(harmonicsMap);
  fCorrectedQnVector = new AliQnCorrectionsQnVector(szCorrectedQnVectorName, nNoOfHarmonics, harmonicsMap);
  fInputQnVector = fDetectorConfiguration->GetPreviousCorrectedQnVector(this);
  delete [] harmonicsMap;
}

/// Asks for support histograms creation
///
/// Allocates the histogram objects and creates the calibration histograms.
///
/// Process concurrency requires Calibration Histograms creation for all
/// concurrent processes but not for Input Histograms so, we delete previously
/// allocated ones.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorAlignment::CreateSupportHistograms(TList *list) {

  TString histoNameAndTitle = TString::Format("%s %s#times%s ",
      szSupportHistogramName,
      fDetectorConfiguration->GetName(),
      fDetectorConfigurationForAlignment->GetName());

  if (fInputHistograms != NULL) delete fInputHistograms;
  fInputHistograms = new AliQnCorrectionsProfileCorrelationComponents((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      fDetectorConfiguration->GetEventClassVariablesSet());
  fInputHistograms->SetNoOfEntriesThreshold(fMinNoOfEntriesToValidate);
  fCalibrationHistograms = new AliQnCorrectionsProfileCorrelationComponents((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      fDetectorConfiguration->GetEventClassVariablesSet());

  fCalibrationHistograms->CreateCorrelationComponentsProfileHistograms(list);
  return kTRUE;
}

/// Attaches the needed input information to the correction step
/// \param list list where the inputs should be found
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsQnVectorAlignment::AttachInput(TList *list) {

  if (fInputHistograms->AttachHistograms(list)) {
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
Bool_t AliQnCorrectionsQnVectorAlignment::CreateQAHistograms(TList *list) {

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
Bool_t AliQnCorrectionsQnVectorAlignment::CreateNveQAHistograms(TList *list) {

  fQANotValidatedBin = new AliQnCorrectionsHistogramSparse(
      TString::Format("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()).Data(),
      TString::Format("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()).Data(),
      fDetectorConfiguration->GetEventClassVariablesSet());
  fQANotValidatedBin->CreateHistogram(list);
  return kTRUE;
}

/// Processes the correction step
///
/// Apply the correction step
/// \return kTRUE if the correction step was applied
Bool_t AliQnCorrectionsQnVectorAlignment::ProcessCorrections(const Float_t *variableContainer) {
  switch (fState) {
  case QCORRSTEP_calibration:
    /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
    /* we have not perform any correction yet */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the correction if the current Qn vector is good enough */
    /* logging */
    AliInfo(TString::Format("Alignment process in detector %s with reference %s: applying correction.",
        fDetectorConfiguration->GetName(),
        fDetectorConfigurationForAlignment->GetName()).Data());
    if (fDetectorConfiguration->GetCurrentQnVector()->IsGoodQuality()) {
      /* we get the properties of the current Qn vector but its name */
      fCorrectedQnVector->Set(fDetectorConfiguration->GetCurrentQnVector(),kFALSE);

      /* let's check the correction histograms */
      Long64_t bin = fInputHistograms->GetBin(variableContainer);
      if (fInputHistograms->BinContentValidated(bin)) {
        /* the bin content is validated so, apply the correction */
        Double_t XX  = fInputHistograms->GetXXBinContent(bin);
        Double_t YY  = fInputHistograms->GetYYBinContent(bin);
        Double_t XY  = fInputHistograms->GetXYBinContent(bin);
        Double_t YX  = fInputHistograms->GetYXBinContent(bin);
        Double_t eXY = fInputHistograms->GetXYBinError(bin);
        Double_t eYX = fInputHistograms->GetYXBinError(bin);

        Double_t deltaPhi = - TMath::ATan2((XY-YX),(XX+YY)) * (1.0 / fHarmonicForAlignment);

        /* significant correction? */
        if (!(TMath::Sqrt((XY-YX)*(XY-YX)/(eXY*eXY+eYX*eYX)) < 2.0)) {
          Int_t harmonic = fDetectorConfiguration->GetCurrentQnVector()->GetFirstHarmonic();
          while (harmonic != -1) {
            fCorrectedQnVector->SetQx(harmonic,
                fDetectorConfiguration->GetCurrentQnVector()->Qx(harmonic) * TMath::Cos(((Double_t) harmonic) * deltaPhi)
                + fDetectorConfiguration->GetCurrentQnVector()->Qy(harmonic) * TMath::Sin (((Double_t) harmonic) * deltaPhi));
            fCorrectedQnVector->SetQy(harmonic,
                fDetectorConfiguration->GetCurrentQnVector()->Qy(harmonic) * TMath::Cos(((Double_t) harmonic) * deltaPhi)
                - fDetectorConfiguration->GetCurrentQnVector()->Qx(harmonic) * TMath::Sin (((Double_t) harmonic) * deltaPhi));
            harmonic = fDetectorConfiguration->GetCurrentQnVector()->GetNextHarmonic(harmonic);
          }
        } /* if the correction is not significant we leave the Q vector untouched */
      } /* if the correction bin is not validated we leave the Q vector untouched */
      else {
        if (fQANotValidatedBin != NULL) fQANotValidatedBin->Fill(variableContainer, 1.0);
      }
    }
    else {
      /* not done! input Q vector with bad quality */
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
/// Collect data for the correction step.
/// \return kTRUE if the correction step was applied
Bool_t AliQnCorrectionsQnVectorAlignment::ProcessDataCollection(const Float_t *variableContainer) {
  switch (fState) {
  case QCORRSTEP_calibration:
    /* logging */
    AliInfo(TString::Format("Alignment process in detector %s with reference %s: collecting data.",
        fDetectorConfiguration->GetName(),
        fDetectorConfigurationForAlignment->GetName()).Data());
    /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
    if ((fInputQnVector->IsGoodQuality()) &&
        (fDetectorConfigurationForAlignment->GetCurrentQnVector()->IsGoodQuality())) {
      fCalibrationHistograms->FillXX(variableContainer,
          fInputQnVector->Qx(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qx(fHarmonicForAlignment) );
      fCalibrationHistograms->FillXY(variableContainer,
          fInputQnVector->Qx(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qy(fHarmonicForAlignment));
      fCalibrationHistograms->FillYX(variableContainer,
          fInputQnVector->Qy(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qx(fHarmonicForAlignment) );
      fCalibrationHistograms->FillYY(variableContainer,
          fInputQnVector->Qy(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qy(fHarmonicForAlignment));
    }
    /* we have not perform any correction yet */
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* logging */
    AliInfo(TString::Format("Alignment process in detector %s with reference %s: collecting data.",
        fDetectorConfiguration->GetName(),
        fDetectorConfigurationForAlignment->GetName()).Data());
    /* collect the data needed to further produce correction parameters if both current Qn vectors are good enough */
    if ((fInputQnVector->IsGoodQuality()) &&
        (fDetectorConfigurationForAlignment->GetCurrentQnVector()->IsGoodQuality())) {
      fCalibrationHistograms->FillXX(variableContainer,
          fInputQnVector->Qx(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qx(fHarmonicForAlignment) );
      fCalibrationHistograms->FillXY(variableContainer,
          fInputQnVector->Qx(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qy(fHarmonicForAlignment));
      fCalibrationHistograms->FillYX(variableContainer,
          fInputQnVector->Qy(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qx(fHarmonicForAlignment) );
      fCalibrationHistograms->FillYY(variableContainer,
          fInputQnVector->Qy(fHarmonicForAlignment)
          * fDetectorConfigurationForAlignment->GetCurrentQnVector()->Qy(fHarmonicForAlignment));
    }
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the correction if the current Qn vector is good enough */
    /* provide QA info if required */
    if (fQAQnAverageHistogram != NULL) {
      Int_t harmonic = fCorrectedQnVector->GetFirstHarmonic();
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
void AliQnCorrectionsQnVectorAlignment::ClearCorrectionStep() {

  fCorrectedQnVector->Reset();
}

/// Reports if the correction step is being applied
/// Returns TRUE if in the proper state for applying the correction step
/// \return TRUE if the correction step is being applied
Bool_t AliQnCorrectionsQnVectorAlignment::IsBeingApplied() const {
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
Bool_t AliQnCorrectionsQnVectorAlignment::ReportUsage(TList *calibrationList, TList *applyList) {
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

