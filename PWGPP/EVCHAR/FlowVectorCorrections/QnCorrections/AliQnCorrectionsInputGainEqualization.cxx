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

/// \file AliQnCorrectionsInputGainEqualization.cxx
/// \brief Implementation of procedures for gain equalization on input data.
#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsProfileChannelizedIngress.h"
#include "AliQnCorrectionsProfileChannelized.h"
#include "AliQnCorrectionsHistogramChannelizedSparse.h"
#include "AliQnCorrectionsDetectorConfigurationChannels.h"
#include "AliLog.h"
#include "AliQnCorrectionsInputGainEqualization.h"

const Float_t  AliQnCorrectionsInputGainEqualization::fMinimumSignificantValue = 1E-6;
const Int_t AliQnCorrectionsInputGainEqualization::fDefaultMinNoOfEntries = 2;
const char *AliQnCorrectionsInputGainEqualization::szCorrectionName = "Gain equalization";
const char *AliQnCorrectionsInputGainEqualization::szKey = "CCCC";
const char *AliQnCorrectionsInputGainEqualization::szSupportHistogramName = "Multiplicity";
const char *AliQnCorrectionsInputGainEqualization::szQAHistogramName = "QA Multiplicity";
const char *AliQnCorrectionsInputGainEqualization::szQANotValidatedHistogramName = "GE NvE";

/// Default value for the shift parameter
#define GAINEQUALIZATION_SHIFTDEFAULT 0.0
/// Default value for the scale parameter
#define GAINEQUALIZATION_SCALEDEFAULT 1.0


/// \cond CLASSIMP
ClassImp(AliQnCorrectionsInputGainEqualization);
/// \endcond

/// Default constructor
/// Passes to the base class the identity data for the Gain equalization correction step
AliQnCorrectionsInputGainEqualization::AliQnCorrectionsInputGainEqualization() :
    AliQnCorrectionsCorrectionOnInputData(szCorrectionName, szKey) {
  fInputHistograms = NULL;
  fCalibrationHistograms = NULL;
  fQAMultiplicityBefore = NULL;
  fQAMultiplicityAfter = NULL;
  fQANotValidatedBin = NULL;
  fEqualizationMethod = GEQUAL_noEqualization;
  fShift = GAINEQUALIZATION_SHIFTDEFAULT;
  fScale = GAINEQUALIZATION_SCALEDEFAULT;
  fUseChannelGroupsWeights = kFALSE;
  fHardCodedWeights = NULL;
  fMinNoOfEntriesToValidate = fDefaultMinNoOfEntries;
}

/// Default destructor
/// Releases the memory taken
AliQnCorrectionsInputGainEqualization::~AliQnCorrectionsInputGainEqualization() {
  if (fInputHistograms != NULL)
    delete fInputHistograms;
  if (fCalibrationHistograms != NULL)
    delete fCalibrationHistograms;
  if (fQAMultiplicityBefore != NULL)
    delete fQAMultiplicityBefore;
  if (fQAMultiplicityAfter != NULL)
    delete fQAMultiplicityAfter;
  if (fQANotValidatedBin != NULL)
    delete fQANotValidatedBin;
}

/// Attaches the needed input information to the correction step
///
/// If the attachment succeeded asks for hard coded group weights to
/// the detector configuration
/// \param list list where the inputs should be found
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsInputGainEqualization::AttachInput(TList *list) {
  AliQnCorrectionsDetectorConfigurationChannels *ownerConfiguration =
      static_cast<AliQnCorrectionsDetectorConfigurationChannels *>(fDetectorConfiguration);
  if (fInputHistograms->AttachHistograms(list,
      ownerConfiguration->GetUsedChannelsMask(), ownerConfiguration->GetChannelsGroups())) {
    fState = QCORRSTEP_applyCollect;
    fHardCodedWeights = ownerConfiguration->GetHardCodedGroupWeights();
    return kTRUE;
  }
  return kFALSE;
}

/// Asks for support data structures creation
///
/// Does nothing for the time being
void AliQnCorrectionsInputGainEqualization::CreateSupportDataStructures() {

}

/// Asks for support histograms creation
///
/// Allocates the histogram objects and creates the calibration histograms.
/// The histograms are constructed with standard deviation error calculation
/// for the proper behavior of the gain equalization.
///
/// Process concurrency requires Calibration Histograms creation for all c
/// concurrent processes but not for Input Histograms so, we delete previously
/// allocated ones.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsInputGainEqualization::CreateSupportHistograms(TList *list) {

  TString histoNameAndTitle = Form("%s %s",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());

AliQnCorrectionsDetectorConfigurationChannels *ownerConfiguration =
      static_cast<AliQnCorrectionsDetectorConfigurationChannels *>(fDetectorConfiguration);
  if (fInputHistograms != NULL) delete fInputHistograms;
  fInputHistograms = new AliQnCorrectionsProfileChannelizedIngress((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      ownerConfiguration->GetEventClassVariablesSet(),ownerConfiguration->GetNoOfChannels(), "s");
  fInputHistograms->SetNoOfEntriesThreshold(fMinNoOfEntriesToValidate);
  fCalibrationHistograms = new AliQnCorrectionsProfileChannelized((const char *) histoNameAndTitle, (const char *) histoNameAndTitle,
      ownerConfiguration->GetEventClassVariablesSet(),ownerConfiguration->GetNoOfChannels(), "s");
  fCalibrationHistograms->CreateProfileHistograms(list,
      ownerConfiguration->GetUsedChannelsMask(), ownerConfiguration->GetChannelsGroups());
  return kTRUE;
}

/// Asks for QA histograms creation
///
/// Allocates the histogram objects and creates the QA histograms.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsInputGainEqualization::CreateQAHistograms(TList *list) {
  TString beforeName = Form("%s %s",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());
  beforeName += "Before";
  TString beforeTitle = Form("%s %s",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());
  beforeTitle += " before gain equalization";
  TString afterName = Form("%s %s",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());
  afterName += "After";
  TString afterTitle = Form("%s %s",
      szSupportHistogramName,
      fDetectorConfiguration->GetName());
  afterTitle += " after gain equalization";
  AliQnCorrectionsDetectorConfigurationChannels *ownerConfiguration =
      static_cast<AliQnCorrectionsDetectorConfigurationChannels *>(fDetectorConfiguration);
  fQAMultiplicityBefore = new AliQnCorrectionsProfileChannelized(
      (const char *) beforeName,
      (const char *) beforeTitle,
      ownerConfiguration->GetEventClassVariablesSet(),ownerConfiguration->GetNoOfChannels());
  fQAMultiplicityBefore->CreateProfileHistograms(list,
      ownerConfiguration->GetUsedChannelsMask(), ownerConfiguration->GetChannelsGroups());
  fQAMultiplicityAfter = new AliQnCorrectionsProfileChannelized(
      (const char *) afterName,
      (const char *) afterTitle,
      ownerConfiguration->GetEventClassVariablesSet(),ownerConfiguration->GetNoOfChannels());
  fQAMultiplicityAfter->CreateProfileHistograms(list,
      ownerConfiguration->GetUsedChannelsMask(), ownerConfiguration->GetChannelsGroups());
  fQANotValidatedBin = new AliQnCorrectionsHistogramChannelizedSparse(
      Form("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()),
      Form("%s %s", szQANotValidatedHistogramName, fDetectorConfiguration->GetName()),
      ownerConfiguration->GetEventClassVariablesSet(),
      ownerConfiguration->GetNoOfChannels());
  fQANotValidatedBin->CreateChannelizedHistogram(list, ownerConfiguration->GetUsedChannelsMask());
  return kTRUE;
}

/// Processes the correction step
///
/// Data are always taken from the data bank from the equalized weights
/// allowing chaining of input data corrections so, caution must be taken to be
/// sure that, on initialising, weight and equalized weight match
/// \return kTRUE if the correction step was applied
Bool_t AliQnCorrectionsInputGainEqualization::Process(const Float_t *variableContainer) {
  switch (fState) {
  case QCORRSTEP_calibration:
    /* collect the data needed to further produce equalization parameters */
    for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
      AliQnCorrectionsDataVectorChannelized *dataVector =
          static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
      fCalibrationHistograms->Fill(variableContainer, dataVector->GetId(), dataVector->EqualizedWeight());
    }
    return kFALSE;
    break;
  case QCORRSTEP_applyCollect:
    /* collect the data needed to further produce equalization parameters */
    for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
      AliQnCorrectionsDataVectorChannelized *dataVector =
          static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
      fCalibrationHistograms->Fill(variableContainer, dataVector->GetId(), dataVector->EqualizedWeight());
    }
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the equalization */
    /* collect QA data if asked */
    if (fQAMultiplicityBefore != NULL) {
      for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
        AliQnCorrectionsDataVectorChannelized *dataVector =
            static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
        fQAMultiplicityBefore->Fill(variableContainer, dataVector->GetId(), dataVector->EqualizedWeight());
      }
    }
    /* store the equalized weights in the data vector bank according to equalization method */
    switch (fEqualizationMethod) {
    case GEQUAL_noEqualization:
      for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
        AliQnCorrectionsDataVectorChannelized *dataVector =
            static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
        dataVector->SetEqualizedWeight(dataVector->EqualizedWeight());
      }
      break;
    case GEQUAL_averageEqualization:
      for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
        AliQnCorrectionsDataVectorChannelized *dataVector =
            static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
        Long64_t bin = fInputHistograms->GetBin(variableContainer, dataVector->GetId());
        if (fInputHistograms->BinContentValidated(bin)) {
          Float_t average = fInputHistograms->GetBinContent(bin);
          /* let's handle the potential group weights usage */
          Float_t groupweight = 1.0;
          if (fUseChannelGroupsWeights) {
            groupweight = fInputHistograms->GetGrpBinContent(fInputHistograms->GetGrpBin(variableContainer, dataVector->GetId()));
          }
          else {
            if (fHardCodedWeights != NULL) {
              groupweight = fHardCodedWeights[dataVector->GetId()];
            }
          }
          if (fMinimumSignificantValue < average)
            dataVector->SetEqualizedWeight((dataVector->EqualizedWeight() / average) * groupweight);
          else
            dataVector->SetEqualizedWeight(0.0);
        }
        else {
          fQANotValidatedBin->Fill(variableContainer, dataVector->GetId(), 1.0);
        }
      }
      break;
    case GEQUAL_widthEqualization:
      for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
        AliQnCorrectionsDataVectorChannelized *dataVector =
            static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
        Long64_t bin = fInputHistograms->GetBin(variableContainer, dataVector->GetId());
        if (fInputHistograms->BinContentValidated(bin)) {
          Float_t average = fInputHistograms->GetBinContent(fInputHistograms->GetBin(variableContainer, dataVector->GetId()));
          Float_t width = fInputHistograms->GetBinError(fInputHistograms->GetBin(variableContainer, dataVector->GetId()));
          /* let's handle the potential group weights usage */
          Float_t groupweight = 1.0;
          if (fUseChannelGroupsWeights) {
            groupweight = fInputHistograms->GetGrpBinContent(fInputHistograms->GetGrpBin(variableContainer, dataVector->GetId()));
          }
          else {
            if (fHardCodedWeights != NULL) {
              groupweight = fHardCodedWeights[dataVector->GetId()];
            }
          }
          if (fMinimumSignificantValue < average)
            dataVector->SetEqualizedWeight((fShift + fScale * (dataVector->EqualizedWeight() - average) / width) * groupweight);
          else
            dataVector->SetEqualizedWeight(0.0);
        }
        else {
          fQANotValidatedBin->Fill(variableContainer, dataVector->GetId(), 1.0);
        }
      }
      break;
    }
    /* collect QA data if asked */
    if (fQAMultiplicityAfter != NULL) {
      for(Int_t ixData = 0; ixData < fDetectorConfiguration->GetInputDataBank()->GetEntriesFast(); ixData++){
        AliQnCorrectionsDataVectorChannelized *dataVector =
            static_cast<AliQnCorrectionsDataVectorChannelized *>(fDetectorConfiguration->GetInputDataBank()->At(ixData));
        fQAMultiplicityAfter->Fill(variableContainer, dataVector->GetId(), dataVector->EqualizedWeight());
      }
    }
    break;
  }
  return kTRUE;
}

/// Report on correction usage
/// Correction step should incorporate its name in calibration
/// list if it is producing information calibration in the ongoing
/// step and in the apply list if it is applying correction in
/// the ongoing step.
/// \param calibrationList list containing the correction steps producing calibration information
/// \param applyList list containing the correction steps applying corrections
/// \return kTRUE if the correction step is being applied
Bool_t AliQnCorrectionsInputGainEqualization::ReportUsage(TList *calibrationList, TList *applyList) {
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
  }
  return kTRUE;
}

