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

/// \file AliQnCorrectionsDetector.cxx
/// \brief Detector class implementation

#include "AliQnCorrectionsDetector.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsDetector);
/// \endcond


/// Default constructor
AliQnCorrectionsDetector::AliQnCorrectionsDetector() : TNamed(),
    fConfigurations(),
    fDataVectorAcceptedConfigurations() {

  fDetectorId = -1;
  fDataVectorAcceptedConfigurations.SetOwner(kFALSE);
  fCorrectionsManager = NULL;
}

/// Normal constructor
/// \param name the name of the detector
/// \param id detector Id
AliQnCorrectionsDetector::AliQnCorrectionsDetector(const char *name, Int_t id) :
    TNamed(name,name),
    fConfigurations(),
    fDataVectorAcceptedConfigurations() {

  fDetectorId = id;
  fDataVectorAcceptedConfigurations.SetOwner(kFALSE);
  fCorrectionsManager = NULL;
}

/// Default destructor
/// The detector class does not own anything
AliQnCorrectionsDetector::~AliQnCorrectionsDetector() {

}

/// Asks for support data structures creation
///
/// The request is transmitted to the attached detector configurations
void AliQnCorrectionsDetector::CreateSupportDataStructures() {

  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->CreateSupportDataStructures();
  }
}

/// Asks for support histograms creation
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetector::CreateSupportHistograms(TList *list) {
  Bool_t retValue = kTRUE;
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    retValue = retValue && (fConfigurations.At(ixConfiguration)->CreateSupportHistograms(list));
  }
  return retValue;
}

/// Asks for QA histograms creation
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetector::CreateQAHistograms(TList *list) {
  Bool_t retValue = kTRUE;
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    retValue = retValue && (fConfigurations.At(ixConfiguration)->CreateQAHistograms(list));
  }
  return retValue;
}

/// Asks for attaching the needed input information to the correction steps
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the input information should be found
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetector::AttachCorrectionInputs(TList *list) {
  /* TODO: do we need to fine tune the list passed according to the detector? */
  Bool_t retValue = kTRUE;
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    Bool_t ret = fConfigurations.At(ixConfiguration)->AttachCorrectionInputs(list);
    retValue = retValue && ret;
  }
  return retValue;
}

/// Stores the framework manager pointer and transmits it to the incorporated detector configurations if any
///
/// \param manager the framework manager
void AliQnCorrectionsDetector::AttachCorrectionsManager(AliQnCorrectionsManager *manager) {
  fCorrectionsManager = manager;

  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntries(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->AttachCorrectionsManager(manager);
  }
}

/// Adds a new detector configuration to the current detector
///
/// Raise an execution error if the configuration detector reference
/// is not empty and if the detector configuration
/// is already incorporated to the detector.
/// \param detectorConfiguration pointer to the configuration to be added
void AliQnCorrectionsDetector::AddDetectorConfiguration(AliQnCorrectionsDetectorConfigurationBase *detectorConfiguration) {
  if (detectorConfiguration->GetDetector() != NULL) {
    AliFatal(Form("You are adding %s detector configuration of detector Id %d to detector Id %d. FIX IT, PLEASE.",
        detectorConfiguration->GetName(),
        detectorConfiguration->GetDetector()->GetId(),
        GetId()));
    return;
  }

  if (fConfigurations.FindObject(detectorConfiguration->GetName())) {
    AliFatal(Form("You are trying to add twice %s detector configuration to detector Id %d. FIX IT, PLEASE.",
        detectorConfiguration->GetName(),
        GetId()));
    return;
  }
  detectorConfiguration->SetDetectorOwner(this);
  detectorConfiguration->AttachCorrectionsManager(fCorrectionsManager);
  fConfigurations.Add(detectorConfiguration);
}

/// Searches for a concrete detector configuration by name
/// \param name the name of the detector configuration to find
/// \return pointer to the found detector configuration (NULL if not found)
AliQnCorrectionsDetectorConfigurationBase *AliQnCorrectionsDetector::FindDetectorConfiguration(const char *name) {
  return (AliQnCorrectionsDetectorConfigurationBase *) fConfigurations.FindObject(name);
}

/// Include the the list of Qn vector associated to the detector
/// into the passed list
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the corrected Qn vector should be added
void AliQnCorrectionsDetector::IncludeQnVectors(TList *list) {

  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->IncludeQnVectors(list);
  }
}

/// Include the name of each detector configuration into the passed list
///
/// \param list the list where to incorporate detector configurations name
void AliQnCorrectionsDetector::FillDetectorConfigurationNameList(TList *list) const {
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    list->Add(new TObjString(fConfigurations.At(ixConfiguration)->GetName()));
  }
}

/// Include the name of the input correction steps on each detector
/// configuration into the passed list
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the corrected Qn vector should be added
void AliQnCorrectionsDetector::FillOverallInputCorrectionStepList(TList *list) const {
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->FillOverallInputCorrectionStepList(list) ;
  }
}

/// Include the name of the Qn vector correction steps on each detector
/// configuration into the passed list
///
/// The request is transmitted to the attached detector configurations
/// \param list list where the corrected Qn vector should be added
void AliQnCorrectionsDetector::FillOverallQnVectorCorrectionStepList(TList *list) const {
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->FillOverallQnVectorCorrectionStepList(list) ;
  }
}


/// Provide information about assigned corrections on each of the detector configurations
///
/// The request is transmitted to the attached detector configurations
/// \param steps list for incorporating the list of assigned correction steps
/// \param calib list for incorporating the list of steps in calibrating status
/// \param apply list for incorporating the list of steps in applying status
void AliQnCorrectionsDetector::ReportOnCorrections(TList *steps, TList *calib, TList *apply) const {
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->ReportOnCorrections(steps, calib, apply) ;
  }
}


