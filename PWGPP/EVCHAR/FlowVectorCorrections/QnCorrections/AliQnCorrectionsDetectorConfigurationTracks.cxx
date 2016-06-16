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

/// \file AliQnCorrectionsDetectorConfigurationTracks.cxx
/// \brief Implementation of the track detector configuration class

#include "AliQnCorrectionsDetectorConfigurationTracks.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsDetectorConfigurationTracks);
/// \endcond

/// Default constructor
AliQnCorrectionsDetectorConfigurationTracks::AliQnCorrectionsDetectorConfigurationTracks() : AliQnCorrectionsDetectorConfigurationBase() {

}

/// Normal constructor
/// Allocates the data vector bank.
/// \param name the name of the detector configuration
/// \param eventClassesVariables the set of event classes variables
/// \param nNoOfHarmonics the number of harmonics that must be handled
/// \param harmonicMap an optional ordered array with the harmonic numbers
AliQnCorrectionsDetectorConfigurationTracks::AliQnCorrectionsDetectorConfigurationTracks(const char *name,
      AliQnCorrectionsEventClassVariablesSet *eventClassesVariables,
      Int_t nNoOfHarmonics,
      Int_t *harmonicMap) :
          AliQnCorrectionsDetectorConfigurationBase(name, eventClassesVariables, nNoOfHarmonics, harmonicMap) {

}

/// Default destructor
/// Memory taken is released by the parent class destructor
AliQnCorrectionsDetectorConfigurationTracks::~AliQnCorrectionsDetectorConfigurationTracks() {

}

/// Stores the framework manager pointer
/// Orders the base class to store the correction manager and informs
/// the Qn vector corrections they are now attached to the framework
/// \param manager the framework manager
void AliQnCorrectionsDetectorConfigurationTracks::AttachCorrectionsManager(AliQnCorrectionsManager *manager) {
  fCorrectionsManager = manager;

  if (manager != NULL) {
    for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
      fQnVectorCorrections.At(ixCorrection)->AttachedToFrameworkManager();
    }
  }
}

/// Asks for support data structures creation
///
/// The input data vector bank is allocated and the request is
/// transmitted to the Q vector corrections.
void AliQnCorrectionsDetectorConfigurationTracks::CreateSupportDataStructures() {

  /* this is executed in the remote node so, allocate the data bank */
  fDataVectorBank = new TClonesArray("AliQnCorrectionsDataVector", INITIALDATAVECTORBANKSIZE);

  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    fQnVectorCorrections.At(ixCorrection)->CreateSupportDataStructures();
  }
}

/// Asks for support histograms creation
///
/// The request is transmitted to the Q vector corrections.
///
/// A new histograms list is created for the detector configuration and incorporated
/// to the passed list. Then the new list is passed to the corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetectorConfigurationTracks::CreateSupportHistograms(TList *list) {
  Bool_t retValue = kTRUE;
  TList *detectorConfigurationList = new TList();
  detectorConfigurationList->SetName(this->GetName());
  detectorConfigurationList->SetOwner(kTRUE);
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    retValue = retValue && (fQnVectorCorrections.At(ixCorrection)->CreateSupportHistograms(detectorConfigurationList));
  }
  /* if list is empty delete it if not incorporate it */
  if (detectorConfigurationList->GetEntries() != 0) {
    list->Add(detectorConfigurationList);
  }
  else {
    delete detectorConfigurationList;
  }
  return retValue;
}

/// Asks for QA histograms creation
///
/// The request is transmitted to the Q vector corrections.
///
/// A new histograms list is created for the detector and incorporated
/// to the passed list. Then the new list is passed to the corrections.
/// \param list list where the histograms should be incorporated for its persistence
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetectorConfigurationTracks::CreateQAHistograms(TList *list) {
  Bool_t retValue = kTRUE;
  TList *detectorConfigurationList = new TList();
  detectorConfigurationList->SetName(this->GetName());
  detectorConfigurationList->SetOwner(kTRUE);
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    retValue = retValue && (fQnVectorCorrections.At(ixCorrection)->CreateQAHistograms(detectorConfigurationList));
  }
  /* if list is empty delete it if not incorporate it */
  if (detectorConfigurationList->GetEntries() != 0) {
    list->Add(detectorConfigurationList);
  }
  else {
    delete detectorConfigurationList;
  }
  return retValue;
}

/// Asks for attaching the needed input information to the correction steps
///
/// The detector list is extracted from the passed list and then
/// the request is transmitted to the Q vector corrections with the found list.
/// \param list list where the input information should be found
/// \return kTRUE if everything went OK
Bool_t AliQnCorrectionsDetectorConfigurationTracks::AttachCorrectionInputs(TList *list) {
  TList *detectorConfigurationList = (TList *) list->FindObject(this->GetName());
  if (detectorConfigurationList != NULL) {
    Bool_t retValue = kTRUE;
    for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
      retValue = retValue && (fQnVectorCorrections.At(ixCorrection)->AttachInput(detectorConfigurationList));
    }
    return retValue;
  }
  return kFALSE;
}

/// Include the the list of Qn vector associated to the detector configuration
/// into the passed list
///
/// A new list is created for the detector configuration and incorporated
/// to the passed list.
///
/// Always includes first the fully corrected Qn vector,
/// and then includes the plain Qn vector and asks to the different correction
/// steps to include their partially corrected Qn vectors.
/// The check if we are already there is because it could be late information
/// about the process name and then the correction histograms could still not
/// be attached and the constructed list does not contain the final Qn vectors.
/// \param list list where the corrected Qn vector should be added
void AliQnCorrectionsDetectorConfigurationTracks::IncludeQnVectors(TList *list) {

  /* we check whether we are already there and if so we clean it and go again */
  Bool_t bAlreadyThere;
  TList *detectorConfigurationList;
  if (list->FindObject(this->GetName())) {
    detectorConfigurationList = (TList*) list->FindObject(this->GetName());
    detectorConfigurationList->Clear();
    bAlreadyThere = kTRUE;
  }
  else {
    detectorConfigurationList = new TList();
    detectorConfigurationList->SetName(this->GetName());
    detectorConfigurationList->SetOwner(kFALSE);
    bAlreadyThere = kFALSE;
  }

  detectorConfigurationList->Add(&fCorrectedQnVector);
  detectorConfigurationList->Add(&fPlainQnVector);
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    fQnVectorCorrections.At(ixCorrection)->IncludeCorrectedQnVector(detectorConfigurationList);
  }
  if (!bAlreadyThere)
    list->Add(detectorConfigurationList);
}

/// Include only one instance of each input correction step
/// in execution order
///
/// There are not input correction so we do nothing
/// \param list list where the correction steps should be incorporated
void AliQnCorrectionsDetectorConfigurationTracks::FillOverallInputCorrectionStepList(TList *list) const {

}

/// Include only one instance of each Qn vector correction step
/// in execution order
///
/// The request is transmitted to the set of Qn vector corrections
/// \param list list where the correction steps should be incorporated
void AliQnCorrectionsDetectorConfigurationTracks::FillOverallQnVectorCorrectionStepList(TList *list) const {

  fQnVectorCorrections.FillOverallCorrectionsList(list);
}

/// Provide information about assigned corrections
///
/// We create three list which items they own, incorporate info from the
/// correction steps and add them to the passed list
/// \param steps list for incorporating the list of assigned correction steps
/// \param calib list for incorporating the list of steps in calibrating status
/// \param apply list for incorporating the list of steps in applying status
void AliQnCorrectionsDetectorConfigurationTracks::ReportOnCorrections(TList *steps, TList *calib, TList *apply) const {
  TList *mysteps = new TList();
  mysteps->SetOwner(kTRUE);
  mysteps->SetName(GetName());
  TList *mycalib = new TList();
  mycalib->SetOwner(kTRUE);
  mycalib->SetName(GetName());
  TList *myapply = new TList();
  myapply->SetOwner(kTRUE);
  myapply->SetName(GetName());

  /* incorporate Qn vector corrections */
  Bool_t keepIncorporating = kTRUE;
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    mysteps->Add(new TObjString(fQnVectorCorrections.At(ixCorrection)->GetName()));
    /* incorporate additional info if the step will be reached */
    if (keepIncorporating) {
      Bool_t keep = fQnVectorCorrections.At(ixCorrection)->ReportUsage(mycalib,myapply);
      keepIncorporating = keepIncorporating && keep;
    }
  }
  steps->Add(mysteps);
  calib->Add(mycalib);
  apply->Add(myapply);
}

