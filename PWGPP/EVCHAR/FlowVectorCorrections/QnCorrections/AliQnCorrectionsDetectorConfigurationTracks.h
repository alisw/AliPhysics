#ifndef ALIQNCORRECTIONS_DETECTORCONFTRACKS_H
#define ALIQNCORRECTIONS_DETECTORCONFTRACKS_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsDetectorConfigurationTracks.h
/// \brief Track detector configuration class for Q vector correction framework
///

#include "AliQnCorrectionsDataVector.h"
#include "AliQnCorrectionsDetectorConfigurationBase.h"

/// \class AliQnCorrectionsDetectorConfigurationTracks
/// \brief Track detector configuration within Q vector correction framework
///
/// A track detector within the Q vector correction framework is defined
/// as one for which its data vectors only involve azimuthal angles and a
/// potential weight. Apart from that no other input data calibration is
/// available.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class AliQnCorrectionsDetectorConfigurationTracks :
    public AliQnCorrectionsDetectorConfigurationBase {
public:
  friend class AliQnCorrectionsCorrectionStepBase;
  friend class AliQnCorrectionsDetector;
  AliQnCorrectionsDetectorConfigurationTracks();
  AliQnCorrectionsDetectorConfigurationTracks(const char *name,
      AliQnCorrectionsEventClassVariablesSet *eventClassesVariables,
      Int_t nNoOfHarmonics,
      Int_t *harmonicMap = NULL);
  virtual ~AliQnCorrectionsDetectorConfigurationTracks();

  virtual void AttachCorrectionsManager(AliQnCorrectionsManager *manager);

  virtual void CreateSupportDataStructures();
  virtual Bool_t CreateSupportHistograms(TList *list);
  virtual Bool_t CreateQAHistograms(TList *list);
  virtual Bool_t AttachCorrectionInputs(TList *list);

  /// Ask for processing corrections for the involved detector configuration
  ///
  /// The request is transmitted to the Q vector correction steps
  /// \return kTRUE if everything went OK
  virtual Bool_t ProcessCorrections(const Float_t *variableContainer);
  virtual Bool_t AddDataVector(const Float_t *variableContainer, Double_t phi, Double_t weight = 1.0, Int_t channelId = -1);

  virtual void BuildQnVector();
  virtual void IncludeQnVectors(TList *list);
  virtual void FillOverallInputCorrectionStepList(TList *list) const;
  virtual void FillOverallQnVectorCorrectionStepList(TList *list) const;
  virtual void ReportOnCorrections(TList *steps, TList *calib, TList *apply) const;

  /// Checks if the current content of the variable bank applies to
  /// the detector configuration
  /// \param variableContainer pointer to the variable content bank
  /// \return kTRUE if the current content applies to the configuration
  virtual Bool_t IsSelected(const Float_t *variableContainer)
    { return ((fCuts != NULL) ? fCuts->IsSelected(variableContainer) : kTRUE); }
  /// wrong call for this class invoke base class behavior
  virtual Bool_t IsSelected(const Float_t *variableContainer, Int_t nChannel)
  { return AliQnCorrectionsDetectorConfigurationBase::IsSelected(variableContainer,nChannel); }

  virtual void ClearConfiguration();

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsDetectorConfigurationTracks, 1);
/// \endcond
};

/// New data vector for the detector configuration.
/// A check is made to see if the current variable bank content passes
/// the associated cuts. If so, the data vector is stored.
/// \param variableContainer pointer to the variable content bank
/// \param phi azimuthal angle
/// \param weight the weight associated to the data vector. For track detector is usually one.
/// \param id the Id associated to the data vector. For track detector configurations could represent the track id.
/// \return kTRUE if the data vector was accepted and stored
inline Bool_t AliQnCorrectionsDetectorConfigurationTracks::AddDataVector(
    const Float_t *variableContainer, Double_t phi, Double_t weight, Int_t id) {
  if (IsSelected(variableContainer)) {
    /// add the data vector to the bank
    new (fDataVectorBank->ConstructedAt(fDataVectorBank->GetEntriesFast()))
        AliQnCorrectionsDataVector(id, phi, weight);
    return kTRUE;
  }
  return kFALSE;
}

/// Clean the configuration to accept a new event
///
/// Transfers the order to the Q vector correction steps and
/// cleans the own Q vector and the input data vector bank
/// for accepting the next event.
inline void AliQnCorrectionsDetectorConfigurationTracks::ClearConfiguration() {
  /* transfer the order to the Q vector corrections */
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    fQnVectorCorrections.At(ixCorrection)->ClearCorrectionStep();
  }
  /* clean the own Q vectors */
  fPlainQnVector.Reset();
  fCorrectedQnVector.Reset();
  /* and now clear the the input data bank */
  fDataVectorBank->Clear("C");
}

/// Builds Qn vectors before Q vector corrections but
/// considering the chosen calibration method.
/// Remember, this configuration does not have a channelized
/// approach so, the built Q vectors are the ones to be used for
/// subsequent corrections.
inline void AliQnCorrectionsDetectorConfigurationTracks::BuildQnVector() {
  fTempQnVector.Reset();

  for(Int_t ixData = 0; ixData < fDataVectorBank->GetEntriesFast(); ixData++){
    AliQnCorrectionsDataVector *dataVector = static_cast<AliQnCorrectionsDataVector *>(fDataVectorBank->At(ixData));
    fTempQnVector.Add(dataVector->Phi(), dataVector->Weight());
  }
  /* check the quality of the Qn vector */
  fTempQnVector.CheckQuality();
  fTempQnVector.Normalize(fQnNormalizationMethod);
  fPlainQnVector.Set(&fTempQnVector, kFALSE);
  fCorrectedQnVector.Set(&fTempQnVector, kFALSE);
}


/// Ask for processing corrections for the involved detector configuration
///
/// The request is transmitted to the Q vector correction steps.
/// The first not applied correction step breaks the loop and kFALSE is returned
/// \return kTRUE if all correction steps were applied
inline Bool_t AliQnCorrectionsDetectorConfigurationTracks::ProcessCorrections(const Float_t *variableContainer) {
  /* first we build the Q vector with the chosen calibration */
  BuildQnVector();

  /* then we transfer the request to the Q vector correction steps */
  /* the loop is broken when a correction step has not been applied */
  for (Int_t ixCorrection = 0; ixCorrection < fQnVectorCorrections.GetEntries(); ixCorrection++) {
    if (fQnVectorCorrections.At(ixCorrection)->Process(variableContainer))
      continue;
    else
      return kFALSE;
  }
  /* all correction steps were applied */
  return kTRUE;
}

#endif // ALIQNCORRECTIONS_DETECTORCONFTRACKS_H
