#ifndef ALIQNCORRECTIONS_DETECTOR_H
#define ALIQNCORRECTIONS_DETECTOR_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsDetector.h
/// \brief Detector and detector configuration classes for Q vector correction framework
///

#include "AliQnCorrectionsDetectorConfigurationBase.h"
#include "AliQnCorrectionsDetectorConfigurationsSet.h"

class AliQnCorrectionsDetectorConfigurationsSet;
class AliQnCorrectionsDetectorConfigurationBase;
class AliQnCorrectionsDetector;

/// \class AliQnCorrectionsDetector
/// \brief Detector class within Q vector correction framework
///
/// The roles of the detector class are: to store its unique name and Id,
/// and to store and handle the list of the different
/// configurations defined for the involved detector.
///
/// The detector Id should only be obtained at creation time and
/// the object, once created, does not allow its modification.
///
/// The detector object is in the path of the whole control flow and
/// as such it should distribute the different commands to the
/// defined detector configurations.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class AliQnCorrectionsDetector : public TNamed {
public:
  AliQnCorrectionsDetector();
  AliQnCorrectionsDetector(const char *name, Int_t id);
  virtual ~AliQnCorrectionsDetector();

  /// Gets the detector Id
  ///
  /// \return detector Id
  Int_t GetId() { return fDetectorId; }

  void CreateSupportDataStructures();
  Bool_t CreateSupportHistograms(TList *list);
  Bool_t CreateQAHistograms(TList *list);
  Bool_t AttachCorrectionInputs(TList *list);
  Bool_t ProcessCorrections(const Float_t *variableContainer);
  void IncludeQnVectors(TList *list);

  /// Gets the name of the detector configuration at index that accepted last data vector
  /// \param index the position in the list of accepted data vector configuration
  /// \return the configuration name
  const char *GetAcceptedDataDetectorConfigurationName(Int_t index) const
  { return fDataVectorAcceptedConfigurations.At(index)->GetName(); }

  void AddDetectorConfiguration(AliQnCorrectionsDetectorConfigurationBase *detectorConfiguration);
  AliQnCorrectionsDetectorConfigurationBase *FindDetectorConfiguration(const char *name);
  void FillDetectorConfigurationNameList(TList *list) const;
  void FillOverallInputCorrectionStepList(TList *list) const;
  void FillOverallQnVectorCorrectionStepList(TList *list) const;
  virtual void ReportOnCorrections(TList *steps, TList *calib, TList *apply) const;

  Int_t AddDataVector(const Float_t *variableContainer, Double_t phi, Double_t weight = 1.0, Int_t channelId = -1);

  virtual void ClearDetector();

private:
  Int_t fDetectorId;            ///< detector Id
  AliQnCorrectionsDetectorConfigurationsSet fConfigurations;  ///< the set of configurations defined for this detector
  AliQnCorrectionsDetectorConfigurationsSet fDataVectorAcceptedConfigurations; ///< the set of configurations that accepted a data vector

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsDetector, 1);
/// \endcond
};

/// New data vector for the detector
/// The request is transmitted to the attached detector configurations.
/// The current content of the variable bank is passed in order to check
/// for optional cuts tha define the detector configurations.
/// \param variableContainer pointer to the variable content bank
/// \param phi azimuthal angle
/// \param weight the weight of the data vector
/// \param channelId the channel Id that originates the data vector
/// \return the number of detector configurations that accepted and stored the data vector
inline Int_t AliQnCorrectionsDetector::AddDataVector(const Float_t *variableContainer, Double_t phi, Double_t weight, Int_t channelId) {
  fDataVectorAcceptedConfigurations.Clear();
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    Bool_t ret = fConfigurations.At(ixConfiguration)->AddDataVector(variableContainer, phi, weight, channelId);
    if (ret) {
      fDataVectorAcceptedConfigurations.Add(fConfigurations.At(ixConfiguration));
    }
  }
  return fDataVectorAcceptedConfigurations.GetEntries();
}

/// Ask for processing corrections for the involved detector
///
/// The request is transmitted to the attached detector configurations
/// \return kTRUE if everything went OK
inline Bool_t AliQnCorrectionsDetector::ProcessCorrections(const Float_t *variableContainer) {
  Bool_t retValue = kTRUE;

  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    Bool_t ret = fConfigurations.At(ixConfiguration)->ProcessCorrections(variableContainer);
    retValue = retValue && ret;
  }
  return retValue;
}

/// Clean the detector to accept a new event
///
/// Transfers the order to the detector configurations
inline void AliQnCorrectionsDetector::ClearDetector() {
  /* transfer the order to the Q vector corrections */
  for (Int_t ixConfiguration = 0; ixConfiguration < fConfigurations.GetEntriesFast(); ixConfiguration++) {
    fConfigurations.At(ixConfiguration)->ClearConfiguration();
  }
}

#endif // ALIQNCORRECTIONS_DETECTOR_H
