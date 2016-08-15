#ifndef ALIQNCORRECTIONS_DETECTORCONFIGBASE_H
#define ALIQNCORRECTIONS_DETECTORCONFIGBASE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsDetectorConfigurationBase.h
/// \brief The base of a concrete detector configuration (sub-detector) within Q vector correction framework
///

#include <TObject.h>
#include <TList.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TH3.h>
#include "AliQnCorrectionsCutsSet.h"
#include "AliQnCorrectionsCorrectionsSetOnInputData.h"
#include "AliQnCorrectionsCorrectionsSetOnQvector.h"
#include "AliQnCorrectionsEventClassVariablesSet.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliQnCorrectionsQnVectorBuild.h"

class AliQnCorrectionsDetectorConfigurationsSet;
class AliQnCorrectionsDetector;
class AliQnCorrectionsManager;

/// \class AliQnCorrectionsDetectorConfigurationBase
/// \brief The base of a concrete detector configuration within Q vector correction framework
///
/// The detector configuration shapes a detector with a concrete
/// set of cuts to make it the target of a Q vector correction process.
///
/// It receives the data input stream and build the corresponding Q vector associated
/// to it for each processing request.
///
/// As such, it incorporates the set of corrections to carry on the input data
/// and the set of corrections to perform on the produced Qn vector. It always stores
/// the plain Qn vector produced after potential input data corrections and the
/// Qn vector that incorporates the latest Qn vector correction step.
///
/// It also incorporates the equivalent support for Q2n vectors which could be the
/// seed for future Q(m,n) support.
///
/// It receives at construction time the set of event classes variables and the
/// detector reference. The reference of the detector should only be obtained at
/// creation time and the detector configuration object, once created, does not
/// allow its modification.
///
/// The class is a base class for further refined detector configurations.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 08, 2016

class AliQnCorrectionsDetectorConfigurationBase : public TNamed {
public:
  friend class AliQnCorrectionsCorrectionStepBase;
  friend class AliQnCorrectionsDetector;
  AliQnCorrectionsDetectorConfigurationBase();
  AliQnCorrectionsDetectorConfigurationBase(const char *name,
      AliQnCorrectionsEventClassVariablesSet *eventClassesVariables,
      Int_t nNoOfHarmonics,
      Int_t *harmonicMap = NULL);
  virtual ~AliQnCorrectionsDetectorConfigurationBase();

  /// Sets the set of cuts for the detector configuration
  /// \param cuts the set of cuts
  void SetCuts(AliQnCorrectionsCutsSet *cuts)
  { fCuts = cuts; }
  /// Sets the normalization method for Q vectors
  /// \param method the Qn vector normalizatio method
  void SetQVectorNormalizationMethod(AliQnCorrectionsQnVector::QnVectorNormalizationMethod method)
  { fQnNormalizationMethod = method; }

public:
  /// Stores the detector reference
  /// \param detector the detector owner
  void SetDetectorOwner(AliQnCorrectionsDetector *detector)
  { fDetector = detector; }
  /// Gets the detector reference
  ///
  /// \return detector pointer
  AliQnCorrectionsDetector *GetDetector() { return fDetector; }
  /// Stores the framework manager pointer
  /// Pure virtual function
  /// \param manager the framework manager
  virtual void AttachCorrectionsManager(AliQnCorrectionsManager *manager) = 0;
public:
  /// Get the input data bank.
  /// Makes it available for input corrections steps.
  /// \return pointer to the input data bank
  TClonesArray *GetInputDataBank()
  { return fDataVectorBank; }
  /// Get the event class variables set
  /// Makes it available for corrections steps
  /// \return pointer to the event class variables set
  AliQnCorrectionsEventClassVariablesSet &GetEventClassVariablesSet()
  { return *fEventClassVariables; }
  /// Get the current Qn vector
  /// Makes it available for subsequent correction steps.
  /// It could have already supported previous correction steps
  /// \return pointer to the current Qn vector instance
  AliQnCorrectionsQnVector *GetCurrentQnVector()
  { return &fCorrectedQnVector; }
  const AliQnCorrectionsQnVector *GetPreviousCorrectedQnVector(AliQnCorrectionsCorrectionOnQvector *correctionOnQn) const;
  Bool_t IsCorrectionStepBeingApplied(const char *step) const;
  /// Get the current Q2n vector
  /// Makes it available for subsequent correction steps.
  /// It could have already supported previous correction steps
  /// \return pointer to the current Q2n vector instance
  AliQnCorrectionsQnVector *GetCurrentQ2nVector()
  { return &fCorrectedQ2nVector; }
  /// Get the plain Qn vector
  /// Makes it available for correction steps which need it.
  /// \return pointer to the plain Qn vector instance
  AliQnCorrectionsQnVector *GetPlainQnVector()
  { return &fPlainQnVector; }
  /// Get the plain Q2n vector
  /// Makes it available for correction steps which need it.
  /// \return pointer to the plain Qn vector instance
  AliQnCorrectionsQnVector *GetPlainQ2nVector()
  { return &fPlainQ2nVector; }
  /// Update the current Qn vector
  /// Update towards what is the latest values of the Qn vector after executing a
  /// correction step to make it available to further steps.
  /// \param newQnVector the new values for the Qn vector
  /// \param changename kTRUE by default to keep track of the subsequent Qn vector corrections
  void UpdateCurrentQnVector(AliQnCorrectionsQnVector *newQnVector, Bool_t changename = kTRUE)
  { fCorrectedQnVector.Set(newQnVector, changename); }
  /// Update the current Q2n vector
  /// Update towards what is the latest values of the Q2n vector after executing a
  /// correction step to make it available to further steps.
  /// \param newQ2nVector the new values for the Q2n vector
  /// \param changename kTRUE by default to keep track of the subsequent Q2n vector corrections
  void UpdateCurrentQ2nVector(AliQnCorrectionsQnVector *newQ2nVector, Bool_t changename = kTRUE)
  { fCorrectedQ2nVector.Set(newQ2nVector, changename); }
  /// Get the number of harmonics handled by the detector configuration
  /// \return the number of handled harmonics
  Int_t GetNoOfHarmonics() const
  { return fCorrectedQnVector.GetNoOfHarmonics(); }
  /// Get the harmonics map handled by the detector configuration
  /// \param store pointer to the memory for storing the harmonics map
  void GetHarmonicMap(Int_t *store) const
  { fCorrectedQnVector.GetHarmonicsMap(store); }
  /// Get the pointer to the framework manager
  /// \return the stored pointer to the corrections framework
  AliQnCorrectionsManager *GetCorrectionsManager() const { return fCorrectionsManager; }
  /// Get if the detector configuration is own by a tracking detector
  /// Pure virtual function
  /// \return TRUE if it is a tracking detector configuration
  virtual Bool_t GetIsTrackingDetector() const = 0;
public:
  /// Asks for support data structures creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  virtual void CreateSupportDataStructures() = 0;

  /// Asks for support histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateSupportHistograms(TList *list) = 0;

  /// Asks for QA histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateQAHistograms(TList *list) = 0;

  /// Asks for non validated entries QA histograms creation
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateNveQAHistograms(TList *list) = 0;

  /// Asks for attaching the needed input information to the correction steps
  ///
  /// The request is transmitted to the different corrections.
  /// Pure virtual function
  /// \param list list where the input information should be found
  /// \return kTRUE if everything went OK
  virtual Bool_t AttachCorrectionInputs(TList *list) = 0;
  /// Perform after calibration histograms attach actions
  /// It is used to inform the different correction step that
  /// all conditions for running the network are in place so
  /// it is time to check if their requirements are satisfied
  ///
  /// Pure virtual function
  virtual void AfterInputsAttachActions() = 0;
  /// Ask for processing corrections for the involved detector configuration
  ///
  /// Pure virtual function.
  /// The request is transmitted to the correction steps
  /// \return kTRUE if everything went OK
  virtual Bool_t ProcessCorrections(const Float_t *variableContainer) = 0;
  /// Ask for processing corrections data collection for the involved detector configuration
  ///
  /// Pure virtual function.
  /// The request is transmitted to the correction steps
  /// \return kTRUE if everything went OK
  virtual Bool_t ProcessDataCollection(const Float_t *variableContainer) = 0;
  virtual void ActivateHarmonic(Int_t harmonic);
  virtual void AddCorrectionOnQnVector(AliQnCorrectionsCorrectionOnQvector *correctionOnQn);
  virtual void AddCorrectionOnInputData(AliQnCorrectionsCorrectionOnInputData *correctionOnInputData);

  /// Builds Qn vector before Q vector corrections but
  /// considering the chosen calibration method.
  /// Pure virtual function
  virtual void BuildQnVector() = 0;
  /// Include the list of associated Qn vectors into the passed list
  ///
  /// Pure virtual function
  /// \param list list where the Qn vectors list should be added
  virtual void IncludeQnVectors(TList *list) = 0;
  /// Include only one instance of each input correction step
  /// in execution order
  ///
  /// Pure virtual function
  /// \param list list where the correction steps should be incorporated
  virtual void FillOverallInputCorrectionStepList(TList *list) const = 0;
  /// Include only one instance of each Qn vector correction step
  /// in execution order
  ///
  /// Pure virtual function
  /// \param list list where the correction steps should be incorporated
  virtual void FillOverallQnVectorCorrectionStepList(TList *list) const = 0;
  /// Provide information about assigned corrections
  ///
  /// Pure virtual function
  /// \param steps list for incorporating the list of assigned correction steps
  /// \param calib list for incorporating the list of steps in calibrating status
  /// \param apply list for incorporating the list of steps in applying status
  virtual void ReportOnCorrections(TList *steps, TList *calib, TList *apply) const = 0;

  /// New data vector for the detector configuration
  /// Pure virtual function
  /// \param variableContainer pointer to the variable content bank
  /// \param phi azimuthal angle
  /// \param weight the weight of the data vector
  /// \param channelId the channel Id that originates the data vector
  /// \return kTRUE if the data vector was accepted and stored
  virtual Bool_t AddDataVector(const Float_t *variableContainer, Double_t phi, Double_t weight = 1.0, Int_t channelId = -1) = 0;

  virtual Bool_t IsSelected(const Float_t *variableContainer);
  virtual Bool_t IsSelected(const Float_t *variableContainer, Int_t nChannel);

  /// Clean the configuration to accept a new event
  /// Pure virtual function
  virtual void ClearConfiguration() = 0;

private:
  AliQnCorrectionsDetector *fDetector;    ///< pointer to the detector that owns the configuration
protected:
  static const char *szPlainQnVectorName; ///< the name of the Qn plain, not corrected Qn vectors
  /// set of cuts that define the detector configuration
  AliQnCorrectionsManager *fCorrectionsManager; /// the framework manager pointer
  AliQnCorrectionsCutsSet *fCuts;         //->
/// The default initial size of data vectors banks
#define INITIALDATAVECTORBANKSIZE 100000
  TClonesArray *fDataVectorBank;        //!<! input data for the current process / event
  AliQnCorrectionsQnVector fPlainQnVector;     ///< Qn vector from the post processed input data
  AliQnCorrectionsQnVector fPlainQ2nVector;     ///< Q2n vector from the post processed input data
  AliQnCorrectionsQnVector fCorrectedQnVector; ///< Qn vector after subsequent correction steps
  AliQnCorrectionsQnVector fCorrectedQ2nVector; ///< Q2n vector after subsequent correction steps
  AliQnCorrectionsQnVectorBuild fTempQnVector; ///< temporary Qn vector for efficient Q vector building
  AliQnCorrectionsQnVectorBuild fTempQ2nVector; ///< temporary Qn vector for efficient Q vector building
  AliQnCorrectionsQnVector::QnVectorNormalizationMethod fQnNormalizationMethod; ///< the method for Q vector normalization
  AliQnCorrectionsCorrectionsSetOnQvector fQnVectorCorrections; ///< set of corrections to apply on Q vectors
  /// set of variables that define event classes
  AliQnCorrectionsEventClassVariablesSet    *fEventClassVariables; //->

private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliQnCorrectionsDetectorConfigurationBase(const AliQnCorrectionsDetectorConfigurationBase &);
  /// Assignment operator
  /// Not allowed. Forced private.
  AliQnCorrectionsDetectorConfigurationBase& operator= (const AliQnCorrectionsDetectorConfigurationBase &);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsDetectorConfigurationBase, 3);
/// \endcond
};

#endif // ALIQNCORRECTIONS_DETECTORCONFIGBASE_H
