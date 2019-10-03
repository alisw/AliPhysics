#ifndef ALIQNCORRECTIONS_MANAGER_H
#define ALIQNCORRECTIONS_MANAGER_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsManager.h
/// \brief Class that orchestrates the Q vector correction framework

/// \class AliQnCorrectionsManager
/// \brief Class orchestrates the Q vector correction framework
///
/// It should be only one instance of it and behaves as the anchor
/// point between the Q vector correction framework and the external
/// run time environment.
///
/// It owns the list of the detectors incorporated to the framework
/// and distributes among them the different commands according to
/// the analysis phase.
///
/// To improve performance a mapping between internal detector address
/// and external detector id is maintained.
///
/// When the framework is in the calibration phase there are no complete
/// calibration information available to fully implement the desired
/// correction on the input data and on the subsequent Q vector. During
/// this phase a list of support histograms is kept to build the
/// needed quantities to produce the intended calibration information.
/// At the same time a list of the available calibration histograms is
/// as well kept to perform the feasible corrections.
///
/// When, finally, the whole calibration information were available the
/// support histograms list will not need to be present and the framework
/// will be performing just the intended corrections on input data and on
/// the subsequent Q vector.
///
/// The Qn correction manager supports an optional list of concurrent
/// processes that will each contribute to the final merged results. Each
/// process has its own name and it is expected to get contribution of
/// different running instances. At merging time, only the contributions
/// from instances of the same process must be merged.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 16, 2016


#include <TObject.h>
#include <TList.h>
#include <TTree.h>
#include "AliQnCorrectionsDetector.h"

class AliQnCorrectionsManager : public TObject {
public:
  AliQnCorrectionsManager();
  virtual ~AliQnCorrectionsManager();

  /// Establishes the list of processes names
  /// \param names an array containing the processes names
  void SetListOfProcessesNames(TObjArray *names) { fProcessesNames = names; }
  void SetCurrentProcessListName(const char *name);
  void SetCalibrationHistogramsList(TFile *calibrationFile);
  /// Enables disables the filling of histograms for building correction parameters
  /// \param enable kTRUE for enabling histograms filling
  void SetShouldFillOutputHistograms(Bool_t enable = kTRUE) { fFillOutputHistograms = enable; }
  /// Enables disables the filling of QA histograms
  /// \param enable kTRUE for enabling QA histograms filling
  void SetShouldFillQAHistograms(Bool_t enable = kTRUE) { fFillQAHistograms = enable; }
  /// Enables disables the filling of non validated entries QA histograms
  /// \param enable kTRUE for enabling non validated entries QA histograms filling
  void SetShouldFillNveQAHistograms(Bool_t enable = kTRUE) { fFillNveQAHistograms = enable; }
  /// Enables disables the output of Qn vector on a TTree structure
  /// \param enable kTRUE for enabling Qn vector output into a TTree
  void SetShouldFillQnVectorTree(Bool_t enable = kTRUE) { fFillQnVectorTree = enable; }

  void AddDetector(AliQnCorrectionsDetector *detector);

  AliQnCorrectionsDetector *FindDetector(const char *name) const;
  AliQnCorrectionsDetector *FindDetector(Int_t id) const;
  AliQnCorrectionsDetectorConfigurationBase *FindDetectorConfiguration(const char *name) const;


  /// Gets a pointer to the data variables bank
  /// \return the pointer to the data container
  Float_t *GetDataContainer() { return fDataContainer; }

  /// Get whether the output histograms should be filled
  /// \return kTRUE if the output histograms should be filled
  Bool_t GetShouldFillOutputHistograms() const { return fFillOutputHistograms; }
  /// Get whether the QA histograms should be filled
  /// \return kTRUE if the QA histograms should be filled
  Bool_t GetShouldFillQAHistograms() const { return fFillQAHistograms; }
  /// Get whether the non validated entries QA histograms should be filled
  /// \return kTRUE if the non validated entries QA histograms should be filled
  Bool_t GetShouldFillNveQAHistograms() const { return fFillNveQAHistograms; }
  /// Get whether the Qn vector tree should be populated
  /// \return kTRUE if the Qn vector should be written into a TTree
  Bool_t GetShouldFillQnVectorTree() const { return fFillQnVectorTree; }
  /// Gets the output histograms list
  /// \return the list of histograms for building correction parameters
  TList *GetOutputHistogramsList() const { return fSupportHistogramsList; }
  /// Gets the QA histograms list
  /// \return the list of QA histograms
  TList *GetQAHistogramsList() const { return fQAHistogramsList; }
  /// Gets the non validated entries QA histograms list
  /// \return the list of QA histograms
  TList *GetNveQAHistogramsList() const { return fNveQAHistogramsList; }
  /// Gets the Qn vector tree
  /// \return the tree of histograms for building correction parameters
  TTree *GetQnVectorTree() const { return fQnVectorTree; }
  /// Gets the Qn vector tree
  /// \return the list of detector configurations Qn vectors
  TList *GetQnVectorList() const { return fQnVectorList; }
  const TList *GetDetectorQnVectorList(const char *subdetector) const;
  const AliQnCorrectionsQnVector *GetDetectorQnVector(const char *subdetector, const char *expectedstep = "latest", const char *altstep = "latest") const;
  /// Gets the name of the calibration histograms container
  /// \return the calibration histograms container name
  const char *GetCalibrationHistogramsContainerName() const
  { return szCalibrationHistogramsKeyName; }
  /// Gets the name of the calibration QA histograms container
  /// \return the calibration QA histograms container name
  const char *GetCalibrationQAHistogramsContainerName() const
  { return szCalibrationQAHistogramsKeyName; }
  /// Gets the name of the non validated calibration entries QA histograms container
  /// \return the calibration QA histograms container name
  const char *GetCalibrationNveQAHistogramsContainerName() const
  { return szCalibrationNveQAHistogramsKeyName; }


  void PrintFrameworkConfiguration() const;
  void InitializeQnCorrectionsFramework();
  Int_t AddDataVector(Int_t detectorId, Double_t phi, Double_t weight = 1.0, Int_t channelId = -1);
  const char *GetAcceptedDataDetectorConfigurationName(Int_t detectorId, Int_t index) const;
  void ProcessEvent();
  void ClearEvent();
  void FinalizeQnCorrectionsFramework();

private:
  static const Int_t nMaxNoOfDetectors;              ///< the highest detector id currently supported by the framework
  static const Int_t nMaxNoOfDataVariables;          ///< the maximum number of variables currently supported by the framework
  static const char *szCalibrationHistogramsKeyName; ///< the name of the key under which calibration histograms lists are stored
  static const char *szCalibrationQAHistogramsKeyName; ///< the name of the key under which calibration QA histograms lists are stored
  static const char *szCalibrationNveQAHistogramsKeyName; ///< the name of the key under which non validated calibration entries QA histograms lists are stored
  static const char *szDummyProcessListName;         ///< accepted temporary name before getting the definitive one
  static const char *szAllProcessesListName;         ///< the name of the list that collects data from all concurrent processes
  TList fDetectorsSet;                  ///< the list of detectors
  AliQnCorrectionsDetector **fDetectorsIdMap; //!<! map between external detector Id and internal detector
  Float_t *fDataContainer;              //!<! the data variables bank
  TList *fCalibrationHistogramsList;    ///< the list of the input calibration histograms
  TList *fSupportHistogramsList;        //!<! the list of the support histograms
  TList *fQAHistogramsList;             //!<! the list of QA histograms
  TList *fNveQAHistogramsList;          //!<! the list of not validated entries QA histograms
  TTree *fQnVectorTree;                 //!<! the tree to out Qn vectors
  TList *fQnVectorList;                 //!<! list that contains the current event corrected Qn vectors
  Bool_t fFillOutputHistograms;         ///< kTRUE if output histograms for building correction parameters must be filled
  Bool_t fFillQAHistograms;             ///< kTRUE if QA histograms must be filled
  Bool_t fFillNveQAHistograms;          ///< kTRUE if non validated entries QA histograms must be filled
  Bool_t fFillQnVectorTree;             ///< kTRUE if Qn vectors must be written in a TTree structure
  TString fProcessListName;             ///< the name of the list associated to the current process
  TObjArray *fProcessesNames;           ///< array with the list of processes names

private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliQnCorrectionsManager(const AliQnCorrectionsManager &);
  /// Assignment operator
  /// Not allowed. Forced private.
  AliQnCorrectionsManager& operator= (const AliQnCorrectionsManager &);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsManager, 5);
/// \endcond
};

/// New data vector for the framework
/// The request is transmitted to the passed detector together with
/// the current content of the variable bank.
/// \param detectorId id of the involved detector
/// \param phi azimuthal angle
/// \param weight the weight of the data vector
/// \param channelId the channel Id that originates the data vector
/// \return the number of detector configurations that accepted and stored the data vector
inline Int_t AliQnCorrectionsManager::AddDataVector(Int_t detectorId, Double_t phi, Double_t weight, Int_t channelId) {
  return fDetectorsIdMap[detectorId]->AddDataVector(fDataContainer, phi, weight, channelId);
}

/// Gets the name of the detector configuration at index that accepted last data vector
/// \param detectorId id of the involved detector
/// \param index the position in the list of accepted data vector configuration
/// \return the configuration name
inline const char *AliQnCorrectionsManager::GetAcceptedDataDetectorConfigurationName(Int_t detectorId, Int_t index) const {
  return fDetectorsIdMap[detectorId]->GetAcceptedDataDetectorConfigurationName(index);
}

/// Process the current event
///
/// The request is transmitted to the different detectors first for applying the different
/// correction steps and then to collect the correction steps data.
///
/// Must be called only when the whole data vectors for the event
/// have been incorporated to the framework.
inline void AliQnCorrectionsManager::ProcessEvent() {
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->ProcessCorrections(fDataContainer);
  }
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->ProcessDataCollection(fDataContainer);
  }
}

/// Clear the current event
///
/// The request is transmitted to the different detectors.
///
/// Must be called only at the end of each event to start processing the next one
inline void AliQnCorrectionsManager::ClearEvent() {
  for (Int_t ixDetector = 0; ixDetector < fDetectorsSet.GetEntries(); ixDetector++) {
    ((AliQnCorrectionsDetector *) fDetectorsSet.At(ixDetector))->ClearDetector();
  }
}

#endif // ALIQNCORRECTIONS_MANAGER_H
