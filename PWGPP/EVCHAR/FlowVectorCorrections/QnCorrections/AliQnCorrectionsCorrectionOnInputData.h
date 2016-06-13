#ifndef ALIQNCORRECTIONS_CORRECTIONONINPUTDATA_H
#define ALIQNCORRECTIONS_CORRECTIONONINPUTDATA_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsCorrectionOnInputData.h
/// \brief Correction steps on input data support within Q vector correction framework
///

#include <TNamed.h>
#include <TList.h>

#include "AliQnCorrectionsCorrectionStepBase.h"

/// \class AliQnCorrectionsCorrectionOnInputData
/// \brief Base class for correction steps applied to input data
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 05, 2016

class AliQnCorrectionsCorrectionOnInputData: public AliQnCorrectionsCorrectionStepBase {
public:
  friend class AliQnCorrectionsDetectorConfigurationChannels;
  AliQnCorrectionsCorrectionOnInputData();
  AliQnCorrectionsCorrectionOnInputData(const char *name, const char *key);
  virtual ~AliQnCorrectionsCorrectionOnInputData();

  /// Attaches the needed input information to the correction step
  ///
  /// Pure virtual function
  /// \param list list where the inputs should be found
  /// \return kTRUE if everything went OK
  virtual Bool_t AttachInput(TList *list) = 0;
  /// Asks for support data structures creation
  ///
  /// Pure virtual function
  virtual void CreateSupportDataStructures() = 0;
  /// Asks for support histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateSupportHistograms(TList *list) = 0;
  /// Processes the correction step
  ///
  /// Pure virtual function
  /// \return kTRUE if everything went OK
  virtual Bool_t Process(const Float_t *variableContainer) = 0;
  /// Include the new corrected Qn vector into the passed list
  ///
  /// Does nothing. Not applicable for corrections on input data
  /// \param list list where the corrected Qn vector should be added
  virtual void IncludeCorrectedQnVector(TList *list)  {}
  /// Clean the correction to accept a new event
  /// Pure virtual function
  virtual void ClearCorrectionStep() = 0;
  /// Report on correction usage
  /// Pure virtual function
  /// Correction step should incorporate its name in calibration
  /// list if it is producing information calibration in the ongoing
  /// step and in the apply list if it is applying correction in
  /// the ongoing step.
  /// \param calibrationList list containing the correction steps producing calibration information
  /// \param applyList list containing the correction steps applying corrections
  /// \return kTRUE if the correction step is being applied
  virtual Bool_t ReportUsage(TList *calibrationList, TList *applyList) = 0;
/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsCorrectionOnInputData, 1);
/// \endcond
};

#endif // ALIQNCORRECTIONS_CORRECTIONONINPUTDATA_H
