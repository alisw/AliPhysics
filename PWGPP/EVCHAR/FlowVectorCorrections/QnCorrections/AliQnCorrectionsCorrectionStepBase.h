#ifndef ALIQNCORRECTIONS_CORRECTIONSTEPBASE_H
#define ALIQNCORRECTIONS_CORRECTIONSTEPBASE_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsCorrectionStepBase.h
/// \brief Base class for the support of the different correction steps within Q vector correction framework
///

#include <TNamed.h>
#include <TList.h>

class AliQnCorrectionsDetectorConfigurationBase;
class AliQnCorrectionsDetectorConfigurationChannels;
class AliQnCorrectionsQnVector;

/// \class AliQnCorrectionsCorrectionStepBase
/// \brief Base class for correction steps
///
/// Each correction has a name and a key. The name identifies it
/// in an open way while the key is used to codify its position
/// in an ordered list of consecutive corrections.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 05, 2016

class AliQnCorrectionsCorrectionStepBase : public TNamed {
public:
  /// \typedef QnCorrectionStepStatus
  /// \brief The class of the id of the correction steps states
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  /// When referring as "data being collected" means that the needed data
  /// for producing new correction parameters are being collected.
  typedef enum {
    QCORRSTEP_calibration,         ///< the correction step is in calibration mode collecting data
    QCORRSTEP_apply,               ///< the correction step is being applied
    QCORRSTEP_applyCollect,        ///< the correction step is being applied and data are being collected
  } QnCorrectionStepStatus;

  friend class AliQnCorrectionsDetectorConfigurationBase;
  AliQnCorrectionsCorrectionStepBase();
  AliQnCorrectionsCorrectionStepBase(const char *name, const char *key);
  virtual ~AliQnCorrectionsCorrectionStepBase();

  /// Gets the correction ordering key
  const char *GetKey() const { return (const char *) fKey; }
  Bool_t Before(const AliQnCorrectionsCorrectionStepBase *correction);

  /// Informs when the detector configuration has been attached to the framework manager
  /// Basically this allows interaction between the different framework sections at configuration time
  /// Pure virtual function
  virtual void AttachedToFrameworkManager() = 0;
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
  /// Asks for QA histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateQAHistograms(TList *list) = 0;
  /// Asks for non validated entries QA histograms creation
  ///
  /// Pure virtual function
  /// \param list list where the histograms should be incorporated for its persistence
  /// \return kTRUE if everything went OK
  virtual Bool_t CreateNveQAHistograms(TList *list) = 0;
  /// Processes the correction step
  ///
  /// Pure virtual function
  /// \return kTRUE if everything went OK
  virtual Bool_t Process(const Float_t *variableContainer) = 0;
  /// Include the new corrected Qn vector into the passed list
  ///
  /// Pure virtual function
  /// \param list list where the corrected Qn vector should be added
  virtual void IncludeCorrectedQnVector(TList *list) = 0;
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
protected:
  /// Stores the detector configuration owner
  /// \param detectorConfiguration the detector configuration owner
  void SetConfigurationOwner(AliQnCorrectionsDetectorConfigurationBase *detectorConfiguration)
  { fDetectorConfiguration = detectorConfiguration; }

  QnCorrectionStepStatus fState;                                  ///< the state in which the correction step is
  AliQnCorrectionsDetectorConfigurationBase *fDetectorConfiguration; ///< pointer to the detector configuration owner
  TString fKey;                                                   ///< the correction key that codifies order information

private:
  /// Copy constructor
  /// Not allowed. Forced private.
  AliQnCorrectionsCorrectionStepBase(AliQnCorrectionsCorrectionStepBase &);
  /// Assignment operator
  /// Not allowed. Forced private.
  AliQnCorrectionsCorrectionStepBase& operator= (const AliQnCorrectionsCorrectionStepBase &);

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsCorrectionStepBase, 1);
/// \endcond
};

#endif // ALIQNCORRECTIONS_CORRECTIONSTEPBASE_H
