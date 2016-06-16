#ifndef ALIQNCORRECTIONS_QNVECTORALIGNMENT_H
#define ALIQNCORRECTIONS_QNVECTORALIGNMENT_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsQnVectorAlignment.h
/// \brief Definition of the class that implements Qn vectors rotations for detector alignment corrections.
///
/// Qn vector alignment is applied on ongoing Qn vector from the involved detector
/// configuration. It is supposed that such a Qn vector has already been built after input data
/// equalization and, potentially, Qn vector recentering, and is with proper normalization for
/// the recentering correction being applied.
///
/// The alignment is applied according to:
/// \f[
///        Q' = \mathcal{R}(\Delta \phi) Q
/// \f]
/// where the rotation angle \f$ \Delta \phi \f$ is given by
/// \f[
///        \Delta \phi_n = - \frac{1}{m} \tan^{-1}
///          \left(\frac{\langle{Q_n}_X{Q^A_m}_Y\rangle - \langle{Q_n}_Y{Q^A_m}_X\rangle}
///                    {\langle{Q_n}_X{Q^A_m}_X\rangle + \langle{Q_n}_Y{Q^A_m}_Y\rangle}\right)
/// \f]
/// with  \f$ \langle \cdots \rangle \f$ as an average over events in a given event class, \f$ A \f$
/// the detector configuration chosen as alignment reference and \f$ m \f$ the harmonic selected
/// for alignment.
///
/// So, options configurable by the user are the detector configuration to use for alignment
/// and the harmonic number also to be used.
///
/// Qn vector rotation(and width equalization) is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the needed correction histograms.
///
/// Correction and data collecting during calibration is performed for all harmonics
/// defined within the involved detector configuration

#include "AliQnCorrectionsCorrectionOnQvector.h"

/// \class AliQnCorrectionsQnVectorAlignment
/// \brief Encapsulates Qn vector rotation for alignment correction
///
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Apr 27, 2016
///
/// Alignment is applied on ongoing Qn vector from the involved detector
/// configuration.
///
/// The alignment is applied according to:
/// \f[
///        Q' = \mathcal{R}(\Delta \phi) Q
/// \f]
/// where the rotation angle \f$ \Delta \phi \f$ is given by
/// \f[
///        \Delta \phi_n = - \frac{1}{m} \tan^{-1}
///          \left(\frac{\langle{Q_n}_X{Q^A_m}_Y\rangle - \langle{Q_n}_Y{Q^A_m}_X\rangle}
///                    {\langle{Q_n}_X{Q^A_m}_X\rangle + \langle{Q_n}_Y{Q^A_m}_Y\rangle}\right)
/// \f]
/// with  \f$ \langle \cdots \rangle \f$ as an average over events in a given event class, \f$ A \f$
/// the detector configuration chosen as alignment reference and \f$ m \f$ the harmonic selected
/// for alignment.
///
/// Alignment is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper correction histograms that will
/// provide the required averages per event class.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the correction histograms.
///
/// Correction and data collecting during calibration is performed for all harmonics
/// defined within the involved detector configuration


class AliQnCorrectionsQnVectorAlignment : public AliQnCorrectionsCorrectionOnQvector {
public:
  AliQnCorrectionsQnVectorAlignment();
  ~AliQnCorrectionsQnVectorAlignment();

  /// Set the harmonic number used for alignment
  /// \param harmonic harmonic number
  void SetHarmonicNumberForAlignment(Int_t harmonic)
  { fHarmonicForAlignment = harmonic; }
  void SetReferenceConfigurationForAlignment(const char *name);
  /// Set the minimum number of entries for calibration histogram bin content validation
  /// \param nNoOfEntries the number of entries threshold
  void SetNoOfEntriesThreshold(Int_t nNoOfEntries) { fMinNoOfEntriesToValidate = nNoOfEntries; }

  virtual void AttachedToFrameworkManager();
  virtual Bool_t AttachInput(TList *list);
  virtual void CreateSupportDataStructures();
  virtual Bool_t CreateSupportHistograms(TList *list);
  virtual Bool_t CreateQAHistograms(TList *list);

  virtual Bool_t Process(const Float_t *variableContainer);
  virtual void ClearCorrectionStep();
  virtual Bool_t ReportUsage(TList *calibrationList, TList *applyList);

private:
  static const Int_t fDefaultMinNoOfEntries;         ///< the minimum number of entries for bin content validation
  static const char *szCorrectionName;               ///< the name of the correction step
  static const char *szKey;                          ///< the key of the correction step for ordering purpose
  static const char *szSupportHistogramName;         ///< the name and title for support histograms
  static const char *szCorrectedQnVectorName;        ///< the name of the Qn vector after applying the correction
  static const char *szQANotValidatedHistogramName;  ///< the name and title for bin not validated QA histograms
  AliQnCorrectionsProfileCorrelationComponents *fInputHistograms; //!<! the histogram with calibration information
  AliQnCorrectionsProfileCorrelationComponents *fCalibrationHistograms; //!<! the histogram for building calibration information
  AliQnCorrectionsHistogram *fQANotValidatedBin;    //!<! the histogram with non validated bin information

  Int_t  fHarmonicForAlignment;              ///< the harmonic number to be used for Qn vector alignment correction
  TString fDetectorConfigurationForAlignmentName; ///< storage for the name of the reference detector configuration for alignment correction
  AliQnCorrectionsDetectorConfigurationBase *fDetectorConfigurationForAlignment; ///< pointer to the detector configuration used as reference for alingment
  Int_t fMinNoOfEntriesToValidate;              ///< number of entries for bin content validation threshold

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsQnVectorAlignment, 2);
/// \endcond
};

#endif // ALIQNCORRECTIONS_QNVECTORALIGNMENT_H
