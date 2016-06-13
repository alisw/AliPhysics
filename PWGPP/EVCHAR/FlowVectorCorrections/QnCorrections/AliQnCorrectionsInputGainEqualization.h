#ifndef ALIQNCORRECTIONS_INPUTGAINEQUALIZATION_H
#define ALIQNCORRECTIONS_INPUTGAINEQUALIZATION_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsInputGainEqualization.h
/// \brief Definition of the class that implements gain equalization of individual channels.
///
/// Gain equalization is applied on raw data from the involved detector.
/// Two procedures are implemented: average gain equalization and width equalization.
///
/// The average gain equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///        \mbox{M}'_{c,i} = \frac{\mbox{M}_{c,i}}{\langle\mbox{M}_{c}\rangle}
/// \f]
/// where  \f$\langle\mbox{M}_{c}\rangle\f$ is an average over events in a given event class
/// \f[
///        \langle\mbox{M}_{c}\rangle = \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}
/// \f]
///
/// The width equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///     \mbox{M}'_{c,i} = \mbox{A} + \mbox{B} \frac{\mbox{M}_{c,i} - \langle\mbox{M}_{c}\rangle}
///                               {\sigma_{{M}_{c}}}
/// \f]
/// with A and B are parameters that are the same for all channels and
/// \f$\sigma_{{M}_{c}}\f$ is the standard deviation of the signals of the channel \f$c\f$
/// for the considered event class
/// \f[
///        \sigma_{{M}_{c}} = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}^2_{c,i} -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}\right)^2}
/// \f]
///
/// The gain equalization is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper calibration histograms that will
/// provide the required averages per event class and channel.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the calibration histograms.

#include "AliQnCorrectionsCorrectionOnInputData.h"

class AliQnCorrectionsProfileChannelizedIngress;
class AliQnCorrectionsProfileChannelized;

/// \class AliQnCorrectionsInputGainEqualization
/// \brief Encapsulates the gain equalization on input data correction step
///
/// Gain equalization is applied on raw data from the involved detector.
/// Two procedures are implemented: average gain equalization and width equalization.
///
/// The average gain equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///        \mbox{M}'_{c,i} = \frac{\mbox{M}_{c,i}}{\langle\mbox{M}_{c}\rangle}
/// \f]
/// where  \f$\langle\mbox{M}_{c}\rangle\f$ is an average over events in a given event class
/// \f[
///        \langle\mbox{M}_{c}\rangle = \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}
/// \f]
///
/// The width equalization is applied for the signal \f$ \mbox{M}_{c,i} \f$ of each detector
/// channel \f$ c \f$ measured for event \f$ i \f$ according to:
/// \f[
///     \mbox{M}'_{c,i} = \mbox{A} + \mbox{B} \frac{\mbox{M}_{c,i} - \langle\mbox{M}_{c}\rangle}
///                               {\sigma_{{M}_{c}}}
/// \f]
/// with A and B are parameters that are the same for all channels and
/// \f$\sigma_{{M}_{c}}\f$ is the standard deviation of the signals of the channel \f$c\f$
/// for the considered event class
/// \f[
///        \sigma_{{M}_{c}} = \sqrt{
///          \frac{1}{\mbox{N}_{ev}} \Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}^2_{c,i} -
///          \frac{1}{\mbox{N}^2_{ev}} \left(\Sigma_{i}^{\mbox{N}_{ev}} \mbox{M}_{c,i}\right)^2}
/// \f]
/// At the class level A is known as the shift and B is known as the scale.
/// The class also provides support for group equalization where a group weight can be
/// extracted from the channels that conform a group or could be passed as hard coded group
/// weights at detector configuration definition.
///
/// The gain equalization is only applied if the class instance
/// is in the correction status. In order to be in that status the instance
/// should have been able to get the proper calibration histograms that will
/// provide the required averages per event class and channel.
/// If the class instance is not in the correction status then, it is
/// in the calibration one, collecting data for producing, once merged in a
/// further phase, the calibration histograms.

class AliQnCorrectionsInputGainEqualization : public AliQnCorrectionsCorrectionOnInputData {
public:
  /// \enum QnGainEqualizationMethod
  /// \brief The class of the id of the supported gain equalization methods
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  enum QnGainEqualizationMethod {
    GEQUAL_noEqualization,         ///< \f$ \mbox{M'} = \mbox{M}\f$
    GEQUAL_averageEqualization,    ///< \f$ \mbox{M}' = \frac{\mbox{M}}{\langle\mbox{M}\rangle} \f$
    GEQUAL_widthEqualization,      ///< \f$ \mbox{M}' = \mbox{A} + \mbox{B} \frac{\mbox{M} - \langle\mbox{M} \rangle}{\sigma_{{M}}} \f$
  };

  AliQnCorrectionsInputGainEqualization();
  ~AliQnCorrectionsInputGainEqualization();

  /// Stores the passed equalization method
  /// \param method the desired equalization method
  void SetEqualizationMethod(QnGainEqualizationMethod method)
  { fEqualizationMethod = method; }

  /// Set the shift (A) width equalization parameter
  /// \param shift the shift parameter value
  void SetShift(Float_t shift)
  { fShift = shift; }
  /// Set the scale (B) equalization parameter
  /// \param scale the scale parameter value
  void SetScale(Float_t scale)
  { fScale = scale; }
  /// Enable or disable the group weights extracted from channel multiplicity
  /// \param enable kTRUE / kFALSE for enable / disable it
  void SetUseChannelGroupsWeights(Bool_t enable)
  { fUseChannelGroupsWeights = enable; }

  virtual Bool_t AttachInput(TList *list);
  virtual void CreateSupportDataStructures();
  virtual Bool_t CreateSupportHistograms(TList *list);
  virtual Bool_t CreateQAHistograms(TList *list);

  virtual Bool_t Process(const Float_t *variableContainer);
  /// Clean the correction to accept a new event
  /// Does nothing for the time being
  virtual void ClearCorrectionStep() {}
  virtual Bool_t ReportUsage(TList *calibrationList, TList *applyList);

private:
  static const Float_t  fMinimumSignificantValue;     ///< the minimum value that will be considered as meaningful for processing
  static const char *szCorrectionName;               ///< the name of the correction step
  static const char *szKey;                          ///< the key of the correction step for ordering purpose
  static const char *szSupportHistogramName;         ///< the name and title for support histograms
  static const char *szQAHistogramName;              ///< the name and title for QA histograms
  AliQnCorrectionsProfileChannelizedIngress *fInputHistograms; //!<! the histogram with calibration information
  AliQnCorrectionsProfileChannelized *fCalibrationHistograms; //!<! the histogram for building calibration information
  AliQnCorrectionsProfileChannelized *fQAMultiplicityBefore;  //!<! the channel multiplicity histogram before gain equalization
  AliQnCorrectionsProfileChannelized *fQAMultiplicityAfter;   //!<! the channel multiplicity histogram after gain equalization
  QnGainEqualizationMethod fEqualizationMethod; ///< the selected equalization method

  Float_t fShift;                               ///< the shift (A) parameter for width equalization
  Float_t fScale;                               ///< the scale (B) parameter for width equalization
  Bool_t fUseChannelGroupsWeights;              ///< use group weights extracted from channel multiplicity
  const Float_t *fHardCodedWeights;             //!<! group hard coded weights stored in the detector configuration

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsInputGainEqualization, 1);
/// \endcond
};

#endif // ALIQNCORRECTIONS_INPUTGAINEQUALIZATION_H
