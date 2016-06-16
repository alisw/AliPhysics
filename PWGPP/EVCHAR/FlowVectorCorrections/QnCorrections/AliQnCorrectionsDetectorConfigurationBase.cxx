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

/// \file AliQnCorrectionsDetectorConfigurationBase.cxx
/// \brief Implementation of the base detector configuration class within Q vector correction framework

#include "AliQnCorrectionsDetectorConfigurationBase.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsDetectorConfigurationBase);
/// \endcond

const char *AliQnCorrectionsDetectorConfigurationBase::szPlainQnVectorName = "plain";


/// Default constructor
AliQnCorrectionsDetectorConfigurationBase::AliQnCorrectionsDetectorConfigurationBase() : TNamed(),
    fPlainQnVector(), fCorrectedQnVector(), fTempQnVector(), fQnVectorCorrections() {
  fDetector = NULL;
  fCorrectionsManager = NULL;
  fCuts = NULL;
  fDataVectorBank = NULL;
  fQnNormalizationMethod = AliQnCorrectionsQnVector::QVNORM_noCalibration;
  fEventClassVariables = NULL;
}

/// Normal constructor
/// \param name the name of the detector configuration
/// \param eventClassesVariables the set of event classes variables
/// \param nNoOfHarmonics the number of harmonics that must be handled
/// \param harmonicMap an optional ordered array with the harmonic numbers
AliQnCorrectionsDetectorConfigurationBase::AliQnCorrectionsDetectorConfigurationBase(const char *name,
      AliQnCorrectionsEventClassVariablesSet *eventClassesVariables,
      Int_t nNoOfHarmonics,
      Int_t *harmonicMap) :
          TNamed(name,name),
          fPlainQnVector(szPlainQnVectorName,nNoOfHarmonics, harmonicMap),
          fCorrectedQnVector(szPlainQnVectorName,nNoOfHarmonics, harmonicMap),
          fTempQnVector("temp",nNoOfHarmonics, harmonicMap),
          fQnVectorCorrections() {

  fDetector = NULL;
  fCorrectionsManager = NULL;
  fCuts = NULL;
  fDataVectorBank = NULL;
  fQnNormalizationMethod = AliQnCorrectionsQnVector::QVNORM_noCalibration;
  fEventClassVariables = eventClassesVariables;
}

/// Default destructor
/// Releases the memory which was taken or passed
AliQnCorrectionsDetectorConfigurationBase::~AliQnCorrectionsDetectorConfigurationBase() {
  if (fDataVectorBank != NULL) {
    delete fDataVectorBank;
  }
  if (fCuts != NULL) {
    delete fCuts;
  }
}

/// Incorporates the passed correction to the set of Q vector corrections
/// \param correctionOnQn the correction to add
void AliQnCorrectionsDetectorConfigurationBase::AddCorrectionOnQnVector(AliQnCorrectionsCorrectionOnQvector *correctionOnQn) {
  correctionOnQn->SetConfigurationOwner(this);
  fQnVectorCorrections.AddCorrection(correctionOnQn);
}

/// Incorporates the passed correction to the set of input data corrections
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param correctionOnInputData the correction to add
void AliQnCorrectionsDetectorConfigurationBase::AddCorrectionOnInputData(AliQnCorrectionsCorrectionOnInputData *correctionOnInputData) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a non channelized detector configuration to calibrate input data. FIX IT, PLEASE.",
      "AliQnCorrectionsDetectorConfigurationBase::AddCorrectionOnInputData()"));
}

/// Checks if the current content of the variable bank applies to
/// the detector configuration
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer pointer to the variable content bank
/// \return kTRUE if the current content applies to the configuration
Bool_t AliQnCorrectionsDetectorConfigurationBase::IsSelected(const Float_t *variableContainer) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a channelized detector configuration without passing the channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsDetectorConfigurationBase::IsSelected()"));
  return kFALSE;
}

/// Checks if the current content of the variable bank applies to
/// the detector configuration for the passed channel.
///
/// Interface declaration function.
/// Default behavior. Base class should not be instantiated.
/// Run time error to support debugging.
///
/// \param variableContainer pointer to the variable content bank
/// \param nChannel the interested external channel number
/// \return kTRUE if the current content applies to the configuration
Bool_t AliQnCorrectionsDetectorConfigurationBase::IsSelected(const Float_t *variableContainer, Int_t nChannel) {
  AliFatal(Form("You have reached base member %s. This means you have instantiated a base class or\n" \
      "you are using a non channelized detector configuration but passing a channel number. FIX IT, PLEASE.",
      "AliQnCorrectionsDetectorConfigurationBase::IsSelected()"));
  return kFALSE;
}

