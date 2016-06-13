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

/// \file AliQnCorrectionsCorrectionOnQvector.cxx
/// \brief Correction steps on Qn vectors base class implementation

#include "AliQnCorrectionsCorrectionOnQvector.h"
#include "AliQnCorrectionsQnVector.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCorrectionOnQvector);
/// \endcond


/// Default constructor
AliQnCorrectionsCorrectionOnQvector::AliQnCorrectionsCorrectionOnQvector() :
    AliQnCorrectionsCorrectionStepBase() {

  fCorrectedQnVector = NULL;
}

/// Normal constructor
/// \param name of the correction step
/// \param key the associated ordering key
AliQnCorrectionsCorrectionOnQvector::AliQnCorrectionsCorrectionOnQvector(const char *name, const char *key) :
    AliQnCorrectionsCorrectionStepBase(name, key) {

  fCorrectedQnVector = NULL;
}

/// Default destructor
AliQnCorrectionsCorrectionOnQvector::~AliQnCorrectionsCorrectionOnQvector() {

  if (fCorrectedQnVector != NULL)
    delete fCorrectedQnVector;
}

/// Include the new corrected Qn vector into the passed list
///
/// Adds the Qn vector to the passed list
/// if the correction step is in correction states.
/// \param list list where the corrected Qn vector should be added
void AliQnCorrectionsCorrectionOnQvector::IncludeCorrectedQnVector(TList *list) {

  switch (fState) {
  case QCORRSTEP_calibration:
    /* collect the data needed to further produce correction parameters */
    break;
  case QCORRSTEP_applyCollect:
    /* collect the data needed to further produce correction parameters */
    /* and proceed to ... */
  case QCORRSTEP_apply: /* apply the correction */
    list->Add(fCorrectedQnVector);
    break;
  }
}


