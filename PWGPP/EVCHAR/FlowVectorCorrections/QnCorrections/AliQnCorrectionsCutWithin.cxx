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
/// \file AliQnCorrectionsCutWithin.cxx
/// \brief Implementation of the within range cut class support for the Q vector correction framework

#include "AliQnCorrectionsCutWithin.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCutWithin);
/// \endcond

/// Default constructor
AliQnCorrectionsCutWithin::AliQnCorrectionsCutWithin() :  AliQnCorrectionsCutsBase() {
  fMinThreshold = 1e10;
  fMaxThreshold = -1e10;
}

/// Copy constructor
/// \param cut the cut object to be cloned
AliQnCorrectionsCutWithin::AliQnCorrectionsCutWithin(const AliQnCorrectionsCutWithin &cut) :
    AliQnCorrectionsCutsBase(cut) {
  fMinThreshold = cut.fMinThreshold;
  fMaxThreshold = cut.fMaxThreshold;
}

/// Normal constructor
/// \param varId external Id for the affected variable
/// \param min the value the variable content should be above
/// \param max the value the variable content should be below
AliQnCorrectionsCutWithin::AliQnCorrectionsCutWithin(Int_t varId, Float_t min, Float_t max) :
    AliQnCorrectionsCutsBase(varId) {
  fMinThreshold = min;
  fMaxThreshold = max;
}

/// Default destructor. Does nothing
AliQnCorrectionsCutWithin::~AliQnCorrectionsCutWithin() {
}


