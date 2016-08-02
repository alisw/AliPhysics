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
/// \file AliQnCorrectionsCutsBase.cxx
/// \brief Implementation of the base class for the classes that model the cuts support for the Q vector correction framework

#include "AliQnCorrectionsCutsBase.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCutsBase);
/// \endcond

const Int_t AliQnCorrectionsCutsBase::nHighestBitNumberSupported = 31;

/// Default constructor
AliQnCorrectionsCutsBase::AliQnCorrectionsCutsBase() : TObject() {
  fVarId = -1;
}

/// Copy constructor
/// \param cut the cut object to clone
AliQnCorrectionsCutsBase::AliQnCorrectionsCutsBase(const AliQnCorrectionsCutsBase &cut) :
    TObject(cut) {
  fVarId = cut.fVarId;
}

/// Normal constructor. Stores the affected variable Id
/// \param varId External Id for the affected variable
AliQnCorrectionsCutsBase::AliQnCorrectionsCutsBase(Int_t varId) :
    TObject() {
  fVarId = varId;
}

/// Default destructor. Does nothing
AliQnCorrectionsCutsBase::~AliQnCorrectionsCutsBase() {
}


