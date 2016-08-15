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
/// \file AliQnCorrectionsCutSetBit.cxx
/// \brief Implementation of the classes that model the cuts support
/// \brief Implementation of the bit setting cut class support for the Q vector correction framework

#include "AliQnCorrectionsCutSetBit.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsCutSetBit);
/// \endcond

/// Default constructor
AliQnCorrectionsCutSetBit::AliQnCorrectionsCutSetBit() : AliQnCorrectionsCutsBase() {
  fBitMask = 0x00000000;
  fExpectedResult = 0xFFFFFFFF;
}

/// Copy constructor
/// \param cut the cut object to be cloned
AliQnCorrectionsCutSetBit::AliQnCorrectionsCutSetBit(const AliQnCorrectionsCutSetBit &cut) :
    AliQnCorrectionsCutsBase(cut) {
  fBitMask = cut.fBitMask;
  fExpectedResult = cut.fExpectedResult;
}

/// Normal constructor
/// \param varId external Id for the affected variable
/// \param bitNo the bit on the variable content to test (from 0 to 31)
/// \param set (kFALSE)kTRUE for cut on the bit (un)set
AliQnCorrectionsCutSetBit::AliQnCorrectionsCutSetBit(Int_t varId, Int_t bitNo, Bool_t set) :
    AliQnCorrectionsCutsBase(varId) {
  if (nHighestBitNumberSupported < bitNo) {
    AliFatal(Form("You requested a cut on bit %d but the highest bit number supported by the framework is currently %d",
        bitNo, nHighestBitNumberSupported));
  }
  fBitMask = 0x00000001 << bitNo;
  fExpectedResult = (set ? fBitMask : 0x00000000);
}

/// Default destructor. Does nothing
AliQnCorrectionsCutSetBit::~AliQnCorrectionsCutSetBit() {
}


