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

/// \file AliQnCorrectionsQnVectorBuild.cxx
/// \brief Implementation of Q vector while building class

#include <Riostream.h>

#include "AliQnCorrectionsQnVectorBuild.h"
#include "AliLog.h"

using std::cout;
using std::endl;

/// \cond CLASSIMP
ClassImp(AliQnCorrectionsQnVectorBuild);
/// \endcond

/// Default constructor
AliQnCorrectionsQnVectorBuild::AliQnCorrectionsQnVectorBuild() : AliQnCorrectionsQnVector() {

  fSumW = 0.0;
}

/// Normal constructor
///
/// Relays on its parent for almost everything
///
/// \param name the name of the Qn vector. Identifies its origin
/// \param nNoOfHarmonics the desired number of harmonics
/// \param harmonicMap ordered array with the external number of the harmonics
AliQnCorrectionsQnVectorBuild::AliQnCorrectionsQnVectorBuild(const char *name, Int_t nNoOfHarmonics, Int_t *harmonicMap) :
    AliQnCorrectionsQnVector(name, nNoOfHarmonics, harmonicMap) {

  fSumW = 0.0;
}

/// Copy constructor from a Q vector
/// \param Qn the Q vector build object to copy after construction
AliQnCorrectionsQnVectorBuild::AliQnCorrectionsQnVectorBuild(const AliQnCorrectionsQnVector &Qn) :
    AliQnCorrectionsQnVector(Qn) {

  fSumW = 0.0;
}

/// Copy constructor
/// \param Qn the Q vector build object to copy after construction
AliQnCorrectionsQnVectorBuild::AliQnCorrectionsQnVectorBuild(const AliQnCorrectionsQnVectorBuild &Qn) :
    AliQnCorrectionsQnVector(Qn) {

  fSumW = Qn.fSumW;
}

/// Default destructor
AliQnCorrectionsQnVectorBuild::~AliQnCorrectionsQnVectorBuild() {

}

/// Sets the X component for the considered harmonic
///
/// It should not be used. Runtime error indication.
void AliQnCorrectionsQnVectorBuild::SetQx(Int_t, Float_t) {

  AliFatal("You are using a forbidden function for a build Q vector");
}

/// Sets the Y component for the considered harmonic
///
/// It should not be used. Runtime error indication.
void AliQnCorrectionsQnVectorBuild::SetQy(Int_t, Float_t) {

  AliFatal("You are using a forbidden function for a build Q vector");
}

/// Copy member function
///
/// The passed Q vector is copied within the current object
/// \param Qn pointer to the Q vector to be copied
void AliQnCorrectionsQnVectorBuild::Set(AliQnCorrectionsQnVectorBuild* Qn) {

  /* the name is not copied from building Qn vectors */
  AliQnCorrectionsQnVector::Set(Qn,kFALSE);
  fSumW = Qn->fSumW;
}

/// Adds a build Q vector
///
/// Warning: the possibility of a different set of harmonics for both
/// build Q vectors is currently not considered
/// \param Qn the build Q vector to add
void AliQnCorrectionsQnVectorBuild::Add(AliQnCorrectionsQnVectorBuild* Qn) {

  for(Int_t h = 1; h < fHighestHarmonic + 1; h++){
    if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
      fQnX[h] += Qn->Qx(h);
      fQnY[h] += Qn->Qy(h);
    }
  }
  fSumW += Qn->GetSumOfWeights();
  fN += Qn->GetN();
}

/// Normalizes the build Q vector for the whole harmonics set
///
/// Normalizes the build Q vector as \f$ Qn = \frac{Qn}{M} \f$.
/// A check for significant value is made. Not passing it
/// does set the Q vector quality as bad
void AliQnCorrectionsQnVectorBuild::NormalizeQoverM() {

  if(fSumW < fMinimumSignificantValue) {
    SetGood(kFALSE);
  }
  else {
    for(Int_t h = 1; h < fHighestHarmonic + 1; h++){
      if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
        fQnX[h] = (fQnX[h] / fSumW);
        fQnY[h] = (fQnY[h] / fSumW);
      }
    }
  }
}

/// Normalizes the build Q vector for the whole harmonics set
///
/// Normalizes the build Q vector as \f$ Qn = \frac{Qn}{\sqrt{M}} \f$.
/// A check for significant value is made. Not passing it
/// does set the Q vector quality as bad
void AliQnCorrectionsQnVectorBuild::NormalizeQoverSquareRootOfM() {

  if(fSumW < fMinimumSignificantValue) {
    SetGood(kFALSE);
  }
  else {
    for(Int_t h = 1; h < fHighestHarmonic + 1; h++){
      if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
        fQnX[h] += fQnX[h] / TMath::Sqrt(fSumW);
        fQnY[h] += fQnY[h] / TMath::Sqrt(fSumW);
      }
    }
  }
}

/// Resets the Q vector values without touching the structure
void AliQnCorrectionsQnVectorBuild::Reset() {

  AliQnCorrectionsQnVector::Reset();
  fSumW = 0.0;
}

/// Print the Qn vector in a readable shape
///
void AliQnCorrectionsQnVectorBuild::Print(Option_t *) const {
  cout <<"OBJ: building Qn vector" << "\t" << "N: " << fN << "\t" << "Sum w: " << fSumW << "\t"
      << "quality: " << ((fGoodQuality) ? "good" : "bad") << endl;
  Int_t harmonic = GetFirstHarmonic();
  while (harmonic != -1) {
    cout << "\t" << "\t" << "harmonic " << harmonic << "\t" << "QX: " << Qx(harmonic) << "\t" << "QY: " << Qy(harmonic) << endl;
    harmonic = GetNextHarmonic(harmonic);
  }
}


