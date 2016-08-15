#ifndef ALIQNCORRECTIONS_QNVECTORSBUILD_H
#define ALIQNCORRECTIONS_QNVECTORSBUILD_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsQnVectorBuild.h
/// \brief Class that models and encapsulates a Q vector set while building it within the Q vector correction framework

#include "AliQnCorrectionsQnVector.h"

/// \class AliQnCorrectionsQnVectorBuild
/// \brief Class that models and encapsulates a Q vector set while building it
///
/// When the Q vector is being built it needs extra support.
/// This class provides such extra support.
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Feb 02, 2016
class AliQnCorrectionsQnVectorBuild : public AliQnCorrectionsQnVector {

public:
  AliQnCorrectionsQnVectorBuild();
  AliQnCorrectionsQnVectorBuild(const char *name, Int_t nNoOfHarmonics, Int_t *harmonicMap = NULL);
  AliQnCorrectionsQnVectorBuild(const AliQnCorrectionsQnVector &Qn);
  AliQnCorrectionsQnVectorBuild(const AliQnCorrectionsQnVectorBuild &Qn);
  virtual ~AliQnCorrectionsQnVectorBuild();

  virtual void SetQx(Int_t harmonic, Float_t qx);
  virtual void SetQy(Int_t harmonic, Float_t qy);

  void Set(AliQnCorrectionsQnVectorBuild* Qn);

  void Add(AliQnCorrectionsQnVectorBuild* qvec);
  void Add(Double_t phi, Double_t weight = 1.0);

  /// Check the quality of the constructed Qn vector
  /// Current criteria is number of contributors should be at least one.
  /// If so happen, sets the good quality flag.
  void CheckQuality() { fGoodQuality = ((0 < fN) ? kTRUE : kFALSE); }
  void Normalize(QnVectorNormalizationMethod method);

  void NormalizeQoverM();
  void NormalizeQoverSquareRootOfM();

  virtual void Reset();

  /// Gets the sum of weights.
  /// \return sum of weights
  Float_t GetSumOfWeights() const { return fSumW; }

  virtual void Print(Option_t *) const;

private:
  /// Assignment operator
  ///
  /// Default implementation to protect against its accidental use.
  /// It will give a compilation error. Don't make it different from
  /// private!
  /// \param Qn the Q vector to assign
  AliQnCorrectionsQnVectorBuild& operator= (const AliQnCorrectionsQnVectorBuild &Qn);

protected:

  Float_t fSumW;   ///< the sum of weights

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsQnVectorBuild, 1);
/// \endcond
};

/// Adds a contribution to the build Q vector
/// A check for weight significant value is made. Not passing it ignores the contribution.
/// The process of incorporating contributions takes into account the harmonic multiplier
/// \param phi azimuthal angle contribution
/// \param weight the weight of the contribution
inline void AliQnCorrectionsQnVectorBuild::Add(Double_t phi, Double_t weight) {

  if (weight < fMinimumSignificantValue) return;
  for(Int_t h = 1; h < fHighestHarmonic + 1; h++){
    if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
      fQnX[h] += (weight * TMath::Cos(h*fHarmonicMultiplier*phi));
      fQnY[h] += (weight * TMath::Sin(h*fHarmonicMultiplier*phi));
    }
  }
  fSumW += weight;
  fN += 1;
}


/// Calibrates the Q vector according to the method passed
/// \param method the method of calibration
inline void AliQnCorrectionsQnVectorBuild::Normalize(QnVectorNormalizationMethod method) {
  switch (method) {
  case QVNORM_noCalibration:
    break;
  case QVNORM_QoverSqrtM:
    NormalizeQoverSquareRootOfM();
    break;
  case QVNORM_QoverM:
    NormalizeQoverM();
    break;
  case QVNORM_QoverQlength:
    AliQnCorrectionsQnVector::Normalize();
    break;
  }
}


#endif /* ALIQNCORRECTIONS_QNVECTORSBUILD_H */
