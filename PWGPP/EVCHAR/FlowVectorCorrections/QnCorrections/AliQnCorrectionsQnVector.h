#ifndef ALIQNCORRECTIONS_QNVECTORS_H
#define ALIQNCORRECTIONS_QNVECTORS_H

/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Víctor González, UCM, victor.gonzalez@cern.ch            *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2016                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/

/// \file AliQnCorrectionsQnVector.h
/// \brief Classes that model Q vectors for different harmonics within the Q vector correction framework

#include <TNamed.h>
#include <TMath.h>

/// The maximum external harmonic number the framework currently support for Q vectors
#define MAXHARMONICNUMBERSUPPORTED 15

/// \class AliQnCorrectionsQnVector
/// \brief Class that models and encapsulates a Q vector set
///
/// \author Jaap Onderwaater <jacobus.onderwaater@cern.ch>, GSI
/// \author Ilya Selyuzhenkov <ilya.selyuzhenkov@gmail.com>, GSI
/// \author Víctor González <victor.gonzalez@cern.ch>, UCM
/// \date Jan 27, 2016
class AliQnCorrectionsQnVector : public TNamed {

public:
  /// \enum QnVectorNormalizationMethod
  /// \brief The class of the id of the supported Q vector normalization methods
  ///
  /// Actually it is not a class because the C++ level of implementation.
  /// But full protection will be reached when were possible declaring it
  /// as a class.
  ///
  /// M is the sum of weights.
  enum QnVectorNormalizationMethod {
    QVNORM_noCalibration, ///< \f$ \mbox{Q'} = \mbox{Q}\f$
    QVNORM_QoverSqrtM,    ///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{\sqrt{\mbox{M}}} \f$
    QVNORM_QoverM,        ///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{\mbox{M}} \f$
    QVNORM_QoverQlength   ///< \f$ \mbox{Q'} = \frac{\mbox{Q}}{|\mbox{Q}|} \f$
  };

  AliQnCorrectionsQnVector();
  AliQnCorrectionsQnVector(const char *name, Int_t nNoOfHarmonics, Int_t *harmonicMap = NULL);
  AliQnCorrectionsQnVector(const AliQnCorrectionsQnVector &Q);
  virtual ~AliQnCorrectionsQnVector();

  void ActivateHarmonic(Int_t harmonic);
  Int_t GetNoOfHarmonics() const;
  void GetHarmonicsMap(Int_t *harmonicMap) const;
  Int_t GetFirstHarmonic() const;
  Int_t GetNextHarmonic(Int_t harmonic) const;

  /// Sets the X component for the considered harmonic
  /// \param harmonic the intended harmonic
  /// \param qx the X component for the Q vector
  virtual void SetQx(Int_t harmonic, Float_t qx) { fQnX[harmonic] = qx; }
  /// Sets the Y component for the considered harmonic
  /// \param harmonic the intended harmonic
  /// \param qy the Y component for the Q vector
  virtual void SetQy(Int_t harmonic, Float_t qy) { fQnY[harmonic] = qy; }
  /// Set the good quality flag
  /// \param good kTRUE  if the quality is good
  virtual void SetGood(Bool_t good) { fGoodQuality = good; }


  void Set(AliQnCorrectionsQnVector* Qn, Bool_t changename);

  void Normalize();
  /// Provides the length of the Q vector for the considered harmonic
  /// \param harmonic the intended harmonic
  /// \return the square root of components square sum
  Float_t Length(Int_t harmonic) const { return  TMath::Sqrt(Qx(harmonic)*Qx(harmonic)+Qy(harmonic)*Qy(harmonic));}
  Float_t QxNorm(Int_t harmonic) const;
  Float_t QyNorm(Int_t harmonic) const;

  virtual void Reset();

  /// Gets the Q vector X component for the considered harmonic
  /// \param harmonic the intended harmonic
  /// \return the Q vector X component
  Float_t Qx(Int_t harmonic) const { return fQnX[harmonic]; }
  /// Gets the Q vector Y component for the considered harmonic
  /// \param harmonic the intended harmonic
  /// \return the Q vector Y component
  Float_t Qy(Int_t harmonic) const { return fQnY[harmonic]; }
  /// Get the Qn vector quality flag
  /// \return Qn vector quality flag
  Bool_t IsGoodQuality() const { return fGoodQuality; }

  /// Gets the number of elements that were used for Q vector building
  /// \return number of elements
  Int_t GetN() const { return fN; }
  Double_t EventPlane(Int_t harmonic) const;

  virtual void Print(Option_t *) const;

private:
  /// Assignment operator
  ///
  /// Default implementation to protect against its accidental use.
  /// It will give a compilation error. Don't make it different from
  /// private!
  /// \param Qn the Q vector to assign
  AliQnCorrectionsQnVector& operator= (const AliQnCorrectionsQnVector &Qn);

protected:

  static const Float_t  fMinimumSignificantValue;     ///< the minimum value that will be considered as meaningful for processing
  static const UInt_t   harmonicNumberMask[];         ///< Mask for each external harmonic number

  Float_t fQnX[MAXHARMONICNUMBERSUPPORTED+1];   ///< the Q vector X component for each harmonic
  Float_t fQnY[MAXHARMONICNUMBERSUPPORTED+1];   ///< the Q vector Y component for each harmonic
  Int_t   fHighestHarmonic;                    ///< the highest harmonic number handled
  UInt_t  fHarmonicMask;                       ///< the mask for the supported harmonics
  Bool_t  fGoodQuality;                        ///< Qn vector good quality flag
  Int_t fN;                                    ///< number of elements used for Qn vector building

/// \cond CLASSIMP
  ClassDef(AliQnCorrectionsQnVector, 2);
/// \endcond
};

/// Get the number of the first harmonic used
/// \return the number of the first harmonic handled by the Q vector, -1 if none
inline Int_t AliQnCorrectionsQnVector::GetFirstHarmonic() const {
  for(Int_t h = 1; h < fHighestHarmonic + 1; h++){
    if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
      return h;
    }
  }
  return -1;
}

/// Get the next harmonic to the one passed as parameter
/// \param harmonic number to find the next one
/// \return the number of the next to the passed harmonic, -1 if none
inline Int_t AliQnCorrectionsQnVector::GetNextHarmonic(Int_t harmonic) const {
  for(Int_t h = harmonic+1; h < fHighestHarmonic + 1; h++){
    if ((fHarmonicMask & harmonicNumberMask[h]) == harmonicNumberMask[h]) {
      return h;
    }
  }
  return -1;
}

/// Provides the X component normalized to one of the Q vector for the considered harmonic
/// A check for Q vector length significant value is made. Not passing it returns a zero value.
/// \param harmonic the intended harmonic
/// \return X component of the normalized Q vector
inline Float_t AliQnCorrectionsQnVector::QxNorm(Int_t harmonic) const {
  if (Length(harmonic) < fMinimumSignificantValue) {
    return 0.0;
  }
  else {
      return  Qx(harmonic)/Length(harmonic);
  }
}

/// Provides the Y component normalized to one of the Q vector for the considered harmonic
/// A check for Q vector length significant value is made. Not passing it returns a zero value.
/// \param harmonic the intended harmonic
/// \return Y component of the normalized Q vector
inline Float_t AliQnCorrectionsQnVector::QyNorm(Int_t harmonic) const {
  if (Length(harmonic) < fMinimumSignificantValue) {
    return 0.0;
  }
  else {
    return  Qy(harmonic)/Length(harmonic);
  }
}

#endif /* ALIQNCORRECTIONS_QNVECTORS_H */
