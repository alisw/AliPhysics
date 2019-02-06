///
/// \file AliFemtoBPLCMS3DCorrFctnKK.h
///
/// \class AliFemtoBPLCMS3DCorrFctnKK
/// \brief A class to calculate 3D correlation for pairs of identical
///        particles.
///
/// It also stored the weighted qinv per bin histogram for the coulomb
/// correction.
/// In analysis the function should be first created in a macro, then
/// added to the analysis, and at the end of the macro the procedure to
/// write out histograms should be called.
/// ML March 2015 I deleted dividing in each event !!! for histo qinv
/// as it was done before
///

#ifndef ALIFEMTOBPLCMS3DCORRFCTNKK_H
#define ALIFEMTOBPLCMS3DCORRFCTNKK_H

#include "AliFemtoCorrFctn.h"
//#include "AliFemtoCoulomb.h"
#include "AliFemtoPairCut.h"
//#include "AliFemtoHisto.h"
#include "TH3D.h"
//#include "AliFemtoSmearPair.h"

class AliFemtoBPLCMS3DCorrFctnKK : public AliFemtoCorrFctn {
public:
  AliFemtoBPLCMS3DCorrFctnKK(const char* title, const int& nbins, const float& QLo, const float& QHi);
  AliFemtoBPLCMS3DCorrFctnKK(const AliFemtoBPLCMS3DCorrFctnKK& aCorrFctn);
  virtual ~AliFemtoBPLCMS3DCorrFctnKK();

  AliFemtoBPLCMS3DCorrFctnKK& operator=(const AliFemtoBPLCMS3DCorrFctnKK& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  TH3D* Numerator();
  TH3D* Denominator();
  TH3D* Ratio();
  TH3D* QinvHisto();

  // here are get and set for the range over which the correlation function
  // is normalized (in Qinv).  The range is set to 0.15..0.18 in the constuctor
  // by default, but the Set's below override this
  void SetNormRangeLo(float qLo);
  void SetNormRangeHi(float qHi);
  float GetNormRangeLo() const;
  float GetNormRangeHi() const;

  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoBPLCMS3DCorrFctnKK(*this); }

  //  void SetCoulombCorrection(AliFemtoCoulomb* Correction);

  void SetUseRPSelection(unsigned short aRPSel);

  //  void SetSmearPair(AliFemtoSmearPair*);
  void SetRout(double guess);
  void SetRside(double guess);
  void SetRlong(double guess);
  void SetLambda(double guess);

private:
  /*   // here are a whole bunch of histos that get filled if we do resolution correction */
  /*   TH3D* fIDNumHisto;        // true pairs numerator   */
  /*   TH3D* fIDDenHisto;        // true pairs denominator */
  /*   TH3D* fIDRatHisto;        // true pairs ratio       */
  /*   //  */
  /*   TH3D* fSMNumHisto;        // mixed pairs numerator   */
  /*   TH3D* fSMDenHisto;     // mixed pairs denominator */
  /*   TH3D* fSMRatHisto;     // mixed pairs ratio       */
  /*   // */
  /*   TH3D* fCorrectionHisto;   // correction histogram */
  /*   TH3D* fCorrCFHisto;       // Corrected CF */

  TH3D* fNumerator;         ///< numerator
  TH3D* fDenominator;       ///< denominator
  //  TH3D* fUncorrectedDenominator;
  TH3D* fRatio;             ///< ratio - the correlation function
  TH3D* fQinvHisto;         ///< Qinv weights

  // for resolution correction
  //  AliFemtoSmearPair* fSmearPair; //!
  double fLambda;           ///< lambda for smearing correction
  double fRout2;            ///< Rout for smearing correction
  double fRside2;           ///< Rside for smearing correction
  double fRlong2;           ///< Rlong for smearing correction

  // upper and lower bounds of Qinv region where to do normalization
  float fQinvNormLo;        ///< Lower bound of Qinv normalization range
  float fQinvNormHi;        ///< Upper bound of Qinv normalization range

  // and here are the number of pairs in that region...
  unsigned long int fNumRealsNorm; ///< pairs in numerator in Qinv normalization range
  unsigned long int fNumMixedNorm; ///< pairs in denominator in Qinv normalization range

protected:
  unsigned short fUseRPSelection;  ///< The pair cut uses RP selection

  //  AliFemtoCoulomb* fCorrection; //!


#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoBPLCMS3DCorrFctnKK, 1);
  /// \endcond
#endif
};

inline  TH3D* AliFemtoBPLCMS3DCorrFctnKK::Numerator()
{
  return fNumerator;
}
inline  TH3D* AliFemtoBPLCMS3DCorrFctnKK::Denominator()
{
  return fDenominator;
}
//inline  TH3D* AliFemtoBPLCMS3DCorrFctnKK::UncorrectedDenominator(){return fUncorrectedDenominator;}
inline  TH3D* AliFemtoBPLCMS3DCorrFctnKK::Ratio()
{
  return fRatio;
}
inline  TH3D* AliFemtoBPLCMS3DCorrFctnKK::QinvHisto()
{
  return fQinvHisto;
}
inline  void AliFemtoBPLCMS3DCorrFctnKK::SetNormRangeLo(float qLo)
{
  fQinvNormLo = qLo;
}
inline  void AliFemtoBPLCMS3DCorrFctnKK::SetNormRangeHi(float qHi)
{
  fQinvNormHi = qHi;
}
inline  float AliFemtoBPLCMS3DCorrFctnKK::GetNormRangeLo() const
{
  return fQinvNormLo;
}
inline  float AliFemtoBPLCMS3DCorrFctnKK::GetNormRangeHi() const
{
  return fQinvNormHi;
}
//inline  void AliFemtoBPLCMS3DCorrFctnKK::SetCoulombCorrection(AliFemtoCoulomb* Correction){fCorrection = Correction;}
//inline  void AliFemtoBPLCMS3DCorrFctnKK::SetSmearPair(AliFemtoSmearPair* sp){fSmearPair = sp;}

inline  void AliFemtoBPLCMS3DCorrFctnKK::SetRout(double r)
{
  fRout2 = r * r;
}
inline  void AliFemtoBPLCMS3DCorrFctnKK::SetRside(double r)
{
  fRside2 = r * r;
}
inline  void AliFemtoBPLCMS3DCorrFctnKK::SetRlong(double r)
{
  fRlong2 = r * r;
}
inline  void AliFemtoBPLCMS3DCorrFctnKK::SetLambda(double l)
{
  fLambda = l;
}

#endif
