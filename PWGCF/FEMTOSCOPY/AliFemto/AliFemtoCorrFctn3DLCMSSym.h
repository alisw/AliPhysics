///
/// \file AliFemtoCorrFctn3DLCMSSym.h
///

#ifndef ALIFEMTOCORRFCTN3DLCMS_H
#define ALIFEMTOCORRFCTN3DLCMS_H

// forward declare classes
class TH3F;
class AliFemtoPairCut;

// include headers
#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrFctn3DLCMSSym
/// \brief A class to calculate 3D correlation functions for pairs of identical
///        particles vs. Bertsh-Pratt coordinates.
///
///
class AliFemtoCorrFctn3DLCMSSym : public AliFemtoCorrFctn {
public:

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of , and q
  ///
  AliFemtoCorrFctn3DLCMSSym(const char* title, const int nbins, const float QHi);
  AliFemtoCorrFctn3DLCMSSym(const char* title, const int nbins, const float QHi, bool enable_errs);

  /// Copy Constructor
  AliFemtoCorrFctn3DLCMSSym(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn);

  /// Deletes histograms
  virtual ~AliFemtoCorrFctn3DLCMSSym();

  AliFemtoCorrFctn3DLCMSSym& operator=(const AliFemtoCorrFctn3DLCMSSym& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  TH3F* Numerator();
  TH3F* Denominator();
  TH3F* NumeratorW();   ///< Return numerator weighed by qinv
  TH3F* DenominatorW(); ///< Return denominator weighed by qinv

  void WriteOutHistos();
  virtual TList* GetOutputList();

  void SetUseLCMS(int);
  int  GetUseLCMS();
  virtual AliFemtoCorrFctn* Clone() const;

private:

  TH3F* fNumerator;     ///< Numerator
  TH3F* fDenominator;   ///< Denominator
  TH3F* fNumeratorW;    ///< Qinv-Weighted numerator
  TH3F* fDenominatorW;  ///< Qinv-Weighted denominator

  int    fUseLCMS;      ///< 0 - Use PRF, 1 - Use LCMS

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctn3DLCMSSym, 1);
  /// \endcond
#endif
};

inline AliFemtoCorrFctn* AliFemtoCorrFctn3DLCMSSym::Clone() const
{
  return new AliFemtoCorrFctn3DLCMSSym(*this);
}

inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Numerator()
{
  return fNumerator;
}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::Denominator()
{
  return fDenominator;
}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::NumeratorW()
{
  return fNumeratorW;
}
inline  TH3F* AliFemtoCorrFctn3DLCMSSym::DenominatorW()
{
  return fDenominatorW;
}
#endif
