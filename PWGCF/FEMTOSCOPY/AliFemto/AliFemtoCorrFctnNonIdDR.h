///
/// \file AliFemtoCorrFctnNonIdDR.h
///


#ifndef ALIFEMTOCORRFCTNNONIDDR_H
#define ALIFEMTOCORRFCTNNONIDDR_H

class TH1D;
class TNtuple;

#include "AliFemtoCorrFctn.h"

/// \class AliFemtoCorrGctnNonIdDR
/// \brief Correlation function for non-identical particles using k* as the
///    function variable.
///
/// Stores the correlation function separately for positive and negative signs
/// of k* projections into out, side and long directions, enabling the
/// calculations of double ratios
///
class AliFemtoCorrFctnNonIdDR : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnNonIdDR(const char* title, const int nbins, const float QinvLo, const float QinvHi);
  AliFemtoCorrFctnNonIdDR(const AliFemtoCorrFctnNonIdDR& aCorrFctn);
  virtual ~AliFemtoCorrFctnNonIdDR();

  AliFemtoCorrFctnNonIdDR& operator=(const AliFemtoCorrFctnNonIdDR& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  virtual TList* GetOutputList();
  void Write();

  virtual AliFemtoCorrFctn* Clone();
  void FillParticleP(bool);

  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctnNonIdDR(*this); }

protected:
  TH1D *fNumOutP;     ///< Numerator for pair with positive k*out
  TH1D *fNumOutN;     ///< Numerator for pair with negative k*out
  TH1D *fNumSideP;    ///< Numerator for pair with positive k*side
  TH1D *fNumSideN;    ///< Numerator for pair with negative k*side
  TH1D *fNumLongP;    ///< Numerator for pair with positive k*long
  TH1D *fNumLongN;    ///< Numerator for pair with negative k*long

  TH1D *fDenOutP;     ///< Denominator for pair with positive k*out
  TH1D *fDenOutN;     ///< Denominator for pair with negative k*out
  TH1D *fDenSideP;    ///< Denominator for pair with positive k*side
  TH1D *fDenSideN;    ///< Denominator for pair with negative k*side
  TH1D *fDenLongP;    ///< Denominator for pair with positive k*long
  TH1D *fDenLongN;    ///< Denominator for pair with negative k*long

  TH1D* fkTMonitor;   ///< Monitor the kT of pairs in the function
  TNtuple *mNtuple;       //momentum and energy of particles

  bool fParticleP;    //kTRUE to fill nTuple which contain px,py,pz and E of particles

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctnNonIdDR, 2);
  /// \endcond
#endif
};

inline AliFemtoCorrFctn* AliFemtoCorrFctnNonIdDR::Clone()
{
  return new AliFemtoCorrFctnNonIdDR(*this);
}

#endif
