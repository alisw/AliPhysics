///
/// \file AliFmeto/AliFemtoCorrFctn3DSpherical.h
///


#ifndef ALIFEMTOCORRFCTN3DSPHERICAL_H
#define ALIFEMTOCORRFCTN3DSPHERICAL_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

/// \class AliFemtoCorrFctn3DSpherical
/// \brief 3D correlation for pairs of identical particles, binned
///        in spherical coordinates (q_inv, phi, cos(theta))
///
class AliFemtoCorrFctn3DSpherical : public AliFemtoCorrFctn {
public:

  AliFemtoCorrFctn3DSpherical(const char* title,
                              const int nqbins,
                              const float QLo,
                              const float QHi,
                              const int nphibins,
                              const int ncthetabins);

  AliFemtoCorrFctn3DSpherical(const AliFemtoCorrFctn3DSpherical& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DSpherical();

  AliFemtoCorrFctn3DSpherical& operator=(const AliFemtoCorrFctn3DSpherical& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctn3DSpherical(*this); }


private:

  /// numerator
  TH3D* fNumerator;

  /// denominator
  TH3D* fDenominator;

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassDef(AliFemtoCorrFctn3DSpherical, 1);
  /// \endcond
#endif
};


#endif
