///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DSpherical: a class to calculate 3D correlation      //
// for pairs of identical particles, binned in spherical coordinates     //
// (q_inv, phi, cos(theta))
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTN3DSPHERICAL_H
#define ALIFEMTOCORRFCTN3DSPHERICAL_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoCorrFctn3DSpherical : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctn3DSpherical(char* title, 
			      const int& nqbins, const float& QLo, const float& QHi,
			      const int& nphibins, const int& ncthetabins);
  AliFemtoCorrFctn3DSpherical(const AliFemtoCorrFctn3DSpherical& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DSpherical();

  AliFemtoCorrFctn3DSpherical& operator=(const AliFemtoCorrFctn3DSpherical& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  void WriteOutHistos();
  virtual TList* GetOutputList();

  void SetSpecificPairCut(AliFemtoPairCut* aCut);

private:
  // here are a whole bunch of histos that get filled if we do resolution correction
  TH3D* fNumerator;         // numerator
  TH3D* fDenominator;       // denominator

  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctn3DSpherical, 1)
#endif
};

inline  void AliFemtoCorrFctn3DSpherical::SetSpecificPairCut(AliFemtoPairCut* pc){fPairCut=pc;}

#endif

