///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoModelCorrFctn3DSpherical: a class to calculate 3D correlation //
// for pairs of identical particles, binned in spherical coordinates     //
// (q_inv, phi, cos(theta))                                              //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOMODELCORRFCTN3DLCMSSPHERICAL_H
#define ALIFEMTOMODELCORRFCTN3DLCMSSPHERICAL_H

#include "AliFemtoModelCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoModelCorrFctn3DLCMSSpherical : public AliFemtoModelCorrFctn {
public:
  AliFemtoModelCorrFctn3DLCMSSpherical(char* title, 
			      const int& nqbins, const float& QLo, const float& QHi,
			      const int& nphibins, const int& ncthetabins);
  AliFemtoModelCorrFctn3DLCMSSpherical(const AliFemtoModelCorrFctn3DLCMSSpherical& aCorrFctn);
  virtual ~AliFemtoModelCorrFctn3DLCMSSpherical();

  AliFemtoModelCorrFctn3DLCMSSpherical& operator=(const AliFemtoModelCorrFctn3DLCMSSpherical& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  void WriteOutHistos();
  virtual TList* GetOutputList();

  void SetSpecificPairCut(AliFemtoPairCut* aCut);

private:
  // here are a whole bunch of histos that get filled if we do resolution correction
  TH3D* fTrueNumeratorSph;         // numerator
  TH3D* fFakeNumeratorSph;         // numerator
  TH3D* fDenominatorSph;       // denominator

  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis

#ifdef __ROOT__
  ClassDef(AliFemtoModelCorrFctn3DLCMSSpherical, 1)
#endif
};

inline  void AliFemtoModelCorrFctn3DLCMSSpherical::SetSpecificPairCut(AliFemtoPairCut* pc){fPairCut=pc;}

#endif

