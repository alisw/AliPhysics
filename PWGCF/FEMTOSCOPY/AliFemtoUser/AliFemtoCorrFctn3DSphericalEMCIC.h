
///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCorrFctn3DSpherical: a class to calculate 3D correlation      //
// for pairs of identical particles, binned in spherical coordinates     //
// (q_inv, phi, cos(theta))
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTN3DSPHERICALEMCIC_H
#define ALIFEMTOCORRFCTN3DSPHERICALEMCIC_H

#include "AliFemtoCorrFctn.h"
#include "AliFemtoPairCut.h"
#include "TH3D.h"

class AliFemtoCorrFctn3DSphericalEMCIC : public AliFemtoCorrFctn{
public:
  AliFemtoCorrFctn3DSphericalEMCIC(const char* title,
			      const int& nqbins, const float& QLo, const float& QHi,
			      const int& nphibins, const int& ncthetabins);
  AliFemtoCorrFctn3DSphericalEMCIC(const AliFemtoCorrFctn3DSphericalEMCIC& aCorrFctn);
  virtual ~AliFemtoCorrFctn3DSphericalEMCIC();

  AliFemtoCorrFctn3DSphericalEMCIC& operator=(const AliFemtoCorrFctn3DSphericalEMCIC& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair( AliFemtoPair* aPair);
  virtual void AddMixedPair( AliFemtoPair* aPair);

  virtual void Finish();

  void WriteOutHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoCorrFctn3DSphericalEMCIC(*this); }

  void SetSpecificPairCut(AliFemtoPairCut* aCut);

 private:

  TH3D* fNumerator;         // numerator
  TH3D* fDenominator;       // denominator
  //EMCIC histograms:
  /*TH3D* fEnergyTotalReal;       // E1+E2 from real pairs
  TH3D* fEnergyMultReal;        // E1*E2
  TH3D* fPzMultReal;            // Pz1*Pz2
  TH3D* fPtMultReal;            // Pt1*Pt2  */
  TH3D* fEnergyTotalMix;       // E1+E2 from mixed pairs
  TH3D* fEnergyMultMix;        // E1*E2
  TH3D* fPzMultMix;            // Pz1*Pz2
  TH3D* fPtMultMix;            // Pt1*Pt2
  AliFemtoPairCut* fPairCut;    //! this is a PairCut specific to THIS CorrFctn, not the Analysis

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctn3DSphericalEMCIC, 1)
#endif
};

inline  void AliFemtoCorrFctn3DSphericalEMCIC::SetSpecificPairCut(AliFemtoPairCut* pc){fPairCut=pc;}

#endif
