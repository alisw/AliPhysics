////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoTPCInnerCorrFctn - A correlation function that saves the         ///
/// distance at the entrance to the TPC between two tracks as a function     ///
/// of qinv                                                                  ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoTPCInnerCorrFctn_hh
#define AliFemtoTPCInnerCorrFctn_hh

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoTPCInnerCorrFctn : public AliFemtoCorrFctn {
public:
  AliFemtoTPCInnerCorrFctn(const char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoTPCInnerCorrFctn(const AliFemtoTPCInnerCorrFctn& aCorrFctn);
  virtual ~AliFemtoTPCInnerCorrFctn();

  AliFemtoTPCInnerCorrFctn& operator=(const AliFemtoTPCInnerCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);
  void SetRadius(double rad);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
  virtual AliFemtoCorrFctn* Clone() const { return new AliFemtoTPCInnerCorrFctn(*this); }

private:

  TH2D *fDTPCNumerator;        ///< Distance at the entrance to the TPC for real pairs
  TH2D *fDTPCDenominator;      ///< Distance at the entrance to tht TPC for mixed pairs
  TH2D *fRadDNumerator;        ///< Distance at the radius for real pairs
  TH2D *fRadDDenominator;      ///< Distance at the radius for mixed pairs
  Double_t fRadius;            ///< Radius at which to calculate the distance

#ifdef __ROOT__
  ClassDef(AliFemtoTPCInnerCorrFctn, 1)
#endif
};


#endif
