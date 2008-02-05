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
  AliFemtoTPCInnerCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoTPCInnerCorrFctn(const AliFemtoTPCInnerCorrFctn& aCorrFctn);
  virtual ~AliFemtoTPCInnerCorrFctn();

  AliFemtoTPCInnerCorrFctn& operator=(const AliFemtoTPCInnerCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fDTPCNumerator;        // Distance at the entrance to the TPC for real pairs
  TH2D *fDTPCDenominator;      // Distance at the entrance to tht TPC for mixed pairs
 
#ifdef __ROOT__
  ClassDef(AliFemtoTPCInnerCorrFctn, 1)
#endif
};


#endif

