////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCorrFctnTPCNcls - A correlation function that saves the correlation//
// function as a function of number of TPC clusters of the track              //
//                                                                            //
// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCORRFCTNTPCNCLS_H
#define ALIFEMTOCORRFCTNTPCNCLS_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoCorrFctnTPCNcls : public AliFemtoCorrFctn {
public:
  AliFemtoCorrFctnTPCNcls(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoCorrFctnTPCNcls(const AliFemtoCorrFctnTPCNcls& aCorrFctn);
  virtual ~AliFemtoCorrFctnTPCNcls();

  AliFemtoCorrFctnTPCNcls& operator=(const AliFemtoCorrFctnTPCNcls& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fNclsTPCMinNumerator;        // Numerator as a function of lower TPC Ncls of the pair 
  TH2D *fNclsTPCMinDenominator;      // Denominator as a function of lower TPC Ncls of the pair

#ifdef __ROOT__
  ClassDef(AliFemtoCorrFctnTPCNcls, 1)
#endif
};


#endif

