////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoShareQualityCorrFctn - A correlation function that saves the     ///
/// amount of sharing and splitting hits per pair as a function of qinv      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef AliFemtoShareQualityCorrFctn_hh
#define AliFemtoShareQualityCorrFctn_hh

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoShareQualityCorrFctn : public AliFemtoCorrFctn {
public:
  AliFemtoShareQualityCorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoShareQualityCorrFctn(const AliFemtoShareQualityCorrFctn& aCorrFctn);
  virtual ~AliFemtoShareQualityCorrFctn();

  AliFemtoShareQualityCorrFctn& operator=(const AliFemtoShareQualityCorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fShareNumerator;        // Share fraction for real pairs
  TH2D *fShareDenominator;      // share fraction for mixed pairs
 
  TH2D *fQualityNumerator;      // quality for real pairs
  TH2D *fQualityDenominator;    // quality for mixed pairs 

  TH2D *fTPCSepNumerator;       // TPCSep for real pairs
  TH2D *fTPCSepDenominator;     // TPCSep for mixed pairs 

#ifdef __ROOT__
  ClassDef(AliFemtoShareQualityCorrFctn, 1)
#endif
};


#endif

