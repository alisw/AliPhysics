////////////////////////////////////////////////////////////////////////////////
///                                                                          ///
/// AliFemtoChi2CorrFctn - A correlation function that saves the correlation ///
/// function as a function of single track quality (chi2/ndof) for its and   ///
/// tpc                                                                      ///
/// Authors: Adam Kisiel kisiel@mps.ohio-state.edu                           ///
///                                                                          ///
////////////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOCHI2CORRFCTN_H
#define ALIFEMTOCHI2CORRFCTN_H

#include "TH1D.h"
#include "TH2D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoChi2CorrFctn : public AliFemtoCorrFctn {
public:
  AliFemtoChi2CorrFctn(char* title, const int& nbins, const float& QinvLo, const float& QinvHi);
  AliFemtoChi2CorrFctn(const AliFemtoChi2CorrFctn& aCorrFctn);
  virtual ~AliFemtoChi2CorrFctn();

  AliFemtoChi2CorrFctn& operator=(const AliFemtoChi2CorrFctn& aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair* aPair);
  virtual void AddMixedPair(AliFemtoPair* aPair);

  virtual void Finish();

  void WriteHistos();
  virtual TList* GetOutputList();
private:
  
  TH2D *fChi2ITSSUMNumerator;        // Numerator as a function of ITS quality sum for the pair
  TH2D *fChi2ITSSUMDenominator;      // Denominator as a function of ITS quality sum for the pair
 
  TH2D *fChi2TPCSUMNumerator;        // Numerator as a function of TPC quality sum for the pair
  TH2D *fChi2TPCSUMDenominator;      // Denominator as a function of TPC quality sum for the pair

  TH2D *fChi2ITSONENumerator;        // Numerator as a function of ITS quality for the worse track
  TH2D *fChi2ITSONEDenominator;      // Denominator as a function of ITS quality for the worse track
 
  TH2D *fChi2TPCONENumerator;        // Numerator as a function of TPC quality for the worse track
  TH2D *fChi2TPCONEDenominator;      // Denominator as a function of TPC quality for the worse track

  TH2D *fSigmaToVertexNumerator;     // Numerator as a function of sigma to vertex
  TH2D *fSigmaToVertexDenominator;   // Numerator as a function of sigma to vertex

#ifdef __ROOT__
  ClassDef(AliFemtoChi2CorrFctn, 1)
#endif
};


#endif

