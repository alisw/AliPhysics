///
/// \file AliFemtoAvgSepCorrFctn.h
/// \author M. Janik, L. Graczykowski, Warsaw University of Technology
///
/// \class AliFemtoAvgSepCorrFctn
/// \brief An average entrance separation correlation function
///

#ifndef ALIFEMTOAVGSEPCORRFCTN_H
#define ALIFEMTOAVGSEPCORRFCTN_H

#include "TH1D.h"
#include "AliFemtoCorrFctn.h"

class AliFemtoAvgSepCorrFctn : public AliFemtoCorrFctn {
public:
  enum PairType {kTracks = 0, kTrackV0 = 1, kV0s = 2};
  typedef enum PairType AliFemtoPairType;

  AliFemtoAvgSepCorrFctn(char *title, const int &nbins, const float &Low, const float &High);
  AliFemtoAvgSepCorrFctn(const AliFemtoAvgSepCorrFctn &aCorrFctn);
  virtual ~AliFemtoAvgSepCorrFctn();

  AliFemtoAvgSepCorrFctn &operator=(const AliFemtoAvgSepCorrFctn &aCorrFctn);

  virtual AliFemtoString Report();
  virtual void AddRealPair(AliFemtoPair *aPair);
  virtual void AddMixedPair(AliFemtoPair *aPair);
  virtual void Finish();

  TH1D *Numerator();
  TH1D *Denominator();
  TH1D *Ratio();

  virtual TList *GetOutputList();
  void Write();
  void SetPairType(AliFemtoPairType pairtype);

private:
  //2 tracks
  TH1D *fNumerator;          // numerator - real pairs
  TH1D *fDenominator;        // denominator - mixed pairs

  //track + V0
  TH1D *fNumeratorPos;          // numerator - real pairs
  TH1D *fDenominatorPos;        // denominator - mixed pairs
  TH1D *fNumeratorNeg;          // numerator - real pairs
  TH1D *fDenominatorNeg;        // denominator - mixed pairs

  //2 V0s
  TH1D *fNumeratorPosPos;          // numerator - real pairs
  TH1D *fDenominatorPosPos;        // denominator - mixed pairs
  TH1D *fNumeratorPosNeg;          // numerator - real pairs
  TH1D *fDenominatorPosNeg;        // denominator - mixed pairs
  TH1D *fNumeratorNegPos;          // numerator - real pairs
  TH1D *fDenominatorNegPos;        // denominator - mixed pairs
  TH1D *fNumeratorNegNeg;          // numerator - real pairs
  TH1D *fDenominatorNegNeg;        // denominator - mixed pairs

  TH1D *fRatio;              // ratio - correlation function
  AliFemtoPairType fPairType;


#ifdef __ROOT__
  ClassDef(AliFemtoAvgSepCorrFctn, 1)
#endif
};

inline  TH1D *AliFemtoAvgSepCorrFctn::Numerator()
{
  return fNumerator;
}
inline  TH1D *AliFemtoAvgSepCorrFctn::Denominator()
{
  return fDenominator;
}
inline  TH1D *AliFemtoAvgSepCorrFctn::Ratio()
{
  return fRatio;
}


#endif

