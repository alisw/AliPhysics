// @(#) $Id$

#ifndef ALIL3HOUGHTEST_H
#define ALIL3HOUGHTEST_H

#include "AliHLTRootTypes.h"

struct AliHLTSimData {
  Int_t fPads[10][10];//maximum 10 pads width
  Int_t fNPads; //number of pads
  Int_t fMinpad; //??
  Int_t fMintime; //??
};

class AliHLTHistogram;
class TH2;
class TH3;
class AliHLTTrackArray;

class AliHLTHoughTest {

 public:
  
  AliHLTHoughTest(); 
  virtual ~AliHLTHoughTest();
  
  Bool_t GenerateTrackData(Double_t pt,Double_t psi,Double_t tgl,Int_t sign,Int_t patch,Int_t minhits);
  void FillImage(TH2 *hist,Int_t row=-1);
  void Transform2Circle(AliHLTHistogram *hist);
  void Transform2CircleC(AliHLTHistogram *hist);
  void Transform2CircleF(AliHLTHistogram *hist);
  void Transform2Line(AliHLTHistogram *hist,Int_t *rowrange);
  void Transform2LineC(AliHLTHistogram *hist,Int_t *rowrange);
  void Transform2Line3D(TH3 *hist,Int_t *rowrange,Float_t *phirange);
  void Transform2LineC3D(TH3 *hist,Int_t *rowrange);
  void TransformLines2Circle(TH3 *hist,AliHLTTrackArray *tracks);
  void Transform2Center(AliHLTHistogram *hist);
  
  void FindAbsMaxima(TH3 *hist,Int_t zsearch,Float_t &maxx,Float_t &maxy,Float_t &maxz,Int_t &maxvalue) const;
  
 private:
  AliHLTSimData *fData; //??
  Int_t fCurrentPatch; //index of the current patch
  
  ClassDef(AliHLTHoughTest,1) //Hough transform base class
};

typedef AliHLTHoughTest AliL3HoughTest; // for backward compatibility

#endif







