// @(#) $Id$

#ifndef ALIL3_HOUGHTEST
#define ALIL3_HOUGHTEST

#include "AliL3RootTypes.h"

struct SimData {
  Int_t pads[10][10];//maximum 10 pads width
  Int_t npads;
  Int_t minpad;
  Int_t mintime;
};

class AliL3Histogram;
class TH2;
class TH3;
class AliL3TrackArray;

class AliL3HoughTest {
  
 private:
  SimData *fData;
  Int_t fCurrentPatch;

 public:
  
  AliL3HoughTest(); 
  virtual ~AliL3HoughTest();
  
  Bool_t GenerateTrackData(Double_t pt,Double_t psi,Double_t tgl,Int_t sign,Int_t patch,Int_t minhits);
  void FillImage(TH2 *hist,Int_t row=-1);
  void Transform2Circle(AliL3Histogram *hist);
  void Transform2CircleC(AliL3Histogram *hist);
  void Transform2CircleF(AliL3Histogram *hist);
  void Transform2Line(AliL3Histogram *hist,Int_t *rowrange);
  void Transform2LineC(AliL3Histogram *hist,Int_t *rowrange);
  void Transform2Line3D(TH3 *hist,Int_t *rowrange,Float_t *phirange);
  void Transform2LineC3D(TH3 *hist,Int_t *rowrange);
  void TransformLines2Circle(TH3 *hist,AliL3TrackArray *tracks);
  void Transform2Center(AliL3Histogram *hist);
  
  void FindAbsMaxima(TH3 *hist,Int_t zsearch,Float_t &max_x,Float_t &max_y,Float_t &max_z,Int_t &maxvalue);
  
  ClassDef(AliL3HoughTest,1) //Hough transform base class
};

#endif







