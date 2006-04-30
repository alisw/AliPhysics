// @(#) $Id$

#ifndef ALIL3HOUGHTRANSFORMERVHDL_H
#define ALIL3HOUGHTRANSFORMERVHDL_H

#include "AliL3Histogram.h"
#include "AliL3HoughTransformerLUT.h"
class AliL3Histogram;

class AliL3HoughTransformerVhdl : public AliL3HoughTransformerLUT 
{

 public:

  AliL3HoughTransformerVhdl(); 
  AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t netasegments,Int_t nits=0);
  virtual ~AliL3HoughTransformerVhdl();

  void CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t pres,Int_t nybin,Float_t psi) {
    AliL3HoughTransformerLUT::CreateHistograms(ptmin,ptmax,pres,nybin,psi);
  }
  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);

  void TransformCircle();
  void TransformCircle(Int_t *row_range,Int_t every) {
    AliL3HoughTransformerLUT::TransformCircle(row_range,every);
  }
  void TransformCircleC(Int_t */*rowrange*/,Int_t /*every*/) {return;}
  
  void Init(Int_t slice=0,Int_t patch=0,Int_t netasegments=100,Int_t nits=-1);
  void Print();
  void PrintVhdl() const;

 protected:
  Float_t fEpsilon;//??
  Float_t fSinEpsilon;//??
  Float_t fCosEpsilon;//??
  Int_t fIts;//??

  Int_t fNxbin;//Number of bins in X
  Float_t fXmin;//Lower limit in X
  Float_t fXmax;//Upper limit in X
  Int_t fNybin;//Number of bins in Y
  Float_t fYmin;//Lower limit in Y
  Float_t fYmax;//Upper limit in Y

  ClassDef(AliL3HoughTransformerVhdl,1) //VHDL Hough transformation class

};

#endif






