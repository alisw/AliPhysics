// @(#) $Id$

#ifndef ALIL3_HOUGHTRANSFORMERVHDL
#define ALIL3_HOUGHTRANSFORMERVDHL

#include "AliL3Histogram.h"
#include "AliL3HoughTransformerLUT.h"
class AliL3Histogram;

class AliL3HoughTransformerVhdl : public AliL3HoughTransformerLUT 
{
 protected:
  Float_t fEpsilon;
  Float_t fSinEpsilon;
  Float_t fCosEpsilon;
  Int_t fIts;

  Int_t fNxbin;
  Float_t fXmin;
  Float_t fXmax;
  Int_t fNybin;
  Float_t fYmin;
  Float_t fYmax;

 public:

  AliL3HoughTransformerVhdl(); 
  AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments,Int_t n_its=0);
  virtual ~AliL3HoughTransformerVhdl();

  void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax);
  void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
			Int_t nybin,Float_t ymin,Float_t ymax);

  void TransformCircle();
  void TransformCircleC(Int_t */*rowrange*/,Int_t /*every*/) {return;}
  
  void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100,Int_t n_its=-1);
  void Print();
  void PrintVhdl();

  ClassDef(AliL3HoughTransformerVhdl,1) //VHDL Hough transformation class

};

#endif






