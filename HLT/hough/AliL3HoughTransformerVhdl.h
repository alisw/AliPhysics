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

 public:

  AliL3HoughTransformerVhdl(); 
  AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments,Int_t n_its=0);
  virtual ~AliL3HoughTransformerVhdl();

  void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax);
  void CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
			Int_t nybin,Double_t ymin,Double_t ymax);

  void TransformCircle();
#if 0


  void Reset();

  void TransformCircleC(Int_t row_range) {STDCERR<<"TransformCircleC is not defined!"<<STDENDL;}
  void TransformLine() {STDCERR<<"TransformLine is not defined!"<<STDENDL;}

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  Double_t GetEta(Int_t eta_index,Int_t slice);


#endif
  void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100,Int_t n_its=-1);
  void Print();

  ClassDef(AliL3HoughTransformerVhdl,1) //Normal Hough transformation class

};

#endif






