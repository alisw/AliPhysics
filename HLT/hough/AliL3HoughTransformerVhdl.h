#ifndef ALIL3_HOUGHTRANSFORMERVHDL
#define ALIL3_HOUGHTRANSFORMERVDHL

#include "AliL3RootTypes.h"
#include "AliL3HoughBaseTransformer.h"

class AliL3Histogram;

class AliL3HoughTransformerVhdl : public AliL3HoughBaseTransformer {
  
 private:
  
  AliL3Histogram **fParamSpace; //!
  void DeleteHistograms();

 public:
  AliL3HoughTransformerVhdl(); 
  AliL3HoughTransformerVhdl(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformerVhdl();
  
  void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax);
  void CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
			Int_t nybin,Double_t ymin,Double_t ymax);
  void Reset();
  void TransformCircle();
  void TransformCircleC(Int_t row_range);
  void TransformLine();

  Int_t GetEtaIndex(Double_t eta);
  AliL3Histogram *GetHistogram(Int_t eta_index);
  

  ClassDef(AliL3HoughTransformerVhdl,1) //Normal Hough transformation class

};

inline AliL3Histogram *AliL3HoughTransformerVhdl::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

#endif




