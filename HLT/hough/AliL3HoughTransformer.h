#ifndef ALIL3_HOUGHTRANSFORMER
#define ALIL3_HOUGHTRANSFORMER

#include "AliL3RootTypes.h"

class AliL3Transform;
class AliL3Histogram;
class AliL3DigitRowData;

class AliL3HoughTransformer : public TObject {
  
 private:

  Int_t fSlice;
  Int_t fPatch;
  Int_t fNEtaSegments;
  Double_t fEtaMin;
  Double_t fEtaMax;
  Int_t fThreshold;
  AliL3Transform *fTransform; //!
  
  //Pointers to histograms
  AliL3Histogram **fParamSpace; //!

  //Data pointers
  UInt_t fNDigitRowData;
  AliL3DigitRowData *fDigitRowData; //!
  
  void DeleteHistograms();

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughTransformer();
  
  void SetInputData(UInt_t ndigits,AliL3DigitRowData *ptr);
  void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax);
  void CreateHistograms(Int_t nxbin=64,Double_t xmin=-0.006,Double_t xmax=0.006,
			Int_t nybin=64,Double_t ymin=-0.26,Double_t ymax=0.26);
  void Reset();
  void TransformCircle();
  void TransformCircleC();
  void TransformLine();

  //Getters
  Int_t GetSlice() {return fSlice;}
  Int_t GetPatch() {return fPatch;}
  Int_t GetNEtaSegments() {return fNEtaSegments;}
  Int_t GetThreshold() {return fThreshold;}
  Double_t GetEtaMin() {return fEtaMin;}
  Double_t GetEtaMax() {return fEtaMax;}
  Double_t GetEtaSlice() {return (fEtaMax - fEtaMin)/fNEtaSegments;}
  void *GetDataPointer() {return (void*)fDigitRowData;}
  AliL3Histogram *GetHistogram(Int_t eta_index);
  
  //setters
  void SetThreshold(Int_t i) {fThreshold = i;}

  ClassDef(AliL3HoughTransformer,1) //Hough transformation class

};

inline AliL3Histogram *AliL3HoughTransformer::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= fNEtaSegments || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

#endif
