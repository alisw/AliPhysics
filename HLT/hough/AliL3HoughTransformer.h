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
  
  //Pointer to histograms
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
  AliL3DigitRowData *UpdateDataPointer(AliL3DigitRowData *tempPt);
  void CreateHistograms(Int_t nxbin=70,Double_t xmin=-0.006,Double_t xmax=0.006,
			Int_t nybin=70,Double_t ymin=-0.26,Double_t ymax=0.26);
  void Transform();
  
  Int_t GetSlice() {return fSlice;}
  Int_t GetPatch() {return fPatch;}
  Int_t GetNEtaSegments() {return fNEtaSegments;}
  Double_t GetEtaMin() {return fEtaMin;}
  Double_t GetEtaMax() {return fEtaMax;}
  void *GetDataPointer() {return (void*)fDigitRowData;}
  AliL3Histogram *GetHistogram(Int_t eta_index);

  ClassDef(AliL3HoughTransformer,1)

};

inline AliL3Histogram *AliL3HoughTransformer::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= fNEtaSegments)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

#endif
