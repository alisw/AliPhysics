#ifndef ALIL3_HOUGHBASETRANSFORMER
#define ALIL3_HOUGHBASETRANSFORMER

#include "AliL3RootTypes.h"

#ifdef do_mc
const UInt_t MaxTrack=15;
struct TrackIndex {
  Int_t fLabel[MaxTrack];
  Int_t fNHits[MaxTrack];
};
typedef struct TrackIndex TrackIndex;
#endif

class AliL3DigitRowData;
class AliL3Histogram;

class AliL3HoughBaseTransformer {
  
 private:

  Int_t fSlice;
  Int_t fPatch;
  Int_t fNEtaSegments;
  Double_t fEtaMin;
  Double_t fEtaMax;
  Int_t fLowerThreshold;
  Int_t fUpperThreshold;
  
  AliL3DigitRowData *fDigitRowData; //!
  
 public:
  AliL3HoughBaseTransformer(); 
  AliL3HoughBaseTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments);
  virtual ~AliL3HoughBaseTransformer();
  
  void SetInputData(UInt_t ndigits,AliL3DigitRowData *ptr) {fDigitRowData = ptr;}
  
  virtual void CreateHistograms(Int_t nxbin,Double_t ptmin,Int_t nybin,Double_t phimin,Double_t phimax) = 0;
  virtual void CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,Int_t nybin,Double_t ymin,Double_t ymax) = 0;
  virtual void Reset() = 0;
  virtual void TransformCircle() = 0;
  virtual void TransformCircleC(Int_t row_range) = 0;
  virtual void TransformLine() = 0;
  
  //Getters
  Int_t GetSlice() {return fSlice;}
  Int_t GetPatch() {return fPatch;}
  Int_t GetNEtaSegments() {return fNEtaSegments;}
  Int_t GetLowerThreshold() {return fLowerThreshold;}
  Int_t GetUpperThreshold() {return fUpperThreshold;}
  Double_t GetEtaMin() {return fEtaMin;}
  Double_t GetEtaMax() {return fEtaMax;}
  
  AliL3DigitRowData *GetDataPointer() {return fDigitRowData;}
 
  virtual Int_t GetEtaIndex(Double_t eta) = 0;
  virtual AliL3Histogram *GetHistogram(Int_t eta_index) = 0;
  virtual Double_t GetEta(Int_t eta_index,Int_t slice) = 0;
  virtual Int_t GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi) {return -1;}
  
  //setters
  virtual void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100);
  void SetLowerThreshold(Int_t i) {fLowerThreshold = i;}
  void SetUpperThreshold(Int_t i) {fUpperThreshold = i;}

  virtual void Print(){};

  ClassDef(AliL3HoughBaseTransformer,1) //Hough transformation base class

};


#endif

