// @(#) $Id$

#ifndef ALIL3_HOUGHBASETRANSFORMER
#define ALIL3_HOUGHBASETRANSFORMER

#include "AliL3RootTypes.h"

#ifdef do_mc
const UInt_t MaxTrack=120;
struct TrackIndex {
  Int_t fLabel[MaxTrack];
  UChar_t fNHits[MaxTrack];
  UChar_t fCurrentRow[MaxTrack];
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

  Float_t fZVertex;

 public:

  AliL3HoughBaseTransformer(); 
  AliL3HoughBaseTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments,Float_t zvertex=0.0);
  virtual ~AliL3HoughBaseTransformer();
  
  void SetInputData(UInt_t /*ndigits*/,AliL3DigitRowData *ptr) {fDigitRowData = ptr;}
  
  //this is for adaptave histograms
  virtual void CreateHistograms(Float_t /*ptmin*/,Float_t /*ptmax*/,Float_t /*pres*/,Int_t /*nybin*/,Float_t /*psi*/)
    {STDCERR<<"Adaptive histograms are not supported  for this Transformer!"<<STDENDL;}

  virtual void CreateHistograms(Int_t nxbin,Float_t ptmin,Int_t nybin,Float_t phimin,Float_t phimax) = 0;
  virtual void CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,Int_t nybin,Float_t ymin,Float_t ymax) = 0;

  virtual void Reset() = 0;
  virtual void TransformCircle()
    {STDCERR<<"TransformCircle is not defined for this transformer!"<<STDENDL;}
  virtual void TransformCircle(Int_t */*row_range*/,Int_t /*every*/)
    {STDCERR<<"TransformCircle is not defined for this transformer!"<<STDENDL;}
  virtual void TransformCircleC(Int_t */*row_range*/,Int_t /*every*/)
    {STDCERR<<"TransformCircleC is not defined for this transformer!"<<STDENDL;}
  virtual void TransformLine(Int_t */*rowrange*/=0,Float_t */*phirange*/=0)
    {STDCERR<<"TransformLine is not defined for this Transformer!"<<STDENDL;}
  virtual void TransformLineC(Int_t */*rowrange*/,Float_t */*phirange*/)
      {STDCERR<<"TransformLineC is not defined for this Transformer!"<<STDENDL;}

  //Getters
  Int_t GetSlice() const {return fSlice;}
  Int_t GetPatch() const {return fPatch;}
  Int_t GetNEtaSegments() const {return fNEtaSegments;}
  Int_t GetLowerThreshold() const {return fLowerThreshold;}
  Int_t GetUpperThreshold() const {return fUpperThreshold;}
  Double_t GetEtaMin() const {return fEtaMin;}
  Double_t GetEtaMax() const {return fEtaMax;}
  Float_t GetZVertex() const {return fZVertex;}

  AliL3DigitRowData *GetDataPointer() {return fDigitRowData;}
 
  virtual Int_t GetEtaIndex(Double_t eta) = 0;
  virtual void GetEtaIndexes(Double_t /*eta*/,Int_t */*indexes*/)
    {STDCERR<<"GetEtaIndexes not implemented for this Transformer class"<<STDENDL;}
  virtual AliL3Histogram *GetHistogram(Int_t eta_index) = 0;
  virtual Double_t GetEta(Int_t eta_index,Int_t slice) const = 0;

  virtual Int_t GetTrackID(Int_t /*eta_index*/,Double_t /*kappa*/,Double_t /*psi*/){
    STDCERR<<"GetTrackID not implemented for this Transformer class"<<STDENDL; 
    return -1;
  }
  
  //setters
  virtual void Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100,Int_t n_segs=-1);
  void SetLowerThreshold(Int_t i) {fLowerThreshold = i;}
  void SetUpperThreshold(Int_t i) {fUpperThreshold = i;}

  virtual void Print(){};

  ClassDef(AliL3HoughBaseTransformer,1) //Hough transformation base class

};


#endif

