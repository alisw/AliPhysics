#ifndef ALIL3_HOUGHTRANSFORMER
#define ALIL3_HOUGHTRANSFORMER

#include "AliL3RootTypes.h"

struct AliL3Digits
{
  UShort_t fCharge;
  UChar_t fPad;
  UShort_t fTime;
  Int_t fIndex;
  AliL3Digits *fNextVolumePixel;
  //AliL3Digits *nextPhiRowPixel;
  AliL3Digits *nextRowPixel;
};

struct AliL3HoughContainer
{
  void *first;
  void *last;
};



class AliL3Histogram;
class AliL3Transform;
class AliL3HoughEvaluate;
class AliL3TrackArray;
class AliL3DigitRowData;

class AliL3HoughTransformer : public TObject {
  
 private:

  friend class AliL3HoughEval;
  
  AliL3Transform *fTransform; //!

  Int_t fNPhiSegments; //Number of patches in phi.
  Float_t fEtaMin;
  Float_t fEtaMax;
  Int_t fNumEtaSegments;
  Int_t fNumOfPadRows;
  Int_t fNRowsInPatch;
  Int_t fBinTableBounds;

  UInt_t fNDigitRowData; //!
  AliL3DigitRowData *fDigitRowData; //!
  
  AliL3HoughContainer *fRowContainer; //!
  AliL3HoughContainer *fPhiRowContainer; //!
  AliL3HoughContainer *fVolume; //!
  
  Int_t fContainerBounds;
  Int_t fNDigits;
  Int_t **fBinTable; //!
  Int_t *fEtaIndex; //!
  UChar_t **fTrackTable; //!
  AliL3Histogram *fHistoPt;
  
  Int_t fSlice;
  Int_t fPatch;

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Float_t *etarange);
  AliL3HoughTransformer(Int_t slice,Int_t patch,Double_t *etarange=0,Int_t n_eta_segments=1);
  virtual ~AliL3HoughTransformer();

  void InitTables();
  void TransformTables(AliL3Histogram **histos,AliL3Histogram **images=0);
  void SetInputData(UInt_t ndigits,AliL3DigitRowData *ptr);
  void WriteTables();
  void SetHistogram(AliL3Histogram *hist) {fHistoPt = hist;}
  /*
    void Transform2Circle(TH2F *hist,Int_t eta_index);
    void Transform2Circle(TH2F **histos,Int_t n_eta_segments,UInt_t ndigits,AliL3DigitRowData *ptr);
    void Transform2Line(TH2F *hist,Int_t ref_row,Int_t *rowrange,Double_t *phirange,TH2F *raw=0);
    void TransformLines2Circle(TH2F *hist,AliL3TrackArray *tracks);
    void GetPixels(Char_t *rootfile,TH2F *hist=0);
    void InitTemplates(TH2F *hist);
  */
  ClassDef(AliL3HoughTransformer,1)

};

#endif
