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

class TH2F;
class AliL3Transform;
class AliL3HoughEvaluate;
class AliL3TrackArray;

class AliL3HoughTransformer : public TObject {
  
 private:

  friend class AliL3HoughEval;
  
  AliL3Transform *fTransform; //!

  Int_t fNPhiSegments; //Number of patches in phi.
  Float_t fEtaMin;
  Float_t fEtaMax;
  Int_t fNumEtaSegments;
  Int_t fNumOfPadRows;
  
  AliL3HoughContainer *fRowContainer; //!
  AliL3HoughContainer *fPhiRowContainer; //!
  AliL3HoughContainer *fVolume; //!
  Int_t fContainerBounds;
  Int_t fNDigits;
  Int_t **fIndex; //!
  
  
  Int_t fSlice;
  Int_t fPatch;

 public:
  AliL3HoughTransformer(); 
  AliL3HoughTransformer(Int_t slice,Int_t patch,Float_t *etarange);
  AliL3HoughTransformer(Int_t slice,Int_t patch,Double_t *etarange=0,Int_t n_eta_segments=1);
  virtual ~AliL3HoughTransformer();
  
  void Transform2Circle(TH2F *hist,Int_t eta_index);
  void Transform2Line(TH2F *hist,Int_t ref_row,Int_t *rowrange,Double_t *phirange,TH2F *raw=0);
  void TransformLines2Circle(TH2F *hist,AliL3TrackArray *tracks);
  void GetPixels(Char_t *rootfile,TH2F *hist=0);
  void InitTemplates(TH2F *hist);
  void CountBins();

  ClassDef(AliL3HoughTransformer,1)

};

#endif
