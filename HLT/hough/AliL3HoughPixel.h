#ifndef ALI_PIXEL
#define ALI_PIXEL

#include <TObject.h>

class AliL3HoughPixel : public TObject {

 private: 
  Short_t fSignal;
  Int_t fPad;
  Int_t fTime;
  /*Int_t fSlice;
  Int_t fPadrow;

  Float_t fX;
  Float_t fY;
  Float_t fZ;
  Double_t fPhi;
  Double_t fEta;
  */
 public:
  
  AliL3HoughPixel();
  //AliL3HoughPixel(Int_t pad,Int_t time,Short_t signal,Int_t slice,Int_t padrow);
  AliL3HoughPixel(Int_t pad,Int_t time,Short_t signal);
  virtual ~AliL3HoughPixel();

  /*Int_t GetPad() {return fPad;}
  Int_t GetTime() {return fTime;}
  Short_t GetSignal() {return fSignal;}
  Double_t GetPhi() {return fPhi;}
  Double_t GetEta() {return fEta;}
  Float_t GetX() {return fX;}
  Float_t GetY() {return fY;}
  Float_t GetZ() {return fZ;}

  void SetCoordinates();
  */
  AliL3HoughPixel *nextRowPixel;

  ClassDef(AliL3HoughPixel,1)
};


#endif
