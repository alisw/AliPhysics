
#include "AliL3Transform.h"
#include "AliL3HoughPixel.h"

ClassImp(AliL3HoughPixel)

AliL3HoughPixel::AliL3HoughPixel()
{

}

AliL3HoughPixel::AliL3HoughPixel(Int_t pad,Int_t time,Short_t signal)
  //				 Int_t slice,Int_t padrow)
{
  //Constructor
  fPad = pad;
  fTime = time;
  fSignal = signal;
  /*fSlice = slice;
  fPadrow = padrow;
  */
}

AliL3HoughPixel::~AliL3HoughPixel()
{

}
/*
void AliL3HoughPixel::SetCoordinates()
{

  AliL3Transform *transform = new AliL3Transform();
  
  Float_t xyz[3];
  Int_t sector,row;
  transform->Slice2Sector(fSlice,fPadrow,sector,row);
  transform->Raw2Global(xyz,sector,row,fPad,fTime);
  fX = xyz[0];
  fY = xyz[1];
  fZ = xyz[2];
  fPhi = transform->GetPhi(xyz);
  fEta = transform->GetEta(xyz);
  
}
*/
