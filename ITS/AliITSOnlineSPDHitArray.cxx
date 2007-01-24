/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// One object for each half stave and step in a scan. It keeps //
// the nr of hits in each pixel.                               //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDHitArray.h"

ClassImp(AliITSOnlineSPDHitArray)

AliITSOnlineSPDHitArray::AliITSOnlineSPDHitArray() {
  // constructor
  for (Int_t i=0; i<81920; i++) {
    fHits[i]=0;
  }
}

AliITSOnlineSPDHitArray* AliITSOnlineSPDHitArray::CloneThis() const {
  // makes a copy of this object and returns it
  AliITSOnlineSPDHitArray* returnpointer = new AliITSOnlineSPDHitArray();
  for (Int_t chip=0; chip<10; chip++) {
    for (Int_t col=0; col<32; col++) {
      for (Int_t row=0; row<256; row++) {
	returnpointer->SetHits(chip,col,row,fHits[GetKey(chip,col,row)]);
      }
    }
  }
  return returnpointer;
}

void   AliITSOnlineSPDHitArray::IncrementHits(UInt_t chip, UInt_t col, UInt_t row) {
  fHits[GetKey(chip,col,row)] ++;
}
void   AliITSOnlineSPDHitArray::SetHits(UInt_t chip, UInt_t col, UInt_t row, UInt_t hits) {
  fHits[GetKey(chip,col,row)] = hits;
}
UInt_t AliITSOnlineSPDHitArray::GetHits(UInt_t chip, UInt_t col, UInt_t row) const {
  return fHits[GetKey(chip,col,row)];
}
UInt_t AliITSOnlineSPDHitArray::GetKey(UInt_t chip, UInt_t col, UInt_t row) const {
  return chip*256*32 + col*256 + row;
}
