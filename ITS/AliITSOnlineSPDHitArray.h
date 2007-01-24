#ifndef ALI_ITS_ONLINESPDHITARRAY_H
#define ALI_ITS_ONLINESPDHITARRAY_H

/////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                      //
// This class is used as a container online.                   //
// One object for each half stave and step in a scan. It keeps //
// the nr of hits in each pixel.                               //
// This class should only be used through the interface of the //
// AliITSOnlineSPDscan class.                                  //
/////////////////////////////////////////////////////////////////

#include <TObject.h>

class AliITSOnlineSPDHitArray : public TObject {
 public:
  AliITSOnlineSPDHitArray();
  virtual ~AliITSOnlineSPDHitArray(){}
  void   IncrementHits(UInt_t chip, UInt_t col, UInt_t row);
  void   SetHits(UInt_t chip, UInt_t col, UInt_t row, UInt_t hits);
  UInt_t GetHits(UInt_t chip, UInt_t col, UInt_t row) const;
  AliITSOnlineSPDHitArray* CloneThis() const;

 private:
  UInt_t fHits[81920]; // nr of hits for each pixel of this half stave
  UInt_t GetKey(UInt_t chip, UInt_t col, UInt_t row) const;
  ClassDef(AliITSOnlineSPDHitArray,1)
    };

    
#endif
