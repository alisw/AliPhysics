#ifndef ALITOFDCSMAPS_H
#define ALITOFDCSMAPS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id:  $ */

///////////////////////////////
//                           //
// AliTOFDCSmaps class       //
// container for HV&&LV maps //
// as found during a run     //
//                           //
///////////////////////////////

#include "TObject.h" 

class AliTOFDCSmaps : public TObject {

 public:
  AliTOFDCSmaps();
  AliTOFDCSmaps(Int_t time, Short_t array[]);
  AliTOFDCSmaps(const AliTOFDCSmaps & digit);
  AliTOFDCSmaps& operator=(const AliTOFDCSmaps & digit) ;
  virtual ~AliTOFDCSmaps();

  Int_t GetTime() const { return fTime; };
  void SetTime(Int_t time) { fTime=time; };
  Short_t * GetArray() { return fArray; };
  Short_t GetCellValue(Int_t index) const { return fArray[index]; };
  void SetCellValue(Int_t index, Short_t val) { fArray[index]=val; };
  void Update(AliTOFDCSmaps *object);
 private: 
  Int_t fTime; // time stamp
  Short_t fArray[91*96*18];  // array

  ClassDef(AliTOFDCSmaps, 1);

};

#endif
