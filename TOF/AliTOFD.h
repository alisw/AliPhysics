#ifndef ALITOFD_H
#define ALITOFD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Digitization classes for set: TOF     //
////////////////////////////////////////////////
 

#include "TObject.h"
#include "TClonesArray.h"
#include "AliTOF.h"
//_______________________________________________________

class AliTOFRoc : public TObject {

 public:
  AliTOFRoc();
  virtual ~AliTOFRoc();
  Int_t   AddItem  (Int_t Fec, Int_t Tdc, Int_t Error, Float_t Charge, Float_t Time);
//  Int_t   AddItem  (Int_t, UInt_t);
// setters for AliTOFRoc object
  void    SetHeader();
  void    SetTime  (UInt_t Item, UInt_t Error, Float_t RealTime);
  void    SetTime  (UInt_t Item, UInt_t tir);
  void    SetCharge(UInt_t Item, UInt_t Fec,UInt_t Tdc,Float_t RealCharge);  
  void    SetCharge(UInt_t Item, UInt_t chr);  
// getters for AliTOFRoc object
  Float_t GetTime  (Int_t Item,UInt_t& Error);
  Float_t GetCharge(Int_t Item);
  Int_t   GetTotPad(Int_t Item);
  UInt_t  GetCheckSum();
  UInt_t  BitCount (UInt_t x);
  UInt_t  SetSize  ();
  
  Int_t  GetSize()           const {return fItems*8+4;}
  Int_t  GetItems()          const {return fItems;}
  UInt_t GetChrgRow(Int_t i) const {return fChrgRow[i];}
  UInt_t GetTimeRow(Int_t i) const {return fTimeRow[i];} 
  void   SetHeader(UInt_t head){fHeader=head;}

 protected:
  Int_t   fItems; // number of items
  Int_t   fSize;  // size
  Int_t   fNRoc;  // Roc number
  UInt_t  fHeader; // Roc header number

/*  class ChargeRow
  {
  public:
    UInt_t RocID:4;
    UInt_t FecID:6;
    UInt_t TdcID:6;
    Int_t  ChADC:16;
  }Charge[1024];

  class TimeRow
  {
  public:
    UInt_t Error:12;
    Int_t  TDC  :24;
  }Time[1024];
*/
  UInt_t fChrgRow[1024]; // adc values
  UInt_t fTimeRow[1024]; // tdc values



  ClassDef(AliTOFRoc,2)
};

//_______________________________________________________
class AliTOFRawDigit : public TObject{
  
public:
  AliTOFRawDigit();
  virtual ~AliTOFRawDigit(){};

protected:
  Int_t fTreeD;     // class under construction
  Int_t fRawDigits; // class under construction

  
  ClassDef(AliTOFRawDigit,2)
};


//_______________________________________________________
class AliTOFRawSector : public TObject{

 public:
  AliTOFRawSector();
  virtual ~AliTOFRawSector();
  void   WriteSector();
  void   ReadSector();
  
  TClonesArray* GetRocData() const {return fRocData;}
  void SetGlobalCS(UInt_t gcs){fGlobalCheckSum=gcs;}
  void SetHeader  (UInt_t hdr){fHeader = hdr;}

 protected:
  TClonesArray* fRocData; // pointer to the TClonesArray of Roc Data
  UInt_t        fHeader;    // RawSector header number
  UInt_t        fGlobalCheckSum; // check flag

  
  ClassDef(AliTOFRawSector,2)
};

#endif /* ALITOFD_H */
