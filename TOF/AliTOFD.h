#ifndef TOFD_H
#define TOFD_H

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
  Int_t   Items;
  Int_t   Size;
  Int_t   NRoc;
  UInt_t  Header;
  
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
  UInt_t ChrgRow[1024];
  UInt_t TimeRow[1024];

 public:
  AliTOFRoc();
  virtual ~AliTOFRoc();
  Int_t   AddItem  (Int_t, Int_t, Int_t, Float_t, Float_t);
//  Int_t   AddItem  (Int_t, UInt_t);
  void    SetHeader();
  void    SetTime  (UInt_t, UInt_t, Float_t);
  void    SetTime  (UInt_t, UInt_t);
  void    SetCharge(UInt_t, UInt_t,UInt_t,Float_t);  
  void    SetCharge(UInt_t, UInt_t);  
  Float_t GetTime  (Int_t,UInt_t&);
  Float_t GetCharge(Int_t);
  Int_t   GetTotPad(Int_t);
  UInt_t  GetCheckSum();
  UInt_t  BitCount (UInt_t);
  UInt_t  SetSize  ();
  
  inline Int_t  GetSize()          {return Items*8+4;}
  inline Int_t  GetItems()         {return Items;}
  inline UInt_t GetChrgRow(Int_t i){return ChrgRow[i];}
  inline UInt_t GetTimeRow(Int_t i){return TimeRow[i];} 
  inline void   SetHeader(UInt_t head){Header=head;}
  ClassDef(AliTOFRoc,2)
};

//_______________________________________________________
class AliTOFRawDigit : public TObject{

public:
  Int_t fTreeD;
  Int_t fRawDigits;
  
public:
  AliTOFRawDigit();
  virtual ~AliTOFRawDigit(){};
  
  ClassDef(AliTOFRawDigit,2)
};


//_______________________________________________________
class AliTOFRawSector : public TObject{

 public:
  TClonesArray  *fRocData;
  UInt_t        Header;
  UInt_t        GlobalCheckSum;

 public:
  AliTOFRawSector();
  virtual ~AliTOFRawSector();
  void   WriteSector();
  void   ReadSector();
  
  inline TClonesArray  *GetRocData() {return fRocData;}
  inline void SetGlobalCS(UInt_t gcs){GlobalCheckSum=gcs;}
  inline void SetHeader  (UInt_t hdr){Header = hdr;}
  
  ClassDef(AliTOFRawSector,2)
};

#endif
