////////////////////////////////////////////////
//  Digitization class for set: TOF           //
//  AliTOFRoc class                           //
//  Interface                                 //
//  Description                               // 
//*-- Authors: Pierella, Seganti, Vicinanza   //
//    (Bologna and Salerno University)        //
////////////////////////////////////////////////

#ifndef ALITOFROC_H
#define ALITOFROC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "TObject.h"
#include "TClonesArray.h"

//_______________________________________________________

class AliTOFRoc : public TObject {

 public:
  AliTOFRoc();
// copy ctor  (required also by RC10 Coding Convention)
  AliTOFRoc(const AliTOFRoc& tofroc);
// assignment operator (required also by RC10 Coding Convention)
  AliTOFRoc& operator = (const AliTOFRoc& tofroc);
// dtor
  virtual ~AliTOFRoc();
  Int_t   AddItem  (Int_t Fec, Int_t Tdc, Int_t Error, Float_t Charge, Float_t Time);
//  Int_t   AddItem  (Int_t, UInt_t);

// setters for AliTOFRoc object
  void    SetHeadVar(Int_t items, Int_t size, Int_t nroc, UInt_t header);
  void    SetHeader();
  void    SetTime  (UInt_t Item, UInt_t Error, Float_t RealTime);
  void    SetTime  (UInt_t Item, UInt_t tir);
  void    SetCharge(UInt_t Item, UInt_t Fec,UInt_t Tdc,Float_t RealCharge);  
  void    SetCharge(UInt_t Item, UInt_t chr);  

// getters for AliTOFRoc object
  Float_t GetTime  (Int_t Item,UInt_t& Error);
  Float_t GetCharge(Int_t Item) const;
  Int_t   GetTotPad(Int_t Item) const;
  UInt_t  GetCheckSum();
  UInt_t  BitCount (UInt_t x) const;
  UInt_t  SetSize  ();
  
  Int_t  GetSize()           const {return fItems*8+4;}
  Int_t  GetItems()          const {return fItems;}
  UInt_t GetChrgRow(Int_t i) const {return fChrgRow[i];}
  UInt_t GetTimeRow(Int_t i) const {return fTimeRow[i];} 
  void   SetHeader(UInt_t head){fHeader=head;}

 protected:
  Int_t   fItems;  // number of items
  Int_t   fSize;   // size
  Int_t   fNRoc;   // Roc number
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

  ClassDef(AliTOFRoc,2) // TOF Read Out Controller class 
};

#endif /* ALITOFROC_H */
