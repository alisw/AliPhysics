#ifndef ALITOFTDCHIT_H
#define ALITOFTDCHIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a definition for TDC hits.          //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"
#include "AliTOFRawDataFormat.h"

class AliTOFTDCHit : public TObject
{
 public:
  AliTOFTDCHit(); //default constructor
  AliTOFTDCHit(const AliTOFTDCHit &source); //copy contructor
  AliTOFTDCHit &operator = (const AliTOFTDCHit &source); //operator =
  //  AliTOFTDCHit &operator - (const AliTOFTDCHit &source); //operator -
  AliTOFTDCHit &operator -= (const AliTOFTDCHit &source); //operator -=
  AliTOFTDCHit &operator << (const AliTOFTDCHit &source); //operator <<
  virtual ~AliTOFTDCHit(); //destructor
  /* getters */
  UInt_t   GetHitTime() const {return fHitTime;}; //get hit time
  UShort_t GetTOTWidth() const {return fTOTWidth;}; //get TOT width
  UShort_t GetChan() const {return fChan;}; //get channel
  UShort_t GetTDCID() const {return fTDCID;}; //get TDC ID
  UShort_t GetEBit() const {return fEBit;}; //get E bit
  UShort_t GetPSBits() const {return fPSBits;}; //get PS bits
  /* setters */
  void SetHitTime(UInt_t HitTime) {fHitTime = HitTime;}; //set hit time
  void SetTOTWidth(UShort_t TOTWidth) {fTOTWidth = TOTWidth;}; //set TOT width
  void SetChan(UShort_t Chan) {fChan = Chan;}; //set channel
  void SetTDCID(UShort_t TDCID) {fTDCID = TDCID;}; //set TDC ID
  void SetEBit(UShort_t EBit) {fEBit = EBit;};
  void SetPSBits(UShort_t PSBits) {fPSBits = PSBits;}; //set PS bits
 private:
  UInt_t   fHitTime; //hit time [24.4 ps]
  UShort_t fTOTWidth; //TOT width [48.8 ps]
  UShort_t fChan; //channel
  UShort_t fTDCID; //TDC ID
  UShort_t fEBit; //E bit
  UShort_t fPSBits; //PS bits

  ClassDef(AliTOFTDCHit, 1);
};

#endif /* ALITOFTDCHIT_H */
