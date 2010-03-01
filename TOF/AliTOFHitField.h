#ifndef ALITOFHITFIELD_H
#define ALITOFHITFIELD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides the minimum-size TOF hit info       //
//                                                           //
//   author: Roberto Preghenella (R+)                        //
//           preghenella@bo.infn.it                          //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TROOT.h"

class AliTOFHitField
{
  
 public:
  
  AliTOFHitField(); // default constructor
  virtual ~AliTOFHitField(); // default destructor
  AliTOFHitField(const AliTOFHitField &source); // copy constructor
  AliTOFHitField &operator=(const AliTOFHitField &source); // operator=

  UInt_t GetIndex() const {return fIndex;}; // get index
  UShort_t GetTimeBin() const {return fTimeBin;}; // get time bin
  UShort_t GetTOTBin() const {return fTOTBin;}; // get TOT bin
  UChar_t GetDeltaBC() const {return fDeltaBC;}; // get delta BC
  UShort_t GetL0L1Latency() const {return fL0L1Latency;}; // get L0-L1 latency

  void SetIndex(UInt_t value) {fIndex = value;}; // set index
  void SetTimeBin(UShort_t value) {fTimeBin = value;}; // set time bin
  void SetTOTBin(UShort_t value) {fTOTBin = value;}; // set TOT bin
  void SetDeltaBC(UChar_t value) {fDeltaBC = value;}; // set delta BC
  void SetL0L1Latency(UShort_t value) {fL0L1Latency = value;}; // set L0-L1 latency
  
 private:
  
  UInt_t fIndex; // channel index
  UShort_t fTimeBin; // time bin [24.4 ps]
  UShort_t fTOTBin; // TOT bin [48.8 ps]
  UChar_t fDeltaBC; // delta BC [BC bins]
  UShort_t fL0L1Latency; // L0-L1 latency [BC bins]
  
  ClassDef(AliTOFHitField, 1);
};

#endif /* ALITOFHITFIELD_H */
