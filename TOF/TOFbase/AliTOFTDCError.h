#ifndef ALITOFTDCERROR_H
#define ALITOFTDCERROR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  author: Roberto Preghenella (R+), preghenella@bo.infn.it
*/

///////////////////////////////////////////////////////////////
//                                                           //
//   This class provides a definition for TDC errors.        //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFTDCError : public TObject
{
 public:
  AliTOFTDCError(); //default constructor
  AliTOFTDCError(const AliTOFTDCError &source); //copy constructor
  AliTOFTDCError &operator = (const AliTOFTDCError &source); //operator =
  virtual ~AliTOFTDCError(); //destructor
  /* getters */
  UShort_t GetErrorFlags() {return fErrorFlags;}; //get error flags
  UShort_t GetTDCID() {return fTDCID;}; //get TDC ID
  /* setters */
  void SetErrorFlags(UShort_t ErrorFlags) {fErrorFlags = ErrorFlags;}; //set error flags
  void SetTDCID(UShort_t TDCID) {fTDCID = TDCID;};
 private:
  UShort_t fErrorFlags; //error flags
  UShort_t fTDCID; //TDC ID
  
  ClassDef(AliTOFTDCError, 1);
};

#endif /* ALITOFTDCERROR_H */
