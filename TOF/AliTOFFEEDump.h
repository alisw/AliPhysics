#ifndef ALITOFFEEDUMP_H
#define ALITOFFEEDUMP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTOFDecoder.h,v 1.2 2007/05/08 11:55:24 arcelli Exp $ */

///////////////////////////////////////////////////////////////
//                                                           //
//   This classes provide the object to store the full dump  //
//   of TOF FEE configuration database.                      //
//                                                           //
///////////////////////////////////////////////////////////////

#include "TObject.h"

class AliTOFFEEDump :
public TObject
{

 public:

  AliTOFFEEDump(); /* default constructor */
  virtual ~AliTOFFEEDump(); /* default destructor */
  Bool_t operator!=(const AliTOFFEEDump &source); /* operator!= */
  Bool_t operator==(const AliTOFFEEDump &source) {return !(*this != source);}; /* operator== */

  UInt_t GetSize() const {return fSize;}; /* get size */
  UChar_t *GetData() {return fData;}; /* get data */
  Bool_t ReadFromFile(const Char_t *filename); /* read from file */
  void DumpData(); /* dump data */

 private:

  AliTOFFEEDump(const AliTOFFEEDump &source); /* copy constructor */
  AliTOFFEEDump &operator=(const AliTOFFEEDump &source); /* operator= */

  UInt_t fSize; // size
  UChar_t *fData; //[fSize] data

  ClassDef(AliTOFFEEDump, 1);

};

#endif /* ALITOFFEEDUMP_H */
