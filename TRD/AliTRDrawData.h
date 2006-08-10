#ifndef ALITRDRAWDATA_H
#define ALITRDRAWDATA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Converts TRD digits into a raw data stream                               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TTree;
class AliTRDdigitsManager;
class AliRawReader;

class AliTRDrawData : public TObject {

 public:

  AliTRDrawData();
  virtual ~AliTRDrawData();

  virtual Bool_t               Digits2Raw(TTree *digits);
  virtual AliTRDdigitsManager* Raw2Digits(AliRawReader* rawReader);

 protected:

  ClassDef(AliTRDrawData,2)             //  TRD raw data class

};
#endif
