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
  AliTRDrawData(const AliTRDrawData &r);
  virtual ~AliTRDrawData();
  AliTRDrawData &operator=(const AliTRDrawData &r);

  virtual void                 Copy(TObject &r);

  virtual Bool_t               Digits2Raw(TTree *digits);
  virtual AliTRDdigitsManager* Raw2Digits(AliRawReader* rawReader);
  virtual void                 SetDebug(Int_t v = 1) { fDebug = v; };

 protected:

  Int_t                fDebug;          //  Debug level

  ClassDef(AliTRDrawData,1)             //  TRD raw data class

};
#endif
