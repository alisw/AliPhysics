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

class AliTRDdigitsManager;

class AliTRDrawData : public TObject {

 public:

  AliTRDrawData();
  AliTRDrawData(const AliTRDrawData &r);
  virtual ~AliTRDrawData();
  AliTRDrawData &operator=(const AliTRDrawData &r);

  virtual void                 Copy(TObject &r);

  virtual Bool_t               OpenInput(const Char_t *name);
  virtual Bool_t               Digit2Raw(const Char_t *name1 = "trd_ldc0.d", 
                                         const Char_t *name2 = "trd_ldc1.d");
  virtual Bool_t               Raw2Digit(const Char_t *name1 = "trd_ldc0.d", 
                                         const Char_t *name2 = "trd_ldc1.d");
  virtual void                 SetDebug(Int_t v = 1) { fDebug = v; };
  virtual AliTRDdigitsManager *GetDigitsManager()    { return fDigitsManager; };

 protected:

  Int_t                fDebug;          //  Debug level
  AliTRDdigitsManager *fDigitsManager;  //! The TRD digits manager

  ClassDef(AliTRDrawData,1)             //  TRD raw data class

};
#endif
