#ifndef ALITRDCALDCSPTR_H
#define ALITRDCALDCSPTR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTR.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD PTR configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTR : public TNamed {

 public:

  AliTRDCalDCSPTR();
  AliTRDCalDCSPTR(const char *name, const char *title);
  virtual ~AliTRDCalDCSPTR() { };

  void    SetDCSid(Int_t dcsid)              { fDCSID               = dcsid; }  
  void    SetConfigID(TString id)            { fConfigID            = id;    }

  Int_t   GetDCSid() const                   { return fDCSID;                }
  TString GetConfigID() const                { return fConfigID;             }

 protected:
  
  Int_t   fDCSID;                     //  ID of the DCS-Board

  TString fConfigID;                  //  Configuration ID

  ClassDef(AliTRDCalDCSPTR,1)         //  TRD calibration class for TRD PTR parameters

};
#endif
