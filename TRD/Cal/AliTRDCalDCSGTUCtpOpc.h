#ifndef AliTRDCALDCSGTUCtpOpc_H
#define AliTRDCALDCSGTUCtpOpc_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUCtpOpc.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSGTUCtpOpc : public TNamed {

 public:

  AliTRDCalDCSGTUCtpOpc();
  AliTRDCalDCSGTUCtpOpc(const char *name, const char *title);
  virtual ~AliTRDCalDCSGTUCtpOpc() { };

  Int_t   GetId()                                     { return fId;                           }
  Int_t   GetOpcode()                                 { return fOpcode;                       }
  Int_t   GetDirection()                              { return fDirection;                    }
  Int_t   GetInverted()                               { return fInverted;                     }
  Int_t   GetDelay()                                  { return fDelay;                        }
  Int_t   GetConnected()                              { return fConnected;                    }

  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetOpcode(Int_t op)                         { fOpcode = op;                         }
  void    SetDirection(Int_t di)                      { fDirection = di;                      }
  void    SetInverted(Int_t in)                       { fInverted = in;                       }
  void    SetDelay(Int_t de)                          { fDelay = de;                          }
  void    SetConnected(Int_t co)                      { fConnected = co;                      }

  protected:
  Int_t   fId;
  Int_t   fOpcode;
  Int_t   fDirection;
  Int_t   fInverted;
  Int_t   fDelay;
  Int_t   fConnected;

  ClassDef(AliTRDCalDCSGTUCtpOpc,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
