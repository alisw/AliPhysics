#ifndef ALITRDCALDCSGTUCTPOPC_H
#define ALITRDCALDCSGTUCTPOPC_H
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

  Int_t   GetId() const                               { return fId;                           }
  Int_t   GetOpcode() const                           { return fOpcode;                       }
  Int_t   GetDirection() const                        { return fDirection;                    }
  Int_t   GetInverted() const                         { return fInverted;                     }
  Int_t   GetDelay() const                            { return fDelay;                        }
  Int_t   GetConnected() const                        { return fConnected;                    }

  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetOpcode(Int_t op)                         { fOpcode = op;                         }
  void    SetDirection(Int_t di)                      { fDirection = di;                      }
  void    SetInverted(Int_t in)                       { fInverted = in;                       }
  void    SetDelay(Int_t de)                          { fDelay = de;                          }
  void    SetConnected(Int_t co)                      { fConnected = co;                      }

  protected:
  Int_t   fId; // value of the attribute named id within the otp_opc tag
  Int_t   fOpcode; // value of the attribute named opcode within the otp_opc tag
  Int_t   fDirection; // value of the attribute named direction within the otp_opc tag
  Int_t   fInverted; // value of the attribute named inverted within the otp_opc tag
  Int_t   fDelay; // value of the attribute named delay within the otp_opc tag
  Int_t   fConnected; // value of the attribute named connected within the otp_opc tag

  ClassDef(AliTRDCalDCSGTUCtpOpc,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
