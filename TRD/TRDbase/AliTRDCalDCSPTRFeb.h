#ifndef ALITRDCALDCSPTRFEB_H
#define ALITRDCALDCSPTRFEB_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSPTRFeb.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSPTRFeb : public TNamed {

 public:

  AliTRDCalDCSPTRFeb();
  AliTRDCalDCSPTRFeb(const char *name, const char *title);
  AliTRDCalDCSPTRFeb(const AliTRDCalDCSPTRFeb &);
  virtual ~AliTRDCalDCSPTRFeb() { };

  TString GetControlBoxSide() const                   { return fSide;                         }
  TString GetDetectorName() const                     { return fDetName;                      }
  Int_t   GetControlBoxPrimary() const                { return fPrimary;                      }

  void    SetControlBoxSide(TString bs)               { fSide = bs;                           }
  void    SetDetectorName(TString bs)                 { fDetName = bs;                        }
  void    SetControlBoxPrimary(Int_t bp)              { fPrimary = bp;                        }
  

 protected:
  TString fSide; // side of the control box, either A, B or C 
  TString fDetName; // Name of the detector eg  TO, V1 etc..
  Int_t   fPrimary; // 1 if its the primary control box, 2 for backup

  ClassDef(AliTRDCalDCSPTRFeb,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
