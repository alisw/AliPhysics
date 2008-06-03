#ifndef AliTRDCALDCSGTU_H
#define AliTRDCALDCSGTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSGTU : public TNamed {

 public:

  AliTRDCalDCSGTU();
  AliTRDCalDCSGTU(const char *name, const char *title);
  virtual ~AliTRDCalDCSGTU() { };

  void    SetDCSid(Int_t dcsid)                    { fDCSID                     = dcsid; }  
  void    SetSMMaskBit(Int_t smid, Int_t mbit)     { fSMMask[smid]               = mbit; }
  void    SetStackMaskBit(Int_t smid, Int_t stid, Int_t mbit)
                                                   { fStackMask[smid][stid]      = mbit; }
  void    SetLinkMaskBit(Int_t smid, Int_t stid, Int_t lkid, Int_t mbit)
                                                   { fLinkMask[smid][stid][lkid] = mbit; }
  Int_t   SetSMMask(const char *smmask);
  Int_t   SetLinkMask(Int_t smid, Int_t stid, const char *lkmask);

  Int_t   GetDCSid() const                         { return fDCSID;                      }
  char    GetSMMaskBit(Int_t smid) const           { return fSMMask[smid];               }
  char    GetStackMaskBit(Int_t smid, Int_t stid) const
                                                   { return fStackMask[smid][stid];      }
  char    GetLinkMaskBit(Int_t smid, Int_t stid, Int_t lkid) const
                                                   { return fLinkMask[smid][stid][lkid]; }

 protected:
  
  Int_t   fDCSID;                  //  ID of the DCS-Board
  Int_t	  fSMMask[18];		   //  supermodule mask [SM-ID]
  Int_t	  fStackMask[18][5];	   //  stack mask [SM-ID][Stack-ID]
  Int_t   fLinkMask[18][5][12];	   //  link mask [SM-ID][Stack-ID][Link-ID]

  ClassDef(AliTRDCalDCSGTU,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
