#ifndef AliTRDCALDCSGTUBoardInfo_H
#define AliTRDCALDCSGTUBoardInfo_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUBoardInfo.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;

class AliTRDCalDCSGTUBoardInfo : public TNamed {

 public:

  AliTRDCalDCSGTUBoardInfo();
  AliTRDCalDCSGTUBoardInfo(const char *name, const char *title);
  virtual ~AliTRDCalDCSGTUBoardInfo() { };

  TString GetId()                                     { return fId;                           }
  Int_t   GetType()                                   { return fType;                         }
  Int_t   GetPciGa()                                  { return fPciGa;                        }

  void    SetId(TString id)                           { fId = id;                             }
  void    SetType(Int_t ty)                           { fType = ty;                           }
  void    SetPciGa(Int_t ga)                          { fPciGa = ga;                          }

  TString GetHwDate()                                 { return fHwDate;                       }
  Int_t   GetHwRev()                                  { return fHwRev;                        }
  Int_t   GetHwClean()                                { return fHwClean;                      }

  void    SetHwDate(TString hd)                       { fHwDate = hd;                         }
  void    SetHwRev(Int_t hr)                          { fHwRev = hr;                          }
  void    SetHwClean(Int_t hc)                        { fHwClean = hc;                        }

  TString GetSwDate()                                 { return fSwDate;                       }
  Int_t   GetSwRev()                                  { return fSwRev;                        }
  Int_t   GetSwClean()                                { return fSwClean;                      }

  void    SetSwDate(TString sd)                       { fSwDate = sd;                         }
  void    SetSwRev(Int_t sr)                          { fSwRev = sr;                          }
  void    SetSwClean(Int_t sc)                        { fSwClean = sc;                        }

  protected:
  TString fId;
  Int_t   fType;
  Int_t   fPciGa;

  TString fHwDate;
  Int_t   fHwRev;
  Int_t   fHwClean;

  TString fSwDate;
  Int_t   fSwRev;
  Int_t   fSwClean;

  ClassDef(AliTRDCalDCSGTUBoardInfo,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
