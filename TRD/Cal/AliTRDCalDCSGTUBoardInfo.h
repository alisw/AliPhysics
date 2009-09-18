#ifndef ALITRDCALDCSGTUBOARDINFO_H
#define ALITRDCALDCSGTUBOARDINFO_H
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

  TString GetId() const                               { return fId;                           }
  Int_t   GetType() const                             { return fType;                         }
  Int_t   GetPciGa() const                            { return fPciGa;                        }

  void    SetId(TString id)                           { fId = id;                             }
  void    SetType(Int_t ty)                           { fType = ty;                           }
  void    SetPciGa(Int_t ga)                          { fPciGa = ga;                          }

  TString GetHwDate() const                           { return fHwDate;                       }
  Int_t   GetHwRev() const                            { return fHwRev;                        }
  Int_t   GetHwClean() const                          { return fHwClean;                      }

  void    SetHwDate(TString hd)                       { fHwDate = hd;                         }
  void    SetHwRev(Int_t hr)                          { fHwRev = hr;                          }
  void    SetHwClean(Int_t hc)                        { fHwClean = hc;                        }

  TString GetSwDate() const                           { return fSwDate;                       }
  Int_t   GetSwRev() const                            { return fSwRev;                        }
  Int_t   GetSwClean() const                          { return fSwClean;                      }

  void    SetSwDate(TString sd)                       { fSwDate = sd;                         }
  void    SetSwRev(Int_t sr)                          { fSwRev = sr;                          }
  void    SetSwClean(Int_t sc)                        { fSwClean = sc;                        }

  protected:
  TString fId; // value from the board_id attribute of the board_info tag within a board_info tag
  Int_t   fType; // value from the design_type attribute of the board_info tag within a board_info tag
  Int_t   fPciGa; // value from the pci_ga attribute of the board_info tag within a board_info tag

  TString fHwDate; // value from the date attribute of the hardware tag within a board_info tag
  Int_t   fHwRev; // value from the rev attribute of the hardware tag within a board_info tag
  Int_t   fHwClean; // value from the clean attribute of the hardware tag within a board_info tag

  TString fSwDate; // value from the date attribute of the software tag within a board_info tag
  Int_t   fSwRev; // value from the rev attribute of the software tag within a board_info tag
  Int_t   fSwClean; // value from the clean attribute of the software tag within a board_info tag

  ClassDef(AliTRDCalDCSGTUBoardInfo,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
