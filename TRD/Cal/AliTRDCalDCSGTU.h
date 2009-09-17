#ifndef ALITRDCALDCSGTU_H
#define ALITRDCALDCSGTU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTU.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"
#include "AliTRDCalDCSGTUTgu.h"

class TString;

class AliTRDCalDCSGTU : public TNamed {

 public:

  AliTRDCalDCSGTU();
  AliTRDCalDCSGTU(const char *name, const char *title);
  AliTRDCalDCSGTU(const AliTRDCalDCSGTU &);
  AliTRDCalDCSGTU& operator=(const AliTRDCalDCSGTU& sh);
  virtual ~AliTRDCalDCSGTU() { };

  Int_t   GetRunNumber()                              { return fRunNumber;                    }
  Int_t   GetSORFlag()                                { return fSORFlag;                      }
  Int_t   GetSerial()                                 { return fSerial;                       }
  Int_t   GetDNR()                                    { return fDNR;                          }

  void    SetRunNumber(Int_t rn)                      { fRunNumber = rn;                      }
  void    SetSORFlag(Int_t fg)                        { fSORFlag = fg;                        }
  void    SetSerial(Int_t se)                         { fSerial = se;                         }
  void    SetDNR(Int_t dn)                            { fDNR = dn;                            }

  TObjArray* GetSegmentArray() const                  { return fSegmentsArr;                  }
  void SetSegmentArray(TObjArray *sa)                 { fSegmentsArr = sa;                    }

  AliTRDCalDCSGTUTgu* GetTgu() const                  { return fTgu;                          }
  void SetTgu(AliTRDCalDCSGTUTgu* tg)                 { fTgu = tg;                            }

 protected:
  Int_t   fRunNumber;
  Int_t   fSORFlag;
  Int_t   fSerial;
  Int_t   fDNR;

  TObjArray *fSegmentsArr;

  AliTRDCalDCSGTUTgu* fTgu;

  ClassDef(AliTRDCalDCSGTU,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
