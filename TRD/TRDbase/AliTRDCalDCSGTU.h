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

class TString;
class TObjArray;
class AliTRDCalDCSGTUTgu;

class AliTRDCalDCSGTU : public TNamed {

 public:

  AliTRDCalDCSGTU();
  AliTRDCalDCSGTU(const char *name, const char *title);
  AliTRDCalDCSGTU(const AliTRDCalDCSGTU &);
  AliTRDCalDCSGTU& operator=(const AliTRDCalDCSGTU& sh);
  virtual ~AliTRDCalDCSGTU();

  Int_t   GetRunNumber() const                        { return fRunNumber;                    }
  Int_t   GetSORFlag() const                          { return fSORFlag;                      }
  Int_t   GetSerial() const                           { return fSerial;                       }
  Int_t   GetDNR() const                              { return fDNR;                          }

  void    SetRunNumber(Int_t rn)                      { fRunNumber = rn;                      }
  void    SetSORFlag(Int_t fg)                        { fSORFlag = fg;                        }
  void    SetSerial(Int_t se)                         { fSerial = se;                         }
  void    SetDNR(Int_t dn)                            { fDNR = dn;                            }

  TObjArray* GetSegmentArray() const                  { return fSegmentsArr;                  }
  void SetSegmentArray(TObjArray * const sa)          { fSegmentsArr = sa;                    }

  AliTRDCalDCSGTUTgu* GetTgu() const                  { return fTgu;                          }
  void SetTgu(AliTRDCalDCSGTUTgu * const tg)          { fTgu = tg;                            }

 protected:
  Int_t   fRunNumber; // contains the number of the run from when this data was saved
  Int_t   fSORFlag; // contains an int indicating whether it was the start(=1) or end(=2) of run
  Int_t   fSerial; // value of the tag named serial
  Int_t   fDNR; // (DNR=does not respond) this indicates whether the GTU responded correctly

  TObjArray *fSegmentsArr; // Contains an array of AliTRDCalDCSGTUSegment objects holding gtu configuration data

  AliTRDCalDCSGTUTgu* fTgu; // this points to an object containing tgu configuration data

  ClassDef(AliTRDCalDCSGTU,2)      //  TRD calibration class for TRD GTU parameters

};
#endif
