#ifndef AliTRDCALDCSGTUTgu_H
#define AliTRDCALDCSGTUTgu_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUTgu.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"
#include "AliTRDCalDCSGTUBoardInfo.h"

class TString;

class AliTRDCalDCSGTUTgu : public TNamed {

 public:

  AliTRDCalDCSGTUTgu();
  AliTRDCalDCSGTUTgu(const char *name, const char *title);
  AliTRDCalDCSGTUTgu(const AliTRDCalDCSGTUTgu&);
  AliTRDCalDCSGTUTgu &operator=(const AliTRDCalDCSGTUTgu &sh);
  virtual ~AliTRDCalDCSGTUTgu() { };
//   ~AliTRDCalDCSGTUTgu() { };

  Int_t   GetFromRunNumber()                          { return fFromRunNum;                   }
  Int_t   GetFromSORFlag()                            { return fFromSORFlag;                  }
  Int_t   GetFromChild()                              { return fFromChild;                    }
  TString GetSegmentMask()                            { return fSegmentMask;                  }
  TString GetBusyMask()                               { return fBusyMask;                     }
  TString GetContribMask()                            { return fContribMask;                  }

  void    SetFromRunNumber(Int_t rn)                  { fFromRunNum = rn;                     }
  void    SetFromSORFlag(Int_t fs)                    { fFromSORFlag = fs;                    }
  void    SetFromChild(Int_t ch)                      { fFromChild = ch;                      }
  void    SetSegmentMask(TString sm)                  { fSegmentMask = sm;                    }
  void    SetBusyMask(TString bm)                     { fBusyMask = bm;                       }
  void    SetContribMask(TString cm)                  { fContribMask = cm;                    }

  AliTRDCalDCSGTUBoardInfo* GetBoardInfo()            { return fBoardInfo;                    }
  void SetBoardInfo(AliTRDCalDCSGTUBoardInfo *bi)     { fBoardInfo = bi;                      }

  TObjArray* GetCtpOpcArray() const                   { return fCtpOpcArr;                    }
  void SetCtpOpcArray(TObjArray *ca)                  { fCtpOpcArr = ca;                      }

  protected:
  Int_t   fFromRunNum;
  Int_t   fFromSORFlag;
  Int_t   fFromChild;
  TString fSegmentMask;
  TString fBusyMask;
  TString fContribMask;

  AliTRDCalDCSGTUBoardInfo *fBoardInfo;

  TObjArray *fCtpOpcArr;
  
  ClassDef(AliTRDCalDCSGTUTgu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
