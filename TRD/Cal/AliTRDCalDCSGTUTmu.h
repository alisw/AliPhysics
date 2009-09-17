#ifndef AliTRDCALDCSGTUTmu_H
#define AliTRDCALDCSGTUTmu_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTU.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "AliTRDCalDCSGTUBoardInfo.h"

class TString;

class AliTRDCalDCSGTUTmu : public TNamed {

 public:

  AliTRDCalDCSGTUTmu();
  AliTRDCalDCSGTUTmu(const char *name, const char *title);
  AliTRDCalDCSGTUTmu(const AliTRDCalDCSGTUTmu&);	
  AliTRDCalDCSGTUTmu &operator=(const AliTRDCalDCSGTUTmu &sh);
  virtual ~AliTRDCalDCSGTUTmu() { };

  TString GetLinkMask()                               { return fLinkMask;                     }
  Int_t   GetId()                                     { return fId;                           }
  Int_t   GetPatternGeneratorEnable()                 { return fPatternGeneratorEnable;       }
  Int_t   GetPatternGeneratorDataWords()              { return fPatternGeneratorDataWords;    }
  Int_t   GetPatternGeneratorTrackletWordsl()         { return fPatternGeneratorTrackletWords;}

  void    SetLinkMask(TString lm)                     { fLinkMask = lm;                       }
  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetPatternGeneratorEnable(Int_t pe)         { fPatternGeneratorEnable = pe;         }
  void    SetPatternGeneratorDataWords(Int_t pw)      { fPatternGeneratorDataWords = pw;      }
  void    SetPatternGeneratorTrackletWords(Int_t pt) { fPatternGeneratorTrackletWords = pt;  }

  AliTRDCalDCSGTUBoardInfo* GetBoardInfo()         { return fBoardInfo;                    }
  void SetBoardInfo(AliTRDCalDCSGTUBoardInfo *bi)  { fBoardInfo = bi;                      }

  protected:
  TString fLinkMask;
  Int_t   fId;
  Int_t   fPatternGeneratorEnable;
  Int_t   fPatternGeneratorDataWords;
  Int_t   fPatternGeneratorTrackletWords;

  AliTRDCalDCSGTUBoardInfo *fBoardInfo;

  ClassDef(AliTRDCalDCSGTUTmu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
