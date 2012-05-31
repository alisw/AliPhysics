#ifndef ALITRDCALDCSGTUTGU_H
#define ALITRDCALDCSGTUTGU_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUTgu.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;
class TObjArray;
class AliTRDCalDCSGTUBoardInfo;

class AliTRDCalDCSGTUTgu : public TNamed {

 public:

  AliTRDCalDCSGTUTgu();
  AliTRDCalDCSGTUTgu(const char *name, const char *title);
  AliTRDCalDCSGTUTgu(const AliTRDCalDCSGTUTgu&);
  AliTRDCalDCSGTUTgu &operator=(const AliTRDCalDCSGTUTgu &sh);
  virtual ~AliTRDCalDCSGTUTgu();

  Int_t   GetFromRunNumber() const                    { return fFromRunNum;                   }
  Int_t   GetFromSORFlag() const                      { return fFromSORFlag;                  }
  Int_t   GetFromChild() const                        { return fFromChild;                    }
  TString GetSegmentMask() const                      { return fSegmentMask;                  }
  TString GetBusyMask() const                         { return fBusyMask;                     }
  TString GetContribMask() const                      { return fContribMask;                  }

  void    SetFromRunNumber(Int_t rn)                  { fFromRunNum = rn;                     }
  void    SetFromSORFlag(Int_t fs)                    { fFromSORFlag = fs;                    }
  void    SetFromChild(Int_t ch)                      { fFromChild = ch;                      }
  void    SetSegmentMask(TString sm)                  { fSegmentMask = sm;                    }
  void    SetBusyMask(TString bm)                     { fBusyMask = bm;                       }
  void    SetContribMask(TString cm)                  { fContribMask = cm;                    }

  AliTRDCalDCSGTUBoardInfo* GetBoardInfo()            { return fBoardInfo;                    }
  void SetBoardInfo(AliTRDCalDCSGTUBoardInfo * const bi) { fBoardInfo = bi;                   }

  TObjArray* GetCtpOpcArray() const                   { return fCtpOpcArr;                    }
  void SetCtpOpcArray(TObjArray * const ca)           { fCtpOpcArr = ca;                      }

  protected:
  Int_t   fFromRunNum; // the run number from when this data was saved
  Int_t   fFromSORFlag; // a flag indicating wether this data was saved from the start(=1) or end(=2) of run
  Int_t   fFromChild; // value of the attribute named child within the from tag
  TString fSegmentMask; // value of the attribute named value within the segment tag
  TString fBusyMask; // value of the attribute named value within the busymask tag
  TString fContribMask; // value of the attribute named value within the contribmask tag

  AliTRDCalDCSGTUBoardInfo *fBoardInfo; // BoardInfo Object holding the information about the tgu

  TObjArray *fCtpOpcArr; // an array of AliTRDCalDCSGTUCtpOpc objects holding their configuration data
  
  ClassDef(AliTRDCalDCSGTUTgu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
