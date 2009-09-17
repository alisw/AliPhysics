#ifndef AliTRDCALDCSGTUSegment_H
#define AliTRDCALDCSGTUSegment_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUSegment.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"
#include "TObjArray.h"
#include "AliTRDCalDCSGTUBoardInfo.h"

class TString;

class AliTRDCalDCSGTUSegment : public TNamed {

 public:

  AliTRDCalDCSGTUSegment();
  AliTRDCalDCSGTUSegment(const char *name, const char *title);
  AliTRDCalDCSGTUSegment(const AliTRDCalDCSGTUSegment &);
  AliTRDCalDCSGTUSegment& operator=(const AliTRDCalDCSGTUSegment& sh);
  virtual ~AliTRDCalDCSGTUSegment() { };

  Int_t   GetId()                                     { return fId;                           }
  Int_t   GetFromRunNumber()                          { return fFromRunNumber;                }
  Int_t   GetFromSORFlag()                            { return fFromSORFlag;                  }
  Int_t   GetFromChild()                                  { return fChild;                        }

  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetFromRunNumber(Int_t rn)                      { fFromRunNumber = rn;                  }
  void    SetFromSORFlag(Int_t fg)                    { fFromSORFlag = fg;                    }
  void    SetFromChild(Int_t ch)                          { fChild = ch;                          }

  TObjArray* GetTmuArray() const                      { return fTmuArr;                       }
  void SetTmuArray(TObjArray *ta)                     { fTmuArr = ta;                         }

  TString GetSmuStackMask()                           { return fSmuStackMask;                 }
  Int_t   GetSmuTracklets()                           { return fSmuTracklets;                 }
  Int_t   GetSmuTracks()                              { return fSmuTracks;                    }
  Int_t   GetSmuIdelay()                              { return fSmuIdelay;                    }
  Int_t   GetSmuTriggerWindowL1Low()                  { return fSmuTriggerWindowL1Low;        }
  Int_t   GetSmuTriggerWindowL1High()                 { return fSmuTriggerWindowL1High;       }
  Int_t   GetSmuTriggerWindowL2Low()                  { return fSmuTriggerWindowL2Low;        }
  Int_t   GetSmuTriggerWindowL2High()                 { return fSmuTriggerWindowL2High;       }
  Int_t   GetSmuTtcEmulatorEnable()                   { return fSmuTtcEmulatorEnable;         }

  void    SetSmuStackMask(TString sm)                 { fSmuStackMask = sm;                   }
  void    SetSmuTracklets(Int_t ts)                   { fSmuTracklets = ts;                   }
  void    SetSmuTracks(Int_t tk)                      { fSmuTracks = tk;                      }
  void    SetSmuIdelay(Int_t id)                      { fSmuIdelay = id;                      }
  void    SetSmuTriggerWindowL1Low(Int_t ll)          { fSmuTriggerWindowL1Low = ll;          }
  void    SetSmuTriggerWindowL1High(Int_t lh)         { fSmuTriggerWindowL1High = lh;         }
  void    SetSmuTriggerWindowL2Low(Int_t ml)          { fSmuTriggerWindowL2Low = ml;          }
  void    SetSmuTriggerWindowL2High(Int_t mh)         { fSmuTriggerWindowL2High = mh;         }
  void    SetSmuTtcEmulatorEnable(Int_t te)           { fSmuTtcEmulatorEnable = te;           }

  AliTRDCalDCSGTUBoardInfo* GetSmuBoardInfo()         { return fSmuBoardInfo;                 }
  void SetSmuBoardInfo(AliTRDCalDCSGTUBoardInfo *bi)  { fSmuBoardInfo = bi;                   }

  protected:
  Int_t   fId;
  Int_t   fFromRunNumber;
  Int_t   fFromSORFlag;
  Int_t   fChild;

  TObjArray *fTmuArr;

  TString fSmuStackMask;
  Int_t   fSmuTracklets;
  Int_t   fSmuTracks;
  Int_t   fSmuIdelay;
  Int_t   fSmuTriggerWindowL1Low;
  Int_t   fSmuTriggerWindowL1High;
  Int_t   fSmuTriggerWindowL2Low;
  Int_t   fSmuTriggerWindowL2High;
  Int_t   fSmuTtcEmulatorEnable;

  AliTRDCalDCSGTUBoardInfo *fSmuBoardInfo;

  ClassDef(AliTRDCalDCSGTUSegment,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
