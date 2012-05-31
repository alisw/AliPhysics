#ifndef ALITRDCALDCSGTUSEGMENT_H
#define ALITRDCALDCSGTUSEGMENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDCalDCSGTUSegment.h 18952 2007-06-08 11:36:12Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for TRD GTU configuration parameters               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class TString;
class TObjArray;
class AliTRDCalDCSGTUBoardInfo;

class AliTRDCalDCSGTUSegment : public TNamed {

 public:

  AliTRDCalDCSGTUSegment();
  AliTRDCalDCSGTUSegment(const char *name, const char *title);
  AliTRDCalDCSGTUSegment(const AliTRDCalDCSGTUSegment &);
  AliTRDCalDCSGTUSegment& operator=(const AliTRDCalDCSGTUSegment& sh);
  virtual ~AliTRDCalDCSGTUSegment() { };

  Int_t   GetId() const                               { return fId;                           }
  Int_t   GetFromRunNumber() const                    { return fFromRunNumber;                }
  Int_t   GetFromSORFlag() const                      { return fFromSORFlag;                  }
  Int_t   GetFromChild() const                        { return fChild;                        }

  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetFromRunNumber(Int_t rn)                  { fFromRunNumber = rn;                  }
  void    SetFromSORFlag(Int_t fg)                    { fFromSORFlag = fg;                    }
  void    SetFromChild(Int_t ch)                      { fChild = ch;                          }

  TObjArray* GetTmuArray() const                      { return fTmuArr;                       }
  void SetTmuArray(TObjArray * const ta)              { fTmuArr = ta;                         }

  TString GetSmuStackMask() const                     { return fSmuStackMask;                 }
  Int_t   GetSmuTracklets() const                     { return fSmuTracklets;                 }
  Int_t   GetSmuTracks() const                        { return fSmuTracks;                    }
  Int_t   GetSmuIdelay() const                        { return fSmuIdelay;                    }
  Int_t   GetSmuTriggerWindowL1Low() const            { return fSmuTriggerWindowL1Low;        }
  Int_t   GetSmuTriggerWindowL1High() const           { return fSmuTriggerWindowL1High;       }
  Int_t   GetSmuTriggerWindowL2Low() const            { return fSmuTriggerWindowL2Low;        }
  Int_t   GetSmuTriggerWindowL2High() const           { return fSmuTriggerWindowL2High;       }
  Int_t   GetSmuTtcEmulatorEnable() const             { return fSmuTtcEmulatorEnable;         }

  void    SetSmuStackMask(TString sm)                 { fSmuStackMask = sm;                   }
  void    SetSmuTracklets(Int_t ts)                   { fSmuTracklets = ts;                   }
  void    SetSmuTracks(Int_t tk)                      { fSmuTracks = tk;                      }
  void    SetSmuIdelay(Int_t id)                      { fSmuIdelay = id;                      }
  void    SetSmuTriggerWindowL1Low(Int_t ll)          { fSmuTriggerWindowL1Low = ll;          }
  void    SetSmuTriggerWindowL1High(Int_t lh)         { fSmuTriggerWindowL1High = lh;         }
  void    SetSmuTriggerWindowL2Low(Int_t ml)          { fSmuTriggerWindowL2Low = ml;          }
  void    SetSmuTriggerWindowL2High(Int_t mh)         { fSmuTriggerWindowL2High = mh;         }
  void    SetSmuTtcEmulatorEnable(Int_t te)           { fSmuTtcEmulatorEnable = te;           }

  AliTRDCalDCSGTUBoardInfo* GetSmuBoardInfo() const   { return fSmuBoardInfo;                 }
  void SetSmuBoardInfo(AliTRDCalDCSGTUBoardInfo * const bi)  { fSmuBoardInfo = bi;            }

  protected:
  Int_t   fId; // this is the number of the segment
  Int_t   fFromRunNumber; // this is the run number from when this configuration data was saved
  Int_t   fFromSORFlag; // this indicates when the data was saved (1 = start of run and 2 = end)
  Int_t   fChild; // this comes from the value of the child attribute of the tag named from

  TObjArray *fTmuArr; // an array of objects holding the segment's tmu information

  TString fSmuStackMask; // value of the attribute named value within the stackmask tag
  Int_t   fSmuTracklets; // value of the attribute named send within the tracklets tag
  Int_t   fSmuTracks; // value of the attribute named send within the tracklets tag
  Int_t   fSmuIdelay; // value of the attribute named value within the idelay tag
  Int_t   fSmuTriggerWindowL1Low; // value of the attribute named l1_low within the trigger_window tag
  Int_t   fSmuTriggerWindowL1High; // value of the attribute named l1_high within the trigger_window tag
  Int_t   fSmuTriggerWindowL2Low; // value of the attribute named l2_low within the trigger_window tag
  Int_t   fSmuTriggerWindowL2High; // value of the attribute named l2_high within the trigger_window tag
  Int_t   fSmuTtcEmulatorEnable; // value of the attribute named enable within the ttc_emulator tag

  AliTRDCalDCSGTUBoardInfo *fSmuBoardInfo; // the boardinfo for the smu

  ClassDef(AliTRDCalDCSGTUSegment,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
