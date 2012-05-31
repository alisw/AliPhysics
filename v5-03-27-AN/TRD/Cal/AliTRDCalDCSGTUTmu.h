#ifndef ALITRDCALDCSGTUTMU_H
#define ALITRDCALDCSGTUTMU_H
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
class AliTRDCalDCSGTUBoardInfo;

class AliTRDCalDCSGTUTmu : public TNamed {

 public:

  AliTRDCalDCSGTUTmu();
  AliTRDCalDCSGTUTmu(const char *name, const char *title);
  AliTRDCalDCSGTUTmu(const AliTRDCalDCSGTUTmu&);	
  AliTRDCalDCSGTUTmu &operator=(const AliTRDCalDCSGTUTmu &sh);
  virtual ~AliTRDCalDCSGTUTmu() { };

  TString GetLinkMask() const                         { return fLinkMask;                     }
  Int_t   GetId() const                               { return fId;                           }
  Int_t   GetPatternGeneratorEnable() const           { return fPatternGeneratorEnable;       }
  Int_t   GetPatternGeneratorDataWords() const        { return fPatternGeneratorDataWords;    }
  Int_t   GetPatternGeneratorTrackletWordsl() const   { return fPatternGeneratorTrackletWords;}

  void    SetLinkMask(TString lm)                     { fLinkMask = lm;                       }
  void    SetId(Int_t id)                             { fId = id;                             }
  void    SetPatternGeneratorEnable(Int_t pe)         { fPatternGeneratorEnable = pe;         }
  void    SetPatternGeneratorDataWords(Int_t pw)      { fPatternGeneratorDataWords = pw;      }
  void    SetPatternGeneratorTrackletWords(Int_t pt)  { fPatternGeneratorTrackletWords = pt;  }

  AliTRDCalDCSGTUBoardInfo* GetBoardInfo() const      { return fBoardInfo;                    }
  void SetBoardInfo(AliTRDCalDCSGTUBoardInfo * const bi) { fBoardInfo = bi;                   }

  protected:
  TString fLinkMask; // value of the attribute named value within the linkmask tag
  Int_t   fId; // the number of the tmu within the segment
  Int_t   fPatternGeneratorEnable; // value of the attribute named enable within the pattern_generator tag
  Int_t   fPatternGeneratorDataWords; // value of the attribute named datawords within the pattern_generator tag
  Int_t   fPatternGeneratorTrackletWords; // value of the attribute named trackletwords within the pattern_generator tag

  AliTRDCalDCSGTUBoardInfo *fBoardInfo; // This contains the board information for this tmu

  ClassDef(AliTRDCalDCSGTUTmu,1)      //  TRD calibration class for TRD GTU parameters

};
#endif
