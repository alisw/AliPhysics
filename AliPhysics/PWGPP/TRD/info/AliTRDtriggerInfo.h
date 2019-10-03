#ifndef ALITRDTRIGGERINFO_H
#define ALITRDTRIGGERINFO_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Chamber Info Incapsulation                                            //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#define kTriggerListSize 1000

class TH1;
class TCollection;
class TObjArray;
class AliTRDtriggerInfo : public TObject
{
public:
  AliTRDtriggerInfo();
  virtual ~AliTRDtriggerInfo();

  Int_t       Add(const Char_t *trigger, Int_t nstat=1, Bool_t select=kFALSE);
  void        Draw(Option_t* option = "");  // *MENU*
  Int_t       GetNTriggers() const;
  const char* GetTrigger(Int_t it) const;
  Int_t       GetTrigger(const char *trigger) const;

  Bool_t      IsSelected(Int_t it) const            { return it>=0&&it<GetNTriggers()?fTriggerSel[it]:kFALSE;}
  Bool_t      IsSelected(const char* trigger) const { Int_t idx(GetTrigger(trigger)); return IsSelected(idx);}

  Long64_t    Merge(TCollection* list);
  void        Print(Option_t *o="") const;  // *MENU*
  void        SetSelectTrigger(Int_t it, Bool_t sel=kTRUE) { if(it>=0&&it<GetNTriggers()) fTriggerSel[it]=sel;}
private:
  AliTRDtriggerInfo(const AliTRDtriggerInfo &ref);
  const AliTRDtriggerInfo& operator=(const AliTRDtriggerInfo &ref);

  Bool_t      fTriggerSel[kTriggerListSize];  // trigger selection corresponding to fTriggerList
  Int_t       fTriggerStat[kTriggerListSize]; // trigger statistics corresponding to fTriggerList
  TObjArray  *fTriggerList;                   // list of trigger names
  TH1        *fHisto;                         //! graphic representation of trigger statistics

  ClassDef(AliTRDtriggerInfo, 1)              // trigger statistics
};

#endif

