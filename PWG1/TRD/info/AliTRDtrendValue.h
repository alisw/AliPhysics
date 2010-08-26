#ifndef ALITRDTRENDVALUE_H
#define ALITRDTRENDVALUE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Trend Value Incapsulation                                             //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TNamed
#include "TNamed.h"
#endif

class AliTRDtrendValue : public TNamed
{
friend class AliTRDtrendingManager; // allow easy access for Manager
public:  
  enum ETRDtrendValue{
    kNlevels = 4
   ,kNnotifiable = 10 
  };
  AliTRDtrendValue();
  AliTRDtrendValue(Char_t *n, Char_t *t);

  Double_t    Get() const { return fValue;}
  Int_t       GetAlarmLevel();
  const char* GetAlarmMessage() const;
  const char* GetClassName() const;
  const char* GetValueName() const;
  const char* GetResponsible(Char_t *n=NULL, Char_t *mail=NULL) const;
  const char* GetNotifiable(Int_t in, Char_t *n=NULL, Char_t *mail=NULL) const;
  Bool_t      IsAlarm() const { return fAlarmLevel>0;}
  void        Print(Option_t *o="") const;
  void        Set(Double_t v) { fValue=v; GetAlarmLevel();}

protected: // only manager can fill these info !! 
  void        SetAlarm(Int_t level, Char_t *m) {memcpy(fAlarmMessage[level], m, 1024*sizeof(Char_t));}
  void        SetLimits(Double_t l[2*(kNlevels+1)]) {memcpy(fLimits, l, 2*(kNlevels+1)*sizeof(Double_t));}
  void        SetResponsible(const Char_t *name, const Char_t *mail);
  void        SetNotifiable(const Char_t *name, const Char_t *mail);

private:
  AliTRDtrendValue(const AliTRDtrendValue &ref);
  const AliTRDtrendValue& operator=(const AliTRDtrendValue &ref);

  UChar_t       fAlarmLevel;// alarm level
  Double_t      fValue;     // current value
  Double_t      fLimits[2*(kNlevels+1)];// limits
  Char_t        fAlarmMessage[kNlevels][1024]; // list of alarm messages

  struct AliTRDtrendValueResponsible{
    AliTRDtrendValueResponsible(Char_t *name=NULL, Char_t *mail=NULL);
    Char_t fNameR[100]; // name of responsible 
    Char_t fMail[200]; // mail of responsible
  };
  AliTRDtrendValueResponsible fResponsible; // responsible person
  Int_t                       fNnotifiable; // number of persons to be notify
  AliTRDtrendValueResponsible fNotifiable[kNnotifiable]; //also notify these persons 

  ClassDef(AliTRDtrendValue, 0) // TRD trending value representation
};

#endif

