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
  AliTRDtrendValue(const Char_t *n, const Char_t *t);
  AliTRDtrendValue(const AliTRDtrendValue &ref);

  Double_t    GetVal() const { return fValue;}
  Double_t    GetErr() const { return fSigma;}
  void        Print(Option_t *o="") const;  // *MENU*
  void        Set(Double_t v, Double_t ve) { fValue=v; fSigma = ve;}

protected: // only manager can fill these info !! 
  AliTRDtrendValue& operator/=(const AliTRDtrendValue &n);
  const char* GetAlarmMessage(Int_t ns) const;
  const char* GetClassName() const;
  const char* GetValueName() const;
  const char* GetResponsible() const;
  const char* GetNotifiable(Int_t in) const;
  void        SetAlarm(Int_t level, Char_t *m);
  void        SetResponsible(const Char_t *name, const Char_t *mail, Option_t *opt="");
  void        SetNotifiable(const Char_t *name, const Char_t *mail);

private:
  const AliTRDtrendValue& operator=(const AliTRDtrendValue &ref);

  Double_t      fValue;     // value; mean for reference, current for running
  Double_t      fSigma;     // error; sigma for reference, current for running
//  Char_t        *fAlarmMessage[kNlevels]; // list of alarm messages
  TNamed       *fResponsible; // responsible person
  TNamed       *fNotifiable[kNnotifiable]; // also notify these persons

  ClassDef(AliTRDtrendValue, 1) // TRD trending value representation
};

#endif

