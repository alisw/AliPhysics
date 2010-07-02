#ifndef ALITRDTRENDINGMANAGER_H
#define ALITRDTRENDINGMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ALITRDTRENDVALUE_H
#include "AliTRDtrendValue.h"
#endif

class TObjArray;
class AliTRDtrendingManager : public TObject
{
public:
  virtual ~AliTRDtrendingManager();
  void              AddValue(Char_t *class_name
                            ,Char_t *value_name
                            ,Char_t *value_title
                            ,Double_t limits[2*(AliTRDtrendValue::kNlevels+1)]
                            ,Char_t *messages[AliTRDtrendValue::kNlevels]
                            ,Char_t *responsible
                            ,Char_t *notifiables=NULL
                            );
  AliTRDtrendValue* GetValue(Char_t *class_name, Char_t *value_name);
  static AliTRDtrendingManager*	Instance();
  Bool_t            ModifyValue(Char_t *class_name
                            ,Char_t *value_name
                            ,Char_t *value_title
                            ,Double_t *limits=NULL
                            ,Char_t **messages=NULL
                            ,Char_t *responsible=NULL
                            ,Char_t *notifiables=NULL
                            );
  void              Print(Option_t *o="") const;
  void              Save();
  void              ResetRunRange(Int_t runStart, Int_t runStop) {fRunRange[0]=runStart; fRunRange[1]=runStop;}
  void              Terminate();

protected:
  AliTRDtrendingManager();
  AliTRDtrendingManager(const AliTRDtrendingManager& ref);
  AliTRDtrendingManager& operator=(const AliTRDtrendingManager& ref);

private:
  static Bool_t	                fgTerminated; // instance terminate flag
	static AliTRDtrendingManager*	fgInstance;	  // instance
  TObjArray        *fEntries;    // list of trending values
  AliTRDtrendValue *fValue;      // current loaded trend value
  Int_t             fRunRange[2];// valability range

  ClassDef(AliTRDtrendingManager, 0) // TRD trending Manager
};

#endif
