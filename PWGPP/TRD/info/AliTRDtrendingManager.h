#ifndef ALITRDTRENDINGMANAGER_H
#define ALITRDTRENDINGMANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Trend Value Manager                                                   //
//                                                                        //
//  Authors:                                                              //
//    Alexandru Bercuci <A.Bercuci@gsi.de>                                //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

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
  void              AddValue(const Char_t *name
                            ,Double_t mean,Double_t sigm=1
                            ,const Char_t *title=NULL
                            ,const Char_t *responsible=NULL
                            ,const Char_t *notifiables=NULL
                            ,Char_t **messages = NULL
                            );
  AliTRDtrendValue* GetValue(const Char_t *name);
  static AliTRDtrendingManager*	Instance();
  void              Load(const char *fn = "$ALICE_ROOT/PWGPP/TRD/data/TRD.Trend.root");
  Bool_t            MakeTrends(const char *fileList);
  Bool_t            ModifyValue(const Char_t *name
                            ,const Char_t *title
                            ,Double_t mean,Double_t sigm
                            ,Char_t **messages=NULL
                            ,const Char_t *responsible=NULL
                            ,const Char_t *notifiables=NULL
                            );
  void              Print(Option_t *o="") const;
//  void              ResetRunRange(Int_t runStart, Int_t runStop) {fRunRange[0]=runStart; fRunRange[1]=runStop;}
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
//  Int_t             fRunRange[2];// valability range

  ClassDef(AliTRDtrendingManager, 0) // TRD trending Manager
};

#endif
