#ifndef ALIRICHHITMAP_H
#define ALIRICHHITMAP_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */


#include <TObject.h>


typedef enum {kEmpty, kUsed, kUnused} FlagType;
const Int_t kMaxNpadx=1200, kMaxNpady=1200;

class AliRICHHitMap :
public TObject {
 public:
    virtual  void    FillHits()                                      =0;
    virtual  void    Clear()                                         =0;
    virtual  void    SetHit(Int_t ix, Int_t iy, Int_t idigit)        =0;
    virtual  void    DeleteHit(Int_t ix, Int_t iy)                   =0;
    virtual Int_t    GetHitIndex(Int_t ix, Int_t iy)                 =0;
    virtual TObject* GetHit(Int_t ix, Int_t iy)                      =0;
    virtual void     FlagHit(Int_t ix, Int_t iy)                     =0;    
    virtual FlagType TestHit(Int_t ix, Int_t iy)                     =0;
    
    ClassDef(AliRICHHitMap,1) //virtual base class for muon HitMap
};
#endif	


