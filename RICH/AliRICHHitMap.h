#ifndef AliRICHHitMap_H
#define AliRICHHitMap_H

#include "AliRICH.h"
#include "TArrayI.h"
typedef enum {empty, used, unused} Flag_t;
const Int_t kMaxNpadx=1200, kMaxNpady=1200;

class AliRICHHitMap :
public TObject {
 public:
    virtual  void  FillHits()                                      =0;
    virtual  void  Clear()                                         =0;
    virtual  void  SetHit(Int_t ix, Int_t iy, Int_t idigit)        =0;
    virtual  void  DeleteHit(Int_t ix, Int_t iy)                   =0;
    virtual Int_t  GetHitIndex(Int_t ix, Int_t iy)                 =0;
    virtual TObject * GetHit(Int_t ix, Int_t iy)                   =0;
    virtual void   FlagHit(Int_t ix, Int_t iy)                     =0;    
    virtual Flag_t TestHit(Int_t ix, Int_t iy)                     =0;
    
    ClassDef(AliRICHHitMap,1) //virtual base class for muon HitMap
};

class AliRICHHitMapA1 :
public AliRICHHitMap 
{
 private:
    AliRICHsegmentation *fSegmentation;
    Int_t fNpx;
    Int_t fNpy;
    TObjArray *fDigits;
    Int_t fNdigits;
    Int_t *fHitMap;
    Int_t fMaxIndex;
    
 public:
    AliRICHHitMapA1(AliRICHsegmentation *seg, TObjArray *dig);
    virtual ~AliRICHHitMapA1();
    virtual  void  FillHits();
    virtual  void  Clear();    
    virtual  void  SetHit(Int_t ix, Int_t iy, Int_t idigit);
    virtual  void  DeleteHit(Int_t ix, Int_t iy);
    virtual Int_t  GetHitIndex(Int_t ix, Int_t iy);
    virtual TObject*  GetHit(Int_t ix, Int_t);
    virtual  void  FlagHit(Int_t ix, Int_t iy);    
    virtual Flag_t TestHit(Int_t ix, Int_t iy);
    private:
    Int_t CheckedIndex(Int_t ix, Int_t iy);
    ClassDef(AliRICHHitMapA1,1) // Implements HitMap as a 2-dim array
};
#endif	


