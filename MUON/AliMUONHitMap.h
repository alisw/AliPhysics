#ifndef AliMUONHitMap_H
#define AliMUONHitMap_H

#include "AliMUON.h"
#include "TArrayI.h"
typedef enum {empty, used, unused} Flag_t;
const Int_t kMaxNpadx=1200, kMaxNpady=1200;

class AliMUONHitMap :
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
    
    ClassDef(AliMUONHitMap,1) //virtual base class for muon HitMap
};

class AliMUONHitMapA1 :
public AliMUONHitMap 
{
 private:
    AliMUONsegmentation *fSegmentation;
    Int_t fNpx;
    Int_t fNpy;
    TObjArray *fDigits;
    Int_t fNdigits;
    Int_t *fHitMap;
    Int_t fMaxIndex;
    
 public:
    AliMUONHitMapA1(AliMUONsegmentation *seg, TObjArray *dig);
    virtual ~AliMUONHitMapA1();
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
    ClassDef(AliMUONHitMapA1,1) // Implements HitMap as a 2-dim array
};
#endif	


