#ifndef AliRICHMap_h
#define AliRICHMap_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHitMap.h"
#include <TMatrix.h>
#include <TMath.h>
#include <TClonesArray.h>

class AliRICHMap : public AliHitMap 
{    
public:
           AliRICHMap(TClonesArray *pDig);
  virtual ~AliRICHMap()                          {delete fMap;}
  void     FillHits();                                                                                             //virtual
  void     Clear(const char *)                   {fMap->Zero();}                                                   //virtual  
  void     DeleteHit(Int_t ix,Int_t iy)          {(*fMap)(ix,iy)=0;}                                               //virtual
  void     SetHit(Int_t ix,Int_t iy,Int_t idigit){(*fMap)(ix,iy)=idigit+1;}                                        //virtual  
  Int_t    GetHitIndex(Int_t ix,Int_t iy)   const{return (Int_t)TMath::Abs((*fMap)(ix, iy))-1;}                    //virtual
  TObject* GetHit(Int_t ix,Int_t iy)        const{Int_t idx=GetHitIndex(ix,iy);return(idx <0)?0:fDigits->At(idx);} //virtual
  void     FlagHit(Int_t ix,Int_t iy)            {(*fMap)(ix, iy)=-TMath::Abs((*fMap)(ix,iy));}                    //virtual
  Bool_t   ValidateHit(Int_t,Int_t)              {return 1;}                                                       //virtual
  inline   FlagType TestHit(Int_t ix,Int_t iy);                                                                    //virtual
  void     Print(const Option_t *)          const{fMap->Print();}                                                  
protected:
  TClonesArray *fDigits;                                 //List of digits
  Int_t         fNdigits;                                //Number of digits
  TMatrix      *fMap;                                    //hit map
  ClassDef(AliRICHMap,0)                                 //Implements map as TMatrix
};
//__________________________________________________________________________________________________
FlagType AliRICHMap::TestHit(Int_t padx,Int_t pady)
{
//Is there a hit for given pad?
  Int_t inf=(Int_t)(*fMap)(padx,pady);
  if(inf<0){//flaged as used
    return kUsed;
  }else if(inf==0){//no hit 
    return kEmpty;
  }else{//index of not yet used hit
    return kUnused;
  }
}//TestHit()  
#endif	
