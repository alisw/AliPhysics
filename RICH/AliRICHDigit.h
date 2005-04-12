#ifndef AliRICHDigit_h
#define AliRICHDigit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliDigit.h>
#include "AliRICHParam.h"

class AliRICHDigit :public AliDigit
{
public:
  AliRICHDigit():AliDigit(),fCFM(0),fChamber(0),fPadX(0),fPadY(0),fQdc(-1){fTracks[0]=fTracks[1]=fTracks[2]=-1;}
  AliRICHDigit(Int_t c,TVector pad,Double_t q,Int_t cfm,Int_t tid0,Int_t tid1,Int_t tid2):fCFM(cfm)  
       {fPadX=(Int_t)pad[0];fPadY=(Int_t)pad[1];fQdc=q;fChamber=10*c+AliRICHParam::Pad2Sec(pad);fTracks[0]=tid0;fTracks[1]=tid1;fTracks[2]=tid2;}
  virtual ~AliRICHDigit() {;}
    
  Int_t    Compare(const TObject *pObj) const
                 {if(Id()==((AliRICHDigit*)pObj)->Id())return 0;else if(Id()>((AliRICHDigit*)pObj)->Id())return 1;else return -1;}  //virtual      
  virtual Bool_t   IsSortable()                 const{return kTRUE;}                              //sort interface
  virtual void     Print(Option_t *option="")   const;                                            //virtual
//private part  
  Int_t    ChFbMi()                     const{return fCFM;}                               //particle mixture for this digit
  Int_t    C()                          const{return fChamber/10;}                        //chamber number
  Int_t    S()                          const{return fChamber-(fChamber/10)*10;}          //sector number
  Int_t    X()                          const{return fPadX;}                              //x position of the pad
  Int_t    Y()                          const{return fPadY;}                              //y postion of the pad
  TVector  Pad()                        const{Float_t v[2]={fPadX,fPadY}; return TVector(2,v);}
  Int_t    Id()                         const{return fChamber*10000000+fPadX*1000+fPadY;} //absolute id of this pad
  Double_t Q()                          const{return fQdc;}                               //charge in terms of ADC channels
  void     AddTidOffset(Int_t offset) 
    {for (Int_t i=0; i<3; i++) if (fTracks[i]>0) fTracks[i]+=offset;};
protected:
  Int_t    fCFM;  //1000000*Ncerenkovs+1000*Nfeedbacks+Nmips  
  Int_t    fChamber;  //10*chamber number+ sector number 
  Int_t    fPadX;     //pad number along X
  Int_t    fPadY;     //pad number along Y
  Double_t fQdc;      //QDC value, fractions are permitted for summable procedure  
  ClassDef(AliRICHDigit,3) //RICH digit class       
};//class AliRICHDigit

#endif
