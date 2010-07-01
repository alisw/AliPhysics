#ifndef AliHMPIDHit_h
#define AliHMPIDHit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//.
// HMPID base class to produce hits
//.
//.
#include <AliHit.h>           //base class
#include <TVector3.h>         //ctor
#include "AliHMPIDDigit.h"    //QdcTot() 

#include <TRandom.h>

class AliHMPIDHit : public AliHit //   TObject-AliHit-AliHMPIDHit
{
public:
  AliHMPIDHit(                                                                           ):AliHit(     ),fCh(-1),fPid(-1 ),fQ(-1),fLx(0),fLy(0),fT(0) {} //default ctor
  AliHMPIDHit(Int_t c,Float_t &e,Int_t pid,Int_t tid,Float_t x,Float_t y, Float_t time, const TVector3 &p):AliHit(0,tid),fCh(c ),fPid(pid),fQ(0 ),fLx(x),fLy(y),fT(time) {e=QdcTot(e,time);fX=p.X();fY=p.Y();fZ=p.Z();}
  AliHMPIDHit(Int_t c,Float_t &e,Int_t pid,Int_t tid,Float_t x,Float_t y, Float_t time                 ):AliHit(     ),fCh(c ),fPid(pid),fQ(0 ),fLx(x),fLy(y),fT(time){e=QdcTot(e,time);fTrack=tid;}//manual ctor 
  AliHMPIDHit(const AliHMPIDHit &h):AliHit(h),fCh(h.fCh),fPid(h.fPid),fQ(h.fQ),fLx(h.fLx),fLy(h.fLy),fT(h.fT) {}//copy ctor
  virtual ~AliHMPIDHit()                                                                                                                                 {}
//framework part
         void    Print(Option_t *opt="")const;                                                    //from TObject to print current status
         void    Draw (Option_t *opt="");                                                         //from TObject to Draw this hit
//private part  
         Int_t   Ch     (                               )const{return fCh;                                    }       //Chamber
         void    Hit2Sdi(TClonesArray *pSdiLst,Int_t n=1)const;                                                       //add sdigits of this hit to the list 
         Float_t LorsX  (                               )const{return fLx;                                    }       //hit X position in LORS, [cm]
         Float_t LorsY  (                               )const{return fLy;                                    }       //hit Y position in LORS, [cm]
         Float_t HitTime(                               )const{return fT;                                     }       //hit formation time, [sec]
         Int_t   Pid    (                               )const{return fPid;                                   }       //PID
         Float_t Q      (                               )const{return fQ;                                     }       //total charge, [QDC]
  inline Float_t QdcTot (Float_t e, Float_t time        );                                                            //calculate total charge of the hit          
         Int_t   Tid    (                               )const{return fTrack;                                 }       //TID
         void    SetQ   (Float_t q                      )     {fQ=q;                                          }       //for debugging...
protected:                                                                     //AliHit has fTrack,fX,fY,fZ 
  Int_t    fCh;                                                                //Chamber
  Int_t    fPid;                                                               //PID
  Float_t  fQ;                                                                 //total charge [QDC]
  Float_t  fLx;                                                                //hit X position in chamber LORS, [cm]
  Float_t  fLy;                                                                //hit Y position in chamber LORS, [cm]
  Float_t  fT;                                                                 //hit formation time, [sec] 
  ClassDef(AliHMPIDHit,5)                                                      //HMPID hit class 
};//class AliHMPIDhit
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
Float_t AliHMPIDHit::QdcTot(Float_t e, Float_t time)
{
// Samples total charge of the hit
// Arguments: e- hit energy [GeV] for mip Eloss for photon Etot   
//   Returns: total QDC
  Int_t pc,px,py;
  AliHMPIDParam::Lors2Pad(fLx,fLy,pc,px,py); 
  if(py<0) fQ=0;
 else {
  Float_t y=AliHMPIDParam::LorsY(pc,py);  
  fLy=((y-fLy)>0)?y-0.2:y+0.2;                                                                       //shift to the nearest anod wire   
  
  Float_t  x=(fLx > 66.6)? fLx-66.6:fLx;                                                             //sagita is for PC (0-64) and not for chamber   
  Float_t  qdcEle=34.06311+0.2337070*x+5.807476e-3*x*x-2.956471e-04*x*x*x+2.310001e-06*x*x*x*x;      //reparametrised from DiMauro
  
  Int_t iNele=Int_t(e/26e-9);  if(iNele<1) iNele = 1;                                                //number of electrons created by hit, if photon e=0 implies iNele=1
  fQ=0;
  for(Int_t i=1;i<=iNele;i++){
    Double_t rnd=gRandom->Rndm(); if(rnd==0) rnd=1e-12;                                              //1e-12 is a protection against 0 from rndm  
    fQ-=qdcEle*TMath::Log(rnd);                
  }
 }
  if(time>1.2e-6) fQ=0;
 
  return fQ;
}  
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
#endif
