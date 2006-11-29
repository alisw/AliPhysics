#ifndef AliHMPIDHit_h
#define AliHMPIDHit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliHit.h>           //base class
#include <TVector3.h>         //ctor

class AliHMPIDHit : public AliHit //   TObject-AliHit-AliHMPIDHit
{
public:
  AliHMPIDHit(                                                                             ):AliHit(     ),fCh(-1),fPid(-1 ),fE(-1),fLorsX(-1),fLorsY(-1) {} //default ctor
  AliHMPIDHit(Int_t c,Float_t e,Int_t pid,Int_t tid,Float_t xl,Float_t yl,const TVector3 &p):AliHit(0,tid),fCh(c ),fPid(pid),fE(e ),fLorsX(xl),fLorsY(yl) {fX=p.X();fY=p.Y();fZ=p.Z();}           
  AliHMPIDHit(Int_t c,Float_t e,Int_t pid,Int_t tid,Float_t xl,Float_t yl                  ):              fCh(c ),fPid(pid),fE(e ),fLorsX(xl),fLorsY(yl) {fTrack=tid;}           
  virtual ~AliHMPIDHit()                                                                                                {}
//framework part
  void     Print(Option_t *option="")const;                                    //from TObject to print current status
//private part  
  Int_t   Ch    ()const{return fCh;                                    }       //Chamber
  Float_t E     ()const{return fE;                                     }       //Eloss for MIP hit or Etot for photon hit, [GeV]
  Float_t LorsX ()const{return fLorsX;                                 }       //hit X position in LORS, [cm]
  Float_t LorsY ()const{return fLorsY;                                 }       //hit Y position in LORS, [cm]
  Int_t   Pid   ()const{return fPid;                                   }       //PID
  Int_t   Tid   ()const{return fTrack;                                 }       //TID
         
protected:                                                                     //AliHit has fTid,fX,fY,fZ 
  Int_t    fCh;                                                                //Chamber
  Int_t    fPid;                                                               //PID
  Float_t  fE;                                                                 //Eloss for MIP or Etot for photon [GeV]
  Float_t  fLorsX;                                                             //hit X position in chamber LORS, [cm]
  Float_t  fLorsY;                                                             //hit Y position in chamber LORS, [cm]
  ClassDef(AliHMPIDHit,4)                                                       //HMPID hit class 
};//class AliHMPIDhit

typedef AliHMPIDHit AliRICHHit; // for backward compatibility
    
#endif
