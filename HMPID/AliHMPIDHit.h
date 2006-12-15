#ifndef AliHMPIDHit_h
#define AliHMPIDHit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliHit.h>           //base class
#include <TVector3.h>         //ctor
#include <TRandom.h>          //QdcTot()

class AliHMPIDHit : public AliHit //   TObject-AliHit-AliHMPIDHit
{
public:
  AliHMPIDHit(                                                                             ):AliHit(     ),fCh(-1),fPid(-1  ),fQ(-1),fLorsX(-1),fLorsY(-1) {} //default ctor
  AliHMPIDHit(Int_t c,Float_t e,Int_t pid,Int_t tid,Float_t xl,Float_t yl,const TVector3 &p):AliHit(0,tid),fCh(c ),fPid(pid ),fQ(-1),fLorsX(xl),fLorsY(yl) {QdcTot(e);fX=p.X();fY=p.Y();fZ=p.Z();}           
  AliHMPIDHit(Int_t c,Int_t q  ,                    Float_t xl,Float_t yl                  ):              fCh(c ),fPid(2212),fQ(q ),fLorsX(xl),fLorsY(yl) {fTrack=Int_t(1000*gRandom->Rndm());}//manual ctor           
  virtual ~AliHMPIDHit()                                                                                                {}
//framework part
  void     Print(Option_t *option="")const;                                    //from TObject to print current status
//private part  
         Int_t   Ch    (         )const{return fCh;                                    }       //Chamber
         Float_t LorsX (         )const{return fLorsX;                                 }       //hit X position in LORS, [cm]
         Float_t LorsY (         )const{return fLorsY;                                 }       //hit Y position in LORS, [cm]
         Int_t   Pid   (         )const{return fPid;                                   }       //PID
         Float_t Q     (         )const{return fQ;                                     }       //Eloss for MIP hit or Etot for photon hit, [GeV]
  inline Float_t QdcTot(Float_t e);                                                            //calculate total charge of the hit          
         Int_t   Tid   (         )const{return fTrack;                                 }       //TID
         
protected:                                                                     //AliHit has fTid,fX,fY,fZ 
  Int_t    fCh;                                                                //Chamber
  Int_t    fPid;                                                               //PID
  Float_t  fQ;                                                                 //total charge [QDC]
  Float_t  fLorsX;                                                             //hit X position in chamber LORS, [cm]
  Float_t  fLorsY;                                                             //hit Y position in chamber LORS, [cm]
  ClassDef(AliHMPIDHit,5)                                                      //HMPID hit class 
};//class AliHMPIDhit
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
Float_t AliHMPIDHit::QdcTot(Float_t e)
{
// Samples total charge of the hit
// Arguments: e- hit energy [GeV]; for mip Eloss for photon Etot   
//   Returns: total QDC
  Float_t  x=(LorsX() > 66.6)? LorsX()-66.6:LorsX();                                                 //sagita is for PC (0-64) and not for chamber   
  Float_t  qdcEle=34.06311+0.2337070*x+5.807476e-3*x*x-2.956471e-04*x*x*x+2.310001e-06*x*x*x*x;      //reparametrised from DiMauro
  
  Int_t iNele=Int_t(e/26e-9);  if(iNele<1) iNele=1;                                                  //number of electrons created by hit
  for(Int_t i=1;i<=iNele;i++) fQ-=qdcEle*TMath::Log(gRandom->Rndm()+1e-6);                           //1e-6 is a protection against 0 from rndm  
  return fQ;
}  
  
  
  
typedef AliHMPIDHit AliRICHHit; // for backward compatibility
    
#endif
