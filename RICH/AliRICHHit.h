#ifndef AliRICHHit_h
#define AliRICHHit_h
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <AliHit.h>
#include <TVector3.h>

//RICH hit container
class AliRICHHit : public AliHit
{
public:
  AliRICHHit()                                                                         :AliHit(     ),fCham(-1) ,fE(-1),fPid(-1 ){fInX3.SetXYZ(0,0,0);fOutX3.SetXYZ(0,0,0);}
  AliRICHHit(Int_t c,Int_t tid,TVector3 in,TVector3 out,Double_t e,Int_t pid)          :AliHit(0,tid),fCham(c ) ,fE(e) ,fPid(pid){fInX3=in; fOutX3=out; fX=out.X();fY=out.Y();fZ=out.Z();}
  AliRICHHit(Int_t tid,Double_t e,Int_t pad,Double_t x,Double_t y,Double_t z,Int_t pid):AliHit(0,tid),fCham(pad),fE(e) ,fPid(pid){fX=x;fY=y;fZ=z;}
           
  virtual ~AliRICHHit()                                                                                            {}
//framework part
  void     Print(Option_t *option="")const;                               //from TObject to print current status
//private part  
  Int_t    C      ()const{return fCham;}                //chamber number 
  Int_t    Chamber()const{return fCham;}                //chamber number 
  Int_t    Pad    ()const{return fCham;}                //absolute pad number, definition in AliRICHParam
  Float_t  Eloss  ()const{return fE;   }                //Eloss for MIP hit or Etot for photon hit
  TVector3 InX3   ()const{return fInX3;}                //track position at the faceplane of the gap 
  TVector3 OutX3  ()const{return fOutX3;}               //track position at the backplane of the gap 
  Double_t Length ()const{return (fOutX3-fInX3).Mag();} //track length inside the amplification gap
protected:
  Int_t     fCham;                         //chamber number or in future absolute pad number
  Double_t  fE;                            //Eloss for MIP or Etot for photon [GeV]
  Int_t     fPid;                          //PID of particle created this hit
  TVector3  fInX3;                         //position at the entrance of the GAP   
  TVector3  fOutX3;                        //position at the exit of the GAP
  ClassDef(AliRICHHit,3)                   //RICH hit class
};//class AliRICHhit
#endif
