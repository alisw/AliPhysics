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
           AliRICHHit()                                                     :AliHit(     ),fChamber(-1),fEloss(-1) {fInX3.SetXYZ(0,0,0);fOutX3.SetXYZ(0,0,0);}
           AliRICHHit(Int_t c,Int_t tid,TVector3 in,TVector3 out,Double_t e):AliHit(0,tid),fChamber(c ),fEloss(e ) {fInX3=in; fOutX3=out; fX=out.X();fY=out.Y();fZ=out.Z();}
  virtual ~AliRICHHit()                                                                                            {}

  Int_t    C()                       const{return fChamber;}              //chamber number 
  Int_t    Chamber()                 const{return fChamber;}              //chamber number 
  Float_t  Eloss()                   const{return fEloss;}                //energy lost by track inside amplification gap  
  TVector3 InX3()                    const{return fInX3;}                 //track position at the faceplane of the gap 
  TVector3 OutX3()                   const{return fOutX3;}                //track position at the backplane of the gap 
  Double_t Length()                  const{return (fOutX3-fInX3).Mag();}  //track length inside the amplification gap
  void     Print(Option_t *option="")const;                               //virtual
protected:
  Int_t     fChamber;                      //chamber number
  Double_t  fEloss;                        //ionisation energy lost in GAP
  TVector3  fInX3;                         //position at the entrance of the GAP   
  TVector3  fOutX3;                        //position at the exit of the GAP
  ClassDef(AliRICHHit,2)                   //RICH hit class
};//class AliRICHhit
#endif
