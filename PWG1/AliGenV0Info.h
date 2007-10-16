#ifndef ALIGENV0INFO_H
#define ALIGENV0INFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenV0Info                               //
//   collect together MC info for comparison purposes - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////
//
// Start of implementation of the class AliTPCdigitRow
//
////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include "AliMCInfo.h"

class TFile;
class AliRunLoader;
class AliStack;
class AliTPCParam;




class AliGenV0Info: public TObject {
public:
  AliGenV0Info();       //
  void Update(Float_t vertex[3]);       
  AliMCInfo &  GetPlus()      {return fMCd;}
  AliMCInfo &  GetMinus()     {return fMCm;}
  TParticle &  GetMopther()   {return fMotherP;}
  Double_t    GetMCDist1() const { return fMCDist1;}
  Double_t    GetMCDist2() const {return fMCDist2;}  
  const Double_t*  GetMCPdr() const {return fMCPdr;}
  const Double_t*  GetMCPd()  const {return fMCPd;}
  const Double_t*  GetMCX()  const {return fMCX;}
  //  const Double_t    fMCXr;
  //
//   Double_t     fMCPm[3];    
//   Double_t     fMCAngle[3]; 
//   Double_t     fMCRr;       
//   Double_t     fMCR;       
//   Int_t        fPdg[2];   
//   Int_t        fLab[2];   
//   //
//   Double_t       fInvMass;  
//   Float_t        fPointAngleFi;
//   Float_t        fPointAngleTh;
//   Float_t        fPointAngle;  

  void SetInfoP(AliMCInfo &plus) {fMCd=plus;}
  void SetInfoM(AliMCInfo &minus){fMCm=minus;}
  void SetMother(TParticle&mother){fMotherP=mother;}
private:
  AliMCInfo   fMCd;       //info about daughter particle - second particle for V0
  AliMCInfo   fMCm;       //info about mother particle   - first particle for V0
  TParticle   fMotherP;   //particle info about mother particle
  Double_t    fMCDist1;    //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;    //info about closest distance parabolic DCA
  //
  Double_t    fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t    fMCPd[4];     //exact momentum from MC info
  Double_t    fMCX[3];      //exact position of the vertex
  Double_t    fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];    //momentum at the vertex mother
  Double_t     fMCAngle[3]; //three angels
  Double_t     fMCRr;       // rec position of the vertex 
  Double_t     fMCR;        //exact r position of the vertex
  Int_t        fPdg[2];   //pdg code of mother and daugter particles
  Int_t        fLab[2];   //MC label of the partecle  
  //
  Double_t       fInvMass;  //reconstructed invariant mass -
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  //
  ClassDef(AliGenV0Info,1)  // container for  
};



#endif
