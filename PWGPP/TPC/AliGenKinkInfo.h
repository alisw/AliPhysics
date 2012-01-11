#ifndef ALIGENKINKINFO_H
#define ALIGENKINKINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliGenInfo                               //
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




class AliGenKinkInfo: public TObject {
public:
  AliGenKinkInfo();          //default cosntructor
  void    Update();          // put some derived info to special field 
  Float_t GetQt();           //
  AliMCInfo &  GetPlus()      {return fMCd;}
  AliMCInfo &  GetMinus()     {return fMCm;}
  void SetInfoDaughter(AliMCInfo &daughter) {fMCd=daughter;}
  void SetInfoMother(AliMCInfo &mother){fMCm=mother;}
private:
  AliMCInfo   fMCd;          //info about daughter particle - second particle for V0
  AliMCInfo   fMCm;          //info about mother particle   - first particle for V0
  Double_t    fMCDist1;      //info about closest distance according closest MC - linear DCA
  Double_t    fMCDist2;      //info about closest distance parabolic DCA
  //
  Double_t     fMCPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t     fMCPd[4];     //exact momentum from MC info
  Double_t     fMCX[3];      //exact position of the vertex
  Double_t     fMCXr[3];     //rec. position according helix
  //
  Double_t     fMCPm[3];     //momentum at the vertex mother
  Double_t     fMCAngle[3];  //three angels
  Double_t     fMCRr;        // rec position of the vertex 
  Double_t     fMCR;         //exact r position of the vertex
  Int_t        fPdg[2];      //pdg code of mother and daugter particles
  Int_t        fLab[2];      //MC label of the partecle
  ClassDef(AliGenKinkInfo,1) // container for  
};

#endif
