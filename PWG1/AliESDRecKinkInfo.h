#ifndef ALIESDRECKINKINFO_H
#define ALIESDRECKINKINFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//////////////////////////////////////////////////////////////////////////////
//                          Class AliRecInfo                                //
//   collect together MC info and Rec info for comparison purposes 
//                                           - effieciency studies and so on//                                                                 //
//   marian.ivanov@cern.ch                                                  //
//////////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "AliESDRecInfo.h"



class AliESDRecKinkInfo: public TObject {
friend class  AliRecInfoMaker;
public:
  void Update();
protected:
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  AliESDkink     fKink;    //kink
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Double_t       fMinR;     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fStatus;       //status -tracks 
  Int_t          fRecStatus;    //kink -status- 0 - not found  1-good -  fake
  Int_t          fMultiple;     // how many times was kink reconstructed
  Int_t          fKinkMultiple; // how many times was kink reconstructed
  ClassDef(AliESDRecKinkInfo,1)   // container for  
};

#endif
