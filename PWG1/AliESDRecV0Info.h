#ifndef ALIESDRECV0INFO_H
#define ALIESDRECV0INFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include "TObject.h"
#include "AliESDRecInfo.h"




/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


class AliESDRecV0Info: public TObject {
  friend class  AliRecInfoMaker;
public:
  void Update(Float_t vertex[3]);
protected:
  AliESDRecInfo  fT1;      //track1
  AliESDRecInfo  fT2;      //track2  
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fRs[2];     // minimum radius in rphi intersection
  Double_t       fDistMinR; // distance at minimal radius
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the partecle
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  Int_t          fV0Status;       // status of the kink
  AliV0*         fV0tpc;           // Vo information from reconsturction according TPC
  AliV0*         fV0its;           // Vo information from reconsturction according ITS
  AliV0*         fV0rec;           // V0 information form the reconstruction
  Int_t          fMultiple;     // how man times V0 was recostructed 
  Int_t          fV0Multiple;   // how man times was V0 reconstucted
  Int_t          fRecStatus;    // status form the reconstuction
  ClassDef(AliESDRecV0Info,2)   // container for  
};



#endif
