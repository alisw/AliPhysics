#ifndef ALIESDV0MI_H
#define ALIESDV0MI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//                          ESD V0 Vertex Class
//          This class is part of the Event Summary Data set of classes
//    Origin: Marian Ivanov marian.ivanov@cern.ch
//-------------------------------------------------------------------------

#include <AliESDv0.h>
#include "AliExternalTrackParam.h"
#include <TPDGCode.h>

class AliESDtrack;

class AliESDV0MI :  public AliESDv0 {
public:
  AliESDV0MI();             //constructor
  //
  void SetP(const AliExternalTrackParam & paramp); 
  void SetM(const AliExternalTrackParam & paramd);
  void UpdatePID(Double_t pidp[5], Double_t pidm[5]);
  Float_t GetEffMass(UInt_t p1, UInt_t p2);
  Float_t GetProb(UInt_t p1, UInt_t p2);
  void Update(Float_t vertex[3]);            //update
  void SetID(Int_t id){fID =id;}
  Int_t GetID(){ return fID;}
  AliExternalTrackParam fParamP;
  AliExternalTrackParam fParamM;
  Float_t        fRP[5];         // combined pid positive
  Float_t        fRM[5];         // combined pid positive
  Int_t          fID;
  Int_t          fLab[2];     //MC label of the partecle
  Int_t          fIndex[2];   //reconstructed labels of the tracks
  //
  //  
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  //
  Double_t       fPP[3];    //momentum  positive   - according approx at DCA
  Double_t       fPM[3];    //momentum negative
  //
  Double_t       fXr[3];      //rec. position according helix
  Double_t       fAngle[3];   //three angles
  Double_t       fRr;         //rec position of the vertex 
  Int_t          fStatus;       //status 
  Int_t          fRow0;         // critical layer
  Int_t          fOrder[3]; //order of the vertex 
  //  quality information
  Double_t       fDistNorm; //normalized  DCA
  Double_t       fDistSigma; //sigma of distance
  Float_t        fChi2Before;   //chi2 of the tracks before V0
  Float_t        fNBefore;      // number of possible points before V0
  Float_t        fChi2After;   // chi2 of the tracks after V0
  Float_t        fNAfter;      // number of possible points after V0
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full



  ClassDef(AliESDV0MI,1)      // ESD V0 vertex
};




#endif

/*

class AliITSRecV0Info: public AliV0vertex {
  friend class AliITStrackerMI;
protected:
  AliITSRecV0Info();
  AliITSRecV0Info(const AliITStrackV2 &neg, const AliITStrackV2 &pos);
  void Update(Float_t vertex[3], Float_t mass1, Float_t mass2);
  Double_t       fDist1;    //info about closest distance according closest MC - linear DCA
  Double_t       fDist2;    //info about closest distance parabolic DCA
  Double_t       fDistNorm; //normalized  DCA
  Double_t       fDistSigma; //sigma of distance
  Double_t       fInvMass;  //reconstructed invariant mass -
  //
  Double_t       fPdr[3];    //momentum at vertex daughter  - according approx at DCA
  Double_t       fXr[3];     //rec. position according helix
  //
  Double_t       fPm[3];    //momentum at the vertex mother
  Double_t       fAngle[3]; //three angles
  Double_t       fRr;       // rec position of the vertex 
  Int_t          fLab[2];   //MC label of the particle
  Int_t          fIndex[2]; //indexes of the tracks
  Int_t          fOrder[3]; //order of the vertex 
  Float_t        fChi2Before;   //chi2 of the tracks before V0
  Float_t        fNBefore;      // number of possible points before V0
  //
  Float_t        fChi2After;   // chi2 of the tracks after V0
  Float_t        fNAfter;      // number of possible points after V0
 
  //
  Float_t        fPointAngleFi; //point angle fi
  Float_t        fPointAngleTh; //point angle theta
  Float_t        fPointAngle;   //point angle full
  ClassDef(AliITSRecV0Info,1)  // container for  
};

*/
