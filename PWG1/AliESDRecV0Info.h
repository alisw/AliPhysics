#ifndef ALIESDRECV0INFO_H
#define ALIESDRECV0INFO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



#include "TObject.h"
#include "AliESDRecInfo.h"
class AliESDVertex;
class AliKFParticle;



/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


class AliESDRecV0Info: public TObject {
  friend class  AliRecInfoMaker;
public:
  AliESDRecV0Info();
  void Reset();   
  void Update(Float_t vertex[3]);
  void UpdateKF(const AliESDVertex &vertex, Int_t pdg0, Int_t pdg1, Float_t mass);
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
  //
  Int_t          fV0Status;       // status of the V0
  AliV0*         fV0tpc;           // Vo information from reconsturction according TPC
  AliV0*         fV0its;           // Vo information from reconsturction according ITS
  AliV0*         fV0rec;           // V0 information form the reconstruction
  AliV0*         fV0recOff;        // V0 information form the reconstruction - OFFLINE
  Int_t          fMultiple;     // how man times V0 was recostructed 
  Int_t          fRecStatus;    // status form the reconstuction - 1 reconstructed - -1 fake
  Int_t          fV0MultipleOn;    // how man times was V0 reconstucted - onfly
  Int_t          fV0MultipleOff;   // how man times was V0 reconstucted - offline
  //
  // AliKF variables - variables to make a selection + resoluton study
  //
  Float_t        fKFrecChi2NC;     //  ONLINE V0 finder non constrained chi2  
  Float_t        fKFrecChi2C;      //  ONLINE V0 finder   constrained chi2 - prim vertex  
  Float_t        fKFrecChi2CM;     //  ONLINE V0 finder   constrained chi2 - prim vertex+mass 
  AliKFParticle* fKFRecNC;         //  non constrained  
  AliKFParticle* fKFRecC;          //  constrained vertex
  AliKFParticle* fKFRecCM;         //  constrained vertex+mass
  //
  Float_t        fKFrecOffChi2NC;  // OFFLINE V0 finder - non constrained chi2  
  Float_t        fKFrecOffChi2C;   // OFFLINE V0 finder -     constrained chi2 - prim vertex  
  Float_t        fKFrecOffChi2CM;  // OFFLINE V0 finder -     constrained chi2 - prim vertex+mass
  AliKFParticle* fKFOffRecNC;       //  non constrained  
  AliKFParticle* fKFOffRecC;        //  constrained vertex
  AliKFParticle* fKFOffRecCM;       //  constrained vertex+mass

 private:
  AliESDRecV0Info(const AliESDRecV0Info&); // Not implemented
  AliESDRecV0Info& operator=(const AliESDRecV0Info&); // Not implemented


  ClassDef(AliESDRecV0Info,2)   // container for  
};



#endif
