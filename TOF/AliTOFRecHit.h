#ifndef ALITOFRECHIT_H
#define ALITOFRECHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliHit.h"

/////////////////////////////////////////////////////////////
//                                            
//  Dummy Hit class for TOF Reconstruction         
//  AliTOFRecHit  class                        
//  (see implementation file for details)
//                               
//-- Authors: Bologna-ITEP-Salerno Group         
/////////////////////////////////////////////////////////////

 
class AliTOFRecHit : public AliHit {
  
public:
  AliTOFRecHit() {}
  AliTOFRecHit(Int_t shunt, Int_t track);
  AliTOFRecHit(const AliTOFRecHit & hit) ;
  virtual ~AliTOFRecHit() {}
       // setters and getters for AliTOFRecHit object

  void SetTrack(Int_t track)     {fTrack=track;}   
  void SetPdgCode(Int_t pdgCode) {fPdgCode=pdgCode;}   
  void SetX(Float_t X)           {fX=X;}
  void SetY(Float_t Y)           {fY=Y;}
  void SetZ(Float_t Z)           {fZ=Z;}
  void SetHit(Int_t track, Int_t pdgCode, Float_t* mrfpos, Float_t mom, Float_t vtxRadius, Int_t isFirstHit);
  void SetP(Float_t p)           {fP=p;}
  void SetVrho(Float_t vrho)     {fVrho=vrho;}
  void SetFirst(Int_t first)     {fFirst=first;}   
  void SetNoise(Int_t noise)     {fNoise=noise;}   
  void SetRmin(Float_t rmin)     {fRmin=rmin;}

  Int_t   GetTrack()   const {return fTrack;}
  Int_t   GetPdgCode() const {return fPdgCode;}
  // getters for fX, fY and fZ implemented in the mother class as X(), Y(), Z()
  Float_t GetP()       const {return fP;}
  Float_t GetVrho()    const {return fVrho;}
  Int_t   GetFirst()   const {return fFirst;}
  Int_t   GetNoise()   const {return fNoise;}
  Float_t GetRmin()    const {return fRmin;}


protected:
  Int_t    fTrack;    //track number of the particle that produced the hit
  Int_t    fPdgCode;  //GEANT code of the particle that produced the hit
  /* they are defined in the mother class
  Float_t  fX;        //x-coordinate of the hit 
  Float_t  fY;        //y-coordinate of the hit 
  Float_t  fZ;        //z-coordinate of the hit
  */
  Float_t  fP;        //momentum
  Float_t  fVrho;     //rho-coordinate of the Vertex
  Int_t    fFirst;    //=1 for the first hit of the track, =0 otherwise
  Int_t    fNoise;    //=1 for the noise hit (Rvtx>200 or second, ... hit), =0 otherwise
  Float_t  fRmin;     //distance to the nearest TOFhit (cm)

  ClassDef(AliTOFRecHit,1)  // Dummy Hit class for TOF Reconstruction
};

#endif /* ALITOFRECHIT_H */
