#ifndef ALITOFTRACK_H
#define ALITOFTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////
//  TOF Reconstructed track
//  AliTOFTrack  class                        
//  (see implementation file for details)
//                               
//-- Authors: Bologna-ITEP-Salerno Group
/////////////////////////////////////////////////////////////


#include "TObject.h"

//_______________________________________________________
class AliTOFTrack : public TObject{

public:
  AliTOFTrack();
  AliTOFTrack(Int_t track, Float_t vtxMom, Int_t pdgcode, Float_t tpcTrackLen, Float_t* tpcpos, Float_t* tpcmom, Float_t* trdpos, Int_t pad, Int_t matchflag, Float_t length, Float_t tof, Float_t massTof);
  ~AliTOFTrack(){};
  void SetTrack(Int_t track, Float_t vtxMom, Int_t pdgcode, Float_t tpcTrackLen, Float_t* tpcpos, Float_t* tpcmom, Float_t* trdpos);
  void SetTrack(Int_t track)   {fTrack=track;}   
  void SetP(Float_t mom)       {fP=mom;}
  void SetPdgCode(Int_t pdgCode) {fPdgCode=pdgCode;}   
  void SetlTPC(Float_t lTPC)   {flTPC=lTPC;}
  void SetRxTPC(Float_t rxTPC) {fRxTPC=rxTPC;}
  void SetRyTPC(Float_t ryTPC) {fRyTPC=ryTPC;}
  void SetRzTPC(Float_t rzTPC) {fRzTPC=rzTPC;}
  void SetPxTPC(Float_t pxTPC) {fPxTPC=pxTPC;}
  void SetPyTPC(Float_t pyTPC) {fPyTPC=pyTPC;}
  void SetPzTPC(Float_t pzTPC) {fPzTPC=pzTPC;}
  void SetRxTRD(Float_t rxTRD) {fRxTRD=rxTRD;}
  void SetRyTRD(Float_t ryTRD) {fRyTRD=ryTRD;}
  void SetPixel(Int_t pad)     {fPad=pad;}
  void SetMatching(Int_t matching) {fMatching=matching;}
  void SetLength(Float_t length)   {fLength=length;}
  void SetTof(Float_t tof)     {fTof=tof;}
  void SetMassTOF(Float_t massTOF) {fMassTOF=massTOF;}

  Int_t   GetTrack()   const {return fTrack;}
  Float_t GetP()       const {return fP;}
  Int_t   GetPdgCode() const {return fPdgCode;}   
  Float_t GetlTPC()    const {return flTPC;}
  Float_t GetRxTPC()   const {return fRxTPC;}
  Float_t GetRyTPC()   const {return fRyTPC;}
  Float_t GetRzTPC()   const {return fRzTPC;}
  Float_t GetPxTPC()   const {return fPxTPC;}
  Float_t GetPyTPC()   const {return fPyTPC;}
  Float_t GetPzTPC()   const {return fPzTPC;}
  Float_t GetRxTRD()   const {return fRxTRD;}
  Float_t GetRyTRD()   const {return fRyTRD;}
  Int_t   GetPad()     const {return fPad;}
  Int_t   GetMatching()const {return fMatching;}
  Float_t GetLength()  const {return fLength;}
  Float_t GetTof()     const {return fTof;}
  Float_t GetMassTOF() const {return fMassTOF;}

private:
  Int_t    fTrack;    // track number
  Float_t  fP;        // vertex momentum
  Int_t    fPdgCode;  // Geant code of particle
  Float_t  flTPC;     // length to TPC
  Float_t  fRxTPC;    // x-coordinate on TPC
  Float_t  fRyTPC;    // y-coordinate on TPC
  Float_t  fRzTPC;    // z-coordinate on TPC
  Float_t  fPxTPC;    // x-momentum on TPC
  Float_t  fPyTPC;    // y-momentum on TPC
  Float_t  fPzTPC;    // z-momentum on TPC
  Float_t  fRxTRD;    // x-coordinate on the last layer of TRD
  Float_t  fRyTRD;    // y-coordinate on the last layer of TRD
//  Float_t  fTof;    // Time of Flight [ns] smearing with RPC resolution
  Int_t    fPad  ;    // pad number
  Int_t    fMatching; // Index of TPC track - TOF pixel matching
  Float_t  fLength  ; // Track length [cm] from the origin to the TOF [cm]
  Float_t  fTof;      // Time [ns] determined by pixel matched with the track
  Float_t  fMassTOF;  // Mass [GeV] determined by fTOF,fLength,fPx,...

  ClassDef(AliTOFTrack,1)   // TOF Reconstructed track
};

#endif /* ALITOFTRACK_H */
