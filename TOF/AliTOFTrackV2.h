#ifndef ALITOFTRACKV2_H
#define ALITOFTRACKV2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////
//  TOF Reconstructed track
//  AliTOFTrackV2  class                        
//  (see implementation file for details)
//                               
//-- Author: F. Pierella
/////////////////////////////////////////////////////////////


#include "TObject.h"

//_______________________________________________________
class AliTOFTrackV2 : public TObject{

 public:
  AliTOFTrackV2();
  AliTOFTrackV2(Int_t trackLabel, Int_t matchingStatus, Float_t tpcMom, Float_t dEdX, Float_t* tpcXYZ, Float_t* tpcPtPz, Float_t* trdXYZ, Float_t* trdPxPyPz);
  ~AliTOFTrackV2(){};

  void AddRecLength(Float_t deltaRecLength)              {fRecTrackLength+=deltaRecLength;}
  void UpdateTrack(Int_t tofDigitTrackLabel, Int_t matching, Float_t tof);
  void UpdateTrack(Int_t pdgCode, Float_t trackLength);

  void SetTrackLabel(Int_t trackLabel)                   {fTrackLabel=trackLabel;}   
  void SetTOFDigitTrackLabel(Int_t tofDigitTrackLabel)   {fTOFDigitTrackLabel=tofDigitTrackLabel;}
  void SetPTPC(Float_t tpcMom)                           {fPTPC=tpcMom;}
  void SetPdgCode(Int_t pdgCode)                         {fPdgCode=pdgCode;}   
  void SetdEdX(Float_t dEdX)                             {fdEdX=dEdX;} 
  void SetxTPC(Float_t xTPC)                             {fxTPC=xTPC;}
  void SetyTPC(Float_t yTPC)                             {fyTPC=yTPC;}
  void SetzTPC(Float_t zTPC)                             {fzTPC=zTPC;}
  void SetPtTPC(Float_t ptTPC)                           {fPtTPC=ptTPC;}
  void SetPzTPC(Float_t pzTPC)                           {fPzTPC=pzTPC;}
  void SetxTRD(Float_t xTRD)                             {fxTRD=xTRD;}
  void SetyTRD(Float_t yTRD)                             {fyTRD=yTRD;}
  void SetzTRD(Float_t zTRD)                             {fzTRD=zTRD;}
  void SetPxTRD(Float_t pxTRD)                           {fPxTRD=pxTRD;}
  void SetPyTRD(Float_t pyTRD)                           {fPyTRD=pyTRD;}
  void SetPzTRD(Float_t pzTRD)                           {fPzTRD=pzTRD;}
  void SetMatchingStatus(Int_t matching)                 {fMatchingStatus=matching;}
  void SetLength(Float_t length)                         {fLength=length;}
  void SetRecLength(Float_t rlength)                     {fRecTrackLength=rlength;}
  void SetTof(Float_t tof)                               {fTof=tof;}
  void SetMassTOF(Float_t massTOF)                       {fMassTOF=massTOF;}
  void SetVertex(Float_t xvtx,Float_t yvtx,Float_t zvtx) {fXRecVtx=xvtx; fYRecVtx=yvtx; fZRecVtx=zvtx;}
  void SetMomVertex(Float_t pxvtx,Float_t pyvtx,Float_t pzvtx) {fPxRecVtx=pxvtx; fPyRecVtx=pyvtx; fPzRecVtx=pzvtx;}

  Int_t    GetTrackLabel()         const {return fTrackLabel;}
  Int_t    GetTOFDigitTrackLabel() const {return fTOFDigitTrackLabel;}
  Float_t  GetPTPC()               const {return fPTPC;}
  Int_t    GetPdgCode()            const {return fPdgCode;}
  Float_t  GetdEdX()               const {return fdEdX;}
  Float_t  GetxTPC()               const {return fxTPC;}
  Float_t  GetyTPC()               const {return fyTPC;}
  Float_t  GetzTPC()               const {return fzTPC;}
  Float_t  GetPtTPC()              const {return fPtTPC;}
  Float_t  GetPzTPC()              const {return fPzTPC;}
  Float_t  GetxTRD()               const {return fxTRD;}
  Float_t  GetyTRD()               const {return fyTRD;}
  Float_t  GetzTRD()               const {return fzTRD;}
  Float_t  GetPxTRD()              const {return fPxTRD;}
  Float_t  GetPyTRD()              const {return fPyTRD;}
  Float_t  GetPzTRD()              const {return fPzTRD;}
  Int_t    GetMatchingStatus()     const {return fMatchingStatus;}
  Float_t  GetLength()             const {return fLength;}
  Float_t  GetRecLength()          const {return fRecTrackLength;}
  Float_t  GetTof()                const {return fTof;}
  Float_t  GetMassTOF()            const {return fMassTOF;}
  Float_t  GetXVtx()               const {return fXRecVtx;}
  Float_t  GetYVtx()               const {return fYRecVtx;}
  Float_t  GetZVtx()               const {return fZRecVtx;}
  Float_t  GetPxVtx()              const {return fPxRecVtx;}
  Float_t  GetPyVtx()              const {return fPyRecVtx;}
  Float_t  GetPzVtx()              const {return fPzRecVtx;}

 private:
  Int_t    fTrackLabel;         // track label (rt->GetLabel()) as coming from TPC reconstruction
  Int_t    fTOFDigitTrackLabel; // track label stored into the TOF digit
                                // assigned to the track
  Float_t  fPTPC;      // momentum as given by reconstruction in TPC
  Int_t    fPdgCode;   // PDG code of the particle (for MC events)
  Float_t  fdEdX;      // total amount of loss energy in TPC and ITS
  Float_t  fxTPC;      // x-coordinate on TPC
  Float_t  fyTPC;      // y-coordinate on TPC
  Float_t  fzTPC;      // z-coordinate on TPC
  Float_t  fPtTPC;     // pt at the end of TPC
  Float_t  fPzTPC;     // pz-momentum at the end of TPC
  Float_t  fxTRD;      // x-coordinate on the last layer of TRD
  Float_t  fyTRD;      // y-coordinate on the last layer of TRD
  Float_t  fzTRD;      // y-coordinate on the last layer of TRD
  Float_t  fPxTRD;     // x-momentum at the end of TRD
  Float_t  fPyTRD;     // y-momentum at the end of TRD
  Float_t  fPzTRD;     // z-momentum at the end of TRD
  Int_t    fMatchingStatus; // matching status (not only for MC events)
                            // see details in the implementation file
  Float_t  fLength  ; // Track length [cm] from the origin to the TOF [cm]
                      // GEANT track length to be compared with that coming from
                      // reconstruction in order to evaluate the resolution
  Float_t  fTof;      // Time [ns] determined by the TOF digit assigned to the track
  Float_t  fMassTOF;  // Mass [GeV] determined by fTOF,fLength, and reconstructed momentum in TPC
  Float_t  fXRecVtx;  // x component of the reconstructed vertex position MRF
  Float_t  fYRecVtx;  // y component of the reconstructed vertex position MRF
  Float_t  fZRecVtx;  // z component of the reconstructed vertex position MRF
  Float_t  fPxRecVtx;  // x component of the reconstructed vertex momentum MRF
  Float_t  fPyRecVtx;  // y component of the reconstructed vertex momentum MRF
  Float_t  fPzRecVtx;  // z component of the reconstructed vertex momentum MRF
  Float_t  fRecTrackLength; // reconstructed track length (coarse)

  ClassDef(AliTOFTrackV2,2)   // TOF Reconstructed track
};

#endif /* ALITOFTRACKV2_H */
