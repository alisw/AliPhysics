#ifndef AliMFTHit_H
#define AliMFTHit_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Hit description for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TLorentzVector.h"
#include "TParticle.h"
#include "AliHit.h"
#include "AliRun.h"
#include "AliMC.h"

//====================================================================================================================================================

class TParticle;

class AliMFTHit : public AliHit {

public:

  AliMFTHit();
//   AliMFTHit(const AliMFTHit &h);
//   AliMFTHit& operator=(const AliMFTHit &h);
  
  virtual ~AliMFTHit() {}
  
  virtual void SetPlane(Int_t plane) { fPlane = plane; }
  virtual void SetDetElemID(Int_t detElemID) { fDetElemID = detElemID; }
  virtual void SetPosition(TLorentzVector &x) { fX =x.X(); fY =x.Y(); fZ =x.Z(); }
  virtual void SetTOF(Double_t time) { fTOF = time; }
  virtual void SetStatus(Int_t status) { fStatus  = status; }
  virtual void SetEloss(Double_t energy) { fEloss = energy; }
  virtual void SetMomentum(TLorentzVector &p) { fPx=p.Px(); fPy=p.Py(); fPz=p.Pz(); }
  
  virtual Int_t GetTrackStatus() const { return fStatus;  }
  virtual Int_t GetPlane()       const { return fPlane; }
  virtual Int_t GetDetElemID()   const { return fDetElemID; }
  virtual Double_t GetEloss()    const { return fEloss; }
  virtual Double_t GetTOF()      const { return fTOF; }
  
  virtual void GetPosition(Double_t &x,Double_t &y,Double_t &z) const { x=fX; y=fY; z=fZ; }
  virtual Double_t GetX() const { return fX; }
  virtual Double_t GetY() const { return fY; }
  virtual Double_t GetZ() const { return fZ; }

  virtual void GetMomentum(Double_t &px,Double_t &py,Double_t &pz) const { px=fPx; py=fPy; pz=fPz; }
  virtual Double_t GetPx() const { return fPx; }
  virtual Double_t GetPy() const { return fPy; }
  virtual Double_t GetPz() const { return fPz; }

  TParticle* GetParticle() const;

  Bool_t IsInside()      const { return (fStatus ==  1); }
  Bool_t IsEntering()    const { return (fStatus ==  2); }
  Bool_t IsExiting()     const { return (fStatus ==  4); }
  Bool_t IsOut()         const { return (fStatus ==  8); }
  Bool_t IsDisappeared() const { return (fStatus == 16); }
  Bool_t IsStopped()     const { return (fStatus == 32); }
  Bool_t IsAlive()       const { return (fStatus == 64); }

protected:

  Int_t   fStatus;   /* The track status flag. This flag indicates the track status
			at the time of creating this hit. 
			It is made up of the following 8 status bits from highest order to lowest order bits 0 :  
			IsTrackAlive(): IsTrackStop(): IsTrackDisappeared(): IsTrackOut(): IsTrackExiting(): IsTrackEntering(): IsTrackInside()     .
			See AliMC for a description of these functions. 
			If the function is true then the bit is set to one, otherwise it is zero. */
  
  Int_t    fPlane;      // Plane number
  Int_t    fDetElemID;  // Detection Element unique ID 
  Double_t fPx;         // PX of particle at the point of the hit
  Double_t fPy;         // PY of particle at the point of the hit
  Double_t fPz;         // PZ of particle at the point of the hit
  Double_t fEloss;      // Energy deposited in the current step
  Double_t fTOF;        // Time of flight at the point of the hit
  
  ClassDef(AliMFTHit,3)
    
};

//====================================================================================================================================================

#endif
