#ifndef ALITOFTRACK_H
#define ALITOFTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------------------//
//                                                                            //
//   Description: class for handling ESD extracted tracks for TOF matching.   //
//                                                                            //
//----------------------------------------------------------------------------//

#include "TMath.h"
#include "TVector2.h"

#include "AliKalmanTrack.h"

#include "AliCluster3D.h"
#include "AliTOFGeometry.h"

class TObject;

class AliESDtrack;

class AliTOFtrack : public AliKalmanTrack {

public:

   AliTOFtrack();
   AliTOFtrack(const AliTOFtrack& t);
   AliTOFtrack(const AliESDtrack& t);
   AliTOFtrack& operator=(const AliTOFtrack &/*source*/); // ass. op.

   Int_t    GetSector() const {
     return Int_t(TVector2::Phi_0_2pi(GetAlpha())/AliTOFGeometry::GetAlpha())%AliTOFGeometry::NSectors();}

   Int_t    GetSeedLabel() const { return fSeedLab; }
   Int_t    GetSeedIndex() const { return fSeedInd; }
   void     SetSeedLabel(Int_t lab) { fSeedLab=lab; }
   void     SetSeedIndex(Int_t ind) { fSeedInd=ind; }
  
   Int_t Compare(const TObject *o) const;

   Double_t GetYat(Double_t xk, Bool_t & skip) const;
   Double_t GetPredictedChi2(const AliCluster3D *c) const;
   Bool_t   PropagateTo(const AliCluster3D *c);
   Bool_t   PropagateTo(Double_t xr, Double_t x0=8.72, Double_t rho=5.86e-3);
   Bool_t   PropagateToInnerTOF();
   Bool_t   PropagateToInnerTOFold();
   Bool_t   Rotate(Double_t angle) {
     return AliExternalTrackParam::Rotate(GetAlpha()+angle);
   }

protected:
   Bool_t Update(const AliCluster */*c*/, Double_t /*chi2*/, Int_t /*idx*/) {
     return 0;
   }
   Double_t GetPredictedChi2(const AliCluster */*c*/) const {return 0.;}
   
   Int_t    fSeedInd;     // ESD seed track index  
   Int_t    fSeedLab;     // track label taken from seeding  

 private:

   void GetPropagationParameters(Double_t *param);
   
   ClassDef(AliTOFtrack,2) // TOF reconstructed tracks

};                     

#endif
