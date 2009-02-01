#ifndef ALIHMPIDTRACK_H
#define ALIHMPIDTRACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Class HMPID track matching in the common tracking framework
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "TMath.h"
#include "TVector2.h"

#include "AliKalmanTrack.h"

#include "AliCluster3D.h"
#include "AliHMPIDCluster.h"

class TObject;
class AliESDtrack;

class AliHMPIDtrack : public AliKalmanTrack {

public:
   AliHMPIDtrack();
   AliHMPIDtrack(const AliHMPIDtrack& t);
   AliHMPIDtrack(const AliESDtrack& t);
   AliHMPIDtrack& operator=(const AliHMPIDtrack &/*source*/); // ass. op.
   Double_t GetPredictedChi2(const AliCluster3D *c) const;                                
   Bool_t   PropagateTo(const AliCluster3D *c);
   Bool_t   PropagateTo(Double_t xr, Double_t x0 = 8.72, Double_t rho = 5.86e-3);                 //Use material definition as for TOF???
   void     Propagate(Double_t len,Double_t x[3],Double_t p[3],Double_t bz) const;                //HMPID method moved from AliExternalTrackParam
   Bool_t   PropagateToR(Double_t r,Double_t step);
   Bool_t   Rotate(Double_t alpha, Bool_t absolute);
   Int_t    GetProlongation(Double_t xk, Double_t &y, Double_t &z);
   Bool_t   Intersect(Double_t pnt[3], Double_t norm[3], Double_t bz) const;                      //HMPID method moved from AliExternalTrackParam
   Bool_t   Intersect(AliHMPIDtrack *pTrk,Double_t pnt[3], Double_t norm[3]) ;                      //just for test 
   Bool_t   Update(const AliHMPIDCluster *pClu, Double_t chi2, Int_t index);
              
protected:
   Bool_t   Update(const AliCluster */*c*/, Double_t /*chi2*/, Int_t /*idx*/) {return 0;}
   Double_t GetPredictedChi2(const AliCluster */*c*/) const {return 0.;}
   
 private:
   ClassDef(AliHMPIDtrack,0) // HMPID reconstructed tracks

};                     

#endif
