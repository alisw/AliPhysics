#ifndef ALIITSTRACKMI_H
#define ALIITSTRACKMI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                       ITS Track Class
//
//        Origin: Marian Ivanov, CERN, Marian.Ivanov@cern.ch 
//     dEdx analysis by: Boris Batyunya, JINR, Boris.Batiounia@cern.ch
//-------------------------------------------------------------------------


/*****************************************************************************
 *                          December 18, 2000                                *
 *  Internal view of the ITS track parametrisation as well as the order of   *
 *           track parameters are subject for possible changes !             *
 *  Use GetExternalParameters() and GetExternalCovariance() to access ITS    *
 *      track information regardless of its internal representation.         *
 * This formation is now fixed in the following way:                         *
 *      external param0:   local Y-coordinate of a track (cm)                *
 *      external param1:   local Z-coordinate of a track (cm)                *
 *      external param2:   local sine of the track momentum azimuthal angle  *
 *      external param3:   tangent of the track momentum dip angle           *
 *      external param4:   1/pt (1/(GeV/c))                                  *
 *****************************************************************************/

#include <AliKalmanTrack.h>

#include "AliITSrecoV2.h"
#include "AliITStrackV2.h"

class AliESDtrack;

//_____________________________________________________________________________
class AliITStrackMI : public AliITStrackV2 {
public:
  AliITStrackMI();
  AliITStrackMI(AliESDtrack& t,Bool_t c=kFALSE) throw (const Char_t *);
  AliITStrackMI(const AliITStrackMI& t);
  Int_t GetProlongationFast(Double_t alpha, Double_t xr,Double_t &y, Double_t &z);
  Bool_t UpdateMI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz, Double_t chi2, Int_t i);  
  Int_t CorrectForMaterial(Double_t d, Double_t x0=21.82);

  void UpdateESDtrack(ULong_t flags);

  void SetReconstructed(Bool_t sr=kTRUE){fReconstructed = sr;}  
  Bool_t GetReconstructed() const {return fReconstructed;}
  void SetChi2MIP(Int_t i,Float_t val){fChi2MIP[i]=val;}
  Float_t GetChi2MIP(Int_t i) const {return fChi2MIP[i];}  
  void IncrementNSkipped(){fNSkipped++;} // increment by 1 the # of skipped cls
  Float_t GetNSkipped() const {return fNSkipped;}
  void SetNSkipped(Float_t n) {fNSkipped=n;}
  void IncrementNUsed(){fNUsed++;} // increment by 1 the # of shared clusters
  Float_t GetNUsed() const {return fNUsed;}
  void SetNUsed(Float_t n) {fNUsed=n;}

  Int_t Compare(const TObject *o) const;
  Double_t GetCov33() const {return GetCovariance()[9];} // cov. matrix el. 3,3
  //Double_t GetCov44() const {return GetCovariance()[15];}// cov. matrix el. 4,4
  Float_t GetDy(Int_t i) const {return fDy[i];}
  Float_t GetDz(Int_t i) const {return fDz[i];}
  Float_t GetD(Int_t i) const {return fD[i];}
  Double_t GetD(Double_t x, Double_t y) const
    {return AliITStrackV2::GetD(x,y);}
  Float_t *GetDP() {return fD;}
  void SetD(Int_t i, Float_t d) {fD[i]=d;}
  Float_t GetDnorm(Int_t i) const {return fDnorm[i];}
  Float_t *GetDnormP() {return fDnorm;}
  void SetDnorm(Int_t i, Float_t d) {fDnorm[i]=d;}
  Float_t GetSigmaY(Int_t i) const {return fSigmaY[i];}
  Float_t GetSigmaZ(Int_t i) const {return fSigmaZ[i];}
  void SetSigmaY(Int_t i, Float_t s) {fSigmaY[i]=s;}
  void SetSigmaZ(Int_t i, Float_t s) {fSigmaZ[i]=s;}
  Float_t GetNDeadZone() const {return fNDeadZone;}
  void SetNDeadZone(Float_t d) {fNDeadZone=d;}
  Int_t* ClIndex() {return fClIndex;}
  Int_t GetClIndex(Int_t i) const {return fClIndex[i];}
  void SetClIndex(Int_t i, Int_t c) {fClIndex[i]=c;}
  Float_t GetNormChi2(Int_t i) const {return fNormChi2[i];}
  void SetNormChi2(Int_t i, Float_t n) {fNormChi2[i]=n;}
  Bool_t GetConstrain() const {return fConstrain;}
  void SetConstrain(Bool_t c) {fConstrain=c;}
  Float_t GetExpQ() const {return fExpQ;}
  void SetExpQ(Float_t f) {fExpQ=f;}
  Float_t GetNormQ(Int_t i) const {return fNormQ[i];}
  void SetNormQ(Int_t i, Float_t q) {fNormQ[i]=q;}
  Float_t GetdEdxMismatch() const {return fdEdxMismatch;}
  void SetdEdxMismatch(Float_t m) {fdEdxMismatch=m;}
  Float_t GetNy(Int_t i) const {return fNy[i];}
  void SetNy(Int_t i, Float_t f) {fNy[i]=f;}
  Float_t GetNz(Int_t i) const {return fNz[i];}
  void SetNz(Int_t i, Float_t f) {fNz[i]=f;}
  Bool_t GetGoldV0() const {return fGoldV0;}
  void SetGoldV0(Bool_t g) {fGoldV0=g;}
  Float_t GetChi22() const {return fChi22;}
  void SetChi22(Float_t c) {fChi22=c;}
  Float_t GetDeadZoneProbability() const {return fDeadZoneProbability;}
  void SetDeadZoneProbability(Float_t d) {fDeadZoneProbability=d;}

  Double_t GetPredictedChi2MI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz) const;
  Bool_t IsGoldPrimary();
protected:

  Float_t fNUsed;                          // number of shared clusters
  Float_t fNSkipped;                       // number of skipped clusters
  Float_t fNDeadZone;                     // number of clusters in dead zone
  Float_t fDeadZoneProbability;          // probability to cross dead zone
  Bool_t  fReconstructed;                 // reconstructed - accepted flag
  Float_t fChi2MIP[12];                   // MIP chi squres 

  Float_t fDy[12];           //dy in layer
  Float_t fDz[12];           //dz in layer
  Float_t fSigmaY[12];       //sigma y 
  Float_t fSigmaZ[12];       //sigma z
  Float_t fNy[6];              //expected size of cluster
  Float_t fNz[6];              //expected size of cluster
  Float_t fD[2];            //distance to the vertex
  Float_t fDnorm[2];        // normalized distance to the vertex
  Float_t fNormQ[6];        // normalized Q
  Float_t fExpQ;            // expected Q
  Float_t fNormChi2[6];     // normalized chi2 
  Float_t fChi22;           // chi22
  Float_t fdEdxMismatch;    
  Bool_t fConstrain;        //indication of the vertex constrain
  Int_t  fClIndex[6];       //cluster Index
  Bool_t fGoldV0;           //corresponding gold V0 found
  ClassDef(AliITStrackMI,1)   //ITS reconstructed track
};

#endif


