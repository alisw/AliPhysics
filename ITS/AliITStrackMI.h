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
  friend class AliITStrackerV2;
  friend class AliITStrackerMI;
public:
  AliITStrackMI();
  AliITStrackMI(AliESDtrack& t,Bool_t c=kFALSE) throw (const Char_t *);
  AliITStrackMI(const AliITStrackMI& t);
  Int_t GetProlongationFast(Double_t alpha, Double_t xr,Double_t &y, Double_t &z);
  Int_t UpdateMI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz, Double_t chi2,UInt_t i);  
  Int_t CorrectForMaterial(Double_t d, Double_t x0=21.82);

  void UpdateESDtrack(ULong_t flags);

  void SetReconstructed(Bool_t sr=kTRUE){fReconstructed = sr;}  
  Bool_t GetReconstructed() const {return fReconstructed;}
  void SetChi2MIP(Int_t i,Float_t val){fChi2MIP[i]=val;}
  Float_t GetChi2MIP(Int_t i) const {return fChi2MIP[i];}  
  void IncrementNSkipped(){fNSkipped++;} // increment by 1 the # of skipped cls
  Float_t GetNSkipped() const {return fNSkipped;}
  void IncrementNUsed(){fNUsed++;} // increment by 1 the # of shared clusters
  Float_t GetNUsed() const {return fNUsed;}

  Int_t Compare(const TObject *o) const;
  Double_t GetCov33() const {return fC33;} // cov. matrix el. 3,3
  Double_t GetCov44() const {return fC44;} // cov. matrix el. 4,4
  Float_t GetDy(Int_t i) const {return fDy[i];}
  Float_t GetDz(Int_t i) const {return fDz[i];}
  Float_t GetSigmaY(Int_t i) const {return fSigmaY[i];}
  Float_t GetSigmaZ(Int_t i) const {return fSigmaZ[i];}

  Double_t GetPredictedChi2MI(Double_t cy, Double_t cz, Double_t cerry, Double_t cerrz) const;
 
protected:

  Float_t fNUsed;                          // number of shared clusters
  Float_t fNSkipped;                       // number of skipped clusters
  Float_t fNDeadZone;                     // number of clusters in dead zone
  Float_t fDeadZoneProbability;          // probability to cross dead zone
  Bool_t fReconstructed;                 // reconstructed - accepted flag
  Float_t fChi2MIP[12];                   // MIP chi squres 

  Float_t fDy[12];           //dy in layer
  Float_t fDz[12];           //dz in layer
  Float_t fSigmaY[12];       //sigma y 
  Float_t fSigmaZ[12];       //sigma z
  Float_t fNy[6];              //expected size of cluster
  Float_t fNz[6];              //expected size of cluster
  Float_t fD[2];            //distance to the vertex
  Float_t fNormQ[6];        // normalized Q
  Float_t fExpQ;            // expected Q
  Float_t fNormChi2[6];     // normalized chi2 
  Float_t fChi22;           // chi22
  Float_t fdEdxMismatch;    
  Bool_t fConstrain;        //indication of the vertex constrain
  Int_t  fClIndex[6];       //cluster Index

  ClassDef(AliITStrackMI,1)   //ITS reconstructed track
};

#endif


