#ifndef AliBarrelTrack_H
#define AliBarrelTrack_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// General class for barrel tracks
// This class contains all the information
// to describe the tracks detected by the barre detectors (ITS, TPC, TRD, TOF)
//

#include "TObject.h"
#include "TMath.h"

class AliBarrelTrack : public TObject {

public:
  
  AliBarrelTrack();
  ~AliBarrelTrack() {}
  
  // Setters

  void SetLabel(Int_t label); 
  void SetX(Double_t x, Double_t alpha);
  void SetRefPlane(Int_t nRefPlane, Int_t isIn);
  void SetNClusters(Int_t nClusters, Double_t chi2);
  void SetTime(Double_t time[5], Double_t length);
  void SetStateVector(Double_t vec[5]);             // external parameters
  void SetCovarianceMatrix(Double_t vec[15]);       // external parameters
 
  void SetNWrongClusters(Int_t n) {fNWrong = n;}
  void SetNRotate(Int_t n) {fNRotate = n;}

  void SetMass(Double_t mass) {fMass = mass;}
  void SetdEdX(Double_t dEdX) {fdEdX = dEdX;}

  // Getters

  // standard
  Int_t GetLabel() const {return fLabel;}
  Int_t GetNClusters() const {return fNClusters;}
  Int_t GetNWrongClusters() const {return fNWrong;}
  Int_t GetNRotate() const {return fNRotate;}
  
  Int_t GetRefPlane() const {return fRefPlane;}

  Double_t GetX() const {return fX;}
  Double_t Y() const {return fY;}
  Double_t Z() const {return fZ;}

  Double_t GetMass() const {return fMass;}

  // track oriented variables
  Double_t Pt() const {return (f1Pt==0)? 0 : 1./f1Pt;}
  Double_t TgLambda() const {return fTgLambda;}
  Double_t Eta() const {return -TMath::Log(TMath::Tan(0.25 *TMath::Pi()-0.5*Lambda()));}
  Double_t Phi() const {return TMath::ASin(fSnPhi);}
  

  // uncertainties
  Double_t DeltaY() const {return TMath::Sqrt(fCy);}
  Double_t DeltaZ() const {return TMath::Sqrt(fCz);}
  Double_t DeltaTgLambda() const {return TMath::Sqrt(fCtg);}
  Double_t DeltaSnPhi() const {return TMath::Sqrt(fCphi);}
  Double_t DeltaPt() const {return Pt() - 1./(Pt() + TMath::Sqrt(fCpt));}


  // reference oriented variables
  Double_t Px() const {return TMath::Cos(Phi()+fAlpha) * Pt();}
  Double_t Py() const {return TMath::Sin(Phi()+fAlpha) * Pt();}
  Double_t Pz() const {return fTgLambda*Pt();}
  Double_t P()  const {return Pt()*(fTgLambda+1);}

  Double_t Lambda() const {return TMath::ATan(fTgLambda);}
  
protected:

  Int_t fLabel;          // kine tree index

  Int_t fRefPlane;       // id of the reference plane
  Int_t fIsIn;           // direction

  Double_t fX;           // Kalman Time
  Double_t fAlpha;       // sector angle

  // state vector
  Double_t fZ;          // Z in global cs
  Double_t fY;          // Y in local cs corresponds to r-phi
  Double_t fTgLambda;   // Tangent of the dip angle
  Double_t fSnPhi;      // Sin 
  Double_t f1Pt;        // inverse of momentum


  // covariance matrix
  Double_t fCz;   // z element of covariance matrix
  Double_t fCy;   // y element of covariance matrix
  Double_t fCtg;  // tangent element of covariance matrix
  Double_t fCphi; // phi element of covariance matrix
  Double_t fCpt;  // pt element of covariance matrix

    
  // track time/length  
  Double_t fTimeHypothesis[5];    // time for all hypoptheses
  Double_t fLength;               // track length

  // performance info
  Int_t fNClusters;         // Number of clusters 
  Int_t fNWrong;            // Number of wrong clusters
  Double_t fChi2;           // Chi 2
  Int_t fNRotate;           // number of rotations / sector crossing
 
  Double_t fMass;           // mass hypothesis
  Double_t fdEdX;           // dE/dX

  ClassDef(AliBarrelTrack,1)
};


#endif
