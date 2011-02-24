/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

/* $Id$ */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Dielectron Pair class. Internally it makes use of AliKFParticle.     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


#include "AliDielectronPair.h"
#include "AliVTrack.h"
#include "AliPID.h"

ClassImp(AliDielectronPair)

AliDielectronPair::AliDielectronPair() :
  fType(-1),
  fLabel(-1),
  fPair(),
  fD1(),
  fD2(),
  fRefD1(),
  fRefD2()
{
  //
  // Default Constructor
  //
  
}

//______________________________________________
AliDielectronPair::AliDielectronPair(AliVTrack * const particle1, Int_t pid1,
                                     AliVTrack * const particle2, Int_t pid2, Char_t type) :
  fType(type),
  fLabel(-1),
  fPair(),
  fD1(),
  fD2(),
  fRefD1(),
  fRefD2()
{
  //
  // Constructor with tracks
  //
  SetTracks(particle1, pid1, particle2, pid2);
}

//______________________________________________
AliDielectronPair::~AliDielectronPair()
{
  //
  // Default Destructor
  //
  
}

//______________________________________________
void AliDielectronPair::SetTracks(AliVTrack * const particle1, Int_t pid1,
                                  AliVTrack * const particle2, Int_t pid2)
{
  //
  // Sort particles by pt, first particle larget Pt
  // set AliKF daughters and pair
  //
  fPair.Initialize();
  fD1.Initialize();
  fD2.Initialize();

  AliKFParticle kf1(*particle1,pid1);
  AliKFParticle kf2(*particle2,pid2);

  fPair.AddDaughter(kf1);
  fPair.AddDaughter(kf2);

  if (particle1->Pt()>particle2->Pt()){
    fRefD1 = particle1;
    fRefD2 = particle2;
    fD1+=kf1;
    fD2+=kf2;
  } else {
    fRefD1 = particle2;
    fRefD2 = particle1;
    fD1+=kf2;
    fD2+=kf1;
  }
}

//______________________________________________
void AliDielectronPair::GetThetaPhiCM(Double_t &thetaHE, Double_t &phiHE, Double_t &thetaCS, Double_t &phiCS) const
{
  //
  // Calculate theta and phi in helicity and Collins-Soper coordinate frame
  //
  const Double_t kBeamEnergy   = 3500.;
  Double_t pxyz1[3]={0,0,0};
  Double_t pxyz2[3]={0,0,0};
  Double_t eleMass=AliPID::ParticleMass(AliPID::kElectron);
  Double_t proMass=AliPID::ParticleMass(AliPID::kProton);
  
  AliVParticle *d1 = static_cast<AliVParticle*>(fRefD1.GetObject());
  AliVParticle *d2 = static_cast<AliVParticle*>(fRefD2.GetObject());

  d1->PxPyPz(pxyz1);
  d2->PxPyPz(pxyz2);
  
  TLorentzVector projMom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  TLorentzVector targMom(0.,0., kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  
  // first & second daughter 4-mom
  TLorentzVector p1Mom(pxyz1[0],pxyz1[1],pxyz1[2],
                       TMath::Sqrt(pxyz1[0]*pxyz1[0]+pxyz1[1]*pxyz1[1]+pxyz1[2]*pxyz1[2]+eleMass*eleMass));
  TLorentzVector p2Mom(pxyz2[0],pxyz2[1],pxyz2[2],
                       TMath::Sqrt(pxyz2[0]*pxyz2[0]+pxyz2[1]*pxyz2[1]+pxyz2[2]*pxyz2[2]+eleMass*eleMass));
  // J/Psi 4-momentum vector
  TLorentzVector motherMom=p1Mom+p2Mom;
  
  // boost all the 4-mom vectors to the mother rest frame
  TVector3 beta = (-1.0/motherMom.E())*motherMom.Vect();
  p1Mom.Boost(beta);
  p2Mom.Boost(beta);
  projMom.Boost(beta);
  targMom.Boost(beta);

    // x,y,z axes
  TVector3 zAxisHE = (motherMom.Vect()).Unit();
  TVector3 zAxisCS = ((projMom.Vect()).Unit()-(targMom.Vect()).Unit()).Unit();
  TVector3 yAxis = ((projMom.Vect()).Cross(targMom.Vect())).Unit();
  TVector3 xAxisHE = (yAxis.Cross(zAxisHE)).Unit();
  TVector3 xAxisCS = (yAxis.Cross(zAxisCS)).Unit();
  
  // fill theta and phi
  if(d1->Charge()>0){
    thetaHE = zAxisHE.Dot((p1Mom.Vect()).Unit());
    thetaCS = zAxisCS.Dot((p1Mom.Vect()).Unit());
    phiHE   = TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxisHE));
    phiCS   = TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxisCS));
  } else {
    thetaHE = zAxisHE.Dot((p2Mom.Vect()).Unit());
    thetaCS = zAxisCS.Dot((p2Mom.Vect()).Unit());
    phiHE   = TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxisHE));
    phiCS   = TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxisCS));
  }
}

//______________________________________________
Double_t AliDielectronPair::ThetaPhiCM(const AliVParticle* d1, const AliVParticle* d2, 
                                       const Bool_t isHE, const Bool_t isTheta)
{
  // The function calculates theta and phi in the mother rest frame with
  // respect to the helicity coordinate system and Collins-Soper coordinate system
  // TO DO: generalize for different decays (only J/Psi->e+e- now)

  // Laboratory frame 4-vectors:
  // projectile beam & target beam 4-mom
  // TODO: need to retrieve the beam energy from somewhere
  const Double_t kBeamEnergy   = 3500.;
  Double_t px1=d1->Px();
  Double_t py1=d1->Py();
  Double_t pz1=d1->Pz();
  Double_t px2=d2->Px();
  Double_t py2=d2->Py();
  Double_t pz2=d2->Pz();
  Double_t eleMass=AliPID::ParticleMass(AliPID::kElectron);
  Double_t proMass=AliPID::ParticleMass(AliPID::kProton);
  
  TLorentzVector projMom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  TLorentzVector targMom(0.,0., kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  
  // first & second daughter 4-mom
  TLorentzVector p1Mom(px1,py1,pz1,TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+eleMass*eleMass));
  TLorentzVector p2Mom(px2,py2,pz2,TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+eleMass*eleMass));
  // J/Psi 4-momentum vector
  TLorentzVector motherMom=p1Mom+p2Mom;
  
  // boost all the 4-mom vectors to the mother rest frame
  TVector3 beta = (-1.0/motherMom.E())*motherMom.Vect();
  p1Mom.Boost(beta);
  p2Mom.Boost(beta);
  projMom.Boost(beta);
  targMom.Boost(beta);
  
  // x,y,z axes 
  TVector3 zAxis;
  if(isHE) zAxis = (motherMom.Vect()).Unit();
  else zAxis = ((projMom.Vect()).Unit()-(targMom.Vect()).Unit()).Unit();
  TVector3 yAxis = ((projMom.Vect()).Cross(targMom.Vect())).Unit();
  TVector3 xAxis = (yAxis.Cross(zAxis)).Unit();
  
  // return either theta or phi
  if(isTheta) {
    if(d1->Charge()>0)
      return zAxis.Dot((p1Mom.Vect()).Unit());
    else 
      return zAxis.Dot((p2Mom.Vect()).Unit());

  }
  else {
    if(d1->Charge()>0)
      return TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxis));
    else
      return TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxis));
  }
}

//______________________________________________
Double_t AliDielectronPair::ThetaPhiCM(const Bool_t isHE, const Bool_t isTheta) const {
  // The function calculates theta and phi in the mother rest frame with 
  // respect to the helicity coordinate system and Collins-Soper coordinate system
  // TO DO: generalize for different decays (only J/Psi->e+e- now)

  // Laboratory frame 4-vectors:
  // projectile beam & target beam 4-mom
  AliVParticle *d1 = static_cast<AliVParticle*>(fRefD1.GetObject());
  AliVParticle *d2 = static_cast<AliVParticle*>(fRefD2.GetObject());
  
  const Double_t kBeamEnergy   = 3500.;
  Double_t px1=d1->Px();
  Double_t py1=d1->Py();
  Double_t pz1=d1->Pz();
  Double_t px2=d2->Px();
  Double_t py2=d2->Py();
  Double_t pz2=d2->Pz();
  Double_t eleMass=AliPID::ParticleMass(AliPID::kElectron);
  Double_t proMass=AliPID::ParticleMass(AliPID::kProton);
  
  TLorentzVector projMom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  TLorentzVector targMom(0.,0., kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+proMass*proMass));
  
  // first & second daughter 4-mom
  // first & second daughter 4-mom
  TLorentzVector p1Mom(px1,py1,pz1,TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+eleMass*eleMass));
  TLorentzVector p2Mom(px2,py2,pz2,TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+eleMass*eleMass));
  // J/Psi 4-momentum vector
  TLorentzVector motherMom=p1Mom+p2Mom;

  // boost all the 4-mom vectors to the mother rest frame
  TVector3 beta = (-1.0/motherMom.E())*motherMom.Vect();
  p1Mom.Boost(beta);
  p2Mom.Boost(beta);
  projMom.Boost(beta);
  targMom.Boost(beta);

  // x,y,z axes 
  TVector3 zAxis;
  if(isHE) zAxis = (motherMom.Vect()).Unit();
  else zAxis = ((projMom.Vect()).Unit()-(targMom.Vect()).Unit()).Unit();
  TVector3 yAxis = ((projMom.Vect()).Cross(targMom.Vect())).Unit();
  TVector3 xAxis = (yAxis.Cross(zAxis)).Unit();

  // return either theta or phi
  if(isTheta) {
    if(fD1.GetQ()>0) 
      return zAxis.Dot((p1Mom.Vect()).Unit());
    else
      return zAxis.Dot((p2Mom.Vect()).Unit());
  }
  else {
    if(fD1.GetQ()>0)
      return TMath::ATan2((p1Mom.Vect()).Dot(yAxis), (p1Mom.Vect()).Dot(xAxis));
    else
      return TMath::ATan2((p2Mom.Vect()).Dot(yAxis), (p2Mom.Vect()).Dot(xAxis));
  }
}
