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
Double_t AliDielectronPair::ThetaPhiCM(const AliVParticle* d1, const AliVParticle* d2, 
				       const Bool_t isHE, const Bool_t isTheta) {
  // The function calculates theta and phi in the mother rest frame with 
  // respect to the helicity coordinate system and Collins-Soper coordinate system
  // TO DO: generalize for different decays (only J/Psi->e+e- now)

  // Laboratory frame 4-vectors:
  // projectile beam & target beam 4-mom
  const Double_t kBeamEnergy   = 3500.;      //TODO: need to retrieve the beam energy from somewhere
  TLorentzVector projMom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))); 
  TLorentzVector targMom(0.,0.,kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))); 
  
  // first & second daughter 4-mom
  TLorentzVector p1Mom(d1->Px(),d1->Py(),d1->Pz(),TMath::Sqrt(d1->Px()*d1->Px()+d1->Py()*d1->Py()+d1->Pz()*d1->Pz()+AliPID::ParticleMass(AliPID::kElectron)*AliPID::ParticleMass(AliPID::kElectron)));
  TLorentzVector p2Mom(d2->Px(),d2->Py(),d2->Pz(),TMath::Sqrt(d2->Px()*d2->Px()+d2->Py()*d2->Py()+d2->Pz()*d2->Pz()+AliPID::ParticleMass(AliPID::kElectron)*AliPID::ParticleMass(AliPID::kElectron)));
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
  const Double_t kBeamEnergy   = 3500.;      //TODO: need to retrieve the beam energy from somewhere
  TLorentzVector projMom(0.,0.,-kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))); 
  TLorentzVector targMom(0.,0.,kBeamEnergy,TMath::Sqrt(kBeamEnergy*kBeamEnergy+AliPID::ParticleMass(AliPID::kProton)*AliPID::ParticleMass(AliPID::kProton))); 

  // first & second daughter 4-mom
  AliVParticle *d1 = dynamic_cast<AliVParticle*>(fRefD1.GetObject());
  AliVParticle *d2 = dynamic_cast<AliVParticle*>(fRefD2.GetObject());
  TLorentzVector p1Mom(d1->Px(),d1->Py(),d1->Pz(),TMath::Sqrt(d1->Px()*d1->Px()+d1->Py()*d1->Py()+d1->Pz()*d1->Pz()+AliPID::ParticleMass(AliPID::kElectron)*AliPID::ParticleMass(AliPID::kElectron)));
  TLorentzVector p2Mom(d2->Px(),d2->Py(),d2->Pz(),TMath::Sqrt(d2->Px()*d2->Px()+d2->Py()*d2->Py()+d2->Pz()*d2->Pz()+AliPID::ParticleMass(AliPID::kElectron)*AliPID::ParticleMass(AliPID::kElectron)));
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
