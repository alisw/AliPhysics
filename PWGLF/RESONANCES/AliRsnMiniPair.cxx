#include "AliRsnMiniParticle.h"
#include "AliRsnMiniEvent.h"
#include "AliRsnMiniPair.h"

#include <AliQnCorrectionsQnVector.h>

ClassImp(AliRsnMiniPair)

//__________________________________________________________________________________________________
void AliRsnMiniPair::Fill
(AliRsnMiniParticle *p1, AliRsnMiniParticle *p2, Double_t m1, Double_t m2, Double_t refMass)
{
//
// Fill this object with data coming
// from arguments
//
   p1->Set4Vector(fP1[0], m1, kFALSE);
   p2->Set4Vector(fP2[0], m2, kFALSE);
   p1->Set4Vector(fP1[1], m1, kTRUE );
   p2->Set4Vector(fP2[1], m2, kTRUE );
   
   fDCA1 = p1->DCA();
   fDCA2 = p2->DCA();  

   fMother = -1;
   fIsFromB = kFALSE;
   fIsQuarkFound = kFALSE;
   fPmother[0] = -1.0;
   fPmother[1] = -1.0;
   fPmother[2] = -1.0;
   if (p1->Mother() == p2->Mother()) {
      fMother = p1->Mother();
      fMotherPDG = p1->MotherPDG();
      fPmother[0] = p1->PmotherX();
      fPmother[1] = p1->PmotherY();
      fPmother[2] = p1->PmotherZ();
      fIsFromB = p1->IsFromB();
      fIsQuarkFound = p1->IsQuarkFound();
   }

   Int_t i;
   for (i = 0; i < 2; i++) {
      fSum[i] = fP1[i] + fP2[i];
      fRef[i].SetXYZM(fSum[i].X(), fSum[i].Y(), fSum[i].Z(), refMass);
   }

   fNSisters=-1;
   if (p1->NTotSisters()==p2->NTotSisters()) fNSisters = p1->NTotSisters();

   fContainsV0Daughter = kFALSE;
   if (p1->IndexV0Pos() == p2->Index()) fContainsV0Daughter = kTRUE;
   if (p1->IndexV0Neg() == p2->Index()) fContainsV0Daughter = kTRUE;
   if (p2->IndexV0Pos() == p1->Index()) fContainsV0Daughter = kTRUE;
   if (p2->IndexV0Neg() == p1->Index()) fContainsV0Daughter = kTRUE;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::CosThetaStar(Bool_t useMC)
{
//
// Return cosine of angle of one daughter to the resonance momentum in its rest frame
//

   TLorentzVector &mother    = fSum[ID(useMC)];
   TLorentzVector &daughter0 = fP1[ID(useMC)];
//    TLorentzVector &daughter1 = fP2[ID(useMC)];
   TVector3 momentumM(mother.Vect());
   TVector3 normal(mother.Y() / momentumM.Mag(), -mother.X() / momentumM.Mag(), 0.0);

   // Computes first the invariant mass of the mother
//    Double_t mass0      = daughter0.M();
//    Double_t mass1      = daughter1.M();
//    Double_t p0         = daughter0.Vect().Mag();
//    Double_t p1         = daughter1.Vect().Mag();
//    Double_t E0         = TMath::Sqrt(mass0 * mass0 + p0 * p0);
//    Double_t E1         = TMath::Sqrt(mass1 * mass1 + p1 * p1);
//    Double_t MotherMass = TMath::Sqrt((E0 + E1) * (E0 + E1) - (p0 * p0 + 2.0 * daughter0.Vect().Dot(daughter1.Vect()) + p1 * p1));
//    MotherMass = mother.M();

   // Computes components of beta
   Double_t betaX = -mother.X() / mother.E();
   Double_t betaY = -mother.Y() / mother.E();
   Double_t betaZ = -mother.Z() / mother.E();

   // Computes Lorentz transformation of the momentum of the first daughter
   // into the rest frame of the mother and theta*
   daughter0.Boost(betaX, betaY, betaZ);
   TVector3 momentumD = daughter0.Vect();

   Double_t cosThetaStar = normal.Dot(momentumD) / momentumD.Mag();

   return cosThetaStar;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::CosThetaJackson(Bool_t useMC)
{
//
// Return cosine of angle of one daughter in the resonance rest frame to beam z-axis (Jackson frame)
//

   TLorentzVector &mother    = fSum[ID(useMC)];
   TLorentzVector daughter = fP1[ID(useMC)];
//    TLorentzVector daughter = fP2[ID(useMC)];
   daughter.Boost(-mother.BoostVector());

//   TVector3 beamAxis(0,0,1);
//   TVector3 momentumD = daughter.Vect();
     // cos(theta) via dot product
//   Double_t cosTheta = momentumD.Dot(beamAxis)/TMath::Sqrt(momentumD.Mag2()*beamAxis.Mag2());
     // Faster way (less computing)
//   Double_t cosTheta = daughter.CosTheta();

   return daughter.CosTheta();
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::CosThetaTransversity(Bool_t useMC)
{
//
// Return cosine of angle of one daughter in the resonance rest frame to normal of
// beam z-axis and resonance production plane (Transversity frame)
//

	TLorentzVector &mother = fSum[ID(useMC)];
	TLorentzVector daughter = fP1[ID(useMC)];
	//    TLorentzVector daughter = fP2[ID(useMC)];
	daughter.Boost(-mother.BoostVector());

	TVector3 beamAxis(0, 0, 1);
	TVector3 transFrame = beamAxis.Cross(mother.Vect());
	TVector3 momentumD = daughter.Vect();

	return momentumD.Dot(transFrame) / TMath::Sqrt((momentumD.Mag2() * transFrame.Mag2()));
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::CosThetaToEventPlane(AliRsnMiniEvent *event, Bool_t useMC)
{
//
// Return cosine of angle of one daughter in the resonance rest frame to Quantization axis.
// Quantization axis - perpendicular between event plane and beam direction
//

  // Get QnVector
  AliQnCorrectionsQnVector *qnVect = event->GetQnVector();
  if (!qnVect) return 0;

  TLorentzVector &mother = fSum[ID(useMC)];
  TLorentzVector daughter = fP1[ID(useMC)];
  //    TLorentzVector daughter = fP2[ID(useMC)];
  daughter.Boost(-mother.BoostVector());

  TVector3 evPlaneVect(qnVect->Qx(1), qnVect->Qy(1), 0);
  TVector3 beamAxis(0, 0, 1);
  TVector3 quantizationAxis = beamAxis.Cross(evPlaneVect);
  TVector3 momentumD = daughter.Vect();

  Double_t cosTheta = momentumD.Dot(quantizationAxis)/TMath::Sqrt(momentumD.Mag2()*quantizationAxis.Mag2());
  return cosTheta;
}

//__________________________________________________________________________________________________
void AliRsnMiniPair::InvertP(Bool_t first)
{
//
// Inverts one 4-momentum and recompute sum
//

   Int_t i;
   for (i = 0; i < 2; i++) {
      if (first) fP1[i].SetVect(fP1[i].Vect() *= -1.0);
      else       fP2[i].SetVect(fP2[i].Vect() *= -1.0);
      fSum[i] = fP1[i] + fP2[i];
      fRef[i].SetXYZM(fSum[i].X(), fSum[i].Y(), fSum[i].Z(), fRef[i].M());
   }
}

//__________________________________________________________________________________________________
void AliRsnMiniPair::FillRef(Double_t mass)
{
//
// Fill ref 4-vectors using the passed mass and the values in 'sum'
//

   Int_t i;
   for (i = 0; i < 2; i++) {
      fRef[i].SetXYZM(fSum[i].X(), fSum[i].Y(), fSum[i].Z(), mass);
   }
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::InvMassRes() const
{
//
// Return invariant mass resolution
//

   if (fSum[1].M() <= 0.0) return 1E20;

   return (fSum[0].M() - fSum[1].M()) / fSum[1].M();
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::InvMassDiff() const
{
//
// Return invariant mass resolution
//

   if (fSum[1].M() <= 0.0) return 1E20;

   return (fSum[0].M() - fSum[1].M());
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::PtRatio(Bool_t mc) const
{
//
// Return ratio of transverse momenta of daughters
//

   Double_t num = TMath::Abs(fP1[ID(mc)].Perp() - fP2[ID(mc)].Perp());
   Double_t den = TMath::Abs(fP1[ID(mc)].Perp() + fP2[ID(mc)].Perp());

   if (den <= 0.0) return 1E20;

   return num / den;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::DipAngle(Bool_t mc) const
{
//
// Opening angle in a Z-T space
//

   const TLorentzVector &p1 = fP1[ID(mc)];
   const TLorentzVector &p2 = fP2[ID(mc)];

   return ((p1.Perp() * p2.Perp() + p1.Z() * p2.Z()) / p1.Mag() / p2.Mag());
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::DeltaCos(Bool_t mc) const
{
//
// Difference of cos(theta) angles
// - alpha : angle between particles of a pair in the case
// when they are daughters of the resonance with the mass M
// - beta : angle between particles of a pair
// More info in Phys.Rev.C65 (2002) 034909

   const TLorentzVector &p1 = fP1[ID(mc)];
   const TLorentzVector &p2 = fP2[ID(mc)];
   const TLorentzVector &mother = fRef[ID(mc)];

   TVector3 p1Vect = p1.Vect();
   TVector3 p2Vect = p2.Vect();

   Double_t magP1P2 = TMath::Sqrt(p1Vect.Mag2()*p2Vect.Mag2());

   Double_t cosA = (p1.E()*p2.E() - 0.5*(mother.M()*mother.M() - p1.M()*p1.M() - p2.M()*p2.M()))/magP1P2;
   Double_t cosB = p1Vect.Dot(p2Vect)/magP1P2;

   return cosB-cosA;
}




//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::PhiV(Bool_t mc)
{
//
// Opening angle between two tracks/decay particles.

  TLorentzVector &p1 = fP1[ID(mc)];
  TLorentzVector &p2 = fP2[ID(mc)];
  
  TVector3 P1 = p1.Vect();
  TVector3 P2 = p2.Vect();
  
  //   Double_t magP1P2 = TMath::Sqrt(p1Vect.Mag2()*p2Vect.Mag2());
  //   Double_t cosB = p1Vect.Dot(p2Vect)/magP1P2;
  //   return cosB;
  
  
  Double_t r(0);
  
  //Randomization of pair ordering (Symmetric Peaks)
  do { r = gRandom -> Uniform (0.0,1.0); } while (r==0.5);
  
  //  if (r < 0.5) { P1.SetXYZ (track1->Px(),track1->Py(),track1->Pz()); P2.SetXYZ (track2->Px(),track2->Py(),track2->Pz()); }
  //    if (r > 0.5) { P1.SetXYZ (track2->Px(),track2->Py(),track2->Pz()); P2.SetXYZ (track1->Px(),track1->Py(),track1->Pz()); }
  
  
  if (r < 0.5) { 
    P1.SetXYZ (fP1[ID(mc)].Px(), fP1[ID(mc)].Py(), fP1[ID(mc)].Pz()); 
     P2.SetXYZ (fP2[ID(mc)].Px(), fP2[ID(mc)].Py(), fP2[ID(mc)].Pz());
  }
  
  
  if (r > 0.5) { 
    P1.SetXYZ (fP2[ID(mc)].Px(), fP2[ID(mc)].Py(), fP2[ID(mc)].Pz());
    P2.SetXYZ (fP1[ID(mc)].Px(), fP1[ID(mc)].Py(), fP1[ID(mc)].Pz()); 
  }
  
  
  
   //PhiV Calculation  === as per note
  
  TVector3 P = P1 + P2;
  TVector3 U ( P.X()/P.Mag(),P.Y()/P.Mag(),P.Z()/P.Mag() );
  TVector3 A ( U.Y()/TMath::Sqrt(U.X()*U.X()+U.Y()*U.Y()),-U.X()/TMath::Sqrt(U.X()*U.X()+ U.Y()*U.Y()),0 );
  TVector3 Vp = P1.Cross(P2);
   TVector3 V (Vp.X()/Vp.Mag(),Vp.Y()/Vp.Mag(),Vp.Z()/Vp.Mag());
   TVector3 W = U.Cross(V);
   
   return (180.0/TMath::Pi())*A.Angle(W);
}




//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::DaughterPt(Int_t daughterId, Bool_t mc)
{
  //returns pt of the <id> daughter 
  // if MC returns generated momenta
  if (daughterId==0)
    return fP1[ID(mc)].Pt();
  else 
    return fP2[ID(mc)].Pt();
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::DaughterDCA(Int_t daughterId)
{   
  //
  //returns dca to Primary Vertex of the <id> daughter 
  //

  if (daughterId==0)
    return fDCA1;
  else 
    return fDCA2;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::DCAProduct()
{
  //
  //returns products of the DCA of the 2 daughters
  //
  
    return fDCA1*fDCA2;
}

//__________________________________________________________________________________________________
void AliRsnMiniPair::DaughterPxPyPz(Int_t daughterId, Bool_t mc, Double_t *pxpypz)
{
  //returns px,py,pz of the <id> daughter by saving them into pxpypz
  // if MC returns generated momenta
  if (!pxpypz)    return;
 
  if (daughterId==0){
    pxpypz[0]=fP1[ID(mc)].Px();
    pxpypz[1]=fP1[ID(mc)].Py();
    pxpypz[2]=fP1[ID(mc)].Pz();
  } else {
    pxpypz[0]=fP2[ID(mc)].Px();
    pxpypz[1]=fP2[ID(mc)].Py();
    pxpypz[2]=fP2[ID(mc)].Pz();
  }
  return;
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::PairPtRes() const
{
//
// Return pair pt resolution
//
   if (Pt(1) <= 0.0) return 1E20;
   return (Pt(0) - Pt(1)) / Pt(1);
}

//__________________________________________________________________________________________________
Double_t AliRsnMiniPair::PairYRes() const
{
//
// Return pair rapidity resolution
//
  if (Y(1) <= 0.0) return 1E20;
  return (Y(0) - Y(1)) / Y(1);
}
