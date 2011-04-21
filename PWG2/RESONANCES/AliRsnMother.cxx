/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

////////////////////////////////////////////////////////////////////////////////
//
//  This class implements a candidate resonance. It has two pointers to its
//  two candidate daughters, whose 4-momenta are combined to obtain the mother
//  invariant mass and other kinematical quantities.
//  This class contains also some methods used to compute kinematical relations
//  between the candidate resonance and other particles.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TVector3.h>

#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"

#include "AliRsnMother.h"

ClassImp(AliRsnMother)

//__________________________________________________________________________________________________
AliRsnMother::AliRsnMother(const AliRsnMother &obj) :
   TObject(obj),
   fRefEvent(obj.fRefEvent),
   fSum(obj.fSum),
   fRef(obj.fRef)
{
//
// Copy constructor.
// Does not duplicate pointers.
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = obj.fDaughter[i];
}

//__________________________________________________________________________________________________
AliRsnMother& AliRsnMother::operator=(const AliRsnMother &obj)
{
//
// Assignment operator.
// Does not duplicate pointers.
//

   fSum = obj.fSum;
   fRef = obj.fRef;
   fRefEvent = obj.fRefEvent;
   fDaughter[0] = obj.fDaughter[0];
   fDaughter[1] = obj.fDaughter[1];
   
   return (*this);
}

//__________________________________________________________________________________________________
AliRsnMother::~AliRsnMother()
{
//
// Desctructor.
// Does nothing, since pointers are not created in this class.
//
}

//_______________________________________________________________________________________________________________________
void AliRsnMother::Reset()
{
//
// Resets the mother, zeroing all data members.
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = 0x0;
   fRefEvent = 0x0;
   fSum.SetXYZM(0.0, 0.0, 0.0, 0.0);
   fRef.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//__________________________________________________________________________________________________
Int_t AliRsnMother::CommonMother() const
{
//
// If MC info is available, checks if the two tracks in the pair have the same mother.
// If the mother label is the same, the function returns the PDG code of mother,
// otherwise it returns 0.
// The two arguments passed by reference contain the GEANT labels of the mother
// of the two particles to which the two daughters point. This is for being able 
// to check if they are really coming from a resonance (indexes >= 0) or not.
//

   Int_t m1  = fDaughter[0]->GetMother();
   Int_t m2  = fDaughter[1]->GetMother();
   Int_t out = 0;
   
   // a true mother makes sense only if both mothers
   // are not-negative and equal
   if (m1 >= 0 && m2 >= 0 && m1 == m2) {
      out = TMath::Abs(fDaughter[0]->GetMotherPDG());
      AliDebugClass(1, Form("Mothers are: %d %d --> EQUAL (PDG = %d)", m1, m2, out));
   } 
   
   return out;
}

//__________________________________________________________________________________________________
Double_t AliRsnMother::AngleToLeading(Bool_t &success)
{
//
// Compute the angle betwee this and the leading particls
// of the reference event (if this was set properly).
// In case one of the two daughters is the leading, return
// a meaningless value, in order to skip this pair.
// if second argument is kTRUE, use MC values.
//

   if (!fRefEvent) {
      success = kFALSE;
      return -99.0;
   }
   
   Int_t id1 = fDaughter[0]->GetID();
   Int_t id2 = fDaughter[1]->GetID();
   Int_t idL = fRefEvent->GetLeadingIndex();
   
   if (id1 == idL || id2 == idL) {
      success = kFALSE;
      return -99.0;
   }
   
   AliRsnDaughter leading = fRefEvent->GetDaughter(idL, kFALSE);
   AliVParticle  *ref     = leading.GetRef();
   Double_t       angle   = ref->Phi() - fSum.Phi();
   
   //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
   while (angle >= 1.5 * TMath::Pi()) angle -= 2 * TMath::Pi();
   while (angle < -0.5 * TMath::Pi()) angle += 2 * TMath::Pi();
   success = kTRUE;
   
   return angle;
}

//__________________________________________________________________________________________________
void AliRsnMother::ComputeSum(Double_t mass0, Double_t mass1, Double_t motherMass)
{
//
// Sets the masses for the 4-momenta of the daughters and then
// sums them, taking into account that the space part is set to
// each of them when the reference object is set (see AliRsnDaughter::SetRef)
//

   fDaughter[0]->FillP(mass0);
   fDaughter[1]->FillP(mass1);

   // sum
   fSum   = fDaughter[0]->Prec() + fDaughter[1]->Prec();
   fSumMC = fDaughter[0]->Psim() + fDaughter[1]->Psim();
   
   // reference
   fRef.SetXYZM(fSum.X(), fSum.Y(), fSum.Z(), motherMass);
   fRefMC.SetXYZM(fSumMC.X(), fSumMC.Y(), fSumMC.Z(), motherMass);
}

//__________________________________________________________________________________________________
Double_t AliRsnMother::CosThetaStar(Bool_t first, Bool_t useMC)
{
//
// Computes the cosine of theta*, which is the angle of one of the daughters
// with respect to the total momentum of the resonance, in its rest frame.
// The arguments are needed to choose which of the daughters one want to use
// and if reconstructed or MC momentum must be used.
// [Contribution from Z. Feckova]
//

   TLorentzVector &mother    = Sum(useMC);
   TLorentzVector &daughter0 = (first ? fDaughter[0]->P(useMC) : fDaughter[1]->P(useMC));
   TLorentzVector &daughter1 = (first ? fDaughter[1]->P(useMC) : fDaughter[0]->P(useMC));
   TVector3 momentumM(mother.Vect());
   TVector3 normal(mother.Y() / momentumM.Mag(), -mother.X() / momentumM.Mag(), 0.0);

   // Computes first the invariant mass of the mother
   Double_t mass0      = fDaughter[0]->P(useMC).M();
   Double_t mass1      = fDaughter[1]->P(useMC).M();
   Double_t p0         = daughter0.Vect().Mag();
   Double_t p1         = daughter1.Vect().Mag();
   Double_t E0         = TMath::Sqrt(mass0 * mass0 + p0 * p0);
   Double_t E1         = TMath::Sqrt(mass1 * mass1 + p1 * p1);
   Double_t MotherMass = TMath::Sqrt((E0 + E1) * (E0 + E1) - (p0 * p0 + 2.0 * daughter0.Vect().Dot(daughter1.Vect()) + p1 * p1));
   MotherMass = fSum.M();

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
void AliRsnMother::PrintInfo(const Option_t * /*option*/) const
{
//
// Print some info of the pair.
// The options are passed to the AliRsnDaughter::Print() method
//

   AliInfo("======== BEGIN PAIR INFO ===========");
   AliInfo("Track #1");
   fDaughter[0]->Print();
   AliInfo("Track #2");
   fDaughter[1]->Print();
   AliInfo("========= END PAIR INFO ===========");
}

//__________________________________________________________________________________________________
Bool_t AliRsnMother::CheckPair(Bool_t checkMC) const
{
//
// Checks that the pair is well initialized:
// - both daughters are good pointers
// - if MC is required, both daughters have a MC reference
//

   if (!fDaughter[0] || !fDaughter[1]) {
      AliError("One of the two tracks is NULL in this pair!");
      return kFALSE;
   }

   if (checkMC) {
      if (fDaughter[0]->GetRefMC() == 0x0 || fDaughter[1]->GetRefMC() == 0x0) {
         AliError("Required MC info but not all MC refs are available");
         return kFALSE;
      }
   }

   return kTRUE;
}
