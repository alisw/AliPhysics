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

//_____________________________________________________________________________
AliRsnMother::AliRsnMother(const AliRsnMother &obj) :
   TObject(obj),
   fRefEvent(obj.fRefEvent),
   fSum(obj.fSum),
   fSumMC(obj.fSumMC)
{
//
// Copy constructor.
// Does not duplicate pointers.
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = obj.fDaughter[i];
}

//_____________________________________________________________________________
AliRsnMother& AliRsnMother::operator=(const AliRsnMother &obj)
{
//
// Assignment operator.
// Does not duplicate pointers.
//

   Int_t i;

   fSum = obj.fSum;
   fSumMC = obj.fSumMC;

   for (i = 0; i < 2; i++) fDaughter[i] = obj.fDaughter[i];
   
   fRefEvent = obj.fRefEvent;

   return (*this);
}

//_____________________________________________________________________________
AliRsnMother::~AliRsnMother()
{
//
// Desctructor.
// Does nothing, since pointers are not created in this class.
//
}

//_____________________________________________________________________________
Int_t AliRsnMother::CommonMother(Int_t &m0, Int_t &m1) const
{
//
// If MC info is available, checks if the two tracks in the pair have the same mother.
// If the mother label is the same, the function returns the PDG code of mother,
// otherwise it returns 0.
// The two arguments passed by reference contain the GEANT labels of the mother
// of the two particles to which the two daughters point. This is for being able 
// to check if they are really coming from a resonance (indexes >= 0) or not.
//

   // if MC info is not available, the check can't be done
   if (!fDaughter[0]->GetRefMC() || !fDaughter[1]->GetRefMC()) {
      AliDebug(AliLog::kDebug, "MC Info absent --> cannot check common mother");
      return 0;
   }

   // check that labels are the same
   m0 = -1;
   m1 = -2;
   if (fDaughter[0]->IsESD() && fDaughter[1]->IsESD()) {
      if (fDaughter[0]->GetRefMCESD() && fDaughter[1]->GetRefMCESD()) {
         m0 = fDaughter[0]->GetRefMCESD()->Particle()->GetFirstMother();
         m1 = fDaughter[1]->GetRefMCESD()->Particle()->GetFirstMother();
      }
   }
   if (fDaughter[0]->IsAOD() && fDaughter[1]->IsAOD()) {
      if (fDaughter[0]->GetRefMCAOD() && fDaughter[1]->GetRefMCAOD()) {
         m0 = fDaughter[0]->GetRefMCAOD()->GetMother();
         m1 = fDaughter[1]->GetRefMCAOD()->GetMother();
      }
   }
   
   // decide the answer depending on the two mother labels
   if (m0 != m1) 
      return 0;
   else
      return TMath::Abs(fDaughter[0]->GetMotherPDG());
}

//_____________________________________________________________________________
void AliRsnMother::ComputeSum(Double_t mass0, Double_t mass1)
{
//
// Sets the masses for the 4-momenta of the daughters and then
// sums them, taking into account that the space part is set to
// each of them when the reference object is set (see AliRsnDaughter::SetRef)
//

   fDaughter[0]->SetMass(mass0);
   fDaughter[1]->SetMass(mass1);

   fSum   = fDaughter[0]->Prec() + fDaughter[1]->Prec();
   fSumMC = fDaughter[0]->Psim() + fDaughter[1]->Psim();
}

//_____________________________________________________________________________
void AliRsnMother::ResetPair()
{
//
// Resets the mother, zeroing all data members.
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = 0x0;
   fRefEvent = 0x0;

   fSum  .SetXYZM(0.0, 0.0, 0.0, 0.0);
   fSumMC.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
Double_t AliRsnMother::AngleTo(AliRsnDaughter *track, Bool_t mc)
{
//
// Compute the angle betwee this and the passed object
// if second argument is kTRUE, use MC values.
//

   TLorentzVector &me = (mc ? fSumMC : fSum);
   TLorentzVector &he = track->P(mc);
   
   return me.Angle(he.Vect());
}

//_____________________________________________________________________________
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
   Int_t idL = fRefEvent->GetLeadingParticleID();
   
   if (id1 == idL || id2 == idL) {
      success = kFALSE;
      return -99.0;
   }
   
   AliRsnDaughter leading = fRefEvent->GetDaughter(idL);
   AliVParticle  *ref     = leading.GetRef();
   Double_t       angle   = ref->Phi() - fSum.Phi();
   
   //return angle w.r.t. leading particle in the range -pi/2, 3/2pi
   while (angle >= 1.5 * TMath::Pi()) angle -= 2 * TMath::Pi();
   while (angle < -0.5 * TMath::Pi()) angle += 2 * TMath::Pi();
   success = kTRUE;
   
   return angle;
}

//_____________________________________________________________________________
Double_t AliRsnMother::CosThetaStar(Bool_t first, Bool_t useMC)
{
//
// Computes the cosine of theta*, which is the angle of one of the daughters
// with respect to the total momentum of the resonance, in its rest frame.
// The arguments are needed to choose which of the daughters one want to use
// and if reconstructed or MC momentum must be used.
// [Contribution from Z. Feckova]
//

   TLorentzVector mother    = (useMC ? fSumMC : fSum);
   TLorentzVector daughter0 = (first ? fDaughter[0]->P(useMC) : fDaughter[1]->P(useMC));
   TLorentzVector daughter1 = (first ? fDaughter[1]->P(useMC) : fDaughter[0]->P(useMC));
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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
