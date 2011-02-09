//
// Class AliRsnMother
//
// Implementation of a pair of tracks, for several purposes
// - computing the total 4-momentum & inv. mass for output histos filling
// - evaluating cut checks on the pair of particles
// - evaluating any kind of kinematic value over their sum
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#include <Riostream.h>
#include <TVector3.h>
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliRsnDaughter.h"
#include "AliRsnPairDef.h"
#include "AliRsnMother.h"

ClassImp(AliRsnMother)

//_____________________________________________________________________________
AliRsnMother::AliRsnMother() :
   fUseMC(kFALSE),
   fSum(),
   fSumMC()
{
//
// Constructor.
// Initializes all variables to meaningless values.
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = 0x0;
}

//_____________________________________________________________________________
AliRsnMother::AliRsnMother(const AliRsnMother &obj) :
   TObject(obj),
   fUseMC(obj.fUseMC),
   fSum(obj.fSum),
   fSumMC(obj.fSumMC)
{
//
// Copy constructor.
// Initializes all variables to copy values.
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
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

   Int_t i;

   fSum = obj.fSum;
   fSumMC = obj.fSumMC;

   for (i = 0; i < 2; i++) fDaughter[i] = obj.fDaughter[i];

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
// Checks if the two tracks in the pair have the same mother.
// This can be known if MC info is present.
// If the mother label is the same, rhe PDG code of the mother is returned,
// otherwise the method returns 0.
// In the two arguments passed by reference, the mothers of the two daghters are stored
//

   // if MC info is not available, the pairs is not true by default
   if (!fDaughter[0]->GetRefMC() || !fDaughter[1]->GetRefMC()) {
      AliInfo("Cannot know if the pairs is true or not because MC Info is not present");
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
   if (m0 != m1) return 0;

   // if we reach this point, the two tracks have the same mother
   // let's check now the PDG code of this common mother
   return TMath::Abs(fDaughter[0]->GetMotherPDG());
}

//_____________________________________________________________________________
void AliRsnMother::SetDaughters
(AliRsnDaughter *d0, Double_t mass0, AliRsnDaughter *d1, Double_t mass1)
{
//
// Sets the pair defined in this usind tso passed daughters and two masses
// which will be assigned to them, in order to recompute their 4-momenta
// and sum them into the datamembers of this object.
//

   if (d0) fDaughter[0] = d0;
   if (d1) fDaughter[1] = d1;

   if (!d0 || !d1) return;

   fDaughter[0]->SetMass(mass0);
   fDaughter[1]->SetMass(mass1);

   fSum   = fDaughter[0]->P(kFALSE) + fDaughter[1]->P(kFALSE);
   fSumMC = fDaughter[0]->P(kTRUE)  + fDaughter[1]->P(kTRUE);
}

//_____________________________________________________________________________
void AliRsnMother::ResetPair()
{
//
// Resets the mother, nullifying all data members
//

   Int_t i;
   for (i = 0; i < 2; i++) fDaughter[i] = 0x0;

   fSum  .SetXYZM(0.0, 0.0, 0.0, 0.0);
   fSumMC.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
Double_t AliRsnMother::CosThetaStar(Bool_t first, Bool_t useMC)
{
   TLorentzVector mother    = (useMC ? fSumMC : fSum);
   TLorentzVector daughter0 = (first ? fDaughter[0]->P() : fDaughter[1]->P());
   TLorentzVector daughter1 = (first ? fDaughter[1]->P() : fDaughter[0]->P());
   TVector3 momentumM(mother.Vect());
   TVector3 normal(mother.Y() / momentumM.Mag(), -mother.X() / momentumM.Mag(), 0.0);

   // Computes first the invariant mass of the mother
   Double_t mass0            = fDaughter[0]->P().M();
   Double_t mass1            = fDaughter[1]->P().M();
   Double_t p0               = daughter0.Vect().Mag();
   Double_t p1               = daughter1.Vect().Mag();
   Double_t E0               = TMath::Sqrt(mass0 * mass0 + p0 * p0);
   Double_t E1               = TMath::Sqrt(mass1 * mass1 + p1 * p1);
   Double_t MotherMass       = TMath::Sqrt((E0 + E1) * (E0 + E1) - (p0 * p0 + 2.0 * daughter0.Vect().Dot(daughter1.Vect()) + p1 * p1));
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
Bool_t AliRsnMother::CheckPair() const
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

   if (fUseMC) {
      if (fDaughter[0]->GetRefMC() == 0x0 || fDaughter[1]->GetRefMC() == 0x0) {
         AliError("Required MC info but not all MC refs are available");
         return kFALSE;
      }
   }

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnMother::MatchesDef(AliRsnPairDef *def)
{
//
// Checks if the daughters, in any order, do match a given decay channel,
// using the specified identification method, which can be the 'true' one
// or the 'realistic' one only.
//

   if (!def) return kFALSE;
   if (!fDaughter[0]->GetRefMC()) return kFALSE;
   if (!fDaughter[1]->GetRefMC()) return kFALSE;

   Bool_t decayMatch = kFALSE;
   Int_t  pdg[2], ref[2];
   pdg[0] = fDaughter[0]->GetPDG();
   pdg[1] = fDaughter[1]->GetPDG();
   ref[0] = TMath::Abs(AliPID::ParticleCode(def->GetPID(0)));
   ref[1] = TMath::Abs(AliPID::ParticleCode(def->GetPID(1)));

   // check #1:
   // if first member of pairDef has same sign as first member of this,
   // daughter[0] perfect PID must match first member of pairDef
   // daughter[1] perfect PID must march second member of pairDef
   if (fDaughter[0]->IsSign(def->GetCharge(0)) && fDaughter[1]->IsSign(def->GetCharge(1))) {
      decayMatch = (pdg[0] == ref[0] && pdg[1] == ref[1]);
   }

   // check #2:
   // if first member of pairDef has same sign as second member of this,
   // daughter[0] perfect PID must match second member of pairDef
   // daughter[1] perfect PID must march first member of pairDef
   if (fDaughter[1]->IsSign(def->GetCharge(0)) && fDaughter[0]->IsSign(def->GetCharge(1))) {
      decayMatch = (pdg[0] == ref[1] && pdg[1] == ref[0]);
   }

   return (decayMatch && (CommonMother() == def->GetMotherPDG()));
}
