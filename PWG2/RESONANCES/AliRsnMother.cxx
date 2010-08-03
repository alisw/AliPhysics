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
#include "AliRsnDaughter.h"
#include "AliRsnPairDef.h"
#include "AliRsnMother.h"

ClassImp(AliRsnMother)

//_____________________________________________________________________________
AliRsnMother::AliRsnMother() : 
  fUseMC(kFALSE),
  fDefaultMass(0.0),
  fSum(),
  fSumMC(),
  fRef(),
  fRefMC()
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
  fDefaultMass(obj.fDefaultMass),
  fSum(obj.fSum),
  fSumMC(obj.fSumMC),
  fRef(obj.fRef),
  fRefMC(obj.fRefMC)
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
  
  fDefaultMass = obj.fDefaultMass;
  fSum = obj.fSum;
  fRef = obj.fRef;
  fSumMC = obj.fSumMC;
  fRefMC = obj.fRefMC;
  
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
Int_t AliRsnMother::CommonMother() const
{
//
// Checks if the two tracks in the pair have the same mother.
// This can be known if MC info is present.
// If the mother label is the same, rhe PDG code of the mother is returned,
// otherwise the method returns 0.
//

  // if MC info is not available, the pairs is not true by default
  if (!fDaughter[0]->GetRefMC() || !fDaughter[1]->GetRefMC()) 
  {
    AliInfo("Cannot know if the pairs is true or not because MC Info is not present");
    return 0;
  }

  // check that labels are the same
  if (fDaughter[0]->GetParticle()->GetFirstMother() != fDaughter[1]->GetParticle()->GetFirstMother())
    return 0;

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
  
  fRef  .SetXYZM(fSum  .X(), fSum  .Y(), fSum  .Z(), fDefaultMass);
  fRefMC.SetXYZM(fSumMC.X(), fSumMC.Y(), fSumMC.Z(), fDefaultMass);
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
  fRef  .SetXYZM(0.0, 0.0, 0.0, 0.0);
  fSumMC.SetXYZM(0.0, 0.0, 0.0, 0.0);
  fRefMC.SetXYZM(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
Double_t AliRsnMother::ThetaStar(Bool_t first, Bool_t useMC)
{
//
// Returns the theta* as the angle of the first daughter
// w.r. to the mother momentum, in its rest frame
//

  TLorentzVector &mother   = (useMC ? fSumMC : fSum);
  TLorentzVector &daughter = (first ? fDaughter[0]->P() : fDaughter[1]->P());

  Double_t beta  = mother.Beta();
  Double_t gamma = 1.0 / TMath::Sqrt(1.0 - beta*beta);
  Double_t angle = daughter.Angle(mother.Vect());
  Double_t pproj = daughter.Mag() * TMath::Cos(angle);
  
  Double_t plstar = gamma * (pproj - beta*daughter.E());
  Double_t ptstar = daughter.Mag() * TMath::Sin(angle);
  
  return TMath::ATan(ptstar / plstar);
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

  if (!fDaughter[0] || !fDaughter[1]) 
  {
    AliError("One of the two tracks is NULL in this pair!");
    return kFALSE;
  }
  
  if (fUseMC)
  {
    if (fDaughter[0]->GetRefMC() == 0x0 || fDaughter[1]->GetRefMC() == 0x0)
    {
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
  pdg[0] = fDaughter[0]->GetRefMC()->Particle()->GetPdgCode();
  pdg[1] = fDaughter[1]->GetRefMC()->Particle()->GetPdgCode();
  ref[0] = TMath::Abs(AliPID::ParticleCode(def->GetPID(0)));
  ref[1] = TMath::Abs(AliPID::ParticleCode(def->GetPID(1)));

  // check #1:
  // if first member of pairDef has same sign as first member of this,
  // daughter[0] perfect PID must match first member of pairDef
  // daughter[1] perfect PID must march second member of pairDef
  if (fDaughter[0]->IsSign(def->GetCharge(0)) && fDaughter[1]->IsSign(def->GetCharge(1))) 
  {
    decayMatch = (pdg[0] == ref[0] && pdg[1] == ref[1]);
  }

  // check #2:
  // if first member of pairDef has same sign as second member of this,
  // daughter[0] perfect PID must match second member of pairDef
  // daughter[1] perfect PID must march first member of pairDef
  if (fDaughter[1]->IsSign(def->GetCharge(0)) && fDaughter[0]->IsSign(def->GetCharge(1))) 
  {
    decayMatch = (pdg[0] == ref[1] && pdg[1] == ref[0]);
  }

  return (decayMatch && (CommonMother() == def->GetMotherPDG()));
}
