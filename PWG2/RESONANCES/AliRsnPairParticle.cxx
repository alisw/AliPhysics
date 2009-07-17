//
// Class AliRsnPairParticle
//
// Implementation of a pair of tracks, for several purposes
// - computing the total 4-momentum & inv. mass for output histos filling
// - evaluating cut checks on the pair of particles
// - evaluating any kind of kinematic value over their sum
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliRsnPairDef.h"
#include "AliRsnPairParticle.h"

ClassImp(AliRsnPairParticle)

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle()
{
//
// Constructor.
// Initializes all variables to meaningless values.
//

  Int_t i, j;

  for (i = 0; i < 3; i++) {
    fPTot[i] = 0.0;
    fPTotMC[i] = 0.0;
    if (i < 2) {
      fMotherLabel[i] = -1;
      fMotherPDG[i] = 0;
      fDaughter[i] = 0x0;
    }
    for (j = 0; j < 2; j++) {
      fPTrack[j][i] = 0.0;
      fPTrackMC[j][i] = 0.0;
    }
  }
}

//_____________________________________________________________________________
AliRsnPairParticle::AliRsnPairParticle(const AliRsnPairParticle &obj) :
    TObject(obj)
{
//
// Copy constructor.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

  Int_t i, j;
  for (i = 0; i < 3; i++) {
    fPTot[i] = obj.fPTot[i];
    fPTotMC[i] = obj.fPTotMC[i];
    if (i < 2) {
      fMotherLabel[i] = obj.fMotherLabel[i];
      fMotherPDG[i] = obj.fMotherPDG[i];
      fDaughter[i] = obj.fDaughter[i];
    }
    for (j = 0; j < 2; j++) {
      fPTrack[j][i] = obj.fPTrack[j][i];
      fPTrackMC[j][i] = obj.fPTrackMC[j][i];
    }
  }
}

//_____________________________________________________________________________
AliRsnPairParticle& AliRsnPairParticle::operator=(const AliRsnPairParticle &obj)
{
//
// Assignment operator.
// Initializes all variables to copy values.
// Does not duplicate pointers.
//

  Int_t i, j;
  for (i = 0; i < 3; i++) {
    fPTot[i] = obj.fPTot[i];
    fPTotMC[i] = obj.fPTotMC[i];
    if (i < 2) {
      fMotherLabel[i] = obj.fMotherLabel[i];
      fMotherPDG[i] = obj.fMotherPDG[i];
      fDaughter[i] = obj.fDaughter[i];
    }
    for (j = 0; j < 2; j++) {
      fPTrack[j][i] = obj.fPTrack[j][i];
      fPTrackMC[j][i] = obj.fPTrackMC[j][i];
    }
  }

  return (*this);
}

//_____________________________________________________________________________
AliRsnPairParticle::~AliRsnPairParticle()
{
//
// Desctructor.
// Does nothing.
//
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMass(Double_t mass0, Double_t mass1)
{
//
// Compute invariant mass using reconstructed values.
// Mass in argument #1 is assigned to first track in the pair (fDaughter[0]),
// mass in argument #2 is assigned to second track in the pair (fDaughter[1]).
// If the third argument is not zero, a rotation of that angle is applied to XY momentum of track #2
// before computing invariant mass.
// Then, the invariant mass of the pair is computed by using their total momentum
// and the sum of their energies computed after assigned masses.
//

  if (!fDaughter[0] || !fDaughter[1]) {
    AliError("One of the two tracks is NULL in this pair!");
    return -1000.0;
  }

  // compute track energies using the shortcut method defined in AliRsnDaughter
  Double_t etot = 0.0;
  etot += fDaughter[0]->E(mass0);
  etot += fDaughter[1]->E(mass1);

  // compute & return invariant mass
  return  TMath::Sqrt(etot * etot - GetP2());
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetInvMassMC(Double_t mass0, Double_t mass1)
{
//
// Compute invariant mass using MC values.
// Mass in argument #1 is assigned to first track in the pair (fDaughter[0]),
// mass in argument #2 is assigned to second track in the pair (fDaughter[1])
// Then, the invariant mass of the pair is computed by using their total momentum
// and the sum of their energies.
//

  if (!fDaughter[0] || !fDaughter[1]) {
    AliError("One of the two tracks is NULL in this pair!");
    return -1000.0;
  }
  if (!fDaughter[0]->GetParticle() || !fDaughter[1]->GetParticle()) {
    AliError("One of the two tracks has a NULL MCInfo in this pair!");
    return -1000.0;
  }

  // compute track energies using the shortcut method defined in AliRsnDaughter
  Double_t etot = 0.0;
  etot += fDaughter[0]->GetMCEnergy(mass0);
  etot += fDaughter[1]->GetMCEnergy(mass1);

  // compute & return invariant mass
  return  TMath::Sqrt(etot * etot - GetP2MC());
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetEtot(Double_t mass0, Double_t mass1) const
{
//
// Compute total pair energy from the sum of single track energies
// with a necessary mass hypothesis (rec. values).
//

  Double_t etot = 0.0;
  etot += fDaughter[0]->E(mass0);
  etot += fDaughter[1]->E(mass1);

  return etot;
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetEtotMC(Double_t mass0, Double_t mass1) const
{
//
// Compute total pair energy from the sum of single track energies
// with a necessary mass hypothesis (MC values).
//
  Double_t etot = 0.0;
  etot += fDaughter[0]->GetMCEnergy(mass0);
  etot += fDaughter[1]->GetMCEnergy(mass1);

  return etot;
}

//_____________________________________________________________________________
Double_t AliRsnPairParticle::GetAngle() const
{
//
// Returns the relative angle between the vector momenta of the tracks
// Return value is in DEGREES.
//

  Double_t dotProd = 0.0;
  dotProd += fDaughter[0]->Px() * fDaughter[1]->Px();
  dotProd += fDaughter[0]->Py() * fDaughter[1]->Py();
  dotProd += fDaughter[0]->Pz() * fDaughter[1]->Pz();

  Double_t cosAngle = dotProd / (fDaughter[0]->P() * fDaughter[1]->P());

  return TMath::ACos(cosAngle) * TMath::RadToDeg();
}

//_____________________________________________________________________________
Bool_t AliRsnPairParticle::IsTruePair(Int_t refPDG)
{
//
// Checks if the two tracks in the pair come from the same resonance.
// This can be known if MC info is present, looking at the GEANT label of mother
// (which should be the same).
// If the argument is 0, the answer is kTRUE whenever the labels of mothers of the
// two tracks is the same. When the argument is not zero, the answer is kTRUE only
// if the mother is the same and its PDG code is equal to the argument.
//

  // if MC info is not available, the pairs is not true by default
  if (!fDaughter[0]->GetParticle() || !fDaughter[1]->GetParticle()) {
    AliInfo("Cannot know if the pairs is true or not because MC Info is not present");
    return kFALSE;
  }

  // check that labels are the same
  if (fDaughter[0]->GetParticle()->GetFirstMother() != fDaughter[1]->GetParticle()->GetFirstMother()) {
    return kFALSE;
  }

  // if we reach this point, the two tracks have the same mother
  // let's check now the PDG code of this common mother
  Int_t motherPDG = TMath::Abs(fDaughter[0]->GetMotherPDG());
  if (refPDG == 0) return kTRUE;
  else return (motherPDG == refPDG);
}

//_____________________________________________________________________________
Int_t AliRsnPairParticle::CommonMother()
{
//
// Checks if the two tracks in the pair have the same mother.
// This can be known if MC info is present.
// If the mother label is the same, rhe PDG code of the mother is returned,
// otherwise the method returns 0.
//

  // if MC info is not available, the pairs is not true by default
  if (!fDaughter[0]->GetParticle() || !fDaughter[1]->GetParticle()) {
    AliInfo("Cannot know if the pairs is true or not because MC Info is not present");
    return 0;
  }

  // check that labels are the same
  if (fDaughter[0]->GetParticle()->GetFirstMother() != fDaughter[1]->GetParticle()->GetFirstMother()) {
    return 0;
  }

  // if we reach this point, the two tracks have the same mother
  // let's check now the PDG code of this common mother
  return TMath::Abs(fDaughter[0]->GetMotherPDG());
}

//_____________________________________________________________________________
void AliRsnPairParticle::SetPair(AliRsnDaughter * const daughter1, AliRsnDaughter * const daughter2)
{
//
// Accepts two AliRsnDaughter's which are the two tracks in the pair,
// fills all data-members which contain their momenta & info,
// and computes the total momentum for REC data and MC if available
//

  Int_t i;

  fDaughter[0] = daughter1;
  fDaughter[1] = daughter2;

  // copy MC info (if available)
  if (fDaughter[0]->GetParticle() && fDaughter[1]->GetParticle()) {
    for (i = 0; i < 2; i++) {
      fPTrackMC[i][0] = fDaughter[i]->GetParticle()->Px();
      fPTrackMC[i][1] = fDaughter[i]->GetParticle()->Py();
      fPTrackMC[i][2] = fDaughter[i]->GetParticle()->Pz();
      fMotherPDG[i] = fDaughter[i]->GetMotherPDG();
    }
    for (i = 0; i < 3; i++) fPTotMC[i] = fPTrackMC[0][i] + fPTrackMC[1][i];
  }

  // copy reconstructed info (always available)
  for (i = 0; i < 2; i++) {
    fPTrack[i][0] = fDaughter[i]->Px();
    fPTrack[i][1] = fDaughter[i]->Py();
    fPTrack[i][2] = fDaughter[i]->Pz();
  }
  for (i = 0; i < 3; i++) fPTot[i] = fPTrack[0][i] + fPTrack[1][i];
}

//_____________________________________________________________________________
void AliRsnPairParticle::ResetPair()
{
//
// Computes the total momentum for REC data and MC if available
//

  Int_t i;

  // copy MC info (if available)
  if (fDaughter[0]->GetParticle() && fDaughter[1]->GetParticle()) {
    for (i = 0; i < 2; i++) {
      fPTrackMC[i][0] = fDaughter[i]->GetParticle()->Px();
      fPTrackMC[i][1] = fDaughter[i]->GetParticle()->Py();
      fPTrackMC[i][2] = fDaughter[i]->GetParticle()->Pz();
      fMotherPDG[i] = fDaughter[i]->GetMotherPDG();
    }
    for (i = 0; i < 3; i++) fPTotMC[i] = fPTrackMC[0][i] + fPTrackMC[1][i];
  }

  // copy reconstructed info (always available)
  for (i = 0; i < 2; i++) {
    fPTrack[i][0] = fDaughter[i]->Px();
    fPTrack[i][1] = fDaughter[i]->Py();
    fPTrack[i][2] = fDaughter[i]->Pz();
  }
  for (i = 0; i < 3; i++) fPTot[i] = fPTrack[0][i] + fPTrack[1][i];
}

//_____________________________________________________________________________
void AliRsnPairParticle::RotateTrack(Int_t idx, Double_t angle, Bool_t isDegrees)
{
//
// Computes the total momentum for REC data and MC if available
//

  if (idx < 0 || idx > 1) return;

  Int_t    i;

  // copy MC info (if available)
  if (fDaughter[0]->GetParticle() && fDaughter[1]->GetParticle()) {
    for (i = 0; i < 2; i++) {
      if (i == idx) {
        fDaughter[i]->RotateP(angle, fPTrackMC[i][0], fPTrackMC[i][1], isDegrees);
      } else {
        fPTrackMC[i][0] = fDaughter[i]->GetParticle()->Px();
        fPTrackMC[i][1] = fDaughter[i]->GetParticle()->Py();
      }
      fPTrackMC[i][2] = fDaughter[i]->GetParticle()->Pz();
    }
    for (i = 0; i < 3; i++) fPTotMC[i] = fPTrackMC[0][i] + fPTrackMC[1][i];
  }

  for (i = 0; i < 2; i++) {
    if (i == idx) {
      fDaughter[i]->RotateP(angle, fPTrack[i][0], fPTrack[i][1], isDegrees);
    } else {
      fPTrack[i][0] = fDaughter[i]->GetParticle()->Px();
      fPTrack[i][1] = fDaughter[i]->GetParticle()->Py();
    }
    fPTrack[i][2] = fDaughter[i]->GetParticle()->Pz();
  }
  for (i = 0; i < 3; i++) fPTot[i] = fPTrack[0][i] + fPTrack[1][i];

}

//_____________________________________________________________________________
void AliRsnPairParticle::PrintInfo(const Option_t *option)
{
//
// Print some info of the pair.
// The options are passed to the AliRsnDaughter::Print() method
//
  option = "ALL";
  AliInfo("======== BEGIN PAIR INFO ===========");
  AliInfo("Track #1");
  fDaughter[0]->Print(option);
  AliInfo("Track #2");
  fDaughter[1]->Print(option);
  AliInfo("========= END PAIR INFO ===========");
}

//_____________________________________________________________________________
Bool_t AliRsnPairParticle::MatchesDef(AliRsnPairDef *def)
{
//
// Checks if the daughters, in any order, do match a given decay channel,
// using the specified identification method, which can be the 'true' one
// or the 'realistic' one only.
//

  if (!def) return kFALSE;

  Bool_t decayMatch = kFALSE;

  // check #1:
  // daughter[0] matches first member of pairDef
  // daughter[1] matches second member of pairDef
  if (fDaughter[0]->IsSign(def->GetCharge(0)) && fDaughter[1]->IsSign(def->GetCharge(1))) {
    decayMatch = (fDaughter[0]->IsPerfectPID(def->GetType(0)) && fDaughter[1]->IsPerfectPID(def->GetType(1)));
    fDaughter[0]->SetRequiredPID(def->GetType(0));
    fDaughter[1]->SetRequiredPID(def->GetType(1));
  }

  // check #2:
  // daughter[1] matches first member of pairDef
  // daughter[0] matches second member of pairDef
  if (fDaughter[1]->IsSign(def->GetCharge(0)) && fDaughter[0]->IsSign(def->GetCharge(1))) {
    decayMatch = (fDaughter[1]->IsPerfectPID(def->GetType(0)) && fDaughter[0]->IsPerfectPID(def->GetType(1)));
    fDaughter[1]->SetRequiredPID(def->GetType(0));
    fDaughter[0]->SetRequiredPID(def->GetType(1));
  }

  return (decayMatch && (CommonMother() == def->GetMotherPDG()));
}
