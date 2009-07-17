//
// Class AliRsnDaughter
//
// Interface to candidate daughters of a resonance (tracks).
// Points to the source of information, which is generally an AliVParticle-derived object
// and contains few internal data-members to store "on fly" some important information
// for the computations required during resonance analysis.
// ---
// Since the package revision, this object is not supposed to be stacked in memory
// but created "on fly" during analysis and used just for computations, as an interface.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#include <Riostream.h>
#include <TParticle.h>

#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

#include "AliRsnPIDDefESD.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

AliRsnDaughter::EPIDMethod AliRsnDaughter::fgPIDMethod = AliRsnDaughter::kRealistic;

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(AliVParticle *ref, TParticle *refMC) :
    fOK((ref != 0)),
    fKinkIndex(0),
    fParticle(refMC),
    fMotherPDG(0),
    fStatus(0),
    fDr(0.0),
    fDz(0.0),
    fReqPID(AliPID::kUnknown),
    fRef(ref)
{
//
// Default constructor.
//
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) :
    TObject(copy),
    fOK(copy.fOK),
    fKinkIndex(copy.fKinkIndex),
    fParticle(copy.fParticle),
    fMotherPDG(copy.fMotherPDG),
    fStatus(copy.fStatus),
    fDr(copy.fDr),
    fDz(copy.fDz),
    fReqPID(copy.fReqPID),
    fRef(copy.fRef)
{
//
// Copy constructor.
// Pointers are NOT duplicated.
//
}

//_____________________________________________________________________________
AliRsnDaughter& AliRsnDaughter::operator=(const AliRsnDaughter &copy)
{
//
// Assignment operator.
//

  (TObject)(*this) = (TObject)copy;

  fOK = copy.fOK;
  fKinkIndex = copy.fKinkIndex;
  fParticle  = copy.fParticle;
  fMotherPDG = copy.fMotherPDG;
  fStatus = copy.fStatus;
  fDr = copy.fDr;
  fDz = copy.fDz;

  fReqPID = copy.fReqPID;

  fRef = copy.fRef;

  return (*this);
}

//_____________________________________________________________________________
AliRsnDaughter::~AliRsnDaughter()
{
//
// Destructor.
// Since pointers do not allocate new objects, nothing is done.
//
}

//_____________________________________________________________________________
void AliRsnDaughter::RotateP
(Double_t angle, Double_t &x, Double_t &y, Bool_t isDegrees)
{
//
// Rotate the transverse momentum by an angle (in DEGREES)
// around Z axis (does not change the Z component).
// Rotated values are stored in the two arguments passed by reference.
//

  if (isDegrees) angle *= TMath::DegToRad();

  Double_t s = TMath::Sin(angle);
  Double_t c = TMath::Cos(angle);

  x = c*Px() - s*Py();
  y = s*Px() + c*Py();
}

//_____________________________________________________________________________
Double_t AliRsnDaughter::AngleTo(AliRsnDaughter d, Bool_t outInDegrees)
{
//
// Compute angle between the vector momentum of this
// and the one of argument.
//

  Double_t arg, dot, ptot2 = P2() * d.P2();

  if (ptot2 <= 0) {
    return 0.0;
  } else {
    dot = Px()*d.Px() + Py()*d.Py() + Pz()*d.Pz();
    arg = dot / TMath::Sqrt(ptot2);
    if (arg >  1.0) arg =  1.0;
    if (arg < -1.0) arg = -1.0;
    if (outInDegrees) return TMath::ACos(arg) * TMath::RadToDeg();
    else return TMath::ACos(arg);
  }
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetID() const
{
//
// Return reference index, using the "GetID" method
// of the possible source object.
//

  AliESDtrack *esd = dynamic_cast<AliESDtrack*>(fRef);
  if (esd) return esd->GetID();

  AliAODTrack *aod = dynamic_cast<AliAODTrack*>(fRef);
  if (aod) return aod->GetID();

  return GetLabel();
}

//_____________________________________________________________________________
AliPID::EParticleType AliRsnDaughter::RealisticPID() const
{
//
// Return the "realistic" PID of this track,
// i.e. the particle species to which corresponds the largest PID probability.
//

  AliPID::EParticleType pid = AliPID::kElectron;
  Double_t prob = fPID[0];

  Int_t i;
  for (i = 1; i < AliPID::kSPECIES; i++) {
    if (fPID[i] > prob) {
      prob = fPID[i];
      pid = (AliPID::EParticleType)i;
    }
  }

  return pid;
}

//_____________________________________________________________________________
AliPID::EParticleType AliRsnDaughter::PerfectPID() const
{
//
// Return the "perfect" PID of this track,
// reading it from the MC information, if available.
//

  if (!fParticle) return AliPID::kUnknown;

  Int_t absPDG = TMath::Abs(fParticle->GetPdgCode());
  switch (absPDG) {
  case   11:
    return AliPID::kElectron;
  case   13:
    return AliPID::kMuon;
  case  211:
    return AliPID::kPion;
  case  321:
    return AliPID::kKaon;
  case 2212:
    return AliPID::kProton;
  default:
    AliDebug(2, Form("PDG code = %d not recognized. Return 'AliPID::kUnknown'", absPDG));
    return AliPID::kUnknown;
  }
}

//_____________________________________________________________________________
AliPID::EParticleType AliRsnDaughter::PIDType(Double_t &prob) const
{
//
// Return the PID type according to the selected method
// in the argument passed by reference, the probability is stored.
// It will be realistic for realistic PID and 1 for perfect PID.
//

  AliPID::EParticleType pid = AssignedPID();

  prob = 1.0;
  if (fgPIDMethod == kRealistic) prob = fPID[(Int_t)pid];

  return pid;
}

//_____________________________________________________________________________
AliPID::EParticleType AliRsnDaughter::AssignedPID() const
{
//
// Return the PID type according to the selected method
// in the argument passed by reference, the probability is stored.
// It will be realistic for realistic PID and 1 for perfect PID.
//

  switch (fgPIDMethod) {
  case kNoPID:
    return AliPID::kUnknown;
  case kPerfect:
    return PerfectPID();
  case kRealistic:
    return RealisticPID();
  default:
    AliWarning("PID method not properly set. Returning realistic PID");
    return RealisticPID();
  }
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::CombineWithPriors(const Double_t *priors, AliRsnPIDDefESD *pidDef)
{
//
// Combine current PID weights (assumed to be them) with prior probs
//

  Int_t       i;
  Double_t    sum = 0.0;

  // get PID weights according to definition
  // if the reference is not ESD or the pidDef is null
  // of it is not null but is requires the ESD pid,
  // the standard PID value is used, otherwise
  if (pidDef && GetRefESD()) {
    pidDef->ComputeWeights(GetRefESD(), fPID);
  } else {
    if (fRef->PID())
      for (i = 0; i < AliPID::kSPECIES; i++) fPID[i] = fRef->PID()[i];
    else if (fParticle) {
      // if the PID() returns 0x0, this is a MC particle
      // and we set the weight corresponding to its species to 1
      // and the others to zero
      for (i = 0; i < AliPID::kSPECIES; i++) fPID[i] = 0.0;
      i = (Int_t)AliRsnDaughter::InternalType(TMath::Abs(fParticle->GetPdgCode()));
      fPID[i] = 1.0;
    }
  }

  // multiply weights and priors
  for (i = 0; i < AliPID::kSPECIES; i++) {
    fPID[i] *= priors[i];
    sum += fPID[i];
  }
  if (sum <= 0.0) {
    if (sum < 0.0) AliError(Form("Sum of weights = %f < 0", sum));
    return kFALSE;
  }

  // normalize
  for (i = 0; i < AliPID::kSPECIES; i++) fPID[i] /= sum;

  return kTRUE;
}

//_____________________________________________________________________________
void AliRsnDaughter::FindMotherPDG(AliStack *const stack)
{
//
// Searches the stack to find the mother and retrieve its PDG code.
//

  if (!stack || !fParticle) return;

  const Int_t mLabel = fParticle->GetFirstMother();
  if (mLabel < 0) {
    fMotherPDG = 0;
  } else {
    TParticle *mum = stack->Particle(mLabel);
    if (mum) fMotherPDG = mum->GetPdgCode();
    else fMotherPDG = 0;
  }
}

//_____________________________________________________________________________
Double_t AliRsnDaughter::GetMCEnergy(Double_t mass)
{
//
// Uses the argument to compute 4-momentum energy
//

  if (!fParticle) return 0.0;

  Double_t p2 = fParticle->Px()*fParticle->Px();
  p2 += fParticle->Py()*fParticle->Py();
  p2 += fParticle->Pz()*fParticle->Pz();

  return TMath::Sqrt(mass*mass + p2);
}

//_____________________________________________________________________________
void AliRsnDaughter::FindKinkIndex(const AliESDtrack *esdTrack)
{
//
// Assign kink index from an ESD track
//

  Int_t i, ik[3];
  for (i = 0; i < 3; i++) ik[i] = esdTrack->GetKinkIndex(i);

  if (ik[0] < 0 || ik[1] < 0 || ik[2] < 0) {
    SetKinkMother();
  } else if (ik[0] > 0 || ik[1] > 0 || ik[2] > 0) {
    SetKinkDaughter();
  } else SetNoKink();
}

//_____________________________________________________________________________
void AliRsnDaughter::FindKinkIndex(AliAODEvent *const event)
{
//
// Assign kink index from an AOD event
//

  Int_t iv, id, nD, nV = event->GetNumberOfVertices();
  for (iv = 0; iv < nV; iv++) {
    AliAODVertex *v = event->GetVertex(iv);
    AliAODVertex::AODVtx_t type = (AliAODVertex::AODVtx_t)v->GetType();
    if (type != AliAODVertex::kKink) continue;
    AliAODTrack *mother = (AliAODTrack*)v->GetParent();
    if (mother == (AliAODTrack*)fRef) {
      SetKinkMother();
      return;
    } else {
      nD = v->GetNDaughters();
      for (id = 0; id < nD; id++) {
        AliAODTrack *son = (AliAODTrack*)v->GetDaughter(id);
        if (son == (AliAODTrack*)fRef) {
          SetKinkDaughter();
          return;
        }
      }
    }
  }

  SetNoKink();
}

//_____________________________________________________________________________
void AliRsnDaughter::Reset()
{
//
// Reset this track to meaningless values
//

  fOK = kFALSE;
  fKinkIndex = 0;
  fParticle = 0x0;
  fMotherPDG = 0;
  fStatus = 0;
  fRef = 0x0;
}

//_____________________________________________________________________________
void AliRsnDaughter::Print(Option_t * const option) const
{
//
// Prints the values of data members, using the options:
// - P --> momentum
// - V --> DCA vertex
// - C --> electric charge
// - F --> flags
// - I --> identification (PID, probability and mass)
// - W --> PID weights
// - M --> Montecarlo
// - L --> index & label
// - A --> angles
// - ALL --> All oprions switched on
//
// Index and label are printed by default.
//

  TString opt(option);
  opt.ToUpper();

  if (opt.Contains("L") || opt.Contains("ALL")) {
    cout << ".......Index            : " << GetID() << endl;
    cout << ".......Label            : " << GetLabel() << endl;
  }
  if (opt.Contains("P") || opt.Contains("ALL")) {
    cout << ".......Px, Py, Pz, Pt   : " << Px() << ' ' << Py() << ' ' << Pz() << ' ' << Pt() << endl;
  }
  if (opt.Contains("A") || opt.Contains("ALL")) {
    cout << ".......Phi, Theta       : " << Phi() << ' ' << Theta() << endl;
  }
  if (opt.Contains("V") || opt.Contains("ALL")) {
    cout << ".......Vx, Vy, Vz       : " << Xv() << ' ' << Yv() << ' ' << Zv() << endl;
  }
  if (opt.Contains("I") || opt.Contains("ALL")) {
    AliPID::EParticleType type;
    Double_t prob;
    type = PIDType(prob);
    cout << ".......PID & prob       : " << AliPID::ParticleName(type) << ' ' << prob << endl;
  }
  if (opt.Contains("C") || opt.Contains("ALL")) {
    cout << ".......Charge           : " << Charge() << endl;
  }
  if (opt.Contains("F") || opt.Contains("ALL")) {
    cout << ".......Flags            : " << fStatus << endl;
  }
  if (opt.Contains("W") || opt.Contains("ALL")) {
    cout << ".......Weights          : ";
    Int_t i;
    for (i = 0; i < AliPID::kSPECIES; i++) cout << fPID[i] << ' ';
    cout << endl;
  }
  if (opt.Contains("M") || opt.Contains("ALL")) {
    if (fParticle) {
      cout << ".......PDG code         : " << fParticle->GetPdgCode() << endl;
      cout << ".......Mother (label)   : " << fParticle->GetFirstMother() << endl;
      cout << ".......Mother (PDG code): " << fMotherPDG << endl;
    } else {
      cout << ".......MC info not present" << endl;
    }
  }
}

//_____________________________________________________________________________
AliPID::EParticleType AliRsnDaughter::InternalType(Int_t pdg)
{
//
// Return the internal enum value corresponding to the PDG
// code passed as argument, if possible.
// Otherwise, returns 'AliPID::kSPECIES' by default.
//

  AliPID::EParticleType value;
  Int_t absPDG = TMath::Abs(pdg);

  switch (absPDG) {
  case 11:
    value = AliPID::kElectron;
    break;
  case 13:
    value = AliPID::kMuon;
    break;
  case 211:
    value = AliPID::kPion;
    break;
  case 321:
    value = AliPID::kKaon;
    break;
  case 2212:
    value = AliPID::kProton;
    break;
  default:
    value = AliPID::kUnknown;
  }
  return value;
}

//_____________________________________________________________________________
const char* AliRsnDaughter::MethodName(EPIDMethod method)
{
//
// Returns a string with the method name
//

  switch (method)
  {
    case kNoPID:
      return "No PID";
    case kRealistic:
      return "Realistic";
    case kPerfect:
      return "Perfect";
    default:
      return "Unknown";
  }
}
