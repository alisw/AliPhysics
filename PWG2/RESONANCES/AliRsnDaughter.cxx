//
// Class AliRsnDaughter
//
// Light-weight AOD object which contains all required track details
// which are used for resonance analysis.
// Provides converters from all kinds of input track type: ESD, AOD and MC.
//
// authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//          M. Vala (martin.vala@cern.ch)
//

#include <Riostream.h>

#include <TParticle.h>
#include <TString.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliMCParticle.h"

#include "AliRsnPIDDefESD.h"
#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

AliRsnDaughter::EPIDMethod AliRsnDaughter::fgPIDMethod = AliRsnDaughter::kRealistic;

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter() :
  AliVParticle(),
  fIndex(-1),
  fLabel(-1),
  fCharge(0),
  fFlags(0),
  fKink(0),
  fMass(0.0),
  fChi2(0.0),
  fNSigmaToVertex(-1.0),
  fITSnum(0),
  fTPCnum(0),
  fRealisticPID(AliRsnPID::kUnknown),
  fMCInfo(0x0)
{
//
// Default constructor.
// Initializes all data-members with meaningless values.
//

  Int_t i;
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    if (i < 3)
    {
      fP[i] = 0.0;
      fV[i] = 0.0;
    }
    fPIDWeight[i] = 0.0;
    fPIDProb[i] = 0.0;
  }
}

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) :
  AliVParticle(copy),
  fIndex(copy.fIndex),
  fLabel(copy.fLabel),
  fCharge(copy.fCharge),
  fFlags(copy.fFlags),
  fKink(copy.fKink),
  fMass(copy.fMass),
  fChi2(copy.fChi2),
  fNSigmaToVertex(copy.fNSigmaToVertex),
  fITSnum(copy.fITSnum),
  fTPCnum(copy.fTPCnum),
  fRealisticPID(copy.fRealisticPID),
  fMCInfo(0x0)
{
//
// Copy constructor.
//

  Int_t i;
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    if (i < 3)
    {
      fP[i] = copy.fP[i];
      fV[i] = copy.fV[i];
    }
    fPIDWeight[i] = copy.fPIDWeight[i];
    fPIDProb[i] = copy.fPIDProb[i];
  }

  // initialize particle object
  // only if it is present in the template object
  if (copy.fMCInfo) fMCInfo = new AliRsnMCInfo(*(copy.fMCInfo));
}

//_____________________________________________________________________________
AliRsnDaughter& AliRsnDaughter::operator=(const AliRsnDaughter &copy)
{
//
// Assignment operator.
// Works like the copy constructor and returns a reference
// to the initialized object for which it is called.
//

  fIndex  = copy.fIndex;
  fLabel  = copy.fLabel;
  fCharge = copy.fCharge;
  fFlags  = copy.fFlags;
  fKink   = copy.fKink;
  fChi2   = copy.fChi2;
  fNSigmaToVertex = copy.fNSigmaToVertex;
  fITSnum = copy.fITSnum;
  fTPCnum = copy.fTPCnum;

  Int_t i;
  for (i = 0; i < AliRsnPID::kSpecies; i++)
  {
    if (i < 3)
    {
      fP[i] = copy.fP[i];
      fV[i] = copy.fV[i];
    }
    fPIDWeight[i] = copy.fPIDWeight[i];
    fPIDProb[i] = copy.fPIDProb[i];
  }

  fMass = copy.fMass;
  fRealisticPID = copy.fRealisticPID;

  // initialize particle object
  // only if it is present in the template object;
  // otherwise, it is just cleared and not replaced with anything
  if (fMCInfo)
  {
    delete fMCInfo;
    fMCInfo = 0x0;
  }
  if (copy.fMCInfo) fMCInfo = new AliRsnMCInfo(*(copy.fMCInfo));

  return (*this);
}

//_____________________________________________________________________________
AliRsnDaughter::~AliRsnDaughter()
{
//
// Destructor
//

  if (fMCInfo)
  {
    delete fMCInfo;
    fMCInfo = 0;
  }
}

//_____________________________________________________________________________
void AliRsnDaughter::RotateP(Double_t angle, Bool_t isDegrees)
{
//
// Rotate the transverse momentum by an angle (in DEGREES)
// around Z axis (does not change the Z component)
//

  if (isDegrees) angle *= TMath::DegToRad();

  Double_t s = TMath::Sin(angle);
  Double_t c = TMath::Cos(angle);
  Double_t xx = fP[0];
  fP[0] = c*xx - s*fP[1];
  fP[1] = s*xx + c*fP[1];
}

//_____________________________________________________________________________
Double_t AliRsnDaughter::AngleTo(AliRsnDaughter *d, Bool_t outInDegrees)
{
//
// Compute angle between the vector momentum of this
// and the one of argument.
//

  Double_t arg, dot, ptot2 = P2() * d->P2();

  if(ptot2 <= 0) {
    return 0.0;
  }
  else {
    dot = Px()*d->Px() + Py()*d->Py() + Pz()*d->Pz();
    arg = dot / TMath::Sqrt(ptot2);
    if (arg >  1.0) arg =  1.0;
    if (arg < -1.0) arg = -1.0;
    if (outInDegrees) return TMath::ACos(arg) * TMath::RadToDeg();
    else return TMath::ACos(arg);
  }
}

//_____________________________________________________________________________
void  AliRsnDaughter::RealisticPID()
{
//
// Assign realistic PID from largest probability
//

  Int_t i, imax = 0;
  Double_t pmax = fPIDProb[0];

  // search for maximum
  for (i = 1; i < AliRsnPID::kSpecies; i++)
  {
    if (fPIDProb[i] > pmax)
    {
      imax = i;
      pmax = fPIDProb[i];
    }
  }

  fRealisticPID = (AliRsnPID::EType)imax;
}

//_____________________________________________________________________________
AliRsnPID::EType AliRsnDaughter::PIDType(Double_t &prob) const
{
//
// Return the PID type according to the selected method
// in the argument passed by reference, the probability is stored.
// It will be realistic for realistic PID and 1 for perfect PID.
//

  switch (fgPIDMethod)
  {
    case kNoPID:
      AliWarning("Requested a PIDtype call in NoPID mode");
      prob = 1.0;
      return AliRsnPID::kUnknown;
    case kPerfect:
      prob = 1.0;
      if (fMCInfo) return AliRsnPID::InternalType(fMCInfo->PDG());
      else return AliRsnPID::kUnknown;
    default:
      if (fRealisticPID >= 0 && fRealisticPID < AliRsnPID::kSpecies)
      {
        prob = fPIDProb[fRealisticPID];
        return fRealisticPID;
      }
      else
      {
        prob = 1.0;
        return AliRsnPID::kUnknown;
      }
  }
}

//_____________________________________________________________________________
void AliRsnDaughter::Print(Option_t *option) const
{
//
// Prints the values of data members, using the options:
// - P --> momentum
// - V --> DCA vertex
// - C --> electric charge
// - F --> flags
// - I --> identification (PID, probability and mass)
// - W --> PID weights
// - M --> Montecarlo (from AliRsnMCInfo)
// - L --> index & label
// - A --> angles
// - ALL --> All oprions switched on
//
// Index and label are printed by default.
//

  TString opt(option);
  opt.ToUpper();

  if (opt.Contains("L") || opt.Contains("ALL"))
  {
    cout << ".......Index            : " << fIndex << endl;
    cout << ".......Label            : " << fLabel << endl;
  }
  if (opt.Contains("P") || opt.Contains("ALL"))
  {
    cout << ".......Px, Py, Pz, Pt   : " << Px() << ' ' << Py() << ' ' << Pz() << ' ' << Pt() << endl;
  }
  if (opt.Contains("A") || opt.Contains("ALL"))
  {
    cout << ".......Phi, Theta       : " << Phi() << ' ' << Theta() << endl;
  }
  if (opt.Contains("V") || opt.Contains("ALL"))
  {
    cout << ".......Vx, Vy, Vz       : " << Xv() << ' ' << Yv() << ' ' << Zv() << endl;
  }
  if (opt.Contains("I") || opt.Contains("ALL"))
  {
    AliRsnPID::EType type;
    Double_t prob;
    type = PIDType(prob);
    cout << ".......PID & prob       : " << AliRsnPID::ParticleName(type) << ' ' << prob << endl;
  }
  if (opt.Contains("C") || opt.Contains("ALL"))
  {
    cout << ".......Charge           : " << fCharge << endl;
  }
  if (opt.Contains("F") || opt.Contains("ALL"))
  {
    cout << ".......Flags            : " << fFlags << endl;
  }
  if (opt.Contains("W") || opt.Contains("ALL"))
  {
    cout << ".......Weights          : ";
    Int_t i;
    for (i = 0; i < AliRsnPID::kSpecies; i++) cout << fPIDWeight[i] << ' ';
    cout << endl;
  }
  if (opt.Contains("M") || opt.Contains("ALL"))
  {
    if (fMCInfo)
    {
      cout << ".......PDG code         : " << fMCInfo->PDG() << endl;
      cout << ".......Mother (label)   : " << fMCInfo->Mother() << endl;
      cout << ".......Mother (PDG code): " << fMCInfo->MotherPDG() << endl;
    }
    else
    {
      cout << ".......MC info not present" << endl;
    }
  }
}

//_____________________________________________________________________________
void AliRsnDaughter::InitMCInfo()
{
//
// Initializes the particle object with default constructor.
//

  fMCInfo = new AliRsnMCInfo;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::InitMCInfo(TParticle *particle)
{
//
// Copies data from an MC particle into the object which
// contains all MC details taken from kinematics info.
// If requested by second argument, momentum and vertex
// of the Particle are copied into the 'fP' and 'fV'
// data members, to simulate a perfect reconstruction.
// If something goes wrong, returns kFALSE,
// otherwise returns kTRUE.
//

  // retrieve the TParticle object pointed by this MC track
  if (!particle)
  {
    AliError("Passed NULL particle object");
    return kFALSE;
  }

  // initialize object if not initialized yet
  if (fMCInfo) delete fMCInfo;
  fMCInfo = new AliRsnMCInfo;
  fMCInfo->Adopt(particle);

  return kTRUE;
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::Compare(const TObject* obj) const
{
//
// Compare two tracks with respect to their transverse momentum.
// Citation from ROOT reference:
// "Must return -1 if this is smaller than obj, 0 if objects are equal
//  and 1 if this is larger than obj".
//

  AliRsnDaughter *that = (AliRsnDaughter*)obj;
  if (Pt() < that->Pt()) return 1;
  else if (Pt() > that->Pt()) return -1;
  else return 0;
}
