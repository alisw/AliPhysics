//
//  This class works as generic interface to each candidate resonance daughter.
//  Its main purpose is to provide a unique reference which includes all the
//  facilities available in the AliVParticle generic base class, plus all info
//  which could be needed during analysis, which are not in AliVParticle but
//  need to be accessed from ESD or AOD objects, usually in different ways.
//  When MC is available, AliRsnDaughter matches each reconstructed object with
//  its corresponding MC particle.
//  
//  Currently, this interface can point to all kinds of single-particle object
//  which one can have in the reconstructed event: charged tracks, V0s and 
//  cascades. It is care of the user to treat each of them in the correct way,
//  regarding cuts, functions to be computed, etc.
//
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//

#include <TParticle.h>
#include <TDatabasePDG.h>

#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter(const AliRsnDaughter &copy) :
   TObject(copy),
   fOK(copy.fOK),
   fLabel(copy.fLabel),
   fMotherPDG(copy.fMotherPDG),
   fRsnID(copy.fRsnID),
   fPrec(copy.fPrec),
   fPsim(copy.fPsim),
   fRef(copy.fRef),
   fRefMC(copy.fRefMC),
   fOwnerEvent(copy.fOwnerEvent)
{
//
// Copy constructor.
// Pointers are NOT duplicated, since they don't come from a 'new'
// statement, but from just referencing something in the data source.
//
}

//_____________________________________________________________________________
AliRsnDaughter& AliRsnDaughter::operator=(const AliRsnDaughter &copy)
{
//
// Assignment operator.
// Pointers are NOT duplicated, since they don't come from a 'new'
// statement, but from just referencing something in the data source.
//

   fOK         = copy.fOK;
   fLabel      = copy.fLabel;
   fMotherPDG  = copy.fMotherPDG;
   fRsnID      = copy.fRsnID;
   fPsim       = copy.fPsim;
   fPrec       = copy.fPrec;
   fRef        = copy.fRef;
   fRefMC      = copy.fRefMC;
   fOwnerEvent = copy.fOwnerEvent;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnDaughter::Reset()
{
//
// Reset this track to meaningless values and to a 'bad' status.
// After this has been done, this object should not be used
// for analysis unless initialized properly.
//

   fOK        = kFALSE;
   fLabel     = -1;
   fMotherPDG =  0;
   fRsnID     = -1;
   
   fPsim.SetXYZT(0.0, 0.0, 0.0, 0.0);
   fPrec.SetXYZT(0.0, 0.0, 0.0, 0.0);

   fRef = fRefMC = 0x0;
   fOwnerEvent = 0x0;
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetPDG()
{
//
// Return the PDG code of the particle from MC ref (if any).
// If argument is kTRUE, returns its absolute value.
//

   if (Match(fRefMC, AliMCParticle::Class()))
      return ((AliMCParticle*)fRefMC)->Particle()->GetPdgCode();
   else if (Match(fRefMC, AliAODMCParticle::Class()))
      return ((AliAODMCParticle*)fRefMC)->GetPdgCode();
   else {
      AliWarning("Cannot retrieve PDG");
      return 0;
   }
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetID()
{
//
// Return reference index, using the "GetID" method
// of the possible source object.
// When this method is unsuccessful (e.g.: V0s), return the label.
//

   // ESD tracks
   AliESDtrack *esd = Ref2ESDtrack();
   if (esd) return esd->GetID();

   // AOD tracks
   AliAODTrack *aod = Ref2AODtrack();
   if (aod) return aod->GetID();

   // whatever else
   return GetLabel();
}

//_____________________________________________________________________________
Int_t AliRsnDaughter::GetMother()
{
//
// Return index of the first mother of the MC reference, if any.
// Otherwise, returns -1 (the same as for primary tracks)
//

   if (!fRefMC) return -1;

   if (fRefMC->InheritsFrom(AliMCParticle::Class())) {
      AliMCParticle *mc = (AliMCParticle*)fRefMC;
      return mc->Particle()->GetFirstMother();
   } else if (fRefMC->InheritsFrom(AliAODMCParticle::Class())) {
      AliAODMCParticle *mc = (AliAODMCParticle*)fRefMC;
      return mc->GetMother();
   }
   else
      return -1;
}
   
   

//______________________________________________________________________________
void AliRsnDaughter::Print(Option_t *opt) const
{
//
// Override of TObject::Print()
//

   AliInfo("=== DAUGHTER INFO ======================================================================");
   AliInfo(Form(" (sim) px,py,pz = %6.2f %6.2f %6.2f", fPsim.X(), fPsim.Y(), fPsim.Z()));
   AliInfo(Form(" (rec) px,py,pz = %6.2f %6.2f %6.2f", fPrec.X(), fPrec.Y(), fPrec.Z()));
   AliInfo(Form(" OK, RsnID, Label, MotherPDG = %s, %5d, %5d, %4d", (fOK ? "true " : "false"), fRsnID, fLabel, fMotherPDG));
   AliInfo("========================================================================================");
}

//______________________________________________________________________________
const char* AliRsnDaughter::SpeciesName(ESpecies species)
{
//
// Return a string with the short name of the particle
//

   switch (species) {
      case kElectron: return "E";
      case kMuon:     return "Mu";
      case kPion:     return "Pi";
      case kKaon:     return "K";
      case kProton:   return "P";
      case kKaon0:    return "K0s";
      case kLambda:   return "Lambda";
      case kXi:       return "Xi";
      case kOmega:    return "Omega";
      default:        return "Undef";
   }
}

//______________________________________________________________________________
Int_t AliRsnDaughter::SpeciesPDG(ESpecies species)
{
//
// Return the PDG code of a particle species (abs value)
//

   switch (species) {
      case kElectron: return 11;
      case kMuon:     return 13;
      case kPion:     return 211;
      case kKaon:     return 321;
      case kProton:   return 2212;
      case kKaon0:    return 310;
      case kLambda:   return 3122;
      case kXi:       return 3312;
      case kOmega:    return 3334;
      default:        return 0;
   }
}

//______________________________________________________________________________
Double_t AliRsnDaughter::SpeciesMass(ESpecies species)
{
//
// Return the mass of a particle species
//

   TDatabasePDG *db = TDatabasePDG::Instance();
   TParticlePDG *part = 0x0;
   
   Int_t pdg = SpeciesPDG(species);
   if (pdg) {
      part = db->GetParticle(pdg);
      return part->Mass();
   }
   else
      return 0.0;
}

//______________________________________________________________________________
EPARTYPE AliRsnDaughter::ToAliPID(ESpecies species)
{
//
// Convert an enum element from this object
// into the enumeration of AliPID.
// If no match are cound 'kUnknown' is returned.
//

   switch (species) {
      case kElectron: return AliPID::kElectron;
      case kMuon:     return AliPID::kMuon;
      case kPion:     return AliPID::kPion;
      case kKaon:     return AliPID::kKaon;
      case kProton:   return AliPID::kProton;
      case kKaon0:    return AliPID::kKaon0;
      default:        return AliPID::kUnknown;
   }
}

//______________________________________________________________________________
AliRsnDaughter::ESpecies AliRsnDaughter::FromAliPID(EPARTYPE pid)
{
//
// Convert an enum element from AliPID
// into the enumeration of this object.
// If no match are cound 'kUnknown' is returned.
//

   switch (pid) {
      case AliPID::kElectron: return kElectron;
      case AliPID::kMuon:     return kMuon;
      case AliPID::kPion:     return kPion;
      case AliPID::kKaon:     return kKaon;
      case AliPID::kProton:   return kProton;
      case AliPID::kKaon0:    return kKaon0;
      default:                return kUnknown;
   }
}
