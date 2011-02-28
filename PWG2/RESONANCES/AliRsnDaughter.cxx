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
////////////////////////////////////////////////////////////////////////////////

#include <TParticle.h>
#include <TDatabasePDG.h>
#include "AliAODVertex.h"

#include "AliRsnDaughter.h"

ClassImp(AliRsnDaughter)

//_____________________________________________________________________________
AliRsnDaughter::AliRsnDaughter() :
   fOK(kFALSE),
   fLabel(-1),
   fMotherPDG(0),
   fRsnID(-1),
   fPrec(0.0, 0.0, 0.0, 0.0),
   fPsim(0.0, 0.0, 0.0, 0.0),
   fRef(0x0),
   fRefMC(0x0)
{
//
// Default constructor.
// Initializes all data members to the same values
// Which will be given by Reset() method.
//
}

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
   fRefMC(copy.fRefMC)
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

   fOK        = copy.fOK;
   fLabel     = copy.fLabel;
   fMotherPDG = copy.fMotherPDG;
   fRsnID     = copy.fRsnID;
   fPrec      = copy.fPrec;
   fPsim      = copy.fPsim;
   fRef       = copy.fRef;
   fRefMC     = copy.fRefMC;

   return (*this);
}

//_____________________________________________________________________________
void AliRsnDaughter::SetRef(AliVParticle *p)
{
//
// Set the pointer to reference VParticle
// and copies its momentum in the 4-vector.
// Mass is assigned by SetMass().
//

   fRef = p;
   if (p) fPrec.SetXYZT(p->Px(), p->Py(), p->Pz(), 0.0);
   else   fPrec.SetXYZT(0.0, 0.0, 0.0, 0.0);
}

//_____________________________________________________________________________
void AliRsnDaughter::SetRefMC(AliVParticle *p)
{
//
// Set the pointer to reference MonteCarlo VParticle
// and copies its momentum in the 4-vector.
// Mass is assigned by SetMass().
//

   fRefMC = p;
   if (p) fPsim.SetXYZT(p->Px(), p->Py(), p->Pz(), 0.0);
   else   fPrec.SetXYZT(0.0, 0.0, 0.0, 0.0);
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

   SetRef(NULL);
   SetRefMC(NULL);
}


//_____________________________________________________________________________
Int_t AliRsnDaughter::GetPDG(Bool_t abs)
{
//
// Return the PDG code of the particle from MC ref (if any).
// If argument is kTRUE, returns its absolute value.
//

   Int_t pdg = 0;
   if (!fRefMC) return pdg;

   // ESD
   AliMCParticle *esd = GetRefMCESD();
   if (esd) pdg = esd->Particle()->GetPdgCode();

   // AOD
   AliAODMCParticle *aod = GetRefMCAOD();
   if (aod) pdg = aod->GetPdgCode();

   // abs value if required
   if (abs) pdg = TMath::Abs(pdg);
   return pdg;
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
   AliESDtrack *esd = GetRefESDtrack();
   if (esd) return esd->GetID();

   // AOD tracks
   AliAODTrack *aod = GetRefAODtrack();
   if (aod) return aod->GetID();

   // whatever else
   return GetLabel();
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::SetMass(Double_t mass)
{
//
// Assign a mass hypothesis to the track.
// This causes the 4-momentum data members to be re-initialized
// using the momenta of referenced tracks/v0s and this mass.
// This step is fundamental for the following of the analysis.
//

   if (mass < 0.) return kFALSE;

   if (fRef)   fPrec.SetXYZM(fRef  ->Px(), fRef  ->Py(), fRef  ->Pz(), mass);
   if (fRefMC) fPsim.SetXYZM(fRefMC->Px(), fRefMC->Py(), fRefMC->Pz(), mass);

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnDaughter::IsKinkDaughter()
{
//
// Checks if this track is a kink daughter.
// this information is important for some cuts, in some cases
// and it is retrieved differently from ESDs and AODs, so
// this is done here in order to have a unique outcome.
//

   AliESDtrack *etrack = GetRefESDtrack();
   AliAODTrack *atrack = GetRefAODtrack();

   if (etrack) {
      return (etrack->GetKinkIndex(0) > 0);
   } else if (atrack) {
      AliAODVertex *vertex = atrack->GetProdVertex();
      if (vertex) if (vertex->GetType() == AliAODVertex::kKink) return kTRUE;
   }

   return kFALSE;
}

//______________________________________________________________________________
AliRsnDaughter::ERefType AliRsnDaughter::RefType(ESpecies species)
{
//
// Returns the expected object type for a candidate daughter
// of the given species.
//

   switch (species) {
      case kElectron:
      case kMuon:
      case kPion:
      case kKaon:
      case kProton:
         return kTrack;
      case kKaon0:
      case kLambda:
         return kV0;
      case kXi:
      case kOmega:
         return kCascade;
      default:
         return kNoType;
   }
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
