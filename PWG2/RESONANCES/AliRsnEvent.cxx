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
//  This class works as generic interface to an event.
//  Its main purpose is to provide a unique reference which includes all the
//  facilities available in the AliVEvent generic base class, plus all info
//  which could be needed during analysis, which are not in AliVEvent but
//  need to be accessed from ESD or AOD objects, usually in different ways.
//  When MC is available, it is properly taken into account.
//  
//  authors: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//           M. Vala (martin.vala@cern.ch)
//
////////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TArrayF.h>

#include "AliGenEventHeader.h"

#include "AliRsnCutSet.h"
#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(AliVEvent *ref, AliVEvent *refMC) :
   fRef(ref),
   fRefMC(refMC),
   fLeading(-1),
   fPID(0x0),
   fAODList(0x0)
{
//
// Default constructor.
//
}

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(const AliRsnEvent &event) :
   TObject(event),
   fRef(event.fRef),
   fRefMC(event.fRefMC),
   fLeading(event.fLeading),
   fPID(event.fPID),
   fAODList(event.fAODList)
{
//
// Copy constructor.
//
}

//_____________________________________________________________________________
AliRsnEvent& AliRsnEvent::operator= (const AliRsnEvent & event)
{
//
// Works in the same way as the copy constructor.
//

   TObject::operator=(event);
   fRef             = event.fRef;
   fRefMC           = event.fRefMC;
   fLeading         = event.fLeading;
   fPID             = event.fPID;
   fAODList         = event.fAODList;

   return (*this);
}

//_____________________________________________________________________________
AliRsnEvent::~AliRsnEvent()
{
//
// Destructor (does nothing since there are not owned pointers)
//
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughter(AliRsnDaughter &out, Int_t index, Bool_t fromMC)
{
//
// Assigns to the first argument the reference to the i-th track in the ref event.
// What assignment method to be used will depend on the index and on the type of input.
// If the last argument is kTRUE and an MC is referenced, then both fRef and fRefMC will
// point to the MC particle (pure MC analysis)
//

   // by default the daughter is reset
   // and assigned the used index
   out.Reset();
   out.SetRsnID(index);
   out.SetOwnerEvent(this);
   
   // check input type
   if (!InputOK()) return;
   Bool_t inputESD = IsESD();
   

      Int_t trueIndex;
      AliRsnDaughter::ERefType type;
      if (!ConvertAbsoluteIndex(index, trueIndex, type)) {
         AliError(Form("Failed to convert absolute index %d", index));
         return;
      }
      switch (type) {
         case AliRsnDaughter::kTrack:
            if (inputESD) {
	        SetDaughterESDtrack(out, trueIndex); 
		
	    } else {
		SetDaughterAODtrack(out, trueIndex);
	    }
            break;
         case AliRsnDaughter::kV0:
            if (inputESD) SetDaughterESDv0(out, trueIndex); else SetDaughterAODv0(out, trueIndex);
            break;
         case AliRsnDaughter::kCascade:
            if (inputESD) SetDaughterESDcascade(out, trueIndex); else SetDaughterAODcascade(out, trueIndex);
            break;
         default:
            AliError("Unrecognized daughter type");
            return;
      }

   // if it is pure MC, the index tells what particle
   // to be read in the stack of MC particles, otherwise
   // it is converted into a real collection index
   if (fromMC) out.SetRef(out.GetRefMC());
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughter(Int_t i, Bool_t fromMC)
{
//
// Returns a daughter set using same criteria as SetDaughter
//

   AliRsnDaughter d;
   SetDaughter(d, i, fromMC);
   return d;
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterESDtrack(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #1: ESD tracks
//

   AliESDEvent *esd = (AliESDEvent*)fRef;
   
   if (i >= 0 && i < esd->GetNumberOfTracks()) {
      AliESDtrack *track = esd->GetTrack(i);
      if (track) {
         out.SetRef(track);
         out.SetGood();
         // if MC is present, assign label and retrieve corresponding particle
         if (fRefMC) {
            out.SetLabel(TMath::Abs(track->GetLabel()));
            if (!SetMCInfoESD(out)) {
               AliWarning("Failed assignment of MC info");
            }
         }
      } else {
         AliWarning("Null track");
      }
   } else {
      AliWarning(Form("Overflow: required index = %d, max = %d", i, esd->GetNumberOfTracks()));
   }
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterAODtrack(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #2: AOD tracks
//

   AliAODEvent *aod = (AliAODEvent*)fRef;

   if (i >= 0 && i < aod->GetNumberOfTracks()) {
      AliAODTrack *track = aod->GetTrack(i);
      if (track) {
         out.SetRef(track);
         out.SetGood();
         // if MC is present, assign label and retrieve corresponding particle
         if (fRefMC) {
            out.SetLabel(TMath::Abs(track->GetLabel()));
            if (!SetMCInfoAOD(out)) {
               AliWarning("Failed assignment of MC info");
            }
         }
      } else {
         AliWarning("Null track");
      }
   } else {
      AliWarning(Form("Overflow: required index = %d, max = %d", i, aod->GetNumberOfTracks()));
   }
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterESDv0(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #3: ESD v0
//

   if (i >= 0 && i < fRef->GetNumberOfV0s()) {
      AliESDEvent *ev = GetRefESD();
      AliESDv0    *v0 = ev->GetV0(i);
      if (v0) {
         out.SetRef(v0);
         out.SetGood();
         // if MC is present, retrieve the label of V0 from those of daughters
         if (fRefMC) {
            AliMCEvent  *mc = (AliMCEvent*)fRefMC;
            AliESDtrack *tp = ev->GetTrack(v0->GetPindex());
            AliESDtrack *tn = ev->GetTrack(v0->GetNindex());
            if (mc && tp && tn) {
               Int_t lp = TMath::Abs(tp->GetLabel());
               Int_t ln = TMath::Abs(tn->GetLabel());
               TParticle *pp = mc->Stack()->Particle(lp);
               TParticle *pn = mc->Stack()->Particle(ln);
               if (pp && pn) {
                  // if their first mothers are the same, the V0 is true
                  // otherwise label remains '-1' --> fake V0
                  if (pp->GetFirstMother() == pn->GetFirstMother() && pp->GetFirstMother() >= 0) {
                     out.SetLabel(pp->GetFirstMother());
                     SetMCInfoESD(out);
                  }
               }
            }
         }
      }
   }
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterAODv0(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #4: AOD v0
//

   if (i >= 0 && i < fRef->GetNumberOfV0s()) {
      AliAODEvent *ev = (AliAODEvent*)fRef;
      AliAODv0    *v0 = ev->GetV0(i);
      if (v0) {
         out.SetRef(v0);
         out.SetGood();
         if (fRefMC) {
            AliAODEvent  *mc = (AliAODEvent*)fRefMC;
            TClonesArray *mcArray = (TClonesArray*)mc->GetList()->FindObject(AliAODMCParticle::StdBranchName());
            AliAODTrack  *tp  = (AliAODTrack*)v0->GetDaughter(0);
            AliAODTrack  *tn  = (AliAODTrack*)v0->GetDaughter(1);
            if (mcArray && tp && tn) {
               Int_t lp = TMath::Abs(tp->GetLabel());
               Int_t ln = TMath::Abs(tn->GetLabel());
               AliAODMCParticle *pp = (AliAODMCParticle*)mcArray->At(lp);
               AliAODMCParticle *pn = (AliAODMCParticle*)mcArray->At(ln);
               if (pp && pn) {
                  // if their first mothers are the same, the V0 is true
                  // otherwise label remains '-1' --> fake V0
                  if (pp->GetMother() == pn->GetMother() && pp->GetMother() >= 0) {
                     out.SetLabel(pp->GetMother());
                     SetMCInfoAOD(out);
                  }
               }
            }
         }
      }
   }
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterESDcascade(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #3: ESD cascade
//

   if (i >= 0 && i < fRef->GetNumberOfCascades()) {
      AliESDEvent   *ev   = GetRefESD();
      AliESDcascade *casc = ev->GetCascade(i);
      if (casc) {
         out.SetRef(casc);
         out.SetGood();
         if (fRefMC) {
         
         }
      }
   }
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughterAODcascade(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #4: AOD cascade
//

   if (i >= 0 && i < fRef->GetNumberOfCascades()) {
      AliAODEvent *ev = GetRefAOD();
      AliAODv0    *casc = ev->GetCascade(i);
      if (casc) {
         out.SetRef(casc);
         out.SetGood();
         if (fRefMC) {
            
         }
      }
   }
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetMCInfoESD(AliRsnDaughter &out)
{
//
// Using the label assigned to the daughter, searches for the MC informations:
// -- MC reference
// -- mother
//

   // if label makes no sense --> failed
   Int_t label = out.GetLabel();
   if (label < 0 || !fRefMC) return kFALSE;
   
   // get number of particles
   Int_t nMC = fRefMC->GetNumberOfTracks();
   
   // if label too large --> failed
   if (label >= nMC) {
      AliWarning(Form("Stack overflow: track label = %d -- stack maximum = %d", label, nMC));
      return kFALSE;
   }

   // retrieve particle
   AliMCEvent    *mc = (AliMCEvent*)fRefMC;
   AliMCParticle *mcPart = (AliMCParticle*)mc->GetTrack(label);
   
   // if particle = NULL --> failed
   if (!mcPart) {
      AliWarning(Form("Stack discontinuity: label %d refers to a NULL object", label));
      return kFALSE;
   }
   // otherwise --> success
   out.SetRefMC(mcPart);

   // if the particle is not primary, find the mother and get its PDG
   Int_t imum = mcPart->Particle()->GetFirstMother();
   if (imum >= 0 && imum < nMC) {
      AliMCParticle *mcMother = (AliMCParticle*)mc->GetTrack(imum);
      if (mcMother) {
         out.SetMotherPDG(TMath::Abs(mcMother->Particle()->GetPdgCode()));
      } else {
         AliWarning(Form("Stack discontinuity: label mother %d refers to a NULL object", imum));
      }
   } else {
      AliWarning(Form("Stack overflow: mother label = %d -- stack maximum = %d", imum, nMC));
   }

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetMCInfoAOD(AliRsnDaughter &out)
{
//
// Using the label assigned to the daughter, searches for the MC informations:
// -- MC reference
// -- mother
//

   // if label makes no sense --> failed
   Int_t label = out.GetLabel();
   if (label < 0 || !fRefMC) return kFALSE;
   
   // retrieve particle
   AliAODEvent  *mc = (AliAODEvent*)fRefMC;
   TClonesArray *mcArray = (TClonesArray*)mc->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   
   // get number of particles
   Int_t nMC = mcArray->GetEntriesFast();
   
   // if label too large --> failed
   if (label >= nMC) {
      AliWarning(Form("Stack overflow: track label = %d -- stack maximum = %d", label, nMC));
      return kFALSE;
   }
   
   // if particle = NULL --> failed
   AliAODMCParticle *mcPart = (AliAODMCParticle*)mcArray->At(label);
   if (!mcPart) {
      AliWarning(Form("Stack discontinuity: label %d refers to a NULL object", label));
      return kFALSE;
   }
   // otherwise --> success
   out.SetRefMC(mcPart);

   // if the particle is not primary, find the mother and get its PDG
   Int_t imum = mcPart->GetMother();
   if (imum >= 0 && imum < nMC) {
      AliAODMCParticle *mcMother = (AliAODMCParticle*)mcArray->At(imum);
      if (mcMother) {
         out.SetMotherPDG(TMath::Abs(mcMother->GetPdgCode()));
      } else {
         AliWarning(Form("Stack discontinuity: label mother %d refers to a NULL object", imum));
      }
   } else if (imum >= nMC) {
      AliWarning(Form("Stack overflow: mother label = %d -- stack maximum = %d", imum, nMC));
   }

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::ConvertAbsoluteIndex(Int_t index, Int_t &realIndex, AliRsnDaughter::ERefType &type)
{
//
// Using the phylosophy of the absolute index, which loops over
// all tracks, V0s and cascades, returns the result of a check
// on it (first argument) based on this criterion:
// 1) if the absolute ID is smaller than number of tracks,
//    return itself and the type 'track'
// 2) if the absolute ID is larger than number of tracks, subtract it
//    and if the result is smaller than number of V0s,
//    return the corresponding V0 index and type
// 3) if the absolute ID is larger than number of tracks + V0s, subtract them
//    and if the result is smaller than number of cascades,
//    return the corresponding cascade index and type
// The results of this check are stored in the reference arguments, while the outcome of
// the function is kTRUE if one of these checks was successful, otherwise it is kFALSE,
// meaning that the absolute index reached the end.
//

   Int_t nTracks   = fRef->GetNumberOfTracks();
   Int_t nV0s      = fRef->GetNumberOfV0s();
   Int_t nCascades = fRef->GetNumberOfCascades();

   if (index < nTracks) {
      realIndex = index;
      type = AliRsnDaughter::kTrack;
      return kTRUE;
   } else if (index >= nTracks && index < nTracks + nV0s) {
      realIndex = index - nTracks;
      type = AliRsnDaughter::kV0;
      return kTRUE;
   } else if (index >= nTracks + nV0s && index < nTracks + nV0s + nCascades) {
      realIndex = index - nTracks - nV0s;
      type = AliRsnDaughter::kCascade;
      return kTRUE;
   }

   realIndex = -1;
   type = AliRsnDaughter::kNoType;
   return kFALSE;
}

//_____________________________________________________________________________
Int_t AliRsnEvent::ConvertRealIndex(Int_t index, AliRsnDaughter::ERefType type)
{
//
// Translates a pair made by index + object type into the corresponding
// absolute index, which is set to -1 in case the real index overflows.
//

   Int_t nTracks   = fRef->GetNumberOfTracks();
   Int_t nV0s      = fRef->GetNumberOfV0s();
   Int_t nCascades = fRef->GetNumberOfCascades();

   switch (type) {
      case AliRsnDaughter::kTrack:
         if (index >= 0 && index < nTracks)
            return index;
         else
            return -1;
      case AliRsnDaughter::kV0:
         if (index >= 0 && index < nV0s)
            return nTracks + index;
         else
            return -1;
      case AliRsnDaughter::kCascade:
         if (index >= 0 && index < nCascades)
            return nTracks + nV0s + index;
         else
            return -1;
      default:
         return -1;
   }
}

//_____________________________________________________________________________
Int_t AliRsnEvent::SelectLeadingParticle(AliRsnCutSet *cuts)
{
//
// Searches the collection of all particles with given PID type and charge,
// and returns the one with largest momentum, provided that it is greater than 1st argument.
// If one specifies AliRsnPID::kUnknown as type or AliRsnDaughter::kNoPID as method,
// the check is done over all particles irrespectively of their PID.
// If the sign argument is '+' or '-', the check is done over the particles of that charge,
// otherwise it is done irrespectively of the charge.
//

   // check input type
   Bool_t inputESD = IsESD();
   if (!inputESD && !IsAOD()) {
      AliError("Need to process ESD or AOD input");
      return -1;
   }

   Double_t ptMax = 0.0;
   Int_t i, nTracks = fRef->GetNumberOfTracks();
   
   fLeading = -1;
   AliRsnDaughter leading;

   for (i = 0; i < nTracks; i++) {
      if (inputESD)
         SetDaughterESDtrack(leading, i);
      else
         SetDaughterAODtrack(leading, i);
      if (!leading.IsOK()) {
         AliDebugClass(1, Form("Failed assignment of track %d", i));
         continue;
      }
      if (cuts && !cuts->IsSelected(&leading)) {
         AliDebugClass(1, Form("Track %d didn't pass cuts", i));
         continue;
      }
      // check if it has largest momentum
      if (leading.GetRef()->Pt() > ptMax) {
         ptMax = leading.GetRef()->Pt();
         fLeading = i;
      }
   }
   
   return fLeading;
}
