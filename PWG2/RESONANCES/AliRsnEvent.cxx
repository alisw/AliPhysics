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

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliAODEvent.h"
#include "AliRsnCutPID.h"
#include "AliESDtrackCuts.h"

#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

AliRsnEvent* AliRsnEvent::fgRsnEvent1 = 0;
AliRsnEvent* AliRsnEvent::fgRsnEvent2 = 0;

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(AliVEvent *ref, AliVEvent *refMC) :
   fRef(ref),
   fRefMC(refMC),
   fLeading(-1),
   fLocalID(-1)
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
   fLocalID(event.fLocalID)
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

   (TObject)(*this) = (TObject)event;
   fRef             = event.fRef;
   fRefMC           = event.fRefMC;
   fLeading         = event.fLeading;
   fLocalID         = event.fLocalID;

   return (*this);
}

//_____________________________________________________________________________
AliRsnEvent::~AliRsnEvent()
{
//
// Destructor.
// Dereferences global pointer, if needed.
//

   //if (gRsnCurrentEvent == this) gRsnCurrentEvent = 0;
   //if (gRsnMixedEvent   == this) gRsnMixedEvent = 0;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughter(AliRsnDaughter &out, Int_t i, AliRsnDaughter::ERefType type)
{
//
// Using the second and third arguments, retrieves the i-th object in the
// appropriate sample (tracks or V0s) and sets the first reference object
// in order to point to that.
// If a MonteCarlo information is provided, sets the useful informations from there,
// and in case of a V0, sets the 'label' data member only when the two daughters
// of the V0 point to the same mother.
// Returns kFALSE whenever the operation fails (out of range, NULL references).
//

   Bool_t ok = kFALSE;

   if (IsESD() && type == AliRsnDaughter::kTrack)   ok = SetDaughterESDtrack  (out, i);
   if (IsAOD() && type == AliRsnDaughter::kTrack)   ok = SetDaughterAODtrack  (out, i);
   if (IsESD() && type == AliRsnDaughter::kV0)      ok = SetDaughterESDv0     (out, i);
   if (IsAOD() && type == AliRsnDaughter::kV0)      ok = SetDaughterAODv0     (out, i);
   if (IsESD() && type == AliRsnDaughter::kCascade) ok = SetDaughterESDcascade(out, i);
   if (IsAOD() && type == AliRsnDaughter::kCascade) ok = SetDaughterAODcascade(out, i);

   return ok;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterAbs(AliRsnDaughter &out, Int_t absIndex)
{
//
// Sets the first argument daughter using the absolute index, which
// runs continuously from tracks, to V0s, to cascades.
// In case the conversion to real index fails, the track is flagged as bad.
// Additionally, sets the daughter internal 'fRsnID' member to this index.
//

   Int_t index;
   AliRsnDaughter::ERefType type;

   if (!ConvertAbsoluteIndex(absIndex, index, type)) {
      out.Reset();
      out.SetBad();
      return kFALSE;
   }
   
   SetDaughter(out, index, type);
   out.SetRsnID(absIndex);
   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterMC(AliRsnDaughter &out, Int_t label)
{
//
// Using the second argument, retrieves the i-th object in the
// MC sample (if present) and assigns the track using only that,
// so that it is considered both as main reference and MC reference.
// (used for MC-only analysis).
//

   if (!fRefMC) {
      out.SetBad();
      return kFALSE;
   }

   // try to convert into both types
   Int_t        imum;
   AliMCEvent  *esd = GetRefMCESD();
   AliAODEvent *aod = GetRefMCAOD();

   // ESD
   if (esd) {
      // if the MC track exists, assign it
      AliMCParticle *track = (AliMCParticle*)fRef->GetTrack(label);
      if (!track) {
         out.SetBad();
         return kFALSE;
      }
      out.SetRef(track);
      out.SetRefMC(track);
      out.SetLabel(label);
      out.SetGood();

      // search for its mother in stack
      imum = track->GetMother();
      if (imum >= 0 && imum < esd->Stack()->GetNtrack()) {
         TParticle *mum = esd->Stack()->Particle(imum);
         if (mum) out.SetMotherPDG(TMath::Abs(mum->GetPdgCode()));
      }
   }

   // AOD
   if (aod) {
      // checks that array of MC particles exists
      TClonesArray *mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
      if (!mcArray) {
         out.SetBad();
         return kFALSE;
      }

      // in this case one has to loop over the sample to find the good one
      TObjArrayIter next(mcArray);
      AliAODMCParticle *part = 0x0;
      while ((part = (AliAODMCParticle*)next())) {
         if (TMath::Abs(part->GetLabel()) == label) {
            // if the MC track exists, assign it
            out.SetRef(part);
            out.SetRefMC(part);
            out.SetLabel(label);
            out.SetGood();

            // search for the mother
            imum = part->GetMother();
            if (imum >= 0 && imum < mcArray->GetEntriesFast()) {
               AliAODMCParticle *mum = (AliAODMCParticle*)mcArray->At(imum);
               if (mum) out.SetMotherPDG(TMath::Abs(mum->GetPdgCode()));
            }
            break;
         }
      }
      return kTRUE;
   }

   return kFALSE;
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughter(Int_t i, AliRsnDaughter::ERefType type)
{
//
// Returns a daughter set using same criteria as SetDaughter
//

   AliRsnDaughter d;
   SetDaughter(d, i, type);
   return d;
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughterAbs(Int_t absIndex)
{
//
// Returns a daughter set using same criteria as SetDaughter
//

   AliRsnDaughter d;
   SetDaughterAbs(d, absIndex);
   return d;
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughterMC(Int_t i)
{
//
// Returns a daughter set using same criteria as SetDaughterMC
//

   AliRsnDaughter d;
   SetDaughterMC(d, i);
   return d;
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetAbsoluteSum()
{
//
// Utility method that returns the sum of all daughter-like objects:
// tracks, V0s and cascades
//

   Int_t total = 0;

   total += fRef->GetNumberOfTracks();
   total += fRef->GetNumberOfV0s();
   total += fRef->GetNumberOfCascades();

   return total;
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
Int_t AliRsnEvent::GetMultiplicity(AliESDtrackCuts *cuts)
{
//
// Returns event multiplicity as the number of tracks.
// If the argument is not NULL, returns instead the
// number of tracks passing the cuts hereby defined.
//

   if (!fRef) return 0;

   AliESDEvent *esd = GetRefESD();
   if (cuts && esd) return cuts->CountAcceptedTracks(esd);
   else return fRef->GetNumberOfTracks();
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetMultiplicityMC()
{
//
// Returns event multiplicity as the number of primary tracks
// counted from the MC sample, if present.
//

   if (!fRefMC) return 0;

   if (IsESD())
      return GetRefMCESD()->Stack()->GetNprimary();
   else if (IsAOD())
      return GetRefMCAOD()->GetNumberOfTracks();
   else
      return 0;
}

//_____________________________________________________________________________
Double_t AliRsnEvent::GetVz()
{
//
// Return Z coord of primary vertex
//

   return fRef->GetPrimaryVertex()->GetZ();
}

//_____________________________________________________________________________
Int_t AliRsnEvent::SelectLeadingParticle
(Double_t ptMin, AliRsnCutPID *cutPID)
{
//
// Searches the collection of all particles with given PID type and charge,
// and returns the one with largest momentum, provided that it is greater than 1st argument.
// If one specifies AliRsnPID::kUnknown as type or AliRsnDaughter::kNoPID as method,
// the check is done over all particles irrespectively of their PID.
// If the sign argument is '+' or '-', the check is done over the particles of that charge,
// otherwise it is done irrespectively of the charge.
//

   Int_t i, nTracks = fRef->GetNumberOfTracks();
   fLeading = -1;
   AliRsnDaughter leading;
   leading.SetBad();

   for (i = 0; i < nTracks; i++) {
      AliRsnDaughter track = GetDaughter(i);
      if (cutPID) if (!cutPID->IsSelected(&track)) continue;
      const AliVParticle *ref = track.GetRef();
      if (ref->Pt() < ptMin) continue;
      //double pt = track.P().Perp();
      //Printf("track %d %g", i, pt);
      if (!leading.IsOK() || ref->Pt() > ptMin) {
         fLeading = i;
         //leading = track;
         ptMin = ref->Pt();
      }
   }
   return fLeading;
}

//_________________________________________________________________________________________________
Double_t AliRsnEvent::GetAverageMomentum(Int_t &count, AliRsnCutPID *cutPID)
{
//
// Loops on the list of tracks and computes average total momentum.
//

   Int_t i, nTracks = fRef->GetNumberOfTracks();
   Double_t pmean = 0.0;

   for (i = 0, count = 0; i < nTracks; i++) {
      AliRsnDaughter track = GetDaughter(i);
      if (cutPID) if (!cutPID->IsSelected(&track)) continue;
      pmean += track.Prec().Mag();
      count++;
   }

   if (count > 0) pmean /= (Double_t)count;
   else pmean = 0.0;

   return pmean;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::GetAngleDistr
(Double_t &angleMean, Double_t &angleRMS, AliRsnDaughter leading)
{
//
// Takes the leading particle and computes the mean and RMS
// of the distribution of directions of all other tracks
// with respect to the direction of leading particle.
//

   if (!leading.IsOK()) return kFALSE;

   Int_t i, count, nTracks = fRef->GetNumberOfTracks();
   Double_t angle, angle2Mean = 0.0;

   angleMean = angle2Mean = 0.0;

   for (i = 0, count = 0; i < nTracks; i++) {
      AliRsnDaughter trk = GetDaughter(i);
      if (trk.GetID() == leading.GetID()) continue;

      angle = leading.Prec().Angle(trk.Prec().Vect());

      angleMean += angle;
      angle2Mean += angle * angle;
      count++;
   }

   if (!count) return kFALSE;

   angleMean /= (Double_t)count;
   angle2Mean /= (Double_t)count;
   angleRMS = TMath::Sqrt(angle2Mean - angleMean * angleMean);

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterESDtrack(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #1: ESD tracks
//

   // check 1: index in good range
   if (i < 0 || i >= fRef->GetNumberOfTracks()) {
      out.SetBad();
      return kFALSE;
   }

   // check 2: not NULL object
   AliVTrack *track = (AliVTrack*)fRef->GetTrack(i);
   if (!track) {
      out.SetBad();
      return kFALSE;
   }

   // assign references of reconstructed track
   Int_t label = TMath::Abs(track->GetLabel());
   out.SetRef(track);
   out.SetLabel(label);
   out.SetGood();

   // assign MC info, if available
   return SetMCInfoESD(out);
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterAODtrack(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #2: AOD tracks
//

   // check 1: index in good range
   if (i < 0 || i >= fRef->GetNumberOfTracks()) {
      out.SetBad();
      return kFALSE;
   }

   // check 2: not NULL object
   AliVTrack *track = (AliVTrack*)fRef->GetTrack(i);
   if (!track) {
      out.SetBad();
      return kFALSE;
   }

   // assign references of reconstructed track
   Int_t label = TMath::Abs(track->GetLabel());
   out.SetRef(track);
   out.SetLabel(label);
   out.SetGood();

   // assign MC info, if available
   return SetMCInfoAOD(out);
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterESDv0(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #3: ESD v0
//

   // check 1: index in good range
   if (i > fRef->GetNumberOfV0s()) {
      out.SetBad();
      return 1;
   }

   // check 2: not NULL object
   AliESDEvent *ev = GetRefESD();
   AliESDv0    *v0 = ev->GetV0(i);
   if (!v0) {
      out.SetBad();
      return 2;
   }

   // assign references of reconstructed track
   out.SetRef(v0);
   out.SetGood();
   out.SetLabel(-1);

   // this time, assigning label is not trivial,
   // it is done only if MC is present and both
   // daughters come from a true particle
   AliMCEvent  *mc = GetRefMCESD();
   AliESDtrack *tp = ev->GetTrack(v0->GetPindex());
   AliESDtrack *tn = ev->GetTrack(v0->GetNindex());
   if (mc && tp && tn) {
      Int_t        lp = TMath::Abs(tp->GetLabel());
      Int_t        ln = TMath::Abs(tn->GetLabel());
      TParticle   *pp = mc->Stack()->Particle(lp);
      TParticle   *pn = mc->Stack()->Particle(ln);
      // if their first mothers are the same, the V0 is true
      // otherwise no label can be assigned
      if (pp->GetFirstMother() == pn->GetFirstMother() && pp->GetFirstMother() >= 0) out.SetLabel(pp->GetFirstMother());
   }

   // assign MC info, if available
   return SetMCInfoESD(out);
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterAODv0(AliRsnDaughter &out, Int_t i)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #4: AOD v0
//

   // check 1: index in good range
   if (i > fRef->GetNumberOfV0s()) {
      out.SetBad();
      return kFALSE;
   }

   // check 2: not NULL object
   AliAODEvent *ev = GetRefAOD();
   AliAODv0    *v0 = ev->GetV0(i);
   if (!v0) {
      out.SetBad();
      return kFALSE;
   }

   TClonesArray *mcArray = (TClonesArray*)ev->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!mcArray) {
      out.SetBad();
      return kFALSE;
   }

   // assign references of reconstructed track
   out.SetRef(v0);
   out.SetGood();
   out.SetLabel(-1);

   // this time, assigning label is not trivial,
   // it is done only if MC is present and both
   // daughters come from a true particle
   AliAODTrack  *tp = ev->GetTrack(v0->GetPosID());
   AliAODTrack  *tn = ev->GetTrack(v0->GetNegID());
   if (mcArray && tp && tn) {
      Int_t        lp = TMath::Abs(tp->GetLabel());
      Int_t        ln = TMath::Abs(tn->GetLabel());
      // loop on array to find MC daughters
      AliAODMCParticle *pp = 0x0, *pn = 0x0;
      TObjArrayIter next(mcArray);
      AliAODMCParticle *part = 0x0;
      while ((part = (AliAODMCParticle*)next())) {
         if (TMath::Abs(part->GetLabel()) == lp) pp = (AliAODMCParticle*)mcArray->IndexOf(part);
         if (TMath::Abs(part->GetLabel()) == ln) pn = (AliAODMCParticle*)mcArray->IndexOf(part);
      }
      // assign a MC reference and a label only to true V0s
      if (pp->GetMother() == pn->GetMother() && pp->GetMother() >= 0) out.SetLabel(pp->GetMother());
   }

   // assign MC info, if available
   return SetMCInfoAOD(out);
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterESDcascade(AliRsnDaughter &, Int_t)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #3: ESD cascade
//

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterAODcascade(AliRsnDaughter &, Int_t)
{
//
// Setup the first argument to the track identified by the index.
// When available, adds the MC information and references.
// ---
// Version #4: AOD cascade
//

   return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetMCInfoESD(AliRsnDaughter &out)
{
//
// Using the label assigned to the daughter, searches for the MC informations:
// -- MC reference
// -- mother
//

   Int_t label = out.GetLabel();

   // if no MC reference is available, exit here (successfully)
   AliMCEvent *mc = GetRefMCESD();
   if (!mc) return kTRUE;
   Int_t nMC = mc->GetNumberOfTracks();
   
   // debug message for fakes
   if (label < 0) {
      AliDebug(AliLog::kDebug + 1, "Fake object (fake track or false V0)");
      return kFALSE;
   }

   // assign MC reference, being aware of eventual
   // overflows in the array (sometimes happened)
   if (label >= nMC) {
      AliWarning(Form("Stack overflow: track label = %d -- stack maximum = %d", label, nMC));
      return kFALSE;
   }
   AliMCParticle *mcPart = (AliMCParticle*)mc->GetTrack(label);
   if (!mcPart) {
      AliWarning(Form("Stack discontinuity: label %d refers to a NULL object", label));
      return kFALSE;
   }
   out.SetRefMC(mcPart);

   // if this is a primary track, exit here (successfully)
   Int_t imum = mcPart->Particle()->GetFirstMother();
   if (imum < 0) {
      out.SetMotherPDG(0);
      return kTRUE;
   }

   // if didn't stop there, search for the mother
   if (imum >= nMC) {
      AliWarning(Form("Stack overflow: track mother label = %d -- stack maximum = %d", label, nMC));
      return kFALSE;
   }
   AliMCParticle *mcMother = (AliMCParticle*)mc->GetTrack(imum);
   if (!mcMother) {
      AliWarning(Form("Stack discontinuity: label mother %d refers to a NULL object", imum));
      return kFALSE;
   }
   out.SetMotherPDG(TMath::Abs(mcMother->Particle()->GetPdgCode()));

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

   Int_t label = out.GetLabel();
   
   // debug message for fakes
   if (label < 0) {
      AliDebug(AliLog::kDebug + 1, "Fake object (fake track or false V0)");
      return kFALSE;
   }

   // if no MC reference is available, exit here (successfully)
   AliAODEvent *mc = GetRefMCAOD();
   if (!mc) return kTRUE;

   // loop on particles using the track label as reference
   // and then assign also the mother type, if found
   TClonesArray *mcArray = (TClonesArray*)mc->GetList()->FindObject(AliAODMCParticle::StdBranchName());
   if (!mcArray) return kFALSE;
   TObjArrayIter next(mcArray);
   AliAODMCParticle *part = 0x0;
   while ((part = (AliAODMCParticle*)next())) {
      if (TMath::Abs(part->GetLabel()) == label) {
         out.SetRefMC(part);
         out.SetMotherPDG(0);
         Int_t imum = part->GetMother();
         if (imum >= 0 && imum <= mcArray->GetEntriesFast()) {
            AliAODMCParticle *mum = (AliAODMCParticle*)mcArray->At(imum);
            if (mum) out.SetMotherPDG(TMath::Abs(mum->GetPdgCode()));
            return kTRUE;
         } else {
            AliWarning(Form("Array overflow: track mother label = %d -- stack maximum = %d", imum, mcArray->GetEntriesFast()));
            return kFALSE;
         }

         // if a MC reference is found, there is no need to go on with the loop
         break;
      }
   }

   return kTRUE;
}
