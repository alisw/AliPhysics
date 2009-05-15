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

//
// *** Class AliRsnEvent ***
//
// A container for a collection of AliRsnDaughter objects from an event.
// Contains also the primary vertex, useful for some cuts.
// In order to retrieve easily the tracks which have been identified
// as a specific type and charge, there is an array of indexes which
// allows to avoid to loop on all tracks and have only the neede ones.
//
// authors: A. Pulvirenti (email: alberto.pulvirenti@ct.infn.it)
//          M. Vala (email: martin.vala@cern.ch)
//

#include <Riostream.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(AliVEvent *ref, AliMCEvent *refMC) :
    fRef(ref),
    fRefMC(refMC)
{
//
// Default constructor
//
}

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(const AliRsnEvent &event) :
    TObject(event),
    fRef(event.fRef),
    fRefMC(event.fRefMC)

{
//
// Copy constructor.
//
}

//_____________________________________________________________________________
AliRsnEvent& AliRsnEvent::operator= (const AliRsnEvent &event)
{
//
// Works in the same way as the copy constructor.
//

  (TObject)(*this) = (TObject)event;
  fRef = event.fRef;
  fRefMC = event.fRefMC;

  return (*this);
}

//_____________________________________________________________________________
AliRsnEvent::~AliRsnEvent()
{
//
// Destructor.
//
}

//_____________________________________________________________________________
void AliRsnEvent::SetDaughter(AliRsnDaughter &out, Int_t i)
{
//
// Return a track stored here in format of AliRsnDaughter.
// and finds in the reference event the informations to set
// the proprietary data members of AliRsnDaughter
//

  // retrieve reference particle from reference event
  AliVParticle *ref = (AliVParticle*)fRef->GetTrack(i);

  if (!ref) return;

  // if MC info is present, retrieve from it
  TParticle *refMC = 0;
  if (fRefMC) {
    Int_t label = TMath::Abs(ref->GetLabel());
    refMC = fRefMC->Stack()->Particle(label);
  }

  // create output object
  out.SetRef(ref);
  out.SetGood();
  out.SetParticle(refMC);
  if (fRefMC)
    out.FindMotherPDG(fRefMC->Stack());

  // retrieve primary vertex and set impact parameters
  Double_t dx = out.Xv(), dy = out.Yv(), dz = out.Zv();
  const AliVVertex *v = fRef->GetPrimaryVertex();
  if (v) {
    dx -= v->GetX();
    dy -= v->GetY();
    dz -= v->GetZ();
  }
  out.SetDr(TMath::Sqrt(dx*dx + dy*dy));
  out.SetDz(dz);

  // dynamic reference to true nature of referenced event
  // to get kink index
  AliESDEvent *esd = dynamic_cast<AliESDEvent*>(fRef);
  AliAODEvent *aod = dynamic_cast<AliAODEvent*>(fRef);

  if (esd) {
    AliESDtrack *esdTrack = esd->GetTrack(i);
    out.FindKinkIndex(esdTrack);
  } else if (aod) {
    out.FindKinkIndex(aod);
  }

  out.SetGood();
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughter(Int_t i)
{
//
// Return an AliRsnDaughter taken from this event,
// with all additional data members well set.
//

  AliRsnDaughter out;
  SetDaughter(out, i);

  return out;
}

//_____________________________________________________________________________
Int_t AliRsnEvent::GetMultiplicity()
{
//
// Returns event multiplicity
//
  AliDebug(AliLog::kDebug+2,"<-");
  if (!fRef) return 0;
  AliDebug(AliLog::kDebug+2,"->");
  return fRef->GetNumberOfTracks();
}

//_____________________________________________________________________________
Double_t AliRsnEvent::GetVz()
{
//
// Return Z coord of primary vertex
//
  AliDebug(AliLog::kDebug+2,"<-");
  return fRef->GetPrimaryVertex()->GetZ();
  AliDebug(AliLog::kDebug+2,"->");
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetLeadingParticle
(Double_t ptMin, AliPID::EParticleType type)
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
  AliRsnDaughter output;

  for (i = 0; i < nTracks; i++) {
    AliRsnDaughter track = GetDaughter(i);
    if (!AcceptTrackPID(&track, type)) continue;
    if (track.Pt() < ptMin) continue;
    if (!output.IsOK() || track.Pt() > output.Pt()) {
      output = track;
      output.SetGood();
    }
  }

  return output;
}

//_________________________________________________________________________________________________
Double_t AliRsnEvent::GetAverageMomentum(Int_t &count, AliPID::EParticleType type)
{
//
// Loops on the list of tracks and computes average total momentum.
//

  Int_t i, nTracks = fRef->GetNumberOfTracks();
  Double_t pmean = 0.0;

  for (i = 0, count = 0; i < nTracks; i++) {
    AliRsnDaughter track = GetDaughter(i);
    if (!AcceptTrackPID(&track, type)) continue;
    pmean += track.P();
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

    angle = leading.AngleTo(trk);

    angleMean += angle;
    angle2Mean += angle * angle;
    count++;
  }

  if (!count) return kFALSE;

  angleMean /= (Double_t)count;
  angle2Mean /= (Double_t)count;
  angleRMS = TMath::Sqrt(angle2Mean - angleMean*angleMean);

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::AcceptTrackPID
(AliRsnDaughter *d, AliPID::EParticleType type)
{
//
// [PRIVATE]
// Checks if the track PID (according to method in use) corresponds
// to the required identification species.
// If the second argument is "kUnknown", answer of this method is always YES.
//

  if (type == AliPID::kUnknown) return kTRUE;

  return (d->AssignedPID() == type);
}
