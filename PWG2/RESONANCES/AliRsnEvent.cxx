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

#include <TArrayF.h>

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliRsnCutPID.h"

#include "AliRsnEvent.h"

ClassImp(AliRsnEvent)

//_____________________________________________________________________________
AliRsnEvent::AliRsnEvent(AliVEvent *ref, AliMCEvent *refMC) :
  fRef(ref),
  fRefMC(refMC),
  fLeading(-1)
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
  fLeading(event.fLeading)
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
Bool_t AliRsnEvent::SetDaughter(AliRsnDaughter &out, Int_t i, AliRsnDaughter::ERefType type)
{
//
// Using the second and third arguments, retrieves the i-th object in the
// appropriate sample (tracks or V0s) and sets the firs reference object
// in order to point to that.
// If a MonteCarlo information is provided, sets the useful informations from there,
// and in case of a V0, sets the 'label' data member only when the two daughters
// of the V0 point to the same mother.
// Returns kFALSE whenever the operation fails (out of range, NULL references).
//

  Int_t label;

  // retrieve reference particle from reference event
  // if it is found, by defaul track can be used (good)
  if (type == AliRsnDaughter::kTrack)
  {
    if (i >= fRef->GetNumberOfTracks())
    {
      out.SetBad();
      return kFALSE;
    }
    AliVTrack *track = (AliVTrack*)fRef->GetTrack(i);
    if (!track)
    {
      out.SetBad();
      return kFALSE;
    }
    else
    {
      label = TMath::Abs(track->GetLabel());
      out.SetRef(track);
      out.SetLabel(label);
      if (fRefMC)
      {
        if (label < fRefMC->GetNumberOfTracks()) 
        {
          AliMCParticle *part = (AliMCParticle*)fRefMC->GetTrack(label);
          out.SetRefMC(part);
        }
      }
      out.SetGood();
    }
  }
  else if (type == AliRsnDaughter::kV0)
  {
    if (i > fRef->GetNumberOfV0s())
    {
      out.SetBad();
      return kFALSE;
    }
    AliESDv0     *esdV = 0x0;
    AliAODv0     *aodV = 0x0;
    Int_t         lp, ln;
    AliVTrack    *tp = 0x0, *tn = 0x0;
    if (IsESD()) esdV = GetRefESD()->GetV0(i);
    if (IsAOD()) aodV = GetRefAOD()->GetV0(i);
    if (!esdV && !aodV)
    {
      out.SetBad();
      return kFALSE;
    }
    else
    {
      if (esdV) out.SetRef(esdV); else out.SetRef(aodV);
      // retrieve the V0 daughters, which must be done differently with ESD and AOD v0s
      if (esdV)
      {
        // get the 2 daughters of the V0
        AliESDEvent *ev = dynamic_cast<AliESDEvent*>(fRef);
        tp = ev->GetTrack(esdV->GetPindex());
        tn = ev->GetTrack(esdV->GetNindex());
      }
      else if (aodV)
      {
        // get the 2 daughters of the V0
        AliAODEvent *ev = dynamic_cast<AliAODEvent*>(fRef);
        tp = ev->GetTrack(aodV->GetPosID());
        tn = ev->GetTrack(aodV->GetNegID());
      }

      // now, if we have a MC, use the two track objects
      // to retrieve the true particle which generated the V0
      // using their labels; by default they are a false V0 with label -1
      label = -1;
      if (tp && tn && fRefMC)
      {
        lp = TMath::Abs(tp->GetLabel());
        ln = TMath::Abs(tn->GetLabel());
        // if labels are meaningful, retrieve corresponding particles
        TParticle *pp = fRefMC->Stack()->Particle(lp);
        TParticle *pn = fRefMC->Stack()->Particle(ln);
        // if their first mothers are the same, the V0 is true
        // otherwise no label can be assigned
        if (pp->GetFirstMother() == pn->GetFirstMother()) label = pp->GetFirstMother();
      }
      out.SetLabel(label);
      out.SetGood();
    }
  }
  
  // finally, in case we have a MC, searches for the mother, in order to store
  // its PDG code into the output AliRsnDaughter
  if (fRefMC)
  {
    label = out.GetLabel();
    AliStack *stack = fRefMC->Stack();
    if (label >= 0 && label < stack->GetNtrack())
    {
      TParticle *part = stack->Particle(label);
      if (part)
      {
        Int_t imum = part->GetFirstMother();
        if (imum >= 0 && imum <= stack->GetNtrack())
        {
          TParticle *mum = stack->Particle(imum);
          if (mum) out.SetMotherPDG(TMath::Abs(mum->GetPdgCode()));
        }
      }
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnEvent::SetDaughterMC(AliRsnDaughter &out, Int_t i)
{
//
// Using the second argument, retrieves the i-th object in the
// appropriate sample and sets the firs reference object
// in order to point to that.
// If a MonteCarlo information is provided, sets the useful informations from there,
// and in case of a V0, sets the 'label' data member only when the two daughters
// of the V0 point to the same mother.
// Returns kFALSE whenever the operation fails (out of range, NULL references).
//

  if (!fRefMC)
  {
    out.SetBad();
    return kFALSE;
  }

  if (i >= fRefMC->GetNumberOfTracks())
  {
    out.SetBad();
    return kFALSE;
  }
  
  AliMCParticle *track = (AliMCParticle*)fRef->GetTrack(i);
  if (!track)
  {
    out.SetBad();
    return kFALSE;
  }
  else
  {
    out.SetRef(track);
    out.SetRefMC(track);
    out.SetLabel(i);
    out.SetGood();
  }
  
  AliStack  *stack = fRefMC->Stack();
  TParticle *part  = track->Particle();
  if (part)
  {
    Int_t imum = part->GetFirstMother();
    if (imum >= 0 && imum <= stack->GetNtrack())
    {
      TParticle *mum = stack->Particle(imum);
      if (mum) out.SetMotherPDG(TMath::Abs(mum->GetPdgCode()));
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughter(Int_t i, AliRsnDaughter::ERefType type)
{
//
// Return an AliRsnDaughter taken from this event,
// with all additional data members well set.
//

  AliRsnDaughter out;
  SetDaughter(out, i, type);

  return out;
}

//_____________________________________________________________________________
AliRsnDaughter AliRsnEvent::GetDaughterMC(Int_t i)
{
//
// Return an AliRsnDaughter taken from this event,
// with all additional data members well set.
//

  AliRsnDaughter out;
  SetDaughterMC(out, i);

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
    AliVParticle *ref = track.GetRef();
    if (ref->Pt() < ptMin) continue;
    //double pt = track.P().Perp();
    //Printf("track %d %g", i, pt);
    if (!leading.IsOK() || ref->Pt() > ptMin)
    {
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
    pmean += track.P().Mag();
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

    angle = leading.P().Angle(trk.P().Vect());

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
