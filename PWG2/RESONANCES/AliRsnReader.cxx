//
// Class AliRsnReader
//
// This is the universal converter from any kind of source event
// (i.e. ESD, standard AOD, MC) into the internal non-standard
// AOD format used by RSN package.
// ---
// This class reads all tracks in an input event and converts them
// into AliRsnDaughters, and computes preliminarily the PID probabilities
// by doing the Bayesian combination of a set of assigned prior probabilities
// with the PID weights defined in each track.
// ---
// When filling the output event (AliRsnEvent), some arrays of indexes
// are created in order to organize tracks according to their PID and charge,
// which will then be used in further analysis steps.
//
// author: A. Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include <TString.h>

#include "AliLog.h"

#include "AliVEvent.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODVertex.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"

#include "AliRsnMCInfo.h"
#include "AliRsnDaughter.h"
#include "AliRsnEvent.h"
#include "AliRsnPIDDefESD.h"

#include "AliRsnReader.h"

ClassImp(AliRsnReader)

//_____________________________________________________________________________
AliRsnReader::AliRsnReader() :
  fCheckSplit(kFALSE),
  fRejectFakes(kFALSE),
  fTPCOnly(kFALSE),
  fPIDDef(),
  fITSClusters(0),
  fTPCClusters(0),
  fTRDClusters(0),
  fTrackRefs(0),
  fTrackRefsITS(0),
  fTrackRefsTPC(0)
{
//
// Constructor.
//
}

//_____________________________________________________________________________
Bool_t AliRsnReader::AreSplitted(AliESDtrack *track1, AliESDtrack *track2)
{
//
// Checks if two tracks are splitted.
// Currently, two split tracks are the ones which
// have equal GEANT labels.
// On real data, this criterion will have to change
// into an "experimental" one.
//

  Int_t lab1 = TMath::Abs(track1->GetLabel());
  Int_t lab2 = TMath::Abs(track2->GetLabel());

  return (lab1 == lab2);
}

//_____________________________________________________________________________
Bool_t AliRsnReader::ResolveSplit(AliESDtrack *track1, AliESDtrack *track2)
{
//
// Resolves a split of two tracks.
// Only two values can be returned:
//  kTRUE  --> accept "track1" and reject "track2"
//  kFALSE --> accept "track2" and reject "track1"
//

  Double_t chiSq1 = track1->GetConstrainedChi2();
  Double_t chiSq2 = track2->GetConstrainedChi2();

  if (chiSq1 < chiSq2) return 1; else return 2;
}

//_____________________________________________________________________________
void AliRsnReader::SetTPCOnly(Bool_t doit)
{
//
// Set the flag for TPC only.
// If this is true, exclude all other detectors from PID
//

  fTPCOnly = doit;

  if (fTPCOnly) {
    fPIDDef.ExcludeAll();
    fPIDDef.IncludeDet(AliRsnPIDDefESD::kTPC);
    fPIDDef.SetDivValue(AliRsnPIDDefESD::kTPC, 0.0);
  }
}

//_____________________________________________________________________________
Bool_t AliRsnReader::Fill
(AliRsnEvent *rsn, AliVEvent *event, AliMCEvent *mc)
{
//
// According to the class type of event and the selected source
// recalls one of the private reading methods to fill the RsnEvent
// passed by reference as first argument.
//

  Bool_t success = kFALSE;
  TString str(event->ClassName());

  if (!str.CompareTo("AliESDEvent")) {
    AliDebug(1, "Reading an ESD event");
    success = FillFromESD(rsn, (AliESDEvent*)event, mc);
  }
  else if (!str.CompareTo("AliAODEvent")) {
    AliDebug(1, "Reading an AOD event");
    success = FillFromAOD(rsn, (AliAODEvent*)event, mc);
  }
  else if (!str.CompareTo("AliMCEvent")) {
    AliDebug(1, "Reading an MC event");
    success = FillFromMC(rsn, (AliMCEvent*)event);
  }
  else {
    AliError(Form("ClassName '%s' not recognized as possible source data: aborting.", str.Data()));
    return kFALSE;
  }

  return success;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromESD(AliRsnEvent *rsn, AliESDEvent *esd, AliMCEvent *mc)
{
//
// Filler from an ESD event.
// Stores all tracks (if a filter is defined, it will store
// only the ones which survive the cuts).
// If a reference MC event is provided, it is used to store
// the MC informations for each track (true PDG code,
// GEANT label of mother, PDG code of mother, if any).
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kESD'.
//

  // retrieve stack (if possible)
  AliStack *stack = 0x0;
  if (mc) stack = mc->Stack();

  // get number of tracks
  Int_t ntracks = esd->GetNumberOfTracks();
  if (!ntracks) return kFALSE;

  // if required with the flag, scans the event
  // and searches all split tracks (= 2 tracks with the same label);
  // for each pair of split tracks, only the better (best chi2) is kept
  // and the other is rejected: this info is stored into a Boolean array
  Int_t i, i1, i2;
  Bool_t *accept = new Bool_t[ntracks];
  for (i = 0; i < ntracks; i++) accept[i] = kTRUE;
  if (fCheckSplit) {
    for (i1 = 0; i1 < ntracks; i1++) {
      AliESDtrack *track1 = esd->GetTrack(i1);
      for (i2 = i1 + 1; i2 < ntracks; i2++) {
        AliESDtrack *track2 = esd->GetTrack(i2);
        if (AreSplitted(track1, track2)) {
          if (ResolveSplit(track1, track2)) accept[i2] = kFALSE;
          else accept[i1] = kFALSE;
        }
      }
    }
  }

  // get primary vertex
  Double_t vertex[3];
  if (!fTPCOnly) {
    // when taking vertex from ESD event there are two options:
    // if a vertex with tracks was successfully reconstructed,
    // it is used for computing DCA;
    // otherwise, the one computed with SPD is used.
    // This is known from the "Status" parameter of the vertex itself.
    const AliESDVertex *v = esd->GetPrimaryVertex();
    if (!v->GetStatus()) v = esd->GetPrimaryVertexSPD();

    // get primary vertex
    vertex[0] = (Double_t)v->GetXv();
    vertex[1] = (Double_t)v->GetYv();
    vertex[2] = (Double_t)v->GetZv();
  }
  else {
    const AliESDVertex *v = esd->GetPrimaryVertexTPC();

    // get primary vertex
    vertex[0] = (Double_t)v->GetXv();
    vertex[1] = (Double_t)v->GetYv();
    vertex[2] = (Double_t)v->GetZv();
  }
  rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);

  // store tracks from ESD
  Float_t p[2], cov[3];
  Int_t   index, label, labmum;
  Bool_t  check;
  AliRsnDaughter temp;
  AliESDtrack *esdTrack = 0x0, esdTrackTmp;
  for (index = 0; index < ntracks; index++) {
    // skip track recognized as the worse one in a splitted pair
    if (!accept[index]) {
      AliDebug(10, Form("Rejecting split track #%d in this event", index));
      continue;
    }
    // get ESD track
    esdTrack = esd->GetTrack(index);
    // check for fake tracks
    label = esdTrack->GetLabel();
    if (fRejectFakes && (label < 0)) continue;
    // check for minimum number of clusters in ITS, TPC and TRD
    if (fITSClusters > 0 && (esdTrack->GetITSclusters(0x0) < fITSClusters)) if (!fTPCOnly) continue;
    if (fTPCClusters > 0 && (esdTrack->GetTPCclusters(0x0) < fTPCClusters)) continue;
    if (fTRDClusters > 0 && (esdTrack->GetTRDclusters(0x0) < fTRDClusters)) if (!fTPCOnly) continue;
    // switch to TPC data if required
    if (fTPCOnly) {
      esdTrack->GetImpactParametersTPC(p, cov);
      if (p[0] == 0. && p[1] == 0.) {
        esdTrack->RelateToVertexTPC(esd->GetPrimaryVertexTPC(), esd->GetMagneticField(), kVeryBig);
      }
      check = esdTrack->FillTPCOnlyTrack(esdTrackTmp);
      if (check) esdTrack = &esdTrackTmp; else continue;
    }
    // try to get information from this track
    // output value tells if this was successful or not
    check = temp.Adopt(esdTrack, fPIDDef);
    if (!check) {
      AliDebug(10, Form("Failed adopting track #%d", index));
      continue;
    }

    // if stack is present, copy MC info
    if (stack) {
      TParticle *part = stack->Particle(TMath::Abs(label));
      if (part) {
        temp.InitMCInfo(part);
        labmum = part->GetFirstMother();
        if (labmum >= 0) {
          TParticle *mum = stack->Particle(labmum);
          temp.GetMCInfo()->SetMotherPDG(mum->GetPdgCode());
        }
      }
    }

    // set index and label
    temp.SetIndex(index);
    temp.SetLabel(label);

    // shifts the track DCA to the found vertex
    temp.ShiftZero(vertex[0], vertex[1], vertex[2]);

    // try to add track to collection and returns an error in case of problems
    AliRsnDaughter *ptr = rsn->AddTrack(temp);
    if (!ptr) AliWarning(Form("Failed storing track#%d", index));
  }

  // compute total multiplicity
  rsn->MakeComputations();
  if (rsn->GetMultiplicity() <= 0) {
      AliDebug(1, "Zero Multiplicity in this event");
      return kFALSE;
  }

  // sort tracks w.r. to Pt (from largest to smallest)
  rsn->SortTracks();
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromAOD(AliRsnEvent *rsn, AliAODEvent *aod, AliMCEvent *mc)
{
//
// Filler from an AOD event.
// Stores all tracks (if a filter is defined, it will store
// only the ones which survive the cuts).
// If a reference MC event is provided, it is used to store
// the MC informations for each track (true PDG code,
// GEANT label of mother, PDG code of mother, if any).
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kAOD'.
//

    // retrieve stack (if possible)
    AliStack *stack = 0x0;
    if (mc) stack = mc->Stack();

    // get number of tracks
    Int_t ntracks = aod->GetNTracks();
    if (!ntracks) return kFALSE;

    // get primary vertex
    Double_t vertex[3];
    vertex[0] = aod->GetPrimaryVertex()->GetX();
    vertex[1] = aod->GetPrimaryVertex()->GetY();
    vertex[2] = aod->GetPrimaryVertex()->GetZ();
    rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);

    // store tracks from ESD
    Int_t  index, label, labmum;
    Bool_t check;
    AliAODTrack *aodTrack = 0;
    AliRsnDaughter temp;
    TObjArrayIter iter(aod->GetTracks());
    while ((aodTrack = (AliAODTrack*)iter.Next()))
    {
        // retrieve index
        index = aod->GetTracks()->IndexOf(aodTrack);
        label = aodTrack->GetLabel();
        if (fRejectFakes && (label < 0)) continue;
        // copy ESD track data into RsnDaughter
        // if unsuccessful, this track is skipped
        check = temp.Adopt(aodTrack);
        if (!check) continue;
        // if stack is present, copy MC info
        if (stack)
        {
            TParticle *part = stack->Particle(TMath::Abs(label));
            if (part)
            {
                temp.InitMCInfo(part);
                labmum = part->GetFirstMother();
                if (labmum >= 0)
                {
                    TParticle *mum = stack->Particle(labmum);
                    temp.GetMCInfo()->SetMotherPDG(mum->GetPdgCode());
                }
            }
        }
        // set index and label and add this object to the output container
        temp.SetIndex(index);
        temp.SetLabel(label);
        AliRsnDaughter *ptr = rsn->AddTrack(temp);
        // if problems occurred while storin, that pointer is NULL
        if (!ptr) AliWarning(Form("Failed storing track#%d"));
    }

    // compute total multiplicity
    rsn->MakeComputations();
    if (rsn->GetMultiplicity() <= 0)
    {
        AliDebug(1, "Zero multiplicity in this event");
        return kFALSE;
    }

    return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromMC(AliRsnEvent *rsn, AliMCEvent *mc)
{
//
// Filler from an ESD event.
// Stores all tracks which generate at least one
// TrackReference (a point in a sensitive volume).
// In this case, the MC info is stored by default and
// perfect particle identification is the unique available.
// When this is used, the 'source' flag of the output
// AliRsnEvent object will be set to 'kMC'.
//

  // get number of tracks
  Int_t ntracks = mc->GetNumberOfTracks();
  if (!ntracks) return kFALSE;

  AliStack *stack = mc->Stack();

  // get primary vertex
  TArrayF fvertex(3);
  Double_t vertex[3];
  mc->GenEventHeader()->PrimaryVertex(fvertex);
  vertex[0] = (Double_t)fvertex[0];
  vertex[1] = (Double_t)fvertex[1];
  vertex[2] = (Double_t)fvertex[2];
  rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);

  // store tracks from MC
  Int_t  i, index, labmum, nHitsITS, nHitsTPC, nRef;
  Bool_t check;
  AliRsnDaughter temp;
  for (index = 0; index < ntracks; index++) {
    // get and check MC track
    AliMCParticle *mcTrack = mc->GetTrack(index);
    // check particle track references
    nRef = mcTrack->GetNumberOfTrackReferences();
    if (fTrackRefs > 0 && nRef < fTrackRefs) continue;
    else if (fTrackRefsITS > 0 || fTrackRefsTPC > 0) {
      nHitsITS = nHitsTPC = 0;
      for (i = 0; i < nRef; i++) {
        AliTrackReference *trackRef = mcTrack->GetTrackReference(i);
        if (trackRef) {
          Int_t detectorId = trackRef->DetectorId();
          switch(detectorId) {
            case AliTrackReference::kITS  : nHitsITS++  ; break ;
            case AliTrackReference::kTPC  : nHitsTPC++  ; break ;
            default : break ;
          }
        }
      }
      if (fTrackRefsITS > 0 && nHitsITS < fTrackRefsITS) continue;
      if (fTrackRefsTPC > 0 && nHitsTPC < fTrackRefsTPC) continue;
    }
    // try to insert in the RsnDaughter its data
    check = temp.Adopt(mcTrack);
    if (!check) continue;
    labmum = temp.GetMCInfo()->Mother();
    if (labmum >= 0)
    {
      TParticle *mum = stack->Particle(labmum);
      temp.GetMCInfo()->SetMotherPDG(mum->GetPdgCode());
    }
    // if successful, set other data and stores it
    temp.SetIndex(index);
    temp.SetLabel(mcTrack->Label());
    AliRsnDaughter *ptr = rsn->AddTrack(temp);
    // if problems occurred while storin, that pointer is NULL
    if (!ptr) AliWarning(Form("Failed storing track#%d", index));
  }

  // compute total multiplicity
  rsn->MakeComputations();
  if (rsn->GetMultiplicity() <= 0)
  {
    AliDebug(1, "Zero multiplicity in this event");
    return kFALSE;
  }

  // sort tracks w.r. to Pt (from largest to smallest)
  rsn->SortTracks();
  return kTRUE;
}
