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
#include "AliRsnPIDWeightsMgr.h"

#include "AliRsnReader.h"

ClassImp(AliRsnReader)

//_____________________________________________________________________________
AliRsnReader::AliRsnReader(AliRsnPIDWeightsMgr *mgr) :
    TNamed("RsnReader", ""),
    fCheckSplit(kFALSE),
    fRejectFakes(kFALSE),
    fWeightsMgr(mgr),
    fCurrentPIDtype(AliRsnDaughter::kEsd),
    fPIDDivValue(0.0),
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
Bool_t AliRsnReader::Fill
(AliRsnEvent *rsn, AliVEvent *event, AliMCEvent *mc, Bool_t useTPCOnly)
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
        success = FillFromESD(rsn, (AliESDEvent*)event, mc, useTPCOnly);
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

    // sort tracks w.r. to Pt (from largest to smallest)
    rsn->SortTracks();
    return success;
}

//_____________________________________________________________________________
Bool_t AliRsnReader::FillFromESD
(AliRsnEvent *rsn, AliESDEvent *esd, AliMCEvent *mc, Bool_t useTPCOnly)
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

    // if required with the flag, scans the event
    // and searches all split tracks (= 2 tracks with the same label);
    // for each pair of split tracks, only the better (best chi2) is kept
    // and the other is rejected: this info is stored into a Boolean array
    Int_t i1, i2, lab1, lab2;
    Bool_t *accept = new Bool_t[ntracks];
    for (i1 = 0; i1 < ntracks; i1++) accept[i1] = kTRUE;
    if (fCheckSplit)
    {
        for (i1 = 0; i1 < ntracks; i1++)
        {
            AliESDtrack *trk1 = esd->GetTrack(i1);
            lab1 = TMath::Abs(trk1->GetLabel());
            for (i2 = i1+1; i2 < ntracks; i2++)
            {
                AliESDtrack *trk2 = esd->GetTrack(i2);
                lab2 = TMath::Abs(trk2->GetLabel());
                // check if labels are equal
                if (lab1 == lab2)
                {
                    if (trk1->GetConstrainedChi2() < trk2->GetConstrainedChi2())
                    {
                        accept[i1] = kTRUE;
                        accept[i2] = kFALSE;
                    }
                    else
                    {
                        accept[i1] = kFALSE;
                        accept[i2] = kTRUE;
                    }
                }
            }
        }
    }

    // get primary vertex
    Double_t vertex[3];
    if (!useTPCOnly)
    {
        vertex[0] = esd->GetVertex()->GetXv();
        vertex[1] = esd->GetVertex()->GetYv();
        vertex[2] = esd->GetVertex()->GetZv();
    }
    else
    {
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
    rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);

    // store tracks from ESD
    Int_t    i, index, label, labmum;
    Bool_t check;
    AliRsnDaughter temp;
    for (index = 0; index < ntracks; index++)
    {
        // skip track recognized as the worse one in a splitted pair
        if (!accept[index])
        {
            AliInfo(Form("Rejecting split track #%d in this event", index));
            continue;
        }
        // get and check ESD track
        AliESDtrack *esdTrack = esd->GetTrack(index);
        label = esdTrack->GetLabel();
        if (fRejectFakes && (label < 0)) continue;
        // copy ESD track data into RsnDaughter
        // if unsuccessful, this track is skipped
//         temp.SetCurrentESDPID(fCurrentPIDtype);

        if (fITSClusters>0)
            if (esdTrack->GetITSclusters(0) < fITSClusters) continue;

        if (fTPCClusters>0)
            if (esdTrack->GetTPCclusters(0) < fTPCClusters) continue;

        if (fTRDClusters>0)
            if (esdTrack->GetTRDclusters(0) < fTRDClusters) continue;

        check = temp.Adopt(esdTrack,fCurrentPIDtype,fPIDDivValue, useTPCOnly);
        if (!check) continue;
        // if the AliRsnWeightsMgr object is initialized
        // this means that ESD PID weights are not used
        // and they are computed according to the Weight manager
        if (fWeightsMgr)
        {
            //AliInfo("Using customized weights");
            //AliInfo(Form("ESD weights = %f %f %f %f %f", temp.PID()[0], temp.PID()[1], temp.PID()[2], temp.PID()[3], temp.PID()[4], temp.PID()[5]));
            esdTrack->GetITSpid(fWeightsMgr->GetWeightArray(AliRsnPIDWeightsMgr::kITS));
            esdTrack->GetTPCpid(fWeightsMgr->GetWeightArray(AliRsnPIDWeightsMgr::kTPC));
            esdTrack->GetTRDpid(fWeightsMgr->GetWeightArray(AliRsnPIDWeightsMgr::kTRD));
            esdTrack->GetTOFpid(fWeightsMgr->GetWeightArray(AliRsnPIDWeightsMgr::kTOF));
            esdTrack->GetHMPIDpid(fWeightsMgr->GetWeightArray(AliRsnPIDWeightsMgr::kHMPID));
            for (i = 0; i < AliRsnPID::kSpecies; i++)
            {
                temp.SetPIDWeight(i, fWeightsMgr->GetWeight((AliRsnPID::EType)i, temp.Pt()));
            }
            //AliInfo(Form("Used weights = %f %f %f %f %f", temp.PID()[0], temp.PID()[1], temp.PID()[2], temp.PID()[3], temp.PID()[4], temp.PID()[5]));
        }
        else
        {
            //AliInfo("Using standard ESD weights");
        }

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
        temp.ShiftZero(vertex[0], vertex[1], vertex[2]);
        AliRsnDaughter *ptr = rsn->AddTrack(temp);
        // if problems occurred while storing, that pointer is NULL
        if (!ptr) AliWarning(Form("Failed storing track#%d", index));
    }

    // compute total multiplicity
    rsn->MakeComputations();
    if (rsn->GetMultiplicity() <= 0)
    {
        AliDebug(1, "Zero Multiplicity in this event");
        return kFALSE;
    }

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
    if (!ntracks)
    {
        AliWarning("No tracks in this event");
        return kFALSE;
    }

    // get primary vertex
    Double_t vertex[3];
    vertex[0] = aod->GetPrimaryVertex()->GetX();
    vertex[1] = aod->GetPrimaryVertex()->GetY();
    vertex[2] = aod->GetPrimaryVertex()->GetZ();
    rsn->SetPrimaryVertex(vertex[0], vertex[1], vertex[2]);

    // store tracks from ESD
    Int_t    index, label, labmum;
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
    if (!ntracks)
    {
        AliWarning("No tracks in this event");
        return kFALSE;
    }

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
    Int_t    i, index, labmum, nHitsITS, nHitsTPC, nRef;
    Bool_t check;
    AliRsnDaughter temp;
    for (index = 0; index < ntracks; index++)
    {
        // get and check MC track
        AliMCParticle *mcTrack = mc->GetTrack(index);
        // check particle track references
        nRef = mcTrack->GetNumberOfTrackReferences();
        if (fTrackRefs > 0 && nRef < fTrackRefs) continue;
        else if (fTrackRefsITS > 0 || fTrackRefsTPC > 0)
        {
          nHitsITS = nHitsTPC = 0;
          for (i = 0; i < nRef; i++) {
            AliTrackReference *trackRef = mcTrack->GetTrackReference(i);
            if(trackRef)
            {
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

    return kTRUE;
}

void AliRsnReader::SetPIDtype(const AliRsnDaughter::EPIDType & theValue, Double_t divValue)
{
    fCurrentPIDtype = theValue;
    fPIDDivValue = divValue;
}


void AliRsnReader::SetITSTPCTRDSectors(const Int_t & its, const Int_t & tpc, const Int_t & trd)
{
    fITSClusters = its;
    fTPCClusters = tpc;
    fTRDClusters = trd;
}
