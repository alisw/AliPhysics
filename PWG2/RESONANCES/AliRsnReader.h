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

#ifndef ALIRSNREADER_H
#define ALIRSNREADER_H

#include <TNamed.h>

#include "AliRsnDaughter.h"
#include "AliRsnPIDDefESD.h"

class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliRsnEvent;

class AliRsnReader : public TObject
{
  public:

    AliRsnReader();
    virtual ~AliRsnReader() {}

    void    SetCheckSplit(Bool_t doit = kTRUE) {fCheckSplit = doit;}
    Bool_t  AreSplitted(AliESDtrack *track1, AliESDtrack *track2);
    Bool_t  ResolveSplit(AliESDtrack *track1, AliESDtrack *track2);
    Bool_t  DoesCheckSplit() {return fCheckSplit;}

    void    SetRejectFakes(Bool_t doit = kTRUE) {fRejectFakes = doit;}
    Bool_t  DoesRejectFakes() {return fRejectFakes;}

    void    SetTPCOnly(Bool_t doit = kTRUE);
    Bool_t  DoesTPCOnly() {return fTPCOnly;}

    AliRsnPIDDefESD& GetPIDDef() {return fPIDDef;}

    void SetTrackRefs(Int_t value) {fTrackRefs = value;}
    void SetTrackRefsITS(Int_t value) {fTrackRefsITS = value;}
    void SetTrackRefsTPC(Int_t value) {fTrackRefsTPC = value;}

    void SetMinITSClusters(Int_t value) {fITSClusters = value;}
    void SetMinTPCClusters(Int_t value) {fTPCClusters = value;}
    void SetMinTRDClusters(Int_t value) {fTRDClusters = value;}
    void SetITSTPCTRDSectors(const Int_t& its = -1, const Int_t& tpc = -1, const Int_t& trd = -1);

    Bool_t  Fill(AliRsnEvent *rsn, AliVEvent *event, AliMCEvent *refMC = 0);
    Bool_t  FillFromESD(AliRsnEvent *rsn, AliESDEvent *event, AliMCEvent *refMC = 0);
    Bool_t  FillFromAOD(AliRsnEvent *rsn, AliAODEvent *event, AliMCEvent *refMC = 0);
    Bool_t  FillFromMC(AliRsnEvent *rsn, AliMCEvent *mc);

  private:

    // dummy copy methods
    AliRsnReader(const AliRsnReader &copy) :
      TObject(copy),fCheckSplit(0),fRejectFakes(0),fTPCOnly(0),fPIDDef(copy.fPIDDef),
      fITSClusters(0),fTPCClusters(0),fTRDClusters(0),
      fTrackRefs(0),fTrackRefsITS(0),fTrackRefsTPC(0) { /*nothing*/ }
    AliRsnReader& operator=(const AliRsnReader&) {return (*this);}

    Bool_t          fCheckSplit;     // flag to check and remove split tracks
    Bool_t          fRejectFakes;    // flag to reject fake tracks (negative label)
    Bool_t          fTPCOnly;        // flag to use only the TPC for reading data

    AliRsnPIDDefESD fPIDDef;         // manager for alternative PID weights (ESD only)

    Int_t           fITSClusters;    // minimum number of ITS clusters to accept a track
    Int_t           fTPCClusters;    // minimum number of TPC clusters to accept a track
    Int_t           fTRDClusters;    // minimum number of TRD clusters to accept a track

    Int_t           fTrackRefs;      // minimum required track references for MC reading
    Int_t           fTrackRefsITS;   // minimum required track references for MC reading (ITS)
    Int_t           fTrackRefsTPC;   // minimum required track references for MC reading (TPC)

    ClassDef(AliRsnReader, 1);
};

#endif
