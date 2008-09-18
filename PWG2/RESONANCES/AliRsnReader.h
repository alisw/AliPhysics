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

class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliRsnEvent;
class AliRsnPIDWeightsMgr;

class AliRsnReader : public TNamed
{
  public:

    AliRsnReader(AliRsnPIDWeightsMgr *mgr = 0x0);
    virtual ~AliRsnReader() {}

    void    SetWeightsMgr(AliRsnPIDWeightsMgr *mgr) {fWeightsMgr = mgr;}
    void    SetCheckSplit(Bool_t doit = kTRUE) {fCheckSplit = doit;}
    void    SetRejectFakes(Bool_t doit = kTRUE) {fRejectFakes = doit;}

    Bool_t  CheckSplit() {return fCheckSplit;}
    Bool_t  RejectFakes() {return fRejectFakes;}

    Bool_t  Fill(AliRsnEvent *rsn, AliVEvent *event, AliMCEvent *refMC = 0, Bool_t useTPCOnly = kFALSE);
    Bool_t  FillFromESD(AliRsnEvent *rsn, AliESDEvent *event, AliMCEvent *refMC = 0, Bool_t useTPCOnly = kFALSE);
    Bool_t  FillFromAOD(AliRsnEvent *rsn, AliAODEvent *event, AliMCEvent *refMC = 0);
    Bool_t  FillFromMC(AliRsnEvent *rsn, AliMCEvent *mc);

    // sets PID TYPE
    void SetPIDtype(const AliRsnDaughter::EPIDType& theValue = AliRsnDaughter::kEsd, Double_t divValue = 0.0);
    AliRsnDaughter::EPIDType GetPIDtype() const { return fCurrentPIDtype; }

    // sets Sectors cut
    void SetITSTPCTRDSectors(const Int_t& its = -1, const Int_t& tpc = -1, const Int_t& trd = -1);

  protected:

    // dummy copy methods
    AliRsnReader(const AliRsnReader &copy) : TNamed(copy),
        fCheckSplit(0),fRejectFakes(0),fWeightsMgr(0x0),fCurrentPIDtype(AliRsnDaughter::kEsd),
        fPIDDivValue(0.),fITSClusters(0),fTPCClusters(0),fTRDClusters(0) { /*nothing*/ }
    AliRsnReader& operator=(const AliRsnReader&) {return (*this);}

    Bool_t                    fCheckSplit;     // flag to check and remove split tracks
    Bool_t                    fRejectFakes;    // flag to reject fake tracks (negative label)

    AliRsnPIDWeightsMgr      *fWeightsMgr;     // manager for alternative PID weights
    AliRsnDaughter::EPIDType  fCurrentPIDtype; // PID type
    Double_t                  fPIDDivValue;    // PID div Value
    Int_t                     fITSClusters;    //
    Int_t                     fTPCClusters;    //
    Int_t                     fTRDClusters;    //

  private:

    ClassDef(AliRsnReader, 1);
};

#endif
