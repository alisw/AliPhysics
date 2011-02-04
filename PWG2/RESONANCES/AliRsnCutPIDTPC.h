//
// *** Class AliRsnCutPIDTPC ***
//
// This class implements all cuts which have to be used for the 2010 runs
// for phi and generic resonance analysis.
// It contains an AliESDtrackCuts object for track quality selection
// and some criteria for particle identification with ITS, TPC and TOF.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTPIDTPC_H
#define ALIRSNCUTPIDTPC_H

#include "AliPID.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"

#include "AliAODpidUtil.h"

#include "AliRsnDaughter.h"
#include "AliRsnCut.h"

class AliRsnCutPIDTPC : public AliRsnCut
{
  public:

    AliRsnCutPIDTPC(const char *name       = "cutTPC",
                    EPARTYPE type          = AliPID::kKaon,
                    Double_t nSigmaMin     = -3.,
                    Double_t nSigmaMax     =  3.,
                    Bool_t   rejectOutside = kTRUE);
                    
    AliRsnCutPIDTPC(const AliRsnCutPIDTPC& copy);
    AliRsnCutPIDTPC& operator=(const AliRsnCutPIDTPC& copy);
    virtual ~AliRsnCutPIDTPC() { }

    AliESDpid*       ESDpid()  {return &fESDpid;}
    AliAODpidUtil*   AODpid()  {return &fAODpid;}

    void             SetRejectOutside(Bool_t yn = kTRUE)           {fRejectOutside = yn;}
    void             SetMomentumRange(Double_t min, Double_t max)  {fMomMin = min; fMomMax = max;}
    void             SetNSigmaRange(Double_t min, Double_t max)    {fMinD = min; fMaxD = max;}
    void             SetRefType(EPARTYPE type)                     {fRefType = type;}
    void             SetBBParam(Double_t *p)                       {SetBBParam(p[0], p[1], p[2], p[3], p[4]);}
    void             SetBBParam(Double_t p0, Double_t p1, Double_t p2, Double_t p3, Double_t p4);

    Bool_t           IsTPC(AliVTrack *vtrack);
    virtual Bool_t   IsSelected(TObject *object);
    virtual void     Print(const Option_t *option = "") const;

  private:

    void Initialize();

    Bool_t          fInitialized;    // a mono-usage flag which initializes the ESD pid object
    Bool_t          fRejectOutside;  // choose if tracks outside momentum range are rejected or not
    Double_t        fMomMin;         // min p in range where this cut is checked
    Double_t        fMomMax;         // max p in range where this cut is checked
    EPARTYPE        fRefType;        // particle type for which PID is checked
    AliESDpid       fESDpid;         // ESD PID object
    AliAODpidUtil   fAODpid;         // AOD PID object

    ClassDef(AliRsnCutPIDTPC, 1)
};

inline Bool_t AliRsnCutPIDTPC::IsTPC(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for a global track
//

  if (!vtrack)
  {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  
  Bool_t isTPCin     = ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);
  Bool_t isITSrefit  = ((vtrack->GetStatus() & AliESDtrack::kITSrefit) != 0);
  Bool_t isITSpid    = ((vtrack->GetStatus() & AliESDtrack::kITSpid) != 0);
  
  return ( isTPCin && isITSrefit && isITSpid );
  
  return kTRUE;
}

#endif
