//
// *** Class AliRsnCutPIDITS ***
//
// This class implements all cuts which have to be used for the 2010 runs
// for phi and generic resonance analysis.
// It contains an AliESDtrackCuts object for track quality selection
// and some criteria for particle identification with ITS, ITS and TOF.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTPIDITS_H
#define ALIRSNCUTPIDITS_H

#include "AliPID.h"
#include "AliVTrack.h"
#include "AliESDpid.h"
#include "AliAODpidUtil.h"
#include "AliRsnCut.h"

class AliRsnCutPIDITS : public AliRsnCut
{
  public:

    AliRsnCutPIDITS(const char *name = "cutPIDITS", AliPID::EParticleType type = AliPID::kKaon, Bool_t fIsMC = kFALSE, Double_t momLimit = 0.0, Double_t cut1 = 3.0, Double_t cut2 = 3.0);
    AliRsnCutPIDITS(const AliRsnCutPIDITS& copy);
    AliRsnCutPIDITS& operator=(const AliRsnCutPIDITS& copy);
    virtual ~AliRsnCutPIDITS() { }

    AliESDpid*       GetESDpid()                         {return &fESDpid;}    
    void             SetPIDType(AliPID::EParticleType t) {fPIDtype = t;}
        
    void             SetMC(Bool_t yn = kTRUE);
    void             SetMomentumLimit(Double_t v)        {fMomentumLimit = v;}
    void             SetLargeCut(Double_t v)             {fLargeCut = v;}
    void             SetSmallCut(Double_t v)             {fSmallCut = v;}
    
    virtual Bool_t   IsSelected(TObject *object);
    virtual void     Print(const Option_t *option = "") const;

  protected:
  
    Bool_t  IsITSTPC (AliVTrack *d);  // check that the track is TPC+ITS
    Bool_t  IsITSSA  (AliVTrack *d);  // check that the track is ITS standalone
  
    AliPID::EParticleType   fPIDtype;          //  PID reference type used for checks
    
    Bool_t                  fIsMC;             //  to know what settings must be used
    Double_t                fMomentumLimit;    //  below this value, large cut is used; above, small one is used (value in GeV/c)
    Double_t                fLargeCut;         //  range for ITS de/dx large cut
    Double_t                fSmallCut;         //  range for ITS de/dx small cut
    
    AliESDpid               fESDpid;           //  ESD PID object
    AliAODpidUtil           fAODpid;           //  AOD PID object

    ClassDef(AliRsnCutPIDITS, 1)
};

inline Bool_t AliRsnCutPIDITS::IsITSTPC(AliVTrack *vtrack)
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

inline Bool_t AliRsnCutPIDITS::IsITSSA(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

  if (!vtrack)
  {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  
  Bool_t isTPCin     = ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);
  Bool_t isITSrefit  = ((vtrack->GetStatus() & AliESDtrack::kITSrefit) != 0);
  Bool_t isITSpureSA = ((vtrack->GetStatus() & AliESDtrack::kITSpureSA) != 0);
  Bool_t isITSpid    = ((vtrack->GetStatus() & AliESDtrack::kITSpid) != 0);
  
  return ( (!isTPCin) && isITSrefit && (!isITSpureSA) && isITSpid );
  
  return kTRUE;
}

#endif
