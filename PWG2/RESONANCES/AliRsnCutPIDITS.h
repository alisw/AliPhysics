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

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDpid.h"

#include "AliAODpidUtil.h"

#include "AliRsnDaughter.h"
#include "AliRsnCut.h"

class AliRsnCutPIDITS : public AliRsnCut {
public:

   AliRsnCutPIDITS(const char *name          = "cutITS",
                   EPARTYPE    ref           = AliPID::kKaon,
                   Double_t    nSigmaMin     = -3.,
                   Double_t    nSigmaMax     =  3.,
                   Bool_t      rejectOutside = kTRUE);

   AliRsnCutPIDITS(const AliRsnCutPIDITS& copy);
   AliRsnCutPIDITS& operator=(const AliRsnCutPIDITS& copy);
   virtual ~AliRsnCutPIDITS() { }

   AliESDpid*       ESDpid()  {return fESDpid;}
   AliAODpidUtil*   AODpid()  {return fAODpid;}

   void             SetMC(Bool_t mc = kTRUE)                      {fIsMC = mc;}
   void             SetRejectOutside(Bool_t yn = kTRUE)           {fRejectOutside = yn;}
   void             SetMomentumRange(Double_t min, Double_t max)  {fMomMin = min; fMomMax = max;}
   void             SetNSigmaRange(Double_t min, Double_t max)    {fMinD = min; fMaxD = max;}
   void             SetRefType(EPARTYPE type)                     {fRefType = type;}

   Bool_t           IsTPC(AliVTrack *vtrack);
   Bool_t           IsSA(AliVTrack *vtrack);
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

private:

   Bool_t          fIsMC;           //  needed to initialize the pid object
   Bool_t          fRejectOutside;  //  choose if tracks outside momentum range are rejected or not
   Double_t        fMomMin;         //  min p in range where this cut is checked
   Double_t        fMomMax;         //  max p in range where this cut is checked
   EPARTYPE        fRefType;        //  particle type for which PID is checked
   AliESDpid      *fESDpid;         //! ESD PID object
   AliAODpidUtil  *fAODpid;         //! AOD PID object

   ClassDef(AliRsnCutPIDITS, 1)
};

inline Bool_t AliRsnCutPIDITS::IsTPC(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for a global track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTPCin = ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);

   return (isTPCin);
}

inline Bool_t AliRsnCutPIDITS::IsSA(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTPCin     = ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);
   Bool_t isITSrefit  = ((vtrack->GetStatus() & AliESDtrack::kITSrefit) != 0);
   Bool_t isITSpureSA = ((vtrack->GetStatus() & AliESDtrack::kITSpureSA) != 0);
   Bool_t isITSpid    = ((vtrack->GetStatus() & AliESDtrack::kITSpid) != 0);

   return ((!isTPCin) && isITSrefit && (!isITSpureSA) && isITSpid);

   return kTRUE;
}

#endif
