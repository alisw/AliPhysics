//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTPIDTOF_H
#define ALIRSNCUTPIDTOF_H

#include "AliPID.h"
#include "AliVTrack.h"

#include "AliRsnCut.h"

class AliRsnCutPIDTOF : public AliRsnCut {
public:

   AliRsnCutPIDTOF(const char *name            = "cutTOF",
                   EPARTYPE    particle        = AliPID::kKaon,
                   Double_t    nSigmaMin       = -3.,
                   Double_t    nSigmaMax       =  3.,
                   Bool_t      rejectUnmatched = kFALSE);

   AliRsnCutPIDTOF(const AliRsnCutPIDTOF& copy);
   AliRsnCutPIDTOF& operator=(const AliRsnCutPIDTOF& copy);
   virtual ~AliRsnCutPIDTOF() { }

   AliESDpid*      ESDpid()  {return fESDpid;}
   AliAODpidUtil*  AODpid()  {return fAODpid;}

   void            SetRejectUnmatched(Bool_t yn = kTRUE)      {fRejectUnmatched = yn;}
   void            SetNSigmaRange(Double_t min, Double_t max) {fMinD = min; fMaxD = max;}
   void            SetRefType(EPARTYPE type)                  {fRefType = type; fRefMass = AliPID::ParticleMass(type);}

   Bool_t          IsMatched(AliVTrack *vtrack);
   virtual Bool_t  IsSelected(TObject *object);
   virtual void    Print(const Option_t *option = "") const;

private:

   void Initialize();

   Bool_t            fRejectUnmatched;  //  decide if non TOF matched tracks pass the cut or not
   EPARTYPE          fRefType;          //  particle type for which PID is checked
   Double_t          fRefMass;          //  reference mass used for computations
   AliESDpid        *fESDpid;           //! PID utility for ESD
   AliAODpidUtil    *fAODpid;           //! PID utility for AOD

   ClassDef(AliRsnCutPIDTOF, 1)
};

inline Bool_t AliRsnCutPIDTOF::IsMatched(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTOFout = ((vtrack->GetStatus() & AliESDtrack::kTOFout) != 0);
   Bool_t isTIME   = ((vtrack->GetStatus() & AliESDtrack::kTIME) != 0);

   // if flags are not set, track is not matched
   if (!isTOFout || !isTIME) return kFALSE;

   // do an additional check on integrated length for ESD tracks
   AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*>(vtrack);
   if (esdTrack) if (esdTrack->GetIntegratedLength() < 350.) return kFALSE;

   // if we are here, flags are OK and length also
   return kTRUE;
}

#endif
