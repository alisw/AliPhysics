#ifndef ALIRSNCUTDAUGHTERKSTAR2010PP_H
#define ALIRSNCUTDAUGHTERKSTAR2010PP_H

//
// Cuts for selecting good pion candidates for K* analysis
// with the data samples from pp runs in 2010.
// Applies track quality selection plus PID selection,
// with different tolerance ranges depending on the momentum.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutDaughterKStar2010PP : public AliRsnCut {

public:

   AliRsnCutDaughterKStar2010PP(const char *name = "", AliPID::EParticleType pid = AliPID::kPion);
   virtual ~AliRsnCutDaughterKStar2010PP() { }
   
   void                   SetNoPID(Bool_t yn = kTRUE)        {fNoPID = yn;}
   
   void                   SetPID(AliPID::EParticleType type) {fPID = type;}
   AliRsnCutTrackQuality *CutQuality()                       {return &fCutQuality;}
   Bool_t                 MatchTOF(const AliVTrack *vtrack);
   virtual Bool_t         IsSelected(TObject *obj);

private:

   Bool_t                fNoPID;            // flag to switch off PID check
   
   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutTrackQuality fCutQuality;       // track quality cut

   ClassDef(AliRsnCutDaughterKStar2010PP,1) // cut definitions for K*

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutDaughterKStar2010PP::MatchTOF(const AliVTrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
   if (!(vtrack->GetStatus() & AliESDtrack::kTIME  )) return kFALSE;

   return kTRUE;
}

#endif
