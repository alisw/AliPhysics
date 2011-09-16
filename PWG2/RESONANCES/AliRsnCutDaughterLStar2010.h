#ifndef ALIRSNCUTDAUGHTERLSTAR2010_H
#define ALIRSNCUTDAUGHTERLSTAR2010_H

//
// Cuts for selecting good proton and kaon candidates for Lambda(1520) analysis
// with the data samples from PbPb runs in 2010.
// Applies track quality selection plus PID selection,
// with different tolerance ranges depending on the momentum.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutDaughterLStar2010 : public AliRsnCut {

public:

   AliRsnCutDaughterLStar2010(const char *name = "", AliPID::EParticleType pid = AliPID::kKaon);
   virtual ~AliRsnCutDaughterLStar2010() { }
   
   void                   SetPID(AliPID::EParticleType type) {fPID = type;}
   AliRsnCutTrackQuality *CutQuality()                       {return &fCutQuality;}
   Bool_t                 MatchTOF(const AliVTrack *vtrack);
   virtual Bool_t         IsSelected(TObject *obj);

private:

   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutTrackQuality fCutQuality;       // track quality cut

   ClassDef(AliRsnCutDaughterLStar2010,1) // cut definitions for L*

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutDaughterLStar2010::MatchTOF(const AliVTrack *vtrack)
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
