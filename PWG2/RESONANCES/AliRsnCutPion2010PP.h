#ifndef ALIRSNCUTPION2010PP_H
#define ALIRSNCUTPION2010PP_H

//
// All cuts for single pions in phi analysis 2010,
// based on track quality and particle identification
// with TPC and TOF.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutPion2010PP : public AliRsnCut {

public:

   AliRsnCutPion2010PP(const char *name = "");
   virtual ~AliRsnCutPion2010PP() { }
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}

private:

   Bool_t MatchTOF(const AliVTrack *vtrack);

   AliRsnCutTrackQuality fCutQuality;  // track quality cut

   ClassDef(AliRsnCutPion2010PP,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutPion2010PP::MatchTOF(const AliVTrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTOFout = ((vtrack->GetStatus() & AliESDtrack::kTOFout) != 0);
   Bool_t isTIME   = ((vtrack->GetStatus() & AliESDtrack::kTIME) != 0);

   return (isTOFout && isTIME);
}

#endif
