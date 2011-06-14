#ifndef ALIRSNCUTPROTON2010PP_H
#define ALIRSNCUTPROTON2010PP_H

//
// All cuts for single Protons in phi analysis 2010,
// based on quality and PID using the TPC and TOF
// detectors, using default definitions for both
// kinds of cuts, for ESD and AOD
// Author: Serguey Kiselev
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutProton2010PP : public AliRsnCut {

public:

   AliRsnCutProton2010PP(const char *name = "");
   virtual ~AliRsnCutProton2010PP();
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}

private:

   Bool_t MatchTOF(const AliVTrack *vtrack);

   AliRsnCutTrackQuality fCutQuality;  // track quality cut

   ClassDef(AliRsnCutProton2010PP,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutProton2010PP::MatchTOF(const AliVTrack *vtrack)
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
