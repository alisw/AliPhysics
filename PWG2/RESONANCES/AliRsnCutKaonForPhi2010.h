#ifndef ALIRSNCUTKaonForPhi2010_H
#define ALIRSNCUTKaonForPhi2010_H

//
// All cuts for single kaons in phi analysis 2010
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutKaonForPhi2010 : public AliRsnCut {

public:

   AliRsnCutKaonForPhi2010(const char *name = "");
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}

private:

   Bool_t MatchTOF(AliVTrack *vtrack);

   AliRsnCutTrackQuality fCutQuality;  // track quality cut

   ClassDef(AliRsnCutKaonForPhi2010,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutKaonForPhi2010::MatchTOF(AliVTrack *vtrack)
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
