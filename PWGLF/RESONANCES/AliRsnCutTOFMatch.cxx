#include "AliRsnCutTOFMatch.h"

ClassImp(AliRsnCutTOFMatch)

AliRsnCutTOFMatch::AliRsnCutTOFMatch() :
   AliRsnCut("cut", AliRsnTarget::kDaughter)
{
   //Default constructor
}

//_________________________________________________________________________________________________
AliRsnCutTOFMatch::AliRsnCutTOFMatch(const char *name) :
   AliRsnCut(name, AliRsnTarget::kDaughter)
{
   //main constructor
}

//_________________________________________________________________________________________________
Bool_t AliRsnCutTOFMatch::IsSelected(TObject *object)
{
//
// Checks if the track has the status flags required for the TPC-TOF matching
//
   Bool_t accept = kFALSE;
   if (!TargetOK(object)) return accept;

   AliVTrack *vtrack = fDaughter->Ref2Vtrack();
   if (!vtrack) return accept;

   if (MatchTOF(vtrack)) accept = kTRUE;
   return accept;
}

//-----------------------------------------------------
inline Bool_t AliRsnCutTOFMatch::MatchTOF(const AliVTrack *vtrack) const
{
//
// Checks if the track has matched the TOF detector
//
   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }
   if ((vtrack->GetStatus() & AliESDtrack::kTOFout) == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kTIME  ) == 0) return kFALSE;
   return kTRUE;
}
