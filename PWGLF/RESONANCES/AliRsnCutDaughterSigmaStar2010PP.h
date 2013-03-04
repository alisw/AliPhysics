#ifndef ALIRSNCUTDAUGHTERSIGMASTAR2010PP_H
#define ALIRSNCUTDAUGHTERSIGMASTAR2010PP_H

//
// Cuts for selecting good pion candidates for Sigma* analysis
// with the data samples from pp runs in 2010.
// Applies track quality selection plus PID selection,
// with different tolerance ranges depending on the momentum.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutDaughterSigmaStar2010PP : public AliRsnCut {

public:

   AliRsnCutDaughterSigmaStar2010PP(const char *name = "", AliPID::EParticleType pid = AliPID::kPion);
   virtual ~AliRsnCutDaughterSigmaStar2010PP() { }

   //void                   SetPID(AliPID::EParticleType type) {fPID = type;}
   AliRsnCutTrackQuality *CutQuality()                       {return &fCutQuality;}
   Bool_t                 MatchTOF(const AliVTrack *vtrack);
   virtual Bool_t         IsSelected(TObject *obj);

   void           SetPIDCut(Double_t value)                  {fPIDCut = value;}
   void           SetMinTPCcluster(Int_t value)              {fMinTPCcluster = value;}
   
private:

   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutTrackQuality fCutQuality;       // track quality cut

   ClassDef(AliRsnCutDaughterSigmaStar2010PP,2) // cut definitions for Sigma*

protected:

   Double_t         fPIDCut;          // nsigmas for pions
   Int_t            fMinTPCcluster;   // min allowed TPC cluster
};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutDaughterSigmaStar2010PP::MatchTOF(const AliVTrack *vtrack)
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
