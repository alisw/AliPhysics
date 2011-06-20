#ifndef ALIRSNCUTKAONFORPHI2010_H
#define ALIRSNCUTKAONFORPHI2010_H

//
// This cut implements all the checks done to accept a track as a Kaon
// for the PbPb analysis using 2010 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutKaonForPhi2010 : public AliRsnCut {

public:

   AliRsnCutKaonForPhi2010(const char *name = "", Double_t nSigmaTPC = 3.0, Double_t nSigmaTOF = 3.0, Double_t tofLimit = 0.8);
   virtual ~AliRsnCutKaonForPhi2010() { }
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality()            {return &fCutQuality;}
   Double_t               GetTOFthreshold() const {return fTOFthreshold;}
   
   void   SetTOFthreshold(Double_t value)   {fTOFthreshold = value;}
   void   SetOnlyQuality(Bool_t yn = kTRUE) {fOnlyQuality = yn;}
   void   SetOnlyTPC(Bool_t yn = kTRUE)     {fOnlyTPC = yn;}
   void   SetOnlyTOF(Bool_t yn = kTRUE)     {fOnlyTOF = yn;}
   void   SetCutTPC(Double_t cut)           {fCutTPC = cut;}
   void   SetCutTOF(Double_t cut)           {fCutTOF = cut;}
   Bool_t MatchTOF(const AliVTrack *vtrack);

private:

   Bool_t                fOnlyQuality;   // switch off PID
   Bool_t                fOnlyTPC;       // use only TPC PID
   Bool_t                fOnlyTOF;       // use only TOF PID
   Double_t              fCutTPC;        // nsigma cut for TPC
   Double_t              fCutTOF;        // nsigma cut fof TOF
   Double_t              fTOFthreshold;  // for Pt above this threshold, TOF is mandatory
   AliRsnCutTrackQuality fCutQuality;    // track quality cut

   ClassDef(AliRsnCutKaonForPhi2010,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutKaonForPhi2010::MatchTOF(const AliVTrack *vtrack)
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
