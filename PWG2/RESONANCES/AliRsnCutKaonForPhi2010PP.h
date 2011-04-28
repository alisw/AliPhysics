#ifndef AliRsnCutKaonForPhi2010PP_H
#define AliRsnCutKaonForPhi2010PP_H

//
// All cuts for single kaons in phi analysis 2010
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutKaonForPhi2010PP : public AliRsnCut {

public:

   AliRsnCutKaonForPhi2010PP(const char *name = "");
   
   void           SetTPCNSigmaLow (Double_t v) {fNSigmaTPCLow  = v;}
   void           SetTPCNSigmaHigh(Double_t v) {fNSigmaTPCHigh = v;}
   void           SetTPCLimit(Double_t v)      {fLimitTPC      = v;}
   void           SetTOFNSigma(Double_t v)     {fNSigmaTOF     = v;}
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}

private:

   Bool_t MatchTOF(AliVTrack *vtrack);
   
   Double_t              fNSigmaTPCLow;   // TPC: nsigma cut below limit
   Double_t              fNSigmaTPCHigh;  // TPC: nsigma cut above limit
   Double_t              fLimitTPC;       // TPC: momentum limit 
   Double_t              fNSigmaTOF;      // TOF: nsigma cut (unique)

   AliRsnCutTrackQuality fCutQuality;     // track quality cut

   ClassDef(AliRsnCutKaonForPhi2010PP,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutKaonForPhi2010PP::MatchTOF(AliVTrack *vtrack)
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
