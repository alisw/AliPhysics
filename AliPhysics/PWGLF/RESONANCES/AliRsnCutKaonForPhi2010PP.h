#ifndef ALIRSNCUTKAONFORPHI2010PP_H
#define ALIRSNCUTKAONFORPHI2010PP_H

//
// This cut implements all the checks done to accept a track as a Kaon
// for the pp analysis using 2010 runs.
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutKaonForPhi2010PP : public AliRsnCut {

public:

   AliRsnCutKaonForPhi2010PP(const char *name = "");
   AliRsnCutKaonForPhi2010PP(const AliRsnCutKaonForPhi2010PP &copy);
   AliRsnCutKaonForPhi2010PP &operator=(const AliRsnCutKaonForPhi2010PP &copy);
   virtual ~AliRsnCutKaonForPhi2010PP() { }

   void           SetTPCNSigmaLow (Double_t v) {fNSigmaTPCLow  = v;}
   void           SetTPCNSigmaHigh(Double_t v) {fNSigmaTPCHigh = v;}
   void           SetTPCLimit(Double_t v)      {fLimitTPC      = v;}
   void           SetTOFNSigma(Double_t v)     {fNSigmaTOF     = v;}

   virtual Bool_t IsSelected(TObject *obj);
   void           InitMyPID(Bool_t isMC, Bool_t isESD);

   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}
   AliPIDResponse        *MyPID()      {return fMyPID;}

private:

   Bool_t MatchTOF(const AliVTrack *vtrack);

   Double_t              fNSigmaTPCLow;   // TPC: nsigma cut below limit
   Double_t              fNSigmaTPCHigh;  // TPC: nsigma cut above limit
   Double_t              fLimitTPC;       // TPC: momentum limit
   Double_t              fNSigmaTOF;      // TOF: nsigma cut (unique)

   AliPIDResponse       *fMyPID;          // PID response object to be configured manyally
   AliRsnCutTrackQuality fCutQuality;     // track quality cut

   ClassDef(AliRsnCutKaonForPhi2010PP,1)

};

//__________________________________________________________________________________________________
inline Bool_t AliRsnCutKaonForPhi2010PP::MatchTOF(const AliVTrack *vtrack)
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
