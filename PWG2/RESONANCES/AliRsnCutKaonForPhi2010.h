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

   enum ECutMode {
      kQuality = 0,
      kOnlyTPC,
      kOnlyTOF,
      kDefaultPID,
      kModes
   };

   AliRsnCutKaonForPhi2010(const char *name = "", Double_t nSigmaTPC = 3.0, Double_t nSigmaTOF = 3.0, Double_t tofLimit = 0.8);
   AliRsnCutKaonForPhi2010(const AliRsnCutKaonForPhi2010 &copy);
   AliRsnCutKaonForPhi2010& operator=(const AliRsnCutKaonForPhi2010 &copy);
   virtual ~AliRsnCutKaonForPhi2010() { }
   
   virtual Bool_t IsSelected(TObject *obj);
   
   AliRsnCutTrackQuality *CutQuality()            {return &fCutQuality;}
   Double_t               GetTOFthreshold() const {return fTOFthreshold;}
   AliPIDResponse        *MyPID()                 {return fMyPID;}
   
   void   SetMode(ECutMode mode)            {fMode = mode;}
   void   SetCutTPC(Double_t cut)           {fCutTPC = cut;}
   void   SetCutTOF(Double_t cut)           {fCutTOF = cut;}
   void   SetTOFthreshold(Double_t value)   {fTOFthreshold = value;}
   
   void   InitMyPID(Bool_t isMC, Bool_t isESD);
   Bool_t MatchTOF(const AliVTrack *vtrack);

private:

   ECutMode              fMode;          // how the cut is applied
   Double_t              fCutTPC;        // nsigma cut for TPC
   Double_t              fCutTOF;        // nsigma cut fof TOF
   Double_t              fTOFthreshold;  // for Pt above this threshold, TOF is mandatory
   AliPIDResponse       *fMyPID;         // PID response object to be configured manyally
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
