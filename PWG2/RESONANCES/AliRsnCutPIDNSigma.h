#ifndef ALIRSNCUTPIDNSIGMA_H
#define ALIRSNCUTPIDNSIGMA_H

//
// Class for generalized n-sigma PID cuts with detectors
//

#include "AliPID.h"

#include "AliRsnCut.h"

class AliPIDResponse;

class AliRsnCutPIDNSigma : public AliRsnCut {
public:

   enum EDetector {
      kITS,
      kTPC,
      kTOF,
      kDetectors
   };

   AliRsnCutPIDNSigma(const char *name = "cutPIDNsigma", AliPID::EParticleType species = AliPID::kUnknown, EDetector det = kDetectors, Double_t nsigma = 3.0);
   AliRsnCutPIDNSigma(const AliRsnCutPIDNSigma& copy);
   AliRsnCutPIDNSigma& operator=(const AliRsnCutPIDNSigma& copy);
   virtual ~AliRsnCutPIDNSigma() { }

   void             SetRejectOutside(Bool_t yn = kTRUE)           {fRejectOutside = yn;}
   void             SetRejectUnmatched(Bool_t yn = kTRUE)         {fRejectUnmatched = yn;}
   void             SetMomentumRange(Double_t min, Double_t max)  {fMomMin = min; fMomMax = max;}
   void             SetNSigmaRange(Double_t min, Double_t max)    {AliRsnCut::SetRangeD(min, max);}
   void             SetSpecies(AliPID::EParticleType type)        {fSpecies = type;}
   
   Bool_t           IsITS();
   Bool_t           IsTPC();
   Bool_t           IsTOF();
   
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

private:

   AliPID::EParticleType   fSpecies;         //  particle species
   EDetector               fDetector;        //  detector used for PID
   Double_t                fMomMin;          //  momentum range (for ITS and TOF it is vertex momentum, for TPC it is inner wall)
   Double_t                fMomMax;          //  momentum range (for ITS and TOF it is vertex momentum, for TPC it is inner wall)
   Bool_t                  fRejectOutside;   //  tracks outside momentum range do pass the cut?
   Bool_t                  fRejectUnmatched; //  tracks not matched to this detector do pass the cut?

   ClassDef(AliRsnCutPIDNSigma, 1)
};

inline Bool_t AliRsnCutPIDNSigma::IsITS()
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   AliVTrack *vtrack = fDaughter->GetRefVtrack();
   
   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isITSin  = ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);
   Bool_t isITSpid = ((vtrack->GetStatus() & AliESDtrack::kITSpid) != 0);

   return (isITSin && isITSpid);
}

inline Bool_t AliRsnCutPIDNSigma::IsTPC()
{
//
// Checks if the track has the status flags required for a TPC track
//

   AliVTrack *vtrack = fDaughter->GetRefVtrack();

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   return ((vtrack->GetStatus() & AliESDtrack::kTPCin) != 0);
}

inline Bool_t AliRsnCutPIDNSigma::IsTOF()
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   AliVTrack *vtrack = fDaughter->GetRefVtrack();

   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }

   Bool_t isTOFout = ((vtrack->GetStatus() & AliESDtrack::kTOFout) != 0);
   Bool_t isTIME   = ((vtrack->GetStatus() & AliESDtrack::kTIME) != 0);

   return (isTOFout && isTIME);
}

#endif
