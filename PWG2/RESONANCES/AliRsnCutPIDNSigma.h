#ifndef ALIRSNCUTPIDNSIGMA_H
#define ALIRSNCUTPIDNSIGMA_H

//
// Class for generalized n-sigma PID cuts with detectors
//

#include "AliPID.h"

#include "AliRsnCut.h"

class AliVTrack;
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

   void             SetRejectUnmatched(Bool_t yn = kTRUE)         {fRejectUnmatched = yn;}
   void             SetMomentumRange(Double_t min, Double_t max)  {fMomMin = min; fMomMax = max;}
   void             SetNSigmaRange(Double_t min, Double_t max)    {AliRsnCut::SetRangeD(min, max);}
   void             SetSpecies(AliPID::EParticleType type)        {fSpecies = type;}
   
   Bool_t           IsITS(AliVTrack *vtrack);
   Bool_t           IsTPC(AliVTrack *vtrack);
   Bool_t           IsTOF(AliVTrack *vtrack);
   
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

private:

   AliPID::EParticleType   fSpecies;         //  particle species
   EDetector               fDetector;        //  detector used for PID
   Double_t                fMomMin;          //  momentum range (for ITS and TOF it is vertex momentum, for TPC it is inner wall)
   Double_t                fMomMax;          //  momentum range (for ITS and TOF it is vertex momentum, for TPC it is inner wall)
   Bool_t                  fRejectUnmatched; //  tracks not matched to this detector do pass the cut?

   ClassDef(AliRsnCutPIDNSigma, 1)
};

inline Bool_t AliRsnCutPIDNSigma::IsITS(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kITSin)  == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kITSpid) == 0) return kFALSE;

   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::IsTPC(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for a TPC track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTPCin) == 0) return kFALSE;
   
   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::IsTOF(AliVTrack *vtrack)
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTOFout) == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kTIME)   == 0) return kFALSE;

   return kTRUE;
}

#endif
