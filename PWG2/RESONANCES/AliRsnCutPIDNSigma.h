#ifndef ALIRSNCUTPIDNSIGMA_H
#define ALIRSNCUTPIDNSIGMA_H

//
// Class for n-sigma PID cuts.
// ---
// Requires:
//
// 1) the used detector, chosen from an enumeration
// 2) the reference charged particle species, chosen from AliPID enumeration
// 3) a momentum range: outside it, the cut is never passed
//

#include <TMath.h>

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

   AliRsnCutPIDNSigma();
   AliRsnCutPIDNSigma(const char *name, AliPID::EParticleType species, EDetector det, Double_t nsigma);
   AliRsnCutPIDNSigma(const AliRsnCutPIDNSigma& copy);
   AliRsnCutPIDNSigma& operator=(const AliRsnCutPIDNSigma& copy);
   virtual ~AliRsnCutPIDNSigma() { }

   void             SetSpecies(AliPID::EParticleType type)        {fSpecies = type;}
   void             SetDetector(EDetector det)                    {fDetector = det;}
   void             SetRejectUnmatched(Bool_t yn = kTRUE)         {fRejectUnmatched = yn;}
   void             SetNSigma(Double_t nsigma)                    {fNSigma = nsigma;}
   void             SetMomentumRange(Double_t min, Double_t max);
   
   Bool_t           IsITS(const AliVTrack *vtrack) const;
   Bool_t           IsTPC(const AliVTrack *vtrack) const;
   Bool_t           IsTOF(const AliVTrack *vtrack) const;
   
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

private:

   AliPID::EParticleType   fSpecies;         //  particle species
   EDetector               fDetector;        //  detector used for PID
   Bool_t                  fRejectUnmatched; //  tracks not matched to this detector do pass the cut?
   Double_t                fMomMin;          //  momentum range (ITS/TOF: vertex momentum, TPC: mom at inner wall)
   Double_t                fMomMax;          //  momentum range (ITS/TOF: vertex momentum, TPC: mom at inner wall)
   Double_t                fNSigma;          //  nsigma range

   ClassDef(AliRsnCutPIDNSigma, 1)
};

inline void AliRsnCutPIDNSigma::SetMomentumRange(Double_t min, Double_t max)
{
//
// Assigns the range in total momentum used for check
// For ITS and TOF, it is used to check momentum at vertex,
// for TPC it is used to check momentum at TPC inner barrel
//

   fMomMin = TMath::Min(min, max);
   fMomMax = TMath::Max(min, max);
}

inline Bool_t AliRsnCutPIDNSigma::IsITS(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kITSin)  == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kITSpid) == 0) return kFALSE;

   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::IsTPC(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for a TPC track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTPCin) == 0) return kFALSE;
   
   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::IsTOF(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTOFout) == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kTIME)   == 0) return kFALSE;

   return kTRUE;
}

#endif
