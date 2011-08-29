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
#include <TClonesArray.h>

#include "AliPID.h"
#include "AliESDtrack.h"

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
   
   //
   // This allows to define several intervals
   //
   class AliRsnPIDRange : public TObject {
   public:
      
      AliRsnPIDRange(Double_t nsigma, Double_t pmin, Double_t pmax)
         : fPMin(pmin), fPMax(pmax), fNSigmaCut(nsigma) { }
         
      Double_t& PMin()      {return fPMin;}
      Double_t& PMax()      {return fPMax;}
      Double_t& NSigmaCut() {return fNSigmaCut;}
      
      Bool_t IsInRange(Double_t mom)  {return (mom >= fPMin && mom <= fPMax);}
      Bool_t CutPass(Double_t nsigma) {return (nsigma <= fNSigmaCut);}
      
   private:
   
      Double_t fPMin;      // lower bound of momentum range
      Double_t fPMax;      // upper bound of momentum range
      Double_t fNSigmaCut; // cut in number of sigmas
      
      ClassDef(AliRsnCutPIDNSigma::AliRsnPIDRange,1)
   };

   AliRsnCutPIDNSigma();
   AliRsnCutPIDNSigma(const char *name, AliPID::EParticleType species, EDetector det);
   AliRsnCutPIDNSigma(const AliRsnCutPIDNSigma& copy);
   AliRsnCutPIDNSigma& operator=(const AliRsnCutPIDNSigma& copy);
   virtual ~AliRsnCutPIDNSigma() { }

   void             SetSpecies(AliPID::EParticleType type)        {fSpecies = type;}
   void             SetDetector(EDetector det)                    {fDetector = det;}
   void             SetRejectUnmatched(Bool_t yn = kTRUE)         {fRejectUnmatched = yn;}
   
   AliPIDResponse  *MyPID()                                       {return fMyPID;}
   void             InitMyPID(Bool_t isMC, Bool_t isESD);
   
   void             SinglePIDRange(Double_t nsigma);
   void             AddPIDRange(Double_t nsigma, Double_t pmin = 0.0, Double_t pmax = 1E20);
   
   Bool_t           MatchITS(const AliVTrack *vtrack) const;
   Bool_t           MatchTPC(const AliVTrack *vtrack) const;
   Bool_t           MatchTOF(const AliVTrack *vtrack) const;
   Bool_t           MatchDetector(const AliVTrack *vtrack) const;
   
   virtual Bool_t   IsSelected(TObject *object);
   virtual void     Print(const Option_t *option = "") const;

private:

   AliPID::EParticleType   fSpecies;         //  particle species
   EDetector               fDetector;        //  detector used for PID
   Bool_t                  fRejectUnmatched; //  tracks not matched to this detector do pass the cut?
   Double_t                fTrackNSigma;     //! tmp track number of sigmas w.r. to chosen detector
   Double_t                fTrackMom;        //! track reference momentum
   AliPIDResponse         *fMyPID;           //  PID response object to be configured manyally
   TClonesArray            fRanges;          //  collection of ranges

   ClassDef(AliRsnCutPIDNSigma, 1)
};

inline Bool_t AliRsnCutPIDNSigma::MatchITS(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kITSin)  == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kITSpid) == 0) return kFALSE;

   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::MatchTPC(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for a TPC track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTPCin) == 0) return kFALSE;
   
   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::MatchTOF(const AliVTrack *vtrack) const
{
//
// Checks if the track has the status flags required for an ITS standalone track
//

   if ((vtrack->GetStatus() & AliESDtrack::kTOFout) == 0) return kFALSE;
   if ((vtrack->GetStatus() & AliESDtrack::kTIME)   == 0) return kFALSE;

   return kTRUE;
}

inline Bool_t AliRsnCutPIDNSigma::MatchDetector(const AliVTrack *vtrack) const
{
//
// Checks if the track has matched the required detector.
// If no valid detector is specified, kFALSE is always returned.
//

   switch (fDetector) {
      case kITS: return MatchITS(vtrack);
      case kTPC: return MatchTPC(vtrack);
      case kTOF: return MatchTOF(vtrack);
      default  : return kFALSE;
   }
}

inline void AliRsnCutPIDNSigma::AddPIDRange(Double_t nsigma, Double_t pmin, Double_t pmax)
{
//
// Add a new slot for checking PID
//

   Int_t n = fRanges.GetEntries();
   
   new (fRanges[n]) AliRsnPIDRange(nsigma, pmin, pmax);
}

inline void AliRsnCutPIDNSigma::SinglePIDRange(Double_t nsigma)
{
//
// Clear all slots and sets a unique one
//

   fRanges.Delete();
   
   new (fRanges[0]) AliRsnPIDRange(nsigma, 0.0, 1E20);
}

#endif
