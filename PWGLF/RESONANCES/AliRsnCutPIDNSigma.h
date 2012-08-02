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
#include "AliRsnPIDRange.h"

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
   AliRsnCutPIDNSigma(const char *name, AliPID::EParticleType species, EDetector det);
   AliRsnCutPIDNSigma(const AliRsnCutPIDNSigma &copy);
   AliRsnCutPIDNSigma &operator=(const AliRsnCutPIDNSigma &copy);
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

#endif
