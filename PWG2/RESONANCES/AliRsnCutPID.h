//
// Class AliRsnCutRange
//
// General implementation of cuts which check a value inside a range.
// This range can be defined by two integers or two doubles.
// A user-friendly enumeration allows to define what is checked.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#ifndef ALIRSNCUTPID_H
#define ALIRSNCUTPID_H

#include "AliPID.h"
#include "AliRsnCut.h"

class AliRsnDaughter;

class AliRsnCutPID : public AliRsnCut {
public:

   enum EDetector {
      kITS,
      kTPC,
      kTRD,
      kTOF,
      kHMPID,
      kDetectors
   };

   AliRsnCutPID();
   AliRsnCutPID(const char *name, AliPID::EParticleType pid, Double_t probMin = 0.0, Bool_t perfectPID = kFALSE);
   virtual ~AliRsnCutPID() {;};

   void           SetDefault(Bool_t yn = kTRUE) {fUseDefault = yn;}
   void           SetPrior(Int_t i, Double_t value) {if (i >= 0 && i < AliPID::kSPECIES) fPrior[i] = value;}
   void           SetPerfect(Bool_t yn = kTRUE) {fPerfect = yn;}

   void           IncludeDetector(EDetector det, Double_t threshold = 0., Bool_t goAbove = kTRUE);
   void           ExcludeDetector(EDetector det) {if (CheckBounds(det)) fUseDetector[det] = kFALSE;}

   Bool_t         ComputeWeights(AliRsnDaughter *daughter);
   Int_t          RealisticPID(AliRsnDaughter * const daughter, Double_t &prob);
   Int_t          PerfectPID(AliRsnDaughter * const daughter);
   Double_t       GetWeight(Int_t i) const {if (i >= 0 && i < AliPID::kSPECIES) return fWeight[i]; return 0.0;}

   virtual Bool_t IsSelected(TObject *object);

protected:

   Bool_t   CheckBounds(EDetector det) const {return (det >= kITS && det < kDetectors);}
   Bool_t   CheckThreshold(EDetector det, Double_t value) const;

   Double_t              fPrior[AliPID::kSPECIES];        // prior probability
   Double_t              fWeight[AliPID::kSPECIES];       // PID weights used for combinations

   Bool_t                fPerfect;                        // choice to use perfect PID
   Bool_t                fUseDefault;                     // choice to use default combined PID weights (or customized)
   Bool_t                fUseDetector[kDetectors];        // flag to include/exclude each single detector
   Double_t              fPtThreshold[kDetectors];        // pT threshold above/below which a detector is considered
   Double_t              fGoAboveThreshold[kDetectors];   // to choose if detector is used balow or above threshold

   ClassDef(AliRsnCutPID, 1)
};

#endif
