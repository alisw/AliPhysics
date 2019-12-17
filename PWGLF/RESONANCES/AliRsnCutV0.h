#ifndef ALIRSNCUTV0_H
#define ALIRSNCUTV0_H

#include <TMath.h>
#include <TString.h>

#include "AliRsnCut.h"
#include "AliPIDResponse.h"

#include "AliRsnCutTrackQuality.h"

class AliESDtrack;
class AliAODTrack;


class AliRsnCutV0 : public AliRsnCut {

public:


   AliRsnCutV0(const char *name = "AliRsnCutV0", Int_t hypothesis = kLambda0, AliPID::EParticleType pid = AliPID::kProton, AliPID::EParticleType pid2 = AliPID::kPion);
   AliRsnCutV0(const AliRsnCutV0 &copy);
   AliRsnCutV0 &operator=(const AliRsnCutV0 &copy);
   virtual ~AliRsnCutV0() { }

   void           SetESDtrackCuts(AliESDtrackCuts *cuts)   {fESDtrackCuts = cuts;}
   void           SetHypothesis(Int_t code);
   void           SetMassTolSigma(Double_t value)          {fMassTolSigma = value;}     
   void           SetpT_Tolerance(Int_t value)          	  {fpT_Tolerance = value;} 
   void           SetTolerance(Double_t value)             {fTolerance = value;}
   void           SetToleranceVeto(Double_t value)         {fToleranceVeto = value;}
   void           SetSwitch(Bool_t value)                  {fSwitch = value;}
   void           SetfLife(Double_t value)                 {fLife = value;}
   void           SetfLowRadius(Double_t value)            {fLowRadius = value;}
   void           SetfHighRadius(Double_t value)           {fHighRadius = value;}    
   void           SetMaxDCAVertex(Double_t value)          {fMaxDCAVertex = value;}
   void           SetMinCosPointingAngle(Double_t value)   {fMinCosPointAngle = value;}
   void           SetMaxDaughtersDCA(Double_t value)       {fMaxDaughtersDCA = value;}
   void           SetMinTPCcluster(Int_t value)            {fMinTPCcluster = value;}
   void           SetMaxRapidity(Double_t value)           {fMaxRapidity = value;}
   void           SetMaxPseudorapidity(Double_t value)           {fMaxPseudorapidity = value;}
   
   void           SetPIDCutProton(Double_t value)          {fPIDCutProton = value;}
   void           SetPIDCutPion(Double_t value)            {fPIDCutPion = value;}
   
   void           SetDifferentDCACutPosNegTrack(Bool_t doDifferentTrackDCACuts){fCustomTrackDCACuts = doDifferentTrackDCACuts;}
   void           SetMinDCAToVtxXYPositiveTrack(Double_t value) {fMinDCAPositiveTrack = value;}
   void           SetMinDCAToVtxXYNegativeTrack(Double_t value) {fMinDCANegativeTrack = value;}
   void           SetCheckOOBPileup(Bool_t value = true)   {fCheckOOBPileup = value;}

   AliRsnCutTrackQuality *CutQuality()                     {return &fCutQuality;}
   void           SetAODTestFilterBit(Int_t value)         {fAODTestFilterBit = value;}
   Int_t          GetAODTestFilterBit()                    {return fAODTestFilterBit;}

   virtual Bool_t IsSelected(TObject *obj);
   virtual void   Print(const Option_t *option = "") const;

protected:

   Bool_t      CheckESD(AliESDv0 *track);
   Bool_t      CheckAOD(AliAODv0 *track);
   Bool_t      TrackPassesOOBPileupCut(AliESDtrack* t, Double_t b);
   
   Int_t            fHypothesis;       // PDG code corresponding to expected V0 hypothesis
   Int_t            fpT_Tolerance=0;     // Switch to set pT dependent Mass Tolerance
   Double_t         fMassTolSigma;      //Sigma cut for pt Dependent Mass Tol Cut    
   Double_t         fMass;             // mass corresponding to hypothesis
   Double_t         fTolerance;        // tolerance in the difference between computed and expected mass
   Double_t         fToleranceVeto;    // Competing V0 Rejection. Read the note in AliRsnCutV0.cxx for more info.
   Bool_t           fSwitch;           // Switch for using Competing V0 Rejection
   Double_t         fLife;             // Lifetime for positive track
   Double_t         fLowRadius;        // Lower Limit on Fiducial Volume
   Double_t         fHighRadius;       // Higher Limit on Fiducial Volume
   Double_t         fMaxDCAVertex;     // max allowed DCA from primary vertex
   Double_t         fMinCosPointAngle; // min allowed cosine of pointing angle
   Double_t         fMaxDaughtersDCA;  // max allowed DCA between the two daughers
   Int_t            fMinTPCcluster;    // min allowed TOC cluster
   Double_t         fMaxRapidity;      // max allowed V0 rapidity
   Double_t         fMaxPseudorapidity; // max allowed V0 pseudorapidity
   Bool_t           fCustomTrackDCACuts; // Use different DCA cuts for positive and negative V0 tracks
   Double_t         fMinDCAPositiveTrack; // DCA of positive V0 track to vertex
   Double_t         fMinDCANegativeTrack; // DCA of negative V0 track to vertex
   Bool_t           fCheckOOBPileup;   // Check out-of-bunch pileup
   
   AliPID::EParticleType fPID;         // PID for track
   AliPID::EParticleType fPID2;        // PID for track

   Double_t         fPIDCutProton;        // nsigmas for protons
   Double_t         fPIDCutPion;          // nsigmas for pions
   
   AliESDtrackCuts *fESDtrackCuts;     // quality cuts for v0 daughters
   
   AliRsnCutTrackQuality fCutQuality;  // track quality cut
   
   Int_t            fAODTestFilterBit; // test filter bit for AODs
   
   ClassDef(AliRsnCutV0, 2)
};

//__________________________________________________________________________________________________
inline void AliRsnCutV0::SetHypothesis(Int_t code)
{
//
// Assign a V0 species hypothesis, which also assign the expected mass
//

   fHypothesis = code;

   switch (fHypothesis) {
      case kLambda0:
         fMass = 1.11568;
         break;
      case kLambda0Bar:
         fMass = 1.11568;
         break;
      case kK0Short:
         fMass = 0.497614;
         break;
      default:
         AliError(Form("You are setting an unexpected hypothesis: %d", fHypothesis));
         fMass = 1E20;
   }
}

#endif
