#ifndef ALIRSNCUTSETDAUGHTERPARTICLE_H
#define ALIRSNCUTSETDAUGHTERPARTICLE_H

//
// Cuts collection for selecting good daughter candidates for rsn analysis
//Requires:
// 1) choice of existing cuts among the enum list
// 2) PID ipothesis for the daughter particle
//
// Author: Francesca Bellini (fbellini@cern.ch)
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutSet.h"
#include "AliRsnCutTrackQuality.h"
#include "AliRsnCutPIDNSigma.h"
#include "AliRsnCutTOFMatch.h"
#include "AliRsnCutPhi.h"

class AliRsnCutSetDaughterParticle : public AliRsnCutSet {

public:

   enum ERsnDaughterCutSet {
      kNoCuts,
      kQualityStd2010,
      kQualityStd2011,    
      kTOFMatch,
      kFastTPCpidNsigma,
      kFastTOFpidNsigma,
      kTPCTOFpidKstarPP2010,
      kTOFpidKstarPbPb2010,
      kTOFTPCmismatchKstarPbPb2010,
      kTOFMatchTRD2010,
      kTOFMatchNoTRD2010,
      kTOFpidKstarPbPbTRD2010,
      kTOFpidKstarPbPbNoTRD2010,
      kTOFMatchTPCpidNsigma,     
      kQualityStd2010TRD,
      kQualityStd2010NoTRD,
      kNDaughterCuts
   };

   AliRsnCutSetDaughterParticle();
   AliRsnCutSetDaughterParticle(const char *name,
                                AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID,
                                AliPID::EParticleType pid,
                                Float_t nsigmaFast,
                                Int_t AODfilterBit);
   AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy);
   AliRsnCutSetDaughterParticle &operator=(const AliRsnCutSetDaughterParticle &copy);
   virtual ~AliRsnCutSetDaughterParticle();

   void           Init();
   void           SetNsigmaForFastTPCpid(Float_t nsigma) {fNsigmaTPC=nsigma; return;};
   void           SetNsigmaForFastTOFpid(Float_t nsigma) {fNsigmaTOF=nsigma; return;};
   void           SetAODTrackCutFilterBit(Int_t ibit) {fAODTrkCutFilterBit=ibit; return;}
   //getters
   const char   *GetAppliedDaughterCutSetName() { return GetName();}
   Int_t   GetAppliedDaughterCutSetId() { return fAppliedCutSetID;}
   const AliRsnCutTrackQuality *GetQualityCut() {return fCutQuality;};

private:

   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet    fAppliedCutSetID;     // ID of applied cut
   Float_t               fNsigmaTPC;         // number of TPC sigmas for fast pid cut only
   Float_t               fNsigmaTOF;         // number of TOF sigmas for fast pid cut only
   AliRsnCutTrackQuality *fCutQuality;       //pointer to quality cut object
   Int_t                 fAODTrkCutFilterBit; //AOD filter bit for track cuts
   ClassDef(AliRsnCutSetDaughterParticle, 2) // cut definitions for K*

};

#endif
