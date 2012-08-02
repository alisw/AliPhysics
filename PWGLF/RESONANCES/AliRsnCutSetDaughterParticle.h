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
     kNDaughterCuts
  };
     
   AliRsnCutSetDaughterParticle();
   AliRsnCutSetDaughterParticle(const char *name, 
				AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID,
				AliPID::EParticleType pid,
				Float_t nsigmaFast);
   AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy);
   AliRsnCutSetDaughterParticle &operator=(const AliRsnCutSetDaughterParticle &copy);
   virtual ~AliRsnCutSetDaughterParticle() { }
   
   void           Init();
   void           SetNsigmaForFastTPCpid(Float_t nsigma){fNsigmaTPC=nsigma; return;};
   void           SetNsigmaForFastTOFpid(Float_t nsigma){fNsigmaTOF=nsigma; return;};
   //getters
   const char *  GetAppliedDaughterCutSetName() { return GetName();}
   const Int_t   GetAppliedDaughterCutSetId() { return fAppliedCutSetID;}

private:

   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet    fAppliedCutSetID;     // ID of applied cut
   Float_t               fNsigmaTPC;         // number of TPC sigmas for fast pid cut only
   Float_t               fNsigmaTOF;         // number of TOF sigmas for fast pid cut only 
   ClassDef(AliRsnCutSetDaughterParticle, 1) // cut definitions for K*

};

#endif
