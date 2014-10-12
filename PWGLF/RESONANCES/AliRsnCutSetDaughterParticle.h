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
    kQualityStd2010,//quality only
    kQualityStd2011,
    kQualityStd2011HighPt,    
    kTOFMatch,  
    kFastTPCpidNsigma, //PbPb analysis
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
    kTOFMatchPPB2011, //pA analysis
    kTPCpidKstarPPB2011,
    kTOFpidKstarPPB2011,
    kTPCTOFpidKstarPPB2011,
    kTPCTOFtightPidKStarPPB2011, //TPC 2.0 (3.0) sigma pid & TOF at 3.0 (5.0) sigma 
    kTPCpidMatchPPB2011, //Match with nsigma = fNsigmaTPC
    kTPCpidTOFveto4s, //TPC n sigma + 4.0 sigma TOF veto
    kTPCpidTOFveto3s, //TPC n sigma + 3.0 sigma TOF veto
    kCombinedPidBestPtDep, 
    kTPCPidPtDep,
    kTOFPidPtDep,
    kTPCRejPtDepTOFNsigma,
    kTPCNsigmaTOFVetoPtDep,
    kTPCTOFpidLstar,//cuts for L* in pA
    kNDaughterCuts
  };

  enum ECustomQualityCuts { 
    kDisableCustom = -1,
    kFilterBitCustom,
    kStdLooserDCAXY, 
    kStdLooserDCAZ, 
    kStdCrossedRows60, 
    kStdCrossedRows80, 
    kStdRowsToCls075, 
    kStdRowsToCls085, 
    kStdCls70, 
    kStdChi2TPCCls35,
    kStdUseTPCNcls,
    kNcustomQualityCuts
  };
  
   AliRsnCutSetDaughterParticle();
   AliRsnCutSetDaughterParticle(const char *name,
                                AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID,
                                AliPID::EParticleType pid,
                                Float_t nsigmaFast,
                                Int_t AODfilterBit,
				Bool_t useTPCCrossedRows);
   AliRsnCutSetDaughterParticle(const char *name, 
				AliRsnCutTrackQuality *rsnTrackQualityCut, 
				AliRsnCutSetDaughterParticle::ERsnDaughterCutSet cutSetID, 
				AliPID::EParticleType pid, 
				Float_t nSigmaFast);
   AliRsnCutSetDaughterParticle(const AliRsnCutSetDaughterParticle &copy);
   AliRsnCutSetDaughterParticle &operator=(const AliRsnCutSetDaughterParticle &copy);
   virtual ~AliRsnCutSetDaughterParticle();

   void           Init();
   void           InitStdQualityCuts(Bool_t useTPCCrossedRows=kTRUE);
   void           SetNsigmaForFastTPCpid(Float_t nsigma) {fNsigmaTPC=nsigma; return;};
   void           SetNsigmaForFastTOFpid(Float_t nsigma) {fNsigmaTOF=nsigma; return;};
   void           SetAODTrackCutFilterBit(Int_t ibit) {fAODTrkCutFilterBit=ibit; return;}
   void           SetUseFilterBitOnly(Bool_t useFilterBitOnly=kTRUE) {fCheckOnlyFilterBit=useFilterBitOnly; return;}   
   void           EnableCustomCuts(Bool_t useCustom=kFALSE) {fUseCustomQualityCuts=useCustom; return;}
   void           SetPtRange(Double_t a, Double_t b)        {fPtRange[0] = TMath::Min(a, b); fPtRange[1] = TMath::Max(a, b); return;}
   void           SetEtaRange(Double_t a, Double_t b)       {fEtaRange[0] = TMath::Min(a, b); fEtaRange[1] = TMath::Max(a, b); return;}
   void           SetUse2011StdQualityCuts(Bool_t use2011=kFALSE) {fIsUse2011stdQualityCuts=use2011; return;}
   void           SetUse2011StdQualityCutsHighPt(Bool_t use2011HighPt=kFALSE) {fIsUse2011stdQualityCutsHighPt=use2011HighPt; return;}
   //getters
   const char   *GetAppliedDaughterCutSetName() { return GetName();}
   Int_t         GetAppliedDaughterCutSetId() { return fAppliedCutSetID;}
   const AliRsnCutTrackQuality *GetQualityCut() {return fCutQuality;};
   void          PrintTrackQualityCuts();

 private:
   
   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutSetDaughterParticle::ERsnDaughterCutSet    fAppliedCutSetID;     // ID of applied cut
   Float_t               fNsigmaTPC;         // number of TPC sigmas for fast pid cut only
   Float_t               fNsigmaTOF;         // number of TOF sigmas for fast pid cut only
   AliRsnCutTrackQuality *fCutQuality;       //pointer to quality cut object
   Int_t                 fAODTrkCutFilterBit; //AOD filter bit for track cuts
   Bool_t                fCheckOnlyFilterBit; //flag to use only filter bit cut
   Bool_t                fUseCustomQualityCuts; //flag to enable the usage of custom quality cuts
   Float_t               fPtRange[2]; //single track pt range (min, max)
   Float_t               fEtaRange[2]; //single track eta range (min, max)
   Bool_t                fIsUse2011stdQualityCuts;//flag to enalble std quality cuts 2011 
   Bool_t                fIsUse2011stdQualityCutsHighPt;//flag to enalble std quality cuts 2011 

   ClassDef(AliRsnCutSetDaughterParticle, 5) // cut definitions for K*

};

#endif

