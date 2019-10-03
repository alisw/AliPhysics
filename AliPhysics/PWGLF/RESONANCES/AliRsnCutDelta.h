#ifndef ALIRSNCUTDELTA_H
#define ALIRSNCUTDELTA_H

//
// This cut implements all the checks done to accept a track as a proton and a pion 
// for the pp analysis using 2010 runs. 
// It is based on standard cuts on track quality and nsigma cuts
// with respect to the TPC and TOF signals for the PID.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"
#include <string.h>




class AliRsnCutDelta : public AliRsnCut {

public:

  AliRsnCutDelta(const char *name = "", AliPID::EParticleType pid = AliPID::kPion, Bool_t TPCMethod=kFALSE, Bool_t usePPb=0);
   virtual ~AliRsnCutDelta() { }
   void           SetPID(AliPID::EParticleType type) {fPID = type;}
   void           SetTPCNSigma       (Double_t v)    {fNSigmaTPC        = v;}
   void           SetTPCLimit        (Double_t v)    {fLimitTPC         = v;}
   void           SetTOFNSigma       (Double_t v)    {fNSigmaTOF        = v;}
   void           SetTPCNSigmaProton (Double_t v)    {fNSigmaTPCProton  = v;}
   void           SetTPCNSigmaPion   (Double_t v)    {fNSigmaTPCPion    = v;}
   void           SetTOFNSigmaProton (Double_t v)    {fNSigmaTOFProton  = v;}
   void           SetTOFNSigmaPion   (Double_t v)    {fNSigmaTOFPion    = v;}
   void           SetTPCMomProton    (Double_t v)    {fTPCMomProton     = v;}
   void           SetTOFMomProton    (Double_t v)    {fTOFMomProton     = v;}
   void           SetEtaRange        (Double_t v)    {fEtaRange         = v;}
   void           SetTPCNCluster     (Int_t    v)    {fTPCNCluster      = v;}
   void           SetPtDepDCASigma   (Double_t v)    {fPtDepDCASigma    = v;}
   
   virtual Bool_t IsSelected(TObject *obj);
   AliRsnCutTrackQuality *CutQuality() {return &fCutQuality;}
   Bool_t MatchTOF(const AliVTrack *vtrack);

   private:


   AliPID::EParticleType fPID;            // PID for track
   Double_t              fNSigmaTPC;      // TPC: nsigma cut below limit
   Double_t              fLimitTPC;       // TPC: momentum limit 
   Double_t              fNSigmaTOF;      // TOF: nsigma cut (unique)
   Bool_t 		 fTPCMethod; 	  // Flag for selection of TPC or TOF Methods  kTRUE->TPC  kFALSE->TOF
   Double_t              fNSigmaTPCProton;// TPC: proton nsigma cut 
   Double_t              fNSigmaTPCPion;  // TPC: pion nsigma cut 
   Double_t              fNSigmaTOFProton;// TOF: proton nsigma cut 
   Double_t              fNSigmaTOFPion;  // TOF: pion nsigma cut 
   Double_t              fTPCMomProton;   // TPC: proton mom cut 
   Double_t              fTOFMomProton;   // TOF: proton mom. cut 
   Double_t              fEtaRange;       // Eta range 
   Int_t                 fTPCNCluster;    // Number of TPC clusters 
   Double_t              fPtDepDCASigma;  // Pt dep. DCA n-sigma

   AliRsnCutTrackQuality fCutQuality;     // track quality cut

   ClassDef(AliRsnCutDelta,1)

};

#endif
