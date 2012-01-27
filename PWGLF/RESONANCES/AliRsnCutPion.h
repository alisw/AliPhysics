#ifndef ALIRSNCUTPION_H
#define ALIRSNCUTPION_H

//
// Cuts for selecting good pion candidates for rsn analysis
// Applies track quality selection plus PID selection,
// with different tolerance ranges depending on the momentum.
//

#include "AliVTrack.h"
#include "AliRsnCut.h"
#include "AliRsnCutTrackQuality.h"

class AliRsnCutPion : public AliRsnCut {

public:
   enum ERsnPionCut {
      kQualityStd2010,
      kTOFMatchStd2010,
      kFastTPCpidNsigma,
      kFastTPCpid1point5sigma,
      kFastTPCpid2sigma,
      kFastTPCpid3sigma,
      kFastTOFpidNsigma,
      kFastTOFpid1point5sigma,
      kFastTOFpid2sigma,
      kFastTOFpid3sigma,
      kTPCTOFpidDefaultKstarPP2010,
      kNpionCuts
   };

   AliRsnCutPion(const char *name = "",  AliRsnCutPion::ERsnPionCut cutID = AliRsnCutPion::kQualityStd2010, AliPID::EParticleType pid = AliPID::kPion);
   virtual ~AliRsnCutPion() { }

   void                   SetNoPID(Bool_t yn = kTRUE)        {fNoPID = yn;}
   AliRsnCutTrackQuality *CutQuality()                       {return &fCutQuality;}
   Bool_t                 MatchTOF(const AliVTrack *vtrack);

   Bool_t         IsSelected(TObject *obj);

   Bool_t         IsSelectedByQualityStd2010(TObject *obj, Bool_t requireTPCin, Bool_t requireTPCrefit,  Bool_t requireITSrefit );
   Bool_t         IsSelectedByTOFMatchStd2010(TObject *obj);

   Bool_t         IsSelectedByFastTPCpidNsigma(TObject *obj,Float_t nSigmaCut);
   Bool_t         IsSelectedByFastTPCpid1point5sigma(TObject *obj)   {return IsSelectedByFastTPCpidNsigma(obj, 1.5);}
   Bool_t         IsSelectedByFastTPCpid2sigma(TObject *obj)         {return IsSelectedByFastTPCpidNsigma(obj, 2.0);}
   Bool_t         IsSelectedByFastTPCpid3sigma(TObject *obj)         {return IsSelectedByFastTPCpidNsigma(obj, 3.0);}

   Bool_t         IsSelectedByFastTOFpidNsigma(TObject *obj,Float_t nSigmaCut);
   Bool_t         IsSelectedByFastTOFpid1point5sigma(TObject *obj)   {return IsSelectedByFastTOFpidNsigma(obj, 1.5);}
   Bool_t         IsSelectedByFastTOFpid2sigma(TObject *obj)         {return IsSelectedByFastTOFpidNsigma(obj, 2.0);}
   Bool_t         IsSelectedByFastTOFpid3sigma(TObject *obj)         {return IsSelectedByFastTOFpidNsigma(obj, 3.0);}

   Bool_t         IsSelectedByTPCTOFpidDefaultKstarPP2010(TObject *obj);

   //setters
   void     SetNsigmaCutTPC(Float_t nMaxSigmaCut=1E6) {fNsigmaCutTPC=nMaxSigmaCut; return;}
   void     SetNsigmaCutTOF(Float_t nMaxSigmaCut=1E6) {fNsigmaCutTOF=nMaxSigmaCut; return;}

   //getters
   const char   *GetAppliedPionCutName() {if (fAppliedCutID>0) { return fgkCutName[fAppliedCutID];} else {return "none";}}
   Int_t   GetAppliedPionCutId() { return fAppliedCutID;}
   const char   *GetPionCutName(AliRsnCutPion::ERsnPionCut cutID) { return fgkCutName[cutID];}

private:

   Bool_t                fNoPID;            // flag to switch off PID check
   AliPID::EParticleType fPID;              // PID for track
   AliRsnCutTrackQuality fCutQuality;       // track quality cut
   Int_t                 fAppliedCutID;     // ID of applied cut
   static const char    *fgkCutName[kNpionCuts]; //array with cuts names
   Float_t               fNsigmaCutTPC;     //  max sigma for TPC dE/dx fast cut
   Float_t               fNsigmaCutTOF;     //  max sigma for TOF t-t0-t_exp fast cut

   ClassDef(AliRsnCutPion,1) // cut definitions for K*

};

#endif
