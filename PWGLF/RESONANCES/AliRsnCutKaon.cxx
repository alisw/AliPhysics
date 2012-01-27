//
// All cuts for single pions
// based on track quality and particle identification
// with TPC and TOF.
// Author: fbellini@cern.ch
//
//

#include <Riostream.h>

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliRsnCutKaon.h"

ClassImp(AliRsnCutKaon)

const char *AliRsnCutKaon::fgkCutName[AliRsnCutKaon::kNkaonCuts]= {
   "kQuality_Std2010",
   "kTOFMatch_Std2010",
   "kFastTPCpid_Nsigma",
   "kFastTPCpid_1point5sigma",
   "kFastTPCpid_2sigma",
   "kFastTPCpid_3sigma",
   "kFastTOFpid_Nsigma",
   "kFastTOFpid_1point5sigma",
   "kFastTOFpid_2sigma",
   "kFastTOFpid_3sigma",
   "kTPCTOFpid_DefaultKstarPP2010"
};

//__________________________________________________________________________________________________
AliRsnCutKaon::AliRsnCutKaon(const char *name, AliRsnCutKaon::ERsnKaonCut cutID, AliPID::EParticleType pid) :
   AliRsnCut(name, AliRsnTarget::kDaughter),
   fNoPID(kFALSE),
   fPID(pid),
   fCutQuality(Form("%sQuality", name)),
   fAppliedCutID(cutID),
   fNsigmaCutTPC(1E6),
   fNsigmaCutTOF(1E6)
{
   //
   // Constructor
   // Initialize track quality cuts to 2010 defaults
   //
   fCutQuality.SetPtRange(0.15, 1E+20);
   fCutQuality.SetEtaRange(-0.8, 0.8);
   fCutQuality.SetDCARPtFormula("0.0182+0.0350/pt^1.01");
   fCutQuality.SetDCAZmax(2.0);
   fCutQuality.SetSPDminNClusters(1);
   fCutQuality.SetITSminNClusters(0);
   fCutQuality.SetITSmaxChi2(1E+20);
   fCutQuality.SetTPCminNClusters(70);
   fCutQuality.SetTPCmaxChi2(4.0);
   fCutQuality.SetRejectKinkDaughters();
   fCutQuality.SetAODTestFilterBit(5);
   AliInfo(Form("Applied cut on kaon candidate: %s", AliRsnCutKaon::fgkCutName[fAppliedCutID]));
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutKaon::MatchTOF(const AliVTrack *vtrack)
{
//
// Checks if the track has matched the TOF detector
//
   if (!vtrack) {
      AliWarning("NULL argument: impossible to check status");
      return kFALSE;
   }
   if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
   if (!(vtrack->GetStatus() & AliESDtrack::kTIME)) return kFALSE;

   return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliRsnCutKaon::IsSelected(TObject *obj)
{
//
// Global check
//
   // coherence check
   if (!TargetOK(obj)) return kFALSE;

   // check track
   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;

   Bool_t isSelected=kFALSE;
   switch (fAppliedCutID)
   {
      case    AliRsnCutKaon::kQualityStd2010:
         isSelected=IsSelectedByQualityStd2010(obj, 1, 1, 1);
         break;

      case    AliRsnCutKaon::kTOFMatchStd2010:
         isSelected=IsSelectedByTOFMatchStd2010(obj);
         break;

      case    AliRsnCutKaon::kFastTPCpidNsigma:
         isSelected=IsSelectedByFastTPCpidNsigma(obj,fNsigmaCutTPC);
         break;

      case    AliRsnCutKaon::kFastTPCpid1point5sigma:
         isSelected=IsSelectedByFastTPCpid1point5sigma(obj);
         break;

      case    AliRsnCutKaon::kFastTPCpid2sigma:
         isSelected=IsSelectedByFastTPCpid2sigma(obj);
         break;

      case    AliRsnCutKaon::kFastTPCpid3sigma:
         isSelected=IsSelectedByFastTPCpid3sigma(obj);
         break;

      case    AliRsnCutKaon::kFastTOFpidNsigma:
         isSelected=IsSelectedByFastTOFpidNsigma(obj,fNsigmaCutTOF);
         break;

      case    AliRsnCutKaon::kFastTOFpid1point5sigma:
         isSelected=IsSelectedByFastTOFpid1point5sigma(obj);
         break;

      case    AliRsnCutKaon::kFastTOFpid2sigma:
         isSelected=IsSelectedByFastTOFpid2sigma(obj);
         break;

      case    AliRsnCutKaon::kFastTOFpid3sigma:
         isSelected=IsSelectedByFastTOFpid3sigma(obj);
         break;

      case    AliRsnCutKaon::kTPCTOFpidDefaultKstarPP2010:
         isSelected=IsSelectedByTPCTOFpidDefaultKstarPP2010(obj);
         break;

      default :
         break;
   }

   return isSelected;
}

//-----------------------------------------------------------
Bool_t AliRsnCutKaon::IsSelectedByQualityStd2010(TObject *obj, Bool_t requireTPCin=kTRUE, Bool_t requireTPCrefit=kTRUE,  Bool_t requireITSrefit=kTRUE)
{
   //Checks if track passes the standard quality cuts
   //defined for analysis in 2010 pp and PbPb data

   if (!obj) {
      AliError("Invalid track object passed to function. Please check!");
      return kFALSE;
   }

   fAppliedCutID= AliRsnCutKaon::kQualityStd2010;
   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;

   //optionally check refit flags (to be used for Phi and kStar analysis)
   if (requireTPCin && ((track->GetStatus() & AliESDtrack::kTPCin) == 0) ) return kFALSE;
   if (requireTPCrefit && ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) ) return kFALSE;
   if (requireITSrefit && ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) ) return kFALSE;

   // quality
   if (!fCutQuality.IsSelected(obj)) return kFALSE;
   return kTRUE;
}

//-----------------------------------------------------------
Bool_t AliRsnCutKaon::IsSelectedByTOFMatchStd2010(TObject *obj)
{
   /*
   fbellini@cern.ch, 05.dec.2011
   definition of TOF match: quality cuts 2010 + kTOFout + kTIME
   (optionally) L > 350. cm cut can be added
   */
   if (!obj) {
      AliError("Invalid track object passed to function. Please check!");
      return kFALSE;
   }

   fAppliedCutID=AliRsnCutKaon::kTOFMatchStd2010;

   if (!IsSelectedByQualityStd2010(obj, kTRUE, kTRUE, kTRUE))
      return kFALSE;

   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;
   return MatchTOF(track);
}

//-----------------------------------------------------------
Bool_t AliRsnCutKaon::IsSelectedByFastTPCpidNsigma(TObject *obj,Float_t nSigmaCut=10.)
{
   /*
   fbellini@cern.ch, 05.dec.2011   Double_t pTPC   = track->GetTPCmomentum();
   definition of fast TPC pid N-sigma: quality cuts 2010 + sigma cut on dE/dx without p dependence
   */
   if (!obj) {
      AliError("Invalid track object passed to function. Please check!");
      return kFALSE;
   }
   fAppliedCutID=AliRsnCutKaon::kFastTPCpidNsigma;

   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   if (!IsSelectedByQualityStd2010(obj, kTRUE, kTRUE, kTRUE))
      return kFALSE;

   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;

   Double_t nsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(track,fPID));
   return (nsTPC <= nSigmaCut);
}
//-----------------------------------------------------------
Bool_t AliRsnCutKaon::IsSelectedByFastTOFpidNsigma(TObject *obj,Float_t nSigmaCut=10.)
{
   /*
   fbellini@cern.ch, 05.dec.2011
   definition of fast TOF pid N-sigma: quality cuts 2010 + sigma cut on t-t_0-t_exp without p dependence
   */
   if (!obj) {
      AliError("Invalid track object passed to function. Please check!");
      return kFALSE;
   }

   fAppliedCutID=AliRsnCutKaon::kFastTOFpidNsigma;

   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   if (!IsSelectedByQualityStd2010(obj, kTRUE, kTRUE, kTRUE))
      return kFALSE;

   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;
   Bool_t isTOF = MatchTOF(track);
   Double_t nsTOF  = isTOF ? TMath::Abs(pid->NumberOfSigmasTOF(track, fPID)) : 1E20;

   return (nsTOF <= nSigmaCut);
}

//-----------------------------------------------------------
Bool_t AliRsnCutKaon::IsSelectedByTPCTOFpidDefaultKstarPP2010(TObject *obj)
{
   /*
   fbellini@cern.ch, 05.dec.2011
   definition of TPC-TOF pid for pp phi and Kstar analysis on 2010 data:
   quality cuts 2010 + sigma cut on dE/dx with pTPC dependence + sigma cut on t-t_0-t_exp with p dependence
   */

   if (!obj) {
      AliError("Invalid track object passed to function. Please check!");
      return kFALSE;
   }
   fAppliedCutID=AliRsnCutKaon::kTPCTOFpidDefaultKstarPP2010;

   // check initialization of PID object
   AliPIDResponse *pid = fEvent->GetPIDResponse();
   if (!pid) {
      AliFatal("NULL PID response");
      return kFALSE;
   }

   if (!IsSelectedByQualityStd2010(obj, kTRUE, kTRUE, kTRUE)) //require ITS+TPCrefit and TPCin
      return kFALSE;

   AliVTrack *track = dynamic_cast<AliVTrack *>(fDaughter->GetRef());
   if (!track) return kFALSE;

   Bool_t   isTOF  = MatchTOF(track);
   Double_t pTPC   = track->GetTPCmomentum();
   Double_t p      = track->P();
   Double_t nsTPC  = TMath::Abs(pid->NumberOfSigmasTPC(track, fPID));
   Double_t nsTOF  = isTOF ? TMath::Abs(pid->NumberOfSigmasTOF(track, fPID)) : 1E20;
   Double_t maxTPC = 1E20;
   Double_t maxTOF = 1E20;

   if (isTOF) {
      // TPC: 5sigma cut for all
      if (nsTPC > 5.0) return kFALSE;
      // TOF: 3sigma below 1.5 GeV, 2sigma above
      if (p < 1.5) maxTOF = 3.0; else maxTOF = 2.0;
      return (nsTOF <= maxTOF);
   } else {
      // TPC:
      // all   below   350         MeV: 5sigma
      // all   between 350 and 500 MeV: 3sigma
      // pios above   500         MeV: 2sigma
      // kaons between 500 and 700 MeV: 2sigma
      // kaons above   700         MeV: rejected
      if (pTPC <= 0.35)
         maxTPC = 5.0;
      else if (pTPC > 0.35 && pTPC <= 0.5)
         maxTPC = 3.0;
      else {
         if (pTPC <= 0.7)  maxTPC = 2.0;
         else
            return kFALSE;
      }
      return (nsTPC <= maxTPC);
   }
}
