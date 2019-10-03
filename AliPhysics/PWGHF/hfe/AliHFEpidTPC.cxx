/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//
// Class for TPC PID
// Implements the abstract base class AliHFEpidBase
// 
// Class contains TPC specific cuts and QA histograms
// Two selection strategies are offered: Selection of certain value
// regions in the TPC dE/dx (by IsSelected), and likelihoods
//
// Authors: 
//
//   Markus Fasel <M.Fasel@gsi.de> 
//   Markus Heide <mheide@uni-muenster.de> 
//  
#include <TF1.h>
#include <TMath.h>
#include <TH2D.h>

#include "AliTPCdEdxInfo.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliPIDResponse.h"

#include "AliHFEpidTPC.h"
#include "AliHFEpidQAmanager.h"

ClassImp(AliHFEpidTPC)

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC() :
// add a list here
AliHFEpidBase()
, fLineCrossingsEnabled(0)
, fkEtaCorrection(NULL)
, fkCentralityCorrection(NULL)
, fkEtaMeanCorrection(NULL)
, fkEtaWidthCorrection(NULL)
, fkPMeanCorrection(NULL)
, fkPWidthCorrection(NULL)
, fkCentralityMeanCorrection(NULL)
, fkCentralityWidthCorrection(NULL)
, fkCentralityEtaCorrectionMeanJpsi(NULL)
, fkCentralityEtaCorrectionWidthJpsi(NULL)
, fHasCutModel(kFALSE)
, fUseOnlyOROC(kFALSE)
, fNsigmaTPC(3)
, fRejectionEnabled(0)
, fUsedEdx(kFALSE)
{
   //
   // default  constructor
   //

   memset(fkUpperSigmaCut, 0, sizeof(const TF1 *) * 12);
   memset(fkLowerSigmaCut, 0, sizeof(const TF1 *) * 12);

   memset(fRejection, 0, sizeof(Float_t) * 4 * AliPID::kSPECIES);
   memset(fLineCrossingSigma, 0, sizeof(Double_t) * AliPID::kSPECIES);
   memset(fPAsigCut, 0, sizeof(Float_t) * 2);
   memset(fNAsigmaTPC, 0, sizeof(Float_t) * 2);

}

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const char* name) :
// add a list here
AliHFEpidBase(name)
, fLineCrossingsEnabled(0)
, fkEtaCorrection(NULL)
, fkCentralityCorrection(NULL)
, fkEtaMeanCorrection(NULL)
, fkEtaWidthCorrection(NULL)
, fkPMeanCorrection(NULL)
, fkPWidthCorrection(NULL)
, fkCentralityMeanCorrection(NULL)
, fkCentralityWidthCorrection(NULL)
, fkCentralityEtaCorrectionMeanJpsi(NULL)
, fkCentralityEtaCorrectionWidthJpsi(NULL)
, fHasCutModel(kFALSE)
, fUseOnlyOROC(kFALSE)
, fNsigmaTPC(3)
, fRejectionEnabled(0)
, fUsedEdx(kFALSE)
{
   //
   // default  constructor
   //
   //
   memset(fkUpperSigmaCut, 0, sizeof(const TF1 *) * 12);
   memset(fkLowerSigmaCut, 0, sizeof(const TF1 *) * 12);

   memset(fRejection, 0, sizeof(Float_t) * 4 * AliPID::kSPECIES);
   memset(fLineCrossingSigma, 0, sizeof(Double_t) * AliPID::kSPECIES);
   memset(fPAsigCut, 0, sizeof(Float_t) * 2);
   memset(fNAsigmaTPC, 0, sizeof(Float_t) * 2);
}

//___________________________________________________________________
AliHFEpidTPC::AliHFEpidTPC(const AliHFEpidTPC &ref) :
AliHFEpidBase("")
, fLineCrossingsEnabled(0)
, fkEtaCorrection(NULL)
, fkCentralityCorrection(NULL)
, fkEtaMeanCorrection(NULL)
, fkEtaWidthCorrection(NULL)
, fkCentralityMeanCorrection(NULL)
, fkCentralityWidthCorrection(NULL)
, fkCentralityEtaCorrectionMeanJpsi(NULL)
, fkCentralityEtaCorrectionWidthJpsi(NULL)
, fHasCutModel(ref.fHasCutModel)
, fUseOnlyOROC(ref.fUseOnlyOROC)
, fNsigmaTPC(2)
, fRejectionEnabled(0)
, fUsedEdx(kFALSE)
{
   //
   // Copy constructor
   //
   ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTPC &AliHFEpidTPC::operator=(const AliHFEpidTPC &ref){
   //
   // Assignment operator
   //
   if(this != &ref){
      ref.Copy(*this);
   }
   return *this;
}
//___________________________________________________________________
void AliHFEpidTPC::Copy(TObject &o) const{
   //
   // Copy function
   // called in copy constructor and assigment operator
   //
   AliHFEpidTPC &target = dynamic_cast<AliHFEpidTPC &>(o);

   target.fkEtaCorrection = fkEtaCorrection;
   target.fkCentralityCorrection = fkCentralityCorrection;
   target.fkEtaMeanCorrection = fkEtaMeanCorrection;
   target.fkEtaWidthCorrection = fkEtaWidthCorrection;
   target.fkPMeanCorrection = fkPMeanCorrection;
   target.fkPWidthCorrection = fkPWidthCorrection;
   target.fkCentralityMeanCorrection = fkCentralityMeanCorrection;
   target.fkCentralityWidthCorrection = fkCentralityWidthCorrection;
   target.fLineCrossingsEnabled = fLineCrossingsEnabled;
   target.fHasCutModel = fHasCutModel;
   target.fUseOnlyOROC = fUseOnlyOROC;
   target.fNsigmaTPC = fNsigmaTPC;
   target.fRejectionEnabled = fRejectionEnabled;

   memcpy(target.fkUpperSigmaCut, fkUpperSigmaCut, sizeof(const TF1 *) * 12);
   memcpy(target.fkLowerSigmaCut, fkLowerSigmaCut, sizeof(const TF1 *) * 12);

   memcpy(target.fLineCrossingSigma, fLineCrossingSigma, sizeof(Double_t) * AliPID::kSPECIES);
   memcpy(target.fPAsigCut, fPAsigCut, sizeof(Float_t) * 2);
   memcpy(target.fNAsigmaTPC, fNAsigmaTPC, sizeof(Float_t) * 2);

   AliHFEpidBase::Copy(target);
}

//___________________________________________________________________
AliHFEpidTPC::~AliHFEpidTPC(){
   //
   // Destructor
   //
}

//___________________________________________________________________
Bool_t AliHFEpidTPC::InitializePID(Int_t /*run*/){
   //
   // Add TPC dE/dx Line crossings
   //
   //AddTPCdEdxLineCrossing(AliPID::kKaon, 0.3, 0.018);
   //AddTPCdEdxLineCrossing(AliPID::kProton, 0.9, 0.054);
   return kTRUE;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const
{
   //
   // For the TPC pid we use the 2-sigma band around the bethe bloch curve
   // for electrons
   // exclusion of the crossing points
   //

   if(!fkPIDResponse) return 0;

   AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;

   // Make clone track in order to be able to apply correction
   const AliVTrack *rectrack;
   AliESDtrack esdtrack;
   AliAODTrack aodtrack;
   Double_t correctedTPCnSigma=0.;
   Bool_t TPCnSigmaCorrected=kFALSE;
   if((fkEtaMeanCorrection&&fkEtaWidthCorrection)|| (fkPMeanCorrection&&fkPWidthCorrection) ||
      (fkCentralityMeanCorrection&&fkCentralityWidthCorrection)){
      TPCnSigmaCorrected=kTRUE;
      correctedTPCnSigma=GetCorrectedTPCnSigma(track->GetRecTrack()->Eta(), track->GetMultiplicity(), fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron), track->GetRecTrack()->P());
   }
   // jpsi
   if((fkCentralityEtaCorrectionMeanJpsi)&&
      (fkCentralityEtaCorrectionWidthJpsi)){
      TPCnSigmaCorrected=kTRUE;
      correctedTPCnSigma=GetCorrectedTPCnSigmaJpsi(track->GetRecTrack()->Eta(), track->GetMultiplicity(), fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron));
   }
   if(fkEtaCorrection || fkCentralityCorrection){
      // Correction available
      // apply it on copy
      if(track->IsESDanalysis()){
         esdtrack.~AliESDtrack();
         new(&esdtrack) AliESDtrack(*(static_cast<const AliESDtrack *>(track->GetRecTrack())));
         if(track->IsPbPb() && HasCentralityCorrection())
            ApplyCentralityCorrection(&esdtrack, track->GetMultiplicity(), anatype);
         if(HasEtaCorrection())
            ApplyEtaCorrection(&esdtrack, anatype);
         rectrack = &esdtrack;
      } else {
         aodtrack.~AliAODTrack();
         new(&aodtrack) AliAODTrack(*(static_cast<const AliAODTrack *>(track->GetRecTrack())));
         if(track->IsPbPb() && HasCentralityCorrection())
            ApplyCentralityCorrection(&aodtrack, track->GetMultiplicity(), anatype);
         if(HasEtaCorrection())
            ApplyEtaCorrection(&aodtrack, anatype);
         rectrack = &aodtrack;
      }
   } else {
      // Correction not available - no need to copy
      rectrack = track->GetRecTrack();
   }
   AliHFEpidObject tpctrack(*track);
   tpctrack.SetRecTrack(rectrack);
   if(TPCnSigmaCorrected)tpctrack.SetCorrectedTPCnSigma(correctedTPCnSigma);

   // QA before selection (after correction)
   if(pidqa) pidqa->ProcessTrack(&tpctrack, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kBeforePID);
   AliDebug(1, "Doing TPC PID based on n-Sigma cut approach");

   // make copy of the track in order to allow for applying the correction
   Float_t nsigma=correctedTPCnSigma;
   if(!TPCnSigmaCorrected)
      nsigma = fUsedEdx ? rectrack->GetTPCsignal() : fkPIDResponse->NumberOfSigmasTPC(rectrack, AliPID::kElectron);
   AliDebug(1, Form("TPC NSigma: %f", nsigma));
   // exclude crossing points:
   // Determine the bethe values for each particle species
   Bool_t isLineCrossing = kFALSE;
   for(Int_t ispecies = 0; ispecies < AliPID::kSPECIES; ispecies++){
      if(ispecies == AliPID::kElectron) continue;
      if(!(fLineCrossingsEnabled & 1 << ispecies)) continue;
      if(TMath::Abs(fkPIDResponse->NumberOfSigmasTPC(rectrack, (AliPID::EParticleType)ispecies)) < fLineCrossingSigma[ispecies] && TMath::Abs(nsigma) < fNsigmaTPC){
         // Point in a line crossing region, no PID possible, but !PID still possible ;-)
         isLineCrossing = kTRUE;
         break;
      }
   }
   if(isLineCrossing) return 0;

   // Check particle rejection
   if(HasParticleRejection()){
      Int_t reject = Reject(rectrack, anatype);
      if(reject != 0) return reject;
   }

   // Check if we have an asymmetric sigma model set
   Int_t pdg = 0;
   if(fHasCutModel){
      pdg = CutSigmaModel(&tpctrack) ? 11 : 0;
   } else {
      // Perform Asymmetric n-sigma cut if required, else perform symmetric TPC sigma cut
      Float_t p = rectrack->P();
      if(HasAsymmetricSigmaCut()) {

         //printf("p %f, fPAsigCut[0] %f, fPAsigCut[1] %f\n",p,fPAsigCut[0],fPAsigCut[1]);
         if(p >= fPAsigCut[0] && p <= fPAsigCut[1]) {
            if(nsigma >= fNAsigmaTPC[0] && nsigma <= fNAsigmaTPC[1]) pdg = 11;
            else pdg = 0;
         }
         else pdg = 0;

      }
      else {
         if(TMath::Abs(nsigma) < fNsigmaTPC ) pdg = 11;
      }
   }
   if(pidqa && pdg != 0) pidqa->ProcessTrack(&tpctrack, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kAfterPID);
   return pdg;

}

//___________________________________________________________________
Bool_t AliHFEpidTPC::CutSigmaModel(const AliHFEpidObject * const track) const {
   //
   // N SigmaCut using parametrization of the cuts
   //
   Bool_t isSelected = kTRUE;
   AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
   Float_t nsigma = fUsedEdx ? track->GetRecTrack()->GetTPCsignal() : fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron);
   Double_t p = GetP(track->GetRecTrack(), anatype);
   Int_t centrality = track->IsPbPb() ? track->GetCentrality() + 1 : 0;
   AliDebug(2, Form("Centrality: %d\n", centrality));
   if(centrality > 11) return kFALSE;
   const TF1 *cutfunction;
   if((cutfunction = fkUpperSigmaCut[centrality]) && nsigma > cutfunction->Eval(p)) isSelected = kFALSE;
   if((cutfunction = fkLowerSigmaCut[centrality]) && nsigma < cutfunction->Eval(p)) isSelected = kFALSE;
   return isSelected;
}

//___________________________________________________________________
Int_t AliHFEpidTPC::Reject(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anaType) const{
   //
   // reject particles based on asymmetric sigma cut
   //
   Int_t pdc[AliPID::kSPECIES] = {11,13,211,321,2212};
   Double_t p = GetP(track, anaType);
   for(Int_t ispec = 0; ispec < AliPID::kSPECIES; ispec++){
      if(!TESTBIT(fRejectionEnabled, ispec)) continue;
      // Particle rejection enabled
      if(p < fRejection[4*ispec] || p > fRejection[4*ispec+2]) continue;
      Double_t sigma = fkPIDResponse->NumberOfSigmasTPC(track, static_cast<AliPID::EParticleType>(ispec));
      if(sigma >= fRejection[4*ispec+1] && sigma <= fRejection[4*ispec+3]) return pdc[ispec] * track->Charge();
   }
   return 0;
}

//___________________________________________________________________
void AliHFEpidTPC::ApplyEtaCorrection(AliVTrack *track, AliHFEpidObject::AnalysisType_t anatype) const{
   //
   // Apply correction for the eta dependence
   // N.B. This correction has to be applied on a copy track
   //
   AliDebug(1, Form("Applying correction function %s\n", fkEtaCorrection->GetName()));
   Double_t original = track->GetTPCsignal(),
   eta = track->Eta();
   if(anatype == AliHFEpidObject::kESDanalysis){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      if(fkEtaCorrection->Eval(eta)>0.0) esdtrack->SetTPCsignal(original/fkEtaCorrection->Eval(eta), esdtrack->GetTPCsignalSigma(), esdtrack->GetTPCsignalN());
   } else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      AliAODPid *pid = aodtrack->GetDetPid();
      if(pid && fkEtaCorrection->Eval(eta)>0.0) pid->SetTPCsignal(original/fkEtaCorrection->Eval(eta));
   }
}

//___________________________________________________________________
void AliHFEpidTPC::ApplyCentralityCorrection(AliVTrack *track, Double_t centralityEstimator, AliHFEpidObject::AnalysisType_t anatype) const{
   //
   // Apply correction for the eta dependence
   // N.B. This correction has to be applied on a copy track
   //
   AliDebug(1, Form("Applying correction function %s\n", fkCentralityCorrection->GetName()));
   Double_t original = track->GetTPCsignal();
   if(anatype == AliHFEpidObject::kESDanalysis){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      if(fkCentralityCorrection->Eval(centralityEstimator)>0.0) esdtrack->SetTPCsignal(original/fkCentralityCorrection->Eval(centralityEstimator), esdtrack->GetTPCsignalSigma(), esdtrack->GetTPCsignalN());
   } else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      AliAODPid *pid = aodtrack->GetDetPid();
      if(pid && fkCentralityCorrection->Eval(centralityEstimator)>0.0) pid->SetTPCsignal(original/fkCentralityCorrection->Eval(centralityEstimator));
   }
}

//___________________________________________________________________
Double_t AliHFEpidTPC::GetCorrectedTPCnSigma(Double_t eta, Double_t centralityEstimator, Double_t tpcNsigma, Double_t mom) const{
   //
   // Apply correction for the eta dependence
   // N.B. This correction has to be applied on a copy track
   //
   Double_t corrtpcNsigma = tpcNsigma;
   if(fkPMeanCorrection&&fkPWidthCorrection){
      if(TMath::Abs(fkPWidthCorrection->Eval(mom))>0.0000001) corrtpcNsigma=(corrtpcNsigma-fkPMeanCorrection->Eval(mom))/fkPWidthCorrection->Eval(mom);
   }
   if(fkEtaMeanCorrection&&fkEtaWidthCorrection){
      if(TMath::Abs(fkEtaWidthCorrection->Eval(eta))>0.0000001) corrtpcNsigma=(corrtpcNsigma-fkEtaMeanCorrection->Eval(eta))/fkEtaWidthCorrection->Eval(eta);
   }
   if(fkCentralityMeanCorrection&&fkCentralityWidthCorrection) {
      if(TMath::Abs(fkCentralityWidthCorrection->Eval(centralityEstimator))>0.0000001) corrtpcNsigma=(corrtpcNsigma-fkCentralityMeanCorrection->Eval(centralityEstimator))/fkCentralityWidthCorrection->Eval(centralityEstimator);
   }
   return corrtpcNsigma;
}

//___________________________________________________________________
Double_t AliHFEpidTPC::GetCorrectedTPCnSigmaJpsi(Double_t eta, Double_t centralityEstimator, Double_t tpcNsigma) const{
   //
   // Apply correction for the eta dependence
   // N.B. This correction has to be applied on a copy track
   //
   Double_t corrtpcNsigma = tpcNsigma;
   if(fkCentralityEtaCorrectionMeanJpsi&&fkCentralityEtaCorrectionWidthJpsi){
      const TAxis *caxis = fkCentralityEtaCorrectionMeanJpsi->GetXaxis();
      const TAxis *eaxis = fkCentralityEtaCorrectionMeanJpsi->GetYaxis();
      Int_t cbin = caxis->FindFixBin(centralityEstimator);
      Int_t ebin = eaxis->FindFixBin(eta);
      //Double_t cbinlowedge = caxis->GetBinLowEdge(cbin);
      //Double_t cbinupedge = caxis->GetBinUpEdge(cbin);
      //Double_t ebinlowedge = eaxis->GetBinLowEdge(ebin);
      //Double_t ebinupedge = eaxis->GetBinUpEdge(ebin);
      Double_t center = fkCentralityEtaCorrectionMeanJpsi->GetBinContent(cbin,ebin);
      Double_t width = fkCentralityEtaCorrectionWidthJpsi->GetBinContent(cbin,ebin);
      //printf("cbin %d, cbinlowe %f, cbinupe %f, centrality %f\n",cbin,cbinlowedge,cbinupedge,centralityEstimator);
      //printf("ebin %d, ebinlowe %f, ebinupe %f, eta %f\n",ebin,ebinlowedge,ebinupedge,eta);
      //printf("mean %f, width %f\n",center,width);
      if(TMath::Abs(width)>0.0000001) corrtpcNsigma=(corrtpcNsigma-center)/width;
   }
   return corrtpcNsigma;
}

//___________________________________________________________________
void AliHFEpidTPC::UseOROC(AliVTrack *track, AliHFEpidObject::AnalysisType_t anatype) const{
   //
   // Use TPC signal from the OROC
   // N.B. This correction has to be applied on a copy track
   //
   //Double_t original = track->GetTPCsignal();

   AliTPCdEdxInfo dEdxInfo;
   if( !track->GetTPCdEdxInfo(dEdxInfo) ) return;

   Double32_t  TPCsignalRegion[4]; // TPC dEdx signal in 4 different regions - 0 - IROC, 1- OROC medium, 2 - OROC long, 3- OROC all, (default truncation used)
   Char_t      TPCsignalNRegion[3]; // number of clusters above threshold used in the dEdx calculation
   Char_t      TPCsignalNRowRegion[3]; // number of crosed rows used in the dEdx calculation - signal below threshold included
   dEdxInfo.GetTPCSignalRegionInfo(TPCsignalRegion,TPCsignalNRegion,TPCsignalNRowRegion);

   if(anatype == AliHFEpidObject::kESDanalysis){
      AliESDtrack *esdtrack = static_cast<AliESDtrack *>(track);
      esdtrack->SetTPCsignal(TPCsignalRegion[3],esdtrack->GetTPCsignalSigma(),(TPCsignalNRegion[1]+TPCsignalNRegion[2])); // the two last are not ok
   } else {
      AliAODTrack *aodtrack = static_cast<AliAODTrack *>(track);
      AliAODPid *pid = aodtrack->GetDetPid();
      if(pid) pid->SetTPCsignal(TPCsignalRegion[3]);
      if(pid) pid->SetTPCsignalN((TPCsignalNRegion[1]+TPCsignalNRegion[2]));
   }

}
//___________________________________________________________________
void AliHFEpidTPC::AddTPCdEdxLineCrossing(Int_t species, Double_t sigma){
   //
   // Add exclusion point for the TPC PID where a dEdx line crosses the electron line
   // Stores line center and line sigma
   //
   if(species >= AliPID::kSPECIES){
      AliError("Species doesn't exist");
      return;
   }
   fLineCrossingsEnabled |= 1 << species;
   fLineCrossingSigma[species] = sigma;
}

//___________________________________________________________________
Double_t AliHFEpidTPC::GetP(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype) const {
   //
   // Get the momentum at the inner wall of the TPC
   //
   Double_t p = -1;
   if(anatype == AliHFEpidObject::kESDanalysis){
      // ESD analysis: Use Inner Params for the momentum estimate
      const AliESDtrack *esdtrack = dynamic_cast<const AliESDtrack *>(track);
      if(esdtrack) p = esdtrack->GetInnerParam() ? esdtrack->GetInnerParam()->GetP() : esdtrack->P();
   }

   if(anatype == AliHFEpidObject::kAODanalysis)
   {
      // AOD analysis: Use TPC momentum stored in the AliAODpid object
      const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
      if(aodtrack) p = aodtrack->GetDetPid() ? aodtrack->GetDetPid()->GetTPCmomentum() : aodtrack->P();
   }
   return p;
}
