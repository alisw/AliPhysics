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
  , fHasCutModel(kFALSE)
  , fNsigmaTPC(3)
  , fRejectionEnabled(0)
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
  , fHasCutModel(kFALSE)
  , fNsigmaTPC(3)
  , fRejectionEnabled(0)
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
  , fHasCutModel(ref.fHasCutModel)
  , fNsigmaTPC(2)
  , fRejectionEnabled(0)
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

  target.fLineCrossingsEnabled = fLineCrossingsEnabled;
  target.fHasCutModel = fHasCutModel;
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
  
  // QA before selection
  if(pidqa) pidqa->ProcessTrack(track, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kBeforePID);
  AliDebug(1, "Doing TPC PID based on n-Sigma cut approach");
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t nsigma = fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron);
  AliDebug(1, Form("TPC NSigma: %f", nsigma));
  // exclude crossing points:
  // Determine the bethe values for each particle species
  Bool_t isLineCrossing = kFALSE;
  for(Int_t ispecies = 0; ispecies < AliPID::kSPECIES; ispecies++){
    if(ispecies == AliPID::kElectron) continue;
    if(!(fLineCrossingsEnabled & 1 << ispecies)) continue;
    if(TMath::Abs(fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), (AliPID::EParticleType)ispecies)) < fLineCrossingSigma[ispecies] && TMath::Abs(nsigma) < fNsigmaTPC){
      // Point in a line crossing region, no PID possible, but !PID still possible ;-)
      isLineCrossing = kTRUE;      
      break;
    }
  }
  if(isLineCrossing) return 0;

  // Check particle rejection
  if(HasParticleRejection()){
    Int_t reject = Reject(track->GetRecTrack(), anatype);
    if(reject != 0) return reject;
  }

  // Check if we have an asymmetric sigma model set
  Int_t pdg = 0;
  if(fHasCutModel){
    pdg = CutSigmaModel(track) ? 11 : 0;
  } else { 
    // Perform Asymmetric n-sigma cut if required, else perform symmetric TPC sigma cut
    Float_t p = 0.;
    if(HasAsymmetricSigmaCut() && (p = track->GetRecTrack()->P()) >= fPAsigCut[0] && p <= fPAsigCut[1]){ 
      if(nsigma >= fNAsigmaTPC[0] && nsigma <= fNAsigmaTPC[1]) pdg = 11; 
    } else {
      if(TMath::Abs(nsigma) < fNsigmaTPC ) pdg = 11;
    }
  }
  if(pidqa && pdg != 0) pidqa->ProcessTrack(track, AliHFEpid::kTPCpid, AliHFEdetPIDqa::kAfterPID);
  return pdg;

}

//___________________________________________________________________
Bool_t AliHFEpidTPC::CutSigmaModel(const AliHFEpidObject * const track) const {
  //
  // N SigmaCut using parametrization of the cuts
  //
  Bool_t isSelected = kTRUE;
  AliHFEpidObject::AnalysisType_t anatype = track->IsESDanalysis() ? AliHFEpidObject::kESDanalysis : AliHFEpidObject::kAODanalysis;
  Float_t nsigma = fkPIDResponse->NumberOfSigmasTPC(track->GetRecTrack(), AliPID::kElectron);
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
  } else { 
    // AOD analysis: Use TPC momentum stored in the AliAODpid object
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(track);
    if(aodtrack) p = aodtrack->GetDetPid() ? aodtrack->GetDetPid()->GetTPCmomentum() : aodtrack->P();
  }
  return p;
}
