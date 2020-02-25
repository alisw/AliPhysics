#include "AliAnalysisNanoAODCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliNanoAODHeader.h"
#include "AliNanoAODTrack.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODConversionPhoton.h"
#include "AliAODv0.h"
#include "AliEventCuts.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisNanoAODTrackCuts)
ClassImp(AliAnalysisNanoAODV0Cuts)
ClassImp(AliAnalysisNanoAODCascadeCuts)
ClassImp(AliAnalysisNanoAODPhotonCuts)
ClassImp(AliAnalysisNanoAODEventCuts)
ClassImp(AliAnalysisNanoAODMCParticleCuts)
ClassImp(AliNanoAODSimpleSetter)


AliAnalysisNanoAODTrackCuts::AliAnalysisNanoAODTrackCuts():
AliAnalysisCuts(), fBitMask(1), fMinPt(0), fMaxEta(10)
{
  // default ctor 
}

Bool_t AliAnalysisNanoAODTrackCuts::IsSelected(TObject* obj)
{
  // Returns true if the track is good!
  AliAODTrack* track = dynamic_cast<AliAODTrack*>(obj);
  
  
  if (!track->TestFilterBit(fBitMask))    return kFALSE;
  if (track->Pt() < fMinPt)               return kFALSE;
  if (TMath::Abs(track->Eta()) > fMaxEta) return kFALSE; 
  
  return kTRUE;  
}

AliAnalysisNanoAODV0Cuts::AliAnalysisNanoAODV0Cuts()
    : AliAnalysisCuts(),
      fSelectOnFly(false),
      fOnFlyStatus(false),
      fv0pTMin(-1),
      fv0EtaMax(-1),
      fTransverseRadiusMin(-1),
      fTransverseRadiusMax(-1),
      fCPAMin(-1),
      fDCADaugv0VtxMax(-1),
      fDCADaugPrimVtxMin(-1),
      fDaughEtaMax(-1),
      fDaugMinClsTPC(-1),
      fLambdaDaugnSigTPCMax(-1),
      fCheckDaughterPileup(kFALSE),
      fCheckDaughterTPCRefit(kFALSE)
{
}

Bool_t AliAnalysisNanoAODV0Cuts::IsSelected(TObject* obj) 
{
  // V0 selection
  
  AliAODv0* v0 = dynamic_cast<AliAODv0*>(obj);
  if (!v0)
    AliFatal("Did not pass a V0");

  if (v0->GetNProngs() != 2)
    return false;

  if (v0->GetNDaughters() != 2)
    return false;

  if (fSelectOnFly && v0->GetOnFlyStatus() == fOnFlyStatus)
    return false;

  if (fv0pTMin > 0 && v0->Pt() < fv0pTMin)
    return false;

  if (fv0EtaMax > 0 && TMath::Abs(v0->Eta()) > fv0EtaMax)
    return false;
  
  if (fTransverseRadiusMin > 0 || fCPAMin > 0) {
    const AliAODEvent* evt = static_cast<const AliAODEvent*>(static_cast<AliAODTrack *>(v0->GetDaughter(0))->GetEvent());
    if (!evt)
      AliFatal("No event but cut on transverse radius or/and pointing angle requested");

    Float_t xvP = evt->GetPrimaryVertex()->GetX();
    Float_t yvP = evt->GetPrimaryVertex()->GetY();
    Float_t zvP = evt->GetPrimaryVertex()->GetZ();
    Double_t vecTarget[3] = { xvP, yvP, zvP };
    Float_t TransverseRadius = v0->DecayLengthXY(vecTarget);
    if ((fTransverseRadiusMin > 0 || fTransverseRadiusMax > 0) && (TransverseRadius > fTransverseRadiusMax || TransverseRadius < fTransverseRadiusMin))
      return false;
    if (fCPAMin > 0 && v0->CosPointingAngle(vecTarget) < fCPAMin)
      return false;
  }
  
  if (fDCADaugv0VtxMax > 0 && v0->DcaV0Daughters() > fDCADaugv0VtxMax)
    return false;
  if (fDCADaugPrimVtxMin > 0 && (v0->DcaPosToPrimVertex() < fDCADaugPrimVtxMin || v0->DcaNegToPrimVertex() < fDCADaugPrimVtxMin))
    return false;

  //TODO: For the preselection it does not matter if the assignment positive and
  //negative is correct, during creation of the Nano AOD this should be checked
  AliAODTrack *pTrack = static_cast<AliAODTrack *>(v0->GetDaughter(0));
  AliAODTrack *nTrack = static_cast<AliAODTrack *>(v0->GetDaughter(1));
  if (fDaughEtaMax > 0 && (TMath::Abs(pTrack->Eta()) > fDaughEtaMax || TMath::Abs(nTrack->Eta()) > fDaughEtaMax))
    return false;

  if (fDaugMinClsTPC > 0 && (pTrack->GetTPCNcls() < fDaugMinClsTPC || nTrack->GetTPCNcls() < fDaugMinClsTPC))
    return false;
  if (fCheckDaughterTPCRefit && !pTrack->IsOn(AliAODTrack::kTPCrefit)){
    return false;
  }
  if (fCheckDaughterTPCRefit && !nTrack->IsOn(AliAODTrack::kTPCrefit)){
    return false;
  }
  if (fLambdaDaugnSigTPCMax > 0) {
    static AliPIDResponse* pidResponse = 0;
    if (!pidResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
      pidResponse = inputHandler->GetPIDResponse();
      if (!pidResponse)
        AliFatal("No PID response but PID selection requested");
    }
    AliPIDResponse::EDetPidStatus statusPosTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, pTrack);
    AliPIDResponse::EDetPidStatus statusNegTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, nTrack);
    if (!(statusPosTPC == AliPIDResponse::kDetPidOk && statusNegTPC == AliPIDResponse::kDetPidOk))
      return false;

    Float_t nSigPosPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kPion);
    Float_t nSigPosProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kProton);
    Float_t nSigNegPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kPion);
    Float_t nSigNegProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kProton);

    // if the daughter tracks are not a proton or a pion within loose cuts, the candidate can be rejected
    if (!((nSigPosPion < fLambdaDaugnSigTPCMax && nSigNegProton < fLambdaDaugnSigTPCMax) ||
          (nSigNegPion < fLambdaDaugnSigTPCMax && nSigPosProton < fLambdaDaugnSigTPCMax)))
      return false;
  }
  
  if (fCheckDaughterPileup) {
    if (!((pTrack->GetTOFBunchCrossing() == 0) || (nTrack->GetTOFBunchCrossing() == 0))) {
      Bool_t PileUpPass = false;
      for (Int_t iLay = 0; iLay < 6; ++iLay) {
        //do not use SDD for Pile Up Rejection since detector is too slow
        if (iLay == 2 || iLay == 3)
          continue;

        if (pTrack->HasPointOnITSLayer(iLay) || nTrack->HasPointOnITSLayer(iLay)) {
          PileUpPass = true;
          break;
        }
      }
      if (!PileUpPass)
        return false;
    }
  }
  
  return kTRUE;
}

AliAnalysisNanoAODV0ParametricCuts::AliAnalysisNanoAODV0ParametricCuts()
: AliAnalysisCuts(),
fSelectOnFly(false),
fOnFlyStatus(false),
fv0pTMin(-1),
fTransverseRadiusMin(-1),
fTransverseRadiusMax(-1),
fCPAMin(-1),
fMaxDCAV0Dau(-1),
fDCADauToPV(-1),
fDaughEtaMax(-1),
fDaugMinClsTPC(-1),
fLambdaDaugnSigTPCMax(-1),
fCheckDaughterPileup(kFALSE),
fUseParametricCosPA(kFALSE),
fParCosPA{0.}
{
}

void AliAnalysisNanoAODV0ParametricCuts::SetupDefaultPbPb2015Cuts()
{
    fSelectOnFly = kFALSE;
    fOnFlyStatus = kFALSE;
    fv0pTMin = -1;
    fTransverseRadiusMin = 4.5;
    fTransverseRadiusMax = 200;
    fCPAMin = 0.95;
    fMaxDCAV0Dau = 1.0;
    fDCADauToPV = 0.05;
    fDaughEtaMax = 0.8;
    fDaugMinClsTPC = -1 ;
    fLambdaDaugnSigTPCMax = -1;
    fCheckDaughterPileup = kFALSE ;
    fUseParametricCosPA = kTRUE;
    fParCosPA[0] = 0.20428;
    fParCosPA[1] = -0.73728;
    fParCosPA[2] = 0.09887;
    fParCosPA[3] = -0.02822;
    fParCosPA[4] = -0.05302;
}

Bool_t AliAnalysisNanoAODV0ParametricCuts::IsSelected(TObject* obj)
{
    // V0 selection
    
    AliAODv0* v0 = dynamic_cast<AliAODv0*>(obj);
    if (!v0)
        AliFatal("Did not pass a V0");
    
    if (v0->GetNProngs() != 2)
        return false;
    
    if (v0->GetNDaughters() != 2)
        return false;
    
    if (fSelectOnFly && v0->GetOnFlyStatus() == fOnFlyStatus)
        return false;
    
    if (fv0pTMin > 0 && v0->Pt() < fv0pTMin)
        return false;
    
    if (fTransverseRadiusMin > 0 || fCPAMin > 0) {
        const AliAODEvent* evt = static_cast<const AliAODEvent*>(static_cast<AliAODTrack *>(v0->GetDaughter(0))->GetEvent());
        if (!evt)
            AliFatal("No event but cut on transverse radius or/and pointing angle requested");
        
        Float_t xvP = evt->GetPrimaryVertex()->GetX();
        Float_t yvP = evt->GetPrimaryVertex()->GetY();
        Float_t zvP = evt->GetPrimaryVertex()->GetZ();
        Double_t vecTarget[3] = { xvP, yvP, zvP };
        Float_t TransverseRadius = v0->DecayLengthXY(vecTarget);
        if (TransverseRadius > fTransverseRadiusMax || TransverseRadius < fTransverseRadiusMin)
            return false;
        
        Float_t lCosPACut = fCPAMin;
        Float_t lVarV0CosPA = TMath::Cos(
                                         fParCosPA[0]*TMath::Exp(fParCosPA[1]*v0->Pt()) +
                                         fParCosPA[2]*TMath::Exp(fParCosPA[3]*v0->Pt()) +
                                         fParCosPA[4]);
        if( fUseParametricCosPA && lCosPACut<lVarV0CosPA )
            lCosPACut=lVarV0CosPA;
        
        if (fCPAMin > 0 && v0->CosPointingAngle(vecTarget) < lCosPACut)
            return false;
    }
    
    if (fMaxDCAV0Dau > 0 && v0->DcaV0Daughters() > fMaxDCAV0Dau)
        return false;
    if (fDCADauToPV > 0 && (v0->DcaPosToPrimVertex() < fDCADauToPV || v0->DcaNegToPrimVertex() < fDCADauToPV))
        return false;
    
    //TODO: For the preselection it does not matter if the assignment positive and
    //negative is correct, during creation of the Nano AOD this should be checked
    AliAODTrack *pTrack = static_cast<AliAODTrack *>(v0->GetDaughter(0));
    AliAODTrack *nTrack = static_cast<AliAODTrack *>(v0->GetDaughter(1));
    if (fDaughEtaMax > 0 && (TMath::Abs(pTrack->Eta()) > fDaughEtaMax || TMath::Abs(nTrack->Eta()) > fDaughEtaMax))
        return false;
    
    if (fDaugMinClsTPC > 0 && (pTrack->GetTPCNcls() < fDaugMinClsTPC || nTrack->GetTPCNcls() < fDaugMinClsTPC))
        return false;
    
    if (fLambdaDaugnSigTPCMax > 0) {
        static AliPIDResponse* pidResponse = 0;
        if (!pidResponse) {
            AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
            AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
            pidResponse = inputHandler->GetPIDResponse();
            if (!pidResponse)
                AliFatal("No PID response but PID selection requested");
        }
        AliPIDResponse::EDetPidStatus statusPosTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, pTrack);
        AliPIDResponse::EDetPidStatus statusNegTPC = pidResponse->CheckPIDStatus(AliPIDResponse::kTPC, nTrack);
        if (!(statusPosTPC == AliPIDResponse::kDetPidOk && statusNegTPC == AliPIDResponse::kDetPidOk))
            return false;
        
        Float_t nSigPosPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kPion);
        Float_t nSigPosProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, pTrack, AliPID::kProton);
        Float_t nSigNegPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kPion);
        Float_t nSigNegProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC, nTrack, AliPID::kProton);
        
        // if the daughter tracks are not a proton or a pion within loose cuts, the candidate can be rejected
        if (!((nSigPosPion < fLambdaDaugnSigTPCMax && nSigNegProton < fLambdaDaugnSigTPCMax) ||
              (nSigNegPion < fLambdaDaugnSigTPCMax && nSigPosProton < fLambdaDaugnSigTPCMax)))
            return false;
    }
    
    if (fCheckDaughterPileup) {
        if (!((pTrack->GetTOFBunchCrossing() == 0) || (nTrack->GetTOFBunchCrossing() == 0))) {
            Bool_t PileUpPass = false;
            for (Int_t iLay = 0; iLay < 6; ++iLay) {
                //do not use SDD for Pile Up Rejection since detector is too slow
                if (iLay == 2 || iLay == 3)
                    continue;
                
                if (pTrack->HasPointOnITSLayer(iLay) || nTrack->HasPointOnITSLayer(iLay)) {
                    PileUpPass = true;
                    break;
                }
            }
            if (!PileUpPass)
                return false;
        }
    }
    
    return kTRUE;
}

AliAnalysisNanoAODCascadeCuts::AliAnalysisNanoAODCascadeCuts()
    : AliAnalysisCuts(),
      fCascpTMin(-99.),
      fDCADaugPrimVtxMin(-99.),
      fCPACascMin(-99.),
      fTransverseRadiusCasc(-99.),
      fCPAv0Min(-99.),
      fTransverseRadiusv0(-99.),
      fDCAv0PrimVtxMin(-99.),
      fDaughEtaMax(-99.),
      fCascDaugnSigTPCMax(-99.),
      fCheckDaughterPileup(false),
      fCheckDaughterTPCRefit(false) {
}

Bool_t AliAnalysisNanoAODCascadeCuts::IsSelected(TObject* obj) {
  AliAODcascade* cascade = dynamic_cast<AliAODcascade*>(obj);
  if (!cascade)
    AliFatal("Did not pass a Cascade");
  TVector3 pCasc(cascade->MomXiX(), cascade->MomXiY(), cascade->MomXiZ());
  if (fCascpTMin > 0. && pCasc.Pt() < fCascpTMin) {
    return false;
  }
  if (fDCADaugPrimVtxMin > 0.
      && cascade->DcaBachToPrimVertex() < fDCADaugPrimVtxMin) {
    return false;
  }

  if (fDCADaugPrimVtxMin > 0.
      && cascade->DcaNegToPrimVertex() < fDCADaugPrimVtxMin) {
    return false;
  }

  if (fDCADaugPrimVtxMin > 0.
      && cascade->DcaPosToPrimVertex() < fDCADaugPrimVtxMin) {
    return false;
  }
  if (fCPACascMin > 0. || fCPAv0Min > 0.) {
    const AliAODEvent* evt =
        static_cast<const AliAODEvent*>(static_cast<AliAODTrack *>(cascade
            ->GetDaughter(0))->GetEvent());
    if (!evt)
      AliFatal("No event but cut pointing angles requested");
    Float_t xvP = evt->GetPrimaryVertex()->GetX();
    Float_t yvP = evt->GetPrimaryVertex()->GetY();
    Float_t zvP = evt->GetPrimaryVertex()->GetZ();
    // Double_t vecTarget[3] = { xvP, yvP, zvP };
    if (fCPACascMin > 0.
        && cascade->CosPointingAngleXi(xvP, yvP, zvP) < fCPACascMin) {
      return false;
    }
    if (fCPAv0Min > 0.
        && cascade->CosPointingAngle(evt->GetPrimaryVertex()) < fCPAv0Min) {
      return false;
    }
  }
  if (fTransverseRadiusCasc > 0.
      && cascade->DecayVertexXiX() > fTransverseRadiusCasc) {
    return false;
  }
  if (fTransverseRadiusCasc > 0.
      && cascade->DecayVertexXiY() > fTransverseRadiusCasc) {
    return false;
  }
  if (fTransverseRadiusv0 > 0.
      && cascade->DecayVertexV0X() > fTransverseRadiusv0) {
    return false;
  }
  if (fTransverseRadiusv0 > 0.
      && cascade->DecayVertexV0Y() > fTransverseRadiusv0) {
    return false;
  }
  if (fDCAv0PrimVtxMin > 0. && cascade->DcaV0ToPrimVertex() < fDCAv0PrimVtxMin) {
    return false;
  }
  AliAODTrack *pTrack = static_cast<AliAODTrack*>(cascade->GetDaughter(0));
  AliAODTrack *nTrack = static_cast<AliAODTrack*>(cascade->GetDaughter(1));
  AliAODTrack *bachTrack = static_cast<AliAODTrack*>(cascade->GetDecayVertexXi()
      ->GetDaughter(0));
  if (fDaughEtaMax > 0.) {
    if (TMath::Abs(pTrack->Eta()) > fDaughEtaMax) {
      return false;
    }
    if (TMath::Abs(nTrack->Eta()) > fDaughEtaMax) {
      return false;
    }
    if (TMath::Abs(bachTrack->Eta()) > fDaughEtaMax) {
      return false;
    }
  }
  if (fCheckDaughterTPCRefit && !pTrack->IsOn(AliAODTrack::kTPCrefit)){
    return false;
  }
  if (fCheckDaughterTPCRefit && !nTrack->IsOn(AliAODTrack::kTPCrefit)){
    return false;
  }
  if (fCheckDaughterTPCRefit && !bachTrack->IsOn(AliAODTrack::kTPCrefit)){
    return false;
  }
  if (fCascDaugnSigTPCMax > 0.) {
    static AliPIDResponse* pidResponse = 0;
    if (!pidResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
          ->GetInputEventHandler());
      pidResponse = inputHandler->GetPIDResponse();
      if (!pidResponse)
        AliFatal("No PID response but PID selection requested");
    }
    AliPIDResponse::EDetPidStatus statusPosTPC = pidResponse->CheckPIDStatus(
        AliPIDResponse::kTPC, pTrack);
    AliPIDResponse::EDetPidStatus statusNegTPC = pidResponse->CheckPIDStatus(
        AliPIDResponse::kTPC, nTrack);
    AliPIDResponse::EDetPidStatus statusBachTPC = pidResponse->CheckPIDStatus(
        AliPIDResponse::kTPC, bachTrack);
    if (!(statusPosTPC == AliPIDResponse::kDetPidOk
        && statusNegTPC == AliPIDResponse::kDetPidOk
        && statusBachTPC == AliPIDResponse::kDetPidOk))
      return false;

    Float_t nSigPosPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                      pTrack, AliPID::kPion);
    Float_t nSigPosProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                        pTrack,
                                                        AliPID::kProton);
    Float_t nSigNegPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                      nTrack, AliPID::kPion);
    Float_t nSigNegProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                        nTrack,
                                                        AliPID::kProton);
    Float_t nSigBachPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                       bachTrack,
                                                        AliPID::kPion);

    Float_t nSigBachKaon = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                       bachTrack,
                                                        AliPID::kKaon);

    // if the Bachelor is not a pion or a kaon, the candidate can go
    if (!((nSigBachPion < fCascDaugnSigTPCMax)||
        (nSigBachKaon < fCascDaugnSigTPCMax))) {
      return false;
    }
    // if the daughter tracks are not a proton or a pion within loose cuts, the candidate can be rejected
    if (!((nSigPosPion < fCascDaugnSigTPCMax
        && nSigNegProton < fCascDaugnSigTPCMax)
        || (nSigNegPion < fCascDaugnSigTPCMax
            && nSigPosProton < fCascDaugnSigTPCMax)))
      return false;
  }

  if (fCheckDaughterPileup) {
    if (!((pTrack->GetTOFBunchCrossing() == 0)
        || (nTrack->GetTOFBunchCrossing() == 0)
        || (bachTrack->GetTOFBunchCrossing() == 0))) {
      Bool_t PileUpPass = false;
      for (Int_t iLay = 0; iLay < 6; ++iLay) {
        //do not use SDD for Pile Up Rejection since detector is too slow
        if (iLay == 2 || iLay == 3)
          continue;

        if (pTrack->HasPointOnITSLayer(iLay) || nTrack->HasPointOnITSLayer(iLay)
            || bachTrack->HasPointOnITSLayer(iLay)) {
          PileUpPass = true;
          break;
        }
      }
      if (!PileUpPass)
        return false;
    }
  }
  return true;
}

AliAnalysisNanoAODCascadeParametricCuts::AliAnalysisNanoAODCascadeParametricCuts()
: AliAnalysisCuts(),
fCascpTMin(-1),
fDCAV0DauToPV(-1),
fDCABachToPV(-1),
fCPACascMin(-1),
fCPAv0Min(-1),
fTransverseRadiusCasc(-1),
fTransverseRadiusv0(0.),
fDCAv0PrimVtxMin(0.),
fDaughEtaMax(0.),
fCascDaugnSigTPCMax(0.),
fCheckDaughterPileup(false),
fUseParametricV0CosPA(kFALSE),
fParV0CosPA{0.},
fUseParametricCascCosPA(kFALSE),
fParCascCosPA{0.}
{
}

void AliAnalysisNanoAODCascadeParametricCuts::SetupDefaultPbPb2015Cuts()
{
    fCascpTMin = 0.6;
    fDCAV0DauToPV = 0.1;
    fDCABachToPV = 0.05;
    fCPACascMin = 0.95;
    fCPAv0Min = 0.95;
    fTransverseRadiusCasc = 0.6;
    fTransverseRadiusv0 = 2;
    fDCAv0PrimVtxMin = 0.05;
    fDaughEtaMax = 0.8;
    fCascDaugnSigTPCMax = 5.;
    fCheckDaughterPileup = kFALSE;
    fUseParametricV0CosPA = kTRUE;
    fParV0CosPA[0] = TMath::Exp(  -1.77429);
    fParV0CosPA[1] =  -0.692453;
    fParV0CosPA[2] = TMath::Exp( -2.01938);
    fParV0CosPA[3] = -0.201574;
    fParV0CosPA[4] = 0.0776465;
    fUseParametricCascCosPA = kTRUE;
    fParCascCosPA[0] = TMath::Exp(  12.8077);
    fParCascCosPA[1] =  -21.2944;
    fParCascCosPA[2] = TMath::Exp( -1.53357);
    fParCascCosPA[3] = -0.920017;
    fParCascCosPA[4] = 0.0262315;
}

Bool_t AliAnalysisNanoAODCascadeParametricCuts::IsSelected(TObject* obj) {
    AliAODcascade* cascade = dynamic_cast<AliAODcascade*>(obj);
    if (!cascade)
        AliFatal("Did not pass a Cascade");
    TVector3 pCasc(cascade->MomXiX(), cascade->MomXiY(), cascade->MomXiZ());
    if (fCascpTMin > 0. && pCasc.Pt() < fCascpTMin) {
        return false;
    }
    if (fDCABachToPV > 0.
        && cascade->DcaBachToPrimVertex() < fDCABachToPV) {
        return false;
    }
    
    if (fDCAV0DauToPV > 0.
        && cascade->DcaNegToPrimVertex() < fDCAV0DauToPV) {
        return false;
    }
    
    if (fDCAV0DauToPV > 0.
        && cascade->DcaPosToPrimVertex() < fDCAV0DauToPV) {
        return false;
    }
    if (fCPACascMin > 0. || fCPAv0Min > 0.) {
        const AliAODEvent* evt =
        static_cast<const AliAODEvent*>(static_cast<AliAODTrack *>(cascade
                                                                   ->GetDaughter(0))->GetEvent());
        if (!evt)
            AliFatal("No event but cut pointing angles requested");
        Float_t xvP = evt->GetPrimaryVertex()->GetX();
        Float_t yvP = evt->GetPrimaryVertex()->GetY();
        Float_t zvP = evt->GetPrimaryVertex()->GetZ();
        // Double_t vecTarget[3] = { xvP, yvP, zvP };
        //======== casc cosPA =========================
        Float_t lCascCosPACut = fCPACascMin;
        Float_t lVarCascCosPA = TMath::Cos(
                                         fParCascCosPA[0]*TMath::Exp(fParCascCosPA[1]*cascade->Pt()) +
                                         fParCascCosPA[2]*TMath::Exp(fParCascCosPA[3]*cascade->Pt()) +
                                         fParCascCosPA[4]);
        if( fUseParametricCascCosPA && lCascCosPACut<lVarCascCosPA )
            lCascCosPACut=lVarCascCosPA;
        if (lCascCosPACut > 0.
            && cascade->CosPointingAngleXi(xvP, yvP, zvP) < lCascCosPACut) {
            return false;
        }
        
        //======== V0 cosPA =========================
        Float_t lCosPACut = fCPAv0Min;
        Float_t lVarV0CosPA = TMath::Cos(
                                         fParV0CosPA[0]*TMath::Exp(fParV0CosPA[1]*cascade->Pt()) +
                                         fParV0CosPA[2]*TMath::Exp(fParV0CosPA[3]*cascade->Pt()) +
                                         fParV0CosPA[4]);
        if( fUseParametricV0CosPA && lCosPACut<lVarV0CosPA )
            lCosPACut=lVarV0CosPA;
        if (lCosPACut > 0.
            && cascade->CosPointingAngle(evt->GetPrimaryVertex()) < lCosPACut) {
            return false;
        }
    }
    if (fTransverseRadiusCasc > 0.
        && TMath::Sqrt(
                       TMath::Power(cascade->DecayVertexXiX(),2)+
                       TMath::Power(cascade->DecayVertexXiY(),2)
                       )< fTransverseRadiusCasc)
    {
        return false;
    }
    
    if (fTransverseRadiusv0 > 0.
        && TMath::Sqrt(
                       TMath::Power(cascade->DecayVertexV0X(),2)+
                       TMath::Power(cascade->DecayVertexV0Y(),2)
                       )< fTransverseRadiusv0)
    {
        return false;
    }
    
    if (fDCAv0PrimVtxMin > 0. && cascade->DcaV0ToPrimVertex() < fDCAv0PrimVtxMin) {
        return false;
    }
    AliAODTrack *pTrack = static_cast<AliAODTrack*>(cascade->GetDaughter(0));
    AliAODTrack *nTrack = static_cast<AliAODTrack*>(cascade->GetDaughter(1));
    AliAODTrack *bachTrack = static_cast<AliAODTrack*>(cascade->GetDecayVertexXi()
                                                       ->GetDaughter(0));
    if (fDaughEtaMax > 0.) {
        if (TMath::Abs(pTrack->Eta()) > fDaughEtaMax) {
            return false;
        }
        if (TMath::Abs(nTrack->Eta()) > fDaughEtaMax) {
            return false;
        }
        if (TMath::Abs(nTrack->Eta()) > fDaughEtaMax) {
            return false;
        }
    }
    
    if (fCascDaugnSigTPCMax > 0) {
        static AliPIDResponse* pidResponse = 0;
        if (!pidResponse) {
            AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
            AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man
                                                                          ->GetInputEventHandler());
            pidResponse = inputHandler->GetPIDResponse();
            if (!pidResponse)
                AliFatal("No PID response but PID selection requested");
        }
        AliPIDResponse::EDetPidStatus statusPosTPC = pidResponse->CheckPIDStatus(
                                                                                 AliPIDResponse::kTPC, pTrack);
        AliPIDResponse::EDetPidStatus statusNegTPC = pidResponse->CheckPIDStatus(
                                                                                 AliPIDResponse::kTPC, nTrack);
        AliPIDResponse::EDetPidStatus statusBachTPC = pidResponse->CheckPIDStatus(
                                                                                  AliPIDResponse::kTPC, bachTrack);
        if (!(statusPosTPC == AliPIDResponse::kDetPidOk
              && statusNegTPC == AliPIDResponse::kDetPidOk
              && statusBachTPC == AliPIDResponse::kDetPidOk))
            return false;
        
        Float_t nSigPosPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                          pTrack, AliPID::kPion);
        Float_t nSigPosProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                            pTrack,
                                                            AliPID::kProton);
        Float_t nSigNegPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                          nTrack, AliPID::kPion);
        Float_t nSigNegProton = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                            nTrack,
                                                            AliPID::kProton);
        Float_t nSigBachPion = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                           nTrack,
                                                           AliPID::kPion);
        Float_t nSigBachKaon = pidResponse->NumberOfSigmas(AliPIDResponse::kTPC,
                                                           nTrack,
                                                           AliPID::kKaon);
        
        // if the Bachelor is not a pion or kaon, the candidate can go
        if (! (nSigBachPion < fCascDaugnSigTPCMax) &&
            ! (nSigBachKaon < fCascDaugnSigTPCMax)) {
            return false;
        }
        // if the daughter tracks are not a proton or a pion within loose cuts, the candidate can be rejected
        if (!((nSigPosPion < fCascDaugnSigTPCMax
               && nSigNegProton < fCascDaugnSigTPCMax)
              || (nSigNegPion < fCascDaugnSigTPCMax
                  && nSigPosProton < fCascDaugnSigTPCMax)))
            return false;
    }
    
    if (fCheckDaughterPileup) {
        if (!((pTrack->GetTOFBunchCrossing() == 0)
              || (nTrack->GetTOFBunchCrossing() == 0)
              || (bachTrack->GetTOFBunchCrossing() == 0))) {
            Bool_t PileUpPass = false;
            for (Int_t iLay = 0; iLay < 6; ++iLay) {
                //do not use SDD for Pile Up Rejection since detector is too slow
                if (iLay == 2 || iLay == 3)
                    continue;
                
                if (pTrack->HasPointOnITSLayer(iLay) || nTrack->HasPointOnITSLayer(iLay)
                    || bachTrack->HasPointOnITSLayer(iLay)) {
                    PileUpPass = true;
                    break;
                }
            }
            if (!PileUpPass)
                return false;
        }
    }
    return true;
}

AliAnalysisNanoAODPhotonCuts::AliAnalysisNanoAODPhotonCuts()
    : AliAnalysisCuts(),
      fPhotonEtaMax(-1.f),
      fPhotonConvRadiusMin(-1.f),
      fPhotonConvRadiusMax(-1.f),
      fPhotonPsiPairMax(-1.f)
{ }

Bool_t AliAnalysisNanoAODPhotonCuts::IsSelected(TObject* obj) {
  // Photon selection

  auto photonCandidate = dynamic_cast<AliAODConversionPhoton *>(obj);
  if (!photonCandidate)
    AliFatal("Did not pass a Photon Candidate");

  if (fPhotonEtaMax > 0
      && TMath::Abs(photonCandidate->GetPhotonEta()) > fPhotonEtaMax)
    return false;

  if (fPhotonConvRadiusMin > 0) {
    const float transRadius = photonCandidate->GetConversionRadius();
    if (transRadius > fPhotonConvRadiusMax
        || transRadius < fPhotonConvRadiusMin)
      return false;
  }

  if (fPhotonPsiPairMax > 0
      && TMath::Abs(photonCandidate->GetPsiPair()) > fPhotonPsiPairMax) {
    return false;
  }

  return kTRUE;
}

AliAnalysisNanoAODEventCuts::AliAnalysisNanoAODEventCuts():
  AliAnalysisCuts(), 
  fTrackCut(0),
  fMinMultiplicity(-1),
  fMaxMultiplicity(-1),
  fEventCuts()
{
  // default ctor   
  
}

Bool_t AliAnalysisNanoAODEventCuts::IsSelected(TObject* obj)
{
  // Returns true if object accepted on the event level
  
  AliVEvent* evt = dynamic_cast<AliVEvent*>(obj);
  
  if (!fEventCuts.AcceptEvent(evt))
    return kFALSE;
  
  if (fTrackCut != 0)
  {
    Int_t trackCount = 0;
    for (Int_t j=0; j<evt->GetNumberOfTracks(); j++)
      if (fTrackCut->IsSelected(evt->GetTrack(j)))
        trackCount++;
      
    if (fMinMultiplicity > 0 && trackCount < fMinMultiplicity)
      return kFALSE;
    if (fMaxMultiplicity > 0 && trackCount > fMaxMultiplicity)
      return kFALSE;
  }
      
  return kTRUE;
}

AliAnalysisNanoAODMCParticleCuts::AliAnalysisNanoAODMCParticleCuts()
    : AliAnalysisCuts(),
      fDoSelectPrimaries(true),
      fDoSelectCharged(true),
      fMinPt(0.f),
      fMaxEta(999.f),
      fPDGToKeep(),
      fPDGV0(),
      fPDGV0Cascade(){
}

bool AliAnalysisNanoAODMCParticleCuts::IsSelected(TObject* obj) {
  // Returns true if the MC particle is to be kept!
  AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(obj);
  if(!part) {
    AliFatal("MC particle missing");
  }

  // if this is activated the specified particles are kept ignoring the cuts
  for (auto it : fPDGToKeep) {
    if (it == TMath::Abs(part->PdgCode())) {
      return true;
    }
  }
  if (part->Pt() < fMinPt) {
    return false;
  }
  if (TMath::Abs(part->Eta()) > fMaxEta) {
    return false;
  }
  if (fDoSelectPrimaries && !part->IsPhysicalPrimary()) {
    return false;
  }
  if (fDoSelectCharged && !part->Charge()) {
    return false;
  }
  return true;
}

void AliNanoAODSimpleSetter::Init(AliNanoAODHeader* head, TString varListHeader)
{
  // Initialize header
  
  TObjArray *vars = varListHeader.Tokenize(",");
  TIter it(vars);
  TObjString *token  = 0;
  Int_t index=0;
  Int_t indexInt=0;
  
  // HACK workaround for bug in initializer of AliNanoAODHeader in AliRoot v5-09-47b
  head->SetNumberOfESDTracksIndex(-1);
  
  std::map<TString,int> cstMap = head->GetMapCstVar();

  while ((token = (TObjString*) it.Next())) {
    TString var = token->GetString().Strip(TString::kBoth, ' ');

    if(var == "FiredTriggerClasses"    ){ head->SetFiredTriggerClassesIndex(indexInt); indexInt++; continue;}
    else if(var == "BunchCrossNumber"  ){ head->SetBunchCrossNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "OrbitNumber"       ){ head->SetOrbitNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "PeriodNumber"      ){ head->SetPeriodNumberIndex(indexInt); indexInt++; continue;}
    else if(var == "CentrV0M"   ) head->SetCentrIndex      (index++);
    else if(var == "CentrTRK"   ) head->SetCentrTRKIndex   (index++);
    else if(var == "CentrCL0"   ) head->SetCentrCL0Index   (index++);
    else if(var == "CentrCL1"   ) head->SetCentrCL1Index   (index++);
    else if(var == "MagField"   ) head->SetMagFieldIndex   (index++);
    else if(var == "OfflineTrigger"   ) head->SetOfflineTriggerIndex   (indexInt++);
    else if(var == "RunNumber"  ) head->SetRunNumberIndex  (indexInt++);
    else if(var == "T0Spread" ) {
      for (int i=0; i<4; i++)
        head->SetT0SpreadIndex(i, index++);
    }
    else if(var == "NumberOfESDTracks"  ) head->SetNumberOfESDTracksIndex  (indexInt++);
    else if(var.BeginsWith("cst") || var.BeginsWith("MultSelection.")) {
      cstMap[var] = index;
      if (var.BeginsWith("MultSelection."))
        fMultMap[var(14, var.Length())] = index;
      std::cout << "ADDING " << index << " " << cstMap[var] << " " << var.Data() << std::endl;
      index++;
    }
    else 
      AliFatal(Form("Invalid var [%s]", var.Data()));
  }

  vars->Delete();
  delete vars;

  head->SetMapCstVar(cstMap);
  
  fInitialized = kTRUE; 
}

void AliNanoAODSimpleSetter::SetNanoAODHeader(const AliAODEvent* event, AliNanoAODHeader* head, TString varListHeader) 
{
  if (!fInitialized)
    Init(head, varListHeader);
  
  AliAODHeader * header = dynamic_cast<AliAODHeader*>(event->GetHeader());
  if (!header) AliFatal("Not a standard AOD");

  // Set custom nano aod vars
  Double_t centrV0M=-1;
  Double_t centrTRK=-1;
  Double_t centrCL1=-1;
  Double_t centrCL0=-1;

  //2015 multiplicity selection
  AliMultSelection* MultSelection = dynamic_cast<AliMultSelection*> (event->FindListObject("MultSelection"));

  if(MultSelection){
    centrV0M = MultSelection->GetMultiplicityPercentile("V0M");
    centrCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    centrCL0 = MultSelection->GetMultiplicityPercentile("CL0");
    centrTRK = MultSelection->GetMultiplicityPercentile("TRK");
    
    for (std::map<TString,int>::iterator it = fMultMap.begin(); it != fMultMap.end(); it++) {
      double value = it->first.EndsWith(".Value") ?  MultSelection->GetEstimator(it->first(0, it->first.Length() - 6))->GetValue() : MultSelection->GetMultiplicityPercentile(it->first);
      head->SetVar(it->second, value);
    }
  }else{
    //2011 
    AliCentrality * centralityObj = header->GetCentralityP();
    centrV0M = centralityObj->GetCentralityPercentile("V0M");
    centrTRK = centralityObj->GetCentralityPercentile("TRK");
    centrCL1 = centralityObj->GetCentralityPercentile("CL1");
    centrCL0 = centralityObj->GetCentralityPercentile("CL0");
  }

  Double_t magfield = header->GetMagneticField();
  UInt_t offlineTrigger = header->GetOfflineTrigger();
  Int_t runNumber = event->GetRunNumber();
  TString firedTriggerClasses = header->GetFiredTriggerClasses();
  UShort_t bunchCrossNumber = header->GetBunchCrossNumber();
  UInt_t orbitNumber = header->GetOrbitNumber();
  UInt_t periodNumber = header->GetPeriodNumber();

  if ((head->GetFiredTriggerClassesIndex())!=-1)  head->SetFiredTriggerClasses(firedTriggerClasses);
  if ((head->GetBunchCrossNumberIndex())!=-1)     head->SetBunchCrossNumber(bunchCrossNumber);
  if ((head->GetOrbitNumberIndex())!=-1)          head->SetOrbitNumber(orbitNumber);
  if ((head->GetPeriodNumberIndex())!=-1)         head->SetPeriodNumber(periodNumber);
  if ((head->GetCentrIndex())!=-1)     head->SetVar(head->GetCentrIndex()    ,           centrV0M );
  if ((head->GetCentrTRKIndex())!=-1)  head->SetVar(head->GetCentrTRKIndex() ,           centrTRK );
  if ((head->GetCentrCL1Index())!=-1)  head->SetVar(head->GetCentrCL1Index() ,           centrCL1 );
  if ((head->GetCentrCL0Index())!=-1)  head->SetVar(head->GetCentrCL0Index() ,           centrCL0 );
  if ((head->GetMagFieldIndex())!=-1)  head->SetVar(head->GetMagFieldIndex() ,           magfield );
  if ((head->GetOfflineTriggerIndex())!=-1)  head->SetVarInt(head->GetOfflineTriggerIndex(), offlineTrigger);
  if ((head->GetRunNumberIndex())!=-1) head->SetRunNumber(runNumber);
  if (head->GetT0SpreadIndex(0) != -1)
    for (int i=0; i<4; i++)
      head->SetVar(head->GetT0SpreadIndex(i), header->GetT0spread(i));
  if ((head->GetNumberOfESDTracksIndex())!=-1) head->SetVarInt(head->GetNumberOfESDTracksIndex(), header->GetNumberOfESDTracks());
}

void AliNanoAODSimpleSetter::SetNanoAODTrack (const AliAODTrack * aodTrack, AliNanoAODTrack * nanoTrack) 
{
  // Set custom variables in the special track

  // PID
  
  // Initialize PID track mapping
  static const Bool_t bPIDFilled = AliNanoAODTrack::InitPIDIndex();
  
  // PID variables
  if (bPIDFilled)
  {
    // Get the PID info
    static AliPIDResponse* pidResponse = 0;
    if (!pidResponse) {
      AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
      AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
      pidResponse = inputHandler->GetPIDResponse();
    }
    
    if (!pidResponse)
      AliFatal("PID response not available but fields requested");

    for (Int_t r = AliNanoAODTrack::kSigmaTPC; r<AliNanoAODTrack::kLAST; r++) {
      for (Int_t p = 0; p<AliPID::kSPECIESC; p++) {
        Int_t index = AliNanoAODTrack::GetPIDIndex((AliNanoAODTrack::ENanoPIDResponse) r, (AliPID::EParticleType) p);
        if (index == -1)
          continue;
        Double_t value = 0;
        if (r == AliNanoAODTrack::kSigmaTPC)
          value = pidResponse->NumberOfSigmasTPC(aodTrack, (AliPID::EParticleType) p);
        else if (r == AliNanoAODTrack::kSigmaTOF)
          value = pidResponse->NumberOfSigmasTOF(aodTrack, (AliPID::EParticleType) p);

        nanoTrack->SetVar(index, value);
      }
    }
  }

  // additional fields which are to be moved to AliNanoAODTrack
  // TPC clusters
  static const Int_t cstTPCClusterInfo21 = AliNanoAODTrackMapping::GetInstance()->GetVarIndex("cstTPCClusterInfo21");
  if (cstTPCClusterInfo21 != -1)
    nanoTrack->SetVar(cstTPCClusterInfo21, aodTrack->GetTPCClusterInfo(2,1));
}
