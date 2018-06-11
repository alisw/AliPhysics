/*
 * AliFemtoDreamTrack.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */
#include "AliAnalysisManager.h"
#include "AliAODPid.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliFemtoDreamTrack.h"
#include "AliLog.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <iostream>
ClassImp(AliFemtoDreamTrack)
AliFemtoDreamTrack::AliFemtoDreamTrack()
:AliFemtoDreamBasePart()
,fPIDResponse(0)
,fstatusTPC(AliPIDResponse::kDetNoParams)
,fstatusTOF(AliPIDResponse::kDetNoParams)
,fFilterMap(0)
,fdcaXY(-99)
,fdcaZ(-99)
,fdcaXYProp(-99)
,fdcaZProp(-99)
,fNClsTPC(0)
,fTPCCrossedRows(0)
,fRatioCR(0)
,fnoSharedClst(0)
,fTPCClsS(0)
,fChi2(0.f)
,fSharedClsITSLayer(0)
,fHasSharedClsITSLayer(false)
,fdEdxTPC(0)
,fbetaTOF(0)
,fHasITSHit(false)
,fITSHit(0)
,fTOFTiming(false)
,fTPCRefit(false)
,fESDTrack(0)
,fAODTrack(0)
,fAODGlobalTrack(0)
,fESDStatus(0)
{
  for (int i=0;i<5;++i) {
    fnSigmaTPC[i]=0;
    fnSigmaTOF[i]=0;
  }
}

AliFemtoDreamTrack::~AliFemtoDreamTrack() {
  // TODO Auto-generated destructor stub
}

void AliFemtoDreamTrack::SetTrack(AliAODTrack *track) {
  this->Reset();
  fAODTrack=track;
  int trackID=fAODTrack->GetID();
  if (trackID<0) {
    if(!fGTI){
      AliFatal("AliFemtoSPTrack::SetTrack No fGTI Set");
      fAODGlobalTrack=NULL;
    }else if(-trackID-1 >= fTrackBufferSize){
//      AliFatal("Buffer Size too small");
      this->fIsSet=false;
      fIsReset=false;
      fAODGlobalTrack=NULL;
    }else if(!CheckGlobalTrack(trackID)){
      fAODGlobalTrack=NULL;
    }else{
      fAODGlobalTrack=fGTI[-trackID-1];
    }
  }else{
    fAODGlobalTrack=track;
  }
  fIsReset=false;
  if (fAODGlobalTrack&&fAODTrack) {
    this->SetIDTracks(fAODGlobalTrack->GetID());
    this->SetAODTrackingInformation();
    this->SetPIDInformation();
    if (fIsMC) {
      this->SetMCInformation();
    }
  } else {
    this->fIsSet=false;
  }
}

void AliFemtoDreamTrack::SetTrack(AliESDtrack *track) {
  this->Reset();
  fESDTrack=track;
  if (fESDTrack) {
    this->fIsReset=false;
    int trackID=fESDTrack->GetID();
    if (trackID < 0) {
      AliWarning("Negative ID for ESD tracks");
    }
    this->SetIDTracks(trackID);
    this->SetESDTrackingInformation();
//    this->SetPIDInformation();
//    if (fIsMC) {
//      this->SetMCInformation();
//    }
  }
}

void AliFemtoDreamTrack::SetESDTrackingInformation() {
  fESDStatus=fESDTrack->GetStatus();
  fESDnClusterITS=fESDTrack->GetITSclusters(0);
  fESDnClustersTPC = fESDTrack->GetTPCclusters(0);
  fESDnCrossedRowsTPC = fESDTrack->GetTPCCrossedRows();
  if (esdTrack->GetTPCNclsF()>0) {
    ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / esdTrack->GetTPCNclsF();
  }


}
void AliFemtoDreamTrack::SetAODTrackingInformation() {
  this->fFilterMap=fAODTrack->GetFilterMap();
  this->SetEta(fAODTrack->Eta());
  this->SetPhi(fAODTrack->Phi());
  this->SetTheta(fAODTrack->Theta());
  this->SetCharge(fAODTrack->Charge());
  this->SetMomentum(fAODTrack->Px(), fAODTrack->Py(), fAODTrack->Pz());
  this->SetMomTPC(fAODGlobalTrack->GetTPCmomentum());
  this->SetPt(fAODTrack->Pt());
  this->fdcaXY=fAODTrack->DCA();
  this->fdcaZ=fAODTrack->ZAtDCA();
  this->fChi2=fAODTrack->Chi2perNDF();
  double dcaVals[2] = {-99., -99.};
  double covar[3]={0.,0.,0.};
  AliAODTrack copy(*fAODGlobalTrack);
  if (copy.PropagateToDCA(copy.GetAODEvent()->GetPrimaryVertex(),
                          copy.GetAODEvent()->GetMagneticField(),
                          10, dcaVals, covar))
  {
    this->fdcaXYProp = dcaVals[0];
    this->fdcaZProp  = dcaVals[1];
  }else{
    this->fdcaXYProp = -99;
    this->fdcaZProp  = -99;
  }
  //loop over the 6 ITS Layrs and check for a hit!
  for (int i=0;i<6;++i) {
    fITSHit.push_back(fAODGlobalTrack->HasPointOnITSLayer(i));
    if (fAODGlobalTrack->HasPointOnITSLayer(i)) {
      this->fHasITSHit=true;
    }
  }
  if (fAODTrack->IsOn(AliAODTrack::kTPCrefit)) {
    fTPCRefit=true;
  }
  if (fAODGlobalTrack->GetTOFBunchCrossing()==0) {
    this->fTOFTiming=true;
  } else {
    this->fTOFTiming=false;
  }

  this->fNClsTPC=fAODTrack->GetTPCNcls();
  //This method was inherited from H. Beck analysis
  // In the documents
  // https://alisoft.cern.ch/AliRoot/trunk/TPC/doc/Definitions/Definitions.pdf
  // TPC people describe the cut strategy for the TPC. It is explicitly
  // stated that a cut on the number of crossed rows and a cut on the
  // number of crossed rows over findable clusters is recommended to
  // remove fakes. In the pdf a cut value of .83 on the ratio
  // is stated, no value for the number of crossed rows. Looking at the
  // AliESDTrackCuts.cxx one sees that exactly this cut is used with
  // 0.8 on the ratio and 70 on the crossed rows.

  // Checked the filter task and AliAODTrack and AliESDtrack and
  // AliESDtrackCuts and the Definitions.pdf:
  // The function to get the findable clusters is GetTPCNclsF()

  // For the number fo crossed rows for ESD tracks, the function
  // GetTPCCrossedRows() usually is used. Looking at the AliESDtrack.cxx
  // one sees that it's just an alias (with additional caching) for
  // GetTPCClusterInfo(2, 1); The identical function exists in the
  // AliAODTrack.cxx
  this->fTPCCrossedRows=fAODTrack->GetTPCClusterInfo(2, 1);
  if (!(fAODTrack->GetTPCNclsF()>0)) {
    this->fRatioCR=0.;
  } else {
    this->fRatioCR=
        fAODTrack->GetTPCClusterInfo(2, 1)/float(fAODTrack->GetTPCNclsF());
  }
  const TBits sharedMap=fAODTrack->GetTPCSharedMap();
  if ((sharedMap.CountBits()) >= 1) {
    // Bad Track, has too many shared clusters!
    this->fnoSharedClst=false;
  }else{
    this->fnoSharedClst=true;
  }
  for (int i=0;i<6;++i) {
    fSharedClsITSLayer.push_back(fAODTrack->HasSharedPointOnITSLayer(i));
    if (fAODTrack->HasSharedPointOnITSLayer(i)) {
      fHasSharedClsITSLayer=true;
    }
  }
  this->fTPCClsS=fAODTrack->GetTPCnclsS();
  if (fIsMC) {
    SetPhiAtRadii();
  }
}
void AliFemtoDreamTrack::SetPhiAtRadii() {
  float TPCradii[9] = {85.,105.,125.,145.,165.,185.,205.,225.,245.};
  float phi0=GetPhi().at(0);
  float pt=GetPt();
  float chg=GetCharge().at(0);
  float bfield=fAODTrack->GetAODEvent()->GetMagneticField();
  this->SetEvtNumber(fAODTrack->GetAODEvent()->GetRunNumber());
  std::vector<float> phiatRadius;
  for(int radius=0;radius<9;radius++)
  {
    phiatRadius.push_back(
        phi0 + TMath::ASin(0.1*chg*bfield*0.3*TPCradii[radius]*0.01/(2.*pt)));
  }
  fPhiAtRadius.push_back(phiatRadius);
  return;
}
void AliFemtoDreamTrack::SetPIDInformation() {
  AliPID::EParticleType particleID[5] = {AliPID::kElectron,AliPID::kMuon,
      AliPID::kPion,AliPID::kKaon,AliPID::kProton};
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if (man) {
    AliInputEventHandler* inputHandler=
        (AliInputEventHandler*)(man->GetInputEventHandler());
    if (inputHandler) {
      fPIDResponse=inputHandler->GetPIDResponse();
      if (!fPIDResponse) {
        AliFatal("No PID Response, did you run your PID Task?");
      }
    } else {
      AliFatal("No Input Handler");
    }
  }else{
    AliFatal("No PID Response");
  }
  AliPIDResponse::EDetPidStatus statusTPC =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,fAODGlobalTrack);
  AliPIDResponse::EDetPidStatus statusTOF =
      fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fAODGlobalTrack);
  this->fstatusTPC=statusTPC;
  this->fstatusTOF=statusTOF;
  this->fdEdxTPC=fAODGlobalTrack->GetTPCsignal();
  this->fbetaTOF=GetBeta(fAODGlobalTrack);
  for (int i=0;i<5;++i) {
    if(statusTPC == AliPIDResponse::kDetPidOk){
      (this->fnSigmaTPC)[i] =
          fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,fAODGlobalTrack,particleID[i]);
    }else{
      (this->fnSigmaTPC)[i] = -999.;
    }
    if(statusTOF == AliPIDResponse::kDetPidOk){
      (this->fnSigmaTOF)[i] =
          fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,fAODGlobalTrack,particleID[i]);
    }else{
      (this->fnSigmaTOF)[i] = -999.;
    }
  }
}

void AliFemtoDreamTrack::SetMCInformation() {
  //Set the phi at radii at the TPC information
  //      SetPhiStar(track,fphiAtRadius);
  TClonesArray *mcarray =
      dynamic_cast<TClonesArray*>(fAODGlobalTrack->GetAODEvent()
          ->FindListObject(AliAODMCParticle::StdBranchName()));
  if (!mcarray) {
    AliError("SPTrack: MC Array not found");
  }
  if (fAODGlobalTrack->GetLabel()>0) {
    AliAODMCParticle * mcPart = (AliAODMCParticle*)mcarray->At(fAODGlobalTrack->GetLabel());
    if (!(mcPart)) {
      this->fIsSet=false;
    } else {
      this->SetMCPhi(mcPart->Phi());
      this->SetMCTheta(mcPart->Theta());
      this->SetMCPDGCode(mcPart->PdgCode());
      this->SetMCPt(mcPart->Pt());
      this->SetMCMomentum(mcPart->Px(),mcPart->Py(),mcPart->Pz());

      //check for secondary and set origin and mother
      if (mcPart->IsPhysicalPrimary() && !mcPart->IsSecondaryFromWeakDecay()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kPhysPrimary);
      } else if(mcPart->IsSecondaryFromWeakDecay() && !mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kWeak);
        this->SetPDGMotherWeak(((AliAODMCParticle*)mcarray->At(mcPart->GetMother()))->PdgCode());
      } else if (mcPart->IsSecondaryFromMaterial()) {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kMaterial);
      } else {
        this->SetParticleOrigin(AliFemtoDreamBasePart::kUnknown);
      }
    }
  } else {
    this->fIsSet =false; //if we don't have MC Information, don't use that track
  }
}

float AliFemtoDreamTrack::GetBeta(AliAODTrack *track) {
  float beta = -999;
  double integratedTimes[9] = {-1.0,-1.0,-1.0,-1.0,-1.0, -1.0, -1.0, -1.0, -1.0};

  track->GetIntegratedTimes(integratedTimes);

  const float c = 2.99792457999999984e-02;
  float p = track->P();
  float l = integratedTimes[0]*c;

  float trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);

  float timeTOF = track->GetTOFsignal()- trackT0;
  if(timeTOF > 0){
    beta  = l/timeTOF/c;
  }
  return beta;
}

bool AliFemtoDreamTrack::CheckGlobalTrack(const Int_t TrackID) {
  //This method was inherited from H. Beck analysis
  //Checks if to the corresponding track a global track exists
  //This is especially useful if one has TPC only tracks and needs the PID information
  bool isGlobal = true;
  if (TMath::Abs(TrackID)<fTrackBufferSize) {
    if (!(fGTI[-TrackID-1])) {
      isGlobal = false;
    }
  }
  return isGlobal;
}
void AliFemtoDreamTrack::Reset() {
  if (!fIsReset) {
    fstatusTPC=AliPIDResponse::kDetNoParams;
    fstatusTOF=AliPIDResponse::kDetNoParams;
    fFilterMap=0;
    fdcaXY=-99;
    fdcaZ=-99;
    fdcaXYProp=-99;
    fdcaZProp=-99;
    fNClsTPC=0;
    fTPCCrossedRows=0;
    fRatioCR=0;
    fnoSharedClst=0;
    fTPCClsS=0;
    fSharedClsITSLayer.clear();
    fHasSharedClsITSLayer=false;
    fdEdxTPC=-999;
    fbetaTOF=1.1;
    for (int i=0;i<5;++i) {
      fnSigmaTPC[i]=99;
      fnSigmaTOF[i]=99;
    }
    fHasITSHit=false;
    fITSHit.clear();
    fTOFTiming=false;
    fTPCRefit=false;
    fP.SetXYZ(0,0,0);
    fMCP.SetXYZ(0,0,0);
    fPt=0;
    fMCPt=0;
    fMCPt=0;
    fP_TPC=0;
    fEta.clear();
    fTheta.clear();
    fMCTheta.clear();
    fPhi.clear();
    fPhiAtRadius.clear();
    fMCPhi.clear();
    fIDTracks.clear();
    fCharge.clear();
    fCPA=0;
    fOrigin=AliFemtoDreamBasePart::kUnknown;
    //we don't want to reset the fPDGCode
    fMCPDGCode=0;
    fPDGMotherWeak=0;
    //we don't want to reset isMC
    fUse=false;
    fIsSet=true;
    fIsReset=true;
  }
}
