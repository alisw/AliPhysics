/*
 * AliLightNTrack.cxx
 *
 *  Created on: Nov 14, 2017
 *      Author: gu74req
 */
#include "AliLightNTrack.h"
ClassImp(AliLightNTrack)
AliLightNTrack::AliLightNTrack()
:AliLightNBasePart()
,fPIDResponse(0)
,fstatusTPC(AliPIDResponse::kDetNoParams)
,fstatusTOF(AliPIDResponse::kDetNoParams)
,fFilterMap(0)
,fdcaXY(-99)
,fdcaZ(-99)
,fdcaXYProp(-99)
,fdcaZProp(-99)
,fNClsTPC(0)
,fNClsITS(0)
,fTPCCrossedRows(0)
,fRatioCR(0)
,fnoSharedClst(0)
,fTPCClsS(0)
,fSharedClsITSLayer(0)
,fHasSharedClsITSLayer(false)
,fdEdxTPC(0)
,fbetaTOF(0)
,fmass2sq(0)
,fHasITSHit(false)
,fITSHit(0)
,fTOFTiming(false)
,fTPCRefit(false)
,fTrack(0)
,fGlobalTrack(0)
{
    for (int i=0;i<6;++i) {
        fnSigmaTPC[i]=0;
        fnSigmaTOF[i]=0;
        fnSigmaITS[i]=0;
    }
}

AliLightNTrack::~AliLightNTrack() {
    // TODO Auto-generated destructor stub
}

void AliLightNTrack::SetTrack(AliAODTrack *track) {
    this->Reset();
    fTrack=track;
    int trackID=fTrack->GetID();
    if (trackID<0) {
        if(!fGTI){
            AliFatal("AliFemtoSPTrack::SetTrack No fGTI Set");
            fGlobalTrack=NULL;
        }else if(-trackID-1 >= fTrackBufferSize){
            AliFatal("Buffer Size too small");
            fGlobalTrack=NULL;
        }else if(!CheckGlobalTrack(trackID)){
            fGlobalTrack=NULL;
        }else{
            fGlobalTrack=fGTI[-trackID-1];
        }
    }else{
        fGlobalTrack=track;
    }
    fIsReset=false;
    if (fGlobalTrack&&fTrack) {
        this->SetIDTracks(fGlobalTrack->GetID());
        this->SetTrackingInformation();
        this->SetPIDInformation();
        if (fIsMC) {
            this->SetMCInformation();
        }
    } else {
        this->fIsSet=false;
    }
}

void AliLightNTrack::SetTrackingInformation() {
    this->fFilterMap=fTrack->GetFilterMap();
    this->SetTrackChi2perNDF(fTrack->Chi2perNDF());
    this->SetEta(fTrack->Eta());
    this->SetPhi(fTrack->Phi());
    this->SetTheta(fTrack->Theta());
    this->SetCharge(fTrack->Charge());
    this->SetMomentum(fTrack->Px(), fTrack->Py(), fTrack->Pz());
    this->SetP(fTrack->P());
    this->SetMomTPC(fGlobalTrack->GetTPCmomentum());
    this->SetPt(fTrack->Pt());
    this->fdcaXY=fTrack->DCA();
    this->fdcaZ=fTrack->ZAtDCA();
    double dcaVals[2] = {-99., -99.};
    double covar[3]={0.,0.,0.};
    AliAODTrack copy(*fGlobalTrack);
    const AliVVertex *AODeventVtx = copy.GetAODEvent()->GetPrimaryVertex();
    Double_t AODeventMagneticF = copy.GetAODEvent()->GetMagneticField();
    if (copy.PropagateToDCA(AODeventVtx, AODeventMagneticF,10, dcaVals, covar)){
        this->fdcaXYProp = dcaVals[0];
        this->fdcaZProp  = dcaVals[1];
    }else{
        this->fdcaXYProp = -99;
        this->fdcaZProp  = -99;
    }
    //loop over the 6 ITS Layrs and check for a hit!
    for (int i=0;i<6;++i) {
        fITSHit.push_back(fGlobalTrack->HasPointOnITSLayer(i));
        if (fGlobalTrack->HasPointOnITSLayer(i)) {
            this->fHasITSHit=true;
        }
    }
    if (fTrack->IsOn(AliAODTrack::kTPCrefit)) {
        fTPCRefit=true;
    }
    if (fGlobalTrack->GetTOFBunchCrossing()==0) {
        this->fTOFTiming=true;
    } else {
        this->fTOFTiming=false;
    }
    
    this->fNClsTPC=fTrack->GetTPCNcls();
    this->fNClsITS=fTrack->GetITSNcls();
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
    this->fTPCCrossedRows=fTrack->GetTPCClusterInfo(2, 1);
    if (!fTrack->GetTPCNclsF()) {
        this->fRatioCR=0.;
    } else {
        this->fRatioCR=
        fTrack->GetTPCClusterInfo(2, 1)/double(fTrack->GetTPCNclsF());
    }
    const TBits sharedMap=fTrack->GetTPCSharedMap();
    if ((sharedMap.CountBits()) >= 1) {
        // Bad Track, has too many shared clusters!
        this->fnoSharedClst=false;
    }else{
        this->fnoSharedClst=true;
    }
    for (int i=0;i<6;++i) {
        fSharedClsITSLayer.push_back(fTrack->HasSharedPointOnITSLayer(i));
        if (fTrack->HasSharedPointOnITSLayer(i)) {
            fHasSharedClsITSLayer=true;
        }
    }
    this->fTPCClsS=fTrack->GetTPCnclsS();
}

void AliLightNTrack::SetPIDInformation() {
    AliPID::EParticleType particleID[6] = {AliPID::kElectron,AliPID::kMuon,
        AliPID::kPion,AliPID::kKaon,AliPID::kProton,AliPID::kDeuteron};
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
    fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,fGlobalTrack);
    AliPIDResponse::EDetPidStatus statusTOF =
    fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,fGlobalTrack);
    AliPIDResponse::EDetPidStatus statusITS =
    fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS,fGlobalTrack);
    this->fstatusTPC=statusTPC;
    this->fstatusTOF=statusTOF;
    this->fstatusITS=statusITS;
    this->fdEdxTPC=fGlobalTrack->GetTPCsignal();
    this->fbetaTOF=GetBeta(fGlobalTrack);
    this->fmass2sq=GetMass2sq(fbetaTOF,fGlobalTrack);
    this->fdEdxITS=fGlobalTrack->GetITSsignal();
    for (int i=0;i<6;++i) {
        
        if(statusTPC == AliPIDResponse::kDetPidOk){
            (this->fnSigmaTPC)[i] =
            fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC,fGlobalTrack,particleID[i]);
        }else{
            (this->fnSigmaTPC)[i] = -999.;
        }
        if(statusTOF == AliPIDResponse::kDetPidOk){
            (this->fnSigmaTOF)[i] =
            fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF,fGlobalTrack,particleID[i]);
        }else{
            (this->fnSigmaTOF)[i] = -999.;
        }
        if(statusITS == AliPIDResponse::kDetPidOk){
            (this->fnSigmaITS)[i] =
            fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS,fGlobalTrack,particleID[i]);
        }else{
            (this->fnSigmaITS)[i] = -999.;
        }
    }
}

void AliLightNTrack::SetMCInformation() {
    //Set the phi at radii at the TPC information
    //      SetPhiStar(track,fphiAtRadius);
    TClonesArray *mcarray =
    dynamic_cast<TClonesArray*>(fGlobalTrack->GetAODEvent()
                                ->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcarray) {
        AliError("SPTrack: MC Array not found");
    }
    if (fGlobalTrack->GetLabel()>=0) {
        AliAODMCParticle * mcPart = (AliAODMCParticle*)mcarray->At(fGlobalTrack->GetLabel());;
        if (!(mcPart)) {
            this->fIsSet=false;
        } else {
            this->SetMCPhi(mcPart->Phi());
            this->SetMCTheta(mcPart->Theta());
            this->SetMCPDGCode(mcPart->GetPdgCode());
            this->SetMCPt(mcPart->Pt());
            this->SetMCMomentum(mcPart->Px(),mcPart->Py(),mcPart->Pz());
            
            //check for secondary and set origin and mother
            if (mcPart->IsPhysicalPrimary() && !mcPart->IsSecondaryFromWeakDecay()) {
                this->SetParticleOrigin(AliLightNBasePart::kPhysPrimary);
            } else if(mcPart->IsSecondaryFromWeakDecay() && !mcPart->IsSecondaryFromMaterial()) {
                this->SetParticleOrigin(AliLightNBasePart::kWeak);
                this->SetPDGMotherWeak(((AliAODMCParticle*)mcarray->At(mcPart->GetMother()))->GetPdgCode());
            } else if (mcPart->IsSecondaryFromMaterial()) {
                this->SetParticleOrigin(AliLightNBasePart::kMaterial);
            } else {
                this->SetParticleOrigin(AliLightNBasePart::kUnknown);
            }
        }
    } else {
        this->fIsSet =false; //if we don't have MC Information, don't use that track
    }
}

float AliLightNTrack::GetBeta(AliAODTrack *track) {
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

Float_t AliLightNTrack::GetMass2sq(Double_t beta,AliAODTrack *track){	
    
    Float_t p = track->P();
    Float_t mass2sq = -999;
    if(!(beta<0)){
        mass2sq = ((1/(beta*beta))-1)*(p*p);
    }
    return mass2sq;
}




bool AliLightNTrack::CheckGlobalTrack(const Int_t TrackID) {
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
void AliLightNTrack::Reset() {
    if (!fIsReset) {
        fstatusTPC=AliPIDResponse::kDetNoParams;
        fstatusTOF=AliPIDResponse::kDetNoParams;
        fFilterMap=0;
        fdcaXY=-99;
        fdcaZ=-99;
        fdcaXYProp=-99;
        fdcaZProp=-99;
        fNClsTPC=0;
        fNClsITS=0;
        fTPCCrossedRows=0;
        fRatioCR=0;
        fnoSharedClst=0;
        fTPCClsS=0;
        fSharedClsITSLayer.clear();
        fHasSharedClsITSLayer=false;
        fdEdxTPC=-999;
        fbetaTOF=1.1;
        fmass2sq=-999;
        for (int i=0;i<6;++i) {
            fnSigmaTPC[i]=99;
            fnSigmaTOF[i]=99;
            fnSigmaITS[i]=99;
        }
        fHasITSHit=false;
        fITSHit.clear();
        fTOFTiming=false;
        fTPCRefit=false;
        fP.SetXYZ(0,0,0);
        fMCP.SetXYZ(0,0,0);
        fTrackChi2perNDF=0;
        fPt=0;
        fMCPt=0;
        fMCPt=0;
        fP_TPC=0;
        fEta.clear();
        fTheta.clear();
        fMCTheta.clear();
        fPhi.clear();
        fMCPhi.clear();
        fIDTracks.clear();
        fCharge.clear();
        fCPA=0;
        fOrigin=AliLightNBasePart::kUnknown;
        //we don't want to reset the fPDGCode
        fMCPDGCode=0;
        fPDGMotherWeak=0;
        //we don't want to reset isMC
        fUse=false;
        fIsSet=true;
        fIsReset=true;
    }
}
