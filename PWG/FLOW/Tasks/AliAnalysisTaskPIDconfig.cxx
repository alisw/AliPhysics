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

/* $Id: AliAnalysisTaskPIDconfig.cxx 43811 2014-10-11 Naghmeh Mohammadi $ */

#include "TChain.h"
#include "TTree.h"
#include "TList.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TROOT.h"
#include "stdio.h"
#include "TCutG.h"


#include "AliTHn.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliPIDCombined.h"
#include "AliAnalysisTask.h"
#include "AliAODHandler.h"
#include <AliInputEventHandler.h>
#include <AliVEventHandler.h>
#include <AliVParticle.h>
#include <AliVTrack.h>
#include <AliTPCPIDResponse.h>
#include <AliTOFPIDResponse.h>
#include "AliAnalysisTaskPIDconfig.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODPid.h"
#include "AliPhysicsSelection.h"
#include "AliCentralitySelectionTask.h"
#include "AliCentrality.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliAODVZERO.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskPIDconfig)
//ClassImp()
//___________________________________________________________________
AliAnalysisTaskPIDconfig::AliAnalysisTaskPIDconfig():
AliAnalysisTaskSE(),
fVevent(0),
fESD(0),
fAOD(0),
fPIDResponse(0),
fTriggerSelection(0),
fCentralityPercentileMin(0.),
fCentralityPercentileMax(5.),
fFilterBit(1),
fDCAxyCut(-1),
fDCAzCut(-1),
fData2011(kFALSE),
fTriggerMB(kTRUE),
fTriggerCentral(kFALSE),
fUseCentrality(kTRUE),
fCutTPCmultiplicityOutliersAOD(kFALSE),
fPIDcuts(kFALSE),
fCentralityEstimator("V0M"),
fContoursFile(0),
fCutContourList(0),
fListQA(0x0),
fListQAtpctof(0x0),
fListQAInfo(0x0),
fhistCentralityPass(0),
fNoEvents(0),
fpVtxZ(0),
fhistDCABefore(0),
fhistDCAAfter(0),
fhistPhiDistBefore(0),
fhistPhiDistAfter(0),
fhistEtaDistBefore(0),
fhistEtaDistAfter(0),
fTPCvsGlobalMultBeforeOutliers(0),
fTPCvsGlobalMultAfterOutliers(0),
fTPCvsGlobalMultAfter(0),
fHistBetavsPTOFbeforePID(0),
fHistdEdxvsPTPCbeforePID(0),
fhistNsigmaP(0),
fhistTPCnSigmavsP(0),
fHistBetavsPTOFafterPID(0),
fHistdEdxvsPTPCafterPID(0),
fHistBetavsPTOFafterPID_2(0),
fHistdEdxvsPTPCafterPID_2(0),
fHistBetavsPTOFafterPIDTPCTOF(0),
fHistdEdxvsPTPCafterPIDTPCTOF(0),
fHistBetavsPTOFafterPIDTPConly(0),
fHistdEdxvsPTPCafterPIDTPConly(0),
fHistPion_BetavsPTOFafterPIDTPCTOF(0),
fHistPion_dEdxvsPTPCafterPIDTPCTOF(0),
fHistKaon_BetavsPTOFafterPIDTPCTOF(0),
fHistKaon_dEdxvsPTPCafterPIDTPCTOF(0),
fHistProton_BetavsPTOFafterPIDTPCTOF(0),
fHistProton_dEdxvsPTPCafterPIDTPCTOF(0),
fhistPionEtaDistAfter(0),
fhistKaonEtaDistAfter(0),
fhistProtonEtaDistAfter(0)
 //fSparseSpecies(0),
//fvalueSpecies(0)
{
    for(int i=0;i<150;i++){
        fCutContour[i]= NULL;
        fCutGraph[i]=NULL;
    }
    //Low momentum nsigma cuts based on Purity>0.7 with TPC info only.
    
    for(int i=0;i<6;i++){
        
        fLowPtPIDTPCnsigLow_Pion[i] = 0;
        fLowPtPIDTPCnsigLow_Kaon[i] = 0;
        fLowPtPIDTPCnsigHigh_Pion[i] =0;
        fLowPtPIDTPCnsigHigh_Kaon[i] =0;
    }

}


//___________________________________________________________________

AliAnalysisTaskPIDconfig::AliAnalysisTaskPIDconfig(const char *name):
AliAnalysisTaskSE(name),
fVevent(0),
fESD(0),
fAOD(0),
fPIDResponse(0),
fTriggerSelection(0),
fCentralityPercentileMin(0.),
fCentralityPercentileMax(5.),
fFilterBit(1),
fDCAxyCut(-1),
fDCAzCut(-1),
fData2011(kFALSE),
fTriggerMB(kTRUE),
fTriggerCentral(kFALSE),
fUseCentrality(kTRUE),
fCutTPCmultiplicityOutliersAOD(kFALSE),
fPIDcuts(kFALSE),
fCentralityEstimator("V0M"),
fContoursFile(0),
fCutContourList(0),
fListQA(0x0),
fListQAtpctof(0x0),
fListQAInfo(0x0),
fhistCentralityPass(0),
fNoEvents(0),
fpVtxZ(0),
fhistDCABefore(0),
fhistDCAAfter(0),
fhistPhiDistBefore(0),
fhistPhiDistAfter(0),
fhistEtaDistBefore(0),
fhistEtaDistAfter(0),
fTPCvsGlobalMultBeforeOutliers(0),
fTPCvsGlobalMultAfterOutliers(0),
fTPCvsGlobalMultAfter(0),
fHistBetavsPTOFbeforePID(0),
fHistdEdxvsPTPCbeforePID(0),
fhistNsigmaP(0),
fhistTPCnSigmavsP(0),
fHistBetavsPTOFafterPID(0),
fHistdEdxvsPTPCafterPID(0),
fHistBetavsPTOFafterPID_2(0),
fHistdEdxvsPTPCafterPID_2(0),
fHistBetavsPTOFafterPIDTPCTOF(0),
fHistdEdxvsPTPCafterPIDTPCTOF(0),
fHistBetavsPTOFafterPIDTPConly(0),
fHistdEdxvsPTPCafterPIDTPConly(0),
fHistPion_BetavsPTOFafterPIDTPCTOF(0),
fHistPion_dEdxvsPTPCafterPIDTPCTOF(0),
fHistKaon_BetavsPTOFafterPIDTPCTOF(0),
fHistKaon_dEdxvsPTPCafterPIDTPCTOF(0),
fHistProton_BetavsPTOFafterPIDTPCTOF(0),
fHistProton_dEdxvsPTPCafterPIDTPCTOF(0),
fhistPionEtaDistAfter(0),
fhistKaonEtaDistAfter(0),
fhistProtonEtaDistAfter(0)
//fSparseSpecies(0),
//fvalueSpecies(0)root
{
    //Default Constructor
    for(int i=0;i<150;i++){
        fCutContour[i]= NULL;
        fCutGraph[i]=NULL;
    }

    //fCutContour[150]=NULL;
    //fCutGraph[150]=NULL;
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
}

//_____________________________________________________________________
AliAnalysisTaskPIDconfig::~AliAnalysisTaskPIDconfig()
{
    //Destructor
    
    fContoursFile->Close();
    for(int i=0;i<150;i++){
        delete fCutContour[i];
        delete fCutGraph[i];
    }
    
    
    //   delete fPID;
    //  delete fPIDqa;
    //   delete fTrackCuts;
    // delete fSparseSpecies;
    //delete []fvalueSpecies;
    
    
}
//______________________________________________________________________
void AliAnalysisTaskPIDconfig::UserCreateOutputObjects()
{
    //
    // Create the output QA objects
    //
    
    AliLog::SetClassDebugLevel("AliAnalysisTaskPIDconfig",10);
    
    //input hander
    AliAnalysisManager *mgr=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(mgr->GetInputEventHandler());
    if (!inputHandler) {
        AliFatal("Input handler needed");
        return;         // to shut up coverity
    }
    
    //pid response object
    fPIDResponse=inputHandler->GetPIDResponse();
    if (!fPIDResponse) AliError("PIDResponse object was not created");
    
    if(fPIDcuts){ GetPIDContours(); cout<<"********** PID cut contours retrieved **********"<<endl;}
    //
    fListQA=new TList;
    fListQA->SetOwner();
    
    fListQAtpctof=new TList;
    fListQAtpctof->SetOwner();
    fListQAtpctof->SetName("PID_TPC_TOF");
    
    fListQAInfo=new TList;
    fListQAInfo->SetOwner();
    fListQAInfo->SetName("Event_Track_Info");
    
    fListQA->Add(fListQAtpctof);
    fListQA->Add(fListQAInfo);
    
    SetupTPCTOFqa();
    SetupEventInfo();
    
    PostData(1,fListQA);
}
//______________________________________________________________________
void AliAnalysisTaskPIDconfig::UserExec(Option_t*){
    //Main loop
    //Called for each event
    
    // create pointer to event
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    
    
    
    if(!(fESD || fAOD)){
        printf("ERROR: fESD & fAOD not available\n");
        return;
    }

    Int_t ntracks=fAOD->GetNumberOfTracks();

    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVevent not available\n");
        return;
    }
    
    Bool_t pass = kFALSE;
    CheckCentrality(fVevent,pass);
    
    if(!pass){ return;}
    
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    
    Double_t pVtxZ = -999;
    pVtxZ = pVtx->GetZ();
    
    if(TMath::Abs(pVtxZ)>10) return;
    
    TH1F *hNoEvents = (TH1F*)fListQAInfo->At(1);
    TH1F *histpVtxZ = (TH1F*)fListQAInfo->At(2);
    
    if(hNoEvents) hNoEvents->Fill(0);
    if(histpVtxZ) histpVtxZ->Fill(pVtxZ);
    
    if(ntracks<2) return;
    
    // if(!pass) return;
    
    TH2F *HistTPCvsGlobalMultBeforeOutliers = (TH2F*)fListQAInfo->At(3);
    
    TH2F *HistTPCvsGlobalMultAfterOutliers = (TH2F*)fListQAInfo->At(4);
    
    Float_t multTPC(0.); // tpc mult estimate
    Float_t multGlobal(0.); // global multiplicity
    
    const Int_t nGoodTracks = fVevent->GetNumberOfTracks();
    if(!fData2011) { // cut on outliers
        
        for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { // fill tpc mult
            AliAODTrack* AODtrack =dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
            if (!AODtrack) continue;
            if (!(AODtrack->TestFilterBit(1))) continue;
            if ((AODtrack->Pt() < .2) || (AODtrack->Pt() > 5.0) || (TMath::Abs(AODtrack->Eta()) > .8) || (AODtrack->GetTPCNcls() < 70)  || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.2)) continue;
            multTPC++;
        }//track loop
        
        for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { // fill global mult
            AliAODTrack *AODtrack=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
            if (!AODtrack) continue;
            if (!(AODtrack->TestFilterBit(16))) continue;
            if ((AODtrack->Pt() < .2) || (AODtrack->Pt() > 5.0) || (TMath::Abs(AODtrack->Eta()) > .8) || (AODtrack->GetTPCNcls() < 70) || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.1)) continue;
            Double_t b[2] = {-99., -99.};
            Double_t bCov[3] = {-99., -99., -99.};
            AliAODTrack copy(*AODtrack);
            if (!(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., b, bCov))) continue;
            if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
            multGlobal++;
        } //track loop
        
        HistTPCvsGlobalMultBeforeOutliers->Fill(multGlobal,multTPC);
        
        if(multTPC < (-40.3+1.22*multGlobal) || multTPC > (32.1+1.59*multGlobal)){ pass = kFALSE;}
        
        if(!pass) return;
        HistTPCvsGlobalMultAfterOutliers->Fill(multGlobal,multTPC);
        
    }
    
    
    if(fData2011) { // cut on outliers
        //Float_t multTPC(0.); // tpc mult estimate
        //Float_t multGlob(0.); // global multiplicity
        for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { // fill tpc mult
            AliAODTrack *AODtrack=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
            if (!AODtrack) continue;
            if (!(AODtrack->TestFilterBit(1))) continue;
            if ((AODtrack->Pt() < .2) || (AODtrack->Pt() > 5.0) || (TMath::Abs(AODtrack->Eta()) > .8) || (AODtrack->GetTPCNcls() < 70)  || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.2)) continue;
            multTPC++;
        }
        for(Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++) { // fill global mult
            AliAODTrack *AODtrack=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(iTrack));
            if (!AODtrack) continue;
            if (!(AODtrack->TestFilterBit(16))) continue;
            if ((AODtrack->Pt() < .2) || (AODtrack->Pt() > 5.0) || (TMath::Abs(AODtrack->Eta()) > .8) || (AODtrack->GetTPCNcls() < 70) || (AODtrack->GetDetPid()->GetTPCsignal() < 10.0) || (AODtrack->Chi2perNDF() < 0.1)) continue;
            Double_t b[2] = {-99., -99.};
            Double_t bCov[3] = {-99., -99., -99.};
            AliAODTrack copy(*AODtrack);
            if (!(copy.PropagateToDCA(fVevent->GetPrimaryVertex(), fVevent->GetMagneticField(), 100., b, bCov))) continue;
            if ((TMath::Abs(b[0]) > 0.3) || (TMath::Abs(b[1]) > 0.3)) continue;
            multGlobal++;
            
        } //track loop
        
        HistTPCvsGlobalMultBeforeOutliers->Fill(multGlobal,multTPC);
        
        if(multTPC < (-36.73 + 1.48*multGlobal) || multTPC > (62.87 + 1.78*multGlobal)){pass = kFALSE;}
        
        if(!pass) return;
        HistTPCvsGlobalMultAfterOutliers->Fill(multGlobal,multTPC);
        
    }
    
    for(Int_t itrack = 0; itrack < ntracks; itrack++){
        
        AliAODTrack *track=dynamic_cast<AliAODTrack*>(fVevent->GetTrack(itrack));
        if(!track) continue;
        
        Float_t dcaXY = track->DCA();
        Float_t dcaZ  = track->ZAtDCA();
        
        TH2F* HistDCAbefore =(TH2F*)fListQAInfo->At(5);
        HistDCAbefore->Fill(dcaZ,dcaXY);
        
        Double_t p = -999, pT = -999, phi = -999, eta = -999, dEdx =-999;
        Double_t length = -999., beta =-999, tofTime = -999., tof = -999.;
        Double_t c = TMath::C()*1.E-9;// m/ns
        
        //cout<<"track->GetFilterMap()= "<<track->GetFilterMap()<<endl;
        if(!track->TestFilterBit(fFilterBit)) continue;
        
        p=track->P();
        pT=track->Pt();
        phi=track->Phi();
        eta=track->Eta();
        dEdx=track->GetDetPid()->GetTPCsignal();
        
        Float_t probMis = fPIDResponse->GetTOFMismatchProbability(track);
        if (probMis < 0.01) { //if u want to reduce mismatch using also TPC
            
            //if ( (track->IsOn(AliAODTrack::kTOFin)) && (track->IsOn(AliAODTrack::kTIME)) &&  (track->IsOn(AliAODTrack::kTOFout))) {
            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                
                tofTime = track->GetTOFsignal();//in ps
                length = track->GetIntegratedLength();
                
                tof = tofTime*1E-3; // ns
                //cout<<"tof = "<<tof<<endl;
                if (tof <= 0) continue;
                //cout<<"length = "<<length<<endl;
                if (length <= 0) continue;
                
                length = length*0.01; // in meters
                tof = tof*c;
                beta = length/tof;
                
                TH2F *HistBetavsPTOFbeforePID = (TH2F*)fListQAInfo->At(6);
                HistBetavsPTOFbeforePID ->Fill(track->P()*track->Charge(),beta);
            }//TOF signal
            
            TH2F *HistdEdxvsPTPCbeforePID = (TH2F*)fListQAInfo->At(7);
            HistdEdxvsPTPCbeforePID -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
        
            //QA plot
            TH1F *HistPhiDistBefore = (TH1F*)fListQAInfo->At(8);
            HistPhiDistBefore->Fill(phi);
            //
            TH1F *HistEtaDistBefore = (TH1F*)fListQAInfo->At(9);
            HistEtaDistBefore->Fill(eta);
        
        
            if(pT<0.1) continue;
            if(TMath::Abs(eta)>0.8) continue;
        
            Int_t TPCNcls = track->GetTPCNcls();
        
            if(TPCNcls<70 || dEdx<10) continue;
        
            // fill QA histograms
        
            TH2F* HistDCAAfter =(TH2F*)fListQAInfo->At(10);
            HistDCAAfter->Fill(dcaZ,dcaXY);
        
            TH1F *HistPhiDistAfter = (TH1F*)fListQAInfo->At(11);
            HistPhiDistAfter->Fill(phi);
        
            TH1F *HistEtaDistAfter = (TH1F*)fListQAInfo->At(12);
            HistEtaDistAfter->Fill(eta);
        
        
            Bool_t pWithinRange = kFALSE;
            Int_t p_bin = -999;
            Double_t pBins[50];
            for(int b=0;b<50;b++){pBins[b] = 0.1*b;}
            for(int i=0;i<50;i++){
                if(p>pBins[i] && p<(pBins[i]+0.1)){
                    pWithinRange = kTRUE;
                    p_bin = i;
                }
            }
        
            for (Int_t ispecie=0; ispecie<AliPID::kSPECIESC; ++ispecie){

                Double_t nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)ispecie);
                Double_t nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)ispecie);
            
                int i = ispecie - AliPID::kPion;
            
                if(fPIDcuts && pWithinRange){// for pions, kaons and protons only
                    if(ispecie==AliPID::kPion || ispecie==AliPID::kKaon || ispecie==AliPID::kProton){
                        int index = 50*i+p_bin;
                       
                        if(fCutContour[index]->IsInside(nSigmaTOF,nSigmaTPC)){//p_bin>7
                            TH3 *hist1 = (TH3*)fListQAtpctof->At(ispecie);
                            if (hist1){
                                hist1->Fill(nSigmaTPC,nSigmaTOF,p);}
                        }
                        if(p_bin>7 && fCutContour[index]->IsInside(nSigmaTOF,nSigmaTPC)){//p_bin>7
                            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                TH2F *HistBetavsPTOFafterPID = (TH2F*)fListQAInfo->At(13);
                                HistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
                            }
                            TH2F *HistdEdxvsPTPCafterPID = (TH2F*)fListQAInfo->At(14);
                            HistdEdxvsPTPCafterPID -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                        }
                        
                        if(p_bin<8 && nSigmaTPC<3 && nSigmaTPC>-3){//p_bin<8
                            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                TH2F *HistBetavsPTOFafterPID = (TH2F*)fListQAInfo->At(13);
                                HistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
                            }
                            TH2F *HistdEdxvsPTPCafterPID = (TH2F*)fListQAInfo->At(14);
                            HistdEdxvsPTPCafterPID -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                        }
                        //====================for low momentum PID based on purity by using only TPC information
                        Double_t LowPtPIDTPCnsigLow_Pion[6] = {-3,-3,-3,-3,-3,-3};
                        Double_t LowPtPIDTPCnsigLow_Kaon[6] = {-3,-2,0,-1.8,-1.2,-0.8}; //for 0.4<Pt<0.5 the purity is lower than 0.7
                        Double_t LowPtPIDTPCnsigHigh_Pion[6] ={2.4,3,3,3,2,1.4};
                        Double_t LowPtPIDTPCnsigHigh_Kaon[6] ={3,2.2,0,-0.2,1,1.8}; //for 0.4<Pt<0.5 the purity is lower than 0.7

                        
                        if(p_bin>7 && fCutContour[index]->IsInside(nSigmaTOF,nSigmaTPC)){//p_bin>7
                            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                TH2F *HistBetavsPTOFafterPID = (TH2F*)fListQAInfo->At(28);
                                HistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
                            }
                            TH2F *HistdEdxvsPTPCafterPID = (TH2F*)fListQAInfo->At(29);
                            HistdEdxvsPTPCafterPID -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                        }

                        if(p_bin<8){//p_bin<8
                            if((ispecie==AliPID::kPion && nSigmaTPC>LowPtPIDTPCnsigLow_Pion[p_bin-2] && nSigmaTPC<LowPtPIDTPCnsigHigh_Pion[p_bin-2]) || (ispecie==AliPID::kKaon && nSigmaTPC>LowPtPIDTPCnsigLow_Kaon[p_bin-2] && nSigmaTPC<LowPtPIDTPCnsigHigh_Kaon[p_bin-2]) || (ispecie==AliPID::kProton && nSigmaTPC>-3 && nSigmaTPC<3)){
                                if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                    TH2F *HistBetavsPTOFafterPID = (TH2F*)fListQAInfo->At(28);
                                    HistBetavsPTOFafterPID ->Fill(track->P()*track->Charge(),beta);
                                }
                                TH2F *HistdEdxvsPTPCafterPID = (TH2F*)fListQAInfo->At(29);
                                HistdEdxvsPTPCafterPID -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                            }
                        }

                        
                        TH2F *hTPCnSigmavsP = (TH2F*)fListQAtpctof->At(ispecie+AliPID::kSPECIESC);
                        if (hTPCnSigmavsP){
                            hTPCnSigmavsP->Fill(track->P()*track->Charge(),nSigmaTPC);}
                        
                        //=======================With TPC+TOF nsigma method Only!==============================
                        if(fCutContour[index]->IsInside(nSigmaTOF,nSigmaTPC)){
                            
                            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                TH2F *HistBetavsPTOFafterPIDTPCTOF = (TH2F*)fListQAInfo->At(15);
                                HistBetavsPTOFafterPIDTPCTOF ->Fill(track->P()*track->Charge(),beta);
                                if(ispecie==AliPID::kPion){
                                    TH2F *HistPion_BetavsPTOFafterPIDTPCTOF = (TH2F*)fListQAInfo->At(19);
                                    HistPion_BetavsPTOFafterPIDTPCTOF ->Fill(track->P()*track->Charge(),beta);
                                }
                                if(ispecie==AliPID::kKaon){
                                    TH2F *HistKaon_BetavsPTOFafterPIDTPCTOF = (TH2F*)fListQAInfo->At(21);
                                    HistKaon_BetavsPTOFafterPIDTPCTOF ->Fill(track->P()*track->Charge(),beta);
                                }
                                if(ispecie==AliPID::kProton){
                                    TH2F *HistProton_BetavsPTOFafterPIDTPCTOF = (TH2F*)fListQAInfo->At(23);
                                    HistProton_BetavsPTOFafterPIDTPCTOF ->Fill(track->P()*track->Charge(),beta);
                                }
                            }
                            TH2F *HistdEdxvsPTPCafterPIDTPCTOF = (TH2F*)fListQAInfo->At(16);
                            HistdEdxvsPTPCafterPIDTPCTOF -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                            if(ispecie==AliPID::kPion){
                                TH2F *HistPion_dEdxvsPTPCafterPIDTPCTOF = (TH2F*)fListQAInfo->At(20);
                                HistPion_dEdxvsPTPCafterPIDTPCTOF -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                                TH1F *HistPionEta = (TH1F*)fListQAInfo->At(25);
                                HistPionEta->Fill(eta);
                            }
                            if(ispecie==AliPID::kKaon){
                                TH2F *HistKaon_dEdxvsPTPCafterPIDTPCTOF = (TH2F*)fListQAInfo->At(22);
                                HistKaon_dEdxvsPTPCafterPIDTPCTOF -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                                TH1F *HistKaonEta = (TH1F*)fListQAInfo->At(26);
                                HistKaonEta->Fill(eta);
                            }
                            if(ispecie==AliPID::kProton){
                                TH2F *HistProton_dEdxvsPTPCafterPIDTPCTOF = (TH2F*)fListQAInfo->At(24);
                                HistProton_dEdxvsPTPCafterPIDTPCTOF -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                                TH1F *HistProtonEta = (TH1F*)fListQAInfo->At(27);
                                HistProtonEta->Fill(eta);

                            }
                            
                        }
                        //======================With TPC nsigma Only!
                        if(nSigmaTPC<3 && nSigmaTPC>-3){
                            if ( (track->IsOn(AliAODTrack::kITSin)) && (track->IsOn(AliAODTrack::kTOFpid)) ) {
                                TH2F *HistBetavsPTOFafterPIDTPConly = (TH2F*)fListQAInfo->At(17);
                                HistBetavsPTOFafterPIDTPConly ->Fill(track->P()*track->Charge(),beta);
                            }
                            TH2F *HistdEdxvsPTPCafterPIDTPConly = (TH2F*)fListQAInfo->At(18);
                            HistdEdxvsPTPCafterPIDTPConly -> Fill(track->P()*track->Charge(),dEdx); //TPC signal
                        }
                        //========================================================================================
                        
                        
                    }
                }
                if(!fPIDcuts){
                    TH3 *hist1 = (TH3*)fListQAtpctof->At(ispecie);
                    if (hist1){
                        hist1->Fill(nSigmaTPC,nSigmaTOF,p);}
                    
                    TH2F *hTPCnSigmavsP = (TH2F*)fListQAtpctof->At(ispecie+AliPID::kSPECIESC);
                    if (hTPCnSigmavsP){
                        hTPCnSigmavsP->Fill(track->P()*track->Charge(),nSigmaTPC);}
                    
                }
            }
        }//probMis
        
    }//track loop
    
    TH2F *HistTPCvsGlobalMultAfter = (TH2F*) fListQAInfo->At(30);
    HistTPCvsGlobalMultAfter->Fill(multGlobal,multTPC);
    
}
//_________________________________________
void AliAnalysisTaskPIDconfig::CheckCentrality(AliVEvent* event, Bool_t &centralitypass)
{
    // Check if event is within the set centrality range. Falls back to V0 centrality determination if no method is set
    if (!fUseCentrality) AliFatal("No centrality method set! FATAL ERROR!");
    Double_t centvalue = event->GetCentrality()->GetCentralityPercentile(fCentralityEstimator);
    //cout << "Centrality evaluated-------------------------: " << centvalue <<endl;
    if ((centvalue >= fCentralityPercentileMin) && (centvalue < fCentralityPercentileMax))
    {
        TH1F *hCentralityPass = (TH1F*)fListQAInfo->At(0);
        hCentralityPass->Fill(centvalue);
        //cout << "--------------Fill pass-------------------------"<<endl;
        centralitypass = kTRUE;
    }
    
}
//______________________________________________________________________________
void AliAnalysisTaskPIDconfig::GetPIDContours()
{
    fContoursFile = new TFile(Form("$ALICE_PHYSICS/PWGCF/FLOW/database/PIDCutContours_%i-%i.root",fCentralityPercentileMin,fCentralityPercentileMax));
    
    fCutContourList=(TDirectory*)fContoursFile->Get("Filterbit1");
    if(!fCutContourList){printf("The contour file is empty"); return;}

    Double_t pBinning[50];
    for(int b=0;b<50;b++){pBinning[b]=b;}
    TString species[3] = {"pion","kaon","proton"};
    
    for(int i=0;i<150;i++){
        int ispecie = i/50;
        int iPbin = i%50;
        TList *Species_contours = (TList*)fCutContourList->Get(species[ispecie]);
        //if(Species_contours){cout<<"Species_contours exists"<<endl;}
        
        TString Graph_Name = "contourlines_";
        Graph_Name += species[ispecie];
        Graph_Name += Form("%.f%.f-%i%icent",pBinning[iPbin],pBinning[iPbin]+1,fCentralityPercentileMin,fCentralityPercentileMax);
        //cout<<Graph_Name<<endl;
        fCutGraph[i] = (TGraph*)Species_contours->FindObject(Graph_Name);
        
        if(!fCutGraph[i]){cout<<"Contour Graph does not exist"<<endl; continue;}
        
        fCutContour[i] = new TCutG(Graph_Name.Data(),fCutGraph[i]->GetN(),fCutGraph[i]->GetX(),fCutGraph[i]->GetY());
        
    }
    
}
//______________________________________________________________________________
void AliAnalysisTaskPIDconfig::SetupTPCTOFqa()
{
    //
    // Create the qa objects for TPC + TOF combination
    
    
    //TPC and TOF signal vs. momentum
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIESC; ++ispecie){
        fhistNsigmaP = new TH3F(Form("NsigmaP_TPC_TOF_%s",AliPID::ParticleName(ispecie)),Form("TPC n#sigma vs. TOF n#sigma %s vs. p ;TPC n#sigma;TOF n#sigma;p [GeV]",AliPID::ParticleName(ispecie)),200,-20,20,200,-20,20,60,0.1,6);
        fListQAtpctof->Add(fhistNsigmaP);
    }
 
    //TPC signal vs. momentum
    for (Int_t ispecie=0; ispecie<AliPID::kSPECIESC; ++ispecie){
        fhistTPCnSigmavsP = new TH2F(Form("NsigmaP_TPC_%s",AliPID::ParticleName(ispecie)),Form("TPC n#sigma %s vs. p ;p [GeV];TPC n#sigma",AliPID::ParticleName(ispecie)),60,0,6,125,-5,20);
        fListQAtpctof->Add(fhistTPCnSigmavsP);
    }
    
}
//______________________________________________________________________________
void AliAnalysisTaskPIDconfig::SetupEventInfo()
{
    //event and track info
    
    fhistCentralityPass = new TH1F("fcentralityPass","centralityPass", 100,0,100);
    fListQAInfo->Add(fhistCentralityPass);
    
    fNoEvents = new TH1F("number of events","no. of events",1,0,1);
    fListQAInfo->Add(fNoEvents);
    
    fpVtxZ = new TH1F("pVtxZ","pVtxZ",100,-20,20);
    fListQAInfo->Add(fpVtxZ);
    
    fTPCvsGlobalMultBeforeOutliers = new TH2F("TPC vs. Global Multiplicity Before","TPC vs. Global Multiplicity Before",500,0,6000,500,0,6000);
    fListQAInfo->Add(fTPCvsGlobalMultBeforeOutliers);
    
    fTPCvsGlobalMultAfterOutliers = new TH2F("TPC vs. Global Multiplicity After outliers","TPC vs. Global Multiplicity After outliers",500,0,6000,500,0,6000);
    fListQAInfo->Add(fTPCvsGlobalMultAfterOutliers);
    
    fhistDCABefore = new TH2F("DCA xy vs z (before)","DCA before",200,0,10,200,0,10);
    fListQAInfo->Add(fhistDCABefore);
    
    fHistBetavsPTOFbeforePID = new TH2F("momentum vs beta before PID","momentum vs beta before PID",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistBetavsPTOFbeforePID);
    
    fHistdEdxvsPTPCbeforePID = new TH2F("momentum vs dEdx before PID","momentum vs dEdx before PID",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistdEdxvsPTPCbeforePID);
    
    fhistPhiDistBefore = new TH1F("Phi Distribution Before Cuts","Phi Distribution Before Cuts",200,0,6.4);
    fListQAInfo->Add(fhistPhiDistBefore);
    
    fhistEtaDistBefore = new TH1F("Eta Distribution Before Cuts","Eta Distribution Before Cuts",100,-2,2);
    fListQAInfo->Add(fhistEtaDistBefore);
    
    fhistDCAAfter = new TH2F("DCA xy vs z (after)","DCA after",200,0,10,200,0,10);
    fListQAInfo->Add(fhistDCAAfter);
    
    fhistPhiDistAfter = new TH1F("Phi Distribution After Cuts","Phi Distribution After Cuts",200,0,6.4);
    fListQAInfo->Add(fhistPhiDistAfter);
    
    fhistEtaDistAfter = new TH1F("Eta Distribution After Cuts","Eta Distribution After Cuts",200,-10,10);
    fListQAInfo->Add(fhistEtaDistAfter);
    
    fHistBetavsPTOFafterPID = new TH2F("momentum vs beta after PID","momentum vs beta after PID",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistBetavsPTOFafterPID);
    
    fHistdEdxvsPTPCafterPID = new TH2F("momentum vs dEdx after PID","momentum vs dEdx after PID",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistdEdxvsPTPCafterPID);
    
    fHistBetavsPTOFafterPIDTPCTOF = new TH2F("momentum vs beta after PID TPC+TOF","momentum vs beta after PID TPC+TOF",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistBetavsPTOFafterPIDTPCTOF);
    
    fHistdEdxvsPTPCafterPIDTPCTOF = new TH2F("momentum vs dEdx after PID TPC+TOF","momentum vs dEdx after PID TPC+TOF",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistdEdxvsPTPCafterPIDTPCTOF);
    
    fHistBetavsPTOFafterPIDTPConly = new TH2F("momentum vs beta after PID TPC only","momentum vs beta after PID TPC only",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistBetavsPTOFafterPIDTPConly);
    
    fHistdEdxvsPTPCafterPIDTPConly = new TH2F("momentum vs dEdx after PID TPC only","momentum vs dEdx after PID TPC only",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistdEdxvsPTPCafterPIDTPConly);
    
    fHistPion_BetavsPTOFafterPIDTPCTOF = new TH2F("Pion momentum vs beta after PID TPC+TOF","Pion momentum vs beta after PID TPC+TOF",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistPion_BetavsPTOFafterPIDTPCTOF);
    
    fHistPion_dEdxvsPTPCafterPIDTPCTOF = new TH2F("Pion momentum vs dEdx after PID TPC+TOF","Pion momentum vs dEdx after PID TPC+TOF",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistPion_dEdxvsPTPCafterPIDTPCTOF);
    
    fHistKaon_BetavsPTOFafterPIDTPCTOF = new TH2F("Kaon momentum vs beta after PID TPC+TOF","Kaon momentum vs beta after PID TPC+TOF",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistKaon_BetavsPTOFafterPIDTPCTOF);
    
    fHistKaon_dEdxvsPTPCafterPIDTPCTOF = new TH2F("Kaon momentum vs dEdx after PID TPC+TOF","Kaon momentum vs dEdx after PID TPC+TOF",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistKaon_dEdxvsPTPCafterPIDTPCTOF);
    
    fHistProton_BetavsPTOFafterPIDTPCTOF = new TH2F("Proton momentum vs beta after PID TPC+TOF","Proton momentum vs beta after PID TPC+TOF",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistProton_BetavsPTOFafterPIDTPCTOF);
    
    fHistProton_dEdxvsPTPCafterPIDTPCTOF = new TH2F("Proton momentum vs dEdx after PID TPC+TOF","Proton momentum vs dEdx after PID TPC+TOF",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistProton_dEdxvsPTPCafterPIDTPCTOF);
    
    fhistPionEtaDistAfter = new TH1F("Pion Eta Distribution After PID Cuts","Pion Eta Distribution After PID Cuts",100,-2,2);
    fListQAInfo->Add(fhistPionEtaDistAfter);
    
    fhistKaonEtaDistAfter = new TH1F("Kaon Eta Distribution After PID Cuts","Kaon Eta Distribution After PID Cuts",100,-2,2);
    fListQAInfo->Add(fhistKaonEtaDistAfter);

    fhistProtonEtaDistAfter = new TH1F("Proton Eta Distribution After PID Cuts","Proton Eta Distribution PID After Cuts",100,-2,2);
    fListQAInfo->Add(fhistProtonEtaDistAfter);
    
    fHistBetavsPTOFafterPID_2 = new TH2F("momentum vs beta after PID (PID in low Pt TPC only with Purity>0.7)","momentum vs beta after PID (PID in low Pt TPC only with Purity>0.7)",1000,-10.,10.,1000,0,1.2);
    fListQAInfo->Add(fHistBetavsPTOFafterPID_2);
    
    fHistdEdxvsPTPCafterPID_2 = new TH2F("momentum vs dEdx after PID (PID in low Pt TPC only with Purity>0.7)","momentum vs dEdx after PID (PID in low Pt TPC only with Purity>0.7)",1000,-10.,10.,1000,0,1000);
    fListQAInfo->Add(fHistdEdxvsPTPCafterPID_2);
    
    fTPCvsGlobalMultAfter = new TH2F("TPC vs. Global Multiplicity After","TPC vs. Global Multiplicity After",500,0,6000,500,0,6000);
    fListQAInfo->Add(fTPCvsGlobalMultAfter);
    
}


