/**************************************************************************
 * Copyright(c) 1998-2012, ALICE Experiment at CERN, All rights reserved. *
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

/*$Id: AliAnalysisTask for HF(DPlus:3Prongs)-Hadron/Kaon/K0 azimuthal correlation
 By: Jitendra Kumar: (jikumar@cern.ch)
 Last edited for MC Reco Corr + QA on 08.10.2016)
 */

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TString.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODPidHF.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliNormalizationCounter.h"
#include "AliVParticle.h"
#include "AliHFAssociatedTrackCuts.h"
#include "AliReducedParticle.h"
#include "AliHFCorrelator.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAODTracklets.h"
#include "AliVertexingHFUtils.h"
#include "AliAnalysisTaskSEDplusCorrelations.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSEDplusCorrelations)

//____________________| Default Constructor
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations():
AliAnalysisTaskSE(),
fSystem(kFALSE),
fReadMC(kFALSE),
fRecoTrk(kFALSE),
fMCParticle(kFALSE),
fMCGenEvType(kFALSE),
farrayMC(0x0),
fMixing(kFALSE),
fAssoParType(0),
fDplusCuts(0),
fAssoCuts(0),
fEffTrack(kFALSE),
fEffDplus(kFALSE),
fCentralityEstimator(0),
fEvalCentrality(kFALSE),
fMinCentrality(0),
fMaxCentrality(100),
fCentrOrMult(-99),
fTCconfig(kFALSE),
fUseBit(kTRUE),
fCorrelator(0x0),
fNPtBins(0),
fHistNEvents(0),
fHistNDplus(0),
fCounter(0x0),
fBinWidth(0),
fPoolByPool(kFALSE),
fWhichPool(0),
fCheckCutDist(kFALSE),
fAODProtection(1),
fCutSuffix(0x0),
fRawCutQA(0x0),
fOutput(0x0),
fOutputCorr(0X0)
{
    // Default constructor
}


//____________________| Specific Constructor
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const Char_t* name,AliRDHFCutsDplustoKpipi* DplusCuts, AliHFAssociatedTrackCuts *AsscCuts) :
AliAnalysisTaskSE(name),
fSystem(kFALSE),
fReadMC(kFALSE),
fRecoTrk(kFALSE),
fMCParticle(kFALSE),
fMCGenEvType(kFALSE),
farrayMC(0x0),
fMixing(kFALSE),
fAssoParType(0),
fDplusCuts(0),
fAssoCuts(AsscCuts),
fEffTrack(kFALSE),
fEffDplus(kFALSE),
fCentralityEstimator(0),
fEvalCentrality(kFALSE),
fMinCentrality(0),
fMaxCentrality(100),
fCentrOrMult(-99.0),
fTCconfig(kFALSE),
fUseBit(kTRUE),
fCorrelator(0x0),
fNPtBins(0),
fHistNEvents(0),
fHistNDplus(0),
fCounter(0x0),
fBinWidth(0.002),
fPoolByPool(kFALSE),
fWhichPool(0),
fCheckCutDist(kFALSE),
fAODProtection(1),
fCutSuffix(0x0),
fRawCutQA(0x0),
fOutput(0x0),
fOutputCorr(0X0)
{
    
    Info("AliAnalysisTaskSEDplusCorrelations","Calling Constructor");
    fNPtBins=DplusCuts->GetNPtBins();
    fDplusCuts=DplusCuts;
    
    DefineInput(0, TChain::Class());
    DefineOutput(1,TList::Class()); // Basic output slot (more needed)
    DefineOutput(2,TList::Class()); // Correlations form Data and MC
    DefineOutput(3,AliRDHFCutsDplustoKpipi::Class());  //D meson cuts
    DefineOutput(4,AliHFAssociatedTrackCuts::Class());  // Associated tracks cuts
    DefineOutput(5,AliNormalizationCounter::Class()); //Norm
    
}


//____________________| Soruce Operator
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const AliAnalysisTaskSEDplusCorrelations &source):
AliAnalysisTaskSE(source),
fSystem(source.fSystem),
fReadMC(source.fReadMC),
fRecoTrk(source.fRecoTrk),
fMCParticle(source.fMCParticle),
fMCGenEvType(source.fMCGenEvType),
farrayMC(source.farrayMC),
fMixing(source.fMixing),
fAssoParType(source.fAssoParType),
fDplusCuts(source.fDplusCuts),
fAssoCuts(source.fAssoCuts),
fEffTrack(source.fEffTrack),
fEffDplus(source.fEffDplus),
fCentralityEstimator(source.fCentralityEstimator),
fEvalCentrality(source.fEvalCentrality),
fMinCentrality(source.fMinCentrality),
fMaxCentrality(source.fMaxCentrality),
fCentrOrMult(source.fCentrOrMult),
fTCconfig(source.fTCconfig),
fUseBit(source.fUseBit),
fCorrelator(source.fCorrelator),
fNPtBins(source.fNPtBins),
fHistNEvents(source.fHistNEvents),
fHistNDplus(source.fHistNDplus),
fCounter(source.fCounter),
fBinWidth(source.fBinWidth),
fPoolByPool(source.fPoolByPool),
fWhichPool(source.fWhichPool),
fCheckCutDist(source.fCheckCutDist),
fAODProtection(source.fAODProtection),
fCutSuffix(source.fCutSuffix),
fRawCutQA(source.fRawCutQA),
fOutput(source.fOutput),
fOutputCorr(source.fOutputCorr)
{
    
}


//____________________| Destructor
AliAnalysisTaskSEDplusCorrelations::~AliAnalysisTaskSEDplusCorrelations()
{
    Info("AliAnalysisTaskSEDplusCorrelations","Calling Destructor");
    
    if(fOutputCorr) {delete fOutputCorr; fOutputCorr = 0;}
    if(fOutput) {delete fOutput; fOutput = 0;}
    if(farrayMC) {delete farrayMC; farrayMC = 0;}
    if(fHistNEvents) {delete fHistNEvents; fHistNEvents = 0;}
    if(fHistNDplus) {delete fHistNDplus; fHistNDplus = 0;}
    if(fCounter) {delete fCounter; fCounter = 0;}
    if(fDplusCuts) {delete fDplusCuts; fDplusCuts = 0;}
    if(fAssoCuts) {delete fAssoCuts; fAssoCuts=0;}
    if(fCorrelator) {delete fCorrelator; fCorrelator = 0;}
    
}

//____________________| Assignment Operator
AliAnalysisTaskSEDplusCorrelations& AliAnalysisTaskSEDplusCorrelations::operator=(const AliAnalysisTaskSEDplusCorrelations& orig)
{
    if (&orig == this) return *this; //if address is the same (same object), returns itself
    AliAnalysisTaskSE::operator=(orig); //Uses the AliAnalysisTaskSE operator to assign the inherited part of the class
    fSystem = orig.fSystem;
    fReadMC = orig.fReadMC;
    fRecoTrk = orig.fRecoTrk;
    fMCParticle = orig.fMCParticle;
    fMCGenEvType = orig.fMCGenEvType;
    farrayMC = orig.farrayMC;
    fMixing = orig.fMixing;
    fAssoParType = orig.fAssoParType;
    fDplusCuts = orig.fDplusCuts;
    fAssoCuts = orig.fAssoCuts;
    fEffTrack = orig.fEffTrack;
    fEffDplus=orig.fEffDplus;
    fCentralityEstimator = orig.fCentralityEstimator;
    fEvalCentrality = orig.fEvalCentrality;
    fMinCentrality = orig.fMinCentrality;
    fMaxCentrality = orig.fMaxCentrality;
    fCentrOrMult = orig.fCentrOrMult;
    fTCconfig=orig.fTCconfig;
    fUseBit=orig.fUseBit;
    fCorrelator=orig.fCorrelator;
    fNPtBins=orig.fNPtBins;
    fHistNEvents=orig.fHistNEvents;
    fHistNDplus=orig.fHistNDplus;
    fCounter=orig.fCounter;
    fBinWidth=orig.fBinWidth;
    fPoolByPool=orig.fPoolByPool;
    fWhichPool=orig.fWhichPool;
    fCheckCutDist=orig.fCheckCutDist;
    fAODProtection=orig.fAODProtection;
    fCutSuffix=orig.fCutSuffix;
    fRawCutQA=orig.fRawCutQA;
    fOutput=orig.fOutput;
    fOutputCorr=orig.fOutputCorr;
    return *this; //returns pointer of the class
    
}


//____________________| Inilizations
void AliAnalysisTaskSEDplusCorrelations::Init()
{
    
    Info("AliAnalysisTaskSEDplusCorrelations","Calling Inilizations");
    
    //Copy of cuts objects
    AliRDHFCutsDplustoKpipi* DplusCutsCopy = new AliRDHFCutsDplustoKpipi(*fDplusCuts);
    if(!DplusCutsCopy)return;
    const char* name=GetOutputSlot(3)->GetContainer()->GetName();
    DplusCutsCopy->SetName(name);
    
    if(fEvalCentrality){
        Bool_t isCentEstimatorOk = kFALSE;
        if(fCentralityEstimator==fDplusCuts->GetUseCentrality())isCentEstimatorOk = kTRUE;
        if(!isCentEstimatorOk) {
            AliFatal(Form("Chosen centrality estimator is conflict with D cuts file"));
            return;
        }
    }
    
    //Posting cuts for D+ and Associated tracks
    PostData(3,DplusCutsCopy);
    PostData(4,fAssoCuts);
    return;
}


//____________________| UserCreateOutputObjects
void AliAnalysisTaskSEDplusCorrelations::UserCreateOutputObjects()
{
    
    Info("AliAnalysisTaskSEDplusCorrelations","Creating UserCreateOutputObjects");
    
    // Category one: Basic O/P
    fOutput = new TList();
    fOutput->SetOwner();
    fOutput->SetName("BasicQAHisto");
    
    fOutputCorr = new TList();
    fOutputCorr->SetOwner();
    fOutputCorr->SetName("DplushCorrHisto");
    
    fHistNEvents = new TH1F("fHistNEvents", "Events Stats", 10, -0.5 , 9.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"All");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"Pileup Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"Cent-Mult Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"Bad-Vtx Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(5,"Trigger Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(6,"Zvtx Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(7,"PS Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(8,"HFCorr Interface Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(9,"Accepted");
    fHistNEvents->GetXaxis()->SetBinLabel(10,"mismatch AOD/dAOD");
    fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistNEvents->SetFillColor(kBlue);
    fHistNEvents->SetFillStyle(3015);
    fHistNEvents->SetMinimum(0);
    //fHistNEvents->Sumw2();
    fOutput->Add(fHistNEvents);
    
    fHistNDplus = new TH1F("fHistNDplus", "Dplus Stats ", 9, -0.5 , 8.5);
    fHistNDplus->GetXaxis()->SetBinLabel(1,"All");
    fHistNDplus->GetXaxis()->SetBinLabel(2,"-pT Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(3,"<2GeV Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(4,"SelBit Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(5,"Topocut Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(6,"TC Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(7,"FiduAcc Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(8,"Accepted");
    fHistNDplus->GetXaxis()->SetBinLabel(9,"D+ failed");
    fHistNDplus->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistNDplus->SetMinimum(0);
    fHistNDplus->SetFillColor(kRed);
    fHistNDplus->SetFillStyle(3015);
    //fHistNDplus->Sumw2();
    fOutput->Add(fHistNDplus);
    
    if(fRawCutQA)fCutSuffix = "BeforeDSel";
    else if(!fRawCutQA)fCutSuffix = "AfterDSel";
    
    HistoNomenclature();// Hist Nomenclature function
    TString normName="NormalizationCounter";
    AliAnalysisDataContainer *cont = GetOutputSlot(5)->GetContainer();
    if(cont)normName=(TString)cont->GetName();
    fCounter = new AliNormalizationCounter(normName.Data());
    fCounter->Init();
    
    //Setting up from Correlator
    Double_t Pi = TMath::Pi();
    fCorrelator = new AliHFCorrelator("Correlator",fAssoCuts,fSystem,fDplusCuts);
    fCorrelator->SetDeltaPhiInterval(-0.5*Pi, 1.5*Pi);
    fCorrelator->SetEventMixing(fMixing);
    fCorrelator->SetAssociatedParticleType(fAssoParType);
    //fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
    fCorrelator->SetUseReco(fRecoTrk);
    fCorrelator->SetUseMC(fReadMC);
    
    Bool_t pooldef = fCorrelator->DefineEventPool();
    if(!pooldef) AliInfo("Warning:: Event pool not defined properly");
    
    //Post Data
    PostData(1,fOutput);
    PostData(2,fOutputCorr);
    PostData(5,fCounter);
    
}



//____________________| UserExec
void AliAnalysisTaskSEDplusCorrelations::UserExec(Option_t *) {
    
    
    Info("AliAnalysisTaskSEDplusCorrelations","Running... UserExec");
    AliAODEvent *aod = dynamic_cast<AliAODEvent*>(InputEvent());
    TClonesArray *array3Prong = 0;
    
    if(fAODProtection>=0){ //New added 05.10
        Int_t matchingAODdeltaAODlevel = AliRDHFCuts::CheckMatchingAODdeltaAODevents();
        if (matchingAODdeltaAODlevel<0 || (matchingAODdeltaAODlevel==0 && fAODProtection==1)){
            fHistNEvents->Fill(9);
            return;
        }
    }
    
    if(!aod && AODEvent() && IsStandardAOD()){
        aod = dynamic_cast<AliAODEvent*> (AODEvent());
        AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        
        if(aodHandler->GetExtensions()) {
            AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();
            array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
        }
    }else if(aod)array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    
    if(!aod || !array3Prong){
        printf("AliAnalysisTaskSEDplusCorrelationselation::UserExec: AOD  Charm3Prong branch not found!\n");
        return;
    }
    
    if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
    
    fCounter->StoreEvent(aod,fDplusCuts,fReadMC);
    fHistNEvents->Fill(0);
    
    //MC Kinemactics
    if(fReadMC){
        
        farrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
        if(!farrayMC){
            printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC particles branch not found!\n");
            return;
        }
        
        AliAODMCHeader  *mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
        if(!mcHeader) {
            printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC header branch not found!\n");
            return;
        }
        
        Double_t zVtxMC = mcHeader->GetVtxZ();
        if(TMath::Abs(zVtxMC)>10) return;
        if(aod->GetTriggerMask()==0 && (aod->GetRunNumber()>=195344 && aod->GetRunNumber()<=195677)) return;
        
        if(fMCGenEvType){
            Bool_t isMCeventgood = kFALSE;
            Int_t      eventType = mcHeader->GetEventType();
            Int_t      NMCevents = fAssoCuts->GetNofMCEventType();
            for(Int_t k=0; k<NMCevents; k++){
                Int_t * MCEventType = fAssoCuts->GetMCEventType();
                if(eventType == MCEventType[k]) isMCeventgood= kTRUE;
            }
            if(NMCevents && !isMCeventgood){
                if(fDebug>2) std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
                return;
            }
        }
        fCorrelator->SetMCArray(farrayMC);
    }
    
    //Event Selection
    if(fRecoTrk && !fDplusCuts->IsEventSelected(aod)){
        //Remarks 1. ON for Data+MC=Reco
        if(fDplusCuts->GetWhyRejection()==1)fHistNEvents->Fill(1);//pileup
        if(fDplusCuts->GetWhyRejection()==5)fHistNEvents->Fill(4);//triggers
        if(fDplusCuts->GetWhyRejection()==6)fHistNEvents->Fill(5);//Zvtx
        if(fDplusCuts->GetWhyRejection()==7)fHistNEvents->Fill(6);//PS
        return;
    }else if(fReadMC && !fRecoTrk){
        //Kinematics-Montecarlo
        AliAODVertex *vertex = (AliAODVertex*)aod->GetPrimaryVertex();
        if(!vertex)return;
        Double_t zVtxMCreco = vertex->GetZ();
        if(TMath::Abs(zVtxMCreco)>10) return;
    }
    
    if(fEvalCentrality){
        fMinCentrality = fDplusCuts->GetMinCentrality();
        fMaxCentrality = fDplusCuts->GetMaxCentrality();
        fCentrOrMult = fDplusCuts->GetCentrality(aod);
        if(fCentrOrMult<fMinCentrality || fCentrOrMult >fMaxCentrality)fHistNEvents->Fill(2);
    }
    else if(!fEvalCentrality){
        Double_t count = -1, mineta = -1.0, maxeta = 1.0;
        AliAODTracklets* tracklets = aod->GetTracklets();
        if(!tracklets)return;
        Int_t nTr=tracklets->GetNumberOfTracklets();
        for(Int_t iTr=0; iTr<nTr; iTr++){
            Double_t theta=tracklets->GetTheta(iTr);
            Double_t eta=-TMath::Log(TMath::Tan(theta/2.));
            if(eta>mineta && eta<maxeta) count++;
        }
        fCentrOrMult = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
        if(!fMixing)((TH2F*)fOutput->FindObject("h2SPDTrkVsTrkMult"))->Fill(aod->GetNumberOfTracks(),fCentrOrMult);
    }
    
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    Bool_t isGoodVtx=kFALSE;
    TString primTitle = vtx1->GetTitle();
    if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)isGoodVtx=kTRUE;
    if(!isGoodVtx){
        fHistNEvents->Fill(3);
        return;
    }
    
    //Pool by Pool Setting
    if(fPoolByPool)fWhichPool = fAssoCuts->GetPoolBin(fCentrOrMult, vtx1->GetZ());
    else if(!fPoolByPool)fWhichPool=0; //integrated pools
    
    //HFCorrelator interface
    fCorrelator->SetPidAssociated();
    fCorrelator->SetAODEvent(aod);
    Bool_t correlatorON = fCorrelator->Initialize(); //
    if(!correlatorON) {
        AliInfo("AliHFCorrelator didn't initialize the pool correctly or processed a bad event");
        fHistNEvents->Fill(7);
        return;
    }
    fHistNEvents->Fill(8);
    
    //D+ Particle loop and Correlation
    Int_t nDplusKpipi = array3Prong->GetEntriesFast();
    Int_t       labDp = -1;
    Int_t pdgDgDplustoKpipi[3]={321,211,211};
    //printf("Number of D+->Kpipi: %d and of tracks: %d\n",nDplusKpipi,aod->GetNumberOfTracks());
    
    AliAnalysisVertexingHF *vHF = new AliAnalysisVertexingHF();
    
    TString isDataOrMC = "";
    if(!fReadMC)isDataOrMC = "Data";
    else if(fReadMC)isDataOrMC = "MC";
    
    
    for (Int_t iDplusKpipi = 0; iDplusKpipi < nDplusKpipi; iDplusKpipi++){
        
        AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(iDplusKpipi); // D+
        if(!d)continue;
        fHistNDplus->Fill(0);
        
        if(!(vHF->FillRecoCand(aod,d))) {//Fill the data members of the candidate only if they are empty.
            fHistNDplus->Fill(8); //monitor how often this fails
            continue;
        }
        
        if(d->Pt()<0.)fHistNDplus->Fill(1);
        
        if(d->Pt()<2.){
            fHistNDplus->Fill(2);
            continue;
        }
        
        if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
            fHistNDplus->Fill(3);
            continue;
        }
        
        if(!fMixing && fRecoTrk && fCheckCutDist && fRawCutQA)DoDplusCutDistFill(d);
        
        //Tight cuts
        Int_t passTightCuts = fDplusCuts->IsSelected(d,AliRDHFCuts::kAll,aod);
        
        //Topological cuts
        Int_t passTopologicalCuts=fDplusCuts->GetIsSelectedCuts();
        if(!passTopologicalCuts) {
            fHistNDplus->Fill(4);
            continue;
        }
        
        //Int_t passPIDCuts=fDplusCuts->GetIsSelectedPID();
        //if(!passPIDCuts) continue; //not needed
        if(fTCconfig){
            if(!passTightCuts){
                fHistNDplus->Fill(5);
                continue;
            }
        }
        
        Bool_t unsetvtx=kFALSE;
        if(!d->GetOwnPrimaryVtx()){
            d->SetOwnPrimaryVtx(vtx1);
            unsetvtx=kTRUE;
        }
        
        Bool_t recVtx=kFALSE;
        AliAODVertex *origownvtx=0x0;
        if(fDplusCuts->GetIsPrimaryWithoutDaughters()){
            if(d->GetOwnPrimaryVtx()) origownvtx=new AliAODVertex(*d->GetOwnPrimaryVtx());
            if(fDplusCuts->RecalcOwnPrimaryVtx(d,aod))recVtx=kTRUE;
            else fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
        }
        
        
        Double_t ptDplusCand = d->Pt();
        Double_t rapDplusCand=d->YDplus();
        Int_t iPtDplus = fDplusCuts->PtBin(ptDplusCand);
        if(iPtDplus<0){
            //Vertexing cleaning...
            if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
            if(unsetvtx)d->UnsetOwnPrimaryVtx();
            continue;
        }
        
        // D+ Fiducial Acceptance in pt and eta
        if(!fMixing && fRecoTrk && fCheckCutDist && !fRawCutQA){
            ((TH2F*)fOutput->FindObject(Form("hDistDpluswoYCutVsPt_Bin%d",iPtDplus)))->Fill(d->Pt(),d->YDplus());
        }
        
        if(!fDplusCuts->IsInFiducialAcceptance(ptDplusCand,rapDplusCand)){
            //Vertexing cleaning...
            if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
            if(unsetvtx) d->UnsetOwnPrimaryVtx();
            fHistNDplus->Fill(6);
            continue;
        }
        
        Int_t fDPlusorig = 0;
        if(fReadMC){
            //Match to MC
            labDp = d->MatchToMC(411,farrayMC,3,pdgDgDplustoKpipi);
            if(labDp>=0)fDPlusorig = CheckOriginPartOfDPlus(farrayMC,(AliAODMCParticle*)farrayMC->At(labDp));
            //cout << "Origin of Particle --> " << fDPlusorig << endl;
        }
        fHistNDplus->Fill(7);
        
        
        //Computing Correlations: Below Funtions
        HadronCorrelations(d,fDPlusorig);
        
        //Vertexing cleaning...
        if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        
    }// loop over D+
    
    
    //EVT MIXING PART
    if(fMixing ) {
        Bool_t updated = fCorrelator->PoolUpdate();
        if(!updated) AliInfo("Pool was not updated");
    }
    
    delete vHF;
    // Posting Slots
    PostData(1,fOutput);
    PostData(2,fOutputCorr);
    PostData(5,fCounter);
    return;
    
}



//____________________| Calculating Correlation
void AliAnalysisTaskSEDplusCorrelations::HadronCorrelations(AliAODRecoDecayHF3Prong* d, Int_t fMCDplusMom){
    
    
    //Switch for String for Hist name
    TString parttype = "", datatype = "";
    if(fAssoParType==1)parttype += "Hdron";
    else if(fAssoParType==2)parttype += "Kaons";
    else if(fAssoParType==3)parttype += "kZero";
    
    if(fRecoTrk){
        if(!fReadMC)datatype += "Data";
        else if(fReadMC)datatype += "MCrc";
    }
    
    TString isDataOrMC = "";
    if(!fReadMC)isDataOrMC = "Data";
    else if(fReadMC)isDataOrMC = "MC";
    
    // D+ Quantities
    Double_t mDplus       = d->InvMassDplus();
    Double_t ptDplusCand  = d->Pt();
    Double_t etaDplusCand = d->YDplus();
    Double_t phiDplusCand = d->Phi();
    
    Double_t effDplus =1.0;
    if(fEffDplus)effDplus = fAssoCuts->GetTrigWeight(ptDplusCand,fCentrOrMult); //Dplus efficiency
    if(effDplus<1.0e-9) effDplus=1; // case of 0 bin content
    
    Int_t iPtBin = fDplusCuts->PtBin(ptDplusCand); // Pt bins
    if(iPtBin<0) return ;
    
    if(!fMixing && fRecoTrk && fCheckCutDist && !fRawCutQA)DoDplusCutDistFill(d);
    
    //Correlation with SE or ME
    Double_t nSparceCorrDplusinfo[1] = {mDplus};
    if(!fMixing){
        ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseMa_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrDplusinfo,1.0/effDplus);
        if(fEffDplus && !fReadMC){
            ((TH1F*)fOutput->FindObject(Form("histgrm1Ma_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(mDplus); //Dplus mass w/o efficiency
        }
    }
    
    Double_t phiDplus = fCorrelator->SetCorrectPhiRange(phiDplusCand);
    fCorrelator->SetTriggerParticleProperties(ptDplusCand,phiDplus,etaDplusCand);
    
    //Pool Setting and Pool events
    Int_t NofEventsinPool = 1; // SE
    Bool_t execPool = fCorrelator->ProcessEventPool();
    if(fMixing){
        if(!execPool){
            AliInfo("Mixed event analysis: pool is not ready");
            NofEventsinPool = 0;
        }else if(execPool){
            NofEventsinPool = fCorrelator->GetNofEventsInPool();
        }
    }
    
    
    //D+ Daughters ids
    Int_t nIDs[3] = {-9999999};
    nIDs[0] = ((AliAODTrack*)d->GetDaughter(0))->GetID();
    nIDs[1] = ((AliAODTrack*)d->GetDaughter(1))->GetID();
    nIDs[2] = ((AliAODTrack*)d->GetDaughter(2))->GetID();
    
    //Correlation with associated tracks
    for (Int_t jMix =0; jMix < NofEventsinPool; jMix++){
        
        Bool_t analyzetracks = fCorrelator->ProcessAssociatedTracks(jMix);
        if(!analyzetracks) {
            AliInfo("AliHFCorrelator::Cannot process the track array");
            continue;
        }
        
        // Leading Particle varibales
        Double_t refpt = 0, effLead = 1, TotaleffLead =1;
        Double_t ptleadHadron=0;
        Double_t DeltaphiLead=0, DeltaetaLead=0;
        
        //Correlation with associated tracks
        Int_t DCount = 0;
        for (Int_t iTrack = 0; iTrack<fCorrelator->GetNofTracks(); iTrack++) {
            
            // Correlation using HFCorrelator
            Bool_t runcorrelation = fCorrelator->Correlate(iTrack);
            if(!runcorrelation) continue;
            
            AliReducedParticle* redpart = fCorrelator->GetAssociatedParticle();
            if(!redpart)continue;
            
            Double_t phiHad=redpart->Phi();
            Double_t ptHad=redpart->Pt();
            Int_t label = redpart->GetLabel();
            Int_t trackid = redpart->GetID();
            
            DCount++;
            if(!fMixing && DCount==1){
                //SE only + on request
                ((TH1F*)fOutput->FindObject(Form("AssoTrkPhi_%s",isDataOrMC.Data())))->Fill(redpart->Phi());
                ((TH1F*)fOutput->FindObject(Form("AssoTrkEta_%s",isDataOrMC.Data())))->Fill(redpart->Eta());
                ((TH1F*)fOutput->FindObject(Form("AssoTrkPt_%s",isDataOrMC.Data())))->Fill(redpart->Pt());
            }
            
            phiHad = fCorrelator->SetCorrectPhiRange(phiHad);
            
            Double_t effTr =1;
            if(fEffTrack)effTr = redpart->GetWeight(); //track efficiency
            if(effTr<1.e-9) effTr=1;
            Double_t effTotal = effTr*effDplus;
            
            if(trackid < 0) continue;
            if (!fMixing)if( trackid == nIDs[0] || trackid == nIDs[1] || trackid == nIDs[2]) continue;
            
            
            // leading particle correlation
            if (ptHad > refpt) {
                refpt = ptHad;
                DeltaphiLead  = fCorrelator->GetDeltaPhi();
                DeltaetaLead  = fCorrelator->GetDeltaEta();
                ptleadHadron  = ptHad;
                effLead = redpart->GetWeight();
                TotaleffLead = effLead*effDplus;
            }
            
            CorrelationNSparsePlots(d, redpart, iPtBin, fMCDplusMom, 1.0/effTotal);
            
        }//track loop end here
        
        //Leading particle correlations
        //pT bins setting
        Double_t ptLim_Sparse = 0.0;
        if(!fMixing){
            if(!fReadMC)ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLa_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
            if(fReadMC)ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLc_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        }else if(fMixing){
            ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLa_%s_%s_Bin%d_evMix",parttype.Data(), datatype.Data(), iPtBin)))->GetAxis(3)->GetXmax();
        }
        
        if(ptleadHadron > ptLim_Sparse) ptleadHadron = ptLim_Sparse-0.01; //filling all above pT in last bin
        Double_t nSparceCorrLeadPart[5] = {mDplus, DeltaphiLead, DeltaetaLead, ptleadHadron, fWhichPool+0.5};
        if(!fMixing){
            if(!fReadMC)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLa_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
            if(fReadMC){
                if(fMCDplusMom==4)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLc_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
                else if(fMCDplusMom==5)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLb_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
                else ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLX_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
            }
        }else if(fMixing){
            ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseLa_%s_%s_Bin%d_evMix",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
        }
    }  // event loop.
    
}


//____________________|Correlation with ThnSparse Histograms
void AliAnalysisTaskSEDplusCorrelations::CorrelationNSparsePlots(AliAODRecoDecayHF3Prong *d, AliReducedParticle* track, Int_t iPtBin, Int_t origDplus, Double_t weightEff) {
    
    iPtBin = fDplusCuts->PtBin(d->Pt());
    Double_t mDplus = d->InvMassDplus();
    Double_t deltaPhi = fCorrelator->GetDeltaPhi();
    Double_t deltaEta = fCorrelator->GetDeltaEta();
    Double_t ptTrack = track->Pt();
    
    TString parttype = "", datatype = "";
    if(fAssoParType==1)parttype += "Hdron";
    else if(fAssoParType==2)parttype += "Kaons";
    else if(fAssoParType==3)parttype += "kZero";
    
    if(fRecoTrk){
        if(!fReadMC)datatype += "Data";
        else if(fReadMC)datatype += "MCrc";
    }
    
    Double_t ptLim_Sparse = 0.0;
    if(!fMixing){
        if(!fReadMC)ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCa_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        if(fReadMC)ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCc_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
    }else if(fMixing){
        ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCa_%s_%s_Bin%d_evMix",parttype.Data(), datatype.Data(), iPtBin)))->GetAxis(3)->GetXmax();
    }
    
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01; //filling all above pT in last bin
    Double_t nSparceCorr[5] = {mDplus, deltaPhi, deltaEta, ptTrack, fWhichPool+0.5};
    if(!fMixing){
        if(!fReadMC)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCa_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
        if(fReadMC){
            if(origDplus==4)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCc_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
            else if(origDplus==5)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCb_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);//FIX
            else ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCX_%s_%s_Bin%d",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
        }
    }else if(fMixing){
        ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseCa_%s_%s_Bin%d_evMix",parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
    }
    
    return;
}



//____________________|  Cut Distributions
void AliAnalysisTaskSEDplusCorrelations::DoDplusCutDistFill(AliAODRecoDecayHF3Prong *d) {
    
    TString isDataOrMC = "";
    if(!fReadMC)isDataOrMC = "Data";
    else if(fReadMC)isDataOrMC = "MC";
    
    Int_t iPtBin = fDplusCuts->PtBin(d->Pt());
    
    Double_t MaxDCA = -9999.;
    for(Int_t i=0; i<3; i++)if(d->GetDCA(i)>MaxDCA)MaxDCA=d->GetDCA(i);
    Double_t d0Square = d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
    
    //DPlus Cuts dist
    ((TH2F*)fOutput->FindObject(Form("hDistPtvsd0_DauPion1_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->PtProng(0), d->Getd0Prong(0));
    ((TH2F*)fOutput->FindObject(Form("hDistPtvsd0_DauKaon_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->PtProng(1), d->Getd0Prong(1));
    ((TH2F*)fOutput->FindObject(Form("hDistPtvsd0_DauPion2_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->PtProng(2), d->Getd0Prong(2));
    ((TH2F*)fOutput->FindObject(Form("hDistCosPAvsDCA_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->CosPointingAngle(), MaxDCA);
    ((TH1F*)fOutput->FindObject(Form("hDistDplusCosPA_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->CosPointingAngle());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusCosPAxy_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->CosPointingAngleXY());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusDCA_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(MaxDCA);
    ((TH1F*)fOutput->FindObject(Form("hDistDplusd0Square_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d0Square);
    ((TH1F*)fOutput->FindObject(Form("hDistDplusSigmaVert_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->GetSigmaVert());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusDecayLen_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->DecayLength());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusDecayLenxy_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->DecayLengthXY());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusNormDecayLenxy_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->NormalizedDecayLengthXY());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusDist12_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->GetDist12toPrim());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusDist23_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->GetDist23toPrim());
    
    //DPlus QA
    ((TH1F*)fOutput->FindObject(Form("hDistDplusRap_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->YDplus());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusEta_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->Eta());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusPhi_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->Phi());
    ((TH1F*)fOutput->FindObject(Form("hDistDplusMass_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->InvMassDplus());
    ((TH2F*)fOutput->FindObject(Form("hDistDplusMassVsPhi_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->InvMassDplus(), d->Phi());
    ((TH2F*)fOutput->FindObject(Form("hDistDplusMassVsEta_Bin%d_%s_%s",iPtBin,isDataOrMC.Data(),fCutSuffix.Data())))->Fill(d->InvMassDplus(), d->Eta());
    
    
}

//____________________|  Histograms Nomenclature
void AliAnalysisTaskSEDplusCorrelations::HistoNomenclature() {
    
    Double_t Pi = TMath::Pi();
    
    //D+ Candidate Vars: InvMass
    Float_t         range  = 0.200;
    Float_t fLowmasslimit  = 1.865 - range;
    Float_t  fUpmasslimit  = 1.865 + range;
    Float_t         width  = fBinWidth;
    Int_t           nbins  = (Int_t)((fUpmasslimit-fLowmasslimit)/width+0.5);
    Int_t     missingbins  = 4 - nbins%4;
    nbins  =  nbins + missingbins;
    width  = (fUpmasslimit-fLowmasslimit)/nbins; // new width
    fBinWidth=width;
    
    Int_t nPoolCorr=1;
    if(fPoolByPool)nPoolCorr= fAssoCuts->GetNZvtxPoolBins()*fAssoCuts->GetNCentPoolBins();
    
    
    //1. Invariant Mass
    Int_t     nBinsDinfo[1] = {nbins};
    Double_t binMinDinfo[1] = {fLowmasslimit};
    Double_t binMaxDinfo[1] = {fUpmasslimit};
    
    //2. SEorME Correlations Particle Vars: InvMass, DeltaPhi, DeltaEta, Centrality
    Int_t           nVarsBins[5] = {nbins,              32,      16,    6,     nPoolCorr};
    Double_t nMinimumEdgeBins[5] = {fLowmasslimit,   -Pi/2,    -1.6,    0.00,  0.};  //is the minimum for all the bins
    Double_t nMaximumEdgeBins[5] = {fUpmasslimit,   3*Pi/2,     1.6,    3.00,  static_cast<double>(nPoolCorr)};  //is the maximum for all the bins
    
    
    //Data
    TString namePlotThnRecoWOeffMass = "histgrm1" ;
    TString namePlotThnRecoMass      = "hnSparse" ;
    TString namePlotThnRecoCorr      = "hnSparse" ;
    TString namePlotThnRecoCorrL     = "hnSparse" ;
    
    //MC reco + origin
    TString namePlotThnRecoMassOrgnAll= "hnSparse" ;
    TString namePlotThnRecoCorrOrgnc  = "hnSparse" ;
    TString namePlotThnRecoCorrOrgnb  = "hnSparse" ;
    TString namePlotThnRecoCorrOrgnX  = "hnSparse" ;
    TString namePlotThnRecoCorrLOrgnc = "hnSparse" ;
    TString namePlotThnRecoCorrLOrgnb = "hnSparse" ;
    TString namePlotThnRecoCorrLOrgnX = "hnSparse" ;
    
    TString nameMasterStting   = "";
    if(fAssoParType==1)nameMasterStting += "_Hdron";
    else  if(fAssoParType==2)nameMasterStting += "_Kaons";
    else  if(fAssoParType==3)nameMasterStting += "_kZero";
    else  nameMasterStting = "_Nulls";
    
    if(!fReadMC){
        nameMasterStting +="_Data";
        namePlotThnRecoWOeffMass  += Form("Ma%s_Bin", nameMasterStting.Data());
        namePlotThnRecoMass  += Form("Ma%s_Bin", nameMasterStting.Data());
        namePlotThnRecoCorr  += Form("Ca%s_Bin", nameMasterStting.Data());
        namePlotThnRecoCorrL += Form("La%s_Bin", nameMasterStting.Data());
    }else if(fReadMC){
        if(fRecoTrk){
            nameMasterStting +="_MCrc";
            namePlotThnRecoMassOrgnAll+=  Form("Ma%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrOrgnc  +=  Form("Cc%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrOrgnb  +=  Form("Cb%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrOrgnX  +=  Form("CX%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrLOrgnc +=  Form("Lc%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrLOrgnb +=  Form("Lb%s_Bin", nameMasterStting.Data());
            namePlotThnRecoCorrLOrgnX +=  Form("LX%s_Bin", nameMasterStting.Data());
        }
    }
    
    TString isDataOrMC = "";
    if(!fReadMC)isDataOrMC = "Data";
    else if(fReadMC)isDataOrMC = "MC";
    
    
    if(!fMixing && fCheckCutDist && fRecoTrk){
        
        for(Int_t i=2;  i<fNPtBins;  i++){
            
            //DPlus related QA
            TH2F* hDistPtvsd0_DauPion1  = new TH2F(Form("hDistPtvsd0_DauPion1_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "#pi_{1}: d_{0} vs. #it{p}_{T} dist.", 100, 0., 10., 50, 0., 5.);
            hDistPtvsd0_DauPion1->GetXaxis()->SetTitle("p_{T}^{#pi_{1}} distribution");
            hDistPtvsd0_DauPion1->GetYaxis()->SetTitle("d_{0}^{#pi_{1}} distribution");
            hDistPtvsd0_DauPion1->Sumw2();
            fOutput->Add(hDistPtvsd0_DauPion1);
            
            TH2F* hDistPtvsd0_DauKaon  = new TH2F(Form("hDistPtvsd0_DauKaon_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "K: d_{0} vs. #it{p}_{T} dist.", 100, 0., 10., 50, 0., 5.);
            hDistPtvsd0_DauKaon->GetXaxis()->SetTitle("p_{T}^{K} distribution");
            hDistPtvsd0_DauKaon->GetYaxis()->SetTitle("d_{0}^{K} distribution");
            hDistPtvsd0_DauKaon->Sumw2();
            fOutput->Add(hDistPtvsd0_DauKaon);
            
            TH2F* hDistPtvsd0_DauPion2  = new TH2F(Form("hDistPtvsd0_DauPion2_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "#pi_{2}: d_{0} vs. #it{p}_{T} dist.", 100, 0., 10., 50, 0., 5.);
            hDistPtvsd0_DauPion2->GetXaxis()->SetTitle("p_{T}^{#pion_{2}} distribution");
            hDistPtvsd0_DauPion2->GetYaxis()->SetTitle("d_{0}^{#pion_{2}} distribution");
            hDistPtvsd0_DauPion2->Sumw2();
            fOutput->Add(hDistPtvsd0_DauPion2);
            
            TH2F* hDistCosPAvsDCA  = new TH2F(Form("hDistCosPAvsDCA_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "cos(#theta) vs. DCA", 500, 0.5, 1., 100, 0., 0.1);
            hDistCosPAvsDCA->GetXaxis()->SetTitle("cos(#theta)");
            hDistCosPAvsDCA->GetYaxis()->SetTitle("DCA");
            hDistCosPAvsDCA->Sumw2();
            fOutput->Add(hDistCosPAvsDCA);
            
            TH1F* hDistDplusCosPA   = new TH1F(Form("hDistDplusCosPA_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "cos(#theta) dist.", 500, 0.5, 1.0);
            hDistDplusCosPA->GetXaxis()->SetTitle("cos(#theta)");
            hDistDplusCosPA->GetYaxis()->SetTitle("#entries");
            hDistDplusCosPA->Sumw2();
            fOutput->Add(hDistDplusCosPA);
            
            TH1F* hDistDplusCosPAxy   = new TH1F(Form("hDistDplusCosPAxy_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "cos(#theta_{xy}) dist.", 500, 0.5, 1.0);
            hDistDplusCosPAxy->GetXaxis()->SetTitle("cos(#theta_{xy})");
            hDistDplusCosPAxy->GetYaxis()->SetTitle("#entries");
            hDistDplusCosPAxy->Sumw2();
            fOutput->Add(hDistDplusCosPAxy);
            
            TH1F* hDistDplusDCA   = new TH1F(Form("hDistDplusDCA_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "DCA dist.", 100, 0.0, 0.10);
            hDistDplusDCA->GetXaxis()->SetTitle("DCA");
            hDistDplusDCA->GetYaxis()->SetTitle("#entries");
            hDistDplusDCA->Sumw2();
            fOutput->Add(hDistDplusDCA);
            
            TH1F* hDistDplusd0Square   = new TH1F(Form("hDistDplusd0Square_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "d_{0}^{2} dist.", 150, 0.0, 1.5);
            hDistDplusd0Square->GetXaxis()->SetTitle("d_{0}^{2}");
            hDistDplusd0Square->GetYaxis()->SetTitle("#entries");
            hDistDplusd0Square->Sumw2();
            fOutput->Add(hDistDplusd0Square);
            
            TH1F* hDistDplusSigmaVert   = new TH1F(Form("hDistDplusSigmaVert_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "#sigma_{vert} dist.", 500, 0.0, 0.10);
            hDistDplusSigmaVert->GetXaxis()->SetTitle("#sigma_{vert}");
            hDistDplusSigmaVert->GetYaxis()->SetTitle("#entries");
            hDistDplusSigmaVert->Sumw2();
            fOutput->Add(hDistDplusSigmaVert);
            
            TH1F* hDistDplusDecayLen   = new TH1F(Form("hDistDplusDecayLen_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Decay length dist.", 1500, 0.0, 1.5);
            hDistDplusDecayLen->GetXaxis()->SetTitle("Decay length");
            hDistDplusDecayLen->GetYaxis()->SetTitle("#entries");
            hDistDplusDecayLen->Sumw2();
            fOutput->Add(hDistDplusDecayLen);
            
            TH1F* hDistDplusDecayLenxy   = new TH1F(Form("hDistDplusDecayLenxy_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Decay length-XY dist.", 1000, 0.0, 1.5);
            hDistDplusDecayLenxy->GetXaxis()->SetTitle("Decay length xy");
            hDistDplusDecayLenxy->GetYaxis()->SetTitle("#entries");
            hDistDplusDecayLenxy->Sumw2();
            fOutput->Add(hDistDplusDecayLenxy);
            
            TH1F* hDistDplusNormDecayLenxy   = new TH1F(Form("hDistDplusNormDecayLenxy_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Norm:Decay length-XY dist.", 1500, 0.0, 15.);
            hDistDplusNormDecayLenxy->GetXaxis()->SetTitle("Norm Decay length xy");
            hDistDplusNormDecayLenxy->GetYaxis()->SetTitle("#entries");
            hDistDplusNormDecayLenxy->Sumw2();
            fOutput->Add(hDistDplusNormDecayLenxy);
            
            TH1F* hDistDplusDist12   = new TH1F(Form("hDistDplusDist12_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "dist1-2 dist.", 5000, 0.0, 50.0);
            hDistDplusDist12->GetXaxis()->SetTitle("dist1-2");
            hDistDplusDist12->GetYaxis()->SetTitle("#entries");
            hDistDplusDist12->Sumw2();
            fOutput->Add(hDistDplusDist12);
            
            TH1F* hDistDplusDist23   = new TH1F(Form("hDistDplusDist23_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "dist2-3 dist.", 5000, 0.0, 50.0);
            hDistDplusDist23->GetXaxis()->SetTitle("dist2-3");
            hDistDplusDist23->GetYaxis()->SetTitle("#entries");
            hDistDplusDist23->Sumw2();
            fOutput->Add(hDistDplusDist23);
            
            //DPlus QA
            TH1F* hDistDplusRap   = new TH1F(Form("hDistDplusRap_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Rapidity dist.", 180, -0.9, 0.9);
            hDistDplusRap->GetXaxis()->SetTitle("Y^{D^{+}}");
            hDistDplusRap->GetYaxis()->SetTitle("#entries");
            hDistDplusRap->Sumw2();
            fOutput->Add(hDistDplusRap);
            
            TH1F* hDistDplusEta   = new TH1F(Form("hDistDplusEta_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Eta dist.", 200, -1.0, 1.0);
            hDistDplusEta->GetXaxis()->SetTitle("#eta^{D^{+}}");
            hDistDplusEta->GetYaxis()->SetTitle("#entries");
            hDistDplusEta->Sumw2();
            fOutput->Add(hDistDplusEta);
            
            TH1F* hDistDplusPhi   = new TH1F(Form("hDistDplusPhi_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "#varphi dist.", 360, 0, 2*Pi);
            hDistDplusPhi->GetXaxis()->SetTitle("#varphi^{D^{+}}");
            hDistDplusPhi->GetYaxis()->SetTitle("#entries");
            hDistDplusPhi->Sumw2();
            fOutput->Add(hDistDplusPhi);
            
            TH1F* hDistDplusMass   = new TH1F(Form("hDistDplusMass_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Inv Mass dist.", nbins,fLowmasslimit,fUpmasslimit);
            hDistDplusMass->GetXaxis()->SetTitle("D^{+} Inv. Mass");
            hDistDplusMass->GetYaxis()->SetTitle("#entries");
            hDistDplusMass->Sumw2();
            fOutput->Add(hDistDplusMass);
            
            TH2F* hDistDplusMassVsPhi = new TH2F(Form("hDistDplusMassVsPhi_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Inv Mass vs. #varphi dist.", nbins,fLowmasslimit,fUpmasslimit,64,0,2*Pi);
            hDistDplusMassVsPhi->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
            hDistDplusMassVsPhi->GetYaxis()->SetTitle("D^{+} #varphi");
            hDistDplusMassVsPhi->Sumw2();
            fOutput->Add(hDistDplusMassVsPhi);
            
            TH2F* hDistDplusMassVsEta = new TH2F(Form("hDistDplusMassVsEta_Bin%d_%s_%s", i, isDataOrMC.Data(), fCutSuffix.Data()), "Inv Mass vs. #eta dist.", nbins,fLowmasslimit,fUpmasslimit,36,-0.9,0.9);
            hDistDplusMassVsEta->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
            hDistDplusMassVsEta->GetYaxis()->SetTitle("D^{+} #eta");
            hDistDplusMassVsEta->Sumw2();
            fOutput->Add(hDistDplusMassVsEta);
            
            TH2F* hDistDpluswoYCutVsPt = new TH2F(Form("hDistDpluswoYCutVsPt_Bin%d", i), "Rap vs. p_{T} dist.", 240, 0, 24.0, 100, -5., 5.);
            hDistDpluswoYCutVsPt->GetXaxis()->SetTitle("D^{+} Y");
            hDistDpluswoYCutVsPt->GetYaxis()->SetTitle("D^{+} p_{T}");
            hDistDpluswoYCutVsPt->Sumw2();
            fOutput->Add(hDistDpluswoYCutVsPt);
            
        }
    }
    
    
    if(!fMixing){
        //Track related QA
        TH1F* hTrackPhi = new TH1F(Form("AssoTrkPhi_%s",isDataOrMC.Data()), "Associated Trk #varphi Dist", 320, 0, 2*Pi);
        hTrackPhi->SetMarkerStyle(29);
        hTrackPhi->SetMarkerSize(0.9);
        hTrackPhi->GetXaxis()->SetTitle("Track #varphi");
        hTrackPhi->GetYaxis()->SetTitle("#entries");
        hTrackPhi->Sumw2();
        fOutput->Add(hTrackPhi);
        
        TH1F* hTrackEta = new TH1F(Form("AssoTrkEta_%s",isDataOrMC.Data()), "Associated Trk #eta Dist", 800,-1.0, 1.0);
        hTrackEta->SetMarkerStyle(29);
        hTrackEta->SetMarkerSize(0.9);
        hTrackEta->GetXaxis()->SetTitle("Track #eta");
        hTrackEta->GetYaxis()->SetTitle("#entries");
        hTrackEta->Sumw2();
        fOutput->Add(hTrackEta);
        
        TH1F* hTrackPt  = new TH1F(Form("AssoTrkPt_%s",isDataOrMC.Data()), "Associated Trk p_{T} Dist", 300,0.0,30);
        hTrackPt->SetMarkerStyle(29);
        hTrackPt->SetMarkerSize(0.9);
        hTrackPt->GetXaxis()->SetTitle("Track p_{T}");
        hTrackPt->GetYaxis()->SetTitle("#entries");
        hTrackPt->Sumw2();
        fOutput->Add(hTrackPt);
        
        Double_t MaxTrklets=100;
        if(fSystem==0)MaxTrklets=150;
        else MaxTrklets=450;
        TH2F* h2SPDTrkVsTrkMult = new TH2F("h2SPDTrkVsTrkMult", "Tracklets Vs Track Mult", 150, 0., MaxTrklets, 150, 0., MaxTrklets);
        h2SPDTrkVsTrkMult->GetXaxis()->SetTitle("SPD tracklets");
        h2SPDTrkVsTrkMult->GetYaxis()->SetTitle("Trk Multiplicity");
        h2SPDTrkVsTrkMult->SetDrawOption("SURF1");
        h2SPDTrkVsTrkMult->Sumw2();
        fOutput->Add(h2SPDTrkVsTrkMult);
    }
    
    
    for(Int_t i=2;  i<fNPtBins;  i++){
        
        if(!fReadMC){
            namePlotThnRecoWOeffMass.Resize(25);
            namePlotThnRecoMass.Resize(25);
            namePlotThnRecoCorr.Resize(25);
            namePlotThnRecoCorrL.Resize(25);
            namePlotThnRecoWOeffMass+=i;
            namePlotThnRecoMass     +=i;
            namePlotThnRecoCorr     +=i;
            namePlotThnRecoCorrL    +=i;
        }else if(fReadMC & !fMixing){
            if(fRecoTrk){
                namePlotThnRecoMassOrgnAll.Resize(25);
                namePlotThnRecoCorrOrgnc.Resize(25);
                namePlotThnRecoCorrOrgnb.Resize(25);
                namePlotThnRecoCorrOrgnX.Resize(25);
                namePlotThnRecoCorrLOrgnc.Resize(25);
                namePlotThnRecoCorrLOrgnb.Resize(25);
                namePlotThnRecoCorrLOrgnX.Resize(25);
                namePlotThnRecoMassOrgnAll+=i;
                namePlotThnRecoCorrOrgnc  +=i;
                namePlotThnRecoCorrOrgnb  +=i;
                namePlotThnRecoCorrOrgnX  +=i;
                namePlotThnRecoCorrLOrgnc +=i;
                namePlotThnRecoCorrLOrgnb +=i;
                namePlotThnRecoCorrLOrgnX +=i;
            }
        }
        
        
        if(!fMixing){
            if(!fReadMC){
                TH1F *hDplusHistoWoEff = new TH1F(namePlotThnRecoWOeffMass.Data(), "Data D^{+} Inv Mass WO Eff", nbins,fLowmasslimit,fUpmasslimit);
                hDplusHistoWoEff->Sumw2();
                fOutput->Add(hDplusHistoWoEff);
                
                THnSparseI *hDplusHistoWEff = new THnSparseI(namePlotThnRecoMass.Data(), "Data_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
                hDplusHistoWEff->Sumw2();
                fOutputCorr->Add(hDplusHistoWEff);
            }else if(fReadMC){
                if(fRecoTrk){
                    THnSparseI *hDplusHistoWoEffFrmc = new THnSparseI(namePlotThnRecoMassOrgnAll.Data(), "Reco_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
                    hDplusHistoWoEffFrmc->Sumw2();
                    fOutputCorr->Add(hDplusHistoWoEffFrmc);
                }
            }
        }
        
        
        
        // Correlations in ThnSparce
        if(!fMixing){
            if(!fReadMC){
                THnSparseI *hPhiSEThnCorr = new THnSparseI(namePlotThnRecoCorr.Data(), "Data_SE_Corr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorr->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorr);
                
                //D-Leading Particles
                THnSparseI *hPhiSEThnCorrL = new THnSparseI(namePlotThnRecoCorrL.Data(), "Data_SE_LCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;", 5, nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrL->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrL);
            }else if(fReadMC){
                THnSparseI *hPhiSEThnCorrFromc = new THnSparseI(namePlotThnRecoCorrOrgnc.Data(), "MCrc_SE_Corr_ThnSparse_Fromc; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrFromc->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrFromc);
                
                THnSparseI *hPhiSEThnCorrFromb = new THnSparseI(namePlotThnRecoCorrOrgnb.Data(), "MCrc_SE_Corr_ThnSparse_Fromb; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrFromb->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrFromb);
                
                THnSparseI *hPhiSEThnCorrFromX = new THnSparseI(namePlotThnRecoCorrOrgnX.Data(), "MCrc_SE_Corr_ThnSparse_FromX; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrFromX->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrFromX);
                
                THnSparseI *hPhiSEThnCorrLFromc = new THnSparseI(namePlotThnRecoCorrLOrgnc.Data(), "MCrc_SE_LCorr_ThnSparse_Fromc; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrLFromc->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrLFromc);
                
                THnSparseI *hPhiSEThnCorrLFromb = new THnSparseI(namePlotThnRecoCorrLOrgnb.Data(), "MCrc_SE_LCorr_ThnSparse_Fromb; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrLFromb->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrLFromb);
                
                THnSparseI *hPhiSEThnCorrLFromX = new THnSparseI(namePlotThnRecoCorrLOrgnX.Data(), "MCrc_SE_LCorr_ThnSparse_FromX; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                hPhiSEThnCorrLFromX->Sumw2();
                fOutputCorr->Add(hPhiSEThnCorrLFromX);
            }
        }else if(fMixing && fRecoTrk){
            // Histogram for Event Mixing
            namePlotThnRecoCorr+="_evMix";
            namePlotThnRecoCorrL+="_evMix";
            
            THnSparseI *hPhiMEThnCorr = new THnSparseI(namePlotThnRecoCorr.Data(), "Data_ME_Corr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
            hPhiMEThnCorr->Sumw2();
            fOutputCorr->Add(hPhiMEThnCorr);
            
            THnSparseI *hPhiMEThnCorrL = new THnSparseI(namePlotThnRecoCorrL.Data(), "Data_ME_LCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt_Pool;", 5, nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
            hPhiMEThnCorrL->Sumw2();
            fOutputCorr->Add(hPhiMEThnCorrL);
            
        }
    }
}



//_________________________________________________________________________________________________
Int_t AliAnalysisTaskSEDplusCorrelations::CheckOriginPartOfDPlus(TClonesArray* arrayMC, AliAODMCParticle *mcDplus) const {
    
    //Checking origin of D+ and retrun the PDF of the mother particle
    printf("AliAnalysisTaskSEDplusCorrelations::CheckOriginPartOfDPlus() \n");
    
    Int_t DpMomPDG = 0;
    Int_t mcDplusMomLabel = mcDplus->GetMother();
    Bool_t isFromB=kFALSE, isQuarkFound=kFALSE;
    
    while (mcDplusMomLabel > 0){
        AliAODMCParticle* mcDplusMom = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mcDplusMomLabel));
        if(mcDplusMom){
            DpMomPDG = TMath::Abs(mcDplusMom->GetPdgCode());
            if ((DpMomPDG > 500 && DpMomPDG < 600) || (DpMomPDG > 5000 && DpMomPDG < 6000))isFromB=kTRUE;
            if(DpMomPDG==4 || DpMomPDG==5) isQuarkFound=kTRUE;
            mcDplusMomLabel = mcDplusMom->GetMother();
        }else{
            break;
        }
    }
    
    if(isQuarkFound) {
        if(isFromB) return 5;
        else return 4;
    }
    else return -1;
}


//____________________| Terminate
void AliAnalysisTaskSEDplusCorrelations::Terminate(Option_t *) {
    
    fOutputCorr = dynamic_cast<TList*> (GetOutputData(2));
    if (!fOutputCorr) {
        printf("ERROR: fOutputCorr not available\n");
        return;
    }
    
    fDplusCuts = dynamic_cast<AliRDHFCutsDplustoKpipi*>(GetOutputData(3));
    if(!fDplusCuts){
        printf("ERROR: fDplusCuts(Dplus Track Cuts) are not available\n");
        return;
    }
    
    fAssoCuts = dynamic_cast<AliHFAssociatedTrackCuts*>(GetOutputData(4));
    if(!fAssoCuts){
        printf("ERROR: fAssoCuts(Associated Track Cuts) are not available\n");
        return;
    }
    
    fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(5));
    if (!fCounter) {
        printf("ERROR: fCounter not available\n");
        return;
    }
    
    return;
}
//EOF Jitendra
