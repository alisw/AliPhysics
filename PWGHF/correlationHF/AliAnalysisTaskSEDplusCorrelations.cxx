

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

/* $Id: AliAnalysisTaskSEDplusCorrelations.cxx 58883 2012-10-02 09:41:01Z prino $ */
/* AliAnalysisTask for HF(DPlus:3Prongs)-Hadron/Kaon/K0 azimuthal correlations
 By: Jitendra Kumar (edited for pool by pool on 31.01.2016)
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
fSelect(0),
fOutput(0x0),
fOutputCorr(0X0),
fReadMC(kFALSE),
fReco(kFALSE),
fMontecarlo(kFALSE),
fMCGenEvType(kFALSE),
fMixing(kFALSE),
farrayMC(0x0),
fSystem(kFALSE),
fUseBit(kTRUE),
fTCconfig(kFALSE),
fHistNEvents(0),
fHistNDplus(0),
fCounter(0x0),
fDplusCuts(0),
fAssoCuts(0),
fCorrelator(0x0),
fNPtBins(0),
fBinWidth(0),
fCentrOrMult(-99),
fMultiplicity(1),
fEffTrack(kFALSE),
fEffDplus(kFALSE),
fCentralityEstimator(0),
fEvalCentrality(kFALSE),
fMinCentrality(0),
fMaxCentrality(100),
fPoolByPool(kFALSE),
fWhichPool(0),
fCheckCutDist(kFALSE)
{
    // Default constructor
}


//____________________| Specific Constructor
AliAnalysisTaskSEDplusCorrelations::AliAnalysisTaskSEDplusCorrelations(const Char_t* name,AliRDHFCutsDplustoKpipi* DplusCuts, AliHFAssociatedTrackCuts *AsscCuts) :
AliAnalysisTaskSE(name),
fSelect(0),
fOutput(0x0),
fOutputCorr(0X0),
fReadMC(kFALSE),
fReco(kFALSE),
fMontecarlo(kFALSE),
fMCGenEvType(kFALSE),
fMixing(kFALSE),
farrayMC(0x0),
fSystem(kFALSE),
fUseBit(kTRUE),
fTCconfig(kFALSE),
fHistNEvents(0),
fHistNDplus(0),
fCounter(0x0),
fDplusCuts(0),
fAssoCuts(AsscCuts),
fCorrelator(0x0),
fNPtBins(0),
fBinWidth(0.002),
fCentrOrMult(-99),
fMultiplicity(1),
fEffTrack(kFALSE),
fEffDplus(kFALSE),
fCentralityEstimator(0),
fEvalCentrality(kFALSE),
fMinCentrality(0),
fMaxCentrality(100),
fPoolByPool(kFALSE),
fWhichPool(0),
fCheckCutDist(kFALSE)
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
fSelect(source.fSelect),
fOutput(source.fOutput),
fOutputCorr(source.fOutputCorr),
fReadMC(source.fReadMC),
fReco(source.fReadMC),
fMontecarlo(source.fMontecarlo),
fMCGenEvType(source.fMCGenEvType),
fMixing(source.fMixing),
farrayMC(source.farrayMC),
fSystem(source.fMixing),
fUseBit(source.fUseBit),
fTCconfig(source.fTCconfig),
fHistNEvents(source.fHistNEvents),
fHistNDplus(source.fHistNDplus),
fCounter(source.fCounter),
fDplusCuts(source.fDplusCuts),
fAssoCuts(source.fAssoCuts),
fCorrelator(source.fCorrelator),
fNPtBins(source.fNPtBins),
fBinWidth(source.fBinWidth),
fCentrOrMult(source.fCentrOrMult),
fMultiplicity(source.fMultiplicity),
fEffTrack(source.fEffTrack),
fEffDplus(source.fEffDplus),
fCentralityEstimator(source.fCentralityEstimator),
fEvalCentrality(source.fEvalCentrality),
fMinCentrality(source.fMinCentrality),
fMaxCentrality(source.fMaxCentrality),
fPoolByPool(source.fPoolByPool),
fWhichPool(source.fWhichPool),
fCheckCutDist(source.fCheckCutDist)
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
    fSelect = orig.fSelect;
    fOutput = orig.fOutput;
    fOutputCorr = orig.fOutputCorr;
    fReadMC = orig.fReadMC;
    fReco = orig.fReco;
    fMontecarlo = orig.fMontecarlo;
    fMCGenEvType = orig.fMCGenEvType;
    fMixing = orig.fMixing;
    farrayMC = orig.farrayMC;
    fSystem = orig.fSystem;
    fUseBit = orig.fUseBit;
    fTCconfig=orig.fTCconfig;
    fHistNEvents = orig.fHistNEvents;
    fHistNDplus = orig.fHistNDplus;
    fCounter = orig.fCounter;
    fDplusCuts = orig.fDplusCuts;
    fAssoCuts = orig.fAssoCuts;
    fCorrelator = orig.fCorrelator;
    fNPtBins = orig.fNPtBins;
    fBinWidth = orig.fBinWidth;
    fCentrOrMult = orig.fCentrOrMult;
    fMultiplicity = orig.fMultiplicity;
    fEffTrack = orig.fEffTrack;
    fEffDplus = orig.fEffDplus;
    fCentralityEstimator = orig.fCentralityEstimator;
    fEvalCentrality=orig.fEvalCentrality;
    fMinCentrality=orig.fMinCentrality;
    fMaxCentrality=orig.fMaxCentrality;
    fPoolByPool=orig.fPoolByPool;
    fWhichPool=orig.fWhichPool;
    fCheckCutDist=orig.fCheckCutDist;
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
        if(fCentralityEstimator==fDplusCuts->GetUseCentrality())isCentEstimatorOk = kTRUE; //
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
    
    fHistNEvents = new TH1F("fHistNEvents", "Events Stats", 9, -0.5 , 8.5);
    fHistNEvents->GetXaxis()->SetBinLabel(1,"All");
    fHistNEvents->GetXaxis()->SetBinLabel(2,"Pileup Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(3,"Cent-Mult Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(4,"Bad-Vtx Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(5,"Trigger Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(6,"Zvtx Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(7,"PS Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(8,"HFCorr Interface Rej");
    fHistNEvents->GetXaxis()->SetBinLabel(9,"Accepted");
    fHistNEvents->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistNEvents->SetFillColor(kBlue);
    fHistNEvents->SetFillStyle(3015);
    fHistNEvents->SetMinimum(0);
    //fHistNEvents->Sumw2();
    fOutput->Add(fHistNEvents);
    
    fHistNDplus = new TH1F("fHistNDplus", "Dplus Stats ", 8, -0.5 , 7.5);
    fHistNDplus->GetXaxis()->SetBinLabel(1,"All");
    fHistNDplus->GetXaxis()->SetBinLabel(2,"-pT Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(3,"<2GeV Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(4,"SelBit Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(5,"Topocut Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(6,"TC Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(7,"FiduAcc Rej");
    fHistNDplus->GetXaxis()->SetBinLabel(8,"Accepted");
    fHistNDplus->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistNDplus->SetMinimum(0);
    fHistNDplus->SetFillColor(kRed);
    fHistNDplus->SetFillStyle(3015);
    //fHistNDplus->Sumw2();
    fOutput->Add(fHistNDplus);
    
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
    fCorrelator->SetAssociatedParticleType(fSelect);
    //fCorrelator->SetApplyDisplacementCut(fDisplacement); //set kFALSE/kTRUE for using the displacement cut
    fCorrelator->SetUseReco(fReco);
    fCorrelator->SetUseMC(fMontecarlo);
    
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
    if(!aod) return;
    
    TClonesArray *array3Prong = 0;
    array3Prong=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    if(!aod && AODEvent() && IsStandardAOD()) {
        
        aod = dynamic_cast<AliAODEvent*> (AODEvent());
        AliAODHandler* aodHandler = (AliAODHandler*)((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
        
        if(aodHandler->GetExtensions()) {
            AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
            AliAODEvent *aodFromExt = ext->GetAOD();
            array3Prong=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
        }
    }
    else if(!aod || !array3Prong){
        printf("AliAnalysisTaskSEDplusCorrelations::UserExec: AOD  Charm3Prong branch not found!\n");
        return;
    }
    
    if(!aod->GetPrimaryVertex()||TMath::Abs(aod->GetMagneticField())<0.001) return;
    
    fCounter->StoreEvent(aod,fDplusCuts,fMontecarlo);
    fHistNEvents->Fill(0);
    
    //MC Kinemactics
    if(fReadMC && fMontecarlo){
        
        farrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
        if(!farrayMC){
            printf("AliAnalysisTaskSEDplusCorrelations::UserExec: MC particles branch not found!\n");
            return;
        }
        fCorrelator->SetMCArray(farrayMC);
        
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
                // ((TH1D*)fOutputBasic->FindObject("EventTypeMC"))->Fill(eventType);
            }
            if(NMCevents && !isMCeventgood){
                if(fDebug>2) std::cout << "The MC event " << eventType << " not interesting for this analysis: skipping" << std::endl;
                return;
            }
        }
    }
    
    //Event Selection
    if(!fReadMC && !fDplusCuts->IsEventSelected(aod)){
        if(fDplusCuts->GetWhyRejection()==1)fHistNEvents->Fill(1);//pileup
        if(fDplusCuts->GetWhyRejection()==5)fHistNEvents->Fill(4);//triggers
        if(fDplusCuts->GetWhyRejection()==6)fHistNEvents->Fill(5);//Zvtx
        if(fDplusCuts->GetWhyRejection()==7)fHistNEvents->Fill(6);//PS
        return;
    }else if(fReadMC){
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
        //fMultiplicity = count;
        //fCentrOrMult = fMultiplicity;
        fCentrOrMult = (Double_t)(AliVertexingHFUtils::GetNumberOfTrackletsInEtaRange(aod,-1.,1.));
        ((TH2F*)fOutput->FindObject("h2SPDTrkVsTrkMult"))->Fill(aod->GetNumberOfTracks(),fCentrOrMult);
        //if(fCentrOrMult<fMinMult || fCentrOrMult >fMaxMult)fHistNEvents->Fill(2);
    }
    
    AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
    Bool_t isGoodVtx=kFALSE;
    TString primTitle = vtx1->GetTitle();
    if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0)isGoodVtx=kTRUE;
    if(!isGoodVtx){
        fHistNEvents->Fill(3);
        return;
    }
    
    
    //Pool by Pool Seetting
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
    Bool_t    isDplus = kFALSE;
    Int_t pdgDgDplustoKpipi[3]={321,211,211};
    
    //printf("Number of D+->Kpipi: %d and of tracks: %d\n",nDplusKpipi,aod->GetNumberOfTracks());
    for (Int_t iDplusKpipi = 0; iDplusKpipi < nDplusKpipi; iDplusKpipi++){
        
        // D+ Primary Vertex Setting
        AliAODRecoDecayHF3Prong *d = (AliAODRecoDecayHF3Prong*)array3Prong->UncheckedAt(iDplusKpipi); // D+
        if(!d)continue;
        fHistNDplus->Fill(0);
        
        if(d->Pt()<0.)fHistNDplus->Fill(1);
        
        if(d->Pt()<2.){
            fHistNDplus->Fill(2);
            continue;
        }
        
        if(fUseBit && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)){
            fHistNDplus->Fill(3);
            continue;
        }
        
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
        if(!fDplusCuts->IsInFiducialAcceptance(ptDplusCand,rapDplusCand)) {
            //Vertexing cleaning...
            if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
            if(unsetvtx) d->UnsetOwnPrimaryVtx();
            fHistNDplus->Fill(6);
            continue;
        }
        if(fMontecarlo){
            labDp = d->MatchToMC(411,farrayMC,3,pdgDgDplustoKpipi);
            if(labDp>=0){
                isDplus = kTRUE;
            }
        }
        fHistNDplus->Fill(7);
        
        //Computing Correlations: Below Funtions
        HadronCorrelations(d,farrayMC,isDplus);
        
        //Vertexing cleaning...
        if(recVtx)fDplusCuts->CleanOwnPrimaryVtx(d,aod,origownvtx);
        if(unsetvtx) d->UnsetOwnPrimaryVtx();
        
    }// loop over D+
    
    
    //EVT MIXING PART
    if(fMixing ) {
        Bool_t updated = fCorrelator->PoolUpdate();
        if(!updated) AliInfo("Pool was not updated");
    }
    
    
    // Posting Slots
    PostData(1,fOutput);
    PostData(2,fOutputCorr);
    PostData(5,fCounter);
    return;
    
}



//____________________| Calculating Correlation
void AliAnalysisTaskSEDplusCorrelations::HadronCorrelations(AliAODRecoDecayHF3Prong* d,  TClonesArray *arrayMC, Bool_t isDplus) {
    
    //Switch for String for Hist name
    TString parttype = "", datatype = "";
    if(fSelect==1)parttype += "Hdron";
    else if(fSelect==2)parttype += "Kaons";
    else if(fSelect==3)parttype += "kZero";
    
    if(fReco){
        if(!fReadMC)datatype += "Data";
        else if(fReadMC)datatype += "Reco";
    }
    
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
    
    
    ((TH2F*)fOutput->FindObject("DplusMVsPhi"))->Fill(d->InvMassDplus(),d->Phi());
    ((TH2F*)fOutput->FindObject("DplusMVsPt"))->Fill(d->InvMassDplus(),d->Pt());
    ((TH2F*)fOutput->FindObject("DplusMVsEta"))->Fill(d->InvMassDplus(),d->Eta());
    
    if(fCheckCutDist && !fMixing){
    	((TH1F*)fOutput->FindObject(Form("hPtDaughterPion1_Bin%d",iPtBin)))->Fill(d->PtProng(0)); //Prong (0)= Pion1 from line 133 of AliAODRecoDecayHF3Prong.cxx
        ((TH1F*)fOutput->FindObject(Form("hPtDaughterKaon_Bin%d",iPtBin)))->Fill(d->PtProng(1)); //Prong (1)= Kaon from line 133 of AliAODRecoDecayHF3Prong.cxx
        ((TH1F*)fOutput->FindObject(Form("hPtDaughterPion2_Bin%d",iPtBin)))->Fill(d->PtProng(2)); //Prong (2)= Pion2 from line 133 of AliAODRecoDecayHF3Prong.cxx
        ((TH1F*)fOutput->FindObject(Form("hd0DaughterPion1_Bin%d",iPtBin)))->Fill(d->Getd0Prong(0)*10.0);
        ((TH1F*)fOutput->FindObject(Form("hd0DaughterKaon_Bin%d",iPtBin)))->Fill(d->Getd0Prong(1)*10.0);
        ((TH1F*)fOutput->FindObject(Form("hd0DaughterPion2_Bin%d",iPtBin)))->Fill(d->Getd0Prong(2)*10.0);
        ((TH1F*)fOutput->FindObject(Form("hdist12_Bin%d",iPtBin)))->Fill(d->GetDist12toPrim()*10.0);
        ((TH1F*)fOutput->FindObject(Form("hdist23_Bin%d",iPtBin)))->Fill(d->GetDist23toPrim()*10.0);
        ((TH1F*)fOutput->FindObject(Form("hCosPA_Bin%d",iPtBin)))->Fill(d->CosPointingAngle());
        ((TH1F*)fOutput->FindObject(Form("hCosPAXY_Bin%d",iPtBin)))->Fill(d->CosPointingAngleXY());
        ((TH1F*)fOutput->FindObject(Form("hNDeacyLenXY_Bin%d",iPtBin)))->Fill(d->NormalizedDecayLengthXY()*10.0);
        ((TH1F*)fOutput->FindObject(Form("hDecayLen_Bin%d",iPtBin)))->Fill(d->DecayLength()*10.0);
        ((TH1F*)fOutput->FindObject(Form("hDecayLenXY_Bin%d",iPtBin)))->Fill(d->DecayLengthXY()*10.0);
        ((TH2F*)fOutput->FindObject(Form("hCosPAvsdPS_Bin%d",iPtBin)))->Fill(d->DecayLength()*10.0,d->CosPointingAngle());
        ((TH2F*)fOutput->FindObject(Form("hd0DaughterKaonvsPtK_Bin%d",iPtBin)))->Fill(d->PtProng(1),d->Getd0Prong(1)*10.0);
        ((TH2F*)fOutput->FindObject(Form("hd0DaughterPion1vsPtPi1_Bin%d",iPtBin)))->Fill(d->PtProng(0),d->Getd0Prong(0)*10.0);
        ((TH2F*)fOutput->FindObject(Form("hd0DaughterPion2vsPtPi2_Bin%d",iPtBin)))->Fill(d->PtProng(2),d->Getd0Prong(2)*10.0);

        Double_t MaxDCA = -9999.;
        for(Int_t i=0; i<3; i++){
            if(d->GetDCA(i)>MaxDCA)MaxDCA=d->GetDCA(i);
        }
        ((TH1F*)fOutput->FindObject(Form("hDCA_Bin%d",iPtBin)))->Fill(MaxDCA*10.0);
        
        Double_t d0Square = d->Getd0Prong(0)*d->Getd0Prong(0)+d->Getd0Prong(1)*d->Getd0Prong(1)+d->Getd0Prong(2)*d->Getd0Prong(2);
        ((TH1F*)fOutput->FindObject(Form("hd0square_Bin%d",iPtBin)))->Fill(d0Square*100.0);
        ((TH1F*)fOutput->FindObject(Form("hSigmaVert_Bin%d",iPtBin)))->Fill(d->GetSigmaVert()*10.0);
        ((TH1F*)fOutput->FindObject(Form("hDrapidty_Bin%d",iPtBin)))->Fill(d->YDplus());
    }
    
    //Correlation with SE or ME
    Double_t nSparceCorrDplusinfo[1] = {mDplus};
    if(!fMixing){
        if(fReco){ // data/MC-Reco
            ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseM_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrDplusinfo,1.0/effDplus);
            if(fEffDplus){
                ((TH1F*)fOutput->FindObject(Form("hDplsMOrG_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(mDplus); //Dplus efficiency
            }
        }
        if(fMontecarlo){
            ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseM_%s_%s_Trth_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrDplusinfo);
        }
    }
    
    //Dplus parameter resetting and storage
    Double_t phiDplus = fCorrelator->SetCorrectPhiRange(phiDplusCand);
    fCorrelator->SetTriggerParticleProperties(ptDplusCand,phiDplus,etaDplusCand);
    
    //Pool Setting and Pool events
    Bool_t execPool = fCorrelator->ProcessEventPool();
    if(fMixing && !execPool) {
        AliInfo("Mixed event analysis: pool is not ready");
    }
    
    Int_t NofEventsinPool = 1; // SE
    if(fMixing) {
        NofEventsinPool = fCorrelator->GetNofEventsInPool();
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
            
            if(!fMixing) {
                //SE only
                ((TH1F*)fOutput->FindObject("AssoTrkPhi"))->Fill(redpart->Phi());
                ((TH1F*)fOutput->FindObject("AssoTrkPt"))->Fill(redpart->Pt());
                ((TH1F*)fOutput->FindObject("AssoTrkEta"))->Fill(redpart->Eta());
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
            
            //filling correlations in ThnSparce format
            if(fReco){
                Bool_t* partSource = NULL;
                CorrelationNSparsePlots(d, redpart, iPtBin, partSource, 1.0/effTotal);
                delete[] partSource;
            }
            
            if(fMontecarlo && isDplus){
                Bool_t* partSource = fAssoCuts->IsMCpartFromHF(label,arrayMC); // check source for 1/2/3
                CorrelationNSparsePlots(d, redpart, iPtBin, partSource, 1.0/effTotal);
                delete[] partSource;
            } // MC Gen Origin ends
        }//track loop end here
        //Leading particle correlations
        
        //pT bins setting
        Double_t ptLim_Sparse = 0.0;
        if(!fMixing) {
            if(fReco)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
            if(fMontecarlo)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        }
        else if(fMixing){
            ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin%d_evMix", parttype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        }
        if(ptleadHadron > ptLim_Sparse) ptleadHadron = ptLim_Sparse-0.01; //filling all above pT in last bin

        
        Double_t nSparceCorrLeadPart[5] = {mDplus, DeltaphiLead, DeltaetaLead, ptleadHadron, fWhichPool+0.5};
        if(!fMixing){
            if(fReco)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
            if(fMontecarlo)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
            //montecarlo temp: need fixs
        }
        
        if(fMixing) {
            ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseL_%s_%s_Bin%d_evMix", parttype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorrLeadPart);
        }
    }  // event loop.
    
    
}


//____________________|Correlation with ThnSparse Histograms
void AliAnalysisTaskSEDplusCorrelations::CorrelationNSparsePlots(AliAODRecoDecayHF3Prong *d, AliReducedParticle* track, Int_t iPtBin, Bool_t *origDplus, Double_t weightEff) {
    
    iPtBin = fDplusCuts->PtBin(d->Pt());
    Double_t mDplus = d->InvMassDplus();
    Double_t deltaPhi = fCorrelator->GetDeltaPhi();
    Double_t deltaEta = fCorrelator->GetDeltaEta();
    Double_t ptTrack = track->Pt();
    
    TString partype = "", datatype = "";
    if(fSelect==1)partype += "Hdron";
    else if(fSelect==2)partype += "Kaons";
    else if(fSelect==3)partype += "kZero";
    
    
    if(fReco){
        if(!fReadMC)datatype += "Data";
        else if(fReadMC)datatype += "Reco";
    }
    
    if(fMontecarlo){
        origDplus = new Bool_t[4];
        if(origDplus[0] && origDplus[2]) datatype += "Frmc"; // is from charm ->D
        else if(origDplus[1] && origDplus[2])datatype += "Frmb"; // is from beauty ->D
        else if(origDplus[1] && origDplus[3])datatype += "FrmB"; // is from beauty ->B
        else   datatype += "Trth";
    }
    
    //pT bins setting
    Double_t ptLim_Sparse = 0.0;
    if(!fMixing) {
        if(fReco)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d",partype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        if(fMontecarlo)ptLim_Sparse = ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d",partype.Data(), datatype.Data(),iPtBin)))->GetAxis(3)->GetXmax();
        //montecarlo temp: need fixs
    }
    else if(fMixing){
        ptLim_Sparse=((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d_evMix",partype.Data(), datatype.Data(), iPtBin)))->GetAxis(3)->GetXmax();
    }
    if(ptTrack > ptLim_Sparse) ptTrack = ptLim_Sparse-0.01; //filling all above pT in last bin
    
    
    
    //Correlations
    Double_t nSparceCorr[5] = {mDplus, deltaPhi, deltaEta, ptTrack, fWhichPool+0.5};
    if(!fMixing) {
        if(fReco)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d",partype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
        if(fMontecarlo)((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d",partype.Data(), datatype.Data(),iPtBin)))->Fill(nSparceCorr,weightEff);
    }
    else if(fMixing) {
        ((THnSparseI*)fOutputCorr->FindObject(Form("hnSparseC_%s_%s_Bin%d_evMix",partype.Data(), datatype.Data(), iPtBin)))->Fill(nSparceCorr);
    }
    
    delete[] origDplus;
    return;
    
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
    
    //Other QA -- more to be added
    TH2F* hDplusMVsPhi = new TH2F("DplusMVsPhi", "D^{+} Inv Mass Dist; D^{+} InvMass; Phi;",nbins,fLowmasslimit,fUpmasslimit,32,0,2*Pi);
    TH2F* hDplusMVsEta = new TH2F("DplusMVsEta", "D^{+} Inv Mass Dist; D^{+} InvMass; Eta;",nbins,fLowmasslimit,fUpmasslimit,40,-2.0,2.0);
    TH2F* hDplusMVsPt  = new TH2F("DplusMVsPt",  "D^{+} Inv Mass Dist; D^{+} InvMass; Pt(GeV/c);",nbins,fLowmasslimit,fUpmasslimit,28,2,16);
    TH1F* hTrackPhi = new TH1F("AssoTrkPhi", "Associated Trk Phi Dist", 32, 0, 2*Pi);
    TH1F* hTrackEta = new TH1F("AssoTrkEta", "Associated Trk Eta Dist", 36,-0.9, 0.9);
    TH1F* hTrackPt  = new TH1F("AssoTrkPt",  "Associated Trk Pt Dist", 80,0.0,20);
    hDplusMVsPhi->Sumw2();
    hDplusMVsEta->Sumw2();
    hDplusMVsPt->Sumw2();
    
    hTrackPhi->SetMarkerStyle(29);
    hTrackPhi->SetMarkerSize(0.9);
    hTrackPhi->Sumw2();
    hTrackEta->SetMarkerStyle(29);
    hTrackEta->SetMarkerSize(0.9);
    hTrackEta->Sumw2();
    hTrackPt->SetMarkerStyle(29);
    hTrackPt->SetMarkerSize(0.9);
    hTrackPt->Sumw2();
    fOutput->Add(hDplusMVsPhi);
    fOutput->Add(hDplusMVsEta);
    fOutput->Add(hDplusMVsPt);
    fOutput->Add(hTrackPhi);
    fOutput->Add(hTrackEta);
    fOutput->Add(hTrackPt);
    
    Double_t MaxTrklets=100;
    if(fSystem==0)MaxTrklets=150;
    else MaxTrklets=2000;
    TH2F* h2SPDTrkVsTrkMult = new TH2F("h2SPDTrkVsTrkMult", "Tracklets Vs Track Mult", 200, 0., MaxTrklets, 200, 0., MaxTrklets);
    h2SPDTrkVsTrkMult->SetDrawOption("SURF1");
    h2SPDTrkVsTrkMult->Sumw2();
    fOutput->Add(h2SPDTrkVsTrkMult);
    
    TString nameMasterStting   = "";
    TString namePlotThnDplus   = "hnSparse" ;
    TString namePlotThnDplusOrg= "hDplsMOr" ;
    TString namePlotThnDplusMC = "hnSparse" ;
    TString namePlotThn        = "hnSparse" ;
    TString namePlotThnL       = "hnSparse" ;
    TString namePlotThnMC      = "hnSparse" ;
    TString namePlotThnLMC     = "hnSparse" ;
    TString namePlotThnMCc     = "hnSparse" ;
    TString namePlotThnMCb     = "hnSparse" ;
    TString namePlotThnMCB     = "hnSparse" ;
    
    if(fSelect==1)nameMasterStting += "_Hdron";
    else  if(fSelect==2)nameMasterStting += "_Kaons";
    else  if(fSelect==3)nameMasterStting += "_kZero";
    else  nameMasterStting = "_Nulls";
    
    if(fReco){
        if(!fReadMC)nameMasterStting +="_Data";
        if(fReadMC)nameMasterStting +="_Reco";
        namePlotThnDplus    += Form("M%s_Bin", nameMasterStting.Data());
        namePlotThnDplusOrg += Form("G%s_Bin", nameMasterStting.Data());
        namePlotThn         += Form("C%s_Bin", nameMasterStting.Data());
        namePlotThnL        += Form("L%s_Bin", nameMasterStting.Data());
    }
    
    if(fMontecarlo){
        namePlotThnDplusMC  +=  Form("M%s_Trth_Bin", nameMasterStting.Data());
        namePlotThnMCc      +=  Form("C%s_Frmc_Bin", nameMasterStting.Data());
        namePlotThnMCb      +=  Form("C%s_Frmb_Bin", nameMasterStting.Data());
        namePlotThnMCB      +=  Form("C%s_FrmB_Bin", nameMasterStting.Data());
        namePlotThnLMC      +=  Form("L%s_Trth_Bin", nameMasterStting.Data());
    }
    
    
    for(Int_t i=0;  i<fNPtBins;  i++){
        
        namePlotThnDplus.Resize(24);
        namePlotThnDplusOrg.Resize(24);
        namePlotThnDplusMC.Resize(24);
        
        namePlotThnDplus    +=  i;
        namePlotThnDplusOrg +=  i;
        namePlotThnDplusMC  +=  i;
        
        if(!fMixing){
            if(fReco){
                TH1F *hDplusHistowoEff = new TH1F(namePlotThnDplusOrg.Data(), "Reco D^{+} Inv Mass WOEff", nbins,fLowmasslimit,fUpmasslimit);
                hDplusHistowoEff->Sumw2();
                fOutput->Add(hDplusHistowoEff);
                
                THnSparseI *DplusHisto = new THnSparseI(namePlotThnDplus.Data(), "Reco_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
                DplusHisto->Sumw2();
                fOutputCorr->Add(DplusHisto);
            }
            if(fMontecarlo){
                THnSparseI *DplusHistoMC = new THnSparseI(namePlotThnDplusMC.Data(), "MCGen_Dplus_Mass_Histo; D^{+} Inv Mass;",1,nBinsDinfo,binMinDinfo,binMaxDinfo);
                DplusHistoMC->Sumw2();
                fOutputCorr->Add(DplusHistoMC);
            }
        }
        
        if(fCheckCutDist){
            TH1F* hPtDauKaon   = new TH1F(Form("hPtDaughterKaon_Bin%d", i), "Dist. Pt Daughter-Kaon (GeV/c)", 100, 0., 50.);
            TH1F* hPtDauPion1  = new TH1F(Form("hPtDaughterPion1_Bin%d", i), "Dist. Pt Daughter-Pion1 (GeV/c)", 100, 0., 50.);
            TH1F* hPtDauPion2  = new TH1F(Form("hPtDaughterPion2_Bin%d", i), "Dist. Pt Daughter-Pion2 (GeV/c)", 100, 0., 50.);
            TH1F* hd0DauKaon   = new TH1F(Form("hd0DaughterKaon_Bin%d", i), "Dist. d0 Daughter-Kaon (mm)", 100, 0., 50.);
            TH1F* hd0DauPion1  = new TH1F(Form("hd0DaughterPion1_Bin%d", i), "Dist. d0 Daughter-Pion1 (mm)", 100, 0., 50.);
            TH1F* hd0DauPion2  = new TH1F(Form("hd0DaughterPion2_Bin%d", i), "Dist. d0 Daughter-Pion2 (mm)", 100, 0., 50.);
            TH1F* hdist12      = new TH1F(Form("hdist12_Bin%d", i), "Dist. Prim-dist12", 100, 0.0, 50.);
            TH1F* hdist23      = new TH1F(Form("hdist23_Bin%d", i), "Dist. Prim-dist23", 100, 0.0, 50.);
            TH1F* hSigmaVert   = new TH1F(Form("hSigmaVert_Bin%d", i), "Dist. Sigma Vtx (mm)", 100, 0., 50.0);
            TH1F* hCosPA       = new TH1F(Form("hCosPA_Bin%d", i), "Dist. CosTheta PA", 100, 0.0, 1.);
            TH1F* hCosPAXY     = new TH1F(Form("hCosPAXY_Bin%d", i), "Dist. CosTheta PAXY", 100, 0.0, 1.);
            TH1F* hNDeacyLen   = new TH1F(Form("hNDeacyLenXY_Bin%d", i), "Dist. Normalized decay length (mm)", 100, -100., 100.);
            TH1F* hDecayLen    = new TH1F(Form("hDecayLen_Bin%d", i), "Dist. Decay length (mm)", 100, 0., 100.0);
            TH1F* hDecayLenxy  = new TH1F(Form("hDecayLenXY_Bin%d", i), "Dist. Decay length XY (mm)", 100, 0., 100.);
            TH1F* hDCA         = new TH1F(Form("hDCA_Bin%d", i), "Dist. Max DCA (mm)", 100, 0.,50.0);
            TH1F* hd0Square    = new TH1F(Form("hd0square_Bin%d", i), "Dist. d0Square length (mm^{2})", 100, 0., 100.0);
            TH1F* hDrapidity   = new TH1F(Form("hDrapidty_Bin%d", i), "Dist. D Rapidity", 100, -5., 5.);
            TH2F* hCosPAvsdPS  = new TH2F(Form("hCosPAvsdPS_Bin%d", i), "Dist. CosPA vs Decay Length", 100, 0., 50.,100, 0.0, 1.0);
            TH2F* hd0KvsPtK    = new TH2F(Form("hd0DaughterKaonvsPtK_Bin%d", i), "Dist. d0 Daughter-Kaon (mm) vs Pt K (GeV/c)", 100, 0., 50.,100, 0.0, 50.0);
            TH2F* hd0Pi1vsPtPi1  = new TH2F(Form("hd0DaughterPion1vsPtPi1_Bin%d", i), "Dist. d0 Daughter-Pion1 (mm) vs Pt Pi1 (GeV/c)", 100, 0., 50.,100, 0.0, 50.0);
            TH2F* hd0Pi2vsPtPi2  = new TH2F(Form("hd0DaughterPion2vsPtPi2_Bin%d", i), "Dist. d0 Daughter-Pion2 (mm) vs Pt Pi2 (GeV/c)", 100, 0., 50.,100, 0.0, 50.0);
            hPtDauKaon->SetMarkerStyle(29);
            hPtDauKaon->SetMarkerSize(0.9);
            hPtDauKaon->Sumw2();
            hPtDauPion1->SetMarkerStyle(29);
            hPtDauPion1->SetMarkerSize(0.9);
            hPtDauPion1->Sumw2();
            hPtDauPion2->SetMarkerStyle(29);
            hPtDauPion2->SetMarkerSize(0.9);
            hPtDauPion2->Sumw2();
            hd0DauKaon->SetMarkerStyle(29);
            hd0DauKaon->SetMarkerSize(0.9);
            hd0DauKaon->Sumw2();
            hd0DauPion1->SetMarkerStyle(29);
            hd0DauPion1->SetMarkerSize(0.9);
            hd0DauPion1->Sumw2();
            hd0DauPion2->SetMarkerStyle(29);
            hd0DauPion2->SetMarkerSize(0.9);
            hd0DauPion2->Sumw2();
            hdist12->SetMarkerStyle(29);
            hdist12->SetMarkerSize(0.9);
            hdist12->Sumw2();
            hdist23->SetMarkerStyle(29);
            hdist23->SetMarkerSize(0.9);
            hdist23->Sumw2();
            hSigmaVert->SetMarkerStyle(29);
            hSigmaVert->SetMarkerSize(0.9);
            hSigmaVert->Sumw2();
            hCosPA->SetMarkerStyle(29);
            hCosPA->SetMarkerSize(0.9);
            hCosPA->Sumw2();
            hCosPAXY->SetMarkerStyle(29);
            hCosPAXY->SetMarkerSize(0.9);
            hCosPAXY->Sumw2();
            hNDeacyLen->SetMarkerStyle(29);
            hNDeacyLen->SetMarkerSize(0.9);
            hNDeacyLen->Sumw2();
            hDecayLen->SetMarkerStyle(29);
            hDecayLen->SetMarkerSize(0.9);
            hDecayLen->Sumw2();
            hDecayLenxy->SetMarkerStyle(29);
            hDecayLenxy->SetMarkerSize(0.9);
            hDecayLenxy->Sumw2();
            hDCA->SetMarkerStyle(29);
            hDCA->SetMarkerSize(0.9);
            hDCA->Sumw2();
            hd0Square->SetMarkerStyle(29);
            hd0Square->SetMarkerSize(0.9);
            hd0Square->Sumw2();
            hDrapidity->SetMarkerStyle(29);
            hDrapidity->SetMarkerSize(0.9);
            hDrapidity->Sumw2();
            hCosPAvsdPS->SetMarkerStyle(29);
            hCosPAvsdPS->SetMarkerSize(0.9);
            hCosPAvsdPS->Sumw2();
            hd0KvsPtK->SetMarkerStyle(29);
            hd0KvsPtK->SetMarkerSize(0.9);
            hd0KvsPtK->Sumw2();
            hd0Pi1vsPtPi1->SetMarkerStyle(29);
            hd0Pi1vsPtPi1->SetMarkerSize(0.9);
            hd0Pi1vsPtPi1->Sumw2();
            hd0Pi2vsPtPi2->SetMarkerStyle(29);
            hd0Pi2vsPtPi2->SetMarkerSize(0.9);
            hd0Pi2vsPtPi2->Sumw2();
            fOutput->Add(hPtDauKaon);
            fOutput->Add(hPtDauPion1);
            fOutput->Add(hPtDauPion2);
            fOutput->Add(hd0DauKaon);
            fOutput->Add(hd0DauPion1);
            fOutput->Add(hd0DauPion2);
            fOutput->Add(hdist12);
            fOutput->Add(hdist23);
            fOutput->Add(hSigmaVert);
            fOutput->Add(hCosPA);
            fOutput->Add(hCosPAXY);
            fOutput->Add(hNDeacyLen);
            fOutput->Add(hDecayLen);
            fOutput->Add(hDecayLenxy);
            fOutput->Add(hDCA);
            fOutput->Add(hd0Square);
            fOutput->Add(hDrapidity);
            fOutput->Add(hCosPAvsdPS);
            fOutput->Add(hd0KvsPtK);
            fOutput->Add(hd0Pi1vsPtPi1);
            fOutput->Add(hd0Pi2vsPtPi2);
        }
        
        namePlotThn.Resize(24);
        namePlotThnL.Resize(24);
        namePlotThnMC.Resize(24);
        namePlotThnLMC.Resize(24);
        namePlotThnMCc.Resize(24);
        namePlotThnMCb.Resize(24);
        namePlotThnMCB.Resize(24);
        
        namePlotThn        +=  i;
        namePlotThnL       +=  i;
        namePlotThnMC      +=  i;
        namePlotThnLMC     +=  i;
        namePlotThnMCc     +=  i;
        namePlotThnMCb     +=  i;
        namePlotThnMCB     +=  i;
        
            // Correlations in ThnSparce
            if(!fMixing) {
                if(fReco){//D-h
                    THnSparseI *hPhiCorr = new THnSparseI(namePlotThn.Data(), "Reco_SE_Corr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hPhiCorr->Sumw2();
                    fOutputCorr->Add(hPhiCorr);
                    
                    //D-Leading Particles
                    THnSparseI *hCorrLead = new THnSparseI(namePlotThnL.Data(), "Reco_SE_LeadingCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;", 5, nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hCorrLead->Sumw2();
                    fOutputCorr->Add(hCorrLead);
                }
                
                if(fMontecarlo){
                    THnSparseI *hPhiCorrMCc = new THnSparseI(namePlotThnMCc.Data(), "MCGen_SE_CorrCharm_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hPhiCorrMCc->Sumw2();
                    fOutputCorr->Add(hPhiCorrMCc);
                    
                    THnSparseI *hPhiCorrMCb = new THnSparseI(namePlotThnMCb.Data(), "MCGen_SE_CorrBeauty_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hPhiCorrMCb->Sumw2();
                    fOutputCorr->Add(hPhiCorrMCb);
                    
                    THnSparseI *hPhiCorrMCB = new THnSparseI(namePlotThnMCB.Data(), "MCGen_SE_CorrBMeson_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hPhiCorrMCB->Sumw2();
                    fOutputCorr->Add(hPhiCorrMCB);
                    
                    THnSparseI *hCorrLeadMC = new THnSparseI(namePlotThnLMC.Data(), "MCGen_SE_LeadingCorr_ThnSparse; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T; Evt Pool;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hCorrLeadMC->Sumw2();
                    fOutputCorr->Add(hCorrLeadMC);
                }
            }
            
            
            
            if(fMixing){  // Histogram for Event Mixing
                if(fReco){
                    namePlotThn+="_evMix";
                    THnSparseI *hPhiCorrMixing = new THnSparseI(namePlotThn.Data(), "Reco_ME_Corr_ThnSparse_EvMix; D^{+} Inv. Mass; #Delta#phi; #Delta#eta; track p_T;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hPhiCorrMixing->Sumw2();
                    fOutputCorr->Add(hPhiCorrMixing);
                    
                    // Correlations 2D for Leading Particles
                    namePlotThnL+="_evMix";
                    THnSparseI *hCorrLead = new THnSparseI(namePlotThnL.Data(), "Reco_ME_LeadingCorr_ThnSparse_EvMix; D^{+} Inv Mass; #Delta#phi; #Delta#eta; track p_T;",5,nVarsBins,nMinimumEdgeBins,nMaximumEdgeBins);
                    hCorrLead->Sumw2();
                    fOutputCorr->Add(hCorrLead);
                }
            }
            
    
    }
}

//____________________| Terminate
void AliAnalysisTaskSEDplusCorrelations::Terminate(Option_t */*option*/) {
    
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
