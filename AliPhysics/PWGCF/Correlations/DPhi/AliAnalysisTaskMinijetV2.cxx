#include <TChain.h>
#include <AliDirList.h>

#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TCanvas.h>
#include "TRandom.h"

#include "AliVEvent.h"
#include "AliVParticle.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliCentrality.h"

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

#include "AliAnalysisTaskMinijetV2.h"

// Analysis task for two-particle correlations using all particles over pt threshold
// pt_trig threshold for trigger particle (event axis) and pt_assoc for possible associated particles.
// Extract mini-jet yield and fragmentation properties via Delta-Phi histograms of these correlations
// post processing of analysis output via macro plot3and2Gaus.C
// Can use ESD or AOD, reconstructed and Monte Carlo data as input
// Author: Eva Sicking, modifications by Emilia Leogrande


ClassImp(AliAnalysisTaskMinijetV2)

//________________________________________________________________________
AliAnalysisTaskMinijetV2::AliAnalysisTaskMinijetV2(const char *name)
    : AliAnalysisTaskSE(name),
    fUseMC(kFALSE),
    fMcOnly(kFALSE),
    fBSign(0),
    fAnalysePrimOnly(kFALSE),// not used
    fPtMin(0.2),
    fPtMax(50.0),
    fCuts(0),
    fTriggerPtCut(0.7),
    fAssociatePtCut(0.4),
    fMode(0),
    fTriggerType(1),
    fFilterBit(128),
    fVertexZCut(10.),
    fEtaCut(0.9),
    fEtaCutSeed(0.9),
    fSelectParticles(1),
    fSelectParticlesAssoc(1),
    fCheckSDD(true),
    fSelOption(1),
    fCorrStrangeness(true),
    fThreeParticleCorr(false),
    fRejectChunks(false),
    fNTPC(10),
    fESDEvent(0),
    fAODEvent(0),
    fNMcPrimAccept(0),
    fNRecAccept(0),
    fNRecAcceptStrangeCorr(0),
    fNMcPrimAcceptTracklet(0),
    fNRecAcceptTracklet(0),
    fVzEvent(0),
    fMeanPtRec(0),
    fLeadingPtRec(0),
    fHists(0),
    fStep(0),
    fEventStat(0),
    fHistPt(0),
    fHistPtMC(0),
    fNContrNtracklets(0),
    fNContrNtracks(0),
    fCorruptedChunks(0),
    fCorruptedChunksAfter(0),
    fNmcNch(0),
    fPNmcNch(0),
    fNmcNchVtx(0),
    fNmcNchVtxStrangeCorr(0),
    fPNmcNchVtx(0),
    fNmcNchTracklet(0),
    fPNmcNchTracklet(0),
    fNmcNchVtxTracklet(0),
    fPNmcNchVtxTracklet(0),
    fChargedPi0(0),
    fVertexCheck(0),
    fPropagateDca(0),
    fAnalysisUtils(0),
    fCentralityMethod("")
{
    
    //Constructor
    
    for(Int_t i = 0;i< 8;i++){
        fMapSingleTrig[i]         =  0;
        fMapPair[i]               =  0;
        fMapEvent[i]              =  0;
        fMapAll[i]                =  0;
        fMapThree[i]              =  0;
        
        fVertexZ[i]               =  0;
        
        fNcharge[i]               =  0;
        fPt[i]                    =  0;
        fEta[i]                   =  0;
        fPhi[i]                   =  0;
        fDcaXY[i]                 =  0;
        fDcaZ[i]                  =  0;
        
        fPtSeed[i]                =  0;
        fEtaSeed[i]               =  0;
        fPhiSeed[i]               =  0;
        
        fPtOthers[i]              =  0;
        fEtaOthers[i]             =  0;
        fPhiOthers[i]             =  0;
        fPtEtaOthers[i]           =  0;
        
        fPhiEta[i]                =  0;
        
        fDPhiDEtaEventAxis[i]     =  0;
        fDPhiDEtaEventAxisSeeds[i]=  0;
        fTriggerNch[i]            =  0;
        
        fTriggerNchSeeds[i]       =  0;
        fTriggerTracklet[i]       =  0;
        
        fNch07Nch[i]              =  0;
        fPNch07Nch[i]             =  0;
        
        fNch07Tracklet[i]         =  0;
        fNchTracklet[i]           =  0;
        fPNch07Tracklet[i]        =  0;
        fDPhiEventAxis[i]         =  0;
        fDPhi1DPhi2[i]            =  0;
    }
    DefineOutput(1, AliDirList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMinijetV2::~AliAnalysisTaskMinijetV2()
{
    // Destructor
    
    if (fHists && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHists;
}

//________________________________________________________________________
void AliAnalysisTaskMinijetV2::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    if(fDebug) Printf("In User Create Output Objects.");
    
    Int_t nbinsCentr = 0;
    Float_t minbinCentr=0, maxbinCentr=0;

    if (fCentralityMethod.Length() > 0)
    {
        nbinsCentr = 105;
        minbinCentr=0;
        maxbinCentr=105;
    }
    else
    {
        nbinsCentr = 101;
        minbinCentr=-0.5;
        maxbinCentr=100.5;
    }

    fStep = new TH1F("fStep", "fStep", 10, -0.5, 9.5);
    fEventStat = new TH1F("fEventStat", "fEventStat", 10, -0.5, 9.5);
    fHistPt = new TH1F("fHistPt", "P_{T} distribution REC", 150, 0.1, 3.1);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);
    fNContrNtracklets = new TH2F ("fNContrNtracklets", ";N_{tracklets};N_{vtx contrib}", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
    fNContrNtracks    = new TH2F ("fNContrNtracks", ";N_{tracks};N_{vtx contrib}", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
    fCorruptedChunks  = new TH2F ("fCorruptedChunks",
                                  ";N_{tracks,TPC};N_{tracks,ITS-TPC}", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
    fCorruptedChunksAfter  = new TH2F ("fCorruptedChunksAfter",
                                       ";N_{tracks,TPC};N_{tracks,ITS-TPC}", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
    

    
    
    if (fUseMC) {
        fHistPtMC = new TH1F("fHistPtMC", "P_{T} distribution MC", 150, 0.1, 3.1);
        fHistPtMC->GetXaxis()->SetTitle("P_{T} (GeV/c)");
        fHistPtMC->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
        fHistPtMC->SetMarkerStyle(kFullCircle);
        
        fNmcNch = new TH2F("fNmcNch", "fNmcNch", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
        fPNmcNch = new TProfile("pNmcNch", "pNmcNch", nbinsCentr,minbinCentr,maxbinCentr);
        fNmcNchVtx = new TH2F("fNmcNchVtx", "fNmcNchVtx", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
        fNmcNchVtxStrangeCorr = new TH2F("fNmcNchVtxStrangeCorr", "fNmcNchVtxStrangeCorr", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
        fPNmcNchVtx = new TProfile("pNmcNchVtx", "pNmcNchVtx", nbinsCentr,minbinCentr,maxbinCentr);
        
        fNmcNchTracklet = new TH2F("fNmcNchTracklet", "fNmcNchTracklet", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
        fPNmcNchTracklet = new TProfile("pNmcNchTracklet", "pNmcNchTracklet", nbinsCentr,minbinCentr,maxbinCentr);
        fNmcNchVtxTracklet = new TH2F("fNmcNchVtxTracklet", "fNmcNchVtxTracklet", nbinsCentr,minbinCentr,maxbinCentr,nbinsCentr,minbinCentr,maxbinCentr);
        fPNmcNchVtxTracklet = new TProfile("pNmcNchVtxTracklet", "pNmcNchVtxTracklet", nbinsCentr,minbinCentr,maxbinCentr);
        

        
    }
    
    fChargedPi0  = new TH2F("fChargedPi0", "fChargedPi0", 200, -0.5, 199.5, 200, -0.5, 199.5);
    fVertexCheck = new TH1F("fVertexCheck", "fVertexCheck", 100, -0.5, 99.5);
    fPropagateDca = new TH1F("fPropagateDca", "fPropagateDca", 2, -0.5, 1.5);
    
    //----------------------
    //bins for pt in THnSpare
    Double_t ptMin = 0.0, ptMax = 100.;
    Int_t nPtBins = 37;
    Double_t binsPt[]  = {0.0, 0.1, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8,
        0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
        10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
    
    //Int_t nPtBinsAll = 100;
    //Double_t ptMinAll = 1.e-2, ptMaxAll = 100.;
    //Double_t *binsPtAll = 0;
    //binsPtAll = (Double_t*)CreateLogAxis(nPtBinsAll,ptMinAll,ptMaxAll);
    
    
    //   Double_t ptMin2 = 0.0, ptMax2 = 100.;
    //   Int_t nPtBins2 = 9;
    //   Double_t binsPt2[]  = {0.1, 0.4, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 100.0};
    
    //  Int_t nPtBins2 = 10;
    //   Double_t ptMin2 = 0.4, ptMax2 = 100.;
    //   Double_t *binsPt2 = 0;
    //   binsPt2 = CreateLogAxis(nPtBins2,ptMin2,ptMax2);
    
    //3 dim matrix
    Int_t binsEffHisto[3]   = {Int_t(fEtaCut*20),  nPtBins,       nbinsCentr };
    Double_t minEffHisto[3] = {-fEtaCut,           ptMin,        minbinCentr };
    Double_t maxEffHisto[3] = {fEtaCut,            ptMax,        maxbinCentr};
    
    //5 dim matrix
    Int_t binsEffHisto5[6]   = {  nPtBins,   nPtBins,    1,                              90,         nbinsCentr ,      2 };
    Double_t minEffHisto5[6] = {  ptMin,     ptMin,     -2*fEtaCut,          -0.5*TMath::Pi(),       minbinCentr ,   -0.5 };
    Double_t maxEffHisto5[6] = {  ptMax,     ptMax,      2*fEtaCut,           1.5*TMath::Pi(),      maxbinCentr ,    1.5 };
    
    
    //4 dim matrix
    Int_t binsEvent[4]   = {   nbinsCentr,        20,   50,  nPtBins };
    Double_t minEvent[4] = {   minbinCentr,       -10,    0,    ptMin };
    Double_t maxEvent[4] = {   maxbinCentr,        10,   10,    ptMax };
    
    //3 dim matrix
    Int_t binsAll[3]   = {Int_t(fEtaCut*20),  nPtBins,       nbinsCentr };
    Double_t minAll[3] = {-fEtaCut,           ptMin,         minbinCentr };
    Double_t maxAll[3] = {fEtaCut,            ptMax,        maxbinCentr };
    
    //3 dim matrix
    Int_t binsThree[3]   = {              90,               90,    nbinsCentr};
    Double_t minThree[3] = {-0.5*TMath::Pi(), -0.5*TMath::Pi(),    minbinCentr};
    Double_t maxThree[3] = { 1.5*TMath::Pi(),  1.5*TMath::Pi(),    maxbinCentr};
    

    
    
    
    //--------------------
    TString dataType[2] ={"ESD", "AOD"};
    TString labels[8]={Form("%sAllAllMcNmc",dataType[fMode].Data()),
        Form("%sTrigAllMcNmc",dataType[fMode].Data()),
        Form("%sTrigVtxMcNmc",dataType[fMode].Data()),
        Form("%sTrigVtxMcNrec",dataType[fMode].Data()),
        Form("%sTrigVtxRecMcPropNrec",dataType[fMode].Data()),
        Form("%sTrigVtxRecNrec",dataType[fMode].Data()),
        Form("%sTrigVtxRecMcPropNrecStrangeCorr",dataType[fMode].Data()),
        Form("%sTrigVtxRecNrecStrangeCorr",dataType[fMode].Data())};
    
    
    for(Int_t i=0;i<8;i++){
        
        fMapSingleTrig[i] = new THnSparseD(Form("fMapSingleTrig%s", labels[i].Data()),"eta:pt:Nrec",
                                           3,binsEffHisto,minEffHisto,maxEffHisto);
        fMapSingleTrig[i]->SetBinEdges(1,binsPt);
        fMapSingleTrig[i]->GetAxis(0)->SetTitle("#eta");
        fMapSingleTrig[i]->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
        fMapSingleTrig[i]->GetAxis(2)->SetTitle("N_{rec}");
        fMapSingleTrig[i]->Sumw2();
        
        fMapPair[i] = new THnSparseD(Form("fMapPair%s", labels[i].Data()),"pt_trig:pt_assoc:DeltaEta:DeltaPhi:Nrec:likesign",
                                     6,binsEffHisto5,minEffHisto5,maxEffHisto5);
        fMapPair[i]->SetBinEdges(0,binsPt);
        fMapPair[i]->SetBinEdges(1,binsPt);
        fMapPair[i]->GetAxis(0)->SetTitle("p_{T, trig} (GeV/c)");
        fMapPair[i]->GetAxis(1)->SetTitle("p_{T, assoc} (GeV/c)");
        fMapPair[i]->GetAxis(2)->SetTitle("#Delta #eta");
        fMapPair[i]->GetAxis(3)->SetTitle("#Delta #phi");
        fMapPair[i]->GetAxis(4)->SetTitle("N_{rec}");
        fMapPair[i]->GetAxis(5)->SetTitle("Like-sign or Unlike-sign");
        fMapPair[i]->Sumw2();
        
        
        fMapEvent[i] = new THnSparseD(Form("fMapEvent%s", labels[i].Data()),"Nrec:vertexZ:meanPt:leadingPt",
                                      4,binsEvent,minEvent,maxEvent);
        fMapEvent[i]->GetAxis(0)->SetTitle("N_{rec}");
        fMapEvent[i]->GetAxis(1)->SetTitle("z_{vertex} (cm)");
        fMapEvent[i]->GetAxis(2)->SetTitle("meanPt");
        fMapEvent[i]->SetBinEdges(3,binsPt);
        fMapEvent[i]->GetAxis(3)->SetTitle("leadingPt");
        fMapEvent[i]->Sumw2();
        
        fMapAll[i] = new THnSparseD(Form("fMapAll%s", labels[i].Data()),"eta:pt:Nrec",
                                    3,binsAll,minAll,maxAll);
        fMapAll[i]->SetBinEdges(1,binsPt);
        fMapAll[i]->GetAxis(0)->SetTitle("#eta");
        fMapAll[i]->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
        fMapAll[i]->GetAxis(2)->SetTitle("N_{rec}");
        fMapAll[i]->Sumw2();
        
        
        fMapThree[i] = new THnSparseD(Form("fMapThree%s", labels[i].Data()),"dphi1:dphi2:Nrec",
                                      3,binsThree,minThree,maxThree);
        fMapThree[i]->GetAxis(0)->SetTitle("#Delta#varphi_{1}");
        fMapThree[i]->GetAxis(1)->SetTitle("#Delta#varphi_{2}");
        fMapThree[i]->GetAxis(2)->SetTitle("N_{rec}");
        fMapThree[i]->Sumw2();
        
        
        fVertexZ[i]                  = new TH1F(Form("fVertexZ%s",labels[i].Data()),
                                                Form("fVertexZ%s",labels[i].Data()) ,
                                                220, -11., 11.);
        fPt[i]                       = new TH1F(Form("fPt%s",labels[i].Data()),
                                                Form("fPt%s",labels[i].Data()) ,
                                                200, 0., 50);
        fNcharge[i]                  = new TH1F(Form("fNcharge%s",labels[i].Data()),
                                                Form("fNcharge%s",labels[i].Data()) ,
                                                250, -0.5, 249.5);
        fEta[i]                      = new TH1F (Form("fEta%s",labels[i].Data()),
                                                 Form("fEta%s",labels[i].Data()) ,
                                                 100, -2., 2);
        fPhi[i]                      = new TH1F(Form("fPhi%s",labels[i].Data()),
                                                Form("fPhi%s",labels[i].Data()) ,
                                                360, 0.,2*TMath::Pi());
        fDcaXY[i]                    = new TH1F(Form("fDcaXY%s",labels[i].Data()),
                                                Form("fDcaXY%s",labels[i].Data()) ,
                                                200, -0.3,0.3);
        fDcaZ[i]                     = new TH1F(Form("fDcaZ%s",labels[i].Data()),
                                                Form("fDcaZ%s",labels[i].Data()) ,
                                                200, -2.2,2.2);
        fPtSeed[i]                       = new TH1F(Form("fPSeedt%s",labels[i].Data()),
                                                    Form("fPtSeed%s",labels[i].Data()) ,
                                                    500, 0., 50);
        fEtaSeed[i]                      = new TH1F (Form("fEtaSeed%s",labels[i].Data()),
                                                     Form("fEtaSeed%s",labels[i].Data()) ,
                                                     100, -2., 2);
        fPhiSeed[i]                      = new TH1F(Form("fPhiSeed%s",labels[i].Data()),
                                                    Form("fPhiSeed%s",labels[i].Data()) ,
                                                    360, 0.,2*TMath::Pi());
        
        fPtOthers[i]                       = new TH1F(Form("fPOtherst%s",labels[i].Data()),
                                                      Form("fPtOthers%s",labels[i].Data()) ,
                                                      500, 0., 50);
        fEtaOthers[i]                      = new TH1F (Form("fEtaOthers%s",labels[i].Data()),
                                                       Form("fEtaOthers%s",labels[i].Data()) ,
                                                       100, -2., 2);
        fPhiOthers[i]                      = new TH1F(Form("fPhiOthers%s",labels[i].Data()),
                                                      Form("fPhiOthers%s",labels[i].Data()) ,
                                                      360, 0.,2*TMath::Pi());
        fPtEtaOthers[i]                      = new TH2F(Form("fPtEtaOthers%s",labels[i].Data()),
                                                        Form("fPtEtaOthers%s",labels[i].Data()) ,
                                                        500, 0., 50, 100, -1., 1);
        fPhiEta[i]                   = new TH2F(Form("fPhiEta%s",labels[i].Data()),
                                                Form("fPhiEta%s",labels[i].Data()) ,
                                                180, 0., 2*TMath::Pi(), 100, -1.,1.);
        fDPhiDEtaEventAxis[i]        = new TH2F(Form("fDPhiDEtaEventAxis%s",labels[i].Data()),
                                                Form("fDPhiDEtaEventAxis%s",labels[i].Data()) ,
                                                180, -0.5* TMath::Pi(), 1.5*TMath::Pi(), 9, -1.8,1.8);
        fTriggerNch[i]               = new TH1F(Form("fTriggerNch%s",labels[i].Data()),
                                                Form("fTriggerNch%s",labels[i].Data()) ,
                                                250, -0.5, 249.5);
        fTriggerNchSeeds[i]          = new TH2F(Form("fTriggerNchSeeds%s",labels[i].Data()),
                                                Form("fTriggerNchSeeds%s",labels[i].Data()) ,
                                                250, -0.5, 249.5, 100, -0.5, 99.5);
        fTriggerTracklet[i]          = new TH1F(Form("fTriggerTracklet%s",labels[i].Data()),
                                                Form("fTriggerTracklet%s",labels[i].Data()) ,
                                                250, -0.5, 249.5);
        fNch07Nch[i]                 = new TH2F(Form("fNch07Nch%s",labels[i].Data()),
                                                Form("fNch07Nch%s",labels[i].Data()) ,
                                                250, -2.5, 247.5,250, -2.5, 247.5);
        fPNch07Nch[i]                 = new TProfile(Form("pNch07Nch%s",labels[i].Data()),
                                                     Form("pNch07Nch%s",labels[i].Data()) ,
                                                     250, -2.5, 247.5);
        fNch07Tracklet[i]            = new TH2F(Form("fNch07Tracklet%s",labels[i].Data()),
                                                Form("fNch07Tracklet%s",labels[i].Data()) ,
                                                250, -2.5, 247.5,250, -2.5, 247.5);
        fNchTracklet[i]              = new TH2F(Form("fNchTracklet%s",labels[i].Data()),
                                                Form("fNchTracklet%s",labels[i].Data()) ,  
                                                250, -2.5, 247.5,250, -2.5, 247.5);
        fPNch07Tracklet[i]            = new TProfile(Form("pNch07Tracklet%s",labels[i].Data()),
                                                     Form("pNch07Tracklet%s",labels[i].Data()) ,  
                                                     250, -2.5, 247.5);
        fDPhiEventAxis[i]          = new TH1F(Form("fDPhiEventAxis%s",
                                                   labels[i].Data()),
                                              Form("fDPhiEventAxis%s",
                                                   labels[i].Data()) ,  
                                              180, -0.5*TMath::Pi(), 1.5*TMath::Pi());
        fDPhi1DPhi2[i] = new TH2F(Form("fDPhi1DPhi2%s",labels[i].Data()),
                                  Form("fDPhi1DPhi2%s",labels[i].Data()),
                                  180, -0.5* TMath::Pi(), 1.5*TMath::Pi(), 
                                  180, -0.5* TMath::Pi(), 1.5*TMath::Pi());   
    }
    
    fHists = new AliDirList();
    fHists->SetOwner();
    
    fHists->Add(fStep);
    fHists->Add(fEventStat);
    fHists->Add(fHistPt);
    fHists->Add(fNContrNtracklets);
    fHists->Add(fNContrNtracks);
    fHists->Add(fCorruptedChunks);
    fHists->Add(fCorruptedChunksAfter);
    
    if(fUseMC){
        fHists->Add(fHistPtMC); 
        
        fHists->Add(fNmcNch); 
        fHists->Add(fPNmcNch); 
        fHists->Add(fNmcNchVtx); 
        fHists->Add(fNmcNchVtxStrangeCorr); 
        fHists->Add(fPNmcNchVtx); 
        
        fHists->Add(fNmcNchTracklet); 
        fHists->Add(fPNmcNchTracklet); 
        fHists->Add(fNmcNchVtxTracklet); 
        fHists->Add(fPNmcNchVtxTracklet); 
    }
    fHists->Add(fChargedPi0);
    fHists->Add(fVertexCheck);
    fHists->Add(fPropagateDca);
    
    for(Int_t i=0;i<8;i++){
        fHists->Add(fMapSingleTrig[i]);
        fHists->Add(fMapPair[i]);
        fHists->Add(fMapEvent[i]);
        fHists->Add(fMapAll[i]);
        fHists->Add(fMapThree[i]);
        fHists->Add(fVertexZ[i]);
        fHists->Add(fPt[i]);
        fHists->Add(fNcharge[i]);
        fHists->Add(fEta[i]);
        fHists->Add(fPhi[i]);
        fHists->Add(fDcaXY[i]);
        fHists->Add(fDcaZ[i]);
        fHists->Add(fPtSeed[i]);
        fHists->Add(fEtaSeed[i]);
        fHists->Add(fPhiSeed[i]);
        fHists->Add(fPtOthers[i]);
        fHists->Add(fEtaOthers[i]);
        fHists->Add(fPhiOthers[i]);
        fHists->Add(fPtEtaOthers[i]);
        fHists->Add(fPhiEta[i]);
        fHists->Add(fDPhiDEtaEventAxis[i]);
        fHists->Add(fTriggerNch[i]);
        fHists->Add(fTriggerNchSeeds[i]);
        fHists->Add(fTriggerTracklet[i]);
        fHists->Add(fNch07Nch[i]);
        fHists->Add(fPNch07Nch[i]);
        fHists->Add(fNch07Tracklet[i]);
        fHists->Add(fNchTracklet[i]);
        fHists->Add(fPNch07Tracklet[i]);
        fHists->Add(fDPhiEventAxis[i]);
        fHists->Add(fDPhi1DPhi2[i]);
    }
    
    PostData(1, fHists);
    
}

//________________________________________________________________________
void AliAnalysisTaskMinijetV2::UserExec(Option_t *)
{
    // Main function, called for each event
    // Kinematics-only, ESD and AOD can be processed.
    // Data is read (ReadEventESD, ReadEventAOD...) and then analysed (Analyse).
    //  - in case of MC with full detector simulation, all correction steps(0-5) can be processed
    //  - for Data, only step 5 is performed
    //  - for kinematics-only, only step 0 is processed
    
    //             - trigger -               - vertex -                       - tracks -                                       - multiplicity -
    // step 7 =  Triggered events, reconstructed accepted vertex,  reconstructed tracks + strangness corr,               reconstructed multiplicity+strangeCorr
    // step 6 =  Triggered events, reconstructed accepted vertex,  reconstructed tracks with MC prop + stangess corr,    reconstructed multiplicity+strangeCorr
    // step 5 =  Triggered events, reconstructed accepted vertex,  reconstructed tracks,                                 reconstructed multiplicity
    // step 4 =  Triggered events, reconstructed accepted vertex,  reconstructed tracks with MC properties,              reconstructed multiplicity
    // step 3 =  Triggered events, reconstructed accepted vertex,  mc primary particles,                                 reconstructed multiplicity
    // step 2 =  Triggered events, reconstructed accepted vertex,  mc primary particles,                                 true multiplicity
    // step 1 =  Triggered events, all                             mc primary particles,                                 true multiplicity
    // step 0 =  All events,       all                             mc primary particles,                                 true multiplicity
    
    
    if(fDebug) Printf("UserExec: Event starts");
    
    AliVEvent *event = InputEvent();
    if (!event) {
        Error("UserExec", "Could not retrieve event");
        return;
    }
    fBSign= event->GetMagneticField();
    
    //get events, either ESD or AOD
    fESDEvent = dynamic_cast<AliESDEvent*> (InputEvent());
    fAODEvent = dynamic_cast<AliAODEvent*> (InputEvent());
    
    vector<Float_t> pt;
    vector<Float_t> eta;
    vector<Float_t> phi;
    vector<Short_t> charge;
    vector<Float_t> strangenessWeight;
    vector<Float_t> px;
    vector<Float_t> py;
    vector<Float_t> pz;
    vector<Float_t> theta;
    
    
    //number of accepted tracks and tracklets
    Double_t ntracks = 0;  //return value for reading functions for ESD and AOD
    //Int_t ntracksRemove = 0;  //return value for reading functions for ESD and AOD
    vector<Double_t> nTracksTracklets; // [0]=nAccepted,1=ntracklets,2=nall(also neutral in case of mc,
    //for real nall=ncharged)
    
    if(!fAODEvent && !fESDEvent)return;
    
    //Centrality
    Double_t centrality = 0;
    if (fCentralityMethod.Length() > 0)
    {
        AliCentrality *centralityObj = 0;
        if (fAODEvent)
            centralityObj = ((AliVAODHeader*)fAODEvent->GetHeader())->GetCentralityP();
        else if (fESDEvent)
            centralityObj = fESDEvent->GetCentrality();
        if (centralityObj)
        {
            centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
            AliInfo(Form("Centrality is %f", centrality));
        }
        else
        {
            Printf("WARNING: Centrality object is 0");
            centrality = -1;
        }
    }
    
    // SDD check for LHC11a
    if (fCheckSDD) {
        
        //ESD
        if(fESDEvent){
            if(fSelOption==0){
                const AliMultiplicity *mul = fESDEvent->GetMultiplicity();
                Int_t nClu3 = mul->GetNumberOfITSClusters(2);
                Int_t nClu4 = mul->GetNumberOfITSClusters(3);
                if(nClu3==0 &&  nClu4==0) return;
            }
            else if (fSelOption==1){
                TString trcl = fESDEvent->GetFiredTriggerClasses().Data();
                if (!(trcl.Contains("CINT1-B-NOPF-ALLNOTRD"))) return;
            }
        }
        
        //AOD
        if(fAODEvent){
            if(fSelOption==0){
                Bool_t useEvent = false;
                Int_t nTracks = fAODEvent->GetNumberOfTracks();
                for(Int_t itrack=0; itrack<nTracks; itrack++) {
                    AliAODTrack * track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(itrack));
                    if(!track) AliFatal("Not a standard AOD");
                    if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
                        useEvent=true;
                        break;
                    }
                }
                if (!useEvent) return;
            }
            else if(fSelOption==1){
                TString trcl = fAODEvent->GetFiredTriggerClasses().Data();
                if (!(trcl.Contains("CINT1-B-NOPF-ALLNOTRD"))) return;
            }
        }
    }
    
    //reset values
    fNMcPrimAccept=0;// number of accepted primaries
    fNRecAccept=0;   // number of accepted tracks
    fNRecAcceptStrangeCorr=0;   // number of accepted tracks + strangeness correction
    fNMcPrimAcceptTracklet=0;// number of accepted primaries (no pt cut)
    fNRecAcceptTracklet=0;   // number of accepted tracklets
    
    // instead of task->SelectCollisionCandidate(mask) in AddTask macro
    Bool_t isSelected = (((AliInputEventHandler*)
                          (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
                         ->IsEventSelected() &  fTriggerType);
    
    
    if(fDebug){
        Printf("IsSelected = %d", isSelected);
        Printf("CheckEvent(true)= %d", CheckEvent(true));
        Printf("CheckEvent(false)= %d", CheckEvent(false));
    }
    
    fEventStat->Fill(0);//all events
    
    //check trigger
    if(isSelected){ // has offline trigger
        
        fEventStat->Fill(1);//triggered event
        
        if(CheckEvent(true)){//step 5 = TrigVtxRecNrec, step 4 = TrigVtxRecMcPropNrec ,step 3 = TrigVtxMcNrec
            
            fEventStat->Fill(2);//triggered event with vertex
            
            if(!fMcOnly){
                //step 5 = TrigVtxRecNrec
                
                // read tracks
                if(fESDEvent)     ntracks = ReadEventESD(pt, eta, phi, charge,strangenessWeight, nTracksTracklets, 5);
                else if(fAODEvent)ntracks = ReadEventAOD(pt, eta, phi, charge,strangenessWeight, nTracksTracklets, 5);
                else AliInfo("Fatal Error");
                
                if (fCentralityMethod.Length() > 0)
                    ntracks = centrality;
                
                // analyse
                if(pt.size()){ //(internally ntracks=fNRecAccept)
                    fEventStat->Fill(3);//triggered event with vertex and one reconstructed track in acceptance
                        Analyse(pt, eta, phi, charge, strangenessWeight, ntracks, nTracksTracklets[1], nTracksTracklets[2], 5);//step 5 = TrigVtxRecNrec
                    }
                
                if(fCorrStrangeness){
                    //step 7 = TrigVtxRecNrecStrangeCorr
                    
                    // read tracks
                    if(fESDEvent)     ntracks = ReadEventESD(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 7);//stagness version not yet implemented
                    else if(fAODEvent)ntracks = ReadEventAOD(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 7);
                    else AliInfo("Fatal Error");
                    
                    if (fCentralityMethod.Length() > 0)
                        ntracks = centrality;
                    
                    // analyse
                    if(pt.size()){ //(internally ntracks=fNRecAccept)
                        Analyse(pt, eta, phi, charge, strangenessWeight, fNRecAcceptStrangeCorr, nTracksTracklets[1], nTracksTracklets[2], 7);//step 7 = TrigVtxRecNrecStrangeCorr
                    }
                }
                
                if(fUseMC){
                    // step 4 = TrigVtxRecMcPropNrec
                    
                    // read tracks
                    if(fESDEvent)       ntracks = ReadEventESDRecMcProp(pt, eta, phi, charge,strangenessWeight, nTracksTracklets, 4);
                    else if(fAODEvent)  ntracks = ReadEventAODRecMcProp(pt, eta, phi, charge,strangenessWeight, nTracksTracklets, 4);
                    else AliInfo("Fatal Error");
                    
                    if (fCentralityMethod.Length() > 0)
                        ntracks = centrality;
                    
                    //analyse
                    if(pt.size()){//(internally ntracks=fNRecAccept)
                        Analyse(pt, eta, phi, charge,strangenessWeight, ntracks, nTracksTracklets[1], nTracksTracklets[2], 4); //step 4 = TrigVtxRecMcPropNrec
                    }
                    
                    
                    if(fCorrStrangeness){
                        // step 6 = TrigVtxRecMcPropNrecStrangeCorr
                        
                        // read tracks
                        if(fESDEvent)       ntracks = ReadEventESDRecMcProp(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 6);//stagness version not yet implemented
                        else if(fAODEvent)  ntracks = ReadEventAODRecMcProp(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 6);
                        else AliInfo("Fatal Error");
                        
                        if (fCentralityMethod.Length() > 0)
                            ntracks = centrality;
                        
                        //analyse
                        if(pt.size()){//(internally ntracks=fNRecAccept)
                            Analyse(pt, eta, phi, charge, strangenessWeight, fNRecAcceptStrangeCorr, nTracksTracklets[1], nTracksTracklets[2], 6); //step 6 = TrigVtxRecMcPropNrecStrangeCorr
                        }
                    }
                    // step 3 = TrigVtxMcNrec
                    
                    // read tracks
                    if(fESDEvent)       ntracks = ReadEventESDMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 3);
                    else if(fAODEvent)  ntracks = ReadEventAODMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 3);
                    else AliInfo("Fatal Error");
                    
                    if (fCentralityMethod.Length() > 0){
                        fNRecAccept = centrality;
                        fNMcPrimAccept = centrality;
                    }
                    
                    // analyse
                    if(pt.size()){
                        Analyse(pt, eta, phi, charge, strangenessWeight, fNRecAccept,    nTracksTracklets[1],nTracksTracklets[2], 3); //step 3 = TrigVtxMcNrec
                        Analyse(pt, eta, phi, charge, strangenessWeight, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 2); //step 2 = TrigVtxMcNmc
                    }
                    
                }
            }
            
        }//check event (true)
        
        
        if(fUseMC && !fMcOnly){
            //reset values
            fNMcPrimAccept=0;// number of accepted primaries
            fNRecAccept=0;   // number of accepted tracks
            fNMcPrimAcceptTracklet=0;// number of accepted primaries (no pt cut)
            fNRecAcceptTracklet=0;   // number of accepted tracklets
            
            if(CheckEvent(false)){// all events, with and without reconstucted vertex
                if(fESDEvent) ntracks       = ReadEventESDMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 1);//read tracks
                else if(fAODEvent) ntracks  = ReadEventAODMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 1);//read tracks
                else AliInfo("Fatal Error");
                
                if (fCentralityMethod.Length() > 0)
                    fNMcPrimAccept = centrality;
                
                
                // analyse
                if(pt.size()){
                    Analyse(pt, eta, phi, charge, strangenessWeight, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 1);  // step 1 = TrigAllMcNmc
                    
                    Analyse(pt, eta, phi, charge, strangenessWeight, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 0);  //first part of step 0 // step 0 = AllAllMcNmc
                }
                
                
            }
        }
    }// triggered event
    
    else { // not selected by physics selection task = not triggered
        if(fUseMC && !fMcOnly){
            if(CheckEvent(false)){
                
                //read tracks
                if(fESDEvent)	   ntracks  = ReadEventESDMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 0);
                else if(fAODEvent) ntracks  = ReadEventAODMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 0);
                else AliInfo("Fatal Error");
                
                if (fCentralityMethod.Length() > 0)
                    fNMcPrimAccept = centrality;
                
                //analyse
                if(pt.size()){
                    Analyse(pt, eta, phi, charge, strangenessWeight, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 0);  //second part of step 0 // step 0 = AllAllMcNmc
                }
            }
        }
    }
    
    if(fMcOnly){
        // read event
        if(fMode==0)       ntracks  = ReadEventESDMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 0);
        else if (fMode==1) ntracks  = ReadEventAODMC(pt, eta, phi, charge, strangenessWeight, nTracksTracklets, 0);
        
        if (fCentralityMethod.Length() > 0)
            fNMcPrimAccept = centrality;
        
        // analyse
        if(pt.size()){
            Analyse(pt, eta, phi, charge, strangenessWeight, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 0); 
        }
    }
}      


//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventESD( vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                           vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
                                           vector<Float_t> &strangeArray,
                                           vector<Double_t>   &nTracksTracklets, const Int_t step)
{
    // gives back the number of esd tracks and pointer to arrays with track
    // properties (pt, eta, phi)
    // Uses TPC tracks with SPD vertex from now on
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    const AliESDVertex*	vtxSPD   = fESDEvent->GetPrimaryVertexSPD(); // uses track or SPD vertexer
    fVertexZ[step]->Fill(vtxSPD->GetZ());
    
    // Retreive the number of all tracks for this event.
    Int_t ntracks = fESDEvent->GetNumberOfTracks();
    if(fDebug>1)  Printf("all ESD tracks: %d", ntracks);
    
    //first loop to check how many tracks are accepted
    //------------------
    Double_t nAcceptedTracks=0;
    //Float_t nAcceptedTracksStrange=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        
        AliESDtrack *esdTrack = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
        if (!esdTrack) {
            Error("ReadEventESD", "Could not receive track %d", iTracks);
            continue;
        }
        
        // use TPC only tracks with non default SPD vertex
        AliESDtrack *track = AliESDtrackCuts::
        GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(fESDEvent),esdTrack->GetID());
        if(!track) continue;
        if(!fCuts->AcceptTrack(track)) {
            delete track;
            continue;// apply TPC track cuts before constrain to SPD vertex
        }
        if(track->Pt()>0.){
            // only constrain tracks above threshold
            AliExternalTrackParam exParam;
            // take the B-field from the ESD, no 3D fieldMap available at this point
            Bool_t relate = false;
            relate = track->RelateToVertexTPC(vtxSPD,fESDEvent->GetMagneticField(),
                                              kVeryBig,&exParam);
            if(!relate){
                delete track;
                continue;
            }
            track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),
                       exParam.GetCovariance());
        }
        else{
            delete track;
            continue;// only if tracks have pt<=0
        }
        
        if (TMath::Abs(track->Eta())<fEtaCut && track->Pt()>fPtMin && track->Pt()<fPtMax) {
            ptArray.push_back(track->Pt());
            etaArray.push_back(track->Eta());
            phiArray.push_back(track->Phi());
            chargeArray.push_back(track->Charge());
            strangeArray.push_back(1);
            ++nAcceptedTracks;
            fHistPt->Fill(track->Pt());
        }
        
        // TPC only track memory needs to be cleaned up
        if(track)
            delete track;
        
    }
    
    //need to be checked
    if(nAcceptedTracks==0) return 0;
    
    //tracklet loop
    Int_t ntrackletsAccept=0;
    AliMultiplicity * mult = (AliMultiplicity*)(fESDEvent->GetMultiplicity());
    Int_t ntracklets = mult->GetNumberOfTracklets();
    for(Int_t i=0;i< ntracklets;i++){
        if(mult->GetDeltaPhi(i)<0.05 && TMath::Abs(mult->GetEta(i))<fEtaCut){
            ntrackletsAccept++;
        }
    }
    nTracksTracklets.push_back(nAcceptedTracks);
    nTracksTracklets.push_back(ntrackletsAccept);
    nTracksTracklets.push_back(nAcceptedTracks);//in order to be compatible with mc analysis
    //where here also neutral particles are counted.
    
    
    fVzEvent=vtxSPD->GetZ(); // needed for correction map
    if(step==5 || step ==2) {
        fNRecAccept=nAcceptedTracks;
        fNRecAcceptTracklet=ntrackletsAccept;
    }
    return fNRecAccept;
    
    
}

//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventESDRecMcProp( vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                                    vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
                                                    vector<Float_t> &strangeArray,
                                                    vector<Double_t> &nTracksTracklets, const Int_t step)
{
    // gives back the number of esd tracks and pointer to arrays with track
    // properties (pt, eta, phi) of mc particles if available
    // Uses TPC tracks with SPD vertex from now on
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    
    AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
    if (!mcEvent) {
        Error("ReadEventESDRecMcProp", "Could not retrieve MC event");
        return 0;
    }
    AliStack* stack = MCEvent()->Stack();
    if(!stack) return 0;
    
    
    // Retreive the number of all tracks for this event.
    Int_t ntracks = fESDEvent->GetNumberOfTracks();
    if(fDebug>1)Printf("all ESD tracks: %d", ntracks);
    
    const AliESDVertex *vtxSPD = fESDEvent->GetPrimaryVertexSPD();
    fVertexZ[step]->Fill(vtxSPD->GetZ());
    
    //track loop
    Double_t nAcceptedTracks=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        
        AliVParticle *vtrack = fESDEvent->GetTrack(iTracks);
        AliESDtrack *esdTrack = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
        if (!esdTrack) {
            Error("ReadEventESDRecMcProp", "Could not receive track %d", iTracks);
            continue;
        }
        
        // use TPC only tracks with non default SPD vertex
        AliESDtrack *track = AliESDtrackCuts::
        GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(fESDEvent),esdTrack->GetID());
        if(!track) continue;
        if(!fCuts->AcceptTrack(track)) {
            delete track;
            continue;// apply TPC track cuts before constrain to SPD vertex
        }
        if(track->Pt()>0.){
            // only constrain tracks above threshold
            AliExternalTrackParam exParam;
            // take the B-field from the ESD, no 3D fieldMap available at this point
            Bool_t relate = false;
            relate = track->RelateToVertexTPC(vtxSPD,fESDEvent->GetMagneticField(),
                                              kVeryBig,&exParam);
            if(!relate){
                delete track;
                continue;
            }
            track->Set(exParam.GetX(),exParam.GetAlpha(),exParam.GetParameter(),
                       exParam.GetCovariance());
        }
        else{
            delete track;
            continue;// only if tracks have pt<=0
        }
        
        //count tracks, if available, use mc particle properties
        if(vtrack->GetLabel()<=0){
            if (TMath::Abs(track->Eta())<fEtaCut && track->Pt()>fPtMin && track->Pt()<fPtMax){
                ptArray.push_back(track->Pt());
                etaArray.push_back(track->Eta());
                phiArray.push_back(track->Phi());
                chargeArray.push_back(track->Charge());
                strangeArray.push_back(1);
                ++nAcceptedTracks;
            }
        }
        else{
            TParticle *partOfTrack = stack->Particle(vtrack->GetLabel());
            if (TMath::Abs(partOfTrack->Eta())<fEtaCut && partOfTrack->Pt()>fPtMin && partOfTrack->Pt()<fPtMax) {
                ptArray.push_back(partOfTrack->Pt());
                etaArray.push_back(partOfTrack->Eta());
                phiArray.push_back(partOfTrack->Phi());
                chargeArray.push_back(vtrack->Charge());
                strangeArray.push_back(1);
                ++nAcceptedTracks;
            }
        }
        
        // TPC only track memory needs to be cleaned up
        if(track)
            delete track;
        
    }
    
    if(nAcceptedTracks==0) return 0;
    
    //tracklet loop
    Int_t ntrackletsAccept=0;
    AliMultiplicity * mult = (AliMultiplicity*)(fESDEvent->GetMultiplicity());
    Int_t ntracklets = mult->GetNumberOfTracklets();
    for(Int_t i=0;i< ntracklets;i++){
        if(mult->GetDeltaPhi(i)<0.05 && TMath::Abs(mult->GetEta(i))<fEtaCut){
            ntrackletsAccept++;
        }
    }
    
    nTracksTracklets.push_back(nAcceptedTracks);
    nTracksTracklets.push_back(ntrackletsAccept);
    nTracksTracklets.push_back(nAcceptedTracks);//in order to be compatible with mc analysis
    //where here also neutral particles are counted.
    
    
    //get mc vertex for correction maps
    AliGenEventHeader*  header = MCEvent()->GenEventHeader();
    TArrayF mcV;
    header->PrimaryVertex(mcV);
    fVzEvent= mcV[2];
    
    return fNRecAccept; // give back reconstructed value
    
    
}




//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventESDMC(vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                             vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
                                             vector<Float_t> &strangeArray,
                                             vector<Double_t> &nTracksTracklets, const Int_t step)
{
    // gives back the number of charged prim MC particle and pointer to arrays
    // with particle properties (pt, eta, phi)
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    fNMcPrimAccept=0;
    
    AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
    if (!mcEvent) {
        Error("ReadEventESDMC", "Could not retrieve MC event");
        return 0;
    }
    
    AliStack* stack = MCEvent()->Stack();
    if(!stack) return 0;
    
    Int_t ntracks = mcEvent->GetNumberOfTracks();
    if(fDebug>1)Printf("MC particles: %d", ntracks);
    
    //vertex
    AliGenEventHeader*  header = MCEvent()->GenEventHeader();
    TArrayF mcV;
    Float_t vzMC=0.;
    if(header){
        header->PrimaryVertex(mcV);
        vzMC = mcV[2];
        if(step==1){
            fVertexZ[0]->Fill(vzMC);
        }
        fVertexZ[step]->Fill(vzMC);
    }
    
    //----------------------------------
    //first track loop to check how many chared primary tracks are available
    Double_t nChargedPrimaries=0;
    Double_t nAllPrimaries=0;
    
    Double_t nPseudoTracklets=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
        if (!track) {
            Error("ReadEventESDMC", "Could not receive track %d", iTracks);
            continue;
        }
        
        
        if(//count also charged particles in case of fSelectParticles==2 (only neutral)
           !SelectParticlePlusCharged(
                                      track->Charge(),
                                      track->PdgCode(),
                                      stack->IsPhysicalPrimary(track->Label())
                                      )
           )
            continue;
        
        //count number of pseudo tracklets
        if(TMath::Abs(track->Eta())<=fEtaCut && track->Pt()>0.0) ++nPseudoTracklets; //0.035
        //same cuts as on ESDtracks
        if(TMath::Abs(track->Eta())>fEtaCut || track->Pt()<fPtMin || track->Pt()>fPtMax) continue;
        
        //count all primaries
        ++nAllPrimaries;
        //count charged primaries
        if (track->Charge() != 0) ++nChargedPrimaries;
        
        if(fDebug>2) Printf("PDG=%d, IsPrim=%d",  track->PdgCode(),stack->IsPhysicalPrimary(track->Label()) );
        
        
    }
    //----------------------------------
    
    if(fDebug>2){
        Printf("All in acceptance=%f",nAllPrimaries);
        Printf("Charged in acceptance =%f",nChargedPrimaries);
    }
    
    fChargedPi0->Fill(nAllPrimaries,nChargedPrimaries);
    
    if(nAllPrimaries==0) return 0;
    if(nChargedPrimaries==0) return 0;
    
    //track loop
    //Int_t iChargedPiK=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
        if (!track) {
            Error("ReadEventESDMC", "Could not receive track %d", iTracks);
            continue;
        }
        
        if(!SelectParticle(
                           track->Charge(),
                           track->PdgCode(),
                           stack->IsPhysicalPrimary(track->Label())
                           )
           ) 
            continue;
        
        
        //same cuts as on ESDtracks
        if(TMath::Abs(track->Eta())>fEtaCut || track->Pt()<fPtMin  || track->Pt()>fPtMax) continue;
        
        if(fDebug>2) Printf("After: PDG=%d, IsPrim=%d",  track->PdgCode(),stack->IsPhysicalPrimary(track->Label()) );
        
        
        fHistPtMC->Fill(track->Pt());
        //fills arrays with track properties
        ptArray.push_back(track->Pt());
        etaArray.push_back(track->Eta());
        phiArray.push_back(track->Phi());
        chargeArray.push_back(track->Charge());
        strangeArray.push_back(1);
        
    } //track loop
    
    nTracksTracklets.push_back(nChargedPrimaries);
    nTracksTracklets.push_back(nPseudoTracklets);
    if(fSelectParticles!=2){
        nTracksTracklets.push_back(nAllPrimaries);
    }
    else{
        nTracksTracklets.push_back(nAllPrimaries-nChargedPrimaries); // only neutral
    }
    
    fNMcPrimAccept = nChargedPrimaries;
    fNMcPrimAcceptTracklet = nPseudoTracklets;
    
    if(step==1){
        fNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
        fPNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
        fNmcNchTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
        fPNmcNchTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
    }
    if(step==3){
        fNmcNchVtx->Fill(fNMcPrimAccept,fNRecAccept);
        fPNmcNchVtx->Fill(fNMcPrimAccept,fNRecAccept);
        fNmcNchVtxTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
        fPNmcNchVtxTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
    }
    
    fVzEvent= vzMC;
    return fNRecAccept;
    
}

//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventAOD( vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                           vector<Float_t> &phiArray,  vector<Short_t> &chargeArray,
                                           vector<Float_t> &strangeArray,
                                           vector<Double_t> &nTracksTracklets, const Int_t step)
{
    // gives back the number of AOD tracks and pointer to arrays with track
    // properties (pt, eta, phi)
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    TClonesArray *mcArray=0x0;
    if(fAnalysePrimOnly || (fCorrStrangeness && fUseMC)){
        mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
    }
    
    
    AliAODVertex*	vertex= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();//GetPrimaryVertex()
    Double_t vzAOD=vertex->GetZ();
    fVertexZ[step]->Fill(vzAOD);
    
    // Retreive the number of tracks for this event.
    Int_t ntracks = fAODEvent->GetNumberOfTracks();
    if(fDebug>1) Printf("AOD tracks: %d", ntracks);
    
    
    Double_t nAcceptedTracks=0;
    Float_t nAcceptedTracksStrange=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTracks));
        if(!track) AliFatal("Not a standard AOD");
        if (!track) {
            Error("ReadEventAOD", "Could not receive track %d", iTracks);
            continue;
        }
        
        AliVParticle *vtrack = fAODEvent->GetTrack(iTracks);
        
        //use only tracks from primaries
        if(fAnalysePrimOnly){
            if(vtrack->GetLabel()<=0)continue;
            if(!(static_cast<AliAODMCParticle*>(mcArray->At(vtrack->GetLabel()))->IsPhysicalPrimary()))continue;
        }
        
        Double_t save= track->Pt();
        Double_t d0rphiz[2],covd0[3];

	AliAODTrack* clone = (AliAODTrack*) track->Clone("trk_clone"); //need clone, in order not to change track parameters
        Bool_t isDca= clone->PropagateToDCA(fAODEvent->GetPrimaryVertex(),fAODEvent->GetMagneticField(),9999.,d0rphiz,covd0);
	delete clone;
        fPropagateDca->Fill(Int_t(isDca));
        if(TMath::Abs(save - track->Pt())>1e-6) Printf("Before pt=%f, After pt=%f",save, track->Pt());
        
        if(track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta())<fEtaCut
           && track->Pt()>fPtMin && track->Pt()<fPtMax){
            
            nAcceptedTracks++;
            
            // Printf("dca= %f", track->DCA());
            //save track properties in vector
            ptArray.push_back(track->Pt());
            etaArray.push_back(track->Eta());
            phiArray.push_back(track->Phi());
            chargeArray.push_back(track->Charge());
            fHistPt->Fill(track->Pt());
            
            
            //correction for underestimation of strangeness in Monte Carlos -> underestimation of contamination
            Float_t factor=1.;
            if(fUseMC && fCorrStrangeness && step==7){
                if(vtrack->GetLabel()>0){
                    AliAODMCParticle* mcprong =(AliAODMCParticle*)mcArray->At(vtrack->GetLabel());
                    if(mcprong){
                        Int_t labmom = mcprong->GetMother();
                        if(labmom>=0){
                            AliAODMCParticle* mcmother=(AliAODMCParticle*)mcArray->At(labmom);
                            Int_t pdgMother=0;
                            if(mcmother) {
                                pdgMother = mcmother->GetPdgCode();
                                if(TMath::Abs(pdgMother)==310 || TMath::Abs(pdgMother)==130 || TMath::Abs(pdgMother)==321){ //K^0_S, K^0_L, K^+-
                                    if(track->Pt()<=1)factor=1./0.7; // values from strangeness publication
                                    else factor=1./0.6;// values from strangeness publication
                                }
                                if(TMath::Abs(pdgMother)==3122) { //Lambda
                                    factor=1./0.25; // values from strangeness publication
                                }
                            }
                        }
                    }
                }
            }
            nAcceptedTracksStrange+=factor;
            strangeArray.push_back(factor);
            fDcaXY[step]->Fill(d0rphiz[0], factor);
            fDcaZ[step]->Fill(d0rphiz[0], factor);
            
        }
    }
    //need to check this option for MC
    if(nAcceptedTracks==0) return 0;
    
    
    //tracklet loop
    Int_t ntrackletsAccept=0;
    AliAODTracklets * mult= (AliAODTracklets*)fAODEvent->GetTracklets();
    for(Int_t i=0;i<mult->GetNumberOfTracklets();++i){
        if(TMath::Abs(mult->GetDeltaPhi(i))<0.05 && TMath::Abs(TMath::Log(TMath::Tan(0.5 * mult->GetTheta(i))))<fEtaCut){
            ++ntrackletsAccept;
        }
    }
    
    
    nTracksTracklets.push_back(nAcceptedTracks);
    nTracksTracklets.push_back(ntrackletsAccept);
    nTracksTracklets.push_back(nAcceptedTracks);//in order to be compatible with mc analysis 
    //where here also neutral particles are counted.
    
    
    fVzEvent= vzAOD;
    if(step==5 || step==2){
        fNRecAccept = nAcceptedTracks; // needed for MC case //step5 = TrigVtxRecNrec
        fNRecAcceptTracklet = ntrackletsAccept; // needed for MC case //step5 = TrigVtxRecNrec
    }
    if(step==7)fNRecAcceptStrangeCorr = nAcceptedTracksStrange;
    
    return fNRecAccept; // at the moment, always return reconstructed multiplicity
    
}   

//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventAODRecMcProp( vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                                    vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
                                                    vector<Float_t> &strangeArray,
                                                    vector<Double_t> &nTracksTracklets, const Int_t step)
{
    // gives back the number of AOD tracks and pointer to arrays with track
    // properties (pt, eta, phi)
    
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    
    // Retreive the number of tracks for this event.
    Int_t ntracks = fAODEvent->GetNumberOfTracks();
    if(fDebug>1) Printf("AOD tracks: %d", ntracks);
    
    
    //get array of mc particles
    TClonesArray *mcArray = (TClonesArray*)fAODEvent->
    FindListObject(AliAODMCParticle::StdBranchName());
    if(!mcArray){
        AliInfo("No MC particle branch found");
        return kFALSE;
    }
    
    AliAODVertex*	vtx= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();//GetPrimaryVertex()
    Double_t vzAOD=vtx->GetZ();
    fVertexZ[step]->Fill(vzAOD);
    
    Double_t nAcceptedTracks=0;
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTracks));
        if(!track) AliFatal("Not a standard AOD");
        
        AliVParticle *vtrack = fAODEvent->GetTrack(iTracks);
        
        if (!track) {
            Error("ReadEventAODRecMcProp", "Could not receive track %d", iTracks);
            continue;
        }
        
        //use only tracks from primaries
        if(fAnalysePrimOnly){
            if(vtrack->GetLabel()<=0)continue;
            if(!(static_cast<AliAODMCParticle*>(mcArray->At(vtrack->GetLabel()))->IsPhysicalPrimary()))continue;
        }
        
        if(track->TestFilterBit(fFilterBit) &&  TMath::Abs(track->Eta())<fEtaCut &&
           track->Pt()>fPtMin && track->Pt()<fPtMax){
            
            nAcceptedTracks++;
            Float_t factor=1.;
            
            //save track properties in vector
            if(vtrack->GetLabel()<=0){ //fake tracks before "label<0", but crash in AOD079 // what is the meaning of label 0
                // 	Printf("Fake track");
                // 	continue;
                ptArray.push_back(track->Pt());
                etaArray.push_back(track->Eta());
                phiArray.push_back(track->Phi());
                chargeArray.push_back(track->Charge());
                
            }
            else{//mc properties
                AliAODMCParticle *partOfTrack = (AliAODMCParticle*)mcArray->At(vtrack->GetLabel());
                if(!partOfTrack) Printf("label=%d", vtrack->GetLabel());
                if(partOfTrack){
                    ptArray.push_back(partOfTrack->Pt());
                    etaArray.push_back(partOfTrack->Eta());
                    phiArray.push_back(partOfTrack->Phi());
                    chargeArray.push_back(vtrack->Charge());//partOfTrack?
                    
                    //correction for underestimation of strangeness in Monte Carlos -> underestimation of contamination
                    if(fUseMC && fCorrStrangeness && step==6){
                        Int_t labmom = partOfTrack->GetMother();
                        if(labmom>=0){
                            AliAODMCParticle* mcmother=(AliAODMCParticle*)mcArray->At(labmom);
                            Int_t pdgMother=0;
                            if(mcmother) {
                                pdgMother = mcmother->GetPdgCode();
                                if(TMath::Abs(pdgMother)==310 || TMath::Abs(pdgMother)==130 || TMath::Abs(pdgMother)==321){ //K^0_S, K^0_L, K^+-
                                    if(track->Pt()<=1)factor=1./0.7; // values from strangeness publication
                                    else factor=1./0.6;// values from strangeness publication
                                }
                                if(TMath::Abs(pdgMother)==3122) { //Lambda
                                    factor=1./0.25; // values from strangeness publication
                                }
                            }
                        }
                    }
                }
            }
            strangeArray.push_back(factor);
            
        }
    }
    //need to check this option for MC
    if(nAcceptedTracks==0) return 0;
    
    //tracklet loop
    Int_t ntrackletsAccept=0;
    AliAODTracklets * mult= (AliAODTracklets*)fAODEvent->GetTracklets();
    for(Int_t i=0;i<mult->GetNumberOfTracklets();++i){
        if(TMath::Abs(mult->GetDeltaPhi(i))<0.05 && TMath::Abs(TMath::Log(TMath::Tan(0.5 * mult->GetTheta(i))))<fEtaCut ){
            ++ntrackletsAccept;
        }
    }
    
    
    nTracksTracklets.push_back(nAcceptedTracks);
    nTracksTracklets.push_back(ntrackletsAccept);
    nTracksTracklets.push_back(nAcceptedTracks);//in order to be compatible with mc analysis
    //where here also neutral particles are counted.
    
    
    //check vertex (mc)
    AliAODMCHeader *aodMCheader = (AliAODMCHeader *) fAODEvent->
    FindListObject(AliAODMCHeader::StdBranchName());
    Float_t vzMC = aodMCheader->GetVtxZ();
    
    fVzEvent= vzMC;
    return fNRecAccept;//this is the rec value from step 5
    
}



//________________________________________________________________________
Double_t AliAnalysisTaskMinijetV2::ReadEventAODMC( vector<Float_t> &ptArray,  vector<Float_t> &etaArray,
                                             vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
                                             vector<Float_t> &strangeArray,
                                             vector<Double_t> &nTracksTracklets, const Int_t step)
{
    // gives back the number of AOD MC particles and pointer to arrays with particle
    // properties (pt, eta, phi)
    
    ptArray.clear();
    etaArray.clear();
    phiArray.clear();
    chargeArray.clear();
    strangeArray.clear();
    nTracksTracklets.clear();
    
    //check vertex
    AliAODMCHeader *aodMCheader = (AliAODMCHeader *) fAODEvent->
    FindListObject(AliAODMCHeader::StdBranchName());
    Float_t vzMC = aodMCheader->GetVtxZ();
    if(step==1){
        fVertexZ[0]->Fill(vzMC);
    }
    fVertexZ[step]->Fill(vzMC);
    
    //retreive MC particles from event
    TClonesArray *mcArray = (TClonesArray*)fAODEvent->
    FindListObject(AliAODMCParticle::StdBranchName());
    if(!mcArray){
        AliInfo("No MC particle branch found");
        return kFALSE;
    }
    
    Int_t ntracks = mcArray->GetEntriesFast();
    if(fDebug>1)Printf("MC particles: %d", ntracks);
    
    
    // Track loop: chek how many particles will be accepted
    //Float_t vzMC=0.;
    Double_t nChargedPrim=0;
    Double_t nAllPrim=0;
    Double_t nPseudoTracklets=0;
    for (Int_t it = 0; it < ntracks; it++) {
        AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
        if (!track) {
            Error("ReadEventAODMC", "Could not receive particle %d", it);
            continue;
        }
        
        if(!SelectParticlePlusCharged(
                                      track->Charge(),
                                      track->PdgCode(),
                                      track->IsPhysicalPrimary()
                                      )
           )
            continue;
        
        if(TMath::Abs(track->Eta())<fEtaCut && track->Pt()>0.0)++nPseudoTracklets; //0.035
        if(TMath::Abs(track->Eta())>fEtaCut || track->Pt()<fPtMin || track->Pt()>fPtMax) continue;
        
        nAllPrim++;
        if(track->Charge()!=0) nChargedPrim++;
        
    }
    
    
    if(nAllPrim==0) return 0;
    if(nChargedPrim==0) return 0;
    
    //generate array with size of number of accepted tracks
    fChargedPi0->Fill(nAllPrim,nChargedPrim);
    
    
    // Track loop: fill arrays for accepted tracks
    // Int_t iChargedPrim=0;
    for (Int_t it = 0; it < ntracks; it++) {
        AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
        if (!track) {
            Error("ReadEventAODMC", "Could not receive particle %d", it);
            continue;
        }
        
        if(!SelectParticle(
                           track->Charge(),
                           track->PdgCode(),
                           track->IsPhysicalPrimary()
                           )
           ) 
            continue;
        
        if(TMath::Abs(track->Eta())>fEtaCut || track->Pt()<fPtMin || track->Pt()>fPtMax) continue;
        
        fHistPtMC->Fill(track->Pt());
        ptArray.push_back(track->Pt());
        etaArray.push_back(track->Eta());
        phiArray.push_back(track->Phi());
        chargeArray.push_back(track->Charge());
        strangeArray.push_back(1);
    }
    
    nTracksTracklets.push_back(nChargedPrim);
    nTracksTracklets.push_back(nPseudoTracklets);
    if(fSelectParticles!=2){
        nTracksTracklets.push_back(nAllPrim);
    }
    else{
        nTracksTracklets.push_back(nAllPrim-nChargedPrim); // only neutral
    }
    
    
    
    fVzEvent= vzMC;
    fNMcPrimAccept = nChargedPrim;
    fNMcPrimAcceptTracklet = nPseudoTracklets;
    
    if(step==1){ // step 1 = Trig All Mc Nmc
        fNmcNch->Fill( fNMcPrimAccept,fNRecAccept);
        fPNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
        fNmcNchTracklet->Fill( fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
        fPNmcNchTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
    }
    if(step==3){ // step 3 = Trig vtx Mc
        fNmcNchVtx->Fill( fNMcPrimAccept,fNRecAccept);
        fNmcNchVtxStrangeCorr->Fill( fNMcPrimAccept,fNRecAcceptStrangeCorr);
        fPNmcNchVtx->Fill(fNMcPrimAccept,fNRecAccept);
        fNmcNchVtxTracklet->Fill( fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
        fPNmcNchVtxTracklet->Fill(fNMcPrimAcceptTracklet,fNRecAcceptTracklet);
    }
    return fNRecAccept; // rec value from step 5 or step 2
    
    
} 

//________________________________________________________________________
void AliAnalysisTaskMinijetV2::Analyse(const vector<Float_t> &pt,
                                     const vector<Float_t> &eta,
                                     const vector<Float_t> &phi,
                                     const vector<Short_t> &charge,
                                     const vector<Float_t> &strangeWeight,
                                     const Double_t ntracksCharged,
                                     const Int_t ntracklets,
                                     const Int_t nAll,
                                     const Int_t step)
{
    
    // analyse track properties (comming from either ESDs or AODs) in order to compute
    // mini jet activity (chared tracks) as function of charged multiplicity
    
    fStep->Fill(step);
    
    if(fDebug){
        Printf("Analysis Step=%d", step);
        if(fDebug>2){
            Printf("nAll=%d",nAll);
            Printf("nCharged=%f",ntracksCharged);
            for (Int_t i = 0; i < nAll; i++) {
                Printf("pt[%d]=%f",i,pt[i]);
            }
        }
    }
    
    //calculation of mean pt for all tracks in case of step==0
    if(step==5 || step ==2){ // rec step
        Double_t meanPt=0.;
        Double_t leadingPt=0.;
        for (Int_t i = 0; i < nAll; i++) {
            if(pt[i]<0.01)continue;
            meanPt+=pt[i];
            if(leadingPt<pt[i])leadingPt=pt[i];
        }
        meanPt=meanPt/nAll;
        fMeanPtRec=meanPt;
        fLeadingPtRec=leadingPt;
    }
    
    Double_t propEvent[4] = {ntracksCharged,fVzEvent,fMeanPtRec, fLeadingPtRec}; //vz: {rec, mc, mc}, meanPt and Nrec is always rec value
    fMapEvent[step]->Fill(propEvent);
    
    fNcharge[step]->Fill(ntracksCharged);
    
    Float_t ptEventAxis=0;  // pt event axis
    Float_t etaEventAxis=0; // eta event axis
    Float_t phiEventAxis=0; // phi event axis
    Short_t chargeEventAxis=0; // charge event axis
    Float_t strangeWeightEventAxis=0;  // strange weight event axis
    
    Float_t ptOthers  = 0; // pt others // for all other tracks around event axis -> see loop
    Float_t etaOthers = 0; // eta others
    Float_t phiOthers = 0; // phi others
    Short_t chargeOthers = 0; // charge others
    Float_t strangeWeightOthers  = 0; // strange weight others
    
    Int_t   *pindexInnerEta  = new Int_t  [nAll+1];
    Float_t *ptInnerEta      = new Float_t[nAll+1];
    
    
    
    for (Int_t i = 0; i < nAll; i++) {
        
        if(pt[i]<0.01)continue;
        
        //fill single particle correction for first step of pair correction
        Double_t propAll[3] = {eta[i],pt[i],ntracksCharged };
        fMapAll[step]->Fill(propAll, strangeWeight[i]);
        
        
        //filling of simple check plots
        if(pt[i]<0.7) continue;
        fPt[step]    -> Fill( pt[i]);
        fEta[step]   -> Fill(eta[i]);
        fPhi[step]   -> Fill(phi[i]);
        fPhiEta[step]-> Fill(phi[i], eta[i]);
        
        pindexInnerEta[i]=0; //set all values to zero
        //fill new array for eta check
        ptInnerEta[i]= pt[i];
        
    }
    
    
    
    // define event axis: leading or random track (with pt>fTriggerPtCut)
    // ---------------------------------------
    Int_t highPtTracks=0;
    Int_t highPtTracksInnerEta=0;
    // Double_t highPtTracksInnerEtaCorr=0;
    Int_t mult09=0;
    
    //count high pt tracks and high pt tracks in acceptance for seeds
    for(Int_t i=0;i<nAll;i++){
        
        if(pt[i]<0.01)continue;
        
        if(TMath::Abs(eta[i])<0.9){
            mult09++;
        }
        
        if(pt[i]>fTriggerPtCut) {
            highPtTracks++;
        }
        
        // seed should be place in middle of acceptance, that complete cone is in acceptance
        if(pt[i]>fTriggerPtCut && TMath::Abs(eta[i])<fEtaCutSeed && charge[i]!=0){
            
            // Printf("eta=%f", eta[i]);
            highPtTracksInnerEta++;
            
        }
        else{
            ptInnerEta[i]=0;
        }
    }
    
    
    //sort array in order to get highest pt tracks first
    //index can be computed with pindexInnerEta[number]
    if(nAll) TMath::Sort(nAll, ptInnerEta, pindexInnerEta, kTRUE);
    
    //     plot of multiplicity distributions
    fNch07Nch[step]->Fill(ntracksCharged, highPtTracksInnerEta);
    fPNch07Nch[step]->Fill(ntracksCharged, highPtTracksInnerEta);
    
    if(ntracklets){
        fNch07Tracklet[step]->Fill(ntracklets, highPtTracksInnerEta);//only counts tracks which can be used as seeds
        fNchTracklet[step]->Fill(ntracklets, ntracksCharged);
        fPNch07Tracklet[step]->Fill(ntracklets, highPtTracksInnerEta);//only counts tracks which can be used as seeds
    }
    
    //analysis can only be performed with event axis, defined by high pt track
    
    
    if(highPtTracks>0 && highPtTracksInnerEta>0){
        
        // build pairs in two track loops
        // loop over all possible trigger particles (defined by pt_trig and eta_acceptance)
        for(Int_t axis=0;(axis<nAll) && (pt[pindexInnerEta[axis]]>=fTriggerPtCut); axis++){
            
            //EventAxisRandom track properties
            ptEventAxis  = pt [pindexInnerEta[axis]];
            etaEventAxis = eta[pindexInnerEta[axis]];
            phiEventAxis = phi[pindexInnerEta[axis]];
            chargeEventAxis = charge[pindexInnerEta[axis]];
            strangeWeightEventAxis = strangeWeight[pindexInnerEta[axis]];
            fPtSeed[step]    -> Fill( ptEventAxis);
            fEtaSeed[step]   -> Fill(etaEventAxis);
            fPhiSeed[step]   -> Fill(phiEventAxis);
            
            
            Double_t prop[3] = {etaEventAxis,ptEventAxis,ntracksCharged };
            fMapSingleTrig[step]->Fill(prop, strangeWeightEventAxis);
            
            //associated tracks
            for (Int_t iTrack = axis+1; iTrack < nAll; iTrack++) {
                
                if(pt[pindexInnerEta[iTrack]]<fAssociatePtCut) continue;
                
                if(fSelectParticlesAssoc==1){
                    if(charge[pindexInnerEta[iTrack]]==0)continue;
                }
                if(fSelectParticlesAssoc==2){
                    if(charge[pindexInnerEta[iTrack]]!=0)continue;
                }
                
                
                ptOthers   = pt [pindexInnerEta[iTrack]];
                etaOthers  = eta[pindexInnerEta[iTrack]];
                phiOthers  = phi[pindexInnerEta[iTrack]];
                chargeOthers = charge[pindexInnerEta[iTrack]];
                strangeWeightOthers = strangeWeight[pindexInnerEta[iTrack]];
                
                
                //plot only properties of tracks with pt>ptassoc
                fPtOthers[step]    -> Fill( ptOthers);
                fEtaOthers[step]   -> Fill(etaOthers);
                fPhiOthers[step]   -> Fill(phiOthers);
                fPtEtaOthers[step]   -> Fill(ptOthers, etaOthers);
                
                //	if(fDebug>2)Printf("%f, %f", pt[pindexInnerEta[axis]], pt[pindexInnerEta[iTrack]]);
                
                Float_t dPhi = phiOthers-phiEventAxis;
                if(dPhi>1.5*TMath::Pi()) dPhi = dPhi-2*TMath::Pi();
                else if(dPhi<-0.5*TMath::Pi())dPhi=dPhi+2*TMath::Pi();
                Float_t dEta=etaOthers-etaEventAxis;
                
                
                fDPhiDEtaEventAxis[step]->Fill(dPhi, dEta, strangeWeightEventAxis*strangeWeightOthers);
                fDPhiEventAxis[step]->Fill(dPhi, strangeWeightEventAxis*strangeWeightOthers);
                
                //check outliers
                if(ptEventAxis< 0.4 || ptEventAxis > 100) AliInfo("particles out of range pt");
                if(ntracksCharged<-1 || ntracksCharged>1500) AliInfo("particles out of range ncharge");
                if(TMath::Abs(dEta)>2*fEtaCut) {
                    AliInfo("particles out of range dEta");
                    AliInfo(Form("eta1=%f, eta2=%f", etaOthers, etaEventAxis));
                    AliInfo(Form("step=%d",step));
                }
                if(dPhi<-0.5*TMath::Pi() || dPhi>1.5*TMath::Pi()){
                    AliInfo(Form("particles out of range dPhi"));
                    AliInfo(Form("phi1=%f, phi2=%f", phiOthers, phiEventAxis));
                }
                
                Bool_t isLikeSign = CheckLikeSign(chargeEventAxis, chargeOthers);
                
                Double_t prop6[6] = {ptEventAxis,ptOthers,dEta,dPhi,ntracksCharged, static_cast<Double_t>(isLikeSign) };
                fMapPair[step]->Fill(prop6, strangeWeightEventAxis*strangeWeightOthers);
                
                //thrid track loop (Andreas: three particle auto-correlations)
                if(fThreeParticleCorr){
                    for (Int_t third = iTrack+1; third < nAll; third++) {
                        if(pt[pindexInnerEta[iTrack]]<fTriggerPtCut) continue;
                        if(pt[pindexInnerEta[third]]<fTriggerPtCut) continue;
                        
                        //dphi1
                        Float_t dPhi1 = phiEventAxis - phiOthers;
                        if(dPhi1>1.5*TMath::Pi()) dPhi1 = dPhi1-2*TMath::Pi();
                        else if(dPhi1<-0.5*TMath::Pi())dPhi1=dPhi1+2*TMath::Pi();
                        
                        Float_t phiThird = phi[pindexInnerEta[third]];
                        Float_t strangeWeightThird = strangeWeight[pindexInnerEta[third]];
                        
                        //dphi2
                        Float_t dPhi2 = phiEventAxis - phiThird;
                        if(dPhi2>1.5*TMath::Pi()) dPhi2 = dPhi2-2*TMath::Pi();
                        else if(dPhi2<-0.5*TMath::Pi())dPhi2=dPhi2+2*TMath::Pi();
                        
                        fDPhi1DPhi2[step]-> Fill(dPhi1, dPhi2, strangeWeightEventAxis*strangeWeightOthers*strangeWeightThird);
                        Double_t propThree[3] = {dPhi1,dPhi2,ntracksCharged}; 
                        fMapThree[step]->Fill(propThree,strangeWeightEventAxis*strangeWeightOthers*strangeWeightThird );
                        
                        
                    }// end of three particle correlation loop
                    
                }// if fThreeParticleCorr is set to true
                
            }// end of inner track loop
            
        } //end of outer track loop
        
        fTriggerNch[step]->Fill(ntracksCharged,highPtTracksInnerEta);
        fTriggerNchSeeds[step]->Fill(ntracksCharged,highPtTracksInnerEta);
        fTriggerTracklet[step]->Fill(ntracklets);
        
        
    }//if there is at least one high pt track
    
    
    if(pindexInnerEta){// clean up array memory used for TMath::Sort
        delete[] pindexInnerEta; 
        pindexInnerEta=0;
    }
    
    if(ptInnerEta){// clean up array memory used for TMath::Sort
        delete[] ptInnerEta; 
        ptInnerEta=0;
    }
    
}



//________________________________________________________________________
void AliAnalysisTaskMinijetV2::Terminate(Option_t*)
{
    //terminate function is called at the end
    //can be used to draw histograms etc.
    
    
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMinijetV2::SelectParticlePlusCharged(const Short_t charge, const Int_t pdg, Bool_t prim)
{
    //selection of mc particle
    //fSelectParticles=0: use charged primaries and pi0 and k0
    //fSelectParticles=1: use only charged primaries
    //fSelectParticles=2: use only pi0 and k0
    
    if(fSelectParticles==0 || fSelectParticles==2){ // in case of 2: need to count also charged particles
        if(charge==0){
            if(!(pdg==111||pdg==130||pdg==310))
                return false;
        }
        else{// charge !=0
            if(!prim)
                return false;
        }
    }
    
    else if(fSelectParticles==1){
        if (charge==0 || !prim){
            return false;
        }
    }
    
    else{
        AliInfo("Error: wrong selection of charged/pi0/k0");
        return 0;
    }
    
    return true;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMinijetV2::SelectParticle(const Short_t charge, const Int_t pdg, const Bool_t prim)
{
    //selection of mc particle
    //fSelectParticles=0: use charged primaries and pi0 and k0
    //fSelectParticles=1: use only charged primaries
    //fSelectParticles=2: use only pi0 and k0
    
    if(fSelectParticles==0){
        if(charge==0){
            if(!(pdg==111||pdg==130||pdg==310))
                return false;
        }
        else{// charge !=0
            if(!prim)
                return false;
        }
    }
    else if (fSelectParticles==1){
        if (charge==0 || !prim){
            return false;
        }
    }
    else if(fSelectParticles==2){
        if(!(pdg==111||pdg==130||pdg==310))
            return false;
    }
    
    return true;
    
}

//________________________________________________________________________
Bool_t AliAnalysisTaskMinijetV2::CheckEvent(const Bool_t recVertex)
{
    // This function tests the quality of an event (ESD/AOD) (rec and/or mc part)
    // recVertex=false:  check if Mc events and stack is there, Nmc>0
    // recVertex=false: " + check if there is a good, reconstructed SPD vertex
    // defined by |z|<fVertexCut(10cm), Contributer>0, no PileUpFromSPD(3,0,8)
    
    if(fMode==0){//esd
        
        //mc
        if(fUseMC){
            
            //mc event
            AliMCEvent *mcEvente = (AliMCEvent*) MCEvent();
            if (!mcEvente) {
                Error("CheckEvent", "Could not retrieve MC event");
                return false;
            }
            
            //stack
            AliStack* stackg = MCEvent()->Stack();
            if(!stackg) return false;
            Int_t ntracksg = mcEvente->GetNumberOfTracks();
            if(ntracksg<0) return false;
            
            //vertex
            AliGenEventHeader*  headerg = MCEvent()->GenEventHeader();
            TArrayF mcVg;
            headerg->PrimaryVertex(mcVg);
            //if(TMath::Abs(mcVg[0])<1e-8 && TMath::Abs(mcVg[0])<1e-8 &&
            //   TMath::Abs(mcVg[0])<1e-8) return false;
            Float_t vzMCg = mcVg[2];
            if(TMath::Abs(vzMCg)>fVertexZCut) return false;
            //hasVtxMc=true;
        }
        
        //rec
        if(recVertex==true){
            if (fAnalysisUtils && fAnalysisUtils->IsPileUpMV(fESDEvent)) {
              fEventStat->Fill(9);
              return false;
            }
            
            //rec vertex
            const AliESDVertex*	vertexESDg   = fESDEvent->GetPrimaryVertex(); // uses track or SPD vertexer
            if(!vertexESDg) return false;
            fVertexCheck->Fill(vertexESDg->GetNContributors());
            if(vertexESDg->GetNContributors()<=0)return false;
            Float_t fVzg= vertexESDg->GetZ();
            if(TMath::Abs(fVzg)>fVertexZCut) return false;
            
            //rec spd vertex
            const AliESDVertex *vtxSPD = fESDEvent->GetPrimaryVertexSPD();
            if(!vtxSPD) return false;
            if(vtxSPD->GetNContributors()<=0)return false;
            Float_t fVzSPD= vtxSPD->GetZ();
            if(TMath::Abs(fVzSPD)>fVertexZCut) return false;
            
        }
        return true;
    }
    
    
    else if(fMode==1){ //aod
        
        if(fUseMC){
            
            //retreive MC particles from event
            TClonesArray *mcArray = (TClonesArray*)fAODEvent->
            FindListObject(AliAODMCParticle::StdBranchName());
            if(!mcArray){
                AliInfo("No MC particle branch found");
                return false;
            }
            
            //mc
            AliAODMCHeader *aodMCheader = (AliAODMCHeader *) fAODEvent->
            FindListObject(AliAODMCHeader::StdBranchName());
            Float_t vzMC = aodMCheader->GetVtxZ();
            if(TMath::Abs(vzMC)>fVertexZCut) return false;
            
            //hasVtxMc=true;
        }
        
        //rec
        if(recVertex==true){
            // pile up
            if (fAnalysisUtils && fAnalysisUtils->IsPileUpMV(fAODEvent)) {
              fEventStat->Fill(9);
              return false;
            }
            
            AliAODVertex*	vertex= (AliAODVertex*)fAODEvent->GetPrimaryVertex();
            if(!vertex) return false;
            TString vtxTitle(vertex->GetTitle());// only allow vertex from tracks, no vertexer z
            // Printf("vtxTitle: %s",vtxTitle.Data());
            //if (!(vtxTitle.Contains("VertexerTracksWithConstraint"))) return false;
            fVertexCheck->Fill(vertex->GetNContributors());
            if(vertex->GetNContributors()<=0) return false;
            Double_t vzAOD=vertex->GetZ();
            // if(TMath::Abs(vzAOD)<1e-9) return false;
            if(TMath::Abs(vzAOD)>fVertexZCut) return false;
            
            AliAODVertex*	vertexSPD= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();
            if(!vertexSPD) return false;
            if(vertexSPD->GetNContributors()<=0) return false;
            Double_t vzSPD=vertexSPD->GetZ();
            //if(TMath::Abs(vzSPD)<1e-9) return false;
            if(TMath::Abs(vzSPD)>fVertexZCut) return false;
            
            //check TPC reconstruction: check for corrupted chunks
            //better: check TPCvertex, but this is not available yet in AODs
            Int_t nAcceptedTracksTPC=0;
            Int_t nAcceptedTracksITSTPC=0;
            for (Int_t iTracks = 0; iTracks < fAODEvent->GetNumberOfTracks(); iTracks++) {
                AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTracks));
                if(!track) AliFatal("Not a standard AOD");
                if (!track) continue;
                if(track->TestFilterBit(128) && TMath::Abs(track->Eta())<fEtaCut &&
                   track->Pt()>fPtMin && track->Pt()<fPtMax)
                    nAcceptedTracksTPC++;
                if(track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta())<fEtaCut &&
                   track->Pt()>fPtMin && track->Pt()<fPtMax)
                    nAcceptedTracksITSTPC++;
            }
            fCorruptedChunks->Fill(nAcceptedTracksTPC,nAcceptedTracksITSTPC);
            if(fRejectChunks){
                if(nAcceptedTracksTPC>fNTPC && nAcceptedTracksITSTPC==0)
                    return false;//most likely corrupted chunk. No ITS tracks are reconstructed
            }
            fCorruptedChunksAfter->Fill(nAcceptedTracksTPC,nAcceptedTracksITSTPC);
            
            //control histograms=================
            //tracklet loop
            Int_t ntrackletsAccept=0;
            AliAODTracklets * mult= (AliAODTracklets*)fAODEvent->GetTracklets();
            for(Int_t i=0;i<mult->GetNumberOfTracklets();++i){
                if(TMath::Abs(mult->GetDeltaPhi(i))<0.05 &&
                   TMath::Abs(TMath::Log(TMath::Tan(0.5 * mult->GetTheta(i))))<fEtaCut) ++ntrackletsAccept;
            }
            Int_t nAcceptedTracks=0;
            for (Int_t iTracks = 0; iTracks < fAODEvent->GetNumberOfTracks(); iTracks++) {
                AliAODTrack *track = dynamic_cast<AliAODTrack*>(fAODEvent->GetTrack(iTracks));
                if(!track) AliFatal("Not a standard AOD");
                if (!track) continue;
                if(track->TestFilterBit(fFilterBit) && TMath::Abs(track->Eta())<fEtaCut
                   && track->Pt()>fPtMin && track->Pt()<fPtMax) nAcceptedTracks++;
            }
            fNContrNtracklets->Fill(ntrackletsAccept,vertexSPD->GetNContributors());
            fNContrNtracks->Fill(nAcceptedTracks,vertexSPD->GetNContributors());
            //====================================
        }
        return true;
        
    }
    
    else {
        Printf("Wrong mode!");
        return false;
    }
    
}

//_____________________________________________________________________________
const Double_t * AliAnalysisTaskMinijetV2::CreateLogAxis(const Int_t nbins,
                                                       const Double_t xmin,
                                                       const Double_t xmax)
{
    // returns pointer to an array with can be used to build a logarithmic axis
    // it is user responsibility to delete the array
    
    Double_t logxmin = TMath::Log10(xmin);
    Double_t logxmax = TMath::Log10(xmax);
    Double_t binwidth = (logxmax-logxmin)/nbins;
    
    Double_t *xbins =  new Double_t[nbins+1];
    
    xbins[0] = xmin;
    for (Int_t i=1;i<=nbins;i++) {
        xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
    }
    
    return xbins;
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskMinijetV2::CheckLikeSign(const Short_t chargeEventAxis, 
                                             const Short_t chargeOthers)
{
    // compute if charge of two particles/tracks has same sign or different sign
    
    if(fDebug>2) Printf("Charge1=%d, Charge2=%d",chargeEventAxis,chargeOthers);
    
    if(chargeEventAxis<0){
        if(chargeOthers<0){
            return true;
        }
        else if(chargeOthers>0){
            return false;
        }
    }
    
    else if(chargeEventAxis>0){
        if(chargeOthers>0){
            return true;
        }
        else if(chargeOthers<0){
            return false;
        }
    }
    
    else{
        AliInfo("Error: Charge not lower nor higher as zero");
        return false;
    }
    
    AliInfo("Error: Check values of Charge");
    return false;
}


