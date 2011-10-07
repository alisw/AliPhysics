#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

#include "TCanvas.h"
#include "TList.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TProfile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

#include "AliMCParticle.h"
#include "AliMCEvent.h"
#include "AliAnalysisTaskEfficiency.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

// Analysis Task for basic QA on the ESD
// Authors: Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing, Veronica Canoa, AM, Eva Sicking

ClassImp(AliAnalysisTaskEfficiency)

//________________________________________________________________________
  AliAnalysisTaskEfficiency::AliAnalysisTaskEfficiency(const char *name) 
    : AliAnalysisTaskSE(name) 
    ,fFieldOn(kTRUE)
    ,fHists(0)
    ,fHistRECpt(0)
    ,fHistMCpt(0)
    ,fHistMCNRpt(0)
    ,fHistFAKEpt(0)
    ,fHistMULTpt(0)
    ,fHistRECeta(0)
    ,fHistMCeta(0)
    ,fHistMCNReta(0)
    ,fHistFAKEeta(0)
    ,fHistMULTeta(0)
    ,fHistRECphi(0)
    ,fHistMCphi(0)
    ,fHistMCNRphi(0)
    ,fHistFAKEphi(0)
    ,fHistMULTphi(0)
    ,fHistRECHPTeta(0)
    ,fHistMCHPTeta(0)
    ,fHistMCNRHPTeta(0)
    ,fHistFAKEHPTeta(0)
    ,fHistRECHPTphi(0)
    ,fHistMCHPTphi(0)
    ,fHistMCNRHPTphi(0)
    ,fHistFAKEHPTphi(0)
    ,fHistRecMult(0)
    ,fHistNCluster(0)  
    ,fh1VtxEff(0) 
    ,fh2VtxTpcSpd(0)
    ,fh2EtaCorrelation(0)      
    ,fh2EtaCorrelationShift(0)  
    ,fh2PhiCorrelation(0)      
    ,fh2PhiCorrelationShift(0)  
    ,fh2PtCorrelation(0)      
    ,fh2PtCorrelationShift(0)
    ,fHistDeltaphiprimaries(0)
    ,fHistDeltaphisecondaries(0)           
    ,fHistDeltaphireject(0) 
    ,fHistTPCITSdeltax(0)
    ,fHistTPCITSdeltay(0)
    ,fHistTPCITSdeltaz(0)
    ,fHistTPCITSdeltar(0)
    ,fHistTPCTRDdeltax(0)
    ,fHistTPCTRDdeltay(0)
    ,fHistTPCTRDdeltaz(0)
    ,fHistTPCTRDdeltar(0)
    ,fHistTPCTRDdeltaphi(0)
    ,fPtPullsPos(0)
    ,fPtPullsNeg(0)
    ,fPtShiftPos(0)
    ,fPtShiftNeg(0)
    ,fTrackType(0)
    ,fCuts(0)
    ,fVertexRvsZ(0)
    ,fVertexRvsZC(0)
    ,fEtaMultiFluc(0)
    ,fEtaMulti(0)
    ,fEtaMultiH(0)
    ,fPhiMultiFluc(0)
    ,fPhiMulti(0)
    ,fPhiMultiH(0)
{
  //
  // Constructor
  //

  for(int i = 0;i< 2;++i){
    fh2MultSpdChips[i] = 0;
  }

  for(int i = 0;i < kTotalAC;++i){
    fh2PhiPadRow[i] = fh2PhiLayer[i] = 0;
  }

  for(int i = 0; i < kTotalVtx; ++i){ 
    fh2VertexCorrelation[i] = 
      fh2VertexCorrelationShift[i] = 0;
    fh1VertexShift[i] =         
      fh1VertexShiftNorm[i] = 0;      
  }

  for(int i = 0;i< 8;++i){
    fHistRECptCharge[i]=0;
    fHistMCptCharge[i]=0;
    fHistMCNRptCharge[i]=0;
    fHistFAKEptCharge[i]=0;

    fHistRECetaCharge[i]=0;
    fHistMCetaCharge[i]=0;
    fHistMCNRetaCharge[i]=0;
    fHistFAKEetaCharge[i]=0;

    fHistRECphiCharge[i]=0;
    fHistMCphiCharge[i]=0;
    fHistMCNRphiCharge[i]=0;
    fHistFAKEphiCharge[i]=0;

  }

  DefineOutput(1,  TList::Class());
}


//________________________________________________________________________
void AliAnalysisTaskEfficiency::UserCreateOutputObjects()
{
  // Create histograms


  fHists = new TList();
  fHistRECpt   = new TH1F("fHistRECpt",  " p_{T} distribution: all reconstructed",        100, 0., 20.);
  fHistMCpt    = new TH1F("fHistMCpt",   " p_{T} distribution: all MC",                   100, 0., 20.);
  fHistMCNRpt  = new TH1F("fHistMCNRpt", " p_{T} distribution: all not-reconstructed MC", 100, 0., 20.);
  fHistFAKEpt  = new TH1F("fHistFAKEpt", " p_{T} distribution: all Fake",                 100, 0., 20.);
  fHistMULTpt  = new TH1F("fHistMULTpt", " p_{T} distribution: multiply rec.",            100, 0., 20.);
    
  fHistRECeta   = new TH1F("fHistRECeta",  " #eta distribution: all reconstructed",        100,-1.0,1.0);
  fHistMCeta    = new TH1F("fHistMCeta",   " #eta distribution: all MC",                   100,-1.0,1.0);
  fHistMCNReta  = new TH1F("fHistMCNReta", " #eta distribution: all not-reconstructed MC", 100,-1.0,1.0);
  fHistFAKEeta  = new TH1F("fHistFAKEeta", " #eta distribution: all Fake",                 100,-1.0,1.0);
  fHistMULTeta  = new TH1F("fHistMULTeta", " #eta distribution: multiply rec.",            100,-1.0,1.0);

  fHistRECphi   = new TH1F("fHistRECphi",  " #phi distribution: all reconstructed",        314, 0., 6.28);
  fHistMCphi    = new TH1F("fHistMCphi",   " #phi distribution: all MC",                   314, 0., 6.28);
  fHistMCNRphi  = new TH1F("fHistMCNRphi", " #phi distribution: all not-reconstructed MC", 314, 0., 6.28);
  fHistFAKEphi  = new TH1F("fHistFAKEphi", " #phi distribution: all Fake",                 314, 0., 6.28);
  fHistMULTphi  = new TH1F("fHistMULTphi", " #phi distribution: multipli rec.",            314, 0., 6.28);

  fHistRECHPTeta   = new TH1F("fHistRECHPTeta",  " #eta distribution: all reconstructed",        100,-1.0,1.0);
  fHistMCHPTeta    = new TH1F("fHistMCHPTeta",   " #eta distribution: all MC",                   100,-1.0,1.0);
  fHistMCNRHPTeta  = new TH1F("fHistMCNRHPTeta", " #eta distribution: all not-reconstructed MC", 100,-1.0,1.0);
  fHistFAKEHPTeta  = new TH1F("fHistFAKEHPTeta", " #eta distribution: all Fake",                 100,-1.0,1.0);

  fHistRECHPTphi   = new TH1F("fHistRECHPTphi",  " #phi distribution: all reconstructed",        314, 0., 6.28);
  fHistMCHPTphi    = new TH1F("fHistMCHPTphi",   " #phi distribution: all MC",                   314, 0., 6.28);
  fHistMCNRHPTphi  = new TH1F("fHistMCNRHPTphi", " #phi distribution: all not-reconstructed MC", 314, 0., 6.28);
  fHistFAKEHPTphi  = new TH1F("fHistFAKEHPTphi", " #phi distribution: all Fake",                 314, 0., 6.28);

  fHistRecMult     = new TH1F("fHistRecMult",  "Multiple reconstructed tracks", 50, 0., 50.);
  fHistNCluster    = new TH1F("fHistNCluster", "Number of clusters for suspicious tracks", 300, 0., 300.);

  // CKB

  fh1VtxEff = new TH1F("fh1VtxEff","VtxEff TPC",4,-0.5,3.5);
  fh1VtxEff->GetXaxis()->SetBinLabel(1,"NO TPC VTX");
  fh1VtxEff->GetXaxis()->SetBinLabel(2,"TPC VTX");
  fh1VtxEff->GetXaxis()->SetBinLabel(3,"NO SPD VTX");
  fh1VtxEff->GetXaxis()->SetBinLabel(4,"SPD VTX");

  TString labels[kTotalAC];
  labels[kNegA]="NegA";
  labels[kPosA]="PosA";
  labels[kNegC]="NegC";
  labels[kPosC]="PosC";

  for(int i = 0;i< kTotalAC;++i){
    fh2PhiPadRow[i] = new TH2F(Form("fh2PhiPadRow_%d",i),Form("Padrow vs phi AC %s;phi;padrow",labels[i].Data()),360,0.,360.,159,-0.5,158.5);
    fh2PhiLayer[i] = new TH2F(Form("fh2PhiLayer_%d",i),Form("layer vs phi AC %s;phi;layer",labels[i].Data()),360,0.,360.,6,-0.5,5.5);
  }
  for(int i = 0;i<2;++i){
    fh2MultSpdChips[i] = new TH2F(Form("fh2ChipsSpdMult_%d",i),"mult SPD vs chips;chips;Mult",201,-1.5,199.5,201,-1.5,199.5);
  }
  fh2VtxTpcSpd  = new TH2F("fh2VtxTpcSpd","SPD vtx vs TPC;tpc-z;spd-z",200,-20.,20.,200,-20,20.);

  TString labelsVtx[kTotalVtx];
  labelsVtx[kTPC] = "TPC";
  labelsVtx[kSPD] = "SPD";
  for(int i = 0; i < kTotalVtx; ++i){ 
    fh2VertexCorrelation[i] = new TH2F(Form("fh2VertexCorrelation_%d",i),Form("VertexCorrelation %s;MC z-vtx;ESD z-vtx",labelsVtx[i].Data()), 120, -30, 30, 120, -30, 30);
    fh2VertexCorrelationShift[i] = new TH2F(Form("fh2VertexCorrelationShift_%d",i), Form("VertexCorrelationShift %s;MC z-vtx;MC z-vtx - ESD z-vtx",labelsVtx[i].Data()), 120, -30, 30, 100, -1, 1); 
    fh1VertexShift[i] =  new TH1F(Form("fh1VertexShift_%d",i), Form("VertexShift %s;(MC z-vtx - ESD z-vtx);Entries",labelsVtx[i].Data()), 201, -2, 2);       
    fh1VertexShiftNorm[i] = new TH1F(Form("fh1VertexShiftNorm_%d",i), Form("VertexShiftNorm %s;(MC z-vtx - ESD z-vtx) / #sigma_{ESD z-vtx};Entries",labelsVtx[i].Data()), 200, -100, 100);      
  }

    
  fh2EtaCorrelation = new TH2F("fh2EtaCorrelation", "EtaCorrelation;MC #eta;ESD #eta", 120, -3, 3, 120, -3, 3);
  fh2EtaCorrelationShift = new TH2F("fh2EtaCorrelationShift", "EtaCorrelationShift;ESD #eta;MC #eta - ESD #eta", 120, -3, 3, 100, -0.1, 0.1);
  fh2PhiCorrelation = new TH2F("fh2PhiCorrelation", "PhiCorrelation;MC #phi;ESD #phi", 120, 0, 360, 120, 0, 360);
  fh2PhiCorrelationShift = new TH2F("fh2PhiCorrelationShift", "PhiCorrelationShift;ESD #phi;MC #phi - ESD #phi", 120, 0, 360, 100, -10, 10);
  fh2PtCorrelation = new TH2F("fh2PtCorrelation", "PtCorrelation;MC p_T;ESD p_T", 200, 0., 200., 200, 0, 200);
  fh2PtCorrelationShift = new TH2F("fh2PtCorrelationShift", "PtCorrelationShift;ESD p_T;MC p_T - ESD #p_T", 120, 0, 10, 100, -2, 2);
  
  fHistDeltaphiprimaries  =new TH1F("fHistDeltaphiprimaries", " #Delta#phi distribution: primaries",     314, -0.2,0.2);               
  fHistDeltaphisecondaries=new TH1F("fHistDeltaphisecondaries", "#Delta#phi distribution: secondaries", 314, -0.2,0.2);              
  fHistDeltaphireject     =new TH1F("fHistDeltaphireject", " #Delta#phi distribution: reject trackled",  314, -0.2,0.2); 
  
  fHistTPCITSdeltax =new TH1F("fHistTPCITSdeltax", "TPC-ITS matching:#Delta x", 100,-20,20); 
  fHistTPCITSdeltay =new TH1F("fHistTPCITSdeltay", "TPC-ITS matching:#Delta y", 100,-20,20); 
  fHistTPCITSdeltaz =new TH1F("fHistTPCITSdeltaz", "TPC-ITS matching:#Delta z", 100,-20,20); 
  fHistTPCITSdeltar =new TH1F("fHistTPCITSdeltar", "TPC-ITS matching:#Delta r", 100,0,20); 
  
  fHistTPCTRDdeltax   = new TH1F("fHistTPCTRDdeltax", "TPC-TRD matching:#Delta x", 400,-400, 400); 
  fHistTPCTRDdeltay   = new TH1F("fHistTPCTRDdeltay", "TPC-TRD matching:#Delta y", 400,-400, 400); 
  fHistTPCTRDdeltaz   = new TH1F("fHistTPCTRDdeltaz", "TPC-TRD matching:#Delta z", 400,-400, 400); 
  fHistTPCTRDdeltar   = new TH1F("fHistTPCTRDdeltar", "TPC-TRD matching:#Delta r", 100,0,20); 
  fHistTPCTRDdeltaphi = new TH1F("fHistTPCTRDdeltarphi", "TPC-TRD matching:#Delta #phi", 128,-3.14 , 3.14); 
  // Pulls
  fPtPullsPos = new TProfile("fPtPullsPos", "Pt Pulls for pos. primaries", 50 , 0., 10., -5., 5.,"S");
  fPtPullsNeg = new TProfile("fPtPullsNeg", "Pt Pulls for neg. primaries", 50 , 0., 10., -5., 5.,"S");
  fPtShiftPos = new TProfile("fPtShiftPos", "Pt Shift for pos. primaries", 50 , 0., 10., -0.2, 0.2,"S");
  fPtShiftNeg = new TProfile("fPtShiftNeg", "Pt Shift for neg. primaries", 50 , 0., 10., -0.2, 0.2,"S");

  fEtaMultiFluc = new TProfile("fEtaMultiFluc", "eta multiplicity fluctuations", 10, -2., 2., 0., 600., "S");
  fEtaMulti     = new TH1F("fEtaMulti",     "eta multiplicity fluctuations", 10, -2., 2.);
  fEtaMultiH    = new TH1F("fEtaMultiH",    "eta multiplicity fluctuations", 600,  0., 600.);

  fPhiMultiFluc = new TProfile("fPhiMultiFluc", "phi multiplicity fluctuations", 64, 0., 6.4, 0., 100., "S");
  fPhiMulti     = new TH1F("fPhiMulti",         "phi multiplicity fluctuations", 64, 0., 6.4);
  fPhiMultiH    = new TH1F("fPhiMultiH",        "phi multiplicity fluctuations", 100, 0., 100.);

  // SPD Vertex
  fVertexRvsZ  = new TH2F("fVertexRvsZ", "SPD Vertex R vs Z", 200, -1., 1., 200, -1., 1.);
  fVertexRvsZC = new TH2F("fVertexRvsZC", "SPD Vertex R vs Z", 200, -1., 1., 200, -1., 1.);
  
  TString charge[8]={"Deu", "AntiDeu", "Tri", "AntiTri", "He3", "AntiHe3", "He4", "AntiHe4"};
  for(Int_t i=0;i<8;i++){
    fHistRECptCharge[i]   = new TH1F(Form("fHistRECptCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all reconstructed",
				     100, 0., 20.);
    fHistMCptCharge[i]    = new TH1F(Form("fHistMCptCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all MC",
				     100, 0., 20.);
    fHistMCNRptCharge[i]  = new TH1F(Form("fHistMCNRptCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all not-reconstructed MC", 
				     100, 0., 20.);
    fHistFAKEptCharge[i]  = new TH1F( Form("fHistFAKEptCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all Fake",                 
				     100, 0., 20.);

    fHistRECetaCharge[i]   = new TH1F(Form("fHistRECetaCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all reconstructed",
				     100, 0., 20.);
    fHistMCetaCharge[i]    = new TH1F(Form("fHistMCetaCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all MC",
				     100, 0., 20.);
    fHistMCNRetaCharge[i]  = new TH1F(Form("fHistMCNRetaCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all not-reconstructed MC", 
				     100, 0., 20.);
    fHistFAKEetaCharge[i]  = new TH1F( Form("fHistFAKEetaCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all Fake",                 
				     100, 0., 20.);

    fHistRECphiCharge[i]   = new TH1F(Form("fHistRECphiCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all reconstructed",
				     100, 0., 20.);
    fHistMCphiCharge[i]    = new TH1F(Form("fHistMCphiCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all MC",
				     100, 0., 20.);
    fHistMCNRphiCharge[i]  = new TH1F(Form("fHistMCNRphiCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all not-reconstructed MC", 
				     100, 0., 20.);
    fHistFAKEphiCharge[i]  = new TH1F( Form("fHistFAKEphiCharge%s",charge[i].Data()), 
				     "p_{T} distribution: all Fake",                 
				     100, 0., 20.);

  }
  // BKC


  fHists->SetOwner();

  fHists->Add(fHistRECpt);
  fHists->Add(fHistMCpt);
  fHists->Add(fHistMCNRpt);
  fHists->Add(fHistFAKEpt);
  fHists->Add(fHistMULTpt);

  fHists->Add(fHistRECeta);
  fHists->Add(fHistMCeta);
  fHists->Add(fHistMCNReta);
  fHists->Add(fHistFAKEeta);
  fHists->Add(fHistMULTeta);

  fHists->Add(fHistRECphi);
  fHists->Add(fHistMCphi);
  fHists->Add(fHistMCNRphi);
  fHists->Add(fHistFAKEphi);
  fHists->Add(fHistMULTphi);

  fHists->Add(fHistRECHPTeta);
  fHists->Add(fHistMCHPTeta);
  fHists->Add(fHistMCNRHPTeta);
  fHists->Add(fHistFAKEHPTeta);

  fHists->Add(fHistRECHPTphi);
  fHists->Add(fHistMCHPTphi);
  fHists->Add(fHistMCNRHPTphi);
  fHists->Add(fHistFAKEHPTphi);

  // CKB

  fHists->Add(fh1VtxEff);
  for(int i = 0;i < kTotalAC;++i){
    fHists->Add(fh2PhiPadRow[i]);
    fHists->Add(fh2PhiLayer[i]);
  }
  for(int i = 0;i<2;++i){
    fHists->Add(fh2MultSpdChips[i]);
  }
  fHists->Add(fh2VtxTpcSpd);
    
  for(int i = 0; i < kTotalVtx; ++i){ 
    fHists->Add(fh2VertexCorrelation[i]);
    fHists->Add(fh2VertexCorrelationShift[i]);
    fHists->Add(fh1VertexShift[i]);      
    fHists->Add(fh1VertexShiftNorm[i]);
  }


  fHists->Add(fh2EtaCorrelation);
  fHists->Add(fh2EtaCorrelationShift);
  fHists->Add(fh2PhiCorrelation);
  fHists->Add(fh2PhiCorrelationShift);
  fHists->Add(fh2PtCorrelation);
  fHists->Add(fh2PtCorrelationShift);

  fHists->Add(fHistDeltaphiprimaries);              
  fHists->Add(fHistDeltaphisecondaries);
  fHists->Add(fHistDeltaphireject);  

  fHists->Add(fHistTPCITSdeltax);  
  fHists->Add(fHistTPCITSdeltay);  
  fHists->Add(fHistTPCITSdeltaz);  
  fHists->Add(fHistTPCITSdeltar);  
    
  fHists->Add(fHistTPCTRDdeltax);  
  fHists->Add(fHistTPCTRDdeltay);  
  fHists->Add(fHistTPCTRDdeltaz);  
  fHists->Add(fHistTPCTRDdeltar);  
  fHists->Add(fHistTPCTRDdeltaphi);  

  fHists->Add(fHistRecMult);
  fHists->Add(fHistNCluster);
 

  fHists->Add(fPtPullsPos);
  fHists->Add(fPtPullsNeg);
  fHists->Add(fPtShiftPos);
  fHists->Add(fPtShiftNeg);

  fHists->Add(fVertexRvsZ);
  fHists->Add(fVertexRvsZC);
  fHists->Add(fEtaMultiFluc);
  fHists->Add(fEtaMulti);
  fHists->Add(fEtaMultiH);
  fHists->Add(fPhiMultiFluc);
  fHists->Add(fPhiMulti);
  fHists->Add(fPhiMultiH);

  for(Int_t i=0;i<8;i++){
    fHists->Add(fHistRECptCharge[i]);
    fHists->Add(fHistMCptCharge[i]);
    fHists->Add(fHistMCNRptCharge[i]);
    fHists->Add(fHistFAKEptCharge[i]);

    fHists->Add(fHistRECetaCharge[i]);
    fHists->Add(fHistMCetaCharge[i]);
    fHists->Add(fHistMCNRetaCharge[i]);
    fHists->Add(fHistFAKEetaCharge[i]);

    fHists->Add(fHistRECphiCharge[i]);
    fHists->Add(fHistMCphiCharge[i]);
    fHists->Add(fHistMCNRphiCharge[i]);
    fHists->Add(fHistFAKEphiCharge[i]);
  }

  for (Int_t i=0; i<fHists->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHists->At(i));
    if (h1){
      h1->Sumw2();
    }
  }
  // BKC


  // Post output data.
  PostData(1, fHists);

}

//________________________________________________________________________
void AliAnalysisTaskEfficiency::UserExec(Option_t *) 
{
  // Event loop
  // Analysis of MC true particles and reconstructed tracks
  // Different track types (Global, TPC, ITS) can be selected

  if (!fInputEvent) {
    Printf("ERROR: fESD not available");
    return;
  }
  static Int_t nfc = 0;

  //AliESDInputHandler* esdI = (AliESDInputHandler*) 
  //(AliAnalysisManager::GetAnalysisManager())->GetInputEventHandler();
  //AliESDEvent* hltEvent = esdI->GetHLTEvent();

  AliStack* stack = MCEvent()->Stack();
  AliESDEvent* esdE = (AliESDEvent*) fInputEvent;
    
  // CKB
  // Fetch the SPD vertex:
  Double_t vtxSPD[3];
  const AliESDVertex* spdVertex = esdE->GetPrimaryVertexSPD();
  if (spdVertex) {
    if(spdVertex->GetNContributors()<=0){
      spdVertex = 0;
    }
    else{
      spdVertex->GetXYZ(vtxSPD);
    }
  }
  
  // Fetch the TPC vertex
  Double_t vtxTPC[3];
  const AliESDVertex* tpcVertex = esdE->GetPrimaryVertexTPC();
  if (tpcVertex) {
    if(tpcVertex->GetNContributors()<=0){
      tpcVertex = 0;
    }
    else{
      tpcVertex->GetXYZ(vtxTPC);
    }
  }
  // SPD Vertex
  if (spdVertex) {
    Double_t x,y,z;
    x = spdVertex->GetX() + 0.07;
    y = spdVertex->GetY() - 0.25;
    z = spdVertex->GetZ();
    if (TMath::Abs(z) < 10.) {
      fVertexRvsZ->Fill(x,y);
      if (TMath::Sqrt(x*x + y*y) > 5. * 0.028) 
	fVertexRvsZC->Fill(x,y);
    }
  }
  // BKC
  //Printf("%s:%d %5d",(char*)__FILE__,__LINE__, esdE->GetNumberOfTracks());
  TArrayI labels(esdE->GetNumberOfTracks());
  Int_t igood = 0;
  // Track loop to fill a pT spectrum
  //printf("ESD track loop \n");

  AliESDtrack *tpcP = 0x0;

  for (Int_t iTracks = 0; iTracks < esdE->GetNumberOfTracks(); iTracks++) {

    //prevent mem leak for TPConly track
    if(fTrackType==2&&tpcP){
      delete tpcP;
      tpcP = 0;
    }

    AliVParticle *track = esdE->GetTrack(iTracks);
    AliESDtrack *esdtrack =  dynamic_cast<AliESDtrack*>(track);
    esdtrack->PropagateToDCA(esdE->GetPrimaryVertex(),
			     esdE->GetMagneticField(), 10000.);

    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    //__________
    // run Task for global tracks or ITS tracks or TPC tracks
    const AliExternalTrackParam *tpcPin = 0x0;
    //Double_t phiIn=0.;
    Float_t phi = 0.;
    Float_t eta = 0.;
    Float_t pt  = 0.;
    Int_t indexAC = GetIndexAC(esdtrack);  

    //select in the steering macro, which track type should be analysed. 
    // Four track types are selectable:
    // 0. Global tracks
    // 1. ITS tracks (SA or Pure SA)
    // 2. TPC only tracks
    // Track selection copied from PWG1/AliAnalysisTaskQAsym

    if(fTrackType==0){
      //Fill all histograms with global tracks
      tpcP = esdtrack;
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      phi = tpcP->Phi();
      eta = tpcP->Eta();
      pt  = tpcP->Pt();
    }
    else if(fTrackType==1){
      //Fill all histograms with ITS tracks
      tpcP = esdtrack;
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      phi = tpcP->Phi();
      eta = tpcP->Eta();
      pt  = tpcP->Pt();
    }
    else if(fTrackType==2){     
      //Fill all histograms with TPC track information
      tpcPin = esdtrack->GetInnerParam();
      if (!tpcPin) continue;
      tpcP = AliESDtrackCuts::
	GetTPCOnlyTrack(dynamic_cast<AliESDEvent*>(esdE),
			esdtrack->GetID());
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      if(tpcP->GetNcls(1)>160)continue;//jacek's track cut
      if(tpcP->GetConstrainedChi2TPC()<0)continue; // jacek's track cut
      phi=tpcPin->Phi();
      eta=tpcPin->Eta();
      pt=tpcPin->Pt();
    }
    else{
      Printf("ERROR: wrong track type \n");
      continue;
    }
    //___________
    //


    labels.AddAt(-1, iTracks);
      

      // Selected
      if (TMath::Abs(eta) > 0.9) {
      } else {

	//Int_t charge=track->Charge();
	// Within acceptance
	Int_t ind = TMath::Abs(esdtrack->GetLabel());
	AliMCParticle* mcPtest = (AliMCParticle*) MCEvent()->GetTrack(ind);
	if (stack->IsPhysicalPrimary(mcPtest->Label())){
	  if(TMath::Abs(mcPtest->PdgCode())>=1000020030){//all helium (+-1000020030,+-1000020040)
	    pt*=2;// reconstruction takes charge=1 -> for particles with charge=2, pT,rec need to be corrected
	  }
	  
	  Int_t index=ConvertePDG(mcPtest->PdgCode());
	  if(fDebug>1)Printf("PDG=%d, index=%d", mcPtest->PdgCode(), index);
	  //fill tracks comming from d,t,3He, 4He
	  if(index<8){
	    fHistRECptCharge [index] ->Fill(pt);
	    fHistRECetaCharge[index]->Fill(eta);
	    fHistRECphiCharge[index]->Fill(phi);
	  }

	}
	fHistRECpt ->Fill(pt);
	fHistRECeta->Fill(eta);
	fHistRECphi->Fill(phi);

	if (pt > 2.) {
	  fHistRECHPTeta->Fill(eta);
	  fHistRECHPTphi->Fill(phi);
	}
	  
	//	Int_t ind = TMath::Abs(esdtrack->GetLabel());
	AliMCParticle* mcP = (AliMCParticle*) MCEvent()->GetTrack(ind);
	if (!(stack->IsPhysicalPrimary(mcP->Label()))) {
	  fHistFAKEpt ->Fill(pt);
	  fHistFAKEeta->Fill(eta);		  
	  fHistFAKEphi->Fill(phi);
	  
	  if (pt > 2.) {
	    fHistFAKEHPTeta->Fill(eta);
	    fHistFAKEHPTphi->Fill(phi);
	  }
	    
	}
	  
	labels.AddAt(TMath::Abs(esdtrack->GetLabel()), iTracks);
	igood++;
	  
	  
	Float_t phiDeg = TMath::RadToDeg()*track->Phi();
	  
	//TPC-ITS matching
	//TPC-TRD matching  
	      
	if (tpcP){
	  Int_t labeltpcits = esdtrack->GetLabel();
	  AliExternalTrackParam * tpcPCopy = new AliExternalTrackParam(*tpcP);
	  Double_t xk = 43.6;  // cm
	  Double_t bz = 5.0;   // kG
	  if(tpcPCopy->PropagateTo(xk,bz)){
	    Double_t alpha=tpcPCopy->GetAlpha();
	    // if(tpcPCopy->Rotate(0.)){ 
	    Float_t xtpc=tpcPCopy->GetX();
	    Float_t ytpc=tpcPCopy->GetY();
	    Float_t ztpc=tpcPCopy->GetZ();
	    Float_t xpr,ypr ; 
	    xpr=xtpc*TMath::Cos(alpha)-ytpc*TMath::Sin(alpha);
	    ypr=xtpc*TMath::Sin(alpha)+ytpc*TMath::Cos(alpha);
		
	    AliMCParticle* mcPart = (AliMCParticle*) MCEvent()->GetTrack(abs(labeltpcits));
	    Int_t ntref = mcPart->GetNumberOfTrackReferences();
		 
	    for (Int_t k = 0; k < ntref; k++) {
	      AliTrackReference* ref = mcPart->GetTrackReference(k);
	      Float_t xhits = ref->X();
	      Float_t yhits = ref->Y();
	      Float_t radio = TMath::Sqrt(xhits*xhits+yhits*yhits);
		
	      if(ref->DetectorId() == 0 && (radio > 42.8 && radio < 43.7)) {
		  
		Float_t xits = ref->X();
		Float_t yits = ref->Y();
		Float_t zits = ref->Z();
		  
		Float_t difx=(xits-xpr);
		fHistTPCITSdeltax->Fill(difx);
		Float_t dify=(yits-ypr);
		fHistTPCITSdeltay->Fill(dify);
		Float_t difz=(zits-ztpc);
		fHistTPCITSdeltaz->Fill(difz);
		Float_t deltar = TMath::Sqrt(difx * difx + dify * dify + difz * difz);
		fHistTPCITSdeltar->Fill(deltar);
		break;
	      }
	      // }
	    } // trackrefs
	  } 
	} //ITS-TPC maching
	  //TPC-TRD maching 
	       
	const AliExternalTrackParam *trd = esdtrack->GetOuterParam();
	  
	if (trd){Int_t labeltpctrd = track->GetLabel();
	  
	  AliExternalTrackParam * trdCopy = new AliExternalTrackParam(*trd);
	  Double_t xktrd=296.0;
	  Double_t bztrd=5.0;
	  if(trdCopy->PropagateTo(xktrd,bztrd)){
	    Float_t xtpc2=trdCopy->GetX();
	    Float_t ytpc2=trdCopy->GetY();
	    Float_t ztpc2=trdCopy->GetZ();
	    Double_t alpha=trdCopy->GetAlpha();
	    Float_t xpr,ypr ; 
	    xpr=xtpc2*TMath::Cos(alpha)-ytpc2*TMath::Sin(alpha);
	    ypr=xtpc2*TMath::Sin(alpha)+ytpc2*TMath::Cos(alpha);
	    AliMCParticle* mcPart = (AliMCParticle*) MCEvent()->GetTrack(abs(labeltpctrd));
	    Int_t ntreftrd = mcPart->GetNumberOfTrackReferences(); 
	    
	    for (Int_t k = 0; k < ntreftrd; k++) {
	      AliTrackReference* ref = mcPart->GetTrackReference(k);
	      if(ref->DetectorId() == 3 && ref->Label() == labeltpctrd){
		
		Float_t xtrd = ref->X();
		Float_t ytrd = ref->Y();
		Float_t ztrd = ref->Z();
		
		Float_t difx=(xtrd-xpr);
		fHistTPCTRDdeltax->Fill(difx);
		Float_t dify=(ytrd-ypr);
		fHistTPCTRDdeltay->Fill(dify);
		Float_t difz=(ztrd-ztpc2);
		fHistTPCTRDdeltaz->Fill(difz);
		Float_t deltar = TMath::Sqrt(difx * difx + dify * dify + difz * difz);
		fHistTPCTRDdeltar->Fill(deltar);
		Float_t phi_tpc = TMath::ATan2(ypr, xpr);
		Float_t phi_trd = TMath::ATan2(ytrd, xtrd);
		Float_t dphi = phi_trd - phi_tpc;
		if (dphi >   TMath::Pi()) dphi -= 2. * TMath::Pi();
  		if (dphi < - TMath::Pi()) dphi += 2. * TMath::Pi();

		fHistTPCTRDdeltaphi->Fill(dphi);
		break;
	      }
	      
	    }
	  }

	}//TRD-TPC maching

	// CKB
	if(pt>2.&&fFieldOn){
	  TBits bits(esdtrack->GetTPCClusterMap());
	  for(unsigned int ip = 0;ip < bits.GetNbits();++ip){
	    if(bits[ip]){
	      fh2PhiPadRow[indexAC]->Fill(phiDeg,ip);
	    }
	  }
	  for(int i = 0;i < 6;++i){
	    if(esdtrack->HasPointOnITSLayer(i))fh2PhiLayer[indexAC]->Fill(phiDeg,i);
	  }
	}
	else if(!fFieldOn){ // field not on all momenta are set to MPV 
	  TBits bits(esdtrack->GetTPCClusterMap());
	  for(unsigned int ip = 0;ip < bits.GetNbits();++ip){
	    if(bits[ip]){
	      fh2PhiPadRow[indexAC]->Fill(phiDeg,ip);
	    }
	    for(int i = 0;i < 6;++i){
	      if(esdtrack->HasPointOnITSLayer(i))fh2PhiLayer[indexAC]->Fill(phiDeg,i);
	    }
	  }      

	}
	// Fill the correlation
	if(MCEvent()){
	  TParticle *part = MCEvent()->Stack()->Particle(TMath::Abs(esdtrack->GetLabel()));
	  if(part){
	    Float_t mcEta = part->Eta();
	    Float_t mcPhi = TMath::RadToDeg()*part->Phi();
	    Float_t mcPt  = part->Pt();
	    //		   if (pt - mcPt > 20.) {
	    fh2EtaCorrelation->Fill(mcEta,eta);
	    fh2EtaCorrelationShift->Fill(eta,mcEta-eta);
	    fh2PhiCorrelation->Fill(mcPhi,phiDeg);
	    fh2PhiCorrelationShift->Fill(phiDeg,mcPhi-phiDeg);
	    fh2PtCorrelation->Fill(mcPt,pt);
	    fh2PtCorrelationShift->Fill(pt,mcPt-pt);
	    //}
	    Double_t sigmaPt = TMath::Sqrt(esdtrack->GetSigma1Pt2());
	    if (track->Charge() > 0.) {
	      fPtPullsPos->Fill(mcPt, (1./pt - 1./mcPt) / sigmaPt); 
	      fPtShiftPos->Fill(mcPt, (1./pt - 1./mcPt));
	    } else {
	      fPtPullsNeg->Fill(mcPt, (1./pt - 1./mcPt) / sigmaPt); 
	      fPtShiftNeg->Fill(mcPt, (1./pt - 1./mcPt));  
	    }
	  }
	  // BKC
	}
      }
   
  } //track loop 

  //prevent mem leak for TPConly track
  if(fTrackType==2&&tpcP){
    delete tpcP;
    tpcP = 0;
  }

  //Printf("%s:%d",(char*)__FILE__,__LINE__);
  // CKB
  // Vertex eff
  if(tpcVertex){
    fh1VtxEff->Fill("TPC VTX",1);
  }
  else{
    fh1VtxEff->Fill("NO TPC VTX",1);
  }

  if(spdVertex){
    fh1VtxEff->Fill("SPD VTX",1);
  }
  else{
    fh1VtxEff->Fill("NO SPD VTX",1);
  }


  // Vertex correlation SPD vs. TPC
  if(tpcVertex&&spdVertex){
    fh2VtxTpcSpd->Fill(vtxTPC[2],vtxSPD[2]);
  }
  //  Printf("%s:%d",(char*)__FILE__,__LINE__);
  // Multiplicity checks in the SPD
  const AliMultiplicity *spdMult = esdE->GetMultiplicity();
  // Multiplicity Analysis
  Int_t nt = spdMult->GetNumberOfTracklets();
  nfc++;
  if (nfc == 520) {
    nfc = 0;
    for (Int_t ib = 1; ib <= 10; ib++) {
      Int_t ic = Int_t(fEtaMulti->GetBinContent(ib));
      Float_t eta = -1.8 + Float_t(ib-1) * 0.4;
      fEtaMultiFluc->Fill(eta, Float_t(ic));
      if (ib == 5) fEtaMultiH->Fill(Float_t(ic));
    }

    for (Int_t ib = 1; ib <= 64; ib++) {
      Int_t ic = Int_t(fPhiMulti->GetBinContent(ib));
      Float_t phi = 0.05 + Float_t(ib-1) * 0.1;
      fPhiMultiFluc->Fill(phi, Float_t(ic));
      if (ib  == 2) fPhiMultiH->Fill(Float_t(ic));
    }

    fEtaMulti->Reset();
    fPhiMulti->Reset();    
  }
	
  for (Int_t j = 0; j < nt; j++) {
    fEtaMulti->Fill(spdMult->GetEta(j));
    fPhiMulti->Fill(spdMult->GetPhi(j));
  }

  Short_t chips0 = spdMult->GetNumberOfFiredChips(0);
  Short_t chips1 = spdMult->GetNumberOfFiredChips(1);
  //Printf("%s:%d",(char*)__FILE__,__LINE__);
  Int_t inputCount = 0;
  
  if(spdVertex){
    // get multiplicity from ITS tracklets
    for (Int_t i=0; i<spdMult->GetNumberOfTracklets(); ++i){
      Float_t deltaPhi = spdMult->GetDeltaPhi(i);
      // prevent values to be shifted by 2 Pi()
      if (deltaPhi < -TMath::Pi())
	deltaPhi += TMath::Pi() * 2;
      if (deltaPhi > TMath::Pi())
	deltaPhi -= TMath::Pi() * 2;      
      if (TMath::Abs(deltaPhi) > 1)
	printf("WARNING: Very high Delta Phi: %d %f %f %f\n", i, spdMult->GetTheta(i), spdMult->GetPhi(i), deltaPhi);
      //trackled Deltaphi
      Int_t label1=spdMult->GetLabel(i,0);
      Int_t label2=spdMult->GetLabel(i,1);
      if ( label1==label2){
	if(stack->IsPhysicalPrimary(label1) == kTRUE)
	  fHistDeltaphiprimaries->Fill(deltaPhi);
     	else
	  fHistDeltaphisecondaries->Fill(deltaPhi);
      }
      else{
      	fHistDeltaphireject->Fill(deltaPhi);
      }
      ++inputCount;
    }
    // Printf("%s:%d",(char*)__FILE__,__LINE__);
    fh2MultSpdChips[0]->Fill(chips0,inputCount);
    fh2MultSpdChips[1]->Fill(chips1,inputCount);
    //Printf("%s:%d",(char*)__FILE__,__LINE__);
  }
  // BKC

  // MC analysis
  Int_t nmc = MCEvent()->GetNumberOfTracks();
  AliGenEventHeader*       header = MCEvent()->GenEventHeader();
  // some test
  AliGenCocktailEventHeader* cH = dynamic_cast<AliGenCocktailEventHeader*> (header);
  AliGenPythiaEventHeader* pH;
  if (cH == 0) 
    {
      header->Dump();
      pH = dynamic_cast<AliGenPythiaEventHeader*>(header);
    } else {
    TObject* entry = cH->GetHeaders()->FindObject("Pythia");
    pH = dynamic_cast<AliGenPythiaEventHeader*>(entry);
  }
  Int_t iproc = -1;

  if (pH) iproc = pH->ProcessType();

  TArrayF mcV;
  header->PrimaryVertex(mcV);
  Float_t vzMC = mcV[2];

  // CKB
  // Fill the correlation with MC
  if(spdVertex&&MCEvent()){
    Float_t diff = vzMC-vtxSPD[2];
    fh2VertexCorrelation[kSPD]->Fill(vzMC,vtxSPD[2]);
    fh2VertexCorrelationShift[kSPD]->Fill(vzMC,diff);
    fh1VertexShift[kSPD]->Fill(diff);
    if(spdVertex->GetZRes()>0)fh1VertexShiftNorm[kSPD]->Fill(diff/spdVertex->GetZRes());
  }
  if(tpcVertex&&MCEvent()){
    Float_t diff = vzMC-vtxTPC[2];
    fh2VertexCorrelation[kTPC]->Fill(vzMC,vtxTPC[2]);
    fh2VertexCorrelationShift[kTPC]->Fill(vzMC,diff);
    fh1VertexShift[kTPC]->Fill(diff);
    if(tpcVertex->GetZRes()>0)fh1VertexShiftNorm[kTPC]->Fill(diff/tpcVertex->GetZRes());
  }
  // BKC


  // MC event loop
  //printf("MC particle loop \n");
  for (Int_t i = 0; i < nmc; i++) {
    AliMCParticle* mcP = (AliMCParticle*) MCEvent()->GetTrack(i);
    //printf("MC particle loop %5d \n", i);
    // Primaries only
    if (!(stack->IsPhysicalPrimary(mcP->Label()))) continue;
    if (mcP->Particle()->GetPDG()->Charge() == 0) continue;
    // Int_t charge=Int_t(mcP->Particle()->GetPDG()->Charge());
    Float_t phi = mcP->Phi();
    Float_t eta = mcP->Eta();
    Float_t pt  = mcP->Pt();
    if (TMath::Abs(eta) > 0.9) continue;
    Int_t ntr = mcP->GetNumberOfTrackReferences();
    Int_t nITS = 0;
    Int_t nTPC = 0;
    Int_t nFRA = 0;
    Float_t  x = 0.;
    Float_t  y = 0.;
      
    for (Int_t j = 0; j < ntr; j++) {
         
      AliTrackReference* ref = mcP->GetTrackReference(j);
	
      if (ref->DetectorId() == 0) nITS++;
      if (ref->DetectorId() == 1) nTPC++;
      if (ref->DetectorId() == 2) nFRA++;
      if (nTPC == 1) {
	x = ref->X();
	y = ref->Y();
	break;
      }
    }

    fHistMCpt ->Fill(pt);
    fHistMCeta->Fill(eta);
    fHistMCphi->Fill(phi);
    Int_t index=ConvertePDG(mcP->PdgCode());
    if(index<8){
      fHistMCptCharge [index] ->Fill(pt);
      fHistMCetaCharge[index]->Fill(eta);
      fHistMCphiCharge[index]->Fill(phi);
    }


    if (pt > 2.) {
      fHistMCHPTeta->Fill(eta);
      fHistMCHPTphi->Fill(phi);
    }

    Bool_t reco = kFALSE;
    Int_t multR =  0;
    Int_t jold  = -1;
      
    for (Int_t j = 0; j < esdE->GetNumberOfTracks(); j++) {
      if (i == labels.At(j)) {
	reco = kTRUE;
	multR++;
	//AliESDtrack* test = esdE->GetTrack(j);
	if (multR > 1) {
	  Int_t nclus = 0;
	  AliESDtrack* track = esdE->GetTrack(jold);
	  nclus = track->GetTPCclusters(0);		 
	  fHistNCluster->Fill(nclus);
	  
		  
	  track = esdE->GetTrack(j);
	  nclus = track->GetTPCclusters(0);
	  fHistNCluster->Fill(nclus);
	  

	}
	jold = j;
      }
    }

    fHistRecMult->Fill(Float_t(multR), 1.);
    if (multR == 0) {
      fHistMCNRpt ->Fill(pt);
      fHistMCNReta->Fill(eta);
      fHistMCNRphi->Fill(phi);
      if(index<8){
	fHistMCNRptCharge [index] ->Fill(pt);
	fHistMCNRetaCharge[index]->Fill(eta);
	fHistMCNRphiCharge[index]->Fill(phi);
      }
     

      if (pt > 2.) {
	fHistMCNRHPTeta->Fill(eta);
	fHistMCNRHPTphi->Fill(phi);
      }
    }
    else if (multR > 1)
      {
	fHistMULTpt ->Fill(pt);
	fHistMULTeta->Fill(eta);
	fHistMULTphi->Fill(phi);
      }	
  } // MC particle loop
  //printf("End of MC particle loop \n");
  
  
}      

//________________________________________________________________________
void AliAnalysisTaskEfficiency::Terminate(Option_t *) 
{

}  


Bool_t AliAnalysisTaskEfficiency::SelectJouri(AliESDtrack* track) 
{

  Bool_t selected = kTRUE;
  AliESDEvent* esdE = (AliESDEvent*) fInputEvent;    
  // > 50 Clusters
  if (track->GetTPCclusters(0) < 50) selected = kFALSE;

  const AliESDVertex *vtx=esdE->GetPrimaryVertexSPD();
  if (!vtx->GetStatus()) selected = kFALSE;
    
  Double_t zv=vtx->GetZv();

    
  const AliExternalTrackParam *ct=track->GetTPCInnerParam();
  if (!ct)  selected = kFALSE;
    
  Double_t xyz[3];
  ct->GetXYZ(xyz);
  Float_t rv = TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
  if (rv > 3.0)                     selected = kFALSE;
  if (TMath::Abs(xyz[2] - zv)>2.5) selected = kFALSE;
  return (selected);
    
}

Int_t AliAnalysisTaskEfficiency::GetIndexAC(AliESDtrack *track){
  if(!track)return -1;

  // Crude selection for A and C side
  // just with eta
  if (track->Eta() > 0) { 
    if (track->Charge() > 0) 
      return kPosA;
    else
      return kNegA;	
  }
  else { // C side
    if (track->Charge() > 0) 
      return kPosC;
    else
      return kNegC;	
  }
  
  return -1;

}


Int_t AliAnalysisTaskEfficiency::ConvertePDG(Int_t pdg)
{

  // function converts the pdg values for nuclei (d, t, 3He, 
  // 4He and their antiparticles) into indices for histo-array
  // filled in UserExec
  
  if      ( pdg ==  1000010020 )return 0;
  else if ( pdg == -1000010020 )return 1;
  else if ( pdg ==  1000010030 )return 2;
  else if ( pdg == -1000010030 )return 3;
  else if ( pdg ==  1000020030 )return 4;
  else if ( pdg == -1000020030 )return 5;
  else if ( pdg ==  1000020040 )return 6;
  else if ( pdg == -1000020040 )return 7;
  else return 9;
  
  
}
