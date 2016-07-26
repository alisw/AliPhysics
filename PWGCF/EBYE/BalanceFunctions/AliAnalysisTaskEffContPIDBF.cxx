#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TParticle.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h" 
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliCentrality.h"
#include "AliGenEventHeader.h"

#include "AliLog.h"
#include "AliAnalysisTaskEffContPIDBF.h"

// ---------------------------------------------------------------------
//
// Task for calculating the efficiency of the Balance Function 
// for single particles and pairs
// 
// Modified By Noor Alam (VECC ,Kolkata) sk.noor.alam@cern.ch
//[ Special thanks to Michael Weber ]
// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEffContPIDBF)

//________________________________________________________________________
AliAnalysisTaskEffContPIDBF::AliAnalysisTaskEffContPIDBF(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAOD(0),
    fPIDType(kNSigmaTPCTOF),
    fNSigmaPID(3.0),
    fArrayMC(0), 
    fQAList(0), 
    fOutputList(0), 
    fHistEventStats(0), 
    fHistCentrality(0),
    fHistVz(0), 
    fHistNSigmaTPCvsPtbeforePID(0),
    fHistNSigmaTPCvsPtafterPID(0),  
    fHistTruthPionPlus(0),
    fHistTruthKaonPlus(0),
    fHistTruthProtonPlus(0),
    fHistTruthPionMinus(0),
    fHistTruthKaonMinus(0),
    fHistTruthProtonMinus(0),
    fHistTruthPion(0),
    fHistTruthKaon(0),
    fHistTruthProton(0),
    HistMCTruthPtAll(0),
    HistMCTruthEtaAll(0),
    HistMCTruthPhiAll(0),
    HistMCTruthPtPion(0),
    HistMCTruthEtaPion(0),
    HistMCTruthPhiPion(0),
    HistMCTruthPtKaon(0),
    HistMCTruthEtaKaon(0),
    HistMCTruthPhiKaon(0),
    HistMCTruthPtProton(0),
    HistMCTruthEtaProton(0),
    HistMCTruthPhiProton(0),
    HistPionContaminationInPt(0),
    HistPionPlusContaminationInPt(0),
    HistPionMinusContaminationInPt(0),
    Hist3dPionContamination(0),
    Hist3dPionPlusContamination(0),
    Hist3dPionMinusContamination(0),
    HistKaonContaminationInPt(0),
    HistKaonPlusContaminationInPt(0),
    HistKaonMinusContaminationInPt(0),
    Hist3dKaonContamination(0),
    Hist3dKaonPlusContamination(0),
    Hist3dKaonMinusContamination(0),
    HistProtonContaminationInPt(0),
    HistProtonPlusContaminationInPt(0),
    HistProtonMinusContaminationInPt(0),
    Hist3dProtonContamination(0),
    Hist3dProtonPlusContamination(0),
    Hist3dProtonMinusContamination(0), 
    HistPionPurityInPt(0),
    HistPionPlusPurityInPt(0),
    HistPionMinusPurityInPt(0),
    Hist3dPionPurity(0),
    Hist3dPionPlusPurity(0),
    Hist3dPionMinusPurity(0),
    HistKaonPurityInPt(0),
    HistKaonPlusPurityInPt(0),
    HistKaonMinusPurityInPt(0),
    Hist3dKaonPurity(0),
    Hist3dKaonPlusPurity(0),
    Hist3dKaonMinusPurity(0),
    HistProtonPurityInPt(0),
    HistProtonPlusPurityInPt(0),
    HistProtonMinusPurityInPt(0),
    Hist3dProtonPurity(0),
    Hist3dProtonPlusPurity(0),
    Hist3dProtonMinusPurity(0),
    fHistSigmaTPCVsTOFPionForPionAfterCut(0),
    fHistSigmaTPCVsTOFProtonForPionAfterCut(0),
    fHistSigmaTPCVsTOFKaonForPionAfterCut(0),
    fHistSigmaTPCVsTOFPionForKaonAfterCut(0),
    fHistSigmaTPCVsTOFProtonForKaonAfterCut(0),
    fHistSigmaTPCVsTOFKaonForKaonAfterCut(0),
    fHistSigmaTPCVsTOFPionForProtonAfterCut(0),
    fHistSigmaTPCVsTOFProtonForProtonAfterCut(0),
    fHistSigmaTPCVsTOFKaonForProtonAfterCut(0),
    h1PionAfterCut(0),
    h1KaonAfterCut(0),
    h1ProtonAfterCut(0),
    h1PionAsNonPion(0),
    h1KaonAsNonKaon(0),
    h1ProtonAsNonProton(0),
    h1PionAsKaon(0),
    h1PionAsProton(0),
    h1KaonAsPion(0),
    h1KaonAsProton(0),
    h1ProtonAsPion(0),
    h1ProtonAsKaon(0),
    fHistMCRecoPionPlus(0),
    fHistMCRecoKaonPlus(0),
    fHistMCRecoProtonPlus(0),
    fHistMCRecoPionMinus(0),
    fHistMCRecoKaonMinus(0),
    fHistMCRecoProtonMinus(0),
    fHistMCRecoPion(0),
    fHistMCRecoKaon(0),
    fHistMCRecoProton(0),
    fHistMCRecoPionAsKaon(0),
    fHistMCRecoPionAsProton(0),
    fHistMCRecoProtonAsKaon(0),
    fHistMCRecoProtonAsPion(0),
    fHistMCRecoKaonAsPion(0),
    fHistMCRecoKaonAsProton(0),
    fUseCentrality(kFALSE),
    fCentralityEstimator("V0M"), 
    fCentralityPercentileMin(0.0), 
    fCentralityPercentileMax(5.0), 
    fInjectedSignals(kFALSE),
    fPIDResponse(0),
    fElectronRejection(kFALSE),
    fElectronOnlyRejection(kFALSE),
    fElectronRejectionNSigma(-1.),
    fElectronRejectionMinPt(0.),
    fElectronRejectionMaxPt(1000.),
    fVxMax(3.0), 
    fVyMax(3.0),
    fVzMax(10.), 
    fAODTrackCutBit(128),
    fMinNumberOfTPCClusters(80),
    fMaxChi2PerTPCCluster(4.0),
    fMaxDCAxy(3.0),
    fMaxDCAz(3.0),
    fMinPt(0.0),
    fMaxPt(20.0),
    fPtTPCMax(0.6),
    fMinEta(-0.8), 
    fMaxEta(0.8),
    fEtaRangeMin(0.0), 
    fEtaRangeMax(1.6), 
    fPtRangeMin(0.0), 
    fPtRangeMax(20.0), 
    fEtaBin(100), //=100 (BF) 16
    fdEtaBin(64), //=64 (BF)  16
    fPtBin(100), //=100 (BF)  36
    fHistdEdxTPC(0),
    fHistBetaTOF(0),
    fParticleType_(kPion)
{ 

for(Int_t ipart=0;ipart<3;ipart++)
    for(Int_t ipid=0;ipid<3;ipid++)
      fnsigmas[ipart][ipid]=999.;


  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskEffContPIDBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  fQAList = new TList();
  fQAList->SetName("QAList");
  fQAList->SetOwner();
  
  fOutputList = new TList();
  fOutputList->SetName("OutputList");
  fOutputList->SetOwner();
  
  //Event stats.
  TString gCutName[4] = {"Total","Offline trigger",
                         "Vertex","Analyzed"};
  fHistEventStats = new TH1F("fHistEventStats",
                             "Event statistics;;N_{events}",
                             4,0.5,4.5);
  for(Int_t i = 1; i <= 4; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fQAList->Add(fHistEventStats);

  //====================================================//
  Int_t ptBin = 40;
  Int_t etaBin = 32;
  Int_t phiBin = 100;

  Double_t nArrayPt[41]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2};
  Double_t nArrayEta[33]={-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8}; 

  Double_t nArrayPhi[phiBin+1];
  for(Int_t iBin = 0; iBin <= phiBin; iBin++) 
    nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;

  //AOD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",1001,-0.5,100.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)

  //Vz addition+++++++++++++++++++++++++++++
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fQAList->Add(fHistVz);

  //Electron cuts -> PID QA
  fHistNSigmaTPCvsPtbeforePID = new TH2F ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore",200, 0, 20, 200, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtbeforePID);

  fHistNSigmaTPCvsPtafterPID = new TH2F ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter",200, 0, 20, 200, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtafterPID);


    HistMCTruthPtAll= new TH1F("HistMCTruthPtAll","Pt Distribution Of all  charged Partilces",ptBin,nArrayPt);
    HistMCTruthEtaAll= new TH1F("HistMCTruthEtaAll","Eta Distribution Of all  charged Partilces",etaBin,nArrayEta);
    HistMCTruthPhiAll= new TH1F("HistMCTruthPhiAll","Phi Distribution Of all  charged Partilces",phiBin,nArrayPhi);

    HistMCTruthPtPion= new TH1F("HistMCTruthPtPion","Pt Distribution Of Pion  Partilces",ptBin,nArrayPt);
    HistMCTruthEtaPion= new TH1F("HistMCTruthEtaPion","Eta Distribution Of Pion   Partilces",etaBin,nArrayEta);
    HistMCTruthPhiPion= new TH1F("HistMCTruthPhiPion","Phi Distribution Of Pion   Partilces",phiBin,nArrayPhi);

    HistMCTruthPtKaon= new TH1F("HistMCTruthPtKaon","Pt Distribution Of Kaon   Partilces",ptBin,nArrayPt);
    HistMCTruthEtaKaon= new TH1F("HistMCTruthEtaKaon","Eta Distribution Of Kaon   Partilces",etaBin,nArrayEta);
    HistMCTruthPhiKaon= new TH1F("HistMCTruthPhiKaon","Phi Distribution Of Kaon   Partilces",phiBin,nArrayPhi);

    HistMCTruthPtProton= new TH1F("HistMCTruthPtProton","Pt Distribution Of Proton   Partilces",ptBin,nArrayPt);
    HistMCTruthEtaProton= new TH1F("HistMCTruthEtaProton","Eta Distribution Of Proton   Partilces",etaBin,nArrayEta);
    HistMCTruthPhiProton= new TH1F("HistMCTruthPhiProton","Phi Distribution Of Proton   Partilces",phiBin,nArrayPhi);

  fOutputList->Add(HistMCTruthPtAll);
  fOutputList->Add(HistMCTruthEtaAll);
  fOutputList->Add(HistMCTruthPhiAll);
  fOutputList->Add(HistMCTruthPtPion);
  fOutputList->Add(HistMCTruthEtaPion);
  fOutputList->Add(HistMCTruthPhiPion);
  fOutputList->Add(HistMCTruthPtKaon);
  fOutputList->Add(HistMCTruthEtaKaon);
  fOutputList->Add(HistMCTruthPhiKaon);
  fOutputList->Add(HistMCTruthPtProton);
  fOutputList->Add(HistMCTruthEtaProton);
  fOutputList->Add(HistMCTruthPhiProton);

  //Contamination and efficiency 
  fHistTruthPionPlus = new TH3D("fHistTruthPionPlus","TruthPionPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPionPlus);

  fHistTruthPionMinus = new TH3D("fHistTruthPionMinus","TruthPionMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPionMinus);

  fHistTruthKaonPlus = new TH3D("fHistTruthKaonPlus","TruthKaonplus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaonPlus);

  fHistTruthKaonMinus = new TH3D("fHistTruthKaonMinus","TruthKaonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaonMinus);

  fHistTruthProtonPlus = new TH3D("fHistTruthProtonPlus","TruthProtonPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProtonPlus);

  fHistTruthProtonMinus = new TH3D("fHistTruthProtonMinus","TruthProtonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProtonMinus);

  fHistTruthPion = new TH3D("fHistTruthPion","TruthPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPion);
  
  fHistTruthKaon = new TH3D("fHistTruthKaon","TruthKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaon);

  fHistTruthProton = new TH3D("fHistTruthProton","TruthProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProton);

 
  //Contamination and Purity Histogram
 
   HistPionContaminationInPt=new TH1F("HistPionContaminationInPt","Pt Distribution of conatmination in Pion", ptBin,nArrayPt);
   HistPionPlusContaminationInPt=new TH1F("HistPionPlusContaminationInPt","Pt Distribution of conatmination in Positive Pion", ptBin,nArrayPt);
   HistPionMinusContaminationInPt=new TH1F("HistPionMinusContaminationInPt","Pt Distribution of conatmination in Negative Pion", ptBin,nArrayPt);
   Hist3dPionContamination =new TH3D("Hist3dPionContamination","Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionPlusContamination =new TH3D("Hist3dPionPlusContamination","Positive Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionMinusContamination =new TH3D("Hist3dPionMinusContamination","Negative Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);

   HistKaonContaminationInPt=new TH1F("HistKaonContaminationInPt","Pt Distribution of conatmination in Kaon", ptBin,nArrayPt);
   HistKaonPlusContaminationInPt=new TH1F("HistKaonPlusContaminationInPt","Pt Distribution of conatmination in Positive Kaon", ptBin,nArrayPt);
   HistKaonMinusContaminationInPt=new TH1F("HistKaonMinusContaminationInPt","Pt Distribution of conatmination in Negative Kaon", ptBin,nArrayPt);
   Hist3dKaonContamination =new TH3D("Hist3dKaonContamination","Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dKaonPlusContamination =new TH3D("Hist3dKaonPlusContamination","Positive Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dKaonMinusContamination =new TH3D("Hist3dKaonMinusContamination","Negative Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);


   HistProtonContaminationInPt=new TH1F("HistProtonContaminationInPt","Pt Distribution of conatmination in Proton", ptBin,nArrayPt);
   HistProtonPlusContaminationInPt=new TH1F("HistProtonPlusContaminationInPt","Pt Distribution of conatmination in Positive Proton", ptBin,nArrayPt);
   HistProtonMinusContaminationInPt=new TH1F("HistProtonMinusContaminationInPt","Pt Distribution of conatmination in Negative Proton", ptBin,nArrayPt);
   Hist3dProtonContamination =new TH3D("Hist3dProtonContamination","Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dProtonPlusContamination =new TH3D("Hist3dProtonPlusContamination","Positive Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dProtonMinusContamination =new TH3D("Hist3dProtonMinusContamination","Negative Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);


   HistPionPurityInPt=new TH1F("HistPionPurityInPt","Pt Distribution of Purity in Pion", ptBin,nArrayPt);
   HistPionPlusPurityInPt=new TH1F("HistPionPlusPurityInPt","Pt Distribution of Purity in Positive Pion", ptBin,nArrayPt);
   HistPionMinusPurityInPt=new TH1F("HistPionMinusPurityInPt","Pt Distribution of Purity in Negative Pion", ptBin,nArrayPt);
   Hist3dPionPurity =new TH3D("Hist3dPionPurity","Pion Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionPlusPurity =new TH3D("Hist3dPionPlusPurity","Positive Pion Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionMinusPurity =new TH3D("Hist3dPionMinusPurity","Negative Pion Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);

   HistKaonPurityInPt=new TH1F("HistKaonPurityInPt","Pt Distribution of Purity in Kaon", ptBin,nArrayPt);
   HistKaonPlusPurityInPt=new TH1F("HistKaonPlusPurityInPt","Pt Distribution of Purity in Positive Kaon", ptBin,nArrayPt);
   HistKaonMinusPurityInPt=new TH1F("HistKaonMinusPurityInPt","Pt Distribution of Purity in Negative Kaon", ptBin,nArrayPt);
   Hist3dKaonPurity =new TH3D("Hist3dKaonPurity","Kaon Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dKaonPlusPurity =new TH3D("Hist3dKaonPlusPurity","Positive Kaon Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dKaonMinusPurity =new TH3D("Hist3dKaonMinusPurity","Negative Kaon Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);


   HistProtonPurityInPt=new TH1F("HistProtonPurityInPt","Pt Distribution of Purity in Proton", ptBin,nArrayPt);
   HistProtonPlusPurityInPt=new TH1F("HistProtonPlusPurityInPt","Pt Distribution of Purity in Positive Proton", ptBin,nArrayPt);
   HistProtonMinusPurityInPt=new TH1F("HistProtonMinusPurityInPt","Pt Distribution of Purity in Negative Proton", ptBin,nArrayPt);
   Hist3dProtonPurity =new TH3D("Hist3dProtonPurity","Proton Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dProtonPlusPurity =new TH3D("Hist3dProtonPlusPurity","Positive Proton Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dProtonMinusPurity =new TH3D("Hist3dProtonMinusPurity","Negative Proton Purity;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);

 fOutputList->Add(HistPionContaminationInPt);
 fOutputList->Add(HistPionPlusContaminationInPt);
 fOutputList->Add(HistPionMinusContaminationInPt);
 fOutputList->Add(Hist3dPionContamination);
 fOutputList->Add(Hist3dPionPlusContamination);
 fOutputList->Add(Hist3dPionMinusContamination);
 fOutputList->Add(HistKaonContaminationInPt);
 fOutputList->Add(HistKaonPlusContaminationInPt);
 fOutputList->Add(HistKaonMinusContaminationInPt);
 fOutputList->Add(Hist3dKaonContamination);
 fOutputList->Add(Hist3dKaonPlusContamination);
 fOutputList->Add(Hist3dKaonMinusContamination);
 fOutputList->Add(HistProtonContaminationInPt);
 fOutputList->Add(HistProtonPlusContaminationInPt);
 fOutputList->Add(HistProtonMinusContaminationInPt);
 fOutputList->Add(Hist3dProtonContamination);
 fOutputList->Add(Hist3dProtonPlusContamination);
 fOutputList->Add(Hist3dProtonMinusContamination); 
 fOutputList->Add(HistPionPurityInPt);
 fOutputList->Add(HistPionPlusPurityInPt);
 fOutputList->Add(HistPionMinusPurityInPt);
 fOutputList->Add(Hist3dPionPurity);
 fOutputList->Add(Hist3dPionPlusPurity);
 fOutputList->Add(Hist3dPionMinusPurity);
 fOutputList->Add(HistKaonPurityInPt);
 fOutputList->Add(HistKaonPlusPurityInPt);
 fOutputList->Add(HistKaonMinusPurityInPt);
 fOutputList->Add(Hist3dKaonPurity);
 fOutputList->Add(Hist3dKaonPlusPurity);
 fOutputList->Add(Hist3dKaonMinusPurity);
 fOutputList->Add(HistProtonPurityInPt);
 fOutputList->Add(HistProtonPlusPurityInPt);
 fOutputList->Add(HistProtonMinusPurityInPt);
 fOutputList->Add(Hist3dProtonPurity);
 fOutputList->Add(Hist3dProtonPlusPurity);
 fOutputList->Add(Hist3dProtonMinusPurity);

// Contamination and Purity  Histogram 

//-----------------------------------------------------------------------------------------------------------------------------------------------------

  fHistMCRecoPionPlus = new TH3D("fHistMCRecoPionPlus","MCRecoPionPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionPlus);

  fHistMCRecoPionMinus = new TH3D("fHistMCRecoPionMinus","MCRecoPionMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionMinus);

  fHistMCRecoKaonPlus = new TH3D("fHistMCRecoKaonPlus","MCRecoKaonplus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonPlus);

  fHistMCRecoKaonMinus = new TH3D("fHistMCRecoKaonMinus","MCRecoKaonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonMinus);

  fHistMCRecoProtonPlus = new TH3D("fHistMCRecoProtonPlus","MCRecoProtonPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonPlus);

  fHistMCRecoProtonMinus = new TH3D("fHistMCRecoProtonMinus","MCRecoProtonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonMinus);

  fHistMCRecoPion = new TH3D("fHistMCRecoPion","MCRecoPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPion);
  
  fHistMCRecoKaon = new TH3D("fHistMCRecoKaon","MCRecoKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaon);

  fHistMCRecoProton = new TH3D("fHistMCRecoProton","MCRecoProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProton);

  fHistMCRecoPionAsKaon = new TH3D("fHistMCRecoPionAsKaon","MCRecoPionAsKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionAsKaon);

  fHistMCRecoPionAsProton = new TH3D("fHistMCRecoPionAsProton","MCRecoPionAsProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionAsProton);

  fHistMCRecoKaonAsProton = new TH3D("fHistMCRecoKaonAsProton","MCRecoKaonAsProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonAsProton);

  fHistMCRecoKaonAsPion = new TH3D("fHistMCRecoKaonAsPion","MCRecoKaonAsPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonAsPion);

  fHistMCRecoProtonAsKaon = new TH3D("fHistMCRecoProtonAsKaon","MCRecoProtonAsKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonAsKaon);

  fHistMCRecoProtonAsPion = new TH3D("fHistMCRecoProtonAsPion","MCRecoProtonAsPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonAsPion);
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------



  fHistdEdxTPC = new TH2F("fHistdEdxTPC", ";p_{T} (GeV/c);dE/dx (au.)",200,0.0,fMaxPt,500, 0., 500.);
fOutputList->Add(fHistdEdxTPC);
fHistBetaTOF = new TH2F(Form("fHistBetaTOF"), ";p_{T} (GeV/c);v/c",200, 0.0, fMaxPt, 500, 0.1, 1.1);
  fOutputList->Add(fHistBetaTOF);

 fHistSigmaTPCVsTOFPionForPionAfterCut=new TH2F("PionForPionAfterCut","Sigma plot for Pion  as a Pion TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFProtonForPionAfterCut=new TH2F("ProtonForPionAfterCut","Sigma plot for Pion  as a Proton TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFKaonForPionAfterCut=new TH2F("KaonForPionAfterCut","Sigma plot for Pion  as a Kaon TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFPionForKaonAfterCut=new TH2F("PionForKaonAfterCut","Sigma plot for Kaon  as a Pion TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFProtonForKaonAfterCut=new TH2F("ProtonForKaonAfterCut","Sigma plot for Kaon  as a Proton TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFKaonForKaonAfterCut=new TH2F("KaonForKaonAfterCut","Sigma plot for Kaon  as a Kaon TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFPionForProtonAfterCut=new TH2F("PionForProtonAfterCut","Sigma plot for Proton   as a Pion TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFProtonForProtonAfterCut=new TH2F("ProtonForProtonAfterCut","Sigma plot for Proton   as a Proton TPC Vs TOF After cut",500,-100,100,500,-100,100);
 fHistSigmaTPCVsTOFKaonForProtonAfterCut=new TH2F("KaonForProtonAfterCut","Sigma plot for Proton   as a Kaon TPC Vs TOF After cut",500,-100,100,500,-100,100);


 fOutputList->Add(fHistSigmaTPCVsTOFPionForPionAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFProtonForPionAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFKaonForPionAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFPionForKaonAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFProtonForKaonAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFKaonForKaonAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFPionForProtonAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFProtonForProtonAfterCut);
 fOutputList->Add(fHistSigmaTPCVsTOFKaonForProtonAfterCut);


h1PionAfterCut=new TH1F("h1PionAFCut","Pion Pt Distribution For getting Pion Entries After Cut",ptBin,nArrayPt);
 h1KaonAfterCut=new TH1F("h1KaonAFCut","Kaon Pt Distribution For getting Kaon Entries After Cut",ptBin,nArrayPt);
 h1ProtonAfterCut=new TH1F("h1ProtonAFCut","Proton Pt Distribution For getting Proton Entries After Cut",ptBin,nArrayPt);

h1PionAsNonPion=new TH1F("PionAsNonPion","Getting Number of Incorrect Pion after cut , Initially it was Pion ",ptBin,nArrayPt);
h1KaonAsNonKaon=new TH1F("KaonAsNonKaon","Getting Number of Incorrect Kaon after cut , Initially it was Kaon ",ptBin,nArrayPt);
h1ProtonAsNonProton=new TH1F("ProtonAsNonProton","Getting Number of Incorrect Proton after cut , Initially it was Proton ",ptBin,nArrayPt);

fOutputList->Add(h1PionAfterCut);
fOutputList->Add(h1KaonAfterCut);
fOutputList->Add(h1ProtonAfterCut);

fOutputList->Add(h1PionAsNonPion);
fOutputList->Add(h1KaonAsNonKaon);
fOutputList->Add(h1ProtonAsNonProton);

h1PionAsKaon=new TH1F("h1PionAsKaon","MC Pion detect as kaon after cut",ptBin,nArrayPt);
h1PionAsProton=new TH1F("h1PionAsProton","MC Pion detect as Proton after cut",ptBin,nArrayPt);
h1KaonAsPion=new TH1F("h1KaonAsPion","MC Kaon detect as Pion after cut",ptBin,nArrayPt);
h1KaonAsProton=new TH1F("h1KaonAsProton","MC Kaon detect as Proton after cut",ptBin,nArrayPt);
h1ProtonAsPion=new TH1F("h1ProtonAsPion","MC Proton detect as Pion after cut",ptBin,nArrayPt);
h1ProtonAsKaon=new TH1F("h1ProtonAsKaon","MC Proton detect as Kaon after cut",ptBin,nArrayPt);

fOutputList->Add(h1PionAsKaon);
fOutputList->Add(h1PionAsProton);
fOutputList->Add(h1ProtonAsKaon);
fOutputList->Add(h1ProtonAsPion);
fOutputList->Add(h1KaonAsProton);
fOutputList->Add(h1KaonAsPion);
 
  PostData(1, fQAList);
  PostData(2, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEffContPIDBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
  
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;
  }

  fArrayMC = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));   
  if (!fArrayMC)  
    AliFatal("No array of MC particles found !!!"); // MW  no AliFatal use return values   
  
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    AliError("ERROR: Could not retrieve MC event");
    return;
  }

  // PID Response task active?
  if(fElectronRejection) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }

  // ==============================================================================================
  // Copy from AliAnalysisTaskPhiCorrelations:
  // For productions with injected signals, figure out above which label to skip particles/tracks
  Int_t skipParticlesAbove = 0;
  if (fInjectedSignals)
  {
    AliGenEventHeader* eventHeader = 0;
    Int_t headers = 0;
    
    // AOD only
    AliAODMCHeader* header = (AliAODMCHeader*) fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!header){
      AliFatal("fInjectedSignals set but no MC header found");
      return;
    }
    
    headers = header->GetNCocktailHeaders();
    eventHeader = header->GetCocktailHeader(0);
    
    
    if (!eventHeader)
      {
	// We avoid AliFatal here, because the AOD productions sometimes have events where the MC header is missing 
	// (due to unreadable Kinematics) and we don't want to loose the whole job because of a few events
	AliError("First event header not found. Skipping this event.");
	return;
      }
    
    skipParticlesAbove = eventHeader->NProduced();
    AliInfo(Form("Injected signals in this event (%d headers). Keeping particles/tracks of %s. Will skip particles/tracks above %d.", headers, eventHeader->ClassName(), skipParticlesAbove)); 
  }
  // ==============================================================================================

  
  

  fHistEventStats->Fill(1); //all events
  
  //Centrality stuff
  Double_t nCentrality = 0;
  if(fUseCentrality) {
    
    AliAODHeader *headerAOD = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
    if (!headerAOD){
      AliFatal("AOD header found");
      return;
    }

    AliCentrality *centrality = headerAOD->GetCentralityP();
    nCentrality =centrality->GetCentralityPercentile(fCentralityEstimator.Data());
    

    if(!centrality->IsEventInCentralityClass(fCentralityPercentileMin,
					     fCentralityPercentileMax,
					     fCentralityEstimator.Data()))
      return;
    else {    
      fHistEventStats->Fill(2); //triggered + centrality
      fHistCentrality->Fill(nCentrality);
    }
  }
  //Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);

  const AliAODVertex *vertex = fAOD->GetPrimaryVertex(); 
  if(vertex) {
    if(vertex->GetNContributors() > 0) {
      Double32_t fCov[6];    
      vertex->GetCovarianceMatrix(fCov);   
      if(fCov[5] != 0) {
	fHistEventStats->Fill(3); //events with a proper vertex
	if(TMath::Abs(vertex->GetX()) < fVxMax) {    // antes Xv
	  //Printf("X Vertex: %lf", vertex->GetX());
	  //Printf("Y Vertex: %lf", vertex->GetY());
	  if(TMath::Abs(vertex->GetY()) < fVyMax) {  // antes Yv
	    if(TMath::Abs(vertex->GetZ()) < fVzMax) {  // antes Zv
	      //Printf("Z Vertex: %lf", vertex->GetZ());
	      
	      fHistEventStats->Fill(4); //analyzed events
	      fHistVz->Fill(vertex->GetZ()); 
	      
	      //+++++++++++++++++ MC Truth Hsito for Pion, Kaon and Proton ++++++++++++++++++//
	      Int_t nMCParticles = fArrayMC->GetEntriesFast();
	      
	      for(Int_t iParticle = 0; iParticle < nMCParticles; iParticle++) {
               
              AliAODMCParticle* currentAODMCParticle = (AliAODMCParticle*) fArrayMC->At(iParticle);
                if (!IsMCParticleCut(currentAODMCParticle)) {continue;}
                if(!currentAODMCParticle->IsPhysicalPrimary()) continue;

          if (fInjectedSignals && currentAODMCParticle->GetLabel() >= skipParticlesAbove) continue;


// PDG Code  from Track label of MC particle 


/*    AliAODMCParticle* pidCodeMC = (AliAODMCParticle*)fArrayMC->At(TMath::Abs(currentAODMCParticle->GetLabel()));
    Int_t pdgCodeMC = pidCodeMC->GetPdgCode();*/

// PDG Code from Track label of MC particle



            Int_t pdgCode=((AliAODMCParticle*)currentAODMCParticle)->GetPdgCode();
//cout<<" PDG CODE MC from label matching is "<<pdgCodeMC<<'\t'<<" PDG Code of MC without label matchiing is "<<pdgCode<<endl; 
            
             if (TMath::Abs(pdgCode)==11) continue;
	
            Short_t gAODmcCharge = currentAODMCParticle->Charge();


            HistMCTruthPtAll->Fill(currentAODMCParticle->Pt());
            HistMCTruthEtaAll->Fill(currentAODMCParticle->Eta());
            HistMCTruthPhiAll->Fill(currentAODMCParticle->Phi());
    

  
           if (currentAODMCParticle->IsPhysicalPrimary()) {
           switch(TMath::Abs(pdgCode)){
           case 211:
           HistMCTruthPtPion->Fill(currentAODMCParticle->Pt());
           HistMCTruthEtaPion->Fill(currentAODMCParticle->Eta());
           HistMCTruthPhiPion->Fill(currentAODMCParticle->Phi());
           fHistTruthPion->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
           if(gAODmcCharge >0){
           fHistTruthPionPlus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
           if(gAODmcCharge < 0){
           fHistTruthPionMinus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
           break; 
           
           case 321:
           HistMCTruthPtKaon->Fill(currentAODMCParticle->Pt());
           HistMCTruthEtaKaon->Fill(currentAODMCParticle->Eta());
           HistMCTruthPhiKaon->Fill(currentAODMCParticle->Phi());
           fHistTruthKaon->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
           if(gAODmcCharge >0){
           fHistTruthKaonPlus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
           if(gAODmcCharge < 0){
           fHistTruthKaonMinus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
           break;
    
          case 2212:
          HistMCTruthPtProton->Fill(currentAODMCParticle->Pt());
          HistMCTruthEtaProton->Fill(currentAODMCParticle->Eta());
          HistMCTruthPhiProton->Fill(currentAODMCParticle->Phi());
           fHistTruthProton->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
           if(gAODmcCharge >0){
           fHistTruthProtonPlus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
           if(gAODmcCharge < 0){
           fHistTruthProtonMinus->Fill(currentAODMCParticle->Eta(),currentAODMCParticle->Pt(),currentAODMCParticle->Phi());
 }
          break;
  
         } // End of switch

     } // If condition end
	
}//loop over tracks
// ++++++++++++++++++++++++++++++++++++++++ MC Truth Hsito for Pion, Kaon and Proton ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	      
//+++++++++++++++++++++++++++++++++++++++++Reconstruced Partcile and PID selection ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	      
	      //AOD track loop
	      Int_t nGoodTracks = fAOD->GetNumberOfTracks();   
	     
	      
	      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));    
		if(!trackAOD) continue;
		
		//track cuts
		if (!trackAOD->TestFilterBit(fAODTrackCutBit)) 
		  continue;

 		Int_t label = TMath::Abs(trackAOD->GetLabel()); 
		if(label > trackAOD->GetLabel()) continue; 
		
		  if(TMath::Abs(trackAOD->Eta()) > fMaxEta) 
		    continue;
		  if((trackAOD->Pt() > fMaxPt)||(trackAOD->Pt() <  fMinPt)) 
		    continue;

                AliAODMCParticle* recoMC = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(trackAOD->GetLabel())));

               if(!recoMC) continue;

                  if (fInjectedSignals)
                  {

                    AliAODMCParticle* mother = recoMC;

                    // find the primary mother (if not already physical primary)
                    while (!((AliAODMCParticle*)mother)->IsPhysicalPrimary())
                      {
                        if (((AliAODMCParticle*)mother)->GetMother() < 0)
                          {
                            mother = 0;
                            break;
                          }

                        mother = (AliAODMCParticle*) fArrayMC->At(((AliAODMCParticle*)mother)->GetMother());
                        if (!mother)
                          break;
                      }


                    if (!mother)
                      {
                        AliError(Form("WARNING: No mother found for particle %d:", recoMC->GetLabel()));
                        continue;
                      }

                    if (mother->GetLabel() >= skipParticlesAbove)
                      {
                        //AliInfo(Form("Remove particle %d (>= %d)",mother->GetLabel(),skipParticlesAbove));
                        continue;
                      }
                  }


                if (((AliAODMCParticle*) recoMC)->IsSecondaryFromWeakDecay()) continue;



		  
		  Short_t gCharge = trackAOD->Charge();
		  Double_t phiRad = trackAOD->Phi();

//========================================================PID (so far only for electron rejection)===================================================================//		    
		  if(fElectronRejection) {
		    
		    // get the electron nsigma
		    Double_t nSigma = fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kElectron);
		    fHistNSigmaTPCvsPtbeforePID->Fill(trackAOD->Pt(),nSigma);		    
		    
		    // check only for given momentum range
		    if( trackAOD->Pt() > fElectronRejectionMinPt && trackAOD->Pt() < fElectronRejectionMaxPt ){
		      
		      //look only at electron nsigma
		      if(!fElectronOnlyRejection){
			
			//Make the decision based on the n-sigma of electrons
			if(TMath::Abs(nSigma) < fElectronRejectionNSigma) continue;
		      }
		      else{
			
			Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kPion));
			Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kKaon));
			Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kProton));
			
			//Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
			if(TMath::Abs(nSigma) < fElectronRejectionNSigma
			   && nSigmaPions   > fElectronRejectionNSigma
			   && nSigmaKaons   > fElectronRejectionNSigma
			   && nSigmaProtons > fElectronRejectionNSigma ) continue;
		      }
		    }
		    
		    fHistNSigmaTPCvsPtafterPID->Fill(trackAOD->Pt(),nSigma);		    

		  }
//=====================================================end of PID (so far only for electron rejection)================================================================//		  
		  
               Int_t pdgCodeReco = ((AliAODMCParticle*)recoMC)->GetPdgCode();

             if (TMath::Abs(pdgCodeReco)==11) continue;

// PID selection start ================================================================================================================================================

 Double_t dEdx   = trackAOD -> GetTPCsignal();
 fHistdEdxTPC->Fill(trackAOD->Pt(),dEdx);
 
  if(IsTOFPID(trackAOD))
{
Double_t beta = Beta(trackAOD);
fHistBetaTOF->Fill(trackAOD->Pt(), beta);
}

if(trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<fMaxPt && (!IsTOFPID(trackAOD))) continue;

Int_t particleMCReco=-999;
particleMCReco=GetParticleSpecies(trackAOD);

if(particleMCReco == kSpUndefined ) continue;


Double_t nsigmaPionTOF,nsigmaKaonTOF,nsigmaProtonTOF;
Double_t nsigmaPionTPC,nsigmaKaonTPC,nsigmaProtonTPC;
Double_t nsigmaPionTPCTOF,nsigmaKaonTPCTOF,nsigmaProtonTPCTOF;

nsigmaPionTOF= fnsigmas[kSpPion][kNSigmaTOF];
nsigmaKaonTOF= fnsigmas[kSpKaon][kNSigmaTOF];
nsigmaProtonTOF= fnsigmas[kSpProton][kNSigmaTOF];


nsigmaPionTPC= fnsigmas[kSpPion][kNSigmaTPC];
nsigmaKaonTPC= fnsigmas[kSpKaon][kNSigmaTPC];
nsigmaProtonTPC= fnsigmas[kSpProton][kNSigmaTPC];

nsigmaPionTPCTOF = fnsigmas[kSpPion][kNSigmaTPCTOF];
nsigmaKaonTPCTOF = fnsigmas[kSpKaon][kNSigmaTPCTOF];
nsigmaProtonTPCTOF = fnsigmas[kSpProton][kNSigmaTPCTOF];



//cout<<" NSigma TOF Pion "<<nsigmaPionTOF<<'\t'<<"NSigma TOF Kaon "<<nsigmaKaonTOF<<'\t'<<"Nsigma TOF Proton "<<nsigmaProtonTOF<<endl;


//Pion

if(TMath::Abs(pdgCodeReco) == 211) {
if(particleMCReco==kSpPion) {
if(IsTOFPID(trackAOD) && trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<fMaxPt && nsigmaPionTOF!=999 && nsigmaPionTPCTOF<3.0){
fHistSigmaTPCVsTOFPionForPionAfterCut->Fill(nsigmaPionTOF, nsigmaPionTPC);
//cout<<"NSigma Pion TOF is "<<nsigmaPionTOF<<endl;
}

h1PionAfterCut->Fill(trackAOD->Pt());
fHistMCRecoPion->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) fHistMCRecoPionPlus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge<0) fHistMCRecoPionMinus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());

}
else {
if(particleMCReco == kSpProton) {
fHistSigmaTPCVsTOFProtonForPionAfterCut->Fill(nsigmaProtonTOF, nsigmaProtonTPC);
h1PionAsProton->Fill(trackAOD->Pt());
fHistMCRecoPionAsProton->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}
else {
if(particleMCReco == kSpKaon){
fHistSigmaTPCVsTOFKaonForPionAfterCut->Fill(nsigmaKaonTOF, nsigmaKaonTPC);
h1PionAsKaon->Fill(trackAOD->Pt());
fHistMCRecoPionAsKaon->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}

}

h1PionAsNonPion->Fill(trackAOD->Pt());

}

}
//Kaon

if(TMath::Abs(pdgCodeReco) == 321) {
if(particleMCReco==kSpKaon) {
if(IsTOFPID(trackAOD) && trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<fMaxPt && nsigmaKaonTOF!=999 && nsigmaKaonTPCTOF<3.0){
fHistSigmaTPCVsTOFKaonForKaonAfterCut->Fill(nsigmaKaonTOF, nsigmaKaonTPC);
//cout<<"NSigma Kaon TOF is "<<nsigmaKaonTOF<<endl;
}
h1KaonAfterCut->Fill(trackAOD->Pt());
fHistMCRecoKaon->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) fHistMCRecoKaonPlus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge<0) fHistMCRecoKaonMinus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());

}
else {
if(particleMCReco == kSpProton) {
fHistSigmaTPCVsTOFProtonForKaonAfterCut->Fill(nsigmaProtonTOF, nsigmaProtonTPC);
h1KaonAsProton->Fill(trackAOD->Pt());
fHistMCRecoKaonAsProton->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}
else {
if(particleMCReco == kSpPion){
fHistSigmaTPCVsTOFPionForKaonAfterCut->Fill(nsigmaPionTOF, nsigmaPionTPC);
h1KaonAsPion->Fill(trackAOD->Pt());
fHistMCRecoKaonAsPion->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}

}
h1KaonAsNonKaon->Fill(trackAOD->Pt());


}

}

//Proton

if(TMath::Abs(pdgCodeReco) == 2212) {
if(particleMCReco==kSpProton) {
if(IsTOFPID(trackAOD) && trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<fMaxPt && nsigmaProtonTOF!=999 && nsigmaProtonTPCTOF<3.0){
fHistSigmaTPCVsTOFProtonForProtonAfterCut->Fill(nsigmaProtonTOF, nsigmaProtonTPC);
//cout<<"NSigma Proton TOF is "<<nsigmaProtonTOF<<endl;
}

h1ProtonAfterCut->Fill(trackAOD->Pt());
fHistMCRecoProton->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) fHistMCRecoProtonPlus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge<0) fHistMCRecoProtonMinus->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());

}
else {
if(particleMCReco == kSpKaon) {
fHistSigmaTPCVsTOFKaonForProtonAfterCut->Fill(nsigmaKaonTOF, nsigmaKaonTPC);
h1ProtonAsKaon->Fill(trackAOD->Pt());
fHistMCRecoProtonAsKaon->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}
else {
if(particleMCReco == kSpPion){
fHistSigmaTPCVsTOFPionForProtonAfterCut->Fill(nsigmaPionTOF, nsigmaPionTPC);
h1ProtonAsPion->Fill(trackAOD->Pt());
fHistMCRecoProtonAsPion->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
}

}
h1ProtonAsNonProton->Fill(trackAOD->Pt());


}

}

// Fill the Histograms for Contamination
if(TMath::Abs(pdgCodeReco) !=211 && particleMCReco == kSpPion) {
HistPionContaminationInPt->Fill(trackAOD->Pt());
Hist3dPionContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) {
Hist3dPionPlusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
HistPionPlusContaminationInPt->Fill(trackAOD->Pt());
}

if(gCharge<0){
 Hist3dPionMinusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
 HistPionMinusContaminationInPt->Fill(trackAOD->Pt());
}

}

if(TMath::Abs(pdgCodeReco)!=321 && particleMCReco == kSpKaon) {
HistKaonContaminationInPt->Fill(trackAOD->Pt());
Hist3dKaonContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) {
Hist3dKaonPlusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
HistKaonPlusContaminationInPt->Fill(trackAOD->Pt());
}

if(gCharge<0){
 Hist3dKaonMinusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
 HistKaonMinusContaminationInPt->Fill(trackAOD->Pt());
}

}

if(TMath::Abs(pdgCodeReco)!=2212 && particleMCReco == kSpProton) {
HistProtonContaminationInPt->Fill(trackAOD->Pt());
Hist3dProtonContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
if(gCharge>0) {
Hist3dProtonPlusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
HistProtonPlusContaminationInPt->Fill(trackAOD->Pt());
}

if(gCharge<0){
 Hist3dProtonMinusContamination->Fill(trackAOD->Eta(),trackAOD->Pt(),trackAOD->Phi());
 HistProtonMinusContaminationInPt->Fill(trackAOD->Pt());
}

}
// Fill the Histograms for Contamination




// PID selection end ================================================================================================================================================

		
	       }//AOD track loop 
	      
	    }//Vz cut
	  }//Vy cut
	}//Vx cut
      }//Vz resolution
    }//number of contributors
  }//valid vertex  
}


Bool_t AliAnalysisTaskEffContPIDBF::IsMCParticleCut(AliAODMCParticle* particle) {


//cout<<"Min Pt ---- "<<fMinPt<<"Max Pt "<<fMaxPt<<"Pt TPC Max "<<fPtTPCMax<<endl;


        if (particle->Charge() == 0) {return kFALSE;}
        if ((particle->Pt() < fMinPt) || (particle->Pt() > fMaxPt)) {return kFALSE;}

        if ((TMath::Abs(particle->Eta()) > fMaxEta)) {return kFALSE;}

        return kTRUE;
}

void AliAnalysisTaskEffContPIDBF::SigmaCalculate(AliAODTrack *trk ){

  Double_t nsigmaTPCkProton = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kProton);
  Double_t nsigmaTPCkKaon   = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kKaon);
  Double_t nsigmaTPCkPion   = fPIDResponse->NumberOfSigmasTPC(trk, AliPID::kPion);




 Double_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;
 Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;

 if(IsTOFPID(trk) && trk->Pt()>fPtTPCMax && trk->Pt()<fMaxPt){
//cout<<" Hi I am now at TPC+TOF  track "<<endl;
    nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kProton);
    nsigmaTOFkKaon   = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kKaon);
    nsigmaTOFkPion   = fPIDResponse->NumberOfSigmasTOF(trk, AliPID::kPion);
    Double_t d2Proton=nsigmaTPCkProton * nsigmaTPCkProton + nsigmaTOFkProton * nsigmaTOFkProton;
    Double_t d2Kaon=nsigmaTPCkKaon * nsigmaTPCkKaon + nsigmaTOFkKaon * nsigmaTOFkKaon;
    Double_t d2Pion=nsigmaTPCkPion * nsigmaTPCkPion + nsigmaTOFkPion * nsigmaTOFkPion;

    nsigmaTPCTOFkProton  =  TMath::Sqrt(d2Proton);
    nsigmaTPCTOFkKaon    =  TMath::Sqrt(d2Kaon);
    nsigmaTPCTOFkPion    =  TMath::Sqrt(d2Pion);

}

else{
 if(IsTPCPID(trk) && trk->Pt()>=fMinPt && trk->Pt()<=fPtTPCMax){
    nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
    nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
    nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
}

  }

  fnsigmas[kSpPion][kNSigmaTPC]=nsigmaTPCkPion;
  fnsigmas[kSpKaon][kNSigmaTPC]=nsigmaTPCkKaon;
  fnsigmas[kSpProton][kNSigmaTPC]=nsigmaTPCkProton;

  fnsigmas[kSpPion][kNSigmaTPCTOF]=nsigmaTPCTOFkPion;
  fnsigmas[kSpKaon][kNSigmaTPCTOF]=nsigmaTPCTOFkKaon;
  fnsigmas[kSpProton][kNSigmaTPCTOF]=nsigmaTPCTOFkProton;
 
  
   fnsigmas[kSpPion][kNSigmaTOF]=nsigmaTOFkPion;
   fnsigmas[kSpKaon][kNSigmaTOF]=nsigmaTOFkKaon;
   fnsigmas[kSpProton][kNSigmaTOF]=nsigmaTOFkProton;

return;

}


Int_t AliAnalysisTaskEffContPIDBF::SigmaCutForParticleSpecies(AliAODTrack *trk ){


Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;

if(trk->Pt()>fPtTPCMax && trk->Pt()<fMaxPt && (!IsTOFPID(trk))) return kSpUndefined;

  switch (fPIDType){
  case kNSigmaTPC:
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPC]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPC])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPC])  ;

   //cout<<" I am using TPC detector "<<endl;
    break;

  case kNSigmaTPCTOF://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTPCTOF]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTPCTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTPCTOF])  ;
   //cout<<" I am using TPCTOF detector "<<endl;
    break;

  case kNSigmaTOF:

    nsigmaProton  =  TMath::Abs(fnsigmas[kSpProton][kNSigmaTOF]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kSpKaon][kNSigmaTOF])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kSpPion][kNSigmaTOF])  ;
   //cout<<" I am using TOF detector "<<endl;
    break;



}

if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return kSpUndefined;

//Kaon

if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < fNSigmaPID)) return kSpKaon;


//Pion

 if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < fNSigmaPID))  return kSpPion;


//Proton

 if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fNSigmaPID))  return kSpProton;

return kSpUndefined;

}


Int_t AliAnalysisTaskEffContPIDBF::GetParticleSpecies(AliAODTrack *trk )
{
  if (!trk) {return kSpUndefined; }

  SigmaCalculate(trk);
  Int_t mypid = -1;
  mypid = SigmaCutForParticleSpecies(trk);
  //Printf(" >>>>>>>>>>>>>>>>> mypid = %d",mypid);
  return mypid;
}


Double_t AliAnalysisTaskEffContPIDBF::Beta(AliAODTrack *track)
{
  Double_t stoptime=track->GetTOFsignal();

  Double_t c=TMath::C()*1.E-9;// m/ns
  Float_t startTime = fPIDResponse->GetTOFResponse().GetStartTime(((AliVTrack*)track)->P());//in ps
  Double_t length= fPIDResponse->GetTOFResponse().GetExpectedSignal(track,AliPID::kElectron)*1E-3*c;
  stoptime -= startTime;      
  Double_t scaleStopTime= stoptime*1E-3;          
  scaleStopTime=scaleStopTime*c;
  return length/scaleStopTime;
}





Bool_t AliAnalysisTaskEffContPIDBF::IsTPCPID(AliAODTrack* track) const
{
  // check PID signal 
  AliPIDResponse::EDetPidStatus statusTPC = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC,track);
  if (statusTPC != AliPIDResponse::kDetPidOk)
    return kFALSE;
  return kTRUE;
}




Bool_t AliAnalysisTaskEffContPIDBF::IsTOFPID(AliAODTrack* track) const
{
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse -> CheckPIDStatus(AliPIDResponse::kTOF,track);
  if(statusTOF != AliPIDResponse::kDetPidOk)
    return kFALSE;
  if(0)
  {
    Int_t startTimeMask = fPIDResponse -> GetTOFResponse().GetStartTimeMask(track->P());
    if (startTimeMask < 0) return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskEffContPIDBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}

//____________________________________________________________________//
Bool_t AliAnalysisTaskEffContPIDBF::IsLabelUsed(TArrayI labelArray, Int_t label) {
  //Checks if the label is used already
  Bool_t status = kFALSE;
  for(Int_t i = 0; i < labelArray.GetSize(); i++) {
    if(labelArray.At(i) == label)
      status = kTRUE;
  }

  return status;
}
