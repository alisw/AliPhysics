#include "TChain.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
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
// By Noor Alam (VECC ,Kolkata) sk.noor.alam@cern.ch and Subhasis Chattopadhyay(sub.chattopadhyay@gmail.com)
//[ Special thanks to Michael Weber  ]

// ---------------------------------------------------------------------

ClassImp(AliAnalysisTaskEffContPIDBF)

//________________________________________________________________________
AliAnalysisTaskEffContPIDBF::AliAnalysisTaskEffContPIDBF(const char *name) 
  : AliAnalysisTaskSE(name), 
    fAOD(0),
    fNSigmaPID(3.0),
    fArrayMC(0), 
    fQAList(0), 
    fOutputList(0), 
    fHistEventStats(0), 
    fHistCentrality(0),
    fHistVz(0), 
    fHistNSigmaTPCvsPtbeforeElectronPID(0),
    fHistNSigmaTPCvsPtafterElectronPID(0),  
    fHistTruthPionPlus(0),
    fHistTruthKaonPlus(0),
    fHistTruthProtonPlus(0),
    fHistTruthPionMinus(0),
    fHistTruthKaonMinus(0),
    fHistTruthProtonMinus(0),
    fHistTruthPion(0),
    fHistTruthKaon(0),
    fHistTruthProton(0),
    Hist3dPionContamination(0),
    Hist3dPionPlusContamination(0),
    Hist3dPionMinusContamination(0),
    Hist3dKaonContamination(0),
    Hist3dKaonPlusContamination(0),
    Hist3dKaonMinusContamination(0),
    Hist3dProtonContamination(0),
    Hist3dProtonPlusContamination(0),
    Hist3dProtonMinusContamination(0),
    fHistMCRecoPionPlus(0),
    fHistMCRecoKaonPlus(0),
    fHistMCRecoProtonPlus(0),
    fHistMCRecoPionMinus(0),
    fHistMCRecoKaonMinus(0),
    fHistMCRecoProtonMinus(0),
    fHistMCRecoPion(0),
    fHistMCRecoKaon(0),
    fHistMCRecoProton(0),
    fHistNsigmaTPCPionBeforePIDCut(0),
    fHistNsigmaTPCKaonBeforePIDCut(0),
    fHistNsigmaTPCProtonBeforePIDCut(0),
    fHistNsigmaTOFPionBeforePIDCut(0),
    fHistNsigmaTOFKaonBeforePIDCut(0),
    fHistNsigmaTOFProtonBeforePIDCut(0),
    fHistNsigmaTPCTOFPionBeforePIDCut(0),
    fHistNsigmaTPCTOFKaonBeforePIDCut(0),
    fHistNsigmaTPCTOFProtonBeforePIDCut(0),
    fHistNsigmaTPCTOFPionAfterPIDCut(0),
    fHistNsigmaTPCTOFKaonAfterPIDCut(0),
    fHistNsigmaTPCTOFProtonAfterPIDCut(0),
    fHistdEdxTPCPionAfterPIDCut(0),
    fHistdEdxTPCKaonAfterPIDCut(0),
    fHistdEdxTPCProtonAfterPIDCut(0),
    fHistBetaTOFPionAfterPIDCut(0),
    fHistBetaTOFKaonAfterPIDCut(0),
    fHistBetaTOFProtonAfterPIDCut(0),
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
    fParticleType_(kPion_),
    fDetectorPID_(kTPCTOFpid_),
    fRapidityInsteadOfEta(kFALSE),
    fMistMatchTOFProb(0.0),
    fTOFMisMatch(kFALSE),
    fZvertexTruthPion(0),
    fZvertexTruthPionPlus(0),
    fZvertexTruthPionMinus(0),
    fZvertexContaminationPion(0),
    fZvertexContaminationPionPlus(0),
    fZvertexContaminationPionMinus(0),
    fZvertexRecoPion(0),
    fZvertexRecoPionPlus(0),
    fZvertexRecoPionMinus(0),
    fZvertexTruthKaon(0),
    fZvertexTruthKaonPlus(0),
    fZvertexTruthKaonMinus(0),
    fZvertexContaminationKaon(0),
    fZvertexContaminationKaonPlus(0),
    fZvertexContaminationKaonMinus(0),
    fZvertexRecoKaon(0),
    fZvertexRecoKaonPlus(0),
    fZvertexRecoKaonMinus(0),
    fZvertexTruthProton(0),
    fZvertexTruthProtonPlus(0),
    fZvertexTruthProtonMinus(0),
    fZvertexContaminationProton(0),
    fZvertexContaminationProtonPlus(0),
    fZvertexContaminationProtonMinus(0),
    fZvertexRecoProton(0),
    fZvertexRecoProtonPlus(0),
    fZvertexRecoProtonMinus(0),
    fHistEtaVertexzTruthPion(0),
    fHistEtaVertexzTruthPionPlus(0),
    fHistEtaVertexzTruthPionMinus(0),
    fHistEtaVertexzTruthKaon(0),
    fHistEtaVertexzTruthKaonPlus(0),
    fHistEtaVertexzTruthKaonMinus(0),
    fHistEtaVertexzTruthProton(0),
    fHistEtaVertexzTruthProtonPlus(0),
    fHistEtaVertexzTruthProtonMinus(0),
    fHistEtaVertexzRecoPion(0),
    fHistEtaVertexzRecoPionPlus(0),
    fHistEtaVertexzRecoPionMinus(0),
    fHistEtaVertexzRecoKaon(0),
    fHistEtaVertexzRecoKaonPlus(0),
    fHistEtaVertexzRecoKaonMinus(0),
    fHistEtaVertexzRecoProton(0),
    fHistEtaVertexzRecoProtonPlus(0),
    fHistEtaVertexzRecoProtonMinus(0),
    fHistEtaVertexzContaminationPion(0),
    fHistEtaVertexzContaminationPionPlus(0),
    fHistEtaVertexzContaminationPionMinus(0),
    fHistEtaVertexzContaminationKaon(0),
    fHistEtaVertexzContaminationKaonPlus(0),
    fHistEtaVertexzContaminationKaonMinus(0),
    fHistEtaVertexzContaminationProton(0),
    fHistEtaVertexzContaminationProtonPlus(0),
    fHistEtaVertexzContaminationProtonMinus(0),
    fHistPhiVertexzTruthPion(0),
    fHistPhiVertexzTruthPionPlus(0),
    fHistPhiVertexzTruthPionMinus(0),
    fHistPhiVertexzTruthKaon(0),
    fHistPhiVertexzTruthKaonPlus(0),
    fHistPhiVertexzTruthKaonMinus(0),
    fHistPhiVertexzTruthProton(0),
    fHistPhiVertexzTruthProtonPlus(0),
    fHistPhiVertexzTruthProtonMinus(0),
    fHistPhiVertexzRecoPion(0),
    fHistPhiVertexzRecoPionPlus(0),
    fHistPhiVertexzRecoPionMinus(0),
    fHistPhiVertexzRecoKaon(0),
    fHistPhiVertexzRecoKaonPlus(0),
    fHistPhiVertexzRecoKaonMinus(0),
    fHistPhiVertexzRecoProton(0),
    fHistPhiVertexzRecoProtonPlus(0),
    fHistPhiVertexzRecoProtonMinus(0),
    fHistPhiVertexzContaminationPion(0),
    fHistPhiVertexzContaminationPionPlus(0),
    fHistPhiVertexzContaminationPionMinus(0),
    fHistPhiVertexzContaminationKaon(0),
    fHistPhiVertexzContaminationKaonPlus(0),
    fHistPhiVertexzContaminationKaonMinus(0),
    fHistPhiVertexzContaminationProton(0),
    fHistPhiVertexzContaminationProtonPlus(0),
    fHistPhiVertexzContaminationProtonMinus(0)

{ 
    for(Int_t ipart=0;ipart<6;ipart++){
   fNsigmaTPC[ipart]=999.0;
   fNsigmaTOF[ipart]=999.0;
   }

  
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
  
  Int_t ptBin = 44;
  Int_t etaBin = 16;
  Int_t phiBin = 100;

  Int_t VertexZBin=10;

   Double_t nArrayPt[45]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};

  Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Double_t nArrayPhi[phiBin+1];
  for(Int_t iBin = 0; iBin <= phiBin; iBin++) 
  nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;


  Double_t nArrayVertexZ[11]={-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};

  //AOD analysis
  fHistCentrality = new TH1F("fHistCentrality",";Centrality bin;Events",1001,-0.5,100.5);
  fQAList->Add(fHistCentrality);
  
  //multiplicity (good MC tracks)

  //Vz addition+++++++++++++++++++++++++++++
  fHistVz = new TH1F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Entries",100,-20.,20.);
  fQAList->Add(fHistVz);

  //Electron cuts -> PID QA
  fHistNSigmaTPCvsPtbeforeElectronPID = new TH2F ("ElectronNSigmaTPCvsPtbefore","Electron NSigmaTPCvsPtbefore PID Cut",1000, -10,10,1000, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtbeforeElectronPID);

  fHistNSigmaTPCvsPtafterElectronPID = new TH2F ("ElectronNSigmaTPCvsPtafter","Electron NSigmaTPCvsPtafter PID Cut",1000, -10, 10, 1000, -10, 10); 
  fQAList->Add(fHistNSigmaTPCvsPtafterElectronPID);

 fHistNsigmaTPCPionBeforePIDCut=new TH2F("HistNsigmaTPCvsPtPionBeforePIDCut","NsigmaTPC vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCKaonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtKaonBeforePIDCut","NsigmaTPC vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCProtonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtProtonBeforePIDCut","NsigmaTPC vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fHistNsigmaTOFPionBeforePIDCut=new TH2F("HistNsigmaTOFvsPtPionBeforePIDCut","NsigmaTOF vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFKaonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtKaonBeforePIDCut","NsigmaTOF vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFProtonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtProtonBeforePIDCut","NsigmaTOF vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);
 

 fHistNsigmaTPCTOFPionBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtPionBeforePIDCut","NsigmaTPCTOF vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFKaonBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtKaonBeforePIDCut","NsigmaTPCTOF vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFProtonBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtProtonBeforePIDCut","NsigmaTPCTOF vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fHistNsigmaTPCTOFPionAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtPionAfterPIDCut","NsigmaTPCTOF vs Pt of  Pion AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFKaonAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtKaonAfterPIDCut","NsigmaTPCTOF vs Pt of  Kaon AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFProtonAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtProtonAfterPIDCut","NsigmaTPCTOF vs Pt of  Proton AfterPIDCut",1000, 0,10,1000, -10, 10);



 fQAList->Add(fHistNsigmaTPCPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCProtonBeforePIDCut);

 fQAList->Add(fHistNsigmaTOFPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFProtonBeforePIDCut);


 fQAList->Add(fHistNsigmaTPCTOFPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCTOFKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCTOFProtonBeforePIDCut);

 fQAList->Add(fHistNsigmaTPCTOFPionAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCTOFKaonAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCTOFProtonAfterPIDCut);

  //Contamination and efficiency 
  fHistTruthPionPlus = new TH3F("fHistTruthPionPlus","TruthPionPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPionPlus);

  fHistTruthPionMinus = new TH3F("fHistTruthPionMinus","TruthPionMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPionMinus);

  fHistTruthKaonPlus = new TH3F("fHistTruthKaonPlus","TruthKaonplus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaonPlus);

  fHistTruthKaonMinus = new TH3F("fHistTruthKaonMinus","TruthKaonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaonMinus);

  fHistTruthProtonPlus = new TH3F("fHistTruthProtonPlus","TruthProtonPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProtonPlus);

  fHistTruthProtonMinus = new TH3F("fHistTruthProtonMinus","TruthProtonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProtonMinus);

  fHistTruthPion = new TH3F("fHistTruthPion","TruthPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthPion);
  
  fHistTruthKaon = new TH3F("fHistTruthKaon","TruthKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthKaon);

  fHistTruthProton = new TH3F("fHistTruthProton","TruthProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistTruthProton);

  fZvertexTruthPion=new TH1F("fZvertexTruthPion","vertex z truth Pion",200,-10,10);
  fOutputList->Add(fZvertexTruthPion);
  fZvertexTruthPionPlus=new TH1F("fZvertexTruthPionPlus","vertex z truth Pion +",200,-10,10);
  fOutputList->Add(fZvertexTruthPionPlus);
  fZvertexTruthPionMinus=new TH1F("fZvertexTruthPionMinus","vertex z truth Pion -",200,-10,10);
  fOutputList->Add(fZvertexTruthPionMinus);

  fZvertexTruthKaon=new TH1F("fZvertexTruthKaon","vertex z truth Kaon",200,-10,10);
  fOutputList->Add(fZvertexTruthKaon);
  fZvertexTruthKaonPlus=new TH1F("fZvertexTruthKaonPlus","vertex z truth Kaon +",200,-10,10);
  fOutputList->Add(fZvertexTruthKaonPlus);
  fZvertexTruthKaonMinus=new TH1F("fZvertexTruthKaonMinus","vertex z truth Kaon -",200,-10,10);
  fOutputList->Add(fZvertexTruthKaonMinus);
 
  fZvertexTruthProton=new TH1F("fZvertexTruthProton","vertex z truth Proton",200,-10,10);
  fOutputList->Add(fZvertexTruthProton);
  fZvertexTruthProtonPlus=new TH1F("fZvertexTruthProtonPlus","vertex z truth Proton +",200,-10,10);
  fOutputList->Add(fZvertexTruthProtonPlus);
  fZvertexTruthProtonMinus=new TH1F("fZvertexTruthProtonMinus","vertex z truth Proton -",200,-10,10);
  fOutputList->Add(fZvertexTruthProtonMinus);
 
  //Contamination and Purity Histogram
 
   Hist3dPionContamination =new TH3F("Hist3dPionContamination","Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionPlusContamination =new TH3F("Hist3dPionPlusContamination","Positive Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
   Hist3dPionMinusContamination =new TH3F("Hist3dPionMinusContamination","Negative Pion Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);

 fOutputList->Add(Hist3dPionContamination);
 fOutputList->Add(Hist3dPionPlusContamination);
 fOutputList->Add(Hist3dPionMinusContamination);
 Hist3dKaonContamination =new TH3F("Hist3dKaonContamination","Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
 Hist3dKaonPlusContamination =new TH3F("Hist3dKaonPlusContamination","Positive Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
 Hist3dKaonMinusContamination =new TH3F("Hist3dKaonMinusContamination","Negative Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);


 Hist3dProtonContamination =new TH3F("Hist3dProtonContamination","Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
 Hist3dProtonPlusContamination =new TH3F("Hist3dProtonPlusContamination","Positive Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
 Hist3dProtonMinusContamination =new TH3F("Hist3dProtonMinusContamination","Negative Proton Contamination;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);


 fOutputList->Add(Hist3dKaonContamination);
 fOutputList->Add(Hist3dKaonPlusContamination);
 fOutputList->Add(Hist3dKaonMinusContamination);


 fOutputList->Add(Hist3dProtonContamination);
 fOutputList->Add(Hist3dProtonPlusContamination);
 fOutputList->Add(Hist3dProtonMinusContamination);



  fZvertexContaminationPion=new TH1F("fZvertexContaminationPion","vertex z Contamination Pion",200,-10,10);
  fOutputList->Add(fZvertexContaminationPion);
  fZvertexContaminationPionPlus=new TH1F("fZvertexContaminationPionPlus","vertex z Contamination Pion +",200,-10,10);
  fOutputList->Add(fZvertexContaminationPionPlus);
  fZvertexContaminationPionMinus=new TH1F("fZvertexContaminationPionMinus","vertex z Contamination Pion -",200,-10,10);
  fOutputList->Add(fZvertexContaminationPionMinus);

  fZvertexContaminationKaon=new TH1F("fZvertexContaminationKaon","vertex z Contamination Kaon",200,-10,10);
  fOutputList->Add(fZvertexContaminationKaon);
  fZvertexContaminationKaonPlus=new TH1F("fZvertexContaminationKaonPlus","vertex z Contamination Kaon +",200,-10,10);
  fOutputList->Add(fZvertexContaminationKaonPlus);
  fZvertexContaminationKaonMinus=new TH1F("fZvertexContaminationKaonMinus","vertex z Contamination Kaon -",200,-10,10);
  fOutputList->Add(fZvertexContaminationKaonMinus);

  fZvertexContaminationProton=new TH1F("fZvertexContaminationProton","vertex z Contamination Proton",200,-10,10);
  fOutputList->Add(fZvertexContaminationProton);
  fZvertexContaminationProtonPlus=new TH1F("fZvertexContaminationProtonPlus","vertex z Contamination Proton +",200,-10,10);
  fOutputList->Add(fZvertexContaminationProtonPlus);
  fZvertexContaminationProtonMinus=new TH1F("fZvertexContaminationProtonMinus","vertex z Contamination Proton -",200,-10,10);
  fOutputList->Add(fZvertexContaminationProtonMinus);


// Contamination and Purity  Histogram 

//-----------------------------------------------------------------------------------------------------------------------------------------------------

  fHistMCRecoPionPlus = new TH3F("fHistMCRecoPionPlus","MCRecoPionPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionPlus);

  fHistMCRecoPionMinus = new TH3F("fHistMCRecoPionMinus","MCRecoPionMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPionMinus);

  fHistMCRecoKaonPlus = new TH3F("fHistMCRecoKaonPlus","MCRecoKaonplus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonPlus);

  fHistMCRecoKaonMinus = new TH3F("fHistMCRecoKaonMinus","MCRecoKaonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaonMinus);

  fHistMCRecoProtonPlus = new TH3F("fHistMCRecoProtonPlus","MCRecoProtonPlus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonPlus);

  fHistMCRecoProtonMinus = new TH3F("fHistMCRecoProtonMinus","MCRecoProtonMinus;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProtonMinus);

  fHistMCRecoPion = new TH3F("fHistMCRecoPion","MCRecoPion;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoPion);
  
  fHistMCRecoKaon = new TH3F("fHistMCRecoKaon","MCRecoKaon;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoKaon);

  fHistMCRecoProton = new TH3F("fHistMCRecoProton","MCRecoProton;#eta;p_{T} (GeV/c);#varphi",etaBin,nArrayEta,ptBin,nArrayPt,phiBin,nArrayPhi);
  fOutputList->Add(fHistMCRecoProton);

  fZvertexRecoPion=new TH1F("fZvertexRecoPion","vertex z Reco Pion",200,-10,10);
  fOutputList->Add(fZvertexRecoPion);
  fZvertexRecoPionPlus=new TH1F("fZvertexRecoPionPlus","vertex z Reco Pion +",200,-10,10);
  fOutputList->Add(fZvertexRecoPionPlus);
  fZvertexRecoPionMinus=new TH1F("fZvertexRecoPionMinus","vertex z Reco Pion -",200,-10,10);
  fOutputList->Add(fZvertexRecoPionMinus);

  fZvertexRecoKaon=new TH1F("fZvertexRecoKaon","vertex z Reco Kaon",200,-10,10);
  fOutputList->Add(fZvertexRecoKaon);
  fZvertexRecoKaonPlus=new TH1F("fZvertexRecoKaonPlus","vertex z Reco Kaon +",200,-10,10);
  fOutputList->Add(fZvertexRecoKaonPlus);
  fZvertexRecoKaonMinus=new TH1F("fZvertexRecoKaonMinus","vertex z Reco Kaon -",200,-10,10);
  fOutputList->Add(fZvertexRecoKaonMinus);

  fZvertexRecoProton=new TH1F("fZvertexRecoProton","vertex z Reco Proton",200,-10,10);
  fOutputList->Add(fZvertexRecoProton);
  fZvertexRecoProtonPlus=new TH1F("fZvertexRecoProtonPlus","vertex z Reco Proton +",200,-10,10);
  fOutputList->Add(fZvertexRecoProtonPlus);
  fZvertexRecoProtonMinus=new TH1F("fZvertexRecoProtonMinus","vertex z Reco Proton -",200,-10,10);
  fOutputList->Add(fZvertexRecoProtonMinus);

//-----------------------------------------------------------------------------------------------------------------------------------------------------------------


// Efficieny in Eta and Phi for Different Vertex-z cut 

    fHistEtaVertexzTruthPion=new TH2F("fHistEtaVertexzTruthPion","Eta-Vertexz TruthPion;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthPionPlus=new TH2F("fHistEtaVertexzTruthPionPlus","Eta-Vertexz TruthPion +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthPionMinus=new TH2F("fHistEtaVertexzTruthPionMinus","Eta-Vertexz TruthPion -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzTruthKaon=new TH2F("fHistEtaVertexzTruthKaon","Eta-Vertexz TruthKaon;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthKaonPlus=new TH2F("fHistEtaVertexzTruthKaonPlus","Eta-Vertexz TruthKaon +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthKaonMinus=new TH2F("fHistEtaVertexzTruthKaonMinus","Eta-Vertexz TruthKaon -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzTruthProton=new TH2F("fHistEtaVertexzTruthProton","Eta-Vertexz TruthProton;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthProtonPlus=new TH2F("fHistEtaVertexzTruthProtonPlus","Eta-Vertexz TruthProton +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzTruthProtonMinus=new TH2F("fHistEtaVertexzTruthProtonMinus","Eta-Vertexz TruthProton -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    
    fHistEtaVertexzRecoPion=new TH2F("fHistEtaVertexzRecoPion","Eta-Vertexz RecoPion;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoPionPlus=new TH2F("fHistEtaVertexzRecoPionPlus","Eta-Vertexz RecoPion +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoPionMinus=new TH2F("fHistEtaVertexzRecoPionMinus","Eta-Vertexz RecoPion -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzRecoKaon=new TH2F("fHistEtaVertexzRecoKaon","Eta-Vertexz RecoKaon;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoKaonPlus=new TH2F("fHistEtaVertexzRecoKaonPlus","Eta-Vertexz RecoKaon +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoKaonMinus=new TH2F("fHistEtaVertexzRecoKaonMinus","Eta-Vertexz RecoKaon -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzRecoProton=new TH2F("fHistEtaVertexzRecoProton","Eta-Vertexz RecoProton;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoProtonPlus=new TH2F("fHistEtaVertexzRecoProtonPlus","Eta-Vertexz RecoProton +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzRecoProtonMinus=new TH2F("fHistEtaVertexzRecoProtonMinus","Eta-Vertexz RecoProton -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzContaminationPion=new TH2F("fHistEtaVertexzContaminationPion","Eta-Vertexz ContaminationPion;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationPionPlus=new TH2F("fHistEtaVertexzContaminationPionPlus","Eta-Vertexz ContaminationPion +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationPionMinus=new TH2F("fHistEtaVertexzContaminationPionMinus","Eta-Vertexz ContaminationPion -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzContaminationKaon=new TH2F("fHistEtaVertexzContaminationKaon","Eta-Vertexz ContaminationKaon;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationKaonPlus=new TH2F("fHistEtaVertexzContaminationKaonPlus","Eta-Vertexz ContaminationKaon +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationKaonMinus=new TH2F("fHistEtaVertexzContaminationKaonMinus","Eta-Vertexz ContaminationKaon -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);

    fHistEtaVertexzContaminationProton=new TH2F("fHistEtaVertexzContaminationProton","Eta-Vertexz ContaminationProton;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationProtonPlus=new TH2F("fHistEtaVertexzContaminationProtonPlus","Eta-Vertexz ContaminationProton +;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);
    fHistEtaVertexzContaminationProtonMinus=new TH2F("fHistEtaVertexzContaminationProtonMinus","Eta-Vertexz ContaminationProton -;#eta;Vz",etaBin,nArrayEta,VertexZBin,nArrayVertexZ);


    fHistPhiVertexzTruthPion=new TH2F("fHistPhiVertexzTruthPion","Phi-Vertexz TruthPion;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthPionPlus=new TH2F("fHistPhiVertexzTruthPionPlus","Phi-Vertexz TruthPion +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthPionMinus=new TH2F("fHistPhiVertexzTruthPionMinus","Phi-Vertexz TruthPion -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzTruthKaon=new TH2F("fHistPhiVertexzTruthKaon","Phi-Vertexz TruthKaon;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthKaonPlus=new TH2F("fHistPhiVertexzTruthKaonPlus","Phi-Vertexz TruthKaon +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthKaonMinus=new TH2F("fHistPhiVertexzTruthKaonMinus","Phi-Vertexz TruthKaon -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzTruthProton=new TH2F("fHistPhiVertexzTruthProton","Phi-Vertexz TruthProton;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthProtonPlus=new TH2F("fHistPhiVertexzTruthProtonPlus","Phi-Vertexz TruthProton +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzTruthProtonMinus=new TH2F("fHistPhiVertexzTruthProtonMinus","Phi-Vertexz TruthProton -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    
    fHistPhiVertexzRecoPion=new TH2F("fHistPhiVertexzRecoPion","Phi-Vertexz RecoPion;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoPionPlus=new TH2F("fHistPhiVertexzRecoPionPlus","Phi-Vertexz RecoPion +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoPionMinus=new TH2F("fHistPhiVertexzRecoPionMinus","Phi-Vertexz RecoPion -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzRecoKaon=new TH2F("fHistPhiVertexzRecoKaon","Phi-Vertexz RecoKaon;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoKaonPlus=new TH2F("fHistPhiVertexzRecoKaonPlus","Phi-Vertexz RecoKaon +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoKaonMinus=new TH2F("fHistPhiVertexzRecoKaonMinus","Phi-Vertexz RecoKaon -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzRecoProton=new TH2F("fHistPhiVertexzRecoProton","Phi-Vertexz RecoProton;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoProtonPlus=new TH2F("fHistPhiVertexzRecoProtonPlus","Phi-Vertexz RecoProton +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzRecoProtonMinus=new TH2F("fHistPhiVertexzRecoProtonMinus","Phi-Vertexz RecoProton -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzContaminationPion=new TH2F("fHistPhiVertexzContaminationPion","Phi-Vertexz ContaminationPion;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationPionPlus=new TH2F("fHistPhiVertexzContaminationPionPlus","Phi-Vertexz ContaminationPion +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationPionMinus=new TH2F("fHistPhiVertexzContaminationPionMinus","Phi-Vertexz ContaminationPion -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzContaminationKaon=new TH2F("fHistPhiVertexzContaminationKaon","Phi-Vertexz ContaminationKaon;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationKaonPlus=new TH2F("fHistPhiVertexzContaminationKaonPlus","Phi-Vertexz ContaminationKaon +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationKaonMinus=new TH2F("fHistPhiVertexzContaminationKaonMinus","Phi-Vertexz ContaminationKaon -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);

    fHistPhiVertexzContaminationProton=new TH2F("fHistPhiVertexzContaminationProton","Phi-Vertexz ContaminationProton;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationProtonPlus=new TH2F("fHistPhiVertexzContaminationProtonPlus","Phi-Vertexz ContaminationProton +;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);
    fHistPhiVertexzContaminationProtonMinus=new TH2F("fHistPhiVertexzContaminationProtonMinus","Phi-Vertexz ContaminationProton -;#phi;Vz",phiBin,nArrayPhi,VertexZBin,nArrayVertexZ);


fOutputList->Add(fHistEtaVertexzTruthPion);
fOutputList->Add(fHistEtaVertexzTruthPionPlus);
fOutputList->Add(fHistEtaVertexzTruthPionMinus);

fOutputList->Add(fHistEtaVertexzTruthKaon);
fOutputList->Add(fHistEtaVertexzTruthKaonPlus);
fOutputList->Add(fHistEtaVertexzTruthKaonMinus);

fOutputList->Add(fHistEtaVertexzTruthProton);
fOutputList->Add(fHistEtaVertexzTruthProtonPlus);
fOutputList->Add(fHistEtaVertexzTruthProtonMinus);

fOutputList->Add(fHistEtaVertexzRecoPion);
fOutputList->Add(fHistEtaVertexzRecoPionPlus);
fOutputList->Add(fHistEtaVertexzRecoPionMinus);

fOutputList->Add(fHistEtaVertexzRecoKaon);
fOutputList->Add(fHistEtaVertexzRecoKaonPlus);
fOutputList->Add(fHistEtaVertexzRecoKaonMinus);

fOutputList->Add(fHistEtaVertexzRecoProton);
fOutputList->Add(fHistEtaVertexzRecoProtonPlus);
fOutputList->Add(fHistEtaVertexzRecoProtonMinus);

fOutputList->Add(fHistEtaVertexzContaminationPion);
fOutputList->Add(fHistEtaVertexzContaminationPionPlus);
fOutputList->Add(fHistEtaVertexzContaminationPionMinus);

fOutputList->Add(fHistEtaVertexzContaminationKaon);
fOutputList->Add(fHistEtaVertexzContaminationKaonPlus);
fOutputList->Add(fHistEtaVertexzContaminationKaonMinus);

fOutputList->Add(fHistEtaVertexzContaminationProton);
fOutputList->Add(fHistEtaVertexzContaminationProtonPlus);
fOutputList->Add(fHistEtaVertexzContaminationProtonMinus);

fOutputList->Add(fHistPhiVertexzTruthPion);
fOutputList->Add(fHistPhiVertexzTruthPionPlus);
fOutputList->Add(fHistPhiVertexzTruthPionMinus);

fOutputList->Add(fHistPhiVertexzTruthKaon);
fOutputList->Add(fHistPhiVertexzTruthKaonPlus);
fOutputList->Add(fHistPhiVertexzTruthKaonMinus);

fOutputList->Add(fHistPhiVertexzTruthProton);
fOutputList->Add(fHistPhiVertexzTruthProtonPlus);
fOutputList->Add(fHistPhiVertexzTruthProtonMinus);

fOutputList->Add(fHistPhiVertexzRecoPion);
fOutputList->Add(fHistPhiVertexzRecoPionPlus);
fOutputList->Add(fHistPhiVertexzRecoPionMinus);

fOutputList->Add(fHistPhiVertexzRecoKaon);
fOutputList->Add(fHistPhiVertexzRecoKaonPlus);
fOutputList->Add(fHistPhiVertexzRecoKaonMinus);

fOutputList->Add(fHistPhiVertexzRecoProton);
fOutputList->Add(fHistPhiVertexzRecoProtonPlus);
fOutputList->Add(fHistPhiVertexzRecoProtonMinus);

fOutputList->Add(fHistPhiVertexzContaminationPion);
fOutputList->Add(fHistPhiVertexzContaminationPionPlus);
fOutputList->Add(fHistPhiVertexzContaminationPionMinus);

fOutputList->Add(fHistPhiVertexzContaminationKaon);
fOutputList->Add(fHistPhiVertexzContaminationKaonPlus);
fOutputList->Add(fHistPhiVertexzContaminationKaonMinus);

fOutputList->Add(fHistPhiVertexzContaminationProton);
fOutputList->Add(fHistPhiVertexzContaminationProtonPlus);
fOutputList->Add(fHistPhiVertexzContaminationProtonMinus);
// Efficieny in Eta and Phi for Different Vertex-z cut 


  fHistdEdxTPC = new TH2F("fHistdEdxTPC", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);
  fOutputList->Add(fHistdEdxTPC);

  fHistBetaTOF = new TH2F("fHistBetaTOF", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);
  fOutputList->Add(fHistBetaTOF);

  fHistdEdxTPCPionAfterPIDCut = new TH2F("fHistdEdxTPCPionAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);
  fHistdEdxTPCKaonAfterPIDCut = new TH2F("fHistdEdxTPCKaonAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);
  fHistdEdxTPCProtonAfterPIDCut = new TH2F("fHistdEdxTPCProtonAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);

  fOutputList->Add(fHistdEdxTPCPionAfterPIDCut);
  fOutputList->Add(fHistdEdxTPCKaonAfterPIDCut);
  fOutputList->Add(fHistdEdxTPCProtonAfterPIDCut);

  fHistBetaTOFPionAfterPIDCut = new TH2F("fHistBetaTOFPionAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);
  fHistBetaTOFKaonAfterPIDCut = new TH2F("fHistBetaTOFKaonAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);
  fHistBetaTOFProtonAfterPIDCut = new TH2F("fHistBetaTOFProtonAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);

  fOutputList->Add(fHistBetaTOFPionAfterPIDCut);
  fOutputList->Add(fHistBetaTOFKaonAfterPIDCut);
  fOutputList->Add(fHistBetaTOFProtonAfterPIDCut);

 
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

  // ==================================================================================================================================================================
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
  // ==================================================================================================================================================================

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

//  Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);

    Double_t MassPion   = 0.139570; // GeV/c2
    Double_t MassKaon   = 0.493677; // GeV/c2
    Double_t MassProton = 0.938272; // GeV/c2

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
	     
//-------------------------------------------------MCAOD----------------------------------------------------------------------------------------------------------
                       Short_t vCharge;
                       Double_t vEta;
                       Double_t vPionY;
                       Double_t vKaonY;
                       Double_t vProtonY;
                       Double_t vPhi;
                       Double_t vPt;
                       Double_t vYCut;
	      
	      //for(Int_t iParticle = 0; iParticle < mcEvent->GetNumberOfTracks(); iParticle++) {
	      for(Int_t iParticle = 0; iParticle < fArrayMC->GetEntriesFast(); iParticle++) {
               
              //AliAODMCParticle* currentAODMCParticle = (AliAODMCParticle*) mcEvent->GetTrack(iParticle);
              AliAODMCParticle* currentAODMCParticle = (AliAODMCParticle*) fArrayMC->At(iParticle);

               if (!currentAODMCParticle) {
               AliError(Form("ERROR: Could not receive track %d (mc loop)", iParticle));
               continue;
               }

                if (!IsMCParticleCut(currentAODMCParticle)) continue;

                if(currentAODMCParticle->IsSecondaryFromWeakDecay()) continue;

                if(!currentAODMCParticle->IsPhysicalPrimary()) continue;

          if (fInjectedSignals && currentAODMCParticle->GetLabel() >= skipParticlesAbove) continue;


            Int_t pdgCode=((AliAODMCParticle*)currentAODMCParticle)->GetPdgCode();
            
             if (TMath::Abs(pdgCode)==11) continue;
	
            Short_t gAODmcCharge = currentAODMCParticle->Charge();

    
      vCharge = currentAODMCParticle->Charge();
      vEta    = currentAODMCParticle->Eta();
      vPhi    = currentAODMCParticle->Phi();// * TMath::RadToDeg();
      vPt     = currentAODMCParticle->Pt();


//cout<<"MCAOD ---- Eta : "<<vEta<<'\t'<<"Phi :"<<vPhi<<'\t'<<"Pt:"<<vPt<<endl;


     vPionY = log( ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(MassPion*MassPion + vPt*vPt) ); // convert eta to y
     vKaonY = log( ( sqrt(MassKaon*MassKaon + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(MassKaon*MassKaon + vPt*vPt) ); // convert eta to y
     vProtonY = log( ( sqrt(MassProton*MassProton + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(MassProton*MassProton + vPt*vPt) ); // convert eta to y

     if(fRapidityInsteadOfEta){


           if(TMath::Abs(pdgCode) == 211) vYCut=vPionY;
           else if(TMath::Abs(pdgCode) == 321) vYCut=vKaonY;
           else if(TMath::Abs(pdgCode) == 2212) vYCut=vProtonY;
          else {
           continue;
          }
 

      if(vYCut < fMinEta || vYCut > fMaxEta) continue;
 }
           if(TMath::Abs(pdgCode) == 211){
//          cout<<"test for Pion at CERN r1"<<endl;
           fZvertexTruthPion->Fill(vertex->GetZ());
           fHistPhiVertexzTruthPion->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthPion->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthPion->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthPion->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPion->Fill(vEta,vertex->GetZ());
           }

           if(gAODmcCharge >0){
           fZvertexTruthPionPlus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthPionPlus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthPionPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthPionPlus->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthPionPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPionPlus->Fill(vEta,vertex->GetZ());
           }
            }

           if(gAODmcCharge < 0){
           fZvertexTruthPionMinus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthPionMinus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta) {
          fHistTruthPionMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthPionMinus->Fill(vYCut,vertex->GetZ());
          }
           else {
           fHistTruthPionMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPionMinus->Fill(vEta,vertex->GetZ());
           }
           }
        
         }        

   
         else if(TMath::Abs(pdgCode) ==321){
           
           fZvertexTruthKaon->Fill(vertex->GetZ());
           fHistPhiVertexzTruthKaon->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthKaon->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthKaon->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthKaon->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaon->Fill(vEta,vertex->GetZ());
           }

           if(gAODmcCharge >0){
           fZvertexTruthKaonPlus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthKaonPlus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthKaonPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthKaonPlus->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthKaonPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaonPlus->Fill(vEta,vertex->GetZ());
           }
            }

           if(gAODmcCharge < 0){
           fZvertexTruthKaonMinus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthKaonMinus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta) {
          fHistTruthKaonMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthKaonMinus->Fill(vYCut,vertex->GetZ());
          }
           else {
           fHistTruthKaonMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaonMinus->Fill(vEta,vertex->GetZ());
           }
           }
        
        }

        else if(TMath::Abs(pdgCode) == 2212){

           fZvertexTruthProton->Fill(vertex->GetZ());
           fHistPhiVertexzTruthProton->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthProton->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthProton->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthProton->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProton->Fill(vEta,vertex->GetZ());
           }

           if(gAODmcCharge >0){
           fZvertexTruthProtonPlus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthProtonPlus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta){
           fHistTruthProtonPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthProtonPlus->Fill(vYCut,vertex->GetZ());
           }
           else {
           fHistTruthProtonPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProtonPlus->Fill(vEta,vertex->GetZ());
           }
            }

           if(gAODmcCharge < 0){
           fZvertexTruthProtonMinus->Fill(vertex->GetZ());
           fHistPhiVertexzTruthProtonMinus->Fill(vPhi,vertex->GetZ());
           if(fRapidityInsteadOfEta) {
          fHistTruthProtonMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthProtonMinus->Fill(vYCut,vertex->GetZ());
          }
           else {
           fHistTruthProtonMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProtonMinus->Fill(vEta,vertex->GetZ());
           }
           }
        
         } 

	
}//loop over tracks

	      
//-------------------------------------------------MCAODrec----------------------------------------------------------------------------------------------------------
	  
                       Short_t vChargeReco;
                       Double_t vEtaReco;
                       Double_t vPionYReco;
                       Double_t vKaonYReco;
                       Double_t vProtonYReco;
                       Double_t vPhiReco;
                       Double_t vPtReco;
                       Double_t vYReco;
   
	      Int_t nGoodTracks = fAOD->GetNumberOfTracks();   
	     
	      
	      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));    
		if(!trackAOD) continue;
		

		//track cuts
		if (!trackAOD->TestFilterBit(fAODTrackCutBit)) 
		  continue;

 		Int_t label = TMath::Abs(trackAOD->GetLabel()); 
		if(label > trackAOD->GetLabel()) continue; 
		
		  if(trackAOD->Eta() <fMinEta || trackAOD->Eta() > fMaxEta) 
		    continue;

		  if((trackAOD->Pt() > fMaxPt)||(trackAOD->Pt() <  fMinPt)) 
		    continue;

// Remove Electron Track : 

        AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);

        if (AODmcTrack){
          if(TMath::Abs(AODmcTrack->GetPdgCode()) == 11) continue;
        }

                AliAODMCParticle* recoMC = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(trackAOD->GetLabel())));

               if(!recoMC) continue;

               if (((AliAODMCParticle*) recoMC)->IsSecondaryFromWeakDecay()) continue;

               if (!recoMC->IsPhysicalPrimary()) continue;
		  
		  Short_t gCharge = trackAOD->Charge();
		  Double_t phiRad = trackAOD->Phi();

// Electron Rejection -------------------------------------------------------------------------------------------------------------------------------------------------		  

                   if(fElectronRejection) {

		    
		    // get the electron nsigma
		    Double_t nSigma = fPIDResponse->NumberOfSigmasTPC(trackAOD,(AliPID::EParticleType)AliPID::kElectron);
		    fHistNSigmaTPCvsPtbeforeElectronPID->Fill(trackAOD->Pt(),nSigma);		    
		    
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
		    
		    fHistNSigmaTPCvsPtafterElectronPID->Fill(trackAOD->Pt(),nSigma);		    

		  }
		  
               Int_t pdgCodeReco = ((AliAODMCParticle*)recoMC)->GetPdgCode();


      vChargeReco = trackAOD->Charge();
      vEtaReco    = trackAOD->Eta();
      vPhiReco    = trackAOD->Phi();// * TMath::RadToDeg();
      vPtReco     = trackAOD->Pt();

     vPionYReco = log( ( sqrt(MassPion*MassPion + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) / sqrt(MassPion*MassPion + vPtReco*vPtReco) ); // convert eta to y
     vKaonYReco = log( ( sqrt(MassKaon*MassKaon + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) / sqrt(MassKaon*MassKaon + vPtReco*vPtReco) ); // convert eta to y
     vProtonYReco = log( ( sqrt(MassProton*MassProton + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) / sqrt(MassProton*MassProton + vPtReco*vPtReco) ); // convert eta to y


//cout<<"MCAODrec ---- Eta : "<<vEtaReco<<'\t'<<"Phi :"<<vPhiReco<<'\t'<<"Pt:"<<vPtReco<<endl;

//----------------------------------------------------------------------Nsigma Method for PID--------------------------------------------------------------------------
Double_t dEdx=-1.0;
Double_t beta=-1.0;

Double_t nsigmaSpecies[6]={999.0};

Double_t nsigmaTPC[6] = {999.0};

Double_t nsigmaTOF[6]={999.0};
Double_t nsigmaTPCTOF[6]={999.0};

for(int iSpecies=0;iSpecies<6;iSpecies++){
nsigmaTPCTOF[iSpecies]  =GetNsigmas(fPIDResponse ,trackAOD,iSpecies);
}

for(int iSpecies=0;iSpecies<6;iSpecies++){

// TPC : Nsigma 
  nsigmaTPC[iSpecies] = fNsigmaTPC[iSpecies];

// TOF : Nsigma
  nsigmaTOF[iSpecies] = fNsigmaTOF[iSpecies];

}

Bool_t IsTPCSignal=kFALSE;
Bool_t IsTOFSignal=kFALSE;

if(trackAOD->GetTPCsignal() >0.0) IsTPCSignal=kTRUE;
if(IsTOFPID(trackAOD) && (Beta(trackAOD) >0.0 && Beta(trackAOD)<=1.0)) IsTOFSignal=kTRUE;



if(IsTPCSignal){
dEdx   = trackAOD->GetTPCsignal();  //dEdX for TPC
fHistdEdxTPC->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);
}


if(IsTOFSignal){
 beta = Beta(trackAOD); // Beta for TOF 
 fHistBetaTOF->Fill(trackAOD->Pt()*trackAOD->Charge(), beta);
}



Double_t MisMatchTOFProb = fPIDResponse->GetTOFMismatchProbability(trackAOD);

if(fTOFMisMatch){


if(MisMatchTOFProb < fMistMatchTOFProb){

if(fDetectorPID_== kTPCTOFpid_){

if(IsTOFSignal){
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
}
else {
continue;
} 
fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);

fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);

}

else if(fDetectorPID_ == kTogether_){
if(IsTOFSignal){
if(trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);

}
}
else if(IsTPCSignal) {
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fPtTPCMax){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTPC[iSpecies]);

fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);

}
}

else {
continue;
}

fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);

}


}
}


else {

if(fDetectorPID_== kTPCTOFpid_){
if(IsTOFSignal){
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
}
else {
continue;
} 
fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);

fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);

}

else if(fDetectorPID_ == kTogether_){
if(IsTOFSignal){
if(trackAOD->Pt()>fPtTPCMax && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);

}
}
else if(IsTPCSignal) {
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fPtTPCMax){
for(int iSpecies=0;iSpecies<6;iSpecies++)
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTPC[iSpecies]);

fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);

}
}

else {
continue;
}

fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);

}

} // end of else 


Int_t MostProbableSpecies = TMath::LocMin(6, nsigmaSpecies);

if(TMath::Abs(nsigmaSpecies[MostProbableSpecies]) > fNSigmaPID) continue ;

if(MostProbableSpecies == 0 || MostProbableSpecies == 1 || MostProbableSpecies == 5 ) continue;


if(MostProbableSpecies == 2) vYReco=vPionYReco;
else if(MostProbableSpecies == 3) vYReco=vKaonYReco;
else if(MostProbableSpecies == 4) vYReco=vProtonYReco;

if(fRapidityInsteadOfEta) {

      if( vYReco < fMinEta || vYReco > fMaxEta)  continue;
}

//Pion

if(TMath::Abs(pdgCodeReco) == 211 && MostProbableSpecies == 2){ 

fHistdEdxTPCPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFPionAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[2]);

fZvertexRecoPion->Fill(vertex->GetZ());
fHistPhiVertexzRecoPion->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexRecoPionPlus->Fill(vertex->GetZ());
fHistPhiVertexzRecoPionPlus->Fill(vPhiReco,vertex->GetZ());
}
if(gCharge<0) {
fZvertexRecoPionMinus->Fill(vertex->GetZ());
fHistPhiVertexzRecoPionMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){
fHistMCRecoPion->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPion->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoPionPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionPlus->Fill(vYReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoPionMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionMinus->Fill(vYReco,vertex->GetZ());
}
}

else{
fHistMCRecoPion->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPion->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoPionPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionPlus->Fill(vEtaReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoPionMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionMinus->Fill(vEtaReco,vertex->GetZ());
}
}

}

//Kaon

if(TMath::Abs(pdgCodeReco) == 321 && MostProbableSpecies == 3){ 

fHistdEdxTPCKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFKaonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[3]);

fZvertexRecoKaon->Fill(vertex->GetZ());
fHistPhiVertexzRecoKaon->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexRecoKaonPlus->Fill(vertex->GetZ());
fHistPhiVertexzRecoKaonPlus->Fill(vPhiReco,vertex->GetZ());
}
if(gCharge<0) {
fZvertexRecoKaonMinus->Fill(vertex->GetZ());
fHistPhiVertexzRecoKaonMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){
fHistMCRecoKaon->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaon->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoKaonPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonPlus->Fill(vYReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoKaonMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonMinus->Fill(vYReco,vertex->GetZ());
}
}

else{
fHistMCRecoKaon->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaon->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoKaonPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonPlus->Fill(vEtaReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoKaonMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonMinus->Fill(vEtaReco,vertex->GetZ());
}
}

}

// Proton

if(TMath::Abs(pdgCodeReco) == 2212 && MostProbableSpecies == 4){ 

fHistdEdxTPCProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFProtonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[4]);

fZvertexRecoProton->Fill(vertex->GetZ());
fHistPhiVertexzRecoProton->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexRecoProtonPlus->Fill(vertex->GetZ());
fHistPhiVertexzRecoProtonPlus->Fill(vPhiReco,vertex->GetZ());
}
if(gCharge<0) {
fZvertexRecoProtonMinus->Fill(vertex->GetZ());
fHistPhiVertexzRecoProtonMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){
fHistMCRecoProton->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProton->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoProtonPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonPlus->Fill(vYReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoProtonMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonMinus->Fill(vYReco,vertex->GetZ());
}
}

else{
fHistMCRecoProton->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProton->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
fHistMCRecoProtonPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonPlus->Fill(vEtaReco,vertex->GetZ());
}
if(gCharge<0) {
fHistMCRecoProtonMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonMinus->Fill(vEtaReco,vertex->GetZ());
}
}

}

// Fill the Histograms for Contamination
if(TMath::Abs(pdgCodeReco) !=211 && MostProbableSpecies == 2) {

fZvertexContaminationPion->Fill(vertex->GetZ());
fHistPhiVertexzContaminationPion->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexContaminationPionPlus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationPionPlus->Fill(vPhiReco,vertex->GetZ());
}

if(gCharge<0) {
fZvertexContaminationPionMinus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationPionMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){

Hist3dPionContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPion->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
Hist3dPionPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionPlus->Fill(vYReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dPionMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionMinus->Fill(vYReco,vertex->GetZ());
}
}

else {
Hist3dPionContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPion->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
Hist3dPionPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionPlus->Fill(vEtaReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dPionMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionMinus->Fill(vEtaReco,vertex->GetZ());
}
}

}


if(TMath::Abs(pdgCodeReco) !=321 && MostProbableSpecies == 3) {

fZvertexContaminationKaon->Fill(vertex->GetZ());
fHistPhiVertexzContaminationKaon->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexContaminationKaonPlus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationKaonPlus->Fill(vPhiReco,vertex->GetZ());
}

if(gCharge<0) {
fZvertexContaminationKaonMinus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationKaonMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){

Hist3dKaonContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaon->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
Hist3dKaonPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonPlus->Fill(vYReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dKaonMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonMinus->Fill(vYReco,vertex->GetZ());
}
}

else {
Hist3dKaonContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaon->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
Hist3dKaonPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonPlus->Fill(vEtaReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dKaonMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonMinus->Fill(vEtaReco,vertex->GetZ());
}
}

}


if(TMath::Abs(pdgCodeReco) !=2212 && MostProbableSpecies == 4) {

fZvertexContaminationProton->Fill(vertex->GetZ());
fHistPhiVertexzContaminationProton->Fill(vPhiReco,vertex->GetZ());

if(gCharge>0) {
fZvertexContaminationProtonPlus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationProtonPlus->Fill(vPhiReco,vertex->GetZ());
}

if(gCharge<0) {
fZvertexContaminationProtonMinus->Fill(vertex->GetZ());
fHistPhiVertexzContaminationProtonMinus->Fill(vPhiReco,vertex->GetZ());
}

if(fRapidityInsteadOfEta){

Hist3dProtonContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProton->Fill(vYReco,vertex->GetZ());
if(gCharge>0) {
Hist3dProtonPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonPlus->Fill(vYReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dProtonMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonMinus->Fill(vYReco,vertex->GetZ());
}
}

else {
Hist3dProtonContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProton->Fill(vEtaReco,vertex->GetZ());
if(gCharge>0) {
Hist3dProtonPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonPlus->Fill(vEtaReco,vertex->GetZ());
}

if(gCharge<0) {
Hist3dProtonMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonMinus->Fill(vEtaReco,vertex->GetZ());
}
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

       
        if (particle->Eta() < fMinEta || particle->Eta() > fMaxEta) {return kFALSE;}

        return kTRUE;
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
  AliPIDResponse::EDetPidStatus statusTOF = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track);
  if(statusTOF != AliPIDResponse::kDetPidOk)
    return kFALSE;
  return kTRUE;
}


Double_t AliAnalysisTaskEffContPIDBF::GetNsigmas(AliPIDResponse* PIDresponse , AliAODTrack* track , Int_t specie)
{
    Double_t nsigmaTPC = 999.0;
    Double_t nsigmaTOF = 999.0;

    AliPIDResponse::EDetPidStatus statusTPC = PIDresponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) specie, nsigmaTPC);
    AliPIDResponse::EDetPidStatus statusTOF = PIDresponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) specie, nsigmaTOF);
    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);/* && trk->IsOn(AliESDtrack::kTPCpid)*/;
    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);

   Bool_t fHasTPCPID=kFALSE;
   Bool_t fHasTOFPID=kFALSE;

  if(track->GetTPCsignal() >0.0)
   fHasTPCPID=kTRUE;


   if(tofIsOk) {
   if(Beta(track) >0.0 && Beta(track)<=1.0)
   fHasTOFPID=kTRUE;
   }

    fNsigmaTPC[specie] = fHasTPCPID? nsigmaTPC: 999.0;
    fNsigmaTOF[specie] = fHasTOFPID? nsigmaTOF: 999.0;
//    if(fHasTOFPID && fNsigmaTOF[specie] == -998) fNsigmaTOF[specie]=999.0;

//cout<<"nsigma TOF :"<<fNsigmaTOF[specie]<<'\t'<<"nsigma TPC: "<<fNsigmaTPC[specie]<<'\t'<<"nsigma TPCTOF:"<<TMath::Hypot(fNsigmaTPC[specie],fNsigmaTOF[specie])<<endl;
    return TMath::Hypot(fNsigmaTPC[specie], fNsigmaTOF[specie]);
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
