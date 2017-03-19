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
#include "AliVEvent.h" 
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
    fQAListTruthReco(0),
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
    fHistPhiVertexzContaminationProtonMinus(0),
    fHistNsigmaTPCPionAfterPIDCut(0),
    fHistNsigmaTPCKaonAfterPIDCut(0),
    fHistNsigmaTPCProtonAfterPIDCut(0),
    fHistNsigmaTOFPionAfterPIDCut(0),
    fHistNsigmaTOFKaonAfterPIDCut(0),
    fHistNsigmaTOFProtonAfterPIDCut(0),
    fHistEtaMCPion(0),
    fHistEtaMCKaon(0),
    fHistEtaMCProton(0),
    fHistPhiMCPion(0),
    fHistEtaMCAll(0),
    fHistPhiMCAll(0),
    fHistPhiMCKaon(0),
    fHistPhiMCProton(0),
    fHistEtaMCRecoPion(0),
    fHistEtaMCRecoKaon(0),
    fHistEtaMCRecoProton(0),
    fHistPhiMCRecoPion(0),
    fHistPhiMCRecoKaon(0),
    fHistPhiMCRecoProton(0),
    fHistPtMCProton(0),
    fHistPhiMCProtonTruthPlus(0),
    fHistPtMCProtonTruthPlus(0),
    fHistEtaMCProtonTruthPlus(0),
    fHistPhiMCProtonTruthMinus(0),
    fHistPtMCProtonTruthMinus(0),
    fHistEtaMCProtonTruthMinus(0),
    fHistPtMCPion(0),
    fHistPhiMCPionTruthPlus(0),
    fHistPtMCPionTruthPlus(0),
    fHistEtaMCPionTruthPlus(0),
    fHistPhiMCPionTruthMinus(0),
    fHistPtMCPionTruthMinus(0),
    fHistEtaMCPionTruthMinus(0),
    fHistPtMCKaon(0),
    fHistPhiMCKaonTruthPlus(0),
    fHistPtMCKaonTruthPlus(0),
    fHistEtaMCKaonTruthPlus(0),
    fHistPhiMCKaonTruthMinus(0),
    fHistPtMCKaonTruthMinus(0),
    fHistEtaMCKaonTruthMinus(0),
    fHistPtMCRecoPion(0),
    fHistPtMCRecoPionPlus(0),
    fHistPhiMCRecoPionPlus(0),
    fHistPtMCRecoPionMinus(0),
    fHistPhiMCRecoPionMinus(0),
    fHistEtaMCRecoPionPlus(0),
    fHistEtaMCRecoPionMinus(0),
    fHistPtMCRecoKaon(0),
    fHistPtMCRecoKaonPlus(0),
    fHistPhiMCRecoKaonPlus(0),
    fHistPtMCRecoKaonMinus(0),
    fHistPhiMCRecoKaonMinus(0),
    fHistEtaMCRecoKaonPlus(0),
    fHistEtaMCRecoKaonMinus(0),
    fHistPtMCRecoProton(0),
    fHistPtMCRecoProtonPlus(0),
    fHistPhiMCRecoProtonPlus(0),
    fHistPtMCRecoProtonMinus(0),
    fHistPhiMCRecoProtonMinus(0),
    fHistEtaMCRecoProtonPlus(0),
    fHistEtaMCRecoProtonMinus(0),
    fHistPhiContaminationPion(0),
    fHistEtaContaminationPion(0),
    fHistPtContaminationPion(0),
    fHistPhiContaminationPionPlus(0),
    fHistEtaContaminationPionPlus(0),
    fHistPtContaminationPionPlus(0),
    fHistPhiContaminationPionMinus(0),
    fHistEtaContaminationPionMinus(0),
    fHistPtContaminationPionMinus(0),
    fHistPhiContaminationKaon(0),
    fHistEtaContaminationKaon(0),
    fHistPtContaminationKaon(0),
    fHistPhiContaminationKaonPlus(0),
    fHistEtaContaminationKaonPlus(0),
    fHistPtContaminationKaonPlus(0),
    fHistPhiContaminationKaonMinus(0),
    fHistEtaContaminationKaonMinus(0),
    fHistPtContaminationKaonMinus(0),
    fHistPhiContaminationProton(0),
    fHistEtaContaminationProton(0),
    fHistPtContaminationProton(0),
    fHistPhiContaminationProtonPlus(0),
    fHistEtaContaminationProtonPlus(0),
    fHistPtContaminationProtonPlus(0),
    fHistPhiContaminationProtonMinus(0),
    fHistEtaContaminationProtonMinus(0),
    fHistPtContaminationProtonMinus(0),
    fVertexZ(0),
    fDCAxyCut(-1),
    fDCAzCut(-1),
    fTPCchi2Cut(-1),
    fNClustersTPCCut(-1)

/*  fHistEtaPtPhiVertxezTruthPion(0),
 fHistEtaPtPhiVertxezTruthPionPlus(0),
 fHistEtaPtPhiVertxezTruthPionMinus(0),
 fHistEtaPtPhiVertxezTruthKaon(0),
 fHistEtaPtPhiVertxezTruthKaonPlus(0),
 fHistEtaPtPhiVertxezTruthKaonMinus(0),
 fHistEtaPtPhiVertxezTruthProton(0),
 fHistEtaPtPhiVertxezTruthProtonPlus(0),
 fHistEtaPtPhiVertxezTruthProtonMinus(0),
 fHistEtaPtPhiVertxezContaminationPion(0),
 fHistEtaPtPhiVertxezContaminationPionPlus(0),
 fHistEtaPtPhiVertxezContaminationPionMinus(0),
 fHistEtaPtPhiVertxezContaminationKaon(0),
 fHistEtaPtPhiVertxezContaminationKaonPlus(0),
 fHistEtaPtPhiVertxezContaminationKaonMinus(0),
 fHistEtaPtPhiVertxezContaminationProton(0),
 fHistEtaPtPhiVertxezContaminationProtonPlus(0),
 fHistEtaPtPhiVertxezContaminationProtonMinus(0),
 fHistEtaPtPhiVertxezRecoPion(0),
 fHistEtaPtPhiVertxezRecoPionPlus(0),
 fHistEtaPtPhiVertxezRecoPionMinus(0),
 fHistEtaPtPhiVertxezRecoKaon(0),
 fHistEtaPtPhiVertxezRecoKaonPlus(0),
 fHistEtaPtPhiVertxezRecoKaonMinus(0),
 fHistEtaPtPhiVertxezRecoProton(0),
 fHistEtaPtPhiVertxezRecoProtonPlus(0),
 fHistEtaPtPhiVertxezRecoProtonMinus(0)*/                                                 
{ 
    for(Int_t ipart=0;ipart<5;ipart++){
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
  DefineOutput(3, TList::Class());
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
 

  fQAListTruthReco = new TList();
  fQAListTruthReco->SetName("QAListTruthReco");
  fQAListTruthReco->SetOwner();
 
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
  
  Int_t ptBin = 36;
 // Int_t ptBin = 44;
  Int_t etaBin = 16;
  Int_t phiBin = 100;

  Int_t VertexZBin=10;

   Double_t nArrayPt[37]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};

//   Double_t nArrayPt[45]={0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.95,2.0,2.05,2.1,2.15,2.2,2.25,2.3,2.35,2.4};

  Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Double_t nArrayPhi[phiBin+1];
  for(Int_t iBin = 0; iBin <= phiBin; iBin++) 
  nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;


//  Double_t nArrayVertexZ[11]={-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
    Double_t nArrayVertexZ[25]={-6,-5.5,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6};

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

    fHistEtaMCPion= new TH1F("fHistEtaMCPion","Eta distribution MC Pion",200,-2,2);
    fHistEtaMCKaon= new TH1F("fHistEtaMCKaon","Eta distribution MC Kaon",200,-2,2);
    fHistEtaMCProton= new TH1F("fHistEtaMCProton","Eta distribution MC Proton",200,-2,2);
    fHistEtaMCAll =new TH1F("fHistEtaMCAll","Eta distribution MC All",200,-2,2);
    fHistPhiMCAll =new TH1F("fHistPhiMCAll","Phi distribution MC All",200,-1,2.5*TMath::Pi());

    fHistPhiMCPion= new TH1F("fHistPhiMCPion","Phi distribution MC Pion",200,-1,2.5*TMath::Pi());
    fHistPhiMCKaon= new TH1F("fHistPhiMCKaon","Phi distribution MC Kaon",200,-1,2.5*TMath::Pi());
    fHistPhiMCProton= new TH1F("fHistPhiMCProton","Phi distribution MC Proton",200,-1,2.5*TMath::Pi());
    
    fHistEtaMCRecoPion= new TH1F("fHistEtaMCRecoPion","Eta distribution MCReco Pion",200,-2,2);
    fHistEtaMCRecoKaon= new TH1F("fHistEtaMCRecoKaon","Eta distribution MCReco Kaon",200,-2,2);
    fHistEtaMCRecoProton= new TH1F("fHistEtaMCRecoProton","Eta distribution MCReco Proton",200,-2,2);
    
    fHistPhiMCRecoPion= new TH1F("fHistPhiMCRecoPion","Phi distribution MCReco Pion",200,-1,2.5*TMath::Pi());
    fHistPhiMCRecoKaon= new TH1F("fHistPhiMCRecoKaon","Phi distribution MCReco Kaon",200,-1,2.5*TMath::Pi());
    fHistPhiMCRecoProton= new TH1F("fHistPhiMCRecoProton","Phi distribution MCReco Proton",200,-1,2.5*TMath::Pi());

    fQAListTruthReco->Add(fHistEtaMCPion);
    fQAListTruthReco->Add(fHistEtaMCKaon);
    fQAListTruthReco->Add(fHistEtaMCProton);
    fQAListTruthReco->Add(fHistPhiMCPion);
    fQAListTruthReco->Add(fHistPhiMCKaon);
    fQAListTruthReco->Add(fHistPhiMCProton);
    fQAListTruthReco->Add(fHistEtaMCRecoPion);
    fQAListTruthReco->Add(fHistEtaMCRecoKaon);
    fQAListTruthReco->Add(fHistEtaMCRecoProton);
    fQAListTruthReco->Add(fHistPhiMCRecoPion);
    fQAListTruthReco->Add(fHistPhiMCRecoKaon);
    fQAListTruthReco->Add(fHistPhiMCRecoProton);


   fQAListTruthReco->Add(fHistEtaMCAll);
   fQAListTruthReco->Add(fHistPhiMCAll);

    fHistPtMCPion=new TH1F("fHistPtMCPion","Pt distribution MC Pion Truth",200,0,2);

    fHistPtMCPionTruthPlus= new TH1F("fHistPtMCPionTruthPlus","Pt distribution MC Pion Truth Plus",200,0,2);
    fHistPtMCPionTruthMinus =new TH1F("fHistPtMCPionTruthMinus","Pt distribution MC Pion Truth Minus",200,0,2);
    fHistPhiMCPionTruthPlus = new TH1F("fHistPhiMCPionTruthPlus","Phi distribution MC Pion Truth Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCPionTruthMinus = new TH1F("fHistPhiMCPionTruthMinus","Phi distribution MC Pion Truth Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCPionTruthPlus = new TH1F("fHistEtaMCPionTruthPlus","Eta distribution MC Pion Truth Plus",200,-2,2);
    fHistEtaMCPionTruthMinus = new TH1F("fHistEtaMCPionTruthMinus","Eta distribution MC Pion Truth Minus",200,-2,2);


    fHistPtMCProton=new TH1F("fHistPtMCProton","Pt distribution MC Proton Truth",200,0,2);

    fHistPtMCProtonTruthPlus= new TH1F("fHistPtMCProtonTruthPlus","Pt distribution MC Proton Truth Plus",200,0,2);
    fHistPtMCProtonTruthMinus= new TH1F("fHistPtMCProtonTruthMinus","Pt distribution MC Proton Truth Minus",200,0,2);
    fHistPhiMCProtonTruthPlus= new TH1F("fHistPhiMCProtonTruthPlus","Phi distribution MC Proton Truth Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCProtonTruthMinus= new TH1F("fHistPhiMCProtonTruthMinus","Phi distribution MC Proton Truth Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCProtonTruthPlus= new TH1F("fHistEtaMCProtonTruthPlus","Eta distribution MC Proton Truth Plus",200,-2,2);
    fHistEtaMCProtonTruthMinus= new TH1F("fHistEtaMCProtonTruthMinus","Eta distribution MC Proton Truth Minus",200,-2,2);

    fHistPtMCKaon=new TH1F("fHistPtMCKaon","Pt distribution MC Kaon Truth",200,0,2);

    fHistPtMCKaonTruthPlus= new TH1F("fHistPtMCKaonTruthPlus","Pt distribution MC Kaon Truth Plus",200,0,2);
    fHistPtMCKaonTruthMinus =new TH1F("fHistPtMCKaonTruthMinus","Pt distribution MC Kaon Truth Minus",200,0,2);
    fHistPhiMCKaonTruthPlus= new TH1F("fHistPhiMCKaonTruthPlus","Phi distribution MC Kaon Truth Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCKaonTruthMinus= new TH1F("fHistPhiMCKaonTruthMinus","Phi distribution MC Kaon Truth Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCKaonTruthPlus= new TH1F("fHistEtaMCKaonTruthPlus","Eta distribution MC Kaon Truth Plus",200,-2,2);
    fHistEtaMCKaonTruthMinus= new TH1F("fHistEtaMCKaonTruthMinus","Eta distribution MC Kaon Truth Minus",200,-2,2);


    fQAListTruthReco ->Add(fHistPtMCProton);
    fQAListTruthReco ->Add(fHistPhiMCProtonTruthPlus);
    fQAListTruthReco ->Add(fHistPtMCProtonTruthPlus);
    fQAListTruthReco ->Add(fHistEtaMCProtonTruthPlus);
    fQAListTruthReco ->Add(fHistPhiMCProtonTruthMinus);
    fQAListTruthReco ->Add(fHistPtMCProtonTruthMinus);
    fQAListTruthReco ->Add(fHistEtaMCProtonTruthMinus);
    fQAListTruthReco ->Add(fHistPtMCPion);
    fQAListTruthReco ->Add(fHistPhiMCPionTruthPlus);
    fQAListTruthReco ->Add(fHistPtMCPionTruthPlus);
    fQAListTruthReco ->Add(fHistEtaMCPionTruthPlus);
    fQAListTruthReco ->Add(fHistPhiMCPionTruthMinus);
    fQAListTruthReco ->Add(fHistPtMCPionTruthMinus);
    fQAListTruthReco ->Add(fHistEtaMCPionTruthMinus);
    fQAListTruthReco ->Add(fHistPtMCKaon);
    fQAListTruthReco ->Add(fHistPhiMCKaonTruthPlus);
    fQAListTruthReco ->Add(fHistPtMCKaonTruthPlus);
    fQAListTruthReco ->Add(fHistEtaMCKaonTruthPlus);
    fQAListTruthReco ->Add(fHistPhiMCKaonTruthMinus);
    fQAListTruthReco ->Add(fHistPtMCKaonTruthMinus);
    fQAListTruthReco ->Add(fHistEtaMCKaonTruthMinus);

    fHistPtMCRecoPion=new TH1F("fHistPtMCRecoPion","Pt distribution MC Reco Pion",200,0,2);
    fHistPtMCRecoPionPlus=new TH1F("fHistPtMCRecoPionPlus","Pt distribution MC Reco Pion Plus",200,0,2);
    fHistPtMCRecoPionMinus=new TH1F("fHistPtMCRecoPionMinus","Pt distribution MC Reco Pion Minus",200,0,2);
    fHistPhiMCRecoPionPlus =new TH1F("fHistPhiMCRecoPionPlus","Phi distribution MC Reco Pion Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCRecoPionMinus =new TH1F("fHistPhiMCRecoPionMinus","Phi distribution MC Reco Pion Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCRecoPionPlus =new TH1F("fHistEtaMCRecoPionPlus","Eta distribution MC Reco Pion Plus",200,-2,2);
    fHistEtaMCRecoPionMinus =new TH1F("fHistEtaMCRecoPionMinus","Eta distribution MC Reco Pion Minus",200,-2,2);

    fHistPtMCRecoKaon=new TH1F("fHistPtMCRecoKaon","Pt distribution MC Reco Kaon",200,0,2);
    fHistPtMCRecoKaonPlus=new TH1F("fHistPtMCRecoKaonPlus","Pt distribution MC Reco Kaon Plus",200,0,2);
    fHistPtMCRecoKaonMinus=new TH1F("fHistPtMCRecoKaonMinus","Pt distribution MC Reco Kaon Minus",200,0,2);
    fHistPhiMCRecoKaonPlus =new TH1F("fHistPhiMCRecoKaonPlus","Phi distribution MC Reco Kaon Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCRecoKaonMinus =new TH1F("fHistPhiMCRecoKaonMinus","Phi distribution MC Reco Kaon Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCRecoKaonPlus =new TH1F("fHistEtaMCRecoKaonPlus","Eta distribution MC Reco Kaon Plus",200,-2,2);
    fHistEtaMCRecoKaonMinus =new TH1F("fHistEtaMCRecoKaonMinus","Eta distribution MC Reco Kaon Minus",200,-2,2);


    fHistPtMCRecoProton=new TH1F("fHistPtMCRecoProton","Pt distribution MC Reco Proton",200,0,2);
    fHistPtMCRecoProtonPlus=new TH1F("fHistPtMCRecoProtonPlus","Pt distribution MC Reco Proton Plus",200,0,2);
    fHistPtMCRecoProtonMinus=new TH1F("fHistPtMCRecoProtonMinus","Pt distribution MC Reco Proton Minus",200,0,2);
    fHistPhiMCRecoProtonPlus =new TH1F("fHistPhiMCRecoProtonPlus","Phi distribution MC Reco Proton Plus",200,-1,2.5*TMath::Pi());
    fHistPhiMCRecoProtonMinus =new TH1F("fHistPhiMCRecoProtonMinus","Phi distribution MC Reco Proton Minus",200,-1,2.5*TMath::Pi());
    fHistEtaMCRecoProtonPlus =new TH1F("fHistEtaMCRecoProtonPlus","Eta distribution MC Reco Proton Plus",200,-2,2);
    fHistEtaMCRecoProtonMinus =new TH1F("fHistEtaMCRecoProtonMinus","Eta distribution MC Reco Proton Minus",200,-2,2);


    fQAListTruthReco->Add(fHistPtMCRecoPion);
    fQAListTruthReco->Add(fHistPtMCRecoPionPlus);
    fQAListTruthReco->Add(fHistPhiMCRecoPionPlus);
    fQAListTruthReco->Add(fHistPtMCRecoPionMinus);
    fQAListTruthReco->Add(fHistPhiMCRecoPionMinus);
    fQAListTruthReco->Add(fHistEtaMCRecoPionPlus);
    fQAListTruthReco->Add(fHistEtaMCRecoPionMinus);

    fQAListTruthReco->Add(fHistPtMCRecoKaon);
    fQAListTruthReco->Add(fHistPtMCRecoKaonPlus);
    fQAListTruthReco->Add(fHistPhiMCRecoKaonPlus);
    fQAListTruthReco->Add(fHistPtMCRecoKaonMinus);
    fQAListTruthReco->Add(fHistPhiMCRecoKaonMinus);
    fQAListTruthReco->Add(fHistEtaMCRecoKaonPlus);
    fQAListTruthReco->Add(fHistEtaMCRecoKaonMinus);

    fQAListTruthReco->Add(fHistPtMCRecoProton);
    fQAListTruthReco->Add(fHistPtMCRecoProtonPlus);
    fQAListTruthReco->Add(fHistPhiMCRecoProtonPlus);
    fQAListTruthReco->Add(fHistPtMCRecoProtonMinus);
    fQAListTruthReco->Add(fHistPhiMCRecoProtonMinus);
    fQAListTruthReco->Add(fHistEtaMCRecoProtonPlus);
    fQAListTruthReco->Add(fHistEtaMCRecoProtonMinus);

    fHistPhiContaminationPion=new TH1F("fHistPhiContaminationPion","Phi distribution of Contamination in Pion",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationPionPlus=new TH1F("fHistPhiContaminationPionPlus","Phi distribution of Contamination in Pion Plus",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationPionMinus=new TH1F("fHistPhiContaminationPionMinus","Phi distribution of Contamination in Pion Minus",200,-1,2.5*TMath::Pi());
    fHistEtaContaminationPion=new TH1F("fHistEtaContaminationPion","Eta distribution of Contamination in Pion",200,-2,2);
    fHistEtaContaminationPionPlus=new TH1F("fHistEtaContaminationPionPlus","Eta distribution of Contamination in Pion Plus",200,-2,2);
    fHistEtaContaminationPionMinus=new TH1F("fHistEtaContaminationPionMinus","Eta distribution of Contamination in Pion Minus",200,-2,2);
    fHistPtContaminationPion=new TH1F("fHistPtContaminationPion","Pt distribution of Contamination in  Pion",200,0,2);
    fHistPtContaminationPionPlus=new TH1F("fHistPtContaminationPionPlus","Pt distribution of Contamination in  Pion Plus",200,0,2);
    fHistPtContaminationPionMinus=new TH1F("fHistPtContaminationPionMinus","Pt distribution of Contamination in  Pion Minus",200,0,2);

    fHistPhiContaminationKaon=new TH1F("fHistPhiContaminationKaon","Phi distribution of Contamination in Kaon",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationKaonPlus=new TH1F("fHistPhiContaminationKaonPlus","Phi distribution of Contamination in Kaon Plus",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationKaonMinus=new TH1F("fHistPhiContaminationKaonMinus","Phi distribution of Contamination in Kaon Minus",200,-1,2.5*TMath::Pi());
    fHistEtaContaminationKaon=new TH1F("fHistEtaContaminationKaon","Eta distribution of Contamination in Kaon",200,-2,2);
    fHistEtaContaminationKaonPlus=new TH1F("fHistEtaContaminationKaonPlus","Eta distribution of Contamination in Kaon Plus",200,-2,2);
    fHistEtaContaminationKaonMinus=new TH1F("fHistEtaContaminationKaonMinus","Eta distribution of Contamination in Kaon Minus",200,-2,2);
    fHistPtContaminationKaon=new TH1F("fHistPtContaminationKaon","Pt distribution of Contamination in  Kaon",200,0,2);
    fHistPtContaminationKaonPlus=new TH1F("fHistPtContaminationKaonPlus","Pt distribution of Contamination in  Kaon Plus",200,0,2);
    fHistPtContaminationKaonMinus=new TH1F("fHistPtContaminationKaonMinus","Pt distribution of Contamination in  Kaon Minus",200,0,2);


    fHistPhiContaminationProton=new TH1F("fHistPhiContaminationProton","Phi distribution of Contamination in Proton",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationProtonPlus=new TH1F("fHistPhiContaminationProtonPlus","Phi distribution of Contamination in Proton Plus",200,-1,2.5*TMath::Pi());
    fHistPhiContaminationProtonMinus=new TH1F("fHistPhiContaminationProtonMinus","Phi distribution of Contamination in Proton Minus",200,-1,2.5*TMath::Pi());
    fHistEtaContaminationProton=new TH1F("fHistEtaContaminationProton","Eta distribution of Contamination in Proton",200,-2,2);
    fHistEtaContaminationProtonPlus=new TH1F("fHistEtaContaminationProtonPlus","Eta distribution of Contamination in Proton Plus",200,-2,2);
    fHistEtaContaminationProtonMinus=new TH1F("fHistEtaContaminationProtonMinus","Eta distribution of Contamination in Proton Minus",200,-2,2);
    fHistPtContaminationProton=new TH1F("fHistPtContaminationProton","Pt distribution of Contamination in  Proton",200,0,2);
    fHistPtContaminationProtonPlus=new TH1F("fHistPtContaminationProtonPlus","Pt distribution of Contamination in  Proton Plus",200,0,2);
    fHistPtContaminationProtonMinus=new TH1F("fHistPtContaminationProtonMinus","Pt distribution of Contamination in  Proton Minus",200,0,2);

    fQAListTruthReco->Add(fHistPhiContaminationPion);
    fQAListTruthReco->Add(fHistEtaContaminationPion);
    fQAListTruthReco->Add(fHistPtContaminationPion);
    fQAListTruthReco->Add(fHistPhiContaminationPionPlus);
    fQAListTruthReco->Add(fHistEtaContaminationPionPlus);
    fQAListTruthReco->Add(fHistPtContaminationPionPlus);
    fQAListTruthReco->Add(fHistPhiContaminationPionMinus);
    fQAListTruthReco->Add(fHistEtaContaminationPionMinus);
    fQAListTruthReco->Add(fHistPtContaminationPionMinus);

    fQAListTruthReco->Add(fHistPhiContaminationKaon);
    fQAListTruthReco->Add(fHistEtaContaminationKaon);
    fQAListTruthReco->Add(fHistPtContaminationKaon);
    fQAListTruthReco->Add(fHistPhiContaminationKaonPlus);
    fQAListTruthReco->Add(fHistEtaContaminationKaonPlus);
    fQAListTruthReco->Add(fHistPtContaminationKaonPlus);
    fQAListTruthReco->Add(fHistPhiContaminationKaonMinus);
    fQAListTruthReco->Add(fHistEtaContaminationKaonMinus);
    fQAListTruthReco->Add(fHistPtContaminationKaonMinus);

    fQAListTruthReco->Add(fHistPhiContaminationProton);
    fQAListTruthReco->Add(fHistEtaContaminationProton);
    fQAListTruthReco->Add(fHistPtContaminationProton);
    fQAListTruthReco->Add(fHistPhiContaminationProtonPlus);
    fQAListTruthReco->Add(fHistEtaContaminationProtonPlus);
    fQAListTruthReco->Add(fHistPtContaminationProtonPlus);
    fQAListTruthReco->Add(fHistPhiContaminationProtonMinus);
    fQAListTruthReco->Add(fHistEtaContaminationProtonMinus);
    fQAListTruthReco->Add(fHistPtContaminationProtonMinus);


if(fDetectorPID_ == kTPC_){
 fHistNsigmaTPCPionBeforePIDCut=new TH2F("HistNsigmaTPCvsPtPionBeforePIDCut","NsigmaTPC vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCKaonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtKaonBeforePIDCut","NsigmaTPC vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCProtonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtProtonBeforePIDCut","NsigmaTPC vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fQAList->Add(fHistNsigmaTPCPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCProtonBeforePIDCut);

 fHistNsigmaTPCPionAfterPIDCut=new TH2F("HistNsigmaTPCvsPtPionAfterPIDCut","NsigmaTPC vs Pt of  Pion AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCKaonAfterPIDCut=new TH2F("HistNsigmaTPCvsPtKaonAfterPIDCut","NsigmaTPC vs Pt of  Kaon AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCProtonAfterPIDCut=new TH2F("HistNsigmaTPCvsPtProtonAfterPIDCut","NsigmaTPC vs Pt of  Proton AfterPIDCut",1000, 0,10,1000, -10, 10);

 fQAList->Add(fHistNsigmaTPCPionAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCKaonAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCProtonAfterPIDCut);
 
}

if(fDetectorPID_ == kTOF_){
 fHistNsigmaTOFPionBeforePIDCut=new TH2F("HistNsigmaTOFvsPtPionBeforePIDCut","NsigmaTOF vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFKaonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtKaonBeforePIDCut","NsigmaTOF vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFProtonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtProtonBeforePIDCut","NsigmaTOF vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);
 
 fQAList->Add(fHistNsigmaTOFPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFProtonBeforePIDCut);

 fHistNsigmaTOFPionAfterPIDCut=new TH2F("HistNsigmaTOFvsPtPionAfterPIDCut","NsigmaTOF vs Pt of  Pion AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFKaonAfterPIDCut=new TH2F("HistNsigmaTOFvsPtKaonAfterPIDCut","NsigmaTOF vs Pt of  Kaon AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFProtonAfterPIDCut=new TH2F("HistNsigmaTOFvsPtProtonAfterPIDCut","NsigmaTOF vs Pt of  Proton AfterPIDCut",1000, 0,10,1000, -10, 10);

 fQAList->Add(fHistNsigmaTOFPionAfterPIDCut);
 fQAList->Add(fHistNsigmaTOFKaonAfterPIDCut);
 fQAList->Add(fHistNsigmaTOFProtonAfterPIDCut);


}

if(fDetectorPID_ == kTPCTOFpid_ || fDetectorPID_ == kTogether_){

 fHistNsigmaTPCTOFPionBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtPionBeforePIDCut","NsigmaTPCTOF vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFKaonBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtKaonBeforePIDCut","NsigmaTPCTOF vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFProtonBeforePIDCut=new TH2F("HistNsigmaTPCTOFvsPtProtonBeforePIDCut","NsigmaTPCTOF vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fHistNsigmaTPCTOFPionAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtPionAfterPIDCut","NsigmaTPCTOF vs Pt of  Pion AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFKaonAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtKaonAfterPIDCut","NsigmaTPCTOF vs Pt of  Kaon AfterPIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCTOFProtonAfterPIDCut=new TH2F("HistNsigmaTPCTOFvsPtProtonAfterPIDCut","NsigmaTPCTOF vs Pt of  Proton AfterPIDCut",1000, 0,10,1000, -10, 10);

 fHistNsigmaTPCPionBeforePIDCut=new TH2F("HistNsigmaTPCvsPtPionBeforePIDCut","NsigmaTPC vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCKaonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtKaonBeforePIDCut","NsigmaTPC vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTPCProtonBeforePIDCut=new TH2F("HistNsigmaTPCvsPtProtonBeforePIDCut","NsigmaTPC vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fHistNsigmaTOFPionBeforePIDCut=new TH2F("HistNsigmaTOFvsPtPionBeforePIDCut","NsigmaTOF vs Pt of  Pion BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFKaonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtKaonBeforePIDCut","NsigmaTOF vs Pt of  Kaon BeforePIDCut",1000, 0,10,1000, -10, 10);
 fHistNsigmaTOFProtonBeforePIDCut=new TH2F("HistNsigmaTOFvsPtProtonBeforePIDCut","NsigmaTOF vs Pt of  Proton BeforePIDCut",1000, 0,10,1000, -10, 10);

 fQAList->Add(fHistNsigmaTPCTOFPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCTOFKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCTOFProtonBeforePIDCut);

 fQAList->Add(fHistNsigmaTPCPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTPCProtonBeforePIDCut);

 fQAList->Add(fHistNsigmaTOFPionBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFKaonBeforePIDCut);
 fQAList->Add(fHistNsigmaTOFProtonBeforePIDCut);

 fQAList->Add(fHistNsigmaTPCTOFPionAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCTOFKaonAfterPIDCut);
 fQAList->Add(fHistNsigmaTPCTOFProtonAfterPIDCut);


}


  //Contamination and efficiency 
  fHistTruthPionPlus = new TH3F("fHistTruthPionPlus","TruthPionPlus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthPionPlus);

  fHistTruthPionMinus = new TH3F("fHistTruthPionMinus","TruthPionMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthPionMinus);

  fHistTruthKaonPlus = new TH3F("fHistTruthKaonPlus","TruthKaonplus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthKaonPlus);

  fHistTruthKaonMinus = new TH3F("fHistTruthKaonMinus","TruthKaonMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthKaonMinus);

  fHistTruthProtonPlus = new TH3F("fHistTruthProtonPlus","TruthProtonPlus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthProtonPlus);

  fHistTruthProtonMinus = new TH3F("fHistTruthProtonMinus","TruthProtonMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthProtonMinus);

  fHistTruthPion = new TH3F("fHistTruthPion","TruthPion;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthPion);
  
  fHistTruthKaon = new TH3F("fHistTruthKaon","TruthKaon;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthKaon);

  fHistTruthProton = new TH3F("fHistTruthProton","TruthProton;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistTruthProton);

  fZvertexTruthPion=new TH1F("fZvertexTruthPion","vertex z truth Pion",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthPion);
  fZvertexTruthPionPlus=new TH1F("fZvertexTruthPionPlus","vertex z truth Pion +",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthPionPlus);
  fZvertexTruthPionMinus=new TH1F("fZvertexTruthPionMinus","vertex z truth Pion -",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthPionMinus);

  fZvertexTruthKaon=new TH1F("fZvertexTruthKaon","vertex z truth Kaon",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthKaon);
  fZvertexTruthKaonPlus=new TH1F("fZvertexTruthKaonPlus","vertex z truth Kaon +",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthKaonPlus);
  fZvertexTruthKaonMinus=new TH1F("fZvertexTruthKaonMinus","vertex z truth Kaon -",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthKaonMinus);
 
  fZvertexTruthProton=new TH1F("fZvertexTruthProton","vertex z truth Proton",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthProton);
  fZvertexTruthProtonPlus=new TH1F("fZvertexTruthProtonPlus","vertex z truth Proton +",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthProtonPlus);
  fZvertexTruthProtonMinus=new TH1F("fZvertexTruthProtonMinus","vertex z truth Proton -",200,-10,10);
  fQAListTruthReco->Add(fZvertexTruthProtonMinus);
 
  //Contamination and Purity Histogram
 
   Hist3dPionContamination =new TH3F("Hist3dPionContamination","Pion Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
   Hist3dPionPlusContamination =new TH3F("Hist3dPionPlusContamination","Positive Pion Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
   Hist3dPionMinusContamination =new TH3F("Hist3dPionMinusContamination","Negative Pion Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());

 fOutputList->Add(Hist3dPionContamination);
 fOutputList->Add(Hist3dPionPlusContamination);
 fOutputList->Add(Hist3dPionMinusContamination);
 Hist3dKaonContamination =new TH3F("Hist3dKaonContamination","Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
 Hist3dKaonPlusContamination =new TH3F("Hist3dKaonPlusContamination","Positive Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
 Hist3dKaonMinusContamination =new TH3F("Hist3dKaonMinusContamination","Negative Kaon Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());


 Hist3dProtonContamination =new TH3F("Hist3dProtonContamination","Proton Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
 Hist3dProtonPlusContamination =new TH3F("Hist3dProtonPlusContamination","Positive Proton Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
 Hist3dProtonMinusContamination =new TH3F("Hist3dProtonMinusContamination","Negative Proton Contamination;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());


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

  fHistMCRecoPionPlus = new TH3F("fHistMCRecoPionPlus","MCRecoPionPlus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoPionPlus);

  fHistMCRecoPionMinus = new TH3F("fHistMCRecoPionMinus","MCRecoPionMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoPionMinus);

  fHistMCRecoKaonPlus = new TH3F("fHistMCRecoKaonPlus","MCRecoKaonplus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoKaonPlus);

  fHistMCRecoKaonMinus = new TH3F("fHistMCRecoKaonMinus","MCRecoKaonMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoKaonMinus);

  fHistMCRecoProtonPlus = new TH3F("fHistMCRecoProtonPlus","MCRecoProtonPlus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoProtonPlus);

  fHistMCRecoProtonMinus = new TH3F("fHistMCRecoProtonMinus","MCRecoProtonMinus;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoProtonMinus);

  fHistMCRecoPion = new TH3F("fHistMCRecoPion","MCRecoPion;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoPion);
  
  fHistMCRecoKaon = new TH3F("fHistMCRecoKaon","MCRecoKaon;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
  fOutputList->Add(fHistMCRecoKaon);

  fHistMCRecoProton = new TH3F("fHistMCRecoProton","MCRecoProton;#eta;p_{T} (GeV/c);#varphi",200,-2,2,200,0.1,2,200,-1,2.5*TMath::Pi());
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

/*fHistEtaPtPhiVertxezTruthPion=new TH4D("fHistEtaPtPhiVertxezTruthPion","EtaPtPhiVzMCTruthPion;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthPionPlus=new TH4D("fHistEtaPtPhiVertxezTruthPionPlus","EtaPtPhiVzMCTruthPionPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthPionMinus=new TH4D("fHistEtaPtPhiVertxezTruthPionMinus","EtaPtPhiVzMCTruthPionMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezTruthKaon=new TH4D("fHistEtaPtPhiVertxezTruthKaon","EtaPtPhiVzMCTruthKaon;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthKaonPlus=new TH4D("fHistEtaPtPhiVertxezTruthKaonPlus","EtaPtPhiVzMCTruthKaonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthKaonMinus=new TH4D("fHistEtaPtPhiVertxezTruthKaonMinus","EtaPtPhiVzMCTruthKaonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezTruthProton=new TH4D("fHistEtaPtPhiVertxezTruthProton","EtaPtPhiVzMCTruthProton;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthProtonPlus=new TH4D("fHistEtaPtPhiVertxezTruthProtonPlus","EtaPtPhiVzMCTruthProtonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezTruthProtonMinus=new TH4D("fHistEtaPtPhiVertxezTruthProtonMinus","EtaPtPhiVzMCTruthProtonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);



fHistEtaPtPhiVertxezContaminationPion=new TH4D("fHistEtaPtPhiVertxezContaminationPion","EtaPtPhiVzMCContaminationPion;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationPionPlus=new TH4D("fHistEtaPtPhiVertxezContaminationPionPlus","EtaPtPhiVzMCContaminationPionPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationPionMinus=new TH4D("fHistEtaPtPhiVertxezContaminationPionMinus","EtaPtPhiVzMCContaminationPionMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezContaminationKaon=new TH4D("fHistEtaPtPhiVertxezContaminationKaon","EtaPtPhiVzMCContaminationKaon;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationKaonPlus=new TH4D("fHistEtaPtPhiVertxezContaminationKaonPlus","EtaPtPhiVzMCContaminationKaonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationKaonMinus=new TH4D("fHistEtaPtPhiVertxezContaminationKaonMinus","EtaPtPhiVzMCContaminationKaonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezContaminationProton=new TH4D("fHistEtaPtPhiVertxezContaminationProton","EtaPtPhiVzMCContaminationProton;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationProtonPlus=new TH4D("fHistEtaPtPhiVertxezContaminationProtonPlus","EtaPtPhiVzMCContaminationProtonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezContaminationProtonMinus=new TH4D("fHistEtaPtPhiVertxezContaminationProtonMinus","EtaPtPhiVzMCContaminationProtonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);


fHistEtaPtPhiVertxezRecoPion=new TH4D("fHistEtaPtPhiVertxezRecoPion","EtaPtPhiVzMCRecoPion;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoPionPlus=new TH4D("fHistEtaPtPhiVertxezRecoPionPlus","EtaPtPhiVzMCRecoPionPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoPionMinus=new TH4D("fHistEtaPtPhiVertxezRecoPionMinus","EtaPtPhiVzMCRecoPionMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezRecoKaon=new TH4D("fHistEtaPtPhiVertxezRecoKaon","EtaPtPhiVzMCRecoKaon;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoKaonPlus=new TH4D("fHistEtaPtPhiVertxezRecoKaonPlus","EtaPtPhiVzMCRecoKaonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoKaonMinus=new TH4D("fHistEtaPtPhiVertxezRecoKaonMinus","EtaPtPhiVzMCRecoKaonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

fHistEtaPtPhiVertxezRecoProton=new TH4D("fHistEtaPtPhiVertxezRecoProton","EtaPtPhiVzMCRecoProton;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoProtonPlus=new TH4D("fHistEtaPtPhiVertxezRecoProtonPlus","EtaPtPhiVzMCRecoProtonPlus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);
fHistEtaPtPhiVertxezRecoProtonMinus=new TH4D("fHistEtaPtPhiVertxezRecoProtonMinus","EtaPtPhiVzMCRecoProtonMinus;#eta;p_{T} (GeV/c);#varphi;V_{z}",200,-2,2,200,0.1,2,200,0,2*TMath::Pi(),24,-6,6);

 fOutputList->Add(fHistEtaPtPhiVertxezTruthPion);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthPionPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthPionMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthKaon);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthKaonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthKaonMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthProton);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthProtonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezTruthProtonMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationPion);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationPionPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationPionMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationKaon);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationKaonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationKaonMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationProton);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationProtonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezContaminationProtonMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoPion);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoPionPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoPionMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoKaon);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoKaonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoKaonMinus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoProton);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoProtonPlus);
 fOutputList->Add(fHistEtaPtPhiVertxezRecoProtonMinus);*/

// Efficieny in Eta and Phi for Different Vertex-z cut 

    fHistEtaVertexzTruthPion=new TH2F("fHistEtaVertexzTruthPion","Eta-Vertexz TruthPion;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthPionPlus=new TH2F("fHistEtaVertexzTruthPionPlus","Eta-Vertexz TruthPion +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthPionMinus=new TH2F("fHistEtaVertexzTruthPionMinus","Eta-Vertexz TruthPion -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzTruthKaon=new TH2F("fHistEtaVertexzTruthKaon","Eta-Vertexz TruthKaon;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthKaonPlus=new TH2F("fHistEtaVertexzTruthKaonPlus","Eta-Vertexz TruthKaon +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthKaonMinus=new TH2F("fHistEtaVertexzTruthKaonMinus","Eta-Vertexz TruthKaon -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzTruthProton=new TH2F("fHistEtaVertexzTruthProton","Eta-Vertexz TruthProton;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthProtonPlus=new TH2F("fHistEtaVertexzTruthProtonPlus","Eta-Vertexz TruthProton +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzTruthProtonMinus=new TH2F("fHistEtaVertexzTruthProtonMinus","Eta-Vertexz TruthProton -;#eta;Vz",200,-2,2,24,-6,6);
    
    fHistEtaVertexzRecoPion=new TH2F("fHistEtaVertexzRecoPion","Eta-Vertexz RecoPion;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoPionPlus=new TH2F("fHistEtaVertexzRecoPionPlus","Eta-Vertexz RecoPion +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoPionMinus=new TH2F("fHistEtaVertexzRecoPionMinus","Eta-Vertexz RecoPion -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzRecoKaon=new TH2F("fHistEtaVertexzRecoKaon","Eta-Vertexz RecoKaon;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoKaonPlus=new TH2F("fHistEtaVertexzRecoKaonPlus","Eta-Vertexz RecoKaon +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoKaonMinus=new TH2F("fHistEtaVertexzRecoKaonMinus","Eta-Vertexz RecoKaon -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzRecoProton=new TH2F("fHistEtaVertexzRecoProton","Eta-Vertexz RecoProton;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoProtonPlus=new TH2F("fHistEtaVertexzRecoProtonPlus","Eta-Vertexz RecoProton +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzRecoProtonMinus=new TH2F("fHistEtaVertexzRecoProtonMinus","Eta-Vertexz RecoProton -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzContaminationPion=new TH2F("fHistEtaVertexzContaminationPion","Eta-Vertexz ContaminationPion;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationPionPlus=new TH2F("fHistEtaVertexzContaminationPionPlus","Eta-Vertexz ContaminationPion +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationPionMinus=new TH2F("fHistEtaVertexzContaminationPionMinus","Eta-Vertexz ContaminationPion -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzContaminationKaon=new TH2F("fHistEtaVertexzContaminationKaon","Eta-Vertexz ContaminationKaon;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationKaonPlus=new TH2F("fHistEtaVertexzContaminationKaonPlus","Eta-Vertexz ContaminationKaon +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationKaonMinus=new TH2F("fHistEtaVertexzContaminationKaonMinus","Eta-Vertexz ContaminationKaon -;#eta;Vz",200,-2,2,24,-6,6);

    fHistEtaVertexzContaminationProton=new TH2F("fHistEtaVertexzContaminationProton","Eta-Vertexz ContaminationProton;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationProtonPlus=new TH2F("fHistEtaVertexzContaminationProtonPlus","Eta-Vertexz ContaminationProton +;#eta;Vz",200,-2,2,24,-6,6);
    fHistEtaVertexzContaminationProtonMinus=new TH2F("fHistEtaVertexzContaminationProtonMinus","Eta-Vertexz ContaminationProton -;#eta;Vz",200,-2,2,24,-6,6);


    fHistPhiVertexzTruthPion=new TH2F("fHistPhiVertexzTruthPion","Phi-Vertexz TruthPion;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthPionPlus=new TH2F("fHistPhiVertexzTruthPionPlus","Phi-Vertexz TruthPion +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthPionMinus=new TH2F("fHistPhiVertexzTruthPionMinus","Phi-Vertexz TruthPion -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzTruthKaon=new TH2F("fHistPhiVertexzTruthKaon","Phi-Vertexz TruthKaon;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthKaonPlus=new TH2F("fHistPhiVertexzTruthKaonPlus","Phi-Vertexz TruthKaon +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthKaonMinus=new TH2F("fHistPhiVertexzTruthKaonMinus","Phi-Vertexz TruthKaon -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzTruthProton=new TH2F("fHistPhiVertexzTruthProton","Phi-Vertexz TruthProton;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthProtonPlus=new TH2F("fHistPhiVertexzTruthProtonPlus","Phi-Vertexz TruthProton +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzTruthProtonMinus=new TH2F("fHistPhiVertexzTruthProtonMinus","Phi-Vertexz TruthProton -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    
    fHistPhiVertexzRecoPion=new TH2F("fHistPhiVertexzRecoPion","Phi-Vertexz RecoPion;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoPionPlus=new TH2F("fHistPhiVertexzRecoPionPlus","Phi-Vertexz RecoPion +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoPionMinus=new TH2F("fHistPhiVertexzRecoPionMinus","Phi-Vertexz RecoPion -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzRecoKaon=new TH2F("fHistPhiVertexzRecoKaon","Phi-Vertexz RecoKaon;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoKaonPlus=new TH2F("fHistPhiVertexzRecoKaonPlus","Phi-Vertexz RecoKaon +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoKaonMinus=new TH2F("fHistPhiVertexzRecoKaonMinus","Phi-Vertexz RecoKaon -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzRecoProton=new TH2F("fHistPhiVertexzRecoProton","Phi-Vertexz RecoProton;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoProtonPlus=new TH2F("fHistPhiVertexzRecoProtonPlus","Phi-Vertexz RecoProton +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzRecoProtonMinus=new TH2F("fHistPhiVertexzRecoProtonMinus","Phi-Vertexz RecoProton -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzContaminationPion=new TH2F("fHistPhiVertexzContaminationPion","Phi-Vertexz ContaminationPion;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationPionPlus=new TH2F("fHistPhiVertexzContaminationPionPlus","Phi-Vertexz ContaminationPion +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationPionMinus=new TH2F("fHistPhiVertexzContaminationPionMinus","Phi-Vertexz ContaminationPion -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzContaminationKaon=new TH2F("fHistPhiVertexzContaminationKaon","Phi-Vertexz ContaminationKaon;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationKaonPlus=new TH2F("fHistPhiVertexzContaminationKaonPlus","Phi-Vertexz ContaminationKaon +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationKaonMinus=new TH2F("fHistPhiVertexzContaminationKaonMinus","Phi-Vertexz ContaminationKaon -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);

    fHistPhiVertexzContaminationProton=new TH2F("fHistPhiVertexzContaminationProton","Phi-Vertexz ContaminationProton;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationProtonPlus=new TH2F("fHistPhiVertexzContaminationProtonPlus","Phi-Vertexz ContaminationProton +;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);
    fHistPhiVertexzContaminationProtonMinus=new TH2F("fHistPhiVertexzContaminationProtonMinus","Phi-Vertexz ContaminationProton -;#phi;Vz",200,-1,2.5*TMath::Pi(),24,-6,6);


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

  if(fDetectorPID_ ==kTPCTOFpid_ || fDetectorPID_ == kTogether_){
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
}

if(fDetectorPID_ == kTPC_) {

  fHistdEdxTPCPionAfterPIDCut = new TH2F("fHistdEdxTPCPionAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);
  fHistdEdxTPCKaonAfterPIDCut = new TH2F("fHistdEdxTPCKaonAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);
  fHistdEdxTPCProtonAfterPIDCut = new TH2F("fHistdEdxTPCProtonAfterPIDCut", ";p_{T} (GeV/c);dE/dx (au.)",1000,-fMaxPt,fMaxPt,1000, 0., 1000.);

  fOutputList->Add(fHistdEdxTPCPionAfterPIDCut);
  fOutputList->Add(fHistdEdxTPCKaonAfterPIDCut);
  fOutputList->Add(fHistdEdxTPCProtonAfterPIDCut);
}

if(fDetectorPID_ == kTOF_){

  fHistBetaTOFPionAfterPIDCut = new TH2F("fHistBetaTOFPionAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);
  fHistBetaTOFKaonAfterPIDCut = new TH2F("fHistBetaTOFKaonAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);
  fHistBetaTOFProtonAfterPIDCut = new TH2F("fHistBetaTOFProtonAfterPIDCut", ";p_{T} (GeV/c);v/c",1000, -fMaxPt, fMaxPt, 1000, 0, 1.2);

  fOutputList->Add(fHistBetaTOFPionAfterPIDCut);
  fOutputList->Add(fHistBetaTOFKaonAfterPIDCut);
  fOutputList->Add(fHistBetaTOFProtonAfterPIDCut);
}

 
  PostData(1, fQAList);
  PostData(2, fOutputList);
  PostData(3, fQAListTruthReco);
}

//________________________________________________________________________
void AliAnalysisTaskEffContPIDBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event
 Double_t lMultiplicityVar     = -999.;
 
  fAOD = dynamic_cast<AliVEvent*>(InputEvent());
//  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
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

 if((lMultiplicityVar = IsEventAccepted(fAOD)) < 0){
    return;
  }


/*

Int_t nMCLabelCounter         = 0;
const Int_t maxMCLabelCounter = 20000;

  Double_t etaMC[maxMCLabelCounter];
  Double_t ptMC[maxMCLabelCounter];
  Double_t phiMC[maxMCLabelCounter];
  Short_t chargeMC[maxMCLabelCounter];
  Int_t pidMC[maxMCLabelCounter];
*/
    

//  Printf("Centrality selection: %lf - %lf",fCentralityPercentileMin,fCentralityPercentileMax);

    Double_t MassPion   = 0.139570; // GeV/c2
    Double_t MassKaon   = 0.493677; // GeV/c2
    Double_t MassProton = 0.938272; // GeV/c2

	     
//-------------------------------------------------MCAOD----------------------------------------------------------------------------------------------------------
                       Short_t vCharge;
                       Double_t vEta;
                       Double_t vPionY;
                       Double_t vKaonY;
                       Double_t vProtonY;
                       Double_t vPhi;
                       Double_t vPt;
                       Double_t vYCut;
                         
                       Int_t nMCParticles = mcEvent->GetNumberOfTracks();
//                       TArrayI labelMCArray(nMCParticles);

	      for(Int_t iParticle = 0; iParticle < mcEvent->GetNumberOfTracks(); iParticle++) {
               
              AliAODMCParticle* currentAODMCParticle = (AliAODMCParticle*) mcEvent->GetTrack(iParticle);

               if (!currentAODMCParticle) {
               AliError(Form("ERROR: Could not receive track %d (mc loop)", iParticle));
               continue;
               }


 
               if (currentAODMCParticle->Charge() == 0) continue;
        if ((currentAODMCParticle->Pt() < fMinPt) || (currentAODMCParticle->Pt() > fMaxPt)) continue;


       if(!fRapidityInsteadOfEta){
 if (currentAODMCParticle->Eta() < fMinEta || currentAODMCParticle->Eta() > fMaxEta) continue;
}




              //  if (!IsMCParticleCut(currentAODMCParticle)) continue;

                if(currentAODMCParticle->IsSecondaryFromWeakDecay()) continue;

                if(!currentAODMCParticle->IsPhysicalPrimary()) continue;




            Int_t pdgCode=((AliAODMCParticle*)currentAODMCParticle)->GetPdgCode();
            
             if (TMath::Abs(pdgCode)==11) continue;
	
            Short_t gAODmcCharge = currentAODMCParticle->Charge();

    
      vCharge = currentAODMCParticle->Charge();
      vEta    = currentAODMCParticle->Eta();
      vPhi    = currentAODMCParticle->Phi();// * TMath::RadToDeg();
      vPt     = currentAODMCParticle->Pt();
//      vYCut   = currentAODMCParticle->Y();


// Eta  and Phi Distribution of All charged Particles

      fHistEtaMCAll->Fill(vEta);
      fHistPhiMCAll->Fill(vPhi);


//cout<<"MCAOD ---- Eta : "<<vEta<<'\t'<<"Phi :"<<vPhi<<'\t'<<"Pt:"<<vPt<<endl;


     //vPionY = log( ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(MassPion*MassPion + vPt*vPt) ); // convert eta to y
     vPionY = 0.5*log( ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) - vPt*sinh(vEta) )); // convert eta to y
     vKaonY = 0.5*log( ( sqrt(MassKaon*MassKaon + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / ( sqrt(MassKaon*MassKaon + vPt*vPt*cosh(vEta)*cosh(vEta)) - vPt*sinh(vEta) )); // convert eta to y
     vProtonY = 0.5*log( ( sqrt(MassProton*MassProton + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / ( sqrt(MassProton*MassProton + vPt*vPt*cosh(vEta)*cosh(vEta)) - vPt*sinh(vEta) ) ); // convert eta to y


     if(fRapidityInsteadOfEta){


           if(TMath::Abs(pdgCode) == 211) vYCut=vPionY;
           else if(TMath::Abs(pdgCode) == 321) vYCut=vKaonY;
           else if(TMath::Abs(pdgCode) == 2212) vYCut=vProtonY;
          else {
           vYCut=999.0;//continue;
          }
 
 }



if(fRapidityInsteadOfEta){

if(fParticleType_ == kPion_){

if( vPionY < fMinEta || vPionY > fMaxEta)  continue;
}

else  if(fParticleType_==kKaon_){
if( vKaonY < fMinEta || vKaonY > fMaxEta)  continue;
}

else  if(fParticleType_==kProton_){
if( vProtonY < fMinEta || vProtonY > fMaxEta)  continue;
}

}

//cout<<"eta fro efiiciein for MCAOd"<<vEta<<endl;
           if(TMath::Abs(pdgCode) == 211){
           fHistPhiMCPion->Fill(vPhi);
           fHistPtMCPion->Fill(vPt);

           fZvertexTruthPion->Fill(fVertexZ);
           fHistPhiVertexzTruthPion->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistTruthPion->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthPion->Fill(vYCut,fVertexZ);
           fHistEtaMCPion->Fill(vYCut);
           }
           else {
           fHistEtaMCPion->Fill(vEta);
           fHistTruthPion->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPion->Fill(vEta,fVertexZ);
           }

           if(gAODmcCharge >0){
           fHistPhiMCPionTruthPlus->Fill(vPhi);
           fHistPtMCPionTruthPlus->Fill(vPt);
           fZvertexTruthPionPlus->Fill(fVertexZ);
           fHistPhiVertexzTruthPionPlus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistEtaMCPionTruthPlus->Fill(vYCut);
           fHistTruthPionPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthPionPlus->Fill(vYCut,fVertexZ);
           }
           else {
           fHistEtaMCPionTruthPlus->Fill(vEta);
           fHistTruthPionPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPionPlus->Fill(vEta,fVertexZ);
           }
            }

           if(gAODmcCharge < 0){
           fHistPhiMCPionTruthMinus->Fill(vPhi);
           fHistPtMCPionTruthMinus->Fill(vPt);
           fZvertexTruthPionMinus->Fill(fVertexZ);
           fHistPhiVertexzTruthPionMinus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta) {
          fHistEtaMCPionTruthMinus->Fill(vYCut);
          fHistTruthPionMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthPionMinus->Fill(vYCut,fVertexZ);
          }
           else {
           fHistEtaMCPionTruthMinus->Fill(vEta);
           fHistTruthPionMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthPionMinus->Fill(vEta,fVertexZ);
           }
           }
        
         } 

        else if(TMath::Abs(pdgCode) == 321){
           fHistPhiMCKaon->Fill(vPhi);
           fHistPtMCKaon->Fill(vPt);

           fZvertexTruthKaon->Fill(fVertexZ);
           fHistPhiVertexzTruthKaon->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistTruthKaon->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthKaon->Fill(vYCut,fVertexZ);
           fHistEtaMCKaon->Fill(vYCut);
           }
           else {
           fHistEtaMCKaon->Fill(vEta);
           fHistTruthKaon->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaon->Fill(vEta,fVertexZ);
           }

           if(gAODmcCharge >0){
           fHistPhiMCKaonTruthPlus->Fill(vPhi);
           fHistPtMCKaonTruthPlus->Fill(vPt);
           fZvertexTruthKaonPlus->Fill(fVertexZ);
           fHistPhiVertexzTruthKaonPlus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistEtaMCKaonTruthPlus->Fill(vYCut);
           fHistTruthKaonPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthKaonPlus->Fill(vYCut,fVertexZ);
           }
           else {
           fHistEtaMCKaonTruthPlus->Fill(vEta);
           fHistTruthKaonPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaonPlus->Fill(vEta,fVertexZ);
           }
            }

           if(gAODmcCharge < 0){
           fHistPhiMCKaonTruthMinus->Fill(vPhi);
           fHistPtMCKaonTruthMinus->Fill(vPt);
           fZvertexTruthKaonMinus->Fill(fVertexZ);
           fHistPhiVertexzTruthKaonMinus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta) {
          fHistEtaMCKaonTruthMinus->Fill(vYCut);
          fHistTruthKaonMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthKaonMinus->Fill(vYCut,fVertexZ);
          }
           else {
           fHistEtaMCKaonTruthMinus->Fill(vEta);
           fHistTruthKaonMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthKaonMinus->Fill(vEta,fVertexZ);
           }
           }
        
         } 




        else if(TMath::Abs(pdgCode) == 2212){
           fHistPhiMCProton->Fill(vPhi);
           fHistPtMCProton->Fill(vPt);

           fZvertexTruthProton->Fill(fVertexZ);
           fHistPhiVertexzTruthProton->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistTruthProton->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthProton->Fill(vYCut,fVertexZ);
           fHistEtaMCProton->Fill(vYCut);
           }
           else {
           fHistEtaMCProton->Fill(vEta);
           fHistTruthProton->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProton->Fill(vEta,fVertexZ);
           }

           if(gAODmcCharge >0){
           fHistPhiMCProtonTruthPlus->Fill(vPhi);
           fHistPtMCProtonTruthPlus->Fill(vPt);
           fZvertexTruthProtonPlus->Fill(fVertexZ);
           fHistPhiVertexzTruthProtonPlus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta){
           fHistEtaMCProtonTruthPlus->Fill(vYCut);
           fHistTruthProtonPlus->Fill(vYCut,vPt,vPhi);
           fHistEtaVertexzTruthProtonPlus->Fill(vYCut,fVertexZ);
           }
           else {
           fHistEtaMCProtonTruthPlus->Fill(vEta);
           fHistTruthProtonPlus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProtonPlus->Fill(vEta,fVertexZ);
           }
            }

           if(gAODmcCharge < 0){
           fHistPhiMCProtonTruthMinus->Fill(vPhi);
           fHistPtMCProtonTruthMinus->Fill(vPt);
           fZvertexTruthProtonMinus->Fill(fVertexZ);
           fHistPhiVertexzTruthProtonMinus->Fill(vPhi,fVertexZ);
           if(fRapidityInsteadOfEta) {
          fHistEtaMCProtonTruthMinus->Fill(vYCut);
          fHistTruthProtonMinus->Fill(vYCut,vPt,vPhi);
          fHistEtaVertexzTruthProtonMinus->Fill(vYCut,fVertexZ);
          }
           else {
           fHistEtaMCProtonTruthMinus->Fill(vEta);
           fHistTruthProtonMinus->Fill(vEta,vPt,vPhi);
           fHistEtaVertexzTruthProtonMinus->Fill(vEta,fVertexZ);
           }
           }
        
         } 




/*          labelMCArray.AddAt(iParticle,nMCLabelCounter);
                  if(nMCLabelCounter >= maxMCLabelCounter){
                    AliWarning(Form("MC Label Counter > Limit (%d) --> stop loop here",maxMCLabelCounter));
                    break;
                  }

 //cout<<"eta fro efiiciencr for MCAOD"<<vEta<<endl;

                if(!fRapidityInsteadOfEta) etaMC[nMCLabelCounter]=vEta;	
                else etaMC[nMCLabelCounter]=vYCut;
                phiMC[nMCLabelCounter]=vPhi;	
                ptMC[nMCLabelCounter]=vPt;
                chargeMC[nMCLabelCounter]=vCharge;
                pidMC[nMCLabelCounter]=pdgCode;	
                nMCLabelCounter += 1;*/

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
//	      TArrayI labelArray(nGoodTracks);
//              Int_t labelCounter = 0;

	      for(Int_t iTracks = 0; iTracks < nGoodTracks; iTracks++) {
		AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));    
		if(!trackAOD) continue;
		

		//track cuts
		if (!trackAOD->TestFilterBit(fAODTrackCutBit)) 
		  continue;

                 Int_t label = TMath::Abs(trackAOD->GetLabel());
//                 Int_t label = trackAOD->GetLabel();
/*                if(IsLabelUsed(labelArray,label)) continue;
                  labelArray.AddAt(label,labelCounter);
                  labelCounter += 1;
           
                  Int_t mcGoods = nMCLabelCounter;
                  for (Int_t kMC = 0; kMC < mcGoods; kMC++) {
                  Int_t mcLabel = labelMCArray.At(kMC);

                  if (mcLabel != TMath::Abs(label)) continue;
                  if(mcLabel != label) continue;
                  if(label > trackAOD->GetLabel()) continue;   */

                  
      vChargeReco = trackAOD->Charge();
      vEtaReco    = trackAOD->Eta();
      vPhiReco    = trackAOD->Phi();// * TMath::RadToDeg();
      vPtReco     = trackAOD->Pt();


      Float_t dcaXY = 0.;
      Float_t DCAZ  = 0.;   // this is the DCA from global track (not exactly what is cut on)

      dcaXY = trackAOD->DCA();      // this is the DCA from global track (not exactly what is cut on)
      DCAZ  = trackAOD->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)

      if((trackAOD->Pt() > fMaxPt)||(trackAOD->Pt() <  fMinPt))
      continue;

      if(!fRapidityInsteadOfEta){
      if(trackAOD->Eta() <fMinEta || trackAOD->Eta() > fMaxEta) continue;
       }



        if( fDCAxyCut != -1 && fDCAzCut != -1){

      Double_t posTrack[3];
      Double_t vertexPos[3];



      const AliVVertex *vertex = fAOD->GetPrimaryVertex();
        vertex->GetXYZ(vertexPos);
        trackAOD->GetXYZ(posTrack);

        Float_t  DCAX = posTrack[0] - vertexPos[0];
        Float_t  DCAY = posTrack[1] - vertexPos[1];
        DCAZ = posTrack[2] - vertexPos[2];


        dcaXY  = TMath::Sqrt(DCAX*DCAX + DCAY*DCAY);


         if (DCAZ     <  -fDCAzCut || DCAZ   > fDCAzCut) continue;
        if( dcaXY    > fDCAxyCut ) continue;

}

       if( fTPCchi2Cut != -1 && trackAOD->Chi2perNDF() > fTPCchi2Cut){
        continue;
      }
      if( fNClustersTPCCut != -1 && trackAOD->GetTPCNcls() < fNClustersTPCCut){
        continue;
      }


//        if(pidMC[kMC]==11) continue;
		  
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
		  
//               Int_t pdgCodeReco = pidMC[kMC];

AliAODMCParticle* recoMC = static_cast<AliAODMCParticle*>(fArrayMC->At(TMath::Abs(trackAOD->GetLabel())));
//AliAODMCParticle* recoMC = static_cast<AliAODMCParticle*>(fArrayMC->At(trackAOD->GetLabel()));
 if(!recoMC) continue;
 if (((AliAODMCParticle*) recoMC)->IsSecondaryFromWeakDecay()) continue;

if (!recoMC->IsPhysicalPrimary()) continue;
Int_t pdgCodeReco = ((AliAODMCParticle*)recoMC)->GetPdgCode();

if(pdgCodeReco == 11) continue;


//cout<<"MCAODReco ---- Eta : "<<vEtaReco<<'\t'<<"Phi :"<<vPhiReco<<'\t'<<"Pt:"<<vPtReco<<endl;


     vPionYReco = 0.5*log( ( sqrt(MassPion*MassPion + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) / ( sqrt(MassPion*MassPion + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) - vPtReco*sinh(vEtaReco) ) ); // convert eta to y
     vKaonYReco = 0.5*log( ( sqrt(MassKaon*MassKaon + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) /  ( sqrt(MassKaon*MassKaon + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) - vPtReco*sinh(vEtaReco) ) ); // convert eta to y
     vProtonYReco = 0.5*log( ( sqrt(MassProton*MassProton + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) + vPtReco*sinh(vEtaReco) ) / ( sqrt(MassProton*MassProton + vPtReco*vPtReco*cosh(vEtaReco)*cosh(vEtaReco)) - vPtReco*sinh(vEtaReco) ) ); // convert eta to y


 /*      if(fRapidityInsteadOfEta){
       if( (vPionYReco < fMinEta || vPionYReco >fMaxEta) || (vKaonYReco < fMinEta || vKaonYReco >fMaxEta) || (vProtonYReco < fMinEta || vProtonYReco >fMaxEta) ) continue;
}*/
/*
 if(fRapidityInsteadOfEta){


           if(TMath::Abs(pdgCodeReco) == 211) vYReco=vPionYReco;
           else if(TMath::Abs(pdgCodeReco) == 321) vYReco=vKaonYReco;
           else if(TMath::Abs(pdgCodeReco) == 2212) vYReco=vProtonYReco;
          else {
           vYReco=999.0;//continue;
          }


      if(vYReco < fMinEta || vYReco > fMaxEta) continue;
 }
*/


if(fRapidityInsteadOfEta){

if(fParticleType_ == kPion_){

if( vPionYReco < fMinEta || vPionYReco > fMaxEta)  continue;
}

else  if(fParticleType_==kKaon_){
if( vKaonYReco < fMinEta || vKaonYReco > fMaxEta)  continue;
}

else  if(fParticleType_==kProton_){
if( vProtonYReco< fMinEta || vProtonYReco > fMaxEta)  continue;
}

}


//cout<<"MCAODreco Eta "<<vEtaReco<<'\t'<<"Pt"<<vPtReco<<'\t'<<"Phi"<<vPhiReco<<endl;



//----------------------------------------------------------------------Nsigma Method for PID--------------------------------------------------------------------------
Double_t dEdx=-1.0;
Double_t beta=-1.0;

Double_t nsigmaSpecies[5]={999.0};

Double_t nsigmaTPC[5] = {999.0};

Double_t nsigmaTOF[5]={999.0};
Double_t nsigmaTPCTOF[5]={999.0};

for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaTPCTOF[iSpecies]  =GetNsigmas(fPIDResponse ,trackAOD,iSpecies);
}

for(int iSpecies=0;iSpecies<5;iSpecies++){

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
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];
//if(nsigmaSpecies[iSpecies] >=999.0) continue;
}

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
else{ continue;}
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
if(trackAOD->Pt()>=fPtTPCMax && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];
//if(nsigmaSpecies[iSpecies] >=999.0) continue;
}
fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
else{ continue;}
}
else if(IsTPCSignal) {
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fPtTPCMax){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTPC[iSpecies]);
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}
fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);
}
else{ continue;}
}
else {
continue;
}
fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);
}

}
else {continue;}

}

else {

if(fDetectorPID_== kTPCTOFpid_){
if(IsTOFSignal){
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
else{ continue;}
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
if(trackAOD->Pt()>=fPtTPCMax && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=nsigmaTPCTOF[iSpecies];
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}

fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);

}
else{ continue;}
}
else if(IsTPCSignal) {
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fPtTPCMax){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTPC[iSpecies]);
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}

fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);

}
else{ continue;}
}

else {
continue;
}

fHistNsigmaTPCTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[2]);
fHistNsigmaTPCTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[3]);
fHistNsigmaTPCTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPCTOF[4]);

}


} // else loop end

if(fDetectorPID_ == kTPC_){
if(IsTPCSignal) {
if(trackAOD->Pt()>=fMinPt && trackAOD->Pt()<=fPtTPCMax){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTPC[iSpecies]);
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}
fHistNsigmaTPCPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[2]);
fHistNsigmaTPCKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[3]);
fHistNsigmaTPCProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTPC[4]);
}
else{ continue;}
}
else {
continue;
}
}

else if(fDetectorPID_ == kTOF_){
if(IsTOFSignal){
if(trackAOD->Pt()>=fPtTPCMax && trackAOD->Pt()<=fMaxPt){
for(int iSpecies=0;iSpecies<5;iSpecies++){
nsigmaSpecies[iSpecies]=TMath::Abs(nsigmaTOF[iSpecies]);
//if(nsigmaSpecies[iSpecies]>=999.0) continue;
}
fHistNsigmaTOFPionBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[2]);
fHistNsigmaTOFKaonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[3]);
fHistNsigmaTOFProtonBeforePIDCut->Fill(trackAOD->Pt(),nsigmaTOF[4]);
}
else{ continue;}
}
else {
continue;
}
}




//Int_t MostProbableSpecies = TMath::LocMin(6, nsigmaSpecies);
Int_t MostProbableSpecies = MinNsigma(5,nsigmaSpecies);

Int_t MostProbableSpecie1 =MostProbableSpecies+2;

if(TMath::Abs(nsigmaSpecies[MostProbableSpecie1]) > fNSigmaPID) continue ;

//cout<<MostProbableSpecies<<'\t'<<nsigmaSpecies[MostProbableSpecies]<<endl;

//if(MostProbableSpecies == 0 || MostProbableSpecies == 1 || MostProbableSpecies == 5 ) continue;


if(MostProbableSpecies == 0) vYReco=vPionYReco;
else if(MostProbableSpecies == 1) vYReco=vKaonYReco;
else if(MostProbableSpecies == 2) vYReco=vProtonYReco;

if(fRapidityInsteadOfEta) {

      if( vYReco <fMinEta || vYReco > fMaxEta)   continue; //cout<<" in of eta range "<<endl;;
}

//Pion

if(TMath::Abs(pdgCodeReco) == 211 && MostProbableSpecies == 0){ 
fHistPhiMCRecoPion->Fill(vPhiReco);
fHistPtMCRecoPion->Fill(vPtReco);

if(fDetectorPID_ == kTPCTOFpid_ || fDetectorPID_ == kTogether_){
fHistdEdxTPCPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFPionAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[2]);
}

if(fDetectorPID_ == kTPC_ ){
fHistdEdxTPCPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);
fHistNsigmaTPCPionAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[2]);
}

if(fDetectorPID_ == kTOF_ ){
fHistBetaTOFPionAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);
fHistNsigmaTOFPionAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[2]);
}

fZvertexRecoPion->Fill(fVertexZ);
fHistPhiVertexzRecoPion->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiMCRecoPionPlus->Fill(vPhiReco);
fHistPtMCRecoPionPlus->Fill(vPtReco);
fZvertexRecoPionPlus->Fill(fVertexZ);
fHistPhiVertexzRecoPionPlus->Fill(vPhiReco,fVertexZ);
}
if(gCharge<0) {
fHistPhiMCRecoPionMinus->Fill(vPhiReco);
fHistPtMCRecoPionMinus->Fill(vPtReco);
fZvertexRecoPionMinus->Fill(fVertexZ);
fHistPhiVertexzRecoPionMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaMCRecoPion->Fill(vYReco);
fHistMCRecoPion->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPion->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoPionPlus->Fill(vYReco);
fHistMCRecoPionPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionPlus->Fill(vYReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoPionMinus->Fill(vYReco);
fHistMCRecoPionMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionMinus->Fill(vYReco,fVertexZ);
}
}

else{
fHistEtaMCRecoPion->Fill(vEtaReco);
fHistMCRecoPion->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPion->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoPionPlus->Fill(vEtaReco);
fHistMCRecoPionPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionPlus->Fill(vEtaReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoPionMinus->Fill(vEtaReco);
fHistMCRecoPionMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoPionMinus->Fill(vEtaReco,fVertexZ);
}
}

}

//Kaon
if(TMath::Abs(pdgCodeReco) == 321 && MostProbableSpecies == 1){ 
fHistPhiMCRecoKaon->Fill(vPhiReco);
fHistPtMCRecoKaon->Fill(vPtReco);

if(fDetectorPID_ == kTPCTOFpid_ || fDetectorPID_ == kTogether_){
fHistdEdxTPCKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFKaonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[3]);
}

if(fDetectorPID_ == kTPC_ ){
fHistdEdxTPCKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);
fHistNsigmaTPCKaonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[3]);
}

if(fDetectorPID_ == kTOF_ ){
fHistBetaTOFKaonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);
fHistNsigmaTOFKaonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[3]);
}

fZvertexRecoKaon->Fill(fVertexZ);
fHistPhiVertexzRecoKaon->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiMCRecoKaonPlus->Fill(vPhiReco);
fHistPtMCRecoKaonPlus->Fill(vPtReco);
fZvertexRecoKaonPlus->Fill(fVertexZ);
fHistPhiVertexzRecoKaonPlus->Fill(vPhiReco,fVertexZ);
}
if(gCharge<0) {
fHistPhiMCRecoKaonMinus->Fill(vPhiReco);
fHistPtMCRecoKaonMinus->Fill(vPtReco);
fZvertexRecoKaonMinus->Fill(fVertexZ);
fHistPhiVertexzRecoKaonMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaMCRecoKaon->Fill(vYReco);
fHistMCRecoKaon->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaon->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoKaonPlus->Fill(vYReco);
fHistMCRecoKaonPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonPlus->Fill(vYReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoKaonMinus->Fill(vYReco);
fHistMCRecoKaonMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonMinus->Fill(vYReco,fVertexZ);
}
}

else{
fHistEtaMCRecoKaon->Fill(vEtaReco);
fHistMCRecoKaon->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaon->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoKaonPlus->Fill(vEtaReco);
fHistMCRecoKaonPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonPlus->Fill(vEtaReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoKaonMinus->Fill(vEtaReco);
fHistMCRecoKaonMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoKaonMinus->Fill(vEtaReco,fVertexZ);
}
}

}


// Proton
if(TMath::Abs(pdgCodeReco) == 2212 && MostProbableSpecies == 2){ 
fHistPhiMCRecoProton->Fill(vPhiReco);
fHistPtMCRecoProton->Fill(vPtReco);

if(fDetectorPID_ == kTPCTOFpid_ || fDetectorPID_ == kTogether_){
fHistdEdxTPCProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);

fHistBetaTOFProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);

fHistNsigmaTPCTOFProtonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[4]);
}

if(fDetectorPID_ == kTPC_ ){
fHistdEdxTPCProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),dEdx);
fHistNsigmaTPCProtonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[4]);
}

if(fDetectorPID_ == kTOF_ ){
fHistBetaTOFProtonAfterPIDCut->Fill(trackAOD->Pt()*trackAOD->Charge(),beta);
fHistNsigmaTOFProtonAfterPIDCut->Fill(trackAOD->Pt(),nsigmaSpecies[4]);
}

fZvertexRecoProton->Fill(fVertexZ);
fHistPhiVertexzRecoProton->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiMCRecoProtonPlus->Fill(vPhiReco);
fHistPtMCRecoProtonPlus->Fill(vPtReco);
fZvertexRecoProtonPlus->Fill(fVertexZ);
fHistPhiVertexzRecoProtonPlus->Fill(vPhiReco,fVertexZ);
}
if(gCharge<0) {
fHistPhiMCRecoProtonMinus->Fill(vPhiReco);
fHistPtMCRecoProtonMinus->Fill(vPtReco);
fZvertexRecoProtonMinus->Fill(fVertexZ);
fHistPhiVertexzRecoProtonMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaMCRecoProton->Fill(vYReco);
fHistMCRecoProton->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProton->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoProtonPlus->Fill(vYReco);
fHistMCRecoProtonPlus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonPlus->Fill(vYReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoProtonMinus->Fill(vYReco);
fHistMCRecoProtonMinus->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonMinus->Fill(vYReco,fVertexZ);
}
}

else{
fHistEtaMCRecoProton->Fill(vEtaReco);
fHistMCRecoProton->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProton->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaMCRecoProtonPlus->Fill(vEtaReco);
fHistMCRecoProtonPlus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonPlus->Fill(vEtaReco,fVertexZ);
}
if(gCharge<0) {
fHistEtaMCRecoProtonMinus->Fill(vEtaReco);
fHistMCRecoProtonMinus->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzRecoProtonMinus->Fill(vEtaReco,fVertexZ);
}
}

}


// Fill the Histograms for Contamination
if(TMath::Abs(pdgCodeReco) !=211 && MostProbableSpecies == 0) {
fHistPhiContaminationPion->Fill(vPhiReco);
fHistPtContaminationPion->Fill(vPtReco);

fZvertexContaminationPion->Fill(fVertexZ);
fHistPhiVertexzContaminationPion->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiContaminationPionPlus->Fill(vPhiReco);
fHistPtContaminationPionPlus->Fill(vPtReco);
fZvertexContaminationPionPlus->Fill(fVertexZ);
fHistPhiVertexzContaminationPionPlus->Fill(vPhiReco,fVertexZ);
}

if(gCharge<0) {
fHistPhiContaminationPionMinus->Fill(vPhiReco);
fHistPtContaminationPionMinus->Fill(vPtReco);
fZvertexContaminationPionMinus->Fill(fVertexZ);
fHistPhiVertexzContaminationPionMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaContaminationPion->Fill(vYReco);

Hist3dPionContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPion->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationPionPlus->Fill(vYReco);
Hist3dPionPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionPlus->Fill(vYReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationPionMinus->Fill(vYReco);
Hist3dPionMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionMinus->Fill(vYReco,fVertexZ);
}
}

else {
fHistEtaContaminationPion->Fill(vEtaReco);
Hist3dPionContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPion->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationPionPlus->Fill(vEtaReco);
Hist3dPionPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionPlus->Fill(vEtaReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationPionMinus->Fill(vEtaReco);
Hist3dPionMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationPionMinus->Fill(vEtaReco,fVertexZ);
}
}

}

if(TMath::Abs(pdgCodeReco) !=321 && MostProbableSpecies == 1) {
fHistPhiContaminationKaon->Fill(vPhiReco);
fHistPtContaminationKaon->Fill(vPtReco);

fZvertexContaminationKaon->Fill(fVertexZ);
fHistPhiVertexzContaminationKaon->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiContaminationKaonPlus->Fill(vPhiReco);
fHistPtContaminationKaonPlus->Fill(vPtReco);
fZvertexContaminationKaonPlus->Fill(fVertexZ);
fHistPhiVertexzContaminationKaonPlus->Fill(vPhiReco,fVertexZ);
}

if(gCharge<0) {
fHistPhiContaminationKaonMinus->Fill(vPhiReco);
fHistPtContaminationKaonMinus->Fill(vPtReco);
fZvertexContaminationKaonMinus->Fill(fVertexZ);
fHistPhiVertexzContaminationKaonMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaContaminationKaon->Fill(vYReco);

Hist3dKaonContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaon->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationKaonPlus->Fill(vYReco);
Hist3dKaonPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonPlus->Fill(vYReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationKaonMinus->Fill(vYReco);
Hist3dKaonMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonMinus->Fill(vYReco,fVertexZ);
}
}

else {
fHistEtaContaminationKaon->Fill(vEtaReco);
Hist3dKaonContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaon->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationKaonPlus->Fill(vEtaReco);
Hist3dKaonPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonPlus->Fill(vEtaReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationKaonMinus->Fill(vEtaReco);
Hist3dKaonMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationKaonMinus->Fill(vEtaReco,fVertexZ);
}
}

}


if(TMath::Abs(pdgCodeReco) !=2212 && MostProbableSpecies == 2) {
fHistPhiContaminationProton->Fill(vPhiReco);
fHistPtContaminationProton->Fill(vPtReco);

fZvertexContaminationProton->Fill(fVertexZ);
fHistPhiVertexzContaminationProton->Fill(vPhiReco,fVertexZ);

if(gCharge>0) {
fHistPhiContaminationProtonPlus->Fill(vPhiReco);
fHistPtContaminationProtonPlus->Fill(vPtReco);
fZvertexContaminationProtonPlus->Fill(fVertexZ);
fHistPhiVertexzContaminationProtonPlus->Fill(vPhiReco,fVertexZ);
}

if(gCharge<0) {
fHistPhiContaminationProtonMinus->Fill(vPhiReco);
fHistPtContaminationProtonMinus->Fill(vPtReco);
fZvertexContaminationProtonMinus->Fill(fVertexZ);
fHistPhiVertexzContaminationProtonMinus->Fill(vPhiReco,fVertexZ);
}

if(fRapidityInsteadOfEta){
fHistEtaContaminationProton->Fill(vYReco);

Hist3dProtonContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProton->Fill(vYReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationProtonPlus->Fill(vYReco);
Hist3dProtonPlusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonPlus->Fill(vYReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationProtonMinus->Fill(vYReco);
Hist3dProtonMinusContamination->Fill(vYReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonMinus->Fill(vYReco,fVertexZ);
}
}

else {
fHistEtaContaminationProton->Fill(vEtaReco);
Hist3dProtonContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProton->Fill(vEtaReco,fVertexZ);
if(gCharge>0) {
fHistEtaContaminationProtonPlus->Fill(vEtaReco);
Hist3dProtonPlusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonPlus->Fill(vEtaReco,fVertexZ);
}

if(gCharge<0) {
fHistEtaContaminationProtonMinus->Fill(vEtaReco);
Hist3dProtonMinusContamination->Fill(vEtaReco,vPtReco,vPhiReco);
fHistEtaVertexzContaminationProtonMinus->Fill(vEtaReco,fVertexZ);
}
}

}


// Fill the Histograms for Contamination


// PID selection end ================================================================================================================================================

//} // MCgood tracks

		
	       }//AOD track loop 

//              labelMCArray.Reset();
//              labelArray.Reset();
	      
}


Bool_t AliAnalysisTaskEffContPIDBF::IsMCParticleCut(AliAODMCParticle* particle) {


//cout<<"Min Pt ---- "<<fMinPt<<"Max Pt "<<fMaxPt<<"Pt TPC Max "<<fPtTPCMax<<endl;


        if (particle->Charge() == 0) {return kFALSE;}
        if ((particle->Pt() < fMinPt) || (particle->Pt() > fMaxPt)) {return kFALSE;}

       
       if(!fRapidityInsteadOfEta){
 if (particle->Eta() < fMinEta || particle->Eta() > fMaxEta) {return kFALSE;}
}

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
/*
    AliPIDResponse::EDetPidStatus statusTPC = PIDresponse->NumberOfSigmas(AliPIDResponse::kTPC, track, (AliPID::EParticleType) specie, nsigmaTPC);
    AliPIDResponse::EDetPidStatus statusTOF = PIDresponse->NumberOfSigmas(AliPIDResponse::kTOF, track, (AliPID::EParticleType) specie, nsigmaTOF);
    Bool_t tpcIsOk = (statusTPC == AliPIDResponse::kDetPidOk);
    Bool_t tofIsOk = (statusTOF == AliPIDResponse::kDetPidOk);*/


     nsigmaTPC= PIDresponse->NumberOfSigmasTPC((AliVTrack*)track,(AliPID::EParticleType)specie);
     nsigmaTOF= PIDresponse->NumberOfSigmasTOF((AliVTrack*)track,(AliPID::EParticleType)specie);

   Bool_t fHasTPCPID=kFALSE;
   Bool_t fHasTOFPID=kFALSE;

  if(track->GetTPCsignal() >0.0)
   fHasTPCPID=kTRUE;

   if(IsTOFPID(track) && (Beta(track) >0.0 && Beta(track)<=1.0))
   fHasTOFPID=kTRUE;
   
    fNsigmaTPC[specie] = fHasTPCPID? nsigmaTPC: 999.0;
    fNsigmaTOF[specie] = fHasTOFPID? nsigmaTOF: 999.0;
//    if(fHasTOFPID && fNsigmaTOF[specie] == -998) fNsigmaTOF[specie]=999.0;

//cout<<"nsigma TOF :"<<fNsigmaTOF[specie]<<'\t'<<"nsigma TPC: "<<fNsigmaTPC[specie]<<'\t'<<"nsigma TPCTOF:"<<TMath::Hypot(fNsigmaTPC[specie],fNsigmaTOF[specie])<<endl;
    return TMath::Hypot(fNsigmaTPC[specie], fNsigmaTOF[specie]);
}


 Int_t AliAnalysisTaskEffContPIDBF::MinNsigma(Int_t n, const Double_t *b)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  if  (n <= 0) return -1;

Double_t a[3];
a[0]=b[2];
a[1]=b[3];
a[2]=b[4];

if(n>3) n=3;

  Double_t xmin = a[0];
  Int_t loc = 0;
  for  (Int_t i = 1; i < n; i++) {
     if (xmin > a[i])  {
         xmin = a[i];
         loc = i;
     }
  }
  return loc;
}


Double_t AliAnalysisTaskEffContPIDBF::IsEventAccepted(AliVEvent *event){

Float_t gRefMultiplicity = -1.;

fHistEventStats->Fill(1);


   const AliVVertex *vertex = event->GetPrimaryVertex();

      if(vertex) {
        Double32_t fCov[6];
        vertex->GetCovarianceMatrix(fCov);
        if(vertex->GetNContributors() > 0) {
          if(fCov[5] != 0) {
            fHistEventStats->Fill(3); //proper vertex
            if(TMath::Abs(vertex->GetX()) < fVxMax) {
              if(TMath::Abs(vertex->GetY()) < fVyMax) {
                if(TMath::Abs(vertex->GetZ()) < fVzMax) {
               gRefMultiplicity = GetRefMultiOrCentrality(event);

                if(fUseCentrality) {
                    if((gRefMultiplicity > fCentralityPercentileMin) && (gRefMultiplicity < fCentralityPercentileMax)){
                     fHistCentrality->Fill(gRefMultiplicity);
                     fHistEventStats->Fill(2);

                     return gRefMultiplicity;
}
}
              fHistEventStats->Fill(4); //analyzed events
              fHistVz->Fill(fVertexZ);
              
              fVertexZ=vertex->GetZ();


             } // vetex z
           } // vertex Y
         } // vertex x

       } // fCov[5]

     } // GetNContributor

   }  // Vetex 

return -1;


}  // Event Loop end


Double_t AliAnalysisTaskEffContPIDBF::GetRefMultiOrCentrality(AliVEvent *event){
    // Checks the Event cuts
    // Fills Event statistics histograms

  Float_t gCentrality = -1.;

  AliAODHeader *header = (AliAODHeader*) event->GetHeader();
      if(header) gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());


 return  gCentrality;

} // Event loop end


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
