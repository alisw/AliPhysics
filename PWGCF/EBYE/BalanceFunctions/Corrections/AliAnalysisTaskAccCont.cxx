#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"

#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliAODInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODMCParticle.h"

#include <AliPID.h>
#include <AliPIDCombined.h>
#include <AliPIDResponse.h>
#include "AliAnalysisUtils.h"
#include "AliESDpid.h"

#include "AliAnalysisTaskAccCont.h"

#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskAccCont)

//________________________________________________________________________
AliAnalysisTaskAccCont::AliAnalysisTaskAccCont(const char *name)
: AliAnalysisTaskSE(name),
  gAOD(0),
  fListQA(0),
  fListResults(0),
  fHistEventStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fCentralityEstimator("V0M"),
  fCentralityPercentileMin(0.),
  fCentralityPercentileMax(5.),
  fHistMultiplicity(0),
  fHistPtCen(0),
  fHistPhiCen(0),
  fHistEtaCen(0),
  fHistDCAToVertex2D(0),
  fHistPt(0),
  fHistPtbin(0),
  fHistCent(0),
  fHistCentbin(0),
  fHistPhi(0),
  fHistEta(0),
  fHistNClustersTPC(0),
  fHistChi2PerClusterTPC(0),
  fHistDCAToVertexZ(0),
  fHistDCAToVertexXY(0),
  fHistEtaPhiCent(0),
  fHistPtEtaCent(0),
  fHistPtPhiCent(0),
  fHistEtaPhiVertexPlus(0),
  fHistEtaPhiVertexMinus(0),
  fHistYPhiVertexPlus(0),
  fHistYPhiVertexMinus(0),
  fHistDCAXYptchargedminus(0),
  fHistDCAXYptchargedplus(0),
  fHistDCAXYptchargedminus_ext(0),
  fHistDCAXYptchargedplus_ext(0),
  fHistGlobalvsESDBeforePileUpCuts(0),
  fHistGlobalvsESDAfterPileUpCuts(0),
  fHistV0MvsTPCoutBeforePileUpCuts(0), 
  fHistV0MvsTPCoutAfterPileUpCuts(0),
  hNSigmaCutApplied(0),
  hBayesProbab(0),
  fUseOfflineTrigger(kFALSE),
  fPbPb(kFALSE),
  fpPb(kFALSE),
  fCheckPileUp(kFALSE),
  fMCrec(kFALSE),
  fArrayMC(0),
  fExcludeSecondariesInMCrec(kFALSE),
  fExcludeElectronsInMCrec(kFALSE),
  fExcludeInjectedSignals(kFALSE),
  fRejectCheckGenName(kFALSE),
  fGenToBeKept("Hijing"),
  fPileupLHC15oSlope(3.38),
  fPileupLHC15oOffset(15000),
  fUseOutOfBunchPileUpCutsLHC15o(kFALSE),
  fUseOutOfBunchPileUpCutsLHC15oJpsi(kFALSE),
  fUsePID(kFALSE),
  fDCAext(kFALSE),
  fUseRapidity(kFALSE),
  fUsePIDnSigmaComb(kFALSE),
  fVxMax(0.5),
  fVyMax(0.5),
  fVzMax(10.),
  fAODtrackCutBit(128),
  fPIDNSigma(3),
  fMassParticleOfInterest(0.13957),
  fParticleOfInterest(AliPID::kPion),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fPtMin(0),
  fPtMax(20),
  fBayesPIDThr(0.8),
  fPIDMomCut(0.7),
  fUtils(0),
  fPIDResponse(0),
  fPIDCombined(0),
  fHistPdg(0),
  fUsePIDNewTrial(kFALSE),
  fHistdEdxVsPTPCbeforePID(0),
  fHistBetavsPTOFbeforePID(0),
  fHistProbTPCvsPtbeforePID(0),
  fHistProbTPCTOFvsPtbeforePID(0),
  fHistNSigmaTPCvsPtbeforePID(0),
  fHistNSigmaTOFvsPtbeforePID(0),
  fHistBetaVsdEdXbeforePID(0),
  fHistNSigmaTPCTOFvsPtbeforePID(0),
  fHistNSigmaTPCTOFPbefPID(0),
  fHistBetavsPTOFafterPID(0),
  fHistdEdxVsPTPCafterPID(0),
  fHistBetaVsdEdXafterPID(0),
  fHistNSigmaTOFvsPtafterPID(0),
  fHistNSigmaTPCvsPtafterPID(0),
  fHistNSigmaTPCTOFvsPtafterPID(0),
  fHistNSigmaTPCTOFPafterPID(0)
{
    // Constructor
    
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 writes into a TH1 container
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskAccCont::~AliAnalysisTaskAccCont() {
    //
}

//________________________________________________________________________
void AliAnalysisTaskAccCont::UserCreateOutputObjects() {
    // Create histograms
    // Called once
    Int_t phiBin = 100;
    Int_t etaBin = 16;
    Int_t vertex_bin = 9;
    
    Double_t nArrayPhi[phiBin+1];
    for(Int_t iBin = 0; iBin <= phiBin; iBin++)
        nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;
    
    Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    Double_t nArrayVertex[10]={-10, -7, -5, -3, -1, 1, 3, 5, 7, 10}; 
    
    fUtils = new AliAnalysisUtils();
    
    //QA list
    fListQA = new TList();
    fListQA->SetName("listQA");
    fListQA->SetOwner();
    
    //Event stats.
    TString gCutName[6] = {"Total","Offline trigger",
        "Vertex","sel. Centrality",
        "No Pile-Up", "Out-of-bunch Pile-Up"};
    fHistEventStats = new TH2F("fHistEventStats",
                               "Event statistics;;Centrality percentile;N_{events}",
                               6,0.5,6.5,110,-5,105);
    for(Int_t i = 1; i <= 6; i++)
        fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
    fListQA->Add(fHistEventStats);
    
    fHistTrackStats = new TH2F("fHistTrackStats","Event statistics;Centrality (%);FilterBit",110,-5,105,1300,0,1300);
    fListQA->Add(fHistTrackStats);
    
    // Vertex distributions
    fHistVx = new TH2F("fHistVx","Primary vertex distribution - x coordinate;Centrality (%);V_{x} (cm);Entries",10,0,100,100,-0.5,0.5);
    fListQA->Add(fHistVx);
    fHistVy = new TH2F("fHistVy","Primary vertex distribution - y coordinate;Centrality (%);V_{y} (cm);Entries",10,0,100,100,-0.5,0.5);
    fListQA->Add(fHistVy);
    fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;Centrality (%);V_{z} (cm);Entries",10,0,100,100,-20.,20.);
    fListQA->Add(fHistVz);
    
    // QA histograms: multiplicities
    fHistMultiplicity = new TH2F("fHistMultiplicity",";Centrality (%);N_{acc.};Counts",10,0,100,5000,-0.5,4999.5);
    fListQA->Add(fHistMultiplicity);
  
    fHistGlobalvsESDBeforePileUpCuts = new TH2F("fHistGlobalvsESDBeforePileUpCuts","Global vs ESD Tracks; ESD tracks; Global tracks;",1000,0,20000,100,0,20000);
    fHistGlobalvsESDAfterPileUpCuts = new TH2F("fHistGlobalvsESDAfterPileUpCuts","Global vs ESD Tracks; ESD tracks; Global tracks;",1000,0,20000,100,0,20000);
    
    fHistV0MvsTPCoutBeforePileUpCuts = new TH2F("fHistV0MvsTPCoutBeforePileUpCuts","V0M amplitude vs TPCout tracks; TPCout tracks; V0M amplitude;",1000,0,20000,1000,0,40000);
    fHistV0MvsTPCoutAfterPileUpCuts = new TH2F("fHistV0MvsTPCoutAfterPileUpCuts","V0M amplitude vs TPCout tracks; TPCout tracks; V0M amplitude;",1000,0,20000,1000,0,40000);
    
    fListQA->Add(fHistGlobalvsESDBeforePileUpCuts);
    fListQA->Add(fHistGlobalvsESDAfterPileUpCuts);

    fListQA->Add(fHistV0MvsTPCoutBeforePileUpCuts);
    fListQA->Add(fHistV0MvsTPCoutAfterPileUpCuts);
    
    //====================================================//
    //Results TList
    fListResults = new TList();
    fListResults->SetName("listResults");
    fListResults->SetOwner();
    
    //Results
    fHistPt = new TH1F("fHistPt","Pt distribution;p_{T} (GeV/c);Counts",100,0,10);
    fListResults->Add(fHistPt);
    fHistCent = new TH1F("fHistCent","Centrality distribution;Centrality;Counts",101,0,101);
    fListResults->Add(fHistCent);
    fHistPtbin = new TH1F("fHistPtbin","Pt distribution;p_{T} (GeV/c);Counts",500,0,10);
    fListResults->Add(fHistPtbin);
    fHistCentbin = new TH1F("fHistCentbin","Centrality distribution;Centrality;Counts",100,0,100);
    fListResults->Add(fHistCentbin);
    fHistPtCen = new TH2F("fHistPtCen",";Centrality (%);p_{T} (GeV/c);Counts",10,0,100,110,-0.5,10.5);
    fListResults->Add(fHistPtCen);
    fHistPhi = new TH1F("fHistPhi","Phi distribution;Phi;Number Of Entries",100,0,2*(TMath::Pi()));
    fListResults->Add(fHistPhi);
    fHistEta = new TH1F("fHistEta","Eta distribution;Eta;Number Of Entries",100,-1,1);
    fListResults->Add(fHistEta);
    fHistPhiCen = new TH2F("fHistPhiCen","Phi vs Centrality;Centrality;Phi;Counts",10,0,100,72,0,2*TMath::Pi());
    fListResults->Add(fHistPhiCen);
    fHistEtaCen = new TH2F("fHistEtaCen","Eta vs Centrality;Centrality;Eta;Counts",10,0,100,100,-1,1);
    fListResults->Add(fHistEtaCen);
    fHistNClustersTPC = new TH1F("fHistNClustersTPC","Number of clusters in TPC;Number of clusters;Counts",100,0,200);
    fListResults->Add(fHistNClustersTPC);
    fHistChi2PerClusterTPC = new TH1F("fHistChi2PerClusterTPC","Chi2 of the fit;Chi^2;Counts",100,0,50);
    fListResults->Add(fHistChi2PerClusterTPC);
    fHistDCAToVertexXY = new TH1F("fHistMaxDCAToVertexXY","DCA_{xy};DCA_{xy};Counts",100,-20,20);
    fListResults->Add(fHistDCAToVertexXY);
    fHistDCAToVertexZ = new TH1F("fHistDCAToVertexZ","DCA_{z};DCA_{z};Counts",100,-20,20);
    fListResults->Add(fHistDCAToVertexZ);
    fHistDCAToVertex2D = new TH2F("fHistDCAToVertex2D","DCA_2D;DCA_{y};DCA_{z};Counts",100,-20,20,100,-20,20);
    fListResults->Add(fHistDCAToVertex2D);
    fHistEtaPhiCent = new TH3F("fHistEtaPhiCent","#eta & #phi vs Centrality;#eta;#phi;Centrality",100,-1,1,100,0,2*(TMath::Pi()),100,0,100);
    fListResults->Add(fHistEtaPhiCent);
    fHistPtEtaCent = new TH3F("fHistPtEtaCent","p_{T} & #eta vs Centrality;p_{T}(GeV/c);#eta;Centrality",100,0,10,100,-1,1,100,0,100);
    fListResults->Add(fHistPtEtaCent);
    fHistPtPhiCent = new TH3F("fHistPtPhiCent","p_{T} & #phi vs Centrality;#eta;#phi;Centrality",100,0,10,100,0,2*(TMath::Pi()),100,0,100);
    fListResults->Add(fHistPtPhiCent);
    
    fHistEtaPhiVertexPlus = new TH3D("fHistEtaPhiVertexPlus",
                                     "Survived positive primaries;#phi;#eta;V_{z} (cm)",
                                     phiBin, nArrayPhi, etaBin, nArrayEta, vertex_bin, nArrayVertex);
    
    fHistEtaPhiVertexMinus = new TH3D("fHistEtaPhiVertexMinus",
                                      "Survived negative primaries;#phi;#eta;V_{z} (cm)",
                                      phiBin, nArrayPhi, etaBin,nArrayEta, vertex_bin, nArrayVertex);

    fHistYPhiVertexPlus = new TH3D("fHistYPhiVertexPlus",
                                     "Survived positive primaries;#phi;y;V_{z} (cm)",
                                     phiBin, nArrayPhi, etaBin, nArrayEta, vertex_bin, nArrayVertex);

    fHistYPhiVertexMinus = new TH3D("fHistYPhiVertexMinus",
                                      "Survived negative primaries;#phi;y;V_{z} (cm)",
                                      phiBin, nArrayPhi, etaBin, nArrayEta, vertex_bin, nArrayVertex);
    
    fHistDCAXYptchargedminus = new TH3F("fHistDCAxychargedminus","DCA_{xy} vs pt for charged particles (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptchargedplus = new TH3F("fHistDCAxychargedplus","DCA_{xy} vs pt for charged particles (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-0.5,0.5);
    fHistDCAXYptchargedminus_ext = new TH3F("fHistDCAxychargedminusext","DCA_{xy} vs pt for charged particles (negative);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    fHistDCAXYptchargedplus_ext = new TH3F("fHistDCAxychargedplusext","DCA_{xy} vs pt for charged particles (positive);p_{T} [GeV/c];#eta;DCA_{xy}",100,0,10,16,-0.8,0.8,1000,-4,4);
    
    fListResults->Add(fHistEtaPhiVertexPlus);
    fListResults->Add(fHistEtaPhiVertexMinus);
    if(fUseRapidity){
      fListResults->Add(fHistYPhiVertexPlus);
      fListResults->Add(fHistYPhiVertexMinus);
    }
    fListResults->Add(fHistDCAXYptchargedminus);
    fListResults->Add(fHistDCAXYptchargedplus);
    if(fDCAext){
    fListResults->Add(fHistDCAXYptchargedminus_ext);
    fListResults->Add(fHistDCAXYptchargedplus_ext);
    }

    hNSigmaCutApplied = new TH3F("hNSigmaCutApplied", "nsigma vs p vs pT for ID particles; n_{#sigma}; p [GeV/c]; p_{T} [GeV/c]", 50, 0, 4, 100, 0, 10, 100, 0, 10);
    fListResults->Add(hNSigmaCutApplied);

    hBayesProbab = new TH3F("hBayesProbab", "BayesProbab vs p vs pT for ID particles; Probability;p [GeV/c]; p_{T} [GeV/c]", 100, 0.5, 1, 100, 0, 10, 100, 0, 10);
    fListResults->Add(hBayesProbab);
  
    fHistPdg  = new TH1F("fHistPdg","Pdg code distribution;pdg code;Entries",6401,-3200.5,3200.5);
    fListResults->Add(fHistPdg);
    
    if(fUsePIDNewTrial) {
        
        fHistdEdxVsPTPCbeforePID = new TH2D ("dEdxVsPTPCbefore","dEdxVsPTPCbefore", 1000, -10.0, 10.0, 1000, 0, 1000);
        fListQA->Add(fHistdEdxVsPTPCbeforePID);
        
        fHistBetavsPTOFbeforePID = new TH2D ("BetavsPTOFbefore","BetavsPTOFbefore", 1000, -10.0, 10., 1000, 0, 1.2);
        fListQA->Add(fHistBetavsPTOFbeforePID);
        
        fHistProbTPCvsPtbeforePID = new TH2D ("ProbTPCvsPtbefore","ProbTPCvsPtbefore", 1000, -10.0,10.0, 1000, 0, 2.0);
        fListQA->Add(fHistProbTPCvsPtbeforePID);
        
        fHistProbTPCTOFvsPtbeforePID =new TH2D ("ProbTPCTOFvsPtbefore","ProbTPCTOFvsPtbefore", 1000, -50, 50, 1000, 0, 2.0);
        fListQA->Add(fHistProbTPCTOFvsPtbeforePID);
        
        fHistNSigmaTPCvsPtbeforePID = new TH2D ("NSigmaTPCvsPtbefore","NSigmaTPCvsPtbefore", 1000, -10, 10, 1000, -25, 25);
        fListQA->Add(fHistNSigmaTPCvsPtbeforePID);
        
        fHistNSigmaTOFvsPtbeforePID = new TH2D ("NSigmaTOFvsPtbefore","NSigmaTOFvsPtbefore", 1000, -10, 10, 1000, -25, 25);
        fListQA->Add(fHistNSigmaTOFvsPtbeforePID);
        
        fHistBetaVsdEdXbeforePID = new TH2D ("BetaVsdEdXbefore","BetaVsdEdXbefore", 1000, 0., 1000, 1000, 0, 1.2);
        fListQA->Add(fHistBetaVsdEdXbeforePID);
        
        fHistNSigmaTPCTOFvsPtbeforePID = new TH2D ("NSigmaTPCTOFvsPtbefore","NSigmaTPCTOFvsPtbefore", 1000, -10., 10., 1000, -25, 25);
        fListQA->Add(fHistNSigmaTPCTOFvsPtbeforePID);
        
        const Int_t pBins = 36;
        Double_t nArrayP[pBins+1]={0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0};
        //nSigma Array
        const Int_t nSigmaBins = 250;
        Double_t nArrayS[nSigmaBins+1];
        for (Int_t i = 0; i <= nSigmaBins; i++){
            nArrayS[i]=i-125; //i+1
            //Printf("nS: %lf - i: %d", nSigmaArray[i], i);
        }
        
        fHistNSigmaTPCTOFPbefPID = new TH3D ("fHistNSigmaTPCTOFPbefPID","fHistNSigmaTPCTOFPbefPID;#sigma_{TPC};#sigma_{TOF};p_{T} (GeV/c)", nSigmaBins, nArrayS, nSigmaBins, nArrayS, pBins,nArrayP);
        fListQA->Add(fHistNSigmaTPCTOFPbefPID);
        
        fHistdEdxVsPTPCafterPID = new TH2D ("dEdxVsPTPCafter","dEdxVsPTPCafter", 1000, -10, 10, 1000, 0, 1000);
        fListQA->Add(fHistdEdxVsPTPCafterPID);
        
        fHistBetavsPTOFafterPID = new TH2D ("BetavsPTOFafter","BetavsPTOFafter", 1000, -10, 10, 1000, 0, 1.2);
        fListQA->Add(fHistBetavsPTOFafterPID);
        
        fHistNSigmaTPCvsPtafterPID = new TH2D ("NSigmaTPCvsPtafter","NSigmaTPCvsPtafter", 1000, -10, 10, 1000, -25, 25);
        fListQA->Add(fHistNSigmaTPCvsPtafterPID);
        
        fHistNSigmaTOFvsPtafterPID = new TH2D ("NSigmaTOFvsPtafter","NSigmaTOFvsPtafter", 1000, -10, 10, 1000, -25, 25);
        fListQA->Add(fHistNSigmaTOFvsPtafterPID);
        
        fHistBetaVsdEdXafterPID = new TH2D ("BetaVsdEdXafter","BetaVsdEdXafter", 1000, 0., 1000, 1000, 0, 1.2);
        fListQA->Add(fHistBetaVsdEdXafterPID);
        
        fHistNSigmaTPCTOFvsPtafterPID = new TH2D ("NSigmaTPCTOFvsPtafter","NSigmaTPCTOFvsPtafter", 1000, -10., 10., 1000, -25, 25);
        fListQA->Add(fHistNSigmaTPCTOFvsPtafterPID);
        
        fHistNSigmaTPCTOFPafterPID = new TH3D ("fHistNSigmaTPCTOFPafterPID","fHistNSigmaTPCTOFPafterPID;#sigma_{TPC};#sigma_{TOF};p_{T} (GeV/c)", nSigmaBins, nArrayS, nSigmaBins, nArrayS, pBins,nArrayP);
        fListQA->Add(fHistNSigmaTPCTOFPafterPID);
        
    }
     
    // Post output data
    PostData(1, fListQA);
    PostData(2, fListResults);
}

//________________________________________________________________________
void AliAnalysisTaskAccCont::UserExec(Option_t *) {
    // Main loop
    // Called for each event
    //AOD analysis (vertex and track cuts also here!!!!)
    gAOD = dynamic_cast<AliAODEvent*>(InputEvent()); // from TaskSE
    if(!gAOD) {
        Printf("ERROR: gAOD not available");
        return;
    }
    
    AliAODHeader *header = dynamic_cast<AliAODHeader *>(gAOD->GetHeader());
    if(!header) {
        Printf("ERROR: AOD header not available");
        return;
    }

    if (fMCrec){
        
        fArrayMC = dynamic_cast<TClonesArray*>(gAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        
        if (!fArrayMC) {
            AliError("No array of MC particles found !!!");
        }
    
    }

    Int_t nAcceptedTracks = 0;
    Float_t gCentrality = -1;
    
    //Centrality
    // if(header)
    //gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());
    
    AliMultSelection *multSelection = 0x0;
    multSelection = (AliMultSelection*) gAOD->FindListObject("MultSelection");
    if (!multSelection){
        AliWarning("AliMultSelection object not found!");
    }
    else{
      gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator, kFALSE);
    }
    // event selection done in AliAnalysisTaskSE::Exec() --> this is not used
 
    if(fUsePID || fUsePIDNewTrial) {
        fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
        if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");

	fPIDCombined=new AliPIDCombined;
	fPIDCombined->SetDefaultTPCPriors();
    
    }

    
    fHistEventStats->Fill(1,gCentrality); //all events
    Bool_t isSelected = kTRUE;
    if(fUseOfflineTrigger)
        isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    if(isSelected) {
        fHistEventStats->Fill(2,gCentrality); //triggered events
        
        const AliVVertex *vertex = gAOD->GetPrimaryVertex();
        if(vertex) {
            Double32_t fCov[6];
            vertex->GetCovarianceMatrix(fCov);
            if(vertex->GetNContributors() > 0) {
                if(fCov[5] != 0) {
                    if(TMath::Abs(vertex->GetX()) < fVxMax) {
                        if(TMath::Abs(vertex->GetY()) < fVyMax) {
                            if(TMath::Abs(vertex->GetZ()) < fVzMax) {
                                fHistEventStats->Fill(3,gCentrality); //analyzed events
                                fHistVx->Fill(gCentrality,vertex->GetX());
                                fHistVy->Fill(gCentrality,vertex->GetY());
                                fHistVz->Fill(gCentrality,vertex->GetZ());
                                
                                if((gCentrality >= fCentralityPercentileMin) &&
                                   (gCentrality < fCentralityPercentileMax)) {
				  
				  fHistEventStats->Fill(4,gCentrality);
                                  
				  // check for pile-up event
				  const Int_t nTracks = gAOD->GetNumberOfTracks(); 
				  Int_t multEsd = ((AliAODHeader*)gAOD->GetHeader())->GetNumberOfESDTracks();
				  fHistGlobalvsESDBeforePileUpCuts->Fill(nTracks,multEsd);
                                   
                                    if(fpPb) {
				    fUtils->SetUseMVPlpSelection(kTRUE);
    				    //fUtils->SetUseOutOfBunchPileUp(kTRUE);
				    }
				    else if (fPbPb) {
				      fUtils->SetUseMVPlpSelection(kFALSE);
				      fUtils->SetUseOutOfBunchPileUp(kFALSE);
				      
				      if (fUseOutOfBunchPileUpCutsLHC15o){
					if (TMath::Abs(multSelection->GetMultiplicityPercentile("V0M") - multSelection->GetMultiplicityPercentile("CL1")) > 7.5) {
					  fHistEventStats->Fill(6,gCentrality);
					  return;
					}
					//const Int_t nTracks = gAOD->GetNumberOfTracks();
					//Int_t multEsd = ((AliAODHeader*)gAOD->GetHeader())->GetNumberOfESDTracks();
					Int_t multTPC = 0;
					for (Int_t it = 0; it < nTracks; it++) {
					  AliAODTrack* AODTrk = (AliAODTrack*)gAOD->GetTrack(it);
					  if (!AODTrk){ delete AODTrk; continue; }
					  if (AODTrk->TestFilterBit(128)) {multTPC++;}
					} // end of for (Int_t it = 0; it < nTracks; it++)
					
					if ((multEsd - fPileupLHC15oSlope*multTPC) > fPileupLHC15oOffset){
					  fHistEventStats->Fill(6,gCentrality);
					  return;
					}
					
					fHistGlobalvsESDAfterPileUpCuts->Fill(nTracks,multEsd);
				      }
				      
				      if (fUseOutOfBunchPileUpCutsLHC15oJpsi){
					if (TMath::Abs(multSelection->GetMultiplicityPercentile("V0M") - multSelection->GetMultiplicityPercentile("CL1")) > 7.5) {
					  fHistEventStats->Fill(6,gCentrality);
					  return;
					}
					Int_t ntrkTPCout = 0;
					for (int it = 0; it < gAOD->GetNumberOfTracks(); it++) {
					  AliAODTrack* AODTrk = (AliAODTrack*)gAOD->GetTrack(it);
					    if ((AODTrk->GetStatus() & AliAODTrack::kTPCout) && AODTrk->GetID() > 0)
					      ntrkTPCout++;
					}
					
					Double_t multVZERO =0; 
					AliVVZERO *vzero = (AliVVZERO*)gAOD->GetVZEROData();
					if(vzero) {
					  for(int ich=0; ich < 64; ich++)
					    multVZERO += vzero->GetMultiplicity(ich);
					}
					
					fHistV0MvsTPCoutBeforePileUpCuts->Fill(ntrkTPCout, multVZERO);
					
					if (multVZERO < (-2200 + 2.5*ntrkTPCout + 1.2e-5*ntrkTPCout*ntrkTPCout))  {
					  fHistEventStats->Fill(6, -1);
					  return;
					}
					fHistV0MvsTPCoutAfterPileUpCuts->Fill(ntrkTPCout, multVZERO);
					
				      }
				    }
				    
 				    if (fCheckPileUp){
				    if(fUtils->IsPileUpEvent(gAOD)){ 
				      fHistEventStats->Fill(6,gCentrality);
				      return;
				    }
				    }

    				    fHistEventStats->Fill(5,gCentrality); 	    
				    

				    fHistCent->Fill(gCentrality);
				    fHistCentbin->Fill(gCentrality);
				    
				    Double_t probTPC[AliPID::kSPECIES]={0.};
				    Double_t probTPCTOF[AliPID::kSPECIES]={0.};
				    Float_t nSigmaTPCOnly=0., nSigmaTPCNsigcomb=0., nSigmaTOFNsigcomb=0.;
				    Float_t combSquaredSigma=0.;
                                
                                    // Printf("There are %d tracks in this event", gAOD->GetNumberOfTracks());
                                    for (Int_t iTracks = 0; iTracks < gAOD->GetNumberOfTracks(); iTracks++) {
                                        AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(iTracks));
                                        if (!aodTrack) {
                                            Printf("ERROR: Could not receive track %d", iTracks);
                                            continue;
                                        }

                                        if(fExcludeSecondariesInMCrec){
                                            
                                            Int_t label = TMath::Abs(aodTrack->GetLabel());
                                            
                                            AliAODMCParticle *AODmcTrack = (AliAODMCParticle*) fArrayMC->At(label);
                                            
                                            if (!AODmcTrack->IsPhysicalPrimary())
                                                continue;
                                        }

					if(fExcludeElectronsInMCrec){

					    Int_t tracklabel = TMath::Abs(aodTrack->GetLabel());
                 			    AliAODMCParticle *trackAODMC = (AliAODMCParticle*) fArrayMC->At(tracklabel);
					    if(TMath::Abs(trackAODMC->GetPdgCode()) == 11) continue;
                 
 					}
					
					if (fExcludeInjectedSignals){
					  if (fRejectCheckGenName){
					    TString generatorName;
					    AliMCEvent* mcevent = dynamic_cast<AliMCEvent*>(MCEvent());
					    Int_t label = TMath::Abs(aodTrack->GetLabel());
					    Bool_t hasGenerator = mcevent->GetCocktailGenerator(label,generatorName);
					 
					    if((!hasGenerator) || (!generatorName.Contains(fGenToBeKept.Data())))
					      continue;
					    
					    //  Printf("mother =%d, generatorName=%s", label, generatorName.Data()); 
					  }
					}
				
					// AOD track cuts
                                        fHistTrackStats->Fill(gCentrality,aodTrack->GetFilterMap());
                                        //Printf("filterbit is: %i",GetFilterMap());
                                        if(!aodTrack->TestFilterBit(fAODtrackCutBit)) continue;
					
					AliAODPid* pidObj = aodTrack->GetDetPid();
                                        
                                        Float_t pt  = aodTrack->Pt();
                                        Float_t eta = aodTrack->Eta();
                                        Float_t phi = aodTrack->Phi();
                                        //Int_t numberofclustersTPC = aodTrack->GetNumberOfTPCClusters();
                                        Double_t chi2 = aodTrack->Chi2perNDF();
                                        Double_t xdca = aodTrack->DCA();
                                        Double_t zdca = aodTrack->ZAtDCA();
                                        Double_t charge = aodTrack->Charge();
 					Double_t y = log( ( sqrt(fMassParticleOfInterest*fMassParticleOfInterest + pt*pt*cosh(eta)*cosh(eta)) + pt*sinh(eta) ) / sqrt(fMassParticleOfInterest*fMassParticleOfInterest + pt*pt) );
					Double_t yfromAODTrack = aodTrack->Y(fMassParticleOfInterest);
					//Printf("y = %f, yfromAODTrack =%f", y, yfromAODTrack);
					
					if (fUseRapidity){
					  if ((y < fEtaMin) || ( y > fEtaMax))
					    continue;
					}
					else{
					  if( eta < fEtaMin || eta > fEtaMax) continue;
					}
					
                                        if( pt < fPtMin || pt > fPtMax) continue;

					Float_t mom = aodTrack->GetTPCmomentum();	
                                        if(fUsePID) {
					  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC); //firts check only TPC
					  UInt_t detUsed = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPC);

					  if (detUsed  == (UInt_t)fPIDCombined->GetDetectorMask()){
					    //Printf("fParticleOfInterest= %d", fParticleOfInterest);

					    nSigmaTPCOnly = fPIDResponse->NumberOfSigmasTPC(aodTrack,fParticleOfInterest);

					    if(pt < fPIDMomCut){
					      
					      if (fUsePIDnSigmaComb){
						if (TMath::Abs(nSigmaTPCOnly)<3.) {
						  if (charge>0)
						    fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
						  if (charge<0)
						    fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ()); 
						}
						else continue;
						hNSigmaCutApplied->Fill(TMath::Abs(nSigmaTPCOnly), mom, pt);
					      }
					      
					      else {
						//Printf(" detUsed = %d, mom = %f, pt =%f, probTPC[%d] =%f", detUsed, mom, pt, fParticleOfInterest, probTPC[fParticleOfInterest]);
						if (probTPC[fParticleOfInterest] > fBayesPIDThr) {
						  if (charge>0)
						    fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
						  if (charge<0)
						    fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ()); 
						}
						else continue;
						hBayesProbab->Fill(probTPC[fParticleOfInterest], mom, pt);
					      }
					      
					    }
					    
					    fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
					    detUsed = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPCTOF);
					    
					    if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()){
					      
					      if(!pidObj || pidObj->GetTOFsignal() > 99999)  continue;
					      
					      nSigmaTPCNsigcomb = fPIDResponse->NumberOfSigmasTPC(aodTrack,fParticleOfInterest);
					      nSigmaTOFNsigcomb = fPIDResponse->NumberOfSigmasTOF(aodTrack,fParticleOfInterest);
					      combSquaredSigma = TMath::Sqrt((nSigmaTPCNsigcomb*nSigmaTPCNsigcomb) + (nSigmaTOFNsigcomb*nSigmaTOFNsigcomb));
					      // printf("highpT : combSquaredSigma =%f", combSquaredSigma);

					      if (pt >= fPIDMomCut){
						if (fParticleOfInterest == (AliPID::kPion)){
						 
						  if (fUsePIDnSigmaComb){
						    if (pt <= 2.5){
						      if (TMath::Abs(combSquaredSigma)< 3){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      }
						      else continue;
						    }
						    else if (pt>2.5){
						      if (TMath::Abs(combSquaredSigma)< 2){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      }
						      else continue;
						    }
						    hNSigmaCutApplied->Fill(combSquaredSigma, mom, pt);
						  }
						  
						  else{
						    if (probTPCTOF[fParticleOfInterest] > fBayesPIDThr) {
						      if (charge>0)
							fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
						      if (charge<0)
							fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						    }
						    else continue;
						    hBayesProbab->Fill(probTPCTOF[fParticleOfInterest], mom, pt);
						  }
						} //end of pions 
						
						if (fParticleOfInterest == (AliPID::kKaon)){
						  
						  if (fUsePIDnSigmaComb){
						    if (pt <= 2.){
						      if (TMath::Abs(combSquaredSigma)<2.5){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      }
						      else continue;
						    }
						    if (pt > 2.) {
						      if (TMath::Abs(combSquaredSigma)< 1.5) {
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      }
						      else continue;
						    }

						    hNSigmaCutApplied->Fill(combSquaredSigma, mom, pt);
						  }
						  else{
						    if (probTPCTOF[fParticleOfInterest] > fBayesPIDThr) {
						      if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
						      if (charge<0)
							fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						    }
						    else continue;
						    hBayesProbab->Fill(probTPCTOF[fParticleOfInterest], mom, pt);
						  }
						} //end of kaons

						if (fParticleOfInterest == (AliPID::kProton)){
						  
						  if (fUsePIDnSigmaComb){
						    if (pt <= 3.){
						      if (TMath::Abs(combSquaredSigma)<3){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      }
						      else continue;
						    }
						    if ((pt > 3.)&&(pt <= 5.)){
						      if (TMath::Abs(combSquaredSigma)< 1.5){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      } else continue;
						    }
						    if (pt > 5.) {
						      if (TMath::Abs(combSquaredSigma)<1){
							if (charge>0)
							  fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
							if (charge<0)
							  fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						      } else continue;
						    }
						    hNSigmaCutApplied->Fill(combSquaredSigma, mom, pt);
						  }
						  else {
						    if (probTPCTOF[fParticleOfInterest] > fBayesPIDThr){
						      if (charge>0)
							fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ()); 
						      if (charge<0)
							fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
						    } else continue;
						    hBayesProbab->Fill(probTPCTOF[fParticleOfInterest], mom, pt);
						  }
						}// end of protons 	
					      } // higher pT PID
					    }
					  }
					  
					}//end of UsePID
                                        
                                        
                      if(fUsePIDNewTrial) {
                          
                          AliAODPid* pidObj = aodTrack->GetDetPid();
                          Bool_t isPartIDselected = kFALSE;
                          
                          Double_t nSigmaTPC = 0.;
                          
                          Double_t nSigmaTPCPions = 0.;
                          Double_t nSigmaTPCKaons = 0.;
                          Double_t nSigmaTPCProtons = 0.;
                          
                          Double_t nSigmaTOFPions = 0.;
                          Double_t nSigmaTOFKaons = 0.;
                          Double_t nSigmaTOFProtons = 0.;
                          
                          Double_t nSigmaTPCTOFPions = 0.;
                          Double_t nSigmaTPCTOFKaons = 0.;
                          Double_t nSigmaTPCTOFProtons = 0.;
                          
                          Double_t tofTime = -999., length = 999., tof = -999.;
                          Double_t c = TMath::C()*1.E-9;// m/ns
                          Double_t beta = -999.;
                          
                          fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTPC); //firts check only TPC
                          UInt_t detUsed = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPC);
                          
                          if (detUsed  == (UInt_t)fPIDCombined->GetDetectorMask()){
                              
                              nSigmaTPC = fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)fParticleOfInterest);
                              
                              nSigmaTPCPions   = fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kPion);
                              nSigmaTPCKaons   = fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kKaon);
                              nSigmaTPCProtons = fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kProton);
                          
                              
                              fHistdEdxVsPTPCbeforePID -> Fill(aodTrack->GetTPCmomentum()*aodTrack->Charge(),aodTrack->GetTPCsignal()); //aodTrack->P()*aodTrack->Charge()
                              fHistProbTPCvsPtbeforePID -> Fill(aodTrack->Pt(),probTPC[fParticleOfInterest]);
                              fHistNSigmaTPCvsPtbeforePID -> Fill(aodTrack->Pt(),nSigmaTPC);
                              
                              fPIDCombined->SetDetectorMask(AliPIDResponse::kDetTOF|AliPIDResponse::kDetTPC);
                              detUsed = fPIDCombined->ComputeProbabilities(aodTrack, fPIDResponse, probTPCTOF);
                              
                              if (detUsed == (UInt_t)fPIDCombined->GetDetectorMask()){
                                  
                                  if(!pidObj || pidObj->GetTOFsignal() > 99999)  continue;
                                  
                                  nSigmaTPCNsigcomb = fPIDResponse->NumberOfSigmasTPC(aodTrack,fParticleOfInterest);
                                  nSigmaTOFNsigcomb = fPIDResponse->NumberOfSigmasTOF(aodTrack,fParticleOfInterest);
                                  combSquaredSigma = TMath::Sqrt((nSigmaTPCNsigcomb*nSigmaTPCNsigcomb) + (nSigmaTOFNsigcomb*nSigmaTOFNsigcomb));
                          
                                  nSigmaTOFPions = fPIDResponse->NumberOfSigmasTOF(aodTrack,(AliPID::EParticleType)AliPID::kPion);
                                  nSigmaTOFKaons = fPIDResponse->NumberOfSigmasTOF(aodTrack,(AliPID::EParticleType)AliPID::kKaon);
                                  nSigmaTOFProtons = fPIDResponse->NumberOfSigmasTOF(aodTrack,(AliPID::EParticleType)AliPID::kProton);
                                  
                                  nSigmaTPCTOFPions = TMath::Sqrt(nSigmaTPCPions*nSigmaTPCPions + nSigmaTOFPions*nSigmaTOFPions);
                                  nSigmaTPCTOFKaons = TMath::Sqrt(nSigmaTPCKaons*nSigmaTPCKaons + nSigmaTOFKaons*nSigmaTOFKaons);
                                  nSigmaTPCTOFProtons = TMath::Sqrt(nSigmaTPCProtons*nSigmaTPCProtons + nSigmaTOFProtons*nSigmaTOFProtons);
                                
                                  
                                  
                                  if ((aodTrack->IsOn(AliAODTrack::kITSin)) && (aodTrack->IsOn(AliAODTrack::kTOFout)) ) {
                                      tofTime = aodTrack->GetTOFsignal();//in ps
                                      length = aodTrack->GetIntegratedLength();
                                      tof = tofTime*1E-3; // ns
                                      if (tof <= 0) {
                                          Printf("WARNING: track with negative TOF time found! Skipping this track for PID checks\n");
                                          continue;
                                      }
                                      if (length <= 0){
                                          // in old productions integrated track length is not stored in AODs -> need workaround
                                          Double_t exptime[10];
                                          aodTrack->GetIntegratedTimes(exptime);
                                          length = exptime[0]*c*1E-3/0.01; //assume electrons are relativistic (and add all multiplication factors)
                                          if (length <= 0){
                                              Printf("WARNING: track with negative length found!Skipping this track for PID checks\n");
                                              continue;
                                          }
                                      }
                                      length = length*0.01; // in meters
                                      tof = tof*c;
                                      beta = length/tof;
                                      
                                      fHistBetavsPTOFbeforePID ->Fill(aodTrack->P()*aodTrack->Charge(),beta);
                                      fHistNSigmaTOFvsPtbeforePID ->Fill(aodTrack->Pt(),nSigmaTOFNsigcomb);
                                      
                                      
                                      fHistProbTPCTOFvsPtbeforePID -> Fill(aodTrack->Pt(),probTPCTOF[fParticleOfInterest]);
                                      fHistBetaVsdEdXbeforePID->Fill(aodTrack->GetTPCsignal(),beta);
                                      fHistNSigmaTPCTOFvsPtbeforePID -> Fill(aodTrack->Pt(),combSquaredSigma);
                                      fHistNSigmaTPCTOFPbefPID ->Fill(nSigmaTPC,nSigmaTOFNsigcomb,aodTrack->P());
                                      
                                  }
                              
                                  
                                      if(pt < fPIDMomCut){
                                          
                                          if (fParticleOfInterest == (AliPID::kPion)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(nSigmaTPCOnly)<2.) && !(TMath::Abs(nSigmaTPCKaons)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)){
                                                      
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                  
                                                  }
                                                  
                                                  else continue;
                                                  
                                              }
                                              
                                          } //end of pions
                                          
                                          if (fParticleOfInterest == (AliPID::kKaon)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(nSigmaTPCOnly)<2.) && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCProtons)<3.)){
                                                      
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                  
                                                  }
                                                  
                                                  else continue;
                                                  
                                              }
                                              
                                          } //end of kaons
                                          
                                          if (fParticleOfInterest == (AliPID::kProton)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(nSigmaTPCOnly)<2.) && !(TMath::Abs(nSigmaTPCPions)<3.) && !(TMath::Abs(nSigmaTPCKaons)<3.)){
                                                      
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                  
                                                  }
                                                  
                                                  else continue;
                                                  
                                              }
                                              
                                          } //end of protons
                                          
                                      }
                                  
                                      if (pt >= fPIDMomCut){
                                          
                                          if (fParticleOfInterest == (AliPID::kPion)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(combSquaredSigma)<2.) && !(TMath::Abs(nSigmaTPCTOFKaons)<3.) && !(TMath::Abs(nSigmaTPCTOFProtons)<3.)){
                                                      
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                      
                                                  }
                                                 
                                                  else continue;
                                                  
                                              }
                                              
                                          } //end of pions
                                          
                                          if (fParticleOfInterest == (AliPID::kKaon)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(combSquaredSigma)<2.) && !(TMath::Abs(nSigmaTPCTOFPions)<3.) && !(TMath::Abs(nSigmaTPCTOFProtons)<3.)){
                                                      
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                  
                                                  }
                        
                                                  else continue;
                                              }
                                              
                                          } //end of kaons
                                          
                                          if (fParticleOfInterest == (AliPID::kProton)){
                                              
                                              if (fUsePIDnSigma){
                                                  
                                                  if ((TMath::Abs(combSquaredSigma)<2.) && !(TMath::Abs(nSigmaTPCTOFPions)<3.) && !(TMath::Abs(nSigmaTPCTOFKaons)<3.)){
                                                  
                                                      isPartIDselected = kTRUE;
                                                      
                                                      if (charge>0)
                                                          fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());
                                                      if (charge<0)
                                                          fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ());
                                                      
                                                  }
                                                  
                                                  else continue;
                                              }
                                              
                                          } //end of protons
                                          
                                      }
                                  
                                          if (fUsePIDnSigma){
                                              fHistNSigmaTOFvsPtafterPID ->Fill(aodTrack->Pt(),nSigmaTOFNsigcomb);
                                              fHistNSigmaTPCvsPtafterPID ->Fill(aodTrack->Pt(),nSigmaTPCOnly);
                                              fHistNSigmaTPCTOFvsPtafterPID ->Fill(aodTrack->Pt(),combSquaredSigma);
                                              fHistNSigmaTPCTOFPafterPID ->Fill(nSigmaTPCOnly,nSigmaTOFNsigcomb,aodTrack->P());  //++++++++++++++
                                          }
                                          
                                          //Fill QA after the PID
                                          fHistBetavsPTOFafterPID ->Fill(aodTrack->P()*aodTrack->Charge(),beta);
                                          fHistdEdxVsPTPCafterPID ->Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
                                          fHistBetaVsdEdXafterPID ->Fill(aodTrack->GetTPCsignal(),beta);
                                      
                                  }
                              }
                              // if no detector flag remove track
                              else {
                                  continue;
                              }
                          
                          if (isPartIDselected == kFALSE) continue;
                      }
                              
					    
					/*Float_t probMis = fPIDResponse->GetTOFMismatchProbability(aodTrack);                                            
					  
					  if (probMis < 0.01) { //if u want to reduce mismatch using also TPC						
					      
					  Double_t nSigmaPionTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kPion));
					  Double_t nSigmaKaonTPC   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kKaon));
					  Double_t nSigmaProtonTPC = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,AliPID::kProton));
					  
					  Double_t nSigmaPionTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kPion));
					  Double_t nSigmaKaonTOF   = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kKaon));
					  Double_t nSigmaProtonTOF = TMath::Abs(fPIDResponse->NumberOfSigmasTOF(aodTrack,AliPID::kProton));
                                          
					  
					  if (fParticleOfInterest == kPion){
                                          
					  if( pt > 0.2 && pt < 0.6 ){
					  if( nSigmaPionTPC > fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
					  continue;
					  }
					  else if(pt > 0.6){
					  if( nSigmaPionTPC > fPIDNSigma )
					  continue;
                                          
					  if( nSigmaPionTOF > fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
					  continue;
					  }
					  } //end of pion case
					  
					  
					  else if(fParticleOfInterest == kKaon){
					  
					  if( pt > 0.2 && pt < 0.4 ){
                                          
					  if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC > fPIDNSigma || nSigmaProtonTPC < fPIDNSigma )
					  continue;
					  }
					  else if(pt >= 0.4  && pt <= 2.5){
                                                    
					  if( nSigmaKaonTPC > fPIDNSigma )
					  continue;
                                          
					  if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF > fPIDNSigma || nSigmaProtonTOF < fPIDNSigma )
					  continue;
					  }
					  } //end of the kaon case
                                          
                                          
					  else if (fParticleOfInterest == kProton){
                    				                            
					  if( pt > 0.2 && pt < 0.6 ){
					  if( nSigmaPionTPC < fPIDNSigma || nSigmaKaonTPC < fPIDNSigma || nSigmaProtonTPC > fPIDNSigma )
					  continue;
					  }
                                          
					  else if(pt > 0.6  && pt < 4.0 ){
					  if( nSigmaProtonTPC > fPIDNSigma )
                                          
					  continue;
                                          
					  if( nSigmaPionTOF < fPIDNSigma || nSigmaKaonTOF < fPIDNSigma || nSigmaProtonTOF > fPIDNSigma )
                                          
					  continue;
					  }
					  }//end of the proton case
                      			  
					  }//end probability check
					*/

				       
					fHistPt->Fill(pt);
					fHistPtbin->Fill(pt);
					fHistPtCen->Fill(gCentrality,pt);
					fHistPhi->Fill(phi);
					fHistPhiCen->Fill(gCentrality,phi);
					fHistEta->Fill(eta);
					fHistEtaCen->Fill(gCentrality,eta);
					//fHistNClustersTPC->Fill(numberofclustersTPC);
					fHistChi2PerClusterTPC->Fill(chi2);
					fHistDCAToVertexXY->Fill(xdca);
					fHistDCAToVertexZ->Fill(zdca);
					fHistDCAToVertex2D->Fill(xdca,zdca);
					fHistEtaPhiCent->Fill(eta,phi,gCentrality);
					fHistPtEtaCent->Fill(pt,eta,gCentrality);
					fHistPtPhiCent->Fill(pt,phi,gCentrality);
					
					Double_t  dca[2] = {0.0,0.0};
					Double_t  cov[3] = {0.0,0.0,0.0};
					
					AliAODTrack copy(*aodTrack);
					
					if (fAODtrackCutBit==768){
					  if (aodTrack->TestFilterBit(256)){
					    if (!copy.PropagateToDCA(vertex,gAOD->GetMagneticField(),300.,dca,cov))
					      continue;
					  }
                                	}
				    
					else {
					  if (!copy.PropagateToDCA(vertex,gAOD->GetMagneticField(),300.,dca,cov))
					    continue;
					}
					
					
					if (charge>0){
					  fHistEtaPhiVertexPlus->Fill(phi,eta,vertex->GetZ());
					  //	  if (fUseRapidity)
					  //   fHistYPhiVertexPlus->Fill(phi,y,vertex->GetZ());  
					  if (fAODtrackCutBit==768){
					    if (fDCAext){
					      if (aodTrack->TestFilterBit(512))
						fHistDCAXYptchargedplus_ext->Fill(pt,eta,aodTrack->DCA());
					      else if (aodTrack->TestFilterBit(256))
						fHistDCAXYptchargedplus_ext->Fill(pt,eta,dca[0]);	
					    }
					    else if (!fDCAext) {
					      if (aodTrack->TestFilterBit(512))
						fHistDCAXYptchargedplus->Fill(pt,eta,aodTrack->DCA());
					      else if (aodTrack->TestFilterBit(256))
						fHistDCAXYptchargedplus->Fill(pt,eta,dca[0]);
					    }
					  }
					  else {
					    fHistDCAXYptchargedplus->Fill(pt,eta,dca[0]);
					    if (fDCAext)
					      fHistDCAXYptchargedplus_ext->Fill(pt,eta,dca[0]);
					  }
					}
                                        else if (charge<0){
					  fHistEtaPhiVertexMinus->Fill(phi,eta,vertex->GetZ());
					  //if (fUseRapidity)
					  // fHistYPhiVertexMinus->Fill(phi,y,vertex->GetZ()); 
					  if (fAODtrackCutBit==768){ 
					    if (fDCAext){
					      if (aodTrack->TestFilterBit(512))	
						fHistDCAXYptchargedminus_ext->Fill(pt,eta,aodTrack->DCA());
					      else if (aodTrack->TestFilterBit(256))	
						fHistDCAXYptchargedminus_ext->Fill(pt,eta,dca[0]);
					    }
					    else if (!fDCAext) {
					      if (aodTrack->TestFilterBit(512))
						fHistDCAXYptchargedminus->Fill(pt,eta,aodTrack->DCA());
					      else if (aodTrack->TestFilterBit(256))
						fHistDCAXYptchargedminus->Fill(pt,eta,dca[0]);
					    }
					  }
					  else {
					    fHistDCAXYptchargedminus->Fill(pt,eta,dca[0]);
					    if (fDCAext)
					      fHistDCAXYptchargedminus_ext->Fill(pt,eta,dca[0]);
						}	
					}

					 if (fMCrec){
                                        Int_t Label = TMath::Abs(aodTrack->GetLabel());
                                        AliAODMCParticle *trackAODMCforpdg = (AliAODMCParticle*) fArrayMC->At(Label);
                                        fHistPdg->Fill(trackAODMCforpdg->GetPdgCode());
                                        }
					
                                        nAcceptedTracks += 1;
                                    } //track loop
                                    fHistMultiplicity->Fill(gCentrality,nAcceptedTracks);
                                } //centrality check
                            }//Vz cut
                        }//Vy cut
                    }//Vx cut
                }//proper vertex resolution
            }//proper number of contributors
        }//vertex object valid
    }//triggered event
}

//________________________________________________________________________
void  AliAnalysisTaskAccCont::FinishTaskOutput(){
    //
}

//________________________________________________________________________
void AliAnalysisTaskAccCont::Terminate(Option_t *) {
    // Draw result to the screen
    // Called once at the end of the query
    
    // not implemented ...
    
}


