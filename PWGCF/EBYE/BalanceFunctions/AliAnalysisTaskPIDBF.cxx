#include "TChain.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h" 
#include "TH2D.h"                  
#include "TH3D.h"
#include "TArrayF.h"
#include "TF1.h"
#include "TRandom.h"
#include "TROOT.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliVTrack.h"
#include "AliAODVZERO.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVParticle.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h" 
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "AliTHn.h"    
#include "AliLog.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"

#include "AliEventPoolManager.h"           

#include "AliPID.h"                
#include "AliPIDResponse.h"        
#include "AliPIDCombined.h"        

#include "AliAnalysisTaskPIDBF.h"
#include "AliPidBFBase.h"
#include "AliAnalysisTaskTriggeredBF.h"
#include "TFile.h"
#include <iostream>
#include "AliEventCuts.h"
#include "AliAODTracklets.h"



// Analysis task for the PID BF code:
// Base Class : AliPidBFBase.cxx
// Noor Alam(VECC, Kolkata) : sk.noor.alam@cern.ch,noor1989phyalam@gmail.com
// Supervisor: Subhasis Chattopadhyay: sub.chattopadhyay@gmail.com
//[Special thanks to Michael Weber(m.weber@cern.ch) and Panos Christakoglou(panos.christakoglou@cern.ch)] 

using std::cout;
using std::endl;


const char * kpidName[]={"TPC","TOF","TPC-TOF"} ;
const char * kparticleName[]={"Pions","Kaons","Protons","Undefined"} ;
const char * kdetectorName[]={"ITS","TPC","TOF"} ;



ClassImp(AliAnalysisTaskPIDBF)

//________________________________________________________________________
AliAnalysisTaskPIDBF::AliAnalysisTaskPIDBF(const char *name) 
: AliAnalysisTaskSE(name),
  fDebugLevel(kFALSE),
  fBalance(0),
  fRunMixing(kFALSE),
  fRunMixingEventPlane(kFALSE),
  fMixingTracks(50000),
  fMixedBalance(0),
  fPoolMgr(0),
  fList(0),
  fListBF(0),
  fListBFM(0),
  fHistListPIDQA(0),
  QA_AliEventCuts(0),
  fEventCuts(0),
  fHistEventStats(0),
  fHistCentStats(0),
  fHistCentStatsUsed(0),
  fHistTriggerStats(0),
  fHistTrackStats(0),
  fHistVx(0),
  fHistVy(0),
  fHistVz(0),
  fHistMixEvents(0),
  fHistMixTracks(0),
  fHistEventPlane(0),
  fHistClus(0),
  fHistDCA(0),
  fHistChi2(0),
  fHistPt(0),
  fHistEta(0),
  fHistEta1D(0),
  fHistPhi1D(0),
  fHistRapidity(0),
  fHistPhi(0),
  fHistEtaPhiPos(0), 	       	 
  fHistEtaPhiNeg(0), 
  fHistPhiPos(0),
  fHistPhiNeg(0),
  fHistV0M(0),
  fHistRefTracks(0),
  fHistdEdxVsPTPCbeforePIDelectron(NULL),
  fHistNSigmaTPCvsPtbeforePIDelectron(NULL),
  fHistdEdxVsPTPCafterPIDelectron(NULL),
  fHistNSigmaTPCvsPtafterPIDelectron(NULL),
  fCentralityArrayBinsForCorrections(kCENTRALITY),
  fPIDResponse(0x0),
  fParticleType_(kPion_),
  fPIDSpeciesHisto(0),
  fUsePID(kFALSE),
  fPIDNSigma(3.0),
  fElectronRejection(kFALSE),
  fElectronOnlyRejection(kFALSE),
  fElectronRejectionNSigma(-1.),
  fElectronRejectionMinPt(0.),
  fElectronRejectionMaxPt(1000.),
  fCentralityEstimator("V0M"),
  fUseCentrality(kFALSE),
  fCentralityPercentileMin(0.), 
  fCentralityPercentileMax(5.),
  fImpactParameterMin(0.),
  fImpactParameterMax(20.),
  fMultiplicityEstimator("V0A"),
  fUseMultiplicity(kFALSE),
  fNumberOfAcceptedTracksMin(0),
  fNumberOfAcceptedTracksMax(10000),
  fHistNumberOfAcceptedTracks(0),
  fHistMultiplicity(0),
  fCheckFirstEventInChunk(kFALSE),
  fCheckPileUp(kFALSE),
  fCheckPrimaryFlagAOD(kFALSE),
  fUseMCforKinematics(kFALSE),
  fVxMax(0.3),
  fVyMax(0.3),
  fVzMax(10.),
  fnAODtrackCutBit(128),
  fPtMin(0.3),
  fPtMax(1.5),
  fEtaMin(-0.8),
  fEtaMax(0.8),
  fPtTOFMin(0.0),
  fPtTOFMax(10.0),
  fPtTPCMin(0.0),
  fPtTPCMax(0.3),
  fPhiMin(0.),
  fPhiMax(360.),
  fDCAxyCut(-1),
  fDCAzCut(-1),
  fTPCchi2Cut(-1),
  fNClustersTPCCut(-1),
  fTPCsharedCut(-1),
  fAcceptanceParameterization(0),
  fDifferentialV2(0),
  fUseFlowAfterBurner(kFALSE),
  fExcludeWeakDecaysInMC(kFALSE),
  fExcludeResonancesInMC(kFALSE),
  fExcludeElectronsInMC(kFALSE),
  fExcludeParticlesExtra(kFALSE),
  fUseMCPdgCode(kFALSE),
  fPDGCodeToBeAnalyzed(-1),
  fExcludeResonancePDGInMC(-1),
  fEventClass("EventPlane"), 
  fCustomBinning(""),
  fHistVZEROAGainEqualizationMap(0),
  fHistVZEROCGainEqualizationMap(0),
  fHistVZEROChannelGainEqualizationMap(0),
  fUtils(0),
  fRapidityInsteadOfEta(kFALSE),
  fTOFMisMatch(kFALSE),
  fMistMatchTOFProb(.01),
  fDetectorPID_(kTPCTOFpid_),
  fUseOfflineTrigger(kFALSE),
  fUseMultSelection(kFALSE),
  fHasTOFPID(kFALSE),
  fLowCut(0),
  fHighCut(0),
  fMultTOFLowCut(0),
  fMultTOFHighCut(0),
  fMultCentLowCut(0),
  fCutMultESDdif(0),
  fHistTOFPid(0),
  fHistdEdxPid(0)  
 {
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain

  //======================================================correction
  for (Int_t i=0; i<kCENTRALITY; i++){
    fHistCorrectionPlus[i] = NULL; 
    fHistCorrectionMinus[i] = NULL; 
    fCentralityArrayForCorrections[i] = -1.;
  }
  //=====================================================correction
  for(Int_t ipart=0;ipart<kProton_+1;ipart++)
    for(Int_t ipid=0;ipid<kProton_+1;ipid++)
      fnsigmas[ipart][ipid]=999.;

   for(Int_t ipart=0;ipart<kProton_+1;ipart++) fHasDoubleCounting[ipart]=kFALSE;


  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
  DefineOutput(4, TList::Class());
  DefineOutput(5, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskPIDBF::~AliAnalysisTaskPIDBF() {

  // delete fBalance; 
  // delete fShuffledBalance; 
  // delete fList;
  // delete fListBF; 
  // delete fListBFS;

  // delete fHistEventStats; 
  // delete fHistTrackStats; 
  // delete fHistVx; 
  // delete fHistVy; 
  // delete fHistVz; 

  // delete fHistClus;
  // delete fHistDCA;
  // delete fHistChi2;
  // delete fHistPt;
  // delete fHistEta;
  // delete fHistPhi;
  // delete fHistEtaPhiPos; 		 	 
  // delete fHistEtaPhiNeg;
  // delete fHistV0M;
}

//________________________________________________________________________
void AliAnalysisTaskPIDBF::UserCreateOutputObjects() {
  // Create histograms
  // Called once

  // global switch disabling the reference 
  // (to avoid "Replacing existing TH1" if several wagons are created in train)
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  if(!fBalance) {
    fBalance = new AliPidBFBase();
    fBalance->SetAnalysisLevel("ESD");
    fBalance->SetEventClass(fEventClass);
    //fBalance->SetNumberOfBins(-1,16);
    //fBalance->SetInterval(-1,-0.8,0.8,16,0.,1.6,15.);
  }
  

 if(fRunMixing) {
    if(!fMixedBalance) {
      fMixedBalance = new AliPidBFBase();
      fMixedBalance->SetAnalysisLevel("ESD");
      fMixedBalance->SetEventClass(fEventClass);
    }
  }

 TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  //QA list
  fList = new TList();
  fList->SetName("listQA");
  fList->SetOwner();

  //Balance Function list
  fListBF = new TList();
  fListBF->SetName("listBF");
  fListBF->SetOwner();


  if(fRunMixing) {
    fListBFM = new TList();
    fListBFM->SetName("listTriggeredBFMixed");
    fListBFM->SetOwner();
  }

  //PID QA list
  if(fUsePID || fElectronRejection) {
    fHistListPIDQA = new TList();
    fHistListPIDQA->SetName("listQAPID");
    fHistListPIDQA->SetOwner();
  }


  if(fUseMultSelection){
  QA_AliEventCuts = new TList();
  QA_AliEventCuts->SetName("QA_AliEventCuts");
  QA_AliEventCuts->SetOwner();
  


 fEventCuts = new AliEventCuts();

 fEventCuts->AddQAplotsToList(QA_AliEventCuts);
 fEventCuts->SetManualMode();
// fEventCuts->SetupLHC15o();

 fEventCuts->fMinVtz = -fVzMax;
 fEventCuts->fMaxVtz = fVxMax;
 fEventCuts->fMinCentrality = fCentralityPercentileMin;
 fEventCuts->fMaxCentrality = fCentralityPercentileMax;

}

// Pile up Function 

    fLowCut = new TF1("fLowCut", "[0]+[1]*x - 5.*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
    fHighCut = new TF1("fHighCut", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);


    fMultTOFLowCut = new TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);

    fMultTOFHighCut = new TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);


    fMultCentLowCut = new TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);


    fLowCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);
    fHighCut->SetParameters(0.0157497, 0.973488, 0.673612, 0.0290718, -0.000546728, 5.82749e-06);

    fMultTOFLowCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    fMultTOFHighCut->SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);

    fMultCentLowCut->SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);


  //Event stats.
  TString gCutName[7] = {"Total","Offline trigger",
                         "Vertex","Analyzed","sel. Centrality","Not1stEvInChunk","No Pile-Up"};
  fHistEventStats = new TH2F("fHistEventStats",
                             "Event statistics;;Centrality percentile;N_{events}",
                             7,0.5,7.5,220,-5,105);
  for(Int_t i = 1; i <= 7; i++)
    fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fList->Add(fHistEventStats);

  TString gCentName[13] = {"V0M","V0A","V0C","FMD","TRK","TKL","CL0","CL1","ZNA","ZPA","V0MvsFMD","TKLvsV0M","ZEMvsZDC"};
  fHistCentStats = new TH2F("fHistCentStats",
                             "Centrality statistics;;Cent percentile",
			    13,-0.5,12.5,220,-5,105);
  for(Int_t i = 1; i <= 13; i++){
    fHistCentStats->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
    //fHistCentStatsUsed->GetXaxis()->SetBinLabel(i,gCentName[i-1].Data());
  }
  fList->Add(fHistCentStats);

  fHistCentStatsUsed = new TH2F("fHistCentStatsUsed","Centrality statistics;;Cent percentile", 1,-0.5,0.5,220,-5,105);
  fHistCentStatsUsed->GetXaxis()->SetBinLabel(1,fCentralityEstimator.Data());
  fList->Add(fHistCentStatsUsed);

  fHistTriggerStats = new TH1F("fHistTriggerStats","Trigger statistics;TriggerBit;N_{events}",1025,0,1025);
  fList->Add(fHistTriggerStats);

  fHistTrackStats = new TH1F("fHistTrackStats","Event statistics;TrackFilterBit;N_{events}",16,0,16);
  fList->Add(fHistTrackStats);

  fHistNumberOfAcceptedTracks = new TH2F("fHistNumberOfAcceptedTracks",";N_{acc.};Centrality percentile;Entries",4001,-0.5,4000.5,220,-5,105);
  fList->Add(fHistNumberOfAcceptedTracks);

  fHistMultiplicity = new TH1F("fHistMultiplicity",";N_{ch.};Entries",30001,-0.5,30000.5);
  fList->Add(fHistMultiplicity);

  // Vertex distributions
  fHistVx = new TH1F("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVx);
  fHistVy = new TH1F("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm);Entries",100,-0.5,0.5);
  fList->Add(fHistVy);
  fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm);Centrality percentile;Entries",100,-20.,20.,220,-5,105);
  fList->Add(fHistVz);

  // Event Mixing
  fHistMixEvents = new TH2F("fHistMixEvents","Number of mixed events;Centrality percentile;N_{mix,evts}",101, 0, 101, 200, 0, 200);
  fList->Add(fHistMixEvents);
  fHistMixTracks = new TH2F("fHistMixTracks","Number of mixed tracks;Centrality percentile;N_{mix,trks}",101, 0, 101, 200, 0, fMixingTracks * 1.5);
  fList->Add(fHistMixTracks);




// Rapidity Histogram 

 if(fUsePID && fRapidityInsteadOfEta){
  fHistRapidity  = new TH2F("fHistRapidity","y distribution;y;Centrality percentile",200,-2,2,220,-5,105);

  fHistListPIDQA->Add(fHistRapidity);
}


// Pt vs NSigma  plot for TPC, TOF and TPC+TOF : MCAODrec

TString typeParticle;
if(fParticleType_ == kPion_) typeParticle ="Pion";
else if( fParticleType_ == kKaon_) typeParticle = "Kaon";
else if( fParticleType_ == kProton_) typeParticle = "Proton";

if(gAnalysisLevel == "MCAODrec" || gAnalysisLevel== "AOD"){
 if(fUsePID){

  fHistTOFPid=new TH2F("fHistTOFPid",Form(" PID signal from TOF of particle %s",typeParticle.Data()),1000,-fPtMax,fPtMax,1000, 0., 1.2);
  fHistdEdxPid=new TH2F("fHistdEdxPid",Form(" PID signal from TPC of particle %s",typeParticle.Data()),1000,-fPtMax,fPtMax, 1000, 0, 1000);

 fHistListPIDQA->Add(fHistTOFPid);
 fHistListPIDQA->Add(fHistdEdxPid);


  for(Int_t ipart=0;ipart<kProton_+1;ipart++){
    for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
      Double_t miny=-30;
      Double_t maxy=30;
      if(ipid==kTPCTOFpid_){miny=0;maxy=50;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigma_%d_%d",ipart,ipid),Form("n#sigma %s %s",kparticleName[ipart],kpidName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kparticleName[ipart],kpidName[ipid]));
      fHistListPIDQA->Add(fHistoNSigma);
    }
  }

   //nsigmaRec plot
  for(Int_t ipart=0;ipart<kProton_+1;ipart++){
    for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==kTPCTOFpid_){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaRec_%d_%d",ipart,ipid),
                                  Form("n#sigma for reconstructed %s %s",kparticleName[ipart],kpidName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kparticleName[ipart],kpidName[ipid]));
      fHistListPIDQA->Add(fHistoNSigma);
    }
  }

// Double Count 

 for(Int_t ipart=0;ipart<kProton_+1;ipart++){
    for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
      Double_t miny=-10;
      Double_t maxy=10;
      if(ipid==kTPCTOFpid_){miny=0;maxy=20;}
      TH2F *fHistoNSigma=new TH2F(Form("NSigmaDC_%d_%d",ipart,ipid),
                                  Form("n#sigma for double counting %s %s",kparticleName[ipart],kpidName[ipid]),200,0,10,500,miny,maxy);
      fHistoNSigma->GetXaxis()->SetTitle("P_{T} (GeV / c)");
      fHistoNSigma->GetYaxis()->SetTitle(Form("n#sigma %s %s",kparticleName[ipart],kpidName[ipid]));
      fHistListPIDQA->Add(fHistoNSigma);
    }
  }

//PID signal plot
  for(Int_t idet=0;idet<kNDetectors_D;idet++){
    for(Int_t ipart=0;ipart<kProton_+1;ipart++){
      Double_t maxy=500;
      if(idet==kTOF_D)maxy=1.1;
      TH2F *fHistoPID=new TH2F(Form("PID_%d_%d",idet,ipart),Form("%s signal - %s",kdetectorName[idet],kparticleName[ipart]),200,0,10,500,-maxy,maxy);
      fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
      fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kdetectorName[idet]));
      fHistListPIDQA->Add(fHistoPID);
    }
  }


  //PID signal plot, before PID cut
  for(Int_t idet=0;idet<kNDetectors_D;idet++){
    Double_t maxy=500;
    if(idet==kTOF_D)maxy=1.1;
    TH2F *fHistoPID=new TH2F(Form("PIDAll_%d",idet),Form("%s signal",kdetectorName[idet]),200,0,10,500,-maxy,maxy);
    fHistoPID->GetXaxis()->SetTitle("P (GeV / c)");
    fHistoPID->GetYaxis()->SetTitle(Form("%s signal",kdetectorName[idet]));
    fHistListPIDQA->Add(fHistoPID);
  }


fPIDSpeciesHisto=new TH1D("fPIDSpeciesHisto","Histogram of PID Species",5,1,6);
   fHistListPIDQA->Add(fPIDSpeciesHisto);

}


}


  //TPC vs VZERO multiplicity


  //Event plane
  fHistEventPlane = new TH2F("fHistEventPlane",";#Psi_{2} [deg.];Centrality percentile;Counts",100,0,360.,220,-5,105);
  fList->Add(fHistEventPlane);

  // QA histograms
  fHistClus = new TH2F("fHistClus","# Cluster (TPC vs. ITS)",10,0,10,200,0,200);
  fList->Add(fHistClus);
  fHistChi2 = new TH2F("fHistChi2","Chi2/NDF distribution;#chi^{2}/ndf;Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistChi2);
  fHistDCA  = new TH2F("fHistDCA","DCA (xy vs. z)",400,-5,5,400,-5,5); 
  fList->Add(fHistDCA);
  fHistPt   = new TH2F("fHistPt","p_{T} distribution;p_{T} (GeV/c);Centrality percentile",200,0,10,220,-5,105);
  fList->Add(fHistPt);
  fHistEta  = new TH2F("fHistEta","#eta distribution;#eta;Centrality percentile",200,-2,2,220,-5,105);
  fList->Add(fHistEta);
  fHistEta1D  = new TH1F("fHistEta1D","#eta distribution 1d ;#eta",200,-2,2);
  fList->Add(fHistEta1D);
  fHistPhi1D  = new TH1F("fHistPhi1D","#phi distribution 1d ;#phi",200,0,2*TMath::Pi());
  fList->Add(fHistPhi1D);

  fHistPhi  = new TH2F("fHistPhi","#phi distribution;#phi (rad);Centrality percentile",200,0.0,2.*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhi);
  fHistEtaPhiPos  = new TH3F("fHistEtaPhiPos","#eta-#phi distribution (+);#eta;#phi (rad);Centrality percentile",40,-1.6,1.6,72,0.,2.*TMath::Pi(),220,-5,105); 		 	 
  fList->Add(fHistEtaPhiPos); 			 
  fHistEtaPhiNeg  = new TH3F("fHistEtaPhiNeg","#eta-#phi distribution (-);#eta;#phi (rad);Centrality percentile",40,-1.6,1.6,72,0.,2.*TMath::Pi(),220,-5,105); 	       	 
  fList->Add(fHistEtaPhiNeg);
//  fHistPhiBefore  = new TH2F("fHistPhiBefore","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
//  fList->Add(fHistPhiBefore);
//  fHistPhiAfter  = new TH2F("fHistPhiAfter","#phi distribution;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
//  fList->Add(fHistPhiAfter);
  fHistPhiPos  = new TH2F("fHistPhiPos","#phi distribution for positive particles;#phi;Centrality percentile",200,0.,2*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiPos);
  fHistPhiNeg  = new TH2F("fHistPhiNeg","#phi distribution for negative particles;#phi;Centrality percentile",200,0.,2.*TMath::Pi(),220,-5,105);
  fList->Add(fHistPhiNeg);
  fHistV0M  = new TH2F("fHistV0M","V0 Multiplicity C vs. A",500, 0, 20000, 500, 0, 20000);
  fList->Add(fHistV0M);
  TString gRefTrackName[6] = {"tracks","tracksPos","tracksNeg","tracksTPConly","clusITS0","clusITS1"};
  fHistRefTracks  = new TH2F("fHistRefTracks","Nr of Ref tracks/event vs. ref track estimator;;Nr of tracks",6, 0, 6, 400, 0, 20000);
  for(Int_t i = 1; i <= 6; i++)
  fHistRefTracks->GetXaxis()->SetBinLabel(i,gRefTrackName[i-1].Data());
  fList->Add(fHistRefTracks);



  // Balance function histograms
  // Initialize histograms if not done yet (including the custom binning)
  if(!fBalance->GetHistNp()){
   AliInfo("Histograms not yet initialized! --> Will be done now");
    fBalance->SetCustomBinning(fCustomBinning);
    fBalance->InitHistograms();
  }


  if(fRunMixing) {
    if(!fMixedBalance->GetHistNp()) {
    AliInfo("Histograms (mixing) not yet initialized! --> Will be done now");
      fMixedBalance->SetCustomBinning(fCustomBinning);
      fMixedBalance->InitHistograms();
    }
  }

  // QA histograms for different cuts
  fList->Add(fBalance->GetQAHistHBTbefore());
  fList->Add(fBalance->GetQAHistHBTafter());
  fList->Add(fBalance->GetQAHistPhiStarHBTbefore());
  fList->Add(fBalance->GetQAHistPhiStarHBTafter());
  fList->Add(fBalance->GetQAHistConversionbefore());
  fList->Add(fBalance->GetQAHistConversionafter());
  fList->Add(fBalance->GetQAHistPsiMinusPhi());
  fList->Add(fBalance->GetQAHistResonancesBefore());
  fList->Add(fBalance->GetQAHistResonancesRho());
  fList->Add(fBalance->GetQAHistResonancesK0());
  fList->Add(fBalance->GetQAHistResonancesLambda());
  fList->Add(fBalance->GetQAHistQbefore());
  fList->Add(fBalance->GetQAHistQafter());

  //for(Int_t a = 0; a < ANALYSIS_TYPES; a++){
  fListBF->Add(fBalance->GetHistNp());
  fListBF->Add(fBalance->GetHistNn());
  fListBF->Add(fBalance->GetHistNpn());
  fListBF->Add(fBalance->GetHistNnn());
  fListBF->Add(fBalance->GetHistNpp());
  fListBF->Add(fBalance->GetHistNnp());


  if(fRunMixing) {
    fListBFM->Add(fMixedBalance->GetHistNp());
    fListBFM->Add(fMixedBalance->GetHistNn());
    fListBFM->Add(fMixedBalance->GetHistNpn());
    fListBFM->Add(fMixedBalance->GetHistNnn());
    fListBFM->Add(fMixedBalance->GetHistNpp());
    fListBFM->Add(fMixedBalance->GetHistNnp());
  }
  //}


  // Event Mixing
  if(fRunMixing){
    Int_t trackDepth = fMixingTracks; 
    Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
    
    // centrality bins
    Double_t* centbins = NULL;
    Int_t nCentralityBins;
    if(fBalance->IsUseVertexBinning()){
      centbins = fBalance->GetBinning(fBalance->GetBinningString(), "centralityVertex", nCentralityBins);
    }
    else{
      centbins = fBalance->GetBinning(fBalance->GetBinningString(), "centrality", nCentralityBins);
    }
    
    // multiplicity bins
    Double_t* multbins = NULL;
    Int_t nMultiplicityBins;
    multbins = fBalance->GetBinning(fBalance->GetBinningString(), "multiplicity", nMultiplicityBins);
    
    // Zvtx bins
    Double_t* vtxbins = NULL; 
    Int_t nVertexBins;
    if(fBalance->IsUseVertexBinning()){
      vtxbins = fBalance->GetBinning(fBalance->GetBinningString(), "vertexVertex", nVertexBins);
    }
    else{
      vtxbins = fBalance->GetBinning(fBalance->GetBinningString(), "vertex", nVertexBins);
    }

    // Event plane angle (Psi) bins
    Double_t* psibins = NULL;
    Int_t nPsiBins; 
    psibins = fBalance->GetBinning(fBalance->GetBinningString(), "eventPlane", nPsiBins);

  
    // run the event mixing also in bins of event plane (statistics!)
    if(fRunMixingEventPlane){
      if(fEventClass=="Multiplicity"){
	if(multbins && vtxbins && psibins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins, nPsiBins, psibins);
	}
      }
      else{
	if(centbins && vtxbins && psibins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins, nPsiBins, psibins);
	}
      }
    }
    else{
      if(fEventClass=="Multiplicity"){
	if(multbins && vtxbins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nMultiplicityBins, multbins, nVertexBins, vtxbins);
	}
      }
      else{
	if(centbins && vtxbins){
	  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centbins, nVertexBins, vtxbins);
	}
      }
    }
    
    if(centbins) delete [] centbins; 
    if(multbins) delete [] multbins; 
    if(vtxbins)  delete [] vtxbins; 
    if(psibins)  delete [] psibins; 

    // set minimum values for track depth, fraction, and number of events
    fPoolMgr->SetTargetValues(fMixingTracks, 0.1, 5);
    
    // check pool manager
    if(!fPoolMgr){
      AliError("Event Mixing required, but Pool Manager not initialized...");
      return;
    }
  }
  

  // for electron rejection only TPC nsigma histograms
  if(fElectronRejection) {
 
    fHistdEdxVsPTPCbeforePIDelectron = new TH2D ("dEdxVsPTPCbeforeelectron","dEdxVsPTPCbeforeelectron", 1000, -10.0, 10.0, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCbeforePIDelectron);
    
    fHistNSigmaTPCvsPtbeforePIDelectron = new TH2D ("NSigmaTPCvsPtbeforeelectron","NSigmaTPCvsPtbeforeelectron", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtbeforePIDelectron);
    
    fHistdEdxVsPTPCafterPIDelectron = new TH2D ("dEdxVsPTPCafterelectron","dEdxVsPTPCafterelectron", 1000, -10, 10, 1000, 0, 1000); 
    fHistListPIDQA->Add(fHistdEdxVsPTPCafterPIDelectron);

    fHistNSigmaTPCvsPtafterPIDelectron = new TH2D ("NSigmaTPCvsPtafterelectron","NSigmaTPCvsPtafterelectron", 1000, -10, 10, 1000, 0, 500); 
    fHistListPIDQA->Add(fHistNSigmaTPCvsPtafterPIDelectron); 
  }
  //====================PID========================//

  // Post output data.
  PostData(1, fList);
  PostData(2, fListBF);
  if(fRunMixing) PostData(3, fListBFM);
  if(fUsePID || fElectronRejection) PostData(4, fHistListPIDQA);       //PID

   AliInfo("Finished setting up the Output");

  TH1::AddDirectory(oldStatus);
  
  fUtils = new AliAnalysisUtils();
}


//________________________________________________________________________
void AliAnalysisTaskPIDBF::SetInputCorrection(TString filename, 
					      Int_t nCentralityBins, 
					      Double_t *centralityArrayForCorrections) 
{
  //Open files that will be used for correction
  fCentralityArrayBinsForCorrections = nCentralityBins;
  for (Int_t i=0; i<nCentralityBins; i++)
    fCentralityArrayForCorrections[i] = centralityArrayForCorrections[i];

  // No file specified -> Abort
  if(!filename.Contains(".root")) {
    AliFatal(Form("No correction file specified (= %s) but correction requested ==> ABORT",filename.Data()));
    return;
  }

  //Open the input file
  TFile *f = TFile::Open(filename);
  if(!f->IsOpen()) {
    AliFatal(Form("File %s not found but correction requested ==> ABORT",filename.Data()));
    return;
  }
    
  //TString listEffName = "";
  for (Int_t iCent = 0; iCent < fCentralityArrayBinsForCorrections-1; iCent++) {    
    //Printf("iCent %d:",iCent);    
    TString histoName = "fHistCorrectionPlus";
    histoName += Form("%d-%d",(Int_t)(fCentralityArrayForCorrections[iCent]),(Int_t)(fCentralityArrayForCorrections[iCent+1]));
    fHistCorrectionPlus[iCent]= dynamic_cast<TH3F *>(f->Get(histoName.Data()));
    if(!fHistCorrectionPlus[iCent]) {
      AliFatal(Form("fHist %s not found but correction requested ==> ABORT",histoName.Data()));
      return;
    }
    
    histoName = "fHistCorrectionMinus";
    histoName += Form("%d-%d",(Int_t)(fCentralityArrayForCorrections[iCent]),(Int_t)(fCentralityArrayForCorrections[iCent+1]));
    fHistCorrectionMinus[iCent] = dynamic_cast<TH3F *>(f->Get(histoName.Data())); 
    if(!fHistCorrectionMinus[iCent]) {
      AliFatal(Form("fHist %s not found but correction requested ==> ABORT",histoName.Data()));
      return; 
    }
  }//loop over centralities: ONLY the PbPb case is covered
}

//________________________________________________________________________
void AliAnalysisTaskPIDBF::UserExec(Option_t *) {
  // Main loop
  // Called for each event


//cout<<" TPC pt min "<<fPtTPCMin<<'\t'<<"TPC Pt Max"<<fPtTPCMax<<'\t'<<"TOF Pt min"<<fPtTOFMin<<'\t'<<"TOF Pt Max "<<fPtTOFMax<<endl;

Float_t gCalculateCentrality =-1;


  TString gAnalysisLevel = fBalance->GetAnalysisLevel();
  Int_t gNumberOfAcceptedTracks = 0;
  Double_t lMultiplicityVar     = -999.; //-1
  Double_t gReactionPlane       = -1.; 
  Float_t bSign = 0.;
  
  // get the event (for generator level: MCEvent())
  AliVEvent* eventMain = NULL;
  if(gAnalysisLevel == "MC") {
    eventMain = dynamic_cast<AliVEvent*>(MCEvent()); 
  }
  else{
    eventMain = dynamic_cast<AliVEvent*>(InputEvent());     
    // for HBT like cuts need magnetic field sign
    bSign = (eventMain->GetMagneticField() > 0) ? 1 : -1;
  }
  if(!eventMain) {
    AliError("eventMain not available");
    return;
  }
 
  // PID Response task active?
  if(fUsePID || fElectronRejection) {
    fPIDResponse = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetPIDResponse();
    if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");
  }


   if (!fEventCuts->AcceptEvent(eventMain)) {
  PostData(5, QA_AliEventCuts);
  return;
}

IsProperVertexPileUp(eventMain);

Float_t vtxX = -999;
Float_t vtxY = -999;
Float_t vtxZ = -999;

  const AliAODVertex *trkVtx = (AliAODVertex*)eventMain->GetPrimaryVertex();

    if (!trkVtx || trkVtx->GetNContributors()<=0)
    return ;

    vtxX = trkVtx->GetX();
    vtxY = trkVtx->GetY();
    vtxZ = trkVtx->GetZ();

    if(vtxX > fVxMax ) return;
    if(vtxY > fVyMax ) return;
    if(vtxZ > fVzMax ) return;
    fHistEventStats->Fill(4,-1);//analyzed events

  if((lMultiplicityVar = IsEventAccepted(eventMain)) < 0){ 
    return;
  }

   fHistVx->Fill(vtxX);
   fHistVy->Fill(vtxY);
   fHistVz->Fill(vtxZ,lMultiplicityVar);



  // get the reaction plane
  if(fEventClass != "Multiplicity" && gAnalysisLevel!="AODnano") {
    gReactionPlane = GetEventPlane(eventMain);
    fHistEventPlane->Fill(gReactionPlane,lMultiplicityVar);
    if(gReactionPlane < 0){
      return;
    }
  }
  
  // get the accepted tracks in main event
  TObjArray *tracksMain = GetAcceptedTracks(eventMain,lMultiplicityVar,gReactionPlane);
  gNumberOfAcceptedTracks = tracksMain->GetEntriesFast();

//cout<<"Number of Accepted tracks "<<gNumberOfAcceptedTracks<<endl;

  //multiplicity cut (used in pp)
  fHistNumberOfAcceptedTracks->Fill(gNumberOfAcceptedTracks,lMultiplicityVar);
  fHistMultiplicity->Fill(gNumberOfAcceptedTracks);

  
  // Event mixing 
  if (fRunMixing)
    {
      // 1. First get an event pool corresponding in mult (cent) and
      //    zvertex to the current event. Once initialized, the pool
      //    should contain nMix (reduced) events. This routine does not
      //    pre-scan the chain. The first several events of every chain
      //    will be skipped until the needed pools are filled to the
      //    specified depth. If the pool categories are not too rare, this
      //    should not be a problem. If they are rare, you could lose`
      //    statistics.
      
      // 2. Collect the whole pool's content of tracks into one TObjArray
      //    (bgTracks), which is effectively a single background super-event.
      
      // 3. The reduced and bgTracks arrays must both be passed into
      //    FillCorrelations(). Also nMix should be passed in, so a weight
      //    of 1./nMix can be applied.
      
      AliEventPool* pool = fPoolMgr->GetEventPool(lMultiplicityVar, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane);
      
      if (!pool){
	AliFatal(Form("No pool found for centrality = %f, zVtx = %f, psi = %f", lMultiplicityVar, eventMain->GetPrimaryVertex()->GetZ(),gReactionPlane));
      }
      else{
	
	//pool->SetDebug(1);

	if (pool->IsReady()){ 
	  
	  
	  Int_t nMix = pool->GetCurrentNEvents();
	  //cout << "nMix = " << nMix << " tracks in pool = " << pool->NTracksInPool() << endl;
	  
	  fHistMixEvents->Fill(lMultiplicityVar, nMix);
	  fHistMixTracks->Fill(lMultiplicityVar, pool->NTracksInPool());

	  // Fill mixed-event histos here  
	  for (Int_t jMix=0; jMix<nMix; jMix++) 
	    {
	      TObjArray* tracksMixed = pool->GetEvent(jMix);
	      fMixedBalance->CalculateBalance(gReactionPlane,tracksMain,tracksMixed,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
	    }
	}
	
	// Update the Event pool
	pool->UpdatePool(tracksMain);
	//pool->PrintInfo();
	
      }//pool NULL check  
    }//run mixing
  
  // calculate balance function
  fBalance->CalculateBalance(gReactionPlane,tracksMain,NULL,bSign,lMultiplicityVar,eventMain->GetPrimaryVertex()->GetZ());
  
  // calculate shuffled balance function
}      


//_______________________________________________________________________  
void AliAnalysisTaskPIDBF::IsProperVertexPileUp(AliVEvent *event){

   //Centrality
    Float_t v0Centr    = -100.;
    Float_t cl1Centr   = -100.;
    Float_t cl0Centr   = -100.;

    AliMultSelection* MultSelection = 0x0;
    MultSelection = (AliMultSelection*) event->FindListObject("MultSelection");
    if( !MultSelection) {
        AliWarning("AliMultSelection object not found!");
        return;
    } else {
        v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
        cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
        cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
    }

    if (v0Centr >= 90. || v0Centr < 0)
        return;


    Int_t nITSClsLy0 = event->GetNumberOfITSClusters(0);
    Int_t nITSClsLy1 = event->GetNumberOfITSClusters(1);
    Int_t nITSCls = nITSClsLy0 + nITSClsLy1;

/*
    AliAODTracklets* aodTrkl = (AliAODTracklets*)event->GetTracklets();
    Int_t nITSTrkls = aodTrkl->GetNumberOfTracklets();
*/

    const Int_t nTracks = event->GetNumberOfTracks();
    Int_t multEsd = ((AliAODHeader*)event->GetHeader())->GetNumberOfESDTracks();

    Int_t multTrk = 0;
    Int_t multTrkBefC = 0;
    Int_t multTrkTOFBefC = 0;
    Int_t multTPC = 0;
    Int_t multTPCout = 0;

    for (Int_t it = 0; it < nTracks; it++) {

     AliAODTrack* aodTrk = (AliAODTrack*)event->GetTrack(it);

        if (!aodTrk){
            delete aodTrk;
            continue;
        }

        if (aodTrk->GetFlags()&AliESDtrack::kTPCout)
            multTPCout++;

        if (aodTrk->TestFilterBit(32)){

            multTrkBefC++;

            if ( TMath::Abs(aodTrk->GetTOFsignalDz()) <= 10. && aodTrk->GetTOFsignal() >= 12000. && aodTrk->GetTOFsignal() <= 25000.)
                multTrkTOFBefC++;

            if ((TMath::Abs(aodTrk->Eta()) < fEtaMax) && (aodTrk->GetTPCNcls() >= fNClustersTPCCut) && (aodTrk->Pt() >= fPtMin) && (aodTrk->Pt() < fPtMax))
                multTrk++;

        }

         if (aodTrk->TestFilterBit(128))
            multTPC++;

} // track loop end 

    Float_t multTPCn = multTPC;
    Float_t multEsdn = multEsd;
    Float_t multESDTPCDif = multEsdn - multTPCn*3.38;

    //cout<<" multTPC"<<multTPC<<""<<"multEsd"<<multEsd<<""<<"multESDTPCDif"<<multESDTPCDif<<endl;

    AliAODVZERO* aodV0 = (AliAODVZERO*) event->GetVZEROData();
//    AliVVZERO* aodV0 = event->GetVZEROData();
    Float_t multV0a = aodV0->GetMTotV0A();
    Float_t multV0c = aodV0->GetMTotV0C();
    Float_t multV0Tot = multV0a + multV0c;

//   cout<<" multV0a"<<multV0a<<""<<"multV0c"<<multV0c<<""<<"multV0Tot"<<multV0Tot<<endl;


 //new vertex selection
      const AliAODVertex* vtTrc = (AliAODVertex*) event->GetPrimaryVertex();
      const AliAODVertex* vtSPD = (AliAODVertex*) event->GetPrimaryVertexSPD();

  //     const AliVVertex *vtTrc = event->GetPrimaryVertex();
 //     const AliVVertex* vtSPD = event->GetPrimaryVertexSPD();

    if (vtTrc->GetNContributors() < 2 || vtSPD->GetNContributors()<1)
        return; // one of vertices is missing

    double covTrc[6], covSPD[6];
    vtTrc->GetCovarianceMatrix(covTrc);
    vtSPD->GetCovarianceMatrix(covSPD);

    double dz = vtTrc->GetZ() - vtSPD->GetZ();

    double errTot = TMath::Sqrt(covTrc[5]+covSPD[5]);
    double errTrc = TMath::Sqrt(covTrc[5]);
    double nsigTot = dz/errTot;
    double nsigTrc = dz/errTrc;

    if (TMath::Abs(dz)>0.2 || TMath::Abs(nsigTot)>10 || TMath::Abs(nsigTrc)>20)
        return; // bad vertexing


    // vertex cut from old selection
    TString vtxTyp = vtSPD->GetTitle();
    Double_t zRes = TMath::Sqrt(covSPD[5]);
    if ((vtxTyp.Contains("vertexer:Z")) && (zRes>0.25) && (vtSPD->GetNContributors() < 20))
        return;


    if(fCheckPileUp){
      
     if (cl0Centr < fLowCut->Eval(v0Centr)) 
        return;

        if (cl0Centr > fHighCut->Eval(v0Centr))
            return;


        if (multESDTPCDif > fCutMultESDdif)
            return;

        if (Float_t(multTrkTOFBefC) < fMultTOFLowCut->Eval(Float_t(multTrkBefC)))
            return;

        if (Float_t(multTrkTOFBefC) > fMultTOFHighCut->Eval(Float_t(multTrkBefC)))
            return;
 
         //Short_t isPileup = event->IsPileupFromSPD(3);
         Short_t isPileup = event->IsPileupFromSPD(5,0.8,3.,2.,.5);
        if (isPileup != 0)
            return; 

        if (((AliAODHeader*)event->GetHeader())->GetRefMultiplicityComb08() < 0)
            return;

        //new function for 2015 to remove incomplete events
        if (event->IsIncompleteDAQ())
            return;


        //new cut to remove outliers
        if (Float_t(multTrk) < fMultCentLowCut->Eval(v0Centr))
            return;
    

    }


}

//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::IsEventAccepted(AliVEvent *event){
  // Checks the Event cuts
  // Fills Event statistics histograms
 
  Bool_t isTriggerselected = kTRUE;
  Float_t gRefMultiplicity = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  AliMCEvent *mcevent = dynamic_cast<AliMCEvent*>(event);

  fHistEventStats->Fill(1,gRefMultiplicity); //all events

  // check first event in chunk (is not needed for new reconstructions)
  if(fCheckFirstEventInChunk){
    if(fUtils->IsFirstEventInChunk(event)) 
      return -1.;
    fHistEventStats->Fill(6,gRefMultiplicity); 
  }

  // Event trigger bits
  fHistTriggerStats->Fill(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
 
   if(fUseOfflineTrigger) {
isTriggerselected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7);
//if(isTriggerselected) cout<<" in side Trigger kINT7 ;"<<endl;
}

else {
 isTriggerselected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral);
}

 
  if(isTriggerselected) {
    fHistEventStats->Fill(2,gRefMultiplicity); //triggered events
      

		    gRefMultiplicity = GetRefMultiOrCentrality(event);

		  
		  // take only events inside centrality class
		  if(fUseCentrality) {
		    if((gRefMultiplicity > fCentralityPercentileMin) && (gRefMultiplicity < fCentralityPercentileMax)){
		      
		      fHistEventStats->Fill(5,gRefMultiplicity); //events with correct centrality
                       
		      return gRefMultiplicity;	
		    }//centrality class
		  }
		  // take events only within the same multiplicity class
		  else if(fUseMultiplicity) {
		    //if(fDebugLevel) 
		    //Printf("N(min): %.0f, N(max): %.0f - N(ref): %.0f",fNumberOfAcceptedTracksMin,
		    //fNumberOfAcceptedTracksMax,gRefMultiplicity);

		    if((gRefMultiplicity > fNumberOfAcceptedTracksMin) && (gRefMultiplicity < fNumberOfAcceptedTracksMax)) {
		      fHistEventStats->Fill(5,gRefMultiplicity); //events with correct multiplicity
		      return gRefMultiplicity;
		    }
		  }//multiplicity range
    }//trigger
  
  // in all other cases return -1 (event not accepted)
  return -1;
}


//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::GetRefMultiOrCentrality(AliVEvent *event){
    // Checks the Event cuts
    // Fills Event statistics histograms

  Float_t gCentrality = -1.;
  TString gAnalysisLevel = fBalance->GetAnalysisLevel();


  // calculate centrality always (not only in centrality mode)
  if(gAnalysisLevel == "AOD"|| gAnalysisLevel == "MCAOD" || gAnalysisLevel == "MCAODrec" ) { //centrality in AOD header  //++++++++++++++

    if(fUseMultSelection){
   AliMultSelection *SelectMult = (AliMultSelection*) event->FindListObject("MultSelection");
   //AliMultSelection *SelectMult = dynamic_cast<AliMultSelection*>( event->FindListObject("MultSelection"));
   if (!SelectMult) AliFatal("MultSelection not found in input event");

      gCentrality = SelectMult->GetMultiplicityPercentile(fCentralityEstimator, kTRUE);
  
       if (gCentrality > fCentralityPercentileMin && gCentrality <fCentralityPercentileMax){
      fHistCentStats->Fill(0.,SelectMult->GetMultiplicityPercentile("V0M", kTRUE));
      fHistCentStats->Fill(1.,SelectMult->GetMultiplicityPercentile("V0A", kTRUE));
      fHistCentStats->Fill(2.,SelectMult->GetMultiplicityPercentile("V0C", kTRUE));
      fHistCentStats->Fill(3.,SelectMult->GetMultiplicityPercentile("FMD", kTRUE));
      fHistCentStats->Fill(4.,SelectMult->GetMultiplicityPercentile("TRK", kTRUE));
      fHistCentStats->Fill(5.,SelectMult->GetMultiplicityPercentile("TKL", kTRUE));
      fHistCentStats->Fill(6.,SelectMult->GetMultiplicityPercentile("CL0", kTRUE));
      fHistCentStats->Fill(7.,SelectMult->GetMultiplicityPercentile("CL1", kTRUE));
      fHistCentStats->Fill(8.,SelectMult->GetMultiplicityPercentile("ZNA", kTRUE));
      fHistCentStats->Fill(9.,SelectMult->GetMultiplicityPercentile("ZPA", kTRUE));
      fHistCentStats->Fill(10.,SelectMult->GetMultiplicityPercentile("V0MvsFMD", kTRUE));
      fHistCentStats->Fill(11.,SelectMult->GetMultiplicityPercentile("TKLvsV0M", kTRUE));
      fHistCentStats->Fill(12.,SelectMult->GetMultiplicityPercentile("ZEMvsZDC", kTRUE));

      // Centrality estimator USED   ++++++++++++++++++++++++++++++
      fHistCentStatsUsed->Fill(0.,SelectMult->GetMultiplicityPercentile(fCentralityEstimator, kTRUE));
      }

       if (gCentrality > 100) gCentrality = -1;

//cout<< "inside LHC150 data "<<endl;

   }

   else {


    AliAODHeader *header = (AliAODHeader*) event->GetHeader();
    if(header){
      gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());

      // QA for centrality estimators
      fHistCentStats->Fill(0.,header->GetCentralityP()->GetCentralityPercentile("V0M"));
      fHistCentStats->Fill(1.,header->GetCentralityP()->GetCentralityPercentile("V0A"));
      fHistCentStats->Fill(2.,header->GetCentralityP()->GetCentralityPercentile("V0C"));
      fHistCentStats->Fill(3.,header->GetCentralityP()->GetCentralityPercentile("FMD"));
      fHistCentStats->Fill(4.,header->GetCentralityP()->GetCentralityPercentile("TRK"));
      fHistCentStats->Fill(5.,header->GetCentralityP()->GetCentralityPercentile("TKL")); 
      fHistCentStats->Fill(6.,header->GetCentralityP()->GetCentralityPercentile("CL0"));
      fHistCentStats->Fill(7.,header->GetCentralityP()->GetCentralityPercentile("CL1"));
      fHistCentStats->Fill(8.,header->GetCentralityP()->GetCentralityPercentile("ZNA"));
      fHistCentStats->Fill(9.,header->GetCentralityP()->GetCentralityPercentile("ZPA"));
      fHistCentStats->Fill(10.,header->GetCentralityP()->GetCentralityPercentile("V0MvsFMD"));
      fHistCentStats->Fill(11.,header->GetCentralityP()->GetCentralityPercentile("TKLvsV0M"));
      fHistCentStats->Fill(12.,header->GetCentralityP()->GetCentralityPercentile("ZEMvsZDC"));
      
      // Centrality estimator USED   ++++++++++++++++++++++++++++++
      fHistCentStatsUsed->Fill(0.,header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data()));
      
      // centrality QA (V0M)
      fHistV0M->Fill(event->GetVZEROData()->GetMTotV0A(), event->GetVZEROData()->GetMTotV0C());
      
      // centrality QA (reference tracks)
      fHistRefTracks->Fill(0.,header->GetRefMultiplicity());
      fHistRefTracks->Fill(1.,header->GetRefMultiplicityPos());
      fHistRefTracks->Fill(2.,header->GetRefMultiplicityNeg());
      fHistRefTracks->Fill(3.,header->GetTPConlyRefMultiplicity());
      fHistRefTracks->Fill(4.,header->GetNumberOfITSClusters(0));
      fHistRefTracks->Fill(5.,header->GetNumberOfITSClusters(1));
      fHistRefTracks->Fill(6.,header->GetNumberOfITSClusters(2));
      fHistRefTracks->Fill(7.,header->GetNumberOfITSClusters(3));
      fHistRefTracks->Fill(8.,header->GetNumberOfITSClusters(4));

    }//AOD header
   } // else 
  }//AOD

  // decide what should be returned only here
  Double_t lReturnVal = -100;
  if(fEventClass=="Centrality"){
    lReturnVal = gCentrality;
  }
  return lReturnVal;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::GetEventPlane(AliVEvent *event){
  // Get the event plane

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  Float_t gVZEROEventPlane    = -10.;
  Float_t gReactionPlane      = -10.;
  Double_t qxTot = 0.0, qyTot = 0.0;

  //MC: from reaction plane
  if(gAnalysisLevel == "MC"){
    if(!event) {
      AliError("mcEvent not available");
      return 0x0;
    }

    AliMCEvent *gMCEvent = dynamic_cast<AliMCEvent*>(event);
    if(gMCEvent){
      AliCollisionGeometry* headerH = dynamic_cast<AliCollisionGeometry*>(gMCEvent->GenEventHeader());    
      if (headerH) {
	gReactionPlane = headerH->ReactionPlaneAngle();
	//gReactionPlane *= TMath::RadToDeg();
      }//MC header
    }//MC event cast
  }//MC
  
  // AOD,ESD,ESDMC: from VZERO Event Plane
  else{
   
    AliEventplane *ep = event->GetEventplane();
    if(ep){ 
      gVZEROEventPlane = ep->CalculateVZEROEventPlane(event,10,2,qxTot,qyTot);
      if(gVZEROEventPlane < 0.) gVZEROEventPlane += TMath::Pi();
      //gReactionPlane = gVZEROEventPlane*TMath::RadToDeg();
      gReactionPlane = gVZEROEventPlane;
    }
  }//AOD,ESD,ESDMC

  return gReactionPlane;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::GetTrackbyTrackCorrectionMatrix( Double_t vEta, 
								Double_t vPhi, 
								Double_t vPt, 
								Short_t vCharge, 
								Double_t gCentrality) {
  // -- Get efficiency correction of particle dependent on (eta, phi, pt, charge, centrality) 

  Double_t correction = 1.;
  Int_t gCentralityInt = -1;

  for (Int_t i=0; i<fCentralityArrayBinsForCorrections-1; i++){
    if((fCentralityArrayForCorrections[i] <= gCentrality)&&(gCentrality <= fCentralityArrayForCorrections[i+1])){
      gCentralityInt = i;
      break;
    }
  }  

  // centrality not in array --> no correction
  if(gCentralityInt < 0){
    correction = 1.;
  }
  else{
    
    //Printf("//=============CENTRALITY=============// %d:",gCentralityInt);
    
    if(fHistCorrectionPlus[gCentralityInt]){
      if (vCharge > 0) {
	correction = fHistCorrectionPlus[gCentralityInt]->GetBinContent(fHistCorrectionPlus[gCentralityInt]->FindBin(vEta,vPt,vPhi));
	//Printf("CORRECTIONplus: %.2f | Centrality %d",correction,gCentralityInt);  
      }
      if (vCharge < 0) {
	correction = fHistCorrectionMinus[gCentralityInt]->GetBinContent(fHistCorrectionMinus[gCentralityInt]->FindBin(vEta,vPt,vPhi));
	//Printf("CORRECTIONminus: %.2f | Centrality %d",correction,gCentralityInt); 
      }
    }
    else {
      correction = 1.;
    }
  }//centrality in array
  
  if (correction == 0.) { 
  //  AliError(Form("Should not happen : bin content = 0. >> eta: %.2f | phi : %.2f | pt : %.2f | cent %d",vEta, vPhi, vPt, gCentralityInt)); 
    correction = 1.; 
  } 
  
  return correction;
}

//________________________________________________________________________
TObjArray* AliAnalysisTaskPIDBF::GetAcceptedTracks(AliVEvent *event, Double_t gCentrality, Double_t gReactionPlane){
  // Returns TObjArray with tracks after all track cuts (only for AOD!)
  // Fills QA histograms

  TString gAnalysisLevel = fBalance->GetAnalysisLevel();

  //output TObjArray holding all good tracks
  TObjArray* tracksAccepted = new TObjArray;
  tracksAccepted->SetOwner(kTRUE);

  Short_t vCharge;
  Double_t vEta;
  Double_t vPhi;
  Double_t vPt;
  Double_t vY=-999;

   
  Double_t MassPID = 0.0;
  Double_t MassPion   = 0.139570; // GeV/c2
  Double_t MassKaon   = 0.493677; // GeV/c2
  Double_t MassProton = 0.938272; // GeV/c2

   if(fUsePID && fRapidityInsteadOfEta){
   if ( fParticleType_ == kPion_ )  MassPID = MassPion;
   else if( fParticleType_ == kKaon_ )  MassPID = MassKaon;
   else if( fParticleType_ == kProton_ )  MassPID = MassProton;
  }

  if(gAnalysisLevel == "AOD") { // handling of TPC only tracks different in AOD and ESD
    // Loop over tracks in event
   

// after PAG meeting on 31.08.2017 

  std::map<int, int> labels;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i = 0; i < event->GetNumberOfTracks(); i++) {
    const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(event->GetTrack(i));
    if (!aodtrack->TestFilterBit(fnAODtrackCutBit)) {
      // Skip TPC-only tracks
      if (aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }

// after PAG meeting on 31.08.2017 

 
    for (Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++) {
      AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(event->GetTrack(iTracks));
      if (!aodTrack) {
	AliError(Form("Could not receive track %d", iTracks));
	continue;
      }
      
      // AOD track cuts
      
      // For ESD Filter Information: ANALYSIS/macros/AddTaskESDfilter.C
      // take only TPC only tracks 
      //fHistTrackStats->Fill(aodTrack->GetFilterMap());
      for(Int_t iTrackBit = 0; iTrackBit < 16; iTrackBit++){
	fHistTrackStats->Fill(iTrackBit,aodTrack->TestFilterBit(1<<iTrackBit));
      }

      if(!aodTrack->TestFilterBit(fnAODtrackCutBit)) continue;


      // additional check on kPrimary flag
      if(fCheckPrimaryFlagAOD){
	if(aodTrack->GetType() != AliAODTrack::kPrimary)
	  continue;
      }

     
      vCharge = aodTrack->Charge();
      vEta    = aodTrack->Eta();
      vPhi    = aodTrack->Phi();// * TMath::RadToDeg();
      vPt     = aodTrack->Pt();
 

     // dcaXY = aodTrack->DCA();      // this is the DCA from global track (not exactly what is cut on)
      //DCAZ  = aodTrack->ZAtDCA();   // this is the DCA from global track (not exactly what is cut on)


    if( vPt < fPtMin || vPt > fPtMax)      continue;

    if( vEta < fEtaMin || vEta > fEtaMax)  continue;
    if (aodTrack->P() == 0.0) continue;
  

    Double_t b[2] = {-99., -99.};
    Double_t bCov[3] = {-99., -99., -99.};

    AliAODTrack* trackAODClone = new AliAODTrack(*aodTrack);
    if (!trackAODClone) {
      AliWarning("Clone of AOD track failed.");
      delete trackAODClone;
      continue;
    }   
    if (!trackAODClone->PropagateToDCA(event->GetPrimaryVertex(), event->GetMagneticField(), 100., b, bCov)){
      delete trackAODClone;
      continue;
    } else {
      delete trackAODClone;
    }
    
    Float_t dcaXY   = b[0];
    Float_t DCAZ   = b[1]; 
 
     if( fDCAxyCut != -1 && fDCAzCut != -1){


        if (DCAZ     <  -fDCAzCut || DCAZ   > fDCAzCut ||  dcaXY    > fDCAxyCut ) continue;

//        if(TMath::Sqrt((dcaXY*dcaXY)/(fDCAxyCut*fDCAxyCut)+(DCAZ*DCAZ)/(fDCAzCut*fDCAzCut)) > 1 )
  //        continue;  // 2D cut
       
}

        Int_t nITSclus = aodTrack->GetITSNcls();
       if(nITSclus < 4 || nITSclus > 100)
        continue;


        // Extra TPC cuts (for systematic studies [!= -1])
      if( fTPCchi2Cut != -1 && aodTrack->Chi2perNDF() > fTPCchi2Cut){
        continue;
      }
      if( fNClustersTPCCut != -1 && (aodTrack->GetTPCNcls() < fNClustersTPCCut || aodTrack->GetTPCNcls() >150 )){
        continue;
      }

      // Extra cut on shared clusters
      if( fTPCsharedCut != -1 && aodTrack->GetTPCnclsS() > fTPCsharedCut){
        continue;
      }

      //===========================PID (so far only for electron rejection)===============================//		    
      if(fElectronRejection) {

	// get the electron nsigma
	Double_t nSigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kElectron));
	
	//Fill QA before the PID
	fHistdEdxVsPTPCbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtbeforePIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	//end of QA-before pid
	
	// check only for given momentum range
	if( vPt > fElectronRejectionMinPt && vPt < fElectronRejectionMaxPt ){
	  	  
	  //look only at electron nsigma
	  if(!fElectronOnlyRejection){
	    
	    //Make the decision based on the n-sigma of electrons
	    if(nSigma < fElectronRejectionNSigma) continue;
	  }
	  else{
	    
	    Double_t nSigmaPions   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kPion));
	    Double_t nSigmaKaons   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kKaon));
	    Double_t nSigmaProtons = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(aodTrack,(AliPID::EParticleType)AliPID::kProton));
	    
	    //Make the decision based on the n-sigma of electrons exclusively ( = track not in nsigma region for other species)
	    if(nSigma < fElectronRejectionNSigma
	       && nSigmaPions   > fElectronRejectionNSigma
	       && nSigmaKaons   > fElectronRejectionNSigma
	       && nSigmaProtons > fElectronRejectionNSigma ) continue;
	  }
	}
  
	//Fill QA after the PID
	fHistdEdxVsPTPCafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),aodTrack->GetTPCsignal());
	fHistNSigmaTPCvsPtafterPIDelectron -> Fill(aodTrack->P()*aodTrack->Charge(),nSigma); 
	
      }
      

//  if(fUsePID && fRapidityInsteadOfEta) vY = log( ( sqrt(MassPID*MassPID + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) / sqrt(MassPID*MassPID + vPt*vPt) ); // convert eta to y

//   vPionYReco = 0.5*log( ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) + vPt*sinh(vEta) ) /  ( sqrt(MassPion*MassPion + vPt*vPt*cosh(vEta)*cosh(vEta)) - vPt*sinh(vEta) ) ); // convert eta to y

//cout<< " vEta :" << vEta<<'\t'<<" vY :"<<vY<<endl;


// ON 4.09.2017 :  For TPC Only tracks we have to copy PID information from corresponding global tracks



   AliVTrack *trackAli=(AliVTrack*)event->GetTrack(iTracks);


    const Int_t pid_track_id = (fnAODtrackCutBit == (1 << 7))
                             ? labels[-1 - trackAli->GetID()]
                             : iTracks;
    AliAODTrack *aodtrackpid = dynamic_cast<AliAODTrack *>(event->GetTrack(pid_track_id));

   
  if(fUsePID && fRapidityInsteadOfEta) vY = log( ( sqrt(MassPID*MassPID + aodtrackpid->Pt()*aodtrackpid->Pt()*cosh(aodtrackpid->Eta())*cosh(aodtrackpid->Eta())) + aodtrackpid->Pt()*sinh(aodtrackpid->Eta()) ) / sqrt(MassPID*MassPID + aodtrackpid->Pt()*aodtrackpid->Pt()) ); // convert eta to y


// ON 4.09.2017 :  For TPC Only tracks we have to copy PID information from corresponding global tracks




if(fUsePID){

Double_t betaParticle= Beta(aodtrackpid);
if(betaParticle >1.0) continue;

Float_t MisMatchTOFProb = fPIDResponse->GetTOFMismatchProbability(aodtrackpid);

if(fTOFMisMatch){
if(MisMatchTOFProb > fMistMatchTOFProb) continue;
}

Int_t ParticleSpecies=GetSpecies(aodtrackpid);

fPIDSpeciesHisto->Fill(ParticleSpecies);

if(fParticleType_ == kPion_){
if(ParticleSpecies !=kPion_) continue;
}

else if(fParticleType_ == kKaon_){
if(ParticleSpecies !=kKaon_) continue;
}

else if(fParticleType_ == kProton_){
if(ParticleSpecies !=kProton_) continue;
}

if(fRapidityInsteadOfEta){
if(vY <fEtaMin | vY > fEtaMax ) continue;
fHistRapidity->Fill(vY,gCentrality);
}

IsTOF(aodtrackpid);
if(fHasTOFPID && aodtrackpid->Pt() > fPtTOFMin) fHistTOFPid->Fill(aodtrackpid->P()*aodtrackpid->Charge(), betaParticle);
fHistdEdxPid->Fill(aodtrackpid->P()*aodtrackpid->Charge(),aodtrackpid->GetTPCsignal());

} // end of PID selection 

// ----------------------------------------------------------------------------Implementation of PID 

      vCharge = aodtrackpid->Charge();
      vEta    = aodtrackpid->Eta();
      vPhi    = aodtrackpid->Phi();// * TMath::RadToDeg();
      vPt     = aodtrackpid->Pt();

      // fill QA histograms
      fHistClus->Fill(aodtrackpid->GetITSNcls(),aodtrackpid->GetTPCNcls());
      fHistDCA->Fill(DCAZ,dcaXY);
      fHistChi2->Fill(aodtrackpid->Chi2perNDF(),gCentrality);
      fHistPt->Fill(vPt,gCentrality);
      fHistEta->Fill(vEta,gCentrality);
      fHistEta1D->Fill(vEta);
      fHistPhi1D->Fill(vPhi);
      if(vCharge > 0) fHistPhiPos->Fill(vPhi,gCentrality);
      else if(vCharge < 0) fHistPhiNeg->Fill(vPhi,gCentrality);
      fHistPhi->Fill(vPhi,gCentrality);
      if(fUsePID && fRapidityInsteadOfEta){
      if(vCharge > 0)      fHistEtaPhiPos->Fill(vY,vPhi,gCentrality); 		 
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vY,vPhi,gCentrality);
      }

      else {
       if(vCharge > 0)      fHistEtaPhiPos->Fill(vEta,vPhi,gCentrality);
      else if(vCharge < 0) fHistEtaPhiNeg->Fill(vEta,vPhi,gCentrality);

     }
      
      //=======================================correction
      Double_t correction=0.0;
      if(fUsePID && fRapidityInsteadOfEta) {
      correction = GetTrackbyTrackCorrectionMatrix(vY,vPhi, vPt, vCharge, gCentrality);  
}

      else {
       correction = GetTrackbyTrackCorrectionMatrix(vEta, vPhi, vPt, vCharge, gCentrality);  
      }



      //Printf("CORRECTIONminus: %.2f | Centrality %lf",correction,gCentrality);
      
      // add the track to the TObjArray
    
 
      if(fUsePID && fRapidityInsteadOfEta){
      tracksAccepted->Add(new AliBFBasicParticle(vY, vPhi, vPt, vCharge, correction));  
} 

else {
      tracksAccepted->Add(new AliBFBasicParticle(vEta, vPhi, vPt, vCharge, correction));  

}


    }//track loop
  } //AOD analysis
 
 
  return tracksAccepted;  

}


//________________________________________________________________________
void  AliAnalysisTaskPIDBF::SetVZEROCalibrationFile(const char* filename,
						    const char* lhcPeriod) {
  //Function to setup the VZERO gain equalization
    //============Get the equilization map============//
  TFile *calibrationFile = TFile::Open(filename);
  if((!calibrationFile)||(!calibrationFile->IsOpen())) {
    Printf("No calibration file found!!!");
    return;
  }

  TList *list = dynamic_cast<TList *>(calibrationFile->Get(lhcPeriod));
  if(!list) {
    Printf("Calibration TList not found!!!");
    return;
  }

  fHistVZEROAGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROAGainEqualizationMap"));
  if(!fHistVZEROAGainEqualizationMap) {
    Printf("VZERO-A calibration object not found!!!");
    return;
  }
  fHistVZEROCGainEqualizationMap = dynamic_cast<TH1F *>(list->FindObject("gHistVZEROCGainEqualizationMap"));
  if(!fHistVZEROCGainEqualizationMap) {
    Printf("VZERO-C calibration object not found!!!");
    return;
  }

  fHistVZEROChannelGainEqualizationMap = dynamic_cast<TH2F *>(list->FindObject("gHistVZEROChannelGainEqualizationMap"));
  if(!fHistVZEROChannelGainEqualizationMap) {
    Printf("VZERO channel calibration object not found!!!");
    return;
  }
}

//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::GetChannelEqualizationFactor(Int_t run, 
							    Int_t channel) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  for(Int_t iBinX = 1; iBinX <= fHistVZEROChannelGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROChannelGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    if(gRunNumber == run)
      return fHistVZEROChannelGainEqualizationMap->GetBinContent(iBinX,channel+1);
  }

  return 1.0;
}

//________________________________________________________________________
Double_t AliAnalysisTaskPIDBF::GetEqualizationFactor(Int_t run, 
						     const char* side) {
  //
  if(!fHistVZEROAGainEqualizationMap) return 1.0;

  TString gVZEROSide = side;
  for(Int_t iBinX = 1; iBinX < fHistVZEROAGainEqualizationMap->GetNbinsX(); iBinX++) {
    Int_t gRunNumber = atoi(fHistVZEROAGainEqualizationMap->GetXaxis()->GetBinLabel(iBinX));
    //cout<<"Looking for run "<<run<<" - current run: "<<gRunNumber<<endl;
    if(gRunNumber == run) {
      if(gVZEROSide == "A") 
	return fHistVZEROAGainEqualizationMap->GetBinContent(iBinX);
      else if(gVZEROSide == "C") 
	return fHistVZEROCGainEqualizationMap->GetBinContent(iBinX);
    }
  }

  return 1.0;
}


//____________________________________________________________________
Bool_t AliAnalysisTaskPIDBF::IsThisAWeakDecayingParticle(TParticle *thisGuy)
{
  // In order to prevent analyzing daughters from weak decays 
  // - AMPT does not only strong decays, so IsPhysicalPrimary does not catch it

 Int_t pdgcode = TMath::Abs( thisGuy->GetPdgCode() );

 Int_t myWeakParticles[7] = { 3322, 3312, 3222, // Xi0 Xi+- Sigma-+
			       3122, 3112, // Lambda0 Sigma+-
			       130, 310 // K_L0 K_S0
 };

 Bool_t found = kFALSE;
 for(Int_t i=0; i!=7; ++i)
   if( myWeakParticles[i] == pdgcode ) {
     found = kTRUE;
     break;
   }
 
 return found;
}

//________________________________________________________________________
void  AliAnalysisTaskPIDBF::FinishTaskOutput(){
  //Printf("END BF");

  if (!fBalance) {
    AliError("fBalance not available");
    return;
  }  

}


void AliAnalysisTaskPIDBF::IsTOF(AliVTrack *track) 
 {

 if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF,track)==0) fHasTOFPID=kFALSE;
  else fHasTOFPID=kTRUE;

  //in addition to TOF status we look at the pt
  if(track->Pt()<fPtTOFMin)fHasTOFPID=kFALSE;

      Int_t startTimeMask = fPIDResponse->GetTOFResponse().GetStartTimeMask(track->P());
      if (startTimeMask < 0)fHasTOFPID=kFALSE;


 } //End of IsTOF


Double_t AliAnalysisTaskPIDBF::Beta(AliVTrack *track)
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



void AliAnalysisTaskPIDBF::GetNsigmas(AliVTrack* trk)
{

  AliVParticle *inEvHMain = dynamic_cast<AliVParticle *>(trk);
  // --- TPC
  Double_t nsigmaTPCkProton = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kProton);
  Double_t nsigmaTPCkKaon   = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kKaon);
  Double_t nsigmaTPCkPion   = fPIDResponse->NumberOfSigmasTPC(inEvHMain, AliPID::kPion);
  // --- TOF
  Double_t nsigmaTOFkProton=999.,nsigmaTOFkKaon=999.,nsigmaTOFkPion=999.;
  Double_t nsigmaTPCTOFkProton=999.,nsigmaTPCTOFkKaon=999.,nsigmaTPCTOFkPion=999.;


   IsTOF(trk);

  if(fHasTOFPID && trk->Pt()>fPtTOFMin){//use TOF information
    nsigmaTOFkProton = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kProton);
    nsigmaTOFkKaon   = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kKaon);
    nsigmaTOFkPion   = fPIDResponse->NumberOfSigmasTOF(inEvHMain, AliPID::kPion);
    Double_t d2Proton=nsigmaTPCkProton * nsigmaTPCkProton + nsigmaTOFkProton * nsigmaTOFkProton;
    Double_t d2Kaon=nsigmaTPCkKaon * nsigmaTPCkKaon + nsigmaTOFkKaon * nsigmaTOFkKaon;
    Double_t d2Pion=nsigmaTPCkPion * nsigmaTPCkPion + nsigmaTOFkPion * nsigmaTOFkPion;

    nsigmaTPCTOFkProton  =  TMath::Sqrt(d2Proton);
    nsigmaTPCTOFkKaon    =  TMath::Sqrt(d2Kaon);
    nsigmaTPCTOFkPion    =  TMath::Sqrt(d2Pion);



  }
   else{
    // --- combined
    // if TOF is missing and below fPtTOFPID only the TPC information is used
    nsigmaTPCTOFkProton = TMath::Abs(nsigmaTPCkProton);
    nsigmaTPCTOFkKaon   = TMath::Abs(nsigmaTPCkKaon);
    nsigmaTPCTOFkPion   = TMath::Abs(nsigmaTPCkPion);
  }

  fnsigmas[kPion_][kTPC_]=nsigmaTPCkPion;
  fnsigmas[kKaon_][kTPC_]=nsigmaTPCkKaon;
  fnsigmas[kProton_][kTPC_]=nsigmaTPCkProton;
  fnsigmas[kPion_][kTOF_]=nsigmaTOFkPion;
  fnsigmas[kKaon_][kTOF_]=nsigmaTOFkKaon;
  fnsigmas[kProton_][kTOF_]=nsigmaTOFkProton;
  fnsigmas[kPion_][kTPCTOFpid_]=nsigmaTPCTOFkPion;
  fnsigmas[kKaon_][kTPCTOFpid_]=nsigmaTPCTOFkKaon;
  fnsigmas[kProton_][kTPCTOFpid_]=nsigmaTPCTOFkProton;

  for(Int_t ipart=0;ipart<kProton_+1;ipart++){
      for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
        if((ipid!=kTPC_) && (!fHasTOFPID) && trk->Pt()<fPtTOFMin)continue;//not filling TOF and combined if no TOF PID
        TH2F *h=GetHistogram2D(Form("NSigma_%d_%d",ipart,ipid));
        h->Fill(trk->Pt(),fnsigmas[ipart][ipid]);
      }
    }



}

// Back to Histogram for Filling : 

TH2F* AliAnalysisTaskPIDBF::GetHistogram2D(const char * name){
  // returns histo named name
  return (TH2F*) fHistListPIDQA->FindObject(name);
}



// Minimum Nsigma : 

 Int_t AliAnalysisTaskPIDBF::MinNsigma(AliVTrack * trk, Bool_t FillHistos)
{
   // Return index of array with the minimum element.
   // If more than one element is minimum returns first found.

  IsTOF(trk);
  if( (!fHasTOFPID) && trk->Pt()>fPtTOFMin)return kSpUndefined;

  //get the identity of the particle with the minimum Nsigma
  Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fDetectorPID_){
  case kTPC_:
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTPC_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTPC_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTPC_])  ;
    break;
  case kTOF_:
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTOF_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTOF_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTOF_])  ;
    break;
  case kTPCTOFpid_://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTPCTOFpid_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTPCTOFpid_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTPCTOFpid_])  ;
    break;
  }


   if( ( nsigmaKaon==nsigmaPion ) && ( nsigmaKaon==nsigmaProton )) return kSpUndefined;

// Pion 

  if( ( nsigmaPion   < nsigmaKaon ) && ( nsigmaPion < nsigmaProton ) && (nsigmaPion   < fPIDNSigma)){

  if(FillHistos){
   for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
        if((ipid!=kTPC_) && (!fHasTOFPID) && (trk->Pt()<fPtTOFMin))continue;//not filling TOF and combined if no TOF PID
        TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kPion_,ipid));
        h->Fill(trk->Pt(),fnsigmas[kPion_][ipid]);
      }
}    

    return kPion_;
  }

 // Kaon

   if( ( nsigmaKaon   < nsigmaPion ) && ( nsigmaKaon < nsigmaProton ) && (nsigmaKaon   < fPIDNSigma)){

    if(FillHistos){

      for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
        if((ipid!=kTPC_) && (!fHasTOFPID) && (trk->Pt()<fPtTOFMin))continue;//not filling TOF and combined if no TOF PID
        TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kKaon_,ipid));
        h->Fill(trk->Pt(),fnsigmas[kKaon_][ipid]);
      }
}    

    return kKaon_;
  }


// Proton

   if( ( nsigmaProton < nsigmaKaon ) && ( nsigmaProton < nsigmaPion ) && (nsigmaProton < fPIDNSigma)){

 if(FillHistos){

   for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
        if((ipid!=kTPC_) && (!fHasTOFPID) && (trk->Pt()<fPtTOFMin))continue;//not filling TOF and combined if no TOF PID
        TH2F *h=GetHistogram2D(Form("NSigmaRec_%d_%d",kProton_,ipid));
        h->Fill(trk->Pt(),fnsigmas[kProton_][ipid]);
      }
}

    return kProton_;
  }

  return kSpUndefined;

}

Bool_t* AliAnalysisTaskPIDBF::GetDoubleCounting(AliVTrack * trk){

for(Int_t ipart=0;ipart<kProton_+1;ipart++)fHasDoubleCounting[ipart]=kFALSE;
 Int_t SigmaMin=MinNsigma(trk,kFALSE);

 IsTOF(trk);

 if(SigmaMin==kSpUndefined)return fHasDoubleCounting;

 Double_t nsigmaPion=999., nsigmaKaon=999., nsigmaProton=999.;
  switch (fDetectorPID_){
  case kTPC_:
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTPC_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTPC_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTPC_])  ;
    break;
  case kTOF_:
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTOF_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTOF_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTOF_])  ;
    break;
  case kTPCTOFpid_://In case of no TOF matching the combined nsigma is the TPC one
    nsigmaProton  =  TMath::Abs(fnsigmas[kProton_][kTPCTOFpid_]);
    nsigmaKaon    =  TMath::Abs(fnsigmas[kKaon_][kTPCTOFpid_])  ;
    nsigmaPion    =  TMath::Abs(fnsigmas[kPion_][kTPCTOFpid_])  ;
    break;
  }

  if(nsigmaPion<fPIDNSigma && SigmaMin!=kPion_)fHasDoubleCounting[kPion_]=kTRUE;
  if(nsigmaKaon<fPIDNSigma && SigmaMin!=kKaon_)fHasDoubleCounting[kKaon_]=kTRUE;
  if(nsigmaProton<fPIDNSigma && SigmaMin!=kProton_)fHasDoubleCounting[kProton_]=kTRUE;

  for(Int_t ipart=0;ipart<kProton_+1;ipart++){
      if(fHasDoubleCounting[ipart]){
        for(Int_t ipid=0;ipid<kTPCTOFpid_+1;ipid++){
          if((ipid!=kTPC_) && (!fHasTOFPID) && (trk->Pt()<fPtTOFMin))continue;//not filling TOF and combined if no TOF PID
          TH2F *h=GetHistogram2D(Form("NSigmaDC_%d_%d",ipart,ipid));
          h->Fill(trk->Pt(),fnsigmas[ipart][ipid]);
        }
      }
    }

return fHasDoubleCounting;

} // End of Double Counting


Int_t AliAnalysisTaskPIDBF:: GetSpecies(AliVTrack *trk ){

Int_t ID=kSpUndefined;
GetNsigmas(trk);

ID=MinNsigma(trk,kTRUE);

Bool_t *HasDC;
HasDC=GetDoubleCounting(trk);
for(Int_t ipart=0;ipart<kProton_+1;ipart++){
if(HasDC[ipart]==kTRUE) ID = kSpUndefined;
}

if(ID != kSpUndefined){
      for(Int_t idet=0;idet<kNDetectors_D;idet++){
        TH2F *h=GetHistogram2D(Form("PID_%d_%d",idet,ID));
        if(idet==kITS_D)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
        if(idet==kTPC_D)h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
        if(idet==kTOF_D && fHasTOFPID)h->Fill(trk->P(),Beta(trk)*trk->Charge());
      }
    }

 for(Int_t idet=0;idet<kNDetectors_D;idet++){
      TH2F *h=GetHistogram2D(Form("PIDAll_%d",idet));
      if(idet==kITS_D)h->Fill(trk->P(),trk->GetITSsignal()*trk->Charge());
      if(idet==kTPC_D) h->Fill(trk->P(),trk->GetTPCsignal()*trk->Charge());
      if(idet==kTOF_D && fHasTOFPID)h->Fill(trk->P(),Beta(trk)*trk->Charge());
    }

//cout<< " particle id :"<<ID<<endl;
return ID;
} // End of the code 


//________________________________________________________________________
void AliAnalysisTaskPIDBF::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}

