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

#include "AliFlowBayesianPID.h"

#include "AliAnalysisUtils.h"
#include "AliESDpid.h"

#include "AliAnalysisTaskNUA.h"

#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskNUA)

//________________________________________________________________________
AliAnalysisTaskNUA::AliAnalysisTaskNUA(const char *name) 
: AliAnalysisTaskSE(name), 
  gAOD(0),
  fListQA(0), fListResults(0),
  fHistEventStats(0),
  fHistTrackStats(0),
  fHistVx(0), fHistVy(0), fHistVz(0),
  fCentralityEstimator("V0M"),
  fCentralityPercentileMin(0.), 
  fCentralityPercentileMax(5.),
  fHistMultiplicity(0),
  fHistMultiplicity_counts(0),
  fHistPt(0),
  fHistCent(0),
  fHistCentbin(0),
  fHistPtbin(0),
  fHistPhi(0),
  fHistEta(0),
  fHistPtCen(0),
  fHistEtaCen(0),
  fHistPhiCen(0),
  fHistNClustersTPC(0),
  fHistChi2PerClusterTPC(0),
  fHistDCAToVertexZ(0),
  fHistDCAToVertexXY(0),
  fHistDCAToVertex2D(0),
  fHistEtaPhiCent(0),
  fHistPtEtaCent(0),
  fHistPtPhiCent(0),
  fUseOfflineTrigger(kFALSE),
  fVxMax(0.3), fVyMax(0.3), fVzMax(10.),
  fAODtrackCutBit(128),
  fEtaMin(-0.8), fEtaMax(0.8),
  fUtils(0),
  fHistEtaPhiVertexPlus(0),
  fHistEtaPhiVertexMinus(0){
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskNUA::~AliAnalysisTaskNUA() {
  //
}

//________________________________________________________________________
void AliAnalysisTaskNUA::UserCreateOutputObjects() {
  // Create histograms
  // Called once
  Int_t phiBin = 100;
  Int_t etaBin = 16;
  Int_t vertex_bin = 4;

  Double_t nArrayPhi[phiBin+1];
    for(Int_t iBin = 0; iBin <= phiBin; iBin++)
        nArrayPhi[iBin] = iBin*TMath::TwoPi()/phiBin;

  Double_t nArrayEta[17]={-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
  Double_t nArrayVertex[5]={-10,-5,0,5,10};

  fUtils = new AliAnalysisUtils();

  //QA list
  fListQA = new TList();
  fListQA->SetName("listQA");
  fListQA->SetOwner();

  //Event stats.
  TString gCutName[5] = {"Total","Offline trigger",
                         "Vertex","sel. Centrality",
			 "No Pile-Up"};
  fHistEventStats = new TH2F("fHistEventStats",
                             "Event statistics;;Centrality percentile;N_{events}",
                             5,0.5,5.5,110,-5,105);
  for(Int_t i = 1; i <= 5; i++)
  fHistEventStats->GetXaxis()->SetBinLabel(i,gCutName[i-1].Data());
  fListQA->Add(fHistEventStats);

  fHistTrackStats = new TH2F("fHistTrackStats","Event statistics;Centrality (%);FilterBit",110,-5,105,1300,0,1300);
  fListQA->Add(fHistTrackStats);

  //LHC11h: Change all the centralioty binning from 10,0,100 to 110,-5,105
  // Vertex distributions
  fHistVx = new TH2F("fHistVx","Primary vertex distribution - x coordinate;Centrality (%);V_{x} (cm);Entries",10,0,100,100,-0.5,0.5);
  fListQA->Add(fHistVx);
  fHistVy = new TH2F("fHistVy","Primary vertex distribution - y coordinate;Centrality (%);V_{y} (cm);Entries",10,0,100,100,-0.5,0.5);
  fListQA->Add(fHistVy);
  fHistVz = new TH2F("fHistVz","Primary vertex distribution - z coordinate;Centrality (%);V_{z} (cm);Entries",10,0,100,100,-20.,20.);
  fListQA->Add(fHistVz);

  // QA histograms: multiplicities
  fHistMultiplicity = new TH2F("fHistMultiplicity",";Centrality (%);N_{acc.};Counts",10,0,100,5000,-0.5,4999.5);
  fHistMultiplicity_counts = new TH1F("fHistMultiplicity_counts",";Multiplicity;Counts",100,0,10000);
  fListQA->Add(fHistMultiplicity);
  fListQA->Add(fHistMultiplicity_counts);
  
  //====================================================//
  //Results TList
  fListResults = new TList();
  fListResults->SetName("listResults");
  fListResults->SetOwner();

  //Result: pT spectra
  fHistPt = new TH1F("fHistPt","Pt distribution;p_{T} (GeV/c);Counts",100,0,10);
  fListResults->Add(fHistPt);
  fHistCent = new TH1F("fHistCent","Centrality distribution;Centrality;Counts",10,0,100);
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
  fHistPhiCen = new TH2F("fHistPhiCen","Phi vs Centrality;Centrality;Phi;Counts",10,0,100,100,-3.4,3.4);
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
                                          "Survived positive particles;#phi;#eta;V_{z} (cm)",
                                          phiBin, nArrayPhi, etaBin, nArrayEta, vertex_bin, nArrayVertex);

  fHistEtaPhiVertexMinus = new TH3D("fHistEtaPhiVertexMinus",
                                          "Survived negative particles;#phi;#eta;V_{z} (cm)",
                                          phiBin, nArrayPhi, etaBin,nArrayEta, vertex_bin, nArrayVertex);

  fListResults->Add(fHistEtaPhiVertexPlus);
  fListResults->Add(fHistEtaPhiVertexMinus);


  // Post output data
  PostData(1, fListQA);
  PostData(2, fListResults);
}

//________________________________________________________________________
void AliAnalysisTaskNUA::UserExec(Option_t *) {
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

  Int_t nAcceptedTracks = 0;
  Float_t gCentrality = -1;

  //Centrality
  if(header)
  gCentrality = header->GetCentralityP()->GetCentralityPercentile(fCentralityEstimator.Data());

  /*  
  AliMultSelection *multSelection = 0x0; 
  multSelection = (AliMultSelection*) gAOD->FindListObject("MultSelection");
  if (!multSelection){
    AliWarning("AliMultSelection object not found!");
  }
  else{
    gCentrality = multSelection->GetMultiplicityPercentile(fCentralityEstimator.Data());
  }

  */ 

  // event selection done in AliAnalysisTaskSE::Exec() --> this is not used

  fHistEventStats->Fill(1,gCentrality); //all events
  Bool_t isSelected = kTRUE;
  if(fUseOfflineTrigger)
    isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if(isSelected) {
    fHistEventStats->Fill(2,gCentrality); //triggered events
    
    const AliAODVertex *vertex = gAOD->GetPrimaryVertex();
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
		
		//if((gCentrality >= fCentralityPercentileMin) && 
		//(gCentrality < fCentralityPercentileMax)) {
		if(gCentrality >= 0) { 
		  fHistEventStats->Fill(4,gCentrality); 
		  
		  // check for pile-up event
		  //if(fCheckPileUp) {
		  //fUtils->SetUseMVPlpSelection(kTRUE);
		  //fUtils->SetUseOutOfBunchPileUp(kTRUE);
		  //if(!fUtils->IsPileUpEvent(gAOD)) {
		  fHistEventStats->Fill(5,gCentrality); 
		  
		 // Printf("There are %d tracks in this event", gAOD->GetNumberOfTracks());
		  for (Int_t iTracks = 0; iTracks < gAOD->GetNumberOfTracks(); iTracks++) {
		    AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(gAOD->GetTrack(iTracks));
		    if (!aodTrack) {
		      Printf("ERROR: Could not receive track %d", iTracks);
		      continue;
		    }
		    
		    // AOD track cuts
		    fHistTrackStats->Fill(gCentrality,aodTrack->GetFilterMap());
	        
            if(!aodTrack->TestFilterBit(fAODtrackCutBit)) continue;
            	    	    
 	    Float_t pt  = aodTrack->Pt();
            Float_t eta = aodTrack->Eta();
            Float_t phi = aodTrack->Phi();
            //Int_t numberofclustersTPC = aodTrack->GetNumberOfTPCClusters();
            Double_t chi2 = aodTrack->Chi2perNDF();
            Double_t xdca = aodTrack->DCA();
            Double_t zdca = aodTrack->ZAtDCA();
		    
		    // Kinematics cuts from ESD track cuts
		    if( eta < fEtaMin || eta > fEtaMax) continue;
		     
            fHistPt->Fill(pt);
            fHistPtbin->Fill(pt);
            fHistCent->Fill(gCentrality);
            fHistCentbin->Fill(gCentrality);
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
	
	    fHistEtaPhiVertexPlus->Fill(phi,eta,vertex->GetZ());
	    fHistEtaPhiVertexMinus->Fill(phi,eta,vertex->GetZ());

              
		    nAcceptedTracks += 1;
		  } //track loop
		  fHistMultiplicity->Fill(gCentrality,nAcceptedTracks);
		if((gCentrality > 20.)&&(gCentrality <= 30.)) 
		fHistMultiplicity_counts->Fill(nAcceptedTracks);
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
void  AliAnalysisTaskNUA::FinishTaskOutput(){
  //
}

//________________________________________________________________________
void AliAnalysisTaskNUA::Terminate(Option_t *) {
  // Draw result to the screen
  // Called once at the end of the query

  // not implemented ...

}

