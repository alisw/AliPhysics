#include <TChain.h>
#include <TList.h>

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

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

#include "AliAnalysisTaskMinijet.h"

// Analysis task for two-particle correlations using all particles over pt threshold
// pt_trig threshold for trigger particle (event axis) and pt_assoc for possible associated particles.
// Extract mini-jet yield and fragmentation properties via Delta-Phi histograms of these correlations
// post processing of analysis output via macro plot3and2Gaus.C
// Can use ESD or AOD, reconstructed and Monte Carlo data as input
// Author: Eva Sicking


ClassImp(AliAnalysisTaskMinijet)

//________________________________________________________________________
  AliAnalysisTaskMinijet::AliAnalysisTaskMinijet(const char *name)
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
      fVertexZCut(10.),
      fEtaCut(0.9),
      fEtaCutSeed(0.9),
      fSelectParticles(1),
      fSelectParticlesAssoc(1),
      fESDEvent(0),
      fAODEvent(0),
      fNMcPrimAccept(0),
      fNRecAccept(0),
      fVzEvent(0),
      fMeanPtRec(0),
      fLeadingPtRec(0),
      fHists(0),
      fStep(0),
      fHistPt(0),
      fHistPtMC(0),
      fNmcNch(0),
      fPNmcNch(0),
      fChargedPi0(0)
{

  //Constructor

  for(Int_t i = 0;i< 6;i++){
    fMapSingleTrig[i]         =  0;
    fMapPair[i]               =  0;
    fMapEvent[i]              =  0;
    fMapAll[i]                =  0;

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
  }
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMinijet::~AliAnalysisTaskMinijet()
{
  // Destructor

  if (fHists && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) delete fHists;
}

//________________________________________________________________________
void AliAnalysisTaskMinijet::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  if(fDebug) Printf("In User Create Output Objects.");
   
  fStep = new TH1F("fStep", "fStep", 10, -0.5, 9.5);
  fHistPt = new TH1F("fHistPt", "P_{T} distribution REC", 150, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  if (fUseMC) {
    fHistPtMC = new TH1F("fHistPtMC", "P_{T} distribution MC", 150, 0.1, 3.1);
    fHistPtMC->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtMC->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtMC->SetMarkerStyle(kFullCircle);
    
    fNmcNch = new TH2F("fNmcNch", "fNmcNch", 100,-0.5,99.5,100,-0.5,99.5);
    fPNmcNch = new TProfile("pNmcNch", "pNmcNch", 100,-0.5,99.5);

  }

  fChargedPi0 = new TH2F("fChargedPi0", "fChargedPi0", 200, -0.5, 199.5, 200, -0.5, 199.5);
    

  //----------------------
  //bins for pt in THnSpare
  Double_t ptMin = 0.0, ptMax = 100.;
  Int_t nPtBins = 39; 
  Double_t binsPt[]  = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 
    			0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 
    			10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
  
  //  Int_t nPtBins = 100;
  //   Double_t ptMin = 1.e-2, ptMax = 100.;
  //   Double_t *binsPt = 0;
  //   binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
  
  
  Double_t ptMin2 = 0.0, ptMax2 = 100.;
  Int_t nPtBins2 = 9; 
  Double_t binsPt2[]  = {0.1, 0.4, 0.7, 1.0, 2.0, 3.0, 4.0, 5.0, 10.0, 100.0};
  
  //  Int_t nPtBins2 = 10;
  //   Double_t ptMin2 = 0.4, ptMax2 = 100.;
  //   Double_t *binsPt2 = 0;
  //   binsPt2 = CreateLogAxis(nPtBins2,ptMin2,ptMax2);
  
  //3 dim matrix
  Int_t binsEffHisto[3]   = {Int_t(fEtaCut*20),  nPtBins,    150   };
  Double_t minEffHisto[3] = {-fEtaCut,           ptMin,        -0.5 };
  Double_t maxEffHisto[3] = {fEtaCut,            ptMax,        149.5 };
  
  //5 dim matrix
  Int_t binsEffHisto5[6]   = {  nPtBins2,  nPtBins2,  Int_t(fEtaCut*10),    180,                   150 ,      2 };
  Double_t minEffHisto5[6] = {  ptMin2,    ptMin2,    -2*fEtaCut,          -0.5*TMath::Pi(),      -0.5 ,   -0.5 };
  Double_t maxEffHisto5[6] = {  ptMax2,    ptMax2,     2*fEtaCut,           1.5*TMath::Pi(),     149.5 ,    1.5 };
   

  //4 dim matrix
  Int_t binsEvent[4]   = {   150,        20,   50,  nPtBins };
  Double_t minEvent[4] = {  -0.5,       -10,    0,    ptMin };
  Double_t maxEvent[4] = { 149.5,        10,   10,    ptMax };

  //3 dim matrix
  Int_t binsAll[3]   = {Int_t(fEtaCut*20),  nPtBins2,       150   };
  Double_t minAll[3] = {-fEtaCut,           ptMin2,        -0.5 };
  Double_t maxAll[3] = {fEtaCut,            ptMax2,        149.5 };

  //--------------------
  TString dataType[2] ={"ESD", "AOD"};
  TString labels[6]={Form("%sAllAllMcNmc",dataType[fMode].Data()),
		     Form("%sTrigAllMcNmc",dataType[fMode].Data()),
		     Form("%sTrigAllMcNrec",dataType[fMode].Data()),
		     Form("%sTrigVtxMcNrec",dataType[fMode].Data()),
		     Form("%sTrigVtxRecMcPropNrec",dataType[fMode].Data()),
		     Form("%sTrigVtxRecNrec",dataType[fMode].Data())};


  for(Int_t i=0;i<6;i++){

    fMapSingleTrig[i] = new THnSparseD(Form("fMapSingleTrig%s", labels[i].Data()),"eta:pt:Nrec",
				       3,binsEffHisto,minEffHisto,maxEffHisto);
    fMapSingleTrig[i]->SetBinEdges(1,binsPt);
    fMapSingleTrig[i]->GetAxis(0)->SetTitle("#eta");
    fMapSingleTrig[i]->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    fMapSingleTrig[i]->GetAxis(2)->SetTitle("N_{rec}");
    fMapSingleTrig[i]->Sumw2(); 
    
    fMapPair[i] = new THnSparseD(Form("fMapPair%s", labels[i].Data()),"pt_trig:pt_assoc:DeltaEta:DeltaPhi:Nrec:likesign",
				 6,binsEffHisto5,minEffHisto5,maxEffHisto5);
    fMapPair[i]->SetBinEdges(0,binsPt2);
    fMapPair[i]->SetBinEdges(1,binsPt2);
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
    fMapAll[i]->SetBinEdges(1,binsPt2);
    fMapAll[i]->GetAxis(0)->SetTitle("#eta");
    fMapAll[i]->GetAxis(1)->SetTitle("p_{T} (GeV/c)");
    fMapAll[i]->GetAxis(2)->SetTitle("N_{rec}");
    fMapAll[i]->Sumw2(); 

  
    fVertexZ[i]                  = new TH1F(Form("fVertexZ%s",labels[i].Data()),
					    Form("fVertexZ%s",labels[i].Data()) ,  
					    200, -10., 10.);
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
					    200, -3,3);
    fDcaZ[i]                     = new TH1F(Form("fDcaZ%s",labels[i].Data()),
					    Form("fDcaZ%s",labels[i].Data()) ,  
					    200, -10,10);
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
    
  }

  fHists = new TList();
  fHists->SetOwner();

  fHists->Add(fStep);
  fHists->Add(fHistPt);

  if(fUseMC){
    fHists->Add(fHistPtMC); 
    fHists->Add(fNmcNch); 
    fHists->Add(fPNmcNch); 
  }
  fHists->Add(fChargedPi0);
  
  for(Int_t i=0;i<6;i++){
    fHists->Add(fMapSingleTrig[i]);
    fHists->Add(fMapPair[i]);
    fHists->Add(fMapEvent[i]);
    fHists->Add(fMapAll[i]);
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
  }

  PostData(1, fHists);
 
}

//________________________________________________________________________
void AliAnalysisTaskMinijet::UserExec(Option_t *)
{
  // Main loop, called for each event
  // Kinematics-only, ESD and AOD can be processed.
  // Data is read (LoopESD, LoopAOD...) and then analysed (Analyse). 
  //  - in case of MC with full detector simulation, all correction steps(0-5) can be processed
  //  - for Data, only step 5 is performed
  //  - for kinematics-only, only step 0 is processed
  // step 5 =  Triggered events, reconstructed accepted vertex, reconstructed tracks,                    reconstructed multiplicity, 
  // step 4 =  Triggered events, reconstructed accepted vertex, reconstructed tracks with MC properties, reconstructed multiplicity
  // step 3 =  Triggered events, reconstructed accepted vertex, mc primary particles,                    reconstructed multiplicity, 
  // step 2 =  Triggered events, all                            mc primary particles,                    reconstructed multiplicity
  // step 1 =  Triggered events, all                            mc primary particles,                    true multiplicity
  // step 0 =  All events,       all                            mc primary particles,                    true multiplicity

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
  vector<Float_t> px;
  vector<Float_t> py;
  vector<Float_t> pz;
  vector<Float_t> theta;


  //number of accepted tracks and tracklets
  Int_t ntracks = 0;  //return value for reading functions for ESD and AOD
  //Int_t ntracksRemove = 0;  //return value for reading functions for ESD and AOD
  vector<Int_t> nTracksTracklets; // [0]=nAccepted,1=ntracklets,2=nall(also neutral in case of mc, 
                                  //for real nall=ncharged) 
                                 
 
  if(!fAODEvent && !fESDEvent)return;

  //=================== AOD ===============

  if(fAODEvent){//AOD loop

    //reset global values
    fNMcPrimAccept=0;// number of accepted primaries
    fNRecAccept=0;   // number of accepted tracks
  
    if(CheckEvent(true)){//step 5 = TrigVtxRecNrec, step 4 = TrigVtxRecMcPropNrec ,step 3 = TrigVtxMcNrec
    
      if(!fMcOnly){ 
	//step 5 = TrigVtxRecNrec
	ntracks = LoopAOD(pt, eta, phi, charge, nTracksTracklets, 5);//read tracks
	if(pt.size()){ //(internally ntracks=fNRecAccept)
	  Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1], nTracksTracklets[2], 5);//analyse
	}
      
	if(fUseMC){
	  // step 4 = TrigVtxRecMcPropNrec
	  ntracks = LoopAODRecMcProp(pt, eta, phi, charge, nTracksTracklets, 4);//read tracks
	  if(pt.size()){//(internally ntracks=fNRecAccept)
	    Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1], nTracksTracklets[2], 4);//analyse
	  }
	}
      }
    
      if(fUseMC){
	// step 3 = TrigVtxMcNrec
	ntracks = LoopAODMC(pt, eta, phi, charge, nTracksTracklets, 3);//read tracks
	if(pt.size()){//(internally ntracks=fNRecAccept)
	  Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1],nTracksTracklets[2], 3);//analyse
	}
      }


    }//check event (true)

    //reset values
    fNMcPrimAccept=0;// number of accepted primaries
    fNRecAccept=0;   // number of accepted tracks

    if(CheckEvent(false)){// all events, with and without reconstucted vertex
      ntracks  = LoopAOD  (pt, eta, phi, charge, nTracksTracklets, 2);//need to compute Nrec once more for all events
      ntracks  = LoopAODMC(pt, eta, phi, charge, nTracksTracklets, 1);//read tracks
      if(pt.size()){
	Analyse(pt, eta, phi, charge, fNRecAccept, nTracksTracklets[1],nTracksTracklets[2], 2); // step 2 = TrigAllMcNrec
	
	Analyse(pt, eta, phi, charge, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 1);  // step 1 = TrigAllMcNmc
      }

      //step 0 (include not triggered events) is not possible with AODs generated before October 2011
      //step 0 can be implemented for new AODs

    }
  }//AOD loop


  //=================== ESD ===============

    
  if(fESDEvent){//ESD loop

    //reset values
    fNMcPrimAccept=0;// number of accepted primaries
    fNRecAccept=0;   // number of accepted tracks
  
    // instead of task->SelectCollisionCandidate(mask) in AddTask macro
    Bool_t isSelected = (((AliInputEventHandler*)
			  (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
			 ->IsEventSelected() &    AliVEvent::kMB);

    
    if(fDebug){
      Printf("IsSelected = %d", isSelected);
      Printf("CheckEvent(true)= %d", CheckEvent(true));
      Printf("CheckEvent(false)= %d", CheckEvent(false));
    }

    //check trigger
    if(isSelected){ // has offline trigger
      
      if(CheckEvent(true)){//step 5 = TrigVtxRecNrec, step 4 = TrigVtxRecMcPropNrec ,step 3 = TrigVtxMcNrec
	
	if(!fMcOnly){ 
	  //step 5 = TrigVtxRecNrec
	  ntracks = LoopESD(pt, eta, phi, charge,nTracksTracklets, 5);//read tracks
	  if(pt.size()){ //(internally ntracks=fNRecAccept)
	    Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1], nTracksTracklets[2], 5);//analyse
	  }
	  
	  if(fUseMC){
	    // step 4 = TrigVtxRecMcPropNrec
	    ntracks = LoopESDRecMcProp(pt, eta, phi, charge, nTracksTracklets, 4);//read tracks
	    if(pt.size()){//(internally ntracks=fNRecAccept)
	      Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1], nTracksTracklets[2], 4);//analyse
	    }
	  }
	}
	
	if(fUseMC){
	  // step 3 = TrigVtxMcNrec
	  ntracks = LoopESDMC(pt, eta, phi, charge, nTracksTracklets, 3);//read tracks
	  if(pt.size()){//(internally ntracks=fNRecAccept)
	    Analyse(pt, eta, phi, charge, ntracks, nTracksTracklets[1],nTracksTracklets[2], 3);//analyse
	  }
	}
	
	
      }//check event (true)
      
      //reset values
      fNMcPrimAccept=0;// number of accepted primaries
      fNRecAccept=0;   // number of accepted tracks
      
      if(CheckEvent(false)){// all events, with and without reconstucted vertex
	ntracks  = LoopESD  (pt, eta, phi, charge, nTracksTracklets, 2);//need to compute Nrec once more for all events
	ntracks  = LoopESDMC(pt, eta, phi, charge, nTracksTracklets, 1);//read tracks
	if(pt.size()){
	  Analyse(pt, eta, phi, charge, fNRecAccept, nTracksTracklets[1],nTracksTracklets[2], 2); // step 2 = TrigAllMcNrec
	  
	  Analyse(pt, eta, phi, charge, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 1);  // step 1 = TrigAllMcNmc

	  Analyse(pt, eta, phi, charge, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 0);  //first part of step 0 // step 0 = AllAllMcNmc
	}
	
	
      }
    }// triggered event
    
    else { // not selected by physics selection task = not triggered
      if(CheckEvent(false)){
	ntracks  = LoopESDMC(pt, eta, phi, charge, nTracksTracklets, 0);//read tracks
	if(pt.size())Analyse(pt, eta, phi, charge, fNMcPrimAccept, nTracksTracklets[1],nTracksTracklets[2], 0);  //second part of step 0 // step 0 = AllAllMcNmc
      }
    }


  }//ESD loop





}      


//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopESD( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
				       vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
				       vector<Int_t> &nTracksTracklets, const Int_t step)
{
  // gives back the number of esd tracks and pointer to arrays with track
  // properties (pt, eta, phi)
  // Uses TPC tracks with SPD vertex from now on
  
  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear(); 
  nTracksTracklets.clear(); 
 
  const AliESDVertex*	vtxSPD   = fESDEvent->GetPrimaryVertexSPD(); // uses track or SPD vertexer
  fVertexZ[step]->Fill(vtxSPD->GetZ());
  
  // Retreive the number of all tracks for this event.
  Int_t ntracks = fESDEvent->GetNumberOfTracks();
  if(fDebug>1)  Printf("all ESD tracks: %d", ntracks);

  //first loop to check how many tracks are accepted
  //------------------
  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
 
    AliESDtrack *esdTrack = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
    if (!esdTrack) {
      Error("LoopESD", "Could not receive track %d", iTracks);
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
  if(step==5 || step ==2) fNRecAccept=nAcceptedTracks;
  return fNRecAccept;


}   

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopESDRecMcProp( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
						vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
						vector<Int_t> &nTracksTracklets, const Int_t step)
{  
  // gives back the number of esd tracks and pointer to arrays with track
  // properties (pt, eta, phi) of mc particles if available
  // Uses TPC tracks with SPD vertex from now on

  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear(); 
  nTracksTracklets.clear(); 

  
  AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
  if (!mcEvent) {
    Error("LoopESDRecMcProp", "Could not retrieve MC event");
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
  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
 
    AliVParticle *vtrack = fESDEvent->GetTrack(iTracks);
    AliESDtrack *esdTrack = (AliESDtrack *)fESDEvent->GetTrack(iTracks);
    if (!esdTrack) {
      Error("LoopESDRecMcProp", "Could not receive track %d", iTracks);
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
    if(vtrack->GetLabel()<0){
      if (TMath::Abs(track->Eta())<fEtaCut && track->Pt()>fPtMin && track->Pt()<fPtMax){
	ptArray.push_back(track->Pt());
	etaArray.push_back(track->Eta());
	phiArray.push_back(track->Phi());
	chargeArray.push_back(track->Charge());
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
Int_t AliAnalysisTaskMinijet::LoopESDMC(vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
					vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
					vector<Int_t> &nTracksTracklets, const Int_t step)
{
  // gives back the number of charged prim MC particle and pointer to arrays 
  // with particle properties (pt, eta, phi)

  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear(); 
  nTracksTracklets.clear(); 

  fNMcPrimAccept=0;

  AliMCEvent *mcEvent = (AliMCEvent*) MCEvent();
  if (!mcEvent) {
    Error("LoopESDMC", "Could not retrieve MC event");
    return 0;
  }

  AliStack* stack = MCEvent()->Stack();
  if(!stack) return 0;
  
  Int_t ntracks = mcEvent->GetNumberOfTracks();
  if(fDebug>1)Printf("MC particles: %d", ntracks);

  //vertex
  AliGenEventHeader*  header = MCEvent()->GenEventHeader();
  TArrayF mcV;
  header->PrimaryVertex(mcV);
  Float_t vzMC = mcV[2];
  fVertexZ[step]->Fill(vzMC);

  //----------------------------------
  //first track loop to check how many chared primary tracks are available
  Int_t nChargedPrimaries=0;
  Int_t nAllPrimaries=0;

  Int_t nPseudoTracklets=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
    if (!track) {
      Error("LoopESDMC", "Could not receive track %d", iTracks);
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
    Printf("All in acceptance=%d",nAllPrimaries);
    Printf("Charged in acceptance =%d",nChargedPrimaries);
  }
  
  fChargedPi0->Fill(nAllPrimaries,nChargedPrimaries);

  if(nAllPrimaries==0) return 0;  
  if(nChargedPrimaries==0) return 0;  

  //track loop
  //Int_t iChargedPiK=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliMCParticle *track = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iTracks));
    if (!track) {
      Error("LoopESDMC", "Could not receive track %d", iTracks);
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
  
    
  } //track loop

  nTracksTracklets.push_back(nChargedPrimaries);
  nTracksTracklets.push_back(nPseudoTracklets);
  if(fSelectParticles!=2){
    nTracksTracklets.push_back(nAllPrimaries);
  }
  else{
    nTracksTracklets.push_back(nAllPrimaries-nChargedPrimaries); // only neutral
  }


  fNMcPrimAccept=nChargedPrimaries;

  if(step==1){
    fNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
    fPNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
  }
  
  fVzEvent= mcV[2];
  return fNRecAccept;
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopAOD( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
				       vector<Float_t> &phiArray,  vector<Short_t> &chargeArray,
				       vector<Int_t> &nTracksTracklets, const Int_t step)
{
  // gives back the number of AOD tracks and pointer to arrays with track 
  // properties (pt, eta, phi)

  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear(); 
  nTracksTracklets.clear(); 

  TClonesArray *mcArray=0x0;
  if(fAnalysePrimOnly){
    mcArray = (TClonesArray*)fAODEvent->FindListObject(AliAODMCParticle::StdBranchName());
  }

  
  AliAODVertex*	vertex= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();//GetPrimaryVertex()
  Double_t vzAOD=vertex->GetZ();
  fVertexZ[step]->Fill(vzAOD);
  
  // Retreive the number of tracks for this event.
  Int_t ntracks = fAODEvent->GetNumberOfTracks();
  if(fDebug>1) Printf("AOD tracks: %d", ntracks);
  

  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(iTracks);
    if (!track) {
      Error("LoopAOD", "Could not receive track %d", iTracks);
      continue;
    }
   
    AliVParticle *vtrack = fAODEvent->GetTrack(iTracks);

    //use only tracks from primaries
    if(fAnalysePrimOnly){
      if(vtrack->GetLabel()<0)continue;
      if(!(static_cast<AliAODMCParticle*>(mcArray->At(vtrack->GetLabel()))->IsPhysicalPrimary()))continue;
    }
    
    if(track->TestFilterBit(128) && TMath::Abs(track->Eta())<fEtaCut  
       && track->Pt()>fPtMin && track->Pt()<fPtMax){
      
      nAcceptedTracks++;

      // Printf("dca= %f", track->DCA());
      //save track properties in vector
      ptArray.push_back(track->Pt());
      etaArray.push_back(track->Eta());
      phiArray.push_back(track->Phi());
      chargeArray.push_back(track->Charge());
      fHistPt->Fill(track->Pt());

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
  if(step==5 || step==2)fNRecAccept = nAcceptedTracks; // needed for MC case //step5 = TrigVtxRecNrec
  return fNRecAccept; // at the moment, always return reconstructed multiplicity

}   

//________________________________________________________________________
Int_t AliAnalysisTaskMinijet::LoopAODRecMcProp( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
						vector<Float_t> &phiArray, vector<Short_t> &chargeArray, 
						vector<Int_t> &nTracksTracklets, const Int_t step)
{
  // gives back the number of AOD tracks and pointer to arrays with track 
  // properties (pt, eta, phi)


  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear(); 
  nTracksTracklets.clear(); 
  

  // Retreive the number of tracks for this event.
  Int_t ntracks = fAODEvent->GetNumberOfTracks();
  if(fDebug>1) Printf("AOD tracks: %d", ntracks);

  
  //get array of mc particles
  TClonesArray *mcArray = (TClonesArray*)fAODEvent->
    FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    Printf("No MC particle branch found");
    return kFALSE;
  }

  AliAODVertex*	vtx= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();//GetPrimaryVertex()  
  Double_t vzAOD=vtx->GetZ();
  fVertexZ[step]->Fill(vzAOD);

  Int_t nAcceptedTracks=0;
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliAODTrack *track = (AliAODTrack *)fAODEvent->GetTrack(iTracks);

    AliVParticle *vtrack = fAODEvent->GetTrack(iTracks);

    if (!track) {
      Error("LoopAODRecMcProp", "Could not receive track %d", iTracks);
      continue;
    }
   
    //use only tracks from primaries
    if(fAnalysePrimOnly){
      if(vtrack->GetLabel()<0)continue;
      if(!(static_cast<AliAODMCParticle*>(mcArray->At(vtrack->GetLabel()))->IsPhysicalPrimary()))continue;
    }

    if(track->TestFilterBit(128) &&  TMath::Abs(track->Eta())<fEtaCut &&
       track->Pt()>fPtMin && track->Pt()<fPtMax){
      
      nAcceptedTracks++;

      //save track properties in vector
      if(vtrack->GetLabel()<0){ //fake tracks
	// 	Printf("Fake track");
	// 	continue;
	ptArray.push_back(track->Pt());
	etaArray.push_back(track->Eta());
	phiArray.push_back(track->Phi());
	chargeArray.push_back(track->Charge());
      }
      else{//mc properties
	AliAODMCParticle *partOfTrack = (AliAODMCParticle*)mcArray->At(vtrack->GetLabel());
	
	ptArray.push_back(partOfTrack->Pt());
	etaArray.push_back(partOfTrack->Eta());
	phiArray.push_back(partOfTrack->Phi());
	chargeArray.push_back(vtrack->Charge());//partOfTrack?
      }

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
Int_t AliAnalysisTaskMinijet::LoopAODMC( vector<Float_t> &ptArray,  vector<Float_t> &etaArray, 
					 vector<Float_t> &phiArray, vector<Short_t> &chargeArray,
					 vector<Int_t> &nTracksTracklets, const Int_t step)
{
  // gives back the number of AOD MC particles and pointer to arrays with particle 
  // properties (pt, eta, phi)
  
  ptArray.clear(); 
  etaArray.clear(); 
  phiArray.clear(); 
  chargeArray.clear();
  nTracksTracklets.clear(); 

  //check vertex
  AliAODMCHeader *aodMCheader = (AliAODMCHeader *) fAODEvent->
    FindListObject(AliAODMCHeader::StdBranchName());
  Float_t vzMC = aodMCheader->GetVtxZ();
  fVertexZ[step]->Fill(vzMC);

  
  //retreive MC particles from event
  TClonesArray *mcArray = (TClonesArray*)fAODEvent->
    FindListObject(AliAODMCParticle::StdBranchName());
  if(!mcArray){
    Printf("No MC particle branch found");
    return kFALSE;
  }
  
  Int_t ntracks = mcArray->GetEntriesFast();
  if(fDebug>1)Printf("MC particles: %d", ntracks);


  // Track loop: chek how many particles will be accepted
  //Float_t vzMC=0.;
  Int_t nChargedPrim=0;
  Int_t nAllPrim=0;
  Int_t nPseudoTracklets=0;
  for (Int_t it = 0; it < ntracks; it++) {
    AliAODMCParticle *track = (AliAODMCParticle*)mcArray->At(it);
    if (!track) {
      Error("LoopAODMC", "Could not receive particle %d", it);
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
      Error("LoopAODMC", "Could not receive particle %d", it);
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
  if(step==1){ // step 1 = Trig All Mc Nmc
    fNmcNch->Fill( fNMcPrimAccept,fNRecAccept);
    fPNmcNch->Fill(fNMcPrimAccept,fNRecAccept);
  }
  return fNRecAccept; // rec value from step 5 or step 2
 

} 

//________________________________________________________________________
void AliAnalysisTaskMinijet::Analyse(const vector<Float_t> &pt, 
				     const vector<Float_t> &eta, 
				     const vector<Float_t> &phi, 
				     const vector<Short_t> &charge, 
				     const Int_t ntracksCharged, 
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
      Printf("nCharged=%d",ntracksCharged); 
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
  
  Float_t ptOthers  = 0; // pt others // for all other tracks around event axis -> see loop
  Float_t etaOthers = 0; // eta others
  Float_t phiOthers = 0; // phi others
  Short_t chargeOthers = 0; // charge others

  Int_t   *pindexInnerEta  = new Int_t  [nAll+1];
  Float_t *ptInnerEta      = new Float_t[nAll+1];
  
 

  for (Int_t i = 0; i < nAll; i++) {

    if(pt[i]<0.01)continue;

    //fill single particle correction for first step of pair correction
    Double_t propAll[3] = {eta[i],pt[i],ntracksCharged }; 
    fMapAll[step]->Fill(propAll);
      

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
      fPtSeed[step]    -> Fill( ptEventAxis);
      fEtaSeed[step]   -> Fill(etaEventAxis);
      fPhiSeed[step]   -> Fill(phiEventAxis);


      Double_t prop[3] = {etaEventAxis,ptEventAxis,ntracksCharged }; 
      fMapSingleTrig[step]->Fill(prop);

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

   
	fDPhiDEtaEventAxis[step]->Fill(dPhi, dEta);
	fDPhiEventAxis[step]->Fill(dPhi);

	//check outliers
	if(ptEventAxis< 0.4 || ptEventAxis > 100) Printf("particles out of range pt");
	if(ntracksCharged<0 || ntracksCharged>150) Printf("particles out of range ncharge");
	if(TMath::Abs(dEta)>2*fEtaCut) {
	  Printf("particles out of range dEta");
	  Printf("eta1=%f, eta2=%f", etaOthers, etaEventAxis);
	  Printf("step=%d",step);
	}
	if(dPhi<-0.5*TMath::Pi() || dPhi>1.5*TMath::Pi()){
	  Printf("particles out of range dPhi");
	  Printf("phi1=%f, phi2=%f", phiOthers, phiEventAxis);
	}
	    
	Bool_t isLikeSign = CheckLikeSign(chargeEventAxis, chargeOthers);
	    
	Double_t prop6[6] = {ptEventAxis,ptOthers,dEta,dPhi,ntracksCharged, isLikeSign }; 
	fMapPair[step]->Fill(prop6);
  
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
void AliAnalysisTaskMinijet::Terminate(Option_t*)
{
  //terminate function is called at the end
  //can be used to draw histograms etc.


}

//________________________________________________________________________
Bool_t AliAnalysisTaskMinijet::SelectParticlePlusCharged(const Short_t charge, const Int_t pdg, Bool_t prim)
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
    Printf("Error: wrong selection of charged/pi0/k0");
    return 0;
  }

  return true;

}

//________________________________________________________________________
Bool_t AliAnalysisTaskMinijet::SelectParticle(const Short_t charge, const Int_t pdg, const Bool_t prim)
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
Bool_t AliAnalysisTaskMinijet::CheckEvent(const Bool_t recVertex)
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
      //AliGenEventHeader*  headerg = MCEvent()->GenEventHeader();
      //TArrayF mcVg;
      //headerg->PrimaryVertex(mcVg);
      //  if(TMath::Abs(mcVg[0])<1e-8 && TMath::Abs(mcVg[0])<1e-8 && 
      // TMath::Abs(mcVg[0])<1e-8) return false;
      // Float_t vzMCg = mcVg[2];
      // if(TMath::Abs(vzMCg)>fVertexZCut) return false;
    }

    //rec
    if(recVertex==true){
      if(fESDEvent->IsPileupFromSPD(3,0,8))return false;
      
      //rec vertex
      // const AliESDVertex*	vertexESDg   = fESDEvent->GetPrimaryVertex(); // uses track or SPD vertexer
      // if(!vertexESDg) return false;
      // if(vertexESDg->GetNContributors()<=0)return false;
      // Float_t fVzg= vertexESDg->GetZ();
      // if(TMath::Abs(fVzg)>fVertexZCut) return false;
      
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
      //mc
      //   AliAODMCHeader *aodMCheader = (AliAODMCHeader *) fAODEvent->
      // 	 FindListObject(AliAODMCHeader::StdBranchName());
      //        Float_t vzMC = aodMCheader->GetVtxZ();
      //        if(TMath::Abs(vzMC)>fVertexZCut) return false;
       
      //retreive MC particles from event
      TClonesArray *mcArray = (TClonesArray*)fAODEvent->
	FindListObject(AliAODMCParticle::StdBranchName());
      if(!mcArray){
	Printf("No MC particle branch found");
	return false;
      }
    }

    //rec
    if(recVertex==true){
      if(fAODEvent->IsPileupFromSPD(3,0.8))return false;
       
      //      AliAODVertex*	vertex= (AliAODVertex*)fAODEvent->GetPrimaryVertex();
      //      if(!vertex) return false;
      //      if(vertex->GetNContributors()<=0) return false;
      //      Double_t vzAOD=vertex->GetZ();
      //      if(TMath::Abs(vzAOD)<1e-9) return false;
      //      if(TMath::Abs(vzAOD)>fVertexZCut) return false;
       
      AliAODVertex*	vertexSPD= (AliAODVertex*)fAODEvent->GetPrimaryVertexSPD();
      if(!vertexSPD) return false;
      if(vertexSPD->GetNContributors()<=0) return false;
      Double_t vzSPD=vertexSPD->GetZ();
      //if(TMath::Abs(vzSPD)<1e-9) return false;
      if(TMath::Abs(vzSPD)>fVertexZCut) return false;
    }
    return true;
   
  }

  else {
    Printf("Wrong mode!");
    return false;
  }
  
}

//_____________________________________________________________________________
const Double_t * AliAnalysisTaskMinijet::CreateLogAxis(const Int_t nbins, 
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
Bool_t AliAnalysisTaskMinijet::CheckLikeSign(const Short_t chargeEventAxis, 
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
    Printf("Error: Charge not lower nor higher as zero");
    return false;
  }
  
  Printf("Error: Check values of Charge");
  return false;
}
 
