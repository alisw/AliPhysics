/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id:$ */

////////////////////////////////////////////////////////////////////////
//
// Analysis class to Correct Underlying Event studies
//
// This class needs as input ESDs.
// The output is an analysis-specific container.
//
// The class is used to get the contamination from secondaries
// from tracks DCA distribution 
// as function of track pT and pseudo-rapidity.
// It provides additional information for the corrections 
// that can not be retrieved by the AliAnalysisTaskLeadingTackUE
// task, which is running on AODs.
// 
////////////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TChain.h>
//#include <TCanvas.h>
//#include <TFile.h>
#include <TList.h>
#include <TMath.h>
//#include <TProfile.h>
#include <TTree.h>
//#include <TVector3.h>
#include <TH3F.h>

#include "AliAnalyseLeadingTrackUE.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskCorrectionsUE.h"
#include "AliAODHandler.h"
#include "AliESDHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODInputHandler.h"
#include "AliESDInputHandler.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisHelperJetTasks.h"

class TCanvas;
class TFile;
class TProfile;
class TVector3; 

ClassImp( AliAnalysisTaskCorrectionsUE)

// Define global pointer
AliAnalysisTaskCorrectionsUE* AliAnalysisTaskCorrectionsUE::fgTaskCorrectionsUE=NULL;

//____________________________________________________________________
AliAnalysisTaskCorrectionsUE:: AliAnalysisTaskCorrectionsUE(const char* name):
AliAnalysisTask(name,""),
fAnalyseUE(0x0),
fDebug(0),
fESDEvent(0x0),
fESDHandler(0x0),
fInputHandler(0x0),
fListOfHistos(0x0),  
fMcEvent(0x0),
fMcHandler(0x0),
fMode(0),
fOutCFcont(0x0),
fhEntries(0x0),
fhFakes(0x0),
fhPtMCAll(0x0),
fhPtMCPrim(0x0),
fhPtMCSec(0x0),
fhPtMCPrimFake(0x0),
fhPtMCSecFake(0x0),
fnTracksVertex(1),  // QA tracks pointing to principal vertex  
fZVertex(10.),
fhVertexContributors(0x0),
fhVertexReso(0x0),
fTrackEtaCut(0.9),
fTrackPtCut(0.),
fEsdTrackCuts(0x0),
fEsdTrackCutsSPD(0x0),
fEsdTrackCutsSDD(0x0),
fEsdTrackCutsDCA(0x0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());

}

//______________________________________________________________
AliAnalysisTaskCorrectionsUE & AliAnalysisTaskCorrectionsUE::operator = (const AliAnalysisTaskCorrectionsUE & /*source*/)
{
  // assignment operator
  return *this;
}


/************** INTERFACE METHODS *****************************/

//______________________________________________________________
Bool_t AliAnalysisTaskCorrectionsUE::Notify()
{
  
  return kTRUE;

}

//____________________________________________________________________
void AliAnalysisTaskCorrectionsUE::ConnectInputData(Option_t* /*option*/)
{
  // Connect the input data  
 
  if (fDebug > 1) AliInfo("ConnectInputData() ");
  
  //Get the input handler
  fESDHandler = (AliESDInputHandler*)AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
  if ( !fESDHandler && fDebug > 0 ) {
  	AliFatal(" No ESD event handler connected !!! ");
	return;
	}

  //Get ESD event
   fESDEvent = (AliESDEvent*)fESDHandler->GetEvent();
   if (!fESDEvent && fDebug > 1) {
  	AliFatal(" No ESD event retrieved !!! ");
        return;
	}

  //Get MC handler 
  fMcHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  
  // Define track cuts
  fEsdTrackCuts = new AliESDtrackCuts("ITSTPC", "ITS+TPC standard 2009 cuts w.o. SPD cluster requirement nor DCA cut");
  // TPC  
  fEsdTrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  fEsdTrackCuts->SetMinNClustersTPC(70);
  //fEsdTrackCuts->SetMinNClustersTPC(90); // ***** TMP *****
  fEsdTrackCuts->SetMaxChi2PerClusterTPC(4);
  fEsdTrackCuts->SetAcceptKinkDaughters(kFALSE);
  fEsdTrackCuts->SetRequireTPCRefit(kTRUE);
  // ITS
  fEsdTrackCuts->SetRequireITSRefit(kTRUE);
  fEsdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  fEsdTrackCuts->SetMaxDCAToVertexZ(2); // new for pile-up !!!

  // Add SPD requirement 
  fEsdTrackCutsSPD = new AliESDtrackCuts("SPD", "Require 1 cluster in SPD");
  fEsdTrackCutsSPD->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
  
  // Add SDD requirement 
  fEsdTrackCutsSDD = new AliESDtrackCuts("SDD", "Require 1 cluster in first layer SDD");
  fEsdTrackCutsSDD->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kFirst);
  
  // Add DCA cuts 
  fEsdTrackCutsDCA = new AliESDtrackCuts("DCA", "pT dependent DCA cut");
  // 7*(0.0050+0.0060/pt^0.9)
  fEsdTrackCutsDCA->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
  
  fEsdTrackCutsDCA->SetMaxDCAToVertexZ(1.e6);
  fEsdTrackCutsDCA->SetDCAToVertex2D(kFALSE);

  // emulates filterbit when getting leading-track from ESD
  fAnalyseUE->DefineESDCuts(16); // any number = standard ITS+TPC + (SPD or SDD)
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::CreateOutputObjects()
{
  // Create the output container
  
  if (fDebug > 1) AliInfo("CreateOutPutData()");
   
  // Create pointer to AliAnalyseLeadingTrackUE, a class implementing the main analysis algorithms
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  if (!fAnalyseUE){
     AliError("UE analysis class not initialized!!!");
     return;
  }
  
  // Create list of output histograms
  if (fListOfHistos != NULL){
	delete fListOfHistos;
        fListOfHistos = NULL;
  	}
  if (!fListOfHistos){
  	fListOfHistos = new TList();
  	fListOfHistos->SetOwner(kTRUE); 
  	}
  
  // Create CF container
  CreateContainer();
  // number of events
  fhEntries = new TH1F("fhEntries","Entries",1,0.,2.); 
  fListOfHistos->Add(fhEntries);
  // tracks contributing to the vertex
  fhVertexContributors = new TH1F("fhVertexContributors","Tracks in vertex",51, -0.5, 50.5);
  fhVertexContributors->GetXaxis()->SetTitle("# tracks in vertex");
  fhVertexContributors->GetYaxis()->SetTitle("Entries");
  fListOfHistos->Add(fhVertexContributors);
  // vertex resolution
  fhVertexReso = new TH3F("fhVertexReso","Vertex resolution",51,-0.5,50.5,100,0.,0.05,100.,0.,0.1); 
  fhVertexContributors->GetXaxis()->SetTitle("# tracks in vertex");
  fhVertexContributors->GetYaxis()->SetTitle("Resolution XY (cm)");
  fhVertexContributors->GetZaxis()->SetTitle("Resolution Z (cm)");
  fListOfHistos->Add(fhVertexReso);
  
  // number of fake tracks
  fhFakes = new TH1F("fhFakes","Fraction of fake tracks",5,-0.5,4.5);
  fhFakes->GetXaxis()->SetBinLabel(1,"No MC");
  fhFakes->GetXaxis()->SetBinLabel(2,"Unique MC primary");
  fhFakes->GetXaxis()->SetBinLabel(3,"Unique MC secondary");
  fhFakes->GetXaxis()->SetBinLabel(4,"Multiple MC");
  fListOfHistos->Add(fhFakes);
  
  //pT distributions
  fhPtMCAll = new TH1F("fhPtMCAll","All MC particles reconstructed",100,0., 20.);
  fListOfHistos->Add(fhPtMCAll);
  fhPtMCPrim = new TH1F("fhPtMCPrim","Primary MC particles reconstructed",100,0., 20.);
  fListOfHistos->Add(fhPtMCPrim);
  fhPtMCSec = new TH1F("fhPtMCSec","Secondary MC particles reconstructed",100,0., 20.);
  fListOfHistos->Add(fhPtMCSec);
  fhPtMCPrimFake = new TH1F("fhPtMCPrimFake","Fake primary MC particles reconstructed",100,0., 20.);
  fListOfHistos->Add(fhPtMCPrimFake);
  fhPtMCSecFake = new TH1F("fhPtMCSecFake","Fake secondary MC particles reconstructed",100,0., 20.);
  fListOfHistos->Add(fhPtMCSecFake);


  // Add task configuration to output list 
  AddSettingsTree();
    

  //Post outputs
  PostData(0,fListOfHistos);
 
}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::Exec(Option_t */*option*/)
{
  
  // Get MC event
  if (fMcHandler){
  	fMcEvent = fMcHandler->MCEvent();
       	if ( fDebug > 3 ) AliInfo( " Processing MC event..." );
  	if (fMode && !fMcEvent) return;
	}
  
  // Do the analysis
  AnalyseCorrectionMode();

  
 PostData(0,fListOfHistos);

}

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  
  if (fDebug >1) AliAnalysisHelperJetTasks::PrintDirectorySize("PWG4_JetTasksOutput.root");
  
  if (!gROOT->IsBatch()){
  	fListOfHistos = dynamic_cast<TList*> (GetOutputData(0));
  	if (!fListOfHistos){
		AliError("Histogram List is not available");
		return;
  		}

  } else {
        AliInfo(" Batch mode, not histograms will be shown...");
  }

  
}


/******************** ANALYSIS METHODS *****************************/

//____________________________________________________________________
void  AliAnalysisTaskCorrectionsUE::AddSettingsTree()
{
  //Write settings to output list
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fnTracksVertex", &fnTracksVertex, "TracksInVertex/D");
  settingsTree->Branch("fZVertex", &fZVertex, "VertexZCut/D");
  settingsTree->Branch("fTrackPtCut", &fTrackPtCut, "TrackPtCut/D");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  

//____________________________________________________________________
void AliAnalysisTaskCorrectionsUE::AnalyseCorrectionMode()
{

  // Analyse the event
  Int_t labelMC  = -1;
  if (fMcHandler && fMcEvent){
  	// Only consider MC events within the vtx-z region used also as cut on the reconstructed vertex
  	if (!fAnalyseUE->VertexSelection(fMcEvent, 0, fZVertex))  
    	return; 
  	// Get MC-true leading particle 
  	TObjArray *ltMC = (TObjArray*)fAnalyseUE->FindLeadingObjects(fMcEvent);
  	AliVParticle* leadingMC = 0;
  	if (ltMC){
    		leadingMC = (AliVParticle*) ltMC->At(0);
		}
	if (!leadingMC)return; 
	labelMC = TMath::Abs(leadingMC->GetLabel());
	}

  // Trigger selection ************************************************
  if (!fAnalyseUE->TriggerSelection(fESDHandler)) return;
  
  // PILEUP-CUT ****** NEW !!!!!!!!! **************
  Bool_t select = kTRUE;
  //select = AliAnalysisHelperJetTasks::TestSelectInfo(AliAnalysisHelperJetTasks::kIsPileUp);
  if (! select) return;
  
  // Vertex selection *************************************************
  
  if(!fAnalyseUE->VertexSelection(fESDEvent, 0, fZVertex)) return;
  AliESDVertex* vertex = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
  Int_t nvtx = vertex->GetNContributors();
  fhVertexContributors->Fill(nvtx);
  if (fMcHandler){
		AliVVertex *vertexMC = (AliVVertex*)fMcEvent->GetPrimaryVertex();
		if (vertexMC){
			Double_t diffx = vertexMC->GetX()-vertex->GetX();
			Double_t diffy = vertexMC->GetY()-vertex->GetY();
			Double_t diffxy = TMath::Sqrt(TMath::Power(diffx,2)+TMath::Power(diffy,2));
			//Double_t diffxy = TMath::Abs(diffx); // **** TMP ****
			//Double_t diffxy = diffy;
			Double_t diffz = TMath::Abs(vertexMC->GetZ()-vertex->GetZ());
		
			fhVertexReso->Fill(nvtx,diffxy,diffz);
		}else if (fDebug>1)Printf("******* NO MC VERTEX ********");
  	}

  if(!fAnalyseUE->VertexSelection(fESDEvent, fnTracksVertex, fZVertex)) return;

  // Get Reconstructed leading particle *******************************
  TObjArray *ltRECO = fAnalyseUE->FindLeadingObjects(fESDEvent);
  if (!ltRECO){
  	delete[] ltRECO;
	return;
 	}
  Int_t labelReco= TMath::Abs(((AliVParticle*)ltRECO->At(0))->GetLabel());
  
  
  // Loop on tracks
  Int_t nTracks = fESDEvent->GetNumberOfTracks();
  if (!nTracks)return;
  // count accepted events
  fhEntries->Fill(1.);

  Int_t npart=0;
  Bool_t *labelsArray = 0; 
  if (fMcHandler){
  	npart = fMcEvent->GetNumberOfTracks();
        labelsArray = new Bool_t[npart];
  	for(Int_t j = 0; j<npart; j++){
		labelsArray[j] = kFALSE;
  		}
	}

  for (Int_t i = 0; i < nTracks; i++){
	AliESDtrack *track = fESDEvent->GetTrack(i);
	// only charged
	if (!track || !track->Charge())continue;
        // apply cuts
	Bool_t cut = fEsdTrackCuts->AcceptTrack(track);
	Bool_t cutSPD = fEsdTrackCutsSPD->AcceptTrack(track);
	Bool_t cutSDD = fEsdTrackCutsSDD->AcceptTrack(track);
	Bool_t cutDCA = fEsdTrackCutsDCA->AcceptTrack(track);
  	
	//Exclude the MC leading track
	Double_t matchLeading = 1.;//no match
	if (fMcHandler){
		if (TMath::Abs(track->GetLabel())==labelMC) matchLeading=0.; //match MC leading
		if (TMath::Abs(track->GetLabel())==labelReco) {
			matchLeading=2.;//match RECO leading
			if (labelMC == labelReco)matchLeading = 3.; // match both (mc = reco leading)
			}
		}
	// Fill step0 (all tracks) 
	FillContainer(track, 0,kFALSE,matchLeading);         
	// Fill step1 (no SPD cluster requirement - no DCA cut )
	if ( cut ) FillContainer(track,1,kFALSE,matchLeading);
	// Fill step2
	if ( cut && cutDCA && (cutSPD || cutSDD)) FillContainer(track,2,kFALSE,matchLeading);
	// Fill step3-step4-step5 
	if ( cut && cutDCA && !cutSPD && !cutSDD) FillContainer(track,3,kTRUE,matchLeading);
	if ( cut && cutDCA && cutSPD) FillContainer(track,4,kTRUE,matchLeading);
	if ( cut && cutDCA && !cutSPD && cutSDD) FillContainer(track,5,kTRUE,matchLeading);
        // Fill step 6 - temporary just to define the standard track cuts  
	if ( cut && cutSPD ) FillContainer(track,6,kFALSE,matchLeading);
	//if ( cut && cutSPD ) FillContainer(track,6,kTRUE,matchLeading); // ***** TMP *****
	if ( cut && (cutSPD || cutSDD) ) FillContainer(track,7,kFALSE,matchLeading);
	// Study contamination form fakes
	if (fMcHandler){
			
		// consider standard ITS+TPC cuts without SPD cluster requirement
		if (cut && cutDCA){
			//check if it points back to any MC
			Int_t label = TMath::Abs(track->GetLabel());
	       		AliVParticle *part = (AliVParticle*)fMcEvent->GetTrack(label);	
			if (!part){
				fhFakes->Fill(0.);
				//Printf("*************** NO MC PARTICLE ************************");
				continue;
				}
				
	                fhPtMCAll->Fill(part->Pt()); 
			// Check if label is not already in array
			if (!labelsArray[label]){
				labelsArray[label]= kTRUE;
				if (fMcEvent->IsPhysicalPrimary(label)){
					fhFakes->Fill(1.);
					fhPtMCPrim->Fill(part->Pt());
				}else{
					fhFakes->Fill(2.);
					fhPtMCSec->Fill(part->Pt());
					}
			}else{
				fhFakes->Fill(3.);
				if (fMcEvent->IsPhysicalPrimary(label))fhPtMCPrimFake->Fill(part->Pt());
				else fhPtMCSecFake->Fill(part->Pt());
				}
		   }	
		}
  	} // end loop on tracks
	if(labelsArray)
	delete[] labelsArray;

}


//____________________________________________________________________
void AliAnalysisTaskCorrectionsUE::CreateContainer()
{
  
  // Create the output CF container
  // relevant variables 
  UInt_t iprim  = 0; // 0: primaries, 1: secondaries from strangness, 2: secondaries form material 
  UInt_t iptmc  = 1;
  UInt_t ipt    = 2;
  UInt_t ieta   = 3;
  UInt_t idcaxy = 4;
  UInt_t idcaz  = 5;
  UInt_t imatch = 6;
  UInt_t icharge = 7;
  // set-up the grid
  UInt_t nstep = 8;
  const Int_t nvar  = 8;
  const Int_t nbin0 = 5;  // prim
  const Int_t nbin1 = 20; // pt resolution
  const Int_t nbin2 = 39; // pt 
  const Int_t nbin3 = 20; // eta
  const Int_t nbin4 = 100; // dca xy
  const Int_t nbin5 = 100; // dca z
  const Int_t nbin6 = 4; // matching with leading track
  const Int_t nbin7 = 2;

  // array for the number of bins in each dimension
  Int_t iBin[nvar];
  iBin[0] = nbin0;
  iBin[1] = nbin1;
  iBin[2] = nbin2;
  iBin[3] = nbin3;
  iBin[4] = nbin4;
  iBin[5] = nbin5;
  iBin[6] = nbin6;
  iBin[7] = nbin7;
  
  // primaries
  Double_t primBins[7] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5};
  // matching with leading-track
  Double_t matchBins[5] = {-0.5,0.5,1.5,2.5,3.5};

  // pT resolution
  Double_t resoBins[nbin1+1];
  for (Int_t i=0; i<=nbin1; i++)
    resoBins[i] = -1.0 + 0.1 * i;
  
  // pT
  Double_t ptBins[nbin2+1] = {0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 100.0};
  
  // eta
  Double_t etaBins[nbin3+1];
  for (Int_t i=0; i<=nbin3; i++)
    etaBins[i] = -1.0 + 0.1 * i;
  
  // dca xy
  Double_t dcaxyBins[nbin4+1];
  for (Int_t i=0; i<=nbin4; i++)
    dcaxyBins[i] = -1.+0.02 * i;
  
  // dca z
  Double_t dcazBins[nbin5+1];
  for (Int_t i=0; i<=nbin5; i++)
    dcazBins[i] = -5.0 + 0.1 * i;
 
  // charge 
  Double_t chargeBins[nbin7+1] = {-1.,0.,1.};

  // create container
  // set variables
  fOutCFcont = new AliCFContainer("fOutCFcont","Output Container",nstep,nvar,iBin);
  fOutCFcont->SetBinLimits(iprim,primBins);
  fOutCFcont->SetVarTitle(iprim, "Particle type");
  fOutCFcont->SetBinLimits(iptmc,resoBins);
  fOutCFcont->SetVarTitle(iptmc, "#Delta p_{T} (DATA-MC) (GeV/c)");
  fOutCFcont->SetBinLimits(ipt,ptBins);
  fOutCFcont->SetVarTitle(ipt, "p_{T} (GeV/c)");
  fOutCFcont->SetBinLimits(ieta,etaBins);
  fOutCFcont->SetVarTitle(ieta, "#eta");
  fOutCFcont->SetBinLimits(idcaxy,dcaxyBins);
  fOutCFcont->SetVarTitle(idcaxy, " DCA_{XY} (cm)");
  fOutCFcont->SetBinLimits(idcaz,dcazBins);
  fOutCFcont->SetVarTitle(idcaz, " DCA_{Z} (cm)");
  fOutCFcont->SetBinLimits(imatch,matchBins);
  fOutCFcont->SetVarTitle(imatch, "Matching with leading-track");
  fOutCFcont->SetBinLimits(icharge,chargeBins);
  fOutCFcont->SetVarTitle(icharge, "Charge");

  // set steps
  fOutCFcont->SetStepTitle(0,"all tracks");
  fOutCFcont->SetStepTitle(1,"ITS+TPC cuts (no SPD cluster requirement and DCA cut)");
  fOutCFcont->SetStepTitle(2,"add DCA cut");
  fOutCFcont->SetStepTitle(3,"NO SPD cluster, NO SDD cluster in first layer");
  fOutCFcont->SetStepTitle(4,"YES SPD cluster, NO SDD cluster in first layer");
  fOutCFcont->SetStepTitle(5,"NO SPD cluster, YES SDD cluster in first layer");
  fOutCFcont->SetStepTitle(6,"ITS+TPC cuts - no DCA cut - SPD cut");
  fOutCFcont->SetStepTitle(7,"ITS+TPC cuts - no DCA cut - SPD or SDD cut");
  fListOfHistos->Add(fOutCFcont);

}


void  AliAnalysisTaskCorrectionsUE::FillContainer(AliESDtrack *track, Int_t step,Bool_t mcvertex, Double_t matchLeading)
{

  // Fill the CF container

  Double_t vars[8];
  Double_t prim = -1.;
  
  if (track->Charge() > 0.) vars[7] = 0.5;
  else vars[7] = -0.5;
  
  if (fMcHandler){
  	// determine if points back to a primary
  	Int_t label = TMath::Abs(track->GetLabel());
  	
	AliMCParticle *part = (AliMCParticle*)fMcEvent->GetTrack(label);  
	for (Int_t i=0; i<=6;i++) vars[i] = -999.;
  	if (part) { //PRIMARY
  		if (fMcEvent->IsPhysicalPrimary(label)){
			prim = 0.;
		}else { //SECONDARY
			// decide if strange
			Int_t labelm = TMath::Abs(part->GetMother());
			AliMCParticle *mother = (AliMCParticle*)fMcEvent->GetTrack(labelm);  
			Int_t code = mother->PdgCode();
			Int_t mfl = Int_t (code/ TMath::Power(10, Int_t(TMath::Log10(code))));
			 if  (mfl == 3) prim = 1.;
				else{   
				        //Printf("***** PROCESS : %d",part->Particle()->GetUniqueID());  
					if (TMath::Abs(code) == 211) prim = 2.; // charged pion decay
					else if (part->Particle()->GetUniqueID() == 13 )prim = 3.; // hadronic interactions
					else if (part->Particle()->GetUniqueID() == 5 )prim = 4.; // photon conversions
					else prim = 5.; // other?
				}
			}
  		vars[1]= part->Pt()-track->Pt();
  		// In step 2 fill MC pT for contamination study
		if (step == 2)vars[2]=part->Pt();
		}
 	}
  vars[0]=prim;

  if (step != 2)vars[2]=track->Pt();
  vars[3]=track->Eta();

  Bool_t dcaControlFlag = kFALSE;
  if (mcvertex && fMcHandler){
          // we want DCA w.r.t. MC vertex
	  AliVVertex *vtxMC = (AliVVertex*)fMcEvent->GetPrimaryVertex();
	  dcaControlFlag = track->RelateToVertex((AliESDVertex*)vtxMC,(Double_t)fESDEvent->GetMagneticField(),10000.);
	  }else{
	  AliESDVertex* vertex = (AliESDVertex*)fESDEvent->GetPrimaryVertex();
	  dcaControlFlag = track->RelateToVertex(vertex,(Double_t)fESDEvent->GetMagneticField(),10000.);
	  }
  if (dcaControlFlag){
  	Float_t dca[2];
  	Float_t dcaCov[2];
  	track->GetImpactParameters(dca,dcaCov);
 	vars[4]=dca[0];
  	vars[5]=dca[1];
	}	  


  vars[6]= matchLeading;
  fOutCFcont->Fill(vars,step);

}


//____________________________________________________________________
AliAnalysisTaskCorrectionsUE* AliAnalysisTaskCorrectionsUE::Instance()
{ 
  
  //Create instance
  if (fgTaskCorrectionsUE) {
	return fgTaskCorrectionsUE;
  } else {
	fgTaskCorrectionsUE = new AliAnalysisTaskCorrectionsUE();
	return fgTaskCorrectionsUE;
  	}
}

