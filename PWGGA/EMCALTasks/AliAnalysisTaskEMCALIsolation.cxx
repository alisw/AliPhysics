#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"

#include "AliVEvent.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"

#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"

#include "AliAnalysisTaskEMCALIsolation.h"

#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

// Task for isolating gammas with EMCAL
// Author: Marco Marquard

ClassImp(AliAnalysisTaskEMCALIsolation)

	//________________________________________________________________________

AliAnalysisTaskEMCALIsolation::AliAnalysisTaskEMCALIsolation()
	: AliAnalysisTaskSE(),
	bVerbose(0),
	bMC(0),
	fESD(0),
	fAOD(0),
	fMC(0),
	fStack(0),
	fOutputList(0),
	fHistGlobalHmap(0),
	fHistGlobalHmap0(0),
	fTreeEvent(0),
	fTreeCluster(0),
	emclus(0),
	prodrad(0),
	contPID(-999),
	mothPID(-999),
	trackmatch(0),
	fGeom(0),
	fGeoName("EMCAL_COMPLETE"),
	fESDCells(0),
	fAODCells(0),
	fEsdClusters(0),
	fAodClusters(0)
{
	//default constructor
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEMCALIsolation::AliAnalysisTaskEMCALIsolation(const char *name) 
	: AliAnalysisTaskSE(name),
	bVerbose(0),
	bMC(0),
	fESD(0),
	fAOD(0),
	fMC(0),
	fStack(0),
	fOutputList(0),
	fHistGlobalHmap(0),
	fHistGlobalHmap0(0),
	fTreeEvent(0),
	fTreeCluster(0),
	emclus(0),
	prodrad(0),
	contPID(-999),
	mothPID(-999),
	trackmatch(0),
	fGeom(0),
	fGeoName("EMCAL_COMPLETE"),
	fESDCells(0),
	fAODCells(0),
	fEsdClusters(0),
	fAodClusters(0)
{
	// Constructor

	// Define input and output slots here
	// Input slot #0 works with a TChain
	DefineInput(0, TChain::Class());
	// Output slot #0 id reserved by the base class for AOD
	// Output slot #1 writes into a TH1 container
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEMCALIsolation::~AliAnalysisTaskEMCALIsolation() 
{
	//Destructor

	delete fESD;
	delete fAOD;
	delete fOutputList;
	delete fHistGlobalHmap;
	delete fHistGlobalHmap0;
	delete fTreeCluster;
	delete fTreeEvent;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsolation::UserCreateOutputObjects()
{
	// Create histograms
	// Called once
	Double_t etagaps[97]= {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};

	Double_t phigaps[125]= {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732,1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

	fOutputList = new TList();

	fHistGlobalHmap = new TH2F("fHistGlobalHmap","Hit map for EMCAL",96,etagaps,124,phigaps);
	fHistGlobalHmap->GetXaxis()->SetTitle("#eta");
	fHistGlobalHmap->GetYaxis()->SetTitle("#phi");
	fOutputList->Add(fHistGlobalHmap);

	fHistGlobalHmap0 = new TH2F("fHistGlobalHmap0","Hit map for EMCAL",96,0,96,120,0,120);
	fHistGlobalHmap0->GetXaxis()->SetTitle("#eta");
	fHistGlobalHmap0->GetYaxis()->SetTitle("#phi");
	fOutputList->Add(fHistGlobalHmap0);
	/*	
		fHist = new TH1F("","",,,,);
		fHist->GetXaxis()->SetTitle("");
		fHist->GetYaxis()->SetTitle("");
		fHist->SetMarkerStyle(kFullCircle);
		fOutputList->Add(fHist);
		*/

	fTreeEvent = new TTree("eventsT","eventsT");
	fTreeEvent->Branch("emclus",	&emclus);

	fTreeCluster = new TTree("clusterT","clusterT");
	fTreeCluster->Branch("emclus",	&emclus);
	fTreeCluster->Branch("prodrad",	&prodrad);
	fTreeCluster->Branch("contPID",	&contPID);
	fTreeCluster->Branch("mothPID",	&mothPID);

	fOutputList->Add(fHistGlobalHmap);
	fOutputList->Add(fHistGlobalHmap0);
	fOutputList->Add(fTreeEvent);
	fOutputList->Add(fTreeCluster);



}

//________________________________________________________________________
void AliAnalysisTaskEMCALIsolation::UserExec(Option_t *) 
{

	AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();


	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
	if (!esdH)
	{
		Printf("ERROR: Could not get ESDInputHandler");
	}

	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
	if (!aodH && !esdH)
	{
		Printf("ERROR: Could not get AODInputHandler");
	}

	if(esdH){
	  fESD = (AliESDEvent*)esdH->GetEvent();
		am->LoadBranch("AliESDRun.");
		am->LoadBranch("AliESDHeader.");
	}
	else if(aodH){
		fAOD = aodH->GetEvent();
	}
	else{
		AliFatal("Neither ESD nor AOD event found");
		return;
	}

	AliMCEventHandler* mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	if (!mcH) {
		Printf("ERROR: Could not retrieve MC event handler");
		return;
	}
	fMC = mcH->MCEvent();
	if (!fMC) {
		Printf("ERROR: Could not retrieve MC event");
		return;
	}
	// stack to get truth infos of particles:
	fStack = fMC->Stack();
	if (!fStack) {
		Printf("ERROR: Could not retrieve stack from MC event");
		return;
	}
	//Printf(" MC particles: %d", fMC->GetNumberOfTracks());


	fGeom = AliEMCALGeometry::GetInstance(fGeoName);
	//AliEMCALEMCGeometry *emc = fGeom->GetEMCGeometry();
	//Double_t phimin = emc->GetArm1PhiMin();
	//Double_t phimax = emc->GetArm1PhiMax();
	//printf("EMCAL coverage in phi from: %f to: %f \n",phimin,phimax);
	if (!fGeom){
		AliFatal("EMCAL(fGeom) geometry is not initialized");
		return;
	}

	AliVCaloCells *cells = 0;

	if (fESD){
		am->LoadBranch("EMCALCells."); 
		fESDCells = fESD->GetEMCALCells();
		cells = fESDCells;
	}
	else{
		am->LoadBranch("emcalCells");
		fAODCells = fAOD->GetEMCALCells();
		cells = fAODCells;
	}




	if (!cells){
		AliFatal("No cells loaded");
		return;
	}

	if(fESD){ 
		am->LoadBranch("CaloClusters");
		TList *l = fESD->GetList();
		fEsdClusters = dynamic_cast<TClonesArray*>(l->FindObject("CaloClusters"));
	}
	else if(fAOD){
		am->LoadBranch("caloClusters");
		TList *l = fAOD->GetList();
		fAodClusters = dynamic_cast<TClonesArray*>(l->FindObject("caloClusters"));
	}
	else{
		AliFatal("Impossible to not have either pointer to ESD or AOD event");
	}

	TObjArray *clusters = fEsdClusters;
	if (!clusters)
		clusters = fAodClusters;
	if (!clusters){
		Printf("ERROR: no clusters node in event!");
		//return;
	}

	//========================================= start of analysis ============================================


	//============> DATA LOOP =================>

	//------------> cell loop ----------------->
	Int_t ncells = cells->GetNumberOfCells();
	//for (Int_t i = 0; i<15000; ++i){
	//	Int_t absID    = i;
	for (Int_t i = 0; i<ncells; ++i){
		Int_t absID    = TMath::Abs(cells->GetCellNumber(i));

		//Check if this absId exists
		if(!(fGeom->CheckAbsCellId(absID))){ 
			continue;
		}
		Int_t nSupMod, nModule, nIphi, nIeta, iphi, ieta;
		// Get from the absId the super module number, the module number and the eta-phi index (0 or 1) in the module
		fGeom->GetCellIndex(absID, nSupMod, nModule, nIphi, nIeta);

		// Get from the  the super module number, the module number and the eta-phi index (0 or 1) in the module the tower row (iphi) and column (ieta)
		fGeom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi, ieta);

		if(bVerbose)
			printf("ID: %d ; supermodule: %d ; module %d ; phi %d ; eta %d \n",absID, nSupMod, nModule, nIphi, nIeta);


		Int_t etainv = 47-ieta;
		//fHSMCellE[nSupMod]->Fill(etainv,phi,i);


		Int_t binx = etainv + 48*((nSupMod+1) % 2);
		Int_t biny = iphi + 24*nSupMod/2;
		if(nSupMod % 2 != 0)
			biny -= 12;

		Double_t cellE = cells->GetAmplitude(i);
		//fHistGlobalHmap0->Fill(binx,biny,cellE);

		Double_t nEtaGlobal, nPhiGlobal; 
		Int_t nSupMod2;

		nSupMod2 = fGeom->GetSuperModuleNumber(absID);
		fGeom->EtaPhiFromIndex(absID,nEtaGlobal,nPhiGlobal);

		if(bVerbose){
			printf("SM number: %d ; global eta: \t%6.4f ; global phi: \t%6.4f \n",nSupMod2,nEtaGlobal,nPhiGlobal);
		}


		//Double_t nPhi;
		//if(nPhiGlobal < 0){
		//	nPhi = TMath::Pi()+TMath::Abs(TMath::Pi()+nPhiGlobal);
		//}
		//else{
		//nPhi = nPhiGlobal;
		//}
		fHistGlobalHmap0->Fill(binx,biny,cellE);
		fHistGlobalHmap->Fill(nEtaGlobal,nPhiGlobal,cellE);

	}

	if(bVerbose){
		printf("\nThere are %d cells in this event\n", ncells);
		printf("There are %d tracks in this event\n", fESD->GetNumberOfTracks());
	}


	//<------------ cell loop <-----------------

	//------------> track loop ----------------->
	for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
		AliESDtrack* track = fESD->GetTrack(iTracks);
		if (!track) {
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}


	}
	//<------------ track loop <-----------------

	//------------> cluster loop ----------------->
	emclus = 0;
	Int_t nclus = clusters->GetEntries();
	//printf("%d clusters\n",nclus);
	//cluster loop 1 for counting EMCal cluster
	for(int i = 0; i < nclus; i++){ 
		AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
		if (!clus){
			//printf("Cluster empty!\n");
			continue;
		}
		if (!clus->IsEMCAL()){
			//printf("Cluster %d is not EMCAL!\n",i);
			//      cout << "at " << i << " cluster not EMCAL" << endl;
			continue;
		}

		emclus++;
	}

	//<------------ cluster loop <-----------------


	//<============ DATA LOOP <=================


	//============> MC LOOP =================>
	if(bMC){

		//------------> cell loop ----------------->
		//<------------ cell loop <-----------------

		//------------> track loop ----------------->
		//<------------ track loop <-----------------

		//------------> cluster loop  ----------------->
		for(int i = 0; i < nclus; i++){
			AliVCluster *clus = static_cast<AliVCluster*>(clusters->At(i));
			if (!clus){
				//printf("Cluster empty!\n");
				continue;
			}
			if (!clus->IsEMCAL()){
				//printf("Cluster %d is not EMCAL!\n",i);
				//      cout << "at " << i << " cluster not EMCAL" << endl;
				continue;
			}
			//MCTruth
			int clslabel = -1;

			// all contributors
			Int_t* mcarr =  clus->GetLabels();
			// number of contributors
			Int_t nl = clus->GetNLabels();




			for(int ip=0;ip<nl;ip++){
				Int_t entry = mcarr[ip];
				AliMCParticle *mcPart = static_cast<AliMCParticle*>(fMC->GetTrack(entry));
				if(!mcPart){
					continue;
				}
				contPID = mcPart->PdgCode();
				if(bVerbose){
					printf("%d \t %d \n",clslabel,mcPart->PdgCode());
				}
				prodrad = mcPart->Particle()->R();
				Int_t indexMoth = mcPart->GetMother();
				if(indexMoth >= 0)
				{
					AliMCParticle* moth = static_cast<AliMCParticle*>(fMC->GetTrack(indexMoth));
					//	//AliMCParticle* moth = (AliMCParticle*)stack->At(indexMoth);
					mothPID = moth->PdgCode();
				}





				//printf("%d \t %d \n",clslabel,nl);

			}
			fTreeCluster->Fill();

		}

		//<------------ cluster loop <-----------------

	}
	//<============ MC LOOP <=================







	//cluster loop 2 for analysis
	//\cluster loop


	fTreeEvent->Fill();


	PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskEMCALIsolation::Terminate(Option_t *) 
{
	// Draw result to the screen
	// Called once at the end of the query

	fOutputList = dynamic_cast<TList*> (GetOutputData(1));
	if (!fOutputList) {
		printf("ERROR: Output list not available\n");
		return;
	}

	fHistGlobalHmap = dynamic_cast<TH2F*> (fOutputList->At(2));
	if (!fHistGlobalHmap) {
		printf("ERROR: fHistGlobalHmap not available\n");
		return;
	}

	fHistGlobalHmap->DrawCopy("colz");

}


const char * AliAnalysisTaskEMCALIsolation::GetParticleName(Int_t pdg)
{
	TParticlePDG * p1 = TDatabasePDG::Instance()->GetParticle(pdg);
	if(p1) return p1->GetName();
	return Form("%d", pdg);
}
