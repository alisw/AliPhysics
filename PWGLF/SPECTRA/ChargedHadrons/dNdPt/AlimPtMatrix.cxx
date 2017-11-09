//------------------------------------------------------------------------------
// AlimPtMatrix 
// 
// generates the correlation matrix between N_acc and N_ch
// creates only one THnSparse: pT_acc:eta_acc:N_acc:N_ch
// 
// made by P. Luettig
// last modified 2011/08/05
// modified my M. Marquard
//------------------------------------------------------------------------------


#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
// #include "THnSparse.h"
#include "TParticle.h"
#include "TRandom3.h"
#include "iostream"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#include "AliHeader.h"

#include "AliTriggerAnalysis.h"
#include "AliPhysicsSelection.h"

#include "AlimPtMatrix.h"

ClassImp(AlimPtMatrix)
	//________________________________________________________________________
	AlimPtMatrix::AlimPtMatrix(const char *name) 
: AliAnalysisTaskSE(name), 
	fESD(0),
	fMC(0),
	fOutputList(0),
	r3(0),
	fNaccCent(0),
	fNchCent(0),
	fNaccNch(0),
	fEvtsCent(0),
	fESDTrackCuts(0),
	fMaxVertexZ(0),
	fUseCentrality(0),
	iNbEvents(0),
	bUseMC(0),
	fCentLimit(0),
	fSubsample(0),
	fSubsampleEvent(1),
	fSubsampleTrack(1),
	fTrigger(AliTriggerAnalysis::kMB1),
	genMinPt(0),
	genMaxPt(0),
	genMinEta(0),
	genMaxEta(0)
{
	// Constructor
	iNbEvents = 0;
	// Define input and output slots here
	// Input slot #0 works with a TChain
	//  DefineInput(0, TChain::Class());
	// Output slot #0 writes into a TList
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AlimPtMatrix::~AlimPtMatrix()
{
	if (r3) delete r3; r3 = 0;
	if (fOutputList) {
		fOutputList->Clear();
		delete fOutputList;
	}
	fOutputList = 0;

	if (fESDTrackCuts) delete fESDTrackCuts; fESDTrackCuts = 0;

}

// //________________________________________________________________________
// void AlimPtMatrix::ConnectInputData(Option_t *) 
// {
//   // Connect ESD or AOD here
//   // Called once
// 
//   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
//   if (!tree) {
//     Printf("ERROR: Could not read chain from input slot 0");
//   } else {
//     // Disable all branches and enable only the needed ones
//     // The next two lines are different when data produced as AliESDEvent is read
//     /*
//        tree->SetBranchStatus("*", kFALSE);
//        tree->SetBranchStatus("fTracks.*", kTRUE);
//        */
// 
//     AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
//     if (!esdH) 
//     {
//       Printf("ERROR: Could not get ESDInputHandler");
//     } 
//     else fESD = esdH->GetEvent();
// 
//     if(bUseMC)
//     {
//       AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
//       if (!mcH) 
//       {
// 	printf("ERROR: Could not get MCInputHandler");
// 	return;
//       }
//       else fMC = mcH->MCEvent(); 
// 
//     }
//   }
// }

//________________________________________________________________________
void AlimPtMatrix::UserCreateOutputObjects()
{
	// Create histograms
	// Called once

	OpenFile(1, "RECREATE");

	//const Int_t multNbins = 3500;
	//const Int_t evtsNbins = 1;
	//  const Int_t ptNbins = 68;
	//  const Int_t etaNbins = 30;
	//const Int_t centNbins = 11;
	//  const Int_t xvNbins = 12;
	//  const Int_t yvNbins = 12;
	//  const Int_t zvNbins = 12;
	//
	//  const Int_t finexvNbins = 100;
	//  const Int_t fineyvNbins = 100;
	//
	//  const Double_t finexvRange = 10.;
	//  const Double_t fineyvRange = 10.;
	//
	//  Double_t fineBinsXv[finexvNbins+1] = {0.};
	//  Double_t fineBinsYv[fineyvNbins+1] = {0.};
	//
	// for(Int_t i = 0; i <= finexvNbins+1; i++) fineBinsXv[i] = -1. * finexvRange + (Double_t)i*(2. * finexvRange)/(Double_t)finexvNbins;
	//  for(Int_t i = 0; i <= fineyvNbins+1; i++) fineBinsYv[i] = -1. * fineyvRange + (Double_t)i*(2. * fineyvRange)/(Double_t)fineyvNbins;

	//Double_t binsMult[multNbins+1] = {-0.5, 0.5 , 1.5 , 2.5 , 3.5 , 4.5 , 5.5 , 6.5 , 7.5 , 8.5,
	//  9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5, 17.5, 18.5,
	//  19.5,20.5, 21.5, 22.5, 23.5, 24.5, 29.5, 149.5};

	//Double_t binsMult[multNbins+1] = {0.};

	//for(Int_t i = 0; i <= multNbins+1; i++) binsMult[i] = -0.5 + (Double_t)i;


	//  Double_t binsPt[ptNbins+1] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
	//    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 
	//    1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 
	//    2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5,
	//    5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0,
	//    14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 
	//    32.0, 34.0, 36.0, 40.0, 45.0, 50.0};

	//  Double_t binsEta[etaNbins+1] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,
	//	                               0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
	//  
	//Double_t binsCent[centNbins+1] = {0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.};
	//Double_t binsEvents[evtsNbins+1] = {0., 2.};
	//  
	//  Double_t binsXv[xvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
	//  Double_t binsYv[yvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
	//  Double_t binsZv[zvNbins+1] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};

	//  cout << endl << endl << genMinPt << " " << genMaxPt << " " << genMinEta << " " << genMaxEta << endl;

	fOutputList = new TList();
	fOutputList->SetOwner();

	// define array of bins
	//Int_t binsNaccNch[3] = {multNbins, multNbins, centNbins};
	//  Int_t binsEvtsCent[2] = {evtsNbins, centNbins};



	// Histograms
	//fNaccNchCent = new THnSparseF("fNaccNchCent","N_acc:N_ch:Centrality",3,binsNaccNch);
	//fNaccNchCent->SetBinEdges(0, binsMult);
	//fNaccNchCent->SetBinEdges(1, binsMult);
	//fNaccNchCent->SetBinEdges(2, binsCent);
	//fNaccNchCent->GetAxis(0)->SetTitle("n_{acc}");
	//fNaccNchCent->GetAxis(1)->SetTitle("n_{ch}");
	//fNaccNchCent->GetAxis(2)->SetTitle("Centrality");

	fNaccCent = new TH2F("fNaccCent","N_acc:Centrality",3500,-0.5,3499.5,12,-1.5,10.5);
	fNaccCent->GetXaxis()->SetTitle("n_{acc}");
	fNaccCent->GetYaxis()->SetTitle("Centrality");

	fNchCent = new TH2F("fNchCent","N_ch:Centrality",3500,-0.5,3499.5,12,-1.5,10.5);
	fNchCent->GetXaxis()->SetTitle("n_{ch}");
	fNchCent->GetYaxis()->SetTitle("Centrality");

	fNaccNch = new TH2F("fNaccNch","Nacc:Nch",3500,-0.5, 3499.5,3500,-0.5,3499.5);
	fNaccNch->GetXaxis()->SetTitle("n_{ch}");
	fNaccNch->GetYaxis()->SetTitle("n_{acc}");


	fEvtsCent = new TH1F("fEvtsCent","Centrality;Number of Events",12,-1.5,10.5);

	// Add Histos, Profiles etc to List
	fOutputList->Add(fNaccCent);
	fOutputList->Add(fNchCent);
	fOutputList->Add(fNaccNch);
	fOutputList->Add(fEvtsCent);

	// Post output data (if histograms are not used later, PostData is at least called here)
	PostData(1, fOutputList);

}

//________________________________________________________________________
void AlimPtMatrix::UserExec(Option_t *) 
{
	// Main loop
	// Called for each event

	//if(bUseMC) Printf("INFO: Running on MC");
	// 	else Printf("INFO: Running on data");

	Double_t centralityVZERO;
	//Int_t centralityVZEROBin;

	//Int_t centBin = 11;
	Int_t centBin = -1;

	Int_t nacc = 0;
	Int_t nch = 0;
	Double_t sumpt = 0.;

	Double_t partPt = 0.;
	Double_t partEta = 0.;
	//Double_t rnum = 0.;

	//   TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
	//   if (!tree) {
	//     Printf("ERROR: Could not read chain from input slot 0");
	//   } else {
	// Disable all branches and enable only the needed ones
	// The next two lines are different when data produced as AliESDEvent is read
	/*
	   tree->SetBranchStatus("*", kFALSE);
	   tree->SetBranchStatus("fTracks.*", kTRUE);
	 */
	switch (fSubsample){
		case 1:{
			       r3 = new TRandom3(0);
			       Double_t randomevent = r3->Rndm();
			       if (fSubsampleEvent<=0 || fSubsampleEvent>1){
				       fSubsampleEvent = 1;
			       }
			       if (randomevent > fSubsampleEvent){
				       return;
			       }
		       }break;
		case 2:{
			       r3 = new TRandom3(0);
		       }break;
		case 3:{
			       r3 = new TRandom3(0);
			       Double_t randomevent = r3->Rndm();
			       if (fSubsampleEvent<=0 || fSubsampleEvent>1){
				       fSubsampleEvent = 1;
			       }
			       if (randomevent > fSubsampleEvent){
				       return;
			       }
		       }break;
	};


	//if (fSubsample != 0){
	//	r3 = new TRandom3(0);
	//}

	//if (fSubsample == 1 || fSubsample == 3){
	//	Double_t randomevent = r3->Rndm();
	//	if (fSubsampleEvent<=0 || fSubsampleEvent>1){
	//		fSubsampleEvent = 1;
	//	}
	//	if (randomevent > fSubsampleEvent){
	//		return;
	//	}
	//}


	AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if (!esdH) 
	{
		Printf("ERROR: Could not get ESDInputHandler");
	} 
	else fESD = (AliESDEvent*)esdH->GetEvent();

	if(bUseMC)
	{
		AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
		if (!mcH) 
		{
			printf("ERROR: Could not get MCInputHandler");
			return;
		}
		else fMC = mcH->MCEvent(); 

	}
	//   }


	if (!fESD) 
	{
		Printf("ERROR: fESD not available");
		return;
	}

	if ( bUseMC && !fMC)
	{
		Printf("ERROR: fMC not available");
		return;
	}

	// trigger selection
	Bool_t isEventTriggered = kFALSE;
	AliPhysicsSelection *physicsSelection = NULL;
	AliTriggerAnalysis* triggerAnalysis = NULL;
	// always MB
	//isEventTriggered = inputHandler->IsEventSelected() & AliVEvent::kMB;
	/*cout << "EventSel: " << esdH->IsEventSelected() << " " << GetCollisionCandidates() << endl;
	  isEventTriggered = esdH->IsEventSelected() & GetCollisionCandidates();
	 */
	physicsSelection = static_cast<AliPhysicsSelection*> (esdH->GetEventSelection());
	if(!physicsSelection) return;
	//SetPhysicsTriggerSelection(physicsSelection);

	//if (isEventTriggered ) {
	// set trigger (V0AND)
	triggerAnalysis = physicsSelection->GetTriggerAnalysis();
	if(!triggerAnalysis) return;
	//cout << isEventTriggered <<" " << GetTrigger() << " ";
	isEventTriggered = triggerAnalysis->IsOfflineTriggerFired(fESD, GetTrigger());
	//cout << isEventTriggered << " " << iNbEvents << endl;
	//}

	if (isEventTriggered == 0) {
		//cout<<"No trigger found" << endl;
		return;
	}


	if( (fESD->GetNumberOfTracks() == 0 ) )  
	{
		//cout << "No tracks" << endl;
		return;
	}

	if(!fESDTrackCuts) 
	{
		Printf("ERROR: No esd track cuts available");
	}

	// track cuts from Jochen
	// const AliESDVertex* vtxESDTPC = fESD->GetPrimaryVertexTPC();
	//  if( vtxESDTPC->GetNContributors() < 1 ) 
	//  {
	//		return;
	//	}

	const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();
	if( vtxESD->GetNContributors() < 1 ) 
	{
		//cout << "No contributors" << endl;
		return;
	}

	// // francesco prino cut
	//  const AliMultiplicity* multESD = fESD->GetMultiplicity();
	//  if( vtxESDTPC->GetNContributors() < (-10.+0.25*multESD->GetNumberOfITSClusters(0)) ) 
	//  {
	//    return;
	//  } 

	//const AliESDVertex* vtxESD = fESD->GetPrimaryVertexTracks();

	// Event cut on the z-Position of the vertex
	if(vtxESD->GetZ() > fMaxVertexZ || vtxESD->GetZ() < (-1.*fMaxVertexZ)) 
	{
		//Printf("VertexZ out of range, Zv = %f",vtxESD->GetZ());
		return;
	}

	if(fUseCentrality != 0)
	{

		AliCentrality *esdCentrality = fESD->GetCentrality();


		//V0
		if (fUseCentrality == 1)
		{

			centBin = esdCentrality->GetCentralityClass10("V0M");	
			//centralityVZERO = esdCentrality->GetCentralityPercentile("V0M");

			//if      ( centralityVZERO >   0. && centralityVZERO <   5.) centBin =  0;
			//else if ( centralityVZERO >=  5. && centralityVZERO <  10.) centBin =  1;
			//else if ( centralityVZERO >= 10. && centralityVZERO <  20.) centBin =  2;
			//else if ( centralityVZERO >= 20. && centralityVZERO <  30.) centBin =  3;
			//else if ( centralityVZERO >= 30. && centralityVZERO <  40.) centBin =  4;
			//else if ( centralityVZERO >= 40. && centralityVZERO <  50.) centBin =  5;
			//else if ( centralityVZERO >= 50. && centralityVZERO <  60.) centBin =  6;
			//else if ( centralityVZERO >= 60. && centralityVZERO <  70.) centBin =  7;
			//else if ( centralityVZERO >= 70. && centralityVZERO <  80.) centBin =  8;
			//else if ( centralityVZERO >= 80. && centralityVZERO <  90.) centBin =  9;
			//else if ( centralityVZERO >= 90. && centralityVZERO <  99.) centBin = 10;
			//else if ( centralityVZERO >= 99. ) centBin = 11;
			//else if ( centralityVZERO <= 0.  ) centBin = 11;

			//cout<<"cent class 10: "<< cent10 <<" cent percentile: "<< centBin <<" diff: "<< cent10-centBin << endl;
		} // end use centrality
	} 
	else
	{
		//         printf("No Centrality!\n");
		centBin = 0; // for pp
		//         return;
	}

	if (centBin < 0 || centBin > fCentLimit){
		//if (centBin >= fCentLimit){
		//cout << "No centrality" << endl;
		return;
	}

	iNbEvents++;

	if (bUseMC)
	{
		// data loop

		// select only primary charged particles
		TArrayF vtxMC(3);
		AliStack *stack = fMC->Stack();
		if(!stack) 
		{
			AliDebug(AliLog::kError, "Stack not available");
			return;
		}
		// get MC event header
		AliHeader *header = fMC->Header();
		if (!header) 
		{
			AliDebug(AliLog::kError, "Header not available");
			return;
		}

		// get MC vertex
		AliGenEventHeader* genHeader = header->GenEventHeader();
		if (!genHeader) 
		{
			AliDebug(AliLog::kError, "Could not retrieve genHeader from Header");
			return;
		}
		genHeader->PrimaryVertex(vtxMC);		

		// MC Track loop to fill a pT spectrum
		for (Int_t iMCTracks = 0; iMCTracks < stack->GetNtrack(); iMCTracks++) 
		{
			switch (fSubsample){
				case 2:{
					       Double_t randomtrack = r3->Rndm();
					       if (fSubsampleTrack<=0 || fSubsampleTrack>1){
						       fSubsampleTrack = 1;
					       }
					       if (randomtrack > fSubsampleTrack){
						       continue;
					       }
				       }break;
				case 3:{
					       Double_t randomtrack = r3->Rndm();
					       if (fSubsampleTrack<=0 || fSubsampleTrack>1){
						       fSubsampleTrack = 1;
					       }
					       if (randomtrack > fSubsampleTrack){
						       continue;
					       }	
				       }break;
			};

			//if (fSubsample == 2 || fSubsample == 3){
			//	Double_t randomtrack = r3->Rndm();
			//	if (fSubsampleTrack<=0 || fSubsampleTrack>1){
			//		fSubsampleTrack = 1;
			//	}
			//	if (randomtrack > fSubsampleTrack){
			//		continue;
			//	}
			//}
			TParticle* particle = stack->Particle(iMCTracks);
			if (!particle) continue;

			// only charged particles
			if(!particle->GetPDG()) continue;
			Double_t charge = particle->GetPDG()->Charge()/3.;
			if ( TMath::Abs(charge) < 0.001 ) continue;

			//Int_t partPdg = TMath::Abs(particle->GetPdgCode());

			// select only pi, K, p
			//if( (partPdg != 211) && (partPdg != 321) && (partPdg != 2212)) continue;
			// physical primary
			Bool_t prim = stack->IsPhysicalPrimary(iMCTracks);
			if(!prim) continue;

			// if(!AcceptParticle(particle)) continue;


			partPt = particle->Pt();
			partEta = particle->Eta();

			//if( partPt < genMinPt) continue;
			//if( partPt > genMaxPt) continue;
			if( TMath::Abs(partEta) > genMaxEta) continue;


			nch++;

		} //MC track loop

		fNchCent->Fill(nch, centBin);
	} // end MC loop

	// loop over all reconstructed tracks
	for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) 
	{
		switch (fSubsample){
			case 2:{
				       Double_t randomtrack = r3->Rndm();
				       if (fSubsampleTrack<=0 || fSubsampleTrack>1){
					       fSubsampleTrack = 1;
				       }
				       if (randomtrack > fSubsampleTrack){
					       continue;
				       }
			       }break;
			case 3:{
				       Double_t randomtrack = r3->Rndm();
				       if (fSubsampleTrack<=0 || fSubsampleTrack>1){
					       fSubsampleTrack = 1;
				       }
				       if (randomtrack > fSubsampleTrack){
					       continue;
				       }	
			       }break;
		};

		//if (fSubsample == 2 || fSubsample == 3){
		//	Double_t randomtrack = r3->Rndm();
		//	if (fSubsampleTrack<=0 || fSubsampleTrack>1){
		//		fSubsampleTrack = 1;
		//	}
		//	if (randomtrack > fSubsampleTrack){
		//		continue;
		//	}
		//}

		AliESDtrack* track = fESD->GetTrack(iTracks);

		if (!track) 
		{
			Printf("ERROR: Could not receive track %d", iTracks);
			continue;
		}

		if(!fESDTrackCuts->AcceptTrack(track)) continue;

		sumpt += track->Pt();

		nacc++;      
	} //track reconstructed loop

	//Double_t mult[3] = {nacc, nch, centBin};
	fNaccCent->Fill(nacc, centBin);

	if(bUseMC) fNaccNch->Fill(nch, nacc);

	//Double_t dNbEvents[2] = {1, centBin};
	fEvtsCent->Fill(centBin);


	PostData(1, fOutputList);

	if (r3) delete r3; r3 = 0;
	}

	//________________________________________________________________________
	void AlimPtMatrix::Terminate(Option_t *) 
	{
		// Draw result to the screen
		// Called once at the end of the query

		fOutputList = dynamic_cast<TList*>(GetOutputData(1));
		if (!fOutputList) {
			Printf("ERROR: fOutputList not available");
			return;
		}

		//cout << "Number of events: " << iNbEvents << endl;
	}

	// //________________________________________________________________________
	// // Cuts for generated Particles
	// Bool_t AlimPtMatrix::AcceptParticle(TParticle *p)
	// {
	//   if(!p) return kFALSE;
	// 
	//   if( (p->Pt() < genMinPt) || 
	//       (p->Pt() > genMaxPt) || 
	//       (p->Eta() < genMinEta) || 
	//       (p->Eta() > genMaxEta) ) return kFALSE;
	// 
	//   return kTRUE;
	// 
	// }
