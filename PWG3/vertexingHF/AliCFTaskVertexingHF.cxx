/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables
// 6 Steps introduced: MC, MC Acc, Reco, Reco Acc, Reco Acc + ITS Cl, 
// Reco Acc + ITS Cl + PPR cuts
// 12 variables used: pt, y, cosThetaStar, ptPi, ptK, ct,
// dca, d0Pi, d0K, d0Pixd0K, cosPointingAngle, phi
//
//-----------------------------------------------------------------------
// Author : C. Zampolli, CERN
//          D. Caffarri, Univ & INFN Padova  caffarri@pd.infn.it
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Base class for HF Unfolding (pt and eta)
// correlation matrix filled at Acceptance and PPR level
// Author: A.Grelli ,  Utrecht - agrelli@uu.nl
//----------------------------------------------------------------------- 
#include <TCanvas.h>
#include <TParticle.h>
#include <TDatabasePDG.h>
#include <TH1I.h>
#include <TStyle.h>
#include <TFile.h>

#include "AliCFTaskVertexingHF.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "AliESDtrackCuts.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliCFVertexingHF2Prong.h"
#include "AliCFVertexingHF3Prong.h"
#include "AliCFVertexingHF.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

//__________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF() :
	AliAnalysisTaskSE(),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCorrelation(0x0),
	fCountMC(0),
	fCountAcc(0),
	fCountVertex(0),
	fCountRefit(0),
	fCountReco(0),
	fCountRecoAcc(0),
	fCountRecoITSClusters(0),
	fCountRecoPPR(0),
	fCountRecoPID(0),
	fEvents(0),
	fDecayChannel(0),
	fFillFromGenerated(kFALSE),
	fOriginDselection(0),
	fAcceptanceUnf(kTRUE),
	fCuts(0),
	fUseWeight(kFALSE),
	fWeight(1.),
	fNvar(0),
	fPartName(""),
	fDauNames(""),
	fSign(2)
{
	//
	//Default ctor
	//
}
//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const Char_t* name, AliRDHFCuts* cuts) :
	AliAnalysisTaskSE(name),
	fCFManager(0x0),
	fHistEventsProcessed(0x0),
	fCorrelation(0x0),
	fCountMC(0),
	fCountAcc(0),
	fCountVertex(0),
	fCountRefit(0),
	fCountReco(0),
	fCountRecoAcc(0),
	fCountRecoITSClusters(0),
	fCountRecoPPR(0),
	fCountRecoPID(0),
	fEvents(0),
	fDecayChannel(0),
	fFillFromGenerated(kFALSE),
	fOriginDselection(0),
	fAcceptanceUnf(kTRUE),
	fCuts(cuts), 
	fUseWeight(kFALSE),
	fWeight(1.),
	fNvar(0),
	fPartName(""),
	fDauNames(""),
	fSign(2)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	/*
	  DefineInput(0) and DefineOutput(0)
	  are taken care of by AliAnalysisTaskSE constructor
	*/
	DefineOutput(1,TH1I::Class());
	DefineOutput(2,AliCFContainer::Class());
	DefineOutput(3,THnSparseD::Class());
	DefineOutput(4,AliRDHFCuts::Class());
	
	fCuts->PrintAll();
}

//___________________________________________________________________________
AliCFTaskVertexingHF& AliCFTaskVertexingHF::operator=(const AliCFTaskVertexingHF& c) 
{
	//
	// Assignment operator
	//
	if (this!=&c) {
		AliAnalysisTaskSE::operator=(c) ;
		fCFManager  = c.fCFManager;
		fHistEventsProcessed = c.fHistEventsProcessed;
		fCuts = c.fCuts;
	}
	return *this;
}

//___________________________________________________________________________
AliCFTaskVertexingHF::AliCFTaskVertexingHF(const AliCFTaskVertexingHF& c) :
	AliAnalysisTaskSE(c),
	fCFManager(c.fCFManager),
	fHistEventsProcessed(c.fHistEventsProcessed),
	fCorrelation(c.fCorrelation),
	fCountMC(c.fCountMC),
	fCountAcc(c.fCountAcc),
	fCountVertex(c.fCountVertex),
	fCountRefit(c.fCountRefit),
	fCountReco(c.fCountReco),
	fCountRecoAcc(c.fCountRecoAcc),
	fCountRecoITSClusters(c.fCountRecoITSClusters),
	fCountRecoPPR(c.fCountRecoPPR),
	fCountRecoPID(c.fCountRecoPID),
	fEvents(c.fEvents),
	fDecayChannel(c.fDecayChannel),
	fFillFromGenerated(c.fFillFromGenerated),
	fOriginDselection(c.fOriginDselection),
	fAcceptanceUnf(c.fAcceptanceUnf),
	fCuts(c.fCuts),
	fUseWeight(c.fUseWeight),
	fWeight(c.fWeight),
	fNvar(c.fNvar),
	fPartName(c.fPartName),
	fDauNames(c.fDauNames),
	fSign(c.fSign)
{
	//
	// Copy Constructor
	//
}

//___________________________________________________________________________
AliCFTaskVertexingHF::~AliCFTaskVertexingHF() 
{
	//
	//destructor
	//
	if (fCFManager)           delete fCFManager ;
	if (fHistEventsProcessed) delete fHistEventsProcessed ;
	if (fCorrelation)	  delete fCorrelation ;
	if (fCuts)                delete fCuts;
}

//_________________________________________________________________________-
void AliCFTaskVertexingHF::Init()
{
	//
	// Initialization
	//
	
	if (fDebug>1) printf("AliCFTaskVertexingHF::Init()");
	AliRDHFCuts *copyfCuts = 0x0;
	if (!fCuts){
		AliFatal("No cuts defined - Exiting...");
		return;
	}

	switch (fDecayChannel){
	case 2:{
		copyfCuts = new AliRDHFCutsD0toKpi(*(dynamic_cast<AliRDHFCutsD0toKpi*>(fCuts)));
		fNvar = 13;
		fPartName="D0";
		fDauNames="K+pi";
		break;
	}
	case 21:{ 
		copyfCuts = new AliRDHFCutsDStartoKpipi(*(dynamic_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
		fNvar = 13;
		fPartName="Dstar";
		fDauNames="K+pi+pi";
		break;
	}
	case 31:{
		copyfCuts = new AliRDHFCutsDplustoKpipi(*(dynamic_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
		fNvar = 12;
		fPartName="Dplus";
		fDauNames="K+pi+pi";
		break;
	}
	case 32:{
		copyfCuts = new AliRDHFCutsLctopKpi(*(dynamic_cast<AliRDHFCutsLctopKpi*>(fCuts)));
		fNvar = 12;
		fPartName="Lambdac";
		fDauNames="p+K+pi";
		break;
	}
	case 33:{
		copyfCuts = new AliRDHFCutsDstoKKpi(*(dynamic_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
		fNvar = 12;
		fPartName="Ds";
		fDauNames="K+K+pi";
		break;
	}
	case 4:{
		copyfCuts = new AliRDHFCutsD0toKpipipi(*(dynamic_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
		fNvar = 13;
		fPartName="D0";
		fDauNames="K+pi+pi+pi";
		break;
	}
	default:
		AliFatal("The decay channel MUST be defined according to AliCFVertexing::DecayChannel - Exiting...");
		break;
	}  
	
	const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
	if (copyfCuts){
		copyfCuts->SetName(nameoutput);
		
		//Post the data
		PostData(4, copyfCuts);
	}
	else{
		AliFatal("Failing initializing AliRDHFCuts object - Exiting...");
	}	
	
	return;
}

//_________________________________________________
void AliCFTaskVertexingHF::UserExec(Option_t *)
{
	//
	// Main loop function
	//
	
	PostData(1,fHistEventsProcessed) ;
	PostData(2,fCFManager->GetParticleContainer()) ;
	PostData(3,fCorrelation) ;
	
	AliESDtrackCuts* trackCuts = fCuts->GetTrackCuts();
	
	if (fFillFromGenerated){
		AliWarning("Flag to fill container with generated value ON ---> dca, d0pi, d0K, d0xd0, cosPointingAngle will be set as dummy!");
	}
	
	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
	
	TClonesArray *arrayBranch=0;
	
	if(!aodEvent && AODEvent() && IsStandardAOD()) {
		// In case there is an AOD handler writing a standard AOD, use the AOD 
		// event in memory rather than the input (ESD) event.    
		aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
		// in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
		// have to taken from the AOD event hold by the AliAODExtension
		AliAODHandler* aodHandler = (AliAODHandler*) 
			((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
		if(aodHandler->GetExtensions()) {
			AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
			AliAODEvent *aodFromExt = ext->GetAOD();
			
			switch (fDecayChannel){
			case 2:{
				arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
				break;
			}
			case 21:{ 
				arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
				break;
			}
			case 31:
			case 32:
			case 33:{
				arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
				break;
			}
			case 4:{
				arrayBranch=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
				break;
			}
			default:
				break;
			}
		}
	} 
	else {
		switch (fDecayChannel){
		case 2:{
			arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
			break;
		}
		case 21:{ 
			arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Dstar");
			break;
		}
		case 31:
		case 32:
		case 33:{
			arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Charm3Prong");
			break;
		}
		case 4:{
			arrayBranch=(TClonesArray*)aodEvent->GetList()->FindObject("Charm4Prong");
			break;
		}
		default:
			break;
		}
	}
	
	AliAODVertex *aodVtx = (AliAODVertex*)aodEvent->GetPrimaryVertex();
	if (!aodVtx) return;
	
	if (!arrayBranch) {
		AliError("Could not find array of HF vertices");
		return;
	}
	
	fEvents++;

	fCFManager->SetRecEventInfo(aodEvent);
	fCFManager->SetMCEventInfo(aodEvent);
	
	//******** DEFINE number of variables of the container***** for now set at 13, in the future in the config macro.
	
	//loop on the MC event
	
	TClonesArray* mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	if (!mcArray) {
		AliError("Could not find Monte-Carlo in AOD");
		return;
	}
	Int_t icountMC = 0;
	Int_t icountAcc = 0;
	Int_t icountReco = 0;
	Int_t icountVertex = 0;
	Int_t icountRefit = 0;
	Int_t icountRecoAcc = 0;
	Int_t icountRecoITSClusters = 0;
	Int_t icountRecoPPR = 0;
	Int_t icountRecoPID = 0;
	Int_t cquarks = 0;
		
	AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	if (!mcHeader) {
		AliError("Could not find MC Header in AOD");
		return;
	}

	Double_t* containerInput = new Double_t[fNvar];
	Double_t* containerInputMC = new Double_t[fNvar]; 
	
       
       	AliCFVertexingHF* cfVtxHF=0x0;
	switch (fDecayChannel){
	case 2:{
		cfVtxHF = new AliCFVertexingHF2Prong(mcArray, fOriginDselection);
		break;
	}
	case 21:{ 
		//	  cfVtxHF = new AliCFVertexingHFCascade(mcArray, originDselection);  // not there yet
		break;
	}
	case 31:
	case 32:
	case 33:{
	  cfVtxHF = new AliCFVertexingHF3Prong(mcArray, fOriginDselection, fDecayChannel); 
		break;
	}
	case 4:{
		//cfVtxHF = new AliCFVertexingHF4Prong(mcArray, originDselection);  // not there yet
		break;
	}
	default:
		break;
	}
	if (!cfVtxHF){
		AliError("No AliCFVertexingHF initialized");
		return;
	}
	
	Double_t zPrimVertex = aodVtx ->GetZ();
	Double_t zMCVertex = mcHeader->GetVtxZ();
	
	//General settings: vertex, feed down and fill reco container with generated values.  			
	cfVtxHF->SetRecoPrimVertex(zPrimVertex);
	cfVtxHF->SetMCPrimaryVertex(zMCVertex);
	cfVtxHF->SetFillFromGenerated(fFillFromGenerated);
	
	for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
		
		AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
		if (!mcPart){
			AliError("Failed casting particle from MC array!, Skipping particle");
			continue;
		}
		// check the MC-level cuts, must be the desidered particle
		if (!fCFManager->CheckParticleCuts(0, mcPart)) continue;  // 0 stands for MC level
		
		cfVtxHF->SetMCCandidateParam(iPart);
		cfVtxHF->SetNVar(fNvar);
		//counting c quarks
		cquarks += cfVtxHF->MCcquarkCounting(mcPart);
		
		if (!(cfVtxHF->SetLabelArray())){
			AliDebug(2,Form("Impossible to set the label array (decaychannel = %d)",fDecayChannel));
			continue;
		}		   

		//check the candiate family at MC level
		if (!(cfVtxHF->CheckMCPartFamily(mcPart, mcArray))) {
			AliDebug(2,Form("Check on the family wrong!!! (decaychannel = %d)",fDecayChannel));
			continue;
		}
		else{
			AliDebug(2,Form("Check on the family OK!!! (decaychannel = %d)",fDecayChannel));
		}
		
		//Fill the MC container
		Bool_t mcContainerFilled = cfVtxHF -> FillMCContainer(containerInputMC);
		if (mcContainerFilled) {
			if (fUseWeight)fWeight = GetWeight(containerInputMC[0]);
			if (!fCuts->IsInFiducialAcceptance(containerInputMC[0],containerInputMC[1])) continue;
			//MC Limited Acceptance
			if (TMath::Abs(containerInputMC[1]) < 0.5) {
				fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGeneratedLimAcc, fWeight);
				AliDebug(3,"MC Lim Acc container filled\n");
			}	    
			
			//MC 
			fCFManager->GetParticleContainer()->Fill(containerInputMC, kStepGenerated, fWeight);
			icountMC++;
			AliDebug(3,"MC cointainer filled \n");
			
			// MC in acceptance
			// check the MC-Acceptance level cuts
			// since standard CF functions are not applicable, using Kine Cuts on daughters
			Bool_t mcAccepStep = cfVtxHF-> MCAcceptanceStep();
			if (mcAccepStep){	
				fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepAcceptance, fWeight);
				AliDebug(3,"MC acceptance cut passed\n");
				icountAcc++;
				
				//MC Vertex step
				if (fCuts->IsEventSelected(aodEvent)){
					// filling the container if the vertex is ok
					fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepVertex, fWeight) ;
					AliDebug(3,"Vertex cut passed and container filled\n");
					icountVertex++;
					
					//mc Refit requirement	
					Bool_t mcRefitStep = cfVtxHF->MCRefitStep(aodEvent, trackCuts);
					if (mcRefitStep){
						fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepRefit, fWeight);
						AliDebug(3,"MC Refit cut passed and container filled\n");
						icountRefit++;
					}
					else{
						AliDebug(3,"MC Refit cut not passed\n");
						continue;
					}					
				}
				else{
				  AliDebug (3, "MC vertex step not passed\n");
				  continue;
				}
			}
			else{
				AliDebug (3, "MC in acceptance step not passed\n");
				continue;
			}			
		}
		else {
			AliDebug (3, "MC container not filled\n");
		}
	}
	
	if (cquarks<2) AliDebug(2,Form("Event with %d c-quarks", cquarks));
	AliDebug(2,Form("Found %i MC particles that are %s!!",icountMC,fPartName.Data()));
	AliDebug(2,Form("Found %i MC particles that are %s and satisfy Acc cuts!!",icountAcc,fPartName.Data()));
	AliDebug(2,Form("Found %i MC particles that are %s and satisfy Vertex cuts!!",icountVertex,fPartName.Data()));
	AliDebug(2,Form("Found %i MC particles that are %s and satisfy Refit cuts!!",icountRefit,fPartName.Data()));

	// Now go to rec level
	fCountMC += icountMC;
	fCountAcc += icountAcc;
	fCountVertex+= icountVertex;
	fCountRefit+= icountRefit;

	AliDebug(2,Form("Found %d vertices for decay channel %d",arrayBranch->GetEntriesFast(),fDecayChannel));
	
	for(Int_t iCandid = 0; iCandid<arrayBranch->GetEntriesFast();iCandid++){
		AliAODRecoDecayHF* charmCandidate=0x0;
		switch (fDecayChannel){
		case 2:{
			charmCandidate = (AliAODRecoDecayHF2Prong*)arrayBranch->At(iCandid);
			break;
		}
		case 21:{ 
			charmCandidate = (AliAODRecoCascadeHF*)arrayBranch->At(iCandid);
			break;
		}
		case 31:
		case 32:
		case 33:{
			charmCandidate = (AliAODRecoDecayHF3Prong*)arrayBranch->At(iCandid);
			break;
		}
		case 4:{
			charmCandidate = (AliAODRecoDecayHF4Prong*)arrayBranch->At(iCandid);
			break;
		}
		default:
			break;
		}
		
		Bool_t unsetvtx=kFALSE;
		if(!charmCandidate->GetOwnPrimaryVtx()) {
			charmCandidate->SetOwnPrimaryVtx(aodVtx); // needed to compute all variables
			unsetvtx=kTRUE;
		}
		
		Bool_t signAssociation = cfVtxHF->SetRecoCandidateParam((AliAODRecoDecayHF*)charmCandidate);
		if (!signAssociation){
			charmCandidate = 0x0;
			continue;
		}

		Int_t isPartOrAntipart = cfVtxHF->CheckReflexion(fSign);
		if (isPartOrAntipart == 0){
			AliDebug(2, Form("The candidate pdg code doesn't match the requirement set in the task (fSign = %d)",fSign));
			continue;
		}

		Bool_t recoContFilled = cfVtxHF->FillRecoContainer(containerInput);
		if (recoContFilled){
			
			if (!fCuts->IsInFiducialAcceptance(containerInput[0],containerInput[1])) continue;	   
			
			//Reco Step
			Bool_t recoStep = cfVtxHF->RecoStep();
			Bool_t vtxCheck = fCuts->IsEventSelected(aodEvent);
			
			if (recoStep && recoContFilled && vtxCheck){
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed, fWeight) ;   
				icountReco++;
				AliDebug(3,"Reco step  passed and container filled\n");
			  			  
				//Reco in the acceptance -- take care of UNFOLDING!!!!
				Bool_t recoAcceptanceStep = cfVtxHF->RecoAcceptStep(trackCuts);
				if (recoAcceptanceStep) {
					fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoAcceptance, fWeight) ;
					icountRecoAcc++; 
					AliDebug(3,"Reco acceptance cut passed and container filled\n");
				  
					if(fAcceptanceUnf){
						Double_t fill[4]; //fill response matrix
						Bool_t bUnfolding = cfVtxHF -> FillUnfoldingMatrix(fill);
						if (bUnfolding) fCorrelation->Fill(fill);
					}
					
					//Number of ITS cluster requirements	
					Int_t recoITSnCluster = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kTracks);
					if (recoITSnCluster){
						fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoITSClusters, fWeight) ;
						icountRecoITSClusters++;   
						AliDebug(3,"Reco n ITS cluster cut passed and container filled\n");
						
						Bool_t iscutsusingpid = fCuts->GetIsUsePID(); 
						Int_t recoAnalysisCuts = -1, recoPidSelection = -1;
						fCuts->SetUsePID(kFALSE);
						recoAnalysisCuts = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kCandidate, aodEvent);

						if (recoAnalysisCuts > 3){ // Ds case, where more possibilities are considered
							if (recoAnalysisCuts >= 8){
								recoAnalysisCuts -= 8; // removing K0star mass
							}
							if (recoAnalysisCuts >= 4){
								recoAnalysisCuts -= 4; // removing Phi mass
							}
						}

						fCuts->SetUsePID(iscutsusingpid); //restoring usage of the PID from the cuts object	
						if (recoAnalysisCuts == 3 || recoAnalysisCuts == isPartOrAntipart){
							fCFManager->GetParticleContainer()->Fill(containerInput, kStepRecoPPR, fWeight);
							icountRecoPPR++;
							AliDebug(3,"Reco Analysis cuts passed and container filled \n");
							//pid selection
							//recoPidSelection = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kPID);
							//if((fCuts->CombineSelectionLevels(3,recoAnalysisCuts,recoPidSelection)==isPartOrAntipart)||(fCuts->CombineSelectionLevels(3,recoAnalysisCuts,recoPidSelection)==3)){
							recoPidSelection = fCuts->IsSelected(charmCandidate, AliRDHFCuts::kCandidate, aodEvent);
							if (recoPidSelection == 3 || recoPidSelection == isPartOrAntipart){
								fCFManager->GetParticleContainer()->Fill(containerInput, kStepRecoPID, fWeight);
								icountRecoPID++;
								AliDebug(3,"Reco PID cuts passed and container filled \n");
								if(!fAcceptanceUnf){
									Double_t fill[4]; //fill response matrix
									Bool_t bUnfolding = cfVtxHF -> FillUnfoldingMatrix(fill);
									if (bUnfolding) fCorrelation->Fill(fill);
								}
							}
							else {
								AliDebug(3, "Analysis Cuts step not passed \n");
								continue;
							}
						}
						else {
							AliDebug(3, "PID selection not passed \n");
							continue;
						}
					}
					else{
						AliDebug(3, "Number of ITS cluster step not passed\n");
						continue;
					}
				}
				else{
					AliDebug(3, "Reco acceptance step not passed\n");
					continue;
				}
			}
			else {
				AliDebug(3, "Reco step not passed\n");
				continue;
			}
		}
		
		if(unsetvtx) charmCandidate->UnsetOwnPrimaryVtx();
	} // end loop on candidate
		
	fCountReco+= icountReco;
       	fCountRecoAcc+= icountRecoAcc;
	fCountRecoITSClusters+= icountRecoITSClusters;
	fCountRecoPPR+= icountRecoPPR;
	fCountRecoPID+= icountRecoPID;
	
	fHistEventsProcessed->Fill(0);

	delete[] containerInput;
	containerInput = 0x0;
	delete[] containerInputMC;
	containerInputMC = 0x0;
	delete cfVtxHF;

}

//___________________________________________________________________________
void AliCFTaskVertexingHF::Terminate(Option_t*)
{
	// The Terminate() function is the last function to be called during
	// a query. It always runs on the client, it can be used to present
	// the results graphically or save the results to file.
	
	AliAnalysisTaskSE::Terminate();

	AliInfo(Form("Found %i MC particles that are %s in MC, in %d events",fCountMC,fPartName.Data(),fEvents));
	AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, in %d events",fCountAcc,fPartName.Data(),fEvents));
	AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, and satisfy Vertex requirement in %d events",fCountVertex,fPartName.Data(),fEvents));
	AliInfo(Form("Found %i MC particles that are %s in MC and satisfy Acc cuts, and satisfy ITS+TPC refit requirementin %d events",fCountRefit,fPartName.Data(),fEvents));
	AliInfo(Form("Found %i reco %s that are decaying in %s, in %d events",fCountReco,fPartName.Data(),fDauNames.Data(),fEvents));
	AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and are in the requested acceptance, in %d events",fCountRecoAcc,fPartName.Data(),fDauNames.Data(),fEvents));
	AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and have at least 5 clusters in ITS, in %d events",fCountRecoITSClusters,fPartName.Data(),fDauNames.Data(),fEvents));
	AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and satisfy PPR cuts, in %d events",fCountRecoPPR,fPartName.Data(),fDauNames.Data(),fEvents));
	AliInfo(Form("Among the above, found %i reco %s that are decaying in %s and satisfy PPR+PID cuts, in %d events",fCountRecoPID,fPartName.Data(),fDauNames.Data(),fEvents));
	
	// draw some example plots....
	
	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
	if(!cont) {
		printf("CONTAINER NOT FOUND\n");
		return;
	}
	// projecting the containers to obtain histograms
	// first argument = variable, second argument = step
	
	TH1D* h[3][12];
	for(Int_t iC=0;iC<12; iC++){ 
	  // MC-level
	  h[0][iC] =   cont->ShowProjection(iC,0);
	  // MC-Acceptance level
	  h[1][iC] =   cont->ShowProjection(iC,1);
	  // Reco-level
	  h[2][iC] =   cont->ShowProjection(iC,4);
	}	
	TString titles[12];
	if(fDecayChannel==31){
	  titles[0]="pT_Dplus (GeV/c)";
	  titles[1]="rapidity";
	  titles[2]="phi (rad)";
	  titles[3]="cT (#mum)";
	  titles[4]="cosPointingAngle";
	  titles[5]="pT_1 (GeV/c)";
	  titles[6]="pT_2 (GeV/c)";
	  titles[7]="pT_3 (GeV/c)";
	  titles[8]="d0_1 (#mum)";
	  titles[9]="d0_2 (#mum)";
	  titles[10]="d0_3 (#mum)";
	  titles[11]="zVertex (cm)";
	}else{
	  titles[0]="pT_D0 (GeV/c)";
	  titles[1]="rapidity";
	  titles[2]="cosThetaStar";
	  titles[3]="pT_pi (GeV/c)";
	  titles[4]="pT_K (Gev/c)";
	  titles[5]="cT (#mum)";
	  titles[6]="dca (#mum)";
	  titles[7]="d0_pi (#mum)";
	  titles[8]="d0_K (#mum)";
	  titles[9]="d0xd0 (#mum^2)";
	  titles[10]="cosPointingAngle";
	  titles[11]="phi (rad)";
	}
	Int_t markers[12]={20,24,21,25,27,28,
			   20,24,21,25,27,28};
	Int_t colors[3]={2,8,4};
	for(Int_t iC=0;iC<12; iC++){ 
	  for(Int_t iStep=0;iStep<3;iStep++){
	    h[iStep][iC]->SetTitle(titles[iC].Data());
	    h[iStep][iC]->GetXaxis()->SetTitle(titles[iC].Data());
	    Double_t maxh=h[iStep][iC]->GetMaximum();
	    h[iStep][iC]->GetYaxis()->SetRangeUser(0,maxh*1.2);
	    h[iStep][iC]->SetMarkerStyle(markers[iC]);
	    h[iStep][iC]->SetMarkerColor(colors[iStep]);	    
	  }
	}

	
	
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);
	
	// drawing in 2 separate canvas for a matter of clearity
	TCanvas * c1 =new TCanvas("c1New","Vars 0,1,2",1100,1600);
	c1->Divide(3,3);
	Int_t iPad=1;
	for(Int_t iVar=0; iVar<3; iVar++){
	  c1->cd(iPad++);
	  h[0][iVar]->Draw("p");
	  c1->cd(iPad++);
	  h[1][iVar]->Draw("p");
	  c1->cd(iPad++);
	  h[2][iVar]->Draw("p");
	}
	
	TCanvas * c2 =new TCanvas("c2New","Vars 3,4,5",1100,1600);
	c2->Divide(3,3);
	iPad=1;
	for(Int_t iVar=3; iVar<6; iVar++){
	  c2->cd(iPad++);
	  h[0][iVar]->Draw("p");
	  c2->cd(iPad++);
	  h[1][iVar]->Draw("p");
	  c2->cd(iPad++);
	  h[2][iVar]->Draw("p");
	}

	
	TCanvas * c3 =new TCanvas("c3New","Vars 6,7,8",1100,1600);
	c3->Divide(3,3);
	iPad=1;
	for(Int_t iVar=6; iVar<9; iVar++){
	  c3->cd(iPad++);
	  h[0][iVar]->Draw("p");
	  c3->cd(iPad++);
	  h[1][iVar]->Draw("p");
	  c3->cd(iPad++);
	  h[2][iVar]->Draw("p");
	}

	
	TCanvas * c4 =new TCanvas("c4New","Vars 9,10,11",1100,1600);
	c4->Divide(3,3);
	iPad=1;
	for(Int_t iVar=9; iVar<11; iVar++){
	  c4->cd(iPad++);
	  h[0][iVar]->Draw("p");
	  c4->cd(iPad++);
	  h[1][iVar]->Draw("p");
	  c4->cd(iPad++);
	  h[2][iVar]->Draw("p");
	}

	
	THnSparseD* hcorr = dynamic_cast<THnSparseD*> (GetOutputData(3));
	
	TH2D* corr1 =hcorr->Projection(0,2);
	TH2D* corr2 = hcorr->Projection(1,3);
	
	TCanvas * c7 =new TCanvas("c7New","",800,400);
	c7->Divide(2,1);
	c7->cd(1);
	corr1->Draw("text");
	c7->cd(2);
	corr2->Draw("text");
	
	
	TFile* file_projection = new TFile("CFtaskHFprojectionNew.root","RECREATE");
	
	corr1->Write();
	corr2->Write();
	for(Int_t iC=0;iC<12; iC++){ 
	  for(Int_t iStep=0;iStep<3;iStep++){
	    h[iStep][iC]->Write(Form("Step%d_%s",iStep,titles[iC].Data()));
	  }
	}
	file_projection->Close();
	
    
	
	/*
	 c1->SaveAs("Plots/pT_rapidity_cosThetaStar.eps");
	 c2->SaveAs("Plots/pTpi_pTK_cT.eps");
	 c3->SaveAs("Plots/dca_d0pi_d0TK.eps");
	 c4->SaveAs("Plots/d0xd0_cosPointingAngle.eps");
	 
	 c1->SaveAs("Plots/pT_rapidity_cosThetaStar.gif");
	 c2->SaveAs("Plots/pTpi_pTK_cT.gif");
	 c3->SaveAs("Plots/dca_d0pi_d0TK.gif");
	 c4->SaveAs("Plots/d0xd0_cosPointingAngle.gif");
	 */	
}

//___________________________________________________________________________
void AliCFTaskVertexingHF::UserCreateOutputObjects() 
{
	//HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
	//TO BE SET BEFORE THE EXECUTION OF THE TASK
	//
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	OpenFile(1);
	fHistEventsProcessed = new TH1I("CFHFchist0","",1,0,1) ;
}


//_________________________________________________________________________
Double_t AliCFTaskVertexingHF::GetWeight(Float_t pt)
{
	//
	// calculating the weight to fill the container
	//
	
	// FNOLL central:
	// p0 = 1.63297e-01 --> 0.322643
	// p1 = 2.96275e+00
	// p2 = 2.30301e+00
	// p3 = 2.50000e+00

	// PYTHIA
	// p0 = 1.85906e-01 --> 0.36609
	// p1 = 1.94635e+00
	// p2 = 1.40463e+00
	// p3 = 2.50000e+00

	Double_t func1[4] = {0.322643,2.96275,2.30301,2.5};
	Double_t func2[4] = {0.36609,1.94635,1.40463,2.5};

	Double_t dndpt_func1 = dNdptFit(pt,func1);
	Double_t dndpt_func2 = dNdptFit(pt,func2);
	AliDebug(2,Form("pt = %f, FONLL = %f, Pythia = %f, ratio = %f",pt,dndpt_func1,dndpt_func2,dndpt_func1/dndpt_func2));
	return dndpt_func1/dndpt_func2;
}

//__________________________________________________________________________________________________
Double_t AliCFTaskVertexingHF::dNdptFit(Float_t pt, Double_t* par)
{	
	// 
	// calculating dNdpt
	//
	
	Double_t denom =  TMath::Power((pt/par[1]), par[3] );
	Double_t dNdpt = par[0]*pt/TMath::Power(1.+denom, par[2]);
	
	return dNdpt;
}
