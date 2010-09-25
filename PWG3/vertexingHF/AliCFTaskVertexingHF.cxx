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
	fNvar(0)
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
	fNvar(0)
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
	fNvar(c.fNvar)
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
	

	switch (fDecayChannel){
	case 2:{
		copyfCuts = new AliRDHFCutsD0toKpi(*(dynamic_cast<AliRDHFCutsD0toKpi*>(fCuts)));
		fNvar = 13;
		break;
	}
	case 21:{ 
		copyfCuts = new AliRDHFCutsDStartoKpipi(*(dynamic_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
		fNvar = 13;
		break;
	}
	case 31:{
		copyfCuts = new AliRDHFCutsDplustoKpipi(*(dynamic_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
		fNvar = 12;
		break;
	}
	case 32:{
		copyfCuts = new AliRDHFCutsLctopKpi(*(dynamic_cast<AliRDHFCutsLctopKpi*>(fCuts)));
		fNvar = 13;
		break;
	}
	case 33:{
		copyfCuts = new AliRDHFCutsDstoKKpi(*(dynamic_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
		fNvar = 13;
		break;
	}
	case 4:{
		copyfCuts = new AliRDHFCutsD0toKpipipi(*(dynamic_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
		fNvar = 13;
		break;
	}
	default:
		AliFatal("The decay channel MUST be defined according to AliCFVertexing::DecayChannel - Exiting...");
		break;
	}  
	
	const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
	copyfCuts->SetName(nameoutput);
	
	//Post the data
	PostData(4, copyfCuts);
	
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
	
	Double_t* containerInput = new Double_t[fNvar];
	Double_t* containerInputMC = new Double_t[fNvar]; 
	
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
	
	Double_t zPrimVertex = aodVtx ->GetZ();
	Double_t zMCVertex = mcHeader->GetVtxZ();
	
	//General settings: vertex, feed down and fill reco container with generated values.  			
	cfVtxHF->SetRecoPrimVertex(zPrimVertex);
	cfVtxHF->SetMCPrimaryVertex(zMCVertex);
	cfVtxHF->SetFillFromGenerated(fFillFromGenerated);
	
	for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
		
		AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
		
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
	AliDebug(2,Form("Found %i MC particles that are D0!!",icountMC));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Acc cuts!!",icountAcc));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Vertex cuts!!",icountVertex));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Refit cuts!!",icountRefit));

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
		
		Int_t isPartOrAntipart = cfVtxHF->CheckReflexion();
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
	
	AliInfo(Form("Found %i MC particles that are D0 in MC, in %d events",fCountMC,fEvents));
	AliInfo(Form("Found %i MC particles that are D0 in MC and satisfy Acc cuts, in %d events",fCountAcc,fEvents));
	AliInfo(Form("Found %i MC particles that are D0 in MC and satisfy Acc cuts, and satisfy Vertex requirement in %d events",fCountVertex,fEvents));
	AliInfo(Form("Found %i MC particles that are D0 in MC and satisfy Acc cuts, and satisfy ITS+TPC refit requirementin %d events",fCountRefit,fEvents));
	AliInfo(Form("Found %i reco D0 that are decaying in K+pi, in %d events",fCountReco,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and are in the requested acceptance, in %d events",fCountRecoAcc,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and have at least 5 clusters in ITS, in %d events",fCountRecoITSClusters,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and satisfy PPR cuts, in %d events",fCountRecoPPR,fEvents));
	
	// draw some example plots....
	
	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
	if(!cont) {
		printf("CONTAINER NOT FOUND\n");
		return;
	}
	// projecting the containers to obtain histograms
	// first argument = variable, second argument = step
	
	// MC-level
	TH1D* h00 =   cont->ShowProjection(0,0) ;   // pt
	TH1D* h10 =   cont->ShowProjection(1,0) ;   // rapidity
	TH1D* h20 =   cont->ShowProjection(2,0) ;   // cosThetaStar
	TH1D* h30 =   cont->ShowProjection(3,0) ;   // pTpi
	TH1D* h40 =   cont->ShowProjection(4,0) ;   // pTK
	TH1D* h50 =   cont->ShowProjection(5,0) ;   // cT
	TH1D* h60 =   cont->ShowProjection(6,0) ;   // dca
	TH1D* h70 =   cont->ShowProjection(7,0) ;   // d0pi
	TH1D* h80 =   cont->ShowProjection(8,0) ;   // d0K
	TH1D* h90 =   cont->ShowProjection(9,0) ;   // d0xd0
	TH1D* h100 =   cont->ShowProjection(10,0) ;   // cosPointingAngle
	TH1D* h110 =   cont->ShowProjection(11,0) ;   // phi
	
	// MC-Acceptance level
	TH1D* h01 =   cont->ShowProjection(0,1) ;   // pt
	TH1D* h11 =   cont->ShowProjection(1,1) ;   // rapidity
	TH1D* h21 =   cont->ShowProjection(2,1) ;   // cosThetaStar
	TH1D* h31 =   cont->ShowProjection(3,1) ;   // pTpi
	TH1D* h41 =   cont->ShowProjection(4,1) ;   // pTK
	TH1D* h51 =   cont->ShowProjection(5,1) ;   // cT
	TH1D* h61 =   cont->ShowProjection(6,1) ;   // dca
	TH1D* h71 =   cont->ShowProjection(7,1) ;   // d0pi
	TH1D* h81 =   cont->ShowProjection(8,1) ;   // d0K
	TH1D* h91 =   cont->ShowProjection(9,1) ;   // d0xd0
	TH1D* h101 =   cont->ShowProjection(10,1) ;   // cosPointingAngle
	TH1D* h111 =   cont->ShowProjection(11,1) ;   // phi
	
	// Reco-level
	TH1D* h04 =   cont->ShowProjection(0,4) ;   // pt
	TH1D* h14 =   cont->ShowProjection(1,4) ;   // rapidity
	TH1D* h24 =   cont->ShowProjection(2,4) ;   // cosThetaStar
	TH1D* h34 =   cont->ShowProjection(3,4) ;   // pTpi
	TH1D* h44 =   cont->ShowProjection(4,4) ;   // pTK
	TH1D* h54 =   cont->ShowProjection(5,4) ;   // cT
	TH1D* h64 =   cont->ShowProjection(6,4) ;   // dca
	TH1D* h74 =   cont->ShowProjection(7,4) ;   // d0pi
	TH1D* h84 =   cont->ShowProjection(8,4) ;   // d0K
	TH1D* h94 =   cont->ShowProjection(9,4) ;   // d0xd0
	TH1D* h104 =   cont->ShowProjection(10,4) ;   // cosPointingAngle
	TH1D* h114 =   cont->ShowProjection(11,4) ;   // phi
	
	h00->SetTitle("pT_D0 (GeV/c)");
	h10->SetTitle("rapidity");
	h20->SetTitle("cosThetaStar");
	h30->SetTitle("pT_pi (GeV/c)");
	h40->SetTitle("pT_K (Gev/c)");
	h50->SetTitle("cT (#mum)");
	h60->SetTitle("dca (#mum)");
	h70->SetTitle("d0_pi (#mum)");
	h80->SetTitle("d0_K (#mum)");
	h90->SetTitle("d0xd0 (#mum^2)");
	h100->SetTitle("cosPointingAngle");
	h100->SetTitle("phi (rad)");
	
	h00->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h10->GetXaxis()->SetTitle("rapidity");
	h20->GetXaxis()->SetTitle("cosThetaStar");
	h30->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h40->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h50->GetXaxis()->SetTitle("cT (#mum)");
	h60->GetXaxis()->SetTitle("dca (#mum)");
	h70->GetXaxis()->SetTitle("d0_pi (#mum)");
	h80->GetXaxis()->SetTitle("d0_K (#mum)");
	h90->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h100->GetXaxis()->SetTitle("cosPointingAngle");
	h110->GetXaxis()->SetTitle("phi (rad)");
	
	h01->SetTitle("pT_D0 (GeV/c)");
	h11->SetTitle("rapidity");
	h21->SetTitle("cosThetaStar");
	h31->SetTitle("pT_pi (GeV/c)");
	h41->SetTitle("pT_K (Gev/c)");
	h51->SetTitle("cT (#mum)");
	h61->SetTitle("dca (#mum)");
	h71->SetTitle("d0_pi (#mum)");
	h81->SetTitle("d0_K (#mum)");
	h91->SetTitle("d0xd0 (#mum^2)");
	h101->SetTitle("cosPointingAngle");
	h111->GetXaxis()->SetTitle("phi (rad)");
	
	h01->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h11->GetXaxis()->SetTitle("rapidity");
	h21->GetXaxis()->SetTitle("cosThetaStar");
	h31->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h41->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h51->GetXaxis()->SetTitle("cT (#mum)");
	h61->GetXaxis()->SetTitle("dca (#mum)");
	h71->GetXaxis()->SetTitle("d0_pi (#mum)");
	h81->GetXaxis()->SetTitle("d0_K (#mum)");
	h91->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h101->GetXaxis()->SetTitle("cosPointingAngle");
	h111->GetXaxis()->SetTitle("phi (rad)");
	
	h04->SetTitle("pT_D0 (GeV/c)");
	h14->SetTitle("rapidity");
	h24->SetTitle("cosThetaStar");
	h34->SetTitle("pT_pi (GeV/c)");
	h44->SetTitle("pT_K (Gev/c)");
	h54->SetTitle("cT (#mum)");
	h64->SetTitle("dca (#mum)");
	h74->SetTitle("d0_pi (#mum)");
	h84->SetTitle("d0_K (#mum)");
	h94->SetTitle("d0xd0 (#mum^2)");
	h104->SetTitle("cosPointingAngle");
	h114->GetXaxis()->SetTitle("phi (rad)");
	
	h04->GetXaxis()->SetTitle("pT_D0 (GeV/c)");
	h14->GetXaxis()->SetTitle("rapidity");
	h24->GetXaxis()->SetTitle("cosThetaStar");
	h34->GetXaxis()->SetTitle("pT_pi (GeV/c)");
	h44->GetXaxis()->SetTitle("pT_K (Gev/c)");
	h54->GetXaxis()->SetTitle("cT (#mum)");
	h64->GetXaxis()->SetTitle("dca (#mum)");
	h74->GetXaxis()->SetTitle("d0_pi (#mum)");
	h84->GetXaxis()->SetTitle("d0_K (#mum)");
	h94->GetXaxis()->SetTitle("d0xd0 (#mum^2)");
	h104->GetXaxis()->SetTitle("cosPointingAngle");
	h114->GetXaxis()->SetTitle("phi (rad)");
	
	Double_t max0 = h00->GetMaximum();
	Double_t max1 = h10->GetMaximum();
	Double_t max2 = h20->GetMaximum();
	Double_t max3 = h30->GetMaximum();
	Double_t max4 = h40->GetMaximum();
	Double_t max5 = h50->GetMaximum();
	Double_t max6 = h60->GetMaximum();
	Double_t max7 = h70->GetMaximum();
	Double_t max8 = h80->GetMaximum();
	Double_t max9 = h90->GetMaximum();
	Double_t max10 = h100->GetMaximum();
	Double_t max11 = h110->GetMaximum();
	
	h00->GetYaxis()->SetRangeUser(0,max0*1.2);
	h10->GetYaxis()->SetRangeUser(0,max1*1.2);
	h20->GetYaxis()->SetRangeUser(0,max2*1.2);
	h30->GetYaxis()->SetRangeUser(0,max3*1.2);
	h40->GetYaxis()->SetRangeUser(0,max4*1.2);
	h50->GetYaxis()->SetRangeUser(0,max5*1.2);
	h60->GetYaxis()->SetRangeUser(0,max6*1.2);
	h70->GetYaxis()->SetRangeUser(0,max7*1.2);
	h80->GetYaxis()->SetRangeUser(0,max8*1.2);
	h90->GetYaxis()->SetRangeUser(0,max9*1.2);
	h100->GetYaxis()->SetRangeUser(0,max10*1.2);
	h110->GetYaxis()->SetRangeUser(0,max11*1.2);
	
	h01->GetYaxis()->SetRangeUser(0,max0*1.2);
	h11->GetYaxis()->SetRangeUser(0,max1*1.2);
	h21->GetYaxis()->SetRangeUser(0,max2*1.2);
	h31->GetYaxis()->SetRangeUser(0,max3*1.2);
	h41->GetYaxis()->SetRangeUser(0,max4*1.2);
	h51->GetYaxis()->SetRangeUser(0,max5*1.2);
	h61->GetYaxis()->SetRangeUser(0,max6*1.2);
	h71->GetYaxis()->SetRangeUser(0,max7*1.2);
	h81->GetYaxis()->SetRangeUser(0,max8*1.2);
	h91->GetYaxis()->SetRangeUser(0,max9*1.2);
	h101->GetYaxis()->SetRangeUser(0,max10*1.2);
	h111->GetYaxis()->SetRangeUser(0,max11*1.2);
	
	h00->SetMarkerStyle(20);
	h10->SetMarkerStyle(24);
	h20->SetMarkerStyle(21);
	h30->SetMarkerStyle(25);
	h40->SetMarkerStyle(27);
	h50->SetMarkerStyle(28);
	h60->SetMarkerStyle(20);
	h70->SetMarkerStyle(24);
	h80->SetMarkerStyle(21);
	h90->SetMarkerStyle(25);
	h100->SetMarkerStyle(27);
	h110->SetMarkerStyle(28);
	
	h00->SetMarkerColor(2);
	h10->SetMarkerColor(2);
	h20->SetMarkerColor(2);
	h30->SetMarkerColor(2);
	h40->SetMarkerColor(2);
	h50->SetMarkerColor(2);
 	h60->SetMarkerColor(2);
	h70->SetMarkerColor(2);
	h80->SetMarkerColor(2);
	h90->SetMarkerColor(2);
	h100->SetMarkerColor(2);
	h110->SetMarkerColor(2);
	
	h01->SetMarkerStyle(20) ;
	h11->SetMarkerStyle(24) ;
	h21->SetMarkerStyle(21) ;
	h31->SetMarkerStyle(25) ;
	h41->SetMarkerStyle(27) ;
	h51->SetMarkerStyle(28) ;
	h61->SetMarkerStyle(20);
	h71->SetMarkerStyle(24);
	h81->SetMarkerStyle(21);
	h91->SetMarkerStyle(25);
	h101->SetMarkerStyle(27);
	h111->SetMarkerStyle(28);
	
	h01->SetMarkerColor(8);
	h11->SetMarkerColor(8);
	h21->SetMarkerColor(8);
	h31->SetMarkerColor(8);
	h41->SetMarkerColor(8);
	h51->SetMarkerColor(8);
 	h61->SetMarkerColor(8);
	h71->SetMarkerColor(8);
	h81->SetMarkerColor(8);
	h91->SetMarkerColor(8);
	h101->SetMarkerColor(8);
	h111->SetMarkerColor(8);
	
	h04->SetMarkerStyle(20) ;
	h14->SetMarkerStyle(24) ;
	h24->SetMarkerStyle(21) ;
	h34->SetMarkerStyle(25) ;
	h44->SetMarkerStyle(27) ;
	h54->SetMarkerStyle(28) ;
	h64->SetMarkerStyle(20);
	h74->SetMarkerStyle(24);
	h84->SetMarkerStyle(21);
	h94->SetMarkerStyle(25);
	h104->SetMarkerStyle(27);
	h114->SetMarkerStyle(28);
	
	h04->SetMarkerColor(4);
	h14->SetMarkerColor(4);
	h24->SetMarkerColor(4);
	h34->SetMarkerColor(4);
	h44->SetMarkerColor(4);
	h54->SetMarkerColor(4);
	h64->SetMarkerColor(4);
	h74->SetMarkerColor(4);
	h84->SetMarkerColor(4);
	h94->SetMarkerColor(4);
	h104->SetMarkerColor(4);
	h114->SetMarkerColor(4);
	
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetStatColor(0);
	
	// drawing in 2 separate canvas for a matter of clearity
	TCanvas * c1 =new TCanvas("c1New","pT, rapidiy, cosThetaStar",1100,1600);
	c1->Divide(3,3);
	
	c1->cd(1);
	h00->Draw("p");
	c1->cd(1);
	c1->cd(2);
	h01->Draw("p");
	c1->cd(2);
	c1->cd(3);
	h04->Draw("p");
	c1->cd(3);
	c1->cd(4);
	h10->Draw("p");
	c1->cd(4);
	c1->cd(5);
	h11->Draw("p");
	c1->cd(5);
	c1->cd(6);
	h14->Draw("p");
	c1->cd(6);
	c1->cd(7);
	h20->Draw("p");
	c1->cd(7);
	c1->cd(8);
	h21->Draw("p");
	c1->cd(8);
	c1->cd(9);
	h24->Draw("p");
	c1->cd(9);
	c1->cd();
	
	TCanvas * c2 =new TCanvas("c2New","pTpi, pTK, cT",1100,1600);
	c2->Divide(3,3);
	c2->cd(1);
	h30->Draw("p");
	c2->cd(1);
	c2->cd(2);
	h31->Draw("p");
	c2->cd(2);
	c2->cd(3);
	h34->Draw("p");
	c2->cd(3);
	c2->cd(4);
	h40->Draw("p");
	c2->cd(4);
	c2->cd(5);
	h41->Draw("p");
	c2->cd(5);
	c2->cd(6);
	h44->Draw("p");
	c2->cd(6);
	c2->cd(7);
	h50->Draw("p");
	c2->cd(7);
	c2->cd(8);
	h51->Draw("p");
	c2->cd(8);
	c2->cd(9);
	h54->Draw("p");
	c2->cd(9);
	c2->cd();
	
	TCanvas * c3 =new TCanvas("c3New","dca, d0pi, d0K",1100,1600);
	c3->Divide(3,3);
	c3->cd(1);
	h60->Draw("p");
	c3->cd(1);
	c3->cd(2);
	h61->Draw("p");
	c3->cd(2);
	c3->cd(3);
	h64->Draw("p");
	c3->cd(3);
	c3->cd(4);
	h70->Draw("p");
	c3->cd(4);
	c3->cd(5);
	h71->Draw("p");
	c3->cd(5);
	c3->cd(6);
	h74->Draw("p");
	c3->cd(6);
	c3->cd(7);
	h80->Draw("p");
	c3->cd(7);
	c3->cd(8);
	h81->Draw("p");
	c3->cd(8);
	c3->cd(9);
	h84->Draw("p");
	c3->cd(9);
	c3->cd();
	
	TCanvas * c4 =new TCanvas("c4New","d0xd0, cosPointingAngle, phi",1100,1600);
	c4->Divide(3,3);
	c4->cd(1);
	h90->Draw("p");
	c4->cd(1);
	c4->cd(2);
	h91->Draw("p");
	c4->cd(2);
	c4->cd(3);
	h94->Draw("p");
	c4->cd(3);
	c4->cd(4);
	h100->Draw("p");
	c4->cd(4);
	c4->cd(5);
	h101->Draw("p");
	c4->cd(5);
	c4->cd(6);
	h104->Draw("p");
	c4->cd(6);
	c4->cd(7);
	h110->Draw("p");
	c4->cd(7);
	c4->cd(8);
	h111->Draw("p");
	c4->cd(8);
	c4->cd(9);
	h114->Draw("p");
	c4->cd(9);
	c4->cd();
	
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
	h00->Write("pT_D0_step0");
	h10->Write("rapidity_step0");
	h20->Write("cosThetaStar_step0");
	h30->Write("pT_pi_step0");
	h40->Write("pT_K_step0");
	h50->Write("cT_step0");
	h60->Write("dca_step0");
	h70->Write("d0_pi_step0");
	h80->Write("d0_K_step0");
	h90->Write("d0xd0_step0");
	h100->Write("cosPointingAngle_step0");
	h110->Write("phi_step0");
	
	h01->Write("pT_D0_step1");
	h11->Write("rapidity_step1");
	h21->Write("cosThetaStar_step1");
	h31->Write("pT_pi_step1");
	h41->Write("pT_K_step1");
	h51->Write("cT_step1");
	h61->Write("dca_step1");
	h71->Write("d0_pi_step1");
	h81->Write("d0_K_step1");
	h91->Write("d0xd0_step1");
	h101->Write("cosPointingAngle_step1");
	h111->Write("phi_step1");
	
	h04->Write("pT_D0_step2");
	h14->Write("rapidity_step2");
	h24->Write("cosThetaStar_step2");
	h34->Write("pT_pi_step2");
	h44->Write("pT_K_step2");
	h54->Write("cT_step2");
	h64->Write("dca_step2");
	h74->Write("d0_pi_step2");
	h80->Write("d0_K_step2");
	h94->Write("d0xd0_step2");
	h104->Write("cosPointingAngle_step2");
	h114->Write("phi_step2");
	
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
