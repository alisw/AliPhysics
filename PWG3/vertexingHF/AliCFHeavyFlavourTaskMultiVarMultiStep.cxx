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

/* $Id$ */

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables
// 6 Steps introduced: MC, MC Acc, Reco, Reco Acc, Reco Acc + ITS Cl, 
// Reco Acc + ITS Cl + PPR cuts
// 13 variables used: pt, y, cosThetaStar, ptPi, ptK, ct,
// dca, d0Pi, d0K, d0Pixd0K, cosPointingAngle, phi, z
//
//-----------------------------------------------------------------------
// Author : C. Zampolli, CERN
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

#include "AliCFHeavyFlavourTaskMultiVarMultiStep.h"
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
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliESDtrack.h"
#include "AliRDHFCutsD0toKpi.h"
#include "TChain.h"
#include "THnSparse.h"
#include "TH2D.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"

//__________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep() :
	AliAnalysisTaskSE(),
	fPDG(0),
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
	fFillFromGenerated(kFALSE),
	fMinITSClusters(5),
        fAcceptanceUnf(kTRUE),
	fKeepD0fromB(kFALSE),
	fKeepD0fromBOnly(kFALSE),
	fCuts(0),
	fUseWeight(kFALSE),
	fWeight(1.),
	fSign(2)
{
	//
	//Default ctor
	//
}
//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep(const Char_t* name, AliRDHFCutsD0toKpi* cuts) :
	AliAnalysisTaskSE(name),
	fPDG(0),
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
	fFillFromGenerated(kFALSE),
	fMinITSClusters(5),
        fAcceptanceUnf(kTRUE),
	fKeepD0fromB(kFALSE),
        fKeepD0fromBOnly(kFALSE),
	fCuts(cuts),
	fUseWeight(kFALSE),
	fWeight(1.),
	fSign(2)
{
	//
	// Constructor. Initialization of Inputs and Outputs
	//
	Info("AliCFHeavyFlavourTaskMultiVarMultiStep","Calling Constructor");
	/*
	  DefineInput(0) and DefineOutput(0)
	  are taken care of by AliAnalysisTaskSE constructor
	*/
	DefineOutput(1,TH1I::Class());
	DefineOutput(2,AliCFContainer::Class());
        DefineOutput(3,THnSparseD::Class());
        DefineOutput(4,AliRDHFCutsD0toKpi::Class());

	fCuts->PrintAll();
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep& AliCFHeavyFlavourTaskMultiVarMultiStep::operator=(const AliCFHeavyFlavourTaskMultiVarMultiStep& c) 
{
	//
	// Assignment operator
	//
	if (this!=&c) {
		AliAnalysisTaskSE::operator=(c) ;
		fPDG      = c.fPDG;
		fCFManager  = c.fCFManager;
		fHistEventsProcessed = c.fHistEventsProcessed;
		fCuts = c.fCuts;
	}
	return *this;
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::AliCFHeavyFlavourTaskMultiVarMultiStep(const AliCFHeavyFlavourTaskMultiVarMultiStep& c) :
	AliAnalysisTaskSE(c),
	fPDG(c.fPDG),
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
	fFillFromGenerated(c.fFillFromGenerated),
	fMinITSClusters(c.fMinITSClusters),
        fAcceptanceUnf(c.fAcceptanceUnf),
	fKeepD0fromB(c.fKeepD0fromB),
        fKeepD0fromBOnly(c.fKeepD0fromBOnly),
	fCuts(c.fCuts),
	fUseWeight(c.fUseWeight),
	fWeight(c.fWeight),
	fSign(c.fSign)
{
	//
	// Copy Constructor
	//
}

//___________________________________________________________________________
AliCFHeavyFlavourTaskMultiVarMultiStep::~AliCFHeavyFlavourTaskMultiVarMultiStep() {
	//
	//destructor
	//
	if (fCFManager)           delete fCFManager ;
	if (fHistEventsProcessed) delete fHistEventsProcessed ;
        if (fCorrelation) delete fCorrelation ;
	//        if (fCuts) delete fCuts ;
}

//________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::Init(){

	//
	// Initialization
	//

	if(fDebug > 1) printf("AliCFHeavyFlavourTaskMultiVarMultiStep::Init() \n");
	
	fMinITSClusters = fCuts->GetTrackCuts()->GetMinNClustersITS();
	AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
	const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
	copyfCuts->SetName(nameoutput);
	// Post the data
	PostData(4,copyfCuts);
	
	return;
}

//_________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::UserExec(Option_t *)
{
	//
	// Main loop function
	//
	
	PostData(1,fHistEventsProcessed) ;
	PostData(2,fCFManager->GetParticleContainer()) ;
        PostData(3,fCorrelation) ;

	AliESDtrackCuts* trackCuts = fCuts->GetTrackCuts(); // track cuts

	if (fFillFromGenerated){
		AliWarning("Flag to fill container with generated value ON ---> dca, d0pi, d0K, d0xd0, cosPointingAngle will be set as dummy!");
	}

	if (!fInputEvent) {
		Error("UserExec","NO EVENT FOUND!");
		return;
	}
	
	// check that the fKeepD0fromB flag is set to true when the fKeepD0fromBOnly flag is true
	if(fKeepD0fromBOnly) { 
	  fKeepD0fromB=true;   
	  if(fEvents<2) AliInfo(Form("Both fKeepD0fromB and fKeepD0fromBOnly flags are true, looking _ONLY_ at D0 FROM B"));
	}

	AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

	TClonesArray *arrayD0toKpi=0;

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
	    arrayD0toKpi=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
	  }
	} else {
	  arrayD0toKpi=(TClonesArray*)aodEvent->GetList()->FindObject("D0toKpi");
	}


	if (!arrayD0toKpi) {
	  AliError("Could not find array of HF vertices");
	  return;
	}

	// fix for temporary bug in ESDfilter 
	// the AODs with null vertex pointer didn't pass the PhysSel
	if(!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField())<0.001) return;

	fEvents++;

	fCFManager->SetRecEventInfo(aodEvent);
	fCFManager->SetMCEventInfo(aodEvent);
	
	// MC-event selection
	Double_t containerInput[13] ;
	Double_t containerInputMC[13] ;
        
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
	
	AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
	if (!mcHeader) {
		AliError("Could not find MC Header in AOD");
		return;
	}

	Int_t cquarks = 0;
		
	// AOD primary vertex
	AliAODVertex *vtx1 = (AliAODVertex*)aodEvent->GetPrimaryVertex();
	if(!vtx1) { 
	  AliError("There is no primary vertex !"); 
	  return; 
	}
	Double_t zPrimVertex = vtx1->GetZ();
	Double_t zMCVertex = mcHeader->GetVtxZ();
	Bool_t vtxFlag = kTRUE;
	TString title=vtx1->GetTitle();
	if(!title.Contains("VertexerTracks")) vtxFlag=kFALSE;

	for (Int_t iPart=0; iPart<mcArray->GetEntriesFast(); iPart++) { 
		AliAODMCParticle* mcPart = dynamic_cast<AliAODMCParticle*>(mcArray->At(iPart));
		if (!mcPart) {
			AliWarning("Particle not found in tree, skipping"); 
			continue;
		} 
		if (mcPart->GetPdgCode() == 4) cquarks++; 
		if (mcPart->GetPdgCode() == -4) cquarks++; 
		
		// check the MC-level cuts

		if (!fCFManager->CheckParticleCuts(0, mcPart)) continue;  // 0 stands for MC level
		Int_t pdgGranma = CheckOrigin(mcPart, mcArray);
		Int_t abspdgGranma = TMath::Abs(pdgGranma);
		if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
			AliDebug(2,Form("Particle has a b-meson, or b-baryon mother (pdg code mother = %d )--> not coming from a c-quark, skipping...", pdgGranma));
			if (!fKeepD0fromB) continue;  // skipping particles that don't come from c quark
		}
		else { if(fKeepD0fromBOnly) continue; } // skipping particles that don't come from b quark
		
		//		if (TMath::Abs(pdgGranma)!=4) {

		// fill the container for Gen-level selection
		Double_t vectorMC[7] = {9999.,9999.,9999.,9999.,9999.,9999.,9999.};

		if (GetGeneratedValuesFromMCParticle(mcPart, mcArray, vectorMC)){
			containerInputMC[0] = vectorMC[0];
			containerInputMC[1] = vectorMC[1] ;
			containerInputMC[2] = vectorMC[2] ;
			containerInputMC[3] = vectorMC[3] ;
			containerInputMC[4] = vectorMC[4] ;
			containerInputMC[5] = vectorMC[5] ;  // in micron
			containerInputMC[6] = 0.;    // dummy value, meaningless in MC, in micron
			containerInputMC[7] = 0.;   // dummy value, meaningless in MC, in micron
			containerInputMC[8] = 0.;   // dummy value, meaningless in MC, in micron
			containerInputMC[9] = -100000.; // dummy value, meaningless in MC, in micron^2
			containerInputMC[10] = 1.01;    // dummy value, meaningless in MC
			containerInputMC[11] = vectorMC[6];    // dummy value, meaningless in MC
			containerInputMC[12] = zMCVertex;    // z of reconstructed of primary vertex
			if (fUseWeight) fWeight = GetWeight(vectorMC[0]); // setting the weight according to the function defined in AliCFHeavyFlavourTaskMultiVarMultiStep::GetWeight(Float_t pt)
			AliDebug(3,Form("weight = %f",fWeight));
			if (!fCuts->IsInFiducialAcceptance(vectorMC[0],vectorMC[1])) continue;
			if (TMath::Abs(vectorMC[1]) < 0.5) {
				fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGeneratedLimAcc,fWeight);
			}
			fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepGenerated,fWeight);
			icountMC++;

			// check the MC-Acceptance level cuts
			// since standard CF functions are not applicable, using Kine Cuts on daughters
			
			Int_t daughter0 = mcPart->GetDaughter(0);
			Int_t daughter1 = mcPart->GetDaughter(1);
			AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
			if (daughter0 == 0 || daughter1 == 0) {
				AliDebug(2, "Error! the D0 MC doesn't have correct daughters!! But it should have, this check was already done...");
			}
			if (TMath::Abs(daughter1 - daughter0) != 1) {
				AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, but it should be, this check was already done...");
			}
			AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
			AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
			if (!mcPartDaughter0 || !mcPartDaughter1) {
				AliWarning("At least one Daughter Particle not found in tree, but it should be, this check was already done..."); 
				continue;
			}
			Double_t eta0 = mcPartDaughter0->Eta();
			Double_t eta1 = mcPartDaughter1->Eta();
			Double_t y0 = mcPartDaughter0->Y();
			Double_t y1 = mcPartDaughter1->Y();
			Double_t pt0 = mcPartDaughter0->Pt();
			Double_t pt1 = mcPartDaughter1->Pt();
			AliDebug(2, Form("Daughter 0: eta = %f, y = %f, pt = %f", eta0, y0, pt0));
			AliDebug(2, Form("Daughter 1: eta = %f, y = %f, pt = %f", eta1, y1, pt1));
			Bool_t daught0inAcceptance = (TMath::Abs(eta0) < 0.9 && pt0 > 0.1); 
			Bool_t daught1inAcceptance = (TMath::Abs(eta1) < 0.9 && pt1 > 0.1); 
			if (daught0inAcceptance && daught1inAcceptance) {
				// checking whether the cuts implemented in the CF are equivalent - simply a cross-check
				AliDebug(2, "Daughter particles in acceptance");
				if (!fCFManager->CheckParticleCuts(1, mcPartDaughter0)) {
					AliDebug(2,"Inconsistency with CF cut for daughter 0!");
				} 
				if (!fCFManager->CheckParticleCuts(1, mcPartDaughter1)) {
					AliDebug(2,"Inconsistency with CF cut for daughter 1!");
				} 
				fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepAcceptance,fWeight);
				icountAcc++;

				// check on the vertex
				if (fCuts->IsEventSelected(aodEvent)){
					AliDebug(2,"Vertex cut passed\n");
					// filling the container if the vertex is ok
					fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepVertex,fWeight) ;
					icountVertex++;
					// check on the kTPCrefit and kITSrefit conditions of the daughters
					Bool_t refitFlag = kTRUE;
					if (trackCuts->GetRequireTPCRefit() || trackCuts->GetRequireITSRefit()){
						Int_t foundDaughters = 0;
						for (Int_t iaod =0; iaod<aodEvent->GetNumberOfTracks(); iaod++){
							AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iaod);
							if(track->GetStatus()&AliESDtrack::kITSpureSA) continue;
								if ((track->GetLabel() == daughter0) || (track->GetLabel() == daughter1)) {
								foundDaughters++;
								if (trackCuts->GetRequireTPCRefit()) {
									    if(!(track->GetStatus()&AliESDtrack::kTPCrefit)){
										    refitFlag = kFALSE;
										    break;
									    }
								}
								if (trackCuts->GetRequireITSRefit()) {
									    if(!(track->GetStatus()&AliESDtrack::kITSrefit)){
										    refitFlag = kFALSE;
										    break;
									    }
								}
							}
							if (foundDaughters == 2){  // both daughters have been checked
								break;
							}
						}
						if (foundDaughters != 2) refitFlag = kFALSE;
					}
					if (refitFlag){
						AliDebug(3,"Refit cut passed\n");
						fCFManager->GetParticleContainer()->Fill(containerInputMC,kStepRefit,fWeight);
						icountRefit++;
					}
					else{
						AliDebug(3,"Refit cut not passed\n");
					}
				}
				else{
					AliDebug(3,"Vertex cut not passed\n");
				}			
			}
			else{
				AliDebug(3,"Acceptance cut not passed\n");
			}
		}
		else {
			AliDebug(3,"Problems in filling the container");
			continue;
		}
 	}

	if (cquarks<2) AliDebug(2,Form("Event found with %d c-quarks", cquarks));
	
	AliDebug(2,Form("Found %i MC particles that are D0!!",icountMC));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Acc cuts!!",icountAcc));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Vertex cuts!!",icountVertex));
	AliDebug(2,Form("Found %i MC particles that are D0 and satisfy Refit cuts!!",icountRefit));

	// Now go to rec level
	fCountMC += icountMC;
	fCountAcc += icountAcc;

	AliDebug(2, Form("Found %d vertices",arrayD0toKpi->GetEntriesFast()));

	Int_t pdgDgD0toKpi[2]={321,211};
	Int_t isD0D0bar=1;// 1 for D0, 2 for D0bar, to be used for the PPR and PID selection steps

	for (Int_t iD0toKpi = 0; iD0toKpi<arrayD0toKpi->GetEntriesFast(); iD0toKpi++) {
		
		AliAODRecoDecayHF2Prong* d0tokpi = (AliAODRecoDecayHF2Prong*)arrayD0toKpi->At(iD0toKpi);
		if(!d0tokpi) continue;
		Bool_t unsetvtx=kFALSE;
		if(!d0tokpi->GetOwnPrimaryVtx()) {
		  d0tokpi->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
		  unsetvtx=kTRUE;
		}

		// find associated MC particle
		Int_t mcLabel = d0tokpi->MatchToMC(421,mcArray,2,pdgDgD0toKpi) ;
		if (mcLabel == -1) 
			{
				AliDebug(2,"No MC particle found");
				continue;
			}
		else {
			AliAODMCParticle* mcVtxHF = (AliAODMCParticle*)mcArray->At(mcLabel);
			if (!mcVtxHF) {
				AliWarning("Could not find associated MC in AOD MC tree");
				continue;
			}

			if (mcVtxHF->GetPdgCode() == 421){  // particle is D0
				if (fSign == 1){ // I ask for D0bar only
					AliDebug(2,"particle is D0, I ask for D0bar only");
					continue;
				}
				else{
					isD0D0bar=1;
				}
			}
			else if (mcVtxHF->GetPdgCode()== -421){ // particle is D0bar
				if (fSign == 0){ // I ask for D0 only
					AliDebug(2,"particle is D0bar, I ask for D0 only");
					continue;
				}
				else{
					isD0D0bar=2;
				}
			} 
			else continue;

			// check whether the daughters have kTPCrefit and kITSrefit set
			AliAODTrack *track0 = (AliAODTrack*)d0tokpi->GetDaughter(0);
			AliAODTrack *track1 = (AliAODTrack*)d0tokpi->GetDaughter(1);
			if( !track0 || !track1 ) {
			  AliWarning("Could not find associated MC daughter tracks in AOD MC tree");
			  continue;
			}
			if ((trackCuts->GetRequireTPCRefit() && (!(track0->GetStatus()&AliESDtrack::kTPCrefit) || !(track1->GetStatus()&AliESDtrack::kTPCrefit))) || 
			    (trackCuts->GetRequireITSRefit() && (!(track0->GetStatus()&AliESDtrack::kITSrefit) || !(track1->GetStatus()&AliESDtrack::kITSrefit)))){
				// skipping if at least one daughter does not have kTPCrefit or kITSrefit, if they were required
				continue;
			}

			// check on the vertex -- could be moved outside the loop on the reconstructed D0... 
			if(!fCuts->IsEventSelected(aodEvent)) {
				// skipping cause vertex is not ok
				continue;
			}

                        const Double_t d0tokpiCuts[9] = {0.3,999999.,1.1,0.,0.,999999.,999999.,999999.,0.};
                        Int_t okD0, okD0bar;
			if (!(d0tokpi->SelectD0(&d0tokpiCuts[0],okD0,okD0bar))){
				// skipping candidate
				continue;
			} 

			// check if associated MC v0 passes the cuts
			if (!fCFManager->CheckParticleCuts(0 ,mcVtxHF)) {  // 0 stands for MC
			        AliDebug(2, "Skipping the particles due to cuts");
				continue; 
			}
			Int_t pdgGranma = CheckOrigin(mcVtxHF, mcArray);
			Int_t abspdgGranma = TMath::Abs(pdgGranma);
			if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
				AliDebug(2,Form("At Reco level, from MC info: Particle has a b-meson, or b-baryon mother (pdg code mother = %d )--> not coming from a c-quark, skipping...", pdgGranma));
				if (!fKeepD0fromB) continue;  // skipping particles that don't come from c quark
			}
			else { if(fKeepD0fromBOnly) continue; } // skipping particles that don't come from b quark

			// fill the container...
			Double_t pt = d0tokpi->Pt();
			Double_t rapidity = d0tokpi->YD0();
			Double_t invMass=0.;
			Double_t cosThetaStar = 9999.;
			Double_t pTpi = 0.;
			Double_t pTK = 0.;
			Double_t dca = d0tokpi->GetDCA();
			Double_t d0pi = 0.;
			Double_t d0K = 0.;
			Double_t d0xd0 = d0tokpi->Prodd0d0();
		        Double_t cosPointingAngle = d0tokpi->CosPointingAngle();
			Double_t phi = d0tokpi->Phi();
			Int_t pdgCode = mcVtxHF->GetPdgCode();
			if (pdgCode > 0){
				cosThetaStar = d0tokpi->CosThetaStarD0();
				pTpi = d0tokpi->PtProng(0);
				pTK = d0tokpi->PtProng(1);
				d0pi = d0tokpi->Getd0Prong(0);
				d0K = d0tokpi->Getd0Prong(1);
				invMass=d0tokpi->InvMassD0();
			}
			else {
				cosThetaStar = d0tokpi->CosThetaStarD0bar();
				pTpi = d0tokpi->PtProng(1);
				pTK = d0tokpi->PtProng(0);
				d0pi = d0tokpi->Getd0Prong(1);
				d0K = d0tokpi->Getd0Prong(0);
				invMass=d0tokpi->InvMassD0bar();
			}

			Double_t cT = d0tokpi->CtD0();

			if (!fFillFromGenerated){
				// ...either with reconstructed values....
				containerInput[0] = pt;
				containerInput[1] = rapidity;
				containerInput[2] = cosThetaStar;
				containerInput[3] = pTpi;
				containerInput[4] = pTK;
				containerInput[5] = cT*1.E4;  // in micron
				containerInput[6] = dca*1.E4;  // in micron
				containerInput[7] = d0pi*1.E4;  // in micron
				containerInput[8] = d0K*1.E4;  // in micron
				containerInput[9] = d0xd0*1.E8;  // in micron^2
				containerInput[10] = cosPointingAngle;  // in micron
				containerInput[11] = phi;  
				containerInput[12] = zPrimVertex;    // z of reconstructed of primary vertex
			}
			else {
				// ... or with generated values				
				Double_t vectorMC[7] = {9999.,9999.,9999.,9999.,9999.,9999.,9999.};
				if (GetGeneratedValuesFromMCParticle(mcVtxHF, mcArray, vectorMC)){
					containerInput[0] = vectorMC[0];
					containerInput[1] = vectorMC[1] ;
					containerInput[2] = vectorMC[2] ;
					containerInput[3] = vectorMC[3] ;
					containerInput[4] = vectorMC[4] ;
					containerInput[5] = vectorMC[5] ;  // in micron
					containerInput[6] = 0.;    // dummy value, meaningless in MC, in micron
					containerInput[7] = 0.;   // dummy value, meaningless in MC, in micron
					containerInput[8] = 0.;   // dummy value, meaningless in MC, in micron
					containerInput[9] = 100000.; // dummy value, meaningless in MC, in micron^2
					containerInput[10] = 1.01;    // dummy value, meaningless in MC
					containerInput[11] = vectorMC[6];   
					containerInput[12] = zMCVertex;    // z of reconstructed of primary vertex
				}
				else {
					AliDebug(3,"Problems in filling the container");
					continue;
				}
			}

			if (!fCuts->IsInFiducialAcceptance(containerInput[0],containerInput[1])) continue; // fiducial region

			AliDebug(2, Form("Filling the container with pt = %f, rapidity = %f, cosThetaStar = %f, pTpi = %f, pTK = %f, cT = %f", containerInput[0], containerInput[1], containerInput[2], containerInput[3], containerInput[4], containerInput[5]));
			icountReco++;
			AliDebug(2,Form("%d: filling RECO step\n",iD0toKpi));
			fCFManager->GetParticleContainer()->Fill(containerInput,kStepReconstructed,fWeight) ;   

			// cut in acceptance
			Float_t etaCutMin, etaCutMax, ptCutMin, ptCutMax;
			trackCuts->GetEtaRange(etaCutMin,etaCutMax);
			trackCuts->GetPtRange(ptCutMin,ptCutMax);
			Bool_t acceptanceProng0 = (d0tokpi->EtaProng(0)>etaCutMin && d0tokpi->EtaProng(0)<etaCutMax && d0tokpi->PtProng(0) > ptCutMin && d0tokpi->PtProng(0) < ptCutMax);
			Bool_t acceptanceProng1 = (d0tokpi->EtaProng(1)>etaCutMin && d0tokpi->EtaProng(1)<etaCutMax && d0tokpi->PtProng(1) > ptCutMin && d0tokpi->PtProng(1) < ptCutMax);
			if (acceptanceProng0 && acceptanceProng1) {
				AliDebug(2,"D0 reco daughters in acceptance");
				fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoAcceptance,fWeight) ;
				icountRecoAcc++; 
				
				if(fAcceptanceUnf){
					
					Double_t fill[4]; //fill response matrix
					
					// dimensions 0&1 : pt,eta (Rec)
					
					fill[0] = pt ;
					fill[1] = rapidity;
					
					// dimensions 2&3 : pt,eta (MC)
					
					fill[2] = mcVtxHF->Pt();
					fill[3] = mcVtxHF->Y();
					
					fCorrelation->Fill(fill);
					
				}  
				
				// cut on the min n. of clusters in ITS
				if (fCuts->IsSelected(d0tokpi,AliRDHFCuts::kTracks)){
					fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoITSClusters) ;
					icountRecoITSClusters++;   
					AliDebug(2,Form("pT = %f, dca = %f, cosThetaStar = %f, pTpi = %f, pTK = %f, d0pi = %f, d0K = %f, d0xd0 = %f, cosPointingAngle = %f", pt, dca, cosThetaStar,pTpi, pTK, d0pi*1E4, d0K*1E4, d0xd0*1E8, cosPointingAngle));
		
					// setting the use of the PID cut when applying the selection to FALSE - whatever it was. Keeping track of the original value
					Bool_t iscutusingpid=fCuts->GetIsUsePID();
					Int_t isselcuts=-1,isselpid=-1;
					fCuts->SetUsePID(kFALSE);	
					//Bool_t origFlag = fCuts->GetIsPrimaryWithoutDaughters();
					//fCuts->SetRemoveDaughtersFromPrim(kFALSE);
					isselcuts = fCuts->IsSelected(d0tokpi,AliRDHFCuts::kCandidate,aodEvent);
					//fCuts->SetRemoveDaughtersFromPrim(origFlag);
					fCuts->SetUsePID(iscutusingpid); // restoring usage of the PID from the cuts object
					if (isselcuts == 3 || isselcuts == isD0D0bar){
						AliDebug(2,"Particle passed PPR cuts (actually cuts for D0 analysis!)");
						fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPPR,fWeight) ;   
						icountRecoPPR++;
						
						if(!fAcceptanceUnf){ // unfolding
							
							Double_t fill[4]; //fill response matrix
							
							// dimensions 0&1 : pt,eta (Rec)
							
							fill[0] = pt ;
							fill[1] = rapidity;
							
							// dimensions 2&3 : pt,eta (MC)
							
							fill[2] = mcVtxHF->Pt();
							fill[3] = mcVtxHF->Y();
							
							fCorrelation->Fill(fill);
							
						}

						isselpid = fCuts->IsSelected(d0tokpi,AliRDHFCuts::kPID);
						if((fCuts->CombineSelectionLevels(3,isselcuts,isselpid)==isD0D0bar)||(fCuts->CombineSelectionLevels(3,isselcuts,isselpid)==3)){
							AliDebug(2,"Particle passed PID cuts");
							fCFManager->GetParticleContainer()->Fill(containerInput,kStepRecoPID,fWeight) ;   
							icountRecoPID++;
						}
					}
				}
			}
		}
		if(unsetvtx) d0tokpi->UnsetOwnPrimaryVtx();
	} // end loop on D0->Kpi

	AliDebug(2, Form("Found %i Reco particles that are D0!!",icountReco));

	fCountReco+= icountReco;
	fCountVertex+= icountVertex;
	fCountRefit+= icountRefit;
	fCountRecoAcc+= icountRecoAcc;
	fCountRecoITSClusters+= icountRecoITSClusters;
	fCountRecoPPR+= icountRecoPPR;
	fCountRecoPID+= icountRecoPID;
	
	fHistEventsProcessed->Fill(0);
	/* PostData(0) is taken care of by AliAnalysisTaskSE */
	//PostData(1,fHistEventsProcessed) ;
	//PostData(2,fCFManager->GetParticleContainer()) ;
        //PostData(3,fCorrelation) ;
}


//___________________________________________________________________________
void AliCFHeavyFlavourTaskMultiVarMultiStep::Terminate(Option_t*)
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
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and have at least %d clusters in ITS, in %d events",fCountRecoITSClusters,fMinITSClusters,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and satisfy PPR cuts, in %d events",fCountRecoPPR,fEvents));
	AliInfo(Form("Among the above, found %i reco D0 that are decaying in K+pi and satisfy PID cuts, in %d events",fCountRecoPID,fEvents));
	
	// draw some example plots....
	
	//	AliCFContainer *cont= dynamic_cast<AliCFContainer*> (GetOutputData(2));
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
	TCanvas * c1 =new TCanvas("c1","pT, rapidiy, cosThetaStar",1100,1600);
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

	TCanvas * c2 =new TCanvas("c2","pTpi, pTK, cT",1100,1600);
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

	TCanvas * c3 =new TCanvas("c3","dca, d0pi, d0K",1100,1600);
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

	TCanvas * c4 =new TCanvas("c4","d0xd0, cosPointingAngle, phi",1100,1600);
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

        TCanvas * c7 =new TCanvas("c7","",800,400);
        c7->Divide(2,1);
        c7->cd(1);
        corr1->Draw("text");
        c7->cd(2);
        corr2->Draw("text");
      

	TString projection_filename="CFtaskHFprojection";
	if(fKeepD0fromBOnly) projection_filename+="_KeepD0fromBOnly";
	else if(fKeepD0fromB) projection_filename+="_KeepD0fromB";
	projection_filename+=".root";
	TFile* file_projection = new TFile(projection_filename.Data(),"RECREATE");

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
void AliCFHeavyFlavourTaskMultiVarMultiStep::UserCreateOutputObjects() {
	//HERE ONE CAN CREATE OUTPUT OBJECTS, IN PARTICULAR IF THE OBJECT PARAMETERS DON'T NEED
	//TO BE SET BEFORE THE EXECUTION OF THE TASK
	//
	Info("UserCreateOutputObjects","CreateOutputObjects of task %s\n", GetName());
	
	//slot #1
	OpenFile(1);
	fHistEventsProcessed = new TH1I("CFHFchist0","",1,0,1) ;
}

//___________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::CosThetaStar(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

	//
	// to calculate cos(ThetaStar) for generated particle
	// using the K, since mcPartDaughter0 and mcPartDaughter1 always correspond to K and pi respectively 
	// (see where the function is called)
	//

	Int_t pdgvtx = mcPart->GetPdgCode();
	/*	if (pdgvtx > 0) { // setting as the first daughter always the kaon, to be used to calculate cos(ThetaStar)
		Int_t pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
		Int_t pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
		AliInfo(Form("D0, with pdgprong0 = %d, pdgprong1 = %d",pdgprong0,pdgprong1));
		AliDebug(2,"This is a D0");
		AliAODMCParticle* mcPartdummy = mcPartDaughter0;
		mcPartDaughter0 = mcPartDaughter1;
		mcPartDaughter1 = mcPartdummy;
	} 
	else{
		AliInfo("D0bar");
	}
	*/
	Int_t pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
	Int_t pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
	if (pdgvtx > 0) { // setting as the first daughter always the kaon, to be used to calculate cos(ThetaStar)
		AliDebug(2,"D0");
	}
	else{
		AliDebug(2,"D0bar");
	}
	if (pdgprong0 == 211){
		AliDebug(2,Form("pdgprong0 = %d, pdgprong1 = %d, switching...",pdgprong0,pdgprong1));
		AliAODMCParticle* mcPartdummy = mcPartDaughter0;
		mcPartDaughter0 = mcPartDaughter1;
		mcPartDaughter1 = mcPartdummy;
		pdgprong0 = TMath::Abs(mcPartDaughter0->GetPdgCode());
		pdgprong1 = TMath::Abs(mcPartDaughter1->GetPdgCode());
	} 

	AliDebug(2,Form("After checking, pdgprong0 = %d, pdgprong1 = %d",pdgprong0,pdgprong1));
	Double_t massvtx = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdgvtx))->Mass();
	Double_t massp[2];
	massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
	massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

	Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);

	Double_t px = mcPartDaughter0->Px()+mcPartDaughter1->Px();
	Double_t py = mcPartDaughter0->Py()+mcPartDaughter1->Py();
	Double_t pz = mcPartDaughter0->Pz()+mcPartDaughter1->Pz();
	Double_t p = TMath::Sqrt(px*px+py*py+pz*pz);
	Double_t e =  TMath::Sqrt(massvtx*massvtx+p*p);
	Double_t beta = p/e;
	Double_t gamma = e/massvtx;
	TVector3 mom(mcPartDaughter0->Px(),mcPartDaughter0->Py(),mcPartDaughter0->Pz());
	TVector3 momTot(mcPartDaughter0->Px()+mcPartDaughter1->Px(),mcPartDaughter0->Py()+mcPartDaughter1->Py(),mcPartDaughter0->Pz()+mcPartDaughter1->Pz());

	Double_t qlprong = mom.Dot(momTot)/momTot.Mag();  // analog to AliAODRecoDecay::QlProng(0)
	
	AliDebug(2,Form("pStar = %f, beta = %f, gamma = %f, qlprong = %f, massp[0] = %f", pStar, beta, gamma, qlprong, massp[0]));
	Double_t cts = (qlprong/gamma-beta*TMath::Sqrt(pStar*pStar+massp[0]*massp[0]))/pStar;
	AliDebug(2,Form("cts = %f", cts));
	return cts;
}
//___________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::CT(AliAODMCParticle* mcPart, AliAODMCParticle* mcPartDaughter0, AliAODMCParticle* mcPartDaughter1) const {

	//
	// to calculate cT for generated particle
	//

	Double_t xmcPart[3] = {0,0,0};
	Double_t xdaughter0[3] = {0,0,0};
	Double_t xdaughter1[3] = {0,0,0};
	mcPart->XvYvZv(xmcPart);  // cm
	mcPartDaughter0->XvYvZv(xdaughter0);  // cm
	mcPartDaughter1->XvYvZv(xdaughter1);  //cm
	Double_t prodVtxD0 = TMath::Sqrt(xmcPart[0]*xmcPart[0]+
					 xmcPart[1]*xmcPart[1]+
					 xmcPart[2]*xmcPart[2]);
	Double_t prodVtxDaughter0 = TMath::Sqrt(xdaughter0[0]*xdaughter0[0]+
						xdaughter0[1]*xdaughter0[1]+
						xdaughter0[2]*xdaughter0[2]);
	Double_t prodVtxDaughter1 = TMath::Sqrt(xdaughter1[0]*xdaughter1[0]+
						xdaughter1[1]*xdaughter1[1]+
						xdaughter1[2]*xdaughter1[2]);
 
	AliDebug(2, Form("D0:        x = %f, y = %f, z = %f, production vertex distance = %f (cm), %f (micron)", xmcPart[0], xmcPart[1], xmcPart[2], prodVtxD0, prodVtxD0*1.E4));
	AliDebug(2, Form("Daughter0: x = %f, y = %f, z = %f, production vertex distance = %f (cm) %f (micron)", xdaughter0[0], xdaughter0[1], xdaughter0[2], prodVtxDaughter0, prodVtxDaughter0*1E4));
	AliDebug(2, Form("Daughter1: x = %f, y = %f, z = %f, production vertex distance = %f (cm) %f (micron)", xdaughter1[0], xdaughter1[1], xdaughter1[2], prodVtxDaughter1, prodVtxDaughter1*1.E4));

	Double_t cT0 = TMath::Sqrt((xdaughter0[0]-xmcPart[0])*(xdaughter0[0]-xmcPart[0])+
				     (xdaughter0[1]-xmcPart[1])*(xdaughter0[1]-xmcPart[1])+
				     (xdaughter0[2]-xmcPart[2])*(xdaughter0[2]-xmcPart[2]));

	Double_t cT1 = TMath::Sqrt((xdaughter1[0]-xmcPart[0])*(xdaughter1[0]-xmcPart[0])+
				     (xdaughter1[1]-xmcPart[1])*(xdaughter1[1]-xmcPart[1])+
				     (xdaughter1[2]-xmcPart[2])*(xdaughter1[2]-xmcPart[2]));

	if ((cT0 - cT1)>1E-5) {
		AliWarning(Form("cT from daughter 0 (%f) different from cT from daughter 1 (%f)! Using cT from daughter 0, but PLEASE, CHECK....",cT0,cT1));
	}

	// calculating cT from cT0

	Double_t mass = TDatabasePDG::Instance()->GetParticle(mcPart->GetPdgCode())->Mass(); // mcPart->GetPdgCode() should return +/- 421 for the D0/D0bar
	Double_t p = mcPart-> P();
	Double_t cT = cT0*mass/p;
	AliDebug(2, Form("cT from daughter 0 = %f (micron)", cT0*1E4)); 
	AliDebug(2, Form("cT from daughter 1 = %f (micron)", cT1*1E4)); 
	AliDebug(2, Form("cT (with mass = %f and p = %f) = %f (micron)", mass, p, cT*1E4));
	return cT;
}
//________________________________________________________________________________
Bool_t AliCFHeavyFlavourTaskMultiVarMultiStep::GetGeneratedValuesFromMCParticle(AliAODMCParticle* mcPart, TClonesArray* mcArray, Double_t* vectorMC) const {

	// 
	// collecting all the necessary info (pt, y, cosThetaStar, ptPi, ptKa, cT) from MC particle
	//

	Bool_t isOk = kFALSE;

	// check whether the D0 decays in pi+K
	// to be added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// could use a cut in the CF, but not clear how to define a TDecayChannel
	// implemented for the time being as a cut on the number of daughters - see later when 
	// getting the daughters

	// getting the daughters
	Int_t daughter0 = mcPart->GetDaughter(0);
	Int_t daughter1 = mcPart->GetDaughter(1);
	AliDebug(2, Form("daughter0 = %d and daughter1 = %d",daughter0,daughter1));
	if (daughter0 == 0 || daughter1 == 0) {
		AliDebug(2, "Error! the D0 MC doesn't have correct daughters!!");
		return isOk;  
	}
	if (TMath::Abs(daughter1 - daughter0) != 1) {
		AliDebug(2, "The D0 MC doesn't come from a 2-prong decay, skipping!!");
		return isOk;  
	}
	AliAODMCParticle* mcPartDaughter0 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter0));
	AliAODMCParticle* mcPartDaughter1 = dynamic_cast<AliAODMCParticle*>(mcArray->At(daughter1));
	if (!mcPartDaughter0 || !mcPartDaughter1) {
		AliWarning("At least one Daughter Particle not found in tree, skipping"); 
		return isOk;  
	}
	if (!(TMath::Abs(mcPartDaughter0->GetPdgCode())==321 &&
	      TMath::Abs(mcPartDaughter1->GetPdgCode())==211) && 
	    !(TMath::Abs(mcPartDaughter0->GetPdgCode())==211 &&
	      TMath::Abs(mcPartDaughter1->GetPdgCode())==321)) {
	  AliDebug(2, "The D0 MC doesn't come from a Kpi decay, skipping!!");
	  return isOk;  
	}

	Double_t vtx1[3] = {0,0,0};   // primary vertex		
	Double_t vtx2daughter0[3] = {0,0,0};   // secondary vertex from daughter 0
	Double_t vtx2daughter1[3] = {0,0,0};   // secondary vertex from daughter 1
	mcPart->XvYvZv(vtx1);  // cm
	// getting vertex from daughters
	mcPartDaughter0->XvYvZv(vtx2daughter0);  // cm
	mcPartDaughter1->XvYvZv(vtx2daughter1);  //cm
	if (vtx2daughter0[0] != vtx2daughter1[0] && vtx2daughter0[1] != vtx2daughter1[1] && vtx2daughter0[2] != vtx2daughter1[2]) {
		AliError("Daughters have different secondary vertex, skipping the track");
		return isOk;
	}
	Int_t nprongs = 2;
	Short_t charge = 0;
	// always instantiate the AliAODRecoDecay with the positive daughter first, the negative second
	AliAODMCParticle* positiveDaugh = mcPartDaughter0;
	AliAODMCParticle* negativeDaugh = mcPartDaughter1;
	if (mcPartDaughter0->GetPdgCode()<0 && mcPartDaughter1->GetPdgCode()>0){
		// inverting in case the positive daughter is the second one
		positiveDaugh = mcPartDaughter1;
		negativeDaugh = mcPartDaughter0;
	}
	// getting the momentum from the daughters
	Double_t px[2] = {positiveDaugh->Px(), negativeDaugh->Px()};		
	Double_t py[2] = {positiveDaugh->Py(), negativeDaugh->Py()};		
	Double_t pz[2] = {positiveDaugh->Pz(), negativeDaugh->Pz()};

	Double_t d0[2] = {0.,0.};		

	AliAODRecoDecayHF* decay = new AliAODRecoDecayHF(vtx1,vtx2daughter0,nprongs,charge,px,py,pz,d0);

	Double_t cosThetaStar = 0.;
	Double_t cosThetaStarD0 = 0.;
	Double_t cosThetaStarD0bar = 0.;
	cosThetaStarD0 = decay->CosThetaStar(1,421,211,321);
	cosThetaStarD0bar = decay->CosThetaStar(0,421,321,211);
	if (mcPart->GetPdgCode() == 421){  // D0
		AliDebug(3, Form("D0, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
		cosThetaStar = cosThetaStarD0;
	}
	else if (mcPart->GetPdgCode() == -421){  // D0bar{
		AliDebug(3, Form("D0bar, with pdgprong0 = %d, pdgprong1 = %d",mcPartDaughter0->GetPdgCode(),mcPartDaughter1->GetPdgCode()));
		cosThetaStar = cosThetaStarD0bar;
	}
	else{
		AliWarning("There are problems!! particle was expected to be either a D0 or a D0bar, check...");
		return vectorMC;
	}
	if (cosThetaStar < -1 || cosThetaStar > 1) {
		AliWarning("Invalid value for cosine Theta star, returning");
		return isOk;
	}

	// calculate cos(Theta*) according to the method implemented herein

	Double_t cts = 9999.;
	cts = CosThetaStar(mcPart, mcPartDaughter0, mcPartDaughter1);
	if (cts < -1 || cts > 1) {
		AliWarning("Invalid value for cosine Theta star from AliCFHeavyFlavourTaskMultiVarMultiStep method");
	}
	if (TMath::Abs(cts - cosThetaStar)>0.001) {
		AliError(Form("cosThetaStar automatically calculated different from that manually calculated!!! cosThetaStar = %f, cosThetaStar = %f", cosThetaStar,cts));
	}
	
	Double_t cT = decay->Ct(421);

	// calculate cT from the method implemented herein
	Double_t cT1 = 0.;
	cT1 = CT(mcPart, mcPartDaughter0, mcPartDaughter1);

	if (TMath::Abs(cT1 - cT)>0.001) {
		AliError(Form("cT automatically calculated different from that manually calculated!!! cT = %f, cT1 = %f",cT,cT1));
	}
	
	// get the pT of the daughters
	
	Double_t pTpi = 0.;
	Double_t pTK = 0.;
	
	if (TMath::Abs(mcPartDaughter0->GetPdgCode()) == 211) {
		pTpi = mcPartDaughter0->Pt();
		pTK = mcPartDaughter1->Pt();
	}
	else {
		pTpi = mcPartDaughter1->Pt();
		pTK = mcPartDaughter0->Pt();
	}

	vectorMC[0] = mcPart->Pt();
	vectorMC[1] = mcPart->Y() ;
	vectorMC[2] = cosThetaStar ;
	vectorMC[3] = pTpi ;
	vectorMC[4] = pTK ;
	vectorMC[5] = cT*1.E4 ;  // in micron
	vectorMC[6] = mcPart->Phi() ;  
	isOk = kTRUE;
	return isOk;
}
//_________________________________________________________________________________________________
Int_t AliCFHeavyFlavourTaskMultiVarMultiStep::CheckOrigin(AliAODMCParticle* mcPart, TClonesArray* mcArray)const{

	//
	// checking whether the very mother of the D0 is a charm or a bottom quark
	//

	Int_t pdgGranma = 0;
	Int_t mother = 0;
	mother = mcPart->GetMother();
	Int_t istep = 0;
	while (mother >0 ){
		istep++;
		AliDebug(2,Form("mother at step %d = %d", istep, mother));
		AliAODMCParticle* mcGranma = dynamic_cast<AliAODMCParticle*>(mcArray->At(mother));
		if(!mcGranma) break;
		pdgGranma = mcGranma->GetPdgCode();
		AliDebug(2,Form("Pdg mother at step %d = %d", istep, pdgGranma));
		Int_t abspdgGranma = TMath::Abs(pdgGranma);
		if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
			break;
		}
		mother = mcGranma->GetMother();
	}
	return pdgGranma;
}
//__________________________________________________________________________________________________
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::GetWeight(Float_t pt){

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
Double_t AliCFHeavyFlavourTaskMultiVarMultiStep::dNdptFit(Float_t pt, Double_t* par){

	// 
	// calculating dNdpt
	//

	Double_t denom =  TMath::Power((pt/par[1]), par[3] );
	Double_t dNdpt = par[0]*pt/TMath::Power(1.+denom, par[2]);
	
	return dNdpt;
}
