/*************** Tree Maker Class **********************
*                                                      *
* Created: 05.10.2016                                  *
* Authors: Aaron Capon      (aaron.capon@cern.ch)      *
*          Sebastian Lehner (sebastian.lehner@cern.ch) *
*                                                      *
*******************************************************/
/**********************************************************************************
*This analysis task is designed to create simple flat structured TTrees which     *
*can then be worked on locally. The motivation for this was to create TTrees that *
*could be easily used within TMVA.                                                *
*                                                                                 *
*The default track cuts are relatively loose, and the PID selects out electrons   *
*within +-4 nSigma in the TPC. This value, as well as cuts on the other detector  *
*response values, can be turned on and off via their relavent setter function.    *
*                                                                                 *
*The GRID PID number for each job is stored with each event along with a simple   *
*event counter ID, so that after merging each event still has a unique label.     *
*                                                                                 *
*By default, the class will return a TTree focused on selecting electrons from    *
*primary decays. There are setter functions which can be used to instead create a *
*TTree filled with tracks originating from V0 decays. In this case the DCA cuts   *
*are removed.  There is a futher option to create a TTree which contains all      *
*generated particles, and, if reconstructed, their corresponding reconstructed    *
*track features.                                                                  *
*TODO: Speed up the creation of the generator TTree (very slow with multiple      *
*loops over the same particles).                                                  *
**********************************************************************************/
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"

#include "AliVEvent.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

#include "AliAODMCParticle.h"

#include "TParticle.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenCocktailEventHeader.h"

#include "AliVVertex.h"

#include "AliESDv0.h"
#include "AliDielectronTrackCuts.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronV0Cuts.h"

#include "TClonesArray.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"

#include "TSystem.h"
#include <string>
#include "AliAnalysisTaskSimpleTreeMaker.h"



ClassImp(AliAnalysisTaskSimpleTreeMaker)

Int_t eventNum = 0;

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker():
    AliAnalysisTaskSE(),
		hasMC(kFALSE),
    fESDtrackCuts(0),
    fPIDResponse(0),
    fTree(0x0),
		primaryVertex{0,0,0},
		multiplicityV0A(0),
		multiplicityV0C(0),
		multiplicityCL1(0),
		runNumber(0),
		event(0),
		pt(0),
		eta(0),
		phi(0),
		nTPCclusters(0),
		nTPCcrossed(0),
		fTPCcrossOverFind(0),
		nTPCfindable(0),
		tpcSharedMap(0),
		nTPCshared(0),
		chi2TPC(0),
		DCA{0,0},
		nITS(0),
		chi2ITS(0),
		fITSshared(0),
		SPDfirst(0),
		charge(0),
		EnSigmaITS(0),
		EnSigmaTPC(0),
		EnSigmaTOF(0),
		PnSigmaTPC(0),
		PnSigmaITS(0),
		PnSigmaTOF(0),
		KnSigmaITS(0),
		KnSigmaTPC(0),
		KnSigmaTOF(0),
		ITSsignal(0),
		TPCsignal(0),
		TOFsignal(0),
		goldenChi2(0),
		mcEta(0),
		mcPhi(0),
		mcPt(0),
		mcVert{0,0,0},
		iPdg(0),
		iPdgMother(0),
		HasMother(0),
		motherLabel(0),
		pointingAngle(0),
		daughtersDCA(0),
		decayLength(0),
		v0mass(0),
		ptArm(0),
		alpha(0),
    fQAhist(0),
    fCentralityPercentileMin(0),
    fCentralityPercentileMax(80), 
    fPtMin(0.2),
    fPtMax(10),
    fEtaMin(-0.8),
    fEtaMax(0.8),
    fESigITSMin(-3.),
    fESigITSMax(3.),
    fESigTPCMin(-4.),
    fESigTPCMax(4.),
    fESigTOFMin(-3.),
    fESigTOFMax(3.),
    fPIDcutITS(kFALSE),
    fPIDcutTOF(kFALSE),
    fPionPIDcutTPC(kFALSE),
    fPSigTPCMin(-99.),
    fPSigTPCMax(-3.),
    fHasSDD(kTRUE),
    fIsV0tree(kFALSE),
    fArmPlot(0),
		fIsEffTree(kFALSE),
    fIsAOD(kTRUE),
    fFilterBit(16),
    fIsGRIDanalysis(kTRUE),
    fGridPID(-1)
{

}

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker(const char *name) :
    AliAnalysisTaskSE(name),
		hasMC(kFALSE),
   	fESDtrackCuts(0),
    fPIDResponse(0),
    fTree(0),
		primaryVertex{0,0,0},
		multiplicityV0A(0),
		multiplicityV0C(0),
		multiplicityCL1(0),
		runNumber(0),
		event(0),
		pt(0),
		eta(0),
		phi(0),
		nTPCclusters(0),
		nTPCcrossed(0),
		fTPCcrossOverFind(0),
		nTPCfindable(0),
		tpcSharedMap(0),
		nTPCshared(0),
		chi2TPC(0),
		DCA{0,0},
		nITS(0),
		chi2ITS(0),
		fITSshared(0),
		SPDfirst(0),
		charge(0),
		EnSigmaITS(0),
		EnSigmaTPC(0),
		EnSigmaTOF(0),
		PnSigmaTPC(0),
		PnSigmaITS(0),
		PnSigmaTOF(0),
		KnSigmaITS(0),
		KnSigmaTPC(0),
		KnSigmaTOF(0),
		ITSsignal(0),
		TPCsignal(0),
		TOFsignal(0),
		goldenChi2(0),
		mcEta(0),
		mcPhi(0),
		mcPt(0),
		mcVert{0,0,0},
		iPdg(0),
		iPdgMother(0),
		HasMother(0),
		motherLabel(0),
		pointingAngle(0),
		daughtersDCA(0),
		decayLength(0),
		v0mass(0),
		ptArm(0),
		alpha(0),
    fQAhist(0),
    fCentralityPercentileMin(0),
    fCentralityPercentileMax(80), 
    fPtMin(0.2),
    fPtMax(10),
    fEtaMin(-0.8),
    fEtaMax(0.8),
    fESigITSMin(-3.),
    fESigITSMax(3.),
    fESigTPCMin(-4.),
    fESigTPCMax(4.),
    fESigTOFMin(-3.),
    fESigTOFMax(3.),
    fPIDcutITS(kFALSE),
    fPIDcutTOF(kFALSE),
    fPionPIDcutTPC(kFALSE),
		fPSigTPCMin(-99.),
    fPSigTPCMax(-3.),
    fHasSDD(kTRUE),
    fIsV0tree(kFALSE),
    fArmPlot(0),
		fIsEffTree(kFALSE),
    fIsAOD(kTRUE),
    fFilterBit(16),
    fIsGRIDanalysis(kTRUE),
    fGridPID(-1)

{
    if(!fIsV0tree && !fIsEffTree){
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
    }
    else if(fIsV0tree && !fIsEffTree){
        fESDtrackCuts = AliESDtrackCuts::GetStandardV0DaughterCuts();
    }

    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    DefineOutput(1, TTree::Class()); //will be connected to fTree
    DefineOutput(2, TH1F::Class());
    DefineOutput(3, TH2F::Class());
}

//________________________________________________________________________

//~ AliAnalysisTaskSimpleTreeMaker::~AliAnalysisTaskSimpleTreeMaker() {

  //~ // Destructor

  //~ // ... not implemented

//~ }


//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserCreateOutputObjects(){
  
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    inputHandler->SetNeedField();
     
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse){
        AliFatal("This task needs the PID response attached to the inputHandler");
      return;
    } 

    fTree = new TTree("tracks", "tracks");
		//Common branches to all variants of class

		//Track variables
		fTree->Branch("pt",                &pt);
		fTree->Branch("eta",               &eta);
		fTree->Branch("phi",               &phi);
		fTree->Branch("nTPCclusters",      &nTPCclusters);
		fTree->Branch("nTPCcrossed",       &nTPCcrossed);
		fTree->Branch("nTPCfindable",      &nTPCfindable);
		fTree->Branch("nTPCshared",        &nTPCshared);
		fTree->Branch("fTPCcrossOverFind", &fTPCcrossOverFind);
		fTree->Branch("chi2TPC",           &chi2TPC);
		fTree->Branch("nITS",              &nITS);
		fTree->Branch("fITSshared",        &fITSshared);
		fTree->Branch("chi2ITS",           &chi2ITS);
		fTree->Branch("SPDfirst",          &SPDfirst);
		fTree->Branch("DCAxy",             &DCA[0]);
		fTree->Branch("DCAz",              &DCA[1]);
		fTree->Branch("goldenChi2",        &goldenChi2);
		fTree->Branch("charge",            &charge);
		fTree->Branch("EsigITS",           &EnSigmaITS);
		fTree->Branch("EsigTPC",           &EnSigmaTPC);
		fTree->Branch("EsigTOF",           &EnSigmaTOF);
		fTree->Branch("PsigITS",           &PnSigmaITS);
		fTree->Branch("PsigTPC",           &PnSigmaTPC);
		fTree->Branch("PsigTOF",           &PnSigmaTOF);
		fTree->Branch("KsigITS",           &KnSigmaITS);
		fTree->Branch("KsigTPC",           &KnSigmaTPC);
		fTree->Branch("KsigTOF",           &KnSigmaTOF);
		fTree->Branch("ITSsignal",         &ITSsignal);
		fTree->Branch("TPCsignal",         &TPCsignal);
		fTree->Branch("TOFsignal",         &TOFsignal);
		//V0 tree variables
		if(fIsV0tree){
				fTree->Branch("pointingAngle", &pointingAngle);
				fTree->Branch("daughtersDCA",  &daughtersDCA);
				fTree->Branch("decayLength",   &decayLength);
				fTree->Branch("v0mass",        &v0mass);
				fTree->Branch("ptArm",         &ptArm);
				fTree->Branch("alpha",         &alpha);
		}
		//MC variables
		if(hasMC){
			fTree->Branch("mcPt",        &mcPt);
			fTree->Branch("mcEta",       &mcEta);
			fTree->Branch("mcPhi",       &mcPhi);
			fTree->Branch("mcVertx",     &mcVert[0]);
			fTree->Branch("mcVerty",     &mcVert[1]);
			fTree->Branch("mcVertz",     &mcVert[2]);
			fTree->Branch("pdg",         &iPdg);
			fTree->Branch("pdgMother",   &iPdgMother);
			fTree->Branch("HasMother",   &HasMother);
			fTree->Branch("motherLabel", &motherLabel);
		}
		//Event variables
		fTree->Branch("vertexX",         &primaryVertex[0]);
		fTree->Branch("vertexY",         &primaryVertex[1]);
		fTree->Branch("vertexZ",         &primaryVertex[2]);
		fTree->Branch("runNumber",       &runNumber);
		fTree->Branch("eventNum",        &eventNum);
		fTree->Branch("multiplicityV0A", &multiplicityV0A);
		fTree->Branch("multiplicityV0C", &multiplicityV0C);
		fTree->Branch("multiplicityCL1", &multiplicityCL1);
		fTree->Branch("gridPID",         &fGridPID);

    //Get grid PID which can be used later to assign unique event numbers
    if(fIsGRIDanalysis){
        const char* gridIDchar = gSystem->Getenv("ALIEN_PROC_ID");
        std::string str(gridIDchar);
        SetGridPID(str);
    }
    else{ 
        fGridPID = -1;
    }


    //Create TH2F for armenteros plot. Filled if creating v0 tree
		//NOTE: Class designed to study electrons with weak EsigTPC cuts
		//applied. Therefore, the number of tracks will not necessarily be twice the
		//number of V0 particles (as one might expect). 
    fArmPlot = new TH2F("ArmPlot", "Armenteros Plot", 100, -1, 1, 100, 0, 0.4);
    if(fIsV0tree){
        fArmPlot->GetXaxis()->SetTitle("#alpha = (p^{+}-p^{-})/(p^{+}+p^{-})");
        fArmPlot->GetYaxis()->SetTitle("p_{T}");
    }  

    fQAhist = new TH1F("h1", "Event and track QA", 8, 0, 8);
		//Fill each bin with nothing to ensure correct ordering
		fQAhist->Fill("Events_check",0);
		fQAhist->Fill("Events_accepted",0);
		fQAhist->Fill("Events_MCcheck",0);
		fQAhist->Fill("Tracks_all",0);
		fQAhist->Fill("Tracks_MCcheck",0);
		fQAhist->Fill("Tracks_KineCuts",0);
		fQAhist->Fill("Tracks_TrackCuts",0);
		fQAhist->Fill("Tracks_PIDcuts",0);
		    
		PostData(1, fTree);
    PostData(2, fQAhist);
    PostData(3, fArmPlot);
}

//________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::UserExec(Option_t *){
		
	//Main loop
	//Called for each event
	AliVEvent* event = 0x0;
	AliESDEvent* esdEvent = 0x0;
	if(!fIsV0tree){
		event = dynamic_cast<AliVEvent*>(InputEvent());
	}else{
		esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
		//event = dynamic_cast<AliESDEvent*>(InputEvent());
		event = esdEvent;
	}
	if(!event){
		AliError("No event");
		return;
	} 
	if(!fIsV0tree){
		fQAhist->Fill("Events_check",1);
	}

	// check event cuts
	if(IsEventAccepted(event) == 0){ 
		return;
	}
	if(!fIsV0tree){
		fQAhist->Fill("Events_accepted",1);
	}

	//Check if running on MC files
	AliMCEvent* mcEvent = MCEvent();
	if(mcEvent){
		hasMC = kTRUE;
		if(!fIsV0tree){
			fQAhist->Fill("Events_MCcheck",1);
		}
	}
	else{
		hasMC = kFALSE;
	}

	eventNum += 1;

    
	//PID Response task active?
	fPIDResponse = (dynamic_cast<AliInputEventHandler*>((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))->GetPIDResponse();

	if(!fPIDResponse){ AliFatal("This task needs the PID response attached to the inputHandler"); }

	AliVVertex* vertex = const_cast<AliVVertex*>(event->GetPrimaryVertex());

	primaryVertex[0] = vertex->GetX();
	primaryVertex[1] = vertex->GetY();
	primaryVertex[2] = vertex->GetZ(); 


	//Get Multiplicity
	AliMultSelection* multSelection = dynamic_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
	if(!multSelection){
		AliWarning("AliMultSelection object not found");
	}
	else{
		multiplicityV0A = multSelection->GetMultiplicityPercentile("V0A");
		multiplicityV0C = multSelection->GetMultiplicityPercentile("V0C");
		multiplicityCL1 = multSelection->GetMultiplicityPercentile("CL1");
	}

	//Check if ESD or AOD analysis
	//Get ClassName from file and set appropriate ESD or AOD flag
	TString className = static_cast<TString>(event->ClassName());
	if(className.Contains("AOD")){
		fIsAOD = kTRUE;
	}
	else if(className.Contains("ESD")){
		fIsAOD = kFALSE;
	}
	else{
		AliError("-----!! Analysis type unknown !!--------");
		return;
	}
    

	runNumber         = event->GetRunNumber();
	Int_t eventTracks = event->GetNumberOfTracks();
	Int_t numV0s      = event->GetNumberOfV0s();

	AliVParticle* mcTrack = 0x0;
	AliVParticle* motherMCtrack = 0x0;

	//Loop over tracks for event
	if(!fIsV0tree && !fIsEffTree){
		for(Int_t iTrack = 0; iTrack < eventTracks; iTrack++){ 

			AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
			if(!track){
				AliError(Form("Could not receive track %d", iTrack));
				continue;
			}

			fQAhist->Fill("Tracks_all",1);

			//Set MC features to dummy values
			mcEta     = -99;
			mcPhi     = -99;
			mcPt      = -99;
			mcVert[0] = {-99};
			mcVert[1] = {-99};
			mcVert[2] = {-99};
			iPdg         = -9999;
			iPdgMother   = -9999;
			HasMother   = kFALSE; 
			motherLabel  = -9999; //Needed to determine whether tracks have same mother
			//Bool_t IsEnhanced = kFALSE;

			//Get MC information
			if(hasMC){
				if(fIsAOD){
					mcTrack = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));

					//Check valid pointer has been returned. If not, disregard track. 
					if(!mcTrack){
						continue;
					}
				}else{
					mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));

					//Check valid pointer has been returned. If not, disregard track. 
					if(!mcTrack){
						continue;
					}
				}
					
				fQAhist->Fill("Tracks_MCcheck", 1);

				iPdg  = mcTrack->PdgCode();
				mcEta = mcTrack->Eta();
				mcPhi = mcTrack->Phi();
				mcPt  = mcTrack->Pt();
				mcTrack->XvYvZv(mcVert);

				Int_t gMotherIndex = mcTrack->GetMother();
				
				if(!(gMotherIndex < 0)){
					if(fIsAOD){
						motherMCtrack = dynamic_cast<AliAODMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}else{
						motherMCtrack = dynamic_cast<AliMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
                    }
						HasMother = kTRUE;
            iPdgMother = motherMCtrack->PdgCode();
            motherLabel = TMath::Abs(motherMCtrack->GetLabel());
				}
					
				//Currently no injected MC productions for Run 2. Hence commented out
				//Now determine whether track comes from an injected MC sample or not
				/*Int_t mcTrackIndex = -9999;
				while(!(mcTrack->GetMother() < 0)){ 
						mcTrackIndex = mcTrack->GetMother();
						mcTrack = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(mcTrack->GetMother()));
				}
				if(!(mcEvent->IsFromBGEvent(TMath::Abs(mcTrackIndex)))){
						IsEnchanced = kTRUE;
				}else{
						IsEnhanced = kFALSE;
				}*/
      }


			//Apply global track filter
			if(!fIsAOD){
				if(!(fESDtrackCuts->AcceptTrack(dynamic_cast<const AliESDtrack*>(track)))){ 
						continue; 
				}
			}
			else{
				if(!((dynamic_cast<AliAODTrack*>(track)))->TestFilterBit(fFilterBit)){
						continue;
				}
			}

			//Apply some track cuts 
			pt = track->Pt();
			if( pt < fPtMin || pt >= fPtMax ){ continue;}
			eta  = track->Eta();
			if( eta < fEtaMin || eta >= fEtaMax ){ continue;} 
			phi  = track->Phi();

			fQAhist->Fill("Tracks_KineCuts", 1);

			//Get TPC information
			//kNclsTPC
			nTPCclusters = track->GetTPCNcls(); 
			if(nTPCclusters < 70){ continue;}
			
			//kNFclsTPCr
			nTPCcrossed = track->GetTPCClusterInfo(2,1);
			if(nTPCcrossed < 60){ continue;}

			fTPCcrossOverFind = 0;
			nTPCfindable = track->GetTPCNclsF();
			//kNFclsTPCfCross
			if(nTPCfindable > 0){
				fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
			}

			if(fTPCcrossOverFind < 0.3 || fTPCcrossOverFind >= 1.1){ continue;}

			tpcSharedMap = 0;
			if(fIsAOD){
				tpcSharedMap = (dynamic_cast<AliAODTrack*>(track))->GetTPCSharedMap();
			}else{
				tpcSharedMap = (dynamic_cast<AliESDtrack*>(track))->GetTPCSharedMap();
			}

			nTPCshared = -1;
			if(fIsAOD){
				nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
			}else{
				nTPCshared = (dynamic_cast<AliESDtrack*>(track))->GetTPCnclsS();
			}

			chi2TPC = track->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs
			
			//Check for refits 
			if((track->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
			if((track->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

			//DCA values
			Float_t  DCAesd[2] = {0.0, 0.0};
			Double_t DCAaod[2] = {0.0, 0.0};
			Double_t DCAcov[2] = {0.0, 0.0};
			if(!fIsAOD){
				//Arguments: xy, z
				track->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
			}
			else{
				GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(track), DCAaod, DCAcov);
			}
			//Final DCA values stored here 
			if(!fIsAOD){
				DCA[0] = static_cast<Double_t>(DCAesd[0]);
				DCA[1] = static_cast<Double_t>(DCAesd[1]);
			}
			else{
				DCA[0] = static_cast<Double_t>(DCAaod[0]);
				DCA[1] = static_cast<Double_t>(DCAaod[1]);
			}

			//DCAxy cut
			//kImpactParXY
			if(DCA[0] < -3.0 || DCA[0] >= 3.0){ continue;}
			//DCAz cut
			//kImpactParZ
			if(DCA[1] < -4.0 || DCA[1] >= 4.0){ continue;}

			//Get ITS information
			//kNclsITS
			nITS = track->GetNcls(0);

			if(fHasSDD){
				if(nITS < 3){ continue;}
			}else{
				if(nITS < 1){ continue;}
			}
			
			chi2ITS = track->GetITSchi2();
			//kITSchi2Cl
			if((chi2ITS/nITS) >= 36){ continue;} 
			
			fITSshared = 0.;
			for(Int_t d = 0; d < 6; d++){
				fITSshared += static_cast<Double_t>(track->HasSharedPointOnITSLayer(d));
			}
			fITSshared /= nITS;

			if(fIsAOD){
				SPDfirst = (dynamic_cast<AliAODTrack*>(track))->HasPointOnITSLayer(0); 
			}else{
				SPDfirst = (dynamic_cast<AliESDtrack*>(track))->HasPointOnITSLayer(0);
			}

			fQAhist->Fill("Tracks_TrackCuts", 1);


			//Get electron nSigma in TPC for cut (inclusive cut)
			EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
			if( EnSigmaTPC >= fESigTPCMax || EnSigmaTPC < fESigTPCMin) { continue;}
				
			EnSigmaITS = -999;
			if(fHasSDD){
				//Get rest of electron nSigma values and apply cuts if requested (inclusive cuts)
				EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron);
				if(fPIDcutITS){
					if(EnSigmaITS < fESigITSMin || EnSigmaITS > fESigITSMax){ continue;}
				}
			}

			EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron);
			if(fPIDcutTOF){
				if(EnSigmaTOF < fESigTOFMin || EnSigmaTOF > fESigTOFMax){ continue;}
			}

			//Get pion nSigma for TPC and apply cut if requested (exclusive cut)
			PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion);
			if(fPionPIDcutTPC){
				if(PnSigmaTPC > fPSigTPCMin && PnSigmaTPC < fPSigTPCMax){ continue;}
			}
		
			fQAhist->Fill("Tracks_PIDcuts",1); 

			//Get rest of nSigma values for pion and kaon
			PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion);
			PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion);

			KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon);
			KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon);
			KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon);


			//Get ITS and TPC signals
			ITSsignal = track->GetITSsignal();
			TPCsignal = track->GetTPCsignal();
			TOFsignal = track->GetTOFsignal();
           

			Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDtrackCuts->kVertexTracks | fESDtrackCuts->kVertexSPD;

			const AliVVertex* vertex = 0;

			if(fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTracks){
				vertex = track->GetEvent()->GetPrimaryVertexTracks();
			}

			if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexSPD){
				vertex = track->GetEvent()->GetPrimaryVertexSPD();
			}

			if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTPC){
				vertex = track->GetEvent()->GetPrimaryVertexTPC();
			}
			
			//Get golden Chi2
			goldenChi2 = -1;
			if(vertex->GetStatus()){
				if(fIsAOD){
					goldenChi2 = dynamic_cast<AliAODTrack*>(track)->GetChi2TPCConstrainedVsGlobal();
				}
				else{
					goldenChi2 = dynamic_cast<AliESDtrack*>(track)->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex*>(vertex));
				}
			}

			charge = track->Charge();

			fTree->Fill();
    } //End loop over tracks
  }//End of "normal" TTree creation
  else if(fIsV0tree && !fIsEffTree){
		for(Int_t iV0 = 0; iV0 < numV0s; iV0++){

			AliESDv0* V0vertex = esdEvent->GetV0(iV0);
	
			if(!V0vertex){
				AliError(Form("Could not receive V0 track %d", iV0));
				continue;
			}

			//fQAhist->Fill("Tracks_all",1);

			//Get V0 daughter tracks
			AliESDtrack* negTrack = esdEvent->GetTrack(V0vertex->GetIndex(0));
			AliESDtrack* posTrack = esdEvent->GetTrack(V0vertex->GetIndex(1));
			if(!negTrack || !posTrack){
				Printf("Daughter track of v0 not found: %p - %p \n",negTrack, posTrack);
				continue;
			}
			//Check for like-sign V0 candidates
			if(negTrack->Charge() == posTrack->Charge()){ continue; }

			//Apply kinematic and PID cuts to both legs 
			if(isV0daughterAccepted(negTrack) != kTRUE){ continue; }
			if(isV0daughterAccepted(posTrack) != kTRUE){ continue; }

			pointingAngle = V0vertex->GetV0CosineOfPointingAngle();
			daughtersDCA  = V0vertex->GetDcaV0Daughters();
			decayLength   = V0vertex->GetRr();
			v0mass        = V0vertex->M();
			//Super loose cuts on V0 topological qualities(stored in Tree to be cut on later)
			if( pointingAngle < 0.8 || daughtersDCA < 0.05 || decayLength < 0.01 ){ continue; }
			
			//fQAhist->Fill("Arm. cuts",1);

			ptArm = V0vertex->PtArmV0();
			alpha = V0vertex->AlphaV0();
			//Armentors plots is not filled until MC checks completed

		

			//TODO: Improve efficiency of MC section
			//Currently: checks neg particle, then pos particle.
			//If both pass, write pos particle features, then get neg features
			//again!
			Int_t label = -9999;
			if(hasMC){
				//-------- Check negative charge particle --------------
				//Only want to check validity of track
				label = negTrack->GetLabel();

				mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
				if(!mcTrack){
					continue;
				}
				
				Int_t gMotherIndex = mcTrack->GetMother();
		
				if(!(gMotherIndex < 0)){
					if(fIsAOD){
						motherMCtrack = dynamic_cast<AliAODMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}else{
						motherMCtrack = dynamic_cast<AliMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}
				}//End loop over motherIndex
				//--------- End negative particle check	 -------------

				//--------- Get postive particle features
				//Set features for writing into TTree
				 
				label = posTrack->GetLabel();

				mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
				if(!mcTrack){
					continue;
				}

				//Fill MC check
				//fQAhist->Fill("Tracks_MCcheck", 1);

				iPdg = mcTrack->PdgCode();
				
				mcEta = mcTrack->Eta();
				mcPhi = mcTrack->Phi();
				mcPt  = mcTrack->Pt();
				mcTrack->XvYvZv(mcVert);

				gMotherIndex = mcTrack->GetMother();
		
				if(!(gMotherIndex < 0)){
					if(fIsAOD){
						motherMCtrack = dynamic_cast<AliAODMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}else{
						motherMCtrack = dynamic_cast<AliMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}
					HasMother   = kTRUE;
					iPdgMother  = motherMCtrack->PdgCode();
					motherLabel = TMath::Abs(motherMCtrack->GetLabel());
				}//End loop over motherIndex
			}//End ifMC section
	
			fArmPlot->Fill(alpha, ptArm);

			//Get positive particle obsevables
			pt  = posTrack->Pt();
			eta = posTrack->Eta();
			phi = posTrack->Phi();

			EnSigmaITS = -999;
			PnSigmaITS = -999;
			KnSigmaITS = -999;

			if(fHasSDD){
				EnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kElectron);
				PnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kPion);
				KnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kKaon);
			}
			EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kElectron);
			EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kElectron);
			PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kPion);
			PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kPion);
			KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kKaon);
			KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kKaon);
	
			
			//DCA values
			//Initialized here incase goto is invoked
			Float_t DCAesd[2]  = {0.0,0.0};
			Double_t DCAaod[2] = {0.0,0.0};
			Double_t DCAcov[2] = {0.0, 0.0};

			if(TMath::Abs(EnSigmaTPC) > fESigTPCMax){ 
				goto negativeTrack;
			}

			//Get TPC information
			nTPCclusters = posTrack->GetTPCNcls(); 
			if(nTPCclusters < 70){ continue;}
			
			nTPCcrossed = posTrack->GetTPCClusterInfo(2,1);
			if(nTPCcrossed < 60){ continue;}

			fTPCcrossOverFind = 0;

			nTPCfindable = posTrack->GetTPCNclsF();
			if(nTPCfindable > 0){
				fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
			}

			if(fTPCcrossOverFind < 0.3 || fTPCcrossOverFind > 1.1){ continue;}

			tpcSharedMap = 0;
			if(fIsAOD){
				tpcSharedMap = (dynamic_cast<AliAODTrack*>(posTrack))->GetTPCSharedMap();
			}else{
				tpcSharedMap = (dynamic_cast<AliESDtrack*>(posTrack))->GetTPCSharedMap();
			}

			nTPCshared = -1;
			if(fIsAOD){
				nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
			}else{
				nTPCshared = (dynamic_cast<AliESDtrack*>(posTrack))->GetTPCnclsS();
			}

			chi2TPC = posTrack->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs
			
			//Check for refits 
			if((posTrack->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
			if((posTrack->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

			if(!fIsAOD){
				posTrack->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
			}
			else{
				GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(posTrack), DCAaod, DCAcov);
			}
			//Final DCA values stored here 
			if(!fIsAOD){
				DCA[0] = static_cast<Double_t>(DCAesd[0]);
				DCA[1] = static_cast<Double_t>(DCAesd[1]);
			}
			else{
				DCA[0] = static_cast<Double_t>(DCAaod[0]);
				DCA[1] = static_cast<Double_t>(DCAaod[1]);
			}

			//ITS clusters and shared clusters
			nITS = posTrack->GetNumberOfITSClusters();
			fITSshared = 0.;
			for(Int_t d = 0; d < 6; d++){
				fITSshared += static_cast<Double_t>(posTrack->HasSharedPointOnITSLayer(d));
			}
			fITSshared /= nITS;

			SPDfirst = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(0);

			charge = posTrack->Charge();

			//Fill with feature from positive track
			fTree->Fill();

			//Skip to this point if positive track failed PID cut
			negativeTrack:

			//Get negative particle obsevables
			pt  = negTrack->Pt();
			eta = negTrack->Eta();
			phi = negTrack->Phi();

			if(fHasSDD){
				EnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kElectron);
				PnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kPion);
				KnSigmaITS = fPIDResponse->NumberOfSigmasITS(negTrack,(AliPID::EParticleType)AliPID::kKaon);
			}
			EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kElectron);
			EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kElectron);
			PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kPion);
			PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kPion);
			KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(negTrack,(AliPID::EParticleType)AliPID::kKaon);
			KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(negTrack,(AliPID::EParticleType)AliPID::kKaon);

			if(TMath::Abs(EnSigmaTPC) > fESigTPCMax){ 
				continue;
			}

			//Get TPC information
			nTPCclusters = negTrack->GetTPCNcls(); 
			if(nTPCclusters < 70){ continue;}
			
			nTPCcrossed = negTrack->GetTPCClusterInfo(2,1);
			if(nTPCcrossed < 60){ continue;}

			fTPCcrossOverFind = 0;

			nTPCfindable = negTrack->GetTPCNclsF();
			if(nTPCfindable > 0){
				fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
			}

			if(fTPCcrossOverFind < 0.3 || fTPCcrossOverFind > 1.1){ continue;}

			tpcSharedMap = 0;
			if(fIsAOD){
				tpcSharedMap = (dynamic_cast<AliAODTrack*>(negTrack))->GetTPCSharedMap();
			}else{
				tpcSharedMap = (dynamic_cast<AliESDtrack*>(negTrack))->GetTPCSharedMap();
			}

			nTPCshared = -1;
			if(fIsAOD){
				nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
			}else{
				nTPCshared = (dynamic_cast<AliESDtrack*>(negTrack))->GetTPCnclsS();
			}

			chi2TPC = negTrack->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs
			
			//Check for refits 
			if((negTrack->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
			if((negTrack->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

			//DCA values
			if(!fIsAOD){
				negTrack->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
			}
			else{
				GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(negTrack), DCAaod, DCAcov);
			}
			//Final DCA values stored here 
			if(!fIsAOD){
				DCA[0] = static_cast<Double_t>(DCAesd[0]);
				DCA[1] = static_cast<Double_t>(DCAesd[1]);
			}
			else{
				DCA[0] = static_cast<Double_t>(DCAaod[0]);
				DCA[1] = static_cast<Double_t>(DCAaod[1]);
			}

			//ITS clusters and shared clusters
			nITS = negTrack->GetNumberOfITSClusters();
			fITSshared = 0.;
			for(Int_t d = 0; d < 6; d++){
				fITSshared += static_cast<Double_t>(negTrack->HasSharedPointOnITSLayer(d));
			}
			fITSshared /= nITS;

			SPDfirst = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(0);

			charge = negTrack->Charge(); 
			//Write negative observales to tree (v0 information written twice. Filter by looking at only pos or neg charge)
			if(hasMC){
				label = negTrack->GetLabel();

				mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
				//Redundant?
				if(!mcTrack){
					continue;
				}

				//Fill MC check
				//fQAhist->Fill("Tracks_MCcheck", 1);

				iPdg = mcTrack->PdgCode();
				
				mcEta = mcTrack->Eta();
				mcPhi = mcTrack->Phi();
				mcPt  = mcTrack->Pt();
				mcTrack->XvYvZv(mcVert);

				Int_t gMotherIndex = mcTrack->GetMother();
	
				if(!(gMotherIndex < 0)){
					if(fIsAOD){
						motherMCtrack = dynamic_cast<AliAODMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}else{
						motherMCtrack = dynamic_cast<AliMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
						//Check for mother particle. 
						if(!motherMCtrack){
							continue;
						}
					}
					HasMother = kTRUE;
					iPdgMother = motherMCtrack->PdgCode();
					motherLabel = TMath::Abs(motherMCtrack->GetLabel());
				}//End loop over motherIndex
			}

			
			//Fill with features from negative track
			fTree->Fill();
		}//End loop over v0's for this event
  }//End V0 code
	else if(!fIsV0tree && fIsEffTree){

		//Vars used to keep track of particle labels 
		std::vector<std::vector<Int_t>> particleList; //List of all generated particle
		std::vector<Int_t> recoInfo(3,0); //Flags for each track: electron, recod, reco. track number
		std::vector<Int_t> labels; //List of generated electrons (generator track labels))

		
		//First loop over all MC particles and get list of electrons
		for(Int_t iMCtrack = 0; iMCtrack < mcEvent->GetNumberOfTracks(); iMCtrack++){

			if(fIsAOD){
				mcTrack = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(iMCtrack));
			}else{
				mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(iMCtrack));
			}
			if(!mcTrack){
				AliWarning(Form("Could not retreive MC particle %d", iMCtrack));
				continue;
			}
  
			//Check if MC particle is electron. If not, skip to next track
			if(TMath::Abs(mcTrack->PdgCode()) != 11){
				recoInfo[0] = 0;
			}
			else{
				recoInfo[0] = 1;
				labels.push_back(iMCtrack);
			}
			particleList.push_back(recoInfo);
		}//End loop over MC tracks
 
		//Next loop over reconstructed tracks and store which tracks were
		//reconstructed
		for(Int_t iTracks = 0; iTracks < event->GetNumberOfTracks(); iTracks++){
			
			AliVTrack* recoTrack = dynamic_cast<AliVTrack*>(event->GetTrack(iTracks));
      if(!recoTrack){
	      AliError(Form("Could not receive reconstructed track %d", iTracks));
	      continue;
      } 
			//Check if electron
			//GetLabel returns MC track number
			//Take absolute value as some tracks recieve a negative value due to poor
			//quality of track (not important for this step)
      if(particleList[TMath::Abs(recoTrack->GetLabel())][0] == 0){
				continue;     
			}
			//Store flag for reconstruction and reconstruction track number
      particleList[TMath::Abs(recoTrack->GetLabel())][1] = 1;
      particleList[TMath::Abs(recoTrack->GetLabel())][2] = iTracks;
		}

		//Finally, loop over generated electrons
		for(UInt_t iTrack = 0; iTrack < labels.size(); iTrack++){
			Int_t elecLabel = labels[iTrack];

			fQAhist->Fill("Tracks_all", 1);
			if(fIsAOD){
				mcTrack = dynamic_cast<AliAODMCParticle*>(mcEvent->GetTrack(elecLabel));
			}else{
				mcTrack = dynamic_cast<AliMCParticle*>(mcEvent->GetTrack(elecLabel));
			}
			if(!mcTrack){
				Printf("Could not get MC track!");
				AliWarning(Form("Could not retreive MC particle %d", elecLabel));
				continue;
			}
  
			//Check if MC particle is electron. If not, issue error and exit
			if(TMath::Abs(mcTrack->PdgCode()) != 11){
				Printf("Second loop over non-electron. Shouldn't happen. Disregard all results and check code");
				return;
			}
			
			fQAhist->Fill("Tracks_MCcheck", 1);
			//Apply kine cuts (to MC)
			pt   = mcTrack->Pt();
			if( pt < fPtMin || pt > fPtMax ){ continue;}
			eta  = mcTrack->Eta();
			if( eta < fEtaMin || eta > fEtaMax ){ continue;} 
			fQAhist->Fill("Tracks_KineCuts", 1);
			phi  = mcTrack->Phi();
			
			//Get MC information
			//Declare MC variables
			mcEta       = -99;
			mcPhi       = -99;
			mcPt        = -99;
			mcVert[0]   = {-99};
			mcVert[1]   = {-99};
			mcVert[2]   = {-99};
			iPdg        = -9999;
			iPdgMother  = -9999;
			HasMother   = kFALSE;
			motherLabel = -9999;

			iPdg  = mcTrack->PdgCode();
			mcEta = mcTrack->Eta();
			mcPhi = mcTrack->Phi();
			mcPt  = mcTrack->Pt();
			mcTrack->XvYvZv(mcVert);

			//Check for mother particle
			//If not found, dummy MC initialisation values will be written
			Int_t gMotherIndex = mcTrack->GetMother();
			if(!(gMotherIndex < 0)){
				if(fIsAOD){
					motherMCtrack = dynamic_cast<AliAODMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
				}else{
					motherMCtrack = dynamic_cast<AliMCParticle*>((mcEvent->GetTrack(gMotherIndex)));
				}
				HasMother   = kTRUE;
				iPdgMother  = motherMCtrack->PdgCode();
				motherLabel = TMath::Abs(motherMCtrack->GetLabel());
			}
		
			//If reconstructed, get track properties.
			//Otherwise branches will be filled with dummy variables.
			if(particleList[elecLabel][1] == 1){

				//Get reconstructed track
				AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(particleList[elecLabel][2]));
				if(!track){
					Printf("Could not retrieve reconstructed track %d", particleList[elecLabel][2]);
					continue;
				}

				nTPCclusters      = track->GetTPCNcls();
				nTPCcrossed       = track->GetTPCClusterInfo(2,1);
				fTPCcrossOverFind = 0;
				nTPCfindable      = track->GetTPCNclsF();
				if(nTPCfindable > 0){
					fTPCcrossOverFind = nTPCcrossed/nTPCfindable;
				}
				tpcSharedMap = 0;
				if(fIsAOD){
					tpcSharedMap = (dynamic_cast<AliAODTrack*>(track))->GetTPCSharedMap();
				}else{
					tpcSharedMap = (dynamic_cast<AliESDtrack*>(track))->GetTPCSharedMap();
				}

				nTPCshared = -1;
				if(fIsAOD){
					nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
				}else{
					nTPCshared = (dynamic_cast<AliESDtrack*>(track))->GetTPCnclsS();
				}

				//chi2TPC = track->GetTPCchi2(); //Function only implemented in ESDs. Returns dummy value for AODs

				//DCA values
				Float_t DCAesd[2]  = {0.0,0.0};
				Double_t DCAaod[2] = {0.0,0.0};
				Double_t DCAcov[2] = {0.0, 0.0};
				if(!fIsAOD){
					track->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
				}
				else{
					GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(track), DCAaod, DCAcov);
				}
				//Final DCA values stored here 
				if(!fIsAOD){
					DCA[0] = static_cast<Double_t>(DCAesd[0]);
					DCA[1] = static_cast<Double_t>(DCAesd[1]);
				}
				else{
					DCA[0] = static_cast<Double_t>(DCAaod[0]);
					DCA[1] = static_cast<Double_t>(DCAaod[1]);
				}

				//Get ITS values
				nITS = track->GetNcls(0);;
				chi2ITS = track->GetITSchi2();
				fITSshared = 0.;
				for(Int_t d = 0; d < 6; d++){
					fITSshared += static_cast<Double_t>(track->HasSharedPointOnITSLayer(d));
				}
				fITSshared /= nITS;
				SPDfirst = kFALSE;
				if(fIsAOD){
					SPDfirst = (dynamic_cast<AliAODTrack*>(track))->HasPointOnITSLayer(0);
				}else{
					SPDfirst = (dynamic_cast<AliESDtrack*>(track))->HasPointOnITSLayer(0);
				}

				//Get electron nSigma in TPC for cut (inclusive cut)
				EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
					
				EnSigmaITS = -999;
				if(fHasSDD){
					//Get rest of electron nSigma values and apply cuts if requested (inclusive cuts)
					EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron);
				}

				EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron);

				//Get pion nSigma for TPC and apply cut if requested (exclusive cut)
				PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion);
			
				//Get rest of nSigma values for pion and kaon
				/* PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion); */
				/* PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion); */
				/* KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon); */
				/* KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon); */
				/* KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon); */

				//Get raw signals
				ITSsignal = track->GetITSsignal();
				TPCsignal = track->GetTPCsignal();
				TOFsignal = track->GetTOFsignal();
						 
				Int_t fCutMaxChi2TPCConstrainedVsGlobalVertexType = fESDtrackCuts->kVertexTracks | fESDtrackCuts->kVertexSPD;

				const AliVVertex* vertex = 0;

				if(fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTracks){
					vertex = track->GetEvent()->GetPrimaryVertexTracks();
				}

				if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexSPD){
					vertex = track->GetEvent()->GetPrimaryVertexSPD();
				}

				if((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & fESDtrackCuts->kVertexTPC){
					vertex = track->GetEvent()->GetPrimaryVertexTPC();
				}
			
				/* //Get golden Chi2 */
				/* Double_t goldenChi2 = -1; */
				/* if(vertex->GetStatus()){ */
				/* 	if(fIsAOD){ */
				/* 		goldenChi2 = dynamic_cast<AliAODTrack*>(track)->GetChi2TPCConstrainedVsGlobal(); */
				/* 	} */
				/* 	else{ */
				/* 		goldenChi2 = dynamic_cast<AliESDtrack*>(track)->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex*>(vertex)); */
				/* 	} */
				/* } */

				charge = -998;
				charge = track->Charge();

				fTree->Fill();
			}
			else{
				//If not reconstructed
				//Reconstructed branches get dummy values
				//MC branches get correct values (generated values)
				pt                = -9999;
				eta               = -9999;
				phi               = -9999;
				nTPCclusters      = -9999;
				nTPCcrossed       = -9999;
				nTPCfindable      = -9999;
				fTPCcrossOverFind = -9999;
				tpcSharedMap      = -9999;
				nTPCshared        = -9999;
				DCA[0]            = -9999;
				DCA[1]            = -9999;
				nITS              = -9999;
				chi2ITS           = -9999;
				fITSshared        = -9999;
				SPDfirst          = kFALSE;
				EnSigmaTPC        = -9999;
				EnSigmaITS        = -9999;
				EnSigmaTOF        = -9999;
				PnSigmaTPC        = -9999;
				PnSigmaITS        = -9999;
				PnSigmaTOF        = -9999;
				KnSigmaITS        = -9999;
				KnSigmaTPC        = -9999;
				KnSigmaTOF        = -9999;
				ITSsignal         = -9999;
				TPCsignal         = -9999;
				TOFsignal         = -9999;
				goldenChi2        = -9999;
				charge            = -9999;

				fTree->Fill();
			}

		}//End loop over tracks
	}//End efficiency loop

}

//~ //________________________________________________________________________

void  AliAnalysisTaskSimpleTreeMaker::FinishTaskOutput(){
    // Finish task output

    // not implemented ...

}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::Terminate(Option_t *){
    // Draw result to the screen

    // Called once at the end of the query

    // not implemented ...

}
//~ 


//________________________________________________________________________

Int_t AliAnalysisTaskSimpleTreeMaker::IsEventAccepted(AliVEvent *event){
    
    if(!fIsV0tree){
			if(TMath::Abs(event->GetPrimaryVertexSPD()->GetZ()) < 10){
				if(event->GetPrimaryVertexSPD()->GetNContributors() >0){ 
						return 1; 
				}
				else{ 
					return 0;
				}   
			}
    }
    else{
			if(TMath::Abs(event->GetPrimaryVertexSPD()->GetZ()) < 10){
				return 1;
			}
    }
    return 0;
}

Bool_t AliAnalysisTaskSimpleTreeMaker::isV0daughterAccepted(AliVTrack* track){
    
    Bool_t answer = kFALSE;
    //Kinematic cuts
    Double_t pt   = track->Pt();
    if( pt < fPtMin || pt > fPtMax ){ return answer; }
    Double_t eta  = track->Eta();
    if( eta < fEtaMin || eta > fEtaMax ){ return answer; } 

    if(!fIsAOD){
        if(!(fESDtrackCuts->AcceptTrack(dynamic_cast<AliESDtrack*>(track)))){
            return answer;
        }
    }
    else{
        if(!(dynamic_cast<AliAODTrack*>(track))->TestFilterBit(fFilterBit)){
            return answer;
        }
    }
    //fQAhist->Fill("Tracks_KineCuts", 1);
		//Do not apply PID cuts
    //fQAhist->Fill("Tracks_PIDcuts",1); 

    answer = kTRUE;
    return answer;
}

Bool_t AliAnalysisTaskSimpleTreeMaker::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0){
// this is a copy of the AliDielectronVarManager

  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    // the covariance matrix is not stored in case of AliAODTrack::kIsDCA
    return kTRUE;
  }

  Bool_t ok = kFALSE;
  if(event){
    AliExternalTrackParam etp; 
		etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if(xstart>3.){
      d0z0[0] = -999.;
      d0z0[1] = -999.;
      return kFALSE;
    }

    const AliAODVertex *vtx =dynamic_cast<const AliAODVertex*>((event->GetPrimaryVertex()));
    Double_t fBzkG = event->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }

  if(!ok){
    d0z0[0] = -999.;
    d0z0[1] = -999.;
  }
  return ok;
}
