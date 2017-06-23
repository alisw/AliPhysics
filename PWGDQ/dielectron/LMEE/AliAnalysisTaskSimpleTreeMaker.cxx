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

/*************** Tree Maker Class **********************
*                                                      *
* Created: 05.10.2016                                  *
* Authors: Aaron Capon      (aaron.capon@cern.ch)      *
*          Sebastian Lehner (sebastian.lehner@cern.ch) *
*                                                      *
*******************************************************/


ClassImp(AliAnalysisTaskSimpleTreeMaker)

Int_t eventNum = 0;

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker():
    AliAnalysisTaskSE(),
    fTree(0),
    fStream(0),
    fQAhist(0),
    fESDtrackCuts(0),
    fPIDResponse(0),
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
    fPSigTPCMin(-99.),
    fPSigTPCMax(-3.),
    fPIDcutITS(kFALSE),
    fPIDcutTOF(kFALSE),
    fPionPIDcutTPC(kFALSE),
    fIsMC(kTRUE),
    fHasSDD(kTRUE),
    fIsGRIDanalysis(kTRUE),
    fIsV0tree(kFALSE),
    fArmPlot(0),
    fIsAOD(kTRUE),
    fFilterBit(4),
    fGridPID(-1)
{

}

AliAnalysisTaskSimpleTreeMaker::AliAnalysisTaskSimpleTreeMaker(const char *name) :
    AliAnalysisTaskSE(name),
    fTree(0),
    fStream(0),
    fQAhist(0),
    fESDtrackCuts(0),
    fPIDResponse(0),
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
    fPSigTPCMin(-99.),
    fPSigTPCMax(-3.),
    fPIDcutITS(kFALSE),
    fPIDcutTOF(kFALSE),
    fPionPIDcutTPC(kFALSE),
    fIsMC(kTRUE),
    fHasSDD(kTRUE),
    fIsGRIDanalysis(kTRUE),
    fIsV0tree(kFALSE),
    fArmPlot(0),
    fIsAOD(kTRUE),
    fFilterBit(4),
    fGridPID(-1)
{
    fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);

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

void AliAnalysisTaskSimpleTreeMaker::UserCreateOutputObjects() {
  
    AliAnalysisManager* man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
    inputHandler->SetNeedField();
     
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse){
        AliFatal("This task needs the PID response attached to the inputHandler");
      return;
    } 

    fStream = new TTreeStream("tracks");
    fTree   = dynamic_cast<TTree*>(fStream->GetTree());

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
    fArmPlot = new TH2F("ArmPlot", "Armenteros Plot", 100, -1, 1, 100, 0, 0.4);
    if(fIsV0tree){
        fArmPlot->GetXaxis()->SetTitle("#alpha = (p^{+}-p^{-})/(p^{+}+p^{-})");
        fArmPlot->GetYaxis()->SetTitle("p_{T}");
    }  

    fQAhist = new TH1F("h1", "Event and track QA", 9, 0, 9);
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

void AliAnalysisTaskSimpleTreeMaker::UserExec(Option_t *) {
    // Main loop
    // Called for each event
    AliVEvent* event = 0x0;
    AliESDEvent* esdEvent = 0x0;
    if(!fIsV0tree){
        event = dynamic_cast<AliVEvent*>(InputEvent());
    }else{
        esdEvent = dynamic_cast<AliESDEvent*>(InputEvent());
        //event = dynamic_cast<AliESDEvent*>(InputEvent());
        event = esdEvent;
    }
    if(!event) {
        AliError("No event");
        return;
    } 
    fQAhist->Fill("Events_check",1);
    
    // check event cuts
    if(IsEventAccepted(event) == 0){ 
        return;
    }
    fQAhist->Fill("Events_accepted",1);

    //Check if running on MC files
    AliMCEvent* mcEvent = MCEvent();
    if(mcEvent){
        fIsMC = kTRUE;
        fQAhist->Fill("Events_MCcheck",1);
    }
    else{
        fIsMC = kFALSE;
    }
    
    eventNum += 1;

    
    // PID Response task active?
    fPIDResponse = (dynamic_cast<AliInputEventHandler*>((AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler())))->GetPIDResponse();

    if(!fPIDResponse){ AliFatal("This task needs the PID response attached to the inputHandler"); }

    AliVVertex* vertex = const_cast<AliVVertex*>(event->GetPrimaryVertex());

    Double_t primaryVertex[3];
    primaryVertex[0] = vertex->GetX();
    primaryVertex[1] = vertex->GetY();
    primaryVertex[2] = vertex->GetZ(); 


    //Get Multiplicity
    Double_t nMultiplicity = -1;
    AliMultSelection* multSelection = dynamic_cast<AliMultSelection*>(event->FindListObject("MultSelection"));
    if(!multSelection){
        AliWarning("AliMultSelection object not found");
    }
    else{
        nMultiplicity = multSelection->GetMultiplicityPercentile("V0M");
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
    

    Int_t eventTracks = event->GetNumberOfTracks();
    Int_t runNumber = event->GetRunNumber();
    Int_t numV0s = event->GetNumberOfV0s();

    AliVParticle* mcTrack = 0x0;
    AliVParticle* motherMCtrack = 0x0;

    //Loop over tracks for event
    if(!fIsV0tree){
        for(Int_t iTrack = 0; iTrack < eventTracks; iTrack++){ 

            AliVTrack* track = dynamic_cast<AliVTrack*>(event->GetTrack(iTrack));
            if(!track){
                AliError(Form("Could not receive track %d", iTrack));
                continue;
            }

            fQAhist->Fill("Tracks_all",1);

            //Declare MC variables
            Double_t mcEta     = -99;
            Double_t mcPhi     = -99;
            Double_t mcPt      = -99;
			Double_t mcVert[3] = {-99,-99,-99};
            Int_t iPdg         = -9999;
            Int_t iPdgMother   = -9999;
			Bool_t HasMother   = kFALSE; 
            Int_t motherLabel  = -9999; //Needed to determine whether tracks have same mother in evet with many ee pairs
            //Bool_t IsEnhanced = kFALSE;

            //Get MC information
            if(fIsMC){
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

                iPdg = mcTrack->PdgCode();
                
                mcEta = mcTrack->Eta();
                mcPhi = mcTrack->Phi();
                mcPt = mcTrack->Pt();
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
            Double_t pt   = track->Pt();
            if( pt < fPtMin || pt > fPtMax ){ continue;}
            Double_t eta  = track->Eta();
            if( eta < fEtaMin || eta > fEtaMax ){ continue;} 

            fQAhist->Fill("Tracks_KineCuts", 1);

            Double_t nTPCclusters = track->GetTPCNcls(); 

            if( nTPCclusters < 70 ){ continue;}
            
            Double_t nTPCcrossed = track->GetTPCClusterInfo(2,1);
            if(nTPCcrossed < 60){ continue;}

            Double_t TPCcrossOverFind = 0;

            Double_t nTPCfindable = track->GetTPCNclsF();
            if(nTPCfindable > 0){
                TPCcrossOverFind = nTPCcrossed/nTPCfindable;
            }

            if(TPCcrossOverFind < 0.3 || TPCcrossOverFind > 1.1){ continue;}

            TBits tpcSharedMap = 0;
            if(fIsAOD){
                tpcSharedMap = (dynamic_cast<AliAODTrack*>(track))->GetTPCSharedMap();
            }else{
                tpcSharedMap = (dynamic_cast<AliESDtrack*>(track))->GetTPCSharedMap();
            }

            Double_t nTPCshared = -1;
            if(fIsAOD){
                nTPCshared = tpcSharedMap.CountBits(0) - tpcSharedMap.CountBits(159);
            }else{
                nTPCshared = (dynamic_cast<AliESDtrack*>(track))->GetTPCnclsS();
            }

            Double_t chi2TPC = track->GetTPCchi2(); //Function only implemented in ESDs. Returns dumym value for AODs
            
            //Check for refits 
            if((track->GetStatus() & AliVTrack::kITSrefit) <= 0){ continue;}
            if((track->GetStatus() & AliVTrack::kTPCrefit) <= 0){ continue;}

            //DCA values
            Float_t DCAesd[2] = {0.0,0.0};
            Double_t DCAaod[2] = {0.0,0.0};
            Double_t DCAcov[2] = {0.0, 0.0};
            if(!fIsAOD){
                track->GetImpactParameters( &DCAesd[0], &DCAesd[1]);
            }
            else{
                GetDCA(const_cast<const AliVEvent*>(event), dynamic_cast<const AliAODTrack*>(track), DCAaod, DCAcov);
            }
            //Final DCA values stored here 
            Double_t DCA[2];
            if(!fIsAOD){
                DCA[0] = static_cast<Double_t>(DCAesd[0]);
                DCA[1] = static_cast<Double_t>(DCAesd[1]);
            }
            else{
                DCA[0] = static_cast<Double_t>(DCAaod[0]);
                DCA[1] = static_cast<Double_t>(DCAaod[1]);
            }

            if(DCA[0] < -1 || DCA[0] > 1){ continue;}
            if(DCA[1] < -3 || DCA[1] > 3){ continue;}

            Int_t nITS = track->GetNcls(0);;

            if(fHasSDD){
                if(nITS < 4){ continue;}
            }else{
                if(nITS < 2){ continue;}
            }
            
            Double_t chi2ITS = track->GetITSchi2();
            if((chi2ITS/nITS) > 36){ continue;} 
            
            Double_t fITS_shared = 0.;
            for(Int_t d = 0; d < 6; d++){
                fITS_shared += static_cast<Double_t>(track->HasSharedPointOnITSLayer(d));
            }
            fITS_shared /= nITS;

            //Store if SPD has hit in first layer
            Bool_t SPDfirst = kFALSE;
            if(fIsAOD){
                SPDfirst = (dynamic_cast<AliAODTrack*>(track))->HasPointOnITSLayer(0); //Method available for ESDs and AODs
            }else{
                SPDfirst = (dynamic_cast<AliESDtrack*>(track))->HasPointOnITSLayer(0);
            }
     
            fQAhist->Fill("Tracks_TrackCuts", 1);

            Double_t phi  = track->Phi();

            //Get electron nSigma in TPC for cut (inclusive cut)
            Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
            if( EnSigmaTPC > fESigTPCMax || EnSigmaTPC < fESigTPCMin) { continue;}
              
            Double_t EnSigmaITS = -999;
            if(fHasSDD){
                //Get rest of electron nSigma values and apply cuts if requested (inclusive cuts)
                EnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kElectron);
                if(fPIDcutITS){
                    if(EnSigmaITS < fESigITSMin || EnSigmaITS > fESigITSMax){ continue;}
                }
            }

            Double_t EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kElectron);
            if(fPIDcutTOF){
                if(EnSigmaTOF < fESigTOFMin || EnSigmaTOF > fESigTOFMax){ continue;}
            }

            //Get pion nSigma for TPC and apply cut if requested (exclusive cut)
            Double_t PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kPion);
            if(fPionPIDcutTPC){
                if(PnSigmaTPC > fPSigTPCMin && PnSigmaTPC < fPSigTPCMax){ continue;}
            }
          
            fQAhist->Fill("Tracks_PIDcuts",1); 

            //Get rest of nSigma values for pion and kaon
            Double_t PnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kPion);
            Double_t PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kPion);

            Double_t KnSigmaITS = fPIDResponse->NumberOfSigmasITS(track,(AliPID::EParticleType)AliPID::kKaon);
            Double_t KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kKaon);
            Double_t KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track,(AliPID::EParticleType)AliPID::kKaon);


            //Get ITS and TPC signals
            Double_t ITSsignal = track->GetITSsignal();
            Double_t TPCsignal = track->GetTPCsignal();
            Double_t TOFsignal = track->GetTOFsignal();
           

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
            Double_t goldenChi2 = -1;
            if(vertex->GetStatus()){
                if(fIsAOD){
                    goldenChi2 = dynamic_cast<AliAODTrack*>(track)->GetChi2TPCConstrainedVsGlobal();
                }
                else{
                    goldenChi2 = dynamic_cast<AliESDtrack*>(track)->GetChi2TPCConstrainedVsGlobal(dynamic_cast<const AliESDVertex*>(vertex));
                }
            }

            Int_t charge = -10;
            charge = track->Charge();

            //Stream values into tree
            if(fIsMC){
                (*fStream)    << "tracks" <<
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 

                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<

                "charge="     << charge <<
                "ITSsignal="  << ITSsignal <<
                "TPCsignal="  << TPCsignal << 
                "TOFsignal="  << TOFsignal <<
                "vertexX="    << primaryVertex[0] <<
                "vertexY="    << primaryVertex[1] <<
                "vertexZ="    << primaryVertex[2] <<

                "nTPCclusters=" << nTPCclusters <<
                "nTPCcrossed="  << nTPCcrossed <<
                "TPCcrossFind=" << TPCcrossOverFind <<
                "nTPCshared="   << nTPCshared <<
                "chi2TPC="      << chi2TPC <<

                "nITS="        << nITS <<
                "fITSshared="  << fITS_shared << 
                "chi2ITS="     << chi2ITS <<
                "SPDfirst="    << SPDfirst <<

                "DCAxy="        << DCA[0] <<
                "DCAz="         << DCA[1] <<
                "goldenChi2="   << goldenChi2 <<
                "multiplicity=" << nMultiplicity << 

                "mcPt="        << mcPt <<
                "mcEta="       << mcEta <<
                "mcPhi="       << mcPhi <<
				"mcVtX="       << mcVert[0] <<
				"mcVtY="       << mcVert[1] <<
				"mcVtZ="       << mcVert[2] <<
                "pdg="         << iPdg <<
				"hasMother="   << HasMother << 
                "pdgMother="   << iPdgMother <<
                "motherLabel=" << motherLabel << 
                "runNumber="   << runNumber << 
                "eventNum="    << eventNum <<
                "gridPID="     << fGridPID <<
                "\n";
            }
            else{
                (*fStream)    << "tracks" <<
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 

                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<

                "charge="     << charge <<
                "ITSsignal="  << ITSsignal <<
                "TPCsignal="  << TPCsignal << 
                "TOFsignal="  << TOFsignal <<
                "vertexX="    << primaryVertex[0] <<
                "vertexY="    << primaryVertex[1] <<
                "vertexZ="    << primaryVertex[2] <<

                "nTPCclusters=" << nTPCclusters <<
                "nTPCcrossed="  << nTPCcrossed <<
                "TPCcrossFind=" << TPCcrossOverFind <<
                "nTPCshared="   << nTPCshared <<
                "chi2TPC="      << chi2TPC <<

                "nITS="        << nITS <<
                "fITSshared="  << fITS_shared << 
                "chi2ITS="     << chi2ITS <<
                "SPDfirst="    << SPDfirst <<

                "DCAxy="        << DCA[0] <<
                "DCAz="         << DCA[1] <<
                "goldenChi2="   << goldenChi2 <<
                "multiplicity=" << nMultiplicity << 

                "runNumber="  << runNumber << 
                "eventNum="   << eventNum <<
                "gridPID="    << fGridPID <<
                "\n";
            }
        } //End loop over tracks
    }
    else{
        for(Int_t iV0 = 0; iV0 < numV0s; iV0++){
     
            AliESDv0* V0vertex = esdEvent->GetV0(iV0);
        
            if(!V0vertex){
                AliError(Form("Could not receive V0 track %d", iV0));
                continue;
            }

            fQAhist->Fill("Tracks_all",1);

            //Get V0 daughter tracks
            AliESDtrack* negTrack = esdEvent->GetTrack(V0vertex->GetIndex(0));
            AliESDtrack* posTrack = esdEvent->GetTrack(V0vertex->GetIndex(1));
            if(!negTrack | !posTrack){
                Printf("Daughter track of v0 not found: %p - %p \n",negTrack, posTrack);
                continue;
            }
            //Check for like-sign V0 candidates
            if(negTrack->Charge() == posTrack->Charge()){ continue; }

            //Apply kinematic and PID cuts to both legs 
            if(isV0daughterAccepted(negTrack) != kTRUE){ continue; }
            if(isV0daughterAccepted(posTrack) != kTRUE){ continue; }

            Double_t pointingAngle = V0vertex->GetV0CosineOfPointingAngle();
            Double_t daughtersDCA = V0vertex->GetDcaV0Daughters();
            Double_t decayLength = V0vertex->GetRr();
            Double_t v0mass = V0vertex->M();
            //Super loose cuts on V0 topological qualities(stored in Tree to be cut on later)
            if( pointingAngle < 0.8 | daughtersDCA < 0.05 | decayLength < 0.01 ){ continue; }

            Double_t ptArm = V0vertex->PtArmV0();
            Double_t alpha = V0vertex->AlphaV0();
            fArmPlot->Fill(alpha, ptArm);

            fQAhist->Fill("Arm. cuts",1);

            //Get positive particle obsevables
            Double_t pt = posTrack->Pt();
            Double_t eta = posTrack->Eta();
            Double_t phi = posTrack->Phi();

            Double_t EnSigmaITS = -999;
            Double_t PnSigmaITS = -999;
            Double_t KnSigmaITS = -999;
            if(fHasSDD){
                EnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kElectron);
                PnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kPion);
                KnSigmaITS = fPIDResponse->NumberOfSigmasITS(posTrack,(AliPID::EParticleType)AliPID::kKaon);
            }
            Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kElectron);
            Double_t EnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kElectron);
            Double_t PnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kPion);
            Double_t PnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kPion);
            Double_t KnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(posTrack,(AliPID::EParticleType)AliPID::kKaon);
            Double_t KnSigmaTOF = fPIDResponse->NumberOfSigmasTOF(posTrack,(AliPID::EParticleType)AliPID::kKaon);
            //DCA values
            Float_t ImpParamXY = 0.;
            Float_t ImpParamZ = 0.;
            posTrack->GetImpactParameters( &ImpParamXY, &ImpParamZ);

            Int_t nITS = 0;
            Double_t fITS_shared = 0;
            if(fHasSDD){
                //ITS clusters and shared clusters
                nITS = posTrack->GetNumberOfITSClusters();
                for(Int_t d = 0; d < 6; d++){
                      fITS_shared += static_cast<Double_t>(posTrack->HasSharedPointOnITSLayer(d));
                }
                fITS_shared /= nITS;
            }
            Int_t daughtCharge = posTrack->Charge();

            //Declare MC variables
            Double_t mcEta   = -99;
            Double_t mcPhi   = -99;
            Double_t mcPt    = -99;
            Int_t iPdg       = 0;
            Int_t iPdgMother = 0;
            Int_t label = -999;

            if(fIsMC){
                label = posTrack->GetLabel();
				AliVEvent* mcEvent = MCEvent();
                AliVParticle* mcTrack = static_cast<AliVParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
                TParticle* mcParticle = static_cast<TParticle*>(mcTrack->Particle());
                mcEta = mcTrack->Eta();
                mcPhi = mcTrack->Phi();
                mcPt = mcTrack->Pt();

                iPdg = mcTrack->PdgCode();

                Int_t gMotherIndex = mcTrack->GetMother();
                AliVParticle* motherTrack = static_cast<AliVParticle*>(mcEvent->GetTrack(gMotherIndex));
                iPdgMother = motherTrack->PdgCode();
                (*fStream)    << "tracks" <<
                //Positive particle obsevables
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 
                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<
                "nITS="       << nITS <<
                "fITSshared=" << fITS_shared << 
                "impParamXY=" << ImpParamXY <<
                "impParamZ="  << ImpParamZ <<
                "charge="     << daughtCharge <<
                "DCA="        << daughtersDCA <<
                "mcEta="      << mcEta <<
                "mcPhi="      << mcPhi <<
                "mcPt="       << mcPt <<
                "pdg="        << iPdg <<
                "pdgMother="  << iPdgMother <<
                //V0 particle observables 
                "v0effMass="  << v0mass <<
                "pointing="   << pointingAngle << 
                "Rlength="    << decayLength <<
                "ptArm="      << ptArm <<
                "alpha="      << alpha <<
                "\n";
            }
            else{
                (*fStream)    << "tracks" <<
                //Positive particle obsevables
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 
                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<
                "nITS="       << nITS <<
                "nITSshared=" << fITS_shared << 
                "impParamXY=" << ImpParamXY <<
                "impParamZ="  << ImpParamZ <<
                "charge="     << daughtCharge <<
                "DCA="        << daughtersDCA <<
                //V0 particle observables 
                "v0effMass="  << v0mass <<
                "pointing="   << pointingAngle << 
                "Rlength="    << decayLength <<
                "ptArm="      << ptArm <<
                "alpha="      << alpha <<
                "\n";
            }
            //Get negative particle obsevables
            pt = negTrack->Pt();
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
            //DCA values
            ImpParamXY = 0.;
            ImpParamZ = 0.;
            negTrack->GetImpactParameters( &ImpParamXY, &ImpParamZ);
  
            if(fHasSDD){
                //ITS clusters andared clusters
                nITS = negTrack->GetNumberOfITSClusters();
                fITS_shared = 0.;
                for(Int_t d = 0; d < 6; d++){
                    fITS_shared += static_cast<Double_t>(negTrack->HasSharedPointOnITSLayer(d));
                }
                fITS_shared /= nITS;
            }
            daughtCharge = negTrack->Charge(); 
            //Write negative observales to tree (v0 information written twice. Filter by looking at only pos or neg charge)
            if(fIsMC){
                label = negTrack->GetLabel();
				AliVEvent* mcEvent = MCEvent();
                AliVParticle* mcTrack = static_cast<AliVParticle*>(mcEvent->GetTrack(TMath::Abs(label)));
                TParticle* mcParticle = static_cast<TParticle*>(mcTrack->Particle());
                mcEta = mcTrack->Eta();
                mcPhi = mcTrack->Phi();
                mcPt = mcTrack->Pt();

                iPdg = mcTrack->PdgCode();

                Int_t gMotherIndex = mcTrack->GetMother();
                AliVParticle* motherTrack = static_cast<AliVParticle*>((mcEvent->GetTrack(gMotherIndex)));
                iPdgMother = motherTrack->PdgCode();
                (*fStream)    << "tracks" <<
                //Positive particle obsevables
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 
                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<
                "nITS="       << nITS <<
                "nITSshared=" << fITS_shared << 
                "impParamXY=" << ImpParamXY <<
                "impParamZ="  << ImpParamZ <<
                "charge="     << daughtCharge <<
                "DCA="        << daughtersDCA <<
                "mcEta="      << mcEta <<
                "mcPhi="      << mcPhi <<
                "mcPt="       << mcPt <<
                "pdg="        << iPdg <<
                "pdgMother="  << iPdgMother <<
                //V0 particle observables 
                "v0effMass="  << v0mass <<
                "pointing="   << pointingAngle << 
                "Rlength="    << decayLength <<
                "ptArm="      << ptArm <<
                "alpha="      << alpha <<
                "\n";
            }
            else{
                (*fStream)    << "tracks" <<
                //Positive particle obsevables
                "pt="         << pt << 
                "eta="        << eta << 
                "phi="        << phi << 
                "EsigITS="    << EnSigmaITS <<
                "EsigTPC="    << EnSigmaTPC <<
                "EsigTOF="    << EnSigmaTOF <<
                "PsigITS="    << PnSigmaITS <<
                "PsigTPC="    << PnSigmaTPC <<
                "PsigTOF="    << PnSigmaTOF <<
                "KsigITS="    << KnSigmaITS <<
                "KsigTPC="    << KnSigmaTPC <<
                "KsigTOF="    << KnSigmaTOF <<
                "nITS="       << nITS <<
                "nITSshared=" << fITS_shared << 
                "impParamXY=" << ImpParamXY <<
                "impParamZ="  << ImpParamZ <<
                "charge="     << daughtCharge <<
                "DCA="        << daughtersDCA <<
                //V0 particle observables 
                "v0effMass="  << v0mass <<
                "pointing="   << pointingAngle << 
                "Rlength="    << decayLength <<
                "ptArm="      << ptArm <<
                "alpha="      << alpha <<
                "\n";
            }

        }//End loop over v0's for this event
    }

}

//~ //________________________________________________________________________

void  AliAnalysisTaskSimpleTreeMaker::FinishTaskOutput(){
    // Finish task output

    // not implemented ...

}
//~ 

//~ //________________________________________________________________________

void AliAnalysisTaskSimpleTreeMaker::Terminate(Option_t *) {
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
        if(!(fESDtrackCuts->AcceptTrack(static_cast<AliESDtrack*>(track)))){
            return answer;
        }
    }
    else{
        if(!(static_cast<AliAODTrack*>(track))->TestFilterBit(fFilterBit)){
            return answer;
        }
    }
    fQAhist->Fill("Tracks_KineCuts", 1);
    //PID cuts
    Double_t EnSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track,(AliPID::EParticleType)AliPID::kElectron);
    if( EnSigmaTPC > fESigTPCMax || EnSigmaTPC < fESigTPCMin) { return answer; }
    fQAhist->Fill("Tracks_PIDcuts",1); 

    answer = kTRUE;
    return answer;
}

Bool_t AliAnalysisTaskSimpleTreeMaker::GetDCA(const AliVEvent* event, const AliAODTrack* track, Double_t* d0z0, Double_t* covd0z0)
// this is a copy of the AliDielectronVarManager
{
  if(track->TestBit(AliAODTrack::kIsDCA)){
    d0z0[0]=track->DCA();
    d0z0[1]=track->ZAtDCA();
    // the covariance matrix is not stored in case of AliAODTrack::kIsDCA
    return kTRUE;
  }

  Bool_t ok=kFALSE;
  if(event) {
    AliExternalTrackParam etp; etp.CopyFromVTrack(track);

    Float_t xstart = etp.GetX();
    if(xstart>3.) {
      d0z0[0]=-999.;
      d0z0[1]=-999.;
      return kFALSE;
    }

    AliAODVertex *vtx =(AliAODVertex*)(event->GetPrimaryVertex());
    Double_t fBzkG = event->GetMagneticField(); // z componenent of field in kG
    ok = etp.PropagateToDCA(vtx,fBzkG,kVeryBig,d0z0,covd0z0);
  }
  if(!ok){
    d0z0[0]=-999.;
    d0z0[1]=-999.;
  }
  return ok;
}
