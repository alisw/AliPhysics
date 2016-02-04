#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include <AliAnalysisDataContainer.h>
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCGenHandler.h"
#include "AliMCEvent.h"
#include "AliFilteredTreeEventCuts.h"
#include "AliFilteredTreeAcceptanceCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include <TParticle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TH1F.h>
#include "AliESDInputHandlerRP.h"
#include "TRandom.h"
#include "AliAnalysisTaskSEFT2Simulation.h"
#include "FT2.h"
#include "AliSysInfo.h"
#include "AliVertexerTracks.h"
#include "AliTracker.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"

ClassImp(AliAnalysisTaskSEFT2Simulation)

//________________________________________________________________________
AliAnalysisTaskSEFT2Simulation::AliAnalysisTaskSEFT2Simulation():
AliAnalysisTaskSE("FT2Simulation")
, fDebug(0)
, fUseMonitorTree(0)
, fesd(0)
, fESDTree(0)
, fNentries(0)
, fCpuEventWatch(0)
, fRealEventWatch(0)
, fCpuSingleTrackWatch(0)
, fRealSingleTrackWatch(0)
, fFilteredTreeEventCuts(0)
, fFilteredTreeAcceptanceCuts(0)
, fTreeSRedirector(0)
, fMCEffTree(0)
, fDummyTrack(0)
, fDummyFT2Track(0)
, fFET(0)
, fStreamLevel(0)
, fTuneOnDataOrMC(0)
, fSimMat(0)
, fUsePID(0)
, fUseKalman(0)
, fAllowDecay(0)
, fAllowAbsorption(0)
, fUseConverisons(0)
, fdNdY(0)
, fMaxStepTGeo(0)
, fTpcParameterizationFile("")
, fCrossSectionFile("")
, fRunNumber(0)
, fStandaloneTune(0)
{
	AliInfo("Default Constructor");
	
	DefineOutput(1,TH1F::Class());	// Monitor
	DefineOutput(2,TH1F::Class());	// timing plot
	DefineOutput(3,TH1F::Class());	// timing plot
	DefineOutput(4,TH1F::Class());	// timing plot
	DefineOutput(5,TH1F::Class());	// timing plot
	DefineOutput(6,TTree::Class());	// Smeared AliESDEvent
	DefineOutput(7,TTree::Class());	// FT2 performance monitor tree
}

//________________________________________________________________________
AliAnalysisTaskSEFT2Simulation::~AliAnalysisTaskSEFT2Simulation(){
	
	AliInfo("Destructor");
	
	if(fesd != NULL){
		delete fesd;
		fesd = NULL;
	}
	if(fESDTree != NULL){
		delete fESDTree;
		fESDTree = NULL;
	}
	if(fNentries != NULL){
		delete fNentries;
		fNentries = NULL;
	}
	if(fCpuEventWatch != NULL){
		delete fCpuEventWatch;
		fCpuEventWatch = NULL;
	}
	if(fRealEventWatch != NULL){
		delete fRealEventWatch;
		fRealEventWatch = NULL;
	}
	if(fCpuSingleTrackWatch != NULL){
		delete fCpuSingleTrackWatch;
		fCpuSingleTrackWatch = NULL;
	}
	if(fRealSingleTrackWatch != NULL){
		delete fRealSingleTrackWatch;
		fRealSingleTrackWatch = NULL;
	}
	if(fFilteredTreeEventCuts != NULL){
		delete fFilteredTreeEventCuts;
		fFilteredTreeEventCuts = NULL;
	}
	if(fFilteredTreeAcceptanceCuts != NULL){
		delete fFilteredTreeAcceptanceCuts;
		fFilteredTreeAcceptanceCuts = NULL;
	}
	if(fUseMonitorTree){
		if(fTreeSRedirector != NULL){
			delete fTreeSRedirector;
			fTreeSRedirector = NULL;
		}
		if(fMCEffTree != NULL){
			delete fMCEffTree;
			fMCEffTree = NULL;
		}
	}
	if(fDummyTrack != NULL){
		delete fDummyTrack;
		fDummyTrack = NULL;
	}
	if(fDummyFT2Track != NULL){
		delete fDummyFT2Track;
		fDummyFT2Track = NULL;
	}
	if(fFET != NULL){
		delete fFET;
		fFET = NULL;
	}
	
}
//________________________________________________________________________
void AliAnalysisTaskSEFT2Simulation::Init(){
	// Initialization
	
	//	if(fDebug > 1)
	AliInfo("Init() \n");
	
}
//________________________________________________________________________
void AliAnalysisTaskSEFT2Simulation::UserCreateOutputObjects(){
	
	// Create the output container
	//
	AliInfo("Creating the output objects");
	AliSysInfo::AddStamp("Start_Task_UserCreateOutputObjects");
	
	if(fDebug > 1)printf("AliAnalysisTaskSEFT2Simulation::UserCreateOutputObjects() \n");
	//
	if(fUseMonitorTree){
		fTreeSRedirector = new TTreeSRedirector();
		//
		// Create trees
		fMCEffTree = ((*fTreeSRedirector)<<"MCEffTree").GetTree();
		PostData(7,fMCEffTree);
	}
	if (!fDummyTrack)  {
		fDummyTrack=new AliESDtrack();
	}
	if (!fDummyFT2Track)  {
		fDummyFT2Track = new FTProbe();
	}
	
	if(fesd == NULL){
		fESDTree = new TTree("esdTree", "Tree with ESD objects");
		fesd =  new AliESDEvent();
		fesd->CreateStdContent();
		fesd->WriteToTree(fESDTree);
	}
	
	fNentries = new TH1F("Monitor", "Monitor",10,0,10);
	fNentries->GetXaxis()->SetBinLabel(1,"# of Events");
	fNentries->GetXaxis()->SetBinLabel(2,"ProcessTrack success");
	fNentries->GetXaxis()->SetBinLabel(3,"ProcessTrack fail");
	fNentries->GetXaxis()->SetBinLabel(4,"HIJING Tracks");
	fNentries->GetXaxis()->SetBinLabel(5,"non-HIJING Tracks");
	fNentries->GetXaxis()->SetBinLabel(6,"Related to prim. Vtx.");
	fNentries->GetXaxis()->SetBinLabel(7,"!Related to prim. Vtx.");
	fNentries->GetXaxis()->SetBinLabel(8,"Track in ESD event");
	
	fCpuEventWatch	= new TH1F("fCpuEventWatch","CPU time of single event processing;t_{cpu} (s);counts",600,0.,60.);
	fRealEventWatch	= new TH1F("fRealEventWatch","Real time of single event processing;t_{real} (s);counts",600,0.,60.);
	fCpuSingleTrackWatch	= new TH1F("fCpuSingleTrackWatch","CPU time of 1000 tracks;t_{cpu} (s);counts",1000,0.,0.1);
	fRealSingleTrackWatch	= new TH1F("fRealSingleTrackWatch","Real time of 1000 tracks;t_{real} (s);counts",1000,0.,0.1);
	
	// initialize fast estimation tool
	fFET = new FT2();
	fFET->SetStreamLevel(fStreamLevel);
	fFET->SetRunNumber(fRunNumber);
	fFET->SetStandaloneTune(fStandaloneTune);
	fFET->InitEnvLocal();
	fFET->InitDetector();
	fFET->SetTuneOnDataOrMC(fTuneOnDataOrMC);
	fFET->SetSimMat(fSimMat);
	fFET->InitTPCParaFile(fTpcParameterizationFile.Data());
	fFET->InitXSectionFile(fCrossSectionFile.Data());
	fFET->InitTPCPIDResponse();
	fFET->SetUsePIDForTracking(fUsePID);
	fFET->SetMaxStepTGeo(fMaxStepTGeo);
	fFET->SetdNdY(fdNdY);
	fFET->SetUseKalmanOut(fUseKalman);
	fFET->SetAllowDecay(fAllowDecay);
	fFET->SetAllowAbsorbtion(fAllowAbsorption);
	fFET->SetUseConversionExtension(fUseConverisons);
	fFET->PrintLayout();
	printf("\n\n\n###################################################\n### FT2 is now setup with:\n### TuneOnDataOrMC?        %i\n### Simulate Material?     %i\n### Use PID for tracking?  %i\n### MaxStepTGeo?           %0.2f\n### dNdY?                  %4.f\n### Use Kalman?            %i\n### Allow Decay?           %i\n### Allow Absorption       %i\n### TPC Param. from?       %s\n### X-Sections from?       %s\n### Runnumber?             %i\n### Magentic Field?        %0.5f\n### Streamer Level?        %i\n### Use conversion algo.?  %i\n### Use standlone tune?    %i\n###################################################\n\n\n",fTuneOnDataOrMC,fSimMat,fUsePID,fMaxStepTGeo,fdNdY,fUseKalman,fAllowDecay,fAllowAbsorption,fTpcParameterizationFile.Data(),fCrossSectionFile.Data(),fRunNumber,fFET->GetMagneticField(),fStreamLevel,fUseConverisons,fStandaloneTune);
	
	// Post the data
	
	PostData(1,fNentries);
	PostData(2,fCpuSingleTrackWatch);
	PostData(3,fRealSingleTrackWatch);
	PostData(4,fCpuEventWatch);
	PostData(5,fRealEventWatch);
	PostData(6,fESDTree);
	
	AliSysInfo::AddStamp("Stop_Task_UserCreateOutputObjects");
	
	return;
}
//________________________________________________________________________
void AliAnalysisTaskSEFT2Simulation::UserExec(Option_t */*option*/){
	//
	// Execute analysis for current event:
	//
	AliInfo("Analyzing Event");
	
	AliSysInfo::AddStamp("Start_Task_UserExec");
	TStopwatch tWatchEvent;
	tWatchEvent.Start();
	
	AliStack* stack=0;
	
	AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	if (!eventHandler) {
		Printf("ERROR: Could not retrieve MC event handler\n");
		tWatchEvent.Stop();
		return;
	}
	AliMCEvent* mcEvent = eventHandler->MCEvent();
	if (!mcEvent) {
		Printf("ERROR: Could not retrieve MC event\n");
		tWatchEvent.Stop();
		return;
	}
	stack = mcEvent->Stack();
	if (!stack) {
		Printf("ERROR: stack not available\n");
		tWatchEvent.Stop();
		return;
	}
	//
	// get selection cuts
	AliFilteredTreeEventCuts *evtCuts = GetEventCuts();
	AliFilteredTreeAcceptanceCuts *accCuts = GetAcceptanceCuts();
	
	if(!evtCuts || !accCuts) {
		AliDebug(AliLog::kError, "Cuts not available");
		return;
	}
	//
	fFET->SetMCTrueTrackMultiplicity(GetMCTrueTrackMult(mcEvent,evtCuts,accCuts));
	Double_t bz = fFET->GetMagneticField();
	fesd->SetMagneticField(bz);
	//
	//
	//
	// MC truth primary vertex
	TArrayF vertex(3);
	AliVertex *primaryVertex = (AliVertex*)mcEvent->GetPrimaryVertex();
	if(fDebug>2)primaryVertex->Print("");
	vertex[0] = primaryVertex->GetX();
	vertex[1] = primaryVertex->GetY();
	vertex[2] = primaryVertex->GetZ();
	//
	AliESDVertex *ESDprimVtx = SmearPrimaryVertex(vertex[0],vertex[1],vertex[2]);
	ESDprimVtx->SetTitle("smearMC");
	//
	fNentries->Fill(0);
	//
	if(fDebug>1) {
		printf("Number of stack entires: %d \n",stack->GetNtrack());
	}
	//
	Int_t stopwatchCounter = 0;
	TStopwatch tWatchSingleTrack;
	
	TObjArray nonHijingESDTracks;
	
	//Loop over particle tracks
	TParticle *part = NULL;
	Bool_t isPrim=kFALSE;
	Bool_t isFromStrangess=kFALSE;
	Bool_t isFromConversion=kFALSE;
	Bool_t isFromMaterial=kFALSE;
	Bool_t isHIJINGparticle = kTRUE;
	Bool_t isFastWeak=kFALSE;
	Double_t deltaOrigin=0.;

	for(Int_t k = 0 ; k < stack->GetNtrack(); k++) {
		if(fDebug>2)printf("Beginning Particle Loop : %i\n",k);
		if(stopwatchCounter==0 || stopwatchCounter==1000){
			tWatchSingleTrack.Start();
		}
		// Create MC truth track
		part = (TParticle*)stack->Particle(k);
		if(!part){
			continue;
		}
		//
		// checking pdg code
		Int_t partAbsPdg = TMath::Abs(part->GetPdgCode());
		if(partAbsPdg!=211 && partAbsPdg!=321 && partAbsPdg!=2212 && partAbsPdg!=13 && partAbsPdg!=11){
			continue;
		}
		//
		// is particle in acceptance?
		if(!accCuts->AcceptTrack(part)) continue;
		// ---------
		// HIJIJNG particle selection
		TString nameGen="";
		GetTrackPrimaryGenerator(k,mcEvent,stack,nameGen);
		if(!nameGen.Contains("ijing")){
			isHIJINGparticle=kFALSE;
		}
		// ---------
		//TVector3 dVertex(part->Vx() - vertex[0],
		//								 part->Vy() - vertex[1],
		//								 part->Vz() - vertex[2]);
		//if (dVertex.Mag() > 0.0001) {continue;}
		//
		//
		//
		isFromStrangess  = IsFromStrangeness(k,stack);
		isFromConversion = IsFromConversion(k,stack);
		isFromMaterial   = IsFromMaterial(k,stack);
		isFastWeak			 = IsFastMcFromWeakDecay(k,stack);
		Bool_t prim = stack->IsPhysicalPrimary(k);
		
		if(isFromMaterial) isFastWeak=kFALSE;

		Double_t deltaX = part->Vx()-primaryVertex->GetX();
		Double_t deltaY = part->Vy()-primaryVertex->GetY();
		
		deltaOrigin = TMath::Sqrt(deltaX*deltaX+deltaY*deltaY);

		//if(!isFromConversion) continue;
		//if(!prim) continue;
		isPrim = prim;
		if(fDebug>2){
			Int_t mfl = (Int_t)(TMath::Abs(part->GetPdgCode())/TMath::Power(10, Int_t(TMath::Log10(TMath::Abs(part->GetPdgCode())))));
			printf("Particle Info\n");
			part->Print();
			printf("Number of daughters: %i\n",part->GetNDaughters());
			printf("Status code: %i\n",part->GetStatusCode());
			printf("Pdg code: %i\n",part->GetPdgCode());
			printf("Physical primary? %i\n",prim);
			printf("This is the generator name: %s\n",nameGen.Data());
			printf("Mother flavour? %i\n",mfl);
			printf("Mother Info\n");
			if(part->GetFirstMother()>0){
				TParticle* partMum = (TParticle*)stack->Particle(part->GetFirstMother());
				partMum->Print();
			}
		}
		//______________________________________________________________________________________
		FTProbe* ft2Track			= NULL;
		AliESDtrack *recTrack = NULL;
		Bool_t isFT2accept		= kFALSE;
		Bool_t isFt2Fake			= kFALSE;
		Bool_t isFT2TCaccept	= kFALSE;
		Double_t ft2DcaXY, ft2DcaZ, ft2NormDcaXY, ft2NormDcaZ = 0;
		
		if(fDebug>2)printf("*** FT2:ProcessTrack ***\n");
		if(fFET->ProcessTrack(part,primaryVertex)){
			if(stopwatchCounter==999){
				tWatchSingleTrack.Stop();
				fCpuSingleTrackWatch->Fill(tWatchSingleTrack.CpuTime());
				fRealSingleTrackWatch->Fill(tWatchSingleTrack.RealTime());
			}
			isFT2accept = kTRUE;
			ft2Track = &fFET->GetProbe();
			//
			const double* dca = fFET->GetDCA();
			const double* cov = fFET->GetDCACov();
			ft2DcaXY			= dca[0];
			ft2NormDcaXY	= dca[0]/TMath::Sqrt(cov[0]);
			ft2DcaZ				= dca[1];
			ft2NormDcaZ		= dca[1]/TMath::Sqrt(cov[2]);
			//
			isFT2TCaccept = ApplyTrackCuts(ft2Track);
			//
			if (fFET->GetNClITSFakes()) { isFt2Fake = kTRUE; }
			//
			fNentries->Fill(1);
			//
			AliESDtrack *trackESD = PrepareESDtrack(ft2Track,part,k);
			if(!trackESD) {
				delete trackESD;
				trackESD = fDummyTrack;
			}
			recTrack = trackESD;
			if(isHIJINGparticle){
				AddTrack(trackESD);
				fNentries->Fill(3);
			}
			else {nonHijingESDTracks.Add(trackESD);}
		}
		else {
			tWatchSingleTrack.Stop();
			fNentries->Fill(2);
			ft2Track = fDummyFT2Track;// &fFET->GetProbe();
			recTrack = fDummyTrack;
		}
		
		// downscale low-pT particles
		if(fUseMonitorTree){
			Double_t scalempt= TMath::Min(part->Pt(),10.);
			Double_t downscaleF = gRandom->Rndm();
			downscaleF *= 200.;
			if(fTreeSRedirector && (TMath::Exp(2*scalempt)>=downscaleF)) {
				(*fTreeSRedirector)<<"MCEffTree"<<									//
				"Bz="<<bz<<																					// magnetic field
				"ft2Track.="<<ft2Track<<														// reconstructed FT2 track
				"esdTrack.="<<recTrack<<														// reconstructed track (only the longest from the loopers)
				"particle.="<<part<<																// MC generated particle
				"isHijing="<<isHIJINGparticle<<											// is particle from HIJING?
				"isPrim="<<isPrim<<																	// is primary particle?
				"isFromStrangess="<<isFromStrangess<<								// is particle from strange decay?
				"isFromConversion="<<isFromConversion<<							// is particle from conversion?
				"isFromMaterial="<<isFromMaterial<<									// is particle from material
				"deltaOrigin="<<deltaOrigin<<												// calculates the distance of the particle production MC vertex to the MC primary vertex in xy
				"isFastWeak="<<isFastWeak<<
				"isFt2Acc="<<isFT2accept<<													// track accepted by fast simulation tool
				"isFt2TCAcc="<<isFT2TCaccept<<											// track accepted by fast simulation tool with additional track cuts
				"isFt2Fake="<<isFt2Fake<<														// reconstructed FT2 fake track
				"ft2DcaXY="<<ft2DcaXY<<														  // DCA XY to primary vertex in FT2
				"ft2NormDcaXY="<<ft2NormDcaXY<<											// normalized DCA XY to primary vertex in FT2
				"ft2DcaZ="<<ft2DcaZ<<																// DCA Z to primary vertex in FT2
				"ft2NormDcaZ="<<ft2NormDcaZ<<												// normalized DCA Z to primary vertex in
				"\n";
			}
		}
		if(stopwatchCounter==999)stopwatchCounter=0;
		else stopwatchCounter++;
	}
	
	fesd->SetDiamond(ESDprimVtx);
	AliVertexerTracks vertexer(bz);
	vertexer.SetITSMode();
	vertexer.SetConstraintOff();
	vertexer.SetVtxStart(ESDprimVtx);
	AliESDVertex *pvertex = vertexer.FindPrimaryVertex(fesd);//,ESDprimVtx);
	pvertex->SetTitle("smearMC");
	//
	fesd->SetPrimaryVertexTracks(pvertex);
	//
	
	for (Int_t i=0; i<fesd->GetNumberOfTracks(); i++) {
		AliESDtrack *trackESD = fesd->GetTrack(i);
		Double_t x[3]; trackESD->GetXYZ(x);
		Double_t b[3]; AliTracker::GetBxByBz(x,b);
		
		if(!trackESD->RelateToVertexBxByBz(pvertex,b,kVeryBig)){
			//	fesd->RemoveTrack(i);
			//	trackESD->Delete();
			fNentries->Fill(6);
		}
		else{fNentries->Fill(5);}
	}
	
	// adding non-HIJING ESD tracks to fesd event without affecting the primary vertex
	Double_t numberNonHIJINGTracks = nonHijingESDTracks.GetEntriesFast();
	if(fDebug>2)AliInfo(Form("There are %.f non-HIJING tracks.",numberNonHIJINGTracks));
	for(Int_t j=0;j<numberNonHIJINGTracks;j++){
		AliESDtrack *trackESDnonHIJING = (AliESDtrack*)nonHijingESDTracks.At(j);
		AddTrack(trackESDnonHIJING);
		nonHijingESDTracks.Remove(trackESDnonHIJING);
		fNentries->Fill(4);
	}
	// check if all tracks are used
	if(nonHijingESDTracks.GetEntriesFast()!=0){
		AliFatal("==> Not all injected tracks were added to the ESD event.");
		return;
	}
	else{
		if(fDebug>2)AliInfo(Form("==> All %.f tracks were added to the ESD event",numberNonHIJINGTracks));
	}
	
	fNentries->Fill(7,fesd->GetNumberOfTracks());
	if(fDebug>2)pvertex->Print();
	fesd->SetRunNumber(fRunNumber);
	
	fESDTree->Fill();
	Reset();
	
	tWatchEvent.Stop();
	fCpuEventWatch->Fill(tWatchEvent.CpuTime());
	fRealEventWatch->Fill(tWatchEvent.RealTime());
	
	// Post the data
	PostData(1,fNentries);
	PostData(2,fCpuSingleTrackWatch);
	PostData(3,fRealSingleTrackWatch);
	PostData(4,fCpuEventWatch);
	PostData(5,fRealEventWatch);
	PostData(6,fESDTree);
	if(fUseMonitorTree){PostData(7,fMCEffTree);}
	
	AliSysInfo::AddStamp("Stop_Task_UserExec");
	
	return;
}
//____________________________________________________________________________
void AliAnalysisTaskSEFT2Simulation::GetTrackPrimaryGenerator(Int_t label,AliMCEvent *mcEvent,AliStack *stack,TString &nameGen){
	// method to check if a track comes from a given generator
	Int_t lab=TMath::Abs(label);
	nameGen=GetGenerator(lab,mcEvent);
	
	while(nameGen.IsWhitespace()){
		TParticle* part = (TParticle*)stack->Particle(lab);
		if(!part){
			printf("AliAnalysisTaskSEFT2Simulation::GetTrackPrimaryGenerator - BREAK: No valid AliESDMCParticle at label %i\n",lab);
			break;
		}
		Int_t mother = part->GetFirstMother();
		if(mother<0){
			printf("AliAnalysisTaskSEFT2Simulation::GetTrackPrimaryGenerator - BREAK: Reached primary particle without valid mother\n");
			break;
		}
		lab=mother;
		nameGen=GetGenerator(mother,mcEvent);
	}
	
	return;
}
//____________________________________________________________________________
TString AliAnalysisTaskSEFT2Simulation::GetGenerator(Int_t label, AliMCEvent* mcEvent){
	//----
	TString genname=mcEvent->GenEventHeader()->ClassName();
	TList* lgen = 0x0;
	Int_t  nsumpart = 0;
	
	if(genname.Contains("CocktailEventHeader")){
		AliGenCocktailEventHeader *cockhead=(AliGenCocktailEventHeader*)mcEvent->GenEventHeader();
		lgen=cockhead->GetHeaders();
		Int_t nh=lgen->GetEntries();
		for(Int_t i=0;i<nh;i++){
			AliGenEventHeader* gh=(AliGenEventHeader*)lgen->At(i);
			genname=gh->GetName();
			Int_t npart=gh->NProduced();
			if(label>=nsumpart && label<(nsumpart+npart)) return genname;
			nsumpart+=npart;
		}
	}
	TString empty="";
	return empty;
}
//____________________________________________________________________________
Bool_t AliAnalysisTaskSEFT2Simulation::ApplyTrackCuts(FTProbe *FETtrack){
	//
	// apply some track quality cuts here, which save pure computing time! can be removed
	//
	if(FETtrack->fProbeNClITS<4)														return kFALSE;
	if(FETtrack->fProbeNClTPC<=0)														return kFALSE;
	if(FETtrack->fProbeChi2TPC/FETtrack->fProbeNClTPC>4.0)	return kFALSE;
	if(FETtrack->fProbeNClTPC<70)														return kFALSE;
	if(FETtrack->Pt()<0.5)																	return kFALSE;
	if(TMath::Abs(FETtrack->Eta())>0.90)										return kFALSE;
	if(!FETtrack->fProbeITSPattern&(0x1<<0) && (!FETtrack->fProbeITSPattern&(0x1<<1) || !FETtrack->fProbeITSPattern&(0x1<<2)))return kFALSE;
	if(FETtrack->fLostInItsTpcMatching!=0) return kFALSE;
	//
	return kTRUE;
}
//____________________________________________________________________________
AliESDtrack * AliAnalysisTaskSEFT2Simulation::PrepareESDtrack(FTProbe *FETtrack,TParticle *part, Int_t iParticle){
	
	AliESDtrack *trackESD = new AliESDtrack(part);
	if(!trackESD)return 0x0;
	
	//
	// this is a flag in the FT2, but not a cut; need to cut here; must not be removed
	if(FETtrack->fLostInItsTpcMatching!=0) return 0x0;
	//
	//
	//
	// apply some track quality cuts here, which save pure computing time! can be removed
	/*
	 if(FETtrack->fProbeNClITS<4)														return 0x0;
	 if(FETtrack->fProbeNClTPC<=0)														return 0x0;
	 if(FETtrack->fProbeChi2TPC/FETtrack->fProbeNClTPC>4.0)	return 0x0;
	 if(FETtrack->fProbeNClTPC<70)														return 0x0;
	 if(FETtrack->Pt()<0.5)																	return 0x0;
	 if(TMath::Abs(FETtrack->Eta())>0.90)										return 0x0;
	 if(!FETtrack->fProbeITSPattern&(0x1<<0) && (!FETtrack->fProbeITSPattern&(0x1<<1) || !FETtrack->fProbeITSPattern&(0x1<<2)))return 0x0;
	 */
	
	const Double_t *covar = FETtrack->GetCovariance();
	const Double_t *param = FETtrack->GetParameter();
	Double_t alpha				= FETtrack->GetAlpha();
	Double_t xref					= FETtrack->GetX();
	
	trackESD->Set(xref,alpha,param,covar);
	
	trackESD->SetLabel(iParticle);
	trackESD->SetStatus(AliESDtrack::kITSrefit);	// -->  4
	trackESD->SetStatus(AliESDtrack::kITSin);			// --> 68
	trackESD->SetStatus(AliESDtrack::kITSout);		// --> 69
	trackESD->SetStatus(AliESDtrack::kITSupg);
	
	trackESD->SetStatus(AliESDtrack::kTPCrefit);	// --> 64
	trackESD->SetStatus(AliESDtrack::kTPCin);
	trackESD->SetStatus(AliESDtrack::kTPCout);
	trackESD->SetStatus(AliESDtrack::kTPCpid);
	
	trackESD->SetStatus(AliESDtrack::kTOFout);
	trackESD->SetStatus(AliESDtrack::kTIME);
	trackESD->SetStatus(AliESDtrack::kTOFpid);
	//	trackESD->SetStatus(AliESDtrack::kTOFmismatch);
	
	trackESD->SetTPCPointsF(fFET->GetNClTPC());
	trackESD->SetTPCsignal(FETtrack->GetTPCsignal(),0.,FETtrack->GetTPCsignalN());
	trackESD->SetTPCsignalTunedOnData(FETtrack->GetTPCsignal());
	//	trackESD->fTPCsignalN = FETtrack->GetTPCsignal();
	//	trackESD->fTPCnclsIter1 = FETtrack->GetTPCsignal();
	
	
	UChar_t hits = fFET->GetITSPattern();
	UChar_t hitsF = fFET->GetITSPatternFakes();
	trackESD->SetITSClusterMap(hits);
	trackESD->SetITSSharedMap(hitsF);
	
	Char_t itsCluster = 0;
	for (int j=7;j--;){
		if ((hits&(0x1<<j)) || (hitsF&(0x1<<j))){
			itsCluster++;
		}
	}
	// new setters
	trackESD->SetTPCNcls(fFET->GetNClTPC());
	trackESD->SetTPCchi2(fFET->GetChi2TPC());
	trackESD->SetITSchi2(fFET->GetChi2ITS());
	trackESD->SetITSNcls(itsCluster);
	
	trackESD->SetTPCClusterMap(fFET->GetTPCHitMap());
	
	Double_t times[5];
	trackESD->GetIntegratedTimes(times);
	
	return trackESD;
}
//____________________________________________________________________________
AliESDVertex *AliAnalysisTaskSEFT2Simulation::SmearPrimaryVertex(Float_t x, Float_t y, Float_t z)
{
	//	ESDprimVtx->SetNContributors(contributors);
	Double_t position[3];
	position[0]=x;
	position[1]=y;
	position[2]=z;
	Double_t covmatrix[6] = {10.,0.,10.,0.,0.,10.};
	Double_t chi2 = 4.;
	Int_t nContributors = 100;
	AliESDVertex *ESDprimVtx = new AliESDVertex(position,covmatrix,chi2,nContributors,"smearMC");
	if(fDebug>3){
		printf("Smeared Vertex:\n");
		ESDprimVtx->Print();
	}
	
	return ESDprimVtx;
}

//____________________________________________________________________________
void AliAnalysisTaskSEFT2Simulation::Terminate(Option_t */*option*/){
	// Terminate analysis
	//
	if(fDebug > 1) printf("AliAnalysisTaskSEFT2Simulation::Terminate() \n");
	return;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskSEFT2Simulation::GetMCTrueTrackMult(AliMCEvent *const mcEvent, AliFilteredTreeEventCuts *const evtCuts, AliFilteredTreeAcceptanceCuts *const accCuts)
{
	//
	// calculate mc event true track multiplicity
	//
	if(!mcEvent) return 0;
	
	AliStack* stack = 0;
	Int_t mult = 0;
	
	// MC particle stack
	stack = mcEvent->Stack();
	if (!stack) return 0;
	//
	Bool_t isEventOK = evtCuts->AcceptMCEvent(mcEvent);
	if(!isEventOK) return 0;
	
	Int_t nPart  = stack->GetNtrack();
	for (Int_t iMc = 0; iMc < nPart; ++iMc)
	{
		TParticle* particle = stack->Particle(iMc);
		if (!particle)
			continue;
		
		// only charged particles
		if(!particle->GetPDG()) continue;
		Double_t charge = particle->GetPDG()->Charge()/3.;
		if (TMath::Abs(charge) < 0.001)
			continue;
		
		// physical primary
		Bool_t prim = stack->IsPhysicalPrimary(iMc);
		if(!prim) continue;
		
		// checked accepted without pt cut
		//if(accCuts->AcceptTrack(particle))
		if( particle->Eta() > accCuts->GetMinEta() && particle->Eta() < accCuts->GetMaxEta() )
		{
			mult++;
		}
	}
	
	return mult;
}
//_____________________________________________________________________________
TParticle *AliAnalysisTaskSEFT2Simulation::GetMother(TParticle *const particle, AliStack *const stack)
{
	if(!particle) return NULL;
	if(!stack) return NULL;
	
	Int_t motherLabel = TMath::Abs(particle->GetMother(0));
	TParticle* mother = NULL;
	Int_t mcStackSize=stack->GetNtrack();
	if (motherLabel>=mcStackSize) return NULL;
	mother = stack->Particle(motherLabel);
	
	return mother;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEFT2Simulation::IsFromConversion(Int_t label, AliStack *const stack)
{
	Bool_t isFromConversion = kFALSE;
	
	if(stack) {
		Int_t mcStackSize=stack->GetNtrack();
		if (label>=mcStackSize) return kFALSE;
		TParticle* particle = stack->Particle(label);
		if (!particle) return kFALSE;
		
		if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0)
		{
			Int_t mech = particle->GetUniqueID(); // production mechanism
			Bool_t isPrim = stack->IsPhysicalPrimary(label);
			
			Int_t motherLabel = TMath::Abs(particle->GetMother(0));
			if (motherLabel>=mcStackSize) return kFALSE;
			TParticle* mother = stack->Particle(motherLabel);
			if(mother) {
				Int_t motherPdg = mother->GetPdgCode();
				
				if(!isPrim && mech==5 && motherPdg==kGamma) {
					isFromConversion=kTRUE;
				}
			}
		}
	}
	
	return isFromConversion;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEFT2Simulation::IsFromMaterial(Int_t label, AliStack *const stack)
{
	Bool_t isFromMaterial = kFALSE;
	
	if(stack) {
		Int_t mcStackSize=stack->GetNtrack();
		if (label>=mcStackSize) return kFALSE;
		TParticle* particle = stack->Particle(label);
		if (!particle) return kFALSE;
		
		if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0)
		{
			Int_t mech = particle->GetUniqueID(); // production mechanism
			Bool_t isPrim = stack->IsPhysicalPrimary(label);
			
			Int_t motherLabel = TMath::Abs(particle->GetMother(0));
			if (motherLabel>=mcStackSize) return kFALSE;
			TParticle* mother = stack->Particle(motherLabel);
			if(mother) {
				if(!isPrim && mech==13) {
					isFromMaterial=kTRUE;
				}
			}
		}
	}
	
	return isFromMaterial;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskSEFT2Simulation::IsFromStrangeness(Int_t label, AliStack *const stack)
{
	Bool_t isFromStrangeness = kFALSE;
	
	if(stack) {
		Int_t mcStackSize=stack->GetNtrack();
		if (label>=mcStackSize) return kFALSE;
		TParticle* particle = stack->Particle(label);
		if (!particle) return kFALSE;
		
		if(particle && particle->GetPDG() && particle->GetPDG()->Charge()!=0)
		{
			Int_t mech = particle->GetUniqueID(); // production mechanism
			Bool_t isPrim = stack->IsPhysicalPrimary(label);
			
			Int_t motherLabel = TMath::Abs(particle->GetMother(0));
			if (motherLabel>=mcStackSize) return kFALSE;
			TParticle* mother = stack->Particle(motherLabel);
			if(mother) {
				Int_t motherPdg = mother->GetPdgCode();
				
				// K+-, lambda, antilambda, K0s decays
				if(!isPrim && mech==4 &&
					 (TMath::Abs(motherPdg)==kKPlus || TMath::Abs(motherPdg)==kLambda0 || motherPdg==kK0Short))
				{
					isFromStrangeness = kTRUE;
				}
			}
		}
	}
	
	return isFromStrangeness;
}
//_______________________________________________________________
Bool_t AliAnalysisTaskSEFT2Simulation::IsFastMcFromWeakDecay(Int_t label, AliStack *const stack){
	
	if(stack) {
		Int_t mcStackSize=stack->GetNtrack();
		if (label>=mcStackSize) return kFALSE;
		
		TParticle* particle = stack->Particle(label);
		if (!particle) return kFALSE;
		Float_t codepart = (Float_t)TMath::Abs(particle->GetPdgCode());
		Int_t motherLabel = particle->GetMother(0);
		if(motherLabel<0) return kFALSE;
		if (motherLabel>=mcStackSize) return kFALSE;
		TParticle* mother = stack->Particle(motherLabel);
		if(mother) {
			Int_t nDgh = mother->GetNDaughters();
			Float_t codemoth = (Float_t)TMath::Abs(mother->GetPdgCode());
			if(codemoth>10000000) return kFALSE;
			Int_t mfl = Int_t (codemoth / TMath::Power(10, Int_t(TMath::Log10(codemoth))));
			Double_t massMum = TDatabasePDG::Instance()->GetParticle(codemoth)->Mass();
			if(mfl>3) {return kFALSE;}
			//	if(mfl!=3) return kFALSE;
			if(codemoth>=331 && codemoth<=337){return kFALSE;}
			for(Int_t i=0;i<nDgh;i++){
				if(mother->GetLastDaughter()>=mcStackSize){return kFALSE;}
				TParticle* mcParticleDgh = stack->Particle(mother->GetLastDaughter()-i);
				if(mcParticleDgh){
					Int_t dfl = 0;
					Float_t codedgh = (Float_t)TMath::Abs(mcParticleDgh->GetPdgCode());
					if(codedgh>10000000 && codepart==13) return kTRUE; // muon decays --> Weak decays ala AliStack
					if(codedgh>10000000 && codepart!=13) return kFALSE;
					dfl = Int_t (codedgh / TMath::Power(10, Int_t(TMath::Log10(codedgh))));
					Double_t massDgh = TDatabasePDG::Instance()->GetParticle(codedgh)->Mass();
					if(codepart==11 && codemoth==13)return kTRUE;
					if(codepart==211 && mfl==dfl && ((codedgh<=3334 && codedgh>=3300) || (codedgh<=4332 && codedgh>=4334) || (codemoth<=3334 && codemoth>=3300) || (codemoth<=4334 && codemoth>=4332))){
						// double/triple strange content
						if(((codedgh<3334 && codedgh>=3300) || (codedgh<=4332 && codedgh>=4334)) && ((codemoth<3334 && codemoth>=3300) || (codemoth<=4334 && codemoth>=4332))) return kFALSE; // double strange particles decay into double strange particles
						if((codedgh<3334 && codedgh>=3300) && (codemoth<=3334))return kTRUE; // triple strange into double strange + weak decay
						return kTRUE; // double strange particles decay into single strange + weak decay
					}
					if ((dfl>=mfl))/* && codepart!=11) || (codepart==11 && mfl!=1 && dfl>=mfl && nDgh==1))*/{return	kFALSE;}
					if(massDgh>massMum){return kFALSE;}
				}
				else return kFALSE;
			}
			return kTRUE;
		}
	}
	return kFALSE;
}
