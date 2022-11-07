#include <iostream>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDTrdTracklet.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliMCParticle.h"
#include "TPDGCode.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHe3EffTree.h"



class AliAnalysisTaskHe3EffTree;   
using namespace std;           

ClassImp(AliAnalysisTaskHe3EffTree) 

AliAnalysisTaskHe3EffTree::AliAnalysisTaskHe3EffTree() 
: AliAnalysisTaskSE(), 
fESDevent(0),
mcEvent(0),
fStack(0),
fInputHandler(0),
fPIDResponse(0),
fEventCuts(),
fHistEvents(0),
fHistdEdx(0),
fOutputList(0),
fTree(0),
fTreeGen(0),
fBetheParamsHe(),
fBetheParamsT(),
fUseExternalSplines(kFALSE),
matching(0),
fMCtrue(0),
fYear(0),
tRunNumber(0),
tTrigMB(-999),			
tTrigHMV0(-999),
tTrigHMSPD(-999),
tTrigHNU(0),
tTrigHQU(0),
tMagField(0),
tTRDvalid(0),
tTRDtrigHNU(0),
tTRDtrigHQU(0),
tTRDPid(0),
tTRDnTracklets(0),
tTRDPt(0),
tTRDLayerMask(0),
tTRDSagitta(-1),
tTRDStack(0),
tTRDSector(0),
tTRDPID0(0),
tTRDPID1(0),
tTRDPID2(0),
tTRDPID3(0),
tTRDPID4(0),
tTRDPID5(0),
tnTPCcluster(0),
tnITScluster(0),
tTPCchi2(0),
tSPDFiredChips0(-999),	
tSPDFiredChips1(-999),
tSPDTracklets(-999),
tSPDCluster(-999),
tV0Multiplicity(-999),
tNTracks(-999),
tMultV0M(-999),	
tMultOfV0M(-999),			
tMultSPDTracklet(-999),	
tMultSPDCluster(-999),	
tMultRef05(-999),			
tMultRef08(-999),	
tPt(-999),				
tCharge(-999),
tY(-999),
tEta(-999),
tPhi(-999),
tP(-999),	
tPx(-999),
tPy(-999),
tPz(-999),
tE(-999),		
tKink(-1),
tTPCrefit(-1),	
tITSrefit(-1),					
tHeDEdx(-999),
tHeSigma(-999),
tTOFSignalHe(-999),
tDcaXY(-999),				
tDcaZ(-999),
tSigmaYX(-999),
tSigmaXYZ(-999),
tSigmaZ(-999),
tMCtrue(-999),			
tPrimary(-999),
tWeak(-999),
tMaterial(-999),
tHypertriton(-999),
tGenCharge(-999),
tGenPt(-999),			
tGenY(-999),
tGenPrimary(-999),			
tGenHypertriton(-999)
{
}
//_____________________________________________________________________________
AliAnalysisTaskHe3EffTree::AliAnalysisTaskHe3EffTree(const char* name) 
: AliAnalysisTaskSE(name),
fESDevent(0),
mcEvent(0),
fStack(0),
fInputHandler(0),
fPIDResponse(0),
fEventCuts(),
fHistEvents(0),
fHistdEdx(0),
fOutputList(0),
fTree(0),
fTreeGen(0),
fBetheParamsHe(),
fBetheParamsT(),
fUseExternalSplines(kFALSE),
matching(0),
fMCtrue(0),
fYear(0),
tRunNumber(0),
tTrigMB(-999),			
tTrigHMV0(-999),
tTrigHMSPD(-999),
tTrigHNU(0),
tTrigHQU(0),
tMagField(0),
tTRDvalid(0),
tTRDtrigHNU(0),
tTRDtrigHQU(0),
tTRDPid(0),
tTRDnTracklets(0),
tTRDPt(0),
tTRDLayerMask(0),
tTRDSagitta(-1),
tTRDStack(0),
tTRDSector(0),
tTRDPID0(0),
tTRDPID1(0),
tTRDPID2(0),
tTRDPID3(0),
tTRDPID4(0),
tTRDPID5(0),
tnTPCcluster(0),
tnITScluster(0),
tTPCchi2(0),
tSPDFiredChips0(-999),	
tSPDFiredChips1(-999),
tSPDTracklets(-999),
tSPDCluster(-999),
tV0Multiplicity(-999),
tNTracks(-999),
tMultV0M(-999),	
tMultOfV0M(-999),			
tMultSPDTracklet(-999),	
tMultSPDCluster(-999),	
tMultRef05(-999),			
tMultRef08(-999),	
tPt(-999),				
tCharge(-999),
tY(-999),
tEta(-999),
tPhi(-999),
tP(-999),	
tPx(-999),
tPy(-999),
tPz(-999),
tE(-999),		
tKink(-1),
tTPCrefit(-1),		
tITSrefit(-1),					
tHeDEdx(-999),
tHeSigma(-999),
tTOFSignalHe(-999),
tDcaXY(-999),				
tDcaZ(-999),
tSigmaYX(-999),
tSigmaXYZ(-999),
tSigmaZ(-999),
tMCtrue(-999),			
tPrimary(-999),
tWeak(-999),
tMaterial(-999),
tHypertriton(-999),
tGenCharge(-999),
tGenPt(-999),			
tGenY(-999),
tGenPrimary(-999),			
tGenHypertriton(-999)

{
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TList::Class());   
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskHe3EffTree::~AliAnalysisTaskHe3EffTree() {
    if(fOutputList) {
        delete fOutputList;    
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::Terminate(Option_t *){
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::UserCreateOutputObjects() {
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	if (man) {
		fInputHandler = dynamic_cast<AliESDInputHandler*> (man->GetInputEventHandler());
		if (fInputHandler)   fPIDResponse = fInputHandler->GetPIDResponse();
	}
	
	matching = new AliTRDonlineTrackMatching(); 

	fOutputList = new TList();         
	fOutputList->SetOwner(kTRUE);     
 
	fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);

	fHistEvents = new TH1F("fHistEvents", "fHistEvents", 6, 0, 6);     
	fHistEvents->GetXaxis()->SetBinLabel(1,"Events");
	fHistEvents->GetXaxis()->SetBinLabel(2,"MB");
	fHistEvents->GetXaxis()->SetBinLabel(3,"HMV0");
	fHistEvents->GetXaxis()->SetBinLabel(4,"HMSPD");
	fHistEvents->GetXaxis()->SetBinLabel(5,"HNU");
	fHistEvents->GetXaxis()->SetBinLabel(6,"HQU");

	fOutputList->Add(fHistEvents);          
	fOutputList->Add(fHistdEdx);
	
		fHistEventsTRD[0] = new TH2D("fHistEventsTRD_bfec", "fHistEventsTRDbeforeEventCuts", 9, 0, 9,10000,0,100);     
	fHistEventsTRD[1] = new TH2D("fHistEventsTRD_afec", "fHistEventsTRDafterEventCuts", 9, 0, 9,10000,0,100);  
	for (int i = 0; i < 2; i++) {   
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(1,"Events");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(2,"MB");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(3,"HNUsimMB");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(4,"HQUsimMB");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(5,"TRD");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(6,"HNU");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(7,"HNUsimTRD");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(8,"HQU");
	  fHistEventsTRD[i]->GetXaxis()->SetBinLabel(9,"HQUsimTRD");
	  fOutputList->Add(fHistEventsTRD[i]);
	}  
	
	fEventCuts.AddQAplotsToList(fOutputList);

	fTree = new TTree("treeHe","fTree");

	fTree->Branch("tRunNumber"      , &tRunNumber      , "tRunNumber/I");
	fTree->Branch("tTrigMB"         , &tTrigMB         , "tTrigMB/I");
	fTree->Branch("tTrigHMV0"       , &tTrigHMV0       , "tTrigHMV0/I");
	fTree->Branch("tTrigHMSPD"      , &tTrigHMSPD      , "tTrigHMSPD/I");
	fTree->Branch("tTrigHNU"        , &tTrigHNU        , "tTrigHNU/I");
	fTree->Branch("tTrigHQU"        , &tTrigHQU        , "tTrigHQU/I");	
	fTree->Branch("tMagField"       , &tMagField       , "tMagField/D");
	fTree->Branch("tTRDvalid"       , &tTRDvalid       , "tTRDvalid/I");
	fTree->Branch("tTRDtrigHNU"     , &tTRDtrigHNU     , "tTRDtrigHNU/I");
	fTree->Branch("tTRDtrigHQU"     , &tTRDtrigHQU     , "tTRDtrigHQU/I");	
	fTree->Branch("tTRDPid"         , &tTRDPid         , "tTRDPid/I");
	fTree->Branch("tTRDnTracklets"  , &tTRDnTracklets  , "tTRDnTracklets/I");
	fTree->Branch("tTRDPt"          , &tTRDPt          , "tTRDPt/I");
	fTree->Branch("tTRDLayerMask"   , &tTRDLayerMask   , "tTRDLayerMask/I");
	fTree->Branch("tTRDSagitta"     , &tTRDSagitta     , "tTRDSagitta/D");	
	fTree->Branch("tTRDStack"     , &tTRDStack     , "tTRDStack/I");	
	fTree->Branch("tTRDSector"     , &tTRDSector     , "tTRDSector/I");	
	fTree->Branch("tTRDPID0"     , &tTRDPID0     , "tTRDPID0/i");	
	fTree->Branch("tTRDPID1"     , &tTRDPID1     , "tTRDPID1/i");	
	fTree->Branch("tTRDPID2"     , &tTRDPID2     , "tTRDPID2/i");	
	fTree->Branch("tTRDPID3"     , &tTRDPID3     , "tTRDPID3/i");	
	fTree->Branch("tTRDPID4"     , &tTRDPID4     , "tTRDPID4/i");
	fTree->Branch("tTRDPID5"     , &tTRDPID5     , "tTRDPID5/i");	
	
	fTree->Branch("tnTPCcluster"    , &tnTPCcluster    , "tnTPCcluster/D");
	fTree->Branch("tnITScluster"    , &tnITScluster    , "tnITScluster/D");
	fTree->Branch("tTPCchi2"        , &tTPCchi2        , "tTPCchi2/D");	
	fTree->Branch("tSPDFiredChips0" , &tSPDFiredChips0 , "tSPDFiredChips0/D");
	fTree->Branch("tSPDFiredChips1" , &tSPDFiredChips1 , "tSPDFiredChips1/D");
	fTree->Branch("tSPDTracklets"   , &tSPDTracklets   , "tSPDTracklets/D");
	fTree->Branch("tSPDCluster"     , &tSPDCluster     , "tSPDCluster/D");
	fTree->Branch("tV0Multiplicity" , &tV0Multiplicity , "tV0Multiplicity/D");
	fTree->Branch("tNTracks"        , &tNTracks        , "tNTracks/D");
	fTree->Branch("tMultV0M"        , &tMultV0M        , "tMultV0M/D");
	fTree->Branch("tMultOfV0M"      , &tMultOfV0M      , "tMultOfV0M/D");
	fTree->Branch("tMultSPDTracklet", &tMultSPDTracklet, "tMultSPDTracklet/D");
	fTree->Branch("tMultSPDCluster" , &tMultSPDCluster , "tMultSPDCluster/D");
	fTree->Branch("tMultRef05"      , &tMultRef05      , "tMultRef05/D");
	fTree->Branch("tMultRef08"      , &tMultRef08      , "tMultRef08/D");
	fTree->Branch("tPt"             , &tPt             , "tPt/D");
	fTree->Branch("tCharge"         , &tCharge         , "tCharge/D");
	fTree->Branch("tY"              , &tY              , "tY/D");
	fTree->Branch("tEta"            , &tEta            , "tEta/D");
	fTree->Branch("tPhi"            , &tPhi            , "tPhi/D");
	fTree->Branch("tP"              , &tP              , "tP/D");
	fTree->Branch("tPx"             , &tPx             , "tPx/D");
	fTree->Branch("tPy"             , &tPy             , "tPy/D");
	fTree->Branch("tPz"             , &tPz             , "tPz/D");
	fTree->Branch("tE"              , &tE              , "tE/D");
	fTree->Branch("tKink"           , &tKink           , "tKink/I");
	fTree->Branch("tTPCrefit"       , &tTPCrefit       , "tTPCrefit/I");
	fTree->Branch("tITSrefit"       , &tITSrefit       , "tITSrefit/I");
	fTree->Branch("tHeDEdx"         , &tHeDEdx         , "tHeDEdx/D");
	fTree->Branch("tHeSigma"        , &tHeSigma        , "tHeSigma/D");
	fTree->Branch("tTOFSignalHe"    , &tTOFSignalHe    , "tTOFSignalHe/D");
	fTree->Branch("tDcaXY"          , &tDcaXY          , "tDcaXY/D");
	fTree->Branch("tDcaZ"           , &tDcaZ           , "tDcaZ/D");
	fTree->Branch("tSigmaYX"        , &tSigmaYX        , "tSigmaYX/D");
	fTree->Branch("tSigmaXYZ"       , &tSigmaXYZ       , "tSigmaXYZ/D");
	fTree->Branch("tSigmaZ"         , &tSigmaZ         , "tSigmaZ/D");
	fTree->Branch("tMCtrue"         , &tMCtrue         , "tMCtrue/I");
	fTree->Branch("tPrimary"        , &tPrimary        , "tPrimary/I");
	fTree->Branch("tWeak"           , &tWeak           , "tWeak/I");
	fTree->Branch("tMaterial"       , &tMaterial       , "tMaterial/I");
	fTree->Branch("tHypertriton"    , &tHypertriton    , "tHypertriton/I");

	fTreeGen = new TTree("treeGenHe","fTreeGen");
	fTreeGen->Branch("tTrigMB"        , &tTrigMB        , "tTrigMB/I");
	fTreeGen->Branch("tTrigHMV0"      , &tTrigHMV0      , "tTrigHMV0/I");
	fTreeGen->Branch("tTrigHMSPD"     , &tTrigHMSPD     , "tTrigHMSPD/I");
	fTreeGen->Branch("tTrigHNU"       , &tTrigHNU       , "tTrigHNU/I");
	fTreeGen->Branch("tTrigHQU"       , &tTrigHQU       , "tTrigHQU/I");
	fTreeGen->Branch("tMagField"      , &tMagField      , "tMagField/D");
	fTreeGen->Branch("tGenCharge"     , &tGenCharge     , "tGenCharge/D");
	fTreeGen->Branch("tGenPt"         , &tGenPt         , "tGenPt/D");
	fTreeGen->Branch("tGenY"          , &tGenY          , "tGenY/D");
	fTreeGen->Branch("tGenPrimary"    , &tGenPrimary    , "tGenPrimary/I");
	fTreeGen->Branch("tGenHypertriton", &tGenHypertriton, "tGenHypertriton/I");

	PostData(1, fOutputList);
	PostData(2, fTree);
	PostData(3, fTreeGen);              
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::UserExec(Option_t *) {
  // MC
	fMCtrue = kTRUE;
	AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>
		(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	if (!mcEventHandler) fMCtrue = kFALSE;
	if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
	if (!mcEvent && fMCtrue) return;
	if (fMCtrue) {
		fStack = mcEvent->Stack();
		if (!fStack) return;
	}

	fESDevent = dynamic_cast<AliESDEvent*>(InputEvent()); 
	fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);


  Int_t runNumber = fESDevent->GetRunNumber();
  if (!fUseExternalSplines) SetBetheBlochParams(runNumber);


	AliCDBManager *cdbMgr = AliCDBManager::Instance();
	if (fMCtrue) {
		cdbMgr->SetDefaultStorage("MC","Full");
	}
	else {
		cdbMgr->SetDefaultStorage (Form("alien://Folder=/alice/data/%d/OCDB", fYear));
	}
	cdbMgr->SetRun(runNumber);
	AliGeomManager::LoadGeometry();
	
	AliESDtrackCuts trackCutsV0("AlitrackCutsV0", "AlitrackCutsV0");
	trackCutsV0.SetEtaRange(-0.9,0.9);
	trackCutsV0.SetAcceptKinkDaughters(kTRUE);
	trackCutsV0.SetRequireTPCRefit(kFALSE);
	trackCutsV0.SetMaxChi2PerClusterTPC(8);
	trackCutsV0.SetMinNClustersTPC(40);

	AliMultSelection *MultSelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
	if (MultSelection) {
		tMultV0M = MultSelection->GetMultiplicityPercentile("V0M");
		tMultOfV0M = MultSelection->GetMultiplicityPercentile("OnlineV0M");
		tMultSPDTracklet = MultSelection->GetMultiplicityPercentile("SPDClusters");
		tMultSPDCluster = MultSelection->GetMultiplicityPercentile("SPDTracklets");
		tMultRef05 = MultSelection->GetMultiplicityPercentile("RefMult05");
		tMultRef08 = MultSelection->GetMultiplicityPercentile("RefMult08");
	}
	
	TRDnEvents();
	
	if(!fEventCuts.AcceptEvent(fESDevent)) return;
	//******************************
	//*   get trigger information  *
	//******************************
	Bool_t MB = kFALSE;		// minimum bias
	Bool_t HMV0 = kFALSE;	// high multiplicity V0
	Bool_t HMSPD = kFALSE;	// high multiplicity SPD
	Int_t HNU = 0;			// TRD nuclei
	Int_t HQU = 0;			// TRD quarkonia

	if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) MB = kTRUE;
	if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0) HMV0 = kTRUE;
	if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) HMSPD = kTRUE;
	
	Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
	if (!fMCtrue){
		// Data: get TRD trigger information from trigger classes 
		TString classes = fESDevent->GetFiredTriggerClasses();   
		if (classes.Contains("HNU")) HNU = 1;
		if (classes.Contains("HQU")) HQU = 1; 
		
	} else {
		// MC: simulate TRD trigger
		Bool_t primaryHeHNU = kFALSE, primaryHeHQU = kFALSE;

		if (nTrdTracks > 0) {
			for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
				AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
				if (!trdTrack) continue;
				
				Int_t label = trdTrack->GetLabel();
				AliMCParticle *particle = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label))->Particle());
			
				// simulate HNU
				if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) || 
					(trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {	
					HNU = 1;
					if (TMath::Abs(particle->PdgCode()) == 1000020030) {
						if (mcEvent->IsPhysicalPrimary(TMath::Abs(label))) primaryHeHNU = kTRUE;
					}
				}
				// simulate HQU
				if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
					trdTrack->GetPID() >= 130 && trdTrack->GetNTracklets() >= 5 && (trdTrack->GetLayerMask() & 1) ){	
					Double_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
					if (sag < 0.2 && sag > -0.2) {
						HQU = 1;
						if (TMath::Abs(particle->PdgCode()) == 1000020030) {
							if (mcEvent->IsPhysicalPrimary(TMath::Abs(label))) primaryHeHQU = kTRUE;
						}
					}
				}
			}
		if (primaryHeHNU) HNU = 2;
		if (primaryHeHQU) HQU = 2;
		}
	}
	// additional multiplicity information
	tNTracks = fESDevent->GetNumberOfTracks();
	
	AliESDVZERO *vzero = fESDevent->GetVZEROData();
	tV0Multiplicity = 0;
	for (Int_t ii = 0; ii < 64; ii++){
		tV0Multiplicity += vzero->GetMultiplicity(ii);
	}
	
	AliMultiplicity *multSPD = fESDevent->GetMultiplicity();
	tSPDCluster	= multSPD->GetNumberOfSPDClusters();
	tSPDTracklets = multSPD->GetNumberOfTracklets();
	tSPDFiredChips0 = multSPD->GetNumberOfFiredChips(0);
	tSPDFiredChips1 = multSPD->GetNumberOfFiredChips(1);


	
	tMagField = fESDevent->GetMagneticField();

	// fill histogram
	fHistEvents->Fill(0);
	if (MB) fHistEvents->Fill(1);
	if (HMV0) fHistEvents->Fill(2);
	if (HMSPD) fHistEvents->Fill(3);
	if (HNU) fHistEvents->Fill(4);
	if (HQU) fHistEvents->Fill(5);
	
		
	//***************************************
	//* loop over all tracks, identify He3  *
	//***************************************
	Int_t iTracks(fESDevent->GetNumberOfTracks());         
	for(Int_t i(0); i < iTracks; i++) {                
		AliESDtrack* track = static_cast<AliESDtrack*>(fESDevent->GetTrack(i));        
		if(!track || !trackCutsV0.AcceptTrack(track)) continue;
		fHistdEdx->Fill(track->GetInnerParam()->GetP() * track->GetSign(), track->GetTPCsignal());   
		
                 
		if (track->GetTPCsignal() > 1500 || track->GetInnerParam()->GetP() > 5) continue;
		tCharge = 2 * track->GetSign();
		if (TMath::Abs(Bethe(*track, AliPID::ParticleMass(AliPID::kHe3), TMath::Abs(tCharge), fBetheParamsHe)) > 5) {
			continue;
			}
		        
		tRunNumber = runNumber;
		tTrigMB = MB;		
		tTrigHMV0 = HMV0;
		tTrigHMSPD = HMSPD;
		tTrigHNU = HNU;
		tTrigHQU = HQU;

		tnTPCcluster = track->GetTPCNcls();
		tTPCchi2 = track->GetTPCchi2() / (Double_t) track->GetTPCclusters(0);
		tnITScluster = track->GetITSNcls();
		
		Double_t momvect[3];
		track->PxPyPz(momvect);
		TLorentzVector He3Vector(momvect[0],momvect[1],momvect[2],0);
		He3Vector.SetVect(TMath::Abs(tCharge) * He3Vector.Vect());
		tE = TMath::Sqrt(AliPID::ParticleMass(AliPID::kHe3) * AliPID::ParticleMass(AliPID::kHe3) + He3Vector.Vect().Mag2());
 		He3Vector.SetE(tE);

		tPx = momvect[0];
		tPy = momvect[1];
		tPz = momvect[2];
		tPt = He3Vector.Pt();
		tY = He3Vector.Rapidity();
		tEta = track->Eta();
		tPhi = track->Phi();
		
		tKink = track->GetKinkIndex(0) > 0;
		tTPCrefit = (track->GetStatus() & AliESDtrack::kTPCrefit) != 0;
		tITSrefit = (track->GetStatus() & AliESDtrack::kITSrefit) != 0;
		
		tP = track->GetInnerParam()->GetP();
		tHeDEdx =  track->GetTPCsignal();
		tHeSigma = Bethe(*track, AliPID::ParticleMass(AliPID::kHe3), TMath::Abs(tCharge), fBetheParamsHe);
		tTOFSignalHe = GetTOFSignalHe3(*track, tP);

		Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
		track->GetImpactParameters(dca, cov);
		tDcaXY = dca[0];
		tDcaZ = dca[1];
		tSigmaYX = cov[0];
		tSigmaXYZ = cov[1];
		tSigmaZ = cov[2];

		TRDtrack(track);

		if (fMCtrue) {
		    Int_t label = track->GetLabel();
			AliMCParticle *particle = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(label))->Particle());	
			Int_t labelMother = mcEvent->GetLabelOfParticleMother(TMath::Abs(label));
			AliMCParticle *particleMother = new AliMCParticle(mcEvent->GetTrack(TMath::Abs(labelMother))->Particle()); 
			tMCtrue = TMath::Abs(particle->PdgCode()) == 1000020030;	
			tPrimary = mcEvent->IsPhysicalPrimary(TMath::Abs(label));
			tWeak = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(label)); 
			tMaterial = mcEvent->IsSecondaryFromMaterial(TMath::Abs(label));
			tHypertriton = TMath::Abs(particleMother->PdgCode()) == 1010010030;		
		}
		fTree->Fill();
	}
	if (fMCtrue) ProcessMCParticles();
	PostData(1, fOutputList);  
	PostData(2, fTree);                             
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::ProcessMCParticles() {
	// fills tree with generated He3
	Int_t nMCprimaries = mcEvent->GetNumberOfPrimaries();

	for (Int_t istack = 0; istack < mcEvent->GetNumberOfTracks(); istack++) {
		AliMCParticle* particle =  (AliMCParticle*)(mcEvent->GetTrack(istack));
		if(!particle){ continue; }
		Long_t pdgCode = particle->PdgCode();
		if (TMath::Abs(pdgCode) == 1000020030) {
			tGenPt = particle->Pt();
			tGenY = particle->Y();
			tGenCharge = -2;
			if (pdgCode == 1000020030) 
			tGenCharge = 2;

			Int_t label = particle->Label();
			tGenPrimary = ((label >= 0  && label <= nMCprimaries) && (particle->IsPhysicalPrimary()));
			tGenHypertriton = mcEvent->IsSecondaryFromWeakDecay(TMath::Abs(istack));	
		fTreeGen->Fill();
		}
	}
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskHe3EffTree::GetInvPtDevFromBC(Int_t b, Int_t c) {
	//returns d(1/Pt) in c/GeV 
	//in case of no gtu simulation -> return maximum 0.5
	if(b==0 && c==0) return 0.5;
	Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
	tmp += (c & 0xfff);
	Double_t invPtDev = tmp * 0.000001;
	return invPtDev;
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::TRDnEvents() {
	//*   get trigger information  *
	Bool_t MB = kFALSE;		// minimum bias
	Int_t HNU = 0;			// TRD nuclei
	Int_t HQU = 0;			// TRD quarkonia
  Int_t HJT = 0;
  
  Bool_t TRD = kFALSE;		// minimum bias
	Int_t HNUTRD = 0;			// TRD nuclei
	Int_t HQUTRD = 0;			// TRD quarkonia
  Int_t HJTTRD = 0;
 	Int_t HNUsimTRD = 0;			// TRD nuclei
	Int_t HQUsimTRD = 0;			// TRD quarkonia
  Int_t HJTsimTRD = 0; 
  
	if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) MB = kTRUE;
	if (fInputHandler->IsEventSelected() & AliVEvent::kTRD) TRD = kTRUE;
	
	if (TRD) {
		TString classes = fESDevent->GetFiredTriggerClasses();   
		if (classes.Contains("HNU")) HNUTRD = 1;
		if (classes.Contains("HQU")) HQUTRD = 1; 
		if (classes.Contains("HJT")) HJTTRD = 1; 
	} 
	
	  Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
		// MC: simulate TRD trigger

		
		Int_t nTracks[90] = { 0 }; // stack-wise counted number of tracks above pt threshold

    
		if (nTrdTracks > 0) {
			for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
				AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
				if (!trdTrack) continue;
	
				// simulate HNU
				if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) || 
					(trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {	
					if (MB) HNU = 1;
					if (TRD) HNUsimTRD = 1;
				}
				// simulate HQU
				if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
					trdTrack->GetPID() >= 130 && trdTrack->GetNTracklets() >= 5 && (trdTrack->GetLayerMask() & 1) ){	
					Double_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
					if (sag < 0.2 && sag > -0.2) {
						if (MB) HQU = 1;
						if (TRD) HQUsimTRD = 1;
					}
				}
      }
		}
  
	// fill histogram
	//fHistEventsTRD[0]->Fill(0,tMultV0M);
	if (MB) fHistEventsTRD[0]->Fill(1,tMultV0M);
	if (HNU) fHistEventsTRD[0]->Fill(2,tMultV0M);
	if (HQU) fHistEventsTRD[0]->Fill(3,tMultV0M);
	if (TRD) fHistEventsTRD[0]->Fill(4,tMultV0M);
	if (HNUTRD) fHistEventsTRD[0]->Fill(5,tMultV0M);
	if (HQUTRD) fHistEventsTRD[0]->Fill(7,tMultV0M);
	if (HNUsimTRD) fHistEventsTRD[0]->Fill(6,tMultV0M);
	if (HQUsimTRD) fHistEventsTRD[0]->Fill(8,tMultV0M);
	
	if (fEventCuts.AcceptEvent(fESDevent)) {
		  //fHistEventsTRD[1]->Fill(0,tMultV0M);
	  if (MB) fHistEventsTRD[1]->Fill(1,tMultV0M);
	  if (HNU) fHistEventsTRD[1]->Fill(2,tMultV0M);
	  if (HQU) fHistEventsTRD[1]->Fill(3,tMultV0M);
	  if (TRD) fHistEventsTRD[1]->Fill(4,tMultV0M);
	  if (HNUTRD) fHistEventsTRD[1]->Fill(5,tMultV0M);
	  if (HQUTRD) fHistEventsTRD[1]->Fill(7,tMultV0M);
	  if (HNUsimTRD) fHistEventsTRD[1]->Fill(6,tMultV0M);
	  if (HQUsimTRD) fHistEventsTRD[1]->Fill(8,tMultV0M);
	}
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskHe3EffTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
	/// Calculates number of sigma deviation from expected dE/dx in TPC
	/// \param track particle track
	/// \param mass hypothesis of particle
	/// \param charge particle charge hypothesis
	/// \param params Parameters of Aleph parametrization of Bethe Energy-loss
	Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
	Double_t sigma = expected*params[5];
	if (TMath::IsNaN(expected)) return -999;
	return (track.GetTPCsignal() - expected) / sigma;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskHe3EffTree::GetTOFSignalHe3(AliESDtrack& trackHe, Double_t tP) {
	// returns the mass calculated from TOF Signal
	Double_t mass = 0, time = -1, beta = 0, gamma = 0, length = 0, time0 = 0;
	length = trackHe.GetIntegratedLength();
	time0 = fPIDResponse->GetTOFResponse().GetStartTime(trackHe.P());
    time = trackHe.GetTOFsignal() - time0;
	if (time > 0 && length > 0) {
	beta = length / (2.99792457999999984e-02 * time);
	gamma = 1 / TMath::Sqrt(1 - beta*beta);
	mass = (tP / TMath::Sqrt(gamma*gamma - 1));
	return mass*mass;
	}
	return -1;
}
//_____________________________________________________________________________
void AliAnalysisTaskHe3EffTree::SetBetheBlochParams(Int_t runNumber) {
	// set Bethe-Bloch parameter
	if (runNumber >= 252235 && runNumber <= 265589) { // 2016 pp data
		fYear = 2016;
    // He3 2016/2018 pass2
    fBetheParamsHe[0] = 4.20995;
    fBetheParamsHe[1] = 10.5007;
    fBetheParamsHe[2] = -0.895979;
    fBetheParamsHe[3] = 2.01748;
    fBetheParamsHe[4] = 0.0798937;
    fBetheParamsHe[5] = 0.06;

    // Triton 2016/2018 pass2
    fBetheParamsT[0] = 12.0774;
    fBetheParamsT[1] = 5.70345;
    fBetheParamsT[2] = 4.764;
    fBetheParamsT[3] = 1.94198;
    fBetheParamsT[4] = -3.03895;
    fBetheParamsT[5] = 0.07;
	}
	if (runNumber > 265589 && runNumber <= 267166) { // 2016 p-Pb data
		fYear = 2016;
		// He3
		fBetheParamsHe[0] = 0.715489;
		fBetheParamsHe[1] = 59.5463;
		fBetheParamsHe[2] = 4.44487e-12;
		fBetheParamsHe[3] = 2.69874;
		fBetheParamsHe[4] = 24.063;
		fBetheParamsHe[5] = 0.04725;
		// Triton
		fBetheParamsT[0] = 0.223948;
		fBetheParamsT[1] = 180.564;
		fBetheParamsT[2] = -3.03884e-10;
		fBetheParamsT[3] = 2.30095;
		fBetheParamsT[4] = 34.2269;
		fBetheParamsT[5] = 0.06517;	
	} 	
	if (runNumber >= 270581 && runNumber <= 282704) { // 2017 pp data
		fYear = 2017;
    // He3 2017 pass2
    fBetheParamsHe[0] = 1.65042;
    fBetheParamsHe[1] = 25.9254;
    fBetheParamsHe[2] = 0.00600469;
    fBetheParamsHe[3] = 2.73841;
    fBetheParamsHe[4] = 10.8988;
    fBetheParamsHe[5] = 0.06;

    // Triton 2017 pass2
    fBetheParamsT[0] = 2.82837;
    fBetheParamsT[1] = 15.4278;
    fBetheParamsT[2] = 1.03545;
    fBetheParamsT[3] = 2.2757;
    fBetheParamsT[4] = 2.7525;
    fBetheParamsT[5] = 0.06;
	}
	if (runNumber >= 285009 && runNumber <= 294925) { // 2018 pp data
		fYear = 2018;
    // He3 2016/2018 pass2
    fBetheParamsHe[0] = 4.20995;
    fBetheParamsHe[1] = 10.5007;
    fBetheParamsHe[2] = -0.895979;
    fBetheParamsHe[3] = 2.01748;
    fBetheParamsHe[4] = 0.0798937;
    fBetheParamsHe[5] = 0.06;

    // Triton 2016/2018 pass2
    fBetheParamsT[0] = 12.0774;
    fBetheParamsT[1] = 5.70345;
    fBetheParamsT[2] = 4.764;
    fBetheParamsT[3] = 1.94198;
    fBetheParamsT[4] = -3.03895;
    fBetheParamsT[5] = 0.07;
	}
}
//_____________________________________________________________________________

Double_t AliAnalysisTaskHe3EffTree::TRDtrack(AliESDtrack* esdTrack) {

		tTRDvalid = 0;
		tTRDtrigHNU = 0;
		tTRDtrigHQU = 0;
		tTRDPid = 0;
		tTRDnTracklets = 0;
		tTRDPt = 0;
		tTRDLayerMask = 0;
		tTRDSagitta = -1;
		tTRDPID0 = 0;
		tTRDStack = 0;
		tTRDSector = 0;
		tTRDPID1 = 0;
		tTRDPID2 = 0;
		tTRDPID3 = 0;
		tTRDPID4 = 0;
		tTRDPID5 = 0;
	
    if(!esdTrack) {
        return 0;
    }

    if(fESDevent->GetNumberOfTrdTracks() == 0) {
        return 0;
    }
    
    AliESDTrdTrack* bestGtuTrack = 0x0;
       
    Double_t esdPt = esdTrack->GetSignedPt();
    Double_t mag = fESDevent->GetMagneticField();
    Double_t currentMatch = 0;
    Double_t bestMatch = 0;

    for (Int_t i = 0; i < fESDevent->GetNumberOfTrdTracks(); i++) {

        AliESDTrdTrack* gtuTrack= fESDevent->GetTrdTrack ( i );
        Double_t gtuPt = gtuTrack->Pt();
        if (mag > 0.) gtuPt = gtuPt * (-1.0);

        Double_t ydist;
        Double_t zdist;

        if (matching->EstimateTrackDistance(esdTrack, gtuTrack, mag, &ydist, &zdist) == 0) {
        	currentMatch = matching->RateTrackMatch(ydist, zdist, esdPt, gtuPt);
				}
				
        if (currentMatch > bestMatch) {
            bestMatch = currentMatch;
            bestGtuTrack = gtuTrack;
        }
    }
    
    if (!bestGtuTrack) {
    	return 0;
    }
    
    tTRDvalid = 1;
		tTRDPid = bestGtuTrack->GetPID();
		tTRDnTracklets = bestGtuTrack->GetNTracklets();
		tTRDPt = (TMath::Abs(bestGtuTrack->GetPt()));
		tTRDLayerMask =  bestGtuTrack->GetLayerMask();
		tTRDSagitta = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());	
		tTRDStack = bestGtuTrack->GetStack();
		tTRDSector = bestGtuTrack->GetSector();

		AliESDTrdTracklet *trk = 0;
		trk = bestGtuTrack->GetTracklet(0);
		if (trk) tTRDPID0 = trk->GetTrackletWord();
		trk = 0;
		trk = bestGtuTrack->GetTracklet(1);
		if (trk) tTRDPID1 = trk->GetTrackletWord();
		trk = 0;
		trk = bestGtuTrack->GetTracklet(2);
		if (trk) tTRDPID2 = trk->GetTrackletWord();
		trk = 0;
		trk = bestGtuTrack->GetTracklet(3);
		if (trk) tTRDPID3 = trk->GetTrackletWord();
		trk = 0;
		trk = bestGtuTrack->GetTracklet(4);
		if (trk) tTRDPID4 = trk->GetTrackletWord();
		trk = 0;
		trk = bestGtuTrack->GetTracklet(5);
		if (trk) tTRDPID5 = trk->GetTrackletWord();		

			
		if((bestGtuTrack->GetPID() >= 255 && bestGtuTrack->GetNTracklets() == 4) || 
			(bestGtuTrack->GetPID() >= 235 && bestGtuTrack->GetNTracklets() > 4)) {
						tTRDtrigHNU = 1;
		}		

			if (TMath::Abs(bestGtuTrack->GetPt()) >= 256 &&
				bestGtuTrack->GetPID() >= 130 && bestGtuTrack->GetNTracklets() >= 5 && (bestGtuTrack->GetLayerMask() & 1) ){	
				Double_t sag = GetInvPtDevFromBC(bestGtuTrack->GetB(), bestGtuTrack->GetC());
					if (sag < 0.2 && sag > -0.2) {
							tTRDtrigHQU = 1;
				 	}
			}	

    return bestMatch;

}
//_____________________________________________________________________________
//_____________________________________________________________________________
