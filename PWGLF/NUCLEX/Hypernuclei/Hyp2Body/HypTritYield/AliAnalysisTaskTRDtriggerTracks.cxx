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
#include "AliESDpid.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliTPCPIDResponse.h"
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
#include "AliEventCuts.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliTOFPIDResponse.h"
#include "AliAnalysisTaskTRDtriggerTracks.h"

class AliAnalysisTaskTRDtriggerTracks;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskTRDtriggerTracks) // classimp: necessary for root

AliAnalysisTaskTRDtriggerTracks::AliAnalysisTaskTRDtriggerTracks() : AliAnalysisTaskSE(), 
fESD(0),    
fPIDResponse(0), 
fOutputList(0),   
fEventCuts(0),   	
fInputHandler(0),   				
fTree(0),      
tTrigMB(0),			
tTrigHMV0(0),
tTrigHMSPD(0),
tTrigHNU(0),
tTrigHQU(0),    
tCharge(0),
tPt(0),
tPx(0),
tPy(0),
tPz(0),			
tY(0),
tP(0),			
tTPCDEdx(0),
tTOFSignal(0),
tDcaXY(0),			
tDcaZ(0),			
tTRDtrigHNU(0),
tTRDtrigHQU(0),			
tTRDPid(0),
tTRDnTracklets(0),
tTRDPt(0),
tTRDLayerMask(0),
tTRDSagitta(0),
tMultV0M(0),
tMultSPDCluster(0),
tSPDTracklets(0),
histMultV0HNU(0), 
histMultSPDHNU(0), 
histMultRefHNU(0),
histMultV0HQU(0),
histMultSPDHQU(0),
histMultRefHQU(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}
//_____________________________________________________________________________
AliAnalysisTaskTRDtriggerTracks::AliAnalysisTaskTRDtriggerTracks(const char* name) : AliAnalysisTaskSE(name),
fESD(0),    
fPIDResponse(0), 
fOutputList(0),   
fEventCuts(0),   	
fInputHandler(0),   				
fTree(0),   
tTrigMB(0),			
tTrigHMV0(0),
tTrigHMSPD(0),
tTrigHNU(0),
tTrigHQU(0),      
tCharge(0),
tPt(0),
tPx(0),
tPy(0),
tPz(0),		
tY(0),
tP(0),			
tTPCDEdx(0),
tTOFSignal(0),
tDcaXY(0),			
tDcaZ(0),			
tTRDtrigHNU(0),
tTRDtrigHQU(0),			
tTRDPid(0),
tTRDnTracklets(0),
tTRDPt(0),
tTRDLayerMask(0),
tTRDSagitta(0),
tMultV0M(0),
tMultSPDCluster(0),
tSPDTracklets(0),
histMultV0HNU(0), 
histMultSPDHNU(0), 
histMultRefHNU(0),
histMultV0HQU(0),
histMultSPDHQU(0),
histMultRefHQU(0)
{
    DefineInput(0, TChain::Class());   
    DefineOutput(1, TList::Class());   
    DefineOutput(2, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskTRDtriggerTracks::~AliAnalysisTaskTRDtriggerTracks()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________
void AliAnalysisTaskTRDtriggerTracks::UserCreateOutputObjects()
{
  fInputHandler = dynamic_cast<AliESDInputHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPIDResponse = fInputHandler->GetPIDResponse();

  fOutputList = new TList();        
  fOutputList->SetOwner(kTRUE);       
	// multiplicity bins
	const Int_t nBins = 13;
	Double_t edges[nBins + 1] = {0.0, 0.01, 0.1, 0.5, 1.0, 5., 10., 15., 20., 30.,40.,50.,70.,100.};
	const Int_t nBins2 = 10;
	Double_t edges2[nBins2 + 1] = {0.0, 1., 5., 10., 15.,20.,30.,40.,50.,70.,100.};

	// create histograms
	histMultV0HNU = new TH1F("histMultV0HNU","HNU - V0M estimator;MultiplicityPercentile;Counts", nBins, edges);
	histMultSPDHNU = new TH1F("histMultSPDHNU","HNU - SPDTracklet estimator;MultiplicityPercentile;Counts",nBins2, edges2);
	histMultRefHNU = new TH1F("histMultRefHNU","HNU - Ref08 estimator;MultiplicityPercentile;Counts",nBins2, edges2);
	histMultV0HQU = new TH1F("histMultV0HQU","HQU - V0M estimator;MultiplicityPercentile;Counts", nBins, edges);
	histMultSPDHQU = new TH1F("histMultSPDHQU","HQU - SPDTracklet estimator;MultiplicityPercentile;Counts",nBins2, edges2);
	histMultRefHQU = new TH1F("histMultRefHQU","HQU - Ref08 estimator;MultiplicityPercentile;Counts",nBins2, edges2);
	fOutputList->Add(histMultV0HNU);
	fOutputList->Add(histMultSPDHNU);
	fOutputList->Add(histMultRefHNU);
	fOutputList->Add(histMultV0HQU);
	fOutputList->Add(histMultSPDHQU);
	fOutputList->Add(histMultRefHQU);                      
    
  fTree = new TTree("tree", "fTree");

	fTree->Branch("tTrigMB"         , &tTrigMB         , "tTrigMB/I");
	fTree->Branch("tTrigHMV0"       , &tTrigHMV0       , "tTrigHMV0/I");
	fTree->Branch("tTrigHMSPD"      , &tTrigHMSPD      , "tTrigHMSPD/I");
	fTree->Branch("tTrigHNU"        , &tTrigHNU        , "tTrigHNU/I");
	fTree->Branch("tTrigHQU"        , &tTrigHQU        , "tTrigHQU/I");
          
	fTree->Branch("tCharge"         , &tCharge         , "tCharge/F");  
	fTree->Branch("tPt"             , &tPt             , "tPt/F");  
	fTree->Branch("tPx"             , &tPx             , "tPx/F");  
	fTree->Branch("tPy"             , &tPy             , "tPy/F");  
	fTree->Branch("tPz"             , &tPz             , "tPz/F");  
  fTree->Branch("tY"              , &tY              , "tY/F");  
  fTree->Branch("tP"              , &tP              , "tP/F");  
  fTree->Branch("tTPCDEdx"        , &tTPCDEdx        , "tTPCDEdx/F");
  fTree->Branch("tTOFSignal"      , &tTOFSignal      , "tTOFSignal/F");
	fTree->Branch("tDcaXY"          , &tDcaXY          , "tDcaXY/F");
	fTree->Branch("tDcaZ"           , &tDcaZ           , "tDcaZ/F");	

	fTree->Branch("tTRDtrigHNU"     , &tTRDtrigHNU     , "tTRDtrigHNU/I");
	fTree->Branch("tTRDtrigHQU"     , &tTRDtrigHQU     , "tTRDtrigHQU/I");
	
	fTree->Branch("tTRDPid"         , &tTRDPid         , "tTRDPid/I");
	fTree->Branch("tTRDnTracklets"  , &tTRDnTracklets  , "tTRDnTracklets/I");
	fTree->Branch("tTRDPt"          , &tTRDPt          , "tTRDPt/I");
	fTree->Branch("tTRDLayerMask"   , &tTRDLayerMask   , "tTRDLayerMask/I");
	fTree->Branch("tTRDSagitta"     , &tTRDSagitta     , "tTRDSagitta/F");
	
	fTree->Branch("tMultV0M"        , &tMultV0M        , "tMultV0M/F");
	fTree->Branch("tMultSPDCluster" , &tMultSPDCluster , "tMultSPDCluster/F");
	fTree->Branch("tSPDTracklets"   , &tSPDTracklets   , "tSPDTracklets/F");
	
	PostData(1, fOutputList);
	PostData(2, fTree);
}
//_____________________________________________________________________________
void AliAnalysisTaskTRDtriggerTracks::UserExec(Option_t *)
{
 
    fESD = dynamic_cast<AliESDEvent*>(InputEvent()); 
    if(!fESD) return;    
  
  AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
	if (!MultSelection) return;
	
		tMultV0M = MultSelection->GetMultiplicityPercentile("V0M");
		tMultSPDCluster = MultSelection->GetMultiplicityPercentile("SPDTracklets");
		Double_t tMultRef08 = MultSelection->GetMultiplicityPercentile("RefMult08");

		TString classes = fESD->GetFiredTriggerClasses();                      
		if (classes.Contains("HNU")) {
			histMultV0HNU->Fill(tMultV0M);
			histMultSPDHNU->Fill(tMultSPDCluster);
			histMultRefHNU->Fill(tMultRef08);		
		}
		if (classes.Contains("HQU")) {
			histMultV0HQU->Fill(tMultV0M);
			histMultSPDHQU->Fill(tMultSPDCluster);
			histMultRefHQU->Fill(tMultRef08);			
		}
	
	tTrigMB = 0;
	tTrigHMV0 = 0;
	tTrigHMSPD = 0;
	tTrigHNU = 0;
	tTrigHQU = 0; 
			
	if (fInputHandler->IsEventSelected() & AliVEvent::kINT7) tTrigMB = 1;
	if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0) tTrigHMV0 = 1;
	if (fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD) tTrigHMSPD = 1;
	if (classes.Contains("HNU")) tTrigHNU = 1;
	if (classes.Contains("HQU")) tTrigHQU = 1; 
		
	AliMultiplicity *multSPD = fESD->GetMultiplicity();
	tSPDTracklets = multSPD->GetNumberOfTracklets();
	
	Int_t nTrdTracks = fESD->GetNumberOfTrdTracks();
	if (nTrdTracks == 0) return;
	
	Int_t nTracks = fESD->GetNumberOfTracks();
	if (nTracks == 0) return;
	
	Int_t fYear = 2018;
  Int_t runNumber = fESD->GetRunNumber();
	if (runNumber >= 252235 && runNumber <= 267166) fYear = 2016;
	if (runNumber >= 270581 && runNumber <= 282704) fYear = 2017;
		
	AliCDBManager *cdbMgr = AliCDBManager::Instance();
	cdbMgr->SetDefaultStorage (Form("alien://Folder=/alice/data/%d/OCDB", fYear));
	cdbMgr->SetRun(runNumber);
	AliGeomManager::LoadGeometry();
	AliTRDonlineTrackMatching *matching = new AliTRDonlineTrackMatching();
	
	AliESDtrackCuts trackCutsV0("AlitrackCutsV0", "AlitrackCutsV0");
	trackCutsV0.SetEtaRange(-0.9,0.9);
	trackCutsV0.SetAcceptKinkDaughters(kTRUE);
	trackCutsV0.SetRequireTPCRefit(kFALSE);
	trackCutsV0.SetMaxChi2PerClusterTPC(6);
	trackCutsV0.SetMinNClustersTPC(60);
		
		
	for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
		tCharge = 0;
		tTRDtrigHNU = 0;
		tTRDtrigHQU = 0;
		tPt = -1;
		tPx = -99;
		tPy = -99;
		tPz = -99;
		tY = -99;
		tTPCDEdx = -1;
		tTOFSignal = -1;
		tDcaXY = -999;
		tDcaZ = -999;
		tTRDPid = 0;
		tTRDnTracklets = 0;
		tTRDPt = 0;
		tTRDLayerMask = 0;
		tTRDSagitta = -1;
		
		AliESDTrdTrack* trdTrack = fESD->GetTrdTrack(iTrack);
		if (!trdTrack) continue;
		tTRDPid = trdTrack->GetPID();
		tTRDnTracklets = trdTrack->GetNTracklets();
		tTRDPt= (TMath::Abs(trdTrack->GetPt()));
		tTRDLayerMask = trdTrack->GetLayerMask();

		Float_t invPtDev = 0;
		Int_t b = trdTrack->GetB();
		Int_t c = trdTrack->GetC();
		if (b==0 && c==0) invPtDev = 0.5;
		else {
			Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
			tmp += (c & 0xfff);
			invPtDev = tmp * 0.000001;
		}
		tTRDSagitta = invPtDev;
		
		// simulate HNU
		if ((tTRDPid >= 255 && tTRDnTracklets == 4) || (tTRDPid >= 235 && tTRDnTracklets > 4)) {	
				tTRDtrigHNU = 1;
		}
		// simulate HQU
		if (tTRDPt >= 256 && tTRDPid >= 130 && tTRDnTracklets >= 5 && (tTRDLayerMask & 1) ) {	
			if (tTRDSagitta < 0.2 && tTRDSagitta > -0.2) {
				tTRDtrigHQU = 1;	
			}
		}
		
		if (!tTRDtrigHQU && !tTRDtrigHNU) continue;
		
		AliESDtrack* bestTrack = 0x0;

    Double_t esdPt = 0;
    Double_t mag = fESD->GetMagneticField();
    Double_t currentMatch = 0;
    Double_t bestMatch = 0;
    
    for(Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
	    AliESDtrack* esdTrack = static_cast<AliESDtrack*>(fESD->GetTrack(i));   
      if(!esdTrack || !trackCutsV0.AcceptTrack(esdTrack)) continue; 
          
				esdPt = esdTrack->GetSignedPt();
        Double_t gtuPt = trdTrack->Pt();

        if(mag > 0.) gtuPt = gtuPt * (-1.0);

        Double_t ydist;
        Double_t zdist;
        
        if (matching->EstimateTrackDistance(esdTrack, trdTrack, mag, &ydist, &zdist) == 0) {
        	currentMatch = matching->RateTrackMatch(ydist, zdist, esdPt, gtuPt);
				}
     
        if(currentMatch > bestMatch) {
            bestMatch = currentMatch;
            bestTrack = esdTrack;
        }
    }
    
    if (bestTrack) {
			Double_t momvect[3];
			bestTrack->PxPyPz(momvect);
			tPx = momvect[0];
			tPy = momvect[1];
			tPz = momvect[2];
		  tCharge = bestTrack->GetSign();
			tPt = bestTrack->Pt();			
			tY  = bestTrack->Y();		
			tP = bestTrack->GetInnerParam()->GetP();			
			tTPCDEdx = bestTrack->GetTPCsignal();
			
			// calculate mass from TOF Signal
			Double_t mass = 0, time = -1, beta = 0, gamma = 0, length = 0, time0 = 0;
			length = bestTrack->GetIntegratedLength();
			
			time0 = fPIDResponse->GetTOFResponse().GetStartTime(bestTrack->P());
    	time = bestTrack->GetTOFsignal() - time0;
			if (time > 0 && length > 0) {
				beta = length / (2.99792457999999984e-02 * time);
				if (1 - beta*beta > 0) {
					gamma = 1 / TMath::Sqrt(1 - beta*beta);
					mass = (tP / TMath::Sqrt(gamma*gamma - 1));
	 				tTOFSignal = mass*mass;
				}
			}
	
			Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
			bestTrack->GetImpactParameters(dca, cov);
			tDcaXY = dca[0];
			tDcaZ = dca[1];
    }
  
		fTree->Fill();	
		
	}
	
		delete matching;
		                                                  
    PostData(1, fOutputList);                          
  	PostData(2, fTree);                                                      
}
//_____________________________________________________________________________
void AliAnalysisTaskTRDtriggerTracks::Terminate(Option_t *)
{
}
//_____________________________________________________________________________
