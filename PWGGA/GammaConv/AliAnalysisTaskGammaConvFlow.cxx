/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Andrea Dubla, Redmer Alexander Bertens, Friederike Bock	      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	      *
 * provided "as is" without express or implied warranty.       		      *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------

// Class used to do analysis on conversion pairs
//---------------------------------------------
///////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaConvFlow.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"

#include "AliFlowCandidateTrack.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEventSimple.h"
#include "AliFlowCommonConstants.h"
#include "AliFlowEvent.h"
#include "AliFlowTrack.h"

class AliFlowTrackCuts;

using namespace std;

ClassImp(AliAnalysisTaskGammaConvFlow)

//________________________________________________________________________
AliAnalysisTaskGammaConvFlow::AliAnalysisTaskGammaConvFlow(): AliAnalysisTaskSE(),
fV0Reader(NULL),
fBGHandler(NULL),
fBGHandlerRP(NULL),
fInputEvent(NULL),

fCutFolder(NULL),
fESDList(NULL),
fBackList(NULL),
fMotherList(NULL),
fPhotonDCAList(NULL),
fMesonDCAList(NULL),
//fTrueList(NULL),
//fMCList(NULL),
fHeaderNameList(NULL),
fOutputContainer(0),
fReaderGammas(NULL),
fGammaCandidates(NULL),
fEventCutArray(NULL),
fEventCuts(NULL),
fCutArray(NULL),
fConversionCuts(NULL),
fMesonCutArray(NULL),
fMesonCuts(NULL),
hESDConvGammaPt(NULL),
hInvMassPair(NULL),
hKappaTPC(NULL),
hKappaTPC_after(NULL),
hESDConvGammaR(NULL),
hESDConvGammaEta(NULL),
//tESDConvGammaPtDcazCat(NULL),
fPtGamma(0),
fDCAzPhoton(0),
fRConvPhoton(0),
fEtaPhoton(0),
iCatPhoton(0),
iPhotonMCInfo(0),
hESDMotherInvMassPt(NULL),
//sESDMotherInvMassPtZM(NULL),
hESDMotherBackInvMassPt(NULL),
//sESDMotherBackInvMassPtZM(NULL),
hESDMotherInvMassEalpha(NULL),
hESDMotherPi0PtY(NULL),
hESDMotherEtaPtY(NULL),
hESDMotherPi0PtAlpha(NULL),
hESDMotherEtaPtAlpha(NULL),
hESDMotherPi0PtOpenAngle(NULL),
hESDMotherEtaPtOpenAngle(NULL),

hNEvents(NULL),
hNGoodESDTracks(NULL),
hNGammaCandidates(NULL),
hNGoodESDTracksVsNGammaCanditates(NULL),
hNV0Tracks(NULL),
hEtaShift(NULL),
tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
fInvMass(0),
fPt(0),
fDCAzGammaMin(0),
fDCAzGammaMax(0),
iFlag(0),
iMesonMCInfo(0),
fEventPlaneAngle(-100),
fRandom(0),
fnGammaCandidates(0),
fUnsmearedPx(NULL),
fUnsmearedPy(NULL),
fUnsmearedPz(NULL),
fUnsmearedE(NULL),
fMCStackPos(NULL),
fMCStackNeg(NULL),
fESDArrayPos(NULL),
fESDArrayNeg(NULL),
fnCuts(0),
fiCut(0),
fMoveParticleAccordingToVertex(kTRUE),
fIsHeavyIon(0),
fDoMesonAnalysis(kTRUE),
fDoMesonQA(0),
fDoPhotonQA(0),
fIsFromMBHeader(kTRUE),
fhistoEPVZ(NULL),
fMinMass(-1),
fMaxMass(10),
fMinKappa(-1),
fMaxKappa(100),
fFilterVariable(1),
fMinFilter(0),
fMaxFilter(0.2),
fDebug(0),
fCutsRP(0),
fNullCuts(0), 
fFlowEvent(0),
fIsMC(0),
fMCEvent(NULL),
fMCStack(NULL)
{
	// DefineOutput(1, TList::Class());
	// DefineOutput(2, AliFlowEventSimple::Class());
}

//________________________________________________________________________
AliAnalysisTaskGammaConvFlow::AliAnalysisTaskGammaConvFlow(const char *name):
AliAnalysisTaskSE(name),
fV0Reader(NULL),
fBGHandler(NULL),
fBGHandlerRP(NULL),
fInputEvent(NULL),
fCutFolder(NULL),
fESDList(NULL),
fBackList(NULL),
fMotherList(NULL),
fPhotonDCAList(NULL),
fMesonDCAList(NULL),
//fTrueList(NULL),
//fMCList(NULL),
fHeaderNameList(NULL),
fOutputContainer(0),
fReaderGammas(NULL),
fGammaCandidates(NULL),
fEventCutArray(NULL),
fEventCuts(NULL),
fCutArray(NULL),
fConversionCuts(NULL),
fMesonCutArray(NULL),
fMesonCuts(NULL),
hESDConvGammaPt(NULL),
hInvMassPair(NULL),
hKappaTPC(NULL),
hKappaTPC_after(NULL),
hESDConvGammaR(NULL),
hESDConvGammaEta(NULL),
//tESDConvGammaPtDcazCat(NULL),
fPtGamma(0),
fDCAzPhoton(0),
fRConvPhoton(0),
fEtaPhoton(0),
iCatPhoton(0),
iPhotonMCInfo(0),
hESDMotherInvMassPt(NULL),
//sESDMotherInvMassPtZM(NULL),
hESDMotherBackInvMassPt(NULL),
//sESDMotherBackInvMassPtZM(NULL),
hESDMotherInvMassEalpha(NULL),
hESDMotherPi0PtY(NULL),
hESDMotherEtaPtY(NULL),
hESDMotherPi0PtAlpha(NULL),
hESDMotherEtaPtAlpha(NULL),
hESDMotherPi0PtOpenAngle(NULL),
hESDMotherEtaPtOpenAngle(NULL),
hNEvents(NULL),
hNGoodESDTracks(NULL),
hNGammaCandidates(NULL),
hNGoodESDTracksVsNGammaCanditates(NULL),
hNV0Tracks(NULL),
hEtaShift(NULL),
tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
fInvMass(0),
fPt(0),
fDCAzGammaMin(0),
fDCAzGammaMax(0),
iFlag(0),
iMesonMCInfo(0),
fEventPlaneAngle(-100),
fRandom(0),
fnGammaCandidates(0),
fUnsmearedPx(NULL),
fUnsmearedPy(NULL),
fUnsmearedPz(NULL),
fUnsmearedE(NULL),
fMCStackPos(NULL),
fMCStackNeg(NULL),
fESDArrayPos(NULL),
fESDArrayNeg(NULL),
fnCuts(0),
fiCut(0),
fMoveParticleAccordingToVertex(kTRUE),
fIsHeavyIon(0),
fDoMesonAnalysis(kTRUE),
fDoMesonQA(0),
fDoPhotonQA(0),
fIsFromMBHeader(kTRUE),
fhistoEPVZ(NULL),
fMinMass(-1),
fMaxMass(10),
fMinKappa(-1),
fMaxKappa(100),
fFilterVariable(1),
fMinFilter(0),
fMaxFilter(0.2),
fDebug(0),
fCutsRP(0), 
fNullCuts(0), 
fFlowEvent(0),
fIsMC(0),
fMCEvent(NULL),
fMCStack(NULL)

{
    // Define output slots here
    DefineOutput(1, TList::Class());
    DefineOutput(2, AliFlowEventSimple::Class());

}
//___________________________________________________________
AliAnalysisTaskGammaConvFlow::~AliAnalysisTaskGammaConvFlow()
{
	if(fGammaCandidates){
		delete fGammaCandidates;
		fGammaCandidates = 0x0;
	}
	if(fBGHandler){
		delete[] fBGHandler;
		fBGHandler = 0x0;
	}
	if(fBGHandlerRP){
		delete[] fBGHandlerRP;
		fBGHandlerRP = 0x0;
	}
	
	if (fFlowEvent) delete fFlowEvent;

}

//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UserCreateOutputObjects(){
	
	// Create histograms
	if(fOutputContainer != NULL){
		delete fOutputContainer;
		fOutputContainer = NULL;
	}
	if(fOutputContainer == NULL){
		fOutputContainer = new TList();
		fOutputContainer->SetOwner(kTRUE);
	}
	
	//========================= again flow setting==========================
	//Create histograms
	//----------hfe initialising begin---------
	fNullCuts = new AliFlowTrackCuts("null_cuts");
	
	AliFlowCommonConstants* cc = AliFlowCommonConstants::GetMaster();
	cc->SetNbinsMult(10000);
	cc->SetMultMin(0);
	cc->SetMultMax(10000);
	
	cc->SetNbinsPt(200);
	cc->SetPtMin(0);
	cc->SetPtMax(20);
	
	cc->SetNbinsPhi(180);
	cc->SetPhiMin(0.0);
	cc->SetPhiMax(TMath::TwoPi());
	
	cc->SetNbinsEta(40);
	cc->SetEtaMin(-2.0);
	cc->SetEtaMax(+2.0);
	
	cc->SetNbinsQ(500);
	cc->SetQMin(0.0);
	cc->SetQMax(3.0);
  
  cc->SetMassMin(fMinFilter);
  cc->SetMassMax(fMaxFilter);
  cc->SetNbinsMass(100);

	// Array of current cut's gammas
	fGammaCandidates = new TList();
	
	fCutFolder = new TList*[fnCuts];
	fESDList = new TList*[fnCuts];
	fBackList = new TList*[fnCuts];
	fMotherList = new TList*[fnCuts];
	hNEvents = new TH1I*[fnCuts];
	hNGoodESDTracks = new TH1I*[fnCuts];
	hNGammaCandidates = new TH1I*[fnCuts];
	hNGoodESDTracksVsNGammaCanditates = new TH2F*[fnCuts];
	hNV0Tracks = new TH1I*[fnCuts];
	hEtaShift = new TProfile*[fnCuts];
	hESDConvGammaPt = new TH1F*[fnCuts];
	hInvMassPair = new TH2F*[fnCuts];
	hKappaTPC = new TH2F*[fnCuts];
  hKappaTPC_after = new TH2F*[fnCuts];
	
	if (fDoPhotonQA == 2){
		fPhotonDCAList = new TList*[fnCuts];
//      tESDConvGammaPtDcazCat = new TTree*[fnCuts];
	}
	if (fDoPhotonQA > 0){
		hESDConvGammaR = new TH1F*[fnCuts];
		hESDConvGammaEta = new TH1F*[fnCuts];
	}
	
	if(fDoMesonAnalysis){
		hESDMotherInvMassPt = new TH2F*[fnCuts];
		hESDMotherBackInvMassPt = new TH2F*[fnCuts];
		hESDMotherInvMassEalpha = new TH2F*[fnCuts];
		if (fDoMesonQA == 2){
			fMesonDCAList = new TList*[fnCuts];
			tESDMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
		}
		if (fDoMesonQA > 0){
			hESDMotherPi0PtY =  new TH2F*[fnCuts];
			hESDMotherEtaPtY =  new TH2F*[fnCuts];
			hESDMotherPi0PtAlpha =  new TH2F*[fnCuts];
			hESDMotherEtaPtAlpha =  new TH2F*[fnCuts];
			hESDMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
			hESDMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
		}
	}
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		
		TString cutstringEvent 	= ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPhoton = ((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson = "NoMesonCut";
		if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
		
		fCutFolder[iCut] = new TList();
		fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
		fCutFolder[iCut]->SetOwner(kTRUE);
		fOutputContainer->Add(fCutFolder[iCut]);
		fESDList[iCut] = new TList();
		fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
		fESDList[iCut]->SetOwner(kTRUE);
		fCutFolder[iCut]->Add(fESDList[iCut]);
		
		hNEvents[iCut] = new TH1I("NEvents","NEvents",12,-0.5,11.5);
		hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
		if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
			TString TriggerNames = "Not Trigger: ";
			TriggerNames = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
			hNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
		} else {
			hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
		}
		hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
		hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
		fESDList[iCut]->Add(hNEvents[iCut]);
		
		if(fIsHeavyIon == 1) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
		else if(fIsHeavyIon == 2) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
		else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
		fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
		if(fIsHeavyIon == 1) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
		else if(fIsHeavyIon == 2) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		else hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
		fESDList[iCut]->Add(hNGammaCandidates[iCut]);
		if(fIsHeavyIon == 1) hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",4000,0,4000,100,0,100);
		else if(fIsHeavyIon == 2) hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",400,0,400,50,0,50);
		else hNGoodESDTracksVsNGammaCanditates[iCut] = new TH2F("GoodESDTracksVsGammaCandidates","GoodESDTracksVsGammaCandidates",200,0,200,50,0,50);
		fESDList[iCut]->Add(hNGoodESDTracksVsNGammaCanditates[iCut]);
		
		
		if(fIsHeavyIon == 1) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
		else if(fIsHeavyIon == 2) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
		else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
		fESDList[iCut]->Add(hNV0Tracks[iCut]);
		hEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
		fESDList[iCut]->Add(hEtaShift[iCut]);
		hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
		fESDList[iCut]->Add(hESDConvGammaPt[iCut]);
		
		hInvMassPair[iCut]= new TH2F("InvMassPair_Pt","Gamma invariant mass vs Pt",200,0,0.2,250,0,25);	
		fESDList[iCut]->Add(hInvMassPair[iCut]);
		
		hKappaTPC[iCut]= new TH2F("KappaTPC_Pt","Gamma KappaTPC vs Pt",200,0,10,250,0,25);
		fESDList[iCut]->Add(hKappaTPC[iCut]);
    
    hKappaTPC_after[iCut]= new TH2F("KappaTPC_Pt_after","Gamma KappaTPC vs Pt after cuts",200,0,10,250,0,25);
    fESDList[iCut]->Add(hKappaTPC_after[iCut]);
    
    //2d histogram filling the cut and value - control check for selections
		
		if (fDoPhotonQA == 2){
			fPhotonDCAList[iCut] = new TList();
			fPhotonDCAList[iCut]->SetName(Form("%s_%s_%s Photon DCA tree",cutstringEvent.Data(),cutstringPhoton.Data(),cutstringMeson.Data()));
			fPhotonDCAList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fPhotonDCAList[iCut]);
			
//            tESDConvGammaPtDcazCat[iCut] = new TTree("ESD_ConvGamma_Pt_Dcaz_R_Eta","ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");
//            tESDConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/F");
//            tESDConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/F");
			//          tESDConvGammaPtDcazCat[iCut]->Branch("R",&fRConvPhoton,"fRConvPhoton/F");
			//          tESDConvGammaPtDcazCat[iCut]->Branch("Eta",&fEtaPhoton,"fEtaPhoton/F");
			
		//    tESDConvGammaPtDcazCat[iCut]->Branch("cat",&iCatPhoton,"iCatPhoton/b");

	//        fPhotonDCAList[iCut]->Add(tESDConvGammaPtDcazCat[iCut]);
		}
		
		if (fDoPhotonQA > 0){
			hESDConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
			fESDList[iCut]->Add(hESDConvGammaR[iCut]);
			hESDConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
			fESDList[iCut]->Add(hESDConvGammaEta[iCut]);
		}
		
		if(fDoMesonAnalysis){
			hESDMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
			hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);
			hESDMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
			fESDList[iCut]->Add(hESDMotherInvMassEalpha[iCut]);
			if (fDoMesonQA == 2){
				fMesonDCAList[iCut] = new TList();
				fMesonDCAList[iCut]->SetName(Form("%s_%s_%s Meson DCA tree",cutstringEvent.Data() ,cutstringPhoton.Data(),cutstringMeson.Data()));
				fMesonDCAList[iCut]->SetOwner(kTRUE);
				fCutFolder[iCut]->Add(fMesonDCAList[iCut]);
				
				tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree("ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag","ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");
				tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
				tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
				tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
				tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
				
				fMesonDCAList[iCut]->Add(tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]);
				
			}
			if (fDoMesonQA > 0 ){
				hESDMotherPi0PtY[iCut] = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
				SetLogBinningXTH2(hESDMotherPi0PtY[iCut]);
				fESDList[iCut]->Add(hESDMotherPi0PtY[iCut]);
				hESDMotherEtaPtY[iCut] = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
				SetLogBinningXTH2(hESDMotherEtaPtY[iCut]);
				fESDList[iCut]->Add(hESDMotherEtaPtY[iCut]);
				hESDMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",150,0.03,15.,100,0,1);
				SetLogBinningXTH2(hESDMotherPi0PtAlpha[iCut]);
				fESDList[iCut]->Add(hESDMotherPi0PtAlpha[iCut]);
				hESDMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",150,0.03,15.,100,0,1);
				SetLogBinningXTH2(hESDMotherEtaPtAlpha[iCut]);
				fESDList[iCut]->Add(hESDMotherEtaPtAlpha[iCut]);
				hESDMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",150,0.03,15.,100,0,TMath::Pi());
				SetLogBinningXTH2(hESDMotherPi0PtOpenAngle[iCut]);
				fESDList[iCut]->Add(hESDMotherPi0PtOpenAngle[iCut]);
				hESDMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",150,0.03,15.,100,0,TMath::Pi());
				SetLogBinningXTH2(hESDMotherEtaPtOpenAngle[iCut]);
				fESDList[iCut]->Add(hESDMotherEtaPtOpenAngle[iCut]);
			}
			
			
		}
		
		
	}
//     if(fDoMesonAnalysis){
//         InitBack(); // Init Background Handler
//     }
	

	
	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
	
	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
		if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
		}
		if(fDoMesonAnalysis){
			if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
				fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
			}
		}
	}
	
	fhistoEPVZ = new TH1D("EPVZ", "EPVZ", 60, -TMath::Pi()/2, TMath::Pi()/2);
    fOutputContainer->Add(fhistoEPVZ);


	
	PostData(1, fOutputContainer);

	fFlowEvent = new AliFlowEvent(10000);
	PostData(2, fFlowEvent);
	
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::Notify()
{
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))->GetDoEtaShift()){
			hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			continue; // No Eta Shift requested, continue
		}
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
			((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
			hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
			continue;
		}
		else{
			printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
				(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber()).Data(),((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift());
			hEtaShift[iCut]->Fill(0.,(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetEtaShift()));
			((AliConvEventCuts*)fEventCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
		}
	}
	return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UserExec(Option_t *)
{
	//
	// Called for each event
	//
	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
	if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
		for(Int_t iCut = 0; iCut<fnCuts; iCut++){
			hNEvents[iCut]->Fill(eventQuality);
		}
		return;
	}
	
	fInputEvent = InputEvent();
  if(fIsMC) fMCEvent = MCEvent();
  if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){ fMCStack = fMCEvent->Stack(); }

	fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
	
	// ------------------- BeginEvent ----------------------------
	
	AliEventplane *EventPlane = fInputEvent->GetEventplane();
	if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
	else fEventPlaneAngle=0.0;
	
	SetNullCuts(fInputEvent);
	PrepareFlowEvent(fInputEvent->GetNumberOfTracks(),fFlowEvent);    //Calculate event plane Qvector and EP resolution for inclusive

	

	
// 	for(Int_t iCut = 0; iCut<fnCuts; iCut++){
		Int_t iCut = 0;
		fiCut = iCut;
		Int_t eventNotAccepted =
		((AliConvEventCuts*)fEventCutArray->At(iCut))
		->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
		if(eventNotAccepted){
			// cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
			hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
			return;
		}
		
		if(eventQuality != 0){// Event Not Accepted
			// cout << "event rejected due to: " <<eventQuality << endl;
			hNEvents[iCut]->Fill(eventQuality);
			return;
		}
				
		hNEvents[iCut]->Fill(eventQuality); // Should be 0 here
		hNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks());
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->IsHeavyIon() == 2)	hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
		else hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
		
		ProcessPhotonCandidates(); // Process this cuts gammas
		ProcessPhotonCandidatesforV2(); // Process this cuts gammas and do v2

		hNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries());
		hNGoodESDTracksVsNGammaCanditates[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks(),fGammaCandidates->GetEntries());

		fGammaCandidates->Clear(); // delete this cuts good gammas
// 	}
	
	fhistoEPVZ->Fill(fEventPlaneAngle);

	PostData(1, fOutputContainer);
	PostData(2, fFlowEvent);

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::ProcessPhotonCandidates()
{
    Int_t nV0 = 0;
    TList *GammaCandidatesStepOne = new TList();
    TList *GammaCandidatesStepTwo = new TList();
    // Loop over Photon Candidates allocated by ReaderV1
    for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
        AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
        if(!PhotonCandidate) continue;
        fIsFromMBHeader = kTRUE;
        
        hKappaTPC[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
        if (PhotonCandidate->GetInvMassPair() < fMinMass || PhotonCandidate->GetInvMassPair() > fMaxMass) continue;
        if (((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent) < fMinKappa || ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent) > fMaxKappa) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
        if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
           !((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
            fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas
            
            if(fIsFromMBHeader){
                hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                if (fDoPhotonQA > 0){
                    hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                    hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                }
            }

            
            if (fIsFromMBHeader && fDoPhotonQA == 2){
                if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                   // tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                  //  tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        } else if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
            ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
            nV0++;
            GammaCandidatesStepOne->Add(PhotonCandidate);
        } else if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
                  ((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
            GammaCandidatesStepTwo->Add(PhotonCandidate);
        }
    }
    
    
    
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
        for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
            AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
            if(!PhotonCandidate) continue;
            fIsFromMBHeader = kTRUE;

            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
                fGammaCandidates->Add(PhotonCandidate);
                if(fIsFromMBHeader){
                    hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                    hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                    hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                    if (fDoPhotonQA > 0){
                        hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                        hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                    }
                }
            }

                
                GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
            
            if (fIsFromMBHeader && fDoPhotonQA == 2){
                if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                //    tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                //    tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        }
    }
    if(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
        for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
            AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
            if(!PhotonCandidate) continue;
            fIsFromMBHeader = kTRUE;

            if(!((AliConversionPhotonCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
            fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
            if(fIsFromMBHeader){
                hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
                hInvMassPair[fiCut]->Fill(PhotonCandidate->GetInvMassPair(),PhotonCandidate->Pt());
                hKappaTPC_after[fiCut]->Fill(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(PhotonCandidate,fInputEvent),PhotonCandidate->Pt());
                if (fDoPhotonQA > 0){
                    hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                    hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
                }
            }

            if (fIsFromMBHeader){
                if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
               //     tESDConvGammaPtDcazCat[fiCut]->Fill();
                } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
                    fPtGamma = PhotonCandidate->Pt();
                    fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
                    fRConvPhoton = PhotonCandidate->GetConversionRadius();
                    fEtaPhoton = PhotonCandidate->GetPhotonEta();
                    iCatPhoton = PhotonCandidate->GetPhotonQuality();
                  //  tESDConvGammaPtDcazCat[fiCut]->Fill();
                }
            }
        }
    }
    
    delete GammaCandidatesStepOne;
    GammaCandidatesStepOne = 0x0;
    delete GammaCandidatesStepTwo;
    GammaCandidatesStepTwo = 0x0;
    
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::UpdateEventByEventData(){
    //see header file for documentation
    if(fGammaCandidates->GetEntries() >0 ){
        if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
            fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
        }
        else{ // means we use #V0s for multiplicity
            fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
        }
    }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::SetLogBinningXTH2(TH2* histoRebin){
    TAxis *axisafter = histoRebin->GetXaxis(); 
    Int_t bins = axisafter->GetNbins();
    Double_t from = axisafter->GetXmin();
    Double_t to = axisafter->GetXmax();
    Double_t *newbins = new Double_t[bins+1];
    newbins[0] = from;
    Double_t factor = TMath::Power(to/from, 1./bins);
    for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
    axisafter->Set(bins, newbins);
    delete [] newbins;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::Terminate(const Option_t *)
{
    
    //fOutputContainer->Print(); // Will crash on GRID
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::ProcessPhotonCandidatesforV2()
{

	// Loop over Photon Candidates allocated by ReaderV1	
	for(Int_t i = 0; i < fGammaCandidates->GetEntries(); i++){
		
		AliAODConversionPhoton *gammaForv2=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(i));
		if (gammaForv2 == NULL) return;
    AliFlowTrack *sTrack = new AliFlowTrack();
		sTrack->SetForRPSelection(kFALSE);
		sTrack->SetForPOISelection(kTRUE);
    
    if(fFilterVariable==1){//using mass for POI selections
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==2){//using Kappa for POI selections
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else if(fFilterVariable==3 && fIsMC){//MC ElectronElectron + mass
      if(!MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==4 && fIsMC){//MC ElectronElectron + kappa
      if(!MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else if(fFilterVariable==5 && fIsMC){//MC Not ElectronElectron + mass
      if(MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(gammaForv2->GetInvMassPair());
    }else if(fFilterVariable==6 && fIsMC){//MC Not ElectronElectron + kappa
      if(MCElectronElectron(gammaForv2)) return;
      sTrack->SetMass(((AliConversionPhotonCuts*)fCutArray->At(fiCut))->GetKappaTPC(gammaForv2,fInputEvent));
    }else{//no additional POI selection
      sTrack->SetMass(123456);
    }
		sTrack->SetPt(gammaForv2->Pt());
		sTrack->SetPhi(gammaForv2->GetPhotonPhi());
		sTrack->SetEta(gammaForv2->GetPhotonEta());

	/*     for(int iRPs=0; iRPs!=fFlowEvent->NumberOfTracks(); ++iRPs)
		{
			//   cout << " no of rps " << iRPs << endl;
			AliFlowTrack *iRP = dynamic_cast<AliFlowTrack*>(fFlowEvent->GetTrack( iRPs ));
			if (!iRP) continue;
			if (!iRP->InRPSelection()) continue;
			if( sTrack->GetID() == iRP->GetID())
			{
				if(fDebug) printf(" was in RP set");
				//  cout << sTrack->GetID() <<"   ==  " << iRP->GetID() << " was in RP set====REMOVED" <<endl;
				iRP->SetForRPSelection(kFALSE);
				// fFlowEvent->SetNumberOfRPs(fFlowEvent->GetNumberOfRPs() - 1);
			}
		} //end of for loop on RPs*/
		fFlowEvent->InsertTrack(((AliFlowTrack*) sTrack));
		fFlowEvent->SetNumberOfPOIs(fFlowEvent->GetNumberOfPOIs()+1);
	}
}

//_____________________________________________________________________________
template <typename T> void AliAnalysisTaskGammaConvFlow::SetNullCuts(T* event)
{
	// Set null cuts
	if (fDebug) cout << " fCutsRP " << fCutsRP << endl;
	fCutsRP->SetEvent(event, MCEvent());
	fNullCuts->SetParamType(AliFlowTrackCuts::kGlobal);
	fNullCuts->SetPtRange(+1, -1); // select nothing QUICK
	fNullCuts->SetEtaRange(+1, -1); // select nothing VZERO
	fNullCuts->SetEvent(event, MCEvent());
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaConvFlow::PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const
{
   //Prepare flow events
    FlowEv->ClearFast();
    FlowEv->Fill(fCutsRP, fNullCuts);
    FlowEv->SetReferenceMultiplicity(iMulti);
    FlowEv->DefineDeadZone(0, 0, 0, 0);
    //  FlowEv->TagSubeventsInEta(-0.7, 0, 0, 0.7);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvFlow::MCElectronElectron( AliAODConversionPhoton *MCPhoton ){
  
  
  TParticle *posDaughter = MCPhoton->GetPositiveMCDaughter(fMCStack);
  TParticle *negDaughter = MCPhoton->GetNegativeMCDaughter(fMCStack);
  if(posDaughter==NULL || negDaughter==NULL) return kFALSE;
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0;
  Bool_t IsItElectronElectron;
  
  if( posDaughter->GetMother(0) != negDaughter->GetMother(0)  || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) {
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) IsItElectronElectron = kTRUE;
  }else{
    IsItElectronElectron = kFALSE;
  }
  
  return IsItElectronElectron;
    
}






