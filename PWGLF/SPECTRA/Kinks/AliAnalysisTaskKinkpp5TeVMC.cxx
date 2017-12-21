//==========================================================//
// Date: 11/5/2017
//Authors: Nur Hussain and Buddhadeb Bhattacharjee, Department of Physics, Gauhati University
//Special thanks to  Martha Spyropoulou-Stassinaki
//purpose:: Kink topology for pp 5TeV for MC data
//=========================================================//

#include "Riostream.h" 
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "TParticle.h"
#include <TVector3.h>
#include "TF1.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
//#include "AliPhysicsSelectionTask.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
//#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTaskKinkpp5TeVMC.h"
#include "AliCentrality.h"
#include "iostream"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"



ClassImp(AliAnalysisTaskKinkpp5TeVMC)

//_________________________________________________________________
AliAnalysisTaskKinkpp5TeVMC::AliAnalysisTaskKinkpp5TeVMC()
: AliAnalysisTaskSE(),fOutputList(0), fHistPt(0),fVtxCut(10.),fMultiplicity(0),fIncompletEvent(0),fMultpileup(0), fMultV0trigger(0),fZvertex(0),fEventVertex(0),
        fRatioCrossedRows(0),fZvXv(0), fZvYv(0), fXvYv(0),fRpr(0),fdcaToVertexXY(0),fdcaToVertexXYafterCut(0),fptAllKink(0),fRatioCrossedRowsKink(0),fPosiKink(0),
        fQtAll(0),fptKink(0),fQtMothP(0),fqT1(0),fEta(0),fqT2(0),fKinkKaonBackg(0),f1(0), f2(0),fPtCut1(0),fAngMotherPi(0),
        fQtInvM(0),fInvMuNuAll(0),fInvMassMuNuPtAll(0),fRadiusPt(0),fKinkRadUp(200.), fKinkRadLow(130.), fLowCluster(30), fLowQt(.12), fRapiK(0.5),
        fAngMotherKC(0),fkaonInvaiant(0),fRadiusNcl(0), fPtKPDG(0),fAngMotherKKinks(0),fPtCut2(0),fPtCut3(0),fTPCSignlMotherK(0),fPtKaon(0), fPtKaonP(0), fPtKaonN(0),
         fTPCSignalP(0),fRadiusNclCln(0),fRadiusPtcln(0),fInvMassMuNuPt(0),fTPCSignlPtpc(0),fMothKinkMomSignl(0),fTPCSignlKinkDau(0),fTPCMomNSigmaAllKaon(0),
        fnSigmaTPC(0),fradiurKink(0),fLenthKink(0),fEtaK(0),frapiKESD(0),fzVertexPositionKinKvsKinkRad(0),fSignPtNcl(0),fSignPtrapiK(0),frapiKNcl(0),fSignPt(0),
        fChi2NclTPC(0),fRatioChi2Ncl(0),flifetime(0),fPtKinkKaon(0),fDCAkink(0),fPtKink(0),fPtKinkPos(0),fPtKinkNeg(0),fPtKinkK0(0),fPtKinkK0P(0),fPtKinkK0N(0),
        fPtKinkGyu(0),fPtKinkGyuP(0),fPtKinkGyuN(0),fKinKRbn(0),fradPtRpDt(0),fAngMomK(0),fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),fPIDResponse(0),fNumberOfEvent(0),

        fMultMC_wo_any_cut(0),fMultMC_incompleteDAQ(0),fMultMC_AfterPileUp(0),fMultTriggerMCAfterV0(0),fTrigSel(AliVEvent::kINT7),
        frapidKMC(0), fptKMC(0),fSignPtGen(0),fPtKPlMC(0),fPtKMnMC(0),flengthTrackRef(0),flifetimeMC(0),flifeSmallMC(0),flifeInt(0),flifeYuri(0),flifetiMCK(0),
        flenYuri(0),flengthMCK(0),flifetime_MCprocess4(0),fradPtRapMC(0), flifetime_kaonpionPDG(0),flifetime_kaonmuonPDG(0),fmaxAngMomKmuMC(0),fradPtRapDC(0),
        fradMC(0),fQtKMuMC(0),fgenPtEtR(0),fgenPtEtRP(0),fMCEtaKaon(0),fSignPtEtaMC(0),fSignPtMC(0),fgenPtEtRN(0),fQtKElMC(0),fQtKPiMC(0),fESDMult(0),
        fFakepipi(0),fFakeKPi(0),fRadiusPtPion(0),fRadiusPtKaon(0),fQtKMu(0),fQtKPi(0),fQtKEl(0),fQtK3PiP(0),fQtK3PiM(0),fHistPtKPDG(0), fHiPtKPDGP(0), fHiPtKPDGN(0),
        fHistEta(0), frapidESDK(0), fHistQt2(0),fAngMomPi(0),fMinvPi(0),fMinvKa(0),fcodeH(0), fZkinkZDau(0),fRadiusPtFake(0),fTPCMomNSgnl(0),fPtPrKink(0),flifTiESDK(0),
        fKinkKaon(0),fkinkKaonP(0),fkinkKaonN(0),fcode2(0),fTPCSgnlPtpc(0),fMothKinkMomSgnlD(0),fKinkKaonBg(0),fMothKinkMomSgnl(0),fcodeDau2(0),fTPCSgnlKinkDau(0),
        fMinvPr(0),fDCAkinkBG(0),fPosiKinKBgXY(0),fPosiKinKBgZY(0),fPosiKinKBgZX(0),fKinKBGP(0),fKinKBGN(0),fdcodeH(0),fcode4(0),
        fNumberOfEvent_cent(0),fESDtrackCuts(0),fEventVsCentrality(0)
{}

//________________________________________________________________________
AliAnalysisTaskKinkpp5TeVMC::AliAnalysisTaskKinkpp5TeVMC(const char *name, Float_t lRadiusKUp,  Float_t lRadiusKLow, Int_t lNCluster, Float_t lLowQtValue, Float_t yRange) 
  : AliAnalysisTaskSE(name),  fOutputList(0), fHistPt(0),fVtxCut(10.),fMultiplicity(0),fIncompletEvent(0),fMultpileup(0), fMultV0trigger(0),fZvertex(0),fEventVertex(0),
	fRatioCrossedRows(0),fZvXv(0), fZvYv(0), fXvYv(0),fRpr(0),fdcaToVertexXY(0),fdcaToVertexXYafterCut(0),fptAllKink(0),fRatioCrossedRowsKink(0),fPosiKink(0),
	fQtAll(0),fptKink(0),fQtMothP(0),fqT1(0),fEta(0),fqT2(0),fKinkKaonBackg(0),f1(0), f2(0),fPtCut1(0),fAngMotherPi(0),
	fQtInvM(0),fInvMuNuAll(0),fInvMassMuNuPtAll(0),fRadiusPt(0),fKinkRadUp(200.), fKinkRadLow(130.), fLowCluster(30), fLowQt(.12), fRapiK(0.5),
	fAngMotherKC(0),fkaonInvaiant(0),fRadiusNcl(0), fPtKPDG(0),fAngMotherKKinks(0),fPtCut2(0),fPtCut3(0),fTPCSignlMotherK(0),fPtKaon(0), fPtKaonP(0), fPtKaonN(0),
	 fTPCSignalP(0),fRadiusNclCln(0),fRadiusPtcln(0),fInvMassMuNuPt(0),fTPCSignlPtpc(0),fMothKinkMomSignl(0),fTPCSignlKinkDau(0),fTPCMomNSigmaAllKaon(0),
	fnSigmaTPC(0),fradiurKink(0),fLenthKink(0),fEtaK(0),frapiKESD(0),fzVertexPositionKinKvsKinkRad(0),fSignPtNcl(0),fSignPtrapiK(0),frapiKNcl(0),fSignPt(0),
	fChi2NclTPC(0),fRatioChi2Ncl(0),flifetime(0),fPtKinkKaon(0),fDCAkink(0),fPtKink(0),fPtKinkPos(0),fPtKinkNeg(0),fPtKinkK0(0),fPtKinkK0P(0),fPtKinkK0N(0),
	fPtKinkGyu(0),fPtKinkGyuP(0),fPtKinkGyuN(0),fKinKRbn(0),fradPtRpDt(0),fAngMomK(0),fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),fPIDResponse(0),fNumberOfEvent(0),
	
	fMultMC_wo_any_cut(0),fMultMC_incompleteDAQ(0),fMultMC_AfterPileUp(0),fMultTriggerMCAfterV0(0),fTrigSel(AliVEvent::kINT7),
	frapidKMC(0), fptKMC(0),fSignPtGen(0),fPtKPlMC(0),fPtKMnMC(0),flengthTrackRef(0),flifetimeMC(0),flifeSmallMC(0),flifeInt(0),flifeYuri(0),flifetiMCK(0),
	flenYuri(0),flengthMCK(0),flifetime_MCprocess4(0),fradPtRapMC(0), flifetime_kaonpionPDG(0),flifetime_kaonmuonPDG(0),fmaxAngMomKmuMC(0),fradPtRapDC(0),
	fradMC(0),fQtKMuMC(0),fgenPtEtR(0),fgenPtEtRP(0),fMCEtaKaon(0),fSignPtEtaMC(0),fSignPtMC(0),fgenPtEtRN(0),fQtKElMC(0),fQtKPiMC(0),fESDMult(0),
	fFakepipi(0),fFakeKPi(0),fRadiusPtPion(0),fRadiusPtKaon(0),fQtKMu(0),fQtKPi(0),fQtKEl(0),fQtK3PiP(0),fQtK3PiM(0),fHistPtKPDG(0), fHiPtKPDGP(0), fHiPtKPDGN(0), 
	fHistEta(0), frapidESDK(0), fHistQt2(0),fAngMomPi(0),fMinvPi(0),fMinvKa(0),fcodeH(0), fZkinkZDau(0),fRadiusPtFake(0),fTPCMomNSgnl(0),fPtPrKink(0),flifTiESDK(0),
	fKinkKaon(0),fkinkKaonP(0),fkinkKaonN(0),fcode2(0),fTPCSgnlPtpc(0),fMothKinkMomSgnlD(0),fKinkKaonBg(0),fMothKinkMomSgnl(0),fcodeDau2(0),fTPCSgnlKinkDau(0),
	fMinvPr(0),fDCAkinkBG(0),fPosiKinKBgXY(0),fPosiKinKBgZY(0),fPosiKinKBgZX(0),fKinKBGP(0),fKinKBGN(0),fdcodeH(0),fcode4(0),
	fNumberOfEvent_cent(0),fESDtrackCuts(0),fEventVsCentrality(0)
	


{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container


	fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  DefineOutput(1, TList::Class());
  //llll
//////////////////  DefineOutput(2, TFile::Class());
}

//________________________________________________________________________
void AliAnalysisTaskKinkpp5TeVMC::UserCreateOutputObjects()
{
  // Create histograms
 
 // Called once

	f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
   	f1->SetParameter(0,0.493677);
   	f1->SetParameter(1,0.9127037);
   	f1->SetParameter(2,TMath::Pi());


   	f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
   	f2->SetParameter(0,0.13957018);
   	f2->SetParameter(1,0.2731374);
   	f2->SetParameter(2,TMath::Pi());



	AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
        AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
        fPIDResponse=inputHandler->GetPIDResponse();
        if (!fPIDResponse) AliFatal("This Task needs the PID response attached to the inputHandler");



  	fOutputList = new TList();
	fOutputList->SetOwner();
        fOutputList->SetName("hihello");
	const Int_t rap_bin1=20;
	const Float_t rap_low1= -0.5;
	const Float_t rap_high1=0.5;
	
	const Int_t kPtBins = 35;
  	const Int_t kMultBins = 11;
	Int_t kDcaBinsTemp = 76;


	const Int_t bin_y1 =10;
        const Int_t bin_y2 =15;
        const Int_t bin_y3 =7;
        const Float_t low_y = -0.5;
        const Float_t high_y = 0.5;
	const Float_t kDcaBinsTPConlyFactor = 5; //need to change binning of DCA plot for tpconly
	// sort pT-bins ..




	 Double_t binsPt2[77] = {0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};
	 Double_t binsPt[77] = { 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0};

  	Double_t binsDca[77] =  {-3,-2.85,-2.7,-2.55,-2.4,-2.25,-2.1,-1.95,-1.8,-1.65,-1.5,-1.35,-1.2,-1.05,-0.9,-0.75,-0.6,-0.45,-0.3,-0.285,-0.27,-0.255,-0.24,-0.225,-0.21,-0.195,-0.18,-0.165,-0.15,-0.135,-0.12,-0.105,-0.09,-0.075,-0.06,-0.045,-0.03,-0.015,0,0.015,0.03,0.045,0.06,0.075,0.09,0.105,0.12,0.135,0.15,0.165,0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3};

	Double_t gPt7Comb[48] = { //pp 7Tev bining
0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0, 1.1, 1.2,
1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.5, 6
 };
   	Double_t gPtpp[57] = { //pp 5Tev bining  //specially for my analysis
0.01,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65, 0.7, 0.75,
0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4,
2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,5.5, 6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0
 };


	// K0 binning from Peter Kalinak,
  	Double_t gPt13K0PKal[45]=    { 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,//  9bins 
                                 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,       //12 bins
                                 2.2,2.4,2.6,2.8,3.0,                                    //5 bins
  	3.3,3.6,3.9,4.2,4.6,5.0,5.4,5.9, 6.5,7.0,7.5,8.0,8.5,9.2,10.,11.,12.,13.5,15.} ;  //  19
	
	// from Gyula , 28/11/2015

        Double_t gPt13HPtGyu[69] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.45,
                                    0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
                                    1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
                                    2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
                                    4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
                                    11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
                                    26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
	Double_t gPt7TOF[47] = { 0.2,0.25, 0.3,0.35,  0.4,0.45,  0.5,0.55,  0.6,0.65,  0.7,0.75,  0.8, 0.85, 0.9, 0.95, 1.0,
		1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0,
                         3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0 };  //  Barbara TOF Kch


  //
  	fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  	fMultiplicity = new TH1F("fMultiplicity", "multiplicity distribution", 200, 0, 1000);
  	fIncompletEvent = new TH1F("fIncompletEvent", "incomplete event", 200, 0, 1000);
  	fMultpileup = new TH1F("fMultpileup", "multiplicity after pile up rejection event", 200, 0, 1000);
  	fMultV0trigger = new TH1F("fMultV0trigger", "multiplicity after selection of trigger particle ", 200, 0, 1000);
  	fZvertex =  new TH1F("fZvertex", "z vertex position  ", 200, -15, 15);
  	fEventVertex =  new TH1F("fEventVertex", " number of event after the vertex cut  ", 10, 0, 10);
	fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);
	
	fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5,60, -15., 15.0);
  	fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
  	fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
	fRpr = new TH1D("fRpr", "rad distribution  PID pr",100,-10.0, 10.0);
	fdcaToVertexXY = new TH1D("fdcaToVertexXY", "dca to vertex xy distribution ",100,-10.0, 10.0);
	fdcaToVertexXYafterCut = new TH1D("fdcaToVertexXYafterCut", "dca to vertex xy distribution after cut",100,-1.0, 1.0);
  	fptAllKink = new TH1F("fptAllKink", "P_{T} distribution", 50, 0.1, 10.);
	fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);
	fPosiKink= new TH2F("fPosiKink", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
	fQtAll= new TH1F ("fQtAll"," Qt distribution of daughter kink",100, 0, 1.);
	fptKink = new TH1F("fptKink", "P_{T}Kaon Kink  bution",300, 0.0,15.0);
	fQtMothP = new TH2F("fQtMothP", " Qt vrs Mother P", 100, 0., 5.0,100, 0.,0.300);
	fqT1 = new TH1F("fqT1", "Q_{T} distribution",100, 0.0,.300);
	fEta = new TH1F("fEta", " eta distribution of all kinks", 26, -1.3, 1.3);
	fqT2 = new TH1F("fqT2", "Q_{T} distribution of kink daughters without low pt cut",100, 0.0,.300);
	fKinkKaonBackg= new TH1F("fKinkKaonBackg", "P_{T}Kaon kinks background",300, 0.0,15.0);
	
	fPtCut1 = new TH1F ("fPtCut1"," pt distribution after first pt cut (<0.200) ",300, 0.0,15.0);
	fAngMotherPi= new TH2F("fAngMotherPi","Decay angle vrs Mother Mom,Pi",100,0.0,5.0,80,0.,80.);
  	fQtInvM= new TH2F("fQtInvM", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300);
	fInvMuNuAll= new TH1F("fInvMuNuAll", " Inv Mass MuNu all kink",600,0.1,0.7); 
	fInvMassMuNuPtAll =new TH2F("fInvMassMuNuPtAll","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 100, 0.0, 10.0  );
	fRadiusPt =new TH2F("fRadiusPt","radius vs pt  ",80, 90.,250.,100, 0.,10.);
	fAngMotherKC= new TH2F("fAngMotherKC","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
	fkaonInvaiant= new TH1F("fkaonInvaiant","Invar m(kaon) from kink->mu+netrino decay",600,0.10, 0.7);	
	fRadiusNcl = new TH2F("fRadiusNcl","kink radius vrs Nclust,K",75,100.,250., 80,0, 160);
	fPtKPDG = new TH1F("fPtKPDG", "P_{T}Kaon distribution",56, gPtpp );
	fAngMotherKKinks = new TH2F("fAngMotherKKinks","Decay angle vrs Mother Mom,Kinks",300,0.0,15.0,100,0.,100.);
	fPtCut2 = new TH1F("fPtCut2", "P_{T}Kaon distribution",56, gPtpp );
	fPtCut3 = new TH1F("fPtCut3", "P_{T}Kaon distribution",56, gPtpp);
	fTPCSignlMotherK = new TH2F("fTPCSignlMotherK","TPC signal de/dx Mom,K",1000,0.0,20.,150, 0.,300.);
	fPtKaon = new TH1F("fPtKaon", "P_{T} distribution of K all charged",56, gPtpp );	
	fPtKaonP = new TH1F("fPtKaonP", "P_{T} of distribution of K^{+}",56, gPtpp);	
	fPtKaonN = new TH1F("fPtKaonN", "P_{T} distribution of K^{-}",56, gPtpp );	
	 fTPCSignalP = new TH2F("fTPCSignalP","TPC signal de/dx Mom,K",1000,0.0,20.0,150,0.,300.);
	fRadiusNclCln = new TH2F("fRadiusNclCln","kink radius vrs Nclust,K Clean ",75,100.,250., 80,0, 160);
	fRadiusPtcln =new TH2F("fRadiusPtcln","radius vs pt clean ",80, 90.,250.,100, 0.,10.);
	fInvMassMuNuPt =new TH2F("fInvMassMuNuPt","Invariant Kaon = mass-munu  vs pt  ",600, 0.10, 0.7, 100, 0.0, 10.0  );
	fTPCSignlPtpc = new TH2F("fTPCSignlPtpc"," TPC signal de/dx vs Mom TPC,K after all the cut; TPC momentum (P) (GeV/c); dE/dx (a.u.) ",300,0.0,15.0,100, 0., 250.    );
	fMothKinkMomSignl = new TH2F("fMothKinkMomSignl","TPC signal de/dx Mom TPC of Kink; Daughter momentum (GeV/c); dE/dx (a.u.)  ",100,0.0,250.0,100, 0., 250.);
	fTPCSignlKinkDau = new TH2F("fTPCSignlKinkDau","TPC signal daughter kink de/dx Mom,K",500,0.0,10.0,400,0.,250.);
	fTPCMomNSigmaAllKaon = new TH2F("fTPCMomNSigmaAllKaon","TPC signal de/dx Mom TPC,all K  ",300,0.0,15.0,20 , -10., 10.);
	fnSigmaTPC = new TH1F("fnSigmaTPC"," n#sigma distribution of tpc detector" , 30 , -7.5, 7.5);
	fradiurKink = new TH1F("fradiurKink", "radius of  Kink generated",100, 0.,1000.0);
	fLenthKink = new TH1F("fLenthKink", "Length of   K generated",100,0.,1000.);
	fEtaK = new TH1F ("fEtaK", "Eta K distribution", 26,-1.3, 1.3);
	frapiKESD=new TH1F("frapiKESD","rapid Kdistribution", 26,-1.3, 1.3);
	fzVertexPositionKinKvsKinkRad = new TH2F("fzVertexPositionKinKvsKinkRad"," z vertex of kink vs kink radius; z vertex of kink; kink radius",100, -300.0,300.0,100,100., 300.);
	fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",80,-4.,4.0,70,20.,160.);
	fSignPtrapiK = new TH2F("fSignPtrapiK","SignPt vrs rapi ,K",80,-4.0,4.0,30,-1.5,1.5);
	frapiKNcl = new TH2F ("frapiKNcl", "rapi K  vs Number of cluster; rapidity K; Number of cluster", 30,-1.5,1.5, 70,20, 160);
	fSignPt= new TH1F("fSignPt","SignPt ,K",80,-4.0,4.0);
	fChi2NclTPC= new TH2F("fChi2NclTPC","Chi2vrs nclust,K",100,0.,500., 70,20, 160);
	fRatioChi2Ncl= new TH1F("fRatioChi2Ncl","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
	flifetime= new TH1F("flifetime", "ct study of K-kinks",100,0.,1000.);
	fPtKinkKaon = new TH1F ("fPtKinkKaon","p_{T} Kaon kinks identied",56, gPtpp );
	fDCAkink = new TH1F("fDCAkink ", "DCA kink vetrex ",50, 0.0,1.0);
	 fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  distribution, counts",56, gPtpp );
	 fPtKinkPos= new TH1F("fPtKinkPos", "Pos P_{T} K^{+} Kink distribution, counts",56, gPtpp);
	 fPtKinkNeg= new TH1F("fPtKinkNeg", "Pos P_{T} K^{-} Kink   distribution, counts",56, gPtpp);
	fPtKinkK0= new TH1F("fPtKinkK0", "P_{T}Kaon Kink  distribution, counts",44, gPt13K0PKal);
	fPtKinkK0P= new TH1F("fPtKinkK0P", "P_{T} KPl Kink  distribution, counts",44, gPt13K0PKal);
  	fPtKinkK0N= new TH1F("fPtKinkK0N", "P_{T} KMn Kink  distribution, counts",44, gPt13K0PKal);
	fPtKinkGyu= new TH1F("fPtKinkGyu", "P_{T}Kaon Kink  distribution, counts",68, gPt13HPtGyu);
  	fPtKinkGyuP= new TH1F("fPtKinkGyuP", "P_{T} KPl Kink  distribution, counts",68, gPt13HPtGyu);
  	fPtKinkGyuN= new TH1F("fPtKinkGyuN", "P_{T} KMn Kink  distribution, counts",68, gPt13HPtGyu);
	fKinKRbn= new TH1F("fKinKRbn", "p_{t}Kaon kinks identi[GeV/c],Entries",46,gPt7TOF);	
	fradPtRpDt=new TH3F("fradPtRpDt","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
	fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
	fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  	fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  	fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);	
	fNumberOfEvent = new TH1F("fNumberOfEvent", "the number of events in this run", 40, 0., 20.);
	
	fMultMC_wo_any_cut = new TH1F ("fMultMC_wo_any_cut"," number of tracks in MC", 100, 0., 50000.);
	fMultMC_incompleteDAQ = new TH1F ("fMultMC_incompleteDAQ"," number of tracks in MC after the cut of Incomplete DAQ", 100, 0., 50000.);
	fMultMC_AfterPileUp = new TH1F ("fMultMC_AfterPileUp"," number of tracks in MC after the pile up cut", 100, 0., 50000.);
	fMultTriggerMCAfterV0 = new TH1F ("fMultTriggerMCAfterV0"," number of tracks in MC after v0 trigger selection", 100, 0., 50000.);
	frapidKMC = new TH1F("frapidKMC ", "rapidity distributon of kaon before any cut  MC  ",26,-1.3, 1.3);
	fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",56, gPtpp);
	fSignPtGen= new TH1F("fSignPtGen","SignPtGen ,K",100,-5.0,5.0);
	fPtKPlMC= new TH1F("fPtKPlMC", "P_{T} Kaon Pos  generated", 56, gPtpp );
	fPtKMnMC= new TH1F("fPtKMnMC", "P_{T} Kaon Minus generated",56, gPtpp );
	flengthTrackRef =new TH1F("flengthTrackRef","lifetime ref K   Decay   ",100,0.,1000.0);
	flifeSmallMC=new TH1F("flifeSmallMC","lifetime ref K   Decay   ",100,0.,1000.0);
  	flifetimeMC =new TH1F("flifetimeMC","lifetime ref K   Decay   ",100,0.,1000.0);
	flifeInt =new TH1F("flifeInt", "lifetime ref K   Decay   ",100,0.,1000.0);
	flifeYuri=new TH1F("flifeYuri","lifetime ref K   Decay   ",100,0.,1000.0);
	flifetiMCK=new TH1F("flifetiMCK", "lifetime ref K   Decay   ",100,0.,1000.0);
	 flenYuri=new TH1F("flenYuri","lifetime ref K   Decay   ",100,0.,1000.0);
	flengthMCK=new TH1F("flengthMCK", "length of K  MCref decay ",100,0.,1000.0);
	flifetime_MCprocess4 =new TH1F("flifetime_MCprocess4", "lifetime ref K   Decay   ",100,0.,1000.0);
	fradPtRapMC=new TH3F("fradPtRapMC","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
	flifetime_kaonpionPDG =new TH1F("flifetime_kaonpionPDG", "lifetime ref K   Decay   ",100,0.,1000.0);
	flifetime_kaonmuonPDG =new TH1F("flifetime_kaonmuonPDG", "lifetime ref K   Decay   ",100,0.,1000.0);
	fmaxAngMomKmuMC= new TH2F("fmaxAngMomKmuMC","Decay angle vrs Mother Mom,Kmu",100,0.0,10.0,120,0.,120.);
	fradPtRapDC=new TH3F("fradPtRapDC","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
	  fradMC= new TH1F("fradMC", "radius  K generated",100,0.,1000.);
	fQtKMuMC= new TH1F("fQtKMuMC", "Q_{T} distribution  K to mu MC",100, 0.0,.500);
	fgenPtEtR= new TH1F("fgenPtEtR", "P_{T}Kaon distribution", 56, gPtpp);
	fgenPtEtRP= new TH1F("fgenPtEtRP", "P_{T}Kaon distribution positive", 56, gPtpp );
	 fMCEtaKaon = new TH1F("fMCEtaKaon"," Hist of Eta K -Kink Selecied",26,-1.3,1.3);
	 fSignPtEtaMC= new TH2F("fSignPtEtaMC","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
	 fSignPtMC= new TH1F("fSignPtMC","SignPt ,K",100,-5.,5.0);
	fgenPtEtRN= new TH1F("fgenPtEtRN", "P_{T}Kaon distribution positive", 56, gPtpp  );
	fQtKElMC= new TH1F("fQtKElMC", "Q_{T} distribution   K to elec MC",100, 0.0,.500);
	  fQtKPiMC= new TH1F("fQtKPiMC", "Q_{T} distribution K to pi MC",100, 0.0,.500);
	 fESDMult= new TH1F("fESDMult", " ESD charge multipliESD",100, 0.0,300.);
	fFakepipi = new TH1F("fFakepipi", "P_{T}fake pipi   ", 56, gPtpp);
        fFakeKPi = new TH1F("fFakeKPi", "P_{T}fake Kpi   ", 56, gPtpp   );
	fRadiusPtPion  =new TH2F("fRadiusPtPion","radius vs pt Pion PDG ",80, 90.,250.,100, 0.,10.);
	fRadiusPtKaon  =new TH2F("fRadiusPtKaon","radius vs pt Kaon PDG ",80, 90.,250.,100, 0.,10.);
	fQtKMu= new TH1F("fQtKMu", "Q_{T} distribution  K to mu ",100, 0.0,.500);
        fQtKPi= new TH1F("fQtKPi", "Q_{T} distribution K to pi",100, 0.0,.500);
        fQtKEl= new TH1F("fQtKEl", "Q_{T} distribution   K to elec",100, 0.0,.500);
	fQtK3PiP= new TH1F("fQtK3PiP", "Q_{T} distribution K to 3pi ",100, 0.0,.500);
  	fQtK3PiM= new TH1F("fQtK3PiM", "Q_{T} distribution K to 3pi ",100, 0.0,.500);
	fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",56, gPtpp );
	fHiPtKPDGP= new TH1F("fHiPtKPDGP", "P_{T}Kaon Pos ESD", 56, gPtpp );
        fHiPtKPDGN= new TH1F("fHiPtKPDGN", "P_{T}Kaon neg ESD", 56, gPtpp );
	fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3);
	frapidESDK= new TH1F("frapidESDK", "rapidity distribution", 26,-1.3, 1.3);
	fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300);
	fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,10.0,80,0.,80.);
	fMinvPi= new TH1F("fMinvPi","Invar m(kaon) from kink-> decay",100,0.0, 1.2);
        fMinvKa= new TH1F("fMinvKa","Invar m(kaon) from kink-> decay",100,0.0, 2.0);
	fcodeH   = new TH2F("fcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fZkinkZDau = new TH2F("fZkinkZDau", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
	fRadiusPtFake  =new TH2F("fRadiusPtFake","radius vs pt Pion Fake ",80, 90.,250.,100, 0.,10.);
	fTPCMomNSgnl = new TH2F("fTPCMomNsgnl","TPC signal de/dx Mom TPC,K  ",100,0.0,8.0,20 , -10., 10.);
	fPtPrKink=new TH1F("fPtPrKink","pt of ESD  kaonKink tracks",56, gPtpp);
	flifTiESDK=new TH1F("flifTiESDK","lifetime ref K   Decay   ",100,0.,1000.0);
	fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi", 56, gPtpp);
	 fkinkKaonP= new TH1F("fKinkKaonP", "P_{T}Kaon distribution positive", 56, gPtpp );
	 fkinkKaonN= new TH1F("fKinkKaonN", "P_{T}Kaon distribution positive", 56, gPtpp);
	fcode2   = new TH2F("fcode2", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fTPCSgnlPtpc = new TH2F("fTPCSgnlPtpc","TPC signal de/dx Mom TPC,K  ",100,0.0,8.0,100, 0., 250.    );
	fMothKinkMomSgnlD = new TH2F("fMothKinkMomSgnlD","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.,100, 0., 250.);
	fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",56, gPtpp);
	fMothKinkMomSgnl  = new TH2F("fMothKinkMomSgnl","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.,100, 0., 250.    );
	 fcodeDau2   = new TH2F("fcodeDau2", "code vrs dcode dist. kinks,K",100,0.,500.,100,0.,500.);
	fTPCSgnlKinkDau = new TH2F("fTPCSgnlKinkDau","TPC signal de/dx Mom,K",100,0.0,8.0,100,0.,250.);
	fMinvPr= new TH1F("fMinvPr","Invar m(kaon) from kink-> decay",100,0.0, 1.2);
	fDCAkinkBG = new TH1F("fDCAkinkBG", "DCA kink vetrex ",50, 0.0,1.0);
	 fPosiKinKBgXY= new TH2F("fPosiKinKBgXY", "Y vrx X kink Vrexbg ",100, -300.0,300.0,100, -300, 300.);
	 fPosiKinKBgZY= new TH2F("fPosiKinKBgZY", "Y vrx Z kink Vrexbg ",100, -300.0,300.0,100, -300, 300.);
	fPosiKinKBgZX= new TH2F("fPosiKinKBgZX", "X vrx Z kink Vrexbg ",100, -20.0,20.0,100, 0., 300.);
	fKinKBGP  = new TH1F("fKinKBGP  ", "P_{T}Kaon Pos ESD", 56, gPtpp);
	fKinKBGN= new TH1F("fKinKBGN", "P_{T}Kaon neg ESD",56, gPtpp );
	fdcodeH = new TH2F("fdcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fcode4   = new TH2F("fcode4", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fNumberOfEvent_cent = new TH1F("fNumberOfEvent_cent", "the number of events in this run", 10, 0., 6.);
	fEventVsCentrality = new TH2F("fEventVsCentrality","centrality vs event", 100, 0., 100., 10, 0.,  5.);

	 fOutputList->Add(fHistPt);
	 fOutputList->Add(fMultiplicity);
	 fOutputList->Add(fIncompletEvent);
	 fOutputList->Add(fMultpileup);
	 fOutputList->Add(fMultV0trigger);
	 fOutputList->Add(fZvertex);
	 fOutputList->Add(fEventVertex);
	 fOutputList->Add(fRatioCrossedRows);
	 fOutputList->Add(fZvXv);
	 fOutputList->Add(fZvYv);
	 fOutputList->Add(fXvYv);
	 fOutputList->Add(fRpr);
	 fOutputList->Add(fdcaToVertexXY);
	 fOutputList->Add(fdcaToVertexXYafterCut);
	 fOutputList->Add(fptAllKink);
	 fOutputList->Add(fRatioCrossedRowsKink);
	 fOutputList->Add(fPosiKink);
	 fOutputList->Add(fQtAll);
	 fOutputList->Add(fptKink);
	 fOutputList->Add(fQtMothP);
	 fOutputList->Add(fqT1);
	 fOutputList->Add(fEta);
	 fOutputList->Add(fqT2);
	 fOutputList->Add(fKinkKaonBackg);
	 fOutputList->Add(fPtCut1);
	 fOutputList->Add(fAngMotherPi);
	 fOutputList->Add(fQtInvM);
	 fOutputList->Add(fInvMuNuAll);
	 fOutputList->Add(fInvMassMuNuPtAll);
	 fOutputList->Add(fRadiusPt);
	 fOutputList->Add(fAngMotherKC);
	 fOutputList->Add(fkaonInvaiant);
	 fOutputList->Add(fRadiusNcl);
	 fOutputList->Add(fPtKPDG);
	 fOutputList->Add(fAngMotherKKinks);
	 fOutputList->Add(fPtCut2);
	 fOutputList->Add(fPtCut3);
	 fOutputList->Add(fTPCSignlMotherK);
	 fOutputList->Add(fPtKaon);
	 fOutputList->Add(fPtKaonP);
	 fOutputList->Add(fPtKaonN);
	 fOutputList->Add(fTPCSignalP);
	 fOutputList->Add(fRadiusNclCln);
	 fOutputList->Add(fRadiusPtcln);
	 fOutputList->Add(fInvMassMuNuPt);
	 fOutputList->Add(fTPCSignlPtpc);
	 fOutputList->Add(fMothKinkMomSignl);
	 fOutputList->Add(fTPCSignlKinkDau);
	 fOutputList->Add(fTPCMomNSigmaAllKaon);
	 fOutputList->Add(fnSigmaTPC);
	 fOutputList->Add(fradiurKink);
	 fOutputList->Add(fLenthKink);
	 fOutputList->Add(fEtaK);
	 fOutputList->Add(frapiKESD);
	 fOutputList->Add(fzVertexPositionKinKvsKinkRad);
	 fOutputList->Add(fSignPtNcl);
	 fOutputList->Add(fSignPtrapiK);
	 fOutputList->Add(frapiKNcl);
	 fOutputList->Add(fSignPt);
	 fOutputList->Add(fChi2NclTPC);
	 fOutputList->Add(fRatioChi2Ncl);
	 fOutputList->Add(flifetime);
	 fOutputList->Add(fPtKinkKaon);
	 fOutputList->Add(fDCAkink);
	 fOutputList->Add(fPtKink);
	 fOutputList->Add(fPtKinkPos);
	 fOutputList->Add(fPtKinkNeg);
	 fOutputList->Add(fPtKinkK0);
	 fOutputList->Add(fPtKinkK0P);
	 fOutputList->Add(fPtKinkK0N);
	 fOutputList->Add(fPtKinkGyu);
	 fOutputList->Add(fPtKinkGyuP);
	 fOutputList->Add(fPtKinkGyuN);
	 fOutputList->Add(fKinKRbn);
	 fOutputList->Add(fradPtRpDt);
	 fOutputList->Add(fAngMomK);
	 fOutputList->Add(fPosiKinkK);
	 fOutputList->Add(fPosiKinKXZ);
	 fOutputList->Add(fPosiKinKYZ);
	 fOutputList->Add(fNumberOfEvent);
	 fOutputList->Add(fMultMC_wo_any_cut);
	 fOutputList->Add(fMultMC_incompleteDAQ);
	 fOutputList->Add(fMultMC_AfterPileUp);
	 fOutputList->Add(fMultTriggerMCAfterV0);
	 fOutputList->Add(frapidKMC);
	 fOutputList->Add(fptKMC);
	 fOutputList->Add(fSignPtGen);
	 fOutputList->Add(fPtKPlMC);
	 fOutputList->Add(fPtKMnMC);
	 fOutputList->Add(flengthTrackRef);
	 fOutputList->Add(flifetimeMC);
	 fOutputList->Add(flifeSmallMC);
	 fOutputList->Add(flifeInt);
	 fOutputList->Add(flifeYuri);
	 fOutputList->Add(flifetiMCK);
	 fOutputList->Add(flenYuri);
	 fOutputList->Add(flengthMCK);
	 fOutputList->Add(flifetime_MCprocess4);
	 fOutputList->Add(fradPtRapMC);
	 fOutputList->Add(flifetime_kaonpionPDG);
	 fOutputList->Add(flifetime_kaonmuonPDG);
	 fOutputList->Add(fmaxAngMomKmuMC);
	 fOutputList->Add(fradPtRapDC);
	 fOutputList->Add(fradMC);
	 fOutputList->Add(fQtKMuMC);
	 fOutputList->Add(fgenPtEtR);
	 fOutputList->Add(fgenPtEtRP);
	 fOutputList->Add(fMCEtaKaon);
	 fOutputList->Add(fSignPtEtaMC);
	 fOutputList->Add(fSignPtMC);
	 fOutputList->Add(fgenPtEtRN);
	 fOutputList->Add(fQtKElMC);
	 fOutputList->Add(fQtKPiMC);
	 fOutputList->Add(fESDMult);
	 fOutputList->Add(fFakepipi);
	 fOutputList->Add(fFakeKPi);
	 fOutputList->Add(fRadiusPtPion);
	 fOutputList->Add(fRadiusPtKaon);
	 fOutputList->Add(fQtKMu);
	 fOutputList->Add(fQtKPi);
	 fOutputList->Add(fQtKEl);
	 fOutputList->Add(fQtK3PiP);
	 fOutputList->Add(fQtK3PiM);
	 fOutputList->Add(fHistPtKPDG);
	 fOutputList->Add(fHiPtKPDGP);
	 fOutputList->Add(fHiPtKPDGN);
	 fOutputList->Add(fHistEta);
	 fOutputList->Add(frapidESDK);
	 fOutputList->Add(fHistQt2);
	 fOutputList->Add(fAngMomPi);
	 fOutputList->Add(fMinvPi);
	 fOutputList->Add(fMinvKa);
	 fOutputList->Add(fcodeH);
	 fOutputList->Add(fZkinkZDau);
	 fOutputList->Add(fRadiusPtFake);
	 fOutputList->Add(fTPCMomNSgnl);
	 fOutputList->Add(fPtPrKink);
	 fOutputList->Add(flifTiESDK);
	 fOutputList->Add(fKinkKaon);
	 fOutputList->Add(fkinkKaonP);
	 fOutputList->Add(fkinkKaonN);
	 fOutputList->Add(fcode2);
	 fOutputList->Add(fTPCSgnlPtpc);
	 fOutputList->Add(fMothKinkMomSgnlD);
	 fOutputList->Add(fKinkKaonBg);
	 fOutputList->Add(fMothKinkMomSgnl);
	 fOutputList->Add(fcodeDau2);
	 fOutputList->Add(fTPCSgnlKinkDau);
	 fOutputList->Add(fMinvPr);
	 fOutputList->Add(fDCAkinkBG);
	 fOutputList->Add(fPosiKinKBgXY);
	 fOutputList->Add(fPosiKinKBgZY);
	 fOutputList->Add(fPosiKinKBgZX);
	 fOutputList->Add(fKinKBGP);
	 fOutputList->Add(fKinKBGN);
	 fOutputList->Add(fdcodeH);
	 fOutputList->Add(fcode4);
	 fOutputList->Add(fNumberOfEvent_cent);
	 fOutputList->Add(fEventVsCentrality);
	

PostData(1, fOutputList);
}


//___________________user exec_____________________________________________________
void AliAnalysisTaskKinkpp5TeVMC::UserExec(Option_t *) 
{
  //initialized values
	Float_t dca1[2], cov1[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
	





	AliVEvent *event = InputEvent();
  	if (!event) {
     	Printf("ERROR: Could not retrieve event");
     	return;
  	}

  	AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  	if (!esd) {
     	Printf("ERROR: Could not retrieve esd");
     	return;
  	}

	AliMCEvent* mcEvent = MCEvent();
  	if (!mcEvent) {
     	Printf("ERROR: Could not retrieve MC event");
     	return;
  	}
	fNumberOfEvent->Fill(1.5);
//	cout<<" hi the MC event is found here***********************"<<endl;
	fMultMC_wo_any_cut->Fill(mcEvent->GetNumberOfTracks() );
	
// Number ESD tracks 
 	Int_t nESDTracks =  esd->GetNumberOfTracks();
        //fMultiplicity->Fill(nESDTracks);
        fMultiplicity->Fill(2);

        fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0105 + 0.0350/pt^1.01");
        //fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
         fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
	fNumberOfEvent->Fill(2.5);
//physics selection
	UInt_t maskIsSelected =
        ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
        Bool_t isSelected = 0;
        isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
         if (!isSelected) {
         ( PostData(1, fOutputList));
            return;
         }
        fMultV0trigger->Fill(2);
	fNumberOfEvent->Fill(3.5);


//multiplicity/ centrality 

        Float_t cent = -999;
        if (esd->GetRunNumber() < 244824) { //OLD multiplicity/centrality class framework
        AliCentrality *centrality = esd->GetCentrality();
        cent = centrality->GetCentralityPercentile("V0M");
        } else { //New multiplicity/centrality class framework
        AliMultSelection *fMultSel = (AliMultSelection *) esd->FindListObject("MultSelection");
        if (!fMultSel) {
  //If you get this warning please check that the AliMultSelectionTask actually ran (before your task) 
        AliWarning("AliMultSelection object not found!");
        } else {
  //Event selection is embedded in the Multiplicity estimator so that the Multiplicity percentiles are well defined and refer to the same sample
        cent = fMultSel->GetMultiplicityPercentile("V0M", kTRUE);
//        if ((cent < 0) || (cent > 100)) return; //Event selection
        }
        }

        fNumberOfEvent->Fill(4.5);
        Double_t mctrack= mcEvent->GetNumberOfTracks();

        fMultTriggerMCAfterV0->Fill(mctrack );
        fMultV0trigger->Fill(2);

// check incomplete event
        if (esd->IsIncompleteDAQ()) return;
        //fIncompletEvent ->Fill(esd->GetNumberOfTracks() );
        fIncompletEvent ->Fill(2 );
        fNumberOfEvent->Fill(5.5);
// check of Pileup   
        if (esd->IsPileupFromSPD()) return;
       //Multpileup->Fill(nESDTracks);
        fMultpileup->Fill(2);
        fNumberOfEvent->Fill(6.5);

//tracklet vs cluster cut
        AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
        Double_t IsCluVstrk = AnalysisUtils->IsSPDClusterVsTrackletBG(esd);
        if(IsCluVstrk)
          return;
	fNumberOfEvent->Fill(7.5);
//TPC or SPD vertex check
        const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
        const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
        Bool_t hasSPD = spdVertex->GetStatus();
        Bool_t hasTrk = trkVertex->GetStatus();
        //Note that AliVertex::GetStatus checks that N_contributors is > 0
        if (!(hasSPD && hasTrk)) return;

 	fNumberOfEvent->Fill(8.5);
// vertex cut
        //const AliESDVertex *vertex= esd->GetPrimaryVertex();    
        //if(!vertex) return;

        const AliESDVertex *vertex=GetEventVertex(esd);
        if(!vertex) return;
	fNumberOfEvent->Fill(9.5);
        Double_t vpos[3];
        vertex->GetXYZ(vpos);
        fZvertex->Fill(vpos[2]);
        if (TMath::Abs( vpos[2] ) > 10. ) return;
        fEventVertex->Fill(2);
        fNumberOfEvent->Fill(10.5);

        fEventVsCentrality->Fill(cent, 2);
	
  	Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
  	AliStack* stack=mcEvent->Stack();
//primary tracks  in MC
        Int_t  nPrimary = stack->GetNprimary();
	// variable to count tracks
  	Int_t nMCTracks = 0;
  	Int_t nMCKinkKs = 0;

//if (cent>=0 && cent<5) {
	fNumberOfEvent_cent->Fill(1);	
	for (Int_t iMc = 0; iMc < mcEvent->GetNumberOfTracks(); iMc++)
  	{

    	TParticle* particle = stack->Particle(iMc);

	//Int_t labelmc = particle->GetLabel();
// keep only primaries
  	if (!stack->IsPhysicalPrimary(iMc)) continue;
    	if (!particle)
    	{
      	continue;
    	}
        Double_t ptK = particle->Pt();
	if( ptK <0.200) continue;  
	Double_t EtaMC  = particle->Eta();
	if ((TMath::Abs(EtaMC)) > 0.8) continue ;
	Float_t charge=0;
	Float_t code = particle->GetPdgCode();
        Int_t  mcProcess=-1011;

	if ((code==321)||(code==-321)){
	 Double_t   etracKMC= TMath::Sqrt(particle->P() *particle->P()  + 0.493677 *0.493677  );
         Double_t rapidiKMC = 0.5 * (TMath::Log(  (etracKMC +particle->Pz())/( etracKMC-particle->Pz() )) )  ;
	if ( (TMath::Abs( rapidiKMC)) > fRapiK ) continue;   // 
         frapidKMC ->Fill(rapidiKMC) ;  
	 if(code > 0 ) charge =1;
         if(code < 0 ) charge =-1;
         Float_t chargePt= charge*ptK;
	 fptKMC->Fill(ptK);
         fSignPtGen->Fill(chargePt);// kaon gensign pt
         if (charge==1 )  fPtKPlMC->Fill( ptK);
         if ( charge==-1) fPtKMnMC->Fill( ptK );

	 Double_t mVz=particle->Vz();
	 TClonesArray* trArray=0;
         TParticle* tempParticle=0;
	TVector3 DecayMomentumK(0,0,0);
         Float_t lengthKMC=0;
	if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) {
        AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
        lengthKMC = MCtrackReference->GetLength();
	}
	flengthTrackRef ->Fill(lengthKMC);
	flifetimeMC ->Fill(  (lengthKMC*0.493667  /particle->P()));  // 19/7
       	if ((lengthKMC>100.)&& (lengthKMC<300.) )  flifeSmallMC->Fill( (lengthKMC*0.493667/particle->P() ) );	

	Int_t nMCKpi =0;
       	Int_t mcProc4 =0;
       	Int_t mcProc13=0;
        Float_t radiusD=0;
   //    Double_t lengthK =0.;
       	Double_t LengthK =0.;
       	Double_t lenYuri =0.;
        Float_t MCQt =0.;
//        Double_t MCQt3[2];
      	Int_t firstD=particle->GetFirstDaughter();
       	Int_t lastD=particle->GetLastDaughter();
	//cout << " last daughter ======"<<lastD<<endl;
	if( (lastD<=0) || (firstD<=0)  ) continue;

	if ( lastD > mcEvent->GetNumberOfTracks() ) continue;
        if (firstD > mcEvent->GetNumberOfTracks() ) continue;
//loop on secondaries
        for (Int_t k=firstD;k<=lastD;k++) {
        if ( k > 0 )    {
        TParticle*    daughter1=stack->Particle(k);     
        Float_t dcode = daughter1->GetPdgCode();
//	    mother momentum trackrefs    and QtMC     
        if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) {
        AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
        DecayMomentumK.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
        }
        const TVector3 MCP3d(daughter1->Px(), daughter1->Py(), daughter1->Pz()); //MC daughter's momentum
        MCQt = MCP3d.Perp(DecayMomentumK); //MC daughter's transverse momentum in mother's frame
	Double_t MCKinkAngle = TMath::ASin(MCQt/daughter1->P() );
        Double_t  MCKinkAngle2= TMath::RadToDeg() * MCKinkAngle; // in degrees 
	mcProcess=daughter1->GetUniqueID();
	//cout <<"m process==============="<<mcProcess<<endl;

//        cout <<"mc prccess unique i======"<<mcProcess<< endl;


	radiusD=daughter1->R();
	Double_t hVz=daughter1->Vz();
	 LengthK = TMath::Sqrt( radiusD*radiusD  + ( mVz-hVz) * (mVz-hVz) ); 
	 lenYuri  = (TMath::Abs( mVz-hVz))* (TMath::Sqrt( 1.+ ( ptK*ptK)/ (particle->Pz() * particle->Pz()) )) ;
	if(mcProcess==13) {
    	mcProc13++;

      	if(mcProc13==1)     flifeInt  ->Fill(  (lengthKMC*0.493667  /particle->P()));  
	}

	if (mcProcess==4) {

    	mcProc4++;
   	if ( mcProc4==1)  {
        flifeYuri ->Fill(  (lenYuri  *0.493667  /particle->P()));  // 19/7
        flifetiMCK->Fill(LengthK);
        flenYuri  ->Fill(lenYuri);
   	flengthMCK->Fill(lengthKMC);  //
        flifetime_MCprocess4 ->Fill(  (lengthKMC*0.493667  /particle->P()));  // 19/7
        fradPtRapMC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
  	} 

	if (  ( ( code==321 )&&( dcode ==-13  ))||( ( code ==-321)&&(dcode== 13) ) || ( ( code==321 )&&( dcode ==-11  )) || ( (code ==-321)&&(dcode== 11))) {
                      flifetime_kaonmuonPDG ->Fill(  (lengthKMC*0.493667  /particle->P()));
	if( (radiusD >fKinkRadLow )&&( radiusD< fKinkRadUp) )
                 fmaxAngMomKmuMC->Fill(particle->P() , MCKinkAngle2);// MC 
              }
	if (( (TMath::Abs(code)==321 )&&(TMath::Abs(dcode)  ==211  ))&& ( mcProc4<2)) flifetime_kaonpionPDG->Fill( lengthKMC *0.493667 /particle->P()) ;

	if(MCKinkAngle2 < 2.) continue;  // as in ESD 
	if (((daughter1->R())> fKinkRadLow )&&((daughter1->R())< fKinkRadUp )&& (MCQt>fLowQt)  ){

	if ( ( code==321 )&&( dcode ==-13  ))   {
                    
	fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
        fradMC->Fill(daughter1->R());
        fQtKMuMC ->Fill(MCQt );//to muon
        fgenPtEtR->Fill( ptK);//to muon
        fgenPtEtRP->Fill( ptK );//to muon
        fMCEtaKaon  ->Fill(rapidiKMC );//to muon
        fSignPtEtaMC->Fill(ptK,rapidiKMC );//to muon
        fSignPtMC->Fill(ptK);//to muon
        flifetiMCK->Fill(lenYuri*0.4933667/particle->P()  );

//  flifetiMCK->Fill(   LengthK*0.4933667/   ptK        
    } //  positive kaon 
	if (  ( code ==-321)&&(dcode== 13)){
	fgenPtEtR->Fill(  ptK );//to muon
        fQtKMuMC ->Fill(MCQt );//to muon
        fgenPtEtRN->Fill(ptK);  //  
        fSignPtEtaMC->Fill(chargePt,rapidiKMC );//to muon
        fMCEtaKaon  ->Fill(rapidiKMC );//to muon
        fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
        fradMC->Fill(daughter1->R());
        fSignPtMC->Fill(chargePt);//to muon
        flifetiMCK->Fill(   lenYuri*0.4933667/particle->P() ) ;


	} // negative kaoni
	if ( ( code==321 )&&( dcode ==-11  ))   {
            fQtKElMC ->Fill(MCQt );//to muon
            fgenPtEtR->Fill( ptK);//to electron
            fgenPtEtRP->Fill( ptK);//to muon
         fMCEtaKaon  ->Fill(rapidiKMC );//to electron
            fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
                  fradMC->Fill(daughter1->R());
         fSignPtEtaMC->Fill(ptK,rapidiKMC );//to electron
         fSignPtMC->Fill(ptK);
                      flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
        //             flifetiMCK->Fill(   LengthK*0.4933667/ ptK      );

    } //  positive kaon

	if (  ( code ==-321)&&(dcode== 11)){
        fgenPtEtR->Fill(   ptK);//to electron
        fQtKElMC ->Fill(MCQt );//to muon
        fgenPtEtRN->Fill(ptK);  //  
       	fSignPtEtaMC->Fill(chargePt,rapidiKMC  );//to electron
       	fMCEtaKaon  ->Fill(rapidiKMC  );//to electron
       	fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
       	fradMC->Fill(daughter1->R());
       	fSignPtMC->Fill(chargePt);//to electron 
       	flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );

    	} //  negative code
	 if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211)))    nMCKpi++ ;
        if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211)))   {
        if ( nMCKpi > 0) {
//              MCQt3[nMCKpi-1] = MCQt ;// k to pipipi 
	}
        }
    nMCKinkKs++;

} // end of the radius and low qt cut
} // end of the mcprocess=4	
} //k>0 loop
} // daughter loop
//cout<<" nMCKpi=========="<<nMCKpi<<endl;

  if ( code > 0) {
              if( nMCKpi == 1) fgenPtEtR->Fill(ptK);  //  k to pipi
              if( nMCKpi == 1) fgenPtEtRP->Fill(ptK);  //  k to pipi
              if( nMCKpi == 1) fSignPtEtaMC->Fill(ptK,rapidiKMC );  //  k to pipi
              if( nMCKpi == 1) fSignPtMC->Fill(ptK);  //  k to pipi
              if( nMCKpi == 1) fMCEtaKaon->Fill(rapidiKMC );  //  k to pipi
              if(nMCKpi==1) fradMC->Fill(radiusD       );
              if(nMCKpi==1) fQtKPiMC->Fill( MCQt   );
              if (nMCKpi==1)  fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
              if(nMCKpi==1)     flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
             //   if(nMCKpi==1)     flifetiMCK->Fill(   LengthK*0.4933667/      ptK      );
   }   //positive kaon 
     //
           if ( code < 0) {

              if( nMCKpi == 1) fgenPtEtR->Fill(   ptK );  //  k to pipi
              if( nMCKpi == 1) fgenPtEtRN->Fill(ptK);  //  k to pipi
              if( nMCKpi == 1) fSignPtEtaMC->Fill(chargePt,rapidiKMC );  //  k to pipi
              if( nMCKpi == 1) fSignPtMC->Fill(chargePt);  //  k to pipi
              if( nMCKpi == 1) fMCEtaKaon->Fill(rapidiKMC );  //  k to pipi
              if(nMCKpi==1) fradMC->Fill(radiusD       );
              if(nMCKpi==1) fQtKPiMC->Fill( MCQt   );
              if( nMCKpi== 1) fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
              if(nMCKpi==1)     flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
         //       if(nMCKpi==1)     flifetiMCK->Fill(   LengthK*0.4933667/      ptK      );

   }   //negative K



} // kaon pid 
} // end of the MC loop


//*************************************************************//
//MC reconstructed
	Double_t nsigmaPion =-100.0;
  //     Double_t nsigmaDau  =-100.0;
       	Double_t dEdxKinkDau =0.0;
	Double_t nsigmall = 100.0;
       	Double_t nsigma = 100.0;
	Double_t dEdxDauMC  =   0.0;
	Double_t vtrack[3], ptrack[3];
	
   	Int_t nGoodTracks =  esd->GetNumberOfTracks();
	fESDMult->Fill(nGoodTracks);




// loop on kink daughters
   	for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {

    	AliESDtrack* trackD = esd->GetTrack(iTrack);
    	if (!trackD) {
      	Printf("ERROR: Could not receive track %d", iTrack);
      	continue;
    	}
//	
        Int_t indexKinkDau=trackD->GetKinkIndex(0);
// daughter kink 
//        AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkDau)-1);
      	if ( indexKinkDau > 0 )    {
    	Int_t labelD = trackD->GetLabel();
    	labelD = TMath::Abs(labelD);
  //       mss 2015 nsigmaPion     = (fPIDResponse->NumberOfSigmasTPC(trackD  , AliPID::kPion));// 26/10 eftihis
   	nsigmaPion     = (fPIDResponse->NumberOfSigmasTPC(trackD  , AliPID::kPion));// 26/10 eftihis
     	dEdxKinkDau =  (trackD->GetTPCsignal()  )  ;  //  daughter kink  dEdx 
//   KinkDauCl=(trackD->GetTPCclusters(0)  )  ;
	fTPCSignlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
        dEdxDauMC   = trackD->GetTPCsignal()     ;  //  daughter kink

 }
   }





	for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {

    	AliESDtrack* track = esd->GetTrack(iTrack);
    	if (!track) {
      	Printf("ERROR: Could not receive track %d", iTrack);
      	continue;
    	}
	if (!fESDtrackCuts->AcceptTrack(track)) continue;

	fHistPt->Fill(track->Pt());

	  //    sigmas
    	 nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
      	nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));

// checking of tpc cros rows

	Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  	Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  	if (track->GetTPCNclsF()>0) {
    	ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
    	fRatioCrossedRows->Fill(ratioCrossedRowsOverFindableClustersTPC);
   	}
// kink index
	Int_t indexKinkPos=track->GetKinkIndex(0);   // kink index 

      	Int_t tpcNCl = track->GetTPCclusters(0);
      	Double_t tpcSign = track->GetSign();

    	Int_t label = track->GetLabel();
    	label = TMath::Abs(label);

	 if(label > mcEvent->GetNumberOfTracks()) continue; //  
	
    	TParticle * part = stack->Particle(label);
    	if (!part) continue;
	// loop only on Primary tracks
     	if (label > nPrimary) continue; /// primary tracks only   
    
	//    pt cut 
      if ( (track->Pt())<.200)continue;

	UInt_t status=track->GetStatus();

    	if((status&AliESDtrack::kITSrefit)==0) continue;
    	if((status&AliESDtrack::kTPCrefit)==0) continue;
      	if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;

      	Double_t extCovPos[15];
      	track->GetExternalCovariance(extCovPos);


    	track->GetXYZ(vtrack);
	fXvYv->Fill(vtrack[0],vtrack[1]);
 	fZvYv->Fill(vtrack[1],vtrack[2]);
 	fZvXv->Fill(vtrack[0],vtrack[2]);

// track momentum, rapidity calculation
     	track->GetPxPyPz(ptrack);

    	TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);

// K-rapidity calcualtion 
         Double_t   etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  );
         Double_t rapiditK = 0.5 * (TMath::Log(  (etracK + ptrack[2]  ) / ( etracK - ptrack[2])  ))  ;

    	Double_t trackEta=trackMom.Eta();
    	Double_t trackPt = track->Pt();

	Float_t bpos[2];
    	Float_t bCovpos[3];
    	track->GetImpactParameters(bpos,bCovpos);

    	if (bCovpos[0]<=0 || bCovpos[2]<=0) {
     	Printf("Estimated b resolution lower or equal zero!");
     	bCovpos[0]=0; bCovpos[2]=0;
    	}
	Float_t dcaToVertexXYpos = bpos[0];
    	Float_t dcaToVertexZpos = bpos[1];
    	fRpr->Fill(dcaToVertexZpos);
	fdcaToVertexXY->Fill(dcaToVertexXYpos);

	//if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5))
	if((TMath::Abs(dcaToVertexZpos)>2.5))
        continue;   //    

//                    if (!fMaxDCAtoVtxCut->AcceptTrack(track)) continue;

    	fdcaToVertexXYafterCut->Fill(dcaToVertexXYpos);
//

// **********************kink analysis*****************************
	
	//  loop on kinks
                if(indexKinkPos<0){     ////mother kink
	fptAllKink->Fill(track->Pt());

	fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);

	 AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);

	 // DCA kink
        Double_t  Dist2 = kink->GetDistance();
	Int_t eSDfLabel1=kink->GetLabel(0);
        TParticle *particle1= stack->Particle(TMath::Abs(eSDfLabel1));
        Int_t code1= particle1->GetPdgCode();
	Int_t eSDfLabeld=kink->GetLabel(1);
          TParticle *particled= stack->Particle(TMath::Abs(eSDfLabeld));
          Int_t dcode1= particled->GetPdgCode();
          Int_t mcProcssDa= particled->GetUniqueID();

	//    loop on MC daugthres for 3Pi    24/9/2010
       	Int_t nESDKpi =0;
       	if(mcProcssDa==4) {
       	Int_t firstDa=particle1->GetFirstDaughter();
       	Int_t lastDa=particle1->GetLastDaughter();

       	if( (lastDa<=0) || (firstDa<=0)  ) continue;

        if ( lastDa > mcEvent->GetNumberOfTracks() ) continue;
        if (firstDa > mcEvent->GetNumberOfTracks() ) continue;
//loop on secondaries
        for (Int_t kk=firstDa;kk<=lastDa;kk++) {
        if ( kk > 0 )    {
        TParticle*    daughter3=stack->Particle(kk);   // 24/9   
        Float_t dcode3= daughter3->GetPdgCode();
       	if (( ( code1==321 )&& ( dcode3==211  ))|| (( code1 == -321 )&& ( dcode3==-211)))    nESDKpi++ ;
        }
        }
        }
	Double_t hVzdau=particled->Vz();

// TPC mother momentum 

         const TVector3 vposKink(kink->GetPosition());
 	fPosiKink ->Fill( vposKink[0], vposKink[1]  );
   	Double_t  dzKink=vpos[2]-vposKink[2];
	Double_t z_vertex_kink = vposKink[2];
        if (z_vertex_kink > 0.4 && z_vertex_kink < 0.65) continue;
	if (TMath::Abs(z_vertex_kink)>225.0) continue;

//   lifitime
        Double_t tanLamda = track->GetTgl();  // 25/6/2010

   	Double_t lifeKink= (TMath::Abs( dzKink ))*( TMath::Sqrt(1.+ tanLamda*tanLamda) ) / (TMath::Abs( tanLamda)) ;
//
       	const TVector3 motherMfromKink(kink->GetMotherP());
       	const TVector3 daughterMKink(kink->GetDaughterP());

        Float_t qT=kink->GetQt();
        Float_t motherPt=motherMfromKink.Pt();
// TPC mother momentun
	Double_t trMomTPC=track->GetTPCmomentum();
	fQtAll->Fill(qT) ;  //  Qt   distr

	fptKink->Fill(motherMfromKink.Pt()); /// pt from kink
           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

 	if(  (TMath::Abs(rapiditK )) > 0.5 ) continue;  
        if ( (track->Pt())<.200)continue;  
//              eta selection 
        if ( TMath::Abs(trackEta) > 0.8) continue;   

        fQtMothP->Fill( track->P(), qT);

        if ( qT> fLowQt )  fqT1  ->Fill(qT) ;  //  Qt   distr

	fEta->Fill(trackEta) ;  //   Eta distr 
        fqT2->Fill(qT);  //      
	fKinkKaonBackg->Fill(motherPt);


	
	
    	if(((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==321))) fFakeKPi->Fill(track->Pt());
    	if(((TMath::Abs(code1)==211)&&(TMath::Abs(dcode1)==211))) fFakepipi->Fill(track->Pt());
	if( (kinkAngle<2.)  ) continue;
// 	cout<<" decay angle =="<<maxDecAngKmu<<endl;

	fPtCut1->Fill(trackPt);
        fAngMotherPi->Fill( track->P(),kinkAngle);


	Float_t signPt= tpcSign*trackPt;
	if  ( (TMath::Abs(code1)==211)&&(TMath::Abs(dcode1)==13))
           fRadiusPtPion->Fill( kink->GetR(), track->Pt()); //
	if( ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) {
           fRadiusPtKaon->Fill( kink->GetR(), track->Pt()); //
              }
	        if((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp )&&(TMath::Abs(rapiditK)<fRapiK)&&(label<nPrimary)) {
//          if((kink->GetR()>120.)&&(kink->GetR()<210.)&&(TMath::Abs(rapiditK)<0.7)&&(label<nPrim)) {
    	if( ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))) fQtKMu->Fill(qT);
    	if     ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))   fQtKEl->Fill(qT);
    	if     ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))   fQtKPi->Fill(qT);
    	if  (( nESDKpi>1) &&    ( ((code1)==321)&&((dcode1)==211)) )   fQtK3PiP->Fill(qT);
    	if  (( nESDKpi>1) &&    ( ((code1)==-321)&&((dcode1)==-211)) )   fQtK3PiM->Fill(qT);
        if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) {
         if(qT>fLowQt )        fHistPtKPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
        if(qT>fLowQt  )    {
      	if(code1>0.)  fHiPtKPDGP->Fill(trackPt ); // 26/feb  // ALL KAONS (pdg) inside ESD  kink sample
      	if(code1<0.)  fHiPtKPDGN->Fill(      trackPt ); // 26/feb  // ALL KAONS (pdg) inside ESD  kink sample
                   }
            fHistEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
            frapidESDK->Fill(rapiditK) ;  //18/feb rapiddistr of PDG kink ESD  kaons
      	if( qT > fLowQt   )  fHistQt2->Fill(qT);  // PDG ESD kaons            
           fRadiusPt->Fill( kink->GetR(), track->Pt()); //
     }
     }

        Double_t maxDecAngKmu=f1->Eval(track->P() ,0.,0.,0.);
        Double_t maxDecAngpimu=f2->Eval(track->P(), 0.,0.,0.);
	if(TMath::Abs(code1)==321) fAngMomK->Fill(track->P(), kinkAngle); //copied from below
        if(TMath::Abs(code1)==211) fAngMomPi->Fill( track->P(), kinkAngle);



// invariant mass of mother track decaying to mu
         Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
        Float_t energyDaughterPi=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.139570*0.139570);
         Float_t energyDaughterKa=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.493677*0.493677); 

	Float_t p1XM= motherMfromKink.Px();
         Float_t p1YM= motherMfromKink.Py();
         Float_t p1ZM= motherMfromKink.Pz();
         Float_t p2XM= daughterMKink.Px();
         Float_t p2YM= daughterMKink.Py();
         Float_t p2ZM= daughterMKink.Pz();
         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());
	 Double_t invariantMassKpi= TMath::Sqrt((energyDaughterPi+p3Daughter)*(energyDaughterPi+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());
         Double_t invariantMassKK = TMath::Sqrt((energyDaughterKa+p3Daughter)*(energyDaughterKa+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());
	fQtInvM -> Fill ( invariantMassKmu,  qT);
           fInvMuNuAll->Fill(invariantMassKmu);
           fInvMassMuNuPtAll ->Fill(invariantMassKmu,  trackPt);
	//fRadiusPt->Fill( kink->GetR(), trackPt); 
	  if (qT> fLowQt )     fSignPtNcl->Fill( signPt  ,   tpcNCl   );   // copied from below

	if((qT> fLowQt )&&((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp ))&&(TMath::Abs(rapiditK )<fRapiK)) {
    // if((qT>0.120)&&((kink->GetR()>120.)&&(kink->GetR()<210.))&&(TMath::Abs(rapiditK )<0.7)) {
	fAngMotherKC->Fill(track->P(), kinkAngle);
	fkaonInvaiant->Fill(invariantMassKmu);
         fMinvPi->Fill(invariantMassKpi);
         fMinvKa->Fill(invariantMassKK);
         fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
                        }

//  tails cleaning
    	if(  ( tpcNCl<fLowCluster) ) continue; 
	// cleaning BG in tails
      	Int_t tpcNClHigh = -42.67+ (11./12.) *( kink->GetR() ) ;
	if (tpcNCl >  tpcNClHigh  )   fcodeH->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
        if (  tpcNCl  >tpcNClHigh)   fZkinkZDau->Fill( vposKink[2],hVzdau              );
           if (tpcNCl > tpcNClHigh) fRadiusPtFake->Fill( kink->GetR(), track->Pt()); //

        if ( tpcNCl > tpcNClHigh) continue;

      	Int_t tpcNClMin  = -78.84 + 0.829  *( kink->GetR() ) ;
        if ( tpcNCl < tpcNClMin ) continue;

	fPtKPDG->Fill(track->Pt());    // ALL  K-candidates until now 


	if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<fRapiK  )&&(invariantMassKmu<0.8)){
        fAngMotherKKinks->Fill(track->P(), kinkAngle);
        fPtCut2->Fill(trackPt );

	//  maximum angles selection with some error cut 
        if( (kinkAngle<maxDecAngpimu*1.2) ) continue;
        if ( (kinkAngle>maxDecAngKmu*0.98) && ( track->P() >1.2 )) continue;  
	fPtCut3->Fill(trackPt);

	//  here the kaons selected by the decay features
        fTPCSignlMotherK->Fill( track->GetInnerParam()->GetP() ,(track->GetTPCsignal() )) ;
//nsigma cut
	if ( nsigma > 3.5) continue;
	fPtKaon->Fill(track->Pt());   //
        if(code1>0. ) fPtKaonP->Fill(track->Pt()) ;   //
        if ( code1 <0.)fPtKaonN->Fill(track->Pt()) ;   //
	fTPCSignalP->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()  ) ) ;
	fTPCMomNSgnl->Fill(track->GetInnerParam()->GetP() ,nsigmall );	
//  loop on kink daughters inside  mother's loop 
        Int_t ESDLabelM   =  0. ;
        Int_t ESDLabelD   =  0. ;
        Double_t dEdxDauMC  =   0.0;
 //       Double_t raDAU=0.;
        Int_t Ikink =0;
        Int_t IRkink =0;
	for (Int_t jTrack = 0; jTrack < esd->GetNumberOfTracks(); jTrack++) {

    	AliESDtrack* trackDau = esd->GetTrack(jTrack);
    	if (!trackDau) {
      	Printf("ERROR: Could not receive track %d", jTrack);
      	continue;
   	}
        Int_t indexKinkDAU =trackDau->GetKinkIndex(0);
        if (indexKinkDAU <0  ){
        AliESDkink *kinkDau=esd->GetKink(TMath::Abs(indexKinkDAU)-1);
   //            raDAU= kinkDau->GetR();
        ESDLabelM=kinkDau->GetLabel(0);   //  mothers's label
        ESDLabelM = TMath::Abs(ESDLabelM);
        ESDLabelD=kinkDau->GetLabel(1);   //  Daughter'slabel
        ESDLabelD = TMath::Abs(ESDLabelD);
        if ( kink->GetR() == kinkDau->GetR() ) IRkink++;
        if ( ESDLabelM == label ) Ikink++  ;
   	}
  //           if (( ESDLabelM== eSDfLabel1))   { 
        if (   (Ikink >0)  && (IRkink>0 )       )   {
 	if((indexKinkDAU      >0))    dEdxDauMC   = trackDau->GetTPCsignal()     ;  //  daughter kink 
   	}
   	}

	fRadiusNclCln->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
        fRadiusPtcln->Fill( kink->GetR(), trackPt); // 
	fInvMassMuNuPt ->Fill(invariantMassKmu,  trackPt);
	fTPCSignlPtpc->Fill(trMomTPC  , (track->GetTPCsignal()  ) ) ;               //trMomTPC
	fMothKinkMomSignl ->Fill(dEdxKinkDau  , (track->GetTPCsignal()  ) ) ;
	fTPCSignlKinkDau->Fill( daughterMKink.Mag(), dEdxKinkDau  ) ;  //  daughter kink
	fTPCMomNSigmaAllKaon->Fill(trMomTPC ,nsigmall );	
	fnSigmaTPC->Fill(nsigmall );
	fradiurKink->Fill(kink->GetR());  // kink 
	fLenthKink->Fill(lifeKink);  // kink
	fEtaK->Fill(trackEta);
	frapiKESD ->Fill(rapiditK);  //  rapidityof kaons
	fzVertexPositionKinKvsKinkRad->Fill( vposKink[2], kink->GetR() );
	
	Float_t signPt= tpcSign*trackPt;
	//fSignPtNcl->Fill( signPt, tpcNCl); // add before
	fSignPtrapiK->Fill( signPt  , rapiditK  );
	frapiKNcl->Fill( rapiditK, tpcNCl);	
	fSignPt->Fill( signPt );
	fChi2NclTPC->Fill( (track->GetTPCchi2()),tpcNCl );
	fRatioChi2Ncl-> Fill (  (track->GetTPCchi2()/track->GetTPCclusters(0)  )) ;	
	flifetime->Fill(( lifeKink*.493667   )  /track->P()   ) ;
	fPtKinkKaon->Fill(track->Pt());		
	//fDCAkink->Fill( Dist2);   //added below
	fPtKink->Fill(track->Pt()); ///  pp 7tev bins 
	//  kaons from k to mun and k to pipi and to e decay 
         if(( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||
           ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||
           ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) {

        if ( label<=nPrimary ) fPtPrKink->Fill(track->Pt());

        flifTiESDK->Fill(  ( lifeKink    /track->P() )*0.493667);



          fKinkKaon->Fill(track->Pt());
          fDCAkink->Fill( Dist2   );
          fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[0], vposKink[2]  );
          fPosiKinKYZ->Fill( vposKink[1], vposKink[2]  );
        if( code1>0.)     fkinkKaonP->Fill(trackPt);  //                  kPtPID kink-kaon
        if( code1<0.)    fkinkKaonN->Fill(trackPt);    //    PID kink-kaon
//     daughters
        if ( TMath::Abs(dEdxDauMC - track->GetTPCsignal() ) <  2)  fcode2->Fill( TMath::Abs(code1), TMath::Abs(dcode1));
  	fTPCSgnlPtpc   ->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 
        fMothKinkMomSgnlD->Fill( dEdxDauMC     , (track->GetTPCsignal()  ) ) ;
                             }
        else {
        fKinkKaonBg->Fill(track->Pt());
        fMothKinkMomSgnl ->Fill( dEdxDauMC     , (track->GetTPCsignal()  ) ) ;
        if ( TMath::Abs(dEdxDauMC - track->GetTPCsignal() ) <  2)  fcodeDau2->Fill( TMath::Abs(code1), TMath::Abs(dcode1));
  	fTPCSgnlKinkDau->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 

        fMinvPr->Fill(invariantMassKmu);

        fDCAkinkBG->Fill( Dist2   );
        fPosiKinKBgXY->Fill( vposKink[0], vposKink[1]  );
        fPosiKinKBgZY->Fill( vposKink[2], vposKink[1]  );
      fPosiKinKBgZX->Fill( vposKink[2], kink->GetR() );  //  31/7/11 
        if( code1>0.)     fKinKBGP  ->Fill(   trackPt);   //all PID kink-kaon
        if( code1<0.)     fKinKBGN  ->Fill( trackPt);   //all PID kink-kaonl
 	fdcodeH->Fill( TMath::Abs(code1), TMath::Abs(dcode1));   // put it here,  22/10/2009
        if(code1<0.)      fcode4->Fill(TMath::Abs(code1), TMath::Abs(dcode1));

          }   // primary and all +BG 



	 if(tpcSign >0.) fPtKinkPos ->Fill( track->Pt()) ;   //K-plus bins Comb 
	if(tpcSign <0.)  fPtKinkNeg ->Fill( track->Pt()) ;   //K-minus bins Comb
	fPtKinkK0->Fill(track->Pt()); ///  K0    bins     
	if(tpcSign >0.)        fPtKinkK0P->Fill( track->Pt() ) ;   //K-plus bins K0 Peter 
        if(tpcSign <0.)        fPtKinkK0N ->Fill( track->Pt() ) ;   //K-minus bins K0 Peter 
	
	fPtKinkGyu->Fill(track->Pt()); ///   K-charged High Pt Gyula
        if(tpcSign >0.)        fPtKinkGyuP->Fill( track->Pt()) ;   //K-plus charged high pt ,Gyula
        if(tpcSign <0.)        fPtKinkGyuN ->Fill( track->Pt()) ;   //K-minus bins High Pt charged Gyula

	fKinKRbn->Fill(track->Pt());       // TOF 	
	fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK);

	
	} // kink cut ends

	} // mother kink lopp ends


	}// end of the track loop
//	} //0-5% centrality

  PostData(1, fOutputList);
}      
/*
//_________________________________________________________________________
Float_t AliAnalysisTaskKinkpp5TeVMC::GetVertex(AliESDEvent* aod) const

{

  Float_t zvtx = -999;

  const AliESDVertex* vtxAOD = aod->GetPrimaryVertex();
  if (!vtxAOD)
    return zvtx;
  if(vtxAOD->GetNContributors()>0)
    zvtx = vtxAOD->GetZ();

  return zvtx;
}

*/
//________________________________________________________________________
void AliAnalysisTaskKinkpp5TeVMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

}

//-------------primary vertex-------

const AliESDVertex* AliAnalysisTaskKinkpp5TeVMC::GetEventVertex(const AliESDEvent* esd) const
  // Get the vertex from the ESD and returns it if the vertex is valid

{
  // Get the vertex 
  	const AliESDVertex* vertex = esd->GetPrimaryVertexTracks();

  //if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
  	if((vertex->GetStatus()==kTRUE)) return vertex;
  	else
  	{
    	 vertex = esd->GetPrimaryVertexSPD();
      	if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>0)) return vertex;
  //   if((vertex->GetStatus()==kTRUE)) return vertex;
     	else
     	return 0;
  	}
}

