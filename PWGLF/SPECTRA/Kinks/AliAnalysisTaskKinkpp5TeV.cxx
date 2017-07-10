//========================================================//
//Date:11/05/2017

//Author :: Nur Hussain and Buddhadeb Bhattacharjee, Gauhati University 
// Thanks to Martha Spyropoulou-Stassinaki for her suggestions for the modification
//purpose::Charged kaon identification using "Kink topology" for pp-5.02 TeV
//========================================================//
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
#include "AliAnalysisTaskKinkpp5TeV.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
//#include "AliMultSelectionCuts.h"

ClassImp(AliAnalysisTaskKinkpp5TeV)

//-------------------------------------------------------------------------

AliAnalysisTaskKinkpp5TeV::AliAnalysisTaskKinkpp5TeV()
: AliAnalysisTaskSE(),fOutputList(0), fHistPt(0),fVtxCut(10.),fMultiplicity(0),fIncompletEvent(0),fMultpileup(0), fMultV0trigger(0),fZvertex(0),fEventVertex(0),
        fRatioCrossedRows(0),fZvXv(0), fZvYv(0), fXvYv(0),fRpr(0),fdcaToVertexXY(0),fdcaToVertexXYafterCut(0),fptAllKink(0),fRatioCrossedRowsKink(0),fPosiKink(0),
        fQtAll(0),fptKink(0),fQtMothP(0),fqT1(0),fEta(0),fqT2(0),fKinkKaonBackg(0),f1(0), f2(0),fPtCut1(0),fAngMotherPi(0),
        fQtInvM(0),fInvMuNuAll(0),fInvMassMuNuPtAll(0),fRadiusPt(0),fKinkRadUp(200.), fKinkRadLow(130.), fLowCluster(30), fLowQt(.12), fRapiK(0.5),
        fAngMotherKC(0),fkaonInvaiant(0),fRadiusNcl(0), fPtKPDG(0),fAngMotherKKinks(0),fPtCut2(0),fPtCut3(0),fTPCSignlMotherK(0),fPtKaon(0), fPtKaonP(0), fPtKaonN(0),
         fTPCSignalP(0),fRadiusNclCln(0),fRadiusPtcln(0),fInvMassMuNuPt(0),fTPCSignlPtpc(0),fMothKinkMomSignl(0),fTPCSignlKinkDau(0),fTPCMomNSigmaAllKaon(0),
        fnSigmaTPC(0),fradiurKink(0),fLenthKink(0),fEtaK(0),frapiKESD(0),fzVertexPositionKinKvsKinkRad(0),fSignPtNcl(0),fSignPtrapiK(0),frapiKNcl(0),fSignPt(0),
        fChi2NclTPC(0),fRatioChi2Ncl(0),flifetime(),fPtKinkKaon(0),fDCAkink(0),fPtKink(0),fPtKinkPos(0),fPtKinkNeg(0),fPtKinkK0(0),fPtKinkK0P(0),fPtKinkK0N(0),
        fPtKinkGyu(0),fPtKinkGyuP(0),fPtKinkGyuN(0),fKinKRbn(0),fradPtRpDt(0),fAngMomK(0),fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),fPIDResponse(0),fNumberOfEvent(0),
        fNumberOfEvent_cent(0),fESDtrackCuts(0),fbgCleaningHigh(0),fTPCSignalPt(0),fqTvsPt(0),fInvMassPt(0),fCent(0),fEventVsCentrality(0)

{}

//________________________________________________________________________
AliAnalysisTaskKinkpp5TeV::AliAnalysisTaskKinkpp5TeV(const char *name, Float_t lRadiusKUp,  Float_t lRadiusKLow, Int_t lNCluster, Float_t lLowQtValue, Float_t yRange) 
  : AliAnalysisTaskSE(name),  fOutputList(0), fHistPt(0),fVtxCut(10.),fMultiplicity(0),fIncompletEvent(0),fMultpileup(0), fMultV0trigger(0),fZvertex(0),fEventVertex(0),
	fRatioCrossedRows(0),fZvXv(0), fZvYv(0), fXvYv(0),fRpr(0),fdcaToVertexXY(0),fdcaToVertexXYafterCut(0),fptAllKink(0),fRatioCrossedRowsKink(0),fPosiKink(0),
	fQtAll(0),fptKink(0),fQtMothP(0),fqT1(0),fEta(0),fqT2(0),fKinkKaonBackg(0),f1(0), f2(0),fPtCut1(0),fAngMotherPi(0),
	fQtInvM(0),fInvMuNuAll(0),fInvMassMuNuPtAll(0),fRadiusPt(0),fKinkRadUp(200.), fKinkRadLow(130.), fLowCluster(30), fLowQt(.12), fRapiK(0.5),
	fAngMotherKC(0),fkaonInvaiant(0),fRadiusNcl(0), fPtKPDG(0),fAngMotherKKinks(0),fPtCut2(0),fPtCut3(0),fTPCSignlMotherK(0),fPtKaon(0), fPtKaonP(0), fPtKaonN(0),
	 fTPCSignalP(0),fRadiusNclCln(0),fRadiusPtcln(0),fInvMassMuNuPt(0),fTPCSignlPtpc(0),fMothKinkMomSignl(0),fTPCSignlKinkDau(0),fTPCMomNSigmaAllKaon(0),
	fnSigmaTPC(0),fradiurKink(0),fLenthKink(0),fEtaK(0),frapiKESD(0),fzVertexPositionKinKvsKinkRad(0),fSignPtNcl(0),fSignPtrapiK(0),frapiKNcl(0),fSignPt(0),
	fChi2NclTPC(0),fRatioChi2Ncl(0),flifetime(),fPtKinkKaon(0),fDCAkink(0),fPtKink(0),fPtKinkPos(0),fPtKinkNeg(0),fPtKinkK0(0),fPtKinkK0P(0),fPtKinkK0N(0),
	fPtKinkGyu(0),fPtKinkGyuP(0),fPtKinkGyuN(0),fKinKRbn(0),fradPtRpDt(0),fAngMomK(0),fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),fPIDResponse(0),fNumberOfEvent(0),
	fNumberOfEvent_cent(0),fESDtrackCuts(0),fbgCleaningHigh(0),fTPCSignalPt(0),fqTvsPt(0),fInvMassPt(0),fCent(0),fEventVsCentrality(0)

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
void AliAnalysisTaskKinkpp5TeV::UserCreateOutputObjects()
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
	Double_t gPtpp[57] = { //pp 5Tev bining
0.01,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65, 0.7, 0.75,
0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4,
2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0,5.5, 6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0
 };

  	Double_t gPt13K0PKal[45]=    { 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,//  9bins 
                                 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,       //12 bins
                                 2.2,2.4,2.6,2.8,3.0,                                    //5 bins
  	3.3,3.6,3.9,4.2,4.6,5.0,5.4,5.9, 6.5,7.0,7.5,8.0,8.5,9.2,10.,11.,12.,13.5,15.} ;  //  19
	
	// from Gyula , 28/11/2015
 	//  const Int_t nPtBins = 68;
        //Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.45,

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
	fkaonInvaiant= new TH1F("fkaonInvaiant","Invar m(kaon) from kink->mu+netrino decay",600,0.10, 0.7); //  23/8/2013	
	fRadiusNcl = new TH2F("fRadiusNcl","kink radius vrs Nclust,K",75,100.,250., 80,0, 160);
	fPtKPDG = new TH1F("fPtKPDG", "P_{T}Kaon distribution",56, gPtpp);
	fAngMotherKKinks = new TH2F("fAngMotherKKinks","Decay angle vrs Mother Mom,Kinks",300,0.0,15.0,100,0.,100.);
	fPtCut2 = new TH1F("fPtCut2", "P_{T}Kaon distribution",56, gPtpp);
	fPtCut3 = new TH1F("fPtCut3", "P_{T}Kaon distribution",56, gPtpp);
	fTPCSignlMotherK = new TH2F("fTPCSignlMotherK","TPC signal de/dx Mom,K",1000,0.0,20.,150, 0.,300.);
	fPtKaon = new TH1F("fPtKaon", "P_{T} distribution of K all charged",56, gPtpp);	
	fPtKaonP = new TH1F("fPtKaonP", "P_{T} of distribution of K^{+}",56, gPtpp);	
	fPtKaonN = new TH1F("fPtKaonN", "P_{T} distribution of K^{-}",56, gPtpp);	
	 fTPCSignalP = new TH2F("fTPCSignalP","TPC signal de/dx Mom,K",1000,0.0,20.0,150,0.,300.);
	fRadiusNclCln = new TH2F("fRadiusNclCln","kink radius vrs Nclust,K Clean ",75,100.,250., 80,0, 160);
	fRadiusPtcln =new TH2F("fRadiusPtcln","radius vs pt clean ",80, 90.,250.,100, 0.,10.);
	fInvMassMuNuPt =new TH2F("fInvMassMuNuPt","Invariant Kaon = mass-munu  vs pt  ",600, 0.10, 0.7, 100, 0.0, 10.0  );
	fTPCSignlPtpc = new TH2F("fTPCSignlPtpc"," TPC signal de/dx vs Mom TPC,K after all the cut; TPC momentum (P) (GeV/c); dE/dx (a.u.) ",300,0.0,15.0,100, 0., 250.    );
	fMothKinkMomSignl = new TH2F("fMothKinkMomSignl","TPC signal de/dx Mom TPC of Kink; Daughter momentum (GeV/c); dE/dx (a.u.)  ",100,0.0,250.0,100, 0., 250.);
	fTPCSignlKinkDau = new TH2F("fTPCSignlKinkDau","TPC signal daughter kink de/dx Mom,K",500,0.0,10.0,400,0.,250.);
	fTPCMomNSigmaAllKaon = new TH2F("fTPCMomNSigmaAllKaon","TPC signal de/dx Mom TPC,all K  ",300,0.0,15.0,40 , -10., 10.);
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
	fPtKinkKaon = new TH1F ("fPtKinkKaon","p_{T} Kaon kinks identied",56, gPtpp);
	fDCAkink = new TH1F("fDCAkink ", "DCA kink vetrex ",50, 0.0,1.0);
	fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  distribution, counts",56, gPtpp);
	fPtKinkPos= new TH1F("fPtKinkPos", "Pos P_{T} K^{+} Kink distribution, counts",56, gPtpp );
	fPtKinkNeg= new TH1F("fPtKinkNeg", "Pos P_{T} K^{-} Kink   distribution, counts",56, gPtpp );
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
	fNumberOfEvent_cent = new TH1F("fNumberOfEvent_cent", "the number of events in this run", 10, 0., 6.);
	fbgCleaningHigh = new TH1F ("fbgCleaningHigh"," BG cleaning histo 1", 300, 0, 300);
	fTPCSignalPt = new TH2F("fTPCSignalPt"," TPC signal de/dx vs Mom pt,K after all the cut; p_{T}) (GeV/c); dE/dx (a.u.) ",300,0.0,15.0,100, 0., 250.    );
        fqTvsPt=new TH2F("fQtMothP", " Qt vrs Mother Pt", 300, 0., 15.0,100, 0.,0.400);
        fInvMassPt =new TH2F("fInvMassPt","Invariant Kaon = mass-munu  vs pt  ",600, 0.10, 0.7, 300, 0.0, 15.0  );
	fCent = new TH1F ("fCent","the centrality distribution",100, 0., 100.);
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
	 fOutputList->Add(fNumberOfEvent_cent);
	 fOutputList->Add(f1);
	 fOutputList->Add(fbgCleaningHigh);
	 fOutputList->Add(f2);
	 fOutputList->Add(fTPCSignalPt);
	 fOutputList->Add(fqTvsPt);
	 fOutputList->Add(fInvMassPt);
	 fOutputList->Add(fCent);
	 fOutputList->Add(fEventVsCentrality);
	

PostData(1, fOutputList);
}


//___________________user exec_____________________________________________________
void AliAnalysisTaskKinkpp5TeV::UserExec(Option_t *) 
{
  //initialized values
	Float_t dca1[2], cov1[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z for the vertex cut
	
	Double_t nsigmaPion =-100.0;
  //     Double_t nsigmaDau  =-100.0;
       	Double_t dEdxKinkDau =0.0;
	Double_t nsigmall = 100.0;
       	Double_t nsigma = 100.0;
	Double_t vtrack[3], ptrack[3];

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
// Number ESD tracks 

	Int_t nESDTracks =  esd->GetNumberOfTracks();
        //fMultiplicity->Fill(nESDTracks);
        fNumberOfEvent->Fill(1.5);
        fMultiplicity->Fill(2);
        //fESDtrackCuts->SetMaxDCAToVertexXYPtDep(0.0105 + 0.0350/ 
        fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0105 + 0.0350/pt^1.01");
//       fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
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
        fCent->Fill(cent);
	fNumberOfEvent->Fill(4.5);

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
 	}
//if((indexKinkDau >0)&& (nsigmaPion>1.2)) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
// Ayto mexri 26/11/2012     if(indexKinkDau >0) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
   	}



//	  if (cent>=0 && cent<5) {
	fNumberOfEvent_cent->Fill(3);

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


    	UInt_t status=track->GetStatus();

    	if((status&AliESDtrack::kITSrefit)==0) continue;
    	if((status&AliESDtrack::kTPCrefit)==0) continue;
      	if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;
//	fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0105 + 0.0350/(pow((track->Pt()),1.1))");
//	fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0105 + 0.0350/pt^1.01");
  //  	fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

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

//	if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5))
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

// TPC mother momentum 

         const TVector3 vposKink(kink->GetPosition());
 	fPosiKink ->Fill( vposKink[0], vposKink[1]  );
   	Double_t  dzKink=vpos[2]-vposKink[2];
	Double_t z_vertex_kink = vposKink[2];
	//cout <<" zvertex kink"<< z_vertex_kink<<endl;
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
      // mss2015  if ( TMath::Abs(trackEta) > 0.8) continue;  // new  NOv   2014 
        if ( TMath::Abs(trackEta) > 0.8) continue;   

        fQtMothP->Fill( track->P(), qT);

//        if ( qT> 0.12 )  fqT1  ->Fill(qT) ;  //  Qt   distr

	fEta->Fill(trackEta) ;  //   Eta distr 
        fqT2->Fill(qT);  //      
	fKinkKaonBackg->Fill(motherPt);


	//          maximum decay angle at a given mother momentum
           //Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
      	Double_t maxDecAngKmu=f1->Eval(track->P() ,0.,0.,0.);
       	Double_t maxDecAngpimu=f2->Eval(track->P(), 0.,0.,0.);
//cout << "kink angle==="<<kinkAngle<<endl;
	if( (kinkAngle<2.)  ) continue;
// 	cout<<" decay angle =="<<maxDecAngKmu<<endl;

	fPtCut1->Fill(trackPt );
        fAngMotherPi->Fill( track->P(),kinkAngle);

// invariant mass of mother track decaying to mu
         Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
         Float_t p1XM= motherMfromKink.Px();
         Float_t p1YM= motherMfromKink.Py();
         Float_t p1ZM= motherMfromKink.Pz();
         Float_t p2XM= daughterMKink.Px();
         Float_t p2YM= daughterMKink.Py();
         Float_t p2ZM= daughterMKink.Pz();
         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());

	fQtInvM -> Fill ( invariantMassKmu,  qT);
        fInvMuNuAll->Fill(invariantMassKmu);
        fInvMassMuNuPtAll ->Fill(invariantMassKmu,  trackPt);
	fRadiusPt->Fill( kink->GetR(), trackPt); 

	 if( ( kink->GetR()> fKinkRadLow ) && ( kink->GetR() <fKinkRadUp   )  )  {
    //  for systematics   if( ( kink->GetR()> 130 ) && ( kink->GetR() < 200 )  )  {
      	if (qT>fLowQt )  fAngMotherKC->Fill(track->P(), kinkAngle);
  //`    if (qT> 0.120  )  fAngMomKC->Fill(track->P(), kinkAngle); 
      // `if (qT> 0.04  )  fAngMomKC->Fill(track->P(), kinkAngle); 
        if ( qT> fLowQt ) fkaonInvaiant->Fill(invariantMassKmu);
          //   if ( qT>  0.12  ) fM1kaon->Fill(invariantMassKmu);
        if ( qT > fLowQt)
        //     if ( qT > 0.12  ) 
         fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
  	}
//  tails cleaning
    	if(  ( tpcNCl<fLowCluster) ) continue; 
	// cleaning BG in tails
  //    	Int_t tpcNClHigh = -31.67+ (11./12.)  *( kink->GetR() ) ;
	Int_t tpcNClHigh = -42.67+ (11./12.)  *( kink->GetR() ) ;
	fbgCleaningHigh->Fill(tpcNClHigh);
        if ( tpcNCl > tpcNClHigh) continue;
	Int_t tpcNClMin  = -78.84 + 0.829  *( kink->GetR() ) ;
//      	Int_t tpcNClMin  = -85.5 + (65./95.)  *( kink->GetR() ) ;
        if ( tpcNCl < tpcNClMin ) continue;

	fPtKPDG->Fill(track->Pt());    // ALL  K-candidates until now 


	if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<fRapiK  )&&(invariantMassKmu<0.8)){
     //          if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.8)){
  // systematics   if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>=130.)&&(kink->GetR()<=200.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.8)){
//
        fAngMotherKKinks->Fill(track->P(), kinkAngle);
        fPtCut2->Fill(trackPt );

	//  maximum angles selection with some error cut 
        if( (kinkAngle<maxDecAngpimu*1.2) ) continue;
       	if ( (kinkAngle>maxDecAngKmu*0.98) && ( track->P() >1.2 )) continue;  
	fPtCut3->Fill(trackPt );

	//  here the kaons selected by the decay features
        fTPCSignlMotherK->Fill( track->GetInnerParam()->GetP() ,(track->GetTPCsignal() )) ;
//nsigma cut
	if ( nsigma > 3.5) continue;
	fqT1  ->Fill(qT) ;
	fqTvsPt->Fill(trackPt, qT);
        fInvMassPt->Fill(invariantMassKmu, trackPt);
	fTPCSignalPt->Fill(trackPt, (track->GetTPCsignal()  ) ) ;
	fPtKaon->Fill(track->Pt());   //
        if(tpcSign >0.) fPtKaonP->Fill(track->Pt()) ;   //
        if ( tpcSign <0.)fPtKaonN->Fill(track->Pt()) ;   //
	fTPCSignalP->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()  ) ) ;
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
	fSignPtNcl->Fill( signPt, tpcNCl); 
	fSignPtrapiK->Fill( signPt  , rapiditK  );
	frapiKNcl->Fill( rapiditK, tpcNCl);	
	fSignPt->Fill( signPt );
	fChi2NclTPC->Fill( (track->GetTPCchi2()),tpcNCl );
	fRatioChi2Ncl-> Fill (  (track->GetTPCchi2()/track->GetTPCclusters(0)  )) ;	
	flifetime->Fill(( lifeKink*.493667   )  /track->P()   ) ;
	fPtKinkKaon->Fill(track->Pt());		
	fDCAkink->Fill( Dist2);
	fPtKink->Fill(track->Pt()); ///  pp 7tev bins 
	 if(tpcSign >0.) fPtKinkPos ->Fill( track->Pt()) ;   //K-plus bins Comb 
	if(tpcSign <0.)  fPtKinkNeg ->Fill( track->Pt()) ;   //K-minus bins Comb
	fPtKinkK0->Fill(track->Pt()); ///  K0    bins     
	if(tpcSign >0.)        fPtKinkK0P->Fill( track->Pt()         ) ;   //K-plus bins K0 Peter 
        if(tpcSign <0.)        fPtKinkK0N ->Fill( track->Pt()         ) ;   //K-minus bins K0 Peter 
	
	fPtKinkGyu->Fill(track->Pt()); ///   K-charged High Pt Gyula
        if(tpcSign >0.)        fPtKinkGyuP->Fill( track->Pt()    ) ;   //K-plus charged high pt ,Gyula
        if(tpcSign <0.)        fPtKinkGyuN ->Fill( track->Pt()  ) ;   //K-minus bins High Pt charged Gyula

	fKinKRbn->Fill(track->Pt());       // TOF 	
	fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK);
	fAngMomK->Fill(track->P(),kinkAngle);
	fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[2], vposKink[0]  );
          fPosiKinKYZ->Fill( vposKink[2], vposKink[1]  );

	
	} // kink cut ends

	} // mother kink lopp ends


	}// end of the track loop
//	} // 0-5% centrality cut

  PostData(1, fOutputList);
}      
/*
//_________________________________________________________________________
Float_t AliAnalysisTaskKinkPbPb::GetVertex(AliESDEvent* aod) const

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
void AliAnalysisTaskKinkpp5TeV::Terminate(Option_t *) 
{
  // Draw result to the screen
}


//-------------primary vertex-------

const AliESDVertex* AliAnalysisTaskKinkpp5TeV::GetEventVertex(const AliESDEvent* esd) const
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





