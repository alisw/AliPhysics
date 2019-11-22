/**************************************************************************
 *  Authors: Martha Spyropoulou-Stassinaki and the  members               * 
 * of the Greek group at Physics Department of Athens University          *
 * Paraskevi Ganoti, Anastasia Belogianni and Filimon Roukoutakis.        *
 * The method is applied in pp and Pb-Pb real data.                       *
 *                                                                        *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisKinkTaskMult13 class
//       Example of an analysis task for kink topology study
//      Kaons from kink topology are 'identified' in this code
//     Nov 2014 : Nominal R->120-210 cm,  Rapidity kaon 0.5, eta< 0.8
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TParticle.h"
#include <TVector3.h>
#include "TF1.h"
#include "AliAnalysisTaskSE.h"
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
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliEventCuts.h"
//#include "AliAnalysisUtils.h"
//#include "AliESDUtils.h"
//#include "AliPPVsMultUtils.h"
//#include "AliOADBContainer.h"
//#include "AliOADBMultSelection.h"
//#include "AliMultEstimator.h"
//#include "AliMultVariable.h"
//#include "AliMultInput.h"
#include "AliMultSelection.h"
/////////#include "AliTENDERSupplies.h"
//#include "AliMultSelectionCuts.h"
#include "AliAnalysisKinkTaskPbPb2018.h"
//#include "AliAnalysisUtils.h"
ClassImp(AliAnalysisKinkTaskPbPb2018)


//________________________________________________________________________
AliAnalysisKinkTaskPbPb2018::AliAnalysisKinkTaskPbPb2018(const char *name) 
  : AliAnalysisTaskSE(name)
  , fEventCuts(kFALSE)
  , fMultiplicityBeforeCuts(0)
  , fIncompletEv(0)
  , fMultiplicityAfterTriggerBit(0)
  , fMultiplicityAfterPileup(0)
  , fZMainVx(0)
  , fMultiplicityAfterVertexCut(0)
  , fTrackPtAfterTrackCuts(0)
  , fTrackPtAll(0)
  , fHistQtAll(0)
  , fQtAfterAcceptance(0)
  , fPtSelectedKaons(0)
  , fKaonKinkPtAfterKinkNclCut(0)
  , fEtaAfterAcceptance(0)
  , fTrackEtaSelectedKaons(0)
  , fPtAllKinks(0)
  , fgenpt(0)
  , fRadiusSelectedKinks(0)
  , fKinKRbn(0)
  , fKinkMomFromMother(0)
  , fInvMassKaonInR(0)
  , fPtFromMotherAllKinks(0)
  , fAngMomK(0)
  , fAngleMomKaonsinTPC(0)
  , fAngleVsMomentumKaonsInR(0)
  , fNclinTPCVsSignedPtSelectedKaons(0)
  , fRapidityVsSignedPtSelectedKons(0)
  , fNclVsRapiditySelectedKAons(0)
  , fSignedPtSelectedKaons(0)
  , fNclVsChi2SelectedKaons(0)
  , fChi2perTPCclusterSelectedKaons(0)
  , fRadiusVsNclInR(0)
  , fTPCSignalVsMomSelectedKaons(0)
  , fTPCSgnlPa(0)
  , fDCAz(0)
  , fDCAxy(0)
  , fnSigmToVx(0)
  , fKinkMothDau(0)
  , fZvXv(0)
  , fZvYv(0)
  , fXvYv(0)
  , fRapiditySelectedKaons(0)
  , fLifetimeSelectedKaons()
  , fKaonLifetimeSelectedKaons(0)
  //, fradPtRpDt(0)
  , fInvMassMuNuKaonTPC(0)
  , fQtInvMassKaonTPC(0)
  , fDCAkinkSelectedKaons(0)
  , fxyKinkPosition(0)
  , fPosiKinkK(0)
  , fPosiKinKXZ(0)
  , fPosiKinKYZ(0)
  , fZKinkProductionVsKinkRadSelectedKaons(0)
  , fQtVsKinkMomAfterAcceptance(0)
  , fdedxMthVsTPCMomSelectedKaons(0)
  , fNsigmaVsTPCmomSelectedKaons(0)
  , fSignalMthVsSignalDaughterSelectedKaons(0)
  , fNsigmaSelectedKaons(0)
  , fSignalDaughterVsDaughterMomSelectedKaons(0)
  , fPtPositiveSelectedKaons(0)
  , fPtNegativeSelectedKaons(0)
  , fNclustersVsRadiusSelectedKaons(0)
  , fRatioCrossedRows(0)
  , fRatioCrossedRowsKink(0)
  , fRadiusVsPtKaonTPC(0)
  , fPtVsRadiusSelectedKaons(0)
  , fPtVsInvMassSelectedKaons(0)
  , fInvMassPtKaonTPC(0)
  , fPtKaonInTPC(0)
  , fPtPreSelectedkinks(0)
  , fPtKinkBeforedEdx(0)
  , fAngleVsMomPreSelectedKinks(0)
  , fPtKinkK0(0)
  , fPtKinkK0P(0)
  , fPtKinkK0N(0)
  , fPtKinkGyu(0)
  , fPtKinkGyuP(0)
  , fPtKinkGyuN(0)
  , f1(0)
  , f2(0)
  , fListOfHistos(0)
  , fKinkRadUp(200.)
  , fKinkRadLow(130.)
  , fLowCluster(30)
  , fLowQt(.12)
  , fRapiK(0.5)
  , fPIDResponse()
  , fMultiBin1ChargedMulti(0)
  , fMultiBin2ChargedMulti(0)
  , fMultiBin3ChargedMulti(0)

  , fMultiBin1KaonKinksPos(0)
  , fMultiBin2KaonKinksPos(0)
  , fMultiBin3KaonKinksPos(0)

  , fMultiBin1KaonKinksNeg(0)
  , fMultiBin2KaonKinksNeg(0)
  , fMultiBin3KaonKinksNeg(0)

  , fMultiBin1Vertex(0)
  , fMultiBin2Vertex(0)
  , fMultiBin3Vertex(0)
 
  , fVertexNet(0)
  , fhPS(0)
  , fhvtxP(0)
  , fESDtrackCuts(0)
  , fpercentile(0)
//  , fUtils(0)
{
  // Constructor

  DefineInput(0, TChain::Class());
  fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts","AliESDtrackCuts");
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
void AliAnalysisKinkTaskPbPb2018::UserCreateOutputObjects() 
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

   Double_t gPt7Comb[61] = {
     0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0, 1.1, 1.2,
     1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
     3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.5, 6, 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10.0, 10.5, 11., 12., 13., 14.
   };  // 25/11/2013 from Francesco

   Double_t gPt7TOF[47] = { 0.2,0.25, 0.3,0.35,  0.4,0.45,  0.5,0.55,  0.6,0.65,  0.7,0.75,  0.8, 0.85, 0.9, 0.95, 1.0,
			    1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
			    2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0, 
			    3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0 };  //  Barbara TOF Kch
   // from Gyula , 28/11/2015
   //  const Int_t nPtBins = 68;
   //Double_t xBins[nPtBins+1] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.45,
   ///*
   Double_t gPt13HPtGyu[69] = {0. ,  0.05, 0.1,  0.15, 0.2,  0.25, 0.3,  0.35, 0.4, 0.45,
			       0.5,  0.55, 0.6,  0.65, 0.7,  0.75, 0.8,  0.85, 0.9,  0.95,
			       1.0,  1.1 , 1.2,  1.3 , 1.4,  1.5 , 1.6,  1.7 , 1.8,  1.9 ,
			       2.0,  2.2 , 2.4,  2.6 , 2.8,  3.0 , 3.2,  3.4 , 3.6,  3.8 ,
			       4.0,  4.5 , 5.0,  5.5 , 6.0,  6.5 , 7.0,  8.0 , 9.0,  10.0,
			       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 18.0, 20.0, 22.0, 24.0,
			       26.0, 28.0, 30.0, 32.0, 34.0, 36.0, 40.0, 45.0, 50.0 };
   //*/                                 
   // K0 binning from Peter Kalinak,  22/1/2016,   einai to idio me to gPt7K0 tou David
   Double_t gPt13K0PKal[45]=    { 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,//  9bins 
				  0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,       //12 bins
				  2.2,2.4,2.6,2.8,3.0,                                    //5 bins
				  3.3,3.6,3.9,4.2,4.6,5.0,5.4,5.9, 6.5,7.0,7.5,8.0,8.5,9.2,10.,11.,12.,13.5,15.} ;  //  19
   
   //
   fMultiplicityBeforeCuts= new TH1F("fMultiplicityBeforeCuts", "charge multiplicity ESD",300, 0.0,300.0); 
   fIncompletEv=new TH1F("fIncompletEv", "charge multiplicity ESD after Incom. cut",100, 0.0,300.0); 
   fMultiplicityAfterTriggerBit=new TH1F("fMultiplicityAfterTriggerBit", "charge multiplity  after Trigger",300, 0.0, 300.0); 
   fMultiplicityAfterPileup=new TH1F("fMultiplicityAfterPileup", "charge multipliESD after Pileup",300, 0.0, 300.0); 
   fZMainVx= new TH1D("fZMainVx", "ESD charge mult. Main Vertex", 60,-15.,15.); 
   fMultiplicityAfterVertexCut= new TH1F("fMultiplicityAfterVertexCut", "ESD charge mult. inside Main Vertx",300, 0.0, 300.0); 
   //
   fTrackPtAfterTrackCuts = new TH1F("fTrackPtAfterTrackCuts", "P_{T} distribution",1500, 0.0,15.0);
   fTrackPtAfterTrackCuts->GetXaxis()->SetTitle("P_{T} (GeV/c)");
   fTrackPtAfterTrackCuts->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
   fTrackPtAfterTrackCuts->SetMarkerStyle(kFullCircle);
   fTrackPtAll = new TH1F("fTrackPtAll", "P_{T} distribution",1500, 0.0,15.0); 
   fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300); 
   fQtAfterAcceptance= new TH1F("fQtAfterAcceptance", "Q_{T} distribution",100, 0.0,.300); 
   fPtSelectedKaons = new TH1F("fPtSelectedKaons", "P_{T}Kaon distribution",1500, 0.0,15.0); 
   fKaonKinkPtAfterKinkNclCut = new TH1F("fKaonKinkPtAfterKinkNclCut", "P_{T}Kaon distribution",1500, 0.0,15.0); 
   fEtaAfterAcceptance= new TH1F("fEtaAfterAcceptance", "Eta distribution", 26,-1.3, 1.3); 
   fTrackEtaSelectedKaons= new TH1F("fTrackEtaSelectedKaons", "EtaK distribution", 26,-1.3, 1.3); 
   fPtAllKinks= new TH1F("fPtAllKinks", "P_{T} all kinks",1500, 0.0,15.0); 
   fgenpt= new TH1F("fgenpt", "genpt   K distribution",300, 0.0,15.0); 
   //fRadiusSelectedKinks= new TH1F("fRadiusSelectedKinks", "radius  K generated",100, 50., 250.0);
   fRadiusSelectedKinks= new TH1F("fRadiusSelectedKinks", "radius  K ",100, 0.,1000.0);
   fKinKRbn= new TH1F("fKinKRbn", "p_{t}Kaon kinks identi[GeV/c],Entries",46,gPt7TOF); 
   fKinkMomFromMother= new TH1F("fKinkMomFromMother", "P_{T}Kaon kinks backgr",1500, 0.0,15.0); 
   //fInvMassKaonInR= new TH1F("fInvMassKaonInR","Invar m(kaon) from kink->mu+netrino decay",180,0.10, 1.0); 
   fInvMassKaonInR= new TH1F("fInvMassKaonInR","Invar m(kaon) from kink->mu+netrino decay",600,0.10, 0.7); //  23/8/2013
   //fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  distribution, counts",44, gPt7K0); 
   fPtFromMotherAllKinks= new TH1F("fPtFromMotherAllKinks", "P_{T}",1500, 0.0,15.0); 
   fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
   fAngleMomKaonsinTPC= new TH2F("fAngleMomKaonsinTPC","Decay angle vrs Mother Mom,Pi",1500, 0.0, 15.0 ,80, 0., 80.);
   fAngleVsMomentumKaonsInR= new TH2F("fAngleVsMomentumKaonsInR","Decay angle vrs Mother Mom,K",1500, 0.0, 15.0, 80, 0., 80.);
   fNclinTPCVsSignedPtSelectedKaons= new TH2F("fNclinTPCVsSignedPtSelectedKaons","SignPt vrs Ncl,K",80,-4.,4.0,70,20.,160.);
   fRapidityVsSignedPtSelectedKons= new TH2F("fRapidityVsSignedPtSelectedKons","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
   fNclVsRapiditySelectedKAons= new TH2F("fNclVsRapiditySelectedKAons","Eta vrs nclust,K",30,-1.5,1.5, 70,20, 160);
   fSignedPtSelectedKaons= new TH1F("fSignedPtSelectedKaons","SignPt ,K",80,-4.0,4.0);
   fNclVsChi2SelectedKaons= new TH2F("fNclVsChi2SelectedKaons","Chi2vrs nclust,K",100,0.,500., 70,20, 160);
   fChi2perTPCclusterSelectedKaons= new TH1F("fChi2perTPCclusterSelectedKaons","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
   fRadiusVsNclInR = new TH2F("fRadiusVsNclInR","kink radius vrs Nclust,K",75,100.,250., 80,0, 160);
   fTPCSignalVsMomSelectedKaons = new TH2F("fTPCSignalVsMomSelectedKaons","TPC signal de/dx Mom,K",1000,0.0,20.0,150,0.,300.);
   fTPCSgnlPa= new TH2F("fTPCSgnlPa","TPC signal de/dx Mom,K",1000,0.0,20.,150, 0.,300.);
   fDCAz = new TH1D("fDCAz", "rad distribution  PID pr",100,-10.0, 10.0);
   fDCAxy = new TH1D("fDCAxy", "dca  distribution PID  ",20,-1.,1.);
   fnSigmToVx = new TH1D("fnSigmToVx", "dca  distribution PID  ",80,0.,8.);
   fKinkMothDau= new TH2F("fKinkMothDau","TPC kink Moth Daugh ,K",50,0.0,2.5,50, 0., 2.5);
   fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-1.5,1.5,60, -15., 15.0);
   fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-1.5,1.5, 60, -15., 15.);
   fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
   fRapiditySelectedKaons=new TH1F("fRapiditySelectedKaons","rapid Kdistribution", 26,-1.3, 1.3); 
   fLifetimeSelectedKaons= new TH1F("fLifetimeSelectedKaons", "ct study of K-kinks",100,0.,1000.); 
   fKaonLifetimeSelectedKaons= new TH1F("fKaonLifetimeSelectedKaons", "Length of   K generated",100,0.,1000.); 
   //fradPtRpDt=new TH3F("fradPtRpDt","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
   //fInvMassMuNuKaonTPC= new TH1F("fInvMassMuNuKaonTPC", " Inv Mass MuNu all kink",180,0.1,1.0); 
   fInvMassMuNuKaonTPC= new TH1F("fInvMassMuNuKaonTPC", " Inv Mass MuNu all kink",600,0.1,0.7); //  23/8/2013
   fQtInvMassKaonTPC= new TH2F("fQtInvMassKaonTPC", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300); 
   fDCAkinkSelectedKaons = new TH1F("fDCAkinkSelectedKaons ", "DCA kink vetrex ",50, 0.0,1.0);
   fxyKinkPosition= new TH2F("fxyKinkPosition", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
   fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fZKinkProductionVsKinkRadSelectedKaons= new TH2F("fZKinkProductionVsKinkRadSelectedKaons", "z vrx kink rad    ",100, -300.0,300.0,100,100., 300.);
   fQtVsKinkMomAfterAcceptance = new TH2F("fQtVsKinkMomAfterAcceptance", " Qt vrs Mother P", 1500, 0., 15.0, 100, 0., 0.300);
   fdedxMthVsTPCMomSelectedKaons = new TH2F("fdedxMthVsTPCMomSelectedKaons","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,100, 0., 250.    );
   fNsigmaVsTPCmomSelectedKaons = new TH2F("fNsigmaVsTPCmomSelectedKaons","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,20 , -10., 10.);
   fSignalMthVsSignalDaughterSelectedKaons  = new TH2F("fSignalMthVsSignalDaughterSelectedKaons","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.0,100, 0., 250.    );
   fNsigmaSelectedKaons = new TH1F("fNsigmaSelectedKaons","TPC Nsigma  de/dx  TPC,K  ", 30 , -7.5, 7.5);
   fSignalDaughterVsDaughterMomSelectedKaons = new TH2F("fSignalDaughterVsDaughterMomSelectedKaons","TPC signal de/dx Mom,K",500,0.0,10.0,100,0.,250.);
   //fPtPositiveSelectedKaons= new TH1F("fPtPositiveSelectedKaons", "Pos P_{T}Kaon Kink  distribution, counts",44, gPt7K0); 
   fPtPositiveSelectedKaons= new TH1F("fPtPositiveSelectedKaons", "Pos P_{T}Kaon Kink  distribution, counts", 1500, 0.0, 15.0); 
   fPtNegativeSelectedKaons= new TH1F("fPtNegativeSelectedKaons", "Neg P_{T}Kaon Kink  distribution, counts", 1500, 0.0, 15.0); 
   fNclustersVsRadiusSelectedKaons = new TH2F("fNclustersVsRadiusSelectedKaons","kink radius vrs Nclust,K Clean ",75,100.,250., 80,0, 160);
   fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);
   fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);
   fRadiusVsPtKaonTPC =new TH2F("fRadiusVsPtKaonTPC","radius vs pt  ",80, 90.,250.,1500, 0.,15.);
   fPtVsRadiusSelectedKaons =new TH2F("fPtVsRadiusSelectedKaons","radius vs pt clean ",80, 90.,250.,1500, 0.,15.);
   //fPtVsInvMassSelectedKaons =new TH2F("fPtVsInvMassSelectedKaons","Invariant mass-munu  vs pt  ",180, 0.10, 1.00, 100, 0.0, 10.0  );
   fPtVsInvMassSelectedKaons =new TH2F("fPtVsInvMassSelectedKaons","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 1500, 0.0, 15.0  );// 23/8/2013
   fInvMassPtKaonTPC =new TH2F("fInvMassPtKaonTPC","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 1500, 0.0, 15.0  );// 23/8/2013
   fPtKaonInTPC = new TH1F("fPtKaonInTPC", "P_{T}Kaon distribution",1500, 0.0,15.0); 
   fPtPreSelectedkinks = new TH1F("fPtPreSelectedkinks", "P_{T}Kaon distribution",1500, 0.0,15.0); 
   fPtKinkBeforedEdx = new TH1F("fPtKinkBeforedEdx", "P_{T}Kaon distribution",1500, 0.0,15.0); 
   fAngleVsMomPreSelectedKinks = new TH2F("fAngleVsMomPreSelectedKinks","Decay angle vrs Mother Mom,Kinks",1500 ,0.0, 15.0,100,0.,100.);
   fPtKinkK0= new TH1F("fPtKinkK0", "P_{T}Kaon Kink  distribution, counts",44, gPt13K0PKal); 
   fPtKinkK0P= new TH1F("fPtKinkK0P", "P_{T} KPl Kink  distribution, counts",44, gPt13K0PKal); 
   fPtKinkK0N= new TH1F("fPtKinkK0N", "P_{T} KMn Kink  distribution, counts",44, gPt13K0PKal); 
   fPtKinkGyu= new TH1F("fPtKinkGyu", "P_{T}Kaon Kink  distribution, counts",68, gPt13HPtGyu); 
   fPtKinkGyuP= new TH1F("fPtKinkGyuP", "P_{T} KPl Kink  distribution, counts",68, gPt13HPtGyu); 
   fPtKinkGyuN= new TH1F("fPtKinkGyuN", "P_{T} KMn Kink  distribution, counts",68, gPt13HPtGyu);

   fMultiBin1ChargedMulti= new TH1F("fMultiBin1ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin2ChargedMulti= new TH1F("fMultiBin2ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin3ChargedMulti= new TH1F("fMultiBin3ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);

   fMultiBin1KaonKinksPos = new TH1F("fMultiBin1KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonKinksPos = new TH1F("fMultiBin2KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonKinksPos = new TH1F("fMultiBin3KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);

   fMultiBin1KaonKinksNeg = new TH1F("fMultiBin1KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonKinksNeg = new TH1F("fMultiBin2KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonKinksNeg = new TH1F("fMultiBin3KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);

   fMultiBin1Vertex= new TH1D("fZMainVxfMultiBin1Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin2Vertex= new TH1D("fZMainVxfMultiBin2Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin3Vertex= new TH1D("fZMainVxfMultiBin3Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);

   fVertexNet= new TH1D("fVertexNet", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fhPS = new TH1F("fhPS", " ",10, 0.0,10.0);
   fhvtxP = new TH1F("fhvtxP", " ",10, 0.0,10.0);
   fpercentile =  new TH1F("fpercentile", " ", 100, 0.0, 100.0);
   
   fListOfHistos=new TList();
   
   fListOfHistos->Add(fMultiplicityBeforeCuts);
   fListOfHistos->Add(fIncompletEv);
   fListOfHistos->Add(fMultiplicityAfterTriggerBit);
   fListOfHistos->Add(fMultiplicityAfterPileup);
   fListOfHistos->Add(fZMainVx);
   fListOfHistos->Add(fMultiplicityAfterVertexCut);
//
   fListOfHistos->Add(fTrackPtAfterTrackCuts);
   fListOfHistos->Add(fTrackPtAll);
   fListOfHistos->Add(fHistQtAll);
   fListOfHistos->Add(fQtAfterAcceptance);
   fListOfHistos->Add(fPtSelectedKaons);
   fListOfHistos->Add(fKaonKinkPtAfterKinkNclCut);
   fListOfHistos->Add(fEtaAfterAcceptance);
   fListOfHistos->Add(fTrackEtaSelectedKaons);
   fListOfHistos->Add(fPtAllKinks);
   fListOfHistos->Add(fgenpt);
   fListOfHistos->Add(fRadiusSelectedKinks);
   fListOfHistos->Add(fKinKRbn);
   fListOfHistos->Add(fKinkMomFromMother);
   fListOfHistos->Add(fInvMassKaonInR);
   fListOfHistos->Add(fPtFromMotherAllKinks);
   fListOfHistos->Add(fAngMomK);
   fListOfHistos->Add(fAngleMomKaonsinTPC);
   fListOfHistos->Add(fAngleVsMomentumKaonsInR);
   fListOfHistos->Add(fNclinTPCVsSignedPtSelectedKaons);
   fListOfHistos->Add(fRapidityVsSignedPtSelectedKons);
   fListOfHistos->Add(fNclVsRapiditySelectedKAons);
   fListOfHistos->Add(fSignedPtSelectedKaons);
   fListOfHistos->Add(fNclVsChi2SelectedKaons);
   fListOfHistos->Add(fChi2perTPCclusterSelectedKaons);
   fListOfHistos->Add(fRadiusVsNclInR);
   fListOfHistos->Add(fTPCSignalVsMomSelectedKaons);
   fListOfHistos->Add(fTPCSgnlPa);
   fListOfHistos->Add(fDCAz);
   fListOfHistos->Add(fDCAxy);
   fListOfHistos->Add(fnSigmToVx);
   fListOfHistos->Add(fKinkMothDau);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fRapiditySelectedKaons);
   fListOfHistos->Add(fLifetimeSelectedKaons);
   fListOfHistos->Add(fKaonLifetimeSelectedKaons);
   //fListOfHistos->Add(fradPtRpDt);
   fListOfHistos->Add(fInvMassMuNuKaonTPC);
   fListOfHistos->Add(fQtInvMassKaonTPC);
   fListOfHistos->Add(fDCAkinkSelectedKaons);
   fListOfHistos->Add(fxyKinkPosition);
   fListOfHistos->Add(fPosiKinkK);
   fListOfHistos->Add(fPosiKinKXZ);
   fListOfHistos->Add(fPosiKinKYZ);
   fListOfHistos->Add(fZKinkProductionVsKinkRadSelectedKaons);
   fListOfHistos->Add(fQtVsKinkMomAfterAcceptance);
   fListOfHistos->Add(fdedxMthVsTPCMomSelectedKaons);
   fListOfHistos->Add(fNsigmaVsTPCmomSelectedKaons);
   fListOfHistos->Add(fSignalMthVsSignalDaughterSelectedKaons);
   fListOfHistos->Add(fNsigmaSelectedKaons);
   fListOfHistos->Add(fSignalDaughterVsDaughterMomSelectedKaons);
   fListOfHistos->Add(fPtPositiveSelectedKaons);
   fListOfHistos->Add(fPtNegativeSelectedKaons);
   fListOfHistos->Add(fNclustersVsRadiusSelectedKaons);
   fListOfHistos->Add(fRatioCrossedRows);
   fListOfHistos->Add(fRatioCrossedRowsKink);
   fListOfHistos->Add(fRadiusVsPtKaonTPC);
   fListOfHistos->Add(fPtVsRadiusSelectedKaons);
   fListOfHistos->Add(fPtVsInvMassSelectedKaons);
   fListOfHistos->Add(fInvMassPtKaonTPC);
   fListOfHistos->Add(fPtKaonInTPC);
   fListOfHistos->Add(fPtPreSelectedkinks);
   fListOfHistos->Add(fPtKinkBeforedEdx);
   fListOfHistos->Add(fAngleVsMomPreSelectedKinks);
   fListOfHistos->Add(fPtKinkK0);
   fListOfHistos->Add(fPtKinkK0P);
   fListOfHistos->Add(fPtKinkK0N);
   fListOfHistos->Add(fPtKinkGyu);
   fListOfHistos->Add(fPtKinkGyuP);
   fListOfHistos->Add(fPtKinkGyuN);
   fListOfHistos->Add(fMultiBin1ChargedMulti);
   fListOfHistos->Add(fMultiBin2ChargedMulti);
   fListOfHistos->Add(fMultiBin3ChargedMulti);

   fListOfHistos->Add(fMultiBin1KaonKinksPos);
   fListOfHistos->Add(fMultiBin2KaonKinksPos);
   fListOfHistos->Add(fMultiBin3KaonKinksPos);

   fListOfHistos->Add(fMultiBin1KaonKinksNeg);
   fListOfHistos->Add(fMultiBin2KaonKinksNeg);
   fListOfHistos->Add(fMultiBin3KaonKinksNeg);

   fListOfHistos->Add(fMultiBin1Vertex);
   fListOfHistos->Add(fMultiBin2Vertex);
   fListOfHistos->Add(fMultiBin3Vertex);

   fListOfHistos->Add(fVertexNet);
   fListOfHistos->Add(fhPS);
   fListOfHistos->Add(fhvtxP);
   fListOfHistos->Add(fpercentile);
   
//   fAnUtils = new AliAnalysisUtils();
//   if(! fUtils ) {
 //       fUtils = new AliAnalysisUtils();
//}
   //DefineOutput(1, TList::Class());
   PostData(1, fListOfHistos);
}
//________________________________________________________________________
void AliAnalysisKinkTaskPbPb2018::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event
  
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

  // Printf("*******new event*****");

  const AliESDVertex * vertex    =    esd->GetPrimaryVertex();
  // Number ESD tracks 
  Int_t nESDTracks =  esd->GetNumberOfTracks();

  //Printf("number of tracks 1 = %f", nESDTracks);
  fMultiplicityBeforeCuts->Fill(nESDTracks);

  // check incomplete events

  fEventCuts.SetManualMode();
  fEventCuts.SetupPbPb2018();
  if (!fEventCuts.AcceptEvent(event)) return;
  
  ///if (esd->IsIncompleteDAQ()) return;
  fIncompletEv ->Fill(esd->GetNumberOfTracks() );

 // check of Pileup   3/2/2016
  ///if (esd->IsPileupFromSPD()) return;
  
  ///if (esd->IsPileupFromSPDInMultBins()) return;
  
  //Printf("number of tracks 2 = %f", nESDTracks);

 // if(fAnUtils->IsSPDClusterVsTrackletBG(esd)) return ;
  fMultiplicityAfterPileup->Fill(nESDTracks);

//  Bool_t SPDvsClustersBG = kFALSE;
  
//  AliAnalysisUtils *AnalysisUtils = new AliAnalysisUtils();
//  if (!AnalysisUtils)
//    {
//      cout<<"------- No AnalysisUtils Object Found --------"<<AnalysisUtils<<endl;
//      return;
//    }
 // else
 //   SPDvsClustersBG = AnalysisUtils->IsSPDClusterVsTrackletBG(fESD); // We want NO BG
//if(fUtils->IsSPDClusterVsTrackletBG(esd)) return ;
//  if (SPDvsClustersBG==kTRUE) return;
  fZMainVx->Fill(vertex->GetZ());

//==================check of Physics selection?

       // Multiplicity selection

      AliMultSelection *fMultSel = (AliMultSelection*) esd -> FindListObject("MultSelection");

      //if(!fMultSel-> IsEventSelected()) return;
      
      if (!fMultSel) { 
  //If you get this warning please check that the AliMultSelectionTask actually ran (before your task) 
  	AliWarning("AliMultSelection object not found!"); 
	}
      Float_t lPercentile= fMultSel->GetMultiplicityPercentile("V0M", kFALSE);
      // if ((lPercentile < 0) || (lPercentile > 100)) return;
      
      //if (!fMultSel->GetThisEventIsNotPileup()) return;
      
      //if (!fMultSel->GetThisEventIsNotPileupInMultBins()) return;
     
      //if (!fMultSel->GetThisEventINELgtZERO()) return; //this cuts...
 
      Bool_t isTracklet = (AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTracklets, 1.0) >= 1);

      //if(isTracklet== kFALSE) return;
   
      fpercentile->Fill(lPercentile);

      fESDtrackCuts->SetMaxDCAToVertexXYPtDep("0.0105 + 0.0350/pt^1.01");
      fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36);
      
      Bool_t isSelected =
	((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7;
	//((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMCEGA;

      //      if ( isSelected ==kFALSE) return;   //  24/6/11 apo MF
//*/


       // end Multiplicity selection
       //Printf("multiplicity percentile = %f", lPercentile);
//
      fMultiplicityAfterTriggerBit->Fill(nESDTracks);

      Float_t mbin=0.;
      if ((lPercentile>0)&&(lPercentile<10)) mbin=0.5;
      if ((lPercentile>=30)&&(lPercentile<50)) mbin=1.5;
      if ((lPercentile>=50)&&(lPercentile<90)) mbin=2.5;


      fhPS->Fill(mbin);
      
      //  Printf("number of tracks 3 = %f", nESDTracks);

 // trigger cut, 14-6-2017
      
     // TString triggerClasses = InputEvent()->GetFiredTriggerClasses();
      //if( triggerClasses.Contains( "EJ1" )) return;

      ///const AliESDVertex * vertex    =    esd->GetPrimaryVertex();
      Bool_t isVtxGood = vertex->GetStatus() &&
				selectVertex2015pp(esd,kTRUE,kFALSE,kTRUE); // requires Tracks and spd vertex, and Zconsistency of 5mm

      //if(isVtxGood==kFALSE) return;
      fhvtxP->Fill(mbin);
      fMultiplicityAfterVertexCut->Fill(nESDTracks);
      ///fVertexNet->Fill(vertex->GetZ());

      //double vertex_z = vertex->GetZ();
      //Bool_t isVtxInZCut=0;
      //if (TMath::Abs(vertex->GetZ())>=10.) return;
      fVertexNet->Fill(vertex->GetZ());
      //if (!selectVertex2015pp(esd,kTRUE,kFALSE,kTRUE)) return;  //new common    
      //const AliESDVertex *vertex=GetEventVertex(esd);    // 22/8
      ///fZMainVx->Fill(vertex->GetZ());

      Double_t vpos[3];
      vertex->GetXYZ(vpos);
      
      Double_t vtrack[3], ptrack[3];

      	     // // Multiplicity bins
	     
	     // if ((lPercentile>0)&&(lPercentile<1)) {
	     //   fMultiBin1ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=1)&&(lPercentile<5)) {
	     //   fMultiBin2ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=5)&&(lPercentile<10)) {
	     //   fMultiBin3ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=10)&&(lPercentile<15)) {
	     //   fMultiBin4ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=15)&&(lPercentile<20)) {
	     //   fMultiBin5ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=20)&&(lPercentile<30)) {
	     //   fMultiBin6ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=30)&&(lPercentile<40)) {
	     //   fMultiBin7ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=40)&&(lPercentile<50)) {
	     //   fMultiBin8ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=50)&&(lPercentile<70)) {
	     //   fMultiBin9ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
	     // if ((lPercentile>=70)&&(lPercentile<100)) {
	     //   fMultiBin10ChargedMulti->Fill(esd->GetNumberOfTracks());
	     // }
      
      Int_t nESDTracK = 0;
      Int_t nMultiTrack = 0;
      // Int_t nESDTrKink = 0;
      
      Double_t nsigmall = 100.0;
      Double_t nsigma = 100.0;
      Double_t nsigmaPion =-100.0;
      //     Double_t nsigmaDau  =-100.0;
      Double_t dEdxKinkDau =0.0;
      //      Double_t KinkDauCl   =0.0;
      //* apo Eftihi 
      if(!fPIDResponse) {
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler =
	  (AliInputEventHandler*)(man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();
      }

      for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
	
	AliESDtrack* track = esd->GetTrack(iTracks);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTracks);
	  continue;
	}
	//Printf("track pt = %f", track->Pt());
	fTrackPtAll->Fill(track->Pt());	 
	
	//    sigmas
	nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	//  nsigmaPion= (fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
	nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	
	//=======================new 
	//*   test back 2015
	Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
	Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
	if (track->GetTPCNclsF()>0) {
	  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
	  fRatioCrossedRows->Fill(ratioCrossedRowsOverFindableClustersTPC);
	}
	//*/
	//_______
	
	Int_t tpcNCl = track->GetTPCclusters(0);  
	Double_t tpcSign = track->GetSign();  
	
	Int_t label = track->GetLabel();
	label = TMath::Abs(label);
	
	UInt_t status=track->GetStatus();
	
	if((status&AliESDtrack::kITSrefit)==0) continue;   
	if((status&AliESDtrack::kTPCrefit)==0) continue;
	if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;  
	
	Double_t extCovPos[15];
	track->GetExternalCovariance(extCovPos);    
	
	
	track->GetXYZ(vtrack);
	// fXvYv->Fill(vtrack[0],vtrack[1]);  
	// fZvYv->Fill(vtrack[1],vtrack[2]);  
	// fZvXv->Fill(vtrack[0],vtrack[2]);  
	
	// track momentum, rapidity calculation
	track->GetPxPyPz(ptrack);
	
	TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
	
	// K-rapidity calculation 
	Double_t etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  ); //assuming Kaon
	Double_t rapiditK = 0.5 * (TMath::Log((etracK + ptrack[2]) / ( etracK - ptrack[2])  ))  ;
	
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
    
	 //fDCAz->Fill(dcaToVertexZpos);
 
//  14/2/13 /================/   
	 //if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;   //    
         if (TMath::Abs(dcaToVertexZpos)>2.5) continue;
	 if (!fESDtrackCuts->AcceptTrack(track)) continue;
	 //fDCAxy->Fill(dcaToVertexXYpos);
//
    
//  track Mult. after selection 
	 nMultiTrack++;        
  //    
//=========================================

      }
      // Printf("n multi track = %d", nMultiTrack);      
            	     // Multiplicity bins
      
      if ((lPercentile>0)&&(lPercentile<10)) {
	fMultiBin1ChargedMulti->Fill(nMultiTrack);
	fMultiBin1Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=30)&&(lPercentile<50)) {
	fMultiBin2ChargedMulti->Fill(nMultiTrack);
	fMultiBin2Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=50)&&(lPercentile<90)) {
	fMultiBin3ChargedMulti->Fill(nMultiTrack);
	fMultiBin3Vertex->Fill(vertex->GetZ());
      }
      
      // for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
	
      // 	AliESDtrack* trackD = esd->GetTrack(iTrack);
      // 	if (!trackD) {
      // 	  Printf("ERROR: Could not receive track %d", iTrack);
      // 	  continue;
      // 	}
      // 	//
      // 	Int_t indexKinkDau=trackD->GetKinkIndex(0);
      // 	// daughter kink 
      // 	//	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkDau)-1);
      // 	if ( indexKinkDau > 0 )    {
      // 	  Int_t labelD = trackD->GetLabel();
      // 	  labelD = TMath::Abs(labelD);
      // 	  //       mss 2015 nsigmaPion     = (fPIDResponse->NumberOfSigmasTPC(trackD  , AliPID::kPion));// 26/10 eftihis
      // 	  nsigmaPion     = (fPIDResponse->NumberOfSigmasTPC(trackD  , AliPID::kPion));// 26/10 eftihis
      // 	  dEdxKinkDau =  (trackD->GetTPCsignal()  )  ;  //  daughter kink  dEdx 
      // 	}
      // 	//if((indexKinkDau >0)&& (nsigmaPion>1.2)) fSignalDaughterVsDaughterMomSelectedKaons->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
      // 	// Ayto mexri 26/11/2012     if(indexKinkDau >0) fSignalDaughterVsDaughterMomSelectedKaons->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
      // } //end of loop on kink daughters
      
      
      // track loop
      //
      for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
	
	AliESDtrack* track = esd->GetTrack(iTracks);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTracks);
	  continue;
	}
	
	//fTrackPtAll->Fill(track->Pt());	 
	
	//    sigmas
	nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	//  nsigmaPion= (fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
	nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	
	//=======================new 
	//*   test back 2015
	Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
	Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
	if (track->GetTPCNclsF()>0) {
	  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
	  //fRatioCrossedRows->Fill(ratioCrossedRowsOverFindableClustersTPC);
	}
	//*/
	//_______
	
	Int_t indexKinkPos=track->GetKinkIndex(0);   // kink index 
	
	Int_t tpcNCl = track->GetTPCclusters(0);  
	Double_t tpcSign = track->GetSign();  
	
	// Int_t label = track->GetLabel();
	// label = TMath::Abs(label);
	
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
	
	// K-rapidity calculation 
	Double_t etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  ); //assuming Kaon
	Double_t rapiditK = 0.5 * (TMath::Log((etracK + ptrack[2]) / ( etracK - ptrack[2])  ))  ;
	
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
    
	 fDCAz->Fill(dcaToVertexZpos);
 
//  14/2/13 /================/   
	 //if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;   //    
	 if (!fESDtrackCuts->AcceptTrack(track)) continue;
	 if (TMath::Abs(dcaToVertexZpos)>2.5) continue;
	 fDCAxy->Fill(dcaToVertexXYpos);
//
    
//  track Mult. after selection 
	 nESDTracK++;        
  //    
//=========================================
	 fTrackPtAfterTrackCuts->Fill(track->Pt());
	 
	 //  select kinks
	 if(indexKinkPos<0) {     ////mother kink
	   
	   fPtAllKinks->Fill(track->Pt());  // Pt from tracks , all kinks
	   
	   fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);
	   
	   // select kink class	
	   
	   AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
	   
	   // DCA kink
	   Double_t  Dist2 = kink->GetDistance();
	   //   if (Dist2 > 0.2) continue; //  systematics 11/8/11 
	   
	   // TPC mother momentum 
	   
	   const TVector3 vposKink(kink->GetPosition()); //reco position
	   fxyKinkPosition ->Fill( vposKink[0], vposKink[1]);
	   Double_t  dzKink=vpos[2]-vposKink[2]; 
	   //
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

           fHistQtAll->Fill(qT) ;  //  Qt   distr
	   
           fPtFromMotherAllKinks->Fill(motherMfromKink.Pt()); /// pt from kink
	   
           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);
	    
	   if((TMath::Abs(rapiditK)) > fRapiK) continue; //  allagh  Nov. 2014 , better acceptance 
	   if ((track->Pt())<.200)continue;  // new Feb 2012

	   if (TMath::Abs(trackEta) > 0.8) continue;  // new  NOv   2014 
	   
	   fQtVsKinkMomAfterAcceptance->Fill(track->P(),qT);
	   fQtAfterAcceptance->Fill(qT); 
	   
	   fEtaAfterAcceptance->Fill(trackEta) ;  //   Eta distr 
	               
	   fKinkMomFromMother->Fill(motherPt);     
	   
	   //          maximum decay angle at a given mother momentum

	   Double_t maxDecAngKmu=f1->Eval(track->P(),0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(track->P(),0.,0.,0.);
	   
	   //  fake kinks are removed 
	   if( (kinkAngle<2.)  ) continue;
	             
	   //  BG  ?????==============
	   if (TMath::Abs(vposKink[2]) >  225.) continue ;
           if (TMath::Abs(vposKink[2]) <  0.5) continue ;
	   
	   fPtKaonInTPC ->Fill(trackPt);     
	   fAngleMomKaonsinTPC->Fill(track->P(), kinkAngle); 
	   //
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
           fQtInvMassKaonTPC -> Fill (invariantMassKmu, qT);
           fInvMassMuNuKaonTPC->Fill(invariantMassKmu);
           fInvMassPtKaonTPC ->Fill(invariantMassKmu, trackPt);
	   //
	   fRadiusVsPtKaonTPC->Fill( kink->GetR(), trackPt); // 
	   //  radius and Minv selection 
	   //   if( ( kink->GetR()> 120 ) && ( kink->GetR() < 210 )  )  {
	   if( ( kink->GetR()> fKinkRadLow ) && ( kink->GetR() <fKinkRadUp   )  )  {
	     //  for systematics   if( ( kink->GetR()> 130 ) && ( kink->GetR() < 200 )  )  {
	     if (qT>fLowQt) {
	       fAngleVsMomentumKaonsInR->Fill(track->P(), kinkAngle); 
	       fInvMassKaonInR->Fill(invariantMassKmu);
	       fRadiusVsNclInR->Fill( (kink->GetR()), (track->GetTPCclusters(0))) ;

	     }
	   }
	   
	   //  tails cleaning
	   if(  (tpcNCl<fLowCluster)) continue;  // test 27 feb 2012 ,, OK

	   // cleaning BG in tails
	   Int_t tpcNClHigh = -31.67+ (11./12.)*(kink->GetR()) ;  
	   if (tpcNCl > tpcNClHigh) continue;   
	   
	   Int_t tpcNClMin  = -85.5 + (65./95.)*(kink->GetR()) ;  
	   if ( tpcNCl < tpcNClMin ) continue;   
	   
	   fKaonKinkPtAfterKinkNclCut->Fill(track->Pt());  // ALL  K-candidates until now                 

	   //Final kaon kinks selection 
	   if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<fRapiK)&&(invariantMassKmu<0.8)){
	  
	     fAngleVsMomPreSelectedKinks->Fill(track->P(), kinkAngle); 
	     fPtPreSelectedkinks->Fill(trackPt);
	     
	     //  maximum angles selection with some error cut
	     //Whyyyyyyyyy????????
	     if ((kinkAngle<maxDecAngpimu*1.2)) continue; 
	     if ((kinkAngle>maxDecAngKmu*.98) && (track->P()>1.2)) continue;  ///5/5/2010
	     
	     fPtKinkBeforedEdx->Fill(trackPt);     
	     //  here the kaons selected by the decay features
	     fTPCSgnlPa->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()));
	     //
	     //  NO dEdx cut test 9/2/13               if ( nsigma               > 3.5) continue;
	     //  mss 2015
	     if (nsigma > 3.) continue; 
	     // 
	     //  next plots for the identified kaons by the kink analysis
	     
	     fTPCSignalVsMomSelectedKaons->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal())) ;
	     fNclustersVsRadiusSelectedKaons->Fill( (kink->GetR()), (track->GetTPCclusters(0))) ;
	     fPtVsRadiusSelectedKaons->Fill(kink->GetR(), trackPt); // 
	     fPtVsInvMassSelectedKaons->Fill(invariantMassKmu, trackPt);
	     
	     fdedxMthVsTPCMomSelectedKaons->Fill(trMomTPC, (track->GetTPCsignal())) ;               //trMomTPC
	     //fSignalMthVsSignalDaughterSelectedKaons ->Fill(trMomTPCKink  , (track->GetTPCsignal()  ) ) ;
	     //fSignalMthVsSignalDaughterSelectedKaons ->Fill(dEdxKinkDau, (track->GetTPCsignal())) ;
	     //
	     //fSignalDaughterVsDaughterMomSelectedKaons->Fill( daughterMKink.Mag(), dEdxKinkDau  ) ;  //  daughter kink 
	     // 
	     fNsigmaVsTPCmomSelectedKaons->Fill(trMomTPC ,nsigmall );     
	     fNsigmaSelectedKaons->Fill(nsigmall );     
	     //
	     fRadiusSelectedKinks->Fill(kink->GetR());  // kink 
	     fKaonLifetimeSelectedKaons->Fill(lifeKink);  // kink 
             fTrackEtaSelectedKaons->Fill(trackEta);
	     fRapiditySelectedKaons->Fill(rapiditK);  //  rapidityof kaons 
	     fZKinkProductionVsKinkRadSelectedKaons->Fill( vposKink[2], kink->GetR() );
	     
	     Float_t signPt= tpcSign*trackPt;
	     fNclinTPCVsSignedPtSelectedKaons->Fill(signPt, tpcNCl);   ///  28/4/2010
	     fRapidityVsSignedPtSelectedKons->Fill(signPt, rapiditK);
	     fNclVsRapiditySelectedKAons->Fill( rapiditK, tpcNCl);
	     fSignedPtSelectedKaons->Fill(signPt);
	     fNclVsChi2SelectedKaons->Fill((track->GetTPCchi2()), tpcNCl);
	     fChi2perTPCclusterSelectedKaons-> Fill((track->GetTPCchi2()/track->GetTPCclusters(0))) ;
	     //   fDCAxy->Fill(dcaToVertexXYpos);
	     //if( dEdxKinkDau> 1.5* (track->GetTPCsignal()   )  )      fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
	     //    if((dEdxKinkDau>  80. ) && (dEdxKinkDau > 4.*nsigmaPion)   )      fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
	     //  mss 2015          if (nsigmaPion>  3.             ) fdedxMthVsTPCMomSelectedKaons->Fill( daughterMKink.Mag(),  dEdxKinkDau    ) ;
	     //if (TMath::Abs(dEdxKinkDau -(track->GetTPCsignal() )> 10. )) fdedxMthVsTPCMomSelectedKaons->Fill( daughterMKink.Mag(),  dEdxKinkDau    ) ;
	     fLifetimeSelectedKaons->Fill((lifeKink*.493667)/track->P());
             fPtSelectedKaons->Fill(track->Pt());        
	     fDCAkinkSelectedKaons->Fill(Dist2);
	     
	     if(tpcSign >0.) fPtPositiveSelectedKaons->Fill(track->Pt());   //K-plus bins Comb 
	     if(tpcSign <0.) fPtNegativeSelectedKaons->Fill(track->Pt());   //K-minus bins Comb   
	     
	     fPtKinkK0->Fill(track->Pt()); ///  K0    bins     
	     if(tpcSign >0.) fPtKinkK0P->Fill(track->Pt());   //K-plus bins K0 Peter 
	     if(tpcSign <0.) fPtKinkK0N ->Fill(track->Pt());   //K-minus bins K0 Peter

	     // Multiplicity bins
	     
	     if ((lPercentile>0)&&(lPercentile<10)) {
	       if (tpcSign >0.) fMultiBin1KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin1KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=30)&&(lPercentile<50)) {
	       if (tpcSign >0.) fMultiBin2KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin2KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=50)&&(lPercentile<90)) {
	       if (tpcSign >0.) fMultiBin3KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin3KaonKinksNeg->Fill(track->Pt());
	     }

	     // End Multiplicity bins
	     
	     fPtKinkGyu->Fill(track->Pt()); ///   K-charged High Pt Gyula
	     if(tpcSign >0.)        fPtKinkGyuP->Fill( track->Pt()    ) ;   //K-plus charged high pt ,Gyula
	     if(tpcSign <0.)        fPtKinkGyuN ->Fill( track->Pt()  ) ;   //K-minus bins High Pt charged Gyula  
	     
             fKinKRbn->Fill(track->Pt());       // TOF      
	     
	     //fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK); 
	     fAngMomK->Fill(    track->P(),        kinkAngle); 
	     fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
	     fPosiKinKXZ->Fill( vposKink[2], vposKink[0]  );
	     fPosiKinKYZ->Fill( vposKink[2], vposKink[1]  );
	     
	   }  //  kink selection 
	   

	 }  //End Kink Information    
  

   } //track loop 
   
   PostData(1, fListOfHistos);
   
}      

//________________________________________________________________________
void AliAnalysisKinkTaskPbPb2018::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}

//____________________________________________________________________//

  // Get the vertex 
const AliESDVertex* AliAnalysisKinkTaskPbPb2018::GetEventVertex(const AliESDEvent* esd) const
  
{
  // Get the vertex 
  
   const AliESDVertex* vertex = esd->GetPrimaryVertexTracks();

  if((vertex->GetStatus()==kTRUE)) return vertex;
  else
  { 
     vertex = esd->GetPrimaryVertexSPD();
      if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>0)) return vertex;
     else
     return 0;
  }
}
//======================================
//
//
Bool_t AliAnalysisKinkTaskPbPb2018::selectVertex2015pp(AliESDEvent *esd,
						       //Bool_t AliAnalysisTaskExtractV0me::selectVertex2015pp(AliESDEvent *esd,
						       Bool_t checkSPDres, //enable check on vtx resolution
						       Bool_t requireSPDandTrk, //ask for both trk and SPD vertex
						       Bool_t checkProximity) //apply cut on relative position of spd and trk 
//verteces
{
  
  if (!esd) return kFALSE;
  
  const AliESDVertex * trkVertex = esd->GetPrimaryVertexTracks();
 const AliESDVertex * spdVertex = esd->GetPrimaryVertexSPD();
 Bool_t hasSPD = spdVertex->GetStatus();
 Bool_t hasTrk = trkVertex->GetStatus();
 //
 //Note that AliVertex::GetStatus checks that N_contributors is > 0
 //reject events if both are explicitly requested and none is available
 if (requireSPDandTrk && !(hasSPD && hasTrk)) return kFALSE;
 
 //reject events if none between the SPD or track verteces are available
 //if no trk vertex, try to fall back to SPD vertex;
 if (!hasTrk) {
   if (!hasSPD) return kFALSE;
   //on demand check the spd vertex resolution and reject if not 
   //satisfied
   if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
   
 } else {
   if (hasSPD) {
     //if enabled check the spd vertex resolution and reject if not 
     ////satisfied
     //       //if enabled, check the proximity between the spd vertex and trak 
     //       //vertex, and reject if not satisfied
     if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
     if ((checkProximity && TMath::Abs(spdVertex->GetZ() - 
				       trkVertex->GetZ())>0.5)) return kFALSE;
   }
 }
 
 //Cut on the vertex z position
 // const  AliESDVertex * vertex = esd->GetPrimaryVertex();
 //if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
 return kTRUE;
}

             //___________________________________________________________________________________ 
//___________________________________________________________________________________

            //Bool_t AliAnalysisTaskExtractV0me::IsGoodSPDvertexRes(const AliESDVertex  * spdVertex)
            Bool_t AliAnalysisKinkTaskPbPb2018::IsGoodSPDvertexRes(const AliESDVertex  * spdVertex)
              {
            if (!spdVertex) return kFALSE;
                if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04
          && spdVertex->GetZRes()<0.25)) return kFALSE;
                              return kTRUE;
                                  }


