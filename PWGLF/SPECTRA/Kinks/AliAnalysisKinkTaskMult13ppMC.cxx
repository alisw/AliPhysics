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
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
//#include "AliAnalysisUtils.h"
//#include "AliESDUtils.h"
//#include "AliPPVsMultUtils.h"
//#include "AliOADBContainer.h"
//#include "AliOADBMultSelection.h"
//#include "AliMultEstimator.h"
//#include "AliMultVariable.h"
//#include "AliMultInput.h"
#include "AliMultSelection.h"
//#include "AliPythia.h"
/////////#include "AliTENDERSupplies.h"
//#include "AliMultSelectionCuts.h"
#include "AliAnalysisKinkTaskMult13ppMC.h"
ClassImp(AliAnalysisKinkTaskMult13ppMC)


//________________________________________________________________________
AliAnalysisKinkTaskMult13ppMC::AliAnalysisKinkTaskMult13ppMC(const char *name) 
  : AliAnalysisTaskSE(name)
  , fMultiplicityBeforeCuts(0)  
  , fIncompletEv(0)   
  , fMultiplicityAfterTriggerBit(0)  
  , fMultiplicityAfterPileup(0)   
  , fZMainVx(0)  
  , fMultiplicityAfterVertexCut(0)   
  , f1(0)
  , f2(0)
  , fListOfHistos(0)
  , fKinkRadUp(210.)
  , fKinkRadLow(120.)
  , fLowCluster(20)
  , fLowQt(.12)
  , fRapiK(0.5)
  , fptKMC(0)  
  , fPtKPlMC(0)  
  , fPtKMnMC(0)   
  , frapidKMC(0)  
  , flenTrRef(0)  
  , flifeSmall(0)   
  , fLifeMCProcess13(0)  
  , fLifeP(0)    
  , fLengthRZ(0)   
  , fLengthP(0)    
  , flengthKMC(0)   
  , fLifeMCProcess4(0)    
  //  , fradPtRapMC(0)
  , fLifeDECAYmuonORelectron(0)
  , fmaxAngleVsMomKmuKel(0)  
  //  , fradiusPtRapidityKplusMu(0)
  , fradiusKplusMu(0)  
  , fQtKplusMu(0)   
  , fPtRapidityKplusMu(0)  
  , fLifetimeKplusMu(0)   
  //  , fradiusPtRapidityKminusMu(0)
  , fradiusKminusMu(0) 
  , fQtKminusMu(0)  
  , fPtRapidityKminusMu(0)  
  , fLifetimeKminusMu(0)
  , fMultiBin1KaonMCPos(0)
  , fMultiBin2KaonMCPos(0)
  , fMultiBin3KaonMCPos(0)
  , fMultiBin4KaonMCPos(0)
  , fMultiBin5KaonMCPos(0)
  , fMultiBin6KaonMCPos(0)
  , fMultiBin7KaonMCPos(0)
  , fMultiBin8KaonMCPos(0)
  , fMultiBin9KaonMCPos(0)
  , fMultiBin10KaonMCPos(0)
  , fMultiBin1KaonMCNeg(0)
  , fMultiBin2KaonMCNeg(0)
  , fMultiBin3KaonMCNeg(0)
  , fMultiBin4KaonMCNeg(0)
  , fMultiBin5KaonMCNeg(0)
  , fMultiBin6KaonMCNeg(0)
  , fMultiBin7KaonMCNeg(0)
  , fMultiBin8KaonMCNeg(0)
  , fMultiBin9KaonMCNeg(0)
  , fMultiBin10KaonMCNeg(0)
  , flifetime(0)
  , fMultiBin1fPtKplusMu(0)
  , fMultiBin1fPtKminusMu(0)
  , fMultiBin1fPtKplusEl(0)
  , fMultiBin1fPtKminusEl(0)
  , fMultiBin1fPtKPiPlus(0)
  , fMultiBin1fPtKPiMinus(0)
  , fMultiBin2fPtKplusMu(0)
  , fMultiBin2fPtKminusMu(0)
  , fMultiBin2fPtKplusEl(0)
  , fMultiBin2fPtKminusEl(0)
  , fMultiBin2fPtKPiPlus(0)
  , fMultiBin2fPtKPiMinus(0)
  , fMultiBin3fPtKplusMu(0)
  , fMultiBin3fPtKminusMu(0)
  , fMultiBin3fPtKplusEl(0)
  , fMultiBin3fPtKminusEl(0)
  , fMultiBin3fPtKPiPlus(0)
  , fMultiBin3fPtKPiMinus(0)
  , fMultiBin4fPtKplusMu(0)
  , fMultiBin4fPtKminusMu(0)
  , fMultiBin4fPtKplusEl(0)
  , fMultiBin4fPtKminusEl(0)
  , fMultiBin4fPtKPiPlus(0)
  , fMultiBin4fPtKPiMinus(0)
  , fMultiBin5fPtKplusMu(0)
  , fMultiBin5fPtKminusMu(0)
  , fMultiBin5fPtKplusEl(0)
  , fMultiBin5fPtKminusEl(0)
  , fMultiBin5fPtKPiPlus(0)
  , fMultiBin5fPtKPiMinus(0)
  , fMultiBin6fPtKplusMu(0)
  , fMultiBin6fPtKminusMu(0)
  , fMultiBin6fPtKplusEl(0)
  , fMultiBin6fPtKminusEl(0)
  , fMultiBin6fPtKPiPlus(0)
  , fMultiBin6fPtKPiMinus(0)
  , fMultiBin7fPtKplusMu(0)
  , fMultiBin7fPtKminusMu(0)
  , fMultiBin7fPtKplusEl(0)
  , fMultiBin7fPtKminusEl(0)
  , fMultiBin7fPtKPiPlus(0)
  , fMultiBin7fPtKPiMinus(0)
  , fMultiBin8fPtKplusMu(0)
  , fMultiBin8fPtKminusMu(0)
  , fMultiBin8fPtKplusEl(0)
  , fMultiBin8fPtKminusEl(0)
  , fMultiBin8fPtKPiPlus(0)
  , fMultiBin8fPtKPiMinus(0)
  , fMultiBin9fPtKplusMu(0)
  , fMultiBin9fPtKminusMu(0)
  , fMultiBin9fPtKplusEl(0)
  , fMultiBin9fPtKminusEl(0)
  , fMultiBin9fPtKPiPlus(0)
  , fMultiBin9fPtKPiMinus(0)
  , fMultiBin10fPtKplusMu(0)
  , fMultiBin10fPtKminusMu(0)
  , fMultiBin10fPtKplusEl(0)
  , fMultiBin10fPtKminusEl(0)
  , fMultiBin10fPtKPiPlus(0)
  , fMultiBin10fPtKPiMinus(0)
  , fradiusKplusEl(0)
  , fQtKplusEl(0)
  , fPtRapidityKplusEl(0)
  , fLifetimeKplusEl(0)
  , fradiusKminusEl(0)
  , fQtKminusEl(0)
  , fPtRapidityKminusEl(0)
  , fLifetimeKminusEl(0)
  , fradiusKPiPlus(0)
  , fQtKPiPlus(0)
  , fPtRapidityKPiPlus(0)
  , fLifetimeKPiPlus(0)
  , fradiusKPiMinus(0)
  , fQtKPiMinus(0)
  , fPtRapidityKPiMinus(0)
  , fLifetimeKPiMinus(0)
  , fPIDResponse()
  , fTrackPtAll(0)
  , fRatioCrossedRows(0)
  , fZvXv(0)
  , fZvYv(0)
  , fXvYv(0)
  , fDCAz(0)
  , fDCAxy(0)
  , fMultiBin1ChargedMulti(0)
  , fMultiBin2ChargedMulti(0)
  , fMultiBin3ChargedMulti(0)
  , fMultiBin4ChargedMulti(0)
  , fMultiBin5ChargedMulti(0)
  , fMultiBin6ChargedMulti(0)
  , fMultiBin7ChargedMulti(0)
  , fMultiBin8ChargedMulti(0)
  , fMultiBin9ChargedMulti(0)
  , fMultiBin10ChargedMulti(0)
  , fMultiBin1Vertex(0)
  , fMultiBin2Vertex(0)
  , fMultiBin3Vertex(0)
  , fMultiBin4Vertex(0)
  , fMultiBin5Vertex(0)
  , fMultiBin6Vertex(0)
  , fMultiBin7Vertex(0)
  , fMultiBin8Vertex(0)
  , fMultiBin9Vertex(0)
  , fMultiBin10Vertex(0)
  , fTrackPtAfterTrackCuts(0)
  , fPtAllKinks(0)
  , fRatioCrossedRowsKink(0)
  , fxyKinkPosition(0)
  , fHistQtAll(0)
  , fPtFromMotherAllKinks(0)
  , fQtVsKinkMomAfterAcceptance(0)
  , fQtAfterAcceptance(0)
  , fEtaAfterAcceptance(0)
  , fKinkMomFromMother(0)
  , fQtBeforeAngleCut(0)
  , fFakepipi(0)
  , fFakeKPi(0)
  , fPtKaonInTPC(0)
  , fAngleMomKaonsinTPC(0)
  , fRadiusPtPion(0)
  , fRadiusPtKaon(0)
  , fQtKMu(0)
  , fQtKPi(0)
  , fQtKEl(0)
  , fQtK3PiP(0)
  , fQtK3PiM(0)
  , fPtKaonPDG(0)
  , fPtKaonPDGPos(0)
  , fPtKaonPDGNeg(0)
  , fKaonPDGEta(0)
  , fKaonPDGrapidity(0)
  , fKaonPDGpTvsRadius(0)
  , fKaonPDGqT(0)
  , fAngMomKaonPDG(0)
  , fAngMomPi(0)
  , fQtInvMassKaonTPC(0)
  , fInvMassPtKaonTPC(0)
  , fRadiusVsPtKaonTPC(0)
  , fSignPtNcl(0)
  , fAngleVsMomentumKaonsInR(0)
  , fInvMassKaonInR(0)
  , fRadiusVsNclInR(0)
  , fkaonToMu(0)
  , fkaonToPi(0)
  , fkaonToKa(0)
  , fRadiusNcl(0)
  , fcodeMotherVsDaughter(0)
  , fZMotherVsDaughter(0)
  , fPtVsRadiusFake(0)
  , fKaonKinkPtAfterKinkNclCut(0)
  , fAngleVsMomPreSelectedKinks(0)
  , fPtPreSelectedkinks(0)
  , fPtKinkBeforedEdx(0)
  , fTPCSgnlPa(0)
  , fTPCSignalVsMomSelectedKaons(0)
  , fNclustersVsRadiusSelectedKaons(0)
  , fPtVsRadiusSelectedKaons(0)
  , fPtVsInvMassSelectedKaons(0)
  , fPtKaonRECpos(0)
  , fPtKaonRECneg(0)
  , fPtKaonREC(0)
  , fdedxMthVsTPCMomSelectedKaons(0)
  , fSignalMthVsSignalDaughterSelectedKaons(0)
  , fSignalDaughterVsDaughterMomSelectedKaons(0)
  , fNsigmaVsTPCmomSelectedKaons(0)
  , fRadiusSelectedKinks(0)
  , fKaonLifetimeSelectedKaons(0)
  , fTrackEtaSelectedKaons(0)
  , fRapiditySelectedKaons(0)
  , fZKinkProductionVsKinkRadSelectedKaons(0)
  , fNclinTPCVsSignedPtSelectedKaons(0)
  , fRapidityVsSignedPtSelectedKons(0)
  , fNclVsRapiditySelectedKAons(0)
  , fNclVsChi2SelectedKaons(0)
  , fChi2perTPCclusterSelectedKaons(0)
  , fLifetimeSelectedKaons()
  , fPtSelectedKaons(0)
  , fDCAkinkSelectedKaons(0)
  , fPtPositiveSelectedKaons(0)
  , fPtNegativeSelectedKaons(0)
  , fkinkKaonPDG(0)
  , fkinkKaonPDGpos(0)
  , fkinkKaonPDGneg(0)
  , fkinkKaonPDGposMulti1(0)
  , fkinkKaonPDGnegMulti1(0)
  , fkinkKaonPDGposMulti2(0)
  , fkinkKaonPDGnegMulti2(0)
  , fkinkKaonPDGposMulti3(0)
  , fkinkKaonPDGnegMulti3(0)
  , fkinkKaonPDGposMulti4(0)
  , fkinkKaonPDGnegMulti4(0)
  , fkinkKaonPDGposMulti5(0)
  , fkinkKaonPDGnegMulti5(0)
  , fkinkKaonPDGposMulti6(0)
  , fkinkKaonPDGnegMulti6(0)
  , fkinkKaonPDGposMulti7(0)
  , fkinkKaonPDGnegMulti7(0)
  , fkinkKaonPDGposMulti8(0)
  , fkinkKaonPDGnegMulti8(0)
  , fkinkKaonPDGposMulti9(0)
  , fkinkKaonPDGnegMulti9(0)
  , fkinkKaonPDGposMulti10(0)
  , fkinkKaonPDGnegMulti10(0)
  , fkinkKaonPDGBkg(0)
  , fkinkKaonPDGBkgMulti1(0)
  , fkinkKaonPDGBkgMulti2(0)
  , fkinkKaonPDGBkgMulti3(0)
  , fkinkKaonPDGBkgMulti4(0)
  , fkinkKaonPDGBkgMulti5(0)
  , fkinkKaonPDGBkgMulti6(0)
  , fkinkKaonPDGBkgMulti7(0)
  , fkinkKaonPDGBkgMulti8(0)
  , fkinkKaonPDGBkgMulti9(0)
  , fkinkKaonPDGBkgMulti10(0)
  , fMultiBin1KaonKinksPos(0)
  , fMultiBin2KaonKinksPos(0)
  , fMultiBin3KaonKinksPos(0)
  , fMultiBin4KaonKinksPos(0)
  , fMultiBin5KaonKinksPos(0)
  , fMultiBin6KaonKinksPos(0)
  , fMultiBin7KaonKinksPos(0)
  , fMultiBin8KaonKinksPos(0)
  , fMultiBin9KaonKinksPos(0)
  , fMultiBin10KaonKinksPos(0)
  , fMultiBin1KaonKinksNeg(0)
  , fMultiBin2KaonKinksNeg(0)
  , fMultiBin3KaonKinksNeg(0)
  , fMultiBin4KaonKinksNeg(0)
  , fMultiBin5KaonKinksNeg(0)
  , fMultiBin6KaonKinksNeg(0)
  , fMultiBin7KaonKinksNeg(0)
  , fMultiBin8KaonKinksNeg(0)
  , fMultiBin9KaonKinksNeg(0)
  , fMultiBin10KaonKinksNeg(0)
  , fAngMomK(0)
  , fPosiKinkK(0)
  , fPosiKinKXZ(0)
  , fPosiKinKYZ(0)
{
  // Constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkTaskMult13ppMC::UserCreateOutputObjects() 
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


   fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated", 1500, 0.0, 15.0); 
   fPtKPlMC= new TH1F("fPtKPlMC", "P_{T}Kaon Pos generated", 1500, 0.0, 15.0);  
   fPtKMnMC= new TH1F("fPtKMnMC", "P_{T}Kaon Minus generated", 1500, 0.0, 15.0);   
   frapidKMC = new TH1F("frapidKMC ", "rapidity distrib MC ", 26, -1.3, 1.3);  
   flenTrRef =new TH1F("flenTrRef","lifetime ref K Decay",100,0.,1000.0);  
   flifeSmall=new TH1F("flifeSmall","lifetime ref K Decay",100,0.,1000.0);  
   fLifeMCProcess13 =new TH1F("fLifeMCProcess13", "lifetime ref K Decay",100,0.,1000.0);  
   fLifeP =new TH1F("fLifeP", "lifetime ref K   Decay",100,0.,1000.0);   
   fLengthRZ =new TH1F("fLengthRZ","lifetime ref K   Decay",100,0.,1000.0); 
   fLengthP =new TH1F("fLengthP","lifetime ref K   Decay",100,0.,1000.0);  
   flengthKMC =new TH1F("flengthKMC","lifetime ref K   Decay",100,0.,1000.0);  
   fLifeMCProcess4 =new TH1F("fLifeMCProcess4", "lifetime ref K Decay",100,0.,1000.0);   
   //   fradPtRapMC=new TH3F("fradPtRapMC","rad pt rap dat",28,100.,240., 1500, 0., 15., 20, -1., 1. );
   fLifeDECAYmuonORelectron =new TH1F("fLifeDECAYmuonORelectron", "lifetime ref K   Decay   ",100,0.,1000.0); 
   fmaxAngleVsMomKmuKel =new TH2F("fmaxAngleVsMomKmuKel","Decay angle vrs Mother Mom,Kmu",1500,0.0,15.0,120,0.,120.);  

   //   fradiusPtRapidityKplusMu=new TH3F("fradiusPtRapidityKplusMu","rad pt rap dat",120,90.,250., 1500, 0., 15., 20, -1., 1.);
   fradiusKplusMu =new TH1F("fradiusKplusMu", "radius  K generated",100,0.,1000.);  
   fQtKplusMu =new TH1F("fQtKplusMu", "Q_{T} distribution  K to mu MC",100, 0.0,.300);  
   fPtRapidityKplusMu =new TH2F("fPtRapidityKplusMu", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5); 
   fLifetimeKplusMu =new TH1F("fLifetimeKplusMu", "lifetime ref K", 100, 0., 1000.0); 

   //   fradiusPtRapidityKminusMu=new TH3F("fradiusPtRapidityKminusMu","rad pt rap dat",120,90.,250., 1500, 0., 15., 20, -1., 1.);
   fradiusKminusMu =new TH1F("fradiusKminusMu", "radius  K generated",100,0.,1000.);  
   fQtKminusMu =new TH1F("fQtKminusMu", "Q_{T} distribution  K to mu MC",100, 0.0,.300);  
   fPtRapidityKminusMu =new TH2F("fPtRapidityKminusMu", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5);
   fLifetimeKminusMu =new TH1F("fLifetimeKminusMu", "lifetime ref K", 100, 0., 1000.0);  

   fMultiBin1KaonMCPos = new TH1F("fMultiBin1KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonMCPos = new TH1F("fMultiBin2KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonMCPos = new TH1F("fMultiBin3KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin4KaonMCPos = new TH1F("fMultiBin4KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin5KaonMCPos = new TH1F("fMultiBin5KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin6KaonMCPos = new TH1F("fMultiBin6KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin7KaonMCPos = new TH1F("fMultiBin7KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin8KaonMCPos = new TH1F("fMultiBin8KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin9KaonMCPos = new TH1F("fMultiBin9KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin10KaonMCPos = new TH1F("fMultiBin10KaonMCPos", "P_{T}Kaon distribution",1500, 0.0,15.0);

   fMultiBin1KaonMCNeg = new TH1F("fMultiBin1KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonMCNeg = new TH1F("fMultiBin2KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonMCNeg = new TH1F("fMultiBin3KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin4KaonMCNeg = new TH1F("fMultiBin4KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin5KaonMCNeg = new TH1F("fMultiBin5KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin6KaonMCNeg = new TH1F("fMultiBin6KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin7KaonMCNeg = new TH1F("fMultiBin7KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin8KaonMCNeg = new TH1F("fMultiBin8KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin9KaonMCNeg = new TH1F("fMultiBin9KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin10KaonMCNeg = new TH1F("fMultiBin10KaonMCNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   flifetime =new TH1F("flifetime","lifetime ref K Decay",100,0.,1000.0);

   fMultiBin1fPtKplusMu= new TH1F("fMultiBin1fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin1fPtKminusMu= new TH1F("fMultiBin1fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin1fPtKplusEl= new TH1F("fMultiBin1fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin1fPtKminusEl= new TH1F("fMultiBin1fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin1fPtKPiPlus= new TH1F("fMultiBin1fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin1fPtKPiMinus= new TH1F("fMultiBin1fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin2fPtKplusMu= new TH1F("fMultiBin2fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin2fPtKminusMu= new TH1F("fMultiBin2fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin2fPtKplusEl= new TH1F("fMultiBin2fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin2fPtKminusEl= new TH1F("fMultiBin2fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin2fPtKPiPlus= new TH1F("fMultiBin2fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin2fPtKPiMinus= new TH1F("fMultiBin2fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin3fPtKplusMu= new TH1F("fMultiBin3fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin3fPtKminusMu= new TH1F("fMultiBin3fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin3fPtKplusEl= new TH1F("fMultiBin3fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin3fPtKminusEl= new TH1F("fMultiBin3fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin3fPtKPiPlus= new TH1F("fMultiBin3fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin3fPtKPiMinus= new TH1F("fMultiBin3fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin4fPtKplusMu= new TH1F("fMultiBin4fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin4fPtKminusMu= new TH1F("fMultiBin4fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin4fPtKplusEl= new TH1F("fMultiBin4fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin4fPtKminusEl= new TH1F("fMultiBin4fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin4fPtKPiPlus= new TH1F("fMultiBin4fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin4fPtKPiMinus= new TH1F("fMultiBin4fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin5fPtKplusMu= new TH1F("fMultiBin5fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin5fPtKminusMu= new TH1F("fMultiBin5fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin5fPtKplusEl= new TH1F("fMultiBin5fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin5fPtKminusEl= new TH1F("fMultiBin5fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin5fPtKPiPlus= new TH1F("fMultiBin5fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin5fPtKPiMinus= new TH1F("fMultiBin5fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin6fPtKplusMu= new TH1F("fMultiBin6fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin6fPtKminusMu= new TH1F("fMultiBin6fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin6fPtKplusEl= new TH1F("fMultiBin6fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin6fPtKminusEl= new TH1F("fMultiBin6fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin6fPtKPiPlus= new TH1F("fMultiBin6fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin6fPtKPiMinus= new TH1F("fMultiBin6fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin7fPtKplusMu= new TH1F("fMultiBin7fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin7fPtKminusMu= new TH1F("fMultiBin7fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin7fPtKplusEl= new TH1F("fMultiBin7fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin7fPtKminusEl= new TH1F("fMultiBin7fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin7fPtKPiPlus= new TH1F("fMultiBin7fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin7fPtKPiMinus= new TH1F("fMultiBin7fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin8fPtKplusMu= new TH1F("fMultiBin8fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin8fPtKminusMu= new TH1F("fMultiBin8fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin8fPtKplusEl= new TH1F("fMultiBin8fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin8fPtKminusEl= new TH1F("fMultiBin8fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin8fPtKPiPlus= new TH1F("fMultiBin8fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin8fPtKPiMinus= new TH1F("fMultiBin8fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin9fPtKplusMu= new TH1F("fMultiBin9fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin9fPtKminusMu= new TH1F("fMultiBin9fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin9fPtKplusEl= new TH1F("fMultiBin9fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin9fPtKminusEl= new TH1F("fMultiBin9fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin9fPtKPiPlus= new TH1F("fMultiBin9fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin9fPtKPiMinus= new TH1F("fMultiBin9fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fMultiBin10fPtKplusMu= new TH1F("fMultiBin10fPtKplusMu", " ", 1500, 0.0, 15.0);
   fMultiBin10fPtKminusMu= new TH1F("fMultiBin10fPtKminusMu", " ", 1500, 0.0, 15.0); 
   fMultiBin10fPtKplusEl= new TH1F("fMultiBin10fPtKplusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin10fPtKminusEl= new TH1F("fMultiBin10fPtKminusEl", " ", 1500, 0.0, 15.0); 
   fMultiBin10fPtKPiPlus= new TH1F("fMultiBin10fPtKPiPlus", " ", 1500, 0.0, 15.0);
   fMultiBin10fPtKPiMinus= new TH1F("fMultiBin10fPtKPiMinus", " ", 1500, 0.0, 15.0);

   fradiusKplusEl =new TH1F("fradiusKplusEl", "radius  K generated",100,0.,1000.);
   fQtKplusEl =new TH1F("fQtKplusEl", "Q_{T} distribution  K to mu MC",100, 0.0,.300); 
   fPtRapidityKplusEl =new TH2F("fPtRapidityKplusEl", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5);
   fLifetimeKplusEl =new TH1F("fLifetimeKplusEl", "lifetime ref K", 100, 0., 1000.0);

   fradiusKminusEl =new TH1F("fradiusKminusEl", "radius  K generated",100,0.,1000.);
   fQtKminusEl =new TH1F("fQtKminusEl", "Q_{T} distribution  K to mu MC",100, 0.0,.300); 
   fPtRapidityKminusEl =new TH2F("fPtRapidityKminusEl", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5);
   fLifetimeKminusEl =new TH1F("fLifetimeKminusEl", "lifetime ref K", 100, 0., 1000.0);
   fradiusKPiPlus =new TH1F("fradiusKPiPlus", "radius  K generated",100,0.,1000.);
   fQtKPiPlus =new TH1F("fQtKPiPlus", "Q_{T} distribution  K to mu MC",100, 0.0,.300); 
   fPtRapidityKPiPlus =new TH2F("fPtRapidityKPiPlus", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5);
   fLifetimeKPiPlus =new TH1F("fLifetimeKPiPlus", "lifetime ref K", 100, 0., 1000.0);

   fradiusKPiMinus =new TH1F("fradiusKPiMinus", "radius  K generated",100,0.,1000.);
   fQtKPiMinus =new TH1F("fQtKPiMinus", "Q_{T} distribution  K to mu MC",100, 0.0,.300); 
   fPtRapidityKPiMinus =new TH2F("fPtRapidityKPiMinus", "Q_{T} distribution  K to mu MC", 1500, 0., 15., 30, -1.5, 1.5);
   fLifetimeKPiMinus =new TH1F("fLifetimeKPiMinus", "lifetime ref K", 100, 0., 1000.0);
   fTrackPtAll = new TH1F("fTrackPtAll", "P_{T} distribution",1500, 0.0,15.0);
   fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);

   fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-1.5,1.5,60, -15., 15.0);
   fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-1.5,1.5, 60, -15., 15.);
   fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
   fDCAz = new TH1D("fDCAz", "rad distribution  PID pr",100,-10.0, 10.0);
   fDCAxy = new TH1D("fDCAxy", "dca  distribution PID  ",20,-1.,1.);
   fMultiBin1ChargedMulti= new TH1F("fMultiBin1ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin2ChargedMulti= new TH1F("fMultiBin2ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin3ChargedMulti= new TH1F("fMultiBin3ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin4ChargedMulti= new TH1F("fMultiBin4ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin5ChargedMulti= new TH1F("fMultiBin5ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin6ChargedMulti= new TH1F("fMultiBin6ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin7ChargedMulti= new TH1F("fMultiBin7ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin8ChargedMulti= new TH1F("fMultiBin8ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin9ChargedMulti= new TH1F("fMultiBin9ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin10ChargedMulti= new TH1F("fMultiBin10ChargedMulti", "charge multiplicity ESD",500, 0.0,500.0);
   fMultiBin1Vertex= new TH1D("fZMainVxfMultiBin1Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin2Vertex= new TH1D("fZMainVxfMultiBin2Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin3Vertex= new TH1D("fZMainVxfMultiBin3Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin4Vertex= new TH1D("fZMainVxfMultiBin4Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin5Vertex= new TH1D("fZMainVxfMultiBin5Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin6Vertex= new TH1D("fZMainVxfMultiBin6Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin7Vertex= new TH1D("fZMainVxfMultiBin7Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin8Vertex= new TH1D("fZMainVxfMultiBin8Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin9Vertex= new TH1D("fZMainVxfMultiBin9Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fMultiBin10Vertex= new TH1D("fZMainVxfMultiBin10Vertex", "ESD charge mult. Main Vertex", 60,-15.,15.);
   fTrackPtAfterTrackCuts = new TH1F("fTrackPtAfterTrackCuts", "P_{T} distribution",1500, 0.0,15.0);
   fPtAllKinks= new TH1F("fPtAllKinks", "P_{T} all kinks",1500, 0.0,15.0);
   fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);
   fxyKinkPosition= new TH2F("fxyKinkPosition", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
   fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300);
   fPtFromMotherAllKinks= new TH1F("fPtFromMotherAllKinks", "P_{T}",1500, 0.0,15.0);
   fQtVsKinkMomAfterAcceptance = new TH2F("fQtVsKinkMomAfterAcceptance", " Qt vrs Mother P", 1500, 0., 15.0, 100, 0., 0.300);
   fQtAfterAcceptance= new TH1F("fQtAfterAcceptance", "Q_{T} distribution",100, 0.0,.300);
   fEtaAfterAcceptance= new TH1F("fEtaAfterAcceptance", "Eta distribution", 26,-1.3, 1.3);
   fKinkMomFromMother= new TH1F("fKinkMomFromMother", "P_{T}Kaon kinks backgr",1500, 0.0,15.0);
   fQtBeforeAngleCut= new TH1F("fQtBeforeAngleCut", "Q_{T} distribution",100, 0.0,.300);
   fFakepipi = new TH1F("fFakepipi", "P_{T}fake pipi   ", 1500, 0.0, 15.0); 
   fFakeKPi = new TH1F("fFakeKPi", "P_{T}fake Kpi   ", 1500, 0.0, 15.0);
   fPtKaonInTPC = new TH1F("fPtKaonInTPC", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fAngleMomKaonsinTPC= new TH2F("fAngleMomKaonsinTPC","Decay angle vrs Mother Mom,Pi",1500, 0.0, 15.0 ,80, 0., 80.);
   fRadiusPtPion =new TH2F("fRadiusPtPion","radius vs pt Pion PDG ",120, 90.,250.,1500, 0.,15.);
   fRadiusPtKaon =new TH2F("fRadiusPtKaon","radius vs pt Kaon PDG ",120, 90.,250.,1500, 0.,15.);
   fQtKMu= new TH1F("fQtKMu", "Q_{T} distribution  K to mu ",100, 0.0,.300); 
   fQtKPi= new TH1F("fQtKPi", "Q_{T} distribution K to pi",100, 0.0,.300); 
   fQtKEl= new TH1F("fQtKEl", "Q_{T} distribution   K to elec",100, 0.0,.300);
   fQtK3PiP= new TH1F("fQtK3PiP", "Q_{T} distribution K to 3pi ",100, 0.0,.300); 
   fQtK3PiM= new TH1F("fQtK3PiM", "Q_{T} distribution K to 3pi ",100, 0.0,.300);
   fPtKaonPDG= new TH1F("fPtKaonPDG", "  ", 1500, 0.0, 15.0);
   fPtKaonPDGPos= new TH1F("fPtKaonPDGPos", "  ", 1500, 0.0, 15.0);
   fPtKaonPDGNeg= new TH1F("fPtKaonPDGNeg", "  ", 1500, 0.0, 15.0);
   fKaonPDGEta= new TH1F("fKaonPDGEta", "Eta distribution", 26,-1.3, 1.3);
   fKaonPDGrapidity= new TH1F("fKaonPDGrapidity", "rapidity distribution", 26,-1.3, 1.3);
   fKaonPDGpTvsRadius =new TH2F("fKaonPDGpTvsRadius","radius vs pt  ",120, 90.,250.,1500, 0.,15.);
   fKaonPDGqT= new TH1F("fKaonPDGqT", "Q_{T} distribution",100, 0.0,.300);
   fAngMomKaonPDG= new TH2F("fAngMomKaonPDG","Decay angle vrs Mother Mom,K",1500,0.0,15.0,80,0.,80.);
   fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",1500,0.0,15.0,80,0.,80.);
   fQtInvMassKaonTPC= new TH2F("fQtInvMassKaonTPC", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300);
   fInvMassPtKaonTPC =new TH2F("fInvMassPtKaonTPC","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 1500, 0.0, 15.0);
   fRadiusVsPtKaonTPC =new TH2F("fRadiusVsPtKaonTPC","radius vs pt  ",120, 90.,250.,1500, 0.,15.);
   fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",3000,-15.,15.0, 70, 0, 160);
   fAngleVsMomentumKaonsInR= new TH2F("fAngleVsMomentumKaonsInR","Decay angle vrs Mother Mom,K",1500, 0.0, 15.0, 80, 0., 80.);
   fInvMassKaonInR= new TH1F("fInvMassKaonInR","Invar m(kaon) from kink->mu+netrino decay",600,0.10, 0.7);
   fRadiusVsNclInR = new TH2F("fRadiusVsNclInR","kink radius vrs Nclust,K",120,90.,250., 80,0, 160);
   fkaonToMu= new TH1F("fkaonToMu","Invar m(kaon) from kink->mu+netrino decay",1000, 0., 1.);
   fkaonToPi= new TH1F("fkaonToPi","Invar m(kaon) from kink-> decay",1000,0.0, 1.0);
   fkaonToKa= new TH1F("fkaonToKa","Invar m(kaon) from kink-> decay",1000,0.0, 1.0);
   fRadiusNcl = new TH2F("fRadiusNcl","KinkRadius Ncl,K",120,90.,250., 80, 0, 160);
   fcodeMotherVsDaughter= new TH2F("fcodeMotherVsDaughter", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
   fZMotherVsDaughter = new TH2F("fZMotherVsDaughter", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fPtVsRadiusFake= new TH2F("fPtVsRadiusFake","radius vs pt Pion Fake ",120, 90.,250.,1500, 0.,15.);
   fKaonKinkPtAfterKinkNclCut = new TH1F("fKaonKinkPtAfterKinkNclCut", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fAngleVsMomPreSelectedKinks = new TH2F("fAngleVsMomPreSelectedKinks","Decay angle vrs Mother Mom,Kinks",1500 ,0.0, 15.0,100,0.,100.);
   fPtPreSelectedkinks = new TH1F("fPtPreSelectedkinks", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fPtKinkBeforedEdx = new TH1F("fPtKinkBeforedEdx", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fTPCSgnlPa= new TH2F("fTPCSgnlPa","TPC signal de/dx Mom,K",1000,0.0,20.,150, 0.,300.);
   fTPCSignalVsMomSelectedKaons = new TH2F("fTPCSignalVsMomSelectedKaons","TPC signal de/dx Mom,K",1000,0.0,20.0,150,0.,300.);
   fNclustersVsRadiusSelectedKaons = new TH2F("fNclustersVsRadiusSelectedKaons","kink radius vrs Nclust,K Clean ",120,90.,250., 80,0, 160);
   fPtVsRadiusSelectedKaons =new TH2F("fPtVsRadiusSelectedKaons","radius vs pt clean ",120, 90.,250.,1500, 0.,15.);
   fPtVsInvMassSelectedKaons =new TH2F("fPtVsInvMassSelectedKaons","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 1500, 0.0, 15.0);
   fPtKaonRECpos= new TH1F("fPtKaonRECpos", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fPtKaonRECneg= new TH1F("fPtKaonRECneg", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fPtKaonREC= new TH1F("fPtKaonREC", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fdedxMthVsTPCMomSelectedKaons = new TH2F("fdedxMthVsTPCMomSelectedKaons","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,100, 0., 250.);
   fSignalMthVsSignalDaughterSelectedKaons  = new TH2F("fSignalMthVsSignalDaughterSelectedKaons","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.0,100, 0., 250.);
   fSignalDaughterVsDaughterMomSelectedKaons = new TH2F("fSignalDaughterVsDaughterMomSelectedKaons","TPC signal de/dx Mom,K",500,0.0,10.0,100,0.,250.);
   fNsigmaVsTPCmomSelectedKaons = new TH2F("fNsigmaVsTPCmomSelectedKaons","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,20 , -10., 10.);
   fRadiusSelectedKinks= new TH1F("fRadiusSelectedKinks", "radius  K ",100, 0.,1000.0);
   fKaonLifetimeSelectedKaons= new TH1F("fKaonLifetimeSelectedKaons", "Length of   K generated",100,0.,1000.);
   fTrackEtaSelectedKaons= new TH1F("fTrackEtaSelectedKaons", "EtaK distribution", 26,-1.3, 1.3);
   fRapiditySelectedKaons=new TH1F("fRapiditySelectedKaons","rapid Kdistribution", 26,-1.3, 1.3);
   fZKinkProductionVsKinkRadSelectedKaons= new TH2F("fZKinkProductionVsKinkRadSelectedKaons", "z vrx kink rad    ",100, -300.0,300.0,100,100., 300.);
   fNclinTPCVsSignedPtSelectedKaons= new TH2F("fNclinTPCVsSignedPtSelectedKaons","SignPt vrs Ncl,K",3000,-15.,15.0,80,0,160);
   fRapidityVsSignedPtSelectedKons= new TH2F("fRapidityVsSignedPtSelectedKons","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
   fNclVsRapiditySelectedKAons= new TH2F("fNclVsRapiditySelectedKAons","Eta vrs nclust,K",30,-1.5,1.5, 80, 0, 160);
   fNclVsChi2SelectedKaons= new TH2F("fNclVsChi2SelectedKaons","Chi2vrs nclust,K",100,0.,500., 80,0, 160);
   fChi2perTPCclusterSelectedKaons= new TH1F("fChi2perTPCclusterSelectedKaons","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
   fLifetimeSelectedKaons= new TH1F("fLifetimeSelectedKaons", "ct study of K-kinks",100,0.,1000.);
   fPtSelectedKaons = new TH1F("fPtSelectedKaons", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fDCAkinkSelectedKaons = new TH1F("fDCAkinkSelectedKaons ", "DCA kink vetrex ",50, 0.0,1.0);
   fPtPositiveSelectedKaons= new TH1F("fPtPositiveSelectedKaons", "Pos P_{T}Kaon Kink  distribution, counts", 1500, 0.0, 15.0); 
   fPtNegativeSelectedKaons= new TH1F("fPtNegativeSelectedKaons", "Neg P_{T}Kaon Kink  distribution, counts", 1500, 0.0, 15.0);
   fkinkKaonPDG= new TH1F("fkinkKaonPDG", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fkinkKaonPDGpos= new TH1F("fkinkKaonPDGpos", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fkinkKaonPDGneg= new TH1F("fkinkKaonPDGneg", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);

   fkinkKaonPDGposMulti1= new TH1F("fkinkKaonPDGposMulti1", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti1= new TH1F("fkinkKaonPDGnegMulti1", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti2= new TH1F("fkinkKaonPDGposMulti2", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti2= new TH1F("fkinkKaonPDGnegMulti2", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti3= new TH1F("fkinkKaonPDGposMulti3", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti3= new TH1F("fkinkKaonPDGnegMulti3", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti4= new TH1F("fkinkKaonPDGposMulti4", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti4= new TH1F("fkinkKaonPDGnegMulti4", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti5= new TH1F("fkinkKaonPDGposMulti5", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti5= new TH1F("fkinkKaonPDGnegMulti5", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti6= new TH1F("fkinkKaonPDGposMulti6", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti6= new TH1F("fkinkKaonPDGnegMulti6", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti7= new TH1F("fkinkKaonPDGposMulti7", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti7= new TH1F("fkinkKaonPDGnegMulti7", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti8= new TH1F("fkinkKaonPDGposMulti8", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti8= new TH1F("fkinkKaonPDGnegMulti8", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti9= new TH1F("fkinkKaonPDGposMulti9", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti9= new TH1F("fkinkKaonPDGnegMulti9", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGposMulti10= new TH1F("fkinkKaonPDGposMulti10", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGnegMulti10= new TH1F("fkinkKaonPDGnegMulti10", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkg= new TH1F("fkinkKaonPDGBkg", "P_{T}Kaon Pos ESD", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti1= new TH1F("fkinkKaonPDGBkgMulti1", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti2= new TH1F("fkinkKaonPDGBkgMulti2", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti3= new TH1F("fkinkKaonPDGBkgMulti3", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti4= new TH1F("fkinkKaonPDGBkgMulti4", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti5= new TH1F("fkinkKaonPDGBkgMulti5", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti6= new TH1F("fkinkKaonPDGBkgMulti6", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti7= new TH1F("fkinkKaonPDGBkgMulti7", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti8= new TH1F("fkinkKaonPDGBkgMulti8", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti9= new TH1F("fkinkKaonPDGBkgMulti9", " ", 1500, 0.0, 15.0);
   fkinkKaonPDGBkgMulti10= new TH1F("fkinkKaonPDGBkgMulti10", " ", 1500, 0.0, 15.0);

   fMultiBin1KaonKinksPos = new TH1F("fMultiBin1KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonKinksPos = new TH1F("fMultiBin2KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonKinksPos = new TH1F("fMultiBin3KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin4KaonKinksPos = new TH1F("fMultiBin4KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin5KaonKinksPos = new TH1F("fMultiBin5KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin6KaonKinksPos = new TH1F("fMultiBin6KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin7KaonKinksPos = new TH1F("fMultiBin7KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin8KaonKinksPos = new TH1F("fMultiBin8KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin9KaonKinksPos = new TH1F("fMultiBin9KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin10KaonKinksPos = new TH1F("fMultiBin10KaonKinksPos", "P_{T}Kaon distribution",1500, 0.0,15.0);

   fMultiBin1KaonKinksNeg = new TH1F("fMultiBin1KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin2KaonKinksNeg = new TH1F("fMultiBin2KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin3KaonKinksNeg = new TH1F("fMultiBin3KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin4KaonKinksNeg = new TH1F("fMultiBin4KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin5KaonKinksNeg = new TH1F("fMultiBin5KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin6KaonKinksNeg = new TH1F("fMultiBin6KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin7KaonKinksNeg = new TH1F("fMultiBin7KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin8KaonKinksNeg = new TH1F("fMultiBin8KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin9KaonKinksNeg = new TH1F("fMultiBin9KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin10KaonKinksNeg = new TH1F("fMultiBin10KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
   fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);

   fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",1500,0.0,15.0,80,0.,80.);
   
   fListOfHistos=new TList();
   
   fListOfHistos->Add(fMultiplicityBeforeCuts);  
   fListOfHistos->Add(fIncompletEv);  
   fListOfHistos->Add(fMultiplicityAfterTriggerBit);   
   fListOfHistos->Add(fMultiplicityAfterPileup);  
   fListOfHistos->Add(fZMainVx);  
   fListOfHistos->Add(fMultiplicityAfterVertexCut);   
//
 
   fListOfHistos->Add(fptKMC);   
   fListOfHistos->Add(fPtKPlMC);
   fListOfHistos->Add(fPtKMnMC);  
   fListOfHistos->Add(frapidKMC);  
   fListOfHistos->Add(flenTrRef);   
   fListOfHistos->Add(flifeSmall);   
   fListOfHistos->Add(fLifeMCProcess13);   
   fListOfHistos->Add(fLifeP);  
   fListOfHistos->Add(fLengthRZ);  
   fListOfHistos->Add(fLengthP);   
   fListOfHistos->Add(flengthKMC);   
   fListOfHistos->Add(fLifeMCProcess4);    
   //   fListOfHistos->Add(fradPtRapMC);
   fListOfHistos->Add(fLifeDECAYmuonORelectron);  
   fListOfHistos->Add(fmaxAngleVsMomKmuKel);  
   
   //  fListOfHistos->Add(fradiusPtRapidityKplusMu);
   fListOfHistos->Add(fradiusKplusMu);  
   fListOfHistos->Add(fQtKplusMu);  
   fListOfHistos->Add(fPtRapidityKplusMu);  
   fListOfHistos->Add(fLifetimeKplusMu);  

   //   fListOfHistos->Add(fradiusPtRapidityKminusMu);
   fListOfHistos->Add(fradiusKminusMu); 
   fListOfHistos->Add(fQtKminusMu);  
   fListOfHistos->Add(fPtRapidityKminusMu);  
   fListOfHistos->Add(fLifetimeKminusMu);

   fListOfHistos->Add(fMultiBin1KaonMCPos);
   fListOfHistos->Add(fMultiBin2KaonMCPos);
   fListOfHistos->Add(fMultiBin3KaonMCPos);
   fListOfHistos->Add(fMultiBin4KaonMCPos);
   fListOfHistos->Add(fMultiBin5KaonMCPos);
   fListOfHistos->Add(fMultiBin6KaonMCPos);
   fListOfHistos->Add(fMultiBin7KaonMCPos);
   fListOfHistos->Add(fMultiBin8KaonMCPos);
   fListOfHistos->Add(fMultiBin9KaonMCPos);
   fListOfHistos->Add(fMultiBin10KaonMCPos);
   fListOfHistos->Add(fMultiBin1KaonMCNeg);
   fListOfHistos->Add(fMultiBin2KaonMCNeg);
   fListOfHistos->Add(fMultiBin3KaonMCNeg);
   fListOfHistos->Add(fMultiBin4KaonMCNeg);
   fListOfHistos->Add(fMultiBin5KaonMCNeg);
   fListOfHistos->Add(fMultiBin6KaonMCNeg);
   fListOfHistos->Add(fMultiBin7KaonMCNeg);
   fListOfHistos->Add(fMultiBin8KaonMCNeg);
   fListOfHistos->Add(fMultiBin9KaonMCNeg);
   fListOfHistos->Add(fMultiBin10KaonMCNeg);
   fListOfHistos->Add(flifetime);

   fListOfHistos->Add(fMultiBin1fPtKplusMu);
   fListOfHistos->Add(fMultiBin1fPtKminusMu);
   fListOfHistos->Add(fMultiBin1fPtKplusEl);
   fListOfHistos->Add(fMultiBin1fPtKminusEl);
   fListOfHistos->Add(fMultiBin1fPtKPiPlus);
   fListOfHistos->Add(fMultiBin1fPtKPiMinus);

   fListOfHistos->Add(fMultiBin2fPtKplusMu);
   fListOfHistos->Add(fMultiBin2fPtKminusMu);
   fListOfHistos->Add(fMultiBin2fPtKplusEl);
   fListOfHistos->Add(fMultiBin2fPtKminusEl);
   fListOfHistos->Add(fMultiBin2fPtKPiPlus);
   fListOfHistos->Add(fMultiBin2fPtKPiMinus);

   fListOfHistos->Add(fMultiBin3fPtKplusMu);
   fListOfHistos->Add(fMultiBin3fPtKminusMu);
   fListOfHistos->Add(fMultiBin3fPtKplusEl);
   fListOfHistos->Add(fMultiBin3fPtKminusEl);
   fListOfHistos->Add(fMultiBin3fPtKPiPlus);
   fListOfHistos->Add(fMultiBin3fPtKPiMinus);

   fListOfHistos->Add(fMultiBin4fPtKplusMu);
   fListOfHistos->Add(fMultiBin4fPtKminusMu);
   fListOfHistos->Add(fMultiBin4fPtKplusEl);
   fListOfHistos->Add(fMultiBin4fPtKminusEl);
   fListOfHistos->Add(fMultiBin4fPtKPiPlus);
   fListOfHistos->Add(fMultiBin4fPtKPiMinus);

   fListOfHistos->Add(fMultiBin5fPtKplusMu);
   fListOfHistos->Add(fMultiBin5fPtKminusMu);
   fListOfHistos->Add(fMultiBin5fPtKplusEl);
   fListOfHistos->Add(fMultiBin5fPtKminusEl);
   fListOfHistos->Add(fMultiBin5fPtKPiPlus);
   fListOfHistos->Add(fMultiBin5fPtKPiMinus);

   fListOfHistos->Add(fMultiBin6fPtKplusMu);
   fListOfHistos->Add(fMultiBin6fPtKminusMu);
   fListOfHistos->Add(fMultiBin6fPtKplusEl);
   fListOfHistos->Add(fMultiBin6fPtKminusEl);
   fListOfHistos->Add(fMultiBin6fPtKPiPlus);
   fListOfHistos->Add(fMultiBin6fPtKPiMinus);

   fListOfHistos->Add(fMultiBin7fPtKplusMu);
   fListOfHistos->Add(fMultiBin7fPtKminusMu);
   fListOfHistos->Add(fMultiBin7fPtKplusEl);
   fListOfHistos->Add(fMultiBin7fPtKminusEl);
   fListOfHistos->Add(fMultiBin7fPtKPiPlus);
   fListOfHistos->Add(fMultiBin7fPtKPiMinus);

   fListOfHistos->Add(fMultiBin8fPtKplusMu);
   fListOfHistos->Add(fMultiBin8fPtKminusMu);
   fListOfHistos->Add(fMultiBin8fPtKplusEl);
   fListOfHistos->Add(fMultiBin8fPtKminusEl);
   fListOfHistos->Add(fMultiBin8fPtKPiPlus);
   fListOfHistos->Add(fMultiBin8fPtKPiMinus);

   fListOfHistos->Add(fMultiBin9fPtKplusMu);
   fListOfHistos->Add(fMultiBin9fPtKminusMu);
   fListOfHistos->Add(fMultiBin9fPtKplusEl);
   fListOfHistos->Add(fMultiBin9fPtKminusEl);
   fListOfHistos->Add(fMultiBin9fPtKPiPlus);
   fListOfHistos->Add(fMultiBin9fPtKPiMinus);

   fListOfHistos->Add(fMultiBin10fPtKplusMu);
   fListOfHistos->Add(fMultiBin10fPtKminusMu);
   fListOfHistos->Add(fMultiBin10fPtKplusEl);
   fListOfHistos->Add(fMultiBin10fPtKminusEl);
   fListOfHistos->Add(fMultiBin10fPtKPiPlus);
   fListOfHistos->Add(fMultiBin10fPtKPiMinus);
   
   fListOfHistos->Add(fradiusKplusEl);
   fListOfHistos->Add(fQtKplusEl);
   fListOfHistos->Add(fPtRapidityKplusEl);
   fListOfHistos->Add(fLifetimeKplusEl);

   fListOfHistos->Add(fradiusKminusEl);
   fListOfHistos->Add(fQtKminusEl);
   fListOfHistos->Add(fPtRapidityKminusEl);
   fListOfHistos->Add(fLifetimeKminusEl);

   fListOfHistos->Add(fradiusKPiPlus);
   fListOfHistos->Add(fQtKPiPlus);
   fListOfHistos->Add(fPtRapidityKPiPlus);
   fListOfHistos->Add(fLifetimeKPiPlus);

   fListOfHistos->Add(fradiusKPiMinus);
   fListOfHistos->Add(fQtKPiMinus);
   fListOfHistos->Add(fPtRapidityKPiMinus);
   fListOfHistos->Add(fLifetimeKPiMinus);
   fListOfHistos->Add(fTrackPtAll);
   fListOfHistos->Add(fRatioCrossedRows);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fDCAz);
   fListOfHistos->Add(fDCAxy);
   fListOfHistos->Add(fMultiBin1ChargedMulti);
   fListOfHistos->Add(fMultiBin2ChargedMulti);
   fListOfHistos->Add(fMultiBin3ChargedMulti);
   fListOfHistos->Add(fMultiBin4ChargedMulti);
   fListOfHistos->Add(fMultiBin5ChargedMulti);
   fListOfHistos->Add(fMultiBin6ChargedMulti);
   fListOfHistos->Add(fMultiBin7ChargedMulti);
   fListOfHistos->Add(fMultiBin8ChargedMulti);
   fListOfHistos->Add(fMultiBin9ChargedMulti);
   fListOfHistos->Add(fMultiBin10ChargedMulti);
   fListOfHistos->Add(fMultiBin1Vertex);
   fListOfHistos->Add(fMultiBin2Vertex);
   fListOfHistos->Add(fMultiBin3Vertex);
   fListOfHistos->Add(fMultiBin4Vertex);
   fListOfHistos->Add(fMultiBin5Vertex);
   fListOfHistos->Add(fMultiBin6Vertex);
   fListOfHistos->Add(fMultiBin7Vertex);
   fListOfHistos->Add(fMultiBin8Vertex);
   fListOfHistos->Add(fMultiBin9Vertex);
   fListOfHistos->Add(fMultiBin10Vertex);
   fListOfHistos->Add(fTrackPtAfterTrackCuts);
   fListOfHistos->Add(fPtAllKinks);
   fListOfHistos->Add(fRatioCrossedRowsKink);
   fListOfHistos->Add(fxyKinkPosition);
   fListOfHistos->Add(fHistQtAll);
   fListOfHistos->Add(fPtFromMotherAllKinks);
   fListOfHistos->Add(fQtVsKinkMomAfterAcceptance);
   fListOfHistos->Add(fQtAfterAcceptance);
   fListOfHistos->Add(fEtaAfterAcceptance);
   fListOfHistos->Add(fKinkMomFromMother);
   fListOfHistos->Add(fQtBeforeAngleCut);
   fListOfHistos->Add(fFakepipi);
   fListOfHistos->Add(fFakeKPi);
   fListOfHistos->Add(fPtKaonInTPC);
   fListOfHistos->Add(fAngleMomKaonsinTPC);
   fListOfHistos->Add(fRadiusPtPion);
   fListOfHistos->Add(fRadiusPtKaon);
   fListOfHistos->Add(fQtKMu);
   fListOfHistos->Add(fQtKPi);
   fListOfHistos->Add(fQtKEl);
   fListOfHistos->Add(fQtK3PiP);
   fListOfHistos->Add(fQtK3PiM);
   fListOfHistos->Add(fPtKaonPDG);
   fListOfHistos->Add(fPtKaonPDGPos);
   fListOfHistos->Add(fPtKaonPDGNeg);
   fListOfHistos->Add(fKaonPDGEta);
   fListOfHistos->Add(fKaonPDGrapidity);
   fListOfHistos->Add(fKaonPDGpTvsRadius);
   fListOfHistos->Add(fKaonPDGqT);
   fListOfHistos->Add(fAngMomKaonPDG);
   fListOfHistos->Add(fAngMomPi);
   fListOfHistos->Add(fQtInvMassKaonTPC);
   fListOfHistos->Add(fInvMassPtKaonTPC);
   fListOfHistos->Add(fRadiusVsPtKaonTPC);
   fListOfHistos->Add(fSignPtNcl);
   fListOfHistos->Add(fAngleVsMomentumKaonsInR);
   fListOfHistos->Add(fInvMassKaonInR);
   fListOfHistos->Add(fRadiusVsNclInR);
   fListOfHistos->Add(fkaonToMu);
   fListOfHistos->Add(fkaonToPi);
   fListOfHistos->Add(fkaonToKa);
   fListOfHistos->Add(fRadiusNcl);
   fListOfHistos->Add(fcodeMotherVsDaughter);
   fListOfHistos->Add(fZMotherVsDaughter);
   fListOfHistos->Add(fPtVsRadiusFake);
   fListOfHistos->Add(fKaonKinkPtAfterKinkNclCut);
   fListOfHistos->Add(fAngleVsMomPreSelectedKinks);
   fListOfHistos->Add(fPtPreSelectedkinks);
   fListOfHistos->Add(fPtKinkBeforedEdx);
   fListOfHistos->Add(fTPCSgnlPa);
   fListOfHistos->Add(fTPCSignalVsMomSelectedKaons);
   fListOfHistos->Add(fNclustersVsRadiusSelectedKaons);
   fListOfHistos->Add(fPtVsRadiusSelectedKaons);
   fListOfHistos->Add(fPtVsInvMassSelectedKaons);
   fListOfHistos->Add(fPtKaonRECpos);
   fListOfHistos->Add(fPtKaonRECneg);
   fListOfHistos->Add(fPtKaonREC);
   fListOfHistos->Add(fdedxMthVsTPCMomSelectedKaons);
   fListOfHistos->Add(fSignalMthVsSignalDaughterSelectedKaons);
   fListOfHistos->Add(fSignalDaughterVsDaughterMomSelectedKaons);
   fListOfHistos->Add(fNsigmaVsTPCmomSelectedKaons);
   fListOfHistos->Add(fRadiusSelectedKinks);
   fListOfHistos->Add(fKaonLifetimeSelectedKaons);
   fListOfHistos->Add(fTrackEtaSelectedKaons);
   fListOfHistos->Add(fRapiditySelectedKaons);
   fListOfHistos->Add(fZKinkProductionVsKinkRadSelectedKaons);
   fListOfHistos->Add(fNclinTPCVsSignedPtSelectedKaons);
   fListOfHistos->Add(fRapidityVsSignedPtSelectedKons);
   fListOfHistos->Add(fNclVsRapiditySelectedKAons);
   fListOfHistos->Add(fNclVsChi2SelectedKaons);
   fListOfHistos->Add(fChi2perTPCclusterSelectedKaons);
   fListOfHistos->Add(fLifetimeSelectedKaons);
   fListOfHistos->Add(fPtSelectedKaons);
   fListOfHistos->Add(fDCAkinkSelectedKaons);
   fListOfHistos->Add(fPtPositiveSelectedKaons);
   fListOfHistos->Add(fPtNegativeSelectedKaons);
   fListOfHistos->Add(fkinkKaonPDG);
   fListOfHistos->Add(fkinkKaonPDGpos);
   fListOfHistos->Add(fkinkKaonPDGneg);

   fListOfHistos->Add(fkinkKaonPDGposMulti1);
   fListOfHistos->Add(fkinkKaonPDGnegMulti1);
   fListOfHistos->Add(fkinkKaonPDGposMulti2);
   fListOfHistos->Add(fkinkKaonPDGnegMulti2);
   fListOfHistos->Add(fkinkKaonPDGposMulti3);
   fListOfHistos->Add(fkinkKaonPDGnegMulti3);
   fListOfHistos->Add(fkinkKaonPDGposMulti4);
   fListOfHistos->Add(fkinkKaonPDGnegMulti4);
   fListOfHistos->Add(fkinkKaonPDGposMulti5);
   fListOfHistos->Add(fkinkKaonPDGnegMulti5);
   fListOfHistos->Add(fkinkKaonPDGposMulti6);
   fListOfHistos->Add(fkinkKaonPDGnegMulti6);
   fListOfHistos->Add(fkinkKaonPDGposMulti7);
   fListOfHistos->Add(fkinkKaonPDGnegMulti7);
   fListOfHistos->Add(fkinkKaonPDGposMulti8);
   fListOfHistos->Add(fkinkKaonPDGnegMulti8);
   fListOfHistos->Add(fkinkKaonPDGposMulti9);
   fListOfHistos->Add(fkinkKaonPDGnegMulti9);
   fListOfHistos->Add(fkinkKaonPDGposMulti10);
   fListOfHistos->Add(fkinkKaonPDGnegMulti10);
   fListOfHistos->Add(fkinkKaonPDGBkg);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti1);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti2);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti3);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti4);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti5);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti6);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti7);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti8);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti9);
   fListOfHistos->Add(fkinkKaonPDGBkgMulti10);
   fListOfHistos->Add(fMultiBin1KaonKinksPos);
   fListOfHistos->Add(fMultiBin2KaonKinksPos);
   fListOfHistos->Add(fMultiBin3KaonKinksPos);
   fListOfHistos->Add(fMultiBin4KaonKinksPos);
   fListOfHistos->Add(fMultiBin5KaonKinksPos);
   fListOfHistos->Add(fMultiBin6KaonKinksPos);
   fListOfHistos->Add(fMultiBin7KaonKinksPos);
   fListOfHistos->Add(fMultiBin8KaonKinksPos);
   fListOfHistos->Add(fMultiBin9KaonKinksPos);
   fListOfHistos->Add(fMultiBin10KaonKinksPos);
   fListOfHistos->Add(fMultiBin1KaonKinksNeg);
   fListOfHistos->Add(fMultiBin2KaonKinksNeg);
   fListOfHistos->Add(fMultiBin3KaonKinksNeg);
   fListOfHistos->Add(fMultiBin4KaonKinksNeg);
   fListOfHistos->Add(fMultiBin5KaonKinksNeg);
   fListOfHistos->Add(fMultiBin6KaonKinksNeg);
   fListOfHistos->Add(fMultiBin7KaonKinksNeg);
   fListOfHistos->Add(fMultiBin8KaonKinksNeg);
   fListOfHistos->Add(fMultiBin9KaonKinksNeg);
   fListOfHistos->Add(fMultiBin10KaonKinksNeg);
   fListOfHistos->Add(fAngMomK);
   fListOfHistos->Add(fPosiKinkK);
   fListOfHistos->Add(fPosiKinKXZ);
   fListOfHistos->Add(fPosiKinKYZ);
   

   //DefineOutput(1, TList::Class());
   PostData(1, fListOfHistos);
}
//________________________________________________________________________
void AliAnalysisKinkTaskMult13ppMC::UserExec(Option_t *) 
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

  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  // Number ESD tracks 
  Int_t nESDTracks =  esd->GetNumberOfTracks();
  fMultiplicityBeforeCuts->Fill(nESDTracks);

  // check incomplete events 
  if (esd->IsIncompleteDAQ()) return;
  fIncompletEv ->Fill(esd->GetNumberOfTracks() );
//
 // check of Pileup   3/2/2016
  if (esd->IsPileupFromSPD()) return;
  if (esd->IsPileupFromSPDInMultBins()) return;
  fMultiplicityAfterPileup->Fill(nESDTracks);

//==================check of Physics selection?

       // Multiplicity selection

      AliMultSelection *fMultSel = (AliMultSelection*) esd -> FindListObject("MultSelection");

      // if(!fMultSel-> IsEventSelected()) return;
      if (!fMultSel) { 
  //If you get this warning please check that the AliMultSelectionTask actually ran (before your task) 
  	AliWarning("AliMultSelection object not found!"); 
	}
      Float_t lPercentile= fMultSel->GetMultiplicityPercentile("V0M", kFALSE);
 //     if ((lPercentile < 0) || (lPercentile > 100)) return; 

      Bool_t isTracklet = (AliESDtrackCuts::GetReferenceMultiplicity(esd, AliESDtrackCuts::kTracklets, 1.0) >= 1);

      if(isTracklet== kFALSE) return;
      
      Bool_t isSelected =
	((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kINT7;
	//((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kEMCEGA;

       if ( isSelected ==kFALSE) return;   //  24/6/11 apo MF
//*/


       // end Multiplicity selection
      Printf("multiplicity percentile = %f", lPercentile);
//
      fMultiplicityAfterTriggerBit->Fill(nESDTracks);

 // trigger cut, 14-6-2017
  
     // TString triggerClasses = InputEvent()->GetFiredTriggerClasses();
      //if( triggerClasses.Contains( "EJ1" )) return;

//
/*
//===============Marek's  multiplicity
   Float_t refmultiplicity=fCutsMul->CountAcceptedTracks(esd);
        if(fLowMulcut>-1)
        {
                if(refmultiplicity<fLowMulcut)
                        return;
        }
        if(fUpMulcut>-1)
        {
                if(refmultiplicity>fUpMulcut)
                        return;
        }



       fMultiplicityAfterPileup->Fill(refmultiplicity);
*/

      if (!selectVertex2015pp(esd,kTRUE,kFALSE,kTRUE)) return;  //new common 

       
      fMultiplicityAfterVertexCut->Fill(nESDTracks);
      
      const AliESDVertex *vertex=GetEventVertex(esd);    // 22/8
      fZMainVx->Fill(vertex->GetZ());

      Double_t vpos[3];
      vertex->GetXYZ(vpos);
      
      Double_t vtrack[3], ptrack[3];

 
      // MC part

      Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

      AliStack* stack=mcEvent->Stack();
      
      //primary tracks  in MC
      Int_t  nPrim = stack->GetNprimary();
      //
      
      // loop over mc particles
      
      // variable to count tracks
      Int_t nMCTracks = 0;
      Int_t nMCKinkKs = 0;
      
      //   15/2/12 OK   for (Int_t iMc = 0; iMc < nPrim; iMc++)
      //for (Int_t iMc = 0; iMc < stack->GetNtrack(); iMc++)
      for (Int_t iMc = 0; iMc < mcEvent->GetNumberOfTracks(); iMc++)   // ??? Is this the reco MC tracks?
	{
	  
	  TParticle* particle = stack->Particle(iMc);
	  
	  if (!particle)
	    {
	      //  AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
	      continue;
	    }
	  // keep only primaries
	  if (!stack->IsPhysicalPrimary(iMc)) continue;  //??? Is this the same with the GetNprimary?

	  Double_t ptK = particle->Pt();
	  
	  if( ptK <0.200) continue;       //    12/2/2012
	  //if( ptK <0.200) return;       //    22/4/2016  //???? Why???????????
	  //
	  Double_t EtaMC  = particle->Eta();
	  // No Eta mss 2015      if ((TMath::Abs(EtaMC)) > 0.8) continue ; //  27/11/2014 
	  if ((TMath::Abs(EtaMC)) > 0.8) continue ; //  27/11/2014 
	  //if ((TMath::Abs(EtaMC)) > 0.8) return ; //  22/4/ 2016 
	  
	  Float_t charg=0;
	  Float_t code = particle->GetPdgCode();
	  Int_t  mcProcess=-1011;
	  //---------------------------------------kaon selection 
	  if ((code==321)||(code==-321)) {	    
	    
	    Double_t   etracKMC= TMath::Sqrt(particle->P() *particle->P()  + 0.493677 *0.493677);
	    Double_t rapidiKMC = 0.5 * (TMath::Log(  (etracKMC +particle->Pz())/( etracKMC-particle->Pz())))  ;

	    if ((TMath::Abs( rapidiKMC)) > fRapiK ) continue;   // 
	    // if ( (TMath::Abs( rapidiKMC)) > fRapiK ) return;     //   22/4/2016 test
            frapidKMC ->Fill(rapidiKMC) ;  //   test for correct  acceptane    april 2016
	    //metefera to cut sthn acceptance edw - Evi, 15/9/2017
	    
	    if(code > 0 ) charg =1;
	    if(code < 0 ) charg =-1;
	    Float_t chargPt= charg*ptK;
	    
	    fptKMC->Fill(ptK);
	    //   fSignPtGen->Fill(chargPt);// kaon gensign pt
	    if (charg==1) fPtKPlMC->Fill(ptK);
	    if (charg==-1) fPtKMnMC->Fill(ptK);

	    if ((lPercentile>0)&&(lPercentile<1)) {
	      if (charg==1) fMultiBin1KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin1KaonMCNeg->Fill(ptK);
	    }
	    if ((lPercentile>=1)&&(lPercentile<5)) {
	      if (charg==1) fMultiBin2KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin2KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=5)&&(lPercentile<10)) {
	      if (charg==1) fMultiBin3KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin3KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=10)&&(lPercentile<15)) {
	      if (charg==1) fMultiBin4KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin4KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=15)&&(lPercentile<20)) {
	      if (charg==1) fMultiBin5KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin5KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=20)&&(lPercentile<30)) {
	      if (charg==1) fMultiBin6KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin6KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=30)&&(lPercentile<40)) {
	      if (charg==1) fMultiBin7KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin7KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=40)&&(lPercentile<50)) {
	      if (charg==1) fMultiBin8KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin8KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=50)&&(lPercentile<70)) {
	      if (charg==1) fMultiBin9KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin9KaonMCNeg->Fill(ptK);
	      
	    }
	    if ((lPercentile>=70)&&(lPercentile<100)) {
	      if (charg==1) fMultiBin10KaonMCPos->Fill(ptK);
	      if (charg==-1) fMultiBin10KaonMCNeg->Fill(ptK);
	      
	    }

	    //
	    //====================================metafore ths rapid gia svsth  acc?===================

	    //
	    // primary   vertex
            //      Double_t mVx=particle->Vx();
            //     Double_t mVy=particle->Vy();
	    Double_t mVz=particle->Vz();
	    // 25/11/2012  ???????back 10/1/2013
	    TClonesArray* trArray=0;
	    TParticle* tempParticle=0;
	    
	    TVector3 DecayMomentumK(0,0,0);  
	    Float_t lengthKMC=0;
	    if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) {
	      AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
	      lengthKMC = MCtrackReference->GetLength();
	      
	      //DecayMomentumK.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
	    }
	    flenTrRef->Fill(lengthKMC);
	    flifetime->Fill((lengthKMC*0.493667/particle->P()));  // 19/7
	    if ((lengthKMC>100.)&& (lengthKMC<300.)) flifeSmall->Fill((lengthKMC*0.493667/particle->P() ) ); 
	    
	    Int_t nMCKpi =0;
	    Int_t mcProc4 =0;
	    Int_t mcProc13=0;
	    Float_t radiusD=0;
	    Double_t LengthRZ =0.;
	    Double_t LengthP =0.;
	    Float_t MCQt =0.;
	    Int_t firstD=particle->GetFirstDaughter();
	    Int_t lastD=particle->GetLastDaughter();
	    
	    if ((lastD<=0)||(firstD<=0)) continue; 
	    
	    if (lastD > mcEvent->GetNumberOfTracks()) continue;
	    if (firstD > mcEvent->GetNumberOfTracks()) continue;
	    TParticle *daughter1=0x0;

	    //loop on secondaries
	    for (Int_t k=firstD;k<=lastD;k++) {
              if (k > 0)    {
		daughter1=stack->Particle(k);   // 27/8   
		Float_t dcode = daughter1->GetPdgCode();
		
		//     mother momentum trackrefs    and QtMC     // 17/9/2010,  test Feb 2011
		if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) { 
		  AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
		  DecayMomentumK.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
		} ///how do I know that this takes mother's mom????????
		const TVector3 MCP3d(daughter1->Px(), daughter1->Py(), daughter1->Pz()); //MC daughter's momentum
		MCQt = MCP3d.Perp(DecayMomentumK); //MC daughter's transverse momentum in mother's frame
		//
		Double_t MCKinkAngle = TMath::ASin(MCQt/daughter1->P() ); 
                Double_t  MCKinkAngle2= TMath::RadToDeg() * MCKinkAngle; // in degrees 
		mcProcess=daughter1->GetUniqueID();
		radiusD=daughter1->R();

		Double_t hVz=daughter1->Vz();
		
		LengthRZ = TMath::Sqrt( radiusD*radiusD  + ( mVz-hVz) * (mVz-hVz) );  //   19/7/2010 mss
		
	        LengthP  = (TMath::Abs( mVz-hVz))* (TMath::Sqrt( 1.+ ( ptK*ptK)/ (particle->Pz() * particle->Pz()) )) ;
		//	??? Poia h diafora???
		if(mcProcess==13) {
		  mcProc13++;
		  
		  if(mcProc13==1) fLifeMCProcess13->Fill((lengthKMC*0.493667/particle->P()));  // 19/7
		  
		}  
		
		
		//
		if (mcProcess==4) {        
		  
		  mcProc4++;
		  if ( mcProc4==1)  {
		    // 10/1/13                if( (radiusD >120.)&&( radiusD< 210.) )  {
		    fLifeP ->Fill((LengthP*0.493667  /particle->P()));  // 19/7
		    fLengthRZ->Fill(LengthRZ);
		    fLengthP ->Fill(LengthP);
		    flengthKMC->Fill(lengthKMC);  //
		    
		    fLifeMCProcess4 ->Fill((lengthKMC*0.493667/particle->P()));  // 19/7
		    //	    fradPtRapMC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		  }

		  //
		  if (((code==321)&&( dcode ==-13))||((code ==-321)&&(dcode== 13))||((code==321 )&&( dcode ==-11))||((code ==-321)&&(dcode== 11))) {  
		    fLifeDECAYmuonORelectron ->Fill((lengthKMC*0.493667/particle->P()));
		    if( (radiusD >fKinkRadLow )&&(radiusD< fKinkRadUp) )  
		      fmaxAngleVsMomKmuKel->Fill(particle->P(), MCKinkAngle2);// MC 
		  } 
		  
		  if (((TMath::Abs(code)==321)&&(TMath::Abs(dcode) ==211))&&(mcProc4<2)) fLifeMCProcess4->Fill(lengthKMC *0.493667 /particle->P()) ;//19/7
		  
		  // test feb 2013                    if ((TMath::Abs(hVz)<0.5) || (TMath::Abs(hVz )>225)) continue;
		  ///   inside radius region ----------------------------------------------
		  if(MCKinkAngle2 < 2.) continue;  // as in ESD 
		  if ((daughter1->R()>fKinkRadLow)&&(daughter1->R()< fKinkRadUp)&&(MCQt>fLowQt)) {
		    
		    if ((code==321)&&(dcode ==-13))   {  
		      //	      fradiusPtRapidityKplusMu->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		      fradiusKplusMu->Fill(daughter1->R());
		      fQtKplusMu ->Fill(MCQt);//to muon
		      fPtRapidityKplusMu->Fill(ptK,rapidiKMC);//to muon
		      fLifetimeKplusMu->Fill(LengthP*0.4933667/particle->P());
		      
		      
		      // Multiplicity bins
		      
		      if ((lPercentile>0)&&(lPercentile<1)) {
		      	fMultiBin1fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=1)&&(lPercentile<5)) {
		      	fMultiBin2fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=5)&&(lPercentile<10)) {
		      	fMultiBin3fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=10)&&(lPercentile<15)) {
		      	fMultiBin4fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=15)&&(lPercentile<20)) {
		      	fMultiBin5fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=20)&&(lPercentile<30)) {
		      	fMultiBin6fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=30)&&(lPercentile<40)) {
		      	fMultiBin7fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=40)&&(lPercentile<50)) {
		      	fMultiBin8fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=50)&&(lPercentile<70)) {
		      	fMultiBin9fPtKplusMu->Fill(ptK);
		      }
		      if ((lPercentile>=70)&&(lPercentile<100)) {
		      	fMultiBin10fPtKplusMu->Fill(ptK);
		      }

	     // End Multiplicity bins
		      
		    } //  positive kaon   
		    if (( code ==-321)&&(dcode== 13)) {
		      //	      fradiusPtRapidityKminusMu->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		      fradiusKminusMu->Fill(daughter1->R());
		      fQtKminusMu ->Fill(MCQt);//to muon
		      fPtRapidityKminusMu->Fill(ptK,rapidiKMC);//to muon
		      fLifetimeKminusMu->Fill(LengthP*0.4933667/particle->P());

		      // 		      // Multiplicity bins
		      
		      if ((lPercentile>0)&&(lPercentile<1)) {
		      	fMultiBin1fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=1)&&(lPercentile<5)) {
		      	fMultiBin2fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=5)&&(lPercentile<10)) {
		      	fMultiBin3fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=10)&&(lPercentile<15)) {
		      	fMultiBin4fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=15)&&(lPercentile<20)) {
		      	fMultiBin5fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=20)&&(lPercentile<30)) {
		      	fMultiBin6fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=30)&&(lPercentile<40)) {
		      	fMultiBin7fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=40)&&(lPercentile<50)) {
		      	fMultiBin8fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=50)&&(lPercentile<70)) {
		      	fMultiBin9fPtKminusMu->Fill(ptK);
		      }
		      if ((lPercentile>=70)&&(lPercentile<100)) {
		      	fMultiBin10fPtKminusMu->Fill(ptK);
		      }
		    } //  negative code
	
		    if ((code==321 )&&( dcode ==-11))   {
		      //	      fradiusPtRapidityKplusEl->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		      fradiusKplusEl->Fill(daughter1->R());
		      fQtKplusEl ->Fill(MCQt);//to muon
		      fPtRapidityKplusEl->Fill(ptK,rapidiKMC);//to muon
		      fLifetimeKplusEl->Fill(LengthP*0.4933667/particle->P());

		      // 		      // Multiplicity bins
		      
		      if ((lPercentile>0)&&(lPercentile<1)) {
		      	fMultiBin1fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=1)&&(lPercentile<5)) {
		      	fMultiBin2fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=5)&&(lPercentile<10)) {
		      	fMultiBin3fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=10)&&(lPercentile<15)) {
		      	fMultiBin4fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=15)&&(lPercentile<20)) {
		      	fMultiBin5fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=20)&&(lPercentile<30)) {
		      	fMultiBin6fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=30)&&(lPercentile<40)) {
		      	fMultiBin7fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=40)&&(lPercentile<50)) {
		      	fMultiBin8fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=50)&&(lPercentile<70)) {
		      	fMultiBin9fPtKplusEl->Fill(ptK);
		      }
		      if ((lPercentile>=70)&&(lPercentile<100)) {
		      	fMultiBin10fPtKplusEl->Fill(ptK);
		      }
		    } //  positive kaon
		    
		    if ((code ==-321)&&(dcode== 11)) {
		      // fradiusPtRapidityKminusEl->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		      fradiusKminusEl->Fill(daughter1->R());
		      fQtKminusEl ->Fill(MCQt);//to muon
		      fPtRapidityKminusEl->Fill(ptK,rapidiKMC);//to muon
		      fLifetimeKminusEl->Fill(LengthP*0.4933667/particle->P());

		      // Multiplicity bins
		      
		      if ((lPercentile>0)&&(lPercentile<1)) {
		      	fMultiBin1fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=1)&&(lPercentile<5)) {
		      	fMultiBin2fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=5)&&(lPercentile<10)) {
		      	fMultiBin3fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=10)&&(lPercentile<15)) {
		      	fMultiBin4fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=15)&&(lPercentile<20)) {
		      	fMultiBin5fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=20)&&(lPercentile<30)) {
		      	fMultiBin6fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=30)&&(lPercentile<40)) {
		      	fMultiBin7fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=40)&&(lPercentile<50)) {
		      	fMultiBin8fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=50)&&(lPercentile<70)) {
		      	fMultiBin9fPtKminusEl->Fill(ptK);
		      }
		      if ((lPercentile>=70)&&(lPercentile<100)) {
		      	fMultiBin10fPtKminusEl->Fill(ptK);
		      }
		    } //  negative code
		    
		    if (((code==321)&&(dcode ==211))||((code == -321)&&(dcode ==-211))) nMCKpi++ ; 
 
		    nMCKinkKs++;
		  } // in radius
		  
		} //    decay, process=4
	      } // positive k
	    }//  loop on daughters
	    	    
	    if ((code > 0)&&(nMCKpi == 1)) {

	      //     fradiusPtRapidityKPiPlus->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
	      fradiusKPiPlus->Fill(daughter1->R());
	      fQtKPiPlus ->Fill(MCQt);//to muon
	      fPtRapidityKPiPlus->Fill(ptK,rapidiKMC);//to muon
	      fLifetimeKPiPlus->Fill(LengthP*0.4933667/particle->P());
	      
	      // // Multiplicity bins
	      
	      if ((lPercentile>0)&&(lPercentile<1)) {
	      	fMultiBin1fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=1)&&(lPercentile<5)) {
	      	fMultiBin2fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=5)&&(lPercentile<10)) {
	      	fMultiBin3fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=10)&&(lPercentile<15)) {
	      	fMultiBin4fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=15)&&(lPercentile<20)) {
	      	fMultiBin5fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=20)&&(lPercentile<30)) {
	      	fMultiBin6fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=30)&&(lPercentile<40)) {
	      	fMultiBin7fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=40)&&(lPercentile<50)) {
	      	fMultiBin8fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=50)&&(lPercentile<70)) {
	      	fMultiBin9fPtKPiPlus->Fill(ptK);
	      }
	      if ((lPercentile>=70)&&(lPercentile<100)) {
	      	fMultiBin10fPtKPiPlus->Fill(ptK);
	      }
	      
	    }   //positive kaon 
	    //
	    if ((code < 0)&&(nMCKpi == 1)) {

	      //fradiusPtRapidityKPiMinus->Fill(radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
	      fradiusKPiMinus->Fill(daughter1->R());
	      fQtKPiMinus ->Fill(MCQt);//to muon
	      fPtRapidityKPiMinus->Fill(ptK,rapidiKMC);//to muon
	      fLifetimeKPiMinus->Fill(LengthP*0.4933667/particle->P());

	      // Multiplicity bins
	      
	      if ((lPercentile>0)&&(lPercentile<1)) {
	      	fMultiBin1fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=1)&&(lPercentile<5)) {
	      	fMultiBin2fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=5)&&(lPercentile<10)) {
	      	fMultiBin3fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=10)&&(lPercentile<15)) {
	      	fMultiBin4fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=15)&&(lPercentile<20)) {
	      	fMultiBin5fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=20)&&(lPercentile<30)) {
	      	fMultiBin6fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=30)&&(lPercentile<40)) {
	      	fMultiBin7fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=40)&&(lPercentile<50)) {
	      	fMultiBin8fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=50)&&(lPercentile<70)) {
	      	fMultiBin9fPtKPiMinus->Fill(ptK);
	      }
	      if ((lPercentile>=70)&&(lPercentile<100)) {
	      	fMultiBin10fPtKPiMinus->Fill(ptK);
	      }
	      
	    }   //negative K
	    
	  }   /// kaons  loop 
	  	  
	  nMCTracks++;
	}// end of mc particle
      
      // end MC part

      if(!fPIDResponse) {
	AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
	AliInputEventHandler* inputHandler =
	  (AliInputEventHandler*)(man->GetInputEventHandler());
	fPIDResponse = inputHandler->GetPIDResponse();
      }
      
      Int_t nESDTracK = 0;
      Int_t nMultiTrack = 0;
      Double_t nsigmall = 100.0;
      Double_t nsigma = 100.0;

      for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
	
	AliESDtrack* track = esd->GetTrack(iTracks);
	if (!track) {
	  Printf("ERROR: Could not receive track %d", iTracks);
	  continue;
	}
	
	fTrackPtAll->Fill(track->Pt());	 
	
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
	
	//	Int_t indexKinkPos=track->GetKinkIndex(0);   // kink index 
	
	Int_t tpcNCl = track->GetTPCclusters(0);  
	Double_t tpcSign = track->GetSign();  
	
	Int_t label = track->GetLabel();
	label = TMath::Abs(label);
	
	if(label > mcEvent->GetNumberOfTracks()) continue;
	if (label > nPrim) continue;
	TParticle * part = stack->Particle(label);
	if (!part) continue;
	
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
	 if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;   //    
	 
	 fDCAxy->Fill(dcaToVertexXYpos);
	 //
	 
	 //  track Mult. after selection 
	 nMultiTrack++;        
	 //
      }

           	     // Multiplicity bins
      
      if ((lPercentile>0)&&(lPercentile<1)) {
	fMultiBin1ChargedMulti->Fill(nMultiTrack);
	fMultiBin1Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=1)&&(lPercentile<5)) {
	fMultiBin2ChargedMulti->Fill(nMultiTrack);
	fMultiBin2Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=5)&&(lPercentile<10)) {
	fMultiBin3ChargedMulti->Fill(nMultiTrack);
	fMultiBin3Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=10)&&(lPercentile<15)) {
	fMultiBin4ChargedMulti->Fill(nMultiTrack);
	fMultiBin4Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=15)&&(lPercentile<20)) {
	fMultiBin5ChargedMulti->Fill(nMultiTrack);
	fMultiBin5Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=20)&&(lPercentile<30)) {
	fMultiBin6ChargedMulti->Fill(nMultiTrack);
	fMultiBin6Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=30)&&(lPercentile<40)) {
	fMultiBin7ChargedMulti->Fill(nMultiTrack);
	fMultiBin7Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=40)&&(lPercentile<50)) {
	fMultiBin8ChargedMulti->Fill(nMultiTrack);
	fMultiBin8Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=50)&&(lPercentile<70)) {
	fMultiBin9ChargedMulti->Fill(nMultiTrack);
	fMultiBin9Vertex->Fill(vertex->GetZ());
      }
      if ((lPercentile>=70)&&(lPercentile<100)) {
	fMultiBin10ChargedMulti->Fill(nMultiTrack);
	fMultiBin10Vertex->Fill(vertex->GetZ());
      }
      //
      for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
	
 	AliESDtrack* track = esd->GetTrack(iTracks);
 	if (!track) {
 	  Printf("ERROR: Could not receive track %d", iTracks);
 	  continue;
 	}

		//    sigmas
	nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
	
	//=======================new 
	//*   test back 2015
	Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
	Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
	if (track->GetTPCNclsF()>0) 
	  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
	
 	Int_t indexKinkPos=track->GetKinkIndex(0);   // kink index 
	
	Int_t tpcNCl = track->GetTPCclusters(0);  
	Double_t tpcSign = track->GetSign();  
	
	Int_t label = track->GetLabel();
	label = TMath::Abs(label);

	if(label > mcEvent->GetNumberOfTracks()) continue;
	if (label > nPrim) continue;
	TParticle * part = stack->Particle(label);
	if (!part) continue;
	
	UInt_t status=track->GetStatus();
	
	if((status&AliESDtrack::kITSrefit)==0) continue;   
	if((status&AliESDtrack::kTPCrefit)==0) continue;
	if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;  
	
	Double_t extCovPos[15];
	track->GetExternalCovariance(extCovPos);    
       	
 	track->GetXYZ(vtrack); 
	
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
 
// //  14/2/13 /================/   
 	 if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;   //    

	 //  track Mult. after selection 
 	 nESDTracK++;        
	 
// //=========================================
 	 fTrackPtAfterTrackCuts->Fill(track->Pt());
	 
// 	 //  select kinks
 	 if(indexKinkPos<0) {     ////mother kink
	   
 	   fPtAllKinks->Fill(track->Pt());  // Pt from tracks , all kinks
	   
 	   fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);
	   
// 	   // select kink class	
	   
 	   AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
	   
// 	   // DCA kink
	   Double_t  Dist2 = kink->GetDistance();
	   
	   Int_t eSDfLabel1=kink->GetLabel(0);
	   Int_t eSDfLabeld=kink->GetLabel(1);
	   
	   TParticle *particle1= stack->Particle(TMath::Abs(eSDfLabel1));
	   Int_t code1= particle1->GetPdgCode();
	   
	   TParticle *particled= stack->Particle(TMath::Abs(eSDfLabeld));
          Int_t dcode1= particled->GetPdgCode();
          Int_t mcProcssDa= particled->GetUniqueID();

// 	  //    loop on MC daugthres for 3Pi    24/9/2010
	  Int_t nESDKpi =0;
	  if(mcProcssDa==4) {
	    Int_t firstDa=particle1->GetFirstDaughter();
	    Int_t lastDa=particle1->GetLastDaughter();
	    
	    if( (lastDa<=0) || (firstDa<=0)  ) continue; 
	    
	    if (lastDa > mcEvent->GetNumberOfTracks()) continue;
	    if (firstDa > mcEvent->GetNumberOfTracks()) continue;
	    //loop on secondaries
	    for (Int_t kk=firstDa;kk<=lastDa;kk++) {
              if ( kk > 0 )    {
		TParticle*    daughter3=stack->Particle(kk);   // 24/9   
		Float_t dcode3= daughter3->GetPdgCode();
		if (((code1==321 )&&(dcode3==211))||((code1 == -321 )&&(dcode3==-211))) nESDKpi++ ; 
             }
	    }
	  }
	  
 	   Double_t hVzdau=particled->Vz();

// 	  //   if (Dist2 > 0.2) continue; //  systematics 11/8/11 
	   
// 	   // TPC mother momentum 
	   
	   const TVector3 vposKink(kink->GetPosition()); //reco position
	   fxyKinkPosition ->Fill( vposKink[0], vposKink[1]);
	   Double_t  dzKink=vpos[2]-vposKink[2]; 
// 	   //
// 	   //   lifitime
	   Double_t tanLamda = track->GetTgl();  // 25/6/2010
	   if (tanLamda ==0 )  continue;
	   
	   Double_t lifeKink= (TMath::Abs(dzKink))*( TMath::Sqrt(1.+ tanLamda*tanLamda) ) / (TMath::Abs( tanLamda)) ;
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
	   
// 	   //          maximum decay angle at a given mother momentum

	   Double_t maxDecAngKmu=f1->Eval(track->P(),0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(track->P(),0.,0.,0.);
	   Float_t signPt= tpcSign*trackPt;

 	   if(( kinkAngle<2.)) fQtBeforeAngleCut->Fill(qT) ;
	   
	   if(((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==321))) fFakeKPi->Fill(track->Pt());
	   if(((TMath::Abs(code1)==211)&&(TMath::Abs(dcode1)==211))) fFakepipi->Fill(track->Pt());
	   
	   //  fake kinks are removed 
	   if( (kinkAngle<2.)) continue;
	             
	   //  BG  ?????==============   ???????Giati einai to parakatw? Den yparxei sto MC all yparxei sta data
	   //  if (TMath::Abs(vposKink[2]) >  225.) continue ;
           //if (TMath::Abs(vposKink[2]) <  0.5) continue ;   XXXXXXXXX???????????
	   
	   fPtKaonInTPC ->Fill(trackPt);     
	   fAngleMomKaonsinTPC->Fill(track->P(), kinkAngle);
	   
	   if ((TMath::Abs(code1)==211)&&(TMath::Abs(dcode1)==13))      
	     fRadiusPtPion->Fill(kink->GetR(), track->Pt()); //
	   
	   if(((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
	      ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))||    
	      ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))) { 
	     fRadiusPtKaon->Fill(kink->GetR(), track->Pt()); //
	   }

	   if((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp )&&(TMath::Abs(rapiditK)<fRapiK)&&(label<nPrim)) {
	     if (((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))) fQtKMu->Fill(qT);
	     if ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11)) fQtKEl->Fill(qT); 
	     if ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211)) fQtKPi->Fill(qT); 
	     if ((nESDKpi>1)&&(((code1)==321)&&((dcode1)==211))) fQtK3PiP->Fill(qT); 
	     if ((nESDKpi>1)&&(((code1)==-321)&&((dcode1)==-211))) fQtK3PiM->Fill(qT); 
	     if (((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
		       ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))||    
		       ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))) {
	       
	       if(qT>fLowQt) fPtKaonPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
	       if(qT>fLowQt) {
		 if(code1>0) fPtKaonPDGPos->Fill(trackPt); // 26/feb  // ALL KAONS (pdg) inside ESD  kink sample
		 if(code1<0) fPtKaonPDGNeg->Fill(trackPt); // 26/feb  // ALL KAONS (pdg) inside ESD  kink sample
	       }
	       fKaonPDGEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
	       fKaonPDGrapidity->Fill(rapiditK) ;  //18/feb rapiddistr of PDG kink ESD  kaons
	               
	       if(qT > fLowQt) fKaonPDGqT->Fill(qT);  // PDG ESD kaons            
	       fKaonPDGpTvsRadius->Fill(kink->GetR(), track->Pt()); //
	     }
	   }

	   if(TMath::Abs(code1)==321) fAngMomKaonPDG->Fill(track->P(), kinkAngle); 
	   if(TMath::Abs(code1)==211) fAngMomPi->Fill(track->P(), kinkAngle); 
	   
	   //
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
	   
           fQtInvMassKaonTPC->Fill(invariantMassKmu, qT);
           fInvMassPtKaonTPC->Fill(invariantMassKmu, trackPt);
	   //
	   fRadiusVsPtKaonTPC->Fill(kink->GetR(), trackPt); //

	   if (qT>fLowQt) fSignPtNcl->Fill(signPt, tpcNCl);

	   if((kink->GetR()>fKinkRadLow)&&(kink->GetR()<fKinkRadUp))  {
	     //  for systematics   if( ( kink->GetR()> 130 ) && ( kink->GetR() < 200 )  )  {
	     if (qT>fLowQt) {
	       fAngleVsMomentumKaonsInR->Fill(track->P(), kinkAngle); 
	       fInvMassKaonInR->Fill(invariantMassKmu);
	       fRadiusVsNclInR->Fill((kink->GetR()), (track->GetTPCclusters(0))) ;
	       if (TMath::Abs(rapiditK)<fRapiK) {
		 fkaonToMu->Fill(invariantMassKmu);
		 fkaonToPi->Fill(invariantMassKpi);
		 fkaonToKa->Fill(invariantMassKK);
		 fRadiusNcl->Fill((kink->GetR()),track->GetTPCclusters(0));
	       }

	     }
	   }
	   
// 	   //  tails cleaning
 	   if((tpcNCl<fLowCluster)) continue;  // test 27 feb 2012 ,, OK

// 	   // cleaning BG in tails
 	   Int_t tpcNClHigh = -31.67+ (11./12.)*(kink->GetR()) ;

	   if (tpcNCl>tpcNClHigh) {
	     fcodeMotherVsDaughter->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
	     fZMotherVsDaughter->Fill( vposKink[2],hVzdau);
	     fPtVsRadiusFake->Fill(kink->GetR(), track->Pt());
	   }	    
	   
	   Int_t tpcNClMin  = -85.5 + (65./95.)*(kink->GetR()) ;
	   if (tpcNCl > tpcNClHigh) continue;  
	   if (tpcNCl < tpcNClMin) continue;   
	   
 	   fKaonKinkPtAfterKinkNclCut->Fill(track->Pt());  // ALL  K-candidates until now                 

// 	   //Final kaon kinks selection 
 	   if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<fRapiK)&&(invariantMassKmu<0.8)) {
	  
	     fAngleVsMomPreSelectedKinks->Fill(track->P(), kinkAngle); 
	     fPtPreSelectedkinks->Fill(trackPt);
	     
// 	     //  maximum angles selection with some error cut
// 	     //Whyyyyyyyyy????????
	     if ((kinkAngle<maxDecAngpimu*1.2)) continue; 
	     if ((kinkAngle>maxDecAngKmu*.98) && (track->P()>1.2)) continue;  ///5/5/2010
	     
 	     fPtKinkBeforedEdx->Fill(trackPt);     
// 	     //  here the kaons selected by the decay features
 	     fTPCSgnlPa->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()));
// 	     //
// 	     //  NO dEdx cut test 9/2/13        
 	     if (nsigma > 3.5) continue; 
// 	     // 
// 	     //  next plots for the identified kaons by the kink analysis
	     
	     fTPCSignalVsMomSelectedKaons->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal())) ;
	     fNclustersVsRadiusSelectedKaons->Fill((kink->GetR()), (track->GetTPCclusters(0))) ;
	     fPtVsRadiusSelectedKaons->Fill(kink->GetR(), trackPt); // 
	     fPtVsInvMassSelectedKaons->Fill(invariantMassKmu, trackPt);

	     Int_t ESDLabelM=0. ;                                      
	     Int_t ESDLabelD= 0. ;                                      
	     Double_t dEdxDauMC=0.0;
	     Int_t Ikink=0;
	     Int_t IRkink=0;
	     
	     for (Int_t jTrack = 0; jTrack < esd->GetNumberOfTracks(); jTrack++) {
	       
	       AliESDtrack* trackDaughter = esd->GetTrack(jTrack);
	       if (!trackDaughter) {
	     	 Printf("ERROR: Could not receive track %d", jTrack);
	     	 continue;
	       }
	       Int_t indexKinkDAU =trackDaughter->GetKinkIndex(0);
                     if (indexKinkDAU <0) {
	     	       AliESDkink *kinkDau=esd->GetKink(TMath::Abs(indexKinkDAU)-1);
	     	       //            raDAU= kinkDau->GetR();
	     	       ESDLabelM=kinkDau->GetLabel(0);   //  mothers's label
	     	       ESDLabelM = TMath::Abs(ESDLabelM);
	     	       ESDLabelD=kinkDau->GetLabel(1);   //  Daughter'slabel
	     	       ESDLabelD = TMath::Abs(ESDLabelD);
	     	       if (kink->GetR() == kinkDau->GetR()) IRkink++;
	     	       if (ESDLabelM == label) Ikink++  ;
	     	     }

	     	     if ((Ikink >0)&&(IRkink>0)) { 
	     	       // daughter kink 

	     	       if((indexKinkDAU>0)) dEdxDauMC = trackDaughter->GetTPCsignal();  //  daughter kink 
	     	     }
	     }

	     fPtKaonREC->Fill(track->Pt());
	     if( code1>0.) fPtKaonRECpos->Fill(track->Pt());   //all PID kink-kaon
	     if( code1<0.) fPtKaonRECneg->Fill(track->Pt());   //all PID kink-kaon
	     
	     fdedxMthVsTPCMomSelectedKaons->Fill(trMomTPC, track->GetTPCsignal()) ;               //trMomTPC
	     fSignalMthVsSignalDaughterSelectedKaons->Fill(dEdxDauMC, track->GetTPCsignal()) ;
	     fSignalDaughterVsDaughterMomSelectedKaons->Fill(daughterMKink.Mag(), dEdxDauMC) ;  //  daughter kink 
// 	     // 
 	     fNsigmaVsTPCmomSelectedKaons->Fill(trMomTPC, nsigmall);     

	     fRadiusSelectedKinks->Fill(kink->GetR());  // kink 
	     fKaonLifetimeSelectedKaons->Fill(lifeKink);  // kink 
             fTrackEtaSelectedKaons->Fill(trackEta);
	     fRapiditySelectedKaons->Fill(rapiditK);  //  rapidityof kaons 
	     fZKinkProductionVsKinkRadSelectedKaons->Fill(vposKink[2], kink->GetR());
	     
	     Float_t signPt= tpcSign*trackPt;
	     fNclinTPCVsSignedPtSelectedKaons->Fill(signPt, tpcNCl);   ///  28/4/2010
	     fRapidityVsSignedPtSelectedKons->Fill(signPt, rapiditK);
	     fNclVsRapiditySelectedKAons->Fill(rapiditK, tpcNCl);
	     fNclVsChi2SelectedKaons->Fill((track->GetTPCchi2()), tpcNCl);
	     fChi2perTPCclusterSelectedKaons-> Fill((track->GetTPCchi2()/track->GetTPCclusters(0))) ;

	     fLifetimeSelectedKaons->Fill((lifeKink*.493667)/track->P());
             fPtSelectedKaons->Fill(track->Pt());        
	     fDCAkinkSelectedKaons->Fill(Dist2);
	     
	     if(tpcSign >0.) fPtPositiveSelectedKaons->Fill(track->Pt());   //K-plus bins Comb 
	     if(tpcSign <0.) fPtNegativeSelectedKaons->Fill(track->Pt());   //K-minus bins Comb
	     

// 	     //  kaons from k to mun and k to pipi and to e decay 
	     if(((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
		 ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))||    
		 ((TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))) {

	        fkinkKaonPDG->Fill(track->Pt());
		if(code1>0.) fkinkKaonPDGpos->Fill(trackPt);  //                  kPtPID kink-kaon
		if(code1<0.) fkinkKaonPDGneg->Fill(trackPt);    //

			     // Multiplicity bins
	     
	     if ((lPercentile>0)&&(lPercentile<1)) {
	       if (code1 >0.) fkinkKaonPDGposMulti1->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti1->Fill(track->Pt());
	     }
	     if ((lPercentile>=1)&&(lPercentile<5)) {
	       if (code1 >0.) fkinkKaonPDGposMulti2->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti2->Fill(track->Pt());
	     }
	     if ((lPercentile>=5)&&(lPercentile<10)) {
	       if (code1 >0.) fkinkKaonPDGposMulti3->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti3->Fill(track->Pt());
	     }
	     if ((lPercentile>=10)&&(lPercentile<15)) {
	       if (code1 >0.) fkinkKaonPDGposMulti4->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti4->Fill(track->Pt());
	     }
	     if ((lPercentile>=15)&&(lPercentile<20)) {
	       if (code1 >0.) fkinkKaonPDGposMulti5->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti5->Fill(track->Pt());
	     }
	     if ((lPercentile>=20)&&(lPercentile<30)) {
	       if (code1 >0.) fkinkKaonPDGposMulti6->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti6->Fill(track->Pt());
	     }
	     if ((lPercentile>=30)&&(lPercentile<40)) {
	       if (code1 >0.) fkinkKaonPDGposMulti7->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti7->Fill(track->Pt());
	     }
	     if ((lPercentile>=40)&&(lPercentile<50)) {
	       if (code1 >0.) fkinkKaonPDGposMulti8->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti8->Fill(track->Pt());
	     }
	     if ((lPercentile>=50)&&(lPercentile<70)) {
	       if (code1 >0.) fkinkKaonPDGposMulti9->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti9->Fill(track->Pt());
	     }
	     if ((lPercentile>=70)&&(lPercentile<100)) {
	       if (code1 >0.) fkinkKaonPDGposMulti10->Fill(track->Pt());
	       if (code1 <0.) fkinkKaonPDGnegMulti10->Fill(track->Pt());
	     }
		
	     }
	     else {
	       fkinkKaonPDGBkg->Fill(trackPt);

	       	      // Multiplicity bins
	      
	      if ((lPercentile>0)&&(lPercentile<1)) {
		fkinkKaonPDGBkgMulti1->Fill(track->Pt());
	      }
	      if ((lPercentile>=1)&&(lPercentile<5)) {
		fkinkKaonPDGBkgMulti2->Fill(track->Pt());
	      }
	      if ((lPercentile>=5)&&(lPercentile<10)) {
		fkinkKaonPDGBkgMulti3->Fill(track->Pt());
	      }
	      if ((lPercentile>=10)&&(lPercentile<15)) {
		fkinkKaonPDGBkgMulti4->Fill(track->Pt());
	      }
	      if ((lPercentile>=15)&&(lPercentile<20)) {
		fkinkKaonPDGBkgMulti5->Fill(track->Pt());
	      }
	      if ((lPercentile>=20)&&(lPercentile<30)) {
		fkinkKaonPDGBkgMulti6->Fill(track->Pt());
	      }
	      if ((lPercentile>=30)&&(lPercentile<40)) {
		fkinkKaonPDGBkgMulti7->Fill(track->Pt());
	      }
	      if ((lPercentile>=40)&&(lPercentile<50)) {
		fkinkKaonPDGBkgMulti8->Fill(track->Pt());
	      }
	      if ((lPercentile>=50)&&(lPercentile<70)) {
		fkinkKaonPDGBkgMulti9->Fill(track->Pt());
	      }
	      if ((lPercentile>=70)&&(lPercentile<100)) {
		fkinkKaonPDGBkgMulti10->Fill(track->Pt());
	      }
	     }

// 	     // Multiplicity bins
	     
	     if ((lPercentile>0)&&(lPercentile<1)) {
	       if (tpcSign >0.) fMultiBin1KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin1KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=1)&&(lPercentile<5)) {
	       if (tpcSign >0.) fMultiBin2KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin2KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=5)&&(lPercentile<10)) {
	       if (tpcSign >0.) fMultiBin3KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin3KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=10)&&(lPercentile<15)) {
	       if (tpcSign >0.) fMultiBin4KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin4KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=15)&&(lPercentile<20)) {
	       if (tpcSign >0.) fMultiBin5KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin5KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=20)&&(lPercentile<30)) {
	       if (tpcSign >0.) fMultiBin6KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin6KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=30)&&(lPercentile<40)) {
	       if (tpcSign >0.) fMultiBin7KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin7KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=40)&&(lPercentile<50)) {
	       if (tpcSign >0.) fMultiBin8KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin8KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=50)&&(lPercentile<70)) {
	       if (tpcSign >0.) fMultiBin9KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin9KaonKinksNeg->Fill(track->Pt());
	     }
	     if ((lPercentile>=70)&&(lPercentile<100)) {
	       if (tpcSign >0.) fMultiBin10KaonKinksPos->Fill(track->Pt());
	       if (tpcSign <0.) fMultiBin10KaonKinksNeg->Fill(track->Pt());
	     }

// 	     // End Multiplicity bins
	     
	     
// 	     //	     fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK); 
 	     fAngMomK->Fill(track->P(), kinkAngle); 
	     fPosiKinkK->Fill(vposKink[0], vposKink[1]);
	     fPosiKinKXZ->Fill(vposKink[2], vposKink[0]);
	     fPosiKinKYZ->Fill(vposKink[2], vposKink[1]);
	     
 	   }  //  kink selection 
	   

 	 }  //End Kink Information    
  

       } //track loop 
   
   PostData(1, fListOfHistos);
 
}      

//________________________________________________________________________
void AliAnalysisKinkTaskMult13ppMC::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}

//____________________________________________________________________//

  // Get the vertex 
const AliESDVertex* AliAnalysisKinkTaskMult13ppMC::GetEventVertex(const AliESDEvent* esd) const
  
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
Bool_t AliAnalysisKinkTaskMult13ppMC::selectVertex2015pp(AliESDEvent *esd,
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
                       const  AliESDVertex * vertex = esd->GetPrimaryVertex();
           if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
              return kTRUE;
                  }

             //___________________________________________________________________________________ 
//___________________________________________________________________________________

            //Bool_t AliAnalysisTaskExtractV0me::IsGoodSPDvertexRes(const AliESDVertex  * spdVertex)
            Bool_t AliAnalysisKinkTaskMult13ppMC::IsGoodSPDvertexRes(const AliESDVertex  * spdVertex)
              {
            if (!spdVertex) return kFALSE;
                if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04
          && spdVertex->GetZRes()<0.25)) return kFALSE;
                              return kTRUE;
                                  }


