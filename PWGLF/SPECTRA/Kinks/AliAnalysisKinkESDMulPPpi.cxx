/*************************************************************************
 *  Authors: Martha Spyropoulou-Stassinaki and the  members               * 
 * of the Greek group at Physics Department of Athens University          *
 * Paraskevi Ganoti, Anastasia Belogianni and Filimon Roukoutakis.        *
 * The method is applied in pp and Pb-Pb real data.                       *
 *        since 2018  is applied on Xe-Xe data , Kaons and Pions          *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDat13 class
//       Example of an analysis task for kink topology study
//      pions from kink topology are 'identified' in this code
//    March2020 : Nominal R->110-205 cm,  Rapidity Pion 0.5, eta< 0.8, back to 0.5 for XeXe 14/4/18
//    27/3/20  sigmaPion  -4 ----1.5  ?      dcaXY < 0.3 and dcaZ < 2.5cm
//    28/3/20  sigmaPion -3.5 ---3.5    NO cut for angle<1 degree , dcaXY < 0.3 and dcaZ < 2.5cm  the same, plots for counting 
//          diorthvsi in Max angle for pions !!!!!
//          1/8/20   sigmaPi <3,  No Angle cut
//          2/3/21   pions from 13 TeV pp  LHC18,  MB 
//          7/3/21 mult range , cent 0-100
//          10/3/21 mult bins 14
//          31/3/ XeXe ->  PPpi
//            14/4   Bins2 and 3 limits 0-.001,  0.001-0.005
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
#include "AliAnalysisKinkESDMulPPpi.h"
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
#include "AliCentrality.h"
#include "iostream"
#include "AliMultiplicity.h"
#include "AliMultSelection.h"
/////////#include "AliTENDERSupplies.h"
ClassImp(AliAnalysisKinkESDMulPPpi)


//________________________________________________________________________
AliAnalysisKinkESDMulPPpi::AliAnalysisKinkESDMulPPpi(const char *name) 
  : AliAnalysisTaskSE(name), fMultiplMC(0),fIncompletEv(0),  fMultMCK(0),   fMultESDK(0),fZpr(0),    fESDMult(0), 
   fZMainVx(0),  fMultVertex(0),                         
   fHistPtESD(0),fHistPt(0),fHistQtAll(0),fHistQt1(0),fHistQt2(0)
  , fHistPtKaon(0),fHistPtKPDG(0),fHistEta(0),fHistEtaK(0),fptKMC(0),fgenpt(0),frad(0),
  fKinkKaon(0),fKinKRbn(0), fKinkKaonBg(0), fM1kaon(0),  fPtKink(0),  fptKink(0),
    fAngMomK(0),fAngMomPi(0), fAngMomKC(0),   
 fSignPtNcl(0), fSignPtEta(0), fEtaNcl(0), fSignPt(0), fChi2NclTPC(0), fRatChi2Ncl(0), fRadiusNcl(0), fTPCSgnlP(0),
   fTPCSgnlPa(0), fRpr(0), fdcatoVxXY(0), fnSigmToVx(0), fKinkMothDau(0),
 fZvXv(0),fZvYv(0), fXvYv(0),  fHistPtKaoP(0), fHistPtKaoN(0),frapiKESD(0), flifetime(), fradLK(0),
    fradPtRpDt(0), fInvMuNuAll(0), fQtInvM(0), 
         fDCAkink(0), fPosiKink(0),  fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0), fPosiKinKBg(0), fQtMothP(0), fTPCSgnlPtpc(0),
       fTPCMomNSgnl(0),  fMothKinkMomSgnl(0), fNSigmTPC(0),  fTPCSgnlKinkDau(0), fPtKinkPos(0), fPtKinkNeg(0),  fRadNclCln(0),
       fRatioCrossedRows(0), fRatioCrossedRowsKink(0), 
      fRadiusPt(0), fRadiusPtcln(0),  fInvMassMuNuPt(0), fInvMassMuNuPtAll(0),fPtCut1(0), fPtCut2(0), 
      fPtCut3(0), fAngMomKKinks(0),fPtKinkK0(0), fPtKinkK0P(0), fPtKinkK0N(0), fPtKinkGyu(0), fPtKinkGyuP(0), fPtKinkGyuN(0), 
 f1(0), f2(0),
  //    fListOfHistos(0), fPIDResponse()
 //   fListOfHistos(0),fLowMulcut(-1),fUpMulcut(-1), fKinkRadUp(210.), fKinkRadLow(120.), fLowCluster(20), fLowQt(.12), fRapiK(0.5),  fCutsMul(0),   fMaxDCAtoVtxCut(0), fPIDResponse()
    // mexri 15/2/18 fListOfHistos(0), fKinkRadUp(210.), fKinkRadLow(120.), fLowCluster(20), fLowQt(.12), fRapiK(0.5),   fPIDResponse(0)
    //  allagh 14/4/18fListOfHistos(0), fKinkRadUp(205.), fKinkRadLow(110.), fLowCluster(20), fLowQt(.12), fRapiK(0.7),   fPIDResponse(0)
    fListOfHistos(0), fKinkRadUp(205.), fKinkRadLow(110.), fLowCluster(20), fLowQt(.04), fRapiK(0.5),   fPIDResponse(0)
 , fNumberOfEvent(0),  fbgCleaningHigh(0), fMultiplicity(0), fEventVertx(0),  fIncompletEvent(0), fMultpileup(0) ,  fcentrality(0), fcentral0To1(0), fcentral0To01(0)
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
  , fMultiBin11ChargedMulti(0)
  , fMultiBin12ChargedMulti(0)
  , fMultiBin13ChargedMulti(0)
  , fMultiBin14ChargedMulti(0)
, fMultiBin15ChargedMulti(0)
  , fMultiBin16ChargedMulti(0)
  , fMultiBin17ChargedMulti(0)
  , fMultiBin18ChargedMulti(0)
  , fMultiBin19ChargedMulti(0)
  , fMultiBin20ChargedMulti(0)
  , fMultiBin21ChargedMulti(0)
  , fMultiBin22ChargedMulti(0)
  , fMultiBin23ChargedMulti(0)
//
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
  , fMultiBin11KaonKinksPos(0)
  , fMultiBin12KaonKinksPos(0)
  , fMultiBin13KaonKinksPos(0)
  , fMultiBin14KaonKinksPos(0)
 , fMultiBin15KaonKinksPos(0)
  , fMultiBin16KaonKinksPos(0)
  , fMultiBin17KaonKinksPos(0)
  , fMultiBin18KaonKinksPos(0)
  , fMultiBin19KaonKinksPos(0)
  , fMultiBin20KaonKinksPos(0)
  , fMultiBin21KaonKinksPos(0)
  , fMultiBin22KaonKinksPos(0)
  , fMultiBin23KaonKinksPos(0)
//
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
  , fMultiBin11KaonKinksNeg(0)
  , fMultiBin12KaonKinksNeg(0)
  , fMultiBin13KaonKinksNeg(0)
  , fMultiBin14KaonKinksNeg(0)
, fMultiBin15KaonKinksNeg(0)
  , fMultiBin16KaonKinksNeg(0)
  , fMultiBin17KaonKinksNeg(0)
  , fMultiBin18KaonKinksNeg(0)
  , fMultiBin19KaonKinksNeg(0)
  , fMultiBin20KaonKinksNeg(0)
  , fMultiBin21KaonKinksNeg(0)
  , fMultiBin22KaonKinksNeg(0)
  , fMultiBin23KaonKinksNeg(0)

//
  , fpercentMul(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
//  DefineOutput(1, TList::Class());
// }
 //   apo Evi
/*{

        // Constructor  
  fClusterZvtx=new TArrayD(1500);
  fNEvClusters=new TArrayD(1500);
  DefineInput (0, TChain::Class());
  DefineOutput(1, TList::Class());
}
*/
/*
//----------------------Marek multiplicity bins 
 fCutsMul=new AliESDtrackCuts("Mul","Mul");
        fCutsMul->SetMinNClustersTPC(70);
        fCutsMul->SetMaxChi2PerClusterTPC(4);
        fCutsMul->SetAcceptKinkDaughters(kFALSE);
        fCutsMul->SetRequireTPCRefit(kTRUE);
        // ITS
        fCutsMul->SetRequireITSRefit(kTRUE);
        fCutsMul->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                                AliESDtrackCuts::kAny);
        fCutsMul->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

        fCutsMul->SetMaxDCAToVertexZ(2);
        fCutsMul->SetDCAToVertex2D(kFALSE);
        fCutsMul->SetRequireSigmaToVertex(kFALSE);

        fCutsMul->SetEtaRange(-0.8,+0.8);
        fCutsMul->SetPtRange(0.15, 1e10);

        fMaxDCAtoVtxCut=new AliESDtrackCuts("fMaxDCA", "fMaxDCA");
       fMaxDCAtoVtxCut->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fMaxDCAtoVtxCut->SetMaxChi2TPCConstrainedGlobal(36);
//        esdTrackCuts->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
  //  esdTrackCuts->SetMaxChi2TPCConstrainedGlobal(36);

*/

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkESDMulPPpi::UserCreateOutputObjects() 
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
   //Open file  1= CAF 
    //OpenFile(1); 
/*
   Double_t gPt7K0[45] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.2, 2.4, 2.6, 2.8,  3.0,   3.3, 3.6, 3.9,   
                         4.2, 4.6,5.0, 5.4, 5.9,  6.5,   7.0,7.5, 8.0,8.5,  9.2, 10., 11., 12., 13.5,15.0 };  // David K0
*/
 Double_t gPt7Comb[48] = {
0.25,0.30,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0, 1.1, 1.2,
1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0, 5.5, 6
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
  fMultiplMC= new TH1F("fMultiplMC", "charge multiplicity ESD",400,10.0,2010.); 
  fIncompletEv=new TH1F("fIncompletEv", "charge multiplicity ESD after Incom. cut",400,10.0,2010.); 
  fMultMCK=new TH1F("fMultMCK", "charge multiplity  after Trigger",400, 0.0,2000.); 
  fMultESDK=new TH1F("fMultESDK", "charge multipliESD after Pileup",400,10.0,2010.); 
  fZpr = new TH1D("fZpr", "z distribution of main vertex",60,-15.,15.);
  fESDMult= new TH1F("fESDMult", "ESD charge mult. inside 10 cm Zvx",400, 0.0,2000.); 
  fZMainVx= new TH1D("fZMainVx", "ESD charge mult. Main Vertex", 60,-15.,15.); 
  fMultVertex= new TH1F("fMultVertex", "ESD charge mult. inside Main Vertx",400, 0.0,2000.); 
//
  fHistPtESD = new TH1F("fHistPtESD", "P_{T} distribution",300, 0.0,15.0);
  fHistPtESD->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtESD->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtESD->SetMarkerStyle(kFullCircle);
  fHistPt = new TH1F("fHistPt", "P_{T} distribution",300, 0.0,15.0); 
  fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300); 
  fHistQt1= new TH1F("fHistQt1", "Q_{T} distribution",100, 0.0,.300); 
  fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300); 
  fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",300, 0.0,15.0); 
  fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",300, 0.0,15.0); 
  fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3); 
  fHistEtaK= new TH1F("fHistEtaK", "EtaK distribution", 26,-1.3, 1.3); 
  fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",300, 0.0,15.0); 
  fgenpt= new TH1F("fgenpt", "genpt   K distribution",300, 0.0,15.0); 
   //frad= new TH1F("frad", "radius  K generated",100, 50., 250.0);
   frad= new TH1F("frad", "radius  K generated",100, 0.,1000.0);
  fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi",300, 0.0,15.0); 
  fKinKRbn= new TH1F("fKinKRbn", "p_{t}Kaon kinks identi[GeV/c],Entries",46,gPt7TOF); 
  fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",300, 0.0,15.0); 
  //fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",180,0.10, 1.0); 
  fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",600,0.10, 0.7); //  23/8/2013
  //fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  distribution, counts",44, gPt7K0); 
  fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  distribution, counts",47, gPt7Comb); 
  fptKink= new TH1F("fptKink", "P_{T}Kaon Kink  bution",300, 0.0,15.0); 
  fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,5.0,80,0.,80.);
  fAngMomKC= new TH2F("fAngMomKC","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",80,-4.,4.0,70,20.,160.);
  fSignPtEta= new TH2F("fSignPtEta","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
  fEtaNcl= new TH2F("fEtaNcl","Eta vrs nclust,K",30,-1.5,1.5, 70,20, 160);
  fSignPt= new TH1F("fSignPt","SignPt ,K",80,-4.0,4.0);
  fChi2NclTPC= new TH2F("fChi2NclTPC","Chi2vrs nclust,K",100,0.,500., 70,20, 160);
  fRatChi2Ncl= new TH1F("fRatChi2Ncl","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
  fRadiusNcl = new TH2F("fRadiusNcl","kink radius vrs Nclust,K",75,100.,250., 80,0, 160);
    fTPCSgnlP = new TH2F("fTPCSgnlP","TPC signal de/dx Mom,K",1000,0.0,20.0,150,0.,300.);
  fTPCSgnlPa= new TH2F("fTPCSgnlPa","TPC signal de/dx Mom,K",1000,0.0,20.,150, 0.,300.);
  fRpr = new TH1D("fRpr", "rad distribution  PID pr",100,-10.0, 10.0);
  fdcatoVxXY = new TH1D("fdcatoVxXY", "dca  distribution PID  ",40,-1.,1.);
  fnSigmToVx = new TH1D("fnSigmToVx", "dca  distribution PID  ",80,0.,8.);
  fKinkMothDau= new TH2F("fKinkMothDau","TPC kink Moth Daugh ,K",50,0.0,2.5,50, 0., 2.5);
  fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5,60, -15., 15.0);
  fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
  fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
  fHistPtKaoP = new TH1F("fHistPtKaoP", "P_{T}KaonP  distribution",300, 0.0,15.0); 
  fHistPtKaoN = new TH1F("fHistPtKaoN", "P_{T}KaonN  distribution",300, 0.0,15.0); 
  frapiKESD=new TH1F("frapiKESD","rapid Kdistribution", 26,-1.3, 1.3); 
  flifetime= new TH1F("flifetime", "ct study of K-kinks",100,0.,1000.); 
  fradLK= new TH1F("fradLK", "Length of   K generated",100,0.,1000.); 
  fradPtRpDt=new TH3F("fradPtRpDt","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
  //fInvMuNuAll= new TH1F("fInvMuNuAll", " Inv Mass MuNu all kink",180,0.1,1.0); 
  fInvMuNuAll= new TH1F("fInvMuNuAll", " Inv Mass MuNu all kink",600,0.1,0.7); //  23/8/2013
  fQtInvM= new TH2F("fQtInvM", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300); 
    fDCAkink = new TH1F("fDCAkink ", "DCA kink vetrex ",40,-1.0,1.0);
  fPosiKink= new TH2F("fPosiKink", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKBg= new TH2F("fPosiKinKBg", "z vrx kink rad    ",100, -300.0,300.0,100,100., 300.);
  fQtMothP = new TH2F("fQtMothP", " Qt vrs Mother P", 100, 0., 5.0,100, 0.,0.300);
    fTPCSgnlPtpc = new TH2F("fTPCSgnlPtpc","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,100, 0., 250.    );
    fTPCMomNSgnl = new TH2F("fTPCMomNsgnl","TPC signal de/dx Mom TPC,K  ",300,0.0,15.0,20 , -10., 10.);
    fMothKinkMomSgnl  = new TH2F("fMothKinkMomSgnl","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.0,100, 0., 250.    );
    fNSigmTPC    = new TH1F("fNSigmTPC","TPC Nsigma  de/dx  TPC,K  ", 30 , -7.5, 7.5);
    fTPCSgnlKinkDau = new TH2F("fTPCSgnlKinkDau","TPC signal de/dx Mom,K",500,0.0,10.0,100,0.,250.);
  //fPtKinkPos= new TH1F("fPtKinkPos", "Pos P_{T}Kaon Kink  distribution, counts",44, gPt7K0); 
  fPtKinkPos= new TH1F("fPtKinkPos", "Pos P_{T}Kaon Kink  distribution, counts",47, gPt7Comb ); 
  fPtKinkNeg= new TH1F("fPtKinkNeg", "Neg P_{T}Kaon Kink  distribution, counts",47, gPt7Comb ); 
  fRadNclCln = new TH2F("fRadNclCln","kink radius vrs Nclust,K Clean ",75,100.,250., 80,0, 160);
  fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);
  fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);
  fRadiusPt =new TH2F("fRadiusPt","radius vs pt  ",80, 90.,250.,100, 0.,10.              );
  fRadiusPtcln =new TH2F("fRadiusPtcln","radius vs pt clean ",80, 90.,250.,100, 0.,10.              );
  //fInvMassMuNuPt =new TH2F("fInvMassMuNuPt","Invariant mass-munu  vs pt  ",180, 0.10, 1.00, 100, 0.0, 10.0  );
  fInvMassMuNuPt =new TH2F("fInvMassMuNuPt","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 100, 0.0, 10.0  );// 23/8/2013
  fInvMassMuNuPtAll =new TH2F("fInvMassMuNuPtAll","Invariant mass-munu  vs pt  ",600, 0.10, 0.7, 100, 0.0, 10.0  );// 23/8/2013
  fPtCut1 = new TH1F("fPtCut1", "P_{T}Kaon distribution",300, 0.0,15.0); 
  fPtCut2 = new TH1F("fPtCut2", "P_{T}Kaon distribution",300, 0.0,15.0); 
  fPtCut3 = new TH1F("fPtCut3", "P_{T}Kaon distribution",300, 0.0,15.0); 
  fAngMomKKinks = new TH2F("fAngMomKKinks","Decay angle vrs Mother Mom,Kinks",300,0.0,15.0,100,0.,100.);
  fPtKinkK0= new TH1F("fPtKinkK0", "P_{T}Kaon Kink  distribution, counts",44, gPt13K0PKal); 
  fPtKinkK0P= new TH1F("fPtKinkK0P", "P_{T} KPl Kink  distribution, counts",44, gPt13K0PKal); 
  fPtKinkK0N= new TH1F("fPtKinkK0N", "P_{T} KMn Kink  distribution, counts",44, gPt13K0PKal); 
  fPtKinkGyu= new TH1F("fPtKinkGyu", "P_{T}Kaon Kink  distribution, counts",68, gPt13HPtGyu); 
  fPtKinkGyuP= new TH1F("fPtKinkGyuP", "P_{T} KPl Kink  distribution, counts",68, gPt13HPtGyu); 
  fPtKinkGyuN= new TH1F("fPtKinkGyuN", "P_{T} KMn Kink  distribution, counts",68, gPt13HPtGyu); 

  fNumberOfEvent =  new TH1F("fNumberOfEvent", "Number of Events in this job", 20,  0.,10.); 
//  fNumberOfEvent_cent  =  new TH1F("fNumberOfEvent_cent", "Number of Events in this job", 10,  0.,6.); 
  fbgCleaningHigh      =  new TH1F("fbgCleaningHigh", "  BG  cleaning histo 1      ",20,0.,10.); 
  fMultiplicity= new TH1F("fMultiplicity", "charge multiplicity ESD",20, 0.,10.0); 
  fEventVertx   =  new TH1F("fEventVertx", "Number of Events  after the vertex cut", 20,  0.,10.); 
  fIncompletEvent=new TH1F("fIncompletEvent", "charge multiplicity ESD after Incomplete cut",20,0.,10.0 ); 
  fMultpileup    =  new TH1F("fMultpileup", "Number of Events  after pileup rejection",20,0.,10.); 
//
fcentrality    =  new TH1F("fcentrality", "centrality distribution  from 0 to 100 %",1000,0.,100. );
  fcentral0To1   =  new TH1F("fcentral0To1", "central0To1 distribution  for HM trigger ",1000,0.,  1. );
  fcentral0To01   =  new TH1F("fcentral0To01", "central0To01 distribution  for HM trigger ",1000,0.,  1. );
//
//
   fMultiBin1ChargedMulti= new TH1F("fMultiBin1ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin2ChargedMulti= new TH1F("fMultiBin2ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin3ChargedMulti= new TH1F("fMultiBin3ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin4ChargedMulti= new TH1F("fMultiBin4ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin5ChargedMulti= new TH1F("fMultiBin5ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin6ChargedMulti= new TH1F("fMultiBin6ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin7ChargedMulti= new TH1F("fMultiBin7ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin8ChargedMulti= new TH1F("fMultiBin8ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin9ChargedMulti= new TH1F("fMultiBin9ChargedMulti", "charge multiplicity ESD",400, 0.,2000.0);
   fMultiBin10ChargedMulti= new TH1F("fMultiBin10ChargedMulti", "charge multiplicity ESD",400,0.0,2000.0);
   fMultiBin11ChargedMulti= new TH1F("fMultiBin11ChargedMulti", "charge multiplicity ESD",400,0.0,2000.0);
   fMultiBin12ChargedMulti= new TH1F("fMultiBin12ChargedMulti", "charge multiplicity ESD",400,0.0,2000.0);
   fMultiBin13ChargedMulti= new TH1F("fMultiBin13ChargedMulti", "charge multiplicity ESD",400,0.0,2000.0);
   fMultiBin14ChargedMulti= new TH1F("fMultiBin14ChargedMulti", "charge multiplicity ESD",400,0.0,2000.0);
fMultiBin15ChargedMulti= new TH1F("fMultiBin15ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin16ChargedMulti= new TH1F("fMultiBin16ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin17ChargedMulti= new TH1F("fMultiBin17ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin18ChargedMulti= new TH1F("fMultiBin18ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin19ChargedMulti= new TH1F("fMultiBin19ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin20ChargedMulti= new TH1F("fMultiBin20ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin21ChargedMulti= new TH1F("fMultiBin21ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin22ChargedMulti= new TH1F("fMultiBin22ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);
   fMultiBin23ChargedMulti= new TH1F("fMultiBin23ChargedMulti", "charge multiplicity ESD",400, 0.0,2000.);


   fMultiBin1KaonKinksPos = new TH1F("fMultiBin1KaonKinksPos", "P_{T}Kaon distr 0-.001",1500, 0.0,15.0);
   fMultiBin2KaonKinksPos = new TH1F("fMultiBin2KaonKinksPos", "P_{T}Kaon distr.001-.005   ",1500, 0.0,15.0);
   fMultiBin3KaonKinksPos = new TH1F("fMultiBin3KaonKinksPos", "P_{T}Kaon distr 0.-.005 ",1500, 0.0,15.0);
   fMultiBin4KaonKinksPos = new TH1F("fMultiBin4KaonKinksPos", "P_{T}Kaon distr .005-0.01% ",1500, 0.0,15.0);
   fMultiBin5KaonKinksPos = new TH1F("fMultiBin5KaonKinksPos", "P_{T}Kaon distr  0.0- 0.01 ",1500, 0.0,15.0);
   fMultiBin6KaonKinksPos = new TH1F("fMultiBin6KaonKinksPos", "P_{T}Kaon distr 0.01 -0.05 ",1500, 0.0,15.0);
   fMultiBin7KaonKinksPos = new TH1F("fMultiBin7KaonKinksPos", "P_{T}Kaon distr 0.05 -0.1 ",1500, 0.0,15.0);
   fMultiBin8KaonKinksPos = new TH1F("fMultiBin8KaonKinksPos", "P_{T}Kaon distr 0.1 - 0.2 ",1500, 0.0,15.0);
   fMultiBin9KaonKinksPos = new TH1F("fMultiBin9KaonKinksPos", "P_{T}Kaon distr 0.2 - 0.5 ",1500, 0.0,15.0);
   fMultiBin10KaonKinksPos = new TH1F("fMultiBin10KaonKinksPos", "P_{T}Kaon distr 0.5-1.",1500, 0.0,15.0);
   fMultiBin11KaonKinksPos = new TH1F("fMultiBin11KaonKinksPos", "P_{T}Kaon distr 0 -1. ",1500, 0.0,15.0);
   fMultiBin12KaonKinksPos = new TH1F("fMultiBin12KaonKinksPos", "P_{T}Kaon distr 1 -5. ",1500, 0.0,15.0);
   fMultiBin13KaonKinksPos = new TH1F("fMultiBin13KaonKinksPos", "P_{T}Kaon distr 5.-10.",1500, 0.0,15.0);
   fMultiBin14KaonKinksPos = new TH1F("fMultiBin14KaonKinksPos", "P_{T}Kaon distr 10-15.",1500, 0.0,15.0);
fMultiBin15KaonKinksPos = new TH1F("fMultiBin15KaonKinksPos", "P_{T}Kaon distr 15.- 20.",1500, 0.0,15.0);
   fMultiBin16KaonKinksPos = new TH1F("fMultiBin16KaonKinksPos", "P_{T}Kaon distr 20.- 30.",1500, 0.0,15.0);
   fMultiBin17KaonKinksPos = new TH1F("fMultiBin17KaonKinksPos", "P_{T}Kaon distr 30. -40.",1500, 0.0,15.0);
   fMultiBin18KaonKinksPos = new TH1F("fMultiBin18KaonKinksPos", "P_{T}Kaon distr 40. -50.",1500, 0.0,15.0);
   fMultiBin19KaonKinksPos = new TH1F("fMultiBin19KaonKinksPos", "P_{T}Kaon distr 50. -60.",1500, 0.0,15.0);
   fMultiBin20KaonKinksPos = new TH1F("fMultiBin20KaonKinksPos", "P_{T}Kaon distr 60. -70.",1500, 0.0,15.0);
   fMultiBin21KaonKinksPos = new TH1F("fMultiBin21KaonKinksPos", "P_{T}Kaon distr 70.-80. ",1500, 0.0,15.0);
   fMultiBin22KaonKinksPos = new TH1F("fMultiBin22KaonKinksPos", "P_{T}Kaon distr 80. -90.",1500, 0.0,15.0);
   fMultiBin23KaonKinksPos = new TH1F("fMultiBin23KaonKinksPos", "P_{T}Kaon distr 90-100. ",1500, 0.0,15.0);

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
   fMultiBin11KaonKinksNeg = new TH1F("fMultiBin11KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin12KaonKinksNeg = new TH1F("fMultiBin12KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin13KaonKinksNeg = new TH1F("fMultiBin13KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin14KaonKinksNeg = new TH1F("fMultiBin14KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);

fMultiBin15KaonKinksNeg = new TH1F("fMultiBin15KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin16KaonKinksNeg = new TH1F("fMultiBin16KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin17KaonKinksNeg = new TH1F("fMultiBin17KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin18KaonKinksNeg = new TH1F("fMultiBin18KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin19KaonKinksNeg = new TH1F("fMultiBin19KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin20KaonKinksNeg = new TH1F("fMultiBin20KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin21KaonKinksNeg = new TH1F("fMultiBin21KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin22KaonKinksNeg = new TH1F("fMultiBin22KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);
   fMultiBin23KaonKinksNeg = new TH1F("fMultiBin23KaonKinksNeg", "P_{T}Kaon distribution",1500, 0.0,15.0);



  fpercentMul= new TH1F("fpercentMul", "V0M  percentile    distrib.",220, 0.0 ,110.); 
//
//
//  fPtKinkHigh= new TH1F("fPtKinkHigh", "Sum  P_{T}Kaon Kink  distrib.in High Pt Binning",47, gPt13High); 
//  fPtKinkPlH= new TH1F("fPtKinkPlH", "Pos P_{T}Kaon Kink  distrib.in High Pt Binning",47, gPt13High); 
//  fPtKinkMnH= new TH1F("fPtKinkMnH", "Neg P_{T}Kaon Kink  distrib.in High Pt Binning",47, gPt13High); 
//
   fListOfHistos=new TList();

   fListOfHistos->Add(fMultiplMC);
   fListOfHistos->Add(fIncompletEv);
   fListOfHistos->Add(fMultMCK);
   fListOfHistos->Add(fMultESDK);
   fListOfHistos->Add(fZpr);
   fListOfHistos->Add(fESDMult);
   fListOfHistos->Add(fZMainVx);
   fListOfHistos->Add(fMultVertex);
//
   fListOfHistos->Add(fHistPtESD);
   fListOfHistos->Add(fHistPt);
   fListOfHistos->Add(fHistQtAll);
   fListOfHistos->Add(fHistQt1);
   fListOfHistos->Add(fHistQt2);
   fListOfHistos->Add(fHistPtKaon);
   fListOfHistos->Add(fHistPtKPDG);
   fListOfHistos->Add(fHistEta);
   fListOfHistos->Add(fHistEtaK);
   fListOfHistos->Add(fptKMC);
   fListOfHistos->Add(fgenpt);
   fListOfHistos->Add(frad);
   fListOfHistos->Add(fKinkKaon);
   fListOfHistos->Add(fKinKRbn);
   fListOfHistos->Add(fKinkKaonBg);
   fListOfHistos->Add(fM1kaon);
   fListOfHistos->Add(fPtKink);
   fListOfHistos->Add(fptKink);
   fListOfHistos->Add(fAngMomK);
   fListOfHistos->Add(fAngMomPi);
   fListOfHistos->Add(fAngMomKC);
   fListOfHistos->Add(fSignPtNcl);
   fListOfHistos->Add(fSignPtEta);
   fListOfHistos->Add(fEtaNcl);
   fListOfHistos->Add(fSignPt);
   fListOfHistos->Add(fChi2NclTPC);
   fListOfHistos->Add(fRatChi2Ncl);
   fListOfHistos->Add(fRadiusNcl);
   fListOfHistos->Add(fTPCSgnlP);
   fListOfHistos->Add(fTPCSgnlPa);
   fListOfHistos->Add(fRpr);
   fListOfHistos->Add(fdcatoVxXY);
   fListOfHistos->Add(fnSigmToVx);
   fListOfHistos->Add(fKinkMothDau);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fHistPtKaoP);
   fListOfHistos->Add(fHistPtKaoN);
   fListOfHistos->Add(frapiKESD);
   fListOfHistos->Add(flifetime);
   fListOfHistos->Add(fradLK);
   fListOfHistos->Add(fradPtRpDt);
   fListOfHistos->Add(fInvMuNuAll);
   fListOfHistos->Add(fQtInvM);
   fListOfHistos->Add(fDCAkink);
   fListOfHistos->Add(fPosiKink);
   fListOfHistos->Add(fPosiKinkK);
   fListOfHistos->Add(fPosiKinKXZ);
   fListOfHistos->Add(fPosiKinKYZ);
   fListOfHistos->Add(fPosiKinKBg);
   fListOfHistos->Add(fQtMothP);
   fListOfHistos->Add(fTPCSgnlPtpc);
   fListOfHistos->Add(fTPCMomNSgnl);
   fListOfHistos->Add(fMothKinkMomSgnl);
   fListOfHistos->Add(fNSigmTPC);
   fListOfHistos->Add(fTPCSgnlKinkDau);
   fListOfHistos->Add(fPtKinkPos);
   fListOfHistos->Add(fPtKinkNeg);
   fListOfHistos->Add(fRadNclCln);
   fListOfHistos->Add(fRatioCrossedRows);
   fListOfHistos->Add(fRatioCrossedRowsKink);
   fListOfHistos->Add(fRadiusPt);
   fListOfHistos->Add(fRadiusPtcln);
   fListOfHistos->Add(fInvMassMuNuPt);
   fListOfHistos->Add(fInvMassMuNuPtAll);
   fListOfHistos->Add(fPtCut1);
   fListOfHistos->Add(fPtCut2);
   fListOfHistos->Add(fPtCut3);
  fListOfHistos->Add(fAngMomKKinks);
  fListOfHistos->Add(fPtKinkK0);
   fListOfHistos->Add(fPtKinkK0P);
   fListOfHistos->Add(fPtKinkK0N);
  fListOfHistos->Add(fPtKinkGyu);
   fListOfHistos->Add(fPtKinkGyuP);
   fListOfHistos->Add(fPtKinkGyuN);
   fListOfHistos->Add(fNumberOfEvent);
//   fListOfHistos->Add(fNumberOfEvent_cent);
   fListOfHistos->Add(fbgCleaningHigh);
//   fListOfHistos->Add(f1);
   fListOfHistos->Add(fMultiplicity);
   fListOfHistos->Add(fEventVertx);
   fListOfHistos->Add(fIncompletEvent);
   fListOfHistos->Add(fMultpileup);
fListOfHistos->Add(fcentrality);
   fListOfHistos->Add(fcentral0To1);
   fListOfHistos->Add(fcentral0To01);
//
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
   fListOfHistos->Add(fMultiBin11ChargedMulti);
   fListOfHistos->Add(fMultiBin12ChargedMulti);
   fListOfHistos->Add(fMultiBin13ChargedMulti);
   fListOfHistos->Add(fMultiBin14ChargedMulti);
fListOfHistos->Add(fMultiBin15ChargedMulti);
   fListOfHistos->Add(fMultiBin16ChargedMulti);
   fListOfHistos->Add(fMultiBin17ChargedMulti);
   fListOfHistos->Add(fMultiBin18ChargedMulti);
   fListOfHistos->Add(fMultiBin19ChargedMulti);
   fListOfHistos->Add(fMultiBin20ChargedMulti);
   fListOfHistos->Add(fMultiBin21ChargedMulti);
   fListOfHistos->Add(fMultiBin22ChargedMulti);
   fListOfHistos->Add(fMultiBin23ChargedMulti);
//
//
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
   fListOfHistos->Add(fMultiBin11KaonKinksPos);
   fListOfHistos->Add(fMultiBin12KaonKinksPos);
   fListOfHistos->Add(fMultiBin13KaonKinksPos);
   fListOfHistos->Add(fMultiBin14KaonKinksPos);
fListOfHistos->Add(fMultiBin15KaonKinksPos);
   fListOfHistos->Add(fMultiBin16KaonKinksPos);
   fListOfHistos->Add(fMultiBin17KaonKinksPos);
   fListOfHistos->Add(fMultiBin18KaonKinksPos);
   fListOfHistos->Add(fMultiBin19KaonKinksPos);
   fListOfHistos->Add(fMultiBin20KaonKinksPos);
   fListOfHistos->Add(fMultiBin21KaonKinksPos);
   fListOfHistos->Add(fMultiBin22KaonKinksPos);
   fListOfHistos->Add(fMultiBin23KaonKinksPos);
//
//

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
   fListOfHistos->Add(fMultiBin11KaonKinksNeg);
   fListOfHistos->Add(fMultiBin12KaonKinksNeg);
   fListOfHistos->Add(fMultiBin13KaonKinksNeg);
   fListOfHistos->Add(fMultiBin14KaonKinksNeg);

fListOfHistos->Add(fMultiBin15KaonKinksNeg);
   fListOfHistos->Add(fMultiBin16KaonKinksNeg);
   fListOfHistos->Add(fMultiBin17KaonKinksNeg);
   fListOfHistos->Add(fMultiBin18KaonKinksNeg);
   fListOfHistos->Add(fMultiBin19KaonKinksNeg);
   fListOfHistos->Add(fMultiBin20KaonKinksNeg);
   fListOfHistos->Add(fMultiBin21KaonKinksNeg);
   fListOfHistos->Add(fMultiBin22KaonKinksNeg);
   fListOfHistos->Add(fMultiBin23KaonKinksNeg);

   fListOfHistos->Add(fpercentMul);

  //DefineOutput(1, TList::Class());
  PostData(1, fListOfHistos);
}
//________________________________________________________________________
void AliAnalysisKinkESDMulPPpi::UserExec(Option_t *) 
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
// Number ESD tracks 
   Int_t nESDTracks =  esd->GetNumberOfTracks();
      fMultiplMC->Fill(nESDTracks);
      fNumberOfEvent->Fill(1.5);
      fMultiplicity->Fill(2);
// check incomplete events 
                if (esd->IsIncompleteDAQ()) return;
                     fIncompletEv ->Fill(esd->GetNumberOfTracks() );
      fNumberOfEvent->Fill(2.5);
                     fIncompletEvent ->Fill(2);
//
 // check of Pileup   3/2/2016
              if (esd->IsPileupFromSPD()) return;
       fMultESDK->Fill(nESDTracks);
       fMultpileup->Fill(2);
      fNumberOfEvent->Fill(3.5);
          //   HM trigger pp 13 TeV 
       Bool_t isSelected =
 ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kHighMultV0;

       if ( isSelected ==kFALSE) return;   //  24/6/11 apo MF
//*/

/*
                  UInt_t maskIsSelected =
     ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
     Bool_t isSelected = 0;
     //   isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
==?? isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kHighMultV0;    //  16/3/21
         if (!isSelected) {
      //   ( PostData(1, fListOfHistos));
            return;
     }
*/
      fMultMCK->Fill(nESDTracks);
//
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
//-----------------------------------------------------------------------------------------------
        cent = fMultSel->GetMultiplicityPercentile("V0M", kTRUE);
        //if ((cent < 0) || (cent > 90)) return; //Event selection
        //if ((cent < 0) || (cent >100)) return; //Event selection
        //   if ((cent < 0) || (cent >1)) return; //Event selection
        //if ((cent < 0) || (cent >0.01)) return; //Event selection     3/4/21
        if ((cent < 0.0) || (cent >0.1)) return; //Event selection    29/9/22
        }
        }

      fNumberOfEvent->Fill(4.5);
       fcentrality->Fill(cent);
          //

             // End Multiplicity bins
//
//*/


       // end Multiplicity selection
//      Printf("multiplicity percentile = %f", lPercentile);
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



       fMultESDK->Fill(refmultiplicity);
*/
/*
       // Multiplicity selection  apo Evi 18 Oct 2017

      AliMultSelection *fMultSel = (AliMultSelection*) esd -> FindListObject("MultSelection");

      // if(!fMultSel-> IsEventSelected()) return;
      if (!fMultSel) {
  //If you get this warning please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
        }
      Float_t lPercentile= fMultSel->GetMultiplicityPercentile("V0M", kTRUE);
    //   Feb 2018 test mss    if ((lPercentile < 0) || (lPercentile > 100)) return;
        if ((lPercentile < 0) || (lPercentile > 100)) return;
*/
//
//    
//       if (!selectVertex2015pp(lESDevent,kTRUE,kFALSE,kTRUE)){ //new common 
      if (!selectVertex2015pp(esd,kTRUE,kFALSE,kTRUE)){ //new common 
//       settings for MB analysis feb 17 2016 (mail David)
//1/marc         ( PostData(1, fListOfHistos));
 //             PostData(1, fList);
                            return;
                                 }
      fMultVertex->Fill(nESDTracks);
       fcentral0To1 ->Fill(cent);

//
  const AliESDVertex *vertex=GetEventVertex(esd);    // 22/8
  if(!vertex) return;
    fZMainVx->Fill(vertex->GetZ());
//-------------
       fpercentMul->Fill(cent);
//
 /* Multiplicitycentrality    bins,  move here 14/4/2018

             if ((cent>0)&&(cent<1)) {
               fMultiBin1ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>1.)&&(cent<5.0)) {
               fMultiBin2ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>5.)&&(cent<10.0 )) {
               fMultiBin3ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>10.0  )&&(cent<15.0 )) {
               fMultiBin4ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>15)&&(cent<20)) {
               fMultiBin5ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>20)&&(cent<30)) {
               fMultiBin6ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>30)&&(cent<40)) {
               fMultiBin7ChargedMulti->Fill(esd->GetNumberOfTracks());
 }
             if ((cent>40)&&(cent<50)) {
               fMultiBin8ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>50)&&(cent<60)) {
               fMultiBin9ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>60)&&(cent<70)) {
               fMultiBin10ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>70)&&(cent<80)) {
               fMultiBin11ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>80)&&(cent<90)) {
               fMultiBin12ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>90)&&(cent<100)) {
               fMultiBin13ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>.5)&&(cent<1.0)) {
               fMultiBin14ChargedMulti->Fill(esd->GetNumberOfTracks());
             }

*/           // End Multiplicity bins
       fcentral0To01->Fill(cent);
// Multiplicity bins
//
          if ((cent>0)&&(cent<0.001)) {
                        fMultiBin1ChargedMulti->Fill(esd->GetNumberOfTracks());
                }
        if ((cent>=0.001)&&(cent<0.005 )) {
             fMultiBin2ChargedMulti->Fill(esd->GetNumberOfTracks());
                   }
                if ((cent>=0.0   )&&(cent<0.005)) {
   fMultiBin3ChargedMulti->Fill(esd->GetNumberOfTracks());
                                          }
             if ((cent>=0.005)&&(cent<0.01 )) {
         fMultiBin4ChargedMulti->Fill(esd->GetNumberOfTracks());
           }
                if ((cent>=0.0 )&&(cent<0.01 )) {
             fMultiBin5ChargedMulti->Fill(esd->GetNumberOfTracks());
                         }
              if ((cent>=0.01 )&&(cent<0.05 )) {
         fMultiBin6ChargedMulti->Fill(esd->GetNumberOfTracks());
                                  }
             if ((cent>=0.050)&&(cent<0.10)) {
          fMultiBin7ChargedMulti->Fill(esd->GetNumberOfTracks());
                                  }
           if ((cent>=0.1)&&(cent<0.2)) {
       fMultiBin8ChargedMulti->Fill(esd->GetNumberOfTracks());
                      }
              if ((cent>=0.2)&&(cent<0.5)) {
         fMultiBin9ChargedMulti->Fill(esd->GetNumberOfTracks());
                       }
          if ((cent>= 0.5)&&(cent< 1.0)) {
      fMultiBin10ChargedMulti->Fill(esd->GetNumberOfTracks());
                    }
       if ((cent>= 0.0)&&(cent< 1.0)) {
               fMultiBin11ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 1.0)&&(cent< 5.0)) {
               fMultiBin12ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 5.0)&&(cent< 10.)) {
               fMultiBin13ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 10.)&&(cent< 15.)) {
               fMultiBin14ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 15.)&&(cent< 20.)) {
               fMultiBin15ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 20.)&&(cent< 30. )) {
               fMultiBin16ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 30.)&&(cent< 40.)) {
               fMultiBin17ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 40.)&&(cent< 50.)) {
               fMultiBin18ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 50.)&&(cent< 60.)) {
               fMultiBin19ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 60.)&&(cent< 70.)) {
               fMultiBin20ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 70.)&&(cent< 80.)) {
               fMultiBin21ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 80.)&&(cent< 90.)) {
               fMultiBin22ChargedMulti->Fill(esd->GetNumberOfTracks());
             }
             if ((cent>= 90.)&&(cent<100.)) {
               fMultiBin23ChargedMulti->Fill(esd->GetNumberOfTracks());
             }                                                             
           // End Multiplicity bins
//
  Double_t vpos[3];
  vertex->GetXYZ(vpos);
    fZpr->Fill(vpos[2]);         
 //     if (TMath::Abs( vpos[2] ) > 10. ) return;   

    
      fNumberOfEvent->Fill(5.5);
      fEventVertx->Fill(2);

  Double_t vtrack[3], ptrack[3];
  
     
 Int_t nESDTracK = 0;
// Int_t nESDTrKink = 0;

   Int_t nGoodTracks =  esd->GetNumberOfTracks();
    fESDMult->Fill(nGoodTracks);
      
       Double_t nsigmall = 100.0;
       Double_t nsigma = 100.0;
       Double_t nsigmaPi = 100.0;
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
//*/  
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
//	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkDau)-1);
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


// track loop
//
   for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {

    AliESDtrack* track = esd->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    fHistPt->Fill(track->Pt());


     //    sigmas
//  mss2015     nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
     nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
  //  nsigmaPion= (fESDpid->NumberOfSigmasTPC(track,AliPID::kPion));
//  mss 2015      nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
  // kaon     nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
   nsigma     = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track  , AliPID::kPion));// 26/10 eftihis
   nsigmaPi   = (fPIDResponse->NumberOfSigmasTPC(track  , AliPID::kPion));//    3/2021  PPpi pions

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

   	    	Int_t indexKinkPos=track->GetKinkIndex(0);   // kink index 

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
 fXvYv->Fill(vtrack[0],vtrack[1]);  
 fZvYv->Fill(vtrack[0],vtrack[2]);  
 fZvXv->Fill(vtrack[1],vtrack[2]);  

// track momentum, rapidity calculation
     track->GetPxPyPz(ptrack);
    
    TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
    
// K-rapidity calcualtion 
          Double_t   etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.13957  *0.13957   );
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
 
//  14/2/13 /================/   if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5))  nESDTrKink++;  //  count of second  23Jul11    
         if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5))
 //     if((TMath::Abs(dcaToVertexXYpos)>0.14)||(TMath::Abs(dcaToVertexZpos)>1.5))
          continue;   //    

//                    if (!fMaxDCAtoVtxCut->AcceptTrack(track)) continue;

    fdcatoVxXY->Fill(dcaToVertexXYpos);
//
    
//  track Mult. after selection 
    nESDTracK++;        
  //    
//=========================================
    fHistPtESD->Fill(track->Pt());

   // Add Kink analysis           =============================
   
// daughter kink 
//if(indexKinkPos >0)fTPCSgnlKinkDau->Fill(track->P(), (track->GetTPCsignal()  ) ) ;  //  daughter kink 
    
//  loop on kinks
		if(indexKinkPos<0){     ////mother kink

             fptKMC   ->Fill(  track->Pt()    );  // Pt from tracks , all kinks

    fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);
//
	// select kink class	
 //
                  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
                        //

         // DCA kink
          Double_t  Dist2 = kink->GetDistance();
         //   fDCAkink->Fill( Dist2   );
            //   if (Dist2 > 0.2) continue; //  systematics 11/8/11 
//
	
// TPC mother momentum 
      
	   const TVector3 vposKink(kink->GetPosition());
 fPosiKink ->Fill( vposKink[0], vposKink[1]  );
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
// Kink  mother momentum 
//     Double_t trMomTPCKink=motherMfromKink.Mag();        
// TPC mother momentun
     Double_t trMomTPC=track->GetTPCmomentum();      
  //     fTPCSgnlKinkDau->Fill( daughterMKink.Mag(), dEdxKinkDau  ) ;  //  daughter kink 
  //   
           fHistQtAll->Fill(qT) ;  //  Qt   distr
                  
           fptKink->Fill(motherMfromKink.Pt()); /// pt from kink

           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

//   rapiditya nd pt selection 
    //       if(  (TMath::Abs(rapiditK )) > 0.7 ) continue;
         //    if(  (TMath::Abs(rapiditK )) > 0.5 ) continue; //  allagh  Nov. 2014 , better acceptance 
//============= if(  (TMath::Abs(rapiditK )) > fRapiK ) continue; //  allagh  Nov. 2014 , better acceptance 
        if ( (track->Pt())<.150)continue;  // new March2020
//              eta selection 
      // mss2015  if ( TMath::Abs(trackEta) > 0.8) continue;  // new  NOv   2014 
        if ( TMath::Abs(trackEta) > 0.8) continue;  // new  NOv   2014 
            fHistEta->Fill(trackEta) ;  //   Eta distr 
//
 if(  (TMath::Abs(rapiditK )) > fRapiK ) continue; //   18/3/2020

                fQtMothP->Fill( track->P(), qT);

        //if ( qT< fLowQt )  fHistQt1  ->Fill(qT) ;  //  Qt   distr

            fKinkKaonBg->Fill(motherPt);     

//          maximum decay angle at a given mother momentum
	   //Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	   Double_t maxDecAngKmu=f1->Eval(track->P()          ,0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(    track->P(),       0.,0.,0.);

//  fake kinks are removed 
// 28/3/20               if( (kinkAngle<1.)  ) continue;

           
      //  BG  ?????==============
             if ( TMath::Abs(vposKink[2]) >  225. ) continue ;
       //      if ( TMath::Abs(vposKink[2]) <  0.5 ) continue ;
            fPtCut1   ->Fill(trackPt );     

                   //    oxi afto 4/4/20  if (vposKink[2]  > 0.3 && vposKink[2]   < 0.7) continue;  //  mss 15/2/2018   test

          fHistQt2->Fill(qT);  //             
//    remove background====================
  if (qT < 0.010) continue; 
            fAngMomPi->Fill( track->P(),           kinkAngle); 

        if ( qT< fLowQt )  fHistQt1  ->Fill(qT) ;  //  Qt   distr
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
           fQtInvM -> Fill ( invariantMassKmu,  qT);
           fInvMuNuAll->Fill(invariantMassKmu);
           fInvMassMuNuPtAll ->Fill(invariantMassKmu,  trackPt);
//
         fRadiusPt->Fill( kink->GetR(), trackPt); // 
 //  radius and Minv selection 
    //   if( ( kink->GetR()> 120 ) && ( kink->GetR() < 210 )  )  {
       if( ( kink->GetR()> fKinkRadLow ) && ( kink->GetR() <fKinkRadUp   )  )  {
    //  for systematics   if( ( kink->GetR()> 130 ) && ( kink->GetR() < 200 )  )  {
if (qT<fLowQt )  fAngMomKC->Fill(track->P(), kinkAngle); 
      //
          if ( qT< fLowQt ) fM1kaon->Fill(invariantMassKmu);
          //   if ( qT>  0.12  ) fM1kaon->Fill(invariantMassKmu);
              if ( qT < fLowQt) 
         fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
  }    
//  tails cleaning
    if(  ( tpcNCl<fLowCluster) ) continue;  // test 27 feb 2012 ,, OK
//              if(  ( tpcNCl<20 ) ) continue;  // test 15 March  13,, OK
// cleaning BG in tails
      //Int_t tpcNClHigh = -31.67+ (11./12.)  *( kink->GetR() ) ;  
      Int_t tpcNClHigh = -42.67+ (11./12.)  *( kink->GetR() ) ;  
               if ( tpcNCl > tpcNClHigh) continue;   
                  
      //Int_t tpcNClMin  = -85.5 + (65./95.)  *( kink->GetR() ) ;  
      Int_t tpcNClMin  = -60.5 + (65./95.)  *( kink->GetR() ) ;  
               if ( tpcNCl < tpcNClMin ) continue;   

   //  back, 20/1/2013   if (ratioCrossedRowsOverFindableClustersTPC< 0.5) continue;// test 14/1/2013 
//
               fHistPtKPDG->Fill(track->Pt());  // ALL  K-candidates until now                 
//  maximum angles selection with some error cut 
        if( (kinkAngle>maxDecAngpimu*1.2) ) continue;   // fpr pions
//-------------------------------------------------------------------------
    //  if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=210.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.6)){
 //  march2016    if((kinkAngle>maxDecAngpimu)&&(qT>0.120)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=210.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.8)){
     //if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.8)){
//      if((kinkAngle>maxDecAngpimu)&&(qT>fLowQt)&&(qT<0.30)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(TMath::Abs(rapiditK)<fRapiK  )&&(invariantMassKmu<0.8)){
if((qT<fLowQt)&&((kink->GetR()>= fKinkRadLow )&&(kink->GetR()<= fKinkRadUp ))&&(invariantMassKmu<0.2)){
  // systematics   if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>=130.)&&(kink->GetR()<=200.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.8)){
//
        fAngMomKKinks->Fill(track->P(), kinkAngle); 
            fPtCut2   ->Fill(trackPt );     

            fPtCut3   ->Fill(trackPt );     
//  here the kaons selected by the decay features
           fTPCSgnlPa->Fill( track->GetInnerParam()->GetP() ,(track->GetTPCsignal()  ) ) ;
//
         //               if ( nsigma               > 3.0) continue;
                   //   mexri 23/3/2020   if ( nsigma               > 3.0) continue; 
                   // if ( nsigma               > 2.0) continue; 
                   //         ????    if ( nsigmaPion< -4.5 || nsigmaPion> 1.5     ) continue;// 27/32020 
                   //if ( nsigmaPi< -3.5 || nsigmaPi> 3.5     ) continue;// 27/32020 
                   if ( nsigmaPi< -3.0 || nsigmaPi> 3.0     ) continue;// 27/32020 
// 
//  next plots for the identified kaons by the kink analysis

                                     fHistPtKaon->Fill(track->Pt());   //
//
              if(tpcSign >0.)        fHistPtKaoP->Fill( track->Pt()         ) ;   //
               if ( tpcSign <0.)    fHistPtKaoN->Fill(  track->Pt()        ) ;   //
          fTPCSgnlP->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()  ) ) ;
         fRadNclCln->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
         fRadiusPtcln->Fill( kink->GetR(), trackPt); // 
           fInvMassMuNuPt ->Fill(invariantMassKmu,  trackPt);
                 
         fTPCSgnlPtpc->Fill(trMomTPC  , (track->GetTPCsignal()  ) ) ;               //trMomTPC
         //fMothKinkMomSgnl ->Fill(trMomTPCKink  , (track->GetTPCsignal()  ) ) ;
         fMothKinkMomSgnl ->Fill(  dEdxKinkDau  , (track->GetTPCsignal()  ) ) ;
//
       fTPCSgnlKinkDau->Fill( daughterMKink.Mag(), dEdxKinkDau  ) ;  //  daughter kink 
// 
         fTPCMomNSgnl->Fill(trMomTPC ,nsigmaPi  );     
      fNSigmTPC   ->Fill(nsigmaPi   ); //new 22/3/20 prosoxi nsigmaPion , nsigma Abs(sigmaPion)    
//
                frad->Fill(kink->GetR());  // kink 
               fradLK->Fill(lifeKink    );  // kink 
             fHistEtaK->Fill(trackEta);
            frapiKESD ->Fill(rapiditK);  //  rapidityof kaons 
          fPosiKinKBg->Fill( vposKink[2], kink->GetR() );

                     Float_t signPt= tpcSign*trackPt;
                  fSignPtNcl->Fill( signPt  ,   tpcNCl   );   ///  28/4/2010
                  fSignPtEta->Fill( signPt  , rapiditK  );
                  fEtaNcl->Fill( rapiditK, tpcNCl    );
                  fSignPt->Fill( signPt );
                  fChi2NclTPC->Fill( (track->GetTPCchi2() ) ,  tpcNCl );
         fRatChi2Ncl-> Fill (  (track->GetTPCchi2()/track->GetTPCclusters(0)  )) ;
 //   fdcatoVxXY->Fill(dcaToVertexXYpos);
           //if( dEdxKinkDau> 1.5* (track->GetTPCsignal()   )  )      fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
       //    if((dEdxKinkDau>  80. ) && (dEdxKinkDau > 4.*nsigmaPion)   )      fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
//  mss 2015          if (nsigmaPion>  3.             ) fTPCSgnlPtpc->Fill( daughterMKink.Mag(),  dEdxKinkDau    ) ;
         //if (TMath::Abs(dEdxKinkDau -(track->GetTPCsignal() )> 10. )) fTPCSgnlPtpc->Fill( daughterMKink.Mag(),  dEdxKinkDau    ) ;
//		        flifetime->Fill(( lifeKink*.493667   )  /track->P()   ) ;
             fKinkKaon ->Fill(track->Pt());        
            fDCAkink->Fill(dcaToVertexXYpos );

               fPtKink->Fill(track->Pt()); ///  Comb  bins     
              if(tpcSign >0.)        fPtKinkPos ->Fill( track->Pt()         ) ;   //K-plus bins Comb 
              if(tpcSign <0.)        fPtKinkNeg ->Fill( track->Pt()         ) ;   //K-minus bins Comb   

               fPtKinkK0->Fill(track->Pt()); ///  K0    bins     
              if(tpcSign >0.)        fPtKinkK0P->Fill( track->Pt()         ) ;   //K-plus bins K0 Peter 
              if(tpcSign <0.)        fPtKinkK0N ->Fill( track->Pt()         ) ;   //K-minus bins K0 Peter 

               fPtKinkGyu->Fill(track->Pt()); ///   K-charged High Pt Gyula
              if(tpcSign >0.)        fPtKinkGyuP->Fill( track->Pt()    ) ;   //K-plus charged high pt ,Gyula
              if(tpcSign <0.)        fPtKinkGyuN ->Fill( track->Pt()  ) ;   //K-minus bins High Pt charged Gyula  

 // Multiplicity bins

             //   wrong 14/4   if ((cent>0)&&(cent<1)) {
             if ((cent>0.)&&(cent<0.001)) {
               if (tpcSign >0.) fMultiBin1KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin1KaonKinksNeg->Fill(track->Pt());
             }
             //  wrong  14/4 if ((cent>0)&&(cent<.001 )) {
             if ((cent>0.001)&&(cent<.005 )) {
               if (tpcSign >0.) fMultiBin2KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin2KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.  )&&(cent<0.005 )) {
               if (tpcSign >0.) fMultiBin3KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin3KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.005 )&&(cent<0.01 )) {
               if (tpcSign >0.) fMultiBin4KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin4KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0. )&&(cent<0.01 )) {
               if (tpcSign >0.) fMultiBin5KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin5KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.01 )&&(cent<0.05 )) {
               if (tpcSign >0.) fMultiBin6KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin6KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.05 )&&(cent<0.1 )) {
               if (tpcSign >0.) fMultiBin7KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin7KaonKinksNeg->Fill(track->Pt());
 }
             if ((cent> 0.1  )&&(cent<0.2 )) {
               if (tpcSign >0.) fMultiBin8KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin8KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent> 0.2  )&&(cent<0.5  )) {
               if (tpcSign >0.) fMultiBin9KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin9KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.5 )&&(cent< 1.0 )) {
               if (tpcSign >0.) fMultiBin10KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin10KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>0.0 )&&(cent< 1.0 )) {
               if (tpcSign >0.) fMultiBin11KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin11KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>1.0 )&&(cent< 5.0 )) {
               if (tpcSign >0.) fMultiBin12KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin12KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent> 5.)&&(cent<10.0)) {
               if (tpcSign >0.) fMultiBin13KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin13KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>10.)&&(cent<15.0)) {
               if (tpcSign >0.) fMultiBin14KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin14KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>15)&&(cent<20)) {
               if (tpcSign >0.) fMultiBin15KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin15KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>20)&&(cent<30)) {
               if (tpcSign >0.) fMultiBin16KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin16KaonKinksNeg->Fill(track->Pt());
             }
//
             if ((cent>30)&&(cent<40)) {
               if (tpcSign >0.) fMultiBin17KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin17KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>40)&&(cent<50)) {
               if (tpcSign >0.) fMultiBin18KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin18KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>50)&&(cent<60)) {
               if (tpcSign >0.) fMultiBin19KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin19KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>60)&&(cent<70)) {
               if (tpcSign >0.) fMultiBin20KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin20KaonKinksNeg->Fill(track->Pt());
             }
//*
             if ((cent>70)&&(cent<80)) {
               if (tpcSign >0.) fMultiBin21KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin21KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>80)&&(cent<90)) {
               if (tpcSign >0.) fMultiBin22KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin22KaonKinksNeg->Fill(track->Pt());
             }
             if ((cent>90)&&(cent<100)) {
               if (tpcSign >0.) fMultiBin23KaonKinksPos->Fill(track->Pt());
               if (tpcSign <0.) fMultiBin23KaonKinksNeg->Fill(track->Pt());
             }

//*/             // End Multiplicity bins
//

             // End Multiplicity bins
//
             fKinKRbn->Fill(track->Pt());       // TOF      

              fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK); 
               fAngMomK->Fill(    track->P(),        kinkAngle); 
       fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[2], vposKink[0]  );
          fPosiKinKYZ->Fill( vposKink[2], vposKink[1]  );

        }  //  kink selection 
                  

	}  //End Kink Information    
  

  } //track loop 
//  } //daughter loop 

 //     fMultiplMC->Fill(nESDTracK );

  PostData(1, fListOfHistos);

}      

//________________________________________________________________________
void AliAnalysisKinkESDMulPPpi::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}

//____________________________________________________________________//

  // Get the vertex 
const AliESDVertex* AliAnalysisKinkESDMulPPpi::GetEventVertex(const AliESDEvent* esd) const
  
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
Bool_t AliAnalysisKinkESDMulPPpi::selectVertex2015pp(AliESDEvent *esd,
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
            Bool_t AliAnalysisKinkESDMulPPpi::IsGoodSPDvertexRes(const AliESDVertex  * spdVertex)
              {
            if (!spdVertex) return kFALSE;
                if (spdVertex->IsFromVertexerZ() && !(spdVertex->GetDispersion()<0.04
          && spdVertex->GetZRes()<0.25)) return kFALSE;
                              return kTRUE;
                                  }

