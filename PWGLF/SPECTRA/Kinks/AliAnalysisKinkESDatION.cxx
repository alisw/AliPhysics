/**************************************************************************
 * Authors: Martha Spyropoulou-Stassinaki and the  members
 * of the Greek group at Physics Department of Athens University
 * Paraskevi Ganoti and Anastasia Belogianni 
 **************************************************************************/
//-----------------------------------------------------------------
//                 AliAnalysisKinkESDatION class
//     Example of an analysis task for kink topology study in Pb-Pb collisions
//      Kaons from kink topology are 'identified' in this code
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
#include "AliAnalysisKinkESDatION.h"
#include "AliStack.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
#include "AliPhysicsSelectionTask.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliPIDResponse.h"
//  #include "AliTPCpidESD.h"

ClassImp(AliAnalysisKinkESDatION)


//________________________________________________________________________
AliAnalysisKinkESDatION::AliAnalysisKinkESDatION(const char *name) 
  : AliAnalysisTaskSE(name), fHistPtESD(0),fHistPt(0),fHistQtAll(0),fHistQt1(0),fHistQt2(0)
  , fHistPtKaon(0),fHistPtKPDG(0),fHistEta(0),fHistEtaK(0),fptKMC(0),fMultiplMC(0),fESDMult(0),frad(0),
  fKinkKaon(0),fKinKRbn(0), fKinkKaonBg(0), fM1kaon(0),  fPtKink(0),  fptKink(0),
    fAngMomK(0),fAngMomPi(0), fAngMomKC(0),  fMultESDK(0), fMultMCK(0),
 fSignPtNcl(0), fSignPtEta(0), fEtaNcl(0), fSignPt(0), fChi2NclTPC(0), fRatChi2Ncl(0), fRadiusNcl(0), fTPCSgnlP(0),
   fTPCSgnlPa(0), fRpr(0),fZpr(0), fdcatoVxXY(0),  fKinkMothDau(0),
 fZvXv(0),fZvYv(0), fXvYv(0), fPtPrKink(0), fHistPtKaoP(0), fHistPtKaoN(0),frapiKESD(0), flifetime(), fradLK(0),
    fradPtRpDt(0), fInvMuNuAll(0), fQtInvM(0), 
         fDCAkink(0), fPosiKink(0),  fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),  fQtMothP(0),
           fKinkKPt05(0), fKinkKPt510(0), fKinkKPt1020(0), fKinkKPt2030(0), fKinkKPt3040(0),
           fKinkKPt4050(0), fKinkKPt5060(0), fKinkKPt6070(0), fKinkKPt7080(0), fKinkKPt8090(0),
           fKinkMul05(0), fKinkMul510(0), fKinkMul1020(0), fKinkMul2030(0), fKinkMul3040(0),
           fKinkMul4050(0), fKinkMul5060(0), fKinkMul6070(0), fKinkMul7080(0), fKinkMul8090(0),
            fRadiusNclAll(0), fRadiusNclK(0), fRadiusNclClean(0), fRatioCrossedRows(0),fRatioCrossedRowsKink(0),
 f1(0), f2(0),
      fListOfHistos(0)   ,fMaxDCAtoVtxCut(0), fPIDResponse(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
 // DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkESDatION::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
                   fMaxDCAtoVtxCut=new AliESDtrackCuts("fMaxDCA", "fMaxDCA");
       fMaxDCAtoVtxCut->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fMaxDCAtoVtxCut->SetMaxChi2TPCConstrainedGlobal(36);

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
 //  Double_t gPt[31] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
   //                     1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
     //                    2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6,3.9, 4.2, 4.5, 4.8};
   Double_t gPt7[43] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0, 
                         3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.8,5.2, 5.6, 6.0,  7.0, 8.0,10.0 };
    Double_t gPtK0[45] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.2, 2.4, 2.6, 2.8,  3.0,   3.3, 3.6, 3.9,
                         4.2, 4.4,5.0, 5.4, 5.9,  6.5,   7.0,7.5, 8.0,8.5,  9.2, 10., 11., 12., 13.5,15.0 };  // David K0

   Double_t gPt7TOF[47] = { 0.2,0.25, 0.3,0.35,  0.4,0.45,  0.5,0.55,  0.6,0.65,  0.7,0.75,  0.8, 0.85, 0.9, 0.95, 1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0,
                         3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0 };//  TOF Barbara

  fHistPtESD = new TH1F("fHistPtESD", "P_{T} distribution",100, 0.0,10.0);
  //fHistPtESD->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  //fHistPtESD->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  //fHistPtESD->SetMarkerStyle(kFullCircle);
  fHistPt = new TH1F("fHistPt", "P_{T} distribution",100, 0.0,10.0); 
  fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300); 
  fHistQt1= new TH1F("fHistQt1", "Q_{T} distribution",100, 0.0,.300); 
  fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300); 
  fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",150, 0.0,15.0); 
  fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",150, 0.0,15.0); 
  fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3); 
  fHistEtaK= new TH1F("fHistEtaK", "EtaK distribution", 26,-1.3, 1.3); 
  fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",100, 0.0,10.0); 
  fMultiplMC= new TH1F("fMultiplMC", "charge multiplicity MC",100, 0.5,4000.5); 
  fESDMult= new TH1F("fESDMult", "charge multipliESD", 100, 0.5,4000.5); 
   frad= new TH1F("frad", "radius  K generated",100, 0.,1000.0);
  fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi",150, 0.0,15.0); 
  fKinKRbn= new TH1F("fKinKRbn", "p_{t}Kaon kinks identi[GeV/c],Entries",44,gPtK0); 
  fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",150, 0.0,15.0); 
  fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",80,0.0, 0.8); 
  fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  bution",150, 0.0,15.0); 
  fptKink= new TH1F("fptKink", "P_{T}Kaon Kink  bution",   42, gPt7  ); 
  fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,10.0,80,0.,80.);
  fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,10.0,80,0.,80.);
  fAngMomKC= new TH2F("fAngMomKC","Decay angle vrs Mother Mom,K",100,0.0,10.0,80,0.,80.);
  fMultESDK=new TH1F("fMultESDK", "charge multipliESD kaons",100, 0.5,4000.5); 
  fMultMCK=new TH1F("fMultMCK", "charge multipli MC kaons",100, 0.5,4000.5); 
  fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",100,-5.,5.0,90,0.,180.);
  fSignPtEta= new TH2F("fSignPtEta","SignPt vrs Eta,K",100,-5.0,5.0,30,-.75,.75);
  fEtaNcl= new TH2F("fEtaNcl","Eta vrs nclust,K",30,-.75,.75, 90,0, 180);
  fSignPt= new TH1F("fSignPt","SignPt ,K",100,-5.0,5.0);
  fChi2NclTPC= new TH2F("fChi2NclTPC","Chi2vrs nclust,K",100,0.,500., 90,0, 180);
  fRatChi2Ncl= new TH1F("fRatChi2Ncl","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
  fRadiusNcl = new TH2F("fRadiusNcl","kink radius vrs Nclust, Dau",75,100.,250., 90,0, 180);
    fTPCSgnlP = new TH2F("fTPCSgnlP","TPC signal de/dx Mom,K",200,0.0,20.0,100,0.,300.);
  fTPCSgnlPa= new TH2F("fTPCSgnlPa","TPC signal de/dx Mom,K",200,0.0,20.,100, 0.,300.);
  fRpr = new TH1D("fRpr", "rad distribution  PID pr",50,0.0, 2.5);
  fZpr = new TH1D("fZpr", "z distribution PID pr  ",80,-20.,20.);
  fdcatoVxXY = new TH1D("fdcatoVxXY", "dca  distribution PID  ",80,-4.,4.);
  fKinkMothDau= new TH2F("fKinkMothDau","TPC kink Moth Daugh ,K",50,0.0,2.5,50, 0., 2.5);
  fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5,60, -15., 15.0);
  fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
  fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
  fPtPrKink=new TH1F("fPtPrKink","pt of ESD  kaonKink tracks",  46, gPt7TOF);
  fHistPtKaoP = new TH1F("fHistPtKaoP", "P_{T}KaonP  distribution",150, 0.0,15.0); 
  fHistPtKaoN = new TH1F("fHistPtKaoN", "P_{T}KaonN  distribution",150, 0.0,15.0); 

  frapiKESD=new TH1F("frapiKESD","rapid Kdistribution", 26,-1.3, 1.3); 
  flifetime= new TH1F("flifetime", "ct study of K-kinks",100,0.,1000.); 
  fradLK= new TH1F("fradLK", "perscentile  in Pb-Pb10",100,0,100); 
  fradPtRpDt=new TH3F("fradPtRpDt","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
  fInvMuNuAll= new TH1F("fInvMuNuAll", " Inv Mass MuNu all kink",80,0.,0.8); 
  fQtInvM= new TH2F("fQtInvM", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300); 
    fDCAkink = new TH1F("fDCAkink ", "DCA kink vetrex ",50, 0.0,1.0);

  fPosiKink= new TH2F("fPosiKink", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fQtMothP = new TH2F("fQtMothP", " Qt vrs Mother P", 100, 0., 5.0,100, 0.,0.300);
  fKinkKPt05 = new TH1F("fKinkKPt05", "P_{T}Kaon kinks identi centra0",150, 0, 15. ); 
  fKinkKPt510  = new TH1F("fKinkKPt510", "P_{T}Kaon kinks identi central5",150, 0, 15. ); 
  fKinkKPt1020 = new TH1F("fKinkKPt1020", "P_{T}Kaon kinks identi centra10",150, 0., 15.0 ); 
  fKinkKPt2030 = new TH1F("fKinkKPt2030", "P_{T}Kaon kinks identi centra20", 150, 0., 15.0 ); 
  fKinkKPt3040 = new TH1F("fKinkKPt3040", "P_{T}Kaon kinks identi centra30",150, 0., 15.0 ); 
  fKinkKPt4050 = new TH1F("fKinkKPt4050", "P_{T}Kaon kinks identi centra40",150, 0., 15.0 ); 
  fKinkKPt5060 = new TH1F("fKinkKPt5060", "P_{T}Kaon kinks identi centra50",150, 0., 15.0 ); 
  fKinkKPt6070 = new TH1F("fKinkKPt6070", "P_{T}Kaon kinks identi centra60",150, 0., 15.0 ); 
  fKinkKPt7080 = new TH1F("fKinkKPt7080", "P_{T}Kaon kinks identi centra70",150, 0., 15.0 ); 
  fKinkKPt8090 = new TH1F("fKinkKPt8090", "P_{T}Kaon kinks identi centra80",150, 0., 15.0 ); 
  fKinkMul05 = new TH1F("fKinkMul05", "charge Multipl  kink-Kaons centr05",100, 2000.5, 12000.5  ); 
  fKinkMul510  = new TH1F("fKinkMul510","charge Multipl  kink-Kaons centr0510",100, 2000.5,12000.5 ); 
  fKinkMul1020 = new TH1F("fKinkMul1020", "charge Multipl  kink-Kaons centr1020",100, 0.5, 10000.5  ); 
  fKinkMul2030 = new TH1F("fKinkMul2030", " charge Multipl  kink-Kaons centr1020",100, 2000.5, 7000.5); 
  fKinkMul3040 = new TH1F("fKinkMul3040", "charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fKinkMul4050 = new TH1F("fKinkMul4050", "charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fKinkMul5060 = new TH1F("fKinkMul5060", " charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fKinkMul6070 = new TH1F("fKinkMul6070", "charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fKinkMul7080 = new TH1F("fKinkMul7080", "charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fKinkMul8090 = new TH1F("fKinkMul8090", "charge Multipl  kink-Kaons centr1020",100, 0.5, 5000.5  ); 
  fRadiusNclAll = new TH2F("fRadiusNclAll","kink radius vrs Nclust,K",75,100.,250., 90,0, 180);
  fRadiusNclK   = new TH2F("fRadiusNclK","kink radius vrs Nclust,K",75,100.,250., 90,0, 180);
  fRadiusNclClean = new TH2F("fRadiusNclClean","kink radius vrs Nclust,K",75,100.,250., 90,0, 180);
fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);
  fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);


   fListOfHistos=new TList();
   fListOfHistos->SetOwner();

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
   fListOfHistos->Add(fMultiplMC);
   fListOfHistos->Add(fESDMult);
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
   fListOfHistos->Add(fMultESDK);
   fListOfHistos->Add(fMultMCK);
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
   fListOfHistos->Add(fZpr);
   fListOfHistos->Add(fdcatoVxXY);
   fListOfHistos->Add(fKinkMothDau);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fPtPrKink);
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
   fListOfHistos->Add(fQtMothP);
   fListOfHistos->Add(fKinkKPt05);
   fListOfHistos->Add(fKinkKPt510);
   fListOfHistos->Add(fKinkKPt1020);
   fListOfHistos->Add(fKinkKPt2030);
   fListOfHistos->Add(fKinkKPt3040);
   fListOfHistos->Add(fKinkKPt4050);
   fListOfHistos->Add(fKinkKPt5060);
   fListOfHistos->Add(fKinkKPt6070);
   fListOfHistos->Add(fKinkKPt7080);
   fListOfHistos->Add(fKinkKPt8090);
   fListOfHistos->Add(fKinkMul05);
   fListOfHistos->Add(fKinkMul510);
   fListOfHistos->Add(fKinkMul1020);
   fListOfHistos->Add(fKinkMul2030);
   fListOfHistos->Add(fKinkMul3040);
   fListOfHistos->Add(fKinkMul4050);
   fListOfHistos->Add(fKinkMul5060);
   fListOfHistos->Add(fKinkMul6070);
   fListOfHistos->Add(fKinkMul7080);
   fListOfHistos->Add(fKinkMul8090);
   fListOfHistos->Add(fRadiusNclAll);
   fListOfHistos->Add(fRadiusNclK);
   fListOfHistos->Add(fRadiusNclClean);
   fListOfHistos->Add(fRatioCrossedRows);
   fListOfHistos->Add(fRatioCrossedRowsKink);


  PostData(1, fListOfHistos);
}

//________________________________________________________________________
void AliAnalysisKinkESDatION::UserExec(Option_t *) 
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

//
    fMultMCK->Fill(esd->GetNumberOfTracks());
//  
///*
//   ==================check of Physics selection?
       Bool_t isSelected =
//    ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB;
((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kCentral;// 10/3/14

        if ( isSelected ==kFALSE) return;   //  

   Int_t nESDTracks =  esd->GetNumberOfTracks();
    fMultiplMC->Fill(nESDTracks);
//*/
//

     AliCentrality *esdCentrality = esd->GetCentrality();
//

//

    Float_t percent=  esdCentrality->GetCentralityPercentile("V0M");
//           Printf("percentntile:%i", percent);
            fradLK->Fill(percent     );  //  
         if ( percent    >= .0  && percent < 5.)     fKinkMul05 -> Fill ( nESDTracks  ); // percent 0-5 
         if ( percent    >= 5.  && percent < 10.) fKinkMul510->Fill( nESDTracks  );  // percent  5-10 
        if ( percent    >= 10. && percent < 20.)  fKinkMul1020->Fill( nESDTracks  );  // percent 10-20
         if ( percent    >= 20. && percent < 30.) fKinkMul2030->Fill( nESDTracks  );  // percent 20-30 
       if ( percent    >= 30. && percent < 40.) fKinkMul3040->Fill( nESDTracks  );  // percent 30-40 
       if ( percent    >= 40. && percent < 50.) fKinkMul4050->Fill( nESDTracks  );  // percent 40-50  
       if ( percent    >= 50. && percent < 60.) fKinkMul5060->Fill( nESDTracks  );  // percent 50-60 
       if ( percent    >= 60. && percent < 70.) fKinkMul6070->Fill( nESDTracks  );  // percent 60-70 
        if ( percent    >= 70. && percent < 80.) fKinkMul7080->Fill( nESDTracks  );  // percent 70-80 
        if ( percent    >= 80. && percent < 90.) fKinkMul8090->Fill( nESDTracks  );  // percent 80-90 

//    if ( centrality >=0 )    fMultESDK->Fill( percent);

//

  const AliESDVertex *vertex=GetEventVertex(esd);    // 22/8
  if(!vertex) return;
//
  Double_t vpos[3];
  vertex->GetXYZ(vpos);
    fZpr->Fill(vpos[2]);         
     if (TMath::Abs( vpos[2] ) > 10. ) return;  // main vertex selection   

    

  Double_t vtrack[3], ptrack[3];
  
     
 // Int_t nESDTracK = 0;

   Int_t nGoodTracks =  esd->GetNumberOfTracks();
    fESDMult->Fill(nGoodTracks);
//================================================================
//               apo Eftihi 
                  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler =
(AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }

 Double_t nsigmall = 100.0;
       Double_t nsigma = 100.0;
 //      Double_t nsigmaPion =-100.0;
//       Double_t nsigmaPi=-100.0;

//    
            Float_t nprimTrk =0.;
// track loop
   for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {

    AliESDtrack* track = esd->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    fHistPt->Fill(track->Pt());


     //   edw ypologizontai ta sigmas
 //   nsigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));



//    sigmas    
             nsigmall  = (fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));
//    nsigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));
 //   nsigmall = (fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));
       nsigma  = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon));

//==================================
  Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
  Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
  if (track->GetTPCNclsF()>0) {
    ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
    fRatioCrossedRows->Fill(ratioCrossedRowsOverFindableClustersTPC);
   }

    fHistPt->Fill(track->Pt());



  //    Int_t tpcNClMoth = track->GetTPCclustersIter1(0);  // 28/4/2010
      Int_t tpcNCl = track->GetTPCclusters(0);  // 25/1/2010
      Double_t tpcSign = track->GetSign();  // 25/1/2010
    
    Int_t label = track->GetLabel();
    label = TMath::Abs(label);


    UInt_t status=track->GetStatus();

    if((status&AliESDtrack::kITSrefit)==0) continue;   
    if((status&AliESDtrack::kTPCrefit)==0) continue;
      if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;  // 

      Double_t extCovPos[15];
      track->GetExternalCovariance(extCovPos);    


    track->GetXYZ(vtrack);
 fXvYv->Fill(vtrack[0],vtrack[1]);  
 fZvYv->Fill(vtrack[0],vtrack[2]);  
 fZvXv->Fill(vtrack[1],vtrack[2]);  

// track momentum
     track->GetPxPyPz(ptrack);
    
    TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
    
          Double_t   etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  );  // kaon energy
         Double_t rapiditK = 0.5 * (TMath::Log(  (etracK + ptrack[2]  ) / ( etracK - ptrack[2])  ))  ; // kaon rapidity 
    
    Double_t trackEta=trackMom.Eta();
//     Double_t trMoment=trackMom.Mag();       
    Double_t trackPt = track->Pt();
    
    
      
    Float_t bpos[2];
    Float_t bCovpos[3];
    track->GetImpactParameters(bpos,bCovpos);
    
    if (bCovpos[0]<=0 || bCovpos[2]<=0) {
     Printf("Estimated b resolution lower or equal zero!");
     bCovpos[0]=0; bCovpos[2]=0;
    }

    Float_t dcaToVertexXYpos = bpos[0];
//    Float_t dcaToVertexZpos = bpos[1];
    
 //   fRpr->Fill(dcaToVertexZpos);
   fRpr      ->Fill(dcaToVertexXYpos);
  // ===   if((dcaToVertexXYpos>0.3)||(dcaToVertexZpos>2.5)) continue;   //    allagi-dokini    3/6                 
    
                       if (!fMaxDCAtoVtxCut->AcceptTrack(track)) continue;


//           count nprimTrk ????

            nprimTrk++ ;
//
 //  cut on eta 
 // xwris tea cut 25/3         if(  (TMath::Abs(trackEta )) > 0.9 ) continue;

    fHistPtESD->Fill(track->Pt());

   // Add Kink analysis
   
   	    	Int_t indexKinkPos=track->GetKinkIndex(0);

              if ( indexKinkPos >0 )   fHistPtESD->Fill(track->Pt());
//  loop on kinks
		if(indexKinkPos<0){     ////mother kink
               fPtKink->Fill(track->Pt()); /// pt from track

	// select kink class	

	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
  fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);

//
         // DCA kink
          Double_t  Dist2 = kink->GetDistance();
          fDCAkink->Fill( Dist2   );
//
    Int_t ESDLabeld = kink->GetLabel(1);
//
	
      
	   const TVector3 vposKink(kink->GetPosition());
 fPosiKink ->Fill( vposKink[0], vposKink[1]  );

//                24/3  
//                  if     ( vposKink[2] < -160. ) continue;   // cut for background 
//
//   Double_t  lengthK = TMath::Sqrt( vposKink[0]*vposKink[0] + vposKink[1]*vposKink[1] + vposKink[2]*vposKink[2] ) ;
   // Double_t dxKink = vpos[0]-vposKink[0], dyKink=vpos[1]-vposKink[1], dzKink=vpos[2]-vposKink[2]; 
   Double_t  dzKink=vpos[2]-vposKink[2]; 
   // Double_t lifeKink= TMath::Sqrt( dxKink*dxKink + dyKink*dyKink + dzKink*dzKink ) ;
//
            Double_t tanLamda = track->GetTgl();  // 25/6/2010

   // Double_t lifeKink= (TMath::Abs( dzKink ))*( TMath::Sqrt(1.+ tanLamda*tanLamda) ) / tanLamda ;
   Double_t lifeKink= (TMath::Abs( dzKink ))*( TMath::Sqrt(1.+ tanLamda*tanLamda) ) / (TMath::Abs( tanLamda)) ;
	   const TVector3 motherMfromKink(kink->GetMotherP());
	   const TVector3 daughterMKink(kink->GetDaughterP());

	   Float_t qT=kink->GetQt();
       //      Float_t motherPt=motherMfromKink.Pt();
       //     Float_t etaMother=motherMfromKink.Eta();

           fHistQtAll->Fill(qT) ;  //  Qt   distr
                  
         //  fptKink->Fill(motherMfromKink.Pt()); /// pt from kink

           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

//     Rapidity Cut
         // if(  (TMath::Abs(etaMother)) > 0.9 ) continue;
          if(  (TMath::Abs(rapiditK )) > 0.5 ) continue;
//                Pt cut 
        if ( (track->Pt())<.200)continue;

//       fake kinks, small Qt and small kink angle
                fQtMothP->Fill( track->P(), qT);

//           if(qT<0.012) continue;  // remove fake kinks
    if ( qT> 0.04)  fHistQt1  ->Fill(qT) ;  //  Qt   distr    , without pions and BG


            fHistEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
      fHistQt2->Fill(qT);  // PDG ESD kaons            

//          maximum decay angle at a given mother momentum
	   Double_t maxDecAngKmu=f1->Eval(track->P()          ,0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(    track->P(),       0.,0.,0.);
//
// 31/1/010         if(  motherMfromKink.Mag() < daughterMKink.Mag() ) continue;
//remove the double taracks 
         //   28/4/2010  if( (kinkAngle<1.)  ) continue;
         if( (kinkAngle<2.)  ) continue;
    //       if(qT<0.012) continue;  // remove fake kinks
           //  BG  ?????==============
              if ( TMath::Abs(vposKink[2]) >  225. ) continue ;
         //     if ( TMath::Abs(vposKink[2]) <  0.5 ) continue ;

            fKinkKaonBg->Fill(track->Pt());     
 fAngMomPi->Fill( track->P(),           kinkAngle); 
         //  fRadiusNclAll->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
       //   2/12  test   fRadiusNclAll->Fill( (kink->GetR()) ,(track->GetNcls(1)         ) ) ;
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
   //      Float_t ptKink=TMath::Sqrt(p1XM*p1XM + p1YM*p1YM);
  
      if( ( kink->GetR()> 110 ) && ( kink->GetR() < 220 )  )  {
      if (qT>0.04)  fAngMomKC->Fill(track->P(), kinkAngle); 
      // 2/10  if (qT>0.12)  fAngMomKC->Fill(track->P(), kinkAngle); 
          if ( qT>0.04) fM1kaon->Fill(invariantMassKmu);
             // if ( qT > 0.12) 
             if ( qT > 0.04) 
           fRadiusNclAll->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ; // mother's clusters 
         // fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
//    try for daughters ncl
            if ( label = ESDLabeld )  fRadiusNcl->Fill( (kink->GetR()) ,(track->GetNcls(1)     )   ); // Dau's clusters;
  }    
          if(  ( tpcNCl<20) ) continue;
      Int_t tpcNClHigh = -51.67+ (11./12.)  *( kink->GetR() ) ;
               if ( tpcNCl > tpcNClHigh) continue;

      Int_t tpcNClMin  = -85.5 + (65./95.)  *( kink->GetR() ) ;
               if ( tpcNCl < tpcNClMin ) continue;

//         if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) > 0.63 ) continue;
  //          if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) < 0.20 ) continue;

//
               fHistPtKPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
         fRadiusNclK->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;

         //fM1kaon->Fill(invariantMassKmu);

          //     frad->Fill(kink->GetR());  // kink 
//  kaon selection from kinks
 // if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=200.))&&(TMath::Abs(etaMother)<0.9)&&(invariantMassKmu<0.6)){
//---------------------------------
    if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>=110.)&&(kink->GetR()<=220.))&&(TMath::Abs(rapiditK)<0.5)&&(invariantMassKmu<0.8)){
//-----------------------
 // 16/10  if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=210.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.6)){
//  if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=133.)&&(kink->GetR()<=179.))&&(TMath::Abs(rapiditK)<0.5)&&(invariantMassKmu<0.6)){      // STAR 

        if( (kinkAngle<maxDecAngpimu*1.1) ) continue; 
                 if ( (kinkAngle>maxDecAngKmu*.98) && ( track->P() >1.1 )) continue;  ///5/5/2010
              //    fAngMomKC->Fill(motherMfromKink.Mag(), kinkAngle);

      //   fAngMomKC->Fill(track->P(), kinkAngle); 

      //   fTPCSgnlPa->Fill( trMoment ,(track->GetTPCsignal()  ) ) ;
                fTPCSgnlPa->Fill( track->GetInnerParam()->GetP() ,(track->GetTPCsignal()  ) ) ;
//

                          if ( nsigma               > 3.5) continue;
                    //      if ( nsigma               > 4.0) continue;
     //  fMultESDK->Fill(nGoodTracks);
                                     fHistPtKaon->Fill(track->Pt());   //all PID kink-kaon
              if(tpcSign >0.)        fHistPtKaoP->Fill( track->Pt()         ) ;   //all PID kink-kaon
               if ( tpcSign <0.)    fHistPtKaoN->Fill(  track->Pt()        ) ;   //all PID kink-kaon
//  percentile centrality
         if ( percent    >= .0  && percent < 5.)     fKinkKPt05 -> Fill ( track->Pt() ); // percent 0-5 
         if ( percent    >= 5.  && percent < 10.) fKinkKPt510->Fill( track->Pt() );  // percent  5-10 
        if ( percent    >= 10. && percent < 20.)  fKinkKPt1020->Fill( track->Pt() );  // percent 10-20
         if ( percent    >= 20. && percent < 30.) fKinkKPt2030->Fill( track->Pt() );  // percent 20-30 
       if ( percent    >= 30. && percent < 40.) fKinkKPt3040->Fill( track->Pt() );  // percent 30-40 
       if ( percent    >= 40. && percent < 50.) fKinkKPt4050->Fill( track->Pt() );  // percent 40-50  
       if ( percent    >= 50. && percent < 60.) fKinkKPt5060->Fill( track->Pt() );  // percent 50-60 
       if ( percent    >= 60. && percent < 70.) fKinkKPt6070->Fill( track->Pt() );  // percent 60-70 
        if ( percent    >= 70. && percent < 80.) fKinkKPt7080->Fill( track->Pt() );  // percent 70-80 
        if ( percent    >= 80. && percent < 90.) fKinkKPt8090->Fill( track->Pt() );  // percent 80-90 
//     
                frad->Fill(kink->GetR());  // kink 
               //frad->Fill(lengthK     );  // kink 
           //    fradLK->Fill(lifeKink    );  // kink 
             fHistEtaK->Fill(trackEta);
            frapiKESD ->Fill(rapiditK);  //  rapidityof kaons 

                     Float_t signPt= tpcSign*trackPt;
                  //fSignPtNcl->Fill( signPt  ,   tpcNCl   );
                  fSignPtNcl->Fill( signPt  ,   tpcNCl   );   ///  28/4/2010
                  fSignPtEta->Fill( signPt  , rapiditK  );
                  fEtaNcl->Fill( rapiditK, tpcNCl    );
                  fSignPt->Fill( signPt );
                  fChi2NclTPC->Fill( (track->GetTPCchi2() ) ,  tpcNCl );
         fRatChi2Ncl-> Fill (  (track->GetTPCchi2()/track->GetTPCclusters(0)  )) ;
    fdcatoVxXY->Fill(dcaToVertexXYpos);
      //if( ( (trMoment<0.4)&&(track->GetTPCsignal()<60) )||( (trMoment<0.6 )&&(track->GetTPCsignal()<50 )) )
        // fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
     // if( ( (trMoment<0.4)&&(track->GetTPCsignal()<60) )||( (trMoment<0.6 )&&(track->GetTPCsignal()<50 )) )
                  fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
            //      fEtaNcl->Fill( trackEta, tpcNCl    );
            //  Double_t tpcdedx=100000.* TMath::Abs( track->GetTPCsignal());
        //  fTPCSgnlP->Fill( track->P() ,tpcdedx  ) ;
  //       fTPCSgnlP->Fill(track->P(), (track->GetTPCsignal()  ) ) ;
          fTPCSgnlP->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()  ) ) ;
         fRadiusNclClean ->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
//   
		        flifetime->Fill(( lifeKink*.493667   )  /track->P()   ) ;
		// 15/7       flifetime->Fill(( lifeKink*.493667   )  /track->Pt()   ) ;
// 
             fKinkKaon->Fill(track->Pt());        
// 

             fKinKRbn->Fill(track->Pt());        
             fptKMC   ->Fill(  track->Pt()    );        
              fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK); 
               fAngMomK->Fill(    track->P(),        kinkAngle); 
       fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[0], vposKink[2]  );
          fPosiKinKYZ->Fill( vposKink[1], vposKink[2]  );

        }  //  kink selection 
                  

	}  //End Kink Information    
  

  } //track loop 
    //   fMultESDK->Fill(nprimTrk);

//         }   // close persent   11 may 2011 

  PostData(1, fListOfHistos);

}      

//________________________________________________________________________
void AliAnalysisKinkESDatION::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}
//*

const AliESDVertex* AliAnalysisKinkESDatION::GetEventVertex(const AliESDEvent* esd) const
  // Get the vertex from the ESD and returns it if the vertex is valid
  
{
  // Get the vertex 
  
//  const AliESDVertex* vertex = esd->GetPrimaryVertex();
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
