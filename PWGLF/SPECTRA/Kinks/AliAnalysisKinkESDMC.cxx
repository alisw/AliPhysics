/**************************************************************************
 * Authors: Martha Spyropoulou-Stassinaki and the  members 
 * of the Greek group at Physics Department of Athens University
 * Paraskevi Ganoti, Anastasia Belogianni and Filimon Roukoutakis 
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDMC class
//       Example of an analysis task for kink topology study
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
#include "AliTrackReference.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisKinkESDMC.h"
#include "AliStack.h"
#include "AliESDpid.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliPhysicsSelectionTask.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
 #include "AliPIDResponse.h"
///#include "AliTPCpidESD.h"

ClassImp(AliAnalysisKinkESDMC)


//________________________________________________________________________
AliAnalysisKinkESDMC::AliAnalysisKinkESDMC(const char *name) 
  : AliAnalysisTaskSE(name), fHistPtESD(0),fHistPt(0),fHistQtAll(0),fHistQt1(0),fHistQt2(0)
  , fHistPtKaon(0),fHistPtKPDG(0),fHistEta(0),fHistEtaK(0),fptKMC(0),fMultiplMC(0),fESDMult(0),frad(0),
  fradMC(0), fKinkKaon(0), fKinkKaonBg(0), fM1kaon(0),  fgenPtEtR(0),fPtKink(0),  
   fcodeH(0), fdcodeH(0), fAngMomK(0),fAngMomPi(0), fAngMomKC(0),  fMultESDK(0), fMultMCK(0),
 fSignPtNcl(0), fSignPtEta(0), fSignPtEtaMC(0), fSignPtMC(0),  fEtaNcl(0), fSignPt(0), fChi2NclTPC(0), fRatChi2Ncl(0),
  fRadiusNcl(0),  fTPCSgnlP(0), fTPCSgnlPa(0),  fSignPtGen(0),
  fRpr(0),fZpr(0), fdcatoVxXY(0),    fMCEtaKaon(0),
 fZvXv(0),fZvYv(0),fXvYv(0),fPtPrKink(0),fgenPtEtRP(0),fgenPtEtRN(0),fkinkKaonP(0),fkinkKaonN(0),
   frapidESDK(0), frapidKMC(0),  fPtKPlMC(0), fPtKMnMC(0),  
     fHistPtKaoP(0), fHistPtKaoN(0), fHiPtKPDGP(0), fHiPtKPDGN(0),fKinKBGP(0),fKinKBGN(0), 
    fQtKMu(0),fQtKPi(0),fQtKEl(0),fFakepipi(0), fFakeKPi(0),
     fDCAkink(0), fDCAkinkBG(0), fPosiKink(0),  fPosiKinkK(0), fPosiKinKXZ(0), fPosiKinKYZ(0),  fPosiKinKBgZY(0), 
    fcode2(0), fcode4(0), fZkinkZDau(0),  
         fQtKMuMC(0),  fQtKElMC(0), fQtKPiMC(0),   fQtK3PiP(0),fQtK3PiM(0),  fmaxAngMomKmu(0),
           fPosiKinKBgZX(0), fPosiKinKBgXY(0),  fMinvPi(0),fMinvKa(0),fMinvPr(0),
                        fTPCSgnlPtpc(0),
       fTPCMomNSgnl(0),  fMothKinkMomSgnl(0), fNSigmTPC(0),  fTPCSgnlKinkDau(0),fcodeDau1(0),fcodeDau2(0), fMothKinkMomSgnlD(0), 
     fInvMassMuNuAll(0),   fInvMassMuNuPt(0), fRatioCrossedRows(0), fRatioCrossedRowsKink(0), fRadiusPt(0), fRadiusPtcln(0),
     fPtCut1(0), fPtCut2(0), fPtCut3(0),  fAngMomKKinks(0),     
  flengthMCK(0), flifetiMCK(0), flifetim2(0), fLHelESDK(0),flifeInt(0), flifeYuri(0), flenYuri(0), flenTrRef(0),flifeSmall(0), flifetime(0),flifTiESDK(0),  
    flifeKink(), flenHelx(0), fradPtRapMC(0), fradPtRapDC(0), fradPtRapESD(0), fRadNclcln(0),
    f1(0), f2(0),
  fListOfHistos(0),fLowMulcut(-1),fUpMulcut(-1), fKinkRadUp(200),fKinkRadLow(130), fCutsMul(0),fMaxDCAtoVtxCut(0), fPIDResponse(0)

{
  // Constructor

  // Define input and output slots here

	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkESDMC::UserCreateOutputObjects() 
{
    
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());liESDtrackCuts("Mul","Mul");
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
//
                   fMaxDCAtoVtxCut=new AliESDtrackCuts("fMaxDCA", "fMaxDCA");
       fMaxDCAtoVtxCut->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");
    fMaxDCAtoVtxCut->SetMaxChi2TPCConstrainedGlobal(36);

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
//
      Double_t gPt7K0[45] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.2, 2.4, 2.6, 2.8,  3.0,   3.3, 3.6, 3.9,
                         4.2, 4.6,5.0, 5.4, 5.9,  6.5,   7.0,7.5, 8.0,8.5,  9.2, 10., 11., 12., 13.5,15.0 };  // David K0
/*
   Double_t gPt7TOF[47] = { 0.2,0.25, 0.3,0.35,  0.4,0.45,  0.5,0.55,  0.6,0.65,  0.7,0.75,  0.8, 0.85, 0.9, 0.95, 1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0,
                         3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8,5.0 };  //  Barbara TOF Kch
*/
	fHistPtESD = new TH1F("fHistPtESD", "P_{T} distribution",200, 0.0,10.0);
	fHistPt = new TH1F("fHistPt", "P_{T} distribution",200, 0.0,10.0); 
	fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300); 
	fHistQt1= new TH1F("fHistQt1", "Q_{T} distribution",100, 0.0,.300); 
	fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300); 
//	fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",200, 0.0,10.0); 
	fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",44,gPt7K0); 
	fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",44, gPt7K0   ); 
	fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3); 
	fHistEtaK= new TH1F("fHistEtaK", "EtaK distribution", 26,-1.3, 1.3); 
	fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",44,  gPt7K0  ); 
	fMultiplMC= new TH1F("fMultiplMC", " charge particle multipl",100, 0.0,2500.);
	fESDMult= new TH1F("fESDMult", "charge multipliESD",100, 0.0,100.); 
	frad= new TH1F("frad", "radius  K ESD recon",100,0.,1000.); 
	fradMC= new TH1F("fradMC", "radius  K generated",100,0.,1000.); 
	fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi", 44, gPt7K0  ); 
	fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",44 , gPt7K0  ); 
	fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",180,0.1, 1.0); 
	fgenPtEtR= new TH1F("fgenPtEtR", "P_{T}Kaon distribution", 44, gPt7K0  ); 
	fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  bution",44, gPt7K0   ); 
	fcodeH   = new TH2F("fcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fdcodeH = new TH2F("fdcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,10.0,80,0.,80.);
	fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,10.0,80,0.,80.);
	fAngMomKC= new TH2F("fAngMomKC","Decay angle vrs Mother Mom,K",100,0.0,10.0,80,0.,80.);
	fMultESDK=new TH1F("fMultESDK", "charge multipliESD kaons",100, 0.,2500.); 
	fMultMCK=new TH1F("fMultMCK", "charge multipli MC kaons",100, 0.,2500.); 
	fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",80,-4.,4.0,70,20.,160.);
	fSignPtEta= new TH2F("fSignPtEta","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
	fSignPtEtaMC= new TH2F("fSignPtEtaMC","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
	fSignPtMC= new TH1F("fSignPtMC","SignPt ,K",100,-5.,5.0);
	fEtaNcl= new TH2F("fEtaNcl","Eta vrs Ncl,K",26,-1.3,1.3,70,20.,160.);
	fSignPt= new TH1F("fSignPt","SignPt ,K",100,-5.,5.0);
	fChi2NclTPC= new TH2F("fChi2NclTPC","Chi2vrs nclust,K",100,0.,500., 70,20, 160);
	fRatChi2Ncl= new TH1F("fRatChi2Ncl","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
	fRadiusNcl = new TH2F("fRadiusNcl","KinkRadius Ncl,K",75,100.,250., 80,0, 160);
	fTPCSgnlP = new TH2F("fTPCSgnlP","Kink TCP de/dx,K",300,0.,15.,150,0.,300);
	fTPCSgnlPa= new TH2F("fTPCSgnlPa","Kink TCP de/dx,a",300,0.,15.,150,0.,300);
	fSignPtGen= new TH1F("fSignPtGen","SignPtGen ,K",100,-5.0,5.0);
	fRpr = new TH1D("fRpr", "rad distribution  PID pr",100,-1.0,1.0);
	fZpr = new TH1D("fZpr", "z distribution PID pr  ",60,-15.,15.);
	fdcatoVxXY = new TH1D("fdcatoVxXY", "dca  distribution PID  ",20,-1.,1.);
	fMCEtaKaon = new TH1F("fMCEtaKaon"," Hist of Eta K -Kink Selecied",26,-1.3,1.3);
	fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
	fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
	fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
	fPtPrKink=new TH1F("fPtPrKink","pt of ESD  kaonKink tracks",300, 0.0,15.0);
	fgenPtEtRP= new TH1F("fgenPtEtRP", "P_{T}Kaon distribution positive", 44, gPt7K0  ); 
	fgenPtEtRN= new TH1F("fgenPtEtRN", "P_{T}Kaon distribution negative", 44, gPt7K0  ); 
	fkinkKaonP= new TH1F("fKinkKaonP", "P_{T}Kaon distribution positive", 44, gPt7K0  ); 
	fkinkKaonN= new TH1F("fKinkKaonN", "P_{T}Kaon distribution negative", 44, gPt7K0  ); 
	frapidESDK= new TH1F("frapidESDK", "rapidity distribution", 26,-1.3, 1.3); 
	frapidKMC = new TH1F("frapidKMC ", "rapidity distri  MC  ",26,-1.3, 1.3); 
	fPtKPlMC= new TH1F("fPtKPlMC", "P_{T}Kaon Pos  generated", 44, gPt7K0  ); 
	fPtKMnMC= new TH1F("fPtKMnMC", "P_{T}Kaon Minusgenerated",44 , gPt7K0  ); 
	fHistPtKaoP= new TH1F("fHistPtKaoP", "P_{T}Kaon Pos ESD", 44, gPt7K0  ); 
	fHistPtKaoN= new TH1F("fHistPtKaoN", "P_{T}Kaon Neg ESD", 44, gPt7K0  ); 
	fHiPtKPDGP= new TH1F("fHiPtKPDGP", "P_{T}Kaon Pos ESD", 44,  gPt7K0 ); 
	fHiPtKPDGN= new TH1F("fHiPtKPDGN", "P_{T}Kaon neg ESD", 44, gPt7K0  ); 
	fKinKBGP  = new TH1F("fKinKBGP  ", "P_{T}Kaon Pos ESD", 44, gPt7K0  ); 
	fKinKBGN= new TH1F("fKinKBGN", "P_{T}Kaon neg ESD", 44, gPt7K0  ); 
	fQtKMu= new TH1F("fQtKMu", "Q_{T} distribution  K to mu ",100, 0.0,.300); 
	fQtKPi= new TH1F("fQtKPi", "Q_{T} distribution K to pi",100, 0.0,.300); 
	fQtKEl= new TH1F("fQtKEl", "Q_{T} distribution   K to elec",100, 0.0,.300); 
	fFakepipi = new TH1F("fFakepipi", "P_{T}fake pipi   ",44 , gPt7K0  ); 
	fFakeKPi = new TH1F("fFakeKPi", "P_{T}fake Kpi   ", 44, gPt7K0  ); 
	fDCAkink = new TH1F("fDCAkink", "DCA kink vetrex ",50, 0.0,1.0); 
	fDCAkinkBG = new TH1F("fDCAkinkBG", "DCA kink vetrex ",50, 0.0,1.0); 
	fPosiKink= new TH2F("fPosiKink", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.); 
	fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.); 
	fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.); 
	fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.); 
	fPosiKinKBgZY= new TH2F("fPosiKinKBgZY", "Y vrx Z kink Vrexbg ",100, -300.0,300.0,100, -300, 300.); 
	fcode2   = new TH2F("fcode2", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fcode4   = new TH2F("fcode4", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
	fZkinkZDau = new TH2F("fZkinkZDau", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.); 
  fQtKMuMC= new TH1F("fQtKMuMC", "Q_{T} distribution  K to mu MC",100, 0.0,.300); 
  fQtKPiMC= new TH1F("fQtKPiMC", "Q_{T} distribution K to pi MC",100, 0.0,.300); 
  fQtKElMC= new TH1F("fQtKElMC", "Q_{T} distribution   K to elec MC",100, 0.0,.300); 
  fQtK3PiP= new TH1F("fQtK3PiP", "Q_{T} distribution K to 3pi ",100, 0.0,.300); 
  fQtK3PiM= new TH1F("fQtK3PiM", "Q_{T} distribution K to 3pi ",100, 0.0,.300); 
  fmaxAngMomKmu= new TH2F("fmaxAngMomKmu","Decay angle vrs Mother Mom,Kmu",100,0.0,10.0,80,0.,80.);
	fPosiKinKBgZX= new TH2F("fPosiKinKBgZX", "X vrx Z kink Vrexbg ",100, -20.0,20.0,100, 0., 300.); 
	fPosiKinKBgXY= new TH2F("fPosiKinKBgXY", "Y vrx X kink Vrexbg ",100, -300.0,300.0,100, -300, 300.); 
	fMinvPi= new TH1F("fMinvPi","Invar m(kaon) from kink-> decay",100,0.0, 1.2); 
	fMinvKa= new TH1F("fMinvKa","Invar m(kaon) from kink-> decay",100,0.0, 2.0); 
	fMinvPr= new TH1F("fMinvPr","Invar m(kaon) from kink-> decay",100,0.0, 1.2); 
                fTPCSgnlPtpc = new TH2F("fTPCSgnlPtpc","TPC signal de/dx Mom TPC,K  ",100,0.0,8.0,100, 0., 250.    );
    fTPCMomNSgnl = new TH2F("fTPCMomNsgnl","TPC signal de/dx Mom TPC,K  ",100,0.0,8.0,20 , -10., 10.);
    fMothKinkMomSgnl  = new TH2F("fMothKinkMomSgnl","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.,100, 0., 250.    );
    fNSigmTPC    = new TH1F("fNSigmTPC","TPC Nsigma  de/dx  TPC,K  ", 30 , -7.5, 7.5);
    fTPCSgnlKinkDau = new TH2F("fTPCSgnlKinkDau","TPC signal de/dx Mom,K",100,0.0,8.0,100,0.,250.);
	fcodeDau1   = new TH2F("fcodeDau1", "code vrs dcode dist. kinks,K",100,0.,500.,100,0.,500.);
	fcodeDau2   = new TH2F("fcodeDau2", "code vrs dcode dist. kinks,K",100,0.,500.,100,0.,500.);
    fMothKinkMomSgnlD = new TH2F("fMothKinkMomSgnlD","TPC signal de/dx Mom TPC,Kink  ",100,0.0,250.,100, 0., 250.    );
	fInvMassMuNuAll = new TH1F("fInvMassMuNuAll","Invar from kink->mu+netrino decay",180,0.1, 1.0); 
	fInvMassMuNuPt  = new TH2F("fInvMassMuNuPt","Invar from kink->mu+netrino decay vs Pt",180,0.1, 1.0, 100, 0. , 10.); 
      fRatioCrossedRows = new TH1F("fRatioCrossedRows","Ratio crossed rows  in TPC",20,0.0,1.0);
  fRatioCrossedRowsKink = new TH1F("fRatioCrossedRowsKink","Ratio crossed rows  in TPC for kinks",20,0.0,1.0);
  fRadiusPt =new TH2F("fRadiusPt","radius vs pt  ",80, 90.,250.,100, 0.,10.              );
  fRadiusPtcln =new TH2F("fRadiusPtcln","radius vs pt clean ",80, 90.,250.,100, 0.,10.              );
  fPtCut1 = new TH1F("fPtCut1", "P_{T}Kaon distribution",300, 0.0,15.0);
  fPtCut2 = new TH1F("fPtCut2", "P_{T}Kaon distribution",300, 0.0,15.0);
  fPtCut3 = new TH1F("fPtCut3", "P_{T}Kaon distribution",300, 0.0,15.0);
  fAngMomKKinks = new TH2F("fAngMomKKinks","Decay angle vrs Mother Mom,Kinks",300,0.0,15.0,100,0.,100.);
   flengthMCK=new TH1F("flengthMCK", "length of K  MCref decay ",100,0.,1000.0);
  flifetiMCK=new TH1F("flifetiMCK", "lifetime ref K   Decay   ",100,0.,1000.0);
  flifetim2 =new TH1F("flifetim2", "lifetime ref K   Decay   ",100,0.,1000.0);
  fLHelESDK =new TH1F("fLHelESDK", "lifetime ref K   Decay   ",100,0.,1000.0);
  flifeInt =new TH1F("flifeInt", "lifetime ref K   Decay   ",100,0.,1000.0);
  flifeYuri=new TH1F("flifeYuri","lifetime ref K   Decay   ",100,0.,1000.0);
  flenYuri=new TH1F("flenYuri","lifetime ref K   Decay   ",100,0.,1000.0);
  flenTrRef =new TH1F("flenTrRef","lifetime ref K   Decay   ",100,0.,1000.0);
  flifeSmall=new TH1F("flifeSmall","lifetime ref K   Decay   ",100,0.,1000.0);
  flifetime =new TH1F("flifetime","lifetime ref K   Decay   ",100,0.,1000.0);
  flifTiESDK=new TH1F("flifTiESDK","lifetime ref K   Decay   ",100,0.,1000.0);
  flifeKink =new TH1F("flifeKink", "lifetime ref K   Decay   ",100,0.,1000.0);
  flenHelx  =new TH1F("flenHelx", "lifetime ref K   Decay   ",100,0.,1000.0);
     fradPtRapMC=new TH3F("fradPtRapMC","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
     fradPtRapDC=new TH3F("fradPtRapDC","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
     fradPtRapESD=new TH3F("fradPtRapESD","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
	fRadNclcln = new TH2F("fRadNclcln","KinkRadius Ncl,K",75,100.,250., 80,0, 160);




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
   fListOfHistos->Add(fradMC);
   fListOfHistos->Add(fKinkKaon);
   fListOfHistos->Add(fKinkKaonBg);
   fListOfHistos->Add(fM1kaon);
   fListOfHistos->Add(fgenPtEtR);
   fListOfHistos->Add(fPtKink);
   fListOfHistos->Add(fcodeH);
   fListOfHistos->Add(fdcodeH);
   fListOfHistos->Add(fAngMomK);
   fListOfHistos->Add(fAngMomPi);
   fListOfHistos->Add(fAngMomKC);
   fListOfHistos->Add(fMultESDK);
   fListOfHistos->Add(fMultMCK);
   fListOfHistos->Add(fSignPtNcl);
   fListOfHistos->Add(fSignPtEta);
   fListOfHistos->Add(fSignPtEtaMC);
   fListOfHistos->Add(fSignPtMC);
   fListOfHistos->Add(fEtaNcl);
   fListOfHistos->Add(fSignPt);
   fListOfHistos->Add(fChi2NclTPC);
   fListOfHistos->Add(fRatChi2Ncl); 
   fListOfHistos->Add(fRadiusNcl); 
   fListOfHistos->Add(fTPCSgnlP); 
   fListOfHistos->Add(fTPCSgnlPa); 
   fListOfHistos->Add(fSignPtGen);
   fListOfHistos->Add(fRpr);
   fListOfHistos->Add(fZpr);
   fListOfHistos->Add(fdcatoVxXY);
   fListOfHistos->Add(fMCEtaKaon);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fPtPrKink);
   fListOfHistos->Add(fgenPtEtRP);
   fListOfHistos->Add(fgenPtEtRN);
   fListOfHistos->Add(fkinkKaonP);
   fListOfHistos->Add(fkinkKaonN);
   fListOfHistos->Add(frapidESDK);
   fListOfHistos->Add(frapidKMC);
   fListOfHistos->Add(fPtKPlMC);
   fListOfHistos->Add(fPtKMnMC);
   fListOfHistos->Add(fHistPtKaoP);
   fListOfHistos->Add(fHistPtKaoN);
   fListOfHistos->Add(fHiPtKPDGP);
   fListOfHistos->Add(fHiPtKPDGN);
   fListOfHistos->Add(fKinKBGP);
   fListOfHistos->Add(fKinKBGN);
   fListOfHistos->Add(fQtKMu);
   fListOfHistos->Add(fQtKPi);
   fListOfHistos->Add(fQtKEl);
   fListOfHistos->Add(fFakepipi);
   fListOfHistos->Add(fFakeKPi);
   fListOfHistos->Add(fDCAkink);
   fListOfHistos->Add(fDCAkinkBG);
   fListOfHistos->Add(fPosiKink);
   fListOfHistos->Add(fPosiKinkK);
   fListOfHistos->Add(fPosiKinKXZ);
   fListOfHistos->Add(fPosiKinKYZ);
   fListOfHistos->Add(fPosiKinKBgZY);
   fListOfHistos->Add(fcode2);
   fListOfHistos->Add(fcode4);
   fListOfHistos->Add(fZkinkZDau);
   fListOfHistos->Add(fQtKMuMC);
   fListOfHistos->Add(fQtKPiMC);
   fListOfHistos->Add(fQtKElMC);
   fListOfHistos->Add(fQtK3PiP);
   fListOfHistos->Add(fQtK3PiM);
   fListOfHistos->Add(fmaxAngMomKmu);
   fListOfHistos->Add(fPosiKinKBgZX);
   fListOfHistos->Add(fPosiKinKBgXY);
   fListOfHistos->Add(fMinvPi);
   fListOfHistos->Add(fMinvKa);
   fListOfHistos->Add(fMinvPr);
   fListOfHistos->Add(fTPCSgnlPtpc);
   fListOfHistos->Add(fTPCMomNSgnl);
   fListOfHistos->Add(fMothKinkMomSgnl);
   fListOfHistos->Add(fNSigmTPC);
   fListOfHistos->Add(fTPCSgnlKinkDau);
   fListOfHistos->Add(fcodeDau1);
   fListOfHistos->Add(fcodeDau2);
   fListOfHistos->Add(fMothKinkMomSgnlD);
   fListOfHistos->Add(fInvMassMuNuAll);
   fListOfHistos->Add(fInvMassMuNuPt);
   fListOfHistos->Add(fRatioCrossedRows);
   fListOfHistos->Add(fRatioCrossedRowsKink);
   fListOfHistos->Add(fRadiusPt);
   fListOfHistos->Add(fRadiusPtcln);
   fListOfHistos->Add(fPtCut1);
   fListOfHistos->Add(fPtCut2);
   fListOfHistos->Add(fPtCut3);
   fListOfHistos->Add(fAngMomKKinks);
   fListOfHistos->Add(flengthMCK);
   fListOfHistos->Add(flifetiMCK);
   fListOfHistos->Add(flifetim2);
   fListOfHistos->Add(fLHelESDK);
   fListOfHistos->Add(flifeInt);
   fListOfHistos->Add(flifeYuri);
   fListOfHistos->Add(flenYuri);
   fListOfHistos->Add(flenTrRef);
   fListOfHistos->Add(flifeSmall); 
   fListOfHistos->Add(flifetime); 
   fListOfHistos->Add(flifTiESDK); 
   fListOfHistos->Add(flifeKink); 
   fListOfHistos->Add(flenHelx); 
   fListOfHistos->Add(fradPtRapMC); 
   fListOfHistos->Add(fradPtRapDC); 
   fListOfHistos->Add(fradPtRapESD); 
   fListOfHistos->Add(fRadNclcln); 

  PostData(1, fListOfHistos);
}

//________________________________________________________________________
void AliAnalysisKinkESDMC::UserExec(Option_t *) 
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
     fMultMCK->Fill(mcEvent->GetNumberOfTracks() );
//


//   multiplicity selection 
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
 for (Int_t iMc = 0; iMc < mcEvent->GetNumberOfTracks(); iMc++)
  {

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
    //  AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }
// keep only primaries
  if (!stack->IsPhysicalPrimary(iMc)) continue;

 //   

            Double_t ptK = particle->Pt();

     if( ptK <0.200) continue;       //    12/2/2012
//

                Float_t charg=0;
      Float_t code = particle->GetPdgCode();
            Int_t  mcProcess=-1011;
//---------------------------------------kaon selection 
      if ((code==321)||(code==-321)){
	    
    
          Double_t   etracKMC= TMath::Sqrt(particle->P() *particle->P()  + 0.493677 *0.493677  );
         Double_t rapidiKMC = 0.5 * (TMath::Log(  (etracKMC +particle->Pz())/( etracKMC-particle->Pz() )) )  ;

     if ( TMath::Abs( rapidiKMC) > 0.7) continue;   // 
            frapidKMC ->Fill(rapidiKMC) ;  //18/feb rapiddistr of PDG kink ESD  kaons
 
	    
//  maximum angle    vrs momentum
  //   Double_t maxDecAnKmu=f1->Eval(particle->P(),      0.,0.,0.);
 //       fmaxAngMomKmu->Fill(particle->P() , maxDecAnKmu);

               if(code > 0 ) charg =1;
               if(code < 0 ) charg =-1;
                         Float_t chargPt= charg*ptK;

	      fptKMC->Fill(ptK);
         fSignPtGen->Fill(chargPt);// kaon gensign pt
	      if (charg==1 )  fPtKPlMC->Fill( ptK );
	       if ( charg==-1) fPtKMnMC->Fill( ptK  );
// primary   vertex
            //      Double_t mVx=particle->Vx();
              //   Double_t mVy=particle->Vy();
                 Double_t mVz=particle->Vz();
// 25/11/2012  ???????back 10/1/2013
						TClonesArray* trArray=0;
						TParticle* tempParticle=0;

                                            TVector3 DecayMomentumK(0,0,0);  
                                              Float_t lengthKMC=0;
                      if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) {
                                                AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
                                                   lengthKMC = MCtrackReference->GetLength();
 
//			DecayMomentumK.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
                                              }
                   flenTrRef ->Fill(lengthKMC);
         flifetime ->Fill(  (lengthKMC*0.493667  /particle->P()));  // 19/7
                 if ((lengthKMC>100.)&& (lengthKMC<300.) )  flifeSmall->Fill( (lengthKMC*0.493667/particle->P() ) ); 

       Int_t nMCKpi =0;
       Int_t mcProc4 =0;
       Int_t mcProc13=0;
        Double_t radiusD=0;
   //    Double_t lengthK =0.;
       Double_t LengthK =0.;
       Double_t lenYuri =0.;
        Double_t MCQt =0.;
        Double_t MCQt3[2];
	    Int_t firstD=particle->GetFirstDaughter();
	    Int_t lastD=particle->GetLastDaughter();

             if( (lastD<=0) || (firstD<=0)  ) continue; 

                     if ( lastD > mcEvent->GetNumberOfTracks() ) continue;
                     if (firstD > mcEvent->GetNumberOfTracks() ) continue;
// 25/112012
//						TClonesArray* trArray=0;
//						TParticle* tempParticle=0;
 //                                           TVector3 DecayMomentumK(0,0,0);  
//loop on secondaries
	     for (Int_t k=firstD;k<=lastD;k++) {
              if ( k > 0 )    {
	     TParticle*    daughter1=stack->Particle(k);   // 27/8   
	     Float_t dcode = daughter1->GetPdgCode();

//     mother momentum trackrefs    and QtMC     // 17/9/2010,  test Feb 2011
						if (mcEvent->GetParticleAndTR(iMc, tempParticle, trArray) != -1) { 
     						AliTrackReference* MCtrackReference = static_cast<AliTrackReference*>(trArray->Last());
			DecayMomentumK.SetXYZ(MCtrackReference->Px(), MCtrackReference->Py(), MCtrackReference->Pz());
                              }
						const TVector3 MCP3d(daughter1->Px(), daughter1->Py(), daughter1->Pz()); //MC daughter's momentum
						 MCQt = MCP3d.Perp(DecayMomentumK); //MC daughter's transverse momentum in mother's frame
//
                  Double_t MCKinkAngle = TMath::ASin(MCQt/daughter1->P() ); 
                Double_t  MCKinkAngle2= TMath::RadToDeg() * MCKinkAngle; // in degrees 
             //    fmaxAngMomKmu->Fill(particle->P() , MCKinkAngle2);// MC 
//
	     mcProcess=daughter1->GetUniqueID();
                  radiusD=daughter1->R();
// secondary vertex
           //       Double_t hVx=daughter1->Vx();
             //   Double_t hVy=daughter1->Vy();
                  Double_t hVz=daughter1->Vz();

           LengthK = TMath::Sqrt( radiusD*radiusD  + ( mVz-hVz) * (mVz-hVz) );  //   19/7/2010 mss

//          lengthK = TMath::Sqrt( (mVx -hVx)*( mVx    -hVx)  + ( mVy    -hVy)* (mVy    -hVy ) + ( mVz   -hVz ) * (mVz     -hVz) );
            lenYuri  = (TMath::Abs( mVz-hVz))* (TMath::Sqrt( 1.+ ( ptK*ptK)/ (particle->Pz() * particle->Pz()) )) ;

       if(mcProcess==13) {
    mcProc13++;
                 
      if(mcProc13==1)     flifeInt  ->Fill(  (lengthKMC*0.493667  /particle->P()));  // 19/7

   //            if( (charg==1)&&(mcProc13==1 )) fradIntKP->Fill(daughter1->R());

     //            if( ( charg ==-1)&&(mcProc13==1))  fradIntKM->Fill(daughter1->R());
  }  


//
     if (mcProcess==4) {        
    
    mcProc4++;
   if ( mcProc4==1)  {
   // 10/1/13                if( (radiusD >120.)&&( radiusD< 210.) )  {
          flifeYuri ->Fill(  (lenYuri  *0.493667  /particle->P()));  // 19/7
                  flifetiMCK->Fill(LengthK);
                  flenYuri  ->Fill(lenYuri);
// 10/1/13                    }

  //                    flengthMCK->Fill(lengthK);  //
  //      flifetime ->Fill(  (lengthK*0.493667  /particle->P()));  // 19/7
                      flengthMCK->Fill(lengthKMC);  //
        flifetime ->Fill(  (lengthKMC*0.493667  /particle->P()));  // 19/7
            fradPtRapMC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
  }

//
    if (  ( ( code==321 )&&( dcode ==-13  ))||( ( code ==-321)&&(dcode== 13) ) || ( ( code==321 )&&( dcode ==-11  )) || ( (code ==-321)&&(dcode== 11))) {  
                      flifetim2 ->Fill(  (lengthKMC*0.493667  /particle->P()));
           //   8/2/2013allgh radius    if( (radiusD >120.)&&( radiusD< 210.) )  
           //   14/2/13 if( (radiusD >130.)&&( radiusD< 200.) )  
           if( (radiusD >fKinkRadLow )&&( radiusD< fKinkRadUp) )  
                 fmaxAngMomKmu->Fill(particle->P() , MCKinkAngle2);// MC 
              } 

      if (( (TMath::Abs(code)==321 )&&(TMath::Abs(dcode)  ==211  ))&& ( mcProc4<2)) flifetim2->Fill( lengthKMC *0.493667 /particle->P()) ;//19/7

///   inside radius region ----------------------------------------------
                       if(MCKinkAngle2 < 2.) continue;  // as in ESD 
		        //          ======  8/2/13 if (((daughter1->R())>120)&&((daughter1->R())<210)&& (MCQt>0.120)  ){
		        if (((daughter1->R())> fKinkRadLow )&&((daughter1->R())< fKinkRadUp )&& (MCQt>0.120)  ){

        if ( ( code==321 )&&( dcode ==-13  ))   {  
            fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		  fradMC->Fill(daughter1->R());
            fQtKMuMC ->Fill(MCQt );//to muon
            fgenPtEtR->Fill( ptK );//to muon
            fgenPtEtRP->Fill( ptK );//to muon
         fMCEtaKaon  ->Fill(rapidiKMC );//to muon
         fSignPtEtaMC->Fill(ptK,rapidiKMC );//to muon
         fSignPtMC->Fill(ptK);//to muon
                 flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
                    //  flifetiMCK->Fill(   LengthK*0.4933667/   ptK        
    } //  positive kaon   
        if (  ( code ==-321)&&(dcode== 13)){
        fgenPtEtR->Fill(  ptK   );//to muon
            fQtKMuMC ->Fill(MCQt );//to muon
               fgenPtEtRN->Fill(ptK);  //  
        fSignPtEtaMC->Fill(chargPt,rapidiKMC );//to muon
        fMCEtaKaon  ->Fill(rapidiKMC );//to muon
            fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		  fradMC->Fill(daughter1->R());
         fSignPtMC->Fill(chargPt);//to muon
                          flifetiMCK->Fill(   lenYuri*0.4933667/particle->P() ) ;
          //           flifetiMCK->Fill(   LengthK*0.4933667/  ptK   ) ;  

    } //  negative code
      //  if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211))) fgenPtEtR->Fill( ptK );//to pion
        if ( ( code==321 )&&( dcode ==-11  ))   {  
            fQtKElMC ->Fill(MCQt );//to muon
            fgenPtEtR->Fill( ptK );//to electron
            fgenPtEtRP->Fill( ptK );//to muon
         fMCEtaKaon  ->Fill(rapidiKMC );//to electron
            fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		  fradMC->Fill(daughter1->R());
         fSignPtEtaMC->Fill(ptK,rapidiKMC );//to electron
         fSignPtMC->Fill(ptK);
                      flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
        //             flifetiMCK->Fill(   LengthK*0.4933667/ ptK      );

    } //  positive kaon   
        if (  ( code ==-321)&&(dcode== 11)){
        fgenPtEtR->Fill(   ptK  );//to electron
            fQtKElMC ->Fill(MCQt );//to muon
               fgenPtEtRN->Fill(ptK);  //  
        fSignPtEtaMC->Fill(chargPt,rapidiKMC  );//to electron
        fMCEtaKaon  ->Fill(rapidiKMC  );//to electron
            fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
		  fradMC->Fill(daughter1->R());
         fSignPtMC->Fill(chargPt);//to electron 
                              flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
      //               flifetiMCK->Fill(   LengthK*0.4933667/   ptK         );

    } //  negative code

        if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211)))    nMCKpi++ ; 
        if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211)))   {                
                 if ( nMCKpi > 0) {
              MCQt3[nMCKpi-1] = MCQt ;// k to pipipi 
}
                       } 
    nMCKinkKs++;
       }

		    }//    decay
         } // positive k
       }//  daughters



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
              if( nMCKpi == 1) fSignPtEtaMC->Fill(chargPt,rapidiKMC );  //  k to pipi
              if( nMCKpi == 1) fSignPtMC->Fill(chargPt);  //  k to pipi
              if( nMCKpi == 1) fMCEtaKaon->Fill(rapidiKMC );  //  k to pipi
		  if(nMCKpi==1) fradMC->Fill(radiusD       );
		  if(nMCKpi==1) fQtKPiMC->Fill( MCQt   );
                             if( nMCKpi== 1) fradPtRapDC->Fill( radiusD, 1./ptK, rapidiKMC);  // systematics 26/8
                 if(nMCKpi==1)     flifetiMCK->Fill(   lenYuri*0.4933667/particle->P()  );
         //       if(nMCKpi==1)     flifetiMCK->Fill(   LengthK*0.4933667/      ptK      );

   }   //negative K

      }   /// kaons  loop 


      
    nMCTracks++;
  }// end of mc particle

//  Phys sel    2012 EFF calculation      
Bool_t isSelected =
((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB;

         if ( isSelected ==kFALSE) return;   //  24/6/11 apo MF
//
       fMultiplMC->Fill(nPrim);
//=======================================================================================
//            main vertex selection 
  const AliESDVertex *vertex=GetEventVertex(esd);
    if(!vertex) return;
///
//       fMultiplMC->Fill(nPrim);
//
/*  / apo Alexander  Feb 2012
    AliESDpid *fESDpid = new AliESDpid();
    if (!fESDpid) fESDpid =
  ((AliESDInputHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->GetESDpid();                
*/       ///========================================================================
//               apo Eftihi 
                  if(!fPIDResponse) {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler =
(AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
  }

  
    Double_t vpos[3];
     vertex->GetXYZ(vpos);
    fZpr->Fill(vpos[2]);         
    if ( TMath::Abs( vpos[2])  > 10. ) return;  ///  it is applied on ESD and generation 

  Double_t vtrack[3], ptrack[3];
  
     
 // Int_t nESDTracK = 0;

   Int_t nGoodTracks =  esd->GetNumberOfTracks();
//
    fESDMult->Fill(nGoodTracks);
     
//
      Double_t nsigmall = 100.0;
       Double_t nsigma = 100.0;
       Double_t nsigmaPion =-100.0;
       Double_t nsigmaPi=-100.0;
  //     Double_t dEdxDauMC  =   0.0;

//
for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {

    AliESDtrack* trackD = esd->GetTrack(iTrack);
    if (!trackD) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }

                Int_t indexKinkDau=trackD->GetKinkIndex(0);
// daughter kink 
          nsigmaPion     = (fPIDResponse->NumberOfSigmasTPC(trackD  , AliPID::kPion));// 26/10 eftihis
 //   nsigmaPion= (fESDpid->NumberOfSigmasTPC(trackD,AliPID::kPion));
 //   22/11/12 if((indexKinkDau >0)&& (nsigmaPion>1.2)) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
 //if((indexKinkDau >0)) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
// if((indexKinkDau >0))    dEdxDauMC   = trackD->GetTPCsignal()     ;  //  daughter kink 
   }

// loop on ESD tracks 

//
   for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {

    AliESDtrack* track = esd->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    


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

    
  Double_t tpcNCl = track->GetTPCclusters(0);  // 
  //Int_t tpcNCl = track->GetTPCclusters(0);  // 
     Double_t tpcSign = track->GetSign();  // 

  Int_t label = track->GetLabel();

    label = TMath::Abs(label);

     if(label > mcEvent->GetNumberOfTracks()) continue; //  

    TParticle * part = stack->Particle(label);
    if (!part) continue;
// loop only on Primary tracks
     if (label > nPrim) continue; /// primary tracks only   ,21/3/10  EFF study

     //    pt cut 
      if ( (track->Pt())<.200)continue;   //12/2/2012     -------------------------------------------

    UInt_t status=track->GetStatus();

     if((status&AliESDtrack::kITSrefit)==0) continue;
    if((status&AliESDtrack::kTPCrefit)==0) continue;   //6 feb
     //  if((track->GetTPCchi2()/track->GetTPCclusters(0))>3.8) continue;
     //if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.5) continue;   // test 10/3/2012
     if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.0) continue;   // test 10/3/2012

      Double_t extCovPos[15];
      track->GetExternalCovariance(extCovPos);    
//   


    track->GetXYZ(vtrack);
 fXvYv->Fill(vtrack[0],vtrack[1]);  
 fZvYv->Fill(vtrack[0],vtrack[2]);  
 fZvXv->Fill(vtrack[1],vtrack[2]);  

// track momentum
     track->GetPxPyPz(ptrack);
    
    TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
    
          Double_t   etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  );
         Double_t rapiditK = 0.5 * (TMath::Log(  (etracK + ptrack[2]  ) / ( etracK - ptrack[2])  ))  ;
    Double_t trackEta=trackMom.Eta();
   // Double_t trMoment= trackMom.Mag();
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
    
    //fRpr->Fill(dcaToVertexZpos);
    fRpr->Fill(dcaToVertexXYpos);

          //if((TMath::Abs(dcaToVertexXYpos) >0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;  //   22/7/11
   //       if((TMath::Abs(dcaToVertexXYpos) >0.24)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;  //   12/2/13
                     if (!fMaxDCAtoVtxCut->AcceptTrack(track)) continue;
    

 //  cut on eta 
 //       if(  (TMath::Abs(trackEta )) > 0.9 ) continue;
        if(  (TMath::Abs(rapiditK  )) > 0.7 ) continue; ////   rapid K, Feb 20
    fHistPtESD->Fill(track->Pt());

   // Add Kink analysis
   
   	    	Int_t indexKinkPos=track->GetKinkIndex(0);
//  loop on mother kinks
		if(indexKinkPos<0){
               fPtKink->Fill(track->Pt()); /// pt from track

                    fRatioCrossedRowsKink->Fill(ratioCrossedRowsOverFindableClustersTPC);

	// select kink class	

	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
//
	
// DCA kink
          Double_t  Dist2 = kink->GetDistance();
  //        fDCAkink->Fill( Dist2   );
//
	  Int_t eSDfLabel1=kink->GetLabel(0);
	  TParticle *particle1= stack->Particle(TMath::Abs(eSDfLabel1));
          Int_t code1= particle1->GetPdgCode();
 //         Int_t mcProcssMo= particle1->GetUniqueID();
	  
	  Int_t eSDfLabeld=kink->GetLabel(1);
	  TParticle *particled= stack->Particle(TMath::Abs(eSDfLabeld));
          Int_t dcode1= particled->GetPdgCode();
          Int_t mcProcssDa= particled->GetUniqueID();
//
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
//---------------------------edw telos  9/2010
                  Double_t hVzdau=particled->Vz();

	   const TVector3 motherMfromKink(kink->GetMotherP());
	   const TVector3 daughterMKink(kink->GetDaughterP());
	   Float_t qT=kink->GetQt();
    //   Float_t motherPt=motherMfromKink.Pt();
// Kink  mother momentum 
     Double_t trMomTPCKink=motherMfromKink.Mag();
// TPC mother momentun
     Double_t trMomTPC=track->GetTPCmomentum();
       //     Float_t etaMother=motherMfromKink.Eta();


           fHistQtAll->Fill(qT) ;  //  Qt   distr
        
	   const TVector3 vposKink(kink->GetPosition());
          fPosiKink ->Fill( vposKink[0], vposKink[1]  );

                   Double_t dxKink = vpos[0]-vposKink[0], dyKink=vpos[1]-vposKink[1], dzKink=vpos[2]-vposKink[2] ;
  //   Double_t  dzKink=vpos[2]-vposKink[2] ;    ///  ??
   Double_t lifeKink= TMath::Sqrt( dxKink*dxKink + dyKink*dyKink + dzKink*dzKink ) ;

             Double_t  tanLamda   = track-> GetTgl()    ;//  ??
                                  if (tanLamda ==0 )  continue;//   ??
        Double_t lenHelx = (TMath::Abs(dzKink    )   ) *(TMath::Sqrt( 1. + tanLamda *tanLamda  ) ) / ( TMath::Abs( tanLamda))  ;// ??




           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

//       fake kinks, small Qt and small kink angle
    if(( kinkAngle<1.))  fHistQt1  ->Fill(qT) ;  //  Qt   distr
//
    if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==321))) fFakeKPi->Fill(track->Pt());
    if(       ( (TMath::Abs(code1)==211)&&(TMath::Abs(dcode1)==211))) fFakepipi->Fill(track->Pt());
//

//remove the double taracks 
           if( (kinkAngle<2.)  ) continue;    // test 15/7 2010 , it removes 3  from   10000 good Kaons 
         //  BG  ?????==============
              if ( TMath::Abs(vposKink[2]) >  225. ) continue ; 
              if ( TMath::Abs(vposKink[2]) <  0.5 ) continue ; 

                    fPtCut1   ->Fill(trackPt );
                  fChi2NclTPC->Fill( (track->GetTPCchi2() ) ,  tpcNCl );

             Float_t signPt= tpcSign*trackPt;
//

         // ======8/1/13 if((kink->GetR()>120.)&&(kink->GetR()<210.)&&(TMath::Abs(rapiditK)<0.7)&&(label<nPrim)) {
         if((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp )&&(TMath::Abs(rapiditK)<0.7)&&(label<nPrim)) {
    if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))) fQtKMu->Fill(qT);
    if     ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))   fQtKEl->Fill(qT); 
    if     ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))   fQtKPi->Fill(qT); 
    if  (( nESDKpi>1) &&    ( ((code1)==321)&&((dcode1)==211)) )   fQtK3PiP->Fill(qT); 
    if  (( nESDKpi>1) &&    ( ((code1)==-321)&&((dcode1)==-211)) )   fQtK3PiM->Fill(qT); 
         if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) { 
         if(qT>0.120)        fHistPtKPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
        if(qT>0.120)    {
      if(code1>0.)  fHiPtKPDGP->Fill(trackPt             ); //  //  positive KAONS (pdg) inside ESD  kink sample
      if(code1<0.)  fHiPtKPDGN->Fill(      trackPt       ); //   // negative  KAONS (pdg) inside ESD  kink sample
                   }
            fHistEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
            frapidESDK->Fill(rapiditK) ;  //18/feb rapiddistr of PDG kink ESD  kaons
      if( qT > 0.120 )  fHistQt2->Fill(qT);  // PDG ESD kaons            
     }
     }


	   //Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	   Double_t maxDecAngKmu=f1->Eval(track->P(),      0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(track->P(),   0.,0.,0.);
//  two dimensional plot 
                if(TMath::Abs(code1)==321) fAngMomK->Fill(track->P(),            kinkAngle); 
                if(TMath::Abs(code1)==211) fAngMomPi->Fill( track->P(), kinkAngle); 
//_______

//
// invariant mass of mother track decaying to mu
	 Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
	 Float_t energyDaughterPi=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.139570*0.139570);
	 Float_t energyDaughterKa=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.493677*0.493677);
//	 Float_t energyDaughterPr=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.938658*0.938658);
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
         fInvMassMuNuAll ->Fill(invariantMassKmu);
           fRadiusPt->Fill( kink->GetR(), track->Pt()); //
  

               if (qT>0.120)     fSignPtNcl->Fill( signPt  ,   tpcNCl   );

//
    //  if((qT>0.12)&&((kink->GetR()>120.)&&(kink->GetR()<210.))&&(TMath::Abs(rapiditK )<0.7)) {
    if((qT>0.12)&&((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp ))&&(TMath::Abs(rapiditK )<0.7)) {
         fM1kaon->Fill(invariantMassKmu);
         fMinvPi->Fill(invariantMassKpi);
         fMinvKa->Fill(invariantMassKK);
         fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
                        }

//
         //  if ( tpcNCl<30  ) continue;
         if ( tpcNCl<20. ) continue;
                             //if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) > 0.63 ) continue;
      Double_t tpcNClHigh = -51.67+ (11./12.)  *( kink->GetR() ) ;
               if ( tpcNCl > tpcNClHigh) continue;

     Double_t tpcNClMin  = -85.5 + (65./95.)  *( kink->GetR() ) ;
               // if ( tpcNClMin < tpcNCl ) continue;   
               if ( tpcNCl < tpcNClMin ) continue;
          //  20/7/2012   if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) < 0.20 ) continue;

//              if( ( ( track->GetTPCclusters(0))  /  ( kink->GetR() ))  > 0.63 ) continue;
  //                     Int_t tpcNClMin  = -87. + (2./3.)  *( kink->GetR() ) ;
               // if ( tpcNClMin < tpcNCl ) continue;   
    //           if ( tpcNCl < tpcNClMin ) continue;

//              if( ( ( track->GetTPCclusters(0))  /  ( kink->GetR() ))  < 0.20 ) continue; //   5feb
       //  back , 20/1/2013             if (ratioCrossedRowsOverFindableClustersTPC< 0.5 )continue;// check for systematics 14/1/2013 
 //

//  kaon selection from kinks
           
   //=====  8/2/13 if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>120.)&&(kink->GetR()<210.))&&(TMath::Abs(rapiditK )<0.7)&&(invariantMassKmu<0.8)) {
   if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()> fKinkRadLow )&&(kink->GetR()< fKinkRadUp ))&&(TMath::Abs(rapiditK )<0.7)&&(invariantMassKmu<0.8)) {
// 29092010     if((kinkAngle>maxDecAngpimu)&&(qT>0.120)&&(qT<0.25)&&((kink->GetR()>120.)&&(kink->GetR()<210.))&&(TMath::Abs(rapiditK )<0.7)&&(invariantMassKmu<0.6)) {
//  if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>133.)&&(kink->GetR()<179.))&&(TMath::Abs(rapiditK )<0.5)&&(invariantMassKmu<0.6)) {   

                fAngMomKKinks->Fill(track->P(), kinkAngle);
            fPtCut2   ->Fill(trackPt );

                 if ( (kinkAngle>maxDecAngKmu*0.98)&& ( track->P() > 1.2 ) ) continue; // maximum angle s
                 if ( (kinkAngle<maxDecAngpimu*1.20)  )  continue; // maximum angle s

                fPtCut3   ->Fill(trackPt );
     //fTPCSgnlPa->Fill(track->P(),track->GetTPCsignal());
     fTPCSgnlPa->Fill(track->GetInnerParam()->GetP(),track->GetTPCsignal());
                //    if( nsigma > 3.5 )      fcode2->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
                  if(nsigma  > 3.5) continue;  // 1/11/12
               //  if(nsigma  > 4.0) continue; // test 17/2/2011  4% or more ? bg? 
//
   fTPCSgnlP->Fill(track->GetInnerParam()->GetP(), (track->GetTPCsignal()  ) ) ;
          
         fInvMassMuNuPt ->Fill(invariantMassKmu,  trackPt);
//  loop on kink daughters inside  mother's loop 
	  Int_t ESDLabelM   =  0. ;                                      
	  Int_t ESDLabelD   =  0. ;                                      
       Double_t dEdxDauMC  =   0.0;
        Double_t raDAU=0.;
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
               raDAU= kinkDau->GetR();
	   ESDLabelM=kinkDau->GetLabel(0);   //  mothers's label
                ESDLabelM = TMath::Abs(ESDLabelM);
	   ESDLabelD=kinkDau->GetLabel(1);   //  Daughter'slabel
                ESDLabelD = TMath::Abs(ESDLabelD);
                     if ( kink->GetR() == kinkDau->GetR() ) IRkink++;
                     if ( ESDLabelM == label ) Ikink++  ;
   }
  //           if (( ESDLabelM== eSDfLabel1))   { 
             if (   (Ikink >0)  && (IRkink>0 )       )   { 
// daughter kink 
     //if(indexKinkDAU >0)     nsigmaPi     = (fPIDResponse->NumberOfSigmasTPC(trackDau,AliPID::kPion));// 26/10 eftihis
     if(indexKinkDAU >0)     nsigmaPi     = (fPIDResponse->NumberOfSigmasTPC(trackDau,AliPID::kKaon));// 26/10 eftihis
 //   nsigmaPion= (fESDpid->NumberOfSigmasTPC(trackD,AliPID::kPion));
 //   22/11/12 if((indexKinkDau >0)&& (nsigmaPion>1.2)) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
 //if((indexKinkDau >0)) fTPCSgnlKinkDau->Fill(trackD->P(), (trackD->GetTPCsignal()  ) ) ;  //  daughter kink 
 if((indexKinkDAU      >0))    dEdxDauMC   = trackDau->GetTPCsignal()     ;  //  daughter kink 
   }
   }
// end internal loop for kink daughters

      //   fTPCSgnlP->Fill(trMomTPC  , (track->GetTPCsignal()  ) ) ;
      //   fTPCSgnlPtpc->Fill(trMomTPC  , (track->GetTPCsignal()  ) ) ;
         //fMothKinkMomSgnl ->Fill(trMomTPCKink  , (track->GetTPCsignal()  ) ) ;
     //    fMothKinkMomSgnl ->Fill( dEdxDauMC     , (track->GetTPCsignal()  ) ) ;
         // fTPCMomNSgnl->Fill(trMomTPC ,pidResponse->NumberOfSigmasTPC(track, AliPID::kKaon)  );     
         fNSigmTPC   ->Fill(nsigmall );
//  daughter selection 
  //fTPCSgnlKinkDau->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 
           // if(  nsigmaPion > 1.0 )        fcode4->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
       //    if( dEdxDauMC > 1.5 *(track->GetTPCsignal()   ) )       fcode4->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
// 
         //fTPCMomNSgnl->Fill(trMomTPC ,nsigmall );
         fTPCMomNSgnl->Fill(track->GetInnerParam()->GetP() ,nsigmall );
//

          fRadNclcln->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
           fRadiusPtcln->Fill( kink->GetR(), trackPt); // 

             fAngMomKC->Fill(track->P()           , kinkAngle);
                                    fHistPtKaon->Fill(track->Pt()         );   //all PID kink-kaon
        if( code1>0.)     fHistPtKaoP->Fill(track->Pt()         );   //all PID kink-kaon
        if( code1<0.)    fHistPtKaoN->Fill(track->Pt()         );   //all PID kink-kaon
// systematics
                 fradPtRapESD->Fill(kink->GetR(), 1./ track->Pt(), rapiditK);
//
                fHistEtaK->Fill(rapiditK );
 		 frad->Fill( kink->GetR() );
               flenHelx->Fill( lenHelx   );  //??
                flifeKink ->Fill(lifeKink      );//??
                        fLHelESDK ->Fill(  ( lenHelx /track->P() )*0.493677);// for all 'PID' kaons  31/7/11// ??
                      flifTiESDK->Fill(  ( lifeKink    /track->P() )*0.493677);  //  ??


 
                  fSignPtNcl->Fill( signPt  ,   tpcNCl   );
                  fSignPtEta->Fill(signPt  , rapiditK  );
                  fEtaNcl->Fill(rapiditK, tpcNCl    );
                  fSignPt->Fill( signPt   );
       fRatChi2Ncl-> Fill((track->GetTPCchi2()/track->GetTPCclusters(0) )) ;
    fdcatoVxXY->Fill(dcaToVertexXYpos);

//  kaons from k to mun and k to pipi and to e decay 
         if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) { 

         if ( label<=nPrim ) fPtPrKink->Fill(track->Pt());

                                  //         flifetim2  ->Fill(  ( lenHelx /track->P() )*0.493667); // to compare with fLHelESDK
                      flifTiESDK->Fill(  ( lifeKink    /track->P() )*0.493667);



              fKinkKaon->Fill(track->Pt());        
          fDCAkink->Fill( Dist2   );
          fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[0], vposKink[2]  );
          fPosiKinKYZ->Fill( vposKink[1], vposKink[2]  );
        if( code1>0.)     fkinkKaonP->Fill(trackPt);  //                  kPtPID kink-kaon
        if( code1<0.)    fkinkKaonN->Fill(trackPt);    //    PID kink-kaon
//     daughters
           if((((nsigmaPi) > 0.)&& ( dEdxDauMC > 70.  ) ))       fcode4->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
          if ( TMath::Abs(dEdxDauMC - track->GetTPCsignal() ) <  2)  fcode2->Fill( TMath::Abs(code1), TMath::Abs(dcode1));  
       //  fTPCSgnlPtpc->Fill(trMomTPC  , (track->GetTPCsignal()  ) ) ;
  fTPCSgnlPtpc   ->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 
         fMothKinkMomSgnlD->Fill( dEdxDauMC     , (track->GetTPCsignal()  ) ) ;
                             }
         else {
              fKinkKaonBg->Fill(track->Pt());     
         fMothKinkMomSgnl ->Fill( dEdxDauMC     , (track->GetTPCsignal()  ) ) ;
//  fTPCSgnlKinkDau->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 
   //        if( dEdxDauMC > 1.5 *(track->GetTPCsignal()   ) )       fcode4->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
      //     if(( (track->P())<1. )&& ( dEdxDauMC > 1.5 *(track->GetTPCsignal()   ) ))   fcodeDau1  ->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
           if( (( nsigmaPi) > 0. ) && ((  dEdxDauMC > 70  ) ))   fcodeDau1  ->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
          if ( TMath::Abs(dEdxDauMC - track->GetTPCsignal() ) <  2)  fcodeDau2->Fill( TMath::Abs(code1), TMath::Abs(dcode1));  
  fTPCSgnlKinkDau->Fill( daughterMKink.Mag() ,     dEdxDauMC  ) ;  //  daughter kink 
//
         fMinvPr->Fill(invariantMassKmu);
//
          fDCAkinkBG->Fill( Dist2   );
          fPosiKinKBgXY->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKBgZY->Fill( vposKink[2], vposKink[1]  );
          fPosiKinKBgZX->Fill( vposKink[2], kink->GetR() );  //  31/7/11 
        if( code1>0.)     fKinKBGP  ->Fill(   trackPt          );   //all PID kink-kaon
        if( code1<0.)     fKinKBGN  ->Fill( trackPt          );   //all PID kink-kaonl
 fdcodeH->Fill( TMath::Abs(code1), TMath::Abs(dcode1));   // put it here,  22/10/2009
         if (eSDfLabel1==eSDfLabeld)   fcodeH->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
         if (eSDfLabeld>nPrim     )   fZkinkZDau->Fill( vposKink[2],hVzdau              );

          }   // primary and all +BG    

        }  //  kink selection 
                  

	}  //End Kink Information    
  

  } //track loop 

//                                       } // vx 10 cm only on  esd
  PostData(1, fListOfHistos);

}      

//________________________________________________________________________
void AliAnalysisKinkESDMC::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}

const AliESDVertex* AliAnalysisKinkESDMC::GetEventVertex(const AliESDEvent* esd) const
  // Get the vertex from the ESD and returns it if the vertex is valid
  
{
  // Get the vertex 
  
// 10/4  const AliESDVertex* vertex = esd->GetPrimaryVertex();
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
