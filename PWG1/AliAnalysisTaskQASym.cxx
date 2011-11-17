#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "TProfile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDVertex.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMultiplicity.h"

#include "AliAnalysisTaskQASym.h"
#include "AliExternalTrackParam.h"
#include "AliTrackReference.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenDPMjetEventHeader.h"

// Analysis Task for basic QA on the ESD

// Authors: Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing,
//          Andreas Morsch, Eva Sicking

ClassImp(AliAnalysisTaskQASym)

  //________________________________________________________________________
  AliAnalysisTaskQASym::AliAnalysisTaskQASym(const char *name) 
    : AliAnalysisTaskSE(name) 
    ,fTrackType(0)
    ,fStandAlone(0)
    ,fLow(0)
    ,fHigh(100)
    ,fFieldOn(kTRUE)
    ,fHists(0)
    ,fHistRECpt(0)
    ,fEta(0)
    ,fEtaWidth(0)
    ,fPhiWidth(0)
    ,fDcaWidth(0)
    ,fPtWidth(0)
    ,fEtaPhi(0)
    ,fEtaPt(0)
    ,fQPt(0)
    ,fDca(0)
    ,fDcaZ(0)
    ,fqRec(0)
    ,fSigmaPtHist(0)
    
    ,fRecPtPos(0)
    ,fRecPtNeg(0)
    ,fRecPhiPos(0)
    ,fRecPhiNeg(0)
    ,fRecEtaPos(0)
    ,fRecEtaNeg(0)
    ,fRecEtaPtPos(0)
    ,fRecEtaPtNeg(0)
    ,fRecDcaPos(0)
    ,fRecDcaNeg(0)
    ,fRecDcaNegInv(0)
    ,fRecDPos(0)
    ,fRecDNeg(0)
    
    
    ,fRecQPtPosEta(0)
    ,fRecQPtNegEta(0)
    ,fRecPtPosEta(0)
    ,fRecPtNegEta(0)
    ,fRecPhiPosEta(0)
    ,fRecPhiNegEta(0)
    ,fRecDcaPosEta(0)
    ,fRecDcaNegEta(0)
    ,fRecDPosEta(0)
    ,fRecDNegEta(0)
    
    ,fRecPtPosVz(0)
    ,fRecPtNegVz(0)
    ,fRecEtaPosVz(0)
    ,fRecEtaNegVz(0)
    ,fRecPhiPosVz(0)
    ,fRecPhiNegVz(0)
    ,fSignedDcaPosVz(0)
    ,fSignedDcaNegVz(0)
    ,fRecQPtPosEtaVz(0)
    ,fRecQPtNegEtaVz(0)
    ,fRecEtaPtPosVz(0)
    ,fRecEtaPtNegVz(0)
    
    
    ,fDeltaPhiAll(0)
    ,fDeltaPhiLeading(0) 
    ,fDiffDcaD(0)
    
    ,fPhiRec(0)
    ,fThetaRec(0)
    ,fNumber(0)
    ,fNumberAfterCut(0)
    ,fVx(0)
    ,fVy(0)
    ,fVz(0)
    ,fNVertexSPD(0)
    ,fNVertexTracks(0)
    ,fRecDcaPosPhi(0)
    ,fRecDcaNegPhi(0)
    ,fRecPtPosPhi(0)
    ,fRecPtNegPhi(0)
    ,fRecEtaPosPhi(0)
    ,fRecEtaNegPhi(0)
    ,fRecQPtPhi(0)
    ,fRecEtaPtPosPhi(0)
    ,fRecEtaPtNegPhi(0)

    ,fRecPtPosEtaPos(0)
    ,fRecPtNegEtaPos(0)
    ,fRecPtPosEtaNeg(0)
    ,fRecPtNegEtaNeg(0)

    ,fRec1PtPosEtaPos(0)
    ,fRec1PtNegEtaPos(0)
    ,fRec1PtPosEtaNeg(0)
    ,fRec1PtNegEtaNeg(0)

    ,fRecPhiPosEtaPos(0)
    ,fRecPhiNegEtaPos(0)
    ,fRecPhiPosEtaNeg(0)
    ,fRecPhiNegEtaNeg(0)

    ,fRecDcaPosPhiEtaPos(0)
    ,fRecDcaNegPhiEtaPos(0) 
    ,fRecDcaPosPhiEtaNeg(0)  
    ,fRecDcaNegPhiEtaNeg(0)  

    ,fRecDcaPosPtEtaPos(0)
    ,fRecDcaNegPtEtaPos(0) 
    ,fRecDcaPosPtEtaNeg(0)  
    ,fRecDcaNegPtEtaNeg(0)  
  
    ,fRecPtPosPhiEtaPos(0)  
    ,fRecPtNegPhiEtaPos(0)  
    ,fRecPtPosPhiEtaNeg(0) 
    ,fRecPtNegPhiEtaNeg(0) 


//    ,fRecDcaPhiPtPosEtaPos(0)
//    ,fRecDcaPhiPtNegEtaPos(0)
//    ,fRecDcaPhiPtPosEtaNeg(0)  
//    ,fRecDcaPhiPtNegEtaNeg(0)  

    ,fEtavPt(0)  
    ,fPhivPt(0) 
 
    ,fCompareTPCparam(0)

    ,fITSlayer(0)
    ,fITSlayerEta(0)
    ,fITSlayerPhi(0)

    ,fCuts(0)

{
  // Constructor
  for(Int_t i = 0;i<4;++i){
    fVertexX[i]=0;
    fVertexY[i]=0;
    fVertexZ[i]=0;
  }

  for(Int_t i = 0;i<18;++i){
    fRecPtTpcSector[i] = 0;
    fRecEtaTpcSector[i] = 0;
    fSignedDcaTpcSector[i] = 0;
    fRecQPtTpcSector[i] = 0;
    fRecEtaPtTpcSector[i] = 0;
  }

  for(Int_t i = 0;i< 7;++i){
    fRecPtPosLadder[i] = 0;
    fRecPtNegLadder[i] = 0;
    fRecPhiPosLadder[i] = 0;
    fRecPhiNegLadder[i] = 0;
    fRecEtaPosLadder[i] = 0;
    fRecEtaNegLadder[i] = 0;
    fSignDcaPos[i] = 0;
    fSignDcaNeg[i] = 0;
    fSignDcaNegInv[i] = 0;
    fPtSigmaPos[i] =0;
    fPtSigmaNeg[i] =0;
    fqPtRec[i] =0;
    fDcaSigmaPos[i] =0;
    fDcaSigmaNeg[i] =0;
  }

  for(Int_t i = 0;i< 3;i++){
    for(Int_t j = 0;j< 2;j++){
      fEtaBinPt[i][j]=0;
      fPhiBinPt[i][j]=0;
      fDcaBinPt[i][j]=0;
      fEtaPhiBinPt[i][j]=0;
    }
  }

  DefineOutput(1,  TList::Class()); 

  
  
}


//________________________________________________________________________
void AliAnalysisTaskQASym::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  Double_t range = 0.3;
  Double_t pt = 20.;

  fHists = new TList();

  fHistRECpt   = new TH1F("fHistRECpt", 
			  " p_{T}",
			  200, 0., pt);
  fEta   = new TH1F("fEta", 
		    " #eta",
		    200, -2., 2.);
  fEtaWidth   = new TH1F("fEtaWidth", 
			 " #eta",
			 200, -2., 2.);
  fPhiWidth   = new TH1F("fPhiWidth", 
			 " #phi",
			 200, 0., 2*TMath::Pi());
  fDcaWidth   = new TH1F("fDcaWidth", 
			 "dca",
			 200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fPtWidth   = new TH1F("fPtWidth", 
			 "p_{T}",
			 200, 0., pt);
  fEtavPt   = new TH2F("fEtavPt", 
		       " #eta -p_{T}",
		       200, -2., 2.,
		       100, -3, 4);
  fPhivPt   = new TH2F("fPhivPt", 
		       " #phi -p_{T}",
		       200, 0, 2*TMath::Pi(),
		       100, -3, 5);
  fCompareTPCparam   = new TH2F("fCompareTPCparam", 
				"fCompareTPCparam",
				100, -1., 1.,100,-5, 5);

  fITSlayer   = new TH1F("fITSlayer", 
			 "fITSlayer",
			 8, -1.5, 6.5);
  fITSlayerEta   = new TH2F("fITSlayerEta", 
			 "fITSlayerEta",
			    8, -1.5, 6.5, 200, -2.,2.);
  fITSlayerPhi   = new TH2F("fITSlayerPhi", 
			    "fITSlayerPhi",
			    8, -1.5, 6.5, 200, 0,2*TMath::Pi());
  fEtaPhi   = new TH2F("fEtaPhi", 
		       " #eta - #phi",
		       200, -2., 2., 128, 0., 2. * TMath::Pi());
  fThetaRec   = new TH1F("fThetaRec", 
			 " #theta",
			 180, 0., TMath::Pi());
  fPhiRec   = new TH1F("fPhiRec", 
		       " #phi",
		       180, 0., 2*TMath::Pi());
  fNumber   = new TH1F("fNumber", 
		       "number of tracks per event",
		       500, -5, 4995);
  fNumberAfterCut   = new TH1F("fNumberAfterCut", 
			       "number of tracks per event after cuts",
			       500, -5, 4995);
  fVx   = new TH1F("fVx", 
		   "X of first track point",
		   100, -1., 1.);
  fVy   = new TH1F("fVy", 
		   "Y of first track point",
		   100, -1., 1.);
  fVz   = new TH1F("fVz", 
		   "Z of first track point",
		   200, -50., 50.);
  fNVertexSPD   = new TH1F("fNVertexSPD", 
			"Number of SPD vertices",
			10, -0.5, 9.5);
  fNVertexTracks   = new TH1F("fNVertexTracks", 
			      "Number of track vertices",
			      10, -0.5, 9.5);
  
  fEtaPt   = new TH1F("fEtaPt", 
		      " #eta/p_{T} ",
		      100, -1., 1.);

  fQPt   = new TH1F("fQPt", 
		    " charge/p_{T} ",
		    100, -1., 1.);

  fDca   = new TH1F("fDca", 
		    " dca ",
		    200,  -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));

  fDcaZ   = new TH1F("fDcaZ", "fDcaZ ",200,  -3, 3);// limits fitting those of track cuts


  fqRec    = new TH1F("fqRec",   
		      " charge all reconstructed particle",
		      21, -9.5, 10.5);
  
  fSigmaPtHist    = new TH1F("fSigmaPtHist",   
			 "Log_{10}(#sigma_{p_{T}})",
			 200, -4., 8.);


  TString lable[4]={"", "SPD", "Track", "TPC"};
  for(Int_t i=0;i<4;i++){
    fVertexX[i]   = new TH1F(Form("fVertexX%s",lable[i].Data()),
			     Form("fVertexX%s",lable[i].Data()),
			     100, -1., 1.);
    fVertexY[i]   = new TH1F(Form("fVertexY%s",lable[i].Data()),
			     Form("fVertexY%s",lable[i].Data()),
			     100, -1., 1.);
    if(i==1 || i==2){
      fVertexZ[i]   = new TH1F(Form("fVertexZ%s",lable[i].Data()),
			       Form("fVertexZ%s",lable[i].Data()),
			       200, -5., 5.);
    }
    else{
      fVertexZ[i]   = new TH1F(Form("fVertexZ%s",lable[i].Data()),
			       Form("fVertexZ%s",lable[i].Data()),
			       200, -50., 50.);
    }
  }

  //------------
  for(Int_t ITSlayer_case=0;ITSlayer_case<7;ITSlayer_case++){

    fSignDcaPos[ITSlayer_case]   = new TH1F(Form("fSignDcaPos%d", ITSlayer_case),  
					    " Signed dca", 
					    200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
    fSignDcaPos[ITSlayer_case]->GetXaxis()->SetTitle("dca");
    fSignDcaPos[ITSlayer_case]->GetYaxis()->SetTitle("");
   
 
    fSignDcaNeg[ITSlayer_case]   = new TH1F(Form("fSignDcaNeg%d", ITSlayer_case),  
					    " Signed dcas",
					    200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
    fSignDcaNeg[ITSlayer_case]->GetXaxis()->SetTitle("dca");
    fSignDcaNeg[ITSlayer_case]->GetYaxis()->SetTitle("");

    fSignDcaNegInv[ITSlayer_case]   = new TH1F(Form("fSignDcaNegInv%d", ITSlayer_case),  
					       " inverse Signed dca ",
					       200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
    fSignDcaNegInv[ITSlayer_case]->GetXaxis()->SetTitle("-dca");
    fSignDcaNegInv[ITSlayer_case]->GetYaxis()->SetTitle("");




    fPtSigmaPos[ITSlayer_case]   = new TH1F(Form("fPtSigmaPos%d", ITSlayer_case),  
					    " #sigma_{pT} ",
					    208, -4., 8.);
    fPtSigmaPos[ITSlayer_case]->GetXaxis()->SetTitle("Log_{10}(#sigma_{pT})");
    fPtSigmaPos[ITSlayer_case]->GetYaxis()->SetTitle("");
    
    
    fPtSigmaNeg[ITSlayer_case]   = new TH1F(Form("fPtSigmaNeg%d",ITSlayer_case),  
					    " #sigma_{pT}",
					    208, -4., 8.);
    fPtSigmaNeg[ITSlayer_case]->GetXaxis()->SetTitle("Log_{10}(#sigma_{pT})");
    fPtSigmaNeg[ITSlayer_case]->GetYaxis()->SetTitle("");





    fqPtRec[ITSlayer_case]   = new TH1F(Form("fqPtRec%d",ITSlayer_case),  
					"q/ p_{T}",
					200, -100., 100.);
    fqPtRec[ITSlayer_case]->GetXaxis()->SetTitle("q_{tr}/p_{T, tr} (GeV/c)");
    fqPtRec[ITSlayer_case]->GetYaxis()->SetTitle("");

  
   


    fDcaSigmaPos[ITSlayer_case]   = new TH2F(Form("fDcaSigmaPos%d", ITSlayer_case),  
					     " p_{T} shift vs #sigma_{pT} ",
					     200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9),200, -4., 4. );
    fDcaSigmaPos[ITSlayer_case]->GetXaxis()->SetTitle("signed DCA)");
    fDcaSigmaPos[ITSlayer_case]->GetYaxis()->SetTitle("log_{10}(#sigma_{pT})");
    
    
    fDcaSigmaNeg[ITSlayer_case]   = new TH2F(Form("fDcaSigmaNeg%d", ITSlayer_case),  
					     " p_{T} shift vs #sigma_{pT} ",
					     200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9),200, -4., 4. );
    fDcaSigmaNeg[ITSlayer_case]->GetXaxis()->SetTitle("signed DCA");
    fDcaSigmaNeg[ITSlayer_case]->GetYaxis()->SetTitle("log_{10}(#sigma_{pT})");


  }
    
 
    
  // YIELDs---------- positive and negative particles
    
  fRecPtPos   = new TH1F("fRecPtPos", 
			 " p_{T}",
			 100, 0.,pt);
  fRecPtPos->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fRecPtNeg   = new TH1F("fRecPtNeg", 
			 " p_{T} ",
			 100, 0., pt);
  fRecPtNeg->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    
  fRecPhiPos   = new TH1F("fRecPhiPos", 
			  "#phi",
			  361, 0., 360.);
  fRecPhiPos->GetXaxis()->SetTitle("#phi (deg)");
  
  fRecPhiNeg   = new TH1F("fRecPhiNeg", 
			  "#phi ",
			  361, 0., 360.);
  fRecPhiNeg->GetXaxis()->SetTitle("#phi (deg)");
    
  fRecEtaPos   = new TH1F("fRecEtaPos", 
			  "#eta",
			  200, -2., 2.);
  fRecEtaPos->GetXaxis()->SetTitle("#eta");

  fRecEtaNeg   = new TH1F("fRecEtaNeg", 
			  "#eta",
			  200, -2., 2.);
  fRecEtaNeg->GetXaxis()->SetTitle("#eta");
    
  fRecEtaPtPos   = new TH1F("fRecEtaPtPos", 
			    "#eta/p_{T}",
			    200, -0.1, .1);
  fRecEtaPtPos->GetXaxis()->SetTitle("#eta/p_{T}");

  fRecEtaPtNeg   = new TH1F("fRecEtaPtNeg", 
			    "#eta/p_{T}",
			    200, -.1, .1);
  fRecEtaPtNeg->GetXaxis()->SetTitle("#eta/p_{T}");

  fRecDcaPos   = new TH1F("fRecDcaPos", 
			  " dca",
			  100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDcaPos->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNeg   = new TH1F("fRecDcaNeg", 
			  " dca",
			  100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDcaNeg->GetXaxis()->SetTitle("dca (cm)");

  fRecDcaNegInv   = new TH1F("fRecDcaNegInv", 
			     " dca",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDcaNegInv->GetXaxis()->SetTitle("dca (cm)");


  fRecDPos   = new TH1F("fRecDPos", 
			" d",
			100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDPos->GetXaxis()->SetTitle("d (cm)");
  fRecDNeg   = new TH1F("fRecDNeg", 
			"d",
			100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDNeg->GetXaxis()->SetTitle("d (cm)");


  //  YIELDs ---------------- positive and negative eta
    
    
  fRecQPtPosEta   = new TH1F("fRecQPtPosEta", 
			     "q/p_{T}",
			     200, -0.5, 0.5);
  fRecQPtPosEta->GetXaxis()->SetTitle("q/p_{T} ");

  fRecQPtNegEta   = new TH1F("fRecQPtNegEta", 
			     "q/p_{T}",
			     200, -0.5, 0.5);
  fRecQPtNegEta->GetXaxis()->SetTitle("q/p_{T}");
    
  fRecPtPosEta   = new TH1F("fRecPtPosEta", 
			    " p_{T} ",
			    100, 0., pt);
  fRecPtPosEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fRecPtNegEta   = new TH1F("fRecPtNegEta", 
			    " p_{T}",
			    100, 0., pt);
  fRecPtNegEta->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    
  fRecPhiPosEta   = new TH1F("fRecPhiPosEta", 
			     "#phi",
			     361, 0., 360);
  fRecPhiPosEta->GetXaxis()->SetTitle("#phi (deg)");

  fRecPhiNegEta   = new TH1F("fRecPhiNegEta", 
			     "#phi ",
			     361, 0, 360);
  fRecPhiNegEta->GetXaxis()->SetTitle("#phi (deg)");

  fRecDcaPosEta   = new TH1F("fRecDcaPosEta", 
			     " dca ",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDcaPosEta->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegEta   = new TH1F("fRecDcaNegEta", 
			     " dca",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDcaNegEta->GetXaxis()->SetTitle("dca (cm)");

  fRecDPosEta   = new TH1F("fRecDPosEta", 
			   " d",
			   100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
  fRecDPosEta->GetXaxis()->SetTitle("d (cm)");
  fRecDNegEta   = new TH1F("fRecDNegEta", 
			   "d",
			   100, -5., 5.);
  fRecDNegEta->GetXaxis()->SetTitle("d (cm)");

  fRecDcaPosPhi   = new TH2F("fRecDcaPosPhi", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaPosPhi->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaPosPhi->GetYaxis()->SetTitle("#phi (rad.)");
  fRecDcaNegPhi   = new TH2F("fRecDcaNegPhi", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaNegPhi->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegPhi->GetYaxis()->SetTitle("#phi (rad.)");

  fRecPtPosPhi   = new TH2F("fRecPtPosPhi", 
			     " log(p_T) vs. phi",
			     100, -2.5, 2., 180, 0, TMath::Pi()*2);
  fRecPtPosPhi->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtPosPhi->GetYaxis()->SetTitle("#phi (rad.)");
  fRecPtNegPhi   = new TH2F("fRecPtNegPhi", 
			     " log(p_T) vs. phi",
			     100,-2.5 , 2., 180, 0, TMath::Pi()*2);
  fRecPtNegPhi->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtNegPhi->GetYaxis()->SetTitle("#phi (rad.)");

  fRecEtaPosPhi   = new TH2F("fRecEtaPosPhi", 
			     " eta vs. phi",
			     100, -1.5, 1.5, 180, 0, TMath::Pi()*2);
  fRecEtaPosPhi->GetXaxis()->SetTitle("#eta");
  fRecEtaPosPhi->GetYaxis()->SetTitle("#phi (rad.)");
  fRecEtaNegPhi   = new TH2F("fRecEtaNegPhi", 
			     " eta vs. phi",
			     100, -1.5, 1.5, 180, 0, TMath::Pi()*2);
  fRecEtaNegPhi->GetXaxis()->SetTitle("#eta");
  fRecEtaNegPhi->GetYaxis()->SetTitle("#phi (rad.)");

  fRecQPtPhi   = new TH2F("fRecQPtPhi", 
			     " charge/p_T vs. phi",
			     100,-1. , 1., 180, 0, TMath::Pi()*2);
  fRecQPtPhi->GetXaxis()->SetTitle("charge/p_{T}");
  fRecQPtPhi->GetYaxis()->SetTitle("#phi (rad.)");

  fRecEtaPtPosPhi   = new TH2F("fRecEtaPtPosPhi", 
			     " eta/p_T vs. phi",
			     100, -5, 5., 180, 0, TMath::Pi()*2);
  fRecEtaPtPosPhi->GetXaxis()->SetTitle("#eta/p_{T}");
  fRecEtaPtPosPhi->GetYaxis()->SetTitle("#phi (rad.)");
  fRecEtaPtNegPhi   = new TH2F("fRecEtaPtNegPhi", 
			     " eta/p_T vs. phi",
			     100,-5 , 5., 180, 0, TMath::Pi()*2);
  fRecEtaPtNegPhi->GetXaxis()->SetTitle("#eta/p_{T}");
  fRecEtaPtNegPhi->GetYaxis()->SetTitle("#phi (rad.)");





  fRecDcaPosPhiEtaPos   = new TH2F("fRecDcaPosPhiEtaPos", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaPosPhiEtaPos->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaPosPhiEtaPos->GetYaxis()->SetTitle("#phi (rad.)");
  fRecDcaNegPhiEtaPos   = new TH2F("fRecDcaNegPhiEtaPos", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaNegPhiEtaPos->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegPhiEtaPos->GetYaxis()->SetTitle("#phi (rad.)");

  fRecPtPosPhiEtaPos   = new TH2F("fRecPtPosPhiEtaPos", 
			     " log(p_T) vs. phi",
			     100, -2.5, 2., 180, 0, TMath::Pi()*2);
  fRecPtPosPhiEtaPos->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtPosPhiEtaPos->GetYaxis()->SetTitle("#phi (rad.)");
  fRecPtNegPhiEtaPos   = new TH2F("fRecPtNegPhiEtaPos", 
			     " log(p_T) vs. phi",
			     100,-2.5 , 2., 180, 0, TMath::Pi()*2);
  fRecPtNegPhiEtaPos->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtNegPhiEtaPos->GetYaxis()->SetTitle("#phi (rad.)");


  fRecDcaPosPhiEtaNeg   = new TH2F("fRecDcaPosPhiEtaNeg", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaPosPhiEtaNeg->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaPosPhiEtaNeg->GetYaxis()->SetTitle("#phi (rad.)");
  fRecDcaNegPhiEtaNeg   = new TH2F("fRecDcaNegPhiEtaNeg", 
			     " dca vs. phi",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 180, 0, TMath::Pi()*2);
  fRecDcaNegPhiEtaNeg->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegPhiEtaNeg->GetYaxis()->SetTitle("#phi (rad.)");

  fRecPtPosPhiEtaNeg   = new TH2F("fRecPtPosPhiEtaNeg", 
			     " log(p_T) vs. phi",
			     100, -2.5, 2., 180, 0, TMath::Pi()*2);
  fRecPtPosPhiEtaNeg->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtPosPhiEtaNeg->GetYaxis()->SetTitle("#phi (rad.)");
  fRecPtNegPhiEtaNeg   = new TH2F("fRecPtNegPhiEtaNeg", 
			     " log(p_T) vs. phi",
			     100,-2.5 , 2., 180, 0, TMath::Pi()*2);
  fRecPtNegPhiEtaNeg->GetXaxis()->SetTitle("log_{10}(p_{T})");
  fRecPtNegPhiEtaNeg->GetYaxis()->SetTitle("#phi (rad.)");
  
  //new
  fRecDcaPosPtEtaPos   = new TH2F("fRecDcaPosPtEtaPos", 
			     " dca vs. pt",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 200, -2, 2);
  fRecDcaPosPtEtaPos->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaPosPtEtaPos->GetYaxis()->SetTitle("log_{10}(p_{T})");

  fRecDcaPosPtEtaNeg   = new TH2F("fRecDcaPosPtEtaNeg", 
			     " dca vs. pt",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 200, -2, 2);
  fRecDcaPosPtEtaNeg->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaPosPtEtaNeg->GetYaxis()->SetTitle("log_{10}(p_{T})");

  fRecDcaNegPtEtaPos   = new TH2F("fRecDcaNegPtEtaPos", 
			     " dca vs. pt",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 200, -2, 2);
  fRecDcaNegPtEtaPos->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegPtEtaPos->GetYaxis()->SetTitle("log_{10}(p_{T})");

  fRecDcaNegPtEtaNeg   = new TH2F("fRecDcaNegPtEtaNeg", 
			     " dca vs. pt",
			     100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 200, -2, 2);
  fRecDcaNegPtEtaNeg->GetXaxis()->SetTitle("dca (cm)");
  fRecDcaNegPtEtaNeg->GetYaxis()->SetTitle("log_{10}(p_{T})");
    


  //  YIELDs ---------------- for TPC sectors
  for(Int_t sector=0; sector<18;sector++){
      

    fRecPtTpcSector[sector]   = new TH1F(Form("fRecPtTpcSector%02d",sector), 
					 Form("p_{T} distribution: TPC sector %d",
					      sector),100, 0., pt);
    fRecPtTpcSector[sector]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    fRecEtaTpcSector[sector]   = new TH1F(Form("fRecEtaTpcSector%02d",sector), 
					  Form("#eta distribution: TPC sector %d",
					       sector),200, -2., 2.);
    fRecEtaTpcSector[sector]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
     

    fSignedDcaTpcSector[sector]   = new TH1F(Form("fSignedDcaTpcSector%02d",sector), 
					     Form("dca distribution: TPC sector %d",
						  sector),200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9) );
    fSignedDcaTpcSector[sector]->GetXaxis()->SetTitle("dca");

    fRecQPtTpcSector[sector]   = new TH1F(Form("fRecQPtTpcSector%02d",sector), 
					  Form("Q/ p_{T} distribution: TPC sector %d",
					       sector),100, -1., 1.);
    fRecQPtTpcSector[sector]->GetXaxis()->SetTitle("Q/p_{T} (GeV/c)");

    fRecEtaPtTpcSector[sector]   = new TH1F(Form("fRecEtaPtTpcSector%02d",sector), 
					    Form("#eta/ p_{T} distribution: TPC sector %d",
						 sector),100, -1., 1.);
    fRecEtaPtTpcSector[sector]->GetXaxis()->SetTitle("#eta/p_{T} (GeV/c)");
 
  }
  // YIELDS ITS ladder
  for(Int_t i=0;i<7;i++){
    fRecPtPosLadder[i]   = new TH1F(Form("fRecPtPosLadder%d", i), 
				    " p_{T} distribution",
				    100, 0., pt);
    fRecPtPosLadder[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fRecPtNegLadder[i]   = new TH1F(Form("fRecPtNegLadder%d",i), 
				    " p_{T} distribution ",
				    100, 0., pt);
    fRecPtNegLadder[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");


    fRecPhiPosLadder[i]   = new TH1F(Form("fRecPhiPosLadder%d",i), 
				     "#phi distribution: all pos eta",
				     361, 0., 360);
    fRecPhiPosLadder[i]->GetXaxis()->SetTitle("#phi (deg)");
      
    fRecPhiNegLadder[i]   = new TH1F(Form("fRecPhiNegLadder%d", i),
				     "#phi distribution: all neg eta",
				     361, 0, 360);
    fRecPhiNegLadder[i]->GetXaxis()->SetTitle("#phi (deg)");



    fRecEtaPosLadder[i]   = new TH1F(Form("fRecEtaPosLadder%d",i), 
				     "#eta distribution",
				     200, -2., 2.);
    fRecEtaPosLadder[i]->GetXaxis()->SetTitle("#eta)");
      
    fRecEtaNegLadder[i]   = new TH1F(Form("fRecEtaNegLadder%d", i),
				     "#eta distribution",
				     200, -2., 2.);
    fRecEtaNegLadder[i]->GetXaxis()->SetTitle("#eta");
  }

  Double_t vzmax = 15.;

  fRecPtPosVz = new TH2F("fRecPtPosVz", 
			 "p_{T} distribution vs Vz()",
			 100, -1., 2., 200,-vzmax,vzmax);
  fRecPtPosVz->GetXaxis()->SetTitle("log_{10}(p_{T})");
    
  fRecPtNegVz = new TH2F("fRecPtNegVz",
			 "p_{T} distribution vs Vz()",
			 100, -1., 2.,200,-vzmax,vzmax);
  fRecPtNegVz->GetXaxis()->SetTitle("Log_{10}(p_{T})");
    
   
  fRecEtaPosVz= new TH2F("fRecEtaPosVz", 
			 "#eta distribution vs Vz()",
			 100, -2., 2., 200,-vzmax,vzmax);
  fRecEtaPosVz->GetXaxis()->SetTitle("#eta");
  fRecEtaNegVz = new TH2F("fRecEtaNegVz",
			  "#eta distribution vs Vz()",
			  100, -2., 2.,200,-vzmax,vzmax);
  fRecEtaNegVz->GetXaxis()->SetTitle("#eta");

  fRecPhiPosVz= new TH2F("fRecPhiPosVz", 
			 "#eta distribution vs Vz()",
			 361, 0., 360., 200,-vzmax,vzmax);
  fRecPhiPosVz->GetXaxis()->SetTitle("#phi (deg)");
  fRecPhiNegVz = new TH2F("fRecPhiNegVz",
			  "dca vs Vz()",
			  361, 0., 360.,200,-vzmax,vzmax);
  fRecPhiNegVz->GetXaxis()->SetTitle("#phi (deg)");

  fSignedDcaPosVz= new TH2F("fSignedDcaPosVz", 
			    "#eta distribution vs Vz()",
			    200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9), 200,-vzmax,vzmax);
  fSignedDcaPosVz->GetXaxis()->SetTitle("dca (cm)");
  fSignedDcaNegVz = new TH2F("fSignedDcaNegVz",
			     "dca vs Vz()",
			     200, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9),200,-vzmax,vzmax);
  fSignedDcaNegVz->GetXaxis()->SetTitle("dca (cm)");

  fRecQPtPosEtaVz= new TH2F("fRecQPtPosEtaVz",
			    " Q/p_{T} distribution vs Vz()",
			    100, -1., 1., 200,-vzmax,vzmax);
  fRecQPtPosEtaVz->GetXaxis()->SetTitle("Q/p_{T}");
  fRecQPtNegEtaVz = new TH2F("fRecQPtNegEtaVz",
			     " Q/p_{T} distribution vs Vz()",
			     100, -1., 1.,200,-vzmax,vzmax);
  fRecQPtNegEtaVz->GetXaxis()->SetTitle("Q/p_{T}");

 
  fRecEtaPtPosVz= new TH2F("fRecEtaPtPosVz",
			   " #eta/p_{T} distribution vs Vz()",
			   100, -1., 1., 200,-vzmax,vzmax);
  fRecEtaPtPosVz->GetXaxis()->SetTitle("#eta/p_{T");
  fRecEtaPtNegVz = new TH2F("fRecEtaPtNegVz",
			    " #eta/p_{T} distribution vs Vz()",
			    100, -1., 1.,200,-vzmax,vzmax);
  fRecEtaPtNegVz->GetXaxis()->SetTitle("#eta/p_{T}");

  //new
  fDeltaPhiAll = new TH1F("fDeltaPhiAll",
			  " #Delta #phi",200,-360,360);
  fDeltaPhiAll->GetXaxis()->SetTitle("#Delta #phi");


  fDeltaPhiLeading = new TH2F("fDeltaPhiLeading",
			      " #Delta #phi",361,-360,360, 361,0, 360);
  fDeltaPhiLeading->GetXaxis()->SetTitle("#Delta #phi (deg.)");
  fDeltaPhiLeading->GetYaxis()->SetTitle("#phi_{leading particle} (deg.)");

  fDiffDcaD    = new TH1F("fDiffDcaD",   
			  "dca-d",
			  200, -1., 1.);

  
  fRecPtPosEtaPos = new TH1F("fRecPtPosEtaPos",
			     "p_{T} distribution",100,0,pt);
  fRecPtPosEtaPos->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fRecPtNegEtaPos = new TH1F("fRecPtNegEtaPos",
			     "p_{T} distribution",100,0,pt);
  fRecPtNegEtaPos->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fRecPtPosEtaNeg = new TH1F("fRecPtPosEtaNeg",
			     "p_{T} distribution",100,0,pt);
  fRecPtPosEtaNeg->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  fRecPtNegEtaNeg = new TH1F("fRecPtNegEtaNeg",
			     "p_{T} distribution",100,0,pt);
  fRecPtNegEtaNeg->GetXaxis()->SetTitle("p_{T} (GeV/c)");



  fRec1PtPosEtaPos = new TH1F("fRec1PtPosEtaPos",
			     "1/p_{T} distribution",100,0,0.5);
  fRec1PtPosEtaPos->GetXaxis()->SetTitle("p_{T} (c/GeV)");

  fRec1PtNegEtaPos = new TH1F("fRec1PtNegEtaPos",
			     "1/p_{T} distribution",100,0,0.5);
  fRec1PtNegEtaPos->GetXaxis()->SetTitle("p_{T} (c/GeV)");

  fRec1PtPosEtaNeg = new TH1F("fRec1PtPosEtaNeg",
			     "1/p_{T} distribution",100,0,0.5);
  fRec1PtPosEtaNeg->GetXaxis()->SetTitle("p_{T} (c/GeV)");

  fRec1PtNegEtaNeg = new TH1F("fRec1PtNegEtaNeg",
			     "1/p_{T} distribution",100,0,0.5);
  fRec1PtNegEtaNeg->GetXaxis()->SetTitle("1/p_{T} (c/GeV)");


 
  fRecPhiPosEtaPos = new TH1F("fRecPhiPosEtaPos",
			     "#phi",180,0,2*TMath::Pi());
  fRecPhiPosEtaPos->GetXaxis()->SetTitle("#phi (rad.)");

  fRecPhiNegEtaPos = new TH1F("fRecPhiNegEtaPos",
			     "#phi",180,0,2*TMath::Pi());
  fRecPhiNegEtaPos->GetXaxis()->SetTitle("#phi (rad.)");

  fRecPhiPosEtaNeg = new TH1F("fRecPhiPosEtaNeg",
			     "#phi",180,0,2*TMath::Pi());
  fRecPhiPosEtaNeg->GetXaxis()->SetTitle("#phi (rad.)");

  fRecPhiNegEtaNeg = new TH1F("fRecPhiNegEtaNeg",
			     "#phi",180,0,2*TMath::Pi());
  fRecPhiNegEtaNeg->GetXaxis()->SetTitle("#phi (rad.)");



//   fRecDcaPhiPtPosEtaPos = new TH3F("fRecDcaPhiPtPosEtaPos",
// 				   "#phi- p_{T} - DCA",
// 				   180,0,2*TMath::Pi(),
// 				   100,0,pt,
// 				   100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
//   fRecDcaPhiPtPosEtaPos->GetXaxis()->SetTitle("#phi (rad.)");
//   fRecDcaPhiPtPosEtaPos->GetYaxis()->SetTitle("p_{T} (GeV/c)");
//   fRecDcaPhiPtPosEtaPos->GetZaxis()->SetTitle("dca (cm)");

//   fRecDcaPhiPtPosEtaNeg = new TH3F("fRecDcaPhiPtPosEtaNeg",
// 				   "#phi- p_{T} - DCA",
// 				   180,0,2*TMath::Pi(),
// 				   100,0,pt,
// 				   100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
//   fRecDcaPhiPtPosEtaNeg->GetZaxis()->SetTitle("dca (cm)");
//   fRecDcaPhiPtPosEtaNeg->GetXaxis()->SetTitle("#phi (rad.)");
//   fRecDcaPhiPtPosEtaNeg->GetYaxis()->SetTitle("p_{T} (GeV/c)");

//   fRecDcaPhiPtNegEtaPos = new TH3F("fRecDcaPhiPtNegEtaPos",
// 				   "#phi- p_{T} - DCA",
// 				   180,0,2*TMath::Pi(),
// 				   100,0,pt,
// 				   100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
//   fRecDcaPhiPtNegEtaPos->GetZaxis()->SetTitle("dca (cm)");
//   fRecDcaPhiPtNegEtaPos->GetXaxis()->SetTitle("#phi (rad.)");
//   fRecDcaPhiPtNegEtaPos->GetYaxis()->SetTitle("p_{T} (GeV/c)");

//   fRecDcaPhiPtNegEtaNeg = new TH3F("fRecDcaPhiPtNegEtaNeg",
// 				   "#phi- p_{T} - DCA",
// 				   180,0,2*TMath::Pi(),
// 				   100,0,pt,
// 				   100, -range*(1+Int_t(fTrackType/2)*9), range*(1+Int_t(fTrackType/2)*9));
//   fRecDcaPhiPtNegEtaNeg->GetZaxis()->SetTitle("dca (cm)");
//   fRecDcaPhiPtNegEtaNeg->GetYaxis()->SetTitle("#phi (rad.)");
//   fRecDcaPhiPtNegEtaNeg->GetXaxis()->SetTitle("p_{T} (GeV/c)");

  TString charge[2];
  charge[0]="Pos";
  charge[1]="Neg";

  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<2;j++){
      fEtaBinPt[i][j]   = new TH1F(Form("fEtaBinPt%d%s", i, charge[j].Data()), 
				   "eta",
				   200, -2., 2.);
      fPhiBinPt[i][j]   = new TH1F(Form("fPhiBinPt%d%s", i,charge[j].Data() ), 
				   "phi",
				   181, 0, 2*TMath::Pi());
      fDcaBinPt[i][j]   = new TH1F(Form("fDcaBinPt%d%s", i, charge[j].Data()), 
				   "DCA",
				   200,-range*(1+Int_t(fTrackType/2)*9),
				   range*(1+Int_t(fTrackType/2)*9) );
      fEtaPhiBinPt[i][j]= new TH2F(Form("fEtaPhiBinPt%d%s", i, charge[j].Data()), 
				      "eta-phi",
				      200, -2., 2., 200, 0.,2*TMath::Pi());
    }
  }


  fHists->SetOwner();

  fHists->Add(fHistRECpt);
  fHists->Add(fEta);
  fHists->Add(fEtaWidth);
  fHists->Add(fPhiWidth);
  fHists->Add(fDcaWidth);
  fHists->Add(fPtWidth);
  fHists->Add(fEtavPt);
  fHists->Add(fPhivPt);
  fHists->Add(fCompareTPCparam);
  fHists->Add(fITSlayer);
  fHists->Add(fITSlayerEta);
  fHists->Add(fITSlayerPhi);
  fHists->Add(fEtaPhi);
  fHists->Add(fThetaRec);
  fHists->Add(fPhiRec);
  fHists->Add(fNumber);
  fHists->Add(fNumberAfterCut);
  fHists->Add(fVx);
  fHists->Add(fVy);
  fHists->Add(fVz);
  fHists->Add(fNVertexSPD);
  fHists->Add(fNVertexTracks);

  fHists->Add(fEtaPt);
  fHists->Add(fQPt);
  fHists->Add(fDca);
  fHists->Add(fDcaZ);

  fHists->Add(fDeltaPhiAll);
  fHists->Add(fDeltaPhiLeading);
  fHists->Add(fDiffDcaD);

  fHists->Add(fqRec);
  fHists->Add(fSigmaPtHist);

  fHists->Add(fRecPtPos);
  fHists->Add(fRecPtNeg);
  fHists->Add(fRecPhiPos);
  fHists->Add(fRecPhiNeg);
  fHists->Add(fRecEtaPos);
  fHists->Add(fRecEtaNeg);
  fHists->Add(fRecEtaPtPos);
  fHists->Add(fRecEtaPtNeg);
  fHists->Add(fRecDcaPos);
  fHists->Add(fRecDcaNeg);
  fHists->Add(fRecDcaNegInv);
  fHists->Add(fRecDPos);
  fHists->Add(fRecDNeg);


  fHists->Add(fRecQPtPosEta);
  fHists->Add(fRecQPtNegEta);
  fHists->Add(fRecPtPosEta);
  fHists->Add(fRecPtNegEta);
  fHists->Add(fRecPhiPosEta);
  fHists->Add(fRecPhiNegEta);
  fHists->Add(fRecDcaPosEta);
  fHists->Add(fRecDcaNegEta);
  fHists->Add(fRecDPosEta);
  fHists->Add(fRecDNegEta);

  for(Int_t i=0;i<4;i++){
    fHists->Add(fVertexX[i]);
    fHists->Add(fVertexY[i]);
    fHists->Add(fVertexZ[i]);
  }
  for(Int_t i=0;i<18;i++){
    fHists->Add(fRecPtTpcSector[i]);
    fHists->Add(fRecEtaTpcSector[i]);
    fHists->Add(fSignedDcaTpcSector[i]);
    fHists->Add(fRecQPtTpcSector[i]);
    fHists->Add(fRecEtaPtTpcSector[i]);
  }

  for(Int_t i=0;i<7;i++){
    fHists->Add(fRecPtPosLadder[i]);
    fHists->Add(fRecPtNegLadder[i]);
    fHists->Add(fRecPhiPosLadder[i]);
    fHists->Add(fRecPhiNegLadder[i]);
    fHists->Add(fRecEtaPosLadder[i]);
    fHists->Add(fRecEtaNegLadder[i]);
  }

  fHists->Add(fRecPtPosVz);
  fHists->Add(fRecPtNegVz);
  fHists->Add(fRecEtaPosVz);
  fHists->Add(fRecEtaNegVz);
  fHists->Add(fRecPhiPosVz);
  fHists->Add(fRecPhiNegVz);
  fHists->Add(fSignedDcaPosVz);
  fHists->Add(fSignedDcaNegVz);
  fHists->Add(fRecQPtPosEtaVz);
  fHists->Add(fRecQPtNegEtaVz);
  fHists->Add(fRecEtaPtPosVz);
  fHists->Add(fRecEtaPtNegVz);


  for(Int_t i=0;i<7;i++){
    fHists->Add(fSignDcaPos[i]);
    fHists->Add(fSignDcaNeg[i]);
    fHists->Add(fSignDcaNegInv[i]);
 
    fHists->Add(fPtSigmaPos[i]);
    fHists->Add(fPtSigmaNeg[i]);
    fHists->Add(fqPtRec[i]);
  
    fHists->Add(fDcaSigmaPos[i]);
    fHists->Add(fDcaSigmaNeg[i]);
 

  } 
  
  fHists->Add(fRecDcaPosPhi);
  fHists->Add(fRecDcaNegPhi);   
  fHists->Add(fRecPtPosPhi);
  fHists->Add(fRecPtNegPhi);   
  fHists->Add(fRecEtaPosPhi);
  fHists->Add(fRecEtaNegPhi);  
  fHists->Add(fRecQPtPhi);   
  fHists->Add(fRecEtaPtPosPhi);
  fHists->Add(fRecEtaPtNegPhi);   

  fHists->Add(fRecPtPosEtaPos);   
  fHists->Add(fRecPtNegEtaPos);   
  fHists->Add(fRecPtPosEtaNeg);   
  fHists->Add(fRecPtNegEtaNeg); 

  fHists->Add(fRec1PtPosEtaPos);   
  fHists->Add(fRec1PtNegEtaPos);   
  fHists->Add(fRec1PtPosEtaNeg);   
  fHists->Add(fRec1PtNegEtaNeg);   
 

  fHists->Add(fRecPhiPosEtaPos);   
  fHists->Add(fRecPhiNegEtaPos);   
  fHists->Add(fRecPhiPosEtaNeg);   
  fHists->Add(fRecPhiNegEtaNeg);   

  fHists->Add(fRecDcaPosPhiEtaPos);
  fHists->Add(fRecDcaNegPhiEtaPos);   
  fHists->Add(fRecPtPosPhiEtaPos);
  fHists->Add(fRecPtNegPhiEtaPos);  
  fHists->Add(fRecDcaPosPhiEtaNeg);
  fHists->Add(fRecDcaNegPhiEtaNeg);   
  fHists->Add(fRecPtPosPhiEtaNeg);
  fHists->Add(fRecPtNegPhiEtaNeg); 

  fHists->Add(fRecDcaPosPtEtaPos);
  fHists->Add(fRecDcaNegPtEtaPos);
  fHists->Add(fRecDcaPosPtEtaNeg);
  fHists->Add(fRecDcaNegPtEtaNeg);

  //  fHists->Add(fRecDcaPhiPtPosEtaPos); 
  //  fHists->Add(fRecDcaPhiPtPosEtaNeg); 
  //  fHists->Add(fRecDcaPhiPtNegEtaPos); 
  //  fHists->Add(fRecDcaPhiPtNegEtaNeg); 

  for(Int_t i=0;i<3;i++){
    for(Int_t j=0;j<2;j++){
      fHists->Add(fEtaBinPt[i][j]);
      fHists->Add(fPhiBinPt[i][j]);
      fHists->Add(fDcaBinPt[i][j]);
      fHists->Add(fEtaPhiBinPt[i][j]);
    }
  }

    
//   for (Int_t i=0; i<fHists->GetEntries(); ++i) {
//     TH1 *h1 = dynamic_cast<TH1*>(fHists->At(i));
//     if (h1){
//     //  Printf("%s ",h1->GetName());
//       h1->Sumw2();
//     }
//   }

  TH1::AddDirectory(oldStatus);
  PostData(1, fHists);

}

//__________________________________________________________

void AliAnalysisTaskQASym::UserExec(Option_t *) 
{
  // QA of global, TPC, ITS and ITS stand alone tracks
  // exploiting basic symmetries

  AliVEvent *event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }


  if(Entry()==0){
    AliESDEvent* esd = static_cast<AliESDEvent*>(event);
    if(esd){
      Printf("We are reading from ESD");
    }
   
  }



  if(fDebug>1)Printf("There are %d tracks in this event", event->GetNumberOfTracks());

  
  Int_t   leadingTrack  =   0;
  Float_t leadingEnergy = -20.;
  Float_t leadingPhi    =   0;//TMath::Pi();


  //check vertices
  AliESDEvent* esd = static_cast<AliESDEvent*>(event);
  Int_t nPileSPDVertices=1+esd->GetNumberOfPileupVerticesSPD(); // also SPD main vertex
  Int_t nPileTrkVertices=esd->GetNumberOfPileupVerticesTracks();
  fNVertexSPD->Fill(nPileSPDVertices);
  fNVertexTracks->Fill(nPileTrkVertices);

  //check primary vertex
  Float_t vx = 0;
  Float_t vy = 0;
  Float_t vz = 0;

  //primary vertex: contribution from different vertexers
  const AliVVertex* vertex = event->GetPrimaryVertex();
  if(!vertex) return;
  vx = vertex->GetX();
  vy = vertex->GetY();
  vz = vertex->GetZ();
  if(vertex->GetNContributors()>0){
    fVertexX[0]->Fill(vx);
    fVertexY[0]->Fill(vy);
    fVertexZ[0]->Fill(vz);     
  }
  

  
  const AliVVertex* vertexSPD = esd->GetPrimaryVertexSPD();
  if(vertexSPD){
    if(vertexSPD->GetNContributors()>0){
      fVertexX[1]->Fill(vertexSPD->GetX());
      fVertexY[1]->Fill(vertexSPD->GetY());
      fVertexZ[1]->Fill(vertexSPD->GetZ());
    }
  }

  const AliVVertex* vertexTrack = esd->GetPrimaryVertexTracks();
  if(vertexTrack){
    if(vertexTrack->GetNContributors()>0){
      fVertexX[2]->Fill(vertexTrack->GetX());
      fVertexY[2]->Fill(vertexTrack->GetY());
      fVertexZ[2]->Fill(vertexTrack->GetZ());
    }
  }

  const AliVVertex* vertexTPC = esd->GetPrimaryVertexTPC();
  if(vertexTPC){
    if(vertexTPC->GetNContributors()>0){
      fVertexX[3]->Fill(vertexTPC->GetX());
      fVertexY[3]->Fill(vertexTPC->GetY());
      fVertexZ[3]->Fill(vertexTPC->GetZ());
    }
  }

  //cuts on general vertex
  if(vertex->GetNContributors()<1) return;
  if (TMath::Abs(vz) > 10.) return;

  fNumber->Fill(event->GetNumberOfTracks());

  AliESDtrack *tpcP = 0x0;
  Int_t fNTracksAccepted=0;
  const Int_t arrSize = event->GetNumberOfTracks();
  Float_t * phiArray      = new Float_t[arrSize];
  Float_t * etaArray      = new Float_t[arrSize];
  Float_t * ptArray       = new Float_t[arrSize];
  Float_t * dcaArray      = new Float_t[arrSize];
  Int_t   * chargeArray   = new Int_t  [arrSize];
  Bool_t  * acceptedArray = new Bool_t [arrSize];

  for (Int_t i = 0; i < event->GetNumberOfTracks(); i++) {
    phiArray[i]     = 0.;
    etaArray[i]     = 0.;
    ptArray[i]      = 0.;
    dcaArray[i]     = 0.;
    chargeArray[i]  = 0;
    acceptedArray[i]= kFALSE;
    
  }



  for (Int_t iTrack = 0; iTrack < event->GetNumberOfTracks(); iTrack++) {
    
    //prevent mem leak for TPConly track
    if(fTrackType==2&&tpcP){
      delete tpcP;
      tpcP = 0;
    }

    AliVParticle *track = event->GetTrack(iTrack);
    AliESDtrack *esdtrack =  static_cast<AliESDtrack*>(track);
    esdtrack->PropagateToDCA(event->GetPrimaryVertex(),
			     event->GetMagneticField(), 10000.);

    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    //__________
    // run Task for global tracks or ITS tracks or TPC tracks
    const AliExternalTrackParam *tpcPin = 0x0;
    Double_t phiIn=0.;

    if(fTrackType==0){
      //Fill all histograms with global tracks
      tpcP = esdtrack;
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      phiIn = tpcP->Phi();
    }
    else if(fTrackType==1){
      //Fill all histograms with ITS tracks
      tpcP = esdtrack;
      phiIn = tpcP->Phi();
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      if(fStandAlone==kTRUE) {
	if(!(tpcP->GetStatus()&AliESDtrack::kITSpureSA))continue;
      }
      else if(fStandAlone==kFALSE){
	if(tpcP->GetStatus()&AliESDtrack::kITSpureSA)continue;
      }
    }
    else if(fTrackType==2){     
      //Fill all histograms with TPC track information
      tpcPin = esdtrack->GetInnerParam();
      if (!tpcPin) continue;
      phiIn=tpcPin->Phi();

      tpcP = AliESDtrackCuts::GetTPCOnlyTrack(static_cast<AliESDEvent*>(event),esdtrack->GetID());
      if (!tpcP) continue;
      if (!fCuts->AcceptTrack(tpcP)) continue;
      if(tpcP->GetNcls(1)>160)continue;//jacek's track cut
      if(tpcP->GetConstrainedChi2TPC()<0)continue; // jacek's track cut
    }
    else{
      Printf("ERROR: wrong track type \n");
      continue;
    }
    //___________
    //
  

    fNTracksAccepted++;
    phiArray[iTrack]     = phiIn;
    etaArray[iTrack]     = tpcP->Eta();
    ptArray[iTrack]      = tpcP->Pt();
    chargeArray[iTrack]  = tpcP->Charge();
    acceptedArray[iTrack]= kTRUE;

 
    if(tpcP->E()>leadingEnergy){
      leadingTrack=iTrack;
      leadingEnergy=tpcP->E();
      leadingPhi=phiIn;
    }
   
    
    fqRec->Fill(tpcP->Charge());
  

    Double_t fSigmaPt = tpcP->GetSigma1Pt2();
    fSigmaPt= sqrt(fSigmaPt);
    fSigmaPt= fSigmaPt *(tpcP->Pt()*tpcP->Pt()); 

    if(TMath::Abs(fSigmaPt) < 1.e-10) continue;

    fSigmaPtHist->Fill(TMath::Log10(fSigmaPt));
 

    // hits in ITS layer
    Int_t cas=-1;
    if(tpcP->HasPointOnITSLayer(0)) 
      cas=0;
    else if(!tpcP->HasPointOnITSLayer(0)
	    &&  tpcP->HasPointOnITSLayer(1)) 
      cas=1;
    else if(!tpcP->HasPointOnITSLayer(0)
	    && !tpcP->HasPointOnITSLayer(1) 
	    &&  tpcP->HasPointOnITSLayer(2)) 
      cas=2;
    else if(!tpcP->HasPointOnITSLayer(0)
	    && !tpcP->HasPointOnITSLayer(1) 
	    && !tpcP->HasPointOnITSLayer(2)
	    &&  tpcP->HasPointOnITSLayer(3)) 
      cas=3;
    else if(!tpcP->HasPointOnITSLayer(0)
	    && !tpcP->HasPointOnITSLayer(1) 
	    && !tpcP->HasPointOnITSLayer(2)
	    && !tpcP->HasPointOnITSLayer(3)
	    &&  tpcP->HasPointOnITSLayer(4)) 
      cas=4;
    else if(   !tpcP->HasPointOnITSLayer(0)
	       && !tpcP->HasPointOnITSLayer(1)
	       && !tpcP->HasPointOnITSLayer(2)
	       && !tpcP->HasPointOnITSLayer(3)
	       && !tpcP->HasPointOnITSLayer(4) 
	       &&  tpcP->HasPointOnITSLayer(5)) 
      cas=5;
    else 
      cas=6;
  
   
   
    //------------------- 
    Float_t fXVertexCor = 0.;
    Float_t fYVertexCor = 0.;

    fXVertexCor = tpcP->Xv() - vertex->GetX(); // coordinate corrected for vertex position
    fYVertexCor = tpcP->Yv() - vertex->GetY(); // "
    Double_t fSignedDca = (tpcP->Py()*fXVertexCor - tpcP->Px()*fYVertexCor)/tpcP->Pt();


    fqPtRec[cas]->Fill(tpcP->Charge()/tpcP->Pt());
    
    

    fHistRECpt->Fill(tpcP->Pt());
    fEta->Fill(tpcP->Eta());
    fEtavPt->Fill(tpcP->Eta(), TMath::Log(tpcP->Pt()));
    fPhivPt->Fill(phiIn, TMath::Log(tpcP->Pt()));
    fEtaPhi->Fill(tpcP->Eta(), phiIn);
    fThetaRec->Fill(tpcP->Theta());
    fPhiRec->Fill(phiIn);
    fVx->Fill(tpcP->Xv());
    fVy->Fill(tpcP->Yv());
    fVz->Fill(tpcP->Zv());
  

    fEtaPt->Fill(tpcP->Eta()/tpcP->Pt());
    fQPt->Fill(tpcP->Charge()/tpcP->Pt());
    fDca->Fill(fSignedDca);
    dcaArray[iTrack]=fSignedDca;
    fRecQPtPhi->Fill(tpcP->Charge()/tpcP->Pt(), phiIn);

    Float_t fXY = 0.;
    Float_t  fZ = 0.;

    tpcP->GetImpactParameters(fXY,fZ);
    fDiffDcaD->Fill(fSignedDca+fXY);
    fDcaZ->Fill(fZ);
    
    if(fTrackType==2) fCompareTPCparam->Fill(fZ,tpcPin->GetTgl());

    if(fTrackType!=2){//for global and ITS tracks
      for(Int_t itsLayer=0;itsLayer<6;itsLayer++){
	if(tpcP->HasPointOnITSLayer(itsLayer)){
	  fITSlayer->Fill(itsLayer);
	  fITSlayerEta->Fill(itsLayer, tpcP->Eta());
	  fITSlayerPhi->Fill(itsLayer, tpcP->Phi());
	}
      }    
    }

    //for positive particles
    
    if(tpcP->Charge()>0){
      fRecPtPos->Fill(tpcP->Pt());
      fRecPtPosLadder[cas]->Fill(tpcP->Pt());
      fRecPtPosVz->Fill(TMath::Log10(tpcP->Pt()),tpcP->Zv());
      fRecPhiPos->Fill(TMath::RadToDeg()*phiIn);
      
      
      fRecPhiPosLadder[cas]->Fill(TMath::RadToDeg()*phiIn);
      fRecPhiPosVz->Fill(TMath::RadToDeg()*phiIn,tpcP->Zv());
      fSignedDcaPosVz->Fill(fSignedDca,tpcP->Zv());

      fRecEtaPos->Fill(tpcP->Eta());
      fRecEtaPosLadder[cas]->Fill(tpcP->Eta());
      fRecEtaPtPos->Fill(tpcP->Eta()/tpcP->Pt());
      fRecEtaPosVz->Fill(tpcP->Eta(),tpcP->Zv());
      fRecEtaPtPosVz->Fill(tpcP->Eta()/tpcP->Pt(),tpcP->Zv());
     
      fRecDcaPos->Fill(fSignedDca);
      fRecDcaPosPhi->Fill(fSignedDca, phiIn);
      fRecPtPosPhi->Fill(TMath::Log10(tpcP->Pt()), phiIn);
      fRecEtaPtPosPhi->Fill(tpcP->Eta()/tpcP->Pt(), phiIn);
      fRecEtaPosPhi->Fill(tpcP->Eta(), phiIn);
      fRecDPos->Fill(fXY);
      fSignDcaPos[cas]->Fill(fSignedDca);
    
     
      fDcaSigmaPos[cas]->Fill(fSignedDca, TMath::Log10(fSigmaPt));
    
      fPtSigmaPos[cas]->Fill(TMath::Log10(fSigmaPt));
      //pos eta
      if(tpcP->Eta()>0){
	fRecPtPosEtaPos->Fill(tpcP->Pt());
	fRec1PtPosEtaPos->Fill(1/tpcP->Pt());
	fRecPhiPosEtaPos->Fill(phiIn);
	fRecDcaPosPhiEtaPos->Fill(fSignedDca, phiIn);
	fRecDcaPosPtEtaPos->Fill(fSignedDca, TMath::Log10(tpcP->Pt()));
	fRecPtPosPhiEtaPos->Fill(TMath::Log10(tpcP->Pt()), phiIn);
	//fRecDcaPhiPtPosEtaPos->Fill(phiIn, tpcP->Pt(), fSignedDca);
      }
      //neg eta
      else{
	fRecPtPosEtaNeg->Fill(tpcP->Pt());
	fRec1PtPosEtaNeg->Fill(1/tpcP->Pt());
	fRecPhiPosEtaNeg->Fill(phiIn);
	fRecDcaPosPhiEtaNeg->Fill(fSignedDca, phiIn);
	fRecDcaPosPtEtaNeg->Fill(fSignedDca, TMath::Log10(tpcP->Pt()));
	fRecPtPosPhiEtaNeg->Fill(TMath::Log10(tpcP->Pt()), phiIn);
	//fRecDcaPhiPtPosEtaNeg->Fill(phiIn, tpcP->Pt(), fSignedDca);
      }
      
    }
    //and negative particles
    else {
      fRecPtNeg->Fill(tpcP->Pt());
      fRecPtNegLadder[cas]->Fill(tpcP->Pt());
      fRecPtNegVz->Fill(TMath::Log10(tpcP->Pt()),tpcP->Zv());
           
      fRecPhiNeg->Fill(TMath::RadToDeg()*phiIn);
      fRecPhiNegLadder[cas]->Fill(TMath::RadToDeg()*phiIn);
      fRecPhiNegVz->Fill(TMath::RadToDeg()*phiIn,tpcP->Zv());
      fSignedDcaNegVz->Fill(fSignedDca,tpcP->Zv());
      fRecEtaPtNegVz->Fill(tpcP->Eta()/tpcP->Pt(),tpcP->Zv());

      fRecEtaNeg->Fill(tpcP->Eta());
      fRecEtaNegLadder[cas]->Fill(tpcP->Eta());
      fRecEtaPtNeg->Fill(tpcP->Eta()/tpcP->Pt());
      fRecEtaNegVz->Fill(tpcP->Eta(),tpcP->Zv());
     
      fRecDcaNeg->Fill(fSignedDca);
      fRecDcaNegInv->Fill(-fSignedDca);
      fRecDcaNegPhi->Fill(fSignedDca, phiIn);
      fRecPtNegPhi->Fill(TMath::Log10(tpcP->Pt()), phiIn);
      fRecEtaNegPhi->Fill(tpcP->Eta(), phiIn);
      fRecEtaPtNegPhi->Fill(tpcP->Eta()/tpcP->Pt(), phiIn);
      fRecDNeg->Fill(fXY);
      fSignDcaNeg[cas]->Fill(fSignedDca);
      fSignDcaNegInv[cas]->Fill(-fSignedDca);
     
     
      fDcaSigmaNeg[cas]->Fill(fSignedDca,TMath::Log10(fSigmaPt));
   
      fPtSigmaNeg[cas]->Fill(TMath::Log10(fSigmaPt));
      
      //pos eta
      if(tpcP->Eta()>0){
	fRecPtNegEtaPos->Fill(tpcP->Pt());
	fRec1PtNegEtaPos->Fill(1/tpcP->Pt());
	fRecPhiNegEtaPos->Fill(phiIn);
	fRecDcaNegPhiEtaPos->Fill(fSignedDca, phiIn);
	fRecDcaNegPtEtaPos->Fill(fSignedDca, TMath::Log10(tpcP->Pt()));
	fRecPtNegPhiEtaPos->Fill(TMath::Log10(tpcP->Pt()), phiIn);
	//fRecDcaPhiPtNegEtaPos->Fill(phiIn, tpcP->Pt(), fSignedDca);
      }
      //neg eta
      else{
	fRecPtNegEtaNeg->Fill(tpcP->Pt());
	fRec1PtNegEtaNeg->Fill(1/tpcP->Pt());
	fRecPhiNegEtaNeg->Fill(phiIn);
	fRecDcaNegPhiEtaNeg->Fill(fSignedDca, phiIn);
	fRecDcaNegPtEtaNeg->Fill(fSignedDca, TMath::Log10(tpcP->Pt()));
	fRecPtNegPhiEtaNeg->Fill(TMath::Log10(tpcP->Pt()), phiIn);
	//fRecDcaPhiPtNegEtaNeg->Fill(phiIn, tpcP->Pt(), fSignedDca);
      }

    }
    


    //all particles with positive eta
    if(tpcP->Eta()>0){
      fRecQPtPosEta->Fill(tpcP->Charge()/tpcP->Pt());
      fRecPtPosEta->Fill(tpcP->Pt());
      fRecPhiPosEta->Fill(TMath::RadToDeg()*phiIn);
      fRecQPtPosEtaVz->Fill(tpcP->Charge()/tpcP->Pt(),tpcP->Zv());
      fRecDcaPosEta->Fill(fSignedDca);
      fRecDPosEta->Fill(fXY);
    }
    //all particles with negative eta (and eta==0)
    else{
      fRecQPtNegEta->Fill(tpcP->Charge()/tpcP->Pt());
      fRecPtNegEta->Fill(tpcP->Pt());
      fRecPhiNegEta->Fill(TMath::RadToDeg()*phiIn);
      fRecQPtNegEtaVz->Fill(tpcP->Charge()/tpcP->Pt(),tpcP->Zv());
      fRecDcaNegEta->Fill(fSignedDca);
      fRecDNegEta->Fill(fXY);

    }


    fRecPtTpcSector[Int_t(phiIn*
			  TMath::RadToDeg()/20)]->Fill(tpcP->Pt());
    fRecEtaTpcSector[Int_t(phiIn*
			   TMath::RadToDeg()/20)]->Fill(tpcP->Eta());
    fSignedDcaTpcSector[Int_t(phiIn*
			      TMath::RadToDeg()/20)]->Fill(fSignedDca); 
    fRecQPtTpcSector[Int_t(phiIn*
			   TMath::RadToDeg()/20)]->Fill(tpcP->Charge()/tpcP->Pt());
    fRecEtaPtTpcSector[Int_t(phiIn*
			     TMath::RadToDeg()/20)]->Fill(tpcP->Eta()/tpcP->Pt());
     


//     // another track loop
//     for (Int_t iTrack2 = 0; iTrack2 < event->GetNumberOfTracks(); iTrack2++) {
      
//       if(LeadingTrack==iTrack2) continue;

//       AliVParticle *track2 = event->GetTrack(iTrack2);
//       AliESDtrack* esdtrack2 =  static_cast<AliESDtrack*>(track2);
//       if (!track2) {
// 	Printf("ERROR: Could not receive track %d", iTrack);
// 	continue;
//       }
//       if (!fCuts->AcceptTrack(esdtrack2)) continue;
//       //propagate to dca
//       esdtrack2->PropagateToDCA(event->GetPrimaryVertex(),
// 				event->GetMagneticField(), 10000.);
 
//       fDeltaPhiLeading->Fill((LeadingPhi-esdtrack2->Phi())*TMath::RadToDeg(),
// 			     LeadingPhi*TMath::RadToDeg() );

     

//     }//second track loop

    // if(fTrackType==2) delete tpcP; // delete in case of TPCOnlyTrack

  }//first track loop

  fNumberAfterCut->Fill(fNTracksAccepted);
  
  //second track loop
 
  for (Int_t iT = 0; iT < event->GetNumberOfTracks(); iT++) {
    if(acceptedArray[iT]){
      if(ptArray[iT]>0.2 && ptArray[iT]<1. ){
	fEtaBinPt[0][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT]);
	fDcaBinPt[0][Bool_t(chargeArray[iT]>0)]->Fill(dcaArray[iT]);
	fPhiBinPt[0][Bool_t(chargeArray[iT]>0)]->Fill(phiArray[iT]);
	fEtaPhiBinPt[0][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT], phiArray[iT]);
      }
      else if(ptArray[iT]>1. && ptArray[iT]<5.){
	fEtaBinPt[1][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT]);
	fDcaBinPt[1][Bool_t(chargeArray[iT]>0)]->Fill(dcaArray[iT]);
	fPhiBinPt[1][Bool_t(chargeArray[iT]>0)]->Fill(phiArray[iT]);
	fEtaPhiBinPt[1][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT], phiArray[iT]);
      }
      else if (ptArray[iT]>5.){
	fEtaBinPt[2][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT]);
	fDcaBinPt[2][Bool_t(chargeArray[iT]>0)]->Fill(dcaArray[iT]);
	fPhiBinPt[2][Bool_t(chargeArray[iT]>0)]->Fill(phiArray[iT]);
	fEtaPhiBinPt[2][Bool_t(chargeArray[iT]>0)]->Fill(etaArray[iT], phiArray[iT]);
      }

      if(fNTracksAccepted>=fLow&&fNTracksAccepted<=fHigh){
	fEtaWidth->Fill(etaArray[iT]);
	fPhiWidth->Fill(phiArray[iT]);
	fDcaWidth->Fill(dcaArray[iT]);
	fPtWidth->Fill(ptArray[iT]);
       }
     }
  }


  //prevent mem leak for TPConly track
  if(fTrackType==2&&tpcP){
    delete tpcP;
    tpcP = 0;
  }

  if(phiArray){
    delete[] phiArray; 
    phiArray=0; 
  }
  
  if(etaArray){
    delete[] etaArray; 
    etaArray=0; 
  }
  
  if(dcaArray){
    delete[] dcaArray; 
    dcaArray=0; 
  }
  
  if(ptArray){
    delete[] ptArray; 
    ptArray=0; 
  }
  
  if(chargeArray){
    delete[] chargeArray; 
    chargeArray=0; 
  }
  
  if(acceptedArray){
    delete[] acceptedArray; 
    acceptedArray=0; 
  }
  
  // Post output data.
  // PostData(1, fHistPt);
  PostData(1, fHists);
}      





//________________________________________________________________________
void AliAnalysisTaskQASym::Terminate(Option_t *) 
{


}  





