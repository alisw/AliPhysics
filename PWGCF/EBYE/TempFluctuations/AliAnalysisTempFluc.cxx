/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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


#include "AliAnalysisTempFluc.h"
#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TAxis.h"
#include "TFile.h"
#include "TString.h"
#include "TGrid.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliAODHeader.h"
#include "AliPID.h"
//-------------------
#include "TParticlePDG.h"                          // for MC only
#include "AliAODMCParticle.h"
#include "TClonesArray.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "TSpline.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"

ClassImp(AliAnalysisTempFluc)

//________________________________________________________________________
AliAnalysisTempFluc::AliAnalysisTempFluc() 
  :AliAnalysisTaskSE(),
   fOutput(0),

fCentralityEstimator("V0M"),
fCentralityCutL(0.0),
fCentralityCutH(5.0),
AnalysisMode("data"),                   // "data" ? "MC"
DataType("AOD"),
mass_pion(0.1395),
fAOD(0x0),
fPIDResponse(0x0),
fHistPtTPC(0),
fHistPtTruthdummy(0),
fHistTruth(0),
fHistRecTPC(0),
fHistRecPrimaryTPC(0),
fHistRecSecWDTPC(0),
fHistRecSecMatTPC(0),
fHistRecMisIdTPC(0),
fHistPM(0),
fHistdEdx(0),
fHistdEdxSigma(0),
fNSigmaTPC(0),
fNSigmaTPCCut(0),
fNSigmaTOF(0),
fNSigmaTPCTOF(0),
fNSigmaTPCTOFCut(0),
fHistPtTPCTOF(0),
fHistRecTPCTOF(0),
fHistRecPrimaryTPCTOF(0),
fHistRecSecWDTPCTOF(0),
fHistRecSecMatTPCTOF(0),
fHistRecMisIdTPCTOF(0),
fHistTOF(0),
mcTruth(0),
fHistPtTruth(0),
fsizearray(0),
//-------------------
fEffptTPC(0),
fEffptTOF(0),
fEffptTPCTOF(0),
fcorrection(kFALSE),
fBinningptL(""),
fBinningptH(""),
fefffilename(""),
//---------------------

DCAxy(0),
DCAz(0),
gEtaL (0),
gEtaH (0),
gRapL (0),
gRapH (0),
fPtL(0),
fPtH(0),
fBit(272),
fCutNSigmaTPC(3.0),
fCutNSigmaTPCTOF(3.0),
   fNEvt(0)
{
    for(Int_t iQA=0;iQA<38;iQA++){
        fQAHist[iQA]=0x0;
    }
    for(Int_t j=0;j<10;j++){
        MeanpTebyeTPC[j]=0;
        Teff_ebyeTPC[j]=0;
        No_pion_evt_TPC[j]=0;
        MeanpTebyeTPCTOF[j]=0;
        Teff_ebyeTPCTOF[j]=0;
        No_pion_evt_TPCTOF[j]=0;
        MeanpTebyeTruth[j]=0;
        Teff_ebyeTruth[j]=0;
        No_pion_evt_Truth[j]=0;
    }


}

//________________________________________________________________________
AliAnalysisTempFluc::AliAnalysisTempFluc(const char *name) 
   :AliAnalysisTaskSE(name),
    fOutput(0),

fCentralityEstimator("V0M"),
fCentralityCutL(0.0),
fCentralityCutH(5.0),
AnalysisMode("mc"),                   // "data" ? "MC"
DataType("AOD"),
mass_pion(0.1395),
fAOD(0x0),
fPIDResponse(0x0),
fHistPtTPC(0),
fHistPtTruthdummy(0),
fHistTruth(0),
fHistRecTPC(0),
fHistRecPrimaryTPC(0),
fHistRecSecWDTPC(0),
fHistRecSecMatTPC(0),
fHistRecMisIdTPC(0),
fHistPM(0),
fHistdEdx(0),
fHistdEdxSigma(0),
fNSigmaTPC(0),
fNSigmaTPCCut(0),
fHistTOF(0),
mcTruth(0),
fNSigmaTOF(0),
fNSigmaTPCTOF(0),
fNSigmaTPCTOFCut(0),
fHistPtTPCTOF(0),
fHistRecTPCTOF(0),
fHistRecPrimaryTPCTOF(0),
fHistRecSecWDTPCTOF(0),
fHistRecSecMatTPCTOF(0),
fHistRecMisIdTPCTOF(0),
fHistPtTruth(0),
fsizearray(0),
//-------------------
fEffptTPC(0),
fEffptTOF(0),
fEffptTPCTOF(0),
fcorrection(kFALSE),
fBinningptL(""),
fBinningptH(""),
fefffilename(""),
//---------------------
DCAxy(0),
DCAz(0),
gEtaL (0),
gEtaH (0),
gRapL (0),
gRapH (0),
fPtL(0),
fPtH(0),
fBit(272),
fCutNSigmaTPC(3.0),
fCutNSigmaTPCTOF(3.0),
fNEvt(0)
{
    for(Int_t iQA=0;iQA<38;iQA++){
        fQAHist[iQA]=0x0;
    }
    for(Int_t j=0;j<10;j++){
        MeanpTebyeTPC[j]=0;
        Teff_ebyeTPC[j]=0;
        No_pion_evt_TPC[j]=0;
        MeanpTebyeTPCTOF[j]=0;
        Teff_ebyeTPCTOF[j]=0;
        No_pion_evt_TPCTOF[j]=0;
        MeanpTebyeTruth[j]=0;
        Teff_ebyeTruth[j]=0;
        No_pion_evt_Truth[j]=0;

    }
  
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTempFluc::~AliAnalysisTempFluc()
{
 
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
}

//________________________________________________________________________
void AliAnalysisTempFluc::UserCreateOutputObjects()
{
   
           
    fOutput = new TList();
    fOutput->SetOwner();  // IMPORTANT!
    
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
    Int_t ptbins = 60;
    fHistPtTPC = new TH1F("fHistPtTPC", "P_{T} distribution for reconstructed TPC", ptbins, fPtL, fPtH);
    fHistPtTPC->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtTPC->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtTPC->Sumw2();

    fHistPtTruthdummy = new TH1F("fHistPtTruthdummy", "P_{T} distribution for Truth ebye", ptbins, fPtL, fPtH);
    fHistPtTruthdummy->Sumw2();
    
    fHistTruth = new TH3F("HistTruth","Rapidity-Phi-Pt for Truth",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistTruth->GetXaxis()->SetTitle("Rapidity");
    fHistTruth->GetYaxis()->SetTitle("#phi");
    fHistTruth->GetZaxis()->SetTitle("p_{T}");
    fHistTruth->Sumw2();
    
    fHistRecTPC = new TH3F("HistRecTPC","Rapidity-Phi-Pt for Reconstructed for TPC Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecTPC->GetXaxis()->SetTitle("Rapidity");
    fHistRecTPC->GetYaxis()->SetTitle("#phi");
    fHistRecTPC->GetZaxis()->SetTitle("p_{T}");
    fHistRecTPC->Sumw2();
    
    fHistRecPrimaryTPC = new TH3F("HistRecTPCPrimary","Rapidity-Phi-Pt for Reconstructed for TPC Primary Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecPrimaryTPC->GetXaxis()->SetTitle("Rapidity");
    fHistRecPrimaryTPC->GetYaxis()->SetTitle("#phi");
    fHistRecPrimaryTPC->GetZaxis()->SetTitle("p_{T}");
    fHistRecPrimaryTPC->Sumw2();
    
    fHistRecSecWDTPC = new TH3F("HistRecTPCSecWD","Rapidity-Phi-Pt for Reconstructed for TPC Secondary WD Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecSecWDTPC->GetXaxis()->SetTitle("Rapidity");
    fHistRecSecWDTPC->GetYaxis()->SetTitle("#phi");
    fHistRecSecWDTPC->GetZaxis()->SetTitle("p_{T}");
    
    
    fHistRecSecMatTPC = new TH3F("HistRecTPCSecMat","Rapidity-Phi-Pt for Reconstructed for TPC Secondary Material Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecSecMatTPC->GetXaxis()->SetTitle("Rapidity");
    fHistRecSecMatTPC->GetYaxis()->SetTitle("#phi");
    fHistRecSecMatTPC->GetZaxis()->SetTitle("p_{T}");

    
    fHistRecMisIdTPC = new TH3F("HistRecTPCMisID","Rapidity-Phi-Pt for Reconstructed for TPC MisIdentified Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecMisIdTPC->GetXaxis()->SetTitle("Rapidity");
    fHistRecMisIdTPC->GetYaxis()->SetTitle("#phi");
    fHistRecMisIdTPC->GetZaxis()->SetTitle("p_{T}");
    
    
    
    
    fHistPtTPCTOF = new TH1F("fHistPtTPCTOF", "P_{T} distribution for reconstructed TPC & TOF", ptbins, fPtL, fPtH);
    fHistPtTPCTOF->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtTPCTOF->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtTPCTOF->Sumw2();
    
    fHistRecTPCTOF = new TH3F("HistRecTPCTOF","Rapidity-Phi-Pt for Reconstructed TPC & TOF",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecTPCTOF->GetXaxis()->SetTitle("Rapidity");
    fHistRecTPCTOF->GetYaxis()->SetTitle("#phi");
    fHistRecTPCTOF->GetZaxis()->SetTitle("p_{T}");
    fHistRecTPCTOF->Sumw2();
    
    fHistRecPrimaryTPCTOF = new TH3F("HistRecTPCTOFPrimary","Rapidity-Phi-Pt for Reconstructed for TPCTOF Primary Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecPrimaryTPCTOF->GetXaxis()->SetTitle("Rapidity");
    fHistRecPrimaryTPCTOF->GetYaxis()->SetTitle("#phi");
    fHistRecPrimaryTPCTOF->GetZaxis()->SetTitle("p_{T}");
    fHistRecPrimaryTPCTOF->Sumw2();
    
    fHistRecSecWDTPCTOF = new TH3F("HistRecTPCTOFSecWD","Rapidity-Phi-Pt for Reconstructed for TPCTOF Secondary WD Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecSecWDTPCTOF->GetXaxis()->SetTitle("Rapidity");
    fHistRecSecWDTPCTOF->GetYaxis()->SetTitle("#phi");
    fHistRecSecWDTPCTOF->GetZaxis()->SetTitle("p_{T}");
    
    
    fHistRecSecMatTPCTOF = new TH3F("HistRecTPCTOFSecMat","Rapidity-Phi-Pt for Reconstructed for TPCTOF Secondary Material Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecSecMatTPCTOF->GetXaxis()->SetTitle("Rapidity");
    fHistRecSecMatTPCTOF->GetYaxis()->SetTitle("#phi");
    fHistRecSecMatTPCTOF->GetZaxis()->SetTitle("p_{T}");
   
    
    fHistRecMisIdTPCTOF = new TH3F("HistRecTPCTOFMisID","Rapidity-Phi-Pt for Reconstructed for TPCTOF MisIdentified Only",100,-0.5,0.5,360,0,6.28,ptbins, fPtL, fPtH);
    fHistRecMisIdTPCTOF->GetXaxis()->SetTitle("Rapidity");
    fHistRecMisIdTPCTOF->GetYaxis()->SetTitle("#phi");
    fHistRecMisIdTPC->GetZaxis()->SetTitle("p_{T}");
    

        
        
    // NEW HISTO should be defined here, with a sensible name,
        
    fHistPM = new TH2F("hHistPlusMinus","Plus ~ minus", 3000,0,3000,3000,0,3000);
    
    fNSigmaTPC=new TH2F("fNSigmaTPC","NSigma - pT distribution for TPC",60,fPtL,fPtH,1000,-10,90);
    fNSigmaTPC->GetXaxis()->SetTitle("#p_{T}");
    fNSigmaTPC->GetYaxis()->SetTitle("N#sigma TPC");
   
    
    fNSigmaTPCCut=new TH2F("fNSigmaTPCCut","NSigma - pT distribution for TPC after N#sigma Cut",60,fPtL,fPtH,1000,-10,90);
    fNSigmaTPCCut->GetXaxis()->SetTitle("#p_{T}");
    fNSigmaTPCCut->GetYaxis()->SetTitle("N#sigma TPC");
   

    fHistTOF=new TH2F("fHistTOF", "TOF_signal momentum Vs #beta",900,fPtL,fPtH,500,0.1,1.1);
    fHistTOF->GetXaxis()->SetTitle("momentum");
    fHistTOF->GetYaxis()->SetTitle("#beta p(GeV/C)");

    
    fNSigmaTOF=new TH2F("fNSigmaTOF","NSigma - pT distribution for TOF",60,fPtL,fPtH,1000,-10,90);
    fNSigmaTPCCut->GetXaxis()->SetTitle("#p_{T}");
    fNSigmaTPCCut->GetYaxis()->SetTitle("N#sigma TOF");
  
    
    fNSigmaTPCTOF=new TH2F("fNSigmaTPCTOF","NSigma - pT distribution for TPC + TOF",60,fPtL,fPtH,100,-2,2);
    fNSigmaTPCTOF->GetXaxis()->SetTitle("#p_{T}");
    fNSigmaTPCTOF->GetYaxis()->SetTitle("N#sigma TPC+TOF");
    
    
    fNSigmaTPCTOFCut=new TH2F("fNSigmaTPCTOFCut","NSigma - pT distribution for TPC + TOF after cut",60,fPtL,fPtH,1000,-10,90);
    fNSigmaTPCTOFCut->GetXaxis()->SetTitle("#p_{T}");
    fNSigmaTPCTOFCut->GetYaxis()->SetTitle("N#sigma TPC+TOF");
    
    
    fHistdEdx = new TH2F("fHistdEdx", "Bethe Bloch for TPC", 1000,-5.0,5.0,1000,-2,500);
    fHistdEdx->GetXaxis()->SetTitle("rigidity");
    fHistdEdx->GetYaxis()->SetTitle("Central dE/dx");
    
   
   fHistdEdxSigma = new TH2F("fHistdEdxSigma", "Bethe Bloch for TPC after NSigma cut", 1000,-5.0,5.0,1000,-2,500);
   fHistdEdxSigma->GetXaxis()->SetTitle("rigidity");
  fHistdEdxSigma->GetYaxis()->SetTitle("Central dE/dx");
    
    fHistPtTruth = new TH1F("fHistPtTruth", "P_{T} distribution for Truth", ptbins, fPtL, fPtH);
    fHistPtTruth->Sumw2();
    
    fQAHist[0] = new TH1F("DCAxy","no of tracks after DCAxy Cut", 300, -15., 15.);               // For TPC only track DCA , for Global track (1024) ---> return nothing
    fQAHist[1] = new TH1F("DCAz","no of tracks after DCAz cut", 300, -15., 15.);
    fQAHist[2] = new TH1F("TPC clusters","no tracks after TPCcluster cut",100,80,180);
    fQAHist[3] = new TH1F("chisq/ndf","no tracks after chisq/ndf cut",500,0,5);
    fQAHist[4] = new TH1F("Event Counter","no of events after cuts",20,-5,15);
    fQAHist[5] = new TH1F("vertex_cut_X","Distribution of X vertex cut",200,-1,1);
    fQAHist[6] = new TH1F("vertex_cut_Y","Distribution of Y vertex cut",200,-1,1);
    fQAHist[7] = new TH1F("vertex_cut_Z","Distribution of Z vertex cut",300,-15,15);
    fQAHist[8] = new TH1F("Eta","Eta Distribution of all charge particles before cut",200,-2.0,2.0);
    fQAHist[9] = new TH1F("Phi","Phi Distribution of all charge particles before cut",1000,-0.5,15.0);
    fQAHist[10] = new TH1F("EtaPostive","Eta Distribution of +ve charge particles before cut",200,-2.0,2.0);
    fQAHist[11] = new TH1F("PhiPositive","Phi Distribution of +ve charge particles before cut",1000,-0.5,15.0);
    fQAHist[12] = new TH1F("EtaNegative","Eta Distribution of -ve charge particles before cut",200,-2.0,2.0);
    fQAHist[13] = new TH1F("PhiNegative","Phi Distribution of -ve charge particles before cut",1000,-0.5,15.0);
    fQAHist[14] = new TH1F("EtaptCut","Eta Distribution of all charge particles after pT cut",200,-2.0,2.0);
    fQAHist[15] = new TH1F("PhiptCut","Phi Distribution of all charge particles after pT cut",1000,-0.5,15.0);
    fQAHist[16] = new TH1F("EtaptCutPositive","Eta Distribution of +ve charge particles after pT cut",200,-2.0,2.0);
    fQAHist[17] = new TH1F("PhiptCutPositive","Phi Distribution of +ve charge particles after pT cut",1000,-0.5,15.0);
    fQAHist[18] = new TH1F("EtaptCutNegative","Eta Distribution of -ve charge particles after pT cut",200,-2.0,2.0);
    fQAHist[19] = new TH1F("PhiptCutNegative","Phi Distribution of -ve charge particles after pT cut",1000,-0.5,15.0);
    fQAHist[20] = new TH1F("EtaptetaCut","Eta Distribution of all charge particles after pT & eta cut",200,-2.0,2.0);
    fQAHist[21] = new TH1F("PhiptetaCut","Phi Distribution of all charge particles after pT & eta cut",1000,-0.5,15.0);
    fQAHist[22] = new TH1F("EtaptetaCutPositive","Eta Distribution of +ve charge particles after pT & eta cut",200,-2.0,2.0);
    fQAHist[23] = new TH1F("PhiptetaCutPositive","Phi Distribution of +ve charge particles after pT & eta cut",1000,-0.5,15.0);
    fQAHist[24] = new TH1F("EtaptetaCutNegative","Eta Distribution of -ve charge particles after pT & eta cut",200,-2.0,2.0);
    fQAHist[25] = new TH1F("PhiptetaCutNegative","Phi Distribution of -ve charge particles after pT & eta cut",1000,-0.5,15.0);

    
    fQAHist[26] = new TH1F("Eta_pion_TPC","Eta Distribution of all Pions from TPC",200,-2.0,2.0);
    fQAHist[27] = new TH1F("Phi_pion_TPC","Phi Distribution of all Pions from TPC",1000,-0.5,15.0);
    fQAHist[28] = new TH1F("Eta_pion_TPCPositive","Eta Distribution of +ve Pions from TPC",200,-2.0,2.0);
    fQAHist[29] = new TH1F("Phi_pion_TPCPositive","Phi Distribution of +ve Pions from TPC",1000,-0.5,15.0);
    fQAHist[30] = new TH1F("Eta_pion_TPCNegative","Eta Distribution of -ve Pions from TPC",200,-2.0,2.0);
    fQAHist[31] = new TH1F("Phi_pion_TPCNegative","Phi Distribution of -ve Pions from TPC",1000,-0.5,15.0);
    
    fQAHist[32] = new TH1F("Eta_pion_TOF","Eta Distribution of all Pions from TPC & TOF",200,-2.0,2.0);
    fQAHist[33] = new TH1F("Phi_pion_TOF","Phi Distribution of all Pions from TPC & TOF",1000,-0.5,15.0);
    fQAHist[34] = new TH1F("Eta_pion_TPCTOFPositive","Eta Distribution of +ve Pions from TPC & TOF",200,-2.0,2.0);
    fQAHist[35] = new TH1F("Phi_pion_TPCTOFPositive","Phi Distribution of +ve Pions from TPC & TOF",1000,-0.5,15.0);
    fQAHist[36] = new TH1F("Eta_pion_TPCTOFNegative","Eta Distribution of -ve Pions from TPC & TOF",200,-2.0,2.0);
    fQAHist[37] = new TH1F("Phi_pion_TPCTOFNegative","Phi Distribution of -ve Pions from TPC & TOF",1000,-0.5,15.0);
    TString Histname4;
    for (Int_t j =0; j<10; j++){
        Histname4 = "Mean_pT_TPC_";
        Histname4+= j;
        MeanpTebyeTPC[j]=new TH1F(Histname4.Data(),"EbyE #LT p_{T} #GT Distribution TPC",1000,0.0,1.0);
        MeanpTebyeTPC[j]->Sumw2();
    }
    TString Histname5;
    for (Int_t j =0; j<10; j++){
        Histname5 = "Teff_TPC_";
        Histname5+= j;
        Teff_ebyeTPC[j]=new TH1F(Histname5.Data(),"EbyE T_{eff} Distribution TPC",1000,0.05,0.55);
        Teff_ebyeTPC[j]->Sumw2();
    }
    

    TString Histname6;
    for (Int_t j =0; j<10; j++){
        Histname6 = "No_pion_evt_TPC_";
        Histname6+= j;
        No_pion_evt_TPC[j]=new TH1F(Histname6.Data(),"EbyE No of pions Distribution TPC",3000,0,3000);
         No_pion_evt_TPC[j]->Sumw2();
    }
    TString Histname7;
    for (Int_t j =0; j<10; j++){
        Histname7 = "Mean_pT_TPCTOF_";
        Histname7+= j;
        MeanpTebyeTPCTOF[j]=new TH1F(Histname7.Data(),"EbyE #LT p_{T} #GT Distribution TPC-TOF",1000,0.0,1.0);
         MeanpTebyeTPCTOF[j]->Sumw2();
    }
    TString Histname8;
    for (Int_t j =0; j<10; j++){
        Histname8 = "Teff_TPCTOF_";
        Histname8+= j;
        Teff_ebyeTPCTOF[j]=new TH1F(Histname8.Data(),"EbyE T_{eff} Distribution TPC-TOF",1000,0.05,0.55);
        Teff_ebyeTPCTOF[j]->Sumw2();
    }
    
    
    TString Histname9;
    for (Int_t j =0; j<10; j++){
        Histname9 = "No_pion_evt_TPCTOF_";
        Histname9+= j;
        No_pion_evt_TPCTOF[j]=new TH1F(Histname9.Data(),"EbyE No of pions Distribution TPC-TOF",3000,0,3000);
         No_pion_evt_TPCTOF[j]->Sumw2();
    }
    TString Histname10;
    for (Int_t j =0; j<10; j++){
        Histname10 = "Mean_pT_Truth_";
        Histname10+= j;
        MeanpTebyeTruth[j]=new TH1F(Histname10.Data(),"EbyE #LT p_{T} #GT Distribution Truth",1000,0.0,1.0);
         MeanpTebyeTruth[j]->Sumw2();
    }
    TString Histname11;
    for (Int_t j =0; j<10; j++){
        Histname11 = "Teff_Truth_";
        Histname11+= j;
        Teff_ebyeTruth[j]=new TH1F(Histname11.Data(),"EbyE T_{eff} Distribution Truth",1000,0.05,0.55);
        Teff_ebyeTruth[j]->Sumw2();
    }
    
    
    TString Histname12;
    for (Int_t j =0; j<10; j++){
        Histname12 = "No_pion_evt_Truth_";
        Histname12+= j;
        No_pion_evt_Truth[j]=new TH1F(Histname12.Data(),"EbyE No of pions Distribution Truth",3000,0,3000);
        No_pion_evt_Truth[j]->Sumw2();
    }

    
//    
    
    
    
    fOutput->Add(fHistPM);
    fOutput->Add(fHistdEdx);
    fOutput->Add(fHistdEdxSigma);
    fOutput->Add(fNSigmaTPC);
    fOutput->Add(fNSigmaTPCCut);
    fOutput->Add(fNSigmaTPCTOF);
    fOutput->Add(fNSigmaTPCTOFCut);
    fOutput->Add(fHistTOF);
    fOutput->Add(fNSigmaTOF);
    fOutput->Add(fHistPtTPC);
    fOutput->Add(fHistRecTPC);
    fOutput->Add(fHistPtTPCTOF);
    fOutput->Add(fHistRecTPCTOF);
    if ( AnalysisMode == "MC" || AnalysisMode == "mc")
    {
        fOutput->Add(fHistRecPrimaryTPC);
        fOutput->Add(fHistRecSecWDTPC);
        fOutput->Add(fHistRecSecMatTPC);
        fOutput->Add(fHistRecMisIdTPC);
        fOutput->Add(fHistRecPrimaryTPCTOF);
        fOutput->Add(fHistRecSecWDTPCTOF);
        fOutput->Add(fHistRecSecMatTPCTOF);
        fOutput->Add(fHistRecMisIdTPCTOF);
        fOutput->Add(fHistTruth);
        fOutput->Add(fHistPtTruth);


    }
    for(Int_t j=0;j<10;j++)
    {
        
        fOutput->Add(MeanpTebyeTPC[j]);
        fOutput->Add(Teff_ebyeTPC[j]);
        fOutput->Add(No_pion_evt_TPC[j]);
        fOutput->Add(MeanpTebyeTPCTOF[j]);
        fOutput->Add(Teff_ebyeTPCTOF[j]);
        fOutput->Add(No_pion_evt_TPCTOF[j]);
        if ( AnalysisMode == "MC" || AnalysisMode == "mc")
        {
        fOutput->Add(MeanpTebyeTruth[j]);
        fOutput->Add(Teff_ebyeTruth[j]);
        fOutput->Add(No_pion_evt_Truth[j]);
        }
    }
    for(Int_t iQA=0;iQA<38;iQA++)
    {
        fOutput->Add(fQAHist[iQA]);
    }
    // make sure Efficiency.root file is in the same folder
    if(fcorrection)
    {
     if (TString(fefffilename).BeginsWith("alien:"))
        TGrid::Connect("alien:");
        TFile *file=TFile::Open(fefffilename);
    fEffptTPC = (TH1F *)file->Get("ptTPC");
    fEffptTPCTOF = (TH1F *)file->Get("ptTPCTOF");
    file->Close();
    }
    
    

    // NEW HISTO added to fOutput here
    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________

void AliAnalysisTempFluc::UserExec(Option_t *) 
{
    
    
    AliVEvent *event = InputEvent();
    if (!event) { Printf("ERROR: Could not retrieve event"); return; }
    
    if (DataType == "AOD" || DataType == "aod")
        
    {
     fAOD = dynamic_cast<AliAODEvent*>(event);
        if (!fAOD) {
             AliError("ERROR: fAOD not available\n");
            return;
        }
    }
      if(!fPIDResponse) return;
    
    if ( AnalysisMode == "MC" || AnalysisMode == "mc")
    {
        //retreive MC particles from event
        mcTruth = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!mcTruth)
        {
            AliInfo("No MC particle branch found");
            return;
        }
    }
    
    fQAHist[4]->Fill(1); //precessed event histo

  if (!AcceptEvent(fAOD)) return;
    
     fQAHist[4]->Fill(3); // Accepeted event histo
    
    AliAODVertex *fVertex = fAOD->GetPrimaryVertex(); //primary vertex
    if (!fVertex) return;
    

  AliCentrality * fCentrality = fAOD->GetCentrality();
  Float_t Centrality = fCentrality->GetCentralityPercentile(fCentralityEstimator);
    
   if (Centrality < fCentralityCutL || Centrality > fCentralityCutH) return;

    fQAHist[4]->Fill(7); // Centrality cut histo
if(mcTruth)
{
     CalEfficiencyMap();
}
  fNEvt++;
  Int_t ntracks = fAOD->GetNumberOfTracks();
  Int_t neg = 0, pos = 0;
  
  for(Int_t i = 0; i < ntracks; i++) {
      AliAODTrack *track =(AliAODTrack *)fAOD->GetTrack(i); // pointer to reconstructed to track
    if(!track) { 
      AliError(Form("ERROR: Could not retrieve aod track %d\n",i));
      continue; 
    }
      
      if(!track->TestFilterBit(fBit)) continue;
      fQAHist[8]->Fill(track->Eta());
      fQAHist[9]->Fill(track->Phi());
      Float_t charge = track->Charge();
      if(charge > 0)
      {
          fQAHist[10]->Fill(track->Eta());
          fQAHist[11]->Fill(track->Phi());
      }
      if(charge < 0)
      {
          fQAHist[12]->Fill(track->Eta());
          fQAHist[13]->Fill(track->Phi());
      }
      if (track->Pt() < fPtL || track->Pt()> fPtH) continue;
      fQAHist[14]->Fill(track->Eta());
      fQAHist[15]->Fill(track->Phi());
      if(charge > 0)
      {
          fQAHist[16]->Fill(track->Eta());
          fQAHist[17]->Fill(track->Phi());
      }
      if(charge < 0)
      {
          fQAHist[18]->Fill(track->Eta());
          fQAHist[19]->Fill(track->Phi());
      }

      
      if (track->Eta() < gEtaL || track->Eta()> gEtaH) continue;
      fQAHist[20]->Fill(track->Eta());
      fQAHist[21]->Fill(track->Phi());
      if(charge > 0)
      {
          fQAHist[22]->Fill(track->Eta());
          fQAHist[23]->Fill(track->Phi());
      }
      if(charge < 0)
      {
          fQAHist[24]->Fill(track->Eta());
          fQAHist[25]->Fill(track->Phi());
      }

      
      
      DCAxy=track->DCA();
      DCAz=track->ZAtDCA();
    
      fQAHist[0]->Fill(DCAxy);
      fQAHist[1]->Fill(DCAz);
    
      
      Double_t nclus = track->GetTPCClusterInfo(2,1);
      Double_t chi2ndf = track->Chi2perNDF();
      
      fQAHist[2]->Fill(nclus);
      fQAHist[3]->Fill(chi2ndf);
      
      
      
      Bool_t TPCPIDSwitch = CheckTPC(track);
      if (!TPCPIDSwitch) continue;
      
          Double_t Ptpc = track->Pt();
          Double_t dEdx = track->GetTPCsignal();
      
          Float_t rigidity = Ptpc*charge;
          fHistdEdx->Fill(rigidity,dEdx);
      
    Double_t  PionNSigmaTPC = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion);
          fNSigmaTPC->Fill(Ptpc,PionNSigmaTPC);
    Double_t rapTPC = CalculateRapidity(track,mass_pion);
     //   -------------------------    TPC NSigma   -----------------------------
      if (TMath::Abs(PionNSigmaTPC)< fCutNSigmaTPC)
           {
               
               if (rapTPC < gRapH && rapTPC> gRapL)
               {
                   fQAHist[26]->Fill(track->Eta());
                   fQAHist[27]->Fill(track->Phi());
                   Float_t charge = track->Charge();
                   if(charge > 0)
                   {
                       fQAHist[28]->Fill(track->Eta());
                       fQAHist[29]->Fill(track->Phi());
                   }
                   if(charge < 0)
                   {
                       fQAHist[30]->Fill(track->Eta());
                       fQAHist[31]->Fill(track->Phi());
                   }

                   
               
               fHistdEdxSigma->Fill(rigidity,dEdx);
               fNSigmaTPCCut->Fill(Ptpc,PionNSigmaTPC);
               fHistPtTPC->Fill(Ptpc);
               fHistRecTPC->Fill(rapTPC,track->Phi(),track->Pt());
               
               
                   if (mcTruth)
                   {
                       Int_t label = TMath::Abs(track->GetLabel());
                       AliAODMCParticle * trackTruth = (AliAODMCParticle*) mcTruth->At(label);
                       if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())==211 && trackTruth->Y() < gRapH && trackTruth->Y()> gRapL){ fHistPtTruthdummy->Fill(trackTruth->Pt());}

                       if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())==211)
                       {fHistRecPrimaryTPC->Fill(rapTPC,track->Phi(),track->Pt());}
                       if(trackTruth->IsSecondaryFromWeakDecay() && TMath::Abs(trackTruth->PdgCode())==211)
                       {fHistRecSecWDTPC->Fill(rapTPC,track->Phi(),track->Pt());}
                       if(trackTruth->IsSecondaryFromMaterial() && TMath::Abs(trackTruth->PdgCode())==211)
                       {fHistRecSecMatTPC->Fill(rapTPC,track->Phi(),track->Pt());}
                       if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())!=211)
                       {fHistRecMisIdTPC ->Fill(rapTPC,track->Phi(),track->Pt());}
                   }

               }
           }
      

    if(charge < 0) neg++;
    if(charge > 0) pos++;
      
      
      // -------------------------    TPC + TOF NSigma   -----------------------------
      if(track->Pt()<0.6)
      {
          
          fNSigmaTPCTOF->Fill(Ptpc,PionNSigmaTPC);
          
          
          if (TMath::Abs(PionNSigmaTPC)<fCutNSigmaTPC){     // 2.0, 3.0 ,4.0
              
              
                  fNSigmaTPCTOFCut->Fill(Ptpc,PionNSigmaTPC);
                  fHistPtTPCTOF->Fill(track->Pt());
                  fHistRecTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());
              
              if (mcTruth)
              {
                  Int_t label = TMath::Abs(track->GetLabel());
                  AliAODMCParticle * trackTruth = (AliAODMCParticle*) mcTruth->At(label);
                  if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecPrimaryTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsSecondaryFromWeakDecay() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecSecWDTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsSecondaryFromMaterial() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecSecMatTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())!=211)
                  {fHistRecMisIdTPCTOF ->Fill(rapTPC,track->Phi(),track->Pt());}
              }


          }
      }
      else
          
      {
          Bool_t TOFPIDSwitch = CheckTOF(track);
          if (!TOFPIDSwitch) continue;
          Double_t P = track->P();
          double_t beta = Getbeta(track);
          fHistTOF->Fill(P,beta);
          
    Double_t PionNSigmaTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion);
          fNSigmaTOF->Fill(Ptpc,PionNSigmaTOF);
          
    Double_t PionNSigmaTPCTOF = TMath::Sqrt(TMath::Power((PionNSigmaTPC/fCutNSigmaTPC),2.0)+ TMath::Power((PionNSigmaTOF/fCutNSigmaTPCTOF),2.0));
      fNSigmaTPCTOF->Fill(Ptpc,PionNSigmaTPCTOF);
          
     if (PionNSigmaTPCTOF<1.0)
          {
              fQAHist[32]->Fill(track->Eta());
              fQAHist[33]->Fill(track->Phi());
              Float_t charge = track->Charge();
              if(charge > 0)
              {
                  fQAHist[34]->Fill(track->Eta());
                  fQAHist[35]->Fill(track->Phi());
              }
              if(charge < 0)
              {
                  fQAHist[36]->Fill(track->Eta());
                  fQAHist[37]->Fill(track->Phi());
              }
              

              
                  fNSigmaTPCTOFCut->Fill(Ptpc,PionNSigmaTPCTOF);
                  fHistPtTPCTOF->Fill(track->Pt());
                  fHistRecTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());
            
              if (mcTruth)
              {
                  Int_t label = TMath::Abs(track->GetLabel());
                  AliAODMCParticle * trackTruth = (AliAODMCParticle*) mcTruth->At(label);
                  if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecPrimaryTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsSecondaryFromWeakDecay() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecSecWDTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsSecondaryFromMaterial() && TMath::Abs(trackTruth->PdgCode())==211)
                  {fHistRecSecMatTPCTOF->Fill(rapTPC,track->Phi(),track->Pt());}
                  if(trackTruth->IsPhysicalPrimary() && TMath::Abs(trackTruth->PdgCode())!=211)
                  {fHistRecMisIdTPCTOF ->Fill(rapTPC,track->Phi(),track->Pt());}
              }
              

          }
              
 }
//-----------------end

    
}//track loop
  
 
  fHistPM->Fill(pos,neg);
    if(fcorrection)
    {
        
        fHistPtTPC->Divide(fEffptTPC);
        fHistPtTPCTOF->Divide(fEffptTPCTOF);
    }

    CalMeanPtEbyETPC(fHistPtTPC);
    CalMeanPtEbyETPCTOF(fHistPtTPCTOF);
    if(mcTruth){CalMeanPtEbyETruth(fHistPtTruthdummy);}

  PostData(1, fOutput);
  
}//event loop


//________________________________________________________________________
void AliAnalysisTempFluc::Terminate(Option_t *) 
{
   
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
       
}

Bool_t
AliAnalysisTempFluc::AcceptEvent(AliAODEvent *event) const
{
    Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB);
    if(!isSelected) return kFALSE;

    
  Bool_t ver = kFALSE;
  const AliAODVertex *vertex = event->GetPrimaryVertex();
  if(vertex) {
    Double32_t fCov[6];
    vertex->GetCovarianceMatrix(fCov);
    if(vertex->GetNContributors() > 0) {
      if(fCov[5] != 0) {
	
	if(TMath::Abs(vertex->GetX()) < 3.0) {
	  if(TMath::Abs(vertex->GetY()) < 3.0) {
	    if(TMath::Abs(vertex->GetZ()) < 10.0) {
	      ver = kTRUE;
            fQAHist[5]->Fill(vertex->GetX());
            fQAHist[6]->Fill(vertex->GetY());
            fQAHist[7]->Fill(vertex->GetZ());
	    }
	  }
	}
      }
    }
  }
  
  AliCentrality *centrality = event->GetCentrality();
  if (centrality->GetQuality() != 0) ver = kFALSE;
  return ver;
}
//_____________________________________
void AliAnalysisTempFluc::CalMeanPtEbyETPC(TH1* ebye_pT_temp)
{
    
    Float_t gPtMin=-1.0;
    Float_t gPtMax=-1.0;
    TString config1(fBinningptL);
    TObjArray* binningptL = config1.Tokenize(",");
    Double_t* PtRangeArrayLow = new Double_t[binningptL->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptL->GetEntriesFast(); j++)
    {PtRangeArrayLow[j]=TString(binningptL->At(j)->GetName()).Atof();}
    
    TString config2(fBinningptH);
    TObjArray* binningptH = config2.Tokenize(",");
    Double_t* PtRangeArrayHigh = new Double_t[binningptH->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptH->GetEntriesFast(); j++)
    {PtRangeArrayHigh[j]=TString(binningptH->At(j)->GetName()).Atof();}
    
      Int_t nPtbins=binningptL->GetEntriesFast();
    Int_t binlow,binhigh;
    Double_t MeanPt,Teff_MeanPt;
    Int_t No_pion;
    
    for(Int_t ibin=0;ibin< nPtbins;ibin++)
    {
        MeanPt=0.0;
        Teff_MeanPt=0.0;
        No_pion=0;
        gPtMin=PtRangeArrayLow[ibin];
        gPtMax=PtRangeArrayHigh[ibin];
        
        binlow=ebye_pT_temp->GetXaxis()->FindBin(gPtMin);
        binhigh=ebye_pT_temp->GetXaxis()->FindBin(gPtMax);
        ebye_pT_temp->GetXaxis()->SetRange(binlow,binhigh);
        MeanPt=ebye_pT_temp->GetMean();
        No_pion=ebye_pT_temp->Integral(binlow,binhigh);
        if(MeanPt<0.22) continue;
        ebye_pT_temp->GetXaxis()->SetRange(0,-1);
        TSpline *spline=func_slope_expo(gPtMin,gPtMax);
        Teff_MeanPt=spline->Eval(MeanPt);
      
        MeanpTebyeTPC[ibin]->Fill(MeanPt);
        Teff_ebyeTPC[ibin]->Fill(Teff_MeanPt);
        No_pion_evt_TPC[ibin]->Fill(No_pion);
        
        
    }
  
    ebye_pT_temp->Reset();
    
}


//_____________________________________

void AliAnalysisTempFluc::CalMeanPtEbyETPCTOF(TH1* ebye_pT_temp)
{
    
    Float_t gPtMin=-1.0;
    Float_t gPtMax=-1.0;
    TString config1(fBinningptL);
    TObjArray* binningptL = config1.Tokenize(",");
    Double_t* PtRangeArrayLow = new Double_t[binningptL->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptL->GetEntriesFast(); j++)
    {PtRangeArrayLow[j]=TString(binningptL->At(j)->GetName()).Atof();}
    
    TString config2(fBinningptH);
    TObjArray* binningptH = config2.Tokenize(",");
    Double_t* PtRangeArrayHigh = new Double_t[binningptH->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptH->GetEntriesFast(); j++)
    {PtRangeArrayHigh[j]=TString(binningptH->At(j)->GetName()).Atof();}
    
    Int_t nPtbins=binningptL->GetEntriesFast();
    Int_t binlow,binhigh;
    Double_t MeanPt,Teff_MeanPt;
    Int_t No_pion;
    
    for(Int_t ibin=0;ibin< nPtbins;ibin++)
    {
        MeanPt=0.0;
        Teff_MeanPt=0.0;
        No_pion=0;
        gPtMin=PtRangeArrayLow[ibin];
        gPtMax=PtRangeArrayHigh[ibin];
        
        binlow=ebye_pT_temp->GetXaxis()->FindBin(gPtMin);
        binhigh=ebye_pT_temp->GetXaxis()->FindBin(gPtMax);
        ebye_pT_temp->GetXaxis()->SetRange(binlow,binhigh);
        MeanPt=ebye_pT_temp->GetMean();
        No_pion=ebye_pT_temp->Integral(binlow,binhigh);
        if(MeanPt<0.22) continue;
        ebye_pT_temp->GetXaxis()->SetRange(0,-1);
        TSpline *spline=func_slope_expo(gPtMin,gPtMax);
        Teff_MeanPt=spline->Eval(MeanPt);
        MeanpTebyeTPCTOF[ibin]->Fill(MeanPt);
        Teff_ebyeTPCTOF[ibin]->Fill(Teff_MeanPt);
        No_pion_evt_TPCTOF[ibin]->Fill(No_pion);
        
        
    }
  
    ebye_pT_temp->Reset();
 }
//----------------------------------------------------------------
void AliAnalysisTempFluc::CalMeanPtEbyETruth(TH1* ebye_pT_temp)
{
    
    Float_t gPtMin=-1.0;
    Float_t gPtMax=-1.0;
    TString config1(fBinningptL);
    TObjArray* binningptL = config1.Tokenize(",");
    Double_t* PtRangeArrayLow = new Double_t[binningptL->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptL->GetEntriesFast(); j++)
    {PtRangeArrayLow[j]=TString(binningptL->At(j)->GetName()).Atof();}
    
    TString config2(fBinningptH);
    TObjArray* binningptH = config2.Tokenize(",");
    Double_t* PtRangeArrayHigh = new Double_t[binningptH->GetEntriesFast()];
    
    for (Int_t j=0; j<binningptH->GetEntriesFast(); j++)
    {PtRangeArrayHigh[j]=TString(binningptH->At(j)->GetName()).Atof();}
     Int_t nPtbins=binningptL->GetEntriesFast();
    Int_t binlow,binhigh;
    Double_t MeanPt,Teff_MeanPt;
    Int_t No_pion;
    
    for(Int_t ibin=0;ibin< nPtbins;ibin++)
    {
        MeanPt=0.0;
        Teff_MeanPt=0.0;
        No_pion=0;
        gPtMin=PtRangeArrayLow[ibin];
        gPtMax=PtRangeArrayHigh[ibin];
        
        binlow=ebye_pT_temp->GetXaxis()->FindBin(gPtMin);
        binhigh=ebye_pT_temp->GetXaxis()->FindBin(gPtMax);
        ebye_pT_temp->GetXaxis()->SetRange(binlow,binhigh);
        MeanPt=ebye_pT_temp->GetMean();
        No_pion=ebye_pT_temp->Integral(binlow,binhigh);
        
        if(MeanPt<0.22) continue;
        ebye_pT_temp->GetXaxis()->SetRange(0,-1);
        TSpline *spline=func_slope_expo(gPtMin,gPtMax);
        Teff_MeanPt=spline->Eval(MeanPt);
         MeanpTebyeTruth[ibin]->Fill(MeanPt);
        Teff_ebyeTruth[ibin]->Fill(Teff_MeanPt);
        No_pion_evt_Truth[ibin]->Fill(No_pion);
        
        
    }
   
    ebye_pT_temp->Reset();
  
}
//----------------------------------------------------------------

TSpline *AliAnalysisTempFluc::func_slope_expo(const Float_t fitmin, const Float_t fitmax)
{
    Double_t x[1000]={0};
    Double_t y[1000]={0};
    
    Float_t pT,T;
    
    for(Int_t i=0; i<1000; i++)
    {
        T=0.0003*i+0.1;
        pT = ((TMath::Power(fitmin,2)*TMath::Exp(-fitmin/T)-TMath::Power(fitmax,2)*TMath::Exp(-fitmax/T))/((fitmin+T)*TMath::Exp(-fitmin/T)-(fitmax+T)*TMath::Exp(-fitmax/T)))+2*T;
        x[i]=pT;
        y[i]=T;
        
    }
    
    TGraph *gr = new TGraph(1000,x,y);
    TSpline3 *sp= new TSpline3("sp",x,y,1000);
    return sp;
    
}

// ============== For TPC PID CHECK ==================

Bool_t AliAnalysisTempFluc::CheckTPC(AliAODTrack *track)
{
    if ((track->GetStatus() & AliAODTrack::kTPCin   ) == 0) return kFALSE;
    if ((track->GetStatus() & AliAODTrack::kTPCrefit) == 0) return kFALSE;
    if ((track->GetStatus() & AliAODTrack::kITSrefit) == 0) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTempFluc::CheckTOF(AliAODTrack * track)
{
    Double_t TOFPtRange =0.6;         // Controling Agent 1
    
    // *************** Check if the particle has TOF Matching ***************
    
    UInt_t status;
    status=track->GetStatus();
    if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)
        return kFALSE;
    // TPC TOF mismatch is be implemented
    Float_t length = track->GetIntegratedLength();
    if (length > 350.)         // in this AOD production it returns 0 ... so not realy implemented for implement replace ">" by "<"
        // menas the min length is reached to TOF
        return kTRUE;
    
    // -------  in addition to KTOFout and kTIME we look at the pt  ------
    
    if(track->Pt()<TOFPtRange) return kFALSE;
    return kTRUE;
}

Double_t AliAnalysisTempFluc::Getbeta(AliAODTrack *track)
{
    Double_t P = track->P();
    Double_t t0 = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(P);
    Double_t tf[10];
    track->GetIntegratedTimes(tf);
    return tf[0]/t0;
    
}
Double_t AliAnalysisTempFluc::CalculateRapidity(AliAODTrack *track , Double_t mass)
{
    Double_t E,rap,pz,pt;
    pt=track->Pt();
    pz = track->Pz();
    E = TMath::Sqrt(pt*pt+pz*pz+mass*mass);
    rap = 0.5 * TMath::Log ((E+pz)/(E-pz));
    return rap;
}
void  AliAnalysisTempFluc::CalEfficiencyMap()
{
    
    Int_t ntracks = mcTruth->GetEntriesFast();                      //number of tracks in event
    
    for (Int_t j = 0; j <ntracks; j++) {
        AliAODMCParticle *trackTruth = (AliAODMCParticle*)mcTruth->At(j);
        if (!trackTruth) continue;                                       //if the position "j" is not empty, then continue
        if(!TMath::Abs(trackTruth->Charge()) || !trackTruth->IsPhysicalPrimary() || trackTruth->Pt()== 0) continue;
        //
        // fill MC histograms
        //
        
        if(TMath::Abs(trackTruth->PdgCode())==211) {
            if (trackTruth->Pt() < fPtL || trackTruth->Pt()> fPtH) continue;
            if (trackTruth->Y() < gRapH && trackTruth->Y()> gRapL){
                fHistTruth->Fill(trackTruth->Y(), trackTruth->Phi(), trackTruth->Pt()); // Fill identical hist as for real data
                 fHistPtTruth->Fill(trackTruth->Pt());
                
            }
        }
    }//track loop mc
    CalMeanPtEbyETPC(fHistPtTruth);
}

