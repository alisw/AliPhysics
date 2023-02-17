#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TList.h"
#include "TObjArray.h"
#include "TString.h"
#include "TParticle.h"
#include "TProfile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"

#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliInputEventHandler.h"

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliPIDResponse.h"
//#include "AliAnalysisUtils.h"

#include "AliAODpidUtil.h"
#include "AliAODHeader.h"

#include "AliAnalysisTaskParticlePurityEfficiency_pik.h"




ClassImp(AliAnalysisTaskParticlePurityEfficiency_pik)
//ClassImp(AliAnalysisTaskParticleEff)

double fV1[3];

//_______________________________________________________

AliAnalysisTaskParticlePurityEfficiency_pik::AliAnalysisTaskParticlePurityEfficiency_pik(const Char_t *partName) :
  AliAnalysisTaskSE(partName), centrality(0), fHistoList(0),  fHistEv(0), fpidResponse(0), fAODpidUtil(0)
{

  for(Int_t i = 0; i < CENTRBINS*PARTTYPES; i++) {
    for(Int_t chg = 0; chg < 2; chg++) {
      fGeneratedMCPrimaries[i][chg] = NULL;
      fMCPrimariesThatAreReconstructed[i][chg] = NULL;
      fMCPrimariesThatAreReconstructedNoNsigma[i][chg] = NULL;
      fReconstructedAfterCuts[i][chg] = NULL;
      fReconstructedNotPrimaries[i][chg] = NULL;
      fReconstructedPrimaries[i][chg] = NULL;
      fContamination[i][chg] = NULL;
    }
  }
  for ( Int_t i = 0; i < 11; i++) { 
    fHistQA[i] = NULL;
    fHistQA2D[i] = NULL;
  }
  for (Int_t i = 0; i < 4; i++) {
    fProfilePTrueReconstructed[i] = NULL;
    fProfileThetaTrueReconstructed[i] = NULL;
    fProfilePhiTrueReconstructed[i] = NULL;
  }
  DefineOutput(1, TList::Class());
}

//_______________________________________________________

AliAnalysisTaskParticlePurityEfficiency_pik::~AliAnalysisTaskParticlePurityEfficiency_pik()
{
  // Destructor
  if (AliAnalysisManager::GetAnalysisManager()->GetAnalysisType() != AliAnalysisManager::kProofAnalysis)
    delete fHistoList;
}

//_______________________________________________________

void AliAnalysisTaskParticlePurityEfficiency_pik::UserCreateOutputObjects()
{

  /* create output */
  fHistoList = new TList();
  fHistoList->SetOwner(kTRUE);
   
  TString hname1, hname2, hname3, hname4, hname5, hname6, hname7, hname8, hname9, hname10, hname11,hname12,hname13,hname14,hname15,hname16,hname17,hname18,hname19,hname20,hname21,hname22,hname23,hname24,hname25,hname26 ;
  
  TString htitle1, htitle2, htitle3, htitle4,htitle6,htitle5,htitle7,htitle8,htitle9,htitle10,htitle11,htitle12,htitle13,htitle14,htitle15,htitle16,htitle17,htitle18,htitle19,htitle20,htitle21,htitle22,htitle23,htitle24,htitle25,htitle26;
  
  TString hname1M, hname2M, hname3M, hname4M, hname5M, hname;
  
  TString htitle1M, htitle2M, htitle3M, htitle4M, htitle5M, htitle;

  TString parttypename = "None";

  for(Int_t j = 0; j < PARTTYPES; j++)  {
    if (j==0) parttypename="All";
    else if (j==1) parttypename="Pion";
    else if (j==2) parttypename="Kaon";
    else if (j==3) parttypename="Proton";

    for(Int_t i = 0; i < CENTRBINS; i++)  {
      hname1  = "hGeneratedMCPrimariesEffM"; hname1+=i; hname1+=parttypename;
      htitle1 = "Kinematic level eta_pT (prim only) M"; htitle1+=i; htitle1+=parttypename;
      fGeneratedMCPrimaries[i*PARTTYPES+j][0] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname1+="Minus";htitle1+="Minus";
      fGeneratedMCPrimaries[i*PARTTYPES+j][1] = new TH2F(hname1.Data(),htitle1.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname3  = "hMCPrimariesThatAreReconstructedNoNsigmaM"; hname3+=i; hname3+=parttypename;
      htitle3 = "Reconstructed level eta_pT (prim only) no Nsigma cut only PDG M"; htitle3+=i; htitle3+=parttypename;
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname3+="Minus";htitle3+="Minus";
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][1] = new TH2F(hname3.Data(),htitle3.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname2  = "hHistoReconstructedAfterCutsM"; hname2+=i; hname2+=parttypename;
      htitle2 = "Total Reconstructed tracks M "; htitle2+=i; htitle2+=parttypename;
      fReconstructedAfterCuts[i*PARTTYPES+j][0] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname2+="Minus";htitle2+="Minus";
      fReconstructedAfterCuts[i*PARTTYPES+j][1] = new TH2F(hname2.Data(),htitle2.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedNotPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (not primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedNotPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedNotPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname4  = "hHistoReconstructedPrimariesM"; hname4+=i; hname4+=parttypename;
      htitle4 = "Reconstructed level eta_pT (primaries) M"; htitle4+=i; htitle4+=parttypename;
      fReconstructedPrimaries[i*PARTTYPES+j][0] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);
      hname4+="Minus";htitle4+="Minus";
      fReconstructedPrimaries[i*PARTTYPES+j][1] = new TH2F(hname4.Data(),htitle4.Data(),50, -1.5, 1.5,1000,0.,10.0);

      hname5  = "hContaminationM"; hname5+=i; hname5+=parttypename;
      htitle5 = "Contamination M"; htitle5+=i; htitle5+=parttypename;
      fContamination[i*PARTTYPES+j][0] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,50,0.,10.0); //50
      hname5+="Minus";htitle5+="Minus";
      fContamination[i*PARTTYPES+j][1] = new TH2F(hname5.Data(),htitle5.Data(),6000, -3000, 3000.,50,0.,10.0); //50

      fReconstructedAfterCuts[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedNotPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fReconstructedPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructedNoNsigma[i*PARTTYPES+j][0]->Sumw2();
      fMCPrimariesThatAreReconstructed[i*PARTTYPES+j][0]->Sumw2();
      fGeneratedMCPrimaries[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][0]->Sumw2();
      fContamination[i*PARTTYPES+j][1]->Sumw2();
    }
    
    hname  = "pidTPCdEdx";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum";
    fHistQAPID[0][j][0] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTime";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum";
    fHistQAPID[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum";
    fHistQAPID[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[2][j][1]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum";
    fHistQAPID[3][j][0] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigma";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma";
    fHistQAPID[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPID[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);

    hname  = "pidTPCdEdxFail";  hname+=parttypename;
    htitle = parttypename + " TPC dEdx vs. momentum Fail";
    fHistQAPIDFail[0][j][0] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[0][j][1] = new TH2F(hname, htitle, 100, 0.0, 5.0, 250, 0.0, 500.0);
    hname  = "pidTOFTimeFail";  hname+=parttypename;
    htitle = parttypename + " TOF Time vs. momentum Fail";
    fHistQAPIDFail[1][j][0] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[1][j][1] = new TH2F(hname, htitle, 100, 0.1, 5.0, 400, -4000.0, 4000.0);
    hname  = "pidTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TOF NSigma vs. momentum Fail";
    fHistQAPIDFail[2][j][0]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[2][j][1]= new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC NSigma vs. momentum Fail";
    fHistQAPIDFail[3][j][0] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[3][j][1] = new TH2F(hname,htitle, 100, 0.0, 5.0, 100, -5.0, 5.0);
    hname  = "pidTPCTOFNSigmaFail";  hname+=parttypename;
    htitle = parttypename + " TPC vs TOF NSigma Fail";
    fHistQAPIDFail[4][j][0] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
    htitle+="Minus"; hname+="Minus";
    fHistQAPIDFail[4][j][1] = new TH2F(hname,htitle, 200, -10.0, 10.0, 200, -10.0, 10.0);
  }

  fHistEv = new TH1F("fHistEv", "Multiplicity", 100, 0, 100);
  fHistoList->Add(fHistEv);

  for(Int_t i = 0; i < CENTRBINS; i++)  {
    hname = "fHistEventCutsM";
    hname+= i;
    
    fHistEvCuts[i] = new TH1F(hname,Form("Event Cuts M%d",i) , 4, 0, 5);
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(1,"All");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(2,"NoVertex");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(3,"PileUp");
    fHistEvCuts[i]->GetXaxis()->SetBinLabel(4,"z-vertex>10");
    fHistoList->Add(fHistEvCuts[i]);

    for(Int_t chg=0;chg<1;chg++){
      hname  = "hMisidentificationM"; hname+=i; 
      htitle = "Misidentification Fraction M"; htitle+=i; 
      fMisidentification[i][0] = new TH3F(hname.Data(),htitle.Data(), 3, 0.5, 3.5, 4 , 0, 4, 200,0.1,2.5);
      fMisidentification[i][0]->GetXaxis()->SetBinLabel(1,"Pions, MC");
      fMisidentification[i][0]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
      fMisidentification[i][0]->GetXaxis()->SetBinLabel(3,"Protons, MC");
      fMisidentification[i][0]->GetYaxis()->SetBinLabel(1,"Pions, Data");
      fMisidentification[i][0]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
      fMisidentification[i][0]->GetYaxis()->SetBinLabel(3,"Protons, Data");
      fMisidentification[i][0]->GetYaxis()->SetBinLabel(4,"Other, Data");
      fHistoList->Add(fMisidentification[i][0]);

    htitle+="Minus"; hname+="Minus";
      fMisidentification[i][1] = new TH3F(hname.Data(),htitle.Data(), 3, 0.5, 3.5, 4 , 0, 4, 200,0.1,2.5);
      fMisidentification[i][1]->GetXaxis()->SetBinLabel(1,"Pions, MC");
      fMisidentification[i][1]->GetXaxis()->SetBinLabel(2,"Kaons, MC");
      fMisidentification[i][1]->GetXaxis()->SetBinLabel(3,"Protons, MC");
      fMisidentification[i][1]->GetYaxis()->SetBinLabel(1,"Pions, Data");
      fMisidentification[i][1]->GetYaxis()->SetBinLabel(2,"Kaons, Data");
      fMisidentification[i][1]->GetYaxis()->SetBinLabel(3,"Protons, Data");
      fMisidentification[i][1]->GetYaxis()->SetBinLabel(4,"Other, Data");
      fHistoList->Add(fMisidentification[i][1]);

    }
  }

  
  fHistQA[0] = new TH1F("fHistVtx", "Z vertex distribution", 100, -15., 15.);
  fHistQA[1] = new TH1F("fHistnTpcCluster", "n TPC Cluster", 100, 0., 200.);
  fHistQA[2] = new TH1F("fHistnTpcClusterF", "n TPC Cluster findable", 100, 0., 200.);
  fHistQA[3] = new TH1F("dcaHistDcaXY1D", "DCA XY", 210, -2.1, 2.1);
  fHistQA[4] = new TH1F("dcaHistDcaZ1D", "DCA Z", 210, -2.1, 2.1);
  fHistQA[5] = new TH1F("fHistChi2Tpc", "Chi2 TPC", 100, 0., 8.);
  fHistQA[6] = new TH1F("fHistpT", "pT distribution",200,0.1,2.5);
  fHistQA[7] = new TH1F("fHistPhi", "Phi distribution" , 100, -TMath::Pi(), TMath::Pi());
  fHistQA[8] = new TH1F("fHistEta", "Eta distribution" , 100, -2, 2);
 
  fHistQA[9] = new TH1F("fHistEventCuts", "Event Cuts" , 4, 0, 5);
  fHistQA[9]->GetXaxis()->SetBinLabel(1,"All");
  fHistQA[9]->GetXaxis()->SetBinLabel(2,"NoVertex");
  fHistQA[9]->GetXaxis()->SetBinLabel(3,"PileUp");
  fHistQA[9]->GetXaxis()->SetBinLabel(4,"z-vertex>10");


  fHistQA[10] = new TH1F("fHistTrackCuts", "Track Cuts" , 7, 0.5, 7.5);
  fHistQA[10]->GetXaxis()->SetBinLabel(1,"AllTracksInEvents");
  fHistQA[10]->GetXaxis()->SetBinLabel(2,"GetTrack");
  fHistQA[10]->GetXaxis()->SetBinLabel(3,"Filter bit");
  fHistQA[10]->GetXaxis()->SetBinLabel(4,"Eta");
  fHistQA[10]->GetXaxis()->SetBinLabel(5,"Pt");
  fHistQA[10]->GetXaxis()->SetBinLabel(6,"DCA");
  fHistQA[10]->GetXaxis()->SetBinLabel(7,"Electron Rejection");

  fHistQA2D[0] = new TH2F("dcaHistDcaXY","DCA XY",50, 0, 5,210, -2.1, 2.1);
  fHistQA2D[1] = new TH2F("dcaHistDcaZ","DCA Z", 50, 0, 5, 210, -2.1, 2.1);
  fHistQA2D[2] = new TH2F("fPhiEta","Eta-Phi",100, -2, 2, 100, -TMath::Pi(), TMath::Pi());

  fHistQA2D[3] = new TH2F("dcaHistDcaXY1DPion", "DCA XY Pion", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[4] = new TH2F("dcaHistDcaZ1DPion", "DCA Z Pion", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[5] = new TH2F("dcaHistDcaXY1DKaon", "DCA XY Kaon", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[6] = new TH2F("dcaHistDcaZ1DKaon", "DCA Z Kaon", 5000, -0.5, 0.5,50, 0, 5);

  fHistQA2D[7] = new TH2F("dcaHistDcaXY1DPionMinus", "DCA XY PionMinus", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[8] = new TH2F("dcaHistDcaZ1DPionMinus", "DCA Z PionMinus", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[9] = new TH2F("dcaHistDcaXY1DKaonMinus", "DCA XY KaonMinus", 500, -0.5, 0.5,50, 0, 5);
  fHistQA2D[10] = new TH2F("dcaHistDcaZ1DKaonMinus", "DCA Z KaonMinus", 5000, -0.5, 0.5,50, 0, 5);

  
  char* particletype = "None";
  for (Int_t i = 0; i < 4; i++) {
    if (i==0) particletype="All";
    else if (i==1) particletype="Pion";
    else if (i==2) particletype="Kaon";
    else if (i==3) particletype="Proton";
    fProfilePTrueReconstructed[i] = new TProfile(Form("fPTrueReconstructed%s", particletype), "Momentum true - reconstructed", 100, 0.0, 3.0);
    fProfileThetaTrueReconstructed[i] = new TProfile(Form("fThetaTrueReconstructed%s", particletype), "Theta angle true - reconstructed", 100, 0.0, 3.0);
    fProfilePhiTrueReconstructed[i] = new TProfile(Form("fPhiTrueReconstructed%s", particletype), "Phi angle true - reconstructed", 100, 0.0, 3.0);
    fHist2DPTrueReconstructed[i] = new TH2D(Form("fHist2DPTrueReconstructed%s", particletype), "Hist Momentum true - reconstructed", 100, 0.0, 3.0, 100, -0.02, 0.02);
    fHist2DThetaTrueReconstructed[i] = new TH2D(Form("fHist2DThetaTrueReconstructed%s", particletype), "Hist Theta true - reconstructed", 100, 0.0, 3.0, 100, -0.02, 0.02);
    fHist2DPhiTrueReconstructed[i] = new TH2D(Form("fHist2DPhiTrueReconstructed%s", particletype), "Hist Phi true - reconstructed", 100, 0.0, 3.0, 100, -0.02, 0.02);
    fHistoList->Add(fProfilePTrueReconstructed[i]);
    fHistoList->Add(fProfileThetaTrueReconstructed[i]);
    fHistoList->Add(fProfilePhiTrueReconstructed[i]);
    fHistoList->Add(fHist2DPTrueReconstructed[i]);
    fHistoList->Add(fHist2DThetaTrueReconstructed[i]);
    fHistoList->Add(fHist2DPhiTrueReconstructed[i]);
  }

  
  for ( Int_t i = 0; i < 11; i++)
    {
      fHistoList->Add(fHistQA[i]);
      fHistoList->Add(fHistQA2D[i]);
      if(i<5) {
	for(Int_t j = 0 ; j<PARTTYPES; j++)
	  for(int chg=0;chg<2;chg++)
	    {
		fHistoList->Add(fHistQAPID[i][j][chg]);
		fHistoList->Add(fHistQAPIDFail[i][j][chg]);
	    }
      }
    }
    
  for (Int_t i = 0; i < CENTRBINS*PARTTYPES; i++){
    for(Int_t chg=0;chg<2;chg++){
      fHistoList->Add(fGeneratedMCPrimaries[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructed[i][chg]);
      fHistoList->Add(fMCPrimariesThatAreReconstructedNoNsigma[i][chg]);
      fHistoList->Add(fReconstructedAfterCuts[i][chg]);
      fHistoList->Add(fReconstructedNotPrimaries[i][chg]);
      fHistoList->Add(fReconstructedPrimaries[i][chg]);
      fHistoList->Add(fContamination[i][chg]);
    }

  }
  
  
  //********** PID ****************

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fpidResponse = inputHandler->GetPIDResponse();
  std::cout<<"*******"<< fpidResponse<<std::endl;




  
  // ************************

  PostData(1, fHistoList);
}


//_____________________________________________________________________

bool IsKaonNSigma(float mom, float nsigmaTPCK, float nsigmaTOFK)
{

  if(mom<0.25)
    {
      if(TMath::Abs(nsigmaTPCK)<2.0)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

  if(mom>=0.25 && mom<0.4)
    {
      if(TMath::Abs(nsigmaTPCK)<2.0)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

  if(mom>=0.4 && mom<0.45)
    {
      if(TMath::Abs(nsigmaTPCK)<1.0)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }

 if(mom>=0.45 && mom<0.5)
    {
      if(TMath::Abs(nsigmaTOFK)<2.0 && TMath::Abs(nsigmaTPCK)<3.0)
	{ 
	  return true;
	} 
      else 
	{
	  return false;
	}
    }
  
  if(mom>=0.5 && mom<0.8)
    {
      if(TMath::Abs(nsigmaTOFK)<2.0 && TMath::Abs(nsigmaTPCK)<3.0) 
	{
	  //     cout<<"500-800 "<<fNsigmaTOF500_800<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }

  if(mom>=0.8 && mom<1.0)
    {
      if(TMath::Abs(nsigmaTOFK)<1.5 && TMath::Abs(nsigmaTPCK)<3.0) 
	{
	  //          cout<<"800-1000 "<<fNsigmaTOF800_1000<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }

 if(mom>=1.0)
    {
      if(TMath::Abs(nsigmaTOFK)<1.0 && TMath::Abs(nsigmaTPCK)<3.0) 
	{
	  //      cout<<">1000 "<<fNsigmaTOFge1000<<" "<<nsigmaTOFK<<endl;
	  return true;
	}
      else
	{
	  return false;
	}
    }
    
    
    
//   if(mom>1.5 || mom<0.15)return false;
  return false;

}

bool IsPionNSigma(float mom, float nsigmaTPCPi, float nsigmaTOFPi)
{
  if (mom > 0.5) {
    //        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
    if (TMath::Hypot( nsigmaTOFPi, nsigmaTPCPi ) < 3.0)
      return true;
  }
  else {
    if (TMath::Abs(nsigmaTPCPi) < 3.0)
      return true;
  }
  return false;
}

bool IsProtonNSigma(float mom, float nsigmaTPCP, float nsigmaTOFP)
{
    if (mom > 0.5) {
        if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP ) < 3.0)
            return true;
	}
    else {
        if (TMath::Abs(nsigmaTPCP) < 3.0)
            return true;
    }

  return false;
}

bool IsElectron(float nsigmaTPCE, float nsigmaTPCPi,float nsigmaTPCK, float nsigmaTPCP)
{
  if(TMath::Abs(nsigmaTPCE)<3 && TMath::Abs(nsigmaTPCPi)>3 && TMath::Abs(nsigmaTPCK)>3 && TMath::Abs(nsigmaTPCP)>3)
    return true;
  else
    return false;
}

//_______________________________________________________

void AliAnalysisTaskParticlePurityEfficiency_pik::UserExec(Option_t *)
{

   AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

   //AliInputEventHandler *aodH = dynamic_cast<AliInputEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) 
  return; 

   AliAODEvent *fAODs = aodH->GetEvent();
  fAODpidUtil = aodH->GetAODpidUtil();
  

  /***Get Event****/
  //AliESDEvent *esdEvent = dynamic_cast<AliESDEvent *>(InputEvent());
  AliAODEvent* aodEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!aodEvent) return;
  AliAODHeader *fAODheader = (AliAODHeader*)aodEvent->GetHeader();
  //Double_t mult = fAODheader->GetRefMultiplicity();

AliMultSelection *MultSelection = (AliMultSelection*)aodEvent->FindListObject("MultSelection");
 Double_t mult = 0.0;
if(!MultSelection) {
        cout << "AliMultSelection object not found!" << endl;
        return;
 }
 else mult = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
    if(mult < 0 || mult >= 80) return;

 cout<<"########fcentrality#######   "<<mult<<endl;
 fHistEv->Fill(mult); 
  


  // EVENT SELECTION ********************

  fHistQA[9]->Fill(1);


  //****** Multiplicity selection *********
  Int_t fcent = -999;  

  if (mult >= 0 && mult < 5) fcent = 0;
  else if (mult >= 5 && mult < 10) fcent = 1;
  else if (mult >= 10 && mult < 20) fcent = 2;
  else if (mult >= 20 && mult < 30) fcent = 3;
  else if (mult >= 30 && mult < 40) fcent = 4;
  else if (mult >= 40 && mult < 50) fcent = 5;
  else if (mult >= 50 && mult < 60) fcent = 6;
  else if (mult >= 60 && mult < 70) fcent = 7;
  else if (mult >= 70 && mult < 80) fcent = 8;
  else return;
 
  if (fcent == 0) fHistEvCuts[0]->Fill(1);
  else if(fcent == 1) fHistEvCuts[1]->Fill(1);
  else if(fcent == 2) fHistEvCuts[2]->Fill(1);
  else if(fcent == 3) fHistEvCuts[3]->Fill(1);
  else if(fcent == 4) fHistEvCuts[4]->Fill(1);
  else if(fcent == 5) fHistEvCuts[5]->Fill(1);
  else if(fcent == 6) fHistEvCuts[6]->Fill(1);
  else if(fcent == 7) fHistEvCuts[7]->Fill(1);
  else if(fcent == 8) fHistEvCuts[8]->Fill(1);

  
  //"ESDs/pass2/AOD049/*AliAOD.root");
  const AliAODVertex* vertex =(AliAODVertex*) aodEvent->GetPrimaryVertex();
  vertex->GetPosition(fV1);
  if (!vertex || vertex->GetNContributors()<=0) return;
  //if (TMath::Abs(vertex->GetZ()) > 10.) return;
    
  fHistQA[9]->Fill(2);
  if(fcent == 0) fHistEvCuts[0]->Fill(2);
  else if(fcent == 1) fHistEvCuts[1]->Fill(2);
  else if(fcent == 2) fHistEvCuts[2]->Fill(2);
  else if(fcent == 3) fHistEvCuts[3]->Fill(2);
  else if(fcent == 4) fHistEvCuts[4]->Fill(2);
  else if(fcent == 5) fHistEvCuts[5]->Fill(2);
  else if(fcent == 6) fHistEvCuts[6]->Fill(2);
  else if(fcent == 7) fHistEvCuts[7]->Fill(2);
  else if(fcent == 8) fHistEvCuts[8]->Fill(2);
  
  AliAnalysisUtils *anaUtil = new AliAnalysisUtils();
    
  Bool_t fpA2013 = kFALSE;
  Bool_t fMVPlp = kFALSE;
  Bool_t fisPileUp = kTRUE;                                                                                                                   
  Int_t fMinPlpContribMV = 0;                                                                                                                 
  Int_t fMinPlpContribSPD = 5;                                                                                                                
  Int_t  fMinPlpZdistSPD=0.8;
  Int_t  fnSigmaPlpZdistSPD=3.;
  Int_t  fnSigmaPlpDiamXYSPD=2.0;
  Int_t  fnSigmaPlpDiamZSPD=5.0;


  if(fpA2013)
    if(anaUtil->IsVertexSelected2013pA(aodEvent)==kFALSE) return;
 
  if(fMVPlp) anaUtil->SetUseMVPlpSelection(kTRUE);
  else anaUtil->SetUseMVPlpSelection(kFALSE);
 
  if(fMinPlpContribMV) anaUtil->SetMinPlpContribMV(fMinPlpContribMV);
  if(fMinPlpContribSPD) anaUtil->SetMinPlpContribSPD(fMinPlpContribSPD);

  if(fisPileUp)
    if(anaUtil->IsPileUpEvent(aodEvent)) return;

  delete anaUtil;   

  fHistQA[9]->Fill(3);
  if(fcent == 0) fHistEvCuts[0]->Fill(3);
  else if(fcent == 1) fHistEvCuts[1]->Fill(3);
  else if(fcent == 2) fHistEvCuts[2]->Fill(3);
  else if(fcent == 3) fHistEvCuts[3]->Fill(3);
  else if(fcent == 4) fHistEvCuts[4]->Fill(3);
  else if(fcent == 5) fHistEvCuts[5]->Fill(3);
  else if(fcent == 6) fHistEvCuts[6]->Fill(3);
  else if(fcent == 7) fHistEvCuts[7]->Fill(3);
  else if(fcent == 8) fHistEvCuts[8]->Fill(3);

  //TString vtxTtl = vertex->GetTitle();
  //if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = vertex->GetZ();
  if (TMath::Abs(zvtx) > 7.0) return;
  fHistQA[0]->Fill(zvtx);
  fHistQA[9]->Fill(4);
  if(fcent == 0) fHistEvCuts[0]->Fill(4);
  else if(fcent == 1) fHistEvCuts[1]->Fill(4);
  else if(fcent == 2) fHistEvCuts[2]->Fill(4);
  else if(fcent == 3) fHistEvCuts[3]->Fill(4);
  else if(fcent == 4) fHistEvCuts[4]->Fill(4);
  else if(fcent == 5) fHistEvCuts[5]->Fill(4);
  else if(fcent == 6) fHistEvCuts[6]->Fill(4);
  else if(fcent == 7) fHistEvCuts[7]->Fill(4);
  else if(fcent == 8) fHistEvCuts[8]->Fill(4);

  //**** getting MC array ******
  TClonesArray  *arrayMC;

  arrayMC = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));


  //get stack 
  //AliStack *mcStack = mcEvent->Stack();
  //if (!mcStack) return;
  //***********************


  // old vertex selection 
  /*const AliESDVertex *vertex = esdEvent->GetPrimaryVertex();
    if (vertex->GetNContributors() < 1) return;

    //z-vertex cut
    if (TMath::Abs(vertex->GetZ()) > 10.) return;
  

    const AliESDVertex *vtxSPD = esdEvent->GetPrimaryVertexSPD();*/
  // Double_t zVertex = vtxSPD->GetZ();

  //std::cout << "Event  Z vtx ==========> " << vertex->GetZ() <<std::endl;
  // centrality selection 
  //AliCentrality *centrality = aodEvent->GetHeader()->GetCentralityP();
  //if (centrality->GetQuality() != 0) return;
  //Double_t cent = centrality->GetCentralityPercentileUnchecked("V0M");
  //if(cent < 0 || cent > 100.) return;


//copying pid information for FB 128
  int labels[20000];
  for (int il=0; il<20000; il++) labels[il] = -1;

  // looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
  for (int i=0;i<aodEvent->GetNumberOfTracks();i++) {
    const AliAODTrack *aodtrack=(AliAODTrack*)aodEvent->GetTrack(i);
    if (!aodtrack->TestFilterBit(7)) {
      if(aodtrack->GetID() < 0) continue;
      labels[aodtrack->GetID()] = i;
    }
  }



  //RECONSTRUCTED TRACKS 

  TObjArray recoParticleArray[PARTTYPES];

  fHistQA[10]->Fill(1,aodEvent->GetNumberOfTracks());
  //loop over AOD tracks 
  for (Int_t iTracks = 0; iTracks < aodEvent->GetNumberOfTracks(); iTracks++) {
    //get track 
    
    //AliESDtrack* track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent*>(esdEvent),iTracks);
    AliAODTrack *track = (AliAODTrack*)aodEvent->GetTrack(iTracks); 
    if (!track)continue;
    fHistQA[10]->Fill(2);

    //UInt_t filterBit = (1 << (0));
    UInt_t filterBit = 7;
    if(!track->TestFilterBit(filterBit))continue;	

    Bool_t fTrackPileUpRemoval = 1;    
    if (fTrackPileUpRemoval) {
      //method which checks if track
      //have at least 1 hit in ITS or TOF.
      bool passTrackPileUp = false;

      // does tof timing exist for our track?
      if (track->GetTOFBunchCrossing() == 0) {
        passTrackPileUp = true;
      }

      // check ITS refit
      if (!(track->GetStatus() & AliESDtrack::kITSrefit)) {
        continue;
      }

      // loop over 2 ITS layers and check for a hit!
      for (int i : {0, 1}) {
        if (track->HasPointOnITSLayer(i)) {
          passTrackPileUp = true;
        }
      }

      if (!passTrackPileUp) {
        continue;
      }
    }

    
    //charge
    //  if(track->Charge() < 0 ) continue;
    Int_t charge = 0;
    if(track->Charge() > 0 ) charge=0;
    else if (track->Charge() < 0 ) charge=1; 
    //if(!track->IsHybridGlobalConstrainedGlobal())continue;
    //if((track->IsHybridGlobalConstrainedGlobal())==false)continue;
    // if(!track->IsHybridTPCConstrainedGlobal())continue;	
    // if(!track->IsTPCConstrained())continue;	
    //if(!track->IsGlobalConstrained())continue;
    //if((track->TestFilterMask(AliAODTrack::kTrkTPCOnly)==false))continue;//cut0_BIT(0)
  
    //   if((track->IsHybridGlobalConstrainedGlobal())==false)
    //  continue;//def_BIT(272)

    //if((track->TestFilterMask(AliAODTrack::kTrkGlobal)==false))continue;//cut1_BIT(5)

    fHistQA[10]->Fill(3);
     
    if(track->Eta() < -0.8 || track->Eta() > 0.8)
      continue; 
    fHistQA[10]->Fill(4);

    if (track->Pt() < 0.1 || track->Pt() > 2.5)
      continue;
    fHistQA[10]->Fill(5);

    //single track cuts
    // if(track->Chi2perNDF() > 4.0) continue;
    if(track->GetTPCNcls() < 70) continue;

    //DCA
    
    Double_t DCAXY;
    Double_t DCAZ;
    //  if(filterBit==(1 << (7))){
    //    DCAXY = -TMath::Abs(track->DCA());
    //DCAZ = -TMath::Abs(track->ZAtDCA());
 
    //if(!(DCAXY==-999 || DCAZ==-999)){
	//if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
	//no DCA cut
	//if(TMath::Abs(DCAXY) > 1000.0) {continue;} //XY
	//if(TMath::Abs(DCAZ) > 1000.0) {continue;} //Z
    //}
    //else {
      // code from Michael and Prabhat from AliAnalysisTaskDptDptCorrelations
      // const AliAODVertex* vertex = (AliAODVertex*) aodEvent->GetPrimaryVertex(); (already defined above)
      float vertexX  = -999.;
      float vertexY  = -999.;
      float vertexZ  = -999.;

      if(vertex) {
	Double32_t fCov[6];
	vertex->GetCovarianceMatrix(fCov);
	if(vertex->GetNContributors() > 0) {
	  if(fCov[5] != 0) {
	    vertexX = vertex->GetX();
	    vertexY = vertex->GetY();
	    vertexZ = vertex->GetZ();

	  }
	}
      }

      Double_t pos[3];
      track->GetXYZ(pos);

      Double_t DCAX = pos[0] - vertexX;
      Double_t DCAY = pos[1] - vertexY;
      DCAZ = pos[2] - vertexZ;
      DCAXY = TMath::Sqrt((DCAX*DCAX) + (DCAY*DCAY));


      if (TMath::Abs(DCAXY) > 0.3)
      continue;
      if (TMath::Abs(DCAZ) > 0.3)
      continue;
      
      //if(TMath::Abs(DCAXY) > 0.0182 + 0.035*TMath::Power(track->Pt(), -1.01)) continue; //XY, Pt dep
      //if(TMath::Abs(impactD) > 0.44 + 0.07*TMath::Power(tPt, -1.94)) continue; //XY, Pt dep
      //no DCA cut
      //if(TMath::Abs(DCAXY) > 1000.0) continue;
      //if(TMath::Abs(DCAZ) > 1000.0) continue;
      //}

    fHistQA[10]->Fill(6);

    AliAODTrack* aodtrackpid;

    //cout << "I am here" << endl;
    //cout << (1 << (7)) << endl;
    //for FB 128 - tpc only tracks
    if(filterBit==(1 << (7))) {
      aodtrackpid =(AliAODTrack*)aodEvent->GetTrack(labels[-1-aodEvent->GetTrack(iTracks)->GetID()]);
      //cout << "I am getting track from aodEvent with new labels" << endl;
    }
    else {
      aodtrackpid = track;
      //cout << "I am getting track from aodEvent with old labels" << endl;
    }
    
    //if (!aodtrackpid) {
    //  cout << "aodtrackpid does not exist!" << endl;
    //}
   
    //Electron rejection
    double nSigmaTPCPi = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kPion);
    double nSigmaTPCK = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kKaon);
    double nSigmaTPCP = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kProton);
    double nSigmaTPCe = fpidResponse->NumberOfSigmasTPC(aodtrackpid,AliPID::kElectron);
    if(IsElectron(nSigmaTPCe,nSigmaTPCPi,nSigmaTPCK,nSigmaTPCP))
      continue;
    /*  fHistQA[10]->Fill(7);*/
    
    fHistQA[1]->Fill(track->GetTPCClusterInfo(2,1)); 
    //fHistQA[2]->Fill(track->GetTPCNclsF());
     fHistQA[3]->Fill(DCAXY);
     fHistQA[4]->Fill(DCAZ);
    Float_t chi2Tpc = track->Chi2perNDF();
    fHistQA[5]->Fill(chi2Tpc);
    fHistQA[6]->Fill(track->Pt());

    float px=track->Px(); float py=track->Py();  float ph=atan2(py,px); //track->Phi()
    float tPt = track->Pt();

    fHistQA[7]->Fill(ph);
    fHistQA[8]->Fill(track->Eta());
    fHistQA2D[2]->Fill(track->Eta(),ph);

    fHistQA2D[0]->Fill(tPt,DCAXY);
    fHistQA2D[1]->Fill(tPt,DCAZ);

    //PID monitors
    double nSigmaTOFPi = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kPion);
    double nSigmaTOFK = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kKaon);
    double nSigmaTOFP = fpidResponse->NumberOfSigmasTOF(aodtrackpid,AliPID::kProton);

    float tdEdx = aodtrackpid->GetTPCsignal();
    float tTofSig = aodtrackpid->GetTOFsignal();
    double pidTime[5]; aodtrackpid->GetIntegratedTimes(pidTime);


    fHistQAPID[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPID[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPID[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPID[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPID[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);

    fHistQAPIDFail[0][0][charge]->Fill(tPt,tdEdx);
    fHistQAPIDFail[1][0][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
    fHistQAPIDFail[2][0][charge]->Fill(tPt,nSigmaTOFPi);
    fHistQAPIDFail[3][0][charge]->Fill(tPt,nSigmaTPCPi);
    fHistQAPIDFail[4][0][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);


    bool isPionNsigma = 0;
    bool isKaonNsigma = 0;
    bool isProtonNsigma  = 0;

    //******** With double counting *******************
    // isPionNsigma = (IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi));
    // isKaonNsigma = (IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK));
    // isProtonNsigma = (IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));

    //******** Without double counting *******************
    // double nSigmaPIDPi = 0, nSigmaPIDK = 0, nSigmaPIDP = 0;
    // if(track->Pt()<0.5){
    //   nSigmaPIDPi = abs(nSigmaTPCPi);
    //   nSigmaPIDK  = abs(nSigmaTPCK);
    //   nSigmaPIDP  = abs(nSigmaTPCP);
    // }
    // else{
    //   nSigmaPIDPi = TMath::Hypot(nSigmaTOFPi, nSigmaTPCPi);
    //   nSigmaPIDK= TMath::Hypot(nSigmaTOFK, nSigmaTPCK);
    //   nSigmaPIDP= TMath::Hypot(nSigmaTOFP, nSigmaTPCP);
    // }
    // if(nSigmaPIDPi<nSigmaPIDK && nSigmaPIDPi<nSigmaPIDP){
    //   isPionNsigma = (IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi));
    // }
    // else if(nSigmaPIDK<nSigmaPIDPi && nSigmaPIDK<nSigmaPIDP){
    //   isKaonNsigma = (IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK));
    // }
    // else if(nSigmaPIDP<nSigmaPIDPi && nSigmaPIDP<nSigmaPIDK){
    //   isProtonNsigma = (IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));
    // }

    //******** Exclusive PID ********************
    //isPionNsigma = (IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi) && !IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));
    //isKaonNsigma = (!IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi)  && IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));
    //isProtonNsigma = (!IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi)  && !IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));


    isPionNsigma = (IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi) && !IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));
    isKaonNsigma = (!IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi)  && IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && !IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));
    isProtonNsigma = (!IsPionNSigma(track->Pt(),nSigmaTPCPi, nSigmaTOFPi)  && !IsKaonNSigma(track->Pt(),nSigmaTPCK, nSigmaTOFK) && IsProtonNSigma(track->Pt(),nSigmaTPCP, nSigmaTOFP));

    if (isPionNsigma){
      //      if (track->Pt() > 1.5) continue;
      fHistQAPID[0][1][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
      fHistQAPID[2][1][charge]->Fill(tPt,nSigmaTOFPi);
      fHistQAPID[3][1][charge]->Fill(tPt,nSigmaTPCPi);
      fHistQAPID[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
    }
    else
      {
	//	if (track->Pt() > 1.5) continue;
	fHistQAPIDFail[0][1][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][1][charge]->Fill(tPt,tTofSig-pidTime[2]);//pion
	fHistQAPIDFail[2][1][charge]->Fill(tPt,nSigmaTOFPi);
	fHistQAPIDFail[3][1][charge]->Fill(tPt,nSigmaTPCPi);
	fHistQAPIDFail[4][1][charge]->Fill(nSigmaTPCPi,nSigmaTOFPi);
      }
    if (isKaonNsigma){
      //if (track->Pt() > 1.5) continue;
      fHistQAPID[0][2][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
      fHistQAPID[2][2][charge]->Fill(tPt,nSigmaTOFK);
      fHistQAPID[3][2][charge]->Fill(tPt,nSigmaTPCK);
      fHistQAPID[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
    }
    else
      {
	//if (track->Pt() > 1.5) continue;	
	fHistQAPIDFail[0][2][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][2][charge]->Fill(tPt,tTofSig-pidTime[3]);//kaon
	fHistQAPIDFail[2][2][charge]->Fill(tPt,nSigmaTOFK);
	fHistQAPIDFail[3][2][charge]->Fill(tPt,nSigmaTPCK);
	fHistQAPIDFail[4][2][charge]->Fill(nSigmaTPCK,nSigmaTOFK);
      }
    if (isProtonNsigma){
      //      if ((track->Pt() < 0.5) ||  (track->Pt() > 4.0)) continue;      
      fHistQAPID[0][3][charge]->Fill(tPt,tdEdx);
      fHistQAPID[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
      fHistQAPID[2][3][charge]->Fill(tPt,nSigmaTOFP);
      fHistQAPID[3][3][charge]->Fill(tPt,nSigmaTPCP);
      fHistQAPID[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
    }
    else
      {
	//if  ((track->Pt() < 0.5) ||  (track->Pt() > 4.0)) continue;      	
	fHistQAPIDFail[0][3][charge]->Fill(tPt,tdEdx);
	fHistQAPIDFail[1][3][charge]->Fill(tPt,tTofSig-pidTime[4]);//proton
	fHistQAPIDFail[2][3][charge]->Fill(tPt,nSigmaTOFP);
	fHistQAPIDFail[3][3][charge]->Fill(tPt,nSigmaTPCP);
	fHistQAPIDFail[4][3][charge]->Fill(nSigmaTPCP,nSigmaTOFP);
      }

      fReconstructedAfterCuts[fcent][charge]->Fill(track->Eta(), track->Pt());//Fills hist. for all reconstructed particles after cuts
 

    if(!arrayMC){
      continue;
    }
    //get coresponding MC particle 
    Int_t label = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label);
    
    fProfilePTrueReconstructed[0]->Fill(track->P(), MCtrk->P() - track->P());
    fProfileThetaTrueReconstructed[0]->Fill(track->P(), MCtrk->Theta() - track->Theta());
    fProfilePhiTrueReconstructed[0]->Fill(track->P(), MCtrk->Phi() - track->Phi());
    
    fHist2DPTrueReconstructed[0]->Fill(track->P(), MCtrk->P() - track->P());
    fHist2DThetaTrueReconstructed[0]->Fill(track->P(), MCtrk->Theta() - track->Theta());
    fHist2DPhiTrueReconstructed[0]->Fill(track->P(), MCtrk->Phi() - track->Phi());

   //getting no. of tracks for each particle species after all the cuts:

    //********* PID - pions ********
     if (isPionNsigma){
       //      if (track->Pt() > 1.5) continue;             
       fReconstructedAfterCuts[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
       fProfilePTrueReconstructed[1]->Fill(track->P(), MCtrk->P() - track->P());
       fProfileThetaTrueReconstructed[1]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fProfilePhiTrueReconstructed[1]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       fHist2DPTrueReconstructed[1]->Fill(track->P(), MCtrk->P() - track->P());
       fHist2DThetaTrueReconstructed[1]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fHist2DPhiTrueReconstructed[1]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       if (!MCtrk) continue;
       recoParticleArray[1].Add(MCtrk);
       }
       //Fills for all identified pions found after cuts (reconstructed) - numerator for Efficiency

     //********* PID - kaons ********
     if (isKaonNsigma){
       //      if (track->Pt() > 1.5) continue;                    
       fReconstructedAfterCuts[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
       fProfilePTrueReconstructed[2]->Fill(track->P(), MCtrk->P() - track->P());
       fProfileThetaTrueReconstructed[2]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fProfilePhiTrueReconstructed[2]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       fHist2DPTrueReconstructed[2]->Fill(track->P(), MCtrk->P() - track->P());
       fHist2DThetaTrueReconstructed[2]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fHist2DPhiTrueReconstructed[2]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       if (!MCtrk) continue;
       recoParticleArray[2].Add(MCtrk);
       }
       //Fills for all identified kaons found after cuts (reconstructed) - numerator for Efficiency

    //********* PID - protons ********
     if (isProtonNsigma){
       //if ((track->Pt() < 0.5) ||  (track->Pt() > 4.0)) continue;                    
       fReconstructedAfterCuts[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
       fProfilePTrueReconstructed[3]->Fill(track->P(), MCtrk->P() - track->P());
       fProfileThetaTrueReconstructed[3]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fProfilePhiTrueReconstructed[3]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       fHist2DPTrueReconstructed[3]->Fill(track->P(), MCtrk->P() - track->P());
       fHist2DThetaTrueReconstructed[3]->Fill(track->P(), MCtrk->Theta() - track->Theta());
       fHist2DPhiTrueReconstructed[3]->Fill(track->P(), MCtrk->Phi() - track->Phi());
       if (!MCtrk) continue;
       recoParticleArray[3].Add(MCtrk);
       }

      //Fills for all identified protos found after cuts (reconstructed) - numerator for Efficiency
   //******************************

     //get coresponding MC particle 
     // Int_t label = TMath::Abs(track->GetLabel()); //moved up
     // if(!label) std::cout<<"no label"<<std::endl;
     //if(label) std::cout<<"label = "<<label<<std::endl;
       
    //AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(label); //moved up
    if (!MCtrk) continue;
    if(MCtrk->Charge()==0){std::cout<<"!!!"<<std::endl; continue;}
    recoParticleArray[0].Add(MCtrk);


    //Fills histogram for particles that are contamination from secondaries:
    if (!MCtrk->IsPhysicalPrimary()) {
      fReconstructedNotPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }
    else{
      fReconstructedPrimaries[PARTTYPES*fcent][charge]->Fill(track->Eta(), track->Pt());
    }

 
    int PDGcode = MCtrk->GetPdgCode();

    if (isPionNsigma && PDGcode==211){
    fHistQA2D[3]->Fill(DCAXY,tPt);
    fHistQA2D[4]->Fill(DCAZ,tPt);
    }

    if (isKaonNsigma && PDGcode==321){
    fHistQA2D[5]->Fill(DCAXY,tPt);
    fHistQA2D[6]->Fill(DCAZ,tPt);
    }

    if (isPionNsigma && PDGcode==-211){
    fHistQA2D[7]->Fill(DCAXY,tPt);
    fHistQA2D[8]->Fill(DCAZ,tPt);
    }

    if (isKaonNsigma && PDGcode==-321){
      fHistQA2D[9]->Fill(DCAXY,tPt);
    fHistQA2D[10]->Fill(DCAZ,tPt);
    }

    
   //And secondaries for different particle species:
    if (!MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) { //secondaries in pions
      fReconstructedNotPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isPionNsigma && abs(PDGcode)==211)) {
      fReconstructedPrimaries[PARTTYPES*fcent+1][charge]->Fill(track->Eta(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) { //secondaries in kaons
      fReconstructedNotPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
    }
    else if(MCtrk->IsPhysicalPrimary() && (isKaonNsigma && abs(PDGcode)==321)) {
      fReconstructedPrimaries[PARTTYPES*fcent+2][charge]->Fill(track->Eta(), track->Pt());
    }

    if (!MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) { //secondaries in protons
      fReconstructedNotPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
    } 
    else if(MCtrk->IsPhysicalPrimary() && (isProtonNsigma && abs(PDGcode)==2212)) {
      fReconstructedPrimaries[PARTTYPES*fcent+3][charge]->Fill(track->Eta(), track->Pt());
    } 

    //fMomentumTrueRecDiff->Fill(track->P(), track->P() - track->

    //step 1, TOF Matching
    UInt_t status;
    status=track->GetStatus();
    if((status&AliVTrack::kTOFout)==0 || (status&AliVTrack::kTIME)==0)
      status=0;
    if(track->Pt()<0.5) status = 1;

    //Misidentification fraction
    //	if(isPionNsigma && (MCtrk->Pt() <= 1.5))
    if(isPionNsigma)
    {
      if(abs(PDGcode)==211)
	   fMisidentification[fcent][charge]-> Fill(1,0.5,track->Pt());
      if(abs(PDGcode)==321)
	  fMisidentification[fcent][charge]-> Fill(1,1.5,track->Pt());
	if(abs(PDGcode)==2212)
	  fMisidentification[fcent][charge]-> Fill(1,2.5,track->Pt());
	if(!(abs(PDGcode)==211) && !(abs(PDGcode)==321) && !(abs(PDGcode)==2212))
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(1,3.5,track->Pt());

      }

    //  else if(isKaonNsigma && (MCtrk->Pt() <= 1.5))
      else if(isKaonNsigma)
      {
	if(abs(PDGcode)==211)
	   fMisidentification[fcent][charge]-> Fill(2,0.5,track->Pt());
	if(abs(PDGcode)==321)
	  fMisidentification[fcent][charge]-> Fill(2,1.5,track->Pt());
	if(abs(PDGcode)==2212)
	  fMisidentification[fcent][charge]-> Fill(2,2.5,track->Pt());
	if(!(abs(PDGcode)==211) && !(abs(PDGcode)==321) && !(abs(PDGcode)==2212))
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(2,3.5,track->Pt());

      }


    //    else if(isProtonNsigma && (MCtrk->Pt() > 0.5)  && (MCtrk->Pt() <= 4.0))
        else if(isProtonNsigma)
      {
	if(abs(PDGcode)==211)
	   fMisidentification[fcent][charge]-> Fill(3,0.5,track->Pt());
	    if(abs(PDGcode)==321)
	  fMisidentification[fcent][charge]-> Fill(3,1.5,track->Pt());
	   if(abs(PDGcode)==2212)
	  fMisidentification[fcent][charge]-> Fill(3,2.5,track->Pt());
	if(!(abs(PDGcode)==211) && !(abs(PDGcode)==321) && !(abs(PDGcode)==2212))
	  if(status)
	    fMisidentification[fcent][charge]-> Fill(3,3.5,track->Pt());

      }



      fContamination[PARTTYPES*fcent][charge]-> Fill(PDGcode,track->Pt()); 
    //Contaminations: "how many pions are in the kaons sample"? etc.
    //Do not use for corrections: using those values will be dependant on i.e. Pi/K ratio in MC
    //Use misidentification fraction instead
    //if(isPionNsigma && (MCtrk->Pt() <= 1.5))
      if(isPionNsigma)
      {
	fContamination[PARTTYPES*fcent+1][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for pions
      }
//    if(isKaonNsigma && (MCtrk->Pt() <= 1.5))
      if(isKaonNsigma)
      {
	fContamination[PARTTYPES*fcent+2][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for kaons
      }
      //if(isProtonNsigma && (MCtrk->Pt() > 0.5) && (MCtrk->Pt() <= 4.0))
      if(isProtonNsigma)
      {
	fContamination[PARTTYPES*fcent+3][charge]-> Fill(PDGcode,track->Pt()); // filling contamination histogram for protons
      }
    
  }
  
  // MONTECARLO PARTICLES 
  if(!arrayMC){
    AliError("Array of MC particles not found");
    return;
  }
  // loop over MC stack 
  for (Int_t ipart = 0; ipart < arrayMC->GetEntries(); ipart++) {
    //std::cout<<"Entered MC loop"<<std::endl;
    
    AliAODMCParticle *MCtrk = (AliAODMCParticle*)arrayMC->At(ipart);

    if (!MCtrk) continue;
    //std::cout<<"particle obtained"<<std::endl;
    
    Int_t PDGcode = TMath::Abs(MCtrk->GetPdgCode()); 

    
    //if(MCtrk->Charge() == 0) continue;
    Int_t charge=0;

    if(MCtrk->Charge() < 0) charge=1;
    else if(MCtrk->Charge() > 0) charge=0;


    if(MCtrk->Charge() == 0)
      {	
	if(MCtrk->GetPdgCode() == 3122) charge = 0;
	else if(MCtrk->GetPdgCode() == -3122) charge = 1;
      }

    
     
    //*** PID - check if pion ***
    //if(PDGcode!=211) continue; //(PDGcode==11 || PDGcode==321 || PDGcode==2212 || PDGcode==13)

      if(MCtrk->Eta() < -0.8 || MCtrk->Eta() > 0.8){
	continue; }
	
      if (MCtrk->Pt() < 0.1 || MCtrk->Pt() > 2.5){
	continue;}


      
      // check physical primary 
      if(MCtrk->IsPhysicalPrimary()) // Not from weak decay!
	{

	// Filling histograms for MC truth particles
	fGeneratedMCPrimaries[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	//if(PDGcode==211 && (MCtrk->Pt() <= 1.5))
	if(abs(PDGcode)==211)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	//else if(PDGcode==321 && (MCtrk->Pt() <= 1.5))
	else if(abs(PDGcode)==321)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	//else if(PDGcode==2212 && (MCtrk->Pt() > 0.5)  && (MCtrk->Pt() <= 4.0))
	else if(abs(PDGcode)==2212)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	else if(abs(PDGcode)==3122)
	  fGeneratedMCPrimaries[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());

	  //Filling data from MC truth particles only for particles that were reconstruced
	if (recoParticleArray[0].Contains(MCtrk)){ //All
	  fMCPrimariesThatAreReconstructed[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());

	  fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  //if(PDGcode==211 && (MCtrk->Pt() <= 1.5))
	  if(abs(PDGcode)==211)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  //if(PDGcode==321 && (MCtrk->Pt() <= 1.5))
	  if(abs(PDGcode)==321)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  //if(PDGcode==2212 && (MCtrk->Pt() > 0.5)&& (MCtrk->Pt() <= 4.0))
	  if(abs(PDGcode)==2212)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	  if(PDGcode==-3212)
	    fMCPrimariesThatAreReconstructedNoNsigma[fcent*PARTTYPES+4][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[1].Contains(MCtrk)){ //Pions
	  //if(PDGcode==211 && (MCtrk->Pt() <= 1.5))
	  if(abs(PDGcode)==211)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+1][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[2].Contains(MCtrk)){ //Kaons
	  //if(PDGcode==321 && (MCtrk->Pt() <= 1.5))
	  if(abs(PDGcode)==321)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+2][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}
	if (recoParticleArray[3].Contains(MCtrk)){ //Protons
	  //if(PDGcode==2212 && (MCtrk->Pt() > 0.5) && (MCtrk->Pt() <= 4.0))
	  if(abs(PDGcode)==2212)
	    fMCPrimariesThatAreReconstructed[fcent*PARTTYPES+3][charge]->Fill(MCtrk->Eta(), MCtrk->Pt());
	}

      }
    
  }
  PostData(1, fHistoList);
}
//-----------------------------------------------------------------

//void AliAnalysisTaskParticleEff::Terminate(Option_t *) 
//{}


