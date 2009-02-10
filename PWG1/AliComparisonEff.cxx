//------------------------------------------------------------------------------
// Implementation of AliComparisonEff class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is 
// a data member of AliComparisonEff.
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_ROOT/PWG1/Macros/LoadMyLibs.C");
  LoadMyLibs();
  TFile f("Output.root");
  //AliComparisonEff * compObj = (AliComparisonEff*)f.Get("AliComparisonEff");
  AliComparisonEff * compObj = (AliComparisonEff*)cOutput->FindObject("AliComparisonEff");

  // Analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderEff" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_Eff.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/


#include <iostream>

#include "TFile.h"
#include "TCint.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph2D.h"
#include "TCanvas.h"
#include "TGraph.h"
// 
#include "AliESDEvent.h"  
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDfriendTrack.h"
#include "AliRecInfoCuts.h" 
#include "AliMCInfoCuts.h" 
#include "AliLog.h" 
//
#include "AliMathBase.h"
#include "AliTreeDraw.h" 
#include "AliMagF.h" 
#include "AliESDVertex.h" 
#include "AliExternalTrackParam.h" 
#include "AliTracker.h" 

#include "AliMCInfo.h" 
#include "AliESDRecInfo.h" 
#include "AliComparisonEff.h" 

using namespace std;


ClassImp(AliComparisonEff)

//_____________________________________________________________________________
AliComparisonEff::AliComparisonEff():
  AliComparisonObject("AliComparisonEff"),

  // histograms
 
  fMCPt(0),
  fMCRecPt(0),
  fMCRecPrimPt(0),
  fMCRecSecPt(0),

  fEffTPCPt(0),
  fEffTPCPtMC(0),
  fEffTPCPtF(0),
  //
  fEffTPCPt_P(0),
  fEffTPCPtMC_P(0),
  fEffTPCPtF_P(0),
  //
  fEffTPCPt_Pi(0),
  fEffTPCPtMC_Pi(0),
  fEffTPCPtF_Pi(0),
  //
  fEffTPCPt_K(0),
  fEffTPCPtMC_K(0),
  fEffTPCPtF_K(0),
 
  fEffTPCTan(0),
  fEffTPCTanMC(0),
  fEffTPCTanF(0),
  //
  fEffTPCPtTan(0),
  fEffTPCPtTanMC(0),
  fEffTPCPtTanF(0),
  // TPC+ITS
  fEffTPCITSPt(0),
  fEffTPCITSTan(0),

  // Cuts 
  fCutsRC(0), 
  fCutsMC(0),

  fVertex(0),

  // histogram folder 
  fAnalysisFolder(0)
{
  // init vertex
  fVertex = new AliESDVertex();
  fVertex->SetXv(0.0); fVertex->SetYv(0.0); fVertex->SetZv(0.0); 

  for(Int_t i=0; i<4; ++i)
  {
    fTPCPtDCASigmaIdeal[i]=0;
    fTPCPtDCASigmaFull[i]=0;
    fTPCPtDCASigmaDay0[i]=0;

    fTPCPtDCAXY[i]=0;
    fTPCPtDCAZ[i]=0;

	fTPCPtDCASigmaIdealPid[i]=0;
	fTPCPtDCASigmaFullPid[i]=0;
	fTPCPtDCASigmaDay0Pid[i]=0;

	fTPCPtDCAXYPid[i]=0;   
	fTPCPtDCAZPid[i]=0; 
  }

  Init();
}

//_____________________________________________________________________________
AliComparisonEff::~AliComparisonEff(){

  // 
  if(fMCPt)  delete  fMCPt; fMCPt=0;
  if(fMCRecPt)  delete  fMCRecPt; fMCRecPt=0;
  if(fMCRecPrimPt)  delete  fMCRecPrimPt; fMCRecPrimPt=0;
  if(fMCRecSecPt)  delete  fMCRecSecPt; fMCRecSecPt=0;

  // 
  if(fEffTPCPt)   delete  fEffTPCPt;   fEffTPCPt=0;
  if(fEffTPCPtMC) delete  fEffTPCPtMC; fEffTPCPtMC=0;
  if(fEffTPCPtF)  delete  fEffTPCPtF;  fEffTPCPtF=0;

  // 
  if(fEffTPCPt_P)   delete  fEffTPCPt_P;   fEffTPCPt_P=0;
  if(fEffTPCPtMC_P) delete  fEffTPCPtMC_P; fEffTPCPtMC_P=0;
  if(fEffTPCPtF_P)  delete  fEffTPCPtF_P;  fEffTPCPtF_P=0;

  // 
  if(fEffTPCPt_Pi)   delete  fEffTPCPt_Pi;   fEffTPCPt_Pi=0;
  if(fEffTPCPtMC_Pi) delete  fEffTPCPtMC_Pi; fEffTPCPtMC_Pi=0;
  if(fEffTPCPtF_Pi)  delete  fEffTPCPtF_Pi;  fEffTPCPtF_Pi=0;

  // 
  if(fEffTPCPt_K)   delete  fEffTPCPt_K;   fEffTPCPt_K=0;
  if(fEffTPCPtMC_K) delete  fEffTPCPtMC_K; fEffTPCPtMC_K=0;
  if(fEffTPCPtF_K)  delete  fEffTPCPtF_K;  fEffTPCPtF_K=0;

  // 
  if(fEffTPCTan)   delete  fEffTPCTan;   fEffTPCTan=0;
  if(fEffTPCTanMC) delete  fEffTPCTanMC; fEffTPCTanMC=0;
  if(fEffTPCTanF)  delete  fEffTPCTanF;  fEffTPCTanF=0;

  //
  if(fEffTPCPtTan)   delete  fEffTPCPtTan;   fEffTPCPtTan=0;
  if(fEffTPCPtTanMC) delete  fEffTPCPtTanMC; fEffTPCPtTanMC=0;
  if(fEffTPCPtTanF)  delete  fEffTPCPtTanF;  fEffTPCPtTanF=0;

  if(fEffTPCITSPt)   delete  fEffTPCITSPt;   fEffTPCITSPt=0;
  if(fEffTPCITSTan)   delete  fEffTPCITSTan;   fEffTPCITSTan=0;

  for(Int_t i=0; i<4; ++i)
  {
    if(fTPCPtDCASigmaIdeal[i]) delete  fTPCPtDCASigmaIdeal[i];  fTPCPtDCASigmaIdeal[i]=0;
    if(fTPCPtDCASigmaFull[i]) delete  fTPCPtDCASigmaFull[i];  fTPCPtDCASigmaFull[i]=0;
    if(fTPCPtDCASigmaDay0[i]) delete  fTPCPtDCASigmaDay0[i];  fTPCPtDCASigmaDay0[i]=0;

    if(fTPCPtDCAXY[i]) delete  fTPCPtDCAXY[i];  fTPCPtDCAXY[i]=0;
    if(fTPCPtDCAZ[i]) delete  fTPCPtDCAZ[i];  fTPCPtDCAZ[i]=0;

 	if(fTPCPtDCASigmaIdealPid[i]) delete  fTPCPtDCASigmaIdealPid[i];  fTPCPtDCASigmaIdealPid[i]=0;
	if(fTPCPtDCASigmaFullPid[i]) delete  fTPCPtDCASigmaFullPid[i];  fTPCPtDCASigmaFullPid[i]=0;
	if(fTPCPtDCASigmaDay0Pid[i]) delete  fTPCPtDCASigmaDay0Pid[i];  fTPCPtDCASigmaDay0Pid[i]=0;

	if(fTPCPtDCAXYPid[i]) delete  fTPCPtDCAXYPid[i]; fTPCPtDCAXYPid[i]=0;
	if(fTPCPtDCAZPid[i]) delete   fTPCPtDCAZPid[i];  fTPCPtDCAZPid[i]=0;
  }

  if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
}

//_____________________________________________________________________________
void AliComparisonEff::Init(){

  // Init histograms
  //
  fMCPt = new TH1F("fMCPt","fMCPt",50,0.1,3);           
  fMCPt->SetXTitle("p_{t}");
  fMCPt->SetYTitle("yield");

  fMCRecPt = new TH1F("fMCRecPt","fMCRecPt",50,0.1,3);           
  fMCRecPt->SetXTitle("p_{t}");
  fMCRecPt->SetYTitle("yield");

  fMCRecPrimPt = new TH1F("fMCRecPrimPt","fMCRecPrimPt",50,0.1,3);           
  fMCRecPrimPt->SetXTitle("p_{t}");
  fMCRecPrimPt->SetYTitle("yield");

  fMCRecSecPt = new TH1F("fMCRecSecPt","fMCRecSecPt",50,0.1,3);           
  fMCRecSecPt->SetXTitle("p_{t}");
  fMCRecSecPt->SetYTitle("yield");

  // Efficiency as function of pt
  fEffTPCPt = new TProfile("Eff_pt","Eff_Pt",50,0.1,3);           
  fEffTPCPt->SetXTitle("p_{t}");
  fEffTPCPt->SetYTitle("TPC Efficiency");

// Efficiency as function of pt (TPC+ITS)
  fEffTPCITSPt = new TProfile("Eff_pt_TPCITS","Eff_Pt_TPCITS",50,0.1,3);           
  fEffTPCITSPt->SetXTitle("p_{t}");
  fEffTPCITSPt->SetYTitle("TPCITS Efficiency");


  fEffTPCPtMC = new TProfile("MC_Eff_pt","MC_Eff_Pt",50,0.1,3);   
  fEffTPCPtMC->SetXTitle("p_{t}");
  fEffTPCPtMC->SetYTitle("MC TPC Efficiency");

  fEffTPCPtF = new TProfile("F_Eff_pt","F_Eff_Pt",50,0.1,3);     
  fEffTPCPtF->SetXTitle("p_{t}");
  fEffTPCPtF->SetYTitle("TPC Findable Efficiency");

  // Efficiency as function of pt protons
  fEffTPCPt_P = new TProfile("Eff_pt_P","Eff_Pt_P",50,0.1,3);           
  fEffTPCPt_P->SetXTitle("p_{t}");
  fEffTPCPt_P->SetYTitle("TPC Efficiency");

  fEffTPCPtMC_P = new TProfile("MC_Eff_pt_P","MC_Eff_Pt_P",50,0.1,3);   
  fEffTPCPtMC_P->SetXTitle("p_{t}");
  fEffTPCPtMC_P->SetYTitle("MC TPC Efficiency");

  fEffTPCPtF_P = new TProfile("F_Eff_pt_P","F_Eff_Pt_P",50,0.1,3);     
  fEffTPCPtF_P->SetXTitle("p_{t}");
  fEffTPCPtF_P->SetYTitle("TPC Findable Efficiency");

  // Efficiency as function of pt pions
  fEffTPCPt_Pi = new TProfile("Eff_pt_Pi","Eff_Pit_Pi",50,0.1,3);           
  fEffTPCPt_Pi->SetXTitle("p_{t}");
  fEffTPCPt_Pi->SetYTitle("TPC Efficiency");

  fEffTPCPtMC_Pi = new TProfile("MC_Eff_pt_Pi","MC_Eff_Pit_Pi",50,0.1,3);   
  fEffTPCPtMC_Pi->SetXTitle("p_{t}");
  fEffTPCPtMC_Pi->SetYTitle("MC TPC Efficiency");

  fEffTPCPtF_Pi = new TProfile("F_Eff_pt_Pi","F_Eff_Pit_Pi",50,0.1,3);     
  fEffTPCPtF_Pi->SetXTitle("p_{t}");
  fEffTPCPtF_Pi->SetYTitle("TPC Findable Efficiency");

  // Efficiency as function of pt kaons
  fEffTPCPt_K = new TProfile("Eff_pt_K","Eff_Kt_K",50,0.1,3);           
  fEffTPCPt_K->SetXTitle("p_{t}");
  fEffTPCPt_K->SetYTitle("TPC Efficiency");

  fEffTPCPtMC_K = new TProfile("MC_Eff_pt_K","MC_Eff_Kt_K",50,0.1,3);   
  fEffTPCPtMC_K->SetXTitle("p_{t}");
  fEffTPCPtMC_K->SetYTitle("MC TPC Efficiency");

  fEffTPCPtF_K = new TProfile("F_Eff_pt_K","F_Eff_Kt_K",50,0.1,3);     
  fEffTPCPtF_K->SetXTitle("p_{t}");
  fEffTPCPtF_K->SetYTitle("TPC Findable Efficiency");

  // Efficiency as function of tan(theta) 
  fEffTPCTan = new TProfile("Eff_tan","Eff_tan",50,-2.5,2.5);         
  fEffTPCTan->SetXTitle("tan(#theta)");
  fEffTPCTan->SetYTitle("TPC Efficiency");

// Efficiency as function of tan(theta) (TPC+ITS)
  fEffTPCITSTan = new TProfile("Eff_tan_TPCITS","Eff_tan_TPCITS",50,-2.5,2.5);         
  fEffTPCITSTan->SetXTitle("tan(#theta)");
  fEffTPCITSTan->SetYTitle("TPCITS Efficiency");


  fEffTPCTanMC = new TProfile("MC_Eff_tan","MC_Eff_tan",50,-2.5,2.5); 
  fEffTPCTanMC->SetXTitle("tan(#theta)");
  fEffTPCTanMC->SetYTitle("MC TPC Efficiency");

  fEffTPCTanF = new TProfile("F_Eff_tan","F_Eff_tan",50,-2.5,2.5);   
  fEffTPCPtF->SetXTitle("tan(#theta)");
  fEffTPCPtF->SetYTitle("TPC Findable Efficiency");

  // Efficiency as function of pt and tan(theta) 
  fEffTPCPtTan = new TProfile2D("Eff_pt_tan","Eff_Pt_tan",10,0.1,3,20,-2.,2.);
  fEffTPCPtTan->SetXTitle("tan(#theta)");
  fEffTPCPtTan->SetYTitle("p_{t}");

  fEffTPCPtTanMC = new TProfile2D("MC_Eff_pt_tan_MC","MC Eff Pt",10,0.1,3,20,-2.,2.);
  fEffTPCPtTanMC->SetXTitle("tan(#theta)");
  fEffTPCPtTanMC->SetYTitle("p_{t}");

  fEffTPCPtTanF = new TProfile2D("MC_Eff_pt_tan_F","MC Eff Pt",10,0.1,3,20,-2.,2.);
  fEffTPCPtTanF->SetXTitle("tan(#theta)");
  fEffTPCPtTanF->SetYTitle("p_{t}");

  char name[256];
  for(Int_t i=0; i<4; ++i)
  {
    sprintf(name, "fTPCPtDCASigmaIdeal_%d",i);
    fTPCPtDCASigmaIdeal[i] = new TH2F(name,name,50,0.1,3,100,0,100);

    sprintf(name, "fTPCPtDCASigmaFull_%d",i);
    fTPCPtDCASigmaFull[i] = new TH2F(name,name,50,0.1,3,100,0,100);

    sprintf(name, "fTPCPtDCASigmaDay0_%d",i);
    fTPCPtDCASigmaDay0[i] = new TH2F(name,name,50,0.1,3,100,0,100);

    sprintf(name, "fTPCPtDCAXY_%d",i);
    fTPCPtDCAXY[i]= new TH2F(name,name,50,0.1,3,100,0,100);
    sprintf(name, "fTPCPtDCAZ_%d",i);
    fTPCPtDCAZ[i]= new TH2F(name,name,50,0.1,3,100,0,100);

    sprintf(name, "fTPCPtDCASigmaIdealPid_%d",i);
   	fTPCPtDCASigmaIdealPid[i] = new TH3F(name,name,50,0.1,3,100,0,100,5,0,5);

    sprintf(name, "fTPCPtDCASigmaFullPid_%d",i);
	fTPCPtDCASigmaFullPid[i]= new TH3F(name,name,50,0.1,3,100,0,100,5,0,5);

    sprintf(name, "fTPCPtDCASigmaDay0Pid_%d",i);
	fTPCPtDCASigmaDay0Pid[i]= new TH3F(name,name,50,0.1,3,100,0,100,5,0,5);

    sprintf(name, "fTPCPtDCAXYPid_%d",i);
	fTPCPtDCAXYPid[i]= new TH3F(name,name,50,0.1,3,100,0,100,5,0,5);

    sprintf(name, "fTPCPtDCAZPid_%d",i);
	fTPCPtDCAZPid[i]= new TH3F(name,name,50,0.1,3,100,0,100,5,0,5);
  }

  // init cuts
  if(!fCutsMC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  if(!fCutsRC) 
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object");

  // init folder
  fAnalysisFolder = CreateFolder("folderEff","Analysis Efficiency Folder");
}

//_____________________________________________________________________________
void AliComparisonEff::Process(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Fill efficiency comparison information
  
  AliExternalTrackParam *track = 0;
  Double_t field      = AliTracker::GetBz(); // nominal Bz field [kG]
  Double_t kMaxD      = 123456.0; // max distance
  Int_t clusterITS[200];

  // systematics
  const Double_t kSigma2Full_xy  = 0.25; // ExB full systematics  [cm]
  const Double_t kSigma2Full_z  =  5.0;  // drift velocity (goofie) [cm] 

  const Double_t kSigma2Day0_xy  = 0.02; //  ExB  [cm]
  const Double_t kSigma2Day0_z  =  0.1;  //  drift velocity (goofie) [cm]  

  //  
  Double_t  DCASigmaIdeal=0;
  Double_t  DCASigmaFull=0;
  Double_t  DCASigmaDay0=0;

  Double_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z

  // distance to Prim. vertex 
  const Double_t* dv = infoMC->GetVDist(); 

  Float_t mcpt = infoMC->GetParticle().Pt();
  Float_t tantheta = TMath::Tan(infoMC->GetParticle().Theta()-TMath::Pi()*0.5);
  Bool_t isPrim = TMath::Sqrt(dv[0]*dv[0] + dv[1]*dv[1])<fCutsMC->GetMaxR() && TMath::Abs(dv[2])<fCutsMC->GetMaxVz();
 
  // calculate and set prim. vertex
  fVertex->SetXv( infoMC->GetParticle().Vx() - dv[0] );
  fVertex->SetYv( infoMC->GetParticle().Vy() - dv[1] );
  fVertex->SetZv( infoMC->GetParticle().Vz() - dv[2] );
  
  // Check selection cuts
  if (fCutsMC->IsPdgParticle(TMath::Abs(infoMC->GetParticle().GetPdgCode())) == kFALSE) return; 

  // transform Pdg to Pid
  // Pdg convension is different for hadrons and leptons 
  // (e.g. K+/K- = 321/-321; e+/e- = -11/11 ) 
  Double_t pid = -1;
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetEM() ) pid = 0; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetMuM() ) pid = 1; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP() ) pid = 2; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP() ) pid = 3; 
  if( TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt() ) pid = 4; 

  //cout << "dv[0] " << dv[0] << " dv[1] " << dv[1]  <<  " dv[2] " << dv[2] << endl; 
  //cout << "v[0] " << fVertex->GetXv()  << " v[1] " << fVertex->GetYv()  <<  " v[2] " << fVertex->GetZv()<< endl; 
  if (TMath::Abs(tantheta)<fCutsRC->GetMaxAbsTanTheta())
  {
    if (infoRC->GetESDtrack() && infoRC->GetESDtrack()->GetTPCInnerParam()) 
    {
      if ((track = new AliExternalTrackParam(*infoRC->GetESDtrack()->GetTPCInnerParam())) != 0 )
      {
        Bool_t bDCAStatus = track->PropagateToDCA(fVertex,field,kMaxD,dca,cov);
		if(bDCAStatus) {
		//
		cov[2] = track->GetCovariance()[2];

	    // Eff = infoRC->GetStatus(1)==3 && isPrim / isPrim;
        // Pt vs ( dca[0]^2/cov[0]^2 + dca[1]^2/cov[2]^2 ) 
        // Pt vs ( dca[0]^2/(cov[0]^2 + kSigma2Full_xy)  + dca[1]^2/( cov[2]^2 + kSigma2Full_z ) 
        // Pt vs ( dca[0]^2/(cov[0]^2 + kSigma2_xy)  + dca[1]^2/( cov[2]^2 + kSigma2_z ) 

         if(cov[0]>0.0 && cov[2]>0.0)
		 {
			 DCASigmaIdeal = TMath::Power(dca[0],2)/cov[0] 
									  + TMath::Power(dca[1],2)/cov[2]; 

			 DCASigmaFull = TMath::Power(dca[0],2)/(cov[0]+kSigma2Full_xy) 
									 + TMath::Power(dca[1],2)/(cov[2]+kSigma2Full_z); 

			 DCASigmaDay0 = TMath::Power(dca[0],2)/(cov[0]+kSigma2Day0_xy) 
									 + TMath::Power(dca[1],2)/(cov[2]+kSigma2Day0_z); 

              //cout << "dca[0] " << dca[0]  << " dca[1] " << dca[1]  << endl; 
              //cout << "cov[0] " << cov[0]  << " cov[2] " << cov[2]  << endl; 
              //cout << "DCASigmaIdeal " << DCASigmaIdeal  << " DCASigmaFull " << DCASigmaFull  << " DCASigmaDay0 "  <<DCASigmaDay0 <<  endl; 
            //cout << " -------------------------------------------------------- "<<  endl; 
		 }

         // MC pt
         fMCPt->Fill(mcpt); 
         if(infoRC->GetStatus(1)==3) fMCRecPt->Fill(mcpt); 
         if(infoRC->GetStatus(1)==3 && isPrim) fMCRecPrimPt->Fill(mcpt); 
         if(infoRC->GetStatus(1)==3 && !isPrim) fMCRecSecPt->Fill(mcpt); 

         if(isPrim)
		 {
           fTPCPtDCASigmaIdeal[0]->Fill(mcpt,DCASigmaIdeal);
           fTPCPtDCASigmaFull[0]->Fill(mcpt,DCASigmaFull);
           fTPCPtDCASigmaDay0[0]->Fill(mcpt,DCASigmaDay0);

           fTPCPtDCAXY[0]->Fill(mcpt,dca[0]);
           fTPCPtDCAZ[0]->Fill(mcpt,dca[1]);

           fTPCPtDCASigmaIdealPid[0]->Fill(mcpt,DCASigmaIdeal,pid);
           fTPCPtDCASigmaFullPid[0]->Fill(mcpt,DCASigmaFull,pid);
           fTPCPtDCASigmaDay0Pid[0]->Fill(mcpt,DCASigmaDay0,pid);

           fTPCPtDCAXYPid[0]->Fill(mcpt,dca[0],pid);
           fTPCPtDCAZPid[0]->Fill(mcpt,dca[1],pid);
           		   
		   if(infoRC->GetStatus(1)==3)
		   {
             fTPCPtDCASigmaIdeal[1]->Fill(mcpt,DCASigmaIdeal);
             fTPCPtDCASigmaFull[1]->Fill(mcpt,DCASigmaFull);
             fTPCPtDCASigmaDay0[1]->Fill(mcpt,DCASigmaDay0);

             fTPCPtDCAXY[1]->Fill(mcpt,dca[0]);
             fTPCPtDCAZ[1]->Fill(mcpt,dca[1]);

             fTPCPtDCASigmaIdealPid[1]->Fill(mcpt,DCASigmaIdeal,pid);
             fTPCPtDCASigmaFullPid[1]->Fill(mcpt,DCASigmaFull,pid);
             fTPCPtDCASigmaDay0Pid[1]->Fill(mcpt,DCASigmaDay0,pid);

             fTPCPtDCAXYPid[1]->Fill(mcpt,dca[0],pid);
             fTPCPtDCAZPid[1]->Fill(mcpt,dca[1],pid);
           }
		 }

        // Cont = infoRC->GetStatus(1)==3 && !isPrim / infoRC->GetStatus(1)==3   
        // Pt vs ( dz[0]^2/cov[0]^2 + dz[1]^2/cov[2]^2 ) 
        // Pt vs ( dz[0]^2/(cov[0]^2 + kSigma2Full_xy)  + dz[1]^2/( cov[2]^2 + kSigma2Full_z ) 
        // Pt vs ( dz[0]^2/(cov[0]^2 + kSigma2_xy)  + dz[1]^2/( cov[2]^2 + kSigma2_z ) 

		 if(infoRC->GetStatus(1)==3) 
		 {
		   fTPCPtDCASigmaIdeal[2]->Fill(mcpt,DCASigmaIdeal);
           fTPCPtDCASigmaFull[2]->Fill(mcpt,DCASigmaFull);
           fTPCPtDCASigmaDay0[2]->Fill(mcpt,DCASigmaDay0);

           fTPCPtDCAXY[2]->Fill(mcpt,dca[0]);
           fTPCPtDCAZ[2]->Fill(mcpt,dca[1]);

           fTPCPtDCASigmaIdealPid[2]->Fill(mcpt,DCASigmaIdeal,pid);
           fTPCPtDCASigmaFullPid[2]->Fill(mcpt,DCASigmaFull,pid);
           fTPCPtDCASigmaDay0Pid[2]->Fill(mcpt,DCASigmaDay0,pid);

           fTPCPtDCAXYPid[2]->Fill(mcpt,dca[0],pid);
           fTPCPtDCAZPid[2]->Fill(mcpt,dca[1],pid);
 
		   if(isPrim==0)
		   {
             fTPCPtDCASigmaIdeal[3]->Fill(mcpt,DCASigmaIdeal);
             fTPCPtDCASigmaFull[3]->Fill(mcpt,DCASigmaFull);
             fTPCPtDCASigmaDay0[3]->Fill(mcpt,DCASigmaDay0);

             fTPCPtDCAXY[3]->Fill(mcpt,dca[0]);
             fTPCPtDCAZ[3]->Fill(mcpt,dca[1]);

             fTPCPtDCASigmaIdealPid[3]->Fill(mcpt,DCASigmaIdeal,pid);
             fTPCPtDCASigmaFullPid[3]->Fill(mcpt,DCASigmaFull,pid);
             fTPCPtDCASigmaDay0Pid[3]->Fill(mcpt,DCASigmaDay0,pid);

             fTPCPtDCAXYPid[3]->Fill(mcpt,dca[0],pid);
             fTPCPtDCAZPid[3]->Fill(mcpt,dca[1],pid);
           }
	     }
	   delete track;
       }
       }
     } 
	 else 
	 {
       if(isPrim)
	   {
	     fTPCPtDCASigmaIdeal[0]->Fill(mcpt,0.0);
         fTPCPtDCASigmaFull[0]->Fill(mcpt,0.0);
         fTPCPtDCASigmaDay0[0]->Fill(mcpt,0.0);

         fTPCPtDCAXY[0]->Fill(mcpt,0.0);
         fTPCPtDCAZ[0]->Fill(mcpt,0.0);

         fTPCPtDCASigmaIdealPid[0]->Fill(mcpt,0.0,pid);
         fTPCPtDCASigmaFullPid[0]->Fill(mcpt,0.0,pid);
         fTPCPtDCASigmaDay0Pid[0]->Fill(mcpt,0.0,pid);

         fTPCPtDCAXYPid[0]->Fill(mcpt,0.0,pid);
         fTPCPtDCAZPid[0]->Fill(mcpt,0.0,pid);
	   }
     }
  }

  // only primary particles
  if (!isPrim) return;

  // pt
  if (TMath::Abs(tantheta)<fCutsRC->GetMaxAbsTanTheta()){

    fEffTPCPt->Fill(mcpt, infoRC->GetStatus(1)==3);
    fEffTPCPtMC->Fill(mcpt, infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());
    if (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()){
      fEffTPCPtF->Fill(mcpt, infoRC->GetStatus(1)==3);
    }
  
    fEffTPCITSPt->Fill(mcpt, infoRC->GetStatus(1)==3 && (infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS()));

    // protons
    if(TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetProt()) { 
	   fEffTPCPt_P->Fill(mcpt, infoRC->GetStatus(1)==3);
       fEffTPCPtMC_P->Fill(mcpt, infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());

       if(infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()) {
         fEffTPCPtF_P->Fill(mcpt, infoRC->GetStatus(1)==3);
       }
	}

    // pions
    if(TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetPiP()) {
	  fEffTPCPt_Pi->Fill(mcpt, infoRC->GetStatus(1)==3);
      fEffTPCPtMC_Pi->Fill(mcpt, infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());

       if(infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()) {
         fEffTPCPtF_Pi->Fill(mcpt, infoRC->GetStatus(1)==3);
       }
	}

	// kaons
    if(TMath::Abs(infoMC->GetParticle().GetPdgCode())==fCutsMC->GetKP()) {
	  fEffTPCPt_K->Fill(mcpt, infoRC->GetStatus(1)==3);
      fEffTPCPtMC_K->Fill(mcpt, infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());

       if(infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()) {
         fEffTPCPtF_K->Fill(mcpt, infoRC->GetStatus(1)==3);
       }
    }
  }

  // theta
  if (TMath::Abs(mcpt)>fCutsRC->GetPtMin()){
    fEffTPCTan->Fill(tantheta, infoRC->GetStatus(1)==3);
    fEffTPCTanMC->Fill(tantheta, infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits());
    if (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()){
      fEffTPCTanF->Fill(tantheta, infoRC->GetStatus(1)==3);
    }

    fEffTPCITSTan->Fill(tantheta, infoRC->GetStatus(1)==3 && (infoRC->GetESDtrack()->GetITSclusters(clusterITS)>fCutsRC->GetMinNClustersITS()));
  }

  // pt-theta
  fEffTPCPtTan->Fill(mcpt,tantheta,infoRC->GetStatus(1)==3);
  fEffTPCPtTanMC->Fill(mcpt,tantheta,infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()); 
  if (infoMC->GetRowsWithDigits()>fCutsMC->GetMinRowsWithDigits()){
    fEffTPCPtTanF->Fill(mcpt,tantheta,infoRC->GetStatus(1)==3); 
  }
}

//_____________________________________________________________________________
void AliComparisonEff::Exec(AliMCInfo* infoMC, AliESDRecInfo *infoRC)
{
  // Process comparison information
  Process(infoMC,infoRC);
}

//_____________________________________________________________________________
Long64_t AliComparisonEff::Merge(TCollection* list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;

  // collection of generated histograms

  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliComparisonEff* entry = dynamic_cast<AliComparisonEff*>(obj);
    if (entry == 0) continue; 
  
    fMCPt->Add(entry->fMCPt);
    fMCRecPt->Add(entry->fMCRecPt);
    fMCRecPrimPt->Add(entry->fMCRecPrimPt);
    fMCRecSecPt->Add(entry->fMCRecSecPt);

    fEffTPCPt->Add(entry->fEffTPCPt);
	fEffTPCPtMC->Add(entry->fEffTPCPtMC);
	fEffTPCPtF->Add(entry->fEffTPCPtF);

    fEffTPCPt_P->Add(entry->fEffTPCPt_P);
	fEffTPCPtMC_P->Add(entry->fEffTPCPtMC_P);
	fEffTPCPtF_P->Add(entry->fEffTPCPtF_P);

    fEffTPCPt_Pi->Add(entry->fEffTPCPt_Pi);
	fEffTPCPtMC_Pi->Add(entry->fEffTPCPtMC_Pi);
	fEffTPCPtF_Pi->Add(entry->fEffTPCPtF_Pi);

    fEffTPCPt_K->Add(entry->fEffTPCPt_K);
	fEffTPCPtMC_K->Add(entry->fEffTPCPtMC_K);
	fEffTPCPtF_K->Add(entry->fEffTPCPtF_K);

	fEffTPCTan->Add(entry->fEffTPCTan);
	fEffTPCTanMC->Add(entry->fEffTPCTanMC);
	fEffTPCTanF->Add(entry->fEffTPCTanF);
	  
	fEffTPCPtTan->Add(entry->fEffTPCPtTan);
	fEffTPCPtTanMC->Add(entry->fEffTPCPtTanMC);
	fEffTPCPtTanF->Add(entry->fEffTPCPtTanF);
    
    fEffTPCITSPt->Add(entry->fEffTPCITSPt);
    fEffTPCITSTan->Add(entry->fEffTPCITSTan);

    for(Int_t i=0; i<4; ++i)
    {
      fTPCPtDCASigmaIdeal[i]->Add(entry->fTPCPtDCASigmaIdeal[i]);
      fTPCPtDCASigmaFull[i]->Add(entry->fTPCPtDCASigmaFull[i]);
      fTPCPtDCASigmaDay0[i]->Add(entry->fTPCPtDCASigmaDay0[i]);

      fTPCPtDCAXY[i]->Add(entry->fTPCPtDCAXY[i]);
      fTPCPtDCAZ[i]->Add(entry->fTPCPtDCAXY[i]);

      fTPCPtDCASigmaIdealPid[i]->Add(entry->fTPCPtDCASigmaIdealPid[i]);
      fTPCPtDCASigmaFullPid[i]->Add(entry->fTPCPtDCASigmaFullPid[i]);
      fTPCPtDCASigmaDay0Pid[i]->Add(entry->fTPCPtDCASigmaDay0Pid[i]);

      fTPCPtDCAXYPid[i]->Add(entry->fTPCPtDCAXYPid[i]);
      fTPCPtDCAZPid[i]->Add(entry->fTPCPtDCAXYPid[i]);
    }

  count++;
  }

return count;
}
 
//_____________________________________________________________________________
void AliComparisonEff::Analyse() 
{
  // Analyse comparison information and store output histograms
  // in the folder "folderEff" 
  //
  
  TH1::AddDirectory(kFALSE);

  AliComparisonEff * comp=this;
  TObjArray *aFolderObj = new TObjArray;

  // calculate efficiency and contamination (4 sigma) 
  TH1 *h_sigmaidealpid[20];
  TH1 *h_sigmafullpid[20];
  TH1 *h_sigmaday0pid[20];

  TH1 *h_sigmaidealpidtot[4];
  TH1 *h_sigmafullpidtot[4];
  TH1 *h_sigmaday0pidtot[4];

  char name[256];
  char name1[256];
  Int_t idx;

  for(Int_t i=0; i<4; ++i)
  {
     //total
     comp->fTPCPtDCASigmaIdealPid[i]->GetYaxis()->SetRange(1,4);
     sprintf(name,"h_sigmaidealpidtot_%d",i);
     h_sigmaidealpidtot[i] = comp->fTPCPtDCASigmaIdealPid[i]->Project3D();
     h_sigmaidealpidtot[i]->SetName(name);

     comp->fTPCPtDCASigmaFullPid[i]->GetYaxis()->SetRange(1,4);
     sprintf(name,"h_sigmafullpidtot_%d",i);
     h_sigmafullpidtot[i] = comp->fTPCPtDCASigmaFullPid[i]->Project3D();
     h_sigmafullpidtot[i]->SetName(name);

     comp->fTPCPtDCASigmaDay0Pid[i]->GetYaxis()->SetRange(1,4);
     sprintf(name,"h_sigmaday0pidtot_%d",i);
     h_sigmaday0pidtot[i] = comp->fTPCPtDCASigmaDay0Pid[i]->Project3D();
     h_sigmaday0pidtot[i]->SetName(name);

     // pid wise
     for(Int_t j=0; j<5; ++j)
	 {
       idx = i*5 + j;

       comp->fTPCPtDCASigmaIdealPid[i]->GetYaxis()->SetRange(1,4);
       comp->fTPCPtDCASigmaIdealPid[i]->GetZaxis()->SetRange(j+1,j+1);

       sprintf(name,"h_sigmaidealpid_%d",idx);
       h_sigmaidealpid[idx] = comp->fTPCPtDCASigmaIdealPid[i]->Project3D();
       h_sigmaidealpid[idx]->SetName(name);
	   

       comp->fTPCPtDCASigmaFullPid[i]->GetYaxis()->SetRange(1,4);
       comp->fTPCPtDCASigmaFullPid[i]->GetZaxis()->SetRange(j+1,j+1);

       sprintf(name,"h_sigmafullpid_%d",idx);
       h_sigmafullpid[idx] = comp->fTPCPtDCASigmaFullPid[i]->Project3D();
       h_sigmafullpid[idx]->SetName(name);
       
       comp->fTPCPtDCASigmaDay0Pid[i]->GetYaxis()->SetRange(1,4);
       comp->fTPCPtDCASigmaDay0Pid[i]->GetZaxis()->SetRange(j+1,j+1);

       sprintf(name,"h_sigmaday0pid_%d",idx);
       h_sigmaday0pid[idx] = comp->fTPCPtDCASigmaDay0Pid[i]->Project3D();
       h_sigmaday0pid[idx]->SetName(name);
	} 
  }

  // calculate efficiency and contamination (all pids)
  h_sigmaidealpidtot[0]->Sumw2();
  h_sigmaidealpidtot[1]->Divide(h_sigmaidealpidtot[0]);
  h_sigmaidealpidtot[2]->Sumw2();
  h_sigmaidealpidtot[3]->Divide(h_sigmaidealpidtot[2]);

  h_sigmafullpidtot[0]->Sumw2();
  h_sigmafullpidtot[1]->Divide(h_sigmafullpidtot[0]);
  h_sigmafullpidtot[2]->Sumw2();
  h_sigmafullpidtot[3]->Divide(h_sigmafullpidtot[2]);

  h_sigmaday0pidtot[0]->Sumw2();
  h_sigmaday0pidtot[1]->Divide(h_sigmaday0pidtot[0]);
  h_sigmaday0pidtot[2]->Sumw2();
  h_sigmaday0pidtot[3]->Divide(h_sigmaday0pidtot[2]);

  // calculate efficiency pid wise
  for(Int_t jdx = 0; jdx<5; jdx++)
  {
    h_sigmaidealpid[jdx]->Sumw2();
    h_sigmaidealpid[jdx+5]->Divide(h_sigmaidealpid[jdx]);

    h_sigmafullpid[jdx]->Sumw2();
    h_sigmafullpid[jdx+5]->Divide(h_sigmafullpid[jdx]);

    h_sigmaday0pid[jdx]->Sumw2();
    h_sigmaday0pid[jdx+5]->Divide(h_sigmaday0pid[jdx]);
  }

  // calculate cont. pid wise
  for(Int_t jdx = 0; jdx<5; jdx++)
  {
    h_sigmaidealpid[jdx+15]->Divide(h_sigmaidealpidtot[2]);
    h_sigmafullpid[jdx+15]->Divide(h_sigmafullpidtot[2]);
    h_sigmaday0pid[jdx+15]->Divide(h_sigmaday0pidtot[2]);
  }

  TCanvas * c = new TCanvas("Efficiency","Track efficiency");
  c->cd();
  c->Divide(1,2);

  //
  c->cd(1);
  comp->fEffTPCTanF->SetXTitle("Tan(#theta)");
  comp->fEffTPCTanF->SetYTitle("eff_{findable}");
  comp->fEffTPCTanF->SetName("EffTanFindable");
  comp->fEffTPCTanF->Draw();
  //
  c->cd(2);
  comp->fEffTPCTan->SetXTitle("Tan(#theta)");
  comp->fEffTPCTan->SetYTitle("eff_{all}");
  comp->fEffTPCTan->SetName("EffTanAll");
  comp->fEffTPCTan->Draw();

  aFolderObj->Add(comp->fEffTPCTanF);
  aFolderObj->Add(comp->fEffTPCTan);

  h_sigmaidealpidtot[1]->SetXTitle("p_{t}");
  h_sigmaidealpidtot[1]->SetYTitle("efficiency");
  h_sigmaidealpidtot[1]->SetTitle("Eff_SigmaIdeal");
  h_sigmaidealpidtot[1]->SetName("Eff_SigmaIdeal");

  h_sigmaidealpidtot[3]->SetXTitle("p_{t}");
  h_sigmaidealpidtot[3]->SetYTitle("contamination");
  h_sigmaidealpidtot[3]->SetTitle("Cont_SigmaIdeal");
  h_sigmaidealpidtot[3]->SetName("Cont_SigmaIdeal");

  aFolderObj->Add(h_sigmaidealpidtot[1]);
  aFolderObj->Add(h_sigmaidealpidtot[3]);

  h_sigmafullpidtot[1]->SetXTitle("p_{t}");
  h_sigmafullpidtot[1]->SetYTitle("efficiency");
  h_sigmafullpidtot[1]->SetTitle("Eff_SigmaFull");
  h_sigmafullpidtot[1]->SetName("Eff_SigmaFull");

  h_sigmafullpidtot[3]->SetXTitle("p_{t}");
  h_sigmafullpidtot[3]->SetYTitle("contamination");
  h_sigmafullpidtot[3]->SetTitle("Cont_SigmaFull");
  h_sigmafullpidtot[3]->SetName("Cont_SigmaFull");

  aFolderObj->Add(h_sigmafullpidtot[1]);
  aFolderObj->Add(h_sigmafullpidtot[3]);

  h_sigmaday0pidtot[1]->SetXTitle("p_{t}");
  h_sigmaday0pidtot[1]->SetYTitle("efficiency");
  h_sigmaday0pidtot[1]->SetTitle("Eff_SigmaDay0");
  h_sigmaday0pidtot[1]->SetName("Eff_SigmaDay0");

  h_sigmaday0pidtot[3]->SetXTitle("p_{t}");
  h_sigmaday0pidtot[3]->SetYTitle("contamination");
  h_sigmaday0pidtot[3]->SetTitle("Cont_SigmaDay0");
  h_sigmaday0pidtot[3]->SetName("Cont_SigmaDay0");

  aFolderObj->Add(h_sigmaday0pidtot[1]);
  aFolderObj->Add(h_sigmaday0pidtot[3]);

  for(Int_t jdx = 0; jdx<5; jdx++)
  {
    sprintf(name,"Eff_SigmaIdeal_%d",jdx);
    sprintf(name1,"Cont_SigmaIdeal_%d",jdx);


    h_sigmaidealpid[jdx+5]->SetXTitle("p_{t}");
    h_sigmaidealpid[jdx+5]->SetYTitle("efficiency");
    h_sigmaidealpid[jdx+5]->SetTitle(name);
    h_sigmaidealpid[jdx+5]->SetName(name);

    h_sigmaidealpid[jdx+15]->SetXTitle("p_{t}");
    h_sigmaidealpid[jdx+15]->SetYTitle("contamination");
    h_sigmaidealpid[jdx+15]->SetTitle(name1);
    h_sigmaidealpid[jdx+15]->SetName(name1);

	aFolderObj->Add(h_sigmaidealpid[jdx+5]);
	aFolderObj->Add(h_sigmaidealpid[jdx+15]);

    sprintf(name,"Eff_SigmaFull_%d",jdx);
    sprintf(name1,"Cont_SigmaFull_%d",jdx);

    h_sigmafullpid[jdx+5]->SetXTitle("p_{t}");
    h_sigmafullpid[jdx+5]->SetYTitle("efficiency");
    h_sigmafullpid[jdx+5]->SetTitle(name);
    h_sigmafullpid[jdx+5]->SetName(name);

    h_sigmafullpid[jdx+15]->SetXTitle("p_{t}");
    h_sigmafullpid[jdx+15]->SetYTitle("contamination");
    h_sigmafullpid[jdx+15]->SetTitle(name1);
    h_sigmafullpid[jdx+15]->SetName(name1);

	aFolderObj->Add(h_sigmafullpid[jdx+5]);
	aFolderObj->Add(h_sigmafullpid[jdx+15]);

    sprintf(name,"Eff_SigmaDay0_%d",jdx);
    sprintf(name1,"Cont_SigmaDay0_%d",jdx);

    h_sigmaday0pid[jdx+5]->SetXTitle("p_{t}");
    h_sigmaday0pid[jdx+5]->SetYTitle("efficiency");
    h_sigmaday0pid[jdx+5]->SetTitle(name);
    h_sigmaday0pid[jdx+5]->SetName(name);

    h_sigmaday0pid[jdx+15]->SetXTitle("p_{t}");
    h_sigmaday0pid[jdx+15]->SetYTitle("contamination");
    h_sigmaday0pid[jdx+15]->SetTitle(name1);
    h_sigmaday0pid[jdx+15]->SetName(name1);

	aFolderObj->Add(h_sigmaday0pid[jdx+5]);
	aFolderObj->Add(h_sigmaday0pid[jdx+15]);
  }

  // export objects to analysis folder
  fAnalysisFolder = ExportToFolder(aFolderObj);

  // delete only TObjArray
  if(aFolderObj) delete aFolderObj;
}

//_____________________________________________________________________________
TFolder* AliComparisonEff::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliComparisonEff * comp=this;
  TFolder *folder = comp->GetAnalysisFolder();

  TString name, title;
  TFolder *newFolder = 0;
  Int_t i = 0;
  Int_t size = array->GetSize();

  if(folder) { 
     // get name and title from old folder
     name = folder->GetName();  
     title = folder->GetTitle();  

	 // delete old one
     delete folder;

	 // create new one
     newFolder = CreateFolder(name.Data(),title.Data());
     newFolder->SetOwner();

	 // add objects to folder
     while(i < size) {
	   newFolder->Add(array->At(i));
	   i++;
	 }
  }

return newFolder;
}


//_____________________________________________________________________________
TFolder* AliComparisonEff::CreateFolder(TString name,TString title) { 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}
