//------------------------------------------------------------------------------
// Implementation of AliPerformanceTPC class. It keeps information from 
// comparison of reconstructed and MC particle tracks. In addtion, 
// it keeps selection cuts used during comparison. The comparison 
// information is stored in the ROOT histograms. Analysis of these 
// histograms can be done by using Analyse() class function. The result of 
// the analysis (histograms/graphs) are stored in the folder which is
// a data member of AliPerformanceTPC.
//
// Author: J.Otwinowski 04/02/2008 
// Changes by M.Knichel 15/10/2010
// Changes by J.Salzwedel 01/10/2014
//------------------------------------------------------------------------------

/*
 
  // after running comparison task, read the file, and get component
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/Macros/LoadMyLibs.C");
  LoadMyLibs();

  TFile f("Output.root");
  AliPerformanceTPC * compObj = (AliPerformanceTPC*)coutput->FindObject("AliPerformanceTPC");
 
  // analyse comparison data
  compObj->Analyse();

  // the output histograms/graphs will be stored in the folder "folderTPC" 
  compObj->GetAnalysisFolder()->ls("*");

  // user can save whole comparison object (or only folder with anlysed histograms) 
  // in the seperate output file (e.g.)
  TFile fout("Analysed_TPC.root","recreate");
  compObj->Write(); // compObj->GetAnalysisFolder()->Write();
  fout.Close();

*/

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TPostScript.h"
#include "TString.h"
#include "TUUID.h"
#include "TTree.h"
#include "TChain.h"
#include "TStopwatch.h"
#include "AliTPCPerformanceSummary.h"
#include "TSystem.h"
#include "AliPerformanceTPC.h"
#include "AliVEvent.h" 
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliESDtrack.h"
#include "AliLog.h" 
#include "AliMCEvent.h" 
#include "AliHeader.h" 
#include "AliGenEventHeader.h" 
#include "AliStack.h" 
#include "AliMCInfoCuts.h" 
#include "AliRecInfoCuts.h" 
#include "AliTracker.h" 
#include "AliTreeDraw.h" 
#include "AliTPCTransform.h" 
#include "AliTPCseed.h" 
#include "AliTPCcalibDB.h" 
#include "AliVfriendEvent.h" 
#include "AliVfriendTrack.h" 
#include "AliTPCclusterMI.h" 

using namespace std;

ClassImp(AliPerformanceTPC)

Bool_t AliPerformanceTPC::fgMergeTHnSparse = kFALSE;
Bool_t AliPerformanceTPC::fgUseMergeTHnSparse = kFALSE;


//Cluster Histograms
TH3D *h_tpc_clust_0_1_2 = 0;

//Event Histograms - Xv:Yv:Zv:mult:multP:multN:vertStatus
TH1D *h_tpc_event_recvertex_0 = 0;
TH1D *h_tpc_event_recvertex_1 = 0;
TH1D *h_tpc_event_recvertex_2 = 0;
TH1D *h_tpc_event_recvertex_3 = 0;
TH1D *h_tpc_event_recvertex_4 = 0;
TH1D *h_tpc_event_recvertex_5 = 0;
TH1D *h_tpc_event_6 = 0;

//Track Histograms - nTPCClust:chi2PerTPCClust:nTPCClustFindRatio:DCAr:DCAz:eta:phi:pt:charge:vertStatus

TH3D* h_tpc_track_pos_recvertex_2_5_6 = 0;
TH3D* h_tpc_track_neg_recvertex_2_5_6 = 0;

TH2D *h_tpc_track_all_recvertex_5_8 = 0;
TH3D *h_tpc_track_all_recvertex_0_5_7 = 0;
TH3D *h_tpc_track_pos_recvertex_0_5_7 = 0;
TH3D *h_tpc_track_neg_recvertex_0_5_7 = 0;

TH3D *h_tpc_track_all_recvertex_1_5_7 = 0;
TH3D *h_tpc_track_all_recvertex_2_5_7 = 0;
TH3D *h_tpc_track_all_recvertex_3_5_7 = 0;
TH3D *h_tpc_track_pos_recvertex_3_5_7 = 0;
TH3D *h_tpc_track_neg_recvertex_3_5_7 = 0;

TH3D *h_tpc_track_all_recvertex_4_5_7 = 0;
TH3D *h_tpc_track_pos_recvertex_4_5_7 = 0;
TH3D *h_tpc_track_neg_recvertex_4_5_7 = 0;
TH3D *h_tpc_track_pos_recvertex_3_5_6 = 0;
TH3D *h_tpc_track_pos_recvertex_4_5_6 = 0;
TH3D *h_tpc_track_neg_recvertex_3_5_6 = 0;
TH3D *h_tpc_track_neg_recvertex_4_5_6 = 0;

//_____________________________________________________________________________
/*
AliPerformanceTPC::AliPerformanceTPC():
  AliPerformanceObject("AliPerformanceTPC"),
  fTPCClustHisto(0),
  fTPCEventHisto(0),
  fTPCTrackHisto(0),
  fFolderObj(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  
  fUseHLT(kFALSE)

{
  Init();
}
*/

//_____________________________________________________________________________
AliPerformanceTPC::AliPerformanceTPC(const Char_t* name, const Char_t* title,Int_t analysisMode,Bool_t hptGenerator, Int_t run, Bool_t highMult, Bool_t useSparse):
  AliPerformanceObject(name,title,run,highMult),
  fTPCClustHisto(0),
  fTPCEventHisto(0),
  fTPCTrackHisto(0),
  fFolderObj(0),

  // Cuts 
  fCutsRC(0),  
  fCutsMC(0),  

  // histogram folder 
  fAnalysisFolder(0),
  fUseHLT(kFALSE)

{

// named constructor
  // 
  SetAnalysisMode(analysisMode);
  SetHptGenerator(hptGenerator);
  fUseSparse = useSparse;
  rtime = 0;
  revent = 0;
  Init();
}


//_____________________________________________________________________________
AliPerformanceTPC::~AliPerformanceTPC()
{
  // destructor
   
    if(fTPCClustHisto) delete fTPCClustHisto; fTPCClustHisto=0;
    if(fTPCEventHisto) delete fTPCEventHisto; fTPCEventHisto=0;
    if(fTPCTrackHisto) delete fTPCTrackHisto; fTPCTrackHisto=0;
    if(fAnalysisFolder) delete fAnalysisFolder; fAnalysisFolder=0;
    if(fFolderObj) delete fFolderObj; fFolderObj=0;
    
    if(h_tpc_clust_0_1_2) delete h_tpc_clust_0_1_2; h_tpc_clust_0_1_2=0;
    if(h_tpc_event_recvertex_0) delete h_tpc_event_recvertex_0; h_tpc_event_recvertex_0=0;
    if(h_tpc_event_recvertex_1) delete h_tpc_event_recvertex_1; h_tpc_event_recvertex_1=0;
    if(h_tpc_event_recvertex_2) delete h_tpc_event_recvertex_2; h_tpc_event_recvertex_2=0;
    if(h_tpc_event_recvertex_3) delete h_tpc_event_recvertex_3; h_tpc_event_recvertex_3=0;
    if(h_tpc_event_recvertex_4) delete h_tpc_event_recvertex_4; h_tpc_event_recvertex_4=0;
    if(h_tpc_event_recvertex_5) delete h_tpc_event_recvertex_5; h_tpc_event_recvertex_5=0;
    if(h_tpc_event_6) delete h_tpc_event_6; h_tpc_event_6=0;
    if(h_tpc_track_all_recvertex_5_8) delete h_tpc_track_all_recvertex_5_8; h_tpc_track_all_recvertex_5_8=0;
    if(h_tpc_track_all_recvertex_0_5_7) delete h_tpc_track_all_recvertex_0_5_7; h_tpc_track_all_recvertex_0_5_7=0;
    if(h_tpc_track_pos_recvertex_0_5_7) delete h_tpc_track_pos_recvertex_0_5_7; h_tpc_track_pos_recvertex_0_5_7=0;
    if(h_tpc_track_neg_recvertex_0_5_7) delete h_tpc_track_neg_recvertex_0_5_7; h_tpc_track_neg_recvertex_0_5_7=0;
    if(h_tpc_track_all_recvertex_1_5_7) delete h_tpc_track_all_recvertex_1_5_7; h_tpc_track_all_recvertex_1_5_7=0;
    if(h_tpc_track_all_recvertex_2_5_7) delete h_tpc_track_all_recvertex_2_5_7; h_tpc_track_all_recvertex_2_5_7=0;
    if(h_tpc_track_all_recvertex_3_5_7) delete h_tpc_track_all_recvertex_3_5_7; h_tpc_track_all_recvertex_3_5_7=0;
    if(h_tpc_track_pos_recvertex_3_5_7) delete h_tpc_track_pos_recvertex_3_5_7; h_tpc_track_pos_recvertex_3_5_7=0;
    if(h_tpc_track_neg_recvertex_3_5_7) delete h_tpc_track_neg_recvertex_3_5_7; h_tpc_track_neg_recvertex_3_5_7=0;
    if(h_tpc_track_all_recvertex_4_5_7) delete h_tpc_track_all_recvertex_4_5_7; h_tpc_track_all_recvertex_4_5_7=0;
    if(h_tpc_track_pos_recvertex_4_5_7) delete h_tpc_track_pos_recvertex_4_5_7; h_tpc_track_pos_recvertex_4_5_7=0;
    if(h_tpc_track_neg_recvertex_4_5_7) delete h_tpc_track_neg_recvertex_4_5_7; h_tpc_track_neg_recvertex_4_5_7=0;
    if(h_tpc_track_pos_recvertex_3_5_6) delete h_tpc_track_pos_recvertex_3_5_6; h_tpc_track_pos_recvertex_3_5_6=0;
    if(h_tpc_track_pos_recvertex_4_5_6) delete h_tpc_track_pos_recvertex_4_5_6; h_tpc_track_pos_recvertex_4_5_6=0;
    if(h_tpc_track_neg_recvertex_3_5_6) delete h_tpc_track_neg_recvertex_3_5_6; h_tpc_track_neg_recvertex_3_5_6=0;
    if(h_tpc_track_neg_recvertex_4_5_6) delete h_tpc_track_neg_recvertex_4_5_6; h_tpc_track_neg_recvertex_4_5_6=0;
    if(h_tpc_track_pos_recvertex_2_5_6) delete h_tpc_track_pos_recvertex_2_5_6; h_tpc_track_pos_recvertex_2_5_6 = 0;
    if(h_tpc_track_neg_recvertex_2_5_6) delete h_tpc_track_neg_recvertex_2_5_6; h_tpc_track_neg_recvertex_2_5_6 = 0;
}


//_____________________________________________________________________________
void AliPerformanceTPC::Init()
{
  //
  // histogram bining
  //
    fAnalysisFolder = CreateFolder("folderTPC","Analysis Resolution Folder");
    if(!fUseSparse) fFolderObj = new TObjArray;
 
    
  // set pt bins
  Int_t nPtBins = 50;
  Double_t ptMin = 1.e-2, ptMax = 20.;

  Double_t *binsPt = 0;

  if (IsHptGenerator())  { 
        ptMax = 100.;
  } 
   binsPt = CreateLogAxis(nPtBins,ptMin,ptMax);
 TString hname, htitle;


  /*
  const Int_t  nCOverPtBins = 80;
  Double_t coverptMin = -10, coverptMax = 10;
  Double_t *binsCOverPtP = 0;
  Double_t *binsCOverPt = new Double_t[nCOverPtBins+1];
  binsCOverPtP = CreateLogAxis(nCOverPtBins/2,0.04,coverptMax-0.04);
  for(Int_t i=0; i < nCOverPtBins/2; i++){
    binsCOverPt[nCOverPtBins - i] = binsCOverPtP[nCOverPtBins/2  - i];
    binsCOverPt[i] = 0 - binsCOverPtP[nCOverPtBins/2  - i];
 }
 */

  /*
  Int_t nPtBins = 31;
  Double_t binsPt[32] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3.,3.5,4.,5.,6.,8.,10.};
  Double_t ptMin = 0., ptMax = 10.; 

  if(IsHptGenerator() == kTRUE) {
    nPtBins = 100;
    ptMin = 0.; ptMax = 100.; 
  }
  */
  //

  //
  //padRow:phi:TPCSide
  Int_t binsTPCClustHisto[3] =   {160,  144,  2};
  Double_t minTPCClustHisto[3] = {0.,   0.,   0.};
  Double_t maxTPCClustHisto[3] = {160., 2.*TMath::Pi(), 2.};

    if(fUseSparse){
        fTPCClustHisto = new THnSparseF("fTPCClustHisto","padRow:phi:TPCSide",3,binsTPCClustHisto,minTPCClustHisto,maxTPCClustHisto);
        fTPCClustHisto->GetAxis(0)->SetTitle("padRow");
        fTPCClustHisto->GetAxis(1)->SetTitle("phi (rad)");
        fTPCClustHisto->GetAxis(2)->SetTitle("TPCSide");
        fTPCClustHisto->Sumw2();
    }
    else{
        h_tpc_clust_0_1_2 = new TH3D("h_tpc_clust_0_1_2","padRow:phi:TPCSide",binsTPCClustHisto[0],minTPCClustHisto[0],maxTPCClustHisto[0],binsTPCClustHisto[1],minTPCClustHisto[1],maxTPCClustHisto[1],binsTPCClustHisto[2],minTPCClustHisto[2],maxTPCClustHisto[2]);
        h_tpc_clust_0_1_2->GetXaxis()->SetTitle("padRow");
        h_tpc_clust_0_1_2->GetYaxis()->SetTitle("phi (rad)");
        h_tpc_clust_0_1_2->GetZaxis()->SetTitle("TPCSide");
    }
    
  Float_t scaleVxy = 1.0;
  if(fAnalysisMode !=0) scaleVxy = 0.1; 

  Int_t maxMult;
  if (fHighMultiplicity) { maxMult = 4001; scaleVxy = 0.1;} else { maxMult = 151; }
  // Xv:Yv:Zv:mult:multP:multN:vertStatus
  Int_t binsTPCEventHisto[7]=  {100,  100,   100,  maxMult,  maxMult,  maxMult, 2   };
  Double_t minTPCEventHisto[7]={-10.*scaleVxy, -10.*scaleVxy, -30., -0.5,  -0.5,  -0.5, -0.5  };
  Double_t maxTPCEventHisto[7]={ 10.*scaleVxy,  10.*scaleVxy,  30.,  maxMult-0.5,  maxMult-0.5, maxMult-0.5, 1.5 };

  if(fUseSparse){
      fTPCEventHisto = new THnSparseF("fTPCEventHisto","Xv:Yv:Zv:mult:multP:multN:vertStatus",7,binsTPCEventHisto,minTPCEventHisto,maxTPCEventHisto);
      fTPCEventHisto->GetAxis(0)->SetTitle("Xv (cm)");
      fTPCEventHisto->GetAxis(1)->SetTitle("Yv (cm)");
      fTPCEventHisto->GetAxis(2)->SetTitle("Zv (cm)");
      fTPCEventHisto->GetAxis(3)->SetTitle("mult");
      fTPCEventHisto->GetAxis(4)->SetTitle("multP");
      fTPCEventHisto->GetAxis(5)->SetTitle("multN");
      fTPCEventHisto->GetAxis(6)->SetTitle("vertStatus");
  }
  else{
      
      h_tpc_event_recvertex_0 = new TH1D("h_tpc_event_recvertex_0","Xv",binsTPCEventHisto[0],minTPCEventHisto[0],maxTPCEventHisto[0]);
      h_tpc_event_recvertex_1 = new TH1D("h_tpc_event_recvertex_1","Yv",binsTPCEventHisto[1],minTPCEventHisto[1],maxTPCEventHisto[1]);
      h_tpc_event_recvertex_2 = new TH1D("h_tpc_event_recvertex_2","Zv",binsTPCEventHisto[2],minTPCEventHisto[2],maxTPCEventHisto[2]);
      h_tpc_event_recvertex_3 = new TH1D("h_tpc_event_recvertex_3","mult",binsTPCEventHisto[3],minTPCEventHisto[3],maxTPCEventHisto[3]);
      h_tpc_event_recvertex_4 = new TH1D("h_tpc_event_recvertex_4","multP",binsTPCEventHisto[4],minTPCEventHisto[4],maxTPCEventHisto[4]);
      h_tpc_event_recvertex_5 = new TH1D("h_tpc_event_recvertex_5","multN",binsTPCEventHisto[5],minTPCEventHisto[5],maxTPCEventHisto[5]);
      h_tpc_event_recvertex_0->GetXaxis()->SetTitle("Xv (cm)");
      h_tpc_event_recvertex_1->GetXaxis()->SetTitle("Yv (cm)");
      h_tpc_event_recvertex_2->GetXaxis()->SetTitle("Zv (cm)");
      h_tpc_event_recvertex_3->GetXaxis()->SetTitle("mult");
      h_tpc_event_recvertex_4->GetXaxis()->SetTitle("multP");
      h_tpc_event_recvertex_5->GetXaxis()->SetTitle("multN");
      h_tpc_event_6 = new TH1D("h_tpc_event_6","vertStatus",binsTPCEventHisto[6],minTPCEventHisto[6],maxTPCEventHisto[6]);
      h_tpc_event_6->SetXTitle("vertStatus");

      fFolderObj->Add(h_tpc_event_recvertex_0);
      fFolderObj->Add(h_tpc_event_recvertex_1);
      fFolderObj->Add(h_tpc_event_recvertex_2);
      fFolderObj->Add(h_tpc_event_recvertex_3);
      fFolderObj->Add(h_tpc_event_recvertex_4);
      fFolderObj->Add(h_tpc_event_recvertex_5);
      fFolderObj->Add(h_tpc_event_6);
      //cout<<"Past make histograms "<<endl;
      //fTPCEventTree = new TNtuple("fTPCEventTree","fTPCEventTree","Xv:Yv:Zv:mult:multP:multN:vertStatus");
  }

  Float_t scaleDCA = 1.0;
  if(fAnalysisMode !=0) scaleDCA = 0.1; 
  // nTPCClust:chi2PerTPCClust:nTPCClustFindRatio:DCAr:DCAz:eta:phi:pt:charge:vertStatus
   //Int_t binsTPCTrackHisto[10]=  { 160,  20,  60,  30, 30,  30,   144,             nPtBins,   nCOverPtBins, 2 };
   //Double_t minTPCTrackHisto[10]={ 0.,   0.,  0., -3*scaleDCA, -3.*scaleDCA, -1.5, 0.,             ptMin,   coverptMin, -0.5 };
   //Double_t maxTPCTrackHisto[10]={ 160., 5., 1.2, 3*scaleDCA,  3.*scaleDCA,  1.5, 2.*TMath::Pi(), ptMax,    coverptMax,  1.5 };
   Int_t binsTPCTrackHisto[10]=  { 160,  20,  60,  30, 30,  30,   144,             nPtBins,   3, 2 };
   Double_t minTPCTrackHisto[10]={ 0.,   0.,  0., -3*scaleDCA, -3.*scaleDCA, -1.5, 0.,             ptMin,  -1.5, -0.5 };
   Double_t maxTPCTrackHisto[10]={ 160., 5., 1.2, 3*scaleDCA,  3.*scaleDCA,  1.5, 2.*TMath::Pi(), ptMax,    1.5,  1.5 };
  
    if(fUseSparse){
        fTPCTrackHisto = new THnSparseF("fTPCTrackHisto","nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge:vertStatus",10,binsTPCTrackHisto,minTPCTrackHisto,maxTPCTrackHisto);
        fTPCTrackHisto->SetBinEdges(7,binsPt);
        fTPCTrackHisto->GetAxis(0)->SetTitle("nClust");
        fTPCTrackHisto->GetAxis(1)->SetTitle("chi2PerClust");
        fTPCTrackHisto->GetAxis(2)->SetTitle("nClust/nFindableClust");
        fTPCTrackHisto->GetAxis(3)->SetTitle("DCAr (cm)");
        fTPCTrackHisto->GetAxis(4)->SetTitle("DCAz (cm)");
        fTPCTrackHisto->GetAxis(5)->SetTitle("#eta");
        fTPCTrackHisto->GetAxis(6)->SetTitle("#phi (rad)");
        fTPCTrackHisto->GetAxis(7)->SetTitle("p_{T} (GeV/c)");
        fTPCTrackHisto->GetAxis(8)->SetTitle("charge");
        fTPCTrackHisto->GetAxis(9)->SetTitle("vertStatus");
    }
    else{

        h_tpc_track_all_recvertex_5_8 = new TH2D("h_tpc_track_all_recvertex_5_8","",binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[8],minTPCTrackHisto[8],maxTPCTrackHisto[8]);

        h_tpc_track_all_recvertex_0_5_7 = new TH3D("h_tpc_track_all_recvertex_0_5_7","",binsTPCTrackHisto[0],minTPCTrackHisto[0],maxTPCTrackHisto[0],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_pos_recvertex_0_5_7 = new TH3D("h_tpc_track_pos_recvertex_0_5_7","",binsTPCTrackHisto[0],minTPCTrackHisto[0],maxTPCTrackHisto[0],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_neg_recvertex_0_5_7 = new TH3D("h_tpc_track_neg_recvertex_0_5_7","",binsTPCTrackHisto[0],minTPCTrackHisto[0],maxTPCTrackHisto[0],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        
        h_tpc_track_all_recvertex_1_5_7 = new TH3D("h_tpc_track_all_recvertex_1_5_7","",binsTPCTrackHisto[1],minTPCTrackHisto[1],maxTPCTrackHisto[1],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_all_recvertex_2_5_7 = new TH3D("h_tpc_track_all_recvertex_2_5_7","",binsTPCTrackHisto[2],minTPCTrackHisto[2],maxTPCTrackHisto[2],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        
        h_tpc_track_all_recvertex_3_5_7 = new TH3D("h_tpc_track_all_recvertex_3_5_7","",binsTPCTrackHisto[3],minTPCTrackHisto[3],maxTPCTrackHisto[3],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_pos_recvertex_3_5_7 = new TH3D("h_tpc_track_pos_recvertex_3_5_7","",binsTPCTrackHisto[3],minTPCTrackHisto[3],maxTPCTrackHisto[3],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_neg_recvertex_3_5_7 = new TH3D("h_tpc_track_neg_recvertex_3_5_7","",binsTPCTrackHisto[3],minTPCTrackHisto[3],maxTPCTrackHisto[3],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        
        h_tpc_track_all_recvertex_4_5_7 = new TH3D("h_tpc_track_all_recvertex_4_5_7","",binsTPCTrackHisto[4],minTPCTrackHisto[4],maxTPCTrackHisto[4],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_pos_recvertex_4_5_7 = new TH3D("h_tpc_track_pos_recvertex_4_5_7","",binsTPCTrackHisto[4],minTPCTrackHisto[4],maxTPCTrackHisto[4],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        h_tpc_track_neg_recvertex_4_5_7 = new TH3D("h_tpc_track_neg_recvertex_4_5_7","",binsTPCTrackHisto[4],minTPCTrackHisto[4],maxTPCTrackHisto[4],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[7],minTPCTrackHisto[7],maxTPCTrackHisto[7]);
        
        h_tpc_track_pos_recvertex_3_5_6 = new TH3D("h_tpc_track_pos_recvertex_3_5_6","",binsTPCTrackHisto[3],minTPCTrackHisto[3],maxTPCTrackHisto[3],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);
        h_tpc_track_pos_recvertex_4_5_6 = new TH3D("h_tpc_track_pos_recvertex_4_5_6","",binsTPCTrackHisto[4],minTPCTrackHisto[4],maxTPCTrackHisto[4],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);
        h_tpc_track_neg_recvertex_3_5_6 = new TH3D("h_tpc_track_neg_recvertex_3_5_6","",binsTPCTrackHisto[3],minTPCTrackHisto[3],maxTPCTrackHisto[3],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);
        h_tpc_track_neg_recvertex_4_5_6 = new TH3D("h_tpc_track_neg_recvertex_4_5_6","",binsTPCTrackHisto[4],minTPCTrackHisto[4],maxTPCTrackHisto[4],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);
        h_tpc_track_pos_recvertex_2_5_6 = new TH3D("h_tpc_track_pos_recvertex_2_5_6","",binsTPCTrackHisto[2],minTPCTrackHisto[2],maxTPCTrackHisto[2],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);
        h_tpc_track_neg_recvertex_2_5_6 = new TH3D("h_tpc_track_neg_recvertex_2_5_6","",binsTPCTrackHisto[2],minTPCTrackHisto[2],maxTPCTrackHisto[2],binsTPCTrackHisto[5],minTPCTrackHisto[5],maxTPCTrackHisto[5],binsTPCTrackHisto[6],minTPCTrackHisto[6],maxTPCTrackHisto[6]);

        
        fFolderObj->Add(h_tpc_track_all_recvertex_5_8);
        fFolderObj->Add(h_tpc_track_all_recvertex_0_5_7);
        fFolderObj->Add(h_tpc_track_pos_recvertex_0_5_7);
        fFolderObj->Add(h_tpc_track_neg_recvertex_0_5_7);
        fFolderObj->Add(h_tpc_track_all_recvertex_1_5_7);
        fFolderObj->Add(h_tpc_track_all_recvertex_2_5_7);
        fFolderObj->Add(h_tpc_track_all_recvertex_3_5_7);
        fFolderObj->Add(h_tpc_track_pos_recvertex_3_5_7);
        fFolderObj->Add(h_tpc_track_neg_recvertex_3_5_7);
        fFolderObj->Add(h_tpc_track_all_recvertex_4_5_7);
        fFolderObj->Add(h_tpc_track_pos_recvertex_4_5_7);
        fFolderObj->Add(h_tpc_track_neg_recvertex_4_5_7);
        fFolderObj->Add(h_tpc_track_pos_recvertex_3_5_6);
        fFolderObj->Add(h_tpc_track_pos_recvertex_4_5_6);
        fFolderObj->Add(h_tpc_track_neg_recvertex_3_5_6);
        fFolderObj->Add(h_tpc_track_neg_recvertex_4_5_6);
        fFolderObj->Add(h_tpc_track_pos_recvertex_2_5_6);
        fFolderObj->Add(h_tpc_track_neg_recvertex_2_5_6);
    }

  // Init cuts 
  if(!fCutsMC) {
    AliDebug(AliLog::kError, "ERROR: Cannot find AliMCInfoCuts object");
  }
  if(!fCutsRC) {
    AliDebug(AliLog::kError, "ERROR: Cannot find AliRecInfoCuts object"); 
  }

  // init folder

  //delete []binsCOverPt;
  
   // save merge status in object
  fMergeTHnSparseObj = fgMergeTHnSparse;
  
}


//_____________________________________________________________________________
void AliPerformanceTPC::ProcessTPC(AliStack* const stack, AliVTrack *const vTrack, AliVEvent *const vEvent, Bool_t vertStatus)
{
//
// fill TPC QA info
//
  if(!vEvent) return;
  if(!vTrack) return;
    
  if(IsUseTOFBunchCrossing()){
    if(vTrack->GetTOFBunchCrossing(vEvent->GetMagneticField())!=0){
      return;
    }
  }
  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    //    Bool_t isOK = esdTrack->RelateToVertexTPCBxByBz(vtxESD, b, kVeryBig);
    Bool_t isOK=kFALSE;
    if(fabs(b[2])>0.000001)
      isOK = vTrack->RelateToVVertexTPCBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;
  }   

  // Fill TPC only resolution comparison information 
  const AliExternalTrackParam *etpTrack = vTrack->GetTPCInnerParam();
  if(!etpTrack) return;

  Float_t dca[2] = {0.,0.};
  Float_t cov[3] = {0.,0.,0.}; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParametersTPC(dca,cov);

  Float_t q = vTrack->Charge();
  Float_t pt = etpTrack->Pt();
  Float_t eta = etpTrack->Eta();
  Float_t phi = etpTrack->Phi();
  UShort_t nClust = vTrack->GetTPCNcls();
  Int_t nFindableClust = vTrack->GetTPCNclsF();

  Float_t chi2PerCluster = 0.;
  if(nClust>0.) chi2PerCluster = vTrack->GetTPCchi2()/Float_t(nClust);

  Float_t clustPerFindClust = 0.;
  if(nFindableClust>0.) clustPerFindClust = Float_t(nClust)/nFindableClust;
  
  // Float_t qpt = 0;
  // if( fabs(pt)>0 ) qpt = q/fabs(pt);

  // filter out noise tracks
  if(vTrack->GetTPCsignal() < 5) return;

  //
  // select primaries
  //
  Double_t dcaToVertex = -1;
  if( fCutsRC->GetDCAToVertex2D() ) 
  {
    dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC->GetMaxDCAToVertexXY()/fCutsRC->GetMaxDCAToVertexXY() + dca[1]*dca[1]/fCutsRC->GetMaxDCAToVertexZ()/fCutsRC->GetMaxDCAToVertexZ()); 
  }
  if(fCutsRC->GetDCAToVertex2D() && dcaToVertex > 1) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC->GetMaxDCAToVertexXY()) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC->GetMaxDCAToVertexZ()) return;
  if(nClust < fCutsRC->GetMinNClustersTPC()) return;

  //Double_t vTPCTrackHisto[10] = {nClust,chi2PerCluster,clustPerFindClust,dca[0],dca[1],eta,phi,pt,qpt,vertStatus};
  Double_t vTPCTrackHisto[10] = {static_cast<Double_t>(nClust),static_cast<Double_t>(chi2PerCluster),static_cast<Double_t>(clustPerFindClust),static_cast<Double_t>(dca[0]),static_cast<Double_t>(dca[1]),static_cast<Double_t>(eta),static_cast<Double_t>(phi),static_cast<Double_t>(pt),static_cast<Double_t>(q),static_cast<Double_t>(vertStatus)};
  
  //nClust:chi2PerClust:nClust/nFindableClust:DCAr:DCAz:eta:phi:pt:charge:vertStatus
    
  if(fUseSparse) fTPCTrackHisto->Fill(vTPCTrackHisto);
    else FillTrackHistogram(vTPCTrackHisto);
    //
  // Fill rec vs MC information
  //
  if(!stack) return;

}


//_____________________________________________________________________________
void AliPerformanceTPC::ProcessTPCITS(AliStack* const stack, AliVTrack *const vTrack, AliVEvent* const vEvent, Bool_t vertStatus)
{
  // Fill comparison information (TPC+ITS) 
  if(!vTrack) return;
  if(!vEvent) return;

  if( IsUseTrackVertex() ) 
  { 
    // Relate TPC inner params to prim. vertex
    const AliVVertex *vVertex = vEvent->GetPrimaryVertexTracks();
    Double_t x[3]; vTrack->GetXYZ(x);
    Double_t b[3]; AliTracker::GetBxByBz(x,b);
    Bool_t isOK = vTrack->RelateToVVertexBxByBz(vVertex, b, kVeryBig);
    if(!isOK) return;

    /*
      // JMT -- recaluclate DCA for HLT if not present
      if ( dca[0] == 0. && dca[1] == 0. ) {
        track->GetDZ( vtxESD->GetX(), vtxESD->GetY(), vtxESD->GetZ(), esdEvent->GetMagneticField(), dca );
      }
    */
  }

  Float_t dca[2] = {0.,0.}; // dca_xy, dca_z
  Float_t cov[3] = {0.,0.,0.}; // sigma_xy, sigma_xy_z, sigma_z
  vTrack->GetImpactParameters(dca,cov);
  if ((vTrack->GetStatus()&AliVTrack::kITSrefit)==0) return; // ITS refit
  if ((vTrack->GetStatus()&AliVTrack::kTPCrefit)==0) return; // TPC refit
  if ((vTrack->HasPointOnITSLayer(0)==kFALSE)&&(vTrack->HasPointOnITSLayer(1)==kFALSE)) return; // at least one SPD
  //if (vTrack->GetITSclusters(0)<fCutsRC->GetMinNClustersITS()) return;  // min. nb. ITS clusters

  Float_t chi2PerCluster = 0.;
  Float_t q = vTrack->Charge();
  Float_t pt = vTrack->Pt();
  Float_t eta = vTrack->Eta();
  Float_t phi = vTrack->Phi();
  UShort_t nClust = vTrack->GetTPCNcls();
  Int_t nFindableClust = vTrack->GetTPCNclsF();
  if(nClust>0.) chi2PerCluster = vTrack->GetTPCchi2()/Float_t(nClust);
  
  Float_t clustPerFindClust = 0.;
  if(nFindableClust>0.) clustPerFindClust = Float_t(nClust)/nFindableClust;
  // Float_t qpt = 0;
  // if( fabs(pt)>0 ) qpt = q/fabs(pt);

  //
  // select primaries
  //
  Double_t dcaToVertex = -1;
  if( fCutsRC->GetDCAToVertex2D() ) 
  {
    dcaToVertex = TMath::Sqrt(dca[0]*dca[0]/fCutsRC->GetMaxDCAToVertexXY()/fCutsRC->GetMaxDCAToVertexXY() + dca[1]*dca[1]/fCutsRC->GetMaxDCAToVertexZ()/fCutsRC->GetMaxDCAToVertexZ()); 
  }
  if(fCutsRC->GetDCAToVertex2D() && dcaToVertex > 1) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[0]) > fCutsRC->GetMaxDCAToVertexXY()) return;
  if(!fCutsRC->GetDCAToVertex2D() && TMath::Abs(dca[1]) > fCutsRC->GetMaxDCAToVertexZ()) return;
  if(nClust < fCutsRC->GetMinNClustersTPC()) return;

  Double_t vTPCTrackHisto[10] = {static_cast<Double_t>(nClust),static_cast<Double_t>(chi2PerCluster),static_cast<Double_t>(clustPerFindClust),static_cast<Double_t>(dca[0]),static_cast<Double_t>(dca[1]),static_cast<Double_t>(eta),static_cast<Double_t>(phi),static_cast<Double_t>(pt),static_cast<Double_t>(q),static_cast<Double_t>(vertStatus)};
    if(fUseSparse) fTPCTrackHisto->Fill(vTPCTrackHisto);
    else FillTrackHistogram(vTPCTrackHisto);
  //
  // Fill rec vs MC information
  //
  if(!stack) return;
}


//_____________________________________________________________________________
void AliPerformanceTPC::ProcessConstrained(AliStack* const /*stack*/, AliVTrack *const /*vTrack*/, AliVEvent* const /*vEvent*/)
{
  // Fill comparison information (constarained parameters) 
  AliDebug(AliLog::kWarning, "Warning: Not implemented");
}


//_____________________________________________________________________________
void AliPerformanceTPC::Exec(AliMCEvent* const mcEvent, AliVEvent *const vEvent, AliVfriendEvent *const vfriendEvent, const Bool_t bUseMC, const Bool_t bUseVfriend)
{
  // Process comparison information 
  //
    
    if(!vEvent)
  {
    Error("Exec","vEvent not available");
    return;
  }

    TStopwatch Watch;
    Watch.Start();

    
  AliHeader* header = 0;
  AliGenEventHeader* genHeader = 0;
  AliStack* stack = 0;
  TArrayF vtxMC(3);
  
  if(bUseMC)
  {
    if(!mcEvent) {
      Error("Exec","mcEvent not available");
      return;
    }
    // get MC event header
    header = mcEvent->Header();
    if (!header) {
      Error("Exec","Header not available");
      return;
    }
    // MC particle stack
    stack = mcEvent->Stack();
    if (!stack) {
      Error("Exec","Stack not available");
      return;
    }
    // get MC vertex
    genHeader = header->GenEventHeader();
    if (!genHeader) {
      Error("Exec","Could not retrieve genHeader from Header");
      return;
    }
    genHeader->PrimaryVertex(vtxMC);
  } 
  
  // check trigger
  
  if(!bUseMC && GetTriggerClass()) {
    Bool_t isEventTriggered = vEvent->IsTriggerClassFired(GetTriggerClass());
    if(!isEventTriggered) {
      printf("ERROR: Could not determine trigger class");
      return;
    }
  }


  // get TPC event vertex
  const AliVVertex *vVertex = NULL;
  if(fUseTrackVertex) {
    vVertex = vEvent->GetPrimaryVertexTracks();
  } else {
    vVertex = vEvent->GetPrimaryVertexTPC();
  }

  if(!vVertex) {
    printf("ERROR: Could not determine primary vertex");
    return;
  }

  //  events with rec. vertex
  Int_t mult=0; Int_t multP=0; Int_t multN=0;
  
  // store vertex status
  Bool_t vertStatus = vVertex->GetStatus(); 
  //  Process events
  for (Int_t iTrack = 0; iTrack < vEvent->GetNumberOfTracks(); iTrack++) 
  {
    
    AliVParticle *particle = vEvent->GetTrack(iTrack);
    if(!particle) continue;
    AliVTrack *vTrack = dynamic_cast<AliVTrack*>(particle);
    if(!vTrack) continue;

    // if not fUseKinkDaughters don't use tracks with kink index > 0
    if(!fUseKinkDaughters && vTrack->GetKinkIndex(0) > 0) continue;
    
    if(bUseVfriend && vfriendEvent && vfriendEvent->TestSkipBit()==kFALSE && iTrack<vfriendEvent->GetNumberOfTracks()) 
    {
      const AliVfriendTrack *friendTrack=vfriendEvent->GetTrack(iTrack);
      if(friendTrack) 
      {
        //
        TObject *calibObject=0;
        AliTPCseed *seed=0;
        for (Int_t j=0;(calibObject=friendTrack->GetCalibObject(j));++j) {
	    if ((seed=dynamic_cast<AliTPCseed*>(calibObject))) {
	    break;
	  }
        }

        // 
	for (Int_t irow=0;irow<159;irow++) {
	if(!seed) continue;
	  
	  AliTPCclusterMI *cluster=seed->GetClusterPointer(irow);
	  if (!cluster) continue;

	     Float_t gclf[3];
	     cluster->GetGlobalXYZ(gclf);

	     //Double_t x[3]={cluster->GetRow(),cluster->GetPad(),cluster->GetTimeBin()};
	     //Int_t i[1]={cluster->GetDetector()};
             //transform->Transform(x,i,0,1);
	     //printf("gx %f gy  %f  gz %f \n", cluster->GetX(), cluster->GetY(),cluster->GetZ());
	     //printf("gclf[0] %f gclf[1]  %f  gclf[2] %f \n", gclf[0], gclf[1],  gclf[2]);
     
             Int_t TPCside; 
	     if(gclf[2]>0.) TPCside=0; // A side 
	     else TPCside=1;

	     //
             //Double_t vTPCClust1[3] = { gclf[0], gclf[1],  TPCside };
             //fTPCClustHisto1->Fill(vTPCClust1);

             //  
	     Double_t phi = TMath::ATan2(gclf[1],gclf[0]);
	     if(phi < 0) phi += 2.*TMath::Pi();
	    
	  //Float_t pad = cluster->GetPad();
	  //Int_t detector = cluster->GetDetector();
	  //Double_t vTPCClust[6] = { irow, phi, TPCside, pad, detector, gclf[2] };
	  Double_t vTPCClust[3] = { static_cast<Double_t>(irow), phi, static_cast<Double_t>(TPCside) };
        if(fUseSparse) fTPCClustHisto->Fill(vTPCClust);
        else{
            h_tpc_clust_0_1_2->Fill(vTPCClust[0],vTPCClust[1],vTPCClust[2]);
        }
    } //end if(bUseVfriend && vfriendEvent && ...)
      }}
    if(GetAnalysisMode() == 0) ProcessTPC(stack,vTrack,vEvent,vertStatus);
    else if(GetAnalysisMode() == 1) ProcessTPCITS(stack,vTrack,vEvent,vertStatus);
    else if(GetAnalysisMode() == 2) ProcessConstrained(stack,vTrack,vEvent);
    else {
      printf("ERROR: AnalysisMode %d \n",fAnalysisMode);
      return;
    }

    // TPC only
    if(!fUseHLT && (GetAnalysisMode() == 0) ){
      AliESDtrack *tpcTrack = AliESDtrackCuts::GetTPCOnlyTrackFromVEvent(vEvent,iTrack);
      if(!tpcTrack) continue;
      // track selection
      if( fCutsRC->AcceptVTrack(tpcTrack) ) { 
	mult++;
	if(tpcTrack->Charge()>0.) multP++;
	if(tpcTrack->Charge()<0.) multN++;
      }
      if(tpcTrack) delete tpcTrack;
    } //end if(!fUseHLT)
    else {
      if( fCutsRC->AcceptVTrack(vTrack) ) { 
	mult++;
	if(vTrack->Charge()>0.) multP++;
	if(vTrack->Charge()<0.) multN++;
      }
    }
  } //end iTrack iteration

    Double_t vtxPosition[3]= {0.,0.,0.};
  vVertex->GetXYZ(vtxPosition);
  Double_t vTPCEvent[7] = {vtxPosition[0],vtxPosition[1],vtxPosition[2],static_cast<Double_t>(mult),static_cast<Double_t>(multP),static_cast<Double_t>(multN),static_cast<Double_t>(vertStatus)};
  if(fUseSparse) fTPCEventHisto->Fill(vTPCEvent);
  else FillEventHistogram(vTPCEvent);
  
    rtime +=Watch.RealTime()*1000;
    revent++;
    cout<<"For event "<<revent<<"TPC Event processing time "<<rtime/revent<<endl;
    Watch.Stop();

}


//_____________________________________________________________________________
void AliPerformanceTPC::Analyse()
{
    //
    // Analyse comparison information and store output histograms
    // in the folder "folderTPC"
    //
    TH1::AddDirectory(kFALSE);
    TH1::SetDefaultSumw2(kFALSE);
    TObjArray *aFolderObj = new TObjArray;

    if(fUseSparse){
        TString selString;
        selString = "all";
        for(Int_t i=0; i <= 2; i++) {
          AddProjection(aFolderObj, "clust", fTPCClustHisto, i, &selString);
        }
        
        //
        // event histograms
        //
        for(Int_t i=0; i<=6; i++) {
          AddProjection(aFolderObj, "event", fTPCEventHisto, i);
        }    
        AddProjection(aFolderObj, "event", fTPCEventHisto, 4, 5);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 0, 1);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 0, 3);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 1, 3);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 2, 3);

        // reconstructed vertex status > 0
        fTPCEventHisto->GetAxis(6)->SetRange(2,2);
        selString = "recVertex";
        for(Int_t i=0; i<=5; i++) {
          AddProjection(aFolderObj, "event", fTPCEventHisto, i, &selString);
        }
        AddProjection(aFolderObj, "event", fTPCEventHisto, 4, 5, &selString);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 0, 1, &selString);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 0, 3, &selString);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 1, 3, &selString);
        AddProjection(aFolderObj, "event", fTPCEventHisto, 2, 3, &selString);

        // reset cuts
        fTPCEventHisto->GetAxis(6)->SetRange(1,2);

        //
        // Track histograms 
        // 
        // all with vertex
        fTPCTrackHisto->GetAxis(8)->SetRangeUser(-1.5,1.5);
        fTPCTrackHisto->GetAxis(9)->SetRangeUser(0.5,1.5);
        selString = "all_recVertex";
        for(Int_t i=0; i <= 9; i++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, &selString);        
        }

        AddProjection(aFolderObj, "track", fTPCTrackHisto, 5, 8, &selString); 

        for(Int_t i=0; i <= 4; i++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, 5, 7, &selString);        
        }    



        // Track histograms (pos with vertex)
        fTPCTrackHisto->GetAxis(8)->SetRangeUser(0,1.5);
        selString = "pos_recVertex";
        for(Int_t i=0; i <= 9; i++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, &selString);
        }
        for(Int_t i=0; i <= 4; i++) { for(Int_t j=5; j <= 5; j++) { for(Int_t k=j+1; k <= 7; k++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, j, k, &selString);
        }  }  }
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 1, 2, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 1, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 2, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 1, 2, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 3, 4, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 5, 6, 7, &selString);
      
        // Track histograms (neg with vertex)
        fTPCTrackHisto->GetAxis(8)->SetRangeUser(-1.5,0);
        selString = "neg_recVertex";
        for(Int_t i=0; i <= 9; i++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, &selString);
        }
        for(Int_t i=0; i <= 4; i++) { for(Int_t j=5; j <= 5; j++) { for(Int_t k=j+1; k <= 7; k++) {
          AddProjection(aFolderObj, "track", fTPCTrackHisto, i, j, k, &selString);
        }  }  }
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 1, 2, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 1, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 0, 2, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 1, 2, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 3, 4, 5, &selString);
        AddProjection(aFolderObj, "track", fTPCTrackHisto, 5, 6, 7, &selString);

        //restore cuts
        fTPCTrackHisto->GetAxis(8)->SetRangeUser(-1.5,1.5);
        fTPCTrackHisto->GetAxis(9)->SetRangeUser(-0.5,1.5);
      
        printf("exportToFolder\n");
        // export objects to analysis folder
        fAnalysisFolder = ExportToFolder(aFolderObj);
        if (fFolderObj) delete fFolderObj;
        fFolderObj = aFolderObj;
        aFolderObj=0;
    }
    else{
        printf("exportToFolder\n");
        fAnalysisFolder = ExportToFolder(fFolderObj);
    }
 
}


//_____________________________________________________________________________
TFolder* AliPerformanceTPC::ExportToFolder(TObjArray * array) 
{
  // recreate folder avery time and export objects to new one
  //
  AliPerformanceTPC * comp=this;
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
Long64_t AliPerformanceTPC::Merge(TCollection* const list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
  return 0;

  if (list->IsEmpty())
  return 1;
  
  Bool_t merge = ((fgUseMergeTHnSparse && fgMergeTHnSparse) || (!fgUseMergeTHnSparse && fMergeTHnSparseObj));

  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  TObjArray* objArrayList = 0;
  objArrayList = new TObjArray();

  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
  {
    AliPerformanceTPC* entry = dynamic_cast<AliPerformanceTPC*>(obj);
    if (entry == 0) continue; 
    if (merge) {
        if ((fTPCClustHisto) && (entry->fTPCClustHisto)) { fTPCClustHisto->Add(entry->fTPCClustHisto); }
        if ((fTPCEventHisto) && (entry->fTPCEventHisto)) { fTPCEventHisto->Add(entry->fTPCEventHisto); }
        if ((fTPCTrackHisto) && (entry->fTPCTrackHisto)) { fTPCTrackHisto->Add(entry->fTPCTrackHisto); }
    }
    // the analysisfolder is only merged if present
    if (entry->fFolderObj) { objArrayList->Add(entry->fFolderObj); }

    count++;
  }
  if (fFolderObj) { fFolderObj->Merge(objArrayList); } 
  // to signal that track histos were not merged: reset
  if (!merge) {
      if(fTPCTrackHisto) fTPCTrackHisto->Reset();
      if(fTPCClustHisto)fTPCClustHisto->Reset();
      if(fTPCEventHisto) fTPCEventHisto->Reset();
  }
  // delete
  if (objArrayList)  delete objArrayList;  objArrayList=0;
return count;
}


//_____________________________________________________________________________
TFolder* AliPerformanceTPC::CreateFolder(TString name, TString title) 
{ 
// create folder for analysed histograms
//
TFolder *folder = 0;
  folder = new TFolder(name.Data(),title.Data());

  return folder;
}

//_____________________________________________________________________________
TTree* AliPerformanceTPC::CreateSummary()
{
    // implementaion removed, switched back to use AliPerformanceSummary (now called in AliPerformanceTask)
    return 0;
}

void AliPerformanceTPC::FillEventHistogram(double *vTPCEvent){

    h_tpc_event_6->Fill(vTPCEvent[6]);
    if(vTPCEvent[6]>0.001){
        h_tpc_event_recvertex_0->Fill(vTPCEvent[0]);
        h_tpc_event_recvertex_1->Fill(vTPCEvent[1]);
        h_tpc_event_recvertex_2->Fill(vTPCEvent[2]);
        h_tpc_event_recvertex_3->Fill(vTPCEvent[3]);
        h_tpc_event_recvertex_4->Fill(vTPCEvent[4]);
        h_tpc_event_recvertex_5->Fill(vTPCEvent[5]);
    }

}

void AliPerformanceTPC::FillTrackHistogram(double *vTPCTrackHisto){

    h_tpc_track_all_recvertex_5_8->Fill(vTPCTrackHisto[5],vTPCTrackHisto[8]);
    h_tpc_track_all_recvertex_1_5_7->Fill(vTPCTrackHisto[1],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    h_tpc_track_all_recvertex_2_5_7->Fill(vTPCTrackHisto[2],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    
    double q = vTPCTrackHisto[8];
    
    h_tpc_track_all_recvertex_0_5_7->Fill(vTPCTrackHisto[0],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    if(q > 0) h_tpc_track_pos_recvertex_0_5_7->Fill(vTPCTrackHisto[0],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    else h_tpc_track_neg_recvertex_0_5_7->Fill(vTPCTrackHisto[0],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    
    h_tpc_track_all_recvertex_3_5_7->Fill(vTPCTrackHisto[3],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    if(q > 0) h_tpc_track_pos_recvertex_3_5_7->Fill(vTPCTrackHisto[3],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    else h_tpc_track_neg_recvertex_3_5_7->Fill(vTPCTrackHisto[3],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    
    h_tpc_track_all_recvertex_4_5_7->Fill(vTPCTrackHisto[4],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    if(q > 0) h_tpc_track_pos_recvertex_4_5_7->Fill(vTPCTrackHisto[4],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    else h_tpc_track_neg_recvertex_4_5_7->Fill(vTPCTrackHisto[4],vTPCTrackHisto[5],vTPCTrackHisto[7]);
    
    if(q > 0)h_tpc_track_pos_recvertex_3_5_6->Fill(vTPCTrackHisto[3],vTPCTrackHisto[5],vTPCTrackHisto[6]);
    else h_tpc_track_neg_recvertex_3_5_6->Fill(vTPCTrackHisto[3],vTPCTrackHisto[5],vTPCTrackHisto[6]);
    if(q > 0)h_tpc_track_pos_recvertex_4_5_6->Fill(vTPCTrackHisto[4],vTPCTrackHisto[5],vTPCTrackHisto[6]);
    else h_tpc_track_neg_recvertex_4_5_6->Fill(vTPCTrackHisto[4],vTPCTrackHisto[5],vTPCTrackHisto[6]);
    
    if(q > 0) h_tpc_track_pos_recvertex_2_5_6->Fill(vTPCTrackHisto[2],vTPCTrackHisto[5],vTPCTrackHisto[6]);
    else h_tpc_track_neg_recvertex_2_5_6->Fill(vTPCTrackHisto[2],vTPCTrackHisto[5],vTPCTrackHisto[6]);

}

