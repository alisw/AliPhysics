#include <Riostream.h>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TList.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TParticle.h"
#include "TParticlePDG.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

#include "AliAnalysisTaskResonanceVsMultiplicityROOT6.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"



using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskResonanceVsMultiplicityROOT6)

AliAnalysisTaskResonanceVsMultiplicityROOT6::AliAnalysisTaskResonanceVsMultiplicityROOT6(): AliAnalysisTaskSE(), fTreeEvent(0), fListHist(0), fQAList(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0), fNch_eta0pt5(0), fQ1_ptmax3(0), fQ2_ptmax3(0), fQ3_ptmax3(0), fQ4_ptmax3(0), fNch_ptmax3(0), fQ1_ptmax2(0), fQ2_ptmax2(0), fQ3_ptmax2(0), fQ4_ptmax2(0), fNch_ptmax2(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist_centrality_beforecut(0), Hist_fTreeTrackVariableDcaXY(0),Hist_fTreeTrackVariableDcaZ(0),Hist_fTreeTrackVariableTpcNCls(0), Hist_fTreeTrackVariableTpcNCrossedRows(0), Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable(0),Hist_fTreeTrackVariableChiSqrPerTpcCls(0), Hist_fTreeTrackVariableChiSqrPerItsCls(0), Hist_Eta(0), Hist_Pt(0),fHistTPConlyVsV0MBeforeAliEventCut(0),fHistCL0VsV0MBeforeAliEventCut(0),fHistTPConlyVsV0MAfterAliEventCut(0),fHistCL0VsV0MAfterAliEventCut(0),fHistTPCtracksVsITStrkltsBeforeAliEventCut(0),fHistTPCtracksVsITStrkltsAfterAliEventCut(0),fHistTPConlyVsV0MBefore(0),fHistCL0VsV0MBefore(0),fHistTPConlyVsV0MAfter(0),fHistCL0VsV0MAfter(0),fHistTPCtracksVsITStrkltsBefore(0),fHistTPCtracksVsITStrkltsAfter(0),fHistTPConlyVsV0MAfterITSTPCcorelationCut(0),fHistCL0VsV0MAfterITSTPCcorelationCut(0),fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut(0),fHistTPCrefitVsITSrefitBeforeAliEventCut(0),fHistTPCrefitVsITSrefitAfterAliEventCut(0),fHistTPCrefitVsITSrefitBefore(0),fHistTPCrefitVsITSrefitAfter(0),fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut(0), fCenCutLowPU(0),fCenCutHighPU(0),fSPDCutPU(0),fV0CutPU(0),fMultCutPU(0), fTreeName(0), fEtaMax(0), fVertexZMax(0), fDCAxyMax(0), fDCAzMax(0), fChi2TPC(0), fChi2ITS(0), fNCrossedRowsTPC(0)
{
  for(int i=0;i<4;i++)
    {
      fQ1_gen[i]=-999;
      fQ2_gen[i]=-999;
      fQ3_gen[i]=-999;
      fQ4_gen[i]=-999;
      fNch_gen[i]=-999;
      fQ1_rec[i]=-999;
      fQ2_rec[i]=-999;
      fQ3_rec[i]=-999;
      fQ4_rec[i]=-999;
      fNch_rec[i]=-999;
    }
}

AliAnalysisTaskResonanceVsMultiplicityROOT6::AliAnalysisTaskResonanceVsMultiplicityROOT6(const char *name): AliAnalysisTaskSE(name), fTreeEvent(0), fListHist(0), fQAList(0), fPIDResponse(0), fESDtrackCuts(0), fEventCuts(0), fTreeVariableCentrality(0), fvertex(0), fNch_eta0pt5(0), fQ1_ptmax3(0), fQ2_ptmax3(0), fQ3_ptmax3(0), fQ4_ptmax3(0), fNch_ptmax3(0), fQ1_ptmax2(0), fQ2_ptmax2(0), fQ3_ptmax2(0), fQ4_ptmax2(0), fNch_ptmax2(0), hist2D_pt_gen_centrality(0), hist2D_pt_rec_centrality(0), Profile_mean_term1(0),Profile_var_term1(0),Profile_var_term2(0),Profile_skewness_term1(0),Profile_skewness_term2(0),Profile_skewness_term3(0),Profile_kurtosis_term1(0),Profile_kurtosis_term2(0),Profile_kurtosis_term3(0),Profile_kurtosis_term4(0), hist_centrality_beforecut(0), Hist_fTreeTrackVariableDcaXY(0),Hist_fTreeTrackVariableDcaZ(0),Hist_fTreeTrackVariableTpcNCls(0), Hist_fTreeTrackVariableTpcNCrossedRows(0), Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable(0), Hist_fTreeTrackVariableChiSqrPerTpcCls(0),Hist_fTreeTrackVariableChiSqrPerItsCls(0), Hist_Eta(0), Hist_Pt(0),fHistTPConlyVsV0MBeforeAliEventCut(0),fHistCL0VsV0MBeforeAliEventCut(0),fHistTPConlyVsV0MAfterAliEventCut(0),fHistCL0VsV0MAfterAliEventCut(0),fHistTPCtracksVsITStrkltsBeforeAliEventCut(0),fHistTPCtracksVsITStrkltsAfterAliEventCut(0),fHistTPConlyVsV0MBefore(0),fHistCL0VsV0MBefore(0),fHistTPConlyVsV0MAfter(0),fHistCL0VsV0MAfter(0),fHistTPCtracksVsITStrkltsBefore(0),fHistTPCtracksVsITStrkltsAfter(0),fHistTPConlyVsV0MAfterITSTPCcorelationCut(0),fHistCL0VsV0MAfterITSTPCcorelationCut(0),fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut(0),fHistTPCrefitVsITSrefitBeforeAliEventCut(0),fHistTPCrefitVsITSrefitAfterAliEventCut(0),fHistTPCrefitVsITSrefitBefore(0),fHistTPCrefitVsITSrefitAfter(0),fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut(0), fCenCutLowPU(0),fCenCutHighPU(0),fSPDCutPU(0),fV0CutPU(0),fMultCutPU(0), fTreeName(0), fEtaMax(0), fVertexZMax(0), fDCAxyMax(0), fDCAzMax(0), fChi2TPC(0), fChi2ITS(0), fNCrossedRowsTPC(0)
{
  for(int i=0;i<4;i++)
    {
      fQ1_gen[i]=-999;
      fQ2_gen[i]=-999;
      fQ3_gen[i]=-999;
      fQ4_gen[i]=-999;
      fNch_gen[i]=-999;
      fQ1_rec[i]=-999;
      fQ2_rec[i]=-999;
      fQ3_rec[i]=-999;
      fQ4_rec[i]=-999;
      fNch_rec[i]=-999;
    }

  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
    
}

AliAnalysisTaskResonanceVsMultiplicityROOT6::~AliAnalysisTaskResonanceVsMultiplicityROOT6()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
  if (fTreeEvent){
    delete fTreeEvent;
    fTreeEvent = 0x0;
  }
  if (fListHist){
        delete fListHist;
        fListHist = 0x0;
    }
  if(fESDtrackCuts) {
    delete fESDtrackCuts;
    fESDtrackCuts=0x0;
  }
}

//________________________________________________________________________
void AliAnalysisTaskResonanceVsMultiplicityROOT6::UserCreateOutputObjects()
{

  //------------------------------------------------
  // Particle Identification Setup
  //------------------------------------------------
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();
  
  //------------------------------------------------
  // track cut
  //------------------------------------------------
  /*
  if(!fESDtrackCuts )
    {
      fESDtrackCuts = new AliESDtrackCuts();
      fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      //fESDtrackCuts->SetPtRange(0.15,8.0);
      fESDtrackCuts->SetEtaRange(-0.8,0.8);

      //++++++++++++++Track parameters to be varied for systematics++++++++++++
      fESDtrackCuts->SetMaxDCAToVertexXY(0.1);
      fESDtrackCuts->SetMaxDCAToVertexZ(1);
      fESDtrackCuts->SetMaxChi2PerClusterTPC(2.5);
      fESDtrackCuts->SetMaxChi2PerClusterITS(36);
      fESDtrackCuts->SetMinNCrossedRowsTPC(70);
      fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36); //new
      //fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
    }
  */
  if(!fESDtrackCuts )
    {
      fESDtrackCuts = new AliESDtrackCuts();
      fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE,1);
      fESDtrackCuts->SetPtRange(0.15,5.0);
      fESDtrackCuts->SetEtaRange(-1.0*fEtaMax,fEtaMax);

      //++++++++++++++Track parameters to be varied for systematics++++++++++++++
      fESDtrackCuts->SetMaxDCAToVertexXY(fDCAxyMax);
      fESDtrackCuts->SetMaxDCAToVertexZ(fDCAzMax);
      fESDtrackCuts->SetMaxChi2PerClusterTPC(fChi2TPC);
      fESDtrackCuts->SetMaxChi2PerClusterITS(fChi2ITS);
      fESDtrackCuts->SetMinNCrossedRowsTPC(fNCrossedRowsTPC);
      fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(36); //new
      //fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
  

  
  
  OpenFile(1);
  fTreeEvent = new TTree(fTreeName,"Event");
  fTreeEvent->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
  fTreeEvent->Branch("fvertex",&fvertex,"fvertex/F");
  //fTreeEvent->Branch("fNch_eta0pt5",&fNch_eta0pt5,"fNch_eta0pt5/F");
  fTreeEvent->Branch("fNch_ptmax3", &fNch_ptmax3, "fNch_ptmax3/F");
  fTreeEvent->Branch("fQ1_ptmax3", &fQ1_ptmax3, "fQ1_ptmax3/F");
  fTreeEvent->Branch("fQ2_ptmax3", &fQ2_ptmax3, "fQ2_ptmax3/F");
  fTreeEvent->Branch("fQ3_ptmax3", &fQ3_ptmax3, "fQ3_ptmax3/F");
  fTreeEvent->Branch("fQ4_ptmax3", &fQ4_ptmax3, "fQ4_ptmax3/F");
  fTreeEvent->Branch("fNch_ptmax2", &fNch_ptmax2, "fNch_ptmax2/F");
  fTreeEvent->Branch("fQ1_ptmax2", &fQ1_ptmax2, "fQ1_ptmax2/F");
  fTreeEvent->Branch("fQ2_ptmax2", &fQ2_ptmax2, "fQ2_ptmax2/F");
  fTreeEvent->Branch("fQ3_ptmax2", &fQ3_ptmax2, "fQ3_ptmax2/F");
  fTreeEvent->Branch("fQ4_ptmax2", &fQ4_ptmax2, "fQ4_ptmax2/F");
  PostData(1, fTreeEvent);
  
  
  OpenFile(2);
  fListHist = new TList();
  fListHist->SetOwner(kTRUE);

  // hist2D_pt_gen_centrality=new TH2F("hist2D_pt_gen_centrality","hist2D_pt_gen_centrality",800,0,8,100,0,100);
  // hist2D_pt_rec_centrality=new TH2F("hist2D_pt_rec_centrality","hist2D_pt_rec_centrality",800,0,8,100,0,100);
  // hist_centrality_beforecut=new TH1D("hist_centrality_beforecut","hist_centrality_beforecut",100,0,100); 

  //QA Histograms
  Hist_fTreeTrackVariableDcaXY = new TH1D("Hist_fTreeTrackVariableDcaXY","Hist_fTreeTrackVariableDcaXY",400,-0.2,+0.2);
  Hist_fTreeTrackVariableDcaZ = new TH1D("Hist_fTreeTrackVariableDcaZ","Hist_fTreeTrackVariableDcaZ",400,-2,+2);
  Hist_fTreeTrackVariableTpcNCls = new TH1D("Hist_fTreeTrackVariableTpcNCls","Hist_fTreeTrackVariableTpcNCls",200,0,200);
  Hist_fTreeTrackVariableTpcNCrossedRows = new TH1D("Hist_fTreeTrackVariableTpcNCrossedRows","Hist_fTreeTrackVariableTpcNCrossedRows",200,0,200);
  Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable = new TH1D("Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable","Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable",200,0,200);
  Hist_fTreeTrackVariableChiSqrPerTpcCls = new TH1D("Hist_fTreeTrackVariableChiSqrPerTpcCls","Hist_fTreeTrackVariableChiSqrPerTpcCls",50,0,5);
  Hist_fTreeTrackVariableChiSqrPerItsCls = new TH1D("Hist_fTreeTrackVariableChiSqrPerItsCls","Hist_fTreeTrackVariableChiSqrPerItsCls",500,0,50);
  Hist_Eta = new TH1D("Hist_Eta","Hist_Eta", 40, -1, +1);
  Hist_Pt = new TH1D("Hist_Pt","Hist_Pt", 500, 0, 5);

  //PileupCut histograms
  fHistTPConlyVsV0MBeforeAliEventCut=new TH2D("fHistTPConlyVsV0MBeforeAliEventCut","fHistTPConlyVsV0MBeforeAliEventCut",100,0,100,250,0,5000);   
  fHistTPConlyVsV0MAfterAliEventCut=new TH2D("fHistTPConlyVsV0MAfterAliEventCut","fHistTPConlyVsV0MAfterAliEventCut",100,0,100,250,0,5000);   
  fHistCL0VsV0MBeforeAliEventCut=new TH2D("fHistCL0VsV0MBeforeAliEventCut","fHistCL0VsV0MBeforeAliEventCut",100,0,100,100,0,100);
  fHistCL0VsV0MAfterAliEventCut=new TH2D("fHistCL0VsV0MAfterAliEventCut","fHistCL0VsV0MAfterAliEventCut",100,0,100,100,0,100);
  fHistTPCtracksVsITStrkltsBeforeAliEventCut=new TH2D("fHistTPCtracksVsITStrkltsBeforeAliEventCut","fHistTPCtracksVsITStrkltsBeforeAliEventCut",250,0,5000,250,0,5000);
  fHistTPCtracksVsITStrkltsAfterAliEventCut=new TH2D("fHistTPCtracksVsITStrkltsAfterAliEventCut","fHistTPCtracksVsITStrkltsAfterAliEventCut",250,0,5000,250,0,5000);
  fHistTPConlyVsV0MBefore=new TH2D("fHistTPConlyVsV0MBefore","fHistTPConlyVsV0MBefore",100,0,100,250,0,5000);   
  fHistTPConlyVsV0MAfter=new TH2D("fHistTPConlyVsV0MAfter","fHistTPConlyVsV0MAfter",100,0,100,250,0,5000);   
  fHistCL0VsV0MBefore=new TH2D("fHistCL0VsV0MBefore","fHistCL0VsV0MBefore",100,0,100,100,0,100);
  fHistCL0VsV0MAfter=new TH2D("fHistCL0VsV0MAfter","fHistCL0VsV0MAfter",100,0,100,100,0,100);
  fHistTPCtracksVsITStrkltsBefore=new TH2D("fHistTPCtracksVsITStrkltsBefore","fHistTPCtracksVsITStrkltsBefore",250,0,5000,250,0,5000);
  fHistTPCtracksVsITStrkltsAfter=new TH2D("fHistTPCtracksVsITStrkltsAfter","fHistTPCtracksVsITStrkltsAfter",250,0,5000,250,0,5000);
  fHistTPConlyVsV0MAfterITSTPCcorelationCut=new TH2D("fHistTPConlyVsV0MAfterITSTPCcorelationCut","fHistTPConlyVsV0MAfterITSTPCcorelationCut",100,0,100,250,0,5000);
  fHistCL0VsV0MAfterITSTPCcorelationCut=new TH2D("fHistCL0VsV0MAfterITSTPCcorelationCut","fHistCL0VsV0MAfterITSTPCcorelationCut",100,0,100,100,0,100);
  fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut=new TH2D("fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut","fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut",250,0,5000,250,0,5000);
  fHistTPCrefitVsITSrefitBeforeAliEventCut=new TH2D("fHistTPCrefitVsITSrefitBeforeAliEventCut","fHistTPCrefitVsITSrefitBeforeAliEventCut",250,0,5000,250,0,5000);
  fHistTPCrefitVsITSrefitAfterAliEventCut=new TH2D("fHistTPCrefitVsITSrefitAfterAliEventCut","fHistTPCrefitVsITSrefitAfterAliEventCut",250,0,5000,250,0,5000);
  fHistTPCrefitVsITSrefitBefore=new TH2D("fHistTPCrefitVsITSrefitBefore","fHistTPCrefitVsITSrefitBefore",250,0,5000,250,0,5000);
  fHistTPCrefitVsITSrefitAfter=new TH2D("fHistTPCrefitVsITSrefitAfter","fHistTPCrefitVsITSrefitAfter",250,0,5000,250,0,5000);
  fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut=new TH2D("fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut","fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut",250,0,5000,250,0,5000);



  // fListHist->Add(hist2D_pt_gen_centrality);
  // fListHist->Add(hist2D_pt_rec_centrality);
  // fListHist->Add(hist_centrality_beforecut);
  fListHist->Add(Hist_fTreeTrackVariableDcaXY);
  fListHist->Add(Hist_fTreeTrackVariableDcaZ);
  fListHist->Add(Hist_fTreeTrackVariableTpcNCls);
  fListHist->Add(Hist_fTreeTrackVariableTpcNCrossedRows);
  fListHist->Add(Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable);
  fListHist->Add(Hist_fTreeTrackVariableChiSqrPerTpcCls);
  fListHist->Add(Hist_fTreeTrackVariableChiSqrPerItsCls);
  fListHist->Add(Hist_Eta);
  fListHist->Add(Hist_Pt);

  fListHist->Add(fHistTPConlyVsV0MBeforeAliEventCut);
  fListHist->Add(fHistTPConlyVsV0MAfterAliEventCut);
  fListHist->Add(fHistCL0VsV0MBeforeAliEventCut);
  fListHist->Add(fHistCL0VsV0MAfterAliEventCut);
  fListHist->Add(fHistTPCtracksVsITStrkltsBeforeAliEventCut);
  fListHist->Add(fHistTPCtracksVsITStrkltsAfterAliEventCut);
  fListHist->Add(fHistTPCrefitVsITSrefitBeforeAliEventCut);
  fListHist->Add(fHistTPCrefitVsITSrefitAfterAliEventCut);
  fListHist->Add(fHistTPConlyVsV0MBefore);
  fListHist->Add(fHistTPConlyVsV0MAfter);
  fListHist->Add(fHistCL0VsV0MBefore);
  fListHist->Add(fHistCL0VsV0MAfter);
  fListHist->Add(fHistTPCtracksVsITStrkltsBefore);
  fListHist->Add(fHistTPCtracksVsITStrkltsAfter);
  fListHist->Add(fHistTPCrefitVsITSrefitBefore);
  fListHist->Add(fHistTPCrefitVsITSrefitAfter);
  fListHist->Add(fHistTPConlyVsV0MAfterITSTPCcorelationCut);
  fListHist->Add(fHistCL0VsV0MAfterITSTPCcorelationCut);
  fListHist->Add(fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut);
  fListHist->Add(fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut);
  PostData(2, fListHist);

  OpenFile(3);
  fQAList = new TList();
  fEventCuts.AddQAplotsToList(fQAList,kTRUE);
  PostData(3, fQAList); 
  fQAList->SetOwner(kTRUE);
  


  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //// PileUp Removal Functions:
  fSPDCutPU = new TF1("fSPDCutPU", "400. + 4.*x", 0, 10000);

  Double_t parV0[8] = {43.8011, 0.822574, 8.49794e-02, 1.34217e+02, 7.09023e+00, 4.99720e-02, -4.99051e-04, 1.55864e-06};
  fV0CutPU = new TF1("fV0CutPU", "[0]+[1]*x - 6.*[2]*([3] + [4]*sqrt(x) + [5]*x + [6]*x*sqrt(x) + [7]*x*x)", 0, 100000);
  fV0CutPU->SetParameters(parV0);

  Double_t parFB32[8] = {2093.36, -66.425, 0.728932, -0.0027611, 1.01801e+02, -5.23083e+00, -1.03792e+00, 5.70399e-03};
  fMultCutPU = new TF1("fMultCutPU", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 6.*([4]+[5]*sqrt(x)+[6]*x+[7]*x*x)", 0, 90);
  fMultCutPU->SetParameters(parFB32);

  Double_t parV0CL0[6] = {0.320462, 0.961793, 1.02278, 0.0330054, -0.000719631, 6.90312e-06};
  fCenCutLowPU = new TF1("fCenCutLowPU", "[0]+[1]*x - 6.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutLowPU->SetParameters(parV0CL0);
  fCenCutHighPU = new TF1("fCenCutHighPU", "[0]+[1]*x + 5.5*([2]+[3]*x+[4]*x*x+[5]*x*x*x)", 0, 100);
  fCenCutHighPU->SetParameters(parV0CL0);
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
}


//________________________________________________________________________
void AliAnalysisTaskResonanceVsMultiplicityROOT6::UserExec(Option_t *)
{



  AliESDEvent *lESDevent = 0x0;

  lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
  if (!lESDevent) {
    AliWarning("ERROR: lESDevent not available \n");
    PostData(1, fTreeEvent);
    PostData(2, fListHist);
    PostData(3, fQAList);
    return;
  }


  cout<<"******************Before trigger calling******************************"<<endl;

  ////tigger/////////////                                                                                                                     
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;

  if ( ! isSelected)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }


  cout<<"******************Before vertex calling******************************"<<endl;
  // primary vertex                                                                                                                          
  //                                                                                                                                         
  const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1)
    {
      // SPD vertex                                                                                                                           
      vertex = lESDevent->GetPrimaryVertexSPD();
      if (vertex->GetNContributors() < 1)
        vertex = 0x0;
    }
  if (!vertex)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }
  if (TMath::Abs(vertex->GetZ()) > fVertexZMax)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }


  //**********************************************************************************************************************************        
  Float_t v0Centr = -100.;
  Float_t cl1Centr = -100.;
  Float_t cl0Centr = -100.;
  AliMultSelection* multSelection = 0x0;
  multSelection = (AliMultSelection*)lESDevent->FindListObject("MultSelection");
  if(!multSelection)
    {
      AliWarning("AliMultSelection object not found!\n\n");
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }
  else
    {
      v0Centr = multSelection->GetMultiplicityPercentile("V0M");
      cl1Centr = multSelection->GetMultiplicityPercentile("CL1");
      cl0Centr = multSelection->GetMultiplicityPercentile("CL0");
    }

  if(v0Centr>=90.||v0Centr<0)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return; //This would have to be adjusted for vs. V0M                                                                                    
    }
  
  AliESDtrackCuts *fESDTPConlyTracks= new AliESDtrackCuts();
  fESDTPConlyTracks = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  Int_t nTPCtracks = lESDevent->GetNumberOfTPCTracks();
  //cout<<"TPC tracks:"<<nTPCtracks<<endl;
  Int_t nITSClsLy0 = lESDevent->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = lESDevent->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
  Int_t nITSTrkls = lESDevent->GetMultiplicity()->GetNumberOfTracklets();
  //cout<<"No. of ITS tracklets:"<<nITSTrkls<<endl;
  const Int_t nTracks = lESDevent->GetNumberOfTracks();
  Int_t multTrk = 0;
  Int_t multTPConlyTrk = 0;
  Int_t nTracksTPC = 0;
  Int_t nTracksITS = 0;
  AliESDtrack *esdTrack;
  for (Int_t it = 0; it < nTracks; it++) {
    esdTrack = (AliESDtrack*)lESDevent->GetTrack(it);
    if(!esdTrack) continue;
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) continue;
    if(fESDtrackCuts->AcceptTrack(esdTrack)) multTrk++;
    if(fESDTPConlyTracks->AcceptTrack(esdTrack)) multTPConlyTrk++;

    if(!esdTrack->GetInnerParam()) continue;
    if(esdTrack->IsOn(AliESDtrack::kTPCrefit)) nTracksTPC++;
    if(esdTrack->IsOn(AliESDtrack::kITSrefit)) nTracksITS++;
  }
  //cout<<"TPC tracks using cut: "<<multTPConlyTrk<<endl;

  //Histograms before AliEventCut imposed
  fHistTPConlyVsV0MBeforeAliEventCut->Fill(v0Centr,multTrk);
  fHistCL0VsV0MBeforeAliEventCut->Fill(v0Centr,cl0Centr);
  fHistTPCtracksVsITStrkltsBeforeAliEventCut->Fill(nITSTrkls,multTPConlyTrk);
  fHistTPCrefitVsITSrefitBeforeAliEventCut->Fill(nTracksTPC,nTracksITS);
 
  //***************************************************************************************************************************************************************

  cout<<"******************Before pileup calling******************************"<<endl;
  //Pileup
  Bool_t EventAccepted;
  fEventCuts.SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE); 
  EventAccepted = fEventCuts.AcceptEvent(lESDevent);                                                                                                                                                     
  if (!EventAccepted)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }

  //Histograms after AliEventCut imposed
  fHistTPConlyVsV0MAfterAliEventCut->Fill(v0Centr,multTrk);
  fHistCL0VsV0MAfterAliEventCut->Fill(v0Centr,cl0Centr);
  fHistTPCtracksVsITStrkltsAfterAliEventCut->Fill(nITSTrkls,multTPConlyTrk);
  fHistTPCrefitVsITSrefitAfterAliEventCut->Fill(nTracksTPC,nTracksITS);

  
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //Strict pileup cut
  Bool_t AcceptEventStatus;
  AcceptEventStatus=AcceptEventAfterPileUpCut(lESDevent);

  if(!AcceptEventStatus)
    {
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  //Histograms after ITS TPC correlation cut
  fHistTPConlyVsV0MAfterITSTPCcorelationCut->Fill(v0Centr,multTrk);
  fHistCL0VsV0MAfterITSTPCcorelationCut->Fill(v0Centr,cl0Centr);
  fHistTPCtracksVsITStrkltsAfterITSTPCcorelationCut->Fill(nITSTrkls,multTPConlyTrk);
  fHistTPCrefitVsITSrefitAfterITSTPCcorelationCut->Fill(nTracksTPC,nTracksITS);


  
  //////centrality selection/////////                                                                                                         
  Float_t lV0M;
  Int_t lEvSelCode = 100;
  AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
  if( !MultSelection)
    {
      AliWarning("AliMultSelection object not found!");
      PostData(1, fTreeEvent);
      PostData(2, fListHist);
      PostData(3, fQAList);
      return;
    }
  else
    {
      lV0M = MultSelection->GetMultiplicityPercentile("V0M");
    }



 
  cout<<"fTreeVariableCentrality:****************** : "<<lV0M<<endl;

  
  //----- Loop on track ----------------------------------------------------------------
  Int_t ntracks = lESDevent->GetNumberOfTracks();

  Float_t mean_term1[2]={-100.0,-100.0};
  Float_t var_term1[2]={-100.0,-100.0};
  Float_t var_term2[2]={-100.0,-100.0};
  Float_t skewness_term1[2]={-100.0,-100.0};
  Float_t skewness_term2[2]={-100.0,-100.0};
  Float_t skewness_term3[2]={-100.0,-100.0};
  Float_t kurtosis_term1[2]={-100.0,-100.0};
  Float_t kurtosis_term2[2]={-100.0,-100.0};
  Float_t kurtosis_term3[2]={-100.0,-100.0};
  Float_t kurtosis_term4[2]={-100.0,-100.0};


  Float_t Q1[2]={0.0,0.0};
  Float_t Q2[2]={0.0,0.0};
  Float_t Q3[2]={0.0,0.0};
  Float_t Q4[2]={0.0,0.0};
  Float_t Nch[2]={0.0,0.0};

  Double_t p[3];
  Int_t trackcharge=0;
  Float_t Track_pt=0.0;
  Float_t Track_eta=0.0;

  Float_t dcaxy     = 0.0;
  Float_t dcaz      = 0.0;
  Int_t nClustersITS = -1;
  Int_t nClustersTPC = -1;
  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  Float_t nCrossedRowsTPC = -1;
  Float_t  ratioCrossedRowsOverFindableClustersTPC = -1;

  Int_t str = 0;


  cout<<"############# No of tracks in event: ########"<<ntracks<<endl;
  for(Int_t itr = 0; itr < ntracks; itr++)
    {
      AliVTrack   *track = (AliVTrack*)lESDevent->GetTrack(itr);
      if(!track)      continue;
      AliESDtrack *esdt  = dynamic_cast<AliESDtrack*>(track);
      if(!esdt)      continue;
      if(!fESDtrackCuts->AcceptTrack(esdt))  continue;
      track->PxPyPz(p); 
      trackcharge=track->Charge();
      Track_eta=track->Eta();
      Track_pt=TMath::Sqrt((p[0]*p[0])+(p[1]*p[1]));
      
      //******************extra track parameters*************************//
      track->GetImpactParameters(dcaxy,dcaz);
      nClustersITS = track->GetITSclusters(0);
      nClustersTPC = -1;
      if(fESDtrackCuts->GetRequireTPCStandAlone()) {
  	  nClustersTPC = esdt->GetTPCNclsIter1();
  	}
  	else {
  	  nClustersTPC = esdt->GetTPCclusters(0);
  	}
  	chi2PerClusterITS = -1;
  	chi2PerClusterTPC = -1;
	
  	if (nClustersITS!=0)
  	  chi2PerClusterITS = esdt->GetITSchi2()/Float_t(nClustersITS);
	
  	if (nClustersTPC!=0) {
  	  if(fESDtrackCuts->GetRequireTPCStandAlone()) {
  	    chi2PerClusterTPC = esdt->GetTPCchi2Iter1()/Float_t(nClustersTPC);
  	  } else {
  	    chi2PerClusterTPC = esdt->GetTPCchi2()/Float_t(nClustersTPC);
  	  }
  	}
  	nCrossedRowsTPC = track->GetTPCCrossedRows();
  	ratioCrossedRowsOverFindableClustersTPC = 1.0;
  	if (track->GetTPCNclsF()>0) {
  	  ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / track->GetTPCNclsF();
  	}

/*
  // 	fTreeTrackVariableDcaXY[str] = dcaxy; 
  // 	fTreeTrackVariableDcaZ[str] = dcaz; 
  // 	fTreeTrackVariableTpcNCls[str] = nClustersTPC; 
  // 	fTreeTrackVariableTpcNCrossedRows[str] = nCrossedRowsTPC; 
  // 	fTreeTrackVariableLeastRatioCrossedRowsOverFindable[str] = ratioCrossedRowsOverFindableClustersTPC; 
  // 	fTreeTrackVariableChiSqrPerTpcCls[str] = chi2PerClusterTPC; 
  // 	fTreeTrackVariableChiSqrPerItsCls[str] = chi2PerClusterITS; 
*/
  	Hist_fTreeTrackVariableDcaXY->Fill(dcaxy);
  	Hist_fTreeTrackVariableDcaZ->Fill(dcaz);
  	Hist_fTreeTrackVariableTpcNCls->Fill(nClustersTPC);
  	Hist_fTreeTrackVariableTpcNCrossedRows->Fill(nCrossedRowsTPC);
  	Hist_fTreeTrackVariableLeastRatioCrossedRowsOverFindable->Fill( ratioCrossedRowsOverFindableClustersTPC);
  	Hist_fTreeTrackVariableChiSqrPerTpcCls->Fill(chi2PerClusterTPC);
  	Hist_fTreeTrackVariableChiSqrPerItsCls->Fill(chi2PerClusterITS);
  	Hist_Eta->Fill(Track_eta);
  	Hist_Pt->Fill(Track_pt);
  	
  	//************************************************************************



      if(TMath::Abs(Track_eta)>0.8)continue;
      if(TMath::Abs(trackcharge)!=1)continue;
      
      if(Track_pt>0.15 && Track_pt<=2.0)
  	{
  	  Q1[0]=Q1[0]+TMath::Power(Track_pt, 1.0);
  	  Q2[0]=Q2[0]+TMath::Power(Track_pt, 2.0);
  	  Q3[0]=Q3[0]+TMath::Power(Track_pt, 3.0);
  	  Q4[0]=Q4[0]+TMath::Power(Track_pt, 4.0);
  	  Nch[0]=Nch[0]+1;
  	}
      if(Track_pt>0.2 && Track_pt<=3.0)
  	{
  	  Q1[1]=Q1[1]+TMath::Power(Track_pt, 1.0);
  	  Q2[1]=Q2[1]+TMath::Power(Track_pt, 2.0);
  	  Q3[1]=Q3[1]+TMath::Power(Track_pt, 3.0);
  	  Q4[1]=Q4[1]+TMath::Power(Track_pt, 4.0);
  	  Nch[1]=Nch[1]+1;
  	}
    }
  //end of track loop

  // //Filling tree variables                                                                                                                    
  fTreeVariableCentrality=lV0M;
  fvertex=vertex->GetZ();
  fQ1_ptmax2=Q1[0];
  fQ2_ptmax2=Q2[0];
  fQ3_ptmax2=Q3[0];
  fQ4_ptmax2=Q4[0];
  fNch_ptmax2=Nch[0];
  fQ1_ptmax3=Q1[1];
  fQ2_ptmax3=Q2[1];
  fQ3_ptmax3=Q3[1];
  fQ4_ptmax3=Q4[1];
  fNch_ptmax3=Nch[1];
  
  fTreeEvent->Fill();
  // Post output data.
  PostData(1, fTreeEvent); 
  PostData(2, fListHist);
  PostData(3, fQAList);
}

//____________________________________________________________________________________________________________________________________
//PILEUP selection

Bool_t  AliAnalysisTaskResonanceVsMultiplicityROOT6::AcceptEventAfterPileUpCut(AliESDEvent* fESD) { //From Alex

  Bool_t BisPileUp=kTRUE;

  Float_t v0Centr = -100.;
  Float_t cl1Centr = -100.;
  Float_t cl0Centr = -100.;
  AliMultSelection* MultSelection = 0x0;
  MultSelection = (AliMultSelection*)fESD->FindListObject("MultSelection");
  if(!MultSelection)
    {
      AliWarning("AliMultSelection object not found!\n\n");
      return kFALSE;
    }
  else
    {
      v0Centr = MultSelection->GetMultiplicityPercentile("V0M");
      cl1Centr = MultSelection->GetMultiplicityPercentile("CL1");
      cl0Centr = MultSelection->GetMultiplicityPercentile("CL0");
    }
  if(v0Centr>=90.||v0Centr<0) return kFALSE; //This would have to be adjusted for vs. V0M

  AliESDtrackCuts *fESDTPConlyTracks= new AliESDtrackCuts();
  fESDTPConlyTracks = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
  Int_t nTPCtracks = fESD->GetNumberOfTPCTracks();
  //cout<<"TPC tracks:"<<nTPCtracks<<endl;
  Int_t nITSClsLy0 = fESD->GetNumberOfITSClusters(0);
  Int_t nITSClsLy1 = fESD->GetNumberOfITSClusters(1);
  Int_t nITSCls = nITSClsLy0 + nITSClsLy1;
  Int_t nITSTrkls = fESD->GetMultiplicity()->GetNumberOfTracklets();
  //cout<<"No. of ITS tracklets:"<<nITSTrkls<<endl;
  const Int_t nTracks = fESD->GetNumberOfTracks();
  Int_t multTrk = 0;
  Int_t multTPConlyTrk = 0;
  Int_t nTracksTPC = 0;
  Int_t nTracksITS = 0;
  AliESDtrack *esdTrack;
  for (Int_t it = 0; it < nTracks; it++) {
    esdTrack = (AliESDtrack*)fESD->GetTrack(it);
    if(!esdTrack) continue;
    if(!fESDtrackCuts->AcceptTrack(esdTrack)) continue;
    if(fESDtrackCuts->AcceptTrack(esdTrack)) multTrk++;
    if(fESDTPConlyTracks->AcceptTrack(esdTrack)) multTPConlyTrk++;

    if(!esdTrack->GetInnerParam()) continue;
    if(esdTrack->IsOn(AliESDtrack::kTPCrefit)) nTracksTPC++;
    if(esdTrack->IsOn(AliESDtrack::kITSrefit)) nTracksITS++;
  }
  //cout<<"TPC tracks using cut: "<<multTPConlyTrk<<endl;
  AliESDVZERO* esdV0 = fESD->GetVZEROData();
  Float_t multV0a = esdV0->GetMTotV0A();
  Float_t multV0c = esdV0->GetMTotV0C();
  Float_t multV0Tot = multV0a + multV0c;
  UShort_t multV0aOn = esdV0->GetTriggerChargeA();
  UShort_t multV0cOn = esdV0->GetTriggerChargeC();
  UShort_t multV0On = multV0aOn + multV0cOn;
  //pile-up cuts
  if(cl0Centr<fCenCutLowPU->Eval(v0Centr)) //return kFALSE;
    BisPileUp=kFALSE;
  if (cl0Centr > fCenCutHighPU->Eval(v0Centr)) //return kFALSE;
    BisPileUp=kFALSE;
  if(Float_t(nITSCls)>fSPDCutPU->Eval(nITSTrkls)) //return kFALSE;
    BisPileUp=kFALSE;
  if(multV0On<fV0CutPU->Eval(multV0Tot)) //return kFALSE; //Problematic for MC for whatever reason? On AODs work perfectly fine
    BisPileUp=kFALSE;
  if(Float_t(multTrk)<fMultCutPU->Eval(v0Centr)) //return kFALSE;
    BisPileUp=kFALSE;
  AliESDtrackCuts::MultEstTrackType estType = fESD->GetPrimaryVertexTracks()->GetStatus() ? AliESDtrackCuts::kTrackletsITSTPC : AliESDtrackCuts::kTracklets;
  if(AliESDtrackCuts::GetReferenceMultiplicity(fESD,estType,0.8) < 0) //return kFALSE;
    BisPileUp=kFALSE;
  if(fESD->IsIncompleteDAQ()) //return kFALSE;
    BisPileUp=kFALSE;

  //Histograms before pile up cut imposed
  fHistTPConlyVsV0MBefore->Fill(v0Centr,multTrk);
  fHistCL0VsV0MBefore->Fill(v0Centr,cl0Centr);
  fHistTPCtracksVsITStrkltsBefore->Fill(nITSTrkls,multTPConlyTrk);
  fHistTPCrefitVsITSrefitBefore->Fill(nTracksTPC,nTracksITS);

  //Histograms after pile up cut imposed
  if(BisPileUp)
    {
      fHistTPConlyVsV0MAfter->Fill(v0Centr,multTrk);
      fHistCL0VsV0MAfter->Fill(v0Centr,cl0Centr);
      fHistTPCtracksVsITStrkltsAfter->Fill(nITSTrkls,multTPConlyTrk);
      fHistTPCrefitVsITSrefitAfter->Fill(nTracksTPC,nTracksITS);
    }



  return BisPileUp;
}

//____________________________________________________________________________________________________________________________________


//________________________________________________________________________
void AliAnalysisTaskResonanceVsMultiplicityROOT6::Terminate(Option_t *)
{
  
}
//----------------------------------------------------------------------------
Double_t AliAnalysisTaskResonanceVsMultiplicityROOT6::MyRapidity(Double_t rE, Double_t rPz) const
{
  // Local calculation for rapidity
  Double_t ReturnValue = -100;
  if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){ 
    ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
  }
  return ReturnValue;
} 
//----------------------------------------------------------------------------
void AliAnalysisTaskResonanceVsMultiplicityROOT6::SetPersonalESDtrackCuts(AliESDtrackCuts* trackcuts){
  if(fESDtrackCuts){ 
    delete fESDtrackCuts;
    fESDtrackCuts = 0x0;
  }
  fESDtrackCuts = trackcuts;
}
//----------------------------------------------------------------------------
Double_t AliAnalysisTaskResonanceVsMultiplicityROOT6::GetTOFBeta(AliVTrack *vtrack)
{
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  if(!esdtrack) return -1;
  const Double_t c = 2.99792457999999984e-02; 
  Double_t p = esdtrack->GetTPCmomentum();
  Double_t l = esdtrack->GetIntegratedLength();
  Double_t trackT0 = fPIDResponse->GetTOFResponse().GetStartTime(p);
  Double_t timeTOF = esdtrack->GetTOFsignal()- trackT0;
  Double_t mass_square=(p*p)*(TMath::Power(c*timeTOF/l,2.0)-1);
  return mass_square;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskResonanceVsMultiplicityROOT6::MatchTOF(AliVTrack *vtrack)
{
  if (!vtrack) {
    AliWarning("NULL argument: impossible to check status");
    return kFALSE;
  }
  if (!(vtrack->GetStatus() & AliESDtrack::kTOFout)) return kFALSE;
  if (!(vtrack->GetStatus() & AliESDtrack::kTIME  )) return kFALSE;

  // if (!(vtrack->GetStatus() & AliESDtrack::kTOFpid)) return kFALSE;
  float probMis = fPIDResponse->GetTOFMismatchProbability(vtrack);
  if(probMis>0.01) return kFALSE;
  
  AliESDtrack *esdtrack  = dynamic_cast<AliESDtrack*>(vtrack);
  Double_t l = esdtrack->GetIntegratedLength();
  if(l<350) return kFALSE;
  
  return kTRUE;
}
