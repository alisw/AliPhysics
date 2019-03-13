#include "AliAnalysisTaskEMCALMesonGGSDMpPb.h"

// ROOT includes
#include <vector>
#include <Riostream.h>
#include <TChain.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

// STEER? includes
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"
#include "AliTrackerBase.h"

// EMCAL includes
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliTrackerBase.h"
#include "AliEMCALCalibTimeDepCorrection.h" // Run dependent
#include "AliEMCALPIDUtils.h"
#include "AliExternalTrackParam.h"

#include "AliCentrality.h"

using std::cout;
using std::endl;


ClassImp(AliAnalysisTaskEMCALMesonGGSDMpPb)

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDMpPb::AliAnalysisTaskEMCALMesonGGSDMpPb() : 
  AliAnalysisTaskSE(),
  fOutput(0),
  fMcMode(0),
  fRecalibrator(0),
  fdRmin_ClustTrack(0),
  fPhimin(0),
  fPhimax(0),
  fEtamin(0),
  fEtamax(0),
  fTrackCuts(0),
  fEsdEv(0),
  fAodEv(0),
  h1_zvtx(0), 
  h1_trigger(0), 
  h1_centrality(0), 
  h2_PhiEtaCluster(0), 
  h2_PhiEtaClusterCut(0), 
  h2_PhiEtaMaxCell(0), 
  h2_PhiEtaMaxCellCut(0), 
  h2_gE_RecTruth(0), 
  h2_eop_E(0),
  h2_eop_pT(0),
  h2_E_time(0),
  h2_Pi0TruthPhiEta(0), 
  h2_PriPi0TruthPhiEta(0), 
  h2_Pi0TruthPhiEtaEmcal(0), 
  h2_PriPi0TruthPhiEtaEmcal(0), 
  h2_Pi0TruthPhiEta_Phi2piEta065(0), 
  h2_Pi0TruthPhiEta_Phi2piEta1(0), 
  h2_TruthPhotonsPhiEta(0),
  h2_PhotonsPhiEtaIsEmcal(0),
  TriggerList(0),
  fHelperClass(0)
{
  // Dummy constructor ALWAYS needed for I/O.
  for(int i=0; i<cent_bins; i++){
    h1_nClusters[i] = 0;
    h1_M[i] = 0;
    h1_M_mix[i] = 0;
    h1_E[i] = 0;
    h1_dR_ClustTrk[i] = 0;
    h1_Pi0TruthPt[i] = 0;
    h1_PriPi0TruthPt[i] = 0;
    h1_Pi0TruthPtEmcal[i] = 0;
    h1_PriPi0TruthPtEmcal[i] = 0;
    h1_Pi0TruthPtPhi2piEta065[i] = 0;
    h1_Pi0TruthPtPhi2piEta1[i] = 0;
    h1_TruthPhotonsEmcal[i] = 0;
    h1_PhotonsEmcal[i] = 0;
    h1_PhotonsNCellsCut[i] = 0;
    h1_PhotonsTrackMatchCut[i] = 0;
    h1_PhotonsAllCut[i] = 0;
    h1_dR_RealMC[i] = 0;
    h1_Chi2[i] = 0;
    h1_nTrkMatch[i] = 0;
    h1_nCells[i] = 0;
    h1_ClusterDisp[i] = 0;
    h2_Ellipse[i] = 0;
    h2_EtaPt[i] = 0;
    h3_MptAsymm[i] = 0;
    h3_MptAsymm_mix[i] = 0;
    h2_dphi_deta[i] = 0;
    h2_dphi_deta_mix[i] = 0;
    h2_DispRes[i] = 0;
    h2_cells_M02[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDMpPb::AliAnalysisTaskEMCALMesonGGSDMpPb(const char *name) :
  AliAnalysisTaskSE(name),
  fOutput(0),
  fMcMode(0),
  fRecalibrator(0),
  fdRmin_ClustTrack(0),
  fPhimin(0),
  fPhimax(0),
  fEtamin(0),
  fEtamax(0),
  fTrackCuts(0),
  fEsdEv(0),
  fAodEv(0),
  h1_zvtx(0), 
  h1_trigger(0), 
  h1_centrality(0), 
  h2_PhiEtaCluster(0), 
  h2_PhiEtaClusterCut(0), 
  h2_PhiEtaMaxCell(0), 
  h2_PhiEtaMaxCellCut(0), 
  h2_gE_RecTruth(0), 
  h2_eop_E(0),
  h2_eop_pT(0),
  h2_E_time(0),
  h2_Pi0TruthPhiEta(0), 
  h2_PriPi0TruthPhiEta(0), 
  h2_Pi0TruthPhiEtaEmcal(0), 
  h2_PriPi0TruthPhiEtaEmcal(0), 
  h2_Pi0TruthPhiEta_Phi2piEta065(0), 
  h2_Pi0TruthPhiEta_Phi2piEta1(0), 
  h2_TruthPhotonsPhiEta(0),
  h2_PhotonsPhiEtaIsEmcal(0),
  TriggerList(0),
  fHelperClass(0)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());                                            // for output list
  for(int i=0; i<cent_bins; i++){
    h1_nClusters[i] = 0;
    h1_M[i] = 0;
    h1_M_mix[i] = 0;
    h1_E[i] = 0;
    h1_dR_ClustTrk[i] = 0;
    h1_Pi0TruthPt[i] = 0;
    h1_PriPi0TruthPt[i] = 0;
    h1_Pi0TruthPtEmcal[i] = 0;
    h1_PriPi0TruthPtEmcal[i] = 0;
    h1_Pi0TruthPtPhi2piEta065[i] = 0;
    h1_Pi0TruthPtPhi2piEta1[i] = 0;
    h1_TruthPhotonsEmcal[i] = 0;
    h1_PhotonsEmcal[i] = 0;
    h1_PhotonsNCellsCut[i] = 0;
    h1_PhotonsTrackMatchCut[i] = 0;
    h1_PhotonsAllCut[i] = 0;
    h1_dR_RealMC[i] = 0;
    h1_Chi2[i] = 0;
    h1_nTrkMatch[i] = 0;
    h1_nCells[i] = 0;
    h1_ClusterDisp[i] = 0;
    h2_Ellipse[i] = 0;
    h2_EtaPt[i] = 0;
    h3_MptAsymm[i] = 0;
    h3_MptAsymm_mix[i] = 0;
    h2_dphi_deta[i] = 0;
    h2_dphi_deta_mix[i] = 0;
    h2_DispRes[i] = 0;
    h2_cells_M02[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDMpPb::~AliAnalysisTaskEMCALMesonGGSDMpPb()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
  delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALMesonGGSDMpPb::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
   
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);

  cout << "__________AliAnalysisTaskEMCALMesonGGSDMpPb: Input settings__________" << endl;
  cout << " fMcMode:             " << fMcMode   << endl;
  cout << " fRecalibrator:       " << fRecalibrator << endl;
  cout << " dRmin_ClustTrack:    " << fdRmin_ClustTrack << endl;
  cout << " phi range:           " << fPhimin << ", " << fPhimax << endl;
  cout << " eta range:           " << fEtamin << ", " << fEtamax << endl;
  cout << " number of zvtx bins: " << zvtx_bins << endl;
  cout << " number of mult bins: " << mult_bins << endl;
  cout << " poolDepth:           " << poolDepth << endl;
  cout << endl;
  
  char saythis1[500];
  char saythis2[500];

  double TotalNBins = 0.0;

  // Create histograms
  Int_t nClustersbins = 501;
  Float_t nClusterslow = -0.5, nClustersup = 500.5;
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_nClusters_%d",i);
    sprintf(saythis2,"# of clusters");
    h1_nClusters[i] = new TH1F(saythis1, saythis2, nClustersbins, nClusterslow, nClustersup);
    h1_nClusters[i]->GetXaxis()->SetTitle("number of clusters/evt");
    h1_nClusters[i]->GetYaxis()->SetTitle("counts");
    h1_nClusters[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=nClustersbins;
  }
  
  Int_t nZvertexbins = 501;
  Float_t Zvertexlow = -50.0, Zvertexup = 50.0;
  h1_zvtx = new TH1F("h1_zvtx", "# of clusters", nZvertexbins, Zvertexlow, Zvertexup);
  h1_zvtx->GetXaxis()->SetTitle("z_{vertex}");
  h1_zvtx->GetYaxis()->SetTitle("counts");
  h1_zvtx->SetMarkerStyle(kFullCircle);
  TotalNBins+=nZvertexbins;

  h1_trigger = new TH1F("h1_trigger", "trigger number returned", 1001,-0.5,1000.5);
  TotalNBins+=1001;

  h1_centrality = new TH1F("h1_centrality", "centrality", 1001,-0.1,100.1);
  TotalNBins+=1001;

  Int_t Mbins = 3000;
  Float_t Mlow = 0.0, Mup = 3.0;
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_M_%d",i);
    sprintf(saythis2,"Invariant Mass");
    h1_M[i] = new TH1F(saythis1, saythis2, Mbins, Mlow, Mup);
    h1_M[i]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    h1_M[i]->GetYaxis()->SetTitle("counts");
    h1_M[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=Mbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_M_mix_%d",i);
    sprintf(saythis2,"Invariant Mass (mixed events)");
    h1_M_mix[i] = new TH1F(saythis1, saythis2, Mbins, Mlow, Mup);
    h1_M_mix[i]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    h1_M_mix[i]->GetYaxis()->SetTitle("counts");
    h1_M_mix[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=Mbins;
  }

  Int_t ptbins = 2000;
  Float_t ptlow = 0.0, ptup = 20.0;  
  Int_t Ebins = 500;
  Float_t Elow = 0.0, Eup = 20.0;
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_E_%d",i);
    sprintf(saythis2,"Cluster Energy in EMCal");
    h1_E[i] = new TH1F(saythis1, saythis2, Ebins, Elow, Eup);
    h1_E[i]->GetXaxis()->SetTitle("E [GeV]");
    h1_E[i]->GetYaxis()->SetTitle("counts");
    h1_E[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=Ebins;
  }

  h2_PhiEtaCluster = new TH2F("h2_PhiEtaCluster", "cluster phi vs eta", 400,1.362,3.178, 300,-0.728,0.728);
  h2_PhiEtaCluster->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhiEtaCluster->GetYaxis()->SetTitle("#eta");
  h2_PhiEtaCluster->GetZaxis()->SetTitle("hits");
  h2_PhiEtaCluster->SetMarkerStyle(kFullCircle);
  TotalNBins+=400*300;

  h2_PhiEtaClusterCut = new TH2F("h2_PhiEtaClusterCut", "cluster phi vs eta (after cuts)", 400,1.362,3.178, 300,-0.728,0.728);
  h2_PhiEtaClusterCut->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhiEtaClusterCut->GetYaxis()->SetTitle("#eta");
  h2_PhiEtaClusterCut->GetZaxis()->SetTitle("hits");
  h2_PhiEtaClusterCut->SetMarkerStyle(kFullCircle);
  TotalNBins+=400*300;

// eta binning
  Double_t EtaBins[97] = {-0.66687,-0.653,-0.63913,-0.62526,-0.61139,-0.59752,-0.58365,-0.56978,-0.55591,-0.54204,-0.52817,-0.5143,-0.50043,-0.48656,-0.47269,-0.45882,-0.44495,-0.43108,-0.41721,-0.40334,-0.38947,-0.3756,-0.36173,-0.34786,-0.33399,-0.32012,-0.30625,-0.29238,-0.27851,-0.26464,-0.25077,-0.2369,-0.22303,-0.20916,-0.19529,-0.18142,-0.16755,-0.15368,-0.13981,-0.12594,-0.11207,-0.0982,-0.08433,-0.07046,-0.05659,-0.04272,-0.02885,-0.01498,-0.00111,0.01276,0.02663,0.0405,0.05437,0.06824,0.08211,0.09598,0.10985,0.12372,0.13759,0.15146,0.16533,0.1792,0.19307,0.20694,0.22081,0.23468,0.24855,0.26242,0.27629,0.29016,0.30403,0.3179,0.33177,0.34564,0.35951,0.37338,0.38725,0.40112,0.41499,0.42886,0.44273,0.4566,0.47047,0.48434,0.49821,0.51208,0.52595,0.53982,0.55369,0.56756,0.58143,0.5953,0.60917,0.62304,0.63691,0.65078,0.66465};
  
  // phi binning
  Double_t PhiBins[125] = {1.408,1.4215,1.435,1.4485,1.462,1.4755,1.489,1.5025,1.516,1.5295,1.543,1.5565,1.57,1.5835,1.597,1.6105,1.624,1.6375,1.651,1.6645,1.678,1.6915,1.705,1.7185,1.732, 1.758,1.7715,1.785,1.7985,1.812,1.8255,1.839,1.8525,1.866,1.8795,1.893,1.9065,1.92,1.9335,1.947,1.9605,1.974,1.9875,2.001,2.0145,2.028,2.0415,2.055,2.0685,2.082,2.108,2.1215,2.135,2.1485,2.162,2.1755,2.189,2.2025,2.216,2.2295,2.243,2.2565,2.27,2.2835,2.297,2.3105,2.324,2.3375,2.351,2.3645,2.378,2.3915,2.405,2.4185,2.432,2.456,2.4695,2.483,2.4965,2.51,2.5235,2.537,2.5505,2.564,2.5775,2.591,2.6045,2.618,2.6315,2.645,2.6585,2.672,2.6855,2.699,2.7125,2.726,2.7395,2.753,2.7665,2.78,2.804,2.8175,2.831,2.8445,2.858,2.8715,2.885,2.8985,2.912,2.9255,2.939,2.9525,2.966,2.9795,2.993,3.0065,3.02,3.0335,3.047,3.0605,3.074,3.0875,3.101,3.1145,3.128};

  h2_PhiEtaMaxCell = new TH2F("h2_PhiEtaMaxCell", "maxcell phi vs eta", 124,PhiBins, 96,EtaBins);
  h2_PhiEtaMaxCell->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhiEtaMaxCell->GetYaxis()->SetTitle("#eta");
  h2_PhiEtaMaxCell->GetZaxis()->SetTitle("hits");
  h2_PhiEtaMaxCell->SetMarkerStyle(kFullCircle);
  TotalNBins+=96*124;

  h2_PhiEtaMaxCellCut = new TH2F("h2_PhiEtaMaxCellCut", "maxcell phi vs eta (after cuts)", 124,PhiBins, 96,EtaBins);
  h2_PhiEtaMaxCellCut->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhiEtaMaxCellCut->GetYaxis()->SetTitle("#eta");
  h2_PhiEtaMaxCellCut->GetZaxis()->SetTitle("hits");
  h2_PhiEtaMaxCellCut->SetMarkerStyle(kFullCircle);
  TotalNBins+=96*124;

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_dR_ClustTrk_%d",i);
    sprintf(saythis2,"Cluster-Track matching");
    h1_dR_ClustTrk[i] = new TH1F(saythis1, saythis2, 5000, -0.01, 5);
    h1_dR_ClustTrk[i]->GetXaxis()->SetTitle("dR [sqrt(d#phi^{2}+d#eta^{2})]");
    h1_dR_ClustTrk[i]->GetYaxis()->SetTitle("N");
    h1_dR_ClustTrk[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=5000;
  }

  h2_gE_RecTruth = new TH2F("h2_gE_RecTruth", "#gamma E_{truth}/E_{clust} vs E_{clust}", Ebins,Elow,Eup, 500,0,2);
  h2_gE_RecTruth->GetXaxis()->SetTitle("E^{rec}_{clust} [GeV]");
  h2_gE_RecTruth->GetYaxis()->SetTitle("E^{rec}_{clust}/E^{truth}_{#gamma}");
  h2_gE_RecTruth->GetZaxis()->SetTitle("counts");
  h2_gE_RecTruth->SetMarkerStyle(kFullCircle);
  TotalNBins+=Ebins*500;
  
  h2_eop_E = new TH2F("h2_eop_E","E/p vs E (using built-in track matching)", Ebins, Elow, Eup, 1200,0,3);
  h2_eop_E->GetXaxis()->SetTitle("cluster Energy [GeV]");
  h2_eop_E->GetYaxis()->SetTitle("E/p");
  TotalNBins+=Ebins*1200;

  h2_eop_pT = new TH2F("h2_eop_pT","E/p vs p_{T} (using built-in track matching)", Ebins, Elow, Eup, 1200,0,3);
  h2_eop_pT->GetXaxis()->SetTitle("cluster Energy [GeV]");
  h2_eop_pT->GetYaxis()->SetTitle("E/p");
  TotalNBins+=Ebins*1200;

  h2_E_time = new TH2F("h2_E_time","cluster energy vs time", Ebins, Elow, Eup, 1000,-1e-6,1e-6);
  h2_E_time->GetXaxis()->SetTitle("cluster Energy [GeV]");
  h2_E_time->GetYaxis()->SetTitle("time [s]");
  TotalNBins+=Ebins*1000;

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_Pi0TruthPt_%d",i);
    sprintf(saythis2,"P_{T} distribution for Truth Pi0's");
    h1_Pi0TruthPt[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_Pi0TruthPt[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_Pi0TruthPt[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_Pi0TruthPt[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PriPi0TruthPt_%d",i);
    sprintf(saythis2,"P_{T} distribution for Truth Primary Pi0's");
    h1_PriPi0TruthPt[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PriPi0TruthPt[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PriPi0TruthPt[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PriPi0TruthPt[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_Pi0TruthPtEmcal_%d",i);
    sprintf(saythis2,"P_{T} distribution for Truth Pi0's (hit EMCal)");
    h1_Pi0TruthPtEmcal[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_Pi0TruthPtEmcal[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_Pi0TruthPtEmcal[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_Pi0TruthPtEmcal[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PriPi0TruthPtEmcal_%d",i);
    sprintf(saythis2,"P_{T} distribution for Truth Primary Pi0's (hit EMCal)");
    h1_PriPi0TruthPtEmcal[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PriPi0TruthPtEmcal[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PriPi0TruthPtEmcal[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PriPi0TruthPtEmcal[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_Pi0TruthPtPhi2piEta065_%d",i);
    sprintf(saythis2,"P_{T} for Truth Pi0's [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]");
    h1_Pi0TruthPtPhi2piEta065[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_Pi0TruthPtPhi2piEta065[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_Pi0TruthPtPhi2piEta065[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_Pi0TruthPtPhi2piEta065[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }
        
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_Pi0TruthPtPhi2piEta1_%d",i);
    sprintf(saythis2,"P_{T} for Truth Pi0's [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]");
    h1_Pi0TruthPtPhi2piEta1[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_Pi0TruthPtPhi2piEta1[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_Pi0TruthPtPhi2piEta1[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_Pi0TruthPtPhi2piEta1[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }
  
  h2_Pi0TruthPhiEta = new TH2F("h2_Pi0TruthPhiEta","Pi0Truth Phi vs Eta ", 380,-0.02,6.30, 200,-10,10);
  h2_Pi0TruthPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEta->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*200;

  h2_PriPi0TruthPhiEta = new TH2F("h2_PriPi0TruthPhiEta","Primary Pi0Truth Phi vs Eta ", 380,-0.02,6.30, 200,-10,10);
  h2_PriPi0TruthPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_PriPi0TruthPhiEta->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*200;

  h2_Pi0TruthPhiEtaEmcal = new TH2F("h2_Pi0TruthPhiEtaEmcal","Pi0Truth Phi vs Eta (in EMCal)", 380,-0.02,6.30, 150,-1.5,1.5);
  h2_Pi0TruthPhiEtaEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEtaEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;

  h2_PriPi0TruthPhiEtaEmcal = new TH2F("h2_PriPi0TruthPhiEtaEmcal","Primary Pi0Truth Phi vs Eta (in EMCal)", 380,-0.02,6.30, 150,-5,5);
  h2_PriPi0TruthPhiEtaEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_PriPi0TruthPhiEtaEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;

  h2_Pi0TruthPhiEta_Phi2piEta065 = new TH2F("h2_Pi0TruthPhiEta_Phi2piEta065",
					    "Pi0Truth Phi vs Eta [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]", 380,-0.02,6.30, 150,-5,5);
  h2_Pi0TruthPhiEta_Phi2piEta065->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEta_Phi2piEta065->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;

  h2_Pi0TruthPhiEta_Phi2piEta1 = new TH2F("h2_Pi0TruthPhiEta_Phi2piEta1",
					    "Pi0Truth Phi vs Eta [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]", 380,-0.02,6.30, 150,-1.5,1.5);
  h2_Pi0TruthPhiEta_Phi2piEta1->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEta_Phi2piEta1->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;
    
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_TruthPhotonsEmcal_%d",i);
    sprintf(saythis2,"P_{T} distribution for photons (in EMCal)");
    h1_TruthPhotonsEmcal[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_TruthPhotonsEmcal[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_TruthPhotonsEmcal[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_TruthPhotonsEmcal[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  h2_TruthPhotonsPhiEta = new TH2F("h2_TruthPhotonsPhiEta", 
				   "Truth Photons Phi vs Eta (pointed at emcal)", 380,-0.02,6.30, 150,-1.5,1.5);
  h2_TruthPhotonsPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_TruthPhotonsPhiEta->GetYaxis()->SetTitle("#eta ");
  h2_TruthPhotonsPhiEta->SetMarkerStyle(kFullCircle);
  TotalNBins+=380*150;

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PhotonsEmcal_%d",i);
    sprintf(saythis2,"P_{T} distribution for photons (in EMCal)");
    h1_PhotonsEmcal[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PhotonsEmcal[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PhotonsEmcal[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PhotonsEmcal[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PhotonsNCellsCut_%d",i);
    sprintf(saythis2,"P_{T} distribution for #gamma's that survive NCells cut");
    h1_PhotonsNCellsCut[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PhotonsNCellsCut[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PhotonsNCellsCut[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PhotonsNCellsCut[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PhotonsTrackMatchCut_%d",i);
    sprintf(saythis2,"P_{T} distribution for #gamma's that survive TrackMatch cut");
    h1_PhotonsTrackMatchCut[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PhotonsTrackMatchCut[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PhotonsTrackMatchCut[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PhotonsTrackMatchCut[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_PhotonsAllCut_%d",i);
    sprintf(saythis2,"P_{T} distribution for #gamma's that survive All cut");
    h1_PhotonsAllCut[i] = new TH1F(saythis1, saythis2, ptbins, ptlow, ptup);
    h1_PhotonsAllCut[i]->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    h1_PhotonsAllCut[i]->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    h1_PhotonsAllCut[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;
  }

  h2_PhotonsPhiEtaIsEmcal = new TH2F("h2_PhotonsPhiEtaIsEmcal",
				     "Photons Phi vs Eta (IsEMCAL()==1)", 380,-0.02,6.30, 100,-1.0,1.0);
  h2_PhotonsPhiEtaIsEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhotonsPhiEtaIsEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*100;
  
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_dR_RealMC_%d",i);
    sprintf(saythis2,"P_{T} distribution for #gamma's that survive All cut");
    h1_dR_RealMC[i] = new TH1F(saythis1, saythis2, 2000, -0.01, 10);
    h1_dR_RealMC[i]->GetXaxis()->SetTitle("dR sqrt(dx^{2}+dy^{2})");
    h1_dR_RealMC[i]->GetYaxis()->SetTitle("N");
    h1_dR_RealMC[i]->SetMarkerStyle(kFullCircle);
    TotalNBins+=2000;
  }

  Int_t chi2bins = 100;
  Float_t chi2low = -2, chi2up = 2;
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_Chi2_%d",i);
    sprintf(saythis2,"#chi^{2} distribution for reconstructed");
    h1_Chi2[i] = new TH1F(saythis1,saythis2,chi2bins, chi2low, chi2up);
    h1_Chi2[i]->GetXaxis()->SetTitle("#chi^{2}");
    h1_Chi2[i]->GetYaxis()->SetTitle("counts");
    TotalNBins+=chi2bins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_nTrkMatch_%d",i);
    sprintf(saythis2,"number of matched tracks");
    h1_nTrkMatch[i] = new TH1F(saythis1,saythis2,14, -1.5, 5.5);
    h1_nTrkMatch[i]->GetXaxis()->SetTitle("nTracksMatched");
    h1_nTrkMatch[i]->GetYaxis()->SetTitle("counts");
    TotalNBins+=14;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h1_ClusterDisp_%d",i);
    sprintf(saythis2,"Dispersion of CaloCluster");
    h1_ClusterDisp[i] = new TH1F(saythis1,saythis2,1000, -1, 3);
    h1_ClusterDisp[i]->GetXaxis()->SetTitle("cluster->GetClusterDisp()");
    h1_ClusterDisp[i]->GetYaxis()->SetTitle("counts");
    TotalNBins+=1000;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_Ellipse_%d",i);
    sprintf(saythis2,"Ellipse axis M20 vs M02");
    h2_Ellipse[i] = new TH2F(saythis1,saythis2,500, -0.01, 1, 500, -0.01, 1);
    h2_Ellipse[i]->GetXaxis()->SetTitle("cluster->GetM20()");
    h2_Ellipse[i]->GetYaxis()->SetTitle("cluster->GetM02()");
    h2_Ellipse[i]->GetZaxis()->SetTitle("counts");
    TotalNBins+=500*500;
  }
       
  Int_t etabins = 150;
  Float_t etalow = -1.5, etaup = 1.5;
  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_EtaPt_%d",i);
    sprintf(saythis2,"Cluster Energy vs ");
    h2_EtaPt[i] = new TH2F(saythis1,saythis2,etabins, etalow, etaup, ptbins, ptlow, ptup);
    h2_EtaPt[i]->GetXaxis()->SetTitle("E [GeV]");
    h2_EtaPt[i]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    TotalNBins+=etabins*ptbins;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h3_MptAsymm_%d",i);
    sprintf(saythis2,"mass vs p_{T} vs Asymm cut");
    h3_MptAsymm[i] = new TH3F(saythis1,saythis2,Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
    h3_MptAsymm[i]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    h3_MptAsymm[i]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    h3_MptAsymm[i]->GetZaxis()->SetTitle("Asymmetry Cut (edges: 0.0, 0.1, 0.7, 1.0)");
    TotalNBins+=Mbins*ptbins*3.0;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h3_MptAsymm_mix_%d",i);
    sprintf(saythis2,"mass vs p_{T} vs Asymm cut (mixed events)");
    h3_MptAsymm_mix[i] = new TH3F(saythis1,saythis2,Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
    h3_MptAsymm_mix[i]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
    h3_MptAsymm_mix[i]->GetYaxis()->SetTitle("p_{T} [GeV/c]");
    h3_MptAsymm_mix[i]->GetZaxis()->SetTitle("Asymmetry Cut (edges: 0.0, 0.1, 0.7, 1.0)");
    TotalNBins+=Mbins*ptbins*3.0;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_dphi_deta_%d",i);
    sprintf(saythis2,"#Delta#phi vs #Delta#eta");
    h2_dphi_deta[i] = new TH2F(saythis1,saythis2, 349,-1.5,5, 400,-2.0,2.0);
    h2_dphi_deta[i]->GetXaxis()->SetTitle("#Delta#phi");
    h2_dphi_deta[i]->GetYaxis()->SetTitle("#Delta#eta");
    TotalNBins+=349*400;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_dphi_deta_mix_%d",i);
    sprintf(saythis2,"#Delta#phi vs #Delta#eta (mixed events)");
    h2_dphi_deta_mix[i] = new TH2F(saythis1,saythis2, 349,-1.5,5, 400,-2.0,2.0);
    h2_dphi_deta_mix[i]->GetXaxis()->SetTitle("#Delta#phi");
    h2_dphi_deta_mix[i]->GetYaxis()->SetTitle("#Delta#eta");
    TotalNBins+=349*400;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_DispRes_%d",i);
    sprintf(saythis2,"zvtx info");
    h2_DispRes[i] = new TH2F(saythis1, saythis2, 500,-0.01,1, 500,-0.1,2);
    h2_DispRes[i]->GetXaxis()->SetTitle("EvtVtx->GetDispersion()");
    h2_DispRes[i]->GetYaxis()->SetTitle("EvtVtx->GetZRes()");
    h2_DispRes[i]->GetZaxis()->SetTitle("counts");
    TotalNBins+=500*500;
  }

  for(int i=0; i<cent_bins; i++){
    sprintf(saythis1,"h2_cells_M02_%d",i);
    sprintf(saythis2,"nCells vs M02");
    h2_cells_M02[i] = new TH2F(saythis1, saythis2, 204,-1.5,100.5, 500,-1,1.5);
    h2_cells_M02[i]->GetXaxis()->SetTitle("nCells");
    h2_cells_M02[i]->GetYaxis()->SetTitle("M02");
    h2_cells_M02[i]->GetZaxis()->SetTitle("counts");
    TotalNBins+=204*500;
  }

  cout << endl << "Total number of bins in booked histograms:  " << TotalNBins << endl << endl;

  // Initialize helper class (for vertex selection & pile up correction)
  fHelperClass = new AliAnalysisUtils();

  //TFile *f = OpenFile(1); 
  //TDirectory::TContext context(f);
    
  fOutput->Add(h1_zvtx);
  fOutput->Add(h1_trigger);
  fOutput->Add(h1_centrality);
  fOutput->Add(h2_PhiEtaCluster);
  fOutput->Add(h2_PhiEtaClusterCut);
  fOutput->Add(h2_PhiEtaMaxCell);
  fOutput->Add(h2_PhiEtaMaxCellCut);
  fOutput->Add(h2_gE_RecTruth);
  fOutput->Add(h2_eop_E);
  fOutput->Add(h2_eop_pT);
  fOutput->Add(h2_E_time);

  for(int i=0; i<cent_bins; i++){
    fOutput->Add(h1_nClusters[i]);
    fOutput->Add(h1_M[i]);
    fOutput->Add(h1_M_mix[i]);
    fOutput->Add(h1_E[i]);
    fOutput->Add(h1_dR_ClustTrk[i]);
    fOutput->Add(h1_Pi0TruthPt[i]);
    fOutput->Add(h1_PriPi0TruthPt[i]);
    fOutput->Add(h1_Pi0TruthPtEmcal[i]);
    fOutput->Add(h1_PriPi0TruthPtEmcal[i]);
    fOutput->Add(h1_Pi0TruthPtPhi2piEta065[i]);
    fOutput->Add(h1_Pi0TruthPtPhi2piEta1[i]);
    fOutput->Add(h1_TruthPhotonsEmcal[i]);
    fOutput->Add(h1_PhotonsEmcal[i]);
    fOutput->Add(h1_PhotonsNCellsCut[i]);
    fOutput->Add(h1_PhotonsTrackMatchCut[i]);
    fOutput->Add(h1_PhotonsAllCut[i]);
    fOutput->Add(h1_dR_RealMC[i]);
    fOutput->Add(h1_Chi2[i]);
    fOutput->Add(h1_nTrkMatch[i]);
    fOutput->Add(h1_ClusterDisp[i]);
    fOutput->Add(h2_Ellipse[i]);
    fOutput->Add(h2_EtaPt[i]);
    fOutput->Add(h3_MptAsymm[i]);
    fOutput->Add(h3_MptAsymm_mix[i]);
    fOutput->Add(h2_dphi_deta[i]);
    fOutput->Add(h2_dphi_deta_mix[i]);
    fOutput->Add(h2_DispRes[i]);
    fOutput->Add(h2_cells_M02[i]);
  }
  fOutput->Add(h2_Pi0TruthPhiEta);
  fOutput->Add(h2_PriPi0TruthPhiEta);
  fOutput->Add(h2_Pi0TruthPhiEtaEmcal);
  fOutput->Add(h2_PriPi0TruthPhiEtaEmcal);
  fOutput->Add(h2_Pi0TruthPhiEta_Phi2piEta065);
  fOutput->Add(h2_Pi0TruthPhiEta_Phi2piEta1);
  fOutput->Add(h2_TruthPhotonsPhiEta);
  fOutput->Add(h2_PhotonsPhiEtaIsEmcal);

  // Post data for ALL output slots >0 here, 
  // To get at least an empty histogram 
  // 1 is the outputnumber of a certain weg of task 1  
  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALMesonGGSDMpPb::UserExec(Option_t *) 
{
  // Main loop Called for each event

  AliMCEvent *mcEvent = MCEvent();  
  Bool_t isMC = bool(mcEvent);//is this the right way to do this? 
  
  TRandom3 randy; randy.SetSeed(0);
  unsigned int iskip = -1;
  TLorentzVector ParentMix;

  double recalScale = 1.0;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();    
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (am->GetInputEventHandler());
  AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*> (am->GetInputEventHandler());
  if (!aodH && !esdH)  Printf("ERROR: Could not get ESD or AODInputHandler");
  
  if(esdH)      fEsdEv = (AliESDEvent*)esdH->GetEvent();    
  else if(aodH) fAodEv = aodH->GetEvent();  
  else{
    AliFatal("Neither ESD nor AOD event found");
    return;
  }

  // get pointer to reconstructed event
  AliVEvent *event = InputEvent();
  if (!event){
    AliError("Pointer == 0, this can not happen!");  return;}
  //AliESDEvent* fEsdEv = dynamic_cast<AliESDEvent*>(event);
  //AliAODEvent* aod = dynamic_cast<AliAODEvent*>(event);
  //if (!fEsdEv){
  //AliError("Cannot get the ESD event");  return;}

  fHelperClass->SetCutOnZVertexSPD(kFALSE);//does the zvtx have to match the spd vertex? 
  fHelperClass->SetMaxVtxZ(1.0e6);//i set this myself later.. 
  // simply makes sure that there is at least 1 contributer to the zvtx determination.
  // this should only remove the *extra* events at zvtx==0.
  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;

  Int_t iTrigger = 0;
  if (fEsdEv)       iTrigger = fEsdEv->GetHeader()->GetL0TriggerInputs();
  else if (fAodEv)  iTrigger = fAodEv->GetHeader()->GetL0TriggerInputs();
  
  char saythis[500];
  Int_t iTriggerBin = 0;
  for(unsigned long j=0; j<TriggerList.size(); j++){
    if(iTrigger==TriggerList[j])
      iTriggerBin=j+1;
  }
  if(iTriggerBin==0){
    TriggerList.push_back(iTrigger);
    iTriggerBin=TriggerList.size();
  }
  
  h1_trigger->SetBinContent(iTriggerBin, h1_trigger->GetBinContent(iTriggerBin)+1);
  sprintf(saythis,"%d",iTrigger);
  h1_trigger->GetXaxis()->SetBinLabel(iTriggerBin, saythis);
  
  Double_t centralityVZERO=0.0;
  Int_t centBin = 0;
  AliCentrality *aliCent=NULL;

  Int_t nclusters=0;
  if(fEsdEv){
    //Int_t evtN      = fEsdEv->GetEventNumberInFile();  
    //Int_t ntracks   = fEsdEv->GetNumberOfTracks();
    nclusters = fEsdEv->GetNumberOfCaloClusters();
    aliCent   = fEsdEv->GetCentrality();
  }
  else if(fAodEv){
    //Int_t evtN      = fAodEv->GetEventNumberInFile();  
    //Int_t ntracks   = fAodEv->GetNumberOfTracks();
    nclusters = fAodEv->GetNumberOfCaloClusters();
    aliCent   = fAodEv->GetCentrality();
  }
  //centBin = aliCent->GetCentralityClass10("V0M");
  //centralityVZERO = aliCent->GetCentralityPercentile("V0M");
  //centBin = aliCent->GetCentralityClass10("V0C");
  //centralityVZERO = aliCent->GetCentralityPercentile("V0C");
  centBin = aliCent->GetCentralityClass10("V0A");
  centralityVZERO = aliCent->GetCentralityPercentile("V0A");

  if     (centralityVZERO<20.0)
    centBin = 0;
  else if(centralityVZERO<40.0)
    centBin = 1;
  else if(centralityVZERO<60.0)
    centBin = 2;
  else
    centBin = 3;

  //cout << "Centrality:  " << centBin << "    " << centralityVZERO << endl;
  
  if (fEsdEv){
    if(!(fEsdEv->GetPrimaryVertex()->GetStatus()))   return;
  }
  //else if (fAodEv){
  //if(!(fAodEv->GetPrimaryVertex()->GetStatus()))   return;
  //}

  Double_t vertDisp=0.0;
  Double_t vertZres=0.0;
  Bool_t vertIsfromZ=0;
  if (fEsdEv){
    vertDisp    = fEsdEv->GetPrimaryVertex()->GetDispersion();
    vertZres    = fEsdEv->GetPrimaryVertex()->GetZRes();
    vertIsfromZ = fEsdEv->GetPrimaryVertex()->IsFromVertexerZ();
  }
  else if (fAodEv){
    vertDisp    = 0;
    vertZres    = 0;
    vertIsfromZ = 0;
  }

  h2_DispRes[centBin]->Fill(vertDisp, vertZres);  
  // if vertex is from spd vertexZ, require more stringent cut
  if (vertIsfromZ) {
    if (vertDisp>0.02 ||  vertZres>0.25 ) 
      return; // bad vertex from VertexerZ
  }
  
  // EMCal cluster loop for reconstructed event
  //numberofclusters set above! 
  TLorentzVector Photon1, Photon2, Parent;
  Double_t vertex[3]; 
  Double_t E1=0.0;
  Double_t vertZ=0.0;
  if (fEsdEv)       vertZ = fEsdEv->GetPrimaryVertex()->GetZ();
  else if (fAodEv)  vertZ = fAodEv->GetPrimaryVertex()->GetZ();    
  
  h1_zvtx->Fill(vertZ);
  //zvertex cut:
  if(fabs(vertZ)>10.0)
    return;
  
  h1_nClusters[centBin]->Fill(nclusters);
  h1_centrality->Fill(centralityVZERO);

  int izvtx = GetZvtxBin(vertZ);
  int imult = GetMultBin(nclusters);

  //cout << iskip << " " << izvtx << " " << imult << endl;  
  //cout << "GetNumberOfVertices(): " << fAodEv->GetNumberOfVertices() << endl;



  //######################### ~~~~~~~~~~~ ##################################
  //######################### STARTING MC ##################################
  //######################### ~~~~~~~~~~~ ##################################
  
  if(isMC){
    int isPrimary = 0;

    if (!mcEvent){
      cout << "no MC event" << endl;
      return;
    }
    
    const AliVVertex *evtVtx = mcEvent->GetPrimaryVertex();
    if (!evtVtx)
      return;
    
    mcEvent->PreReadAll();    
    
    Int_t nTracksMC  = mcEvent->GetNumberOfTracks();
    Int_t nPTracksMC = mcEvent->GetNumberOfPrimaries();
    //cout << "We have  " << nPTracksMC << "  primaries of  " << nTracksMC << "  total tracks." << endl;
    
    for (Int_t iTrack = 0; iTrack<nTracksMC; ++iTrack) {
      AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
      if (!mcP)
	continue;
      

      // it's a pion !! 
      if(mcP->PdgCode() != 111)
	continue;
      
      /*
      // primary particle
      Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
                                (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
      if(dR <= 0.01)  isPrimary = 1;
      else            isPrimary = 0;
      */
      
      if(iTrack<nPTracksMC)  isPrimary = 1;
      else                   isPrimary = 0;
            
      h1_Pi0TruthPt    [centBin]->Fill(mcP->Pt());
      h2_Pi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());

      if(isPrimary==1){
	h1_PriPi0TruthPt    [centBin]->Fill(mcP->Pt());
	h2_PriPi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());
      }
      
      if(mcP->Eta()<-1.0 || mcP->Eta()>1.0)
	continue;
      
      h1_Pi0TruthPtPhi2piEta1    [centBin]->Fill(mcP->Pt());
      h2_Pi0TruthPhiEta_Phi2piEta1->Fill(mcP->Phi(),mcP->Eta());      
      
      if(mcP->Eta()>fEtamin && mcP->Eta()<fEtamax){
	h1_Pi0TruthPtPhi2piEta065    [centBin]->Fill(mcP->Pt());
	h2_Pi0TruthPhiEta_Phi2piEta065->Fill(mcP->Phi(),mcP->Eta());      	
      }
      
      
      Int_t d1 = mcP->GetDaughterFirst();
      Int_t d2 = mcP->GetDaughterLast();
      
      if (d1<0)  continue;
      if (d2<0)  d2=d1;      
      if (d2-d1 != 1)  continue;
      
      bool bacc = true;
      bool binp = true;
      for (Int_t i=d1;i<=d2;++i){
        const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(i));
        Double_t eta_d = dmc->Eta();
        Double_t phi_d = dmc->Phi();
        if(!(dmc->PdgCode()==22)){
	  binp = false;
        }
        if(!(dmc->PdgCode()==22 && eta_d>fEtamin && eta_d<fEtamax && phi_d>fPhimin && phi_d<fPhimax)){
	  bacc = false;
        }	
      }

      if(binp && bacc){// 2 Photons hit the EMCAL! 
	
	for (Int_t j=d1;j<=d2;++j){//both truth photons.
	  
	  const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(j));
	  Double_t eta_d = dmc->Eta();
	  Double_t phi_d = dmc->Phi();
	  
	  if( dmc->PdgCode()==22 && 
	      dmc->Eta()>fEtamin && dmc->Eta()<fEtamax && 
	      dmc->Phi()>fPhimin && dmc->Phi()<fPhimax ){
	    h1_TruthPhotonsEmcal[centBin]->Fill(dmc->Pt());
	    h2_TruthPhotonsPhiEta->Fill(dmc->Phi(),dmc->Eta());
	  }

	  for(int i=0; i<nclusters; i++) {
	    
	    Bool_t matches_pion_photon = 0;
	    
	    AliESDCaloCluster* esdCluster=NULL;
	    AliAODCaloCluster* aodCluster=NULL;
	    if (fEsdEv)       esdCluster = fEsdEv->GetCaloCluster(i); // pointer to EMCal cluster
	    else if (fAodEv)  aodCluster = fAodEv->GetCaloCluster(i); // pointer to EMCal cluster
	    
	    Double_t clustMC_phi, clustMC_eta;
	    
	    if(fEsdEv){
	      
	      if(esdCluster->IsEMCAL()){
		
		Float_t pos[3] = {0,0,0};
		esdCluster->GetPosition(pos);
		TVector3 vpos(pos);
		//h1_Phi->Fill(vpos.Phi());
		clustMC_phi = vpos.Phi();
		clustMC_eta = vpos.Eta();
		
		Double_t dR = TMath::Sqrt((eta_d-clustMC_eta)*(eta_d-clustMC_eta) + 
					  (phi_d-clustMC_phi)*(phi_d-clustMC_phi));
		h1_dR_RealMC[centBin]->Fill(dR);
		if(dR<=0.04) matches_pion_photon = 1;
		
		vpos.Delete();
	      }
	      if(matches_pion_photon){		
		if(esdCluster->IsEMCAL()){
		  h1_PhotonsEmcal[centBin]->Fill(esdCluster->E());
		  h2_PhotonsPhiEtaIsEmcal->Fill(clustMC_phi,clustMC_eta);
		}
		if(esdCluster->IsEMCAL() && esdCluster->GetNCells()>=2)
		  h1_PhotonsNCellsCut[centBin]->Fill(esdCluster->E());
		if(esdCluster->IsEMCAL() && esdCluster->GetNTracksMatched()==0)
		  h1_PhotonsTrackMatchCut[centBin]->Fill(esdCluster->E());
		if(esdCluster->IsEMCAL() && esdCluster->GetNCells()>=2 && esdCluster->GetNTracksMatched()==0)
		  h1_PhotonsAllCut[centBin]->Fill(esdCluster->E());		  
	      }//if(matches_pion_photon)
		
	    }//if(fEsdEv)
	    else if(fAodEv){
	      
	      if(aodCluster->IsEMCAL()){
		
		Float_t pos[3] = {0,0,0};
		aodCluster->GetPosition(pos);  
		TVector3 vpos(pos); 
		//h1_Phi->Fill(vpos.Phi());
		clustMC_phi = vpos.Phi();
		clustMC_eta = vpos.Eta();
		
		Double_t dR = TMath::Sqrt((eta_d-clustMC_eta)*(eta_d-clustMC_eta) + 
					  (phi_d-clustMC_phi)*(phi_d-clustMC_phi));
		h1_dR_RealMC[centBin]->Fill(dR);
		if(dR<=0.04) matches_pion_photon = 1;
		
		vpos.Delete();
	      }
	      if(matches_pion_photon){		
		if(aodCluster->IsEMCAL()){
		  h1_PhotonsEmcal[centBin]->Fill(aodCluster->E());
		  h2_PhotonsPhiEtaIsEmcal->Fill(clustMC_phi,clustMC_eta);
		}
		if(aodCluster->IsEMCAL() && aodCluster->GetNCells()>=2)
		  h1_PhotonsNCellsCut[centBin]->Fill(aodCluster->E());
		if(aodCluster->IsEMCAL() && aodCluster->GetNTracksMatched()==0)
		  h1_PhotonsTrackMatchCut[centBin]->Fill(aodCluster->E());
		if(aodCluster->IsEMCAL() && aodCluster->GetNCells()>=2 && aodCluster->GetNTracksMatched()==0)
		  h1_PhotonsAllCut[centBin]->Fill(aodCluster->E());
		
	      }//if(matches_pion_photon)
	      
	    }//if(fAodEv)
	    
	  }//loop over nclusters. 
	  
	}//both truth photons.
	
      }// 2 Photons hit the EMCAL! 
      
      
      if(binp && bacc){// 2 Photons hit the EMCAL! 
	h1_Pi0TruthPtEmcal    [centBin]->Fill(mcP->Pt());
	h2_Pi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());	
	
	if(isPrimary==1){
	  h1_PriPi0TruthPtEmcal    [centBin]->Fill(mcP->Pt());
	  h2_PriPi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());
	}
	
      }//2 photons hit the EMCAL! 
      
    }//for(nTracksMC)    
    
  }//if(isMC)
  
  //######################### ~~~~~~~~~~~~ ##################################
  //######################### DONE WITH MC ##################################
  //######################### ~~~~~~~~~~~~ ##################################


  for(int i=0; i<nclusters; i++) {

    AliESDCaloCluster* esdCluster=NULL;
    AliAODCaloCluster* aodCluster=NULL;
    if (fEsdEv)       esdCluster = fEsdEv->GetCaloCluster(i); // pointer to EMCal cluster
    else if (fAodEv)  aodCluster = fAodEv->GetCaloCluster(i); // pointer to EMCal cluster
    if(!esdCluster && !aodCluster) { 
      AliError(Form("ERROR: Could not retrieve any (ESD or AOD) Cluster %d",i)); 
      continue; 
    }
    
    if(fEsdEv){

      recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
      
      //uncomment this to do the track matching (1 of 3 lines, esd part)!! 
      //Bool_t MatchesToTrack = 0;
      if(esdCluster->IsEMCAL()){
	
	Float_t pos[3] = {0,0,0};
	Short_t maxCellID = -1;
	Float_t celleta, cellphi;
	esdCluster->GetPosition(pos);
	TVector3 clusterPosition(pos);
	h2_PhiEtaCluster->Fill(clusterPosition.Phi(),clusterPosition.Eta());
	GetMaxCellEnergy(esdCluster, maxCellID);
	AliEMCALGeometry *fGeom = AliEMCALGeometry::GetInstance();
	fGeom->EtaPhiFromIndex(maxCellID,celleta,cellphi);
	h2_PhiEtaMaxCell->Fill(cellphi,celleta);
	
	// _______________Track loop for reconstructed event_____________
	for(Int_t itrk = 0; itrk < fEsdEv->GetNumberOfTracks(); itrk++) {
	  AliESDtrack* esdTrack = fEsdEv->GetTrack(itrk); // pointer to reconstructed to track
	  if(!esdTrack) { 
	    AliError(Form("ERROR: Could not retrieve any (ESD) track %d",itrk)); 
	    continue; 
	  }
	  
	  Double_t posTrk[3] = {0,0,0};
	  esdTrack->GetXYZ(posTrk);
	  TVector3 vposTrk(posTrk);
	  
	  Double_t fMass          = 0.139;
	  Double_t fStepSurface   = 20.;
	  Float_t etaproj, phiproj, pttrackproj;

	  AliExternalTrackParam *trackParam =  const_cast<AliExternalTrackParam*>(esdTrack->GetInnerParam());
	  if(!trackParam) continue;
	  AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(trackParam, 440., fMass, fStepSurface, etaproj, phiproj, pttrackproj);
	  
	  double dR_clusttrk = sqrt((phiproj-clusterPosition.Phi())*(phiproj-clusterPosition.Phi()) + 
				    (etaproj-clusterPosition.Eta())*(etaproj-clusterPosition.Eta()) );
	  
	  h1_dR_ClustTrk[centBin]->Fill(dR_clusttrk);
	  
	  //uncomment this to do the track matching (2 of 3 lines)!! 
	  //if(dR_clusttrk<fdRmin_ClustTrack)
	  //MatchesToTrack = 1;
	  
	}//_____________________________nTracks__________________________
		
	h2_cells_M02  [centBin]->Fill(esdCluster->GetNCells(),esdCluster->GetM02());
	h2_Ellipse    [centBin]->Fill(esdCluster->GetM20(),esdCluster->GetM02());
	h1_Chi2       [centBin]->Fill(esdCluster->Chi2());//always -1. 
	h1_nTrkMatch  [centBin]->Fill(esdCluster->GetNTracksMatched());
	h1_ClusterDisp[centBin]->Fill(esdCluster->GetDispersion());
	h2_E_time              ->Fill(esdCluster->E(),esdCluster->GetTOF());

	TArrayI *TrackLabels = esdCluster->GetTracksMatched();
	if(TrackLabels){
	  if(TrackLabels->GetSize()>0){
	    Int_t trackindex = TrackLabels->At(0);
	    AliESDtrack* matchingT = fEsdEv->GetTrack(trackindex); // pointer to reconstructed to track
	  
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h2_eop_E ->Fill(esdCluster->E()*recalScale, esdCluster->E()*recalScale/matchingT->P());
	    h2_eop_pT->Fill(matchingT->Pt(),            esdCluster->E()*recalScale/matchingT->P());
	  }
	}

	//uncomment this to do the track matching (3 of 3 lines)!! 
	//if(isGoodEsdCluster(esdCluster) && !MatchesToTrack){
	if(isGoodEsdCluster(esdCluster)){
	  recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	  E1 = esdCluster->E()*recalScale;// TOTAL HACK - JJ
	  fEsdEv->GetVertex()->GetXYZ(vertex);
	  esdCluster->GetMomentum(Photon1,vertex);
	  Photon1.SetPx(Photon1.Px()*recalScale);// TOTAL HACK - JJ
	  Photon1.SetPy(Photon1.Py()*recalScale);// TOTAL HACK - JJ
	  Photon1.SetPz(Photon1.Pz()*recalScale);// TOTAL HACK - JJ
	  Photons[0][izvtx][imult].push_back( TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) );
	  h1_E[centBin]->Fill(E1);
	  h2_PhiEtaClusterCut->Fill(clusterPosition.Phi(),clusterPosition.Eta());
	  h2_PhiEtaMaxCellCut->Fill(cellphi,celleta);
	}
	clusterPosition.Delete();	
      }//if(esdCluster->isEMCAL())
    }//if(fEsdEv)
    else if(fAodEv){
      
      recalScale = PrivateEnergyRecal(aodCluster->E(), fRecalibrator);
      
      //uncomment this to do the track matching (1 of 3 lines, aod part)!! 
      //Bool_t MatchesToTrack = 0;
      if(aodCluster->IsEMCAL()){

	Float_t pos[3] = {0,0,0};
	Short_t maxCellID = -1;
	Float_t celleta, cellphi;
	aodCluster->GetPosition(pos);  
	TVector3 clusterPosition(pos); 
	h2_PhiEtaCluster->Fill(clusterPosition.Phi(),clusterPosition.Eta());
	GetMaxCellEnergy(aodCluster, maxCellID);
	AliEMCALGeometry *fGeom = AliEMCALGeometry::GetInstance();
	fGeom->EtaPhiFromIndex(maxCellID,celleta,cellphi);
	h2_PhiEtaMaxCell->Fill(cellphi,celleta);

	// _______________Track loop for reconstructed event_____________
	for(Int_t itrk = 0; itrk < fAodEv->GetNumberOfTracks(); itrk++) {
	  AliAODTrack* aodTrack = dynamic_cast<AliAODTrack*>(fAodEv->GetTrack(itrk));
	  if(!aodTrack) AliFatal("Not a standard AOD"); // pointer to reconstructed to track
	  if(!aodTrack) { 
	    AliError(Form("ERROR: Could not retrieve any (AOD) track %d",itrk)); 
	    continue; 
	  }

	  Double_t posTrk[3] = {0,0,0};
	  aodTrack->GetXYZ(posTrk);
	  TVector3 vposTrk(posTrk);
	  
	  Double_t fMass          = 0.139;
	  Double_t fStepSurface   = 20.;
	  Float_t etaproj, phiproj, pttrackproj;
	  
	  AliExternalTrackParam *trackParam =  const_cast<AliExternalTrackParam*>(aodTrack->GetInnerParam());
	  if(!trackParam) continue;
	  AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(trackParam, 440., fMass, fStepSurface, etaproj, phiproj, pttrackproj);
	  
	  double dR_clusttrk = sqrt((phiproj-clusterPosition.Phi())*(phiproj-clusterPosition.Phi()) + 
				    (etaproj-clusterPosition.Eta())*(etaproj-clusterPosition.Eta()) );
	  
	  h1_dR_ClustTrk[centBin]->Fill(dR_clusttrk);
	  
	  //uncomment this to do the track matching (2 of 3 lines, aod part)!! 
	  //if(dR_clusttrk<fdRmin_ClustTrack)
	  //MatchesToTrack = 1;


	}//_____________________________nTracks__________________________

	h2_cells_M02  [centBin]->Fill(aodCluster->GetNCells(),aodCluster->GetM02());
	h2_Ellipse    [centBin]->Fill(aodCluster->GetM20(),aodCluster->GetM02());
	h1_Chi2       [centBin]->Fill(aodCluster->Chi2());//always -1. 
	h1_nTrkMatch  [centBin]->Fill(aodCluster->GetNTracksMatched());
	h1_ClusterDisp[centBin]->Fill(aodCluster->GetDispersion());
	h2_E_time              ->Fill(aodCluster->E(),aodCluster->GetTOF());

	// #################################################
	// track matching eop histograms are handled here... 
	// #################################################
      
	//uncomment this to do the track matching (3 of 3 lines, aod part)!! 
	//if(isGoodAodCluster(aodCluster) && !MatchesToTrack){
	if(isGoodAodCluster(aodCluster)){
	  recalScale = PrivateEnergyRecal(aodCluster->E(), fRecalibrator);
	  E1 = aodCluster->E()*recalScale;// TOTAL HACK - JJ
	  fAodEv->GetVertex(0)->GetXYZ(vertex);
	  aodCluster->GetMomentum(Photon1,vertex);
	  Photon1.SetPx(Photon1.Px()*recalScale);// TOTAL HACK - JJ
	  Photon1.SetPy(Photon1.Py()*recalScale);// TOTAL HACK - JJ
	  Photon1.SetPz(Photon1.Pz()*recalScale);// TOTAL HACK - JJ
	  Photons[0][izvtx][imult].push_back( TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) );
	  h1_E[centBin]->Fill(E1);
	  h2_PhiEtaClusterCut->Fill(clusterPosition.Phi(),clusterPosition.Eta());
	  h2_PhiEtaMaxCellCut->Fill(cellphi,celleta);	  
	}
	clusterPosition.Delete();
      }//if(aodCluster->IsEMCAL())
    }//if(fAodEv)
    
  }//loop over nclusters. 
  
  //Make same event pions... 
  for(unsigned int i=0; i<Photons[0][izvtx][imult].size(); i++){
    for(unsigned int j=i+1; j<Photons[0][izvtx][imult].size(); j++){
      Parent = Photons[0][izvtx][imult][i] + Photons[0][izvtx][imult][j];
      Double_t deltaphi = getDeltaPhi(Photons[0][izvtx][imult][i],Photons[0][izvtx][imult][j]);
      Double_t deltaeta = getDeltaEta(Photons[0][izvtx][imult][i],Photons[0][izvtx][imult][j]);
      Double_t pairasym = fabs(Photons[0][izvtx][imult][i].Pt()-Photons[0][izvtx][imult][j].Pt())/
	                      (Photons[0][izvtx][imult][i].Pt()+Photons[0][izvtx][imult][j].Pt());
      Int_t asymCut = 0;
      if     (pairasym<0.1)  asymCut = 1;
      else if(pairasym<0.7)  asymCut = 2;
      else                   asymCut = 3;
      
      h1_M        [centBin]->Fill(Parent.M());
      h3_MptAsymm [centBin]->Fill(Parent.M(),Parent.Pt(),asymCut);
      h2_dphi_deta[centBin]->Fill(deltaphi,deltaeta);
    }
  }
  
  //Make mixed event...
  for(unsigned int i=0; i<Photons[0][izvtx][imult].size(); i++){
    for(unsigned int ipool=1; ipool<poolDepth; ipool++){
      for(unsigned int j=0; j<Photons[ipool][izvtx][imult].size(); j++){
	iskip = randy.Integer(Photons[0][izvtx][imult].size());
	if(j==iskip) continue;
	Parent = Photons[0][izvtx][imult][i]+Photons[ipool][izvtx][imult][j];
	Double_t deltaphi = getDeltaPhi(Photons[0][izvtx][imult][i],Photons[ipool][izvtx][imult][j]);
	Double_t deltaeta = getDeltaEta(Photons[0][izvtx][imult][i],Photons[ipool][izvtx][imult][j]);
	Double_t pairasym = fabs(Photons[0][izvtx][imult][i].Pt()-Photons[ipool][izvtx][imult][j].Pt())/
	                        (Photons[0][izvtx][imult][i].Pt()+Photons[ipool][izvtx][imult][j].Pt());
	Int_t asymCut = 0;
	if     (pairasym<0.1)  asymCut = 1;
	else if(pairasym<0.7)  asymCut = 2;
	else                   asymCut = 3;

	h1_M_mix        [centBin]->Fill(Parent.M());
	h3_MptAsymm_mix [centBin]->Fill(Parent.M(),Parent.Pt(),asymCut);
	h2_dphi_deta_mix[centBin]->Fill(deltaphi,deltaeta);
      }
    }
  } 
    
  for(int ipool=poolDepth-1; ipool>0; ipool--){
    Photons[ipool][izvtx][imult].clear();
    for(unsigned int i=0; i<Photons[ipool-1][izvtx][imult].size(); i++)
      Photons[ipool][izvtx][imult].push_back(Photons[ipool-1][izvtx][imult][i]);     
  }
  Photons[0][izvtx][imult].clear();
    

  
  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
  }

//________________________________________________________________________
void AliAnalysisTaskEMCALMesonGGSDMpPb::Terminate(Option_t *) //specify what you want to have done
{
  // Called once at the end of the query.
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::GetZvtxBin(Double_t vertZ)
{
  
  int izvtx = -1;
  
  if     (vertZ<-35)
    izvtx=0;
  else if(vertZ<-30)
    izvtx=1;
  else if(vertZ<-25)
    izvtx=2;
  else if(vertZ<-20)
    izvtx=3;
  else if(vertZ<-15)
    izvtx=4;
  else if(vertZ<-10)
    izvtx=5;
  else if(vertZ< -5)
    izvtx=6;
  else if(vertZ<  0)
    izvtx=7;
  else if(vertZ<  5)
    izvtx=8;
  else if(vertZ< 10)
    izvtx=9;
  else if(vertZ< 15)
    izvtx=10;
  else if(vertZ< 20)
    izvtx=11;
  else if(vertZ< 25)
    izvtx=12;
  else if(vertZ< 30)
    izvtx=13;
  else if(vertZ< 35)
    izvtx=14;
  else
    izvtx=15;
  
  return izvtx;  
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::GetMultBin(Int_t mult){

  int imult = -1;
  
  if     (mult<2)
    imult=0;
  else if(mult<25)
    imult=mult-2;
  else
    imult=24;
  
  return imult;  
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::isGoodEsdCluster(AliESDCaloCluster* esdclust){

  int pass = 1;
  int nMinCells  = 2;
  double MinE    = 0.4;
  //double MinErat = 0;
  //double MinEcc  = 0;
  
  if (!esdclust)
    pass = 0;    
  if (!esdclust->IsEMCAL()) 
    pass = 0;
  if (esdclust->E()<MinE)
    pass = 0;
  if (esdclust->GetNCells()<nMinCells)
    pass = 0;
  //if (GetMaxCellEnergy(esdclust)/esdclust->E()<MinErat)
  //pass = 0;
  //if (esdclust->Chi2()<MinEcc) // eccentricity cut
  //pass = 0;//this is always -1.
    
  //if(esdclust->GetM02()<0.1)
  //  pass = 0;
  //if(esdclust->GetM02()>0.5)
  //  pass = 0;

  Float_t pos[3] = {0,0,0};
  esdclust->GetPosition(pos);
  TVector3 clusterPosition(pos);
  if(clusterPosition.Eta()<fEtamin || clusterPosition.Eta()>fEtamax || 
     clusterPosition.Phi()<fPhimin || clusterPosition.Phi()>fPhimax  )
    pass = 0;
  clusterPosition.Delete();
  
  //DOING THIS BY HAND NOW... 
  //if(!esdclust->GetNTracksMatched()==0)
  //pass = 0;
  
  return pass;
}

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::isGoodAodCluster(AliAODCaloCluster* aodclust){

  int pass = 1;
  int nMinCells  = 2;
  double MinE    = 0.4;
  //double MinErat = 0;
  //double MinEcc  = 0;
  
  if (!aodclust)
    pass = 0;    
  if (!aodclust->IsEMCAL()) 
    pass = 0;
  if (aodclust->E()<MinE)
    pass = 0;
  if (aodclust->GetNCells()<nMinCells)
    pass = 0;
  //if (GetMaxCellEnergy(aodclust)/aodclust->E()<MinErat)
  //pass = 0;
  //if (aodclust->Chi2()<MinEcc) // eccentricity cut
  //pass = 0;//this is always -1.
    
  //if(aodclust->GetM02()<0.1)
  //pass = 0;
  //if(aodclust->GetM02()>0.5)
  //pass = 0;

  Float_t pos[3] = {0,0,0};
  aodclust->GetPosition(pos);
  TVector3 clusterPosition(pos);
  if(clusterPosition.Eta()<fEtamin || clusterPosition.Eta()>fEtamax || 
     clusterPosition.Phi()<fPhimin || clusterPosition.Phi()>fPhimax  )
    pass = 0;
  clusterPosition.Delete();
  
  //DOING THIS BY HAND NOW... 
  //if(!aodclust->GetNTracksMatched()==0)
  //pass = 0;
  
  return pass;
}
 
//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDMpPb::getDeltaPhi(TLorentzVector p1, TLorentzVector p2){

  double dphi = p1.Phi() - p2.Phi();

  if(dphi<0.5*TMath::Pi())  
    dphi = dphi + 2.0*TMath::Pi();

  if(dphi>1.5*TMath::Pi())  
    dphi = dphi - 2.0*TMath::Pi();

  return dphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDMpPb::getDeltaEta(TLorentzVector p1, TLorentzVector p2){

  double deta = p1.PseudoRapidity() - p2.PseudoRapidity();

  return deta;
}


//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDMpPb::PrivateEnergyRecal(Double_t energy, Int_t iCalib){

  double recalibfactor = 0.0;

  if(iCalib==0){// no recalibration! 
    recalibfactor = 1.0;
  }
  else if(iCalib==1){// just a scale factor: 
    recalibfactor = 0.984;
  }
  else if(iCalib==2){// Symmetric Decay Fit - corrects data to uncorrected MC. 
    Double_t p[3] = {0.96968, -2.68720, -0.831607};
    recalibfactor = p[0] + exp(p[1] + p[2]*energy*2.0);
  }
  else if(iCalib==3){// Jason's fit to the LHC12f1a MC single photons - 04 Aug 2013 (call it kPi0MCv4??)
    Double_t p[7] = {1.00000e+00, 3.04925e-02, 4.69043e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.00046e+00};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==4){// Jason's fit to the test beam data - 04 Aug 2013(call it kBTCv3??)
    Double_t p[7] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==5){// Based on kSDM/kTBCv3 (call it kPi0MCv4??)
    Double_t p[10] = {9.78672e-01, 2.39745e-01, 6.41199e-01, 9.13538e-02, 1.46058e+02, 1.99469e+01, 9.72716e-01, 0.96968, -2.68720, -0.831607};
    recalibfactor = ( (p[6]/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) ) / ( p[7] + exp(p[8] + p[9]*energy/2.0) );
  }
  else if(iCalib==6){// kBeamTestCorrectedv2 - in AliROOT! 
    Double_t p[7] = {9.83504e-01, 2.10106e-01, 8.97274e-01, 8.29064e-02, 1.52299e+02, 3.15028e+01, 0.968};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==7){// kPi0MCv3 - in AliROOT! 
    Double_t p[7] = {9.81039e-01, 1.13508e-01, 1.00173e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.0};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==8){// Jason's fit to the noNL MC/data- based on kSDM and kPi0MCv5 - 28 Oct 2013 (call it... ??)
    Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.964, -3.132, -0.435};
    //Double_t p[10] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286, 0.96968, -2.68720, -0.831607};//same SDM piece as iCalib==2
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5]))))) * (p[7] + exp(p[8]+p[9]*energy*2.0));
  }
  else if(iCalib==9){// Jason's fit to the LHC12f1a/b MC single photons (above 400MeV), including conversions - 28 Oct 2013 (call it kPi0MCv5??)
    Double_t p[7] = {1.0, 6.64778e-02, 1.57000e+00, 9.67998e-02, 2.19381e+02, 6.31604e+01, 1.01286};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==10){// Jason played with test beam data
    Double_t p[7] = {1.0, 0.237767, 0.651203, 0.183741, 155.427, 17.0335, 0.987054};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }
  else if(iCalib==11){// Jason played with test beam MC
    Double_t p[7] = {1.0, 0.0797873, 1.68322, 0.0806098, 244.586, 116.938, 1.00437};
    recalibfactor = ((p[6])/(p[0]*(1./(1.+p[1]*exp(-energy/p[2]))*1./(1.+p[3]*exp((energy-p[4])/p[5])))));
  }

  return recalibfactor;
}


//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDMpPb::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;
  AliVCaloCells *fVCells=NULL;
  if(fEsdEv)      fVCells = fEsdEv->GetEMCALCells();
  else if(fAodEv) fVCells = fAodEv->GetEMCALCells();
  if(!fVCells)
    return 0;
  
  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = fVCells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}


//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::IsPhysPrimJ(AliMCEvent *mcEvent, Int_t iTrack){

  AliMCParticle *mcP  = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
  
  Int_t nPTracks= mcEvent->GetNumberOfPrimaries();
  
  Int_t isPhysPrimary   = 1;
  Int_t ismHF           = 0;
  Int_t ismLongLivedOrK = 0;

  if(mcP->GetMother()<0)//if it has no mother... 
    return isPhysPrimary;
  
  Int_t imTrack = mcP->GetMother();
  AliMCParticle *mcPm = static_cast<AliMCParticle*>(mcEvent->GetTrack(imTrack));
  
  if( TMath::Abs(mcPm->PdgCode())<10 )//if mother is a single quark...
    return isPhysPrimary;
  

  //############################################
  //get the PDG digits.... 
  int num = mcPm->PdgCode();
  int RevDigits[10] = {0};
  int nDigits = 0;  
  while (num >= 1){
    RevDigits[nDigits++] = num%10;
    num = num / 10;
  }
  //##############################################


  if(RevDigits[3]>3)//Baryons
    ismHF = 1;
  else if(RevDigits[2]>3)//Mesons
    ismHF = 1;
  
  ismLongLivedOrK = IsLongLivedOrK(mcPm->PdgCode());
  
  if(!ismHF && ismLongLivedOrK)
    isPhysPrimary = 0;
  else{ // check grandmother, greatgrandmothers, etc... 
    while(imTrack >= nPTracks){

      if(mcPm->GetMother()<0)//if it has no mother... 
	break;
      
      if( TMath::Abs(mcPm->PdgCode()<10) )//if mother is a single quark...
	return isPhysPrimary;
      
      imTrack = mcPm->GetMother();
      mcPm = static_cast<AliMCParticle*>(mcEvent->GetTrack(imTrack));      
      
      //############################################
      //get the PDG digits.... 
      num = mcPm->PdgCode();
      for(int i=0; i<10; i++)  RevDigits[i] = 0;
      nDigits = 0;  
      while (num >= 1){
	RevDigits[nDigits++] = num%10;
	num = num / 10;
      }
      //##############################################
      if(RevDigits[3]>3)//Baryons
	ismHF = 1;
      else if(RevDigits[2]>3)//Mesons
	ismHF = 1;
      
      ismLongLivedOrK = IsLongLivedOrK(mcPm->PdgCode());
      
      if(!ismHF && ismLongLivedOrK)
	isPhysPrimary = 0;
      
    }//while( >=nPTracks)
  }
  
  return isPhysPrimary;
}


//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDMpPb::IsLongLivedOrK(Int_t MyPDGcode){

  Int_t MyFlag = 0;

  if(
     (TMath::Abs(MyPDGcode) == 22  ) ||        // Photon
     (TMath::Abs(MyPDGcode) == 11  ) ||        // Electron
     (TMath::Abs(MyPDGcode) == 13  ) ||        // Muon(-) 
     (TMath::Abs(MyPDGcode) == 211 ) ||        // Pion
     (TMath::Abs(MyPDGcode) == 321 ) ||        // Kaon
     (TMath::Abs(MyPDGcode) == 310 ) ||        // K0s
     (TMath::Abs(MyPDGcode) == 130 ) ||        // K0l
     (TMath::Abs(MyPDGcode) == 2212) ||        // Proton 
     (TMath::Abs(MyPDGcode) == 2112) ||        // Neutron
     (TMath::Abs(MyPDGcode) == 3122) ||        // Lambda_0
     (TMath::Abs(MyPDGcode) == 3112) ||        // Sigma Minus
     (TMath::Abs(MyPDGcode) == 3222) ||        // Sigma Plus
     (TMath::Abs(MyPDGcode) == 3312) ||        // Xsi Minus 
     (TMath::Abs(MyPDGcode) == 3322) ||        // Xsi 
     (TMath::Abs(MyPDGcode) == 3334) ||        // Omega
     (TMath::Abs(MyPDGcode) == 12  ) ||        // Electron Neutrino 
     (TMath::Abs(MyPDGcode) == 14  ) ||        // Muon Neutrino
     (TMath::Abs(MyPDGcode) == 16  )   )       // Tau Neutrino
    MyFlag = 1;

  return MyFlag; 
}
