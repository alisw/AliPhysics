#include "AliAnalysisTaskEMCALMesonGGSDM.h"

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

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliEMCALGeometry.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"

#include "AliEMCALRecoUtils.h"
#include "AliExternalTrackParam.h"

// ROOT includes
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoBBox.h>
#include <TH2F.h>
#include <TArrayI.h>
#include <TArrayF.h>
#include <TObjArray.h>

// STEER includes
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
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

#include "AliGenCocktailEventHeader.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALMesonGGSDM)

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDM::AliAnalysisTaskEMCALMesonGGSDM() : 
  AliAnalysisTaskSE(),
  fOutput(0),
  fMcMode(0),
  fMyMCType(0),
  fRecalibrator(0),
  fdRmin_ClustTrack(0),
  fPhimin(0),
  fPhimax(0),
  fEtamin(0),
  fEtamax(0),
  fTrackCuts(0),
  fEsdEv(0),
  fAodEv(0),
  h1_nClusters(0), 
  h1_zvtx(0), 
  h1_trigger(0), 
  h1_M(0), 
  h1_M_mix(0), 
  h1_E(0), 
  h2_PhiEtaCluster(0), 
  h2_PhiEtaClusterCut(0), 
  h2_PhiEtaMaxCell(0), 
  h2_PhiEtaMaxCellCut(0), 
  h1_dR_ClustTrk(0),
  h2_gE_RecTruth(0), 
  h2_eop_E(0),
  h2_eop_pT(0),
  h2_E_time(0),
  h1_Pi0TruthPt(0), 
  h1_K0Pi0TruthPt(0), 
  h1_PriPi0TruthPt(0), 
  h1_PhysPi0TruthPt(0), 
  h1_Pi0TruthPtEmcal(0), 
  h1_K0Pi0TruthPtEmcal(0), 
  h1_PriPi0TruthPtEmcal(0), 
  h1_PhysPi0TruthPtEmcal(0), 
  h1_Pi0TruthPtPhi2piEta065(0), 
  h1_K0Pi0TruthPtPhi2piEta065(0), 
  h1_PriPi0TruthPtPhi2piEta065(0), 
  h1_PhysPi0TruthPtPhi2piEta065(0), 
  h1_Pi0TruthPtPhi2piEta1(0), 
  h1_K0Pi0TruthPtPhi2piEta1(0), 
  h1_PriPi0TruthPtPhi2piEta1(0), 
  h1_PhysPi0TruthPtPhi2piEta1(0), 
  h2_Pi0TruthPhiEta(0), 
  h2_PriPi0TruthPhiEta(0), 
  h2_Pi0TruthPhiEtaEmcal(0), 
  h2_PriPi0TruthPhiEtaEmcal(0), 
  h1_TruthPhotonsEmcal(0), 
  h2_TruthPhotonsPhiEta(0),
  h1_PhotonsEmcal(0), 
  h1_PhotonsNCellsCut(0), 
  h1_PhotonsTrackMatchCut(0), 
  h1_PhotonsAllCut(0), 
  h2_PhotonsPhiEtaIsEmcal(0),
  h1_dR_RealMC(0),
  h2_Mpt_Pri(0),
  h2_Mpt_Sec(0),
  h3_MptR_Sec(0),
  h3_MptR_K0s(0),
  h3_MptR_Mat(0),
  h2_PtR_MatM(0),
  h2_Mpt_Pri_conv(0),
  h2_Mpt_Sec_conv(0),
  h3_MptR_Sec_conv(0),
  h3_MptR_K0s_conv(0),
  h3_MptR_Mat_conv(0),
  h1_eConversionR(0),
  h1_PriPi0Mother(0),
  h1_SecPi0Mother(0),
  h1_Chi2(0),
  h1_nTrkMatch(0),
  h1_nCells(0),
  h1_ClusterDisp(0),
  h2_Ellipse(0),
  h2_EtaPt(0),
//h2_Mpt(0), 
  h3_MptAsymm(0), 
//h2_Mpt_mix(0), 
  h3_MptAsymm_mix(0), 
  h2_dphi_deta(0), 
  h2_dphi_deta_mix(0), 
  h2_DispRes(0),
  h2_cells_M02(0),
  TriggerList(0),
  fHelperClass(0)
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDM::AliAnalysisTaskEMCALMesonGGSDM(const char *name) :
  AliAnalysisTaskSE(name),
  fOutput(0),
  fMcMode(0),
  fMyMCType(0),
  fRecalibrator(0),
  fdRmin_ClustTrack(0),
  fPhimin(0),
  fPhimax(0),
  fEtamin(0),
  fEtamax(0),
  fTrackCuts(0),
  fEsdEv(0),
  fAodEv(0),
  h1_nClusters(0), 
  h1_zvtx(0), 
  h1_trigger(0), 
  h1_M(0), 
  h1_M_mix(0), 
  h1_E(0), 
  h2_PhiEtaCluster(0), 
  h2_PhiEtaClusterCut(0), 
  h2_PhiEtaMaxCell(0), 
  h2_PhiEtaMaxCellCut(0), 
  h1_dR_ClustTrk(0),
  h2_gE_RecTruth(0), 
  h2_eop_E(0),
  h2_eop_pT(0),
  h2_E_time(0),
  h1_Pi0TruthPt(0), 
  h1_K0Pi0TruthPt(0),
  h1_PriPi0TruthPt(0), 
  h1_PhysPi0TruthPt(0), 
  h1_Pi0TruthPtEmcal(0), 
  h1_K0Pi0TruthPtEmcal(0), 
  h1_PriPi0TruthPtEmcal(0), 
  h1_PhysPi0TruthPtEmcal(0), 
  h1_Pi0TruthPtPhi2piEta065(0), 
  h1_K0Pi0TruthPtPhi2piEta065(0), 
  h1_PriPi0TruthPtPhi2piEta065(0), 
  h1_PhysPi0TruthPtPhi2piEta065(0), 
  h1_Pi0TruthPtPhi2piEta1(0), 
  h1_K0Pi0TruthPtPhi2piEta1(0), 
  h1_PriPi0TruthPtPhi2piEta1(0), 
  h1_PhysPi0TruthPtPhi2piEta1(0), 
  h2_Pi0TruthPhiEta(0), 
  h2_PriPi0TruthPhiEta(0), 
  h2_Pi0TruthPhiEtaEmcal(0), 
  h2_PriPi0TruthPhiEtaEmcal(0), 
  h1_TruthPhotonsEmcal(0), 
  h2_TruthPhotonsPhiEta(0),
  h1_PhotonsEmcal(0), 
  h1_PhotonsNCellsCut(0), 
  h1_PhotonsTrackMatchCut(0), 
  h1_PhotonsAllCut(0), 
  h2_PhotonsPhiEtaIsEmcal(0),
  h1_dR_RealMC(0),
  h2_Mpt_Pri(0),
  h2_Mpt_Sec(0),
  h3_MptR_Sec(0),
  h3_MptR_K0s(0),
  h3_MptR_Mat(0),
  h2_PtR_MatM(0),
  h2_Mpt_Pri_conv(0),
  h2_Mpt_Sec_conv(0),
  h3_MptR_Sec_conv(0),
  h3_MptR_K0s_conv(0),
  h3_MptR_Mat_conv(0),
  h1_eConversionR(0),
  h1_PriPi0Mother(0),
  h1_SecPi0Mother(0),
  h1_Chi2(0),
  h1_nTrkMatch(0),
  h1_nCells(0),
  h1_ClusterDisp(0),
  h2_Ellipse(0),
  h2_EtaPt(0),
//h2_Mpt(0), 
  h3_MptAsymm(0), 
//h2_Mpt_mix(0), 
  h3_MptAsymm_mix(0), 
  h2_dphi_deta(0), 
  h2_dphi_deta_mix(0), 
  h2_DispRes(0), 
  h2_cells_M02(0),
  TriggerList(0),
  fHelperClass(0)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container

  DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskEMCALMesonGGSDM::~AliAnalysisTaskEMCALMesonGGSDM()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
  delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskEMCALMesonGGSDM::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
   
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);

  cout << "__________AliAnalysisTaskEMCALMesonGGSDM: Input settings__________" << endl;
  cout << " fMcMode:             " << fMcMode       << endl;
  cout << " fMyMCType:           " << fMyMCType     << endl;
  cout << " fRecalibrator:       " << fRecalibrator << endl;
  cout << " dRmin_ClustTrack:    " << fdRmin_ClustTrack << endl;
  cout << " phi range:           " << fPhimin << ", " << fPhimax << endl;
  cout << " eta range:           " << fEtamin << ", " << fEtamax << endl;
  cout << " number of zvtx bins: " << zvtx_bins     << endl;
  cout << " number of mult bins: " << mult_bins     << endl;
  cout << " poolDepth:           " << poolDepth     << endl;
  cout << endl;
  

  //AliAnalysisManager  *man = AliAnalysisManager::GetAnalysisManager();
  //AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
  //fPIDResponse = (AliPIDResponse*)inputHandler->GetPIDResponse();
  //
  //fPIDCombined = new AliPIDCombined();
  //fPIDCombined->SetSelectedSpecies(AliPID::kSPECIES);
  //fPIDCombined->SetDetectorMask(AliPIDResponse::kDetEMCAL);
  //fPIDCombined->SetEnablePriors(kFALSE);
  
  double TotalNBins = 0.0;

  // Create histograms
  Int_t nClustersbins = 501;
  Float_t nClusterslow = -0.5, nClustersup = 500.5;
  h1_nClusters = new TH1F("h1_nClusters", "# of clusters", nClustersbins, nClusterslow, nClustersup);
  h1_nClusters->GetXaxis()->SetTitle("number of clusters/evt");
  h1_nClusters->GetYaxis()->SetTitle("counts");
  h1_nClusters->SetMarkerStyle(kFullCircle);
  TotalNBins+=nClustersbins;

  Int_t nZvertexbins = 501;
  Float_t Zvertexlow = -50.0, Zvertexup = 50.0;
  h1_zvtx = new TH1F("h1_zvtx", "# of clusters", nZvertexbins, Zvertexlow, Zvertexup);
  h1_zvtx->GetXaxis()->SetTitle("z_{vertex}");
  h1_zvtx->GetYaxis()->SetTitle("counts");
  h1_zvtx->SetMarkerStyle(kFullCircle);
  TotalNBins+=nZvertexbins;

  h1_trigger = new TH1F("h1_trigger", "trigger number returned", 1001,-0.5,1000.5);
  TotalNBins+=1001;

  Int_t Mbins = 3000;
  Float_t Mlow = 0.0, Mup = 3.0;
  h1_M = new TH1F("h1_M", "Invariant Mass", Mbins, Mlow, Mup);
  h1_M->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h1_M->GetYaxis()->SetTitle("counts");
  h1_M->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins;

  h1_M_mix = new TH1F("h1_M_mix", "Invariant Mass (mixed events)", Mbins, Mlow, Mup);
  h1_M_mix->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h1_M_mix->GetYaxis()->SetTitle("counts");
  h1_M_mix->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins;

  Int_t ptbins = 2000;
  Float_t ptlow = 0.0, ptup = 20.0;
  Int_t Ebins = 1000;
  Float_t Elow = 0.0, Eup = 20.0;
  h1_E = new TH1F("h1_E", "Cluster Energy in EMCal", Ebins, Elow, Eup);
  h1_E->GetXaxis()->SetTitle("E [GeV]");
  h1_E->GetYaxis()->SetTitle("counts");
  h1_E->SetMarkerStyle(kFullCircle);
  TotalNBins+=Ebins;

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

  h1_dR_ClustTrk = new TH1F("h1_dR_ClustTrk", "Cluster-Track matching", 5000, -0.01, 5);
  h1_dR_ClustTrk->GetXaxis()->SetTitle("dR [sqrt(d#phi^{2}+d#eta^{2})]");
  h1_dR_ClustTrk->GetYaxis()->SetTitle("N");
  h1_dR_ClustTrk->SetMarkerStyle(kFullCircle);
  TotalNBins+=5000;

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

  h1_Pi0TruthPt = new TH1F("h1_Pi0TruthPt", "P_{T} distribution for Truth Pi0's", ptbins, ptlow, ptup);
  h1_Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_K0Pi0TruthPt = new TH1F("h1_K0Pi0TruthPt", "P_{T} distribution for Truth Pi0's from K^{0}_{s} decays", ptbins, ptlow, ptup);
  h1_K0Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_K0Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_K0Pi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h1_PriPi0TruthPt = new TH1F("h1_PriPi0TruthPt", "P_{T} distribution for Truth Primary Pi0's", ptbins, ptlow, ptup);
  h1_PriPi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PhysPi0TruthPt = new TH1F("h1_PhysPi0TruthPt", "P_{T} distribution for Truth Physical Primary Pi0's", ptbins, ptlow, ptup);
  h1_PhysPi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhysPi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhysPi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_Pi0TruthPtEmcal = new TH1F("h1_Pi0TruthPtEmcal", "P_{T} distribution for Truth Pi0's (hit EMCal)", ptbins, ptlow, ptup);
  h1_Pi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_K0Pi0TruthPtEmcal = new TH1F("h1_K0Pi0TruthPtEmcal", "P_{T} distribution for Truth Pi0's from K^{0}_{s} decays (hit EMCal)", ptbins, ptlow, ptup);
  h1_K0Pi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_K0Pi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_K0Pi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PriPi0TruthPtEmcal = new TH1F("h1_PriPi0TruthPtEmcal", "P_{T} distribution for Truth Primary Pi0's (hit EMCal)", ptbins, ptlow, ptup);
  h1_PriPi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PhysPi0TruthPtEmcal = new TH1F("h1_PhysPi0TruthPtEmcal", "P_{T} distribution for Truth Physical Primary Pi0's (hit EMCal)", ptbins, ptlow, ptup);
  h1_PhysPi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhysPi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhysPi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_Pi0TruthPtPhi2piEta065 = new TH1F("h1_Pi0TruthPtPhi2piEta065", 
				       "P_{T} for Truth Pi0's [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_Pi0TruthPtPhi2piEta065->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPtPhi2piEta065->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPtPhi2piEta065->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
        
  h1_K0Pi0TruthPtPhi2piEta065 = new TH1F("h1_K0Pi0TruthPtPhi2piEta065", 
					 "P_{T} for Truth Pi0's (K^{0}_{s} decays) [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_K0Pi0TruthPtPhi2piEta065->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_K0Pi0TruthPtPhi2piEta065->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_K0Pi0TruthPtPhi2piEta065->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
        
  h1_PriPi0TruthPtPhi2piEta065 = new TH1F("h1_PriPi0TruthPtPhi2piEta065",
					  "P_{T} for Primary Truth Pi0's [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_PriPi0TruthPtPhi2piEta065->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPtPhi2piEta065->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPtPhi2piEta065->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
        
  h1_PhysPi0TruthPtPhi2piEta065 = new TH1F("h1_PhysPi0TruthPtPhi2piEta065", 
					   "P_{T} for Truth Pi0's (not from material) [|#eta_{#pi^{0}}|<0.65 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_PhysPi0TruthPtPhi2piEta065->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhysPi0TruthPtPhi2piEta065->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhysPi0TruthPtPhi2piEta065->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h1_Pi0TruthPtPhi2piEta1 = new TH1F("h1_Pi0TruthPtPhi2piEta1", 
				     "P_{T} for Truth Pi0's [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_Pi0TruthPtPhi2piEta1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPtPhi2piEta1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPtPhi2piEta1->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h1_K0Pi0TruthPtPhi2piEta1 = new TH1F("h1_K0Pi0TruthPtPhi2piEta1", 
				       "P_{T} for Truth Pi0's (k^{0}_{s} decays) [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_K0Pi0TruthPtPhi2piEta1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_K0Pi0TruthPtPhi2piEta1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_K0Pi0TruthPtPhi2piEta1->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h1_PriPi0TruthPtPhi2piEta1 = new TH1F("h1_PriPi0TruthPtPhi2piEta1", 
					"P_{T} for Primary Truth Pi0's [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_PriPi0TruthPtPhi2piEta1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPtPhi2piEta1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPtPhi2piEta1->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h1_PhysPi0TruthPtPhi2piEta1 = new TH1F("h1_PhysPi0TruthPtPhi2piEta1", 
				     "P_{T} for Truth Pi0's (not from material) [|#eta_{#pi^{0}}|<1.0 && 0<#phi_{#pi^{0}}<2#pi]", ptbins, ptlow, ptup);
  h1_PhysPi0TruthPtPhi2piEta1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhysPi0TruthPtPhi2piEta1->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhysPi0TruthPtPhi2piEta1->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;
  
  h2_Pi0TruthPhiEta = new TH2F("h2_Pi0TruthPhiEta","Pi0Truth Phi vs Eta ", 380,-0.02,6.30, 200,-10,10);
  h2_Pi0TruthPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEta->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*200;

  h2_PriPi0TruthPhiEta = new TH2F("h2_PriPi0TruthPhiEta","Primary Pi0Truth Phi vs Eta ", 380,-0.02,6.30, 200,-10,10);
  h2_PriPi0TruthPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_PriPi0TruthPhiEta->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*200;

  h2_Pi0TruthPhiEtaEmcal = new TH2F("h2_Pi0TruthPhiEtaEmcal","Pi0Truth Phi vs Eta (in EMCal)", 380,-0.02,6.30, 150,-5,5);
  h2_Pi0TruthPhiEtaEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_Pi0TruthPhiEtaEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;

  h2_PriPi0TruthPhiEtaEmcal = new TH2F("h2_PriPi0TruthPhiEtaEmcal","Primary Pi0Truth Phi vs Eta (in EMCal)", 380,-0.02,6.30, 150,-5,5);
  h2_PriPi0TruthPhiEtaEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_PriPi0TruthPhiEtaEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;

  h1_TruthPhotonsEmcal = new TH1F("h1_TruthPhotonsEmcal", "P_{T} distribution for photons (in EMCal)", ptbins, ptlow, ptup);
  h1_TruthPhotonsEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_TruthPhotonsEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_TruthPhotonsEmcal->SetMarkerStyle(kFullCircle);
    TotalNBins+=ptbins;

  h2_TruthPhotonsPhiEta = new TH2F("h2_TruthPhotonsPhiEta", 
				   "Truth Photons Phi vs Eta (pointed at emcal)", 380,-0.02,6.30, 150,-1.5,1.5);
  h2_TruthPhotonsPhiEta->GetXaxis()->SetTitle("#phi [rad]");
  h2_TruthPhotonsPhiEta->GetYaxis()->SetTitle("#eta ");
  h2_TruthPhotonsPhiEta->SetMarkerStyle(kFullCircle);
  TotalNBins+=380*150;

  h1_PhotonsEmcal = new TH1F("h1_PhotonsEmcal", "P_{T} distribution for photons (in EMCal)", ptbins, ptlow, ptup);
  h1_PhotonsEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhotonsEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhotonsEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PhotonsNCellsCut = new TH1F("h1_PhotonsNCellsCut", "P_{T} distribution for #gamma's that survive NCells cut", ptbins, ptlow, ptup);
  h1_PhotonsNCellsCut->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhotonsNCellsCut->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhotonsNCellsCut->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PhotonsTrackMatchCut = new TH1F("h1_PhotonsTrackMatchCut", "P_{T} distribution for #gamma's that survive TrackMatch cut", ptbins, ptlow, ptup);
  h1_PhotonsTrackMatchCut->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhotonsTrackMatchCut->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhotonsTrackMatchCut->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PhotonsAllCut = new TH1F("h1_PhotonsAllCut", "P_{T} distribution for #gamma's that survive All cut", ptbins, ptlow, ptup);
  h1_PhotonsAllCut->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PhotonsAllCut->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PhotonsAllCut->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h2_PhotonsPhiEtaIsEmcal = new TH2F("h2_PhotonsPhiEtaIsEmcal",
				     "Photons Phi vs Eta (IsEMCAL()==1)", 380,-0.02,6.30, 100,-1.0,1.0);
  h2_PhotonsPhiEtaIsEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhotonsPhiEtaIsEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*100;
  
  h1_dR_RealMC = new TH1F("h1_dR_RealMC", "P_{T} distribution for #gamma's that survive All cut", 2000, -0.01, 10);
  h1_dR_RealMC->GetXaxis()->SetTitle("dR sqrt(dx^{2}+dy^{2})");
  h1_dR_RealMC->GetYaxis()->SetTitle("N");
  h1_dR_RealMC->SetMarkerStyle(kFullCircle);
  TotalNBins+=2000;

  h2_Mpt_Pri = new TH2F("h2_Mpt_Pri", "mass vs pT for primary pions", Mbins, Mlow, Mup, ptbins, ptlow, ptup);
  h2_Mpt_Pri->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h2_Mpt_Pri->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_Mpt_Pri->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins*ptbins;

  h2_Mpt_Sec = new TH2F("h2_Mpt_Sec", "mass vs pT for secondary pions", Mbins, Mlow, Mup, ptbins, ptlow, ptup);
  h2_Mpt_Sec->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h2_Mpt_Sec->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_Mpt_Sec->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins*ptbins;

  h3_MptR_Sec = new TH3F("h3_MptR_Sec", "mass vs pT vs production radius for secondary pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_Sec->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_Sec->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_Sec->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_Sec->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h3_MptR_K0s = new TH3F("h3_MptR_K0s", "mass vs pT vs production radius for K0s pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_K0s->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_K0s->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_K0s->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_K0s->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h3_MptR_Mat = new TH3F("h3_MptR_Mat", "mass vs pT vs production radius for material pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_Mat->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_Mat->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_Mat->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_Mat->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h2_PtR_MatM = new TH2F("h2_PtR_MatM", "pT vs production radius for merged material pions (pi mass assumed)", 100,0,20, 300,0,600);
  h2_PtR_MatM->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  h2_PtR_MatM->GetYaxis()->SetTitle("production radius [cm]");
  h2_PtR_MatM->SetMarkerStyle(kFullCircle);
  TotalNBins+=100*300;

  h2_Mpt_Pri_conv = new TH2F("h2_Mpt_Pri_conv", "mass vs pT for primary pions", Mbins, Mlow, Mup, ptbins, ptlow, ptup);
  h2_Mpt_Pri_conv->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h2_Mpt_Pri_conv->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_Mpt_Pri_conv->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins*ptbins;

  h2_Mpt_Sec_conv = new TH2F("h2_Mpt_Sec_conv", "mass vs pT for secondary pions", Mbins, Mlow, Mup, ptbins, ptlow, ptup);
  h2_Mpt_Sec_conv->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h2_Mpt_Sec_conv->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h2_Mpt_Sec_conv->SetMarkerStyle(kFullCircle);
  TotalNBins+=Mbins*ptbins;

  h3_MptR_Sec_conv = new TH3F("h3_MptR_Sec_conv", "mass vs pT vs production radius for secondary pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_Sec_conv->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_Sec_conv->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_Sec_conv->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_Sec_conv->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h3_MptR_K0s_conv = new TH3F("h3_MptR_K0s_conv", "mass vs pT vs production radius for K0s pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_K0s_conv->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_K0s_conv->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_K0s_conv->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_K0s_conv->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h3_MptR_Mat_conv = new TH3F("h3_MptR_Mat_conv", "mass vs pT vs production radius for material pions", 500,0,0.5, 100,0,20, 300,0,600);
  h3_MptR_Mat_conv->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptR_Mat_conv->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptR_Mat_conv->GetZaxis()->SetTitle("production radius [cm]");
  h3_MptR_Mat_conv->SetMarkerStyle(kFullCircle);
  TotalNBins+=500*100*300;

  h1_eConversionR = new TH1F("h1_eConversionR", "conversion point (radius)", 600,0,600);
  h1_eConversionR->GetXaxis()->SetTitle("production radius [cm]");
  h1_eConversionR->SetMarkerStyle(kFullCircle);
  TotalNBins+=600;

  h1_PriPi0Mother = new TH1F("h1_PriPi0Mother", "primary pi0 mother ID", 12001,-6000.5,6000.5);
  h1_PriPi0Mother->GetXaxis()->SetTitle("#pi^{0} mother ID");
  h1_PriPi0Mother->SetMarkerStyle(kFullCircle);
  TotalNBins+=12001;

  h1_SecPi0Mother = new TH1F("h1_SecPi0Mother", "secondray pi0 mother ID", 12001,-6000.5,6000.5);
  h1_SecPi0Mother->GetXaxis()->SetTitle("#pi^{0} mother ID");
  h1_SecPi0Mother->SetMarkerStyle(kFullCircle);
  TotalNBins+=12001;

  Int_t chi2bins = 100;
  Float_t chi2low = -2, chi2up = 2;
  h1_Chi2 = new TH1F("h1_Chi2","#chi^{2} distribution for reconstructed",chi2bins, chi2low, chi2up);
  h1_Chi2->GetXaxis()->SetTitle("#chi^{2}");
  h1_Chi2->GetYaxis()->SetTitle("counts");
  TotalNBins+=chi2bins;

  h1_nTrkMatch = new TH1F("h1_nTrkMatch","number of matched tracks",14, -1.5, 5.5);
  h1_nTrkMatch->GetXaxis()->SetTitle("nTracksMatched");
  h1_nTrkMatch->GetYaxis()->SetTitle("counts");
  TotalNBins+=14;
       
  h1_ClusterDisp = new TH1F("h1_ClusterDisp","Dispersion of CaloCluster",1000, -1, 3);
  h1_ClusterDisp->GetXaxis()->SetTitle("cluster->GetClusterDisp()");
  h1_ClusterDisp->GetYaxis()->SetTitle("counts");
  TotalNBins+=1000;
       
  h2_Ellipse = new TH2F("h2_Ellipse","Ellipse axis M20 vs M02",500, -0.01, 1, 500, -0.01, 1);
  h2_Ellipse->GetXaxis()->SetTitle("cluster->GetM20()");
  h2_Ellipse->GetYaxis()->SetTitle("cluster->GetM02()");
  h2_Ellipse->GetZaxis()->SetTitle("counts");
  TotalNBins+=500*500;

  Int_t etabins = 150;
  Float_t etalow = -1.5, etaup = 1.5;
  h2_EtaPt = new TH2F("h2_EtaPt","Cluster Energy vs ",etabins, etalow, etaup, ptbins, ptlow, ptup);
  h2_EtaPt->GetXaxis()->SetTitle("E [GeV]");
  h2_EtaPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  TotalNBins+=etabins*ptbins;

  h3_MptAsymm = new TH3F("h3_MptAsymm","mass vs p_{T} vs Asymm cut",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
  h3_MptAsymm->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptAsymm->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptAsymm->GetZaxis()->SetTitle("Asymmetry Cut (edges: 0.0, 0.1, 0.7, 1.0)");
  TotalNBins+=Mbins*ptbins*3.0;

  h3_MptAsymm_mix = new TH3F("h3_MptAsymm_mix","mass vs p_{T} vs Asymm cut (mixed events)",Mbins,Mlow,Mup, ptbins,ptlow,ptup, 3,0.5,3.5);
  h3_MptAsymm_mix->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
  h3_MptAsymm_mix->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  h3_MptAsymm_mix->GetZaxis()->SetTitle("Asymmetry Cut (edges: 0.0, 0.1, 0.7, 1.0)");
  TotalNBins+=Mbins*ptbins*3.0;

  h2_dphi_deta = new TH2F("h2_dphi_deta","#Delta#phi vs #Delta#eta", 349,-1.5,5, 400,-2.0,2.0);
  h2_dphi_deta->GetXaxis()->SetTitle("#Delta#phi");
  h2_dphi_deta->GetYaxis()->SetTitle("#Delta#eta");
  TotalNBins+=349*400;
  
  h2_dphi_deta_mix = new TH2F("h2_dphi_deta_mix","#Delta#phi vs #Delta#eta (mixed events)", 349,-1.5,5, 400,-2.0,2.0);
  h2_dphi_deta_mix->GetXaxis()->SetTitle("#Delta#phi");
  h2_dphi_deta_mix->GetYaxis()->SetTitle("#Delta#eta");
  TotalNBins+=349*400;

  h2_DispRes = new TH2F("h2_DispRes", "zvtx info", 500,-0.01,1, 500,-0.1,2);
  h2_DispRes->GetXaxis()->SetTitle("EvtVtx->GetDispersion()");
  h2_DispRes->GetYaxis()->SetTitle("EvtVtx->GetZRes()");
  h2_DispRes->GetZaxis()->SetTitle("counts");
  TotalNBins+=500*500;

  h2_cells_M02 = new TH2F("h2_cells_M02", "nCells vs M02", 204,-1.5,100.5, 500,-1,1.5);
  h2_cells_M02->GetXaxis()->SetTitle("nCells");
  h2_cells_M02->GetYaxis()->SetTitle("M02");
  h2_cells_M02->GetZaxis()->SetTitle("counts");
  TotalNBins+=204*500;

  cout << endl << "Total number of bins in booked histograms:  " << TotalNBins << endl << endl;

  // Initialize helper class (for vertex selection & pile up correction)
  fHelperClass = new AliAnalysisUtils();

  //TFile *f = OpenFile(1); 
  //TDirectory::TContext context(f);
    
  fOutput->Add(h1_nClusters);
  fOutput->Add(h1_zvtx);
  fOutput->Add(h1_trigger);
  fOutput->Add(h1_M);
  fOutput->Add(h1_M_mix);
  fOutput->Add(h1_E);
  fOutput->Add(h2_PhiEtaCluster);
  fOutput->Add(h2_PhiEtaClusterCut);
  fOutput->Add(h2_PhiEtaMaxCell);
  fOutput->Add(h2_PhiEtaMaxCellCut);
  fOutput->Add(h1_dR_ClustTrk);
  fOutput->Add(h2_gE_RecTruth);
  fOutput->Add(h2_eop_E);
  fOutput->Add(h2_eop_pT);
  fOutput->Add(h2_E_time);
  fOutput->Add(h1_Pi0TruthPt);
  fOutput->Add(h1_K0Pi0TruthPt);
  fOutput->Add(h1_PriPi0TruthPt);
  fOutput->Add(h1_PhysPi0TruthPt);
  fOutput->Add(h1_Pi0TruthPtEmcal);
  fOutput->Add(h1_K0Pi0TruthPtEmcal);
  fOutput->Add(h1_PriPi0TruthPtEmcal);
  fOutput->Add(h1_PhysPi0TruthPtEmcal);
  fOutput->Add(h1_Pi0TruthPtPhi2piEta065);
  fOutput->Add(h1_K0Pi0TruthPtPhi2piEta065);
  fOutput->Add(h1_PriPi0TruthPtPhi2piEta065);
  fOutput->Add(h1_PhysPi0TruthPtPhi2piEta065);
  fOutput->Add(h1_Pi0TruthPtPhi2piEta1);
  fOutput->Add(h1_K0Pi0TruthPtPhi2piEta1);
  fOutput->Add(h1_PriPi0TruthPtPhi2piEta1);
  fOutput->Add(h1_PhysPi0TruthPtPhi2piEta1);
  fOutput->Add(h2_Pi0TruthPhiEta);
  fOutput->Add(h2_PriPi0TruthPhiEta);
  fOutput->Add(h2_Pi0TruthPhiEtaEmcal);
  fOutput->Add(h2_PriPi0TruthPhiEtaEmcal);
  fOutput->Add(h1_TruthPhotonsEmcal);
  fOutput->Add(h2_TruthPhotonsPhiEta);
  fOutput->Add(h1_PhotonsEmcal);
  fOutput->Add(h1_PhotonsNCellsCut);
  fOutput->Add(h1_PhotonsTrackMatchCut);
  fOutput->Add(h1_PhotonsAllCut);
  fOutput->Add(h2_PhotonsPhiEtaIsEmcal);
  fOutput->Add(h1_dR_RealMC);
  fOutput->Add(h2_Mpt_Pri);
  fOutput->Add(h2_Mpt_Sec);
  fOutput->Add(h3_MptR_Sec);
  fOutput->Add(h3_MptR_K0s);
  fOutput->Add(h3_MptR_Mat);
  fOutput->Add(h2_PtR_MatM);
  fOutput->Add(h2_Mpt_Pri_conv);
  fOutput->Add(h2_Mpt_Sec_conv);
  fOutput->Add(h3_MptR_Sec_conv);
  fOutput->Add(h3_MptR_K0s_conv);
  fOutput->Add(h3_MptR_Mat_conv);
  fOutput->Add(h1_eConversionR);
  fOutput->Add(h1_PriPi0Mother);
  fOutput->Add(h1_SecPi0Mother);
  fOutput->Add(h1_Chi2);
  fOutput->Add(h1_nTrkMatch);
  fOutput->Add(h1_ClusterDisp);
  fOutput->Add(h2_Ellipse);
  fOutput->Add(h2_EtaPt);
  fOutput->Add(h3_MptAsymm);
  fOutput->Add(h3_MptAsymm_mix);
  fOutput->Add(h2_dphi_deta);
  fOutput->Add(h2_dphi_deta_mix);
  fOutput->Add(h2_DispRes);
  fOutput->Add(h2_cells_M02);

  // Post data for ALL output slots >0 here, 
  // To get at least an empty histogram 
  // 1 is the outputnumber of a certain weg of task 1  
  PostData(1, fOutput); 
}

//________________________________________________________________________
void AliAnalysisTaskEMCALMesonGGSDM::UserExec(Option_t *) 
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

  /*
  fHelperClass->SetCutOnZVertexSPD(kFALSE);//does the zvtx have to match the spd vertex? 
  fHelperClass->SetMaxVtxZ(1.0e6);//i set this myself later.. 
  // simply makes sure that there is at least 1 contributer to the zvtx determination.
  // this should only remove the *extra* events at zvtx==0.
  if(!fHelperClass->IsVertexSelected2013pA(event))
    return;
  */

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
  
  if(fEsdEv){
    TString trigClasses = fEsdEv->GetFiredTriggerClasses();
    // remove "fast cluster events": 
    if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL"))
      return;
  }
  else if(fAodEv){
    TString trigClasses = fAodEv->GetFiredTriggerClasses();
    // remove "fast cluster events": 
    if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL"))
      return;
  }
  
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

  h2_DispRes->Fill(vertDisp, vertZres);  
  // if vertex is from spd vertexZ, require more stringent cut
  if (vertIsfromZ) {
    if (vertDisp>0.02 ||  vertZres>0.25 ) 
      return; // bad vertex from VertexerZ
  }

  
  Int_t nclusters=0;
  if(fEsdEv){
    //Int_t evtN      = fEsdEv->GetEventNumberInFile();  
    //Int_t ntracks   = fEsdEv->GetNumberOfTracks();
    nclusters = fEsdEv->GetNumberOfCaloClusters();
  }
  else if(fAodEv){
    //Int_t evtN      = fAodEv->GetEventNumberInFile();  
    //Int_t ntracks   = fAodEv->GetNumberOfTracks();
    nclusters = fAodEv->GetNumberOfCaloClusters();
  }

  // EMCal cluster loop for reconstructed event
  //numberofclusters set above! 
  TLorentzVector Photon1, Photon2, Parent;
  Double_t vertex[3]; 
  Double_t E1=0.0;
  Double_t E2=0.0;
  Double_t vertZ=0.0;
  if (fEsdEv)       vertZ = fEsdEv->GetPrimaryVertex()->GetZ();
  else if (fAodEv)  vertZ = fAodEv->GetPrimaryVertex()->GetZ();    

  h1_zvtx->Fill(vertZ);
  //zvertex cut:
  if(fabs(vertZ)>10.0)
    return;
  
  h1_nClusters->Fill(nclusters);

  int izvtx = GetZvtxBin(vertZ);
  int imult = GetMultBin(nclusters);
  
  //cout << iskip << " " << izvtx << " " << imult << endl;  
  //cout << "GetNumberOfVertices(): " << fAodEv->GetNumberOfVertices() << endl;



  //######################### ~~~~~~~~~~~ ##################################
  //######################### STARTING MC ##################################
  //######################### ~~~~~~~~~~~ ##################################
  
  if(isMC){
    int isPrimary     = 0;
    int isMaterialSec    = 0;
    int isK0sDecay    = 0;

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
    Int_t MyEMCPion = -1;
    for (Int_t iTrack = 0; iTrack<nTracksMC; ++iTrack) {
      AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
      if (!mcP)
	continue;
      
      if(iTrack<nPTracksMC)  isPrimary = 1;
      else                   isPrimary = 0;
      
      isK0sDecay = 0;
      if(mcP->GetMother()>-1){
	if( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  310 ||
	    ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -310  )
	  isK0sDecay = 1;
      }      
      
      // it's a pion !! 
      if(mcP->PdgCode() != 111)
	continue;
       
      MyEMCPion = 0;
      if(strcmp(fMyMCType,"ANY")==0)
	MyEMCPion = 1;
      else if(IsMyMCHeaderType(iTrack, fMyMCType, mcEvent)){
	//cout << "evtN: " << fEsdEv->GetEventNumberInFile() << "   i: " << iTrack << "    pdg: " << mcP->PdgCode() << "   pT: " << mcP->Pt() << endl;
	//cout << "iTrack: " << iTrack << "   nPrimaryMC: " << nPTracksMC << "   nTracksMC: " << nTracksMC << endl;
	MyEMCPion = 1;
      }

      //if(MyEMCPion)
      //cout << "evtN: " << fEsdEv->GetEventNumberInFile() << "   i: " << iTrack << "    pdg: " << mcP->PdgCode() << "   pT: " << mcP->Pt() << endl;
      //cout << "evtN: " << fEsdEv->GetEventNumberInFile() << "   i: " << iTrack << "    pdg: " << static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack+1))->PdgCode() << "   pT: " << mcP->Pt() << endl;
      
      if(MyEMCPion!=1 && isPrimary==1)
	continue;
      
      
      if(isPrimary==1 && mcP->GetMother()>-1)
	h1_PriPi0Mother->Fill( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() );
      else if(isPrimary==0 && mcP->GetMother()>-1)
	h1_SecPi0Mother->Fill( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() );
      
      Int_t daughter[2] = {-1,-1};
      daughter[0] = mcP->GetDaughterFirst();
      daughter[1] = mcP->GetDaughterLast();
      
      if (daughter[0]<0)  continue;
      if (daughter[1]<0)  daughter[1]=daughter[0];      
      if (daughter[1]-daughter[0] != 1)  continue;
      
      Int_t eIndexofConvertedPhoton[2] = {-1,-1};

      bool bacc = true;
      bool binp = true;
      Double_t eta_d[2] = {0.0,0.0};
      Double_t phi_d[2] = {0.0,0.0};
      for (Int_t daughter_index=0; daughter_index<2; daughter_index++){
        const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(daughter[daughter_index]));
	eta_d[daughter_index] = dmc->Eta();
	phi_d[daughter_index] = dmc->Phi();
        if(!(dmc->PdgCode()==22))	  binp = false;
        if(!(dmc->PdgCode()==22 && 
	     eta_d[daughter_index]>fEtamin && eta_d[daughter_index]<fEtamax && 
	     phi_d[daughter_index]>fPhimin && phi_d[daughter_index]<fPhimax))   bacc = false;	

	if(dmc->GetDaughterFirst()>0 && dmc->GetDaughterLast()>0) {
	  // get the photons's daughters... 
	  const AliMCParticle *dmcd1 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetDaughterFirst()));
	  const AliMCParticle *dmcd2 = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dmc->GetDaughterLast()));
	  Double_t productionR1 = TMath::Sqrt(dmcd1->Xv()*dmcd1->Xv() + dmcd1->Yv()*dmcd1->Yv());
	  if(bacc)  h1_eConversionR->Fill(productionR1);
	  // check if this is a conversion... 
	  if( (dmcd1->PdgCode()== -1.0*dmcd2->PdgCode()) &&
	      (dmcd1->PdgCode()==11 || dmcd1->PdgCode()==-11) &&
	      productionR1<440.0){
	    //find the conv e with highest energy, assign it to be that photon decay product.
	    if( dmcd1->E() > dmcd2->E() )
	      eIndexofConvertedPhoton[daughter_index] = dmc->GetDaughterFirst();
	    else
	      eIndexofConvertedPhoton[daughter_index] = dmc->GetDaughterLast();
	  }
	}
      }

      if(binp!=true)
	continue;

      // primary particle
      //Double_t dR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
      //                          (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
      //if(dR <= 0.01)  isPrimary = 1;
      //else            isPrimary = 0;

      isMaterialSec = 0;
      if(isPrimary!=1){
	if(((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2212 || //proton
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2212 || //anti-proton
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  2112 || //neutron
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -2112 || //anti-neutron
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  321  || //K+
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -321  || //K-
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  211  || //pi+
	   ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -211     //pi-
	   )
	  isMaterialSec = 1;
      }

      h1_Pi0TruthPt                  ->Fill(mcP->Pt());
      if(isK0sDecay)  h1_K0Pi0TruthPt->Fill(mcP->Pt());
      h2_Pi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());
      
      if(isPrimary==1){
	h1_PriPi0TruthPt    ->Fill(mcP->Pt());
	h2_PriPi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());
      }
      if(isPrimary!=1 && isMaterialSec!=1)   h1_PhysPi0TruthPt->Fill(mcP->Pt());
     
      if(mcP->Eta()<-1.0 || mcP->Eta()>1.0)
	continue;
      
      h1_Pi0TruthPtPhi2piEta1         ->Fill(mcP->Pt());
      if(isPrimary==1)
	h1_PriPi0TruthPtPhi2piEta1    ->Fill(mcP->Pt());
      if(isK0sDecay)      
	h1_K0Pi0TruthPtPhi2piEta1     ->Fill(mcP->Pt());
      if(isPrimary!=1 && isMaterialSec!=1)
	h1_PhysPi0TruthPtPhi2piEta1   ->Fill(mcP->Pt());
      
      if(mcP->Eta()>fEtamin && mcP->Eta()<fEtamax){
	h1_Pi0TruthPtPhi2piEta065       ->Fill(mcP->Pt());
	if(isPrimary==1)
	  h1_PriPi0TruthPtPhi2piEta065  ->Fill(mcP->Pt());
	if(isK0sDecay)      
	  h1_K0Pi0TruthPtPhi2piEta065   ->Fill(mcP->Pt());
	if(isPrimary!=1 && isMaterialSec!=1)
	  h1_PhysPi0TruthPtPhi2piEta065 ->Fill(mcP->Pt());
      }      
      
      
      if(binp && bacc){// 2 Photons hit the EMCAL! 
	
	Int_t Nfoundphotons = 0;
	Int_t iFoundphotons[10] = {0,0,0,0,0,0,0,0,0,0};
	Int_t Nfoundelectrons = 0;
	Int_t iFoundelectrons[10] = {0,0,0,0,0,0,0,0,0,0};
	for (Int_t daughter_index=0; daughter_index<2; daughter_index++){//both truth photons. (also includes conversions..)
	  
	  const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(daughter[daughter_index]));
	  
	  h1_TruthPhotonsEmcal->Fill(dmc->Pt());
	  h2_TruthPhotonsPhiEta->Fill(dmc->Phi(),dmc->Eta());
	  
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
		
		Double_t dR = TMath::Sqrt((eta_d[daughter_index]-clustMC_eta)*(eta_d[daughter_index]-clustMC_eta) + 
					  (phi_d[daughter_index]-clustMC_phi)*(phi_d[daughter_index]-clustMC_phi));
		h1_dR_RealMC->Fill(dR);
		matches_pion_photon = 0;
		//if(dR<=0.04) matches_pion_photon = 1;
		

		TArrayI *TruthLabelsA = esdCluster->GetLabelsArray();
		if(TruthLabelsA){
		  Int_t trackindex = TruthLabelsA->At(0);
		  if( trackindex==daughter[daughter_index] ){
		    matches_pion_photon = 1;
		    iFoundphotons[Nfoundphotons] = i;
		    Nfoundphotons++;
		  }
		  else if( trackindex==eIndexofConvertedPhoton[daughter_index] ){
		    iFoundelectrons[Nfoundelectrons] = i;
		    Nfoundelectrons++;
		  }
		  AliMCParticle *truthP = (AliMCParticle*)(mcEvent->GetTrack(trackindex));
		  
		  if(matches_pion_photon){
		    
		    h1_PhotonsEmcal->Fill(esdCluster->E());
		    h2_PhotonsPhiEtaIsEmcal->Fill(clustMC_phi,clustMC_eta);
		    if(esdCluster->GetNCells()>=2)
		      h1_PhotonsNCellsCut->Fill(esdCluster->E());
		    if(esdCluster->GetNTracksMatched()==0)
		      h1_PhotonsTrackMatchCut->Fill(esdCluster->E());
		    if(esdCluster->GetNCells()>=2 && esdCluster->GetNTracksMatched()==0)
		      h1_PhotonsAllCut->Fill(esdCluster->E());		  
		    
		    if(esdCluster->GetNCells()>=2){
		      recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);			    
		      h2_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()));
		    }
		  }//if(matches_pion_photon)
		  
		}//if Truthlabels exists
		vpos.Delete();
	      }//if(IsEMCAL())
	      
	    }//if(fEsdEv)
	    else if(fAodEv){
	      
	      if(aodCluster->IsEMCAL()){
		
		Float_t pos[3] = {0,0,0};
		aodCluster->GetPosition(pos);  
		TVector3 vpos(pos); 
		//h1_Phi->Fill(vpos.Phi());
		clustMC_phi = vpos.Phi();
		clustMC_eta = vpos.Eta();
		
		Double_t dR = TMath::Sqrt((eta_d[daughter_index]-clustMC_eta)*(eta_d[daughter_index]-clustMC_eta) + 
					  (phi_d[daughter_index]-clustMC_phi)*(phi_d[daughter_index]-clustMC_phi));
		h1_dR_RealMC->Fill(dR);
		matches_pion_photon = 0;
		if(dR<=0.04) matches_pion_photon = 1;
		
		//TArrayI *TruthLabelsA = aodCluster->GetLabelsArray();
		//if(TruthLabelsA){
		//  Int_t trackindex = TruthLabelsA->At(0);
		//  if( trackindex==daughter[daughter_index] )
		//    matches_pion_photon = 1;
		//  AliMCParticle *truthP = (AliMCParticle*)(mcEvent->GetTrack(trackindex));
		  		  		  
		if(matches_pion_photon){		
		  
		  h1_PhotonsEmcal->Fill(aodCluster->E());
		  h2_PhotonsPhiEtaIsEmcal->Fill(clustMC_phi,clustMC_eta);
		  if(aodCluster->GetNCells()>=2)
		    h1_PhotonsNCellsCut->Fill(aodCluster->E());
		  if(aodCluster->GetNTracksMatched()==0)
		    h1_PhotonsTrackMatchCut->Fill(aodCluster->E());
		  if(aodCluster->GetNCells()>=2 && aodCluster->GetNTracksMatched()==0)
		    h1_PhotonsAllCut->Fill(aodCluster->E());		  
		  
		  //if(aodCluster->GetNCells()>=2){
		  //recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);			    
		  //h2_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()));
		  //}
		}//if(matches_pion_photon)
		
		//}//if Truthlabels exists
		vpos.Delete();
	      }//if(IsEMCAL())	      
	      
	    }//if(fAodEv)
	    
	  }//loop over nclusters. 
	  
	}//both truth photons.

	if(Nfoundphotons>1){
	  AliESDCaloCluster* esdCluster1 = fEsdEv->GetCaloCluster(iFoundphotons[0]); // pointer to EMCal cluster
	  AliESDCaloCluster* esdCluster2 = fEsdEv->GetCaloCluster(iFoundphotons[1]); // pointer to EMCal cluster

	  if( isGoodEsdCluster(esdCluster1) && isGoodEsdCluster(esdCluster2) ){

	    recalScale = PrivateEnergyRecal(esdCluster1->E(), fRecalibrator);
	    E1 = esdCluster1->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster1->GetMomentum(Photon1,vertex);
	    Photon1.SetPx(Photon1.Px()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPy(Photon1.Py()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPz(Photon1.Pz()*recalScale);// TOTAL HACK - JJ

	    recalScale = PrivateEnergyRecal(esdCluster2->E(), fRecalibrator);
	    E2 = esdCluster2->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster2->GetMomentum(Photon2,vertex);
	    Photon2.SetPx(Photon2.Px()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPy(Photon2.Py()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPz(Photon2.Pz()*recalScale);// TOTAL HACK - JJ

	    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
	    	  
	    //double productionR = TMath::Sqrt( mcP->Xv()*mcP->Xv() + mcP->Yv()*mcP->Yv() );
	    Double_t productionR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
					       (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
	    if(isPrimary==1){
	      //cout << "Primary production vertex: " << productionR << endl;
	      h2_Mpt_Pri->Fill(Parent.M(),Parent.Pt());
	    }
	    else{
	      //cout << "Secondary production vertex: " << productionR << endl;
	      h2_Mpt_Sec ->Fill(Parent.M(),Parent.Pt());	    
	      h3_MptR_Sec->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isK0sDecay)
		h3_MptR_K0s->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isMaterialSec)
		h3_MptR_Mat->Fill(Parent.M(),Parent.Pt(),productionR);
	    }
	  }//both good clusters
	}//found 2 photons.
	else if(Nfoundphotons==1){
	  int mergedPion = 0;
	  AliESDCaloCluster* esdCluster1 = fEsdEv->GetCaloCluster(iFoundphotons[0]); // pointer to EMCal cluster
	  
	  TArrayI *TruthLabelsA = esdCluster1->GetLabelsArray();
	  if(TruthLabelsA){
	    if(TruthLabelsA->GetSize()>1){
	      Int_t trackindex[2];
	      trackindex[0] = TruthLabelsA->At(0);
	      trackindex[1] = TruthLabelsA->At(1);
	      if( (trackindex[0]==daughter[0] && trackindex[1]==daughter[1]) ||
		  (trackindex[0]==daughter[1] && trackindex[1]==daughter[0]) ){
		mergedPion = 1;
	      }
	      if(mergedPion==1){
		recalScale = PrivateEnergyRecal(esdCluster1->E(), fRecalibrator);
		E1 = esdCluster1->E()*recalScale;// TOTAL HACK - JJ
		Double_t productionR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
						   (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
		h2_PtR_MatM->Fill(E1,productionR);
	      }//if merged pion.
	    }//truthlabel.size > 1
	  }//if truthlabels
	}// Nfoundphotons==1

	if(Nfoundphotons==1 && Nfoundelectrons==1){
	  AliESDCaloCluster* esdCluster1 = fEsdEv->GetCaloCluster(iFoundphotons[0]); // pointer to EMCal cluster
	  AliESDCaloCluster* esdCluster2 = fEsdEv->GetCaloCluster(iFoundelectrons[0]); // pointer to EMCal cluster

	  if( isGoodEsdCluster(esdCluster1) && isGoodEsdCluster(esdCluster2) ){

	    recalScale = PrivateEnergyRecal(esdCluster1->E(), fRecalibrator);
	    E1 = esdCluster1->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster1->GetMomentum(Photon1,vertex);
	    Photon1.SetPx(Photon1.Px()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPy(Photon1.Py()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPz(Photon1.Pz()*recalScale);// TOTAL HACK - JJ

	    recalScale = PrivateEnergyRecal(esdCluster2->E(), fRecalibrator);
	    E2 = esdCluster2->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster2->GetMomentum(Photon2,vertex);
	    Photon2.SetPx(Photon2.Px()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPy(Photon2.Py()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPz(Photon2.Pz()*recalScale);// TOTAL HACK - JJ

	    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
	    	  
	    //double productionR = TMath::Sqrt( mcP->Xv()*mcP->Xv() + mcP->Yv()*mcP->Yv() );
	    Double_t productionR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
					       (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
	    if(isPrimary==1){
	      //cout << "Primary production vertex: " << productionR << endl;
	      h2_Mpt_Pri_conv->Fill(Parent.M(),Parent.Pt());
	    }
	    else{
	      //cout << "Secondary production vertex: " << productionR << endl;
	      h2_Mpt_Sec_conv ->Fill(Parent.M(),Parent.Pt());	    
	      h3_MptR_Sec_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isK0sDecay)
		h3_MptR_K0s_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isMaterialSec)
		h3_MptR_Mat_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	    }
	  }//both good clusters
	}// Nfoundphotons==1 && Nfoundelectrons==1
	else if(Nfoundelectrons==2){
	  AliESDCaloCluster* esdCluster1 = fEsdEv->GetCaloCluster(iFoundelectrons[0]); // pointer to EMCal cluster
	  AliESDCaloCluster* esdCluster2 = fEsdEv->GetCaloCluster(iFoundelectrons[1]); // pointer to EMCal cluster

	  if( isGoodEsdCluster(esdCluster1) && isGoodEsdCluster(esdCluster2) ){

	    recalScale = PrivateEnergyRecal(esdCluster1->E(), fRecalibrator);
	    E1 = esdCluster1->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster1->GetMomentum(Photon1,vertex);
	    Photon1.SetPx(Photon1.Px()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPy(Photon1.Py()*recalScale);// TOTAL HACK - JJ
	    Photon1.SetPz(Photon1.Pz()*recalScale);// TOTAL HACK - JJ

	    recalScale = PrivateEnergyRecal(esdCluster2->E(), fRecalibrator);
	    E2 = esdCluster2->E()*recalScale;// TOTAL HACK - JJ
	    fEsdEv->GetVertex()->GetXYZ(vertex);
	    esdCluster2->GetMomentum(Photon2,vertex);
	    Photon2.SetPx(Photon2.Px()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPy(Photon2.Py()*recalScale);// TOTAL HACK - JJ
	    Photon2.SetPz(Photon2.Pz()*recalScale);// TOTAL HACK - JJ

	    Parent =  TLorentzVector(Photon1.Px(),Photon1.Py(),Photon1.Pz(),E1) + TLorentzVector(Photon2.Px(),Photon2.Py(),Photon2.Pz(),E2);
	    	  
	    //double productionR = TMath::Sqrt( mcP->Xv()*mcP->Xv() + mcP->Yv()*mcP->Yv() );
	    Double_t productionR = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
					       (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
	    if(isPrimary==1){
	      //cout << "Primary production vertex: " << productionR << endl;
	      h2_Mpt_Pri_conv->Fill(Parent.M(),Parent.Pt());
	    }
	    else{
	      //cout << "Secondary production vertex: " << productionR << endl;
	      h2_Mpt_Sec_conv ->Fill(Parent.M(),Parent.Pt());	    
	      h3_MptR_Sec_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isK0sDecay)
		h3_MptR_K0s_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	      if(isMaterialSec)
		h3_MptR_Mat_conv->Fill(Parent.M(),Parent.Pt(),productionR);
	    }
	  }//both good clusters
	}// Nfoundelectrons==2
	
	h1_Pi0TruthPtEmcal    ->Fill(mcP->Pt());
	if(isK0sDecay)
	  h1_K0Pi0TruthPtEmcal    ->Fill(mcP->Pt());
	h2_Pi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());	
	
	if(isPrimary==1){
	  h1_PriPi0TruthPtEmcal    ->Fill(mcP->Pt());
	  h2_PriPi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());
	}
	if(isPrimary!=1 && isMaterialSec!=1)   h1_PhysPi0TruthPtEmcal->Fill(mcP->Pt());
	
	
      }// 2 Photons hit the EMCAL! 
      
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
	  
	  h1_dR_ClustTrk->Fill(dR_clusttrk);
	  
	  //uncomment this to do the track matching (2 of 3 lines, esd part)!! 
	  //if(dR_clusttrk<fdRmin_ClustTrack)
	  //MatchesToTrack = 1;
	  
	}//_____________________________nTracks__________________________
	
	h2_cells_M02  ->Fill(esdCluster->GetNCells(),esdCluster->GetM02());
	h2_Ellipse    ->Fill(esdCluster->GetM20(),esdCluster->GetM02());
	h1_Chi2       ->Fill(esdCluster->Chi2());//always -1. 
	h1_nTrkMatch  ->Fill(esdCluster->GetNTracksMatched());
	h1_ClusterDisp->Fill(esdCluster->GetDispersion());
	h2_E_time     ->Fill(esdCluster->E(),esdCluster->GetTOF());
	
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

	//uncomment this to do the track matching (2 of 3 lines, esd part)!! 
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
	  h1_E->Fill(E1);
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
	  Double_t momTrk[3] = {0,0,0};
	  aodTrack->GetXYZ(posTrk);
	  aodTrack->GetPxPyPz(momTrk);
	  //TVector3 vposTrk(posTrk);
	  
	  //####################################################################################################	  
	  //
	  // commented all this stuff just to satisfy aliroot warnings. 
	  // but I may need it again if I want to do the track matching for aods. 
	  /*
	  Double_t fMass          = 0.139;
	  Double_t fStepSurface   = 20.;
	  Float_t etaproj=0.0;
	  Float_t phiproj=0.0;
	  Float_t pttrackproj=0.0;

	  Double_t cv[21] = {0.0};	  
	  aodTrack->GetCovarianceXYZPxPyPz(cv);
	  AliExternalTrackParam *trackParam = new AliExternalTrackParam(posTrk,momTrk,cv,aodTrack->Charge());	  
	  //AliExternalTrackParam emcalParam(*trackParam);
	  //AliExternalTrackParam *trackParam =  const_cast<AliExternalTrackParam*>(aodTrack->GetInnerParam());
	  if(!trackParam) continue;
	  ////AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(trackParam, 440., fMass, fStepSurface, etaproj, phiproj, pttrackproj);
	  //AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcalParam, 440., fMass, fStepSurface, etaproj, phiproj, pttrackproj);
	  delete trackParam;

	  //Constantin's implementation... gives funny result. 
	  //AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(aodTrack,440.0);
	  //phiproj = aodTrack->GetTrackPhiOnEMCal();
	  //etaproj = aodTrack->GetTrackPhiOnEMCal();
	  
	  double dR_clusttrk = sqrt((phiproj-clusterPosition.Phi())*(phiproj-clusterPosition.Phi()) + 
				    (etaproj-clusterPosition.Eta())*(etaproj-clusterPosition.Eta()) );

	  h1_dR_ClustTrk->Fill(dR_clusttrk);
	  */
	  //####################################################################################################

	  //uncomment this to do the track matching (2 of 3 lines, aod part)!! 
	  //if(dR_clusttrk<fdRmin_ClustTrack)
	  //MatchesToTrack = 1;
	  	  

	}//_____________________________nTracks__________________________

	h2_cells_M02  ->Fill(aodCluster->GetNCells(),aodCluster->GetM02());
	h2_Ellipse    ->Fill(aodCluster->GetM20(),aodCluster->GetM02());
	h1_Chi2       ->Fill(aodCluster->Chi2());//always -1. 
	h1_nTrkMatch  ->Fill(aodCluster->GetNTracksMatched());
	h1_ClusterDisp->Fill(aodCluster->GetDispersion());
	h2_E_time     ->Fill(aodCluster->E(),aodCluster->GetTOF());

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
	  h1_E->Fill(E1);
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
      	  
      h1_M        ->Fill(Parent.M());
      h3_MptAsymm ->Fill(Parent.M(),Parent.Pt(),asymCut);
      h2_dphi_deta->Fill(deltaphi,deltaeta);
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

	h1_M_mix        ->Fill(Parent.M());
	h3_MptAsymm_mix ->Fill(Parent.M(),Parent.Pt(),asymCut);
	h2_dphi_deta_mix->Fill(deltaphi,deltaeta);
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
void AliAnalysisTaskEMCALMesonGGSDM::Terminate(Option_t *) //specify what you want to have done
{
  // Called once at the end of the query.
  
}


// //________________________________________________________________________
// Int_t AliAnalysisTaskEMCALMesonGGSDM::GetZvtxBin(Double_t vertZ)
// {
//   
//   int izvtx = -1;
//   
//   if     (vertZ<-35)
//     izvtx=0;
//   else if(vertZ<-30)
//     izvtx=1;
//   else if(vertZ<-25)
//     izvtx=2;
//   else if(vertZ<-20)
//     izvtx=3;
//   else if(vertZ<-15)
//     izvtx=4;
//   else if(vertZ<-10)
//     izvtx=5;
//   else if(vertZ< -5)
//     izvtx=6;
//   else if(vertZ<  0)
//     izvtx=7;
//   else if(vertZ<  5)
//     izvtx=8;
//   else if(vertZ< 10)
//     izvtx=9;
//   else if(vertZ< 15)
//     izvtx=10;
//   else if(vertZ< 20)
//     izvtx=11;
//   else if(vertZ< 25)
//     izvtx=12;
//   else if(vertZ< 30)
//     izvtx=13;
//   else if(vertZ< 35)
//     izvtx=14;
//   else
//     izvtx=15;
//   
//   return izvtx;  
// }

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDM::GetZvtxBin(Double_t vertZ)
{
  
  int izvtx = -1;
  
  if     (vertZ<-3.375)
    izvtx=0;
  else if(vertZ<-1.605)
    izvtx=1;
  else if(vertZ<-0.225)
    izvtx=2;
  else if(vertZ<1.065)
    izvtx=3;
  else if(vertZ<-2.445)
    izvtx=4;
  else if(vertZ<-4.245)
    izvtx=5;
  else
    izvtx=6;
  
  return izvtx;  
}


// //________________________________________________________________________
// Int_t AliAnalysisTaskEMCALMesonGGSDM::GetMultBin(Int_t mult){
// 
//   int imult = -1;
//   
//   if     (mult<2)
//     imult=0;
//   else if(mult<25)
//     imult=mult-2;
//   else
//     imult=24;
//   
//   return imult;  
// }

//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDM::GetMultBin(Int_t mult){

  int imult = -1;
  
  if     (mult<2)
    imult=0;
  else if(mult<3)
    imult=1;
  else if(mult<4)
	imult=2;  
  else if(mult<8)	  
    imult=3;
  else if(mult<15)	  
    imult=4;
  else
    imult=5;
  
  return imult;  
}


//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDM::isGoodEsdCluster(AliESDCaloCluster* esdclust){

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
Int_t AliAnalysisTaskEMCALMesonGGSDM::isGoodAodCluster(AliAODCaloCluster* aodclust){

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
Double_t AliAnalysisTaskEMCALMesonGGSDM::getDeltaPhi(TLorentzVector p1, TLorentzVector p2){

  double dphi = p1.Phi() - p2.Phi();

  if(dphi<0.5*TMath::Pi())  
    dphi = dphi + 2.0*TMath::Pi();

  if(dphi>1.5*TMath::Pi())  
    dphi = dphi - 2.0*TMath::Pi();

  return dphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDM::getDeltaEta(TLorentzVector p1, TLorentzVector p2){

  double deta = p1.PseudoRapidity() - p2.PseudoRapidity();

  return deta;
}


//________________________________________________________________________
Double_t AliAnalysisTaskEMCALMesonGGSDM::PrivateEnergyRecal(Double_t energy, Int_t iCalib){

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
Double_t AliAnalysisTaskEMCALMesonGGSDM::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
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
Int_t AliAnalysisTaskEMCALMesonGGSDM::IsPhysPrimJ(AliMCEvent *mcEvent, Int_t iTrack){

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
Int_t AliAnalysisTaskEMCALMesonGGSDM::IsLongLivedOrK(Int_t MyPDGcode){

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


//________________________________________________________________________
Int_t AliAnalysisTaskEMCALMesonGGSDM::IsMyMCHeaderType(Int_t iTrack, char *MyType, AliMCEvent *mcEvent) const
{

  Int_t isMyType = 0;

  AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(mcEvent->GenEventHeader());
  if(!cocktail)
    return 0;

  TList *genHeaders = cocktail->GetHeaders();
  
  Int_t nGenerators = genHeaders->GetEntries();
  Int_t indexMyType = -1;
  Int_t startParticle=0;

  for(Int_t igen = 0; igen < nGenerators; igen++){
    AliGenEventHeader* eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
    TString name = eventHeader2->GetName();
    startParticle += eventHeader2->NProduced();
    //cout << name << endl;
    if (name.Contains(MyType,TString::kIgnoreCase)){
      indexMyType = igen;
      startParticle -= eventHeader2->NProduced();
      break;
    }
  }

  AliGenEventHeader *addedPi0Header = (AliGenEventHeader*)genHeaders->At(indexMyType);
  Int_t ipi0min = startParticle;
  Int_t ipi0max = ipi0min+addedPi0Header->NProduced()-1;
  if(iTrack >= ipi0min && iTrack <= ipi0max)
    isMyType = 1;
  
  return isMyType; 
}


