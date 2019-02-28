#include "AliAnalysisTaskSDMGammaMC.h"
 
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

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSDMGammaMC)

//________________________________________________________________________
AliAnalysisTaskSDMGammaMC::AliAnalysisTaskSDMGammaMC() : 
  AliAnalysisTaskSE(),
  fOutput(0),
  fMcMode(0),
  fRecalibrator(0),
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
  h1_E(0), 
  h1_Phi(0), 
  h2_PiMotherID(0), 
  h2_GaMotherID(0), 
  h3_gE_RecTruth(0), 
  h3_gE_RecTruth_ncellscut(0), 
  h1_Pi0TruthPt(0), 
  h1_PriPi0TruthPt(0), 
  h1_Pi0TruthPtEmcal(0), 
  h1_PriPi0TruthPtEmcal(0), 
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
  h1_Eta(0),
  h1_Chi2(0),
  h1_nTrkMatch(0),
  h1_nCells(0),
  h1_ClusterDisp(0),
  h2_Ellipse(0),
  h2_EtaPt(0),
  h2_dphi_deta(0), 
  h2_dphi_deta_mix(0), 
  h2_DispRes(0),
  h2_cells_M02(0),
  TriggerList(0)
{
  // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskSDMGammaMC::AliAnalysisTaskSDMGammaMC(const char *name) :
  AliAnalysisTaskSE(name),
  fOutput(0),
  fMcMode(0),
  fRecalibrator(0),
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
  h1_E(0), 
  h1_Phi(0), 
  h2_PiMotherID(0), 
  h2_GaMotherID(0), 
  h3_gE_RecTruth(0), 
  h3_gE_RecTruth_ncellscut(0), 
  h1_Pi0TruthPt(0), 
  h1_PriPi0TruthPt(0), 
  h1_Pi0TruthPtEmcal(0), 
  h1_PriPi0TruthPtEmcal(0), 
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
  h1_Eta(0),
  h1_Chi2(0),
  h1_nTrkMatch(0),
  h1_nCells(0),
  h1_ClusterDisp(0),
  h2_Ellipse(0),
  h2_EtaPt(0),
  h2_dphi_deta(0), 
  h2_dphi_deta_mix(0), 
  h2_DispRes(0), 
  h2_cells_M02(0),
  TriggerList(0)
{
  // Constructor
  // Define input and output slots here (never in the dummy constructor)
  // Input slot #0 works with a TChain - it is connected to the default input container
  // Output slot #1 writes into a TH1 container


  DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskSDMGammaMC::~AliAnalysisTaskSDMGammaMC()
{
  // Destructor. Clean-up the output list, but not the histograms that are put inside
  // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutput;
  }
  delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskSDMGammaMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once (on the worker node)

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
   
  fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);

  cout << "__________AliAnalysisTaskSDMGammaMC: Input settings__________" << endl;
  cout << " fMcMode:             " << fMcMode   << endl;
  cout << " fRecalibrator:       " << fRecalibrator << endl;
  cout << " phi range:           " << fPhimin << ", " << fPhimax << endl;
  cout << " eta range:           " << fEtamin << ", " << fEtamax << endl;
  cout << " number of zvtx bins: " << zvtx_bins << endl;
  cout << " number of mult bins: " << mult_bins << endl;
  cout << " poolDepth:           " << poolDepth << endl;
  cout << endl;
  

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

  Int_t ptbins = 2000;
  Float_t ptlow = 0.0, ptup = 20.0;
  Int_t Ebins = 1000;
  Float_t Elow = 0.0, Eup = 20.0;
  h1_E = new TH1F("h1_E", "Cluster Energy in EMCal", Ebins, Elow, Eup);
  h1_E->GetXaxis()->SetTitle("E [GeV]");
  h1_E->GetYaxis()->SetTitle("counts");
  h1_E->SetMarkerStyle(kFullCircle);
  TotalNBins+=Ebins;

  h1_Phi = new TH1F("h1_Phi", "phi distribution", 1000, -7, 7);
  h1_Phi->GetXaxis()->SetTitle("#phi [rad]");
  h1_Phi->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Phi->SetMarkerStyle(kFullCircle);
  TotalNBins+=1000;

  h2_PiMotherID = new TH2F("h2_PiMotherID", "Mother ID for Truth Pi0's", 100001, -0.5,100000.5, 2, 0.5,2.5);
  h2_PiMotherID->GetXaxis()->SetTitle("#pi^{0} Mother Particle ID");
  h2_PiMotherID->GetYaxis()->SetTitle("primary or non-primary");
  h2_PiMotherID->GetZaxis()->SetTitle("counts");
  h2_PiMotherID->SetMarkerStyle(kFullCircle);
  TotalNBins+=2*100001;

  h2_GaMotherID = new TH2F("h2_GaMotherID", "Mother ID for Truth #gamma's", 100001, -0.5,100000.5, 2, 0.5,2.5);
  h2_GaMotherID->GetXaxis()->SetTitle("#gamma Mother Particle ID");
  h2_GaMotherID->GetYaxis()->SetTitle("primary or non-primary");
  h2_GaMotherID->GetZaxis()->SetTitle("counts");
  h2_GaMotherID->SetMarkerStyle(kFullCircle);
  TotalNBins+=2*100001;

  h3_gE_RecTruth = new TH3F("h3_gE_RecTruth", "#gamma E_{truth}/E_{clust} vs E_{clust}", Ebins,Elow,Eup, 1000,0,4, 4,0.5,4.5);
  h3_gE_RecTruth->GetXaxis()->SetTitle("E^{rec}_{clust} [GeV]");
  h3_gE_RecTruth->GetYaxis()->SetTitle("E^{rec}_{clust}/E^{truth}_{#gamma}");
  h3_gE_RecTruth->GetZaxis()->SetTitle("category");
  h3_gE_RecTruth->SetMarkerStyle(kFullCircle);
  TotalNBins+=Ebins*1000*4;
  
  h3_gE_RecTruth_ncellscut = new TH3F("h3_gE_RecTruth_ncellscut", "#gamma E_{truth}/E_{clust} vs E_{clust} (for nCells>1)", Ebins,Elow,Eup, 1000,0,4, 4,0.5,4.5);
  h3_gE_RecTruth_ncellscut->GetXaxis()->SetTitle("E^{rec}_{clust} [GeV]");
  h3_gE_RecTruth_ncellscut->GetYaxis()->SetTitle("E^{rec}_{clust}/E^{truth}_{#gamma}");
  h3_gE_RecTruth_ncellscut->GetZaxis()->SetTitle("category");
  h3_gE_RecTruth_ncellscut->SetMarkerStyle(kFullCircle);
  TotalNBins+=Ebins*1000*4;
  
  h1_Pi0TruthPt = new TH1F("h1_Pi0TruthPt", "P_{T} distribution for Truth Pi0's", ptbins, ptlow, ptup);
  h1_Pi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PriPi0TruthPt = new TH1F("h1_PriPi0TruthPt", "P_{T} distribution for Truth Primary Pi0's", ptbins, ptlow, ptup);
  h1_PriPi0TruthPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPt->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_Pi0TruthPtEmcal = new TH1F("h1_Pi0TruthPtEmcal", "P_{T} distribution for Truth Pi0's (hit EMCal)", ptbins, ptlow, ptup);
  h1_Pi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_Pi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_Pi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
  TotalNBins+=ptbins;

  h1_PriPi0TruthPtEmcal = new TH1F("h1_PriPi0TruthPtEmcal", "P_{T} distribution for Truth Primary Pi0's (hit EMCal)", ptbins, ptlow, ptup);
  h1_PriPi0TruthPtEmcal->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  h1_PriPi0TruthPtEmcal->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  h1_PriPi0TruthPtEmcal->SetMarkerStyle(kFullCircle);
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
				     "Photons Phi vs Eta (IsEMCAL()==1)", 380,-0.02,6.30, 150,-1.5,1.5);
  h2_PhotonsPhiEtaIsEmcal->GetXaxis()->SetTitle("#phi [rad]");
  h2_PhotonsPhiEtaIsEmcal->GetYaxis()->SetTitle("#eta ");
  TotalNBins+=380*150;
  
  h1_dR_RealMC = new TH1F("h1_dR_RealMC", "P_{T} distribution for #gamma's that survive All cut", 2000, -0.01, 10);
  h1_dR_RealMC->GetXaxis()->SetTitle("dR sqrt(dx^{2}+dy^{2})");
  h1_dR_RealMC->GetYaxis()->SetTitle("N");
  h1_dR_RealMC->SetMarkerStyle(kFullCircle);
  TotalNBins+=2000;

  Int_t etabins = 150;
  Float_t etalow = -1.5, etaup = 1.5;
  h1_Eta = new TH1F("h1_Eta","#eta distribution for reconstructed",etabins, etalow, etaup);
  h1_Eta->GetXaxis()->SetTitle("#eta");
  h1_Eta->GetYaxis()->SetTitle("counts");
  TotalNBins+=etabins;
  
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
       
  h2_EtaPt = new TH2F("h2_EtaPt","Cluster Energy vs ",etabins, etalow, etaup, ptbins, ptlow, ptup);
  h2_EtaPt->GetXaxis()->SetTitle("E [GeV]");
  h2_EtaPt->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  TotalNBins+=etabins*ptbins;

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

  //TFile *f = OpenFile(1); 
  //TDirectory::TContext context(f);
    
  fOutput->Add(h1_nClusters);
  fOutput->Add(h1_zvtx);
  fOutput->Add(h1_trigger);
  fOutput->Add(h1_E);
  fOutput->Add(h1_Phi);
  fOutput->Add(h2_PiMotherID);
  fOutput->Add(h2_GaMotherID);
  fOutput->Add(h3_gE_RecTruth);
  fOutput->Add(h3_gE_RecTruth_ncellscut);
  fOutput->Add(h1_Pi0TruthPt);
  fOutput->Add(h1_PriPi0TruthPt);
  fOutput->Add(h1_Pi0TruthPtEmcal);
  fOutput->Add(h1_PriPi0TruthPtEmcal);
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
  fOutput->Add(h1_Eta);
  fOutput->Add(h1_Chi2);
  fOutput->Add(h1_nTrkMatch);
  fOutput->Add(h1_ClusterDisp);
  fOutput->Add(h2_Ellipse);
  fOutput->Add(h2_EtaPt);
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
void AliAnalysisTaskSDMGammaMC::UserExec(Option_t *) 
{
  // Main loop Called for each event

  AliMCEvent *mcEvent = MCEvent();  
  Bool_t isMC = bool(mcEvent);//is this the right way to do this? 
  if (!mcEvent){
    cout << "no MC event" << endl;
    return;
  }
  
  TRandom3 randy; randy.SetSeed(0);

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
  
  Int_t iTrigger = 0;
  if (fEsdEv)       iTrigger = fEsdEv->GetHeader()->GetL0TriggerInputs();
  else if (fAodEv)  iTrigger = fAodEv->GetHeader()->GetL0TriggerInputs();
  //h1_trigger->Fill(iTrigger);
  
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
  Double_t vertZ=0.0;
  if (fEsdEv)       vertZ = fEsdEv->GetPrimaryVertex()->GetZ();
  else if (fAodEv)  vertZ = fAodEv->GetPrimaryVertex()->GetZ();    

  h1_zvtx->Fill(vertZ);
  //zvertex cut:
  if(fabs(vertZ)>10.0)
    return;
  
  h1_nClusters->Fill(nclusters);

  //cout << iskip << " " << izvtx << " " << imult << endl;  
  //cout << "GetNumberOfVertices(): " << fAodEv->GetNumberOfVertices() << endl;



  //######################### ~~~~~~~~~~~ ##################################
  //######################### STARTING MC ##################################
  //######################### ~~~~~~~~~~~ ##################################
  
  if(isMC){
    int isPrimary  = 0;
    int isK0sDecay = 0;

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

    for (Int_t iTrack = 0; iTrack<nTracksMC; ++iTrack) {
      AliMCParticle *mcP = static_cast<AliMCParticle*>(mcEvent->GetTrack(iTrack));
      if (!mcP)
	continue;            
      
      if(iTrack<nPTracksMC)  isPrimary = 1;
      else                   isPrimary = 0;
      
      if(mcP->PdgCode() == 22){
	if(isPrimary==1){
	  if(mcP->GetMother()>-1)
	    h2_GaMotherID->Fill(( (AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()) )->PdgCode(), 1);
	  else
	    h2_GaMotherID->Fill(0.0,1);
	}
	else
	  h2_GaMotherID->Fill(( (AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()) )->PdgCode(), 2);	
      }
      
      // it's a pion !! 
      if(mcP->PdgCode() != 111)
	continue;
      
      isK0sDecay = 0;
      if(mcP->GetMother()>-1){
	if( ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() ==  310 ||
	    ((AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()))->PdgCode() == -310  )
	  isK0sDecay = 1;
      }      
      
      // primary particle
      //Double_t dR_vtx = TMath::Sqrt((mcP->Xv()-evtVtx->GetX())*(mcP->Xv()-evtVtx->GetX()) + 
      //			    (mcP->Yv()-evtVtx->GetY())*(mcP->Yv()-evtVtx->GetY()));
      //if(dR_vtx <= 0.01)  isPrimary = 1;
      //else            isPrimary = 0;
      
      
      if(isPrimary==1){
	if(mcP->GetMother()>-1)
	  h2_PiMotherID->Fill(( (AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()) )->PdgCode(), 1);
	else
	  h2_PiMotherID->Fill(0.0,1);
      }
      else
	h2_PiMotherID->Fill(( (AliMCParticle*)mcEvent->GetTrack(mcP->GetMother()) )->PdgCode(), 2);
      
      h1_Pi0TruthPt    ->Fill(mcP->Pt());
      h2_Pi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());
      
      if(isPrimary==1){
	h1_PriPi0TruthPt    ->Fill(mcP->Pt());
	h2_PriPi0TruthPhiEta->Fill(mcP->Phi(),mcP->Eta());
      }
      
      if(mcP->Eta()<-1.0 || mcP->Eta()>1.0)
	continue;
      
      Int_t DecayPhotonLabel[2] = {mcP->GetDaughterFirst(),
				   mcP->GetDaughterLast() };
      
      if (DecayPhotonLabel[0]<0)  continue;
      if (DecayPhotonLabel[1]<0)  DecayPhotonLabel[1]=DecayPhotonLabel[0];      
      if (DecayPhotonLabel[1]-DecayPhotonLabel[0] != 1)  continue;
      
      bool bacc = true;
      bool binp = true;
      bool isConv[2] = {1,1};
      Int_t convIndices[2][2] = { {-1,-1},{-1,-1} };
      Double_t eta_d[2] = {0.0,0.0};
      Double_t phi_d[2] = {0.0,0.0};
      Int_t daughter_index = -1;
      for (Int_t iPhoton=DecayPhotonLabel[0];iPhoton<=DecayPhotonLabel[1];++iPhoton){
	if(iPhoton==DecayPhotonLabel[0]) daughter_index=0;
	else                             daughter_index=1;
        const AliMCParticle *dmc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(iPhoton));
	eta_d[daughter_index] = dmc->Eta();
	phi_d[daughter_index] = dmc->Phi();
        if(!(dmc->PdgCode()==22))             binp = false;
        if(!(eta_d[daughter_index]>fEtamin && eta_d[daughter_index]<fEtamax && 
	     phi_d[daughter_index]>fPhimin && phi_d[daughter_index]<fPhimax   ))   bacc = false;	
	
	if( ((TParticle*)dmc->Particle())->GetNDaughters() != 2 )  isConv[daughter_index] = 0;
	else{//if photon has 2 daughters. 
	  
	  Int_t dd1 = dmc->GetDaughterFirst();
	  Int_t dd2 = dmc->GetDaughterLast();
	  if (dd2-dd1 != 1)  cout << "How can this happen???? " << endl;
	  const AliMCParticle *dd1mc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dd1));
	  const AliMCParticle *dd2mc = static_cast<const AliMCParticle *>(mcEvent->GetTrack(dd2));
	  if( dd1mc->PdgCode() != -dd2mc->PdgCode() )
	    isConv[daughter_index] = 0;
	  else if( TMath::Abs(dd1mc->PdgCode())!=11 )
	    isConv[daughter_index] = 0;
	  if(isConv[daughter_index]==1){//store the e+e- indices...
	    convIndices[daughter_index][0] = dd1;
	    convIndices[daughter_index][1] = dd2;
	  }//if this photon converted. 
	}//close else-if photon has 2 daughters.
      }//loop over 2 decay photons (iPhoton)
      
      if(binp && bacc){// 2 Photons hit the EMCAL! 

	h1_Pi0TruthPtEmcal    ->Fill(mcP->Pt());
	h2_Pi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());		
	
	if(isPrimary==1){
	  h1_PriPi0TruthPtEmcal    ->Fill(mcP->Pt());
	  h2_PriPi0TruthPhiEtaEmcal->Fill(mcP->Phi(),mcP->Eta());
	}
	
	Int_t PhotonClusterMatch[2][3]  = { {0,-1,-1},
					    {0,-1,-1} };
	Int_t PhotonElectronMatch[2][6] = { {0,-1,-1,-1,-1,-1},
					    {0,-1,-1,-1,-1,-1} };
	
	for(int iCluster=0; iCluster<nclusters; iCluster++) {	    
	  
	  AliESDCaloCluster* esdCluster=NULL;
	  AliAODCaloCluster* aodCluster=NULL;
	  if (fEsdEv)       esdCluster = fEsdEv->GetCaloCluster(iCluster); // pointer to EMCal cluster
	  else if (fAodEv)  aodCluster = fAodEv->GetCaloCluster(iCluster); // pointer to EMCal cluster
	  
	  Double_t clustMC_phi, clustMC_eta;	  
	  if(fEsdEv){	    
	    if(esdCluster->IsEMCAL()){		

	      if(!isGoodEsdCluster(esdCluster))
		continue;
	      
	      Float_t pos[3] = {0,0,0};
	      esdCluster->GetPosition(pos);
	      TVector3 vpos(pos);
	      h1_Phi->Fill(vpos.Phi());
	      clustMC_phi = vpos.Phi();
	      clustMC_eta = vpos.Eta();
	      
	      Double_t dR = TMath::Sqrt((eta_d[daughter_index]-clustMC_eta)*(eta_d[daughter_index]-clustMC_eta) + 
					(phi_d[daughter_index]-clustMC_phi)*(phi_d[daughter_index]-clustMC_phi));
	      h1_dR_RealMC->Fill(dR);
	      //matches_pion_photon = 0;
	      //if(dR<=0.04) matches_pion_photon = 1;
	      
	      TArrayI *TruthLabelsA = esdCluster->GetLabelsArray();
	      if(TruthLabelsA){
		for(int itl=0; itl<TruthLabelsA->GetSize(); itl++){
		  
		  for(int iPhoton=0; iPhoton<2; iPhoton++){
		    if(TruthLabelsA->At(itl)==DecayPhotonLabel[iPhoton]){
		      PhotonClusterMatch[iPhoton][0] = 1;
		      PhotonClusterMatch[iPhoton][1] = DecayPhotonLabel[iPhoton];
		      PhotonClusterMatch[iPhoton][2] = iCluster;
		    }
		  }//loop over truth labels.
		  
		  AliMCParticle *elecCandidate = (AliMCParticle*)(mcEvent->GetTrack(TruthLabelsA->At(itl)));
		  if(TMath::Abs(elecCandidate->PdgCode())==11){//if we have an electron...
		    Int_t elecMother_index = elecCandidate->GetMother();
		    if(elecMother_index>1 && elecMother_index<nTracksMC){
		      AliMCParticle *elecMother   = (AliMCParticle*)(mcEvent->GetTrack(elecMother_index));
		      if( TMath::Abs(elecMother->PdgCode())==22 ){//if the e's mother is a photon...
			Int_t elecGrandMother_index = elecMother->GetMother();
			if(elecGrandMother_index==iTrack){//if the e's gMother is THE pi0 in question...
			  AliMCParticle *elecGrandMother = (AliMCParticle*)(mcEvent->GetTrack(elecGrandMother_index));
			  if( TMath::Abs(elecGrandMother->PdgCode())!=111 ) cout << "|| This can't happen!!  A pion is a pion is a pion is a pion... ||" << endl;
			  
			  for(int iPhoton=0; iPhoton<2; iPhoton++){
			    //if(convIndices[iPhoton][0]==elecMother_index){
			    if(convIndices[iPhoton][0]==TruthLabelsA->At(itl) || convIndices[iPhoton][1]==TruthLabelsA->At(itl)){
			      if(PhotonElectronMatch[iPhoton][1] == DecayPhotonLabel[iPhoton]) PhotonElectronMatch[iPhoton][0] = 2;
			      else                                                             PhotonElectronMatch[iPhoton][0] = 1;
			      PhotonElectronMatch[iPhoton][1] = DecayPhotonLabel[iPhoton];
			      if(PhotonElectronMatch[iPhoton][2]==-1) PhotonElectronMatch[iPhoton][2] = iCluster;//first cluster
			      else                                    PhotonElectronMatch[iPhoton][3] = iCluster;//second cluster
			      if(PhotonElectronMatch[iPhoton][2]==-1) PhotonElectronMatch[iPhoton][4] = elecCandidate->PdgCode();
			      else                                    PhotonElectronMatch[iPhoton][5] = elecCandidate->PdgCode();
			    }			    
			  }//loop over both decay photons
			  
			}//if it's THE pi0.
		      }//if we have a photon.	 
		    }//if we have an electron.
		  }//if the candidate has a real mother.
		  
		}//itl (TruthLabel loop)
	      }//if(TruthLabelsA exists)
	    }//if(isEMCal)
	  }//if(esdEv)
	}//loop over clusters. 

	
	for(int iPhoton=0; iPhoton<2; iPhoton++){
	  AliMCParticle *truthP = (AliMCParticle*)(mcEvent->GetTrack(DecayPhotonLabel[iPhoton]));
	  if(!truthP)
	    continue;
	  if(PhotonClusterMatch[iPhoton][0]==1){
	    AliESDCaloCluster *esdCluster = fEsdEv->GetCaloCluster(PhotonClusterMatch[iPhoton][2]);
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h3_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 1);
	    if(esdCluster->GetNCells()>=2)
	      h3_gE_RecTruth_ncellscut->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 1);
	  }
	  else if(PhotonElectronMatch[iPhoton][0]==2 && PhotonElectronMatch[iPhoton][2] == PhotonElectronMatch[iPhoton][3]){//merged conv photon
	    AliESDCaloCluster *esdCluster = fEsdEv->GetCaloCluster(PhotonElectronMatch[iPhoton][2]);
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h3_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 2);
	    if(esdCluster->GetNCells()>=2)
	      h3_gE_RecTruth_ncellscut->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 2);
	  }
	  else if(PhotonElectronMatch[iPhoton][0]==2){//non-merged conv photon (but both hit emcal)
	    AliESDCaloCluster *esdCluster = fEsdEv->GetCaloCluster(PhotonElectronMatch[iPhoton][2]);
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h3_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 4);
	    if(esdCluster->GetNCells()>=2)
	      h3_gE_RecTruth_ncellscut->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 4);
	                       esdCluster = fEsdEv->GetCaloCluster(PhotonElectronMatch[iPhoton][3]);
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h3_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 4);
	    if(esdCluster->GetNCells()>=2)
	      h3_gE_RecTruth_ncellscut->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 4);
	  }
	  else if(PhotonElectronMatch[iPhoton][0]==1){//non-merged conv photon (one missed emcal)
	    AliESDCaloCluster *esdCluster = fEsdEv->GetCaloCluster(PhotonElectronMatch[iPhoton][2]);
	    recalScale = PrivateEnergyRecal(esdCluster->E(), fRecalibrator);
	    h3_gE_RecTruth->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 3);
	    if(esdCluster->GetNCells()>=2)
	      h3_gE_RecTruth_ncellscut->Fill(recalScale*esdCluster->E(), truthP->E()/(recalScale*esdCluster->E()), 3);
	  }
	}//loop over decay photons (iPhoton).	  
	
      }// 2 Photons pointed at the EMCAL!       
    }//for(nTracksMC) ie. Truth Pion loop. 
    
  }//if(isMC)
  
  //######################### ~~~~~~~~~~~~ ##################################
  //######################### DONE WITH MC ##################################
  //######################### ~~~~~~~~~~~~ ##################################
  
  
  
  // NEW HISTO should be filled before this point, as PostData puts the
  // information for this iteration of the UserExec in the container
  PostData(1, fOutput);
  }

//________________________________________________________________________
void AliAnalysisTaskSDMGammaMC::Terminate(Option_t *) //specify what you want to have done
{
  // Called once at the end of the query.
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskSDMGammaMC::GetZvtxBin(Double_t vertZ)
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
Int_t AliAnalysisTaskSDMGammaMC::GetMultBin(Int_t mult){

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
Int_t AliAnalysisTaskSDMGammaMC::isGoodEsdCluster(AliESDCaloCluster* esdclust){

  int pass = 1;
  int nMinCells  = 1;
  double MinE    = 0.4;
  //double MinErat = 0;
  //double MinEcc  = 0;
  
  if (!esdclust)
    pass = 0;    
  if (!esdclust->IsEMCAL()) 
    pass = 0;//removes ~70% of clusters.
  if (esdclust->E()<MinE)
    pass = 0;//does nothing
  if (esdclust->GetNCells()<nMinCells)
    pass = 0;//does nothing
  //if (GetMaxCellEnergy(esdclust)/esdclust->E()<MinErat)
  //pass = 0;
  //if (esdclust->Chi2()<MinEcc) // eccentricity cut
  //pass = 0;//this is always -1.
    
  /*
  //This cuts out more than just 1 cell clusters
  //and drains the statistics badly.  
  //haven't figured out what it does yet. 
  if(esdclust->GetM20()<0.02)
  pass = 0;
  */
  //if(esdclust->GetM02()<0.1)
  //  pass = 0;
  //if(esdclust->GetM02()>0.5)
  //  pass = 0;
  //if(esdclust->GetNCells()<2)
  //  pass = 0;    

  Float_t pos[3] = {0,0,0};
  esdclust->GetPosition(pos);
  TVector3 clusterPosition(pos);
  if(clusterPosition.Eta()<fEtamin || clusterPosition.Eta()>fEtamax || 
     clusterPosition.Phi()<fPhimin || clusterPosition.Phi()>fPhimax  )
    pass = 0;
  clusterPosition.Delete();

  //doing this by hand now... 
  //if(!esdclust->GetNTracksMatched()==0)
  //pass = 0;
  
  return pass;
}

//________________________________________________________________________
Int_t AliAnalysisTaskSDMGammaMC::isGoodAodCluster(AliAODCaloCluster* aodclust){

  int pass = 1;
  int nMinCells  = 1;
  double MinE    = 0.4;
  //double MinErat = 0;
  //double MinEcc  = 0;
  
  if (!aodclust)
    pass = 0;    
  if (!aodclust->IsEMCAL()) 
    pass = 0;//removes ~70% of clusters.
  if (aodclust->E()<MinE)
    pass = 0;//does nothing
  if (aodclust->GetNCells()<nMinCells)
    pass = 0;//does nothing
  //if (GetMaxCellEnergy(aodclust)/aodclust->E()<MinErat)
  //pass = 0;
  //if (aodclust->Chi2()<MinEcc) // eccentricity cut
  //pass = 0;//this is always -1.
    
  /*
  //This cuts out more than just 1 cell clusters
  //and drains the statistics badly.  
  //haven't figured out what it does yet. 
  if(aodclust->GetM20()<0.02)
  pass = 0;
  if(aodclust->GetM02()<0.02)
  pass = 0;
  */
  //if(aodclust->GetM02()<0.1)
  //  pass = 0;
  //if(aodclust->GetM02()>0.5)
  //  pass = 0;
  //if(aodclust->GetNCells()<2)
  //  pass = 0;    

  Float_t pos[3] = {0,0,0};
  aodclust->GetPosition(pos);
  TVector3 clusterPosition(pos);
  if(clusterPosition.Eta()<fEtamin || clusterPosition.Eta()>fEtamax || 
     clusterPosition.Phi()<fPhimin || clusterPosition.Phi()>fPhimax  )
    pass = 0;
  clusterPosition.Delete();

  //if(!aodclust->GetNTracksMatched()==0)
  //pass = 0;
  
  return pass;
}
 
//________________________________________________________________________
Double_t AliAnalysisTaskSDMGammaMC::getDeltaPhi(TLorentzVector p1, TLorentzVector p2){

  double dphi = p1.Phi() - p2.Phi();

  if(dphi<0.5*TMath::Pi())  
    dphi = dphi + 2.0*TMath::Pi();

  if(dphi>1.5*TMath::Pi())  
    dphi = dphi - 2.0*TMath::Pi();

  return dphi;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSDMGammaMC::getDeltaEta(TLorentzVector p1, TLorentzVector p2){

  double deta = p1.PseudoRapidity() - p2.PseudoRapidity();

  return deta;
}


//________________________________________________________________________
Double_t AliAnalysisTaskSDMGammaMC::PrivateEnergyRecal(Double_t energy, Int_t iCalib){
  
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

