#define ana_sgl_cxx
#include "ana_sgl.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <deque>
#include <cstdlib> 
#include <TRandom.h>

using namespace std;

using std::vector;

class etrk;

//Buffer for event mixing
static const int NBUF=100; //depth of buffer
static const int NMix=10; //# of events mixed (for +-)
//static const int NMix=2; //# of events mixed (for +-)


static const int NZBIN=10;
static const int NCENT=10;
int d_ibuf[NZBIN][NCENT];
vector<etrk> d_vep[NBUF][NZBIN][NCENT];
vector<etrk> d_vem[NBUF][NZBIN][NCENT];
  
static const unsigned int MAXPOOL=150;
//static const unsigned int MAXPOOL=50;
static const int MAX_TRY=3;
deque<etrk> d_poolp[NZBIN][NCENT];
deque<etrk> d_poolm[NZBIN][NCENT]; 
  
void reshuffle_buffer(vector<etrk> &ve,
		      deque<etrk> &pool){
  //If there is not enough electron in the pool, give up
  unsigned int ne = ve.size();
  unsigned int poolsize = pool.size();
  if(poolsize < ne) {
    cout <<" pool size="<<poolsize<<" ne"<<ne<<endl;
    return;
  }
  for(unsigned int ie=0; ie < ne; ie++) {
    int j = rand()%poolsize;
    ve[ie] = pool[j];
  }
} 

void ana_sgl::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L ana_sgl.C
//      Root > ana_sgl t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      ana_event(ientry, jentry);
      // if (Cut(ientry) < 0) continue;
   }
}

//___________________________________________
  void ana_sgl::ana_init(char *outname){
  fout  = new TFile(outname,"recreate");
  char name[100];

  hdedx_pt = new TH2D("hdedx_pt","dedx vs. pT",100,0, 10, 100, 0, 200);
  hdedx_tof_elec_pt = new TH2D("hdedx_tof_elec_pt","dedx vs. pT with tof veto",100,0, 10, 100, 0, 200);
  hdedx_tof_all_pt = new TH2D("hdedx_tof_all_pt","dedx vs. pT with tof veto",100,0, 10, 100, 0, 200);
  hdedx_tof_elec_emc_pt = new TH2D("hdedx_tof_elec_emc_pt","dedx vs. pT with tof veto and EMC",100,0, 10, 100, 0, 200);
  hdedx_tof_all_emc_pt = new TH2D("hdedx_tof_all_emc_pt","dedx vs. pT with tof veto and EMC",100,0, 10, 100, 0, 200);
  hdedx_emc_pt = new TH2D("hdedx_emc_pt","dedx vs. pT with EMC",100,0, 10, 100, 0, 200);


  hbetatof_pt = new TH2D("hbetatof_pt","beta vs. pT",100,0, 10, 100, 0, 1);
  hbetatof_tof_elec_pt = new TH2D("hbetatof_tof_elec_pt","beta vs. pT with tof cut",100,0, 10, 100, 0, 1);
  hbetatof_tof_all_pt = new TH2D("hbetatof_tof_all_pt","beta vs. pT with tof cut",100,0, 10, 100, 0, 1);
  hbetatof_tof_elec_emc_pt = new TH2D("hbetatof_tof_elec_emc_pt","beta vs. pT with tof cut",100,0, 10, 100, 0, 1);
  hbetatof_tof_all_emc_pt = new TH2D("hbetatof_tof_all_emc_pt","beta vs. pT with tof cut and emc E/p",100,0, 10, 100, 0, 1);
  hbetatof_emc_pt = new TH2D("hbetatof_emc_pt","dedx vs. pT with EMC",100,0, 10, 100, 0, 200);

  hCentrality = new TH1F("hCentrality","hCentrality",100,0,100);
  hV0AC = new TH2F("hV0AC","V0AC",5000,0,15000,5000,0,15000);
  hV0AC_Ntrk = new TH2F("hV0AC_Ntrk","V0AC vs. Ntrk",5000,0,30000,5000,0,15000);
  hV0AC_NaccTrcklts  = new TH2F("hV0AC_NaccTrcklts","V0AC vs. Ntrklets",5000,0,30000,2500,0,5000);
  


  fBinWidth = 0.05;

  nHistos = (int)(10/fBinWidth);

  for(int i=0;i<nHistos;i++){
    sprintf(name,"hdedx_p%d",i);
    hdedx[i] = new TH1D(name, name, 200,0,200);

    sprintf(name,"hdedx_tof_elec_p%d",i);
    hdedx_tof_elec[i] = new TH1D(name, name, 200,0,200);

    sprintf(name,"hdedx_tof_all_p%d",i);
    hdedx_tof_all[i] = new TH1D(name, name, 200,0,200);


    sprintf(name,"hdedx_tof_elec_emc_p%d",i);
    hdedx_tof_elec_emc[i] = new TH1D(name, name, 200,0,200);

    sprintf(name,"hdedx_tof_all_emc_p%d",i);
    hdedx_tof_all_emc[i] = new TH1D(name, name, 200,0,200);

    sprintf(name,"hdedx_emc_p%d",i);
    hdedx_emc[i] = new TH1D(name, name, 200,0,200);

  }


  int nbinx=400;
  float max_x=20;
  float min_x=0.2;
  float binw = (TMath::Log(max_x)-TMath::Log(min_x))/nbinx;
  double xbin[401];
  for(int ii=0;ii<nbinx+1;ii++){
    xbin[ii] = TMath::Exp(TMath::Log(min_x) + 0.5*binw+binw*ii);
  }


  //// my histo 
  fEventStat = new TH1D("EventStat ","EventStat ",7,0,7);
  fEvent = new TH1D("Event","Number of Events",   60,0,60);
  fdEdXvsPt = new TH2D("dEdXvsPt","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  fdEdXnSigmaElecvsPt = new TH2D("fdEdXnSigmaElecvsPt","dE/dX normalized to electron vs. pT of TPC",
                                 nbinx, xbin, 2000, -10, 10);
  fTOFbetavsPt = new TH2D("fTOFbetavsPt","TOF beta vs. p", 400, 0, 20, 1200, 0, 1.2);
  fTOFnSigmaElecvsPt = new TH2D("fTOFnSigmaElecvsPt","TOF nsigma for electron", 400, 0, 20, 2000, -10, 10);

  fdEdXvsPtWithCut = new TH2D("dEdXvsPtWithCuts","dE/dX vs. PT of TPC", nbinx, xbin, 2000,0,200);
  

  fPtElec[0] = new TH1D("hPtElec0","elctron pt",250,0,5);
  fPtElec[1] = new TH1D("hPtElec1","elctron pt",250,0,5);


  fNelc_pos = new TH2D("hNelc_pos","# of electrons and positrons",40, -0.5, 39.5, 40, -0.5, 39.5);
  fNelc_all = new TH1D("hNelc_all","# of electrons and positrons",50, -0.5, 49.5);
  fNelc_all_pT =  new TH2D("hNelc_all_pT","# of electrons and positrons in each pT bin",10, 0, 5, 50, -0.5, 49.5);  

  for(int ic=0;ic<10;ic++){
    sprintf(name,"hNelc_pos_cent%d",ic);
    fNelc_pos_cent[ic] = new TH2D(name, name, 40, -0.5, 39.5, 40, -0.5, 39.5);
    sprintf(name,"hNelc_all_cent%d",ic);
    fNelc_all_cent[ic] = new TH1D(name, name, 50, -0.5, 49.5);
  }


  //// tree branch
  d_tree = new TTree("d_tree","single tree");
  d_tree->Branch("evt", &d_evt, "evt/I");
  d_tree->Branch("cent", &d_cent, "cent/D");
  d_tree->Branch("ntrk", &d_ntrk, "ntrk/D");
  d_tree->Branch("xvprim", &d_xvprim, "xvprim/D");
  d_tree->Branch("yvprim", &d_yvprim, "yvprim/D");
  d_tree->Branch("zvprim", &d_zvprim, "zvprim/D");
  d_tree->Branch("nacctrklets", &d_nacctrklets, "nacctrklets/D");
  d_tree->Branch("xres", &d_xres, "xres/D");
  d_tree->Branch("yres", &d_yres, "yres/D");
  d_tree->Branch("zres", &d_zres, "zres/D");
  d_tree->Branch("nelc", &d_nelc, "nelc/I");
  d_tree->Branch("px", d_px, "px[nelc]/D");
  d_tree->Branch("py", d_py, "py[nelc]/D");
  d_tree->Branch("pz", d_pz, "pz[nelc]/D");
  d_tree->Branch("p", d_p, "p[nelc]/D");
  d_tree->Branch("pt", d_pt, "pt[nelc]/D");
  d_tree->Branch("xv", d_xv, "xv[nelc]/D");
  d_tree->Branch("yv", d_yv, "yv[nelc]/D");
  d_tree->Branch("zv", d_zv, "zv[nelc]/D");
  d_tree->Branch("phi", d_phi, "phi[nelc]/D");
  d_tree->Branch("theta", d_theta, "theta[nelc]/D");
  d_tree->Branch("eta", d_eta, "eta[nelc]/D");
  d_tree->Branch("c", d_c, "c[nelc]/D");
  d_tree->Branch("nclusITS", d_nclusITS, "nclusITS[nelc]/D");
  d_tree->Branch("nclusTPC", d_nclusTPC, "nclusTPC[nelc]/D");
  d_tree->Branch("nclusTPCiter", d_nclusTPCiter, "nclusTPCiter[nelc]/D");
  d_tree->Branch("nfclusTPC", d_nfclusTPC, "nfclusTPC[nelc]/D");
  d_tree->Branch("nfclusTPCr", d_nfclusTPCr, "nfclusTPCr[nelc]/D");
  d_tree->Branch("nfclusTPCrFrac", d_nfclusTPCrFrac, "nfclusTPCrFrac[nelc]/D");
  d_tree->Branch("TPCsignalN", d_TPCsignalN, "TPCsignalN[nelc]/D");
  d_tree->Branch("TPCsignalNfrac", d_TPCsignalNfrac, "TPCsignalNfrac[nelc]/D");
  d_tree->Branch("TPCchi2cl", d_TPCchi2cl, "TPCchi2cl[nelc]/D");
  d_tree->Branch("trkstat", d_trkstat, "trkstat[nelc]/D");
  d_tree->Branch("nclsTRD", d_nclsTRD, "nclsTRD[nelc]/D");
  d_tree->Branch("TRDntracklets", d_TRDntracklets, "TRDntracklets[nelc]/D");
  d_tree->Branch("TRDpidquality", d_TRDpidquality, "TRDpidquality[nelc]/D");
  d_tree->Branch("TRDprobEle", d_TRDprobEle, "TRDprobEle[nelc]/D");
  d_tree->Branch("TRDprobPio", d_TRDprobPio, "TRDprobPio[nelc]/D");
  d_tree->Branch("impactXY", d_impactXY, "impactXY[nelc]/D");
  d_tree->Branch("impactZ", d_impactZ, "impactZ[nelc]/D");
  d_tree->Branch("tracklength", d_tracklength, "tracklength[nelc]/D");
  d_tree->Branch("ITSsignal", d_ITSsignal, "ITSsignal[nelc]/D");
  d_tree->Branch("ITSnsigmaEle", d_ITSnsigmaEle, "ITSnsigmaEle[nelc]/D");
  d_tree->Branch("ITSnsigmaPio", d_ITSnsigmaPio, "ITSnsigmaPio[nelc]/D");
  d_tree->Branch("ITSnsigmaMuo", d_ITSnsigmaMuo, "ITSnsigmaMuo[nelc]/D");
  d_tree->Branch("ITSnsigmaKao", d_ITSnsigmaKao, "ITSnsigmaKao[nelc]/D");
  d_tree->Branch("ITSnsigmaPro", d_ITSnsigmaPro, "ITSnsigmaPro[nelc]/D");
  d_tree->Branch("PIn", d_PIn, "PIn[nelc]/D");
  d_tree->Branch("TPCsignal", d_TPCsignal, "TPCsignal[nelc]/D");
  d_tree->Branch("TOFsignal", d_TOFsignal, "TOFsignal[nelc]/D");
  d_tree->Branch("TOFbeta", d_TOFbeta, "TOFbeta[nelc]/D");
  d_tree->Branch("TPCnSigmaEle", d_TPCnSigmaEle, "TPCnSigmaEle[nelc]/D");
  d_tree->Branch("TPCnSigmaPio", d_TPCnSigmaPio, "TPCnSigmaPio[nelc]/D");
  d_tree->Branch("TPCnSigmaMuo", d_TPCnSigmaMuo, "TPCnSigmaMuo[nelc]/D");
  d_tree->Branch("TPCnSigmaKao", d_TPCnSigmaKao, "TPCnSigmaKao[nelc]/D");
  d_tree->Branch("TPCnSigmaPro", d_TPCnSigmaPro, "TPCnSigmaPro[nelc]/D");
  d_tree->Branch("TOFnSigmaEle", d_TOFnSigmaEle, "TOFnSigmaEle[nelc]/D");
  d_tree->Branch("TOFnSigmaPio", d_TOFnSigmaPio, "TOFnSigmaPio[nelc]/D");
  d_tree->Branch("TOFnSigmaMuo", d_TOFnSigmaMuo, "TOFnSigmaMuo[nelc]/D");
  d_tree->Branch("TOFnSigmaKao", d_TOFnSigmaKao, "TOFnSigmaKao[nelc]/D");
  d_tree->Branch("TOFnSigmaPro", d_TOFnSigmaPro, "TOFnSigmaPro[nelc]/D");
  d_tree->Branch("E", d_E, "E[nelc]/D");
  d_tree->Branch("dphi", d_phi, "dphi[nelc]/D");
  d_tree->Branch("deta", d_eta, "deta[nelc]/D");

  d_tree->Branch("chi2ndf", d_chi2ndf, "chi2ndf[nelc]/D");

  cout<<"aaaaaaaaaaa"<<endl;
  d_ntpair = new TTree("ntpair","pair");
  d_ntpair->Branch("run", &d_run, "run/D");
  d_ntpair->Branch("event", &d_event, "event/I");
  d_ntpair->Branch("centrality",&d_centrality,"centrality/D");
  d_ntpair->Branch("prim_xv",&d_prim_xv,"prim_xv/D");
  d_ntpair->Branch("prim_yv",&d_prim_yv,"prim_yv/D");
  d_ntpair->Branch("prim_zv",&d_prim_zv,"prim_zv/D");
  d_ntpair->Branch("mass",&d_mass,"mass/D");
  d_ntpair->Branch("pxpair",&d_pxpair,"pxpair/D");
  d_ntpair->Branch("pypair",&d_pypair,"pypair/D");
  d_ntpair->Branch("pzpair",&d_pzpair,"pzpair/D");
  d_ntpair->Branch("ptpair",&d_ptpair,"ptpair/D");
  d_ntpair->Branch("epair",&d_epair,"epair/D");
  d_ntpair->Branch("eta",&d_etapair,"eta/D");
  d_ntpair->Branch("phi",&d_phipair,"phi/D");
  d_ntpair->Branch("cos",&d_cos,"cos/D");
  d_ntpair->Branch("phiv",&d_phiv,"phiv/D");
  d_ntpair->Branch("psi",&d_psi,"psi/D");
  d_ntpair->Branch("pairtype",&d_pairtype,"pairtype/I");
  d_ntpair->Branch("cent1",&d_cent1,"cent1/D");
  d_ntpair->Branch("xv1",&d_xv1,"xv1/D");
  d_ntpair->Branch("yv1",&d_yv1,"yv1/D");
  d_ntpair->Branch("zv1",&d_zv1,"zv1/D");
  d_ntpair->Branch("px1",&d_px1,"px1/D");
  d_ntpair->Branch("py1",&d_py1,"py1/D");
  d_ntpair->Branch("pz1",&d_pz1,"pz1/D");
  d_ntpair->Branch("pt1",&d_pt1,"pt1/D");
  d_ntpair->Branch("eta1",&d_eta1,"eta1/D");
  d_ntpair->Branch("phi1",&d_phi1,"phi1/D");
  d_ntpair->Branch("theta1",&d_theta1,"theta1/D");
  d_ntpair->Branch("tpc1",&d_tpc1,"tpc1/D");
  d_ntpair->Branch("ntpc_ele1",&d_ntpc_ele1,"ntpc_ele1/D");
  d_ntpair->Branch("ntpc_pio1",&d_ntpc_pio1,"ntpc_pio1/D");
  d_ntpair->Branch("ntpc_kao1",&d_ntpc_kao1,"ntpc_kao1/D");
  d_ntpair->Branch("ntpc_pro1",&d_ntpc_pro1,"ntpc_pro1/D");
  d_ntpair->Branch("beta1",&d_beta1,"beta1/D");
  d_ntpair->Branch("ntof_ele1",&d_ntof_ele1,"ntof_ele1/D");
  d_ntpair->Branch("ntof_pio1",&d_ntof_pio1,"ntof_pio1/D");
  d_ntpair->Branch("ntof_kao1",&d_ntof_kao1,"ntof_kao1/D");
  d_ntpair->Branch("ntof_pro1",&d_ntof_pro1,"ntof_pro1/D");
  d_ntpair->Branch("its1",&d_its1,"its1/D");
  d_ntpair->Branch("nits1",&d_nits1,"nits1/D");
  d_ntpair->Branch("ntpc1",&d_ntpc1,"ntpc1/D");
  d_ntpair->Branch("e1",&d_e1,"e1/D");
  d_ntpair->Branch("dphi1",&d_dphi1,"dphi1/D");
  d_ntpair->Branch("deta1",&d_deta1,"deta1/D");
  d_ntpair->Branch("dcaxy1",&d_dcaxy1,"dcaxy1/D");
  d_ntpair->Branch("dcaz1",&d_dcaz1,"dcaz1/D");
  d_ntpair->Branch("conv1",&d_conv1,"conv1/I");




  d_ntpair->Branch("cent2",&d_cent2,"cent2/D");
  d_ntpair->Branch("xv2",&d_xv2,"xv2/D");
  d_ntpair->Branch("yv2",&d_yv2,"yv2/D");
  d_ntpair->Branch("zv2",&d_zv2,"zv2/D");
  d_ntpair->Branch("px2",&d_px2,"px2/D");
  d_ntpair->Branch("py2",&d_py2,"py2/D");
  d_ntpair->Branch("pz2",&d_pz2,"pz2/D");
  d_ntpair->Branch("pt2",&d_pt2,"pt2/D");
  d_ntpair->Branch("eta2",&d_eta2,"eta2/D");
  d_ntpair->Branch("phi2",&d_phi2,"phi2/D");
  d_ntpair->Branch("theta2",&d_theta2,"theta2/D");
  d_ntpair->Branch("tpc2",&d_tpc2,"tpc2/D");
  d_ntpair->Branch("ntpc_ele2",&d_ntpc_ele2,"ntpc_ele2/D");
  d_ntpair->Branch("ntpc_pio2",&d_ntpc_pio2,"ntpc_pio2/D");
  d_ntpair->Branch("ntpc_kao2",&d_ntpc_kao2,"ntpc_kao2/D");
  d_ntpair->Branch("ntpc_pro2",&d_ntpc_pro2,"ntpc_pro2/D");
  d_ntpair->Branch("beta2",&d_beta2,"beta2/D");
  d_ntpair->Branch("ntof_ele2",&d_ntof_ele2,"ntof_ele2/D");
  d_ntpair->Branch("ntof_pio2",&d_ntof_pio2,"ntof_pio2/D");
  d_ntpair->Branch("ntof_kao2",&d_ntof_kao2,"ntof_kao2/D");
  d_ntpair->Branch("ntof_pro2",&d_ntof_pro2,"ntof_pro2/D");
  d_ntpair->Branch("its2",&d_its2,"its2/D");
  d_ntpair->Branch("nits2",&d_nits2,"nits2/D");
  d_ntpair->Branch("ntpc2",&d_ntpc2,"ntpc2/D");
  d_ntpair->Branch("e2",&d_e2,"e2/D");
  d_ntpair->Branch("dphi2",&d_dphi2,"dphi2/D");
  d_ntpair->Branch("deta2",&d_deta2,"deta2/D");
  d_ntpair->Branch("dcaxy2",&d_dcaxy2,"dcaxy2/D");
  d_ntpair->Branch("dcaz2",&d_dcaz2,"dcaz2/D");
  d_ntpair->Branch("conv2",&d_conv2,"conv2/I");

  cout<<"aaaaaaaaaaa"<<endl;
  for(int i=0;i<7;i++){
    for(int j=0;j<11;j++){
      sprintf(name,"hmasspt_cent%d_pair%d",j,i);
      hmasspt[i][j] = new TH2D(name, name, 500, 0, 5, 500, 0, 5);
      sprintf(name,"hmasspt_weight_cent%d_pair%d",j,i);
      hmasspt_weight[i][j] = new TH2D(name, name, 500, 0, 5, 500, 0, 5);
    }
  }

  vep.clear();  
  vem.clear();
  vep_tmp.clear();  
  vem_tmp.clear();

  d_evt = 0;
  d_event = 0;

  simflag = false;
  d_conv_flag = false;
  
  for(int i=0;i<10;i++){
    nelec_pos[i] = 0;
  }

  magnetic_field_mm = true; 

  d_flag_tof_cut = false;
  d_flag_emc_cut =false;
  d_flag_phiv = false;
  d_tpc_dedx_low = 75;
  d_tpc_dedx_high = 90;
  d_tof_low = -3;
  d_tof_high = 3;
  d_emc_low = 0.7;
  d_emc_high = 1.3;
  d_phiv_cut = 0.6;

  d_flag_kaon_veto = false;
  d_flag_proton_veto = false;
  d_dedx_kaon_veto_low = -2;
  d_dedx_kaon_veto_high = 2;
  d_dedx_proton_veto_low = -2;
  d_dedx_proton_veto_high = 2;


  cout<<"ana_init end:"<<endl;
}

//___________________________________________
void ana_sgl::ana_end(void){

  fout->cd();
  /*
  hdedx_pt->Write();
  hdedx_tof_elec_pt->Write();
  hdedx_tof_all_pt->Write();
  hdedx_tof_elec_emc_pt->Write();
  hdedx_tof_all_emc_pt->Write();
  hdedx_emc_pt->Write();

  hbetatof_pt->Write();
  hbetatof_tof_elec_pt->Write();
  hbetatof_tof_all_pt->Write();
  hbetatof_tof_elec_emc_pt->Write();
  hbetatof_tof_all_emc_pt->Write();
  hbetatof_emc_pt->Write();

  for(int i=0;i<nHistos;i++){
    hdedx[i]->Write();
    hdedx_tof_elec[i]->Write();
    hdedx_tof_all[i]->Write();
    hdedx_tof_elec_emc[i]->Write();
    hdedx_tof_all_emc[i]->Write();
    hdedx_emc[i]->Write();
  }

  */

  if(simflag==false){
    fEventStat->Write();
    fEvent->Write();
    fdEdXvsPt->Write();
    fdEdXnSigmaElecvsPt->Write();
    fTOFbetavsPt->Write();
    fTOFnSigmaElecvsPt->Write();
    fdEdXvsPtWithCut->Write();
    fPtElec[0]->Write();
    fPtElec[1]->Write();
  }
  fNelc_pos->Write();
  fNelc_all->Write();
  for(int i=0;i<10;i++){
    fNelc_pos_cent[i]->Write();
    fNelc_all_cent[i]->Write();
  }
  fNelc_all_pT->Write();
  //d_tree->Write();
  hCentrality->Write();
  hV0AC->Write();
  hV0AC_Ntrk->Write();
  hV0AC_NaccTrcklts->Write();
  for(int j=0;j<11;j++){
    for(int i=0;i<7;i++){
      hmasspt[i][j]->Write();
      hmasspt_weight[i][j]->Write();
    }
  }

  //d_ntpair->Write();
  fout->Close();
}

//____________________________________________
void ana_sgl::loop_a_file(char *file){

  TFile *treefile = TFile::Open(file);
  TDirectory *d = (TDirectory*)treefile->Get("PWG3_dielectron");
  if(d==0){
    cout<<" PWG3_dielectron is not found "<<endl;
    return ; 
  }

  //TTree *tree = (TTree*)d->Get("t");
  TTree *tree = (TTree*)d->Get("tree_MultiDie_CENT1");
  if(tree == 0) {
    cout <<"tree is not found in "<<file<<endl;
    treefile->Close();
    return;
  }
  cout << file <<" is opened"<<endl;
  if(simflag==false){
    //add_histograms(treefile);
  }
  Init(tree);
  Loop();
  tree->Clear();
  d->Clear();
  tree->Delete();
  d->Delete();
  delete d;
  treefile->Close();
  delete treefile;



  cout <<"one file processed"<<endl;
}

//____________________________________________
void ana_sgl::ana_event(int ientry, int jentry){

  //select trigger class: 
  if(sel_trigger==2){
    if(kTriggerCent<100){
      return ;
    }
  }else if(sel_trigger==1){
    if(kTriggerCent%100<10){
      return ;
    }
  }else if(sel_trigger==0){
    if(kTriggerCent%10!=1){
      return ;
    }
  }

  //  if(fkRunNumber>138350){
  /*
  if(kMag>0){
    magnetic_field_mm = false;
  }else{
    magnetic_field_mm = true;
  }
  */
  if(fkRunNumber>169591){
    magnetic_field_mm = false;
  }else{
    magnetic_field_mm = true;
  }

  if(ientry%1000==0){
    cout<<" event processing "<<ientry<<" / "<<d_evt<<" / "<<d_event<<" : trigger "<<kTriggerCent<<endl;
  }
  
  d_nelc = 0;
  for(int i=0;i<10;i++){
    nelec_pos[i] = 0;
  }

  hCentrality->Fill(fkCentrality);
  hV0AC->Fill(fkV0C, fkV0A);
  hV0AC_Ntrk->Fill(fkNTrk, fkV0C+fkV0A);
  hV0AC_NaccTrcklts->Fill(fkNaccTrcklts, fkV0C+fkV0A);
  fill_to_tree_variables();
  
  for(int i=0;i<fkNPar;i++){

    if(simflag==false){
      
      //fdEdXvsPtWithCut->Fill(kPIn[i], kTPCsignal[i]);

      if(GlobalTrackcut(i)==false){
	continue;
      }
      fill_histograms(i);
    }

    
    fill_to_tree_track_variables(i);

    //if(PairTrackcut(i)==false){
    //continue;
    //}

    if(fkCentrality>0 && fkCentrality<10){
      int iptbin = (int)(kP[i]/0.5);
      if(iptbin>=10){
	iptbin=9;
      }
      nelec_pos[iptbin]++;
    }

    if(kCharge[i]>0){
      etrk e = etrk(
		    fkCentrality, fkXvPrim, fkYvPrim, fkZvPrim,
		    kXv[i], kYv[i], kZv[i], 
		    kPx[i], kPy[i], kPz[i], kPt[i],
		    kEta[i], kPhi[i], kTheta[i],
		    kTPCsignal[i], kTOFbeta[i], 
		    kE[i], kDeltaPhi[i], kDeltaEta[i],
		    kTPCnSigmaEle[i], kTPCnSigmaPio[i], kTPCnSigmaKao[i], kTPCnSigmaPro[i],
		    kTOFnSigmaEle[i], kTOFnSigmaPio[i], kTOFnSigmaKao[i], kTOFnSigmaPro[i],
		    kITSsignal[i], kNclsITS[i], kNclsTPC[i], kLegDistXY[i], kLegDist[i]
		  );
      vep_tmp.push_back(e);
    }else{
      etrk e = etrk(
		    fkCentrality, fkXvPrim, fkYvPrim, fkZvPrim,
		    kXv[i], kYv[i], kZv[i], 
		    kPx[i], kPy[i], kPz[i], kPt[i],
		    kEta[i], kPhi[i], kTheta[i],
		    kTPCsignal[i], kTOFbeta[i], 
		    kE[i], kDeltaPhi[i], kDeltaEta[i],
		    kTPCnSigmaEle[i], kTPCnSigmaPio[i], kTPCnSigmaKao[i], kTPCnSigmaPro[i],
		    kTOFnSigmaEle[i], kTOFnSigmaPio[i], kTOFnSigmaKao[i], kTOFnSigmaPro[i],
		    kITSsignal[i], kNclsITS[i], kNclsTPC[i], kLegDistXY[i], kLegDist[i]
		    );
      vem_tmp.push_back(e);
    }
  }
  
  //////// fill to the tree //////////////////
  //d_tree->Fill();
  


  if(d_conv_flag==true){
    check_conversion_pairs(vep_tmp, vem_tmp);
  }

  check_ghost_pairs(vep_tmp);
  check_ghost_pairs(vem_tmp);
  randomize_pool(vep_tmp, vem_tmp);
  if(d_conv_flag==false){
    check_conversion_pairs(vep, vem);
  }

  fNelc_pos->Fill(vep.size(), vem.size());
  fNelc_all->Fill(vep.size()+vem.size());

  if(fkCentrality>0 && fkCentrality<10){
    for(int i=0;i<10;i++){
      fNelc_all_pT->Fill(0.25+0.5*i, nelec_pos[i]);
    }
  }


  calc_pair(vep, vem);

  int icent = (int)(fkCentrality/10.0);
  int izbin = (int)((fkZvPrim+10)/2.0);
  if(icent<0) icent=0;
  if(icent>=NCENT) icent=NCENT-1;
  if(izbin<0) izbin=0;
  if(izbin>=NZBIN) izbin=NZBIN-1;
  
  fNelc_pos_cent[icent]->Fill(vep.size(), vem.size());
  fNelc_all_cent[icent]->Fill(vep.size()+vem.size());



  d_vep[d_ibuf[izbin][icent]][izbin][icent].clear();
  vector<etrk>::iterator iep;
  for(iep = vep.begin();iep != vep.end();++iep) {
    d_vep[d_ibuf[izbin][icent]][izbin][icent].push_back(*iep);
    d_poolp[izbin][icent].push_back(*iep);
    if(d_poolp[izbin][icent].size()>MAXPOOL) {
      d_poolp[izbin][icent].pop_front();
    }
  }
  d_vem[d_ibuf[izbin][icent]][izbin][icent].clear();
  vector<etrk>::iterator iem;
  for(iem = vem.begin();iem != vem.end();++iem) {
    d_vem[d_ibuf[izbin][icent]][izbin][icent].push_back(*iem);
    d_poolm[izbin][icent].push_back(*iem);
    if(d_poolm[izbin][icent].size()>MAXPOOL) {
      d_poolm[izbin][icent].pop_front();
    }
  }
  // Update the buffer pointer
  d_ibuf[izbin][icent]++;
  if(d_ibuf[izbin][icent]>= NBUF) d_ibuf[izbin][icent]=0; 

  ///////////////////////////////////////////
  vem.clear();
  vep.clear();
  vem_tmp.clear();
  vep_tmp.clear();
  d_evt ++;
  d_event ++;
  
}  


//____________________________________________
bool ana_sgl::kTOFcut(int itrk){

  if(kTOFsignal[itrk]!=9999 && kTOFbeta[itrk]>0.3 && 
     (kTOFnSigmaEle[itrk]<3 && kTOFnSigmaEle[itrk]>-3)
     /*
     (kTOFnSigmaKao>3 || kTOFnSigmaKao<-3) &&
     (kTOFnSigmaPro>3 || kTOFnSigmaPro<-3) 
     */

     ){

    return true;
  }else{
    return false;
  }
}

//____________________________________________
bool ana_sgl::GlobalTrackcut(int itrk){


  if((
     (d_flag_kaon_veto == true && (kTPCnSigmaKao[itrk]<d_dedx_kaon_veto_low || kTPCnSigmaKao[itrk]>d_dedx_kaon_veto_high)) ||
     (d_flag_kaon_veto == false)
     )
     && (
	 (d_flag_proton_veto == true && (kTPCnSigmaPro[itrk]<d_dedx_proton_veto_low || kTPCnSigmaPro[itrk]>d_dedx_proton_veto_high)) ||
	 (d_flag_proton_veto == false)
	 )
     ){
    fdEdXvsPtWithCut->Fill(kPIn[itrk], kTPCsignal[itrk]);
  }

  if( kNclsTPC[itrk]>120 
      && kTPCsignal[itrk]>d_tpc_dedx_low && kTPCsignal[itrk]<d_tpc_dedx_high
      && (
	  (d_flag_tof_cut==true && kTOFnSigmaEle[itrk]>d_tof_low && kTOFnSigmaEle[itrk]<d_tof_high) ||
	  (d_flag_tof_cut==false)
	  )
      && (
	  (d_flag_pt_cut == true && kPt[itrk]>d_pt_cut_low && kPt[itrk]<d_pt_cut_high) ||
	  (d_flag_pt_cut == false)
	  )
      && (
	  (d_flag_kaon_veto == true && (kTPCnSigmaKao[itrk]<d_dedx_kaon_veto_low || kTPCnSigmaKao[itrk]>d_dedx_kaon_veto_high)) ||
	  (d_flag_kaon_veto == false)
	  )
      && (
	  (d_flag_proton_veto == true && (kTPCnSigmaPro[itrk]<d_dedx_proton_veto_low || kTPCnSigmaPro[itrk]>d_dedx_proton_veto_high)) ||
	  (d_flag_proton_veto == false)
	  )
      ){
    return true;
  }else{
    return false;
  }
}

/*
//____________________________________________
bool ana_sgl::PairTrackcut(int itrk){

  if( (kTOFnSigmaEle[itrk]<3 && kTOFnSigmaEle[itrk]>-3) &&
      (kTOFnSigmaKao[itrk]>3||kTOFnSigmaKao[itrk]<-3) &&
      (kTOFnSigmaPro[itrk]>3||kTOFnSigmaPro[itrk]<-3) //&&
      //(kTPCnSigmaEle[itrk]-1.65226*exp(-kPIn[itrk]*kPIn[itrk]*1.60890)+0.838)<1 &&
      //(kTPCnSigmaEle[itrk]-1.65226*exp(-kPIn[itrk]*kPIn[itrk]*1.60890)+0.838)>-1
      ){
    return true;
  }else{
    return false;
  }


//  if(
//     !(kTPCsignal>75-210*(kPIn-0.5) && 
//       kTPCsignal<100-190*(kPIn-0.5)) &&
//     !(kTPCsignal>130-120*(kPIn-0.5) && 
//       kTPCsignal<135-80*(kPIn-0.5)) &&
//     (kTPCnSigmaEle-1.65226*exp(-kPIn*kPIn*1.60890)+0.838)<1 &&
//      (kTPCnSigmaEle-1.65226*exp(-kPIn*kPIn*1.60890)+0.838)>-1
//     ){
//    return true;
//  }else{
//   return false;
//  }


}
*/

//____________________________________________
void ana_sgl::fill_histograms(int itrk){
  /*
  hdedx_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
  hbetatof_pt->Fill(kPIn[itrk], kTOFbeta[itrk]);

  if(kTOFcut(itrk)==true){
    hdedx_tof_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
    hbetatof_tof_pt->Fill(kPIn[itrk], kTOFbeta[itrk]);
  }

  int iptbin = (int)((kPIn[itrk])/fBinWidth);
  if(iptbin>=0 && iptbin<nHistos){
    hdedx[iptbin]->Fill(kTPCsignal[itrk]);
    if(kTOFcut(itrk)==true)hdedx_tof[iptbin]->Fill(kTPCsignal[itrk]);
  }
  */

  if(kCharge[itrk]>0){
    fPtElec[1]->Fill(kPIn[itrk]);
  }else{
    fPtElec[0]->Fill(kPIn[itrk]);
  }
  
  hdedx_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
  hbetatof_pt->Fill(kPIn[itrk], kTOFbeta[itrk]);

  if(kTOFnSigmaEle[itrk]>-3 && kTOFnSigmaEle[itrk]<3){
    hdedx_tof_elec_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
    hbetatof_tof_elec_pt->Fill(kPIn[itrk], kTOFbeta[itrk]);  
  
    if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
      hdedx_tof_elec_emc_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
      hbetatof_emc_pt->Fill(kPIn[itrk], kTOFbeta[itrk]);  
    }      

    if( (kTOFnSigmaKao[itrk]>3 || kTOFnSigmaKao[itrk]<-3) &&
	(kTOFnSigmaPro[itrk]>3 || kTOFnSigmaPro[itrk]<-3) 
	){
      hdedx_tof_all_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
      hbetatof_tof_all_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
      if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
	hdedx_tof_all_emc_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
	hbetatof_tof_all_emc_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
      }      
    }
  }

  if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
    hdedx_emc_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
    hbetatof_emc_pt->Fill(kPIn[itrk], kTPCsignal[itrk]);
  }

  
  

  int iptbin = (int)((kPIn[itrk])/fBinWidth);
  if(iptbin>=0 && iptbin<nHistos){
    hdedx[iptbin]->Fill(kTPCsignal[itrk]);
    //if(kTOFcut(itrk)==true)hdedx_tof[iptbin]->Fill(kTPCsignal[itrk]);

    if(kTOFnSigmaEle[itrk]>-3 && kTOFnSigmaEle[itrk]<3){
      hdedx_tof_elec[iptbin]->Fill(kTPCsignal[itrk]);

      if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
	hdedx_tof_elec_emc[iptbin]->Fill(kTPCsignal[itrk]);
      }
      if( (kTOFnSigmaKao[itrk]>3 || kTOFnSigmaKao[itrk]<-3) &&
	  (kTOFnSigmaPro[itrk]>3 || kTOFnSigmaPro[itrk]<-3) 
	  ){
	hdedx_tof_all[iptbin]->Fill(kTPCsignal[itrk]);
	if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
	  hdedx_tof_all_emc[iptbin]->Fill(kTPCsignal[itrk]);
	}
      }
    }
    if(kE[itrk]/kPIn[itrk]>0.7 && kE[itrk]/kPIn[itrk]<1.3){
      hdedx_emc[iptbin]->Fill(kTPCsignal[itrk]);
    }
  }

}

//____________________________________________
void ana_sgl::add_histograms(TFile *fin){

  TDirectory *d = (TDirectory*)fin->Get("PWG3_dielectron");
  if(d==0){
    cout<<" PWG3_dielectron is not found "<<endl;
    return ; 
  }
  fEventStat->Add((TH1D*)d->Get("hEventStat_MultiDie_CENT1"));
  TList *list = (TList*)d->Get("jpsi_QA_CENT1");
  TList *listQA = (TList*)list->FindObject("QAElectron");
  fEvent->Add((TH1D*)listQA->FindObject("Event"));
  fdEdXvsPt->Add((TH1D*)listQA->FindObject("dEdXvsPt"));
  fdEdXnSigmaElecvsPt->Add((TH1D*)listQA->FindObject("fdEdXnSigmaElecvsPt"));
  fTOFbetavsPt->Add((TH1D*)listQA->FindObject("fTOFbetavsPt"));
  fTOFnSigmaElecvsPt->Add((TH1D*)listQA->FindObject("fTOFnSigmaElecvsPt"));

  cout<<" ana_sgl::add_histograms "<<fin->GetName()<<" done "<<endl;

}

//_____________________________________________
void ana_sgl::fill_to_tree_variables(void){

  d_cent =  fkCentrality;
  d_ntrk= fkNTrk;
  d_xvprim= fkXvPrim;
  d_yvprim= fkYvPrim;
  d_zvprim= fkZvPrim;
  d_nacctrklets= fkNaccTrcklts;
  d_xres= fkXRes;
  d_yres= fkYRes;
  d_zres= fkZRes;

}
//_____________________________________________
void ana_sgl::fill_to_tree_track_variables(int itrk){
  

  d_px[d_nelc]= kPx[itrk];
  d_py[d_nelc]= kPy[itrk];
  d_pz[d_nelc]= kPz[itrk];
  d_p[d_nelc]= kP[itrk];
  d_pt[d_nelc]= kPt[itrk];
  d_xv[d_nelc]= kXv[itrk];
  d_yv[d_nelc]= kYv[itrk];
  d_zv[d_nelc]= kZv[itrk];
  d_phi[d_nelc]= kPhi[itrk];
  d_theta[d_nelc]= kTheta[itrk];
  d_eta[d_nelc]= kEta[itrk];
  d_c[d_nelc]= kCharge[itrk];
  d_nclusITS[d_nelc]= kNclsITS[itrk];
  d_nclusTPC[d_nelc]= kNclsTPC[itrk];
  d_nclusTPCiter[d_nelc]= kNclsTPCiter1[itrk];
  d_nfclusTPC[d_nelc]= kNFclsTPC[itrk];
  d_nfclusTPCr[d_nelc]= kNFclsTPCr[itrk];
  d_nfclusTPCrFrac[d_nelc]= kNFclsTPCrFrac[itrk];
  d_TPCsignalN[d_nelc]= kTPCsignalN[itrk];
  d_TPCsignalNfrac[d_nelc]= kTPCsignalNfrac[itrk];
  d_TPCchi2cl[d_nelc]= kTPCchi2Cl[itrk];
  d_trkstat[d_nelc]= kTrackStatus[itrk];
  d_nclsTRD[d_nelc]= kNclsTRD[itrk];
  d_TRDntracklets[d_nelc]= kTRDntracklets[itrk];
  d_TRDpidquality[d_nelc]= kTRDpidQuality[itrk];
  d_TRDprobEle[d_nelc]= kTRDprobEle[itrk];
  d_TRDprobPio[d_nelc]= kTRDprobPio[itrk];
  d_impactXY[d_nelc]= kImpactParXY[itrk];
  d_impactZ[d_nelc]= kImpactParZ[itrk];
  d_tracklength[d_nelc]= kTrackLength[itrk];
  d_ITSsignal[d_nelc]= kITSsignal[itrk];
  d_ITSnsigmaEle[d_nelc]= kITSnSigmaEle[itrk];
  d_ITSnsigmaPio[d_nelc]= kITSnSigmaPio[itrk];
  d_ITSnsigmaMuo[d_nelc]= kITSnSigmaMuo[itrk];
  d_ITSnsigmaKao[d_nelc]= kITSnSigmaKao[itrk];
  d_ITSnsigmaPro[d_nelc]= kITSnSigmaPro[itrk];
  d_PIn[d_nelc]= kPIn[itrk];
  d_TPCsignal[d_nelc]= kTPCsignal[itrk];
  d_TOFsignal[d_nelc]= kTOFsignal[itrk];
  d_TOFbeta[d_nelc]= kTOFbeta[itrk];
  d_TPCnSigmaEle[d_nelc]= kTPCnSigmaEle[itrk];
  d_TPCnSigmaPio[d_nelc]= kTPCnSigmaPio[itrk];
  d_TPCnSigmaMuo[d_nelc]= kTPCnSigmaMuo[itrk];
  d_TPCnSigmaKao[d_nelc]= kTPCnSigmaKao[itrk];
  d_TPCnSigmaPro[d_nelc]= kTPCnSigmaPro[itrk];
  d_TOFnSigmaEle[d_nelc]= kTOFnSigmaEle[itrk];
  d_TOFnSigmaPio[d_nelc]= kTOFnSigmaPio[itrk];
  d_TOFnSigmaMuo[d_nelc]= kTOFnSigmaMuo[itrk];
  d_TOFnSigmaKao[d_nelc]= kTOFnSigmaKao[itrk];
  d_TOFnSigmaPro[d_nelc]= kTOFnSigmaPro[itrk];

  d_chi2ndf[d_nelc]= kChi2NDF[itrk];
  d_E[d_nelc] = kE[itrk];
  d_dphi[d_nelc] = kDeltaPhi[itrk];
  d_deta[d_nelc] = kDeltaEta[itrk];

  d_nelc++;

}

//____________________________________________________________________
void ana_sgl::randomize_pool(vector<etrk> e1, vector<etrk> e2){
  
  int size1 = e1.size();
  int used_index[1000];
  for(int i=0;i<1000;i++){
    used_index[i] = -1;
  }
  for(int i=0;i<size1;i++){
    used_index[i] = 0;
  }

  for(unsigned int i=0;i<size1;i++){
    int j = (int)(gRandom->Uniform(0,size1));
    while(used_index[j]==1){
      j = (int)(gRandom->Uniform(0,size1));
    }
    if( (e1[j].ghost_flag==1) &&
	(
	 (d_conv_flag==true && e1[j].conv_flag==1) ||
	 (d_conv_flag==false)
	 )
	){
      vep.push_back(e1[j]);
    }
    used_index[j] = 1;
  }
  

  int size2 = e2.size();
  for(int i=0;i<1000;i++){
    used_index[i] = -1;
  }
  for(int i=0;i<size2;i++){
    used_index[i] = 0;
  }

  for(unsigned int i=0;i<size2;i++){
    int j = (int)(gRandom->Uniform(0,size2));
    while(used_index[j]==1){
      j = (int)(gRandom->Uniform(0,size2));
    }
    if( (e2[j].ghost_flag==1) &&
      (
       (d_conv_flag==true && e2[j].conv_flag==1) ||
       (d_conv_flag==false)
       )
	){
      vem.push_back(e2[j]);
    }
    used_index[j] = 1;
  }
}

//____________________________________________________________________
void ana_sgl::calc_pair(vector<etrk> vep, vector<etrk> vem){
  
  //vector<etrk>::iterator iep, iem; 

  vector<etrk>::iterator iep;
  vector<etrk>::iterator iem;
  
  //cout<<vep.size()<<" "<<vem.size()<<endl;

  ///unlike pairs
  for(iep=vep.begin(); iep!=vep.end(); ++iep){
    for(iem=vem.begin(); iem!=vem.end(); ++iem){
      if(PairTrackcut(iep, iem)==false) continue;
      fill_pair(iep, iem, 0);  
    }                                                                                                                                                      
  }                                                                                                                                                                                                    
  for(iep=vep.begin(); iep!=vep.end(); ++iep){
    vector<etrk>::iterator iep2=iep;
    ++iep2;
    for(iem=iep2; iem!=vep.end(); ++iem){
      if(PairTrackcut(iep, iem)==false) continue;
      fill_pair(iep, iem, 1);
    }
  }                                              
                                                                                                                                                         
  for(iep=vem.begin(); iep!=vem.end(); ++iep){
    vector<etrk>::iterator iep2=iep;
    ++iep2;
    for(iem=iep2; iem!=vem.end(); ++iem){
      if(PairTrackcut(iep, iem)==false) continue;
      fill_pair(iep, iem, 2);   
    }
  }

  int icent = (int)(fkCentrality/10.0);
  int izbin = (int)((fkZvPrim+10)/2.0);
  if(icent<0) icent=0;
  if(icent>=NCENT) icent=NCENT-1;
  if(izbin<0) izbin=0;
  if(izbin>=NZBIN) izbin=NZBIN-1;


  int nmixed;
  if(vep.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
	reshuffle_buffer(d_vem[ibuf][izbin][icent],d_poolm[izbin][icent]);
	ntry++;
      }
      for(iep=vep.begin(); iep!=vep.end(); ++iep){
	for(iem=d_vem[ibuf][izbin][icent].begin(); 
	    iem!=d_vem[ibuf][izbin][icent].end(); ++iem){
	  if(PairTrackcut(iep, iem)==false) continue;
	  fill_pair(iep,iem,3);
	}
      }
      ++nmixed;
    }//for(ibuf)
  }

  if(vem.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
	reshuffle_buffer(d_vep[ibuf][izbin][icent],d_poolp[izbin][icent]);
	ntry++;
      }
      for(iem=vem.begin(); iem!=vem.end(); ++iem){
	for(iep=d_vep[ibuf][izbin][icent].begin(); 
	    iep!=d_vep[ibuf][izbin][icent].end(); ++iep){
	  if(PairTrackcut(iep, iem)==false) continue;
	  fill_pair(iep,iem,4);
	}
      }
      ++nmixed;
    }//for(ibuf)
  }


  if(vep.size()>0) {
    //
    // Now mixed event for ++ pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
	reshuffle_buffer(d_vep[ibuf][izbin][icent],d_poolp[izbin][icent]);
	ntry++;
      }
      for(iep=vep.begin(); iep!=vep.end(); ++iep){
	for(iem=d_vep[ibuf][izbin][icent].begin(); 
	    iem!=d_vep[ibuf][izbin][icent].end(); ++iem){
	  if(PairTrackcut(iep, iem)==false) continue;
	  fill_pair(iep,iem,5);
	}
      }
      ++nmixed;
    }//for(ibuf)
  }

  if(vem.size()>0) {
    //
    // Now mixed event for +- pairs
    //
    nmixed = 0;
    for(int ibuf=0;(nmixed<NMix);ibuf++) {
      int ntry = 0;
      while(ntry<MAX_TRY) {
	reshuffle_buffer(d_vem[ibuf][izbin][icent],d_poolm[izbin][icent]);
	ntry++;
      }
      for(iem=vem.begin(); iem!=vem.end(); ++iem){
	for(iep=d_vem[ibuf][izbin][icent].begin(); 
	    iep!=d_vem[ibuf][izbin][icent].end(); ++iep){
	  if(PairTrackcut(iep, iem)==false) continue;
	  fill_pair(iep,iem,6);
	}
      }
      ++nmixed;
    }//for(ibuf)
  }
}


//___________________________________________________________________________
//void ana_sgl::fill_pair(etrk* iep, etrk* iem, int type){
//void ana_sgl::fill_pair(etrk iep, etrk iem, int type){
void ana_sgl::fill_pair(vector<etrk>::iterator iep, vector<etrk>::iterator iem, int type){


  d_pairtype = type;

  calc_vars(iep, iem, d_mass, d_phiv, d_pxpair, d_pypair, d_pzpair, 
	    d_ptpair, d_epair, d_phipair, d_etapair, d_cos, d_psi);

  if(type==0||type==1||type==2||type==3||type==5){
    d_centrality = iep->cent;
    d_prim_xv = iep->pxv;
    d_prim_yv = iep->pyv;
    d_prim_zv = iep->pzv;
  }else if(type==4 || type==6){
    d_centrality = iem->cent;
    d_prim_xv = iem->pxv;
    d_prim_yv = iem->pyv;
    d_prim_zv = iem->pzv;

  }

  
  //cout<<iep->px<<" "<<iem->px<<" "<<d_centrality<<endl;  
  
  d_cent1 = iep->cent;
  d_xv1 = iep->xv;
  d_yv1 = iep->yv;
  d_zv1 = iep->zv;
  d_px1 = iep->px;
  d_py1 = iep->py;
  d_pz1 = iep->pz;
  d_pt1 = iep->pt;
  d_eta1 = iep->eta;
  d_phi1 = iep->phi;
  d_theta1 = iep->theta;
  d_tpc1 = iep->tpc;
  d_ntpc_ele1 = iep->ntpc_ele;
  d_ntpc_pio1 = iep->ntpc_pio;
  d_ntpc_kao1 = iep->ntpc_kao;
  d_ntpc_pro1 = iep->ntpc_pro;
  d_ntof_ele1 = iep->ntof_ele;
  d_ntof_pio1 = iep->ntof_pio;
  d_ntof_kao1 = iep->ntof_kao;
  d_ntof_pro1 = iep->ntof_pro;
  d_its1 = iep->its;
  d_nits1 = iep->nits;
  d_ntpc1 = iep->ntpc;
  d_e1 = iep->e;
  d_dphi1 = iep->dphi;
  d_deta1 = iep->deta;
  d_dcaxy1 = iep->dcaxy;
  d_dcaz1 = iep->dcaz;
  d_conv1 = iep->conv_flag;

  d_cent2 = iem->cent;
  d_xv2 = iem->xv;
  d_yv2 = iem->yv;
  d_zv2 = iem->zv;
  d_px2 = iem->px;
  d_py2 = iem->py;
  d_pz2 = iem->pz;
  d_pt2 = iem->pt;
  d_eta2 = iem->eta;
  d_phi2 = iem->phi;
  d_theta2 = iem->theta;
  d_tpc2 = iem->tpc;
  d_ntpc_ele2 = iem->ntpc_ele;
  d_ntpc_pio2 = iem->ntpc_pio;
  d_ntpc_kao2 = iem->ntpc_kao;
  d_ntpc_pro2 = iem->ntpc_pro;
  d_ntof_ele2 = iem->ntof_ele;
  d_ntof_pio2 = iem->ntof_pio;
  d_ntof_kao2 = iem->ntof_kao;
  d_ntof_pro2 = iem->ntof_pro;
  d_its2 = iem->its;
  d_nits2 = iem->nits;
  d_ntpc2 = iem->ntpc;
  d_e2 = iem->e;
  d_dphi2 = iem->dphi;
  d_deta2 = iem->deta;
  d_dcaxy2 = iem->dcaxy;
  d_dcaz2 = iem->dcaz;
  d_conv2 = iem->conv_flag;
  
  d_run = fkRunNumber;
  //  d_event = d_evt;
  

  int icent = (int)(fkCentrality/10.0);
  if(icent<0) icent=0;
  if(icent>=NCENT) icent=NCENT-1;
  if(pair_cut()==true){
    hmasspt[type][icent]->Fill(d_mass, d_ptpair);
    hmasspt[type][10]->Fill(d_mass, d_ptpair);
  }

  ////////////// pt weighting for mixed event
  float weight = 1;
  float k = 10000;
  if(type==0 || type==1 || type==2){
    weight = 1;
  }else if(type==3 || type==5 || type==6){
    k = 18.7*exp(-2.4*d_pt2)+1.24;
    weight = 1+1/k;
  }else if(type==4){
    k = 18.7*exp(-2.4*d_pt1)+1.24;
    weight = 1+1/k;
  }
  if(pair_cut()==true){
    hmasspt_weight[type][icent]->Fill(d_mass, d_ptpair, weight);
    hmasspt_weight[type][10]->Fill(d_mass, d_ptpair, weight);
  }


  //d_ntpair->Fill();

}

//____________________________________________
bool ana_sgl::PairTrackcut(vector<etrk>::iterator e1, vector<etrk>::iterator e2){
  /*
  double p1 = sqrt(pow(e1.px,2)+
		   pow(e1.py,2)+
		   pow(e1.pz,2));
		   
  double p2 = sqrt(pow(e2.px,2)+
		   pow(e2.py,2)+
		   pow(e2.pz,2));
		   
  if( 
     (e1.e/p1>0.7 && e1.e/p1<1.3) ||
     (e2.e/p1>0.7 && e2.e/p2<1.3)
     ){
    return true;
  }else{
    return false;
  }
o  */
  return true;

}

//____________________________________________
void ana_sgl::check_conversion_pairs(vector<etrk> &e1, vector<etrk> &e2){
  vector<etrk>::iterator iep;
  vector<etrk>::iterator iem;
  bool reject = false;
  if(e1.size()>0 && e2.size()>0){
    for(iep = e1.begin(); iep != e1.end(); ++iep){
      reject = false;
      for(iem = e2.begin(); iem != e2.end(); ++iem){
	double mass, phiv, px, py, pz, pt, e, phi, eta, cos, psi;
	calc_vars(iep, iem, mass, phiv, px, py, pz, pt, e, phi, eta, cos, psi);
	if(magnetic_field_mm==true){ //pp
	  if(phiv<0.6 && iep->phi-iem->phi<0){ // this depends on the magntic field 
	    reject = true;
	    iem->conv_flag = 0;
	  }
	}else{
	  if(phiv>acos(-1.0)-0.6 && iep->phi-iem->phi>0){ // this depends on the magntic field 
	    reject = true;
	    iem->conv_flag = 0;
	  }
	}
      }
      if(reject==true) iep->conv_flag=0;
    }
  }
}

//____________________________________________                                                                                                                                                                                                                                                                 
void ana_sgl::check_ghost_pairs(vector<etrk> &e1){
  vector<etrk>::iterator iep;
  vector<etrk>::iterator iem;
  bool reject = false;
  if(e1.size()>1){
    for(iep = e1.begin(); iep != e1.end(); ++iep){
      reject = false;
      vector<etrk>::iterator iep2=iep;
      ++iep2;
      for(iem = iep2; iem != e1.end(); ++iem){
        double mass, phiv, px, py, pz, pt, e, phi, eta, cos, psi;
        calc_vars(iep, iem, mass, phiv, px, py, pz, pt, e, phi, eta, cos, psi);
        if(mass<0.01){
          reject = true;
          iem->ghost_flag = 0;
        }
      }
      if(reject==true) iep->ghost_flag=0;
    }
  }
}

//____________________________________________
//void ana_sgl::calc_vars(etrk iep, etrk iem, double &mass, double &phiv, double &px, double &py, double&pz,
void ana_sgl::calc_vars(vector<etrk>::iterator iep, vector<etrk>::iterator iem, double &mass, double &phiv, double &px, double &py, double&pz,
			double &pt, double &e, double &phi, double &eta, double &cos, double &psi){
  
  px = iep->px+iem->px;
  py = iep->py+iem->py;
  pz = iep->pz+iem->pz;
  pt = sqrt(px*px+py*py);
  double d_ppair = sqrt(pt*pt+pz*pz);
  static const double me=0.0005109989;
  e = sqrt(me*me+iep->px*iep->px+iep->py*iep->py+iep->pz*iep->pz)
    + sqrt(me*me+iem->px*iem->px+iem->py*iem->py+iem->pz*iem->pz);
  
  mass =  e*e-px*px-py*py-pz*pz;
  if(mass<0){
    mass = mass;
  }else{
    mass = sqrt(mass);
  }
   
  
  
  phi = atan2(py, px);
  eta = -0.5*TMath::Log((d_ppair+pz)/(d_ppair-pz));
  double p1 = sqrt(pow(iep->px,2)+pow(iep->py,2)+pow(iep->pz,2));
  double p2 = sqrt(pow(iem->px,2)+pow(iem->py,2)+pow(iem->pz,2));
  cos = acos((iep->px*iem->px+iep->py*iem->py+iep->pz*iem->pz)/(p1*p2));


  double dtheta = iep->theta-iem->theta;
  psi = asin(dtheta/cos);


  //unit vector of (pep+pem) 
  float pl = d_ppair;
  float ux = px/pl;
  float uy = py/pl;
  float uz = pz/pl;
  float ax = uy/sqrt(ux*ux+uy*uy);
  float ay = -ux/sqrt(ux*ux+uy*uy); 
  
  //momentum of e+ and e- in (ax,ay,az) axis. Note that az=0 by 
  //definition. 
  float ptep = iep->px*ax + iep->py*ay; 
  float ptem = iem->px*ax + iem->py*ay; 
  
  float pxep = iep->px;
  float pyep = iep->py;
  float pzep = iep->pz;
  float pxem = iem->px;
  float pyem = iem->py;
  float pzem = iem->pz;
  
  
  //vector product of pep X pem 
  float vpx = pyep*pzem - pzep*pyem; 
  float vpy = pzep*pxem - pxep*pzem; 
  float vpz = pxep*pyem - pyep*pxem; 
  float vp = sqrt(vpx*vpx+vpy*vpy+vpz*vpz); 
  float thev = acos(vpz/vp); 
  
  //unit vector of pep X pem 
  float vx = vpx/vp; 
  float vy = vpy/vp; 
  float vz = vpz/vp; 
  
  //The third axis defined by vector product (ux,uy,uz)X(vx,vy,vz) 
  float wx = uy*vz - uz*vy; 
  float wy = uz*vx - ux*vz; 
  float wz = ux*vy - uy*vx; 
  float wl = sqrt(wx*wx+wy*wy+wz*wz); 
  // by construction, (wx,wy,wz) must be a unit vector. 
  if(fabs(wl - 1.0) > 0.00001) cout << "Calculation error in W vector"<<endl; 
  // measure angle between (wx,wy,wz) and (ax,ay,0). The angle between them 
  // should be small if the pair is conversion 
  //
  float cosPhiV = wx*ax + wy*ay; 
  phiv = acos(cosPhiV); 
  
}

//____________________________________________
bool ana_sgl::pair_cut(void){
  bool ret = true;


  if(d_flag_phiv==true){
    if(d_run>169591 && d_phiv>acos(-1.0)-d_phiv_cut){ 
      ret = false;
    }
    if(d_run<169591 && d_phiv<d_phiv_cut){ 
      ret = false;
    }
  }


  if(d_flag_emc_cut==true){
    if( !( 
	  (d_e1/d_pt1>d_emc_low && d_e1/d_pt1<d_emc_high) ||
	  (d_e2/d_pt2>d_emc_low && d_e2/d_pt2<d_emc_high)
	  )
	){
      ret = false;
    }
  }


  if(d_flag_pt_cut==true){
    if(!(d_pt1>d_pt_cut_low && d_pt2>d_pt_cut_low && 
	 d_pt1<d_pt_cut_high && d_pt2<d_pt_cut_high)
       ){
      ret = false;
    }
  }

  return ret;

}

//____________________________________________ 
void ana_sgl::print_cuts(void){

  cout<<" ********** lists of cuts ********** "<<endl;
  cout<<d_tpc_dedx_low<<" < TPC dE/dx < "<<d_tpc_dedx_high<<endl;


  if(d_flag_kaon_veto==true){
    cout<<d_dedx_kaon_veto_low<<" < TPC veto for Kaon < "<<d_dedx_kaon_veto_high<<endl;
  }else{
    cout<<"  No TPC Kaon Veto "<<endl;
  }

  if(d_flag_proton_veto==true){
    cout<<d_dedx_proton_veto_low<<" < TPC veto for proton < "<<d_dedx_proton_veto_high<<endl;
  }else{
    cout<<" No TPC Proton Veto "<<endl;
  }

  if(d_flag_tof_cut==true){
    cout<<d_tof_low<<" < kTOFnSigmaEle < "<<d_tof_high<<endl;
  }else{
    cout<<" No TOF cuts "<<endl;
  }
  if(d_flag_emc_cut==true){
    cout<<d_emc_low<<" < EMC E/p (or for pairs) < "<<d_emc_high<<endl;
  }else{
    cout<<" No EMC cuts "<<endl;
  }
  if(d_flag_phiv==true){
    cout<<" Phiv cuts "<<d_phiv_cut<<endl;
  }else{
    cout<<" No phiv cuts "<<endl;
  }

  if(d_flag_pt_cut==true){
    cout<<d_pt_cut_low<<" < pT of electrons < "<<d_pt_cut_high<<endl;
  }else{
    cout<<" No electron pt cuts "<<endl;
  }

  if(d_conv_flag==true){
    cout<<" remove converson pairs from electron/positron lists "<<endl;
  }else{
    cout<<" conversion candidates are in the pool."<<endl;
  }
  cout<<" ********** end of lists of cuts ********** "<<endl;

  if(sel_trigger==0){
    cout<<" MB "<<endl;
  }else if(sel_trigger==1){
    cout<<" SemiCentral "<<endl;
  }else if(sel_trigger==2){
    cout<<" Central "<<endl;
  }

}
  
  
