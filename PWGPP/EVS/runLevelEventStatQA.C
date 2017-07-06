#ifndef __CINT__
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TPad.h"
#endif
#include "triggerInfo.C"

// TODO read number of bits from AliVEvent?
#define NBITS 29
#define NMAXCLASSES 100
TString bitNames[NBITS] = {
    "kINT1",
    "kINT7",
    "kMUON",
    "kHighMultSPD",
    "kEMC1",
    "kINT5",
    "kINT7inMUON",
    "kMuonSingleHighPt7",
    "kMuonLikeLowPt7",
    "kMuonUnlikeLowPt7",
    "kEMC7",
    "kMuonSingleLowPt7",
    "kPHI1",
    "kPHI78",
    "kEMCEJE",
    "kEMCEGA",
    "kHighMultV0",
    "kSemiCentral",
    "kDG",
    "kZED",
    "kSPI78",
    "kINT8",
    "kMuonSingleLowPt8",
    "kMuonSingleHighPt8",
    "kMuonLikeLowPt8",
    "kMuonUnlikeLowPt8",
    "kMuonUnlikeLowPt0",
    "kUserDefined",
    "kTRD"
};

void writeTree(TFile* fout, TTree* t){
  fout->cd();
  t->Fill();
  t->Write();
  fout->Close();
}

//Int_t runLevelEventStatQA(TString qafilename="event_stat.root", Int_t run=254422, TString ocdbStorage = "raw://"){
//Int_t runLevelEventStatQA(TString qafilename="event_stat.root", Int_t run=255042, TString ocdbStorage = "raw://"){
//Int_t runLevelEventStatQA(TString qafilename="EventStat_temp.root", Int_t run=256676, TString ocdbStorage = "local:///cvmfs/alice.cern.ch/calibration/data/2016/OCDB"){
Int_t runLevelEventStatQA(TString qafilename="event_stat_new.root", Int_t run=272399, TString ocdbStorage = "local:///cvmfs/alice.cern.ch/calibration/data/2017/OCDB"){
  
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetLegendBorderSize(0);
  // tree variables
  TObjArray classes = TObjArray();
  TObjString activeDetectors = TObjString();
  TObjString objPartition = TObjString();
  TObjString objLhcState = TObjString();
  TObjString objLhcPeriod = TObjString();
  Int_t fill               = 0;
  Double_t run_duration    = 0;
  Int_t nBCsPerOrbit       = 0;
  Double_t refCounts       = 0;
  Double_t mu              = 0;
  Double_t lumi_seen       = 0;
  Double_t interactionRate = 0;
  ULong64_t class_lMb[NMAXCLASSES]         = {0};
  ULong64_t class_lMa[NMAXCLASSES]         = {0};
  ULong64_t class_l0b[NMAXCLASSES]         = {0};
  ULong64_t class_l0a[NMAXCLASSES]         = {0};
  ULong64_t class_l1b[NMAXCLASSES]         = {0};
  ULong64_t class_l1a[NMAXCLASSES]         = {0};
  ULong64_t class_l2b[NMAXCLASSES]         = {0};
  ULong64_t class_l2a[NMAXCLASSES]         = {0};
  ULong64_t class_duration[NMAXCLASSES]    = {0};
  Double_t  class_lifetime[NMAXCLASSES]    = {0};
  Double_t  class_lumi[NMAXCLASSES]        = {0};
  Double_t  class_ds[NMAXCLASSES]          = {0};
  ULong64_t alias_recorded[NBITS]          = {0};
  ULong64_t alias_reconstructed[NBITS]     = {0};
  ULong64_t alias_accepted[NBITS]          = {0};
  ULong64_t alias_acc_step1[NBITS]         = {0};
  ULong64_t alias_acc_step2[NBITS]         = {0};
  ULong64_t alias_acc_step3[NBITS]         = {0};
  ULong64_t alias_acc_step4[NBITS]         = {0};
  ULong64_t alias_acc_step5[NBITS]         = {0};
  ULong64_t alias_acc_step6[NBITS]         = {0};
  ULong64_t alias_acc_step7[NBITS]         = {0};
  ULong64_t alias_acc_step8[NBITS]         = {0};
  ULong64_t alias_acc_step9[NBITS]         = {0};
  Double_t alias_l0b_rate[NBITS]           = {0};
  Double_t alias_lifetime[NBITS]           = {0};
  Double_t alias_lumi_recorded[NBITS]      = {0};
  Double_t alias_lumi_reconstructed[NBITS] = {0};
  Double_t alias_lumi_accepted[NBITS]      = {0};
  Int_t timeStart = 0;
  Int_t timeEnd = 0;
  Double_t meanV0MOn = 0;
  Double_t meanV0MOf = 0;
  Double_t meanOFO = 0;
  Double_t meanTKL = 0;
  Double_t meanErrV0MOn = 0;
  Double_t meanErrV0MOf = 0;
  Double_t meanErrOFO = 0;
  Double_t meanErrTKL = 0;
  Double_t meanV0MOnHM = 0;
  Double_t meanV0MOfHM = 0;
  Double_t meanOFOHM = 0;
  Double_t meanTKLHM = 0;
  Double_t meanErrV0MOnHM = 0;
  Double_t meanErrV0MOfHM = 0;
  Double_t meanErrOFOHM = 0;
  Double_t meanErrTKLHM = 0;
  Double_t thresholdV0M = 0;
  Double_t aV0MOnVsOfVal = 0;
  Double_t bV0MOnVsOfVal = 0;
  Double_t aV0MOnVsOfErr = 0;
  Double_t bV0MOnVsOfErr = 0;
  Double_t aSPDOnVsOfVal = 0;
  Double_t bSPDOnVsOfVal = 0;
  Double_t aSPDOnVsOfErr = 0;
  Double_t bSPDOnVsOfErr = 0;
  Double_t aV0MOnVsOfOADB = 0;
  Double_t bV0MOnVsOfOADB = 0;
  Double_t aSPDOnVsOfOADB = 0;
  Double_t bSPDOnVsOfOADB = 0;

  TH2F* hHistStat = new TH2F();

  TFile* fout = new TFile("trending.root","recreate");
  TTree* t = new TTree("trending","tree of trending variables");
  t->Branch("run",&run);
  t->Branch("fill",&fill);
  t->Branch("bcs",&nBCsPerOrbit);
  t->Branch("run_duration",&run_duration);
  t->Branch("mu",&mu);
  t->Branch("interactionRate",&interactionRate);
  t->Branch("refCounts",&refCounts);
  t->Branch("lumi_seen",&lumi_seen);
  t->Branch("classes",&classes);
  t->Branch("class_lMb",&class_lMb,Form("class_lMb[%i]/l",NMAXCLASSES));
  t->Branch("class_lMa",&class_lMa,Form("class_lMa[%i]/l",NMAXCLASSES));
  t->Branch("class_l0b",&class_l0b,Form("class_l0b[%i]/l",NMAXCLASSES));
  t->Branch("class_l0a",&class_l0a,Form("class_l0a[%i]/l",NMAXCLASSES));
  t->Branch("class_l1b",&class_l1b,Form("class_l1b[%i]/l",NMAXCLASSES));
  t->Branch("class_l1a",&class_l1a,Form("class_l1a[%i]/l",NMAXCLASSES));
  t->Branch("class_l2b",&class_l2b,Form("class_l2b[%i]/l",NMAXCLASSES));
  t->Branch("class_l2a",&class_l2a,Form("class_l2a[%i]/l",NMAXCLASSES));
  t->Branch("class_duration",&class_duration,Form("class_duration[%i]/l",NMAXCLASSES));
  t->Branch("class_lifetime",&class_lifetime,Form("class_lifetime[%i]/D",NMAXCLASSES));
  t->Branch("class_lumi",&class_lumi,Form("class_lumi[%i]/D",NMAXCLASSES));
  t->Branch("class_ds",&class_ds,Form("class_ds[%i]/D",NMAXCLASSES));
  t->Branch("alias_recorded",&alias_recorded,Form("alias_recorded[%i]/l",NBITS));
  t->Branch("alias_reconstructed",&alias_reconstructed,Form("alias_reconstructed[%i]/l",NBITS));
  t->Branch("alias_accepted",&alias_accepted,Form("alias_accepted[%i]/l",NBITS));
  t->Branch("alias_acc_step1",&alias_acc_step1,Form("alias_acc_step1[%i]/l",NBITS));
  t->Branch("alias_acc_step2",&alias_acc_step2,Form("alias_acc_step2[%i]/l",NBITS));
  t->Branch("alias_acc_step3",&alias_acc_step3,Form("alias_acc_step3[%i]/l",NBITS));
  t->Branch("alias_acc_step4",&alias_acc_step4,Form("alias_acc_step4[%i]/l",NBITS));
  t->Branch("alias_acc_step5",&alias_acc_step5,Form("alias_acc_step5[%i]/l",NBITS));
  t->Branch("alias_acc_step6",&alias_acc_step6,Form("alias_acc_step6[%i]/l",NBITS));
  t->Branch("alias_acc_step7",&alias_acc_step7,Form("alias_acc_step7[%i]/l",NBITS));
  t->Branch("alias_acc_step8",&alias_acc_step6,Form("alias_acc_step8[%i]/l",NBITS));
  t->Branch("alias_acc_step9",&alias_acc_step7,Form("alias_acc_step9[%i]/l",NBITS));
  t->Branch("alias_l0b_rate",&alias_lifetime,Form("alias_l0b_rate[%i]/D",NBITS));
  t->Branch("alias_lifetime",&alias_lifetime,Form("alias_lifetime[%i]/D",NBITS));
  t->Branch("alias_lumi_recorded",&alias_lumi_recorded,Form("alias_lumi_recorded[%i]/D",NBITS));
  t->Branch("alias_lumi_reconstructed",&alias_lumi_reconstructed,Form("alias_lumi_reconstructed[%i]/D",NBITS));
  t->Branch("alias_lumi_accepted",&alias_lumi_accepted,Form("alias_lumi_accepted[%i]/D",NBITS));
  t->Branch("activeDetectors",&activeDetectors);
  t->Branch("lhcState",&objLhcState);
  t->Branch("lhcPeriod",&objLhcPeriod);
  t->Branch("partition",&objPartition);
  t->Branch("timeStart",&timeStart);
  t->Branch("timeEnd",&timeEnd);
  t->Branch("meanV0MOn",&meanV0MOn);
  t->Branch("meanV0MOf",&meanV0MOf);
  t->Branch("meanOFO",&meanOFO);
  t->Branch("meanTKL",&meanTKL);
  t->Branch("meanErrV0MOn",&meanErrV0MOn);
  t->Branch("meanErrV0MOf",&meanErrV0MOf);
  t->Branch("meanErrOFO",&meanErrOFO);
  t->Branch("meanErrTKL",&meanErrTKL);
  t->Branch("meanV0MOnHM",&meanV0MOnHM);
  t->Branch("meanV0MOfHM",&meanV0MOfHM);
  t->Branch("meanOFOHM",&meanOFOHM);
  t->Branch("meanTKLHM",&meanTKLHM);
  t->Branch("meanErrV0MOnHM",&meanErrV0MOnHM);
  t->Branch("meanErrV0MOfHM",&meanErrV0MOfHM);
  t->Branch("meanErrOFOHM",&meanErrOFOHM);
  t->Branch("meanErrTKLHM",&meanErrTKLHM);
  t->Branch("thresholdV0M",&thresholdV0M);
  t->Branch("aV0MOnVsOfVal",&aV0MOnVsOfVal);
  t->Branch("bV0MOnVsOfVal",&bV0MOnVsOfVal);
  t->Branch("aV0MOnVsOfErr",&aV0MOnVsOfErr);
  t->Branch("bV0MOnVsOfErr",&bV0MOnVsOfErr);
  t->Branch("aSPDOnVsOfVal",&aSPDOnVsOfVal);
  t->Branch("bSPDOnVsOfVal",&bSPDOnVsOfVal);
  t->Branch("aSPDOnVsOfErr",&aSPDOnVsOfErr);
  t->Branch("bSPDOnVsOfErr",&bSPDOnVsOfErr);
  t->Branch("aV0MOnVsOfOADB",&aV0MOnVsOfOADB);
  t->Branch("bV0MOnVsOfOADB",&bV0MOnVsOfOADB);
  t->Branch("aSPDOnVsOfOADB",&aSPDOnVsOfOADB);
  t->Branch("bSPDOnVsOfOADB",&bSPDOnVsOfOADB);
  
  t->Branch("hHistStat",&hHistStat);
  
  TString refClass="";
  Double_t refSigma=-1;
  Double_t refEff = 1.;
  Double_t refMu = 1.e-11;
  if      (               run<=118501) { refSigma=  62.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=118502 && run<=118561) { refSigma=  47.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_0.90: 47mb=52 mb *0.91=sigma(INEL)*R(INT1/INEL) (arxiv: 1208.4968, fig.10 + table 3)
  else if (run>=118903 && run<=120829) { refSigma=  62.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=121039 && run<=121040) { refSigma=  47.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_0.90: 47mb=52 mb *0.91=sigma(INEL)*R(INT1/INEL) (arxiv: 1208.4968, fig.10 + table 3)
  else if (run>=121041 && run<=126437) { refSigma=  62.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=126438 && run<=127718) { refSigma=  62.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=127719 && run<=127730) { refSigma=  62.; refEff = 1.00; refClass = "CINT1B-ABCE-NOPF-ALL";   } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=127731 && run<=136848) { refSigma=  62.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_7.00: 62mb=54.3mb*1.15=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=136849 && run<=139316) { refSigma=5970.; refEff = 0.78; refClass = "C0SMH-B-NOPF-ALL";       } // PbPb_2.76: (Oyama,2011-05-20,RunCond), sigma_hardronic = 7.64 b
  else if (run>=139328 && run<=139517) { refSigma=5970.; refEff = 0.78; refClass = "C0SMH-B-NOPF-ALLNOTRD";  } // PbPb_2.76: (Oyama,2011-05-20,RunCond), sigma_hardronic = 7.64 b
  else if (run>=145289 && run<=146860) { refSigma=  57.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146808 && run<=146814) { refSigma=  57.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALL";       } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=145815 && run<=146856) { refSigma=  57.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146857 && run<=146857) { refSigma=  57.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALL";       } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=146858 && run<=146860) { refSigma=  57.; refEff = 1.00; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_2.76: 57mb=47.7mb*1.20=sigma(VBAND)*R(INT1/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=148370 && run<=157078) { refSigma=  54.; refEff = 1.00; refClass = "CVBAND-B-NOPF-ALLNOTRD"; } // pp_7.00: 54.3mb (Martino,2012-03-12,RunCond)
  else if (run>=157079 && run<=165746) { refSigma=  24.; refEff = 0.44; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_7.00: 24mb=54.3mb*0.44=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=166477 && run<=170593) { refSigma=4100.; refEff = 0.54; refClass = "CVLN-B-NOPF-ALLNOTRD";   } // PbPb_2.76: (Martino,2013-03-15,RunCond)
  else if (run>=176658 && run<=177143) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177146 && run<=177147) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177148 && run<=177149) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177150 && run<=177165) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177166 && run<=177166) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177167 && run<=177167) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177168 && run<=177168) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177169 && run<=177172) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177173 && run<=177173) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177174 && run<=177506) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=177507 && run<=178017) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=178018 && run<=178029) { refSigma=  67.; refEff = 0.90; refClass = "CINT1-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), CINT1/C0TVX=2.7 from 178052
  else if (run>=178030 && run<=178053) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=178055 && run<=178062) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALL";       } // pp_8.00: (Artem, 2013-10-04,RunCond), vdM
  else if (run>=178062 && run<=178220) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=179444 && run<=180715) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-S-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=180716 && run<=180720) { refSigma=  56.; refEff = 0.75; refClass = "CINT7-S-NOPF-ALLNOTRD";  } // no C0TVX in these runs, taking VBAND cross section
  else if (run>=180721 && run<=184844) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-S-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=184845 && run<=184990) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=184991 && run<=188229) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-S-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=188230 && run<=188366) { refSigma=1590.; refEff = 0.76; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pPb_5.02: pilot. arxiv:1405.1849
  else if (run>=188367 && run<=193692) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-S-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=193693 && run<=193766) { refSigma=  25.; refEff = 0.33; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_8.00: (Artem, 2013-10-04,RunCond), TOTEM INEL = 74.7+/-1.7 mb
  else if (run>=195344 && run<=197388) { refSigma=1590.; refEff = 0.76; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pPb_5.02: arxiv:1405.1849
  else if (run>=197470 && run<=197692) { refSigma=  18.; refEff = 0.39; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // pp_2.76: 18mb=47.7mb*0.39=sigma(VBAND)*R(0TVX/VBAND) (Martino,2012-03-12,RunCond)
  else if (run>=221835 && run<=223669) { refSigma= 16.8; refEff = 0.32; refClass = "CADAND-B-NOPF-ALLNOTRD"; } // estimates from Martino
  else if (run>=221670 && run<=223983) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // estimates from Martino and MC
  else if (run>=223984 && run<=223984) { refSigma= 50.0; refEff = 0.66; refClass = "CADAND-B-NOPF-ALLNOTRD"; } // estimates from Martino and MC
  else if (run>=223985 && run<=226110) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // estimates from Martino and MC
  else if (run>=226111 && run<=226115) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=226116 && run<=228909) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // estimates from Martino and MC
  else if (run>=228910 && run<=229376) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=229386 && run<=229398) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-MUON";      } // estimates from Martino and MC
  else if (run>=229409 && run<=229410) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=229416 && run<=229893) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-MUON";      } // estimates from Martino and MC
  else if (run>=229894 && run<=229899) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // estimates from Martino and MC
  else if (run>=229942 && run<=231321) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-MUON";      } // estimates from Martino and MC
  else if (run>=232914 && run<=233858) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENT";      } // estimates from Martino and MC
  else if (run>=233910 && run<=234050) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-ALLNOTRD";  } // estimates from Martino and MC
  else if (run>=234051 && run<=238669) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENT";      } // estimates from Martino and MC
  else if (run>=238670 && run<=240150) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=240151 && run<=240151) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-MUON";      } // estimates from Martino and MC
  else if (run>=240152 && run<=243373) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=243374 && run<=243398) { refSigma= 21.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=243399 && run<=243984) { refSigma=6700.; refEff = 0.90; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=243985 && run<=244912) { refSigma= 21.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=244913 && run<=246994) { refSigma=4600.; refEff = 0.60; refClass = "C0V0M-B-NOPF-CENTNOTRD"; } // estimates from Cvetan and Alberica
  else if (run>=246995 && run<=256147) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=256148 && run<=256157) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENT";      } // estimates from Martino and MC
  else if (run>=256158 && run<=265301) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Martino and MC
  else if (run>=265304 && run<=265586) { refSigma=1590.; refEff = 0.76; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // taken from old p-Pb
  else if (run>=265587 && run<=267131) { refSigma=1715.; refEff = 0.82; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // T0 efficiency estimate from MC, refSigma assuming cs(INEL) = 2092
  else if (run>=267132 && run<=267166) { refSigma=1590.; refEff = 0.76; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // taken from old p-Pb
  else if (run>=267167               ) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // back to pp at 13TeV

  
  Double_t orbitRate = 11245.;
  TString partition;
  TString lhcState;
  TString lhcPeriod;
  TString activeDetectorsString;
  //Int_t run, TString ocdbStorage, TString &partition, TString &activeDetectorsString, Double_t& run_duration, 
  classes = GetClasses(run,ocdbStorage,partition,activeDetectorsString,run_duration,class_lMb,class_lMa,class_l0b,class_l0a,class_l1b,class_l1a,class_l2b,class_l2a,class_duration);
  activeDetectors.SetString(activeDetectorsString.Data());
  objPartition.SetString(partition.Data());
  Int_t status = triggerInfo(run,ocdbStorage,lhcPeriod,lhcState,fill,nBCsPerOrbit,timeStart,timeEnd);
  if (status>0) { writeTree(fout,t); return 1; }
  objLhcState.SetString(lhcState.Data());
  objLhcPeriod.SetString(lhcPeriod.Data());

  AliTriggerClass* refClassObject = (AliTriggerClass*) classes.FindObject(refClass);
  if (!refClassObject) { writeTree(fout,t); return 2; }
  Int_t refId = classes.IndexOf(refClassObject);
  TString refCluster = refClassObject->GetCluster()->GetName();
  refCounts = (activeDetectorsString.Contains("TRD") && (refCluster.EqualTo("CENT") || refCluster.EqualTo("ALL") || refCluster.EqualTo("FAST"))) ? class_lMb[refId] : class_l0b[refId];
  if (refClass.Contains("1B-ABCE-")){
    Int_t emptyClassId = classes.IndexOf(classes.FindObject("CBEAMB-ABCE-NOPF-ALL"));
    Int_t emptyL0B = class_l0b[emptyClassId];
    if (nBCsPerOrbit<0 && emptyL0B>0) nBCsPerOrbit = Double_t(emptyL0B)/orbitRate/run_duration;
  } else {
    nBCsPerOrbit= refClassObject->GetBCMask()->GetNUnmaskedBCs();
  }
  
  Double_t totalBCs = orbitRate*run_duration*nBCsPerOrbit;
  if (totalBCs<1 || refCounts<1) { writeTree(fout,t); return 3; }
  refMu = -TMath::Log(1-Double_t(refCounts)/totalBCs); // mu
  Double_t refInteractionRate = (run_duration>1e-10) ? (refMu>1.e-10 ? refMu/(1-TMath::Exp(-refMu)):1)*refCounts/run_duration : 0;
  interactionRate = refInteractionRate/refEff; // INEL interaction rate 
  mu              = refMu/refEff;              // INEL mu value
  if (refSigma>1.e-10) lumi_seen = run_duration*refInteractionRate/refSigma/1000; //[ub-1]

  for (Int_t i=0;i<classes.GetEntriesFast();i++){
    // printf("%30s %12lli %10lli %10lli %10lli %10lli %10lli\n",classes.At(i)->GetName(),class_l0b[i],class_l0a[i],class_l1b[i],class_l1a[i],class_l2b[i],class_l2a[i]);
    class_lifetime[i] = class_lMb[i]>0 ? Double_t(class_lMa[i])/class_lMb[i]: 0;
    class_lifetime[i]*= class_l0b[i]>0 ? Double_t(class_l0a[i])/class_l0b[i]: 0;
    class_lifetime[i]*= class_l1b[i]>0 ? Double_t(class_l1a[i])/class_l1b[i]: 0;
    class_lifetime[i]*= class_l2b[i]>0 ? Double_t(class_l2a[i])/class_l2b[i]: 0;
    class_lumi[i] = lumi_seen*class_lifetime[i];
    AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
    cl->GetDownscaleFactor(class_ds[i]);
  }

  // special treatment for Pb-Pb
  if (run>=244917 && run<=246994) {
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
      TObjArray* tokens = TString(cl->GetName()).Tokenize("-");
      TString cluster = tokens->At(3)->GetName();
      TString lifetimeClassName = Form("C0VHM-B-NOPF-%s",cluster.Data());
      if (run<245256) lifetimeClassName = Form("C0V0M-B-NOPF-%s",cluster.Data());
      TObject* lifetimeClass = classes.FindObject(lifetimeClassName.Data());
      if (lifetimeClass) {
        Int_t index = classes.IndexOf(lifetimeClass);
        if (class_ds[index]>0) {
          printf("%s %f %f\n",lifetimeClassName.Data(),class_ds[i],class_ds[index]);
          Float_t ds_ratio = class_ds[i]/class_ds[index];
          class_lifetime[i]=class_lifetime[index]*ds_ratio;
          class_lumi[i]=class_lumi[index]*ds_ratio;
        }
      }
      tokens->Delete();
      delete tokens;
    }
  }

  // special treatment for INT11
  if (run>=252235 && run<=255011 && run!=253434 && run!=253437) {
    AliTriggerClass* clINT11 = (AliTriggerClass*) classes.FindObject("CINT11-B-NOPF-CENTNOTRD");
    AliTriggerClass* cl0TVX  = (AliTriggerClass*) classes.FindObject("C0TVX-B-NOPF-CENTNOTRD");
    if (clINT11 && cl0TVX) {
      Int_t iCINT11 = classes.IndexOf(clINT11);
      Int_t iC0TVX  = classes.IndexOf(cl0TVX);
      class_lifetime[iCINT11] = class_lifetime[iC0TVX];
      class_lumi[iCINT11]     = class_lumi[iC0TVX];
    }
  }
  
  // special treatment for CUP2
  if (run>=255276 && run<=256418) {
    Int_t iCUP2=-1;
    Int_t iCTRUE_NOPF=-1;
    Int_t iCTRUE_SPD1=-1;
    Int_t iCVHMV0M=-1;
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
      if (TString(cl->GetName()).EqualTo("CCUP2-B-SPD1-CENTNOTRD")   ) iCUP2=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-SPD1-CENTNOTRD")   ) iCTRUE_SPD1=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-NOPF-CENTNOTRD")   ) iCTRUE_NOPF=i;
      if (TString(cl->GetName()).EqualTo("CVHMV0M-B-NOPF-CENTNOTRD") ) iCVHMV0M=i;
    }
    Int_t iREF = iCVHMV0M>=0 ? iCVHMV0M : iCTRUE_NOPF;
    if (iCUP2>=0 && iCTRUE_NOPF>=0 && iCTRUE_SPD1>=0 && iREF>=0) {
      class_lifetime[iCUP2] = (class_ds[iREF]>1e-7 && class_l0b[iREF]>0) ? (Double_t(class_l0a[iREF])/class_l0b[iREF]/class_ds[iREF]) : (Double_t(class_l0a[iCTRUE_NOPF])/class_l0b[iCTRUE_NOPF]/class_ds[iCTRUE_NOPF]);
      class_lifetime[iCUP2]*=class_l2a[iCTRUE_SPD1];
      class_lifetime[iCUP2]*=class_l2a[iCTRUE_NOPF]>0 ? 1./class_l2a[iCTRUE_NOPF] : 0;
      class_lifetime[iCUP2]*=class_ds[iCUP2];
      class_lumi[iCUP2] = lumi_seen*class_lifetime[iCUP2];
    }
    if (iCUP2>=0) printf("lumi=%f\n",class_lumi[iCUP2]);
  }
  
  // special treatment for CUP13
  if (run>=256504) {
    Int_t iCUP13=-1;
    Int_t iCTRUE_NOPF=-1;
    Int_t iCTRUE_SPD1=-1;
    Int_t iCEMC7=-1;
    Int_t iCPHI7=-1;
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
      if (TString(cl->GetName()).EqualTo("CCUP13-B-SPD1-CENTNOTRD")  ) iCUP13=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-SPD1-CENTNOTRD")   ) iCTRUE_SPD1=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-NOPF-CENTNOTRD")   ) iCTRUE_NOPF=i;
      if (TString(cl->GetName()).EqualTo("CEMC7EG1-B-NOPF-CENTNOTRD")) iCEMC7=i;
      if (TString(cl->GetName()).EqualTo("CPHI7-B-NOPF-CENTNOTRD"))    iCPHI7=i;
    }
    Int_t iREF = iCEMC7>=0 ? iCEMC7 : (iCPHI7>=0 ? iCPHI7 : iCTRUE_NOPF);
    if (iCUP13>=0 && iCTRUE_NOPF>=0 && iCTRUE_SPD1>=0 && iREF>=0) {
      class_lifetime[iCUP13] = (class_ds[iREF]>1e-7 && class_l0b[iREF]>0) ? (Double_t(class_l0a[iREF])/class_l0b[iREF]/class_ds[iREF]) : (Double_t(class_l0a[iCTRUE_NOPF])/class_l0b[iCTRUE_NOPF]/class_ds[iCTRUE_NOPF]);
      class_lifetime[iCUP13]*=class_l2a[iCTRUE_SPD1]>0 ? class_l2a[iCTRUE_SPD1] : 1;
      class_lifetime[iCUP13]*=class_l2a[iCTRUE_NOPF]>0 ? 1./class_l2a[iCTRUE_NOPF] : 0;
      class_lifetime[iCUP13]*=1-TMath::Power(1-class_ds[iCUP13],4.);
      class_lumi[iCUP13] = lumi_seen*class_lifetime[iCUP13];
    }
  }
  
  // special treatment for CUP25
  if (run>=272151) {
    Int_t iCUP25=-1;
    Int_t iCTRUE_NOPF=-1;
    Int_t iCTRUE_SPD1=-1;
    Int_t iCEMC7=-1;
    Int_t iCPHI7=-1;
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
      if (TString(cl->GetName()).EqualTo("CCUP25-B-SPD1-CENTNOTRD")  ) iCUP25=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-SPD1-CENTNOTRD")   ) iCTRUE_SPD1=i;
      if (TString(cl->GetName()).EqualTo("CTRUE-B-NOPF-CENTNOTRD")   ) iCTRUE_NOPF=i;
      if (TString(cl->GetName()).EqualTo("CEMC7EG1-B-NOPF-CENTNOTRD")) iCEMC7=i;
      if (TString(cl->GetName()).EqualTo("CPHI7-B-NOPF-CENTNOTRD"))    iCPHI7=i;
    }
    Int_t iREF = iCEMC7>=0 ? iCEMC7 : (iCPHI7>=0 ? iCPHI7 : iCTRUE_NOPF);
    if (iCUP25>=0 && iCTRUE_NOPF>=0 && iCTRUE_SPD1>=0 && iREF>=0) {
      class_lifetime[iCUP25] = (class_ds[iREF]>1e-7 && class_l0b[iREF]>0) ? (Double_t(class_l0a[iREF])/class_l0b[iREF]/class_ds[iREF]) : (Double_t(class_l0a[iCTRUE_NOPF])/class_l0b[iCTRUE_NOPF]/class_ds[iCTRUE_NOPF]);
      class_lifetime[iCUP25]*=class_l2a[iCTRUE_SPD1]>0 ? class_l2a[iCTRUE_SPD1] : 1;
      class_lifetime[iCUP25]*=class_l2a[iCTRUE_NOPF]>0 ? 1./class_l2a[iCTRUE_NOPF] : 0;
      class_lifetime[iCUP25]*=class_ds[iCUP25];
      class_lumi[iCUP25] = lumi_seen*class_lifetime[iCUP25];
    }
  }
  
  // special treatment for UPC central barrel in p-Pb
  if (run>=265589 && run<=267166) {
    for (Int_t i=0;i<classes.GetEntriesFast();i++){
      AliTriggerClass* cl = (AliTriggerClass*) classes.At(i);
      if (!cl) continue;
      TString name = cl->GetName();
      if (!(name.Contains("CCUP") && name.Contains("-CENTNOTRD"))) continue;
      Bool_t pf = name.Contains("-SPD2-");
      TObject* emc = classes.FindObject(pf ? "CEMC7EG1-B-SPD2-CENTNOTRD" : "CEMC7EG1-B-NOPF-CENTNOTRD");
      TObject* phs = classes.FindObject(pf ? "CPHI7PHL-B-SPD2-CENTNOTRD" : "CPHI7PHL-B-NOPF-CENTNOTRD");
      if (!emc && !phs) continue;
      Int_t ilt = emc!=0 ? classes.IndexOf(emc) : classes.IndexOf(phs);
      if (ilt<0) continue;
      class_lifetime[i] = class_lifetime[ilt];
      class_lumi[i] = lumi_seen*class_lifetime[i]*class_ds[i];
    }
  }

  
  TFile* fin = new TFile(qafilename);
  if (!fin) {
    printf("qa file not found");
    writeTree(fout,t);
    return 4;
  }
  
  AliPhysicsSelection* ps = 0;
  
  if (qafilename.Contains("event_stat")) {
    hHistStat = (TH2F*) fin->Get("fHistStat");
  } else {
    TList* statsout = (TList*) fin->Get("cstatsout");
    ps = statsout ? (AliPhysicsSelection*) statsout->FindObject("AliPhysicsSelection") : 0;
    hHistStat = ps ? (TH2F*) ps->GetStatistics("") : 0;
  }
  
  if (!hHistStat) {
     printf("fHistStat not found\n");
     writeTree(fout,t);
     return 5; 
  }
  
  fout->cd();

  meanV0MOn = 0;
  meanV0MOf = 0;
  meanOFO   = 0;
  meanTKL   = 0;
  meanErrV0MOn = 0;
  meanErrV0MOf = 0;
  meanErrOFO   = 0;
  meanErrTKL   = 0;
  meanV0MOnHM = 0;
  meanV0MOfHM = 0;
  meanOFOHM   = 0;
  meanTKLHM   = 0;
  meanErrV0MOnHM = 0;
  meanErrV0MOfHM = 0;
  meanErrOFOHM   = 0;
  meanErrTKLHM   = 0;
  thresholdV0M   = 0;
  Double_t* abV0M;
  Double_t* abSPD;

  for (Int_t j=1;j<=hHistStat->GetNbinsY();j++){
    TString label = hHistStat->GetYaxis()->GetBinLabel(j);
    // kINT7
    TList* list = NULL;
    if (qafilename.Contains("event_stat")) {
      list = (TList*) fin->Get(Form("trigger_histograms_%s/histos",label.Data()));
    } else {
      AliTriggerAnalysis* ta = ps->GetTriggerAnalysis(j-1);
      list = ta->GetHistList();
    }
    if (!list) continue;
    TH1F* hV0MOnAcc = (TH1F*) list->FindObject("fHistV0MOnAcc");
    TH1F* hV0MOfAcc = (TH1F*) list->FindObject("fHistV0MOfAcc");
    TH1F* hOFOAcc   = (TH1F*) list->FindObject("fHistOFOAcc");
    TH1F* hTKLAcc   = (TH1F*) list->FindObject("fHistTKLAcc");
    TH2F* hV0MOnVsOfCln = (TH2F*) list->FindObject("fHistV0MOnVsOfCln");
    TH2F* hSPDOnVsOfCln = (TH2F*) list->FindObject("fHistSPDOnVsOfCln");
    if (label.Contains(" &2 ")) {
      meanV0MOn    = hV0MOnAcc ? hV0MOnAcc->GetMean()      : 0;
      meanV0MOf    = hV0MOfAcc ? hV0MOfAcc->GetMean()      : 0;
      meanOFO      = hOFOAcc   ? hOFOAcc->GetMean()        : 0;
      meanTKL      = hTKLAcc   ? hTKLAcc->GetMean()        : 0;
      meanErrV0MOn = hV0MOnAcc ? hV0MOnAcc->GetMeanError() : 0;
      meanErrV0MOf = hV0MOfAcc ? hV0MOfAcc->GetMeanError() : 0;
      meanErrOFO   = hOFOAcc   ? hOFOAcc->GetMeanError()   : 0;
      meanErrTKL   = hTKLAcc   ? hTKLAcc->GetMeanError()   : 0;
      aV0MOnVsOfOADB = ((TF1*) hV0MOnVsOfCln->GetListOfFunctions()->FindObject("fFuncV0MOnVsOf"))->GetParameter(0);
      bV0MOnVsOfOADB = ((TF1*) hV0MOnVsOfCln->GetListOfFunctions()->FindObject("fFuncV0MOnVsOf"))->GetParameter(1);
      aSPDOnVsOfOADB = ((TF1*) hSPDOnVsOfCln->GetListOfFunctions()->FindObject("fFuncSPDOnVsOf"))->GetParameter(0);
      bSPDOnVsOfOADB = ((TF1*) hSPDOnVsOfCln->GetListOfFunctions()->FindObject("fFuncSPDOnVsOf"))->GetParameter(1);
      
      abV0M = getCutParametersOnlineVsOffline(hV0MOnVsOfCln, "V0M");
      aV0MOnVsOfVal = abV0M ? abV0M[0] : 0;
      aV0MOnVsOfErr = abV0M ? abV0M[1] : 0;
      bV0MOnVsOfVal = abV0M ? abV0M[2] : 0;
      bV0MOnVsOfErr = abV0M ? abV0M[3] : 0;
      abSPD = getCutParametersOnlineVsOffline(hSPDOnVsOfCln, "SPD");
      aSPDOnVsOfVal = abSPD ? abSPD[0] : 0;
      aSPDOnVsOfErr = abSPD ? abSPD[1] : 0;
      bSPDOnVsOfVal = abSPD ? abSPD[2] : 0;
      bSPDOnVsOfErr = abSPD ? abSPD[3] : 0;
    } else if (label.Contains(" &65536 ")) {
      meanV0MOnHM    = hV0MOnAcc ? hV0MOnAcc->GetMean()      : 0;
      meanV0MOfHM    = hV0MOfAcc ? hV0MOfAcc->GetMean()      : 0;
      meanOFOHM      = hOFOAcc   ? hOFOAcc->GetMean()        : 0;
      meanTKLHM      = hTKLAcc   ? hTKLAcc->GetMean()        : 0;
      meanErrV0MOnHM = hV0MOnAcc ? hV0MOnAcc->GetMeanError() : 0;
      meanErrV0MOfHM = hV0MOfAcc ? hV0MOfAcc->GetMeanError() : 0;
      meanErrOFOHM   = hOFOAcc   ? hOFOAcc->GetMeanError()   : 0;
      meanErrTKLHM   = hTKLAcc   ? hTKLAcc->GetMeanError()   : 0;
      thresholdV0M   = hV0MOnAcc->GetBinLowEdge(hV0MOnAcc->FindFirstBinAbove(9));
      printf("thresholdV0M=%f\n",thresholdV0M);
    }
  }

  for (Int_t j=1;j<=hHistStat->GetNbinsY();j++){
    TString label = hHistStat->GetYaxis()->GetBinLabel(j);
    if (label.Length()<2) continue;
    // skip background triggers
    // TODO introduce identifier to filter-out background triggers
    if      (label.Contains("-A-"))      continue;
    else if (label.Contains("-C-"))      continue;
    else if (label.Contains("-E-"))      continue;
    else if (label.Contains("-AC-"))     continue;
    else if (label.Contains("-ACE-"))    continue;
    else if (label.Contains("-GA-"))     continue;
    else if (label.Contains("-GC-"))     continue;
    else if (label.Contains("1A-ABCE-")) continue;
    else if (label.Contains("1C-ABCE-")) continue;
    else if (label.Contains("C0LSR-ABCE-")) continue;

    // Read mask
    // TODO think how to propagate mask with TBit aliases
    UInt_t mask = 0;
    TString classList = ""; // list of classes for given PS bit
    TObjArray* array = label.Tokenize(" ");
    for (Int_t itoken=0;itoken<array->GetEntries();itoken++){
      TString token = array->At(itoken)->GetName();
      if (itoken==0) classList = token; 
      if (token[0]!='&') continue;
      token.Remove(0,1);
      mask = token.Atoi();
      break;
    }
    array->Delete();
    delete array;
    printf("%s\n",label.Data());
    printf("%i\n",mask);
    if (!mask) continue;
    // Fill all and accepted counters for the most significant bit
    Int_t ibit = TMath::Nint(TMath::Log2(mask));

    // FIXME
    // if (!bitNames[ibit].EqualTo("kINT7") && !bitNames[ibit].EqualTo("kINT7inMUON") && !bitNames[ibit].EqualTo("kMuonUnlikeLowPt7")) continue;

    if (ibit>=NBITS) continue;
//    if (alias_recorded[ibit]) break; 

    alias_reconstructed[ibit] = Int_t(hHistStat->GetBinContent(1,j));
    alias_accepted[ibit]      = Int_t(hHistStat->GetBinContent(2,j));
    if (hHistStat->GetNbinsX()!=19) { 
      // old stat histo without ZDC background monitoring
      // setting ZDCBG step equal to V0AND and take other steps shifted by 1
      alias_acc_step1[ibit]     = Int_t(hHistStat->GetBinContent(3,j));
      alias_acc_step2[ibit]     = Int_t(hHistStat->GetBinContent(3,j)); 
      alias_acc_step3[ibit]     = Int_t(hHistStat->GetBinContent(4,j));
      alias_acc_step4[ibit]     = Int_t(hHistStat->GetBinContent(5,j));
      alias_acc_step5[ibit]     = Int_t(hHistStat->GetBinContent(6,j));
      alias_acc_step6[ibit]     = Int_t(hHistStat->GetBinContent(7,j));
      alias_acc_step7[ibit]     = Int_t(hHistStat->GetBinContent(8,j));
      alias_acc_step8[ibit]     = Int_t(hHistStat->GetBinContent(9,j));
      alias_acc_step9[ibit]     = Int_t(hHistStat->GetBinContent(10,j));
    } else {
      alias_acc_step1[ibit]     = Int_t(hHistStat->GetBinContent(3,j));
      alias_acc_step2[ibit]     = Int_t(hHistStat->GetBinContent(4,j));
      alias_acc_step3[ibit]     = Int_t(hHistStat->GetBinContent(5,j));
      alias_acc_step4[ibit]     = Int_t(hHistStat->GetBinContent(6,j));
      alias_acc_step5[ibit]     = Int_t(hHistStat->GetBinContent(7,j));
      alias_acc_step6[ibit]     = Int_t(hHistStat->GetBinContent(8,j));
      alias_acc_step7[ibit]     = Int_t(hHistStat->GetBinContent(9,j));
      alias_acc_step8[ibit]     = Int_t(hHistStat->GetBinContent(10,j));
      alias_acc_step9[ibit]     = Int_t(hHistStat->GetBinContent(11,j));
    }
    //printf("%4i %8i %8i\n",ibit,alias_reconstructed[ibit],alias_accepted[ibit]);
    
    classList.Remove(0,1); // remove +
    array = classList.Tokenize(",");
    // if trigger bit corresponds to several active classes, just take the last one
    // example: kTRD
    // TODO think about more elegant solution
    for (Int_t i=0;i<array->GetEntriesFast();i++){
      TString token = array->At(i)->GetName();
      AliTriggerClass* cl = (AliTriggerClass*) classes.FindObject(token.Data());
      if (!cl) continue;
      Int_t iclass = classes.IndexOf(cl);
      printf(" %30s",token.Data());
      printf(" %12lli",class_l0b[iclass]);
      printf(" %12lli",class_l0a[iclass]);
      printf(" %12lli",class_l1b[iclass]);
      printf(" %12lli",class_l1a[iclass]);
      printf(" %12lli",class_l2b[iclass]);
      printf(" %12lli",class_l2a[iclass]);
      printf("\n");
      alias_recorded[ibit]      = class_l2a[iclass];
      alias_lifetime[ibit]      = class_lifetime[iclass];
      alias_lumi_recorded[ibit] = class_lumi[iclass];
      if (!alias_recorded[ibit]) continue;
      alias_lumi_reconstructed[ibit] = alias_lumi_recorded[ibit]/alias_recorded[ibit]*alias_reconstructed[ibit];
      alias_lumi_accepted[ibit]      = alias_lumi_recorded[ibit]/alias_recorded[ibit]*alias_accepted[ibit];
    }
    array->Delete();
    delete array;

    // Fill run QA histograms
    const char* bitName = bitNames[ibit].Data();
    TList* list = NULL;
    if (qafilename.Contains("event_stat")) {
      list = (TList*) fin->Get(Form("trigger_histograms_%s/histos",label.Data()));
    } else {
      AliTriggerAnalysis* ta = ps->GetTriggerAnalysis(j-1);
      list = ta->GetHistList();
    }
    if (!list) continue;
    
    gSystem->mkdir(bitName);
    gSystem->cd(bitName);
    TH1F* hV0AAll          = (TH1F*) list->FindObject("fHistV0AAll");
    TH1F* hV0AAcc          = (TH1F*) list->FindObject("fHistV0AAcc");
    TH1F* hV0CAll          = (TH1F*) list->FindObject("fHistV0CAll");
    TH1F* hV0CAcc          = (TH1F*) list->FindObject("fHistV0CAcc");
    TH1F* hADAAll          = (TH1F*) list->FindObject("fHistADAAll");
    TH1F* hADAAcc          = (TH1F*) list->FindObject("fHistADAAcc");
    TH1F* hADCAll          = (TH1F*) list->FindObject("fHistADCAll");
    TH1F* hADCAcc          = (TH1F*) list->FindObject("fHistADCAcc");
    TH1F* hBBAflagsAll     = (TH1F*) list->FindObject("fHistBBAflagsAll");
    TH1F* hBBAflagsAcc     = (TH1F*) list->FindObject("fHistBBAflagsAcc");
    TH1F* hBBCflagsAll     = (TH1F*) list->FindObject("fHistBBCflagsAll");
    TH1F* hBBCflagsAcc     = (TH1F*) list->FindObject("fHistBBCflagsAcc");
    TH1F* hBGAflagsAll     = (TH1F*) list->FindObject("fHistBGAflagsAll");
    TH1F* hBGAflagsAcc     = (TH1F*) list->FindObject("fHistBGAflagsAcc");
    TH1F* hBGCflagsAll     = (TH1F*) list->FindObject("fHistBGCflagsAll");
    TH1F* hBGCflagsAcc     = (TH1F*) list->FindObject("fHistBGCflagsAcc");
    TH1F* hV0MOnAll        = (TH1F*) list->FindObject("fHistV0MOnAll");
    TH1F* hV0MOnAcc        = (TH1F*) list->FindObject("fHistV0MOnAcc");
    TH1F* hV0MOnVHM        = (TH1F*) list->FindObject("fHistV0MOnVHM");
    TH1F* hV0MOfAll        = (TH1F*) list->FindObject("fHistV0MOfAll");
    TH1F* hV0MOfAcc        = (TH1F*) list->FindObject("fHistV0MOfAcc");
    TH1F* hFiredBitsSPD    = (TH1F*) list->FindObject("fHistFiredBitsSPD");
    TH2F* hSPDClsVsTklAll  = (TH2F*) list->FindObject("fHistSPDClsVsTklAll");
    TH2F* hSPDClsVsTklCln  = (TH2F*) list->FindObject("fHistSPDClsVsTklCln");
    TH2F* hV0C012vsTklAll  = (TH2F*) list->FindObject("fHistV0C012vsTklAll");
    TH2F* hV0C012vsTklCln  = (TH2F*) list->FindObject("fHistV0C012vsTklCln");
    TH2F* hV0MOnVsOfAll    = (TH2F*) list->FindObject("fHistV0MOnVsOfAll");
    TH2F* hV0MOnVsOfCln    = (TH2F*) list->FindObject("fHistV0MOnVsOfCln");
    TH2F* hSPDOnVsOfAll    = (TH2F*) list->FindObject("fHistSPDOnVsOfAll");
    TH2F* hSPDOnVsOfCln    = (TH2F*) list->FindObject("fHistSPDOnVsOfCln");
    TH2F* hV0C3vs012All    = (TH2F*) list->FindObject("fHistV0C3vs012All");
    TH2F* hV0C3vs012Cln    = (TH2F*) list->FindObject("fHistV0C3vs012Cln");
    TH1F* hSPDVtxPileupAll = (TH1F*) list->FindObject("fHistSPDVtxPileupAll");
    TH1F* hSPDVtxPileupCln = (TH1F*) list->FindObject("fHistSPDVtxPileupCln");
    TH2F* hVIRvsBCmod4pup  = (TH2F*) list->FindObject("fHistVIRvsBCmod4pup");
    TH2F* hVIRvsBCmod4acc  = (TH2F*) list->FindObject("fHistVIRvsBCmod4acc");
    TH1F* hOFOAll          = (TH1F*) list->FindObject("fHistOFOAll");
    TH1F* hOFOAcc          = (TH1F*) list->FindObject("fHistOFOAcc");
    TH1F* hOFOVHM          = (TH1F*) list->FindObject("fHistOFOVHM");
    TH1F* hTKLAll          = (TH1F*) list->FindObject("fHistTKLAll");
    TH1F* hTKLAcc          = (TH1F*) list->FindObject("fHistTKLAcc");
    TH2F* hAD              = (TH2F*) list->FindObject("fHistAD");
    TH1F* hTimeZNA         = (TH1F*) list->FindObject("fHistTimeZNA");
    TH1F* hTimeZNC         = (TH1F*) list->FindObject("fHistTimeZNC");
    TH2F* hTimeZNSumVsDif  = (TH2F*) list->FindObject("fHistTimeZNSumVsDif");
    TH2F* hTimeCorrZDC     = (TH2F*) list->FindObject("fHistTimeCorrZDC");
    TH2F* hOFOvsTKLAcc     = (TH2F*) list->FindObject("fHistOFOvsTKLAcc");
    TH2F* hV0MOnVsOfAcc    = (TH2F*) list->FindObject("fHistV0MOnVsOfAcc");
        
    DrawV0MOnVsOf(hV0MOnVsOfAll,hV0MOnVsOfCln,bitName,abV0M);
    DrawSPDOnVsOf(hSPDOnVsOfAll,hSPDOnVsOfCln,bitName,abSPD);
    DrawSPDClsVsTkl(hSPDClsVsTklAll,hSPDClsVsTklCln,bitName);
    DrawV0C012vsTkl(hV0C012vsTklAll,hV0C012vsTklCln,bitName);
    DrawVIR(hVIRvsBCmod4pup,hVIRvsBCmod4acc,bitName);
    DrawV0C3vs012(hV0C3vs012All,hV0C3vs012Cln,bitName);

    DrawTiming(hV0AAll,hV0AAcc,hV0CAll,hV0CAcc,"V0",bitName);
    DrawTiming(hADAAll,hADAAcc,hADCAll,hADCAcc,"AD",bitName);
    DrawFlags(hBBAflagsAll,hBBAflagsAcc,hBBCflagsAll,hBBCflagsAcc,hBGAflagsAll,hBGAflagsAcc,hBGCflagsAll,hBGCflagsAcc,bitName);
    DrawMultiplicity(hV0MOnAll,hV0MOnAcc,hV0MOnVHM,hV0MOfAll,hV0MOfAcc,bitName,"V0M","V0M",meanV0MOn,meanV0MOf);
    DrawMultiplicity(hOFOAll,hOFOAcc,hOFOVHM,hTKLAll,hTKLAcc,bitName,"OFO","TKL",meanOFO,meanTKL);
    DrawZDC(hTimeZNA,hTimeZNC,hTimeZNSumVsDif,hTimeCorrZDC,bitName);
    if (bitNames[ibit].EqualTo("kINT7")) DrawEfficiency(hOFOvsTKLAcc,"OFO",bitName);
    if (bitNames[ibit].EqualTo("kINT7")) DrawEfficiency(hV0MOnVsOfAcc,"V0M",bitName);
//
//    if (hFiredBitsSPD) {
//      TCanvas* cFiredBitsSPD = new TCanvas(Form("cFiredBitsSPD_%s",bitName),Form("cFiredBitsSPD_%s",bitName),1800,500);
//      gPad->SetLogy();
//      gPad->SetMargin(0.05,0.01,0.12,0.06);
//      hFiredBitsSPD->SetTitle(Form("%s: hardware FO",bitName));
//      hFiredBitsSPD->SetTitleFont(43);
//      hFiredBitsSPD->SetTitleSize(25);
//      hFiredBitsSPD->GetYaxis()->SetTitleFont(43);
//      hFiredBitsSPD->GetXaxis()->SetLabelFont(43);
//      hFiredBitsSPD->GetYaxis()->SetLabelFont(43);
//      hFiredBitsSPD->GetYaxis()->SetTitleSize(25);
//      hFiredBitsSPD->GetXaxis()->SetLabelSize(25);
//      hFiredBitsSPD->GetYaxis()->SetLabelSize(25);
//      hFiredBitsSPD->GetYaxis()->SetTickLength(0.01);
//      hFiredBitsSPD->GetYaxis()->SetTitleOffset(0.5);
//      hFiredBitsSPD->GetYaxis()->SetDecimals(1);
//      hFiredBitsSPD->SetLineWidth(2);
//      hFiredBitsSPD->SetLineColor(kBlue);
//      hFiredBitsSPD->Draw();
//      gPad->Print(Form("%s_FiredBitsSPD.png",bitName));
//      hFiredBitsSPD->Write(Form("%s_FiredBitsSPD",bitName));
//    } else printf("QA histogram not found\n"); 
//
//    if (hBitsSPD) {
//      TCanvas* cBitsSPD = new TCanvas(Form("cBitsSPD_%s",bitName),Form("cBitsSPD_%s",bitName),800,800);
//      gPad->SetLogz();
//      gPad->SetMargin(0.12,0.12,0.10,0.06);
//      hBitsSPD->SetTitle(Form("%s: hardware FO vs offline FO",bitName));
//      hBitsSPD->GetXaxis()->SetTitleOffset(1.3);
//      hBitsSPD->GetYaxis()->SetTitleOffset(1.6);
//      hBitsSPD->Draw("colz");
//      gPad->Print(Form("%s_BitsSPD.png",bitName));
//      hBitsSPD->Write(Form("%s_BitsSPD",bitName));
//    } else printf("QA histogram not found\n"); 
    gSystem->cd("..");
  }
  
  writeTree(fout,t);
  return 0;
}

void DrawEfficiency(TH2F* h2D, const char* multOn, const char* bitName){
  if (!h2D) {
    printf("h2D histo for HM trigger efficiency monitoring not found\n");
    return;
  }
  TCanvas* c = new TCanvas(Form("c_%s_%s",multOn,bitName),Form("c_%s_%s",multOn,bitName),1000,1000);
  c->Divide(2,2,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  h2D->SetTitle(bitName);
  h2D->GetYaxis()->SetTitleOffset(1.5);
  h2D->Draw("colz");
  TH2F* h2Dcum = (TH2F*) h2D->Clone();
  h2Dcum->Clear();
  TH1F** proj = new TH1F*[h2D->GetNbinsY()];
  TH1F** eff  = new TH1F*[h2D->GetNbinsY()];
  TH1F* hEff90VsThreshold = new TH1F(Form("hEff90VsThreshold_%s",multOn),Form(";Online %s threshold",multOn),h2D->GetNbinsY(),0,h2D->GetYaxis()->GetXmax());
  TH1F* hEff95VsThreshold = new TH1F(Form("hEff95VsThreshold_%s",multOn),Form(";Online %s threshold",multOn),h2D->GetNbinsY(),0,h2D->GetYaxis()->GetXmax());
  
  TLegend* leg = new TLegend(0.6,0.55,0.98,0.85);

  for (Int_t on=1;on<=h2D->GetNbinsY();on++){
    for (Int_t of=1;of<=h2D->GetNbinsX();of++){
      h2Dcum->SetBinContent(of,on,h2D->Integral(of,of,on,h2D->GetNbinsY()+1));
    }
    proj[on-1] = (TH1F*) h2Dcum->ProjectionX(Form("hProj_%i",on),on,on);
    proj[on-1]->SetLineWidth(2);
    eff[on-1] = (TH1F*) proj[on-1]->Clone(Form("hEff_%i",on));
    eff[on-1]->Divide(proj[0]);
    Int_t bin90 = eff[on-1]->FindFirstBinAbove(0.90);
    Int_t bin95 = eff[on-1]->FindFirstBinAbove(0.95);
    Float_t all     = proj[on-1]->Integral(1,proj[on-1]->GetNbinsX()+1);
    Float_t above90 = proj[on-1]->Integral(bin90,proj[on-1]->GetNbinsX()+1);
    Float_t above95 = proj[on-1]->Integral(bin95,proj[on-1]->GetNbinsX()+1);
    hEff90VsThreshold->SetBinContent(on,all>0 ? above90/all : 0);
    hEff95VsThreshold->SetBinContent(on,all>0 ? above95/all : 0);
    if (on==1){
      c->cd(2);
      proj[0]->SetTitle(bitName);
      proj[0]->SetLineColor(kBlack);
      proj[0]->DrawCopy();
    }
    if ((TString(multOn).Contains("OFO") && (on==60 || on==70  || on== 80 || on== 90)) ||
        (TString(multOn).Contains("V0M") && (on==60 || on==120 || on==180 || on==240))) {
      proj[on-1]->SetLineColor(TString(multOn).Contains("OFO") ? (on-60)/10+2 : on/60+1);
      eff[on-1]->SetLineColor(proj[on-1]->GetLineColor());
      eff[on-1]->SetTitle(bitName);
      c->cd(2);
      proj[on-1]->DrawCopy("same");
      leg->AddEntry(proj[on-1],Form("%s>=%.0f\n",multOn,h2Dcum->GetYaxis()->GetBinUpEdge(on)));
      c->cd(3);
      eff[on-1]->DrawCopy(on==60 ? "": "same");
    }
  }
  
  c->cd(2);
  gPad->SetLogy();
  leg->Draw();
  c->cd(3);
  gPad->SetLogy();
  leg->Draw();
  
  c->cd(4);
  hEff90VsThreshold->SetLineWidth(2);
  hEff95VsThreshold->SetLineWidth(2);
  hEff90VsThreshold->SetLineColor(kBlue);
  hEff95VsThreshold->SetLineColor(kRed);
//  hEff90VsThreshold->GetXaxis()->SetRangeUser(0,150);
  hEff90VsThreshold->DrawCopy();
  hEff95VsThreshold->DrawCopy("same");
  
  TLegend* legEff = new TLegend(0.25,0.75,0.98,0.89);
  legEff->AddEntry(hEff90VsThreshold,"Purity at 90% efficiency");
  legEff->AddEntry(hEff95VsThreshold,"Purity at 95% efficiency");
  legEff->Draw();

  c->Print(Form("%s_efficiency.png",multOn));
}

/*
void DrawV0MOnVsOf(TH2F* hAll,TH2F* hCln,const char* bitName){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }

  TCanvas* c = new TCanvas(Form("cV0MOnVsOf_%s",bitName),Form("cV0MOnVsOf_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->SetTitle(Form("%s",bitName));
  hAll->GetYaxis()->SetTitleOffset(1.7);
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  hCln->SetTitle(Form("%s",bitName));
  hCln->GetYaxis()->SetTitleOffset(1.7);
  hCln->Draw("colz");
  c->Print("V0MOnVsOf.png");
//  TH1F* hFitResults = (TH1F*) hCln->ProjectionX("hFitResults");
//  TF1* fExp = new TF1("fExp","exp(-([0]-x)*[1])",0,10000);
//  TF1* fPol1 = new TF1("fPol1","pol1",0,hCln->GetXaxis()->GetXmax());
//  fExp->SetLineColor(kMagenta);
//  hFitResults->Reset();
//
//  hCln->Sumw2();
//  TCanvas* cFit = new TCanvas(Form("cFit_%s",bitName),Form("cFit_%s",bitName),1800,800);
//  cFit->Divide(2,1,0.001,0.001);
//  cFit->cd(1);
//  gPad->SetLogy();
//  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
//  for (Int_t i=1;i<=hCln->GetNbinsX();i++){
//    printf("%i\n",i);
//    TH1D* hProj = hCln->ProjectionY(Form("proj%i",i),i,i);
//    hProj->GetXaxis()->SetTitleOffset(1.0);
//    hProj->SetTitle(Form("%s;Online V0M;",bitName));
//    Float_t integ = hProj->Integral();
//    gPad->Print("fit.png");
//    if (integ<100) continue;
//    Double_t mMax = hProj->GetBinLowEdge(hProj->GetMaximumBin()-1);
//    Double_t mMin = hProj->GetBinLowEdge(hProj->GetMaximumBin()-6);
//    fExp->SetRange(mMin,mMax);
//    hProj->Fit("fExp","QN","",mMin,mMax);
//    if (i==20 || i==40 || i==60 || i==80 || i==100) {
//      hProj->SetLineColor(i/20);
//      hProj->DrawCopy(i==20 ? "" : "same");
//      fExp->DrawCopy("lsame");
//
//      leg->AddEntry(hProj,Form("Offline V0M = %.0f",hCln->GetXaxis()->GetBinCenter(i)));
//    }
//    
//    Double_t x0     = fExp->GetParameter(0);
//    Double_t lambda = 1/fExp->GetParameter(1);
//    Double_t x0err  = fExp->GetParError(0);
//    Double_t lambdaerr = lambda*fExp->GetParError(1)/fExp->GetParameter(1);
//    hFitResults->SetBinContent(i,x0-7*lambda);
//    hFitResults->SetBinError(i,x0err+7*lambdaerr);
//    gPad->Print("fit.png");
//  }
//  cFit->WaitPrimitive();
//  return;
//  leg->Draw();
//  
//  cFit->cd(2);
//  gPad->SetLogz();
//  gPad->SetLeftMargin(0.11);
//  gPad->SetRightMargin(0.11);
//  hCln->GetListOfFunctions()->Clear();
//  hCln->DrawCopy("colz");
//  hFitResults->SetLineWidth(2);
//  hFitResults->SetLineColor(kMagenta+1);
//  hFitResults->SetMarkerColor(kMagenta+1);
//  hFitResults->SetMarkerStyle(kOpenCircle);
//  hFitResults->DrawCopy("psame");
//  hFitResults->Fit(fPol1,"QN","",50,hCln->GetXaxis()->GetXmax());
//  printf("par[0]=%f\n",fPol1->GetParameter(0));
//  printf("par[1]=%f\n",fPol1->GetParameter(1));
//  fPol1->Draw("same");
//  cFit->Print(Form("%s_Fit.png",bitName));
}
*/


void DrawMultiplicity(TH1F* hOnAll, TH1F* hOnAcc, TH1F* hOnVHM, TH1F* hOfAll,  TH1F* hOfAcc, const char* bitName, TString multOn, TString multOf,
 Float_t meanOn, Float_t meanOf
){
  if (!hOnAll || !hOfAll || !hOnAcc || !hOfAcc) {
    printf("QA histogram not found\n");
    printf("%p %p %p %p\n",hOnAll,hOfAll,hOnAcc,hOfAcc);
    return;
  }

  hOnAll->SetLineColor(kBlack);
  hOfAll->SetLineColor(kBlack);
  hOnAcc->SetLineColor(kBlue);
  hOfAcc->SetLineColor(kBlue);
  hOnAcc->SetLineWidth(2);
  hOnAll->SetLineWidth(2);
  hOfAll->SetLineWidth(2);
  hOfAcc->SetLineWidth(2);
  hOfAcc->SetMarkerColor(kBlue);
  hOnAcc->SetMarkerColor(kBlue);
  hOnAcc->SetFillColor(kBlue);
  hOfAcc->SetFillColor(kBlue);

  hOnAll->SetNameTitle(Form("%s_%s_All",bitName,multOn.Data()),bitName);
  hOfAll->SetNameTitle(Form("%s_%s_All",bitName,multOf.Data()),bitName);

  TLegend* leg1 = new TLegend(0.40,0.80,0.97,0.90);
  leg1->AddEntry(hOnAll,Form("All: %.0f",hOnAll->Integral(0,hOnAll->GetNbinsX()+1)));
  leg1->AddEntry(hOnAcc,Form("Accepted: %.0f, <%s>=%.1f",hOnAcc->Integral(0,hOnAcc->GetNbinsX()+1),hOnAcc->GetMean(),multOn.Data()),"l");

  TLegend* leg2 = new TLegend(0.40,0.80,0.97,0.90);
  leg2->AddEntry(hOfAll,Form("All: %.0f",hOfAll->Integral(0,hOfAll->GetNbinsX()+1)));
  leg2->AddEntry(hOfAcc,Form("Accepted: %.0f, <%s>=%.1f",hOfAcc->Integral(0,hOfAcc->GetNbinsX()+1),hOfAcc->GetMean(),multOf.Data()),"l");

  TCanvas* c = new TCanvas(Form("%s_%s",bitName,multOn.Data()),Form("c_%s",bitName),1800,900);
  c->Divide(2,1,0.001,0.001);

  c->cd(1);
  gPad->SetLogy();
  hOnAll->SetMinimum(0.9);
  hOnAll->Draw();
  hOnAcc->Draw("same");
  if (hOnVHM) hOnVHM->Draw("same");
  leg1->Draw();
  
  c->cd(2);
  gPad->SetLogy();
  hOfAll->SetMinimum(0.9);
  hOfAll->Draw();
  hOfAcc->Draw("same");
  leg2->Draw();

  c->Print(Form("%s.png",multOn.Data()));
  
  TCanvas* cNorm = new TCanvas(Form("cNorm_%s_%s",bitName,multOn.Data()),Form("cNorm_%s",bitName),1800,900);
  cNorm->Divide(2,1,0.001,0.001);
  
  TH1F* hOnNorm = new TH1F(Form("hOnNorm_%s_%s",bitName,multOn.Data()),Form("%s;%s/<%s>;",bitName,multOn.Data(),multOn.Data()),50,0,10);
  TH1F* hOfNorm = new TH1F(Form("hOfNorm_%s_%s",bitName,multOf.Data()),Form("%s;%s/<%s>;",bitName,multOf.Data(),multOf.Data()),50,0,10);

  for (Int_t i=1;i<=hOnAcc->GetNbinsX();i++){
    Float_t binCenter  = hOnAcc->GetXaxis()->GetBinCenter(i);
    Float_t binContent = hOnAcc->GetBinContent(i);
    hOnNorm->Fill(meanOn>1e-10 ? binCenter/meanOn : 0,binContent);
  }
  for (Int_t i=1;i<=hOfAcc->GetNbinsX();i++){
    Float_t binCenter  = hOfAcc->GetXaxis()->GetBinCenter(i);
    Float_t binContent = hOfAcc->GetBinContent(i);
    hOfNorm->Fill(meanOf>1e-10 ? binCenter/meanOf : 0,binContent);
  }

  hOnNorm->SetLineColor(kBlue);
  hOfNorm->SetLineColor(kBlue);
  hOnNorm->SetLineWidth(2);
  hOfNorm->SetLineWidth(2);

  cNorm->cd(1);
  gPad->SetLogy();
  hOnNorm->DrawCopy();

  cNorm->cd(2);
  gPad->SetLogy();
  hOfNorm->DrawCopy();

  cNorm->Print(Form("%s_norm.png",multOn.Data()));

  TFile* fout = new TFile(Form("fout_%s.root",multOn.Data()),"recreate");
  hOfNorm->Write();
  hOnNorm->Write();
  fout->Close();

  TH1F* hOnAllCum  = (TH1F*) hOnAll->GetCumulative(kFALSE);
  TH1F* hOnAccCum  = (TH1F*) hOnAcc->GetCumulative(kFALSE);
  TH1F* hOfAllCum  = (TH1F*) hOfAll->GetCumulative(kFALSE);
  TH1F* hOfAccCum  = (TH1F*) hOfAcc->GetCumulative(kFALSE);
  TH1F* hOnNormCum = (TH1F*) hOnNorm->GetCumulative(kFALSE);
  TH1F* hOfNormCum = (TH1F*) hOfNorm->GetCumulative(kFALSE);
  TH1F* hOnVHMCum = hOnVHM ? (TH1F*) hOnVHM->GetCumulative(kFALSE) : 0;
  
  hOnAllCum->Scale(hOnAllCum->GetBinContent(1)>1e-10 ? 1./hOnAllCum->GetBinContent(1) : 1);
  hOnAccCum->Scale(hOnAccCum->GetBinContent(1)>1e-10 ? 1./hOnAccCum->GetBinContent(1) : 1);
  hOfAllCum->Scale(hOfAllCum->GetBinContent(1)>1e-10 ? 1./hOfAllCum->GetBinContent(1) : 1);
  hOfAccCum->Scale(hOfAccCum->GetBinContent(1)>1e-10 ? 1./hOfAccCum->GetBinContent(1) : 1);
  hOnNormCum->Scale(hOnNormCum->GetBinContent(1)>1e-10 ? 1./hOnNormCum->GetBinContent(1) : 1);
  hOfNormCum->Scale(hOfNormCum->GetBinContent(1)>1e-10 ? 1./hOfNormCum->GetBinContent(1) : 1);
  if (hOnVHMCum) hOnVHMCum->Scale(hOnVHMCum->GetBinContent(1)>1e-10 ? 1./hOnVHMCum->GetBinContent(1) : 1);
  
  hOnAccCum->SetFillStyle(0);
  hOfAccCum->SetFillStyle(0);
  hOnAllCum->SetNameTitle(Form("%s_%s_All_Cum",bitName,multOn.Data()),Form("%s: cumulative",bitName));
  hOfAllCum->SetNameTitle(Form("%s_%s_All_Cum",bitName,multOf.Data()),Form("%s: cumulative",bitName));
  hOnNormCum->SetNameTitle(Form("%s_%s_Norm_Cum",bitName,multOn.Data()),Form("%s: cumulative",bitName));
  hOfNormCum->SetNameTitle(Form("%s_%s_Norm_Cum",bitName,multOf.Data()),Form("%s: cumulative",bitName));
  
  TCanvas* ccum = new TCanvas(Form("cCum_%s_%s",bitName,multOn.Data()),Form("cCum_%s",bitName),1800,1800);
  ccum->Divide(2,2,0.001,0.001);
  
  ccum->cd(1);
  gPad->SetLogy();
  hOnAllCum->DrawCopy();
  hOnAccCum->DrawCopy("same");
  if (hOnVHMCum) hOnVHMCum->DrawCopy("same");
  leg1->Draw();
  
  ccum->cd(2);
  gPad->SetLogy();
  hOfAllCum->DrawCopy();
  hOfAccCum->DrawCopy("same");
  leg2->Draw();
  
  ccum->cd(3);
  gPad->SetLogy();
  hOnNormCum->DrawCopy();
  
  ccum->cd(4);
  gPad->SetLogy();
  hOfNormCum->DrawCopy();
  
  ccum->Print(Form("%s_cum.png",multOn.Data()));
  
}


void DrawSPDClsVsTkl(TH2F* hAll,TH2F* hCln, const char* bitName){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }
  hAll->SetTitle(bitName);
  hCln->SetTitle(bitName);
  hAll->GetYaxis()->SetTitleOffset(1.6);
  hCln->GetYaxis()->SetTitleOffset(1.6);

  TCanvas* c = new TCanvas(Form("cSPDClsVsTkl_%s",bitName),Form("cSPDClsVsTkl_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hCln->Draw("colz");
  c->Print("ClsVsTkl.png");
}


void DrawV0C012vsTkl(TH2F* hAll,TH2F* hCln, const char* bitName){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }
  hAll->SetTitle(bitName);
  hCln->SetTitle(bitName);
  hAll->GetYaxis()->SetTitleOffset(1.3);
  hCln->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* c = new TCanvas(Form("cV0C012vsTkl_%s",bitName),Form("cV0C012vsTkl_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hCln->Draw("colz");
  c->Print("V0C012vsTkl.png");
} 

void DrawV0C3vs012(TH2F* hAll,TH2F* hCln, const char* bitName){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }
  hAll->SetTitle(bitName);
  hCln->SetTitle(bitName);
  hAll->GetYaxis()->SetTitleOffset(1.3);
  hCln->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* c = new TCanvas(Form("cV0C3vs012_%s",bitName),Form("cV0C3vs012_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hCln->Draw("colz");
  c->Print("V0C3vs012.png");
} 



void DrawVIR(TH2F* hPup,TH2F* hAcc, const char* bitName){
  if (!hPup || !hAcc) {
    printf("QA histogram not found\n");
    return;
  }
  hPup->SetTitle(bitName);
  hAcc->SetTitle(bitName);
  hPup->GetYaxis()->SetTitleOffset(1.3);
  hAcc->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* c = new TCanvas(Form("cVIR_%s",bitName),Form("cVIR_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hPup->Draw("colz");
  c->cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAcc->Draw("colz");
  c->Print("VIR.png");
} 


void DrawTiming(TH1F* hAall,TH1F* hAacc,TH1F* hCall,TH1F* hCacc, TString cname, const char* bitName){
  if (!hAall || !hAacc || !hCall || !hCacc) {
    printf("QA histogram not found\n");
    return;
  }
  hCall->SetTitle(bitName);
  hAall->SetTitle(bitName);
  hAall->SetLineWidth(2);
  hCall->SetLineWidth(2);
  hAall->SetLineColor(kBlack);
  hCall->SetLineColor(kBlack);
  hAacc->SetLineColor(kBlue);
  hCacc->SetLineColor(kBlue);
  hAacc->SetFillColor(kBlue);
  hCacc->SetFillColor(kBlue);

  TLegend* leg = new TLegend(0.15,0.8,0.45,0.9);
  leg->AddEntry(hAall,"All events");
  leg->AddEntry(hAacc,"Accepted events");
  
  TCanvas* c = new TCanvas(Form("%s_timing_%s",cname.Data(),bitName),Form("%s_timing_%s",cname.Data(),bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  
  c->cd(1);
  gPad->SetLogy();
  hAall->Draw();
  hAacc->Draw("same");
  leg->Draw();

  c->cd(2);
  gPad->SetLogy();
  hCall->Draw();
  hCacc->Draw("same");
  leg->Draw();
  
  c->Print(Form("%s.png",cname.Data()));
}

void DrawFlags(TH1F* hBBAflagsAll,TH1F* hBBAflagsAcc,TH1F* hBBCflagsAll,TH1F* hBBCflagsAcc,
               TH1F* hBGAflagsAll,TH1F* hBGAflagsAcc,TH1F* hBGCflagsAll,TH1F* hBGCflagsAcc, const char* bitName){
  if (!hBBAflagsAll || !hBBAflagsAcc || !hBBCflagsAll || !hBBCflagsAcc || !hBGAflagsAll || !hBGAflagsAcc || !hBGCflagsAll || !hBGCflagsAcc) {
    printf("QA histogram not found\n");
    return;
  }
  hBBAflagsAll->SetMinimum(0.9);
  hBBCflagsAll->SetMinimum(0.9);
  hBGAflagsAll->SetMinimum(0.9);
  hBGCflagsAll->SetMinimum(0.9);
  hBBAflagsAll->SetTitle(bitName);
  hBBCflagsAll->SetTitle(bitName);
  hBGAflagsAll->SetTitle(bitName);
  hBGCflagsAll->SetTitle(bitName);
  hBBAflagsAll->SetLineColor(kBlack);
  hBBCflagsAll->SetLineColor(kBlack);
  hBGAflagsAll->SetLineColor(kBlack);
  hBGCflagsAll->SetLineColor(kBlack);
  hBBAflagsAcc->SetLineColor(kBlue);
  hBBCflagsAcc->SetLineColor(kBlue);
  hBGAflagsAcc->SetLineColor(kBlue);
  hBGCflagsAcc->SetLineColor(kBlue);
  hBBAflagsAcc->SetFillColor(kBlue);
  hBBCflagsAcc->SetFillColor(kBlue);
  hBGAflagsAcc->SetFillColor(kBlue);
  hBGCflagsAcc->SetFillColor(kBlue);

  TCanvas* c = new TCanvas(Form("c_V0flags_%s",bitName),Form("c_%s",bitName),1000,900);
  c->Divide(2,2,0.001,0.001);
  c->cd(1);
  gPad->SetLogy();
  hBBAflagsAll->Draw();
  hBBAflagsAcc->Draw("same");
  c->cd(2);
  gPad->SetLogy();
  hBBCflagsAll->Draw();
  hBBCflagsAcc->Draw("same");
  c->cd(3);
  gPad->SetLogy();
  hBGAflagsAll->Draw();
  hBGAflagsAcc->Draw("same");
  c->cd(4);
  gPad->SetLogy();
  hBGCflagsAll->Draw();
  hBGCflagsAcc->Draw("same");
  c->Print("V0flags.png");
}

void DrawZDC(TH1F* hTimeZNA, TH1F* hTimeZNC, TH2F* hTimeZNSumVsDif, TH2F* hTimeCorrZDC,const char* bitName){
  if (!hTimeZNA || !hTimeZNC || !hTimeZNSumVsDif || !hTimeCorrZDC) {
    printf("ZDC histogram not found: ");
    printf("%p %p %p %p\n",hTimeZNA,hTimeZNC,hTimeZNSumVsDif,hTimeCorrZDC);
    return;
  }
  hTimeZNA       ->SetTitle(bitName);
  hTimeZNC       ->SetTitle(bitName);
  hTimeCorrZDC   ->SetTitle(bitName);
  hTimeZNSumVsDif->SetTitle(bitName);
  hTimeCorrZDC->GetYaxis()->SetTitleOffset(1.3);
  hTimeZNSumVsDif->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* cZDC = new TCanvas(Form("cZDC_%s",bitName),Form("cZDC_%s",bitName),1000,1000);
  cZDC->Divide(2,2,0.001,0.001);
  
  cZDC->cd(1);
  gPad->SetLogy();
  hTimeZNA->Draw();
  
  cZDC->cd(2);
  gPad->SetLogy();
  hTimeZNC->Draw();
  
  cZDC->cd(3);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();

  hTimeCorrZDC->Draw("colz");
  
  cZDC->cd(4);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hTimeZNSumVsDif->Draw("colz");
  
  cZDC->Print("ZDC.png");
}


/*
Moving forwards (backwards) from the maximum value which is the first bin to have a value below threshold*max

Parameters
----------
h : TH1
    The histogram in question
threshold : float
    The threshold factor below which the bin content has to be
direction : Int_t
    {-1, 1}. If -1, move backwards from max, else move forwards

Returns
-------
Int_t :
    Bin number fullfilling the requirement or the first (last) bin of the histogram with content
*/
Int_t findBinNtimesLessThanMax(const TH1& h, Float_t threshold, Int_t direction) {
  Int_t ibin = h.GetMaximumBin();
  Float_t max_val = h.GetBinContent(ibin);
  while  (true) {
    ibin += direction;
    if (ibin < 1 || ibin > h.GetNbinsX() || h.GetBinContent(ibin) < 1) {
      return ibin + direction; // we went to far...
    }
    if (h.GetBinContent(ibin) <= max_val*threshold) {
      return ibin;
    }
  }
}

/*
Helper function to set a point and its yerror in one go
*/
void setPointWithError(TGraphErrors *gr, Int_t idx, Float_t x, Float_t y, Float_t yerr) {
  gr->SetPoint(idx, x, y);
  gr->SetPointError(idx, 0, yerr / 2.0); // not clear if this is the total error or the error in one direction...
}

/* 
Extract Online vs. Offline cut values.

The logic for extracting these values depends on the detector (SPD or V0). 

Parameters
----------
h2d : TH2F*
    Online vs. Offline histogram with previous cuts already applied ("clean")
detector : std::string
    Name of the detector. {SPD, V0M}

Returns
-------
Double_t[4] :
    [offset, offset_err, slope, slope_err] or NULL if no parameters could be extracted
 */
Double_t* getCutParametersOnlineVsOffline(TH2* h2d, const std::string detector) {
  // Fit ranges
  Int_t low_bin, up_bin;
  Double_t low_edge, up_edge;
  // N*width of the cut in terms of sigma (SPD) lambda (V0M)
  Int_t cut_width_factor;
  // Function used to fit individual slice
  TF1 slice_fit;
  // Minimum number of entries per slice
  Int_t min_entries_per_slice;
  
  if (detector == "SPD") {
    slice_fit = TF1("slice_fit", "[2]*exp(-0.5*((x-[0])/[1])**2)", 0, h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX()));
    min_entries_per_slice = 100;
    cut_width_factor = 5;
  }
  else if (detector == "V0M") {
    // Note that parameter 2 will later be fixed. prior to the fit
    slice_fit = TF1("slice_fit","[2]*exp(-([0]-x)/[1])", 0, h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX()));
    min_entries_per_slice = 200;
    cut_width_factor = 7;
  }
  else {
    return NULL;
  }
  // Set default parameters and restrict the width to values > 0
  slice_fit.SetParameters(1, 1, 1);
  // Ensure that the width parameter > 0
  slice_fit.SetParLimits(1, 0, 1e4);
  TF1 fit_lin_mean = TF1("linear_fit_means","pol1", 0, h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX()));
  TF1 fit_lin_width = TF1("linear_fit_widths","pol1", 0, h2d->GetXaxis()->GetBinUpEdge(h2d->GetNbinsX()));
  TGraphErrors widths = TGraphErrors();
  TGraphErrors means = TGraphErrors();
  TGraphErrors nevents = TGraphErrors();
  Int_t valid_fits = 0;
  // Loop through each slice along the offline axis.  Fit each slice
  // with the appropriate function and add this point to a graph which
  // is latter fitted with a linear regression
  for (Int_t ibin = 1; ibin <= h2d->GetNbinsX(); ibin++) {
    TH1* h_py = h2d->ProjectionY("slice", ibin, ibin+1);
    if (h_py->Integral() < min_entries_per_slice) {
      // It seems like ~100 events is about the thresholds where the fitting errors skyrocket
      continue;
    }
    if (detector == "V0M") {
      // Alice and Evgeny move to one lower bin for the upper edge! ie h_py.GetMaximumBin() - 1
      up_bin = h_py->GetMaximumBin() - 1;
      low_bin = findBinNtimesLessThanMax(*h_py, 0.1, -1);
      if (up_bin - low_bin < 5) {
    continue;
      }
      up_edge = h_py->GetXaxis()->GetBinCenter(up_bin);
      low_edge = TMath::Min(0.85 * up_edge, h_py->GetXaxis()->GetBinCenter(h_py->GetMaximumBin() -3));
      // Fix parameter [2]
      Float_t scalling = h_py->GetBinContent(up_bin);
      slice_fit.SetParameters(up_edge, up_edge/20.0, scalling);
      slice_fit.FixParameter(2, scalling);
    }
    if (detector =="SPD") {
      up_bin = findBinNtimesLessThanMax(*h_py, 0, 1); // last filled bin
      low_edge = h_py->GetXaxis()->GetBinLowEdge(h_py->GetMaximumBin() - 1);
      up_edge = h_py->GetXaxis()->GetBinUpEdge(up_bin);
      slice_fit.SetParameter(2, h_py->GetBinContent(h_py->GetMaximumBin()));
    }
    // returns negative value in case of an error not connected with the minimization procedure
    Int_t fit_res = h_py->Fit(&slice_fit, "QRSLL", "", low_edge, up_edge);
    if (fit_res != 0 || slice_fit.GetNDF() == 0) {
      continue;
    }
    setPointWithError(&means, valid_fits,
             h2d->GetXaxis()->GetBinCenter(ibin), slice_fit.GetParameter(0),
             slice_fit.GetParError(0));
    setPointWithError(&widths, valid_fits,
             h2d->GetXaxis()->GetBinCenter(ibin), slice_fit.GetParameter(1),
             slice_fit.GetParError(1));
    nevents.SetPoint(valid_fits, h2d->GetXaxis()->GetBinCenter(ibin), h_py->Integral());
    valid_fits++;
  }
  // Quality assurance after each slice has been tried to fit
  // Require a minimum number of successful fits in order to perform a linear regression
  // Technically, only 2 points are necessary, of cause.
  if (valid_fits < 3) {
    return NULL;
  }
  means.Fit(&fit_lin_mean, "QS");
  widths.Fit(&fit_lin_width, "QRS", "", 0, widths.GetXaxis()->GetXmax());
  Double_t *ret = new Double_t[4];
  ret[0] = fit_lin_mean.GetParameter(0) - cut_width_factor * fit_lin_width.GetParameter(0);
  ret[1] = fit_lin_mean.GetParError(0) + cut_width_factor*fit_lin_width.GetParError(0);
  ret[2] = fit_lin_mean.GetParameter(1) - cut_width_factor*fit_lin_width.GetParameter(1);
  ret[3] = fit_lin_mean.GetParError(1) + cut_width_factor*fit_lin_width.GetParError(1);
  return ret;
}

void DrawV0MOnVsOf(TH2F* hAll,TH2F* hCln,const char* bitName, Double_t* ab=0){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }

  TCanvas* c = new TCanvas(Form("cV0MOnVsOf_%s",bitName),Form("cV0MOnVsOf_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->SetTitle(Form("%s",bitName));
  hAll->GetYaxis()->SetTitleOffset(1.7);
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLogz();
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  hCln->SetTitle(Form("%s",bitName));
  hCln->GetYaxis()->SetTitleOffset(1.7);
  hCln->Draw("colz");
  if (ab) {
    TF1 pol1("fpol1","pol1",hCln->GetXaxis()->GetXmin(),hCln->GetXaxis()->GetXmax());
    pol1.SetLineColor(kMagenta);
    pol1.SetParameter(0,ab[0]);
    pol1.SetParameter(1,ab[2]);
    pol1.Draw("same");
  }
  c->Print("V0MOnVsOf.png");

}


void DrawSPDOnVsOf(TH2F* hAll,TH2F* hCln, const char* bitName, Double_t* ab=0){
  if (!hAll || !hCln) {
    printf("QA histogram not found\n");
    return;
  }
  
  hAll->SetTitle(bitName);
  hCln->SetTitle(bitName);
  hAll->GetYaxis()->SetTitleOffset(1.3);
  hCln->GetYaxis()->SetTitleOffset(1.3);

  TCanvas* c = new TCanvas(Form("cSPDOnVsOf_%s",bitName),Form("cSPDOnVsOf_%s",bitName),1800,800);
  c->Divide(2,1,0.001,0.001);
  c->cd(1);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hAll->Draw("colz");
  c->cd(2);
  gPad->SetLeftMargin(0.11);
  gPad->SetRightMargin(0.11);
  gPad->SetLogz();
  hCln->Draw("colz");
  if (ab) {
    TF1 pol1("fpol1","pol1",hCln->GetXaxis()->GetXmin(),hCln->GetXaxis()->GetXmax());
    pol1.SetLineColor(kMagenta);
    pol1.SetParameter(0,ab[0]);
    pol1.SetParameter(1,ab[2]);
    pol1.Draw("same");
  }
  c->Print("SPDOnVsOf.png");
}


