#ifndef __CINT__
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#endif
#include "triggerInfo.C"

// TODO read number of bits from AliVEvent?
#define NBITS 29
#define NMAXCLASSES 100
TString bitNames[NBITS] = {
"kINT1",
"kINT7",
"kMUON",
"kHighMult",
"kEMC1",
"kCINT5",
"kCMUS5",
"kMUSH7",
"kMUL7",
"kMUU7",
"kEMC7",
"kMUS7",
"kPHI1",
"kPHI78",
"kEMCEJE",
"kEMCEGA",
"kCentral",
"kSemiCentral",
"kDG5",
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

Int_t runLevelEventStatQA(TString qafilename="event_stat.root", Int_t run=254394, TString ocdbStorage = "raw://"){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetPadLeftMargin(0.07);

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
  Double_t  class_lifetime[NMAXCLASSES]    = {0};
  Double_t  class_lumi[NMAXCLASSES]        = {0};
  Double_t  class_ds[NMAXCLASSES]          = {0};
  ULong64_t alias_recorded[NBITS]          = {0};
  ULong64_t alias_reconstructed[NBITS]     = {0};
  ULong64_t alias_accepted[NBITS]          = {0};
  Double_t alias_l0b_rate[NBITS]           = {0};
  Double_t alias_lifetime[NBITS]           = {0};
  Double_t alias_lumi_recorded[NBITS]      = {0};
  Double_t alias_lumi_reconstructed[NBITS] = {0};
  Double_t alias_lumi_accepted[NBITS]      = {0};
  Int_t timeStart = 0;
  Int_t timeEnd = 0;

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
  t->Branch("class_lifetime",&class_lifetime,Form("class_lifetime[%i]/D",NMAXCLASSES));
  t->Branch("class_lumi",&class_lumi,Form("class_lumi[%i]/D",NMAXCLASSES));
  t->Branch("class_ds",&class_ds,Form("class_ds[%i]/D",NMAXCLASSES));
  t->Branch("alias_recorded",&alias_recorded,Form("alias_recorded[%i]/l",NBITS));
  t->Branch("alias_reconstructed",&alias_reconstructed,Form("alias_reconstructed[%i]/l",NBITS));
  t->Branch("alias_accepted",&alias_accepted,Form("alias_accepted[%i]/l",NBITS));
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
  else if (run>=246995               ) { refSigma= 30.0; refEff = 0.40; refClass = "C0TVX-B-NOPF-CENTNOTRD"; } // estimates from Cvetan and Alberica
  Double_t orbitRate = 11245.;
  TString partition;
  TString lhcState;
  TString lhcPeriod;
  TString activeDetectorsString;
  //Int_t run, TString ocdbStorage, TString &partition, TString &activeDetectorsString, Double_t& run_duration, 
  classes = GetClasses(run,ocdbStorage,partition,activeDetectorsString,run_duration,class_lMb,class_lMa,class_l0b,class_l0a,class_l1b,class_l1a,class_l2b,class_l2a);
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

  if (run>=244917) {
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
  
  
  
  TFile* fin = new TFile(qafilename);
  TH2D* h = fin ? (TH2D*) fin->Get("fHistStat") : 0;
  if (!fin || !h) {
    printf("fHistStatistics not found\n");
    writeTree(fout,t);
    return 4; 
  }

  fout->cd();
  
  for (Int_t j=1;j<=h->GetNbinsY();j++){
    if (j!=3) continue;
    TString label = h->GetYaxis()->GetBinLabel(j);
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
    if (ibit>=NBITS) continue;
//    if (alias_recorded[ibit]) break; 

    alias_reconstructed[ibit] = Int_t(h->GetBinContent(1,j));
    alias_accepted[ibit]      = Int_t(h->GetBinContent(2,j));
    
//    printf("%4i %8i %8i\n",ibit,alias_reconstructed[ibit],alias_accepted[ibit]);
    
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
    TList* list = (TList*) fin->Get(Form("trigger_histograms_%s/histos",label.Data()));
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
    TH1F* hV0MOfAll        = (TH1F*) list->FindObject("fHistV0MOfAll");
    TH1F* hV0MOfAcc        = (TH1F*) list->FindObject("fHistV0MOfAcc");
    TH1F* hFiredBitsSPD    = (TH1F*) list->FindObject("fHistFiredBitsSPD");
    TH2F* hSPDClsVsTklAll  = (TH2F*) list->FindObject("fHistSPDClsVsTklAll");
    TH2F* hSPDClsVsTklCln  = (TH2F*) list->FindObject("fHistSPDClsVsTklCln");
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
    TH1F* hTKLAll          = (TH1F*) list->FindObject("fHistTKLAll");
    TH1F* hTKLAcc          = (TH1F*) list->FindObject("fHistTKLAcc");
    TH2F* hAD              = (TH2F*) list->FindObject("fHistAD");
    TH1F* hTimeZNA         = (TH1F*) list->FindObject("fHistTimeZNA");
    TH1F* hTimeZNC         = (TH1F*) list->FindObject("fHistTimeZNC");
    TH2F* hTimeZNSumVsDif  = (TH2F*) list->FindObject("fHistTimeZNSumVsDif");
    TH2F* hTimeCorrZDC     = (TH2F*) list->FindObject("fHistTimeCorrZDC");
    
    TH1F* hFitResults = (TH1F*) hV0MOnVsOfCln->ProjectionX("hFitResults");
    TF1* fExp = new TF1("fExp","exp(-([0]-x)*[1])",0,10000);
    TF1* fPol1 = new TF1("fPol1","pol1",0,hV0MOnVsOfCln->GetXaxis()->GetXmax());
    fExp->SetLineColor(kMagenta);
    hFitResults->Reset();
    if (0 && hV0MOnVsOfAll && hV0MOnVsOfCln) {
      TCanvas* cV0MOnVsOf = new TCanvas(Form("cV0MOnVsOf_%s",bitName),Form("cV0MOnVsOf_%s",bitName),1800,800);
      cV0MOnVsOf->Divide(2,1,0.001,0.001);
      cV0MOnVsOf->cd(1);
      gPad->SetLeftMargin(0.11);
      gPad->SetRightMargin(0.11);
      gPad->SetLogz();
      hV0MOnVsOfAll->SetTitle(Form("%s",bitName));
      hV0MOnVsOfAll->GetYaxis()->SetTitleOffset(1.7);
      hV0MOnVsOfAll->Draw("colz");
      cV0MOnVsOf->cd(2);
      gPad->SetLogz();
      gPad->SetLeftMargin(0.11);
      gPad->SetRightMargin(0.11);
      hV0MOnVsOfCln->SetTitle(Form("%s",bitName));
      hV0MOnVsOfCln->GetYaxis()->SetTitleOffset(1.7);
      hV0MOnVsOfCln->Draw("colz");
      cV0MOnVsOf->Print(Form("%s_V0MOnVsOf.pdf",bitName));
      hV0MOnVsOfCln->Sumw2();
      TCanvas* cV0MOnVsOfFit = new TCanvas(Form("cV0MOnVsOfFit_%s",bitName),Form("cV0MOnVsOfFit_%s",bitName),1800,800);
      cV0MOnVsOfFit->Divide(2,1,0.001,0.001);
      cV0MOnVsOfFit->cd(1);
      gPad->SetLogy();
      TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
      leg->SetBorderSize(0);
      for (Int_t i=1;i<=hV0MOnVsOfCln->GetNbinsX();i++){
        TH1D* hProj = hV0MOnVsOfCln->ProjectionY(Form("projV0MOnVsOf%i",i),i,i);
        hProj->GetXaxis()->SetTitleOffset(1.0);
        hProj->SetTitle(Form("%s;Online V0M;",bitName));
        Float_t integ = hProj->Integral();
        if (integ<200) continue;
        Double_t mMax = hProj->GetBinLowEdge(hProj->GetMaximumBin()-1);
        Double_t mMin = hProj->GetBinLowEdge(hProj->GetMaximumBin()-6);
        fExp->SetRange(mMin,mMax);
        hProj->Fit("fExp","QN","",mMin,mMax);
        if (i==20 || i==40 || i==60 || i==80 || i==100) {
          hProj->SetLineColor(i/20);
          hProj->DrawCopy(i==20 ? "" : "same");
          fExp->DrawCopy("lsame");

          leg->AddEntry(hProj,Form("Offline V0M = %.0f",hV0MOnVsOfCln->GetXaxis()->GetBinCenter(i)));
        }
        
        Double_t x0     = fExp->GetParameter(0);
        Double_t lambda = 1/fExp->GetParameter(1);
        Double_t x0err  = fExp->GetParError(0);
        Double_t lambdaerr = lambda*fExp->GetParError(1)/fExp->GetParameter(1);
        hFitResults->SetBinContent(i,x0-7*lambda);
        hFitResults->SetBinError(i,x0err+7*lambdaerr);
//        gPad->Print("fit.png");
      }
      leg->Draw();
      
      cV0MOnVsOfFit->cd(2);
      gPad->SetLogz();
      gPad->SetLeftMargin(0.11);
      gPad->SetRightMargin(0.11);
      hV0MOnVsOfCln->GetListOfFunctions()->Clear();
      hV0MOnVsOfCln->DrawCopy("colz");
      hFitResults->SetLineWidth(2);
      hFitResults->SetLineColor(kMagenta+1);
      hFitResults->SetMarkerColor(kMagenta+1);
      hFitResults->SetMarkerStyle(kOpenCircle);
      hFitResults->DrawCopy("psame");
      hFitResults->Fit(fPol1,"QN","",50,hV0MOnVsOfCln->GetXaxis()->GetXmax());
      printf("par[0]=%f\n",fPol1->GetParameter(0));
      printf("par[1]=%f\n",fPol1->GetParameter(1));
      fPol1->Draw("same");
      cV0MOnVsOfFit->Print(Form("%s_V0MOnVsOfFit.pdf",bitName));
      
    } else printf("QA histogram not found\n");
    
    if (hSPDClsVsTklAll && hSPDClsVsTklCln) {
      TCanvas* cClsVsTkl = new TCanvas(Form("cClsVsTkl_%s",bitName),Form("cClsVsTkl_%s",bitName),1800,800);
      cClsVsTkl->Divide(2,1,0.001,0.001);
      cClsVsTkl->cd(1);
      gPad->SetLogz();
      hSPDClsVsTklAll->SetTitle(Form("%s",bitName));
      hSPDClsVsTklAll->Draw("colz");
      cClsVsTkl->cd(2);
      gPad->SetLogz();
      hSPDClsVsTklCln->SetTitle(Form("%s",bitName));
      hSPDClsVsTklCln->Draw("colz");
      cClsVsTkl->Print(Form("%s_ClsVsTkl.pdf",bitName));
    } else printf("QA histogram not found\n");
    
    if (hSPDOnVsOfAll && hSPDOnVsOfCln) {
      TCanvas* cSPDOnVsOf = new TCanvas(Form("cSPDOnVsOf_%s",bitName),Form("cSPDOnVsOf_%s",bitName),1800,800);
      cSPDOnVsOf->Divide(2,1,0.001,0.001);
      cSPDOnVsOf->cd(1);
      gPad->SetLogz();
      hSPDOnVsOfAll->SetTitle(Form("%s",bitName));
      hSPDOnVsOfAll->Draw("colz");
      cSPDOnVsOf->cd(2);
      gPad->SetLogz();
      hSPDOnVsOfCln->SetTitle(Form("%s",bitName));
      hSPDOnVsOfCln->Draw("colz");
      cSPDOnVsOf->Print(Form("%s_SPDOnVsOf.pdf",bitName));
    } else printf("QA histogram not found\n");

    
    if (hBBAflagsAll && hBBAflagsAcc && hBBCflagsAll && hBBCflagsAcc && hBGAflagsAll && hBGAflagsAcc && hBGCflagsAll && hBGCflagsAcc) {
      TCanvas* cV0flags = new TCanvas(Form("V0flags_%s",bitName),Form("cV0flags_%s",bitName),1000,900);
      cV0flags->Divide(2,2,0.001,0.001);
      cV0flags->cd(1);
      gPad->SetLogy();
      hBBAflagsAll->SetTitle(Form("%s",bitName));
      hBBAflagsAll->SetMinimum(0.9);
      hBBAflagsAll->SetLineColor(kBlack);
      hBBAflagsAcc->SetLineColor(kBlue);
      hBBAflagsAcc->SetFillColor(kBlue);
      hBBAflagsAll->Draw();
      hBBAflagsAcc->Draw("same");
      cV0flags->cd(2);
      gPad->SetLogy();
      hBBCflagsAll->SetTitle(Form("%s",bitName));
      hBBCflagsAll->SetMinimum(0.9);
      hBBCflagsAll->SetLineColor(kBlack);
      hBBCflagsAcc->SetLineColor(kBlue);
      hBBCflagsAcc->SetFillColor(kBlue);
      hBBCflagsAll->Draw();
      hBBCflagsAcc->Draw("same");
      cV0flags->cd(3);
      gPad->SetLogy();
      hBGAflagsAll->SetTitle(Form("%s",bitName));
      hBGAflagsAll->SetMinimum(0.9);
      hBGAflagsAll->SetLineColor(kBlack);
      hBGAflagsAcc->SetLineColor(kBlue);
      hBGAflagsAcc->SetFillColor(kBlue);
      hBGAflagsAll->Draw();
      hBGAflagsAcc->Draw("same");
      cV0flags->cd(4);
      gPad->SetLogy();
      hBGCflagsAll->SetTitle(Form("%s",bitName));
      hBGCflagsAll->SetMinimum(0.9);
      hBGCflagsAll->SetLineColor(kBlack);
      hBGCflagsAcc->SetLineColor(kBlue);
      hBGCflagsAcc->SetFillColor(kBlue);
      hBGCflagsAll->Draw();
      hBGCflagsAcc->Draw("same");
      cV0flags->Print(Form("%s_V0flags.pdf",bitName));
    }

    if (hV0AAll && hV0AAcc && hV0CAll && hV0CAcc) {
      TCanvas* cV0 = new TCanvas(Form("V0_%s",bitName),Form("cV0_%s",bitName),1800,800);
      cV0->Divide(2,1,0.001,0.001);
      cV0->cd(1);
      gPad->SetLogy();
      hV0AAll->SetTitle(Form("%s",bitName));
      hV0AAll->SetLineWidth(2);
      hV0AAll->SetLineColor(kBlack);
      hV0AAll->Draw();
      hV0AAcc->SetLineColor(kBlue);
      hV0AAcc->SetFillColor(kBlue);
      hV0AAcc->Draw("same");
      TLegend* leg = new TLegend(0.15,0.8,0.45,0.9);
      leg->SetBorderSize(0);
      leg->AddEntry(hV0AAll,"All events");
      leg->AddEntry(hV0AAcc,"Accepted events");
      leg->Draw();
      cV0->cd(2);
      gPad->SetLogy();
      hV0CAll->SetTitle(Form("%s",bitName));
      hV0CAll->SetLineWidth(2);
      hV0CAll->SetLineColor(kBlack);
      hV0CAll->Draw();
      hV0CAcc->SetLineColor(kBlue);
      hV0CAcc->SetFillColor(kBlue);
      hV0CAcc->Draw("same");
      leg->Draw();
      cV0->Print(Form("%s_V0.pdf",bitName));
    } else printf("QA histogram not found\n");

    if (hADAAll && hADAAcc && hADCAll && hADCAcc) {
      TCanvas* cAD = new TCanvas(Form("AD_%s",bitName),Form("cAD_%s",bitName),1800,800);
      cAD->Divide(2,1,0.001,0.001);
      cAD->cd(1);
      gPad->SetLogy();
      hADAAll->SetTitle(Form("%s",bitName));
      hADAAll->SetLineWidth(2);
      hADAAll->SetLineColor(kBlack);
      hADAAll->Draw();
      hADAAcc->SetLineColor(kBlue);
      hADAAcc->SetFillColor(kBlue);
      hADAAcc->Draw("same");
      TLegend* leg = new TLegend(0.15,0.8,0.45,0.9);
      leg->SetBorderSize(0);
      leg->AddEntry(hADAAll,"All events");
      leg->AddEntry(hADAAcc,"Accepted events");
      leg->Draw();
      cAD->cd(2);
      gPad->SetLogy();
      hADCAll->SetTitle(Form("%s",bitName));
      hADCAll->SetLineWidth(2);
      hADCAll->SetLineColor(kBlack);
      hADCAll->Draw();
      hADCAcc->SetLineColor(kBlue);
      hADCAcc->SetFillColor(kBlue);
      hADCAcc->Draw("same");
      leg->Draw();
      cAD->Print(Form("%s_AD.pdf",bitName));
    } else printf("QA histogram not found\n"); 

    TH1F* hV0MOnNorm = new TH1F("hV0MOnNorm",";Online V0M/<V0M>;",50,0,10);
    TH1F* hV0MOfNorm = new TH1F("hV0MOfNorm",";Offline V0M/<V0M>;",50,0,10);
    hV0MOnNorm->SetLineWidth(2);
    hV0MOfNorm->SetLineWidth(2);
    hV0MOnNorm->SetLineColor(kBlue);
    hV0MOfNorm->SetLineWidth(kBlue);

    if (hV0MOnAll && hV0MOfAll && hV0MOnAcc && hV0MOfAcc) {
      TCanvas* cV0M = new TCanvas(Form("V0M_%s",bitName),Form("cV0M_%s",bitName),1800,1100);
      cV0M->Divide(2,2,0.001,0.001);
      
      cV0M->cd(1);
      gPad->SetLogy();
      hV0MOnAll->SetTitle(Form("%s",bitName));
      hV0MOnAll->SetLineWidth(2);
      hV0MOnAll->SetLineColor(kBlack);
      hV0MOnAll->Draw();
      hV0MOnAcc->SetLineWidth(2);
      hV0MOnAcc->SetMarkerColor(kBlue);
      hV0MOnAcc->SetLineColor(kBlue);
      hV0MOnAcc->SetFillColor(kBlue);
      hV0MOnAcc->Draw("same");
      TLegend* leg1 = new TLegend(0.40,0.8,0.97,0.9);
      leg1->SetBorderSize(0);
      leg1->AddEntry(hV0MOnAll,Form("All events: %.0f",hV0MOnAll->Integral(0,hV0MOnAll->GetNbinsX()+1)));
      leg1->AddEntry(hV0MOnAcc,Form("Accepted events: %.0f, <V0M>=%.1f",hV0MOnAcc->Integral(0,hV0MOnAcc->GetNbinsX()+1),hV0MOnAcc->GetMean()),"l");
      leg1->Draw();
      
      cV0M->cd(2);
      gPad->SetLogy();
      hV0MOfAll->SetTitle(Form("%s",bitName));
      hV0MOfAll->SetLineWidth(2);
      hV0MOfAll->SetLineColor(kBlack);
      hV0MOfAll->Draw();
      hV0MOfAcc->SetLineWidth(2);
      hV0MOfAcc->SetMarkerColor(kBlue);
      hV0MOfAcc->SetLineColor(kBlue);
      hV0MOfAcc->SetFillColor(kBlue);
      hV0MOfAcc->Draw("same");
      TLegend* leg2 = new TLegend(0.40,0.8,0.97,0.9);
      leg2->SetBorderSize(0);
      leg2->AddEntry(hV0MOfAll,Form("All events: %.0f",hV0MOfAll->Integral(0,hV0MOfAll->GetNbinsX()+1)));
      leg2->AddEntry(hV0MOfAcc,Form("Accepted events: %.0f, <V0M>=%.1f",hV0MOfAcc->Integral(0,hV0MOfAcc->GetNbinsX()+1),hV0MOfAcc->GetMean()),"l");
      leg2->Draw();
      
      cV0M->cd(3);
      gPad->SetLogy();
      hV0MOnNorm->SetTitle(Form("%s",bitName));
      Float_t meanV0MOn = hV0MOnAcc->GetMean();
      for (Int_t i=1;i<=hV0MOnAcc->GetNbinsX();i++){
        Float_t binCenter  = hV0MOnAcc->GetXaxis()->GetBinCenter(i);
        Float_t binContent = hV0MOnAcc->GetBinContent(i);
        hV0MOnNorm->Fill(binCenter/meanV0MOn,binContent);
      }
      hV0MOnNorm->DrawCopy();
      cV0M->cd(4);
      gPad->SetLogy();
      hV0MOfNorm->SetTitle(Form("%s",bitName));
      Float_t meanV0MOf = hV0MOfAcc->GetMean();
      for (Int_t i=1;i<=hV0MOfAcc->GetNbinsX();i++){
        Float_t binCenter  = hV0MOfAcc->GetXaxis()->GetBinCenter(i);
        Float_t binContent = hV0MOfAcc->GetBinContent(i);
        hV0MOfNorm->Fill(binCenter/meanV0MOf,binContent);
      }
      hV0MOfNorm->DrawCopy();
      
      cV0M->Print(Form("%s_V0M.pdf(",bitName));
      
      TCanvas* cV0Mcum = new TCanvas(Form("V0Mcum_%s",bitName),Form("cV0Mcum_%s",bitName),1800,1100);
      cV0Mcum->Divide(2,2,0.001,0.001);
      cV0Mcum->cd(1);
      gPad->SetLogy();
      TH1F* hV0MOnAllCum = (TH1F*) hV0MOnAll->GetCumulative(kFALSE);
      TH1F* hV0MOnAccCum = (TH1F*) hV0MOnAcc->GetCumulative(kFALSE);
      hV0MOnAllCum->Scale(1./hV0MOnAllCum->GetBinContent(1));
      hV0MOnAccCum->Scale(1./hV0MOnAccCum->GetBinContent(1));
      hV0MOnAllCum->SetTitle(Form("%s: cumulative",bitName));
      hV0MOnAllCum->DrawCopy();
      hV0MOnAccCum->SetFillStyle(0);
      hV0MOnAccCum->DrawCopy("same");
      leg1->Draw();
      cV0Mcum->cd(2);
      gPad->SetLogy();
      TH1F* hV0MOfAllCum = hV0MOfAll->GetCumulative(kFALSE);
      TH1F* hV0MOfAccCum = hV0MOfAcc->GetCumulative(kFALSE);
      hV0MOfAllCum->Scale(1./hV0MOfAllCum->GetBinContent(1));
      hV0MOfAccCum->Scale(1./hV0MOfAccCum->GetBinContent(1));
      hV0MOfAllCum->SetTitle(Form("%s: cumulative",bitName));
      hV0MOfAllCum->DrawCopy();
      hV0MOfAccCum->SetFillStyle(0);
      hV0MOfAccCum->DrawCopy("same");
      leg2->Draw();
      cV0Mcum->cd(3);
      gPad->SetLogy();
      TH1F* hV0MOnNormCum = hV0MOnNorm->GetCumulative(kFALSE);
      hV0MOnNormCum->Scale(1./hV0MOnNormCum->GetBinContent(1));
      hV0MOnNormCum->SetTitle(Form("%s: cumulative",bitName));
      hV0MOnNormCum->DrawCopy();
      cV0Mcum->cd(4);
      gPad->SetLogy();
      TH1F* hV0MOfNormCum = hV0MOfNorm->GetCumulative(kFALSE);
      hV0MOfNormCum->Scale(1./hV0MOfNormCum->GetBinContent(1));
      hV0MOfNormCum->SetTitle(Form("%s: cumulative",bitName));
      hV0MOfNormCum->DrawCopy();
      cV0Mcum->Print(Form("%s_V0M.pdf)",bitName));
    } else printf("QA histogram not found\n");
    
    TH1F* hOFONorm = new TH1F("hOFONorm",";Outer FO/<Outer FO>;",50,0,10);
    TH1F* hTKLNorm = new TH1F("hTKLNorm",";TKL/<TKL>;",50,0,10);
    hOFONorm->SetLineWidth(2);
    hTKLNorm->SetLineWidth(2);
    hOFONorm->SetLineColor(kBlue);
    hTKLNorm->SetLineWidth(kBlue);

    if (hOFOAll && hTKLAll && hOFOAcc && hTKLAcc) {
      TCanvas* cSPD = new TCanvas(Form("SPD_%s",bitName),Form("cSPD_%s",bitName),1800,1100);
      cSPD->Divide(2,2,0.001,0.001);
      
      cSPD->cd(1);
      gPad->SetLogy();
      hOFOAll->SetTitle(Form("%s",bitName));
      hOFOAll->SetLineWidth(2);
      hOFOAll->SetLineColor(kBlack);
      hOFOAll->Draw();
      hOFOAcc->SetLineWidth(2);
      hOFOAcc->SetMarkerColor(kBlue);
      hOFOAcc->SetLineColor(kBlue);
      hOFOAcc->SetFillColor(kBlue);
      hOFOAcc->Draw("same");
      TLegend* leg1 = new TLegend(0.40,0.8,0.97,0.9);
      leg1->SetBorderSize(0);
      leg1->AddEntry(hOFOAll,Form("All events: %.0f",hOFOAll->Integral(0,hOFOAll->GetNbinsX()+1)));
      leg1->AddEntry(hOFOAcc,Form("Accepted events: %.0f, <OFO>=%.1f",hOFOAcc->Integral(0,hOFOAcc->GetNbinsX()+1),hOFOAcc->GetMean()),"l");
      leg1->Draw();
      
      cSPD->cd(2);
      gPad->SetLogy();
      hTKLAll->SetTitle(Form("%s",bitName));
      hTKLAll->SetLineWidth(2);
      hTKLAll->SetLineColor(kBlack);
      hTKLAll->Draw();
      hTKLAcc->SetLineWidth(2);
      hTKLAcc->SetMarkerColor(kBlue);
      hTKLAcc->SetLineColor(kBlue);
      hTKLAcc->SetFillColor(kBlue);
      hTKLAcc->Draw("same");
      TLegend* leg2 = new TLegend(0.40,0.8,0.97,0.9);
      leg2->SetBorderSize(0);
      leg2->AddEntry(hTKLAll,Form("All events: %.0f",hTKLAll->Integral(0,hTKLAll->GetNbinsX()+1)));
      leg2->AddEntry(hTKLAcc,Form("Accepted events: %.0f, <TKL>=%.1f",hTKLAcc->Integral(0,hTKLAcc->GetNbinsX()+1),hTKLAcc->GetMean()),"l");
      leg2->Draw();
      
      cSPD->cd(3);
      gPad->SetLogy();
      hOFONorm->SetTitle(Form("%s",bitName));
      Float_t meanOFO = hOFOAcc->GetMean();
      for (Int_t i=1;i<=hOFOAcc->GetNbinsX();i++){
        Float_t binCenter  = hOFOAcc->GetXaxis()->GetBinCenter(i);
        Float_t binContent = hOFOAcc->GetBinContent(i);
        hOFONorm->Fill(binCenter/meanOFO,binContent);
      }
      hOFONorm->DrawCopy();
      cSPD->cd(4);
      gPad->SetLogy();
      hTKLNorm->SetTitle(Form("%s",bitName));
      Float_t meanTKL = hTKLAcc->GetMean();
      for (Int_t i=1;i<=hTKLAcc->GetNbinsX();i++){
        Float_t binCenter  = hTKLAcc->GetXaxis()->GetBinCenter(i);
        Float_t binContent = hTKLAcc->GetBinContent(i);
        hTKLNorm->Fill(binCenter/meanTKL,binContent);
      }
      hTKLNorm->DrawCopy();
      
      cSPD->Print(Form("%s_SPD.pdf(",bitName));
      
      TCanvas* cSPDcum = new TCanvas(Form("SPDcum_%s",bitName),Form("cSPDcum_%s",bitName),1800,1100);
      cSPDcum->Divide(2,2,0.001,0.001);
      cSPDcum->cd(1);
      gPad->SetLogy();
      TH1F* hOFOAllCum = hOFOAll->GetCumulative(kFALSE);
      TH1F* hOFOAccCum = hOFOAcc->GetCumulative(kFALSE);
      hOFOAllCum->Scale(1./hOFOAllCum->GetBinContent(1));
      hOFOAccCum->Scale(1./hOFOAccCum->GetBinContent(1));
      hOFOAllCum->SetTitle(Form("%s: cumulative",bitName));
      hOFOAllCum->DrawCopy();
      hOFOAccCum->SetFillStyle(0);
      hOFOAccCum->DrawCopy("same");
      leg1->Draw();
      cSPDcum->cd(2);
      gPad->SetLogy();
      TH1F* hTKLAllCum = hTKLAll->GetCumulative(kFALSE);
      TH1F* hTKLAccCum = hTKLAcc->GetCumulative(kFALSE);
      hTKLAllCum->Scale(1./hTKLAllCum->GetBinContent(1));
      hTKLAccCum->Scale(1./hTKLAccCum->GetBinContent(1));
      hTKLAllCum->SetTitle(Form("%s: cumulative",bitName));
      hTKLAllCum->DrawCopy();
      hTKLAccCum->SetFillStyle(0);
      hTKLAccCum->DrawCopy("same");
      leg2->Draw();
      cSPDcum->cd(3);
      gPad->SetLogy();
      TH1F* hOFONormCum = hOFONorm->GetCumulative(kFALSE);
      hOFONormCum->Scale(1./hOFONormCum->GetBinContent(1));
      hOFONormCum->SetTitle(Form("%s: cumulative",bitName));
      hOFONormCum->DrawCopy();
      cSPDcum->cd(4);
      gPad->SetLogy();
      TH1F* hTKLNormCum = hTKLNorm->GetCumulative(kFALSE);
      hTKLNormCum->Scale(1./hTKLNormCum->GetBinContent(1));
      hTKLNormCum->SetTitle(Form("%s: cumulative",bitName));
      hTKLNormCum->DrawCopy();
      cSPDcum->Print(Form("%s_SPD.pdf)",bitName));
    } else printf("QA histogram not found\n");

    if (hTimeZNA && hTimeZNC && hTimeZNSumVsDif && hTimeCorrZDC) {
      TCanvas* cZDC = new TCanvas(Form("cZDC_%s",bitName),Form("cZDC_%s",bitName),1800,800);
      cZDC->Divide(2,2,0.001,0.001);
      cZDC->cd(1);
      gPad->SetLogy();
      hTimeZNA->SetTitle(Form("%s",bitName));
      hTimeZNA->Draw();
      cZDC->cd(2);
      gPad->SetLogy();
      hTimeZNC->SetTitle(Form("%s",bitName));
      hTimeZNC->Draw();
      cZDC->Print(Form("%s_ZDC.pdf",bitName));
      cZDC->cd(3);
      gPad->SetLogz();
      hTimeCorrZDC->SetTitle(Form("%s",bitName));
      hTimeCorrZDC->Draw("colz");
      cZDC->cd(4);
      gPad->SetLogz();
      hTimeZNSumVsDif->SetTitle(Form("%s",bitName));
      hTimeZNSumVsDif->Draw("colz");
      cZDC->Print(Form("%s_ZDC.pdf",bitName));
    } else printf("QA histogram not found\n");
    
    //    TH2F* hAD           = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistAD"          ,label.Data()));
    //    TH1F* hADA          = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistADA"         ,label.Data()));
    //    TH1F* hADC          = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistADC"         ,label.Data()));
    //    TH1F* hV0A          = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistV0A"         ,label.Data()));
    //    TH1F* hV0C          = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistV0C"         ,label.Data()));
    //    TH1F* hFiredBitsSPD = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistFiredBitsSPD",label.Data()));
    //    TH2F* hBitsSPD      = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistBitsSPD"     ,label.Data()));
    //    TH1F* hTDCZDC       = (TH1F*) fin->Get(Form("trigger_histograms_%s/fHistTDCZDC"      ,label.Data()));
    //    TH2F* hTimeZDC      = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistTimeZDC"     ,label.Data()));
    //    TH2F* hTimeCorrZDC  = (TH2F*) fin->Get(Form("trigger_histograms_%s/fHistTimeCorrZDC" ,label.Data()));
    
//    if (hAD) {
//      TCanvas* cAD = new TCanvas(Form("cAD_%s",bitName),Form("cAD_%s",bitName),800,800);
//      gPad->SetLogz();
//      gPad->SetMargin(0.12,0.12,0.10,0.06);
//      hAD->SetTitle(Form("%s: AD timing;AD timing C-A;AD timing C+A",bitName));
//      hAD->GetXaxis()->SetTitleOffset(1.3);
//      hAD->GetYaxis()->SetTitleOffset(1.6);
//      hAD->Draw("colz");
//      gPad->Print(Form("%s_AD.pdf",bitName));
//      hAD->Write(Form("%s_AD",bitName));
//    } else printf("QA histogram not found\n"); 
//    
//    if (hADA) {
//      TCanvas* cADA = new TCanvas(Form("cADA_%s",bitName),Form("cADA_%s",bitName),1000,800);
//      gPad->SetLogy();
//      hADA->SetTitle(Form("%s: ADA",bitName));
//      hADA->SetLineWidth(2);
//      hADA->SetLineColor(kBlue);
//      hADA->Draw();
//      gPad->Print(Form("%s_ADA.pdf",bitName));
//      hADA->Write(Form("%s_ADA",bitName));
//    } else printf("QA histogram not found\n"); 
//
//    if (hADC) {
//      TCanvas* cADC = new TCanvas(Form("cADC_%s",bitName),Form("cADC_%s",bitName),1000,800);
//      gPad->SetLogy();
//      hADC->SetTitle(Form("%s: ADC",bitName));
//      hADC->SetLineWidth(2);
//      hADC->SetLineColor(kBlue);
//      hADC->Draw();
//      gPad->Print(Form("%s_ADC.pdf",bitName));
//      hADC->Write(Form("%s_ADC",bitName));
//    } else printf("QA histogram not found\n"); 
//
//    if (hV0A) {
//      TCanvas* cV0A = new TCanvas(Form("cV0A_%s",bitName),Form("cV0A_%s",bitName),1000,800);
//      gPad->SetLogy();
//      hV0A->SetTitle(Form("%s: V0A",bitName));
//      hV0A->SetLineWidth(2);
//      hV0A->SetLineColor(kBlue);
//      hV0A->Draw();
//      gPad->Print(Form("%s_V0A.pdf",bitName));
//      hV0A->Write(Form("%s_V0A",bitName));
//    } else printf("QA histogram not found\n"); 
//    
//    if (hV0C) {
//      TCanvas* cV0C = new TCanvas(Form("cV0C_%s",bitName),Form("cV0C_%s",bitName),1000,800);
//      gPad->SetLogy();
//      hV0C->SetTitle(Form("%s: V0C",bitName));
//      hV0C->SetLineWidth(2);
//      hV0C->SetLineColor(kBlue);
//      hV0C->Draw();
//      gPad->Print(Form("%s_V0C.pdf",bitName));
//      hV0C->Write(Form("%s_V0C",bitName));
//    } else printf("QA histogram not found\n"); 
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
//      gPad->Print(Form("%s_FiredBitsSPD.pdf",bitName));
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
//      gPad->Print(Form("%s_BitsSPD.pdf",bitName));
//      hBitsSPD->Write(Form("%s_BitsSPD",bitName));
//    } else printf("QA histogram not found\n"); 
//    
//    if (hTimeZDC) {
//      TCanvas* cTimeZDC = new TCanvas(Form("cTimeZDC_%s",bitName),Form("cTimeZDC_%s",bitName),800,800);
//      gPad->SetLogz();
//      gPad->SetMargin(0.12,0.12,0.10,0.06);
//      hTimeZDC->SetTitle(Form("%s: ZDC timing;TDC timing C-A;TDC timing C+A",bitName));
//      hTimeZDC->GetXaxis()->SetTitleOffset(1.3);
//      hTimeZDC->GetYaxis()->SetTitleOffset(1.6);
//      hTimeZDC->Draw("colz");
//      gPad->Print(Form("%s_TimeZDC.pdf",bitName));
//      hTimeZDC->Write(Form("%s_TimeZDC",bitName));
//    } else printf("QA histogram not found\n"); 
//
//    if (hTimeCorrZDC) {
//      TCanvas* cTimeCorrZDC = new TCanvas(Form("cTimeCorrZDC_%s",bitName),Form("cTimeCorrZDC_%s",bitName),800,800);
//      gPad->SetLogz();
//      gPad->SetMargin(0.12,0.12,0.10,0.06);
//      hTimeCorrZDC->SetTitle(Form("%s: corrected ZDC timing;TDC timing C-A;TDC timing C+A",bitName));
//      hTimeCorrZDC->GetXaxis()->SetTitleOffset(1.3);
//      hTimeCorrZDC->GetYaxis()->SetTitleOffset(1.6);
//      hTimeCorrZDC->Draw("colz");
//      gPad->Print(Form("%s_TimeCorrZDC.pdf",bitName));
//      hTimeCorrZDC->Write(Form("%s_TimeCorrZDC",bitName));
//    } else printf("QA histogram not found\n"); 
    gSystem->cd("..");
  }
  
  t->Fill();
  t->Write();
  fout->Close();
  return 0;
}


