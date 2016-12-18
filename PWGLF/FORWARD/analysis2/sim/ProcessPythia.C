/**
 * @file   ProcessPythia.C
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Fri Mar 20 12:16:33 2015
 * 
 * @brief  Example of processing output from FastSim.C 
 * 
 * 
 */

/** 
 * 
 * 
 * @param max 
 * @param type 
 * @param ppEG 
 * @param sNN 
 */
void Short(Long64_t    max,
	   Int_t      type,
	   const char* ppEG,
	   UShort_t    sNN)
{
  ULong_t runNo = 0;
  switch (sNN) {
  case 900:   runNo  = 118506; break;
  case 2760:  runNo  = 197669; break;
  case 7000:  runNo  = 126407; break;
  case 8000:  runNo  = 192708; break;
  default:
    Printf("No run number for sqrt{s}=%dGeV", sNN);
    return;
  }
  TString pp(Form("%s_%09lu_100k.root", ppEG, runNo));
  TString data(Form("ppExport/pp_%04hu_%s_coarse_empirical.root", sNN,
		    (type==0 ? "INEL" : "V0AND")));
  TString tit(Form("pp %s @ #sqrt{#it{s}}=",(type==0 ? "INEL" : "NSD")));
  if      (sNN < 1000)        tit.Append(Form("%dGeV", sNN));
  else if ((sNN % 1000) == 0) tit.Append(Form("%dTeV",sNN/1000));
  else                        tit.Append(Form("%.2fTeV",Float_t(sNN)/1000));

  ProcessFast(max, type, pp, data, tit);
}

void PythiaMult(Long64_t max=-1) {
  ProcessFast(max,2,"Pythia_000126407_PP_07000_100k.root",0,"Pythia_by_mult");
}
void Pythia0900(bool inel, Long64_t max=-1) { Short(max,inel,"Pythia", 900);}
void Pythia2760(bool inel, Long64_t max=-1) { Short(max,inel,"Pythia",2760);}
void Pythia7000(bool inel, Long64_t max=-1) { Short(max,inel,"Pythia",7000);}
void Pythia8000(bool inel, Long64_t max=-1) { Short(max,inel,"Pythia",8000);}

void Pythia()
{
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/FORWARD/analysis2/sim/ProcessFast.C+");
  Pythia0900(true);
  Pythia2760(true);
  Pythia7000(true);
  Pythia8000(true);
  Pythia0900(false);
  Pythia2760(false);
  Pythia7000(false);
  Pythia8000(false);
  
}

//
//  EOF
//  
	  
