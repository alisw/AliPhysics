#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjString.h>

#include "runCorrelationsStudiesConfigMacro.H"

void load2010bTestRunNumber();
void load2010cTestRunNumber();
void load2010bRunNumbers();
void loadLocal2010hMCTestRunNumber();
void load2010hTestRunNumber();
void load2010hAODTestRunNumber();
void load2010hRunNumbers();
void load2010hBMMRunNumbers();
void load2010hBPPRunNumbers();
void load2010hMCRunNumbers();
void load2010hMCTestRunNumber();
void loadMC2010h11a10a_bisBMMRunNumbers();
void loadMC2010h11a10a_bisBPPRunNumbers();
void loadMC2010h11a10a_bisAODRunNumbers(Int_t);
void loadMCAMPT2010h12A11RunNumbers();
void loadMCAMPT2010hCentrality(const char *cent);
void loadMC2013b2_efixRunNumbers();
void load2015oMCTestRunNumber();
void load2015oLIRpass2RunNumbers();
void load2015oLIRpass3RunNumbers();
void load2015oLIRpass4RunNumbers();
void load2015oHIRRunNumbers();
void load2015oHIRTestRunNumber();
void load2016kTestRunNumber();
void load2013bTestRunNumber();
void load2013bAODRunNumbers();
void load2013cAODRunNumbers();
void load2013bAODTestRunNumber();
void loadAMPT2015oHIRTestRunNumber();
void load2017nRunNumbers();
void load2017nMCRunNumbers();
void loadAMPT2760RunNumbers();
void load2018qRunNumbers();
void load2018rRunNumbers();
void load2018rMCRunNumbers();
void load2018rESDRunNumbers();

void runCorrelationsStudiesConfigMacro() {

  bTrainScope                = kFALSE;
  bGRIDPlugin                = kTRUE && !bTrainScope;
  bMC                        = kFALSE;
  bMConlyTruth               = kFALSE;

  /* tasks to use */
  /* only valid if not under train scope */
  bUsePhysicsSelection       = kTRUE;
  bUsePIDResponse            = kTRUE;
  // Centrality
  bUseCentralityTask         = kTRUE;
  bUseMultiplicityTask       = kTRUE;

  /* data input */
  isHeavyIon                 = kTRUE;

  /* reconstruction pass default */
  szpass = "2";

  /* Running conditions */
  szAliPhysicsVersion = "vAN-20181128_ROOT6-1";

  /* the number of files we want to test */
  nNoOfInputFiles = 30;
  nNoOfTestFiles = 4;

  /* load the run numbers */
  // load2010bTestRunNumber();
  // loadLocal2010hMCTestRunNumber();
  // load2010hTestRunNumber();
  // load2010hMCTestRunNumber();
  // load2010hAODTestRunNumber();
  // loadMC2010h11a10a_bisAODRunNumbers(-1);
  // load2015oMCTestRunNumber();
  // load2015oHIRTestRunNumber();
  // load2015oLIRpass3RunNumbers();
  // load2010hBPPRunNumbers();
  // load2010hMCRunNumbers();
  // load2016kTestRunNumber();
  // load2013bAODTestRunNumber();
  // load2013bAODRunNumbers();
  // loadMC2013b2_efixRunNumbers();
  // loadMCAMPT2010hCentrality("0-5");
  // loadAMPT2015oHIRTestRunNumber();
  // load2017nRunNumbers();
  // loadAMPT2760RunNumbers();
  // load2018qRunNumbers();
  load2018rESDRunNumbers();

  szRunPrefix = bMC ? "" : "000";

  /* local file list */
  szLocalFileList = (bMC? "filelist_mc.txt" : "filelist.txt");
}

void loadLocal2010hMCTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  bTrainScope                = kFALSE;
  bGRIDPlugin                = kFALSE && !bTrainScope;
  bMC                        = kTRUE;
}



/* 2010h period */
static const int nNoOf2010hBmmRuns = 45;
static const char *bmm2010hRunNumbers[nNoOf2010hBmmRuns] = {
    "138275", "138225", "138201", "138197", "138192", "138190", "137848", "137844", "137752", "137751",
    "137724", "137722", "137718", "137704", "137693", "137692", "137691", "137686", "137685", "137639",
    "137638", "137608", "137595", "137549", "137546", "137544", "137541", "137539", "137531", "137530",
    "137443", "137441", "137440", "137439", "137434", "137432", "137431", "137366", "137243", "137236",
    "137235", "137232", "137231", "137162", "137161"
};

static const int nNoOf2010hBppRuns = 45;
static const char *bpp2010hRunNumbers[nNoOf2010hBppRuns] = {
    "139510", "139507", "139505", "139503", "139465", "139438", "139437", "139360", "139329", "139328",
    "139314", "139310", "139309", "139173", "139107", "139105", "139038", "139037", "139036", "139029",
    "139028", "138872", "138871", "138870", "138837", "138732", "138730", "138666", "138662", "138653",
    "138652", "138638", "138624", "138621", "138583", "138582", "138579", "138578", "138534", "138469",
    "138442", "138439", "138438", "138396", "138364"
};

void load2010bTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010bTest";

  /* 2010b */
  szDataDir = "/alice/data/2010/LHC10b";

  /* 2010b */
  szDataPattern = "*pass4/*/AliESDs.root";

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("114918"));
}

void load2010cTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010cTest";

  /* 2010b */
  szDataDir = "/alice/data/2010/LHC10c";

  /* 2010b */
  szDataPattern = "*pass4/*/AliESDs.root";

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("120616"));
}

void load2010bRunNumbers() {
  const int nRuns = 45;
  const char *runnumber[nRuns] = {
      "114786","114798","114918","114920","114924","114930","114931","115186","115193","115310",
      "115318","115322","115328","115335","115345","115393","115399","115401","115414","115521",
      "116079","116081","116102","116288","116402","116403","116562","116571","116574","116643",
      "116645","117048","117050","117052","117053","117059","117060","117063","117092","117099",
      "117109","117112","117116","117220","117222"
  };

  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010bTest";

  /* 2010b */
  szDataDir = "/alice/data/2010/LHC10b";

  /* 2010b */
  szDataPattern = "*pass4/*/AliESDs.root";

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nRuns; run++)
    listOfActiveRuns.Add(new TObjString(runnumber[run]));
}




void load2010hTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hTest";

  /* 2010h */
  szDataDir = "/alice/data/2010/LHC10h";

  /* 2010h */
  szDataPattern = "*ESDs/pass2/*/AliESDs.root";

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("138534"));
}

void load2010hAODTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hTest";

  /* 2010h */
  szDataDir = "/alice/data/2010/LHC10h";

  /* 2010h */
  szDataPattern = "*ESDs/pass2/AOD160/*/AliAOD.root";

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("138534"));
}

void load2010hMCTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hMCTest";

  szDataDir = "/alice/sim/LHC11a10b_bis";
  szDataPattern = "*/AliESDs.root";

  bMC = kTRUE;

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("138534"));
}

void load2010hRunNumbers() {

  load2010hBMMRunNumbers();
  load2010hBPPRunNumbers();

  /* modify the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010h";
}

void load2010hBMMRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hBMM";

  /* 2010h */
  szDataDir = "/alice/data/2010/LHC10h";

  /* 2010h */
  szDataPattern = "*ESDs/pass2/*/AliESDs.root";

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nNoOf2010hBmmRuns; run++)
    listOfActiveRuns.Add(new TObjString(bmm2010hRunNumbers[run]));
}
void load2010hBPPRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hBPP";

  /* 2010h */
  szDataDir = "/alice/data/2010/LHC10h";

  /* 2010h */
  szDataPattern = "*ESDs/pass2/*/AliESDs.root";

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nNoOf2010hBppRuns; run++)
    listOfActiveRuns.Add(new TObjString(bpp2010hRunNumbers[run]));
}

void load2010hMCRunNumbers() {

  loadMC2010h11a10a_bisBMMRunNumbers();
  loadMC2010h11a10a_bisBPPRunNumbers();

  /* modify the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hMC";
}

void loadMC2010h11a10a_bisAODRunNumbers(Int_t step) {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hAODMC";

  szDataDir = "/alice/sim/LHC11a10a_bis";
  szDataPattern = "*/AOD162/*/AliAOD.root";

  bMC = kTRUE;

#define NOOFRUNSPERSTEP 22

  /* the list of runs to analyze */
  Int_t nNoOfRuns;
  Int_t nFirstRun;
  const char **runlist;
  switch (step) {
  case -1:
    nNoOfRuns = nNoOf2010hBmmRuns;
    nFirstRun = 0;
    runlist = bmm2010hRunNumbers;
    break;
  case -2:
    nNoOfRuns = nNoOf2010hBppRuns;
    nFirstRun = 0;
    runlist = bpp2010hRunNumbers;
    break;
  case 0:
    nNoOfRuns = NOOFRUNSPERSTEP;
    nFirstRun = 0;
    runlist = bmm2010hRunNumbers;
    break;
  case 1:
    nNoOfRuns = nNoOf2010hBmmRuns;
    nFirstRun = NOOFRUNSPERSTEP;
    runlist = bmm2010hRunNumbers;
    break;
  case 2:
    nNoOfRuns = NOOFRUNSPERSTEP;
    nFirstRun = 0;
    runlist = bpp2010hRunNumbers;
    break;
  case 3:
    nNoOfRuns = nNoOf2010hBppRuns;
    nFirstRun = NOOFRUNSPERSTEP;
    runlist = bpp2010hRunNumbers;
    break;
  default:
    Fatal("loadMC2010h11a10a_bisAODRunNumbers", "ERROR - Wrong step code %d - ABORTING", step);
  }
  for (Int_t run = nFirstRun; run < nNoOfRuns; run++)
    listOfActiveRuns.Add(new TObjString(runlist[run]));
}



void loadMC2010h11a10a_bisBMMRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hBMMMC";

  szDataDir = "/alice/sim/LHC11a10a_bis";
  szDataPattern = "*/AliESDs.root";

  bMC = kTRUE;

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nNoOf2010hBmmRuns; run++)
    listOfActiveRuns.Add(new TObjString(bmm2010hRunNumbers[run]));
}

void loadMC2010h11a10a_bisBPPRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hBPPMC";

  szDataDir = "/alice/sim/LHC11a10a_bis";
  szDataPattern = "*/AliESDs.root";

  bMC = kTRUE;

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nNoOf2010hBppRuns; run++)
    listOfActiveRuns.Add(new TObjString(bpp2010hRunNumbers[run]));
}

void load2013bTestRunNumber() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* teh GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2013b";

  /* 2013b */
  szDataDir = "/alice/data/2013/LHC13b";

  /* 2013B */
  szDataPattern = "*ESDs/pass3/*/AliESDs.root";

  listOfActiveRuns.Add(new TObjString("195344"));
}

void load2013bAODTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* teh GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2013b";

  /* 2013b */
  szDataDir = "/alice/data/2013/LHC13b";

  /* 2013B */
  szDataPattern = "*ESDs/pass3/*/AliAOD.root";

  listOfActiveRuns.Add(new TObjString("195344"));
}


void load2013bAODRunNumbers() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* teh GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2013bc";

  /* 2013b */
  szDataDir = "/alice/data/2013/LHC13b";

  /* 2013B */
  szDataPattern = "*pass4/AOD/*/AliAOD.root";

  const Int_t nRuns = 8;
  const char *runNumbers[nRuns] = {
    "195344",
    "195351",
    "195389",
    "195391",
    "195479",
    "195480",
    "195482",
    "195483"
  };

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nRuns; run++)
    listOfActiveRuns.Add(new TObjString(runNumbers[run]));
}

void load2013cAODRunNumbers() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* teh GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2013bc";

  /* 2013b */
  szDataDir = "/alice/data/2013/LHC13c";

  /* 2013B */
  szDataPattern = "*pass4/AOD/*/AliAOD.root";

  const Int_t nRuns = 14;
  const char *runNumbers[nRuns] = {
    "195529",
    "195531",
    "195566",
    "195567",
    "195568",
    "195592",
    "195593",
    "195596",
    "195633",
    "195635",
    "195644",
    "195673",
    "195675",
    "195677"
  };

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nRuns; run++)
    listOfActiveRuns.Add(new TObjString(runNumbers[run]));
}

void loadMC2013b2_efixRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2013b2";

  szDataDir = "/alice/sim/2013/LHC13b2_efix_p1";
  szDataPattern = "*/AliESDs.root";

  bMC = kTRUE;

  const Int_t nRuns = 23;
  const char *runNumbers[nRuns] = {
    "195344",
    "195351",
    "195389",
    "195391",
    "195478",
    "195479",
    "195480",
    "195482",
    "195483",
    "195529",
    "195531",
    "195566",
    "195567",
    "195568",
    "195592",
    "195593",
    "195596",
    "195633",
    "195635",
    "195644",
    "195673",
    "195675",
    "195677"
  };

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nRuns; run++)
    listOfActiveRuns.Add(new TObjString(runNumbers[run]));

}
void loadMCAMPT2010h12A11RunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2010hAMPT";

  szDataDir = "/alice/sim/2012/LHC12a11";
  szDataPattern = "*/AliESDs.root";

  bMC = kTRUE;

  /* the list of runs to analyze */
  listOfActiveRuns.Add(new TObjString("137686"));
  listOfActiveRuns.Add(new TObjString("138534"));
  listOfActiveRuns.Add(new TObjString("138653"));
  listOfActiveRuns.Add(new TObjString("139038"));
  listOfActiveRuns.Add(new TObjString("139437"));
}

void loadMCAMPT2010hCentrality(const char *cent) {

  loadMCAMPT2010h12A11RunNumbers();
  bUseMultiplicityTask       = kFALSE;

  if (TString(cent).Contains("0-5")) {
    szDataDir = szDataDir + "a";
  }
  else if (TString(cent).Contains("30-40")) {
    szDataDir = szDataDir + "e";
  }
  else if (TString(cent).Contains("60-70")) {
    szDataDir = szDataDir + "h";
  }
}

void load2015oLIRRunNumbers() {
  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2015oLIR";

  /* 2015o */
  szDataDir = "/alice/data/2015/LHC15o";

  /* 2015o LIR */
  szDataPattern = "*/pass3_lowIR_pidfix/*/AliESDs.root";

  /* reco pass */
  szpass = "3";

  /* the list of runs to analyze 2015o LIR*/
  /* (B: + +) */
  listOfActiveRuns.Add(new TObjString("244918"));
  listOfActiveRuns.Add(new TObjString("244975"));
  listOfActiveRuns.Add(new TObjString("244982"));
  listOfActiveRuns.Add(new TObjString("244983"));
  listOfActiveRuns.Add(new TObjString("245064"));
  listOfActiveRuns.Add(new TObjString("245066"));
  listOfActiveRuns.Add(new TObjString("245068"));
  /* (B: - -) */
  listOfActiveRuns.Add(new TObjString("246390"));
  listOfActiveRuns.Add(new TObjString("246391"));
  listOfActiveRuns.Add(new TObjString("246392"));
}

void load2015oLIRpass2RunNumbers() {

  /* load the 2015oLIR run numbers */
  load2015oLIRRunNumbers();

  /* reco pass */
  szpass = "2";

  /* 2015o LIR pass2 filtering */
  szDataPattern = "*/pass2_lowIR/*/AliESDs.root";
}

void load2015oLIRpass3RunNumbers() {

  /* load the 2015oLIR run numbers */
  load2015oLIRRunNumbers();

  /* reco pass */
  szpass = "3";

  /* 2015o LIR pass3 filtering */
  szDataPattern = "*/pass3_lowIR_pidfix/*/AliESDs.root";
}

void load2015oLIRpass4RunNumbers() {

  /* load the 2015oLIR run numbers */
  load2015oLIRRunNumbers();

  /* reco pass */
  szpass = "4";

  /* 2015o LIR pass4 filtering */
  szDataPattern = "*/pass4_lowIR_pidfix_cookdedx/*/AliESDs.root";
}


void load2015oMCTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2015oMC";

  bMC = kTRUE;

  /* reco pass */
  szpass = "1";

  /* 2015o */
  szDataDir = "/alice/sim/2016/LHC16g1";
  /* 2015o HIR */
  szDataPattern = "*/AliAOD.root";

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("246994"));
}

void load2015oHIRTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_2015oHIR";

  /* reco pass */
  szpass = "1";

  /* 2015o */
  szDataDir = "/alice/data/2015/LHC15o";
  /* 2015o HIR */
  szDataPattern = "*/pass1/*/AliAOD.root";

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("246087"));
}

void loadAMPT2015oHIRTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_AMPT2015oHIR";

  bMC = kTRUE;

  /* reco pass */
  szpass = "1";

  /* 2017i2 */
  szDataDir = "/alice/sim/2017/LHC17i2";
  /* 2017i2 */
  szDataPattern = "*/AOD/*/AliAOD.root";

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("246087"));
}

void load2016kTestRunNumber() {
  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "Test_pp_Studies_2016k";

  /* 2016ko */
  szDataDir = "/alice/data/2016/LHC16k";
  /* 2015o HIR */
  szDataPattern = "*/pass1/*/AliAOD.root";

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("258454"));
}

void load2017nRunNumbers() {
  /* the Xe-Xe dataset */
  bUseESD                   = kFALSE;
  bUseAOD                   = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudiesXeXe_2017n";

  /* reco pass */
  szpass = "1";

  /* 2017n */
  szDataDir = "/alice/data/2017/LHC17n";
  /* 2015o HIR */
  szDataPattern = "*/pass1/AOD/*/AliAOD.root";

  /* heavy data files */
  nNoOfInputFiles = 20;
  nNoOfTestFiles = 1;

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("280234"));
  listOfActiveRuns.Add(new TObjString("280235"));
}

void load2017nMCRunNumbers() {
  /* the Xe-Xe dataset */
  bUseESD                   = kFALSE;
  bUseAOD                   = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudiesXeXe_2017nMC";

  bMC = kTRUE;

  /* reco pass */
  szpass = "1";

  /* 2017n */
  szDataDir = "/alice/sim/2017/LHC17j7";
  /* 2015o HIR */
  szDataPattern = "*/AOD/*/AliAOD.root";

  /* heavy data files */
  nNoOfInputFiles = 20;
  nNoOfTestFiles = 1;

  /* the list of runs to analyze 2015o HIR*/
  listOfActiveRuns.Add(new TObjString("280234"));
  listOfActiveRuns.Add(new TObjString("280235"));
}

/* PbPb, AMPT, fast generation, 2.76TeV (min. bias), String melting ON, rescattering OFF, ID #48823 */
static const int nNoOf2013f3aRuns = 15;
static const char *s2013f3aRuns[nNoOf2013f3aRuns] = {
    "1", "2", "3", "4", "5", "6", "7", "8", "9",
    "11", "12", "13", "14", "15", "16"
};

void loadAMPT2760RunNumbers() {
  /* the fast AMPT dataset */
  bUseESD                   = kTRUE;
  bUseAOD                   = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudiesAMPTTruth";

  bMC                       = kTRUE;
  bMConlyTruth              = kTRUE;

  /* reco pass */
  szpass = "1";

  /* 2013f3 */
  szDataDir = "/alice/sim/2013/LHC13f3a";
  /*  fast AMPT */
  szDataPattern = "*/galice.root";

  /* heavy data files */
  nNoOfInputFiles = 20;
  nNoOfTestFiles = 1;

  /* the list of runs to analyze */
  for (Int_t run = 0; run < nNoOf2013f3aRuns; run++)
    listOfActiveRuns.Add(new TObjString(s2013f3aRuns[run]));

}

void load2018qRunNumbers() {
  const Int_t nruns = 14;
  const char *runs[nruns] = {
      "296068",
      "296066",
      "296065",
      "296063",
      "296062",
      "296061",
      "296060",
      "295589",
      "295588",
      "295587",
      "295586",
      "295585",
      "295584",
      "295581"
  };

  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_LHC2018q";

  bMC = kFALSE;

  /* reco pass */
  szpass = "1";

  /* 2018q */
  szDataDir = "/alice/data/2018/LHC18q";
  /* 2018q uncalibrated pass 1 */
  szDataPattern = "*/pass1_uncalibrated/*/AliESDs.root";

  /* the list of runs to analyze 2015o HIR*/
  for (Int_t i = 0; i < nruns; i++)
    listOfActiveRuns.Add(new TObjString(runs[i]));

}

void load2018rMCRunNumbers() {
  const Int_t nruns = 241;
  char buffer[256];
  int runsno[nruns] = {
      295585, 295586, 295587, 295588, 295589, 295610, 295611, 295612, 295615, 295665,
      295666, 295667, 295668, 295671, 295673, 295675, 295676, 295677, 295712, 295714,
      295717, 295718, 295719, 295721, 295723, 295725, 295753, 295754, 295755, 295756,
      295758, 295759, 295762, 295763, 295786, 295788, 295791, 295816, 295818, 295819,
      295822, 295825, 295826, 295829, 295831, 295853, 295854, 295855, 295856, 295859,
      295860, 295861, 295881, 295908, 295909, 295910, 295913, 295936, 295937, 295941,
      295942, 295943, 295945, 295947, 296016, 296060, 296061, 296062, 296063, 296065,
      296066, 296068, 296074, 296123, 296132, 296133, 296134, 296135, 296142, 296143,
      296191, 296192, 296194, 296195, 296196, 296197, 296198, 296240, 296241, 296242,
      296243, 296244, 296246, 296247, 296269, 296270, 296273, 296279, 296280, 296303,
      296304, 296307, 296309, 296312, 296375, 296376, 296377, 296378, 296379, 296380,
      296381, 296383, 296414, 296415, 296419, 296420, 296423, 296424, 296433, 296472,
      296509, 296510, 296511, 296512, 296514, 296516, 296547, 296548, 296549, 296550,
      296551, 296552, 296553, 296594, 296615, 296616, 296618, 296619, 296621, 296622,
      296623, 296690, 296691, 296693, 296694, 296749, 296750, 296752, 296781, 296784,
      296785, 296786, 296787, 296790, 296793, 296794, 296799, 296835, 296836, 296838,
      296839, 296848, 296849, 296850, 296851, 296852, 296890, 296894, 296899, 296900,
      296903, 296930, 296931, 296932, 296934, 296935, 296938, 296941, 296966, 297029,
      297031, 297035, 297085, 297117, 297118, 297119, 297123, 297124, 297128, 297129,
      297132, 297133, 297193, 297194, 297195, 297196, 297218, 297219, 297221, 297222,
      297278, 297310, 297311, 297312, 297315, 297317, 297332, 297333, 297335, 297336,
      297363, 297366, 297367, 297372, 297379, 297380, 297405, 297406, 297413, 297414,
      297415, 297441, 297442, 297446, 297450, 297451, 297452, 297479, 297481, 297483,
      297512, 297537, 297540, 297541, 297542, 297544, 297558, 297588, 297590, 297595,
      297624
  };

  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_LHC2018r";

  bMC = kTRUE;

  /* reco pass */
  szpass = "1";

  /* 2018q */
  szDataDir = "/alice/sim/2018/LHC18l8a4";
  /* 2018q uncalibrated pass 1 */
  szDataPattern = "*/AOD/*/AliAOD.root";

  /* the list of runs to analyze 2018l8a4*/
  for (Int_t i = 0; i < 2 /*nruns*/; i++) {
    sprintf(buffer,"%d",runsno[i]);
    listOfActiveRuns.Add(new TObjString(buffer));
  }
}

void load2018rRunNumbers() {
  const Int_t nruns = 92;
  char buffer[256];
  int runsno[nruns] = {
      297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512,
      297483, 297481, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415,
      297414, 297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363,
      297336, 297335, 297333, 297332, 297317, 297315, 297312, 297311, 297310, 297278,
      297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128,
      297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941,
      296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899, 296894,
      296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794,
      296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693,
      296691, 296690
  };

  bUseESD                    = kFALSE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_LHC2018r";

  bMC = kFALSE;

  /* reco pass */
  szpass = "1";

  /* 2018q */
  szDataDir = "/alice/data/2018/LHC18r";
  /* 2018q uncalibrated pass 1 */
  szDataPattern = "*/pass1/AOD/*/AliAOD.root";

  /* the list of runs to analyze 2018r */
  for (Int_t i = 0; i < 2 /*nruns*/; i++) {
    sprintf(buffer,"%d",runsno[i]);
    listOfActiveRuns.Add(new TObjString(buffer));
  }
}

void load2018rESDRunNumbers() {
  const Int_t nruns = 92;
  char buffer[256];
  int runsno[nruns] = {
      297595, 297590, 297588, 297558, 297544, 297542, 297541, 297540, 297537, 297512,
      297483, 297481, 297479, 297452, 297451, 297450, 297446, 297442, 297441, 297415,
      297414, 297413, 297406, 297405, 297380, 297379, 297372, 297367, 297366, 297363,
      297336, 297335, 297333, 297332, 297317, 297315, 297312, 297311, 297310, 297278,
      297222, 297221, 297218, 297196, 297195, 297193, 297133, 297132, 297129, 297128,
      297124, 297123, 297119, 297118, 297117, 297085, 297035, 297031, 296966, 296941,
      296938, 296935, 296934, 296932, 296931, 296930, 296903, 296900, 296899, 296894,
      296852, 296851, 296850, 296848, 296839, 296838, 296836, 296835, 296799, 296794,
      296793, 296790, 296787, 296786, 296785, 296784, 296781, 296752, 296694, 296693,
      296691, 296690
  };

  bUseESD                    = kTRUE;
  bUseAOD                    = !bUseESD;

  /* the GRID working directory */
  szGridWorkingDir = "CorrelationStudies_LHC2018r";

  bMC = kFALSE;

  /* reco pass */
  szpass = "1";

  /* 2018q */
  szDataDir = "/alice/data/2018/LHC18r";
  /* 2018q uncalibrated pass 1 */
  szDataPattern = "*/pass1/*/AliESDs.root";

  /* the list of runs to analyze 2018r */
  for (Int_t i = 0; i < 2 /*nruns*/; i++) {
    sprintf(buffer,"%d",runsno[i]);
    listOfActiveRuns.Add(new TObjString(buffer));
  }
}

