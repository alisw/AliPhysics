#ifdef __ECLIPSE_IDE
//  few includes and external declarations just for the IDE
#include <Riostream.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TObjString.h>

#endif // #ifdef __ECLIPSE_IDE

#include "runCorrelationsStudiesConfigMacro.H"

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
void load2015oMCTestRunNumber();
void load2015oLIRpass2RunNumbers();
void load2015oLIRpass3RunNumbers();
void load2015oLIRpass4RunNumbers();
void load2015oHIRRunNumbers();
void load2015oHIRTestRunNumber();
void load2016kTestRunNumber();
void load2013bTestRunNumber();
void loadAMPT2015oHIRTestRunNumber();
void load2017nRunNumbers();
void load2017nMCRunNumbers();
void loadAMPT2760RunNumbers();

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
  szAliPhysicsVersion = "vAN-20171222-1";

  /* the number of files we want to test */
  nNoOfInputFiles = 30;
  nNoOfTestFiles = 2;

  /* load the run numbers */
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
  load2013bTestRunNumber();
  // loadMCAMPT2010hCentrality("0-5");
  // loadAMPT2015oHIRTestRunNumber();
  // load2017nRunNumbers();
  // loadAMPT2760RunNumbers();

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
