#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>
#include <AliPIDResponse.h>
#include <Riostream.h>
#include <AliLog.h>
#include <AliAnalysisGrid.h>
#include <AliAnalysisAlien.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisTaskPIDqa.h>
#include <AliAnalysisTaskTOFSpectra.h>
#include <AliMultSelectionTask.h>
#include <AliESDInputHandler.h>
#include <AliAnalysisTaskPIDResponse.h>
#include <AliMCEventHandler.h>
#include <AliAnalysisDataContainer.h>
#include <AddTaskTOFSpectra.C>
#include <TSystem.h>
#define MODE "COMPILED"
#else
class AliAnalysisGrid;
class AliAnalysisAlien;
#define MODE "INTERPRETED"
#endif

fstream logFile;
TString startitme, optStr;

Bool_t optTree = kTRUE;//Option to use the Tree
const Bool_t optHeavyIon = kTRUE;//Option to use Pb--Pb data instead of pp
const Bool_t optMismatchrun = kFALSE;//option to run with tof channel
const Bool_t optCutVar = optTree ? kTRUE : kFALSE;//Cut variation for the tree!!
const Bool_t optSimpleCutVar = kTRUE;//Cut variation for the list!!
const Bool_t optUseHijing = kTRUE;//Option to use HIJING instead of DPMJET
const Bool_t optJDLMerge = optTree ? kFALSE : kTRUE;//Option to merge using JDL
const Bool_t optRemoveCut = kFALSE;//Option to remove one cut at a time, the cut for tracks

const Int_t optTestMC = -1;//Option to test the analysis on the new MC productions
const Bool_t optHiInteractionRate = kFALSE;//Option to use the HI interaction rage

//Data Periods
TString period = "";
TString pass = "";

AliAnalysisGrid* CreateAlienHandler(const TString runmode, const Bool_t isMC, const TString fname);

void RunTaskTOFSpectra(const TString runmode = "test", const Bool_t isMC = kFALSE, const TString fname = "test", const Int_t hi = -1){
  cout<<"Running in "<<MODE<<" Mode"<<endl;
  AliLog::SetGlobalDebugLevel(AliLog::kInfo);
  
  if(!runmode.EqualTo("test") && !runmode.EqualTo("terminate") && !runmode.EqualTo("full")){
    cout<<"Uknown run option, Aborting!"<<endl;    
    return;
  }
  
  ///////////////////////////////////////////////////////Log
  TString logStr= Form("***Runmode(%s) MC(%i) AnalName(%s)***",runmode.Data(),isMC,fname.Data());
  optStr=Form("OPT: R_%s M_%i",runmode.Data(),isMC);
  gSystem->Exec("date >> logRunAnalysis.txt");
  startitme=gSystem->GetFromPipe("date");
  logFile.open("logRunAnalysis.txt",ios::in | ios::out | ios::ate);
  if (logFile.is_open()){
    logFile << logStr<<endl;
    logFile<<"DATE: "<<startitme<<endl;
    logFile<<optStr<<endl;    
  }
  else cout << "Unable to open file";
  //////////////////////////////////////////////////////
  
  //Include path
  //   gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/ANALYSIS/ -I${ALICE_ROOT}/STEER/ESD -I${ALICE_ROOT}/include ");
  //   gSystem->SetIncludePath("-I$ROOTSYS/include -I${ALICE_ROOT}/include -I${ALICE_ROOT}/ANALYSIS/macros/ -I${ALICE_PHYSICS}/include/ -I${ALICE_PHYSICS}/OADB/macros ");
  
  //__________________________________________________________________________
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/TOF ");
  
  cout<<"Include path: "<<gSystem->GetIncludePath()<<endl;
  
  // Load common libraries
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libProof");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libRAWDatarec.so");
  gSystem->Load("libRAWDatasim.so");
  gSystem->Load("libTOFrec.so");
  gSystem->Load("libASImage.so");
  gSystem->Load("libHistPainter.so");  
  gSystem->Load("libProofPlayer.so");
  gSystem->Load("libTreePlayer");
  gSystem->Load("libHist.so");
  //   gSystem->Load("libTENDER.so"); 
  //   gSystem->Load("libTENDERSupplies.so");
  
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  
  if(isMC) optTree = kFALSE;
  if(optMismatchrun) optTree = kFALSE;
  
  // Create and configure the alien handler plugin
  AliAnalysisGrid *alienHandler = CreateAlienHandler(runmode,isMC,fname);  
  if (!alienHandler) return;
  mgr->SetGridHandler(alienHandler);
  logFile.close();
  
  AliESDInputHandler* esdHanldler = new AliESDInputHandler();
  esdHanldler->SetNeedField(kTRUE);
  mgr->SetInputEventHandler(esdHanldler);
  
  //__________________________________________________________________________
  //Connect to OCDB
  //   cout<<"Loading TaskOCDBConnect"<<endl;
  //   gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
  //   AddTaskCDBconnect();
  
  //__________________________________________________________________________
  // Apply the event selection
  cout<<"Loading TaskPhysicsSelection"<<endl;
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
  AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(isMC,kTRUE);
  
  if(isMC){
    AliMCEventHandler *mcHandler = new AliMCEventHandler();
    mgr->SetMCtruthEventHandler(mcHandler);
    physSelTask->GetPhysicsSelection()->SetAnalyzeMC(); 
  }
  
  //__________________________________________________________________________
  
  //DEPRECATED
  //   cout<<"Loading TaskCentrality"<<endl;
  //   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  //   AliCentralitySelectionTask *centsSelTask = AddTaskCentrality(kTRUE,kFALSE);
  
  if(optHeavyIon){
    cout<<"Loading TaskMultSelection"<<endl;
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask * multSelTask = AddTaskMultSelection(kFALSE); // user mode:
    //   AliMultSelectionTask * multSelTask = AddTaskMultSelection(kTRUE); // calibration mode:, To run in calibration mode, the task must run AFTER the Physics Selection
    //In case you want to change the trigger class used for calibration, you need to call 
    //   multSelTask->SetSelectedTriggerClass(AliVEvent::kINT7); // kINT7 is default, this is OK for Run2; in LHC10h you need kMB
    //In case you want to use the default calibration for runs which have not yet been calibrated, you should call (for data and MC respectively)
    //   if(!isMC) multSelTask->SetUseDefaultCalib(kTRUE); // data
    //   else multSelTask->SetUseDefaultMCCalib(kTRUE); // MC 
  }
  
  
  //__________________________________________________________________________
  
  cout<<"Loading TaskPIDResponse"<<endl;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *setupTask = AddTaskPIDResponse(isMC);
  //   AliAnalysisTaskPIDResponse *setupTask = AddTaskPIDResponse(isMC,kTRUE,kTRUE, pass.Contains("pass1") ? 1 : 2);
  
  //__________________________________________________________________________
  
  cout<<"Loading TaskPIDqa"<<endl;
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
  AliAnalysisTaskPIDqa *pidQA = AddTaskPIDqa();
  
  //__________________________________________________________________________
  
  //gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/TenderSupplies/AddTaskTender.C");
  //AliAnalysisTask* tender=0x0;
  //if(mc==0)
  //  {
  //    tender = AddTaskTender(kTRUE,kTRUE,kFALSE,kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE);
  // tender->SetDebugLevel(10);
  // }
  //else
  //  {
  //   tender = AddTaskTender(kFALSE,kTRUE,kFALSE,kTRUE,kTRUE,kTRUE,kTRUE,kFALSE,kFALSE);
  // tender->SetDebugLevel(10);
  //}
  
  //__________________________________________________________________________
  cout<<"Loading TaskTOFSpectra"<<endl;
  gROOT->LoadMacro("./AliAnTOFtrack.cxx++g");
  gROOT->LoadMacro("./AliAnalysisTaskTOFSpectra.cxx++g");
  gROOT->LoadMacro("./AddTaskTOFSpectra.C");
  gROOT->LoadMacro("./runutils.h");
  const TString tname = "TOFSpectra";
  AliAnalysisTaskTOFSpectra *taskTOFSpectra = AddTaskTOFSpectra(optTree, isMC, optHeavyIon, optMismatchrun, optCutVar, -1, "", tname);//Standard task
  if(optSimpleCutVar){//Cut variation
    if(1){//TPC Rows
      const Int_t offset = kTPCSetL;
      const UInt_t cut = kTPCrows;
      Int_t index = 0;
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_TPCRows[CutIndex[cut]-1];//Remove the standard case
      for(Int_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the Std value 
        
        const TString prefix = Form("%s", CutVarsName[offset + index].Data());
        taskTOFSpectra_CutVar_TPCRows[index] = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, offset + index, prefix, Form("%s_%s%i", tname.Data(), Cuts[cut].Data(), index));
        index++;
      }
      
    }
    
    if(1){//Track Chi2 
      const Int_t offset = kTPCChi2SetL;
      const UInt_t cut = kTrkChi2;
      Int_t index = 0;
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_TrackChi2[CutIndex[cut]-1];//Remove the standard case
      for(Int_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the Std value 
        
        const TString prefix = Form("%s", CutVarsName[offset + index].Data());
        taskTOFSpectra_CutVar_TrackChi2[index] = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, offset + index, prefix, Form("%s_%s%i", tname.Data(), Cuts[cut].Data(), index));
        index++;
      }
      
    }
    
    if(1){//DCAz 
      const Int_t offset = kDCAzSetL;
      const UInt_t cut = kDCAz;
      Int_t index = 0;
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_DCAz[CutIndex[cut]-1];//Remove the standard case
      for(Int_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the Std value 
        
        const TString prefix = Form("%s", CutVarsName[offset + index].Data());
        taskTOFSpectra_CutVar_DCAz[index] = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, offset + index, prefix, Form("%s_%s%i", tname.Data(), Cuts[cut].Data(), index));
        index++;
      }
      
    }
    
    if(1){//DCAxy 
      const Int_t offset = kPrimSetL;
      const UInt_t cut = kDCAxy;
      Int_t index = 0;
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_DCAxy[CutIndex[cut]-1];//Remove the standard case
      for(Int_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the Std value 
        
        const TString prefix = Form("%s", CutVarsName[offset + index].Data());
        taskTOFSpectra_CutVar_DCAxy[index] = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, offset + index, prefix, Form("%s_%s%i", tname.Data(), Cuts[cut].Data(), index));
        index++;
      }
      
    }
    
    if(1){//GeoCut
      const Int_t offset = kGeoCutSet1;
      const UInt_t cut = kGeo;
      Int_t index = 0;
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_Geo[CutIndex[cut]-1];//Remove the standard case
      for(Int_t i = 0; i < CutIndex[cut]; i++){
        if(i == CutStdIndex[cut]) continue;//Skip the Std value 
        
        const TString prefix = Form("%s", CutVarsName[offset + index].Data());
        taskTOFSpectra_CutVar_Geo[index] = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, offset + index, prefix, Form("%s_%s%i", tname.Data(), Cuts[cut].Data(), index));
        index++;
      }
      
    }
    
  }
  
  if(optRemoveCut){//Cut removal
    if(1){//TPC Rows
      const UInt_t cut = kTPCrows;
      const TString prefix = Form("No%s", Cuts[cut].Data());
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoTPCRow = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_No%s", tname.Data(), Cuts[cut].Data()));
      //       taskTOFSpectra_CutVar_NoTPCRow->SetTPCRowsCut(0);
      
      AliESDtrackCuts* trkcus = (AliESDtrackCuts*) AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE, 1);//WARNING KEEP THE DCA SO THAT YOU CAN USE THE SECONDARIES
      trkcus->SetCutGeoNcrNcl(3., 130., 1.5, 0.85, 0.7);
      trkcus->SetMinNCrossedRowsTPC(0);
      trkcus->SetName("asdasdolad");
      taskTOFSpectra_CutVar_NoTPCRow->fESDtrackCuts = trkcus;
      taskTOFSpectra_CutVar_NoTPCRow->asd = 123;
    }
    
    if(0){//Track Chi2 
      const UInt_t cut = kTrkChi2;
      const TString prefix = Form("No%s", Cuts[cut].Data());
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoTrackChi2 = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_No%s", tname.Data(), Cuts[cut].Data()));
      taskTOFSpectra_CutVar_NoTrackChi2->SetTrkChi2Cut(0);
    }
    
    if(0){//Track Chi2 in ITS
      const TString prefix = Form("NoITSChi2");
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoTrackITSChi2 = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_%s", tname.Data(), prefix.Data()));
      taskTOFSpectra_CutVar_NoTrackITSChi2->SetTrkChi2CutITS(0);
    }
    
    if(0){//Ratio of crossed rows and findable clusters
      const TString prefix = Form("NoCrsClsRatio");
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoCrsClsRatio = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_%s", tname.Data(), prefix.Data()));
      taskTOFSpectra_CutVar_NoCrsClsRatio->SetRatioCrossedRowsFindableCls(0);
    }
    
    if(0){//DCAz 
      const UInt_t cut = kDCAz;
      const TString prefix = Form("No%s", Cuts[cut].Data());
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoDCAz = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_No%s", tname.Data(), Cuts[cut].Data()));
      taskTOFSpectra_CutVar_NoDCAz->SetDCAzCut(10);
    }
    
    if(0){//DCAxy 
      const UInt_t cut = kDCAxy;
      const TString prefix = Form("No%s", Cuts[cut].Data());
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoDCAxy = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_No%s", tname.Data(), Cuts[cut].Data()));
      taskTOFSpectra_CutVar_NoDCAxy->SetDCAxyCut(100);
    }
    
    if(0){//GeoCut 
      const UInt_t cut = kGeo;
      const TString prefix = Form("No%s", Cuts[cut].Data());
      AliAnalysisTaskTOFSpectra *taskTOFSpectra_CutVar_NoGeo = AddTaskTOFSpectra(kFALSE, isMC, optHeavyIon, optMismatchrun, kFALSE, -1, prefix, Form("%s_No%s", tname.Data(), Cuts[cut].Data()));
      taskTOFSpectra_CutVar_NoGeo->SetGeoCut(0, 0, 0, 0, 0);
    }
    
  }
  
  //__________________________________________________________________________
  
  
  mgr->InitAnalysis();
  mgr->PrintStatus();
  // Start analysis in grid.
  mgr->StartAnalysis("grid");
};


AliAnalysisGrid* CreateAlienHandler(const TString runmode, const Bool_t isMC, const TString fname) {
  // Check if user has a valid token, otherwise make one. This has limitations
  // One can always follow the standard procedure of calling alien-token-init then
  //   source /tmp/gclient_env_$UID in the current shell.
  //if (!AliAnalysisGrid::CreateToken()) return NULL;
  AliAnalysisAlien *plugin = new AliAnalysisAlien();
  // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
  plugin->SetRunMode(runmode.Data());
  plugin->SetNtestFiles(1);
  // Set versions of used packages
  plugin->SetAPIVersion("V1.1x");
  const TString swversion = "vAN-20161110-1";
  plugin->SetAliPhysicsVersion(swversion);
  logFile<<swversion<<endl;
  
  TString griddata = "/alice/";
  if(isMC) griddata += "sim/";
  else griddata += "data/";
  
  const TString dataprefix = isMC ? "" : "000";
  const TString mcsuffix = isMC ? "_MC" : "_DATA";
  const TString taskname = (optHeavyIon ? "taskSpectraPbPb5" : "taskSpectrapp5") + mcsuffix;
  const TString outdir = optHeavyIon ? "SpectraPbPb/SpectraPbPb2015/" : "Spectrapp/Spectrapp2015/";
  const Int_t nsplit = isMC ? 40 : 30;  
  
  TString partialmc = "";//Partial list for partial productions
  
  if(period.EqualTo("")){//If not set set it!
    if(isMC){//MC
      if(optHeavyIon){//Hi
        
        //LHC16g1, LHC16g1a, LHC16g1b, LHC16g1c, LHC16g2, LHC16g3, LHC16g4
        if(optTestMC > 0){//Test of general purpose
          switch (optTestMC) {
            case 0:
            period = "LHC16g1";//Hijing min.bias
            partialmc = " 245064,";
            break;
            case 1:
            period = "LHC16g1a";//Hijing 0-10%
            partialmc = " 245064,";
            break;
            case 2:
            period = "LHC16g1b";//Hijing 10-50%
            partialmc = " 245064,";
            break;
            case 3:
            period = "LHC16g1c";//Hijing 50-90%
            partialmc = " 245064,";
            break;
            case 4:
            period = "LHC16g2";//EPOS-LHC
            partialmc = " ";
            break;
            case 5:
            period = "LHC16g3";//DPMJET 1%
            partialmc = " ";
            break;
            case 6:
            period = "LHC16g4";//Hijing high-pt 
            partialmc = " 244917, 244918, 244975, 244980,";
            break;
            case 7:
            period = "LHC16j7a";//Hijing general purpose 
            partialmc = " 244982, 245064, 246392,";
            break;
            case 8:
            period = "LHC16j7b";//Hijing general purpose
            partialmc = " 244982, 245064, 246392,";
            break;
            default:
            Fatal("CreateAlienHandler", "mc option %i not recognized", optTestMC);
            break;
          }
        }
        else{//Standard
          if(optUseHijing) period = "LHC15k1a1";//HIJING
          else period = "LHC15k1b1";//DPMJET
        }
        
      }
      else{//pp
        period = "LHC15l1a2";//PYTHIA8
      }
    }
    else{//Data
      if(optHeavyIon){//Hi
        period = "LHC15o";
        if(optHiInteractionRate) pass = "pass1";
        // else pass = "pass2_lowIR";
        else pass = "pass4_lowIR_pidfix_cookdedx";
      }
      else{//pp
        period = "LHC15n";
        pass = "pass2";
      }
    }
  }
  
  if(!period.EqualTo("")){ //Setting the data periods
    
    if(period.Contains("LHC15")) griddata += "2015/";
    else if(period.Contains("LHC16")) griddata += "2016/";
    else Fatal("CreateAlienHandler", "Cannot find year in %s", period.Data());
    
    griddata += period;
    
    plugin->SetGridDataDir(griddata);//DPMJET
    plugin->SetDataPattern(Form("*%s%s/*AliESDs.root", pass.EqualTo("") ? "" : "/", pass.Data()));
    plugin->SetRunPrefix(dataprefix);     
    plugin->SetGridWorkingDir(Form("%s%s%s", outdir.Data(), fname.Data(), mcsuffix.Data()));
    plugin->SetExecutable(Form("%s%s.sh", fname.Data(), mcsuffix.Data()));
    plugin->SetJDLName(Form("%s%s.jdl", taskname.Data(), mcsuffix.Data()));
    plugin->SetAnalysisMacro(Form("%s%s.C", taskname.Data(), mcsuffix.Data()));
    plugin->SetSplitMaxInputFileNumber(nsplit);
    
  }  
  else Fatal("CreateAlienHandler", "No data period provided");
  
  //Set number of runs
  Int_t* runs = 0x0;
  Int_t nruns = 0;
  if(isMC){//MC
    if(optHeavyIon){
      if(optHiInteractionRate) nruns = 161;
      else nruns = 12;
    }
    else nruns = 7;
  }
  else{//Data
    if(optHeavyIon){
      if(optHiInteractionRate) nruns = 161;
      else nruns = pass.EqualTo("pass2_lowIR") ? 13 : 12;
    }
    else nruns = 25;
  }
  
  //Default
  Int_t ntorun = nruns;
  Int_t ntoskip = 0;
  Int_t nrunssplit = 1;
  
  runs = new Int_t[nruns];
  
  //Set runlist
  Int_t counter = 0;
  if(isMC){//MC
    if(optHeavyIon){//Hi
      if(optHiInteractionRate){
        runs[counter++] = 245145;
        runs[counter++] = 245146;
        runs[counter++] = 245148;
        runs[counter++] = 245151;
        runs[counter++] = 245152;
        runs[counter++] = 245231;
        runs[counter++] = 245232;
        runs[counter++] = 245233;
        runs[counter++] = 245259;
        runs[counter++] = 245343;
        runs[counter++] = 245345;
        runs[counter++] = 245346;
        runs[counter++] = 245347;
        runs[counter++] = 245349;
        runs[counter++] = 245353;
        runs[counter++] = 245396;
        runs[counter++] = 245397;
        runs[counter++] = 245401;
        runs[counter++] = 245407;
        runs[counter++] = 245409;
        runs[counter++] = 245410;
        runs[counter++] = 245411;
        runs[counter++] = 245439;
        runs[counter++] = 245441;
        runs[counter++] = 245446;
        runs[counter++] = 245450;
        runs[counter++] = 245452;
        runs[counter++] = 245453;
        runs[counter++] = 245454;
        runs[counter++] = 245496;
        runs[counter++] = 245497;
        runs[counter++] = 245501;
        runs[counter++] = 245504;
        runs[counter++] = 245505;
        runs[counter++] = 245507;
        runs[counter++] = 245535;
        runs[counter++] = 245540;
        runs[counter++] = 245542;
        runs[counter++] = 245543;
        runs[counter++] = 245544;
        runs[counter++] = 245545;
        runs[counter++] = 245554;
        runs[counter++] = 245683;
        runs[counter++] = 245692;
        runs[counter++] = 245700;
        runs[counter++] = 245702;
        runs[counter++] = 245705;
        runs[counter++] = 245729;
        runs[counter++] = 245731;
        runs[counter++] = 245738;
        runs[counter++] = 245752;
        runs[counter++] = 245759;
        runs[counter++] = 245766;
        runs[counter++] = 245775;
        runs[counter++] = 245785;
        runs[counter++] = 245793;
        runs[counter++] = 245829;
        runs[counter++] = 245831;
        runs[counter++] = 245833;
        runs[counter++] = 245923;
        runs[counter++] = 245949;
        runs[counter++] = 245952;
        runs[counter++] = 245954;
        runs[counter++] = 245963;
        runs[counter++] = 245985;
        runs[counter++] = 245996;
        runs[counter++] = 246001;
        runs[counter++] = 246003;
        runs[counter++] = 246012;
        runs[counter++] = 246036;
        runs[counter++] = 246037;
        runs[counter++] = 246042;
        runs[counter++] = 246048;
        runs[counter++] = 246049;
        runs[counter++] = 246052;
        runs[counter++] = 246053;
        runs[counter++] = 246087;
        runs[counter++] = 246089;
        runs[counter++] = 246113;
        runs[counter++] = 246115;
        runs[counter++] = 246148;
        runs[counter++] = 246151;
        runs[counter++] = 246152;
        runs[counter++] = 246153;
        runs[counter++] = 246178;
        runs[counter++] = 246180;
        runs[counter++] = 246181;
        runs[counter++] = 246182;
        runs[counter++] = 246185;
        runs[counter++] = 246187;
        runs[counter++] = 246217;
        runs[counter++] = 246222;
        runs[counter++] = 246225;
        runs[counter++] = 246271;
        runs[counter++] = 246272;
        runs[counter++] = 246275;
        runs[counter++] = 246276;
        runs[counter++] = 246424;
        runs[counter++] = 246428;
        runs[counter++] = 246431;
        runs[counter++] = 246433;
        runs[counter++] = 246434;
        runs[counter++] = 246487;
        runs[counter++] = 246488;
        runs[counter++] = 246493;
        runs[counter++] = 246495;
        runs[counter++] = 246540;
        runs[counter++] = 246543;
        runs[counter++] = 246553;
        runs[counter++] = 246567;
        runs[counter++] = 246568;
        runs[counter++] = 246575;
        runs[counter++] = 246583;
        runs[counter++] = 246639;
        runs[counter++] = 246648;
        runs[counter++] = 246671;
        runs[counter++] = 246675;
        runs[counter++] = 246676;
        runs[counter++] = 246750;
        runs[counter++] = 246751;
        runs[counter++] = 246755;
        runs[counter++] = 246757;
        runs[counter++] = 246758;
        runs[counter++] = 246759;
        runs[counter++] = 246760;
        runs[counter++] = 246763;
        runs[counter++] = 246765;
        runs[counter++] = 246766;
        runs[counter++] = 246804;
        runs[counter++] = 246805;
        runs[counter++] = 246806;
        runs[counter++] = 246807;
        runs[counter++] = 246808;
        runs[counter++] = 246809;
        runs[counter++] = 246810;
        runs[counter++] = 246844;
        runs[counter++] = 246845;
        runs[counter++] = 246846;
        runs[counter++] = 246847;
        runs[counter++] = 246851;
        runs[counter++] = 246855;
        runs[counter++] = 246858;
        runs[counter++] = 246859;
        runs[counter++] = 246864;
        runs[counter++] = 246865;
        runs[counter++] = 246867;
        runs[counter++] = 246870;
        runs[counter++] = 246871;
        runs[counter++] = 246928;
        runs[counter++] = 246930;
        runs[counter++] = 246937;
        runs[counter++] = 246942;
        runs[counter++] = 246945;
        runs[counter++] = 246948;
        runs[counter++] = 246949;
        runs[counter++] = 246980;
        runs[counter++] = 246982;
        runs[counter++] = 246984;
        runs[counter++] = 246989;
        runs[counter++] = 246991;
        runs[counter++] = 246994;
      }
      else{
        runs[counter++] = 244918;
        runs[counter++] = 244975;
        runs[counter++] = 244982;
        runs[counter++] = 244983;
        runs[counter++] = 245064;
        runs[counter++] = 245066;
        runs[counter++] = 245068;
        runs[counter++] = 246392;
        runs[counter++] = 246391;
        runs[counter++] = 246390;
        runs[counter++] = 245145;
        runs[counter++] = 244917;
      }
    }
    else{//pp
      runs[counter++] = 244377;
      runs[counter++] = 244364;
      runs[counter++] = 244359;
      runs[counter++] = 244355;
      runs[counter++] = 244351;
      runs[counter++] = 244343;
      runs[counter++] = 244340;
    }
  }
  else{//Data
    if(optHeavyIon){//Hi
      if(optHiInteractionRate){
        runs[counter++] = 245145;
        runs[counter++] = 245146;
        runs[counter++] = 245148;
        runs[counter++] = 245151;
        runs[counter++] = 245152;
        runs[counter++] = 245231;
        runs[counter++] = 245232;
        runs[counter++] = 245233;
        runs[counter++] = 245259;
        runs[counter++] = 245343;
        runs[counter++] = 245345;
        runs[counter++] = 245346;
        runs[counter++] = 245347;
        runs[counter++] = 245349;
        runs[counter++] = 245353;
        runs[counter++] = 245396;
        runs[counter++] = 245397;
        runs[counter++] = 245401;
        runs[counter++] = 245407;
        runs[counter++] = 245409;
        runs[counter++] = 245410;
        runs[counter++] = 245411;
        runs[counter++] = 245439;
        runs[counter++] = 245441;
        runs[counter++] = 245446;
        runs[counter++] = 245450;
        runs[counter++] = 245452;
        runs[counter++] = 245453;
        runs[counter++] = 245454;
        runs[counter++] = 245496;
        runs[counter++] = 245497;
        runs[counter++] = 245501;
        runs[counter++] = 245504;
        runs[counter++] = 245505;
        runs[counter++] = 245507;
        runs[counter++] = 245535;
        runs[counter++] = 245540;
        runs[counter++] = 245542;
        runs[counter++] = 245543;
        runs[counter++] = 245544;
        runs[counter++] = 245545;
        runs[counter++] = 245554;
        runs[counter++] = 245683;
        runs[counter++] = 245692;
        runs[counter++] = 245700;
        runs[counter++] = 245702;
        runs[counter++] = 245705;
        runs[counter++] = 245729;
        runs[counter++] = 245731;
        runs[counter++] = 245738;
        runs[counter++] = 245752;
        runs[counter++] = 245759;
        runs[counter++] = 245766;
        runs[counter++] = 245775;
        runs[counter++] = 245785;
        runs[counter++] = 245793;
        runs[counter++] = 245829;
        runs[counter++] = 245831;
        runs[counter++] = 245833;
        runs[counter++] = 245923;
        runs[counter++] = 245949;
        runs[counter++] = 245952;
        runs[counter++] = 245954;
        runs[counter++] = 245963;
        runs[counter++] = 245985;
        runs[counter++] = 245996;
        runs[counter++] = 246001;
        runs[counter++] = 246003;
        runs[counter++] = 246012;
        runs[counter++] = 246036;
        runs[counter++] = 246037;
        runs[counter++] = 246042;
        runs[counter++] = 246048;
        runs[counter++] = 246049;
        runs[counter++] = 246052;
        runs[counter++] = 246053;
        runs[counter++] = 246087;
        runs[counter++] = 246089;
        runs[counter++] = 246113;
        runs[counter++] = 246115;
        runs[counter++] = 246148;
        runs[counter++] = 246151;
        runs[counter++] = 246152;
        runs[counter++] = 246153;
        runs[counter++] = 246178;
        runs[counter++] = 246180;
        runs[counter++] = 246181;
        runs[counter++] = 246182;
        runs[counter++] = 246185;
        runs[counter++] = 246187;
        runs[counter++] = 246217;
        runs[counter++] = 246222;
        runs[counter++] = 246225;
        runs[counter++] = 246271;
        runs[counter++] = 246272;
        runs[counter++] = 246275;
        runs[counter++] = 246276;
        runs[counter++] = 246424;
        runs[counter++] = 246428;
        runs[counter++] = 246431;
        runs[counter++] = 246433;
        runs[counter++] = 246434;
        runs[counter++] = 246487;
        runs[counter++] = 246488;
        runs[counter++] = 246493;
        runs[counter++] = 246495;
        runs[counter++] = 246540;
        runs[counter++] = 246543;
        runs[counter++] = 246553;
        runs[counter++] = 246567;
        runs[counter++] = 246568;
        runs[counter++] = 246575;
        runs[counter++] = 246583;
        runs[counter++] = 246639;
        runs[counter++] = 246648;
        runs[counter++] = 246671;
        runs[counter++] = 246675;
        runs[counter++] = 246676;
        runs[counter++] = 246750;
        runs[counter++] = 246751;
        runs[counter++] = 246755;
        runs[counter++] = 246757;
        runs[counter++] = 246758;
        runs[counter++] = 246759;
        runs[counter++] = 246760;
        runs[counter++] = 246763;
        runs[counter++] = 246765;
        runs[counter++] = 246766;
        runs[counter++] = 246804;
        runs[counter++] = 246805;
        runs[counter++] = 246806;
        runs[counter++] = 246807;
        runs[counter++] = 246808;
        runs[counter++] = 246809;
        runs[counter++] = 246810;
        runs[counter++] = 246844;
        runs[counter++] = 246845;
        runs[counter++] = 246846;
        runs[counter++] = 246847;
        runs[counter++] = 246851;
        runs[counter++] = 246855;
        runs[counter++] = 246858;
        runs[counter++] = 246859;
        runs[counter++] = 246864;
        runs[counter++] = 246865;
        runs[counter++] = 246867;
        runs[counter++] = 246870;
        runs[counter++] = 246871;
        runs[counter++] = 246928;
        runs[counter++] = 246930;
        runs[counter++] = 246937;
        runs[counter++] = 246942;
        runs[counter++] = 246945;
        runs[counter++] = 246948;
        runs[counter++] = 246949;
        runs[counter++] = 246980;
        runs[counter++] = 246982;
        runs[counter++] = 246984;
        runs[counter++] = 246989;
        runs[counter++] = 246991;
        runs[counter++] = 246994;
      }
      else{
        runs[counter++] = 244918;
        if(pass.EqualTo("pass2_lowIR")) runs[counter++] = 244917;
        runs[counter++] = 244975;
        runs[counter++] = 244980;
        runs[counter++] = 244982;
        runs[counter++] = 244983;
        runs[counter++] = 245061;
        runs[counter++] = 245064;
        runs[counter++] = 245066;
        runs[counter++] = 245068; 
        runs[counter++] = 246390; 
        runs[counter++] = 246391; 
        runs[counter++] = 246392;
      }
    }
    else{//pp
      runs[counter++] = 244628;
      runs[counter++] = 244627;
      runs[counter++] = 244626;
      runs[counter++] = 244619;
      runs[counter++] = 244618;
      runs[counter++] = 244617;
      runs[counter++] = 244542;
      runs[counter++] = 244540;
      runs[counter++] = 244531;
      runs[counter++] = 244484;
      runs[counter++] = 244483;
      runs[counter++] = 244482;
      runs[counter++] = 244481;
      runs[counter++] = 244480;
      runs[counter++] = 244456;
      runs[counter++] = 244453;
      runs[counter++] = 244421;
      runs[counter++] = 244416;
      runs[counter++] = 244377;
      runs[counter++] = 244364;
      runs[counter++] = 244359;
      runs[counter++] = 244355;
      runs[counter++] = 244351;
      runs[counter++] = 244343;
      runs[counter++] = 244340;
      
    }
  }    
  
  if(counter != nruns){
    cout<<"Different number of runs has been used!! Aborting!!"<<endl;
    return 0x0;
  }
  
  for (Int_t i = 0; i < nruns; i++) for (Int_t j = 0; j < nruns; j++){
    if(i == j) continue;
    if(runs[i] == runs[j]){
      cout<<"run is present twice!! Aborting!!"<<endl;
      return 0x0;
    }
  }
  
  //Set on which run to run
  if(optHeavyIon){//Hi
    if(isMC){
      nrunssplit = 1;
      ntorun = 2;
      ntoskip = 2;
    }
    else{
      ntorun = 13;
      ntoskip = 0;
    }
  }
  else{//pp
    if(isMC){
      nrunssplit = 1;
      ntorun = 7;
      ntoskip = 0;
    }
    else{
      nrunssplit = 1;
      ntorun = 15;
      ntoskip = 0;
    }
  }
  
  if(fname.EqualTo("test")){
    ntorun = 1;
    ntoskip = 1;
  }
  
  const TString useonly = optTestMC < 0 ? "" : partialmc;
  
  if(useonly.EqualTo("")){//Standard way
    for(int i=ntoskip; i<ntorun+ntoskip && i<nruns; i++){
      if(i>=nruns){
        ::Fatal(" runs selection : ","Wrong definition of Start-Stop");
        break;
      }    
      cout<<"Adding run number "<<runs[i]<<" element "<<i+1<<" of "<<nruns<<endl;
      plugin->AddRunNumber(runs[i]);
      /////////////////////////////////////Log
      logFile<<Form("**(%i/%i) Run #%i#",i+1,nruns,runs[i])<<endl;
      ////////////////////////////////////////
    }
  }
  else{
    cout<<"Required custom list: "<<useonly<<endl;
    for(int i = 0; i < nruns; i++){
      cout<<"Try from custom list run number "<<runs[i]<<" element "<<i+1<<" of "<<nruns<<endl;
      if(!useonly.Contains(Form(" %i,", runs[i]))) continue;
      cout<<"Added !"<<endl;
      
      plugin->AddRunNumber(runs[i]);
      /////////////////////////////////////Log
      logFile<<Form("**(%i/%i) Run #%i#",i+1,nruns,runs[i])<<endl;
      ////////////////////////////////////////
    }
  }
  
  plugin->SetOutputToRunNo();
  plugin ->SetKeepLogs(kTRUE);
  plugin->SetNrunsPerMaster(nrunssplit);
  
  plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
  //   plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/ITS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/CONTAINERS  -I$ALICE_ROOT/STEER/STEER -I$ALICE_ROOT/STEER/STEERBase -I$ALICE_ROOT/STEER/ESD -I$ALICE_ROOT/STEER/CDB -I$ALICE_ROOT/TRD -I$ALICE_ROOT/macros -I$ALICE_ROOT/ANALYSIS -g");
  plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include -g");
  plugin->SetAnalysisSource("AliAnTOFtrack.cxx AliAnalysisTaskTOFSpectra.cxx"); 
  plugin->SetAdditionalLibs("AliAnalysisTaskTOFSpectra.h AliAnalysisTaskTOFSpectra.cxx UtilTOFParams.h AliAnTOFtrack.h AliAnTOFtrack.cxx libSTEERBase.so libESD.so libTOFrec.so libCDB.so  libRAWDatarec.so libRAWDatasim.so libTOFrec.so libANALYSIS.so libANALYSISalice.so");
  
  // Declare the output file names separated by blancs.
  // (can be like: file.root or file.root@ALICE::Niham::File)
  // To only save certain files, use SetDefaultOutputs(kFALSE), and then
  // SetOutputFiles("list.root other.filename") to choose which files to save
  plugin->SetDefaultOutputs(kFALSE);
  
  TString outputFileNames = "TListTOF.root";
  if(optTree){
    outputFileNames += " TreeTOF.root";
  }
  if(optMismatchrun) outputFileNames = "TListTOF_Mismatch.root";
  
  if(isMC){
    outputFileNames = "TListTOF_MC.root";
    
    if(optTree){
      outputFileNames += " TreeTOF_MC.root";
    }
    
    if(optMismatchrun) outputFileNames = "TListTOF_MC_Mismatch.root";
  }
  
  outputFileNames += " AnalysisResults.root";
  
  
  plugin->SetOutputFiles(outputFileNames);
  
  // Optionally define the files to be archived.
  //   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
  plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
  // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
  plugin->SetMergeViaJDL(optJDLMerge);
  if(optJDLMerge){
    plugin->SetOneStageMerging(kFALSE);
    //     plugin->SetMaxMergeFiles(50);
    plugin->SetMaxMergeStages(2);
  }
  //  plugin->SetMaxInitFailed(5);
  // Optionally resubmit threshold.
  plugin->SetMasterResubmitThreshold(90);
  // Optionally set time to live (default 30000 sec)
  plugin->SetTTL(30000);//was set to 20000
  // Optionally set input format (default xml-single)
  plugin->SetInputFormat("xml-single");
  // Optionally modify the name of the generated JDL (default analysis.jdl)
  // Optionally modify job price (default 1)
  plugin->SetPrice(1);      
  // Optionally modify split mode (default 'se')    
  plugin->SetSplitMode("se");
  delete runs;
  return (AliAnalysisGrid*) plugin;
}
