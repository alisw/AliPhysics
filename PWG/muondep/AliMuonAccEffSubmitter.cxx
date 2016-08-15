/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
// AliMuonAccEffSubmitter : a class to help submit Acc x Eff simulations
// anchored to real runs for J/psi, upsilon, single muons, etc...
//
// This class is dealing with 3 different directories :
//
// - template directory ($ALICE_ROOT/PWG/muondep/AccEffTemplates) containing the
//   basic template files to be used for a simuation. A template can contain
//   some variables that will be replaced during during the copy from template
//   to local dir
//
// - local directory, where the files from the template directory, are copied
//   once the class has been configured properly (i.e. using the various Set, Use,
//   etc... methods). Some other files (e.g. JDL ones) are generated from
//   scratch and also copied into this directory.
//   At this point one could(should) check the files, as they are the ones
//   to be copied to the remote directory for the production
//
// - remote directory, the alien directory where the files will be copied
//   (from the local directory) before the actual submission
//
// ==========================================================
//
// Basic usage
//
// AliMuonAccEffSubmitter a; // (1)
// a.UseOCDBSnapshots(kFALSE);
// a.SetRemoteDir("/alice/cern.ch/user/l/laphecet/Analysis/LHC13d/simjpsi/pp503z0");
// a.ShouldOverwriteFiles(true);
// a.MakeNofEventsPropToTriggerCount("CMUL7-B-NOPF-MUON");
// a.SetVar("VAR_GENLIB_PARNAME","\"pp 5.03\"");
// a.SetRunList(195682);
// a.Print();
// a.Run("test"); // will do everything but the submit
// a.Submit(false); // actual submission
//
// author: Laurent Aphecetche (Subatech)
//

#include "AliMuonAccEffSubmitter.h"

#include "AliAnalysisTriggerScalers.h"
#include "AliLog.h"
#include "TFile.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TMap.h"
#include "TMath.h"
#include "TObjString.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
using std::ifstream;
namespace
{
  Int_t splitLevel=10;
}

ClassImp(AliMuonAccEffSubmitter)

//______________________________________________________________________________
AliMuonAccEffSubmitter::AliMuonAccEffSubmitter(const char* generator, Bool_t localOnly,
                                               const char* generatorVersion)
: AliMuonGridSubmitter(AliMuonGridSubmitter::kAccEff,localOnly),
fRatio(-1.0),
fFixedNofEvents(10000),
fMaxEventsPerChunk(5000),
fOCDBPath(""),
fSplitMaxInputFileNumber(20),
fLogOutToKeep(""),
fRootOutToKeep(""),
fExternalConfig(""),
fUseOCDBSnapshots(kFALSE),
fSnapshotDir(""),
fUseAODMerging(kFALSE)
{
  // ctor
  //
  // if generator contains "pythia8" and generatorVersion is given then
  // the pythiaversion must represent the integer part XXX of the
  // include directory $ALICE_ROOT/PYTHI8/pythiaXXX/include where the file
  // Analysis.h is to be found.
  //
  // if generator contains "pythia6" then generatorVersion should be the
  // X.YY part of libpythia6.X.YY.so
  //
  
  SetCompactMode(1);

  AddIncludePath("-I$ALICE_ROOT/include");
  
  TString ocdbPath("raw://");
  
  if (localOnly) {
    ocdbPath = "local://$ALICE_ROOT/OCDB";
  }
  
  SetOCDBPath(ocdbPath.Data());
  
  SetLocalDirectory("Snapshot",LocalDir());
  
  SetVar("VAR_OCDB_PATH",Form("\"%s\"",ocdbPath.Data()));
  SetVar("VAR_AOD_MERGE_FILES","\"AliAOD.root,AliAOD.Muons.root\"");

  SetVar("VAR_GENPARAM_INCLUDE","AliGenMUONlib.h");
  SetVar("VAR_GENPARAM_NPART","1");
  SetVar("VAR_GENPARAM_GENLIB_TYPE","AliGenMUONlib::kJpsi");
  SetVar("VAR_GENPARAM_GENLIB_PARNAME","\"pPb 5.03\"");

  SetVar("VAR_GENCORRHF_QUARK","5");
  SetVar("VAR_GENCORRHF_ENERGY","5");

  // some default values for J/psi
  SetVar("VAR_GENPARAMCUSTOM_PDGPARTICLECODE","443");

  // default values below are from J/psi p+Pb (from muon_calo pass)
  SetVar("VAR_GENPARAMCUSTOM_Y_P0","4.08E5");
  SetVar("VAR_GENPARAMCUSTOM_Y_P1","7.1E4");
  
  SetVar("VAR_GENPARAMCUSTOM_PT_P0","1.13E9");
  SetVar("VAR_GENPARAMCUSTOM_PT_P1","18.05");
  SetVar("VAR_GENPARAMCUSTOM_PT_P2","2.05");
  SetVar("VAR_GENPARAMCUSTOM_PT_P3","3.34");

  // some default values for single muons
  SetVar("VAR_GENPARAMCUSTOMSINGLE_PTMIN","0.35");
  
  SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P0","4.05962");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P1","1.0");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P2","2.46187");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_PT_P3","2.08644");

  SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P0","0.729545");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P1","0.53837");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P2","0.141776");
  SetVar("VAR_GENPARAMCUSTOMSINGLE_Y_P3","0.0130173");

  // some default values for GenBox
  
  SetVar("VAR_GENMUBOX_PTMIN","0");
  SetVar("VAR_GENMUBOX_PTMAX","20");
  SetVar("VAR_GENMUBOX_YMIN","-4.1");
  SetVar("VAR_GENMUBOX_YMAX","-2.4");

  SetVar("VAR_PYTHIA8_CMS_ENERGY","8000");
  SetVar("VAR_PYTHIA6_CMS_ENERGY","8000");
  
  SetVar("VAR_PURELY_LOCAL",Form("%d",localOnly));

  SetVar("VAR_USE_RAW_ALIGN","1");

  SetVar("VAR_SIM_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Ideal\"");
  
  SetVar("VAR_REC_ALIGNDATA","\"alien://folder=/alice/simulation/2008/v4-15-Release/Residual\"");
  
  SetVar("VAR_USE_ITS_RECO","0");
  SetVar("VAR_USE_MC_VERTEX","1");
  SetVar("VAR_VERTEX_SIGMA_X","0.0025");
  SetVar("VAR_VERTEX_SIGMA_Y","0.0029");
  
  UseOCDBSnapshots(fUseOCDBSnapshots);

  SetVar("VAR_TRIGGER_CONFIGURATION","");
  
  SetVar("VAR_LHAPDF","liblhapdf");
  SetVar("VAR_MUONMCMODE","1");

  SetVar("VAR_PYTHIA8_SETENV","");
  SetVar("VAR_PYTHIA6_SETENV","");
  SetVar("VAR_NEEDS_PYTHIA6", "0");
  SetVar("VAR_NEEDS_PYTHIA8", "0");

  if ( TString(generator).Contains("pythia8",TString::kIgnoreCase) )
  {
    fMaxEventsPerChunk =  500; // 5000 is not reasonable with Pythia8 (and ITS+MUON...)
    
    SetCompactMode(2); // keep AOD as for the time being the filtering driven from AODtrain.C cannot
    // add SPD tracklets to muon AODs.
    
    SetVar("VAR_USE_ITS_RECO","1");

    SetupPythia8(generatorVersion);
  
//    TString p8env;
//    
//    p8env += Form("  gSystem->Setenv(\"PYTHIA8DATA\", gSystem->ExpandPathName(\"$ALICE_ROOT/PYTHIA8/pythia%s/xmldoc\"));\n",generatorVersion);
//    
//    p8env += "  gSystem->Setenv(\"LHAPDF\",gSystem->ExpandPathName(\"$ALICE_ROOT/LHAPDF\"));\n";
//    
//    p8env +=  "  gSystem->Setenv(\"LHAPATH\",gSystem->ExpandPathName(\"$ALICE_ROOT/LHAPDF/PDFsets\"));\n";
//    
//    SetVar("VAR_PYTHIA8_SETENV",p8env.Data());
//
//    SetVar("VAR_PYTHIA8_SETUP_STRINGS","\"\"");

    SetVar("VAR_TRIGGER_CONFIGURATION","p-p");
  }
  
  if ( TString(generator).Contains("pythia6",TString::kIgnoreCase) )
  {
    fMaxEventsPerChunk =  500; // 5000 is not reasonable with Pythia6 (and ITS+MUON...)

    SetCompactMode(2); // keep AOD as for the time being the filtering driven from AODtrain.C cannot
    // add SPD tracklets to muon AODs.

    SetupPythia6(generatorVersion);
//    TString p6env;
//    
//    p6env += Form("gSystem->Load(\"libpythia6_%s\");",generatorVersion);
//    
//    SetVar("VAR_PYTHIA6_SETENV",p6env.Data());

    SetVar("VAR_USE_ITS_RECO","1");
    
    SetVar("VAR_TRIGGER_CONFIGURATION","p-p");
  }

  SetVar("VAR_EVENTS_PER_JOB",Form("%i",fMaxEventsPerChunk));

  SetGenerator(generator);
  
  if (localOnly)
  {
    MakeNofEventsFixed(10);
  }
  else
  {
    MakeNofEventsPropToTriggerCount();
  }
  
  AddToTemplateFileList("CheckESD.C");
  AddToTemplateFileList("CheckAOD.C");
  AddToTemplateFileList("AODtrainsim.C");
//  AddToTemplateFileList("validation.sh");
  
  AddToTemplateFileList("Config.C");
  AddToTemplateFileList("rec.C");
  AddToTemplateFileList("sim.C");
  AddToTemplateFileList("simrun.sh");
  AddToTemplateFileList(RunJDLName().Data());
  
  UseExternalConfig(fExternalConfig);
}

//______________________________________________________________________________
AliMuonAccEffSubmitter::~AliMuonAccEffSubmitter()
{
  // dtor
}

///______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::Generate(const char* jdlname) const
{
  if ( TString(jdlname).Contains("merge",TString::kIgnoreCase) )
  {
    return GenerateMergeJDL(jdlname);
  }
  else
  {
    return GenerateRunJDL(jdlname);
  }
}

///______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::GenerateMergeJDL(const char* name) const
{
  /// Create the JDL for merging jobs
  /// FIXME: not checked !
  
  AliDebug(1,"");

  std::ostream* os = CreateJDLFile(name);
  
  if (!os)
  {
    return kFALSE;
  }
  
  Bool_t final = TString(name).Contains("merge",TString::kIgnoreCase);

  (*os) << "# Generated merging jdl (production mode)" << std::endl
  << "# $1 = run number" << std::endl
  << "# $2 = merging stage" << std::endl
  << "# Stage_<n>.xml made via: find <OutputDir> *Stage<n-1>/*root_archive.zip" << std::endl;

  OutputToJDL(*os,"Packages",GetMapValue("AliPhysics"));
  
  OutputToJDL(*os,"Executable","AOD_merge.sh");
  
  OutputToJDL(*os,"Price","1");

  if ( final )
  {
    OutputToJDL(*os,"Jobtag","comment: AliMuonAccEffSubmitter final merging");
  }
  else
  {
    OutputToJDL(*os,"Jobtag","comment: AliMuonAccEffSubmitter merging stage $2");
  }
  
  OutputToJDL(*os,"Workdirectorysize","5000MB");
  
  OutputToJDL(*os,"Validationcommand",Form("%s/validation_merge.sh",RemoteDir().Data()));
  
  OutputToJDL(*os,"TTL","14400");
  
  OutputToJDL(*os,"OutputArchive",
    "log_archive.zip:stderr,stdout@disk=1",
    "root_archive.zip:AliAOD.root,AliAOD.Muons.root,AnalysisResults.root@disk=3"
         );
  
  OutputToJDL(*os,"Arguments",(final ? "2":"1")); // for AOD_merge.sh, 1 means intermediate merging stage, 2 means final merging
  
  if ( !final )
  {
    OutputToJDL(*os,"InputFile",Form("LF:%s/AODtrainsim.C",RemoteDir().Data()));
    OutputToJDL(*os,"OutputDir",Form("%s/$1/Stage_$2/#alien_counter_03i#",RemoteDir().Data()));
    OutputToJDL(*os,"InputDataCollection",Form("%s/$1/Stage_$2.xml,nodownload",RemoteDir().Data()));
    OutputToJDL(*os,"split","se");
    OutputToJDL(*os,"SplitMaxInputFileNumber",GetSplitMaxInputFileNumber());
    OutputToJDL(*os,"InputDataListFormat","xml-single");
    OutputToJDL(*os,"InputDataList","wn.xml");
  }
  else
  {
    OutputToJDL(*os,"InputFile",Form("LF:%s/AODtrainsim.C",RemoteDir().Data()),
           Form("LF:%s/$1/wn.xml",RemoteDir().Data()));
    OutputToJDL(*os,"OutputDir",Form("%s/$1",RemoteDir().Data()));
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::GenerateRunJDL(const char* name) const
{
  /// Generate (locally) the JDL to perform the simulation+reco+aod filtering
  /// (to be then copied to the grid and finally submitted)
  
  AliDebug(1,"");

  std::ostream* os = CreateJDLFile(name);
  
  if (!os)
  {
    return kFALSE;
  }
  
  OutputToJDL(*os,"Packages",GetMapValue("Generator"),GetMapValue("AliPhysics"));
              
  OutputToJDL(*os,"Jobtag","comment: AliMuonAccEffSubmitter RUN $1");

  OutputToJDL(*os,"split","production:$2-$3");

  OutputToJDL(*os,"Price","1");
  
  OutputToJDL(*os,"OutputDir",Form("%s/$1/#alien_counter_03i#",RemoteDir().Data()));

  OutputToJDL(*os,"Executable","/alice/bin/aliroot_new");
  
  TObjArray files;
  files.SetOwner(kTRUE);
  TIter next(LocalFileList());
  TObjString* file;
  
  while ( ( file = static_cast<TObjString*>(next())) )
  {
    if ( !file->String().Contains("jdl",TString::kIgnoreCase) &&
         !file->String().Contains("OCDB_") )
    {
      files.Add(new TObjString(Form("LF:%s/%s",RemoteDir().Data(),file->String().Data())));
    }
  }
  
  if ( fUseOCDBSnapshots )
  {
    files.Add(new TObjString(Form("LF:%s/OCDB/$1/OCDB_sim.root",RemoteDir().Data())));
    files.Add(new TObjString(Form("LF:%s/OCDB/$1/OCDB_rec.root",RemoteDir().Data())));
  }
  
  OutputToJDL(*os,"InputFile",files);
  
  if ( fLogOutToKeep.IsNull() || fRootOutToKeep.IsNull() ) {
    AliError(Form("Output files not correctly set. Log: %s   Root: %s",fLogOutToKeep.Data(),fRootOutToKeep.Data()));
    delete os;
    return kFALSE;
  }
  else {
    OutputToJDL(*os,"OutputArchive",fLogOutToKeep.Data(),fRootOutToKeep.Data());
  }
  
  OutputToJDL(*os,"splitarguments","--run $1 --event #alien_counter# --eventsPerJob $4");
  
  OutputToJDL(*os,"Workdirectorysize","5000MB");
  
  OutputToJDL(*os,"JDLVariables","Packages","OutputDir");

  OutputToJDL(*os,"Validationcommand","/alice/validation/validation.sh");

  if ( GetVar("VAR_GENERATOR").Contains("pythia",TString::kIgnoreCase) )
  {
    OutputToJDL(*os,"TTL","36000");
  }
  else
  {
    OutputToJDL(*os,"TTL","14400");
  }
  
  return kTRUE;
}

//______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::MakeOCDBSnapshots()
{
  /// Run sim.C and rec.C in a special mode to generate OCDB snapshots
  /// Can only be done after the templates have been copied locally
  
  if (!IsValid()) return kFALSE;

  if (!fUseOCDBSnapshots) return kTRUE;
  
  if (!NofRuns()) return kFALSE;
  
  AliDebug(1,"");

  Bool_t ok(kTRUE);
  
  const std::vector<int>& runs = RunList();
  
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    Int_t runNumber = runs[i];

    TString ocdbSim(Form("%s/OCDB/%d/OCDB_sim.root",SnapshotDir().Data(),runNumber));
    TString ocdbRec(Form("%s/OCDB/%d/OCDB_rec.root",SnapshotDir().Data(),runNumber));

    if ( !gSystem->AccessPathName(ocdbSim.Data()) &&
         !gSystem->AccessPathName(ocdbRec.Data()) )
    {
      AliWarning(Form("Local OCDB snapshots already there for run %d. Will not redo them. If you want to force them, delete them by hand !",runNumber));
      continue;
    }
    else
    {
      gSystem->Exec(Form("simrun.sh --run %d --snapshot",runNumber));
    
      if ( gSystem->AccessPathName(ocdbSim.Data()) )
      {
        AliError(Form("Could not create OCDB snapshot for simulation"));
        ok = kFALSE;
      }

      if ( gSystem->AccessPathName(ocdbRec.Data()) )
      {
        AliError(Form("Could not create OCDB snapshot for reconstruction"));
        ok = kFALSE;
      }
    }
    
    AddToLocalFileList(ocdbSim);
    AddToLocalFileList(ocdbRec);
  }
  
  return ok;
}

//______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::Merge(Int_t stage, Bool_t dryRun)
{
  /// Submit multiple merging jobs with the format "submit AOD_merge(_final).jdl run# (stage#)".
  /// Also produce the xml collection before sending jobs
  /// Initial AODs will be taken from fRemoteDir/[RUNNUMBER] while the merged
  /// ones will be put into fMergedDir/AODs/[RUNNUMBER]
  ///
  /// Example:
  /// - inDir = "/alice/sim/2012/LHC12a10_bis" (where to find the data to merge)
  ///         = 0x0 --> inDir = homeDir/outDir/resDir
  /// - outDir = "Sim/LHC11h/embedding/AODs" (where to store merged results)
  /// - runList.txt must contains the list of run number
  /// - stage=0 --> final merging / stage>0 --> intermediate merging i
  ///
  
  if (!RemoteDirectoryExists(MergedDir().Data())) {
    AliError(Form("directory %s does not exist", MergedDir().Data()));
    return kFALSE;
  }
  
  gGrid->Cd(MergedDir().Data());
  
  TString jdl = MergeJDLName(stage==0);
  
  if (!RemoteFileExists(jdl.Data()))
  {
    AliError(Form("file %s does not exist in %s\n", jdl.Data(), RemoteDir().Data()));
    return kFALSE;
  }
  
  const std::vector<int>& runs = RunList();
  
  if (runs.empty())
  {
    AliError("No run to work with");
    return 0;
  }

  TString currRun;
  TString reply = "";
  gSystem->Exec("rm -f __failed__");
  Bool_t failedRun = kFALSE;
  
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    Int_t run = runs[i];
    AliInfo(Form("\n --- processing run %d ---\n", run));
    
    TString runDir = Form("%s/%d", MergedDir().Data(), run);
    
    if (!RemoteDirectoryExists(runDir.Data()))
    {
      AliInfo(Form(" - creating output directory %s\n", runDir.Data()));
      gSystem->Exec(Form("alien_mkdir -p %s", runDir.Data()));
    }
    
    if (RemoteFileExists(Form("%s/root_archive.zip", runDir.Data())))
    {
      AliWarning(" ! final merging already done");
      continue;
    }
    
    Int_t lastStage = GetLastStage(runDir.Data());
    
    if (stage > 0 && stage != lastStage+1)
    {
      AliError(Form(" ! lastest merging stage = %d. Next must be stage %d or final stage\n", lastStage, lastStage+1));
      continue;
    }
    
    TString wn = (stage > 0) ? Form("Stage_%d.xml", stage) : "wn.xml";
    TString find = (lastStage == 0) ?
    Form("alien_find -x %s %s/%d *root_archive.zip", wn.Data(), RemoteDir().Data(), run) :
    Form("alien_find -x %s %s/%d/Stage_%d *root_archive.zip", wn.Data(), RemoteDir().Data(), run, lastStage);
    gSystem->Exec(Form("%s 1> %s 2>/dev/null", find.Data(), wn.Data()));
    gSystem->Exec(Form("grep -c /event %s > __nfiles__", wn.Data()));
    ifstream f2("__nfiles__");
    TString nFiles;
    nFiles.ReadLine(f2,kTRUE);
    f2.close();
    gSystem->Exec("rm -f __nfiles__");
    printf(" - number of files to merge = %d\n", nFiles.Atoi());
    if (nFiles.Atoi() == 0) {
      printf(" ! collection of files to merge is empty\n");
      gSystem->Exec(Form("rm -f %s", wn.Data()));
      continue;
    } else if (stage > 0 && nFiles.Atoi() <= splitLevel && !reply.BeginsWith("y")) {
      if (!reply.BeginsWith("n")) {
        printf(" ! number of files to merge <= split level (%d). Continue? [Y/n] ", splitLevel);
        fflush(stdout);
        reply.Gets(stdin,kTRUE);
        reply.ToLower();
      }
      if (reply.BeginsWith("n")) {
        gSystem->Exec(Form("rm -f %s", wn.Data()));
        continue;
      } else reply = "y";
    }
    
    if (!dryRun)
    {
      TString dirwn = Form("%s/%s", runDir.Data(), wn.Data());
      if (RemoteFileExists(dirwn.Data())) gGrid->Rm(dirwn.Data());
      gSystem->Exec(Form("alien_cp file:%s alien://%s", wn.Data(), dirwn.Data()));
      gSystem->Exec(Form("rm -f %s", wn.Data()));
    }
    
    TString query;
    if (stage > 0) query = Form("submit %s %d %d", jdl.Data(), run, stage);
    else query = Form("submit %s %d", jdl.Data(), run);
    printf(" - %s ...", query.Data());
    fflush(stdout);
    
    if (dryRun)
    {
      AliInfo(" dry run");
      continue;
    }
    
    Bool_t done = kFALSE;
    TGridResult *res = gGrid->Command(query);
    if (res)
    {
      TString cjobId1 = res->GetKey(0,"jobId");
      if (!cjobId1.IsDec())
      {
        AliError(" FAILED");
        gGrid->Stdout();
        gGrid->Stderr();
      }
      else
      {
        AliInfo(Form(" DONE\n   --> the job Id is: %s \n", cjobId1.Data()));
        done = kTRUE;
      }
      delete res;
    }
    else
    {
      AliError(" FAILED");
    }
    
    if (!done)
    {
      gSystem->Exec(Form("echo %d >> __failed__", run));
      failedRun = kTRUE;
    }
    
  }
  
  if (failedRun)
  {
    AliInfo("\n--------------------\n");
    AliInfo("list of failed runs:\n");
    gSystem->Exec("cat __failed__");
    gSystem->Exec("rm -f __failed__");
    return kFALSE;
  }
  
  return kTRUE;
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::Print(Option_t* opt) const
{
  /// Printout
  
  AliMuonGridSubmitter::Print(opt);

  if ( fRatio > 0 )
  {
    std::cout << std::endl << Form("-- For each run, will generate %5.2f times the number of real events for trigger %s",
                      fRatio,ReferenceTrigger().Data()) << std::endl;
  }
  else
  {
    std::cout << std::endl <<  Form("-- For each run, will generate %10d events",fFixedNofEvents) << std::endl;
  }
  
  std::cout << "-- MaxEventsPerChunk = " << fMaxEventsPerChunk << std::endl;
  
  std::cout << "-- Will" << (fUseOCDBSnapshots ? "" : " NOT") << " use OCDB snaphosts" << std::endl;
}

//______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::Run(const char* mode)
{
  /// mode can be one of (case insensitive)
  ///
  /// LOCAL : copy the template files from the template directory to the local one
  /// UPLOAD : copy the local files to the grid (requires LOCAL)
  /// OCDB : make ocdb snapshots (requires LOCAL)
  /// SUBMIT : submit the jobs (requires LOCAL + UPLOAD)
  /// FULL : all of the above (requires all of the above)
  ///
  /// TEST : as SUBMIT, but in dry mode (does not actually submit the jobs)
  ///
  /// LOCALTEST : completely local test (including execution)
  
  if (!IsValid()) return kFALSE;
  
  if ( fRootOutToKeep.Contains("AliAOD.Muons.root") && ! fRootOutToKeep.Contains("AliAOD.root") ) SetVar("VAR_AOD_MERGE_FILES","\"AliAOD.Muons.root\"");

  TString smode(mode);
  smode.ToUpper();
  
  if ( smode == "FULL")
  {
    return  ( Run("LOCAL") && Run("OCDB") && Run("UPLOAD") && Run("SUBMIT") );
  }
  
  if ( smode == "LOCAL")
  {
    return CopyTemplateFilesToLocal();
  }
  
  if ( smode == "UPLOAD" )
  {
    return (CopyLocalFilesToRemote());
  }
  
  if ( smode == "OCDB" )
  {
    Bool_t ok = Run("LOCAL");
    if (ok)
    {
      ok = MakeOCDBSnapshots();
    }
    return ok;
  }
  
  if ( smode == "TEST" )
  {
    Bool_t ok = Run("LOCAL") && Run("OCDB") && Run("UPLOAD");
    if ( ok )
    {
      ok = (Submit(kTRUE)>0);
    }
    return ok;
  }
  
  if ( smode == "FULL" )
  {
    Bool_t ok = Run("LOCAL")  && Run("OCDB") && Run("UPLOAD");
    if ( ok )
    {
      ok = (Submit(kFALSE)>0);
    }
    return ok;
  }

  if( smode == "SUBMIT" )
  {
    return (Submit(kFALSE)>0);
  }
  
  if ( smode == "LOCALTEST" )
  {
    Bool_t ok = Run("LOCAL");
    if ( ok )
    {
      ok = LocalTest();
    }
    return ok;
  }
  
  return kFALSE;
}

//______________________________________________________________________________
Bool_t AliMuonAccEffSubmitter::SetGenerator(const char* generator)
{
  // set the variable to select the generator macro in Config.C

  Invalidate();
  
  TString generatorFile(Form("%s/%s.C",TemplateDir().Data(),generator));
  
  Int_t nofMissingVariables(0);
  
  // first check we indeed have such a macro
  if (!gSystem->AccessPathName(generatorFile.Data()))
  {
    TObjArray* variables = GetVariables(generatorFile.Data());
    
    TIter next(variables);
    TObjString* var;
    
    while ( ( var = static_cast<TObjString*>(next())) )
    {
      if ( !Vars()->GetValue(var->String()) )
      {
        ++nofMissingVariables;
        AliError(Form("file %s expect the variable %s to be defined, but we've not defined it !",generatorFile.Data(),var->String().Data()));
      }
    }
    
    delete variables;
    
    if ( !nofMissingVariables )
    {
      if (CheckCompilation(generatorFile.Data()))
      {
        Validate();
        SetVar("VAR_GENERATOR",Form("%s",generator));        
        AddToTemplateFileList(Form("%s.C",generator));
        return kTRUE;
      }
    }
    else
    {
      return kFALSE;
    }
  }
  else
  {
    AliError(Form("Can not work with the macro %s",generatorFile.Data()));
  }
  return kFALSE;
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::SetOCDBPath(const char* ocdbPath)
{
  /// Sets the OCDB path to be used
  
  SetMapKeyValue("OCDBPath",ocdbPath);
}


//______________________________________________________________________________
void AliMuonAccEffSubmitter::SetOCDBSnapshotDir(const char* dir)
{
  // change the directory used for snapshot
  
  if (gSystem->AccessPathName(Form("%s/OCDB",dir)))
  {
    AliError(Form("Snapshot top directory (%s) should contain an OCDB subdir with runnumbers in there",dir));
  }
  else
  {
    SetMapKeyValue("OCDBSnapshot",dir);
  }
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::MakeNofEventsPropToTriggerCount(const char* trigger, Float_t ratio)
{
  SetMapKeyValue("ReferenceTrigger",trigger);
  fRatio = ratio;
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::MakeNofEventsFixed(Int_t nevents)
{
  fFixedNofEvents = nevents;
  fRatio=0.0;
  SetMapKeyValue("ReferenceTrigger","");
}

//______________________________________________________________________________
Int_t AliMuonAccEffSubmitter::LocalTest()
{
  /// Generate a local macro (simrun.sh) to execute locally a full scale test
  /// Can only be used with a fixed number of events (and runnumber is fixed to zero)
  
  if ( fRatio > 0 )
  {
    AliError("Can only work in local test with a fixed number of events");
    return 0;
  }
  
  if ( fFixedNofEvents <= 0 )
  {
    AliError("Please fix the number of input events using MakeNofEventsFixed()");
    return 0;
  }
  
  const std::vector<int>& runs = RunList();

  if ( runs.empty() )
  {
    AliError("No run to work with");
    return 0;
  }
  
//  std::cout << "Generating script to execute : ./simrun.sh" << std::endl;
//
//  std::ofstream out("simrun.sh");
//
//  out << "#!/bin/bash" << std::endl;
////  root.exe -b -q simrun.C  --run <x> --chunk <y> --event <n>
//  out << "root.exe -b -q simrun.C --run "<< runs[0] <<" --event " << fFixedNofEvents << std::endl;
//
//  out.close();

  gSystem->Exec("chmod +x simrun.sh");
  gSystem->Exec("alien_cp alien:///alice/bin/aliroot_new file:");
  gSystem->Exec("chmod u+x aliroot_new");

  std::cout << "Cleaning up left-over files from previous simulation/reconstructions" << std::endl;

  gSystem->Exec("rm -rf TrackRefs.root *.SDigits*.root Kinematics.root *.Hits.root geometry.root gphysi.dat Run*.tag.root HLT*.root *.ps *.Digits.root *.RecPoints.root galice.root *QA*.root Trigger.root *.log AliESD* AliAOD* *.d *.so *.stat");

  TString command = Form("./aliroot_new --run %i --event 1 --eventsPerJob %i", runs[0], fFixedNofEvents);

  std::cout << "Executing the script : " << command.Data() << std::endl;


  gSystem->Exec(command.Data());
  
  return 1;
}

namespace  {

  void OutputRunList(const char* filename, const std::vector<int>& runlist)
  {
    /// output a runlist to ASCII file
    
    std::ofstream out(filename);

    for ( std::vector<int>::size_type j = 0; j < runlist.size(); ++j )
    {
      out << runlist[j] << std::endl;
    }
  }
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::SetCompactMode ( Int_t mode )
{
  /// Set the compact mode:
  /// 0 -> keep all root files + all logs
  /// 1 -> keep only AOD.Muons + all logs
  /// 2 -> keep only AODs and AOD.Muons + all logs
  /// 10 -> keep all root files + stout,stderr
  /// 11 -> keep only AOD.Muons + stout,stderr
  /// 12 -> keep only AODs and AOD.Muons + stout,stderr

  const char* allLogs = "stderr,stdout,*.log";
  const char* minLogs = "stderr,stdout";
  const char* allRoot = "galice*.root,Kinematics*.root,TrackRefs*.root,AliESDs.root,AliAOD.root,AliAOD.Muons.root,Merged.QA.Data.root,Run*.root";
  const char* muonAodRoot = "AliAOD.Muons.root,Merged.QA.Data.root";
  const char* aodRoot = "AliAOD.root,Merged.QA.Data.root";

  fLogOutToKeep = "";
  fRootOutToKeep = "";

  switch (mode) {
    case 0:
      fLogOutToKeep = allLogs;
      fRootOutToKeep = allRoot;
      break;
    case 1:
      fLogOutToKeep = allLogs;
      fRootOutToKeep = muonAodRoot;
      break;
    case 2:
      fLogOutToKeep = allLogs;
      fRootOutToKeep = aodRoot;
      break;
    case 10:
      fLogOutToKeep = minLogs;
      fRootOutToKeep = allRoot;
      break;
    case 11:
      fLogOutToKeep = minLogs;
      fRootOutToKeep = muonAodRoot;
      break;
    case 12:
      fLogOutToKeep = minLogs;
      fRootOutToKeep = aodRoot;
      break;
    default:
      AliError(Form("Unknown CompactMode %i",mode));
      break;
  }

  if ( ! fLogOutToKeep.IsNull() ) {
    fLogOutToKeep.Prepend("log_archive.zip:");
    fLogOutToKeep.Append("@disk=1");
  }
  if ( ! fRootOutToKeep.IsNull() ) {
    fRootOutToKeep.Prepend("root_archive.zip:");
    fRootOutToKeep.Append("@disk=2");
  }
}


///______________________________________________________________________________
void AliMuonAccEffSubmitter::SetupPythia6 ( const char *version )
{
  /// Setup pythia 6
  SetVar("VAR_NEEDS_PYTHIA6", "1");

  TString p6env = Form("gSystem->Load(\"libpythia6_%s\");",version);
  SetVar("VAR_PYTHIA6_SETENV",p6env.Data());
}

///______________________________________________________________________________
void AliMuonAccEffSubmitter::SetupPythia8 ( const char *version, const char* configStrings )
{
  /// Setup pythia 6
  SetVar("VAR_NEEDS_PYTHIA8", "1");

  TString p8env = Form("  gSystem->Setenv(\"PYTHIA8DATA\", gSystem->ExpandPathName(\"$ALICE_ROOT/PYTHIA8/pythia%s/xmldoc\"));\n",version);
  SetVar("VAR_PYTHIA8_SETENV",p8env.Data());

  SetVar("VAR_PYTHIA8_SETUP_STRINGS",Form("\"%s\"",configStrings));
}

//______________________________________________________________________________
Int_t AliMuonAccEffSubmitter::SplitRunList(const char* inputList, int maxJobs)
{
  /// In order to be able to submit, split a given runlist into chunks that will
  /// fit within maxJobs (1500 for a typical user)

  std::vector<int> runs;
  
  AliAnalysisTriggerScalers tmp(inputList);
  runs = tmp.GetRunList();
  
  AliAnalysisTriggerScalers* ts(0x0);
  std::vector<int> currentRunList;
  
  int nJobs(0);
  int nTotalJobs(0);
  int nEvts(0);
  int nFiles(0);
  
  for (std::vector<int>::size_type i=0; i < runs.size(); ++i)
  {
    Int_t runNumber = runs[i];
  
    Int_t nEvtRun(fFixedNofEvents);
    
    if ( fRatio > 0 )
    {
      if (!ts)
      {
        AliInfo(Form("Creating AliAnalysisTriggerScalers from OCDB=%s",OCDBPath().Data()));
        ts = new AliAnalysisTriggerScalers(runs,OCDBPath().Data());
      }
      
      AliAnalysisTriggerScalerItem* trigger = ts->GetTriggerScaler(runNumber, "L2A", ReferenceTrigger().Data());
      
      if (!trigger)
      {
        AliError(Form("Could not get trigger %s for run %09d",ReferenceTrigger().Data(),runNumber));
        continue;
      }
      nEvtRun = TMath::Nint(fRatio * trigger->Value());
    }
    
    Int_t nChunk = 1;
    
    while (nEvtRun/nChunk+0.5 > MaxEventsPerChunk())
    {
      ++nChunk;
    }
    
    Int_t nEvtChunk = TMath::Nint(nEvtRun/nChunk + 0.5);
    
    nJobs += nChunk;
    
    nTotalJobs += nChunk;
    
    nEvts += nChunk*nEvtChunk;

    if ( nJobs > maxJobs )
    {
      ++nFiles;
      
      OutputRunList(Form("%s.%d",inputList,nFiles),currentRunList);
      nJobs = 0;
      currentRunList.clear();
    }
    
    
    currentRunList.push_back(runNumber);
    
  }
  
  if ( !currentRunList.empty() )
  {
    ++nFiles;
    OutputRunList(Form("%s.%d",inputList,nFiles),currentRunList);

  }
  
  delete ts;
  
  std::cout << Form("input run list was split into %d files. Total number of jobs %d. Total number of events %d",
                    nFiles,nTotalJobs,nEvts) << std::endl;
  
  return nFiles;
}


//______________________________________________________________________________
Int_t AliMuonAccEffSubmitter::Submit(Bool_t dryRun)
{
  /// Submit multiple production jobs with the format "submit jdl 000run#.xml 000run#".
  ///
  /// Return the number of submitted (master) jobs
  ///
  /// Example:
  /// - outputDir = "/alice/cern.ch/user/p/ppillot/Sim/LHC10h/JPsiPbPb276/AlignRawVtxRaw/ESDs"
  /// - runList must contains the list of run number
  /// - trigger is the (fully qualified) trigger name used to compute the base number of events
  /// - mult is the factor to apply to the number of trigger to get the number of events to be generated
  ///   (# generated events = # triggers x mult
  
  if (!IsValid()) return 0;
  
  AliDebug(1,"");

  gGrid->Cd(RemoteDir());
  
  if (!RemoteFileExists(RunJDLName()))
  {
    AliError(Form("file %s does not exist in %s", RunJDLName().Data(), RemoteDir().Data()));
    return 0;
  }
  
  if ( !NofRuns() )
  {
    AliError("No run list set. Use SetRunList");
    return 0;
  }
  const std::vector<int>& runs = RunList();
  
  if (runs.empty())
  {
    AliError("No run to work with");
    return 0;
  }
  
  //  cout << "total number of selected MB events = " << totEvt << endl;
  //  cout << "required number of generated events = " << nGenEvents << endl;
  //  cout << "number of generated events per MB event = " << ratio << endl;
  //  cout << endl;
  
  std::cout << "run\tfirstChunk\tlastChunk\teventsPerJob" << std::endl;
  std::cout << "----------------------" << std::endl;
  
  Int_t nJobs(0);
  Int_t nEvts(0);
  
  AliAnalysisTriggerScalers* ts(0x0);
  
  for (std::vector<int>::size_type i=0; i < runs.size(); ++i)
  {
    Int_t runNumber = runs[i];
    
    Int_t nEvtRun(fFixedNofEvents);
    
    if ( fRatio > 0 )
    {
      if (!ts)
      {
        AliInfo(Form("Creating AliAnalysisTriggerScalers from OCDB=%s",OCDBPath().Data()));
        ts = new AliAnalysisTriggerScalers(runs,OCDBPath().Data());
      }
      
      AliAnalysisTriggerScalerItem* trigger = ts->GetTriggerScaler(runNumber, "L2A", ReferenceTrigger().Data());
    
      if (!trigger)
      {
        AliError(Form("Could not get trigger %s for run %09d",ReferenceTrigger().Data(),runNumber));
        continue;
      }
      nEvtRun = TMath::Nint(fRatio * trigger->Value());
    }

    Int_t nChunk = nEvtRun/MaxEventsPerChunk();
    if ( nChunk == 0 ) nChunk++;
    Int_t nEvtChunk = 0, delta = MaxEventsPerChunk();
    for ( Int_t tmpnChunk=nChunk; tmpnChunk<=nChunk+1; tmpnChunk++ ) {
      Int_t tmpnEventChunk = TMath::Min(nEvtRun/tmpnChunk,MaxEventsPerChunk());
      Int_t tmpDelta = TMath::Abs(tmpnChunk*tmpnEventChunk-nEvtRun);
      if ( tmpDelta < delta ) {
        nChunk = tmpnChunk;
        nEvtChunk = tmpnEventChunk;
        delta = tmpDelta;
      }
    }

    nJobs += nChunk;
    
    nEvts += nChunk*nEvtChunk;
    
    std::cout << runNumber << "\t1\t" << nChunk << "\t" << nEvtChunk << std::endl;
    
    TString query(Form("submit %s %d 1 %d %d", RunJDLName().Data(), runNumber, nChunk, nEvtChunk));
    
    std::cout << query.Data() << " ..." << std::flush;
    
    TGridResult* res = 0x0;
    
    if (!dryRun)
    {
      res = gGrid->Command(query);
    }
    
    if (res)
    {
      TString cjobId1 = res->GetKey(0,"jobId");
      
      if (!cjobId1.Length())
      {
        std::cout << " FAILED" << std::endl << std::endl;
        gGrid->Stdout();
        gGrid->Stderr();
      }
      else
      {
        std::cout << "DONE" << std::endl;
        std::cout << Form("   --> the job Id is: %s",cjobId1.Data()) << std::endl << std::endl;
      }
    }
    else
    {
      std::cout << " FAILED" << std::endl << std::endl;
    }
    
    delete res;
  }
  
  std::cout << std::endl
  << "total number of jobs = " << nJobs << std::endl
  << "total number of generated events = " << nEvts << std::endl
  << std::endl;
  
  delete ts;
  
  return nJobs;
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::UpdateLocalFileList(Bool_t clearSnapshots)
{
  /// Update the list of local files
  
  AliMuonGridSubmitter::UpdateLocalFileList();
  
  if (!NofRuns()) return;
  
  if ( clearSnapshots )
  {
    TIter next(LocalFileList());
    TObjString* file;
    
    while ( ( file = static_cast<TObjString*>(next())) )
    {
      if ( file->String().Contains("OCDB_") )
      {
        LocalFileList()->Remove(file);
      }
    }
    LocalFileList()->Compress();
  }

  const char* type[] = { "sim","rec" };
  
  const std::vector<int>& runs = RunList();
  
  for ( std::vector<int>::size_type i = 0; i < runs.size(); ++i )
  {
    Int_t runNumber = runs[i];
    
    for ( Int_t t = 0; t < 2; ++t )
    {
      TString snapshot(Form("%s/OCDB/%d/OCDB_%s.root",SnapshotDir().Data(),runNumber,type[t]));
      
      if ( !gSystem->AccessPathName(snapshot.Data()) )
      {
        AddToLocalFileList(snapshot);
      }
    }
  }
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::UseOCDBSnapshots(Bool_t flag)
{
  /// Whether or not to use OCDB snapshots
  /// Using OCDB snapshots will speed-up both the sim and reco initialization
  /// phases on each worker node, but takes time to produce...
  /// So using them is not always a win-win...
  
  fUseOCDBSnapshots = flag;
  if ( flag )
  {
    SetVar("VAR_OCDB_SNAPSHOT","kTRUE");
    
    // for some reason must include ITS objects in the snapshot
    // (to be able to instantiante the vertexer later on ?)
    
    SetVar("VAR_USE_ITS_RECO","1");
  }
  else
  {
    SetVar("VAR_OCDB_SNAPSHOT","kFALSE");
  }
  
  UpdateLocalFileList();
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::UseAODMerging(Bool_t flag)
{
  /// whether or not we should generate JDL for merging AODs
  
  fUseAODMerging = flag;
  
  AddToTemplateFileList(MergeJDLName(kFALSE).Data());
  AddToTemplateFileList(MergeJDLName(kTRUE).Data());
  AddToTemplateFileList("AOD_merge.sh");
  AddToTemplateFileList("validation_merge.sh");
}

//______________________________________________________________________________
void AliMuonAccEffSubmitter::UseExternalConfig(const char* externalConfigFullFilePath)
{
  // use an external config (or the default Config.C if externalConfigFullFilePath="")
  
  fExternalConfig = externalConfigFullFilePath;
  if ( fExternalConfig.Length() > 0 )
  {
    AddToTemplateFileList(fExternalConfig);
  }
  else
  {
    AddToTemplateFileList("Config.C");
  }
}
