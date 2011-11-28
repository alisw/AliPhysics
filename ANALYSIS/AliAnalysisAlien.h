#ifndef ALIANALYSISALIEN_H
#define ALIANALYSISALIEN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Mihaela Gheata, 01/09/2008

//==============================================================================
//   AliAnalysisAlien - AliEn utility class. Provides interface for creating
// a personalized JDL, finding and creating a dataset.
//==============================================================================

#ifndef ALIANALYSISGRID_H
#include "AliAnalysisGrid.h"
#endif

#ifndef ROOT_TString
#include <TString.h>
#endif

#ifndef ROOT_TMath
#include <TMath.h>
#endif

#ifndef ROOT_TMap
#include <TMap.h>
#endif

class AliAnalysisManager;
class AliAnalysisTaskCfg;
class TGridJDL;

class AliAnalysisAlien : public AliAnalysisGrid {

public:

   AliAnalysisAlien();
   AliAnalysisAlien(const char *name);
   virtual ~AliAnalysisAlien();
   AliAnalysisAlien(const AliAnalysisAlien& other); 
   AliAnalysisAlien& operator=(const AliAnalysisAlien& other);
// Setters   
   virtual void        AddAdditionalLibrary(const char *name);
   virtual void        AddIncludePath(const char *path);
   virtual void        AddRunNumber(Int_t run);
   virtual void        AddRunNumber(const char *run);
   virtual void        AddRunList(const char *runList);
   virtual void        AddDataFile(const char *lfn);
   virtual void        AddExternalPackage(const char *name);
   virtual void        SetExecutable(const char *name="analysis.sh")     {fExecutable = name;}
   virtual void        SetExecutableCommand(const char *command="root -b -q") {fExecutableCommand = command;}
   virtual void        SetArguments(const char *name="")                 {fArguments = name;}
   virtual void        SetExecutableArgs(const char *name="")            {fExecutableArgs = name;}
   virtual void        SetAnalysisMacro(const char *name="myAnalysis.C") {fAnalysisMacro = name;}
   virtual void        SetAnalysisSource(const char *name="myAnalysisClass.cxx") {fAnalysisSource = name;}
   virtual void        SetValidationScript(const char *name="validation.sh") {fValidationScript = name;}
   virtual void        SetAdditionalLibs(const char *list)               {fAdditionalLibs = list;}
   virtual void        SetAdditionalRootLibs(const char *list)           {fAdditionalRootLibs = list;}
   virtual void        SetPrice(Int_t price=1)                           {fPrice = price;}
   virtual void        SetRunRange(Int_t min, Int_t max)                 {fRunRange[0] = min; fRunRange[1] = max;}
   virtual void        SetJobTag(const char *tag="")                     {fJobTag = tag;}
   virtual void        SetNrunsPerMaster(Int_t nruns=1)                  {fNrunsPerMaster = nruns;}
   virtual void        SetMaxMergeFiles(Int_t nfiles)                    {fMaxMergeFiles = nfiles;}
   virtual void        SetMaxMergeStages(Int_t nstages)                  {fMaxMergeStages = nstages;}
   virtual void        SetSplitMode(const char *type="se")               {fSplitMode = type;}
   virtual void        SetSplitMaxInputFileNumber(Int_t nfiles=100)      {fSplitMaxInputFileNumber = nfiles;}
   virtual void        SetAPIVersion(const char *version)                {fAPIVersion = version;}
   virtual void        SetROOTVersion(const char *version)               {fROOTVersion = version;}
   virtual void        SetAliROOTVersion(const char *version)            {fAliROOTVersion=version;}
   virtual void        SetUser(const char *user)                         {fUser = user;}
   virtual void        SetTTL(Int_t ttl=30000)                           {fTTL = ttl;}
   virtual void        SetGridWorkingDir(const char *name="workdir")     {fGridWorkingDir = name;}
   virtual void        SetGridDataDir(const char *name)                  {fGridDataDir = name;}
   virtual void        SetDataPattern(const char *pattern="*AliESDs.root") {fDataPattern = pattern;}
   virtual void        SetFriendChainName(const char *name="")           {fFriendChainName = name;}
   virtual void        SetDefaultOutputs(Bool_t flag);
   virtual void        SetGridOutputDir(const char *name="output")       {fGridOutputDir = name;}
   virtual void        SetOutputArchive(const char *list="log_archive.zip:std*@disk=1 root_archive.zip:*.root@disk=2");
   virtual void        SetOutputFiles(const char *list);
   virtual void        SetOutputToRunNo(Int_t mode=1)                    {fOutputToRunNo = mode;}
   virtual void        SetInputFormat(const char *format="xml-single")   {fInputFormat = format;}
   virtual void        SetMaxInitFailed(Int_t nfail=5)                   {fMaxInitFailed = nfail;}
   virtual void        SetTerminateFiles(const char *list)               {fTerminateFiles = list;}
   virtual void        SetMergeExcludes(const char *list)                {fMergeExcludes = list;};
   virtual void        SetMergeViaJDL(Bool_t on=kTRUE)                   {fMergeViaJDL = on ? 1 : 0;}
   virtual void        SetMergeDirName(const char *name)                 {fMergeDirName = name;}
   virtual void        SetMasterResubmitThreshold(Int_t percentage)      {fMasterResubmitThreshold = percentage;}
   virtual void        SetNtestFiles(Int_t nfiles)                       {fNtestFiles = nfiles;}
   virtual void        SetNumberOfReplicas(Int_t ncopies)                {fNreplicas = TMath::Min(ncopies,4);}
   virtual void        SetJDLName(const char *name="analysis.jdl")       {fJDLName = name;}
   virtual void        SetPreferedSE(const char *se);
   virtual void        SetProductionMode(Int_t mode=1)                   {fProductionMode = mode;}
   virtual void        SetRegisterExcludes(const char *list)             {fRegisterExcludes = list;}
   virtual void        SetRunPrefix(const char *prefix);
   virtual void        SetOutputSingleFolder(const char *folder)         {fOutputSingle = folder; fSplitMode="file"; fSplitMaxInputFileNumber=1;}
   virtual void        SetFastReadOption(Bool_t on=kTRUE)                {fFastReadOption = on ? 1 : 0;}
   virtual void        SetOverwriteMode(Bool_t on=kTRUE)                 {fOverwriteMode = on ? 1 : 0;}

   TGridJDL           *GetGridJDL() const {return fGridJDL;}
   TGridJDL           *GetMergingJDL() const {return fMergingJDL;}
   const char         *GetGridOutputDir() const                          {return fGridOutputDir;}
//Utilities
   void                AddModule(AliAnalysisTaskCfg *module);
   void                AddModules(TObjArray *list);
   AliAnalysisManager *CreateAnalysisManager(const char *name, const char *filename="");
   Int_t               GetNmodules() const;
   AliAnalysisTaskCfg *GetModule(const char *name);
   Bool_t              LoadModules();
   Bool_t              GenerateTest(const char *name, const char *modname="");
   Bool_t              GenerateTrain(const char *name);
   virtual Bool_t      CreateDataset(const char *pattern);
   Bool_t              CopyLocalDataset(const char *griddir, const char *pattern, Int_t nfiles, const char *output="data.txt", const char *archivefile="", const char *outputdir="data");
   virtual Bool_t      CreateJDL();
   virtual void        EnablePackage(const char *package);
   static Bool_t       DirectoryExists(const char *lfn);
   static Bool_t       FileExists(const char *lfn);
   static const char  *GetJobStatus(Int_t jobidstart, Int_t lastid, Int_t &nrunning, Int_t &nwaiting, Int_t &nerror, Int_t &ndone);
   const char         *GetListOfFiles(const char *type);
   Bool_t              CheckMergedFiles(const char *filename, const char *aliendir, Int_t nperchunk, const char *jdl="");
   static AliAnalysisManager *LoadAnalysisManager(const char *fname);
   static Bool_t       MergeOutput(const char *output, const char *basedir, Int_t nmaxmerge, Int_t stage=0);
   virtual Bool_t      MergeOutputs();
   virtual void        Print(Option_t *option="") const;
   virtual Bool_t      StartAnalysis(Long64_t nentries=123456789, Long64_t firstentry=0);
   static Bool_t       SetupPar(const char *package);
   virtual Bool_t      Submit();
   virtual Bool_t      SubmitMerging();
   static Int_t        SubmitSingleJob(const char *query);
   virtual void        WriteAnalysisFile();
   virtual void        WriteAnalysisMacro();
   virtual void        WriteMergingMacro();
   virtual void        WriteMergeExecutable();
   virtual void        WriteExecutable();
   virtual Bool_t      WriteJDL(Bool_t copy);
   virtual void        WriteProductionFile(const char *filename) const;
   virtual void        WriteValidationScript(Bool_t merge=kFALSE);

// PROOF mode
   virtual void        SetProofCluster(const char *cluster)              {fProofCluster = cluster;}
   virtual void        SetProofDataSet(const char *dataset)              {fProofDataSet = dataset;}
   virtual const char *GetProofDataSet() const                           {return fProofDataSet.Data();}
   virtual void        SetProofParameter(const char *pname, const char *value);
   const char         *GetProofParameter(const char *pname) const;
   virtual void        SetProofReset(Int_t mode)                         {fProofReset = mode;}
   virtual void        SetNproofWorkers(Int_t nworkers)                  {fNproofWorkers = nworkers;}
   virtual void        SetNproofWorkersPerSlave(Int_t nworkers)          {fNproofWorkersPerSlave = nworkers;}
   virtual void        SetRootVersionForProof(const char *version)       {fRootVersionForProof = version;}
   virtual void        SetAliRootMode(const char *mode)                  {fAliRootMode = mode;}
   // .txt file containing the list of files to be chained in test mode
   virtual void        SetFileForTestMode(const char *filename)          {fFileForTestMode = filename;}
   virtual TChain     *GetChainForTestMode(const char *treeName) const;

protected:
   void                CdWork();
   Bool_t              CheckInputData();
   void                CheckDataType(const char *lfn, Bool_t &is_collection, Bool_t &is_xml, Bool_t &use_tags);
   virtual Bool_t      Connect();
   virtual void        SetDefaults();  
   Bool_t              SubmitNext();

   Bool_t              IsCollection(const char *lfn) const;
   virtual Bool_t      IsSingleOutput() const;
   Bool_t              IsUsingTags() const {return TObject::TestBit(AliAnalysisGrid::kUseTags);}
   Bool_t              LoadModule(AliAnalysisTaskCfg *mod);
   Bool_t              CheckDependencies();
   Bool_t              CheckFileCopy(const char *alienpath);

private:
   TGridJDL        *fGridJDL;         //! JDL maker
   TGridJDL        *fMergingJDL;      //! JDL maker
   Int_t            fPrice;           // Grid price for the job;
   Int_t            fTTL;             // Time to live.
   Int_t            fSplitMaxInputFileNumber; // Maximum number of files to be processed per subjob
   Int_t            fMaxInitFailed;   // Maximum initial consecutive subjobs accepted to fail
   Int_t            fMasterResubmitThreshold; // Failed jobs will be resubmitted until this DONE ratio
   Int_t            fNtestFiles;      // Number of files used in the testing case
   Int_t            fRunRange[2];     // Run range
   Int_t            fNrunsPerMaster;  // Number of runs per masterjob
   Int_t            fMaxMergeFiles;   // Maximum number of files to be merged in one chunk
   Int_t            fMaxMergeStages;  // Maximum number of merging stages
   Int_t            fNsubmitted;      // Number of jobs submitted
   Int_t            fProductionMode;  // Production mode (0-off, 1-on)
   Int_t            fOutputToRunNo;   // Use run number as output directory
   Int_t            fMergeViaJDL;     // Enable merging via automatic JDL
   Int_t            fFastReadOption;  // Use xrootd tweaks to reduce timeouts in file access
   Int_t            fOverwriteMode;   // Overwrite existing files if any
   Int_t            fNreplicas;       // Number of replicas for the output files
   Int_t            fNproofWorkers;   // Number of workers in proof mode
   Int_t            fNproofWorkersPerSlave; // Max number of workers per slave in proof mode
   Int_t            fProofReset;      // Proof reset mode: 0=no reset, 1=soft, 2=hard
   TString          fRunNumbers;      // List of runs to be processed
   TString          fExecutable;      // Executable script for AliEn job
   TString          fExecutableCommand;  // Command(s) to be executed in the executable script
   TString          fArguments;       // Arguments for the executable script
   TString          fExecutableArgs;  // arguments added to the executable script after the analysis macro
   TString          fAnalysisMacro;   // Root macro steering the analysis
   TString          fAnalysisSource;  // User analysis implementation (.cxx) file(s)
   TString          fValidationScript; // Name of the validation script
   TString          fAdditionalRootLibs;  // List (separated by blacs) of additional libraries needed for/before analysis libs/par file compilation
   TString          fAdditionalLibs;  // List (separated by blacs) of additional libraries needed for the analysis loaded AFTER all par files
   TString          fSplitMode;       // Job split mode
   TString          fAPIVersion;      // API version
   TString          fROOTVersion;     // ROOT version
   TString          fAliROOTVersion;  // AliROOT version
   TString          fExternalPackages; // External packages
   TString          fUser;            // AliEn user name
   TString          fGridWorkingDir;  // AliEn directory containing the input packages
   TString          fGridDataDir;     // AliEn data production directory
   TString          fDataPattern;     // Data pattern for 'find' command
   TString          fGridOutputDir;   // AliEn directory (wrt work dir) where the output should be written
   TString          fOutputArchive;   // List of output archives separated by blancs
   TString          fOutputFiles;     // List of output files separated by blancs
   TString          fInputFormat;     // Input format (xml-single)
   TString          fDatasetName;     // Dataset xml file to be created
   TString          fJDLName;         // JDL file to be generated
   TString          fTerminateFiles;  // List of output files produced during Terminate
   TString          fMergeExcludes;   // List of output files excluded from merging
   TString          fRegisterExcludes; // List of liles not to be registered/merged
   TString          fIncludePath;     // Include path
   TString          fCloseSE;         // Preffered storage element. Taken from alien_CLOSE_SE environment.
   TString          fFriendChainName; // File name to construct friend chain (for AOD)
   TString          fJobTag;          // Job tag
   TString          fOutputSingle;    // Directory name for the output when split is per file
   TString          fRunPrefix;       // Run prefix to be applied to run numbers
   TString          fProofCluster;    // Proof cluster name
   TString          fProofDataSet;    // Proof dataset to be used
   TString          fFileForTestMode; // .txt file for the chain to be used in PROOF test mode
   TString          fRootVersionForProof; // ROOT version to be used in PROOF mode. The default one taken if empty.
   TString          fAliRootMode;     // AliRoot mode among the list supported by the proof cluster
   TString          fMergeDirName;    // Name of the directory that should be added to the output directory
   TObjArray       *fInputFiles;      // List of input files to be processed by the job
   TObjArray       *fPackages;        // List of packages to be used
   TObjArray       *fModules;         // List of AliAnalysisTaskCfg modules
   TMap             fProofParam;      // Key-value pairs for proof mode
   
   ClassDef(AliAnalysisAlien, 20)   // Class providing some AliEn utilities
};
#endif
