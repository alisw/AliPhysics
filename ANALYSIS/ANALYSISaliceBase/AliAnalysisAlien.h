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
   const char         *GetExecutable() const                             {return fExecutable;}
   virtual void        SetExecutableCommand(const char *command="root -b -q") {fExecutableCommand = command;}
   const char         *GetExecutableCommand() const                      {return fExecutableCommand;}
   virtual void        SetArguments(const char *name="")                 {fArguments = name;}
   const char         *GetArguments() const                              {return fArguments;}
   virtual void        SetExecutableArgs(const char *name="")            {fExecutableArgs = name;}
   const char         *GetExecutableArgs() const                         {return fExecutableArgs;}
   virtual void        SetAnalysisMacro(const char *name="myAnalysis.C") {fAnalysisMacro = name;}
   const char         *GetAnalysisMacro() const                          {return fAnalysisMacro;}
   virtual void        SetAnalysisSource(const char *name="myAnalysisClass.cxx") {fAnalysisSource = name;}
   const char         *GetAnalysisSource() const                         {return fAnalysisSource;}
   virtual void        SetValidationScript(const char *name="validation.sh") {fValidationScript = name;}
   const char         *GetValidationScript() const                       {return fValidationScript;}
   virtual void        SetAdditionalLibs(const char *list)               {fAdditionalLibs = list;}
   const char         *GetAdditionalLibs() const                         {return fAdditionalLibs;}
   virtual void        SetAdditionalRootLibs(const char *list)           {fAdditionalRootLibs = list;}
   const char         *GetAdditionalRootLibs() const                     {return fAdditionalRootLibs;}
   virtual void        SetPrice(Int_t price=1)                           {fPrice = price;}
   Int_t               GetPrice() const                                  {return fPrice;}
   virtual void        SetRunRange(Int_t min, Int_t max)                 {fRunRange[0] = min; fRunRange[1] = max;}
   void                GetRunRange(Int_t &min, Int_t &max)               {min = fRunRange[0]; max = fRunRange[1];}
   virtual void        SetJobTag(const char *tag="")                     {fJobTag = tag;}
   const char         *GetJobTag() const                                 {return fJobTag;}
   virtual void        SetNrunsPerMaster(Int_t nruns=1)                  {fNrunsPerMaster = nruns;}
   Int_t               GetNrunsPerMaster() const                         {return fNrunsPerMaster;}
   virtual void        SetMaxMergeFiles(Int_t nfiles)                    {fMaxMergeFiles = nfiles;}
   Int_t               GetMaxMergeFiles()  const                         {return fMaxMergeFiles;}
   virtual void        SetMaxMergeStages(Int_t nstages)                  {fMaxMergeStages = nstages;}
   Int_t               GetMaxMergeStages() const                         {return fMaxMergeStages;}
   virtual void        SetSplitMode(const char *type="se")               {fSplitMode = type;}
   const char         *GetSplitMode() const                              {return fSplitMode;}
   virtual void        SetSplitMaxInputFileNumber(Int_t nfiles=100)      {fSplitMaxInputFileNumber = nfiles;}
   Int_t               GetSplitMaxInputFileNumber() const                {return fSplitMaxInputFileNumber;}
   virtual void        SetAPIVersion(const char *version)                {fAPIVersion = version;}
   const char         *GetAPIVersion() const                             {return fAPIVersion;}
   virtual void        SetROOTVersion(const char *version)               {fROOTVersion = version;}
   const char         *GetROOTVersion() const                            {return fROOTVersion;}
   virtual void        SetAliROOTVersion(const char *version)            {fAliROOTVersion=version;}
   const char         *GetAliROOTVersion() const                         {return fAliROOTVersion;}
   virtual void        SetAliPhysicsVersion(const char *version)         {fAliPhysicsVersion=version;}
   const char         *GetAliPhysicsVersion() const                      {return fAliPhysicsVersion;}
   virtual void        SetUser(const char *user)                         {fUser = user;}
   const char         *GetUser() const                                   {return fUser;}
   virtual void        SetTTL(Int_t ttl=30000)                           {fTTL = ttl;}
   Int_t               GetTTL() const                                    {return fTTL;}
   virtual void        SetGridWorkingDir(const char *name="workdir")     {fGridWorkingDir = name;}
   const char         *GetGridWorkingDir() const                         {return fGridWorkingDir;}
   virtual void        SetGridDataDir(const char *name)                  {fGridDataDir = name;}
   const char         *GetGridDataDir() const                            {return fGridDataDir;}
   void                SetGeneratorLibs(const char *libs)                {fGeneratorLibs = libs;}
   const char         *GetGeneratorLibs() const                          {return fGeneratorLibs;}
   virtual void        SetDataPattern(const char *pattern="*AliESDs.root") {fDataPattern = pattern;}
   const char         *GetDataPattern() const                            {return fDataPattern;}
   virtual void        SetFriendChainName(const char *name="", const char *libnames="");
   const char         *GetFriendChainName() const                        {return fFriendChainName;}
   virtual void        SetDefaultOutputs(Bool_t flag);
   virtual void        SetGridOutputDir(const char *name="output")       {fGridOutputDir = name;}
   const char         *GetGridOutputDir() const                          {return fGridOutputDir;}
   virtual void        SetOutputArchive(const char *list="log_archive.zip:std*@disk=1 root_archive.zip:*.root@disk=2");
   const char         *GetOutputArchive() const                          {return fOutputArchive;}
   virtual void        SetOutputFiles(const char *list);
   const char         *GetOutputFiles() const                            {return fOutputFiles;}
   virtual void        SetOutputToRunNo(Int_t mode=1)                    {fOutputToRunNo = mode;}
   Int_t               GetOutputToRunNoMode() const                      {return fOutputToRunNo;}
   virtual void        SetInputFormat(const char *format="xml-single")   {fInputFormat = format;}
   const char         *GetInputFormat() const                            {return fInputFormat;}   
   virtual void        SetMaxInitFailed(Int_t nfail=5)                   {fMaxInitFailed = nfail;}
   Int_t               GetMaxInitFailed() const                          {return fMaxInitFailed;}
   virtual void        SetTerminateFiles(const char *list)               {fTerminateFiles = list;}
   const char         *GetTerminateFiles() const                         {return fTerminateFiles;}
   virtual void        SetMergeExcludes(const char *list)                {fMergeExcludes = list; fMergeExcludes.ReplaceAll(",", " "); }
   const char         *GetMergeExcludes() const                          {return fMergeExcludes;}
   virtual void        SetMergeViaJDL(Bool_t on=kTRUE)                   {fMergeViaJDL = on ? 1 : 0;}
   Bool_t              IsMergeViaJDL() const                             {return fMergeViaJDL;}
   virtual void        SetMergeDirName(const char *name)                 {fMergeDirName = name;}
   const char         *GetMergeDirName() const                           {return fMergeDirName;}
   virtual void        SetMasterResubmitThreshold(Int_t percentage)      {fMasterResubmitThreshold = percentage;}
   Int_t               GetMasterResubmitThreshold() const                {return fMasterResubmitThreshold;}
   void                SetMCLoop(Bool_t flag=kTRUE)                      {fMCLoop = flag;}
   virtual void        SetNtestFiles(Int_t nfiles)                       {fNtestFiles = nfiles;}
   Int_t               GetNtestFiles() const                             {return fNtestFiles;}
   virtual void        SetNumberOfReplicas(Int_t ncopies)                {fNreplicas = TMath::Min(ncopies,4);}
   Int_t               GetNumberOfReplicas() const                       {return fNreplicas;}
   virtual void        SetJDLName(const char *name="analysis.jdl")       {fJDLName = name;}
   const char         *GetJDLName() const                                {return fJDLName;}
   virtual void        SetProductionMode(Int_t mode=1)                   {fProductionMode = mode;}
   Int_t               GetProductionMode() const                         {return fProductionMode;}
   virtual void        SetRegisterExcludes(const char *list)             {fRegisterExcludes = list; fRegisterExcludes.ReplaceAll(",", " "); }
   const char         *GetRegisterExcludes() const                       {return fRegisterExcludes;}
   virtual void        SetRunPrefix(const char *prefix);
   const char         *GetRunPrefix() const                              {return fRunPrefix;}
   virtual void        SetOutputSingleFolder(const char *folder)         {fOutputSingle = folder; fSplitMode="file"; fSplitMaxInputFileNumber=1;}
   const char         *GetOutputSingleFolder() const                     {return fOutputSingle;}
   virtual void        SetFastReadOption(Bool_t on=kTRUE)                {fFastReadOption = on ? 1 : 0;}
   Bool_t              IsFastReadOption() const                          {return fFastReadOption;}
   virtual void        SetOverwriteMode(Bool_t on=kTRUE)                 {fOverwriteMode = on ? 1 : 0;}
   Bool_t              IsOverwriteMode() const                           {return fOverwriteMode;}
   virtual void        SetDropToShell(Bool_t drop=true)                  {fDropToShell = drop;}
   Bool_t              IsDropToShell() const                             {return fDropToShell;}
   virtual void        SetTreeName(const char *name)                     {fTreeName = name;}
   const char         *GetTreeName() const                               {return fTreeName;}

   TGridJDL           *GetGridJDL() const                                {return fGridJDL;}
   TGridJDL           *GetMergingJDL() const                             {return fMergingJDL;}
   Int_t               GetNMCevents() const                              {return fNMCevents;}
   Int_t               GetNMCjobs() const                                {return fNMCjobs;}
   void                SetNMCevents(Int_t nevents)                       {fNMCevents = nevents;}
   void                SetNMCjobs(Int_t njobs)                           {fNMCjobs = njobs;}
//Utilities
   void                AddModule(AliAnalysisTaskCfg *module);
   void                AddModules(TObjArray *list);
   AliAnalysisManager *CreateAnalysisManager(const char *name, const char *filename="");
   Int_t               GetNmodules() const;
   AliAnalysisTaskCfg *GetModule(const char *name);
   Bool_t              LoadModules();
   Bool_t              LoadFriendLibs() const;
   Bool_t              GenerateTest(const char *name, const char *modname="");
   Bool_t              GenerateTrain(const char *name);
   virtual Bool_t      CreateDataset(const char *pattern);
   Int_t               CopyLocalDataset(const char *griddir, const char *pattern, Int_t nfiles, const char *output="data.txt", const char *archivefile="", const char *outputdir="data");
   virtual Bool_t      CreateJDL();
   virtual void        EnablePackage(const char *package);
   static Bool_t       DirectoryExists(const char *lfn);
   static Bool_t       FileExists(const char *lfn);
   static const char  *GetJobStatus(Int_t jobidstart, Int_t lastid, Int_t &nrunning, Int_t &nwaiting, Int_t &nerror, Int_t &ndone);
   const char         *GetListOfFiles(const char *type);
   Bool_t              CheckMergedFiles(const char *filename, const char *aliendir, Int_t nperchunk, const char *jdl="");
   static AliAnalysisManager *LoadAnalysisManager(const char *fname);
   static Bool_t       MergeInfo(const char *output, const char *collection);
   static Bool_t       MergeOutput(const char *output, const char *basedir, Int_t nmaxmerge, Int_t stage=0);
   virtual Bool_t      MergeOutputs();
   virtual void        Print(Option_t *option="") const;
   static Long64_t     RunMacroAndExtractLibs(const char* macro, const char *args, TString &libs);
   virtual Bool_t      StartAnalysis(Long64_t nentries=123456789, Long64_t firstentry=0);
   static Bool_t       SetupPar(const char *package);
   virtual Bool_t      Submit();
   virtual Bool_t      SubmitMerging();
   static Int_t        SubmitSingleJob(const char *query);
   virtual void        WriteAnalysisFile();
   virtual void        WriteAnalysisMacro(Long64_t nentries=123456789, Long64_t firstentry=0);
   virtual void        WriteMergingMacro();
   virtual void        WriteMergeExecutable();
   virtual void        WriteExecutable();
   virtual Bool_t      WriteJDL(Bool_t copy);
   virtual void        WriteProductionFile(const char *filename) const;
   virtual void        WriteValidationScript(Bool_t merge=kFALSE);

// PROOF mode
   virtual void        SetProofCluster(const char *cluster)              {fProofCluster = cluster;}
   virtual void        SetProofDataSet(const char *dataset)              {fProofDataSet = dataset;}
   virtual const char *GetProofDataSet() const                           {return fProofDataSet;}
   virtual void        SetProofParameter(const char *pname, const char *value);
   const char         *GetProofParameter(const char *pname) const;
   virtual void        SetProofReset(Int_t mode)                         {fProofReset = mode;}
   virtual void        SetNproofWorkers(Int_t nworkers)                  {fNproofWorkers = nworkers;}
   virtual void        SetNproofWorkersPerSlave(Int_t nworkers)          {fNproofWorkersPerSlave = nworkers;}
   virtual void        SetRootVersionForProof(const char *version);
   virtual void        SetAliRootMode(const char *mode)                  {fAliRootMode = mode;}
   virtual void        SetProofProcessOpt(const char *proofOpt="")       {fProofProcessOpt = proofOpt;}
   virtual TString     GetProofProcessOpt()                              {return fProofProcessOpt;}
   // .txt file containing the list of files to be chained in test mode
   virtual void        SetFileForTestMode(const char *filename)          {fFileForTestMode = filename;}
   const char         *GetFileForTestMode() const                        {return fFileForTestMode;}
   virtual TChain     *GetChainForTestMode(const char *treeName) const;
   virtual const TString& GetGridJobIDs() const { return fGridJobIDs; }
   virtual const TString& GetGridStages() const { return fGridStages; }
protected:
   void                CdWork();
   Bool_t              CheckInputData();
   void                CheckDataType(const char *lfn, Bool_t &is_collection, Bool_t &is_xml, Bool_t &use_tags);
   virtual Bool_t      Connect();
   virtual void        SetDefaults();  
   Bool_t              SubmitNext();

   Bool_t              IsCollection(const char *lfn) const;
   Bool_t              IsMCLoop() const {return fMCLoop;}
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
   Int_t            fNMCevents;       // Number of MC events in MC loop mode
   Int_t            fNMCjobs;         // Number of MC jobs in MC loop mode
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
   TString          fGeneratorLibs;   // Extra libraries needed by the generator
   TString          fSplitMode;       // Job split mode
   TString          fAPIVersion;      // API version
   TString          fROOTVersion;     // ROOT version
   TString          fAliROOTVersion;  // AliROOT version
   TString          fAliPhysicsVersion; // AliPhysics version
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
   TString          fAliRootMode;     // AliRoot mode among the list supported by the proof cluster
   TString          fProofProcessOpt; // Option passed to proof process
   TString          fMergeDirName;    // Name of the directory that should be added to the output directory
   TObjArray       *fInputFiles;      // List of input files to be processed by the job
   TObjArray       *fPackages;        // List of packages to be used
   TObjArray       *fModules;         // List of AliAnalysisTaskCfg modules
   TMap             fProofParam;      // Key-value pairs for proof mode
   Bool_t           fDropToShell;     // If true, execute aliensh on start
   Bool_t           fMCLoop;          // MC loop flag
   TString          fGridJobIDs;      // List of last committed jobs
   TString          fGridStages;      // List of last committed jobs
   TString          fFriendLibs;      // List of libs (separated by blacs) needed for friends processing
   TString          fTreeName;        // Name of the tree to be analyzed

   ClassDef(AliAnalysisAlien, 27)   // Class providing some AliEn utilities
};
#endif
