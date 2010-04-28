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

class TGridJDL;

class AliAnalysisAlien : public AliAnalysisGrid {

public:

   AliAnalysisAlien();
   AliAnalysisAlien(const char *name);
   virtual ~AliAnalysisAlien();
   AliAnalysisAlien(const AliAnalysisAlien& other); 
   AliAnalysisAlien& operator=(const AliAnalysisAlien& other);
// Setters   
   virtual void        AddIncludePath(const char *path);
   virtual void        AddRunNumber(Int_t run);
   virtual void        AddRunNumber(const char *run);
   virtual void        AddDataFile(const char *lfn);
   virtual void        AddExternalPackage(const char *name);
   virtual void        SetExecutable(const char *name="analysis.sh")     {fExecutable = name;}
   virtual void        SetExecutableCommand(const char *command="root -b -q") {fExecutableCommand = command;}
   virtual void        SetArguments(const char *name="")                 {fArguments = name;}
   virtual void        SetExecutableArgs(const char *name="")            {fExecutableArgs = name;}
   virtual void        SetAnalysisMacro(const char *name="myAnalysis.C") {fAnalysisMacro = name;}
   virtual void        SetAnalysisSource(const char *name="myAnalysisClass.cxx") {fAnalysisSource = name;}
   virtual void        SetAdditionalLibs(const char *list)               {fAdditionalLibs = list;}
   virtual void        SetAdditionalRootLibs(const char *list)           {fAdditionalRootLibs = list;}
   virtual void        SetPrice(Int_t price=1)                           {fPrice = price;}
   virtual void        SetRunRange(Int_t min, Int_t max)                 {fRunRange[0] = min; fRunRange[1] = max;}
   virtual void        SetJobTag(const char *tag="")                     {fJobTag = tag;}
   virtual void        SetNrunsPerMaster(Int_t nruns=1)                  {fNrunsPerMaster = nruns;}
   virtual void        SetMaxMergeFiles(Int_t nfiles)                    {fMaxMergeFiles = nfiles;}
   virtual void        SetSplitMode(const char *type="se")               {fSplitMode = type;}
   virtual void        SetSplitMaxInputFileNumber(Int_t nfiles=100)      {fSplitMaxInputFileNumber = nfiles;}
   virtual void        SetAPIVersion(const char *version="V2.4") {fAPIVersion = version;}
   virtual void        SetROOTVersion(const char *version="v5-21-01-alice") {fROOTVersion = version;}
   virtual void        SetAliROOTVersion(const char *version="v4-14-Rev-02") {fAliROOTVersion=version;}
   virtual void        SetUser(const char *user)                         {fUser = user;}
   virtual void        SetTTL(Int_t ttl=30000)                           {fTTL = ttl;}
   virtual void        SetGridWorkingDir(const char *name="workdir")     {fGridWorkingDir = name;}
   virtual void        SetGridDataDir(const char *name)                  {fGridDataDir = name;}
   virtual void        SetDataPattern(const char *pattern="*AliESDs.root") {fDataPattern = pattern;}
   virtual void        SetFriendChainName(const char *name="")           {fFriendChainName = name;}
   virtual void        SetDefaultOutputs(Bool_t flag);
   virtual void        SetGridOutputDir(const char *name="output")       {fGridOutputDir = name;}
   virtual void        SetOutputArchive(const char *list="log_archive.zip:stdout,stderr root_archive.zip:*.root") {fOutputArchive = list;}
   virtual void        SetOutputFiles(const char *list)                  {fOutputFiles = list;}
   virtual void        SetOutputToRunNo(Int_t mode=1)                    {fOutputToRunNo = mode;}
   virtual void        SetInputFormat(const char *format="xml-single")   {fInputFormat = format;}
   virtual void        SetMaxInitFailed(Int_t nfail=5)                   {fMaxInitFailed = nfail;}
   virtual void        SetMergeExcludes(const char *list)                {fMergeExcludes = list;};
   virtual void        SetMergeViaJDL(Bool_t on=kTRUE)                   {fMergeViaJDL = on ? 1 : 0;}
   virtual void        SetMasterResubmitThreshold(Int_t percentage)      {fMasterResubmitThreshold = percentage;}
   virtual void        SetNtestFiles(Int_t nfiles)                       {fNtestFiles = nfiles;}
   virtual void        SetJDLName(const char *name="analysis.jdl")       {fJDLName = name;}
   virtual void        SetPreferedSE(const char *se)                     {fCloseSE = se;}
   virtual void        SetProductionMode(Int_t mode=1)                   {fProductionMode = mode;}
   virtual void        SetRunPrefix(const char *prefix)                  {fRunPrefix = prefix;}
   virtual void        SetOutputSingleFolder(const char *folder)         {fOutputSingle = folder; fSplitMode="file"; fSplitMaxInputFileNumber=1;}
   virtual void        SetFastReadOption(Bool_t on=kTRUE)                {fFastReadOption = on ? 1 : 0;}
   virtual void        SetOverwriteMode(Bool_t on=kTRUE)                 {fOverwriteMode = on ? 1 : 0;}

   TGridJDL           *GetGridJDL() const {return fGridJDL;}
   TGridJDL           *GetMergingJDL() const {return fMergingJDL;}
   const char         *GetGridOutputDir() const                          {return fGridOutputDir;}
//Utilities
   virtual Bool_t      CreateDataset(const char *pattern);
   virtual Bool_t      CreateJDL();
   virtual void        EnablePackage(const char *package);
   static Bool_t       DirectoryExists(const char *lfn);
   static Bool_t       FileExists(const char *lfn);
   static const char  *GetJobStatus(Int_t jobidstart, Int_t lastid, Int_t &nrunning, Int_t &nwaiting, Int_t &nerror, Int_t &ndone);
   static Bool_t       MergeOutput(const char *output, const char *basedir, Int_t nmaxmerge);
   virtual Bool_t      MergeOutputs();
   virtual void        Print(Option_t *option="") const;
   virtual Bool_t      StartAnalysis(Long64_t nentries=123456789, Long64_t firstentry=0);
   static Bool_t       SetupPar(const char *package);
   virtual Bool_t      Submit();
   virtual Bool_t      SubmitMerging();
   virtual void        WriteAnalysisFile();
   virtual void        WriteAnalysisMacro();
   virtual void        WriteMergingMacro();
   virtual void        WriteMergeExecutable();
   virtual void        WriteExecutable();
   virtual Bool_t      WriteJDL(Bool_t copy);
   virtual void        WriteProductionFile(const char *filename) const;
   virtual void        WriteValidationScript(Bool_t merge=kFALSE);

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
   Int_t            fNsubmitted;      // Number of jobs submitted
   Int_t            fProductionMode;  // Production mode (0-off, 1-on)
   Int_t            fOutputToRunNo;   // Use run number as output directory
   Int_t            fMergeViaJDL;     // Enable merging via automatic JDL
   Int_t            fFastReadOption;  // Use xrootd tweaks to reduce timeouts in file access
   Int_t            fOverwriteMode;   // Overwrite existing files if any
   TString          fRunNumbers;      // List of runs to be processed
   TString          fExecutable;      // Executable script for AliEn job
   TString          fExecutableCommand;  // Command(s) to be executed in the executable script
   TString          fArguments;       // Arguments for the executable script
   TString          fExecutableArgs;  // arguments added to the executable script after the analysis macro
   TString          fAnalysisMacro;   // Root macro steering the analysis
   TString          fAnalysisSource;  // User analysis implementation (.cxx) file(s)
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
   TString          fMergeExcludes;   // List of output files excluded from merging
   TString          fIncludePath;     // Include path
   TString          fCloseSE;         // Preffered storage element. Taken from alien_CLOSE_SE environment.
   TString          fFriendChainName; // File name to construct friend chain (for AOD)
   TString          fJobTag;          // Job tag
   TString          fOutputSingle;    // Directory name for the output when split is per file
   TString          fRunPrefix;       // Run prefix to be applied to run numbers
   TObjArray       *fInputFiles;      // List of input files to be processed by the job
   TObjArray       *fPackages;        // List of packages to be used
   
   ClassDef(AliAnalysisAlien, 12)   // Class providing some AliEn utilities
};
#endif
