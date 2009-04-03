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
   virtual void        AddDataFile(const char *lfn);
   virtual void        SetExecutable(const char *name="analysis.sh")     {fExecutable = name;}
   virtual void        SetArguments(const char *name="")                 {fArguments = name;}
   virtual void        SetAnalysisMacro(const char *name="myAnalysis.C") {fAnalysisMacro = name;}
   virtual void        SetAnalysisSource(const char *name="myAnalysisClass.cxx") {fAnalysisSource = name;}
   virtual void        SetAdditionalLibs(const char *list)               {fAdditionalLibs = list;}
   virtual void        SetPrice(Int_t price=1)                           {fPrice = price;}
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
   virtual void        SetDefaultOutputs(Bool_t flag);
   virtual void        SetGridOutputDir(const char *name="output")       {fGridOutputDir = name;}
   virtual void        SetOutputArchive(const char *list="log_archive.zip:stdout,stderr root_archive.zip:*.root") {fOutputArchive = list;}
   virtual void        SetOutputFiles(const char *list)                  {fOutputFiles = list;}
   virtual void        SetInputFormat(const char *format="xml-single")   {fInputFormat = format;}
   virtual void        SetMaxInitFailed(Int_t nfail=5)                   {fMaxInitFailed = nfail;}
   virtual void        SetMergeExcludes(const char *list)                {fMergeExcludes = list;};
   virtual void        SetMasterResubmitThreshold(Int_t percentage)      {fMasterResubmitThreshold = percentage;}
   virtual void        SetNtestFiles(Int_t nfiles)                       {fNtestFiles = nfiles;}
   virtual void        SetJDLName(const char *name="analysis.jdl")       {fJDLName = name;}
   virtual void        SetPreferedSE(const char *se)                     {fCloseSE = se;}

   TGridJDL           *GetGridJDL() {return fGridJDL;}
//Utilities
   virtual Bool_t      CreateDataset(const char *pattern);
   virtual Bool_t      CreateJDL();
   virtual void        EnablePackage(const char *package);
   virtual Bool_t      MergeOutputs();
   virtual void        StartAnalysis(Long64_t nentries=123456789, Long64_t firstentry=0);
   virtual void        WriteAnalysisFile();
   virtual void        WriteAnalysisMacro();
   virtual void        WriteExecutable();
   virtual void        WriteValidationScript();

protected:
   void                CdWork();
   Bool_t              CheckInputData();
   void                CheckDataType(const char *lfn, Bool_t &is_collection, Bool_t &is_xml, Bool_t &use_tags);
   virtual Bool_t      Connect();
   virtual void        SetDefaults();  

   Bool_t              FileExists(const char *lfn) const;
   Bool_t              IsCollection(const char *lfn) const;
   Bool_t              IsUsingTags() const {return TObject::TestBit(AliAnalysisGrid::kUseTags);}

private:
   TGridJDL        *fGridJDL;         //! JDL maker
   Int_t            fPrice;           // Grid price for the job;
   Int_t            fTTL;             // Time to live.
   Int_t            fSplitMaxInputFileNumber; // Maximum number of files to be processed per subjob
   Int_t            fMaxInitFailed;   // Maximum initial consecutive subjobs accepted to fail
   Int_t            fMasterResubmitThreshold; // Failed jobs will be resubmitted until this DONE ratio
   Int_t            fNtestFiles;      // Number of files used in the testing case
   TString          fRunNumbers;      // List of runs to be processed
   TString          fExecutable;      // Executable script for AliEn job
   TString          fArguments;       // Arguments for the executable script
   TString          fAnalysisMacro;   // Root macro steering the analysis
   TString          fAnalysisSource;  // User analysis implementation (.cxx) file(s)
   TString          fAdditionalLibs;  // List (separated by blacs) of additional libraries needed for the analysis
   TString          fSplitMode;       // Job split mode
   TString          fAPIVersion;      // API version
   TString          fROOTVersion;     // ROOT version
   TString          fAliROOTVersion;  // AliROOT version
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
   TObjArray       *fInputFiles;      // List of input files to be processed by the job
   TObjArray       *fPackages;        // List of packages to be used
   
   ClassDef(AliAnalysisAlien, 2)   // Class providing some AliEn utilities
};
#endif
