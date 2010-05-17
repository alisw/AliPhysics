#ifndef ALIANALYSISGRID_H
#define ALIANALYSISGRID_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Mihaela Gheata, 01/09/2008

//==============================================================================
//   AliAnalysisGrid - Base grid utility class. Provides interface for creating
// a personalized JDL, finding and creating a dataset.
//==============================================================================

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif

class AliAnalysisGrid : public TNamed {

public:

enum EPluginRunMode {
   kFull     = 0,
   kTest     = BIT(14),
   kOffline  = BIT(15),
   kSubmit   = BIT(16),
   kMerge    = BIT(17),
   kUseTags  = BIT(18),
   kUseESD   = BIT(19),
   kUseAOD   = BIT(20),
   kUseMC    = BIT(21),
   kUsePars  = BIT(22),
   kDefaultOutputs = BIT(23)
};   
enum EPluginBits {
   kBitMask32  = 0xffffffff,
   kUseCopy    = BIT(0),
   kCheckCopy  = BIT(1)
};

   AliAnalysisGrid() : TNamed(), fSpecialBits(0) {}
   AliAnalysisGrid(const char *name) : TNamed(name,""), fSpecialBits(0) {}
   virtual ~AliAnalysisGrid() {}
   AliAnalysisGrid(const AliAnalysisGrid& other); 
   AliAnalysisGrid& operator=(const AliAnalysisGrid& other);
// Getters
   virtual EPluginRunMode GetRunMode() const;
// Setters   
   virtual void        AddIncludePath(const char *path)                  = 0;
   virtual void        AddRunNumber(Int_t run)                           = 0;
   virtual void        AddRunNumber(const char *run)                     = 0;
   virtual void        AddDataFile(const char *lfn)                      = 0;
   virtual Bool_t      IsSingleOutput() const                            = 0;
   virtual void        SetExecutable(const char *name="analysis.sh")     = 0;
   virtual void        SetArguments(const char *name="")                 = 0;
   virtual void        SetAnalysisMacro(const char *name="myAnalysis.C") = 0;
   virtual void        SetAnalysisSource(const char *name="myAnalysisClass.cxx") = 0;
   virtual void        SetAdditionalLibs(const char *list)               = 0;
   virtual void        SetPrice(Int_t price=1)                           = 0;
   virtual void        SetJobTag(const char *tag="")                     = 0;
   virtual void        SetSplitMode(const char *type="se")               = 0;
   virtual void        SetSplitMaxInputFileNumber(Int_t nfiles=100)      = 0;
   virtual void        SetAPIVersion(const char *version)                = 0;
   virtual void        SetROOTVersion(const char *version)               = 0;
   virtual void        SetAliROOTVersion(const char *version)            = 0;
   virtual void        SetUser(const char *user)                         = 0;
   virtual void        SetTTL(Int_t ttl=30000)                           = 0;
   virtual void        SetGridWorkingDir(const char *name="workdir")     = 0;
   virtual void        SetGridDataDir(const char *name)                  = 0;
   virtual void        SetDataPattern(const char *pattern)               = 0;
   virtual void        SetDefaultOutputs(Bool_t flag=kTRUE)              = 0;
   virtual void        SetGridOutputDir(const char *name="output")       = 0;
   virtual void        SetOutputArchive(const char *list="log_archive.zip:stdout,stderr root_archive.zip:*.root") = 0;
   virtual void        SetOutputFiles(const char *list)                  = 0;
   virtual void        SetInputFormat(const char *format="xml-single")   = 0;
   virtual void        SetMaxInitFailed(Int_t nfail=5)                   = 0;
   virtual void        SetMergeExcludes(const char *list)                = 0;
   virtual void        SetMergeViaJDL(Bool_t on=kTRUE)                   = 0;
   virtual void        SetMasterResubmitThreshold(Int_t percentage)      = 0;
   virtual void        SetNtestFiles(Int_t nfiles)                       = 0;
   virtual void        SetJDLName(const char *name="analysis.jdl")       = 0;
   virtual void        SetPreferedSE(const char *se)                     = 0;
   virtual void        SetProductionMode(Int_t mode=1)                   = 0;
   virtual void        SetRunPrefix(const char *prefix)                  = 0;
   virtual void        SetOutputSingleFolder(const char *folder)         = 0;
   virtual void        SetFastReadOption(Bool_t on=kTRUE)                = 0;
   virtual void        SetOverwriteMode(Bool_t on=kTRUE)                 = 0;
   
 // Set run mode.  Can be "full", "test", "offline", "submit" or "merge"
   virtual void        SetRunMode(const char *mode="full");
//Utilities
   static  Bool_t      CreateToken(const char *username=0);
   virtual Bool_t      CreateDataset(const char *pattern)                = 0;
   virtual Bool_t      CreateJDL()                                       = 0;
   virtual void        EnablePackage(const char *package)                = 0;
   virtual Bool_t      MergeOutputs()                                    = 0;
   virtual Bool_t      StartAnalysis(Long64_t nentries=123456789, Long64_t firstentry=0) = 0;
   virtual void        WriteAnalysisFile()                               = 0;
   virtual void        WriteAnalysisMacro()                              = 0;
   virtual void        WriteExecutable()                                 = 0;
   virtual void        WriteValidationScript(Bool_t merge=kFALSE)        = 0;

// Flags
   Bool_t              IsUseCopy() const {return TestSpecialBit(kUseCopy);}
   void                SetUseCopy(Bool_t flag=kTRUE) {SetSpecialBit(kUseCopy,flag);}
   Bool_t              IsCheckCopy() const {return TestSpecialBit(kCheckCopy);}
   void                SetCheckCopy(Bool_t flag=kTRUE) {SetSpecialBit(kCheckCopy,flag);}

protected:
// Methods
   virtual Bool_t      Connect()                                         = 0;
   virtual void        SetDefaults()                                     = 0;
   void     SetSpecialBit(UInt_t f) { fSpecialBits |= f & kBitMask32; }
   void     ResetSpecialBit(UInt_t f) { fSpecialBits &= ~(f & kBitMask32); }
   void     SetSpecialBit(UInt_t f, Bool_t set) {(set)?SetSpecialBit(f):ResetSpecialBit(f);}
   Bool_t   TestSpecialBit(UInt_t f) const { return (Bool_t) ((fSpecialBits & f) != 0); }
   Int_t    TestSpecialBits(UInt_t f) const { return (Int_t) (fSpecialBits & f); }
   void     InvertSpecialBit(UInt_t f) { fSpecialBits ^= f & kBitMask32; }

protected:
   UInt_t              fSpecialBits; // special bits
  

   ClassDef(AliAnalysisGrid, 2)   // Base class for GRID utilities
};
#endif
