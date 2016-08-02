#ifndef ALIMUONACCEFFSUBMITTER_H
#define ALIMUONACCEFFSUBMITTER_H

//
// AliMuonAccEffSubmitter : a class to help submit Acc x Eff simulations
// anchored to real runs for J/psi, upsilon, single muons, etc...
//
// author: Laurent Aphecetche (Subatech)
//

#include "AliMuonGridSubmitter.h"

class AliMuonAccEffSubmitter : public AliMuonGridSubmitter
{
public:
  AliMuonAccEffSubmitter(const char* generator="GenParamCustom",
                         Bool_t localOnly=kFALSE,
                         const char* generatorVersion="8125");

  virtual Bool_t Generate(const char* jdlname) const;
  virtual Bool_t Run(const char* mode);

  virtual ~AliMuonAccEffSubmitter();

  void MakeNofEventsPropToTriggerCount(const char* trigger="CMUL7-B-NOPF-MUON", Float_t ratio=1.0);
  
  void MakeNofEventsFixed(Int_t nevents);
  
  void UseOCDBSnapshots(Bool_t flag);
  
  void UseExternalConfig(const char* externalConfigFullFilePath);
  
  void UseAODMerging(Bool_t flag);

  void SetupPythia6 ( const char *version );
  void SetupPythia8 ( const char *version, const char* configStrings = "" );
  
  Bool_t Merge(Int_t stage, Bool_t dryRun=kTRUE);

  Int_t Submit(Bool_t dryRun=kTRUE);
  
  Int_t LocalTest();
  
  TString RunJDLName() const { return "JDL"; }

  TString MergeJDLName(Bool_t final) const { return (final ? "AOD_merge_final.jdl" : "AOD_merge.jdl"); }

  virtual void Print(Option_t* opt="") const;

  void SetSplitMaxInputFileNumber(Int_t n) { fSplitMaxInputFileNumber = n; }
  
  Int_t GetSplitMaxInputFileNumber() const { return fSplitMaxInputFileNumber; }
  
  void SetCompactMode(Int_t mode);
  
  /// Set the output files to be kept
  void SetCustomOutFiles(const char* logFiles, const char* rootFiles) { fLogOutToKeep = logFiles; fRootOutToKeep = rootFiles; }
  
  Bool_t MakeOCDBSnapshots();
  
  void SetOCDBPath(const char* ocdbPath);

  void SetOCDBSnapshotDir(const char* dir);

  Bool_t SetGenerator(const char* generator);

  Int_t MaxEventsPerChunk() const { return fMaxEventsPerChunk; }
  void SetMaxEventsPerChunk(Int_t n) { fMaxEventsPerChunk = n; SetVar("VAR_EVENTS_PER_JOB", Form("%i",n)); }

  TString OCDBPath() const { return GetMapValue("OCDBPath"); }
  
  TString SnapshotDir() const { return GetMapValue("OCDBsnapshot"); }

  Int_t SplitRunList(const char* inputList, int maxJobs=1500);
  
private:
  
  Bool_t GenerateRunJDL(const char* name) const;
  
  Bool_t GenerateMergeJDL(const char* name) const;
  
  TString ReferenceTrigger() const { return GetMapValue("ReferenceTrigger"); }
  
  void UpdateLocalFileList(Bool_t clearSnapshots=kFALSE);

private:
  AliMuonAccEffSubmitter(const AliMuonAccEffSubmitter& rhs);
  AliMuonAccEffSubmitter& operator=(const AliMuonAccEffSubmitter& rhs);
  
private:
  Float_t fRatio; // ratio simulated events vs real events
  Int_t fFixedNofEvents; // fixed number of events to be used per run
  Int_t fMaxEventsPerChunk; // max events to generate per subjob
  TString fOCDBPath; // OCDB path
  Int_t fSplitMaxInputFileNumber; // used for merging jdl
  TString fLogOutToKeep; // specify the log files to be kept
  TString fRootOutToKeep; // specify the output files to be kept
  TString fExternalConfig; // path to an (optional) external config file
  Bool_t fUseOCDBSnapshots; // whether to use OCDB snapshots or not
  TString fSnapshotDir; // directory for OCDB snapshots
  Bool_t fUseAODMerging; // whether or not to perform (aod) merging
  
  ClassDef(AliMuonAccEffSubmitter,3) // Helper class to submit AccxEff single particle simulations
};

#endif

