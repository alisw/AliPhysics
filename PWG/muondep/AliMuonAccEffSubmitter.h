#ifndef ALIMUONACCEFFSUBMITTER_H
#define ALIMUONACCEFFSUBMITTER_H

//
// AliMuonAccEffSubmitter : a class to help submit Acc x Eff simulations
// anchored to real runs for J/psi, upsilon, single muons, etc...
//
// author: Laurent Aphecetche (Subatech)
//

#include "TObject.h"
#include "TString.h"
#include "Riostream.h"

class AliAnalysisTriggerScalers;
class TMap;

class AliMuonAccEffSubmitter : public TObject
{
public:
  AliMuonAccEffSubmitter(const char* generator="GenParamCustom");

  virtual ~AliMuonAccEffSubmitter();

  Bool_t SetRemoteDir(const char* dir, Bool_t create = kTRUE);
  void SetLocalDir(const char* localdir) { fLocalDir = localdir; }
  Bool_t SetMergedDir(const char* dir, Bool_t create = kTRUE);

  void UseOCDBSnapshots(Bool_t flag);
  
  void UseExternalConfig(const char* externalConfigFullFilePath) { fExternalConfig = externalConfigFullFilePath; }
  
  Bool_t CheckLocal() const;
  Bool_t CheckRemote() const;
  
  void CleanLocal(Bool_t cleanSnapshots=kFALSE) const;
  void CleanRemote() const;

  TString MergedDir() const { return fMergedDir; }  
  TString RemoteDir() const { return fRemoteDir; }
  TString LocalDir() const { return fLocalDir; }
  
  TString FilePath(const char* what) const;

  Int_t MaxEventsPerChunk() const { return fMaxEventsPerChunk; }
  void SetMaxEventsPerChunk(Int_t n) { fMaxEventsPerChunk = n; }

  UInt_t NofRuns() const;
  
  void SetRatio(Float_t ratio) { fRatio = ratio; }
  void SetOCDBPath(const char* ocdbPath);
  
  void MakeNofEventsPropToTriggerCount(const char* trigger="CMUL8-S-NOPF-MUON", Float_t ratio=1.0) { fReferenceTrigger = trigger; fRatio = ratio; }
  void MakeNofEventsFixed(Int_t nevents) { fFixedNofEvents = nevents; fReferenceTrigger=""; fRatio=0.0; }
  
  void SetRunList(const char* runlist);
  void SetRunList(int runNumber);
  
  TString ReferenceTrigger() const { return fReferenceTrigger; }
  
  Bool_t Run(const char* mode);
  
  Bool_t Merge(Int_t stage, Bool_t dryRun=kTRUE);

  Int_t Submit(Bool_t dryRun=kTRUE);
  
  TString RunJDLName() const { return "run.jdl"; }

  TString MergeJDLName(Bool_t final) const { return (final ? "AOD_merge_final.jdl" : "AOD_merge.jdl"); }

  virtual void Print(Option_t* opt="") const;

  Bool_t CopyLocalFilesToRemote();

  Bool_t CopyTemplateFilesToLocal();

  Bool_t GenerateRunJDL(const char* name);

  Bool_t GenerateMergeJDL(const char* name);

  void SetPackages(const char* aliroot, const char* root, const char* geant3,
                   const char* api="VO_ALICE@APISCONFIG::V1.1x");
  
  void SetSplitMaxInputFileNumber(Int_t n) { fSplitMaxInputFileNumber = n; }
  
  Int_t GetSplitMaxInputFileNumber() const { return fSplitMaxInputFileNumber; }
  
  Int_t CompactMode() const { return fCompactMode; }
  
  void SetCompactMode(Int_t mode) { fCompactMode=mode; }
  
  Bool_t ShouldOverwriteFiles() const { return fShouldOverwriteFiles; }

  void ShouldOverwriteFiles(Bool_t flag) { fShouldOverwriteFiles = flag; }

  Bool_t SetVar(const char* varname, const char* value);
  
  Bool_t MakeOCDBSnapshots();
  
  void SetOCDBSnapshotDir(const char* dir);

  Bool_t SetGenerator(const char* generator);
  
  TObjArray* GetVariables(const char* file) const;
  
  Bool_t IsValid() const { return fIsValid; }
  
private:

  TString SnapshotDir() const { return fSnapshotDir; }
  
  TString GetRemoteDir(const char* dir, Bool_t create);

  std::ostream* CreateJDLFile(const char* name) const;

  Bool_t CheckRemoteDir() const;

  Bool_t CopyFile(const char* localFile);
  
  Bool_t GetLastStage(const char* remoteDir) const;

  Bool_t RemoteDirectoryExists(const char *dirname) const;
  Bool_t RemoteFileExists(const char *lfn) const;

  void Output(std::ostream& out, const char* key, const char* v1,
              const char* v2="", const char* v3="", const char* v4="", const char* v5="",
              const char* v6="", const char* v7="", const char* v8="", const char* v9="") const;
  
  void Output(std::ostream& out, const char* key, const TObjArray& values) const;

  Bool_t ReplaceVars(const char* file) const;

  TObjArray* TemplateFileList() const;

  TObjArray* LocalFileList() const;
  
  Bool_t HasVars(const char* localFile) const;

  void UpdateLocalFileList(Bool_t clearSnapshot=kFALSE);
  
  Bool_t CheckCompilation(const char* file) const;
  
private:
  AliMuonAccEffSubmitter(const AliMuonAccEffSubmitter& rhs);
  AliMuonAccEffSubmitter& operator=(const AliMuonAccEffSubmitter& rhs);
  
private:
  AliAnalysisTriggerScalers* fScalers; // helper class used to handle the runlist and the scalers
  TString fRemoteDir; // remote directory to used in alien
  TString fReferenceTrigger; // reference trigger (if any) to be used to get the number of events to be used per run
  Float_t fRatio; // ratio simulated events vs real events
  Int_t fFixedNofEvents; // fixed number of events to be used per run
  Int_t fMaxEventsPerChunk; // max events to generate per subjob
  TString fLocalDir; // local directory
  TString fOCDBPath; // OCDB path
  TString fTemplateDir; // template directory
  TString fPackageAliroot; // which aliroot package to use
  TString fPackageGeant3; // which geant3 package to use
  TString fPackageRoot; // which root package to use (for valid root,geant3,aliroot combinations see http://alimonitor.cern.ch/packages/)
  TString fPackageApi; // which API package to use
  TString fMergedDir; // merge directory
  Int_t fSplitMaxInputFileNumber; // used for merging jdl
  Int_t fCompactMode; // controls which outputs are kept (0=everything, 1=only aods)
  Bool_t fShouldOverwriteFiles; // whether any copy (of template to local) is allowed to overwrite existing files
  TMap* fVars; // map of the variables we can replace in template files
  TString fExternalConfig; // path to an (optional) external config file
  Bool_t fUseOCDBSnapshots; // whether to use OCDB snapshots or not
  Bool_t fIsValid; // whether this object is valid (i.e. properly configured)
  mutable TObjArray* fTemplateFileList; // list of template files
  mutable TObjArray* fLocalFileList; // list of local files
  TString fSnapshotDir; // directory for OCDB snapshots
  
  ClassDef(AliMuonAccEffSubmitter,1) // Helper class to submit AccxEff single particle simulations
};

#endif

