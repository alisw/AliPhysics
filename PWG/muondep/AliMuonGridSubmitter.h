#ifndef ALIMUONGRIDSUBMITTER_H
#define ALIMUONGRIDSUBMITTER_H

//
// AliMuonGridSubmitter : a base class to help submit some
// mini-productions for muon studies (e.g. Acc x Eff, qa merging, etc...)
//
// author: Laurent Aphecetche (Subatech)
//

#include "TObject.h"
#include "TString.h"
#include "Riostream.h"
#include <vector>

class TMap;

class AliMuonGridSubmitter : public TObject
{
public:
  
  enum EJobType
  {
    kAccEff=0,
    kQAMerge=1
  };
  
  AliMuonGridSubmitter(AliMuonGridSubmitter::EJobType jobType, Bool_t localOnly=kFALSE);
  virtual ~AliMuonGridSubmitter();

  virtual Bool_t Generate(const char* jdlname) const = 0;
  virtual Bool_t Run(const char* mode) = 0;
  
  TString JobTypeName(AliMuonGridSubmitter::EJobType jobType) const;
    
  Bool_t SetLocalDir(const char* localdir) { return SetLocalDirectory("Local",localdir); }
  Bool_t SetMergedDir(const char* dir) { return SetRemoteDirectory("Merged",dir); }
  Bool_t SetRemoteDir(const char* dir) { return SetRemoteDirectory("Remote",dir); }
  Bool_t SetTemplateDir(const char* templatedir) { return SetLocalDirectory("Template",templatedir); }
  
  Bool_t CheckLocal() const;
  Bool_t CheckRemote() const;
  
  void CleanLocal(const char* excludeList="") const;
  
  void Invalidate() { fIsValid = kFALSE; }

  void Validate() { fIsValid = kTRUE; }

  TString MergedDir() const { return GetMapValue("Merged"); }
  TString RemoteDir() const { return GetMapValue("Remote"); }
  TString LocalDir() const { return GetMapValue("Local"); }
  TString TemplateDir() const { return GetMapValue("Template"); }
  
  //  TString FilePath(const char* what) const; // Not implemented
  
  UInt_t NofRuns() const;
  
  const std::vector<Int_t>& RunList() const;
  
  void SetRunList(const char* runlist);
  void SetRunList(int runNumber);
  
  virtual void Print(Option_t* opt="") const;

  Bool_t CopyLocalFilesToRemote();

  Bool_t CopyTemplateFilesToLocal();


  void SetAliPhysicsVersion(const char* aliphysics);

  void SetAliRootVersion(const char* aliroot);
  
  void SetPackages(const char* aliroot, const char* root="", const char* geant3="",
                   const char* api="VO_ALICE@APISCONFIG::V1.1x");

  void SetGeneratorPackage(const char* generator);
  
  Bool_t ShouldOverwriteFiles() const { return fShouldOverwriteFiles; }

  void ShouldOverwriteFiles(Bool_t flag) { fShouldOverwriteFiles = flag; }

  Bool_t IsValid() const { return fIsValid; }
  
  void OutputToJDL(std::ostream& out, const char* key, const char* v1,
              const char* v2="", const char* v3="", const char* v4="", const char* v5="",
              const char* v6="", const char* v7="", const char* v8="", const char* v9="") const;
  
  void OutputToJDL(std::ostream& out, const char* key, const TObjArray& values) const;
  
  TString GetRemoteDir(const char* dir, Bool_t create=kTRUE);
  
  Bool_t RemoteDirectoryExists(const char *dirname) const;
  Bool_t RemoteFileExists(const char *lfn);
  
  //  Bool_t CopyLocalFilesToRemote(const TObjArray& localFiles); // Not implemented

  Bool_t CopyFile(const char* localFile);
  
  Int_t GetLastStage(const char* remoteDir);
  
  std::ostream* CreateJDLFile(const char* name) const;

  Bool_t CheckRemoteDir() const;

  Bool_t ReplaceVars(const char* file) const;
  
  Bool_t HasVars(const char* localFile) const;

  Bool_t SetVar(const char* varname, const char* value);
  
  TObjArray* GetVariables(const char* file) const;
  
  Bool_t CheckCompilation(const char* file) const;

  TObjArray* LocalFileList() const;
  
  TObjArray* TemplateFileList() const;
  
  void AddToTemplateFileList(const char* filename);

  void AddToLocalFileList(const char* filename);

  void AddIncludePath(const char* pathList) const;

  TString GetVar(const char* key) const;

protected:
  
  TObjArray* OrderKeys(const TMap& map) const;

  std::ostream* CreateJDLFile(const char* name);

  TString GetMapValue(const char* key) const;
  
  TMap* InternalMap() const;
  TMap* Vars() const;
  
  void SetMapKeyValue(const char* key, const char* value);
  
  Bool_t SetLocalDirectory(const char* type, const char* path);

  Bool_t SetRemoteDirectory(const char* type, const char* path);

  void UpdateLocalFileList();
  
private:
  AliMuonGridSubmitter(const AliMuonGridSubmitter& rhs);
  AliMuonGridSubmitter& operator=(const AliMuonGridSubmitter& rhs);
  
private:
  mutable TMap* fInternalMap; // map of directory paths and packages versions
  mutable TMap* fVars; // map of the variables we can replace in template files
  Bool_t fIsValid; // whether this object is valid (i.e. properly configured)
  Bool_t fShouldOverwriteFiles; // whether or not to overwrite the local files each time we run
  mutable TObjArray* fTemplateFileList; // list of template files
  mutable TObjArray* fLocalFileList; // list of local files
  std::vector<Int_t> fRunList; // run list to process

  ClassDef(AliMuonGridSubmitter,0) // Helper class to submit some muon jobs
};

#endif

