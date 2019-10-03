#ifndef ALIOFFLINETRIGGER_H
#define ALIOFFLINETRIGGER_H

/// \class AliOfflineTrigger
/// \brief This class provides fucntionality for OFFLINE Trigger raw data selection
///        and  conistency checks
///
/// Related task: https://alice.its.cern.ch/jira/browse/PWGPP-6 and PWGPP-134


/// \author Marian Ivanov   - marian.ivanov@cern.ch
/// \author Mesut Arslandok, Mikolaj  - older versions (rawmerege.C)


class AliOfflineTrigger : public TNamed {
public:
  AliOfflineTrigger(const char *triggerName, Int_t timeOut=30, Int_t cacheSize=500000000);
  virtual ~AliOfflineTrigger(){;}
  void     DumpGIDRAWReader(const char *rawFile="raw.root");
  void     DumpGIDRAWTree(const char *rawFile="raw.root");
  void     DumpGIDESD(const char * chinput="AliESDs.root", const char *trigger="1", const char *choutput="gidesd.list");
  void     TestDiffGIDList();
  //  ESD trigger setup
  void     AddESDAlias(const char *aliasName, const char *aliasValue );
  void     SetTriggerAlias(TTree* tree, const char * trigger);
  static   TTree*    MakeDiffTree(const char *refTree, const char *friendTrees);
  //  Filtered trees trigger 
  void  ExtractSelected(const char *rawList, const char * triggerList, const char * outputName, Long_t maxCounter,  Int_t verbose=1);
  void  LoadTriggerList(const char * triggerList);
  Int_t LoadMapFromRawData(const char *rawFile="raw.root", Int_t verbose=1);
  void  ExtractSelected(const char *rawFile="raw.root", Int_t verbose=1);
public:
  //private: 
  std::map<ULong64_t, TString>  fTrgGIDChunkName;  /// GID -> ChunkName 
  std::map<ULong64_t, TString>  fTrgGIDTrigger;    /// GID -> Trigger type 
  std::map<ULong64_t, UInt_t>   fTrgGIDEventNr;    /// GID -> EventNumber
  std::map<ULong64_t, UInt_t>   fTrgGIDTimeStamp;  /// triger map GID -> TimeStamp map
  //
  Int_t fCounterFileOutput;                        ///  input file counter
  Int_t fCounterEventOutput;                       ///  input event counter
  Int_t fCounterFileInput;                         ///  input file counter
  Int_t fCounterEventInput;                        ///  input event counter
  TString fRawName;                                ///  name of the output file
  TFile * fRawTriggerFile;                         ///! pointer to ouput  raw trigger files
  TTree * fRawTriggerTree;                         ///! pointer to output raw trigger tree
  std::map<ULong64_t, TString>  fRAWGIDChunkName;  ///  GID -> ChunkName 
  std::map<ULong64_t, UInt_t>   fRAWGIDEventNr;    ///  GID -> EventNumber
  std::map<UInt_t,ULong64_t>    fRAWEventNrGID;    ///  EventNumber -> GID
  std::map<ULong64_t, UInt_t>   fRAWGIDTimeStamp;  ///  triger map GID -> TimeStamp map
  //
  Int_t   fDefaultTimeOut;         // default time out in seconds for reading and file opening
  Int_t   fDefaultTreeCache;       // default tree cache size
  TObjArray *fESDTriggerList;      // esd trigger list (set of aliases)
  ClassDef(AliOfflineTrigger, 1);  // interface to various trending trees
};

#endif // ALIOFFLINETRIGGER_H
