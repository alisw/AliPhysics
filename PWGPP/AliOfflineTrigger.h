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
  //  Filtered trees trigger setup
public:
  Int_t   fDefaultTimeOut;         // default time out in seconds for reading and file opening
  Int_t   fDefaultTreeCache;       // default tree cache size
  TObjArray *fESDTriggerList;      // esd trigger list (set of aliases)
  ClassDef(AliOfflineTrigger, 1);  // interface to various trending trees
};

#endif // ALIOFFLINETRIGGER_H
