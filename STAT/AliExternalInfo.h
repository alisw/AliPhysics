#ifndef ALIEXTERNALINFO_H
#define ALIEXTERNALINFO_H

/// \class AliExternalInfo
/// \brief This class gives you an interface to different trees of information spread throughout ALICE
///
/// Related task: https://alice.its.cern.ch/jira/browse/ATO-46
/// With this class you can easily download and use information from different
/// resources in ALICE, like the RCT, logbook, QA...
///
/// For the correct usage of this class you need to have a certificate which does not need typing in your password.
/**

Examples:

1) Get the tree from a QA tree

AliExternalInfo b(".", "config.cfg"); // first parameter is the working directory. Second parameter is the configuration file.
b.Cache("QA.TPC", "LHC15f", "pass2"); // Actually not needed because "GetTree(...)" chechs if the resource has been if this is not the case it will automatically call Cache(...)
TTree* treeTPCQA15fPass2 = b.GetTree("QA.TPC", "LHC15f", "pass2"); // Returns a tree with the information of the resource
// TTree* treeTPCQA15fPass2 = b.GetTreeDataQA("TPC", "LHC15f", "pass2"); // Different method to get the tree. Calls internally the "Cache(...)" function
treeTPCQA15fPass2->Print(); // Do something with the tree

2a) Create a chain
b.Cache("QA.EVS", "LHC10b", "pass4"); Downloads the resource
b.Cache("QA.EVS", "LHC10c", "pass4");
b.Cache("QA.EVS", "LHC10d", "pass4");
TChain* chain10evs = b.GetChain("QA.EVS", "LHC10*", "pass4"); // Creates a chain of all downloaded resources. Note the '*' which is used like in ls
chain10evs->Draw("interactionRate:run", "", "*"); // Draws a histogram interaction rate vs run

2b) Create a chain
TTree * treeProdCycle7076 = b.GetTreeProdCycleByID("7076");
TTree * treeProdCycle7236 = b.GetTreeProdCycleByID("7236");
TTree * treeProdCycle7214 = b.GetTreeProdCycleByID("7214");
TChain* chainProdCycleIDs = b.GetChainProdCycleByID("[0-9]*"); // needs number wildcards because a ProdCycle.root is also available
chainProdCycleIDs->Print();

3) Using the friends tree
TTree* treeTPCQA15fPass1 = b.GetTree("QA.TPC", "LHC15f", "pass2");
TTree* treeEVSQA15fPass1 = b.GetTree("QA.EVS", "LHC15f", "pass2");
// TTree* treeTPCQA15hPass1 = b.GetTree("QA.TPC", "LHC15h", "pass1"); // Only one single tree of type can be added to the friends tree
// TTree* treeEVSQA15hPass1 = b.GetTree("QA.EVS", "LHC15h", "pass1"); // Only one single tree of type can be added to the friends tree
b.GetFriendsTree()->Draw("TPC.tpcItsMatchA:EVS.interactionRate","TPC.meanMult>0", "*");
c1.SaveAs("c1.png");
*/
/// \author Carsten Klein <Carsten.Klein@cern.ch>, Goethe Universität Frankfurt
/// \date Jan 18, 2016

#include <map>

#include "TString.h"

// forward declarations
class TTree;
class TChain;

class AliExternalInfo : public TObject {
public:
  AliExternalInfo (TString localStorageDirectory = ".", TString configLocation = "$ALICE_ROOT/STAT/Macros/AliExternalInfo.cfg"/*, Bool_t copyToLocalStorage = kTRUE*/);
  virtual ~AliExternalInfo();
  void ReadConfig();
  void PrintConfig();
  Bool_t Cache(TString type="", TString period = "", TString pass=""); // Downloads the tree in the working directory
  Bool_t CacheMC()                                                      {return Cache("MonALISA.MC", "", "");}
  Bool_t CacheRCT(TString period, TString pass)                         {return Cache("MonALISA.RCT", period, pass);}
  Bool_t CacheDataQA(TString detector, TString period, TString pass)    {return Cache("QA." + detector, period, pass);}
  Bool_t CacheLogbook(TString period)                                   {return Cache("Logbook", period, "");}
  Bool_t CacheTriggerClasses(TString period)                            {return Cache("TriggerClasses", period, "");}
  Bool_t CacheProdCycle()                                               {return Cache("MonALISA.ProductionCycle", "", "");}
  Bool_t CacheCPass()                                                   {return Cache("MonALISA.ProductionCPass", "", "");}
  Bool_t CacheProdCycleByID(TString ID)                                 {return Cache("MonALISA.ProductionCycleID", ID, "");}

  TTree* GetTree(TString type, TString period, TString pass);
  TTree* GetTree(TString type, TString period, TString pass, TString friendList);
  TTree* GetTreeMC()                                                    {return GetTree("MonALISA.MC", "", "");}
  // TTree* GetTreeMC(TString period = "", TString anchorYear = "", TString productionTag = "") {return GetTree("MonALISA.MC", "", "");} // deprecated; not supported anymore
  TTree* GetTreeRCT(TString period, TString pass)                       {return GetTree("MonALISA.RCT", period, pass);}
  TTree* GetTreeDataQA(TString detector, TString period, TString pass)  {return GetTree("QA." + detector, period, pass);}
  TTree* GetTreeLogbook(TString period)                                 {return GetTree("Logbook", period, "");}
  TTree* GetTreeTriggerClasses(TString period)                          {return GetTree("TriggerClasses", period, "");}
  TTree* GetTreeProdCycle()                                             {return GetTree("MonALISA.ProductionCycle", "", "");}
  TTree* GetTreeCPass()                                                 {return GetTree("MonALISA.ProductionCPass", "", "");}
  TTree* GetTreeProdCycleByID(TString ID)                               {return GetTree("MonALISA.ProductionCycleID", ID, "");}
  TTree*  GetCPassTree(const char * period, const  char *pass); 

  TChain* GetChain(TString type, TString period, TString pass);
  TChain* GetChainMC()                                                  {return GetChain("MonALISA.MC", "", "");}
  TChain* GetChainRCT(TString period, TString pass)                     {return GetChain("MonALISA.RCT", period, pass);}
  TChain* GetChainDataQA(TString detector, TString period, TString pass){return GetChain("QA." + detector, period, pass);}
  TChain* GetChainLogbook(TString period)                               {return GetChain("Logbook", period, "");}
  TChain* GetChainTriggerClasses(TString period)                        {return GetChain("TriggerClasses", period, "");}
  TChain* GetChainProdCycle()                                           {return GetChain("MonALISA.ProductionCycle", "", "");}
  TChain* GetChainProdCycleByID(TString ID)                             {return GetChain("MonALISA.ProductionCycleID", ID, "");}

  TTree* GetFriendsTree() const {return fTree;}
  TChain* GetFriendsChain() const {return fChain;} // _Not_ working properly!!!

  void     SetMaxCacheSize(Long64_t size) { fMaxCacheSize=size;   }
  Long64_t GetMaxCacheSize() const        { return fMaxCacheSize; }

  static const TString& GetDefaultConfig() { return fgkDefaultConfig; }
  static void BuildHashIndex(TTree* tree, const char *chbranchName,  const char *chindexName);
private:
  Bool_t AddTree(TTree* tree, TString type);
  Bool_t AddChain(TString type, TString period, TString pass);
  void SetupVariables(TString& internalFilename, TString& internalLocation, Bool_t& resourceIsTree, TString& pathStructure, \
                      TString& detector, TString& rootFileName, TString& treeName, const TString& type, const TString& period, const TString& pass, TString &indexName);
  const TString GetYearFromPeriod(const TString& period);
  const TString Wget(TString& mifFilePath, const TString& internalLocation, TString rootFileName, const TString& externalLocation);
  const TString CreatePath(TString type, TString period, TString pass);
  Bool_t IsDownloadNeeded(TString file, TString type);

  // Bool_t fCopyDataToLocalStorage;
  TString fConfigLocation;                          ///< location of the config file
  TString fLocalStorageDirectory;                   ///< location of the local cache directory
  std::map<TString, TString> fLocationTimeOutMap;   ///< map with configuration parameters
  TTree* fTree;                                     ///< master tree with friends
  TChain* fChain;                                   ///< master chain with friends
  std::map<TString, TChain*> fChainMap;             ///< map of chains
  Long64_t fMaxCacheSize;                           ///< maximum chache size for trees and chains

  static const TString fgkDefaultConfig;            ///< default config file

  ClassDef(AliExternalInfo, 0);  // interface to various trending trees

};

#endif // ALIEXTERNALINFO_H
