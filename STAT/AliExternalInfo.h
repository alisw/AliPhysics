#ifndef ALIEXTERNALINFO_H
#define ALIEXTERNALINFO_H

#include <vector>
#include <map>
#include "TTree.h"
#include "TObject.h"

class TChain;
class TString;
//class Period;

/*
  This class gives an interface to different trees.

  To get the pointer to a tree use GetTree*(<id>).
  Cache*(<id>) downloads the necessary information.

  For the correct usage of this class you need to have a certificate which does not need typing in your password.
*/

class AliExternalInfo
{
public:
  AliExternalInfo();
  AliExternalInfo(TString GlobalPath);
  ~AliExternalInfo();

  Bool_t CacheMC();
  Bool_t CacheRCT(TString period, TString pass); // Calls Cache(...) with the corresponding type
  Bool_t CacheLogbook(TString period);
  Bool_t CacheTriggerClasses(TString year);
  Bool_t CacheDataQA(TString detector, TString period, TString pass);
  Bool_t CacheProdCycle(); // downloads the master page
  Bool_t CacheProdPasses();
  Bool_t CacheProdCycle(TString id); // can be used, if someone knows the id at mona lisa
  // Bool_t CacheProdCycle(TString period, TString pass); // NOT YET IMPLEMENTED
  // Bool_t CacheAll(); // NOT YET IMPLEMENTED


  TTree* GetTreeMC(TString period = "", TString anchorYear = "", TString productionTag = "");
  TTree* GetTreeRCT(TString period, TString pass); // Calls GetTree(...) with the corresponding type
  TTree* GetTreeLogbook(TString period);
  TTree* GetTreeTriggerClasses(TString year);
  TTree* GetTreeDataQA(TString detector, TString period, TString pass);
  TTree* GetTreeProdCycle(); // downloads master page
  TTree* GetTreeProdPasses(); // downloads master page
  TTree* GetTreeProdCycle(TString id); // can be used, if someone knows the id at mona lisa
  // TTree* GetTreeProdCycle(TString period, TString pass);  // NOT YET IMPLEMENTED


  TChain* GetChainMC(TString period = "", TString anchorYear = "", TString productionTag = "");
  TChain* GetChainRCT(TString period, TString pass);
  TChain* GetChainLogbook(TString period);
  TChain* GetChainTriggerClasses(TString year);
  TChain* GetChainDataQA(TString detector, TString period, TString pass);
  TChain* GetChainProdCycle(TString id); // can be used, if someone knows the id at mona lisa
  // TChain* GetChainProdCycle(TString period, TString pass);  // NOT YET IMPLEMENTED

  TString GetLocalLocation(){return fGlobalDir; }
  void SetLocalLocation(TString location){fGlobalDir = location; }

  void SetGlobalTimeLimit(long int t){ fTimeLimit = t;};
  Double_t GetGlobalTimeLimit(){return fTimeLimit;}

  void CreateMapWithPassesAndPeriods();

  TChain* GetChain(TString period, TString pass, TString anchorYear, TString productionTag, Int_t type);
  TTree* GetTree(TString period, TString pass, TString anchorYear, TString productionTag, Int_t type);
  Bool_t Cache(TString period, TString pass, Int_t type); // Downloads the tree in the working directory

  enum type{kRCT, kMC, kLogbook, kTriggerClasses, kTPC, kEVS, kTOF, kProdCycle, kProdCycleCpasses, kRawQATPC, kAll};

private:
  Bool_t IsDownloadNeeded(TString file);
  Bool_t CheckPeriod(TString period);
  Bool_t CheckPass(TString pass);

  void SetUpForPeriodPass(const TString& period, const TString& pass, TString& localpath, TString& file, const TString& filename, const TString& filetype, Int_t type);

  TString wget(const TString& localpath, const TString& filename, const TString& filetype, const TString& remotepath);
  TString GetYearFromPeriod(TString period);

  TString fCertificate;
  TString fPrivateKey;

  TString fGlobalDir;
  TString fFormMCFilename;
  TString fFormMCDir;
  TString fFormMCFiletype;
  TString fFormRCTFilename;
  TString fFormRCTDir;
  TString fFormRCTFiletype;
  TString fFormLogbookFilename;
  TString fFormLogbookDir;
  TString fFormLogbookFiletype;
  TString fFormTriggerClassesFilename;
  TString fFormTriggerClassesDir;
  TString fFormTriggerClassesFiletype;
  TString fFormTrendingFilename;
  TString fFormTrendingDir;
  TString fFormTrendingFiletype;
  TString fFormProdCycleFilename;
  TString fFormProdCycleDir;
  TString fFormProdCycleFiletype;
  TString fFormProdPassesFilename;
  TString fFormProdPassesDir;
  TString fFormProdPassesFiletype;

  long int fTimeLimit; // in seconds

  TTree* fTree;

  class Period
  {
    public:
      Period(){};
      Period(TString name) : fName(name){};
      ~Period(){};

      void SetPeriod(TString name){
        fName = name;
      }
      TString GetPeriod() {
        return fName;
      }
      void AddPass(TString pass){
        fPasses.push_back(pass);
      }
      std::vector<TString> GetPasses(){
        return fPasses;
      }
      void PrintPeriodPass(){
        std::cout << "Period: " << GetPeriod() << "   Passes: ";
        for (unsigned int i = 0; i < fPasses.size(); ++i){
          std::cout << fPasses.at(i) << " ";
        }
        std::cout << std::endl;
      }
    private:
      TString fName;
      std::vector<TString> fPasses;
  };

  std::vector<Period> fPeriods;

  ClassDef(AliExternalInfo, 0);  // interface to various trending trees
};

#endif // ALIEXTERNALINFO_H


