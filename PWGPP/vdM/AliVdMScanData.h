// -*- C++ -*-

#ifndef PWGPP_VDM_ALIVDMSCANDATA_H
#define PWGPP_VDM_ALIVDMSCANDATA_H

#include <vector>
#include <map>

#include <TObject.h>
#include <TTree.h>
#include <TString.h>

#include "AliVdMTree.h"
#include "AliVdMMetaData.h"

// -----------------------------------------------------------------------------
// this class holds data for all scans and all trigger combinations
class AliVdMScanData : public TObject {
public:
  typedef std::map<std::string, AliVdMTree> map_type;
  typedef std::vector<map_type> vector_type;

  AliVdMScanData()
    : TObject()
    , fData() {}

  virtual ~AliVdMScanData() {}

  TTree* GetTTree(Int_t iScan, const char* trigName) {
    return GetMap(iScan)[trigName].GetTTree();
  }

  Int_t GetNScans() const { return fData.size(); }

  AliVdMScanData& FillDefaultBranches(const AliVdMMetaData& vdmMetaData, TTree *t,
                                      const std::vector<std::string>& triggerNames);
  AliVdMScanData& FillDefaultBranchesFromCTPScalers(const AliVdMMetaData& vdmMetaData, TTree *t,
                                                    const std::vector<std::string>& triggerNames);

  map_type& GetMap(std::size_t iScan) {
    for (; iScan >= fData.size();)
      fData.emplace_back(map_type());
    return fData[iScan];
  }

protected:
  TTree* CreateDefaultBranches(Int_t iScan, std::string trigName) {
    const TString treeName = MakeTreeName(iScan, trigName);
    return GetMap(iScan)[trigName].CreateDefaultBranches(treeName);
  }

  static TString MakeTreeName(Int_t iScan, TString trigName) {
    return TString::Format("T%s_scan%d", trigName.Data(), iScan);
  }
private:
  vector_type fData; // vector<map<string, VdMTree> >

  ClassDef(AliVdMScanData, 1);
} ;

#endif // PWGPP_VDM_ALIVDMSCANDATA_H
