// -*- C++ -*-

#ifndef PWGPP_VDM_ALIVDMSCANDATA_H
#define PWGPP_VDM_ALIVDMSCANDATA_H

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

  enum { MaxNumScans = 6 };

  AliVdMScanData()
    : TObject()
    , fData() {}

  virtual ~AliVdMScanData() {}

  TTree* GetTTree(Int_t iScan, const char* trigName) {
    return fData[iScan][trigName].GetTTree();
  }

  void FillDefaultBranches(const AliVdMMetaData& vdmMetaData, TTree *t,
                           const std::vector<std::string>& triggerNames);

  map_type& GetMap(Int_t iScan) { return fData[iScan]; }

protected:
  TTree* CreateDefaultBranches(Int_t iScan, std::string trigName) {
    const TString treeName = MakeTreeName(iScan, trigName);
    return fData[iScan][trigName].CreateDefaultBranches(treeName);
  }

  static TString MakeTreeName(Int_t iScan, TString trigName) {
    return TString::Format("T%s_scan%d", trigName.Data(), iScan);
  }
private:
  map_type fData[MaxNumScans]; // max number of scans = MaxNumScans

  ClassDef(AliVdMScanData, 1);
} ;

#endif // PWGPP_VDM_ALIVDMSCANDATA_H
