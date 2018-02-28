// -*- C++ -*-

#include <sstream>

#include "AliVdMScanData.h"

ClassImp(AliVdMScanData);

void AliVdMScanData::FillDefaultBranches(const AliVdMMetaData& vdmMetaData,
                                         TTree *VdM,
                                         const std::vector<std::string>& triggerNames)
{
  const Double_t deltaT_tree = 2.0;
  const AliTriggerBCMask& bcMask = vdmMetaData.GetTriggerBCMask();

  Double_t timeOfCounters=0;
  VdM->SetBranchAddress("time", &timeOfCounters);
  std::vector<TVectorD> counters     (triggerNames.size(), TVectorD(3564));
  std::vector<TVectorD> sumOfCounters(triggerNames.size(), TVectorD(3564));

  VdM->SetBranchStatus("*", 0);
  VdM->SetBranchStatus("time", 1);
  for (Int_t i=0, n=triggerNames.size(); i<n; ++i) {
    VdM->SetBranchStatus(triggerNames[i].c_str(), 1);
    VdM->SetBranchAddress(triggerNames[i].c_str(), counters[i].GetMatrixArray());
  }

  Int_t counter = 0;
  for (auto it=vdmMetaData.GetScansBegin(), iend=vdmMetaData.GetScansEnd(); it!=iend; ++it) {
    const AliXMLEngine::Node xmlNode(*it);
    if (TString(xmlNode.GetAttr("type").GetData()) == "LSC")
      continue;
    Printf("%s %d",
           AliVdMMetaData::GetScanName(xmlNode).Data(),
           AliVdMMetaData::GetScanType(xmlNode));
    Printf("counter=%d", counter);
    const Int_t iScanType = AliVdMMetaData::GetScanType(xmlNode);

    for (const std::string& triggerName : triggerNames)
      CreateDefaultBranches(counter, triggerName);

    TTree tSep;
    std::istringstream iss(xmlNode.GetData());
    tSep.ReadStream(iss);

    AliVdMTree::DefaultBranchData defBranchData;
    defBranchData.fSep(1-iScanType) = 1e-3*atof(xmlNode.GetAttr("offset_um").GetData()); // mu -> mm
    Double_t timeStart=0, timeEnd=0;
    tSep.SetBranchAddress("timeStart", &timeStart);
    tSep.SetBranchAddress("timeEnd",   &timeEnd);
    tSep.SetBranchAddress("sep",       defBranchData.fSep.GetMatrixArray() + iScanType);
    for (Int_t i=0, n=tSep.GetEntries(); i<n; ++i) {
      tSep.GetEntry(i);
      const Long64_t nSelected = VdM->Draw("Entry$", TString::Format("time>%f && time<%f", timeStart, timeEnd), "GOFF");
      const Double_t      *idx = VdM->GetV1();
      if (nSelected < 2)
        continue;
      std::fill(sumOfCounters.begin(),  sumOfCounters.end(), 0);
      defBranchData.fTime = 0;
      for (Int_t j=idx[0], m=idx[nSelected-1]; j<m; ++j) {
        VdM->GetEntry(j);
        if (j == idx[0])
          defBranchData.StartTime() = timeOfCounters - 0.5*deltaT_tree;
        defBranchData.EndTime() = timeOfCounters + 0.5*deltaT_tree;
        std::transform(counters.begin(), counters.end(), sumOfCounters.begin(), sumOfCounters.begin(), std::plus<TVectorD>());
      }

      for (Int_t j=0, m=triggerNames.size(); j<m; ++j) {
        AliVdMTree& t = GetMap(counter)[triggerNames[j]];
        for (Int_t bc=0; bc<3564; ++bc) {
          if (bcMask.GetMask(bc))
            continue;
          defBranchData.BCID() = bc;
          for (Int_t deltaBC=0; deltaBC<30; ++deltaBC) {
            const Int_t _bc = ((3564+bc-deltaBC) % 3564);
            defBranchData.Counter(deltaBC) = sumOfCounters[j][_bc];
          }
          t.FillDefaultBranches(defBranchData);
        }
      }
    }
    ++counter;
  }
}
