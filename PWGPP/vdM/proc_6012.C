// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>

#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"
#include "AliVdMPileup.h"

#include "proc.h"
#include "proc_pileup.h"
#include "proc_bkgd.h"

// -----------------------------------------------------------------------------
//  analysis for fill 6012

AliVdMScanData vdmScanData;

void proc_6012()
{
  const std::vector<std::string> triggerNames{
    "c2TVX",
    "c2VBAandVBC",
    "c2UBAandUBC",

    "c2T0AandNotT0C",
    "c2T0CandNotT0A",

    "c2VBAandNotVBC",
    "c2VBCandNotVBA",

    "c2UBAandNotUBC",
    "c2UBCandNotUBA"
  };

  const AliVdMMetaData vdmMetaData(AliVdMMetaData::GetFileName("6012/6012.xml"));

  TFile::Open(AliVdMMetaData::GetFileName("6012/vdm_time_6012_5m11.5_11p17.5_1_v3.root"), "READ");
  gDirectory->ls();
  TTree * VdM = nullptr;
  gFile->GetObject("VdM", VdM);
  VdM->AddFriend("DDL2", AliVdMMetaData::GetFileName("6012/vdm_DDL2_6012-5m11.5_11p17.5_1_v3.root"));

  vdmScanData.FillDefaultBranches(vdmMetaData, VdM, triggerNames);
  proc(vdmMetaData, vdmScanData, triggerNames);

#if 1
  // (3) determine pile-up
  proc_pileup(vdmMetaData, vdmScanData,
              "c2UBAandUBC", "c2UBAandNotUBC", "c2UBCandNotUBA",
              {0.49,0.17, 1.5e-5,1.5e-5});

  proc_pileup(vdmMetaData, vdmScanData,
              "c2VBAandVBC", "c2VBAandNotVBC", "c2VBCandNotVBA",
              {0.070,0.086, 4.2e-5,4.6e-5});

  proc_pileup(vdmMetaData, vdmScanData,
              "c2TVX", "c2T0AandNotT0C", "c2T0CandNotT0A",
              {0.57,0.63, 2.2e-5,3.8e-5});
#endif

//  proc_bkgd(vdmMetaData, vdmScanData, "c2UBAandUBC");
}

