// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>

#include "AliVdMMetaData.h"
#include "AliVdMScanData.h"
#include "AliVdMPileup.h"

#include "proc.h"
#include "proc_pileup.h"

// -----------------------------------------------------------------------------
//  analysis for fill 4937

AliVdMScanData vdmScanData;

void proc_4937()
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

  const AliVdMMetaData vdmMetaData(AliVdMMetaData::GetFileName("4937/4937.xml"));

  TFile::Open(AliVdMMetaData::GetFileName("4937/vdm_time_4937_6m11_12p17_1_v3.root"), "READ");
  TTree * VdM = nullptr;
  gFile->GetObject("VdM", VdM);
  VdM->AddFriend("DDL2", AliVdMMetaData::GetFileName("4937/vdm_DDL2_4937.root"));

  vdmScanData.FillDefaultBranches(vdmMetaData, VdM, triggerNames);
  proc(vdmMetaData, vdmScanData, triggerNames);

  // determine pile-up
  proc_pileup(vdmMetaData, vdmScanData,
              "c2UBAandUBC", "c2UBAandNotUBC", "c2UBCandNotUBA",
              {0.24,0.13, 4.4e-5,2.3e-5});

  proc_pileup(vdmMetaData, vdmScanData,
              "c2VBAandVBC", "c2VBAandNotVBC", "c2VBCandNotVBA",
              {0.061,0.076, 14e-5,3.5e-5});

  proc_pileup(vdmMetaData, vdmScanData,
              "c2TVX", "c2T0AandNotT0C", "c2T0CandNotT0A",
              {0.39,0.44, 13e-5,3.3e-5});
}

