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
//  analysis for fill 6380

AliVdMScanData vdmScanData;

void proc_4634()
{
  const std::vector<std::string> triggerNames{
    "C0TVX",
    "CINT7",
    "CADAND",

    "CT0ANOTC",
    "CT0CNOTA",

    "CVBANOTC",
    "CVBCNOTA",

    "CADANOTC",
    "CADCNOTA"};

  const AliVdMMetaData vdmMetaData(AliVdMMetaData::GetFileName("4634/4634.xml"));

  TFile::Open(AliVdMMetaData::GetFileName("4634/vdm_time_4634_5.5m11.5_11.5p17.5_1_v3x.root"), "READ");
  gDirectory->ls();
  TTree * VdM = nullptr;
  gFile->GetObject("VdM", VdM);
  // VdM->AddFriend("DDL2", AliVdMMetaData::GetFileName("4634/vdm_DDL2_4634-4m12_10p18.root"));

  TFile::Open(AliVdMMetaData::GetFileName("4634/scalers_244369.root"), "READ");
  TTree *TS = nullptr;
  gFile->GetObject("TS", TS);

  vdmScanData.FillDefaultBranchesFromCTPScalers(vdmMetaData, TS, triggerNames);
  const Bool_t computeBkgd = kFALSE;
  proc(vdmMetaData, vdmScanData, triggerNames, computeBkgd);

#if 1
  // (3) determine pile-up
  proc_pileup(vdmMetaData, vdmScanData,
              "CADAND", "CADANOTC", "CADCNOTA",
              {0.2,0.2, 1.5e-5,1.5e-5}, computeBkgd);

  proc_pileup(vdmMetaData, vdmScanData,
              "CINT7", "CVBANOTC", "CVBCNOTA",
              {0.070,0.086, 4.2e-5,4.6e-5}, computeBkgd);

  proc_pileup(vdmMetaData, vdmScanData,
              "C0TVX", "CT0ANOTC", "CT0CNOTA",
              {0.57,0.63, 2.2e-5,3.8e-5}, computeBkgd);
#endif

//  proc_bkgd(vdmMetaData, vdmScanData, "c2UBAandUBC");
}

