/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "TSystem.h"

#include "AliProdInfo.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliOADBContainer.h"

#include "AliTimeRangeCut.h"

ClassImp(AliTimeRangeCut)

//______________________________________________________________________________
void AliTimeRangeCut::InitFromEvent(const AliVEvent* event)
{
  const Int_t run = event->GetRunNumber();

  InitFromRunNumber(run);
}

//______________________________________________________________________________
void AliTimeRangeCut::InitFromRunNumber(const Int_t run)
{
  if (run == fLastRun) return;
  fLastRun = run;

  // analysis manager
  const AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();

  // ===| set OADB path |===
  if (fOADBPath.IsNull()) {
    if (mgr) {
      fOADBPath = mgr->GetOADBPath();
    }
    else {
      fOADBPath = gSystem->Getenv("ALICE_PHYSICS");
      fOADBPath += "/OADB";
    }
  }

  // ===| set pass name |===
  TString passName = "pass1";
  if (mgr) {
    AliVEventHandler *inputHandler = mgr->GetInputEventHandler();
    if (inputHandler) {
      TList *uiList = inputHandler->GetUserInfo();
      AliProdInfo prodInfo(uiList);

      passName = prodInfo.GetRecoPassName();

      AliInfoF("Setting pass name from ProdInfo: %s", passName.Data());
    }
  }
  printf("pass: %s\n", passName.Data());

  // ===| Get the AliTimeRangeMasking object |===
  AliOADBContainer cont("TimeRangeMasking");
  cont.InitFromFile(Form("%s/COMMON/PHYSICSSELECTION/data/TimeRangeMasking.root", fOADBPath.Data()), "TimeRangeMasking");

  fTimeRangeMasking = (AliTimeRangeMasking<ULong64_t, UShort_t>*)cont.GetObject(run, "", passName);

}

//______________________________________________________________________________
UShort_t AliTimeRangeCut::GetMask(const AliVEvent* event) const
{
  if (!fTimeRangeMasking) return 0;
  if (!event) return 0;

  AliVHeader* header = event->GetHeader();
  if (!header) return 0;

  const ULong64_t gid = header->GetEventIdAsLong();
  return GetMask(gid);
}

//______________________________________________________________________________
UShort_t AliTimeRangeCut::GetMask(const ULong64_t gid) const
{
  AliTimeRangeMask<ULong64_t, UShort_t>* range = fTimeRangeMasking->FindTimeRangeMask(gid);

  if (!range) return 0;
  return range->GetMaskReasons();
}

//______________________________________________________________________________
Bool_t AliTimeRangeCut::CutEvent(const AliVEvent* event, const UShort_t mask/* = 0*/) const
{
  const UShort_t maskReasons = GetMask(event);

  if (maskReasons == 0) return kFALSE;
  if (mask == 0) return kTRUE;

  return (maskReasons & mask);
}

//______________________________________________________________________________
Bool_t AliTimeRangeCut::CutEvent(const ULong64_t gid, const UShort_t mask/* = 0*/) const
{
  const UShort_t maskReasons = GetMask(gid);

  if (maskReasons == 0) return kFALSE;
  if (mask == 0) return kTRUE;

  return (maskReasons & mask);
}
