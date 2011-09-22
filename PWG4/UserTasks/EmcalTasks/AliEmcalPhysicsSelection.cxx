/* $Id$ */

#include "AliEmcalPhysicsSelection.h"
#include "AliESDEvent.h"
#include <AliLog.h>

ClassImp(AliEmcalPhysicsSelection)

AliEmcalPhysicsSelection::AliEmcalPhysicsSelection() : 
  AliPhysicsSelection(), 
  fMarkFastOnly(0), 
  fMarkLedEvent(0),
  fSkipFastOnly(0), 
  fSkipLedEvent(0),
  fCellMinE(2), 
  fClusMinE(5), 
  fIsFastOnly(0),
  fIsLedEvent(0),
  fIsGoodEvent(0),
  fCellMaxE(0), 
  fClusMaxE(0) 
{
  // Default constructor.
}

//__________________________________________________________________________________________________
UInt_t AliEmcalPhysicsSelection::GetSelectionMask(const TObject* obj) 
{ 
  const AliESDEvent *ev = static_cast<const AliESDEvent*>(obj);
  UInt_t res = IsCollisionCandidate(ev); 

  fIsFastOnly  = kFALSE;
  fIsGoodEvent = kFALSE;
  fIsLedEvent  = kFALSE;
  fCellMaxE    = 0;
  fClusMaxE    = 0;

  if ((res & AliVEvent::kAnyINT) || 
      (res & AliVEvent::kEMC1)   || 
      (res & AliVEvent::kEMC7))
    fIsGoodEvent = kTRUE;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  am->LoadBranch("EMCALCells.");
  am->LoadBranch("CaloClusters");

  AliESDCaloCells *cells = ev->GetEMCALCells();
  const Short_t nCells   = cells->GetNumberOfCells();
    
  // mark LHC11a fast only partition if requested 
  if (res & AliVEvent::kFastOnly) {
    fIsFastOnly = kTRUE;
    if (fMarkFastOnly||fSkipFastOnly) {
      if (nCells>0) {
        AliFatal(Form("Number of cells %d, even though EMCAL should not be in fast only partition.",nCells));
      }
      fIsGoodEvent = kFALSE;
    }
  }

  // count cells above threshold
  Int_t nCellCount[10] = {0,0,0,0,0,0,0,0,0,0};
  for(Int_t iCell=0; iCell<nCells; ++iCell) {
    Short_t cellId = cells->GetCellNumber(iCell);
    Double_t cellE = cells->GetCellAmplitude(cellId);
    Int_t sm       = cellId / (24*48);
    if (cellE>0.1)
      ++nCellCount[sm];
    if (cellE>fCellMaxE)
      fCellMaxE = cellE;
  }
  
  const Int_t nCaloClusters = ev->GetNumberOfCaloClusters();
  for(Int_t iClus = 0; iClus<nCaloClusters; ++iClus) {
    AliESDCaloCluster *cl = ev->GetCaloCluster(iClus);
    if (!cl->IsEMCAL()) 
      continue;
    Double_t e = cl->E();
    if (e>fClusMaxE)
      fClusMaxE = e;
  }

  // bad cell criterion for LHC11a from 
  // https://indico.cern.ch/materialDisplay.py?contribId=4&materialId=slides&confId=147067
  const Int_t runN = ev->GetRunNumber();
  if ((runN>=144871) && (runN<=146860)) { 
    if (nCellCount[4] > 100)
      fIsLedEvent = kTRUE;
    else {
      if ((runN>=146858) && (runN<=146860)) {
        if ((res&AliVEvent::kMB) && (nCellCount[3]>=21))
          fIsLedEvent = kTRUE;
        else if ((res&AliVEvent::kEMC1) && (nCellCount[3]>=35))
          fIsLedEvent = kTRUE;
      }
    }
    if (fIsLedEvent) {
      fIsGoodEvent = kFALSE;
    }
  }

  if (fCellMaxE>=fCellMinE)
    res |= kEmcalHC;

  if (fClusMaxE>=fClusMinE)
    res |= kEmcalHT;

  if (fIsGoodEvent)
    res |= kEmcalOk;

  if ((fSkipLedEvent && fIsLedEvent) ||
      (fSkipFastOnly && fIsFastOnly))
    res = 0;

  return res;
}
