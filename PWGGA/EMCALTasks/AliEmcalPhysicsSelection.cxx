// $Id$
//
// Emcal physics selection class.
//
// Author: C.Loizides

#include "AliEmcalPhysicsSelection.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliLog.h"

ClassImp(AliEmcalPhysicsSelection)

AliEmcalPhysicsSelection::AliEmcalPhysicsSelection() : 
  AliPhysicsSelection(), 
  fMarkFastOnly(0), 
  fMarkLedEvent(0),
  fSkipFastOnly(0), 
  fSkipLedEvent(0),
  fCellMinE(-1), 
  fClusMinE(-1),
  fTrackMinPt(-1), 
  fTriggers(0),
  fIsFastOnly(0),
  fIsLedEvent(0),
  fIsGoodEvent(0),
  fCellMaxE(0), 
  fClusMaxE(0),
  fTrackMaxPt(0)
{
  // Default constructor.
}

//__________________________________________________________________________________________________
UInt_t AliEmcalPhysicsSelection::GetSelectionMask(const TObject* obj) 
{ 
  // Calculate selection mask.
  
  const AliVEvent *ev   = dynamic_cast<const AliVEvent*>(obj);
  if (!ev) {
    return 0;
  }

  UInt_t res = 0;
  const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(obj);
  const AliAODEvent *aev = 0;
  if (eev) {
    res = IsCollisionCandidate(eev); 
  } else {
    aev = dynamic_cast<const AliAODEvent*>(obj);
    res = aev->GetHeader()->GetOfflineTrigger();
  }

  if (fTriggers) { // only process given triggers
    if (res & fTriggers == 0)
      return res;
  }

  fIsFastOnly  = kFALSE;
  fIsGoodEvent = kFALSE;
  fIsLedEvent  = kFALSE;
  fCellMaxE    = 0;
  fClusMaxE    = 0;
  fTrackMaxPt  = 0;

  if ((res & AliVEvent::kAnyINT) || 
      (res & AliVEvent::kEMC1)   || 
      (res & AliVEvent::kEMC7)   || 
      (res & AliVEvent::kEMCEJE) || 
      (res & AliVEvent::kEMCEGA))
    fIsGoodEvent = kTRUE;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  AliVCaloCells *cells   = ev->GetEMCALCells();
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

  if (fCellMinE>0) {
    if (eev)
      am->LoadBranch("EMCALCells.");
    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = cells->GetCellNumber(iCell);
      Double_t cellE = cells->GetCellAmplitude(cellId);
      if (cellE>fCellMaxE)
        fCellMaxE = cellE;
    }
  }

  if (fClusMinE>0) {
    if (eev)
      am->LoadBranch("CaloClusters");

    const Int_t nCaloClusters = ev->GetNumberOfCaloClusters();
    for(Int_t iClus = 0; iClus<nCaloClusters; ++iClus) {
      AliVCluster *cl = ev->GetCaloCluster(iClus);
      if (!cl->IsEMCAL()) 
        continue;
      Double_t e = cl->E();
      if (e>fClusMaxE)
        fClusMaxE = e;
    }
  }

  if (fTrackMinPt>0) {
    TClonesArray *trks = 0;
    if (eev) {
      am->LoadBranch("PicoTracks");
      trks = dynamic_cast<TClonesArray*>(eev->FindListObject("PicoTracks"));
      if (!trks) {
        am->LoadBranch("Tracks");
        trks = dynamic_cast<TClonesArray*>(eev->FindListObject("Tracks"));
      }
    } else {
      trks = dynamic_cast<TClonesArray*>(aev->FindListObject("tracks"));
    }
    const Int_t Ntracks = trks->GetEntriesFast();
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVTrack *track = dynamic_cast<AliVTrack*>(trks->At(iTracks));
      if (!track)
        continue;
      Double_t pt = track->Pt();
      if (pt>fTrackMaxPt)
        fTrackMaxPt = pt;
    }
  }

  // bad cell criterion for LHC11a from 
  // https://indico.cern.ch/materialDisplay.py?contribId=4&materialId=slides&confId=147067
  const Int_t runN = ev->GetRunNumber();
  if ((runN>=144871) && (runN<=146860)) { 

    if (eev)
      am->LoadBranch("EMCALCells.");

    // count cells above threshold
    Int_t nCellCount[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
    for(Int_t iCell=0; iCell<nCells; ++iCell) {
      Short_t cellId = cells->GetCellNumber(iCell);
      Double_t cellE = cells->GetCellAmplitude(cellId);
      Int_t sm       = cellId / (24*48);
      if (cellE>0.1)
        ++nCellCount[sm];
    }

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

  if ((fClusMaxE>=fClusMinE) || (fTrackMaxPt>=fTrackMinPt))
    res |= kEmcalHT;

  if (fIsGoodEvent)
    res |= kEmcalOk;

  if ((fSkipLedEvent && fIsLedEvent) ||
      (fSkipFastOnly && fIsFastOnly))
    res = 0;

  return res;
}
