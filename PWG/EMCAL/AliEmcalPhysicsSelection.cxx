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
  fZvertex(-1),
  fZvertexDiff(0),
  fCentMin(-1),
  fCentMax(-1),
  fMinCellTrackScale(-1),
  fMaxCellTrackScale(-1),
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
  if (!ev)
    return 0;

  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();

  UInt_t res = 0;
  const AliESDEvent *eev = dynamic_cast<const AliESDEvent*>(obj);
  const AliAODEvent *aev = 0;
  if (eev) {
    am->LoadBranch("AliESDHeader.");
    am->LoadBranch("AliESDRun.");
    TString title(eev->GetHeader()->GetTitle());
    if (1&&(title.Length()>0)) {
      res = ((AliVAODHeader*)eev->GetHeader())->GetUniqueID();
      res &= 0x4FFFFFFF;
    } else {
      res = IsCollisionCandidate(eev); 
    }
  } else {
    aev = dynamic_cast<const AliAODEvent*>(obj);
    res = ((AliVAODHeader*)aev->GetHeader())->GetOfflineTrigger();
  }

  // return 0, if 0 found
  if (res==0)
    return 0;

  // return 0, if ptrs are not set
  if ((eev==0) && (aev==0))
    return 0;

  if (fTriggers) { // only process given triggers
    if ((res & fTriggers) == 0)
      return res;
  }

  fIsFastOnly  = kFALSE;
  fIsGoodEvent = kFALSE;
  fIsLedEvent  = kFALSE;
  fCellMaxE    = -1;
  fClusMaxE    = -1;
  fTrackMaxPt  = -1;

  if ((res & AliVEvent::kAnyINT)      || 
      (res & AliVEvent::kSemiCentral) || 
      (res & AliVEvent::kCentral)     || 
      (res & AliVEvent::kEMC1)        || 
      (res & AliVEvent::kEMC7)        || 
      (res & AliVEvent::kEMCEJE)      || 
      (res & AliVEvent::kEMCEGA))
    fIsGoodEvent = kTRUE;
  else {
    return 0;
  }

  if (fZvertexDiff || (fZvertex>0)) {
    Double_t vzPRI = +999;
    Double_t vzSPD = -999;
    const AliVVertex *pv = 0;
    if (eev) {
      am->LoadBranch("PrimaryVertex.");
      pv = eev->GetPrimaryVertex();
    } else {
      pv = aev->GetPrimaryVertex();
    }
    if (pv && pv->GetNContributors()>0) {
      vzPRI = pv->GetZ();
    }
    const AliVVertex *sv = 0;
    if (eev) {
      am->LoadBranch("SPDVertex.");
      sv = eev->GetPrimaryVertexSPD();
    } else {
      sv = aev->GetPrimaryVertexSPD();
    }
    if (sv && sv->GetNContributors()>0) {
      vzSPD = sv->GetZ();
    }
    Double_t  dvertex = TMath::Abs(vzPRI-vzSPD);
    // skip events with dvertex>1mm if requested
    // https://indico.cern.ch/getFile.py/access?contribId=4&resId=0&materialId=slides&confId=189624
    // also check on vertex z if requested
    if (fZvertexDiff && (dvertex>0.1))
      fIsGoodEvent = kFALSE;
    if ((fZvertex>0) && (TMath::Abs(vzPRI)>fZvertex))
      fIsGoodEvent = kFALSE;
  }

  if ((fCentMin>-1) && (fCentMax>-1)) {
    Double_t v0mcent = -1;
    AliCentrality *centin = 0;
    if (eev) {
      TTree *tree = am->GetTree();
      if (tree->GetBranch("Centrality."))
        am->LoadBranch("Centrality.");
      centin = dynamic_cast<AliCentrality*>(eev->FindListObject("Centrality"));
    } else {
      centin = const_cast<AliAODEvent*>(aev)->GetCentrality();
    }
    if (centin)
      v0mcent = centin->GetCentralityPercentileUnchecked("V0M");
    if ((v0mcent<fCentMin) || (v0mcent>fCentMax))
      fIsGoodEvent = kFALSE;
  }

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

  TClonesArray *trks = 0;
  Int_t Ntracks = 0;

  if (fTrackMinPt>0 || fMinCellTrackScale > 0 || fMaxCellTrackScale > 0) {
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
    if (trks)
      Ntracks = trks->GetEntriesFast();
  }
  
  Int_t nAccTracks = 0;

  if (fTrackMinPt>0 || fMinCellTrackScale > 0 || fMaxCellTrackScale > 0) {
    for (Int_t iTracks = 0; iTracks < Ntracks; ++iTracks) {
      AliVTrack *track = static_cast<AliVTrack*>(trks->At(iTracks));
      if (!track)
	continue;
      if (aev) { // a bit ugly since cuts are hard coded for now
	AliAODTrack *aodtrack = static_cast<AliAODTrack*>(track);
	if (!aodtrack->TestFilterBit(256) && !aodtrack->TestFilterBit(512))
	  continue;
      }
      nAccTracks++;
      Double_t pt = track->Pt();
      if (pt>fTrackMaxPt)
	fTrackMaxPt = pt;
    }
  }

  if (fMinCellTrackScale > 0 || fMaxCellTrackScale > 0) {
    if (nCells < fMinCellTrackScale * nAccTracks || nCells > fMaxCellTrackScale * nAccTracks)
      fIsGoodEvent = kFALSE;
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

  if (fIsGoodEvent) {
    if (fCellMaxE>fCellMinE)
      res |= kEmcalHC;
    if ((fClusMaxE>fClusMinE) || (fTrackMaxPt>fTrackMinPt))
      res |= kEmcalHT;
    res |= kEmcalOk;
  }

  if ((fSkipLedEvent && fIsLedEvent) ||
      (fSkipFastOnly && fIsFastOnly))
    res = 0;

  return res;
}
