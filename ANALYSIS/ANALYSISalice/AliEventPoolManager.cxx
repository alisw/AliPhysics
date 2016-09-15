#include "AliEventPoolManager.h"
#include "TList.h"
#include "TRandom.h"
#include <iostream>

using std::cout;
using std::endl;
ClassImp(AliEventPool)

void AliEventPool::PrintInfo() const
{
  cout << Form("%20s: %d events", "Pool capacity", fMixDepth) << endl;
  cout << Form("%20s: %d events, %d tracks", "Current size", 
	       GetCurrentNEvents(), NTracksInPool()) << endl;
  cout << Form("%20s: %.1f to %.1f", "Sub-event mult.", fMultMin, fMultMax) << endl;
  cout << Form("%20s: %.1f to %.1f", "Z-vtx range", fZvtxMin, fZvtxMax) << endl;
  cout << Form("%20s: %.1f to %.1f", "Psi range", fPsiMin, fPsiMax) << endl;
  cout << Form("%20s: %.1f to %.1f", "Pt range", fPtMin, fPtMax) << endl;

  return;
}

Bool_t AliEventPool::EventMatchesBin(Int_t mult, Double_t zvtx, Double_t psi, Double_t pt) const
{
  return EventMatchesBin((Double_t) mult, zvtx, psi, pt);
}

Bool_t AliEventPool::EventMatchesBin(Double_t mult, Double_t zvtx, Double_t psi, Double_t pt) const
{
  // Lower bin limit included; upper limit excluded.

  Bool_t multOK = (mult >= fMultMin && mult < fMultMax);
  Bool_t zvtxOK = (zvtx >= fZvtxMin && zvtx < fZvtxMax);
  Bool_t psiOK  = (psi >= fPsiMin   && psi  < fPsiMax);
  Bool_t ptOK   = (pt >= fPtMin     && pt   < fPtMax);

  return (multOK && zvtxOK && psiOK && ptOK);
}

Int_t AliEventPool::NTracksInPool() const
{
  // Number of tracks for this cent, zvtx bin; possibly includes many
  // events.

  Int_t ntrk=0;
  for (Int_t i=0; i<(Int_t)fEvents.size(); ++i) {
    ntrk += fNTracksInEvent.at(i);
  }
  return ntrk;
}

Int_t AliEventPool::SetEventMultRange(Int_t multMin, Int_t multMax)
{
  fMultMin = (Double_t)multMin;
  fMultMax = (Double_t)multMax;
  return 0;
}

Int_t AliEventPool::SetEventMultRange(Double_t multMin, Double_t multMax)
{
  fMultMin = multMin;
  fMultMax = multMax;
  return 0;
}

Int_t AliEventPool::SetEventZvtxRange(Double_t zvtxMin, Double_t zvtxMax)
{
  fZvtxMin = zvtxMin;
  fZvtxMax = zvtxMax;
  return 0;
}

Int_t AliEventPool::SetEventPsiRange(Double_t psiMin, Double_t psiMax)
{
  fPsiMin = psiMin;
  fPsiMax = psiMax;
  return 0;
}

Int_t AliEventPool::SetEventPtRange(Double_t ptMin, Double_t ptMax)
{
  fPtMin = ptMin;
  fPtMax = ptMax;
  return 0;
}

Int_t AliEventPool::GlobalEventIndex(Int_t j) const
{
  // Index returned from passing local pool event index.

  if (j < 0 || j >= (Int_t)fEventIndex.size()) {
    cout << "ERROR in AliEventPool::GlobalEventIndex(): "
	 << " Invalid index " << j << endl;
    return -99;
  }
  return fEventIndex.at(j);
}

Int_t AliEventPool::UpdatePool(TObjArray *trk)
{
  // A rolling buffer (a double-ended queue) is updated by removing
  // the oldest event, and appending the newest.
  //
  // the ownership of <trk> is delegated to this class

  if(fLockFlag)
  {
    AliFatal("Tried to fill a locked AliEventPool.");
    return fEvents.size();
  }

  static Int_t iEvent = -1; 
  iEvent++;

  Int_t mult = trk->GetEntries();
  
  // Do not fill empty events
  if (mult == 0)
    return fEvents.size();
  
  Int_t nTrk = NTracksInPool();

  if (!IsReady() && IsReady(nTrk + mult, GetCurrentNEvents() + 1))
    fNTimes++;

  // remove 0th element before appending this event
  Bool_t removeFirstEvent = 0;
  if (nTrk>fTargetTrackDepth) {
    Int_t nTrksFirstEvent= fNTracksInEvent.front();
    Int_t diff = nTrk - nTrksFirstEvent + mult;
    if (diff>fTargetTrackDepth)
      removeFirstEvent = 1;
  }
  if (removeFirstEvent) {
    TObjArray *fa = fEvents.front();
    delete fa;
    fEvents.pop_front();         // remove first track array 
    fNTracksInEvent.pop_front(); // remove first int
    fEventIndex.pop_front();
  }

  fNTracksInEvent.push_back(mult);
  fEvents.push_back(trk);
  fEventIndex.push_back(iEvent);

  if (fNTimes==1) {
    fFirstFilled = kTRUE;
    if (AliEventPool::fDebug) {
      cout << "\nPool " << MultBinIndex() << ", " << ZvtxBinIndex() 
           << " ready at event "<< iEvent;
      PrintInfo();
      cout << endl;
    }
    fNTimes++; // See this message exactly once/pool
  } else {
    fFirstFilled = kFALSE;
  }

  fWasUpdated = true;

  if (AliEventPool::fDebug) {
    cout << " Event " << fEventIndex.back() << endl;
    cout << " PoolDepth = " << GetCurrentNEvents() << endl;
    cout << " NTracksInCurrentEvent = " << NTracksInCurrentEvent() << endl;
  }

  return fEvents.size();
}

Long64_t AliEventPool::Merge(TCollection* hlist)
{
  if (!hlist)
  	return 0;

  Bool_t origLock = fLockFlag;
  fLockFlag = kFALSE; // temporary deactivate lockflag to allow filling
  AliEventPool* tmpObj = 0;
  TIter objIter(hlist);
  // Iterate through all objects to be merged
  while ( (tmpObj = static_cast<AliEventPool*>(objIter())) )
  {
    // Update this pool (it won't get fuller than demanded)
    for(Int_t i=0; i<tmpObj->fEvents.size(); i++)
      UpdatePool(tmpObj->fEvents.at(i));
  }
  fLockFlag = origLock;
  return hlist->GetEntries() + 1;
}

void AliEventPool::Clear()
{
  // Clear the pool without deleting the object
  // Don't touch lock or save flag here to be fully flexible
  fEvents.clear();
  fNTracksInEvent.clear();
  fEventIndex.clear();
  fWasUpdated = 0;
  fFirstFilled = 0;
  fWasUpdated = 0;
  fFirstFilled = 0;
  fNTimes = 0;
}

TObject* AliEventPool::GetRandomTrack() const
{
  // Get any random track from the pool, sampled with uniform probability.

  UInt_t ranEvt = gRandom->Integer(fEvents.size());
  TObjArray *tca = fEvents.at(ranEvt);
  UInt_t ranTrk = gRandom->Integer(tca->GetEntries());
  TObject *trk = (TObject*)tca->At(ranTrk);
  return trk;
}

TObjArray* AliEventPool::GetEvent(Int_t i) const
{
  if (i<0 || i>=(Int_t)fEvents.size()) {
    cout << "AliEventPool::GetEvent(" 
	 << i << "): Invalid index" << endl;
    return 0x0;
  }

  TObjArray *tca = fEvents.at(i);
  return tca;
}

TObjArray* AliEventPool::GetRandomEvent() const
{
  UInt_t ranEvt = gRandom->Integer(fEvents.size());
  TObjArray *tca = fEvents.at(ranEvt);
  return tca;
}

Int_t AliEventPool::NTracksInEvent(Int_t iEvent) const
{
  // Return number of tracks in iEvent, which is the local pool index.

  Int_t n = -1;
  Int_t curEvent = fEventIndex.back();
  Int_t offset = curEvent - iEvent;
  Int_t pos = fEventIndex.size() - offset - 1;

  if (offset==0)
    n = fNTracksInEvent.back();
  else if (offset < 0 || iEvent < 0) {
    n = 0;
  }
  else if (offset > 0 && offset <= (int)fEventIndex.size()) {
    n = fNTracksInEvent.at(pos);
  }
  else
    cout << "Event info no longer in memory" << endl;
  return n;
}

ClassImp(AliEventPoolManager)

AliEventPoolManager::AliEventPoolManager(Int_t depth,     Int_t minNTracks,
					 Int_t nMultBins, Double_t *multbins,
					 Int_t nZvtxBins, Double_t *zvtxbins) :
fDebug(0), fNMultBins(0), fNZvtxBins(0), fNPsiBins(0), fNPtBins(0), fMultBins(), fZvtxBins(), fPsiBins(), fPtBins(), fEvPool(0), fTargetTrackDepth(minNTracks) 
{
  // Constructor.
  // without Event plane bins or pt bins
  Int_t nPsiBins = 1;
  Double_t psibins[2] = {-999.,999.};
  Int_t nPtBins = 1;
  Double_t ptbins[2] = {-9999.,9999.};

  InitEventPools(depth, nMultBins, multbins, nZvtxBins, zvtxbins, nPsiBins, psibins, nPtBins, ptbins);
  cout << "AliEventPoolManager initialized." << endl;
}

AliEventPoolManager::AliEventPoolManager(Int_t depth,     Int_t minNTracks,
					 Int_t nMultBins, Double_t *multbins,
					 Int_t nZvtxBins, Double_t *zvtxbins,
					 Int_t nPsiBins, Double_t *psibins) :
fDebug(0), fNMultBins(0), fNZvtxBins(0), fNPsiBins(0), fNPtBins(0), fMultBins(), fZvtxBins(), fPsiBins(), fPtBins(), fEvPool(0), fTargetTrackDepth(minNTracks) 
{
  // Constructor.
  // without pt bins
  Int_t nPtBins = 1;
  Double_t ptbins[2] = {-9999.,9999.};

  InitEventPools(depth, nMultBins, multbins, nZvtxBins, zvtxbins, nPsiBins, psibins, nPtBins, ptbins);
  cout << "AliEventPoolManager initialized." << endl;
}

AliEventPoolManager::AliEventPoolManager(Int_t depth,     Int_t minNTracks,
					 Int_t nMultBins, Double_t *multbins,
					 Int_t nZvtxBins, Double_t *zvtxbins,
					 Int_t nPsiBins, Double_t *psibins,
           Int_t nPtBins, Double_t *ptbins) :
fDebug(0), fNMultBins(0), fNZvtxBins(0), fNPsiBins(0), fNPtBins(0), fMultBins(), fZvtxBins(), fPsiBins(), fPtBins(), fEvPool(0), fTargetTrackDepth(minNTracks) 
{
  // Constructor.

  InitEventPools(depth, nMultBins, multbins, nZvtxBins, zvtxbins, nPsiBins, psibins, nPtBins, ptbins);
  cout << "AliEventPoolManager initialized." << endl;
}

AliEventPoolManager::AliEventPoolManager(Int_t depth,     Int_t minNTracks, const char* binning) :
fDebug(0), fNMultBins(0), fNZvtxBins(0), fNPsiBins(0), fNPtBins(0), fMultBins(), fZvtxBins(), fPsiBins(), fPtBins(), fEvPool(0), fTargetTrackDepth(minNTracks) 
{
  Double_t psidummy[2] = {-999.,999.};
  Double_t ptdummy[2] = {-9999.,9999.};

  // Constructor.
  Double_t* multbins = GetBinning(binning, "multiplicity", fNMultBins);
  Double_t* zvtxbins = GetBinning(binning, "vertex", fNZvtxBins);
  
  Double_t* psibins = GetBinning(binning, "psi", fNPsiBins); //optional
  if (!psibins) {
    psibins = psidummy;
    fNPsiBins = 1;
  }
  
  Double_t* ptbins = GetBinning(binning, "pt", fNPtBins); //optional
  if (!ptbins) {
    ptbins = ptdummy;
    fNPtBins = 1;
  }

  InitEventPools(depth, fNMultBins, multbins, fNZvtxBins, zvtxbins, fNPsiBins, psibins, fNPtBins, ptbins);

  cout << "AliEventPoolManager initialized." << endl;
}

Int_t AliEventPoolManager::InitEventPools(Int_t depth, 
					  Int_t nMultBins, Double_t *multbin, 
					  Int_t nZvtxBins, Double_t *zvtxbin, 
					  Int_t nPsiBins, Double_t *psibin,
            Int_t nPtBins, Double_t *ptbin)
{
  // Assign AliEventPoolManager members. (with Event plane + pt)

  fNMultBins = nMultBins;
  fNZvtxBins = nZvtxBins;
  fNPsiBins  = nPsiBins;
  fNPtBins   = nPtBins;

  fMultBins.assign(multbin, multbin+nMultBins+1);
  fZvtxBins.assign(zvtxbin, zvtxbin+nZvtxBins+1);
  fPsiBins.assign(psibin, psibin+nPsiBins+1);
  fPtBins.assign(ptbin, ptbin+nPtBins+1);
  
  for (Int_t iM=0; iM<nMultBins; iM++) {
    for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
      for (Int_t iP=0; iP<nPsiBins; iP++) {
        for (Int_t iPt=0; iPt<nPtBins; iPt++) {

          fEvPool.push_back(new AliEventPool(depth, 
             multbin[iM], multbin[iM+1], 
             zvtxbin[iZ], zvtxbin[iZ+1],
             psibin[iP], psibin[iP+1],
             ptbin[iPt], ptbin[iPt+1] ));
        }
      }
    }
  }
  
  
  for (Int_t iM=0; iM<nMultBins; iM++) {
    for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
      for (Int_t iP=0; iP<nPsiBins; iP++) {
        for (Int_t iPt=0; iPt<nPtBins; iPt++) {
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetMultBinIndex(iM);
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetZvtxBinIndex(iZ);
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetPsiBinIndex(iP);
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetPtBinIndex(iPt);
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetTargetTrackDepth(fTargetTrackDepth);
        }
      }
    }
  }
    
  if (0) {
    cout << "fEvPool outer size: " << fEvPool.size() << endl;
    for (Int_t iM=0; iM<nMultBins; iM++) {
      for (Int_t iZ=0; iZ<nZvtxBins; iZ++) {
        for (Int_t iP=0; iP<nPsiBins; iP++) {
          for (Int_t iPt=0; iPt<nPtBins; iPt++) {
            if(fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))) {
              cout << "multiplicity bin: " << iM;
              cout << ", z-vertex bin: " << iZ;
              cout << ", psi bin: " << iP;
              cout << ", pt bin: " << iPt;
              fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->PrintInfo();
            }
          }
        }
      }
    }
  }
  
  return fEvPool.size();
}

Long64_t AliEventPoolManager::Merge(TCollection* hlist)
{
  if (!hlist)
  	return 0;
  	
  AliEventPoolManager* tmpObj = 0;
  TIter objIter(hlist);

  // Iterate through all objects to be merged
  while ( (tmpObj = static_cast<AliEventPoolManager*>(objIter())) )
  {
    for(Int_t i = 0; i<GetNumberOfMultBins(); i++)
      for(Int_t j = 0; j<GetNumberOfZVtxBins(); j++)
        for(Int_t k = 0; k<GetNumberOfPsiBins(); k++)
          for(Int_t l = 0; l<GetNumberOfPtBins(); l++)
          {
            TList* poolList = new TList();
            AliEventPool* objPool = tmpObj->GetEventPool(i,j,k,l);
            AliEventPool* pool    = GetEventPool(i,j,k,l);

            poolList->Add(objPool);
            pool->Merge(poolList);
            delete poolList;
          }
  }
  return hlist->GetEntries() + 1;
}

void AliEventPoolManager::SetTargetValues(Int_t trackDepth, Float_t fraction, Int_t events)
{
  // sets target values (when a pool becomes ready) in all event pools
  
  fTargetTrackDepth = trackDepth;
  
  for (Int_t iM=0; iM<fNMultBins; iM++) {
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++) {
      for (Int_t iP=0; iP<fNPsiBins; iP++) {
        for (Int_t iPt=0; iPt<fNPtBins; iPt++) {
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetTargetTrackDepth(trackDepth, fraction);
          fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->SetTargetEvents(events);
        }
      }
    }
  }
}

void AliEventPoolManager::ClearPools()
{
  // Clear the pools that are not marked to be saved
  // Those marked to be saved are now flagged as locked and save deflagged
  // to serve a valid input when importing this pool
  // Call this function in FinishTaskOutput() of your class

  for(Int_t i = 0; i<GetNumberOfMultBins(); i++)
    for(Int_t j = 0; j<GetNumberOfZVtxBins(); j++)
      for(Int_t k = 0; k<GetNumberOfPsiBins(); k++)
        for(Int_t l = 0; l<GetNumberOfPtBins(); l++)
        {
          AliEventPool* pool    = GetEventPool(i,j,k,l);
          if(!pool->GetSaveFlag())
            pool->Clear();
          else
          {
            pool->SetLockFlag(kTRUE);
            pool->SetSaveFlag(kFALSE);
          }
        }
}

void AliEventPoolManager::ClearPools(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi, Double_t maxPsi, Double_t minPt, Double_t maxPt)
{
  // Clear some pools, given by the ranges
  for(Int_t i = 0; i<GetNumberOfMultBins(); i++)
    for(Int_t j = 0; j<GetNumberOfZVtxBins(); j++)
      for(Int_t k = 0; k<GetNumberOfPsiBins(); k++)
        for(Int_t l = 0; l<GetNumberOfPtBins(); l++)
        {
          AliEventPool* pool    = GetEventPool(i,j,0,l);
          if( (minCent < pool->GetMultMax()) && (maxCent > pool->GetMultMin()) &&
              (minZvtx < pool->GetZvtxMax()) && (maxZvtx > pool->GetZvtxMin()) &&
              (minPsi  < pool->GetPsiMax())  && (maxPsi  > pool->GetPsiMin()) &&
              (minPt   < pool->GetPtMax())   && (maxPt   > pool->GetPtMin()) )
          {
            pool->SetLockFlag(kFALSE); // unlock pool
            pool->Clear(); //clear pool
          }
        }
}

void AliEventPoolManager::SetSaveFlag(Double_t minCent, Double_t maxCent, Double_t minZvtx, Double_t maxZvtx, Double_t minPsi, Double_t maxPsi, Double_t minPt, Double_t maxPt)
{
  // set save flag on the pools in range
  for(Int_t i = 0; i<GetNumberOfMultBins(); i++)
    for(Int_t j = 0; j<GetNumberOfZVtxBins(); j++)
      for(Int_t k = 0; k<GetNumberOfPsiBins(); k++)
        for(Int_t l = 0; l<GetNumberOfPtBins(); l++)
        {
          AliEventPool* pool    = GetEventPool(i,j,k,l);
          if( (minCent < pool->GetMultMax()) && (maxCent > pool->GetMultMin()) &&
              (minZvtx < pool->GetZvtxMax()) && (maxZvtx > pool->GetZvtxMin()) &&
              (minPsi  < pool->GetPsiMax())  && (maxPsi  > pool->GetPsiMin()) &&
              (minPt   < pool->GetPtMax())   && (maxPt   > pool->GetPtMin()) )
          {
            if(pool->GetLockFlag())
              AliWarning("A pool that was already imported is flagged to be saved. Is this really intended?");
            pool->SetSaveFlag(kTRUE);
          }
        }  
}

void AliEventPoolManager::Validate()
{
  std::cout << "############## AliEventPoolManager ##############\n"; 
  std::cout << "== Binning ==\n"; 
  TString tmpStr = "cent/mult: ";
  for(Int_t i = 0; i<fMultBins.size(); i++)
    tmpStr += Form("%4.3f ", fMultBins[i]);
  tmpStr += " vertex: ";
  for(Int_t i = 0; i<fZvtxBins.size(); i++)
    tmpStr += Form("%4.3f ", fZvtxBins[i]);
  tmpStr += " psi: ";
  for(Int_t i = 0; i<fPsiBins.size(); i++)
    tmpStr += Form("%4.3f ", fPsiBins[i]);
  tmpStr += " pt: ";
  for(Int_t i = 0; i<fPtBins.size(); i++)
    tmpStr += Form("%4.3f ", fPtBins[i]);
  std::cout << tmpStr.Data() << std::endl;
  std::cout << "== Event pools ==\n";

  std::cout << "Note: Locked pools won't be filled. They are intended to serve as external input.\n";
  std::cout << "      Pools with save flag: Those pools are intended to be written to the output file.\n";

  for (Int_t iM=0; iM<fNMultBins; iM++)
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++)
      for (Int_t iP=0; iP<fNPsiBins; iP++)
        for (Int_t iPt=0; iPt<fNPtBins; iPt++) 
        {
          AliEventPool* pool = GetEventPool(iM, iZ, iP, iPt);
          if(!pool)
            AliFatal(Form("Pool (%i,%i,%i,%i) is not correctly initialized!", iM, iZ, iP, iPt));
          if(pool->GetLockFlag())
            std::cout << Form(" Pool (mult=%4.3f-%4.3f, zvertex=%4.3f-%4.3f, psi=%4.3f-%4.3f, pt=%4.3f-%4.3f) locked", fMultBins[iM], fMultBins[iM+1], fZvtxBins[iZ], fZvtxBins[iZ+1], fPsiBins[iP], fPsiBins[iP+1], fPtBins[iPt], fPtBins[iPt+1]) << std::endl;
        }


  for (Int_t iM=0; iM<fNMultBins; iM++)
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++)
      for (Int_t iP=0; iP<fNPsiBins; iP++)
        for (Int_t iPt=0; iPt<fNPtBins; iPt++) 
        {
          AliEventPool* pool = GetEventPool(iM, iZ, iP, iPt);
          if(pool->GetSaveFlag())
            std::cout << Form(" Pool (mult=%4.3f-%4.3f, zvertex=%3.3f-%4.3f, psi=%4.3f-%4.3f, pt=%4.3f-%4.3f) will be saved", fMultBins[iM], fMultBins[iM+1], fZvtxBins[iZ], fZvtxBins[iZ+1], fPsiBins[iP], fPsiBins[iP+1], fPtBins[iPt], fPtBins[iPt+1]) << std::endl;
        }

  std::cout << "############## AliEventPoolManager ##############\n"; 

}


AliEventPool *AliEventPoolManager::GetEventPool(Int_t iMult, Int_t iZvtx, Int_t iPsi, Int_t iPt) const
{
  if (iMult < 0 || iMult >= fNMultBins) 
  {
    AliError(Form("Mult bin %i exceeds maximum of %i",iMult, fNMultBins));
    return 0x0;
  }
  if (iZvtx < 0 || iZvtx >= fNZvtxBins) 
  {
    AliError(Form("Zvtx bin %i exceeds maximum of %i",iZvtx, fNZvtxBins));
    return 0x0;
  }
  if (iPsi < 0 || iPsi >= fNPsiBins) 
  {
    AliError(Form("Psi bin %i exceeds maximum of %i",iPsi, fNPsiBins));
    return 0x0;
  }
  if (iPt < 0 || iPt >= fNPtBins) 
  {
    AliError(Form("Pt bin %i exceeds maximum of %i",iPt, fNPtBins));
    return 0x0;
  }

  if(fEvPool.at(GetBinIndex(iMult, iZvtx, iPsi, iPt)))
    return fEvPool.at(GetBinIndex(iMult, iZvtx, iPsi, iPt));
  else
    AliFatal("Pool does not exist");

  return 0;
}

AliEventPool *AliEventPoolManager::GetEventPool(Int_t centVal, Double_t zVtxVal, Double_t psiVal, Int_t iPt) const
{
  return GetEventPool((Double_t)centVal, zVtxVal, psiVal, iPt);
}

AliEventPool *AliEventPoolManager::GetEventPool(Double_t centVal, Double_t zVtxVal, Double_t psiVal, Int_t iPt) const
{
  // Return appropriate pool for this centrality and z-vertex value.

  for (Int_t iM=0; iM<fNMultBins; iM++) {
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++) {
      for (Int_t iP=0; iP<fNPsiBins; iP++) {
        AliEventPool* pool = GetEventPool(iM, iZ, iP, iPt);
        if (pool->EventMatchesBin(centVal, zVtxVal, psiVal, (pool->GetPtMin()+pool->GetPtMax())/2 ))
          return pool;
      }
    }
  }
  return 0x0;
}

Int_t AliEventPoolManager::UpdatePools(TObjArray *trk)
{
  // Call UpdatePool for all bins.

  for (Int_t iM=0; iM<fNMultBins; iM++) {
    for (Int_t iZ=0; iZ<fNZvtxBins; iZ++) {
      for (Int_t iP=0; iP<fNPsiBins; iP++) {
        for (Int_t iPt=0; iPt<fNPtBins; iPt++) {
          if (fEvPool.at(GetBinIndex(iM, iZ, iP, iPt))->UpdatePool(trk) > -1)
            break;
        }
      }
    }
  }  
  return 0;
}

Double_t* AliEventPoolManager::GetBinning(const char* configuration, const char* tag, Int_t& nBins) const
{
  // Same as in AliUEHist
  
  TString config(configuration);
  TObjArray* lines = config.Tokenize("\n");
  for (Int_t i=0; i<lines->GetEntriesFast(); i++)
  {
    TString line(lines->At(i)->GetName());
    if (line.BeginsWith(TString(tag) + ":"))
    {
      line.Remove(0, strlen(tag) + 1);
      line.ReplaceAll(" ", "");
      TObjArray* binning = line.Tokenize(",");
      Double_t* bins = new Double_t[binning->GetEntriesFast()];
      for (Int_t j=0; j<binning->GetEntriesFast(); j++)
        bins[j] = TString(binning->At(j)->GetName()).Atof();
      
      nBins = binning->GetEntriesFast() - 1;

      delete binning;
      delete lines;
      return bins;
    }
  }
  
  delete lines;
  return 0;
}

