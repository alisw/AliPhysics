// ////////////////////////////////////////////////////////////////////////////
//
// CEP event buffer
//
// structure to hold event information
//
// Authors:
// P. Buehler, paul.buehler@oeaw.ac.at       27.06.2016
//  major revision                           27.01.2016
//
//
// ----------------------------------------------------------------------------
#include "CEPEventBuffer.h"

ClassImp(CEPEventBuffer)

// ----------------------------------------------------------------------------
CEPEventBuffer::CEPEventBuffer()
  : TObject()
  , fRunNumber(CEPTrackBuffer::kdumval)
  , fEventNumber(CEPTrackBuffer::kdumval)
  , fPeriodNumber(AliCEPBase::kdumval)
  , fOrbitNumber(AliCEPBase::kdumval)
  , fBunchCrossNumber(AliCEPBase::kdumval)
  , fCollissionType(AliCEPBase::kBinEventUnknown)
  , fMagnField(-999.9)
  , fPhysel(kFALSE)
  , fEvHandlersel(kFALSE)
  , fisPileup(kFALSE)
  , fisClusterCut(kFALSE)
  , fisMBOR(kFALSE)
  , fisMBAND(kFALSE)
  , fGapCondition(AliCEPBase::kBitBaseLine)
  , fCEPTracks(new TList())
  , fnTracklets(0)
  , fnTracks(0)
  , fnTracksCombined(0)
  , fnTracksITSpure(0)
  , fnResiduals(0)
  , fnMSelection(0)
  , fnV0(0)
  , fVtxType(-1)
  , fVtxPos(TVector3(-999.9,-999.9,-999.9))
  , fMCProcessType(AliCEPBase::kdumval)
  , fMCGenerator("")
  , fMCVtxPos(TVector3(-999.9,-999.9,-999.9))
{

		// the CEPTrackBuffer objects in fCEPTracks belong to fCEPTracks
    // they are normally created in an application and added to fCEPTracks
    // with the method CEPEventBuffer::AddTrack
    // they nedd however to be deleted with CEPEventBuffer::Reset()
    fCEPTracks->SetOwner(kTRUE);
  // printf("A CEPEventBuffer was created...\n");


}

// ----------------------------------------------------------------------------
CEPEventBuffer::~CEPEventBuffer()
{

  // clear the tracks and delete fCEPTracks
	if (fCEPTracks) {
		fCEPTracks->Clear();
		delete fCEPTracks;
		fCEPTracks = 0x0;
	}

}

// ----------------------------------------------------------------------------
void CEPEventBuffer::Reset()
{

  // reset all private variables
  // event header information
  fRunNumber       = AliCEPBase::kdumval;
  fEventNumber     = AliCEPBase::kdumval;
  fPeriodNumber    = AliCEPBase::kdumval;
  fOrbitNumber     = AliCEPBase::kdumval;
  fBunchCrossNumber= AliCEPBase::kdumval;
  fCollissionType  = AliCEPBase::kBinEventUnknown;
  fMagnField       = -999.9;

  // general event features
  fPhysel          = kFALSE;
  fEvHandlersel    = kFALSE;
  fisPileup        = kFALSE;
  fisClusterCut    = kFALSE;
  fisMBOR          = kFALSE;
  fisMBAND         = kFALSE;

  fGapCondition  = AliCEPBase::kBitBaseLine;  
  
  // summary track information
  fCEPTracks->Clear();
  fnTracklets     = 0;
  fnTracks        = 0;
  fnTracksCombined= 0;
  fnTracksITSpure = 0;
  fnResiduals     = 0;
  fnMSelection    = 0;
  fnV0            = 0;
  fVtxType        = -1;
  fVtxPos         = TVector3(-999.9,-999.9,-999.9);

  // Monte Carlo information
  fMCProcessType = AliCEPBase::kdumval;
  fMCGenerator   = "";
  fMCVtxPos      = TVector3(-999.9,-999.9,-999.9);
    
 }

// ----------------------------------------------------------------------------
void CEPEventBuffer::AddTrack(CEPTrackBuffer* trk)
{
  
  // add track to next element
  fCEPTracks->Add(trk);
  fnTracks++;
  
  // update track counters
  if (trk->TTisSet(AliCEPBase::kTTITSpure)) {
    fnTracksITSpure++;
  } else {
    fnTracksCombined++;
  }
  
}

// ----------------------------------------------------------------------------
CEPTrackBuffer* CEPEventBuffer::GetTrack(Int_t ind)
{
    
  // initialize the result track
  CEPTrackBuffer *trk = NULL;

  if (fCEPTracks->GetEntries() > ind) {
    trk = (CEPTrackBuffer*) fCEPTracks->At(ind);
  }
  
  return trk;
  
}

// ----------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracks(Int_t mask, Int_t pattern)
{
  
  // initialisations
  Int_t ngood = 0;
  
  for (Int_t ii=0; ii<GetnTracks(); ii++) {
    if ((GetTrack(ii)->GetTrackStatus() & mask) == pattern) ngood++;
  }
  
  return ngood;

}

// ----------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracks(Int_t mask, Int_t pattern, TArrayI *indices)
{
  
  // initialisations
  Int_t ngood = 0;
  indices->Reset(0);
  indices->Set(GetnTracks());
  
  for (Int_t ii=0; ii<GetnTracks(); ii++) {
    if ((GetTrack(ii)->GetTrackStatus() & mask) == pattern) {
    
      // increment the counter
      ngood++;
      
      // update the indices
      indices->SetAt(ii,ngood-1);
      
    }
  }
  indices->Set(ngood);
  
  return ngood;

}

// ----------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracks(TArrayI *masks, TArrayI *patterns)
{

  // checks status words with masks and patterns
  // makes logical OR of all tests
  
  // initialisations
  Bool_t ktmp;
  Int_t ngood = 0;

  // number of masks and patterns must agree
  Int_t ntests = masks->GetSize();
  if (ntests != patterns->GetSize()) return ngood;

  for (Int_t ii=0; ii<GetnTracks(); ii++) {
    
    ktmp = kFALSE;
    for (Int_t jj=0; jj<ntests;jj++) {
      
      if ((GetTrack(ii)->GetTrackStatus() & masks->At(jj)) ==
        patterns->At(jj)) {
        ktmp = kTRUE;
        break;
      }
    }
      
    if (ktmp) ngood++;

  }
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracks(TArrayI *masks, TArrayI *patterns,
  TArrayI* indices)
{

  // checks status words with masks and patterns
  // makes logical OR of all tests
  
  // initialisations
  Bool_t ktmp;
  Int_t ngood = 0;
  indices->Reset(0);

  // number of masks and patterns must agree
  Int_t ntests = masks->GetSize();
  if (ntests != patterns->GetSize()) return ngood;
  indices->Set(GetnTracks());

  for (Int_t ii=0; ii<GetnTracks(); ii++) {
    
    ktmp = kFALSE;
    for (Int_t jj=0; jj<ntests;jj++) {
      
      if ((GetTrack(ii)->GetTrackStatus() & masks->At(jj)) ==
        patterns->At(jj)) {
        ktmp = kTRUE;
        break;
      }
    }
      
    if (ktmp) {
      // increment the counter
      ngood++;
    
      // update the indices
      indices->SetAt(ii,ngood-1);
      
    }

  }
  indices->Set(ngood);
  
  return ngood;

}

// ------------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracksisSet(Int_t TTest)
{
  
  // determine number of tracks with a TrackStatus & TTest == TTest
  // all set bits in TTest are also set in TrackStatus
  Int_t nTrksSel = 0;
  
  // loop over all tracks contained in the TrackBuffer
  for (Int_t ii=0; ii<this->GetnTracks();ii++) {
    CEPTrackBuffer* trk = GetTrack(ii);
    if (!trk) continue;
      
    // increment counter
    if (trk->TTisSet(TTest)) nTrksSel++;

  }

  return nTrksSel;
  
}

// ----------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracksisEqual(Int_t TTest)
{
  
  // determine number of tracks with a TrackStatus == TTest
  Int_t nTrksSel = 0;
  
  // loop over all tracks contained in the TrackBuffer
  for (Int_t ii=0; ii<this->GetnTracks();ii++) {
    CEPTrackBuffer* trk = GetTrack(ii);
    if (!trk) continue;
      
    // increment counter
    if (trk->GetTrackStatus() == TTest) nTrksSel++;

  }

  return nTrksSel;
  
}

// ----------------------------------------------------------------------------


