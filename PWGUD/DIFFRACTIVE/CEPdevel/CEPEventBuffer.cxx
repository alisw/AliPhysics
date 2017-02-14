// ////////////////////////////////////////////////////////////////////////////
//
// CEP event buffer
//
// structure to hold event information
//
// Authors:
// P. Buehler, paul.buehler@oeaw.ac.at       27.06.2016
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
  , fFiredTriggerClasses(TString(""))
  , fMagnField(-999.9)
  , fPhysel(kFALSE)
  , fEventCutsel(kFALSE)
  , fisPileup(kFALSE)
  , fisClusterCut(kFALSE)
  , fisDGTrigger(kFALSE)
  , fGapCondition(AliCEPBase::kBitBaseLine)
  , fCEPTracks(new TObjArray())
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
  , fVertexPos(TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval))
  , fMCVertexPos(TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval))
{

  // printf("A CEPEventBuffer was created...\n");

}

// ----------------------------------------------------------------------------
CEPEventBuffer::~CEPEventBuffer()
{

  //this->Reset();
  // printf("A CEPEventBuffer was reset...\n");
  //for (Int_t ii=0; ii<ftrb.GetEntries(); ii++)  {
  //  ftrb[ii]->Delete();
  //}
	if (fCEPTracks) {
		fCEPTracks->SetOwner(kTRUE);
		fCEPTracks->Clear();
		delete fCEPTracks;
		fCEPTracks = 0x0;
	}

  //fCEPTracks->Delete();
  //delete fCEPTracks;
  //fCEPTracks = 0x0;
  
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
  fFiredTriggerClasses = TString("");

  // general event features
  fEventCutsel     = kFALSE;
  fPhysel          = kFALSE;
  fisPileup        = kFALSE;
  fisClusterCut    = kFALSE;
  fisDGTrigger     = kFALSE;

  // clear the track list
  fCEPTracks->Clear();
  fnTracks        = 0;
  fnTracksITSpure = 0;
  fnTracksCombined= 0;
  fnTracklets     = 0;
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
  // printf("<I - CEPEventBuffer::AddTrack> %i\n",fnTracks);
  
  // update track counter
  if (trk->GetisSoft()) {
    fnumSoftTracks++;
  } else {
    fnumTracks++;
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
Int_t CEPEventBuffer::GetnTracks(UInt_t mask, UInt_t pattern)
{
  
  // initialisations
  Int_t ngood = 0;
  
  for (Int_t ii=0; ii<GetnTracks(); ii++) {
    if ((GetTrack(ii)->GetTrackStatus() & mask) == pattern) ngood++;
  }
  
  return ngood;

}

// ----------------------------------------------------------------------------
Int_t CEPEventBuffer::GetnTracks(UInt_t mask, UInt_t pattern, TArrayI *indices)
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
Int_t CEPEventBuffer::GetnTracksisSet(UInt_t TTest)
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
Int_t CEPEventBuffer::GetnTracksisEqual(UInt_t TTest)
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
