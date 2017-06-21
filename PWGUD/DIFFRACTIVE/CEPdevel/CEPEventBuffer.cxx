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
  , fMagnField(AliCEPBase::kdumval)
  , fFiredTriggerClasses(TString(""))
  , fPFBBFlagV0()
  , fPFBGFlagV0()
  , fPFBBFlagAD()
  , fPFBGFlagAD()
  , fEventCutsel(kFALSE)
  , fPhysel(kFALSE)
  , fisPileup(kFALSE)
  , fisClusterCut(kFALSE)
  , fisDGTrigger(kFALSE)
  , fEventCondition(AliCEPBase::kETBaseLine)
  , fnTracklets(0)
  , fnTracks(0)
  , fnTracksCombined(0)
  , fnTracksITSpure(0)
  , fnResiduals(0)
  , fnMSelection(0)
  , fnV0(0)
  , fVtxType(-1)
  , fVtxPos(TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval))
  , fMCProcessType(AliCEPBase::kdumval)
  , fMCGenerator("")
  , fMCVtxPos(TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval))
  , fCEPTracks(new TObjArray())
{

}

// ----------------------------------------------------------------------------
CEPEventBuffer::~CEPEventBuffer()
{

	// delete fCEPTracks and all the tracks it contains
  if (fCEPTracks) {
		fCEPTracks->SetOwner(kTRUE);
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
  fMagnField       = AliCEPBase::kdumval;
  fFiredTriggerClasses = TString("");

  // general event features
  fEventCutsel     = kFALSE;
  fPhysel          = kFALSE;
  fisPileup        = kFALSE;
  fisClusterCut    = kFALSE;
  fisDGTrigger     = kFALSE;

  fnTracklets     = 0;
  fnTracks        = 0;
  fnTracksCombined= 0;
  fnTracksITSpure = 0;
  fnResiduals     = 0;
  fnMSelection    = 0;
  fnV0            = 0;
  fVtxType        = -1;
  fVtxPos         = TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval);

  // Monte Carlo information
  fMCProcessType = AliCEPBase::kdumval;
  fMCGenerator   = "";
  fMCVtxPos      = TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval);
    
  // clear the track list
  fCEPTracks->SetOwner(kTRUE);
  fCEPTracks->Clear();
 }

// ----------------------------------------------------------------------------
void CEPEventBuffer::AddTrack(CEPTrackBuffer* trk)
{
  
  // add track to next element
  fCEPTracks->Add(trk);
  fnTracks++;
  
  // update track counters
  if (trk->GetTrackStatus() & AliCEPBase::kTTITSpure) {
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
Bool_t CEPEventBuffer::RemoveTrack(Int_t ind)
{
    
  // initialize the result track
  //printf("Removing track %i of %i ...",ind,fCEPTracks->GetEntries());
  Bool_t done = kFALSE;
  CEPTrackBuffer *trk = NULL;

  if (fCEPTracks->GetEntries() > ind) {
    trk = (CEPTrackBuffer*) fCEPTracks->RemoveAt(ind);
    fCEPTracks->Compress();
    //printf("ntracks left %i ...",fCEPTracks->GetEntries());

    // update track counters
    fnTracks--;
    if (trk->GetTrackStatus() & AliCEPBase::kTTITSpure) {
      fnTracksITSpure--;
    } else {
      fnTracksCombined--;
    }
    
    done = kTRUE;
  }
  //printf("... done %i!\n",done);
  
  return done;
  
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
void CEPEventBuffer::SetPFFlags(AliVEvent *Event)
{

  // V0
  // there are 64 channels - 32 channels per side (A, C) - and 21 bc values
  AliVVZERO* esdV0 = Event->GetVZEROData();
  if (esdV0) {
    // side A (32-63) and C (0-31)
    for (Int_t bc=0; bc<21; bc++) {
      fPFBBFlagV0[bc] = kFALSE;
      fPFBGFlagV0[bc] = kFALSE;
      for (Int_t ch=0; ch<64; ch++) {
        fPFBBFlagV0[bc] |= (esdV0->GetPFBBFlag(ch,bc)>0);
        fPFBGFlagV0[bc] |= (esdV0->GetPFBGFlag(ch,bc)>0);
      }
    }
  }
  
  // AD
  // there are 16 channels - 8 channels per side (A, C) - and 21 bc values
  AliESDAD* esdAD = (AliESDAD*)Event->GetADData();
  if (esdAD) {
    // side A (8-15) and C (0-7)
    for (Int_t bc=0; bc<21; bc++) {
      fPFBBFlagAD[bc] = kFALSE;
      fPFBGFlagAD[bc] = kFALSE;
      for (Int_t ch=0; ch<8; ch++) {
        fPFBBFlagAD[bc] |= (esdAD->GetPFBBFlag(ch,bc)>0);
        fPFBGFlagAD[bc] |= (esdAD->GetPFBGFlag(ch,bc)>0);
      }
    }
  }

}

// ------------------------------------------------------------------------------

