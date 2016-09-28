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
  , fEventNumber(0)
  , fnumTracks(0)
  , fnumSoftTracks(0)
  , fnumResiduals(0)
  , fMCProcessType(-1)
  , fGapCondition(0)
  , fCEPTracks(new TClonesArray("CEPTrackBuffer"))
  , ftrb(*fCEPTracks)
{

}

// ----------------------------------------------------------------------------
CEPEventBuffer::~CEPEventBuffer()
{

  this->Reset();
  
}

// ----------------------------------------------------------------------------
void CEPEventBuffer::Reset()
{

  // reset all counters
  fRunNumber     = CEPTrackBuffer::kdumval;
  fEventNumber   = 0;
  fnumResiduals  = 0;
  fGapCondition  = 0;  
  fMCProcessType = -1;

  // clear the track list
  ftrb.Clear();
  fnumTracks     = 0;
  fnumSoftTracks = 0;

}

// ----------------------------------------------------------------------------
void CEPEventBuffer::AddTrack(CEPTrackBuffer* trk)
{
  
  // add track to next element
  ftrb[ftrb.GetEntries()] = trk;
  
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
