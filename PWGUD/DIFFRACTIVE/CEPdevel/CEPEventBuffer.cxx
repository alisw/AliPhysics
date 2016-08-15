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
  , fRunNumber(0x0)
  , fEventNumber(0)
  , fnumTracks(0)
  , fnumSoftTracks(0)
  , fnumResiduals(0)
  , fTracks(0x0)
  , fGapCondition(0)
{

  fRunNumber = CEPTrackBuffer::kdumval;

}

// ----------------------------------------------------------------------------
CEPEventBuffer::~CEPEventBuffer()
{

  // delete list of tracks
  this->Reset();

}

// ----------------------------------------------------------------------------
void CEPEventBuffer::Reset()
{

  // reset all counters
  fRunNumber     = CEPTrackBuffer::kdumval;
  fEventNumber   = 0;
  fnumTracks     = 0;
  fnumSoftTracks = 0;
  fnumResiduals  = 0;
  fGapCondition  = 0;  

  // clear the track list
  //printf("I am reseting the ftracks TList!!\n");
  if (fTracks) fTracks->Clear();  

}

// ----------------------------------------------------------------------------
void CEPEventBuffer::AddTrack(CEPTrackBuffer* trk)
{
  
  // make sure that fTracks is defined
  if (!fTracks) fTracks = new TList();

  // add new track
  fTracks->Add((TObject*) trk);
  
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
    
  CEPTrackBuffer *trk = NULL;
  
  // get track at position ind - if there are at least ind tracks in the buffer
  if (fTracks) {
    if (fTracks->GetEntries() > ind) {
      trk = new CEPTrackBuffer();
      trk = (CEPTrackBuffer*) fTracks->At(ind);
    }
  }
  
  return trk;
  
}

// ----------------------------------------------------------------------------












