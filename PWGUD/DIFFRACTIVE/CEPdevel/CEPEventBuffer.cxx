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
  , fisPileup(kFALSE)
  , fnumTracks(0)
  , fnumSoftTracks(0)
  , fnumResiduals(0)
  , fGapCondition(0)
  , fCEPTracks(new TClonesArray("CEPTrackBuffer"))
  , ftrb(*fCEPTracks)
  , fMCProcessType(-1)
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

  // reset all counters
  fRunNumber     = CEPTrackBuffer::kdumval;
  fEventNumber   = CEPTrackBuffer::kdumval;
  fisPileup      = kFALSE;
  fGapCondition  = 0;  
  fMCProcessType = CEPTrackBuffer::kdumval;
  fMCGenerator   = "";
  fPhysel        = kFALSE;
  fisPileup      = kFALSE;

  // clear the track list
  fCEPTracks->Clear();
  fnumTracks     = 0;
  fnumSoftTracks = 0;
  fnumResiduals  = 0;
  
  fVertexPos     = TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval);
  fMCVertexPos   = TVector3(CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval,CEPTrackBuffer::kdumval);
    
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
