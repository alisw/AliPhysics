// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Track.h"
#include "AliL3TrackArray.h"
#include "AliL3Transform.h"

//_____________________________________________________________
// Merging Hough tracks across slices

ClassImp(AliL3HoughGlobalMerger)

AliL3HoughGlobalMerger::AliL3HoughGlobalMerger()
{
  fTracks = 0;
}

AliL3HoughGlobalMerger::AliL3HoughGlobalMerger(Int_t first,Int_t last)
{
  fNSlices = last-first+1;
  fTracks = new AliL3TrackArray*[fNSlices];
  for(Int_t i=0; i<fNSlices; i++)
    fTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
  
  
}
AliL3HoughGlobalMerger::~AliL3HoughGlobalMerger()
{
  if(fTracks)
    {
      for(Int_t i=0; i<fNSlices; i++)
	{
	  if(!fTracks[i])
	    continue;
	  delete fTracks;
	}
      delete [] fTracks;
    }
  
}

void AliL3HoughGlobalMerger::FillTracks(AliL3TrackArray *tracks,Int_t slice)
{
  
  fTracks[slice]->AddTracks(tracks,kTRUE,slice);
  
}

void AliL3HoughGlobalMerger::Merge()
{
  for(Int_t slice=0; slice<fNSlices; slice++)
    {
      if(slice+1 == fNSlices) continue;
      AliL3TrackArray *t1 = fTracks[slice];
      //AliL3TrackArray *t2 = fTracks[slice+1];
      Float_t angle = AliL3Transform::Pi()/18;
      AliL3Transform::Local2GlobalAngle(&angle,slice);
      
      for(Int_t i=0; i<t1->GetNTracks(); i++)
	{
	  
	}
    }
  
}
