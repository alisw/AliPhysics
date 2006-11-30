// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTHoughGlobalMerger.h"
#include "AliHLTTrack.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTransform.h"

//_____________________________________________________________
// Merging Hough tracks across slices

ClassImp(AliHLTHoughGlobalMerger)

AliHLTHoughGlobalMerger::AliHLTHoughGlobalMerger()
{
  fTracks = 0;
}

AliHLTHoughGlobalMerger::AliHLTHoughGlobalMerger(Int_t first,Int_t last)
{
  fNSlices = last-first+1;
  fTracks = new AliHLTTrackArray*[fNSlices];
  for(Int_t i=0; i<fNSlices; i++)
    fTracks[i] = new AliHLTTrackArray("AliHLTHoughTrack");
  
  
}
AliHLTHoughGlobalMerger::~AliHLTHoughGlobalMerger()
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

void AliHLTHoughGlobalMerger::FillTracks(AliHLTTrackArray *tracks,Int_t slice)
{
  
  fTracks[slice]->AddTracks(tracks,kTRUE,slice);
  
}

void AliHLTHoughGlobalMerger::Merge()
{
  for(Int_t slice=0; slice<fNSlices; slice++)
    {
      if(slice+1 == fNSlices) continue;
      AliHLTTrackArray *t1 = fTracks[slice];
      //AliHLTTrackArray *t2 = fTracks[slice+1];
      Float_t angle = AliHLTTransform::Pi()/18;
      AliHLTTransform::Local2GlobalAngle(&angle,slice);
      
      for(Int_t i=0; i<t1->GetNTracks(); i++)
	{
	  
	}
    }
  
}
