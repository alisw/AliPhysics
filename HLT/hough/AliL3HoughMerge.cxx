#include <TH2.h>

#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughMerge.h"
#include "AliL3HoughTransformer.h"

ClassImp(AliL3HoughMerge)

  
AliL3HoughMerge::AliL3HoughMerge()
{
  //Default constructor

  fInTracks = NULL;
  fOutTracks = NULL;
  fNPatches = 5;
}


AliL3HoughMerge::AliL3HoughMerge(Int_t slice,Int_t row_patches)
{
  //Constructor

  //fInTracks = (AliL3TrackArray**)new Byte_t[row_patches*sizeof(AliL3TrackArray*)];
  fInTracks = new AliL3TrackArray*[row_patches];
  fNPatches = row_patches;
  for(Int_t i=0; i<row_patches; i++)
    fInTracks[i] = new AliL3TrackArray("AliL3HoughTrack");

  fOutTracks = new AliL3TrackArray("AliL3HoughTrack");
}


AliL3HoughMerge::~AliL3HoughMerge()
{
  //Destructor
  if(fInTracks)
    delete fInTracks;
  if(fOutTracks)
    delete fOutTracks;
  
}

void AliL3HoughMerge::FillTracks(AliL3TrackArray *tracks,Int_t patch)
{
  fInTracks[patch]->AddTracks(tracks); //copies tracks to new trackarray. Does not delete the track objects.
}

void AliL3HoughMerge::FillHisto(TH2F *merge_hist)
{
  
  for(Int_t pat=0; pat < fNPatches; pat++)
    {
      for(Int_t t=0; t<fInTracks[pat]->GetNTracks(); t++)
	{
	  AliL3HoughTrack *tr = (AliL3HoughTrack*)fInTracks[pat]->GetCheckedTrack(t);
	  if(!tr) {printf("AliL3HoughMerge NO TRACK\n"); continue;}
	  merge_hist->Fill(tr->GetKappa(),tr->GetPhi0(),tr->GetNHits());
	}
    }
}

/*
void AliL3HoughMerge::MergeLines()
{

  

}
*/
