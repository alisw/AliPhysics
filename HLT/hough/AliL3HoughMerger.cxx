
#include "AliL3Logging.h"
#include "AliL3Defs.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughMerger.h"
#include "AliL3HoughTransformer.h"

ClassImp(AliL3HoughMerger)

  
AliL3HoughMerger::AliL3HoughMerger()
{
  //Default constructor

}


AliL3HoughMerger::AliL3HoughMerger(Int_t nsubsectors)
{
  //Constructor
  fNIn = nsubsectors;
  SetArray();
}


AliL3HoughMerger::~AliL3HoughMerger()
{
  DeleteArray();
  if(fOutTracks)
    delete fOutTracks;
}

void AliL3HoughMerger::DeleteArray()
{
  if(!fInTracks)
    return;
  for(Int_t i=0; i<fNIn; i++)
    {
      if(!fInTracks[i]) continue;
      delete fInTracks;
    }
  //delete [] fInTracks;
}

void AliL3HoughMerger::SetArray()
{
  fInTracks = new AliL3TrackArray*[fNIn];
  for(Int_t i=0; i<fNIn; i++)
    fInTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
  fOutTracks = new AliL3TrackArray("AliL3HoughTrack");
}

void AliL3HoughMerger::FillTracks(AliL3TrackArray *tracks,Int_t patch)
{
  if(tracks->GetNTracks()==0)
    LOG(AliL3Log::kWarning,"AliL3HoughMerger::FillTracks","Track Array")
      <<"Adding empty track array"<<ENDLOG;
  fInTracks[patch]->AddTracks(tracks,kFALSE); //Copy tracks
}

void AliL3HoughMerger::MergePatches()
{
  
  
}

void AliL3HoughMerger::MergeEtaSlices(Int_t patch)
{
  AliL3TrackArray *tracks = fInTracks[patch];
  Int_t ntracks = tracks->GetNTracks();
  printf("Number of tracks to merging %d\n",ntracks);
  AliL3HoughTrack *tr1,*tr2;
  Double_t ptdiff=0.01;
  Double_t phi0diff=0.01;
  for(Int_t i=0; i<ntracks; i++)
    {
      tr1 = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr1) continue;
      for(Int_t j=0; j<ntracks; j++)
	{
	  tr2 = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  if(!tr2) continue;
	  if(tr1==tr2) continue;
	  if(tr1->GetEtaIndex() == tr2->GetEtaIndex()) continue;
	  if(tr1->GetEtaIndex() == tr2->GetEtaIndex()-1 || tr1->GetEtaIndex() == tr2->GetEtaIndex()+1)
	    if(fabs(tr1->GetPt()-tr2->GetPt())<ptdiff && fabs(tr1->GetPhi0()-tr2->GetPhi0())<phi0diff) 
	      MergeTracks(tracks,i,j);
	}
    }
  tracks->Compress();
  
  printf("Number of tracks that was not merged: %d\n",tracks->GetNTracks());
  //Add rest of the tracks, which was not merged
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      tr1 = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr1) continue;
      tr2 = (AliL3HoughTrack*)fOutTracks->NextTrack();
      tr2->Set(tr1);
    }
  printf("Total number of tracks after merging: %d\n",fOutTracks->GetNTracks());
}

void AliL3HoughMerger::MergeTracks(AliL3TrackArray *intracks,Int_t i,Int_t j)
{
  AliL3HoughTrack *newtrack = (AliL3HoughTrack*)fOutTracks->NextTrack();
  AliL3HoughTrack *mergetrack = (AliL3HoughTrack*)intracks->GetCheckedTrack(i);
  AliL3HoughTrack *mergetrack2 = (AliL3HoughTrack*)intracks->GetCheckedTrack(j);
  if(!mergetrack || !mergetrack2)
    {
      printf("\nALiL3HoughMerger::MergeTracks : NO TRACK!!!\n");
      return;
    }
  Int_t w = mergetrack->GetWeight() + mergetrack2->GetWeight();
  newtrack->SetTrackParameters(mergetrack->GetKappa(),mergetrack->GetPhi0(),w);
  newtrack->SetEtaIndex(mergetrack->GetEtaIndex());
  newtrack->SetEta(mergetrack->GetEta());
  
  intracks->Remove(i);
  intracks->Remove(j);
  
}
