
//Author:        Uli Frankenfeld
//Last Modified: 06.12.2000

#include <math.h>
#include <iostream.h>
#include "AliL3Logging.h"
#include "AliL3InterMerger.h"
#include "AliL3Track.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

//_____________________________________________________________
//
// The L3 track segment merger
//

ClassImp(AliL3InterMerger)

AliL3InterMerger::AliL3InterMerger():AliL3Merger(1){
  //Default constructor
  Is2Global(kFALSE);
  SetParameter(2,2,.3,.3,.3);
  fRowMax = fRowMin = 0;
}


AliL3InterMerger::~AliL3InterMerger(){
  //Destructor
  
}

void AliL3InterMerger::SlowMerge(){
  Int_t nrow= fRowMax-fRowMin+1;
  void *ntuple=GetNtuple();
  AliL3TrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliL3Track *tr[2];
  Bool_t merge = kTRUE;
  for(Int_t in=0;in<kNIn;in++)
    tracks->GetCheckedTrack(in)->CalculateHelix();
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNIn;out++){
    AliL3Track *outertrack=tracks->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        if(in==out) continue;
        AliL3Track *innertrack=tracks->GetCheckedTrack(in);
        if(!innertrack) continue;
        if(outertrack->GetNHits()+innertrack->GetNHits()>nrow) continue;

        Double_t diff = TrackDiff(innertrack,outertrack);

        if(diff>=0&&diff<min){
          min=diff;
          inmin=in;
          outmin=out; 
        }
      } 
    }
    if(inmin>=0&&outmin>=0){
      AliL3Track *outertrack=tracks->GetTrack(outmin);
      AliL3Track *innertrack=tracks->GetTrack(inmin);
      tr[0]=innertrack;
      tr[1]=outertrack;
      SortTracks(tr,2);
      MultiMerge(tracks,tr,2);
      outertrack->CalculatePoint(tr[0]->GetLastPointX());
      innertrack->CalculatePoint(tr[0]->GetLastPointX());
      PrintDiff(innertrack,outertrack);
      FillNtuple(ntuple,innertrack,outertrack);
      tracks->Remove(outmin);
      tracks->Remove(inmin);
    }
    else merge = kFALSE;
  }
  LOG(AliL3Log::kInformational,"AliL3InterMerger::SlowMerge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "<<tracks->GetNTracks()-kNIn<<ENDLOG;

  char name[256];
  sprintf(name,"ntuple_i_%d.root",fPatch);
  WriteNtuple(name,ntuple);
}

void AliL3InterMerger::MMerge(){
  while(Merge());
  GetOutTracks()->AddTracks(GetInTracks(0));
}

Int_t AliL3InterMerger::Merge(){
  Int_t nrow= fRowMax-fRowMin+1;
  Double_t xval =fTransformer->Row2X((fRowMax+fRowMin)/2);
  AliL3TrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliL3Track *tr[2];
  for(Int_t in=0;in<kNIn;in++){
    AliL3Track *t = tracks->GetCheckedTrack(in);
    if(t){
      t->CalculateHelix();
      t->CalculatePoint(xval);
    }
  }
  for(Int_t out=0;out<kNIn;out++){
  AliL3Track *outertrack=tracks->GetCheckedTrack(out);
  if(!outertrack) continue;
    for(Int_t in=0;in<kNIn;in++){
      if(in==out) continue;
      AliL3Track *innertrack=tracks->GetCheckedTrack(in);
      if(!innertrack) continue;
      if(outertrack->GetNHits()+innertrack->GetNHits()>nrow) continue;

      if(IsTrack(innertrack,outertrack)){
        tr[0]=innertrack;
        tr[1]=outertrack;
        SortTracks(tr,2);
        if(tr[0]->GetLastPointX()<tr[1]->GetFirstPointX()){
          MultiMerge(tracks,tr,2);
          tracks->Remove(out);
          tracks->Remove(in);
          break;
        }
      }
    } 
  }
  Int_t nmerged = tracks->GetNTracks()-kNIn; 
  LOG(AliL3Log::kInformational,"AliL3InterMerger::Merge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  //add in tracks
//  GetOutTracks()->AddTracks(GetInTracks(0)); 

  return nmerged;
}


