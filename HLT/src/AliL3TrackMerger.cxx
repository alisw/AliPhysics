
//Author:        Uli Frankenfeld
//Last Modified: 06.12.2000

#include <math.h>
#include <iostream.h>
#include "AliL3Logging.h"
#include "AliL3TrackMerger.h"
#include "AliL3Track.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

//_____________________________________________________________
//
// The L3 track segment merger
//

ClassImp(AliL3TrackMerger)

AliL3TrackMerger::AliL3TrackMerger(){
  //Default constructor
  Is2Global(kFALSE);
  fSlow = kFALSE;
  SetParameter();
//  fRowMin = 0;
//  fRowMax = 0;
}


AliL3TrackMerger::AliL3TrackMerger(Int_t nsubsectors):AliL3Merger(nsubsectors){
  //Constructor.
  fNSubSector = nsubsectors;
  Is2Global(kFALSE);
  fSlow = kFALSE;
  SetParameter();
  //fRowMin = new Int_t(nsubsectors);
  //fRowMax = new Int_t(nsubsectors);
}

AliL3TrackMerger::~AliL3TrackMerger(){
  //Destructor
}

void AliL3TrackMerger::SetRows(Int_t *row){
  for(Int_t i =0;i<fNSubSector;i++){
    fRowMin[i]=*(row+(2*i));
    fRowMax[i]=*(row+(2*i+1));
  }
}

void AliL3TrackMerger::InitSector(Int_t slice,Int_t subsector){
  // 
  // Select Sector and subsector. The following FillTracks call will 
  // fill this subsector
  //
  fSlice = slice;
  fSubSector = subsector;
  fCurrentTracks = fSubSector;
}

void AliL3TrackMerger::SlowMerge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout,Double_t xval){
  void *ntuple=GetNtuple();
  const Int_t  kNOut=tracksout->GetNTracks();
  const Int_t  kNIn =tracksin->GetNTracks();
  const Int_t  kNMerged =mergedtrack->GetNTracks();
  AliL3Track *tracks[2];
  Bool_t merge = kTRUE;
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNOut;out++){
    AliL3Track *outertrack=tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        AliL3Track *innertrack=tracksin->GetCheckedTrack(in);
        if(!innertrack) continue;
        Double_t diff = TrackDiff(innertrack,outertrack);
        if(diff>=0&&diff<min){
          min=diff;
          inmin=in;
          outmin=out; 
        }
      } 
    }
    if(inmin>=0&&outmin>=0){
      AliL3Track *outertrack=tracksout->GetTrack(outmin);
      AliL3Track *innertrack=tracksin->GetTrack(inmin);
      tracks[0]=innertrack;
      tracks[1]=outertrack;
      SortTracks(tracks,2);
      MultiMerge(mergedtrack,tracks,2);
      outertrack->CalculatePoint(xval);
      innertrack->CalculatePoint(xval);
      PrintDiff(innertrack,outertrack);
      FillNtuple(ntuple,innertrack,outertrack);
      tracksout->Remove(outmin);
      tracksin->Remove(inmin);
//      tracksout->Compress();
//      tracksin->Compress(); 
    }
    else merge = kFALSE;
  }
  LOG(AliL3Log::kInformational,"AliL3TrackMerger::SlowMerge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "
  <<mergedtrack->GetNTracks()-kNMerged<<ENDLOG;
  char name[256] = "ntuple_t.root";
  for(Int_t i=0;i<4;i++)
    if(tracksin==GetInTracks(i))
      sprintf(name,"ntuple_t_%d.root",i);
  WriteNtuple(name,ntuple);
}

void AliL3TrackMerger::SlowMerge(){
  fSlow = kTRUE;
  Merge();
}

void AliL3TrackMerger::InterMerge(){

  for(Int_t patch=0;patch< GetNIn();patch++){
    AliL3TrackArray * tracks = GetInTracks(patch);
    Double_t xval = fTransformer->Row2X((fRowMax[patch]+fRowMin[patch])/2);
    Int_t nrow= fRowMax[patch]-fRowMin[patch]+1;
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
    LOG(AliL3Log::kInformational,"AliL3TrackMerger::InterMerge","Result")
    <<AliL3Log::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  }
}

void AliL3TrackMerger::Merge(){
  //Loop over tracks and pass them to the track merger.
  Double_t edge0 = PI/18.;
  Double_t edge1 = 2*PI - edge0;
  AliL3TrackArray *ttt = GetOutTracks();
  if(fNSubSector==1) {
    GetOutTracks()->AddTracks(GetInTracks(0)); 
    LOG(AliL3Log::kInformational,"AliL3TrackMerger::Merge","Result")
    <<AliL3Log::kDec<<"Total Copied Tracks: "<<GetOutTracks()->GetNPresent()
    <<ENDLOG;
    return;
  }
  Int_t subsec = fNSubSector -2; 
  for(Int_t i=subsec;i>=0;i--){
    AliL3TrackArray *tout = GetOutTracks();
    if(i==subsec) tout = GetInTracks(subsec+1);
    AliL3TrackArray *tin = GetInTracks(i);
    Double_t xval = fTransformer->Row2X(fRowMax[i]);
    Double_t xmax = fTransformer->Row2X(fRowMax[i+1]);
    Double_t ymax = xval*tan(edge0);
    for(Int_t out=0;out<tout->GetNTracks();out++){
      AliL3Track *outtrack=tout->GetCheckedTrack(out);
      if(!outtrack) continue;
      outtrack->CalculateHelix();
      outtrack->CalculatePoint(xval);
      if(outtrack->IsPoint()&&fabs(outtrack->GetPointY())>ymax){
        if(outtrack->GetNHits()<10)
          tout->Remove(out);
      }
    }
//    tout->Compress();
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliL3Track *intrack=(AliL3Track*)tin->GetTrack(in);
      intrack->CalculateHelix();
      intrack->CalculatePoint(xval);
    }
    tin->QSort();
    tout->QSort();

    if(fSlow) SlowMerge(ttt,tin,tout,xval);
    else Merge(ttt,tin,tout);
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliL3Track *intrack=(AliL3Track*)tin->GetCheckedTrack(in);
      if(!intrack) continue;
      if(intrack->CalculateEdgePoint(edge0)){
        if(intrack->GetPointX()<xmax ){
          AddTrack(ttt,intrack);
          tin->Remove(in);
        }
      } 
      else if(intrack->CalculateEdgePoint(edge1)){
        if(intrack->GetPointX()<xmax ){
          AddTrack(ttt,intrack);
          tin->Remove(in);
        }
      }
    }
/*
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliL3Track *intrack=(AliL3Track*)tin->GetCheckedTrack(in);
      if(!intrack) continue;
      if(intrack->GetNHits()<10) continue;
      AddTrack(ttt,intrack);
      tin->Remove(in);
    }
*/
  } // end subsector loop
  LOG(AliL3Log::kInformational,"AliL3TrackMerger::Merge","Result")
  <<AliL3Log::kDec<<"Total Merged Tracks: "<<GetOutTracks()->GetNPresent()
  <<ENDLOG;
}

Int_t AliL3TrackMerger::Merge(AliL3TrackArray* mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout){
  //Loop over tracks and pass them to the track merger.
  AliL3Track *tracks[2];

  const Int_t  kNOut=tracksout->GetNTracks();
  const Int_t  kNIn =tracksin->GetNTracks();
  const Int_t  kNMerged =mergedtrack->GetNTracks();

  Bool_t *ismatchedin  = new Bool_t[kNIn];
  for(Int_t in =0;in<kNIn;in++)
    ismatchedin[in]=kFALSE;
  Bool_t *ismatchedout = new Bool_t[kNOut];
  for(Int_t out =0;out<kNOut;out++)
    ismatchedout[out] = kFALSE;
  for(Int_t out =0;out<kNOut;out++){
    AliL3Track *outertrack=(AliL3Track*)tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
    for(Int_t in =0;in<kNIn;in++){
      if(ismatchedin[in]) continue;
      AliL3Track *innertrack=(AliL3Track*)tracksin->GetCheckedTrack(in);
      if(!innertrack) continue;
      if(outertrack==innertrack) continue;
      if(outertrack->GetCharge()!=innertrack->GetCharge()) continue;
      if(IsTrack(innertrack,outertrack)){
        tracks[0]=innertrack; tracks[1]=outertrack; 
        SortTracks(tracks,2);  
        if(tracks[0]->GetLastPointX()<tracks[1]->GetFirstPointX()){
          MultiMerge(mergedtrack,tracks,2);
          tracksout->Remove(out);
          tracksin->Remove(in);
          ismatchedin[in]=kTRUE;
          ismatchedout[out]=kTRUE;
          break;
        }
      }
    }
  }

  Int_t nmerged = mergedtrack->GetNTracks()-kNMerged;
  LOG(AliL3Log::kInformational,"AliL3TrackMerger::Merge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  delete[] ismatchedin;
  delete[] ismatchedout;
  return nmerged;
}


