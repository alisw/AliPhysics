//$Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy Uli 

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3GlobalMerger.h"
#include "AliL3Track.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

//_____________________________________________________________
// AliL3GlobalMerger
//
// The L3 Slice merger
//

ClassImp(AliL3GlobalMerger)

AliL3GlobalMerger::AliL3GlobalMerger(){
  //Default constructor
  Is2Global(kTRUE);
  SetParameter(2,2,0.001,0.05,0.1);
}


AliL3GlobalMerger::AliL3GlobalMerger(Int_t first,Int_t last):AliL3Merger(last-first+1){
  //Constructor.
  fNSlices = last-first+1;
  fFirst = first;
  fLast = last;
  Is2Global(kTRUE);
  SetParameter(2,2,0.001,0.05,0.1);
}

AliL3GlobalMerger::~AliL3GlobalMerger(){
  //Destructor
}

void AliL3GlobalMerger::InitSlice(Int_t slice){
  // 
  // Select Sector The following FillTracks call will 
  // fill this Sector
  //
  fSlice = slice;
  fCurrentTracks = fSlice - fFirst; 
}

Double_t AliL3GlobalMerger::CheckTracks(AliL3Track *innertrack,AliL3Track *outertrack,Int_t slice)
{
  //Compare the tracks by propagating the outermost track to the last and first point plane
  //of the innermost track. This plane is defined by the padrow plane where these points
  //are.
  
  
  Float_t angle = PI/36;//5 degrees = perpendicular to padrowplane (in local system)
  AliL3Transform::Local2GlobalAngle(&angle,slice);
  Double_t dx[2],dy[2],dz[2];
  Double_t diff =-1;
  AliL3Track *tracks[2];
  tracks[0] = innertrack;
  tracks[1] = outertrack;
  SortGlobalTracks(tracks,2);
  innertrack = tracks[0]; 
  outertrack = tracks[1];
  
  Float_t point[3];
  
  point[0]=innertrack->GetLastPointX();
  point[1]=innertrack->GetLastPointY();
  point[2]=innertrack->GetLastPointZ();
  AliL3Transform::Global2Local(point,slice,kTRUE);
  
  outertrack->CalculateReferencePoint(angle,point[0]);//local x = global distance to padrowplane
  if(!outertrack->IsPoint()) return diff;
  dx[0] = fabs(outertrack->GetPointX()-innertrack->GetLastPointX());
  dy[0] = fabs(outertrack->GetPointY()-innertrack->GetLastPointY());
  dz[0] = fabs(outertrack->GetPointZ()-innertrack->GetLastPointZ());
  
  point[0]=innertrack->GetFirstPointX();
  point[1]=innertrack->GetFirstPointY();
  point[2]=innertrack->GetFirstPointZ();
  AliL3Transform::Global2Local(point,slice,kTRUE);
  
  outertrack->CalculateReferencePoint(angle,point[0]);//local x = global distance to padrowplane
  if(!outertrack->IsPoint()) return diff;
  dx[1] = fabs(outertrack->GetPointX()-innertrack->GetFirstPointX());
  dy[1] = fabs(outertrack->GetPointY()-innertrack->GetFirstPointY());
  dz[1] = fabs(outertrack->GetPointZ()-innertrack->GetFirstPointZ());
  
  for(Int_t i=0; i<2; i++)
    diff += sqrt(dx[i]*dx[i] + dy[i]*dy[i] + dz[i]*dz[i]);
  return diff;
}

void AliL3GlobalMerger::SlowMerge(){
  void* ntuple=GetNtuple();
  AliL3Track *track[2];
  AliL3TrackArray *tout = GetOutTracks();
  if(fNSlices<2){
    LOG(AliL3Log::kWarning,"AliL3GlobalMerger::SlowMerge","Slice Number")
    <<"Need more than one Slice!"<<ENDLOG;
    return;
  }
  for(Int_t i=0; i<fNSlices; i++){
    if(fNSlices!=18 && i+1 == fNSlices) continue; //full cicle == 18 slices
    Int_t slice = fFirst + i;
    AliL3TrackArray *ttt0=GetInTracks(i);
    Int_t slice2 = i+1;
    if(slice2==fNSlices) slice2 =0; 
    AliL3TrackArray *ttt1=GetInTracks(slice2);
    Float_t angle = PI/18.; //10 degrees -> the border of the slices
    AliL3Transform::Local2GlobalAngle(&angle,slice);
    if(i==0)
      ttt0->QSort();
    ttt1->QSort();
    for(Int_t s0=0;s0<ttt0->GetNTracks();s0++){
      AliL3Track *track0=ttt0->GetCheckedTrack(s0);
      if(!track0) continue;
      track0->CalculateHelix();
      track0->CalculateEdgePoint(angle);
      //      if(track0->IsPoint()) AddTrack(tout,track0);
    }
    for(Int_t s1=0;s1<ttt1->GetNTracks();s1++){
      AliL3Track *track1=ttt1->GetCheckedTrack(s1);
      if(!track1) continue;
      track1->CalculateHelix();
      track1->CalculateEdgePoint(angle);
      //      if(track1->IsPoint())  AddTrack(tout,track1); 
    }
    Bool_t merge = kTRUE;
    while(merge){
      Int_t min0=-1,min1=-1;
      Double_t min=10;
      for(Int_t s0=0;s0<ttt0->GetNTracks();s0++){
        AliL3Track *track0=ttt0->GetCheckedTrack(s0);
        if(!track0) continue;
        if(!track0->IsPoint()) continue;
        for(Int_t s1=0;s1<ttt1->GetNTracks();s1++){
          AliL3Track *track1=ttt1->GetCheckedTrack(s1);
          if(!track1) continue;
          if(!track1->IsPoint()) continue;
	  
          //Double_t diff = TrackDiff(track0,track1,angle);
	  Double_t diff = CheckTracks(track0,track1,slice);
	  //PrintDiff(track0,track1);
	  if(diff>=0&&diff<min){
            min=diff;
            min0=s0;
            min1=s1;
          }
        }
      }
      if(min0>=0&&min1>=0){
        AliL3Track *track0=ttt0->GetTrack(min0);
        AliL3Track *track1=ttt1->GetTrack(min1);
	track[0] = track0;
        track[1] = track1;
        SortGlobalTracks(track,2);
	track1->CalculateEdgePoint((angle+PI/9));
	if(track1->IsPoint())//Check if the track will cross the boundary of yet another slice.
	  MultiMerge(ttt1,track,2);
	else
	  MultiMerge(tout,track,2); 
        track0->CalculateReferencePoint(angle);
        track1->CalculateReferencePoint(angle);
	PrintDiff(track0,track1);
        FillNtuple(ntuple,track0,track1);
        ttt0->Remove(min0);
        ttt1->Remove(min1);
	
      }
      else merge = kFALSE;
    }
    ttt0->Compress();
    ttt1->Compress();
  LOG(AliL3Log::kInformational,"AliL3GlobalMerger::SlowMerge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "<<tout->GetNTracks()<<" at:"
  <<angle<<ENDLOG;
  }
  WriteNtuple("ntuple_s.root",ntuple);
}

void AliL3GlobalMerger::Merge(){
  AliL3Track *track[2];
  AliL3TrackArray *tout = GetOutTracks();
  if(fNSlices<2){
  LOG(AliL3Log::kWarning,"AliL3GlobalMerger::Merge","Slice Number")
    <<"Need more than one Slice!"<<ENDLOG;
    return;
  }
  for(Int_t i=0; i<fNSlices; i++){
    if(fNSlices!=18 && i+1 == fNSlices) continue; //full cicle == 18 slices
    Int_t slice = fFirst + i;
    AliL3TrackArray *ttt0=GetInTracks(i);
    Int_t slice2 = i+1;
    if(slice2==fNSlices) slice2 =0;
    AliL3TrackArray *ttt1=GetInTracks(slice2);
    Float_t angle = PI/18.; //10 degrees -> the border of the slices
    AliL3Transform::Local2GlobalAngle(&angle,slice);
    if(i==0)
      ttt0->QSort();
    ttt1->QSort();
    Bool_t *ismatched0  = new Bool_t[ttt0->GetNTracks()];
    Bool_t *ismatched1  = new Bool_t[ttt1->GetNTracks()];
    Int_t n0=0,n1=0;
    for(Int_t s0=0;s0<ttt0->GetNTracks();s0++){
      ismatched0[s0]=kFALSE;
      AliL3Track *track0=ttt0->GetCheckedTrack(s0);
      if(!track0) continue;
      track0->CalculateHelix();
      track0->CalculateEdgePoint(angle);
      if(track0->IsPoint()) {n0++;track0->CalculateReferencePoint(angle);}
    }
    for(Int_t s1=0;s1<ttt1->GetNTracks();s1++){
      ismatched1[s1]=kFALSE;
      AliL3Track *track1=ttt1->GetCheckedTrack(s1);
      if(!track1) continue;
      track1->CalculateHelix();
      track1->CalculateEdgePoint(angle);
      if(track1->IsPoint()) {n1++;track1->CalculateReferencePoint(angle);}
    }
    for(Int_t s0=0;s0<ttt0->GetNTracks();s0++){
      AliL3Track *track0=ttt0->GetCheckedTrack(s0);
      if(!track0) continue;
      if(!track0->IsPoint()) continue;
      for(Int_t s1=0;s1<ttt1->GetNTracks();s1++){
        if(ismatched1[s1]) continue;
        AliL3Track *track1=ttt1->GetCheckedTrack(s1);
        if(!track1) continue;
        if(!track1->IsPoint()) continue;
        if(IsRTrack(track0,track1)){
          track[0] = track0;
          track[1] = track1;
          SortGlobalTracks(track,2);
          Double_t r0 = pow(track[0]->GetLastPointX(),2)+
                        pow(track[0]->GetLastPointY(),2);
          Double_t r1 = pow(track[1]->GetFirstPointX(),2)+
                        pow(track[1]->GetFirstPointY(),2);
          if(r0<r1){
            MultiMerge(tout,track,2); 
            ismatched0[s0]=kTRUE;
            ismatched1[s1]=kTRUE;
            ttt0->Remove(s0);
            ttt1->Remove(s1);
          }
        }
      }
    }
  LOG(AliL3Log::kInformational,"AliL3GlobalMerger::Merge","Result")
  <<AliL3Log::kDec<<"slice0: "<<n0<<" slice1: "<<n1
  <<" Merged Tracks: "<<tout->GetNTracks()<<ENDLOG;
  }
}
