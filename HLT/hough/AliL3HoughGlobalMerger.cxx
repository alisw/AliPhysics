// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 

#include <math.h>
#include <iostream.h>
#include "AliL3Logging.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Track.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

//_____________________________________________________________
// Merging Hough tracks across slices

ClassImp(AliL3HoughGlobalMerger)

AliL3HoughGlobalMerger::AliL3HoughGlobalMerger(){
  //Default constructor
  Is2Global(kTRUE);
  SetParameter(2,2,0.001,0.05,0.1);
}


AliL3HoughGlobalMerger::AliL3HoughGlobalMerger(Int_t first,Int_t last) : AliL3Merger(last-first+1,"AliL3HoughTrack")
{
  //Constructor.
  fNSlices = last-first+1;
  fFirst = first;
  fLast = last;
  Is2Global(kTRUE);
  SetParameter(2,2,0.001,0.05,0.1);
}

AliL3HoughGlobalMerger::~AliL3HoughGlobalMerger(){
  //Destructor
}


void AliL3HoughGlobalMerger::FillTracks(AliL3TrackArray *tracks,Int_t slice)
{
  fSlice = slice;
  fCurrentTracks = fSlice - fFirst; 
  if(tracks->GetNTracks()==0)
    LOG(AliL3Log::kWarning,"AliL3HoughGlobalMerger::FillTracks","Track Array")
      <<AliL3Log::kDec<<"Adding empty track array in slice "<<fSlice<<ENDLOG;
  
  GetInTracks(fCurrentTracks)->AddTracks(tracks,kFALSE,fSlice);//Copy tracks, and rotate them to global coordinates
}

Bool_t AliL3HoughGlobalMerger::IsTrack(AliL3Track *innertrack,AliL3Track *outertrack)
{
  //Check if the tracks can be merged, called by the track merger
  
  AliL3HoughTrack *tr1 = (AliL3HoughTrack*)innertrack;
  AliL3HoughTrack *tr2 = (AliL3HoughTrack*)outertrack;
  
  if( (!tr1->IsPoint()) || (!tr2->IsPoint()) )  return kFALSE; 
  if(abs(tr1->GetEtaIndex() - tr2->GetEtaIndex()) > 1) return kFALSE;
  if(tr1->GetCharge() != tr2->GetCharge()) return kFALSE;
  if(fabs(tr1->GetPhi0() - tr2->GetPhi0()) > fMaxPhi0) return kFALSE;
  if(fabs(tr1->GetKappa() - tr2->GetKappa()) > fMaxKappa) return kFALSE;
  
  return kTRUE;
}

AliL3Track *AliL3HoughGlobalMerger::MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack)
{
  //Called by the track merger

  AliL3HoughTrack *newtrack = (AliL3HoughTrack*)mergedtrack->NextTrack();
  AliL3HoughTrack **trs = (AliL3HoughTrack**)tracks;
  Int_t weight=0;

  //Sum up the total weight:
  for(Int_t i=ntrack-1; i>=0; i--)
    weight += trs[i]->GetWeight();
  
  AliL3HoughTrack *tpt=trs[0];//This is the innermost track
  AliL3HoughTrack *tpl=trs[ntrack-1];
  newtrack->SetTrackParameters(tpt->GetKappa(),tpt->GetPhi0(),weight);
  newtrack->SetEtaIndex(tpt->GetEtaIndex());
  newtrack->SetEta(tpt->GetEta());
  newtrack->SetPsi(tpt->GetPsi());
  newtrack->SetCenterX(tpt->GetCenterX());
  newtrack->SetCenterY(tpt->GetCenterY());
  newtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
  newtrack->SetLastPoint(tpl->GetLastPointX(),tpl->GetLastPointY(),tpl->GetLastPointZ());
  newtrack->SetCharge(tpt->GetCharge());
  newtrack->SetRowRange(tpt->GetFirstRow(),tpl->GetLastRow());
  
  return (AliL3Track*)newtrack;

}

void AliL3HoughGlobalMerger::SlowMerge(){
  void* ntuple=GetNtuple();
  AliL3Track *track[2];
  AliL3TrackArray *tout = GetOutTracks();
  if(fNSlices<2){
    LOG(AliL3Log::kWarning,"AliL3HoughGlobalMerger::SlowMerge","Slice Number")
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
    fTransformer->Local2GlobalAngle(&angle,slice);
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
          Double_t diff = TrackDiff(track0,track1);
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
  LOG(AliL3Log::kInformational,"AliL3HoughGlobalMerger::SlowMerge","Result")
  <<AliL3Log::kDec<<"Merged Tracks: "<<tout->GetNTracks()<<" at:"
  <<angle<<ENDLOG;
  }
  WriteNtuple("ntuple_s.root",ntuple);
}

void AliL3HoughGlobalMerger::Merge()
{
  AliL3Track *track[2];
  AliL3TrackArray *tout = GetOutTracks();
  if(fNSlices<2){
    LOG(AliL3Log::kWarning,"AliL3HoughGlobalMerger::Merge","Slice Number")
      <<"Need more than one Slice!"<<ENDLOG;
    return;
  }
  Bool_t *ismatched0=0;
  Bool_t *ismatched1=0;
  for(Int_t i=0; i<fNSlices; i++)
    {
      if(fNSlices!=18 && i+1 == fNSlices) continue; //full cicle == 18 slices
      Int_t slice = fFirst + i;
      AliL3TrackArray *ttt0=GetInTracks(i);
      Int_t slice2 = i+1;
      if(slice2==fNSlices) slice2 =0;
      AliL3TrackArray *ttt1=GetInTracks(slice2);
      Float_t angle = PI/18.; //10 degrees -> the border of the slices
      fTransformer->Local2GlobalAngle(&angle,slice);
      if(i==0)
	ttt0->QSort();
      ttt1->QSort();
      if(ismatched0) delete [] ismatched0;
      if(ismatched1) delete [] ismatched1;
      ismatched0  = new Bool_t[ttt0->GetNTracks()];
      ismatched1  = new Bool_t[ttt1->GetNTracks()];
      Int_t n0=0,n1=0;
      for(Int_t s0=0;s0<ttt0->GetNTracks();s0++)
	{
	  ismatched0[s0]=kFALSE;
	  AliL3Track *track0=ttt0->GetCheckedTrack(s0);
	  if(!track0) continue;
	  track0->CalculateHelix();
	  track0->CalculateEdgePoint(angle);
	  if(track0->IsPoint()) {n0++;track0->CalculateReferencePoint(angle);
	  }
	}
      for(Int_t s1=0;s1<ttt1->GetNTracks();s1++)
	{
	  ismatched1[s1]=kFALSE;
	  AliL3Track *track1=ttt1->GetCheckedTrack(s1);
	  if(!track1) continue;
	  track1->CalculateHelix();
	  track1->CalculateEdgePoint(angle);
	  if(track1->IsPoint()) {n1++;track1->CalculateReferencePoint(angle);
	  }
	}
      for(Int_t s0=0;s0<ttt0->GetNTracks();s0++)
	{
	  AliL3Track *track0=ttt0->GetCheckedTrack(s0);
	  if(!track0) continue;
	  if(!track0->IsPoint()) continue;
	  for(Int_t s1=0;s1<ttt1->GetNTracks();s1++)
	    {
	      if(ismatched1[s1]) continue;
	      AliL3Track *track1=ttt1->GetCheckedTrack(s1);
	      if(!track1) continue;
	      if(!track1->IsPoint()) continue;
	      if(IsTrack(track0,track1))
		{
		  track[0] = track0;
		  track[1] = track1;
		  SortGlobalTracks(track,2);
		  Double_t r0 = pow(track[0]->GetLastPointX(),2) + pow(track[0]->GetLastPointY(),2);
		  Double_t r1 = pow(track[1]->GetFirstPointX(),2) + pow(track[1]->GetFirstPointY(),2);
		  if(r0<r1)
		    {
		      MultiMerge(tout,track,2); 
		      ismatched0[s0]=kTRUE;
		      ismatched1[s1]=kTRUE;
		      ttt0->Remove(s0);
		      ttt1->Remove(s1);
		    }
		}
	    }
	}
      LOG(AliL3Log::kInformational,"AliL3HoughGlobalMerger::Merge","Result")
	<<AliL3Log::kDec<<"slice0: "<<n0<<" slice1: "<<n1
	<<" Merged Tracks: "<<tout->GetNTracks()<<ENDLOG;
    }
  if(ismatched0) delete [] ismatched0;
  if(ismatched1) delete [] ismatched1;
}
