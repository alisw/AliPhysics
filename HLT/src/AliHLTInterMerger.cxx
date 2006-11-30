// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTInterMerger.h"
#include "AliHLTTrack.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTTransform.h"
#include "AliHLTTrackArray.h"

/** \class AliHLTInterMerger
<pre>
//_____________________________________________________________
// AliHLTInterMerger
//
// The L3 track segment merger
//
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTInterMerger)

AliHLTInterMerger::AliHLTInterMerger()
{
  //Default constructor
  InitMerger(1);
  Is2Global(kFALSE);
//  SetParameter(2,2,.3,.3,.3);
  SetParameter(1,0.5,0.0005,0.05,0.1);
  fRowMax = fRowMin = 0;
}


AliHLTInterMerger::~AliHLTInterMerger(){
  //Destructor
  
}

void AliHLTInterMerger::SlowMerge(){
  Int_t nrow= fRowMax-fRowMin+1;
  void *ntuple=GetNtuple();
  AliHLTTrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliHLTTrack *tr[2];
  Bool_t merge = kTRUE;
  for(Int_t in=0;in<kNIn;in++)
    tracks->GetCheckedTrack(in)->CalculateHelix();
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNIn;out++){
    AliHLTTrack *outertrack=tracks->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        if(in==out) continue;
        AliHLTTrack *innertrack=tracks->GetCheckedTrack(in);
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
      AliHLTTrack *outertrack=tracks->GetTrack(outmin);
      AliHLTTrack *innertrack=tracks->GetTrack(inmin);
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
  LOG(AliHLTLog::kInformational,"AliHLTInterMerger::SlowMerge","Result")
  <<AliHLTLog::kDec<<"Merged Tracks: "<<tracks->GetNTracks()-kNIn<<ENDLOG;

  char name[256];
  sprintf(name,"ntuple_i_%d.root",fPatch);
  WriteNtuple(name,ntuple);
}

void AliHLTInterMerger::MMerge(){
  while(Merge());
  GetOutTracks()->AddTracks(GetInTracks(0));
}

Int_t AliHLTInterMerger::Merge(){
  Int_t nrow= fRowMax-fRowMin+1;
  Double_t xval =AliHLTTransform::Row2X((fRowMax+fRowMin)/2);
  AliHLTTrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliHLTTrack *tr[2];
  for(Int_t in=0;in<kNIn;in++){
    AliHLTTrack *t = tracks->GetCheckedTrack(in);
    if(t){
      t->CalculateHelix();
      t->CalculatePoint(xval);
    }
  }
  for(Int_t out=0;out<kNIn;out++){
  AliHLTTrack *outertrack=tracks->GetCheckedTrack(out);
  if(!outertrack) continue;
    for(Int_t in=0;in<kNIn;in++){
      if(in==out) continue;
      AliHLTTrack *innertrack=tracks->GetCheckedTrack(in);
      if(!innertrack) continue;
      if(outertrack->GetNHits()+innertrack->GetNHits()>nrow) continue;

      if(IsTrack(innertrack,outertrack)){
        tr[0]=innertrack;
        tr[1]=outertrack;
        SortGlobalTracks(tr,2);

        Double_t r0 = pow(tr[0]->GetLastPointX(),2)+
                      pow(tr[0]->GetLastPointY(),2);
        Double_t r1 = pow(tr[1]->GetFirstPointX(),2)+
                      pow(tr[1]->GetFirstPointY(),2);
        if(r0<r1){
          MultiMerge(tracks,tr,2);
          tracks->Remove(out);
          tracks->Remove(in);
          break;
        }
      }
    } 
  }
  Int_t nmerged = tracks->GetNTracks()-kNIn; 
  LOG(AliHLTLog::kInformational,"AliHLTInterMerger::Merge","Result")
  <<AliHLTLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  //add in tracks
//  GetOutTracks()->AddTracks(GetInTracks(0)); 

  return nmerged;
}


