// @(#) $Id$
// Original: AliHLTTrackMerger.cxx,v 1.12 2005/06/14 10:55:21 cvetan 

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

/** \class AliHLTTPCTrackMerger
<pre>
//_____________________________________________________________
// AliHLTTPCTrackMerger
//
// The HLTTPC track segment merger
//
</pre
*/

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTrackMerger.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackArray.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCTrackMerger)

AliHLTTPCTrackMerger::AliHLTTPCTrackMerger(){
  //Default constructor
  Is2Global(kFALSE);
  fSlow = kFALSE;
  SetParameter();
  fRowMin = 0;
  fRowMax = 0;
}


AliHLTTPCTrackMerger::AliHLTTPCTrackMerger(Int_t nsubsectors) : AliHLTTPCMerger()
{
  //Constructor.
  InitMerger(nsubsectors);
  fNSubSector = nsubsectors;
  Is2Global(kFALSE);
  fSlow = kFALSE;
  SetParameter();
  fRowMin = new Int_t[nsubsectors];
  fRowMax = new Int_t[nsubsectors];
  
}

AliHLTTPCTrackMerger::~AliHLTTPCTrackMerger(){
  //Destructor
}

void AliHLTTPCTrackMerger::SetRows(Int_t *row){
  //Set the indeces of the first and last
  //TPC padrows
  //
  for(Int_t i=0;i<fNSubSector;i++){
    fRowMin[i]=*(row+(2*i));
    fRowMax[i]=*(row+(2*i+1));
  }
}

void AliHLTTPCTrackMerger::InitSector(Int_t slice,Int_t subsector){
  // 
  // Select Sector and subsector. The following FillTracks call will 
  // fill this subsector
  //
  fSlice = slice;
  fSubSector = subsector;
  fCurrentTracks = fSubSector;
}

void AliHLTTPCTrackMerger::SlowMerge(AliHLTTPCTrackArray *mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout,Double_t xval){
  // 
  // Slow merging of two AliHLTTPCTrackArrays
  // at reference plane x=xval
  //
  void *ntuple=GetNtuple();
  const Int_t  kNOut=tracksout->GetNTracks();
  const Int_t  kNIn =tracksin->GetNTracks();
  const Int_t  kNMerged =mergedtrack->GetNTracks();
  AliHLTTPCTrack *tracks[2];
  Bool_t merge = kTRUE;
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNOut;out++){
    AliHLTTPCTrack *outertrack=tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        AliHLTTPCTrack *innertrack=tracksin->GetCheckedTrack(in);
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
      AliHLTTPCTrack *outertrack=tracksout->GetTrack(outmin);
      AliHLTTPCTrack *innertrack=tracksin->GetTrack(inmin);
      tracks[0]=innertrack;
      tracks[1]=outertrack;
      SortTracks(tracks,2);
      MultiMerge(mergedtrack,tracks,2);
      outertrack->CalculatePoint(xval);
      innertrack->CalculatePoint(xval);
      PrintDiff(innertrack,outertrack);
      //FillNtuple(ntuple,innertrack,outertrack);
      tracksout->Remove(outmin);
      tracksin->Remove(inmin);
//      tracksout->Compress();
//      tracksin->Compress(); 
    }
    else merge = kFALSE;
  }
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrackMerger::SlowMerge","Result")
  <<AliHLTTPCLog::kDec<<"Merged Tracks: "
  <<mergedtrack->GetNTracks()-kNMerged<<ENDLOG;
  char name[256] = "ntuple_t.root";
  for(Int_t i=0;i<4;i++)
    if(tracksin==GetInTracks(i))
      sprintf(name,"ntuple_t_%d.root",i);
  WriteNtuple(name,ntuple);
}

void AliHLTTPCTrackMerger::SlowMerge(){
  fSlow = kTRUE;
  Merge();
}

void AliHLTTPCTrackMerger::InterMerge(){
  // 
  // Merging of the tracks
  // between readout patches
  //
  for(Int_t patch=0;patch< GetNIn();patch++){
    AliHLTTPCTrackArray * tracks = GetInTracks(patch);
    Double_t xval = AliHLTTPCTransform::Row2X((fRowMax[patch]+fRowMin[patch])/2);
    Int_t nrow= fRowMax[patch]-fRowMin[patch]+1;
    const Int_t  kNIn =tracks->GetNTracks();
    AliHLTTPCTrack *tr[2];
    for(Int_t in=0;in<kNIn;in++){
      AliHLTTPCTrack *t = tracks->GetCheckedTrack(in);
      if(t){
        t->CalculateHelix();
        t->CalculatePoint(xval);
      }
    }
    for(Int_t out=0;out<kNIn;out++){
    AliHLTTPCTrack *outertrack=tracks->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        if(in==out) continue;
        AliHLTTPCTrack *innertrack=tracks->GetCheckedTrack(in);
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
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrackMerger::InterMerge","Result")
    <<AliHLTTPCLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  }
}

void AliHLTTPCTrackMerger::Merge(){
  //Loop over tracks and pass them to the track merger.
  Double_t edge0 = AliHLTTPCTransform::Pi()/18;
  Double_t edge1 = AliHLTTPCTransform::TwoPi() - edge0;
  AliHLTTPCTrackArray *ttt = GetOutTracks();
  if(fNSubSector==1) {
    GetOutTracks()->AddTracks(GetInTracks(0)); 
    LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrackMerger::Merge","Result")
    <<AliHLTTPCLog::kDec<<"Total Copied Tracks: "<<GetOutTracks()->GetNPresent()
    <<ENDLOG;
    return;
  }
  Int_t subsec = fNSubSector -2; 
  for(Int_t i=subsec;i>=0;i--){
    AliHLTTPCTrackArray *tout = GetOutTracks();
    if(i==subsec) tout = GetInTracks(subsec+1);
    AliHLTTPCTrackArray *tin = GetInTracks(i);
    Double_t xval = AliHLTTPCTransform::Row2X(fRowMax[i]);
    Double_t xmax = AliHLTTPCTransform::Row2X(fRowMax[i+1]);
    Double_t ymax = xval*tan(edge0);
    for(Int_t out=0;out<tout->GetNTracks();out++){
      AliHLTTPCTrack *outtrack=tout->GetCheckedTrack(out);
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
      AliHLTTPCTrack *intrack=(AliHLTTPCTrack*)tin->GetTrack(in);
      intrack->CalculateHelix();
      intrack->CalculatePoint(xval);
    }
    tin->QSort();
    tout->QSort();

    if(fSlow) SlowMerge(ttt,tin,tout,xval);
    else Merge(ttt,tin,tout);
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliHLTTPCTrack *intrack=(AliHLTTPCTrack*)tin->GetCheckedTrack(in);
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
      AliHLTTPCTrack *intrack=(AliHLTTPCTrack*)tin->GetCheckedTrack(in);
      if(!intrack) continue;
      if(intrack->GetNHits()<10) continue;
      AddTrack(ttt,intrack);
      tin->Remove(in);
    }
*/
  } // end subsector loop
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrackMerger::Merge","Result")
  <<AliHLTTPCLog::kDec<<"Total Merged Tracks: "<<GetOutTracks()->GetNPresent()
  <<ENDLOG;
}

Int_t AliHLTTPCTrackMerger::Merge(AliHLTTPCTrackArray* mergedtrack,AliHLTTPCTrackArray *tracksin,AliHLTTPCTrackArray *tracksout){
  //Loop over tracks and pass them to the track merger.
  AliHLTTPCTrack *tracks[2];

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
    AliHLTTPCTrack *outertrack=(AliHLTTPCTrack*)tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
    for(Int_t in =0;in<kNIn;in++){
      if(ismatchedin[in]) continue;
      AliHLTTPCTrack *innertrack=(AliHLTTPCTrack*)tracksin->GetCheckedTrack(in);
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
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCTrackMerger::Merge","Result")
  <<AliHLTTPCLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  delete[] ismatchedin;
  delete[] ismatchedout;
  return nmerged;
}


