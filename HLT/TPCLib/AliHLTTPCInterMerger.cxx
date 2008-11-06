// @(#) $Id$
// Original: AliHLTInterMerger.cxx,v 1.8 2005/06/14 10:55:21 cvetan

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Uli Frankenfeld, maintained by                          *
 *                  Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCInterMerger.cxx
    @author Uli Frankenfeld, maintained by Matthias Richter
    @date   
    @brief  The HLT TPC track segment merger
*/

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCInterMerger.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackSegmentData.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackArray.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTTPCInterMerger)

AliHLTTPCInterMerger::AliHLTTPCInterMerger()
  :
  fPatch(0),
  fRowMin(0),
  fRowMax(0)
{
  //Default constructor
  InitMerger(1);
  Is2Global(kFALSE);
//  SetParameter(2,2,.3,.3,.3);
  SetParameter(1,0.5,0.0005,0.05,0.1);
}


AliHLTTPCInterMerger::~AliHLTTPCInterMerger(){
  //Destructor
  
}

void AliHLTTPCInterMerger::SlowMerge(){
  Int_t nrow= fRowMax-fRowMin+1;
  void *ntuple=GetNtuple();
  AliHLTTPCTrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliHLTTPCTrack *tr[2];
  Bool_t merge = kTRUE;
  for(Int_t in=0;in<kNIn;in++)
    tracks->GetCheckedTrack(in)->CalculateHelix();
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNIn;out++){
    AliHLTTPCTrack *outertrack=tracks->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        if(in==out) continue;
        AliHLTTPCTrack *innertrack=tracks->GetCheckedTrack(in);
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
      AliHLTTPCTrack *outertrack=tracks->GetTrack(outmin);
      AliHLTTPCTrack *innertrack=tracks->GetTrack(inmin);
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
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCInterMerger::SlowMerge","Result")
  <<AliHLTTPCLog::kDec<<"Merged Tracks: "<<tracks->GetNTracks()-kNIn<<ENDLOG;

  char name[256];
  sprintf(name,"ntuple_i_%d.root",fPatch);
  WriteNtuple(name,ntuple);
}

void AliHLTTPCInterMerger::MMerge(){
  while(Merge()) {/*empty body*/};
  GetOutTracks()->AddTracks(GetInTracks(0));
}

Int_t AliHLTTPCInterMerger::Merge(){
  Int_t nrow= fRowMax-fRowMin+1;
  Double_t xval =AliHLTTPCTransform::Row2X((fRowMax+fRowMin)/2);
  AliHLTTPCTrackArray * tracks = GetInTracks(0);
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
  LOG(AliHLTTPCLog::kInformational,"AliHLTTPCInterMerger::Merge","Result")
  <<AliHLTTPCLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  //add in tracks
//  GetOutTracks()->AddTracks(GetInTracks(0)); 

  return nmerged;
}


