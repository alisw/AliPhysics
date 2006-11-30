// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTTrackArray.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTHoughMerger.h"
#include "AliHLTHoughTransformer.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// AliHLTHoughMerger
//
// Patch merging class for Hough tracklets

ClassImp(AliHLTHoughMerger)

  
AliHLTHoughMerger::AliHLTHoughMerger()
{
  //Default constructor
}


AliHLTHoughMerger::AliHLTHoughMerger(Int_t nsubsectors) 
{
  //Constructor
  InitMerger(nsubsectors,"AliHLTHoughTrack");
  Is2Global(kFALSE);
  SetParameters(0.001,0.1,0.05);
}


AliHLTHoughMerger::~AliHLTHoughMerger()
{
  //dtor 
}

void AliHLTHoughMerger::FillTracks(AliHLTTrackArray *tracks,Int_t patch)
{
  //Fills tracks into merger
  if(tracks->GetNTracks()==0)
    LOG(AliHLTLog::kWarning,"AliHLTHoughMerger::FillTracks","Track Array")
      <<"Adding empty track array"<<ENDLOG;
  
  GetInTracks(patch)->AddTracks(tracks,kFALSE);//Copy tracks
  printf("Filling %d tracks to merger\n",tracks->GetNTracks());
}

void AliHLTHoughMerger::SetParameters(Double_t maxkappa,Double_t maxpsi,Double_t maxphi0)
{
  //Set merger params
  fMaxKappa = maxkappa;
  fMaxPsi = maxpsi;
  fMaxPhi0 = maxphi0;
}

Bool_t AliHLTHoughMerger::IsTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack)
{
  //Check if the tracks can be merged, called by the track merger
  
  AliHLTHoughTrack *tr1 = (AliHLTHoughTrack*)innertrack;
  AliHLTHoughTrack *tr2 = (AliHLTHoughTrack*)outertrack;
  
  if( (!tr1->IsPoint()) || (!tr2->IsPoint()) )  return kFALSE; 
  if(abs(tr1->GetEtaIndex() - tr2->GetEtaIndex()) > 1) return kFALSE;
  if(tr1->GetCharge() != tr2->GetCharge()) return kFALSE;
  if(fabs(tr1->GetPhi0() - tr2->GetPhi0()) > fMaxPhi0) return kFALSE;
  if(fabs(tr1->GetKappa() - tr2->GetKappa()) > fMaxKappa) return kFALSE;
  
  /*
    if( (!tr1->IsPoint()) || (!tr2->IsPoint()) )  return kFALSE; 
    if(fabs(innertrack->GetPointY()-outertrack->GetPointY()) >fMaxY) return kFALSE;
    if(fabs(innertrack->GetPointZ()-outertrack->GetPointZ()) >fMaxZ) return kFALSE;
    if(fabs(innertrack->GetKappa()-outertrack->GetKappa())   >fMaxKappa) return kFALSE;
    if(GetAngle(innertrack->GetPointPsi(),outertrack->GetPointPsi()) >fMaxPsi) return kFALSE;
    if(fabs(innertrack->GetTgl()-outertrack->GetTgl()) >fMaxTgl) return kFALSE;
  */
  
  return kTRUE;//Tracks could be merged
}

void AliHLTHoughMerger::AddTrack(AliHLTTrackArray *mergedtrack,AliHLTTrack *track)
{
  //Adds track to an already merged one
  AliHLTTrack *t[1];
  t[0] = track;
  MultiMerge(mergedtrack,t,1);
}

AliHLTTrack *AliHLTHoughMerger::MultiMerge(AliHLTTrackArray *mergedtrack,AliHLTTrack **tracks, Int_t ntrack)
{
  //Called by the track merger

  AliHLTHoughTrack *newtrack = (AliHLTHoughTrack*)mergedtrack->NextTrack();
  AliHLTHoughTrack **trs = (AliHLTHoughTrack**)tracks;
  Int_t weight=0;

  //Sum up the total weight:
  for(Int_t i=ntrack-1; i>=0; i--)
    weight += trs[i]->GetWeight();
  
  AliHLTHoughTrack *tpt=trs[0];//This is the innermost track
  AliHLTHoughTrack *tpl=trs[ntrack-1];
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
  
  return (AliHLTTrack*)newtrack;

}

void AliHLTHoughMerger::MergePatches(Bool_t slow)
{
  //Merge tracks from across the patches.
  
  fSlow = slow;
  AliHLTTrackArray *tracks;
  AliHLTHoughTrack *track;
  for(Int_t i=0; i<GetNIn(); i++)
    {
      tracks = GetInTracks(i);
      for(Int_t j=0; j<tracks->GetNTracks(); j++)
	{
	  track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(j);
	  if(!track) continue;
	  track->UpdateToFirstRow();
	}
    }
  Merge();
  
}

void AliHLTHoughMerger::Merge()
{
  //Merging of tracks
  Double_t edge0 = AliHLTTransform::Pi()/18;
  //Double_t edge1 = 2*PI - edge0;
  AliHLTTrackArray *ttt = GetOutTracks();
  
  Int_t subsec = GetNIn() - 2; 
  for(Int_t i=subsec;i>=0;i--){
    AliHLTTrackArray *tout = GetOutTracks();
    if(i==subsec) tout = GetInTracks(subsec+1);
    AliHLTTrackArray *tin = GetInTracks(i);
    Double_t xval = AliHLTTransform::Row2X(AliHLTTransform::GetLastRow(i));
    //Double_t xmax = AliHLTTransform::Row2X(AliHLTTransform::GetLastRow(i+1));
    Double_t ymax = xval*tan(edge0);
    for(Int_t out=0;out<tout->GetNTracks();out++){
      AliHLTTrack *outtrack=tout->GetCheckedTrack(out);
      if(!outtrack) continue;
      //outtrack->CalculateHelix();
      outtrack->CalculatePoint(xval);
      if(outtrack->IsPoint()&&fabs(outtrack->GetPointY())>ymax){
	tout->Remove(out);
      }
    }
    //    tout->Compress();
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliHLTTrack *intrack=(AliHLTTrack*)tin->GetTrack(in);
      //intrack->CalculateHelix();
      intrack->CalculatePoint(xval);
    }
    tin->QSort();
    tout->QSort();

    if(fSlow) SlowMerge(ttt,tin,tout,xval);
    else Merge(ttt,tin,tout);

    /*
    //Add the tracks that cross the sector boundary:
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliHLTTrack *intrack=(AliHLTTrack*)tin->GetCheckedTrack(in);
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
    */
  } // end subsector loop
  LOG(AliHLTLog::kInformational,"AliHLTHoughMerger::Merge","Result")
    <<AliHLTLog::kDec<<"Total Merged Tracks: "<<GetOutTracks()->GetNPresent()
    <<ENDLOG;
}

Int_t AliHLTHoughMerger::Merge(AliHLTTrackArray* mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout)
{
  //Merging of tracks  
  AliHLTTrack *tracks[2];
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
    AliHLTTrack *outertrack=(AliHLTTrack*)tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
    for(Int_t in =0;in<kNIn;in++){
      if(ismatchedin[in]) continue;
      AliHLTTrack *innertrack=(AliHLTTrack*)tracksin->GetCheckedTrack(in);
      if(!innertrack) continue;
      if(outertrack==innertrack) continue;
      
      if(IsTrack(innertrack,outertrack)) //They can be merged
	{
	  tracks[0]=innertrack; tracks[1]=outertrack; 
	  SortTracks(tracks,2); //Sort the tracks according to minimum x-point
	  //if(tracks[0]->GetLastPointX()<tracks[1]->GetFirstPointX()){
	  MultiMerge(mergedtrack,tracks,2);
	  tracksout->Remove(out);
	  tracksin->Remove(in);
	  ismatchedin[in]=kTRUE;
	  ismatchedout[out]=kTRUE;
	  break;
	  // }
	}
    }
  }
  
  Int_t nmerged = mergedtrack->GetNTracks()-kNMerged;
  LOG(AliHLTLog::kInformational,"AliHLTHoughMerger::Merge","Result")
    <<AliHLTLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  delete[] ismatchedin;
  delete[] ismatchedout;
  return nmerged;
}

void AliHLTHoughMerger::SlowMerge(AliHLTTrackArray *mergedtrack,AliHLTTrackArray *tracksin,AliHLTTrackArray *tracksout,Double_t xval)
{
  //Slow merging of tracks??
  void *ntuple=GetNtuple();
  const Int_t  kNOut=tracksout->GetNTracks();
  const Int_t  kNIn =tracksin->GetNTracks();
  const Int_t  kNMerged =mergedtrack->GetNTracks();
  AliHLTTrack *tracks[2];
  Bool_t merge = kTRUE;
  while(merge){
    Int_t inmin=-1,outmin=-1;
    Double_t min=10;
    for(Int_t out=0;out<kNOut;out++){
    AliHLTTrack *outertrack=tracksout->GetCheckedTrack(out);
    if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++){
        AliHLTTrack *innertrack=tracksin->GetCheckedTrack(in);
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
      AliHLTTrack *outertrack=tracksout->GetTrack(outmin);
      AliHLTTrack *innertrack=tracksin->GetTrack(inmin);
      tracks[0]=innertrack;
      tracks[1]=outertrack;
      SortTracks(tracks,2);
      Print(tracks);
      MultiMerge(mergedtrack,tracks,2);
      outertrack->CalculatePoint(xval);
      innertrack->CalculatePoint(xval);
      FillNtuple(ntuple,innertrack,outertrack);
      tracksout->Remove(outmin);
      tracksin->Remove(inmin);
      //      tracksout->Compress();
      //      tracksin->Compress(); 
    }
    else merge = kFALSE;
  }
  LOG(AliHLTLog::kInformational,"AliHLTHoughMerger::SlowMerge","Result")
    <<AliHLTLog::kDec<<"Merged Tracks: "
    <<mergedtrack->GetNTracks()-kNMerged<<ENDLOG;
  char name[256] = "ntuple_t.root";
  for(Int_t i=0;i<GetNIn();i++)
    if(tracksin==GetInTracks(i))
      sprintf(name,"ntuple_t_%d.root",i);
  WriteNtuple(name,ntuple);
}

void AliHLTHoughMerger::Print(AliHLTTrack **tracks)
{
  //Print merging results
  AliHLTHoughTrack *tr1 = (AliHLTHoughTrack*)tracks[0];
  AliHLTHoughTrack *tr2 = (AliHLTHoughTrack*)tracks[1];
  Double_t kappadiff = fabs(tr1->GetKappa()-tr2->GetKappa());
  Double_t phi0diff = fabs(tr1->GetPhi0()-tr2->GetPhi0());
  cout << "---------Difference in merged tracks---------"<<endl;
  cout << "Kappa: "<<kappadiff<<" Phi0 : "<<phi0diff<<endl;
  
}
