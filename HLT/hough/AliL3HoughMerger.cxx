// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughMerger.h"
#include "AliL3HoughTransformer.h"

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughMerger
//
// Patch merging class for Hough tracklets

ClassImp(AliL3HoughMerger)

  
AliL3HoughMerger::AliL3HoughMerger()
{
  //Default constructor
}


AliL3HoughMerger::AliL3HoughMerger(Int_t nsubsectors) 
{
  //Constructor
  InitMerger(nsubsectors,"AliL3HoughTrack");
  Is2Global(kFALSE);
  SetParameters(0.001,0.1,0.05);
}


AliL3HoughMerger::~AliL3HoughMerger()
{
  //dtor 
}

void AliL3HoughMerger::FillTracks(AliL3TrackArray *tracks,Int_t patch)
{
  //Fills tracks into merger
  if(tracks->GetNTracks()==0)
    LOG(AliL3Log::kWarning,"AliL3HoughMerger::FillTracks","Track Array")
      <<"Adding empty track array"<<ENDLOG;
  
  GetInTracks(patch)->AddTracks(tracks,kFALSE);//Copy tracks
  printf("Filling %d tracks to merger\n",tracks->GetNTracks());
}

void AliL3HoughMerger::SetParameters(Double_t maxkappa,Double_t maxpsi,Double_t maxphi0)
{
  //Set merger params
  fMaxKappa = maxkappa;
  fMaxPsi = maxpsi;
  fMaxPhi0 = maxphi0;
}

Bool_t AliL3HoughMerger::IsTrack(AliL3Track *innertrack,AliL3Track *outertrack)
{
  //Check if the tracks can be merged, called by the track merger
  
  AliL3HoughTrack *tr1 = (AliL3HoughTrack*)innertrack;
  AliL3HoughTrack *tr2 = (AliL3HoughTrack*)outertrack;
  
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

void AliL3HoughMerger::AddTrack(AliL3TrackArray *mergedtrack,AliL3Track *track)
{
  //Adds track to an already merged one
  AliL3Track *t[1];
  t[0] = track;
  MultiMerge(mergedtrack,t,1);
}

AliL3Track *AliL3HoughMerger::MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t ntrack)
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

void AliL3HoughMerger::MergePatches(Bool_t slow)
{
  //Merge tracks from across the patches.
  
  fSlow = slow;
  AliL3TrackArray *tracks;
  AliL3HoughTrack *track;
  for(Int_t i=0; i<GetNIn(); i++)
    {
      tracks = GetInTracks(i);
      for(Int_t j=0; j<tracks->GetNTracks(); j++)
	{
	  track = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  if(!track) continue;
	  track->UpdateToFirstRow();
	}
    }
  Merge();
  
}

void AliL3HoughMerger::Merge()
{
  //Merging of tracks
  Double_t edge0 = AliL3Transform::Pi()/18;
  //Double_t edge1 = 2*PI - edge0;
  AliL3TrackArray *ttt = GetOutTracks();
  
  Int_t subsec = GetNIn() - 2; 
  for(Int_t i=subsec;i>=0;i--){
    AliL3TrackArray *tout = GetOutTracks();
    if(i==subsec) tout = GetInTracks(subsec+1);
    AliL3TrackArray *tin = GetInTracks(i);
    Double_t xval = AliL3Transform::Row2X(AliL3Transform::GetLastRow(i));
    //Double_t xmax = AliL3Transform::Row2X(AliL3Transform::GetLastRow(i+1));
    Double_t ymax = xval*tan(edge0);
    for(Int_t out=0;out<tout->GetNTracks();out++){
      AliL3Track *outtrack=tout->GetCheckedTrack(out);
      if(!outtrack) continue;
      //outtrack->CalculateHelix();
      outtrack->CalculatePoint(xval);
      if(outtrack->IsPoint()&&fabs(outtrack->GetPointY())>ymax){
	tout->Remove(out);
      }
    }
    //    tout->Compress();
    for(Int_t in=0;in<tin->GetNTracks();in++){
      AliL3Track *intrack=(AliL3Track*)tin->GetTrack(in);
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
    */
  } // end subsector loop
  LOG(AliL3Log::kInformational,"AliL3HoughMerger::Merge","Result")
    <<AliL3Log::kDec<<"Total Merged Tracks: "<<GetOutTracks()->GetNPresent()
    <<ENDLOG;
}

Int_t AliL3HoughMerger::Merge(AliL3TrackArray* mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout)
{
  //Merging of tracks  
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
  LOG(AliL3Log::kInformational,"AliL3HoughMerger::Merge","Result")
    <<AliL3Log::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  delete[] ismatchedin;
  delete[] ismatchedout;
  return nmerged;
}

void AliL3HoughMerger::SlowMerge(AliL3TrackArray *mergedtrack,AliL3TrackArray *tracksin,AliL3TrackArray *tracksout,Double_t xval)
{
  //Slow merging of tracks??
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
  LOG(AliL3Log::kInformational,"AliL3HoughMerger::SlowMerge","Result")
    <<AliL3Log::kDec<<"Merged Tracks: "
    <<mergedtrack->GetNTracks()-kNMerged<<ENDLOG;
  char name[256] = "ntuple_t.root";
  for(Int_t i=0;i<GetNIn();i++)
    if(tracksin==GetInTracks(i))
      sprintf(name,"ntuple_t_%d.root",i);
  WriteNtuple(name,ntuple);
}

void AliL3HoughMerger::Print(AliL3Track **tracks)
{
  //Print merging results
  AliL3HoughTrack *tr1 = (AliL3HoughTrack*)tracks[0];
  AliL3HoughTrack *tr2 = (AliL3HoughTrack*)tracks[1];
  Double_t kappadiff = fabs(tr1->GetKappa()-tr2->GetKappa());
  Double_t phi0diff = fabs(tr1->GetPhi0()-tr2->GetPhi0());
  cout << "---------Difference in merged tracks---------"<<endl;
  cout << "Kappa: "<<kappadiff<<" Phi0 : "<<phi0diff<<endl;
  
}
