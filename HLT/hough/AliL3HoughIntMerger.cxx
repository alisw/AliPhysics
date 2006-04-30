// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughIntMerger.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackSegmentData.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// Hough Inter merger
// Merging of multiple reconstructed tracks

ClassImp(AliL3HoughIntMerger)

AliL3HoughIntMerger::AliL3HoughIntMerger()
{
  //Default constructor
  InitMerger(1,"AliL3HoughTrack");
  fRowMax = fRowMin = 0;
  SetParameters(0.001,0.05,10);
  Is2Global(kFALSE);
}


AliL3HoughIntMerger::~AliL3HoughIntMerger()
{
  //Destructor
  
}

void AliL3HoughIntMerger::SetParameters(Double_t maxkappa, Double_t maxphi0, Double_t maxtgl)
{
  //Set merger params
  fMaxKappa = maxkappa;
  fMaxPhi0 = maxphi0;
  fMaxTgl = maxtgl;
}

void AliL3HoughIntMerger::FillTracks(AliL3TrackArray *tracks)
{
  //Fills tracks into merger
  if(tracks->GetNTracks()==0)
    LOG(AliL3Log::kWarning,"AliL3HoughIntMerger::FillTracks","Track Array")
      <<"Adding empty track array"<<ENDLOG;
  
  GetInTracks(0)->AddTracks(tracks,kFALSE);//Copy tracks
  printf("Filling %d tracks to intermerger\n",tracks->GetNTracks());
}

Bool_t AliL3HoughIntMerger::IsTrack(AliL3Track *innertrack,AliL3Track *outertrack)
{
  //Check if the tracks can be merged, called by the track merger
  
  AliL3HoughTrack *tr1 = (AliL3HoughTrack*)innertrack;
  AliL3HoughTrack *tr2 = (AliL3HoughTrack*)outertrack;
  
  if(abs(tr1->GetEtaIndex() - tr2->GetEtaIndex()) > 1) return kFALSE;
  if(tr1->GetCharge()!=tr2->GetCharge()) return kFALSE;
  if(fabs(tr1->GetKappa()-tr2->GetKappa())   >fMaxKappa) return kFALSE;
  if(fabs(tr1->GetPhi0()-tr2->GetPhi0()) > fMaxPhi0) return kFALSE;

  return kTRUE;//Tracks could be merged
}

AliL3Track *AliL3HoughIntMerger::MultiMerge(AliL3TrackArray *mergedtrack,AliL3Track **tracks, Int_t /*ntrack*/)
{
  //Called by the track merger

  AliL3HoughTrack *newtrack = (AliL3HoughTrack*)mergedtrack->NextTrack();
  AliL3HoughTrack **trs = (AliL3HoughTrack**)tracks;
      
  AliL3HoughTrack *tpt=trs[0];//this is the "best" track
  //AliL3HoughTrack *tpl=trs[ntrack-1];
  newtrack->Set(tpt);
  return (AliL3Track*)newtrack;
}



void AliL3HoughIntMerger::MMerge()
{
  //Track merging??
  GetInTracks(0)->QSort();
  while(Merge());
  GetOutTracks()->AddTracks(GetInTracks(0));
}

Int_t AliL3HoughIntMerger::Merge()
{
  //Track merging??  
  AliL3TrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliL3Track *tr[2];

  for(Int_t out=0;out<kNIn;out++)
    {
      AliL3HoughTrack *outertrack=(AliL3HoughTrack*)tracks->GetCheckedTrack(out);
      if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++)
	{
	  if(in==out) continue;
	  AliL3HoughTrack *innertrack=(AliL3HoughTrack*)tracks->GetCheckedTrack(in);
	  if(!innertrack) continue;
	  if(IsTrack(innertrack,outertrack))
	    {
	      tr[0]=innertrack;
	      tr[1]=outertrack;
	      SortTracks(tr,2);
	      Print(tr);
	      MultiMerge(tracks,tr,2);
	      tracks->Remove(out);
	      tracks->Remove(in);
	      break;
	    }
	} 
    }
  Int_t nmerged = tracks->GetNTracks()-kNIn; 
  LOG(AliL3Log::kInformational,"AliL3HoughIntMerger::Merge","Result")
    <<AliL3Log::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  
  //add in tracks
  //GetOutTracks()->AddTracks(GetInTracks(0)); 
  
  return nmerged;
}

void AliL3HoughIntMerger::Print(AliL3Track **tracks)
{
  //Prints merger results
  AliL3HoughTrack *tr1 = (AliL3HoughTrack*)tracks[0];
  AliL3HoughTrack *tr2 = (AliL3HoughTrack*)tracks[1];
  Double_t kappadiff = fabs(tr1->GetKappa()-tr2->GetKappa());
  Double_t phi0diff = fabs(tr1->GetPhi0()-tr2->GetPhi0());
  cout << "---------Difference in intermerged tracks---------"<<endl;
  cout << "Kappa: "<<kappadiff<<" Phi0 : "<<phi0diff<<endl;
  
}

void AliL3HoughIntMerger::SortTracks(AliL3Track **trs, Int_t ntrack) const
{
  //Sort the tracks according to their weight

  AliL3HoughTrack **tracks = (AliL3HoughTrack**)trs;
  AliL3HoughTrack **tmp = new  AliL3HoughTrack*[ntrack];
  for(Int_t i=0;i<ntrack;i++) tmp[i] = (AliL3HoughTrack*)tracks[i];
  Int_t *t = new Int_t[ntrack];
  for(Int_t i=0;i<ntrack;i++) t[i]=-1;
  
  for(Int_t j=0;j<ntrack;j++)
    {
      Double_t maxw=0;
      Int_t    maxi=0;
      for(Int_t i=0;i<ntrack;i++)
	{
	  if(!tracks[i]) continue;
	  if(tracks[i]->GetWeight() > maxw)
	    {
	      maxw=tracks[i]->GetWeight();
	      maxi=i;
	    }     
	}
      t[j]=maxi;  
      tracks[maxi]=0;
    }
  for(Int_t i=0;i<ntrack;i++) tracks[i] = tmp[t[i]];
  delete[] t;
  delete[] tmp;
}
