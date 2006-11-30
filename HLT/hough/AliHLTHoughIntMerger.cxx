// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTHoughIntMerger.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTTransform.h"
#include "AliHLTTrackArray.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// Hough Inter merger
// Merging of multiple reconstructed tracks

ClassImp(AliHLTHoughIntMerger)

AliHLTHoughIntMerger::AliHLTHoughIntMerger()
{
  //Default constructor
  InitMerger(1,"AliHLTHoughTrack");
  fRowMax = fRowMin = 0;
  SetParameters(0.001,0.05,10);
  Is2Global(kFALSE);
}


AliHLTHoughIntMerger::~AliHLTHoughIntMerger()
{
  //Destructor
  
}

void AliHLTHoughIntMerger::SetParameters(Double_t maxkappa, Double_t maxphi0, Double_t maxtgl)
{
  //Set merger params
  fMaxKappa = maxkappa;
  fMaxPhi0 = maxphi0;
  fMaxTgl = maxtgl;
}

void AliHLTHoughIntMerger::FillTracks(AliHLTTrackArray *tracks)
{
  //Fills tracks into merger
  if(tracks->GetNTracks()==0)
    LOG(AliHLTLog::kWarning,"AliHLTHoughIntMerger::FillTracks","Track Array")
      <<"Adding empty track array"<<ENDLOG;
  
  GetInTracks(0)->AddTracks(tracks,kFALSE);//Copy tracks
  printf("Filling %d tracks to intermerger\n",tracks->GetNTracks());
}

Bool_t AliHLTHoughIntMerger::IsTrack(AliHLTTrack *innertrack,AliHLTTrack *outertrack)
{
  //Check if the tracks can be merged, called by the track merger
  
  AliHLTHoughTrack *tr1 = (AliHLTHoughTrack*)innertrack;
  AliHLTHoughTrack *tr2 = (AliHLTHoughTrack*)outertrack;
  
  if(abs(tr1->GetEtaIndex() - tr2->GetEtaIndex()) > 1) return kFALSE;
  if(tr1->GetCharge()!=tr2->GetCharge()) return kFALSE;
  if(fabs(tr1->GetKappa()-tr2->GetKappa())   >fMaxKappa) return kFALSE;
  if(fabs(tr1->GetPhi0()-tr2->GetPhi0()) > fMaxPhi0) return kFALSE;

  return kTRUE;//Tracks could be merged
}

AliHLTTrack *AliHLTHoughIntMerger::MultiMerge(AliHLTTrackArray *mergedtrack,AliHLTTrack **tracks, Int_t /*ntrack*/)
{
  //Called by the track merger

  AliHLTHoughTrack *newtrack = (AliHLTHoughTrack*)mergedtrack->NextTrack();
  AliHLTHoughTrack **trs = (AliHLTHoughTrack**)tracks;
      
  AliHLTHoughTrack *tpt=trs[0];//this is the "best" track
  //AliHLTHoughTrack *tpl=trs[ntrack-1];
  newtrack->Set(tpt);
  return (AliHLTTrack*)newtrack;
}



void AliHLTHoughIntMerger::MMerge()
{
  //Track merging??
  GetInTracks(0)->QSort();
  while(Merge());
  GetOutTracks()->AddTracks(GetInTracks(0));
}

Int_t AliHLTHoughIntMerger::Merge()
{
  //Track merging??  
  AliHLTTrackArray * tracks = GetInTracks(0);
  const Int_t  kNIn =tracks->GetNTracks();
  AliHLTTrack *tr[2];

  for(Int_t out=0;out<kNIn;out++)
    {
      AliHLTHoughTrack *outertrack=(AliHLTHoughTrack*)tracks->GetCheckedTrack(out);
      if(!outertrack) continue;
      for(Int_t in=0;in<kNIn;in++)
	{
	  if(in==out) continue;
	  AliHLTHoughTrack *innertrack=(AliHLTHoughTrack*)tracks->GetCheckedTrack(in);
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
  LOG(AliHLTLog::kInformational,"AliHLTHoughIntMerger::Merge","Result")
    <<AliHLTLog::kDec<<"Merged Tracks: "<<nmerged<<ENDLOG;
  
  //add in tracks
  //GetOutTracks()->AddTracks(GetInTracks(0)); 
  
  return nmerged;
}

void AliHLTHoughIntMerger::Print(AliHLTTrack **tracks)
{
  //Prints merger results
  AliHLTHoughTrack *tr1 = (AliHLTHoughTrack*)tracks[0];
  AliHLTHoughTrack *tr2 = (AliHLTHoughTrack*)tracks[1];
  Double_t kappadiff = fabs(tr1->GetKappa()-tr2->GetKappa());
  Double_t phi0diff = fabs(tr1->GetPhi0()-tr2->GetPhi0());
  cout << "---------Difference in intermerged tracks---------"<<endl;
  cout << "Kappa: "<<kappadiff<<" Phi0 : "<<phi0diff<<endl;
  
}

void AliHLTHoughIntMerger::SortTracks(AliHLTTrack **trs, Int_t ntrack) const
{
  //Sort the tracks according to their weight

  AliHLTHoughTrack **tracks = (AliHLTHoughTrack**)trs;
  AliHLTHoughTrack **tmp = new  AliHLTHoughTrack*[ntrack];
  for(Int_t i=0;i<ntrack;i++) tmp[i] = (AliHLTHoughTrack*)tracks[i];
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
