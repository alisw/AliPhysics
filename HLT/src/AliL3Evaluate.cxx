// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPC.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersArray.h>
#include <AliTPCClustersRow.h>
#include <AliTPCParam.h>
#include <AliComplexCluster.h>
#include <AliStack.h>

#if __GNUC__ == 3
#include <fstream>
#include <iosfwd>
#else
#include <fstream.h>
#endif

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3SpacePointData.h"
#include "AliL3Track.h"
#include "AliL3FileHandler.h"
#include "AliL3TrackArray.h"
#include "AliL3Evaluate.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3Evaluate
<pre>
//_____________________________________________________________
// AliL3Evaluate
//
// Evaluation class for tracking. Plots efficiencies etc..
//
</pre>
*/

ClassImp(AliL3Evaluate)

AliL3Evaluate::AliL3Evaluate()
{
  fTracks = 0;
  fNFastPoints = 0;
  fMcIndex = 0;
  fMcId = 0;
  fMinSlice=0;
  fMaxSlice=0;
  fGoodTracks=0;
  fNtupleRes=0;
  fDigitsTree=0;
  fPtRes=0;
  fNGoodTracksPt=0;
  fNFoundTracksPt=0;
  fNFakeTracksPt=0;
  fTrackEffPt=0;
  fFakeTrackEffPt=0;
  fNGoodTracksEta=0;
  fNFoundTracksEta=0;
  fNFakeTracksEta=0;
  fTrackEffEta=0;
  fFakeTrackEffEta=0;
  fMcIndex=0;
  fMcId=0;
  fNtuppel=0;
  fStandardComparison=kTRUE;
}

AliL3Evaluate::AliL3Evaluate(Char_t *datapath,Int_t min_clusters,Int_t minhits,Double_t minpt,Double_t maxpt,Int_t *slice)
{

  if(slice)
    {
      fMinSlice=slice[0];
      fMaxSlice=slice[1];
    }
  else
    {
      fMinSlice=0;
      fMaxSlice=35;
    }
  sprintf(fPath,"%s",datapath);
  fMaxFalseClusters = 0.1;
  fGoodFound = 0;
  fGoodGen = 0;
  fMinPointsOnTrack = min_clusters;
  fMinHitsFromParticle = minhits;
  fMinGoodPt = minpt;
  fMaxGoodPt = maxpt;
  memset(fClusters,0,36*6*sizeof(AliL3SpacePointData*));
  fTracks=0;
  fGoodTracks=0;
  fNtupleRes=0;
  fDigitsTree=0;
  fPtRes=0;
  fNGoodTracksPt=0;
  fNFoundTracksPt=0;
  fNFakeTracksPt=0;
  fTrackEffPt=0;
  fFakeTrackEffPt=0;
  fNGoodTracksEta=0;
  fNFoundTracksEta=0;
  fNFakeTracksEta=0;
  fTrackEffEta=0;
  fFakeTrackEffEta=0;
  fMcIndex=0;
  fMcId=0;
  fNtuppel=0;
  fStandardComparison=kTRUE;
}

void AliL3Evaluate::LoadData(Int_t event,Bool_t sp)
{
  Char_t fname[1024];
  AliL3FileHandler *clusterfile[36][6];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliL3Transform::GetNPatches(); p++)
	{
	  Int_t patch;
	  if(sp==kTRUE)
	    patch=-1;
	  else
	    patch=p;
	  
	  if(fClusters[s][p])
	    delete fClusters[s][p];
	  fClusters[s][p] = 0;
	  clusterfile[s][p] = new AliL3FileHandler();
	  if(event<0)
	    sprintf(fname,"%s/points_%d_%d.raw",fPath,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",fPath,event,s,patch);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
              delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliL3SpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  if(sp==kTRUE)
	    break;
	}
    }
   
  sprintf(fname,"%s/tracks_%d.raw",fPath,event);
  AliL3FileHandler *tfile = new AliL3FileHandler();
  if(!tfile->SetBinaryInput(fname)){
    LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
      <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
    return;
  }
  if(fTracks)
    delete fTracks;
  fTracks = new AliL3TrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
  fTracks->QSort();
}


AliL3Evaluate::~AliL3Evaluate()
{
  if(fGoodTracks) delete fGoodTracks;
  if(fDigitsTree) fDigitsTree->Delete();
  if(fTracks) delete fTracks;
  if(fPtRes) delete fPtRes;
  if(fNGoodTracksPt) delete fNGoodTracksPt;
  if(fNFoundTracksPt) delete fNFoundTracksPt;
  if(fNFakeTracksPt) delete fNFakeTracksPt;
  if(fTrackEffPt) delete fTrackEffPt;
  if(fFakeTrackEffPt) delete fFakeTrackEffPt;
  if(fNGoodTracksEta) delete fNGoodTracksEta;
  if(fNFoundTracksEta) delete fNFoundTracksEta;
  if(fNFakeTracksEta) delete fNFakeTracksEta;
  if(fTrackEffEta) delete fTrackEffEta;
  if(fFakeTrackEffEta) delete fFakeTrackEffEta;
  if(fMcIndex) delete [] fMcIndex;
  if(fMcId)    delete [] fMcId;
  if(fNtuppel) delete fNtuppel;
  if(fNtupleRes) delete fNtupleRes;
  for(Int_t s=0; s<=35; s++)
    for(Int_t p=0; p<6; p++)
      if(fClusters[s][p]) delete fClusters[s][p];
}

void AliL3Evaluate::AssignPIDs()
{
  fTracks->QSort();
  LOG(AliL3Log::kDebug,"AliL3Evaluate::AssignPIDs","Track Loop")
    <<"Assigning pid to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      if(track->GetNumberOfPoints() < fMinPointsOnTrack)
	track->SetPID(0);
      else {
	Float_t pid = GetTrackPID(track);
	track->SetPID(pid);
      }
    }
}
  
void AliL3Evaluate::AssignIDs()
{
  //Assign MC id to the tracks.
#ifndef do_mc
  cerr<<"AliL3Evaluate::AssignIDs() : You need to compile with the do_mc flag!"<<endl;
  return;
#endif
  fGoodFound=0;
  fTracks->QSort();
  LOG(AliL3Log::kDebug,"AliL3Evaluate::AssignIDs","Track Loop")
    <<"Assigning MC id to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      if(track->GetNumberOfPoints() < fMinPointsOnTrack) break;
      
      fGoodFound++;
      Int_t tID = GetMCTrackLabel(track);
      track->SetMCid(tID);
    }
  //cout<<"Found "<<fGoodFound<<" good tracks "<<endl;
}

Float_t AliL3Evaluate::GetTrackPID(AliL3Track *track)
{
  track->CalculateHelix();
  // Track dEdx
  Int_t nc=track->GetNHits();
  UInt_t *hits = track->GetHitNumbers();
  Float_t sampleDEdx[159];
  for (Int_t iHit = 0; iHit < nc; iHit++) {
    UInt_t hitID = hits[iHit];
    Int_t iSector = (hitID>>25) & 0x7f;
    Int_t patch = (hitID>>22) & 0x7;
    UInt_t position = hitID&0x3fffff;
    AliL3SpacePointData *points = fClusters[iSector][patch];
    if(!points) continue; 
    if(position>=fNcl[iSector][patch]) 
      {
	LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	  <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	continue;
      }
    UChar_t padrow = points[position].fPadRow;
    Float_t pWidth = AliL3Transform::GetPadPitchWidthLow();
    if (padrow>63)
      pWidth = AliL3Transform::GetPadPitchWidthUp(); 
    Float_t corr=1.; if (padrow>63) corr=0.67;
    sampleDEdx[iHit] = points[position].fCharge/pWidth*corr;
    Double_t crossingangle = track->GetCrossingAngle(padrow,iSector);
    Double_t s = sin(crossingangle);
    Double_t t = track->GetTgl();
    sampleDEdx[iHit] *= sqrt((1-s*s)/(1+t*t));
  }

  /* Cook dEdx */
  Int_t i;
  Int_t swap;//stupid sorting
  do {
    swap=0;
    for (i=0; i<nc-1; i++) {
      if (sampleDEdx[i]<=sampleDEdx[i+1]) continue;
      Float_t tmp=sampleDEdx[i];
      sampleDEdx[i]=sampleDEdx[i+1]; sampleDEdx[i+1]=tmp;
      swap++;
    }
  } while (swap);

  Double_t low=0.05; Double_t up=0.7;
  Int_t nl=Int_t(low*nc), nu=Int_t(up*nc);
  Float_t trackDEdx=0;
  for (i=nl; i<=nu; i++) trackDEdx += sampleDEdx[i];
  trackDEdx /= (nu-nl+1);

  //  cout<<" PID: "<<nc<<" "<<nl<<" "<<nu<<" "<<trackDEdx<<" "<<track->GetPt()<<endl;
  return trackDEdx;
}

struct S {Int_t lab; Int_t max;};
Int_t AliL3Evaluate::GetMCTrackLabel(AliL3Track *track){ 
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.
  //Definitions are identical to offline.
  //Fake track means:
  // - more than 10 percent of clusters are assigned incorrectly
  // - more than half of the innermost 10% of clusters were assigned incorrectly.
  
  
#ifdef do_mc
  Int_t num_of_clusters = track->GetNumberOfPoints();
  S *s=new S[num_of_clusters];
  Int_t i;
  for (i=0; i<num_of_clusters; i++) s[i].lab=s[i].max=0;
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;
    
  Int_t lab=123456789;
  for (i=0; i<num_of_clusters; i++) 
    {
      //Tricks to get the clusters belonging to this track:
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      //Get the label of the cluster:
      lab=points[pos].fTrackID[0];

      Int_t j;
      for (j=0; j<num_of_clusters; j++)
        if (s[j].lab==lab || s[j].max==0) break;
      s[j].lab=lab;
      s[j].max++;
    }
  
  Int_t max=0;
  for (i=0; i<num_of_clusters; i++) 
    if (s[i].max>max) {max=s[i].max; lab=s[i].lab;}
  
  if(lab == -1)
    return -1; //If most clusters is -1, this is a noise track.
  if(lab < 0)
    cerr<<"AliL3Evaluate::GetMCTrackLabel : Track label negative :"<<lab<<endl;
  
  delete[] s;
  
  for (i=0; i<num_of_clusters; i++) 
    {
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      if (abs(points[pos].fTrackID[1]) == lab || 
	  abs(points[pos].fTrackID[2]) == lab ) max++;
    }
  
  
  //Check if more than 10% of the clusters were assigned incorrectly:
  if (1.-Float_t(max)/num_of_clusters > fMaxFalseClusters) 
    {
      return -lab;
    }
  else //Check if at least half of the 10% innermost clusters are assigned correctly.
    {
      Int_t tail=Int_t(0.10*num_of_clusters);
      max=0;
      for (i=1; i<=tail; i++) 
	{
	  id = hitnum[num_of_clusters - i];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  
	  AliL3SpacePointData *points = fClusters[slice][patch];
	  if(lab == abs(points[pos].fTrackID[0]) ||
	     lab == abs(points[pos].fTrackID[1]) ||
	     lab == abs(points[pos].fTrackID[2])) max++;
	}
      if (max < Int_t(0.5*tail)) return -lab;
    }

  return lab;
#else //If we are running with mc_ids or not
  return 0;
#endif

}

void AliL3Evaluate::GetFastClusterIDs(Char_t *path)
{
  //Get the MC id of space points in case of using the fast simulator. 
  char fname[256];
  sprintf(fname,"%s/point_mc.dat",path);
  FILE *infile = fopen(fname,"r");
  if(!infile) return;
  Int_t hitid,hitmc,i;
  
  for(i=0; ; i++)
    if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
  rewind(infile);
  fNFastPoints = i;
  fMcId = new Int_t[fNFastPoints];
  fMcIndex = new UInt_t[fNFastPoints];
  
  for(i=0; i<fNFastPoints; i++)
    {
      if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
      fMcId[i] = hitmc;
      fMcIndex[i] = hitid;
    }
  fclose(infile);
}

void AliL3Evaluate::CreateHistos(Int_t nbin,Float_t xlow,Float_t xup)
{
  //Create the histograms 
  
  LOG(AliL3Log::kInformational,"AliL3Evaluate::CreateHistos","Allocating")
    <<"Creating histograms..."<<ENDLOG;
  
  fNtuppel = new TNtuple("fNtuppel","Pt resolution","pt_gen:pt_found:nHits");
  fNtuppel->SetDirectory(0);
  fPtRes = new TH1F("fPtRes","Relative Pt resolution",30,-10.,10.); 
  fNGoodTracksPt = new TH1F("fNGoodTracksPt","Good tracks vs pt",nbin,xlow,xup);    
  fNFoundTracksPt = new TH1F("fNFoundTracksPt","Found tracks vs pt",nbin,xlow,xup);
  fNFakeTracksPt = new TH1F("fNFakeTracksPt","Fake tracks vs pt",nbin,xlow,xup);
  fTrackEffPt = new TH1F("fTrackEffPt","Tracking efficiency vs pt",nbin,xlow,xup);
  fFakeTrackEffPt = new TH1F("fFakeTrackEffPt","Efficiency for fake tracks vs pt",nbin,xlow,xup);
  
  fNGoodTracksEta = new TH1F("fNGoodTracksEta","Good tracks vs eta",20,-50,50);
  fNFoundTracksEta = new TH1F("fNFoundTracksEta","Found tracks vs eta",20,-50,50);
  fNFakeTracksEta = new TH1F("fNFakeTracksEta","Fake tracks vs eta",20,-50,50);
  fTrackEffEta = new TH1F("fTrackEffEta","Tracking efficienct vs eta",20,-50,50);
  fFakeTrackEffEta = new TH1F("fFakeTrackEffEta","Efficiency for fake tracks vs eta",20,-50,50);

}

void AliL3Evaluate::GetGoodParticles(Char_t *path,Int_t event,Int_t *padrowrange)
{
  //Read the good particles from file. This file should already have been
  //generated by macro AliTPCComparison.C.
  
  Char_t filename[1024];
  if(event<0 && !padrowrange)
    sprintf(filename,"%s/good_tracks_tpc",path);
  else if(event>=0 && !padrowrange)
    sprintf(filename,"%s/good_tracks_tpc_%d",path,event);
  else
    sprintf(filename,"%s/good_tracks_tpc_%d_%d_%d",path,event,padrowrange[0],padrowrange[1]);
  ifstream in(filename);
  if(!in)
    {
      cerr<<"AliL3Evaluate::GetGoodParticles : Problems opening file :"<<filename<<endl;
      return;
    }
  Int_t MaxTracks=20000;
  if(fGoodTracks)
    delete [] fGoodTracks;
  fGoodGen=0;
  fGoodTracks = new GoodTrack[MaxTracks];
  
  if(fStandardComparison){
    while (in>>fGoodTracks[fGoodGen].label>>fGoodTracks[fGoodGen].code>>
	   fGoodTracks[fGoodGen].px>>fGoodTracks[fGoodGen].py>>fGoodTracks[fGoodGen].pz>>
	   fGoodTracks[fGoodGen].x>>fGoodTracks[fGoodGen].y>>fGoodTracks[fGoodGen].z)
      {
        fGoodTracks[fGoodGen].nhits=-1;
        fGoodTracks[fGoodGen].sector=-1; 
	fGoodGen++;
	if (fGoodGen==MaxTracks) 
	  {
	    cerr<<"AliL3Evaluate::GetGoodParticles : Too many good tracks !\n";
	    break;
	  }
      }
  } else {
    while (in>>fGoodTracks[fGoodGen].label>>fGoodTracks[fGoodGen].code>>
	   fGoodTracks[fGoodGen].px>>fGoodTracks[fGoodGen].py>>fGoodTracks[fGoodGen].pz>>
	   fGoodTracks[fGoodGen].x>>fGoodTracks[fGoodGen].y >>fGoodTracks[fGoodGen].z>>fGoodTracks[fGoodGen].nhits>>fGoodTracks[fGoodGen].sector) 
      {
	fGoodGen++;
	if (fGoodGen==MaxTracks) 
	  {
	    cerr<<"AliL3Evaluate::GetGoodParticles : Too many good tracks !\n";
	    break;
	  }
      }
  }
}

//has to be modified for fakes.

void AliL3Evaluate::FillEffHistos()
{  
  if(!fGoodTracks)
    {
      cerr<<"AliL3Evaluate::FillEffHistos : No good tracks"<<endl;
      return;
    }
  //cout<<"Comparing "<<fGoodGen<<" good tracks ..."<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      //cout<<"Checking particle "<<i<<endl;
      if(!fStandardComparison) 
	if(fGoodTracks[i].nhits < fMinHitsFromParticle) continue;
      Float_t ptg = TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py);
      if(ptg < fMinGoodPt || ptg > fMaxGoodPt) continue;
      Float_t pzg=fGoodTracks[i].pz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();
      
      //If we are only considering tracks on one side of the TPC:
      if(fMaxSlice <= 17)
	if(dipangle < 0)
	  continue;

      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(dipangle);
      Int_t found = 0;
      
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliL3Track *track = fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  Int_t nHits = track->GetNumberOfPoints();
	  if(nHits < fMinPointsOnTrack) break;
	  Int_t tracklabel;
	  tracklabel = track->GetMCid();
	  
	  if(TMath::Abs(tracklabel) != fGoodTracks[i].label) continue;
	  found=1;
	  Float_t pt=track->GetPt();
	  if(tracklabel == fGoodTracks[i].label) 
	    {
	      fNFoundTracksPt->Fill(ptg); 
	      fNFoundTracksEta->Fill(dipangle);
	      fNtuppel->Fill(ptg,pt,nHits);
	      fPtRes->Fill((pt-ptg)/ptg*100.);
	    }
	  else 
	    {
	      fNFakeTracksPt->Fill(ptg); 
	      fNFakeTracksEta->Fill(dipangle);
	    }
	  //fPtRes->Fill((pt-ptg)/ptg*100.);
	  //fNtuppel->Fill(ptg,pt,nHits);
	  break;
	  
	}
      //if(!found)
      //cout<<"Track "<<fGoodTracks[i].label<<" was not found"<<endl;
    }
}

void AliL3Evaluate::FillEffHistosNAIVE()
{  
  //Fill the efficiency histograms.
  
  cout<<endl<<"Note: Doing NAIVE evaluation "<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      if(!fStandardComparison) 
	if(fGoodTracks[i].nhits < fMinHitsFromParticle) continue;
      Double_t ptg=TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py);
      if(ptg < fMinGoodPt || ptg > fMaxGoodPt) continue;
      Double_t pzg=fGoodTracks[i].pz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();
      //printf("filling particle with pt %f and dipangle %f\n",ptg,dipangle);
      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(dipangle);
      
    }
  
  for(Int_t k=0; k<fTracks->GetNTracks(); k++)
    {
      AliL3Track *track = fTracks->GetCheckedTrack(k);
      if(!track) continue;
      Int_t nHits = track->GetNumberOfPoints();
      if(nHits < fMinPointsOnTrack) break;
      if(track->GetPt()<fMinGoodPt || track->GetPt() > fMaxGoodPt) continue;
      if(fabs(track->GetPseudoRapidity())>0.9) continue;

      fNFoundTracksPt->Fill(track->GetPt()); fNFoundTracksEta->Fill(track->GetPseudoRapidity());
      //Float_t pt=track->GetPt();
      //fPtRes->Fill((pt-ptg)/ptg*100.);
      //fNtuppel->Fill(ptg,pt,nHits);
            
    }
}

void AliL3Evaluate::CalcEffHistos()
{  

  Stat_t ngood=fNGoodTracksPt->GetEntries();
  Stat_t nfound=fNFoundTracksPt->GetEntries();
  Stat_t nfake=fNFakeTracksPt->GetEntries();

  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"There was "<<ngood<<" generated good tracks"<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Found "<<nfound<<" tracks"<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Integral efficiency is about "<<nfound/ngood*100<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Fake tracks relative is about "<<nfake/ngood*100<<ENDLOG;
  //LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
  //<<AliL3Log::kDec<<"Naive efficiency "<<(Double_t)fGoodFound/(Double_t)fGoodGen<<ENDLOG;

  fNFoundTracksPt->Sumw2(); fNGoodTracksPt->Sumw2();
  fTrackEffPt->Divide(fNFoundTracksPt,fNGoodTracksPt,1,1.,"b");
  fFakeTrackEffPt->Divide(fNFakeTracksPt,fNGoodTracksPt,1,1.,"b");
  fTrackEffPt->SetMaximum(1.4);
  fTrackEffPt->SetXTitle("P_{T} [GeV]");
  fTrackEffPt->SetLineWidth(2);
  fFakeTrackEffPt->SetFillStyle(3013);
  fTrackEffPt->SetLineColor(4);
  fFakeTrackEffPt->SetFillColor(2);

  fNFoundTracksEta->Sumw2(); fNGoodTracksEta->Sumw2();
  fTrackEffEta->Divide(fNFoundTracksEta,fNGoodTracksEta,1,1.,"b");
  fFakeTrackEffEta->Divide(fNFakeTracksEta,fNGoodTracksEta,1,1.,"b");
  fTrackEffEta->SetMaximum(1.4);
  fTrackEffEta->SetXTitle("#lambda [degrees]");
  fTrackEffEta->SetLineWidth(2);
  fFakeTrackEffEta->SetFillStyle(3013);
  fTrackEffEta->SetLineColor(4);
  fFakeTrackEffEta->SetFillColor(2);
       
}

void AliL3Evaluate::Write2File(Char_t *outputfile)
{
  //Write histograms to file:
  
  TFile *of = TFile::Open(outputfile,"RECREATE");
  if(!of->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::Write2File","File Open")
	<<"Problems opening rootfile"<<ENDLOG;
      return;
    }
  
  of->cd();
  fNtuppel->Write();
  fPtRes->Write();
  fNGoodTracksPt->Write();
  fNFoundTracksPt->Write();
  fNFakeTracksPt->Write();
  fTrackEffPt->Write();
  fFakeTrackEffPt->Write();
  fNGoodTracksEta->Write();
  fNFoundTracksEta->Write();
  fNFakeTracksEta->Write();
  fTrackEffEta->Write();
  fFakeTrackEffEta->Write();
  
  of->Close();

}

TNtuple *AliL3Evaluate::GetNtuple()
{
  if(!fNtupleRes)
    {
      fNtupleRes = new TNtuple("ntuppel","Residuals","residual_trans:residual_long:zHit:pt:dipangle:beta:padrow:nHits");
      fNtupleRes->SetDirectory(0);//Bug in older version of root.
    }
  return fNtupleRes;
}

void AliL3Evaluate::CalculateResiduals()
{

  TNtuple *ntuppel = GetNtuple();
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits() < fMinPointsOnTrack) break;
      
      track->CalculateHelix();
      UInt_t *hitnum = track->GetHitNumbers();
      UInt_t id;
      
      Float_t xyz[3];
      Int_t padrow;
      for(Int_t j=0; j<track->GetNumberOfPoints()-1; j++)
	{
	  id = hitnum[j];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      

	  //if(slice<18) continue;
	  
	  AliL3SpacePointData *points = fClusters[slice][patch];
	  
	  if(!points) 
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluate::CalculateResiduals","Clusterarray")
		<<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	      continue;
	    }
	  if(pos>=fNcl[slice][patch]) 
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluate::CalculateResiduals","Clusterarray")
		<<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	      continue;
	    }
	  
	  xyz[0] = points[pos].fX;
	  xyz[1] = points[pos].fY;
	  xyz[2] = points[pos].fZ;
	  padrow = points[pos].fPadRow;
	  //AliL3Transform::Global2Local(xyz,slice,kTRUE);
	  AliL3Transform::Global2LocHLT(xyz,slice);
	  
	  Float_t angle = 0;
	  AliL3Transform::Local2GlobalAngle(&angle,slice);
	  if(!track->CalculateReferencePoint(angle,AliL3Transform::Row2X(padrow)))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluate::CalculateResiduals","Crossing point")
		<<"Track does not crossing padrow "<<padrow<<" in slice "<<slice<<ENDLOG;
	      continue;
	    }
	  
	  Float_t xyz_cross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  //AliL3Transform::Global2Local(xyz_cross,slice,kTRUE);	  
	  AliL3Transform::Global2LocHLT(xyz_cross,slice);
	  
	  Double_t beta = track->GetCrossingAngle(padrow,slice);
	  
	  Double_t yres = xyz_cross[1] - xyz[1];
	  Double_t zres = xyz_cross[2] - xyz[2];
	  Double_t dipangle = atan(track->GetTgl());
	  ntuppel->Fill(yres,zres,xyz_cross[2],track->GetPt(),dipangle,beta,padrow,track->GetNumberOfPoints());
	  
	}
    }
}

enum tagprimary {kPrimaryCharged = 0x4000};
void AliL3Evaluate::EvaluatePoints(Char_t *rootfile,Char_t *exactfile,Char_t *tofile,Int_t nevent,Bool_t offline,Bool_t sp)
{
  //Compare points to the exact crossing points of track and padrows.
  //The input file to this function, contains the exact clusters calculated
  //in AliTPC::Hits2ExactClusters.
  
#ifndef do_mc
  cerr<<"AliL3Evaluate::EvaluatePoints : Compile with do_mc flag!"<<endl;
  return;
#else
  cout<<"Evaluating points"<<endl;
  TNtuple *ntuppel = new TNtuple("ntuppel_res","Cluster properties",
				 "slice:padrow:charge:resy:resz:zHit:pt:beta:sigmaY2:sigmaZ2:psigmaY2:psigmaZ2");
  ntuppel->SetDirectory(0);
  
  TNtuple *ntuppel2 = new TNtuple("ntuppel_eff","Efficiency","slice:padrow:nfound:ngen");
  ntuppel2->SetDirectory(0);

  TFile *exfile = TFile::Open(rootfile);
  if(!exfile)
    {
      cerr<<"Error opening rootfile "<<rootfile<<endl;
      return;
    }
  gAlice = (AliRun*)exfile->Get("gAlice");
  if (!gAlice) 
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::InitMC","gAlice")
	<<"AliRun object non existing on file"<<ENDLOG;
      return;
    }

  AliTPCParam *param = (AliTPCParam*)exfile->Get(AliL3Transform::GetParamName());
  
  TFile *exact = TFile::Open(exactfile);
  if(!exact)
    {
      cerr<<"AliL3Evaluate::EvaluatePoints : Problems opening file :"<<exactfile<<endl;
      return;
    }
  
  AliStack *astack=gAlice->Stack();

  AliTPCClustersArray *arr=0;
  for(Int_t event=0; event<nevent; event++)
    {
      LoadData(event,sp);   
      exfile->cd();
      if(arr)
	delete arr;
      Int_t nparticles = gAlice->GetEvent(event);
      Int_t nprimaries = 0;//FindPrimaries(nparticles);
      cout<<"Event "<<event<<" had "<<nparticles<<" particles and "<<nprimaries<<" primaries"<<endl;
      exact->cd();
      
      //Get the exact clusters from file:
      AliTPCClustersArray *arr = new AliTPCClustersArray;
      arr->Setup(param);
      arr->SetClusterType("AliComplexCluster");
      char treeName[500];
      sprintf(treeName,"TreeCExact_%s_%d",param->GetTitle(),event);
      Bool_t clusterok = arr->ConnectTree(treeName);//Segment Tree (for offline clusters)
      if(!clusterok) {printf("AliL3Evaluate::EvaluatePoints : Error in clusterloading\n"); return;}
      
      //cout<<"Entering loop with "<<(Int_t)arr->GetTree()->GetEntries()<<endl;
      for(Int_t i=0; i<arr->GetTree()->GetEntries(); i++)
	{
	  //Get the exact clusters for this row:
	  Int_t cursec,currow;
	  AliSegmentID *s = arr->LoadEntry(i);
	  param->AdjustSectorRow(s->GetID(),cursec,currow);
	  
	  AliTPCClustersRow *ro = (AliTPCClustersRow *)arr->GetRow(cursec,currow);
	  TClonesArray *clusters = ro->GetArray();
	  int num_of_offline=clusters->GetEntriesFast();
	  
	  //Get the found clusters:
	  Int_t slice,padrow;
	  AliL3Transform::Sector2Slice(slice,padrow,cursec,currow);
	  if(slice < fMinSlice) continue;
	  if(slice > fMaxSlice) break;
	  
	  Int_t patch = AliL3Transform::GetPatch(padrow);
	  if(sp)
	    patch=0;
	  AliL3SpacePointData *points = fClusters[slice][patch];
	  if(!points)
	    continue;
	  
	  //cout<<"Slice "<<slice<<" padrow "<<padrow<<" has "<<num_of_offline<<" clusters "<<endl;
	  Int_t clustercount=0;
	  Int_t crosscount=0;
	  for(Int_t m=0; m<num_of_offline; m++)
	    {
	      AliComplexCluster *cluster = (AliComplexCluster *)clusters->UncheckedAt(m);
#ifdef use_newio
	      Int_t mcId = cluster->GetTrack(0);
#else
	      Int_t mcId = cluster->fTracks[0];
#endif	      
	      if(mcId <0) continue;
	
#ifdef use_newio      
	      if(cluster->GetY() < 1 || cluster->GetY() > AliL3Transform::GetNPads(padrow) - 2 ||
		 cluster->GetX() < 1 || cluster->GetX() > AliL3Transform::GetNTimeBins() - 2)
		continue;
#else
	      if(cluster->fY < 1 || cluster->fY > AliL3Transform::GetNPads(padrow) - 2 ||
		 cluster->fX < 1 || cluster->fX > AliL3Transform::GetNTimeBins() - 2)
		continue;
#endif	      
	      Float_t xyz_ex[3];
	      
#ifdef use_newio
	      AliL3Transform::Raw2Local(xyz_ex,cursec,currow,cluster->GetY(),cluster->GetX());
#else	      
	      AliL3Transform::Raw2Local(xyz_ex,cursec,currow,cluster->fY,cluster->fX);
#endif	      
	      //In function AliTPC::Hits2ExactClusters the time offset is not included,
	      //so we have to substract it again here.
	      if(slice<18)
		xyz_ex[2]-=AliL3Transform::GetZOffset();
	      else
		xyz_ex[2]+=AliL3Transform::GetZOffset();
	      
	      //Outside our cone:
	      if(param->GetPadRowRadii(cursec,currow)<230./250.*fabs(xyz_ex[2]))
		continue;
	      
	      TParticle *part = astack->Particle(mcId);
	      crosscount++;
	      
	      if(part->Pt() < fMinGoodPt) continue;
	      
	      //Dont take secondaries, because in width calculation we assume primaries:
	      //if(!(part->TestBit(kPrimaryCharged))) continue;
	      if(part->GetFirstMother()>=0) continue;
	      
	      Int_t tempcount=0;
	      for(UInt_t c=0; c<fNcl[slice][patch]; c++)
		{
		  if((Int_t)points[c].fPadRow!=padrow) continue;
		  Float_t xyz_cl[3] = {points[c].fX,points[c].fY,points[c].fZ};
		  
		  if(!offline)
		    AliL3Transform::Global2Local(xyz_cl,cursec);
		  tempcount++;
		  
		  if(points[c].fTrackID[0] != mcId &&
		     points[c].fTrackID[1] != mcId &&
		     points[c].fTrackID[2] != mcId)
		    continue;
		  
		  //Residuals:
		  Float_t resy = xyz_cl[1] - xyz_ex[1];
		  Float_t resz = xyz_cl[2] - xyz_ex[2];
		  
		  //Cluster shape
		  Int_t charge = (Int_t)points[c].fCharge;
		  Float_t beta = GetCrossingAngle(part,slice,padrow,xyz_ex);
		  Double_t tanl = xyz_ex[2]/sqrt(xyz_ex[0]*xyz_ex[0]+xyz_ex[1]*xyz_ex[1]);
		  Float_t psigmaY2 = AliL3Transform::GetParSigmaY2(padrow,xyz_ex[2],beta);
		  Float_t psigmaZ2 = AliL3Transform::GetParSigmaZ2(padrow,xyz_ex[2],tanl);
		  Float_t sigmaY2 = points[c].fSigmaY2;
		  Float_t sigmaZ2 = points[c].fSigmaZ2;
		  ntuppel->Fill(slice,padrow,charge,resy,resz,xyz_ex[2],part->Pt(),beta,sigmaY2,sigmaZ2,psigmaY2,psigmaZ2);
		}
	      clustercount=tempcount;
	    }
	  ntuppel2->Fill(slice,padrow,clustercount,crosscount);
	  arr->ClearRow(cursec,currow);
	}
    }
  exfile->Close();
  exact->Close();

  TFile *ofile = TFile::Open(tofile,"RECREATE");
  ntuppel->Write();
  ntuppel2->Write();
  ofile->Close();
  
#endif
}

void AliL3Evaluate::GetCFeff(Char_t *path,Char_t *outfile,Int_t nevent,Bool_t sp)
{
  //Evaluate the cluster finder efficiency.
  
#ifndef do_mc
  cerr<<"AliL3Evaluate::GetCFeff : Compile with do_mc flag"<<endl;
  return;
#else
  TNtuple *ntuppel = new TNtuple("ntuppel","Cluster finder efficiency","slice:row:ncrossings:nclusters");
  ntuppel->SetDirectory(0);
  
  Char_t filename[1024];
  sprintf(filename,"%s/alirunfile.root",path);
  TFile *rfile = TFile::Open(filename);
  gAlice = (AliRun*)rfile->Get("gAlice");

  AliStack *astack=gAlice->Stack();
  
  AliTPCParam *param = (AliTPCParam*)rfile->Get(AliL3Transform::GetParamName());
      
  Int_t zero=param->GetZeroSup();

  sprintf(filename,"%s/digitfile.root",path);
  TFile *dfile = TFile::Open(filename);
  
  for(Int_t event=0; event<nevent; event++)
    {
      LoadData(event,sp);
      rfile->cd();
      gAlice->GetEvent(event);
      Int_t np = astack->GetNtrack();
      cout<<"Processing event "<<event<<" with "<<np<<" particles "<<endl;
      dfile->cd();
      sprintf(filename,"TreeD_75x40_100x60_150x60_%d",event);
      TTree *TD=(TTree*)gDirectory->Get(filename);
      AliSimDigits da, *digits=&da;
      TD->GetBranch("Segment")->SetAddress(&digits);
      
      Int_t crossed=0,recs=0;
      Int_t *count = new Int_t[np]; //np number of particles.
      Int_t i;
      Float_t xyz[3];
      for (i=0; i<np; i++) count[i]=0;
      
      
      Int_t sec,row,sl,sr;
      for(Int_t i=0; i<(Int_t)TD->GetEntries(); i++)
	{
	  crossed=recs=0;
	  if (!TD->GetEvent(i)) continue;
	  param->AdjustSectorRow(digits->GetID(),sec,row);
	  AliL3Transform::Sector2Slice(sl,sr,sec,row);
	  if(sl < fMinSlice) continue;
	  if(sl > fMaxSlice) break;
	  cout<<"Processing slice "<<sl<<" row "<<sr<<endl;
	  digits->First();
	  do {
	    Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
	    Short_t dig = digits->GetDigit(it,ip);
	    
	    if(dig<=param->GetZeroSup()) continue;
	    AliL3Transform::Raw2Local(xyz,sec,row,ip,it);
	    if(param->GetPadRowRadii(sec,row)<230./250.*fabs(xyz[2]))
	      continue;
	    
	    Int_t idx0=digits->GetTrackID(it,ip,0); 
	    Int_t idx1=digits->GetTrackID(it,ip,1);
	    Int_t idx2=digits->GetTrackID(it,ip,2);
	    
	    if (idx0>=0 && dig>=zero) count[idx0]+=1;
	    if (idx1>=0 && dig>=zero) count[idx1]+=1;
	    if (idx2>=0 && dig>=zero) count[idx2]+=1;
	  } while (digits->Next());
	  for (Int_t j=0; j<np; j++) 
	    {
	      TParticle *part = astack->Particle(j);
	      if(part->Pt() < fMinGoodPt) continue;
	      if(part->GetFirstMother() >= 0) continue;
	      if (count[j]>1) //at least two digits at this padrow 
		{
		  crossed++;
		  count[j]=0;
		}
	    }

	  Int_t patch = AliL3Transform::GetPatch(sr);
	  if(sp==kTRUE)
	    patch=0;
	  AliL3SpacePointData *points = fClusters[sl][patch];
	  if(!points)
	    continue;
	  for(UInt_t k=0; k<fNcl[sl][patch]; k++)
	    {
	      if(points[k].fPadRow!=sr) continue;
	      recs++;
	    }
	  ntuppel->Fill(sl,sr,crossed,recs);
	}
      
      TD->Delete();
      delete[] count;
    }
  TFile *file = TFile::Open(outfile,"RECREATE");
  ntuppel->Write();
  file->Close();
  
  rfile->Close();
  dfile->Close();
#endif
}

Float_t AliL3Evaluate::GetCrossingAngle(TParticle *part,Int_t slice,Int_t /*padrow*/,Float_t *xyz)
{
  //Calculate the padrow crossing angle of the particle
  
  Double_t kappa = AliL3Transform::GetBField()*AliL3Transform::GetBFact()/part->Pt();
  
  Double_t radius = 1/fabs(kappa);
  if(part->GetPdgCode() > 0) kappa = -kappa;

  Float_t angl[1] = {part->Phi()};
  
  AliL3Transform::Global2LocalAngle(angl,slice);
  
  Double_t charge = -1.*kappa;

  Double_t trackPhi0 = angl[0] + charge*0.5*AliL3Transform::Pi()/fabs(charge);
  
  Double_t x0=0;
  Double_t y0=0;
  Double_t xc = x0 - radius * cos(trackPhi0);
  Double_t yc = y0 - radius * sin(trackPhi0);

  Double_t tangent[2];
  tangent[0] = -1.*(xyz[1] - yc)/radius;
  tangent[1] = (xyz[0] - xc)/radius;
  
  Double_t perp_padrow[2] = {1,0}; //locally in slice
  
  Double_t cos_beta = fabs(tangent[0]*perp_padrow[0] + tangent[1]*perp_padrow[1]);
  if(cos_beta > 1) cos_beta=1;
  return acos(cos_beta);
}

Int_t AliL3Evaluate::FindPrimaries(Int_t nparticles)
{
  // cuts:
  Double_t vertcut = 0.001;
  Double_t decacut = 3.;
  Double_t timecut = 0.;
  Int_t nprch1=0;
  AliStack *astack=gAlice->Stack();
  TParticle * part = astack->Particle(0);
  Double_t xori = part->Vx();
  Double_t yori = part->Vy();
  Double_t zori = part->Vz();
  for(Int_t iprim = 0; iprim<nparticles; iprim++){   //loop on  tracks
    
    part = astack->Particle(iprim);
    const char * xxx=strstr(part->GetName(),"XXX");
    if(xxx)continue;
    
    TParticlePDG *ppdg = part->GetPDG();
    if(TMath::Abs(ppdg->Charge())!=3)continue;  // only charged (no quarks)
    
    Double_t dist=TMath::Sqrt((part->Vx()-xori)*(part->Vx()-xori)+(part->Vy()-yori)*(part->Vy()-yori)+(part->Vz()-zori)*(part->Vz()-zori));
    if(dist>vertcut)continue;  // cut on the vertex
    
    if(part->T()>timecut)continue;
    
    Double_t ptot=TMath::Sqrt(part->Px()*part->Px()+part->Py()*part->Py()+part->Pz()*part->Pz());
    if(ptot==(TMath::Abs(part->Pz())))continue; // no beam particles
    
    Bool_t prmch = kTRUE;   // candidate primary track
    Int_t fidau=part->GetFirstDaughter();  // cut on daughters
    Int_t lasdau=0;
    Int_t ndau=0;
    if(fidau>=0){
      lasdau=part->GetLastDaughter();
      ndau=lasdau-fidau+1;
    }
    if(ndau>0){
      for(Int_t j=fidau;j<=lasdau;j++){
        TParticle *dau=astack->Particle(j);
        Double_t distd=TMath::Sqrt((dau->Vx()-xori)*(dau->Vx()-xori)+(dau->Vy()-yori)*(dau->Vy()-yori)+(dau->Vz()-zori)*(dau->Vz()-zori));
        if(distd<decacut)prmch=kFALSE;  // eliminate if the decay is near the vertex
      }
    }
    
    if(prmch){
      nprch1++;
      part->SetBit(kPrimaryCharged);
    }
  }

  return nprch1;
}
