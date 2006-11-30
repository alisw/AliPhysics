// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

/** \class AliHLTEvaluate
<pre>
//_____________________________________________________________
// AliHLTEvaluate
//
// Evaluation class for tracking; plots, efficiencies etc..
//
</pre>
*/

#include <TObject.h>
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TNtuple.h>

#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPC.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersArray.h>
#include <AliTPCClustersRow.h>
#include <AliTPCParam.h>
#include <AliComplexCluster.h>
#include <AliStack.h>

#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTTrack.h"
#include "AliHLTFileHandler.h"
#include "AliHLTTrackArray.h"
#include "AliHLTEvaluate.h"

#if __GNUC__ == 3
#include <iosfwd>
using namespace std;
#endif

ClassImp(AliHLTEvaluate)

AliHLTEvaluate::AliHLTEvaluate()
{ 
  //constructor
  Clear();
}

AliHLTEvaluate::AliHLTEvaluate(Char_t *datapath,Int_t minclusters,Int_t minhits,Double_t minpt,Double_t maxpt,Int_t *slice)
{ 
  //constructor
  Clear();

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
  fMinPointsOnTrack = minclusters;
  fMinHitsFromParticle = minhits;
  fMinGoodPt = minpt;
  fMaxGoodPt = maxpt;
}

AliHLTEvaluate::~AliHLTEvaluate()
{ 
  //destructor
  if(fGoodTracks) delete[] fGoodTracks;
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
      {
	//cout << s << " " << p <<endl;
	if(fClustersFile[s][p]) delete fClustersFile[s][p];
	fClusters[s][p]=0;
      }
}

void AliHLTEvaluate::Clear()
{ 
  //clear
  fTracks = 0;
  fMinSlice=0;
  fMaxSlice=0;
  fNFastPoints = 0;
  fMcIndex = 0;
  fMcId = 0;
  fGoodFound = 0;
  fGoodGen = 0;
  fMinGoodPt = 0;
  fMaxGoodPt = 0;
  fMinPointsOnTrack = 0;  
  fMinHitsFromParticle = 0;
  fGoodTracks = 0;
  fMaxFalseClusters = 0.1;

  fNtuppel=0;
  fPtRes=0;
  fNGoodTracksPt = 0;
  fNFoundTracksPt = 0;
  fNFakeTracksPt = 0;
  fTrackEffPt = 0;
  fFakeTrackEffPt = 0;
  fNGoodTracksEta = 0;
  fNFoundTracksEta = 0;
  fNFakeTracksEta = 0;
  fTrackEffEta = 0;
  fFakeTrackEffEta = 0;
  fNtupleRes = 0;

  fStandardComparison=kTRUE;
  for(Int_t s=0; s<=35; s++)
    for(Int_t p=0; p<6; p++){
      fClusters[s][p]=0;
      fClustersFile[s][p]=0;
    }
  sprintf(fPath,"./");
}

void AliHLTEvaluate::LoadData(Int_t event,Bool_t sp)
{ 
  //load cluster points
  Char_t fname[1024];

  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliHLTTransform::GetNPatches(); p++)
	{
	  Int_t patch;
	  if(sp==kTRUE)
	    patch=-1;
	  else
	    patch=p;
	  
	  if(fClustersFile[s][p])
	    delete fClustersFile[s][p];
	  fClusters[s][p] = 0;
	  fClustersFile[s][p] = new AliHLTFileHandler();
	  if(event<0)
	    sprintf(fname,"%s/points_%d_%d.raw",fPath,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",fPath,event,s,patch);
	  if(!fClustersFile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliHLTLog::kError,"AliHLTEvaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
              delete fClustersFile[s][p];
              fClustersFile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliHLTSpacePointData*)fClustersFile[s][p]->Allocate();
	  fClustersFile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  fClustersFile[s][p]->CloseBinaryInput();
	  if(sp==kTRUE)
	    break;
	}
    }
   
  sprintf(fname,"%s/tracks_%d.raw",fPath,event);
  AliHLTFileHandler *tfile = new AliHLTFileHandler();
  if(!tfile->SetBinaryInput(fname)){
    LOG(AliHLTLog::kError,"AliHLTEvaluation::Setup","File Open")
      <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
    return;
  }
  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  fTracks->QSort();
  delete tfile;
}

void AliHLTEvaluate::AssignPIDs()
{ 
  //assign pid 
  if(!fTracks) return;
  fTracks->QSort();
  LOG(AliHLTLog::kDebug,"AliHLTEvaluate::AssignPIDs","Track Loop")
    <<"Assigning pid to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTTrack *track = (AliHLTTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      if(track->GetNumberOfPoints() < fMinPointsOnTrack)
	track->SetPID(0);
      else {
	Float_t pid = GetTrackPID(track);
	track->SetPID(pid);
      }
    }
}
  
void AliHLTEvaluate::AssignIDs()
{ 
  //Assign MC id to the tracks.
#ifndef do_mc
  cerr<<"AliHLTEvaluate::AssignIDs() : You need to compile with the do_mc flag!"<<endl;
  return;
#else
  if(!fTracks) return;
  fGoodFound=0;
  fTracks->QSort();
  LOG(AliHLTLog::kDebug,"AliHLTEvaluate::AssignIDs","Track Loop")
    <<"Assigning MC id to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTTrack *track = (AliHLTTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      if(track->GetNumberOfPoints() < fMinPointsOnTrack) break;
      
      fGoodFound++;
      Int_t tID = GetMCTrackLabel(track);
      track->SetMCid(tID);
    }
  //cout<<"Found "<<fGoodFound<<" good tracks "<<endl;
#endif
}

Float_t AliHLTEvaluate::GetTrackPID(AliHLTTrack *track)
{ 
  //get track pid
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
    AliHLTSpacePointData *points = fClusters[iSector][patch];
    if(!points) continue; 
    if(position>=fNcl[iSector][patch]) 
      {
	LOG(AliHLTLog::kError,"AliHLTEvaluate::GetMCTrackLabel","Clusterarray")
	  <<AliHLTLog::kDec<<"ERROR"<<ENDLOG;
	continue;
      }
    UChar_t padrow = points[position].fPadRow;
    Float_t pWidth = AliHLTTransform::GetPadPitchWidthLow();
    if (padrow>63)
      pWidth = AliHLTTransform::GetPadPitchWidthUp(); 
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

#ifdef do_mc
Int_t AliHLTEvaluate::GetMCTrackLabel(AliHLTTrack *track)
{ 
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.
  //Definitions are identical to offline.
  //Fake track means:
  // - more than 10 percent of clusters are assigned incorrectly
  // - more than half of the innermost 10% of clusters were assigned incorrectly.
  
  
  Int_t numofclusters = track->GetNumberOfPoints();
  AliS *s=new AliS[numofclusters];
  Int_t i;
  for (i=0; i<numofclusters; i++) s[i].flab=s[i].fmax=0;
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;
    
  Int_t lab=123456789;
  for (i=0; i<numofclusters; i++) 
    {
      //Tricks to get the clusters belonging to this track:
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliHLTSpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliHLTLog::kError,"AliHLTEvaluate::GetMCTrackLabel","Clusterarray")
	    <<AliHLTLog::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      //Get the label of the cluster:
      lab=points[pos].fTrackID[0];

      Int_t j;
      for (j=0; j<numofclusters; j++)
        if (s[j].flab==lab || s[j].fmax==0) break;
      s[j].flab=lab;
      s[j].fmax++;
    }
  
  Int_t max=0;
  for (i=0; i<numofclusters; i++) 
    if (s[i].fmax>max) {max=s[i].fmax; lab=s[i].flab;}
  
  if(lab == -1)
    return -1; //If most clusters is -1, this is a noise track.
  if(lab < 0)
    cerr<<"AliHLTEvaluate::GetMCTrackLabel : Track label negative :"<<lab<<endl;
  
  delete[] s;
  
  for (i=0; i<numofclusters; i++) 
    {
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliHLTSpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliHLTLog::kError,"AliHLTEvaluate::GetMCTrackLabel","Clusterarray")
	    <<AliHLTLog::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      if (abs(points[pos].fTrackID[1]) == lab || 
	  abs(points[pos].fTrackID[2]) == lab ) max++;
    }
  
  
  //Check if more than 10% of the clusters were assigned incorrectly:
  if (1.-Float_t(max)/numofclusters > fMaxFalseClusters) 
    {
      return -lab;
    }
  else //Check if at least half of the 10% innermost clusters are assigned correctly.
    {
      Int_t tail=Int_t(0.10*numofclusters);
      max=0;
      for (i=1; i<=tail; i++) 
	{
	  id = hitnum[numofclusters - i];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  
	  AliHLTSpacePointData *points = fClusters[slice][patch];
	  if(!points) continue;
	  if(lab == abs(points[pos].fTrackID[0]) ||
	     lab == abs(points[pos].fTrackID[1]) ||
	     lab == abs(points[pos].fTrackID[2])) max++;
	}
      if (max < Int_t(0.5*tail)) return -lab;
    }

  return lab;
#else //If we are running with mc_ids or not
  Int_t AliHLTEvaluate::GetMCTrackLabel(AliHLTTrack */*track*/)
{ 
  // Does nothing if do_mc undefined
  return 0;
#endif

}

void AliHLTEvaluate::GetFastClusterIDs(Char_t *path)
{
  //Get the MC id of space points in case of using the fast simulator. 
  char name[256];
  sprintf(name,"%s/point_mc.dat",path);
  FILE *infile = fopen(name,"r");
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

void AliHLTEvaluate::CreateHistos(Int_t nbin,Float_t xlow,Float_t xup)
{
  //Create the histograms 
  
  LOG(AliHLTLog::kInformational,"AliHLTEvaluate::CreateHistos","Allocating")
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

void AliHLTEvaluate::GetGoodParticles(Char_t *path,Int_t event,Int_t *padrowrange)
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
      cerr<<"AliHLTEvaluate::GetGoodParticles : Problems opening file :"<<filename<<endl;
      return;
    }
  Int_t MaxTracks=20000;
  if(fGoodTracks)
    delete [] fGoodTracks;
  fGoodGen=0;
  fGoodTracks = new AliGoodTrack[MaxTracks];
  
  if(fStandardComparison){
    while (in>>fGoodTracks[fGoodGen].flabel>>fGoodTracks[fGoodGen].fcode>>
	   fGoodTracks[fGoodGen].fpx>>fGoodTracks[fGoodGen].fpy>>fGoodTracks[fGoodGen].fpz>>
	   fGoodTracks[fGoodGen].fx>>fGoodTracks[fGoodGen].fy>>fGoodTracks[fGoodGen].fz)
      {
        fGoodTracks[fGoodGen].fnhits=-1;
        fGoodTracks[fGoodGen].fsector=-1; 
	fGoodGen++;
	if (fGoodGen==MaxTracks) 
	  {
	    cerr<<"AliHLTEvaluate::GetGoodParticles : Too many good tracks !\n";
	    break;
	  }
      }
  } else {
    while (in>>fGoodTracks[fGoodGen].flabel>>fGoodTracks[fGoodGen].fcode>>
	   fGoodTracks[fGoodGen].fpx>>fGoodTracks[fGoodGen].fpy>>fGoodTracks[fGoodGen].fpz>>
	   fGoodTracks[fGoodGen].fx>>fGoodTracks[fGoodGen].fy >>fGoodTracks[fGoodGen].fz>>
	   fGoodTracks[fGoodGen].fnhits>>fGoodTracks[fGoodGen].fsector) 
      {
	fGoodGen++;
	if (fGoodGen==MaxTracks) 
	  {
	    cerr<<"AliHLTEvaluate::GetGoodParticles : Too many good tracks !\n";
	    break;
	  }
      }
  }
}

void AliHLTEvaluate::FillEffHistos()
{ 
  //has to be modified for fakes.

  if(!fGoodTracks)
    {
      cerr<<"AliHLTEvaluate::FillEffHistos : No good tracks"<<endl;
      return;
    }
  if(!fTracks) return;

  //cout<<"Comparing "<<fGoodGen<<" good tracks ..."<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      //cout<<"Checking particle "<<i<<endl;
      if(!fStandardComparison) 
	if(fGoodTracks[i].fnhits < fMinHitsFromParticle) continue;
      Float_t ptg = TMath::Sqrt(fGoodTracks[i].fpx*fGoodTracks[i].fpx + fGoodTracks[i].fpy*fGoodTracks[i].fpy);
      if(ptg < fMinGoodPt || ptg > fMaxGoodPt) continue;
      Float_t pzg=fGoodTracks[i].fpz;
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
	  AliHLTTrack *track = fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  Int_t nHits = track->GetNumberOfPoints();
	  if(nHits < fMinPointsOnTrack) break;
	  Int_t tracklabel;
	  tracklabel = track->GetMCid();
	  
	  if(TMath::Abs(tracklabel) != fGoodTracks[i].flabel) continue;
	  found=1;
	  Float_t pt=track->GetPt();
	  if(tracklabel == fGoodTracks[i].flabel) 
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

void AliHLTEvaluate::FillEffHistosNAIVE()
{  
  //Fill the efficiency histograms.
  
  cout<<endl<<"Note: Doing NAIVE evaluation "<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      if(!fStandardComparison) 
	if(fGoodTracks[i].fnhits < fMinHitsFromParticle) continue;
      Double_t ptg=TMath::Sqrt(fGoodTracks[i].fpx*fGoodTracks[i].fpx + fGoodTracks[i].fpy*fGoodTracks[i].fpy);
      if(ptg < fMinGoodPt || ptg > fMaxGoodPt) continue;
      Double_t pzg=fGoodTracks[i].fpz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();
      //printf("filling particle with pt %f and dipangle %f\n",ptg,dipangle);
      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(dipangle);
      
    }
  
  for(Int_t k=0; k<fTracks->GetNTracks(); k++)
    {
      AliHLTTrack *track = fTracks->GetCheckedTrack(k);
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

void AliHLTEvaluate::CalcEffHistos()
{ 
  //calc eff histos

  Stat_t ngood=fNGoodTracksPt->GetEntries();
  Stat_t nfound=fNFoundTracksPt->GetEntries();
  Stat_t nfake=fNFakeTracksPt->GetEntries();

  LOG(AliHLTLog::kInformational,"AliHLTEvaluate::FillEffHistos","Efficiency")
    <<AliHLTLog::kDec<<"There was "<<ngood<<" generated good tracks"<<ENDLOG;
  LOG(AliHLTLog::kInformational,"AliHLTEvaluate::FillEffHistos","Efficiency")
    <<AliHLTLog::kDec<<"Found "<<nfound<<" tracks"<<ENDLOG;
  LOG(AliHLTLog::kInformational,"AliHLTEvaluate::FillEffHistos","Efficiency")
    <<AliHLTLog::kDec<<"Integral efficiency is about "<<nfound/ngood*100<<ENDLOG;
  LOG(AliHLTLog::kInformational,"AliHLTEvaluate::FillEffHistos","Efficiency")
    <<AliHLTLog::kDec<<"Fake tracks relative is about "<<nfake/ngood*100<<ENDLOG;
  //LOG(AliHLTLog::kInformational,"AliHLTEvaluate::FillEffHistos","Efficiency")
  //<<AliHLTLog::kDec<<"Naive efficiency "<<(Double_t)fGoodFound/(Double_t)fGoodGen<<ENDLOG;

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

void AliHLTEvaluate::Write2File(Char_t *outputfile)
{
  //Write histograms to file:
  
  TFile *of = TFile::Open(outputfile,"RECREATE");
  if(!of->IsOpen())
    {
      LOG(AliHLTLog::kError,"AliHLTEvaluate::Write2File","File Open")
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

TNtuple *AliHLTEvaluate::GetNtuple()
{ 
  //get ntuple
  if(!fNtupleRes)
    {
      fNtupleRes = new TNtuple("ntuppel","Residuals","residual_trans:residual_long:zHit:pt:dipangle:beta:padrow:nHits");
      fNtupleRes->SetDirectory(0);//Bug in older version of root.
    }
  return fNtupleRes;
}

void AliHLTEvaluate::CalculateResiduals()
{ 
  //calculate residuals
  TNtuple *ntuppel = GetNtuple();
  
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      
      AliHLTTrack *track = (AliHLTTrack*)fTracks->GetCheckedTrack(i);
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
	  
	  AliHLTSpacePointData *points = fClusters[slice][patch];
	  if(!points) 
	    {
	      LOG(AliHLTLog::kError,"AliHLTEvaluate::CalculateResiduals","Clusterarray")
		<<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	      continue;
	    }
	  if(pos>=fNcl[slice][patch]) 
	    {
	      LOG(AliHLTLog::kError,"AliHLTEvaluate::CalculateResiduals","Clusterarray")
		<<AliHLTLog::kDec<<"ERROR"<<ENDLOG;
	      continue;
	    }
	  
	  xyz[0] = points[pos].fX;
	  xyz[1] = points[pos].fY;
	  xyz[2] = points[pos].fZ;
	  padrow = points[pos].fPadRow;
	  //AliHLTTransform::Global2Local(xyz,slice,kTRUE);
	  AliHLTTransform::Global2LocHLT(xyz,slice);
	  
	  Float_t angle = 0;
	  AliHLTTransform::Local2GlobalAngle(&angle,slice);
	  if(!track->CalculateReferencePoint(angle,AliHLTTransform::Row2X(padrow)))
	    {
	      LOG(AliHLTLog::kError,"AliHLTEvaluate::CalculateResiduals","Crossing point")
		<<"Track does not crossing padrow "<<padrow<<" in slice "<<slice<<ENDLOG;
	      continue;
	    }
	  
	  Float_t xyzcross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	  //AliHLTTransform::Global2Local(xyzcross,slice,kTRUE);	  
	  AliHLTTransform::Global2LocHLT(xyzcross,slice);
	  
	  Double_t beta = track->GetCrossingAngle(padrow,slice);
	  
	  Double_t yres = xyzcross[1] - xyz[1];
	  Double_t zres = xyzcross[2] - xyz[2];
	  Double_t dipangle = atan(track->GetTgl());
	  ntuppel->Fill(yres,zres,xyzcross[2],track->GetPt(),dipangle,beta,padrow,track->GetNumberOfPoints());
	  
	}
    }
}

enum tagprimary {kPrimaryCharged = 0x4000};
#ifndef do_mc
void AliHLTEvaluate::EvaluatePoints(Char_t */*rootfile*/,Char_t */*exactfile*/,Char_t */*tofile*/,Int_t /*nevent*/,Bool_t /*offline*/,Bool_t /*sp*/)
{
  // Does nothing if do_mc undefined
  
  cerr<<"AliHLTEvaluate::EvaluatePoints : Compile with do_mc flag!"<<endl;
  return;
#else
void AliHLTEvaluate::EvaluatePoints(Char_t *rootfile,Char_t *exactfile,Char_t *tofile,Int_t nevent,Bool_t offline,Bool_t sp)
{
  //Compare points to the exact crossing points of track and padrows.
  //The input file to this function, contains the exact clusters calculated
  //in AliTPC::Hits2ExactClusters.
  
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
      LOG(AliHLTLog::kError,"AliHLTEvaluate::InitMC","gAlice")
	<<"AliRun object non existing on file"<<ENDLOG;
      return;
    }

  AliTPCParam *param = (AliTPCParam*)exfile->Get(AliHLTTransform::GetParamName());
  
  TFile *exact = TFile::Open(exactfile);
  if(!exact)
    {
      cerr<<"AliHLTEvaluate::EvaluatePoints : Problems opening file :"<<exactfile<<endl;
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
      if(!clusterok) {printf("AliHLTEvaluate::EvaluatePoints : Error in clusterloading\n"); return;}
      
      //cout<<"Entering loop with "<<(Int_t)arr->GetTree()->GetEntries()<<endl;
      for(Int_t i=0; i<arr->GetTree()->GetEntries(); i++)
	{
	  //Get the exact clusters for this row:
	  Int_t cursec,currow;
	  AliSegmentID *s = arr->LoadEntry(i);
	  param->AdjustSectorRow(s->GetID(),cursec,currow);
	  
	  AliTPCClustersRow *ro = (AliTPCClustersRow *)arr->GetRow(cursec,currow);
	  TClonesArray *clusters = ro->GetArray();
	  int numofoffline=clusters->GetEntriesFast();
	  
	  //Get the found clusters:
	  Int_t slice,padrow;
	  AliHLTTransform::Sector2Slice(slice,padrow,cursec,currow);
	  if(slice < fMinSlice) continue;
	  if(slice > fMaxSlice) break;
	  
	  Int_t patch = AliHLTTransform::GetPatch(padrow);
	  if(sp)
	    patch=0;
	  AliHLTSpacePointData *points = fClusters[slice][patch];
	  if(!points)
	    continue;
	  
	  //cout<<"Slice "<<slice<<" padrow "<<padrow<<" has "<<numofoffline<<" clusters "<<endl;
	  Int_t clustercount=0;
	  Int_t crosscount=0;
	  for(Int_t m=0; m<numofoffline; m++)
	    {
	      AliComplexCluster *cluster = (AliComplexCluster *)clusters->UncheckedAt(m);
#ifdef use_newio
	      Int_t mcId = cluster->GetTrack(0);
#else
	      Int_t mcId = cluster->fTracks[0];
#endif	      
	      if(mcId <0) continue;
	
#ifdef use_newio      
	      if(cluster->GetY() < 1 || cluster->GetY() > AliHLTTransform::GetNPads(padrow) - 2 ||
		 cluster->GetX() < 1 || cluster->GetX() > AliHLTTransform::GetNTimeBins() - 2)
		continue;
#else
	      if(cluster->fY < 1 || cluster->fY > AliHLTTransform::GetNPads(padrow) - 2 ||
		 cluster->fX < 1 || cluster->fX > AliHLTTransform::GetNTimeBins() - 2)
		continue;
#endif	      
	      Float_t xyzex[3];
	      
#ifdef use_newio
	      AliHLTTransform::Raw2Local(xyzex,cursec,currow,cluster->GetY(),cluster->GetX());
#else	      
	      AliHLTTransform::Raw2Local(xyzex,cursec,currow,cluster->fY,cluster->fX);
#endif	      
	      //In function AliTPC::Hits2ExactClusters the time offset is not included,
	      //so we have to substract it again here.
	      if(slice<18)
		xyzex[2]-=AliHLTTransform::GetZOffset();
	      else
		xyzex[2]+=AliHLTTransform::GetZOffset();
	      
	      //Outside our cone:
	      if(param->GetPadRowRadii(cursec,currow)<230./250.*fabs(xyzex[2]))
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
		  Float_t xyzcl[3] = {points[c].fX,points[c].fY,points[c].fZ};
		  
		  if(!offline)
		    AliHLTTransform::Global2Local(xyzcl,cursec);
		  tempcount++;
		  
		  if(points[c].fTrackID[0] != mcId &&
		     points[c].fTrackID[1] != mcId &&
		     points[c].fTrackID[2] != mcId)
		    continue;
		  
		  //Residuals:
		  Float_t resy = xyzcl[1] - xyzex[1];
		  Float_t resz = xyzcl[2] - xyzex[2];
		  
		  //Cluster shape
		  Int_t charge = (Int_t)points[c].fCharge;
		  Float_t beta = GetCrossingAngle(part,slice,padrow,xyzex);
		  Double_t tanl = xyzex[2]/sqrt(xyzex[0]*xyzex[0]+xyzex[1]*xyzex[1]);
		  Float_t psigmaY2 = AliHLTTransform::GetParSigmaY2(padrow,xyzex[2],beta);
		  Float_t psigmaZ2 = AliHLTTransform::GetParSigmaZ2(padrow,xyzex[2],tanl);
		  Float_t sigmaY2 = points[c].fSigmaY2;
		  Float_t sigmaZ2 = points[c].fSigmaZ2;
		  ntuppel->Fill(slice,padrow,charge,resy,resz,xyzex[2],part->Pt(),beta,sigmaY2,sigmaZ2,psigmaY2,psigmaZ2);
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

#ifndef do_mc
void AliHLTEvaluate::GetCFeff(Char_t */*path*/,Char_t */*outfile*/,Int_t /*nevent*/,Bool_t /*sp*/)
{
  // Does nothing if do_mc undefined
  
  cerr<<"AliHLTEvaluate::GetCFeff : Compile with do_mc flag"<<endl;
  return;
#else
void AliHLTEvaluate::GetCFeff(Char_t *path,Char_t *outfile,Int_t nevent,Bool_t sp)
{
  //Evaluate the cluster finder efficiency.
  
  TNtuple *ntuppel = new TNtuple("ntuppel","Cluster finder efficiency","slice:row:ncrossings:nclusters");
  ntuppel->SetDirectory(0);
  
  Char_t filename[1024];
  sprintf(filename,"%s/alirunfile.root",path);
  TFile *rfile = TFile::Open(filename);
  gAlice = (AliRun*)rfile->Get("gAlice");

  AliStack *astack=gAlice->Stack();
  
  AliTPCParam *param = (AliTPCParam*)rfile->Get(AliHLTTransform::GetParamName());
      
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
	  AliHLTTransform::Sector2Slice(sl,sr,sec,row);
	  if(sl < fMinSlice) continue;
	  if(sl > fMaxSlice) break;
	  cout<<"Processing slice "<<sl<<" row "<<sr<<endl;
	  digits->First();
	  do {
	    Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
	    Short_t dig = digits->GetDigit(it,ip);
	    
	    if(dig<=param->GetZeroSup()) continue;
	    AliHLTTransform::Raw2Local(xyz,sec,row,ip,it);
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

	  Int_t patch = AliHLTTransform::GetPatch(sr);
	  if(sp==kTRUE)
	    patch=0;
	  AliHLTSpacePointData *points = fClusters[sl][patch];
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

Float_t AliHLTEvaluate::GetCrossingAngle(TParticle *part,Int_t slice,Int_t /*padrow*/,Float_t *xyz)
{
  //Calculate the padrow crossing angle of the particle
  
  Double_t kappa = AliHLTTransform::GetBField()*AliHLTTransform::GetBFact()/part->Pt();
  
  Double_t radius = 1/fabs(kappa);
  if(part->GetPdgCode() > 0) kappa = -kappa;

  Float_t angl[1] = {part->Phi()};
  
  AliHLTTransform::Global2LocalAngle(angl,slice);
  
  Double_t charge = -1.*kappa;

  Double_t trackPhi0 = angl[0] + charge*0.5*AliHLTTransform::Pi()/fabs(charge);
  
  Double_t x0=0;
  Double_t y0=0;
  Double_t xc = x0 - radius * cos(trackPhi0);
  Double_t yc = y0 - radius * sin(trackPhi0);

  Double_t tangent[2];
  tangent[0] = -1.*(xyz[1] - yc)/radius;
  tangent[1] = (xyz[0] - xc)/radius;
  
  Double_t perppadrow[2] = {1,0}; //locally in slice
  
  Double_t cosbeta = fabs(tangent[0]*perppadrow[0] + tangent[1]*perppadrow[1]);
  if(cosbeta > 1) cosbeta=1;
  return acos(cosbeta);
}

Int_t AliHLTEvaluate::FindPrimaries(Int_t nparticles)
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
