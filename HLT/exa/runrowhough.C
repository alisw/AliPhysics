// $Id$

#ifndef __CINT__
#include "AliL3Logger.h"
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"
#include "AliL3Hough.h"
#include "AliL3TrackArray.h"
#include "AliL3Track.h"
#include "AliL3HoughTrack.h"
#include "AliL3Fitter.h"
#include "AliL3ClusterFitter.h"
#include "AliL3Vertex.h"
#include "AliL3Benchmark.h"
#include <AliRunLoader.h>
#include <AliStack.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <TBenchmark.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#endif

Int_t runrowhough(Char_t *path="./",Char_t *outpath="./fitter",int s1=0,int s2=35,int nevent=1,Bool_t skip=kTRUE)
{
  Bool_t isinit=AliL3Transform::Init(path,kTRUE);
  if(!isinit){
    cerr << "Could not create transform settings, please check log for error messages!" << endl;
    return 1;
  }
  Float_t ptmin = 0.1*AliL3Transform::GetSolenoidField();
  Float_t zvertex;

  {
    AliRunLoader *rl = AliRunLoader::Open("galice.root");
    rl->LoadHeader();
    rl->LoadKinematics();
    AliStack* stack = rl->Stack();
    TParticle *orig = (TParticle*)stack->Particle(0);
    Float_t xori = orig->Vx();
    Float_t yori = orig->Vy();
    zvertex = orig->Vz();
    cout<<" Primary vertex at ("<<xori<<","<<yori<<","<<zvertex<<")"<<endl; 
    if (rl->LoadgAlice()) {
      cerr<<"Error occured while loading gAlice"<<endl;
      return 1;
    }
    delete rl;
  }

  cout<<" Hough Tranform will run with ptmin="<<ptmin<<" and zvertex="<<zvertex<<endl;
  
  AliL3Benchmark *fBenchmark = new AliL3Benchmark();
  AliL3Hough *hough = new AliL3Hough();
  hough->SetThreshold(4);
  hough->SetTransformerParams(140,76,ptmin,-1);
  hough->SetPeakThreshold(50,-1);
  hough->Init(path, kFALSE, 100, kFALSE,4,0,0,zvertex);
  hough->SetAddHistograms();

  for(int ev=0; ev<nevent; ev++)
    {
      for(int slice=s1; slice<=s2; slice++)
	{
	  cout<<"Processing slice "<<slice<<endl;
	  hough->ReadData(slice,ev);
	  hough->Transform();
	  hough->AddAllHistogramsRows();
	  hough->FindTrackCandidates();
	  hough->AddTracks();
	}
      hough->WriteTracks(outpath);
      Char_t bname[100];
      sprintf(bname,"rowhough_%d",ev);
      hough->DoBench(bname);

      if(!skip) {
	// Run cluster fitter
	AliL3ClusterFitter *fitter = new AliL3ClusterFitter(path);

	// Set debug flag for the cluster fitter
	//  fitter->Debug();

	// Setting fitter parameters
	fitter->SetInnerWidthFactor(1,1.5);
	fitter->SetOuterWidthFactor(1,1.5);
	fitter->SetNmaxOverlaps(5);
  
	//  fitter->SetChiSqMax(5,kFALSE); //isolated clusters
	fitter->SetChiSqMax(5,kTRUE);  //overlapping clusters

	Int_t rowrange[2] = {0,AliL3Transform::GetNRows()-1};

	// Takes input from global hough tracks produced by HT
	fitter->LoadSeeds(rowrange,kFALSE,ev,zvertex);

	UInt_t ndigits;

	for(int slice=s1; slice<=s2; slice++)
	  {
	    for(Int_t ipatch = 0; ipatch < AliL3Transform::GetNPatches(); ipatch++)
	      {
		// Read digits
		hough->GetMemHandler(ipatch)->Free();
		hough->GetMemHandler(ipatch)->Init(slice,ipatch);
		AliL3DigitRowData *digits = (AliL3DigitRowData *)hough->GetMemHandler(ipatch)->AliAltroDigits2Memory(ndigits,ev);

		fBenchmark->Start("Fitter Init");
		fitter->Init(slice,ipatch);
		fBenchmark->Stop("Fitter Init");
		fitter->SetInputData(digits);
		fBenchmark->Start("Fitter cluster finder");
		fitter->FindClusters();
		fBenchmark->Stop("Fitter cluster finder");
		fitter->WriteClusters();
	      }
	  }

	// Refit of the clusters
	AliL3Vertex vertex;
	//The seeds are the input tracks from circle HT
	AliL3TrackArray *tracks = fitter->GetSeeds();
	AliL3Fitter *ft = new AliL3Fitter(&vertex,1);

	ft->LoadClusters("./fitter/",ev,kFALSE);
	fBenchmark->Start("Track fitter");
	for(Int_t i=0; i<tracks->GetNTracks(); i++)
	  {
	    AliL3Track *track = tracks->GetCheckedTrack(i);
	    if(!track) continue;
	    if(track->GetNHits() < 20) continue;
	    ft->SortTrackClusters(track);
	    ft->FitHelix(track);
	    track->UpdateToFirstPoint();
	  }
	fBenchmark->Stop("Track fitter");
	delete ft;
        
	//Write the final tracks
	fitter->WriteTracks(20);

	delete fitter;
      }
    }

  if(!skip) {
    fBenchmark->Analyze("fitter");
  }

  hough->DoBench("rowhough");
  delete hough;
  delete fBenchmark;
  return 0;
}
