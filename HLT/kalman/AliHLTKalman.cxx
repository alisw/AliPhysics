// $Id$

#include "AliHLTStandardIncludes.h"
#include "AliHLTRootTypes.h"
#include <sys/time.h>
#include <TNtuple.h>
#include <TTimer.h>
#include "AliHLTTrack.h"
#include "AliHLTKalmanTrack.h"
#include "AliHLTBenchmark.h"
#include "AliHLTMemHandler.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDataHandler.h"
#include "AliHLTTransform.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTDigitData.h"
#include "AliHLTLogging.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTInterMerger.h"
#include "AliHLTTrackMerger.h"
#include "AliHLTKalman.h"

/*
  AliHLTKalman
*/

ClassImp(AliHLTKalman)

AliHLTKalman::AliHLTKalman(Char_t *datapath, Int_t *slice, Int_t min_clusters = 0){
  // Constructor
  if (slice)
    {
      fMinSlice = slice[0];
      fMaxSlice = slice[1];
    }
  else
    {
      fMinSlice = 0;
      fMaxSlice = 35;
    }
  
  sprintf(fPath,"%s",datapath);
  fMinPointsOnTrack = 0; //min_clusters;
  fTracks = 0;
  fKalmanTracks = 0;

  // NB! fNrow under only for single-patch, must include other possibilities 
  // later on. ?? Maybe better also to put it in an Init-function
  fRow[0][0] = 0;
  fRow[0][1] = AliHLTTransform::GetLastRow(-1);
  fWriteOut = kTRUE;
  fBenchmark = 0;

}  

AliHLTKalman::~AliHLTKalman()
{
  // Destructor
  if (fBenchmark) delete fBenchmark;
}

void AliHLTKalman::Init()
{
  fBenchmark = new AliHLTBenchmark();
}

void AliHLTKalman::LoadTracks(Int_t event, Bool_t sp)
{
  // Load space points and tracks from conformal tracker
  // Must also be possible to take seeds (and clusters) from Hough-transform??

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Load tracks");

  // Load space points
  Char_t fname[1024];
  AliHLTFileHandler *clusterfile[36][6];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliHLTTransform::GetNPatches(); p++)
        {
          Int_t patch;
          if(sp==kTRUE)
            patch=-1;
          else
            patch=p;
	  
          fClusters[s][p] = 0;
          clusterfile[s][p] = new AliHLTFileHandler();
          if(event<0)
            sprintf(fname,"%s/points_%d_%d.raw",fPath,s,patch);
          else
            sprintf(fname,"%s/points_%d_%d_%d.raw",fPath,event,s,patch);
          if(!clusterfile[s][p]->SetBinaryInput(fname))
            {
              LOG(AliHLTLog::kError,"AliHLTEvaluation::Setup","File Open")
                <<"Inputfile "<<fname<<" does not exist"<<ENDLOG;
              delete clusterfile[s][p];
              clusterfile[s][p] = 0;
              continue;
            }
          fClusters[s][p] = (AliHLTSpacePointData*)clusterfile[s][p]->Allocate();
          clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
          clusterfile[s][p]->CloseBinaryInput();
	  if(sp==kTRUE)
            break;
	  delete clusterfile[s][p];
        }
    }

  // Load tracks
  //sprintf(fname,"%s/kalmantracks_%d.raw",fPath,event);
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
  delete tfile;
  fBenchmark->Stop("Load tracks");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliHLTLog::kInformational,"AliHLTKalman::LoadTracks()","Timing")
    <<"Loading tracks in "<<cpuTime*1000<<" ms"<<ENDLOG;

}

void AliHLTKalman::ProcessTracks()
{
  // Run the Kalman filter algorithm on the loaded tracks. 
  // If the track is OK, the loaded track is saved in file kalmantracks_0.raw
  // The kalman filter variables (that is the state vector, covariance matrix 
  // and chi2) is written to root-file kalmantracks.root. 
  // Should be extended to correct the track parameters of the loaded track
  // or save the kalmantrack instead of the loaded track.??
  UInt_t fEvent = 0;

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();

  fBenchmark->Start("Process tracks");

  fTracks->QSort();

  fKalmanTracks = new AliHLTTrackArray();

  // Make a ntuple to store state vector, covariance matrix and chisquare
  // Will eventually not need a TTree??
  TNtuple *kalmanTree = new TNtuple("kalmanTree","kalmantracks","x0:x1:x2:x3:x4:c0:c1:c2:c3:c4:c5:c6:c7:c8:c9:c10:c11:c12:c13:c14:chisq");
  Float_t meas[21];

  // Go through the tracks from conformal or hough tracker
  for (Int_t iTrack = 0; iTrack < fTracks->GetNTracks(); iTrack++)
    {
      /*Double_t initTime,cpuTime;
      initTime = GetCpuTime();
      fBenchmark->Start("Process tracks");*/

      AliHLTTrack *track = (AliHLTTrack*)fTracks->GetCheckedTrack(iTrack);
      if (!track) continue;
      if (track->GetNumberOfPoints() < fMinPointsOnTrack) continue;    

      AliHLTKalmanTrack *kalmantrack = new AliHLTKalmanTrack();

      Bool_t save = kTRUE;

      /*if (InitKalmanTrack(kalmantrack, track) == 0)
	{
	  save = kFALSE;
	  continue;
	  }*/
      
      if (MakeKalmanSeed(kalmantrack,track) ==0)
	{
	  save = kFALSE;
	  continue;
	}

      if (Propagate(kalmantrack, track) == 0) 
	{
	  save = kFALSE;
	}

      if (save) {// cout << track->GetPt() << endl;
	Float_t x[5]; 
	kalmantrack->GetStateVector(x);
	Float_t c[15]; 
	kalmantrack->GetCovariance(c);
	Float_t chisq = kalmantrack->GetChisq();
	meas[0]  = x[0];
	meas[1]  = x[1];
	meas[2]  = x[2];
	meas[3]  = x[3];
	meas[4]  = x[4];
	meas[5]  = c[0];
	meas[6]  = c[1];
	meas[7]  = c[2];
	meas[8]  = c[3];
	meas[9]  = c[4];
	meas[10] = c[5];
	meas[11] = c[6];
	meas[12] = c[7];
	meas[13] = c[8];
	meas[14] = c[9];
	meas[15] = c[10];
	meas[16] = c[11];
	meas[17] = c[12];
	meas[18] = c[13];
	meas[19] = c[14];
	meas[20] = chisq;

	// Add the track to the trackarray	
	AliHLTTrack *outtrack = (AliHLTTrack*)fKalmanTracks->NextTrack();
	outtrack->Set(track);
	// SET THE PARAMETERS ACCORDING TO KALMAN FILTER
	outtrack->SetTgl(x[3]);
	// The factor 2 in the expression for Pt is not included in the similar offline expression. However
	// it should be like this if I use a factor 1/2 in the calculation of par4??
	//outtrack->SetPt(1/(2*TMath::Abs(1e-9*TMath::Abs(x[4])/x[4] + x[4])/(0.0029980*AliHLTTransform::GetBField())));
	//outtrack->SetPt(1/(TMath::Abs(x[4])/0.0029980*AliHLTTransform::GetBField()));
	outtrack->SetPsi(x[2]);
	//outtrack->Set(track);

	// Fill the ntuple with the state vector, covariance matrix and
	// chisquare
	kalmanTree->Fill(meas);
      }

      delete track;
      delete kalmantrack;

      /*fBenchmark->Stop("Process tracks");
      cpuTime = GetCpuTime() - initTime;
      LOG(AliHLTLog::kInformational,"AliHLTKalman::ProcessTracks()","Timing")
      <<"Processed track "<<iTrack<<" in "<<cpuTime*1000<<" ms"<<ENDLOG;*/
    }

  fBenchmark->Stop("Process tracks");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliHLTLog::kInformational,"AliHLTKalman::ProcessTracks()","Timing")
    <<"Process tracks in "<<cpuTime*1000<<" ms"<<ENDLOG;
  
  if (fWriteOut)
    {
      Char_t tname[80];
      sprintf(tname,"%s/kalmantracks_%d.raw",fWriteOutPath,fEvent);
      AliHLTMemHandler *mem = new AliHLTMemHandler();
      mem->SetBinaryOutput(tname);
      mem->TrackArray2Binary(fKalmanTracks);
      mem->CloseBinaryOutput();
      delete mem;
    }
  
  // This will be removed??
  TFile *out = new TFile("kalmantracks.root","recreate");      
  kalmanTree->Write();
  out->Close();

  delete kalmanTree;
  
}

Int_t AliHLTKalman::MakeKalmanSeed(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track)
{
  Int_t num_of_clusters = track->GetNumberOfPoints();

  UInt_t *hitnum = track->GetHitNumbers();
  UInt_t id;

  id = hitnum[0];
  Int_t slice0 = (id>>25) & 0x7f;
  Int_t patch0 = (id>>22) & 0x7;	
  UInt_t pos0 = id&0x3fffff;
  AliHLTSpacePointData *points0 = fClusters[slice0][patch0];

  id = hitnum[Int_t(num_of_clusters/2)];
  Int_t slice1 = (id>>25) & 0x7f;
  Int_t patch1 = (id>>22) & 0x7;	
  UInt_t pos1 = id&0x3fffff;
  AliHLTSpacePointData *points1 = fClusters[slice1][patch1];

  id = hitnum[num_of_clusters-1];
  Int_t slice2 = (id>>25) & 0x7f;
  Int_t patch2 = (id>>22) & 0x7;	
  UInt_t pos2 = id&0x3fffff;
  AliHLTSpacePointData *points2 = fClusters[slice2][patch2];

  return kalmantrack->MakeSeed(track, points0, pos0, slice0, points1, pos1, slice1, points2, pos2, slice2);
}

Int_t AliHLTKalman::InitKalmanTrack(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track)
{
  UInt_t *hitnum = track->GetHitNumbers();
  UInt_t id;

  id = hitnum[0];
  Int_t slice = (id>>25) & 0x7f;
  Int_t patch = (id>>22) & 0x7;	
  UInt_t pos = id&0x3fffff;
  AliHLTSpacePointData *points = fClusters[slice][patch];

  return kalmantrack->Init(track, points, pos, slice);
}

Int_t AliHLTKalman::Propagate(AliHLTKalmanTrack *kalmantrack, AliHLTTrack *track)
{
  // This function propagtes the kalmantrack to the next cluster of the loaded 
  // track 
  Int_t num_of_clusters = track->GetNumberOfPoints();
  UInt_t *hitnum = track->GetHitNumbers();
  UInt_t id;
  UInt_t badpoint = 0;

  for (Int_t icl = 1; icl < num_of_clusters; icl++)
    {

      id = hitnum[icl];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;

      AliHLTSpacePointData *points = fClusters[slice][patch];
      if (!points) continue;
      if (kalmantrack->Propagate(points,pos,slice) == 0) 
	{
	  badpoint++;
	  continue;
	}
    }
  
  // If too many clusters are missing, the track is no good
  if (badpoint >= UInt_t(num_of_clusters*0.8)) return 0;
  return 1;
}

Double_t AliHLTKalman::GetCpuTime()
{
  //Return the Cputime in seconds.
 struct timeval tv;
 gettimeofday( &tv, NULL );
 return tv.tv_sec+(((Double_t)tv.tv_usec)/1000000.);
 //return (Double_t)(clock()) / CLOCKS_PER_SEC;
}
