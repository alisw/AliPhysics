#include "AliL3StandardIncludes.h"
#include <sys/time.h>

#include "AliL3Benchmark.h"
#ifdef use_aliroot
#include "AliL3FileHandler.h"
#else
#include "AliL3MemHandler.h"
#endif
#include "AliL3DataHandler.h"
#include "AliL3Transform.h"
#include "AliL3SpacePointData.h"
#include "AliL3DigitData.h"
#include "AliL3Kalman.h"
#include "AliL3Logging.h"
#include "AliL3TrackArray.h"
#include "AliL3Track.h"
#include "AliL3KalmanTrack.h"

#include "TNtuple.h"

/*
  AliL3Kalman
*/

ClassImp(AliL3Kalman)

AliL3Kalman::AliL3Kalman(Char_t *datapath, Int_t *slice, Int_t min_clusters = 0)
{
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
  fRow[0][1] = AliL3Transform::GetLastRow(-1);

  fBenchmark = 0;
}  

AliL3Kalman::~AliL3Kalman()
{
  // Destructor
  if (fBenchmark) delete fBenchmark;
}

void AliL3Kalman::Init()
{
  fBenchmark = new AliL3Benchmark();
}

void AliL3Kalman::LoadTracks(Int_t event, Bool_t sp)
{
  // Load tracks from conformal tracker
  // Must also be possible to take seeds (and clusters) from Hough-transform??

  // Load spacepoints into clusterfile

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Load tracks");

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
	  delete clusterfile[s][p];
        }
    }

  // Load tracks
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
  cout << "Number of loaded tracks " << fTracks->GetNTracks() << endl;
  fBenchmark->Stop("Load tracks");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Kalman::LoadTracks()","Timing")
    <<"Loading tracks in "<<cpuTime*1000<<" ms"<<ENDLOG;

}

void AliL3Kalman::ProcessTracks()
{

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();

  fBenchmark->Start("Process tracks");

  Int_t fEvent = 0;
  fTracks->QSort();

  fKalmanTracks = new AliL3TrackArray();

  // Make a ntuple to store state vector, covariance matrix and chisquare
  TNtuple *kalmanTree = new TNtuple("kalmanTree","kalmantracks","x0:x1:x2:x3:x4:c0:c1:c2:c3:c4:c5:c6:c7:c8:c9:c10:c11:c12:c13:c14:chisq");
  Float_t meas[21];

  // Go through the tracks from conformal or hough tracker
  for (Int_t iTrack = 0; iTrack < fTracks->GetNTracks(); iTrack++)
    {
      AliL3KalmanTrack *kalmantrack = new AliL3KalmanTrack();
      kalmantrack->Init();
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(iTrack);
      if (!track) continue;
      if (track->GetNumberOfPoints() < fMinPointsOnTrack) continue;    

      Bool_t save = kTRUE;

      if (MakeSeed(kalmantrack, track) == 0)
	{
	  save = kFALSE;
	  continue;
	}

      if (Propagate(kalmantrack, track) == 0) 
	{
	  save = kFALSE;
	}

      if (save) {
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
	//cout << "Chisq = " << chisq << endl;
	/*cout << "track " << iTrack << endl; 
	cout << "x = " << meas[0] << ", " << meas[1] << ", " << meas[2] 
	     << ", " << meas[3] << ", " << meas[4] << ", " << endl;
	cout << "c = " << endl 
	     << meas[5] << endl
	     << meas[6] << " " << meas[7] << endl
	     << meas[8] << " " << meas[9] << " " << meas[10] << endl 
	     << meas[11] << " " << meas[12] << " " << meas[13] 
	     << " " << meas[14] << endl
	     << meas[15] << " " << meas[16] << " " << meas[17] << " " 
	     << meas[18] << " " << meas[19] << endl;
	     cout << "chisq = " << meas[20] << endl;*/ 
	//cout << endl;
	// Add the track to the trackarray	
	fKalmanTracks->AddLast(kalmantrack);
	
	// Fill the ntuple with the state vector, covariance matrix and
	// chisquare
	kalmanTree->Fill(meas);
      }

      delete track;
      delete kalmantrack;
    }

  fBenchmark->Stop("Process tracks");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Kalman::ProcessTracks()","Timing")
    <<"Process tracks in "<<cpuTime*1000<<" ms"<<ENDLOG;

  // Write tracks to binary file
  if (fWriteOut)
    {
      Char_t tname[80];
      sprintf(tname,"%s/kalmantracks_%d.raw",fWriteOutPath,fEvent);
      AliL3MemHandler *mem = new AliL3MemHandler();
      mem->SetBinaryOutput(tname);
      mem->TrackArray2Binary(fKalmanTracks);
      mem->CloseBinaryOutput();
      delete mem;
    }
  
  TFile *out = new TFile("out.root","recreate");      
  kalmanTree->Write();
  out->Close();

  delete kalmanTree;
}

Int_t AliL3Kalman::MakeSeed(AliL3KalmanTrack *kalmantrack, AliL3Track *track)
{  
  // Makes a rough state vector and covariance matrix based on the first
  // space point in the outmost row
 
  /*Double_t initTime,cpuTime;
  initTime = GetCpuTime();

  fBenchmark->Start("Make seed");*/

  UInt_t *hitnum = track->GetHitNumbers();
  UInt_t id1, id2, id3;

  // Should do something to make sure that points1 2 and 3 are really not empty
  id1 = hitnum[0];
  Int_t slice1 = (id1>>25) & 0x7f;
  Int_t patch1 = (id1>>22) & 0x7;	
  UInt_t pos1 = id1&0x3fffff;
  AliL3SpacePointData *points1 = fClusters[slice1][patch1];
  
  id2 = hitnum[1];
  Int_t slice2 = (id2>>25) & 0x7f;
  Int_t patch2 = (id2>>22) & 0x7;	
  UInt_t pos2 = id2&0x3fffff;
  AliL3SpacePointData *points2 = fClusters[slice2][patch2];
  if (!points2) return 0;

  id3 = hitnum[2];
  Int_t slice3 = (id3>>25) & 0x7f;
  Int_t patch3 = (id3>>22) & 0x7;	
  UInt_t pos3 = id3&0x3fffff;
  AliL3SpacePointData *points3 = fClusters[slice3][patch3];

  /*fBenchmark->Stop("Make seed");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Kalman::MakeSeed()","Timing")
  <<"Make seeds in "<<cpuTime*1000<<" ms"<<ENDLOG;*/

  return   kalmantrack->MakeTrackSeed(points1,pos1,points2,pos2,points3,pos3);

}

Int_t AliL3Kalman::Propagate(AliL3KalmanTrack *kalmantrack, AliL3Track *track)
{
  // This function should propagte the track to the next layer thera's a 
  // cluster.
  // Must do the following steps,
  // 1. Go to next spacepoint
  // 2. Must find the layer position (old and new), predict new state vector,
  //    and covariance matrix (See AliTPCtrack::PropagateTo
  // 3. Call Update function
  /*Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Propagate");*/

  Int_t num_of_clusters = track->GetNumberOfPoints();
  
  UInt_t *hitnum = track->GetHitNumbers();
  UInt_t id;

  for (Int_t icl = 1; icl < num_of_clusters; icl++)
    {
      id = hitnum[icl];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;	
      UInt_t pos = id&0x3fffff;
      AliL3SpacePointData *points = fClusters[slice][patch];
      if (!points) continue;
      if (kalmantrack->Propagate(points,pos) == 0) return 0;
    }

  /*fBenchmark->Stop("Propagate");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Kalman::Propagate()","Timing")
    <<"Propagate in "<<cpuTime*1000<<" ms"<<ENDLOG;*/

  return 1;
}

void AliL3Kalman::WriteKalmanTrack(AliL3KalmanTrack *kalman)
{
  /*  AliL3MemHandler *memory = new AliL3MemHandler();
  memory->SetBinaryOutput(filename);
  if(opt=='a'||opt=='i'){  //add intracks
    for(Int_t i=0;i<merger->GetNIn();i++){
      AliL3TrackArray *tr=merger->GetInTracks(i);
      memory->TrackArray2Binary(tr);
    }
  }

  if(opt=='o'||opt=='a'){
    AliL3TrackArray *tr=merger->GetOutTracks();
    memory->TrackArray2Binary(tr);
  }

  memory->CloseBinaryOutput();

  return 1;*/
}

void Display()
{
  
}

Double_t AliL3Kalman::GetCpuTime()
{
  //Return the Cputime in seconds.
 struct timeval tv;
 gettimeofday( &tv, NULL );
 return tv.tv_sec+(((Double_t)tv.tv_usec)/1000000.);
 //return (Double_t)(clock()) / CLOCKS_PER_SEC;
}
