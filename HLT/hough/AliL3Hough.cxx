//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 


#include <string.h>
#include <sys/time.h>

#include "AliL3HoughMerger.h"
#include "AliL3HoughIntMerger.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Logging.h"
#include "AliL3Histogram.h"
#include "AliL3Hough.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughTransformerVhdl.h"
#include "AliL3HoughMaxFinder.h"
#ifdef use_aliroot
#include "AliL3FileHandler.h"
#else
#include "AliL3MemHandler.h"
#endif
#include "AliL3DataHandler.h"
#include "AliL3DigitData.h"
#include "AliL3HoughEval.h"
#include "AliL3Transform.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"


//_____________________________________________________________
// AliL3Hough
//
// Interface class for the Hough transform
//
// Example how to use:
//
// AliL3Hough *hough = new AliL3Hough(path,kTRUE,NumberOfEtaSegments);
// hough->ReadData(slice);
// hough->Transform();
// hough->FindTrackCandidates();
// 
// AliL3TrackArray *tracks = hough->GetTracks(patch);

ClassImp(AliL3Hough)

AliL3Hough::AliL3Hough()
{
  //Constructor
  
  fBinary        = kFALSE;
  fAddHistograms = kFALSE;
  fDoIterative   = kFALSE; 
  fWriteDigits   = kFALSE;
  fUse8bits      = kFALSE;

  fMemHandler       = 0;
  fHoughTransformer = 0;
  fEval             = 0;
  fPeakFinder       = 0;
  fTracks           = 0;
  fMerger           = 0;
  fInterMerger      = 0;
  fGlobalMerger     = 0;

  fNEtaSegments     = 0;
  fNPatches         = 0;
  fVersion          = 0;
  fCurrentSlice     = 0;

  SetTransformerParams();
  SetThreshold();
}

AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments,Bool_t bit8=kFALSE,Int_t tv=0)
{
  //Default ctor.

  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments  = n_eta_segments;
  fAddHistograms = kFALSE;
  fDoIterative   = kFALSE; 
  fWriteDigits   = kFALSE;
  fUse8bits      = bit8;
  fVersion       = tv;
}

AliL3Hough::~AliL3Hough()
{
  //dtor

  CleanUp();
  if(fMerger)
    delete fMerger;
  if(fInterMerger)
    delete fInterMerger;
  if(fPeakFinder)
    delete fPeakFinder;
  if(fGlobalMerger)
    delete fGlobalMerger;
}

void AliL3Hough::CleanUp()
{
  //Cleanup memory
  
  for(Int_t i=0; i<fNPatches; i++)
    {
      if(fTracks[i]) delete fTracks[i];
      if(fEval[i]) delete fEval[i];
      if(fHoughTransformer[i]) delete fHoughTransformer[i];
      if(fMemHandler[i]) delete fMemHandler[i];
    }
  
  /*
  if(fTracks) delete [] fTracks;
  if(fEval) delete [] fEval;
  if(fHoughTransformer) delete [] fHoughTransformer;
  if(fMemHandler) delete [] fMemHandler;
  */
}

void AliL3Hough::Init(Char_t *path,Bool_t binary,Int_t n_eta_segments,Bool_t bit8=kFALSE,Int_t tv=0)
{
  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments = n_eta_segments;
  fWriteDigits  = kFALSE;
  fUse8bits     = bit8;
  fVersion      = tv;

  Init(); //do the rest
}

void AliL3Hough::Init(Bool_t doit=kFALSE, Bool_t addhists=kFALSE){
  fDoIterative   = doit; 
  fAddHistograms = addhists;

  AliL3Transform::Init(fPath);
  fNPatches = AliL3Transform::GetNPatches();

  fHoughTransformer = new AliL3HoughBaseTransformer*[fNPatches];
  fMemHandler = new AliL3MemHandler*[fNPatches];
  fTracks = new AliL3TrackArray*[fNPatches];
  fEval = new AliL3HoughEval*[fNPatches];

  for(Int_t i=0; i<fNPatches; i++)
    {
      switch (fVersion){ //choose Transformer
      case 1: 
	fHoughTransformer[i] = new AliL3HoughTransformerVhdl(1,i,fNEtaSegments);
	break;
      default:
	fHoughTransformer[i] = new AliL3HoughTransformer(1,i,fNEtaSegments);
      }

      fHoughTransformer[i]->CreateHistograms(fNBinX,fLowPt,fNBinY,-fPhi,fPhi);
      fHoughTransformer[i]->SetLowerThreshold(fThreshold);

      LOG(AliL3Log::kInformational,"AliL3Hough::Init","Version")
	<<"Initializing Hough transformer version "<<fVersion<<ENDLOG;

      fEval[i] = new AliL3HoughEval();
      fTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
      if(fUse8bits)
	fMemHandler[i] = new AliL3DataHandler();
      else
#ifdef use_aliroot
      	{
	  fMemHandler[i] = new AliL3FileHandler();
	  if(!fBinary)
	    {
	      Char_t filename[100];
	      sprintf(filename,"%s/digitfile",fPath);
	      fMemHandler[i]->SetAliInput(filename);
	    }
	}
#else
      fMemHandler[i] = new AliL3MemHandler();
#endif
    }

  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi",100);
  fMerger = new AliL3HoughMerger(fNPatches);
  fInterMerger = new AliL3HoughIntMerger();
  fGlobalMerger = 0;
}

void AliL3Hough::Process(Int_t minslice,Int_t maxslice)
{
  //Process all slices [minslice,maxslice].
  fGlobalMerger = new AliL3HoughGlobalMerger(minslice,maxslice);
  
  for(Int_t i=minslice; i<=maxslice; i++)
    {
      ReadData(i);
      Transform();
      if(fAddHistograms)
	AddAllHistograms();
      FindTrackCandidates();
      Evaluate();
      fGlobalMerger->FillTracks(fTracks[0],i);
    }
}

void AliL3Hough::ReadData(Int_t slice,Int_t eventnr=0)
{
  //Read data from files, binary or root.
  
  fCurrentSlice = slice;
  for(Int_t i=0; i<fNPatches; i++)
    {
      fMemHandler[i]->Free();
      UInt_t ndigits=0;
      AliL3DigitRowData *digits =0;
      Char_t name[256];
      fMemHandler[i]->Init(slice,i);
      if(fBinary)//take input data from binary files
	{
	  if(fUse8bits)
	    sprintf(name,"%sdigits_c8_%d_%d.raw",fPath,slice,i);
	  else
	    sprintf(name,"%sdigits_%d_%d.raw",fPath,slice,i);
	  fMemHandler[i]->SetBinaryInput(name);
	  digits = (AliL3DigitRowData *)fMemHandler[i]->CompBinary2Memory(ndigits);
	  fMemHandler[i]->CloseBinaryInput();
	}
      else //read data from root file
	{
#ifdef use_aliroot
	  digits=(AliL3DigitRowData *)fMemHandler[i]->AliDigits2Memory(ndigits,eventnr);
	  fMemHandler[i]->FreeDigitsTree();
#else
	  cerr<<"You cannot read from rootfile now"<<endl;
#endif
	}
      fHoughTransformer[i]->SetInputData(ndigits,digits);
    }
}

void AliL3Hough::Transform(Int_t row_range=-1)
{
  //Transform all data given to the transformer within the given slice
  //(after ReadData(slice))
  
  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  for(Int_t i=0; i<fNPatches; i++)
    {
      fHoughTransformer[i]->Reset();//Reset the histograms
      if(row_range < 0)
	fHoughTransformer[i]->TransformCircle();
      else
	fHoughTransformer[i]->TransformCircleC(row_range);
    }
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::Transform()","Timing")
    <<"Transform done in average per patch of "<<cpuTime*1000/fNPatches<<" ms"<<ENDLOG;
}

void AliL3Hough::MergePatches()
{
  if(fAddHistograms) //Nothing to merge here
    return;
  fMerger->MergePatches(kTRUE);
}

void AliL3Hough::MergeInternally()
{
  if(fAddHistograms)
    fInterMerger->FillTracks(fTracks[0]);
  else
    fInterMerger->FillTracks(fMerger->GetOutTracks());
  
  fInterMerger->MMerge();
}

void AliL3Hough::ProcessSliceIter()
{
  //Process current slice (after ReadData(slice)) iteratively.
  
  for(Int_t i=0; i<fNPatches; i++)
    {
      ProcessPatchIter(i);
      fMerger->FillTracks(fTracks[i],i); //Copy tracks to merger
    }
}

void AliL3Hough::ProcessPatchIter(Int_t patch)
{
  //Process patch in a iterative way. 
  //transform + peakfinding + evaluation + transform +...

  Int_t num_of_tries = 10;
  AliL3HoughBaseTransformer *tr = fHoughTransformer[patch];
  AliL3TrackArray *tracks = fTracks[patch];
  tracks->Reset();
  AliL3HoughEval *ev = fEval[patch];
  ev->InitTransformer(tr);
  ev->RemoveFoundTracks();
  ev->SetNumOfRowsToMiss(2);
  ev->SetNumOfPadsToLook(2);
  AliL3Histogram *hist;
  for(Int_t t=0; t<num_of_tries; t++)
    {
      tr->Reset();
      tr->TransformCircle();
      for(Int_t i=0; i<fNEtaSegments; i++)
	{
	  hist = tr->GetHistogram(i);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->SetHistogram(hist);
	  //Int_t n=1;
	  Float_t x,y;
	  //fPeakFinder->FindAbsMaxima(*x,*y);
	  fPeakFinder->FindPeak(3,0.95,5,x,y);
	  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();
	  track->SetTrackParameters(x,y,1);
	  if(!ev->LookInsideRoad(track,i))
	    {	
	      tracks->Remove(tracks->GetNTracks()-1);
	      tracks->Compress();
	    }
	}
    }
  LOG(AliL3Log::kInformational,"AliL3Hough::ProcessPatch","NTracks")
    <<AliL3Log::kDec<<"Found "<<tracks->GetNTracks()<<" tracks in patch "<<patch<<ENDLOG;
}

void AliL3Hough::AddAllHistograms()
{
  //Add the histograms within one etaslice.
  //Resulting histogram are in patch=0.

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(i);
      for(Int_t j=1; j<fNPatches; j++)
	{
	  AliL3Histogram *hist = fHoughTransformer[j]->GetHistogram(i);
	  hist0->Add(hist);
	}
    }
  fAddHistograms = kTRUE;
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::AddAllHistograms()","Timing")
    <<"Adding histograms in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliL3Hough::FindTrackCandidates()
{
  //Look for peaks in histograms, and find the track candidates
  
  Int_t n_patches;
  if(fAddHistograms)
    n_patches = 1; //Histograms have been added.
  else
    n_patches = fNPatches;
  
  Double_t initTime,cpuTime;
  initTime = GetCpuTime();

  for(Int_t i=0; i<n_patches; i++)
    {
      AliL3HoughBaseTransformer *tr = fHoughTransformer[i];
      Double_t eta_slice = (tr->GetEtaMax()-tr->GetEtaMin())/tr->GetNEtaSegments();
      fTracks[i]->Reset();

      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->Reset();
	  fPeakFinder->SetHistogram(hist);
	  fPeakFinder->FindMaxima(0,0); //Simple maxima finder
	  //fPeakFinder->FindAbsMaxima();

	  for(Int_t k=0; k<fPeakFinder->GetEntries(); k++)
	    {
	      if(fPeakFinder->GetWeight(k) == 0) continue;
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[i]->NextTrack();
	      track->SetTrackParameters(fPeakFinder->GetXPeak(k),fPeakFinder->GetYPeak(k),fPeakFinder->GetWeight(k));
	      track->SetEtaIndex(j);
	      Double_t eta = (Double_t)((j+0.5)*eta_slice);
	      if(fCurrentSlice > 17) eta*=-1;
	      track->SetEta(eta);
	      track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));
	    }
	}
      fTracks[i]->QSort();
    }
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::FindTrackCandidates()","Timing")
    <<"Maxima finding done in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliL3Hough::InitEvaluate()
{
  //Pass the transformer objects to the AliL3HoughEval objects:
  //This will provide the evaluation objects with all the necessary
  //data and parameters it needs.
  
  for(Int_t i=0; i<fNPatches; i++) 
    fEval[i]->InitTransformer(fHoughTransformer[i]);
}

Int_t AliL3Hough::Evaluate(Int_t road_width,Int_t nrowstomiss)
{
  //Evaluate the tracks, by looking along the road in the raw data.
  //If track does not cross all padrows - rows2miss, it is removed from the arrray.
  //If histograms were not added, the check is done locally in patch,
  //meaning that nrowstomiss is the number of padrows the road can miss with respect
  //to the number of rows in the patch.
  //If the histograms were added, the comparison is done globally in the _slice_, 
  //meaing that nrowstomiss is the number of padrows the road can miss with
  //respect to the total number of padrows in the slice.
  //
  //Return value = number of tracks which were removed (only in case of fAddHistograms)
  
  if(!fTracks[0])
    {
      LOG(AliL3Log::kError,"AliL3Hough::Evaluate","Track Array")
	<<"No tracks to work with..."<<ENDLOG;
      return 0;
    }
  
  InitEvaluate();
  
  Int_t removed_tracks=0;
  AliL3TrackArray *tracks=0;
  Int_t *total_rows=0;
  if(fAddHistograms)
    {    
      tracks = fTracks[0];
      total_rows = new Int_t[tracks->GetNTracks()];
      for(Int_t i=0; i<tracks->GetNTracks(); i++)
	total_rows[i]=0;
    }

  for(Int_t i=0; i<fNPatches; i++)
    {
      fEval[i]->InitTransformer(fHoughTransformer[i]);
      fEval[i]->SetNumOfPadsToLook(road_width);
      fEval[i]->SetNumOfRowsToMiss(nrowstomiss);
      if(!fAddHistograms)
	tracks = fTracks[i];
      
      for(Int_t j=0; j<tracks->GetNTracks(); j++)
	{
	  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  if(!track)
	    {
	      LOG(AliL3Log::kWarning,"AliL3Hough::Evaluate","Track array")
		<<"Track object missing!"<<ENDLOG;
	      continue;
	    }
	  
	  Bool_t result = fEval[i]->LookInsideRoad(track,total_rows[j]);
	  if(!fAddHistograms)//the track crossed too few good padrows (padrows with signal) in the patch, so remove it
	    {
	      if(result == kFALSE)
		tracks->Remove(j);
	    }
	}
      if(!fAddHistograms)
	{
	  tracks->Compress();
	  tracks->QSort(); 
	  fMerger->FillTracks(tracks,i); //Copy tracks to the track merger
	}
    }
  
  if(fAddHistograms) //Here we check the tracks globally; how many good rows (padrows with signal) did it cross in the slice
    {
      for(Int_t j=0; j<tracks->GetNTracks(); j++)
	{
	  if(total_rows[j] < AliL3Transform::GetNRows() - nrowstomiss)
	    {
	      tracks->Remove(j);
	      removed_tracks++;
	    }
	}
      tracks->Compress();
      tracks->QSort();
    }
  
  if(total_rows)
    delete [] total_rows;
  
  return removed_tracks;
}

void AliL3Hough::EvaluateWithEta()
{

  
  if(!fTracks[0])
    {
      printf("AliL3Hough::EvaluateWithEta: NO TRACKS\n");
      return;
    }
  printf("Number of tracks before evaluation %d\n",fTracks[0]->GetNTracks());
 
  for(Int_t i=0; i<fNPatches; i++)
    {
      fEval[i]->InitTransformer(fHoughTransformer[i]);
      fEval[i]->FindEta(fTracks[0]);
    }
  fMerger->FillTracks(fTracks[0],0);
}

void AliL3Hough::WriteTracks(Int_t slice,Char_t *path)
{
  //Write the tracks in slice
  
  AliL3MemHandler *mem = new AliL3MemHandler();
  Char_t fname[100];
  if(fAddHistograms)
    {
      sprintf(fname,"%s/tracks_ho_%d.raw",path,slice);
      mem->SetBinaryOutput(fname);
      mem->TrackArray2Binary(fTracks[0]);
      mem->CloseBinaryOutput();
    }
  else 
    {
      for(Int_t i=0; i<fNPatches; i++)
	{
	  sprintf(fname,"%s/tracks_ho_%d_%d.raw",path,slice,i);
	  mem->SetBinaryOutput(fname);
	  mem->TrackArray2Binary(fTracks[i]);
	  mem->CloseBinaryOutput();
	}
    }
  delete mem;
  
}

void AliL3Hough::WriteDigits(Char_t *outfile)
{
#ifdef use_aliroot  
  //Write the current data to a new rootfile.

  for(Int_t i=0; i<fNPatches; i++)
    {
      AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fHoughTransformer[i]->GetDataPointer();
      fMemHandler[i]->AliDigits2RootFile(tempPt,outfile);
    }
#else
  cerr<<"AliL3Hough::WriteDigits : You need to compile with AliROOT!"<<endl;
  return;
#endif  
}

Double_t AliL3Hough::GetCpuTime()
{
  //Return the Cputime in seconds.
 struct timeval tv;
 gettimeofday( &tv, NULL );
 return tv.tv_sec+(((Double_t)tv.tv_usec)/1000000.);
 //return (Double_t)(clock()) / CLOCKS_PER_SEC;
}

