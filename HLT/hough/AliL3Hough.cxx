//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 


#include <string.h>

#include "AliL3HoughMerger.h"
#include "AliL3HoughIntMerger.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Logging.h"
#include "AliL3Histogram.h"
#include "AliL3Hough.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughMaxFinder.h"
#ifdef use_aliroot
#include "AliL3FileHandler.h"
#else
#include "AliL3MemHandler.h"
#endif
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
  
  fBinary = kFALSE;
  fNEtaSegments = 0;
  fAddHistograms = kFALSE;
  fDoIterative = kFALSE; 
  fWriteDigits=kFALSE;
  fNPatches=0;
  fMemHandler = 0;
  fHoughTransformer = 0;
  fEval = 0;
  fPeakFinder = 0;
  fTracks = 0;
  fMerger = 0;
  fInterMerger = 0;
  fGlobalMerger = 0;
}


AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments)
{
  //Default ctor.

  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments = n_eta_segments;
  fAddHistograms = kFALSE;
  fDoIterative = kFALSE; 
  fWriteDigits = kFALSE;
  Init();
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

void AliL3Hough::Init()
{
  AliL3Transform::Init(fPath);
  fPeakThreshold = 0;
  fNPatches = AliL3Transform::GetNPatches();
  fHoughTransformer = new AliL3HoughBaseTransformer*[fNPatches];
#ifdef use_aliroot
  fMemHandler = new AliL3FileHandler*[fNPatches];
#else
  fMemHandler = new AliL3MemHandler*[fNPatches];
#endif
  fTracks = new AliL3TrackArray*[fNPatches];
  fEval = new AliL3HoughEval*[fNPatches];
  for(Int_t i=0; i<fNPatches; i++)
    {
      fHoughTransformer[i] = new AliL3HoughTransformer(1,i,fNEtaSegments);
      //fHoughTransformer[i]->CreateHistograms(64,-0.003,0.003,64,-0.26,0.26);
      fHoughTransformer[i]->CreateHistograms(64,0.1,64,-30,30);
      fHoughTransformer[i]->SetThreshold(3);
      fEval[i] = new AliL3HoughEval();
      fTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
#ifdef use_aliroot
      fMemHandler[i] = new AliL3FileHandler();
      if(!fBinary)
	fMemHandler[i]->SetAliInput(fPath);
#else
      fMemHandler[i] = new AliL3MemHandler();
#endif
      
    }
  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi");
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

void AliL3Hough::ReadData(Int_t slice)
{
  //Read data from files, binary or root.

  for(Int_t i=0; i<fNPatches; i++)
    {
      fMemHandler[i]->Free();
      UInt_t ndigits=0;
      AliL3DigitRowData *digits =0;
      Char_t name[256];
      if(fBinary)//take input data from binary files
	{
	  sprintf(name,"%sdigits_%d_%d.raw",fPath,slice,i);
	  fMemHandler[i]->SetBinaryInput(name);
	  digits = (AliL3DigitRowData *)fMemHandler[i]->CompBinary2Memory(ndigits);
	  fMemHandler[i]->CloseBinaryInput();
	}
      else //read data from root file
	{
#ifdef use_aliroot
	  const Int_t rows[2] = {AliL3Transform::GetFirstRow(i),AliL3Transform::GetLastRow(i)};
	  fMemHandler[i]->Init(slice,i,rows);
	  digits=(AliL3DigitRowData *)fMemHandler[i]->AliDigits2Memory(ndigits); 
#else
	  cerr<<"You cannot read from rootfile now"<<endl;
#endif
	}
      fHoughTransformer[i]->SetInputData(ndigits,digits);
    }
}

void AliL3Hough::Transform(Int_t row_range)
{
  //Transform all data given to the transformer within the given slice
  //(after ReadData(slice))

  for(Int_t i=0; i<fNPatches; i++)
    {
      fHoughTransformer[i]->Reset();//Reset the histograms
      if(row_range < 0)
	fHoughTransformer[i]->TransformCircle();
      else
	fHoughTransformer[i]->TransformCircleC(row_range);
    }
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
}

void AliL3Hough::FindTrackCandidates()
{
  //Look for peaks in histograms, and find the track candidates
  
  Int_t n_patches;
  if(fAddHistograms)
    n_patches = 1; //Histograms has been added.
  else
    n_patches = fNPatches;

  
  for(Int_t i=0; i<n_patches; i++)
    {
      AliL3HoughBaseTransformer *tr = fHoughTransformer[i];
      Double_t eta_slice = (tr->GetEtaMax()-tr->GetEtaMin()/tr->GetNEtaSegments());
      fTracks[i]->Reset();
      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->SetHistogram(hist);
	  fPeakFinder->SetThreshold(fPeakThreshold);
	  Int_t n=20;
	  Float_t x[n];
	  Float_t y[n];
	  Int_t weight[n];
	  //fPeakFinder->FindPeak1(x,y,weight,n,2,1);
	  fPeakFinder->FindMaxima(x,y,weight,n);
	  for(Int_t k=0; k<n; k++)
	    {
	      if(weight[k] == 0) continue;
	      
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[i]->NextTrack();
	      track->SetTrackParameters(x[k],y[k],weight[k]);
	      track->SetEtaIndex(j);
	      track->SetEta((Double_t)(j*eta_slice));
	      track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));
	    }
	}
      fTracks[i]->QSort();
    }
}

void AliL3Hough::Evaluate(Int_t road_width)
{
  //Evaluate the tracks, by looking along the road in the raw data.
  

  if(!fTracks[0])
    {
      LOG(AliL3Log::kError,"AliL3Hough::Evaluate","Track Array")
	<<"No tracks to work with..."<<ENDLOG;
      return;
    }
  
  printf("Number of tracks before evaluation %d\n",fTracks[0]->GetNTracks());
  AliL3TrackArray *tracks;
  for(Int_t i=0; i<fNPatches; i++)
    {
      fEval[i]->InitTransformer(fHoughTransformer[i]);
      continue;
      fEval[i]->SetNumOfRowsToMiss(2);
      fEval[i]->SetNumOfPadsToLook(road_width);
      if(fAddHistograms)
	tracks = fTracks[0];
      else
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
	  
	  if(!fEval[i]->LookInsideRoad(track,track->GetEtaIndex()))
	    tracks->Remove(j);
	  if(fAddHistograms)
	    track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));//All rows included
	}
      tracks->Compress();
      tracks->QSort(); //Sort the tracks according to weight
      
      if(!fAddHistograms)
	fMerger->FillTracks(tracks,i); //Copy tracks to the track merger
    }
  
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

void AliL3Hough::WriteTracks(Char_t *path)
{
  AliL3MemHandler *mem = new AliL3MemHandler();
  Char_t fname[100];
  if(fAddHistograms)
    {
      sprintf(fname,"%s/tracks.raw",path);
      mem->SetBinaryOutput(fname);
      mem->TrackArray2Binary(fTracks[0]);
      mem->CloseBinaryOutput();
    }
  else 
    {
      for(Int_t i=0; i<fNPatches; i++)
	{
	  sprintf(fname,"%s/tracks_%d.raw",path,i);
	  mem->SetBinaryOutput(fname);
	  mem->TrackArray2Binary(fTracks[i]);
	  mem->CloseBinaryOutput();
	}
    }
  delete mem;
  
}
#ifdef use_aliroot
void AliL3Hough::WriteDigits(Char_t *outfile)
{
  //Write the current data to a new rootfile.

  for(Int_t i=0; i<fNPatches; i++)
    {
      AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fHoughTransformer[i]->GetDataPointer();
      fMemHandler[i]->AliDigits2RootFile(tempPt,outfile);
    }
  
}
#endif
