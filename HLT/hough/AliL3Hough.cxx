//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 


#include <string.h>
#include <TCanvas.h>
#include <TFile.h>

#include "AliL3HoughMerger.h"
#include "AliL3HoughIntMerger.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Logging.h"
#include "AliL3Histogram.h"
#include "AliL3Hough.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughMaxFinder.h"
#include "AliL3FileHandler.h"
#include "AliL3DigitData.h"
#include "AliL3HoughEval.h"
#include "AliL3Transform.h"
#include "AliL3Defs.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3Benchmark.h"

//_____________________________________________________________
// AliL3Hough
//
// Base class for the Hough transform
//


ClassImp(AliL3Hough)

AliL3Hough::AliL3Hough()
{
  //Constructor
  
  fBinary = kFALSE;
  fNEtaSegments = 0;
  fAddHistograms = kFALSE;
  fDoIterative = kFALSE; 
  fWriteDigits=kFALSE;
}


AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments)
{
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
  CleanUp();
  if(fMerger)
    delete fMerger;
  if(fInterMerger)
    delete fInterMerger;
  if(fPeakFinder)
    delete fPeakFinder;
}

void AliL3Hough::CleanUp()
{
  //Cleanup memory
  
  for(Int_t i=0; i<NPatches; i++)
    {
      if(fTracks[i]) delete fTracks[i];
      if(fEval[i]) delete fEval[i];
      if(fHoughTransformer[i]) delete fHoughTransformer[i];
      if(fMemHandler[i]) delete fMemHandler[i];
    }
  
  /*Shitty compiler doesn't allow this:
    
  if(fTracks) delete [] fTracks;
  if(fEval) delete [] fEval;
  if(fHoughTransformer) delete [] fHoughTransformer;
  if(fMemHandler) delete [] fMemHandler;
  */
}

void AliL3Hough::Init()
{
  fHoughTransformer = new AliL3HoughTransformer*[NPatches];
  fMemHandler = new AliL3FileHandler*[NPatches];
  fTracks = new AliL3TrackArray*[NPatches];
  fEval = new AliL3HoughEval*[NPatches];
  for(Int_t i=0; i<NPatches; i++)
    {
      fHoughTransformer[i] = new AliL3HoughTransformer(1,i,fNEtaSegments);
      fHoughTransformer[i]->CreateHistograms(64,-0.003,0.003,64,-0.26,0.26);
      fHoughTransformer[i]->SetThreshold(3);
      fEval[i] = new AliL3HoughEval();
      fTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
      fMemHandler[i] = new AliL3FileHandler();
      if(!fBinary)
	fMemHandler[i]->SetAliInput(fPath);
    }
  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi");
  fMerger = new AliL3HoughMerger(NPatches);
  fInterMerger = new AliL3HoughIntMerger();
}

void AliL3Hough::Process(Int_t minslice,Int_t maxslice)
{
  //Process all slices [minslice,maxslice].
  
  for(Int_t i=minslice; i<=maxslice; i++)
    {
      ReadData(i);
      Transform();
      if(fAddHistograms)
	AddAllHistograms();
      FindTrackCandidates();
      Evaluate();
      if(fWriteDigits)
	WriteDigits();
    }
}

void AliL3Hough::ReadData(Int_t slice)
{
  //Read data from files, binary or root.

  for(Int_t i=0; i<NPatches; i++)
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
	  fMemHandler[i]->Init(slice,i,NRows[i]);
	  digits=(AliL3DigitRowData *)fMemHandler[i]->AliDigits2Memory(ndigits); 
	}
      fHoughTransformer[i]->SetInputData(ndigits,digits);
    }
}

void AliL3Hough::Transform()
{
  //Transform all data given to the transformer within the given slice
  //(after ReadData(slice))

  Double_t initTime,finalTime;
  for(Int_t i=0; i<NPatches; i++)
    {
      fHoughTransformer[i]->Reset();//Reset the histograms
      initTime = AliL3Benchmark::GetCpuTime();
      fHoughTransformer[i]->TransformCircle();
      finalTime = AliL3Benchmark::GetCpuTime();
      LOG(AliL3Log::kInformational,"AliL3Hough::Transform","Timing")
	<<AliL3Log::kDec<<"Transform finished in "<<(finalTime-initTime)*1000<<"ms"<<ENDLOG;
    }
}

void AliL3Hough::MergePatches()
{
  if(fAddHistograms) //Nothing to merge here
    return;
  AliL3Transform *tr = new AliL3Transform();
  fMerger->SetTransformer(tr);
  fMerger->MergePatches(kTRUE);
  delete tr;
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
  
  for(Int_t i=0; i<NPatches; i++)
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
  AliL3HoughTransformer *tr = fHoughTransformer[patch];
  AliL3TrackArray *tracks = fTracks[patch];
  tracks->Reset();
  AliL3HoughEval *ev = fEval[patch];
  ev->InitTransformer(tr);
  ev->RemoveFoundTracks();
  ev->SetNumOfRowsToMiss(2);
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
	  Int_t n=1;
	  Int_t x[n],y[n];
	  fPeakFinder->FindAbsMaxima(*x,*y);
	  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();
	  track->SetTrackParameters(hist->GetBinCenterX(*x),hist->GetBinCenterY(*y),1);
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

AliL3Histogram *AliL3Hough::AddHistograms(Int_t eta_index)
{

  AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(eta_index);
  for(Int_t i=1; i<NPatches; i++)
    {
      AliL3Histogram *hist = fHoughTransformer[i]->GetHistogram(eta_index);
      hist0->Add(hist);
    }
  
  return hist0;
}

void AliL3Hough::AddAllHistograms()
{
  //Add the histograms within one etaslice.
  //Resulting histogram are in patch=0.

  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(i);
      for(Int_t j=1; j<NPatches; j++)
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
    n_patches = NPatches;
  
  for(Int_t i=0; i<n_patches; i++)
    {
      AliL3HoughTransformer *tr = fHoughTransformer[i];
      fTracks[i]->Reset();
      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->SetHistogram(hist);
	  Int_t n=10;
	  Float_t x[10];
	  Float_t y[10];
	  Int_t weight[10];
	  fPeakFinder->FindPeak1(x,y,weight,n);
	  for(Int_t k=0; k<n; k++)
	    {
	      if(weight[k] == 0) continue;
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[i]->NextTrack();
	      track->SetTrackParameters(x[k],y[k],weight[k]);
	      track->SetEtaIndex(j);
	      track->SetEta((Double_t)(j*tr->GetEtaSlice()));
	      track->SetRowRange(NRows[0][0],NRows[5][1]);
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
  for(Int_t i=0; i<NPatches; i++)
    {
      fEval[i]->InitTransformer(fHoughTransformer[i]);
      fEval[i]->SetNumOfRowsToMiss(1);
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
	    track->SetRowRange(NRows[0][0],NRows[5][1]);//All rows included
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
 
  for(Int_t i=0; i<NPatches; i++)
    {
      fEval[i]->InitTransformer(fHoughTransformer[i]);
      fEval[i]->FindEta(fTracks[0]);
    }
  fMerger->FillTracks(fTracks[0],0);
  fMerger->MergeEtaSlices(0);
}

void AliL3Hough::WriteDigits(Char_t *outfile)
{
  //Write the current data to a new rootfile.

  for(Int_t i=0; i<NPatches; i++)
    {
      AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fHoughTransformer[i]->GetDataPointer();
      fMemHandler[i]->AliDigits2RootFile(tempPt,outfile);
    }
  
}
