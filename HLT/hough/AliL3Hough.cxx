//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include <string.h>
#include <TCanvas.h>
#include <TFile.h>

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

ClassImp(AliL3Hough)

AliL3Hough::AliL3Hough()
{
  fBinary = kFALSE;
  fNEtaSegments = 0;
  fAddHistograms = kFALSE;
  fRemoveFoundTracks = kFALSE; 
  fWriteDigits=kFALSE;
}


AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments)
{
  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments = n_eta_segments;
  fAddHistograms = kFALSE;
  fRemoveFoundTracks = kFALSE; 
  fWriteDigits = kFALSE;
  Init();
}


AliL3Hough::~AliL3Hough()
{
  if(fMemHandler)
    DeleteMemory();
  if(fHoughTransformer)
    DeleteTransformers();
  if(fEval)
    DeleteEval();
  if(fPeakFinder)
    delete fPeakFinder;
  if(fTracks)
    delete fTracks;
  if(fRootFile)
    {
      fRootFile->Close();
      delete fRootFile;
    }
}

void AliL3Hough::DeleteEval()
{
  for(Int_t i=0; i<NPatches; i++)
    {
      if(!fEval[i]) continue;
      delete fEval[i];
    }
  delete [] fEval;
}

void AliL3Hough::DeleteTransformers()
{
  for(Int_t i=0; i<NPatches; i++)
    {
      if(!fHoughTransformer[i]) continue;
      delete fHoughTransformer[i];
    }
  delete [] fHoughTransformer;
}

void AliL3Hough::DeleteMemory()
{
  for(Int_t i=0; i<NPatches; i++)
    {
      if(!fMemHandler[i]) continue;
      delete fMemHandler[i];
    }
  delete [] fMemHandler;
}

void AliL3Hough::Init()
{
  fHoughTransformer = new AliL3HoughTransformer*[NPatches];
  fMemHandler = new AliL3FileHandler*[NPatches];
  for(Int_t i=0; i<NPatches; i++)
    {
      fHoughTransformer[i] = new AliL3HoughTransformer(1,i,fNEtaSegments);
      fHoughTransformer[i]->CreateHistograms(64,-0.003,0.003,64,-0.26,0.26);
      fHoughTransformer[i]->SetThreshold(3);
      fMemHandler[i] = new AliL3FileHandler();
      if(!fBinary)
	fMemHandler[i]->SetAliInput(fPath);
    }
  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi");
}

void AliL3Hough::Process(Int_t minslice,Int_t maxslice)
{
  //Process all slices [minslice,maxslice].

  for(Int_t i=minslice; i<=maxslice; i++)
    {
      TransformSlice(i);
      if(fAddHistograms)
	AddAllHistograms();
      FindTrackCandidates();
      Evaluate(fRemoveFoundTracks);
      if(fWriteDigits)
	WriteDigits();
    }
}

void AliL3Hough::TransformSlice(Int_t slice)
{
  
  for(Int_t i=0; i<NPatches; i++)
    {
      //Reset memories
      fHoughTransformer[i]->Reset();      
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
      fHoughTransformer[i]->TransformCircle();
    }
  
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
  
  LOG(AliL3Log::kDebug,"AliL3Hough::AddAllHistograms","Progress")
    <<"Adding all histograms"<<ENDLOG;
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(i);
      for(Int_t j=1; j<NPatches; j++)
	{
	  AliL3Histogram *hist = fHoughTransformer[j]->GetHistogram(i);
	  hist0->Add(hist);
	}
    }
}

void AliL3Hough::FindTrackCandidates()
{
  //Look for peaks in histograms, and find the track candidates
  
  if(fTracks)
  {
    LOG(AliL3Log::kDebug,"AliL3Hough::FindTrackCandidates","Track array")
      <<"Deleting old track array"<<ENDLOG;
    delete fTracks;
  }
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
  
  Int_t n_patches;
  if(fAddHistograms)
    n_patches = 1; //Histograms has been added.
  else
    n_patches = NPatches;

  for(Int_t i=0; i<n_patches; i++)
    {
      AliL3HoughTransformer *tr = fHoughTransformer[i];
      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  fPeakFinder->SetHistogram(hist);
	  Int_t n=10;
	  Float_t x[10];
	  Float_t y[10];
	  fPeakFinder->FindPeak1(x,y,n);
	  for(Int_t k=0; k<n; k++)
	    {
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks->NextTrack();
	      track->SetTrackParameters(x[k],y[k],1);
	      track->SetEtaIndex(j);
	    }
	}
    }
  
}

void AliL3Hough::Evaluate(Bool_t remove)
{
  //Evaluate the tracks, by looking along the road in the raw data.
  //You may choose to remove the found tracks from the image.
  
  if(!fTracks)
    {
      LOG(AliL3Log::kError,"AliL3Hough::Evaluate","Track array")
	<<AliL3Log::kHex<<"No tracks to work on "<<(Int_t)fTracks<<ENDLOG;
      return;
    }
  
  if(fEval)
    {
      LOG(AliL3Log::kDebug,"AliL3Hough::Evaluate","Evaluate object")
	<<"Deleting old AliL3HoughEval objects"<<ENDLOG;
      DeleteEval();
    }
  
  fEval = new AliL3HoughEval*[NPatches];
  for(Int_t i=0; i<NPatches; i++)
    {
      fEval[i] = new AliL3HoughEval(fHoughTransformer[i]);
      if(remove)
	fEval[i]->RemoveFoundTracks();
      for(Int_t j=0; j<fTracks->GetNTracks(); j++)
	{
	  AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks->GetCheckedTrack(j);
	  if(!track)
	    {
	      printf("AliL3Hough::Evaluate : Missing track object...\n");
	      continue;
	    }
	  //fEval[i]->LookInsideRoad(track,track->GetEtaIndex());
	  if(!fEval[i]->LookInsideRoad(track,track->GetEtaIndex()))
	    fTracks->Remove(j);
	}
      fTracks->Compress();
    }
  
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
