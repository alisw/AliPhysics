//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include <string.h>
#include <TCanvas.h>
#include <TFile.h>

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



}


AliL3Hough::AliL3Hough(Int_t n_eta_segments,Int_t xbin,Double_t *xrange,Int_t ybin,Double_t *yrange)
{
  
  fNEtaSegments = n_eta_segments;
  fNxbin = xbin;
  fNybin = ybin;
  fXmin = xrange[0];
  fXmax = xrange[1];
  fYmin = yrange[0];
  fYmax = yrange[1];

  fMemHandler = new AliL3FileHandler();
  fMaxFinder = new AliL3HoughMaxFinder("KappaPhi");
  fEval = new AliL3HoughEval();
  fTransform = new AliL3Transform();
  fDeleteTrack = kTRUE;
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
}


AliL3Hough::~AliL3Hough()
{
  
  if(fHoughTransformer)
    delete fHoughTransformer;
  if(fMemHandler)
    delete fMemHandler;
  if(fMaxFinder)
    delete fMaxFinder;
  if(fEval)
    delete fEval;
  if(fTransform)
    delete fTransform;
}

void AliL3Hough::SetInput(Char_t *input,Bool_t binary)
{
  if(binary)
    {
      strcpy(fPath,input);
      fUseBinary = kTRUE;
    }
  else
    {
      TFile *file = new TFile(input);
      fMemHandler->SetAliInput(file);
      fUseBinary = kFALSE;
    }
  
}

void AliL3Hough::ProcessSlice(Int_t slice)
{
  
}

void AliL3Hough::ProcessPatch(Int_t slice,Int_t patch)
{
  
  Char_t histname[50];
  Int_t i;
  
  if(fHoughTransformer)
    delete fHoughTransformer;
  fHoughTransformer = new AliL3HoughTransformer(slice,patch);//,0,fNEtaSegments);
  
  fHistos = new AliL3Histogram*[fNEtaSegments];
  printf("Allocating %d bytes to histograms\n",fNEtaSegments*sizeof(AliL3Histogram));
  for(i=0; i<fNEtaSegments; i++)
    {
      sprintf(histname,"hist%d",i);
      fHistos[i] = new AliL3Histogram(histname,"",fNxbin,fXmin,fXmax,fNybin,fYmin,fYmax);
    }

  Char_t name[256];
  
  UInt_t ndigits=0;
  AliL3DigitRowData *digits =0;
  // fMemHandler->Init(slice,patch,NRows[patch]);
  //fMemHandler->Init(fTransform);
  if(fUseBinary)
    {
      fMemHandler->Free();
      sprintf(name,"%sdigits_%d_%d.raw",fPath,slice,patch);
      fMemHandler->SetBinaryInput(name);
      digits = (AliL3DigitRowData *)fMemHandler->CompBinary2Memory(ndigits);
      fMemHandler->CloseBinaryInput();
    }
  else
    {
      digits=(AliL3DigitRowData *)fMemHandler->AliDigits2Memory(ndigits); 
    }
  printf("Setting up tables\n");
  fHoughTransformer->SetHistogram(fHistos[0]);
  fHoughTransformer->InitTables();
  fHoughTransformer->SetInputData(ndigits,digits);
  fEval->SetTransformer(fHoughTransformer);
  
  AliL3HoughTrack *track;
  Int_t good_count;
  while(1)
    {
      fHoughTransformer->TransformTables(fHistos);
      
      good_count=0;
      for(Int_t e=0; e<fNEtaSegments; e++)
	{
	  fMaxFinder->SetHistogram(fHistos[e]);
	  track = (AliL3HoughTrack*)fMaxFinder->FindPeak(3,0.95,5);
	  if(fEval->LookInsideRoad(track,e,fDeleteTrack))
	    {
	      //Found a good track here
	      fTracks->AddLast(track);
	      good_count++;
	    }
	}
      break;
      if(good_count==0)
	break;
    }
  printf("good_count %d\n",good_count);
  
}

