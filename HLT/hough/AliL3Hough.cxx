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


AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t n_eta_segments)
{
  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments = n_eta_segments;
  Init();
}


AliL3Hough::~AliL3Hough()
{
  if(fMemHandler)
    DeleteMemory();
  if(fHoughTransformer)
    DeleteTransformers();
  if(fRootFile)
    {
      fRootFile->Close();
      delete fRootFile;
    }
}

void AliL3Hough::DeleteTransformers()
{
  for(Int_t i=0; i<NPatches; i++)
    {
      if(!fHoughTransformer[i]) continue;
      delete fHoughTransformer;
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
      fMemHandler[i] = new AliL3FileHandler();
    }
  if(!fBinary)
    fRootFile = new TFile(fPath);
}

void AliL3Hough::TransformSlice(Int_t slice)
{
  
  for(Int_t i=0; i<NPatches; i++)
    {
      fHoughTransformer[i]->CreateHistograms(64,-0.006,0.006,64,-0.26,0.26);
      fHoughTransformer[i]->SetThreshold(3);
      fMemHandler[i]->Free();
      UInt_t ndigits=0;
      AliL3DigitRowData *digits =0;
      Char_t name[256];
      if(fBinary)
	{
	  sprintf(name,"%sdigits_%d_%d.raw",fPath,slice,i);
	  fMemHandler[i]->SetBinaryInput(name);
	  digits = (AliL3DigitRowData *)fMemHandler[i]->CompBinary2Memory(ndigits);
	  fMemHandler[i]->CloseBinaryInput();
	}
      else //read data from root file
	{
	  fMemHandler[i]->SetAliInput(fRootFile);
	  fMemHandler[i]->Init(slice,i,NRows[i]);
	  digits=(AliL3DigitRowData *)fMemHandler[i]->AliDigits2Memory(ndigits); 
	}
      fHoughTransformer[i]->SetInputData(ndigits,digits);
      fHoughTransformer[i]->TransformCircle();
    }
  
}

AliL3Histogram *AliL3Hough::AddHistograms()
{
  AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(0);
  for(Int_t i=1; i<NPatches; i++)
    {
      AliL3Histogram *hist = fHoughTransformer[i]->GetHistogram(0);
      hist0->Add(hist);
    }
  
  return hist0;
}

void AliL3Hough::Evaluate(AliL3Histogram *hist)
{
  
  AliL3HoughEval **eval = new AliL3HoughEval*[NPatches];
  for(Int_t i=0; i<NPatches; i++)
    {
      eval[i] = new AliL3HoughEval(fHoughTransformer[i]);
      eval[i]->DisplayEtaSlice(0,hist);
      delete eval[i];
    }
  
  delete [] eval;

}
