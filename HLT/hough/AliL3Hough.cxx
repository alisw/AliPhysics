// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"
#include <sys/time.h>

#include "AliL3Logging.h"
#include "AliL3HoughMerger.h"
#include "AliL3HoughIntMerger.h"
#include "AliL3HoughGlobalMerger.h"
#include "AliL3Histogram.h"
#include "AliL3Hough.h"
#include "AliL3HoughTransformer.h"
#include "AliL3HoughClusterTransformer.h"
#include "AliL3HoughTransformerLUT.h"
#include "AliL3HoughTransformerVhdl.h"
#include "AliL3HoughTransformerRow.h"
#include "AliL3HoughMaxFinder.h"
#include "AliL3Benchmark.h"
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
#include "AliL3DDLDataFileHandler.h"

#include "TThread.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** /class AliL3Hough
//<pre>
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
//
//</pre>
*/

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
  fGlobalTracks     = 0;
  fMerger           = 0;
  fInterMerger      = 0;
  fGlobalMerger     = 0;
  fBenchmark        = 0;
  
  fNEtaSegments     = 0;
  fNPatches         = 0;
  fLastPatch        =-1;
  fVersion          = 0;
  fCurrentSlice     = 0;
  fEvent            = 0;
  
  fKappaSpread = 6;
  fPeakRatio   = 0.5;
  fInputFile   = 0;
  fInputPtr    = 0;
  fRawEvent    = 0;
  
  SetTransformerParams();
  SetThreshold();
  SetNSaveIterations();
  SetPeakThreshold();
#ifdef use_aliroot
  //just be sure that index is empty for new event
    AliL3FileHandler::CleanStaticIndex(); 
#ifdef use_newio
    fRunLoader = 0;
#endif
#endif
  fThread = 0;
}

AliL3Hough::AliL3Hough(Char_t *path,Bool_t binary,Int_t netasegments,Bool_t bit8,Int_t tv,Char_t *infile,Char_t *ptr)
{
  //Normal constructor
  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments  = netasegments;
  fAddHistograms = kFALSE;
  fDoIterative   = kFALSE; 
  fWriteDigits   = kFALSE;
  fUse8bits      = bit8;
  fVersion       = tv;
  fKappaSpread=6;
  fPeakRatio=0.5;
  if(!fBinary) {
    if(infile) {
      fInputFile = infile;
      fInputPtr = 0;
    }
    else {
      fInputFile = 0;
      fInputPtr = ptr;
    }
  }
  else {
    fInputFile = 0;
    fInputPtr = 0;
  }
#ifdef use_aliroot
  //just be sure that index is empty for new event
    AliL3FileHandler::CleanStaticIndex(); 
#ifdef use_newio
    fRunLoader = 0;
#endif
#endif
  fThread = 0;
}

AliL3Hough::~AliL3Hough()
{
  //dtor

  CleanUp();
  if(fMerger)
    delete fMerger;
  //cout << "Cleaned class merger " << endl;
  if(fInterMerger)
    delete fInterMerger;
  //cout << "Cleaned class inter " << endl;
  if(fPeakFinder)
    delete fPeakFinder;
  //cout << "Cleaned class peak " << endl;
  if(fGlobalMerger)
    delete fGlobalMerger;
  //cout << "Cleaned class global " << endl;
  if(fBenchmark)
    delete fBenchmark;
  //cout << "Cleaned class bench " << endl;
  if(fGlobalTracks)
    delete fGlobalTracks;
  //cout << "Cleaned class globaltracks " << endl;
  if(fThread) {
    //    fThread->Delete();
    delete fThread;
    fThread = 0;
  }
}

void AliL3Hough::CleanUp()
{
  //Cleanup memory
  
  for(Int_t i=0; i<fNPatches; i++)
    {
      if(fTracks[i]) delete fTracks[i];
      //cout << "Cleaned tracks " << i << endl;
      if(fEval[i]) delete fEval[i];
      //cout << "Cleaned eval " << i << endl;
      if(fHoughTransformer[i]) delete fHoughTransformer[i];
      //cout << "Cleaned traf " << i << endl;
      if(fMemHandler[i]) delete fMemHandler[i];
      //cout << "Cleaned mem " << i << endl;
    }
  
  if(fTracks) delete [] fTracks;
  //cout << "Cleaned class tracks " << endl;
  if(fEval) delete [] fEval;
  //cout << "Cleaned class eval " << endl;
  if(fHoughTransformer) delete [] fHoughTransformer;
  //cout << "Cleaned cleass trafo " << endl;
  if(fMemHandler) delete [] fMemHandler;
  //cout << "Cleaned class mem " << endl;
}

void AliL3Hough::Init(Int_t netasegments,Int_t tv,AliRawEvent *rawevent,Float_t zvertex)
{
  //Normal constructor
  fNEtaSegments  = netasegments;
  fVersion       = tv;
  fRawEvent      = rawevent;
  fZVertex       = zvertex;

  Init();
}

void AliL3Hough::Init(Char_t *path,Bool_t binary,Int_t netasegments,Bool_t bit8,Int_t tv,Char_t *infile,Char_t *ptr,Float_t zvertex)
{
  //Normal init of the AliL3Hough
  fBinary = binary;
  strcpy(fPath,path);
  fNEtaSegments = netasegments;
  fWriteDigits  = kFALSE;
  fUse8bits     = bit8;
  fVersion      = tv;
  if(!fBinary) {
    if(infile) {
      fInputFile = infile;
      fInputPtr = 0;
    }
    else {
      fInputFile = 0;
      fInputPtr = ptr;
    }
  }
  else {
    fInputFile = 0;
    fInputPtr = 0;
  }
  fZVertex = zvertex;

  Init(); //do the rest
}

void AliL3Hough::Init(Bool_t doit, Bool_t addhists)
{
  // Init
  fDoIterative   = doit; 
  fAddHistograms = addhists;

  fNPatches = AliL3Transform::GetNPatches();
  fHoughTransformer = new AliL3HoughBaseTransformer*[fNPatches];
  fMemHandler = new AliL3MemHandler*[fNPatches];

  fTracks = new AliL3TrackArray*[fNPatches];
  fEval = new AliL3HoughEval*[fNPatches];
  
  fGlobalTracks = new AliL3TrackArray("AliL3HoughTrack");
  
  AliL3HoughBaseTransformer *lasttransformer = 0;

  for(Int_t i=0; i<fNPatches; i++)
    {
      switch (fVersion){ //choose Transformer
      case 1: 
	fHoughTransformer[i] = new AliL3HoughTransformerLUT(0,i,fNEtaSegments);
	break;
      case 2:
	fHoughTransformer[i] = new AliL3HoughClusterTransformer(0,i,fNEtaSegments);
	break;
      case 3:
	fHoughTransformer[i] = new AliL3HoughTransformerVhdl(0,i,fNEtaSegments,fNSaveIterations);
	break;
      case 4:
	fHoughTransformer[i] = new AliL3HoughTransformerRow(0,i,fNEtaSegments,kFALSE,fZVertex);
	break;
      default:
	fHoughTransformer[i] = new AliL3HoughTransformer(0,i,fNEtaSegments,kFALSE,kFALSE);
      }

      fHoughTransformer[i]->SetLastTransformer(lasttransformer);
      lasttransformer = fHoughTransformer[i];
      //      fHoughTransformer[i]->CreateHistograms(fNBinX[i],fLowPt[i],fNBinY[i],-fPhi[i],fPhi[i]);
      fHoughTransformer[i]->CreateHistograms(fNBinX[i],-fLowPt[i],fLowPt[i],fNBinY[i],-fPhi[i],fPhi[i]);
      //fHoughTransformer[i]->CreateHistograms(fLowPt[i],fUpperPt[i],fPtRes[i],fNBinY[i],fPhi[i]);

      fHoughTransformer[i]->SetLowerThreshold(fThreshold[i]);
      fHoughTransformer[i]->SetUpperThreshold(100);

      LOG(AliL3Log::kInformational,"AliL3Hough::Init","Version")
	<<"Initializing Hough transformer version "<<fVersion<<ENDLOG;
      
      fEval[i] = new AliL3HoughEval();
      fTracks[i] = new AliL3TrackArray("AliL3HoughTrack");
      if(fUse8bits)
	fMemHandler[i] = new AliL3DataHandler();
      else
#ifdef use_aliroot
      	{
	  if(!fRawEvent) {
	    if(!fInputFile) {
	      if(!fInputPtr) {
		/* In case of reading digits file */
		fMemHandler[i] = new AliL3FileHandler(kTRUE); //use static index
		if(!fBinary) {
#if use_newio
		  if(!fRunLoader) {
#endif
		    Char_t filename[1024];
		    sprintf(filename,"%s/digitfile.root",fPath);
		    fMemHandler[i]->SetAliInput(filename);
#if use_newio
		  }
		  else {
		    fMemHandler[i]->SetAliInput(fRunLoader);
		  }
#endif
		}
	      }
	      else {
		/* In case of reading from DATE */
		fMemHandler[i] = new AliL3DDLDataFileHandler();
		fMemHandler[i]->SetReaderInput(fInputPtr,-1);
	      }
	    }
	    else {
	      /* In case of reading rawdata from ROOT file */
	      fMemHandler[i] = new AliL3DDLDataFileHandler();
	      fMemHandler[i]->SetReaderInput(fInputFile);
	    }
	  }
	  else {
	    /* In case of reading rawdata using AliRawEvent */
	    fMemHandler[i] = new AliL3DDLDataFileHandler();
	    fMemHandler[i]->SetReaderInput(fRawEvent);
	  }
	}
#else
      fMemHandler[i] = new AliL3MemHandler();
#endif
    }

  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi",50000);
  fMerger = new AliL3HoughMerger(fNPatches);
  fInterMerger = new AliL3HoughIntMerger();
  fGlobalMerger = 0;
  fBenchmark = new AliL3Benchmark();
}

void AliL3Hough::SetTransformerParams(Float_t ptres,Float_t ptmin,Float_t ptmax,Int_t ny,Int_t patch)
{
  // Setup the parameters for the Hough Transformer
  Int_t mrow;
  Float_t psi=0;
  if(patch==-1)
    mrow = 80;
  else
    mrow = AliL3Transform::GetLastRow(patch);
  if(ptmin)
    {
      Double_t lineradius = sqrt(pow(AliL3Transform::Row2X(mrow),2) + pow(AliL3Transform::GetMaxY(mrow),2));
      Double_t kappa = -1*AliL3Transform::GetBField()*AliL3Transform::GetBFact()/ptmin;
      psi = AliL3Transform::Deg2Rad(10) - asin(lineradius*kappa/2);
      cout<<"Calculated psi range "<<psi<<" in patch "<<patch<<endl;
    }

  if(patch==-1)
    {
      Int_t i=0;
      while(i < 6)
	{
	  fPtRes[i] = ptres;
	  fLowPt[i] = ptmin;
	  fUpperPt[i] = ptmax;
	  fNBinY[i] = ny;
	  fPhi[i] = psi;
	  fNBinX[i]=0;
	  i++;
	}
      return;
    }

  fPtRes[patch] = ptres;
  fLowPt[patch] = ptmin;
  fUpperPt[patch] = ptmax;
  fNBinY[patch] = ny;
  fPhi[patch] = psi;
}
/*
void AliL3Hough::SetTransformerParams(Int_t nx,Int_t ny,Float_t ptmin,Int_t patch)
{
  // Setup the parameters for the Hough Transformer

  Int_t mrow=80;
  Double_t lineradius = sqrt(pow(AliL3Transform::Row2X(mrow),2) + pow(AliL3Transform::GetMaxY(mrow),2));
  Double_t kappa = -1*AliL3Transform::GetBField()*AliL3Transform::GetBFact()/ptmin;
  Double_t psi = AliL3Transform::Deg2Rad(10) - asin(lineradius*kappa/2);
  cout<<"Calculated psi range "<<psi<<" in patch "<<patch<<endl;
  
  Int_t i=0;
  while(i < 6)
    {
      fLowPt[i] = ptmin;
      fNBinY[i] = ny;
      fNBinX[i] = nx;
      fPhi[i] = psi;
      i++;
    }
}
*/
void AliL3Hough::SetTransformerParams(Int_t nx,Int_t ny,Float_t ptmin,Int_t /*patch*/)
{
  // Setup the parameters for the Hough Transformer


  Int_t mrow=79;
  Double_t lineradius = sqrt(pow(AliL3Transform::Row2X(mrow),2) + pow(AliL3Transform::GetMaxY(mrow),2));
  Double_t alpha1 = AliL3Transform::GetMaxY(mrow)/pow(lineradius,2);
  Double_t kappa = 1*AliL3Transform::GetBField()*AliL3Transform::GetBFact()/ptmin;
  Double_t psi = AliL3Transform::Deg2Rad(10) - asin(lineradius*kappa/2);
  //  cout<<"Calculated psi range "<<psi<<" in patch "<<patch<<endl;
  AliL3HoughTrack track;
  track.SetTrackParameters(kappa,psi,1);
  Float_t hit[3];
  Int_t mrow2 = 158;
  track.GetCrossingPoint(mrow2,hit);
  Double_t lineradius2 = sqrt(pow(AliL3Transform::Row2X(mrow2),2) + pow(AliL3Transform::GetMaxY(mrow2),2));
  Double_t alpha2 = hit[1]/pow(lineradius2,2);
  //  cout<<"Calculated alphas range "<<alpha1<<" "<<alpha2<<" in patch "<<patch<<endl;

  Int_t i=0;
  while(i < 6)
    {
      fLowPt[i] = 1.15*alpha1;
      fNBinY[i] = ny;
      fNBinX[i] = nx;
      fPhi[i] = 1.15*alpha2;
      i++;
    }
}

void AliL3Hough::SetTransformerParams(Int_t nx,Int_t ny,Float_t lpt,Float_t phi)
{
  Int_t i=0;
  while(i < 6)
    {
      fLowPt[i] = lpt;
      fNBinY[i] = ny;
      fNBinX[i] = nx;
      fPhi[i] = phi;
      i++;
    }
}

void AliL3Hough::SetThreshold(Int_t t3,Int_t patch)
{
  // Set digits threshold
  if(patch==-1)
    {
      Int_t i=0;
      while(i < 6)
	fThreshold[i++]=t3;
      return;
    }
  fThreshold[patch]=t3;
}

void AliL3Hough::SetPeakThreshold(Int_t threshold,Int_t patch)
{
  // Set Peak Finder threshold
  if(patch==-1)
    {
      Int_t i=0;
      while(i < 6)
	fPeakThreshold[i++]=threshold;
      return;
    }
  fPeakThreshold[patch]=threshold;
}

void AliL3Hough::DoBench(Char_t *name)
{
  fBenchmark->Analyze(name);
}

void AliL3Hough::Process(Int_t minslice,Int_t maxslice)
{
  //Process all slices [minslice,maxslice].
  fGlobalMerger = new AliL3HoughGlobalMerger(minslice,maxslice);
  
  for(Int_t i=minslice; i<=maxslice; i++)
    {
      ReadData(i);
      Transform();
      if(fAddHistograms) {
	if(fVersion != 4)
	  AddAllHistograms();
	else
	  AddAllHistogramsRows();
      }
      FindTrackCandidates();
      //Evaluate();
      //fGlobalMerger->FillTracks(fTracks[0],i);
    }
}

void AliL3Hough::ReadData(Int_t slice,Int_t eventnr)
{
  //Read data from files, binary or root.
  
#ifdef use_aliroot
  if(fEvent!=eventnr) //just be sure that index is empty for new event
    AliL3FileHandler::CleanStaticIndex(); 
#endif
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
	    sprintf(name,"%s/binaries/digits_c8_%d_%d_%d.raw",fPath,eventnr,slice,i);
	  else
	    sprintf(name,"%s/binaries/digits_%d_%d_%d.raw",fPath,eventnr,slice,i);

	  fMemHandler[i]->SetBinaryInput(name);
	  digits = (AliL3DigitRowData *)fMemHandler[i]->CompBinary2Memory(ndigits);
	  fMemHandler[i]->CloseBinaryInput();
	}
      else //read data from root file
	{
#ifdef use_aliroot
	  if(fEvent!=eventnr)
	    fMemHandler[i]->FreeDigitsTree();//or else the new event is not loaded
	  digits=(AliL3DigitRowData *)fMemHandler[i]->AliAltroDigits2Memory(ndigits,eventnr);
#else
	  cerr<<"You cannot read from rootfile now"<<endl;
#endif
	}

      //Set the pointer to the TPCRawStream in case of fast raw data reading
      fHoughTransformer[i]->SetTPCRawStream(fMemHandler[i]->GetTPCRawStream());

      //set input data and init transformer
      fHoughTransformer[i]->SetInputData(ndigits,digits);
      fHoughTransformer[i]->Init(slice,i,fNEtaSegments);
    }

  fEvent=eventnr;
}

void AliL3Hough::Transform(Int_t *rowrange)
{
  //Transform all data given to the transformer within the given slice
  //(after ReadData(slice))
  
  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  Int_t patchorder[6] = {5,2,0,1,3,4}; //The order in which patches are processed
  //  Int_t patchorder[6] = {0,1,2,3,4,5}; //The order in which patches are processed
  //  Int_t patchorder[6] = {5,4,3,2,1,0}; //The order in which patches are processed
  //  Int_t patchorder[6] = {5,2,4,3,1,0}; //The order in which patches are processed
  fLastPatch=-1;
  for(Int_t i=0; i<fNPatches; i++)
    {
      // In case of Row transformer reset the arrays only once
      if((fVersion != 4) || (i == 0)) {
	fBenchmark->Start("Hough Reset");
	fHoughTransformer[0]->Reset();//Reset the histograms
	fBenchmark->Stop("Hough Reset");
      }
      fBenchmark->Start("Hough Transform");
      PrepareForNextPatch(patchorder[i]);
      if(!rowrange) {
	char buf[256];
	sprintf(buf,"Patch %d",patchorder[i]);
	fBenchmark->Start(buf);
	fHoughTransformer[patchorder[i]]->SetLastPatch(fLastPatch);
	fHoughTransformer[patchorder[i]]->TransformCircle();
	fBenchmark->Stop(buf);
      }
      else
	fHoughTransformer[i]->TransformCircleC(rowrange,1);
      fBenchmark->Stop("Hough Transform");
      fLastPatch=patchorder[i];
    }
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::Transform()","Timing")
    <<"Transform done in average per patch of "<<cpuTime*1000/fNPatches<<" ms"<<ENDLOG;
}

void AliL3Hough::MergePatches()
{
  // Merge patches if they are not summed
  if(fAddHistograms) //Nothing to merge here
    return;
  fMerger->MergePatches(kTRUE);
}

void AliL3Hough::MergeInternally()
{
  // Merge patches internally
  if(fAddHistograms)
    fInterMerger->FillTracks(fTracks[0]);
  else
    fInterMerger->FillTracks(fMerger->GetOutTracks());
  
  fInterMerger->MMerge();
}

void AliL3Hough::ProcessSliceIter()
{
  //Process current slice (after ReadData(slice)) iteratively.
  
  if(!fAddHistograms)
    {
      for(Int_t i=0; i<fNPatches; i++)
	{
	  ProcessPatchIter(i);
	  fMerger->FillTracks(fTracks[i],i); //Copy tracks to merger
	}
    }
  else
    {
      for(Int_t i=0; i<10; i++)
	{
	  Transform();
	  AddAllHistograms();
	  InitEvaluate();
	  AliL3HoughBaseTransformer *tr = fHoughTransformer[0];
	  for(Int_t j=0; j<fNEtaSegments; j++)
	    {
	      AliL3Histogram *hist = tr->GetHistogram(j);
	      if(hist->GetNEntries()==0) continue;
	      fPeakFinder->Reset();
	      fPeakFinder->SetHistogram(hist);
	      fPeakFinder->FindAbsMaxima();
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[0]->NextTrack();
	      track->SetTrackParameters(fPeakFinder->GetXPeak(0),fPeakFinder->GetYPeak(0),fPeakFinder->GetWeight(0));
	      track->SetEtaIndex(j);
	      track->SetEta(tr->GetEta(j,fCurrentSlice));
	      for(Int_t k=0; k<fNPatches; k++)
		{
		  fEval[i]->SetNumOfPadsToLook(2);
		  fEval[i]->SetNumOfRowsToMiss(2);
		  fEval[i]->RemoveFoundTracks();
		  /*
		  Int_t nrows=0;
		  if(!fEval[i]->LookInsideRoad(track,nrows))
		    {
		      fTracks[0]->Remove(fTracks[0]->GetNTracks()-1);
		      fTracks[0]->Compress();
		    }
		  */
		}
	    }
	  
	}
      
    }
}

void AliL3Hough::ProcessPatchIter(Int_t patch)
{
  //Process patch in a iterative way. 
  //transform + peakfinding + evaluation + transform +...

  Int_t numoftries = 5;
  AliL3HoughBaseTransformer *tr = fHoughTransformer[patch];
  AliL3TrackArray *tracks = fTracks[patch];
  tracks->Reset();
  AliL3HoughEval *ev = fEval[patch];
  ev->InitTransformer(tr);
  //ev->RemoveFoundTracks();
  ev->SetNumOfRowsToMiss(3);
  ev->SetNumOfPadsToLook(2);
  AliL3Histogram *hist;
  for(Int_t t=0; t<numoftries; t++)
    {
      tr->Reset();
      tr->TransformCircle();
      for(Int_t i=0; i<fNEtaSegments; i++)
	{
	  hist = tr->GetHistogram(i);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->Reset();
	  fPeakFinder->SetHistogram(hist);
	  fPeakFinder->FindAbsMaxima();
	  //fPeakFinder->FindPeak1();
	  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();
	  track->SetTrackParameters(fPeakFinder->GetXPeak(0),fPeakFinder->GetYPeak(0),fPeakFinder->GetWeight(0));
	  track->SetEtaIndex(i);
	  track->SetEta(tr->GetEta(i,fCurrentSlice));
	  /*
	  Int_t nrows=0;
	  if(!ev->LookInsideRoad(track,nrows))
	    {	
	      tracks->Remove(tracks->GetNTracks()-1);
	      tracks->Compress();
	    }
	  */
	}
    }
  fTracks[0]->QSort();
  LOG(AliL3Log::kInformational,"AliL3Hough::ProcessPatch","NTracks")
    <<AliL3Log::kDec<<"Found "<<tracks->GetNTracks()<<" tracks in patch "<<patch<<ENDLOG;
}

void AliL3Hough::AddAllHistograms()
{
  //Add the histograms within one etaslice.
  //Resulting histogram are in patch=0.

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Add Histograms");
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      AliL3Histogram *hist0 = fHoughTransformer[0]->GetHistogram(i);
      for(Int_t j=1; j<fNPatches; j++)
	{
	  AliL3Histogram *hist = fHoughTransformer[j]->GetHistogram(i);
	  hist0->Add(hist);
	}
    }
  fBenchmark->Stop("Add Histograms");
  fAddHistograms = kTRUE;
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::AddAllHistograms()","Timing")
    <<"Adding histograms in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliL3Hough::AddAllHistogramsRows()
{
  //Add the histograms within one etaslice.
  //Resulting histogram are in patch=0.

  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Add HistogramsRows");

  UChar_t lastpatchlastrow = AliL3Transform::GetLastRowOnDDL(fLastPatch)+1;

  UChar_t *tracklastrow = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetTrackLastRow();

  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      UChar_t *gapcount = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetGapCount(i);
      UChar_t *currentrowcount = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetCurrentRowCount(i);

      AliL3Histogram *hist = fHoughTransformer[0]->GetHistogram(i);
      Int_t xmin = hist->GetFirstXbin();
      Int_t xmax = hist->GetLastXbin();
      Int_t ymin = hist->GetFirstYbin();
      Int_t ymax = hist->GetLastYbin();
      Int_t nxbins = hist->GetNbinsX()+2;

      for(Int_t ybin=ymin; ybin<=ymax; ybin++)
	{
	  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
	    {
	      Int_t bin = xbin + ybin*nxbins; //Int_t bin = hist->GetBin(xbin,ybin);
	      if(gapcount[bin] < MAX_N_GAPS) {
		if(tracklastrow[bin] > lastpatchlastrow) {
		  if(lastpatchlastrow > currentrowcount[bin])
		    gapcount[bin] += (lastpatchlastrow-currentrowcount[bin]-1);
		}
		else {
		  if(tracklastrow[bin] > currentrowcount[bin])
		    gapcount[bin] += (tracklastrow[bin]-currentrowcount[bin]-1);
		}
		if(gapcount[bin] < MAX_N_GAPS)
		  hist->AddBinContent(bin,(159-gapcount[bin]));
	      }
	    }
	}
    }

  fBenchmark->Stop("Add HistogramsRows");
  fAddHistograms = kTRUE;
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::AddAllHistogramsRows()","Timing")
    <<"Adding histograms in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliL3Hough::PrepareForNextPatch(Int_t nextpatch)
{
  char buf[256];
  sprintf(buf,"Prepare For Patch %d",nextpatch);
  fBenchmark->Start(buf);

  UChar_t lastpatchlastrow;
  if(fLastPatch == -1)
    lastpatchlastrow = 0;
  else
    lastpatchlastrow = AliL3Transform::GetLastRowOnDDL(fLastPatch)+1;
  UChar_t nextpatchfirstrow;
  if(nextpatch==0)
    nextpatchfirstrow = 0;
  else
    nextpatchfirstrow = AliL3Transform::GetFirstRowOnDDL(nextpatch)-1;

  UChar_t *trackfirstrow = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetTrackFirstRow();
  UChar_t *tracklastrow = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetTrackLastRow();

  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      UChar_t *gapcount = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetGapCount(i);
      UChar_t *currentrowcount = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetCurrentRowCount(i);
      UChar_t *prevbin = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetPrevBin(i);
      UChar_t *nextbin = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetNextBin(i);
      UChar_t *nextrow = ((AliL3HoughTransformerRow *)fHoughTransformer[0])->GetNextRow(i);

      AliL3Histogram *hist = fHoughTransformer[0]->GetHistogram(i);
      Int_t xmin = hist->GetFirstXbin();
      Int_t xmax = hist->GetLastXbin();
      Int_t ymin = hist->GetFirstYbin();
      Int_t ymax = hist->GetLastYbin();
      Int_t nxbins = hist->GetNbinsX()+2;

      if(fLastPatch != -1) {
	UChar_t lastyvalue = 0;
	Int_t endybin = ymin - 1;
	for(Int_t ybin=nextrow[ymin]; ybin<=ymax; ybin = nextrow[++ybin])
	  {
	    UChar_t lastxvalue = 0;
	    UChar_t maxvalue = 0;
	    Int_t endxbin = xmin - 1;
	    for(Int_t xbin=xmin; xbin<=xmax; xbin++)
	      {
		Int_t bin = xbin + ybin*nxbins;
		UChar_t value = 0;
		if(gapcount[bin] < MAX_N_GAPS) {
		  if(tracklastrow[bin] > lastpatchlastrow) {
		    if(lastpatchlastrow > currentrowcount[bin])
		      gapcount[bin] += (lastpatchlastrow-currentrowcount[bin]-1);
		  }
		  else {
		    if(tracklastrow[bin] > currentrowcount[bin])
		      gapcount[bin] += (tracklastrow[bin]-currentrowcount[bin]-1);
		  }
		  if(gapcount[bin] < MAX_N_GAPS) {
		    value = 1;
		    maxvalue = 1;
		    if(trackfirstrow[bin] < nextpatchfirstrow)
		      currentrowcount[bin] = nextpatchfirstrow;
		    else
		      currentrowcount[bin] = trackfirstrow[bin];
		  }
		}
		if(value > 0)
		  {
		    nextbin[xbin + ybin*nxbins] = (UChar_t)xbin;
		    prevbin[xbin + ybin*nxbins] = (UChar_t)xbin;
		    if(value > lastxvalue)
		      {
			UChar_t *tempnextbin = nextbin + endxbin + 1 + ybin*nxbins;
			memset(tempnextbin,(UChar_t)xbin,xbin-endxbin-1);
		      }
		    endxbin = xbin;
		  }
		else
		  {
		    prevbin[xbin + ybin*nxbins] = (UChar_t)endxbin;
		  }
		lastxvalue = value;
	      }
	    UChar_t *tempnextbin = nextbin + endxbin + 1 + ybin*nxbins;
	    memset(tempnextbin,(UChar_t)(xmax+1),xmax-endxbin);
	    if(maxvalue > 0)
	      {
		nextrow[ybin] = (UChar_t)ybin;
		if(maxvalue > lastyvalue)
		  {
		    UChar_t *tempnextrow = nextrow + endybin + 1;
		    memset(tempnextrow,(UChar_t)ybin,ybin-endybin-1);
		  }
		endybin = ybin;
	      }
	    lastyvalue = maxvalue;
	  }
	UChar_t *tempnextrow = nextrow + endybin + 1;
	memset(tempnextrow,(UChar_t)(ymax+1),ymax-endybin);
      }
      else {
	UChar_t lastyvalue = 0;
	Int_t endybin = ymin - 1;
	for(Int_t ybin=ymin; ybin<=ymax; ybin++)
	  {
	    UChar_t maxvalue = 0;
	    for(Int_t xbin=xmin; xbin<=xmax; xbin++)
	      {
		Int_t bin = xbin + ybin*nxbins;
		if(gapcount[bin] < MAX_N_GAPS) {
		  maxvalue = 1;
		  if(trackfirstrow[bin] < nextpatchfirstrow)
		    currentrowcount[bin] = nextpatchfirstrow;
		  else
		    currentrowcount[bin] = trackfirstrow[bin];
		}
	      }
	    if(maxvalue > 0)
	      {
		nextrow[ybin] = (UChar_t)ybin;
		if(maxvalue > lastyvalue)
		  {
		    UChar_t *tempnextrow = nextrow + endybin + 1;
		    memset(tempnextrow,(UChar_t)ybin,ybin-endybin-1);
		  }
		endybin = ybin;
	      }
	    lastyvalue = maxvalue;
	  }
	UChar_t *tempnextrow = nextrow + endybin + 1;
	memset(tempnextrow,(UChar_t)(ymax+1),ymax-endybin);
      }
    }

  fBenchmark->Stop(buf);
}

void AliL3Hough::AddTracks()
{
  // Add current slice slice tracks to the global list of found tracks
  if(!fTracks[0])
    {
      cerr<<"AliL3Hough::AddTracks : No tracks"<<endl;
      return;
    }
  AliL3TrackArray *tracks = fTracks[0];
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3Track *track = tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(track->GetNHits()!=1) cerr<<"NHITS "<<track->GetNHits()<<endl;
      UInt_t *ids = track->GetHitNumbers();
      ids[0] = (fCurrentSlice&0x7f)<<25;
    }
  
  fGlobalTracks->AddTracks(fTracks[0],0,fCurrentSlice);
}

void AliL3Hough::FindTrackCandidatesRow()
{
  // Find AliL3HoughTransformerRow track candidates
  if(fVersion != 4) {
    LOG(AliL3Log::kError,"AliL3Hough::FindTrackCandidatesRow()","")
      <<"Incompatible Peak Finder version!"<<ENDLOG;
    return;
  }

  //Look for peaks in histograms, and find the track candidates
  Int_t npatches;
  if(fAddHistograms)
    npatches = 1; //Histograms have been added.
  else
    npatches = fNPatches;
  
  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Find Maxima");
  for(Int_t i=0; i<npatches; i++)
    {
      AliL3HoughBaseTransformer *tr = fHoughTransformer[i];
      fTracks[i]->Reset();
      fPeakFinder->Reset();
      
      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->SetHistogram(hist);
	  fPeakFinder->SetEtaSlice(j);
	  fPeakFinder->SetTrackLUTs(((AliL3HoughTransformerRow *)tr)->GetTrackNRows(),((AliL3HoughTransformerRow *)tr)->GetTrackFirstRow(),((AliL3HoughTransformerRow *)tr)->GetTrackLastRow());
#ifdef do_mc
	  LOG(AliL3Log::kInformational,"AliL3Hough::FindTrackCandidates()","")
	    <<"Starting "<<j<<" etaslice"<<ENDLOG;
#endif
	  fPeakFinder->SetThreshold(fPeakThreshold[i]);
	  fPeakFinder->FindAdaptedRowPeaks(1,0,0);//Maxima finder for HoughTransformerRow

	  //fPeakFinder->FindMaxima(fPeakThreshold[i]); //Simple maxima finder
	}
  
      for(Int_t k=0; k<fPeakFinder->GetEntries(); k++)
	{
	  //	  if(fPeakFinder->GetWeight(k) < 0) continue;
	  AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[i]->NextTrack();
	  Float_t psi = atan((fPeakFinder->GetXPeak(k)-fPeakFinder->GetYPeak(k))/(AliL3HoughTransformerRow::GetBeta1()-AliL3HoughTransformerRow::GetBeta2()));
	  Float_t kappa = 2.0*(fPeakFinder->GetXPeak(k)*cos(psi)-AliL3HoughTransformerRow::GetBeta1()*sin(psi));
	  //	      track->SetTrackParameters(fPeakFinder->GetXPeak(k),fPeakFinder->GetYPeak(k),fPeakFinder->GetWeight(k));
	  track->SetTrackParameters(kappa,psi,fPeakFinder->GetWeight(k));
	  track->SetBinXY(fPeakFinder->GetXPeak(k),fPeakFinder->GetYPeak(k),fPeakFinder->GetXPeakSize(k),fPeakFinder->GetYPeakSize(k));
	  Int_t etaindex = (fPeakFinder->GetStartEta(k)+fPeakFinder->GetEndEta(k))/2;
	  track->SetEtaIndex(etaindex);
	  Float_t starteta = tr->GetEta(fPeakFinder->GetStartEta(k),fCurrentSlice);
	  Float_t endeta = tr->GetEta(fPeakFinder->GetEndEta(k),fCurrentSlice);
	  track->SetEta((starteta+endeta)/2.0);
	  track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));
	  track->SetSector(fCurrentSlice);
	  track->SetSlice(fCurrentSlice);
#ifdef do_mc
	  Int_t label = tr->GetTrackID(etaindex,fPeakFinder->GetXPeak(k),fPeakFinder->GetYPeak(k));
	  track->SetMCid(label);
	  //	  cout<<"Track found with label "<<label<<" at "<<fPeakFinder->GetXPeak(k)<<" "<<fPeakFinder->GetYPeak(k)<<" with weight "<<fPeakFinder->GetWeight(k)<<endl; 
#endif
	}
      LOG(AliL3Log::kInformational,"AliL3Hough::FindTrackCandidates()","")
	<<"Found "<<fTracks[i]->GetNTracks()<<" tracks in slice "<<fCurrentSlice<<ENDLOG;
      fTracks[i]->QSort();
    }
  fBenchmark->Stop("Find Maxima");
  cpuTime = GetCpuTime() - initTime;
  LOG(AliL3Log::kInformational,"AliL3Hough::FindTrackCandidates()","Timing")
    <<"Maxima finding done in "<<cpuTime*1000<<" ms"<<ENDLOG;
}

void AliL3Hough::FindTrackCandidates()
{
  // Find AliL3HoughTransformer track candidates
  if(fVersion == 4) {
    LOG(AliL3Log::kError,"AliL3Hough::FindTrackCandidatesRow()","")
      <<"Incompatible Peak Finder version!"<<ENDLOG;
    return;
  }

  Int_t npatches;
  if(fAddHistograms)
    npatches = 1; //Histograms have been added.
  else
    npatches = fNPatches;
  
  Double_t initTime,cpuTime;
  initTime = GetCpuTime();
  fBenchmark->Start("Find Maxima");
  for(Int_t i=0; i<npatches; i++)
    {
      AliL3HoughBaseTransformer *tr = fHoughTransformer[i];
      fTracks[i]->Reset();
      
      for(Int_t j=0; j<fNEtaSegments; j++)
	{
	  AliL3Histogram *hist = tr->GetHistogram(j);
	  if(hist->GetNEntries()==0) continue;
	  fPeakFinder->Reset();
	  fPeakFinder->SetHistogram(hist);
#ifdef do_mc
	  cout<<"Starting "<<j<<" etaslice"<<endl;
#endif
	  fPeakFinder->SetThreshold(fPeakThreshold[i]);
	  fPeakFinder->FindAdaptedPeaks(fKappaSpread,fPeakRatio);
	  
	  for(Int_t k=0; k<fPeakFinder->GetEntries(); k++)
	    {
	      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks[i]->NextTrack();
	      track->SetTrackParameters(fPeakFinder->GetXPeak(k),fPeakFinder->GetYPeak(k),fPeakFinder->GetWeight(k));
	      track->SetEtaIndex(j);
	      track->SetEta(tr->GetEta(j,fCurrentSlice));
	      track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));
	    }
	}
      cout<<"Found "<<fTracks[i]->GetNTracks()<<" tracks in patch "<<i<<endl;
      fTracks[i]->QSort();
    }
  fBenchmark->Stop("Find Maxima");
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

Int_t AliL3Hough::Evaluate(Int_t roadwidth,Int_t nrowstomiss)
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
  
  Int_t removedtracks=0;
  AliL3TrackArray *tracks=0;

  if(fAddHistograms)
    {
      tracks = fTracks[0];
      for(Int_t i=0; i<tracks->GetNTracks(); i++)
	{
	  AliL3Track *track = tracks->GetCheckedTrack(i);
	  if(!track) continue;
	  track->SetNHits(0);
	}
    }
  
  for(Int_t i=0; i<fNPatches; i++)
    EvaluatePatch(i,roadwidth,nrowstomiss);
  
  //Here we check the tracks globally; 
  //how many good rows (padrows with signal) 
  //did it cross in the slice
  if(fAddHistograms) 
    {
      for(Int_t j=0; j<tracks->GetNTracks(); j++)
	{
	  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  
	  if(track->GetNHits() < AliL3Transform::GetNRows() - nrowstomiss)
	    {
	      tracks->Remove(j);
	      removedtracks++;
	    }
	}
      tracks->Compress();
      tracks->QSort();
    }
    
  return removedtracks;
}

void AliL3Hough::EvaluatePatch(Int_t i,Int_t roadwidth,Int_t nrowstomiss)
{
  //Evaluate patch i.
  
  fEval[i]->InitTransformer(fHoughTransformer[i]);
  fEval[i]->SetNumOfPadsToLook(roadwidth);
  fEval[i]->SetNumOfRowsToMiss(nrowstomiss);
  //fEval[i]->RemoveFoundTracks();
  
  AliL3TrackArray *tracks=0;
  
  if(!fAddHistograms)
    tracks = fTracks[i];
  else
    tracks = fTracks[0];
  
  Int_t nrows=0;
  for(Int_t j=0; j<tracks->GetNTracks(); j++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
      if(!track)
	{
	  LOG(AliL3Log::kWarning,"AliL3Hough::EvaluatePatch","Track array")
	    <<"Track object missing!"<<ENDLOG;
	  continue;
	} 
      nrows=0;
      Int_t rowrange[2] = {AliL3Transform::GetFirstRow(i),AliL3Transform::GetLastRow(i)};
      Bool_t result = fEval[i]->LookInsideRoad(track,nrows,rowrange);
      if(fAddHistograms)
	{
	  Int_t pre=track->GetNHits();
	  track->SetNHits(pre+nrows);
	}
      else//the track crossed too few good padrows (padrows with signal) in the patch, so remove it
	{
	  if(result == kFALSE)
	    tracks->Remove(j);
	}
    }
  
  tracks->Compress();

}

void AliL3Hough::MergeEtaSlices()
{
  //Merge tracks found in neighbouring eta slices.
  //Removes the track with the lower weight.
  
  fBenchmark->Start("Merge Eta-slices");
  AliL3TrackArray *tracks = fTracks[0];
  if(!tracks)
    {
      cerr<<"AliL3Hough::MergeEtaSlices : No tracks "<<endl;
      return;
    }
  for(Int_t j=0; j<tracks->GetNTracks(); j++)
    {
      AliL3HoughTrack *track1 = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
      if(!track1) continue;
      for(Int_t k=j+1; k<tracks->GetNTracks(); k++)
	{
	  AliL3HoughTrack *track2 = (AliL3HoughTrack*)tracks->GetCheckedTrack(k);
	  if(!track2) continue;
	  if(abs(track1->GetEtaIndex() - track2->GetEtaIndex()) != 1) continue;
	  if(fabs(track1->GetKappa()-track2->GetKappa()) < 0.006 && 
	     fabs(track1->GetPsi()- track2->GetPsi()) < 0.1)
	    {
	      //cout<<"Merging track in slices "<<track1->GetEtaIndex()<<" "<<track2->GetEtaIndex()<<endl;
	      if(track1->GetWeight() > track2->GetWeight())
		tracks->Remove(k);
	      else
		tracks->Remove(j);
	    }
	}
    }
  fBenchmark->Stop("Merge Eta-slices");
  tracks->Compress();
}

void AliL3Hough::WriteTracks(Char_t *path)
{
  // Write found tracks into file
  //cout<<"AliL3Hough::WriteTracks : Sorting the tracsk"<<endl;
  //fGlobalTracks->QSort();
  
  Char_t filename[1024];
  sprintf(filename,"%s/tracks_%d.raw",path,fEvent);
  AliL3MemHandler mem;
  mem.SetBinaryOutput(filename);
  mem.TrackArray2Binary(fGlobalTracks);
  mem.CloseBinaryOutput();
  fGlobalTracks->Reset();
}

void AliL3Hough::WriteTracks(Int_t slice,Char_t *path)
{
  // Write found tracks slice by slice into file
  
  AliL3MemHandler mem;
  Char_t fname[100];
  if(fAddHistograms)
    {
      sprintf(fname,"%s/tracks_ho_%d_%d.raw",path,fEvent,slice);
      mem.SetBinaryOutput(fname);
      mem.TrackArray2Binary(fTracks[0]);
      mem.CloseBinaryOutput();
    }
  else 
    {
      for(Int_t i=0; i<fNPatches; i++)
	{
	  sprintf(fname,"%s/tracks_ho_%d_%d_%d.raw",path,fEvent,slice,i);
	  mem.SetBinaryOutput(fname);
	  mem.TrackArray2Binary(fTracks[i]);
	  mem.CloseBinaryOutput();
	}
    }
}

#ifdef use_aliroot
Int_t AliL3Hough::FillESD(AliESD *esd)
{
  if(!fGlobalTracks) return 0;
  Int_t nglobaltracks = 0;
  for(Int_t i=0; i<fGlobalTracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *tpt = (AliL3HoughTrack *)fGlobalTracks->GetCheckedTrack(i);
      if(!tpt) continue; 
      
      AliESDHLTtrack *esdtrack = new AliESDHLTtrack(); 

      esdtrack->SetRowRange(tpt->GetFirstRow(),tpt->GetLastRow());
      esdtrack->SetNHits(tpt->GetNHits());
      esdtrack->SetFirstPoint(tpt->GetFirstPointX(),tpt->GetFirstPointY(),tpt->GetFirstPointZ());
      esdtrack->SetLastPoint(tpt->GetLastPointX(),tpt->GetLastPointY(),tpt->GetLastPointZ());
      esdtrack->SetPt(tpt->GetPt());
      esdtrack->SetPsi(tpt->GetPsi());
      esdtrack->SetTgl(tpt->GetTgl());
      esdtrack->SetCharge(tpt->GetCharge());
      esdtrack->SetMCid(tpt->GetMCid());
      esdtrack->SetWeight(tpt->GetWeight());
      esdtrack->SetSector(tpt->GetSector());
      esdtrack->SetBinXY(tpt->GetBinX(),tpt->GetBinY(),tpt->GetSizeX(),tpt->GetSizeY());
      esdtrack->SetPID(tpt->GetPID());
      esdtrack->ComesFromMainVertex(tpt->ComesFromMainVertex());

      esd->AddHLTHoughTrack(esdtrack);
      nglobaltracks++;
      delete esdtrack;
    }
  return nglobaltracks;
}
#endif

void AliL3Hough::WriteDigits(Char_t *outfile)
{
  //Write the current data to a new rootfile.
#ifdef use_aliroot  

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
}

void *AliL3Hough::ProcessInThread(void *args)
{
  AliL3Hough *instance = (AliL3Hough *)args;
  Int_t minslice = instance->GetMinSlice();
  Int_t maxslice = instance->GetMaxSlice();
  for(Int_t i=minslice; i<=maxslice; i++)
    {
      instance->ReadData(i,0);
      instance->Transform();
      instance->AddAllHistogramsRows();
      instance->FindTrackCandidatesRow();
      instance->AddTracks();
    }
  return (void *)0;
}

void AliL3Hough::StartProcessInThread(Int_t minslice,Int_t maxslice)
{
  if(!fThread) {
    char buf[255];
    sprintf(buf,"houghtrans_%d_%d",minslice,maxslice);
    SetMinMaxSlices(minslice,maxslice);
    //    fThread = new TThread(buf,(void (*) (void *))&ProcessInThread,(void *)this);
    fThread = new TThread(buf,&ProcessInThread,(void *)this);
    fThread->Run();
  }
  return;
}

Int_t AliL3Hough::WaitForThreadFinish()
{
#if ROOT_VERSION_CODE < 262403
  return TThread::Join(fThread->GetId());
#else
  return fThread->Join(fThread->GetId());
#endif
}
