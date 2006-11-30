//$Id$

/* Example of how to use the AliHLTClusterFitter to fit clusters
   to the track candidates given as a AliHLTTrackArray.

   The path "path" should contain the link to the digitsfile, 
   and the directories called fitter (for the results) and hough
   (for the track files). For NEWIO, make sure that the file 
   TPC.Digits.root is in the path (symlink)!

   Also provide the neccessary parameters in SetFitParameters.C.

   RUN with ALIROOT (not ROOT), no other method is 
   supported right now.
*/

#ifndef __CINT__
#include "AliHLTLogger.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTHough.h"
#include "AliHLTFitter.h"
#include "AliHLTClusterFitter.h"
#include "AliHLTVertex.h"
#include "AliHLTTrackArray.h"
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#endif

void deconvclusters(Char_t *path,Int_t minslice=0,Int_t maxslice=35,Int_t nevent=1)
{
  
  AliHLTTransform::Init(path,kTRUE);
  
  Char_t filename[1024];
  AliHLTFileHandler *file = new AliHLTFileHandler(kTRUE); //static index
  UInt_t ndigits=0;
  
  sprintf(filename,"%s/digitfile.root",path);
  file->SetAliInput(filename);
  AliHLTClusterFitter *fitter = new AliHLTClusterFitter(path);
  
#ifdef __CINT__
  Char_t macroname[1024];
  gROOT->LoadMacro("SetFitParameters.C");
  SetFitParameters(fitter);
#else /*compiled version */
  fitter->SetInnerWidthFactor(1,1);
  fitter->SetOuterWidthFactor(1,1);
  fitter->SetNmaxOverlaps(5);

  fitter->SetChiSqMax(5,0);
  fitter->SetChiSqMax(5,1);
  fitter->SetChiSqMax(5,2);
#endif
  
  TStopwatch tloader;tloader.Stop();
  TStopwatch tfinder;tfinder.Stop();
  TStopwatch trefitter;trefitter.Stop();

  Int_t patch=-1;
  Int_t rowrange[2] = {0,AliHLTTransform::GetNRows()-1};  
  for(Int_t ev=0; ev<nevent; ev++)
    {
      AliHLTFileHandler::LoadStaticIndex(0,ev);
      fitter->LoadSeeds(rowrange,kFALSE,ev); //Takes input from global hough tracks
      
      for(Int_t slice=minslice; slice<=maxslice; slice++)
	{
	  tloader.Start(0);
	  file->Init(slice,-1);
	  cout<<"Processing event "<<ev<<" slice "<<slice<<" patch "<<patch<<endl;
	  AliHLTDigitRowData *digits = (AliHLTDigitRowData*)file->AliAltroDigits2Memory(ndigits,ev);
	  fitter->Init(slice,patch);
	  fitter->SetInputData(digits);
	  tloader.Stop();
	  
	  tfinder.Start(0);
	  fitter->FindClusters();
	  fitter->WriteClusters();
	  tfinder.Stop();

	  tloader.Start(0);
	  file->Free();
	  tloader.Stop();
	}
      
      //If you want a refit of the clusters;-------------------------
      tloader.Start(0);
      AliHLTVertex vertex;
      AliHLTTrackArray *tracks = fitter->GetSeeds(); //The seeds are the 
                                                    //input tracks from circle HT
      AliHLTFitter *ft = new AliHLTFitter(&vertex,1);
      sprintf(filename,"%s/fitter/",path);
      ft->LoadClusters(filename,0,kTRUE);
      tloader.Stop();

      trefitter.Start(0);
      for(Int_t i=0; i<tracks->GetNTracks(); i++)
	{
	  AliHLTTrack *track = tracks->GetCheckedTrack(i);
	  if(!track) continue;
	  if(track->GetNHits() < 40) continue;
	  ft->SortTrackClusters(track);
	  ft->FitHelix(track);
	  track->UpdateToFirstPoint();
	}
      trefitter.Stop();
      //-------------------------------------------------------------
      tloader.Start(0);
      delete ft;
      fitter->WriteTracks(5); //Write the final tracks
      file->FreeDigitsTree();
      tloader.Stop();
    }
  cout << " --- Timing values --- " << endl;
  cout << "Data  Loading:        "; tloader.Print("m");
  cout << "Cluster Deconvolution "; tfinder.Print("m");
  cout << "Track ReFitter        "; trefitter.Print("m");

  delete fitter;
}

void deconvlocally(Char_t *path,Int_t minslice=0,Int_t maxslice=17)
{
  
  AliHLTTransform::Init(path,kTRUE);
  
  Char_t filename[1024];
  AliHLTFileHandler *file = new AliHLTFileHandler(kTRUE);
  UInt_t ndigits=0;
  
  sprintf(filename,"%s/digitfile.root",path);
  file->SetAliInput(filename);
  AliHLTClusterFitter *fitter = new AliHLTClusterFitter(path);
  
#ifdef __CINT__
  Char_t macroname[1024];
  sprintf(macroname,"%s/SetFitParameters.C",path);
  gROOT->LoadMacro(macroname);
  SetFitParameters(fitter);
#else /*compiled version */
  fitter->SetInnerWidthFactor(1,1);
  fitter->SetOuterWidthFactor(1,1);
  fitter->SetNmaxOverlaps(5);

  fitter->SetChiSqMax(5,0);
  fitter->SetChiSqMax(5,1);
  fitter->SetChiSqMax(5,2);
#endif
  
  Int_t patch=-1;
  
  for(Int_t slice=minslice; slice<=maxslice; slice++)
    {
      file->Init(slice,patch);
      cout<<"Processing slice "<<slice<<" patch "<<patch<<endl;
      AliHLTDigitRowData *digits = (AliHLTDigitRowData*)file->AliAltroDigits2Memory(ndigits);
      
      fitter->Init(slice,patch);
      fitter->LoadLocalSegments();
      
      fitter->SetInputData(digits);
      
      fitter->FindClusters();
      fitter->WriteClusters();
      
      file->Free();
    
    }
  
  fitter->WriteTracks(5);
  delete fitter;
}

