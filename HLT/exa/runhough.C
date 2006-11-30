//$Id$

/**
   Run this macro for Hough track candidate finder
   (see steering class AliHLTHough).
   In argument path, you have to provide the path to 
   the directory where the data files should be located. 
   In case of reading from a rootfile, you have to
   make a symbolic link "digitfile.root", which points 
   to the rootfile containing AliROOT digits tree 
   and a symbolic link "alirunfile.root" pointing to a file 
   containing the ALIROOT geometry (TPC param).
   For NEWIO, make sure that the file TPC.Digits.root is in the 
   path (or make a symlink to it)!
 
   Also provide the neccessary parameters in SetHoughParameters.C.

   RUN with ALIROOT (not ROOT) if using root files.

*/

#ifndef __CINT__
#include "AliHLTLogger.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include "AliHLTHough.h"
#include "AliHLTTrackArray.h"
#include "AliHLTTrack.h"
#include "AliHLTHoughTrack.h"
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#endif

void runhough(Char_t *path,Char_t *outpath,Int_t s1=0,Int_t s2=35,Int_t nevent=1)
{

  Bool_t isinit=AliHLTTransform::Init(path,kTRUE);
  if(!isinit){
    cerr << "Could not create transform settings, please check log for error messages!" << endl;
    return;
  }

  Int_t tversion=1; //0 = normal transformer
                    //1 = LUT transformer

  AliHLTHough *hough = new AliHLTHough();
#ifdef __CINT__
  Char_t macroname[1024];
  sprintf(macroname,"SetHoughParameters.C");
  gROOT->LoadMacro(macroname);
  SetHoughParameters(hough,path,tversion);
#else /*compiled version*/
  Bool_t binary = kFALSE;    //binary files input
  Int_t n_eta_segments=100;
  Double_t histptmin = 0.5; //mininum pt to find (controls the histogram range) 
  Int_t threshold=6000;  //peak threshold
  //Int_t threshold=5000;  //peak threshold
  Int_t nxbins = 140;
  Int_t nybins = 150;
  //Int_t nxbins = 190;
  //Int_t nybins = 200;
  Int_t patch=-1; //-1 -> Hough transform on slices (means adding histograms)
  hough->SetThreshold(4); //noise threshold on single digits
  hough->SetTransformerParams(nxbins,nybins,histptmin,patch);
  hough->SetPeakThreshold(threshold,patch);
  hough->Init(path,binary,n_eta_segments,kFALSE,tversion); 
#endif

  TStopwatch tloader;tloader.Stop();
  TStopwatch ttransform;ttransform.Stop();
  TStopwatch tfinder;tfinder.Stop();

  for(Int_t ev=0; ev<nevent; ev++)
    {
      AliHLTFileHandler::LoadStaticIndex(0,ev);
      for(Int_t slice=s1; slice<=s2; slice++)
	{
	  cout<<"Processing slice "<<slice<<endl;
	  tloader.Start(0);hough->ReadData(slice,ev);tloader.Stop();
	  ttransform.Start(0);hough->Transform();ttransform.Stop();
	  tfinder.Start(0);
	  hough->AddAllHistograms();
	  hough->FindTrackCandidates();
	  hough->AddTracks();
	  tfinder.Stop();
#if 0 /*print track list */
	  AliHLTTrackArray *tracks = (AliHLTTrackArray*)hough->GetTracks(0);
	  tracks->QSort();
	  for(int i=0; i<tracks->GetNTracks(); i++)
	    {
	      AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
	      if(!track) continue;
	      cout<<"pt "<<track->GetPt()<<" psi "<<track->GetPsi()<<" eta "<<track->GetEta()<<" etaindex "<<track->GetEtaIndex()<<" weight "<<track->GetWeight()<<endl;
	    }
#endif
	}
      hough->WriteTracks(outpath);
      AliHLTFileHandler::SaveStaticIndex(0,ev);
    }

  cout << " --- Timing values --- " << endl;
  cout << "Data  Loading:     "; tloader.Print("m");
  cout << "Hough Transforming "; ttransform.Print("m");
  cout << "Track Finding      "; tfinder.Print("m");
  
  delete hough;
}

