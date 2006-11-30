// $Id$

/**
   Run this macro for cluster finder and track follower 
   (see steering class AliLevel3).
   In argument path, you have to provide the path to the directory 
   where the data files should be located. In case of reading from a rootfile, you have to
   make a symbolic link "digitfile.root", which points to the rootfile containing AliROOT 
   digits tree and a symbolic link "alirunfile.root" pointing to a file containing
   the ALIROOT geometry (TPC param). For NEWIO, make sure that the 
   file TPC.Digits.root is in the path (symlink)!

   RUN with ALIROOT (not ROOT) if using root files.
*/

#ifndef __CINT__
#include "AliHLTLogger.h"
#include "AliHLTFileHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"
#include "AliLevel3.h"
#include <TNtuple.h>
#include <TRandom.h>
#include <TSystem.h>
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#endif

void runtracker(Int_t minslice=0,Int_t maxslice=35,Char_t* path="./",Int_t nevent=1,Char_t *opath="./tracker/")
{
  //Set your configuration here:
  AliLevel3::EFileType filetype=AliLevel3::kRoot; //Input is RLE binary files or rootfile.
  Bool_t pileup=kFALSE; //Assume input is pileup event = non RLE binary files.
  Bool_t nonvertex = kFALSE; //Set this to true if a second nonvertex contrained tracking pass should be performed.
  Int_t npatches = 1;   //Options; 1, 2 and 6.
#ifdef __CINT__
  Char_t trackparams[] = "SetTrackingParameters_4000bf04.C"; //Set this to correspond 
                                                             //with mult. and BField
#endif
  
  //for aliroot the path should point to a file 
  //containing the tpc geometry called alirunfile.root
  Bool_t isinit=AliHLTTransform::Init(path,(filetype!=AliLevel3::kBinary));
  if(!isinit){
    cerr << "Could not create transform settings, please check log for error messages!" << endl;
    return;
  }

  for(Int_t ev=0; ev<nevent; ev++)
    {
      AliLevel3 *a;
      if(filetype==AliLevel3::kBinary)
	a = new AliLevel3();
      else 
	{
	  Char_t fname[1024];
	  if(filetype==AliLevel3::kRaw)
           sprintf(fname,"%s/raw.root",path);
	  else
           sprintf(fname,"%s/digitfile.root",path);
          a = new AliLevel3(fname);
	}
      
      a->Init(path,filetype,npatches);
      
#ifdef __CINT__
      gROOT->LoadMacro(trackparams);
      SetTrackingParameters(a);
#else /* compiled for 4000 and 0.4 */
      Int_t phi_segments,eta_segments,trackletlength,tracklength;
      Int_t rowscopetracklet,rowscopetrack;
      Double_t min_pt_fit,maxangle,goodDist,hitChi2Cut,xyerror,zerror;
      Double_t goodHitChi2,trackChi2Cut,maxphi,maxeta;

      phi_segments = 50;    //devide the space into phi_segments and eta_segments
      eta_segments = 100;   //to access the search of points to that area!
      trackletlength = 3;   //number of hits a tracklet has to have
      tracklength = 10;     //number of hits a track has to have 
      rowscopetracklet = 2; //search range of rows for a tracklet
      rowscopetrack = 10;   //search range of rows for a track
      min_pt_fit = 0;      
      maxangle = 0.1745;   //AliHLTTransform::Deg2Rad(10);
                           //maximum angle for the three point look ahead
      goodDist = 5;        //threshold distance between two hits when building tracklets
      maxphi=0.1;          //maximum phi difference for neighboring hits
      maxeta=0.1;          //maximum eta difference for neighboring hits
      hitChi2Cut = 20;     //maximum chi2 of added hit to track
      goodHitChi2 = 5;     //stop looking for next hit to add if chi2 is less than goodHitChi2
      trackChi2Cut = 10;   //maximum chi2 for track after final fit
      xyerror = -1;
      zerror =  -1;
  
      a->SetClusterFinderParam(xyerror,zerror,kTRUE);
      a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
			   rowscopetracklet,rowscopetrack,
			   min_pt_fit,maxangle,goodDist,hitChi2Cut,
			   goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kTRUE);
      
      if(nonvertex)
	{
	  //Set parameters for nonvertextracking
	  a->SetTrackerParam(phi_segments,eta_segments,trackletlength,tracklength,
			     rowscopetracklet,rowscopetrack,
			     min_pt_fit,maxangle,goodDist,hitChi2Cut,
			     goodHitChi2,trackChi2Cut,50,maxphi,maxeta,kFALSE);
	}
#endif

      if(pileup)
	a->DoPileup();
      //a->DoRoi();    /*do region of interest*/
      //a->DoMc();     /*do monte carlo identification*/
      
      if(nonvertex)
	a->DoNonVertexTracking(); /*2 tracking passes, last without vertex contraint.*/
      
      a->WriteFiles(opath); /*enable output*/
      a->ProcessEvent(minslice,maxslice);
      Char_t bname[100];
      sprintf(bname,"benchmark_tracker_%d",ev);
      a->DoBench(bname);
      delete a;
    } // event loop
}


