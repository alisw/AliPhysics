/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */
//_________________________________________________________________________
//Basic Memory Leak utility.    
//     You can use this tiny class to *see* if your program is leaking.
//     Usage:
//     AliMemoryWatcher memwatcher;
//     some program loop on events here {
//       if ( nevents % x == 0 ) 
//       {
//       // take a sample every x events
//         memwatcher.Watch(nevents);
//       }
//     }
//     TFile f("out.root","RECREATE");
//     memwatcher.Write();
//     f.Close();
//     In the output root file you'll get 3 graphs representing
//     the evolAliPHOSon, as a function of the number of events, of :
//     - VSIZE is the virtual size (in KBytes) of your program, that is sort of
//     the total memory used
//     - RSSIZE is the resident size (in KBytes), that is, the part of your 
//     program which is really in physical memory.
//     - TIME is an estimate of time per event (really it's the time elasped
//     between two calls to watch method)
//     WARNING: this is far from a bulletproof memory report (it's basically 
//     using UNIX command ps -h -p [PID] -o vsize,rssize to do its job).
//     It has only been tested on Linux so far.    
//     But by fitting the VSIZE by a pol1 under ROOT, you'll see right away
//     by how much your program is leaking.          
//*-- Author: Laurent Aphecetche(SUBATECH)
// --- std system ---
#include <cassert> 
#ifdef __APPLE__
#include <stdlib.h>
#else
#include <malloc.h>
#endif
// --- AliRoot header files ---
#include "AliLog.h"
#include "AliMemoryWatcher.h"
// --- ROOT system ---
#include "TSystem.h"
#include "TGraph.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TError.h"

ClassImp(AliMemoryWatcher)

//_____________________________________________________________________________
AliMemoryWatcher::AliMemoryWatcher(UInt_t maxsize) :
  TObject(),
  fUseMallinfo(kTRUE),
  fPID(gSystem->GetPid()),
  fMAXSIZE(maxsize),
  fSize(0),
  fX(new Int_t[fMAXSIZE]),
  fVSIZE(new Int_t[fMAXSIZE]),
  fRSSIZE(new Int_t[fMAXSIZE]),
  fTIME(new Double_t[fMAXSIZE]),
  fTimer(0),
  fDisabled(kFALSE)
{
  //
  //ctor
  //
  sprintf(fCmd,"ps -h -p %d -o vsize,rssize",fPID);
}

//_____________________________________________________________________________
AliMemoryWatcher::AliMemoryWatcher(const AliMemoryWatcher& mw):
  TObject(mw),
  fUseMallinfo(mw.fUseMallinfo),
  fPID(mw.fPID),
  fMAXSIZE(mw.fMAXSIZE),
  fSize(0),
  fX(new Int_t[fMAXSIZE]),
  fVSIZE(new Int_t[fMAXSIZE]),
  fRSSIZE(new Int_t[fMAXSIZE]),
  fTIME(new Double_t[fMAXSIZE]),
  fTimer(0),
  fDisabled(kFALSE)
{
  //copy ctor
  strcpy(fCmd, mw.fCmd) ; 
}

//_____________________________________________________________________________
AliMemoryWatcher::~AliMemoryWatcher()
{
  // dtor
  delete[] fVSIZE;
  delete[] fRSSIZE;
  delete[] fX;
  delete[] fTIME;
  delete fTimer;
}
//_____________________________________________________________________________
void AliMemoryWatcher::Watch(Int_t x)
{
  // Sets the point where CPU parameters have to be monitored
  if ( !fDisabled && fSize < fMAXSIZE ) {
    if ( fSize==0 ) {
      assert(fTimer==0);
      fTimer = new TStopwatch;
      fTimer->Start(true);
      fTimer->Stop();
    }
    if ( fUseMallinfo ) {
#ifdef __linux
      static struct mallinfo meminfo;
      meminfo = mallinfo();
      fX[fSize] = x ;
      fVSIZE[fSize] = (meminfo.hblkhd + meminfo.uordblks) / 1024;
      fRSSIZE[fSize] =  meminfo.uordblks / 1024;
      fTIME[fSize] = fTimer->CpuTime();
      fSize++;
#else
      AliFatal("Please SetUseMallinfo to kFALSE on this system");
#endif
    } else {
      static Int_t vsize, rssize;
      static FILE* pipe = 0;
      pipe = popen(fCmd,"r");
      if ( pipe ) {
    
	fscanf(pipe,"%d %d",&vsize,&rssize);
      
	fX[fSize] = x ;
	fVSIZE[fSize] = vsize ;
	fRSSIZE[fSize] = rssize ;
	fTIME[fSize] = fTimer->CpuTime();
	fSize++;
      }
      assert(pclose(pipe)!=-1);
    }
    fTimer->Start(true);
  }
  else {
    fDisabled=true;
    AliError("I'm full !" ) ;
  }
}
//_____________________________________________________________________________
TGraph*
AliMemoryWatcher::GraphVSIZE(void)
{
  // Fills the graph with the virtual memory sized used
  TGraph* g = 0;
  if ( Size() )
    {
      g = new TGraph(Size());
      Int_t i ; 
      for (i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),VSIZE(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TGraph*
AliMemoryWatcher::GraphRSSIZE(void)
{
  // Fills the graph with the real memory sized used
  TGraph* g = 0;
  if ( Size() ) 
    {
      g = new TGraph(Size());
      Int_t i ; 
      for (i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),RSSIZE(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TGraph*
AliMemoryWatcher::GraphTIME(void)
{
  // Fills the raph with the used CPU time
  TGraph* g = 0;
  if ( Size() ) 
    {
      g = new TGraph(Size());
      Int_t i ; 
      for (i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),TIME(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TH2*
AliMemoryWatcher::Frame(void) const
{
  //creates the frame histo in which the graphs will be plotted 
  Double_t xmin=1E30;
  Double_t xmax=0;
  Double_t ymin=1;
  Double_t ymax=0;
  UInt_t i ; 
  for (i=0; i < Size() ; i++ ) {
    if ( X(i) < xmin ) xmin = X(i);
    if ( X(i) > xmax ) xmax = X(i);
    Double_t y = VSIZE(i)+RSSIZE(i);
    if ( y > ymax ) ymax = y;
    if ( VSIZE(i) < ymin ) ymin = VSIZE(i);
    if ( RSSIZE(i) < ymin ) ymin = RSSIZE(i);
  }
  TH2F* h = new TH2F("frame","",10,xmin,xmax,10,ymin*0.8,ymax*1.2);
  return h;
}
//_____________________________________________________________________________
Int_t
AliMemoryWatcher::WriteToFile()
{
  // Stores the graphs in a file 
  if ( GraphVSIZE() ) GraphVSIZE()->Write("VSIZE",TObject::kOverwrite);
  if ( GraphRSSIZE() ) GraphRSSIZE() ->Write("RSSIZE",TObject::kOverwrite);
  if ( GraphTIME() ) GraphTIME()->Write("TIME",TObject::kOverwrite);
  return 0;
}
