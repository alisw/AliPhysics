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
/*Basic Memory Leak utility.
    
    You can use this tiny class to *see* if your program is leaking.
    Usage:
    AliPHOSMemoryWatcher memwatcher;
    some program loop on events here {
      if ( nevents % x == 0 ) 
      {
      // take a sample every x events
        memwatcher.watch(nevents);
      }
    }
    TFile f("out.root","RECREATE");
    memwatcher.write();
    f.Close();
    In the output root file you'll get 3 graphs representing
    the evolAliPHOSon, as a function of the number of events, of :
    - VSIZE is the virtual size (in KBytes) of your program, that is sort of
    the total memory used
    - RSSIZE is the resident size (in KBytes), that is, the part of your 
    program which is really in physical memory.
    - TIME is an estimate of time per event (really it's the time elasped
    between two calls to watch method)
    WARNING: this is far from a bulletproof memory report (it's basically 
    using UNIX command ps -h -p [PID] -o vsize,rssize to do its job).
    It has only been tested on Linux so far.
    
    But by fitting the VSIZE by a pol1 under ROOT, you'll see right away
    by how much your program is leaking.
*/             
//*-- Author: Laurent Aphecetche(SUBATECH)
// --- std system ---
#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// --- AliRoot header files ---
#include "AliPHOSMemoryWatcher.h"
// --- ROOT system ---
#include "TSystem.h"
#include "TGraph.h"
#include "TH2.h"
#include "TStopwatch.h"
//_____________________________________________________________________________
AliPHOSMemoryWatcher::AliPHOSMemoryWatcher(unsigned int maxsize)
{
  fMAXSIZE=maxsize;
  fPID = gSystem->GetPid();
  sprintf(fCmd,"ps -h -p %d -o vsize,rssize",fPID);
  fX = new int[fMAXSIZE];
  fVSIZE = new int[fMAXSIZE];
  fRSSIZE = new int[fMAXSIZE];
  fTIME = new double[fMAXSIZE];
  fSize=0;
  fDisabled=false;
  fTimer=0;
}
//_____________________________________________________________________________
AliPHOSMemoryWatcher::~AliPHOSMemoryWatcher()
{
  delete[] fVSIZE;
  delete[] fRSSIZE;
  delete[] fX;
  delete[] fTIME;
  delete fTimer;
}
//_____________________________________________________________________________
void AliPHOSMemoryWatcher::watch(int x)
{
  if ( !fDisabled && fSize < fMAXSIZE ) {
    if ( fSize==0 ) {
      assert(fTimer==0);
      fTimer = new TStopwatch;
      fTimer->Start(true);
      fTimer->Stop();
    }
    static int vsize, rssize;
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
    fTimer->Start(true);
  }
  else {
    fDisabled=true;
    cerr << "AliPHOSMemoryWatcher::watch : I'm full !" << endl;
  }
}
//_____________________________________________________________________________
TGraph*
AliPHOSMemoryWatcher::graphVSIZE(void)
{
  TGraph* g = 0;
  if ( size() )
    {
      g = new TGraph(size());
      for (int i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),VSIZE(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TGraph*
AliPHOSMemoryWatcher::graphRSSIZE(void)
{
  TGraph* g = 0;
  if ( size() ) 
    {
      g = new TGraph(size());
      for (int i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),RSSIZE(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TGraph*
AliPHOSMemoryWatcher::graphTIME(void)
{
  TGraph* g = 0;
  if ( size() ) 
    {
      g = new TGraph(size());
      for (int i=0; i < g->GetN(); i++ ) {
        g->SetPoint(i,X(i),TIME(i));
      }
    }
  return g;
}
//_____________________________________________________________________________
TH2*
AliPHOSMemoryWatcher::frame(void)
{
  double xmin=1E30;
  double xmax=0;
  double ymin=1;
  double ymax=0;
  for (unsigned int i=0; i < size() ; i++ ) {
    if ( X(i) < xmin ) xmin = X(i);
    if ( X(i) > xmax ) xmax = X(i);
    double y = VSIZE(i)+RSSIZE(i);
    if ( y > ymax ) ymax = y;
    if ( VSIZE(i) < ymin ) ymin = VSIZE(i);
    if ( RSSIZE(i) < ymin ) ymin = RSSIZE(i);
  }
  TH2F* h = new TH2F("frame","",10,xmin,xmax,10,ymin*0.8,ymax*1.2);
  return h;
}
//_____________________________________________________________________________
void 
AliPHOSMemoryWatcher::write(void)
{
  if ( graphVSIZE() ) graphVSIZE()->Write("VSIZE",TObject::kOverwrite);
  if ( graphRSSIZE() ) graphRSSIZE()->Write("RSSIZE",TObject::kOverwrite);
  if ( graphTIME() ) graphTIME()->Write("TIME",TObject::kOverwrite);
}
