/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#ifndef __CINT__
  #include <Riostream.h>
  #include "AliTPCtracker.h"

  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliTPCrefit() {
  
  cout << "Propagating tracks inward through the TPC..." << endl;
  
  AliKalmanTrack::SetConvConst(1000/0.299792458/4.);
  
  TFile *in = new TFile("AliTPCBackTracks.root", "READ");
  TFile *out= new TFile("AliTPCrefited.root", "RECREATE");   
  TFile *file = new TFile("AliTPCclusters.root");
  
  AliTPCParam *param=(AliTPCParam*)file->Get("75x40_100x60_150x60");
  if (!param) {cerr<<"Can't get TPC parameters !\n"; return 4;}
  
  TStopwatch timer;
  AliTPCtracker *tracker = new AliTPCtracker(param);
  Int_t rc=tracker->RefitInward(in,out);
  delete tracker;
  timer.Stop(); timer.Print();
   
  file->Close();
  in->Close();
  out->Close();
  
  return rc;
}
