////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorEventNumber - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorEventNumber.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorEventNumber::AliFemtoCutMonitorEventNumber():
  fNevent(0)
{
  // Default constructor
}

AliFemtoCutMonitorEventNumber::AliFemtoCutMonitorEventNumber(const char *title, const int &zvtxBins=10, const double& zvtxmin=-10, const double& zvtxmax=10, const int &multBins=5, const int& multmin=0, const int& multmax=100):
  AliFemtoCutMonitor(),
  fNevent(0)
{
  // Normal constructor
 
  // set up numerator
  char tTitNevent[101] = "Nevent";
  strncat(tTitNevent,title, 100);
  fNevent = new TH2F(tTitNevent,title,multBins,multmin,multmax,zvtxBins,zvtxmin,zvtxmax);
  
  // to enable error bar calculation...
  fNevent->Sumw2();
}

AliFemtoCutMonitorEventNumber::AliFemtoCutMonitorEventNumber(const AliFemtoCutMonitorEventNumber &aCut):
  AliFemtoCutMonitor(),
  fNevent(0)
{
  // copy constructor
  fNevent= new TH2F(*aCut.fNevent);
  
}

AliFemtoCutMonitorEventNumber::~AliFemtoCutMonitorEventNumber()
{
  // Destructor
  delete fNevent;
}

AliFemtoCutMonitorEventNumber& AliFemtoCutMonitorEventNumber::operator=(const AliFemtoCutMonitorEventNumber& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fNevent) delete fNevent;
  fNevent = new TH2F(*aCut.fNevent);
    
  return *this;
}

AliFemtoString AliFemtoCutMonitorEventNumber::Report(){ 
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorEventNumber report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

void AliFemtoCutMonitorEventNumber::Fill(const AliFemtoEvent* aEvent) 
{
  // Fill in the monitor histograms with the values from the current track
  
  double mult = -1;
  double zvtx = -11;

  mult = aEvent->UncorrectedNumberOfPrimaries();

  AliFemtoThreeVector vtx = aEvent->PrimVertPos();
  zvtx = vtx.z();
  //cout<<mult<<" "<<zvtx<<endl;
  fNevent->Fill(mult,zvtx);  
 
}


void AliFemtoCutMonitorEventNumber::Write()
{
  // Write out the relevant histograms
  fNevent->Write();
}

TList *AliFemtoCutMonitorEventNumber::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fNevent);

  return tOutputList;
}
