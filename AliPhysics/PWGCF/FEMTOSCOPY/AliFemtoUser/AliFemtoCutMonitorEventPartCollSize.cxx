////////////////////////////////////////////////////////////////////////////////
// AliFemtoCutMonitorEventPartCollSize                                        //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorEventVertex.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

#include "AliFemtoCutMonitorEventPartCollSize.h"

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCutMonitorEventPartCollSize);
  /// \endcond
#endif

//___________________________________________________
AliFemtoCutMonitorEventPartCollSize::AliFemtoCutMonitorEventPartCollSize():
  AliFemtoCutMonitor(),
  fCollSizePerEvent1(0),
  fCollSizePerEvent2(0),
  fCollSizePerEvent1vsEvent2(0)
{
  //Default constructor
  fCollSizePerEvent1 = new TH1D("MultPerEventPart1", "Part1 Mult. Per Event", 100, 0, 1000);
  fCollSizePerEvent2 = new TH1D("MultPerEventPart2", "Part2 Mult. Per Event", 100, 0, 1000);
  fCollSizePerEvent1vsEvent2 = new TH2D("MultPerEventPart1vsPart2", "Part1 vs Part2 Mult. Per Event", 100, 0, 1000, 100, 0, 1000);

  fCollSizePerEvent1->Sumw2();
  fCollSizePerEvent2->Sumw2();
  fCollSizePerEvent1vsEvent2->Sumw2();
}


//___________________________________________________
AliFemtoCutMonitorEventPartCollSize::AliFemtoCutMonitorEventPartCollSize(const char *aName1, const int& aNBins1, const int& aNLow1, const int& aNHigh1, const char *aName2, const int& aNBins2, const int& aNLow2, const int& aNHigh2):
  AliFemtoCutMonitor(),
  fCollSizePerEvent1(0),
  fCollSizePerEvent2(0),
  fCollSizePerEvent1vsEvent2(0)
{
  //Normal constructor
  char name[200];
  snprintf(name, 200, "MultPerEvent%s", aName1);
  fCollSizePerEvent1 = new TH1D(name, "Part1 Mult. Per Event", aNBins1, aNLow1, aNHigh1);
  snprintf(name, 200, "MultPerEvent%s", aName2);
  fCollSizePerEvent2 = new TH1D(name, "Part2 Mult. Per Event", aNBins2, aNLow2, aNHigh2);
  snprintf(name, 200, "MultPerEvent%svs%s", aName1, aName2);
  fCollSizePerEvent1vsEvent2 = new TH2D(name, "Part1 vs Part2 Mult. Per Event", aNBins1, aNLow1, aNHigh1, aNBins2, aNLow2, aNHigh2);

  fCollSizePerEvent1->Sumw2();
  fCollSizePerEvent2->Sumw2();
  fCollSizePerEvent1vsEvent2->Sumw2();
}

//___________________________________________________
AliFemtoCutMonitorEventPartCollSize::~AliFemtoCutMonitorEventPartCollSize()
{
  //Destructor
  delete fCollSizePerEvent1;
  delete fCollSizePerEvent2;
  delete fCollSizePerEvent1vsEvent2;
}

//___________________________________________________
AliFemtoCutMonitorEventPartCollSize::AliFemtoCutMonitorEventPartCollSize(const AliFemtoCutMonitorEventPartCollSize& aCut):
  AliFemtoCutMonitor()
{
  //copy constructor
  if(fCollSizePerEvent1) fCollSizePerEvent1 = new TH1D(*aCut.fCollSizePerEvent1);
  else fCollSizePerEvent1 = 0;

  if(fCollSizePerEvent2) fCollSizePerEvent2 = new TH1D(*aCut.fCollSizePerEvent2);
  else fCollSizePerEvent2 = 0;

  if(fCollSizePerEvent1vsEvent2) fCollSizePerEvent1vsEvent2 = new TH2D(*aCut.fCollSizePerEvent1vsEvent2);
  else fCollSizePerEvent1vsEvent2 = 0;
}

//___________________________________________________
AliFemtoCutMonitorEventPartCollSize& AliFemtoCutMonitorEventPartCollSize::operator=(const AliFemtoCutMonitorEventPartCollSize& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  if(fCollSizePerEvent1) fCollSizePerEvent1 = new TH1D(*aCut.fCollSizePerEvent1);
  else fCollSizePerEvent1 = 0;

  if(fCollSizePerEvent2) fCollSizePerEvent2 = new TH1D(*aCut.fCollSizePerEvent2);
  else fCollSizePerEvent2 = 0;

  if(fCollSizePerEvent1vsEvent2) fCollSizePerEvent1vsEvent2 = new TH2D(*aCut.fCollSizePerEvent1vsEvent2);
  else fCollSizePerEvent1vsEvent2 = 0;

  return *this;
}

//___________________________________________________
AliFemtoString AliFemtoCutMonitorEventPartCollSize::Report()
{
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorEventPartCollSize report"; 
  AliFemtoString returnThis = stemp;
  return returnThis; 
}

//___________________________________________________
void AliFemtoCutMonitorEventPartCollSize::Fill(const AliFemtoParticleCollection *aCollection1, const AliFemtoParticleCollection *aCollection2)
{
  fCollSizePerEvent1->Fill(aCollection1->size());
  fCollSizePerEvent2->Fill(aCollection2->size());
  fCollSizePerEvent1vsEvent2->Fill(aCollection1->size(),aCollection2->size());
}

//___________________________________________________
void AliFemtoCutMonitorEventPartCollSize::Write()
{
  fCollSizePerEvent1->Write();
  fCollSizePerEvent2->Write();
  fCollSizePerEvent1vsEvent2->Write();
}

//___________________________________________________
TList* AliFemtoCutMonitorEventPartCollSize::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fCollSizePerEvent1);
  tOutputList->Add(fCollSizePerEvent2);
  tOutputList->Add(fCollSizePerEvent1vsEvent2);

  return tOutputList;
}
  
