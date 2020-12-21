///
/// \file AliFemto/AliFemtoCutMonitorCollections.cxx
///

#include "AliFemtoCutMonitorCollections.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TString.h>

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections():
  AliFemtoCutMonitorCollections("")  // construct with empty 'name'
{
}

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections(const char *aName):
  AliFemtoCutMonitor(),
  fCollection1Mult(nullptr),
  fCollection2Mult(nullptr)
{
  // Normal constructor
  fCollection1Mult = new TH1D(TString("Coll1Mult") + aName, "Collection 1 Multiplicity", 5001, -0.5, 5000.5);
  fCollection2Mult = new TH1D(TString("Coll2Mult") + aName, "Collection 2 Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections(const AliFemtoCutMonitorCollections &aCutMonitor):
  AliFemtoCutMonitor(),
  fCollection1Mult(nullptr),
  fCollection2Mult(nullptr)
{
  // copy constructor
  fCollection1Mult = new TH1D(*aCutMonitor.fCollection1Mult);
  fCollection2Mult = new TH1D(*aCutMonitor.fCollection2Mult);
}

AliFemtoCutMonitorCollections::~AliFemtoCutMonitorCollections()
{
  // Destructor
  delete fCollection1Mult;
  delete fCollection2Mult;
}

AliFemtoCutMonitorCollections& AliFemtoCutMonitorCollections::operator=(const AliFemtoCutMonitorCollections& aCutMonitor)
{
  // assignment operator
  if (this != &aCutMonitor) {
    *fCollection1Mult = *aCutMonitor.fCollection1Mult;
    *fCollection2Mult = *aCutMonitor.fCollection2Mult;
  }
  return *this;
}

AliFemtoString AliFemtoCutMonitorCollections::Report()
{
  // Prepare report from the execution
  string stemp = "*** AliFemtoCutMonitorCollections report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}

void AliFemtoCutMonitorCollections::Fill(const AliFemtoParticleCollection* aCollection1,const AliFemtoParticleCollection* aCollection2)
{
  // Fill in the monitor histograms with the values from the current event
  //cout<<"Monitor collection sizes: "<<aCollection1->size()<<" "<<aCollection2->size()<<endl;
  fCollection1Mult->Fill(aCollection1->size());
  fCollection2Mult->Fill(aCollection2->size());
}

void AliFemtoCutMonitorCollections::Write()
{
  // Write out the relevant histograms
  fCollection1Mult->Write();
  fCollection2Mult->Write();
}

TList *AliFemtoCutMonitorCollections::GetOutputList()
{
  TList *tOutputList = new TList();
  tOutputList->Add(fCollection1Mult);
  tOutputList->Add(fCollection2Mult);

  return tOutputList;
}


