////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorCollections - the cut monitor for particles to study    //
// the difference between reconstructed and true momentum                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorCollections.h"
#include "AliFemtoModelHiddenInfo.h"
#include "AliFemtoEvent.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections():
  fCollection1Mult(0),
  fCollection2Mult(0)
{
  // Default constructor
  fCollection1Mult = new TH1D("Coll1Mult", "Collection 1 Multiplicity", 5001, -0.5, 5000.5);
  fCollection2Mult = new TH1D("Coll2Mult","Collection 2 Multiplicity",5001,-0.5,5000.5);
}

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections(const char *aName):
  fCollection1Mult(0),
  fCollection2Mult(0)
{
  // Normal constructor
  char name[200];
  snprintf(name, 200, "Coll1Mult%s", aName);
  fCollection1Mult = new TH1D(name, "Collection 1 Multiplicity", 5001, -0.5, 5000.5);

  snprintf(name, 200, "Coll2Mult%s", aName);
  fCollection2Mult = new TH1D(name, "Collection 2 Multiplicity", 5001, -0.5, 5000.5);
}

AliFemtoCutMonitorCollections::AliFemtoCutMonitorCollections(const AliFemtoCutMonitorCollections &aCut):
  fCollection1Mult(0),
  fCollection2Mult(0)
{
  // copy constructor
  if (fCollection1Mult) delete fCollection1Mult;
  fCollection1Mult = new TH1D(*aCut.fCollection1Mult);

  if (fCollection2Mult) delete fCollection2Mult;
  fCollection2Mult = new TH1D(*aCut.fCollection2Mult);
}

AliFemtoCutMonitorCollections::~AliFemtoCutMonitorCollections()
{
  // Destructor
  delete fCollection1Mult;
  delete fCollection2Mult;
}

AliFemtoCutMonitorCollections& AliFemtoCutMonitorCollections::operator=(const AliFemtoCutMonitorCollections& aCut)
{
  // assignment operator
  if (this == &aCut) 
    return *this;

  if (fCollection1Mult) delete fCollection1Mult;
  fCollection1Mult = new TH1D(*aCut.fCollection1Mult);

  if (fCollection2Mult) delete fCollection2Mult;
  fCollection2Mult = new TH1D(*aCut.fCollection2Mult);

  return *this;
}

AliFemtoString AliFemtoCutMonitorCollections::Report(){ 
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


