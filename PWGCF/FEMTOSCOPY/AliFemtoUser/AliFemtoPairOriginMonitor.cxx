#include "AliFemtoPairOriginMonitor.h"
#include <TH1D.h>
#include <TH3D.h>
#include <TList.h>
#include "AliFemtoModelHiddenInfo.h"

AliFemtoPairOriginMonitor::AliFemtoPairOriginMonitor():
  fParticle1Origin(0),
  fParticle2Origin(0),
  fParticle1Id(0),
  fParticle2Id(0)
{
  // Default constructor

  fParticle1Origin =  new TH1D("POriginFirst", "Mothers PDG Codes", 6000, 0.0, 6000.0);
  fParticle1Id =  new TH1D("PIdFirst", "Particle PDG Codes", 6000, 0.0, 6000.0);
  fParticle2Origin =  new TH1D("POriginFirst", "Mothers PDG Codes", 6000, 0.0, 6000.0);
  fParticle2Id =  new TH1D("PIdFirst", "Particle PDG Codes", 6000, 0.0, 6000.0);

}
AliFemtoPairOriginMonitor::AliFemtoPairOriginMonitor(const char *aName):
  AliFemtoCutMonitor(),
  fParticle1Origin(0),
  fParticle2Origin(0),
  fParticle1Id(0),
  fParticle2Id(0)
{
  // Normal constructor
  char name[200];

  snprintf(name, 200, "POriginFirst%s", aName);
  fParticle1Origin =  new TH1D(name, "Mothers PDG Codes", 6000, 0.0, 6000.0);
  snprintf(name, 200, "PIdFirst%s", aName);
  fParticle1Id =  new TH1D(name, "Particle PDG Codes", 6000, 0.0, 6000.0);

  snprintf(name, 200, "POriginSecond%s", aName);
  fParticle2Origin =  new TH1D(name, "Mothers PDG Codes", 6000, 0.0, 6000.0);
  snprintf(name, 200, "PIdSecond%s", aName);
  fParticle2Id =  new TH1D(name, "Particle PDG Codes", 6000, 0.0, 6000.0);

}

AliFemtoPairOriginMonitor::AliFemtoPairOriginMonitor(const AliFemtoPairOriginMonitor &aCut):
  AliFemtoCutMonitor(),
  fParticle1Origin(0),
  fParticle2Origin(0),
  fParticle1Id(0),
  fParticle2Id(0)
{
  // copy constructor

  if (fParticle1Origin) delete fParticle1Origin;
  fParticle1Origin= new TH1D(*aCut.fParticle1Origin);
  if (fParticle1Id) delete fParticle1Id;
  fParticle1Id= new TH1D(*aCut.fParticle1Id);

  if (fParticle1Origin) delete fParticle2Origin;
  fParticle1Origin= new TH1D(*aCut.fParticle2Origin);
  if (fParticle1Id) delete fParticle2Id;
  fParticle1Id= new TH1D(*aCut.fParticle2Id);

}

void AliFemtoPairOriginMonitor::Fill(const AliFemtoPair* aPair) {

  AliFemtoParticle *first = (AliFemtoParticle*)aPair->Track1();
  AliFemtoParticle *second = (AliFemtoParticle*)aPair->Track2();

  if(first!=NULL && second!=NULL) {
    
    AliFemtoModelHiddenInfo *tInfo1 = (AliFemtoModelHiddenInfo*)first->GetHiddenInfo();
    AliFemtoModelHiddenInfo *tInfo2 = (AliFemtoModelHiddenInfo*)second->GetHiddenInfo();
    if(tInfo1!=NULL && tInfo2!=NULL) {
      Int_t part1ID = TMath::Abs(tInfo1->GetPDGPid());
      Int_t mother1ID = TMath::Abs(tInfo1->GetMotherPdgCode());
      
      Int_t part2ID = TMath::Abs(tInfo2->GetPDGPid());
      Int_t mother2ID = TMath::Abs(tInfo2->GetMotherPdgCode());      
      
      fParticle1Id->Fill(part1ID);
      fParticle1Origin->Fill(mother1ID);
      fParticle2Id->Fill(part2ID);
      fParticle2Origin->Fill(mother2ID);
    }
  }
}

AliFemtoPairOriginMonitor::~AliFemtoPairOriginMonitor()
{
  // Destructor

  delete fParticle1Origin;
  delete fParticle1Id;
  delete fParticle2Origin;
  delete fParticle2Id;
}

AliFemtoPairOriginMonitor& AliFemtoPairOriginMonitor::operator=(const AliFemtoPairOriginMonitor& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;


  if (fParticle1Origin) delete fParticle1Origin;
  fParticle1Origin= new TH1D(*aCut.fParticle1Origin);
  if (fParticle1Id) delete fParticle1Id;
  fParticle1Id= new TH1D(*aCut.fParticle1Id);

  if (fParticle2Origin) delete fParticle2Origin;
  fParticle2Origin= new TH1D(*aCut.fParticle2Origin);
  if (fParticle2Id) delete fParticle2Id;
  fParticle2Id= new TH1D(*aCut.fParticle2Id);


  return *this;
}

AliFemtoString AliFemtoPairOriginMonitor::Report(){
  // Prepare report from the execution
  string stemp = "*** AliFemtoPairOriginMonitor report";
  AliFemtoString returnThis = stemp;
  return returnThis;
}


void AliFemtoPairOriginMonitor::Write()
{
  // Write out the relevant histograms

  fParticle1Id->Write();
  fParticle1Origin->Write();
  fParticle2Id->Write();
  fParticle2Origin->Write();

}

TList *AliFemtoPairOriginMonitor::GetOutputList()
{
  // Get the list of histograms to write
  TList *tOutputList = new TList();

  tOutputList->Add(fParticle1Id);
  tOutputList->Add(fParticle1Origin);
  tOutputList->Add(fParticle2Id);
  tOutputList->Add(fParticle2Origin);

  return tOutputList;
}
