#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TParticle.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliHLTTriggerDecision.h"

#include "AliESDCaloCluster.h"
#endif

void CreateInput(const char* filename, Int_t numOfTracks, Double_t minPt, Double_t maxPt)
{
  gRandom->SetSeed(123);
  TFile* file = new TFile(filename, "RECREATE");
  AliESDEvent event;
  event.CreateStdContent();

  AliESDCaloCluster cluster;
  cluster = new AliESDCaloCluster();
  
  Double_t et = gRandom->Rndm() * (maxPt - minPt) + minPt;
  
  cluster.SetE(et);
  cluster.SetClusterType(AliESDCaloCluster::kPHOSCluster);
  //cluster.SetClusterType(nuOfmTracks);

  event.AddCaloCluster(&cluster);

  event.Write();
  delete file;
}

bool CheckIfOutputIsOk()
{
  const char* filename = "PhosClusterEnergyTriggerTestOutput.root";
  TFile file(filename, "READ");
  AliHLTTriggerDecision* decision1 = dynamic_cast<AliHLTTriggerDecision*>(file.Get("PhosClusterEnergyTrigger;2"));
  AliHLTTriggerDecision* decision2 = dynamic_cast<AliHLTTriggerDecision*>(file.Get("PhosClusterEnergyTrigger;3"));
  if (decision1 == NULL)
  {
    cerr << "ERROR: 'PhosClusterEnergyTrigger;1' AliHLTTriggerDecision object not found in file " << filename << endl;
    return false;
  }
  if (decision2 == NULL)
  {
    cerr << "ERROR: 'PhosClusterEnergyTrigger;2' AliHLTTriggerDecision object not found in file " << filename << endl;
    return false;
  }
  cout << "============================== Decision 1 ==============================" << endl;
  decision1->Print();
  if (decision1->Result() != 0)
  {
    cout << "FAILED result check. Expected a result of 0 for trigger decision 1 but received: " << decision1->Result() << endl;
    return false;
  }
  cout << "============================== Decision 2 ==============================" << endl;
  decision2->Print();
  if (decision2->Result() != 1)
  {
    cout << "FAILED result check. Expected a result of 1 for trigger decision 2 but received: " << decision2->Result() << endl;
    return false;
  }

  return true;
}

bool testPhosClusterEnergyTrigger()
{
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
  gSystem->Load("libAliHLTUtil.so");
  gSystem->Load("libAliHLTMUON.so");
  gSystem->Load("libAliHLTTRD.so");
  gSystem->Load("libAliHLTTrigger.so");
  CreateInput("PhosClusterEnergyTriggerTestInput1.root", -2, 0.1, 0.99);
  CreateInput("PhosClusterEnergyTriggerTestInput2.root", 0, 2.1, 4.0);
  AliHLTSystem sys;
  sys.LoadComponentLibraries("libAliHLTUtil.so");
  sys.LoadComponentLibraries("libAliHLTTrigger.so");
  const char* cmdline = " -datatype ROOTTOBJ 'HLT ' -datafile PhosClusterEnergyTriggerTestInput1.root -nextevent -datafile PhosClusterEnergyTriggerTestInput2.root";
  AliHLTConfiguration pub("pub", "ROOTFilePublisher", NULL, cmdline);
  AliHLTConfiguration proc("proc", "PhosClusterEnergyTrigger", "pub", "");
  AliHLTConfiguration sink("sink", "ROOTFileWriter", "proc", "-datafile PhosClusterEnergyTriggerTestOutput.root -concatenate-events");
  sys.BuildTaskList("sink");
  sys.Run(2);
  return CheckIfOutputIsOk();
}
#ifndef __MAKECINT__
int main(int /*argc*/, const char** /*argv*/)
{
  bool resultOk = testPhosClusterEnergyTrigger();
  if (not resultOk) return 1;
  return 0;
}
#endif // __MAKECINT__
