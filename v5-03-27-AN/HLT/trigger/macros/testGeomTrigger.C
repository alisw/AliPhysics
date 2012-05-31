#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Riostream.h"
#include "TFile.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TParticle.h"
#include "AliESDEvent.h"
#include "AliESDMuonTrack.h"
#include "AliHLTSystem.h"
#include "AliHLTConfiguration.h"
#include "AliHLTTriggerDecision.h"

#endif

void CreateTrack(AliESDMuonTrack& track, Double_t minPt, Double_t maxPt, Double_t charge)
{
  Double_t pz = -100;
  Double_t pt = gRandom->Rndm() * (maxPt - minPt) + minPt;
  Double_t angle = gRandom->Rndm() * 2 * TMath::Pi();
  Double_t px = pt * TMath::Cos(angle);
  Double_t py = pt * TMath::Sin(angle);
  track.SetThetaX(TMath::ATan2(px, pz));
  track.SetThetaY(TMath::ATan2(py, pz));
  Double_t invP = 1. / TMath::Sqrt(py*py + pz*pz) * charge;
  track.SetInverseBendingMomentum(invP);
}

TParticle CreateParticle(Double_t minPt, Double_t maxPt)
{
  Double_t mass = 0.1056583699; // [GeV/c] muon mass.
  Double_t pz = 0;
  Double_t pt = gRandom->Rndm() * (maxPt - minPt) + minPt;
  Double_t angle = gRandom->Rndm() * 2 * TMath::Pi();
  //  Double_t angle = 270*TMath::DegToRad();
  Double_t px = pt * TMath::Cos(angle);
  Double_t py = pt * TMath::Sin(angle);
  Double_t etot = TMath::Sqrt(pt*pt + pz*pz + mass*mass);
  return TParticle(13, 0, 0, 0, 0, 0, px, py, pz, etot, 0, 0, 0, 0);
}

void CreateInput(
    const char* filename,
    Int_t numOfBarrelTracks, Double_t minBarrelPt, Double_t maxBarrelPt,
    Int_t numOfMuonTracks, Double_t minMuonPt, Double_t maxMuonPt
  )
{
  TFile* file = new TFile(filename, "RECREATE");
  AliESDEvent event;
  event.CreateStdContent();
  for (int i = 0; i < numOfBarrelTracks; i++)
  {
    TParticle p = CreateParticle(minBarrelPt, maxBarrelPt);
    AliESDtrack track(&p);
    event.AddTrack(&track);
  }
  for (int i = 0; i < numOfMuonTracks; i++)
  {
    AliESDMuonTrack *mutrack = event.NewMuonTrack();
    CreateTrack(*mutrack, minMuonPt, maxMuonPt, (i%2==0) ? -1 : 1);
  }
  event.Write();
  delete file;
}

bool CheckIfOutputIsOk()
{
  const char* filename = "FullTriggerTestOutput.root";
  TFile file(filename, "READ");
  int expectedResult[4] = {0, 1, 0, 1};
  AliHLTTriggerDecision* decision[4];
  for (int i = 0; i < 4; i++)
  {
    char name[1024];
    sprintf(name, "HLTGlobalTrigger;%d", i+1);
    decision[i] = dynamic_cast<AliHLTTriggerDecision*>(file.Get(name));
    if (decision[i] == NULL)
    {
      cerr << "ERROR: '" << name << "' AliHLTTriggerDecision object not found in file " << filename << endl;
      return false;
    }
    cout << "===============================================================================" << endl;
    cout << "============================== Global Decision " << i+1 << " ==============================" << endl;
    cout << "===============================================================================" << endl;
    decision[i]->Print();
    if (decision[i]->Result() != expectedResult[i])
    {
      cout << "FAILED result check. Expected a result of " << expectedResult[i]
           << " for trigger decision " << i+1
           << " but received: " << decision[i]->Result()
           << endl;
      return false;
    }
  }
  return true;
}

bool testGeomTrigger()
{
  gRandom->SetSeed(123);
  gSystem->Load("libAliHLTTrigger.so");
  CreateInput("FullTriggerTestInput1.root", 5, 0.1, 1.9, 0, 0.1, 0.9);
  CreateInput("FullTriggerTestInput2.root", 3, 2.1, 4.0, 0, 0.1, 0.9);
  CreateInput("FullTriggerTestInput3.root", 6, 0.1, 1.9, 0, 1.1, 2.0);
  CreateInput("FullTriggerTestInput4.root", 2, 2.1, 4.0, 0, 1.1, 2.0);
  AliHLTSystem sys;
  sys.LoadComponentLibraries("libAliHLTUtil.so");
  sys.LoadComponentLibraries("libAliHLTTrigger.so");
  const char* cmdline = " -datatype ROOTTOBJ 'HLT ' -datafile FullTriggerTestInput1.root"
    " -nextevent -datafile FullTriggerTestInput2.root"
    " -nextevent -datafile FullTriggerTestInput3.root"
    " -nextevent -datafile FullTriggerTestInput4.root";
  AliHLTConfiguration pub("pub", "ROOTFilePublisher", NULL, cmdline);
  AliHLTConfiguration trg1("trg1", "BarrelGeomMultiplicityTrigger", "pub", "-triggername PHOSgeomtrigger -geomfile PHOSgeomtrigger.root");
  AliHLTConfiguration glob("global", "HLTGlobalTrigger", "trg1", "-config testGeomTriggerConfig.C");
  AliHLTConfiguration sink("sink", "ROOTFileWriter", "global", "-datafile FullTriggerTestOutput.root -concatenate-events");
  sys.BuildTaskList("sink");
  sys.Run(4);
  return CheckIfOutputIsOk();
}
