#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TH1F.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TFile.h>
#include <Riostream.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TArrayF.h>
#include <TStyle.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TParticlePDG.h>
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliITSRecPointContainer.h"
#include "AliITSgeomTGeo.h"
#include "AliITSRecPoint.h"
#include "AliCDBManager.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
//

#include "../HLT/ITS/trackingSAP/AliITSSAPAux.h"
#include "../HLT/ITS/trackingSAP/AliITSSAPLayer.h"
#include "../HLT/ITS/trackingSAP/AliITSSAPTracker.h"
#include "TProfile.h"
#include "TStopwatch.h"
#endif


AliITSSAPTracker* tracker=0;

TTree* treeInp = 0;
AliRunLoader* runLoader = 0;
AliESDEvent* esd=0;

void ProcessEvent(int iev);
void Process(const char* path);
void ProcChunk(const char* path);
void TestTracker(TTree* tRP, const AliESDVertex* vtx);
void LoadClusters(TTree* tRP);
Double_t* DefLogAx(double xMn,double xMx, int nbin);
void InitTracker(int runNumber);

#ifdef _TIMING_
TProfile* hTiming[AliITSSAPTracker::kNSW];
#endif

void Process(const char* inpData)
{
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root") || inpDtStr.EndsWith(".zip#")) {
    ProcChunk(inpDtStr.Data());
  }
  else {
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return;
    }
    //
    inpDtStr.ReadLine(inpf);
    while ( !inpDtStr.IsNull() ) {
      inpDtStr = inpDtStr.Strip(TString::kBoth,' ');
      if (inpDtStr.BeginsWith("//") || inpDtStr.BeginsWith("#")) {inpDtStr.ReadLine(inpf); continue;}
      inpDtStr = inpDtStr.Strip(TString::kBoth,',');
      inpDtStr = inpDtStr.Strip(TString::kBoth,'"');
      ProcChunk(inpDtStr.Data());
      inpDtStr.ReadLine(inpf);
    }
  }
#ifdef _CONTROLH_
  tracker->SaveHistos();
#endif
}

void ProcChunk(const char* path)
{
  //
  TString dir=path;
  //
  printf("Processing %s\n",dir.Data());
  //
  if (dir.EndsWith("AliESDs.root")) {
    dir.ReplaceAll("AliESDs.root","");
  }
  //
  esd = new AliESDEvent();
  //
  runLoader = AliRunLoader::Open(Form("%sgalice.root",dir.Data()));
  if (!runLoader) {
    printf("galice not found\n");
    return;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints("ITS");
  // ESD
  TFile* esdFile = TFile::Open(Form("%sAliESDs.root",dir.Data()));
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file\n");
    runLoader->UnloadKinematics();
    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }
  treeInp = (TTree*) esdFile->Get("esdTree");
  if (!treeInp) {
    printf("Error: no ESD tree found\n");
    runLoader->UnloadKinematics();
    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }
  esd->ReadFromTree(treeInp);
  //
  for (int iEv=0;iEv<runLoader->GetNumberOfEvents(); iEv++) {
    printf("ev %d\n",iEv);
    ProcessEvent(iEv);
  }
  runLoader->UnloadRecPoints("all");
  runLoader->UnloadKinematics();
  runLoader->UnloadHeader();
  runLoader->UnloadgAlice();
  delete runLoader; 
  runLoader = 0;
  esdFile->Close();
  delete esdFile;
}

//_______________________________________________
void ProcessEvent(int iEv)
{
  runLoader->GetEvent(iEv);
  treeInp->GetEvent(iEv);
  //
  if (!tracker) InitTracker(esd->GetRunNumber());
  //
  AliStack *stack = runLoader->Stack(); 
  Int_t nParticles = stack ? stack->GetNtrack() : 0;
  Int_t nTracks= esd->GetNumberOfTracks();
  printf("Ntr: %d NPartMC: %d\n",nTracks,nParticles);
  const AliESDVertex* vtx = esd->GetPrimaryVertexSPD();
  if (!vtx || !vtx->GetStatus()) return;
  //
  TTree* treeITSRP = runLoader->GetTreeR("ITS",kFALSE);
  if (!treeITSRP) {
    printf("Failed to fetch ITS recpoints\n");
    exit(1);
  }
  //
#ifdef _TIMING_
  static Double_t cpuTime[AliITSSAPTracker::kNSW] = {0};
#endif
  printf("\n\n\nEvent: %d\n",iEv);
  TestTracker(treeITSRP, vtx);
  //
#ifdef _TIMING_
  if (vtx->GetStatus()==1) {
    double mlt = esd->GetMultiplicity()->GetNumberOfTracklets();
    for (int i=0;i<AliITSSAPTracker::kNSW;i++) {
      TStopwatch& sw = (TStopwatch&)tracker->GetStopwatch(i);
      double tm = sw.CpuTime();
      hTiming[i]->Fill(mlt, tm - cpuTime[i]);
      cpuTime[i] = tm;
    }
  }
#endif
  //
  //esd->Reset();
  //
}

//_________________________________________________
void TestTracker(TTree* tRP, const AliESDVertex* vtx)
{
  tracker->Clear(); 
  LoadClusters(tRP); 
  tracker->SetSPDVertex(vtx);
  tracker->SetBz(esd->GetMagneticField());
  tracker->ProcessEvent();
  //
  tracker->PrintTracklets();
  tracker->PrintTracks();  
  //
  //
  //  vtx->Print();
  //  esd->GetMultiplicity()->Print("t");
  //
  /*
  AliHeader* hd = runLoader->GetHeader();
  AliGenEventHeader* hdmc;
  if (tracker->GetTrackVertex().GetStatus()==1 && hd && (hdmc=hd->GenEventHeader()) ) {
    TArrayF vtxMC;
    hdmc->PrimaryVertex(vtxMC);
    printf("MCvtx: %f %f %f\n",vtxMC[0],vtxMC[1],vtxMC[2]);
  }
  */
}

//_________________________________________________
void LoadClusters(TTree* tRP)
{
  AliITSRecPointContainer* rpcont = AliITSRecPointContainer::Instance();
  TClonesArray* itsClusters = rpcont->FetchClusters(0,tRP);
  if(!rpcont->IsSPDActive()){
    printf("No SPD rec points found, multiplicity not calculated\n");
    tRP->Print();
    return;
  }
  int nMod = AliITSgeomTGeo::GetNModules();
  for (int imd=0;imd<nMod;imd++) {
    itsClusters = rpcont->UncheckedGetClusters(imd);
    int nClusters = itsClusters->GetEntriesFast();
    if (!nClusters) continue;
    while(nClusters--) {
      AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nClusters);
      if (!cluster) continue;
      tracker->AddCluster(cluster);
    }
  }
  //
}

//______________________________________________
Double_t* DefLogAx(double xMn,double xMx, int nbin)
{
  // get array for log axis
  if (xMn<=0 || xMx<=xMn || nbin<2) {
    printf("Wrong axis request: %f %f %d\n",xMn,xMx,nbin);
    return 0;
  }
  double dx = log(xMx/xMn)/nbin;
  double *xax = new Double_t[nbin+1];
  for (int i=0;i<=nbin;i++) xax[i]= xMn*exp(dx*i);
  return xax;
}


//______________________________________________
void InitTracker(int runNumber)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetRun(runNumber);
  tracker = new AliITSSAPTracker(); 
  tracker->Init();
  //
  //
#ifdef _TIMING_
  int nbMlt = 100;
  double *mltAx = DefLogAx(1,6000,nbMlt);
  for (int i=0;i<AliITSSAPTracker::kNSW;i++) {
    hTiming[i] = new TProfile(tracker->GetStopwatchName(i),tracker->GetStopwatchName(i),nbMlt,mltAx);
  }
#endif
  //  
}
