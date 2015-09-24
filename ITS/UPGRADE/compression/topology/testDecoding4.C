#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TBits.h"
#include "AliITSUClusterPix.h"
#include "AliITSURecoLayer.h"
#include "AliITSURecoDet.h"
#include "AliITSUHit.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliGeomManager.h"
#include "AliStack.h"
#include "AliLoader.h"
#include "AliCDBManager.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGeoMatrix.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "./MyClasses/Topology.h"
#include "./MyClasses/TopDatabase.h"

#endif

TopDatabase DB;

TStopwatch Timer;
Int_t nClashes=0;

TObjArray histoArr;
enum {kNPixAll=0,kNPixSPL=1,kDR=0,kDTXodd,kDTXeven,kDTZ, kDTXoddSPL,kDTXevenSPL,kDTZSPL};

void LoadDB(const char* fname);
void Top2Word(const TBits* top, UChar_t* Word);
UInt_t FuncMurmurHash2(const void * key, Int_t len);
std::ostream& printTop(TBits top, std::ostream &out);
void testDecoding4(int nev=-1, int nRepetintions=100);
Int_t Top2Int(const TBits* top);

TBits clTop;

void testDecoding4(int nev, int nRepetintions){
  const int kSplit=0x1<<22;
  const int kSplCheck=0x1<<23;
  //
  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeSim");
  gSystem->Load("libITSUpgradeRec");
  gROOT->SetStyle("Plain");
  gStyle->SetStripDecimals("kFalse");

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Align/Data",
			  Form("local://%s",gSystem->pwd()));
  man->SetSpecificStorage("ITS/Calib/RecoParam",
			  Form("local://%s",gSystem->pwd()));
  man->SetRun(0);

  gAlice=NULL;
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  runLoader->LoadgAlice();

  gAlice = runLoader->GetAliRun();

  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadRecPoints();
  runLoader->LoadSDigits();
  runLoader->LoadHits();

  AliLoader *dl = runLoader->GetDetectorLoader("ITS");

  AliGeomManager::LoadGeometry("geometry.root");
  TObjArray algITS;
  AliGeomManager::LoadAlignObjsFromCDBSingleDet("ITS",algITS);
  AliGeomManager::ApplyAlignObjsToGeom(algITS);
  //
  AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE);
  AliITSUClusterPix::SetGeom(gm);
  //
  AliITSURecoDet *its = new AliITSURecoDet(gm, "ITSinterface");
  its->CreateClusterArrays();
  //
  Double_t xg1,yg1,zg1=0.,xg0,yg0,zg0=0.,tg0;
  Double_t xExit,yExit,zExit,xEnt,yEnt,zEnt,tof1;
  //
  TTree * cluTree = 0x0;
  TTree *hitTree = 0x0;
  TClonesArray *hitList=new TClonesArray("AliITSUHit");
  //
  TObjArray arrMCTracks; // array of hit arrays for each particle
  //
  Float_t xyzClGloF[3];
  Double_t xyzClGlo[3],xyzClTr[3];
  Int_t labels[3];
  int nLab = 0;
  int nlr=its->GetNLayersActive();
  int ntotev = (Int_t)runLoader->GetNumberOfEvents();
  printf("N Events : %i \n",ntotev);
  if (nev>0) ntotev = TMath::Min(nev,ntotev);
  //
  // load topology database
  LoadDB("TopologyDatabase.root");
  DB.PrintDB("check.txt");
  Int_t nPatterns = DB.GetN();

  TCanvas* c = new TCanvas("c","cInterpolation");
  c->Divide(2,1);
  
  TH1F* RealInter = new TH1F("RealInter","Real time with interpolation search",50,4e-7,8e-7);
  RealInter->SetDirectory(0);
  RealInter->GetXaxis()->SetTitle("t (s)");
  RealInter->SetFillColor(kBlue);
  RealInter->SetFillStyle(3005);
  RealInter->SetNdivisions(505,"X");

  TH1F* CPUInter = new TH1F("CPUInter","CPU time with interpolation search",50,1e-7,8e-7);
  CPUInter->SetDirectory(0);
  CPUInter->GetXaxis()->SetTitle("t (s)");
  CPUInter->SetFillColor(kRed);
  CPUInter->SetFillStyle(3005);


  for(Int_t irep=0; irep<nRepetintions; irep++){
    printf("%d / %d\n", irep, nRepetintions);  
    Int_t totClusters=0;
    for (Int_t iEvent = 0; iEvent < ntotev; iEvent++) {
      //printf("\n Event %i \n",iEvent);
      runLoader->GetEvent(iEvent);
      AliStack *stack = runLoader->Stack();
      cluTree=dl->TreeR();
      hitTree=dl->TreeH();
      hitTree->SetBranchAddress("ITS",&hitList);
      // 
      // read clusters
      for (int ilr=nlr;ilr--;) {
	TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",ilr));
	if (!br) {printf("Did not find cluster branch for lr %d\n",ilr); exit(1);}
	br->SetAddress(its->GetLayerActive(ilr)->GetClustersAddress());
      }
      cluTree->GetEntry(0); 
      its->ProcessClusters();
      //
      //printf(" tree entries: %lld\n",cluTree->GetEntries());
      //
      for (int ilr=0;ilr<nlr;ilr++) {
	AliITSURecoLayer* lr = its->GetLayerActive(ilr);
	TClonesArray* clr = lr->GetClusters();
	int nClu = clr->GetEntries();
	//printf("Layer %d : %d clusters\n",ilr,nClu);
	//
	for (int icl=0;icl<nClu;icl++){
	  //if(icl%1000==0)printf("%d / %d\n", icl, nClu);
	  AliITSUClusterPix *cl = (AliITSUClusterPix*)clr->At(icl);
	  Timer.Start(!totClusters);
	  Int_t num = DB.FromCluster2GroupID(*cl);
	  Timer.Stop();
	  totClusters++;
	}
      }
      //layerClus.Clear();
      //
      arrMCTracks.Delete();
      Double_t CPUTime = (Timer.CpuTime())/totClusters;
      Double_t RealTime = (Timer.RealTime())/totClusters;
      CPUInter->Fill(CPUTime);
      RealInter->Fill(RealTime);    
    }//event loop
  }
  c->cd(1);
  RealInter->Draw();
  c->cd(2);
  CPUInter->Draw();
  c->Print("DecodingTime.pdf");
}
 
void LoadDB(const char* fname)
{
  // load database
  TFile* fl = TFile::Open(fname);
  DB = *((TopDatabase*)fl->Get("DB"));
}

void Top2Word(const TBits* top, UChar_t* Word){
  UInt_t UID = top->GetUniqueID();
  Int_t rs = UID>>16;
  Int_t cs = UID&0xffff;
  Word[0]=rs;
  Word[1]=cs;
  UChar_t tempChar=0;
  Int_t index=2;
  Int_t BitCounter=7;
  for(Int_t ir=0; ir<rs; ir++){
    for(Int_t ic=0; ic<cs; ic++){
      if(BitCounter<0) {
	Word[index]=tempChar;
	tempChar=0;
	BitCounter=7;
	index++;
      }	
      if(top->TestBitNumber(ir*cs+ic)) tempChar+=(1<<BitCounter);
      BitCounter--;
    }
  }
  Word[index]=tempChar;
}

UInt_t FuncMurmurHash2(const void* key, Int_t len){
  // 'm' and 'r' are mixing constants generated offline.
  const UInt_t m =0x5bd1e995;
  const Int_t r = 24;
  // Initialize the hash
  UInt_t h = 0;
  // Mix 4 bytes at a time into the hash
  const UChar_t* data = (const UChar_t *)key;
  /*________________DEBUG____________________
    TBits verifica;
    Word2Top(data,verifica);
    printf("Topology in the hash scope:\n\n");
    printTop(verifica,cout);
  *///_________________________________________
  //Int_t recIndex=0;
  while(len >= 4){
    UInt_t k = *(UInt_t*)data;
    //if(recIndex==0) printf("first vale of k: %u\n", k);
    k *= m;
    k ^= k >> r;
    k *= m;
    h *= m;
    h ^= k;
    data += 4;
    len -= 4;
    //recIndex++;
  }
  // Handle the last few bytes of the input array
  switch(len){
  case 3: h ^= data[2] << 16;
  case 2: h ^= data[1] << 8;
  case 1: h ^= data[0];
    h *= m;
  };
  // Do a few final mixes of the hash to ensure the last few
  // bytes are well-incorporated.
  h ^= h >> 13;
  h *= m;
  h ^= h >> 15;
  return h;
}

std::ostream& printTop(TBits top, std::ostream &out){
  UInt_t UID = top.GetUniqueID();
  Int_t rs = UID>>16;
  Int_t cs = UID&0xffff;
  for (Int_t ir=0;ir<rs;ir++){
    out << "|"; 
    for (Int_t ic=0; ic<cs; ic++) {
      out << Form("%c",top.TestBitNumber(ir*cs+ic) ? '+':' ');
    }
    out << ("|\n");
  }
  out<< endl;
}

Int_t Top2Int(const TBits* top){
  Int_t output=0;
  for(Int_t i=0; i<32; i++){
    if(top->TestBitNumber(i)) output+=(1<<(31-i));
  }
  return output;
}
