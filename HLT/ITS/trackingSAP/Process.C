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
#include "AliITStrackV2.h"
#include "AliESDUtils.h"
#include "AliVertexerTracks.h"
#include "AliCDBEntry.h"
#include "AliGRPRecoParam.h"
//

#include "HLT/ITS/trackingSAP/AliITSSAPAux.h"
#include "HLT/ITS/trackingSAP/AliITSSAPLayer.h"
#include "HLT/ITS/trackingSAP/AliITSSAPTracker.h"
#include "TProfile.h"
#include "TStopwatch.h"
#endif


AliITSSAPTracker* tracker=0;

TTree* treeInp = 0;
AliRunLoader* runLoader = 0;
AliESDEvent* esd=0;
AliESDEvent* esdF=0;
Int_t maxEvProc = -1, nEvProc = 0;

typedef struct {
  Bool_t  validSPD;
  Bool_t  validTrc;
  Bool_t  isPrim;
  Char_t  itsPatt;
  Char_t  ntlC;
  Char_t  ntlF;
  Char_t  mtlC;
  Char_t  mtlF;
  Char_t  ntrC;
  Char_t  ntrF;
  Int_t   pdg;
  Int_t   evID;
  Int_t   multSPDESD;
  Int_t   multSPD;
  Int_t   multTrc;
  Int_t   multSAESD;
  Int_t   multGLOESD;
  Int_t   trId;
  Float_t ptSA;
  Float_t ptGL;
  Float_t pt;
  Float_t ptMC;
  Float_t etaMC;
  Float_t phiMC;
  int     nvSPDHLTCont;
  int     nvTrcHLTCont;
  int     nvSPDESDCont;
  int     nvTrcESDCont;
  int     nvTrcHLTrefCont;
  int     nvTrcESDrefCont;
  Float_t vtxSPD[3];
  Float_t vtxTrcHLT[3];
  Float_t vtxTrcESD[3];
  Float_t vtxSPDESD[3];
  Float_t vtxMC[3];
  Float_t vtxTrcHLTref[3];
  Float_t vtxTrcESDref[3];
  //
  Float_t dcaTrc[2];
  Float_t dcaSA[2];
  Float_t dcaGL[2];
  //
  Float_t spdDPhi;
  Float_t spdDTht;
  Float_t spdChi2;
} tres_t;

tres_t tres;
TFile* flOut = 0;
TTree* trOut = 0;


void ProcessEvent(int iev);
void Process(const char* path, int maxEv=-1);
void ProcChunk(const char* path);
void TestTracker(TTree* tRP, const AliESDVertex* vtx);
void LoadClusters(TTree* tRP);
Double_t* DefLogAx(double xMn,double xMx, int nbin);
void CheckRecStatus();
void InitOutTree(const char* fname="rcInfo.root");
void SaveOutTree();
void InitTracker(int runNumber);

void CreateFakeEvent();

#ifdef _TIMING_
TProfile* hTiming[AliITSSAPTracker::kNSW];
#endif

void Process(const char* inpData, int maxEv)
{
  //
  maxEvProc = maxEv;
  InitOutTree();
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
      if (maxEvProc>0 && nEvProc>=maxEvProc) break;
    }
  }
#ifdef _CONTROLH_
  tracker->SaveHistos();
#endif
  SaveOutTree();
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
    nEvProc++;
    if (maxEvProc>0 && nEvProc>=maxEvProc) break;
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
  CheckRecStatus();
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
  //  tracker->PrintTracklets();
  //  tracker->PrintTracks();  
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
  //  AliITSRecPointContainer* rpcont = AliITSRecPointContainer::Instance();
  //  TClonesArray* itsClusters = rpcont->FetchClusters(0,tRP);
  //  if(!rpcont->IsSPDActive()){
  //    printf("No SPD rec points found, multiplicity not calculated\n");
  //    tRP->Print();
  //    return;
  //  }
  //  else printf("NSP0: %d\n",itsClusters->GetEntriesFast());
  static TClonesArray* itsModAdd[2198] = {0};
  static Bool_t first = kTRUE;
  if (first) {
    first = 0;
    for (int i=0;i<2198;i++) itsModAdd[i] = new TClonesArray("AliITSRecPoint");
  }
  int nMod = AliITSgeomTGeo::GetNModules();
  for (int imd=0;imd<nMod;imd++) {
    TClonesArray* itsClusters = itsModAdd[imd];
    tRP->SetBranchAddress("ITSRecPoints",&itsClusters);    
    tRP->GetEntry(imd);
    //    itsClusters = rpcont->UncheckedGetClusters(imd);
    
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
void CheckRecStatus()
{
  //
  static int nev = -1;
  static TBits patternMC;
  static vector<char> ntrCorr;
  static vector<char> ntrFake;
  static vector<char> mtlCorr;
  static vector<char> mtlFake;
  static vector<char> ntlCorr;
  static vector<char> ntlFake;
  static vector<float> spdDPhi;
  static vector<float> spdDTht;
  static vector<float> spdChi2;
  static vector<int> refMC2Tr;
  static vector<int> refMC2SA;
  static vector<int> refMC2GL;
  //
  AliStack* stack = 0;
  AliRunLoader* rl = AliRunLoader::Instance();
  if (!rl || !(stack=rl->Stack())) return;
  nev++;
  //
  CreateFakeEvent();
  //
  double bz = esd->GetMagneticField();
  AliGenEventHeader* mcGenH = rl->GetHeader()->GenEventHeader();
  TArrayF vtmc(3);
  mcGenH->PrimaryVertex(vtmc);
  for (int i=3;i--;) tres.vtxMC[i] = vtmc[i];
  //
  enum {kIsPrim=AliITSSAPTracker::kNLrActive,kValidTracklet,kValidTrack,kRecDone,kBitPerTrack};
  int nTrkMC = stack->GetNtrack();
  patternMC.SetBitNumber(nTrkMC*kBitPerTrack,0);
  patternMC.ResetAllBits();
  //
  ntrCorr.clear();
  ntrFake.clear();
  ntlCorr.clear();
  ntlFake.clear();
  mtlCorr.clear();
  mtlFake.clear();
  //
  spdDPhi.clear();
  spdDTht.clear();
  spdChi2.clear();
  //  
  ntrCorr.resize(nTrkMC,0);
  ntrFake.resize(nTrkMC,0);
  ntlCorr.resize(nTrkMC,0);
  ntlFake.resize(nTrkMC,0);
  mtlCorr.resize(nTrkMC,0);
  mtlFake.resize(nTrkMC,0);
  spdDPhi.resize(nTrkMC,0);
  spdDTht.resize(nTrkMC,0);
  spdChi2.resize(nTrkMC,0);
  refMC2Tr.resize(nTrkMC,0);
  refMC2SA.resize(nTrkMC,0);
  refMC2GL.resize(nTrkMC,0);
  //
  // fill MC track patterns
  for (int ilr=6;ilr--;) {
    AliITSSAPLayer *lr = tracker->GetLayer(ilr);
    int ncl = lr->GetNClusters();
    for (int icl=ncl;icl--;) {
      AliITSRecPoint* cl = lr->GetClusterUnSorted(icl);
      for (int j=0;j<3;j++) {
	int lb = cl->GetLabel(j);
	if (lb<0 || lb>=nTrkMC) break;
	patternMC.SetBitNumber(lb*kBitPerTrack+ilr,kTRUE);
      }
    }
  }
  //
  int nTrk = tracker->GetNTracklets();
  if (!nTrk) return;
  for (int itr=0;itr<nTrk;itr++) {
    const AliITSSAPTracker::SPDtracklet_t& trlet = tracker->GetTracklet(itr);
    //    PrintTracklet(itr);
    //
    int lbl = trlet.label;
    if (lbl==-3141593) continue;
    int lblA = TMath::Abs(lbl);
    if (lblA==lbl) ntlCorr[lblA]++;
    else           ntlFake[lblA]++;
    //
    if (spdChi2[lblA]==0 || spdChi2[lblA]<trlet.chi2) {
      spdChi2[lblA] = trlet.chi2;
      spdDPhi[lblA] = trlet.dphi;
      spdDTht[lblA] = trlet.dtht;
    }
  }
  //
  AliITSSAPLayer* lr0 = tracker->GetLayer(0);
  AliITSSAPLayer* lr1 = tracker->GetLayer(1);
  for (int itrm=0;itrm<nTrkMC;itrm++) {
    refMC2Tr[itrm] = -1;
    refMC2SA[itrm] = -1;
    refMC2GL[itrm] = -1;
    if (ntlCorr[itrm]+ntlFake[itrm]<2) continue;
    printf("\nExtra for tr %d nC:%d nF:%d\n",itrm,ntlCorr[itrm],ntlFake[itrm]);
    //
    int cnt = 0;
    for (int itr=0;itr<nTrk;itr++) {
      const AliITSSAPTracker::SPDtracklet_t& trlet = tracker->GetTracklet(itr);
      if (TMath::Abs(trlet.label)!=itrm) continue;
      AliITSSAPLayer::ClsInfo_t* clinf0 = lr0->GetClusterInfo(trlet.id1);
      AliITSSAPLayer::ClsInfo_t* clinf1 = lr1->GetClusterInfo(trlet.id2);
      printf("#%2d%s%4d chi:%.2f [%4d/%3d] [%4d/%3d]\n",cnt++,trlet.label<0 ? "-":"+",itr,trlet.chi2,
             trlet.id1,clinf0->detid,
             trlet.id2,clinf1->detid);
    }
  }
  //
  tres.multSAESD = 0;
  for (int it=esd->GetNumberOfTracks();it--;) {
    AliESDtrack* trc = esd->GetTrack(it);
    int lb  = trc->GetLabel();
    int lbA = TMath::Abs(lb);
    //
    if (trc->IsOn(AliESDtrack::kITSpureSA) && trc->IsOn(AliESDtrack::kITSrefit)) {
      tres.multSAESD++;
      if (lbA<nTrkMC) refMC2SA[lbA] = it;
    }
    else if (trc->IsOn(AliESDtrack::kTPCrefit) && trc->IsOn(AliESDtrack::kITSrefit))  {
      tres.multGLOESD++;
      if (lbA<nTrkMC) refMC2GL[lbA] = it;
    }
  }
  
  const AliMultiplicity* mltESD = esd->GetMultiplicity();
  nTrk = mltESD->GetNumberOfTracklets();
  tres.multSPDESD = nTrk;
  for (int itr=0;itr<nTrk;itr++) {
    int lb0 = mltESD->GetLabel(itr,0);
    int lb1 = mltESD->GetLabel(itr,1);    
    if (lb0==lb1) mtlCorr[lb1]++;
    else          mtlFake[lb1]++;
  }
  //
  nTrk = tracker->GetNTracks();
  for (int itr=0;itr<nTrk;itr++) {
    const AliITSSAPTracker::ITStrack_t &track = tracker->GetTrack(itr);
    //
    int lbl = track.label;
    if (lbl==-3141593) continue;
    int lblA = TMath::Abs(lbl);
    if (lblA==lbl) ntrCorr[lblA]++;
    else           ntrFake[lblA]++;
    refMC2Tr[lblA] = itr;
  }
  //
  // set reconstructability
  //
  const AliESDVertex* vtxSPD = tracker->GetSPDVertex();
  if (vtxSPD->GetStatus()>0) {
    tres.vtxSPD[0] = vtxSPD->GetX();
    tres.vtxSPD[1] = vtxSPD->GetY();
    tres.vtxSPD[2] = vtxSPD->GetZ();
  }
  else {
    tres.vtxSPD[0] = -999;
    tres.vtxSPD[1] = -999;
    tres.vtxSPD[2] = -999;
  }
  tres.nvSPDHLTCont = vtxSPD->GetNContributors();
  //
  const AliESDVertex* vtxTRC = &tracker->GetTrackVertex();
  if (vtxTRC->GetStatus()>0) {
    tres.vtxTrcHLT[0] = vtxTRC->GetX();
    tres.vtxTrcHLT[1] = vtxTRC->GetY();
    tres.vtxTrcHLT[2] = vtxTRC->GetZ();
  }
  else {
    tres.vtxTrcHLT[0] = -999;
    tres.vtxTrcHLT[1] = -999;
    tres.vtxTrcHLT[2] = -999;
  }
  tres.nvTrcHLTCont = vtxTRC->GetNContributors();
  //
  const AliESDVertex* vtxTRCESD = esd->GetPrimaryVertexTracks();
  if (vtxTRCESD->GetStatus()>0) {
    tres.vtxTrcESD[0] = vtxTRCESD->GetX();
    tres.vtxTrcESD[1] = vtxTRCESD->GetY();
    tres.vtxTrcESD[2] = vtxTRCESD->GetZ();
  }
  else {
    tres.vtxTrcESD[0] = -999;
    tres.vtxTrcESD[1] = -999;
    tres.vtxTrcESD[2] = -999;
  }
  tres.nvTrcESDCont = vtxTRCESD->GetNContributors();
  //
  const AliESDVertex* vtxTRCEXT = esdF->GetPrimaryVertexTracks();
  if (vtxTRCEXT->GetStatus()>0) {
    tres.vtxTrcHLTref[0] = vtxTRCEXT->GetX();
    tres.vtxTrcHLTref[1] = vtxTRCEXT->GetY();
    tres.vtxTrcHLTref[2] = vtxTRCEXT->GetZ();
  }
  else {
    tres.vtxTrcHLTref[0] = -999;
    tres.vtxTrcHLTref[1] = -999;
    tres.vtxTrcHLTref[2] = -999;
  }
  tres.nvTrcHLTrefCont = vtxTRCEXT->GetNContributors();
  //
  const AliESDVertex* vtxTrcESDref = esdF->GetPrimaryVertexSPD();
  if (vtxTrcESDref->GetStatus()>0) {
    tres.vtxTrcESDref[0] = vtxTrcESDref->GetX();
    tres.vtxTrcESDref[1] = vtxTrcESDref->GetY();
    tres.vtxTrcESDref[2] = vtxTrcESDref->GetZ();
  }
  else {
    tres.vtxTrcESDref[0] = -999;
    tres.vtxTrcESDref[1] = -999;
    tres.vtxTrcESDref[2] = -999;
  }
  tres.nvTrcESDrefCont = vtxTrcESDref->GetNContributors();
  //
  //
  const AliESDVertex* vtxSPDESD = esd->GetPrimaryVertexSPD();
  if (vtxSPDESD->GetStatus()>0) {
    tres.vtxSPDESD[0] = vtxSPDESD->GetX();
    tres.vtxSPDESD[1] = vtxSPDESD->GetY();
    tres.vtxSPDESD[2] = vtxSPDESD->GetZ();
  }
  else {
    tres.vtxSPDESD[0] = -999;
    tres.vtxSPDESD[1] = -999;
    tres.vtxSPDESD[2] = -999;
  }
  tres.nvSPDESDCont = vtxSPDESD->GetNContributors();
  //
  tres.trId = 0;
  tres.evID = nev;
  tres.multSPD = tracker->GetNTracklets();
  tres.multTrc = tracker->GetNTracks();
  //
  for (int itr=0;itr<nTrkMC;itr++) {
    tres.itsPatt = 0;
    int bitoffs = itr*kBitPerTrack;
    //
    if (patternMC.TestBitNumber(bitoffs+AliITSSAPTracker::kALrSPD1)) tres.itsPatt |= 0x1;
    if (patternMC.TestBitNumber(bitoffs+AliITSSAPTracker::kALrSPD2)) tres.itsPatt |= 0x2;
    tres.validSPD = tres.itsPatt == 0x3; // both spd layers
    //
    int nmiss = 0;
    for (int il=AliITSSAPTracker::kALrSDD1;il<=AliITSSAPTracker::kALrSSD2;il++) {
      Bool_t hit = patternMC.TestBitNumber(bitoffs+il);
      if (hit) tres.itsPatt |= 0x1<<il;
      else if (!tracker->GetSkipLayer(il)) nmiss++;
    }
    tres.validTrc = tres.validSPD && (nmiss<=tracker->GetMaxMissedLayers());
    //
    int trID = refMC2Tr[itr];
    int trSA = refMC2SA[itr];
    int trGL = refMC2GL[itr];
    //
    TParticle* mctr = stack->Particle(itr);
    tres.etaMC = mctr->Eta();
    if ( !tres.validSPD && !tres.validTrc && 
	 (ntlCorr[itr]+ntlFake[itr])==0 && 
	 (ntrCorr[itr]+ntrFake[itr])==0 &&
	 (mtlCorr[itr]+mtlFake[itr])==0 &&
	 trSA<0 && trGL<0  
			//	 && TMath::Abs(tres.etaMC)>1.3
   ) continue;
    //
    tres.pt = -1.0;
    tres.dcaTrc[0] = tres.dcaTrc[1] = -999;
    if (trID>=0) {
      //      printf("%d -> %d (of %d)\n",itr, trID,nTrk);
      const AliITSSAPTracker::ITStrack_t &track = tracker->GetTrack(trID);
      const AliExternalTrackParam& etp = track.paramInw;
      tres.pt = etp.Pt();
      etp.GetDZ(vtmc[0],vtmc[1],vtmc[2], bz, tres.dcaTrc);
    }
    tres.ptSA = -1.0;
    tres.dcaSA[0] = tres.dcaSA[1] = -999;
    if (trSA>=0) {
      AliESDtrack* trcESD = esd->GetTrack(trSA);
      tres.ptSA = trcESD->Pt();
      trcESD->GetDZ(vtmc[0],vtmc[1],vtmc[2], bz, tres.dcaSA);
    }
    tres.ptGL = -1.0;
    tres.dcaGL[0] = tres.dcaGL[1] = -999;
    if (trGL>=0) {
      AliESDtrack* trcESD = esd->GetTrack(trGL);
      tres.ptGL = trcESD->Pt();
      trcESD->GetDZ(vtmc[0],vtmc[1],vtmc[2], bz, tres.dcaGL);
    }

    tres.isPrim = stack->IsPhysicalPrimary(itr);
    tres.pdg = mctr->GetPdgCode();
    tres.ptMC =  mctr->Pt();
    tres.phiMC = mctr->Phi();
    //
    tres.ntlC = ntlCorr[itr];
    tres.ntlF = ntlFake[itr];
    tres.ntrC = ntrCorr[itr];
    tres.ntrF = ntrFake[itr];
    //
    tres.mtlC = mtlCorr[itr];
    tres.mtlF = mtlFake[itr];
    //
    tres.spdDPhi = spdDPhi[itr];
    tres.spdDTht = spdDTht[itr];
    tres.spdChi2 = spdChi2[itr];
    //
    trOut->Fill();
    tres.trId++;
  }
  //
}

//___________________________________________________
void InitOutTree(const char* fname)
{
  // output tree
  flOut = TFile::Open(fname,"recreate");
  trOut = new TTree("rcinfo","rcinfo");
  trOut->Branch("res",&tres);
}

//___________________________________________________
void SaveOutTree()
{
  if (!trOut) return;
  flOut->cd();
  trOut->Write();
  delete trOut;
  flOut->Close();
  delete flOut;
  trOut = 0;
  flOut = 0;
}


//______________________________________________
void InitTracker(int runNumber)
{
  AliCDBManager* man = AliCDBManager::Instance();
  man->SetRun(runNumber);
  tracker = new AliITSSAPTracker(); 
  tracker->Init();
  //  tracker->SetMaxMissedLayers(2);
  //  tracker->SetSkipLayer(2);
  //  tracker->SetSkipLayer(3);
  //  tracker->SetMaxMissedLayers(0);
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

void CreateFakeEvent()
{
  if (!esdF) {
    esdF = new AliESDEvent();
    esdF->CreateStdContent();
  }
  esdF->ResetStdContent();
  AliESDRun* rn = (AliESDRun*)esdF->GetESDRun();
  *rn = *(AliESDRun*)esd->GetESDRun();
  AliESDHeader* hd = (AliESDHeader*)esdF->GetHeader();
  *hd = *(AliESDHeader*)esd->GetHeader();
  static AliESDVertex vtFake;
  Double_t diamondxyzD[3]={esdF->GetDiamondX(),esdF->GetDiamondY(),0};    
  Double_t diamondcovxyD[6]={1,0,1,0,0,36};
  vtFake.SetXYZ(diamondxyzD);
  vtFake.SetCovarianceMatrix(diamondcovxyD);
  esdF->SetDiamond(&vtFake);
  //
  static AliESDtrack trESD;
  static AliITStrackV2 trITS;
  int nTrk = tracker->GetNTracklets(); 
  for (int itr = nTrk;itr--;) {
    const AliITSSAPTracker::ITStrack_t &track = tracker->GetTrack(itr);
    if (track.paramOut.TestBit(AliITSSAPTracker::kInvalidBit) ||
	track.paramInw.TestBit(AliITSSAPTracker::kInvalidBit)) continue;
    trITS.AliExternalTrackParam::operator=(track.paramInw);
    trITS.SetNumberOfClusters(track.ncl);
    trITS.SetChi2(track.chi2);
    trITS.SetLabel(track.label);
    trESD.UpdateTrackParams(&trITS,AliESDtrack::kITSrefit);
    trESD.SetFlag(AliESDtrack::kITSin | AliESDtrack::kITSout | AliESDtrack::kITSrefit);
    esdF->AddTrack(&trESD);
  }
  //
  static AliVertexerTracks* vtFinder = 0;
  if (!vtFinder) { // create vertexer
    vtFinder = new AliVertexerTracks(esdF->GetMagneticField());
    AliGRPRecoParam *grpRecoParam=0x0;
    AliCDBEntry* ent = AliCDBManager::Instance()->Get("GRP/Calib/RecoParam");
    grpRecoParam= (AliGRPRecoParam*) ((TObjArray*)ent->GetObject())->At(1);
    Int_t nCutsVertexer = grpRecoParam->GetVertexerTracksNCuts();
    Double_t *cutsVertexer = new Double_t[nCutsVertexer];
    grpRecoParam->GetVertexerTracksCutsITS(cutsVertexer,nCutsVertexer);
    vtFinder->SetITSMode();
    vtFinder->SetCuts(cutsVertexer,nCutsVertexer);
    vtFinder->SetConstraintOff();
    printf("Initialized vertexer for VertexTracks refit with field %f kG\n",esdF->GetMagneticField());
    //
  }
  /*
  const int kNCuts=21;
  const double vtCuts[kNCuts] = 
    {1.00e-01,1.00e-01,5.00e-01,3.00e+00,1.00e+00,3.00e+00,1.00e+02,
     1.00e+03,3.00e+00,3.00e+01,
     6.00e+00, // algo
     4.00e+00,7.00e+00,1.00e+03,
     5.00e+00,5.00e-02,1.00e-03,2.00e+00,1.00e+01,1.00e+00,5.00e+01};
  //
  */
  int ntrESD = esdF->GetNumberOfTracks();
  //  Bool_t res = AliESDUtils::RefitESDVertexTracks(esdF); //,kNCuts,vtCuts);
  AliESDVertex* vtxNCESD = vtFinder->FindPrimaryVertex(esd);
  if (vtxNCESD) {
    AliESDVertex* vtspd = (AliESDVertex*)esdF->GetPrimaryVertexSPD();
    *vtspd = *vtxNCESD;
    delete vtxNCESD;
  }
  //
  AliESDVertex* vtxNC = vtFinder->FindPrimaryVertex(esdF);
  if (vtxNC) {
     AliESDVertex* vttrc = (AliESDVertex*)esdF->GetPrimaryVertexTracks();
    *vttrc = *vtxNC;
    delete vtxNC;
  }
  
}

