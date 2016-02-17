/*
This macro loops over all the cluster (output of compClusHits.C) and creates a database containing the topologies of the cluster and infromation about:
pattern in TBits format, coordinate of the centre of the pixel contining the COG wrt down-left corner of the osculating rectangle (fraction of the pixel pithc),
distance between COG and centre of the pixel (fraction), number of fired pixels, number of colums, number of rows.
*/

#if !defined(__CInt_t__) || defined(__MAKECInt_t__)
#include "AliESDEvent.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TBits.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "../ITS/UPGRADE/AliITSUClusterPix.h"
#include "../ITS/UPGRADE/AliITSURecoLayer.h"
#include "../ITS/UPGRADE/AliITSURecoDet.h"
#include "../ITS/UPGRADE/AliITSUHit.h"
#include "../ITS/UPGRADE/AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliGeomManager.h"

#endif

TTree* treeInp = 0;
TTree* cluTree = 0;
AliRunLoader* runLoader = 0;
AliESDEvent* esd=0;

AliITSUGeomTGeo* gm=0;
Int_t nLayers = 0;
AliITSURecoDet *its = 0;
TObjArray patterns;
TArrayF pattCount;
Int_t nPatterns=0,nClusters=0;
TString outDBFileName = "clusterTopology.root";

void ProcessEvent(Int_t iev);
void ProcChunk(const char* path);
void BuildClTopDB(const char* inpData="AliESDs.root", const char* outDBFile="clusterTopology.root");
void AccountTopology(TBits& patt);
void CreateDB();

void BuildClTopDB(const char* inpData, const char* outDBFile)
{
  //
  // parse input list of single root file and process it
  outDBFileName = outDBFile;
  if (outDBFileName.IsNull()) outDBFileName = "clusterTopology.root";
  //
  AliGeomManager::LoadGeometry("geometry.root");
  gm = new AliITSUGeomTGeo(kTRUE);
  AliITSUClusterPix::SetGeom(gm);
  nLayers = gm->GetNLayers();
  its = new AliITSURecoDet(gm, "ITSInt_terface");
  its->CreateClusterArrays();
  esd = new AliESDEvent();
  //
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
  //
  CreateDB();
  //
}

void CreateDB()
{
  // sort patterns according to frequencies and store them
  printf("Storing %d patterns in %s\n",nPatterns, outDBFileName.Data());
  if (nPatterns<1) {
    printf("No patterns found\n");
    return;
  }
  TArrayI sortIndex;
  sortIndex.Set(nPatterns);
  TMath::Sort(nPatterns,pattCount.GetArray(),sortIndex.GetArray());
  //
  //create a file.txt containing topologies and corresponding frequecies
  ofstream top("topologies.txt");

  TObjArray arrStore(nPatterns);
  TVectorF arrFreq(nPatterns);
  TVectorF arrzCOGPix(nPatterns); //It's the cordinate of the pixel containing along the rows direction
  TVectorF arrxCOGPix(nPatterns); //It's the COG cordinate along the columns direction
  TVectorF arrzCOGshift(nPatterns); //It's the distance between COG and the centre of the pixel containing COG
  TVectorF arrxCOGshift(nPatterns);
  TVectorF arrNpix(nPatterns);
  TVectorF arrNcol(nPatterns);
  TVectorF arrNrow(nPatterns);

  TCanvas* cnv = new TCanvas("cnv","Freq",900,600);
  TH1F* hfreq = new TH1F("hfreq","Pattern frequency",300/2,-0.5, 300/2-0.5);
  hfreq->SetStats(0);
  hfreq->GetXaxis()->SetTitle("pattern ID");
  hfreq->SetFillColor(kRed);
  hfreq->SetFillStyle(3005);
  TLegend* leg = new TLegend(70, 300000, 140,450000 ,"","");
  leg->AddEntry((TObject*)0, Form("Number of patterns: %d",nPatterns), "");
  leg->AddEntry((TObject*)0, "Number of events: 50", "");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  
  for (Int_t i=0;i<nPatterns;i++) {
    Int_t ind = sortIndex[i];
    Float_t fr = arrFreq[i] = Float_t(pattCount[ind])/nClusters;
    arrStore.AddAt(patterns[ind],i);
    UInt_t idd = patterns[ind]->GetUniqueID();

    hfreq->Fill(i,pattCount[ind]);

    Int_t firedPixels=0;
    Int_t tempxCOG=0;
    Int_t tempzCOG=0;

    Int_t rc = idd>>16;//It's the rows span
    Int_t cl = idd&0xffff;//It's the columns span
    top << Form("\nPatt %4d R:%d C: %d Freq:%f \n\n",i,rc,cl,fr);
    for (Int_t ir=rc;ir--;) {
      top << "|"; 
      for (Int_t ic=0; ic<cl; ic++) {
        top << Form("%c",((TBits*)patterns[ind])->TestBitNumber(ir*cl+ic) ? '+':' ');
        if(((TBits*)patterns[ind])->TestBitNumber(ir*cl+ic)){
          firedPixels++;
          tempxCOG+=ir;
          tempzCOG+=ic;
        }
      }
      top << ("|\n");
    }
    Float_t xsh=Float_t((tempxCOG%firedPixels))/firedPixels; //distance between COG end centre of the pixel containing COG
    Float_t zsh=Float_t((tempzCOG%firedPixels))/firedPixels;
    tempxCOG/=firedPixels;
    tempzCOG/=firedPixels;
    if(xsh>0.5){
      tempxCOG+=1;
      xsh-=1;
    }
    if(zsh>0.5){
      tempzCOG+=1;
      zsh-=1;
    }
    Float_t xc=(Float_t) tempxCOG+0.5;
    Float_t zc=(Float_t) tempzCOG+0.5;
    arrxCOGPix[i]=xc;
    arrxCOGshift[i]=xsh;
    arrzCOGPix[i]=zc;
    arrzCOGshift[i]=zsh;
    arrNpix[i]=firedPixels;
    arrNrow[i]=rc;
    arrNcol[i]=cl;
    
    top << "\nxCentre: " <<  xc << " + " << xsh
    <<" zCentre: " << zc << " + " << zsh <<"\n\n...............................................\n";
  }

  cnv->cd();
  hfreq->Draw();
  leg->Draw();
  cnv->Print("Pattern_freq.jpg");

  top.close();

  printf("\nWe have %d clusters corresponding to %d patterns!!!\nEnjoy!!\n\n", nClusters, nPatterns);

  /*
  //Controlling if frequencies have been stored correctly
  ofstream a("frequency.txt");
  a << setw(15) << "i" << setw(20) << "Frequency" << endl;
  for(Int_t i=0; i<nPatterns; i++) a << setw(15) << i << setw(20) << arrFreq[i] << endl;
  a.close();
  */
  //
  TFile* flDB = TFile::Open(outDBFileName.Data(), "recreate");
  flDB->WriteObject(&arrStore,"TopDB","kSingleKey");
  flDB->WriteObject(&arrFreq, "TopFreq","kSingleKey");
  flDB->WriteObject(&arrxCOGPix,"xCOG","kSingleKey");
  flDB->WriteObject(&arrxCOGshift,"xShift","kSingleKey");
  flDB->WriteObject(&arrzCOGPix,"zCOG","kSingleKey");
  flDB->WriteObject(&arrzCOGshift,"zShift","kSingleKey");
  flDB->WriteObject(&arrNpix,"NPix","kSingleKey");
  flDB->WriteObject(&arrNrow,"NRow","kSingleKey");
  flDB->WriteObject(&arrNcol,"NCol","kSingleKey");
  //
  flDB->Close();
  delete flDB;
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
  //
  runLoader = AliRunLoader::Open(Form("%sgalice.root",dir.Data()));
  if (!runLoader) {
    printf("galice not found\n");
    return;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  //  runLoader->LoadHeader();
  //  runLoader->LoadKinematics();
  runLoader->LoadRecPoints("ITS");
  // ESD
  TFile* esdFile = TFile::Open(Form("%sAliESDs.root",dir.Data()));
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file\n");
    //    runLoader->UnloadKinematics();
    //    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }

  treeInp = (TTree*) esdFile->Get("esdTree");
  if (!treeInp) {
    printf("Error: no ESD tree found\n");
    //    runLoader->UnloadKinematics();
    //    runLoader->UnloadHeader();
    runLoader->UnloadgAlice();
    delete runLoader;
    return;
  }
  esd->ReadFromTree(treeInp);
  //
  for (Int_t iEv=0;iEv<runLoader->GetNumberOfEvents(); iEv++) {
    ProcessEvent(iEv);
  }
  runLoader->UnloadRecPoints("all");
  //  runLoader->UnloadKinematics();
  //  runLoader->UnloadHeader();
  runLoader->UnloadgAlice();
  delete runLoader; 
  runLoader = 0;
  esdFile->Close();
  delete esdFile;
}


void ProcessEvent(Int_t iEv) 
{
  // process single event clusters
  runLoader->GetEvent(iEv);
  treeInp->GetEvent(iEv);
  cluTree = runLoader->GetTreeR("ITS",kFALSE);
  if (!cluTree) {
    printf("Did not get cluster tree for event %d\n",iEv);
    return;
  }
  //
  for (Int_t ilr=0;ilr<nLayers;ilr++) {
    TBranch* br = cluTree->GetBranch(Form("ITSRecPoints%d",ilr));
    if (!br) {printf("Did not find cluster branch for lr %d\n",ilr); exit(1);}
    br->SetAddress(its->GetLayerActive(ilr)->GetClustersAddress());
  }
  cluTree->GetEntry(0); 
  its->ProcessClusters();
  //
  TBits currTop;
  //
  for (Int_t ilr=0;ilr<nLayers;ilr++) {
    AliITSURecoLayer* lr = its->GetLayerActive(ilr);
    TClonesArray* clr = lr->GetClusters();
    Int_t nClu = clr->GetEntries();
    printf("Ev%d Lr:%d Ncl:%d\n",iEv,ilr,nClu);
    //    printf("Layer %d : %d clusters\n",ilr,nClu);
    //
    for (Int_t icl=0;icl<nClu;icl++) {
      AliITSUClusterPix *cl = (AliITSUClusterPix*)clr->At(icl);
      currTop.Clear();
      UShort_t rs = (UShort_t)cl->GetPatternRowSpan();
      UShort_t cs = (UShort_t)cl->GetPatternColSpan();
      for (Int_t ir=0;ir<rs;ir++)
	     for (Int_t ic=0;ic<cs;ic++) 
	       if (cl->TestPixel(ir,ic)) currTop.SetBitNumber(ir*cs+ic);
      //
      currTop.SetUniqueID( (rs<<16) + (cs));
      //
      AccountTopology(currTop);
    }
  }
  //
}

//___________________________________________
void AccountTopology(TBits& patt)
{
  // account topology in the database
  Bool_t newPatt = kTRUE;
  nClusters++;
  for (Int_t ip=0;ip<nPatterns;ip++) {
    TBits* pattOld = (TBits*)patterns.At(ip);
    if (*pattOld == patt && pattOld->GetUniqueID()==patt.GetUniqueID()) {
      // require that the matrix size and the bit patterns are equal
      newPatt = kFALSE;
      pattCount[ip]++;
      break;
    }
  }
  //
  if (newPatt) {
    TBits* pt = new TBits(patt);
    patterns.AddLast(pt);
    if (pattCount.GetSize()<nPatterns+100) pattCount.Set(nPatterns+100);
    pattCount[nPatterns++] = 1;
  }
  //
}