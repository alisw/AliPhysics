// $Id$

#ifndef __CINT__
#include <stdio.h>
#include <Riostream.h>
#include <time.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TFile.h>
#include <TMath.h>
#include <TChain.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <AliStack.h>
#include <AliRunLoader.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliJetParticle.h>
#include <AliJetParticlesReader.h>
#include <AliJetParticlesReaderKine.h>
#include <AliJetParticlesReaderKineGoodTPC.h>
#include <AliJetParticlesReaderESD.h>
#include <AliJetParticlesReaderHLT.h>
#include <AliJetEventParticles.h>
#endif

#if 0
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kPHOSpid=0x10000, kHMPIDpid=0x20000, kEMCALpid=0x40000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
#endif

//#define constrained
#define useHoughITS

void createEvents(Char_t *dir, Int_t files,Char_t *sdir=".");
void createEvents(Int_t esd, TObjArray *dirs, Float_t ptcut, Char_t *savefile);
void createEvents(Int_t esd,Char_t *dir="./input",Float_t ptcut=0.05,Char_t *savefile=0);
void createEvents(Int_t esd,Char_t *dirfile,Char_t *dir,Float_t ptcut=0.05,Char_t *savefile=0);
void createMixedEvents(Char_t *fname,Char_t *hname,Char_t *savefile=0,Int_t nMaxEvents=-1);
void createDirs(Char_t *dir,Int_t files,Char_t *savefile);
void testEvents(Char_t *sdir);
void readEvents(Char_t *filename,Int_t nMaxEvents=-1);
void displayMixedEvents(Char_t *fname,Char_t *hname,Float_t ptcut=2.0,Int_t nMaxEvents=-1);

//-----------------------

void createEvents(Char_t *dir, Int_t files, Char_t *sdir)
{
  createDirs(dir,files,"./dirlist.txt");
  Char_t filename[1000];
  sprintf(filename,"%s/aliev-type0.root",sdir);
  createEvents(0,"./dirlist.txt",dir,0.4,filename);
  sprintf(filename,"%s/aliev-type1.root",sdir);
  createEvents(1,"./dirlist.txt",dir,0.4,filename);
  sprintf(filename,"%s/aliev-type2.root",sdir);
  createEvents(2,"./dirlist.txt",dir,0.4,filename);
  sprintf(filename,"%s/aliev-type3.root",sdir);
  createEvents(3,"./dirlist.txt",dir,0.4,filename);
#ifdef useHoughITS
  sprintf(filename,"%s/aliev-type4.root",sdir);
  createEvents(4,"./dirlist.txt",dir,0.4,filename);
#endif
  sprintf(filename,"%s/aliev-type10.root",sdir);
  createEvents(10,"./dirlist.txt",dir,0.4,filename);
  sprintf(filename,"%s/aliev-type11.root",sdir);
  createEvents(11,"./dirlist.txt",dir,0.4,filename);
  sprintf(filename,"%s/aliev-type12.root",sdir);
  createEvents(12,"./dirlist.txt",dir,0.4,filename);
}

void createEvents(Int_t esd,Char_t *dirfile,Char_t *dir,
                  Float_t ptcut,Char_t *savefile)
{
  TObjArray *dirs=new TObjArray(10000);

  FILE *c=fopen(dirfile,"r");
  if(!c) return;

  Int_t saveErrIgLevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kFatal;
  while(!feof(c)){
    Char_t adddir[1024];
    fscanf(c,"%s\n",adddir);
    dirs->Add(new TObjString(adddir));
    //cout << adddir << endl;
  }
  gErrorIgnoreLevel=saveErrIgLevel;

  Char_t buffer[8096];
  if(savefile)
    sprintf(buffer,"%s",savefile);
  else
    sprintf(buffer,"%s/aliev-type%d.root",dir,esd);

  createEvents(esd,dirs,ptcut,buffer);
}

void createEvents(Int_t esd,Char_t *dir,
		  Float_t ptcut,Char_t *savefile)
{
  TObjArray *dirs=new TObjArray(1);
  dirs->Add(new TObjString(dir));

  Char_t buffer[8096];
  if(savefile)
    sprintf(buffer,"%s",savefile);
  else
    sprintf(buffer,"%s/aliev-type%d.root",dir,esd);

  createEvents(esd,dirs,ptcut,buffer);
}

void createEvents(Int_t esd, TObjArray *dirs, Float_t ptcut, Char_t *savefile)
{
  AliJetParticlesReader *reader=0;
  TString desc="";
  if (esd==0) {
    reader=new AliJetParticlesReaderKine(dirs);
    //((AliJetParticlesReaderKine*)reader)->SetUseTracks(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetCharged(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetEM(kFALSE);
    ((AliJetParticlesReaderKine*)reader)->SetNeutral(kFALSE);
    desc+="Kine";
  } else if(esd==1){
#ifdef constrained
    reader=new AliJetParticlesReaderESD(1,dirs);
    ((AliJetParticlesReaderESD*)reader)->SetCompareFlag(
	                   AliESDtrack::kITSrefit+AliESDtrack::kTPCrefit+AliESDtrack::kTRDrefit);
#else
    reader=new AliJetParticlesReaderESD(0,dirs);
    ((AliJetParticlesReaderESD*)reader)->SetCompareFlag(AliESDtrack::kTPCin);
#endif
    desc+="ESD";
  } else if (esd==2) {
    reader=new AliJetParticlesReaderHLT(kTRUE,dirs);
    //((AliJetParticlesReaderHLT*)reader)->SetMinHits(30);
    desc+="HLTConf";
  } else if (esd==3) {
    reader=new AliJetParticlesReaderHLT(kFALSE,dirs);
    desc+="HLTHough";
  } else if (esd==4) {
    reader=new AliJetParticlesReaderHLT(kFALSE,dirs,"AliESDits.root");
    //((AliJetParticlesReaderHLT*)reader)->SetMinHits(4);
    desc+="HLTITSHough";
  } else if (esd==10) {
    reader=new AliJetParticlesReaderKine(dirs);
    //((AliJetParticlesReaderKine*)reader)->SetUseTracks(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetCharged(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetEM(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetNeutral(kTRUE);
    desc+="Kine-All";
  } else if (esd==11) {
    reader=new AliJetParticlesReaderKine(dirs);
    //((AliJetParticlesReaderKine*)reader)->SetUseTracks(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetCharged(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetEM(kTRUE);
    ((AliJetParticlesReaderKine*)reader)->SetNeutral(kFALSE);
    desc+="Kine-EM";
  } else if (esd==12) {
    reader=new AliJetParticlesReaderKineGoodTPC(dirs);
    desc+="GoodTPC";
  } else {
    cout << "Parameter settings: " << endl;
    cout << "0  == Kine (ch)"  << endl;
    cout << "1  == TPC"  << endl;
    cout << "2  == Conformal"  << endl;
    cout << "3  == Hough"  << endl;
    cout << "4  == ITSHough"  << endl;
    cout << "10 == Kine (all)"  << endl;
    cout << "11 == Kine (ch+em)"  << endl;
    cout << "12 == Kine (good tpc)"  << endl;
    return;
  }
  //reader->ReadEventsFromTo(1,0);
  reader->SetPtCut(ptcut,100);
  reader->SetEtaCut(-0.9,0.9);

  desc+=" Ptcut: ";
  desc+=ptcut;

  TFile* file = new TFile(savefile,"RECREATE");
  TTree* tree = new TTree("AJEPtree","AliJetEventParticles Tree");
  reader->SetTree(tree);
  Int_t counter=0;
  while(reader->Next())
    {
      AliJetEventParticles *ev=new AliJetEventParticles(*reader->GetEventParticles());
      if(gErrorIgnoreLevel<=kWarning){
	cout << "Read event: " << ++counter << endl;
	ev->Print();
      }
      delete ev;
    }
  file->Write();
  file->Close();
  delete file;
  delete reader;
}

void createDirs(Char_t *dir,Int_t files,Char_t *savefile)
{
  Char_t buffer[8096];
  if(savefile)
    sprintf(buffer,"%s",savefile);
  else
    sprintf(buffer,"%s/dirinput.txt",dir);

  FILE *c=fopen(buffer,"w");
  if(!c) return;

  Int_t saveErrIgLevel=gErrorIgnoreLevel;
  gErrorIgnoreLevel=kFatal;
  for(Int_t i=0;i<files;i++){
    Char_t adddir[1024];
    //sprintf(adddir,"%s/%05d",dir,i);
    sprintf(adddir,"%s/%d",dir,i);
    Char_t fname[1024];
    sprintf(fname,"%s/galice.root",adddir);
    TFile f(fname);
    if(!f.IsOpen()) continue;
    //cout << f.GetErrno() << endl;
    f.Close();
    AliRunLoader *r = AliRunLoader::Open(fname); 
    if(!r) continue;
    if(r->GetNumberOfEvents() <= 0){
      delete r;
      continue;
    }
    if(r->LoadKinematics()){
      delete r;
      continue;
    }
    delete r;
    sprintf(fname,"%s/AliESDs.root",adddir);
    TFile g(fname);
    if(!g.IsOpen()) continue;

    TTree *tree = dynamic_cast<TTree*>(g.Get("esdTree"));
    if(!tree || tree->GetEntries()<=0){
      g.Close();
      continue;
    }
    g.Close();
#ifdef useHoughITS
    sprintf(fname,"%s/AliESDits.root",adddir);
    TFile h(fname);
    if(!h.IsOpen()) continue;

    tree = dynamic_cast<TTree*>(h.Get("esdTree"));
    if(!tree || tree->GetEntries()<=0){
      h.Close();
      continue;
    }
    h.Close();
#endif
   
    cout << "Added " << adddir << endl;
    fprintf(c,"%s\n",adddir);
  }
  gErrorIgnoreLevel=saveErrIgLevel;
  fclose(c);
}

void createMixedEvents(Char_t *fname,Char_t *hname,Char_t *savefile,Int_t nMaxEvents)
{
  //connect to events
  TChain *theTree = new TChain("AJEPtree");
  theTree->Add(fname);
  AliJetEventParticles *ev=new AliJetEventParticles();
  theTree->SetBranchAddress("particles",&ev);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;
  //cout << "Found " << nMaxEvents << " in " << fname << endl;

  TChain *backtheTree = new TChain("AJEPtree");
  backtheTree->Add(hname);
  AliJetEventParticles *backev=new AliJetEventParticles();
  backtheTree->SetBranchAddress("particles",&backev);

  Int_t backtreeentries=(Int_t)backtheTree->GetEntries();
  Int_t nPerBackground=nMaxEvents/backtreeentries;
  if(nPerBackground==0) nPerBackground=1;

  Char_t buffer[8096];
  if(savefile)
    sprintf(buffer,"%s",savefile);
  else
    sprintf(buffer,"./aliev-mixed.root");
  TFile* file = new TFile(buffer,"RECREATE");
  TTree* tree = new TTree("AJEPtree","AliJetEventParticles Tree");
  AliJetEventParticles *mixedev=new AliJetEventParticles(0);
  tree->Branch("particles","AliJetEventParticles",&mixedev,32000,1);

  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent = 0;
  Int_t nEventHijing = -1;
  Int_t nEventHijingCounter = nPerBackground;
  while(nEvent<nMaxEvents){

    //get the event
    theTree->GetEvent(nEvent);

    //load background if needed
    if(nEventHijingCounter==nPerBackground){    
      backev->Reset();
      nEventHijing++;
      backtheTree->GetEvent(nEventHijing);
      if(nEventHijing==backtreeentries) nEventHijing=0;
      nEventHijingCounter=0;
    }

    backev->AddSignal(*ev);
    TString dummy="Counter: ";
    dummy+=nEvent;
    dummy+=" ";
    //dummy+=ev->GetHeader();
    dummy+="(Pythia ";dummy+=nEvent;
    dummy+=" ";
    //dummy+=backev->GetHeader();
    dummy+=", Hijing ";dummy+=nEventHijing;
    dummy+=")";
    mixedev->Set(*backev);
    mixedev->SetHeader(dummy);
    if(gErrorIgnoreLevel<=kWarning){
      cout << "Read event: " << nEvent << " " << nEventHijing << endl;
      mixedev->Print();
    }

    tree->Fill();

    nEvent++;
    nEventHijingCounter++;
    ev->Reset();
    mixedev->Reset();
  } //end of nev loop

  file->Write();
  file->Close();
  delete file;
  delete ev;
  delete theTree;
  delete backtheTree;
  delete backev;
}

void testEvents(Char_t *sdir)
{
  Char_t filename[1000];
  sprintf(filename,"%s/aliev-type0.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type1.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type2.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type3.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type10.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type11.root",sdir);
  readEvents(filename,0);
  sprintf(filename,"%s/aliev-type12.root",sdir);
}

void readEvents(Char_t *filename,Int_t nMaxEvents)
{
  //connect to events
  TChain *theTree = new TChain("AJEPtree");
  theTree->Add(filename);
  AliJetEventParticles *ev=new AliJetEventParticles();
  theTree->SetBranchAddress("particles",&ev);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;

  cout << "Found " << treeentries << ", but use " << nMaxEvents << " in " << filename << endl;

  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent = 0;
  while(nEvent<nMaxEvents){
    if ((nEvent % 100) == 0) {
      cout << "Reading event " << nEvent << endl;
    }

    //connect the event
    theTree->GetEvent(nEvent);

    ev->Print();
#if 0
    const TClonesArray *p=ev->GetParticles();
    for(Int_t i=0;i<p->GetEntriesFast();i++) {
      cout << i << endl;
      ((AliJetParticle*)p->At(i))->Print();
    }
#endif
    ev->Clear();
    nEvent++;
  } //end of nev loop

  delete ev;
  delete theTree;
}

void displayMixedEvents(Char_t *fname,Char_t *hname,Float_t ptcut,Int_t nMaxEvents)
{
  //connect to events
  TChain *theTree = new TChain("AJEPtree");
  theTree->Add(fname);
  AliJetEventParticles *ev=new AliJetEventParticles();
  theTree->SetBranchAddress("particles",&ev);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;
  //cout << "Found " << nMaxEvents << " in " << fname << endl;

  TChain *backtheTree = new TChain("AJEPtree");
  backtheTree->Add(hname);
  AliJetEventParticles *backev=new AliJetEventParticles();
  backtheTree->SetBranchAddress("particles",&backev);

  Int_t backtreeentries=(Int_t)backtheTree->GetEntries();
  Int_t nPerBackground=nMaxEvents/backtreeentries;
  if(nPerBackground==0) nPerBackground=1;

  AliJetEventParticles *mixedev=new AliJetEventParticles(0);

  Int_t neta=20;
  Int_t nphi=62;

  TCanvas *c1=new TCanvas("c1","",1000,500);
  c1->Divide(2,1);
  c1->Update();
  gStyle->SetLabelSize(0.03,"XYZ");
  gStyle->SetTitleOffset(1.0,"XYZ");
  gStyle->SetTitleOffset(1.3,"Y");
  gROOT->ForceStyle();

  TH2F *h2pp=new TH2F("h2pp","E_{T}^{ch} in PP [GeV];#phi;#eta",nphi,0,TMath::TwoPi(),neta,-1,1);
  h2pp->SetStats(0);
  TH2F *h2pb=new TH2F("h2pb","E_{T}^{ch} in PbPb [GeV];#phi;#eta;",nphi,0,TMath::TwoPi(),neta,-1,1);
  h2pb->SetStats(0);

  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent = 0;
  Int_t nEventHijing = -1;
  Int_t nEventHijingCounter = nPerBackground;
  while(nEvent<nMaxEvents){

    //get the event
    theTree->GetEvent(nEvent);

    //load background if needed
    if(nEventHijingCounter==nPerBackground){    
      backev->Reset();
      nEventHijing++;
      backtheTree->GetEvent(nEventHijing);
      if(nEventHijing==backtreeentries) nEventHijing=0;
      nEventHijingCounter=0;
    }
    backev->AddSignal(*ev);
    TString dummy="Counter: ";
    dummy+=nEvent;
    dummy+=" ";
    //dummy+=ev->GetHeader();
    dummy+="(Pythia ";dummy+=nEvent;
    dummy+=" ";
    //dummy+=backev->GetHeader();
    dummy+=", Hijing ";dummy+=nEventHijing;
    dummy+=")";
    mixedev->Set(*backev);
    mixedev->SetHeader(dummy);
    if(gErrorIgnoreLevel<=kWarning){
      cout << "Read event: " << nEvent << " " << nEventHijing << endl;
      mixedev->Print();
    }

    // loop over all particles...
    c1->cd(1);
    h2pp->Reset();
    AliJetParticle *aliparticle = NULL;
    TIterator *iter = ev->GetParticles()->MakeIterator();
    while ((aliparticle = (AliJetParticle *) iter->Next()) != NULL) {
      Float_t pt=aliparticle->Pt();
      if(pt<ptcut) continue;
      h2pp->Fill(aliparticle->Phi(),aliparticle->Eta(),pt);
    }
    delete iter;
    h2pp->Draw("Lego");
    c1->cd(2);
    h2pb->Reset();
    iter = mixedev->GetParticles()->MakeIterator();
    while ((aliparticle = (AliJetParticle *) iter->Next()) != NULL) {
      Float_t pt=aliparticle->Pt();
      if(pt<ptcut) continue;
      h2pb->Fill(aliparticle->Phi(),aliparticle->Eta(),pt);
    }
    delete iter;

    h2pb->Draw("Lego");
    c1->Update();
    Char_t c;cin >> c; 
    if(c=='s'){
      Char_t sfilenamehist[1024];
      sprintf(sfilenamehist,"./display-planehists-%d-%d.eps",nEvent,nEventHijingCounter);
      c1->SaveAs(sfilenamehist);
    }
    else if(c=='q') break;

    nEvent++;
    nEventHijingCounter++;
    ev->Reset();
    mixedev->Reset();
  } //end of nev loop


  delete ev;
  delete theTree;
  delete backtheTree;
  delete backev;
}
