// $Id$

#ifndef __CINT__
#include <stdio.h>
#include <iostream.h>
#include <time.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TStopwatch.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <AliStack.h>
#include <AliRunLoader.h>
#include <AliHeader.h>
#include <AliGenPythiaEventHeader.h>
#include <AliJetParticle.h>
#include <AliJetParticlesReader.h>
#include <AliJetParticlesReaderKine.h>
#include <AliJetParticlesReaderESD.h>
#include <AliJetParticlesReaderHLT.h>
#include <AliJetEventParticles.h>
#include <AliTkConeJetEvent.h>
#include <AliTkConeJetFinderV2.h>
#include <AliTkChargedJetFinder.h>
#endif

void findJets(Char_t *filename,Float_t ptcut=1.0,Float_t radius=0.3,
	      Int_t seedpt=-1,Char_t *savename=0,Int_t nMaxEvents=-1)
{
  //connect to jets
  TChain *theTree = new TChain("AJEPtree");
  theTree->Add(filename);
  AliJetEventParticles *ev=new AliJetEventParticles();
  theTree->SetBranchAddress("particles",&ev);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;
  //cout << "Found " << nMaxEvents << " in " << filename << endl;

  // create the jet finder
  Char_t buffer[8096];
  AliTkConeJetFinderV2 *ConeFinder = new AliTkConeJetFinderV2();
  ConeFinder->defaultSettings();
  //ConeFinder->setSettings(120,40);
  if(gErrorIgnoreLevel<=kWarning)
    ConeFinder->setOutput(kTRUE);
  if(seedpt<5)
    ConeFinder->setEtMinJet(5.);
  else 
    ConeFinder->setEtMinJet(seedpt);
  ConeFinder->setPtCut(ptcut);
  ConeFinder->setEtCut(seedpt);
  ConeFinder->setRadius(radius);
  if(savename) sprintf(buffer,"%s",savename);
  else {
    Char_t buffer2[8096];
    strncpy(buffer2,filename,strlen(filename)-5);
    buffer2[strlen(filename)-5]='\0';
    sprintf(buffer,"%s-thresh-%.1f-rad-%.1f.cone.evout.root",buffer2,ptcut,radius);
  }
  ConeFinder->setEvOutFilename(buffer);
  ConeFinder->init();

  //=========================================================================
  // start the event loop
  //=========================================================================
  TStopwatch *stopwatch=new TStopwatch();
  Int_t nEvent = 0;
  while(nEvent<nMaxEvents){

    //get the event
    theTree->GetEvent(nEvent);

    //const AliJetEventParticles *ev=reader->GetEventParticles();
    if(gErrorIgnoreLevel<=kWarning){
      cout << "Read event: " << nEvent << endl;
      ev->Print();
    }
    TString dummy="Counter: ";
    dummy+=nEvent;

    stopwatch->Reset();
    stopwatch->Start();
    ConeFinder->initEvent(ev,dummy);
    ConeFinder->run();
    ConeFinder->finishEvent();
    if(gErrorIgnoreLevel<=kWarning){
      cout << "JetFinding done:  CPU time " << stopwatch->CpuTime()
	   << " Real Time " << stopwatch->RealTime() << endl;
    }
    nEvent++;
    ev->Reset();
  } //end of nev loop

  ConeFinder->finish();

  delete stopwatch;
  delete ev;
  delete theTree;
  delete ConeFinder;
}

void findChargedJets(Char_t *filename,Float_t ptcut=1.0,Float_t radius=0.3,
		     Int_t seedpt=10,Char_t *savename=0,Int_t nMaxEvents=-1)
{
  //connect to jets
  TChain *theTree = new TChain("AJEPtree");
  theTree->Add(filename);
  AliJetEventParticles *ev=new AliJetEventParticles();
  theTree->SetBranchAddress("particles",&ev);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;
  //cout << "Found " << nMaxEvents << " in " << filename << endl;

  // create the jet finder
  Char_t buffer[8096];
  AliTkChargedJetFinder *ChFinder = new AliTkChargedJetFinder();
  ChFinder->defaultSettings();
  //  if(gErrorIgnoreLevel<=kWarning)
    ChFinder->setOutput(kTRUE);
  ChFinder->setEtMinJet(seedpt);
  ChFinder->setPtCut(ptcut);
  ChFinder->setPtSeed(seedpt);
  ChFinder->setRadius(radius);
  if(savename) sprintf(buffer,"%s",savename);
  else {
    Char_t buffer2[8096];
    strncpy(buffer2,filename,strlen(filename)-5);
    buffer2[strlen(filename)-5]='\0';
    sprintf(buffer,"%s-thresh-%.1f-rad-%.1f.charged.evout.root",buffer2,ptcut,radius);
  }
  ChFinder->setEvOutFilename(buffer);
  ChFinder->init();

  //=========================================================================
  // start the event loop
  //=========================================================================
  TStopwatch *stopwatch=new TStopwatch();
  Int_t nEvent = 0;
  while(nEvent<nMaxEvents){

    //get the event
    theTree->GetEvent(nEvent);

    //  const AliJetEventParticles *ev=reader->GetEventParticles();
    if(gErrorIgnoreLevel<=kWarning){
      cout << "Read event: " << nEvent << endl;
      ev->Print();
    }
    TString dummy="Counter: ";
    dummy+=nEvent;

    stopwatch->Reset();
    stopwatch->Start();
    ChFinder->initEvent(ev,dummy);
    ChFinder->run();
    ChFinder->finishEvent();
    stopwatch->Stop();
    if(gErrorIgnoreLevel<=kWarning){
      cout << "JetFinding done: CPU time:" << stopwatch->CpuTime()
	   << " Real Time:" << stopwatch->RealTime() << endl;
    }
    nEvent++;
    ev->Reset();
  } //end of nev loop

  ChFinder->finish();

  delete stopwatch;
  delete ev;
  delete theTree;
  delete ChFinder;
}
