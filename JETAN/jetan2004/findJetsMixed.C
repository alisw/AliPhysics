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
#include <TStopwatch.h>
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
#endif

void findJetsMixed(Char_t *fname,Char_t *hname,Float_t ptcut=2.0,Float_t radius=0.3,Char_t *savename=0,Int_t nMaxEvents=-1)
{
  //connect to jets
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

  // create the jet finder
  Char_t buffer[8096];
  AliTkConeJetFinderV2 *ConeFinder = new AliTkConeJetFinderV2();
  ConeFinder->defaultSettings();
  //ConeFinder->setSettings(120,40);
  if(gErrorIgnoreLevel<=kWarning)
    ConeFinder->setOutput(kTRUE);
  ConeFinder->setEtMinJet(10);
  ConeFinder->setPtCut(ptcut);
  ConeFinder->setEtCut(ptcut);
  ConeFinder->setRadius(radius);
  if(savename) sprintf(buffer,"%s",savename);
  else {
    Char_t buffer2[8096];
    strncpy(buffer2,fname,strlen(fname)-5);
    buffer2[strlen(fname)-5]='\0';
    sprintf(buffer,"%s-mixed-thresh-%.1f-rad-%.1f.cone.evout.root",buffer2,ptcut,radius);
  }
  ConeFinder->setEvOutFilename(buffer);
  ConeFinder->init();

  //=========================================================================
  // start the event loop
  //=========================================================================
  TStopwatch *stopwatch=new TStopwatch();
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
    backev->SetHeader(dummy);
    if(gErrorIgnoreLevel<=kWarning){
      cout << "Read event: " << nEvent << " " << nEventHijing << endl;
      backev->Print();
    }

    stopwatch->Reset();
    stopwatch->Start();
    ConeFinder->initEvent(backev,dummy);
    ConeFinder->run();
    ConeFinder->finishEvent();
    if(gErrorIgnoreLevel<=kWarning){
      cout << "JetFinding done:  CPU time " << stopwatch->CpuTime()
	   << " Real Time " << stopwatch->RealTime() << endl;
    }

    nEvent++;
    nEventHijingCounter++;
    ev->Reset();
  } //end of nev loop

  ConeFinder->finish();

  delete ev;
  delete theTree;
  delete backev;
  delete backtheTree;
  delete ConeFinder;
}

