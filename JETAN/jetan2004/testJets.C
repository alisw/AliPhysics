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
#include <AliTkConeJetEvent.h>
#include <AliTkConeJetFinderV2.h>
#include <AliJetParticle.h>
#include <AliJetParticlesReader.h>
#include <AliJetParticlesReaderKine.h>
#include <AliJetParticlesReaderESD.h>
#include <AliJetParticlesReaderHLT.h>
#include <AliJetEventParticles.h>
#endif

#if 0
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kPHOSpid=0x10000, kHMPIDpid=0x20000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
#endif

void testJets(Char_t *input,Int_t esd=0,Float_t ptcut=2.0,Float_t radius=0.3)
{

  TObjArray *dirs=new TObjArray(10);
  dirs->Add(new TObjString(input));

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
    reader=new AliJetParticlesReaderESD(0,dirs);
    ((AliJetParticlesReaderESD*)reader)->SetCompareFlag(0x0030);
    desc+="ESD";
  } else if (esd==2) {
    reader=new AliJetParticlesReaderHLT(kTRUE,dirs);
    desc+="HLTConf";
  } else if (esd==3) {
    reader=new AliJetParticlesReaderHLT(kFALSE,dirs);
    desc+="HLTHough";
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
  } else {
    cout << "Parameter settings: " << endl;
    cout << "0  == Kine (ch)"  << endl;
    cout << "1  == TPC"  << endl;
    cout << "2  == Conformal"  << endl;
    cout << "3  == Hough"  << endl;
    cout << "10 == Kine (all)"  << endl;
    cout << "11 == Kine (ch+em)"  << endl;
    return;
  }
  //reader->ReadEventsFromTo(0,10);
  reader->SetPtCut(ptcut,100);
  reader->SetEtaCut(-0.9,0.9);

  desc+=" Ptcut: ";
  desc+=ptcut;

  // create the jet finder
  Char_t buffer[8096];
  AliTkConeJetFinderV2 *ConeFinder = new AliTkConeJetFinderV2();
  ConeFinder->defaultSettings();
  //ConeFinder->setSettings(120,40);
  if(gErrorIgnoreLevel<=kWarning)
    ConeFinder->setOutput(kTRUE);
  ConeFinder->setEtMinJet(1.);
  ConeFinder->setEtCut(-1);
  ConeFinder->setRadius(radius);
  sprintf(buffer,"./conefinder-type%d.evout.root",esd);
  //cout << "Cone Jet Finder init: " << endl;
  ConeFinder->setEvOutFilename(buffer);
  ConeFinder->init();

  Int_t counter=0;
  while(reader->Next())
    {
      const AliJetEventParticles *ev=reader->GetEventParticles();
      if(gErrorIgnoreLevel<=kWarning){
	cout << "Read event: " << counter << endl;
	ev->Print();
      }
      TString dummy="Counter: ";
      dummy+=counter;
      dummy+=" ";
      dummy+=desc;
      ConeFinder->initEvent(ev,dummy);
      ConeFinder->run();
      ConeFinder->finishEvent();
    }
  ConeFinder->finish();

  delete ConeFinder;
  delete reader;
}

void readJets(Int_t nMaxEvents,Char_t *filename)
{
  //connect to jets
  TChain *theTree = new TChain("jets");
  theTree->Add(filename);
  AliTkConeJetEvent *event = new AliTkConeJetEvent();
  theTree->SetBranchAddress("ConeFinder",&event);

  Int_t treeentries=(Int_t)theTree->GetEntries();
  if((nMaxEvents<0) || (nMaxEvents>treeentries))
    nMaxEvents=treeentries;

  cout << "Found " << nMaxEvents << " in " << filename << endl;

  //=========================================================================
  // start the event loop
  //=========================================================================
  Int_t nEvent = 0;
  while(nEvent<nMaxEvents){
    if ((nEvent % 100) == 0) {
      cout << "Analysing event " << nEvent << endl;
    }

    //connect the cone jets
    theTree->GetEvent(nEvent);

    AliJetEventParticles *jetparts=event->getJetParticles();
    jetparts->Print();
#if 0
    const TClonesArray *p=jetparts->GetParticles();
    for(Int_t i=0;i<p->GetEntriesFast();i++) {
      cout << i << endl;
      ((AliJetParticle*)p->At(i))->Print();
    }
#endif

    TClonesArray *tkjets=event->getJets();
    if(!tkjets){
      cerr << "No Cone jet found in event " << nEvent << ", exiting..." <<endl;
      continue;
    }

    event->Print("");
    nEvent++;
    event->Clear();
  } //end of nev loop

  delete event;
  delete theTree;
}
