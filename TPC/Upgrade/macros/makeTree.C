#include <iostream> 
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include "AliToyMCEvent.h"
#include "AliToyMCEventGenerator.h"
#include "AliToyMCEventGeneratorSimple.h"
#include <TRandom.h>

//to use: root -l loadlibs.C
//    .L makeTree.C+
//    makeTree(50,10,10) for 50 kHz average frequency, 10MHz bunch crossing and 10 events. 


void makeTree(Double_t collFreq/*kHz*/, Double_t bunchFreq/*MHz*/, Int_t nEvents, Bool_t fixedFreq = kFALSE, Bool_t noVertSpread = kFALSE) {

  
  TFile* outFile = new TFile(Form("toymcevents%2.1fkHz_%d%s%s.root",collFreq,nEvents,fixedFreq?"_fixedfreq":"",noVertSpread?"_novertspread":"" ),"recreate");
  TTree* outTree = new TTree("AliToyMC","AliToyMC");
  AliToyMCEvent* event = 0x0;
  outTree->Bronch("AliToyMCEvents","AliToyMCEvent",&event);


  Double_t collProb = (collFreq*1000)/(bunchFreq*1000000);
  Int_t generatedEvents = 0;
  AliToyMCEventGeneratorSimple* evGen = new AliToyMCEventGeneratorSimple();
  if(noVertSpread)evGen->SetParameters("files/params.root",0.,0.); //mean, sigma of vertex
  else evGen->SetParameters("files/params.root",0,0.7); //mean, sigma of vertex
  //generate events
  Double_t time = 0.;
  while(generatedEvents < nEvents) {

    //draw number of collissions in crossing
    Int_t nColls = gRandom->Poisson(collProb);
    if(fixedFreq) nColls =1;
    for (Int_t iColl = 0; iColl<nColls; iColl++)       {
      
      event = evGen->Generate(time);
      outTree->Fill();
      generatedEvents++;
      delete event;
    }
    if(!fixedFreq) time += 1./(bunchFreq); //microseconds
    else time += 1./(collFreq/1000);
    ///*if(generatedEvents%10==0)*/ std::cout << generatedEvents << " time: " << time << " " << collProb << std::endl;
  }
  outFile->Write();
  outFile->Close();
}  
 
