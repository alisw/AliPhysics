#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "TText.h"
#include "TFile.h"
#include "TBenchmark.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliTrackReference.h"

#include "AliT0MultiplicityTask.h"

ClassImp(AliT0MultiplicityTask)


AliT0MultiplicityTask::AliT0MultiplicityTask()
  : AliAnalysisTaskSE("AliT0MultiplicityTask"),
    fListOfHistos(0),
    fOrA(0),
    fOrC(0),
    fMean(0),
    fVertex(0),
    fTime(0),
    fAmp(0),
    fTotalMult(0),
    fMultRecTotal(0),
    fMultRecRealA(0),
    fMultRecRealC(0),
    fPrim(0),
    fRef(0),
    fRec(0)
{
  // Default constructor
  AliInfo("Default constructor AliT0MultiplicityTask");
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}


AliT0MultiplicityTask::AliT0MultiplicityTask(const char* name)
  : AliAnalysisTaskSE(name),
    fListOfHistos(0),
    fOrA(0),
    fOrC(0),
    fMean(0),
    fVertex(0),
    fTime(0),
    fAmp(0),   
    fTotalMult(0),
    fMultRecTotal(0),
    fMultRecRealA(0),
    fMultRecRealC(0),
    fPrim(0),
    fRef(0),
    fRec(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #1 TList
  DefineOutput(1, TList::Class());
}


void AliT0MultiplicityTask::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  // Create output container
  fListOfHistos = new TList();
  
 
  fOrA = new TH1F("hOrA"," T0 A", 100,2450,2700);
  fOrC = new TH1F("hOrC"," T0 C", 100,2450,2700);
  fMean= new TH1F("hMean"," T0 ",100, 12000, 13000);
  fVertex = new TH1F("hVertex","Z position of vertex",   100,-30, 30);
   
  fTime = new TH1F("Time", " Amp vs Time",100, 12000, 13000);
  fAmp = new TH1F("fAmp", " Amplitude", 100, 0, 200);

  //  fHighMult  = new TH1F("fHighMult", " events with high amp ",100, 0, 10);
  
  fTotalMult = new TH1F("fTotalMult","total multiplicity",500,0,5000);
  
  fMultRecRealA = new TH2F("fMultRecRealA"," ",100,0.,200,100,0,200); 
  fMultRecRealC = new TH2F("fMultRecRealC"," ",100,0,200,100,0,200); 
  fMultRecTotal = new TH2F("fMultRecTotal"," ",100,0,200,100,0,200); 
  
  fPrim = new TH1F("fPrim", " primary",100, 0, 200);
  fRef  = new TH1F("fRef", " from TR ",100, 0, 200);
  fRec  = new TH1F("fRec", " in ESD ",100, 0, 200);
   
  fListOfHistos->Add(fOrA);
  fListOfHistos->Add(fOrC);
  fListOfHistos->Add(fMean);
  fListOfHistos->Add(fVertex);
  fListOfHistos->Add(fAmp);
  fListOfHistos->Add(fTime);
  fListOfHistos->Add(fTotalMult);
  fListOfHistos->Add(fMultRecRealA);
  fListOfHistos->Add(fMultRecRealC);
  fListOfHistos->Add(fMultRecTotal);
  fListOfHistos->Add(fPrim);
  fListOfHistos->Add(fRef);
  fListOfHistos->Add(fRec);
}  


void AliT0MultiplicityTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

    
  // MC information
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
    Printf("ERROR: Could not retrieve MC event");
    return;
  }
  
  
  Int_t primaryA=0;
  Int_t primaryC=0;
  Int_t numPrim = 0;
  Int_t refT0=0;

  AliStack* stack = mcEvent->Stack();
  // printf("AliT0MultiplicityAnalysis: Number of tracks on stack %5d\n", stack->GetNtrack());
  Int_t nTracks  = mcEvent->GetNumberOfTracks();
  for (Int_t ipart=0; ipart<nTracks; ipart++)
    {
      //    TParticle* particle = stack->Particle(ipart);
      AliMCParticle* track = mcEvent->GetTrack(ipart);
      if (!track) {
	Printf("ERROR: Could not receive track %d (mc loop)", ipart);
	continue;
      }
      Int_t label = track->GetLabel();
      if(stack->IsPhysicalPrimary(label) == kFALSE)
	continue;
      
      numPrim++;
      Double_t eta=track->Eta();
      
      if(eta<-2.97 && eta>-3.28)  primaryC++;
      if (eta >4.61 && eta<4.92)  primaryA++;
    
       // Loop over Track References
       AliTrackReference* trackRef = 0;

      for (Int_t iTrackRef = 0; iTrackRef  < track->GetNumberOfTrackReferences(); iTrackRef++) {
	trackRef = track->GetTrackReference(iTrackRef);
	if(trackRef) {
	  Int_t detectorId = trackRef->DetectorId(); 
	  if (detectorId == 7)  refT0++;
	}      
      }
    }
  fTotalMult->Fill(numPrim);
  fPrim->Fill(primaryC+primaryA);
  fRef->Fill(refT0);
  //     printf (" tracks %i primaries %d \n",nTracks, numPrim);
 

  // ESD information  
  AliVEvent* event = InputEvent();
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    return;
  }
    
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);

  Bool_t eventTriggered = false;
  ULong64_t triggerMask = esd->GetTriggerMask();
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
    eventTriggered == true;
 
    //= AliPWG0Helper::IsEventTriggered(esd, AliPWG0Helper::kMB1);
  //  printf("!!!!! eventTriggered %i",eventTriggered);
  // if(!eventTriggered) return;
   Double_t besttimeC=99999.;
  Double_t besttimeA=99999.;
  
//  Float_t coefA = 0.891466;
//  Float_t coefC = 0.922772;
  Float_t coefA = 1;
  Float_t coefC = 1;
  Float_t sumampA=0;
  Float_t sumampC=0;
  Int_t highA=0, highC=0;
  const Double_t *amp = esd->GetT0amplitude();
  const Double_t *time = esd->GetT0time();
  Float_t vertex = esd->GetT0zVertex();
  if(vertex<999) fVertex->Fill(vertex);
  Float_t   mean = esd->GetT0();
  fMean->Fill(mean);
  
  for (Int_t i=0; i<12; i++) {
    sumampC += amp[i];  
    fTime->Fill(amp[i],time[i]); 
    fAmp->Fill(amp[i]); 
    if(time[i]<besttimeC && time[i]>0) besttimeC=time[i]; //timeC
  }
  fOrC->Fill(besttimeC);
  
  for (Int_t i=12; i<24; i++){
    sumampA += amp[i];  
    fAmp->Fill(amp[i]); 
    fTime->Fill(time[i]);
    if(time[i]<besttimeA && time[i]>0) besttimeA=time[i]; //timeC
  }

  fOrA->Fill(besttimeA);
  
  fMultRecRealA ->Fill(primaryA, sumampA*coefA);
  fMultRecRealC ->Fill(primaryC, sumampC*coefC);
  fMultRecTotal->Fill(numPrim, sumampA*coefA +  sumampC*coefC);
  fRec->Fill(Int_t(sumampA+sumampC));
 
  // Post output data.
  PostData(1, fListOfHistos);
}      


void AliT0MultiplicityTask::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query 
  //  printf(" AliT0MultiplicityTask::Terminate ");
  
  fListOfHistos = dynamic_cast<TList*>(GetOutputData(1));
  if (!fListOfHistos) {
    Printf("ERROR: fListOfHistos not available");
    return;
  }
  //  printf(" before cast ");
  fOrA = dynamic_cast<TH1F*>(fListOfHistos->At(0));
  fOrC = dynamic_cast<TH1F*>(fListOfHistos->At(1));
  fMean  = dynamic_cast<TH1F*>(fListOfHistos->At(2));  
  fVertex  = dynamic_cast<TH1F*>(fListOfHistos->At(3));
  fAmp = dynamic_cast<TH1F*>(fListOfHistos->At(4));
  fTime = dynamic_cast<TH1F*>(fListOfHistos->At(5));
  fTotalMult = dynamic_cast<TH1F*>(fListOfHistos->At(6));
  fMultRecTotal = dynamic_cast<TH2F*>(fListOfHistos->At(7));
  fMultRecRealA = dynamic_cast<TH2F*>(fListOfHistos->At(8));
  fMultRecRealC= dynamic_cast<TH2F*>(fListOfHistos->At(9));
  fPrim= dynamic_cast<TH1F*>(fListOfHistos->At(10));
  fRef= dynamic_cast<TH1F*>(fListOfHistos->At(11));
  fRec= dynamic_cast<TH1F*>(fListOfHistos->At(12));

  TFile fc("MultHist.root\n", "RECREATE");
  //  printf(" File MultHist.root recreated\n");
  fOrC->Write();
  //  printf("1\n");
  fOrA->Write();
  //  printf("2\n");
  fMean->Write();
  //  printf("3\n");
  fVertex->Write();
  //  printf("4\n");
  fAmp->Write();
  //  printf("5\n");
  fTime->Write();
 fTotalMult->Write();
 //  printf("6\n");
 fMultRecTotal->Write();
 //  printf("7\n");
 fMultRecRealA->Write();
 // printf("8\n");
  fMultRecRealC->Write();
  // printf("9\n");
  fPrim->Write();
  fRef->Write();
  fRec->Write();
  
  fc.Close();
 printf(" fc.Close()\n");

}
