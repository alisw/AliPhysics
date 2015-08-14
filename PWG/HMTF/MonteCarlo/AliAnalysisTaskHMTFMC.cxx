#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAODEvent.h"

#include "AliAnalysisTaskHMTFMC.h"
#include "TGraphErrors.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHepMCEventHeader.h"

#include <iostream>
#include "AliPDG.h"
#include "AliGenDPMjetEventHeader.h"
#include "TH2F.h"

using namespace std;

ClassImp(AliAnalysisTaskHMTFMC)

AliAnalysisTaskHMTFMC::AliAnalysisTaskHMTFMC() 
: AliAnalysisTaskSE(), fHistEta(0), fHistNch(0), fHistNchUnweighted(0), fHistRawMult(0x0), fHistIev(0), fMyOut(0),
  fPrimaryPDGs(), fMotherPDGs()
{

}

//________________________________________________________________________
AliAnalysisTaskHMTFMC::AliAnalysisTaskHMTFMC(const char *name) 
  : AliAnalysisTaskSE(name), fHistEta(0), fHistNch(0), fHistNchUnweighted(0), fHistRawMult(0x0), fHistIev(0), fMyOut(0),
  fPrimaryPDGs(), fMotherPDGs()
{

  AliPDG::AddParticlesToPdgDataBase();

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fMyOut = new TList();

  fHistEta  = new TH1F("fHistdNdetaMCInel", "dN/d#eta inel;#eta;dN/d#eta", 100, -2., 2.);
  fHistEta->SetMarkerStyle(kFullCircle);
  fHistEta->Sumw2();
  fMyOut->Add(fHistEta);

  fHistNch = new TH1F ("fHistNch", "Multiplicity distribution, #||{#eta} < 1.0;N;counts", 100, -0.5, 99.5);
  fHistNch->Sumw2();
  fMyOut->Add(fHistNch);

  fHistNchUnweighted = new TH1F ("fHistNchUnweighted", "Multiplicity distribution, #||{#eta} < 1.0 (unweighted);N;counts", 100, -0.5, 99.5);
  fHistNchUnweighted->Sumw2();
  fMyOut->Add(fHistNchUnweighted);

  fHistRawMult = new TH1F ("fHistRawMult", "Raw mult for downscaling;N;counts", 100, -0.5, 999.5);
  fHistRawMult->Sumw2();
  fMyOut->Add(fHistRawMult);

  // One needs to use a mergeable object to count events. In this case
  // we use a TH1D (careful: a smaller type, e.g. TH1I could
  // saturate). If you split events in different classes (inel, nsd,
  // INEL>0, ...) you can add more bins here
  //
  // Here we count 2 classes of events: all processed events, and the
  // equivalent generate statistics
  fHistIev = new TH1D("fHistIev", "event statistics;;counts", 2, -0.5, 1.5);
  fHistIev->GetXaxis()->SetBinLabel(1, "events");
  fHistIev->GetXaxis()->SetBinLabel(2, "weight");
  fMyOut->Add(fHistIev);

  fMyOut->SetOwner();
  
  // Suppress annoying printout
  AliLog::SetGlobalLogLevel(AliLog::kError);
  PostData(1, fMyOut);
}

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // also a AliEvent...
  AliMCEvent* mcEvent = MCEvent();

  // Here we need some generator-dependent logic, in case we need to extract useful information from the headers.
  AliGenPythiaEventHeader * headPy  = 0;
  AliGenDPMjetEventHeader * headPho = 0;
  AliGenHepMCEventHeader * headHepMC = 0x0;
  AliGenEventHeader * htmp = mcEvent->GenEventHeader();
  if(!htmp) {
    AliError("Cannot Get MC Header!!");
    return;
  }
  if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
    headPy =  (AliGenPythiaEventHeader*) htmp;
  } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
    headPho = (AliGenDPMjetEventHeader*) htmp;
  } else if (TString(htmp->IsA()->GetName()) == "AliGenHepMCEventHeader") {
    headHepMC = (AliGenHepMCEventHeader*) htmp;
  } else {
    AliWarning("Unknown header");
  }

  AliHeader* header = mcEvent->Header();
  AliStack* stack = header->Stack();

  // plot multiplicity used for downscaling
  // fHistRawMult->Fill(htmp->NProduced());
  fHistRawMult->Fill(stack->GetNtransported());

  // In the downscaled productions, every event is given a weight which corresponds to the downscaling factor applied. This is needed to fill natural-looking histos
  Float_t evWeight=htmp->EventWeight();

  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  // This is left as an example, not really used in the present task. Could be useful to distinguish different classes of events
  Bool_t isSD = kFALSE;
  Bool_t isND = kFALSE;
  if(headPy)   {
    //    cout << "Process: " << headPy->ProcessType() << endl;
    if(headPy->ProcessType() == 92 || headPy->ProcessType() == 93) {
      isSD = kTRUE; // is single difractive
    }
    if(headPy->ProcessType() != 92 && headPy->ProcessType() != 93 && headPy->ProcessType() != 94) {     
      isND = kTRUE; // is non-diffractive
    }

  } else if (headPho) {
    if(headPho->ProcessType() == 5 || headPho->ProcessType() == 6 ) {
      isSD = kTRUE;
    }       
    if(headPho->ProcessType() != 5 && headPho->ProcessType() != 6  && headPho->ProcessType() != 7 ) {
      isND = kTRUE;
    }       
  }

  
  // Track loop
  Int_t nchMidEta = 0;
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {    
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    Bool_t isPrimary = mcEvent->Stack()->IsPhysicalPrimary(iTrack);
    if (isPrimary && track->Charge() != 0){

      if(TMath::Abs(track->Eta()) < 1.0) nchMidEta++;
      
      fHistEta->Fill(track->Eta(),evWeight);
      fPrimaryPDGs[abs(track->PdgCode())]++;
      fMotherPDGs [abs(mcEvent->GetTrack(track->GetMother())->PdgCode())]++;
    }
  } //track loop 
  if(nchMidEta == 0) return; //INEL > 0 request
  fHistNch->Fill(nchMidEta,evWeight);
  fHistNchUnweighted->Fill(nchMidEta);

  fHistIev->Fill(0);           // number of processed events
  fHistIev->Fill(1, evWeight); // sum of weights

  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskHMTFMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query   
  //  return;
  //  fHistPt = dynamic_cast<TH1F*> (GetOutputData(1));
  fMyOut  = dynamic_cast<TList*> (GetOutputData(1));
  
  fHistIev  = (TH1D*) fMyOut->FindObject("fHistIev");

  // Draw maps
  TH1I * hPrimaryPDG = new TH1I ("hPrimaryPDG", "hPrimaryPDG", fPrimaryPDGs.size(), 0, fPrimaryPDGs.size());
  Int_t ibin = 0;
  for (std::map<int,int>::iterator it=fPrimaryPDGs.begin(); it!=fPrimaryPDGs.end(); ++it){
    
    ibin++;
    std::cout << ibin <<  " -> " << it->first <<":" << it->second << std::endl;
    hPrimaryPDG->SetBinContent(ibin, it->second);
    hPrimaryPDG->GetXaxis()->SetBinLabel(ibin, Form("%d", it->first));
  }

  TH1I * hMotherPDG = new TH1I ("hMotherPDG", "hMotherPDG", fMotherPDGs.size(), 0, fMotherPDGs.size());
  ibin = 0;
  for (std::map<int,int>::iterator it=fMotherPDGs.begin(); it!=fMotherPDGs.end(); ++it){
    ibin++;
    hMotherPDG->SetBinContent(ibin, it->second);
    hMotherPDG->GetXaxis()->SetBinLabel(ibin, Form("%d", it->first));
  }

  TCanvas * c = new TCanvas("cPDG", "cPDG");
  c->Divide(2,1);
  c->cd(1);
  hPrimaryPDG->Draw();
  c->cd(2);
  hMotherPDG->Draw();

  std::cout << "Processed " << fHistIev->GetBinContent(1) << " events, equivalent to "<< fHistIev->GetBinContent(2) << " generated events."  << std::endl;
  
}
