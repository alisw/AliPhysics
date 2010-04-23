
//----------------------------------------------------------------
//      Implementation of Class AliAnalysisTaskdNdetaMC
//
// Task used to analize simulations at generation level (i.e. only
// needs galice.root and Kinematics.root).
// 
// The tasks produces multiplicity, eta and pt histograms for
// different event classes (listed in the enum in the header file).
// The multiplicity histogram are further produced for different eta
// ranges, defined in the constructor.
//
// Author: Michele Floris, CERN
//
//----------------------------------------------------------------


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

#include "AliAnalysisTaskdNdetaMC.h"
#include "TGraphErrors.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include <iostream>
#include "AliPDG.h"
#include "AliGenDPMjetEventHeader.h"
#include "TFile.h"

using namespace std;

ClassImp(AliAnalysisTaskdNdetaMC)

Float_t      AliAnalysisTaskdNdetaMC::fEtaMax = 0.5;
Int_t        AliAnalysisTaskdNdetaMC::fPDGCodes[]  = {211,2212,321,-211,-2212,-321} ;
const char * AliAnalysisTaskdNdetaMC::fPartNames[] = {"PionPos", "ProtonPos", "KaonPos", "PionNeg", "ProtonNeg", "KaonNeg"} ;


AliAnalysisTaskdNdetaMC::AliAnalysisTaskdNdetaMC() 
  : AliAnalysisTaskSE(), fNchDens(0), fMyOut(0), 
    fHistIev(0), fHistNParticlesAtMidRapidity(0), fSkipNormalization(0)
{

  // default constructor
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    fHistEta[ihist]  = 0;
    fHistPt[ihist]  = 0;
    for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++){
      fHistMult[ihist][ihist2] = 0;    
    }
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      fHistPtID[ihist][ipart] = 0;
    }
  }  

  


  fEtaBins[kEta05] = 0.5;
  fEtaBins[kEta10] = 1.0;
  fEtaBins[kEta14] = 1.3;

  fEtaMax=0.5;

}

//________________________________________________________________________
AliAnalysisTaskdNdetaMC::AliAnalysisTaskdNdetaMC(const char *name) 
  : AliAnalysisTaskSE(name), fNchDens(0), fMyOut(0),
    fHistIev(0), fHistNParticlesAtMidRapidity(0), fSkipNormalization(0)
{
  // constructor

  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    fHistEta[ihist]  = 0;
    fHistPt[ihist]  = 0;
    for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++){
      fHistMult[ihist][ihist2] = 0;    
    }
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      fHistPtID[ihist][ipart] = 0;
    }
  }
  fEtaBins[kEta05] = 0.5;
  fEtaBins[kEta10] = 1.0;
  fEtaBins[kEta14] = 1.3;


  fEtaMax=0.5;

  AliPDG::AddParticlesToPdgDataBase();

  DefineOutput(1, TList::Class());
}

AliAnalysisTaskdNdetaMC::~AliAnalysisTaskdNdetaMC() {

  // destructor

  if(fMyOut) {
    // fMyOut owns the histos
    delete fMyOut;
    fMyOut = 0;
  }

}


AliAnalysisTaskdNdetaMC::AliAnalysisTaskdNdetaMC(const char *name, const char * fname) 
  : AliAnalysisTaskSE(name), fNchDens(0), fMyOut(0),
    fHistIev(0), fHistNParticlesAtMidRapidity(0), fSkipNormalization(0)
{
  // This constructor open list from a dndeta file (useful to finalize after merging)
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    fHistEta[ihist]  = 0;
    fHistPt[ihist]  = 0;
    for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++){
      fHistMult[ihist][ihist2] = 0;    
    }
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      fHistPtID[ihist][ipart] = 0;
    }
  }
  fEtaBins[kEta05] = 0.5;
  fEtaBins[kEta10] = 1.0;
  fEtaBins[kEta14] = 1.3;


  fEtaMax=0.5;

  AliPDG::AddParticlesToPdgDataBase();

  TFile * f =new TFile(fname);
  if (!f) AliFatal(Form("Cannot open file %s!",fname));
  fMyOut  = (TList*) f->Get("coutput");
  if (!fMyOut) AliFatal(Form("Cannot get output from file %s!",fname));

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskdNdetaMC::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fMyOut = new TList();

  fHistEta[kHistINEL] = BookHetaHist("fHistdNdetaMCInel", "MC dN/deta distribution Inel");
  fHistEta[kHistNSD]  = BookHetaHist("fHistdNdetaMCNSD",  "MC dN/deta distribution NSD");
  fHistEta[kHistSiD]  = BookHetaHist("fHistdNdetaMCSiD",  "MC dN/deta distribution SiD");
  fHistEta[kHistND]   = BookHetaHist("fHistdNdetaMCND",   "MC dN/deta distribution Non-Diffractive");

  fHistEta[kHistHL]   = BookHetaHist("fHistdNdetaMCHL",     "MC dN/deta distribution, at least 1 particle |eta| < 1");


  fHistPt[kHistINEL] = BookHptHist("fHistdNdptMCInel", "MC dN/dpt distribution Inel (|#eta| < 0.8)");
  fHistPt[kHistNSD]  = BookHptHist("fHistdNdptMCNSD",  "MC dN/dpt distribution NSD (|#eta| < 0.8)");
  fHistPt[kHistSiD]  = BookHptHist("fHistdNdptMCSiD",  "MC dN/dpt distribution SiD (|#eta| < 0.8)");
  fHistPt[kHistND]   = BookHptHist("fHistdNdptMCND",   "MC dN/dpt distribution Non-Diffractive (|#eta| < 0.8)");
  fHistPt[kHistHL]   = BookHptHist("fHistdNdptMCHL",   "MC dN/dpt distribution at least 1 particle |eta| < 1 (|#eta| < 0.8)");


  const char * labelType[] = {"INEL", "NSD", "SiD", "ND", "HL"};
  for(Int_t ihist = 0; ihist < kNHist; ihist++){ // type
    for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++) { // eta range
      fHistMult[ihist][ihist2] = BookMultHisto(Form("fHistMult_%s_%1.1f",labelType[ihist],fEtaBins[ihist2]),
					       Form("(dN/dN_{ch})_{|#eta_{max}|<%1.1f} (%s)",fEtaBins[ihist2],labelType[ihist]));    
    }
  }
  for(Int_t ihist = 0; ihist < kNHist; ihist++){ // type
    for(Int_t ipart = 0; ipart < kNPart; ipart++){ // particle
      fHistPtID[ihist][ipart] = BookHptHist(Form("fHistPtID_%s_%s",labelType[ihist],fPartNames[ipart]),
					    Form("fHistPtID (%s), %s - |y| < 0.5",labelType[ihist],fPartNames[ipart]));
    }
  }
  
  
  fNchDens = new TGraphErrors();
  fNchDens -> SetName  ("fNchDens");
  fNchDens -> SetTitle ("Charged tracks density at mid-rapidity (|#eta| < fEtaMax)");
  fNchDens->SetMarkerStyle(kFullCircle);
  fMyOut->Add(fNchDens);

  fHistNParticlesAtMidRapidity = new TH1I("fHistNParticlesAtMidRapidity","Number of particles at midrapidity", kNHist, -0.5, kNHist-0.5); 
  fMyOut->Add(fHistNParticlesAtMidRapidity);

  fHistIev = new TH1I("fHistIev","Number of particles at midrapidity", kNHist, -0.5, kNHist-0.5); 
  fMyOut->Add(fHistIev);


  fMyOut->SetOwner();

  // Suppress annoying printout
  AliLog::SetGlobalLogLevel(AliLog::kError);

}

TH1F* AliAnalysisTaskdNdetaMC::BookHetaHist(const char * name, const char * title) {

  // Book Eta histos
  TH1F * h = new TH1F(name, title, 200, -10., 10.);
  h->GetXaxis()->SetTitle("#eta");
  h->GetYaxis()->SetTitle("dN/d#eta");
  h->SetMarkerStyle(kFullCircle);
  h->Sumw2();
  fMyOut->Add(h);
  return h;

}

TH1F* AliAnalysisTaskdNdetaMC::BookHptHist(const char * name, const char * title) {

  // Book Pt histos
  TH1F * h = new TH1F(name, title, 200, 0., 20.);
  h->GetXaxis()->SetTitle("p_{T}");
  h->GetYaxis()->SetTitle("dN/dp_{T}");
  h->SetMarkerStyle(kFullCircle);
  h->Sumw2();
  fMyOut->Add(h);
  return h;

}

TH1F* AliAnalysisTaskdNdetaMC::BookMultHisto(const char * name, const char * title) {

  // Book multiplicity histos
  Int_t maxmult = 100;
  TH1F * h = new TH1F(name, title, maxmult, -0.5, maxmult-0.5);
  h->GetXaxis()->SetTitle("N_{ch}");
  h->GetYaxis()->SetTitle("dN/dN_{ch}");
  h->SetMarkerStyle(kFullCircle);
  h->Sumw2();
  fMyOut->Add(h);
  return h;

}

//________________________________________________________________________
void AliAnalysisTaskdNdetaMC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // also a AliEvent...
  //  AliVEvent* mcEvent = MCEvent();
  AliMCEvent* mcEvent = MCEvent();
  AliGenPythiaEventHeader * headPy  = 0;
  AliGenDPMjetEventHeader * headPho = 0;
  AliGenEventHeader * htmp = mcEvent->GenEventHeader();
  if(!htmp) {
    AliError("Cannot Get MC Header!!");
    return;
  }
  if( TString(htmp->IsA()->GetName()) == "AliGenPythiaEventHeader") {
    headPy =  (AliGenPythiaEventHeader*) htmp;
  } else if (TString(htmp->IsA()->GetName()) == "AliGenDPMjetEventHeader") {
    headPho = (AliGenDPMjetEventHeader*) htmp;
  } else {
    cout << "Unknown header" << endl;
    
  }

  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  //  Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
  //  Check if the evend is single diffractive

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

  // HL definition: is there at least one particle in |eta|<1?
  Bool_t isThereOneCentralPart = kFALSE;
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    Bool_t isPrimary = mcEvent->Stack()->IsPhysicalPrimary(iTrack);
    if (isPrimary && track->Charge() != 0){
      if (track->Eta() > -1 && track->Eta() < 1) {
	isThereOneCentralPart = kTRUE;
	break;
      }
    }
  }

  

  fHistIev->Fill(kHistINEL);
  if (!isSD)                 fHistIev->Fill(kHistNSD);
  if (isSD)                  fHistIev->Fill(kHistSiD);
  if (isND)                  fHistIev->Fill(kHistND);
  if (isThereOneCentralPart) fHistIev->Fill(kHistHL);
  if(!(Int_t(fHistIev->GetBinContent(fHistIev->FindBin(kHistINEL)))%500)) 
    cout << "Event " << Int_t(fHistIev->GetBinContent(fHistIev->FindBin(kHistINEL))) << endl;
  

  static const Float_t ymax =0.5;
  // Track loop
  Int_t multiplicity[kNHist][kNEtaHist] = {{0}};  
  for (Int_t iTrack = 0; iTrack < mcEvent->GetNumberOfTracks(); iTrack++) {
    AliMCParticle *track = (AliMCParticle*)mcEvent->GetTrack(iTrack);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTrack);
      continue;
    }
    Bool_t isPrimary = mcEvent->Stack()->IsPhysicalPrimary(iTrack);
    if (isPrimary && track->Charge() != 0){
      Bool_t isEtaLess08 = (track->Eta() > -0.8 && track->Eta() < 0.8);

      fHistEta[kHistINEL]->Fill(track->Eta());            
      if (isEtaLess08) fHistPt[kHistINEL] ->Fill(track->Pt());            
      if(track->Eta() > -fEtaMax && track->Eta() < fEtaMax) fHistNParticlesAtMidRapidity->Fill(kHistINEL);

      for(Int_t ipart = 0; ipart < kNPart; ipart++){
	if(track->Y() > -ymax && track->Y() < ymax && track->PdgCode() == fPDGCodes[ipart])  fHistPtID[kHistINEL][ipart]->Fill(track->Pt());
      }
      


      if(!isSD) {
	fHistEta[kHistNSD]->Fill(track->Eta());      
	if (isEtaLess08) fHistPt [kHistNSD]->Fill(track->Pt());      
	if(track->Eta() > -fEtaMax && track->Eta() < fEtaMax) fHistNParticlesAtMidRapidity->Fill(kHistNSD);
	for(Int_t ipart = 0; ipart < kNPart; ipart++){
	  if(track->Y() > -ymax && track->Y() < ymax && track->PdgCode() == fPDGCodes[ipart])  fHistPtID[kHistNSD][ipart]->Fill(track->Pt());
	}
      }
      if(isSD) {
	fHistEta[kHistSiD]->Fill(track->Eta());      
	if (isEtaLess08) fHistPt [kHistSiD]->Fill(track->Pt());      
	if(track->Eta() > -fEtaMax && track->Eta() < fEtaMax) fHistNParticlesAtMidRapidity->Fill(kHistSiD);
	for(Int_t ipart = 0; ipart < kNPart; ipart++){
	  if(track->Y() > -ymax && track->Y() < ymax && track->PdgCode() == fPDGCodes[ipart])  fHistPtID[kHistSiD][ipart]->Fill(track->Pt());
	}
      }
      if (isND) {
	fHistEta[kHistND]->Fill(track->Eta());      
	if (isEtaLess08) fHistPt [kHistND]->Fill(track->Pt());      
	if(track->Eta() > -fEtaMax && track->Eta() < fEtaMax) fHistNParticlesAtMidRapidity->Fill(kHistND);
		for(Int_t ipart = 0; ipart < kNPart; ipart++){
	  if(track->Y() > -ymax && track->Y() < ymax && track->PdgCode() == fPDGCodes[ipart])  fHistPtID[kHistND][ipart]->Fill(track->Pt());
	}

      }
      if (isThereOneCentralPart) {
	fHistEta[kHistHL]->Fill(track->Eta());      
	if (isEtaLess08) fHistPt [kHistHL]->Fill(track->Pt());      
	if(track->Eta() > -fEtaMax && track->Eta() < fEtaMax) fHistNParticlesAtMidRapidity->Fill(kHistHL);
	for(Int_t ipart = 0; ipart < kNPart; ipart++){
	  if(track->Y() > -ymax && track->Y() < ymax && track->PdgCode() == fPDGCodes[ipart])  fHistPtID[kHistHL][ipart]->Fill(track->Pt());
	}

      }

      // fill array of multiplicity for different classes of events
      // and in different eta ranges
      for(Int_t ihist = 0; ihist < kNEtaHist; ihist++){
	if(track->Eta() > -fEtaBins[ihist] && track->Eta() < fEtaBins[ihist]) {
	  multiplicity[kHistINEL][ihist]++;
	  if(!isSD)                  multiplicity[kHistNSD][ihist]++;
	  if(isSD)                   multiplicity[kHistSiD][ihist]++;
	  if(isND)                   multiplicity[kHistND] [ihist]++;
	  if(isThereOneCentralPart)  multiplicity[kHistHL] [ihist]++;
  
	}
      }
    }
  } //track loop 



  // Fill multiplicity histos
  for(Int_t ihist = 0; ihist < kNEtaHist; ihist++){
    fHistMult[kHistINEL][ihist] ->Fill(multiplicity[kHistINEL][ihist]);
    if(!isSD)                  fHistMult[kHistNSD][ihist]->Fill(multiplicity[kHistNSD][ihist]);
    if(isSD)                   fHistMult[kHistSiD][ihist]->Fill(multiplicity[kHistSiD][ihist]);
    if(isND)                   fHistMult[kHistND] [ihist]->Fill(multiplicity[kHistND] [ihist]);
    if(isThereOneCentralPart)  fHistMult[kHistHL] [ihist]->Fill(multiplicity[kHistHL] [ihist]);    
  }

  // Post output data.
  PostData(1, fMyOut);
}      

//________________________________________________________________________
void AliAnalysisTaskdNdetaMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query   

  //  fHistPt = dynamic_cast<TH1F*> (GetOutputData(1));
  fMyOut  = dynamic_cast<TList*> (GetOutputData(1));
  
  Finalize();

}



void AliAnalysisTaskdNdetaMC::Finalize() {

  // Scale histos and computes dNdeta

  fHistEta[kHistINEL] = (TH1F*) fMyOut->FindObject("fHistdNdetaMCInel");
  fHistEta[kHistNSD]  = (TH1F*) fMyOut->FindObject("fHistdNdetaMCNSD" );
  fHistEta[kHistSiD]  = (TH1F*) fMyOut->FindObject("fHistdNdetaMCSiD" );
  fHistEta[kHistND]   = (TH1F*) fMyOut->FindObject("fHistdNdetaMCND"  );
  fHistEta[kHistHL]   = (TH1F*) fMyOut->FindObject("fHistdNdetaMCHL"  );
  fHistPt[kHistINEL] = (TH1F*) fMyOut->FindObject("fHistdNdptMCInel");
  fHistPt[kHistNSD]  = (TH1F*) fMyOut->FindObject("fHistdNdptMCNSD" );
  fHistPt[kHistSiD]  = (TH1F*) fMyOut->FindObject("fHistdNdptMCSiD" );
  fHistPt[kHistND]   = (TH1F*) fMyOut->FindObject("fHistdNdptMCND"  );
  fHistPt[kHistHL]   = (TH1F*) fMyOut->FindObject("fHistdNdptMCHL"  );
  fNchDens            = (TGraphErrors*) fMyOut->FindObject("fNchDens");

  const char * labelType[] = {"INEL", "NSD", "SiD", "ND", "HL"};
  for(Int_t ihist = 0; ihist < kNHist; ihist++){
    for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++){
      fHistMult[ihist][ihist2] = (TH1F*) fMyOut->FindObject(Form("fHistMult_%s_%1.1f",labelType[ihist],fEtaBins[ihist2]));      
      if (!fHistMult[ihist][ihist2]) cout << "Cannot get histo " << Form("fHistMult_%s_%1.1f",labelType[ihist],fEtaBins[ihist2]) << endl;
    }
    for(Int_t ipart = 0; ipart < kNPart; ipart++){
      fHistPtID[ihist][ipart] = (TH1F*) fMyOut->FindObject(Form("fHistPtID_%s_%s",labelType[ihist],fPartNames[ipart]));      
      if(!fHistPtID[ihist][ipart]) AliWarning(Form("Cannot get histo fHistPtID_%s_%s",labelType[ihist],fPartNames[ipart]));
    }
  }

  fHistIev  = (TH1I*) fMyOut->FindObject("fHistIev"  );
  fHistNParticlesAtMidRapidity  = (TH1I*) fMyOut->FindObject("fHistNParticlesAtMidRapidity"  );
    

//   fHistIev->Draw();

  cout << "Eta Max: " << fEtaMax << endl;


  for(Int_t ihist = 0; ihist < kNHist; ihist++){

    Int_t iev = (Int_t) fHistIev->GetBinContent(fHistIev->FindBin(ihist));
    Int_t npm = (Int_t) fHistNParticlesAtMidRapidity->GetBinContent(fHistNParticlesAtMidRapidity->FindBin(ihist));

    // Renormalize dNdeta to the number of events
    //    cout << fHistEta[ihist] << " " << iev << endl;
    
    // compute density at midrapidity (Delta_y = 1);
    Int_t firstbin =  fHistEta[ihist]->FindBin(-fEtaMax);
    Int_t lastbin =  fHistEta[ihist]->FindBin(fEtaMax);

    if (!fSkipNormalization) {
      fHistEta[ihist]->Scale(1./iev, "width");       
      fHistPt[ihist] ->Scale(1./iev, "width");       
      for(Int_t ihist2 = 0; ihist2 < kNEtaHist; ihist2++){
	fHistMult[ihist][ihist2]->Scale(1./iev);
      }
      for(Int_t ipart = 0; ipart < kNPart; ipart++){
	fHistPtID[ihist][ipart]->Scale(1./iev, "width");
      }
      
    }


    Double_t meaneta = 0;
    Double_t meanerr = 0;
    Double_t sumweight  = 0;
    for(Int_t ibin = firstbin; ibin <=lastbin ; ibin++){
      Double_t x    = fHistEta[ihist]->GetBinContent(ibin);
      Double_t xerr = fHistEta[ihist]->GetBinError(ibin);
      
      
      Double_t xerr2 = xerr*xerr;
      if(xerr2){
	//      cout << "xe2 " << xerr2 << endl;
	Double_t weight = 1. / xerr2;
	sumweight += weight;
	meaneta += weight * x;
      }
      
    }
    
    if(sumweight){
      meaneta /= sumweight;
      meanerr = TMath::Sqrt(1./ sumweight);
    }
    else {
      meaneta = 0;
      meanerr = 0;
    }


    Double_t range = 2*fEtaMax;

//     meaneta /= fHistEta[ihist]->GetBinWidth(firstbin);
//     meanerr /= fHistEta[ihist]->GetBinWidth(firstbin);
    cout << "Histo: " << fHistEta[ihist]->GetName() << endl;
    cout << " Evts: " << iev << endl;
    
    cout << " Density at midrapidity:       " << meaneta << "+-" << meanerr << endl;

    // Direct computation
    Double_t errdir = TMath::Sqrt(npm) / iev;
    Double_t etadir = Double_t(npm) / iev;
    
    cout << " Density at midrapidity: (dir) " << etadir/range << "+-" << errdir/range << endl;
  
    fNchDens->SetPoint     (ihist, ihist+1, etadir);
    fNchDens->SetPointError(ihist, 0, errdir);

    cout << " Density at midrapidity: (TH1, Eta<0.5) " << fHistMult[ihist][0]->GetMean() << "+-" << fHistMult[ihist][0]->GetMeanError() << endl;

  }
  
  // Draw fHistEtaAll
  TCanvas *c1 = new TCanvas("AliAnalysisTaskdNdetaMC","dNdetaMC",10,10,800,400);
  //  c1->cd(1)->SetLogy();
  c1->Divide(2,1);
  c1->cd(1);
  fHistEta[0]->DrawCopy("");
  fHistEta[1]->SetLineColor(kRed);
  fHistEta[1]->SetMarkerColor(kRed);
  fHistEta[1]->SetMarkerStyle(kOpenSquare);
  fHistEta[1]->DrawCopy("same");
  c1->cd(2);
  fNchDens->Draw("AP");

}
