
#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TClonesArray.h"

#include "AliModule.h"
#include "AliFRAMEv1.h"

#include "AliTRDv1.h"
#include "AliTRDhit.h"
#include "AliTRDsim.h"
#include "AliTRDsimple.h"
#include "AliTRDsimpleGen.h"
#include "AliTRDsimpleMC.h"
#include "AliTRDgeometry.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDdataArrayI.h"

#endif

void AliTRDrunSimple()
{

  ///////////////////////////////////////////////////////////////////////////////
  //
  //  Macro to run the simplified simulator
  //
  ///////////////////////////////////////////////////////////////////////////////

  //_____________________________________________________________________________ 
  //
  // Initialization part
  //_____________________________________________________________________________ 
  //

  // The simple simulator
  AliTRDsimple *simple = new AliTRDsimple();
  simple->Init();

  // Initialize a dummy frame so that the TRD works
  new AliFRAMEv1("FRAME","Space Frame");             

  // Initialize the TRD detector
  AliTRDv1     *trd    = new AliTRDv1("TRD","TRD slow simulator");
  trd->SetHitTypeStandard();

  // Needed for some material properties 
  trd->CreateMaterials();

  // Select the gas mixture (0: 97% Xe + 3% isobutane, 1: 85% Xe + 15% CO2)
  trd->SetGasMix(1);
 
  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");

  parameter->SetADCthreshold(0);
  parameter->SetNTimeBin(30);        // The number of timebins
  parameter->SetExpandTimeBin(5,15); // The additional timebins
    
  // Switch on TR production
  AliTRDsim *trdsim = trd->CreateTR();
  trdsim->SetNFoils(100);
  trdsim->SetFoilThick(0.0013);
  trdsim->SetGapThick(0.0500);

  // Initialize the TRD object
  trd->Init();

  // The digitizer
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("TRDdigitizer","Digitizer class");

  digitizer->SetParameter(parameter);
  digitizer->SetSimple();
  digitizer->InitDetector();
                                    
  // The event generator
  AliTRDsimpleGen *generator = simple->GetGenerator();         
  generator->SetMomentum(3.0,3.0);
  generator->SetPdg(11);                             // Electrons 

  //_____________________________________________________________________________ 
  //
  // Histograms
  //_____________________________________________________________________________ 
  //

  Int_t timeMax  = parameter->GetTimeTotal();   
  Int_t adcRange = ((Int_t) parameter->GetADCoutRange()); 

  TH1F     *hQ       = new TH1F("hQ"    ,"Charge per hit (all)" ,100,0.0,1000.0);
  TH1F     *hQdedx   = new TH1F("hQdedx","Charge per hit (dedx)",100,0.0,1000.0);
  TH1F     *hQtr     = new TH1F("hQtr  ","Charge per hit (tr)"  ,100,0.0,1000.0);

  TProfile *hQX      = new TProfile("hQX"    ,"Charge per hit vs X (all)" ,35,0.0,3.5);
  TProfile *hQXdedx  = new TProfile("hQXdedx","Charge per hit vs X (dedx)",35,0.0,3.5);
  TProfile *hQXtr    = new TProfile("hQXtr  ","Charge per hit vs X (tr)"  ,35,0.0,3.5);

  TH1F     *hNstep   = new TH1F("hNstep","Number of steps / cm"    , 151,-0.5, 150.5);
  TH1F     *hNel     = new TH1F("hNel"  ,"Number of electrons / cm",1001,-0.5,1000.5);

  TH1F     *hAmp     = new TH1F("hAmp","Amplitude of the digits"
                                      ,adcRange+1,-0.5,((Float_t) adcRange)+0.5);         
  TProfile *hAmpTime = new TProfile("hAmpTime","Amplitude vs timebin"
				              ,timeMax,-0.5,((Float_t) timeMax)-0.5);

  //_____________________________________________________________________________ 
  //
  // Event loop
  //_____________________________________________________________________________ 
  //

  // Number of events (tracks)
  Int_t nEvent = 10000;

  Float_t x0 = parameter->GetTime0(0) - AliTRDgeometry::DrThick(); 

  TClonesArray *hitsArray = trd->Hits();

  for (Int_t iEvent = 0; iEvent < nEvent; iEvent++) {

    if (!(iEvent % 100) && (iEvent)) {
      printf("Event no. %d\n",iEvent);
    }

    // Generate the hits for one track
    simple->ProcessEvent(iEvent);

    // Analyze the hits
    Float_t nElectrons = 0.0;
    for (Int_t iHit = 0; iHit < hitsArray->GetEntries(); iHit++) {

      AliTRDhit *hit = (AliTRDhit *) hitsArray->UncheckedAt(iHit);
      Int_t   charge = TMath::Abs(hit->GetCharge());
      Float_t x      = hit->X() - x0;
      hQ->Fill(charge);
      hQX->Fill(x,charge);
      if (hit->FromDrift() ||
          hit->FromAmplification()) {
        hQdedx->Fill(charge);
        hQXdedx->Fill(x,charge);
        nElectrons += charge;
      }
      if (hit->FromTRphoton()) {
        hQtr->Fill(charge);
        hQXtr->Fill(x,charge);
      }

    }

    hNstep->Fill(((AliTRDsimpleMC *) gMC)->GetNStep() / 3.5);
    hNel->Fill(nElectrons / 3.5);

    // Digitize the hits
    digitizer->MakeDigits();

    // Analyze the digits
    Int_t row =  8;
    Int_t col = 48;
    AliTRDdigitsManager *digitsManager = digitizer->Digits();
    AliTRDdataArrayI    *digitsArray   = digitsManager->GetDigits(12);
    for (Int_t time = 0; time < timeMax; time++) {  

      Int_t amp = digitsArray->GetDataUnchecked(row,col,time);
      if (amp > 0) {
        hAmp->Fill((Double_t) amp);
        hAmpTime->Fill((Double_t) time, (Double_t) amp);
      }

    }
    
  }

  //_____________________________________________________________________________ 
  //
  // Plot the spectra
  //_____________________________________________________________________________ 
  //

  TCanvas *cHit = new TCanvas("cHit","Hits",50,50,800,600);
  cHit->Divide(3,2);

  cHit->cd(1);
  gPad->SetLogy();
  hQ->SetLineColor(4);
  hQ->SetXTitle("Q (a.u.)");
  hQ->SetYTitle("Entries");
  hQ->Draw();
  hQdedx->SetLineColor(3);
  hQdedx->Draw("SAME");
  hQtr->SetLineColor(2);
  hQtr->Draw("SAME"); 

  cHit->cd(2); 
  gPad->SetLogy();
  hQX->SetLineColor(4);
  hQX->SetXTitle("x (cm)");
  hQX->SetYTitle("<Q> (a.u.)");
  hQX->Draw("HIST");
  hQXdedx->SetLineColor(3);
  hQXdedx->Draw("SAMEHIST");

  cHit->cd(3); 
  gPad->SetLogy();
  hQXtr->SetLineColor(2);
  hQXtr->SetXTitle("x (cm)");
  hQXtr->SetYTitle("<Q> (a.u.)");
  hQXtr->Draw("HIST"); 

  cHit->cd(4);
  hNstep->SetLineColor(4);
  hNstep->SetXTitle("N_{step}  / cm");
  hNstep->SetYTitle("Entries");
  hNstep->Draw();

  cHit->cd(5);
  hNel->SetLineColor(4);
  hNel->SetXTitle("N_{electron} / cm");
  hNel->SetYTitle("Entries");
  hNel->Draw();

  TCanvas *cDigit = new TCanvas("cDigit","Digits",50,50,600,400);
  cDigit->Divide(2,1);

  cDigit->cd(1);
  gPad->SetLogy();
  hAmp->SetLineColor(2);
  hAmp->SetXTitle("ADC channels");
  hAmp->SetYTitle("Entries");
  hAmp->Draw();

  cDigit->cd(2);
  hAmpTime->SetLineColor(2);
  hAmpTime->SetXTitle("Timebin");
  hAmpTime->SetYTitle("<ADC channels>");
  hAmpTime->Draw("HIST");

  //_____________________________________________________________________________ 
  //
  // Save the histograms
  //_____________________________________________________________________________ 
  //

  TFile *fileOut = new TFile("simple.root","RECREATE");

  hQ->Write();
  hQdedx->Write();
  hQtr->Write();

  hQX->Write();
  hQXdedx->Write();
  hQXtr->Write();

  hNstep->Write();
  hNel->Write();

  hAmp->Write();
  hAmpTime->Write();

  fileOut->Close();

}
