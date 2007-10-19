#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliRawReaderRoot.h"
#include "AliCaloRawStream.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1S.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "iostream.h"

#endif

/* $Id$ */

void PedAmp(const char *rawFile, const Int_t debug=0)
{
  // Read raw data, decode it to samples,
  // calculate pedestals from presamples, 
  // evaluate the signal amplitude as a maximum sample, 
  // and fill histograms with pedestals and amplitudes
  // This script should be compiled to speed up the data processing:
  // .L PedAmp.C++
  //___
  // Yuri Kharlov. 6 September 2007

  TStopwatch stopwatch;
  stopwatch.Start();
  
  AliRawReaderRoot* rf = new AliRawReaderRoot(rawFile);
  AliCaloRawStream in(rf,"PHOS");
  in.SetOldRCUFormat(kTRUE);
  const Int_t nPresample = 10;

  TString baseNamePed ="hPed";
  TString baseTitlePed="Ped in cell (";
  TString baseNameLed ="hLED";
  TString baseTitleLed="LED in cell (";
  const char* sgain[2]={"high","low"};

  const Int_t gainMax=2,modMax=5,colMax=56,rowMax=64;
  TH1F *hPed[2][56][64];
  TH1F *hLed[2][56][64];
  for (Int_t gain=0; gain<gainMax; gain++) {
    for (Int_t mod=0; mod<modMax; mod++) {
      for (Int_t col=0; col<colMax; col++) {
	for (Int_t row=0; row<rowMax; row++) {
	  hPed[gain][col][row] = 0;
	  hLed[gain][col][row] = 0;
	}
      }
    }
  }
  TH1F *hPedHiMean1 = new TH1F("hPedHiMean1","Mean pedestals, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1  = new TH1F("hPedHiRMS1" ,"RMS pedestals, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1 = new TH1F("hPedLoMean1","Mean pedestals, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1  = new TH1F("hPedLoRMS1" ,"RMS pedestals, low gain"   ,100,0.,50.);
  hPedHiMean1->Sumw2();
  hPedHiRMS1 ->Sumw2();
  hPedLoMean1->Sumw2();
  hPedLoRMS1 ->Sumw2();

  TH2F *hPedHiMean  = new TH2F("hPedHiMean","Mean pedestals, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hPedHiRMS   = new TH2F("hPedHiRMS" ,"R.M.S. of pedestals, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hPedHiNum   = new TH2F("hPedHiNum" ,"Number of pedestals, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hPedLoMean  = new TH2F("hPedLoMean","Mean pedestals, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hPedLoRMS   = new TH2F("hPedLoRMS" ,"R.M.S. of pedestals, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hPedLoNum   = new TH2F("hPedLoNum" ,"Number of pedestals, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  
  TH2F *hLedHiMean  = new TH2F("hLedHiMean","Mean LED, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hLedHiRMS   = new TH2F("hLedHiRMS" ,"R.M.S. of LED, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hLedHiNum   = new TH2F("hLedHiNum" ,"Number of LED, high gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hLedLoMean  = new TH2F("hLedLoMean","Mean LED, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hLedLoRMS   = new TH2F("hLedLoRMS" ,"R.M.S. of LED, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);
  TH2F *hLedLoNum   = new TH2F("hLedLoNum" ,"Number of LED, low gain",
			       rowMax,0.,rowMax, colMax,0.,colMax);

  Int_t iEvent=0;
  Int_t iBin=0;
  Int_t gain=-1,mod=-1,col=-1,row=-1;
  Int_t runNum = 0;
  Double_t maxAmp = 0;
//   Double_t meanPed,rmsPed;
  Float_t pedMean = 0;
  Int_t   nPed = 0;
  while (rf->NextEvent()) {
    cout << "\n\n\tEvent "<< iEvent << " of type " << rf->GetType() << endl;
    runNum = rf->GetRunNumber();
    while ( in.Next() ) {
      if (!in.IsNewHWAddress() && debug > 1) {
	cout << " time=" << in.GetTime() << "/"<< in.GetTimeLength()
	     << " amp=" <<  in.GetSignal() <<endl;
      }
      if (in.IsNewHWAddress()) {
	iBin   = 0;
	if (gain!=-1 && col!=-1 && row!=-1) {
	  if (nPed < 1) continue;
	  maxAmp -= (Double_t)(pedMean/nPed); // "pedestal subtraction"
	  if (maxAmp > 5) hLed[gain][col][row]->Fill(maxAmp);
	}
	maxAmp  = 0;
	pedMean = 0;
	nPed    = 0;
	gain = in.IsLowGain(); // 0-high, 1-low
	mod  = in.GetModule();
	col  = in.GetColumn();
	row  = in.GetRow();
	if (debug > 0) {
	  cout << "\n\tNew HW address: "<<in.GetHWAddress()
	       << ", gain="<< gain
	       << ", module="<<mod
	       << ", xz=("<<col<<","<<row<<")"<< endl;
	}

	// Create a new histo for the new cell, if it does not exist yet
	if (!hPed[gain][col][row]) {
	  TString name  = baseNamePed;
	  TString title = baseTitlePed;
	  name +="_g"; name +=gain;
	  name +="_m"; name +=mod;
	  name +="_z"; name +=col;
	  name +="_x"; name +=row;

	  title +=mod; title +=",";
	  title +=col; title +=",";
	  title +=row; title +="), ";
	  title +=sgain[gain]; title +=" gain";

	  hPed[gain][col][row] = new TH1F(name,title,100,0.,100.);
	  hPed[gain][col][row]->Sumw2();
	  hPed[gain][col][row]->SetMarkerStyle(20);
	  hPed[gain][col][row]->SetOption("eph");
	}
	if (!hLed[gain][col][row]) {
	  TString name  = baseNameLed;
	  TString title = baseTitleLed;
	  name +="_g"; name +=gain;
	  name +="_m"; name +=mod;
	  name +="_z"; name +=col;
	  name +="_x"; name +=row;

	  title +=mod; title +=",";
	  title +=col; title +=",";
	  title +=row; title +="), ";
	  title +=sgain[gain]; title +=" gain";

	  hLed[gain][col][row] = new TH1F(name,title,1023,0.,1023.);
	  hLed[gain][col][row]->Sumw2();
	  hLed[gain][col][row]->SetMarkerStyle(20);
	  hLed[gain][col][row]->SetOption("eph");
	}
      }

      iBin++;
      // Assume that first nPresample are pedestals
      if (in.GetTimeLength()-iBin < nPresample) {
	hPed[gain][col][row]->Fill(in.GetSignal());
	pedMean += in.GetSignal();
	nPed++;
      }
      if (in.GetSignal() > maxAmp) maxAmp = in.GetSignal();
    } // end of next signal
    iEvent++;
  } // end of next event

  // Fill 2-dim histograms for mean, rms and n pedestals

  gain = 0;
  for (Int_t col=0; col<colMax; col++) {
    for (Int_t row=0; row<rowMax; row++) {
      if (hPed[gain][col][row] != 0) {
	hPedHiMean1->Fill( hPed[gain][col][row]->GetMean());
	hPedHiRMS1 ->Fill( hPed[gain][col][row]->GetRMS() );
	hPedHiMean ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetMean()    );
	hPedHiRMS  ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetRMS()     );
	hPedHiNum  ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetEntries() );
      }
      if (hLed[gain][col][row] != 0) {
	hLedHiMean ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetMean()    );
	hLedHiRMS  ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetRMS()     );
	hLedHiNum  ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetEntries() );
      }
    }
  }
  gain = 1;
  for (Int_t col=0; col<colMax; col++) {
    for (Int_t row=0; row<rowMax; row++) {
      if (hPed[gain][col][row] != 0) {
	hPedLoMean1->Fill( hPed[gain][col][row]->GetMean());
	hPedLoRMS1 ->Fill( hPed[gain][col][row]->GetRMS() );
	hPedLoMean ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetMean()    );
	hPedLoRMS  ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetRMS()     );
	hPedLoNum  ->Fill( (Double_t)row, (Double_t)col, hPed[gain][col][row]->GetEntries() );
      }
      if (hLed[gain][col][row] != 0) {
	hLedLoMean ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetMean()    );
	hLedLoRMS  ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetRMS()     );
	hLedLoNum  ->Fill( (Double_t)row, (Double_t)col, hLed[gain][col][row]->GetEntries() );
      }
    }
  }

  // Write existing histograms to a root file

  TString fileName = "ped";
  fileName += runNum;
  fileName += ".root";
  TFile *file = new TFile(fileName,"RECREATE");

  for (Int_t gain=0; gain<gainMax; gain++) {
    for (Int_t col=0; col<colMax; col++) {
      for (Int_t row=0; row<rowMax; row++) {
	if (hPed[gain][col][row] != 0)
	  hPed[gain][col][row]->Write();
	if (hLed[gain][col][row] != 0)
	  hLed[gain][col][row]->Write();
      }
    }
  }
  hPedHiMean1->Write();
  hPedHiRMS1 ->Write();
  hPedLoMean1->Write();
  hPedLoRMS1 ->Write();
  hPedHiMean ->Write();
  hPedHiRMS  ->Write();
  hPedHiNum  ->Write();
  hPedLoMean ->Write();
  hPedLoRMS  ->Write();
  hPedLoNum  ->Write();
  
  hLedHiMean ->Write();
  hLedHiRMS  ->Write();
  hLedHiNum  ->Write();
  hLedLoMean ->Write();
  hLedLoRMS  ->Write();
  hLedLoNum  ->Write();
  
  file->Close();
  stopwatch.Print();
}
