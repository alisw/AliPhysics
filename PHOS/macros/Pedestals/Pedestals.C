#if !defined(__CINT__) || defined(__MAKECINT__)
using namespace std;
#include "iostream"

#include <TStopwatch.h>
#include <TStyle.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TGrid.h>
#include <TMath.h>

#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliLog.h"
#include "AliCentralTrigger.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerClass.h"
#include "AliCDBManager.h"

#endif

//-----------------------------------------------------------------------------
static Bool_t              firstEvent = kTRUE;
static Int_t               runNum;
static UInt_t              period;
static UInt_t              orbitID;
static UInt_t              bcID;
static AliRawReader       *reader;
TString GetTriggerClass(ULong64_t);

//-----------------------------------------------------------------------------
void Pedestals(const TString rawFile=0, const char *selectTrigger="CPHI")
{
  // Read raw data, decode it to samples,
  // calculate pedestals from presamples, 
  // evaluate the signal amplitude as a maximum sample, 
  // and fill histograms with pedestals and amplitudes
  // This script should be compiled to speed up the data processing:
  // .L Pedestals.C++
  //___
  // Yuri Kharlov. 6 September 2007

  TStopwatch stopwatch;
  stopwatch.Start();
  
  if (rawFile.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("raw://");

  reader = AliRawReader::Create(rawFile);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,"PHOS");

  TString baseNamePed ="hPed";
  TString baseTitlePed="Ped in cell (";
  const char* sgain[3]={"LG","HG", "TRU"};

  const Int_t caloFlagMax=2,modMax=5,cellXMax=64,cellZMax=56;
  Int_t module,caloFlag,cellX,cellZ;
  TH1F *hPed[5][2][64][56];
  for (module=0; module<modMax; module++) {
    for (caloFlag=0; caloFlag<caloFlagMax; caloFlag++) {
      for (cellX=0; cellX<cellXMax; cellX++) {
	for (cellZ=0; cellZ<cellZMax; cellZ++) {
	  hPed[module][caloFlag][cellX][cellZ] = 0;
	}
      }
    }
  }
  TH1F *hPedHiMean1m2 = new TH1F("hPedHiMean1m2","Mean pedestals in module 2, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m2  = new TH1F("hPedHiRMS1m2" ,"RMS pedestals in module 2, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m2 = new TH1F("hPedLoMean1m2","Mean pedestals in module 2, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m2  = new TH1F("hPedLoRMS1m2" ,"RMS pedestals in module 2, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m2 = new TH1F("hPedTRUMean1m2","Mean pedestals in module 2, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m2  = new TH1F("hPedTRURMS1m2" ,"RMS pedestals in module 2, TRU"      ,100,0.,50.);

  TH1F *hPedHiMean1m3 = new TH1F("hPedHiMean1m3","Mean pedestals in module 3, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m3  = new TH1F("hPedHiRMS1m3" ,"RMS pedestals in module 3, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m3 = new TH1F("hPedLoMean1m3","Mean pedestals in module 3, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m3  = new TH1F("hPedLoRMS1m3" ,"RMS pedestals in module 3, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m3 = new TH1F("hPedTRUMean1m3","Mean pedestals in module 3, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m3  = new TH1F("hPedTRURMS1m3" ,"RMS pedestals in module 3, TRU"      ,100,0.,50.);

  TH1F *hPedHiMean1m4 = new TH1F("hPedHiMean1m4","Mean pedestals in module 4, high gain" ,100,0.,100.);
  TH1F *hPedHiRMS1m4  = new TH1F("hPedHiRMS1m4" ,"RMS pedestals in module 4, high gain"  ,100,0.,50.);
  TH1F *hPedLoMean1m4 = new TH1F("hPedLoMean1m4","Mean pedestals in module 4, low gain"  ,100,0.,100.);
  TH1F *hPedLoRMS1m4  = new TH1F("hPedLoRMS1m4" ,"RMS pedestals in module 4, low gain"   ,100,0.,50.);
  TH1F *hPedTRUMean1m4 = new TH1F("hPedTRUMean1m4","Mean pedestals in module 4, TRU"     ,1000,0.,1000.);
  TH1F *hPedTRURMS1m4  = new TH1F("hPedTRURMS1m4" ,"RMS pedestals in module 4, TRU"      ,100,0.,50.);

  hPedHiMean1m2->Sumw2();
  hPedHiRMS1m2 ->Sumw2();
  hPedLoMean1m2->Sumw2();
  hPedLoRMS1m2 ->Sumw2();
  hPedTRUMean1m2->Sumw2();
  hPedTRURMS1m2 ->Sumw2();
  hPedHiMean1m3->Sumw2();
  hPedHiRMS1m3 ->Sumw2();
  hPedLoMean1m3->Sumw2();
  hPedLoRMS1m3 ->Sumw2();
  hPedTRUMean1m3->Sumw2();
  hPedTRURMS1m3 ->Sumw2();
  hPedHiMean1m4->Sumw2();
  hPedHiRMS1m4 ->Sumw2();
  hPedLoMean1m4->Sumw2();
  hPedLoRMS1m4 ->Sumw2();
  hPedTRUMean1m4->Sumw2();
  hPedTRURMS1m4 ->Sumw2();

  TH2F *hPedHiMeanm2  = new TH2F("hPedHiMeanm2","Mean pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm2   = new TH2F("hPedHiRMSm2" ,"R.M.S. of pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm2   = new TH2F("hPedHiNumm2" ,"Number of pedestals in module 2, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm2  = new TH2F("hPedLoMeanm2","Mean pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm2   = new TH2F("hPedLoRMSm2" ,"R.M.S. of pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm2   = new TH2F("hPedLoNumm2" ,"Number of pedestals in module 2, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH2F *hPedHiMeanm3  = new TH2F("hPedHiMeanm3","Mean pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm3   = new TH2F("hPedHiRMSm3" ,"R.M.S. of pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm3   = new TH2F("hPedHiNumm3" ,"Number of pedestals in module 3, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm3  = new TH2F("hPedLoMeanm3","Mean pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm3   = new TH2F("hPedLoRMSm3" ,"R.M.S. of pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm3   = new TH2F("hPedLoNumm3" ,"Number of pedestals in module 3, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH2F *hPedHiMeanm4  = new TH2F("hPedHiMeanm4","Mean pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiRMSm4   = new TH2F("hPedHiRMSm4" ,"R.M.S. of pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedHiNumm4   = new TH2F("hPedHiNumm4" ,"Number of pedestals in module 4, high gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoMeanm4  = new TH2F("hPedLoMeanm4","Mean pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoRMSm4   = new TH2F("hPedLoRMSm4" ,"R.M.S. of pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);
  TH2F *hPedLoNumm4   = new TH2F("hPedLoNumm4" ,"Number of pedestals in module 4, low gain",
			       cellXMax,0.,cellXMax, cellZMax,0.,cellZMax);

  TH1I *hNBunches = new TH1I("hNBunches","Number of bunches",10,0,10);
  TH2I *hHWaddr   = new TH2I("hHWaddr","DDL is vs HW addr",216,0,216,4096,0,4096);
  TH1I *hModule   = new TH1I("hModule" ,"Module number", 5,0.,5);

  Int_t iEvent=0;
  runNum=0;

  while (reader->NextEvent()) {
    if (firstEvent) {
      firstEvent = kFALSE;
      runNum = reader->GetRunNumber();
      man = AliCDBManager::Instance();
      man ->SetRun(runNum);
    }
    ULong64_t triggerMask  = reader->GetClassMask();
    TString trclasses = GetTriggerClass(triggerMask);
    
    period  = reader->GetPeriod();
    orbitID = reader->GetOrbitID();
    bcID    = reader->GetBCID();
    iEvent++;
    if (!trclasses.Contains(selectTrigger)) continue;
    AliInfoGeneral("",Form("Reading event %d of type %d, time %d, trig.class \"%s\"",
			   iEvent,reader->GetType(), reader->GetTimestamp(), trclasses.Data()));
    while (stream->NextDDL()) {
      while (stream->NextChannel()) {
	module   =   stream->GetModule();   // [0-4]
	cellX    =   stream->GetCellX();    // [0-63]
	cellZ    =   stream->GetCellZ();    // [0-55]
	caloFlag =   stream->GetCaloFlag(); // [0-3]
 	if (caloFlag!=0 && caloFlag!=1) continue;
	if (module<0 || module>=modMax  || 
	    cellX<0  || cellX>=cellXMax || 
	    cellZ<0  || cellZ>=cellZMax) {
	  AliInfoGeneral("",Form("Wrong cell ID (m,x,z)=(%d,%d,%d)",module,cellX,cellZ));
	  break;
	}

	hHWaddr->Fill(stream->GetDDLNumber(),stream->GetHWAddress());
	hModule->Fill(module);
	if (hPed[module][caloFlag][cellX][cellZ] == 0) {
	  TString name  = baseNamePed;
	  TString title = baseTitlePed;
	  name +="_g"; name +=caloFlag;
	  name +="_m"; name +=module;
	  name +="_x"; name +=cellX;
	  name +="_z"; name +=cellZ;

	  title +=module; title +=",";
	  title +=cellX ; title +=",";
	  title +=cellZ ; title +="), ";
	  title +=sgain[caloFlag];

	  Int_t nx,xmin,xmax;
	  if (caloFlag==0 || caloFlag==1) {
	    nx=100;
	    xmin=0.;
	    xmax=100.;
	  }
	  else {
	    nx=1000;
	    xmin=0.;
	    xmax=1000.;
	  }
	  hPed[module][caloFlag][cellX][cellZ] = new TH1F(name,title,100,0.,100.);
	  hPed[module][caloFlag][cellX][cellZ]->Sumw2();
	  hPed[module][caloFlag][cellX][cellZ]->SetMarkerStyle(20);
	  hPed[module][caloFlag][cellX][cellZ]->SetOption("eph");
	}

	Int_t nBunches = 0;
	while (stream->NextBunch()) {
	  nBunches++;
	  const UShort_t *sig = stream->GetSignals();
	  Int_t sigLength = stream->GetBunchLength();
// 	  for (Int_t i = sigLength-70; i < sigLength; i++) {
	  for (Int_t i = 0; i < sigLength; i++) {
	    hPed[module][caloFlag][cellX][cellZ]->Fill(sig[i]);
	  }
	}
	hNBunches->Fill(nBunches);
      } // end of NextChannel()

    } // end of NextDDL()
  } // end of nextEvent()

  // Fill 2-dim histograms for mean, rms and n pedestals

  for (module=0; module<modMax; module++) {
    for (caloFlag=0; caloFlag<2; caloFlag++) {
      for (cellX=0; cellX<cellXMax; cellX++) {
	for (cellZ=0; cellZ<cellZMax; cellZ++) {
	  if (hPed[module][caloFlag][cellX][cellZ] != 0) {
	    if      (caloFlag == 0) {
	      if (module==2) {
		hPedLoMean1m2->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m2 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedLoMeanm2 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedLoRMSm2  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedLoNumm2  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	      else if (module==3) {
		hPedLoMean1m3->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m3 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedLoMeanm3 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedLoRMSm3  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedLoNumm3  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	      else if (module==4) {
		hPedLoMean1m4->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedLoRMS1m4 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedLoMeanm4 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedLoRMSm4  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedLoNumm4  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	    }
	    else if (caloFlag == 1) {
	      if (module==2) {
		hPedHiMean1m2->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m2 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm2 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedHiRMSm2  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedHiNumm2  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	      if (module==3) {
		hPedHiMean1m3->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m3 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm3 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedHiRMSm3  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedHiNumm3  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	      if (module==4) {
		hPedHiMean1m4->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedHiRMS1m4 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
		hPedHiMeanm4 ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetMean()    );
		hPedHiRMSm4  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetRMS()     );
		hPedHiNumm4  ->Fill( cellX, cellZ, hPed[module][caloFlag][cellX][cellZ]->GetEntries() );
	      }
	    }
	    else if (caloFlag == 2) {
	      if (module==2) {
		hPedTRUMean1m2->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m2 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	      if (module==3) {
		hPedTRUMean1m3->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m3 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	      if (module==4) {
		hPedTRUMean1m4->Fill( hPed[module][caloFlag][cellX][cellZ]->GetMean());
		hPedTRURMS1m4 ->Fill( hPed[module][caloFlag][cellX][cellZ]->GetRMS() );
	      }
	    }
	  }
	}
      }
    }
  }

  // Write existing histograms to a root file

  TString fileName = "ped";
  fileName += runNum;
  fileName += ".root";
  TFile *file = new TFile(fileName,"RECREATE");

  for (module=0; module<modMax; module++) {
    for (caloFlag=0; caloFlag<caloFlagMax; caloFlag++) {
      for (cellX=0; cellX<cellXMax; cellX++) {
	for (cellZ=0; cellZ<cellZMax; cellZ++) {
	  if (hPed[module][caloFlag][cellX][cellZ] != 0)
	    hPed[module][caloFlag][cellX][cellZ]->Write();
	}
      }
    }
  }

  hPedHiMean1m2->Write();
  hPedHiRMS1m2 ->Write();
  hPedLoMean1m2->Write();
  hPedLoRMS1m2 ->Write();
  hPedHiMeanm2 ->Write();
  hPedHiRMSm2  ->Write();
  hPedHiNumm2  ->Write();
  hPedLoMeanm2 ->Write();
  hPedLoRMSm2  ->Write();
  hPedLoNumm2  ->Write();
  hPedTRUMean1m2->Write();
  hPedTRURMS1m2 ->Write();
  
  hPedHiMean1m3->Write();
  hPedHiRMS1m3 ->Write();
  hPedLoMean1m3->Write();
  hPedLoRMS1m3 ->Write();
  hPedHiMeanm3 ->Write();
  hPedHiRMSm3  ->Write();
  hPedHiNumm3  ->Write();
  hPedLoMeanm3 ->Write();
  hPedLoRMSm3  ->Write();
  hPedLoNumm3  ->Write();
  hPedTRUMean1m3->Write();
  hPedTRURMS1m3 ->Write();

  hPedHiMean1m4->Write();
  hPedHiRMS1m4 ->Write();
  hPedLoMean1m4->Write();
  hPedLoRMS1m4 ->Write();
  hPedHiMeanm4 ->Write();
  hPedHiRMSm4  ->Write();
  hPedHiNumm4  ->Write();
  hPedLoMeanm4 ->Write();
  hPedLoRMSm4  ->Write();
  hPedLoNumm4  ->Write();
  hPedTRUMean1m4->Write();
  hPedTRURMS1m4 ->Write();

  hNBunches  ->Write();
  hHWaddr    ->Write();
  hModule    ->Write();
  
  file->Close();
  stopwatch.Print();
}
//-----------------------------------------------------------------------------

TString GetTriggerClass(ULong64_t triggerMask)
{
  // Convert a trigger mask to a trigger class
  
  AliCentralTrigger aCTP;
  TString configstr("");
  TString trclasses;
  if (!aCTP.LoadConfiguration(configstr)) { // Load CTP config from OCDB
    AliInfoGeneral("","No trigger configuration found in OCDB! The trigger configuration information will not be used!");
    return trclasses;
  }
  aCTP.SetClassMask(triggerMask);
  AliTriggerConfiguration *config = aCTP.GetConfiguration();
  const TObjArray& classesArray = config->GetClasses();
  Int_t nclasses = classesArray.GetEntriesFast();
  for( Int_t iclass=0; iclass < nclasses; iclass++ ) {
    AliTriggerClass* trclass = (AliTriggerClass*)classesArray.At(iclass);
    if (trclass && trclass->GetMask()>0) {
      Int_t trindex = TMath::Nint(TMath::Log2(trclass->GetMask()));
      reader->LoadTriggerClass(trclass->GetName(),trindex);
      if (triggerMask & (1ull << trindex)) {
	trclasses += " ";
	trclasses += trclass->GetName();
	trclasses += " ";
      }
    }
  }
  return trclasses;
}
