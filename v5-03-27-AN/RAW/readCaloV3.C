#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TStopwatch.h>
  #include <TStyle.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TString.h>
  #include <TCanvas.h>
  #include "AliRawReader.h"
  #include "AliCaloRawStreamV3.h"
  #include "AliLog.h"
#endif


void readCaloV3(const char *fileName = "./", const TString calo="PHOS")
{
  TH1D *hAmplHG  = new TH1D("hAmplHG" ,"HG amplitude" ,1024,-0.5,1023.5);
  TH1D *hAmplLG  = new TH1D("hAmplLG" ,"LG amplitude" ,1024,-0.5,1023.5);
  TH1D *hAmplTRU = new TH1D("hAmplTRU","TRU amplitude",1024,-0.5,1023.5);
  TH1F *hModule  = new TH1F("hModule" ,"Module number",10,-0.5,9.5);
  TH2F *hXZHG    = new TH2F("hXZHG"   ,"XZ HG cells",64,-0.5,63.5,64,-0.5,63.5);
  TH2F *hXZLG    = new TH2F("hXZLG"   ,"XZ LG cells",64,-0.5,63.5,64,-0.5,63.5);
  TH2F *hHWaddr  = new TH2F("hHWaddr","DDL is vs HW addr",216,-0.5,215.5,4096,-0.5,4095.5);

  //  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  //  AliLog::SetPrintRepetitions(kFALSE);

  AliRawReader *reader = AliRawReader::Create(fileName);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliCaloRawStreamV3 *stream = new AliCaloRawStreamV3(reader,calo);

  Int_t iev = 0;

  while (reader->NextEvent()) {
    AliInfoGeneral("",Form("Reading event %d\n",iev++));
    while (stream->NextDDL()) {
      while (stream->NextChannel()) {
// 	AliInfoGeneral("",Form("New channel: HW=%d, module=%d, row=%d, col=%d, caloflag=%d",
// 			    stream->GetHWAddress(),
// 			    stream->GetModule(),
// 			    stream->GetRow(),
// 			    stream->GetColumn(),
// 			    stream->GetCaloFlag()));
	hHWaddr->Fill(stream->GetDDLNumber(),stream->GetHWAddress());
	hModule->Fill(stream->GetModule());
	if (stream->IsHighGain())
	  hXZHG->Fill(stream->GetRow(),stream->GetColumn());
	if (stream->IsLowGain())
	  hXZLG->Fill(stream->GetRow(),stream->GetColumn());

	while (stream->NextBunch()) {
	  const UShort_t *sig = stream->GetSignals();
	  Int_t startBin = stream->GetStartTimeBin();
	  for (Int_t i = 0; i < stream->GetBunchLength(); i++) {
	    if (stream->IsHighGain())
	      hAmplHG->Fill(sig[i]);
	    if (stream->IsLowGain())
	      hAmplLG->Fill(sig[i]);
	    if (stream->IsTRUData())
	      hAmplTRU->Fill(sig[i]);
	  }
	}
      }
    }
    stream->Reset();
  }

  timer.Stop();
  timer.Print();

  gStyle->SetPalette(1);
  TCanvas *c1 = new TCanvas("c1","",0,0,900,600);
  c1->Divide(3,1);
  c1->cd(1);
  hAmplHG->Draw();
  c1->cd(2);
  hAmplLG->Draw();
  c1->cd(3);
  hAmplTRU->Draw();

  TCanvas *c2 = new TCanvas("c2","",0,0,900,600);
  c2->Divide(2,1);
  c2->cd(1);
  hModule->Draw();
  c2->cd(2);
  hHWaddr->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","",0,0,900,600);
  c3->Divide(2,1);
  c3->cd(1);
  hXZHG->Draw("colz");
  c3->cd(2);
  hXZLG->Draw("colz");
}
