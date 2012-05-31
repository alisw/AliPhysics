#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TStopwatch.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TCanvas.h>
  #include "AliRawReader.h"
  #include "AliAltroRawStreamV3.h"
  #include "AliLog.h"
#endif


void readAltroV3(const char *fileName = "./",Int_t ddlId = -1)
{
  TH1D *hsignal = new TH1D("signalall","All ADC counts",1024,-0.5,1023.5);
  TH1D *htime = new TH1D("timeall","All Time-Bins",1024,-0.5,1023.5);
  TH2F *hhwaddr = new TH2F("hhwaddr","DDL is vs HW addr",216,-0.5,215.5,4096,-0.5,4095.5);

  //  AliLog::SetGlobalLogLevel(AliLog::kFatal);
  //  AliLog::SetPrintRepetitions(kFALSE);

  AliRawReader *reader = AliRawReader::Create(fileName);
  reader->Reset();

  TStopwatch timer;
  timer.Start();

  AliAltroRawStreamV3 *stream = new AliAltroRawStreamV3(reader);
  if (ddlId >= 0)
    reader->Select("TPC",ddlId,ddlId);
  else
    reader->Select("TPC");

  //  stream->SelectRawData("TPC");

  Int_t iev = 0;

  while (reader->NextEvent()) {
    AliInfoGeneral("",Form("Reading event %d\n",iev++));
    while (stream->NextDDL()) {
      while (stream->NextChannel()) {
	while (stream->NextBunch()) {
	  const UShort_t *sig = stream->GetSignals();
	  Int_t startBin = stream->GetStartTimeBin();
	  for (Int_t i = 0; i < stream->GetBunchLength(); i++) {
	    hsignal->Fill(sig[i]);
	    htime->Fill(startBin--);
	    hhwaddr->Fill(stream->GetDDLNumber(),stream->GetHWAddress());
	  }
	}
      }
    }
    stream->Reset();
  }

  timer.Stop();
  timer.Print();

  TCanvas *c1 = new TCanvas("c1","",0,0,850,850);
  hsignal->Draw();
  TCanvas *c2 = new TCanvas("c2","",0,0,850,850);
  htime->Draw();
  TCanvas *c3 = new TCanvas("c3","",0,0,850,850);
  hhwaddr->Draw("colz");

}
