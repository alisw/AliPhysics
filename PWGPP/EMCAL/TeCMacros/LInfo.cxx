#include "LInfo.h"
#include <TDatime.h>
#include <TFile.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TMap.h>
#include <TNtuple.h>

void LInfo::Print(Option_t *option) const 
{ 
  cout << "Runno: " << fRunNo << endl;
}

void LInfo::Compute()  
{ 
  char id[100];
  char title[100]; 
  const char *sideStr[] = {"A", "C"};
  const Int_t kNCol     = NCol();
  const Int_t kNRow     = NRow();
  const Int_t kNStrip   = NStrip();

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2; 
    for (Int_t igain=0; igain<2; ++igain) {
      if (!fhAmpOverMon[iSM][igain]) {
	sprintf(id, "hAmpOverMon%02d%d", iSM, igain);
	sprintf(title, "LED amplitude over LEDMON: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
	fhAmpOverMon[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5);   
	fhAmpOverMon[iSM][igain]->SetDirectory(0);
      } else 
	fhAmpOverMon[iSM][igain]->Reset();
    }
  }

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2; 
    for (Int_t igain=0; igain<2; ++igain) {
      for (Int_t col=0;col<kNCol;++col) {
	Int_t strip = col / 2;
	Int_t mbin = fhStrip[iSM][igain]->FindBin(strip);
	Double_t ledMonAmp = fhStrip[iSM][igain]->GetBinContent(mbin);
	for (Int_t row=0;row<kNRow;++row) {
	  Int_t lbin = fhLed[iSM][igain]->FindBin(col,row);
	  Double_t ledAmp = fhLed[iSM][igain]->GetBinContent(lbin);
	  Double_t weightf   = 0;
	  if (ledMonAmp!=0)
	    weightf = ledAmp/ledMonAmp;  
	  fhAmpOverMon[iSM][igain]->Fill(col, row, weightf);
	}
      }
    }
  }
}

void LInfo::CreateHistograms()  
{ 
  // book histograms
  char id[100];
  char title[100]; 
  const char *sideStr[] = {"A", "C"};
  const Int_t kNCol     = NCol();
  const Int_t kNRow     = NRow();
  const Int_t kNStrip   = NStrip();

  for (Int_t iSM=0; iSM<kNSM; ++iSM) {
    Int_t isector = iSM/2;
    Int_t iside = iSM%2; 
    for (Int_t igain=0; igain<2; ++igain) {
      sprintf(id, "hStrip%02d%d", iSM, igain);
      sprintf(title, "LEDMon Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhStrip[iSM][igain] = new TH1F(id, title, kNStrip, -0.5, kNStrip-0.5);
      fhStrip[iSM][igain]->SetDirectory(0);

      sprintf(id, "hStripCount%02d%d", iSM, igain);
      sprintf(title, "LEDMon Entries: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhStripCount[iSM][igain] = new TH1F(id, title, kNStrip, -0.5, kNStrip-0.5);
      fhStripCount[iSM][igain]->SetDirectory(0);

      sprintf(id, "hLed%02d%d", iSM, igain);
      sprintf(title, "Led Amplitude: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhLed[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5);
      fhLed[iSM][igain]->SetDirectory(0);

      sprintf(id, "hLedCount%02d%d", iSM, igain);
      sprintf(title, "Led Entries: SM %d (%1s%d) gain %d", iSM, sideStr[iside], isector, igain);
      fhLedCount[iSM][igain] = new TH2F(id, title, kNCol, -0.5, kNCol-0.5, kNRow, -0.5, kNRow - 0.5);
      fhLedCount[iSM][igain]->SetDirectory(0);

      fhAmpOverMon[iSM][igain] = 0;
    }
  }
}

void LInfo::FillStrip(Int_t mod, Int_t gain, Int_t strip, Double_t amp, Double_t rms)
{
  fhStrip[mod][gain]->Fill(strip, amp);
  fhStrip[mod][gain]->Fill(strip, amp);
  fhStripCount[mod][gain]->Fill(strip, rms);
  fhStripCount[mod][gain]->Fill(strip, rms);
}

void LInfo::FillLed(Int_t mod,Int_t gain, Int_t col, Int_t row, Double_t amp, Double_t rms)
{
  fhLed[mod][gain]->Fill(col, row, amp);  
  fhLedCount[mod][gain]->Fill(col, row, rms);
}


