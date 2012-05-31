#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "Riostream.h"

#include "AliITSgeomTGeo.h"
#include "AliITSGainSSDv2.h"
#include "AliITSBadChannelsSSDv2.h"
#include "AliITSNoiseSSDv2.h"
#include "AliITSPedestalSSDv2.h"
#include "AliITSGainSSD.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSNoiseSSD.h"
#include "AliITSPedestalSSD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#endif


//====================================================================//
void Noise(AliCDBManager * man);
void Pedestal(AliCDBManager * man);
void GainCalibration(AliCDBManager * man);
void ReadOldSSDPedestal(TObjArray *array, AliITSPedestalSSDv2 *pedestalSSD);
void ReadOldSSDNoise(TObjArray *array, AliITSNoiseSSDv2 *noiseSSD);
void ReadOldSSDBadChannels(TObjArray *array, AliITSBadChannelsSSDv2 *badChannelsSSD);
void ReadOldSSDGain(TObjArray *array, AliITSGainSSDv2 *gainSSD);
void drawNoiseDistributions(Int_t runNumber);
void drawPedestalDistributions(Int_t runNumber);
//====================================================================//

//_____________________________________________________________________//
void correlateOCDBforSSD(const char* type = "alien", Int_t runNumber = 0) {
  //This macro allows to read the pedestal values from the reference 
  //directory and correlate them with the bad channel list.
  //It also allows to correlate the noise values in the OCDB with 
  //the bad channel list.
  TString gType = type;
  
  AliCDBManager *man1 = AliCDBManager::Instance();
  
  if(gType == "alien") {
    man1->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB/");
    man1->SetSpecificStorage("ITS/Ref/PedestalSSD",
			     "alien://folder=/alice/data/2009/Reference/");
  }
  else if(gType == "local") {
    man1->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    man1->SetSpecificStorage("ITS/Ref/PedestalSSD",
			     "local://$ALICE_ROOT/");
  }
  else {
    cout<<"Allowed types: local or alien!!!"<<endl;
    abort();
  }
  
  man1->SetRun(runNumber);
    
  Pedestal(man1);
  Noise(man1);
}

//_____________________________________________________________________//
void drawNoiseDistributions(Int_t runNumber) {
  //Draws the noise distributions for each side and layer
  TString filename = "noiseDistributionsSSD."; filename += runNumber;
  filename += ".root";
  TFile *f = TFile::Open(filename.Data());
  TH1F *gHistNoisePSideLayer5 = dynamic_cast<TH1F *>(f->Get("gHistNoisePSideLayer5"));
  TH1F *gHistNoiseNSideLayer5 = dynamic_cast<TH1F *>(f->Get("gHistNoiseNSideLayer5"));
  TH1F *gHistNoisePSideLayer6 = dynamic_cast<TH1F *>(f->Get("gHistNoisePSideLayer6"));
  TH1F *gHistNoiseNSideLayer6 = dynamic_cast<TH1F *>(f->Get("gHistNoiseNSideLayer6"));

  TCanvas *c1 = new TCanvas("c1","Noise distribution (P-side, Layer 5)",
			    0,0,400,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->SetLogy();
  gHistNoisePSideLayer5->SetStats(kFALSE); 
  gHistNoisePSideLayer5->GetXaxis()->SetRangeUser(0.0,20.);
  gHistNoisePSideLayer5->Draw();
  c1->SaveAs("Noise-PSide-Layer5.eps");

  TCanvas *c2 = new TCanvas("c2","Noise distribution (N-side, Layer 5)",
			    400,0,400,400);
  c2->SetFillColor(10); c2->SetHighLightColor(10); c2->SetLogy();
  gHistNoiseNSideLayer5->SetStats(kFALSE); 
  gHistNoiseNSideLayer5->GetXaxis()->SetRangeUser(0.0,20.);
  gHistNoiseNSideLayer5->Draw();
  c2->SaveAs("Noise-NSide-Layer5.eps");

  TCanvas *c3 = new TCanvas("c3","Noise distribution (P-side, Layer 6)",
			    0,400,400,400);
  c3->SetFillColor(10); c3->SetHighLightColor(10); c3->SetLogy();
  gHistNoisePSideLayer6->SetStats(kFALSE); 
  gHistNoisePSideLayer6->GetXaxis()->SetRangeUser(0.0,20.);
  gHistNoisePSideLayer6->Draw();
  c3->SaveAs("Noise-PSide-Layer6.eps");

  TCanvas *c4 = new TCanvas("c4","Noise distribution (N-side, Layer 6)",
			    400,400,400,400);
  c4->SetFillColor(10); c4->SetHighLightColor(10); c4->SetLogy();
  gHistNoiseNSideLayer6->SetStats(kFALSE); 
  gHistNoiseNSideLayer6->GetXaxis()->SetRangeUser(0.0,20.);
  gHistNoiseNSideLayer6->Draw();
  c4->SaveAs("Noise-NSide-Layer6.eps");
}

//_____________________________________________________________________//
void drawPedestalDistributions(Int_t runNumber) {
  //Draws the pedestal distributions for each side and layer
  TString filename = "pedestalDistributionsSSD."; filename += runNumber;
  filename += ".root";
  TFile *f = TFile::Open(filename.Data());
  TH1F *gHistPedestalPSideLayer5 = dynamic_cast<TH1F *>(f->Get("gHistPedestalPSideLayer5"));
  TH1F *gHistPedestalNSideLayer5 = dynamic_cast<TH1F *>(f->Get("gHistPedestalNSideLayer5"));
  TH1F *gHistPedestalPSideLayer6 = dynamic_cast<TH1F *>(f->Get("gHistPedestalPSideLayer6"));
  TH1F *gHistPedestalNSideLayer6 = dynamic_cast<TH1F *>(f->Get("gHistPedestalNSideLayer6"));

  TCanvas *c1 = new TCanvas("c1","Pedestal distribution (P-side, Layer 5)",
			    0,0,400,400);
  c1->SetFillColor(10); c1->SetHighLightColor(10); c1->SetLogy();
  gHistPedestalPSideLayer5->SetStats(kFALSE); gHistPedestalPSideLayer5->Draw();
  c1->SaveAs("Pedestal-PSide-Layer5.eps");

  TCanvas *c2 = new TCanvas("c2","Pedestal distribution (N-side, Layer 5)",
			    400,0,400,400);
  c2->SetFillColor(10); c2->SetHighLightColor(10); c2->SetLogy();
  gHistPedestalNSideLayer5->SetStats(kFALSE); gHistPedestalNSideLayer5->Draw();
  c2->SaveAs("Pedestal-NSide-Layer5.eps");
  
  TCanvas *c3 = new TCanvas("c3","Pedestal distribution (P-side, Layer 6)",
			    0,400,400,400);
  c3->SetFillColor(10); c3->SetHighLightColor(10); c3->SetLogy();
  gHistPedestalPSideLayer6->SetStats(kFALSE); gHistPedestalPSideLayer6->Draw();
  c3->SaveAs("Pedestal-PSide-Layer6.eps");

  TCanvas *c4 = new TCanvas("c4","Pedestal distribution (N-side, Layer 6)",
			    400,400,400,400);
  c4->SetFillColor(10); c4->SetHighLightColor(10); c4->SetLogy();
  gHistPedestalNSideLayer6->SetStats(kFALSE); gHistPedestalNSideLayer6->Draw();
  c4->SaveAs("Pedestal-NSide-Layer6.eps");
}

//_____________________________________________________________________//
void Pedestal(AliCDBManager * man) {
  //Reads the noise OCDB file
  const Int_t fgkSSDMODULES = 1698;
  static const Int_t fgkDefaultNStripsSSD = 768;

  Int_t runNumber = man->GetRun();

  //=========================================================//
  //pedestal histograms
  TH1F *gHistPedestalPSideLayer5 = new TH1F("gHistPedestalPSideLayer5",
					    "Pedestal values (P-side, Layer5); ADC counts; Entries;",
					    1000,-500,500);
  TH1F *gHistPedestalNSideLayer5 = new TH1F("gHistPedestalNSideLayer5",
					    "Pedestal values (N-side, Layer5); ADC counts; Entries;",
					    1000,-500,500);
  TH1F *gHistPedestalPSideLayer6 = new TH1F("gHistPedestalPSideLayer6",
					    "Pedestal values (P-side, Layer6); ADC counts; Entries;",
					    1000,-500,500);
  TH1F *gHistPedestalNSideLayer6 = new TH1F("gHistPedestalNSideLayer6",
					    "Pedestal values (N-side, Layer6); ADC counts; Entries;",
					    1000,-500,500);

  //=========================================================//
  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  
  //=========================================================//
  AliITSBadChannelsSSDv2 *badChannelsSSD = new AliITSBadChannelsSSDv2();
  AliCDBEntry *entryBadChannelsSSD = man->Get("ITS/Calib/BadChannelsSSD");
  TObject *emptyBadChannel = (TObject *)entryBadChannelsSSD->GetObject();
  TString objectnameBadChannel = emptyBadChannel->GetName();
  if(objectnameBadChannel=="TObjArray") {
    TObjArray *badChannelsSSDOld = (TObjArray *)entryBadChannelsSSD->GetObject();
    ReadOldSSDBadChannels(badChannelsSSDOld, badChannelsSSD);
  }
  else if(objectnameBadChannel=="AliITSBadChannelsSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    badChannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
  }

  //=========================================================//
  AliITSPedestalSSDv2 *pedestalSSD = new AliITSPedestalSSDv2();
  AliCDBEntry *entryPedestalSSD = man->Get("ITS/Ref/PedestalSSD");
  TObject *emptyPedestal = (TObject *)entryPedestalSSD->GetObject();
  TString objectnamePedestal = emptyPedestal->GetName();
  if(objectnamePedestal=="TObjArray") {
    TObjArray *pedestalSSDOld = (TObjArray *)entryPedestalSSD->GetObject();
    ReadOldSSDPedestal(pedestalSSDOld, pedestalSSD);
  }
  else if(objectnamePedestal=="AliITSPedestalSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    pedestalSSD = (AliITSPedestalSSDv2 *)entryPedestalSSD->GetObject();
  }

  Double_t pedestalPSide = 0.0, pedestalNSide = 0.0;
  Int_t badChannelPSide = 0, badChannelNSide = 0;
  for (Int_t i = 0; i < fgkSSDMODULES; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule);      
    //cout<<"Pedestal for module: "<<i+500<<" - Layer: "<<fLayer<<endl;
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      badChannelPSide = (Int_t)(badChannelsSSD->GetBadChannelP(i,j));
      pedestalPSide = pedestalSSD->GetPedestalP(i,j);
      if(badChannelPSide == 0) {
	//Printf("Pedestal value: %lf",pedestal);
	if(fLayer == 5) 
	  gHistPedestalPSideLayer5->Fill(pedestalPSide);
	if(fLayer == 6) 
	  gHistPedestalPSideLayer6->Fill(pedestalPSide);
      }
      
      badChannelNSide = (Int_t)(badChannelsSSD->GetBadChannelN(i,j));
      pedestalNSide = pedestalSSD->GetPedestalN(i,j);
      if(badChannelNSide == 0) {
	if(fLayer == 5) 
	  gHistPedestalNSideLayer5->Fill(pedestalNSide);
	if(fLayer == 6) 
	  gHistPedestalNSideLayer6->Fill(pedestalNSide);
      }
    }//loop over strips
  }//loop over modules

  TString output = "pedestalDistributionsSSD."; output += runNumber; 
  output += ".root";
  TFile *f = TFile::Open(output.Data(),"recreate");
  gHistPedestalPSideLayer5->Write();
  gHistPedestalNSideLayer5->Write();
  gHistPedestalPSideLayer6->Write();
  gHistPedestalNSideLayer6->Write();
  f->Close();
}

//_____________________________________________________________________//
void Noise(AliCDBManager * man) {
  //Reads the noise OCDB file
  const Int_t fgkSSDMODULES = 1698;
  const Int_t fgkSSDSTRIPSPERMODULE = 1536;
  static const Int_t fgkDefaultNStripsSSD = 768;

  Int_t runNumber = man->GetRun();

  //noise histograms
  TH1F *gHistNoisePSideLayer5 = new TH1F("gHistNoisePSideLayer5",
					 "Noise values (P-side, Layer5); ADC counts; Entries;",
					 1000,0,100);
  TH1F *gHistNoiseNSideLayer5 = new TH1F("gHistNoiseNSideLayer5",
					 "Noise values (N-side, Layer5); ADC counts; Entries;",
					 1000,0,100);
  TH1F *gHistNoisePSideLayer6 = new TH1F("gHistNoisePSideLayer6",
					 "Noise values (P-side, Layer6); ADC counts; Entries;",
					 1000,0,100);
  TH1F *gHistNoiseNSideLayer6 = new TH1F("gHistNoiseNSideLayer6",
					 "Noise values (N-side, Layer6); ADC counts; Entries;",
					 1000,0,100);

  Int_t fLayer = 0,fLadder = 0, fModule = 0;
  Int_t fHistCounter = 0;
  TString fTitle;
  TObjArray *array = new TObjArray();
  TH1D *hNoiseModule[fgkSSDMODULES];   
  for(Int_t i = 500; i < fgkSSDMODULES + 500; i++) {
    AliITSgeomTGeo::GetModuleId(i,fLayer,fLadder,fModule);
    fTitle = "SSD_Noise_Layer"; fTitle += fLayer;
    fTitle += "_Ladder"; fTitle += fLadder;
    fTitle += "_Module"; fTitle += fModule;
    
    hNoiseModule[fHistCounter] = new TH1D(fTitle.Data(),fTitle.Data(),1540,0,1540);
    hNoiseModule[fHistCounter]->GetXaxis()->SetTitleColor(1);
    hNoiseModule[fHistCounter]->GetXaxis()->SetTitle("Strip number");
    hNoiseModule[fHistCounter]->GetYaxis()->SetTitle("Noise");
    array->AddLast(hNoiseModule[fHistCounter]);
    fHistCounter += 1;
  }
  
  //=========================================================//
  AliITSBadChannelsSSDv2 *badChannelsSSD = new AliITSBadChannelsSSDv2();
  AliCDBEntry *entryBadChannelsSSD = man->Get("ITS/Calib/BadChannelsSSD");
  TObject *emptyBadChannel = (TObject *)entryBadChannelsSSD->GetObject();
  TString objectnameBadChannel = emptyBadChannel->GetName();
  if(objectnameBadChannel=="TObjArray") {
    TObjArray *badChannelsSSDOld = (TObjArray *)entryBadChannelsSSD->GetObject();
    ReadOldSSDBadChannels(badChannelsSSDOld, badChannelsSSD);
  }
  else if(objectnameBadChannel=="AliITSBadChannelsSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    badChannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
  }

  //=========================================================//
  AliITSNoiseSSDv2 *noiseSSD = new AliITSNoiseSSDv2();
  AliCDBEntry *entryNoiseSSD = man->Get("ITS/Calib/NoiseSSD");
  TObject *emptyNoise = (TObject *)entryNoiseSSD->GetObject();
  TString objectnameNoise = emptyNoise->GetName();
  if(objectnameNoise=="TObjArray") {
    TObjArray *noiseSSDOld = (TObjArray *)entryNoiseSSD->GetObject();
    ReadOldSSDNoise(noiseSSDOld, noiseSSD);
  }
  else if(objectnameNoise=="AliITSNoiseSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    noiseSSD = (AliITSNoiseSSDv2 *)entryNoiseSSD->GetObject();
  }
  
  Double_t noisePSide = 0.0, noiseNSide = 0.0;
  Int_t badChannelPSide = 0, badChannelNSide = 0;
  for (Int_t i = 0; i < fgkSSDMODULES; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,fLayer,fLadder,fModule);      
    //cout<<"Noise for module: "<<i+500<<" - Layer: "<<fLayer<<endl;
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      badChannelPSide = (Int_t)(badChannelsSSD->GetBadChannelP(i,j));
      noisePSide = noiseSSD->GetNoiseP(i,j);
      hNoiseModule[i]->SetBinContent(j+1,noisePSide);
      
      badChannelNSide = (Int_t)(badChannelsSSD->GetBadChannelN(i,j));
      noiseNSide = noiseSSD->GetNoiseN(i,j);
      hNoiseModule[i]->SetBinContent(fgkSSDSTRIPSPERMODULE-j,noiseNSide);

      if(badChannelPSide == 0) {
	if(fLayer == 5) 
	  gHistNoisePSideLayer5->Fill(noisePSide);
	if(fLayer == 6) 
	  gHistNoisePSideLayer6->Fill(noisePSide);
      }
      if(badChannelNSide == 0) {
	if(fLayer == 5) 
	  gHistNoiseNSideLayer5->Fill(noiseNSide);
	if(fLayer == 6) 
	  gHistNoiseNSideLayer6->Fill(noiseNSide);
      }
    }//loop over strips
  }//loop over modules
  
  TString output1 = "noiseSSD."; output1 += runNumber; output1 += ".root";
  TFile *f1 = TFile::Open(output1.Data(),"recreate");
  array->Write();
  f1->Close();

  TString output2 = "noiseDistributionsSSD."; output2 += runNumber; 
  output2 += ".root";
  TFile *f2 = TFile::Open(output2.Data(),"recreate");
  gHistNoisePSideLayer5->Write();
  gHistNoiseNSideLayer5->Write();
  gHistNoisePSideLayer6->Write();
  gHistNoiseNSideLayer6->Write();
  f2->Close();
}

//_____________________________________________________________________//
void GainCalibration(AliCDBManager * man) {
  const Int_t fgkSSDMODULES = 1698;
  static const Int_t fgkDefaultNStripsSSD = 768;

  TH2F *fHistGainMapLayer5 = new TH2F("fHistGainMapLayer5",
				      "Layer 5;N_{strip};N_{module}",
				      1537,0,1537,
				      750,499,1249);
  fHistGainMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistGainMapLayer5->SetStats(kFALSE);
  fHistGainMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  TH2F *fHistGainMapLayer6 = new TH2F("fHistGainMapLayer6",
				      "Layer 6;N_{strip};N_{module}",
				      1537,0,1537,
				      952,1249,2199);
  fHistGainMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistGainMapLayer6->SetStats(kFALSE);
  fHistGainMapLayer6->GetYaxis()->SetTitleOffset(1.8);

  AliITSGainSSDv2 *gainSSD = new AliITSGainSSDv2();
  AliCDBEntry *entryGainSSD = man->Get("ITS/Calib/GainSSD");
  TObject *empty = (TObject *)entryGainSSD->GetObject();
  TString objectname = empty->GetName();
  if(objectname=="Gain") {
    TObjArray *gainSSDOld = (TObjArray *)entryGainSSD->GetObject();
    ReadOldSSDGain(gainSSDOld, gainSSD);
  }
  else if(objectname=="AliITSGainSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    gainSSD = (AliITSGainSSDv2 *)entryGainSSD->GetObject();
  }
  
  Int_t layer = 0, ladder = 0, module = 0;
  Double_t gain = 0;
  for (Int_t i = 0; i < fgkSSDMODULES; i++) {
  //for (Int_t i = 0; i < 1; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      gain = gainSSD->GetGainP(i,j);
      //cout<<"GainP: "<<gain<<endl;
      if(layer == 5)
	fHistGainMapLayer5->Fill(j,i+500,gain);
      if(layer == 6)
	fHistGainMapLayer6->Fill(j,i+500,gain);
      
      gain = gainSSD->GetGainN(i,j);
      //cout<<"GainN: "<<gain<<endl;
      if(layer == 5)
	fHistGainMapLayer5->Fill(fgkDefaultNStripsSSD+j,i+500,gain);
      if(layer == 6)
	fHistGainMapLayer6->Fill(fgkDefaultNStripsSSD+j,i+500,gain);
    }//strip loop
  }//module loop

  TCanvas *cGain = new TCanvas("cGain","Gain calibration map",0,300,600,300);
  cGain->SetHighLightColor(10); cGain->SetFillColor(10); cGain->Divide(2,1);
  
  cGain->cd(1)->SetBottomMargin(.2); cGain->cd(1)->SetLeftMargin(.15);
  cGain->cd(1); fHistGainMapLayer5->Draw("colz");
  cGain->cd(2)->SetBottomMargin(.2); cGain->cd(2)->SetLeftMargin(.15);
  cGain->cd(2); fHistGainMapLayer6->Draw("colz");
}

//_____________________________________________________________________//
void ReadOldSSDNoise(TObjArray *array, 
		     AliITSNoiseSSDv2 *noiseSSD) {
  const Int_t fgkSSDSTRIPSPERMODULE = 1536;
  const Int_t fgkSSDPSIDESTRIPSPERMODULE = 768;

  Int_t fNMod = array->GetEntries();
  cout<<"Converting old calibration object for noise..."<<endl;

  //NOISE
  Double_t noise = 0.0;
  for (Int_t iModule = 0; iModule < fNMod; iModule++) {
    AliITSNoiseSSD *noiseModule = (AliITSNoiseSSD*) (array->At(iModule));
    for(Int_t iStrip = 0; iStrip < fgkSSDSTRIPSPERMODULE; iStrip++) {
      noise = (iStrip < fgkSSDPSIDESTRIPSPERMODULE) ? noiseModule->GetNoiseP(iStrip) : noiseModule->GetNoiseN(1535 - iStrip);
      if(iStrip < fgkSSDPSIDESTRIPSPERMODULE)
	noiseSSD->AddNoiseP(iModule,iStrip,noise);
      if(iStrip >= fgkSSDPSIDESTRIPSPERMODULE)
	noiseSSD->AddNoiseN(iModule,1535 - iStrip,noise);
    }//loop over strips
  }//loop over modules      
}

//_____________________________________________________________________//
void ReadOldSSDPedestal(TObjArray *array, 
			AliITSPedestalSSDv2 *pedestalSSD) {
  const Int_t fgkSSDSTRIPSPERMODULE = 1536;
  const Int_t fgkSSDPSIDESTRIPSPERMODULE = 768;

  Int_t fNMod = array->GetEntries();
  cout<<"Converting old calibration object for pedestal..."<<endl;

  //PEDESTAL
  Double_t pedestal = 0.0;
  for (Int_t iModule = 0; iModule < fNMod; iModule++) {
    AliITSPedestalSSD *pedestalModule = (AliITSPedestalSSD*) (array->At(iModule));
    for(Int_t iStrip = 0; iStrip < fgkSSDSTRIPSPERMODULE; iStrip++) {
      pedestal = (iStrip < fgkSSDPSIDESTRIPSPERMODULE) ? pedestalModule->GetPedestalP(iStrip) : pedestalModule->GetPedestalN(1535 - iStrip);
      if(iStrip < fgkSSDPSIDESTRIPSPERMODULE)
	pedestalSSD->AddPedestalP(iModule,iStrip,pedestal);
      if(iStrip >= fgkSSDPSIDESTRIPSPERMODULE)
	pedestalSSD->AddPedestalN(iModule,1535 - iStrip,pedestal);
    }//loop over strips
  }//loop over modules      
}

//_____________________________________________________________________//
void ReadOldSSDBadChannels(TObjArray *array, 
			   AliITSBadChannelsSSDv2 *badChannelsSSD) {
  Int_t fNMod = array->GetEntries();
  cout<<"Converting old calibration object for bad channels..."<<endl;
  for (Int_t iModule = 0; iModule < fNMod; iModule++) {
    //for (Int_t iModule = 0; iModule < 1; iModule++) {
    AliITSBadChannelsSSD *bad = (AliITSBadChannelsSSD*) (array->At(iModule));
    TArrayI arrayPSide = bad->GetBadPChannelsList();
    for(Int_t iPCounter = 0; iPCounter < arrayPSide.GetSize(); iPCounter++) 
      badChannelsSSD->AddBadChannelP(iModule,
				     iPCounter,
				     (Char_t)arrayPSide.At(iPCounter));
        
    TArrayI arrayNSide = bad->GetBadNChannelsList();
    for(Int_t iNCounter = 0; iNCounter < arrayNSide.GetSize(); iNCounter++) 
      badChannelsSSD->AddBadChannelN(iModule,
				     iNCounter,
				     (Char_t)arrayNSide.At(iNCounter));
    
  }//loop over modules      
}

//_____________________________________________________________________//
void ReadOldSSDGain(TObjArray *array, 
		    AliITSGainSSDv2 *gainSSD) {
  Int_t fNMod = array->GetEntries();
  cout<<"Converting old calibration object for gain..."<<endl;

  //GAIN
  for (Int_t iModule = 0; iModule < fNMod; iModule++) {
    AliITSGainSSD *gainModule = (AliITSGainSSD*) (array->At(iModule));
    TArrayF arrayPSide = gainModule->GetGainP();
    for(Int_t iPCounter = 0; iPCounter < arrayPSide.GetSize(); iPCounter++)
      gainSSD->AddGainP(iModule,
			iPCounter,
			arrayPSide.At(iPCounter));
    TArrayF arrayNSide = gainModule->GetGainN();
    for(Int_t iNCounter = 0; iNCounter < arrayNSide.GetSize(); iNCounter++)
      gainSSD->AddGainN(iModule,
			iNCounter,
			arrayNSide.At(iNCounter));
  }//loop over modules 
}
