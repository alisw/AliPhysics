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
#include "AliITSGainSSD.h"
#include "AliITSBadChannelsSSD.h"
#include "AliITSNoiseSSD.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#endif

/*  $Id$    */



//====================================================================//
void Noise(AliCDBManager * man, Int_t runNumber);
void BadChannelMap(AliCDBManager * man);
void GainCalibration(AliCDBManager * man);
void ReadOldSSDNoise(TObjArray *array, AliITSNoiseSSDv2 *noiseSSD);
void ReadOldSSDBadChannels(TObjArray *array, AliITSBadChannelsSSDv2 *badChannelsSSD);
void ReadOldSSDGain(TObjArray *array, AliITSGainSSDv2 *gainSSD);
//====================================================================//

//_____________________________________________________________________//
void readSSDOCDBEntry(const char* type = "alien", Int_t runNumber = 0) {
  //This macro allows to visualize the bad channels in the OCDB
  //The type can be either "local" or "alien" (where the OCDB file comes from)
  //The run nmber is the pedestal one
  gStyle->SetPalette(1,0);
  TString gType = type;
  
  AliCDBManager * man = AliCDBManager::Instance();
  
  if(gType == "alien") {
    man->SetDefaultStorage("alien://folder=/alice/data/2009/OCDB/");
  }
  else if(gType == "local") 
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  else {
    cout<<"Allowed types: local or alien!!!"<<endl;
    abort();
  }
  
  man->SetRun(runNumber);
    
  Noise(man,runNumber);
  BadChannelMap(man);
  GainCalibration(man);
}

//_____________________________________________________________________//
void Noise(AliCDBManager * man, Int_t runNumber) {
  const Int_t fgkSSDMODULES = 1698;
  const Int_t fgkSSDSTRIPSPERMODULE = 1536;
  static const Int_t fgkDefaultNStripsSSD = 768;

  //noise histograms
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
  
  AliITSNoiseSSDv2 *noiseSSD = new AliITSNoiseSSDv2();
  AliCDBEntry *entryNoiseSSD = man->Get("ITS/Calib/NoiseSSD");
  TObject *empty = (TObject *)entryNoiseSSD->GetObject();
  TString objectname = empty->GetName();
  if(objectname=="TObjArray") {
    TObjArray *noiseSSDOld = (TObjArray *)entryNoiseSSD->GetObject();
    ReadOldSSDNoise(noiseSSDOld, noiseSSD);
  }
  else if(objectname=="AliITSNoiseSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    noiseSSD = (AliITSNoiseSSDv2 *)entryNoiseSSD->GetObject();
  }
  
  Double_t noise = 0.0;
  for (Int_t i = 0; i < fgkSSDMODULES; i++) {
    //cout<<"Noise for module: "<<i+1<<endl;
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      noise = noiseSSD->GetNoiseP(i,j);
      hNoiseModule[i]->SetBinContent(j+1,noise);
      noise = noiseSSD->GetNoiseN(i,j);
      hNoiseModule[i]->SetBinContent(fgkSSDSTRIPSPERMODULE-j,noise);
    }//loop over strips
  }//loop over modules
  
  TString output = "noiseSSD."; output += runNumber; output += ".root";
  TFile *f = TFile::Open(output.Data(),"recreate");
  array->Write();
  f->Close();
}

//_____________________________________________________________________//
void BadChannelMap(AliCDBManager * man) {
  const Int_t fgkSSDMODULES = 1698;
  static const Int_t fgkDefaultNStripsSSD = 768;

  //_____________________________________________________________________________//
  TH2F *fHistPSideBadChannelMapLayer5 = new TH2F("fHistPSideBadChannelMapLayer5",
						 "Layer 5;N_{module};N_{ladder}",
						 22,1,23,
						 34,500,534);
  fHistPSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistPSideBadChannelMapLayer5->GetZaxis()->SetRangeUser(0.,100.);
  fHistPSideBadChannelMapLayer5->SetStats(kFALSE);
  fHistPSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistPSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistPSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistPSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistPSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
  TH2F *fHistNSideBadChannelMapLayer5 = new TH2F("fHistNSideBadChannelMapLayer5",
						 "Layer 5;N_{module};N_{ladder}",
						 22,1,23,
						 34,500,534);
  fHistNSideBadChannelMapLayer5->GetXaxis()->SetTitleColor(1);
  fHistNSideBadChannelMapLayer5->GetZaxis()->SetRangeUser(0.,100.);
  fHistNSideBadChannelMapLayer5->SetStats(kFALSE);
  fHistNSideBadChannelMapLayer5->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideBadChannelMapLayer5->GetXaxis()->SetNdivisions(22);
  fHistNSideBadChannelMapLayer5->GetYaxis()->SetNdivisions(34);
  fHistNSideBadChannelMapLayer5->GetXaxis()->SetLabelSize(0.03);
  fHistNSideBadChannelMapLayer5->GetYaxis()->SetLabelSize(0.03);
  fHistNSideBadChannelMapLayer5->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideBadChannelMapLayer5->GetZaxis()->SetTitle("Bad channels (n-side)[%]");

  TH2F *fHistPSideBadChannelMapLayer6 = new TH2F("fHistPSideBadChannelMapLayer6",
					    "Layer 6;N_{module};N_{ladder}",
					    25,1,26,
					    38,600,638);
  fHistPSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistPSideBadChannelMapLayer6->GetZaxis()->SetRangeUser(0.,100.);
  fHistPSideBadChannelMapLayer6->SetStats(kFALSE);
  fHistPSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistPSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistPSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistPSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistPSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistPSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistPSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (p-side)[%]");
  TH2F *fHistNSideBadChannelMapLayer6 = new TH2F("fHistNSideBadChannelMapLayer6",
					    "Layer 6;N_{module};N_{ladder}",
					    25,1,26,
					    38,600,638);
  fHistNSideBadChannelMapLayer6->GetXaxis()->SetTitleColor(1);
  fHistNSideBadChannelMapLayer6->GetZaxis()->SetRangeUser(0.,100.);
  fHistNSideBadChannelMapLayer6->SetStats(kFALSE);
  fHistNSideBadChannelMapLayer6->GetYaxis()->SetTitleOffset(1.8);
  fHistNSideBadChannelMapLayer6->GetXaxis()->SetNdivisions(25);
  fHistNSideBadChannelMapLayer6->GetYaxis()->SetNdivisions(38);
  fHistNSideBadChannelMapLayer6->GetXaxis()->SetLabelSize(0.03);
  fHistNSideBadChannelMapLayer6->GetYaxis()->SetLabelSize(0.03);
  fHistNSideBadChannelMapLayer6->GetZaxis()->SetTitleOffset(1.6);
  fHistNSideBadChannelMapLayer6->GetZaxis()->SetTitle("Bad channels (n-side)[%]");
  //_____________________________________________________________________________//
  
  //_____________________________________________________________________________//
  AliITSBadChannelsSSDv2 *badChannelsSSD = new AliITSBadChannelsSSDv2();
  AliCDBEntry *entryBadChannelsSSD = man->Get("ITS/Calib/BadChannelsSSD");
  TObject *empty = (TObject *)entryBadChannelsSSD->GetObject();
  TString objectname = empty->GetName();
  if(objectname=="TObjArray") {
    TObjArray *badChannelsSSDOld = (TObjArray *)entryBadChannelsSSD->GetObject();
    ReadOldSSDBadChannels(badChannelsSSDOld, badChannelsSSD);
  }
  else if(objectname=="AliITSBadChannelsSSDv2") {
    cout<<"Reading the new format of the calibration file..."<<endl;
    badChannelsSSD = (AliITSBadChannelsSSDv2 *)entryBadChannelsSSD->GetObject();
  }
  //_____________________________________________________________________________//

  //_____________________________________________________________________________//
  Int_t nPSideChannelsTotal = 0, nNSideChannelsTotal = 0;
  Int_t nBadPSideChannelsTotal = 0, nBadNSideChannelsTotal = 0;
  Int_t nBadPSideChannels = 0, nBadNSideChannels = 0;
  Int_t layer = 0, ladder = 0, module = 0;
  Int_t nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
  Int_t nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;
  //_____________________________________________________________________________//

  for(Int_t i = 0; i < fgkSSDMODULES; i++) {
    //for(Int_t i = 0; i < 1; i++) {
    AliITSgeomTGeo::GetModuleId(i+500,layer,ladder,module);
    nBadPSideChannels = 0, nBadNSideChannels = 0;
    nPSideChannelsLayer5 = 0, nNSideChannelsLayer5 = 0;
    nPSideChannelsLayer6 = 0, nNSideChannelsLayer6 = 0;

    Int_t badChannel = 0;
    for(Int_t j = 0; j < fgkDefaultNStripsSSD; j++) {
      badChannel = (Int_t)(badChannelsSSD->GetBadChannelP(i,j));
      //cout<<"Module: "<<i+500<< " Strip: "<<j<<" - "<<badChannel<<endl;
      if(badChannel != 0) {
	if(layer == 5)
	  nPSideChannelsLayer5 += 1;
	if(layer == 6)
	  nPSideChannelsLayer6 += 1;
	nBadPSideChannels += 1;
      }
      badChannel = (Int_t)(badChannelsSSD->GetBadChannelN(i,j));
      //cout<<"Module: "<<i+500<< " Strip: "<<fgkDefaultNStripsSSD+j+1<<" - "<<badChannel<<endl;
      if(badChannel != 0) {
	if(layer == 5)                                                    
	  nNSideChannelsLayer5 += 1;
	if(layer == 6)
	  nNSideChannelsLayer6 += 1;
	nBadNSideChannels += 1;
      }
    }
    if(layer == 5) {
      if(nPSideChannelsLayer5 > 0)
	fHistPSideBadChannelMapLayer5->Fill(module,499+ladder,
					    100.*nPSideChannelsLayer5/fgkDefaultNStripsSSD);
      else fHistPSideBadChannelMapLayer5->Fill(module,499+ladder,0.0001);
      if(nNSideChannelsLayer5 > 0)
	fHistNSideBadChannelMapLayer5->Fill(module,499+ladder,
					    100.*nNSideChannelsLayer5/fgkDefaultNStripsSSD);
      else fHistNSideBadChannelMapLayer5->Fill(module,499+ladder,0.0001);
    }//layer 5
    if(layer == 6) {
      if(nPSideChannelsLayer6 > 0) 
	fHistPSideBadChannelMapLayer6->Fill(module,599+ladder,
					    100.*nPSideChannelsLayer6/fgkDefaultNStripsSSD);
      else fHistPSideBadChannelMapLayer6->Fill(module,599+ladder,0.0001);
      if(nNSideChannelsLayer6 > 0) 
	fHistNSideBadChannelMapLayer6->Fill(module,599+ladder,
					    100.*nNSideChannelsLayer6/fgkDefaultNStripsSSD);
      else fHistNSideBadChannelMapLayer6->Fill(module,599+ladder,0.0001);
    }//layer 6
      
    nBadPSideChannelsTotal += nBadPSideChannels;
    nBadNSideChannelsTotal += nBadNSideChannels;
    nPSideChannelsTotal += fgkDefaultNStripsSSD;
    nNSideChannelsTotal += fgkDefaultNStripsSSD;
  }

  cout<<"================================="<<endl;
  cout<<"Bad p-Side channels: "<<100.*nBadPSideChannelsTotal/nPSideChannelsTotal<<endl;
  cout<<"Bad n-Side channels: "<<100.*nBadNSideChannelsTotal/nNSideChannelsTotal<<endl;
  cout<<"================================="<<endl;

  TCanvas *cBadChannel = new TCanvas("cBadChannel",
				     "Bad channel list",0,0,900,900);
  cBadChannel->SetHighLightColor(10); cBadChannel->SetFillColor(10); 
  cBadChannel->Divide(2,2);

  cBadChannel->cd(1)->SetBottomMargin(.2); 
  cBadChannel->cd(1)->SetLeftMargin(.15);
  cBadChannel->cd(1)->SetRightMargin(.2);
  cBadChannel->cd(1)->SetGridx(); cBadChannel->cd(1)->SetGridy();
  cBadChannel->cd(1); fHistPSideBadChannelMapLayer5->Draw("colz"); 
  cBadChannel->cd(2)->SetBottomMargin(.2); 
  cBadChannel->cd(2)->SetLeftMargin(.15);
  cBadChannel->cd(2)->SetRightMargin(.2);
  cBadChannel->cd(2)->SetGridx(); cBadChannel->cd(2)->SetGridy();
  cBadChannel->cd(2); fHistPSideBadChannelMapLayer6->Draw("colz");
  cBadChannel->cd(3)->SetBottomMargin(.2); 
  cBadChannel->cd(3)->SetLeftMargin(.15);
  cBadChannel->cd(3)->SetRightMargin(.2);
  cBadChannel->cd(3)->SetGridx(); cBadChannel->cd(3)->SetGridy();
  cBadChannel->cd(3); fHistNSideBadChannelMapLayer5->Draw("colz"); 
  cBadChannel->cd(4)->SetBottomMargin(.2); 
  cBadChannel->cd(4)->SetLeftMargin(.15);
  cBadChannel->cd(4)->SetRightMargin(.2);
  cBadChannel->cd(4)->SetGridx(); cBadChannel->cd(4)->SetGridy();
  cBadChannel->cd(4); fHistNSideBadChannelMapLayer6->Draw("colz");
  cBadChannel->SaveAs("Run-BadChannels.gif");
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
