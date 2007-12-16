/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors:                                                       *
 *   Artur Szostak <artursz@iafrica.com>                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

///
///  @file   AliHLTMUONTrackerCalibratorComponent.cxx
///  @author Artur Szostak <artursz@iafrica.com>
///  @date   
///  @brief  Implementation of the AliHLTMUONTrackerCalibratorComponent class.
///

#include "AliHLTMUONTrackerCalibratorComponent.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include <cassert>

ClassImp(AliHLTMUONTrackerCalibratorComponent);


///////////////////////////////////////////////////////////////////////////////
// The code from here on was copied from MUONTRKda.cxx and addapted to the
// HLT framework.
//TODO: test that any of this actually works and clean up the code.
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
using namespace std;

#include <cstdio>
#include <cstdlib>

//AliRoot
#include "AliRawReaderMemory.h"
#include "AliMUONRawStreamTracker.h"
#include "AliMUONDspHeader.h"
#include "AliMUONBlockHeader.h"
#include "AliMUONBusStruct.h"
#include "AliMUONDDLTracker.h"
#include "AliMUONVStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"
#include "AliRawReaderDate.h"

//ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TTimeStamp.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TFitter.h"

#define  NFITPARAMS 4

namespace
{

// global variables
const Int_t gkNChannels = AliMpConstants::ManuNofChannels();
const Int_t gkADCMax    = 4095;

AliMUONVStore* gPedestalStore = NULL;

Int_t  gNManu       = 0;
Int_t  gNChannel    = 0;
UInt_t gRunNumber   = 0;
Int_t  gNEvents     = 0;
Int_t  gNDateEvents = 0;
Int_t  gPrintLevel  = 1;  // global printout variable
Int_t  gPlotLevel  = 0;  // global plot variable

TH1F*  gPedMeanHisto  = 0x0;
TH1F*  gPedSigmaHisto = 0x0;
Char_t gHistoFileName[256];

// used by makegain 
TString gHistoFileName_gain = "MuonTrkDA_data.root";
TString gRootFileName = "MuonTrkDA_gain.root";
Char_t filenam[256] = "MuonTrkDA_gain.param"; // if gPrintLevel  = 2

TString gCommand("ped");

}; // end of namespace


//________________
Double_t funcLin(Double_t *x, Double_t *par)
{
  return par[0] + par[1]*x[0];
}

//________________
Double_t funcParabolic(Double_t *x, Double_t *par)
{
  return par[0]*x[0]*x[0];
}

//________________
Double_t funcCalib(Double_t *x, Double_t *par)
{
  Double_t xLim= par[3];

  if(x[0] <= xLim) return par[0] + par[1]*x[0];

  Double_t yLim = par[0]+ par[1]*xLim;
  return yLim + par[1]*(x[0] - xLim) + par[2]*(x[0] - xLim)*(x[0] - xLim);
}

//__________
void MakePed(Int_t busPatchId, Int_t manuId, Int_t channelId, Int_t charge)
{

    AliMUONVCalibParam* ped = 
	static_cast<AliMUONVCalibParam*>(gPedestalStore->FindObject(busPatchId, manuId));

    if (!ped) {
      gNManu++;
      ped = new AliMUONCalibParamND(2, gkNChannels,busPatchId, manuId, -1.); // put default wise -1, not connected channel
      gPedestalStore->Add(ped);	
    }

    if (gNEvents == 1) {
      ped->SetValueAsDouble(channelId, 0, 0.);
      ped->SetValueAsDouble(channelId, 1, 0.);
    }

    Double_t pedMean  = ped->ValueAsDouble(channelId, 0) + charge;
    Double_t pedSigma = ped->ValueAsDouble(channelId, 1) + charge*charge;

    ped->SetValueAsDouble(channelId, 0, pedMean);
    ped->SetValueAsDouble(channelId, 1, pedSigma);

}

//________________
void MakePedStore(TString flatOutputFile = "")
{
  TTimeStamp date;
  Double_t pedMean;
  Double_t pedSigma;
  ofstream fileout;
  Int_t busPatchId;
  Int_t manuId;
  Int_t channelId;

 // histo
//   sprintf(gHistoFileName,"mutrkped-%d.root",gRunNumber);
  sprintf(gHistoFileName,"MuonTrkDA_%d_ped.root",gRunNumber);
  TFile*  histoFile = new TFile(gHistoFileName,"RECREATE","MUON Tracking pedestals");

  Char_t name[255];
  Char_t title[255];
  sprintf(name,"pedmean_allch");
  sprintf(title,"Pedestal mean all channels");
  Int_t nx = 1000;
  Int_t xmin = 0;
  Int_t xmax = 1000; 
  gPedMeanHisto = new TH1F(name,title,nx,xmin,xmax);
  gPedMeanHisto->SetDirectory(histoFile);

  sprintf(name,"pedsigma_allch");
  sprintf(title,"Pedestal sigma all channels");
  nx = 200;
  xmin = 0;
  xmax = 50; 
  gPedSigmaHisto = new TH1F(name,title,nx,xmin,xmax);
  gPedSigmaHisto->SetDirectory(histoFile);
    
  TTree* tree = new TTree("t","Pedestal tree");
  tree->Branch("bp",&busPatchId,"bp/I");
  tree->Branch("manu",&manuId,",manu/I");
  tree->Branch("channel",&channelId,",channel/I");
  tree->Branch("pedMean",&pedMean,",pedMean/D");
  tree->Branch("pedSigma",&pedSigma,",pedSigma/D");

  if (!flatOutputFile.IsNull()) {
    fileout.open(flatOutputFile.Data());
    fileout<<"//===========================================================================" << endl;
    fileout<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
    fileout<<"//===========================================================================" << endl;
    fileout<<"//       * Run           : " << gRunNumber << endl; 
    fileout<<"//       * Date          : " << date.AsString("l") <<endl;
    fileout<<"//       * Statictics    : " << gNEvents << endl;
    fileout<<"//       * # of MANUS    : " << gNManu << endl;
    fileout<<"//       * # of channels : " << gNChannel << endl;
    fileout<<"//"<<endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;
//     fileout<<"//format : BUS_PATCH MANU_ID CHANNEL MEAN SIGMA"<<endl;
    fileout<<"//      BP     MANU     CH.      MEAN    SIGMA"<<endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;

  }

  // iterator over pedestal
  TIter next(gPedestalStore->CreateIterator());
  AliMUONVCalibParam* ped;
  
  while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    busPatchId              = ped->ID0();
    manuId                  = ped->ID1();

    for (channelId = 0; channelId < ped->Size() ; ++channelId) {

      pedMean  = ped->ValueAsDouble(channelId, 0);

      if (pedMean > 0) { // connected channels

	ped->SetValueAsDouble(channelId, 0, pedMean/(Double_t)gNEvents);

	pedMean  = ped->ValueAsDouble(channelId, 0);
	pedSigma = ped->ValueAsDouble(channelId, 1);

	ped->SetValueAsDouble(channelId, 1, TMath::Sqrt(TMath::Abs(pedSigma/(Double_t)gNEvents - pedMean*pedMean)));

	pedMean  = ped->ValueAsDouble(channelId, 0) + 0.5 ;
// 	pedMean  = ped->ValueAsDouble(channelId, 0) ;
	pedSigma = ped->ValueAsDouble(channelId, 1);


	if (!flatOutputFile.IsNull()) {
	  fileout << "\t" << busPatchId << "\t" << manuId <<"\t"<< channelId << "\t"
		  << pedMean <<"\t"<< pedSigma << endl;
	}

	gPedMeanHisto->Fill(pedMean);
	gPedSigmaHisto->Fill(pedSigma);

	tree->Fill();
      }
    }
  }

  // file outputs
  if (!flatOutputFile.IsNull()) 
      fileout.close();

  histoFile->Write();
  histoFile->Close();

//  delete tree;


// not usable need an update of DATE ->5.28, but it compiles
// uncomment when DATE version ok.
// setenv DATE_FES_PATH
// setenv DATE_RUN_NUMBER
// setenv DATE_ROLE_NAME
// setenv DATE_DETECTOR_CODE

//   if (!flatOutputFile.IsNull()) {

//     flatOutputFile.Prepend("./");

//     status = daqDA_FES_storeFile(flatOutputFile.Data(),"MUONTRKda_results");
//     if (status) {
//       printf("Failed to export file : %d\n",status);
//       return -1;
//     }
//  }

}

//________________
void MakePedStoreForGain(Int_t injCharge)
{
    // store pedestal map in root file

//     Int_t injCharge = 200;

    TTree* tree = 0x0;

    // compute and store pedestals
    MakePedStore();

    // store in root file
//     sprintf(gHistoFileName,"mutrkgain.root");
    
    TString mode("UPDATE");

    if (gCommand.Contains("cre")) {
	mode = "RECREATE";
    }
    TFile* histoFile = new TFile(gHistoFileName_gain, mode.Data(), "MUON Tracking Gains");

    // second argument should be the injected charge, taken from config crocus file
    // put also info about run number could be usefull
    AliMpIntPair* pair   = new AliMpIntPair(gRunNumber, injCharge);

    if (mode.CompareTo("UPDATE") == 0) {
      tree = (TTree*)histoFile->Get("t");
      tree->SetBranchAddress("run",&pair);
      tree->SetBranchAddress("ped",&gPedestalStore);

    } else {
      tree = new TTree("t","Pedestal tree");
      tree->Branch("run", "AliMpIntPair",&pair);
      tree->Branch("ped", "AliMUON2DMap",&gPedestalStore);
      tree->SetBranchAddress("run",&pair);
      tree->SetBranchAddress("ped",&gPedestalStore);

    }

    tree->Fill();
    tree->Write("t", TObject::kOverwrite); // overwrite the tree
    histoFile->Close();

    delete pair;
}

//________________
void MakeGainStore(TString flatOutputFile)
{
    Double_t goodA1Min =  0.5;
    Double_t goodA1Max =  2.;
    Double_t goodA2Min = -0.5E-03;
    Double_t goodA2Max =  1.E-03;

    // open file mutrkgain.root
    // read again the pedestal for the calibration runs (9 runs ?)
    // need the injection charge from config file (to be done)
    // For each channel make a TGraphErrors (mean, sigma) vs injected charge
    // Fit with a polynomial fct
    // store the result in a flat file.

    Double_t pedMean[11];
    Double_t pedSigma[11];
    Double_t injCharge[11];
    Double_t injChargeErr[11];

    ofstream fileout;
    TTimeStamp date;

// full print out 

    // why 2 files ? (Ch. F.)
    FILE *pfilen = 0;
    if(gPrintLevel==2)
    {
//       Char_t filenam[256]="MuonTrkDA_gain.param";
      cout << " fit parameter file             = " << filenam << "\n";
      pfilen = fopen (filenam,"w");

      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//  BP MANU CH. a0     a1       a2      xlim  P(chi2) P(chi2)2    Q\n");
      fprintf(pfilen,"//===================================================================\n");
    }

    FILE *pfilew=0;
    if (!flatOutputFile.IsNull()) 
    {
      pfilew = fopen (flatOutputFile.Data(),"w");

      fprintf(pfilew,"//=================================================\n");
      fprintf(pfilew,"//  Calibration file calculated by MUONTRKda \n");
      fprintf(pfilew,"//=================================================\n");
      fprintf(pfilew,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilew,"//   * Date          : %s \n",date.AsString("l"));
      fprintf(pfilew,"//   * Statictics    : %d \n",gNEvents);
      fprintf(pfilew,"//   * # of MANUS    : %d \n",gNManu);
      fprintf(pfilew,"//   * # of channels : %d \n",gNChannel);
      fprintf(pfilew,"//-------------------------------------------------\n");
      fprintf(pfilew,"//=======================================\n");
      fprintf(pfilew,"// BP MANU CH.   a1      a2     thres. Q\n");
      fprintf(pfilew,"//=======================================\n");
    }



//     sprintf(gHistoFileName,"mutrkgain.root");
    TFile*  histoFile = new TFile(gHistoFileName_gain);

    AliMUON2DMap* map[11];
    AliMUONVCalibParam* ped[11];
    AliMpIntPair* run[11];

//  plot out 

    TFile* gainFile = 0x0;
//     sprintf(gRootFileName,"makegain.root");
    gainFile = new TFile(gRootFileName,"RECREATE");

    Double_t chi2    = 0.;
    Double_t chi2P2  = 0.;
    Double_t prChi2  = 0; 
    Double_t prChi2P2 =0;
    Double_t a0,a1,a2;
    Int_t busPatchId ;
    Int_t manuId     ;
    Int_t channelId ;
    Int_t threshold = 0;
    Int_t Q = 0;

    TTree *tg = new TTree("tg","TTree avec class Manu_DiMu");

    tg->Branch("bp",&busPatchId, "busPatchId/I");
    tg->Branch("manu",&manuId, "manuId/I");
    tg->Branch("channel",&channelId, "channelId/I");

    tg->Branch("a0",&a0, "a0/D");
    tg->Branch("a1",&a1, "a1/D");
    tg->Branch("a2",&a2, "a2/D");
    tg->Branch("Pchi2",&prChi2, "prChi2/D");
    tg->Branch("Pchi2_2",&prChi2P2, "prChi2P2/D");
    tg->Branch("Threshold",&threshold, "threshold/I");
    tg->Branch("Q",&Q, "Q/I");

    //read back from root file
    TTree* tree = (TTree*)histoFile->Get("t");
    Int_t nEntries = tree->GetEntries();
//     tree->Branch("a1",&a1, "a1/D");

    if(gPrintLevel) cout << " nEntries = " << nEntries << " DAC values \n" << endl; 

    // read back info
    for (Int_t i = 0; i < nEntries; ++i) {
      map[i] = 0x0;
      run[i] = 0x0;
      tree->SetBranchAddress("ped",&map[i]);
      tree->SetBranchAddress("run",&run[i]);
      tree->GetEvent(i);

//       std::cout << map[i] << " " << run[i] << std::endl;
    }
      
    // Q = f(ADC)
    TF1 *f1 = new TF1("f1",funcLin,0.,gkADCMax,2);
    TF1 *f2 = new TF1("f2",funcParabolic,0.,gkADCMax,1);

    char graphName[256];

    // iterates over the first pedestal run
    TIter next(map[0]->CreateIterator());
    AliMUONVCalibParam* p;

    Int_t    nmanu         = 0;
    Int_t    nBadChannel   = 0;
    Double_t sumProbChi2   = 0.;
    Double_t sumA1         = 0.;
    Double_t sumProbChi2P2 = 0.;
    Double_t sumA2         = 0.;

    Double_t x[11], xErr[11], y[11], yErr[11];

    while ( ( p = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
    {
      ped[0]  = p;

//       Int_t busPatchId = p->ID0();
//       Int_t manuId     = p->ID1();
      busPatchId = p->ID0();
      manuId     = p->ID1();

      // read back pedestal from the other runs for the given (bupatch, manu)
      for (Int_t i = 1; i < nEntries; ++i) {
	ped[i] = static_cast<AliMUONVCalibParam*>(map[i]->FindObject(busPatchId, manuId));
      }

      // compute for each channel the gain parameters
//       for ( Int_t channelId = 0; channelId < ped[0]->Size() ; ++channelId ) {
      for ( channelId = 0; channelId < ped[0]->Size() ; ++channelId ) {

	Int_t n = 0;
	for (Int_t i = 0; i < nEntries; ++i) {

	  if (!ped[i]) continue; //shouldn't happen.
	  pedMean[i]      = ped[i]->ValueAsDouble(channelId, 0);
	  pedSigma[i]     = ped[i]->ValueAsDouble(channelId, 1);
	  injCharge[i]    = (Double_t)run[i]->GetSecond();
	  injChargeErr[i] = 0.01*injCharge[i];
	  if(injChargeErr[i] <= 1.) injChargeErr[i]=1.;

// 	  if(n<2)cout << nEntries << " " << i << " " << injCharge[i] << endl;

// 	  cout << busPatchId << "\t" << manuId <<"\t"<< channelId << "\t" << n << " " << pedMean[i] << " " << pedSigma[i] << " " << injCharge[i] << " " << injChargeErr[i] << endl;

	  if (pedMean[i] < 0) continue; // not connected

	  if (pedSigma[i] <= 0) pedSigma[i] = 1.; // should not happen.
	  n++;
	}

	// makegain (JLC)


	// Fit Method:  Linear fit over 6 points + fit parabolic function  over 3  points) 

	// 1. - linear fit over 6 points

	Double_t par[4] = {0.,0.,0.,gkADCMax};

	Int_t nInit = 1;
	Int_t nbs   = nEntries - nInit;
	Int_t nbpf1 = 6; // linear fit over nbf1 points

	for (Int_t j = 0; j < nbs; ++j)
	{
	  Int_t k = j + nInit;
	  x[j]    = pedMean[k];
	  xErr[j] = pedSigma[k];
	  y[j]    = injCharge[k];
	  yErr[j] = injChargeErr[k];

	}

	TGraph *graphErr = new TGraphErrors(nbpf1, x, y, xErr, yErr);

	f1->SetParameters(0,0);

	graphErr->Fit("f1","RQ");

	chi2 = f1->GetChisquare();
	f1->GetParameters(par);

	prChi2 = TMath::Prob(chi2, nbpf1 - 2);

	Double_t xLim = pedMean[nInit + nbpf1 - 1];
	Double_t yLim = par[0]+par[1] * xLim;

	a0 = par[0];
	a1 = par[1];

	delete graphErr;


	// 2. - Translation : new origin (xLim, yLim) + parabolic fit over nbf2 points

	Int_t nbpf2 = nEntries - (nInit + nbpf1) + 1;

	if(nbpf2 > 1)
	{
	  for (Int_t j = 0; j < nbpf2; ++j)
	  {
	    Int_t k  = j + (nInit + nbpf1) - 1;
	    x[j]    = pedMean[k] - xLim;
	    xErr[j] = pedSigma[k];

	    y[j]    = injCharge[k] - yLim - par[1]*x[j];
	    yErr[j] = injChargeErr[k];

	  }

	  TGraph *graphErr = new TGraphErrors(nbpf2, x, y, xErr, yErr);

	  graphErr->Fit(f2,"RQ");
	  chi2P2 = f2->GetChisquare();
	  f2->GetParameters(par);

	  prChi2P2 = TMath::Prob(chi2P2, nbpf2-1);

	  a2 = par[0];

	  par[0] = a0;
	  par[1] = a1;
	  par[2] = a2;
	  par[3] = xLim;

	  delete graphErr;

	}

	// Prints

	Int_t p1 = TMath::Nint(ceil(prChi2*15));
	Int_t p2 = TMath::Nint(ceil(prChi2P2*15));
	Q  = p1*16 + p2;

	Double_t x0 = -par[0]/par[1]; // value of x corresponding to à 0 fC 
	threshold = TMath::Nint(ceil(par[3]-x0)); // linear if x < threshold

	if(gPrintLevel==2)
	{
	  fprintf(pfilen,"%4i %4i %2i",busPatchId,manuId,channelId);
	  fprintf(pfilen," %6.2f %6.4f %10.3e %4.2f  %5.3f   %5.3f    %x\n",
		  par[0], par[1], par[2], par[3], prChi2, prChi2P2, Q);
	}

	// some tests
 
	if(par[1]< goodA1Min ||  par[1]> goodA1Max )
	{ 
	  if (gPrintLevel) 
	  {
	    cout << " !!!!!!!!!!!!! Bad Calib.: BP= " << busPatchId << " Manu_Id= " << manuId << 
		" Ch.= " << channelId << ":";
	    cout << "  a1 = " << par[1] << "    out of limit : [" <<  goodA1Min << "," << goodA1Max << 
		"]" << endl;
	  }
	  Q=0;
	  nBadChannel++;
	}
	else if(par[2]< goodA2Min ||  par[2]> goodA2Max )
	{ 
	  if (gPrintLevel) 
	  {
	    cout << " !!!!!!!!!!!!! Bad Calib.: BP= " << busPatchId << " Manu_Id= " << manuId 
		 << " Ch.= " << channelId << ":";
	    cout << "  a2 = " << par[2] << "    out of limit : [" <<  goodA2Min << "," << goodA2Max 
		 << "]" << endl;
	  }
	  Q=0;
	  nBadChannel++;
	}
	else 
	{
	  sumProbChi2   += prChi2;
	  sumA1         += par[1];
	  sumProbChi2P2 += prChi2P2;
	  sumA2         += par[2];

	  if (gPrintLevel) 
	  {
	    if(!p1)
	    { 
	      cout << " ** Warning ** Bad Fit   : BP= " << busPatchId << " Manu_Id= " << manuId 
		   << "Ch.= " << channelId << ":";
	      cout << "  a1 = " << par[1] << "  P(chi2)_lin = " << prChi2  << endl ;
	    }
	    if(!p2)
	    { 
	      cout << " ** Warning ** Bad Fit   : BP= " << busPatchId << " Manu_Id= " << manuId 
		   << " Ch.= " << channelId << ":";
	    cout << "  a2 = " << par[2] << "  P_(chi2)_parab = " <<  prChi2P2  << endl;
	    }
	  }
	}

	tg->Fill();

	if (!flatOutputFile.IsNull()) 
	  {
	    fprintf(pfilew,"%4i %5i %2i %7.4f %10.3e %4i %2x\n",busPatchId,manuId,channelId,par[1],par[2],threshold,Q);
	  }

	// Plots

	if(gPlotLevel){
	  TF1 *f2Calib = new TF1("f2Calib",funcCalib,0.,gkADCMax,NFITPARAMS);

	  graphErr = new TGraphErrors(nEntries,pedMean,injCharge,pedSigma,injChargeErr);

	  sprintf(graphName,"BusPatch_%d_Manu_%d_Ch_%d",busPatchId, manuId,channelId);

	  graphErr->SetTitle(graphName);
	  graphErr->SetMarkerColor(3);
	  graphErr->SetMarkerStyle(12);
	  graphErr->Write(graphName);

	  sprintf(graphName,"f2_BusPatch_%d_Manu_%d_Ch_%d",busPatchId, manuId,channelId);
	  f2Calib->SetTitle(graphName);
	  f2Calib->SetLineColor(4);
	  f2Calib->SetParameters(par);
	  f2Calib->Write(graphName);

	  delete graphErr;
	  delete f2Calib;//LA
	}
      }
      nmanu++;
    }

    // file outputs for gain
    if (!flatOutputFile.IsNull()) 
      fileout.close();

    tg->Write();
    histoFile->Close();

    delete f1;
    delete f2;

    //OutPut
    if (gPrintLevel) 
    {
      cout << "\n Nb of Manu in raw data       = " << nmanu << " (" << nmanu*64 << " channels)" <<  endl;
      cout << "\n Nb of   calibrated channels  = " << nmanu*64 - nBadChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
	   << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
      cout << "\n Nb of UNcalibrated channels  = " << nBadChannel << " (a1 or a2 out of range)\n" << endl;

      Double_t meanA1         = sumA1/(nmanu*64 - nBadChannel);
      Double_t meanProbChi2   = sumProbChi2/(nmanu*64 - nBadChannel);
      Double_t meanA2         = sumA2/(nmanu*64 - nBadChannel);
      Double_t meanProbChi2P2 = sumProbChi2P2/(nmanu*64 - nBadChannel);

      Double_t capaManu = 0.2; // pF
      cout << "\n linear fit   : <a1> = " << meanA1 << "\t  <gain>  = " <<  1./(meanA1*capaManu) 
	   << " mV/fC (capa= " << capaManu << " pF)" << endl;
      cout <<   "        Prob(chi2)>  = " <<  meanProbChi2 << endl;
      cout << "\n parabolic fit: <a2> = " << meanA2  << endl;
      cout <<   "        Prob(chi2)>  = " <<  meanProbChi2P2 << "\n" << endl;
    }
}


AliHLTMUONTrackerCalibratorComponent::AliHLTMUONTrackerCalibratorComponent() :
	AliHLTCalibrationProcessor(),
	fFitter(NULL),
	fSkipEvents(0),
	fMaxEvents(1000000),
	fFlatOutputFile(),
	fCrocusOutputFile(),
	fCrocusConfigFile(),
	fInjCharge(0),
	fNoSigma(3),
	fThreshold(-1)
{
	/// Default contructor.
}


AliHLTMUONTrackerCalibratorComponent::~AliHLTMUONTrackerCalibratorComponent()
{
	/// Default destructor.
	
	assert( fFitter == NULL );
}


const char* AliHLTMUONTrackerCalibratorComponent::GetComponentID()
{
	/// Inherited from AliHLTComponent.
	/// Returns the component ID string for this component type.
	
	return AliHLTMUONConstants::TrackerCalibratorId();
}

void AliHLTMUONTrackerCalibratorComponent::GetInputDataTypes(
		vector<AliHLTComponentDataType>& list
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns the list of input block types expected by this component.
	
	list.clear();
	list.push_back( AliHLTMUONConstants::DDLRawDataType() );
}

AliHLTComponentDataType AliHLTMUONTrackerCalibratorComponent::GetOutputDataType()
{
	/// Inherited from AliHLTComponent.
	/// Returns the type of output block generated by this component.

	//TODO: fix.
	return AliHLTMUONConstants::DDLRawDataType();
}

void AliHLTMUONTrackerCalibratorComponent::GetOutputDataSize(
		unsigned long& constBase, double& inputMultiplier
	)
{
	/// Inherited from AliHLTComponent.
	/// Returns an estimate of the expected output data size.

	constBase = 0;
	inputMultiplier = 2.;  //TODO: is this the correct estimate.
}

AliHLTComponent* AliHLTMUONTrackerCalibratorComponent::Spawn()
{
	/// Inherited from AliHLTComponent.
	/// Creates a new instance of AliHLTMUONTrackerCalibratorComponent.

	return new AliHLTMUONTrackerCalibratorComponent();
}


Int_t AliHLTMUONTrackerCalibratorComponent::ScanArgument(int argc, const char** argv)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Parses the command line parameters.
	
	TString arg = argv[0];
	
	if (arg.CompareTo("-h") == 0)
	{
		HLTInfo("******************* usage **********************");
		HLTInfo("The available options are :");
		HLTInfo("-h help                   (this screen)");
		HLTInfo("");
		HLTInfo(" Input");
		HLTInfo("-c <Crocus config. file>  (default = %s)", fCrocusConfigFile.Data());
		HLTInfo("");
		HLTInfo(" Output");
		HLTInfo("-a <Flat ASCII file>      (default = %s)", fFlatOutputFile.Data());
		HLTInfo("-o <CROUCUS cmd file>     (default = %s)", fCrocusOutputFile.Data());
		HLTInfo("");
		HLTInfo(" Options");
		HLTInfo("-d <print level>          (default = %d)", gPrintLevel);
		HLTInfo("-g <plot level>           (default = %d)", gPlotLevel);
		HLTInfo("-l <DAC level>            (default = %d)", fInjCharge);
		HLTInfo("-s <skip events>          (default = %d)", fSkipEvents);
		HLTInfo("-n <max events>           (default = %d)", fMaxEvents);
		HLTInfo("-p <n sigmas>             (default = %f)", fNoSigma);
		HLTInfo("-r root file data for gain(default = %s)", gHistoFileName_gain.Data());
		HLTInfo("-t <threshold (-1 = no)>  (default = %d)", fThreshold);
		HLTInfo("-e <execute ped/gain>     (default = %s)", gCommand.Data());
		HLTInfo("-e <gain create>           make gain & create a new root file");
		HLTInfo("-e <gain>                  make gain & update root file");
		HLTInfo("-e <gain compute>          make gain & compute gains");
		return 0;  // Zero parameters parsed.
	}
	
	if (arg.CompareTo("-a") == 0)
	{
		if (argc < 2) return -EPROTO;
		fFlatOutputFile = argv[1];
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-o") == 0)
	{
		if (argc < 2) return -EPROTO;
		fCrocusOutputFile = argv[1];
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-c") == 0)
	{
		if (argc < 2) return -EPROTO;
		fCrocusConfigFile = argv[1];
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-e") == 0)
	{
		if (argc < 2) return -EPROTO;
		gCommand = argv[1];
		gCommand.ToLower();  // set command to lower case
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-d") == 0)
	{
		if (argc < 2) return -EPROTO;
		gPrintLevel = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-g") == 0)
	{
		if (argc < 2) return -EPROTO;
		gPlotLevel = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-s") == 0)
	{
		if (argc < 2) return -EPROTO;
		fSkipEvents = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-l") == 0)
	{
		if (argc < 2) return -EPROTO;
		fInjCharge = atoi(argv[1]);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-n") == 0)
	{
		if (argc < 2) return -EPROTO;
		sscanf(argv[1],"%d",&fMaxEvents);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-p") == 0)
	{
		if (argc < 2) return -EPROTO;
		sscanf(argv[1],"%lf",&fNoSigma);
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-r") == 0)
	{
		if (argc < 2) return -EPROTO;
		gHistoFileName_gain = argv[1];
		return 1;  // 1 parameter parsed.
	}
	if (arg.CompareTo("-t") == 0)
	{
		if (argc < 2) return -EPROTO;
		sscanf(argv[1],"%d",&fThreshold);
		return 1;  // 1 parameter parsed.
	}
	
	// Do not know what this argument is so return an error code.
	HLTError("Bad argument %s (please check with -h)", arg.Data());
	return -EINVAL;
}


Int_t AliHLTMUONTrackerCalibratorComponent::InitCalibration()
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Initialise the calibration component.
	
	TFitter* fFitter = new TFitter(NFITPARAMS);
	TVirtualFitter::SetFitter(fFitter);
	
	Int_t status;
	
	gPedestalStore = new AliMUON2DMap(kFALSE);
	
	gPedMeanHisto = 0x0;
	gPedSigmaHisto = 0x0;
	
	// once we have a configuration file in db
	// copy locally a file from daq detector config db 
	// The current detector is identified by detector code in variable
	// DATE_DETECTOR_CODE. It must be defined.
	// If environment variable DAQDA_TEST_DIR is defined, files are copied from DAQDA_TEST_DIR
	// instead of the database. The usual environment variables are not needed.
	if (!fCrocusConfigFile.IsNull()) {
		//status = daqDA_DB_getFile("myconfig", fCrocusConfigFile.Data());
		status = 0;
		if (status) {
			printf("Failed to get config file : %d\n",status);
			return -1;
		}
	}
	
	return 0;
}


Int_t AliHLTMUONTrackerCalibratorComponent::DeinitCalibration()
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Cleanup the calibration component releasing allocated memory.
	
	delete gPedestalStore;
	gPedestalStore = NULL;
	
	delete fFitter;
	fFitter = NULL;
	TVirtualFitter::SetFitter(0);
	
	return 0;
}


Int_t AliHLTMUONTrackerCalibratorComponent::ProcessCalibration(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Perform a calibration procedure on the new event.
	
	// Skip Events if needed
	if (fSkipEvents > 0)
	{
		fSkipEvents--;
		return 0;
	}
	
	// Do not process more than fMaxEvents.
	if (gNEvents >= fMaxEvents) return 0;
	if (gNEvents && gNEvents % 100 == 0)
		HLTInfo("Cumulated events %d", gNEvents);
	
	gNEvents++;
	
	gRunNumber = GetRunNo();
	
	/*
	Int_t eventType = GetRunType();
	if (eventType != 7)  // PHYSICS_EVENT - from event.h (DATE software)
	{
		HLTWarning("This event is not a physics event");
		return 0;
	}
	*/
	
	Int_t status;
	Int_t busPatchId;
	UShort_t manuId;
	UChar_t channelId;
	UShort_t charge;
	
	const AliHLTComponentBlockData* iter = NULL;
	
	// Loop over all DDL raw data input blocks and decode the event.
	iter = GetFirstInputBlock( AliHLTMUONConstants::DDLRawDataType() );
	while (iter != NULL)
	{
		// Make sure we have the correct muon tracker DDL type.
		if (not AliHLTMUONUtils::IsTrackerDDL(iter->fSpecification))
		{
			iter = GetNextInputBlock();
			continue;
		}

		// decoding rawdata headers
		AliRawReaderMemory* rawReader = new AliRawReaderMemory(
				reinterpret_cast<UChar_t*>(iter->fPtr), iter->fSize
			);
		rawReader->SetEquipmentID(AliHLTMUONUtils::SpecToEquipId(iter->fSpecification));
		
		// decoding MUON payload
		AliMUONRawStreamTracker* rawStream  = new AliMUONRawStreamTracker(rawReader);
		
		// loops over DDL 
		rawStream->First();
		while( (status = rawStream->Next(busPatchId, manuId, channelId, charge)) )
		{
			if (gNEvents == 1)
				gNChannel++;
			
			//       if (gPrintLevel) printf("manuId: %d, channelId: %d charge: %d\n", manuId, 
			// 			     channelId, charge);
			
			MakePed(busPatchId, (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		} // Next digit
		
		delete rawReader;
		delete rawStream;
	}
	
	if (gCommand.CompareTo("ped") == 0)
	{
		Char_t flatFile[256];
		sprintf(flatFile,"MuonTrkDA_%d_ped.ped",gRunNumber);
		if(fFlatOutputFile.IsNull())fFlatOutputFile=flatFile;
		MakePedStore(fFlatOutputFile);
	}
	else
	{
		if(fFlatOutputFile.IsNull())fFlatOutputFile="MuonTrkDA_gain.par";
	}
	
	// option gain -> update root file with pedestal results
	// gain + create -> recreate root file
	// gain + comp -> update root file and compute gain parameters
	
	if (gCommand.Contains("gain"))
		MakePedStoreForGain(fInjCharge);
	
	if (gCommand.Contains("comp"))
		MakeGainStore(fFlatOutputFile);
	
	if (gCommand.CompareTo("comp") != 0)
	{
		HLTInfo("MUONTRKda : Nb of events used     = %d", gNEvents);
	}
	if (gCommand.CompareTo("ped") == 0)
	{
		if (!(fCrocusConfigFile.IsNull()))
			HLTInfo("MUONTRKda : CROCUS command file generated : %s", fCrocusOutputFile.Data());
		else
			HLTInfo("MUONTRKda : WARNING no CROCUS command file generated");
		HLTInfo("MUONTRKda : Histo file generated for pedestal : %s", gHistoFileName);
	}
	else
	{
		HLTInfo("MUONTRKda : Histo file generated for gain     : %s", gHistoFileName_gain.Data());
		HLTInfo("MUONTRKda : Root file generated               : %s", gRootFileName.Data());
	}
	
	HLTInfo("MUONTRKda : Flat ASCII file generated         : %s", fFlatOutputFile.Data());
	
	//TODO:
	// PushBack data to shared memory ...
	//PushBack(.....);

	return 0;
}


Int_t AliHLTMUONTrackerCalibratorComponent::ShipDataToFXS(
		const AliHLTComponentEventData& /*evtData*/,
		AliHLTComponentTriggerData& /*trigData*/
	)
{
	/// Inherited from AliHLTCalibrationProcessor.
	/// Push the data to the FXS to ship off to the offline CDB.
	
	//TODO:
	// PushBack data to FXS ...
	//PushToFXS( ..... ) ;
	
	return 0;
}
