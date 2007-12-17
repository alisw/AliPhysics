/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
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

/*-------------------------------------------------------------------------------
/* 14/12/07 New version: MUONTRKda.cxx,v 1.10
/*-------------------------------------------------------------------------------

Version for MUONTRKda MUON tracking
(A. Baldisseri, J.-L. Charvet & Ch. Finck)


Rem:  AliMUON2DMap stores all channels, even those which are not connected
      if pedMean == -1, channel not connected to a pad  


*/
extern "C" {
#include <daqDA.h>
}

#include "event.h"
#include "monitor.h"

#include <Riostream.h>
#include <stdio.h>
#include <stdlib.h>

//AliRoot
#include "AliMUONLogger.h"
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
#include "TSystem.h"
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

// global variables
const Int_t gkNChannels = AliMpConstants::ManuNofChannels();
const Int_t gkADCMax    = 4095;

AliMUONVStore* gPedestalStore =  new AliMUON2DMap(kFALSE);

Int_t  gNManu       = 0;
Int_t  gNChannel    = 0;
UInt_t gRunNumber   = 0;
Int_t  gNEvents     = 0;
Int_t  gNDateEvents = 0;
Int_t  gPrintLevel  = 1;  // global printout variable (others: 2 and 3)
Int_t  gPlotLevel  = 0;  // global plot variable

TH1F*  gPedMeanHisto  = 0x0;
TH1F*  gPedSigmaHisto = 0x0;
Char_t gHistoFileName[256];

// used by makegain 
Char_t gHistoFileName_gain[256]="MUONTRKda_gain_data.root";
Char_t gRootFileName[256];
Char_t gOutFolder[256]=".";
Char_t filename[256];
Char_t filenam[256]="MUONTRKda_gain"; 
Char_t flatFile[256];

ofstream filcout;

TString flatOutputFile;
TString logOutputFile;
TString gCommand("ped");
TTimeStamp date;

// funtions


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

    if (gNEvents == 0) {
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
//   TTimeStamp date;
  Double_t pedMean;
  Double_t pedSigma;
  ofstream fileout;
  Int_t busPatchId;
  Int_t manuId;
  Int_t channelId;

 // histo
  TFile*  histoFile = 0;
  TTree* tree = 0;
  if (gCommand.CompareTo("ped") == 0)
    {
      sprintf(gHistoFileName,"%s/MUONTRKda_ped_%d.root",gOutFolder,gRunNumber);
      histoFile = new TFile(gHistoFileName,"RECREATE","MUON Tracking pedestals");

      Char_t name[255];
      Char_t title[255];
      sprintf(name,"pedmean_allch");
      sprintf(title,"Pedestal mean all channels");
      Int_t nx = 4096;
      Int_t xmin = 0;
      Int_t xmax = 4095; 
      gPedMeanHisto = new TH1F(name,title,nx,xmin,xmax);
      gPedMeanHisto->SetDirectory(histoFile);

      sprintf(name,"pedsigma_allch");
      sprintf(title,"Pedestal sigma all channels");
      nx = 201;
      xmin = 0;
      xmax = 200; 
      gPedSigmaHisto = new TH1F(name,title,nx,xmin,xmax);
      gPedSigmaHisto->SetDirectory(histoFile);
    
      tree = new TTree("t","Pedestal tree");
      tree->Branch("bp",&busPatchId,"bp/I");
      tree->Branch("manu",&manuId,",manu/I");
      tree->Branch("channel",&channelId,",channel/I");
      tree->Branch("pedMean",&pedMean,",pedMean/D");
      tree->Branch("pedSigma",&pedSigma,",pedSigma/D");
    }

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

	if (gCommand.CompareTo("ped") == 0)
	  {
	    gPedMeanHisto->Fill(pedMean);
	    gPedSigmaHisto->Fill(pedSigma);

	    tree->Fill();
	  }
      }
    }
  }

  // file outputs
  if (!flatOutputFile.IsNull())  fileout.close();

  if (gCommand.CompareTo("ped") == 0)
    {
      histoFile->Write();
      histoFile->Close();
    }

//  delete tree;

}

//________________
void MakePedStoreForGain(Int_t injCharge)
{
    // store pedestal map in root file

//     Int_t injCharge = 200;

    TTree* tree = 0x0;

    // compute and store pedestals
    sprintf(flatFile,"%s/%s_%d_DAC_%d.ped",gOutFolder,filenam,gRunNumber,injCharge);
    cout << "\nMUONTRKda : Flat file  generated             : " << flatFile << "\n";
    MakePedStore(flatFile);
    
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
// void MakeGainStore(TString flatOutputFile)
void MakeGainStore()
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


    TFile*  histoFile = new TFile(gHistoFileName_gain);

    AliMUON2DMap* map[11];
    AliMUONVCalibParam* ped[11];
    AliMpIntPair* run[11];

    //read back from root file
    TTree* tree = (TTree*)histoFile->Get("t");
    Int_t nEntries = tree->GetEntries();

    // read back info
    for (Int_t i = 0; i < nEntries; ++i) {
      map[i] = 0x0;
      run[i] = 0x0;
      tree->SetBranchAddress("ped",&map[i]);
      tree->SetBranchAddress("run",&run[i]);
      tree->GetEvent(i);
//       std::cout << map[i] << " " << run[i] << std::endl;
    }
    gRunNumber=(UInt_t)run[0]->GetFirst();

    // some print
    cout<<"\n ********  MUONTRKda for Gain computing (Run = " << gRunNumber << ")\n" << endl;
    cout<<" * Date          : " << date.AsString("l") << "\n" << endl;
    cout << " Entries = " << nEntries << " DAC values \n" << endl; 
    for (Int_t i = 0; i < nEntries; ++i) {
      cout<< " Run = " << (Double_t)run[i]->GetFirst() << "    DAC = " << (Double_t)run[i]->GetSecond() << endl;
    }
    cout << "" << endl;


    Double_t pedMean[11];
    Double_t pedSigma[11];
    Double_t injCharge[11];
    Double_t injChargeErr[11];

// full print out 

    sprintf(filename,"%s/%s_%d.log",gOutFolder,filenam,gRunNumber);
    logOutputFile=filename;

    filcout.open(logOutputFile.Data());
    filcout<<"//====================================================" << endl;
    filcout<<"//        MUONTRKda for Gain computing (Run = " << gRunNumber << ")" << endl;
    filcout<<"//====================================================" << endl;
    filcout<<"//   * Date          : " << date.AsString("l") << "\n" << endl;




    // why 2 files ? (Ch. F.)
    FILE *pfilen = 0;
    FILE *pfilef = 0;
    if(gPrintLevel>=2)
    {
      sprintf(filename,"%s/%s_%d.param",gOutFolder,filenam,gRunNumber);
      cout << " fit parameter file               = " << filename << "\n";
      pfilen = fopen (filename,"w");

      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//  BP MANU CH. a0     a1       a2      xlim  P(chi2) P(chi2)2    Q\n");
      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilen,"//===================================================================\n");

      sprintf(filename,"%s/%s_%d.bad",gOutFolder,filenam,gRunNumber);
      cout << " Bad channel file                 = " << filename << "\n";
      pfilef = fopen (filename,"w");

      fprintf(pfilef,"//=================================================\n");
      fprintf(pfilef,"//  Bad Channel file calculated by MUONTRKda \n");
      fprintf(pfilef,"//=================================================\n");
      fprintf(pfilef,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilef,"//   * Date          : %s \n",date.AsString("l"));
      fprintf(pfilef,"//=======================================\n");
      fprintf(pfilef,"// BP MANU CH.   a1      a2     thres. Q\n");
      fprintf(pfilef,"//=======================================\n");
    }

    FILE *pfilew=0;
    if(flatOutputFile.IsNull())
      {
	sprintf(filename,"%s_%d.par",filenam,gRunNumber);
	flatOutputFile=filename;
      }
    if(!flatOutputFile.IsNull())
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

    FILE *pfilep = 0;
    if(gPrintLevel==3)
    {
      sprintf(filename,"%s/%s_%d.peak",gOutFolder,filenam,gRunNumber);
      cout << " File containing Peak mean values = " << filename << "\n";
      pfilep = fopen (filename,"w");

      fprintf(pfilep,"//===============================================================================================================================\n");
      fprintf(pfilep,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilep,"//===============================================================================================================================\n");
      fprintf(pfilep,"// BP  MANU  CH.    Ped.     <0>      <1>      <2>      <3>      <4>      <5>      <6>      <7>      <8>      <9>     <10> \n"); 
//       fprintf(pfilep,"//                      %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f  fC\n",level_fC[0],level_fC[1],level_fC[2],level_fC[3],level_fC[4],level_fC[5],level_fC[6],level_fC[7],level_fC[8],level_fC[9],level_fC[10]);
//       fprintf(pfilep,"//                      %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f\n",level_err[0],level_err[1],level_err[2],level_err[3],level_err[4],level_err[5],level_err[6],level_err[7],level_err[8],level_err[9],level_err[10]);
      fprintf(pfilep,"//===============================================================================================================================\n");
    }



//  plot out 

    TFile* gainFile = 0x0;
    sprintf(gRootFileName,"%s/%s_%d.root",gOutFolder,filenam,gRunNumber);
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
    Int_t p1 ;
    Int_t p2 ;
    Double_t gain; 
    Double_t capa=0.2; // internal capacitor (pF)

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
    tg->Branch("p1",&p1, "p1/I");
    tg->Branch("p2",&p2, "p2/I");
    tg->Branch("gain",&gain, "gain/D");

// bad BusPatch and manu
    Int_t num_tot_BP=800;
    Int_t num_tot_Manu=1500; 
//     Int_t bad_channel[num_tot_BP][num_tot_Manu];
    Int_t bad_channel[800][1500];
    for ( Int_t i = 0; i < num_tot_BP ; i++ ) 
      { for ( Int_t j = 0; j < num_tot_Manu ; j++ )   bad_channel[i][j]=0;}

      
    char graphName[256];

    // iterates over the first pedestal run
    TIter next(map[0]->CreateIterator());
    AliMUONVCalibParam* p;

    Int_t    nmanu         = 0;
    Int_t    nGoodChannel   = 0;
    Int_t    nGoodChannel_a1   = 0;
    Int_t    nBadChannel   = 0;
    Int_t    nBadChannel_a1   = 0;
    Int_t    nBadChannel_a2   = 0;
    Int_t    nplot=0;
    Double_t sumProbChi2   = 0.;
    Double_t sumA1         = 0.;
    Double_t sumProbChi2P2 = 0.;
    Double_t sumA2         = 0.;

    Double_t x[11], xErr[11], y[11], yErr[11];
    Double_t xp[11], xpErr[11], yp[11], ypErr[11];

    while ( ( p = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
    {
      ped[0]  = p;

      busPatchId = p->ID0();
      manuId     = p->ID1();

      // read back pedestal from the other runs for the given (bupatch, manu)
      for (Int_t i = 1; i < nEntries; ++i) {
	ped[i] = static_cast<AliMUONVCalibParam*>(map[i]->FindObject(busPatchId, manuId));
      }

      // compute for each channel the gain parameters
      for ( channelId = 0; channelId < ped[0]->Size() ; ++channelId ) {

	gain=0.4;

	Int_t n = 0;
	for (Int_t i = 0; i < nEntries; ++i) {

	  if (!ped[i]) continue; //shouldn't happen.
	  pedMean[i]      = ped[i]->ValueAsDouble(channelId, 0);
	  pedSigma[i]     = ped[i]->ValueAsDouble(channelId, 1);
	  injCharge[i]    = (Double_t)run[i]->GetSecond();
	  injChargeErr[i] = 0.01*injCharge[i];
	  if(injChargeErr[i] <= 1.) injChargeErr[i]=1.;

	  if (pedMean[i] < 0) continue; // not connected

	  if (pedSigma[i] <= 0) pedSigma[i] = 1.; // should not happen.
	  n++;
	}


	// print_peak_mean_values
	if(gPrintLevel==3)
	  {

	    fprintf(pfilep,"%4i%5i%5i%10.3f",busPatchId,manuId,channelId,0.);
	    fprintf(pfilep,"%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f  mV\n",pedMean[0],pedMean[1],pedMean[2],pedMean[3],pedMean[4],pedMean[5],pedMean[6],pedMean[7],pedMean[8],pedMean[9],pedMean[10]);
	    fprintf(pfilep,"                        %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f \n",pedSigma[0],pedSigma[1],pedSigma[2],pedSigma[3],pedSigma[4],pedSigma[5],pedSigma[6],pedSigma[7],pedSigma[8],pedSigma[9],pedSigma[10]);
	  }

	// makegain 


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

	TF1 *f1 = new TF1("f1",funcLin,0.,gkADCMax,2);
	TGraphErrors *graphErr = new TGraphErrors(nbpf1, x, y, xErr, yErr);

	f1->SetParameters(0,0);

	graphErr->Fit("f1","RQ");

	chi2 = f1->GetChisquare();
	f1->GetParameters(par);

	delete graphErr;
	graphErr=0;
	delete f1;

	prChi2 = TMath::Prob(chi2, nbpf1 - 2);

	Double_t xLim = pedMean[nInit + nbpf1 - 1];
	Double_t yLim = par[0]+par[1] * xLim;

	a0 = par[0];
	a1 = par[1];

	// 2. - Translation : new origin (xLim, yLim) + parabolic fit over nbf2 points

	Int_t nbpf2 = nEntries - (nInit + nbpf1) + 1;

	if(nbpf2 > 1)
	{
	  for (Int_t j = 0; j < nbpf2; j++)
	  {
	    Int_t k  = j + (nInit + nbpf1) - 1;
	    xp[j]    = pedMean[k] - xLim;
	    xpErr[j] = pedSigma[k];

	    yp[j]    = injCharge[k] - yLim - par[1]*xp[j];
	    ypErr[j] = injChargeErr[k];

	  }

	  TF1 *f2 = new TF1("f2",funcParabolic,0.,gkADCMax,1);
	  TGraphErrors *graphErr = new TGraphErrors(nbpf2, xp, yp, xpErr, ypErr);

	  graphErr->Fit(f2,"RQ");
	  chi2P2 = f2->GetChisquare();
	  f2->GetParameters(par);

	  delete graphErr;
	  graphErr=0;
	  delete f2;

	  prChi2P2 = TMath::Prob(chi2P2, nbpf2-1);


	      // ------------- print out in log file
// 	  if (busPatchId == 6 && manuId == 116 && ( channelId >= 17 && channelId <= 20) ) 
// 	    {
// 	      filcout << " \n ********! Print_out.: BP= " << busPatchId << " Manu_Id= " << manuId 
// 			<< " Ch.= " << channelId << ":" << endl;

// 	      for (Int_t j = 0; j < nbpf1; ++j)
// 		{filcout << j << " " << x[j] << " " << xErr[j] << " " << y[j] << " " << yErr[j] << endl;}
// 	      filcout << "  a0,a1 = " << a0 << " , " << a1 << " pr_chi2 = " <<  prChi2 << endl ;

// 	      for (Int_t j = 0; j < nbpf2; ++j)
// 		{filcout << j << " " << xp[j] << " " << xpErr[j] << " " << yp[j] << " " << ypErr[j] << endl;}
// 	      filcout << "  a2 = " << par[0] << " pr_chi2_2 = " <<  prChi2P2 << endl;
	      
// 	    }
	// ------------------------------------------



	  a2 = par[0];

	  par[0] = a0;
	  par[1] = a1;
	  par[2] = a2;
	  par[3] = xLim;

// 	  delete graphErr;

	}

	// Prints

	p1 = TMath::Nint(ceil(prChi2*14))+1;
	p2 = TMath::Nint(ceil(prChi2P2*14))+1;

	Double_t x0 = -par[0]/par[1]; // value of x corresponding to à 0 fC 
	threshold = TMath::Nint(ceil(par[3]-x0)); // linear if x < threshold

	if(gPrintLevel>=2)
	{
	  fprintf(pfilen,"%4i %4i %2i",busPatchId,manuId,channelId);
	  fprintf(pfilen," %6.2f %6.4f %10.3e %4.2f  %5.3f   %x  %5.3f  %x\n",
		  par[0], par[1], par[2], par[3], prChi2, p1, prChi2P2, p2);
	}

	// some tests
 
	if(par[1]< goodA1Min ||  par[1]> goodA1Max)
	{ 
	  p1=0;
	  nBadChannel_a1++;
	  if (gPrintLevel && nBadChannel_a1 < 1) 
	  {
	    cout << " !!!!! " << nBadChannel_a1 << " !!!!!!!! Bad Calib.: BP= " << busPatchId << " Manu_Id= " << manuId << 
		" Ch.= " << channelId << ":";
	    cout << "  a1 = " << par[1] << "    out of limit : [" <<  goodA1Min << "," << goodA1Max << 
		"]" << endl;
	  }
	}

	if(par[2]< goodA2Min ||  par[2]> goodA2Max)
	{ 
	  p2=0;
	  nBadChannel_a2++;
	  if (gPrintLevel && nBadChannel_a2 < 1) 
	  {
	    cout << " !!!!! " << nBadChannel_a2 << " !!!!!!!! Bad Calib.: BP= " << busPatchId << " Manu_Id= " << manuId 
		 << " Ch.= " << channelId << ":";
	    cout << "  a2 = " << par[2] << "    out of limit : [" <<  goodA2Min << "," << goodA2Max 
		 << "]" << endl;

	    for (Int_t j = 0; j < nbpf2; ++j)
	      {cout << j << " " << x[j] << " " << xErr[j] << " " << y[j] << " " << yErr[j] << endl;}
	  }
	}

	Q  = p1*16 + p2;  // fit quality 
	if(p1==0)Q=0;  // bad linear fit <=> bad calibration
 
	if(p1>0 && p2>0) 
	  {
	    nGoodChannel++;
	    sumProbChi2P2   += prChi2P2;
	    sumA2         += par[2];
	  }
	else
	  {
	    nBadChannel++;
	    if(busPatchId < num_tot_BP  && manuId < num_tot_Manu)  bad_channel[busPatchId][manuId]++;
	    else{cout << " Warning : busPatch = " << busPatchId << " Manu = " << manuId << endl;}
	    if(gPrintLevel>=2)fprintf(pfilef,"%4i %5i %2i %7.4f %10.3e %4i %2x\n",busPatchId,manuId,channelId,par[1],par[2],threshold,Q);
	  }


	if(p1>0)
	  {
	    nGoodChannel_a1++;
	    sumProbChi2   += prChi2;
	    sumA1         += par[1];
	    gain=1./(par[1]*capa);
	  }


	tg->Fill();

	if (!flatOutputFile.IsNull()) 
	  {
	    fprintf(pfilew,"%4i %5i %2i %7.4f %10.3e %4i %2x\n",busPatchId,manuId,channelId,par[1],par[2],threshold,Q);
	  }

	// Plots

	if(gPlotLevel){
	  if(Q==0  and  nplot < 100)
// 	  if(p1>1 && p2==0  and  nplot < 100)
// 	  if(p1>1 && p2>1  and  nplot < 100)
	    {
	      nplot++;
	      TF1 *f2Calib = new TF1("f2Calib",funcCalib,0.,gkADCMax,NFITPARAMS);

	      TGraphErrors *graphErr = new TGraphErrors(nEntries,pedMean,injCharge,pedSigma,injChargeErr);

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
	      delete f2Calib;
	    }
	}
      }
      nmanu++;
      if(fmod(nmanu,100)==0)std::cout << " Nb manu = " << nmanu << std::endl;
    }

    // file outputs for gain
    if (!flatOutputFile.IsNull())  fclose(pfilew);
    if(gPrintLevel==2){ fclose(pfilen); fclose(pfilef);}
    if(gPrintLevel==3)  fclose(pfilep); 

    tg->Write();
    histoFile->Close();

    //OutPut
    if (gPrintLevel) 
    {
      filcout << "\n List of problematic BusPatch and Manu " << endl;
      filcout << " ========================================" << endl;
      filcout << "        BP       Manu        Nb Channel  " << endl ;
      filcout << " ========================================" << endl;
      for ( Int_t i = 0 ; i < num_tot_BP ; i++ )
	{ for ( Int_t j = 0 ; j < num_tot_Manu ; j++ )
	    if (bad_channel[i][j] != 0 ) filcout << "\t" << i << "\t " << j << "\t\t" << bad_channel[i][j] << endl;}
      filcout << " ========================================" << endl;


      filcout << "\n Nb of channels in raw data     = " << nmanu*64 << " (" << nmanu << " Manu)" <<  endl;
      filcout << "\n Nb of fully calibrated channel = " << nGoodChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
	   << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
      filcout << "\n Nb of Bad channel              = " << nBadChannel << endl;

      filcout << "\n Nb of Good a1 channels  = " << nGoodChannel_a1 << " (" << goodA1Min << "<a1<" << goodA1Max <<  ") " << endl;

      Double_t meanA1         = sumA1/(nGoodChannel_a1);
      Double_t meanProbChi2   = sumProbChi2/(nGoodChannel_a1);
      Double_t meanA2         = sumA2/(nGoodChannel);
      Double_t meanProbChi2P2 = sumProbChi2P2/(nGoodChannel);

      Double_t capaManu = 0.2; // pF
      filcout << "\n linear fit   : <a1> = " << meanA1 << "\t  <gain>  = " <<  1./(meanA1*capaManu) 
	   << " mV/fC (capa= " << capaManu << " pF)" << endl;
      filcout <<   "        Prob(chi2)>  = " <<  meanProbChi2 << endl;
      filcout << "\n parabolic fit: <a2> = " << meanA2  << endl;
      filcout <<   "        Prob(chi2)>  = " <<  meanProbChi2P2 << "\n" << endl;

    }  


    return  ;

}

//*************************************************************//

// main routine
int main(Int_t argc, Char_t **argv) 
{
  
    // needed for streamer application
    gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
					  "*",
					  "TStreamerInfo",
					  "RIO",
					  "TStreamerInfo()"); 

    TFitter *minuitFit = new TFitter(NFITPARAMS);
    TVirtualFitter::SetFitter(minuitFit);

    Int_t skipEvents = 0;
    Int_t maxEvents  = 1000000;
    Int_t MaxDateEvents  = 1000000;
    Int_t injCharge = 0;
    Double_t nSigma = 3;
    Int_t threshold = -1;
    Char_t inputFile[256];

    Int_t gGlitchErrors= 0;
    Int_t gParityErrors= 0;
    Int_t gPaddingErrors= 0;

    TString fesOutputFile;
    TString crocusOutputFile;
    TString crocusConfigFile;

// option handler

   // decode the input line
  for (Int_t i = 1; i < argc; i++) // argument 0 is the executable name
  {
      Char_t* arg;
      
      arg = argv[i];
      if (arg[0] != '-') continue;
      switch (arg[1])
	{
	case 'f' : 
	  i++;
	  sprintf(inputFile,argv[i]);
	  break;
	case 'a' : 
	  i++;
	  flatOutputFile = argv[i];
	  break;
	case 'b' : 
	  i++;
	  sprintf(gOutFolder,argv[i]);
	  break;
	case 'o' : 
	  i++;
	  crocusOutputFile = argv[i];
	  break;
	case 'c' : 
	  i++;
	  crocusConfigFile = argv[i];
	  break;
	case 'e' : 
	  i++;
	  gCommand = argv[i];
	  break;
	case 'd' :
	  i++; 
	  gPrintLevel=atoi(argv[i]);
	  break;
	case 'g' :
	  i++; 
	  gPlotLevel=atoi(argv[i]);
	  break;
	case 's' :
	  i++; 
	  skipEvents=atoi(argv[i]);
	  break;
	case 'l' :
	  i++; 
	  injCharge=atoi(argv[i]); 
	  break;
	case 'm' :
	  i++; 
	  sscanf(argv[i],"%d",&MaxDateEvents);
	  break;
	case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
	case 'p' :
	  i++; 
	  sscanf(argv[i],"%lf",&nSigma);
	  break;
	case 'r' : 
	  i++;
	  sprintf(gHistoFileName_gain,argv[i]);
	  break;
	case 't' :
	  i++; 
	  sscanf(argv[i],"%d",&threshold);
	  break;
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                   (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>        (default = %s)",inputFile); 
	  printf("\n-c <Crocus config. file>  (default = %s)",crocusConfigFile.Data()); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n-a <Flat ASCII file>      (default = %s)",flatOutputFile.Data()); 
	  printf("\n-o <CROCUS cmd file>      (default = %s)",crocusOutputFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-b <output directory>     (default = %s)",gOutFolder);
	  printf("\n-d <print level>          (default = %d)",gPrintLevel);
	  printf("\n-g <plot level>           (default = %d)",gPlotLevel);
	  printf("\n-l <DAC level>            (default = %d)",injCharge);
	  printf("\n-m <max date events>      (default = %d)",MaxDateEvents);
	  printf("\n-s <skip events>          (default = %d)",skipEvents);
	  printf("\n-n <max events>           (default = %d)",maxEvents);
	  printf("\n-p <n sigmas>             (default = %f)",nSigma);
	  printf("\n-r root file data for gain(default = %s)",gHistoFileName_gain); 
	  printf("\n-t <threshold (-1 = no)>  (default = %d)",threshold);
	  printf("\n-e <execute ped/gain>     (default = %s)",gCommand.Data());
	  printf("\n-e <gain create>           make gain & create a new root file");
	  printf("\n-e <gain>                  make gain & update root file");
	  printf("\n-e <gain compute>          make gain & compute gains");

	  printf("\n\n");
	  exit(-1);
	default :
	  printf("%s : bad argument %s (please check %s -h)\n",argv[0],argv[i],argv[0]);
	  argc = 2; exit(-1); // exit if error
	} // end of switch  
    } // end of for i  

  // set gCommand to lower case
  gCommand.ToLower();


  // decoding the events
  
  Int_t status;
  void* event;

  gPedMeanHisto = 0x0;
  gPedSigmaHisto = 0x0;

  TStopwatch timers;

  timers.Start(kTRUE); 

  // once we have a configuration file in db
  // copy locally a file from daq detector config db 
  // The current detector is identified by detector code in variable
  // DATE_DETECTOR_CODE. It must be defined.
  // If environment variable DAQDA_TEST_DIR is defined, files are copied from DAQDA_TEST_DIR
  // instead of the database. The usual environment variables are not needed.
  if (!crocusConfigFile.IsNull()) {
    status = daqDA_DB_getFile("myconfig", crocusConfigFile.Data());
    if (status) {
      printf("Failed to get config file : %d\n",status);
      return -1;
    }
  }


  status = monitorSetDataSource(inputFile);
  if (status) {
    cerr << "ERROR : monitorSetDataSource status (hex) = " << hex << status
	      << " " << monitorDecodeError(status) << endl;
    return -1;
  }
  status = monitorDeclareMp("MUON Tracking monitoring");
  if (status) {
    cerr << "ERROR : monitorDeclareMp status (hex) = " << hex << status
	      << " " << monitorDecodeError(status) << endl;
    return -1;
  }

  Int_t busPatchId;
  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;
  TString key("MUONTRKda :");

  AliMUONRawStreamTracker* rawStream  = 0;

  
  if (gCommand.CompareTo("comp") != 0)
    {
      cout << "\nMUONTRKda : Reading data from file " << inputFile <<endl;

      while(1) 
	{
	  if (gNDateEvents >= MaxDateEvents) break;
	  if (gNEvents >= maxEvents) break;
	  if (gNEvents && gNEvents % 100 == 0) 	
	    cout<<"Cumulated:  DATE events = " << gNDateEvents << "   Used events = " << gNEvents << endl;

	  // check shutdown condition 
	  if (daqDA_checkShutdown()) 
	    break;

	  // Skip Events if needed
	  while (skipEvents) {
	    status = monitorGetEventDynamic(&event);
	    skipEvents--;
	  }

	  // starts reading
	  status = monitorGetEventDynamic(&event);
	  if (status < 0)  {
	    cout<<"EOF found"<<endl;
	    break;
	  }

	  // decoding rawdata headers
	  AliRawReader *rawReader = new AliRawReaderDate(event);
 
	  Int_t eventType = rawReader->GetType();
	  gRunNumber = rawReader->GetRunNumber();

	  // Output log file initialisations

	  if(gNDateEvents==0)
	    {
	      if (gCommand.CompareTo("ped") == 0){
		sprintf(flatFile,"%s/MUONTRKda_ped_%d.log",gOutFolder,gRunNumber);
		logOutputFile=flatFile;

		filcout.open(logOutputFile.Data());
		filcout<<"//=================================================" << endl;
		filcout<<"//        MUONTRKda for Pedestal run = "   << gRunNumber << endl;
		cout<<"\n ********  MUONTRKda for Pedestal run = " << gRunNumber << "\n" << endl;
	      }

	      if (gCommand.Contains("gain")){
		sprintf(flatFile,"%s/%s_%d_DAC_%d.log",gOutFolder,filenam,gRunNumber,injCharge);
		logOutputFile=flatFile;

		filcout.open(logOutputFile.Data());
		filcout<<"//=================================================" << endl;
		filcout<<"//        MUONTRKda for Gain run = " << gRunNumber << "  (DAC=" << injCharge << ")" << endl;
		cout<<"\n ********  MUONTRKda for Gain run = " << gRunNumber << "  (DAC=" << injCharge << ")\n" << endl;
	      }

	      filcout<<"//=================================================" << endl;
	      filcout<<"//   * Date          : " << date.AsString("l") << "\n" << endl;
	      cout<<" * Date          : " << date.AsString("l") << "\n" << endl;

	    }

	  gNDateEvents++;



	  if (eventType != PHYSICS_EVENT)
	    continue; // for the moment

	  // decoding MUON payload
// 	  AliMUONRawStreamTracker* rawStream  = new AliMUONRawStreamTracker(rawReader);
	  rawStream  = new AliMUONRawStreamTracker(rawReader);
          rawStream->DisableWarnings();

// 	  // loops over DDL 
// 	  rawStream->First();  // if GlitchError ? what we are doing ?
// 	  while( (status = rawStream->Next(busPatchId, manuId, channelId, charge)) ) 
// 	    {
  
// 	      if (gNEvents == 0) gNChannel++;
            
// 	      MakePed(busPatchId, (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		  
// 	    } // Next digit

//           if (!rawStream->IsErrorMessage()) {
//             gNEvents++;
//           }
          
	  // loops over DDL to find good events  (Alberto 11/12/07)
	  rawStream->First();  // if GlitchError ? what we are doing ?
	  while( (status = rawStream->Next(busPatchId, manuId, channelId, charge)) ) {
	  } // Next digit

          if (!rawStream->IsErrorMessage()) {
	    // loops over DDL to find good events
	    rawStream->First();  // if GlitchError ? what we are doing ?
	    while( (status = rawStream->Next(busPatchId, manuId, channelId, charge)) ) {
	      
	      if (gNEvents == 0)
		gNChannel++;
	      	      
	      MakePed(busPatchId, (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
	    } // Next digit
	    gNEvents++;
	  }
	  else
	    {
	      filcout<<"Event # "<<*(rawReader->GetEventId())<<" rejected"<<endl;
	    }
          if (rawStream->GetPayLoad()->GetGlitchErrors())  gGlitchErrors++;
          if (rawStream->GetPayLoad()->GetParityErrors())  gParityErrors++;
          if (rawStream->GetPayLoad()->GetPaddingErrors()) gPaddingErrors++;

          AliMUONLogger* log = rawStream->GetPayLoad()->GetErrorLogger();
          log->Print(key, filcout);

	  delete rawReader;
	  delete rawStream;

	} // while (1)
    }



    if (gCommand.CompareTo("ped") == 0)
      {
	sprintf(flatFile,"MUONTRKda_ped_%d.ped",gRunNumber);
	if(flatOutputFile.IsNull())flatOutputFile=flatFile;
	MakePedStore(flatOutputFile);
      }

  // option gain -> update root file with pedestal results
  // gain + create -> recreate root file
  // gain + comp -> update root file and compute gain parameters

    if (gCommand.Contains("gain")) 
      {
	MakePedStoreForGain(injCharge);
      }
  
    if (gCommand.Contains("comp")) 
      {

// 	if(flatOutputFile.IsNull())flatOutputFile="MUONTRKda_gain.par";
// 	MakeGainStore(flatOutputFile);
	MakeGainStore();
      }
  

  delete gPedestalStore;

  delete minuitFit;
  TVirtualFitter::SetFitter(0);

  timers.Stop();

  if (gCommand.CompareTo("comp") != 0)
    {
      cout << "\nMUONTRKda : Nb of DATE events     = "         << gNDateEvents    << endl;
      cout << "MUONTRKda : Nb of Glitch errors   = "         << gGlitchErrors  << endl;
      cout << "MUONTRKda : Nb of Parity errors   = "         << gParityErrors  << endl;
      cout << "MUONTRKda : Nb of Padding errors  = "         << gPaddingErrors << endl;
      cout << "MUONTRKda : Nb of events used     = "         << gNEvents        << endl;

      filcout << "\nMUONTRKda : Nb of DATE events     = "         << gNDateEvents    << endl;
      filcout << "MUONTRKda : Nb of Glitch errors   = "         << gGlitchErrors << endl;
      filcout << "MUONTRKda : Nb of Parity errors   = "         << gParityErrors << endl;
      filcout << "MUONTRKda : Nb of Padding errors  = "         << gPaddingErrors << endl;
      filcout << "MUONTRKda : Nb of events used     = "         << gNEvents        << endl;

    }


  cout << "\nMUONTRKda : Output logfile generated         : " << logOutputFile  << endl;

  if (gCommand.CompareTo("ped") == 0)
    {
      if (!(crocusConfigFile.IsNull()))
	cout << "MUONTRKda : CROCUS command file generated    : " << crocusOutputFile.Data() << endl;
      else
	cout << "MUONTRKda : WARNING no CROCUS command file generated" << endl;
      cout << "MUONTRKda : Pedestal Histo file              : " << gHistoFileName  << endl;
      cout << "MUONTRKda : Flat pedestal file  (to SHUTTLE) : " << flatOutputFile << endl;   
    }
  else
    {
      cout << "MUONTRKda : Data file for gain calculation   : " << gHistoFileName_gain  << endl;
    }

  if (gCommand.CompareTo("comp") == 0)
    {
      cout << "MUONTRKda : Root Histo. file generated       : " << gRootFileName  << endl;
      cout << "MUONTRKda : Flat gain file (to SHUTTLE)      : " << flatOutputFile << endl;   
    }

  

  // Store IN FES

  if (gCommand.CompareTo("comp") == 0 || gCommand.CompareTo("ped") == 0)
    {
      printf("\n *****  STORE FILE in FES ****** \n");

      // to be sure that env variable is set
//       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

      if (!flatOutputFile.IsNull()) 
	{

	  //       flatOutputFile.Prepend("./");
	if (gCommand.CompareTo("ped") == 0)
	  status = daqDA_FES_storeFile(flatOutputFile.Data(),"PEDESTALS");
        else 
          status = daqDA_FES_storeFile(flatOutputFile.Data(),"GAINS");
	  
	if (status) 
	    {
	      printf(" Failed to export file : %d\n",status);
	    }
	  else if(gPrintLevel) printf("Export file: %s\n",flatOutputFile.Data());
	}
    }

  filcout.close();
  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return status;
}