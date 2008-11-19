/*
Contact: Jean-Luc Charvet <jean-luc.charvet@cea.fr>
	Link: http://aliceinfo.cern.ch/static/Offline/dimuon/muon_html/README_Mchda
Run Type: PEDESTAL, CALIBRATION
	DA Type: LDC
	Number of events needed: 400 events for pedestal and each calibration run
	Input Files: Rawdata file (DATE format)
	Output Files: local dir (not persistent) -> MUONTRKda_ped_<run#>.ped , MUONTRKda_gain_<run#>.par
	FXS -> run<#>_MCH_<ldc>_PEDESTALS, run<#>_MCH_<ldc>_GAINS
	Trigger types used:
*/

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

/*
	-------------------------------------------------------------------------
	2008-11-14 New version: MUONTRKda.cxx,v 1.15
	-------------------------------------------------------------------------

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
#include <math.h> 

//AliRoot
#include "AliMUONLogger.h"
#include "AliMUONRawStreamTrackerHP.h"
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
#include "AliRawDataErrorLog.h"

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
#include "THashTable.h"
#include <THashList.h>

#define  NFITPARAMS 4

// global variables
const Int_t gkNChannels = AliMpConstants::ManuNofChannels();
const Int_t gkADCMax    = 4095;

AliMUONVStore* gPedestalStore =  new AliMUON2DMap(kFALSE);

Int_t  gNManu       = 0;
Int_t  gNChannel    = 0;
UInt_t gRunNumber   = 0;
Int_t  gNEvents     = 0;
Int_t  gNEventsRecovered = 0;
Int_t  gNDateEvents = 0;
Int_t  gPrintLevel  = 1;  // global printout variable (others: 2 and 3)
Int_t  gPlotLevel  = 0;  // global plot variable
Int_t  gFES  = 1;  // by default FES is used

TH1F*  gPedMeanHisto  = 0x0;
TH1F*  gPedSigmaHisto = 0x0;
Char_t gHistoFileName[256];

// used for computing gain parameters 
Int_t nbpf1 = 6; // linear fit over nbf1 points

Char_t gHistoFileName_gain[256]="MUONTRKda_gain.data";
Char_t gRootFileName[256];
Char_t gOutFolder[256]=".";
Char_t filename[256];
Char_t filenam[256]="MUONTRKda_gain"; 
Char_t flatFile[256]="";


ofstream filcout;

TString flatOutputFile;
TString logOutputFile;
TString logOutputFile_comp;
TString gCommand("ped");
TTimeStamp date;

class ErrorCounter : public TNamed
{
public :
  ErrorCounter(Int_t bp = 0, Int_t manu = 0, Int_t ev = 1) : busPatch(bp), manuId(manu), events(ev) {}
  void Increment() {events++;}
  Int_t BusPatch() {return busPatch;}
  Int_t ManuId() {return manuId;}
  Int_t Events() {return events;}
  Int_t Compare(const TObject*) const;

  void Print(Option_t* option="") const
  {
    TNamed::Print(option);
    cout<<"bp "<<busPatch<<" events "<<events<<endl;
  }
  void Print_uncal(Option_t* option="") const
  {
    TNamed::Print(option);
    cout<<"bp =  "<<busPatch<< "  manu = " << manuId << " uncal = "<< events <<endl;
  }

private :
  Int_t busPatch; // Buspath ID
  Int_t manuId;   // Manu ID
  Int_t events;   // Events with error in this buspatch
};

Int_t ErrorCounter::Compare(const TObject* obj) const
{
  Int_t patch1, patch2, manu1, manu2;
  patch1 = busPatch;
  manu1 = manuId;
  patch2 = ((ErrorCounter*)obj)->BusPatch();
  manu2 = ((ErrorCounter*)obj)->ManuId();

  if (patch1 == patch2)
    {
      if (manu1 == manu2)
	{
	  return 0;
	}
      else
	return (manu1 >= manu2) ? 1 : -1;
    }
  else
    return (patch1 >= patch2) ? 1 : -1;
};


// Table for buspatches with parity errors 
THashTable* gErrorBuspatchTable = new THashTable(100,2);

// Table for uncalibrated  buspatches and manus
// THashTable* gUncalBuspatchManuTable = new THashTable(1000,2);
 THashList* gUncalBuspatchManuTable = new THashList(1000,2);


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

	// if (gNEvents == 0) {
	// 	ped->SetValueAsDouble(channelId, 0, 0.);
	// 	ped->SetValueAsDouble(channelId, 1, 0.);
	// }
	
	// Initialization for the first value
	if (ped->ValueAsDouble(channelId, 0) == -1) ped->SetValueAsDouble(channelId, 0, 0.);
	if (ped->ValueAsDouble(channelId, 1) == -1) ped->SetValueAsDouble(channelId, 1, 0.);

	Double_t pedMean  = ped->ValueAsDouble(channelId, 0) + (Double_t) charge;
	Double_t pedSigma = ped->ValueAsDouble(channelId, 1) + (Double_t) charge*charge;

	ped->SetValueAsDouble(channelId, 0, pedMean);
	ped->SetValueAsDouble(channelId, 1, pedSigma);

}

//________________
void MakePedStore(TString flatOutputFile_1 = "")
//void MakePedStore()
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

	if (!flatOutputFile_1.IsNull()) {
		fileout.open(flatOutputFile_1.Data());
		fileout<<"//===========================================================================" << endl;
		fileout<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
		fileout<<"//===========================================================================" << endl;
		fileout<<"//       * Run           : " << gRunNumber << endl; 
		fileout<<"//       * Date          : " << date.AsString("l") <<endl;
		fileout<<"//       * Statictics    : " << gNEvents << endl;
		fileout<<"//       * # of MANUS    : " << gNManu << endl;
		fileout<<"//       * # of channels : " << gNChannel << endl;
		if (gErrorBuspatchTable->GetSize())
		{
			fileout<<"//"<<endl;
			fileout<<"//       * Buspatches with less statistics (due to parity errors)"<<endl;		
			TIterator* iter = gErrorBuspatchTable->MakeIterator();
			ErrorCounter* parityerror;
			while((parityerror = (ErrorCounter*) iter->Next()))
			{
				fileout<<"//         bp "<<parityerror->BusPatch()<<" events used "<<gNEvents-parityerror->Events()<<endl;
			}
	
		}	
		fileout<<"//"<<endl;
		fileout<<"//---------------------------------------------------------------------------" << endl;
		fileout<<"//---------------------------------------------------------------------------" << endl;
		fileout<<"//      BP     MANU     CH.      MEAN    SIGMA"<<endl;
		fileout<<"//---------------------------------------------------------------------------" << endl;

	}
	// print in logfile
	if (gErrorBuspatchTable->GetSize())
	  {
	    cout<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;		
	    filcout<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;		
	    TIterator* iter = gErrorBuspatchTable->MakeIterator();
	    ErrorCounter* parityerror;
	    while((parityerror = (ErrorCounter*) iter->Next()))
	      {
		cout<<"  bp "<<parityerror->BusPatch()<<": events used = "<<gNEvents-parityerror->Events()<<endl;
		filcout<<"  bp "<<parityerror->BusPatch()<<": events used = "<<gNEvents-parityerror->Events()<<endl;
	      }
	
	  }

	
// iterator over pedestal
	TIter next(gPedestalStore->CreateIterator());
	AliMUONVCalibParam* ped;

	while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
	{
		busPatchId              = ped->ID0();
		manuId                  = ped->ID1();
		Int_t eventCounter;

		// Correct the number of events for buspatch with errors
		char bpname[256];
		ErrorCounter* errorCounter;
		sprintf(bpname,"bp%d",busPatchId);						
		if ((errorCounter = (ErrorCounter*)gErrorBuspatchTable->FindObject(bpname)))
		{
			eventCounter = gNEvents - errorCounter->Events();
		}
		else
		{
			eventCounter = gNEvents;
		}			

		for (channelId = 0; channelId < ped->Size() ; ++channelId) {
			pedMean  = ped->ValueAsDouble(channelId, 0);
			
			if (pedMean > 0) { // connected channels

  			ped->SetValueAsDouble(channelId, 0, pedMean/(Double_t)eventCounter);

			pedMean  = ped->ValueAsDouble(channelId, 0);
			pedSigma = ped->ValueAsDouble(channelId, 1);

			ped->SetValueAsDouble(channelId, 1, TMath::Sqrt(TMath::Abs(pedSigma/(Double_t)eventCounter - pedMean*pedMean)));

			pedMean  = ped->ValueAsDouble(channelId, 0);
			pedSigma = ped->ValueAsDouble(channelId, 1);


			if (!flatOutputFile_1.IsNull()) {
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
if (!flatOutputFile_1.IsNull())  fileout.close();

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

	FILE *pfilew=0;
	if (gCommand.Contains("gain") && !gCommand.Contains("comp")) {
		if(flatOutputFile.IsNull())
		{
			sprintf(filename,"%s_%d_DAC_%d.par",filenam,gRunNumber,injCharge);
			flatOutputFile=filename;
		}
		if(!flatOutputFile.IsNull())
		{
			pfilew = fopen (flatOutputFile.Data(),"w");

			fprintf(pfilew,"//DUMMY FILE (to prevent Shuttle failure)\n");
			fprintf(pfilew,"//================================================\n");
			fprintf(pfilew,"//       MUONTRKda: Calibration run  \n");
			fprintf(pfilew,"//=================================================\n");
			fprintf(pfilew,"//   * Run           : %d \n",gRunNumber); 
			fprintf(pfilew,"//   * Date          : %s \n",date.AsString("l"));
			fprintf(pfilew,"//   * DAC           : %d \n",injCharge);
			fprintf(pfilew,"//-------------------------------------------------\n");
			fclose(pfilew);
		}
	}

	if(gPrintLevel==2)
	{
		// compute and store pedestals
		sprintf(flatFile,"%s/%s_%d_DAC_%d.ped",gOutFolder,filenam,gRunNumber,injCharge);
		cout << "\nMUONTRKda : Flat file  generated  : " << flatFile << "\n";
		MakePedStore(flatFile);
	}
	else
		MakePedStore();

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
void MakeGainStore()
{
  ofstream filcouc;

  Int_t nInit = 1; // DAC=0 excluded from fit procedure
  Double_t goodA1Min =  0.5;
  Double_t goodA1Max =  2.;
  //     Double_t goodA1Min =  0.7;
  //     Double_t goodA1Max =  1.7;
  Double_t goodA2Min = -0.5E-03;
  Double_t goodA2Max =  1.E-03;

  Int_t num_RUN[15];

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
  Int_t nbpf2 = nEntries - (nInit + nbpf1) + 1; // nb pts used for 2nd order fit

  // read back info
  for (Int_t i = 0; i < nEntries; ++i) {
    map[i] = 0x0;
    run[i] = 0x0;
    tree->SetBranchAddress("ped",&map[i]);
    tree->SetBranchAddress("run",&run[i]);
    tree->GetEvent(i);
    //        std::cout << map[i] << " " << run[i] << std::endl;
  }
  //jlc_feb_08  modif:   gRunNumber=(UInt_t)run[0]->GetFirst();
  gRunNumber=(UInt_t)run[nEntries-1]->GetFirst();
  //     sscanf(getenv("DATE_RUN_NUMBER"),"%d",&gRunNumber);

  Double_t pedMean[11];
  Double_t pedSigma[11];
  for ( Int_t i=0 ; i<11 ; i++) {pedMean[i]=0.;pedSigma[i]=1.;};
  Double_t injCharge[11];
  Double_t injChargeErr[11];
  for ( Int_t i=0 ; i<11 ; i++) {injCharge[i]=0.;injChargeErr[i]=1.;};

  // some print
  cout<<"\n ********  MUONTRKda for Gain computing (Run = " << gRunNumber << ")\n" << endl;
  cout<<" * Date          : " << date.AsString("l") << "\n" << endl;
  cout << " Entries = " << nEntries << " DAC values \n" << endl; 
  for (Int_t i = 0; i < nEntries; ++i) {
    cout<< " Run = " << run[i]->GetFirst() << "    DAC = " << run[i]->GetSecond() << endl;
    num_RUN[i] = run[i]->GetFirst();
    injCharge[i] = run[i]->GetSecond();
    injChargeErr[i] = 0.01*injCharge[i];
    if(injChargeErr[i] <= 1.) injChargeErr[i]=1.;
  }
  cout << "" << endl;

  // full print out 

  sprintf(filename,"%s/%s_%d.log",gOutFolder,filenam,gRunNumber);
  logOutputFile_comp=filename;

  filcouc.open(logOutputFile_comp.Data());
  filcouc<<"//====================================================" << endl;
  filcouc<<"//        MUONTRKda for Gain computing (Run = " << gRunNumber << ")" << endl;
  filcouc<<"//====================================================" << endl;
  filcouc<<"//   * Date          : " << date.AsString("l") << "\n" << endl;



  // why 2 files ? (Ch. F.)
  FILE *pfilen = 0;
  if(gPrintLevel==2)
    {
      sprintf(filename,"%s/%s_%d.param",gOutFolder,filenam,gRunNumber);
      cout << " fit parameter file               = " << filename << "\n";
      pfilen = fopen (filename,"w");

      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//  BP MANU CH. par[0]     [1]     [2]     [3]      xlim          P(chi2) p1        P(chi2)2  p2\n");
      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilen,"//===================================================================\n");
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

      fprintf(pfilew,"//================================================\n");
      fprintf(pfilew,"//  Calibration file calculated by MUONTRKda \n");
      fprintf(pfilew,"//=================================================\n");
      fprintf(pfilew,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilew,"//   * Date          : %s \n",date.AsString("l"));
      fprintf(pfilew,"//   * Statictics    : %d \n",gNEvents);
      fprintf(pfilew,"//   * # of MANUS    : %d \n",gNManu);
      fprintf(pfilew,"//   * # of channels : %d \n",gNChannel);
      fprintf(pfilew,"//-------------------------------------------------\n");
      if(nInit==0)
	fprintf(pfilew,"//   %d DAC values  fit:  %d pts (1st order) %d pts (2nd order) \n",nEntries,nbpf1,nbpf2);
      if(nInit==1)
	fprintf(pfilew,"//   %d DAC values  fit: %d pts (1st order) %d pts (2nd order) DAC=0 excluded\n",nEntries,nbpf1,nbpf2);
      fprintf(pfilew,"//   RUN     DAC   \n");
      fprintf(pfilew,"//-----------------\n");
      for (Int_t i = 0; i < nEntries; ++i) {
	tree->SetBranchAddress("run",&run[i]);
	fprintf(pfilew,"//   %d    %5.0f \n",num_RUN[i],injCharge[i]);
      }
      fprintf(pfilew,"//=======================================\n");
      fprintf(pfilew,"// BP MANU CH.   a1      a2     thres. Q\n");
      fprintf(pfilew,"//=======================================\n");
    }

  FILE *pfilep = 0;
  if(gPrintLevel==2)
    {
      sprintf(filename,"%s/%s_%d.peak",gOutFolder,filenam,gRunNumber);
      cout << " File containing Peak mean values = " << filename << "\n";
      pfilep = fopen (filename,"w");

      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//   * Run           : %d \n",gRunNumber); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"// BP  MANU  CH.    Ped.     <0>      <1>      <2>      <3>      <4>      <5>      <6>      <7>      <8>      <9>     <10> \n"); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//                 DAC= %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f  fC\n",injCharge[0],injCharge[1],injCharge[2],injCharge[3],injCharge[4],injCharge[5],injCharge[6],injCharge[7],injCharge[8],injCharge[9],injCharge[10]);
      fprintf(pfilep,"//                      %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f\n",injChargeErr[0],injChargeErr[1],injChargeErr[2],injChargeErr[3],injChargeErr[4],injChargeErr[5],injChargeErr[6],injChargeErr[7],injChargeErr[8],injChargeErr[9],injChargeErr[10]);
      fprintf(pfilep,"//==============================================================================================================================\n");
    }



  //  plot out 

  TFile* gainFile = 0x0;
  sprintf(gRootFileName,"%s/%s_%d.root",gOutFolder,filenam,gRunNumber);
  gainFile = new TFile(gRootFileName,"RECREATE");

  Double_t chi2    = 0.;
  Double_t chi2P2  = 0.;
  Double_t prChi2  = 0; 
  Double_t prChi2P2 =0;
  Double_t a0=0.,a1=1.,a2=0.;
  Int_t busPatchId ;
  Int_t manuId     ;
  Int_t channelId ;
  Int_t threshold = 0;
  Int_t Q = 0;
  Int_t p1 =0;
  Int_t p2 =0;
  Double_t gain=0; 
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

  char graphName[256];

  // iterates over the first pedestal run
  TIter next(map[0]->CreateIterator());
  AliMUONVCalibParam* p;

  Int_t    nmanu         = 0;
  Int_t    nGoodChannel   = 0;
  Int_t    nBadChannel   = 0;
  Int_t    noFitChannel   = 0;
  Int_t    nplot=0;
  Double_t sumProbChi2   = 0.;
  Double_t sumA1         = 0.;
  Double_t sumProbChi2P2 = 0.;
  Double_t sumA2         = 0.;

  Double_t x[11], xErr[11], y[11], yErr[11];
  Double_t xp[11], xpErr[11], yp[11], ypErr[11];

  Int_t uncalcountertotal=0 ;

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
      for ( channelId = 0; channelId < ped[0]->Size() ; ++channelId ) 
	{

	  Int_t n = 0;
	  for (Int_t i = 0; i < nEntries; ++i) {

	    if (!ped[i]) continue; //shouldn't happen.
	    pedMean[i]      = ped[i]->ValueAsDouble(channelId, 0);
	    pedSigma[i]     = ped[i]->ValueAsDouble(channelId, 1);

	    if (pedMean[i] < 0) continue; // not connected

	    if (pedSigma[i] <= 0) pedSigma[i] = 1.; // should not happen.
	    n++;
	  }


	  // print_peak_mean_values
	  if(gPrintLevel==2)
	    {

	      fprintf(pfilep,"%4i%5i%5i%10.3f",busPatchId,manuId,channelId,0.);
	      fprintf(pfilep,"%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f \n",pedMean[0],pedMean[1],pedMean[2],pedMean[3],pedMean[4],pedMean[5],pedMean[6],pedMean[7],pedMean[8],pedMean[9],pedMean[10]);
	      fprintf(pfilep,"                   sig= %9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f \n",pedSigma[0],pedSigma[1],pedSigma[2],pedSigma[3],pedSigma[4],pedSigma[5],pedSigma[6],pedSigma[7],pedSigma[8],pedSigma[9],pedSigma[10]);
	    }

	  // makegain 


	  // Fit Method:  Linear fit over nbpf1 points + parabolic fit  over nbpf2  points) 
	  // nInit=1 : 1st pt DAC=0 excluded

	  // 1. - linear fit over nbpf1 points

	  Double_t par[4] = {0.,0.5,0.,gkADCMax};
	  Int_t nbs   = nEntries - nInit;
	  if(nbs < nbpf1)nbpf1=nbs;

	  Int_t FitProceed=1;
	  for (Int_t j = 0; j < nbs; ++j)
	    {
	      Int_t k = j + nInit;
	      x[j]    = pedMean[k];
	      if(x[j]==0.)FitProceed=0;
	      xErr[j] = pedSigma[k];
	      y[j]    = injCharge[k];
	      yErr[j] = injChargeErr[k];

	    }

	  TGraphErrors *graphErr;
	  if(!FitProceed) { p1=0; p2=0; noFitChannel++;}

	  if(FitProceed)
	    {
		      
	      TF1 *f1 = new TF1("f1",funcLin,0.,gkADCMax,2);
	      graphErr = new TGraphErrors(nbpf1, x, y, xErr, yErr);

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
		  graphErr = new TGraphErrors(nbpf2, xp, yp, xpErr, ypErr);

		  graphErr->Fit(f2,"RQ");
		  chi2P2 = f2->GetChisquare();
		  f2->GetParameters(par);

		  delete graphErr;
		  graphErr=0;
		  delete f2;

		  prChi2P2 = TMath::Prob(chi2P2, nbpf2-1);
		  a2 = par[0];

		  // 	  delete graphErr;

		}

	      par[0] = a0;
	      par[1] = a1;
	      par[2] = a2;
	      par[3] = xLim;

// 	      p1 = TMath::Nint(ceil(prChi2*14))+1;      // round up value : ceil (2.2)=3.
// 	      p2 = TMath::Nint(ceil(prChi2P2*14))+1;
	      if(prChi2>0.999999)prChi2=0.999999 ; if(prChi2P2>0.999999)prChi2P2=0.9999999; // avoiding Pr(Chi2)=1 value
	      p1 = TMath::Nint(floor(prChi2*15))+1;    // round down value : floor(2.8)=2.
	      p2 = TMath::Nint(floor(prChi2P2*15))+1;
	      Q  = p1*16 + p2;  // fit quality 

	      Double_t x0 = -par[0]/par[1]; // value of x corresponding to à 0 fC 
	      threshold = TMath::Nint(ceil(par[3]-x0)); // linear if x < threshold

	      if(gPrintLevel==2)
		{
		  fprintf(pfilen,"%4i %4i %2i",busPatchId,manuId,channelId);
		  fprintf(pfilen," %6.2f %6.4f %10.3e %4.2f %4i          %8.6f %8.6f   %x          %8.6f  %8.6f   %x\n",
			  par[0], par[1], par[2], par[3], threshold, prChi2, floor(prChi2*15), p1,  prChi2P2, floor(prChi2P2*15),p2);
		}
	      // tests
	      if(par[1]< goodA1Min ||  par[1]> goodA1Max) p1=0;
	      if(par[2]< goodA2Min ||  par[2]> goodA2Max) p2=0;

	    } // FitProceed

	  if(FitProceed && p1>0 && p2>0) 
	    {
	      nGoodChannel++;
	      sumProbChi2   += prChi2;
	      sumA1         += par[1];
	      gain=1./(par[1]*capa);
	      sumProbChi2P2   += prChi2P2;
	      sumA2         += par[2];
	    }
	  else // bad calibration
	    {
	      nBadChannel++;
	      Q=0;  
	      par[1]=0.5; a1=0.5; p1=0;
	      par[2]=0.;  a2=0.;  p2=0;
	      threshold=gkADCMax;	

	      char bpmanuname[256];
	      ErrorCounter* uncalcounter;

	      sprintf(bpmanuname,"bp%dmanu%d",busPatchId,manuId);
	      if (!(uncalcounter = (ErrorCounter*)gUncalBuspatchManuTable->FindObject(bpmanuname)))
		{
		  // New buspatch_manu name
		  uncalcounter= new ErrorCounter (busPatchId,manuId);
		  uncalcounter->SetName(bpmanuname);
		  gUncalBuspatchManuTable->Add(uncalcounter);
		}
	      else
		{
		  // Existing buspatch_manu name
		  uncalcounter->Increment();
		}
	      //			    uncalcounter->Print_uncal()
	      uncalcountertotal ++;
	    }

	  if(gPlotLevel){
	    //		      if(Q==0  and  nplot < 100)
	    // 	  if(p1>1 && p2==0  and  nplot < 100)
	    //	    if(p1>1 && p2>1  and  nplot < 100)
	      //	if(p1>=1 and p1<=2  and  nplot < 100)
	    if((p1==1 || p2==1) and  nplot < 100)
	      {
		nplot++;
		// 	      cout << " nplot = " << nplot << endl;
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
		graphErr=0;
		delete f2Calib;
	      }
	  }


	  tg->Fill();

	  if (!flatOutputFile.IsNull()) 
	    {
	      fprintf(pfilew,"%4i %5i %2i %7.4f %10.3e %4i %2x\n",busPatchId,manuId,channelId,par[1],par[2],threshold,Q);
	    }

	}
      nmanu++;
      if(fmod(nmanu,500)==0)std::cout << " Nb manu = " << nmanu << std::endl;
    }

  // file outputs for gain
  if (!flatOutputFile.IsNull())  fclose(pfilew);
  if(gPrintLevel==2){ fclose(pfilen); fclose(pfilep); }

  tg->Write();
  histoFile->Close();

  //OutPut
  if (gPrintLevel) 
    {
      // print in logfile
      if (gUncalBuspatchManuTable->GetSize())
	{
	  gUncalBuspatchManuTable->Sort();  // use compare
	  TIterator* iter = gUncalBuspatchManuTable->MakeIterator();
	  ErrorCounter* uncalcounter;
	  filcouc << "\n List of problematic BusPatch and Manu " << endl;
	  filcouc << " ========================================" << endl;
	  filcouc << "        BP       Manu        Nb Channel  " << endl ;
	  filcouc << " ========================================" << endl;
	  while((uncalcounter = (ErrorCounter*) iter->Next()))
	    {
	      filcouc << "\t" << uncalcounter->BusPatch() << "\t " << uncalcounter->ManuId() << "\t\t"   << uncalcounter->Events() << endl;
	    }
	  filcouc << " ========================================" << endl;

	  filcouc << " Number of bad calibrated Manu    = " << gUncalBuspatchManuTable->GetSize() << endl ;
	  filcouc << " Number of bad calibrated channel = " << uncalcountertotal << endl;
	
	}


      filcouc << "\n Nb of channels in raw data = " << nmanu*64 << " (" << nmanu << " Manu)" <<  endl;
      filcouc << " Nb of calibrated channel   = " << nGoodChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
	      << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
      filcouc << " Nb of uncalibrated channel = " << nBadChannel << " (" << noFitChannel << " unfitted)" << endl;

      cout << "\n Nb of channels in raw data = " << nmanu*64 << " (" << nmanu << " Manu)" <<  endl;
      cout << " Nb of calibrated channel   = " << nGoodChannel << " (" << goodA1Min << "<a1<" << goodA1Max 
	   << " and " << goodA2Min << "<a2<" << goodA2Max << ") " << endl;
      cout << " Nb of uncalibrated channel = " << nBadChannel << " (" << noFitChannel << " unfitted)" << endl;

      Double_t meanA1         = sumA1/(nGoodChannel);
      Double_t meanProbChi2   = sumProbChi2/(nGoodChannel);
      Double_t meanA2         = sumA2/(nGoodChannel);
      Double_t meanProbChi2P2 = sumProbChi2P2/(nGoodChannel);

      Double_t capaManu = 0.2; // pF
      filcouc << "\n linear fit   : <a1> = " << meanA1 << "\t  <gain>  = " <<  1./(meanA1*capaManu) 
	      << " mV/fC (capa= " << capaManu << " pF)" << endl;
      filcouc <<   "        Prob(chi2)>  = " <<  meanProbChi2 << endl;
      filcouc << "\n parabolic fit: <a2> = " << meanA2  << endl;
      filcouc <<   "        Prob(chi2)>  = " <<  meanProbChi2P2 << "\n" << endl;

      cout << "\n  <gain>  = " <<  1./(meanA1*capaManu) 
	   << " mV/fC (capa= " << capaManu << " pF)" 
	   <<  "  Prob(chi2)>  = " <<  meanProbChi2 << endl;
    }  

  filcouc.close();

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

  // 	ofstream filcout;

  Int_t skipEvents = 0;
  Int_t maxEvents  = 1000000;
  Int_t MaxDateEvents  = 1000000;
  Int_t injCharge = 0;
  Char_t inputFile[256]="";

  Int_t gGlitchErrors= 0;
  Int_t gParityErrors= 0;
  Int_t gPaddingErrors= 0;
  Int_t recoverParityErrors = 1;

  TString fesOutputFile;

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
	case 'c' : 
	  i++;
	  gFES=atoi(argv[i]);
	  break;
	case 'd' :
	  i++; 
	  gPrintLevel=atoi(argv[i]);
	  break;
	case 'e' : 
	  i++;
	  gCommand = argv[i];
	  break;
	case 'g' :
	  i++; 
	  gPlotLevel=atoi(argv[i]);
	  break;
	case 'i' :
	  i++; 
	  nbpf1=atoi(argv[i]);
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
	case 'r' : 
	  i++;
	  sprintf(gHistoFileName_gain,argv[i]);
	  break;
	case 'p' : 
	  i++;
	  sscanf(argv[i],"%d",&recoverParityErrors);
	  break;
	case 'h' :
	  i++;
	  printf("\n******************* %s usage **********************",argv[0]);
	  printf("\n%s -options, the available options are :",argv[0]);
	  printf("\n-h help                    (this screen)");
	  printf("\n");
	  printf("\n Input");
	  printf("\n-f <raw data file>         (default = %s)",inputFile); 
	  printf("\n");
	  printf("\n Output");
	  printf("\n-a <Flat ASCII file>       (default = %s)",flatOutputFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-b <output directory>      (default = %s)",gOutFolder);
	  printf("\n-c <FES switch>            (default = %d)",gFES);
	  printf("\n-d <print level>           (default = %d)",gPrintLevel);
	  printf("\n-g <plot level>            (default = %d)",gPlotLevel);
	  printf("\n-i <nb linear points>      (default = %d)",nbpf1);
	  printf("\n-l <DAC level>             (default = %d)",injCharge);
	  printf("\n-m <max date events>       (default = %d)",MaxDateEvents);
	  printf("\n-s <skip events>           (default = %d)",skipEvents);
	  printf("\n-n <max events>            (default = %d)",maxEvents);
	  printf("\n-r root file data for gain (default = %s)",gHistoFileName_gain); 
	  printf("\n-e <execute ped/gain>      (default = %s)",gCommand.Data());
	  printf("\n-e <gain create>           make gain & create a new root file");
	  printf("\n-e <gain>                  make gain & update root file");
	  printf("\n-e <gain compute>          make gain & compute gains");
	  printf("\n-p <Recover parity errors> (default = %d)",recoverParityErrors);

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

  Int_t status=0;
  //  void* event;

  gPedMeanHisto = 0x0;
  gPedSigmaHisto = 0x0;

  TStopwatch timers;

  timers.Start(kTRUE); 

  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;
  TString key("MUONTRKda :");

  // AliMUONRawStreamTrackerHP* rawStream  = 0;

  if (gCommand.CompareTo("comp") != 0)
    {
      
      // Rawdeader, RawStreamHP
      AliRawReader* rawReader = AliRawReader::Create(inputFile);
      AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
      rawStream->DisableWarnings();
      rawStream->EnabbleErrorLogger();

      cout << "\nMUONTRKda : Reading data from file " << inputFile  << endl;

      while (rawReader->NextEvent())
	{
	  if (gNDateEvents >= MaxDateEvents) break;
	  if (gNEvents >= maxEvents) break;
	  if (gNDateEvents>0 &&  gNDateEvents % 100 == 0) 	
	    cout<<"Cumulated:  DATE events = " << gNDateEvents << "   Used events = " << gNEvents << endl;

	  // check shutdown condition 
	  if (daqDA_checkShutdown()) 
	    break;

	  //Skip events
	  while (skipEvents)
	    {
	      rawReader->NextEvent();
	      skipEvents--;
	    }

	  // starts reading
	  // 	      status = monitorGetEventDynamic(&event);
	  // 	      if (status < 0)  {
	  // cout<<"EOF found"<<endl;
	  // break;
	  // 	      }

	  // decoding rawdata headers
	  // AliRawReader *rawReader = new AliRawReaderDate(event);

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
	  // rawStream  = new AliMUONRawStreamTrackerHP(rawReader);
	  // rawStream->DisableWarnings();
	  // rawStream->EnabbleErrorLogger();

	  // First lopp over DDL's to find good events
	  // Error counters per event (counters in the decoding lib are for each DDL)
	  Bool_t eventIsErrorMessage = kFALSE;
	  int eventGlitchErrors = 0;
	  int eventParityErrors = 0;
	  int eventPaddingErrors = 0;
	  rawStream->First();
	  do
	    {
	      if (rawStream->IsErrorMessage()) eventIsErrorMessage = kTRUE;
	      eventGlitchErrors += rawStream->GetGlitchErrors();
	      eventParityErrors += rawStream->GetParityErrors();
	      eventPaddingErrors += rawStream->GetPaddingErrors();
	    } while(rawStream->NextDDL()); 

	  AliMUONRawStreamTrackerHP::AliBusPatch* busPatch;
	  if (!eventIsErrorMessage) 
	    {
	      // Good events (no error) -> compute pedestal for all channels
	      rawStream->First(); 
	      while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
		{
		  for(int i = 0; i < busPatch->GetLength(); ++i)
		    {
		      if (gNEvents == 0) gNChannel++;
		      busPatch->GetData(i, manuId, channelId, charge);
		      MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		    }
		}
	      gNEvents++;
	    }
	  else
	    {
	      // Events with errors
	      if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
		{
		  // Recover parity errors -> compute pedestal for all good buspatches
		  if ( TEST_SYSTEM_ATTRIBUTE( rawReader->GetAttributes(),
					      ATTR_ORBIT_BC )) 
		    {
		      filcout <<"Event recovered -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		    } 
		  else 
		    {
		      filcout <<"Event recovered -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
			      <<" burstNb:"<<EVENT_ID_GET_BURST_NB( rawReader->GetEventId() )
			      <<" nbInBurst:"<<EVENT_ID_GET_NB_IN_BURST( rawReader->GetEventId() )<<endl;
		    }
		  rawStream->First();
		  while( (busPatch = (AliMUONRawStreamTrackerHP::AliBusPatch*) rawStream->Next())) 
		    {
		      // Check the buspatch -> if error not use it in the pedestal calculation
		      int errorCount = 0;
		      for(int i = 0; i < busPatch->GetLength(); ++i)
			{
			  if (!busPatch->IsParityOk(i)) errorCount++;
			}
		      if (!errorCount) 
			{
			  // Good buspatch
			  for(int i = 0; i < busPatch->GetLength(); ++i)
			    {
			      if (gNEvents == 0) gNChannel++;
			      busPatch->GetData(i, manuId, channelId, charge);
			      // if (busPatch->GetBusPatchId()==1719 && manuId == 1 && channelId == 0) cout <<"Recovered charge "<<charge<<endl;
			      MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
			    }
			}
		      else
			{
			  char bpname[256];
			  ErrorCounter* errorCounter;
			  // Bad buspatch -> not used (just print)
			  filcout<<"bpId "<<busPatch->GetBusPatchId()<<" words "<<busPatch->GetLength()
				 <<" parity errors "<<errorCount<<endl;
			  // Number of events where this buspatch is missing
			  sprintf(bpname,"bp%d",busPatch->GetBusPatchId());						
			  if (!(errorCounter = (ErrorCounter*)gErrorBuspatchTable->FindObject(bpname)))
			    {
			      // New buspatch
			      errorCounter = new ErrorCounter(busPatch->GetBusPatchId());
			      errorCounter->SetName(bpname);
			      gErrorBuspatchTable->Add(errorCounter);
			    }
			  else
			    {
			      // Existing buspatch
			      errorCounter->Increment();
			    }	
			  // errorCounter->Print();						
			} // end of if (!errorCount)
		    } // end of while( (busPatch = (AliMUONRawStreamTrackerHP ...
		  gNEvents++;
		  gNEventsRecovered++;
		} //end of if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	      else
		{
		  // Fatal errors reject the event
		  if ( TEST_SYSTEM_ATTRIBUTE( rawReader->GetAttributes(),
					      ATTR_ORBIT_BC )) 
		    {
		      filcout <<"Event rejected -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		    } 
		  else 
		    {
		      filcout <<"Event rejected -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
			      <<" burstNb:"<<EVENT_ID_GET_BURST_NB( rawReader->GetEventId() )
			      <<" nbInBurst:"<<EVENT_ID_GET_NB_IN_BURST( rawReader->GetEventId() )<<endl;

		    }
		} // end of if (!rawStream->GetGlitchErrors() && !rawStream->GetPaddingErrors() ...
	      filcout<<"Number of errors : Glitch "<<eventGlitchErrors
		     <<" Parity "<<eventParityErrors
		     <<" Padding "<<eventPaddingErrors<<endl;
	      filcout<<endl;			
	    } // end of if (!rawStream->IsErrorMessage())

	  if (eventGlitchErrors)  gGlitchErrors++;
	  if (eventParityErrors)  gParityErrors++;
	  if (eventPaddingErrors) gPaddingErrors++;

	} // while (rawReader->NextEvent())
      delete rawReader;
      delete rawStream;


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


      delete gPedestalStore;

      delete minuitFit;
      TVirtualFitter::SetFitter(0);

      timers.Stop();

      cout << "\nMUONTRKda : Nb of DATE events           = " << gNDateEvents    << endl;
      cout << "MUONTRKda : Nb of Glitch errors         = "   << gGlitchErrors  << endl;
      cout << "MUONTRKda : Nb of Parity errors         = "   << gParityErrors  << endl;
      cout << "MUONTRKda : Nb of Padding errors        = "   << gPaddingErrors << endl;		
      cout << "MUONTRKda : Nb of events recovered      = "   << gNEventsRecovered<< endl;
      cout << "MUONTRKda : Nb of events without errors = "   << gNEvents-gNEventsRecovered<< endl;
      cout << "MUONTRKda : Nb of events used           = "   << gNEvents        << endl;

      filcout << "\nMUONTRKda : Nb of DATE events           = " << gNDateEvents    << endl;
      filcout << "MUONTRKda : Nb of Glitch errors         = "   << gGlitchErrors << endl;
      filcout << "MUONTRKda : Nb of Parity errors         = "   << gParityErrors << endl;
      filcout << "MUONTRKda : Nb of Padding errors        = "   << gPaddingErrors << endl;
      filcout << "MUONTRKda : Nb of events recovered      = "   << gNEventsRecovered<< endl;	
      filcout << "MUONTRKda : Nb of events without errors = "   << gNEvents-gNEventsRecovered<< endl;
      filcout << "MUONTRKda : Nb of events used           = "   << gNEvents        << endl;

      if (gCommand.CompareTo("ped") == 0)
	{
          cout << "\nMUONTRKda : Output logfile             : " << logOutputFile  << endl;
	  cout << "MUONTRKda : Pedestal Histo file        : " << gHistoFileName  << endl;
	  cout << "MUONTRKda : Pedestal file (to SHUTTLE) : " << flatOutputFile << endl;   
	}
      else
	{
          cout << "\nMUONTRKda : Output logfile          : " << logOutputFile  << endl;
	  cout << "MUONTRKda : DAC data (root file)    : " << gHistoFileName_gain  << endl;
	  cout << "MUONTRKda : Dummy file (to SHUTTLE) : " << flatOutputFile << endl;   
	}

    }

  // Compute gain parameters


  if (gCommand.Contains("comp")) 
    {
      flatOutputFile="";

      MakeGainStore();

      cout << "\nMUONTRKda : Output logfile          : " << logOutputFile_comp  << endl;
      cout << "MUONTRKda : Root Histo. file        : " << gRootFileName  << endl;
      cout << "MUONTRKda : Gain file (to SHUTTLE)  : " << flatOutputFile << endl;   
    }


  if(gFES) // Store IN FES
    {
      printf("\n *****  STORE FILE in FES ****** \n");

      // be sure that env variable DAQDALIB_PATH is set in script file
      //       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

      if (!flatOutputFile.IsNull()) 
	{
	  if (gCommand.CompareTo("ped") == 0)
	    status = daqDA_FES_storeFile(flatOutputFile.Data(),"PEDESTALS");
	  else
	    status = daqDA_FES_storeFile(flatOutputFile.Data(),"GAINS");

	  if (status) 
	    {
	      printf(" Failed to export file : %d\n",status);
	    }
	  else if(gPrintLevel) printf(" %s successfully exported to FES  \n",flatOutputFile.Data());
	}
    }

  filcout.close();

  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return status;
}
