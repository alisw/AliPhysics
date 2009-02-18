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
#include <sstream>
#include <math.h> 

//AliRoot
// #include "AliMUONLogger.h"
#include "AliMUONRawStreamTrackerHP.h"
// #include "AliMUONDspHeader.h"
// #include "AliMUONBlockHeader.h"
// #include "AliMUONBusStruct.h"
// #include "AliMUONDDLTracker.h"
#include "AliRawReader.h"
#include "AliMUONVStore.h"
#include "AliMUON2DMap.h"
#include "AliMUONCalibParamND.h"
#include "AliMpIntPair.h"
#include "AliMpConstants.h"
#include "AliRawDataErrorLog.h"

#include "AliMUONTrackerIO.h"

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
#include "TObjString.h"
#include "THashTable.h"
#include <THashList.h>
//
//AMORE
//
#ifdef ALI_AMORE
#include <AmoreDA.h>
#endif


#define  NFITPARAMS 4

// global variables
const Int_t kNChannels = AliMpConstants::ManuNofChannels();
const Int_t kADCMax    = 4095;

AliMUONVStore* gAliPedestalStore =  new AliMUON2DMap(kFALSE);

Int_t  gAliNManu       = 0;
Int_t  gAliNChannel    = 0;
UInt_t gAliRunNumber   = 0;
Int_t  gAliNEvents     = 0;
Int_t  gAliNEventsRecovered = 0;
Int_t  gAliPrintLevel  = 1;  // global printout variable (others: 2 and 3)
Int_t  gAliPlotLevel  = 0;  // global plot variable

Char_t gAliHistoFileName[256];

// used for computing gain parameters 
Int_t gAlinbpf1 = 6; // linear fit over nbf1 points

Char_t gAliHistoFileNamegain[256]="MUONTRKda_gain.data";
Char_t gAliOutFolder[256]=".";
Char_t gAlifilename[256];
Char_t gAlifilenam[256]="MUONTRKda_gain"; 
Char_t gAliflatFile[256]="";

ofstream gAlifilcout;

TString gAliOutputFile;
TString gAliCommand("ped");
TTimeStamp gAlidate;

class ErrorCounter : public TNamed
{
public :
  ErrorCounter(Int_t bp = 0, Int_t manu = 0, Int_t ev = 1) : busPatch(bp), manuId(manu), events(ev) {}
  void Increment() {events++;}
  Int_t BusPatch() const {return busPatch;}
  Int_t ManuId() const {return manuId;}
  Int_t Events() const {return events;}
  Int_t Compare(const TObject* obj) const
	{
		/// Compare function
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
	}
	
  void Print(const Option_t* option="") const
  {
    TNamed::Print(option);
    cout<<"bp "<<busPatch<<" events "<<events<<endl;
  }
  void Print_uncal(const Option_t* option="") const
  {
    TNamed::Print(option);
    cout<<"bp =  "<<busPatch<< "  manu = " << manuId << " uncal = "<< events <<endl;
  }

private :
  Int_t busPatch; // Buspath ID
  Int_t manuId;   // Manu ID
  Int_t events;   // Events with error in this buspatch
};

// Table for buspatches with parity errors 
THashTable* gAliErrorBuspatchTable = new THashTable(100,2);

// functions


//________________
Double_t funcLin (const Double_t *x, const Double_t *par)
{
	// Linear function
	return par[0] + par[1]*x[0];
}

//________________
Double_t funcParabolic (const Double_t *x, const Double_t *par)
{
	/// Parabolic function
	return par[0]*x[0]*x[0];
}

//________________
Double_t funcCalib (const Double_t *x, const Double_t *par)  
{
	/// Calibration function
	Double_t xLim= par[3];

	if(x[0] <= xLim) return par[0] + par[1]*x[0];

	Double_t yLim = par[0]+ par[1]*xLim;
	return yLim + par[1]*(x[0] - xLim) + par[2]*(x[0] - xLim)*(x[0] - xLim);
}


//__________
void MakePed(Int_t busPatchId, Int_t manuId, Int_t channelId, Int_t charge)
{
	/// Compute pedestals values
	AliMUONVCalibParam* ped = 
		static_cast<AliMUONVCalibParam*>(gAliPedestalStore->FindObject(busPatchId, manuId));

	if (!ped) {
		gAliNManu++;
		ped = new AliMUONCalibParamND(2, kNChannels,busPatchId, manuId, -1.); // put default wise -1, not connected channel
		gAliPedestalStore->Add(ped);	
	}

	// if (gAliNEvents == 0) {
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
TString WritePedHeader(void) 
{
	ostringstream stream;
	stream<<"//===========================================================================" << endl;
	stream<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
	stream<<"//===========================================================================" << endl;
	stream<<"//       * Run           : " << gAliRunNumber << endl; 
	stream<<"//       * Date          : " << gAlidate.AsString("l") <<endl;
	stream<<"//       * Statictics    : " << gAliNEvents << endl;
	stream<<"//       * # of MANUS    : " << gAliNManu << endl;
	stream<<"//       * # of channels : " << gAliNChannel << endl;
	if (gAliErrorBuspatchTable->GetSize())
	{
		stream<<"//"<<endl;
		stream<<"//       * Buspatches with less statistics (due to parity errors)"<<endl;		
		TIterator* iter = gAliErrorBuspatchTable->MakeIterator();
		ErrorCounter* parityerror;
		while((parityerror = (ErrorCounter*) iter->Next()))
		{
			stream<<"//         bp "<<parityerror->BusPatch()<<" events used "<<gAliNEvents-parityerror->Events()<<endl;
		}
			}	
	stream<<"//"<<endl;
	stream<<"//---------------------------------------------------------------------------" << endl;
	stream<<"//---------------------------------------------------------------------------" << endl;
	stream<<"//      BP     MANU     CH.      MEAN    SIGMA"<<endl;
	stream<<"//---------------------------------------------------------------------------" << endl;
	
	return TString(stream.str().c_str());
}

//________________
TString WritePedData(Int_t BP, Int_t Manu, Int_t ch, Double_t pedMean, Double_t pedSigma) 
{
	ostringstream stream("");
	stream << "\t" << BP << "\t" << Manu <<"\t"<< ch << "\t"
	       << pedMean <<"\t"<< pedSigma << endl;
	return TString(stream.str().c_str());

}

//________________
void MakePedStore(TString gAliOutputFile_1 = "")
{
	
	/// Store pedestals in ASCII files
	Double_t pedMean;
	Double_t pedSigma;
	ofstream fileout;
#ifdef ALI_AMORE
	ostringstream stringout; // String to be sent to AMORE_DB
#endif
	TString tempstring;	
	Int_t busPatchId;
	Int_t manuId;
	Int_t channelId;

// histo
	TFile*  histoFile = 0;
	TTree* tree = 0;
	TH1F* pedMeanHisto = 0;
	TH1F* pedSigmaHisto = 0;
	if (gAliCommand.CompareTo("ped") == 0)
	{
		sprintf(gAliHistoFileName,"%s/MUONTRKda_ped_%d.root",gAliOutFolder,gAliRunNumber);
		histoFile = new TFile(gAliHistoFileName,"RECREATE","MUON Tracking pedestals");

		Char_t name[255];
		Char_t title[255];
		sprintf(name,"pedmean_allch");
		sprintf(title,"Pedestal mean all channels");
		Int_t nx = 4096;
		Int_t xmin = 0;
		Int_t xmax = 4095; 
		pedMeanHisto = new TH1F(name,title,nx,xmin,xmax);
		pedMeanHisto->SetDirectory(histoFile);

		sprintf(name,"pedsigma_allch");
		sprintf(title,"Pedestal sigma all channels");
		nx = 201;
		xmin = 0;
		xmax = 200; 
		pedSigmaHisto = new TH1F(name,title,nx,xmin,xmax);
		pedSigmaHisto->SetDirectory(histoFile);

		tree = new TTree("t","Pedestal tree");
		tree->Branch("bp",&busPatchId,"bp/I");
		tree->Branch("manu",&manuId,",manu/I");
		tree->Branch("channel",&channelId,",channel/I");
		tree->Branch("pedMean",&pedMean,",pedMean/D");
		tree->Branch("pedSigma",&pedSigma,",pedSigma/D");
	}

	if (!gAliOutputFile_1.IsNull()) {
		fileout.open(gAliOutputFile_1.Data());
		tempstring = WritePedHeader();
		fileout << tempstring;
#ifdef ALI_AMORE
		stringout << tempstring;
#endif
	}
	// print in logfile
	if (gAliErrorBuspatchTable->GetSize())
	  {
	    cout<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;		
	    gAlifilcout<<"\n* Buspatches with less statistics (due to parity errors)"<<endl;		
	    TIterator* iter = gAliErrorBuspatchTable->MakeIterator();
	    ErrorCounter* parityerror;
	    while((parityerror = (ErrorCounter*) iter->Next()))
	      {
		cout<<"  bp "<<parityerror->BusPatch()<<": events used = "<<gAliNEvents-parityerror->Events()<<endl;
		gAlifilcout<<"  bp "<<parityerror->BusPatch()<<": events used = "<<gAliNEvents-parityerror->Events()<<endl;
	      }
	
	  }

	
// iterator over pedestal
	TIter next(gAliPedestalStore->CreateIterator());
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
		if ((errorCounter = (ErrorCounter*)gAliErrorBuspatchTable->FindObject(bpname)))
		{
			eventCounter = gAliNEvents - errorCounter->Events();
		}
		else
		{
			eventCounter = gAliNEvents;
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


			if (!gAliOutputFile_1.IsNull()) {
				tempstring = WritePedData(busPatchId,manuId,channelId,pedMean,pedSigma);
				fileout << tempstring;
#ifdef ALI_AMORE
				stringout << tempstring;
#endif
			}

			if (gAliCommand.CompareTo("ped") == 0)
			{
				pedMeanHisto->Fill(pedMean);
				pedSigmaHisto->Fill(pedSigma);

				tree->Fill();
			}
		}
	}
}

// file outputs
if (!gAliOutputFile_1.IsNull())  fileout.close();

// Outputs to root file and AMORE DB
if (gAliCommand.CompareTo("ped") == 0)
{
#ifdef ALI_AMORE
  //
  //Send objects to the AMORE DB
  //
  const char *role=gSystem->Getenv("AMORE_DA_NAME");
  if ( role ){
    amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
//     TObjString peddata(stringout.str().c_str());
    TObjString peddata(stringout.str().c_str());
    Int_t status =0;
    status = amoreDA.Send("Pedestals",&peddata);
    if ( status )
      cout << "Warning: Failed to write Pedestals in the AMORE database : " << status << endl;
    // reset env var
    if (amoreDANameorig) gSystem->Setenv("AMORE_DA_NAME",amoreDANameorig);
    } 
  else {
    cout << "Warning: environment variable 'AMORE_DA_NAME' not set. Cannot write to the AMORE database" << endl;
    }
#endif
	histoFile->Write();
	histoFile->Close();
}

//  delete tree;

}

//________________
void MakePedStoreForGain(Int_t injCharge)
{
	/// Store pedestal map in root file
	TTree* tree = 0x0;

	FILE *pfilew=0;
	if (gAliCommand.Contains("gain") && !gAliCommand.Contains("comp")) {
		if(gAliOutputFile.IsNull())
		{
			sprintf(gAlifilename,"%s_%d_DAC_%d.par",gAlifilenam,gAliRunNumber,injCharge);
			gAliOutputFile=gAlifilename;
		}
		if(!gAliOutputFile.IsNull())
		{
			pfilew = fopen (gAliOutputFile.Data(),"w");

			fprintf(pfilew,"//DUMMY FILE (to prevent Shuttle failure)\n");
			fprintf(pfilew,"//================================================\n");
			fprintf(pfilew,"//       MUONTRKda: Calibration run  \n");
			fprintf(pfilew,"//=================================================\n");
			fprintf(pfilew,"//   * Run           : %d \n",gAliRunNumber); 
			fprintf(pfilew,"//   * Date          : %s \n",gAlidate.AsString("l"));
			fprintf(pfilew,"//   * DAC           : %d \n",injCharge);
			fprintf(pfilew,"//-------------------------------------------------\n");
			fclose(pfilew);
		}
	}

	if(gAliPrintLevel==2)
	{
		// compute and store pedestals
		sprintf(gAliflatFile,"%s/%s_%d_DAC_%d.ped",gAliOutFolder,gAlifilenam,gAliRunNumber,injCharge);
		cout << "\nMUONTRKda : Flat file  generated  : " << gAliflatFile << "\n";
		MakePedStore(gAliflatFile);
	}
	else
		MakePedStore();

	TString mode("UPDATE");

	if (gAliCommand.Contains("cre")) {
		mode = "RECREATE";
	}
	TFile* histoFile = new TFile(gAliHistoFileNamegain, mode.Data(), "MUON Tracking Gains");

	// second argument should be the injected charge, taken from config crocus file
	// put also info about run number could be usefull
	AliMpIntPair* pair   = new AliMpIntPair(gAliRunNumber, injCharge);

	if (mode.CompareTo("UPDATE") == 0) {
		tree = (TTree*)histoFile->Get("t");
		tree->SetBranchAddress("run",&pair);
		tree->SetBranchAddress("ped",&gAliPedestalStore);

	} else {
		tree = new TTree("t","Pedestal tree");
		tree->Branch("run", "AliMpIntPair",&pair);
		tree->Branch("ped", "AliMUON2DMap",&gAliPedestalStore);
		tree->SetBranchAddress("run",&pair);
		tree->SetBranchAddress("ped",&gAliPedestalStore);

	}

	tree->Fill();
	tree->Write("t", TObject::kOverwrite); // overwrite the tree
	histoFile->Close();

	delete pair;
}

//________________
TString WriteGainHeader(Int_t nInit, Int_t nEntries, Int_t nbpf2, Int_t *numrun, Double_t *injCharge) 
{
      ostringstream stream;
      stream << "//================================================" << endl;
      stream << "//  Calibration file calculated by MUONTRKda " << endl;
      stream << "//=================================================" << endl;
      stream << "//   * Run           : " << gAliRunNumber << endl; 
      stream << "//   * Date          : " << gAlidate.AsString("l") << endl;
      stream << "//   * Statictics    : " << gAliNEvents << endl;
      stream << "//   * # of MANUS    : " << gAliNManu << endl;
      stream << "//   * # of channels : " << gAliNChannel << endl;
      stream << "//-------------------------------------------------" << endl;
      if(nInit==0)
	stream << "//   " << nEntries << " DAC values  fit:  " << gAlinbpf1 << " pts (1st order) " << nbpf2 << " pts (2nd order)" << endl;
      if(nInit==1)
	stream << "//   " << nEntries << " DAC values  fit: " << gAlinbpf1 << " pts (1st order) " << nbpf2 << " pts (2nd order) DAC=0 excluded" << endl;
      stream << "//   RUN     DAC   " << endl;
      stream << "//-----------------" << endl;
      for (Int_t i = 0; i < nEntries; ++i) stream << Form("//   %d   %5.0f",numrun[i],injCharge[i]) << endl;
      stream << "//=======================================" << endl;
      stream << "// BP MANU CH.   a1      a2     thres. q" << endl;
      stream << "//=======================================" << endl;
      return TString(stream.str().c_str());
}

//________________
TString WriteGainData(Int_t busPatchId, Int_t manuId, Int_t channelId, Double_t par1, Double_t par2, Int_t threshold, Int_t q)
{
	ostringstream stream("");
	stream << Form("%4i %5i %2i %7.4f %10.3e %4i %2x",busPatchId,manuId,channelId,par1,par2,threshold,q) << endl;
	return TString(stream.str().c_str());
}

//________________
void MakeGainStore()
{
	/// Store gains in ASCII files
  ofstream filcouc;
  TString tempstring;

  Int_t nInit = 1; // DAC=0 excluded from fit procedure
  Double_t goodA1Min =  0.5;
  Double_t goodA1Max =  2.;
  //     Double_t goodA1Min =  0.7;
  //     Double_t goodA1Max =  1.7;
  Double_t goodA2Min = -0.5E-03;
  Double_t goodA2Max =  1.E-03;
	Char_t rootFileName[256];
	TString logOutputFilecomp;
	// Table for uncalibrated  buspatches and manus
 	THashList* uncalBuspatchManuTable = new THashList(1000,2);

  Int_t numrun[15];

  // open file mutrkgain.root
  // read again the pedestal for the calibration runs (9 runs ?)
  // need the injection charge from config file (to be done)
  // For each channel make a TGraphErrors (mean, sigma) vs injected charge
  // Fit with a polynomial fct
  // store the result in a flat file.


  TFile*  histoFile = new TFile(gAliHistoFileNamegain);

  AliMUON2DMap* map[11];
  AliMUONVCalibParam* ped[11];
  AliMpIntPair* run[11];

  //read back from root file
  TTree* tree = (TTree*)histoFile->Get("t");
  Int_t nEntries = tree->GetEntries();
  Int_t nbpf2 = nEntries - (nInit + gAlinbpf1) + 1; // nb pts used for 2nd order fit

  // read back info
  for (Int_t i = 0; i < nEntries; ++i) {
    map[i] = 0x0;
    run[i] = 0x0;
    tree->SetBranchAddress("ped",&map[i]);
    tree->SetBranchAddress("run",&run[i]);
    tree->GetEvent(i);
    //        std::cout << map[i] << " " << run[i] << std::endl;
  }
  //jlc_feb_08  modif:   gAliRunNumber=(UInt_t)run[0]->GetFirst();
  gAliRunNumber=(UInt_t)run[nEntries-1]->GetFirst();
  //     sscanf(getenv("DATE_RUN_NUMBER"),"%d",&gAliRunNumber);

  Double_t pedMean[11];
  Double_t pedSigma[11];
  for ( Int_t i=0 ; i<11 ; i++) {pedMean[i]=0.;pedSigma[i]=1.;};
  Double_t injCharge[11];
  Double_t injChargeErr[11];
  for ( Int_t i=0 ; i<11 ; i++) {injCharge[i]=0.;injChargeErr[i]=1.;};

  // some print
  cout<<"\n ********  MUONTRKda for Gain computing (Run = " << gAliRunNumber << ")\n" << endl;
  cout<<" * Date          : " << gAlidate.AsString("l") << "\n" << endl;
  cout << " Entries = " << nEntries << " DAC values \n" << endl; 
  for (Int_t i = 0; i < nEntries; ++i) {
    cout<< " Run = " << run[i]->GetFirst() << "    DAC = " << run[i]->GetSecond() << endl;
    numrun[i] = run[i]->GetFirst();
    injCharge[i] = run[i]->GetSecond();
    injChargeErr[i] = 0.01*injCharge[i];
    if(injChargeErr[i] <= 1.) injChargeErr[i]=1.;
  }
  cout << "" << endl;

  // full print out 

  sprintf(gAlifilename,"%s/%s_%d.log",gAliOutFolder,gAlifilenam,gAliRunNumber);
  logOutputFilecomp=gAlifilename;

  filcouc.open(logOutputFilecomp.Data());
  filcouc<<"//====================================================" << endl;
  filcouc<<"//        MUONTRKda for Gain computing (Run = " << gAliRunNumber << ")" << endl;
  filcouc<<"//====================================================" << endl;
  filcouc<<"//   * Date          : " << gAlidate.AsString("l") << "\n" << endl;



  // why 2 files ? (Ch. F.)
  FILE *pfilen = 0;
  if(gAliPrintLevel==2)
    {
      sprintf(gAlifilename,"%s/%s_%d.param",gAliOutFolder,gAlifilenam,gAliRunNumber);
      cout << " fit parameter file               = " << gAlifilename << "\n";
      pfilen = fopen (gAlifilename,"w");

      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//  BP MANU CH. par[0]     [1]     [2]     [3]      xlim          P(chi2) p1        P(chi2)2  p2\n");
      fprintf(pfilen,"//===================================================================\n");
      fprintf(pfilen,"//   * Run           : %d \n",gAliRunNumber); 
      fprintf(pfilen,"//===================================================================\n");
    }

  ofstream pfilew;
#ifdef ALI_AMORE
  ostringstream pstringw;
#endif
  if(gAliOutputFile.IsNull())
    {
      sprintf(gAlifilename,"%s_%d.par",gAlifilenam,gAliRunNumber);
      gAliOutputFile=gAlifilename;
    }
  if(!gAliOutputFile.IsNull())
    {    
      pfilew.open(gAliOutputFile.Data());
      pfilew << WriteGainHeader(nInit,nEntries,nbpf2,numrun,injCharge);
#ifdef ALI_AMORE
      pstringw << WriteGainHeader(nInit,nEntries,nbpf2,numrun,injCharge);
#endif
      for (Int_t i = 0; i < nEntries; ++i) {
	tree->SetBranchAddress("run",&run[i]);
      }
    }

  FILE *pfilep = 0;
  if(gAliPrintLevel==2)
    {
      sprintf(gAlifilename,"%s/%s_%d.peak",gAliOutFolder,gAlifilenam,gAliRunNumber);
      cout << " File containing Peak mean values = " << gAlifilename << "\n";
      pfilep = fopen (gAlifilename,"w");

      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//   * Run           : %d \n",gAliRunNumber); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"// BP  MANU  CH.    Ped.     <0>      <1>      <2>      <3>      <4>      <5>      <6>      <7>      <8>      <9>     <10> \n"); 
      fprintf(pfilep,"//==============================================================================================================================\n");
      fprintf(pfilep,"//                 DAC= %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f  fC\n",injCharge[0],injCharge[1],injCharge[2],injCharge[3],injCharge[4],injCharge[5],injCharge[6],injCharge[7],injCharge[8],injCharge[9],injCharge[10]);
      fprintf(pfilep,"//                      %9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f\n",injChargeErr[0],injChargeErr[1],injChargeErr[2],injChargeErr[3],injChargeErr[4],injChargeErr[5],injChargeErr[6],injChargeErr[7],injChargeErr[8],injChargeErr[9],injChargeErr[10]);
      fprintf(pfilep,"//==============================================================================================================================\n");
    }



  //  plot out 

  TFile* gainFile = 0x0;
  sprintf(rootFileName,"%s/%s_%d.root",gAliOutFolder,gAlifilenam,gAliRunNumber);
  gainFile = new TFile(rootFileName,"RECREATE");

  Double_t chi2    = 0.;
  Double_t chi2P2  = 0.;
  Double_t prChi2  = 0; 
  Double_t prChi2P2 =0;
  Double_t a0=0.,a1=1.,a2=0.;
  Int_t busPatchId ;
  Int_t manuId     ;
  Int_t channelId ;
  Int_t threshold = 0;
  Int_t q = 0;
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
  tg->Branch("q",&q, "q/I");
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
	  if(gAliPrintLevel==2)
	    {

	      fprintf(pfilep,"%4i%5i%5i%10.3f",busPatchId,manuId,channelId,0.);
	      fprintf(pfilep,"%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f%9.1f \n",pedMean[0],pedMean[1],pedMean[2],pedMean[3],pedMean[4],pedMean[5],pedMean[6],pedMean[7],pedMean[8],pedMean[9],pedMean[10]);
	      fprintf(pfilep,"                   sig= %9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f \n",pedSigma[0],pedSigma[1],pedSigma[2],pedSigma[3],pedSigma[4],pedSigma[5],pedSigma[6],pedSigma[7],pedSigma[8],pedSigma[9],pedSigma[10]);
	    }

	  // makegain 


	  // Fit Method:  Linear fit over gAlinbpf1 points + parabolic fit  over nbpf2  points) 
	  // nInit=1 : 1st pt DAC=0 excluded

	  // 1. - linear fit over gAlinbpf1 points

	  Double_t par[4] = {0.,0.5,0.,kADCMax};
	  Int_t nbs   = nEntries - nInit;
	  if(nbs < gAlinbpf1)gAlinbpf1=nbs;

	  Int_t fitproceed=1;
	  for (Int_t j = 0; j < nbs; ++j)
	    {
	      Int_t k = j + nInit;
	      x[j]    = pedMean[k];
	      if(x[j]==0.)fitproceed=0;
	      xErr[j] = pedSigma[k];
	      y[j]    = injCharge[k];
	      yErr[j] = injChargeErr[k];

	    }

	  TGraphErrors *graphErr;
	  if(!fitproceed) { p1=0; p2=0; noFitChannel++;}

	  if(fitproceed)
	    {
		      
	      TF1 *f1 = new TF1("f1",funcLin,0.,kADCMax,2);
	      graphErr = new TGraphErrors(gAlinbpf1, x, y, xErr, yErr);

	      f1->SetParameters(0,0);

	      graphErr->Fit("f1","RQ");

	      chi2 = f1->GetChisquare();
	      f1->GetParameters(par);

	      delete graphErr;
	      graphErr=0;
	      delete f1;

	      prChi2 = TMath::Prob(chi2, gAlinbpf1 - 2);

	      Double_t xLim = pedMean[nInit + gAlinbpf1 - 1];
	      Double_t yLim = par[0]+par[1] * xLim;

	      a0 = par[0];
	      a1 = par[1];

	      // 2. - Translation : new origin (xLim, yLim) + parabolic fit over nbf2 points

	      if(nbpf2 > 1)
		{
		  for (Int_t j = 0; j < nbpf2; j++)
		    {
		      Int_t k  = j + (nInit + gAlinbpf1) - 1;
		      xp[j]    = pedMean[k] - xLim;
		      xpErr[j] = pedSigma[k];

		      yp[j]    = injCharge[k] - yLim - par[1]*xp[j];
		      ypErr[j] = injChargeErr[k];
		    }

		  TF1 *f2 = new TF1("f2",funcParabolic,0.,kADCMax,1);
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
	      q  = p1*16 + p2;  // fit quality 

	      Double_t x0 = -par[0]/par[1]; // value of x corresponding to à 0 fC 
	      threshold = TMath::Nint(ceil(par[3]-x0)); // linear if x < threshold

	      if(gAliPrintLevel==2)
		{
		  fprintf(pfilen,"%4i %4i %2i",busPatchId,manuId,channelId);
		  fprintf(pfilen," %6.2f %6.4f %10.3e %4.2f %4i          %8.6f %8.6f   %x          %8.6f  %8.6f   %x\n",
			  par[0], par[1], par[2], par[3], threshold, prChi2, floor(prChi2*15), p1,  prChi2P2, floor(prChi2P2*15),p2);
		}
	      // tests
	      if(par[1]< goodA1Min ||  par[1]> goodA1Max) p1=0;
	      if(par[2]< goodA2Min ||  par[2]> goodA2Max) p2=0;

	    } // fitproceed

	  if(fitproceed && p1>0 && p2>0) 
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
	      q=0;  
	      par[1]=0.5; a1=0.5; p1=0;
	      par[2]=0.;  a2=0.;  p2=0;
	      threshold=kADCMax;	

	      char bpmanuname[256];
	      ErrorCounter* uncalcounter;

	      sprintf(bpmanuname,"bp%dmanu%d",busPatchId,manuId);
	      if (!(uncalcounter = (ErrorCounter*)uncalBuspatchManuTable->FindObject(bpmanuname)))
		{
		  // New buspatch_manu name
		  uncalcounter= new ErrorCounter (busPatchId,manuId);
		  uncalcounter->SetName(bpmanuname);
		  uncalBuspatchManuTable->Add(uncalcounter);
		}
	      else
		{
		  // Existing buspatch_manu name
		  uncalcounter->Increment();
		}
	      //			    uncalcounter->Print_uncal()
	      uncalcountertotal ++;
	    }

	  if(gAliPlotLevel){
	    //		      if(q==0  and  nplot < 100)
	    // 	  if(p1>1 && p2==0  and  nplot < 100)
	    //	    if(p1>1 && p2>1  and  nplot < 100)
	      //	if(p1>=1 and p1<=2  and  nplot < 100)
	    if((p1==1 || p2==1) and  nplot < 100)
	      {
		nplot++;
		// 	      cout << " nplot = " << nplot << endl;
		TF1 *f2Calib = new TF1("f2Calib",funcCalib,0.,kADCMax,NFITPARAMS);

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

	  if (!gAliOutputFile.IsNull()) 
	    {
	      tempstring = WriteGainData(busPatchId,manuId,channelId,par[1],par[2],threshold,q);
	      if(manuId && (busPatchId!=1621)) {// Add a protection to avoid a future crash in Amore due to manuId = 0 (bug not understood/fixed yet)
	      	pfilew << tempstring;
#ifdef ALI_AMORE
	      	pstringw << tempstring;
#endif
	      }
	    }

	}
      nmanu++;
      if(fmod(nmanu,500)==0)std::cout << " Nb manu = " << nmanu << std::endl;
    }

  // outputs for gain (file + AMORE DB)
  if (!gAliOutputFile.IsNull())  {
  	pfilew.close();
#ifdef ALI_AMORE
  //
  //Send objects to the AMORE DB
  //
  const char *role=gSystem->Getenv("AMORE_DA_NAME");
  if ( role ){
    	amore::da::AmoreDA amoreDA(amore::da::AmoreDA::kSender);
	TObjString gaindata(pstringw.str().c_str());
    	if ( amoreDA.Send("Gains",&gaindata) )
      	   cout << "Warning: Failed to write Gains to the AMORE database" << endl;
    // reset env var
    	} 
  else {
	cout << "Warning: environment variable 'AMORE_DA_NAME' not set. Cannot write to the AMORE database" << endl;
    	}
#endif
  }
  if(gAliPrintLevel==2){ fclose(pfilen); fclose(pfilep); }

  tg->Write();
  histoFile->Close();

  //OutPut
  if (gAliPrintLevel) 
    {
      // print in logfile
      if (uncalBuspatchManuTable->GetSize())
	{
	  uncalBuspatchManuTable->Sort();  // use compare
	  TIterator* iter = uncalBuspatchManuTable->MakeIterator();
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

	  filcouc << " Number of bad calibrated Manu    = " << uncalBuspatchManuTable->GetSize() << endl ;
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
  cout << "\nMUONTRKda : Output logfile          : " << logOutputFilecomp  << endl;
  cout << "MUONTRKda : Root Histo. file        : " << rootFileName  << endl;
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

  // 	ofstream gAlifilcout;

  Int_t fes  = 1;  // by default FES is used
  Int_t skipEvents = 0;
  Int_t maxEvents  = 1000000;
  Int_t maxDateEvents  = 1000000;
  Int_t injCharge = 0;
  Char_t inputFile[256]="";

	Int_t  nDateEvents = 0;
  Int_t gGlitchErrors= 0;
  Int_t gParityErrors= 0;
  Int_t gPaddingErrors= 0;
  Int_t recoverParityErrors = 1;

  TString fesOutputFile;
	TString logOutputFile;

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
	  gAliOutputFile = argv[i];
	  break;
	case 'b' : 
	  i++;
	  sprintf(gAliOutFolder,argv[i]);
	  break;
	case 'c' : 
	  i++;
	  fes=atoi(argv[i]);
	  break;
	case 'd' :
	  i++; 
	  gAliPrintLevel=atoi(argv[i]);
	  break;
	case 'e' : 
	  i++;
	  gAliCommand = argv[i];
	  break;
	case 'g' :
	  i++; 
	  gAliPlotLevel=atoi(argv[i]);
	  break;
	case 'i' :
	  i++; 
	  gAlinbpf1=atoi(argv[i]);
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
	  sscanf(argv[i],"%d",&maxDateEvents);
	  break;
	case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
	case 'r' : 
	  i++;
	  sprintf(gAliHistoFileNamegain,argv[i]);
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
	  printf("\n-a <Flat ASCII file>       (default = %s)",gAliOutputFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-b <output directory>      (default = %s)",gAliOutFolder);
	  printf("\n-c <FES switch>            (default = %d)",fes);
	  printf("\n-d <print level>           (default = %d)",gAliPrintLevel);
	  printf("\n-g <plot level>            (default = %d)",gAliPlotLevel);
	  printf("\n-i <nb linear points>      (default = %d)",gAlinbpf1);
	  printf("\n-l <DAC level>             (default = %d)",injCharge);
	  printf("\n-m <max date events>       (default = %d)",maxDateEvents);
	  printf("\n-s <skip events>           (default = %d)",skipEvents);
	  printf("\n-n <max events>            (default = %d)",maxEvents);
	  printf("\n-r root file data for gain (default = %s)",gAliHistoFileNamegain); 
	  printf("\n-e <execute ped/gain>      (default = %s)",gAliCommand.Data());
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

  // set gAliCommand to lower case
  gAliCommand.ToLower();


  // decoding the events

  Int_t status=0;
  //  void* event;

  // gAliPedMeanHisto = 0x0;
  // gAliPedSigmaHisto = 0x0;

  TStopwatch timers;

  timers.Start(kTRUE); 

  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;
  TString key("MUONTRKda :");

  // AliMUONRawStreamTrackerHP* rawStream  = 0;

  if (gAliCommand.CompareTo("comp") != 0)
    {
      
      // Rawdeader, RawStreamHP
      AliRawReader* rawReader = AliRawReader::Create(inputFile);
      AliMUONRawStreamTrackerHP* rawStream  = new AliMUONRawStreamTrackerHP(rawReader);    
      rawStream->DisableWarnings();
      rawStream->EnabbleErrorLogger();

      cout << "\nMUONTRKda : Reading data from file " << inputFile  << endl;

      while (rawReader->NextEvent())
	{
	  if (nDateEvents >= maxDateEvents) break;
	  if (gAliNEvents >= maxEvents) break;
	  if (nDateEvents>0 &&  nDateEvents % 100 == 0) 	
	    cout<<"Cumulated:  DATE events = " << nDateEvents << "   Used events = " << gAliNEvents << endl;

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
	  gAliRunNumber = rawReader->GetRunNumber();

	  // Output log file initialisations

	  if(nDateEvents==0)
	    {
	      if (gAliCommand.CompareTo("ped") == 0){
		sprintf(gAliflatFile,"%s/MUONTRKda_ped_%d.log",gAliOutFolder,gAliRunNumber);
		logOutputFile=gAliflatFile;

		gAlifilcout.open(logOutputFile.Data());
		gAlifilcout<<"//=================================================" << endl;
		gAlifilcout<<"//        MUONTRKda for Pedestal run = "   << gAliRunNumber << endl;
		cout<<"\n ********  MUONTRKda for Pedestal run = " << gAliRunNumber << "\n" << endl;
	      }

	      if (gAliCommand.Contains("gain")){
		sprintf(gAliflatFile,"%s/%s_%d_DAC_%d.log",gAliOutFolder,gAlifilenam,gAliRunNumber,injCharge);
		logOutputFile=gAliflatFile;

		gAlifilcout.open(logOutputFile.Data());
		gAlifilcout<<"//=================================================" << endl;
		gAlifilcout<<"//        MUONTRKda for Gain run = " << gAliRunNumber << "  (DAC=" << injCharge << ")" << endl;
		cout<<"\n ********  MUONTRKda for Gain run = " << gAliRunNumber << "  (DAC=" << injCharge << ")\n" << endl;
	      }

	      gAlifilcout<<"//=================================================" << endl;
	      gAlifilcout<<"//   * Date          : " << gAlidate.AsString("l") << "\n" << endl;
	      cout<<" * Date          : " << gAlidate.AsString("l") << "\n" << endl;

	    }

	  nDateEvents++;

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
		      if (gAliNEvents == 0) gAliNChannel++;
		      busPatch->GetData(i, manuId, channelId, charge);
		      MakePed(busPatch->GetBusPatchId(), (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		    }
		}
	      gAliNEvents++;
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
		      gAlifilcout <<"Event recovered -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		    } 
		  else 
		    {
		      gAlifilcout <<"Event recovered -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
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
			      if (gAliNEvents == 0) gAliNChannel++;
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
			  gAlifilcout<<"bpId "<<busPatch->GetBusPatchId()<<" words "<<busPatch->GetLength()
				 <<" parity errors "<<errorCount<<endl;
			  // Number of events where this buspatch is missing
			  sprintf(bpname,"bp%d",busPatch->GetBusPatchId());						
			  if (!(errorCounter = (ErrorCounter*)gAliErrorBuspatchTable->FindObject(bpname)))
			    {
			      // New buspatch
			      errorCounter = new ErrorCounter(busPatch->GetBusPatchId());
			      errorCounter->SetName(bpname);
			      gAliErrorBuspatchTable->Add(errorCounter);
			    }
			  else
			    {
			      // Existing buspatch
			      errorCounter->Increment();
			    }	
			  // errorCounter->Print();						
			} // end of if (!errorCount)
		    } // end of while( (busPatch = (AliMUONRawStreamTrackerHP ...
		  gAliNEvents++;
		  gAliNEventsRecovered++;
		} //end of if (recoverParityErrors && eventParityErrors && !eventGlitchErrors&& !eventPaddingErrors)
	      else
		{
		  // Fatal errors reject the event
		  if ( TEST_SYSTEM_ATTRIBUTE( rawReader->GetAttributes(),
					      ATTR_ORBIT_BC )) 
		    {
		      gAlifilcout <<"Event rejected -> Period:"<<EVENT_ID_GET_PERIOD( rawReader->GetEventId() )
			      <<" Orbit:"<<EVENT_ID_GET_ORBIT( rawReader->GetEventId() )
			      <<" BunchCrossing:"<<EVENT_ID_GET_BUNCH_CROSSING( rawReader->GetEventId() )<<endl;				
		    } 
		  else 
		    {
		      gAlifilcout <<"Event rejected -> nbInRun:"<<EVENT_ID_GET_NB_IN_RUN( rawReader->GetEventId() )
			      <<" burstNb:"<<EVENT_ID_GET_BURST_NB( rawReader->GetEventId() )
			      <<" nbInBurst:"<<EVENT_ID_GET_NB_IN_BURST( rawReader->GetEventId() )<<endl;

		    }
		} // end of if (!rawStream->GetGlitchErrors() && !rawStream->GetPaddingErrors() ...
	      gAlifilcout<<"Number of errors : Glitch "<<eventGlitchErrors
		     <<" Parity "<<eventParityErrors
		     <<" Padding "<<eventPaddingErrors<<endl;
	      gAlifilcout<<endl;			
	    } // end of if (!rawStream->IsErrorMessage())

	  if (eventGlitchErrors)  gGlitchErrors++;
	  if (eventParityErrors)  gParityErrors++;
	  if (eventPaddingErrors) gPaddingErrors++;

	} // while (rawReader->NextEvent())
      delete rawReader;
      delete rawStream;


      if (gAliCommand.CompareTo("ped") == 0)
	{
	  sprintf(gAliflatFile,"MUONTRKda_ped_%d.ped",gAliRunNumber);
	  if(gAliOutputFile.IsNull())gAliOutputFile=gAliflatFile;
	  MakePedStore(gAliOutputFile);
	}

      // option gain -> update root file with pedestal results
      // gain + create -> recreate root file
      // gain + comp -> update root file and compute gain parameters

      if (gAliCommand.Contains("gain")) 
	{
	  MakePedStoreForGain(injCharge);
	}


      delete gAliPedestalStore;

      delete minuitFit;
      TVirtualFitter::SetFitter(0);

      timers.Stop();

      cout << "\nMUONTRKda : Nb of DATE events           = " << nDateEvents    << endl;
      cout << "MUONTRKda : Nb of Glitch errors         = "   << gGlitchErrors  << endl;
      cout << "MUONTRKda : Nb of Parity errors         = "   << gParityErrors  << endl;
      cout << "MUONTRKda : Nb of Padding errors        = "   << gPaddingErrors << endl;		
      cout << "MUONTRKda : Nb of events recovered      = "   << gAliNEventsRecovered<< endl;
      cout << "MUONTRKda : Nb of events without errors = "   << gAliNEvents-gAliNEventsRecovered<< endl;
      cout << "MUONTRKda : Nb of events used           = "   << gAliNEvents        << endl;

      gAlifilcout << "\nMUONTRKda : Nb of DATE events           = " << nDateEvents    << endl;
      gAlifilcout << "MUONTRKda : Nb of Glitch errors         = "   << gGlitchErrors << endl;
      gAlifilcout << "MUONTRKda : Nb of Parity errors         = "   << gParityErrors << endl;
      gAlifilcout << "MUONTRKda : Nb of Padding errors        = "   << gPaddingErrors << endl;
      gAlifilcout << "MUONTRKda : Nb of events recovered      = "   << gAliNEventsRecovered<< endl;	
      gAlifilcout << "MUONTRKda : Nb of events without errors = "   << gAliNEvents-gAliNEventsRecovered<< endl;
      gAlifilcout << "MUONTRKda : Nb of events used           = "   << gAliNEvents        << endl;

      if (gAliCommand.CompareTo("ped") == 0)
	{
          cout << "\nMUONTRKda : Output logfile             : " << logOutputFile  << endl;
	  cout << "MUONTRKda : Pedestal Histo file        : " << gAliHistoFileName  << endl;
	  cout << "MUONTRKda : Pedestal file (to SHUTTLE) : " << gAliOutputFile << endl;   
	}
      else
	{
          cout << "\nMUONTRKda : Output logfile          : " << logOutputFile  << endl;
	  cout << "MUONTRKda : DAC data (root file)    : " << gAliHistoFileNamegain  << endl;
	  cout << "MUONTRKda : Dummy file (to SHUTTLE) : " << gAliOutputFile << endl;   
	}

    }

  // Compute gain parameters


  if (gAliCommand.Contains("comp")) 
    {
      gAliOutputFile="";

      MakeGainStore();
      cout << "MUONTRKda : Gain file (to SHUTTLE)  : " << gAliOutputFile << endl;   
    }


  if(fes) // Store IN FES
    {
      printf("\n *****  STORE FILE in FES ****** \n");

      // be sure that env variable DAQDALIB_PATH is set in script file
      //       gSystem->Setenv("DAQDALIB_PATH", "$DATE_SITE/infoLogger");

      if (!gAliOutputFile.IsNull()) 
	{
	  if (gAliCommand.CompareTo("ped") == 0)
	    status = daqDA_FES_storeFile(gAliOutputFile.Data(),"PEDESTALS");
	  else
	    status = daqDA_FES_storeFile(gAliOutputFile.Data(),"GAINS");

	  if (status) 
	    {
	      printf(" Failed to export file : %d\n",status);
	    }
	  else if(gAliPrintLevel) printf(" %s successfully exported to FES  \n",gAliOutputFile.Data());
	}
    }

  gAlifilcout.close();

  printf("\nExecution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return status;
}

