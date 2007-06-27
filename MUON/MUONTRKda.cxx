/*

Version 1 for MUONTRKda MUON tracking
Working version for pedestal
Framework for gain computing
(Ch. Finck)


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


// global variables
AliMUONVStore* pedestalStore =  new AliMUON2DMap(kFALSE);
const Int_t kNchannels = AliMpConstants::ManuNofChannels();
Int_t nManu = 0;
Int_t nChannel = 0;
UInt_t runNumber = 0;
Int_t nEvents = 0;
TH1F* pedMeanHisto = 0x0;
TH1F* pedSigmaHisto = 0x0;
Char_t histoFileName[256];
TString command("ped");

// funtions

//________________
Double_t fitFunc(Double_t *x, Double_t *par)
{
    //fit function for gains

    Double_t xx = x[0];
    Double_t f;

    if (xx < par[3])
	f = par[0] + par[1]*xx;
    else
	f= par[0] + par[1]*xx + par[2]*xx*xx;

    return f;
}



//__________
void MakePed(Int_t busPatchId, Int_t manuId, Int_t channelId, Int_t charge)
{

    AliMUONVCalibParam* ped = 
	static_cast<AliMUONVCalibParam*>(pedestalStore->FindObject(busPatchId, manuId));

    if (!ped) {
      nManu++;
      ped = new AliMUONCalibParamND(2, kNchannels,busPatchId, manuId, -1.); // put default wise -1, not connected channel
      pedestalStore->Add(ped);	
    }

    if (nEvents == 1) {
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
  sprintf(histoFileName,"mutrkped-%d.root",runNumber);
  TFile*  histoFile = new TFile(histoFileName,"RECREATE","MUON Tracking pedestals");

  Char_t name[255];
  Char_t title[255];
  sprintf(name,"pedmean_allch");
  sprintf(title,"Pedestal mean all channels");
  Int_t nx = 1000;
  Int_t xmin = 0;
  Int_t xmax = 1000; 
  pedMeanHisto = new TH1F(name,title,nx,xmin,xmax);
  pedMeanHisto->SetDirectory(histoFile);

  sprintf(name,"pedsigma_allch");
  sprintf(title,"Pedestal sigma all channels");
  nx = 200;
  xmin = 0;
  xmax = 50; 
  pedSigmaHisto = new TH1F(name,title,nx,xmin,xmax);
  pedSigmaHisto->SetDirectory(histoFile);
    
  TTree* tree = new TTree("t","Pedestal tree");
  tree->Branch("bp",&busPatchId,"bp/I");
  tree->Branch("manu",&manuId,",manu/I");
  tree->Branch("channel",&channelId,",channel/I");

  if (!flatOutputFile.IsNull()) {
    fileout.open(flatOutputFile.Data());
    fileout<<"//===========================================================================" << endl;
    fileout<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
    fileout<<"//===========================================================================" << endl;
    fileout<<"//       * Run           : " << runNumber << endl; 
    fileout<<"//       * Date          : " << date.AsString("l") <<endl;
    fileout<<"//       * Statictics    : " << nEvents << endl;
    fileout<<"//       * # of MANUS    : " << nManu << endl;
    fileout<<"//       * # of channels : " << nChannel << endl;
    fileout<<"//"<<endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;
    fileout<<"//format : BUS_PATCH MANU_ID CHANNEL MEAN SIGMA"<<endl;
    fileout<<"//---------------------------------------------------------------------------" << endl;

  }

  // iterator over pedestal
  TIter next(pedestalStore->CreateIterator());
  AliMUONVCalibParam* ped;
  
  while ( ( ped = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
  {
    busPatchId              = ped->ID0();
    manuId                  = ped->ID1();

    for (channelId = 0; channelId < ped->Size() ; ++channelId) {

      pedMean  = ped->ValueAsDouble(channelId, 0);

      if (pedMean > 0) { // connected channels

	ped->SetValueAsDouble(channelId, 0, pedMean/(Double_t)nEvents);

	pedMean  = ped->ValueAsDouble(channelId, 0);
	pedSigma = ped->ValueAsDouble(channelId, 1);

	ped->SetValueAsDouble(channelId, 1, TMath::Sqrt(TMath::Abs(pedSigma/(Double_t)nEvents - pedMean*pedMean)));

	pedMean  = ped->ValueAsDouble(channelId, 0) + 0.5 ;
	pedSigma = ped->ValueAsDouble(channelId, 1);


	if (!flatOutputFile.IsNull()) {
	  fileout << "\t" << busPatchId << "\t" << manuId <<"\t"<< channelId << "\t"
		  << pedMean <<"\t"<< pedSigma << endl;
	}

	pedMeanHisto->Fill(pedMean);
	pedSigmaHisto->Fill(pedSigma);

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
void MakePedStoreForGain()
{
    // store pedestal map in root file

    Int_t injCharge = 200;

    TTree* tree = 0x0;

    // compute and store pedestals
    MakePedStore();

    // store in root file
    sprintf(histoFileName,"mutrkgain.root");
    
    TString mode("UPDATE");

    if (command.Contains("cre")) {
	mode = "RECREATE";
    }
    TFile* histoFile = new TFile(histoFileName, mode.Data(), "MUON Tracking Gains");

    // second argument should be the injected charge, taken from config crocus file
    // put also info about run number could be usefull
    AliMpIntPair* pair   = new AliMpIntPair(runNumber, injCharge);

    if (mode.CompareTo("UPDATE") == 0) {
      tree = (TTree*)histoFile->Get("t");
      tree->SetBranchAddress("run",&pair);
      tree->SetBranchAddress("ped",&pedestalStore);

    } else {
      tree = new TTree("t","Pedestal tree");
      tree->Branch("run", "AliMpIntPair",&pair);
      tree->Branch("ped", "AliMUON2DMap",&pedestalStore);
      tree->SetBranchAddress("run",&pair);
      tree->SetBranchAddress("ped",&pedestalStore);

    }

    tree->Fill();
    tree->Write("t", TObject::kOverwrite); // overwrite the tree
    histoFile->Close();

    delete pair;
}

//________________
void MakeGainStore(TString flatOutputFile)
{
    
    // open file mutrkgain.root
    // read again the pedestal for the calibration runs (9 runs ?)
    // need the injection charge from config file (to be done)
    // For each channel make a TGraphErrors (mean, sigma) vs injected charge
    // Fit with a polynomial fct
    // store the result in a flat file.

    Double_t pedMean[10];
    Double_t pedSigma[10];
    Double_t injCharge[10];
    Double_t injChargeErr[10];

    ofstream fileout;
    TTimeStamp date;

    if (!flatOutputFile.IsNull()) {
      fileout.open(flatOutputFile.Data());
      fileout<<"//===========================================================================" << endl;
      fileout<<"//                       Pedestal file calculated by MUONTRKda"<<endl;
      fileout<<"//===========================================================================" << endl;
      fileout<<"//       * Run           : " << runNumber << endl; 
      fileout<<"//       * Date          : " << date.AsString("l") <<endl;
      fileout<<"//       * Statictics    : " << nEvents << endl;
      fileout<<"//       * # of MANUS    : " << nManu << endl;
      fileout<<"//       * # of channels : " << nChannel << endl;
      fileout<<"//"<<endl;
      fileout<<"//---------------------------------------------------------------------------" << endl;
      fileout<<"//---------------------------------------------------------------------------" << endl;
      fileout<<"//format : BUS_PATCH MANU_ID CHANNEL GAIN0 GAIN1 GAIN2 XLIM CHI2"<<endl;
      fileout<<"//---------------------------------------------------------------------------" << endl;

    }


    sprintf(histoFileName,"mutrkgain.root");
    TFile*  histoFile = new TFile(histoFileName);

    AliMUON2DMap* map[10];
    AliMUONVCalibParam* ped[10];
    //   AliMpIntPair* pair;
    AliMpIntPair* run[10];

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
    }
      
    // Q = f(ADC)
    TF1* func = new TF1("func",fitFunc, 0., 2500., 4);
    // TF1* func = new TF1("func","pol1");

    // iterates over the first pedestal run
    TIter next(map[0]->CreateIterator());
    AliMUONVCalibParam* p;

    while ( ( p = dynamic_cast<AliMUONVCalibParam*>(next() ) ) )
    {
      ped[0]  = p;

      Int_t busPatchId = p->ID0();
      Int_t manuId     = p->ID1();

      // read back pedestal from the other runs for the given (bupatch, manu)
      for (Int_t i = 1; i < nEntries; ++i) {
	ped[i] = static_cast<AliMUONVCalibParam*>(map[i]->FindObject(busPatchId, manuId));
      }

      // compute for each channel the gain parameters
      for ( Int_t channelId = 0; channelId < ped[0]->Size() ; ++channelId ) {

	Int_t n = 0;
	for (Int_t i = 0; i < nEntries; ++i) {

	  if (!ped[i]) continue; //shouldn't happen.
	  pedMean[i]      = ped[i]->ValueAsDouble(channelId, 0);
	  pedSigma[i]     = ped[i]->ValueAsDouble(channelId, 1);
	  injCharge[i]    = (Double_t)run[i]->GetSecond();
	  injChargeErr[i] = 1.;

	  if (pedMean[i] < 0) continue; // not connected

	  if (pedSigma[i] <= 0) pedSigma[i] = 1.; // should not happen.
	  n++;
	}

	if (n > 4) {
	  // if (n > 1) {
	  //fit 
	  TGraph *gain = new TGraphErrors(n, pedMean, injCharge, pedSigma, injChargeErr);
	  //should set some initial parameters
	  func->SetParameter(0,-300);  // a0
	  func->SetParameter(1, 1.);   // a1
	  func->SetParameter(2, 0.00001);// a2
	  func->SetParameter(3, 1100.); // xlim in ADC

	  gain->Fit("func","q");
	  cout	<< setw(8) << func->GetParameter(0)  <<"\t"<< setw(8) << func->GetParameter(1) << endl;

	  if (!flatOutputFile.IsNull()) {
	    fileout << "\t" << busPatchId << "\t" << manuId <<"\t"<< channelId << "\t"
		    << setw(8) << func->GetParameter(0)  <<"\t"<< setw(8) << func->GetParameter(1) 
		    << setw(8) << func->GetParameter(2)  <<"\t"<< setw(5) << func->GetParameter(3) 
		    << "\t" << func->GetChisquare() << endl;
	  }
	  delete gain;
	}

      }
    }

    // file outputs for gain
    if (!flatOutputFile.IsNull()) 
	fileout.close();

    delete func;

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


    Int_t printLevel = 0;  // Global variable defined as extern in the others .cxx files
    Int_t skipEvents = 0;
    Int_t maxEvents  = 1000000;
    Double_t nSigma = 3;
    Int_t threshold = -1;
    Char_t inputFile[256];
    TString flatOutputFile;
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
	  command = argv[i];
	  break;
	case 'd' :
	  i++; 
	  printLevel=atoi(argv[i]);
	  break;
	case 's' :
	  i++; 
	  skipEvents=atoi(argv[i]);
	  break;
	case 'n' :
	  i++; 
	  sscanf(argv[i],"%d",&maxEvents);
	  break;
	case 'p' :
	  i++; 
	  sscanf(argv[i],"%lf",&nSigma);
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
	  printf("\n-o <CROUCUS cmd file>     (default = %s)",crocusOutputFile.Data()); 
	  printf("\n");
	  printf("\n Options");
	  printf("\n-d <print level>          (default = %d)",printLevel);
	  printf("\n-s <skip events>          (default = %d)",skipEvents);
	  printf("\n-n <max events>           (default = %d)",maxEvents);
	  printf("\n-p <n sigmas>             (default = %f)",nSigma);
	  printf("\n-t <threshold (-1 = no)>  (default = %d)",threshold);
	  printf("\n-e <execute ped/gain>     (default = %s)",command.Data());
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

  // set command to lower case
  command.ToLower();

  // decoding the events
  
  Int_t status;
  Int_t nDateEvents = 0;

  void* event;

  pedMeanHisto = 0x0;
  pedSigmaHisto = 0x0;

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

  cout << "MUONTRKda : Reading data from file " << inputFile <<endl;

  Int_t busPatchId;
  UShort_t manuId;  
  UChar_t channelId;
  UShort_t charge;

  while(1) 
  {
    if (nEvents >= maxEvents) break;
    if (nEvents && nEvents % 100 == 0) 	
	cout<<"Cumulated events " << nEvents << endl;

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

    nDateEvents++;

    // decoding rawdata headers
    AliRawReader *rawReader = new AliRawReaderDate(event);
 
    Int_t eventType = rawReader->GetType();
    runNumber = rawReader->GetRunNumber();

    if (eventType != PHYSICS_EVENT)
	continue; // for the moment

    nEvents++;

    // decoding MUON payload
    AliMUONRawStreamTracker* rawStream  = new AliMUONRawStreamTracker(rawReader);

    // loops over DDL 
    rawStream->First();
    while( (status = rawStream->Next(busPatchId, manuId, channelId, charge)) ) {
  
      if (nEvents == 1)
	  nChannel++;
      
      if (printLevel) printf("manuId: %d, channelId: %d charge: %d\n", manuId, 
			     channelId, charge);
      
      MakePed(busPatchId, (Int_t)manuId, (Int_t)channelId, (Int_t)charge);
		  
    } // Next digit
     
    delete rawReader;
    delete rawStream;

  } // while (1)


  if (command.CompareTo("ped") == 0)
      MakePedStore(flatOutputFile);

  // option gain -> update root file with pedestal results
  // gain + create -> recreate root file
  // gain + comp -> update root file and compute gain parameters

  if (command.Contains("gain")) 
      MakePedStoreForGain();
  
  if (command.Contains("comp")) 
      MakeGainStore(flatOutputFile);
  

  delete pedestalStore;

  timers.Stop();

  if (!(crocusConfigFile.IsNull()))
      cout << "MUONTRKda : CROCUS command file generated : " << crocusOutputFile.Data() << endl;
  else
      cout << "MUONTRKda : WARNING no CROCUS command file generated" << endl;

  cout << "MUONTRKda : Flat ASCII file generated     : " << flatOutputFile << endl;
  cout << "MUONTRKda : Histo file generated          : " << histoFileName  << endl;
  cout << "MUONTRKda : Nb of DATE events     = "         << nDateEvents    << endl;
  cout << "MUONTRKda : Nb of events used     = "         << nEvents        << endl;

  printf("Execution time : R:%7.2fs C:%7.2fs\n", timers.RealTime(), timers.CpuTime());

  return status;
}
