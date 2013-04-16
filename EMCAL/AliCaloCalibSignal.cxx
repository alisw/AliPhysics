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
/* $Id: AliCaloCalibSignal.cxx $ */

//________________________________________________________________________
//
// A help class for monitoring and calibration tools: MOOD, AMORE etc.,
// It can be created and used a la (ctor):
/*
  //Create the object for making the histograms
  fSignals = new AliCaloCalibSignal( fDetType );
  // AliCaloCalibSignal knows how many modules we have for PHOS or EMCAL
  fNumModules = fSignals->GetModules();
*/
// fed an event:
//  fSignals->ProcessEvent(fCaloRawStream,fRawEventHeaderBase);
// get some info:
//  fSignals->GetXXX..()
// etc.
//________________________________________________________________________
#include <string>
#include <sstream>
#include <fstream>

#include "TProfile.h"
#include "TFile.h"

#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"

#include "AliCaloConstants.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliCaloRawAnalyzer.h"
#include "AliCaloRawAnalyzerFactory.h"

//The include file
#include "AliCaloCalibSignal.h"

ClassImp(AliCaloCalibSignal)

using namespace std;

// variables for TTree filling; not sure if they should be static or not
static int fChannelNum = 0; // for regular towers
static int fRefNum = 0; // for LED
static double fAmp = 0;
static double fAvgAmp = 0;
static double fRMS = 0;

// ctor; initialize everything in order to avoid compiler warnings
// put some reasonable defaults
AliCaloCalibSignal::AliCaloCalibSignal(kDetType detectorType) :  
  TObject(),
  fDetType(kNone),
  fColumns(0),
  fRows(0),
  fLEDRefs(0),
  fModules(0),
  fCaloString(),
  fMapping(NULL),
  fFittingAlgorithm(0),  
  fRawAnalyzer(0),
  fRunNumber(-1),
  fStartTime(0),
  fAmpCut(40), // min. 40 ADC counts as default
  fReqFractionAboveAmpCutVal(0.6), // 60% in a strip, per default
  fReqFractionAboveAmp(kTRUE),
  fAmpCutLEDRef(100), // min. 100 ADC counts as default
  fReqLEDRefAboveAmpCutVal(kTRUE),
  fHour(0),
  fLatestHour(0),
  fUseAverage(kTRUE),
  fSecInAverage(1800), 
  fDownscale(10), 
  fNEvents(0),
  fNAcceptedEvents(0),
  fTreeAmpVsTime(NULL),
  fTreeAvgAmpVsTime(NULL),
  fTreeLEDAmpVsTime(NULL),
  fTreeLEDAvgAmpVsTime(NULL)
{
  //Default constructor. First we set the detector-type related constants.
  if (detectorType == kPhos) {
    fColumns = fgkPhosCols;
    fRows = fgkPhosRows;
    fLEDRefs = fgkPhosLEDRefs;
    fModules = fgkPhosModules;
    fCaloString = "PHOS";
  } 
  else {
    //We'll just trust the enum to keep everything in line, so that if detectorType
    //isn't kPhos then it is kEmCal. Note, however, that this is not necessarily the
    //case, if someone intentionally gives another number
    fColumns = AliEMCALGeoParams::fgkEMCALCols;
    fRows = AliEMCALGeoParams::fgkEMCALRows;
    fLEDRefs = AliEMCALGeoParams::fgkEMCALLEDRefs;
    fModules = AliEMCALGeoParams::fgkEMCALModules;
    fCaloString = "EMCAL";
  }

  fDetType = detectorType;
  SetFittingAlgorithm(Algo::kStandard);
  ResetInfo(); // trees and counters
} 

// dtor
//_____________________________________________________________________
AliCaloCalibSignal::~AliCaloCalibSignal()
{
  DeleteTrees();
}

//_____________________________________________________________________
void AliCaloCalibSignal::DeleteTrees()
{
  // delete what was created in the ctor (TTrees)
  if (fTreeAmpVsTime) delete fTreeAmpVsTime;
  if (fTreeAvgAmpVsTime) delete fTreeAvgAmpVsTime;
  if (fTreeLEDAmpVsTime) delete fTreeLEDAmpVsTime;
  if (fTreeLEDAvgAmpVsTime) delete fTreeLEDAvgAmpVsTime;
  // and reset pointers
  fTreeAmpVsTime = NULL;
  fTreeAvgAmpVsTime = NULL;
  fTreeLEDAmpVsTime = NULL;
  fTreeLEDAvgAmpVsTime = NULL;

  return;
}

// copy ctor
//_____________________________________________________________________
//AliCaloCalibSignal::AliCaloCalibSignal(const AliCaloCalibSignal &sig) :
//  TObject(sig),
//  fDetType(sig.GetDetectorType()),
//  fColumns(sig.GetColumns()),
//  fRows(sig.GetRows()),
//  fLEDRefs(sig.GetLEDRefs()),
//  fModules(sig.GetModules()),
//  fCaloString(sig.GetCaloString()),
//  fMapping(), //! note that we are not copying the map info
//  fRunNumber(sig.GetRunNumber()),
//  fStartTime(sig.GetStartTime()),
//  fAmpCut(sig.GetAmpCut()),
//  fReqFractionAboveAmpCutVal(sig.GetReqFractionAboveAmpCutVal()),
//  fReqFractionAboveAmp(sig.GetReqFractionAboveAmp()),
//  fAmpCutLEDRef(sig.GetAmpCutLEDRef()),
//  fReqLEDRefAboveAmpCutVal(sig.GetReqLEDRefAboveAmpCutVal()),
//  fHour(sig.GetHour()),
//  fLatestHour(sig.GetLatestHour()),
//  fUseAverage(sig.GetUseAverage()),
//  fSecInAverage(sig.GetSecInAverage()),
//  fDownscale(sig.GetDownscale()),
//  fNEvents(sig.GetNEvents()),
//  fNAcceptedEvents(sig.GetNAcceptedEvents()),
//  fTreeAmpVsTime(),
//  fTreeAvgAmpVsTime(),
//  fTreeLEDAmpVsTime(),
//  fTreeLEDAvgAmpVsTime()
//{
//  // also the TTree contents
//  AddInfo(&sig);
//  for (Int_t i = 0; i<fgkMaxTowers; i++) {
//      fNHighGain[i] = sig.fNHighGain[i];
//      fNLowGain[i]  = sig.fNLowGain[i]; 
//  }
//  for (Int_t i = 0; i<(2*fgkMaxRefs); i++) {
//    fNRef[i] = sig.fNRef[i]; 
//  }
//  
//  
//}
//
// assignment operator; use copy ctor to make life easy..
//_____________________________________________________________________
//AliCaloCalibSignal& AliCaloCalibSignal::operator = (const AliCaloCalibSignal &source)
//{
//  // assignment operator; use copy ctor
//  if (&source == this) return *this;
//
//  new (this) AliCaloCalibSignal(source);
//  return *this;
//}

//_____________________________________________________________________
void AliCaloCalibSignal::CreateTrees()
{
  // initialize trees
  // first, regular version
  fTreeAmpVsTime = new TTree("fTreeAmpVsTime","Amplitude vs. Time Tree Variables");

  fTreeAmpVsTime->Branch("fChannelNum", &fChannelNum, "fChannelNum/I");
  fTreeAmpVsTime->Branch("fHour", &fHour, "fHour/D");
  fTreeAmpVsTime->Branch("fAmp", &fAmp, "fAmp/D");

  // then, average version
  fTreeAvgAmpVsTime = new TTree("fTreeAvgAmpVsTime","Average Amplitude vs. Time Tree Variables");

  fTreeAvgAmpVsTime->Branch("fChannelNum", &fChannelNum, "fChannelNum/I");
  fTreeAvgAmpVsTime->Branch("fHour", &fHour, "fHour/D");
  fTreeAvgAmpVsTime->Branch("fAvgAmp", &fAvgAmp, "fAvgAmp/D");
  fTreeAvgAmpVsTime->Branch("fRMS", &fRMS, "fRMS/D");

  // then same for LED..
  fTreeLEDAmpVsTime = new TTree("fTreeLEDAmpVsTime","LED Amplitude vs. Time Tree Variables");
  fTreeLEDAmpVsTime->Branch("fRefNum", &fRefNum, "fRefNum/I");
  fTreeLEDAmpVsTime->Branch("fHour", &fHour, "fHour/D");
  fTreeLEDAmpVsTime->Branch("fAmp", &fAmp, "fAmp/D");

  fTreeLEDAvgAmpVsTime = new TTree("fTreeLEDAvgAmpVsTime","Average LED Amplitude vs. Time Tree Variables");
  fTreeLEDAvgAmpVsTime->Branch("fRefNum", &fRefNum, "fRefNum/I");
  fTreeLEDAvgAmpVsTime->Branch("fHour", &fHour, "fHour/D");
  fTreeLEDAvgAmpVsTime->Branch("fAvgAmp", &fAvgAmp, "fAvgAmp/D");
  fTreeLEDAvgAmpVsTime->Branch("fRMS", &fRMS, "fRMS/D");

  return;
}

//_____________________________________________________________________
void AliCaloCalibSignal::ResetInfo()
{ // reset trees and counters
  Zero(); // set all counters to 0
  DeleteTrees(); // delete previous stuff
  CreateTrees(); // and create some new ones
  return;
}

//_____________________________________________________________________
void AliCaloCalibSignal::Zero()
{
  // set all counters to 0; not cuts etc. though
  fHour = 0;
  fLatestHour = 0;
  fNEvents = 0;
  fNAcceptedEvents = 0;

  // Set the number of points for each tower: Amp vs. Time
  memset(fNHighGain, 0, sizeof(fNHighGain));
  memset(fNLowGain, 0, sizeof(fNLowGain));
  // and LED reference
  memset(fNRef, 0, sizeof(fNRef));

  return;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::CheckFractionAboveAmp(const int *iAmpVal, 
						 int resultArray[]) const
{ // check fraction of towers, per column, that are above amplitude cut
  Bool_t returnCode = false;
    
  int iTowerNum = 0;
  double fraction = 0;
  for (int i = 0; i<fModules; i++) {
    for (int j = 0; j<fColumns; j++) {
      int nAbove = 0;
      for (int k = 0; k<fRows; k++) {
	iTowerNum = GetTowerNum(i,j,k);
	if (iAmpVal[iTowerNum] > fAmpCut) { 
	  nAbove++;
	}
      }
      resultArray[i*fColumns +j] = 0; // init. to denied
      if (nAbove > 0) {
	fraction = (1.0*nAbove) / fRows;
	/*
	printf("DS mod %d col %d nAbove %d fraction %3.2f\n",
	       i, j, nAbove, fraction);
	*/
	if (fraction > fReqFractionAboveAmpCutVal) {
	  resultArray[i*fColumns + j] = nAbove;
	  returnCode = true;
	}  
      }
    }
  } // modules loop
  
  return returnCode;
}


//_____________________________________________________________________
Bool_t AliCaloCalibSignal::CheckLEDRefAboveAmp(const int *iAmpVal, 
					       int resultArray[]) const
{ // check which LEDRef/Mon strips are above amplitude cut
  Bool_t returnCode = false;
    
  int iRefNum = 0;
  int gain = 1; // look at high gain; this should be rather saturated usually..
  for (int i = 0; i<fModules; i++) {
    for (int j = 0; j<fLEDRefs; j++) {
      iRefNum = GetRefNum(i, j, gain);
      if (iAmpVal[iRefNum] > fAmpCutLEDRef) { 
	resultArray[i*fLEDRefs +j] = 1; // enough signal
	returnCode = true;
      }
      else {
	resultArray[i*fLEDRefs +j] = 0; // not enough signal
      }
      
      /*
      printf("DS mod %d LEDRef %d ampVal %d\n",
	     i, j, iAmpVal[iRefNum]);
      */
    } // LEDRefs
  } // modules loop
  
  return returnCode;
}

// Parameter/cut handling
//_____________________________________________________________________
void AliCaloCalibSignal::SetParametersFromFile(const char *parameterFile)
{ // set parameters from file
  static const string delimitor("::");
	
  // open, check input file
  ifstream in( parameterFile );
  if( !in ) {
    printf("in AliCaloCalibSignal::SetParametersFromFile - Using default/run_time parameters.\n");
    return;
  } 

  // Note: this method is a bit more complicated than it really has to be
  // - allowing for multiple entries per line, arbitrary order of the
  // different variables etc. But I wanted to try and do this in as
  // correct a C++ way as I could (as an exercise).

  // read in
  char readline[1024];
  while ((in.rdstate() & ios::failbit) == 0 ) {
    
    // Read into the raw char array and then construct a string
    // to do the searching
    in.getline(readline, 1024);
    istringstream s(readline);		
		
    while ( ( s.rdstate() & ios::failbit ) == 0 ) {
			
      string keyValue; 
      s >> keyValue;
      
      // check stream status
      if( ( s.rdstate() & ios::failbit ) == ios::failbit ) break;
			
      // skip rest of line if comments found
      if( keyValue.substr( 0, 2 ) == "//" ) break;
			
      // look for "::" in keyValue pair
      size_t position = keyValue.find( delimitor );
      if( position == string::npos ) {
	printf("wrong format for key::value pair: %s\n", keyValue.c_str());
      }
				
      // split keyValue pair
      string key( keyValue.substr( 0, position ) );
      string value( keyValue.substr( position+delimitor.size(), 
				      keyValue.size()-delimitor.size() ) );
			
      // check value does not contain a new delimitor
      if( value.find( delimitor ) != string::npos ) {
	printf("wrong format for key::value pair: %s\n", keyValue.c_str());
      }
      
      // debug: check key value pair
      // printf("AliCaloCalibSignal::SetParametersFromFile - key %s value %s\n", key.c_str(), value.c_str());

      // if the key matches with something we expect, we assign the new value
      if ( (key == "fAmpCut") || (key == "fReqFractionAboveAmpCutVal") ||
	   (key == "fAmpCutLEDRef") || (key == "fSecInAverage") || 
	   (key == "fFittingAlgorithm") || (key == "fDownscale") ) {
	istringstream iss(value);
	printf("AliCaloCalibSignal::SetParametersFromFile - key %s value %s\n", key.c_str(), value.c_str());

	if (key == "fAmpCut") { 
	  iss >> fAmpCut; 
	}
	else if (key == "fReqFractionAboveAmpCutVal") { 
	  iss >> fReqFractionAboveAmpCutVal; 
	}
	else if (key == "fAmpCutLEDRef") { 
	  iss >> fAmpCutLEDRef; 
	}
	else if (key == "fSecInAverage") { 
	  iss >> fSecInAverage; 
	}
	else if (key == "fFittingAlgorithm") { 
	  iss >> fFittingAlgorithm;
	  SetFittingAlgorithm( fFittingAlgorithm );
	}
	else if (key == "fDownscale") { 
	  iss >> fDownscale; 
	}
      } // some match found/expected

    }		
  }

  in.close();
  return;	
}

//_____________________________________________________________________
void AliCaloCalibSignal::WriteParametersToFile(const char *parameterFile)
{ // write parameters to file
  static const string delimitor("::");
  ofstream out( parameterFile );
  out << "// " << parameterFile << endl;
  out << "fAmpCut" << "::" << fAmpCut << endl;
  out << "fReqFractionAboveAmpCutVal" << "::" << fReqFractionAboveAmpCutVal << endl;
  out << "fAmpCutLEDRef" << "::" << fAmpCutLEDRef << endl;
  out << "fSecInAverage" << "::" << fSecInAverage << endl;
  out << "fFittingAlgorithm" << "::" << fFittingAlgorithm << endl;
  out << "fDownscale" << "::" << fDownscale << endl;

  out.close();
  return;
}

//_____________________________________________________________________
void AliCaloCalibSignal::SetFittingAlgorithm(Int_t fitAlgo)              
{ // select which fitting algo should be used
  fFittingAlgorithm = fitAlgo;
  delete fRawAnalyzer; // delete doesn't do anything if the pointer is 0x0
  fRawAnalyzer = AliCaloRawAnalyzerFactory::CreateAnalyzer( fitAlgo );
  fRawAnalyzer->SetIsZeroSuppressed(true); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::AddInfo(const AliCaloCalibSignal *sig)
{ 
  // note/FIXME: we are not yet adding correctly the info for fN{HighGain,LowGain,Ref} here - but consider this a feature for now (20080905): we'll do Analyze() unless entries were found for a tower in this original object.

  // add info from sig's TTrees to ours..
  TTree *sigAmp = sig->GetTreeAmpVsTime();
  TTree *sigAvgAmp = sig->GetTreeAvgAmpVsTime();

  // we could try some merging via TList or what also as a more elegant approach
  // but I wanted with the stupid/simple and hopefully safe approach of looping
  // over what we want to add..

  // associate variables for sigAmp and sigAvgAmp:
  sigAmp->SetBranchAddress("fChannelNum",&fChannelNum);
  sigAmp->SetBranchAddress("fHour",&fHour);
  sigAmp->SetBranchAddress("fAmp",&fAmp);

  // loop over the trees.. note that since we use the same variables we should not need
  // to do any assignments between the getting and filling
  for (int i=0; i<sigAmp->GetEntries(); i++) {
    sigAmp->GetEntry(i);
    fTreeAmpVsTime->Fill();
  }

  sigAvgAmp->SetBranchAddress("fChannelNum",&fChannelNum);
  sigAvgAmp->SetBranchAddress("fHour",&fHour);
  sigAvgAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  sigAvgAmp->SetBranchAddress("fRMS",&fRMS);

  for (int i=0; i<sigAvgAmp->GetEntries(); i++) {
    sigAvgAmp->GetEntry(i);
    fTreeAvgAmpVsTime->Fill();
  }

  // also LED.. 
  TTree *sigLEDAmp = sig->GetTreeLEDAmpVsTime();
  TTree *sigLEDAvgAmp = sig->GetTreeLEDAvgAmpVsTime();

  // associate variables for sigAmp and sigAvgAmp:
  sigLEDAmp->SetBranchAddress("fRefNum",&fRefNum);
  sigLEDAmp->SetBranchAddress("fHour",&fHour);
  sigLEDAmp->SetBranchAddress("fAmp",&fAmp);

  // loop over the trees.. note that since we use the same variables we should not need
  // to do any assignments between the getting and filling
  for (int i=0; i<sigLEDAmp->GetEntries(); i++) {
    sigLEDAmp->GetEntry(i);
    fTreeLEDAmpVsTime->Fill();
  }

  sigLEDAvgAmp->SetBranchAddress("fRefNum",&fRefNum);
  sigLEDAvgAmp->SetBranchAddress("fHour",&fHour);
  sigLEDAvgAmp->SetBranchAddress("fAvgAmp",&fAvgAmp);
  sigLEDAvgAmp->SetBranchAddress("fRMS",&fRMS);

  for (int i=0; i<sigLEDAvgAmp->GetEntries(); i++) {
    sigLEDAvgAmp->GetEntry(i);
    fTreeLEDAvgAmpVsTime->Fill();
  }

  // We should also copy other pieces of info: counters and parameters 
  // (not number of columns and rows etc which should be the same)
  // note that I just assign them here rather than Add them, but we
  // normally just Add (e.g. in Preprocessor) one object so this should be fine.
  fRunNumber = sig->GetRunNumber();
  fStartTime = sig->GetStartTime();
  fAmpCut = sig->GetAmpCut();
  fReqFractionAboveAmpCutVal = sig->GetReqFractionAboveAmpCutVal();
  fReqFractionAboveAmp = sig->GetReqFractionAboveAmp();
  fAmpCutLEDRef = sig->GetAmpCutLEDRef();
  fReqLEDRefAboveAmpCutVal = sig->GetReqLEDRefAboveAmpCutVal();
  fHour = sig->GetHour();
  fLatestHour = sig->GetLatestHour();
  fUseAverage = sig->GetUseAverage();
  fSecInAverage = sig->GetSecInAverage();
  fDownscale = sig->GetDownscale();
  fNEvents = sig->GetNEvents();
  fNAcceptedEvents = sig->GetNAcceptedEvents();

  return kTRUE;//We hopefully succesfully added info from the supplied object
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::ProcessEvent(AliRawReader *rawReader)
{
  // if fMapping is NULL the rawstream will crate its own mapping
  AliCaloRawStreamV3 rawStream(rawReader, fCaloString, (AliAltroMapping**)fMapping);  
  if (fDetType == kEmCal) {
    rawReader->Select("EMCAL", 0, AliEMCALGeoParams::fgkLastAltroDDL) ; //select EMCAL DDL range 
  }

  return ProcessEvent( &rawStream, rawReader->GetTimestamp() );
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::ProcessEvent(AliCaloRawStreamV3 *in, UInt_t Timestamp)
{ 
  // Method to process=analyze one event in the data stream
  if (!in) return kFALSE; //Return right away if there's a null pointer
  
  fNEvents++; // one more event

  if ( (fNEvents%fDownscale)!=0 ) return kFALSE; // mechanism to skip some of the input events, if we want

  // use maximum numbers to set array sizes
  int iAmpValHighGain[fgkMaxTowers];
  int iAmpValLowGain[fgkMaxTowers];
  memset(iAmpValHighGain, 0, sizeof(iAmpValHighGain));
  memset(iAmpValLowGain, 0, sizeof(iAmpValLowGain));

  // also for LED reference
  int iLEDAmpVal[fgkMaxRefs * 2]; // factor 2 is for the two gain values
  memset(iLEDAmpVal, 0, sizeof(iLEDAmpVal));

  int gain = 0; // high or low gain
  
  // Number of Low and High gain, and LED Ref, channels for this event:
  int nLowChan = 0; 
  int nHighChan = 0; 
  int nLEDRefChan = 0;

  int iTowerNum = 0; // array index for regular towers
  int iRefNum = 0; // array index for LED references

  // loop first to get the fraction of channels with amplitudes above cut

  while (in->NextDDL()) {
    while (in->NextChannel()) {

      vector<AliCaloBunchInfo> bunchlist; 
      while (in->NextBunch()) {
	bunchlist.push_back( AliCaloBunchInfo(in->GetStartTimeBin(), in->GetBunchLength(), in->GetSignals() ) );
      } 
      if (bunchlist.size() == 0) continue;

      gain = -1; // init to not valid value
      //If we're here then we're done with this tower
      if ( in->IsLowGain() ) {
	gain = 0;
      }
      else if ( in->IsHighGain() ) {
	gain = 1;
      }
      else if ( in->IsLEDMonData() ) {
	gain = in->GetRow(); // gain coded in (in RCU/Altro mapping) as Row info for LED refs..
      }
      else { continue; } // don't try to fit TRU..

      // it should be enough to check the SuperModule info for each DDL really, but let's keep it here for now
      int arrayPos = in->GetModule(); //The modules are numbered starting from 0
      //Debug
      if (arrayPos < 0 || arrayPos >= fModules) {
	printf("AliCaloCalibSignal::ProcessEvent = Oh no: arrayPos = %i.\n", arrayPos); 
	return kFALSE;
      }

      AliCaloFitResults res =  fRawAnalyzer->Evaluate( bunchlist, in->GetAltroCFG1(), in->GetAltroCFG2());  
      if ( in->IsHighGain() || in->IsLowGain() ) { // regular tower
	// get tower number for AmpVal array
	iTowerNum = GetTowerNum(arrayPos, in->GetColumn(), in->GetRow()); 

	if (gain == 0) {
	  // fill amplitude into the array	   
	  iAmpValLowGain[iTowerNum]  = (int) res.GetAmp();
	  nLowChan++;
	} 
	else if (gain==1) {//fill the high gain ones
	  // fill amplitude into the array
	  iAmpValHighGain[iTowerNum] = (int) res.GetAmp();
	  nHighChan++;
	}//end if gain
      } // regular tower
      else if ( in->IsLEDMonData() ) { // LED ref.; 
	// strip # is coded is 'column' in the channel maps 
	iRefNum = GetRefNum(arrayPos, in->GetColumn(), gain); 
	iLEDAmpVal[iRefNum] = (int) res.GetAmp();
	nLEDRefChan++;
      } // end of LED ref

    } // end while over channel 
   
  }//end while over DDL's, of input stream
  
  in->Reset(); // just in case the next customer forgets to check if the stream was reset..

  // now check if it was an LED event, using the LED Reference info per strip

  // by default all columns are accepted (init check to > 0)
  int checkResultArray[AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALCols];
  for (int ia=0; ia<(AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALCols); ia++) { 
    checkResultArray[ia] = 1; 
  }
  if (fReqFractionAboveAmp) {
    bool ok = false;
    if (nHighChan > 0) { 
      ok = CheckFractionAboveAmp(iAmpValHighGain, checkResultArray); 
    }
    if (!ok) return false; // skip event
  }

  // by default all columns are accepted (init check to > 0)
  int checkResultArrayLEDRef[AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALLEDRefs];
  for (int ia=0; ia<(AliEMCALGeoParams::fgkEMCALModules * AliEMCALGeoParams::fgkEMCALLEDRefs); ia++) { 
    checkResultArrayLEDRef[ia] = 1; 
  }
  if (fReqLEDRefAboveAmpCutVal) {
    bool ok = false;
    if (nLEDRefChan > 0) { 
      ok = CheckLEDRefAboveAmp(iLEDAmpVal, checkResultArrayLEDRef); 
    }
    if (!ok) return false; // skip event
  }

  fNAcceptedEvents++; // one more event accepted

  if (fStartTime == 0) { // if start-timestamp wasn't set,we'll pick it up from the first event we encounter
    fStartTime = Timestamp;
  }

  fHour = (Timestamp - fStartTime)/(double)fgkNumSecInHr;
  if (fLatestHour < fHour) {
    fLatestHour = fHour; 
  }
  
  // it is a led event, now fill TTree
  // We also do the activity check for LEDRefs/Strips, but need to translate between column
  // and strip indices for that; based on these relations: 
  // iStrip = AliEMCALGeoParams::GetStripModule(iSM, iCol);
  // iStrip = (iSM%2==0) ? iCol/2 : AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - iCol/2;
  // which leads to
  // iColFirst = (iSM%2==0) ? iStrip*2 : (AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - iStrip)*2;

  for(int i=0; i<fModules; i++){
    for(int j=0; j<fColumns; j++) {
      int iStrip = (i%2==0) ? j/2 : AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - j/2;
      if (checkResultArray[i*fColumns + j]>0  && checkResultArrayLEDRef[i*fLEDRefs + iStrip]>0) { // column passed check 
      for(int k=0; k<fRows; k++){
	
	iTowerNum = GetTowerNum(i, j, k); 

	if(iAmpValHighGain[iTowerNum]) {
	  fAmp = iAmpValHighGain[iTowerNum];
	  fChannelNum = GetChannelNum(i,j,k,1);
	  fTreeAmpVsTime->Fill();//fChannelNum,fHour,AmpValHighGain[iTowerNum]);
	  fNHighGain[iTowerNum]++;
	}
	if(iAmpValLowGain[iTowerNum]) {
	  fAmp = iAmpValLowGain[iTowerNum];
	  fChannelNum = GetChannelNum(i,j,k,0);
	  fTreeAmpVsTime->Fill();//fChannelNum,fHour,AmpValLowGain[iTowerNum]);
	  fNLowGain[iTowerNum]++;
	}
      } // rows
      } // column passed check, and LED Ref for strip passed check (if any)
    } // columns

    // also LED refs
    for(int j=0; j<fLEDRefs; j++){
      int iColFirst = (i%2==0) ? j*2 : (AliEMCALGeoParams::fgkEMCALLEDRefs - 1 - j)*2; //CHECKME!!!
      if ( ((checkResultArray[i*fColumns + iColFirst]>0) || (checkResultArray[i*fColumns + iColFirst + 1]>0)) && // at least one column in strip passed check 
	   (checkResultArrayLEDRef[i*fLEDRefs + j]>0) ) { // and LED Ref passed checks
	for (gain=0; gain<2; gain++) {
	  fRefNum = GetRefNum(i, j, gain); 
	  if (iLEDAmpVal[fRefNum]) {
	    fAmp = iLEDAmpVal[fRefNum];
	    fTreeLEDAmpVsTime->Fill();//fRefNum,fHour,fAmp);
	    fNRef[fRefNum]++;
	  }
	} // gain
      } // at least one column in strip passed check, and LED Ref passed check (if any) 
    }

  } // modules
  
  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::Save(TString fileName)
{
  //Saves all the TTrees to the designated file
  
  TFile destFile(fileName, "recreate");
  
  if (destFile.IsZombie()) {
    return kFALSE;
  }
  
  destFile.cd();

  // save the trees
  fTreeAmpVsTime->Write();
  fTreeLEDAmpVsTime->Write();
  if (fUseAverage) { 
    Analyze(); // get the latest and greatest averages
    fTreeAvgAmpVsTime->Write();
    fTreeLEDAvgAmpVsTime->Write();
  }

  destFile.Close();
  
  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::Analyze()
{
  // Fill the tree holding the average values
  if (!fUseAverage) { return kFALSE; }

  // Reset the average TTree if Analyze has already been called earlier,
  // meaning that the TTree could have been partially filled
  if (fTreeAvgAmpVsTime->GetEntries() > 0) {
    fTreeAvgAmpVsTime->Reset();
  }

  //0: setup variables for the TProfile plots that we'll use to do the averages
  int numProfBins = 0;
  double timeMin = 0;
  double timeMax = 0;
  if (fSecInAverage > 0) {
    numProfBins = (int)( (fLatestHour*fgkNumSecInHr)/fSecInAverage + 1 ); // round-off
  }
  numProfBins += 2; // add extra buffer : first and last
  double binSize = 1.0*fSecInAverage / fgkNumSecInHr;
  timeMin = - binSize;
  timeMax = timeMin + numProfBins*binSize;

  //1: set up TProfiles for the towers that had data
  TProfile * profile[fgkMaxTowers*2]; // *2 is since we include both high and low gains
  memset(profile, 0, sizeof(profile));
  const Int_t buffersize = 200;
  char name[buffersize]; // for profile id and title
  int iTowerNum = 0;

  for (int i = 0; i<fModules; i++) {
    for (int ic=0; ic<fColumns; ic++){
      for (int ir=0; ir<fRows; ir++) {

	iTowerNum = GetTowerNum(i, ic, ir);
	// high gain
	if (fNHighGain[iTowerNum] > 0) {
	  fChannelNum = GetChannelNum(i, ic, ir, 1); 
	  snprintf(name,buffersize,"profileChan%d", fChannelNum);
	  profile[fChannelNum] = new TProfile(name, name, numProfBins, timeMin, timeMax, "s");
	}

	// same for low gain
	if (fNLowGain[iTowerNum] > 0) {
	  fChannelNum = GetChannelNum(i, ic, ir, 0); 
	  snprintf(name,buffersize,"profileChan%d", fChannelNum);
	  profile[fChannelNum] = new TProfile(name, name, numProfBins, timeMin, timeMax, "s");
	}

      } // rows
    } // columns
  } // modules

  //2: fill profiles by looping over tree
  // Set addresses for tree-readback also
  fTreeAmpVsTime->SetBranchAddress("fChannelNum", &fChannelNum);
  fTreeAmpVsTime->SetBranchAddress("fHour", &fHour);
  fTreeAmpVsTime->SetBranchAddress("fAmp", &fAmp);

  for (int ient=0; ient<fTreeAmpVsTime->GetEntries(); ient++) {
    fTreeAmpVsTime->GetEntry(ient);
    if (profile[fChannelNum]) { 
      // profile should always have been created above, for active channels
      profile[fChannelNum]->Fill(fHour, fAmp);
    }
  }

  // re-associating the branch addresses here seems to be needed for OK 'average' storage	    
  fTreeAvgAmpVsTime->SetBranchAddress("fChannelNum", &fChannelNum);
  fTreeAvgAmpVsTime->SetBranchAddress("fHour", &fHour);
  fTreeAvgAmpVsTime->SetBranchAddress("fAvgAmp", &fAvgAmp);
  fTreeAvgAmpVsTime->SetBranchAddress("fRMS", &fRMS);

  //3: fill avg tree by looping over the profiles
  for (fChannelNum = 0; fChannelNum<(fgkMaxTowers*2); fChannelNum++) {
    if (profile[fChannelNum]) { // profile was created
      if (profile[fChannelNum]->GetEntries() > 0) { // profile had some entries
	for(int it=0; it<numProfBins; it++) {
	  if (profile[fChannelNum]->GetBinEntries(it+1) > 0) {
	    fAvgAmp = profile[fChannelNum]->GetBinContent(it+1);
	    fHour = profile[fChannelNum]->GetBinCenter(it+1);
	    fRMS = profile[fChannelNum]->GetBinError(it+1);
	    fTreeAvgAmpVsTime->Fill();
	  } // some entries for this bin
	} // loop over bins
      } // some entries for this profile
    } // profile exists  
  } // loop over all possible channels


  // and finally, go through same exercise for LED also.. 

  //1: set up TProfiles for the towers that had data
  TProfile * profileLED[fgkMaxRefs*2]; // *2 is since we include both high and low gains
  memset(profileLED, 0, sizeof(profileLED));

  for (int i = 0; i<fModules; i++) {
    for(int j=0; j<fLEDRefs; j++){
      for (int gain=0; gain<2; gain++) {
	fRefNum = GetRefNum(i, j, gain);
	if (fNRef[fRefNum] > 0) { 
	  snprintf(name, buffersize, "profileLEDRef%d", fRefNum);
	  profileLED[fRefNum] = new TProfile(name, name, numProfBins, timeMin, timeMax, "s");
	} 
      }// gain
    } 
  } // modules

  //2: fill profiles by looping over tree
  // Set addresses for tree-readback also
  fTreeLEDAmpVsTime->SetBranchAddress("fRefNum", &fRefNum);
  fTreeLEDAmpVsTime->SetBranchAddress("fHour", &fHour);
  fTreeLEDAmpVsTime->SetBranchAddress("fAmp", &fAmp);

  for (int ient=0; ient<fTreeLEDAmpVsTime->GetEntries(); ient++) {
    fTreeLEDAmpVsTime->GetEntry(ient);
    if (profileLED[fRefNum]) { 
      // profile should always have been created above, for active channels
      profileLED[fRefNum]->Fill(fHour, fAmp);
    }
  }

  // re-associating the branch addresses here seems to be needed for OK 'average' storage	    
  fTreeLEDAvgAmpVsTime->SetBranchAddress("fRefNum", &fRefNum);
  fTreeLEDAvgAmpVsTime->SetBranchAddress("fHour", &fHour);
  fTreeLEDAvgAmpVsTime->SetBranchAddress("fAvgAmp", &fAvgAmp);
  fTreeLEDAvgAmpVsTime->SetBranchAddress("fRMS", &fRMS);

  //3: fill avg tree by looping over the profiles
  for (fRefNum = 0; fRefNum<(fgkMaxRefs*2); fRefNum++) {
    if (profileLED[fRefNum]) { // profile was created
      if (profileLED[fRefNum]->GetEntries() > 0) { // profile had some entries
	for(int it=0; it<numProfBins; it++) {
	  if (profileLED[fRefNum]->GetBinEntries(it+1) > 0) {
	    fAvgAmp = profileLED[fRefNum]->GetBinContent(it+1);
	    fHour = profileLED[fRefNum]->GetBinCenter(it+1);
	    fRMS = profileLED[fRefNum]->GetBinError(it+1);
	    fTreeLEDAvgAmpVsTime->Fill();
	  } // some entries for this bin
	} // loop over bins
      } // some entries for this profile
    } // profile exists  
  } // loop over all possible channels

  // OK, we're done..

  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::DecodeChannelNum(const int chanId, 
					    int *imod, int *icol, int *irow, int *igain) const  
{ // return the module, column, row, and gain for a given channel number
  *igain = chanId/(fModules*fColumns*fRows);
  *imod = (chanId/(fColumns*fRows)) % fModules;
  *icol = (chanId/fRows) % fColumns;
  *irow = chanId % fRows;
  return kTRUE;
} 

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::DecodeRefNum(const int refId, 
					int *imod, int *istripMod, int *igain) const 
{ // return the module, stripModule, and gain for a given reference number
  *igain = refId/(fModules*fLEDRefs);
  *imod = (refId/(fLEDRefs)) % fModules;
  *istripMod = refId % fLEDRefs;
  return kTRUE;
} 
