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
//* $Id$ */

//________________________________________________________________________
//
// A help class for monitoring and calibration tools: MOOD, AMORE etc.,
// It can be created and used a la (ctor):
/*
  //Create the object for making the histograms
  fPedestals = new AliCaloCalibPedestal( fDetType );
  // AliCaloCalibPedestal knows how many modules we have for PHOS or EMCAL
  fNumModules = fPedestals->GetModules();
*/
// fed an event:
//  fPedestals->ProcessEvent(fCaloRawStream);
// asked to draw histograms:
//  fPedestals->GetDeadMap(i)->Draw("col");
// or
//  fPedestals->GetPeakProfileHighGainRatio((i < fNumModules) ? i : fVisibleModule)->Draw("colz");
// etc.
// The pseudo-code examples above were from the first implementation in MOOD (summer 2007).
//________________________________________________________________________

//#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <cmath>

#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"

//The include file
#include "AliCaloCalibPedestal.h"

ClassImp(AliCaloCalibPedestal)

using namespace std;

// ctor; initialize everything in order to avoid compiler warnings
AliCaloCalibPedestal::AliCaloCalibPedestal(kDetType detectorType) :  
  TObject(),
  fPedestalLowGain(),
  fPedestalHighGain(),
  fPedestalLEDRefLowGain(),
  fPedestalLEDRefHighGain(),
  fPeakMinusPedLowGain(),
  fPeakMinusPedHighGain(),
  fPeakMinusPedHighGainHisto(),
  fPedestalLowGainDiff(),
  fPedestalHighGainDiff(),
  fPedestalLEDRefLowGainDiff(),
  fPedestalLEDRefHighGainDiff(),
  fPeakMinusPedLowGainDiff(),
  fPeakMinusPedHighGainDiff(),
  fPedestalLowGainRatio(),
  fPedestalHighGainRatio(),
  fPedestalLEDRefLowGainRatio(),
  fPedestalLEDRefHighGainRatio(),
  fPeakMinusPedLowGainRatio(),
  fPeakMinusPedHighGainRatio(),
  fDeadMap(),
  fNEvents(0),
  fNChanFills(0),
  fDeadTowers(0),
  fNewDeadTowers(0),
  fResurrectedTowers(0),
  fReference(0),
  fDetType(kNone),
  fColumns(0),
  fRows(0),
  fLEDRefs(0),
  fModules(0),
  fRowMin(0),
  fRowMax(0),
  fRowMultiplier(0),
  fCaloString(),
  fMapping(NULL),
  fRunNumber(-1),
  fSelectPedestalSamples(kTRUE), 
  fFirstPedestalSample(0),
  fLastPedestalSample(15),
  fDeadThreshold(5),
  fWarningThreshold(50),
  fWarningFraction(0.002),
  fHotSigma(5)
{
  //Default constructor. First we set the detector-type related constants.
  if (detectorType == kPhos) {
    fColumns = fgkPhosCols;
    fRows = fgkPhosRows;
    fLEDRefs = fgkPhosLEDRefs;
    fModules = fgkPhosModules;
    fCaloString = "PHOS";
    fRowMin = -1*fRows;
    fRowMax = 0;
    fRowMultiplier = -1;
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
    fRowMin = 0;
    fRowMax = fRows;
    fRowMultiplier = 1;
  } 
  fDetType = detectorType;
 
  //Then, loop for the requested number of modules
  TString title, name;
  for (int i = 0; i < fModules; i++) {
    //Pedestals, low gain
    name = "hPedlowgain";
    name += i;
    title = "Pedestals, low gain, module ";
    title += i; 
    fPedestalLowGain.Add(new TProfile2D(name, title,
					fColumns, 0.0, fColumns, 
					fRows, fRowMin, fRowMax,"s"));
  
    //Pedestals, high gain
    name = "hPedhighgain";
    name += i;
    title = "Pedestals, high gain, module ";
    title += i; 
    fPedestalHighGain.Add(new TProfile2D(name, title,
					 fColumns, 0.0, fColumns, 
					 fRows, fRowMin, fRowMax,"s"));

    //LED Ref/Mon pedestals, low gain
    name = "hPedestalLEDReflowgain";
    name += i;
    title = "Pedestal LEDRef, low gain, module ";
    title += i; 
    fPedestalLEDRefLowGain.Add(new TProfile(name, title,
					    fLEDRefs, 0.0, fLEDRefs, "s"));
    
    //LED Ref/Mon pedestals, high gain
    name = "hPedestalLEDRefhighgain";
    name += i;
    title = "Pedestal LEDRef, high gain, module ";
    title += i; 
    fPedestalLEDRefHighGain.Add(new TProfile(name, title,
					     fLEDRefs, 0.0, fLEDRefs, "s"));
  
    //Peak-Pedestals, low gain
    name = "hPeakMinusPedlowgain";
    name += i;
    title = "Peak-Pedestal, low gain, module ";
    title += i; 
    fPeakMinusPedLowGain.Add(new TProfile2D(name, title,
					    fColumns, 0.0, fColumns, 
					    fRows, fRowMin, fRowMax,"s"));
  
    //Peak-Pedestals, high gain
    name = "hPeakMinusPedhighgain";
    name += i;
    title = "Peak-Pedestal, high gain, module ";
    title += i; 
    fPeakMinusPedHighGain.Add(new TProfile2D(name, title,
					     fColumns, 0.0, fColumns, 
					     fRows, fRowMin, fRowMax,"s"));

    //Peak-Pedestals, high gain - TH2F histo
    name = "hPeakMinusPedhighgainHisto";
    name += i;
    title = "Peak-Pedestal, high gain, module ";
    title += i; 
    fPeakMinusPedHighGainHisto.Add(new TH2F(name, title,
					    fColumns*fRows, 0.0, fColumns*fRows, 
					    100, 0, 1000));
 
    name = "hDeadMap";
    name += i;
    title = "Dead map, module ";
    title += i;
    fDeadMap.Add(new TH2D(name, title, fColumns, 0.0, fColumns, 
			  fRows, fRowMin, fRowMax));
  
  }//end for nModules create the histograms
 
  //Compress the arrays, in order to remove the empty objects (a 16 slot array is created by default)
  fPedestalLowGain.Compress();
  fPedestalHighGain.Compress();
  fPedestalLEDRefLowGain.Compress();
  fPedestalLEDRefHighGain.Compress();
  fPeakMinusPedLowGain.Compress();
  fPeakMinusPedHighGain.Compress();
  fPeakMinusPedHighGainHisto.Compress();
  fDeadMap.Compress();

}

// dtor
//_____________________________________________________________________
AliCaloCalibPedestal::~AliCaloCalibPedestal()
{
  if (fReference) delete fReference;//Delete the reference object, if it has been loaded
  //TObjArray will delete the histos/profiles when it is deleted.
}

// copy ctor
//_____________________________________________________________________
AliCaloCalibPedestal::AliCaloCalibPedestal(const AliCaloCalibPedestal &ped) :
  TObject(ped),
  fPedestalLowGain(),
  fPedestalHighGain(),
  fPedestalLEDRefLowGain(),
  fPedestalLEDRefHighGain(),
  fPeakMinusPedLowGain(),
  fPeakMinusPedHighGain(),
  fPeakMinusPedHighGainHisto(),
  fPedestalLowGainDiff(),
  fPedestalHighGainDiff(),
  fPedestalLEDRefLowGainDiff(),
  fPedestalLEDRefHighGainDiff(),
  fPeakMinusPedLowGainDiff(),
  fPeakMinusPedHighGainDiff(),
  fPedestalLowGainRatio(),
  fPedestalHighGainRatio(),
  fPedestalLEDRefLowGainRatio(),
  fPedestalLEDRefHighGainRatio(),
  fPeakMinusPedLowGainRatio(),
  fPeakMinusPedHighGainRatio(),
  fDeadMap(),
  fNEvents(ped.GetNEvents()),
  fNChanFills(ped.GetNChanFills()),
  fDeadTowers(ped.GetDeadTowerCount()),
  fNewDeadTowers(ped.GetDeadTowerNew()),
  fResurrectedTowers(ped.GetDeadTowerResurrected()),
  fReference( 0 ), //! note that we do not try to copy the reference info here
  fDetType(ped.GetDetectorType()),
  fColumns(ped.GetColumns()),
  fRows(ped.GetRows()),
  fLEDRefs(ped.GetLEDRefs()),
  fModules(ped.GetModules()),
  fRowMin(ped.GetRowMin()),
  fRowMax(ped.GetRowMax()),
  fRowMultiplier(ped.GetRowMultiplier()),
  fCaloString(ped.GetCaloString()),
  fMapping(NULL), //! note that we are not copying the map info
  fRunNumber(ped.GetRunNumber()),
  fSelectPedestalSamples(ped.GetSelectPedestalSamples()),
  fFirstPedestalSample(ped.GetFirstPedestalSample()),
  fLastPedestalSample(ped.GetLastPedestalSample()),
  fDeadThreshold(ped.GetDeadThreshold()),
  fWarningThreshold(ped.GetWarningThreshold()),
  fWarningFraction(ped.GetWarningFraction()),
  fHotSigma(ped.GetHotSigma())
{
  // Then the ObjArray ones; we add the histograms rather than trying TObjArray = assignment
  //DS: this has not really been tested yet..
  for (int i = 0; i < fModules; i++) {
    fPedestalLowGain.Add( ped.GetPedProfileLowGain(i) );
    fPedestalHighGain.Add( ped.GetPedProfileHighGain(i) );
    fPedestalLEDRefLowGain.Add( ped.GetPedLEDRefProfileLowGain(i) );
    fPedestalLEDRefHighGain.Add( ped.GetPedLEDRefProfileHighGain(i) );
    fPeakMinusPedLowGain.Add( ped.GetPeakProfileLowGain(i) );
    fPeakMinusPedHighGain.Add( ped.GetPeakProfileHighGain(i) );
    fPeakMinusPedHighGainHisto.Add( ped.GetPeakHighGainHisto(i) );

    fDeadMap.Add( ped.GetDeadMap(i) );  
  }//end for nModules 
 
  //Compress the arrays, in order to remove the empty objects (a 16 slot array is created by default)
  fPedestalLowGain.Compress();
  fPedestalHighGain.Compress();
  fPedestalLEDRefLowGain.Compress();
  fPedestalLEDRefHighGain.Compress();
  fPeakMinusPedLowGain.Compress();
  fPeakMinusPedHighGain.Compress();
  fPeakMinusPedHighGainHisto.Compress();

  fDeadMap.Compress();
}

// assignment operator; use copy ctor to make life easy..
//_____________________________________________________________________
AliCaloCalibPedestal& AliCaloCalibPedestal::operator = (const AliCaloCalibPedestal &source)
{
  // assignment operator; use copy ctor
  if (&source == this) return *this;

  new (this) AliCaloCalibPedestal(source);
  return *this;
}

//_____________________________________________________________________
void AliCaloCalibPedestal::Reset()
{
  // Reset all arrays/histograms
  for (int i = 0; i < fModules; i++) {
    GetPedProfileLowGain(i)->Reset();
    GetPedProfileHighGain(i)->Reset();
    GetPedLEDRefProfileLowGain(i)->Reset();
    GetPedLEDRefProfileHighGain(i)->Reset();
    GetPeakProfileLowGain(i)->Reset();
    GetPeakProfileHighGain(i)->Reset();
    GetPeakHighGainHisto(i)->Reset();
    GetDeadMap(i)->Reset();
    
    if (!fPedestalLowGainDiff.IsEmpty()) {
      //This means that the comparison profiles have been created.
  
      GetPedProfileLowGainDiff(i)->Reset();
      GetPedProfileHighGainDiff(i)->Reset();
      GetPedLEDRefProfileLowGainDiff(i)->Reset();
      GetPedLEDRefProfileHighGainDiff(i)->Reset();
      GetPeakProfileLowGainDiff(i)->Reset();
      GetPeakProfileHighGainDiff(i)->Reset();
      
      GetPedProfileLowGainRatio(i)->Reset();
      GetPedProfileHighGainRatio(i)->Reset();
      GetPedLEDRefProfileLowGainRatio(i)->Reset();
      GetPedLEDRefProfileHighGainRatio(i)->Reset();
      GetPeakProfileLowGainRatio(i)->Reset();
      GetPeakProfileHighGainRatio(i)->Reset();
    }
  }
  fNEvents = 0;
  fNChanFills = 0;
  fDeadTowers = 0;
  fNewDeadTowers = 0;
  fResurrectedTowers = 0;
 
  //To think about: should fReference be deleted too?... let's not do it this time, at least...
}

// Parameter/cut handling
//_____________________________________________________________________
void AliCaloCalibPedestal::SetParametersFromFile(const char *parameterFile)
{  
  // Note: this method is a bit more complicated than it really has to be
  // - allowing for multiple entries per line, arbitrary order of the
  // different variables etc. But I wanted to try and do this in as
  // correct a C++ way as I could (as an exercise).

  static const string delimitor("::");
	
  // open, check input file
  ifstream in( parameterFile );
  if( !in ) {
    printf("in AliCaloCalibPedestal::SetParametersFromFile - Using default/run_time parameters.\n");
    return;
  } 


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
      if( s.rdstate() & ios::failbit ) break;
			
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
      // printf("AliCaloCalibPedestal::SetParametersFromFile - key %s value %s\n", key.c_str(), value.c_str());

      // if the key matches with something we expect, we assign the new value
      istringstream iss(value);
      // the comparison strings defined at the beginning of this method
      if ( (key == "fFirstPedestalSample") || (key == "fLastPedestalSample") || (key == "fDeadThreshold") || (key == "fWarningThreshold") || (key == "fWarningFraction") || (key == "fHotSigma") ) {
	printf("AliCaloCalibPedestal::SetParametersFromFile - key %s value %s\n", key.c_str(), value.c_str());

	if (key == "fFirstPedestalSample") { 
	  iss >> fFirstPedestalSample; 
	}
	else if (key == "fLastPedestalSample") { 
	  iss >> fLastPedestalSample; 
	}
	else if (key == "fDeadThreshold") { 
	  iss >> fDeadThreshold; 
	}
	else if (key == "fWarningThreshold") { 
	  iss >> fWarningThreshold; 
	}
	else if (key == "fWarningFraction") { 
	  iss >> fWarningFraction; 
	}
	else if (key == "fHotSigma") { 
	  iss >> fHotSigma; 
	}

      } // some match

    }		
  }

  in.close();
  return;
	
}

//_____________________________________________________________________
void AliCaloCalibPedestal::WriteParametersToFile(const char *parameterFile)
{
  //Write parameters in file.
	
  static const string delimitor("::");
  ofstream out( parameterFile );
  out << "// " << parameterFile << endl;
  out << "fFirstPedestalSample" << "::" << fFirstPedestalSample << endl;
  out << "fLastPedestalSample" << "::" << fLastPedestalSample << endl;
  out << "fDeadThreshold" << "::" << fDeadThreshold << endl;
  out << "fWarningThreshold" << "::" << fWarningThreshold << endl;
  out << "fWarningFraction" << "::" << fWarningFraction << endl;
  out << "fHotSigma" << "::" << fHotSigma << endl;

  out.close();
  return;
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::AddInfo(const AliCaloCalibPedestal *ped)
{
  // just do this for the basic histograms/profiles that get filled in ProcessEvent
  // may not have data for all modules, but let's just Add everything..
  for (int i = 0; i < fModules; i++) {
    GetPedProfileLowGain(i)->Add( ped->GetPedProfileLowGain(i) );
    GetPedProfileHighGain(i)->Add( ped->GetPedProfileHighGain(i) );
    GetPeakProfileLowGain(i)->Add( ped->GetPeakProfileLowGain(i) );
    GetPeakProfileHighGain(i)->Add( ped->GetPeakProfileHighGain(i) );
    GetPeakHighGainHisto(i)->Add( ped->GetPeakHighGainHisto(i) );

  }//end for nModules 

  // DeadMap; Diff profiles etc would need to be redone after this operation

  return kTRUE;//We succesfully added info from the supplied object
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::ProcessEvent(AliRawReader *rawReader)
{ 
  // if fMapping is NULL the rawstream will crate its own mapping
  AliCaloRawStreamV3 rawStream(rawReader, fCaloString, (AliAltroMapping**)fMapping);
  if (fDetType == kEmCal) {
    rawReader->Select("EMCAL", 0, AliEMCALGeoParams::fgkLastAltroDDL) ; //select EMCAL DDL range 
  }
  return ProcessEvent(&rawStream);
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::ProcessEvent(AliCaloRawStreamV3 *in)
{ 
  // Method to process=analyze one event in the data stream
  if (!in) return kFALSE; //Return right away if there's a null pointer
  fNEvents++; // one more event
  
  // indices for the reading
  int sample = 0;
  int time = 0;
  int i = 0; // sample counter
  int startBin = 0;

  // start loop over input stream 
  while (in->NextDDL()) {
    while (in->NextChannel()) {

      // counters
      int max = AliEMCALGeoParams::fgkSampleMin, min = AliEMCALGeoParams::fgkSampleMax; // min and max sample values
      int nsamples = 0;

      // pedestal samples
      int nPed = 0;
      vector<int> pedSamples; 

      while (in->NextBunch()) {
	const UShort_t *sig = in->GetSignals();
	startBin = in->GetStartTimeBin();
	nsamples += in->GetBunchLength();
	for (i = 0; i < in->GetBunchLength(); i++) {
	  sample = sig[i];
	  time = startBin--;

	  // check if it's a min or max value
	  if (sample < min) min = sample;
	  if (sample > max) max = sample;
	  
	  // should we add it for the pedestal calculation?
	  if ( (fFirstPedestalSample<=time && time<=fLastPedestalSample) || // sample time in range
	       !fSelectPedestalSamples ) { // or we don't restrict the sample range.. - then we'll take all 
	    pedSamples.push_back( sig[i] );
	    nPed++;
	  }
	  
	} // loop over samples in bunch
      } // loop over bunches

      if (nsamples > 0) { // this check is needed for when we have zero-supp. on, but not sparse readout

      // it should be enough to check the SuperModule info for each DDL really, but let's keep it here for now
      int arrayPos = in->GetModule(); //The modules are numbered starting from 0
      if (arrayPos >= fModules) {
	//TODO: return an error message, if appopriate (perhaps if debug>0?)
	return kFALSE;
      }     
      //Debug
      if (arrayPos < 0 || arrayPos >= fModules) {
	printf("Oh no: arrayPos = %i.\n", arrayPos); 
      }
      
      fNChanFills++; // one more channel found, and profile to be filled
      //NOTE: coordinates are (column, row) for the profiles
      if ( in->IsLowGain() ) {
	//fill the low gain histograms
	((TProfile2D*)fPeakMinusPedLowGain[arrayPos])->Fill(in->GetColumn(), fRowMultiplier*in->GetRow(), max - min);
	if (nPed>0) { // only fill pedestal info in case it could be calculated
	  for ( i=0; i<nPed; i++) {
	    ((TProfile2D*)fPedestalLowGain[arrayPos])->Fill(in->GetColumn(), fRowMultiplier*in->GetRow(), pedSamples[i]); 
	  }
	}
      } 
      else if ( in->IsHighGain() ) {	
      	//fill the high gain ones
	((TProfile2D*)fPeakMinusPedHighGain[arrayPos])->Fill(in->GetColumn(), fRowMultiplier*in->GetRow(), max - min);
	if (nPed>0) { // only fill pedestal info in case it could be calculated
	  for ( i=0; i<nPed; i++) {
	    ((TProfile2D*)fPedestalHighGain[arrayPos])->Fill(in->GetColumn(), fRowMultiplier*in->GetRow(), pedSamples[i]); 
	  }	  
	}
	// for warning checks
	int idx = in->GetRow() + fRows * in->GetColumn();
	((TH2F*)fPeakMinusPedHighGainHisto[arrayPos])->Fill(idx, max - min);
      } 
      else if ( in->IsLEDMonData() ) {
	// for LED Mon data, the mapping class holds the gain info in the Row variable
	// and the Strip number in the Column..
	int gain = in->GetRow(); 
	int stripId = in->GetColumn();
	if (nPed>0 && stripId<fLEDRefs) {
	  if (gain == 0) {
	    for ( i=0; i<nPed; i++) {
	      ((TProfile*)fPedestalLEDRefLowGain[arrayPos])->Fill(stripId, pedSamples[i]);
	    }
	  }
	  else {
	    for ( i=0; i<nPed; i++) {
	      ((TProfile*)fPedestalLEDRefHighGain[arrayPos])->Fill(stripId, pedSamples[i]);
	    }
	  }
	}
      }

      } // nsamples>0 check, some data found for this channel; not only trailer/header
    }// end while over channel   
  }//end while over DDL's, of input stream 

  in->Reset(); // just in case the next customer forgets to check if the stream was reset..
 
  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::SaveHistograms(TString fileName, Bool_t saveEmptyHistos)
{
  //Saves all the histograms (or profiles, to be accurate) to the designated file
  
  TFile destFile(fileName, "recreate");
  
  if (destFile.IsZombie()) {
    return kFALSE;
  }
  
  destFile.cd();
  
  for (int i = 0; i < fModules; i++) {
    if( ((TProfile2D *)fPeakMinusPedLowGain[i])->GetEntries() || saveEmptyHistos) {
      fPeakMinusPedLowGain[i]->Write();
    }
    if( ((TProfile2D *)fPeakMinusPedHighGain[i])->GetEntries() || saveEmptyHistos) { 
      fPeakMinusPedHighGain[i]->Write();
    }
    if( ((TProfile2D *)fPedestalLowGain[i])->GetEntries() || saveEmptyHistos) {
      fPedestalLowGain[i]->Write();
    }
    if( ((TProfile2D *)fPedestalHighGain[i])->GetEntries() || saveEmptyHistos) {
      fPedestalHighGain[i]->Write();
    }
    if( ((TProfile *)fPedestalLEDRefLowGain[i])->GetEntries() || saveEmptyHistos) {
      fPedestalLEDRefLowGain[i]->Write();
    }
    if( ((TProfile *)fPedestalLEDRefHighGain[i])->GetEntries() || saveEmptyHistos) {
      fPedestalLEDRefHighGain[i]->Write();
    }
    if( ((TH2F *)fPeakMinusPedHighGainHisto[i])->GetEntries() || saveEmptyHistos) { 
      fPeakMinusPedHighGainHisto[i]->Write();
    }

  } 
  
  destFile.Close();
  
  return kTRUE;
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::LoadReferenceCalib(TString fileName, TString objectName)
{
  
  //Make sure that the histograms created when loading the object are not destroyed as the file object is destroyed
  TH1::AddDirectory(kFALSE);
  
  TFile *sourceFile = new TFile(fileName);
  if (sourceFile->IsZombie()) {
    return kFALSE;//We couldn't load the reference
  }

  if (fReference) delete fReference;//Delete the reference object, if it already exists
  fReference = 0;
  
  fReference = (AliCaloCalibPedestal*)sourceFile->Get(objectName);
 
  if (!fReference || !(fReference->InheritsFrom(AliCaloCalibPedestal::Class())) || (fReference->GetDetectorType() != fDetType)) {
    if (fReference) delete fReference;//Delete the object, in case we had an object of the wrong type
    fReference = 0;
    return kFALSE;
  }
	
  delete sourceFile;

  //Reset the histogram ownership behaviour. NOTE: a better workaround would be good, since this may accidentally set AddDirectory to true, even
  //if we are called by someone who has set it to false...
  TH1::AddDirectory(kTRUE);
 
  return kTRUE;//We succesfully loaded the object
}


//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::SetReference(AliCaloCalibPedestal *ref)
{
  if (fReference) delete fReference;//Delete the reference object, if it already exists
  fReference = 0;
  
  fReference = ref;
 
  if (!fReference || (fReference->GetDetectorType() != fDetType)) {
    if (fReference) delete fReference;//Delete the object, in case we had an object of the wrong type
    fReference = 0;
    return kFALSE;
  }

  return kTRUE;//We succesfully loaded the object
}

//_____________________________________________________________________
void AliCaloCalibPedestal::ValidateComparisonProfiles()
{
  //Make sure the comparison histos exist
  if (!fPedestalLowGainDiff.IsEmpty()) return; //The profiles already exist. We just check one, because they're all created at
  //the same time
						
						
  //Then, loop for the requested number of modules
  TString title, name;
  for (int i = 0; i < fModules; i++) {
    //Pedestals, low gain
    name = "hPedlowgainDiff";
    name += i;
    title = "Pedestals difference, low gain, module ";
    title += i; 
    fPedestalLowGainDiff.Add(new TProfile2D(name, title,
					    fColumns, 0.0, fColumns, 
					    fRows, fRowMin, fRowMax,"s"));
  
    //Pedestals, high gain
    name = "hPedhighgainDiff";
    name += i;
    title = "Pedestals difference, high gain, module ";
    title += i; 
    fPedestalHighGainDiff.Add(new TProfile2D(name, title,
					     fColumns, 0.0, fColumns, 
					     fRows, fRowMin, fRowMax,"s"));

    //LED Ref/Mon pedestals, low gain
    name = "hPedestalLEDReflowgainDiff";
    name += i;
    title = "Pedestal difference LEDRef, low gain, module ";
    title += i; 
    fPedestalLEDRefLowGainDiff.Add(new TProfile(name, title,
						fLEDRefs, 0.0, fLEDRefs, "s"));
    
    //LED Ref/Mon pedestals, high gain
    name = "hPedestalLEDRefhighgainDiff";
    name += i;
    title = "Pedestal difference LEDRef, high gain, module ";
    title += i; 
    fPedestalLEDRefHighGainDiff.Add(new TProfile(name, title,
						 fLEDRefs, 0.0, fLEDRefs, "s"));

    //Peak-Pedestals, high gain
    name = "hPeakMinusPedhighgainDiff";
    name += i;
    title = "Peak-Pedestal difference, high gain, module ";
    title += i; 
    fPeakMinusPedHighGainDiff.Add(new TProfile2D(name, title,
						 fColumns, 0.0, fColumns, 
						 fRows, fRowMin, fRowMax,"s"));

    //Peak-Pedestals, low gain
    name = "hPeakMinusPedlowgainDiff";
    name += i;
    title = "Peak-Pedestal difference, low gain, module ";
    title += i; 
    fPeakMinusPedLowGainDiff.Add(new TProfile2D(name, title,
						fColumns, 0.0, fColumns, 
						fRows, fRowMin, fRowMax,"s"));
  
    //Pedestals, low gain
    name = "hPedlowgainRatio";
    name += i;
    title = "Pedestals ratio, low gain, module ";
    title += i; 
    fPedestalLowGainRatio.Add(new TProfile2D(name, title,
					     fColumns, 0.0, fColumns, 
					     fRows, fRowMin, fRowMax,"s"));
  
    //Pedestals, high gain
    name = "hPedhighgainRatio";
    name += i;
    title = "Pedestals ratio, high gain, module ";
    title += i; 
    fPedestalHighGainRatio.Add(new TProfile2D(name, title,
					      fColumns, 0.0, fColumns, 
					      fRows, fRowMin, fRowMax,"s"));

    //LED Ref/Mon pedestals, low gain
    name = "hPedestalLEDReflowgainRatio";
    name += i;
    title = "Pedestal ratio LEDRef, low gain, module ";
    title += i; 
    fPedestalLEDRefLowGainRatio.Add(new TProfile(name, title,
						 fLEDRefs, 0.0, fLEDRefs, "s"));
    
    //LED Ref/Mon pedestals, high gain
    name = "hPedestalLEDRefhighgainRatio";
    name += i;
    title = "Pedestal ratio LEDRef, high gain, module ";
    title += i; 
    fPedestalLEDRefHighGainRatio.Add(new TProfile(name, title,
						  fLEDRefs, 0.0, fLEDRefs, "s"));
  
    //Peak-Pedestals, low gain
    name = "hPeakMinusPedlowgainRatio";
    name += i;
    title = "Peak-Pedestal ratio, low gain, module ";
    title += i; 
    fPeakMinusPedLowGainRatio.Add(new TProfile2D(name, title,
						 fColumns, 0.0, fColumns, 
						 fRows, fRowMin, fRowMax,"s"));
  
    //Peak-Pedestals, high gain
    name = "hPeakMinusPedhighgainRatio";
    name += i;
    title = "Peak-Pedestal ratio, high gain, module ";
    title += i; 
    fPeakMinusPedHighGainRatio.Add(new TProfile2D(name, title,
						  fColumns, 0.0, fColumns, 
						  fRows, fRowMin, fRowMax,"s"));
    
  }//end for nModules create the histograms
}

//_____________________________________________________________________
void AliCaloCalibPedestal::ComputeDiffAndRatio()
{
  // calculate differences and ratios relative to a reference
  ValidateComparisonProfiles();//Make sure the comparison histos exist
 
  if (!fReference) {
    return;//Return if the reference object isn't loaded
  }

  int bin = 0;
  double diff = 0;
  double ratio = 1;
  for (int i = 0; i < fModules; i++) {
    //For computing the difference, we cannot simply do TProfile2D->Add(), because that subtracts the sum of all entries,
    //which means that the mean of the new profile will not be the difference of the means. So do it by hand:
    for (int j = 0; j < fColumns; j++) {
      for (int k = 0; k < fRows; k++) {
	bin = ((TProfile2D*)fPeakMinusPedHighGainDiff[i])->GetBin(j+1, k+1);//Note that we assume here that all histos have the same structure...

	if (fReference->GetPeakProfileHighGain(i)->GetBinContent(bin) > 0) {
	  diff = GetPeakProfileHighGain(i)->GetBinContent(bin) - fReference->GetPeakProfileHighGain(i)->GetBinContent(bin);
	  ((TProfile2D*)fPeakMinusPedHighGainDiff[i])->SetBinContent(bin, diff);
	  ((TProfile2D*)fPeakMinusPedHighGainDiff[i])->SetBinEntries(bin, 1);
	  ratio = GetPeakProfileHighGain(i)->GetBinContent(bin) / fReference->GetPeakProfileHighGain(i)->GetBinContent(bin);  
	  ((TProfile2D*)fPeakMinusPedHighGainRatio[i])->SetBinContent(bin, ratio);
	  ((TProfile2D*)fPeakMinusPedHighGainRatio[i])->SetBinEntries(bin, 1);
	}

	if (fReference->GetPeakProfileLowGain(i)->GetBinContent(bin) > 0) {
	  diff = GetPeakProfileLowGain(i)->GetBinContent(bin) - fReference->GetPeakProfileLowGain(i)->GetBinContent(bin);
	  ((TProfile2D*)fPeakMinusPedLowGainDiff[i])->SetBinContent(bin, diff);
	  ((TProfile2D*)fPeakMinusPedLowGainDiff[i])->SetBinEntries(bin, 1);
	  ratio = GetPeakProfileLowGain(i)->GetBinContent(bin) / fReference->GetPeakProfileLowGain(i)->GetBinContent(bin);  
	  ((TProfile2D*)fPeakMinusPedLowGainRatio[i])->SetBinContent(bin, ratio);
	  ((TProfile2D*)fPeakMinusPedLowGainRatio[i])->SetBinEntries(bin, 1);
	}

	if (fReference->GetPedProfileHighGain(i)->GetBinContent(bin) > 0) {
	  diff = GetPedProfileHighGain(i)->GetBinContent(bin) - fReference->GetPedProfileHighGain(i)->GetBinContent(bin);
	  ((TProfile2D*)fPedestalHighGainDiff[i])->SetBinContent(bin, diff);
	  ((TProfile2D*)fPedestalHighGainDiff[i])->SetBinEntries(bin, 1);
	  ratio = GetPedProfileHighGain(i)->GetBinContent(bin) / fReference->GetPedProfileHighGain(i)->GetBinContent(bin);  
	  ((TProfile2D*)fPedestalHighGainRatio[i])->SetBinContent(bin, ratio);
	  ((TProfile2D*)fPedestalHighGainRatio[i])->SetBinEntries(bin, 1);
	}

	if (fReference->GetPedProfileLowGain(i)->GetBinContent(bin) > 0) {
	  diff = GetPedProfileLowGain(i)->GetBinContent(bin) - fReference->GetPedProfileLowGain(i)->GetBinContent(bin);
	  ((TProfile2D*)fPedestalLowGainDiff[i])->SetBinContent(bin, diff);
	  ((TProfile2D*)fPedestalLowGainDiff[i])->SetBinEntries(bin, 1);
	  ratio = GetPedProfileLowGain(i)->GetBinContent(bin) / fReference->GetPedProfileLowGain(i)->GetBinContent(bin);  
	  ((TProfile2D*)fPedestalLowGainRatio[i])->SetBinContent(bin, ratio);
	  ((TProfile2D*)fPedestalLowGainRatio[i])->SetBinEntries(bin, 1);
	}

      } // rows
    } // columns

    // same for LED Ref/Mon channels
    for (int j = 0; j <= fLEDRefs; j++) {    
      bin = j+1;//Note that we assume here that all histos have the same structure...

      if (fReference->GetPedLEDRefProfileHighGain(i)->GetBinContent(bin) > 0) {
	diff = GetPedLEDRefProfileHighGain(i)->GetBinContent(bin) - fReference->GetPedLEDRefProfileHighGain(i)->GetBinContent(bin);
	((TProfile*)fPedestalLEDRefHighGainDiff[i])->SetBinContent(bin, diff);
	((TProfile*)fPedestalLEDRefHighGainDiff[i])->SetBinEntries(bin, 1);
	ratio = GetPedLEDRefProfileHighGain(i)->GetBinContent(bin) / fReference->GetPedLEDRefProfileHighGain(i)->GetBinContent(bin);  
	((TProfile*)fPedestalLEDRefHighGainRatio[i])->SetBinContent(bin, ratio);
	((TProfile*)fPedestalLEDRefHighGainRatio[i])->SetBinEntries(bin, 1);
      }

      if (fReference->GetPedLEDRefProfileLowGain(i)->GetBinContent(bin) > 0) {
	diff = GetPedLEDRefProfileLowGain(i)->GetBinContent(bin) - fReference->GetPedLEDRefProfileLowGain(i)->GetBinContent(bin);
	((TProfile*)fPedestalLEDRefLowGainDiff[i])->SetBinContent(bin, diff);
	((TProfile*)fPedestalLEDRefLowGainDiff[i])->SetBinEntries(bin, 1);
	ratio = GetPedLEDRefProfileLowGain(i)->GetBinContent(bin) / fReference->GetPedLEDRefProfileLowGain(i)->GetBinContent(bin);  
	((TProfile*)fPedestalLEDRefLowGainRatio[i])->SetBinContent(bin, ratio);
	((TProfile*)fPedestalLEDRefLowGainRatio[i])->SetBinEntries(bin, 1);
      } 
     
    }

  } // modules
 
}

//_____________________________________________________________________
void AliCaloCalibPedestal::ComputeHotAndWarningTowers(const char * hotMapFile)
{ // look for hot/noisy towers
  ofstream * fout = 0;
  char name[512];//Quite a long temp buffer, just in case the filename includes a path

  if (hotMapFile) {
    snprintf(name, 512, "%s.txt", hotMapFile);
    fout = new ofstream(name);
    if (!fout->is_open()) {
      delete fout;
      fout = 0;//Set the pointer to empty if the file was not opened
    }
  }
 
  for(int i = 0; i < fModules; i++){
		
    //first we compute the peak-pedestal distribution for each supermodule...
    if( GetPeakHighGainHisto(i)->GetEntries() > 0 ) {
      double min = GetPeakProfileHighGain(i)->GetBinContent(GetPeakProfileHighGain(i)->GetMinimumBin());
      double max = GetPeakProfileHighGain(i)->GetBinContent(GetPeakProfileHighGain(i)->GetMaximumBin());
      TH1D *hPeakFit = new TH1D(Form("hFit_%d", i), Form("hFit_%d", i), (int)((max-min)*10), min-1, max+1);

      for (int j = 1; j <= fColumns; j++) {
	for (int k = 1; k <= fRows; k++) {	  
	  hPeakFit->Fill(GetPeakProfileHighGain(i)->GetBinContent(j, k));
	}
      }

      //...and then we try to fit it
      double mean  = hPeakFit->GetMean();
      double sigma = hPeakFit->GetRMS();
      try {
	hPeakFit->Fit("gaus", "OQ", "",  mean - 3*sigma, mean + 3*sigma);
	mean  = hPeakFit->GetFunction("gaus")->GetParameter(1);
	sigma = hPeakFit->GetFunction("gaus")->GetParameter(2);
      }
      catch (const std::exception & e) {
	printf("AliCaloCalibPedestal: TH1D PeakFit exception %s", e.what()); 
      }      
      //hPeakFit->Draw();

      delete hPeakFit;

      //Then we look for warm/hot towers
      TH2F * hPeak2D = GetPeakHighGainHisto(i);
      hPeak2D->GetYaxis()->SetRangeUser( fWarningThreshold, hPeak2D->GetYaxis()->GetBinUpEdge(hPeak2D->GetNbinsY()) );

      int idx = 0 ;
      int warnCounter = 0;
      for (int j = 1; j <= fColumns; j++) {
	for (int k = 1; k <= fRows; k++) {				
	  //we start looking for warm/warning towers...
	  // histogram x-axis index
	  idx = k-1 + fRows*(j-1); // this is what is used in the Fill call
	  hPeak2D->GetXaxis()->SetRangeUser(idx, idx);
	  warnCounter = (int) hPeak2D->Integral();
	  if(warnCounter > fNEvents * fWarningFraction) {
	    ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kWarning); 	    
	    /* printf("mod %d col %d row %d warnCounter %d - status %d\n", 
	       i, j-1, k-1, warnCounter, (int) (kWarning)); */  
	  }
	  //...then we look for hot ones (towers whose values are greater than mean + X*sigma)	
	  if(GetPeakProfileHighGain(i)->GetBinContent(j, k) > mean + fHotSigma*sigma ) {
	    ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kHot); 
	    /* printf("mod %d col %d row %d  binc %d - status %d\n", 
	       i, j-1, k-1, (int)(GetPeakProfileHighGain(i)->GetBinContent(j, k)), (int) (kHot)); */	  
	  }

	  //Write the status to the hot/warm map file, if the file is open.
	  // module - column - row - status (1=dead, 2= warm/warning , 3 = hot, see .h file enum)
	  if (fout && ((TH2D*)fDeadMap[i])->GetBinContent(j, k) > 1) {
	    
	    (*fout) << i << " " 
		    << (j - 1) << " " 
		    << (k - 1) << " " 
		    << ((TH2D*)fDeadMap[i])->GetBinContent(j, k) << endl;	    	}
	  
	}
      }

    } 
  }
  return;
}

//_____________________________________________________________________
void AliCaloCalibPedestal::ComputeDeadTowers(const char * deadMapFile)
{
  //Computes the number of dead towers etc etc into memory, after this you can call the GetDead... -functions
  int countTot = 0;
  int countNew = 0;
  int countRes = 0;
  ofstream * fout = 0;
  ofstream * diff = 0;
  char name[512];//Quite a long temp buffer, just in case the filename includes a path
  
  if (deadMapFile) {
    snprintf(name, 512, "%s.txt", deadMapFile);
    fout = new ofstream(name);
    snprintf(name, 512, "%sdiff.txt", deadMapFile);
    diff = new ofstream(name);
    if (!fout->is_open()) {
      delete fout;
      fout = 0;//Set the pointer to empty if the file was not opened
    }
    if (!diff->is_open()) {
      delete diff;
      fout = 0;//Set the pointer to empty if the file was not opened
    }
  }
 
  for (int i = 0; i < fModules; i++) {
    if (GetPeakProfileHighGain(i)->GetEntries() > 0) { //don't care about empty histos
      for (int j = 1; j <= fColumns; j++) {
	for (int k = 1; k <= fRows; k++) {

	  if (GetPeakProfileHighGain(i)->GetBinContent(j, k) < fDeadThreshold) {//It's dead
	    countTot++;//One more dead total
	    if (fout) {
	      (*fout) << i << " " 
		      << (j - 1) << " " 
		      << (k - 1) << " " 
		      << "1" << " " 
		      << "0" << endl;//Write the status to the deadmap file, if the file is open.
	    }
	    
	    if (fReference && fReference->GetPeakProfileHighGain(i)->GetBinContent(j, k) >= fDeadThreshold) {
	      ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kRecentlyDeceased); 
	      countNew++;//This tower wasn't dead before!
	      if (diff) {
		( *diff) << i << " " 
			 << (j - 1) << " " 
			 << (k - 1) << " " 
			 << "1" << " " 
			 << "0" << endl;//Write the status to the deadmap difference file, if the file is open.
	      }
	    } 
	    else {
	      ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kDead);//This has been dead before. Nothing new		
	    }
	  } 
	  else { //It's ALIVE!!
	    //Don't bother with writing the live ones.
	    if (fReference && fReference->GetPeakProfileHighGain(i)->GetBinContent(j, k) < fDeadThreshold) {
	      ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kResurrected);
	      countRes++; //This tower was dead before => it's a miracle! :P
	      if (diff) {
		(*diff) << i << " " 
			<< (j - 1) << " " 
			<< (k - 1) << " " 
			<< "1" << " " 
			<< "1" << endl;//Write the status to the deadmap difference file, if the file is open.
	      }
	    } 
	    else {
	      ((TH2D*)fDeadMap[i])->SetBinContent(j, k, kAlive);
	    }
	  }
	    
	}//end for k/rows
      }//end for j/columns
    }//end if GetEntries >= 0
  
  }//end for modules
 
 if (fout) {
   fout->close();
   delete fout;
 }
 
 fDeadTowers = countTot;
 fNewDeadTowers = countNew;
 fResurrectedTowers = countRes;
}

//_____________________________________________________________________
Bool_t AliCaloCalibPedestal::IsBadChannel(int imod, int icol, int irow) const
{
  //Check if channel is dead or hot.  
  Int_t status =  (Int_t) ( ((TH2D*)fDeadMap[imod])->GetBinContent(icol,irow) );
  if(status == kAlive)
    return kFALSE;
  else 
    return kTRUE;
  
}

//_____________________________________________________________________
void AliCaloCalibPedestal::SetChannelStatus(int imod, int icol, int irow, int status)
{
  //Set status of channel dead, hot, alive ...  
  ((TH2D*)fDeadMap[imod])->SetBinContent(icol, irow, status);	
}
