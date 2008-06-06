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
// asked to draw graphs or profiles:
//  fSignals->GetGraphAmpVsTimeHighGain(module,column,row)->Draw("ap");
// or
//  fSignals->GetProfAmpVsTimeHighGain(module,column,row)->Draw();
// etc.
//________________________________________________________________________

#include "TFile.h"

#include "AliRawReader.h"
#include "AliRawEventHeaderBase.h"
#include "AliCaloRawStream.h"

//The include file
#include "AliCaloCalibSignal.h"

ClassImp(AliCaloCalibSignal)

using namespace std;

// ctor; initialize everything in order to avoid compiler warnings
// put some reasonable defaults
AliCaloCalibSignal::AliCaloCalibSignal(kDetType detectorType) :  
  TObject(),
  fDetType(kNone),
  fColumns(0),
  fRows(0),
  fModules(0),
  fCaloString(),
  fMapping(NULL),
  fRunNumber(-1),
  fStartTime(0),
  fAmpCut(50),
  fReqFractionAboveAmpCutVal(0.8),
  fReqFractionAboveAmp(kTRUE),
  fHour(0),
  fLatestHour(0),
  fUseAverage(kTRUE),
  fSecInAverage(1800), 
  fNEvents(0),
  fNAcceptedEvents(0)
{
  //Default constructor. First we set the detector-type related constants.
  if (detectorType == kPhos) {
    fColumns = fgkPhosCols;
    fRows = fgkPhosRows;
    fModules = fgkPhosModules;
    fCaloString = "PHOS";
  } 
  else {
    //We'll just trust the enum to keep everything in line, so that if detectorType
    //isn't kPhos then it is kEmCal. Note, however, that this is not necessarily the
    //case, if someone intentionally gives another number
    fColumns = fgkEmCalCols;
    fRows = fgkEmCalRows;
    fModules = fgkEmCalModules;
    fCaloString = "EMCAL";
  }

  fDetType = detectorType;

  // Set the number of points for each Amp vs. Time graph to 0
  memset(fNHighGain, 0, sizeof(fNHighGain));
  memset(fNLowGain, 0, sizeof(fNLowGain));

  CreateGraphs(); // set up the TGraphs

  // init TProfiles to NULL=0 also
  memset(fProfAmpVsTimeHighGain, 0, sizeof(fProfAmpVsTimeHighGain));
  memset(fProfAmpVsTimeLowGain, 0, sizeof(fProfAmpVsTimeLowGain));
} 

// dtor
//_____________________________________________________________________
AliCaloCalibSignal::~AliCaloCalibSignal()
{
  ClearObjects();
}

//_____________________________________________________________________
void AliCaloCalibSignal::ClearObjects()
{
  // delete what was created in the ctor (TGraphs), and possible later (TProfiles)
  for (int i=0; i<fgkMaxTowers; i++) {
    if ( fGraphAmpVsTimeHighGain[i] ) { delete fGraphAmpVsTimeHighGain[i]; }
    if ( fGraphAmpVsTimeLowGain[i] ) { delete fGraphAmpVsTimeLowGain[i]; }
    if ( fProfAmpVsTimeHighGain[i] ) { delete fProfAmpVsTimeHighGain[i]; }
    if ( fProfAmpVsTimeLowGain[i] ) { delete fProfAmpVsTimeLowGain[i]; }
  }
  // set pointers
  memset(fGraphAmpVsTimeHighGain, 0, sizeof(fGraphAmpVsTimeHighGain));
  memset(fGraphAmpVsTimeLowGain, 0, sizeof(fGraphAmpVsTimeLowGain));
  memset(fProfAmpVsTimeHighGain, 0, sizeof(fProfAmpVsTimeHighGain));
  memset(fProfAmpVsTimeLowGain, 0, sizeof(fProfAmpVsTimeLowGain));

  return;
}

// copy ctor
//_____________________________________________________________________
AliCaloCalibSignal::AliCaloCalibSignal(const AliCaloCalibSignal &sig) :
  TObject(sig),
  fDetType(sig.GetDetectorType()),
  fColumns(sig.GetColumns()),
  fRows(sig.GetRows()),
  fModules(sig.GetModules()),
  fCaloString(sig.GetCaloString()),
  fMapping(NULL), //! note that we are not copying the map info
  fRunNumber(sig.GetRunNumber()),
  fStartTime(sig.GetStartTime()),
  fAmpCut(sig.GetAmpCut()),
  fReqFractionAboveAmpCutVal(sig.GetReqFractionAboveAmpCutVal()),
  fReqFractionAboveAmp(sig.GetReqFractionAboveAmp()),
  fHour(sig.GetHour()),
  fLatestHour(sig.GetLatestHour()),
  fUseAverage(sig.GetUseAverage()),
  fSecInAverage(sig.GetSecInAverage()),
  fNEvents(sig.GetNEvents()),
  fNAcceptedEvents(sig.GetNAcceptedEvents())
{
  // also the TGraph contents
  AddInfo(&sig);
}

// assignment operator; use copy ctor to make life easy..
//_____________________________________________________________________
AliCaloCalibSignal& AliCaloCalibSignal::operator = (const AliCaloCalibSignal &source)
{
  // assignment operator; use copy ctor
  if (&source == this) return *this;

  new (this) AliCaloCalibSignal(source);
  return *this;
}

//_____________________________________________________________________
void AliCaloCalibSignal::CreateGraphs()
{
  //Then, loop for the requested number of modules
  TString title, name;
  for (int i = 0; i < fModules; i++) {
    
    // Amplitude vs. Time graph for each channel
    for(int ic=0;ic < fColumns;ic++){
      for(int ir=0;ir < fRows;ir++){
	
	int id = GetTowerNum(i, ic, ir);
	
	// high gain graph
	name = "fGraphAmpVsTimeHighGain_"; name += i;
	name += "_"; name += ic;
	name += "_"; name += ir;
	title = "Amp vs. Time High Gain Mod "; title += i;
	title += " Col "; title += ic;
	title += " Row "; title += ir;
	
	fGraphAmpVsTimeHighGain[id] = new TGraph();
	fGraphAmpVsTimeHighGain[id]->SetName(name);
	fGraphAmpVsTimeHighGain[id]->SetTitle(title);
	fGraphAmpVsTimeHighGain[id]->SetMarkerStyle(20);
	
	// Low Gain
	name = "fGraphAmpVsTimeLowGain_"; name += i;
	name += "_"; name += ic;
	name += "_"; name += ir;
	title = "Amp vs. Time Low Gain Mod "; title += i;
	title += " Col "; title += ic;
	title += " Row "; title += ir;
	
	fGraphAmpVsTimeLowGain[id] = new TGraph();
	fGraphAmpVsTimeLowGain[id]->SetName(name);
	fGraphAmpVsTimeLowGain[id]->SetTitle(title);
	fGraphAmpVsTimeLowGain[id]->SetMarkerStyle(20);
	
      }
    }

  }//end for nModules 
}

//_____________________________________________________________________
void AliCaloCalibSignal::Reset()
{
  Zero(); // set all counters to 0
  ClearObjects(); // delete previous TGraphs and TProfiles
  CreateGraphs(); // and create some new ones
  return;
}

//_____________________________________________________________________
void AliCaloCalibSignal::Zero()
{
  // set all counters to 0; not cuts etc.though
  fHour = 0;
  fLatestHour = 0;
  fNEvents = 0;
  fNAcceptedEvents = 0;
  return;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::CheckFractionAboveAmp(int *AmpVal, int nTotChan)
{
  int nAbove = 0;
    
  int TowerNum = 0;
  for (int i = 0; i<fModules; i++) {
    for (int j = 0; j<fColumns; j++) {
      for (int k = 0; k<fRows; k++) {
	TowerNum = GetTowerNum(i,j,k);
	if (AmpVal[TowerNum] > fAmpCut) { 
	  nAbove++;
	}
      }
    }
  }
  
  double fraction = (1.0*nAbove) / nTotChan;
  
  if (fraction > fReqFractionAboveAmpCutVal) {  
    return true;
  }
  else return false;
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::AddInfo(const AliCaloCalibSignal *sig)
{
  // just do this for the basic graphs/profiles that get filled in ProcessEvent
  // may not have data for all channels, but let's just Add everything..
  // Note: this method will run into problems with TProfile adding if the binning of
  // the local profiles is not the same as those provided by the argument *sig..
  int numGraphPoints = 0;
  int id = 0;
  int ip = 0;
  for (int i = 0; i < fModules; i++) {
    for (int j = 0; j < fColumns; j++) {
      for (int k = 0; k < fRows; k++) {
	
	id = GetTowerNum(i,j,k);

	if(fUseAverage){ // add to Profiles
 	  if (sig->GetProfAmpVsTimeHighGain(id)) {
	    GetProfAmpVsTimeHighGain(id)->Add(sig->GetProfAmpVsTimeHighGain(id));
	  }
 	  if (sig->GetProfAmpVsTimeLowGain(id)) {
	    GetProfAmpVsTimeLowGain(id)->Add(sig->GetProfAmpVsTimeLowGain(id));
	  }
	}
	else{ // add to Graphs	  
	  // high gain
	  numGraphPoints= sig->GetGraphAmpVsTimeHighGain(id)->GetN();
	  if (numGraphPoints > 0) {
	    // get the values
	    double *graphX = sig->GetGraphAmpVsTimeHighGain(id)->GetX();
	    double *graphY = sig->GetGraphAmpVsTimeHighGain(id)->GetY();
	    for(ip=0; ip < numGraphPoints; ip++){
	      fGraphAmpVsTimeHighGain[id]->SetPoint(fNHighGain[id]++,graphX[ip],graphY[ip]);
	    }
	  }
	  // low gain
	  numGraphPoints= sig->GetGraphAmpVsTimeLowGain(id)->GetN();
	  if (numGraphPoints > 0) {
	    // get the values
	    double *graphX = sig->GetGraphAmpVsTimeLowGain(id)->GetX();
	    double *graphY = sig->GetGraphAmpVsTimeLowGain(id)->GetY();
	    for(ip=0; ip < numGraphPoints; ip++){
	      fGraphAmpVsTimeLowGain[id]->SetPoint(fNLowGain[id]++,graphX[ip],graphY[ip]);
	    }
	  }

	}

      }//end for nModules 
    }//end for nColumns
  }//end for nRows

  return kTRUE;//We succesfully added info from the supplied object
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::ProcessEvent(AliRawReader *rawReader)
{
  // if fMapping is NULL the rawstream will crate its own mapping
  AliCaloRawStream rawStream(rawReader, fCaloString, (AliAltroMapping**)fMapping);

  return ProcessEvent( &rawStream, (AliRawEventHeaderBase*)rawReader->GetEventHeader() );
}

//_____________________________________________________________________
Bool_t AliCaloCalibSignal::ProcessEvent(AliCaloRawStream *in, AliRawEventHeaderBase *aliHeader)
{ 
  // Method to process=analyze one event in the data stream
  if (!in) return kFALSE; //Return right away if there's a null pointer
  
  fNEvents++; // one more event

  // PHOS has more towers than EMCAL, so use PHOS numbers to set array sizes
  int AmpValHighGain[fgkMaxTowers];
  int AmpValLowGain[fgkMaxTowers];

  memset(AmpValHighGain, 0, sizeof(AmpValHighGain));
  memset(AmpValLowGain, 0, sizeof(AmpValLowGain));

  int sample, i = 0; //The sample temp, and the sample number in current event.
  int max = fgkSampleMin, min = fgkSampleMax;//Use these for picking the signal
  int gain = 0;
  
  // Number of Low and High gain channels for this event:
  int nLowChan = 0; 
  int nHighChan = 0; 

  int TowerNum = 0; // array index for TGraphs etc.

  // loop first to get the fraction of channels with amplitudes above cut
  while (in->Next()) {
    sample = in->GetSignal(); //Get the adc signal
    if (sample < min) min = sample;
    if (sample > max) max = sample;
    i++;
    if ( i >= in->GetTimeLength()) {
      //If we're here then we're done with this tower
      gain = 1 - in->IsLowGain();
      
      int arrayPos = in->GetModule(); //The modules are numbered starting from 0
      if (arrayPos >= fModules) {
	//TODO: return an error message, if appopriate (perhaps if debug>0?)
	return kFALSE;
      } 
      
      //Debug
      if (arrayPos < 0 || arrayPos >= fModules) {
	printf("Oh no: arrayPos = %i.\n", arrayPos); 
      }

      // get tower number for AmpVal array
      TowerNum = GetTowerNum(arrayPos, in->GetColumn(), in->GetRow()); 

      if (gain == 0) {
	// fill amplitude into the array	   
        AmpValLowGain[TowerNum] = max - min;
	nLowChan++;
      } 
      else if (gain==1) {//fill the high gain ones
	// fill amplitude into the array
	AmpValHighGain[TowerNum] = max - min;
	nHighChan++;
      }//end if gain

      
      max = fgkSampleMin; min = fgkSampleMax;
      i = 0;
      
    }//End if end of tower
   
  }//end while, of stream
  
  // now check if it was a led event, only use high gain (that should be sufficient)
  if (fReqFractionAboveAmp) {
    bool ok = false;
    if (nHighChan > 0) { 
      ok = CheckFractionAboveAmp(AmpValHighGain, nHighChan); 
    }
    if (!ok) return false; // skip event
  }

  fNAcceptedEvents++; // one more event accepted

  if (fStartTime == 0) { // if start-timestamp wasn't set,we'll pick it up from the first event we encounter
    fStartTime = aliHeader->Get("Timestamp");
  }

  fHour = (aliHeader->Get("Timestamp")-fStartTime)/(double)fgkNumSecInHr;
  if (fLatestHour < fHour) {
    fLatestHour = fHour; 
  }
  
  // it is a led event, now fill graphs (maybe profiles later)
  for(int i=0;i<fModules;i++){
    for(int j=0;j<fColumns;j++){
      for(int k=0;k<fRows;k++){
	
	TowerNum = GetTowerNum(i, j, k); 

	if(AmpValHighGain[TowerNum]) {
	  fGraphAmpVsTimeHighGain[TowerNum]->SetPoint(fNHighGain[TowerNum]++,fHour,AmpValHighGain[TowerNum]);
	}
	if(AmpValLowGain[TowerNum]) {
 	  fGraphAmpVsTimeLowGain[TowerNum]->SetPoint(fNLowGain[TowerNum]++,fHour,AmpValLowGain[TowerNum]);
	}
      }
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________
void AliCaloCalibSignal::CreateProfile(int imod, int ic, int ir, int towerId, int gain,
				       int nbins, double min, double max)
{ //! create/setup a TProfile
  TString title, name;   
  if (gain == 0) { 
    name = "fProfAmpVsTimeLowGain_";   
    title = "Amp vs. Time Low Gain Mod "; 
  } 
  else if (gain == 1) { 
    name = "fProfAmpVsTimeHighGain_"; 
    title = "Amp vs. Time High Gain Mod "; 
  } 
  name += imod;
  name += "_"; name += ic;
  name += "_"; name += ir;
  title += imod;
  title += " Col "; title += ic;
  title += " Row "; title += ir;
	    
  // use "s" option for RMS
  if (gain==0) { 
    fProfAmpVsTimeLowGain[towerId] = new TProfile(name,title, nbins, min, max,"s");
  }
  else if (gain==1) {
    fProfAmpVsTimeHighGain[towerId] = new TProfile(name,title, nbins, min, max,"s");
  }

  return;
}
//_____________________________________________________________________
Bool_t AliCaloCalibSignal::Save(TString fileName, Bool_t saveEmptyGraphs)
{
  //Saves all the histograms (or profiles, to be accurate) to the designated file
  
  TFile destFile(fileName, "recreate");
  
  if (destFile.IsZombie()) {
    return kFALSE;
  }
  
  destFile.cd();

  // setup variables for the TProfile plot
  int numProfBins = 0;
  double timeMin = 0;
  double timeMax = 0;
  if (fUseAverage) {
    if (fSecInAverage > 0) {
      numProfBins = (int)( (fLatestHour*fgkNumSecInHr)/fSecInAverage + 1 ); // round-off
    }
    numProfBins += 2; // add extra buffer : first and last
    double binSize = 1.0*fSecInAverage / fgkNumSecInHr;
    timeMin = - binSize;
    timeMax = timeMin + numProfBins*binSize;
  }

  int numGraphPoints= 0;
  int TowerNum = 0;    
  for (int i = 0; i < fModules; i++) {
    
    for(int ic=0;ic < fColumns;ic++){
      for(int ir=0;ir < fRows;ir++){

	TowerNum = GetTowerNum(i, ic, ir); 

	// 1st: high gain
	numGraphPoints= fGraphAmpVsTimeHighGain[TowerNum]->GetN();
	if( numGraphPoints>0 || saveEmptyGraphs) {
	  
	  // average the graphs points over time if requested and put them in a profile plot
	  if(fUseAverage && numGraphPoints>0) {
	    
	    // get the values
	    double *graphX = fGraphAmpVsTimeHighGain[TowerNum]->GetX();
	    double *graphY = fGraphAmpVsTimeHighGain[TowerNum]->GetY();

	    // create the TProfile: 1 is for High gain	    	    
	    CreateProfile(i, ic, ir, TowerNum, 1,
			  numProfBins, timeMin, timeMax);

	    // loop over graph points and fill profile
	    for(int ip=0; ip < numGraphPoints; ip++){
	      fProfAmpVsTimeHighGain[TowerNum]->Fill(graphX[ip],graphY[ip]);
	    }
	    
	    fProfAmpVsTimeHighGain[TowerNum]->GetXaxis()->SetTitle("Hours");
	    fProfAmpVsTimeHighGain[TowerNum]->GetYaxis()->SetTitle("MaxAmplitude - Pedestal");
	    fProfAmpVsTimeHighGain[TowerNum]->Write();

	  }
	   else{
	     //otherwise, just save the graphs and forget the profiling
	     fGraphAmpVsTimeHighGain[TowerNum]->GetXaxis()->SetTitle("Hours");
	     fGraphAmpVsTimeHighGain[TowerNum]->GetYaxis()->SetTitle("MaxAmplitude - Pedestal");
	     fGraphAmpVsTimeHighGain[TowerNum]->Write();
	   }
	  
	} // low gain graph info should be saved in one form or another
	
	// 2nd: now go to the low gain case
	numGraphPoints= fGraphAmpVsTimeLowGain[TowerNum]->GetN();
	if( numGraphPoints>0 || saveEmptyGraphs) {
	  
	  // average the graphs points over time if requested and put them in a profile plot
	  if(fUseAverage && numGraphPoints>0) {
	    
	    double *graphX = fGraphAmpVsTimeLowGain[TowerNum]->GetX();
	    double *graphY = fGraphAmpVsTimeLowGain[TowerNum]->GetY();
	    
	    // create the TProfile: 0 is for Low gain	    
	    CreateProfile(i, ic, ir, TowerNum, 0,
			  numProfBins, timeMin, timeMax);

	    // loop over graph points and fill profile
	    for(int ip=0; ip < numGraphPoints; ip++){
	      fProfAmpVsTimeLowGain[TowerNum]->Fill(graphX[ip],graphY[ip]);
	    }
	    
	    fProfAmpVsTimeLowGain[TowerNum]->GetXaxis()->SetTitle("Hours");
	    fProfAmpVsTimeLowGain[TowerNum]->GetYaxis()->SetTitle("MaxAmplitude - Pedestal");
	    fProfAmpVsTimeLowGain[TowerNum]->Write();

	  }
	  else{
	     //otherwise, just save the graphs and forget the profiling
	    fGraphAmpVsTimeLowGain[TowerNum]->GetXaxis()->SetTitle("Hours");
	    fGraphAmpVsTimeLowGain[TowerNum]->GetYaxis()->SetTitle("MaxAmplitude - Pedestal");
	    fGraphAmpVsTimeLowGain[TowerNum]->Write();
	  }
	  
	} // low gain graph info should be saved in one form or another

      } // fRows
    } // fColumns

  } // fModules
  destFile.Close();
  
  return kTRUE;
}
