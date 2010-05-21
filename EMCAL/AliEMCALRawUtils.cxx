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

//_________________________________________________________________________
//  Utility Class for handling Raw data
//  Does all transitions from Digits to Raw and vice versa, 
//  for simu and reconstruction
//
//  Note: the current version is still simplified. Only 
//    one raw signal per digit is generated; either high-gain or low-gain
//    Need to add concurrent high and low-gain info in the future
//    No pedestal is added to the raw signal.
//*-- Author: Marco van Leeuwen (LBL)

#include "AliEMCALRawUtils.h"
#include <stdexcept>
  
#include "TF1.h"
#include "TGraph.h"
#include <TRandom.h>
class TSystem;
  
class AliLog;
#include "AliRun.h"
#include "AliRunLoader.h"
class AliCaloAltroMapping;
#include "AliAltroBuffer.h"
#include "AliRawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliDAQ.h"
  
#include "AliEMCALRecParam.h"
#include "AliEMCALLoader.h"
#include "AliEMCALGeometry.h"
class AliEMCALDigitizer;
#include "AliEMCALDigit.h"
#include "AliEMCALRawDigit.h"
#include "AliEMCAL.h"
#include "AliCaloCalibPedestal.h"  
#include "AliCaloFastAltroFitv0.h"
#include "AliCaloNeuralFit.h"
#include "AliCaloBunchInfo.h"
#include "AliCaloFitResults.h"
#include "AliCaloRawAnalyzerFastFit.h"
#include "AliCaloRawAnalyzerNN.h"
#include "AliCaloRawAnalyzerLMS.h"
#include "AliCaloRawAnalyzerPeakFinder.h"
#include "AliCaloRawAnalyzerCrude.h"

ClassImp(AliEMCALRawUtils)
  
// Signal shape parameters
Int_t    AliEMCALRawUtils::fgTimeBins = 256; // number of sampling bins of the raw RO signal (we typically use 15-50; theoretical max is 1k+) 
Double_t AliEMCALRawUtils::fgTimeBinWidth  = 100E-9 ; // each sample is 100 ns
Double_t AliEMCALRawUtils::fgTimeTrigger = 1.5E-6 ;   // 15 time bins ~ 1.5 musec

// some digitization constants
Int_t    AliEMCALRawUtils::fgThreshold = 1;
Int_t    AliEMCALRawUtils::fgDDLPerSuperModule = 2;  // 2 ddls per SuperModule
Int_t    AliEMCALRawUtils::fgPedestalValue = 0;     // pedestal value for digits2raw, default generate ZS data
Double_t AliEMCALRawUtils::fgFEENoise = 3.;          // 3 ADC channels of noise (sampled)

AliEMCALRawUtils::AliEMCALRawUtils(fitAlgorithm fitAlgo)
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(0), fOption(""),
    fRemoveBadChannels(kTRUE),fFittingAlgorithm(0),  
    fTimeMin(-1.),fTimeMax(1.),
    fUseFALTRO(kFALSE),fRawAnalyzer(0)
{

  //These are default parameters.  
  //Can be re-set from without with setter functions
  //Already set in the OCDB and passed via setter in the AliEMCALReconstructor
  fHighLowGainFactor = 16. ;   // Adjusted for a low gain range of 82 GeV (10 bits) 
  fOrder             = 2;      // Order of gamma fn
  fTau               = 2.35;   // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold    = 3;      // 3 ADC counts is approx. noise level
  fNPedSamples       = 4;      // Less than this value => likely pedestal samples
  fRemoveBadChannels = kFALSE; // Do not remove bad channels before fitting
  fUseFALTRO         = kTRUE;  // Get the trigger FALTRO information and pass it to digits.
  SetFittingAlgorithm(fitAlgo);

  //Get Mapping RCU files from the AliEMCALRecParam                                 
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  //To make sure we match with the geometry in a simulation file,
  //let's try to get it first.  If not, take the default geometry
  AliRunLoader *rl = AliRunLoader::Instance();
  if (rl && rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  } else {
    AliDebug(1, Form("Using default geometry in raw reco"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

}

//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(AliEMCALGeometry *pGeometry, fitAlgorithm fitAlgo)
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(pGeometry), fOption(""),
    fRemoveBadChannels(kTRUE),fFittingAlgorithm(0),
    fTimeMin(-1.),fTimeMax(1.),
    fUseFALTRO(kFALSE),fRawAnalyzer()
{
  //
  // Initialize with the given geometry - constructor required by HLT
  // HLT does not use/support AliRunLoader(s) instances
  // This is a minimum intervention solution
  // Comment by MPloskon@lbl.gov
  //

  //These are default parameters. 
  //Can be re-set from without with setter functions 
  //Already set in the OCDB and passed via setter in the AliEMCALReconstructor
  fHighLowGainFactor = 16. ;   // adjusted for a low gain range of 82 GeV (10 bits)
  fOrder             = 2;      // order of gamma fn
  fTau               = 2.35;   // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold    = 3;      // 3 ADC counts is approx. noise level
  fNPedSamples       = 4;      // Less than this value => likely pedestal samples
  fRemoveBadChannels = kFALSE; // Do not remove bad channels before fitting
  fUseFALTRO         = kTRUE;  // Get the trigger FALTRO information and pass it to digits.
  SetFittingAlgorithm(fitAlgo);

  //Get Mapping RCU files from the AliEMCALRecParam
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

}

//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(const AliEMCALRawUtils& rawU)
  : TObject(),
    fHighLowGainFactor(rawU.fHighLowGainFactor), 
    fOrder(rawU.fOrder),
    fTau(rawU.fTau),
    fNoiseThreshold(rawU.fNoiseThreshold),
    fNPedSamples(rawU.fNPedSamples),
    fGeom(rawU.fGeom), 
    fOption(rawU.fOption),
    fRemoveBadChannels(rawU.fRemoveBadChannels),
    fFittingAlgorithm(rawU.fFittingAlgorithm),
    fTimeMin(rawU.fTimeMin),fTimeMax(rawU.fTimeMax),
	fUseFALTRO(rawU.fUseFALTRO),
    fRawAnalyzer(rawU.fRawAnalyzer)
{
  //copy ctor
  fMapping[0] = rawU.fMapping[0];
  fMapping[1] = rawU.fMapping[1];
  fMapping[2] = rawU.fMapping[2];
  fMapping[3] = rawU.fMapping[3];
}

//____________________________________________________________________________
AliEMCALRawUtils& AliEMCALRawUtils::operator =(const AliEMCALRawUtils &rawU)
{
  //assignment operator

  if(this != &rawU) {
    fHighLowGainFactor = rawU.fHighLowGainFactor;
    fOrder             = rawU.fOrder;
    fTau               = rawU.fTau;
    fNoiseThreshold    = rawU.fNoiseThreshold;
    fNPedSamples       = rawU.fNPedSamples;
    fGeom              = rawU.fGeom;
    fOption            = rawU.fOption;
    fRemoveBadChannels = rawU.fRemoveBadChannels;
    fFittingAlgorithm  = rawU.fFittingAlgorithm;
	fTimeMin           = rawU.fTimeMin;
	fTimeMax           = rawU.fTimeMax;
    fUseFALTRO         = rawU.fUseFALTRO;
    fRawAnalyzer       = rawU.fRawAnalyzer;
    fMapping[0]        = rawU.fMapping[0];
    fMapping[1]        = rawU.fMapping[1];
    fMapping[2]        = rawU.fMapping[2];
    fMapping[3]        = rawU.fMapping[3];
  }

  return *this;

}

//____________________________________________________________________________
AliEMCALRawUtils::~AliEMCALRawUtils() {
  //dtor

}

//____________________________________________________________________________
void AliEMCALRawUtils::Digits2Raw()
{
  // convert digits of the current event to raw data
  
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *loader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  // get the digits
  loader->LoadDigits("EMCAL");
  loader->GetEvent();
  TClonesArray* digits = loader->Digits() ;
  
  if (!digits) {
    Warning("Digits2Raw", "no digits found !");
    return;
  }

  static const Int_t nDDL = 12*2; // 12 SM hardcoded for now. Buffers allocated dynamically, when needed, so just need an upper limit here
  AliAltroBuffer* buffers[nDDL];
  for (Int_t i=0; i < nDDL; i++)
    buffers[i] = 0;

  TArrayI adcValuesLow(fgTimeBins);
  TArrayI adcValuesHigh(fgTimeBins);

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliEMCALDigit* digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit)) ;
    if (digit->GetAmplitude() < fgThreshold) 
      continue;

    //get cell indices
    Int_t nSM = 0;
    Int_t nIphi = 0;
    Int_t nIeta = 0;
    Int_t iphi = 0;
    Int_t ieta = 0;
    Int_t nModule = 0;
    fGeom->GetCellIndex(digit->GetId(), nSM, nModule, nIphi, nIeta);
    fGeom->GetCellPhiEtaIndexInSModule(nSM, nModule, nIphi, nIeta,iphi, ieta) ;
    
    //Check which is the RCU, 0 or 1, of the cell.
    Int_t iRCU = -111;
    //RCU0
    if (0<=iphi&&iphi<8) iRCU=0; // first cable row
    else if (8<=iphi&&iphi<16 && 0<=ieta&&ieta<24) iRCU=0; // first half; 
    //second cable row
    //RCU1
    else if(8<=iphi&&iphi<16 && 24<=ieta&&ieta<48) iRCU=1; // second half; 
    //second cable row
    else if(16<=iphi&&iphi<24) iRCU=1; // third cable row

    if (nSM%2==1) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same

    if (iRCU<0) 
      Fatal("Digits2Raw()","Non-existent RCU number: %d", iRCU);
    
    //Which DDL?
    Int_t iDDL = fgDDLPerSuperModule* nSM + iRCU;
    if (iDDL >= nDDL)
      Fatal("Digits2Raw()","Non-existent DDL board number: %d", iDDL);

    if (buffers[iDDL] == 0) {      
      // open new file and write dummy header
      TString fileName = AliDAQ::DdlFileName("EMCAL",iDDL);
      //Select mapping file RCU0A, RCU0C, RCU1A, RCU1C
      Int_t iRCUside=iRCU+(nSM%2)*2;
      //iRCU=0 and even (0) SM -> RCU0A.data   0
      //iRCU=1 and even (0) SM -> RCU1A.data   1
      //iRCU=0 and odd  (1) SM -> RCU0C.data   2
      //iRCU=1 and odd  (1) SM -> RCU1C.data   3
      //cout<<" nSM "<<nSM<<"; iRCU "<<iRCU<<"; iRCUside "<<iRCUside<<endl;
      buffers[iDDL] = new AliAltroBuffer(fileName.Data(),fMapping[iRCUside]);
      buffers[iDDL]->WriteDataHeader(kTRUE, kFALSE);  //Dummy;
    }
    
    // out of time range signal (?)
    if (digit->GetTimeR() > GetRawFormatTimeMax() ) {
      AliInfo("Signal is out of time range.\n");
      buffers[iDDL]->FillBuffer((Int_t)digit->GetAmplitude());
      buffers[iDDL]->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffers[iDDL]->FillBuffer(3);          // bunch length      
      buffers[iDDL]->WriteTrailer(3, ieta, iphi, nSM);  // trailer
      // calculate the time response function
    } else {
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), digit->GetAmplitude(), adcValuesHigh.GetArray(), adcValuesLow.GetArray()) ; 
      if (lowgain) 
	buffers[iDDL]->WriteChannel(ieta, iphi, 0, GetRawFormatTimeBins(), adcValuesLow.GetArray(), fgThreshold);
      else 
	buffers[iDDL]->WriteChannel(ieta,iphi, 1, GetRawFormatTimeBins(), adcValuesHigh.GetArray(), fgThreshold);
    }
  }
  
  // write headers and close files
  for (Int_t i=0; i < nDDL; i++) {
    if (buffers[i]) {
      buffers[i]->Flush();
      buffers[i]->WriteDataHeader(kFALSE, kFALSE);
      delete buffers[i];
    }
  }

  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr, const AliCaloCalibPedestal* pedbadmap, TClonesArray *digitsTRG)
{
  // convert raw data of the current event to digits                                                                                     

  digitsArr->Clear(); 

  if (!digitsArr) {
    Error("Raw2Digits", "no digits found !");
    return;
  }
  if (!reader) {
    Error("Raw2Digits", "no raw reader found !");
    return;
  }

  AliCaloRawStreamV3 in(reader,"EMCAL",fMapping);
  // Select EMCAL DDL's;
  reader->Select("EMCAL",0,43); // 43 = AliEMCALGeoParams::fgkLastAltroDDL

  // fRawAnalyzer setup
  fRawAnalyzer->SetNsampleCut(5); // requirement for fits to be done
  fRawAnalyzer->SetAmpCut(fNoiseThreshold);
  fRawAnalyzer->SetFitArrayCut(fNoiseThreshold);
  fRawAnalyzer->SetIsZeroSuppressed(true); // TMP - should use stream->IsZeroSuppressed(), or altro cfg registers later

  // channel info parameters
  Int_t lowGain  = 0;
  Int_t caloFlag = 0; // low, high gain, or TRU, or LED ref.

  // start loop over input stream 
  while (in.NextDDL()) {
	  
//    if ( in.GetDDLNumber() != 0 && in.GetDDLNumber() != 2 ) continue;

    while (in.NextChannel()) {

/*
	  Int_t    hhwAdd    = in.GetHWAddress();
	  UShort_t iiBranch  = ( hhwAdd >> 11 ) & 0x1; // 0/1
	  UShort_t iiFEC     = ( hhwAdd >>  7 ) & 0xF;
	  UShort_t iiChip    = ( hhwAdd >>  4 ) & 0x7;
	  UShort_t iiChannel =   hhwAdd         & 0xF;
		 
	  if ( !( iiBranch == 0 && iiFEC == 1 && iiChip == 3 && ( iiChannel >= 8 && iiChannel <= 15 ) ) && !( iiBranch == 1 && iiFEC == 0 && in.GetColumn() == 0 ) ) continue;
*/
		
      //Check if the signal  is high or low gain and then do the fit, 
      //if it  is from TRU or LEDMon do not fit
      caloFlag = in.GetCaloFlag();
//		if (caloFlag != 0 && caloFlag != 1) continue; 
	  if (caloFlag > 2) continue; // Work with ALTRO and FALTRO 
		
      //Do not fit bad channels of ALTRO
      if(caloFlag < 2 && fRemoveBadChannels && pedbadmap->IsBadChannel(in.GetModule(),in.GetColumn(),in.GetRow())) {
	//printf("Tower from SM %d, column %d, row %d is BAD!!! Skip \n", in.GetModule(),in.GetColumn(),in.GetRow());
	continue;
      }  

      vector<AliCaloBunchInfo> bunchlist; 
      while (in.NextBunch()) {
	bunchlist.push_back( AliCaloBunchInfo(in.GetStartTimeBin(), in.GetBunchLength(), in.GetSignals() ) );
      } // loop over bunches

   
      if ( caloFlag < 2 ){ // ALTRO
		
	Float_t time = 0; 
	Float_t amp  = 0; 
	short timeEstimate  = 0;
	Float_t ampEstimate = 0;
	Bool_t fitDone = kFALSE;
		
      if ( fFittingAlgorithm == kFastFit || fFittingAlgorithm == kNeuralNet || fFittingAlgorithm == kLMS || fFittingAlgorithm == kPeakFinder || fFittingAlgorithm == kCrude) {
	// all functionality to determine amp and time etc is encapsulated inside the Evaluate call for these methods 
	AliCaloFitResults fitResults = fRawAnalyzer->Evaluate( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2()); 

	amp          = fitResults.GetAmp();
	time         = fitResults.GetTime();
	timeEstimate = fitResults.GetMaxTimebin();
	ampEstimate  = fitResults.GetMaxSig();
	if (fitResults.GetStatus() == AliCaloFitResults::kFitPar) {
	  fitDone = kTRUE;
	} 
      }
      else { // for the other methods we for now use the functionality of 
	// AliCaloRawAnalyzer as well, to select samples and prepare for fits, 
	// if it looks like there is something to fit

	// parameters init.
	Float_t pedEstimate  = 0;
	short maxADC = 0;
	Int_t first = 0;
	Int_t last = 0;
	Int_t bunchIndex = 0;
	//
	// The PreFitEvaluateSamples + later call to FitRaw will hopefully 
	// be replaced by a single Evaluate call or so soon, like for the other
	// methods, but this should be good enough for evaluation of 
	// the methods for now (Jan. 2010)
	//
	int nsamples = fRawAnalyzer->PreFitEvaluateSamples( bunchlist, in.GetAltroCFG1(), in.GetAltroCFG2(), bunchIndex, ampEstimate, maxADC, timeEstimate, pedEstimate, first, last); 
	
	if (ampEstimate >= fNoiseThreshold) { // something worth looking at
	  
	  time = timeEstimate; // maxrev in AliCaloRawAnalyzer speak; comes with an offset w.r.t. real timebin
	  Int_t timebinOffset = bunchlist.at(bunchIndex).GetStartBin() - (bunchlist.at(bunchIndex).GetLength()-1); 
	  amp = ampEstimate; 
	  
	  if ( nsamples > 1 ) { // possibly something to fit
	    FitRaw(first, last, amp, time, fitDone);
	    time += timebinOffset;
	    timeEstimate += timebinOffset;
	  }
	  
	} // ampEstimate check
      } // method selection

      if ( fitDone ) { // brief sanity check of fit results	    
	Float_t ampAsymm = (amp - ampEstimate)/(amp + ampEstimate);
	Float_t timeDiff = time - timeEstimate;
	if ( (TMath::Abs(ampAsymm) > 0.1) || (TMath::Abs(timeDiff) > 2) ) {
	  // AliDebug(2,Form("Fit results amp %f time %f not consistent with expectations amp %f time %d", amp, time, ampEstimate, timeEstimate));
	  
	  // for now just overwrite the fit results with the simple/initial estimate
	  amp     = ampEstimate;
	  time    = timeEstimate; 
	  fitDone = kFALSE;
	} 
      } // fitDone
    
      if (amp >= fNoiseThreshold  && amp<fgkRawSignalOverflow) { // something to be stored
	if ( ! fitDone) { // smear ADC with +- 0.5 uniform (avoid discrete effects)
	  amp += (0.5 - gRandom->Rndm()); // Rndm generates a number in ]0,1]
	}

	Int_t id = fGeom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
	lowGain  = in.IsLowGain();

	// go from time-bin units to physical time fgtimetrigger
	time = time * GetRawFormatTimeBinWidth(); // skip subtraction of fgTimeTrigger?
	// subtract RCU L1 phase (L1Phase is in seconds) w.r.t. L0:
	time -= in.GetL1Phase();

	AliDebug(2,Form("id %d lowGain %d amp %g", id, lowGain, amp));
	// printf("Added tower: SM %d, row %d, column %d, amp %3.2f\n",in.GetModule(), in.GetRow(), in.GetColumn(),amp);
	AddDigit(digitsArr, id, lowGain, amp, time); 
      }
      
	}//ALTRO
	else if(fUseFALTRO)
	{// Fake ALTRO
		//		if (maxTimeBin && gSig->GetN() > maxTimeBin + 10) gSig->Set(maxTimeBin + 10); // set actual max size of TGraph
		Int_t    hwAdd    = in.GetHWAddress();
		UShort_t iRCU     = in.GetDDLNumber() % 2; // 0/1
		UShort_t iBranch  = ( hwAdd >> 11 ) & 0x1; // 0/1
		
		// Now find TRU number
		Int_t itru = 3 * in.GetModule() + ( (iRCU << 1) | iBranch ) - 1;
		
		AliDebug(1,Form("Found TRG digit in TRU: %2d ADC: %2d",itru,in.GetColumn()));
		
		Int_t idtrg;
		
		Bool_t isOK = fGeom->GetAbsFastORIndexFromTRU(itru, in.GetColumn(), idtrg);
		
		Int_t timeSamples[256]; for (Int_t j=0;j<256;j++) timeSamples[j] = 0;
		Int_t nSamples = 0;
		
		for (std::vector<AliCaloBunchInfo>::iterator itVectorData = bunchlist.begin(); itVectorData != bunchlist.end(); itVectorData++)
		{
			AliCaloBunchInfo bunch = *(itVectorData);
			
			const UShort_t* sig = bunch.GetData();
			Int_t startBin = bunch.GetStartBin();
			
			for (Int_t iS = 0; iS < bunch.GetLength(); iS++) 
			{
		  		Int_t time = startBin--;
				Int_t amp  = sig[iS];
				
				if ( amp ) timeSamples[nSamples++] = ( ( time << 12 ) & 0xFF000 ) | ( amp & 0xFFF );
			}
		}
		
		if (nSamples && isOK) AddDigit(digitsTRG, idtrg, timeSamples, nSamples);
	}//Fake ALTRO
   } // end while over channel   
  } //end while over DDL's, of input stream 

  TrimDigits(digitsArr);
	
  return ; 
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t timeSamples[], Int_t nSamples) 
{
  //Add raw sample to raw digit 
  new((*digitsArr)[digitsArr->GetEntriesFast()]) AliEMCALRawDigit(id, timeSamples, nSamples);	
  
  //	Int_t idx = digitsArr->GetEntriesFast()-1;
  //	AliEMCALRawDigit* d = (AliEMCALRawDigit*)digitsArr->At(idx);
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Float_t amp, Float_t time) {
  //
  // Add a new digit. 
  // This routine checks whether a digit exists already for this tower 
  // and then decides whether to use the high or low gain info
  //
  // Called by Raw2Digits
  
  AliEMCALDigit *digit = 0, *tmpdigit = 0;
  TIter nextdigit(digitsArr);
  while (digit == 0 && (tmpdigit = (AliEMCALDigit*) nextdigit())) {
    if (tmpdigit->GetId() == id) digit = tmpdigit;
  }

  if (!digit) { // no digit existed for this tower; create one
		Int_t type = AliEMCALDigit::kHG; // use enum in AliEMCALDigit
		if (lowGain) { 
			amp *= fHighLowGainFactor;
			type = AliEMCALDigit::kLGnoHG;
		} 
		Int_t idigit = digitsArr->GetEntries();
		new((*digitsArr)[idigit]) AliEMCALDigit( -1, -1, id, amp, time, type, idigit) ; 
		AliDebug(2,Form("Add digit Id %d for the first time, type %d", id, type));
  }//digit added first time
  else { // a digit already exists, check range 
		// (use high gain if signal < cut value, otherwise low gain)
		if (lowGain) { // new digit is low gain
			if (digit->GetAmplitude() > fgkOverflowCut) {  // use if previously stored (HG) digit is out of range
				digit->SetAmplitude(fHighLowGainFactor * amp);
				digit->SetTime(time);
				digit->SetType(AliEMCALDigit::kLG);
				AliDebug(2,Form("Add LG digit ID %d for the second time, type %d", digit->GetId(), digit->GetType()));
			}
		}//new low gain digit
		else { // new digit is high gain 
			if (amp < fgkOverflowCut) { // new digit is high gain; use if not out of range
				digit->SetAmplitude(amp);
				digit->SetTime(time);
				digit->SetType(AliEMCALDigit::kHG);
				AliDebug(2,Form("Add HG digit ID %d for the second time, type %d", digit->GetId(), digit->GetType()));
			}
			else { // HG out of range, just change flag value to show that HG did exist
				digit->SetType(AliEMCALDigit::kLG);
				AliDebug(2,Form("Change LG digit to HG, ID %d, type %d", digit->GetId(), digit->GetType()));
			}
		}//new high gain digit
  }//digit existed replace it
  
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::TrimDigits(TClonesArray *digitsArr) 
{
  // Remove digits with only low gain and large time
  
  AliEMCALDigit *digit = 0;
  Int_t n = 0;
  Int_t nDigits = digitsArr->GetEntriesFast();
  TIter nextdigit(digitsArr);
  while ((digit = (AliEMCALDigit*) nextdigit())) {
    
    //Check if only LG existed, remove if so
    if (digit->GetType() == AliEMCALDigit::kLGnoHG) {
      AliDebug(1,Form("Remove digit with id %d, LGnoHG",digit->GetId()));
      digitsArr->Remove(digit);
    }
    //Check if time if too large or too small, remove if so
    else if(fTimeMin > digit->GetTime() || fTimeMax < digit->GetTime()) {
      digitsArr->Remove(digit);
      AliDebug(1,Form("Remove digit with id %d, Bad Time %e",digit->GetId(), digit->GetTime()));
    }
    //Good digit, just reassign the index of the digit in case there was a previous removal
    else {
      digit->SetIndexInList(n);	
      n++;
    }    
  }//while
  
  digitsArr->Compress();
  AliDebug(1,Form("N Digits before trimming : %d; after array compression %d",nDigits,digitsArr->GetEntriesFast()));
	   
}
	
//____________________________________________________________________________ 
void AliEMCALRawUtils::FitRaw(const Int_t firstTimeBin, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Bool_t & fitDone) const 
{ // Fits the raw signal time distribution
  
  //--------------------------------------------------
  //Do the fit, different fitting algorithms available
  //--------------------------------------------------
  int nsamples = lastTimeBin - firstTimeBin + 1;
  fitDone = kFALSE;

  switch(fFittingAlgorithm) {
  case kStandard:
    {
      if (nsamples < 3) { return; } // nothing much to fit
      //printf("Standard fitter \n");

      // Create Graph to hold data we will fit 
      TGraph *gSig =  new TGraph( nsamples); 
      for (int i=0; i<nsamples; i++) {
	Int_t timebin = firstTimeBin + i;    
	gSig->SetPoint(i, timebin, fRawAnalyzer->GetReversed(timebin)); 
      }

      TF1 * signalF = new TF1("signal", RawResponseFunction, 0, GetRawFormatTimeBins(), 5);
      signalF->SetParameters(10.,5.,fTau,fOrder,0.); //set all defaults once, just to be safe
      signalF->SetParNames("amp","t0","tau","N","ped");
      signalF->FixParameter(2,fTau); // tau in units of time bin
      signalF->FixParameter(3,fOrder); // order
      signalF->FixParameter(4, 0); // pedestal should be subtracted when we get here 
      signalF->SetParameter(1, time);
      signalF->SetParameter(0, amp);
      // set rather loose parameter limits
      signalF->SetParLimits(0, 0.5*amp, 2*amp );
      signalF->SetParLimits(1, time - 4, time + 4); 

      try {			
	gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points
	// assign fit results
	amp  = signalF->GetParameter(0); 
	time = signalF->GetParameter(1);

	// cross-check with ParabolaFit to see if the results make sense
	FitParabola(gSig, amp); // amp is possibly updated
	fitDone = kTRUE;
      }
      catch (const std::exception & e) {
	AliError( Form("TGraph Fit exception %s", e.what()) ); 
	// stay with default amp and time in case of exception, i.e. no special action required
	fitDone = kFALSE;
      }
      delete signalF;

      //printf("Std   : Amp %f, time %g\n",amp, time);
      delete gSig; // delete TGraph
				
      break;
    }//kStandard Fitter
    //----------------------------
  case kLogFit:
    {
      if (nsamples < 3) { return; } // nothing much to fit
      //printf("LogFit \n");

      // Create Graph to hold data we will fit 
      TGraph *gSigLog =  new TGraph( nsamples); 
      for (int i=0; i<nsamples; i++) {
	Int_t timebin = firstTimeBin + i;    
	gSigLog->SetPoint(timebin, timebin, TMath::Log(fRawAnalyzer->GetReversed(timebin) ) ); 
      }

      TF1 * signalFLog = new TF1("signalLog", RawResponseFunctionLog, 0, GetRawFormatTimeBins(), 5);
      signalFLog->SetParameters(2.3, 5.,fTau,fOrder,0.); //set all defaults once, just to be safe
      signalFLog->SetParNames("amplog","t0","tau","N","ped");
      signalFLog->FixParameter(2,fTau); // tau in units of time bin
      signalFLog->FixParameter(3,fOrder); // order
      signalFLog->FixParameter(4, 0); // pedestal should be subtracted when we get here 
      signalFLog->SetParameter(1, time);
      if (amp>=1) {
	signalFLog->SetParameter(0, TMath::Log(amp));
      }
	
      gSigLog->Fit(signalFLog, "QROW"); // Note option 'W': equal errors on all points
				
      // assign fit results
      Double_t amplog = signalFLog->GetParameter(0); //Not Amp, but Log of Amp
      amp = TMath::Exp(amplog);
      time = signalFLog->GetParameter(1);
      fitDone = kTRUE;

      delete signalFLog;
      //printf("LogFit: Amp %f, time %g\n",amp, time);
      delete gSigLog; 
      break;
    } //kLogFit 
    //----------------------------	
    
    //----------------------------
  }//switch fitting algorithms

  return;
}

//__________________________________________________________________
void AliEMCALRawUtils::FitParabola(const TGraph *gSig, Float_t & amp) const 
{
  //BEG YS alternative methods to calculate the amplitude
  Double_t * ymx = gSig->GetX() ; 
  Double_t * ymy = gSig->GetY() ; 
  const Int_t kN = 3 ; 
  Double_t ymMaxX[kN] = {0., 0., 0.} ; 
  Double_t ymMaxY[kN] = {0., 0., 0.} ; 
  Double_t ymax = 0. ; 
  // find the maximum amplitude
  Int_t ymiMax = 0 ;  
  for (Int_t ymi = 0; ymi < gSig->GetN(); ymi++) {
    if (ymy[ymi] > ymMaxY[0] ) {
      ymMaxY[0] = ymy[ymi] ; //<========== This is the maximum amplitude
      ymMaxX[0] = ymx[ymi] ;
      ymiMax = ymi ; 
    }
  }
  // find the maximum by fitting a parabola through the max and the two adjacent samples
  if ( ymiMax < gSig->GetN()-1 && ymiMax > 0) {
    ymMaxY[1] = ymy[ymiMax+1] ;
    ymMaxY[2] = ymy[ymiMax-1] ; 
    ymMaxX[1] = ymx[ymiMax+1] ;
    ymMaxX[2] = ymx[ymiMax-1] ; 
    if (ymMaxY[0]*ymMaxY[1]*ymMaxY[2] > 0) {
      //fit a parabola through the 3 points y= a+bx+x*x*x
      Double_t sy = 0 ; 
      Double_t sx = 0 ; 
      Double_t sx2 = 0 ; 
      Double_t sx3 = 0 ; 
      Double_t sx4 = 0 ; 
      Double_t sxy = 0 ; 
      Double_t sx2y = 0 ; 
      for (Int_t i = 0; i < kN ; i++) {
	sy += ymMaxY[i] ; 
	sx += ymMaxX[i] ; 		
	sx2 += ymMaxX[i]*ymMaxX[i] ; 
	sx3 += ymMaxX[i]*ymMaxX[i]*ymMaxX[i] ; 
	sx4 += ymMaxX[i]*ymMaxX[i]*ymMaxX[i]*ymMaxX[i] ; 
	sxy += ymMaxX[i]*ymMaxY[i] ; 
	sx2y += ymMaxX[i]*ymMaxX[i]*ymMaxY[i] ; 
      }
      Double_t cN = (sx2y*kN-sy*sx2)*(sx3*sx-sx2*sx2)-(sx2y*sx-sxy*sx2)*(sx3*kN-sx*sx2); 
      Double_t cD = (sx4*kN-sx2*sx2)*(sx3*sx-sx2*sx2)-(sx4*sx-sx3*sx2)*(sx3*kN-sx*sx2) ;
      Double_t c  = cN / cD ; 
      Double_t b  = ((sx2y*kN-sy*sx2)-c*(sx4*kN-sx2*sx2))/(sx3*kN-sx*sx2) ;
      Double_t a  = (sy-b*sx-c*sx2)/kN  ;
      Double_t xmax = -b/(2*c) ; 
      ymax = a + b*xmax + c*xmax*xmax ;//<========== This is the maximum amplitude
      amp = ymax;
    }
  }
  
  Double_t diff = TMath::Abs(1-ymMaxY[0]/amp) ; 
  if (diff > 0.1) 
    amp = ymMaxY[0] ; 
  //printf("Yves   : Amp %f, time %g\n",amp, time);
  //END YS
  return;
}

//__________________________________________________________________
Double_t AliEMCALRawUtils::RawResponseFunction(Double_t *x, Double_t *par)
{
  // Matches version used in 2007 beam test
  //
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  //
  // xx = (t - t0 + tau) / tau  [xx is just a convenient help variable]
  // F = A * (xx**N * exp( N * ( 1 - xx) )   for xx >= 0
  // F = 0                                   for xx < 0 
  //
  // parameters:
  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]
  // tau: par[2]
  // N:   par[3]
  // ped: par[4]
  //
  Double_t signal ;
  Double_t tau =par[2];
  Double_t n =par[3];
  Double_t ped = par[4];
  Double_t xx = ( x[0] - par[1] + tau ) / tau ;

  if (xx <= 0) 
    signal = ped ;  
  else {  
    signal = ped + par[0] * TMath::Power(xx , n) * TMath::Exp(n * (1 - xx )) ; 
  }
  return signal ;  
}

//__________________________________________________________________
Double_t AliEMCALRawUtils::RawResponseFunctionLog(Double_t *x, Double_t *par)
{
  // Matches version used in 2007 beam test
  //
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  //
  // xx = (t - t0 + tau) / tau  [xx is just a convenient help variable]
  // F = A * (xx**N * exp( N * ( 1 - xx) )   for xx >= 0
  // F = 0                                   for xx < 0 
  //
  // parameters:
  // Log[A]:   par[0]   // Amplitude = peak value
  // t0:  par[1]
  // tau: par[2]
  // N:   par[3]
  // ped: par[4]
  //
  Double_t signal ;
  Double_t tau =par[2];
  Double_t n =par[3];
  //Double_t ped = par[4]; // not used
  Double_t xx = ( x[0] - par[1] + tau ) / tau ;

  if (xx < 0) 
    signal = par[0] - n*TMath::Log(TMath::Abs(xx)) + n * (1 - xx ) ;  
  else {  
    signal = par[0] + n*TMath::Log(xx) + n * (1 - xx ) ; 
  }
  return signal ;  
}

//__________________________________________________________________
Bool_t AliEMCALRawUtils::RawSampledResponse(const Double_t dtime, const Double_t damp, Int_t * adcH, Int_t * adcL, const Int_t keyErr) const 
{
  // for a start time dtime and an amplitude damp given by digit, 
  // calculates the raw sampled response AliEMCAL::RawResponseFunction

  Bool_t lowGain = kFALSE ; 

  // A:   par[0]   // Amplitude = peak value
  // t0:  par[1]                            
  // tau: par[2]                            
  // N:   par[3]                            
  // ped: par[4]

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeBins(), 5);
  signalF.SetParameter(0, damp) ; 
  signalF.SetParameter(1, (dtime + fgTimeTrigger)/fgTimeBinWidth) ; 
  signalF.SetParameter(2, fTau) ; 
  signalF.SetParameter(3, fOrder);
  signalF.SetParameter(4, fgPedestalValue);
	
  Double_t signal=0.0, noise=0.0;
  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    signal = signalF.Eval(iTime) ;  
	
    // Next lines commeted for the moment but in principle it is not necessary to add
    // extra noise since noise already added at the digits level.	

    //According to Terry Awes, 13-Apr-2008
    //add gaussian noise in quadrature to each sample
    //Double_t noise = gRandom->Gaus(0.,fgFEENoise);
    //signal = sqrt(signal*signal + noise*noise);

    // March 17,09 for fast fit simulations by Alexei Pavlinov.
    // Get from PHOS analysis. In some sense it is open questions.
	if(keyErr>0) {
		noise = gRandom->Gaus(0.,fgFEENoise);
		signal += noise; 
	}
	  
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcH[iTime] > fgkRawSignalOverflow ){  // larger than 10 bits 
      adcH[iTime] = fgkRawSignalOverflow ;
      lowGain = kTRUE ; 
    }

    signal /= fHighLowGainFactor;

    adcL[iTime] =  static_cast<Int_t>(signal + 0.5) ;
    if ( adcL[iTime] > fgkRawSignalOverflow)  // larger than 10 bits 
      adcL[iTime] = fgkRawSignalOverflow ;
  }
  return lowGain ; 
}

//__________________________________________________________________
void AliEMCALRawUtils::CalculateChi2(const Double_t* t, const Double_t* y, const Int_t nPoints, 
const Double_t sig, const Double_t tau, const Double_t amp, const Double_t t0, Double_t &chi2)
{
  //   Input:
  //   t[]   - array of time bins
  //   y[]   - array of amplitudes after pedestal subtractions;
  //   nPoints  - number of points 
  //   sig   - error of amplitude measurement (one value for all channels)
  //           if sig<0 that mean sig=1.
  //   tau   - filter time response (in timebin units)
  //   amp   - amplitude at t0;
  //   t0    - time of max amplitude; 
  // Output:
  //   chi2 - chi2
  //   ndf = nPoints - 2 when tau fixed 
  //   ndf = nPoints - 3 when tau free
  static Double_t par[5]={0.0, 0.0, 0.0, 2.0, 0.0};

  par[0] = amp;
  par[1] = t0;
  par[2] = tau;
  // par[3]=n=2.; par[4]=ped=0.0

  Double_t dy = 0.0, x = 0.0, f=0.0;
  for(Int_t i=0; i<nPoints; i++){
    x     = t[i];
    f     = RawResponseFunction(&x, par);
    dy    = y[i] - f;
    chi2 += dy*dy;
    printf(" AliEMCALRawUtils::CalculateChi2 : %i : y %f -> f %f : dy %f \n", i, y[i], f, dy); 
  }
  if(sig>0.0) chi2 /= (sig*sig);
}

//__________________________________________________________________
void AliEMCALRawUtils::SetFittingAlgorithm(Int_t fitAlgo)              
{
	//Set fitting algorithm and initialize it if this same algorithm was not set before.
	//printf("**** Set Algorithm , number %d ****\n",fitAlgo);

	if(fitAlgo == fFittingAlgorithm && fRawAnalyzer) {
		//Do nothing, this same algorithm already set before.
		//printf("**** Algorithm already set before, number %d, %s ****\n",fitAlgo, fRawAnalyzer->GetName());
		return;
	}
	//Initialize the requested algorithm
	if(fitAlgo != fFittingAlgorithm || !fRawAnalyzer) {
		//printf("**** Init Algorithm , number %d ****\n",fitAlgo);
		
		fFittingAlgorithm = fitAlgo; 
		if (fRawAnalyzer) delete fRawAnalyzer;  // delete prev. analyzer if existed.
		
		if (fitAlgo == kFastFit) {
			fRawAnalyzer = new AliCaloRawAnalyzerFastFit();
		}
		else if (fitAlgo == kNeuralNet) {
			fRawAnalyzer = new AliCaloRawAnalyzerNN();
		}
		else if (fitAlgo == kLMS) {
			fRawAnalyzer = new AliCaloRawAnalyzerLMS();
		}
		else if (fitAlgo == kPeakFinder) {
			fRawAnalyzer = new AliCaloRawAnalyzerPeakFinder();
		}
		else if (fitAlgo == kCrude) {
			fRawAnalyzer = new AliCaloRawAnalyzerCrude();
		}
		else {
			fRawAnalyzer = new AliCaloRawAnalyzer();
		}
	}
	
}


