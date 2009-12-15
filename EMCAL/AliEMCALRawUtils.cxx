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
  
#include "TF1.h"
#include "TGraph.h"
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
#include "AliEMCAL.h"
#include "AliCaloCalibPedestal.h"  
  
ClassImp(AliEMCALRawUtils)
  
// Signal shape parameters
Int_t    AliEMCALRawUtils::fgTimeBins = 256; // number of sampling bins of the raw RO signal (we typically use 15-50; theoretical max is 1k+) 
Double_t AliEMCALRawUtils::fgTimeBinWidth  = 100E-9 ; // each sample is 100 ns
Double_t AliEMCALRawUtils::fgTimeTrigger = 1.5E-6 ;   // 15 time bins ~ 1.5 musec

// some digitization constants
Int_t    AliEMCALRawUtils::fgThreshold = 1;
Int_t    AliEMCALRawUtils::fgDDLPerSuperModule = 2;  // 2 ddls per SuperModule
Int_t    AliEMCALRawUtils::fgPedestalValue = 32;     // pedestal value for digits2raw
Double_t AliEMCALRawUtils::fgFEENoise = 3.;          // 3 ADC channels of noise (sampled)

AliEMCALRawUtils::AliEMCALRawUtils()
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(0), fOption("")
{

  //These are default parameters.  
  //Can be re-set from without with setter functions
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits) 
  fOrder = 2;                         // order of gamma fn
  fTau = 2.35;                        // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold = 3; // 3 ADC counts is approx. noise level
  fNPedSamples = 4;    // less than this value => likely pedestal samples

  //Get Mapping RCU files from the AliEMCALRecParam                                 
  const TObjArray* maps = AliEMCALRecParam::GetMappings();
  if(!maps) AliFatal("Cannot retrieve ALTRO mappings!!");

  for(Int_t i = 0; i < 4; i++) {
    fMapping[i] = (AliAltroMapping*)maps->At(i);
  }

  //To make sure we match with the geometry in a simulation file,
  //let's try to get it first.  If not, take the default geometry
  AliRunLoader *rl = AliRunLoader::Instance();
  if(!rl) AliError("Cannot find RunLoader!");
  if (rl->GetAliRun() && rl->GetAliRun()->GetDetector("EMCAL")) {
    fGeom = dynamic_cast<AliEMCAL*>(rl->GetAliRun()->GetDetector("EMCAL"))->GetGeometry();
  } else {
    AliInfo(Form("Using default geometry in raw reco"));
    fGeom =  AliEMCALGeometry::GetInstance(AliEMCALGeometry::GetDefaultGeometryName());
  }

  if(!fGeom) AliFatal(Form("Could not get geometry!"));

}

//____________________________________________________________________________
AliEMCALRawUtils::AliEMCALRawUtils(AliEMCALGeometry *pGeometry)
  : fHighLowGainFactor(0.), fOrder(0), fTau(0.), fNoiseThreshold(0),
    fNPedSamples(0), fGeom(pGeometry), fOption("")
{
  //
  // Initialize with the given geometry - constructor required by HLT
  // HLT does not use/support AliRunLoader(s) instances
  // This is a minimum intervention solution
  // Comment by MPloskon@lbl.gov
  //

  //These are default parameters. 
  //Can be re-set from without with setter functions 
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits)
  fOrder = 2;                         // order of gamma fn
  fTau = 2.35;                        // in units of timebin, from CERN 2007 testbeam
  fNoiseThreshold = 3; // 3 ADC counts is approx. noise level
  fNPedSamples = 4;    // less than this value => likely pedestal samples

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
    fOption(rawU.fOption)
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
    fOrder = rawU.fOrder;
    fTau = rawU.fTau;
    fNoiseThreshold = rawU.fNoiseThreshold;
    fNPedSamples = rawU.fNPedSamples;
    fGeom = rawU.fGeom;
    fOption = rawU.fOption;
    fMapping[0] = rawU.fMapping[0];
    fMapping[1] = rawU.fMapping[1];
    fMapping[2] = rawU.fMapping[2];
    fMapping[3] = rawU.fMapping[3];
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
    if (digit->GetAmp() < fgThreshold) 
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
      buffers[iDDL]->FillBuffer((Int_t)digit->GetAmp());
      buffers[iDDL]->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffers[iDDL]->FillBuffer(3);          // bunch length      
      buffers[iDDL]->WriteTrailer(3, ieta, iphi, nSM);  // trailer
      // calculate the time response function
    } else {
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), digit->GetAmp(), adcValuesHigh.GetArray(), adcValuesLow.GetArray()) ; 
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
void AliEMCALRawUtils::Raw2Digits(AliRawReader* reader,TClonesArray *digitsArr, AliCaloCalibPedestal* pedbadmap)
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

  //Updated fitting routine from 2007 beam test takes into account
  //possibility of two peaks in data and selects first one for fitting
  //Also sets some of the starting parameters based on the shape of the
  //given raw signal being fit

  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, GetRawFormatTimeBins(), 5);
  signalF->SetParameters(10.,5.,fTau,fOrder,0.); //set all defaults once, just to be safe
  signalF->SetParNames("amp","t0","tau","N","ped");
  signalF->FixParameter(2,fTau); // tau in units of time bin
  signalF->FixParameter(3,fOrder); // order
  
  Int_t id =  -1;
  Float_t time = 0. ; 
  Float_t amp = 0. ; 
  Float_t ped = 0. ;
  Float_t ampEstimate  = 0;
  Float_t timeEstimate = 0;
  Float_t pedEstimate = 0;
  Int_t i = 0;
  Int_t startBin = 0;

  //Graph to hold data we will fit (should be converted to an array
  //later to speed up processing
  TGraph * gSig = new TGraph(GetRawFormatTimeBins()); 

  Int_t lowGain = 0;
  Int_t caloFlag = 0; // low, high gain, or TRU, or LED ref.

  // start loop over input stream 
  while (in.NextDDL()) {
    while (in.NextChannel()) {

      //Check if the signal  is high or low gain and then do the fit, 
      //if it  is from TRU do not fit
      caloFlag = in.GetCaloFlag();
      if (caloFlag != 0 && caloFlag != 1) continue; 
	      
      //Do not fit bad channels
      if(pedbadmap->IsBadChannel(in.GetModule(),in.GetColumn(),in.GetRow())) {
	//printf("Tower from SM %d, column %d, row %d is BAD!!! Skip \n", in.GetModule(),in.GetColumn(),in.GetRow());
	continue;
      }  

      // There can be zero-suppression in the raw data, 
      // so set up the TGraph in advance
      for (i=0; i < GetRawFormatTimeBins(); i++) {
	gSig->SetPoint(i, i , -1); // init to out-of-range values
      }

      Int_t maxTimeBin = 0;
      Int_t min = 0x3ff; // init to 10-bit max
      Int_t max = 0; // init to 10-bit min
      while (in.NextBunch()) {

	const UShort_t *sig = in.GetSignals();
	startBin = in.GetStartTimeBin();
	if (maxTimeBin < startBin) {
	  maxTimeBin = startBin; // timebins come in reverse order
	}	
	if (maxTimeBin < 0 || maxTimeBin >= GetRawFormatTimeBins()) {
	  AliWarning(Form("Invalid time bin %d",maxTimeBin));
	  maxTimeBin = GetRawFormatTimeBins();
	}
	
	for (i = 0; i < in.GetBunchLength(); i++) {
	  time = startBin--;
	  gSig->SetPoint((Int_t)time, time, (Double_t) sig[i]) ;
	  if (max < sig[i]) max = sig[i];
	  if (min > sig[i]) min = sig[i];
	  
	}
      } // loop over bunches

      gSig->Set(maxTimeBin+1); // set actual max size of TGraph
      
      //Initialize the variables, do not keep previous values.
      // not really necessary to reset all of them (only amp and time at the moment), but better safe than sorry
      amp  = -1 ;
      time = -1 ;
      ped = -1;
      ampEstimate  = -1 ;
      timeEstimate = -1 ;
      pedEstimate = -1;
      if ( (max - min) > fNoiseThreshold) {
	FitRaw(gSig, signalF, maxTimeBin, amp, time, ped,
	       ampEstimate, timeEstimate, pedEstimate);
      }
           
      if ( amp>0 && amp<2000 && time>0 && time<(maxTimeBin*GetRawFormatTimeBinWidth()) ) {  //check both high and low end of amplitude result, and time
	//2000 is somewhat arbitrary - not nice with magic numbers in the code..
	id =  fGeom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
	lowGain = in.IsLowGain();

	// check fit results: should be consistent with initial estimates
	// more magic numbers, but very loose cuts, for now..
	// We have checked that amp an time values are positive so division for assymmetry
	// calculation should be OK/safe
	Float_t ampAsymm = (amp - ampEstimate)/(amp + ampEstimate);
	if ( (TMath::Abs(ampAsymm) > 0.1) ||
	     (TMath::Abs(time - timeEstimate) > 2*GetRawFormatTimeBinWidth()) ) {
	  AliDebug(2,Form("Fit results ped %f amp %f time %f not consistent with expectations ped %f max-ped %f time %d",
		      ped, amp, time, pedEstimate, ampEstimate, timeEstimate));

	  // what should do we do then? skip this channel or assign the simple estimate? 
	  // for now just overwrite the fit results with the simple estimate
	  amp = ampEstimate;
	  time = timeEstimate; 
	}

	AliDebug(2,Form("id %d lowGain %d amp %g", id, lowGain, amp));
	// printf("Added tower: SM %d, row %d, column %d, amp %3.2f\n",in.GetModule(), in.GetRow(), in.GetColumn(),amp);
	// round off amplitude value to nearest integer
	AddDigit(digitsArr, id, lowGain, TMath::Nint(amp), time); 
      }
      
      // Reset graph
      for (Int_t index = 0; index < gSig->GetN(); index++) {
	gSig->SetPoint(index, index, -1) ;  
      } 
      // Reset starting parameters for fit function
      signalF->SetParameters(10.,5.,fTau,fOrder,0.); //reset all defaults just to be safe

   } // end while over channel   
  } //end while over DDL's, of input stream 
  
  delete signalF ; 
  delete gSig;
  
  return ; 
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::AddDigit(TClonesArray *digitsArr, Int_t id, Int_t lowGain, Int_t amp, Float_t time) {
  //
  // Add a new digit. 
  // This routine checks whether a digit exists already for this tower 
  // and then decides whether to use the high or low gain info
  //
  // Called by Raw2Digits
  
  AliEMCALDigit *digit = 0, *tmpdigit = 0;
  TIter nextdigit(digitsArr);
  while (digit == 0 && (tmpdigit = (AliEMCALDigit*) nextdigit())) {
    if (tmpdigit->GetId() == id)
      digit = tmpdigit;
  }

  if (!digit) { // no digit existed for this tower; create one
    if (lowGain && amp > fgkOverflowCut) 
      amp = Int_t(fHighLowGainFactor * amp); 
    Int_t idigit = digitsArr->GetEntries();
    new((*digitsArr)[idigit]) AliEMCALDigit( -1, -1, id, amp, time, idigit) ;	
  }
  else { // a digit already exists, check range 
         // (use high gain if signal < cut value, otherwise low gain)
    if (lowGain) { // new digit is low gain
      if (digit->GetAmp() > fgkOverflowCut) {  // use if stored digit is out of range
	digit->SetAmp(Int_t(fHighLowGainFactor * amp));
	digit->SetTime(time);
      }
    }
    else if (amp < fgkOverflowCut) { // new digit is high gain; use if not out of range
      digit->SetAmp(amp);
      digit->SetTime(time);
    }
  }
}

//____________________________________________________________________________ 
void AliEMCALRawUtils::FitRaw(TGraph * gSig, TF1* signalF, const Int_t lastTimeBin, Float_t & amp, Float_t & time, Float_t & ped, Float_t & ampEstimate, Float_t & timeEstimate, Float_t & pedEstimate, const Float_t cut) const 
{
  // Fits the raw signal time distribution; from AliEMCALGetter 
  // last argument: Float_t cut = 0.0; // indicating how much of amplitude w.r.t. max value fit should be above noise and pedestal 

  // initialize return values
  amp = 0; 
  time = 0; 
  ped = 0;
  ampEstimate = 0;
  timeEstimate = 0;
  pedEstimate = 0;

  // 0th step: remove plateau / overflow candidates
  // before trying to estimate amplitude, search for maxima etc.
  //
  Int_t nOrig = gSig->GetN(); // number of samples before we remove any overflows
  // Values for readback from input graph
  Double_t ttime = 0;
  Double_t signal = 0;

  /*
  // start: tmp dump of all values
  for (Int_t i=0; i<gSig->GetN(); i++) {
    gSig->GetPoint(i, ttime, signal) ; // get values
    printf("orig: i %d, time %f, signal %f\n",i, ttime, signal);
  }
  // end: tmp dump of all values
  */

  // start from back of TGraph since RemovePoint will downshift indices
  for (Int_t i=nOrig-1; i>=0; i--) {
    gSig->GetPoint(i, ttime, signal) ; // get values
    if (signal >= (pedEstimate + fgkOverflowCut) ) {
      gSig->RemovePoint(i);
    }
  }

  // 1st step: we try to estimate the pedestal value
  Int_t nPed = 0;
  for (Int_t index = 0; index < gSig->GetN(); index++) {
    gSig->GetPoint(index, ttime, signal) ; 
    // ttime < fNPedsamples used for pedestal estimate; 
    // ttime >= fNPedSamples used for signal checks
    if (signal >= 0 && ttime<fNPedSamples) { // valid value
      pedEstimate += signal;
      nPed++;
    }
  }

  if (nPed > 0)
    pedEstimate /= nPed;
  else {
    //AliWarning("Could not determine pedestal");	  
    AliDebug(1,"Could not determine pedestal");
    pedEstimate = 0; // good estimate for ZeroSupp data (non ZS data should have no problem with pedestal estimate)
  }

  // 2nd step: we look through the rest of the time-bins/ADC values and
  // see if we have something that looks like a signal.
  // We look for a first local maxima, as well as for a global maxima 
  Int_t locMaxFound = 0;
  Int_t locMaxId = 0; // time-bin index at first local max
  Float_t locMaxSig = -1; // actual local max value
  Int_t globMaxId = 0; // time-bin index at global max
  Float_t globMaxSig = -1; // actual global max value
  // We will also look for any values that look like they are in overflow region
  for (Int_t i=0; i<gSig->GetN(); i++) {
    gSig->GetPoint(i, ttime, signal) ; // get values

    // ttime < fNPedsamples used for pedestal estimate; 
    // ttime >= fNPedSamples used for signal checks
    if (ttime >= fNPedSamples) { 

      // look for first local maximum signal=ADC value
      if (!locMaxFound && signal > locMaxSig) {
	locMaxId = i;
	locMaxSig = signal;
      }
      else if ( locMaxSig > (pedEstimate + fNoiseThreshold) ) { 
	// we enter this condition after signal<=max, but previous
	// max value was large enough. I.e. at least a significant local 
	// maxima has been found (just before)
	locMaxFound = 1;
      }

      // also check for global maximum..
      if (signal > globMaxSig) {
	globMaxId = i;
	globMaxSig = signal;
      }
    } // ttime check
  } // end for-loop over samples after pedestal

  // OK, we have looked through the signal spectra, let's see if we should try to make the fit
  ampEstimate = locMaxSig - pedEstimate; // estimate using first local maxima 
  if ( ampEstimate > fNoiseThreshold ) { // else it's just noise 

    //Check that the local maximum we will use is not at the end or beginning of time sample range
    Double_t timeMax = -1;
    Int_t iMax = locMaxId;
    gSig->GetPoint(locMaxId, timeMax, signal) ;
    if (timeMax < 2 || timeMax > lastTimeBin-1) { // lastTimeBin is the lowest kept time-sample; current (Dec 2009) case
      //    if (timeMax < 2 || timeMax > lastTimeBin-2) { // for when lastTimeBin is the lowest read-out time-sample, future (2010) case
      AliDebug(1,Form("Skip fit, maximum of the sample close to the edges : timeMax %3.2f, ampEstimate %3.2f",timeMax, ampEstimate));
      return;
    }

    // Check if the local and global maximum disagree
    if (locMaxId != globMaxId) {
      AliDebug(1,Form("Warning, local first maximum %d does not agree with global maximum %d\n", locMaxId, globMaxId));
      return;
    }
    
    // Get the maximum and find the lowest timebin (tailmin) where the ADC value is not 
    // significantly different from the pedestal
    // first lower times edge a.k.a. tailmin
    Int_t tailMin = 0;
    Double_t tmptime = 0;
    for (Int_t i=iMax-1; i > 0; i--) {
      gSig->GetPoint(i, tmptime, signal) ;
      if((signal-pedEstimate) < fNoiseThreshold){
	tailMin = i;
	break;
      }
    }
    // then same exercise for the higher times edge a.k.a. tailmax
    Int_t tailMax = lastTimeBin;
    for (Int_t i=iMax+1; i < gSig->GetN(); i++) {
      gSig->GetPoint(i, tmptime, signal) ;
      if ((signal-pedEstimate) <= (ampEstimate*cut + fNoiseThreshold)) { // stop fit at cut-fraction of amplitude above noise-threshold (cut>0 would mean avoid the pulse shape falling edge)
	tailMax = i;
	break;
      }
    }

    // remove all points which are not in the distribution around maximum
    // i.e. up to tailmin, and from tailmax
    if ( tailMax != (gSig->GetN()-1) ){ // else nothing to remove
      nOrig = gSig->GetN(); // can't use GetN call in for loop below since gSig size changes..
      for(int j = tailMax; j < nOrig; j++) gSig->RemovePoint(tailMax);
    }
    for(int j = 0; j<=tailMin; j++) gSig->RemovePoint(0);

    if(gSig->GetN() < 3) {
      AliDebug(2,Form("Skip fit, number of entries in sample smaller than number of fitting parameters: in sample %d, fitting param 3", 
		      gSig->GetN() ));
      return;
    }

    timeEstimate = timeMax * GetRawFormatTimeBinWidth();

    // determine what the valid fit range is
    Double_t minFit = 9999;
    Double_t maxFit = 0;
    for (Int_t i=0; i < gSig->GetN(); i++) {
      gSig->GetPoint(i, ttime, signal); 
      if (minFit > ttime) minFit=ttime;
      if (maxFit < ttime) maxFit=ttime;
      //debug: printf("no tail: i %d, time %f, signal %f\n",i, ttime, signal); 
    } 
    signalF->SetRange(minFit, maxFit);

    signalF->FixParameter(4, pedEstimate) ; 
    signalF->SetParameter(1, timeMax);
    signalF->SetParameter(0, ampEstimate);
    
    gSig->Fit(signalF, "QROW"); // Note option 'W': equal errors on all points

    // assign fit results
    amp = signalF->GetParameter(0); 
    time = signalF->GetParameter(1) * GetRawFormatTimeBinWidth(); // skip subtraction of fgTimeTrigger?
    ped = signalF->GetParameter(4); 

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
      }
    }

    Double_t diff = TMath::Abs(1-ymMaxY[0]/amp) ; 
    if (diff > 0.1) 
      amp = ymMaxY[0] ; 

      //END YS

  } // ampEstimate > fNoiseThreshold
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
Bool_t AliEMCALRawUtils::RawSampledResponse(
const Double_t dtime, const Double_t damp, Int_t * adcH, Int_t * adcL) const 
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

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    Double_t signal = signalF.Eval(iTime) ;     

    // Next lines commeted for the moment but in principle it is not necessary to add
    // extra noise since noise already added at the digits level.	

    //According to Terry Awes, 13-Apr-2008
    //add gaussian noise in quadrature to each sample
    //Double_t noise = gRandom->Gaus(0.,fgFEENoise);
    //signal = sqrt(signal*signal + noise*noise);

    // March 17,09 for fast fit simulations by Alexei Pavlinov.
    // Get from PHOS analysis. In some sense it is open questions.
    //Double_t noise = gRandom->Gaus(0.,fgFEENoise);
    //signal += noise; 

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
