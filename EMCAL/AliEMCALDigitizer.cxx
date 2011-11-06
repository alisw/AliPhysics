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
// 
//////////////////////////////////////////////////////////////////////////////
// Class performs digitization of Summable digits from simulated data
//  
// In addition it performs mixing/embedding of summable digits from different events.
//
// For each event 3 branches are created in TreeD:
//   "EMCAL" - list of digits
//   "EMCALTRG" - list of trigger digits
//   "AliEMCALDigitizer" - AliEMCALDigitizer with all parameters used in digitization
//
//
////////////////////////////////////////////////////////////////////////////////////
//
//*-- Author: Sahal Yacoob (LBL)
// based on : AliEMCALDigitizer
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new IO (a la PHOS)
//  November 2003 Aleksei Pavlinov : adopted for Shish-Kebab geometry 
//  July 2011 GCB: Digitizer modified to accomodate embedding. 
//                 Time calibration added. Decalibration possibility of energy and time added
//_________________________________________________________________________________

// --- ROOT system ---
#include <TROOT.h>
#include <TTree.h>
#include <TSystem.h>
#include <TBenchmark.h>
#include <TBrowser.h>
#include <TObjectTable.h>
#include <TRandom.h>
#include <TF1.h>
#include <cassert>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRun.h"
#include "AliDigitizationInput.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALDigit.h"
#include "AliEMCAL.h"
#include "AliEMCALLoader.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALCalibData.h"
#include "AliEMCALSimParam.h"
#include "AliEMCALRawDigit.h"

  namespace
  {
    Double_t HeavisideTheta(Double_t x)
    {
      Double_t signal = 0.;
      
      if (x > 0.) signal = 1.;  
      
      return signal;  
    }
    
    Double_t AnalogFastORFunction(Double_t *x, Double_t *par)
    {
      Double_t v0 = par[0];
      Double_t t0 = par[1];
      Double_t tr = par[2];
      
      Double_t R1 = 1000.;
      Double_t C1 = 33e-12;
      Double_t R2 = 1800;
      Double_t C2 = 22e-12;
      
      Double_t t  =   x[0];
      
      return (((0.8*(-((TMath::Power(C1,2)*C2*TMath::Power(TMath::E(),(-t + t0)/(C1*R1))*
			TMath::Power(R1,2)*R2)/(C1*R1 - C2*R2)) + 
		     C1*C2*R1*R2*(1 - (C2*TMath::Power(TMath::E(),(-t + t0)/(C2*R2))*R2)/(-(C1*R1) + C2*R2)))*v0*
		HeavisideTheta(t - t0))/tr 
	       - (0.8*(C1*C2*R1*R2 - 
		       (TMath::Power(C1,2)*C2*TMath::Power(TMath::E(),(-1.*t + t0 + 1.25*tr)/(C1*R1))*
			TMath::Power(R1,2)*R2)/(C1*R1 - C2*R2) + 
		       (C1*TMath::Power(C2,2)*TMath::Power(TMath::E(),(-1.*t + t0 + 1.25*tr)/(C2*R2))*
			R1*TMath::Power(R2,2))/(C1*R1 - C2*R2))*v0*
		  HeavisideTheta(t - t0 - 1.25*tr))/tr)/(C2*R1));
    }
  }

ClassImp(AliEMCALDigitizer)
  
  
//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer()
  : AliDigitizer("",""),
    fDefaultInit(kTRUE),
    fDigitsInRun(0),
    fInit(0),
    fInput(0),
    fInputFileNames(0x0),
    fEventNames(0x0),
    fDigitThreshold(0),
    fMeanPhotonElectron(0),
    fGainFluctuations(0),
    fPinNoise(0),
    fTimeNoise(0),
    fTimeDelay(0),
    fTimeResolutionPar0(0),
    fTimeResolutionPar1(0),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fADCchannelECDecal(0),
    fTimeChannel(0),
    fTimeChannelDecal(0),
    fNADCEC(0),
    fEventFolderName(""),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0),
    fSDigitizer(0x0)
{
  // ctor
  InitParameters() ; 
  fDigInput = 0 ;                     // We work in the standalone mode
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(TString alirunFileName, TString eventFolderName)
  : AliDigitizer("EMCALDigitizer", alirunFileName),
    fDefaultInit(kFALSE),
    fDigitsInRun(0),
    fInit(0),
    fInput(0),
    fInputFileNames(0), 
    fEventNames(0), 
    fDigitThreshold(0),
    fMeanPhotonElectron(0),
    fGainFluctuations(0),
    fPinNoise(0),
    fTimeNoise(0),
    fTimeDelay(0),
    fTimeResolutionPar0(0),
    fTimeResolutionPar1(0),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fADCchannelECDecal(0),
    fTimeChannel(0),
    fTimeChannelDecal(0),
    fNADCEC(0),
    fEventFolderName(eventFolderName),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0),
    fSDigitizer(0x0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fDigInput = 0 ;                     // We work in the standalone mode
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(const AliEMCALDigitizer & d) 
  : AliDigitizer(d.GetName(),d.GetTitle()),
    fDefaultInit(d.fDefaultInit),
    fDigitsInRun(d.fDigitsInRun),
    fInit(d.fInit),
    fInput(d.fInput),
    fInputFileNames(d.fInputFileNames),
    fEventNames(d.fEventNames),
    fDigitThreshold(d.fDigitThreshold),
    fMeanPhotonElectron(d.fMeanPhotonElectron),
    fGainFluctuations(d.fGainFluctuations),
    fPinNoise(d.fPinNoise),
    fTimeNoise(d.fTimeNoise),
    fTimeDelay(d.fTimeDelay),
    fTimeResolutionPar0(d.fTimeResolutionPar0),
    fTimeResolutionPar1(d.fTimeResolutionPar1),
    fADCchannelEC(d.fADCchannelEC),
    fADCpedestalEC(d.fADCpedestalEC),
    fADCchannelECDecal(d.fADCchannelECDecal),
    fTimeChannel(d.fTimeChannel), fTimeChannelDecal(d.fTimeChannelDecal),
    fNADCEC(d.fNADCEC),
    fEventFolderName(d.fEventFolderName),
    fFirstEvent(d.fFirstEvent),
    fLastEvent(d.fLastEvent),
    fCalibData(d.fCalibData),
    fSDigitizer(d.fSDigitizer ? new AliEMCALSDigitizer(*d.fSDigitizer) : 0)
{
  // copyy ctor 
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(AliDigitizationInput * rd)
  : AliDigitizer(rd,"EMCALDigitizer"),
    fDefaultInit(kFALSE),
    fDigitsInRun(0),
    fInit(0),
    fInput(0),
    fInputFileNames(0),
    fEventNames(0),
    fDigitThreshold(0),
    fMeanPhotonElectron(0),
    fGainFluctuations(0),
    fPinNoise(0.),
    fTimeNoise(0.),
    fTimeDelay(0.),
    fTimeResolutionPar0(0.),
    fTimeResolutionPar1(0.),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fADCchannelECDecal(0),
    fTimeChannel(0),
    fTimeChannelDecal(0),    
    fNADCEC(0),
    fEventFolderName(0),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0),
    fSDigitizer(0x0)
{
  // ctor Init() is called by RunDigitizer
  fDigInput = rd ; 
  fEventFolderName = fDigInput->GetInputFolderName(0) ;
  SetTitle(dynamic_cast<AliStream*>(fDigInput->GetInputStream(0))->GetFileName(0));
  InitParameters() ; 
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  //dtor
  delete [] fInputFileNames ; 
  delete [] fEventNames ; 
  if (fSDigitizer) delete fSDigitizer;
}

//____________________________________________________________________________
void AliEMCALDigitizer::Digitize(Int_t event) 
{ 
  
  // Makes the digitization of the collected summable digits
  // for this it first creates the array of all EMCAL modules
  // filled with noise and after that adds contributions from 
  // SDigits. This design helps to avoid scanning over the 
  // list of digits to add  contribution of any new SDigit.
  //
  // JLK 26-Jun-2008
  // Note that SDigit energy info is stored as an amplitude, so we
  // must call the Calibrate() method of the SDigitizer to convert it
  // back to an energy in GeV before adding it to the Digit
  //
  static int nEMC=0; //max number of digits possible
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  
  if(emcalLoader){
    Int_t readEvent = event ; 
    // fDigInput is data member from AliDigitizer
    if (fDigInput) 
      readEvent = dynamic_cast<AliStream*>(fDigInput->GetInputStream(0))->GetCurrentEventNumber() ; 
    AliDebug(1,Form("Adding event %d from input stream 0 %s %s", 
                    readEvent, GetTitle(), fEventFolderName.Data())) ; 
    rl->GetEvent(readEvent);
    
    TClonesArray * digits = emcalLoader->Digits() ; 
    digits->Delete() ;  //correct way to clear array when memory is
    //allocated by objects stored in array
    
    // Load Geometry
    AliEMCALGeometry *geom = 0;
    if (!rl->GetAliRun()) {
      AliFatal("Could not get AliRun from runLoader");
    }
    else{
      AliEMCAL * emcal  = (AliEMCAL*)rl->GetAliRun()->GetDetector("EMCAL");
      geom = emcal->GetGeometry();
      nEMC = geom->GetNCells();
      AliDebug(1,Form("nEMC %i (number cells in EMCAL) | %s \n", nEMC, geom->GetName()));
    }
    
    Int_t absID = -1 ;
    
    digits->Expand(nEMC) ;
    
    // RS create a digitizer on the fly
    if (!fSDigitizer) fSDigitizer = new AliEMCALSDigitizer(rl->GetFileName().Data());
    fSDigitizer->SetEventRange(0, -1) ;
    //
    //take all the inputs to add together and load the SDigits
    TObjArray * sdigArray = new TObjArray(fInput) ;
    sdigArray->AddAt(emcalLoader->SDigits(), 0) ;
    Int_t i ;
    Int_t embed = kFALSE;
    for(i = 1 ; i < fInput ; i++){
      TString tempo(fEventNames[i]) ; 
      tempo += i ;
        
      AliRunLoader *  rl2 = AliRunLoader::GetRunLoader(tempo) ; 
        
      if (rl2==0) 
	rl2 = AliRunLoader::Open(fInputFileNames[i], tempo) ; 
        
      if (fDigInput) 
	readEvent = dynamic_cast<AliStream*>(fDigInput->GetInputStream(i))->GetCurrentEventNumber() ; 
      Info("Digitize", "Adding event %d from input stream %d %s %s", readEvent, i, fInputFileNames[i].Data(), tempo.Data()) ; 
      rl2->LoadSDigits();
      //     rl2->LoadDigits();
      rl2->GetEvent(readEvent);
      if(rl2->GetDetectorLoader("EMCAL"))
        {
          AliEMCALLoader *emcalLoader2 = dynamic_cast<AliEMCALLoader*>(rl2->GetDetectorLoader("EMCAL"));
 
          if(emcalLoader2) {
            //Merge 2 simulated sdigits
            if(emcalLoader2->SDigits()){
              TClonesArray* sdigits2 = emcalLoader2->SDigits();
              sdigArray->AddAt(sdigits2, i) ;
              // Check if first sdigit is of embedded type, if so, handle the sdigits differently:
              // do not smear energy of embedded, do not add noise to any sdigits
              if(sdigits2->GetEntriesFast()>0){
		//printf("Merged digit type: %d\n",dynamic_cast<AliEMCALDigit*> (sdigits2->At(0))->GetType());
		AliEMCALDigit * digit2 = dynamic_cast<AliEMCALDigit*> (sdigits2->At(0));
		if(digit2 && digit2->GetType()==AliEMCALDigit::kEmbedded){
		  embed = kTRUE;
		}
              }
            }
            
          }//loader2
          else  AliFatal("EMCAL Loader of second event not available!");
          
        }
      else Fatal("Digitize", "Loader of second input not found");
    }// input loop
      
    //Find the first tower with signal
    Int_t nextSig = nEMC + 1 ; 
    TClonesArray * sdigits ;  
    for(i = 0 ; i < fInput ; i++){
      if(i > 0 && embed) continue;
      sdigits = dynamic_cast<TClonesArray *>(sdigArray->At(i)) ;
      if(sdigits){
	if (sdigits->GetEntriesFast() ){
	  AliEMCALDigit *sd = dynamic_cast<AliEMCALDigit *>(sdigits->At(0));
	  if(sd){
	    Int_t curNext = sd->GetId() ;
	    if(curNext < nextSig) 
	      nextSig = curNext ;
	    AliDebug(1,Form("input %i : #sdigits %i \n",
			    i, sdigits->GetEntriesFast()));
              
	  }//First entry is not NULL
	  else {
	    Info("Digitize", "NULL sdigit pointer");
	    continue;
	  }
	}//There is at least one entry
	else {
	  Info("Digitize", "NULL sdigits array");
	  continue;
	}
      }// SDigits array exists
      else {
	Info("Digitizer","Null sdigit pointer");
	continue;
      }
    }// input loop
    AliDebug(1,Form("FIRST tower with signal %i \n", nextSig));
      
    TArrayI index(fInput) ;
    index.Reset() ;  //Set all indexes to zero
      
    AliEMCALDigit * digit ;
    AliEMCALDigit * curSDigit ;
      
    //Put Noise contribution, smear time and energy     
    Float_t timeResolution = 0;
    for(absID = 0; absID < nEMC; absID++){ // Nov 30, 2006 by PAI; was from 1 to nEMC
      Float_t energy = 0 ;    
        
      // amplitude set to zero, noise will be added later
      Float_t noiseTime = 0.;
      if(!embed) noiseTime = TimeOfNoise(); //No need for embedded events?
      new((*digits)[absID]) AliEMCALDigit( -1, -1, absID, 0., noiseTime,kFALSE); // absID-1->absID
      //look if we have to add signal?
      digit = dynamic_cast<AliEMCALDigit *>(digits->At(absID)); // absID-1->absID
        
      if (digit) {
          
	if(absID==nextSig){

	  // Calculate time as time of the largest digit
	  Float_t time = digit->GetTime() ;
	  Float_t aTime= digit->GetAmplitude() ;
            
	  // loop over input
	  Int_t nInputs = fInput;
	  if(embed) nInputs = 1; // In case of embedding, merge later real digits, do not smear energy and time of data
	  for(i = 0; i< nInputs ; i++){  //loop over (possible) merge sources
	    TClonesArray* sdtclarr = dynamic_cast<TClonesArray *>(sdigArray->At(i));
	    if(sdtclarr) {
	      Int_t sDigitEntries = sdtclarr->GetEntriesFast();
                
	      if(sDigitEntries > index[i] )
		curSDigit = dynamic_cast<AliEMCALDigit*>(sdtclarr->At(index[i])) ; 	
	      else
		curSDigit = 0 ;
                
	      //May be several digits will contribute from the same input
	      while(curSDigit && (curSDigit->GetId() == absID)){	   
		//Shift primary to separate primaries belonging different inputs
		Int_t primaryoffset ;
		if(fDigInput)
		  primaryoffset = fDigInput->GetMask(i) ; 
		else
		  primaryoffset = i ;
		curSDigit->ShiftPrimary(primaryoffset) ;
                  
		if(curSDigit->GetAmplitude()>aTime) {
		  aTime = curSDigit->GetAmplitude();
		  time  = curSDigit->GetTime();
		}
                  
		*digit = *digit + *curSDigit ;  //adds amplitudes of each digit
                  
		index[i]++ ;
                  
		if( sDigitEntries > index[i] )
		  curSDigit = dynamic_cast<AliEMCALDigit*>(sdtclarr->At(index[i])) ;
		else
		  curSDigit = 0 ;
	      }// while
	    }// source exists
	  }// loop over merging sources
            
            
	  //Here we convert the summed amplitude to an energy in GeV only for simulation or mixing of simulations
	  energy = fSDigitizer->Calibrate(digit->GetAmplitude()) ; // GeV
            
	  // add fluctuations for photo-electron creation
    // corrected fluctuations after comparison with beam test, Paraskevi Ganoti (06/11/2011)
    Float_t fluct = static_cast<Float_t>((energy*fMeanPhotonElectron)/fGainFluctuations);
    energy       *= static_cast<Float_t>(gRandom->Poisson(fluct)) / fluct ;    
    
	  //calculate and set time
	  digit->SetTime(time) ;
            
	  //Find next signal module
	  nextSig = nEMC + 1 ;
	  for(i = 0 ; i < nInputs ; i++){
	    sdigits = dynamic_cast<TClonesArray *>(sdigArray->At(i)) ;
              
	    if(sdigits){
	      Int_t curNext = nextSig ;
	      if(sdigits->GetEntriesFast() > index[i])
                {
                  AliEMCALDigit * tmpdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(index[i]));
                  if(tmpdigit)
		    {
		      curNext = tmpdigit->GetId() ;	  
		    }
                }
	      if(curNext < nextSig) nextSig = curNext ;
	    }// sdigits exist
	  } // input loop       
            
	}//absID==nextSig
          
	// add the noise now, no need for embedded events with real data
	if(!embed)
	  energy += TMath::Abs(gRandom->Gaus(0., fPinNoise)) ;
          

	// JLK 26-June-2008
	//Now digitize the energy according to the fSDigitizer method,
	//which merely converts the energy to an integer.  Later we will
	//check that the stored value matches our allowed dynamic ranges
	digit->SetAmplitude(fSDigitizer->Digitize(energy)) ; 
          
	//Set time resolution with final energy
	timeResolution = GetTimeResolution(energy);
	digit->SetTime(gRandom->Gaus(digit->GetTime(),timeResolution) ) ;
	AliDebug(10,Form(" absID %5i energy %f nextSig %5i\n",
			 absID, energy, nextSig));
	//Add delay to time
	digit->SetTime(digit->GetTime()+fTimeDelay) ;
          
      }// Digit pointer exists
      else{
	Info("Digitizer","Digit pointer is null");
      }
    } // for(absID = 0; absID < nEMC; absID++)
      
      //ticks->Delete() ;
      //delete ticks ;

      //Embed simulated amplitude (and time?) to real data digits
    if(embed){
      for(Int_t input = 1; input<fInput; input++){
	TClonesArray *realDigits = dynamic_cast<TClonesArray*> (sdigArray->At(input));
	if(!realDigits){
	  AliDebug(1,"No sdigits to merge\n");
	  continue;
	}
	for(Int_t i2 = 0 ; i2 < realDigits->GetEntriesFast() ; i2++){
	  AliEMCALDigit * digit2 = dynamic_cast<AliEMCALDigit*>( realDigits->At(i2) ) ;
	  if(digit2){
	    digit = dynamic_cast<AliEMCALDigit*>( digits->At(digit2->GetId()) ) ;
	    if(digit){
                
	      // Put the embedded cell energy in same units as simulated sdigits -> transform to energy units
	      // and multiply to get a big int.
	      Float_t time2 = digit2->GetTime(); 
	      Float_t e2    = digit2->GetAmplitude();
	      CalibrateADCTime(e2,time2,digit2->GetId());
	      Float_t amp2 = fSDigitizer->Digitize(e2);
                
	      Float_t e0    = digit ->GetAmplitude();
	      if(e0>0.01)
		AliDebug(1,Form("digit 1: Abs ID %d, amp %f, type %d, time %e; digit2: Abs ID %d, amp %f, type %d, time %e\n",
				digit ->GetId(),digit ->GetAmplitude(), digit ->GetType(), digit->GetTime(),
				digit2->GetId(),amp2,                   digit2->GetType(), time2           ));
                
	      // Sum amplitudes, change digit amplitude
	      digit->SetAmplitude( digit->GetAmplitude() + amp2);
	      // Rough assumption, assign time of the largest of the 2 digits, 
	      // data or signal to the final digit.
	      if(amp2 > digit->GetAmplitude())  digit->SetTime(time2);
	      //Mark the new digit as embedded
	      digit->SetType(AliEMCALDigit::kEmbedded);
                
	      if(digit2->GetAmplitude()>0.01 && e0> 0.01 )
		AliDebug(1,Form("Embedded digit: Abs ID %d, amp %f, type %d\n",
				digit->GetId(), digit->GetAmplitude(), digit->GetType()));
	    }//digit2 
	  }//digit2
	}//loop digit 2
      }//input loop
    }//embed
      
    //JLK is it better to call Clear() here?
    delete sdigArray ; //We should not delete its contents
      
    //remove digits below thresholds
    // until 10-02-2010 remove digits with energy smaller than fDigitThreshold 3*fPinNoise
    // now, remove digits with Digitized ADC smaller than fDigitThreshold = 3,
    // before merge in the same loop real data digits if available
    Float_t energy = 0;
    Float_t time   = 0;
    for(i = 0 ; i < nEMC ; i++){
      digit = dynamic_cast<AliEMCALDigit*>( digits->At(i) ) ;
      if(digit){
          
	//Then get the energy in GeV units.
	energy = fSDigitizer->Calibrate(digit->GetAmplitude()) ;
	//Then digitize using the calibration constants of the ocdb
	Float_t ampADC = energy;
	DigitizeEnergyTime(ampADC, time, digit->GetId())  ; 	  
	if(ampADC < fDigitThreshold)
	  digits->RemoveAt(i) ;
          
      }// digit exists
    } // digit loop
      
    digits->Compress() ;  
      
    Int_t ndigits = digits->GetEntriesFast() ; 

    //JLK 26-June-2008
    //After we have done the summing and digitizing to create the
    //digits, now we want to calibrate the resulting amplitude to match
    //the dynamic range of our real data.  
    for (i = 0 ; i < ndigits ; i++) { 
      digit = dynamic_cast<AliEMCALDigit *>( digits->At(i) ) ; 
      if(digit){
	digit->SetIndexInList(i) ; 
	energy = fSDigitizer->Calibrate(digit->GetAmplitude()) ;
	time   = digit->GetTime();
	Float_t ampADC = energy;
	DigitizeEnergyTime(ampADC, time, digit->GetId());
	digit->SetAmplitude(ampADC) ;
	digit->SetTime(time);
	// printf("digit amplitude set at end: i %d, amp %f\n",i,digit->GetAmplitude());
      }// digit exists
    }//Digit loop
      
  }
  else AliFatal("EMCALLoader is NULL!") ;
  
}

//_____________________________________________________________________
void AliEMCALDigitizer::DigitizeEnergyTime(Float_t & energy, Float_t & time, const Int_t absId)
{ 
  // JLK 26-June-2008
  // Returns digitized value of the energy in a cell absId
  // using the calibration constants stored in the OCDB
  // or default values if no CalibData object is found.
  // This effectively converts everything to match the dynamic range
  // of the real data we will collect
  //
  // Load Geometry
  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  
  if (geom==0){
    AliFatal("Did not get geometry from EMCALLoader");
    return;
  }
  
  Int_t iSupMod = -1;
  Int_t nModule = -1;
  Int_t nIphi   = -1;
  Int_t nIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Bool_t bCell = geom->GetCellIndex(absId, iSupMod, nModule, nIphi, nIeta) ;
  if(!bCell)
    Error("DigitizeEnergyTime","Wrong cell id number : absId %i ", absId) ;
  geom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
  
  if(fCalibData) {
    fADCpedestalEC     = fCalibData->GetADCpedestal     (iSupMod,ieta,iphi);
    fADCchannelEC      = fCalibData->GetADCchannel      (iSupMod,ieta,iphi);
    fADCchannelECDecal = fCalibData->GetADCchannelDecal (iSupMod,ieta,iphi);
    fTimeChannel       = fCalibData->GetTimeChannel     (iSupMod,ieta,iphi);
    fTimeChannelDecal  = fCalibData->GetTimeChannelDecal(iSupMod,ieta,iphi);      
  }
  
  //Apply calibration to get ADC counts and partial decalibration as especified in OCDB
  energy = (energy + fADCpedestalEC)/fADCchannelEC/fADCchannelECDecal   ;
  time  += fTimeChannel-fTimeChannelDecal;
  
  if(energy > fNADCEC ) 
    energy =  fNADCEC ; 
}

//_____________________________________________________________________
void AliEMCALDigitizer::CalibrateADCTime(Float_t & adc, Float_t & time, const Int_t absId)
{ 
  // Returns the energy in a cell absId with a given adc value
  // using the calibration constants stored in the OCDB. Time also corrected from parameter in OCDB
  // Used in case of embedding, transform ADC counts from real event into energy
  // so that we can add the energy of the simulated sdigits which are in energy
  // units.
  // Same as in AliEMCALClusterizer::Calibrate() but here we do not reject channels being marked as hot
  // or with time out of window
  
  // Load Geometry
  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  
  if (!geom){
    AliFatal("Did not get geometry from EMCALLoader");
    return;
  }
  
  Int_t iSupMod = -1;
  Int_t nModule = -1;
  Int_t nIphi   = -1;
  Int_t nIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Bool_t bCell = geom->GetCellIndex(absId, iSupMod, nModule, nIphi, nIeta) ;
  if(!bCell)
    Error("CalibrateADCTime","Wrong cell id number : absId %i ", absId) ;
  geom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
  
  if(fCalibData) {
    fADCpedestalEC = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);
    fADCchannelEC  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
    fTimeChannel   = fCalibData->GetTimeChannel(iSupMod,ieta,iphi);
  }
  
  adc   = adc * fADCchannelEC - fADCpedestalEC;
  time -= fTimeChannel;
  
  
}


//____________________________________________________________________________
void AliEMCALDigitizer::Digitize(Option_t *option) 
{ 
  // Steering method to process digitization for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1, then process events until the end.
  // by default fLastEvent = fFirstEvent (process only one event)

  if (!fInit) { // to prevent overwrite existing file
    Error( "Digitize", "Give a version name different from %s", fEventFolderName.Data() ) ;
    return ;
  } 

  if (strstr(option,"print")) {

    Print();
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALDigitizer");

  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  Int_t nEvents = 0;
  if(!emcalLoader){
    AliFatal("Did not get the  Loader");
  }
  else{
    
    if (fLastEvent == -1)  {
      fLastEvent = rl->GetNumberOfEvents() - 1 ;
    }
    else if (fDigInput) 
      fLastEvent = fFirstEvent ; // what is this ??
    
    nEvents = fLastEvent - fFirstEvent + 1;
    Int_t ievent  = -1;

    TClonesArray* digitsTMP = new TClonesArray("AliEMCALDigit",    32*96);
    TClonesArray* digitsTRG = new TClonesArray("AliEMCALRawDigit", 32*96);

    rl->LoadSDigits("EMCAL");
    for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
      
      rl->GetEvent(ievent);
      
      Digitize(ievent) ; //Add prepared SDigits to digits and add the noise
      
      WriteDigits() ;
      
      //Trigger Digits
      //-------------------------------------
     

      Digits2FastOR(digitsTMP, digitsTRG);  
      
      WriteDigits(digitsTRG);
      
      (emcalLoader->TreeD())->Fill();
      
      emcalLoader->WriteDigits(   "OVERWRITE");
      
      Unload();
      
      digitsTRG  ->Delete();
      digitsTMP  ->Delete();
      
      //-------------------------------------
      
      if(strstr(option,"deb"))
        PrintDigits(option);
      if(strstr(option,"table")) gObjectTable->Print();
      
      //increment the total number of Digits per run 
      fDigitsInRun += emcalLoader->Digits()->GetEntriesFast() ;  
    }//loop
       
  }//loader exists
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALDigitizer");
    AliInfo(Form("Digitize: took %f seconds for Digitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALDigitizer"), gBenchmark->GetCpuTime("EMCALDigitizer")/nEvents )) ;
  } 
}

//__________________________________________________________________
Float_t AliEMCALDigitizer::GetTimeResolution(Float_t energy) const
{  
  // From F. Blanco
  Float_t res = -1;
  if (energy > 0) {
    res = TMath::Sqrt(fTimeResolutionPar0 + 
		      fTimeResolutionPar1/(energy*energy) );
  }
  // parametrization above is for ns. Convert to seconds:
  return res*1e-9;
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::Digits2FastOR(TClonesArray* digitsTMP, TClonesArray* digitsTRG)
{
  // FEE digits afterburner to produce TRG digits 
  // we are only interested in the FEE digit deposited energy
  // to be converted later into a voltage value
  
  // push the FEE digit to its associated FastOR (numbered from 0:95)
  // TRU is in charge of summing module digits
  
  AliRunLoader *runLoader = AliRunLoader::Instance();
  
  AliRun* run = runLoader->GetAliRun();
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(runLoader->GetDetectorLoader("EMCAL"));
  if(!emcalLoader){
    AliFatal("Did not get the  Loader");
  }
  else {
    AliEMCAL* emcal = dynamic_cast<AliEMCAL*>(run->GetDetector("EMCAL"));
    if(emcal){
      AliEMCALGeometry* geom = emcal->GetGeometry();
      
      // build FOR from simulated digits
      // and xfer to the corresponding TRU input (mapping)
      
      TClonesArray* digits = emcalLoader->Digits();
      
      TIter NextDigit(digits);
      while (AliEMCALDigit* digit = (AliEMCALDigit*)NextDigit())
      {
        Int_t id = digit->GetId();
        
        Int_t trgid;
        if (geom && geom->GetFastORIndexFromCellIndex(id , trgid)) 
        {
          AliDebug(1,Form("trigger digit id: %d from cell id: %d\n",trgid,id));
          
          AliEMCALDigit* d = static_cast<AliEMCALDigit*>(digitsTMP->At(trgid));
          
          if (!d)
          {
            new((*digitsTMP)[trgid]) AliEMCALDigit(*digit);
            d = (AliEMCALDigit*)digitsTMP->At(trgid);
            d->SetId(trgid);
          }	
          else
          {
            *d = *d + *digit;
          }
        }
      }
      
      if (AliDebugLevel()) printf("Number of TRG digits: %d\n",digitsTMP->GetEntriesFast());
      
      Int_t    nSamples = 32;
      Int_t *timeSamples = new Int_t[nSamples];
      
      NextDigit = TIter(digitsTMP);
      while (AliEMCALDigit* digit = (AliEMCALDigit*)NextDigit())
      {
        if (digit)
        {
          Int_t     id = digit->GetId();
          Float_t time = digit->GetTime();
          
          Double_t depositedEnergy = 0.;
          for (Int_t j = 1; j <= digit->GetNprimary(); j++) depositedEnergy += digit->GetDEPrimary(j);
          
          if (AliDebugLevel()) printf("Deposited Energy: %f\n", depositedEnergy);
          
          // FIXME: Check digit time!
          if (depositedEnergy)
          {
            DigitalFastOR(time, depositedEnergy, timeSamples, nSamples);
            
            for (Int_t j=0;j<nSamples;j++) 
            {
              timeSamples[j] = ((j << 12) & 0xFF000) | (timeSamples[j] & 0xFFF);
            }
            
            new((*digitsTRG)[digitsTRG->GetEntriesFast()]) AliEMCALRawDigit(id, timeSamples, nSamples);
          }
        }
      }
      
      delete [] timeSamples;
    }// AliEMCAL exists
    else AliFatal("Could not get AliEMCAL");
  }// loader exists
  
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::DigitalFastOR( Double_t time, Double_t dE, Int_t timeSamples[], Int_t nSamples )
{
	// parameters:	
	// id: 0..95
	const Int_t    reso = 11;      // 11-bit resolution ADC
	const Double_t vFSR = 1;       // Full scale input voltage range
	const Double_t dNe   = 125;     // signal of the APD per MeV of energy deposit in a tower: 125 photo-e-/MeV @ M=30
	const Double_t vA   = .136e-6; // CSP output range: 0.136uV/e-
	const Double_t rise = 40e-9;   // rise time (10-90%) of the FastOR signal before shaping
	
	const Double_t kTimeBinWidth = 25E-9; // sampling frequency (40MHz)
	
	Double_t vV = 1000. * dE * dNe * vA; // GeV 2 MeV
	
	TF1 signalF("signal", AnalogFastORFunction, 0, nSamples * kTimeBinWidth, 3);
	signalF.SetParameter( 0,   vV ); 
	signalF.SetParameter( 1, time ); // FIXME: when does the signal arrive? Might account for cable lengths
	signalF.SetParameter( 2, rise );
	
	for (Int_t iTime=0; iTime<nSamples; iTime++) 
	{
		// FIXME: add noise (probably not simply Gaussian) according to DA measurements
		// probably plan an access to OCDB
		
		timeSamples[iTime] = int((TMath::Power(2, reso) / vFSR) * signalF.Eval(iTime * kTimeBinWidth) + 0.5);
	}
}

//____________________________________________________________________________ 
Bool_t AliEMCALDigitizer::Init()
{
  // Makes all memory allocations
  fInit = kTRUE ; 
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  
  if ( emcalLoader == 0 ) {
    Fatal("Init", "Could not obtain the AliEMCALLoader");  
    return kFALSE;
  } 
  
  fFirstEvent = 0 ; 
  fLastEvent = fFirstEvent ; 
  
  if(fDigInput)
    fInput = fDigInput->GetNinputs() ; 
  else 
    fInput           = 1 ;

  fInputFileNames  = new TString[fInput] ; 
  fEventNames      = new TString[fInput] ; 
  fInputFileNames[0] = GetTitle() ; 
  fEventNames[0]     = fEventFolderName.Data() ; 
  Int_t index ; 
  for (index = 1 ; index < fInput ; index++) {
    fInputFileNames[index] = dynamic_cast<AliStream*>(fDigInput->GetInputStream(index))->GetFileName(0); 
    TString tempo = fDigInput->GetInputFolderName(index) ;
    fEventNames[index] = tempo.Remove(tempo.Length()-1) ; // strip of the stream number added bt fDigInput 
  }
    
  //Calibration instance
  fCalibData = emcalLoader->CalibData();
  return fInit ;    
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::InitParameters()
{ 
  // Parameter initialization for digitizer
  
  // Get the parameters from the OCDB via the loader
  AliRunLoader *rl = AliRunLoader::Instance();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  AliEMCALSimParam * simParam = 0x0;
  if(emcalLoader) simParam = emcalLoader->SimulationParameters();
	
  if(!simParam){
    simParam = AliEMCALSimParam::GetInstance();
    AliWarning("Simulation Parameters not available in OCDB?");
  }
  
  fMeanPhotonElectron = simParam->GetMeanPhotonElectron() ; // 4400;  // electrons per GeV 
  fGainFluctuations   = simParam->GetGainFluctuations()   ; // 15.0; 
  
  fPinNoise           = simParam->GetPinNoise();//0.012; // pin noise in GeV from analysis test beam data 
  if (fPinNoise < 0.0001 ) 
    Warning("InitParameters", "No noise added\n") ; 
  fTimeNoise          = simParam->GetTimeNoise(); // 1.28E-5 s
  fDigitThreshold     = simParam->GetDigitThreshold(); //fPinNoise * 3; // 3 * sigma
  fTimeResolutionPar0 = simParam->GetTimeResolutionPar0(); 
  fTimeResolutionPar1 = simParam->GetTimeResolutionPar1(); 
  fTimeDelay          = simParam->GetTimeDelay(); //600e-9 ; // 600 ns
  
  // These defaults are normally not used. 
  // Values are read from calibration database instead
  fADCchannelEC       = 0.0153; // Update 24 Apr 2007: 250./16/1024 - width of one ADC channel in GeV
  fADCpedestalEC      = 0.0 ;   // GeV
  fADCchannelECDecal  = 1.0;    // No decalibration by default
  fTimeChannel        = 0.0;    // No time calibration by default
  fTimeChannelDecal   = 0.0;    // No time decalibration by default

  fNADCEC             = simParam->GetNADCEC();//(Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC - 65536
  
  AliDebug(2,Form("Mean Photon Electron %d, Gain Fluct. %2.1f; Noise: APD %f, Time %f; Digit Threshold %d,Time Resolution Par0 %g Par1 %g,NADCEC %d",
                  fMeanPhotonElectron, fGainFluctuations, fPinNoise,fTimeNoise, fDigitThreshold,fTimeResolutionPar0,fTimeResolutionPar1,fNADCEC));
  
}

//__________________________________________________________________
void AliEMCALDigitizer::Print1(Option_t * option)
{ // 19-nov-04 - just for convenience
  Print(); 
  PrintDigits(option);
}

//__________________________________________________________________
void AliEMCALDigitizer::Print (Option_t * ) const
{
  // Print Digitizer's parameters
  printf("Print: \n------------------- %s -------------", GetName() ) ; 
  if( strcmp(fEventFolderName.Data(), "") != 0 ){
    printf(" Writing Digits to branch with title  %s\n", fEventFolderName.Data()) ;
    
    Int_t nStreams ; 
    if (fDigInput) 
      nStreams =  GetNInputStreams() ;
    else 
      nStreams = fInput ; 
    
    AliRunLoader *rl=0;
    
    Int_t index = 0 ;  
    for (index = 0 ; index < nStreams ; index++) {  
      TString tempo(fEventNames[index]) ; 
      tempo += index ;
      if ((rl = AliRunLoader::GetRunLoader(tempo)) == 0)
        rl = AliRunLoader::Open(fInputFileNames[index], tempo) ; 
      AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
      if(emcalLoader){
        TString fileName( emcalLoader->GetSDigitsFileName() ) ; 
        if ( fEventNames[index] != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
          fileName = fileName.ReplaceAll(".root", "") + "_" + fEventNames[index]  + ".root" ;
        printf ("Adding SDigits from %s %s\n", fInputFileNames[index].Data(), fileName.Data()) ; 
      }//loader
    }
    
    AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
    
    if(emcalLoader) printf("\nWriting digits to %s", emcalLoader->GetDigitsFileName().Data()) ;
    else printf("\nNULL LOADER");
    
    printf("\nWith following parameters:\n") ;
    printf("    Electronics noise in EMC, APD (fPinNoise) = %f, Time = %f \n", fPinNoise, fTimeNoise) ;
    printf("    Threshold  in Tower  (fDigitThreshold) = %d\n", fDigitThreshold)  ;
    printf("---------------------------------------------------\n")  ;
  }
  else
    printf("Print: AliEMCALDigitizer not initialized") ; 
}

//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option)
{
  //utility method for printing digit information
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  if(emcalLoader){
    TClonesArray * digits  = emcalLoader->Digits() ;
    TClonesArray * sdigits = emcalLoader->SDigits() ;
    
    printf("\n #Digits: %d : sdigits %d ", digits->GetEntriesFast(), sdigits->GetEntriesFast()) ; 
    printf("\n event %d", emcalLoader->GetRunLoader()->GetEventNumber());
    
    if(strstr(option,"all")){  
      //loop over digits
      AliEMCALDigit * digit;
      printf("\nEMC digits (with primaries):\n")  ;
      printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
      Int_t index ;
      for (index = 0 ; index < digits->GetEntries() ; index++) {
        digit = dynamic_cast<AliEMCALDigit *>(digits->At(index)) ;
        if(digit){
          printf("\n%6d  %8f    %6.5e %4d      %2d : ",
                 digit->GetId(), digit->GetAmplitude(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
          Int_t iprimary;
          for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
            printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
          }
        }// digit exists
      }// loop
    }
    printf("\n");
  }// loader exists
  else printf("NULL LOADER, cannot print\n");
}

//__________________________________________________________________
Float_t AliEMCALDigitizer::TimeOfNoise(void)
{  
  // Calculates the time signal generated by noise
  //printf("Time noise %e\n",fTimeNoise);
  return gRandom->Rndm() * fTimeNoise;
}

//__________________________________________________________________
void AliEMCALDigitizer::Unload() 
{  
  // Unloads the SDigits and Digits
  AliRunLoader *rl=0;
    
  Int_t i ; 
  for(i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i;
    if ((rl = AliRunLoader::GetRunLoader(tempo))) 
      rl->GetDetectorLoader("EMCAL")->UnloadSDigits() ; 
  }
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  if(emcalLoader)emcalLoader->UnloadDigits() ; 
  else AliFatal("Did not get the loader");
}

//_________________________________________________________________________________________
void AliEMCALDigitizer::WriteDigits()
{ // Makes TreeD in the output file. 
  // Check if branch already exists: 
  //   if yes, exit without writing: ROOT TTree does not support overwriting/updating of 
  //      already existing branches. 
  //   else creates branch with Digits, named "EMCAL", title "...",
  //      and branch "AliEMCALDigitizer", with the same title to keep all the parameters
  //      and names of files, from which digits are made.
  
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  
  if(emcalLoader){
    const TClonesArray * digits = emcalLoader->Digits() ; 
    TTree * treeD = emcalLoader->TreeD(); 
    if ( !treeD ) {
      emcalLoader->MakeDigitsContainer();
      treeD = emcalLoader->TreeD(); 
    }
    
    // -- create Digits branch
    Int_t bufferSize = 32000 ;    
    TBranch * digitsBranch = 0;
    if ((digitsBranch = treeD->GetBranch("EMCAL"))) {
      digitsBranch->SetAddress(&digits);
      AliWarning("Digits Branch already exists. Not all digits will be visible");
    }
    else
      treeD->Branch("EMCAL","TClonesArray",&digits,bufferSize);
    //digitsBranch->SetTitle(fEventFolderName);
    
    //	treeD->Fill() ;
    /*  
     emcalLoader->WriteDigits("OVERWRITE");
     emcalLoader->WriteDigitizer("OVERWRITE");
     
     Unload() ; 
     */
    
  }// loader exists
  else AliFatal("Loader not available");
}

//__________________________________________________________________
void AliEMCALDigitizer::WriteDigits(TClonesArray* digits, const char* branchName)
{ // overloaded method
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  if(emcalLoader){
    
    TTree* treeD = emcalLoader->TreeD(); 
    if (!treeD) 
      {
	emcalLoader->MakeDigitsContainer();
	treeD = emcalLoader->TreeD(); 
      }
    
    // -- create Digits branch
    Int_t bufferSize = 32000;
    
    if (TBranch* triggerBranch = treeD->GetBranch(branchName)) 
      {
	triggerBranch->SetAddress(&digits);
      }
    else
      {
	treeD->Branch(branchName,"TClonesArray",&digits,bufferSize);
      }
    
    //	treeD->Fill();
  }// loader exists
  else AliFatal("Loader not available");
}

