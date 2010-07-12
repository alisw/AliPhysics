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
// In addition it performs mixing of summable digits from different events.
//
// For each event two branches are created in TreeD:
//   "EMCAL" - list of digits
//   "AliEMCALDigitizer" - AliEMCALDigitizer with all parameters used in digitization
//
// Note, that one cset title for new digits branch, and repeat digitization with
// another set of parameters.
//
// Examples of use:
// root[0] AliEMCALDigitizer * d = new AliEMCALDigitizer() ;
// root[1] d->ExecuteTask()             
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       //Digitizes SDigitis in all events found in file galice.root 
//
// root[2] AliEMCALDigitizer * d1 = new AliEMCALDigitizer("galice1.root") ;  
//                       // Will read sdigits from galice1.root
// root[3] d1->MixWith("galice2.root")       
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       // Reads another portion of sdigits from galice2.root
// root[3] d1->MixWith("galice3.root")       
//                       // Reads another portion of sdigits from galice3.root
// root[4] d->ExecuteTask("deb timing")    
//                       // Reads SDigits from files galice1.root, galice2.root ....
//                       // mixes them and stores produced Digits in file galice1.root          
//                       // deb - prints number of produced digits
//                       // deb all - prints list of produced digits
//                       // timing  - prints time used for digitization
////////////////////////////////////////////////////////////////////////////////////
//
//*-- Author: Sahal Yacoob (LBL)
// based on : AliEMCALDigitizer
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//  November 2003 Aleksei Pavlinov : adopted for Shish-Kebab geometry 
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
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliEMCALDigit.h"
#include "AliEMCAL.h"
#include "AliEMCALLoader.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
//#include "AliEMCALTick.h"
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
//    fPedestal(0), //Not used, remove?
//    fSlope(0),    //Not used, remove?
    fPinNoise(0),
    fTimeDelay(0),
    fTimeResolution(0),
//    fTimeThreshold(0),    //Not used, remove?
//    fTimeSignalLength(0), //Not used, remove?
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(""),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0)
{
  // ctor
  InitParameters() ; 
  fManager = 0 ;                     // We work in the standalone mode
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(TString alirunFileName, TString eventFolderName)
  : AliDigitizer("EMCAL"+AliConfig::Instance()->GetDigitizerTaskName(), alirunFileName),
    fDefaultInit(kFALSE),
    fDigitsInRun(0),
    fInit(0),
    fInput(0),
    fInputFileNames(0), 
    fEventNames(0), 
    fDigitThreshold(0),
    fMeanPhotonElectron(0),
//    fPedestal(0),//Not used, remove?
//    fSlope(0),   //Not used, remove?
    fPinNoise(0),
	fTimeDelay(0),
    fTimeResolution(0),
//    fTimeThreshold(0),    //Not used, remove?
//    fTimeSignalLength(0), //Not used, remove?
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(eventFolderName),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fManager = 0 ;                     // We work in the standalone mode
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
//    fPedestal(d.fPedestal), //Not used, remove?
//    fSlope(d.fSlope),       //Not used, remove?
    fPinNoise(d.fPinNoise),
    fTimeDelay(d.fTimeDelay),
    fTimeResolution(d.fTimeResolution),
//    fTimeThreshold(d.fTimeThreshold),       //Not used, remove?
//    fTimeSignalLength(d.fTimeSignalLength), //Not used, remove?
    fADCchannelEC(d.fADCchannelEC),
    fADCpedestalEC(d.fADCpedestalEC),
    fNADCEC(d.fNADCEC),
    fEventFolderName(d.fEventFolderName),
    fFirstEvent(d.fFirstEvent),
    fLastEvent(d.fLastEvent),
    fCalibData(d.fCalibData)
{
  // copyy ctor 
 }

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(AliRunDigitizer * rd)
  : AliDigitizer(rd,"EMCAL"+AliConfig::Instance()->GetDigitizerTaskName()),
    fDefaultInit(kFALSE),
    fDigitsInRun(0),
    fInit(0),
    fInput(0),
    fInputFileNames(0),
    fEventNames(0),
    fDigitThreshold(0),
    fMeanPhotonElectron(0),
//    fPedestal(0), //Not used, remove?
//    fSlope(0.),   //Not used, remove?
    fPinNoise(0.),
    fTimeDelay(0.),
    fTimeResolution(0.),
//    fTimeThreshold(0),    //Not used, remove?
//    fTimeSignalLength(0), //Not used, remove?
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(0),
    fFirstEvent(0),
    fLastEvent(0),
    fCalibData(0x0)
{
  // ctor Init() is called by RunDigitizer
  fManager = rd ; 
  fEventFolderName = fManager->GetInputFolderName(0) ;
  SetTitle(dynamic_cast<AliStream*>(fManager->GetInputStream(0))->GetFileName(0));
  InitParameters() ; 
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  //dtor
  if (AliRunLoader::Instance()) {
    AliLoader *emcalLoader=0;
    if ((emcalLoader = AliRunLoader::Instance()->GetDetectorLoader("EMCAL")))
      emcalLoader->CleanDigitizer();
  }
  else
    AliDebug(1," no runloader present");
  delete [] fInputFileNames ; 
  delete [] fEventNames ; 

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
  Int_t readEvent = event ; 
  // fManager is data member from AliDigitizer
  if (fManager) 
    readEvent = dynamic_cast<AliStream*>(fManager->GetInputStream(0))->GetCurrentEventNumber() ; 
  AliDebug(1,Form("Adding event %d from input stream 0 %s %s", 
		  readEvent, GetTitle(), fEventFolderName.Data())) ; 
  rl->GetEvent(readEvent);

  TClonesArray * digits = emcalLoader->Digits() ; 
  digits->Delete() ;  //correct way to clear array when memory is
		      //allocated by objects stored in array

  // Load Geometry
  AliEMCALGeometry *geom = 0;
  if (rl->GetAliRun()) {
    AliEMCAL * emcal  = (AliEMCAL*)rl->GetAliRun()->GetDetector("EMCAL");
    geom = emcal->GetGeometry();
  }
  else 
    AliFatal("Could not get AliRun from runLoader");

  nEMC = geom->GetNCells();
  AliDebug(1,Form("nEMC %i (number cells in EMCAL) | %s \n", nEMC, geom->GetName()));
  
  Int_t absID ;

  digits->Expand(nEMC) ;

  // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  if ( !emcalLoader->SDigitizer() ) 
    emcalLoader->LoadSDigitizer();
  AliEMCALSDigitizer * sDigitizer = dynamic_cast<AliEMCALSDigitizer*>(emcalLoader->SDigitizer()); 
  
  if ( !sDigitizer )
    Fatal("Digitize", "SDigitizer with name %s %s not found", fEventFolderName.Data(), GetTitle() ) ; 

  //take all the inputs to add together and load the SDigits
  TObjArray * sdigArray = new TObjArray(fInput) ;
  sdigArray->AddAt(emcalLoader->SDigits(), 0) ;
  Int_t i ;

  for(i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i ;

    AliRunLoader *  rl2 = AliRunLoader::GetRunLoader(tempo) ; 

    if (rl2==0) 
      rl2 = AliRunLoader::Open(fInputFileNames[i], tempo) ; 

    if (fManager) 
      readEvent = dynamic_cast<AliStream*>(fManager->GetInputStream(i))->GetCurrentEventNumber() ; 
    Info("Digitize", "Adding event %d from input stream %d %s %s", readEvent, i, fInputFileNames[i].Data(), tempo.Data()) ; 
    rl2->LoadSDigits();
    rl2->GetEvent(readEvent);
    AliEMCALLoader *emcalLoader2 = dynamic_cast<AliEMCALLoader*>(rl2->GetDetectorLoader("EMCAL"));
    sdigArray->AddAt(emcalLoader2->SDigits(), i) ;
  }
  
  //Find the first tower with signal
  Int_t nextSig = nEMC + 1 ; 
  TClonesArray * sdigits ;  
  for(i = 0 ; i < fInput ; i++){
    sdigits = dynamic_cast<TClonesArray *>(sdigArray->At(i)) ;
    if ( !sdigits->GetEntriesFast() )
      continue ; 
    Int_t curNext = dynamic_cast<AliEMCALDigit *>(sdigits->At(0))->GetId() ;
     if(curNext < nextSig) 
       nextSig = curNext ;
     AliDebug(1,Form("input %i : #sdigits %i \n",
		     i, sdigits->GetEntriesFast()));
  }
  AliDebug(1,Form("FIRST tower with signal %i \n", nextSig));

  TArrayI index(fInput) ;
  index.Reset() ;  //Set all indexes to zero

  AliEMCALDigit * digit ;
  AliEMCALDigit * curSDigit ;

  //  TClonesArray * ticks = new TClonesArray("AliEMCALTick",1000) ;

  //Put Noise contribution
  for(absID = 0; absID < nEMC; absID++){ // Nov 30, 2006 by PAI; was from 1 to nEMC
    Float_t energy = 0 ;
    // amplitude set to zero, noise will be added later
    new((*digits)[absID]) AliEMCALDigit( -1, -1, absID, 0., TimeOfNoise(),kFALSE); // absID-1->absID
    //look if we have to add signal?
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(absID)); // absID-1->absID
    
    if(absID==nextSig){
      //Add SDigits from all inputs    
      //      ticks->Clear() ;
      //Int_t contrib = 0 ;

      //Follow PHOS and comment out this timing model til a better one
      //can be developed - JLK 28-Apr-2008

      //Float_t a = digit->GetAmplitude() ;
      //Float_t b = TMath::Abs( a /fTimeSignalLength) ;
      //Mark the beginning of the signal
      //new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime(),0, b);  
      //Mark the end of the signal     
      //new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime()+fTimeSignalLength, -a, -b);

      // Calculate time as time of the largest digit
      Float_t time = digit->GetTime() ;
      Float_t aTime= digit->GetAmplitude() ;
      
      // loop over input
      for(i = 0; i< fInput ; i++){  //loop over (possible) merge sources
    	if(dynamic_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	  curSDigit = dynamic_cast<AliEMCALDigit*>(dynamic_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ; 	
	else
	  curSDigit = 0 ;
	//May be several digits will contribute from the same input
	while(curSDigit && (curSDigit->GetId() == absID)){	   
	  //Shift primary to separate primaries belonging different inputs
	  Int_t primaryoffset ;
	  if(fManager)
	    primaryoffset = fManager->GetMask(i) ; 
	  else
	    primaryoffset = i ;
	  curSDigit->ShiftPrimary(primaryoffset) ;

	  //Remove old timing model - JLK 28-April-2008
	  //a = curSDigit->GetAmplitude() ;
	  //b = a /fTimeSignalLength ;
	  //new((*ticks)[contrib++]) AliEMCALTick(curSDigit->GetTime(),0, b);  
	  //new((*ticks)[contrib++]) AliEMCALTick(curSDigit->GetTime()+fTimeSignalLength, -a, -b); 
	  if(curSDigit->GetAmplitude()>aTime) {
	    aTime = curSDigit->GetAmplitude();
	    time = curSDigit->GetTime();
	  }

	  *digit = *digit + *curSDigit ;  //adds amplitudes of each digit

	  index[i]++ ;
	  if( dynamic_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	    curSDigit = dynamic_cast<AliEMCALDigit*>(dynamic_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ;
	  else
	    curSDigit = 0 ;
	}
      }
      //Here we convert the summed amplitude to an energy in GeV
      energy = sDigitizer->Calibrate(digit->GetAmplitude()) ; // GeV
      // add fluctuations for photo-electron creation
      energy *= static_cast<Float_t>(gRandom->Poisson(fMeanPhotonElectron)) / static_cast<Float_t>(fMeanPhotonElectron) ;
  
      //calculate and set time
      //New timing model needed - JLK 28-April-2008
      //Float_t time = FrontEdgeTime(ticks) ;
      digit->SetTime(time) ;

      //Find next signal module
      nextSig = nEMC + 1 ;
      for(i = 0 ; i < fInput ; i++){
	sdigits = dynamic_cast<TClonesArray *>(sdigArray->At(i)) ;
	Int_t curNext = nextSig ;
	if(sdigits->GetEntriesFast() > index[i] ){
	  curNext = dynamic_cast<AliEMCALDigit *>(sdigits->At(index[i]))->GetId() ;	  
	}
	if(curNext < nextSig) nextSig = curNext ;
      }
    }
    // add the noise now
    energy += TMath::Abs(gRandom->Gaus(0., fPinNoise)) ;
    // JLK 26-June-2008
    //Now digitize the energy according to the sDigitizer method,
    //which merely converts the energy to an integer.  Later we will
    //check that the stored value matches our allowed dynamic ranges
    digit->SetAmplitude(sDigitizer->Digitize(energy)) ;  
    AliDebug(10,Form(" absID %5i energy %f nextSig %5i\n",
		     absID, energy, nextSig));
  } // for(absID = 0; absID < nEMC; absID++)
  
  //ticks->Delete() ;
  //delete ticks ;

  //JLK is it better to call Clear() here?
  delete sdigArray ; //We should not delete its contents

  //remove digits below thresholds
  // until 10-02-2010 remove digits with energy smaller than fDigitThreshold 3*fPinNoise
  // now, remove digits with Digitized ADC smaller than fDigitThreshold = 3
  Float_t energy=0;
  for(i = 0 ; i < nEMC ; i++){
    digit = dynamic_cast<AliEMCALDigit*>( digits->At(i) ) ;
    //First get the energy in GeV units.
    energy = sDigitizer->Calibrate(digit->GetAmplitude()) ;
    //Then digitize using the calibration constants of the ocdb
    Float_t ampADC = DigitizeEnergy(energy, digit->GetId())  ; 	  
    //if(ampADC>2)printf("Digit energy %f, amp %d, amp cal %d, threshold %d\n",energy,digit->GetAmplitude(),ampADC,fDigitThreshold);
    if(ampADC < fDigitThreshold)
      digits->RemoveAt(i) ;
    else 
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
  }
  
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ; 
  
  //JLK 26-June-2008
  //After we have done the summing and digitizing to create the
  //digits, now we want to calibrate the resulting amplitude to match
  //the dynamic range of our real data.  
  for (i = 0 ; i < ndigits ; i++) { 
    digit = dynamic_cast<AliEMCALDigit *>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ; 
    energy = sDigitizer->Calibrate(digit->GetAmplitude()) ;
    digit->SetAmplitude(DigitizeEnergy(energy, digit->GetId()) ) ;
	//Add delay to time
	digit->SetTime(digit->GetTime()+fTimeDelay) ;
	 // printf("digit amplitude set at end: i %d, amp %f\n",i,digit->GetAmplitude());
  }

}

// //_____________________________________________________________________
Float_t AliEMCALDigitizer::DigitizeEnergy(Float_t energy, Int_t AbsId)
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

  if (geom==0)
    AliFatal("Did not get geometry from EMCALLoader");

  Int_t iSupMod = -1;
  Int_t nModule  = -1;
  Int_t nIphi   = -1;
  Int_t nIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Float_t channel = -999; 

  Bool_t bCell = geom->GetCellIndex(AbsId, iSupMod, nModule, nIphi, nIeta) ;
  if(!bCell)
    Error("DigitizeEnergy","Wrong cell id number : AbsId %i ", AbsId) ;
  geom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);
  
  if(fCalibData) {
    fADCpedestalEC = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);
    fADCchannelEC = fCalibData->GetADCchannel(iSupMod,ieta,iphi);
  }
  
  //channel = static_cast<Int_t> (TMath::Floor( (energy + fADCpedestalEC)/fADCchannelEC ))  ;
  channel = (energy + fADCpedestalEC)/fADCchannelEC   ;

  if(channel > fNADCEC ) 
    channel =  fNADCEC ; 
  return channel ;
  
}

//____________________________________________________________________________
void AliEMCALDigitizer::Exec(Option_t *option) 
{ 
  // Steering method to process digitization for events
  // in the range from fFirstEvent to fLastEvent.
  // This range is optionally set by SetEventRange().
  // if fLastEvent=-1, then process events until the end.
  // by default fLastEvent = fFirstEvent (process only one event)

  if (!fInit) { // to prevent overwrite existing file
    Error( "Exec", "Give a version name different from %s", fEventFolderName.Data() ) ;
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

  // Post Digitizer to the white board
  emcalLoader->PostDigitizer(this) ;
  
  if (fLastEvent == -1)  {
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  }
  else if (fManager) 
    fLastEvent = fFirstEvent ; // what is this ??

  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  Int_t ievent;

  TClonesArray* digitsTRG = new TClonesArray("AliEMCALRawDigit", 32 * 96);
  TClonesArray* digitsTMP = new TClonesArray("AliEMCALDigit",    32 * 96);
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
	emcalLoader->WriteDigitizer("OVERWRITE");
	  
	Unload();
	  
	digitsTRG->Clear();
	digitsTMP->Clear();
	//-------------------------------------

    if(strstr(option,"deb"))
      PrintDigits(option);
    if(strstr(option,"table")) gObjectTable->Print();

    //increment the total number of Digits per run 
    fDigitsInRun += emcalLoader->Digits()->GetEntriesFast() ;  
  }
  
  emcalLoader->CleanDigitizer() ;

  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALDigitizer");
    AliInfo(Form("Exec: took %f seconds for Digitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALDigitizer"), gBenchmark->GetCpuTime("EMCALDigitizer")/nEvents )) ;
  } 
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
	
	AliEMCALGeometry* geom = dynamic_cast<AliEMCAL*>(run->GetDetector("EMCAL"))->GetGeometry();
	
	// build FOR from simulated digits
	// and xfer to the corresponding TRU input (mapping)
	
	TClonesArray* digits = emcalLoader->Digits();
	
	TIter NextDigit(digits);
	while (AliEMCALDigit* digit = (AliEMCALDigit*)NextDigit())
	{
		Int_t id = digit->GetId();
		
		Int_t iSupMod, nModule, nIphi, nIeta, iphi, ieta, iphim, ietam;
		
		geom->GetCellIndex(              id, iSupMod, nModule, nIphi, nIeta );
		geom->GetModulePhiEtaIndexInSModule( iSupMod, nModule, iphim, ietam );		
		geom->GetCellPhiEtaIndexInSModule(   iSupMod, nModule, nIphi, nIeta, iphi, ieta); 
		
		// identify to which TRU this FEE digit belong
		//Int_t itru = (iSupMod < 11) ? iphim / 4 + 3 * iSupMod : 31;
		Int_t itru = -1;
		if (iSupMod < 11)
			itru = (iSupMod % 2) ? (2 - int(iphim / 4)) + 3 * iSupMod : iphim / 4 + 3 * iSupMod;
		else 
			itru = 31;
		
		//---------
		//
		// FIXME: bad numbering solution to deal w/ the last 2 SM which have only 1 TRU each
		// using the AliEMCALGeometry official numbering
		// only 1 TRU/SM in SM 10 & SM 11
		//
		//---------
		if ((itru == 31 && iphim < 2) || (itru == 30 && iphim > 5)) continue;
		
		// to be compliant with %4 per TRU
		if (itru == 31) iphim -= 2;
		
		Int_t trgid;
		Bool_t isOK = geom->GetAbsFastORIndexFromPositionInTRU(itru, ietam, iphim % 4, trgid);
		
		AliDebug(2,Form("trigger digit id: %d itru: %d isOK: %d\n",trgid,itru,isOK));
		
		if (isOK) 
		{
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
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::DigitalFastOR( Double_t time, Double_t dE, Int_t timeSamples[], Int_t nSamples )
{
	// parameters:	
	// id: 0..95
	const Int_t    reso = 11;      // 11-bit resolution ADC
	const Double_t vFSR = 1;       // Full scale input voltage range
	const Double_t Ne   = 125;     // signal of the APD per MeV of energy deposit in a tower: 125 photo-e-/MeV @ M=30
	const Double_t vA   = .136e-6; // CSP output range: 0.136uV/e-
	const Double_t rise = 40e-9;   // rise time (10-90%) of the FastOR signal before shaping
	
	const Double_t kTimeBinWidth = 25E-9; // sampling frequency (40MHz)
	
	Double_t vV = 1000. * dE * Ne * vA; // GeV 2 MeV
	
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
//Float_t AliEMCALDigitizer::FrontEdgeTime(TClonesArray * ticks) 
//{ 
//  //  Returns the shortest time among all time ticks
//
//  ticks->Sort() ; //Sort in accordance with times of ticks
//  TIter it(ticks) ;
//  AliEMCALTick * ctick = (AliEMCALTick *) it.Next() ;
//  Float_t time = ctick->CrossingTime(fTimeThreshold) ;    
//  
//  AliEMCALTick * t ;  
//  while((t=(AliEMCALTick*) it.Next())){
//    if(t->GetTime() < time)  //This tick starts before crossing
//      *ctick+=*t ;
//    else
//      return time ;
//    
//    time = ctick->CrossingTime(fTimeThreshold) ;    
//  }
//  return time ;
//}
//

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
  
  if(fManager)
    fInput = fManager->GetNinputs() ; 
  else 
    fInput           = 1 ;

  fInputFileNames  = new TString[fInput] ; 
  fEventNames      = new TString[fInput] ; 
  fInputFileNames[0] = GetTitle() ; 
  fEventNames[0]     = fEventFolderName.Data() ; 
  Int_t index ; 
  for (index = 1 ; index < fInput ; index++) {
    fInputFileNames[index] = dynamic_cast<AliStream*>(fManager->GetInputStream(index))->GetFileName(0); 
    TString tempo = fManager->GetInputFolderName(index) ;
    fEventNames[index] = tempo.Remove(tempo.Length()-1) ; // strip of the stream number added bt fManager 
  }
  
  //to prevent cleaning of this object while GetEvent is called
  emcalLoader->GetDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);

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
	
  fMeanPhotonElectron = simParam->GetMeanPhotonElectron();//4400;  // electrons per GeV 
  fPinNoise           = simParam->GetPinNoise();//0.012; // pin noise in GeV from analysis test beam data 
  if (fPinNoise < 0.0001 ) 
    Warning("InitParameters", "No noise added\n") ; 
  fDigitThreshold     = simParam->GetDigitThreshold(); //fPinNoise * 3; // 3 * sigma
  fTimeResolution     = simParam->GetTimeResolution(); //0.6e-9 ; // 600 pc
  fTimeDelay          = simParam->GetTimeDelay(); //600e-9 ; // 600 ns

  // These defaults are normally not used. 
  // Values are read from calibration database instead
  fADCchannelEC       = 0.0153; // Update 24 Apr 2007: 250./16/1024 - width of one ADC channel in GeV
  fADCpedestalEC      = 0.0 ;  // GeV

  fNADCEC          = simParam->GetNADCEC();//(Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC - 65536

  AliDebug(2,Form("Mean Photon Electron %d, Noise %f, Digit Threshold %d,Time Resolution %g,NADCEC %d",
		fMeanPhotonElectron,fPinNoise,fDigitThreshold,fTimeResolution,fNADCEC));

  // Not used anymore, remove?
  // fTimeSignalLength   = 1.0e-9 ;
  // fTimeThreshold      = 0.001*10000000 ; // Means 1 MeV in terms of SDigits amplitude ??

}

//__________________________________________________________________
void AliEMCALDigitizer::MixWith(TString alirunFileName, TString eventFolderName)
{
  // Allows to produce digits by superimposing background and signal event.
  // It is assumed, that headers file with SIGNAL events is opened in 
  // the constructor. 
  // Sets the BACKGROUND event, with which the SIGNAL event is to be mixed 
  // Thus we avoid writing (changing) huge and expensive 
  // backgound files: all output will be writen into SIGNAL, i.e. 
  // opened in constructor file. 
  //
  // One can open as many files to mix with as one needs.
  // However only Sdigits with the same name (i.e. constructed with the same SDigitizer)
  // can be mixed.

  if( strcmp(GetName(), "") == 0 )
    Init() ;
  
  if(fManager){
    Error("MixWith", "Cannot use this method under AliRunDigitizer") ;
    return ;
  } 
  // looking for file which contains AliRun
  if (gSystem->AccessPathName(alirunFileName)) {// file does not exist
    Error("MixWith", "File %s does not exist!", alirunFileName.Data()) ;
    return ; 
  }
  // looking for the file which contains SDigits
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
  TString fileName( emcalLoader->GetSDigitsFileName() ) ; 
    if ( eventFolderName != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
      fileName = fileName.ReplaceAll(".root", "") + "_" + eventFolderName + ".root" ;
    if ( (gSystem->AccessPathName(fileName)) ) { 
      Error("MixWith", "The file %s does not exist!", fileName.Data()) ;
      return ;
    }
    // need to increase the arrays
    // MvL: This code only works when fInput == 1, I think.
    TString tempo = fInputFileNames[fInput-1] ; 
    delete [] fInputFileNames ; 
    fInputFileNames = new TString[fInput+1] ; 
    fInputFileNames[fInput-1] = tempo ; 
 
    tempo = fEventNames[fInput-1] ; 
    delete [] fEventNames ; 
    fEventNames = new TString[fInput+1] ; 
    fEventNames[fInput-1] = tempo ; 

    fInputFileNames[fInput] = alirunFileName ; 
    fEventNames[fInput]     = eventFolderName ;
    fInput++ ;
}  

//__________________________________________________________________
void AliEMCALDigitizer::Print1(Option_t * option)
{ // 19-nov-04 - just for convinience
  Print(); 
  PrintDigits(option);
}

//__________________________________________________________________
void AliEMCALDigitizer::Print(Option_t*)const 
{
  // Print Digitizer's parameters
  printf("Print: \n------------------- %s -------------", GetName() ) ; 
  if( strcmp(fEventFolderName.Data(), "") != 0 ){
    printf(" Writing Digits to branch with title  %s\n", fEventFolderName.Data()) ;
    
    Int_t nStreams ; 
    if (fManager) 
      nStreams =  GetNInputStreams() ;
    else 
      nStreams = fInput ; 
    
    Int_t index = 0 ;  

    AliRunLoader *rl=0;

    for (index = 0 ; index < nStreams ; index++) {  
      TString tempo(fEventNames[index]) ; 
      tempo += index ;
      if ((rl = AliRunLoader::GetRunLoader(tempo)) == 0)
	rl = AliRunLoader::Open(fInputFileNames[index], tempo) ; 
      AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
      TString fileName( emcalLoader->GetSDigitsFileName() ) ; 
      if ( fEventNames[index] != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
	fileName = fileName.ReplaceAll(".root", "") + "_" + fEventNames[index]  + ".root" ;
      printf ("Adding SDigits from %s %s\n", fInputFileNames[index].Data(), fileName.Data()) ; 
    }

    AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));

    printf("\nWriting digits to %s", emcalLoader->GetDigitsFileName().Data()) ;
    
    printf("\nWith following parameters:\n") ;
    
    printf("    Electronics noise in EMC (fPinNoise) = %f\n", fPinNoise) ;
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
      printf("\n%6d  %8f    %6.5e %4d      %2d : ",
	      digit->GetId(), digit->GetAmplitude(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	printf("%d ",digit->GetPrimary(iprimary+1) ) ; 
      }
    }   
  }
  printf("\n");
}

//__________________________________________________________________
Float_t AliEMCALDigitizer::TimeOfNoise(void)
{  
  // Calculates the time signal generated by noise
  //PH  Info("TimeOfNoise", "Change me") ; 
  return gRandom->Rndm() * 1.28E-5;
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
  emcalLoader->UnloadDigits() ; 
}

//_________________________________________________________________________________________
void AliEMCALDigitizer::WriteDigits()
{

  // Makes TreeD in the output file. 
  // Check if branch already exists: 
  //   if yes, exit without writing: ROOT TTree does not support overwriting/updating of 
  //      already existing branches. 
  //   else creates branch with Digits, named "EMCAL", title "...",
  //      and branch "AliEMCALDigitizer", with the same title to keep all the parameters
  //      and names of files, from which digits are made.

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));

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
}

//__________________________________________________________________
void AliEMCALDigitizer::WriteDigits(TClonesArray* digits, const char* branchName)
{
	//
	AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::Instance()->GetDetectorLoader("EMCAL"));
	
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
}

//__________________________________________________________________
void AliEMCALDigitizer::Browse(TBrowser* b)
{
  TTask::Browse(b);
}
