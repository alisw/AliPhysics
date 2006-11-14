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
// Class performs digitization of Summable digits 
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
#include <TList.h>
#include <TH1.h>
#include <TBrowser.h>
#include <TObjectTable.h>
#include <TRandom.h>

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRun.h"
#include "AliRunDigitizer.h"
#include "AliRunLoader.h"
#include "AliEMCALDigit.h"
#include "AliEMCAL.h"
#include "AliEMCALLoader.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTick.h"
#include "AliEMCALHistoUtilities.h"

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
    fPedestal(0),
    fSlope(0),
    fPinNoise(0),
    fTimeResolution(0),
    fTimeThreshold(0),    
    fTimeSignalLength(0),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(""),
    fFirstEvent(0),
    fLastEvent(0),
    fControlHists(0),
    fHists(0)
{
  // ctor
  InitParameters() ; 
  fManager = 0 ;                     // We work in the standalong mode
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
    fPedestal(0),
    fSlope(0),
    fPinNoise(0),
    fTimeResolution(0),
    fTimeThreshold(0),
    fTimeSignalLength(0),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(eventFolderName),
    fFirstEvent(0),
    fLastEvent(0),
    fControlHists(0),
    fHists(0)
{
  // ctor
  InitParameters() ; 
  Init() ;
  fManager = 0 ;                     // We work in the standalong mode
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
    fPedestal(d.fPedestal),
    fSlope(d.fSlope),
    fPinNoise(d.fPinNoise),
    fTimeResolution(d.fTimeResolution),
    fTimeThreshold(d.fTimeThreshold),
    fTimeSignalLength(d.fTimeSignalLength),
    fADCchannelEC(d.fADCchannelEC),
    fADCpedestalEC(d.fADCpedestalEC),
    fNADCEC(d.fNADCEC),
    fEventFolderName(d.fEventFolderName),
    fFirstEvent(d.fFirstEvent),
    fLastEvent(d.fLastEvent),
    fControlHists(d.fControlHists),
    fHists(d.fHists)
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
    fDigitThreshold(0.),
    fMeanPhotonElectron(0),
    fPedestal(0),
    fSlope(0.),
    fPinNoise(0),
    fTimeResolution(0.),
    fTimeThreshold(0),
    fTimeSignalLength(0),
    fADCchannelEC(0),
    fADCpedestalEC(0),
    fNADCEC(0),
    fEventFolderName(0),
    fFirstEvent(0),
    fLastEvent(0),
    fControlHists(0),
    fHists(0)
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
  if (AliRunLoader::GetRunLoader()) {
    AliLoader *emcalLoader=0;
    if ((emcalLoader = AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL")))
      emcalLoader->CleanDigitizer();
  }
  else
    AliDebug(1," no runloader present");
  delete [] fInputFileNames ; 
  delete [] fEventNames ; 

  if(fHists) delete fHists;
}

//____________________________________________________________________________
void AliEMCALDigitizer::Digitize(Int_t event) 
{ 

  // Makes the digitization of the collected summable digits
  // for this it first creates the array of all EMCAL modules
  // filled with noise (different for EMC, CPV and PPSD) and
  // after that adds contributions from SDigits. This design 
  // helps to avoid scanning over the list of digits to add 
  // contribution of any new SDigit.
  static int isTrd1Geom = -1; // -1 - mean undefined 
  static int nEMC=0; //max number of digits possible

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));
  Int_t readEvent = event ; 
  // fManager is data member from AliDigitizer
  if (fManager) 
    readEvent = dynamic_cast<AliStream*>(fManager->GetInputStream(0))->GetCurrentEventNumber() ; 
  AliDebug(1,Form("Adding event %d from input stream 0 %s %s", 
		  readEvent, GetTitle(), fEventFolderName.Data())) ; 
  rl->GetEvent(readEvent);

  TClonesArray * digits = emcalLoader->Digits() ; 
  digits->Clear() ;

  // Load Geometry
  // const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  rl->LoadgAlice(); 
  AliRun * gAlice = rl->GetAliRun(); 
  AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
  AliEMCALGeometry * geom = emcal->GetGeometry();

  if(isTrd1Geom < 0) { 
    TString ng(geom->GetName());
    isTrd1Geom = 0;
    if(ng.Contains("SHISH") &&  ng.Contains("TRD1")) isTrd1Geom = 1;

    if(isTrd1Geom == 0) nEMC = geom->GetNPhi()*geom->GetNZ();
    else                nEMC = geom->GetNCells();
    AliDebug(1,Form("nEMC %i (number cells in EMCAL) | %s | isTrd1Geom %i\n", nEMC, geom->GetName(), isTrd1Geom));
  }
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

  TClonesArray * ticks = new TClonesArray("AliEMCALTick",1000) ;

  //Put Noise contribution
  for(absID = 1; absID <= nEMC; absID++){
    Float_t amp = 0 ;
    // amplitude set to zero, noise will be added later
    new((*digits)[absID-1]) AliEMCALDigit( -1, -1, absID, 0, TimeOfNoise() ) ;
    //look if we have to add signal?
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(absID-1)) ;
    
    if(absID==nextSig){
      //Add SDigits from all inputs    
      ticks->Clear() ;
      Int_t contrib = 0 ;
      Float_t a = digit->GetAmp() ;
      Float_t b = TMath::Abs( a /fTimeSignalLength) ;
      //Mark the beginning of the signal
      new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime(),0, b);  
      //Mark the end of the signal     
      new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime()+fTimeSignalLength, -a, -b);
      
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
	  
	  a = curSDigit->GetAmp() ;
	  b = a /fTimeSignalLength ;
	  new((*ticks)[contrib++]) AliEMCALTick(curSDigit->GetTime(),0, b);  
	  new((*ticks)[contrib++]) AliEMCALTick(curSDigit->GetTime()+fTimeSignalLength, -a, -b); 

	  *digit = *digit + *curSDigit ;  //add energies

	  index[i]++ ;
	  if( dynamic_cast<TClonesArray *>(sdigArray->At(i))->GetEntriesFast() > index[i] )
	    curSDigit = dynamic_cast<AliEMCALDigit*>(dynamic_cast<TClonesArray *>(sdigArray->At(i))->At(index[i])) ;
	  else
	    curSDigit = 0 ;
	}
      }
      // add fluctuations for photo-electron creation
      amp = sDigitizer->Calibrate(digit->GetAmp()) ; // GeV
      amp *= static_cast<Float_t>(gRandom->Poisson(fMeanPhotonElectron)) / static_cast<Float_t>(fMeanPhotonElectron) ;
  
      //calculate and set time
      Float_t time = FrontEdgeTime(ticks) ;
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
    amp += TMath::Abs(gRandom->Gaus(0., fPinNoise)) ;
    digit->SetAmp(sDigitizer->Digitize(amp)) ;  
    AliDebug(10,Form(" absID %5i amp %f nextSig %5i\n",
		     absID, amp, nextSig));
  } // for(absID = 1; absID <= nEMC; absID++)
  
  ticks->Delete() ;
  delete ticks ;

  delete sdigArray ; //We should not delete its contents

  //remove digits below thresholds
  for(i = 0 ; i < nEMC ; i++){
    digit = dynamic_cast<AliEMCALDigit*>( digits->At(i) ) ;
    Float_t threshold = fDigitThreshold ; 
    if(sDigitizer->Calibrate( digit->GetAmp() ) < threshold)
      digits->RemoveAt(i) ;
    else 
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
  }
  
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ; 
  digits->Expand(ndigits) ;
  
  //Set indexes in list of digits and fill hists.
  AliEMCALHistoUtilities::FillH1(fHists, 0, Double_t(ndigits));
  Float_t energy=0., esum=0.;
  for (i = 0 ; i < ndigits ; i++) { 
    digit = dynamic_cast<AliEMCALDigit *>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ; 
    energy = sDigitizer->Calibrate(digit->GetAmp()) ;
    esum += energy;
    digit->SetAmp(DigitizeEnergy(energy, digit->GetId()) ) ; // for what ??
    AliEMCALHistoUtilities::FillH1(fHists, 2, double(digit->GetAmp()));
    AliEMCALHistoUtilities::FillH1(fHists, 3, double(energy));
    AliEMCALHistoUtilities::FillH1(fHists, 4, double(digit->GetId()));
  }
  AliEMCALHistoUtilities::FillH1(fHists, 1, esum);
}

// //_____________________________________________________________________
Int_t AliEMCALDigitizer::DigitizeEnergy(Float_t energy, Int_t AbsId)
{ 
  // Returns digitized value of the energy in a cell absId
  // Loader
  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>
    (rl->GetDetectorLoader("EMCAL"));
  
  // Load EMCAL Geometry
  rl->LoadgAlice(); 
  AliRun * gAlice = rl->GetAliRun(); 
  AliEMCAL * emcal  = (AliEMCAL*)gAlice->GetDetector("EMCAL");
  AliEMCALGeometry * geom = emcal->GetGeometry();

  if (geom==0)
    AliFatal("Did not get geometry from EMCALLoader") ;

  Int_t iSupMod = -1;
  Int_t nTower  = -1;
  Int_t nIphi   = -1;
  Int_t nIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t channel = -999; 

  Bool_t bCell = geom->GetCellIndex(AbsId, iSupMod, nTower, nIphi, nIeta) ;
  if(!bCell)
    Error("DigitizeEnergy","Wrong cell id number") ;
  geom->GetCellPhiEtaIndexInSModule(iSupMod,nTower,nIphi, nIeta,iphi,ieta);
  
  if(emcalLoader->CalibData()) {
    fADCpedestalEC = emcalLoader->CalibData()
      ->GetADCpedestal(iSupMod,ieta,iphi);
    fADCchannelEC = emcalLoader->CalibData()
      ->GetADCchannel(iSupMod,ieta,iphi);
  }
  
  channel = static_cast<Int_t> (TMath::Ceil( (energy + fADCpedestalEC)/fADCchannelEC ))  ;
  
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

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  // Post Digitizer to the white board
  emcalLoader->PostDigitizer(this) ;
  
  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  else if (fManager) 
    fLastEvent = fFirstEvent ; // what is this ??

  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  Int_t ievent;

  rl->LoadSDigits("EMCAL");
  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    
    rl->GetEvent(ievent);

    Digitize(ievent) ; //Add prepared SDigits to digits and add the noise

    WriteDigits() ;

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
Float_t AliEMCALDigitizer::FrontEdgeTime(TClonesArray * ticks) 
{ 
  //  Returns the shortest time among all time ticks

  ticks->Sort() ; //Sort in accordance with times of ticks
  TIter it(ticks) ;
  AliEMCALTick * ctick = (AliEMCALTick *) it.Next() ;
  Float_t time = ctick->CrossingTime(fTimeThreshold) ;    
  
  AliEMCALTick * t ;  
  while((t=(AliEMCALTick*) it.Next())){
    if(t->GetTime() < time)  //This tick starts before crossing
      *ctick+=*t ;
    else
      return time ;
    
    time = ctick->CrossingTime(fTimeThreshold) ;    
  }
  return time ;
}

//____________________________________________________________________________ 
Bool_t AliEMCALDigitizer::Init()
{
  // Makes all memory allocations
  fInit = kTRUE ; 
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

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

  //PH  Print();
  
  return fInit ;    
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::InitParameters()
{ 
  //parameter initialization for digitizer
  // Tune parameters - 24-nov-04

  fMeanPhotonElectron = 3300 ; // electrons per GeV 
  fPinNoise           = 0.004; 
  if (fPinNoise == 0. ) 
    Warning("InitParameters", "No noise added\n") ; 
  fDigitThreshold     = fPinNoise * 3; // 3 * sigma
  fTimeResolution     = 0.3e-9 ; // 300 psc
  fTimeSignalLength   = 1.0e-9 ;

  fADCchannelEC    = 0.00305; // 200./65536 - width of one ADC channel in GeV
  fADCpedestalEC   = 0.009 ;  // GeV
  fNADCEC          = (Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC - 65536

  fTimeThreshold      = 0.001*10000000 ; // Means 1 MeV in terms of SDigits amplitude ??
  // hists. for control; no hists on default
  fControlHists = 0;
  fHists        = 0;
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
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
  TString fileName( emcalLoader->GetSDigitsFileName() ) ; 
    if ( eventFolderName != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
      fileName = fileName.ReplaceAll(".root", "") + "_" + eventFolderName + ".root" ;
    if ( (gSystem->AccessPathName(fileName)) ) { 
      Error("MixWith", "The file %s does not exist!", fileName.Data()) ;
      return ;
    }
    // need to increase the arrays
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

    AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

    printf("\nWriting digits to %s", emcalLoader->GetDigitsFileName().Data()) ;
    
    printf("\nWith following parameters:\n") ;
    
    printf("    Electronics noise in EMC (fPinNoise) = %f\n", fPinNoise) ;
    printf("    Threshold  in EMC  (fDigitThreshold) = %f\n", fDigitThreshold)  ;
    printf("---------------------------------------------------\n")  ;
  }
  else
    printf("Print: AliEMCALDigitizer not initialized") ; 
}

//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option)
{
  //utility method for printing digit information

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
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
      printf("\n%6d  %8d    %6.5e %4d      %2d : ",
	      digit->GetId(), digit->GetAmp(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
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
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
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

  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

  const TClonesArray * digits = emcalLoader->Digits() ; 
  TTree * treeD = emcalLoader->TreeD(); 
  if ( !treeD ) {
    emcalLoader->MakeDigitsContainer();
    treeD = emcalLoader->TreeD(); 
  }

  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = 0;
  if ((digitsBranch = treeD->GetBranch("EMCAL")))
    digitsBranch->SetAddress(&digits);
  else
    treeD->Branch("EMCAL","TClonesArray",&digits,bufferSize);
  //digitsBranch->SetTitle(fEventFolderName);
  treeD->Fill() ;
  
  emcalLoader->WriteDigits("OVERWRITE");
  emcalLoader->WriteDigitizer("OVERWRITE");

  Unload() ; 

}

void AliEMCALDigitizer::Browse(TBrowser* b)
{
  if(fHists) b->Add(fHists);
  TTask::Browse(b);
}

TList *AliEMCALDigitizer::BookControlHists(int var)
{ 
  // 22-nov-04
  // histograms for monitoring digitizer performance

  Info("BookControlHists"," started ");
  gROOT->cd();
  const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance();
  if(var>=1){
    new TH1F("hDigiN",  "#EMCAL digits with fAmp > fDigitThreshold", 
    fNADCEC+1, -0.5, Double_t(fNADCEC));
    new TH1F("HDigiSumEnergy","Sum.EMCAL energy from digi", 1000, 0.0, 200.);
    new TH1F("hDigiAmp",  "EMCAL digital amplitude", fNADCEC+1, -0.5, Double_t(fNADCEC));
    new TH1F("hDigiEnergy","EMCAL cell energy", 2000, 0.0, 200.);
    new TH1F("hDigiAbsId","EMCAL absId cells with fAmp > fDigitThreshold ",
    geom->GetNCells(), 0.5, Double_t(geom->GetNCells())+0.5);
  }

  fHists = AliEMCALHistoUtilities::MoveHistsToList("EmcalDigiControlHists", kFALSE);
  fHists = 0; //huh? JLK 03-Mar-2006
  return fHists;
}

void AliEMCALDigitizer::SaveHists(const char* name, Bool_t kSingleKey, const char* opt)
{
  AliEMCALHistoUtilities::SaveListOfHists(fHists, name, kSingleKey, opt); 
}
