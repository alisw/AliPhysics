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
///_________________________________________________________________________________

// --- ROOT system ---
#include "TTree.h"
#include "TSystem.h"
#include "TBenchmark.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliRunDigitizer.h"
#include "AliEMCALDigit.h"
#include "AliEMCAL.h"
#include "AliEMCALGetter.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTick.h"

ClassImp(AliEMCALDigitizer)


//____________________________________________________________________________ 
  AliEMCALDigitizer::AliEMCALDigitizer():AliDigitizer("",""),
				       fInput(0),
				       fInputFileNames(0x0),
				       fEventNames(0x0)
{
  // ctor
  InitParameters() ; 
  fDefaultInit = kTRUE ; 
  fManager = 0 ;                     // We work in the standalong mode
  fEventFolderName = "" ; 
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(TString alirunFileName, TString eventFolderName):
  AliDigitizer("EMCAL"+AliConfig::Instance()->GetDigitizerTaskName(), alirunFileName),
  fInputFileNames(0), fEventNames(0), fEventFolderName(eventFolderName)
{
  // ctor

  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
  fManager = 0 ;                     // We work in the standalong mode
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(const AliEMCALDigitizer & d) : AliDigitizer(d)
{
  // copyy ctor 

  SetName(d.GetName()) ; 
  SetTitle(d.GetTitle()) ; 
  fDigitThreshold = d.fDigitThreshold ; 
  fMeanPhotonElectron = d.fMeanPhotonElectron ; 
  fPedestal           = d.fPedestal ; 
  fSlope              = d.fSlope ; 
  fPinNoise           = d.fPinNoise ; 
  fTimeResolution     = d.fTimeResolution ; 
  fTimeThreshold      = d.fTimeThreshold ; 
  fTimeSignalLength   = d.fTimeSignalLength ; 
  fADCchannelEC       = d.fADCchannelEC ; 
  fADCpedestalEC      = d.fADCpedestalEC ; 
  fNADCEC             = d.fNADCEC ;
  fEventFolderName    = d.fEventFolderName;
 }

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(AliRunDigitizer * rd):
 AliDigitizer(rd,"EMCAL"+AliConfig::Instance()->GetDigitizerTaskName()),
 fEventFolderName(0)
{
  // ctor Init() is called by RunDigitizer
  fManager = rd ; 
  fEventFolderName = fManager->GetInputFolderName(0) ;
  SetTitle(dynamic_cast<AliStream*>(fManager->GetInputStream(0))->GetFileName(0));
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  AliEMCALGetter * gime =AliEMCALGetter::Instance(GetTitle(),fEventFolderName);
  gime->EmcalLoader()->CleanDigitizer();

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

  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ; 
  Int_t ReadEvent = event ; 
  if (fManager) 
    ReadEvent = dynamic_cast<AliStream*>(fManager->GetInputStream(0))->GetCurrentEventNumber() ; 
  Info("Digitize", "Adding event %d from input stream 0 %s %s", ReadEvent, GetTitle(), fEventFolderName.Data()) ; 
  gime->Event(ReadEvent, "S") ;
  TClonesArray * digits = gime->Digits() ; 
  digits->Clear() ;

  const AliEMCALGeometry *geom = gime->EMCALGeometry() ; 

  Int_t nEMC = geom->GetNPhi()*geom->GetNZ(); //max number of digits possible
  
  Int_t absID ;

  digits->Expand(nEMC) ;

  // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  if ( !gime->SDigitizer() ) 
    gime->LoadSDigitizer();
  AliEMCALSDigitizer * sDigitizer = gime->SDigitizer(); 
  
  if ( !sDigitizer )
    Fatal("Digitize", "SDigitizer with name %s %s not found", fEventFolderName.Data(), GetTitle() ) ; 

  //take all the inputs to add together and load the SDigits
  TObjArray * sdigArray = new TObjArray(fInput) ;
  sdigArray->AddAt(gime->SDigits(), 0) ;
  Int_t i ;
  for(i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i ;
    AliEMCALGetter * gime = AliEMCALGetter::Instance(fInputFileNames[i], tempo) ; 
    if (fManager) 
      ReadEvent = dynamic_cast<AliStream*>(fManager->GetInputStream(i))->GetCurrentEventNumber() ; 
    Info("Digitize", "Adding event %d from input stream %d %s %s", ReadEvent, i, fInputFileNames[i].Data(), tempo.Data()) ; 
    gime->Event(ReadEvent,"S");
    sdigArray->AddAt(gime->SDigits(), i) ;
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
  }

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
  }
  
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
  
  //Set indexes in list of digits
  for (i = 0 ; i < ndigits ; i++) { 
    digit = dynamic_cast<AliEMCALDigit *>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ; 
    Float_t energy = sDigitizer->Calibrate(digit->GetAmp()) ;
    digit->SetAmp(DigitizeEnergy(energy) ) ;
  }
}

//____________________________________________________________________________

Int_t AliEMCALDigitizer::DigitizeEnergy(Float_t energy)
{ 
  Int_t channel = -999; 
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
  
  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle()) ;
   
  if (fLastEvent == -1) 
    fLastEvent = gime->MaxEvent() - 1 ;
  else if (fManager) 
    fLastEvent = fFirstEvent ; 

  Int_t nEvents   = fLastEvent - fFirstEvent + 1;
  
  Int_t ievent ;

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
  
    gime->Event(ievent,"S") ; 

    Digitize(ievent) ; //Add prepared SDigits to digits and add the noise

    WriteDigits() ;

    if(strstr(option,"deb"))
      PrintDigits(option);

    //increment the total number of Digits per run 
    fDigitsInRun += gime->Digits()->GetEntriesFast() ;  
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALDigitizer");
    printf("Exec: took %f seconds for Digitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALDigitizer"), gBenchmark->GetCpuTime("EMCALDigitizer")/nEvents ) ;
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
  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ;  
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object for file %s and event %s !", GetTitle(), fEventFolderName.Data()) ;   
    return kFALSE;
  } 
  TString opt("Digits") ; 
  if(gime->VersionExists(opt) ) { 
    Error( "Init", "Give a version name different from %s", fEventFolderName.Data() ) ;
    fInit = kFALSE ; 
  }

  // Post Digitizer to the white board
  gime->PostDigitizer(this) ;
  
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
  gime->EmcalLoader()->GetDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
  return fInit ;    
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::InitParameters()
{
  fMeanPhotonElectron = 18200 ; // electrons per GeV
  fPinNoise           = 0.003 ; // noise equivalent GeV (random choice)
  if (fPinNoise == 0. ) 
    Warning("InitParameters", "No noise added\n") ; 
  fDigitThreshold     = fPinNoise * 3; //2 sigma
  fTimeResolution     = 1.0e-9 ;
  fTimeSignalLength   = 1.0e-9 ;

  fADCchannelEC    = 0.00305;                     // width of one ADC channel in GeV - HG fix so that we see 200 GeV gammas
  fADCpedestalEC   = 0.005 ;                       // GeV
  fNADCEC          = (Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC

  fTimeThreshold      = 0.001*10000000 ; //Means 1 MeV in terms of SDigits amplitude
 
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
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  TString fileName( gime->GetSDigitsFileName() ) ; 
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
 
//__________________________________________________________________
void AliEMCALDigitizer::Print()const 
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
    for (index = 0 ; index < nStreams ; index++) {  
      TString tempo(fEventNames[index]) ; 
      tempo += index ;
      AliEMCALGetter * gime = AliEMCALGetter::Instance(fInputFileNames[index], tempo) ; 
      TString fileName( gime->GetSDigitsFileName() ) ; 
      if ( fEventNames[index] != AliConfig::GetDefaultEventFolderName()) // only if not the default folder name 
	fileName = fileName.ReplaceAll(".root", "") + "_" + fEventNames[index]  + ".root" ;
      printf ("Adding SDigits from %s %s\n", fInputFileNames[index].Data(), fileName.Data()) ; 
    }
    AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ; 
    printf("\nWriting digits to %s", gime->GetDigitsFileName().Data()) ;
    
    printf("\nWith following parameters:\n") ;
    
    printf("    Electronics noise in EMC (fPinNoise) = %f\n", fPinNoise) ;
    printf("    Threshold  in EMC  (fDigitThreshold) = %f\n", fDigitThreshold)  ;
    printf("---------------------------------------------------\n")  ;
  }
  else
    printf("Print: AliEMCALDigitizer not initialized") ; 
}

//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option){
    
  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ; 
  TClonesArray * digits = gime->Digits() ;
  
  printf("PrintDigits: %d", digits->GetEntriesFast()) ; 
  printf("\nevent %d", gAlice->GetEvNumber()) ;
  printf("\n       Number of entries in Digits list %d", digits->GetEntriesFast() )  ;  
  
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
}

//__________________________________________________________________
Float_t AliEMCALDigitizer::TimeOfNoise(void)
{  // Calculates the time signal generated by noise
  //to be rewritten, now returns just big number
  return 1. ;

}

//__________________________________________________________________
void AliEMCALDigitizer::Unload() 
{  
  // Unloads the SDigits and Digits
  Int_t i ; 
  for(i = 1 ; i < fInput ; i++){
    TString tempo(fEventNames[i]) ; 
    tempo += i ;
    AliEMCALGetter * gime = AliEMCALGetter::Instance(fInputFileNames[i], tempo) ; 
    gime->EmcalLoader()->UnloadSDigits() ; 
  }
  
  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ; 
  gime->EmcalLoader()->UnloadDigits() ; 
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

  AliEMCALGetter * gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName) ; 
  const TClonesArray * digits = gime->Digits() ; 
  TTree * treeD = gime->TreeD(); 

  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = treeD->Branch("EMCAL",&digits,bufferSize);
  digitsBranch->SetTitle(fEventFolderName);
  digitsBranch->Fill() ;
  
  gime->WriteDigits("OVERWRITE");
  gime->WriteDigitizer("OVERWRITE");

  Unload() ; 

}

