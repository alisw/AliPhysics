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
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TGeometry.h"
#include "TBenchmark.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliHeader.h"
#include "AliStream.h"
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
  AliEMCALDigitizer::AliEMCALDigitizer()
{
  // ctor
  //InitParameters() ; 
  fDefaultInit = kTRUE ; 
  fManager = 0 ;                     // We work in the standalong mode
 
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(const char *headerFile, const char *name, const Bool_t toSplit)
{
  // ctor

  SetTitle(headerFile) ;
  SetName(name) ;
  fManager = 0 ;                     // We work in the standalong mode
  fSplitFile= 0 ; 
  fToSplit = toSplit ;
  Init() ;
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(AliRunDigitizer * ard):AliDigitizer(ard)
{
  // ctor
  SetTitle(ard->GetInputFileName(0,0)) ;
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
  fSplitFile   = 0 ; 
 
  if (ard->GetOutputFile()) {
    SetName(ard->GetOutputFile().Data());
    fToSplit = kTRUE ;
  } else {
    SetName("Default") ;
    fToSplit = kFALSE ;
  }
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  // dtor

    fSplitFile = 0 ; 
}

//____________________________________________________________________________
void AliEMCALDigitizer::Digitize(const Int_t event) 
{ 

  // Makes the digitization of the collected summable digits
  // for this it first creates the array of all EMCAL modules
  // filled with noise (different for EMC, CPV and PPSD) and
  // after that adds contributions from SDigits. This design 
  // helps to avoid scanning over the list of digits to add 
  // contribution of any new SDigit.

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits(GetName()) ; 
  
  digits->Clear() ;

  const AliEMCALGeometry *geom = gime->EMCALGeometry() ; 
  //Making digits with noise, first EMC
  Int_t nEMC = 0 ; 
  if (geom->GetNHCLayers() > 0 )
    nEMC = 3*geom->GetNPhi()*geom->GetNZ(); //max number of digits possible (Preshower, ECAL, HCAL)
  else 
    nEMC = 2*geom->GetNPhi()*geom->GetNZ(); //max number of digits possible (Preshower, ECAL)
 
 Int_t absID ;
  TString name      =  geom->GetName() ;

 // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  const AliEMCALSDigitizer * sDigitizer = gime->SDigitizer(GetName()); 
  if ( !sDigitizer) {
    Fatal("Digitize", "SDigitizer with name %s not found", GetName() ); 
  }

// loop through the sdigits posted to the White Board and add them to the noise
  TCollection * folderslist = gime->SDigitsFolder()->GetListOfFolders() ; 
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TClonesArray * sdigits = 0 ;
  Int_t input = 0 ;
  TObjArray * sdigArray = new TObjArray(2) ;
  while ( (folder = (TFolder*)next()) ) {
    if ( (sdigits = (TClonesArray*)folder->FindObject(GetName()) ) ) {
      TString fileName(folder->GetName()) ;
      fileName.ReplaceAll("_","/") ;
      if (gDebug)
	Info("Digitize", "Adding SDigits %s from %s", GetName(), fileName.Data()) ; 
      sdigArray->AddAt(sdigits, input) ;
      input++ ;
    }
  }

  //Find the first tower with signal
  Int_t nextSig = nEMC + 1 ; 
  Int_t i;
  for(i=0; i<input; i++){
    sdigits = (TClonesArray *)sdigArray->At(i) ;
    if ( !sdigits->GetEntriesFast() )
      continue ; 
    Int_t curNext = ((AliEMCALDigit *)sdigits->At(0))->GetId() ;
     if(curNext < nextSig) 
       nextSig = curNext ;
  }

  TArrayI index(input) ;
  index.Reset() ;  //Set all indexes to zero

  AliEMCALDigit * digit ;
  AliEMCALDigit * curSDigit ;

  TClonesArray * ticks = new TClonesArray("AliEMCALTick",1000) ;

  //Put Noise contribution
  for(absID = 1; absID <= nEMC; absID++){
    Float_t amp = 0 ; 
    //    Float_t noise = TMath::Abs(gRandom->Gaus(0., fPinNoise)); 
    //new((*digits)[absID-1]) AliEMCALDigit( -1, -1, absID, sDigitizer->Digitize(noise), TimeOfNoise() ) ;
    // amplitude set to zero, noise will be added later
    new((*digits)[absID-1]) AliEMCALDigit( -1, -1, absID, 0, TimeOfNoise() ) ;
    //look if we have to add signal?
    digit = (AliEMCALDigit *) digits->At(absID-1) ;
    if(absID==nextSig){
      //Add SDigits from all inputs    
      ticks->Clear() ;
      Int_t contrib = 0 ;
      Float_t a = digit->GetAmp() ;
      Float_t b = TMath::Abs( a /fTimeSignalLength) ;
      //Mark the beginnign of the signal
      new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime(),0, b);  
      //Mark the end of the ignal     
      new((*ticks)[contrib++]) AliEMCALTick(digit->GetTime()+fTimeSignalLength, -a, -b);
      
 // loop over input
      for(i = 0; i< input ; i++){  //loop over (possible) merge sources
    	if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	  curSDigit = (AliEMCALDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
	else
	  curSDigit = 0 ;
	//May be several digits will contribute from the same input
	while(curSDigit && curSDigit->GetId() == absID){	   
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
	  if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	    curSDigit = (AliEMCALDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
	  else
	    curSDigit = 0 ;
	}
      }
      // add fluctuations for photo-electron creation
      amp = sDigitizer->Calibrate(digit->GetAmp()) ; // GeV
      if (gDebug==1) 
	Info("Digitize", "id = %d BEFORE amp = %f\n", digit->GetId(),amp) ; 
      amp *= static_cast<Float_t>(gRandom->Poisson(fMeanPhotonElectron)) / static_cast<Float_t>(fMeanPhotonElectron) ;
      if (gDebug==1) 
	Info("Digitize", "After fluctuation amp = %f\n", amp) ; 
      //calculate and set time
      Float_t time = FrontEdgeTime(ticks) ;
      digit->SetTime(time) ;

      //Find next signal module
      nextSig = nEMC + 1 ;
      for(i=0; i<input; i++){
	sdigits = ((TClonesArray *)sdigArray->At(i)) ;
	Int_t curNext = nextSig ;
	if(sdigits->GetEntriesFast() > index[i] ){
	  curNext = ((AliEMCALDigit *) sdigits->At(index[i]))->GetId() ;
	  
	}
	if(curNext < nextSig) nextSig = curNext ;
      }
    }
    // add the noise now
    
    if (geom->IsInECAL(digit->GetId())) 
      amp += TMath::Abs(gRandom->Gaus(0., fPinNoise)) ;
    else if (geom->IsInPRE(digit->GetId())) 
      amp += TMath::Abs(gRandom->Gaus(0., fPinNoise/100.)) ; // arbitrarely divide by 100, assuming that the gain of APD will be higher
    else if (geom->IsInHCAL(digit->GetId())) 
      amp += TMath::Abs(gRandom->Gaus(0., fPinNoise/10.)) ;  // arbitrarely divide by 10, assuming that the gain of APD will be higher
   if (gDebug==1) 
      Info("Digitize", "After noise amp = %f \n", amp) ; 
    digit->SetAmp(sDigitizer->Digitize(amp)) ;  
  }
  
  ticks->Delete() ;
  delete ticks ;

  delete sdigArray ; //We should not delete its contents

  //remove digits below thresholds
  for(absID = 0; absID < nEMC; absID++){
    digit = dynamic_cast<AliEMCALDigit*>( digits->At(absID) ) ;
    Float_t threshold = 0 ; 

    if (geom->IsInECAL(digit->GetId())) 
      threshold = fDigitThreshold ; 
    else  if (geom->IsInPRE(digit->GetId()))
      threshold = fDigitThreshold / 100. ; // arbitrary see before when noise is added
    else  if (geom->IsInHCAL(digit->GetId()))
      threshold = fDigitThreshold / 10. ; // arbitrary see before when noise is added    

    if(sDigitizer->Calibrate( digit->GetAmp() ) <= threshold)
      digits->RemoveAt(absID) ;
    else {
      if (gDebug == 1)
	Info("Digitize", "id = %d, amp = %f, noise = %f", absID, sDigitizer->Calibrate(digit->GetAmp()), fPinNoise) ; 
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
    }
  }
  
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ; 
  //  digits->Expand(ndigits) ;
  
  //Set indexes in list of digits
  for (i = 0 ; i < ndigits ; i++) { 
    digit = dynamic_cast<AliEMCALDigit *>( digits->At(i) ) ; 
    digit->SetIndexInList(i) ; 
    Float_t energy = sDigitizer->Calibrate(digit->GetAmp()) ;
    digit->SetAmp(DigitizeEnergy(energy,digit->GetId()) ) ;
  }
}

//____________________________________________________________________________

Int_t AliEMCALDigitizer::DigitizeEnergy(Float_t energy, Int_t absId)
{ 
  Int_t channel = -999;
  AliEMCALGeometry * geom = AliEMCALGetter::GetInstance()->EMCALGeometry() ; 
  
  if(geom->IsInPRE(absId)){        //digitize as PRE section
    channel =  static_cast<Int_t>(TMath::Ceil( (energy + fADCpedestalPRE)/fADCchannelPRE ))  ;
    if(channel > fNADCPRE ) 
      channel =  fNADCPRE ;
  }
  else if(geom->IsInECAL(absId)){  //digitize as ECAL section
    channel = static_cast<Int_t> (TMath::Ceil( (energy + fADCpedestalEC)/fADCchannelEC ))  ;
    if(channel > fNADCEC ) 
      channel =  fNADCEC ;
  } 
  else if(geom->IsInHCAL(absId)){  //digitize as HCAL section
    channel = static_cast<Int_t> (TMath::Ceil( (energy + fADCpedestalHC)/fADCchannelHC ))  ;
    if(channel > fNADCHC ) 
      channel =  fNADCHC ;
  }
  
  return channel ;
}

//____________________________________________________________________________
void AliEMCALDigitizer::Exec(Option_t *option) 
{ 
  // Managing method

  if(strcmp(GetName(), "") == 0 )   
    Init() ;
  if (strstr(option,"print")) {
    Print("");
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALDigitizer");

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  
  Int_t nevents ;
  
  TTree * treeD ;
  
  if(fManager){
    treeD = fManager->GetTreeD() ;
    nevents = 1 ;    // Will process only one event
    
    //Check, if this branch already exits
    if (treeD) { 
      TObjArray * lob = (TObjArray*)treeD->GetListOfBranches() ;
      TIter next(lob) ; 
      TBranch * branch = 0 ;  
      Bool_t emcalfound = kFALSE, digitizerfound = kFALSE ; 
      
      while ( (branch = (TBranch*)next()) && (!emcalfound || !digitizerfound) ) {
	if ( (strcmp(branch->GetName(), "EMCAL")==0) && 
	     (strcmp(branch->GetTitle(), GetName())==0) ) 
	  emcalfound = kTRUE ;
	
	else if ( (strcmp(branch->GetName(), "AliEMCALDigitizer")==0) && 
		  (strcmp(branch->GetTitle(), GetName())==0) ) 
	  digitizerfound = kTRUE ; 
      }
      
      if ( emcalfound ) {
	Error( "Exec", "Digits branch with name %s already exits", GetName() ) ;
	return ; 
      }   
      if ( digitizerfound ) {
	Error( "Exec", "Digitizer branch with name %s already exit", GetName() ) ;
	       return ; 
      }
    }   
  }
  else { //EMCAL standalone
    if(gime->BranchExists("Digits") ) 
      return ;
    nevents=gime->MaxEvent() ;
  }
  
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){
    
    if(fManager){

      Int_t input ;
      for(input = 0 ; input < fManager->GetNinputs(); input ++){
  	TTree * treeS = fManager->GetInputTreeS(input) ;
	if(!treeS){
	  Error( "Exec", "No Input") ;
	  return ;
	}
	gime->ReadTreeS(treeS,input) ;
      }

    }
    else 
      gime->Event(ievent,"S") ; 
    
    Digitize(ievent) ; //Add prepared SDigits to digits and add the noise
    
    WriteDigits(ievent) ;
    
    if(strstr(option,"deb"))
      PrintDigits(option);
    
    //increment the total number of Digits per run 
    fDigitsInRun += gime->Digits()->GetEntriesFast() ;  
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALDigitizer");
    Info("Exec", "took %f seconds for Digitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALDigitizer"), gBenchmark->GetCpuTime("EMCALDigitizer")/nevents ) ;
  }
  
}

//____________________________________________________________________________ 
Float_t AliEMCALDigitizer::FrontEdgeTime(TClonesArray * ticks) 
{ // 
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

  if( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(), GetName(), fToSplit) ; 
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object !" ) ; 
    return kFALSE;
  } 
  
  //const AliEMCALGeometry * geom = gime->EMCALGeometry() ;

  //fEmcCrystals = geom->GetNModules() *  geom->GetNCristalsInModule() ;
  
  // Post Digits to the white board
  gime->PostDigits(GetName() ) ;   
  
  // Post Digitizer to the white board
  gime->PostDigitizer(this) ;
 
  fSplitFile = 0 ;
  if(fToSplit){
    // construct the name of the file as /path/EMCAL.SDigits.root
    //First - extract full path if necessary
    TString digitsFileName(GetTitle()) ;
    Ssiz_t islash = digitsFileName.Last('/') ;
    if(islash<digitsFileName.Length())
      digitsFileName.Remove(islash+1,digitsFileName.Length()) ;
    else
      digitsFileName="" ;
    // Next - append the file name 
    digitsFileName+="EMCAL.Digits." ;
    if((strcmp(GetName(),"Default")!=0)&&(strcmp(GetName(),"")!=0)){
      digitsFileName+=GetName() ;
      digitsFileName+="." ;
    }
    digitsFileName+="root" ;
    // Finally - check if the file already opened or open the file
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(digitsFileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(digitsFileName.Data(),"update") ;
  }
  
  //Mark that we will use current header file
  if(!fManager){
    gime->PostSDigits(GetName(),GetTitle()) ;
    gime->PostSDigitizer(GetName(),GetTitle()) ;
  }
  return kTRUE ;    
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::InitParameters()
{
  fMeanPhotonElectron = 1250 ; // electrons per GeV
  
  fPinNoise           = 0.001 ; // noise equivalent GeV (random choice)
  fDigitThreshold     = fPinNoise * 3; //2 sigma
  fTimeResolution     = 0.5e-9 ;
  fTimeSignalLength   = 1.0e-9 ;

  fADCchannelEC    = 0.000220;                     // width of one ADC channel in GeV
  fADCpedestalEC   = 0.005 ;                       // GeV
  fNADCEC          = (Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC

  fADCchannelHC    = 0.000220;                     // width of one ADC channel in GeV
  fADCpedestalHC   = 0.005 ;                       // GeV
  fNADCHC          = (Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC

  fADCchannelPRE   = 0.0000300;                    // width of one ADC channel in Pre Shower
  fADCpedestalPRE  = 0.005 ;                       // GeV 
  fNADCPRE         = (Int_t) TMath::Power(2,12);   // number of channels in Pre ShowerADC

  fTimeThreshold      = 0.001*10000000 ; //Means 1 MeV in terms of SDigits amplitude
 
}

//__________________________________________________________________
void AliEMCALDigitizer::MixWith(char* headerFile)
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
  
  // check if the specified SDigits do not already exist on the White Board:
  // //Folders/RunMC/Event/Data/EMCAL/SDigits/headerFile/sdigitsname

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TString path = gime->SDigitsFolder()->GetName() ; 

  // before it was ???? "Folders/RunMC/Event/Data/EMCAL/SDigits" ; 
  path += headerFile ; 
  path += "/" ; 
  path += GetName() ;
  if ( gROOT->FindObjectAny(path.Data()) ) {
    Error("MixWith", "Entry already exists, do not add" ) ;
    return;
  }

  gime->PostSDigits(GetName(),headerFile) ;
  
  // check if the requested file is already open or exist and if SDigits Branch exist
  TFile * file = (TFile*)gROOT->FindObject(headerFile); 
  if ( !file ) { 
    file = new TFile(headerFile, "READ") ; 
    if (!file) { 
      Error("MixWith", "File %s does not exist!", headerFile) ; 
      return ; 
    }
  }
  
}

//__________________________________________________________________
void AliEMCALDigitizer::Print(Option_t* option)const {

  TString message("\n") ; 

  if( strcmp(GetName(), "") != 0) {
    message += "------------------- " ; 
    message += GetName() ; 
    message += " -------------" ;
    const Int_t nStreams = GetNInputStreams() ; 
    if (nStreams) {
      Int_t index = 0 ;  
      for (index = 0 ; index < nStreams ; index++) {  
	message += "\nAdding SDigits " ; 
	message += GetName() ; 
	message += " from " ; 
	message += fManager->GetInputFileName(index, 0) ;
      } 
     
      message += "\nWriting digits to " ; 
      message += fManager->GetInputFileName(0, 0) ;   
    } else { 
      message += "\nWriting digits to " ;
      message += GetTitle() ;
    }       
    message += "\nWith following parameters: " ;
    message += "\n     Electronics noise in EMC (fPinNoise) = " ; 
    message += fPinNoise ;
    message += "\n  Threshold  in EMC  (fDigitThreshold) = " ; 
    message += fDigitThreshold  ;
    message += "\n---------------------------------------------------"  ;
  }
  else
    message += "\nAliEMCALDigitizer not initialized " ;

  Info("Print", message.Data() ) ; 
}

//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option){
    
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TClonesArray * fDigits = gime->Digits() ;

  TString message("\n") ; 
  message += "       Number of entries in Digits list " ; 
  message += fDigits->GetEntriesFast() ;

  if(strstr(option,"all")){
    
    //loop over digits
    AliEMCALDigit * digit;
    message += "\n   Id  Amplitude    Time          Index Nprim: Primaries list \n" ;    
    Int_t index ;
    char * tempo = new char[8192]; 
    for (index = 0 ; index < fDigits->GetEntries() ; index++) {
      digit = (AliEMCALDigit * )  fDigits->At(index) ;
      sprintf(tempo, "\n%6d  %8d    %6.5e %4d      %2d : ",
	      digit->GetId(), digit->GetAmp(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      message += tempo ; 
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	sprintf(tempo, "%d ",digit->GetPrimary(iprimary+1) ) ; 
	message += tempo ; 
      }
    }   
    delete tempo ; 
  }
  Info("PrintDigits", message.Data() ) ; 
}

//__________________________________________________________________
Float_t AliEMCALDigitizer::TimeOfNoise(void)
{  // Calculates the time signal generated by noise
  //to be rewritten, now returns just big number
  return 1. ;

}

//_________________________________________________________________________________________
void AliEMCALDigitizer::WriteDigits(Int_t event)
{

  // Makes TreeD in the output file. 
  // Check if branch already exists: 
  //   if yes, exit without writing: ROOT TTree does not support overwriting/updating of 
  //      already existing branches. 
  //   else creates branch with Digits, named "EMCAL", title "...",
  //      and branch "AliEMCALDigitizer", with the same title to keep all the parameters
  //      and names of files, from which digits are made.

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  const TClonesArray * digits = gime->Digits(GetName()) ; 
 TTree * treeD ;

 if(fManager)
   treeD = fManager->GetTreeD() ;
 else {
   if ((gAlice->TreeD() == 0) || (fSplitFile)) // we should not create TreeD if it is already here 
     gAlice->MakeTree("D", fSplitFile); // We overwrite TreeD in split file in the case of second reconstruction
   if(fSplitFile)
     fSplitFile->cd() ;
   treeD = gAlice->TreeD();
 }
 
 // -- create Digits branch
 Int_t bufferSize = 32000 ;    
 TBranch * digitsBranch = treeD->Branch("EMCAL",&digits,bufferSize);
 digitsBranch->SetTitle(GetName());
 
 // -- Create Digitizer branch
 Int_t splitlevel = 0 ;
 const  AliEMCALDigitizer * d = gime->Digitizer(GetName()) ;
 TBranch * digitizerBranch = treeD->Branch("AliEMCALDigitizer", "AliEMCALDigitizer", &d,bufferSize,splitlevel); 
 digitizerBranch->SetTitle(GetName());
 
 digitsBranch->Fill() ;
 digitizerBranch->Fill() ; 
 treeD->AutoSave() ;  
 
}

