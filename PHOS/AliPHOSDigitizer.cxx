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
//*-- Author :  Dmitri Peressounko (SUBATECH & Kurchatov Institute) 
//////////////////////////////////////////////////////////////////////////////
// Class performs digitization of Summable digits (in the PHOS case this is just
// sum of contributions of all primary particles into given cell). 
// In addition it performs mixing of summable digits from different events.
// Examples of use:
// root[0] AliPHOSDigitizer * d = new AliPHOSDigitizer() ;
// root[1] d->ExecuteTask()             
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       //Digitizes SDigitis in all events found in file galice.root 
//                       //Depending on variable "CONFIG_SPLIT_FILE" reads branches stored in galice.root
//                       //or in PHOS.SDigits.root
// root[2] AliPHOSDigitizer * d1 = new AliPHOSDigitizer("galice1.root") ;  
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
//
// For each event two branches are created in TreeD:
//   "PHOS" - list of digits
//   "AliPHOSDigitizer" - AliPHOSDigitizer with all parameters used in digitization
//
// Note, that one can specify new file name for digits branch, and repeat digitization with
// another set of parameters.

// --- ROOT system ---
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TObjString.h"
#include "TBenchmark.h"
// --- Standard library ---
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSDigit.h"
#include "AliPHOSHit.h"
#include "AliPHOSv1.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSDigitizer)


//____________________________________________________________________________ 
  AliPHOSDigitizer::AliPHOSDigitizer():TTask("AliPHOSDigitizer","") 
{
  // ctor

  fSDigitizer = 0 ;
  fNinputs = 1 ;
  fPinNoise = 0.01 ;
  fEMCDigitThreshold = 0.01 ;
  fCPVNoise = 0.01;
  fCPVDigitThreshold = 0.09 ;
  fPPSDNoise = 0.0000001;
  fPPSDDigitThreshold = 0.0000002 ;
  fInitialized = kFALSE ;

  fHeaderFiles = 0;
  fSDigitsTitles = 0;
  fSDigits  = 0 ;
  fDigits = 0;

}
//____________________________________________________________________________ 
void AliPHOSDigitizer::Init(){
// Mades all memory allocations and defiles, 
// whether first (default) file will be output file (isOutFile !=0) 

  if(!fInitialized){
    
    cout << "In Init" << endl ;

    fHeaderFiles  = new TClonesArray("TObjString",1) ;
    new((*fHeaderFiles)[0]) TObjString("galice.root") ;
    
    //Test, if this file already open
    
    TFile *file = (TFile*) gROOT->GetFile(((TObjString *) fHeaderFiles->At(0))->GetString() ) ;
    
    if(file == 0){
      file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }
    else
      file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString()) ;
    
    file->cd() ;
    
    fSDigitsTitles = new TClonesArray("TObjString",1);
    new((*fSDigitsTitles)[0]) TObjString("") ;   
    
    fSDigits      = new TClonesArray("TClonesArray",1) ;
    new((*fSDigits)[0]) TClonesArray("AliPHOSDigit",1000) ;

    fDigitsTitle = "" ;
    
    fDigits = new TClonesArray("AliPHOSDigit",200000) ;
    
    fIevent    = new TArrayI(1) ;
    fIevent->AddAt(-1,0 ) ; 
    fIeventMax = new TArrayI(1) ;
    
    fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
    
    // add Task to //root/Tasks folder
    TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
    roottasks->Add(this) ; 
    
    fInitialized = kTRUE ;
  }
  
}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(const char *HeaderFile,const char *sDigitsTitle):
  TTask("AliPHOSDigitizer","")
{
  // ctor
  fHeaderFiles  = new TClonesArray("TObjString",1) ;          
  new((*fHeaderFiles)[0]) TObjString(HeaderFile) ;
  
  // Header file, where result will be stored
  TFile * file = (TFile*) gROOT->GetFile(((TObjString *) fHeaderFiles->At(0))->GetString() ) ;
  if(file==0){
    file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;      
    gAlice = (AliRun *) file->Get("gAlice") ;  //If not read yet
  }
  
  file->cd() ;
  
  fSDigitsTitles = new TClonesArray("TObjString",1);         // Title name of the SDigits branch
  new((*fSDigitsTitles)[0]) TObjString(sDigitsTitle) ;  
    
  fSDigits      = new TClonesArray("TClonesArray",1) ;      // here list of SDigits wil be stored
  new((*fSDigits)[0]) TClonesArray("AliPHOSDigit",1000) ;
    
  fDigits = new TClonesArray("AliPHOSDigit",200000) ;
  
  fDigitsTitle = "" ; 
  
  fIevent    = new TArrayI(1) ;
  fIevent->AddAt(-1,0 ) ; 
  fIeventMax = new TArrayI(1) ;
  
  // Get number of events to process
  fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
  
  fNinputs = 1 ;
  
  fPinNoise = 0.01 ;
  fEMCDigitThreshold = 0.01 ;
  fCPVNoise = 0.01;
  fCPVDigitThreshold = 0.09 ;
  fPPSDNoise = 0.0000001;
  fPPSDDigitThreshold = 0.0000002 ;  
  
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
  fInitialized = kTRUE ;
  
}

//____________________________________________________________________________ 
  AliPHOSDigitizer::~AliPHOSDigitizer()
{
  // dtor

  if(fHeaderFiles)  delete fHeaderFiles ;
  if(fSDigitsTitles) delete fSDigitsTitles ;
  if(fSDigits)      delete fSDigits ;
  if(fDigits)       delete fDigits ;
}
//____________________________________________________________________________
void AliPHOSDigitizer::Reset() { 
  //sets current event number to the beginning

  if(!fInitialized)
    Init() ;

  Int_t inputs ;
  for(inputs = 0; inputs < fNinputs ;inputs++)
      fIevent->AddAt(-1, inputs ) ;
  
}
//____________________________________________________________________________
Bool_t AliPHOSDigitizer::Combinator() { 

  //Makes all desirable combinations Signal+Background,
  // returns kFALSE when all combinations are made
  // May be useful to introduce some options like "One-to-One", "All-to-One" and "All-to-All" ?

  //realizing "One-to-One" option...

  if(!fInitialized)
    Init() ;

  Int_t inputs ;
  Bool_t endNotReached = kTRUE ;

  for(inputs = 0; (inputs < fNinputs) && endNotReached ;inputs++){
    if(fIevent->At(inputs)+1 < fIeventMax->At(inputs))
      fIevent->AddAt(fIevent->At(inputs)+1, inputs ) ;
    else
      if(inputs == 0)
	endNotReached = kFALSE ;
      else //for inputs other than base one start from the beginning
	fIevent->AddAt(0, inputs ) ;
    
  }
  return endNotReached ;
  
}

//____________________________________________________________________________
void AliPHOSDigitizer::Digitize(Option_t *option) { 

  //Makes the digitization of the collected summable digits

  if(!fInitialized)
    Init() ;

  fDigits->Clear() ;

  AliPHOS * PHOS = (AliPHOS *) gAlice->GetDetector("PHOS") ;   
  AliPHOSGeometry *geom = AliPHOSGeometry::GetInstance( PHOS->GetGeometry()->GetName(), PHOS->GetGeometry()->GetTitle() );

  //Making digits with noise, first EMC
  Int_t nEMC = geom->GetNModules()*geom->GetNPhi()*geom->GetNZ();
  
  Int_t nCPV ;
  Int_t nPPSD ;
  Int_t absID ;
  TString name      =  geom->GetName() ;
  
  if ( name == "IHEP" || name == "MIXT" )    
    nCPV =nEMC + geom->GetNumberOfCPVPadsZ()*geom->GetNumberOfCPVPadsPhi()*
      geom->GetNCPVModules()*geom->GetNumberOfCPVLayers() ;
  else
    nCPV = nEMC; 
  
  if ( name == "GPS2" || name == "MIXT" )    
    nPPSD =nCPV+2*geom->GetNPPSDModules()*geom->GetNumberOfModulesPhi()*geom->GetNumberOfModulesZ()*
      geom->GetNumberOfPadsPhi()*geom->GetNumberOfPadsZ() ;
  else
    nPPSD = nCPV; 


  fDigits->Expand(nPPSD) ;

  
  for(absID = 1; absID <= nEMC; absID++){
    Float_t noise = gRandom->Gaus(0., fPinNoise) ; 
    new((*fDigits)[absID-1]) AliPHOSDigit( -1,absID,fSDigitizer->Digitize(noise) ) ;
  }
  
  for(absID = nEMC+1; absID <= nCPV; absID++){
    Float_t noise = gRandom->Gaus(0., fCPVNoise) ; 
    new((*fDigits)[absID-1]) AliPHOSDigit( -1,absID,fSDigitizer->Digitize(noise) ) ;
  }
  
  for(absID = nCPV+1; absID <= nPPSD; absID++){
    Float_t noise = gRandom->Gaus(0., fPPSDNoise) ; 
    new((*fDigits)[absID-1]) AliPHOSDigit( -1,absID,fSDigitizer->Digitize(noise) ) ;
  }
  

  // Now look throught (unsorted) list of SDigits and add corresponding digits  
  AliPHOSDigit *curSDigit ;
  AliPHOSDigit *digit ;
    
  Int_t inputs;
  for(inputs = 0; inputs< fNinputs ; inputs++){  //loop over (possible) merge sources
    
    TClonesArray * sdigits= (TClonesArray *)fSDigits->At(inputs) ;
    Int_t isdigit ;

    Int_t nSDigits = sdigits->GetEntries() ;     
    for(isdigit=0;isdigit< nSDigits; isdigit++){
      curSDigit = (AliPHOSDigit *)sdigits->At(isdigit) ;
      if(inputs)                                       //Shift primaries for non-background sdigits
	curSDigit->ShiftPrimary(inputs) ;
      digit = (AliPHOSDigit *)fDigits->At(curSDigit->GetId() - 1);
      *digit = *digit + *curSDigit ;
    }  
  }


  //remove digits below thresholds
  for(absID = 0; absID < nEMC ; absID++)
    if(fSDigitizer->Calibrate(((AliPHOSDigit*)fDigits->At(absID))->GetAmp()) < fEMCDigitThreshold)
      fDigits->RemoveAt(absID) ;
  for(absID = nEMC; absID < nCPV ; absID++)
    if(fSDigitizer->Calibrate(((AliPHOSDigit*)fDigits->At(absID))->GetAmp()) < fCPVDigitThreshold)
      fDigits->RemoveAt(absID) ;
  for(absID = nCPV; absID < nPPSD ; absID++)
    if(fSDigitizer->Calibrate(((AliPHOSDigit *)fDigits->At(absID))->GetAmp()) < fPPSDDigitThreshold)
      fDigits->RemoveAt(absID) ;
  
  fDigits->Compress() ;  
  
  Int_t ndigits = fDigits->GetEntriesFast() ;

  fDigits->Expand(ndigits) ;


  //Set indexes in list of digits
  Int_t i ;
  for (i = 0 ; i < ndigits ; i++) { 
    AliPHOSDigit * digit = (AliPHOSDigit *) fDigits->At(i) ; 
    digit->SetIndexInList(i) ;     
  }
}
//____________________________________________________________________________
void AliPHOSDigitizer::WriteDigits(){

  //Made TreeD in the output file if necessary and writes digiths there.

  gAlice->GetEvent(fIevent->At(0)) ;  // Suitable only for One-To-One mixing
  gAlice->SetEvent(fIevent->At(0)) ;  // for all-to-all will produce a lot of branches in TreeD

  if(gAlice->TreeD()==0)
    gAlice->MakeTree("D") ;

  
  //Check, if this branch already exits?
  TBranch * digitsBranch = 0;
  TBranch * digitizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t phosNotFound = kTRUE ;
  Bool_t digitizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(phosNotFound){
      digitsBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp("PHOS",digitsBranch->GetName())==0 ) &&
	  (fDigitsTitle.CompareTo(digitsBranch->GetTitle()) == 0) )
	phosNotFound = kFALSE ;
    }
    if(digitizerNotFound){
      digitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) &&
	  (fDigitsTitle.CompareTo(digitizerBranch->GetTitle()) == 0))
	digitizerNotFound = kFALSE ;
    }
  }
  
  
  if(!(digitizerNotFound && phosNotFound)){ 
    cout << "AliPHOSDigitizer error: " << endl ;
    cout << "       can not update/overwrite existing branches "<< endl ;
    cout << "       do not write " << endl ;
    return ;
  }

  // create new branches

  //First generate file name
  char * file =0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
    file = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(file,"%s/PHOS.Digits.root",gAlice->GetBaseFile()) ;
  }
  
  TDirectory *cwd = gDirectory;
  
  //First create list of sdigits
  Int_t bufferSize = 32000 ;    
  digitsBranch = gAlice->TreeD()->Branch("PHOS",&fDigits,bufferSize);
  digitsBranch->SetTitle(fDigitsTitle.Data());
  if (file) {
    digitsBranch->SetFile(file);
    TIter next( digitsBranch->GetListOfBranches());
    while ((digitsBranch=(TBranch*)next())) {
      digitsBranch->SetFile(file);
    }   
    cwd->cd();
  } 
    
  //second - create Digitizer
  Int_t splitlevel = 0 ;
  AliPHOSDigitizer * d = this ;
  digitizerBranch = gAlice->TreeD()->Branch("AliPHOSDigitizer","AliPHOSDigitizer",
					    &d,bufferSize,splitlevel); 
  digitizerBranch->SetTitle(fDigitsTitle.Data());
  if (file) {
    digitizerBranch->SetFile(file);
    TIter next( digitizerBranch->GetListOfBranches());
    while ((digitizerBranch=(TBranch*)next())) {
      digitizerBranch->SetFile(file);
    }   
    cwd->cd();
  }

  gAlice->TreeD()->Fill() ;
  
  gAlice->TreeD()->Write(0,kOverwrite) ;  
}

//____________________________________________________________________________
void AliPHOSDigitizer::Exec(Option_t *option) { 
  //manager
  if(!fInitialized)    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSDigitizer");

  //reset events numbers to start from the beginnig
  Reset() ;
  
  while(Combinator()){  
    
    if(!ReadSDigits()) //read sdigits event(s) evaluated by Combinator() from file(s)
      return ;    
    
    Digitize(option) ; //Add prepared SDigits to digits and add the noise
    WriteDigits() ;
    
    if(strstr(option,"deb"))
      PrintDigits(option);

  }

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSDigitizer");
    cout << "AliPHOSDigitizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("PHOSDigitizer")/(fIeventMax->At(0)) << " seconds per event " << endl ;
    cout << endl ;
  }
  
}

//__________________________________________________________________
Bool_t AliPHOSDigitizer::ReadSDigits(){
// Reads summable digits from the opened files for the particular set of events given by fIevent

  if(!fInitialized)    Init() ;

  Int_t inputs ;
  for(inputs = fNinputs-1; inputs >= 0; inputs --){

    Int_t event = fIevent->At(inputs) ;

    TFile * file = (TFile*) gROOT->GetFile(((TObjString *) fHeaderFiles->At(inputs))->GetString() ) ;
    file->cd() ;

    // Get SDigits Tree header from file
    char treeName[20]; 
    sprintf(treeName,"TreeS%d",event);
    TTree * treeS = (TTree*)file->Get(treeName);
   
    if(treeS==0){
      cout << "Error at AliPHOSDigitizer: no "<<treeName << "   in file " << file->GetName() << endl ;
      cout << "Do nothing " << endl ;
      return kFALSE ;
    }

    TBranch * sdigitsBranch = 0;
    TBranch * sdigitizerBranch = 0;

    TObjArray * branches = treeS->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t sdigitizerNotFound = kTRUE ;
  
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
            
      if(phosNotFound){
	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
	if(( strcmp("PHOS",sdigitsBranch->GetName())==0 ) &&
	   ((TObjString*) fSDigitsTitles->At(inputs))->GetString().CompareTo(sdigitsBranch->GetTitle())== 0 )
	      phosNotFound = kFALSE ;
	
      }
      
      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if(( strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0) &&
	   ((TObjString*) fSDigitsTitles->At(inputs))->GetString().CompareTo(sdigitizerBranch->GetTitle())== 0 )
	      sdigitizerNotFound = kFALSE ;
	
      }
    }
    
    if(sdigitizerNotFound || phosNotFound){
      cout << "Can't find Branch with sdigits or SDigitizer in the file " ;
      if( ((TObjString*)fSDigitsTitles->At(inputs))->GetString().IsNull() )
	cout << file->GetName() << endl ;	
      else
	cout << ((TObjString*)fSDigitsTitles->At(inputs))->GetString().Data() << endl ;
      cout << "Do nothing" <<endl  ;
      return kFALSE ;
    }
    
    TClonesArray * sdigits = (TClonesArray*) fSDigits->At(inputs) ;  
    sdigitsBranch->SetAddress(&sdigits) ;
    
    AliPHOSSDigitizer *sDigitizer = new AliPHOSSDigitizer();
    sdigitizerBranch->SetAddress(&sDigitizer) ;
    treeS->GetEvent(0) ;
    
    if(fSDigitizer == 0)
      fSDigitizer = sDigitizer ;
    else
      if(!((*fSDigitizer)==(*sDigitizer)) ){
	cout << "ERROR: you are using sdigits made with different SDigitizers" << endl ;
	fSDigitizer->Print("") ;
	sDigitizer->Print("") ;
	cout << "Do Nothing " << endl ;
	return kFALSE ;
      }
    
  }
  fPedestal = fSDigitizer->GetPedestalParameter() ;
  fSlope    = fSDigitizer->GetCalibrationParameter() ;
  
  return kTRUE ;

}
//__________________________________________________________________
void AliPHOSDigitizer::MixWith(char* HeaderFile, char* sDigitsTitle){
//

  if(!fInitialized)
    Init() ;


  if(HeaderFile == 0){
    cout << "Specify at least header file to merge"<< endl ;
    return ;
  }
  
  Int_t inputs ;
  for(inputs = 0; inputs < fNinputs ; inputs++){
    if(strcmp(((TObjString *)fHeaderFiles->At(inputs))->GetString(),HeaderFile) == 0 ){
      if(sDigitsTitle == 0){ 
	if(((TObjString*)fSDigitsTitles->At(inputs))->GetString().CompareTo("")  == 0){
	  cout << "Entry already exists, do not add" << endl ;
	  return ;
	}
      }
      else
	if(((TObjString*)fSDigitsTitles->At(inputs))->GetString().CompareTo(sDigitsTitle)){
	  cout << "Entry already exists, do not add" << endl ;
	  return;
	}
    }	
  }  
  
  fHeaderFiles->Expand(fNinputs+1) ;
  new((*fHeaderFiles)[fNinputs]) TObjString(HeaderFile) ;
  
  
  TFile * file = new TFile(((TObjString *) fHeaderFiles->At(fNinputs))->GetString()) ;  
  
  file->cd() ;
  
  fSDigitsTitles->Expand(fNinputs+1) ;
  new((*fSDigitsTitles)[fNinputs]) TObjString(sDigitsTitle) ;
  
  fSDigits->Expand(fNinputs+1) ;
  new((*fSDigits)[fNinputs]) TClonesArray("AliPHOSDigit",1000) ;
  
  fIevent->Set(fNinputs+1) ;
  fIevent->AddAt(-1, fNinputs) ;
  
  fIeventMax->Set(fNinputs+1) ;  
  
  TTree * te = (TTree *) file->Get("TE") ;
  fIeventMax->AddAt((Int_t) te->GetEntries(), fNinputs );
  
  fNinputs++ ;
  
}
//__________________________________________________________________
void AliPHOSDigitizer::Print(Option_t* option)const {
  
  if(fInitialized){
    
    cout << "------------------- "<< GetName() << " -------------" << endl ;
    cout << "Digitizing sDigits from file(s): " <<endl ;
    Int_t input ;
    for(input = 0; input < fNinputs ; input++) {
      cout << "          " << ((TObjString *) fHeaderFiles->At(input))->GetString() << 
	"   Branch title:" << ((TObjString *) fSDigitsTitles->At(input))->GetString() << endl ;
    }
    cout << endl ;
    cout << "Writing digits to " << ((TObjString *) fHeaderFiles->At(0))->GetString() << endl ;
    
    cout << endl ;
    cout << "With following parameters: " << endl ;
    cout << "     Electronics noise in EMC (fPinNoise) = " << fPinNoise << endl ;
    cout << "  Threshold  in EMC  (fEMCDigitThreshold) = " << fEMCDigitThreshold  << endl ; ;
    cout << "                 Noise in CPV (fCPVNoise) = " << fCPVNoise << endl ; 
    cout << "    Threshold in CPV (fCPVDigitThreshold) = " << fCPVDigitThreshold << endl ; 
    cout << "               Noise in PPSD (fPPSDNoise) = " << fPPSDNoise << endl ;
    cout << "  Threshold in PPSD (fPPSDDigitThreshold) = " << fPPSDDigitThreshold << endl ;
    cout << "---------------------------------------------------" << endl ;
  }
  else
    cout << "AliPHOSDigitizer not initialized " << endl ;
  
}
//__________________________________________________________________
void AliPHOSDigitizer::PrintDigits(Option_t * option){
    
  cout << "AliPHOSDigitiser:"<< endl ;
  cout << "       Number of entries in Digits list " << fDigits->GetEntriesFast() << endl ;
  cout << endl ;
  if(strstr(option,"all")){
    
    //loop over digits
    AliPHOSDigit * digit;
    cout << "Digit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;      
    Int_t index ;
    for (index = 0 ; index < fDigits->GetEntries() ; index++) {
      digit = (AliPHOSDigit * )  fDigits->At(index) ;
      cout << setw(8)  <<  digit->GetId() << " "  << 	setw(3)  <<  digit->GetAmp() <<   "  "  
	   << setw(6)  <<  digit->GetIndexInList() << "  "   
	   << setw(5)  <<  digit->GetNprimary() <<"  ";
      
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++)
	cout << setw(5)  <<  digit->GetPrimary(iprimary+1) << " ";
      cout << endl;  	 
    }
    
  }
}
//__________________________________________________________________
void AliPHOSDigitizer::SetSDigitsBranch(const char* title){
  // we set name of the SDigits branch file in the first! header file
  if(!fInitialized)    Init() ;

  ((TObjString*) fSDigitsTitles->At(0) )->SetString((char*)title) ;

}
//__________________________________________________________________
void AliPHOSDigitizer::SetDigitsBranch(const char* title){
  //Sets the name of the file to which Digits branch will be diverted 
  if(!fInitialized)    Init() ;
  
  fDigitsTitle = title ;

}
