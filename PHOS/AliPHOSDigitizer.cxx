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
// This is a TTask that constructs SDigits out of Hits
// A Summable Digits is the sum of all hits in a cell
// A threshold is applied 
//
//*-- Author :  Dmitri Peressounko (SUBATECH & Kurchatov Institute) 
//////////////////////////////////////////////////////////////////////////////
// Class performs digitization of Summable digits (in the PHOS case this is just
// sum of contributions of all primary particles into give cell). 
// In addition it performs mixing of summable digits from different events.
// Examples of use:
// root[0] AliPHOSDigitizer * d = new AliPHOSDigitizer() ;
// root[1] d->ExecuteTask()             
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       //Digitizes SDigitis in all events found in file galice.root 
//                       //Depending on variable "CONFIG_SPLIT_FILE" reads branches stored in galice.root
//                       //or in PHOS.SDigits.root
// root[2] AliPHOSDigitizer * d1 = new AliPHOSDigitizer("galice1.root") ;  // Will read sdigits from galice1.root
// root[3] d1->MixWith("galice2.root",1)       // Reads another portion of sdigits from galice2.root
//                                             // says, that this will be output file
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root[3] d1->MixWith("galice3.root",1)       // Reads another portion of sdigits from galice3.root
//                                             // overwrides previous definition of output file
// root[4] d->ExecuteTask()    // Reads SDigits from files galice1.root, galice2.root ....
//                             // mixes them and stores produced Digits in file galice3.root          
//
//
// 

// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSHit.h"
#include "AliPHOSv1.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "TROOT.h"
#include "TFolder.h"

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
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 

}
//____________________________________________________________________________ 
void AliPHOSDigitizer::Init(Int_t isOutFile){
// Mades all memory allocations and defiles, 
// whether first (default) file will be output file (isOutFile !=0) 

  if(!fInitialized){

    fHeaderFiles  = new TClonesArray("TObjString",1) ;
    new((*fHeaderFiles)[0]) TObjString("galice.root") ;
    TFile * file ;

    if(isOutFile)
      file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;
    else
      file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString()) ;

    file->cd() ;
  
    fSDigitsFiles = new TClonesArray("TObjString",1);
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")) 
      new((*fSDigitsFiles)[0]) TObjString("./PHOS.SDigits.root") ;   
    else
      new((*fSDigitsFiles)[0]) TObjString("") ;   
    
    fSDigits      = new TClonesArray("TClonesArray",1) ;
    new((*fSDigits)[0]) TClonesArray("AliPHOSDigit",1000) ;
    
    fDigits = new TClonesArray("AliPHOSDigit",200000) ;
    
    fIevent    = new TArrayI(1) ;
    fIevent->AddAt(-1,0 ) ; 
    fIeventMax = new TArrayI(1) ;

    //Store digits in this file
    if(isOutFile){
      gAlice = (AliRun *) file->Get("gAlice") ;
      fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
      fOutFileNumber = 0 ;
    }
    else{
      TTree * te = (TTree *) file->Get("TE") ;
      fIeventMax->AddAt((Int_t) te->GetEntries(), 0 );
      fOutFileNumber = -1 ;
    }

    fInitialized = kTRUE ;
  }

}


//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(char *HeaderFile,char *DigitsFile = 0):TTask("AliPHOSDigitizer","")
{
  // ctor
  fHeaderFiles  = new TClonesArray("TFile",1) ;          
  new((*fHeaderFiles)[0]) TObjString(HeaderFile) ;
  TFile * file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;      // Header file, where result will be stored

  file->cd() ;
  
  fSDigitsFiles = new TClonesArray("TObjString",1);         // File name of the SDigits branch
  if(DigitsFile)
    new((*fSDigitsFiles)[0]) TObjString(DigitsFile) ;   
  else
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")) 
      new((*fSDigitsFiles)[0]) TObjString("./PHOS.SDigits.root") ;   
    else
      new((*fSDigitsFiles)[0]) TObjString("") ;   
    
  fSDigits      = new TClonesArray("TClonesArray",1) ;      // here list of SDigits wil be stored
  new((*fSDigits)[0]) TClonesArray("AliPHOSDigit",1000) ;
    
  fDigits = new TClonesArray("AliPHOSDigit",200000) ;
  fDigitsFile="PHOS.Digits" ; 
    
  fIevent    = new TArrayI(1) ;
  fIevent->AddAt(-1,0 ) ; 
  fIeventMax = new TArrayI(1) ;
  //Should be check whether gAlice in memory is the same as in file
  //However, there is no such method (?) ==> we are forced to read it
  // if(gAlice->TreeE()==0)

  gAlice = (AliRun *) file->Get("gAlice") ;  //If not read yet

  // Get number of events to process
  fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
  fOutFileNumber = 0 ;

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
  delete fHeaderFiles ;
  delete fSDigitsFiles ;
  delete fSDigits ;
  delete fDigits ;
}
//____________________________________________________________________________
Bool_t AliPHOSDigitizer::Combinator() { 

  //Makes all desirable combinations Signal+Background,
  // returns kFALSE when all combinations are made
  // May be useful to introduce some options like "One-to-One", "All-to-One" and "All-to-All" ?

  //realizing "One-to-One" option...

  if(!fInitialized)
    Init(1) ;

  Int_t inputs ;
  Bool_t endNotReached = kTRUE ;

  for(inputs = 0; (inputs < fNinputs) && endNotReached ;inputs++){
    if(fIevent->At(inputs)+1 < fIeventMax->At(inputs))
      fIevent->AddAt(fIevent->At(inputs)+1, inputs ) ;
    else
      endNotReached = kFALSE ;
  }
  return endNotReached ;

}

//____________________________________________________________________________
void AliPHOSDigitizer::Digitize(Option_t *option) { 

  //Makes the digitization of the collected summable digits

  if(!fInitialized)
    Init(1) ;

  //Collects all hits in the same active volume into digit
  //if(option == "raw")    // add simulated data to row data -- to be implemented


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
  
  Int_t ndigits = fDigits->GetEntries() ;
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

  gAlice->GetEvent(fIevent->At(fOutFileNumber)) ;  // Suitable only for One-To-One mixing
  gAlice->SetEvent(fIevent->At(fOutFileNumber)) ;  // for all-to-all will produce a lot of branches in TreeD

  if(gAlice->TreeD()==0)
    gAlice->MakeTree("D") ;
  
  //Make branches in TreeD for digits and Digitizer
  char branchname[20];
  sprintf(branchname,"PHOS");  

  Int_t bufferSize = 16000 ;
  char * filename = 0;
  if(!fDigitsFile.IsNull())
    filename = (char*) fDigitsFile.Data() ; //ievent ;
  else
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){ //generating file name
      filename = new char[30] ;
      //	sprintf(filename,"PHOS.Digits%d.root",ievent) ;
      sprintf(filename,"PHOS.Digits.root") ;
    }
    else
      filename = 0 ;
  
  //Link digits  
  gAlice->MakeBranchInTree(gAlice->TreeD(),branchname,&fDigits,bufferSize,filename);  
  //Link Digitizer
  AliPHOSDigitizer * d = this ;
  Int_t splitlevel = 0 ;
  sprintf(branchname,"AliPHOSDigitizer");   
  gAlice->MakeBranchInTree(gAlice->TreeD(),branchname,"AliPHOSDigitizer",&d, bufferSize, splitlevel,filename); 

  gAlice->TreeD()->Fill() ;
   
  gAlice->TreeD()->Write(0,kOverwrite) ;  
}

//____________________________________________________________________________
void AliPHOSDigitizer::Exec(Option_t *option) { 
  //manager

  if(!fInitialized)    Init(1) ;
  
  while(Combinator()){  
    
    if(!ReadSDigits()) //read sdigits event(s) evaluated by Combinator() from file(s)
      return ;    
    
    Digitize(option) ; //Add prepared SDigits to digits and add the noise
    WriteDigits() ;
  }

//   //Close all opened files
//   Int_t input ;
//   for(input = 0; input < fNinputs ; input ++){
//     TFile * file = (TFile*) gROOT->GetFile(((TObjString *) fHeaderFiles->At(input))->GetString() ) ;
//     file->Close() ;
//   }

  
}

//__________________________________________________________________
Bool_t AliPHOSDigitizer::ReadSDigits(){
// Reads summable digits from the opened files for the particular set of events given by fIevent

  if(!fInitialized)    Init(1) ;

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
	if( ((TObjString*)fSDigitsFiles->At(inputs))->GetString().CompareTo(sdigitsBranch->GetFileName())==0 ){
	  if( strcmp(sdigitsBranch->GetName(),"PHOS") == 0) {
	    phosNotFound = kFALSE ;
	  }
	}
      }

      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if( ((TObjString*)fSDigitsFiles->At(inputs))->GetString().CompareTo(sdigitizerBranch->GetFileName()) == 0){
	  if( strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0) {
	    sdigitizerNotFound = kFALSE ;
	  }
	}
      }
      
    }

    if(sdigitizerNotFound || phosNotFound){
      cout << "Can't find Branch with sdigits or SDigitizer in the file " ;
      if( ((TObjString*)fSDigitsFiles->At(inputs))->GetString().IsNull() )
	cout << file->GetName() << endl ;	
      else
	cout << ((TObjString*)fSDigitsFiles->At(inputs))->GetString().Data() << endl ;
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
  

  return kTRUE ;

}
//__________________________________________________________________
void AliPHOSDigitizer::MixWith(char* HeaderFile,Int_t isOutFile = 1, char* SDigitsFile =0){
//

  if(!fInitialized)
    if(isOutFile)
      Init(0) ;     //Do not read gAlice from Background file
    else
      Init(1) ;     //read gAlice from background file

  if(HeaderFile == 0){
    cout << "Specify at least header file to merge"<< endl ;
    return ;
  }
  
  Int_t inputs ;
  for(inputs = 0; inputs < fNinputs ; inputs++){
    if(strcmp(((TObjString *)fHeaderFiles->At(inputs))->GetString(),HeaderFile) == 0 ){
      if(SDigitsFile == 0){ 
	if(((TObjString*)fSDigitsFiles->At(inputs))->GetString().CompareTo("")  == 0){
	  cout << "Entry already exists, do not add" << endl ;
	  return ;
	}
      }
      else
	if(((TObjString*)fSDigitsFiles->At(inputs))->GetString().CompareTo(SDigitsFile) == 0){
	cout << "Entry already exists, do not add" << endl ;
	return;
      }
    }	
  }  
  
  fHeaderFiles->Expand(fNinputs+1) ;
  new((*fHeaderFiles)[fNinputs]) TObjString(HeaderFile) ;

  
  TFile * file ;
  if(isOutFile)
    file = new TFile(((TObjString *) fHeaderFiles->At(fNinputs))->GetString(),"update") ;  
  else
    file = new TFile(((TObjString *) fHeaderFiles->At(fNinputs))->GetString()) ;  

  file->cd() ;

  fSDigitsFiles->Expand(fNinputs+1) ;
  new((*fSDigitsFiles)[fNinputs]) TObjString(SDigitsFile) ;

  fSDigits->Expand(fNinputs+1) ;
  new((*fSDigits)[fNinputs]) TClonesArray("AliPHOSDigit",1000) ;

  fIevent->Set(fNinputs+1) ;
  fIevent->AddAt(-1, fNinputs) ;

  fIeventMax->Set(fNinputs+1) ;  

  if(isOutFile){
    gAlice = (AliRun*) file->Get("gAlice") ;
    fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), fNinputs );
    fOutFileNumber = fNinputs ;
  }
  else{
    TTree * te = (TTree *) file->Get("TE") ;
    fIeventMax->AddAt((Int_t) te->GetEntries(), fNinputs );
  }

  fNinputs++ ;

}
//__________________________________________________________________
void AliPHOSDigitizer::Print(Option_t* option){

  if(!fInitialized)    Init(1) ;

  cout << "------------------- "<< GetName() << " -------------" << endl ;
  cout << "Digitizing sDigits from file(s): " <<endl ;
  Int_t input ;
  for(input = 0; input < fNinputs ; input++) {
    cout << "          " << ((TObjString *) fHeaderFiles->At(input))->GetString() << 
      "   Branch: " << ((TObjString *) fSDigitsFiles->At(input))->GetString() << endl ;
  }
  cout << endl ;
  cout << "Writing digits to " << ((TObjString *) fHeaderFiles->At(fOutFileNumber))->GetString() << endl ;

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
