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
// Class performs digitization of Summable digits (in the EMCAL case this is just
// sum of contributions of all primary particles into given cell). 
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
//

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
#include "AliEMCALDigit.h"
#include "AliEMCALHit.h"
#include "AliEMCALv1.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"

ClassImp(AliEMCALDigitizer)


//____________________________________________________________________________ 
  AliEMCALDigitizer::AliEMCALDigitizer():TTask("AliEMCALDigitizer","") 
{
  // ctor

  fSDigitizer = 0 ;
  fNinputs = 1 ;
  fPinNoise = 0.01 ;
  fEMCDigitThreshold = 0.01 ;
  fInitialized = kFALSE ;

  fHeaderFiles = 0;
  fSDigitsTitles = 0;
  fSDigits  = 0 ;
  fDigits = 0;

}
//____________________________________________________________________________ 
void AliEMCALDigitizer::Init(){
// Makes all memory allocations

  if(!fInitialized){
    
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
    new((*fSDigits)[0]) TClonesArray("AliEMCALDigit",1000) ;

    fSDigitizer = 0 ;
    
    fDigitsTitle = "" ;
    
    fDigits = new TClonesArray("AliEMCALDigit",200000) ;
    
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
AliEMCALDigitizer::AliEMCALDigitizer(const char *headerFile,const char *sDigitsTitle):
  TTask("AliEMCALDigitizer","")
{
  // ctor
  fHeaderFiles  = new TClonesArray("TObjString",1) ;          
  new((*fHeaderFiles)[0]) TObjString(headerFile) ;
  
  // Header file, where result will be stored
  TFile * file = (TFile*) gROOT->GetFile(((TObjString *) fHeaderFiles->At(0))->GetString() ) ;
  if(file==0){
      if(((TObjString *) fHeaderFiles->At(0))->GetString().Contains("rfio"))
	file =	TFile::Open(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;
      else
	file = new TFile(((TObjString *) fHeaderFiles->At(0))->GetString(),"update") ;      
    gAlice = (AliRun *) file->Get("gAlice") ;  //If not read yet
  }
  
  file->cd() ;
  
  fSDigitsTitles = new TClonesArray("TObjString",1);         // Title name of the SDigits branch
  new((*fSDigitsTitles)[0]) TObjString(sDigitsTitle) ;  
    
  fSDigits      = new TClonesArray("TClonesArray",1) ;      // here list of SDigits wil be stored
  new((*fSDigits)[0]) TClonesArray("AliEMCALDigit",1000) ;
    
  fDigits = new TClonesArray("AliEMCALDigit",200000) ;
  
  fDigitsTitle = "" ; 

  
  fSDigitizer = 0 ;
  
  fIevent    = new TArrayI(1) ;
  fIevent->AddAt(-1,0 ) ; 
  fIeventMax = new TArrayI(1) ;
  
  // Get number of events to process
  fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
  
  fNinputs = 1 ;
  
  fPinNoise = 0.01 ;
  fEMCDigitThreshold = 0.01 ;
  
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
  fInitialized = kTRUE ;
  
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  // dtor

  if(fHeaderFiles)  delete fHeaderFiles ;
  if(fSDigitsTitles) delete fSDigitsTitles ;
  if(fSDigits)      delete fSDigits ;
  if(fDigits)       delete fDigits ;
}
//____________________________________________________________________________
void AliEMCALDigitizer::Reset() { 
  //sets current event number to the first simulated event

  if(!fInitialized)
    Init() ;

  Int_t inputs ;
  for(inputs = 0; inputs < fNinputs ;inputs++)
      fIevent->AddAt(-1, inputs ) ;
  
}
//____________________________________________________________________________
Bool_t AliEMCALDigitizer::Combinator() { 

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
void AliEMCALDigitizer::Digitize(Option_t *option) { 

  // Makes the digitization of the collected summable digits
  // for this it first creates the array of all EMCAL modules
  // filled with noise (different for EMC, CPV and PPSD) and
  // after that adds contributions from SDigits. This design 
  // helps to avoid scanning over the list of digits to add 
  // contribution of any new SDigit.

  if(!fInitialized)
    Init() ;

  fDigits->Clear() ;

  AliEMCAL * EMCAL = (AliEMCAL *) gAlice->GetDetector("EMCAL") ;   
  AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance( EMCAL->GetGeometry()->GetName(), EMCAL->GetGeometry()->GetTitle() );

  //Making digits with noise, first EMC
  Int_t nEMC = geom->GetNPhi()*geom->GetNZ();
  Int_t absID ;
  TString name      =  geom->GetName() ;
  
  for(absID = 1; absID <= nEMC; absID++){
    Float_t noise = gRandom->Gaus(0., fPinNoise) ; 
    new((*fDigits)[absID-1]) AliEMCALDigit( -1,-1, absID,fSDigitizer->Digitize(noise) ) ;
  }
  

  // Now look throught (unsorted) list of SDigits and add corresponding digits  
  AliEMCALDigit *curSDigit ;
  AliEMCALDigit *digit ;
    
  Int_t inputs;
  for(inputs = 0; inputs< fNinputs ; inputs++){  //loop over (possible) merge sources
    
    TClonesArray * sdigits= (TClonesArray *)fSDigits->At(inputs) ;
    Int_t isdigit ;

    Int_t nSDigits = sdigits->GetEntries() ;     
    for(isdigit=0;isdigit< nSDigits; isdigit++){
      curSDigit = (AliEMCALDigit *)sdigits->At(isdigit) ;
      if(inputs)                                       //Shift primaries for non-background sdigits
	curSDigit->ShiftPrimary(inputs) ;
      digit = (AliEMCALDigit *)fDigits->At(curSDigit->GetId() - 1);
      *digit = *digit + *curSDigit ;
    }  
  }


  //remove digits below thresholds
  for(absID = 0; absID < nEMC ; absID++)
    if(fSDigitizer->Calibrate(((AliEMCALDigit*)fDigits->At(absID))->GetAmp()) < fEMCDigitThreshold)
      fDigits->RemoveAt(absID) ;
  
  fDigits->Compress() ;  
  
  Int_t ndigits = fDigits->GetEntriesFast() ;

  fDigits->Expand(ndigits) ;


  //Set indexes in list of digits
  Int_t i ;
  for (i = 0 ; i < ndigits ; i++) { 
    AliEMCALDigit * digit = (AliEMCALDigit *) fDigits->At(i) ; 
    digit->SetIndexInList(i) ;     
  }
}
//____________________________________________________________________________
void AliEMCALDigitizer::WriteDigits(){

  // Made TreeD in the output file. Check if branch already exists: if yes, exits 
  // without writing: ROOT TTree does not suppert overwriting/updating of the 
  // already existing branches. Creates branch with Digits, named "EMCAL", title "...",
  // and branch "AliEMCALDigitizer", with the same title to keep all the parameters
  // and names of files, from which digits are made.

  gAlice->GetEvent(fIevent->At(0)) ;  // Suitable only for One-To-One mixing
  gAlice->SetEvent(fIevent->At(0)) ;  // for all-to-all will produce a lot of branches in TreeD

  if(gAlice->TreeD()==0)
    gAlice->MakeTree("D") ;  

  //Check, if this branch already exits?
  TBranch * digitsBranch = 0;
  TBranch * digitizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t emcalNotFound = kTRUE ;
  Bool_t digitizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(emcalNotFound){
      digitsBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp("EMCAL",digitsBranch->GetName())==0 ) &&
	  (fDigitsTitle.CompareTo(digitsBranch->GetTitle()) == 0) )
	emcalNotFound = kFALSE ;
    }
    if(digitizerNotFound){
      digitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (strcmp(digitizerBranch->GetName(),"AliEMCALDigitizer") == 0) &&
	  (fDigitsTitle.CompareTo(digitizerBranch->GetTitle()) == 0))
	digitizerNotFound = kFALSE ;
    }
  }
  
  
  if(!(digitizerNotFound && emcalNotFound)){ 
    cout << "AliEMCALDigitizer error: " << endl ;
    cout << "       can not update/overwrite existing branches "<< endl ;
    cout << "       do not write " << endl ;
    return ;
  }

  // create new branches

  //First generate file name
  char * file =0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
    file = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(file,"%s/EMCAL.Digits.root",gAlice->GetBaseFile()) ;
  }
  
  TDirectory *cwd = gDirectory;
  
  //First create list of sdigits
  Int_t bufferSize = 32000 ;    
  digitsBranch = gAlice->TreeD()->Branch("EMCAL",&fDigits,bufferSize);
  digitsBranch->SetTitle(fDigitsTitle.Data());
  if (file) {
    digitsBranch->SetFile(file);
    TIter next( digitsBranch->GetListOfBranches());
    TBranch * sbr ;
    while ((sbr=(TBranch*)next())) {
      sbr->SetFile(file);
    }   
    cwd->cd();
  } 
    
  //second - create Digitizer
  Int_t splitlevel = 0 ;
  AliEMCALDigitizer * d = this ;
  digitizerBranch = gAlice->TreeD()->Branch("AliEMCALDigitizer","AliEMCALDigitizer",
					    &d,bufferSize,splitlevel); 
  digitizerBranch->SetTitle(fDigitsTitle.Data());
  if (file) {
    digitizerBranch->SetFile(file);
    TIter next( digitizerBranch->GetListOfBranches());
    TBranch * sbr;
    while ((sbr=(TBranch*)next())) {
      sbr->SetFile(file);
    }   
    cwd->cd();
  }

  digitsBranch->Fill() ;
  digitizerBranch->Fill() ;
  
  gAlice->TreeD()->Write(0,kOverwrite) ;  

  //remove fSDigitizer before new event.  
  if(fSDigitizer){
    delete fSDigitizer ;
    fSDigitizer = 0 ;
  }


}

//____________________________________________________________________________
void AliEMCALDigitizer::Exec(Option_t *option) { 
  // Managing method

  if(!fInitialized)    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALDigitizer");

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
    gBenchmark->Stop("EMCALDigitizer");
    cout << "AliEMCALDigitizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("EMCALDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("EMCALDigitizer")/(fIeventMax->At(0)) << " seconds per event " << endl ;
    cout << endl ;
  }
  
}

//__________________________________________________________________
Bool_t AliEMCALDigitizer::ReadSDigits(){
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
      cout << "Error at AliEMCALDigitizer: no "<<treeName << "   in file " << file->GetName() << endl ;
      cout << "Do nothing " << endl ;
      return kFALSE ;
    }

    TBranch * sdigitsBranch = 0;
    TBranch * sdigitizerBranch = 0;

    TObjArray * branches = treeS->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t emcalNotFound = kTRUE ;
    Bool_t sdigitizerNotFound = kTRUE ;
  
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
            
      if(emcalNotFound){
	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
	if(( strcmp("EMCAL",sdigitsBranch->GetName())==0 ) &&
	   ((TObjString*) fSDigitsTitles->At(inputs))->GetString().CompareTo(sdigitsBranch->GetTitle())== 0 )
	      emcalNotFound = kFALSE ;
	
      }
      
      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if(( strcmp(sdigitizerBranch->GetName(),"AliEMCALSDigitizer") == 0) &&
	   ((TObjString*) fSDigitsTitles->At(inputs))->GetString().CompareTo(sdigitizerBranch->GetTitle())== 0 )
	      sdigitizerNotFound = kFALSE ;
	
      }
    }
    
    if(sdigitizerNotFound || emcalNotFound){
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
    
    AliEMCALSDigitizer *sDigitizer = new AliEMCALSDigitizer();
    sdigitizerBranch->SetAddress(&sDigitizer) ;

    sdigitsBranch->GetEntry(0) ;
    sdigitizerBranch->GetEntry(0) ;
    
    if(fSDigitizer == 0)
      fSDigitizer = sDigitizer ;
    else
      if(!((*fSDigitizer)==(*sDigitizer)) ){
	cout << "AliEMCALDigitizer ERROR:" << endl ;
	cout << "       you are using sdigits made with different SDigitizers" << endl ;
	cout << "fSD " << fSDigitizer << "  SD" << sDigitizer << endl ;
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
void AliEMCALDigitizer::MixWith(char* HeaderFile, char* sDigitsTitle){
  // Alows produce digits by superimposing background and signal event.
  // It is assumed, that headers file with SIGNAL events is opened in 
  // constructor, and now we set the BACKGROUND event, with which we 
  // will mix. Thus we avoid writing (changing) huge and expencive 
  // backgound files: all output will be writen into SIGNAL, i.e. 
  // opened in constructor file. 
  //
  // One can open as many files to mix with as one wants.


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
  new((*fSDigits)[fNinputs]) TClonesArray("AliEMCALDigit",1000) ;
  
  fIevent->Set(fNinputs+1) ;
  fIevent->AddAt(-1, fNinputs) ;
  
  fIeventMax->Set(fNinputs+1) ;  
  
  TTree * te = (TTree *) file->Get("TE") ;
  fIeventMax->AddAt((Int_t) te->GetEntries(), fNinputs );
  
  fNinputs++ ;
  
}
//__________________________________________________________________
void AliEMCALDigitizer::Print(Option_t* option)const {
  
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
    cout << "---------------------------------------------------" << endl ;
  }
  else
    cout << "AliEMCALDigitizer not initialized " << endl ;
  
}
//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option){
    
  cout << "AliEMCALDigitiser:"<< endl ;
  cout << "       Number of entries in Digits list " << fDigits->GetEntriesFast() << endl ;
  cout << endl ;
  if(strstr(option,"all")){
    
    //loop over digits
    AliEMCALDigit * digit;
    cout << "Digit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;      
    Int_t index ;
    for (index = 0 ; index < fDigits->GetEntries() ; index++) {
      digit = (AliEMCALDigit * )  fDigits->At(index) ;
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
void AliEMCALDigitizer::SetSDigitsBranch(const char* title){
  // we set title (comment) of the SDigits branch in the first! header file
  if(!fInitialized)    Init() ;

  ((TObjString*) fSDigitsTitles->At(0) )->SetString((char*)title) ;

}
//__________________________________________________________________
void AliEMCALDigitizer::SetDigitsBranch(const char* title){
  //Sets the title (comment) of the branch to which Digits branch
  if(!fInitialized)    Init() ;
  
  fDigitsTitle = title ;

}
