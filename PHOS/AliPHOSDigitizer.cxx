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
// This TTask performs digitization of Summable digits (in the PHOS case it is just
// the sum of contributions from all primary particles into a given cell). 
// In addition it performs mixing of summable digits from different events.
// The name of the TTask is also the title of the branch that will contain 
// the created SDigits
// The title of the TTAsk is the name of the file that contains the hits from
// which the SDigits are created
//
// For each event two branches are created in TreeD:
//   "PHOS" - list of digits
//   "AliPHOSDigitizer" - AliPHOSDigitizer with all parameters used in digitization
//
// Note, that one can set a title for new digits branch, and repeat digitization with
// another set of parameters.
//
// Use case:
// root[0] AliPHOSDigitizer * d = new AliPHOSDigitizer() ;
// root[1] d->ExecuteTask()             
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       //Digitizes SDigitis in all events found in file galice.root 
//
// root[2] AliPHOSDigitizer * d1 = new AliPHOSDigitizer("galice1.root") ;  
//                       // Will read sdigits from galice1.root
// root[3] d1->MixWith("galice2.root")       
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//                       // Reads another set of sdigits from galice2.root
// root[3] d1->MixWith("galice3.root")       
//                       // Reads another set of sdigits from galice3.root
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
#include "AliPHOSDigit.h"
#include "AliPHOS.h"
#include "AliPHOSGetter.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSDigitizer)


//____________________________________________________________________________ 
  AliPHOSDigitizer::AliPHOSDigitizer():TTask("","") 
{
  // ctor

  fPinNoise           = 0.01 ;
  fEMCDigitThreshold  = 0.01 ;
  fCPVNoise           = 0.01;
  fCPVDigitThreshold  = 0.09 ;
  fPPSDNoise          = 0.0000001;
  fPPSDDigitThreshold = 0.0000002 ;  

}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(const char *headerFile,const char * name):
  TTask(name, headerFile)
{
  // ctor
   
  fPinNoise           = 0.01 ;
  fEMCDigitThreshold  = 0.01 ;
  fCPVNoise           = 0.01;
  fCPVDigitThreshold  = 0.09 ;
  fPPSDNoise          = 0.0000001;
  fPPSDDigitThreshold = 0.0000002 ;  

  Init() ;
  
}

//____________________________________________________________________________ 
  AliPHOSDigitizer::~AliPHOSDigitizer()
{
  // dtor


}

//____________________________________________________________________________
void AliPHOSDigitizer::Digitize(const Int_t event) 
{ 
  
  // Makes the digitization of the collected summable digits.
  //  It first creates the array of all PHOS modules
  //  filled with noise (different for EMC, CPV and PPSD) and
  //  then adds contributions from SDigits. 
  // This design avoids scanning over the list of digits to add 
  // contribution to new SDigits only.

  if( strcmp(GetName(), "") == 0 )
    Init() ;

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits() ; 

  digits->Clear() ;

  const AliPHOSGeometry *geom = gime->PHOSGeometry() ; 

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


  digits->Expand(nPPSD) ;


  // sdigitize random gaussian noise and add it to all cells (EMCA+CPV+PPSD) 
  // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  const AliPHOSSDigitizer * sDigitizer = gime->SDigitizer(GetName()); 
  if ( !sDigitizer) {
    cerr << "ERROR: AliPHOSDigitizer::Digitize -> SDigitizer with name " << GetName() << " not found " << endl ; 
    abort() ; 
  }
  for(absID = 1; absID <= nEMC; absID++){
    Float_t noise = gRandom->Gaus(0., fPinNoise) ; 
    new((*digits)[absID-1]) AliPHOSDigit( -1,absID,sDigitizer->Digitize(noise) ) ;
  }
  
  for(absID = nEMC+1; absID <= nCPV; absID++){
    Float_t noise = gRandom->Gaus(0., fCPVNoise) ; 
    new((*digits)[absID-1]) AliPHOSDigit( -1,absID,sDigitizer->Digitize(noise) ) ;
  }
  
  for(absID = nCPV+1; absID <= nPPSD; absID++){
    Float_t noise = gRandom->Gaus(0., fPPSDNoise) ; 
    new((*digits)[absID-1]) AliPHOSDigit( -1,absID,sDigitizer->Digitize(noise) ) ;
  }
  
  // loop through the sdigits posted to the White Board and add them to the noise
  TCollection * folderslist = ((TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/SDigits/PHOS"))->GetListOfFolders() ; 
  TIter next(folderslist) ; 
  TFolder * folder = 0 ; 
  TClonesArray * sdigits = 0 ;  
  TString eventS ; 
  eventS += event ;
  while ( (folder = (TFolder*)next()) ) {
   if ( (strcmp(folder->GetTitle(), eventS.Data()) == 0) || (strcmp(folder->GetTitle(), "") == 0) ) {
      Int_t numberoffiles = 0 ; 
      if ( (sdigits = (TClonesArray*)folder->FindObject(GetName()) ) ) {
	cout << "INFO: AliPHOSDigitizer::Exec -> Adding SDigits " << GetName() << " from " << folder->GetName() << endl ; 
	numberoffiles++ ; 
	Int_t index ; 
	AliPHOSDigit * curSDigit ; 
	AliPHOSDigit * digit ; 
	for ( index = 0 ; index < sdigits->GetEntriesFast(); index++) { 
	  curSDigit = (AliPHOSDigit*)sdigits->At(index) ; 
	  curSDigit->ShiftPrimary(numberoffiles) ;
	  digit = (AliPHOSDigit*)digits->At(curSDigit->GetId() - 1 ) ; 
	  *digit = *digit + *curSDigit ; 
	}
      }
    }
  }
  //remove digits below thresholds
  for(absID = 0; absID < nEMC ; absID++)
    if(sDigitizer->Calibrate(((AliPHOSDigit*)digits->At(absID))->GetAmp()) < fEMCDigitThreshold)
      digits->RemoveAt(absID) ;
  
  for(absID = nEMC; absID < nCPV ; absID++)
    if(sDigitizer->Calibrate(((AliPHOSDigit*)digits->At(absID))->GetAmp()) < fCPVDigitThreshold)
      digits->RemoveAt(absID) ;
  
  for(absID = nCPV; absID < nPPSD ; absID++)
    if(sDigitizer->Calibrate(((AliPHOSDigit *)digits->At(absID))->GetAmp()) < fPPSDDigitThreshold)
      digits->RemoveAt(absID) ;
  
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ;

  digits->Expand(ndigits) ;


  //Set indexes in list of digits
  Int_t i ;
  for (i = 0 ; i < ndigits ; i++) { 
    AliPHOSDigit * digit = (AliPHOSDigit *) digits->At(i) ; 
    digit->SetIndexInList(i) ;     
  }

}

//____________________________________________________________________________
void AliPHOSDigitizer::Exec(Option_t *option) 
{ 
  // Managing method

  if( strcmp(GetName(), "") == 0 )    
    Init() ;
  
  if (strstr(option,"print")) {
    Print("");
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSDigitizer");

  //Check, if this branch already exits
  TObjArray * lob = (TObjArray*)gAlice->TreeD()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  Bool_t phosfound = kFALSE, digitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !digitizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) 
      phosfound = kTRUE ;
    
    else if ( (strcmp(branch->GetName(), "AliPHOSDigitizer")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) 
      digitizerfound = kTRUE ; 
  }

  if ( phosfound || digitizerfound ) {
    cerr << "WARNING: AliPHOSDigitizer::WriteDigits -> Digits and/or Digitizer branch with name " << GetName() 
	 << " already exits" << endl ;
    return ; 
  }   

  Int_t nevents = (Int_t) gAlice->TreeE()->GetEntries() ;
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){

    if(!ReadSDigits(ievent)) //read sdigits event(s) evaluated by Combinator() from file(s)
      return ;    
    
    Digitize(ievent) ; //Add prepared SDigits to digits and add the noise
    
    WriteDigits(ievent) ;
  }
   
  if(strstr(option,"deb"))
    PrintDigits(option);
  

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSDigitizer");
    cout << "AliPHOSDigitizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSDigitizer") << " seconds for SDigitizing " 
	 <<  gBenchmark->GetCpuTime("PHOSDigitizer") << " seconds per event " << endl ;
    //	 <<  gBenchmark->GetCpuTime("PHOSDigitizer")/(fIeventMax->At(0)) << " seconds per event " << endl ;
    cout << endl ;
  }
  
}

//____________________________________________________________________________ 
void AliPHOSDigitizer::Init()
{
  // Makes all memory allocations
  // Adds Digitizer task to the folder of PHOS tasks
   //============================================================= YS
  //  The initialisation is now done by AliPHOSGetter
    
  if( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
   
  // the SDigits name is stored by AliPHOSGetter as the name of the TClones Array 
  // //YSAlice/WhiteBoard/SDigits/PHOS/headerFile/branchname and has branchTitle as title.    
    
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), GetName()) ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSDigitizer::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
   
//   fIevent    = new TArrayI(1) ;
//   fIevent->AddAt(-1,0 ) ; 
//   fIeventMax = new TArrayI(1) ;
  
//   fIeventMax->AddAt((Int_t) gAlice->TreeE()->GetEntries(), 0 );
  
  //add Task to //YSAlice/tasks/Digitizer/PHOS
  TTask * aliceSD  = (TTask*)gROOT->FindObjectAny("YSAlice/tasks/Digitizer") ; 
  TTask * phosSD   = (TTask*)aliceSD->GetListOfTasks()->FindObject("PHOS") ;
  phosSD->Add(this) ; 
  // create a folder on the white board //YSAlice/WhiteBoard/Digits/PHOS/headerFile/digitsTitle
  gime->Post(GetTitle(), "D",  GetName() ) ;   

}

//__________________________________________________________________
void AliPHOSDigitizer::MixWith(const char* headerFile)
{
  // Alows to produce digits by superimposing background and signal event.
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

  const char* sDigitsTitle = GetName() ; 
  
  // check if the specified SDigits do not already exist on the White Board:
  // //YSAlice/WhiteBoard/SDigits/PHOS/headerFile/sDigitsTitle

  TString path = "YSAlice/WhiteBoard/SDigits/PHOS/" ; 
  path += headerFile ; 
  path += "/" ; 
  path += sDigitsTitle ;
  if ( gROOT->FindObjectAny(path.Data()) ) {
    cerr << "WARNING: AliPHOSDigitizer::MixWith -> Entry already exists, do not add" << endl ;
    return;
  }
  // check if the requested file is already open or exist and if SDigits Branch exist
  TFile * file = (TFile*)gROOT->FindObject(headerFile); 
  if ( !file ) { 
    file = new TFile(headerFile, "READ") ; 
    if (!file) { 
      cerr << "ERROR: AliPHOSDigitizer::MixWith -> File " << headerFile << " does not exist!" << endl ; 
      return ; 
    }
  }
  Int_t nevent = (Int_t)((TTree*)file->Get("TE"))->GetEntries() ;
  Int_t ievent ; 
  for (ievent = 0; ievent < nevent; ievent++) {
    TString tsname("TreeS") ; 
    tsname += ievent ; 
    TTree * ts = (TTree*)file->Get(tsname.Data()) ;
    if ( !ts ) {
      cerr << "ERROR: AliPHOSDigitizer::MixWith -> TreeS0 " << " does not exist in " << headerFile << endl ; 
      return ;
    }
    
    TObjArray * lob = (TObjArray*)ts->GetListOfBranches() ;
    TIter next(lob) ; 
    TBranch * branch = 0 ; 
    TBranch * sdigitsbranch = 0 ; 
    TBranch * sdigitizerbranch = 0 ; 
    Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
    
    while ( (branch = (TBranch*)next()) && (!phosfound || !sdigitizerfound) ) {
      if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), sDigitsTitle)==0) ) {
	sdigitsbranch = branch ; 
	phosfound = kTRUE ;
      }
      else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && (strcmp(branch->GetTitle(), sDigitsTitle)==0) ) {
	sdigitizerbranch = branch ; 
	sdigitizerfound = kTRUE ; 
      }
    }
    
    if ( !phosfound || !sdigitizerfound ) {
      cerr << "WARNING: AliPHOSDigitizer::MixWith -> Cannot find SDigits and/or SDigitizer with name " << sDigitsTitle << endl ;
      return ; 
    }   
    
    // post the new SDigits to the White Board
    AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
    gime->Post(headerFile, "S", sDigitsTitle, ievent ) ; 
    TClonesArray * sdigits = gime->SDigits(sDigitsTitle, headerFile) ;  
    sdigitsbranch->SetAddress(&sdigits) ;
    sdigitsbranch->GetEvent(0) ;
  }
}

//__________________________________________________________________
void AliPHOSDigitizer::Print(Option_t* option)const {
  // Print Digitizer's parameters
  if( strcmp(GetName(), "") != 0 ){
    
    cout << "------------------- "<< GetName() << " -------------" << endl ;
    cout << "Digitizing sDigits from file(s): " <<endl ;
    
     TCollection * folderslist = ((TFolder*)gROOT->FindObjectAny("YSAlice/WhiteBoard/SDigits/PHOS"))->GetListOfFolders() ; 
    TIter next(folderslist) ; 
    TFolder * folder = 0 ; 
    
    while ( (folder = (TFolder*)next()) ) {
      if ( folder->FindObject(GetName())  ) 
	cout << "Adding SDigits " << GetName() << " from " << folder->GetName() << endl ; 
    }
    cout << endl ;
    cout << "Writing digits to " << GetTitle() << endl ;
    
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
  // Print a table of digits

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits() ; 

  cout << "AliPHOSDigitiser:"<< endl ;
  cout << "       Number of entries in Digits list " << digits->GetEntriesFast() << endl ;
  cout << endl ;
  if(strstr(option,"all")){
    
    //loop over digits
    AliPHOSDigit * digit;
    cout << "Digit Id " << " Amplitude " <<  " Index "  <<  " Nprim " << " Primaries list " <<  endl;      
    Int_t index ;
    for (index = 0 ; index < digits->GetEntries() ; index++) {
      digit = (AliPHOSDigit * )  digits->At(index) ;
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
void AliPHOSDigitizer::SetSDigitsBranch(const char* title)
{
  // we set title (comment) of the SDigits branch in the first! header file
  if( strcmp(GetName(), "") == 0 )
    Init() ;

  AliPHOSGetter::GetInstance()->SDigits()->SetName(title) ; 
 
}

//__________________________________________________________________
Bool_t AliPHOSDigitizer::ReadSDigits(Int_t event)
{
  // Reads summable digits from the opened files for the particular set of events given by fIevent

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  
  TFile * file = (TFile*)gROOT->GetFile(GetTitle()); 
  file->cd() ;

// Get SDigits Tree header from file
  TString treeName("TreeS") ;
  treeName += event ; 
  TTree * treeS = (TTree*)file->Get(treeName.Data());
   
  if(treeS==0){
    cerr << "ERROR: AliPHOSDigitizer::ReadSDigits There is no SDigit Tree" << endl;
    return kFALSE;
  }
 
  //set address of the SDigits and SDigitizer
  TBranch * sdigitsBranch = 0;
  TBranch * sdigitizerBranch = 0;
  TObjArray * lob = (TObjArray*)gAlice->TreeS()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !sdigitizerfound) ) {
   if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) {
      phosfound = kTRUE ;
      sdigitsBranch = branch ; 
    }
    
    else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && (strcmp(branch->GetTitle(), GetName())==0) ) {
      sdigitizerfound = kTRUE ; 
      sdigitizerBranch = branch ;
    }
  }
  if ( !phosfound || !sdigitizerfound ) {
    cerr << "WARNING: AliPHOSDigitizer::ReadSDigits -> Digits and/or Digitizer branch with name " << GetName() 
	 << " not found" << endl ;
    return kFALSE ; 
  }   
  
  
  TClonesArray * sdigits = gime->SDigits() ; 
  sdigitsBranch->SetAddress(&sdigits) ;
  
  AliPHOSSDigitizer * sdigitizer = gime->SDigitizer() ; 
  sdigitizerBranch->SetAddress(&sdigitizer) ;

  sdigitsBranch->GetEntry(0) ;
  sdigitizerBranch->GetEntry(0) ;
 
  fPedestal = sdigitizer->GetPedestalParameter() ;
  fSlope    = sdigitizer->GetCalibrationParameter() ;
  
   return kTRUE ;

}

//____________________________________________________________________________
void AliPHOSDigitizer::Reset() 
{ 
  // sets current event number to the first simulated event

  if( strcmp(GetName(), "") == 0 )
    Init() ;

 //  Int_t inputs ;
//   for(inputs = 0; inputs < fNinputs ;inputs++)
//       fIevent->AddAt(-1, inputs ) ;
  
}

//____________________________________________________________________________
void AliPHOSDigitizer::WriteDigits(Int_t event)
{

  // Makes TreeD in the output file. 
  // Check if branch already exists: 
  //   if yes, exit without writing: ROOT TTree does not support overwriting/updating of 
  //      already existing branches. 
  //   else creates branch with Digits, named "PHOS", title "...",
  //      and branch "AliPHOSDigitizer", with the same title to keep all the parameters
  //      and names of files, from which digits are made.

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits() ; 


  gAlice->GetEvent(event) ; 
  if(gAlice->TreeD()==0)
    gAlice->MakeTree("D") ;  

  // create new branches
  // -- generate file name if necessary
  char * file =0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
    file = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(file,"%s/PHOS.Digits.root",gAlice->GetBaseFile()) ;
  }

  TDirectory *cwd = gDirectory;
  
  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = gAlice->TreeD()->Branch("PHOS",&digits,bufferSize);
  digitsBranch->SetTitle(GetName());
  if (file) {
    digitsBranch->SetFile(file);
    TIter next( digitsBranch->GetListOfBranches());
    TBranch * sbr ;
    while ((sbr=(TBranch*)next())) {
      sbr->SetFile(file);
    }   
    cwd->cd();
  } 
    
  // -- Create Digitizer branch
  Int_t splitlevel = 0 ;
  AliPHOSDigitizer * d = gime->Digitizer(GetName()) ;
  TBranch * digitizerBranch = gAlice->TreeD()->Branch("AliPHOSDigitizer", "AliPHOSDigitizer", &d,bufferSize,splitlevel); 
  digitizerBranch->SetTitle(GetName());
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

}


