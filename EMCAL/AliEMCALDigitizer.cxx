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
//_________________________________________________________________________________

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
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliHeader.h"
#include "AliEMCALDigit.h"
#include "AliEMCALHit.h"
#include "AliEMCALTick.h"
#include "AliEMCALv1.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGetter.h"
#include "AliRunDigitizer.h"
ClassImp(AliEMCALDigitizer)


//____________________________________________________________________________ 
  AliEMCALDigitizer::AliEMCALDigitizer()
{
  // ctor

  InitParameters() ; 
  fDefaultInit = kTRUE ; 

}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(const char *headerFile,const char *name)
{
  SetName(name) ;
  SetTitle(headerFile) ;
  fManager = 0 ;                     // We work in the standalong mode
  fSplitFile= 0 ; 
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 

}

//____________________________________________________________________________ 
AliEMCALDigitizer::AliEMCALDigitizer(AliRunDigitizer * ard):AliDigitizer(ard)
{
  // ctor
  SetName("Default");    
  SetTitle("aliroot") ;  
}

//____________________________________________________________________________ 
  AliEMCALDigitizer::~AliEMCALDigitizer()
{
  // dtor
  // fDefaultInit = kTRUE if Digitizer created by default ctor (to get just the parameters)
  
  if (!fDefaultInit) {
    AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
    
    // remove the task from the folder list
    gime->RemoveTask("S",GetName()) ;
    gime->RemoveTask("D",GetName()) ;
    
    // remove the Digits from the folder list
    gime->RemoveObjects("D", GetName()) ;
    
    // remove the SDigits from the folder list
    gime->RemoveSDigits() ;
    
    // Delete gAlice
    gime->CloseFile() ; 
    
    fSplitFile = 0 ; 
  }
}

//____________________________________________________________________________ 
void AliEMCALDigitizer::InitParameters()
{
  fSDigitizer = 0 ;
  fNinputs = 1 ;
  fPinNoise = 0.00001 ;
  fTowerDigitThreshold = 0.001 ;
  fTimeResolution     = 0.5e-9 ;
  fTimeSignalLength   = 1.0e-9 ;
  fPreShowerDigitThreshold = fTowerDigitThreshold/25. ;
  fInitialized = kFALSE ;
  fADCchannelTower = 0.000220;       // width of one ADC channel in GeV
  fADCpedestalTower = 0.005 ;      // GeV
  fNADCTower = (Int_t) TMath::Power(2,16) ;  // number of channels in Tower ADC

  fADCchannelPreSho = 0.0000300;          // width of one ADC channel in Pre Shower
  fADCpedestalPreSho = 0.005 ;         // 
  fNADCPreSho = (Int_t) TMath::Power(2,12);      // number of channels in Pre ShowerADC

  fTimeThreshold = 0.001*10000000 ; //Means 1 MeV in terms of SDigits amplitude
 
}

//____________________________________________________________________________ 
Bool_t AliEMCALDigitizer::Init()
{
  // Makes all memory allocations

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(), GetName(), "update") ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliEMCALDigitizer::Init -> Could not obtain the Getter object !" << endl ; 
    return kFALSE;
  } 
  
  //const AliEMCALGeometry * geom = gime->EMCALGeometry() ;
  //fEmcCrystals = geom->GetNModules() *  geom->GetNCristalsInModule() ;
  
  // Post Digits to the white board
  gime->PostDigits(GetName() ) ;   
  
  // Post Digitizer to the white board
  gime->PostDigitizer(this) ;
  
  //Mark that we will use current header file
  if(!fManager){
    gime->PostSDigits(GetName(),GetTitle()) ;
    gime->PostSDigitizer(GetName(),GetTitle()) ;
  }
  return kTRUE ;
    
}

//____________________________________________________________________________
void AliEMCALDigitizer::Reset() { 
  //sets current event number to the first simulated event
if( strcmp(GetName(), "") == 0 )
  Init() ;

      // Int_t inputs ;
      // for(inputs = 0; inputs < fNinputs ;inputs++)
      //  fIevent->AddAt(-1, inputs ) ;
  
}

//____________________________________________________________________________
void AliEMCALDigitizer::Digitize(const Int_t event) { 

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
  Int_t nEMC = 2*geom->GetNPhi()*geom->GetNZ();
  Int_t absID ;
  TString name      =  geom->GetName() ;

 // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  const AliEMCALSDigitizer * sDigitizer = gime->SDigitizer(GetName()); 
  if ( !sDigitizer) {
    cerr << "ERROR: AliEMCALDigitizer::Digitize -> SDigitizer with name " << GetName() << " not found " << endl ; 
    abort() ; 
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
      cout << "INFO: AliEMCALDigitizer::Digitize -> Adding SDigits " 
	   << GetName() << " from " << fileName << endl ; 
      sdigArray->AddAt(sdigits, input) ;
      input++ ;
    }
  }

  //Find the first tower with signal
  Int_t nextSig = 200000 ; 
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

  AliEMCALDigit * digit = 0 ;
  AliEMCALDigit * curSDigit = 0 ;

  TClonesArray * ticks = new TClonesArray("AliEMCALTick",1000) ;

  //Put Noise contribution
  for(absID = 1; absID <= nEMC; absID++){
    Float_t noise = gRandom->Gaus(0., fPinNoise); 
    new((*digits)[absID-1]) AliEMCALDigit( -1, -1, absID,sDigitizer->Digitize(noise), TimeOfNoise() ) ;
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

//calculate and set time
      Float_t time = FrontEdgeTime(ticks) ;
      digit->SetTime(time) ;

      //Find next signal module
      nextSig = 200000 ;
      for(i=0; i<input; i++){
	sdigits = ((TClonesArray *)sdigArray->At(i)) ;
	Int_t curNext = nextSig ;
	if(sdigits->GetEntriesFast() > index[i] ){
	  curNext = ((AliEMCALDigit *) sdigits->At(index[i]))->GetId() ;
	  
	}
	if(curNext < nextSig) nextSig = curNext ;
      }
    }
  }
  
  ticks->Delete() ;
  delete ticks ;




  //remove digits below thresholds
 
  for(absID = 0; absID < nEMC/2 ; absID++){
   if(sDigitizer->Calibrate(((AliEMCALDigit*)digits->At(absID))->GetAmp()) < fTowerDigitThreshold)
      digits->RemoveAt(absID) ;
    else
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
  }
  
  for(absID = nEMC/2; absID < nEMC ; absID++){

    if(sDigitizer->Calibrate(((AliEMCALDigit*)digits->At(absID))->GetAmp()) < fPreShowerDigitThreshold)
      digits->RemoveAt(absID) ;
    else
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
  }
  
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ;
  
  digits->Expand(ndigits) ;
  
  
  //Set indexes in list of digits
  //Int_t i ;
 for (i = 0 ; i < ndigits ; i++) { 
    AliEMCALDigit * digit = (AliEMCALDigit *) digits->At(i) ; 
    digit->SetIndexInList(i) ; 
    Float_t energy = sDigitizer->Calibrate(digit->GetAmp()) ;
    digit->SetAmp(DigitizeEnergy(energy,digit->GetId()) ) ;
  }
}

//____________________________________________________________________________

Int_t AliEMCALDigitizer::DigitizeEnergy(Float_t energy, Int_t absId)
{ 
  Int_t channel = -999;
  Int_t nphi = AliEMCALGetter::GetInstance()->EMCALGeometry()->GetNPhi() ; 
  Int_t nz   = AliEMCALGetter::GetInstance()->EMCALGeometry()->GetNZ() ;
  
  if(absId <= nphi*nz){  //digitize as tower
    channel = static_cast<Int_t> (TMath::Ceil( (energy + fADCpedestalTower)/fADCchannelTower ))  ;
  if(channel > fNADCTower ) 
    channel =  fNADCTower ;
  } else {
    channel =  static_cast<Int_t>(TMath::Ceil( (energy + fADCpedestalPreSho)/fADCchannelPreSho ))  ;
  if(channel > fNADCPreSho ) 
    channel =  fNADCPreSho ;
  }
  
  return channel ;
}

//____________________________________________________________________________
void AliEMCALDigitizer::Exec(Option_t *option) { 
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
  }
  else {
    gAlice->GetEvent(0) ;
    nevents = (Int_t) gAlice->TreeE()->GetEntries() ;
    treeD=gAlice->TreeD() ;
  }
 

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
      cerr << "WARNING: AliEMCALDigitizer -> Digits branch with name " << GetName() 
	   << " already exits" << endl ;
      return ; 
    }   
    if ( digitizerfound ) {
      cerr << "WARNING: AliEMCALDigitizer -> Digitizer branch with name " << GetName() 
	   << " already exits" << endl ;
      return ; 
    }   
  }
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){
    
    if(fManager){
      Int_t input ;
      for(input = 0 ; input < fManager->GetNinputs(); input ++){
  	TTree * treeS = fManager->GetInputTreeS(input) ;
	if(!treeS){
	  cerr << "AliEMCALDigitizer -> No Input " << endl ;
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
    cout << "AliEMCALDigitizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("EMCALDigitizer") << " seconds for Digitizing " 
	 <<  gBenchmark->GetCpuTime("EMCALDigitizer")/nevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
}



//__________________________________________________________________
void AliEMCALDigitizer::MixWith(char* headerFile){
  // Alows produce digits by superimposing background and signal event.
  // It is assumed, that headers file with SIGNAL events is opened in 
  // constructor, and now we set the BACKGROUND event, with which we 
  // will mix. Thus we avoid writing (changing) huge and expencive 
  // backgound files: all output will be writen into SIGNAL, i.e. 
  // opened in constructor file. 
  //
  // One can open as many files to mix with as one wants.

if( strcmp(GetName(), "") == 0 )
    Init() ;
  
 if(fManager){
    cout << "Can not use this method under AliRunDigitizer " << endl ;
    return ;
  } // check if the specified SDigits do not already exist on the White Board:
  // //Folders/RunMC/Event/Data/EMCAL/SDigits/headerFile/sdigitsname

  TString path = "Folders/RunMC/Event/Data/EMCAL/SDigits" ; 
  path += headerFile ; 
  path += "/" ; 
  path += GetName() ;
  if ( gROOT->FindObjectAny(path.Data()) ) {
    cerr << "WARNING: AliEMCALDigitizer::MixWith -> Entry already exists, do not add" << endl ;
    return;
  }

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  gime->PostSDigits(GetName(),headerFile) ;
  
  // check if the requested file is already open or exist and if SDigits Branch exist
  TFile * file = (TFile*)gROOT->FindObject(headerFile); 
  if ( !file ) { 
    file = new TFile(headerFile, "READ") ; 
    if (!file) { 
      cerr << "ERROR: AliEMCALDigitizer::MixWith -> File " << headerFile << " does not exist!" << endl ; 
      return ; 
    }
  }
  
}
 
//__________________________________________________________________
void AliEMCALDigitizer::Print(Option_t* option)const {
  if( strcmp(GetName(), "") != 0) {
    
    cout << "------------------- "<< GetName() << " -------------" << endl ;
    cout << "Digitizing sDigits from file(s): " <<endl ;
    
  TCollection * folderslist = ((TFolder*)gROOT->FindObjectAny("Folders/RunMC/Event/Data/EMCAL/SDigits"))->GetListOfFolders() ; 
    TIter next(folderslist) ; 
    TFolder * folder = 0 ;
    while ( (folder = (TFolder*)next()) ) 
      if ( folder->FindObject(GetName())  ) 
	{
	cout << "Adding SDigits " << GetName() << " from " << folder->GetName() << endl ; 
      cout << endl ;
      cout << "Writing digits to " << GetTitle() << endl ;
      
      cout << endl ;
      cout << "With following parameters: " << endl ;
      cout << "     Electronics noise in EMC (fPinNoise) = " << fPinNoise << endl ;
      cout << "  Threshold  in EMC  (fTowerDigitThreshold) = " << fTowerDigitThreshold  << endl;
      cout << "  Threshold  in PreShower  (fPreShowerDigitThreshold) = " << fPreShowerDigitThreshold  << endl ; ;
      cout << "---------------------------------------------------" << endl ;
    }
    else
      cout << "AliEMCALDigitizer not initialized " << endl ;
    }
}

//__________________________________________________________________
void AliEMCALDigitizer::PrintDigits(Option_t * option){
    
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TClonesArray * fDigits = gime->Digits() ;

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
   if (!gAlice->TreeD() ) 
     gAlice->MakeTree("D",fSplitFile);
   treeD = gAlice->TreeD();
 }
 
 // -- create Digits branch
 Int_t bufferSize = 32000 ;    
 TBranch * digitsBranch = treeD->Branch("EMCAL",&digits,bufferSize);
 digitsBranch->SetTitle(GetName());
 
 // -- Create Digitizer branch
 Int_t splitlevel = 0 ;
 AliEMCALDigitizer * d = gime->Digitizer(GetName()) ;
 TBranch * digitizerBranch = treeD->Branch("AliEMCALDigitizer", "AliEMCALDigitizer", &d,bufferSize,splitlevel); 
 digitizerBranch->SetTitle(GetName());
 
 digitsBranch->Fill() ;
 digitizerBranch->Fill() ; 
 treeD->AutoSave() ;  
 
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
Float_t AliEMCALDigitizer::TimeOfNoise(void)
{  // Calculates the time signal generated by noise
  //to be rewritten, now returns just big number
  return 1. ;

}
//____________________________________________________________________________
void AliEMCALDigitizer::SetSDigitsBranch(const char* title)
{
  // we set title (comment) of the SDigits branch in the first! header file
  if( strcmp(GetName(), "") == 0 )
    Init() ;

  AliEMCALGetter::GetInstance()->SDigits()->SetName(title) ; 
}

//__________________________________________________________________
void AliEMCALDigitizer::SetSplitFile(const TString splitFileName)
{
  // Diverts the Digits in a file separate from the hits file
  
  // I guess it is not going to work if we do merging
  if (fManager) {
    cerr << "ERROR: AliEMCALDigitizer::SetSplitFile -> Not yet available in case of merging activated " << endl ;  
    return ; 
  }

  TDirectory * cwd = gDirectory ;
  if ( !(gAlice->GetTreeDFileName() == splitFileName) ) {
    if (gAlice->GetTreeDFile() )  
      gAlice->GetTreeDFile()->Close() ; 
  }

  fSplitFile = gAlice->InitTreeFile("D",splitFileName.Data());
  fSplitFile->cd() ; 
  gAlice->Write(0, TObject::kOverwrite);
  
  TTree *treeE  = gAlice->TreeE();
  if (!treeE) {
    cerr << "ERROR: AliEMCALDigitizer::SetSplitFile -> No TreeE found "<<endl;
    abort() ;
  }      
  
  // copy TreeE
  AliHeader *header = new AliHeader();
  treeE->SetBranchAddress("Header", &header);
  treeE->SetBranchStatus("*",1);
  TTree *treeENew =  treeE->CloneTree();
  treeENew->Write(0, TObject::kOverwrite);

  // copy AliceGeom
  TGeometry *AliceGeom = static_cast<TGeometry*>(cwd->Get("AliceGeom"));
  if (!AliceGeom) {
    cerr << "ERROR: AliEMCALDigitizer::SetSplitFile -> AliceGeom was not found in the input file "<<endl;
    abort() ;
  }
  AliceGeom->Write(0, TObject::kOverwrite);
  
  gAlice->MakeTree("D",fSplitFile);
  cwd->cd() ; 
  cout << "INFO: AliEMCALDigitizer::SetSPlitMode -> Digits will be stored in " << splitFileName.Data() << endl ; 
}
