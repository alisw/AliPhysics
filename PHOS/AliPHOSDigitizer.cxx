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
#include "TGeometry.h"
#include "TBenchmark.h"

// --- Standard library ---
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliHeader.h"
#include "AliStream.h"
#include "AliRunDigitizer.h"
#include "AliPHOSDigit.h"
#include "AliPHOS.h"
#include "AliPHOSGetter.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSTick.h"

ClassImp(AliPHOSDigitizer)


//____________________________________________________________________________ 
  AliPHOSDigitizer::AliPHOSDigitizer() 
{
  // ctor

  InitParameters() ; 
  fDefaultInit = kTRUE ;
 
  fHitsFileName    = "" ;
  fSDigitsFileName = "" ; 
 }

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(const char *headerFile,const char * name)
{
  // ctor

  SetTitle(headerFile) ;
  SetName(name) ;
  fManager = 0 ;                     // We work in the standalong mode
  fSplitFile= 0 ; 
  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 
  fSDigitsFileName = headerFile ; 
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  gime->Event(0, "S") ; 
  fHitsFileName = gime->SDigitizer()->GetTitle() ; 
}

//____________________________________________________________________________ 
AliPHOSDigitizer::AliPHOSDigitizer(AliRunDigitizer * ard):AliDigitizer(ard)
{
  // ctor
  SetTitle("aliroot") ;
  SetName("Default") ;
  InitParameters() ; 
  fDefaultInit = kTRUE ; 
  
  fSDigitsFileName = fManager->GetInputFileName(0, 0) ;
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(fSDigitsFileName, GetName()) ; 
  gime->Event(0,"S") ; 
  fHitsFileName = gime->SDigitizer()->GetTitle() ; 
  
}

//____________________________________________________________________________ 
  AliPHOSDigitizer::~AliPHOSDigitizer()
{
  // dtor
  // fDefaultInit = kTRUE if Digitizer created by default ctor (to get just the parameters)
  
  if (!fDefaultInit) {
    AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
 
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
void AliPHOSDigitizer::Digitize(const Int_t event) 
{ 
  
  // Makes the digitization of the collected summable digits.
  //  It first creates the array of all PHOS modules
  //  filled with noise (different for EMC, CPV and PPSD) and
  //  then adds contributions from SDigits. 
  // This design avoids scanning over the list of digits to add 
  // contribution to new SDigits only.

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TClonesArray * digits = gime->Digits(GetName()) ; 

  digits->Clear() ;

  const AliPHOSGeometry *geom = gime->PHOSGeometry() ; 
  //Making digits with noise, first EMC
  Int_t nEMC = geom->GetNModules()*geom->GetNPhi()*geom->GetNZ();
  
  Int_t nCPV ;
  Int_t absID ;
  TString name      =  geom->GetName() ;
  
  nCPV = nEMC + geom->GetNumberOfCPVPadsZ()*geom->GetNumberOfCPVPadsPhi()*
    geom->GetNModules() ;

  digits->Expand(nCPV) ;

  // get first the sdigitizer from the tasks list (must have same name as the digitizer)
  const AliPHOSSDigitizer * sDigitizer = gime->SDigitizer(GetName()); 
  if ( !sDigitizer) {
    cerr << "ERROR: AliPHOSDigitizer::Digitize -> SDigitizer with name " << GetName() << " not found " << endl ; 
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
      cout << "INFO: AliPHOSDigitizer::Digitize -> Adding SDigits " 
	   << GetName() << " from " << fileName << endl ; 
      sdigArray->AddAt(sdigits, input) ;
      input++ ;
    }
  }

  //Find the first crystall with signal
  Int_t nextSig = 200000 ; 
  Int_t i;
  for(i=0; i<input; i++){
    sdigits = (TClonesArray *)sdigArray->At(i) ;
    if ( !sdigits->GetEntriesFast() )
      continue ; 
    Int_t curNext = ((AliPHOSDigit *)sdigits->At(0))->GetId() ;
     if(curNext < nextSig) 
       nextSig = curNext ;
  }

  TArrayI index(input) ;
  index.Reset() ;  //Set all indexes to zero

  AliPHOSDigit * digit ;
  AliPHOSDigit * curSDigit ;

  TClonesArray * ticks = new TClonesArray("AliPHOSTick",1000) ;

  //Put Noise contribution
  for(absID = 1; absID <= nEMC; absID++){
    Float_t noise = gRandom->Gaus(0., fPinNoise) ; 
    new((*digits)[absID-1]) AliPHOSDigit( -1,absID,sDigitizer->Digitize(noise), TimeOfNoise() ) ;
    //look if we have to add signal?
    digit = (AliPHOSDigit *) digits->At(absID-1) ;
 
    if(absID==nextSig){
      //Add SDigits from all inputs 
      ticks->Clear() ;
      Int_t contrib = 0 ;
      Float_t a = digit->GetAmp() ;
      Float_t b = TMath::Abs( a /fTimeSignalLength) ;
      //Mark the beginnign of the signal
      new((*ticks)[contrib++]) AliPHOSTick(digit->GetTime(),0, b);  
      //Mark the end of the ignal     
      new((*ticks)[contrib++]) AliPHOSTick(digit->GetTime()+fTimeSignalLength, -a, -b); 

      //loop over inputs
      for(i=0; i<input; i++){
	if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	  curSDigit = (AliPHOSDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
	else
	  curSDigit = 0 ;
	//May be several digits will contribute from the same input
	while(curSDigit && curSDigit->GetId() == absID){	   
	  //Shift primary to separate primaries belonging different inputs
	  Int_t primaryoffset ;
	  if(fManager)
	    primaryoffset = fManager->GetMask(i) ; 
	  else
	    primaryoffset = 10000000*i ;
	  curSDigit->ShiftPrimary(primaryoffset) ;
	  
	  a = curSDigit->GetAmp() ;
	  b = a /fTimeSignalLength ;
	  new((*ticks)[contrib++]) AliPHOSTick(curSDigit->GetTime(),0, b);  
	  new((*ticks)[contrib++]) AliPHOSTick(curSDigit->GetTime()+fTimeSignalLength, -a, -b); 

	  *digit = *digit + *curSDigit ;  //add energies

	  index[i]++ ;
	  if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	    curSDigit = (AliPHOSDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
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
	  curNext = ((AliPHOSDigit *) sdigits->At(index[i]))->GetId() ;
	  
	}
	if(curNext < nextSig) nextSig = curNext ;
      }
    }
  }
  
  ticks->Delete() ;
  delete ticks ;


  //Now CPV digits (different noise and no timing)
  for(absID = nEMC+1; absID <= nCPV; absID++){
    Float_t noise = gRandom->Gaus(0., fCPVNoise) ; 
    new((*digits)[absID-1]) AliPHOSDigit( -1,absID,sDigitizer->Digitize(noise), TimeOfNoise() ) ;
    //look if we have to add signal?
    if(absID==nextSig){
      digit = (AliPHOSDigit *) digits->At(absID-1) ;
      //Add SDigits from all inputs
      for(i=0; i<input; i++){
	if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	  curSDigit = (AliPHOSDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
	else
	  curSDigit = 0 ;

	//May be several digits will contribute from the same input
	while(curSDigit && curSDigit->GetId() == absID){	   
	  //Shift primary to separate primaries belonging different inputs
	  Int_t primaryoffset ;
	  if(fManager)
	    primaryoffset = fManager->GetMask(i) ; 
	  else
	    primaryoffset = 10000000*i ;
	  curSDigit->ShiftPrimary(primaryoffset) ;

	  //add energies
	  *digit = *digit + *curSDigit ;  
	  index[i]++ ;
	  if(((TClonesArray *)sdigArray->At(i))->GetEntriesFast() > index[i] )
	    curSDigit = (AliPHOSDigit*)((TClonesArray *)sdigArray->At(i))->At(index[i]) ; 	
	  else
	    curSDigit = 0 ;
	}
      }

      //Find next signal module
      nextSig = 200000 ;
      for(i=0; i<input; i++){
	sdigits = (TClonesArray *)sdigArray->At(i) ;
	Int_t curNext = nextSig ;
	if(sdigits->GetEntriesFast() > index[i] )
	  curNext = ((AliPHOSDigit *) sdigits->At(index[i]))->GetId() ;
	if(curNext < nextSig) nextSig = curNext ;
      }
      
    }
  }
  delete sdigArray ; //We should not delete its contents
  
  
  
  //remove digits below thresholds
  for(i = 0; i < nEMC ; i++){
    digit = (AliPHOSDigit*) digits->At(i) ;
    if(sDigitizer->Calibrate( digit->GetAmp() ) < fEMCDigitThreshold)
      digits->RemoveAt(i) ;
    else
      digit->SetTime(gRandom->Gaus(digit->GetTime(),fTimeResolution) ) ;
  }


  for(i = nEMC; i < nCPV ; i++)
    if(sDigitizer->Calibrate(((AliPHOSDigit*)digits->At(i))->GetAmp()) < fCPVDigitThreshold)
      digits->RemoveAt(i) ;
    
  digits->Compress() ;  
  
  Int_t ndigits = digits->GetEntriesFast() ;
  digits->Expand(ndigits) ;

  //Set indexes in list of digits and make true digitization of the energy
  for (i = 0 ; i < ndigits ; i++) { 
    AliPHOSDigit * digit = (AliPHOSDigit *) digits->At(i) ; 
    digit->SetIndexInList(i) ;     
    Float_t energy = sDigitizer->Calibrate(digit->GetAmp()) ;
    digit->SetAmp(DigitizeEnergy(energy,digit->GetId()) ) ;
  }

}
//____________________________________________________________________________
Int_t AliPHOSDigitizer::DigitizeEnergy(Float_t energy, Int_t absId)
{
  Int_t chanel ;
  if(absId <= fEmcCrystals){ //digitize as EMC 
    chanel = (Int_t) TMath::Ceil((energy - fADCpedestalEmc)/fADCchanelEmc) ;       
    if(chanel > fNADCemc ) chanel =  fNADCemc ;
  }
  else{ //Digitize as CPV
    chanel = (Int_t) TMath::Ceil((energy - fADCpedestalCpv)/fADCchanelCpv) ;       
    if(chanel > fNADCcpv ) chanel =  fNADCcpv ;
  }
  return chanel ;
}
//____________________________________________________________________________
void AliPHOSDigitizer::Exec(Option_t *option) 
{ 
  // Managing method

  if(strcmp(GetName(), "") == 0 )   
    Init() ;
  if (strstr(option,"print")) {
    Print("");
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSDigitizer");

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  
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
    Bool_t phosfound = kFALSE, digitizerfound = kFALSE ; 
    
    while ( (branch = (TBranch*)next()) && (!phosfound || !digitizerfound) ) {
      if ( (strcmp(branch->GetName(), "PHOS")==0) && 
	   (strcmp(branch->GetTitle(), GetName())==0) ) 
	phosfound = kTRUE ;
      
      else if ( (strcmp(branch->GetName(), "AliPHOSDigitizer")==0) && 
		(strcmp(branch->GetTitle(), GetName())==0) ) 
	digitizerfound = kTRUE ; 
    }
    
    if ( phosfound ) {
      cerr << "WARNING: AliPHOSDigitizer -> Digits branch with name " << GetName() 
	   << " already exits" << endl ;
      return ; 
    }   
    if ( digitizerfound ) {
      cerr << "WARNING: AliPHOSDigitizer -> Digitizer branch with name " << GetName() 
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
	  cerr << "AliPHOSDigitizer::Exec -> No Input " << endl ;
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
    gBenchmark->Stop("PHOSDigitizer");
    cout << "AliPHOSDigitizer:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSDigitizer") << " seconds for Digitizing " 
	 <<  gBenchmark->GetCpuTime("PHOSDigitizer")/nevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
}

//____________________________________________________________________________ 
Float_t AliPHOSDigitizer::FrontEdgeTime(TClonesArray * ticks) 
{ // 
  ticks->Sort() ; //Sort in accordance with times of ticks
  TIter it(ticks) ;
  AliPHOSTick * ctick = (AliPHOSTick *) it.Next() ;
  Float_t time = ctick->CrossingTime(fTimeThreshold) ;    

  AliPHOSTick * t ;  
  while((t=(AliPHOSTick*) it.Next())){
    if(t->GetTime() < time)  //This tick starts before crossing
      *ctick+=*t ;
    else
      return time ;

    time = ctick->CrossingTime(fTimeThreshold) ;    
  }
  return time ;
}

//____________________________________________________________________________ 
Bool_t AliPHOSDigitizer::Init()
{

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), GetName()) ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSDigitizer::Init -> Could not obtain the Getter object !" << endl ; 
    return kFALSE;
  } 
  
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ;

  fEmcCrystals = geom->GetNModules() *  geom->GetNCristalsInModule() ;
  
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
void AliPHOSDigitizer::InitParameters()
{
  fPinNoise           = 0.004 ;
  fEMCDigitThreshold  = 0.012 ;
  fCPVNoise           = 0.01;
  fCPVDigitThreshold  = 0.09 ;
  fTimeResolution     = 0.5e-9 ;
  fTimeSignalLength   = 1.0e-9 ;
  fDigitsInRun  = 0 ; 
  fADCchanelEmc = 0.0015;        // width of one ADC channel in GeV
  fADCpedestalEmc = 0.005 ;      //
  fNADCemc = (Int_t) TMath::Power(2,16) ;  // number of channels in EMC ADC

  fADCchanelCpv = 0.0012 ;          // width of one ADC channel in CPV 'popugais'
  fADCpedestalCpv = 0.012 ;         // 
  fNADCcpv = (Int_t) TMath::Power(2,12);      // number of channels in CPV ADC

  fTimeThreshold = 0.001*10000000 ; //Means 1 MeV in terms of SDigits amplitude
    
}

//__________________________________________________________________
void AliPHOSDigitizer::MixWith(const char* headerFile)
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
    cout << "Can not use this method under AliRunDigitizer " << endl ;
    return ;
  }
  
  // check if the specified SDigits do not already exist on the White Board:
  // //Folders/RunMC/Event/Data/PHOS/SDigits/headerFile/sdigitsname

  TString path = "Folders/RunMC/Event/Data/PHOS/SDigits" ; 
  path += headerFile ; 
  path += "/" ; 
  path += GetName() ;
  if ( gROOT->FindObjectAny(path.Data()) ) {
    cerr << "WARNING: AliPHOSDigitizer::MixWith -> Entry already exists, do not add" << endl ;
    return;
  }

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  gime->PostSDigits(GetName(),headerFile) ;
  
  // check if the requested file is already open or exist and if SDigits Branch exist
  TFile * file = (TFile*)gROOT->FindObject(headerFile); 
  if ( !file ) { 
    file = new TFile(headerFile, "READ") ; 
    if (!file) { 
      cerr << "ERROR: AliPHOSDigitizer::MixWith -> File " << headerFile << " does not exist!" << endl ; 
      return ; 
    }
  }
  
}

//__________________________________________________________________
void AliPHOSDigitizer::SetSplitFile(const TString splitFileName)
{
  // Diverts the Digits in a file separate from the hits file
  
  // I guess it is not going to work if we do merging
//   if (fManager) {
//     cerr << "ERROR: AliPHOSDigitizer::SetSplitFile -> Not yet available in case of merging activated " << endl ;  
//     return ; 
//   }

  cout << "AliPHOSDigitizer::SetSplitFile " << gDirectory->GetName() << endl ;  
  cout << "AliPHOSDigitizer::SetSplitFile " << gAlice->GetTreeDFileName() << endl ;  
  cout << "AliPHOSDigitizer::SetSplitFile " <<  splitFileName.Data() << endl ;
  
  SetTitle(splitFileName) ; 

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
    cerr << "ERROR: AliPHOSDigitizer::SetSplitFile -> No TreeE found "<<endl;
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
    cerr << "ERROR: AliPHOSDigitizer::SetSplitFile -> AliceGeom was not found in the input file "<<endl;
    abort() ;
    }
  AliceGeom->Write(0, TObject::kOverwrite);
  
  gAlice->MakeTree("D",fSplitFile);
  cwd->cd() ; 
  cout << "INFO: AliPHOSDigitizer::SetSPlitMode -> Digits will be stored in " << splitFileName.Data() << endl ; 
}

//__________________________________________________________________
void AliPHOSDigitizer::Print(Option_t* option)const {
  // Print Digitizer's parameters
  if( strcmp(GetName(), "") != 0 ){
    
    cout << "------------------- "<< GetName() << " -------------" << endl ;

    const Int_t nStreams = GetNInputStreams() ; 
    if (nStreams) {
      Int_t index = 0 ;  
      for (index = 0 ; index < nStreams ; index++)  
	cout << "Adding SDigits " << GetName() << " from " <<  fManager->GetInputFileName(index, 0) << endl ; 
      
      cout << endl ;
      cout << "Writing digits to " <<   fManager->GetInputFileName(0, 0) << endl ;   
    } else { 
//       AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;  
//       gime->Folder("sdigits")  ;
//       cout << "Digitizing sDigits from file(s): " <<endl ;
//       TCollection * folderslist = gime->Folder("sdigits")->GetListOfFolders() ; 
//       TIter next(folderslist) ; 
//       TFolder * folder = 0 ; 
      
//       while ( (folder = (TFolder*)next()) ) {
// 	if ( folder->FindObject(GetName())  ) 
      cout << "Adding SDigits " << GetName() << " from " << GetSDigitsFileName() << endl ; 
//      }
      cout << endl ;
      cout << "Writing digits to " << GetTitle() << endl ;
    }    
    cout << endl ;
    cout << "With following parameters: " << endl ;
    cout << "     Electronics noise in EMC (fPinNoise) = " << fPinNoise << endl ;
    cout << "  Threshold  in EMC  (fEMCDigitThreshold) = " << fEMCDigitThreshold  << endl ; ;
    cout << "                 Noise in CPV (fCPVNoise) = " << fCPVNoise << endl ; 
    cout << "    Threshold in CPV (fCPVDigitThreshold) = " << fCPVDigitThreshold << endl ; 
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

  cout << "AliPHOSDigitiser: event " << gAlice->GetEvNumber() << endl ;
  cout << "       Number of entries in Digits list " << digits->GetEntriesFast() << endl ;
  cout << endl ;
  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliPHOSDigit * digit;
    cout << "EMC digits (with primaries): " << endl ;
    cout << "Digit Id     Amplitude   Index       Nprim  Primaries list " <<  endl;      
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; (index < digits->GetEntriesFast()) && 
	 (((AliPHOSDigit * )  digits->At(index))->GetId() <= maxEmc) ; index++) {
      digit = (AliPHOSDigit * )  digits->At(index) ;
      if(digit->GetNprimary() == 0) continue;
      cout << setw(6)  <<  digit->GetId() << "   "  << 	setw(10)  <<  digit->GetAmp() <<   "    "  
	   << setw(6)  <<  digit->GetIndexInList() << "    "   
	   << setw(5)  <<  digit->GetNprimary() <<"    ";
      
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++)
	cout << setw(5)  <<  digit->GetPrimary(iprimary+1) << "    ";
      cout << endl;  	 
    }    
    cout << endl;
  }

  if(strstr(option,"all")||strstr(option,"CPV")){
    
    //loop over CPV digits
    AliPHOSDigit * digit;
    cout << "CPV digits: " << endl ;
    cout << "Digit Id     Amplitude   Index       Nprim  Primaries list " <<  endl;      
    Int_t maxEmc = gime->PHOSGeometry()->GetNModules()*gime->PHOSGeometry()->GetNCristalsInModule() ;
    Int_t index ;
    for (index = 0 ; index < digits->GetEntriesFast(); index++) {
      digit = (AliPHOSDigit * )  digits->At(index) ;
      if(digit->GetId() > maxEmc){
	cout << setw(6)  <<  digit->GetId() << "   "  << 	setw(10)  <<  digit->GetAmp() <<   "    "  
	     << setw(6)  <<  digit->GetIndexInList() << "    "   
	     << setw(5)  <<  digit->GetNprimary() <<"    ";
	
	Int_t iprimary;
	for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++)
	  cout << setw(5)  <<  digit->GetPrimary(iprimary+1) << "    ";
	cout << endl;  	 
      }    
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
Float_t AliPHOSDigitizer::TimeOfNoise(void)
{  // Calculates the time signal generated by noise
  //to be rewritten, now returns just big number
  return 1. ;

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
  const TClonesArray * digits = gime->Digits(GetName()) ; 
  TTree * treeD ;

 
 if(fManager) 
   treeD = fManager->GetTreeD() ;
 else {
   if (!gAlice->TreeD() ) 
     gAlice->MakeTree("D", fSplitFile);
   treeD = gAlice->TreeD();
 }


  // -- create Digits branch
  Int_t bufferSize = 32000 ;    
  TBranch * digitsBranch = treeD->Branch("PHOS",&digits,bufferSize);
  digitsBranch->SetTitle(GetName());
    
  // -- Create Digitizer branch
  Int_t splitlevel = 0 ;
  AliPHOSDigitizer * d = gime->Digitizer(GetName()) ;
  TBranch * digitizerBranch = treeD->Branch("AliPHOSDigitizer", "AliPHOSDigitizer", &d,bufferSize,splitlevel); 
  digitizerBranch->SetTitle(GetName());

  digitsBranch->Fill() ;
  digitizerBranch->Fill() ; 
  treeD->AutoSave() ; 
 
}
