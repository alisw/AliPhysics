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


//_________________________________________________________________________
// This is a TTask that makes TOF-Digits out of TOF-Hits.
// A TOF Digit is essentially a TOF hit with the smearing for strip 
// time resolution and simulation of the ADC correlation signal. 
// Digits are written to TreeD in branch "TOF".
// AliTOFDigitizer with all current parameters used for digitization 
// (time resolution and ADC parameter) is written 
// to TreeD branch "AliTOFDigitizer".
// Both branches have the same title. If necessary one can produce 
// another set of Digits with different parameters. Two versions
// can be distunguished using titles of the branches.
// User case:
//  root [0] AliTOFDigitizer * s = new AliTOFDigitizer("galice.root")
//  Warning in <AliITSgeomSPD425Long::Default Creator>: Detector size may not be write.
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] s->ExecuteTask()
//             // Makes Digits for all events stored in galice.root
//  root [2] s->SetTimeRes(100.)
//             // One can change parameters of digitization
//  root [3] s->SetDigitsBranch("Time Resolution 100. ps")
//             // and write them into the new branch
//  root [4] s->ExecuteTask("deb all tim")
//             // available parameters:
//             deb - print number of produced Digits
//             deb all  - print number and list of produced Digits
//             tim - print benchmarking information
//
// -- Author :  F. Pierella (Bologna University) pierella@bo.infn.it
//////////////////////////////////////////////////////////////////////////////

#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TBenchmark.h"
#include "TRandom.h"

#include <iomanip.h>

#include "AliRun.h"
#include "AliTOFdigit.h"
#include "AliTOFhit.h"
#include "AliTOFDigitizer.h"


ClassImp(AliTOFDigitizer)

           
//____________________________________________________________________________ 
  AliTOFDigitizer::AliTOFDigitizer():TTask("AliTOFDigitizer","") 
{
  // ctor
  fTimeRes       = 100;  // ps
  fChrgRes       = 100.; // pC
  fNevents       = 0 ;      
  fDigits        = 0 ;
  fHits          = 0 ;
  fIsInitialized = kFALSE ;

}

//____________________________________________________________________________ 
AliTOFDigitizer::AliTOFDigitizer(const char* headerFile, const char *digitsTitle):TTask("AliTOFDigitizer","")
{
  // ctor
  fTimeRes       = 100;  // ps
  fChrgRes       = 100.; // pC
  fNevents       = 0 ;      
  fDigitsTitle   = digitsTitle ;
  fHeadersFile   = headerFile ;
  fIsInitialized = kFALSE ;
  Init();
}

//____________________________________________________________________________ 
AliTOFDigitizer::~AliTOFDigitizer()
{
  // dtor
  if(fDigits)
    delete fDigits ;
  if(fHits)
    delete fHits ;
}
//____________________________________________________________________________ 
void AliTOFDigitizer::Init()
{
  // Initialization: open root-file, allocate arrays for hits and digits,
  // attach task Digitizer to the list of TOF tasks
  // 
  // Initialization can not be done in the default constructor

  if(!fIsInitialized){

    if(fHeadersFile.IsNull())
      fHeadersFile="galice.root" ;

    TFile * file = (TFile*) gROOT->GetFile(fHeadersFile.Data() ) ;
    
    //if file was not opened yet, read gAlice
    if(file == 0){
      if(fHeadersFile.Contains("rfio"))
	file =	TFile::Open(fHeadersFile,"update") ;
      else
	file = new TFile(fHeadersFile.Data(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }
    fHits    = new TClonesArray("AliTOFhit"  , 405);
    fDigits  = new TClonesArray("AliTOFdigit", 405);
    
    //add Task to //FPAlice/tasks/(S)Digitizer/TOF
    TFolder * alice  = (TFolder*)gROOT->GetListOfBrowsables()->FindObject("FPAlice") ;
    TTask * aliceSD  = (TTask*)alice->FindObject("tasks/(S)Digitizer") ;
    TTask * tofD   = (TTask*)aliceSD->GetListOfTasks()->FindObject("TOF") ;
    tofD->Add(this) ;

    fIsInitialized = kTRUE ;
  }
}
//____________________________________________________________________________
void AliTOFDigitizer::Exec(Option_t *option) 
{ 
  Int_t    tracks[3];    // track info
  Int_t    vol[5];       // dummy location for digit
  Float_t  digit[2];     // TOF digit variables

  TRandom* rnd = new TRandom();

  // Collects all hits in the same active volume into digit
  
  if(!fIsInitialized)
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("TOFDigitizer");
  
  fNevents = (Int_t) gAlice->TreeE()->GetEntries() ; 
  
  Int_t ievent ;
  // start loop on events
  for(ievent = 0; ievent < fNevents; ievent++){
    gAlice->GetEvent(ievent) ;
    gAlice->SetEvent(ievent) ;
    
    if(gAlice->TreeH()==0){
      cout << "AliTOFDigitizer: There is no Hit Tree" << endl;
      return ;
    }
    
    //set address of the hits 
    TBranch * branch = gAlice->TreeH()->GetBranch("TOF");
    if (branch) 
      branch->SetAddress(&fHits);
    else{
      cout << "ERROR in AliTOFDigitizer: "<< endl ;
      cout << "      no branch TOF in TreeH"<< endl ;
      cout << "      do nothing " << endl ;
      return ;
    }
    
    fDigits->Clear();
    Int_t ndigits = 0 ;
    
    //Now made Digits from hits, for TOF it is the same a part for the tof smearing
    // and some missing MC variables     
    Int_t itrack ;
    for (itrack=0; itrack < gAlice->GetNtrack(); itrack++){

      //=========== Get the TOF branch from Hits Tree for the Primary track itrack
      branch->GetEntry(itrack,0);
      
      Int_t i;
      for ( i = 0 ; i < fHits->GetEntries() ; i++ ) {

	AliTOFhit * hit = (AliTOFhit*)fHits->At(i) ;

        vol[0] = hit->GetSector();
        vol[1] = hit->GetPlate();
        vol[2] = hit->GetPadx();
        vol[3] = hit->GetPadz();
        vol[4] = hit->GetStrip();
        // 95% of efficiency to be inserted here
        // edge effect to be inserted here
        // cross talk  to be inserted here
        // simulation of the detector response
        Float_t idealtime = hit->GetTof(); // unit s
        idealtime *= 1.E+12;  // conversion from s to ps  
                              // fTimeRes is given usually in ps
        Float_t tdctime   = rnd->Gaus(idealtime, fTimeRes);
        digit[0] = tdctime;
        // typical Landau Distribution to be inserted here
        Float_t idealcharge = hit->GetEdep();
        Float_t adccharge = rnd->Gaus(idealcharge, fChrgRes);
        digit[1] = adccharge;
        // to be added a check for overlaps
        new((*fDigits)[ndigits]) AliTOFdigit(tracks, vol, digit);
	ndigits++ ;  
      } 
      
    } // end loop over tracks

    delete rnd;
    rnd = 0;
    
    ndigits = fDigits->GetEntriesFast() ;
    printf("AliTOFDigitizer: Total number of digits %d\n",ndigits);

    if(gAlice->TreeD() == 0)
      gAlice->MakeTree("D") ;
    
    //check, if this branch already exits
    TBranch * digitsBranch = 0;
    TBranch * digitizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t tofNotFound = kTRUE ;
    Bool_t digitizerNotFound = kTRUE ;
    
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
      
      if(tofNotFound){
	digitsBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp("TOF",digitsBranch->GetName())==0 ) &&
	    (fDigitsTitle.CompareTo(digitsBranch->GetTitle()) == 0) )
	  tofNotFound = kFALSE ;
      }
      if(digitizerNotFound){
	digitizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(digitizerBranch->GetName(),"AliTOFDigitizer") == 0)&&
	    (fDigitsTitle.CompareTo(digitizerBranch->GetTitle()) == 0) )
	  digitizerNotFound = kFALSE ;
      }
    }

    if(!(digitizerNotFound && tofNotFound)){
      cout << "AliTOFdigitizer error:" << endl ;
      cout << "Can not overwrite existing branches: do not write" << endl ;
      return ;
    }
    
    //Make (if necessary) branches    
    char * file =0;
    if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
      file = new char[strlen(gAlice->GetBaseFile())+20] ;
      sprintf(file,"%s/TOF.Digits.root",gAlice->GetBaseFile()) ;
    }
    
    TDirectory *cwd = gDirectory;
    
    //First list of digits
    Int_t bufferSize = 32000 ;    
    digitsBranch = gAlice->TreeD()->Branch("TOF",&fDigits,bufferSize);
    digitsBranch->SetTitle(fDigitsTitle.Data());
    if (file) {
      digitsBranch->SetFile(file);
      TIter next( digitsBranch->GetListOfBranches());
      TBranch * subbr;
      while ((subbr=(TBranch*)next())) {
	subbr->SetFile(file);
      }   
      cwd->cd();
    } 
      
    //second - Digitizer
    Int_t splitlevel = 0 ;
    AliTOFDigitizer * digtz = this ;
    digitizerBranch = gAlice->TreeD()->Branch("AliTOFDigitizer","AliTOFDigitizer",
					       &digtz,bufferSize,splitlevel); 
    digitizerBranch->SetTitle(fDigitsTitle.Data());
    if (file) {
      digitizerBranch->SetFile(file);
      TIter next( digitizerBranch->GetListOfBranches());
      TBranch * subbr ;
      while ((subbr=(TBranch*)next())) {
	subbr->SetFile(file);
      }   
      cwd->cd();
      delete file;
    }

    digitsBranch->Fill();
    digitizerBranch->Fill();
    gAlice->TreeD()->Write(0,TObject::kOverwrite) ;
    
    if(strstr(option,"deb"))
      PrintDigits(option) ;
    
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("TOFDigitizer");
    cout << "AliTOFDigitizer:" << endl ;
    cout << "   took " << gBenchmark->GetCpuTime("TOFDigitizer") << " seconds for Digitizing " 
	 <<  gBenchmark->GetCpuTime("TOFDigitizer")/fNevents << " seconds per event " << endl ;
    cout << endl ;
  }
  
  
}
//__________________________________________________________________
void AliTOFDigitizer::SetDigitsBranch(const char* title )
{
  // Setting title to branch Digits 
  if(!fDigitsTitle.IsNull())
    cout << "AliTOFdigitizer: changing Digits file from " <<fDigitsTitle.Data() << " to " << title << endl ;
  fDigitsTitle=title ;
}
//__________________________________________________________________
void AliTOFDigitizer::Print(Option_t* option)const
{
  // Prints parameters of Digitizer
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  cout << "   Writing Digits to branch with title  " << fDigitsTitle.Data() << endl ;
  cout << "   with digitization parameters Time resolution        = " << fTimeRes << endl ;
  cout << "                                Adc smearing parameter = " << fChrgRes << endl ;
  cout << "---------------------------------------------------"<<endl ;
  
}
//__________________________________________________________________
Bool_t AliTOFDigitizer::operator==( AliTOFDigitizer const &digtz )const
{
  // Equal operator.
  // Dititizers are equal if their time resolution and Adc 
  // smearing parameter are equal

  if( (fTimeRes==digtz.fTimeRes)&&(fChrgRes==digtz.fChrgRes))
    return kTRUE ;
  else
    return kFALSE ;
}
//__________________________________________________________________
void AliTOFDigitizer::PrintDigits(Option_t * option)
{
  // Prints list of digits produced in the current pass of AliTOFDigitizer
  
  cout << "AliTOFDigitizer: " << endl ;
  cout << "       Number of entries in Digits list  " << fDigits->GetEntriesFast() << endl ;
  cout << endl ;
  
  if(strstr(option,"all")){// print all digits
    
    //loop over digits
    AliTOFdigit * digit;
    cout << "Digit # " << " Time of flight(ps) " << 
            "  ADC(pC)"  <<  " Sector #" << " Plate #" << 
            " Strip #" << " Pad # (x) " << " Pad # (z) " << endl;    
    Int_t index ;
    for (index = 0 ; index < fDigits->GetEntries() ; index++) {
      digit = (AliTOFdigit * )  fDigits->At(index) ;

      // set the width of the output with the setw(int width) manipulator

      cout << setw(7)  <<  index << " " 
           << setw(13) <<  digit->GetTdc()    << "  "  
	   << setw(13) <<  digit->GetAdc()    << "  "   
	   << setw(6)  <<  digit->GetSector() << "  "   
	   << setw(6)  <<  digit->GetPlate()  << "  "   
	   << setw(6)  <<  digit->GetStrip()  << "  "   
	   << setw(7)  <<  digit->GetPadx()   << "  "   
	   << setw(7)  <<  digit->GetPadz()   << "  " << endl;   
      
    } // end loop on digits
    
  } // close if "all" option selected
}
