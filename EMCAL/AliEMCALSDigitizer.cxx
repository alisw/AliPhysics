/*************************************************************************
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
// This is a TTask that makes SDigits out of Hits
// A Summable Digits is the sum of all hits originating 
// from one in one tower of the EMCAL 
// A threshold for assignment of the primary to SDigit is applied 
// SDigits are written to TreeS, branch "EMCAL"
// AliEMCALSDigitizer with all current parameters is written 
// to TreeS branch "AliEMCALSDigitizer".
// Both branches have the same title. If necessary one can produce 
// another set of SDigits with different parameters. Two versions
// can be distunguished using titles of the branches.
// User case:
// root [0] AliEMCALSDigitizer * s = new AliEMCALSDigitizer("galice.root")
// Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
// root [1] s->ExecuteTask()
//             // Makes SDigitis for all events stored in galice.root
// root [2] s->SetPedestalParameter(0.001)
//             // One can change parameters of digitization
// root [3] s->SetSDigitsBranch("Redestal 0.001")
//             // and write them into the new branch
// root [4] s->ExeciteTask("deb all tim")
//             // available parameters:
//             deb - print # of produced SDigitis
//             deb all  - print # and list of produced SDigits
//             tim - print benchmarking information
//
//*-- Author : Sahal Yacoob (LBL)
// based on  : AliPHOSSDigitzer 
// Modif: 
//  August 2002 Yves Schutz: clone PHOS as closely as possible and intoduction
//                           of new  IO (à la PHOS)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include <TBenchmark.h>
#include <TH1.h>
#include <TBrowser.h>
#include <Riostream.h>
#include <TMath.h>

// --- Standard library ---
#include "stdlib.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliEMCALDigit.h"
#include "AliEMCALLoader.h"
#include "AliEMCALHit.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHistoUtilities.h"

ClassImp(AliEMCALSDigitizer)
           
//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer()
  : TTask("",""),
    fA(0.),fB(0.),fECPrimThreshold(0.),
    fDefaultInit(kTRUE),
    fEventFolderName(0),
    fInit(0),
    fSDigitsInRun(0),
    fFirstEvent(0),
    fLastEvent(0),
    fSampling(0.),
    fControlHists(0),
    fHists(0)
{
  // ctor
  InitParameters();
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const char * alirunFileName, 
				       const char * eventFolderName)
  : TTask("EMCAL"+AliConfig::Instance()->GetSDigitizerTaskName(), alirunFileName),
    fA(0.),fB(0.),fECPrimThreshold(0.),
    fDefaultInit(kFALSE),
    fEventFolderName(eventFolderName),
    fInit(0),
    fSDigitsInRun(0),
    fFirstEvent(0),
    fLastEvent(0),
    fSampling(0.),
    fControlHists(1),
    fHists(0)
{
  // ctor
  Init();
  InitParameters() ; 
  if(fControlHists) BookControlHists(1);
}


//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const AliEMCALSDigitizer & sd) 
  : TTask(sd.GetName(),sd.GetTitle()),
    fA(sd.fA),
    fB(sd.fB),
    fECPrimThreshold(sd.fECPrimThreshold),
    fDefaultInit(sd.fDefaultInit),
    fEventFolderName(sd.fEventFolderName),
    fInit(sd.fInit),
    fSDigitsInRun(sd.fSDigitsInRun),
    fFirstEvent(sd.fFirstEvent),
    fLastEvent(sd.fLastEvent),
    fSampling(sd.fSampling),
    fControlHists(sd.fControlHists),
    fHists(sd.fHists)
{
  //cpy ctor 
}


//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer() {
  //dtor
  AliLoader *emcalLoader=0;
  if ((emcalLoader = AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL")))
      emcalLoader->CleanSDigitizer();
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init(){
  // Initialization: open root-file, allocate arrays for hits and sdigits,
  // attach task SDigitizer to the list of EMCAL tasks
  // 
  // Initialization can not be done in the default constructor
  //============================================================= YS
  //  The initialisation is now done by the getter

  fInit = kTRUE; 
   
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));

  if ( emcalLoader == 0 ) {
    Fatal("Init", "Could not obtain the AliEMCALLoader");
    return ;
  } 
  
  emcalLoader->PostSDigitizer(this);
  emcalLoader->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::InitParameters()
{
  //initialize parameters for sdigitization

  const AliEMCALGeometry * geom = AliEMCALGeometry::GetInstance();
  if (geom->GetSampling() == 0.) {
    Fatal("InitParameters", "Sampling factor not set !") ; 
  }
  // Initializes parameters
  fA         = 0;
  fB         = 1.e+7;
  fSampling  = geom->GetSampling();

 // threshold for deposit energy of hit
  fECPrimThreshold  = 0.; // 24-nov-04 - was 1.e-6;
  Print("");
}

//____________________________________________________________________________
void AliEMCALSDigitizer::Exec(Option_t *option) 
{ 
  // Collects all hit of the same tower into digits
  TString o(option); o.ToUpper();
  if (strstr(option, "print") ) {
    Print() ; 
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALSDigitizer");

  AliRunLoader *rl = AliRunLoader::GetRunLoader();
  AliEMCALLoader *emcalLoader = dynamic_cast<AliEMCALLoader*>(rl->GetDetectorLoader("EMCAL"));

  //switch off reloading of this task while getting event
  if (!fInit) { // to prevent overwrite existing file
    AliError( Form("Give a version name different from %s", fEventFolderName.Data()) ) ;
    return ;
    }

  if (fLastEvent == -1) 
    fLastEvent = rl->GetNumberOfEvents() - 1 ;
  else {
    fLastEvent = TMath::Min(fLastEvent, rl->GetNumberOfEvents()-1);
  }
  Int_t nEvents   = fLastEvent - fFirstEvent + 1;

  Int_t ievent;
  Float_t energy=0.; // de * fSampling - 23-nov-04
  rl->LoadKinematics();
  rl->LoadHits("EMCAL");

  for (ievent = fFirstEvent; ievent <= fLastEvent; ievent++) {
    rl->GetEvent(ievent);
    TTree * treeS = emcalLoader->TreeS();
    if ( !treeS ) { 
      emcalLoader->MakeSDigitsContainer();
      treeS = emcalLoader->TreeS();
    }
    TClonesArray * hits = emcalLoader->Hits() ; 
    TClonesArray * sdigits = emcalLoader->SDigits() ;
    sdigits->Clear();

    Int_t nSdigits = 0 ;
    Int_t i;
    AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance(); 
    for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
      AliEMCALHit * hit = dynamic_cast<AliEMCALHit*>(hits->At(i)) ;
      AliEMCALDigit * curSDigit = 0 ;
      AliEMCALDigit * sdigit = 0 ;
      Bool_t newsdigit = kTRUE; 
      
      // hit->GetId() - Absolute Id number EMCAL segment
      if(geom->CheckAbsCellId(hit->GetId())) { // was IsInECA(hit->GetId())
	energy = hit->GetEnergy() * fSampling; // 23-nov-04
	if(energy >  fECPrimThreshold )
	  // Assign primary number only if deposited energy is significant
	  curSDigit =  new AliEMCALDigit( hit->GetPrimary(),
					  hit->GetIparent(), hit->GetId(), 
					  Digitize(energy), hit->GetTime() ) ;
	  else
	    curSDigit =  new AliEMCALDigit( -1               , 
					    -1               ,
					    hit->GetId(), 
					    Digitize(energy), hit->GetTime() ) ;
      } else {
	Warning("Exec"," abs id %i is bad \n", hit->GetId());
	newsdigit = kFALSE;
	curSDigit = 0;
      }
      
      if(curSDigit != 0){
	for(Int_t check= 0; check < nSdigits ; check++) {
	  sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(check)) ;
	  
	  if( sdigit->GetId() == curSDigit->GetId()) { // Are we in the same ECAL tower ?              
	    *sdigit = *sdigit + *curSDigit;
	    newsdigit = kFALSE;
	  }
	}
      }
      if (newsdigit) {
	new((*sdigits)[nSdigits])  AliEMCALDigit(*curSDigit);
	nSdigits++ ;  
      }
      delete curSDigit ; 
    }  // loop over all hits (hit = deposited energy/entering particle)
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntriesFast() ;
    fSDigitsInRun += nSdigits ;  
    
    Double_t e=0.,esum=0.;
    AliEMCALHistoUtilities::FillH1(fHists, 0, double(sdigits->GetEntriesFast()));
    for (i = 0 ; i < sdigits->GetEntriesFast() ; i++) { 
      AliEMCALDigit * sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(i)) ;
      sdigit->SetIndexInList(i) ;
      
      AliEMCALHistoUtilities::FillH1(fHists, 2, double(sdigit->GetAmp()));
      e = double(Calibrate(sdigit->GetAmp()));
      esum += e;
      AliEMCALHistoUtilities::FillH1(fHists, 3, e);
      AliEMCALHistoUtilities::FillH1(fHists, 4, double(sdigit->GetId()));
    }
    if(esum>0.) AliEMCALHistoUtilities::FillH1(fHists, 1, esum);
    
    // Now write SDigits    
    
    Int_t bufferSize = 32000 ;    
    TBranch * sdigitsBranch = treeS->GetBranch("EMCAL");
    if (sdigitsBranch)
      sdigitsBranch->SetAddress(&sdigits);
    else
      treeS->Branch("EMCAL",&sdigits,bufferSize);
    
    treeS->Fill();
    
    emcalLoader->WriteSDigits("OVERWRITE");
    
    //NEXT - SDigitizer
    emcalLoader->WriteSDigitizer("OVERWRITE");  // why in event cycle ?
    
    if(strstr(option,"deb"))
      PrintSDigits(option) ;  
  }
  
  Unload();
  
  emcalLoader->GetSDigitsDataLoader()->GetBaseTaskLoader()->SetDoNotReload(kTRUE);
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALSDigitizer"); 
    printf("\n Exec: took %f seconds for SDigitizing %f seconds per event\n", 
	   gBenchmark->GetCpuTime("EMCALSDigitizer"), gBenchmark->GetCpuTime("EMCALSDigitizer")/nEvents ) ; 
  }
}

//__________________________________________________________________

Int_t AliEMCALSDigitizer::Digitize(Float_t energy)const {
  // Digitize the energy
    Double_t aSignal = fA + energy*fB;
    if (TMath::Abs(aSignal)>2147483647.0) { 
      //PH 2147483647 is the max. integer
      //PH This apparently is a problem which needs investigation
      AliWarning(Form("Too big or too small energy %f",aSignal));
      aSignal = TMath::Sign((Double_t)2147483647,aSignal);
    }
    return (Int_t ) aSignal;
  }
 

//__________________________________________________________________

void AliEMCALSDigitizer::Print1(Option_t * option)
{
  Print(); 
  PrintSDigits(option);
}

//__________________________________________________________________
void AliEMCALSDigitizer::Print(Option_t *option) const
{ 
  // Prints parameters of SDigitizer
  printf("Print: \n------------------- %s ------------- option %s\n", GetName() , option) ; 
  printf("   fInit                                 %i\n", int(fInit));
  printf("   fFirstEvent                           %i\n", fFirstEvent);
  printf("   fLastEvent                            %i\n", fLastEvent);
  printf("   Writing SDigits to branch with title  %s\n", fEventFolderName.Data()) ;
  printf("   with digitization parameters       A = %f\n", fA) ; 
  printf("                                      B = %f\n", fB) ;
  printf("   Threshold for EC Primary assignment  = %f\n", fECPrimThreshold)  ;
  printf("   Sampling                             = %f\n", fSampling);
  printf("---------------------------------------------------\n") ;
}

//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal
  if( (fA==sd.fA)&&(fB==sd.fB)&&
      (fECPrimThreshold==sd.fECPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option)
{
  //Prints list of digits produced at the current pass of AliEMCALDigitizer
  
  
  // AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
  const TClonesArray * sdigits = rl->SDigits() ; 
  
  printf("\n") ;  
  printf("event %i", rl->GetRunLoader()->GetEventNumber());
  printf(" Number of entries in SDigits list %i", sdigits->GetEntriesFast()); 
  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliEMCALDigit * digit;
    printf("\n   Id  Amplitude    Time          Index Nprim: Primaries list \n") ;    
    Int_t index, isum=0;
    char * tempo = new char[8192]; 
    for (index = 0 ; index < sdigits->GetEntries() ; index++) {
      digit = dynamic_cast<AliEMCALDigit *>( sdigits->At(index) ) ;
      sprintf(tempo, "\n%6d  %8d    %6.5e %4d      %2d :",
	      digit->GetId(), digit->GetAmp(), digit->GetTime(), digit->GetIndexInList(), digit->GetNprimary()) ;  
      printf(tempo);
      isum += digit->GetAmp();
      
      Int_t iprimary;
      for (iprimary=0; iprimary<digit->GetNprimary(); iprimary++) {
	sprintf(tempo, "%d ",digit->GetPrimary(iprimary+1) ) ; 
	printf(tempo); 
      }  	 
    }
    delete tempo ;
    printf("\n** Sum %i : %10.3f GeV/c **\n ", isum, double(isum)*1.e-6);
  } else printf("\n");
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Unload() const
{
  // Unload Hits and SDigits from the folder
  AliEMCALLoader *rl = dynamic_cast<AliEMCALLoader*>(AliRunLoader::GetRunLoader()->GetDetectorLoader("EMCAL"));
  rl->UnloadHits() ; 
  rl->UnloadSDigits() ; 
}

void AliEMCALSDigitizer::Browse(TBrowser* b)
{
  if(fHists) b->Add(fHists);
  TTask::Browse(b);
}

TList *AliEMCALSDigitizer::BookControlHists(int var)
{ 
  //book histograms for monitoring sdigitization
  // 22-nov-04
  gROOT->cd();
  const AliEMCALGeometry *geom = AliEMCALGeometry::GetInstance() ;
  if(var>=1){
    printf(" AliEMCALSDigitizer::BookControlHists() in action \n");
    new TH1F("HSDigiN",  "#EMCAL  sdigits ", 1001, -0.5, 1000.5);
    new TH1F("HSDigiSumEnergy","Sum.EMCAL energy", 1000, 0.0, 100.);
    new TH1F("HSDigiAmp",  "EMCAL sdigits amplitude", 1000, 0., 2.e+9);
    new TH1F("HSDigiEnergy","EMCAL cell energy", 1000, 0.0, 100.);
    new TH1F("HSDigiAbsId","EMCAL absID for sdigits",
    geom->GetNCells(), 0.5, Double_t(geom->GetNCells())+0.5);
  }

  fHists = AliEMCALHistoUtilities::MoveHistsToList("EmcalSDigiControlHists", kFALSE);
  fHists = 0;

  return fHists;
}

void AliEMCALSDigitizer::SaveHists(const char* name, Bool_t kSingleKey, const char* opt)
{
  AliEMCALHistoUtilities::SaveListOfHists(fHists, name, kSingleKey, opt); 
}
