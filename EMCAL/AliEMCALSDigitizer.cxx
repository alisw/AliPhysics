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
#include "TFile.h"
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFolder.h"
#include "TBenchmark.h"
#include "TGeometry.h"

// --- Standard library ---

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliHeader.h"
#include "AliEMCALDigit.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGetter.h"
#include "AliEMCALHit.h"
#include "AliEMCALSDigitizer.h"

ClassImp(AliEMCALSDigitizer)
           
//____________________________________________________________________________ 
  AliEMCALSDigitizer::AliEMCALSDigitizer():TTask("AliEMCALSDigitizer","") 
{
  // ctor
  InitParameters() ; 
  fSplitFile           = 0 ; 
  fToSplit             = kFALSE ;
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::AliEMCALSDigitizer(const char* headerFile, const char *sDigitsTitle, const Bool_t toSplit):
TTask(sDigitsTitle, headerFile)
{
  // ctor
  fToSplit = toSplit ;
  Init();
  InitParameters() ; 
  fDefaultInit = kFALSE ; 
}

//____________________________________________________________________________ 
AliEMCALSDigitizer::~AliEMCALSDigitizer()
{
  // dtor
  fSplitFile = 0 ; 
}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::Init(){
  // Initialization: open root-file, allocate arrays for hits and sdigits,
  // attach task SDigitizer to the list of EMCAL tasks
  // 
  // Initialization can not be done in the default constructor
  //============================================================= YS
  //  The initialisation is now done by the getter

  if( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
   
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(), GetName(), fToSplit) ; 
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object !" ) ;  
    return ;
  } 
  
  gime->PostSDigits( GetName(), GetTitle() ) ; 
  
  fSplitFile = 0 ;
 
  if(fToSplit){
    // construct the name of the file as /path/EMCAL.SDigits.root
    // First - extract full path if necessary
    TString sDigitsFileName(GetTitle()) ;
    Ssiz_t islash = sDigitsFileName.Last('/') ;
    if(islash<sDigitsFileName.Length())
      sDigitsFileName.Remove(islash+1,sDigitsFileName.Length()) ;
    else
      sDigitsFileName="" ;

    // Next - append the file name 
    sDigitsFileName+="EMCAL.SDigits." ;
    if((strcmp(GetName(),"Default")!=0)&&(strcmp(GetName(),"")!=0)){
      sDigitsFileName+=GetName() ;
      sDigitsFileName+="." ;
    }
    sDigitsFileName+="root" ;

    // Finally - check if the file already opened or open the file
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(sDigitsFileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(sDigitsFileName.Data(),"update") ;
  }

  TString sdname(GetName() );
  sdname.Append(":") ;
  sdname.Append(GetTitle() ) ;
  SetName(sdname) ;
  gime->PostSDigitizer(this) ;

}

//____________________________________________________________________________ 
void AliEMCALSDigitizer::InitParameters()
{
  fA                      = 0;
  fB                      = 10000000.;

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  const AliEMCALGeometry * geom = gime->EMCALGeometry() ; 
  if (geom->GetSampling() == 0.) {
    Error("InitParameters", "Sampling factor not set !") ; 
    abort() ;
  }
  else
    Info("InitParameters", "Sampling factor set to %f\n", geom->GetSampling()) ; 
  
  // this threshold corresponds approximately to 100 MeV
  fECPrimThreshold     = 100E-3 / ( geom->GetSampling() * ( geom->GetNPRLayers() + geom->GetNECLayers()) ) * geom->GetNECLayers() ;
  fPREPrimThreshold    = 100E-3 / ( geom->GetSampling() * ( geom->GetNPRLayers() + geom->GetNECLayers()) ) * geom->GetNPRLayers() ; 
  fHCPrimThreshold     = fECPrimThreshold/5. ; // 5 is totally arbitrary

}

//____________________________________________________________________________
void AliEMCALSDigitizer::Exec(Option_t *option) 
{ 
  // Collects all hits in the section (PRE/ECAL/HCAL) of the same tower into digit

  if( strcmp(GetName(), "") == 0 )
    Init() ;
  
  if (strstr(option, "print") ) {
    Print("") ; 
    return ; 
  }
  
  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALSDigitizer");

  //Check, if this branch already exits
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ;
  if(gime->BranchExists("SDigits") ) 
    return;   

  TString sdname(GetName()) ;
  sdname.Remove(sdname.Index(GetTitle())-1) ;

  Int_t nevents = gime->MaxEvent() ; 
  Int_t ievent ;   
  for(ievent = 0; ievent < nevents; ievent++){     
    gime->Event(ievent,"H") ;  
    const TClonesArray * hits = gime->Hits() ; 
    TClonesArray * sdigits = gime->SDigits(sdname.Data()) ;
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    
    //Now make SDigits from hits, for EMCAL it is the same, so just copy    
    Int_t nPrim =  static_cast<Int_t>((gAlice->TreeH())->GetEntries()) ; 
    // Attention nPrim is the number of primaries tracked by Geant 
    // and this number could be different to the number of Primaries in TreeK;
    Int_t iprim ;
    for ( iprim = 0 ; iprim < nPrim ; iprim++ ) { 
      //=========== Get the EMCAL branch from Hits Tree for the Primary iprim
      gime->Track(iprim) ;
      Int_t i;
      for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
	AliEMCALHit * hit = dynamic_cast<AliEMCALHit*>(hits->At(i)) ;
	AliEMCALDigit * curSDigit = 0 ;
	AliEMCALDigit * sdigit = 0 ;
	Bool_t newsdigit = kTRUE; 

	// Assign primary number only if deposited energy is significant

	const AliEMCALGeometry * geom = gime->EMCALGeometry() ; 

	if (gDebug) 
	  Info("Exec", "id = %d energy = %f thresholdPRE = %f thresholdEC = %f thresholdHC = %f \n", 
	       hit->GetId(), hit->GetEnergy(), fPREPrimThreshold, fECPrimThreshold, fHCPrimThreshold) ;    
	
	if( geom->IsInPRE(hit->GetId()) )  
	  if( hit->GetEnergy() > fPREPrimThreshold )
	    curSDigit =  new AliEMCALDigit( hit->GetPrimary(),
					    hit->GetIparent(), hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	  else
	    curSDigit =  new AliEMCALDigit( -1               , 
					    -1               ,
					    hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	else if( geom->IsInECAL(hit->GetId()) )
	  if( hit->GetEnergy() >  fECPrimThreshold )
	    curSDigit =  new AliEMCALDigit( hit->GetPrimary(),
					    hit->GetIparent(), hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	  else
	    curSDigit =  new AliEMCALDigit( -1               , 
					    -1               ,
					    hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	else if( geom->IsInHCAL(hit->GetId()) )
	  if( hit->GetEnergy() >  fHCPrimThreshold )
	    
	    curSDigit =  new AliEMCALDigit( hit->GetPrimary(),
					    hit->GetIparent(), hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	  else
	    curSDigit =  new AliEMCALDigit( -1               , 
					    -1               ,
					    hit->GetId(), 
					    Digitize(hit->GetEnergy()), hit->GetTime() ) ;
	
	Int_t check = 0 ;
	for(check= 0; check < nSdigits ; check++) {
	  sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(check)) ;
	  if( sdigit->GetId() == curSDigit->GetId()) { // Are we in the same ECAL/HCAL/preshower tower ?              
	    *sdigit = *sdigit + *curSDigit;
	    newsdigit = kFALSE;
	  }
	}
	if (newsdigit) { 
	  new((*sdigits)[nSdigits])  AliEMCALDigit(*curSDigit);
	  nSdigits++ ;  
	}
	delete curSDigit ; 
      }  // loop over all hits (hit = deposited energy/layer/entering particle)
    } // loop over iprim
    
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntriesFast() ;
    fSDigitsInRun += nSdigits ;  
    sdigits->Expand(nSdigits) ;

    Int_t NPrimarymax = -1 ; 
    Int_t i ;
    for (i = 0 ; i < sdigits->GetEntriesFast() ; i++) { 
      AliEMCALDigit * sdigit = dynamic_cast<AliEMCALDigit *>(sdigits->At(i)) ;
      sdigit->SetIndexInList(i) ;
    }
    
    for (i = 0 ; i < sdigits->GetEntriesFast() ; i++) {   
      if (((dynamic_cast<AliEMCALDigit *>(sdigits->At(i)))->GetNprimary()) > NPrimarymax)
	NPrimarymax = ((dynamic_cast<AliEMCALDigit *>(sdigits->At(i)))->GetNprimary()) ;
    }
    
    // Now write SDigits
    
    if(gAlice->TreeS() == 0 || (fSplitFile))  //<--- To be checked: we should not create TreeS if it is already here
      gAlice->MakeTree("S",fSplitFile);
  
    if(fSplitFile)
      fSplitFile->cd() ;
    
    //First list of sdigits
    Int_t bufferSize = 32000 ;    
    TBranch * sdigitsBranch = gAlice->TreeS()->Branch("EMCAL",&sdigits,bufferSize);
    sdigitsBranch->SetTitle(sdname);
    
    //NEXT - SDigitizer
    Int_t splitlevel = 0 ;
    AliEMCALSDigitizer * sd = this ;
    TBranch * sdigitizerBranch = gAlice->TreeS()->Branch("AliEMCALSDigitizer","AliEMCALSDigitizer",
							 &sd,bufferSize,splitlevel); 
    sdigitizerBranch->SetTitle(sdname);
    
    sdigitsBranch->Fill() ; 
    sdigitizerBranch->Fill() ; 
  
    gAlice->TreeS()->AutoSave() ;
    
    if(strstr(option,"deb"))
      PrintSDigits(option) ;  
  }
   
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALSDigitizer"); 
    Info("Exec", "took %f seconds for SDigitizing %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALSDigitizer"), gBenchmark->GetCpuTime("EMCALSDigitizer") ) ; 
  }
}

//__________________________________________________________________
void AliEMCALSDigitizer::SetSDigitsBranch(const char * title ){
 
  // Setting title to branch SDigits 

  TString stitle(title) ;

  // check if branch with title already exists
  TBranch * sdigitsBranch    = 
    static_cast<TBranch*>(gAlice->TreeS()->GetListOfBranches()->FindObject("EMCAL")) ; 
  TBranch * sdigitizerBranch =  
    static_cast<TBranch*>(gAlice->TreeS()->GetListOfBranches()->FindObject("AliEMCALSDigitizer")) ;
  const char * sdigitsTitle    = sdigitsBranch ->GetTitle() ;  
  const char * sdigitizerTitle = sdigitizerBranch ->GetTitle() ;
  if ( stitle.CompareTo(sdigitsTitle)==0 || stitle.CompareTo(sdigitizerTitle)==0 ){
    Error("SetSDigitsBranch", "Cannot overwrite existing branch with title %s", title) ;
    return ;
  }
  
  Info("SetSDigitsBranch", "Changing SDigits file from %s to %s", GetName(), title) ;

  SetName(title) ; 
    
  // Post to the WhiteBoard
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  gime->PostSDigits( title, GetTitle()) ; 
}


//__________________________________________________________________
void AliEMCALSDigitizer::Print(Option_t* option)const
{ 
  // Prints parameters of SDigitizer

  TString message("\n") ; 
  message += "------------------- "; 
  message += GetName() ; 
  message += " -------------\n" ;
  message += "   Writing SDigitis to branch with title  " ; 
  message += GetName() ;
  message += "\n   with digitization parameters  A               = " ; 
  message += fA ;
  message += "\n                                 B               = " ; 
  message += fB ; 
  message += "\n   Threshold for Primary assignment in PreShower = " ; 
  message += fPREPrimThreshold ; 
  message += "\n   Threshold for Primary assignment in EC section= " ; 
  message += fECPrimThreshold ; 
  message += "\n   Threshold for Primary assignment in HC section= " ; 
  message += fHCPrimThreshold ; 
  message += "\n---------------------------------------------------" ;
  
  Info("Print", message.Data() ) ; 
}

//__________________________________________________________________
Bool_t AliEMCALSDigitizer::operator==( AliEMCALSDigitizer const &sd )const
{
  // Equal operator.
  // SDititizers are equal if their pedestal, slope and threshold are equal

  if( (fA==sd.fA)&&(fB==sd.fB)&&
      (fECPrimThreshold==sd.fECPrimThreshold) &&
      (fHCPrimThreshold==sd.fHCPrimThreshold) &&
      (fPREPrimThreshold==sd.fPREPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}

//__________________________________________________________________
void AliEMCALSDigitizer::PrintSDigits(Option_t * option){
  //Prints list of digits produced at the current pass of AliEMCALDigitizer
  
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TString sdname(GetName()) ;
  sdname.Remove(sdname.Index(GetTitle())-1) ;
  const TClonesArray * sdigits = gime->SDigits(sdname.Data()) ; 
  
  TString message("\n") ;  
  message += "event " ; 
  message += gAlice->GetEvNumber() ;
  message += "\n      Number of entries in SDigits list " ;
  message += sdigits->GetEntriesFast() ; 
  if(strstr(option,"all")||strstr(option,"EMC")){
    
    //loop over digits
    AliEMCALDigit * digit;
    message += "\n   Id  Amplitude    Time          Index Nprim: Primaries list \n" ;    
    Int_t index ;
    char * tempo = new char[8192]; 
    for (index = 0 ; index < sdigits->GetEntries() ; index++) {
      digit = dynamic_cast<AliEMCALDigit *>( sdigits->At(index) ) ;
      sprintf(tempo, "\n%6d  %8d    %6.5e %4d      %2d :",
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
  Info("PrintSDigits", message.Data() ) ; 
}

//_______________________________________________________________________________________
// void AliEMCALSDigitizer::TestTowerID(void)
// {
//   Int_t j;

//   Bool_t preshower = kFALSE;
//   for (j = 0 ; j < 10 ; j++){  // loop over hit id
//     Int_t i;
//    for (i = 0 ; i <= 2 ; i++){  // loop over 
//      Int_t k = i*96*144+j*144+1;
//       Info("TestTowerID", " Hit Index = %d  %d   TOWERID = %d", k, j*10, Layer2TowerID(k, preshower) ) ;
//     }
//   }
// }

//____________________________________________________________________________ 
void AliEMCALSDigitizer::UseHitsFrom(const char * filename)
{
  SetTitle(filename) ; 
  Init() ; 
}
