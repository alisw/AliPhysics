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
//*-- Author :  Dmitri Peressounko (SUBATECH & KI) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSHit.h"
#include "AliPHOSv1.h"
#include "AliPHOSSDigitizer.h"

#include "TROOT.h"
#include "TFolder.h"

ClassImp(AliPHOSSDigitizer)

           
//____________________________________________________________________________ 
  AliPHOSSDigitizer::AliPHOSSDigitizer():TTask("AliPHOSSDigitizer","") 
{
  // ctor
  fA = 0;
  fB = 10000000. ;
  fPrimThreshold = 0.01 ;
  fNevents = 0 ;     // Number of events to digitize, 0 means all evens in current file
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 

}
//____________________________________________________________________________ 
AliPHOSSDigitizer::AliPHOSSDigitizer(char* HeaderFile, char *SDigitsFile):TTask("AliPHOSSDigitizer","")
{
  // ctor
  fA = 0;
  fB = 10000000.;
  fPrimThreshold = 0.01 ;
  fNevents = 0 ;         // Number of events to digitize, 0 means all events in current file
  fSDigitsFile = SDigitsFile ;
  fHeadersFile = HeaderFile ;
  //add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
    
}

//____________________________________________________________________________ 
  AliPHOSSDigitizer::~AliPHOSSDigitizer()
{
  // dtor
}


//____________________________________________________________________________
void AliPHOSSDigitizer::Exec(Option_t *option) { 
  //Collects all hits in the same active volume into digit

  TFile * file = 0;

  // if(gAlice->TreeE()==0)        //If gAlice not yet read/constructed
  if(fHeadersFile.IsNull())
    file = new TFile("galice.root", "update") ;
  else
    file = new TFile(fHeadersFile.Data(),"update") ;

  gAlice = (AliRun *) file->Get("gAlice") ;
  
  
  
  TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1000) ;

  AliPHOS * PHOS = (AliPHOS *) gAlice->GetDetector("PHOS") ;
  
    
  if(fNevents == 0) 
    fNevents = (Int_t) gAlice->TreeE()->GetEntries() ; 
  
  Int_t ievent ;
  for(ievent = 0; ievent < fNevents; ievent++){
    gAlice->GetEvent(ievent) ;
    gAlice->SetEvent(ievent) ;
    
    if(gAlice->TreeH()==0){
      cout << "AliPHOSSDigitizer: There is no Hit Tree" << endl;
      return ;
    }
    
    if(gAlice->TreeS() == 0)
      gAlice->MakeTree("S") ;
    
    TClonesArray * hits = PHOS->Hits() ;
    
    sdigits->Clear();
    Int_t nSdigits = 0 ;
    
    //Make branches
    char branchname[20];
    sprintf(branchname,"%s",PHOS->GetName());  
    
    Int_t bufferSize = 16000 ;
    char * file =0;
    if(!fSDigitsFile.IsNull())
      file = (char*) fSDigitsFile.Data() ; //ievent ;
    else
      if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
	file = new char[30] ;
	//	sprintf(file,"PHOS.SDigits%d.root",ievent) ;
	sprintf(file,"PHOS.SDigits.root") ;
      }
      else
	file = 0 ;
    
    gAlice->MakeBranchInTree(gAlice->TreeS(),branchname,&sdigits,bufferSize,file);  


    Int_t splitlevel = 0 ;
    sprintf(branchname,"AliPHOSSDigitizer");   
    AliPHOSSDigitizer * sd = this ;
    gAlice->MakeBranchInTree(gAlice->TreeS(),branchname,"AliPHOSSDigitizer",&sd, bufferSize, splitlevel,file); 
    

    //Now made SDigits from hits, for PHOS it is the same

    Int_t itrack ;
    for (itrack=0; itrack<gAlice->GetNtrack(); itrack++){
      
      //=========== Get the Hits Tree for the Primary track itrack
      gAlice->ResetHits();    
      gAlice->TreeH()->GetEvent(itrack);
      
      Int_t i;
      for ( i = 0 ; i < hits->GetEntries() ; i++ ) {
	AliPHOSHit * hit = (AliPHOSHit*)hits->At(i) ;
	AliPHOSDigit * newdigit ;

	// Assign primary number only if contribution is significant
	if( hit->GetEnergy() > fPrimThreshold)
	  newdigit = new AliPHOSDigit( hit->GetPrimary(), hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
	else
	  newdigit = new AliPHOSDigit( -1               , hit->GetId(), Digitize( hit->GetEnergy() ) ) ;
	
	new((*sdigits)[nSdigits]) AliPHOSDigit(* newdigit) ;
	nSdigits++ ;  
	
	delete newdigit ;    
      } 
      
    } // loop over tracks
    
    sdigits->Sort() ;
    
    nSdigits = sdigits->GetEntries() ;
    sdigits->Expand(nSdigits) ;
    
    Int_t i ;
    for (i = 0 ; i < nSdigits ; i++) { 
      AliPHOSDigit * digit = (AliPHOSDigit *) sdigits->At(i) ; 
      digit->SetIndexInList(i) ;     
    }
    
    gAlice->TreeS()->Fill() ;
    gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
  }

  delete sdigits ;
  if(file)
    file->Close() ;

}
//__________________________________________________________________
void AliPHOSSDigitizer::SetSDigitsFile(char * file ){
  if(!fSDigitsFile.IsNull())
    cout << "Changing SDigits file from " <<(char *)fSDigitsFile.Data() << " to " << file << endl ;
  fSDigitsFile=file ;
}
//__________________________________________________________________
void AliPHOSSDigitizer::Print(Option_t* option)const{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  if(fSDigitsFile.IsNull())
    cout << " Writing SDigitis to file galice.root "<< endl ;
  else
    cout << "    Writing SDigitis to file  " << (char*) fSDigitsFile.Data() << endl ;
  cout << "   with digitization parameters A = " << fA << endl ;
  cout << "                                B = " << fB << endl ;
  cout << "Threshold for Primary assignment  = " << fPrimThreshold << endl ; 
  cout << "---------------------------------------------------"<<endl ;

}
//__________________________________________________________________
Bool_t AliPHOSSDigitizer::operator==( AliPHOSSDigitizer const &sd )const{
  if( (fA==sd.fA)&&(fB==sd.fB)&&(fPrimThreshold==sd.fPrimThreshold))
    return kTRUE ;
  else
    return kFALSE ;
}
