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
// This is a TTask that constructs SDigits out of Hits
// A Summable Digits is the sum of all hits in a cell
// A threshold is applied 
//
//-- Author: Alla Maevskaia(INR)
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TTask.h"
#include "TTree.h"
#include "TSystem.h"
#include "TFile.h"
// --- Standard library ---

// --- AliRoot header files ---

#include "AliFMDdigit.h"
#include "AliFMDhit.h"
#include "AliFMD.h"
#include "AliFMDv1.h"
#include "AliFMDSDigitizer.h"
#include "AliRun.h"
#include "AliDetector.h"
#include "AliMC.h"

#include "TROOT.h"
#include "TFolder.h"
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

ClassImp(AliFMDSDigitizer)

           
//____________________________________________________________________________ 
  AliFMDSDigitizer::AliFMDSDigitizer():TTask("AliFMDSDigitizer","") 
{
  fNevents = 0 ;     // Number of events to digitize, 0 means all evens in current file
  // add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
}
//____________________________________________________________________________ 
AliFMDSDigitizer::AliFMDSDigitizer(char* HeaderFile, char *SDigitsFile):TTask("AliFMDSDigitizer","")
{
  // ctor
  fNevents = 0 ;    // Number of events to digitize, 0 means all events in current file
  fSDigitsFile = SDigitsFile ;
  fHeadersFile = HeaderFile ;
  //add Task to //root/Tasks folder
  TTask * roottasks = (TTask*)gROOT->GetRootFolder()->FindObject("Tasks") ; 
  roottasks->Add(this) ; 
    
}

//____________________________________________________________________________ 
  AliFMDSDigitizer::~AliFMDSDigitizer()
{
  // dtor
}


//____________________________________________________________________________
void AliFMDSDigitizer::Exec(Option_t *option) { 
  //Collects all hits in the same active volume into digit
  TClonesArray * sdigits = new TClonesArray("AliFMDdigit",1000) ;
  TFile * file = 0;

  AliFMD * FMD = (AliFMD *) gAlice->GetDetector("FMD") ;
  
  if(fNevents == 0) 
    fNevents = (Int_t) gAlice->TreeE()->GetEntries() ; 

  cout<<"AliFMDSDigitizer-> Nevents"<<fNevents<<endl;
  for(Int_t ievent = 0; ievent < fNevents; ievent++){
    gAlice->GetEvent(ievent) ;
    if(gAlice->TreeH()==0) return ;

    if(gAlice->TreeS() == 0)    	
      gAlice->MakeTree("S") ;
    
    TClonesArray * FMDhits = FMD->Hits() ;

    

    Int_t nSdigits = 0 ;
    
    //Make branches
    char branchname[20];
    sprintf(branchname,"%s",FMD->GetName());  
    
    Int_t bufferSize = 16000 ;
    char * file =0;
    if(!fSDigitsFile.IsNull())
      file = (char*) fSDigitsFile.Data() ; //ievent ;
    else
      if(gSystem->Getenv("CONFIG_SPLIT_FILE")){ //generating file name
	file = new char[30] ;
	sprintf(file,"FMD.SDigits.root") ;
	cout<<"CONFIG_SPLIT_FILE "<<file<<endl; 
      }
      else{
	file = 0 ;
	cout<<" FILE =0 "<<endl;}
	cout<<"After CONFIG_SPLIT_FILE "<<file<<endl; 
    
    gAlice->MakeBranchInTree(gAlice->TreeS(),branchname,&sdigits,bufferSize,file);  
    /*    
    Int_t splitlevel = 0 ;
    sprintf(branchname,"AliFMDSDigitizer");   
    AliFMDSDigitizer * sd = this ;
    gAlice->MakeBranchInTree(gAlice->TreeS(),branchname,"AliFMDSDigitizer",&sd, bufferSize, splitlevel,file); 
    */

    //Now made SDigits from hits, for PHOS it is the same
  Int_t volume,sector,ring,charge;
  Float_t e;
  Float_t de[10][20][150];
  Int_t ivol,isec,iring;
  Int_t hit,nbytes;
  TParticle *particle;
  AliFMDhit  *fmdHit;

 
 // Event ------------------------- LOOP  
 
  for (ivol=1; ivol<=6; ivol++)
    for(isec=0; isec<16; isec++)
      for(iring=0; iring<128; iring++)
	de[ivol][isec][iring]=0;
      
  if (FMD){
    FMDhits   = FMD->Hits();
    TTree *TH = gAlice->TreeH();
    Stat_t ntracks    = TH->GetEntries();
 
   for (Int_t track=0; track<ntracks;track++) {
     gAlice->ResetHits();
     nbytes += TH->GetEvent(track);
     particle=gAlice->Particle(track);
     Int_t nhits =FMDhits->GetEntriesFast();
     
     for (hit=0;hit<nhits;hit++) {
       fmdHit   = (AliFMDhit*)FMDhits->UncheckedAt(hit);
       
       volume = fmdHit->Volume();
       sector = fmdHit->NumberOfSector();
       ring   = fmdHit->NumberOfRing();
       e = fmdHit->Edep();
       de[volume][sector][ring]=de[volume][sector][ring]+e;

     } //hit loop
   } //track loop
  }//if FMD
  Int_t digit[5];
  Float_t I=1.664*0.04*2.33/22400; // = 0.69e-6;
  for(ivol=1; ivol<=6; ivol++){
    for(isec=0; isec<16; isec++){
      for(iring=0; iring<128; iring++){
	if(de[ivol][isec][iring]>0.){
	  digit[0]=ivol;
	  digit[1]=isec;
	  digit[2]=iring;
	  charge=Int_t (de[ivol][isec][iring]/I);
	  digit[3]=charge;

	  //dinamic diapason from MIP(0.155MeV) to 30MIP(4.65MeV)
	  //1024 ADC channels 
	  Float_t channelWidth=(22400*30)/1024;
	  digit[4]=Int_t(digit[3]/channelWidth);

	  new((*sdigits)[nSdigits++]) AliFMDdigit(digit) ;
	  
	} //de >threshold
      }// iring loop
    }//sector loop
  } // volume loop


  
  
    gAlice->TreeS()->Fill() ;
    TFile *f1 = new TFile(file,"RECREATE");
    f1->cd();
    gAlice->TreeS()->Write(0,TObject::kOverwrite) ;
  }

  delete sdigits ;
  if(file)
    file->Close() ;

}
//__________________________________________________________________
void AliFMDSDigitizer::SetSDigitsFile(char * file ){
  if(!fSDigitsFile.IsNull())
    cout << "Changing SDigits file from " <<(char *)fSDigitsFile.Data() << " to " << file << endl ;
  fSDigitsFile=file ;
}
//__________________________________________________________________
void AliFMDSDigitizer::Print(Option_t* option)const{
  cout << "------------------- "<< GetName() << " -------------" << endl ;
  if(fSDigitsFile.IsNull())
    cout << " Writing SDigitis to file galice.root "<< endl ;
  else
    cout << "    Writing SDigitis to file  " << (char*) fSDigitsFile.Data() << endl ;

}
