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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                          V-Zero   Detector                            //
//  This class contains the base procedures for the VZERO  detector      //
//  Default geometry of November 2003 :   V0R box is 4.4 cm thick        //
//                                  scintillators are 2 cm thick         //
//  All comments should be sent to Brigitte CHEYNIS :                    //
//                                 b.cheynis@ipnl.in2p3.fr               //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


// --- Standard libraries ---
#include <Riostream.h>
#include <stdlib.h>

// --- ROOT libraries ---
#include <TNamed.h>
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliVZERO.h"
#include "AliVZEROLoader.h"
#include "AliVZERODigitizer.h"
#include "AliVZEROBuffer.h"
#include "AliRunDigitizer.h"
#include "AliVZEROdigit.h"
#include "AliDAQ.h"
#include "AliRawReader.h"
#include "AliVZERORawStream.h"

ClassImp(AliVZERO)
 //__________________________________________________________________
AliVZERO::AliVZERO(): AliDetector(),
          fIdSens1(0),
          fThickness(0.),
	  fThickness1(0.),
	  fMaxStepQua(0.),
	  fMaxStepAlu(0.),
	  fMaxDestepQua(0.),
	  fMaxDestepAlu(0.)
{
/// Default Constructor
    
    AliDebug(1,Form("default (empty) ctor this=%p",this));
    fIshunt          = 0;	  
}
//_____________________________________________________________________________
AliVZERO::AliVZERO(const char *name, const char *title)
       : AliDetector(name,title),
         fIdSens1(0),
         fThickness(4.4),
	 fThickness1(2.0),
	 fMaxStepQua(0.05),
	 fMaxStepAlu(0.01),
	 fMaxDestepQua(-1.0),
	 fMaxDestepAlu(-1.0)
{
  
  // Standard constructor for VZERO Detector
  
  AliDebug(1,Form("ctor this=%p",this));
  
  //  fIshunt       =  1;  // All hits are associated with primary particles  
   
  fHits         =  new TClonesArray("AliVZEROhit", 400);
  fDigits       =  new TClonesArray("AliVZEROdigit",400); 
   
  gAlice->GetMCApp()->AddHitList(fHits);

//   fThickness    =  4.4;   // total thickness of the V0R box in cm
//   fThickness1   =  2.0;   // thickness of scintillating cells in cm
//   
//   fMaxStepQua   =  0.05; 
//   fMaxStepAlu   =  0.01; 
//   
//   fMaxDestepQua =  -1.0;
//   fMaxDestepAlu =  -1.0;
  
  
}

//_____________________________________________________________________________
AliVZERO::~AliVZERO()
{
  //
  // Default destructor for VZERO Detector
  //
  
    if (fHits) {
        fHits->Delete();
        delete fHits;
	fHits=0; }
    
    if (fDigits) {
        fDigits->Delete();
        delete fDigits;
        fDigits=0; }
}

//_____________________________________________________________________________
void AliVZERO::CreateGeometry()
{
  //
  // Builds simple Geant3 geometry 
  //
}
//_____________________________________________________________________________
void AliVZERO::CreateMaterials()
{
  //
  // Creates materials used for Geant3 geometry 
  //
}

//_____________________________________________________________________________
void AliVZERO::Init()
{
  //
  // Initialises the VZERO  class after it has been built
  //
}


//_____________________________________________________________________________
void AliVZERO::SetMaxStepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials
  //
     fMaxStepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxStepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum foils (not used in 
  // version v2)
  //
    fMaxStepAlu = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials (quartz)
  //
    fMaxDestepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum (not used in 
  // version v2)
  //
    fMaxDestepAlu = p1;
}

//_____________________________________________________________________________
AliLoader* AliVZERO::MakeLoader(const char* topfoldername)
{ 
  //
  // Builds VZEROgetter (AliLoader type)
  // if detector wants to use customized getter, it must overload this method
  //
//  Info("MakeLoader","Creating AliVZEROLoader. Top folder is %s.",topfoldername); 
 
  AliDebug(1,Form("Creating AliVZEROLoader, Top folder is %s ",topfoldername));
  fLoader = new AliVZEROLoader(GetName(),topfoldername);
  return fLoader;
}

//_____________________________________________________________________________
void AliVZERO::SetTreeAddress()
{
  //
  // Sets tree address for hits.
  //
  if (fLoader->TreeH() && (fHits == 0x0))
    fHits = new  TClonesArray("AliVZEROhit", 400);

  AliDetector::SetTreeAddress();
}

//_____________________________________________________________________________
AliDigitizer* AliVZERO::CreateDigitizer(AliRunDigitizer* manager) const
{
  //
  // Creates a digitizer for VZERO
  //
  return new AliVZERODigitizer(manager);
}

//_____________________________________________________________________________
void AliVZERO::Hits2Digits(){
  //
  // Converts hits to digits of the current event
  //
  // Inputs file name
  const char *alifile = "galice.root";   

  // Create the run digitizer 
  AliRunDigitizer* manager = new AliRunDigitizer(1, 1);
  manager->SetInputStream(0, alifile);
  manager->SetOutputFile("H2Dfile");

  // Creates the VZERO digitizer 
  AliVZERODigitizer* dig = new AliVZERODigitizer(manager);

  // Creates the digits
  dig->Exec("");

}
//_____________________________________________________________________________
void AliVZERO::Digits2Raw()
{
   //
   //  Converts digits of the current event to raw data
   //
  
   AliVZERO *fVZERO = (AliVZERO*)gAlice->GetDetector("VZERO");
   fLoader->LoadDigits();
   TTree* digits = fLoader->TreeD();
   if (!digits) {
      Error("Digits2Raw", "no digits tree");
      return;
   }
   TClonesArray * VZEROdigits = new TClonesArray("AliVZEROdigit",1000);
   fVZERO->SetTreeAddress();  		
   digits->GetBranch("VZERODigit")->SetAddress(&VZEROdigits); 
  
   const char *fileName    = AliDAQ::DdlFileName("VZERO",0);
   AliVZEROBuffer* buffer  = new AliVZEROBuffer(fileName);
  
   //  Verbose level
   //  0: Silent
   //  1: cout messages
   //  2: txt files with digits 
   //  BE CAREFUL, verbose level 2 MUST be used only for debugging and
   //  it is highly suggested to use this mode only for debugging digits files
   //  reasonably small, because otherwise the size of the txt files can reach
   //  quickly several MB wasting time and disk space.
  
   ofstream ftxt;
   buffer->SetVerbose(0);
   Int_t verbose = buffer->GetVerbose();

   // Get Trigger information first
   // Read trigger inputs from trigger-detector object
   AliDataLoader * dataLoader = fLoader->GetDigitsDataLoader();
   if( !dataLoader->IsFileOpen() ) 
        dataLoader->OpenFile( "READ" );
   AliTriggerDetector* trgdet = (AliTriggerDetector*)dataLoader->GetDirectory()->Get( "Trigger" );
   UInt_t triggerInfo = 0;
   if(trgdet) {
      triggerInfo = trgdet->GetMask() & 0xffff;
   }
   else {
      AliError(Form("There is no trigger object for %s",fLoader->GetName()));
   }

   buffer->WriteTriggerInfo((UInt_t)triggerInfo); 
   buffer->WriteTriggerScalers(); 
   buffer->WriteBunchNumbers(); 
  
   // Now retrieve the channel information: charge+time and 
   // dump it into ADC and Time arrays
   // We assume here an ordered (by PMNumber) array of
   // digits!!

   Int_t nEntries = Int_t(digits->GetEntries());
   UInt_t ADC[64];
   Int_t PMNumber[64];
   UInt_t Time[64];
  
   for (Int_t i = 0; i < nEntries; i++) {
     fVZERO->ResetDigits();
     digits->GetEvent(i);
     Int_t ndig = VZEROdigits->GetEntriesFast(); 
   
     if(ndig == 0) continue;
     if(verbose == 2) {ftxt.open("VZEROdigits.txt",ios::app);}
     for(Int_t k=0; k<ndig; k++){
         AliVZEROdigit* fVZERODigit = (AliVZEROdigit*) VZEROdigits->At(k);
         // Convert aliroot channel k into FEE channel iChannel before writing data
	 Int_t iChannel = buffer->GetOnlineChannel(k);		
	 ADC[iChannel]       = (UInt_t) fVZERODigit->ADC();
	 PMNumber[iChannel]  = fVZERODigit->PMNumber();
	 Time[iChannel]      = (UInt_t) fVZERODigit->Time();
         if(verbose == 1) { cout <<"DDL: "<<fileName<< "\tdigit number: "<< k<<"\tPM number: "
	                    <<PMNumber[k]<<"\tADC: "<< ADC[k] << "\tTime: "<< Time[k] << endl;} 
	 if(verbose == 2) {
	      ftxt<<"DDL: "<<fileName<< "\tdigit number: "<< k<<"\tPM number: "
	                   <<PMNumber[k]<<"\tADC: "<< ADC[k] << "\tTime: "<< Time[k] << endl;}	      
//	 printf("DDL: %s, channel: %d, PM: %d, ADC: %d, Time: %d \n", 
//	            fileName,k,PMNumber[k],ADC[k],Time[k]); 
     }        
   if(verbose==2) ftxt.close();
   }

   // Now fill raw data
          
   for (Int_t  iCIU = 0; iCIU < 8; iCIU++) { 
 
   // decoding of one Channel Interface Unit numbered iCIU - there are 8 channels per CIU (and 8 CIUs) :
  
      for(Int_t iChannel_Offset = iCIU*8; iChannel_Offset < (iCIU*8)+8; iChannel_Offset=iChannel_Offset+4) { 
         for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
             buffer->WriteChannel(iChannel, ADC[iChannel], Time[iChannel]);       
         }
         buffer->WriteBeamFlags(); 
         buffer->WriteMBInfo(); 
         buffer->WriteMBFlags();   
         buffer->WriteBeamScalers(); 
      } 
      for(Int_t iChannel=0; iChannel < 8; iChannel++) {
          buffer->WriteTiming(iChannel, ADC[iChannel], Time[iChannel]); 
      }

    // End of decoding of one CIU card
    
  } // end of decoding the eight CIUs
     
  delete buffer;
  fLoader->UnloadDigits();  
}

//_____________________________________________________________________________
Bool_t AliVZERO::Raw2SDigits(AliRawReader* rawReader){
  // Converts the VZERO raw data into digits
  // The method is used for merging simulated and
  // real data events
  TStopwatch timer;
  timer.Start();

  if(!fLoader) {
    AliError("no VZERO loader found");
    return kFALSE; }

  TTree* treeD  = fLoader->TreeD();
  if(!treeD) {
      fLoader->MakeTree("D");
      treeD = fLoader->TreeD(); }
        
  AliVZEROdigit  digit;
  AliVZEROdigit* pdigit = &digit;
  const Int_t kBufferSize = 4000;
   
  treeD->Branch("VZERO", "AliVZEROdigit",  &pdigit, kBufferSize);

  rawReader->Reset();
  AliVZERORawStream* rawStream  = new AliVZERORawStream(rawReader);    
     
  if (!rawStream->Next()) return kFALSE; // No VZERO data found
  
  for(Int_t i=0; i<64; i++) {
      new(pdigit) AliVZEROdigit(i, (Int_t)rawStream->GetADC(i), (Int_t)rawStream->GetTime(i)); 
      treeD->Fill();
  }
 
// Checks if everything is OK by printing results 

//   for(int i=0;i<64;i++) {
// 	printf("Channel %d : %d %d \n",i,rawStream->GetADC(i),rawStream->GetTime(i)); }
//   treeD->Print(); printf(" \n"); 
   	
  fLoader->WriteDigits("OVERWRITE");
  fLoader->UnloadDigits();	
	
  delete rawStream;

  timer.Stop();
  timer.Print();
  return kTRUE;
}


