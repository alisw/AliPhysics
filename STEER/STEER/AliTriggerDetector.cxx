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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Base Class for Detector specific Trigger                                 //                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TNamed.h>
#include <TString.h>
#include <TObjArray.h>
#include <TRandom.h>


#include "AliLog.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliLoader.h"

#include "AliTriggerInput.h"
#include "AliTriggerDetector.h"


using std::endl;
using std::cout;
using std::hex;
using std::dec;
ClassImp( AliTriggerDetector )

//_____________________________________________________________________________
AliTriggerDetector::AliTriggerDetector() :
  TNamed(),
  fMask(0),
  fInputs()
{
   // Default constructor
}

//_____________________________________________________________________________
AliTriggerDetector::AliTriggerDetector(const AliTriggerDetector & de ):
  TNamed(de),
  fMask(de.fMask),
  fInputs(de.fInputs)
{
  // Copy constructor
  
}

//_____________________________________________________________________________
AliTriggerDetector::~AliTriggerDetector()
{
  // Destructor
  // Delete only inputs that are not
  // associated to the CTP
  Int_t ninputs = fInputs.GetEntriesFast();
  for(Int_t i = 0; i < ninputs; i++) {
    AliTriggerInput *inp = (AliTriggerInput *)fInputs.At(i);
    if (inp->GetSignature() == -1 &&
	inp->GetMask() == 0)
      delete fInputs.RemoveAt(i);
  }
}

//_____________________________________________________________________________
void AliTriggerDetector::AssignInputs(const TObjArray &inputs)
{
  // Cross-check the trigger inputs provided by the detector trigger
  // processor (TP) and the inputs defined in CTP
  // a) If the input is defined in the TP, but not in CTP it
  // will be generated but not used by CTP. It will be possibly stored
  // in the detector raw data if the hardware allows this.
  // b) If hte input is not defined in the TP, but is defined in CTP
  // then a warning is issued and the CTP simulation is working
  // with this input disabled.
 
   // Check if we have to create the inputs first
   if( fInputs.GetEntriesFast() == 0 ) {
     // Create the inputs that the detector can provide
     CreateInputs();
   }

   TString name = GetName();

   // Pointer to the available inputs provided by the detector
   TObjArray* availInputs = &fInputs;

   Int_t ninputs = inputs.GetEntriesFast();
   for( Int_t j=0; j<ninputs; j++ ) {
     AliTriggerInput *inp = (AliTriggerInput*)inputs.At(j);
     if ( name.CompareTo(inp->GetModule().Data()) ) continue;
     AliTriggerInput *tempObj = (AliTriggerInput*)availInputs->FindObject(inp->GetInputName().Data());
     if ( tempObj ) {
       Int_t tempIndex = availInputs->IndexOf(tempObj);
       if (tempObj->GetSignature()!=-1) {
	 AliErrorF("Ignore attempt to overwrite input%d %s [Sign:%d, Mask:%llu] by [Sign:%d, Mask:%llu]",
		   tempIndex,tempObj->GetInputName().Data(),
		   tempObj->GetSignature(),tempObj->GetMask(),inp->GetSignature(),inp->GetMask());
       }
       else {
	 delete availInputs->Remove(tempObj);
	 fInputs.AddAt( inp, tempIndex );
	 inp->Enable();
	 AliDebug(1,Form("Trigger input (%s) is found in the CTP configuration. Therefore it is enabled for trigger detector (%s)",
			 inp->GetInputName().Data(),name.Data()));
       }
     }
     else {
       AliWarning(Form("Trigger Input (%s) is not implemented for the trigger detector (%s) ! It will be disabled !",
		       inp->GetInputName().Data(),name.Data()));
     }
   }

   for( Int_t j=0; j<fInputs.GetEntriesFast(); j++ ) {
     AliTriggerInput *inp = (AliTriggerInput *)fInputs.At(j);
     if (inp->GetSignature() == -1 &&
	 inp->GetMask() == 0) {
       inp->Enable();
       AliDebug(1,Form("Trigger input (%s) was not found in the CTP configuration. Therefore it will be run in a stand-alone mode",
		    inp->GetInputName().Data()));
     }
   }

   fInputs.SetOwner(kFALSE);
}

//_____________________________________________________________________________
void AliTriggerDetector::CreateInputs()
{
   // Define the inputs to the Central Trigger Processor
   // This is a dummy version 
   
   // Do not create inputs again!!
   if( fInputs.GetEntriesFast() > 0 ) return;
   
   fInputs.AddLast( new AliTriggerInput( "TEST1_L0", "TEST", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "TEST2_L0", "TEST", 0 ) );
   fInputs.AddLast( new AliTriggerInput( "TEST3_L0", "TEST", 0 ) );
}

//_____________________________________________________________________________
void AliTriggerDetector::Trigger()
{
   // This is a dummy version set all inputs in a random way

   AliWarning( Form( "Triggering dummy detector %s", GetName() ) );

   //  ********** Get Digits for the current event **********
   AliRunLoader* runLoader = AliRunLoader::Instance();
   AliInfo( Form( "Event %d", runLoader->GetEventNumber() ) );
   
   TString loadername = GetName(); 
   loadername.Append( "Loader" );
   AliLoader * loader = runLoader->GetLoader( loadername.Data() );
   if( loader ) {
      loader->LoadDigits( "READ" );
      TTree* digits = loader->TreeD();
      // Do something with the digits !!!
      if( digits ) { 
//         digits->Print();
      }   
      loader->UnloadDigits();
   }
   //  ******************************************************

   // set all inputs in a random way
   Int_t nInputs = fInputs.GetEntriesFast();
   for( Int_t j=0; j<nInputs; j++ ) {
      AliTriggerInput* in = (AliTriggerInput*)fInputs.At( j );
      if( gRandom->Rndm() > 0.5 ) {
         in->Set();
         fMask |= in->GetValue();
      }
   }
}

//_____________________________________________________________________________
void AliTriggerDetector::SetInput( TString& name )
{
   // Set Input by name
   SetInput( name.Data() );
}

//_____________________________________________________________________________
void AliTriggerDetector::SetInput( const char * name )
{
   // Set Input by name
   AliTriggerInput* in = ((AliTriggerInput*)fInputs.FindObject( name ));
   if( in ) {
      in->Set();
      fMask |= in->GetValue();
   } else
      AliError( Form( "There is not input named %s", name ) );
}

//_____________________________________________________________________________
void AliTriggerDetector::Print( const Option_t* opt ) const
{
   // Print
   cout << "Trigger Detector : " << GetName() << endl;
   cout << "  Trigger Class Mask: 0x" << hex << GetMask() << dec << endl;
   Int_t nInputs = fInputs.GetEntriesFast();
   for( Int_t j=0; j<nInputs; j++ ) {
      AliTriggerInput* in = (AliTriggerInput*)fInputs.At( j );
      in->Print( opt );
   }
}
